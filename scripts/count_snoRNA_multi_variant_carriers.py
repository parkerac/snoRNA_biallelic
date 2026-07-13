#!/usr/bin/env python3
"""Count snoRNA-overlapping variant sites per participant from VCF shards."""

import argparse
import csv
import gzip
import os
from collections import defaultdict


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def parse_gtf(path, feature_types):
    genes = []
    with open_text(path) as fh:
        for line in fh:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 9:
                continue
            chrom, _, _, start, end, _, _, _, attrs = parts[:9]
            attr_map = {}
            for item in attrs.split(";"):
                item = item.strip()
                if not item or " " not in item:
                    continue
                key, value = item.split(" ", 1)
                attr_map[key] = value.strip().strip('"')
            gene_type = attr_map.get("gene_type") or attr_map.get("gene_biotype")
            if gene_type not in feature_types:
                continue
            genes.append(
                {
                    "chrom": chrom,
                    "start": int(start),
                    "end": int(end),
                    "rna_class": gene_type,
                    "gene_name": attr_map.get("gene_name", "UNKNOWN"),
                    "gene_id": attr_map.get("gene_id", "UNKNOWN"),
                }
            )
    return sorted(genes, key=lambda g: (g["chrom"], g["start"], g["end"], g["gene_name"], g["gene_id"]))


def parse_shard_bed(path):
    shard_bed = defaultdict(list)
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 8:
                continue
            chrom = parts[0]
            start = int(parts[1]) + 1
            end = int(parts[2])
            shard = parts[4]
            subshard = parts[5]
            vcf_path = parts[6]
            shard_bed[chrom].append((start, end, shard, subshard, vcf_path))
    for chrom in shard_bed:
        shard_bed[chrom].sort()
    return shard_bed


def gene_overlaps_shard(gene, shard):
    return not (shard[1] < gene["start"] or shard[0] > gene["end"])


def gene_overlaps_variant(gene, pos):
    return gene["start"] <= pos <= gene["end"]


def find_overlapping_vcfs_for_gene(gene, shard_bed, vcf_root=None, vcf_template="shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz"):
    hits = set()
    for shard in shard_bed.get(gene["chrom"], []):
        if not gene_overlaps_shard(gene, shard):
            continue
        path = os.path.join(vcf_root, vcf_template.format(shard=shard[2], subshard=shard[3])) if vcf_root else shard[4]
        hits.add(path)
    return sorted(hits)


def is_carrier(gt):
    return any(allele not in {".", "0"} for allele in gt.replace("|", "/").split("/"))


def read_participants(path, id_col, group_col):
    rows = {}
    groups = defaultdict(int)
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if id_col not in reader.fieldnames:
            raise ValueError(f"Missing participant id column: {id_col}")
        if group_col and group_col not in reader.fieldnames:
            raise ValueError(f"Missing group column: {group_col}")
        for row in reader:
            pid = row[id_col]
            rows[pid] = row
            if group_col:
                groups[row[group_col]] += 1
    return rows, groups


def scan_gene_vcfs(gene, vcf_paths):
    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
    site_rows = []
    variant_count = 0
    for path in vcf_paths:
        with open_text(path) as fh:
            samples = []
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10 or not samples:
                    continue
                if parts[0] != gene["chrom"]:
                    continue
                pos = int(parts[1])
                if not gene_overlaps_variant(gene, pos):
                    continue
                ref, alts = parts[3], parts[4]
                fmt = parts[8].split(":")
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")
                variant_id = f"{parts[0]}:{pos}:{ref}:{alts}"
                variant_count += 1
                for sample, sample_field in zip(samples, parts[9:]):
                    fields = sample_field.split(":")
                    if len(fields) <= gt_idx:
                        continue
                    gt = fields[gt_idx]
                    if gt in {".", "./.", ".|."} or not is_carrier(gt):
                        continue
                    participants[sample]["variants"].add(variant_id)
                    participants[sample]["genes"].add((gene["gene_name"], gene["gene_id"]))
                    site_rows.append(
                        {
                            "gene_name": gene["gene_name"],
                            "gene_id": gene["gene_id"],
                            "rna_class": gene["rna_class"],
                            "chrom": parts[0],
                            "pos": pos,
                            "ref": ref,
                            "alt": alts,
                            "variant_id": variant_id,
                            "participant_id": sample,
                        }
                    )
    return site_rows, participants, variant_count


def scan_all_vcfs(vcf_dir, genes):
    by_chrom = defaultdict(list)
    for gene in genes:
        by_chrom[gene["chrom"]].append(gene)
    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
    site_rows = []
    variant_count = 0
    for name in sorted(os.listdir(vcf_dir)):
        if not name.endswith((".vcf", ".vcf.gz")):
            continue
        path = os.path.join(vcf_dir, name)
        with open_text(path) as fh:
            samples = []
            for line in fh:
                if line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.rstrip("\n").split("\t")[9:]
                    continue
                parts = line.rstrip("\n").split("\t")
                if len(parts) < 10 or not samples or parts[0] not in by_chrom:
                    continue
                pos = int(parts[1])
                ref, alts = parts[3], parts[4]
                hits = [gene for gene in by_chrom[parts[0]] if gene_overlaps_variant(gene, pos)]
                if not hits:
                    continue
                fmt = parts[8].split(":")
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")
                variant_id = f"{parts[0]}:{pos}:{ref}:{alts}"
                variant_count += 1
                for sample, sample_field in zip(samples, parts[9:]):
                    fields = sample_field.split(":")
                    if len(fields) <= gt_idx:
                        continue
                    gt = fields[gt_idx]
                    if gt in {".", "./.", ".|."} or not is_carrier(gt):
                        continue
                    participants[sample]["variants"].add(variant_id)
                    for gene in hits:
                        participants[sample]["genes"].add((gene["gene_name"], gene["gene_id"]))
                        site_rows.append(
                            {
                                "gene_name": gene["gene_name"],
                                "gene_id": gene["gene_id"],
                                "rna_class": gene["rna_class"],
                                "chrom": parts[0],
                                "pos": pos,
                                "ref": ref,
                                "alt": alts,
                                "variant_id": variant_id,
                                "participant_id": sample,
                            }
                        )
    return site_rows, participants, variant_count


def write_outputs(out_prefix, participants, participant_rows, group_col, gene_stats, gene_out_path=None):
    if gene_out_path is None:
        gene_out_path = f"{out_prefix}.gene_summary.tsv"
    print(f"Writing gene summary to {gene_out_path}", flush=True)
    with open(gene_out_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["gene_name", "gene_id", "rna_class", "chrom", "start", "end", "n_vcfs", "n_variants", "n_participants"],
        )
        writer.writeheader()
        for row in gene_stats:
            writer.writerow(row)

    summary_path = f"{out_prefix}.summary.tsv"
    print(f"Writing participant summary to {summary_path}", flush=True)
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["group", "n_participants", "n_multiple_variant_carriers", "n_unique_variants"],
        )
        writer.writeheader()
        if group_col:
            by_group = defaultdict(lambda: {"n_participants": 0, "n_multiple_variant_carriers": 0, "n_unique_variants": 0})
            for pid, row in participant_rows.items():
                group = row.get(group_col, "")
                by_group[group]["n_participants"] += 1
                n_variants = len(participants.get(pid, {"variants": set()})["variants"])
                by_group[group]["n_unique_variants"] += n_variants
                if n_variants >= 2:
                    by_group[group]["n_multiple_variant_carriers"] += 1
            for group in sorted(by_group):
                writer.writerow({"group": group, **by_group[group]})
        else:
            writer.writerow(
                {
                    "group": "all",
                    "n_participants": len(participant_rows),
                    "n_multiple_variant_carriers": sum(
                        1 for pid in participant_rows if len(participants.get(pid, {"variants": set()})["variants"]) >= 2
                    ),
                    "n_unique_variants": sum(len(participants.get(pid, {"variants": set()})["variants"]) for pid in participant_rows),
                }
            )


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--vcf-dir")
    parser.add_argument("--shard-bed")
    parser.add_argument("--vcf-root")
    parser.add_argument("--vcf-template", default="shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz")
    parser.add_argument("--participant-tsv", required=True)
    parser.add_argument("--participant-id-col", required=True)
    parser.add_argument("--group-col")
    parser.add_argument("--feature-type", action="append", default=["snoRNA"])
    parser.add_argument("--out-prefix", required=True)
    args = parser.parse_args()

    feature_types = set(args.feature_type)
    print(f"Loading genes from {args.gtf}", flush=True)
    genes = parse_gtf(args.gtf, feature_types)
    print(f"Loaded {len(genes)} genes", flush=True)
    participant_rows, _ = read_participants(args.participant_tsv, args.participant_id_col, args.group_col)

    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
    gene_stats = []
    detail_path = f"{args.out_prefix}.variants.tsv"
    print(f"Writing detailed rows to {detail_path}", flush=True)
    detail_fh = open(detail_path, "w", newline="")
    detail_writer = csv.DictWriter(
        detail_fh,
        delimiter="\t",
        fieldnames=["gene_name", "gene_id", "rna_class", "chrom", "pos", "ref", "alt", "variant_id", "participant_id"],
    )
    detail_writer.writeheader()

    try:
        if args.shard_bed:
            print(f"Loading shard BED from {args.shard_bed}", flush=True)
            shard_bed = parse_shard_bed(args.shard_bed)
            for i, gene in enumerate(genes, start=1):
                print(f"[{i}/{len(genes)}] {gene['gene_name']} ({gene['gene_id']})", flush=True)
                vcf_paths = find_overlapping_vcfs_for_gene(gene, shard_bed, vcf_root=args.vcf_root, vcf_template=args.vcf_template)
                print(f"  overlapping VCFs: {len(vcf_paths)}", flush=True)
                if not vcf_paths:
                    gene_stats.append(
                        {
                            "gene_name": gene["gene_name"],
                            "gene_id": gene["gene_id"],
                            "rna_class": gene["rna_class"],
                            "chrom": gene["chrom"],
                            "start": gene["start"],
                            "end": gene["end"],
                            "n_vcfs": 0,
                            "n_variants": 0,
                            "n_participants": 0,
                        }
                    )
                    continue
                gene_rows, gene_participants, n_variants = scan_gene_vcfs(gene, vcf_paths)
                print(f"  variants found: {n_variants}, carrier rows: {len(gene_rows)}", flush=True)
                for row in gene_rows:
                    detail_writer.writerow(row)
                detail_fh.flush()
                for pid, data in gene_participants.items():
                    participants[pid]["variants"].update(data["variants"])
                    participants[pid]["genes"].update(data["genes"])
                gene_stats.append(
                    {
                        "gene_name": gene["gene_name"],
                        "gene_id": gene["gene_id"],
                        "rna_class": gene["rna_class"],
                        "chrom": gene["chrom"],
                        "start": gene["start"],
                        "end": gene["end"],
                        "n_vcfs": len(vcf_paths),
                        "n_variants": n_variants,
                        "n_participants": len(gene_participants),
                    }
                )
        else:
            if not args.vcf_dir:
                raise SystemExit("--vcf-dir is required when --shard-bed is not provided")
            print(f"Scanning VCF directory {args.vcf_dir}", flush=True)
            gene_rows, gene_participants, _ = scan_all_vcfs(args.vcf_dir, genes)
            participants = gene_participants
            for row in gene_rows:
                detail_writer.writerow(row)
            detail_fh.flush()
            for gene in genes:
                gene_stats.append(
                    {
                        "gene_name": gene["gene_name"],
                        "gene_id": gene["gene_id"],
                        "rna_class": gene["rna_class"],
                        "chrom": gene["chrom"],
                        "start": gene["start"],
                        "end": gene["end"],
                        "n_vcfs": 0,
                        "n_variants": sum(1 for row in gene_rows if row["gene_id"] == gene["gene_id"]),
                        "n_participants": len({row["participant_id"] for row in gene_rows if row["gene_id"] == gene["gene_id"]}),
                    }
                )
    finally:
        detail_fh.close()
    write_outputs(args.out_prefix, participants, participant_rows, args.group_col, gene_stats)


if __name__ == "__main__":
    main()
