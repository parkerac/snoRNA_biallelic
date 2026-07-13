#!/usr/bin/env python3
"""Count snoRNA-overlapping variants gene by gene from AGGV3 shard BEDs."""

import argparse
import csv
import gzip
import os
import re
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed
from concurrent.futures import ThreadPoolExecutor


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def safe_name(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return value.strip("._") or "unknown"


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
            if len(parts) < 6:
                continue
            chrom = parts[0]
            start = int(parts[1]) + 1
            end = int(parts[2])
            shard = parts[4] if len(parts) > 4 else ""
            subshard = parts[5] if len(parts) > 5 else ""
            vcf_path = parts[6] if len(parts) > 6 else ""
            shard_bed[chrom].append((start, end, shard, subshard, vcf_path))
    for chrom in shard_bed:
        shard_bed[chrom].sort()
    return shard_bed


def gene_overlaps_shard(gene, shard):
    return not (shard[1] < gene["start"] or shard[0] > gene["end"])


def gene_overlaps_variant(gene, pos):
    return gene["start"] <= pos <= gene["end"]


def find_overlapping_vcfs_for_gene(gene, shard_bed, vcf_root=None, vcf_template="shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz"):
    hits = []
    for shard in shard_bed.get(gene["chrom"], []):
        if not gene_overlaps_shard(gene, shard):
            continue
        if vcf_root:
            path = os.path.join(vcf_root, vcf_template.format(shard=shard[2], subshard=shard[3]))
        else:
            path = shard[4]
        if path:
            hits.append(path)
    return sorted(set(hits))


def is_carrier(gt):
    return any(allele not in {".", "0"} for allele in gt.replace("|", "/").split("/"))


def read_participants(path, id_col, group_col):
    rows = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if id_col not in reader.fieldnames:
            raise ValueError(f"Missing participant id column: {id_col}")
        if group_col and group_col not in reader.fieldnames:
            raise ValueError(f"Missing group column: {group_col}")
        for row in reader:
            rows[row[id_col]] = row
    return rows


def detail_path_for_gene(detail_dir, gene_index, gene):
    return os.path.join(
        detail_dir,
        f"{gene_index:04d}_{safe_name(gene['gene_name'])}_{safe_name(gene['gene_id'])}.tsv",
    )


def process_gene(task):
    gene = task["gene"]
    vcf_paths = task["vcf_paths"]
    detail_path = task["detail_path"]
    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
    seen_variants = set()

    with open(detail_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["gene_name", "gene_id", "rna_class", "chrom", "pos", "ref", "alt", "variant_id", "participant_id"],
        )
        writer.writeheader()
        for path in vcf_paths:
            with open_text(path) as vcf:
                samples = []
                for line in vcf:
                    if line.startswith("##"):
                        continue
                    if line.startswith("#CHROM"):
                        samples = line.rstrip("\n").split("\t")[9:]
                        continue
                    parts = line.rstrip("\n").split("\t")
                    if len(parts) < 10 or not samples or parts[0] != gene["chrom"]:
                        continue
                    pos = int(parts[1])
                    if not gene_overlaps_variant(gene, pos):
                        continue
                    fmt = parts[8].split(":")
                    if "GT" not in fmt:
                        continue
                    variant_id = f"{parts[0]}:{pos}:{parts[3]}:{parts[4]}"
                    if variant_id in seen_variants:
                        continue
                    gt_idx = fmt.index("GT")
                    carrier_rows = []
                    for sample, sample_field in zip(samples, parts[9:]):
                        fields = sample_field.split(":")
                        if len(fields) <= gt_idx:
                            continue
                        gt = fields[gt_idx]
                        if gt in {".", "./.", ".|."} or not is_carrier(gt):
                            continue
                        carrier_rows.append(sample)
                    if not carrier_rows:
                        seen_variants.add(variant_id)
                        continue
                    for sample in carrier_rows:
                        participants[sample]["variants"].add(variant_id)
                        participants[sample]["genes"].add((gene["gene_name"], gene["gene_id"]))
                        writer.writerow(
                            {
                                "gene_name": gene["gene_name"],
                                "gene_id": gene["gene_id"],
                                "rna_class": gene["rna_class"],
                                "chrom": parts[0],
                                "pos": pos,
                                "ref": parts[3],
                                "alt": parts[4],
                                "variant_id": variant_id,
                                "participant_id": sample,
                            }
                        )
                    seen_variants.add(variant_id)

    return {
        "gene_index": task["gene_index"],
        "gene_name": gene["gene_name"],
        "gene_id": gene["gene_id"],
        "rna_class": gene["rna_class"],
        "chrom": gene["chrom"],
        "start": gene["start"],
        "end": gene["end"],
        "detail_path": detail_path,
        "n_vcfs": len(vcf_paths),
        "n_variants": len(seen_variants),
        "n_participants": len(participants),
        "participants": {pid: {"variants": data["variants"], "genes": data["genes"]} for pid, data in participants.items()},
    }


def merge_participants(dest, src):
    for pid, data in src.items():
        dest[pid]["variants"].update(data["variants"])
        dest[pid]["genes"].update(data["genes"])


def write_gene_summary(out_prefix, gene_stats):
    path = f"{out_prefix}.gene_summary.tsv"
    print(f"Writing gene summary to {path}", flush=True)
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "gene_index",
                "gene_name",
                "gene_id",
                "rna_class",
                "chrom",
                "start",
                "end",
                "n_vcfs",
                "n_variants",
                "n_participants",
                "detail_path",
            ],
        )
        writer.writeheader()
        for row in sorted(gene_stats, key=lambda r: r["gene_index"]):
            writer.writerow(row)


def write_participant_outputs(out_prefix, participants, participant_rows, group_col):
    part_path = f"{out_prefix}.participants.tsv"
    summary_path = f"{out_prefix}.summary.tsv"
    print(f"Writing participant details to {part_path}", flush=True)
    with open(part_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["participant_id", "group", "n_variants", "n_genes", "variants", "genes"],
        )
        writer.writeheader()
        for pid in sorted(participant_rows):
            row = participant_rows[pid]
            data = participants.get(pid, {"variants": set(), "genes": set()})
            writer.writerow(
                {
                    "participant_id": pid,
                    "group": row.get(group_col, "") if group_col else "",
                    "n_variants": len(data["variants"]),
                    "n_genes": len(data["genes"]),
                    "variants": ";".join(sorted(data["variants"])),
                    "genes": ";".join(sorted(f"{name}|{gene_id}" for name, gene_id in data["genes"])),
                }
            )

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
                n_variants = len(participants.get(pid, {"variants": set()})["variants"])
                by_group[group]["n_participants"] += 1
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


def execute_tasks(tasks, workers, executor_cls):
    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
    gene_stats = []
    with executor_cls(max_workers=workers) as pool:
        futures = [pool.submit(process_gene, task) for task in tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            merge_participants(participants, result["participants"])
            gene_stats.append({k: result[k] for k in ("gene_index", "gene_name", "gene_id", "rna_class", "chrom", "start", "end", "n_vcfs", "n_variants", "n_participants", "detail_path")})
            print(
                f"[done {done}/{len(tasks)}] {result['gene_name']} ({result['gene_id']}): "
                f"{result['n_variants']} variants, {result['n_participants']} participants",
                flush=True,
            )
    return participants, gene_stats


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True)
    parser.add_argument("--shard-bed", required=True)
    parser.add_argument("--vcf-root", required=True)
    parser.add_argument("--vcf-template", default="shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz")
    parser.add_argument("--participant-tsv", required=True)
    parser.add_argument("--participant-id-col", required=True)
    parser.add_argument("--group-col")
    parser.add_argument("--feature-type", action="append", default=["snoRNA"])
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--cpus", type=int, default=16)
    args = parser.parse_args()

    feature_types = set(args.feature_type)
    print(f"Loading genes from {args.gtf}", flush=True)
    genes = parse_gtf(args.gtf, feature_types)
    print(f"Loaded {len(genes)} genes", flush=True)
    print(f"Loading shard BED from {args.shard_bed}", flush=True)
    shard_bed = parse_shard_bed(args.shard_bed)
    participant_rows = read_participants(args.participant_tsv, args.participant_id_col, args.group_col)

    out_prefix_dir = os.path.dirname(args.out_prefix)
    if out_prefix_dir:
        os.makedirs(out_prefix_dir, exist_ok=True)
    detail_dir = f"{args.out_prefix}.genes"
    os.makedirs(detail_dir, exist_ok=True)

    tasks = []
    for gene_index, gene in enumerate(genes, start=1):
        vcf_paths = find_overlapping_vcfs_for_gene(gene, shard_bed, vcf_root=args.vcf_root, vcf_template=args.vcf_template)
        detail_path = detail_path_for_gene(detail_dir, gene_index, gene)
        print(f"[queued {gene_index}/{len(genes)}] {gene['gene_name']} ({gene['gene_id']}): {len(vcf_paths)} VCFs", flush=True)
        tasks.append({"gene_index": gene_index, "gene": gene, "vcf_paths": vcf_paths, "detail_path": detail_path})

    workers = max(1, min(args.cpus, len(tasks)))
    print(f"Processing {len(tasks)} genes with {workers} workers", flush=True)

    try:
        participants, gene_stats = execute_tasks(tasks, workers, ProcessPoolExecutor)
    except (PermissionError, OSError) as exc:
        print(f"Process pool unavailable ({exc}); falling back to threads", flush=True)
        participants, gene_stats = execute_tasks(tasks, workers, ThreadPoolExecutor)

    write_gene_summary(args.out_prefix, gene_stats)
    write_participant_outputs(args.out_prefix, participants, participant_rows, args.group_col)


if __name__ == "__main__":
    main()
