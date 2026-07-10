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
    regions = defaultdict(list)
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
            regions[chrom].append((int(start), int(end), attr_map.get("gene_name", "UNKNOWN"), attr_map.get("gene_id", "UNKNOWN")))
    for chrom in regions:
        regions[chrom].sort()
    return regions


def regions_to_bed(regions):
    bed = []
    for chrom, intervals in regions.items():
        for start, end, gene_name, gene_id in intervals:
            bed.append((chrom, start - 1, end, gene_name, gene_id))
    bed.sort()
    return bed


def find_overlapping_shards(regions, shard_bed, vcf_root=None, vcf_template="shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz"):
    hits = set()
    for chrom, start, end, *_ in regions_to_bed(regions):
        for row in shard_bed.get(chrom, []):
            shard_start, shard_end, shard, subshard = row[:4]
            if shard_end < start or shard_start > end:
                continue
            if vcf_root:
                path = os.path.join(vcf_root, vcf_template.format(shard=shard, subshard=subshard))
            else:
                path = row[4]
            hits.add(path)
    return sorted(hits)


def is_carrier(gt):
    return any(allele not in {".", "0"} for allele in gt.replace("|", "/").split("/"))


def overlaps(chrom, pos, regions):
    hits = []
    if chrom not in regions:
        return hits
    for start, end, gene_name, gene_id in regions[chrom]:
        if start <= pos <= end:
            hits.append((gene_name, gene_id))
    return hits


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


def scan_vcfs(vcf_dir, regions):
    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
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
                if len(parts) < 10 or not samples:
                    continue
                chrom, pos, ref, alts = parts[0], int(parts[1]), parts[3], parts[4]
                hits = overlaps(chrom, pos, regions)
                if not hits:
                    continue
                fmt = parts[8].split(":")
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")
                variant_id = f"{chrom}:{pos}:{ref}:{alts}"
                for sample, sample_field in zip(samples, parts[9:]):
                    fields = sample_field.split(":")
                    if len(fields) <= gt_idx:
                        continue
                    gt = fields[gt_idx]
                    if gt in {".", "./.", ".|."} or not is_carrier(gt):
                        continue
                    participants[sample]["variants"].add(variant_id)
                    for gene_name, gene_id in hits:
                        participants[sample]["genes"].add((gene_name, gene_id))
    return participants


def scan_vcf_paths(vcf_paths, regions):
    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
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
                chrom, pos, ref, alts = parts[0], int(parts[1]), parts[3], parts[4]
                hits = overlaps(chrom, pos, regions)
                if not hits:
                    continue
                fmt = parts[8].split(":")
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")
                variant_id = f"{chrom}:{pos}:{ref}:{alts}"
                for sample, sample_field in zip(samples, parts[9:]):
                    fields = sample_field.split(":")
                    if len(fields) <= gt_idx:
                        continue
                    gt = fields[gt_idx]
                    if gt in {".", "./.", ".|."} or not is_carrier(gt):
                        continue
                    participants[sample]["variants"].add(variant_id)
                    for gene_name, gene_id in hits:
                        participants[sample]["genes"].add((gene_name, gene_id))
    return participants


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
    return shard_bed


def write_outputs(out_prefix, participants, participant_rows, group_col):
    part_path = f"{out_prefix}.participants.tsv"
    summary_path = f"{out_prefix}.summary.tsv"
    with open(part_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["participant_id", "group", "n_variants", "n_genes", "variants", "genes"],
        )
        writer.writeheader()
        for pid in sorted(set(participant_rows) | set(participants)):
            data = participants.get(pid, {"variants": set(), "genes": set()})
            row = participant_rows.get(pid, {})
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
    with open(summary_path, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=["group", "n_participants", "n_multiple_variant_carriers"],
        )
        writer.writeheader()
        if group_col:
            by_group = defaultdict(lambda: {"n_participants": 0, "n_multiple_variant_carriers": 0})
            for pid, row in participant_rows.items():
                group = row.get(group_col, "")
                by_group[group]["n_participants"] += 1
                if len(participants.get(pid, {"variants": set()})["variants"]) >= 2:
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

    regions = parse_gtf(args.gtf, set(args.feature_type))
    participant_rows, _ = read_participants(args.participant_tsv, args.participant_id_col, args.group_col)
    if args.shard_bed:
        shard_bed = parse_shard_bed(args.shard_bed)
        vcf_paths = find_overlapping_shards(regions, shard_bed, vcf_root=args.vcf_root, vcf_template=args.vcf_template)
        participants = scan_vcf_paths(vcf_paths, regions)
    else:
        if not args.vcf_dir:
            raise SystemExit("--vcf-dir is required when --shard-bed is not provided")
        participants = scan_vcfs(args.vcf_dir, regions)
    write_outputs(args.out_prefix, participants, participant_rows, args.group_col)


if __name__ == "__main__":
    main()
