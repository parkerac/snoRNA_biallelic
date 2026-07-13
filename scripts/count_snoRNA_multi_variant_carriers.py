#!/usr/bin/env python3
"""Count snoRNA-overlapping variants gene by gene from AGGV3 shard BEDs."""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


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
    shards = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 7:
                continue
            shards.append(
                {
                    "chrom": parts[0],
                    "start": int(parts[1]) + 1,
                    "end": int(parts[2]),
                    "shard": parts[4] if len(parts) > 4 else "",
                    "subshard": parts[5] if len(parts) > 5 else "",
                    "vcf_path": parts[6] if len(parts) > 6 else "",
                }
            )
    return sorted(shards, key=lambda s: (s["chrom"], s["start"], s["end"], s["shard"], s["subshard"]))


def gene_overlaps_interval(gene, start, end):
    return not (gene["end"] < start or gene["start"] > end)


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


def safe_tabix_regions(regions):
    return [f"{chrom}:{start}-{end}" for chrom, start, end in regions]


def iter_vcf_lines(vcf_path, regions, region_access):
    regions = list(dict.fromkeys(safe_tabix_regions(regions)))
    use_tabix = region_access in {"auto", "tabix"} and shutil.which("tabix")
    if use_tabix:
        cmd = ["tabix", "-h", vcf_path, *regions]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        try:
            for line in proc.stdout:
                yield line.rstrip("\n")
        finally:
            if proc.stdout:
                proc.stdout.close()
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f"tabix failed for {vcf_path} with exit code {rc}")
        return

    with open_text(vcf_path) as fh:
        for line in fh:
            yield line.rstrip("\n")


def make_gene_detail_path(detail_dir, gene_index, gene):
    return os.path.join(
        detail_dir,
        f"{gene_index:04d}_{safe_name(gene['gene_name'])}_{safe_name(gene['gene_id'])}.tsv",
    )


def build_shard_tasks(genes, shards, vcf_root, vcf_template):
    by_chrom = defaultdict(list)
    for gene in genes:
        by_chrom[gene["chrom"]].append(gene)

    tasks = []
    for shard_index, shard in enumerate(shards, start=1):
        vcf_path = os.path.join(vcf_root, vcf_template.format(shard=shard["shard"], subshard=shard["subshard"]))
        shard_genes = [gene for gene in by_chrom.get(shard["chrom"], []) if gene_overlaps_interval(gene, shard["start"], shard["end"])]
        if not shard_genes:
            continue
        tasks.append(
            {
                "shard_index": shard_index,
                "shard": shard,
                "vcf_path": vcf_path,
                "genes": shard_genes,
                "regions": [(gene["chrom"], gene["start"], gene["end"]) for gene in shard_genes],
            }
        )
    return tasks


def process_shard(task):
    shard = task["shard"]
    genes = task["genes"]
    vcf_path = task["vcf_path"]
    region_access = task["region_access"]
    gene_by_index = {gene["gene_index"]: gene for gene in genes}
    seen_variants = set()
    samples = []
    rows_by_gene = defaultdict(list)
    gene_variant_ids = defaultdict(set)
    gene_participants = defaultdict(set)

    for line in iter_vcf_lines(vcf_path, task["regions"], region_access):
        if not line or line.startswith("##"):
            continue
        if line.startswith("#CHROM"):
            samples = line.split("\t")[9:]
            continue
        parts = line.split("\t")
        if len(parts) < 10 or not samples or parts[0] != shard["chrom"]:
            continue
        pos = int(parts[1])
        overlapping_genes = [gene for gene in genes if gene["start"] <= pos <= gene["end"]]
        if not overlapping_genes:
            continue
        fmt = parts[8].split(":")
        if "GT" not in fmt:
            continue
        gt_idx = fmt.index("GT")
        variant_id = f"{parts[0]}:{pos}:{parts[3]}:{parts[4]}"
        if variant_id in seen_variants:
            continue
        carrier_rows = []
        for sample, sample_field in zip(samples, parts[9:]):
            fields = sample_field.split(":")
            if len(fields) <= gt_idx:
                continue
            gt = fields[gt_idx]
            if gt in {".", "./.", ".|."} or not is_carrier(gt):
                continue
            carrier_rows.append((sample, gt))
        if not carrier_rows:
            seen_variants.add(variant_id)
            continue
        for gene in overlapping_genes:
            gene_key = gene["gene_index"]
            gene_variant_ids[gene_key].add(variant_id)
            for sample, gt in carrier_rows:
                gene_participants[gene_key].add(sample)
                rows_by_gene[gene_key].append(
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
                        "genotype": gt,
                    }
                )
        seen_variants.add(variant_id)

    gene_stats = {}
    for gene_index, gene in gene_by_index.items():
        gene_stats[gene_index] = {
            "gene_index": gene_index,
            "gene_name": gene["gene_name"],
            "gene_id": gene["gene_id"],
            "rna_class": gene["rna_class"],
            "chrom": gene["chrom"],
            "start": gene["start"],
            "end": gene["end"],
            "n_vcfs": 1,
            "n_variants": len(gene_variant_ids.get(gene_index, set())),
            "n_participants": len(gene_participants.get(gene_index, set())),
        }

    return {
        "shard_index": task["shard_index"],
        "shard": shard,
        "vcf_path": vcf_path,
        "rows_by_gene": rows_by_gene,
        "gene_stats": gene_stats,
        "gene_participants": gene_participants,
        "gene_variant_ids": gene_variant_ids,
    }


def write_gene_summary(out_prefix, gene_agg, gene_paths):
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
        for gene_index in sorted(gene_agg):
            row = gene_agg[gene_index]
            writer.writerow(
                {
                    "gene_index": gene_index,
                    "gene_name": row["gene_name"],
                    "gene_id": row["gene_id"],
                    "rna_class": row["rna_class"],
                    "chrom": row["chrom"],
                    "start": row["start"],
                    "end": row["end"],
                    "n_vcfs": row["n_vcfs"],
                    "n_variants": len(row["variant_ids"]),
                    "n_participants": len(row["participants"]),
                    "detail_path": gene_paths[gene_index],
                }
            )


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
    parser.add_argument("--region-access", choices=["auto", "tabix", "scan"], default="auto")
    args = parser.parse_args()

    print(f"Loading genes from {args.gtf}", flush=True)
    genes = parse_gtf(args.gtf, set(args.feature_type))
    print(f"Loaded {len(genes)} genes", flush=True)
    print(f"Loading shard BED from {args.shard_bed}", flush=True)
    shards = parse_shard_bed(args.shard_bed)
    participant_rows = read_participants(args.participant_tsv, args.participant_id_col, args.group_col)

    out_prefix_dir = os.path.dirname(args.out_prefix)
    if out_prefix_dir:
        os.makedirs(out_prefix_dir, exist_ok=True)
    detail_dir = f"{args.out_prefix}.genes"
    os.makedirs(detail_dir, exist_ok=True)

    gene_paths = {}
    gene_agg = {}
    for gene_index, gene in enumerate(genes, start=1):
        gene["gene_index"] = gene_index
        gene_paths[gene_index] = make_gene_detail_path(detail_dir, gene_index, gene)
        gene_agg[gene_index] = {
            "gene_name": gene["gene_name"],
            "gene_id": gene["gene_id"],
            "rna_class": gene["rna_class"],
            "chrom": gene["chrom"],
            "start": gene["start"],
            "end": gene["end"],
            "n_vcfs": 0,
            "variant_ids": set(),
            "participants": set(),
        }
        with open(gene_paths[gene_index], "w", newline="") as fh:
            writer = csv.DictWriter(
                fh,
                delimiter="\t",
                fieldnames=["gene_name", "gene_id", "rna_class", "chrom", "pos", "ref", "alt", "variant_id", "participant_id", "genotype"],
            )
            writer.writeheader()

    shard_tasks = build_shard_tasks(genes, shards, args.vcf_root, args.vcf_template)
    for task in shard_tasks:
        task["region_access"] = args.region_access
        print(
            f"[queued shard {task['shard_index']}] {task['shard']['shard']}/{task['shard']['subshard']}: "
            f"{len(task['genes'])} genes, {task['vcf_path']}",
            flush=True,
        )

    participants = defaultdict(lambda: {"variants": set(), "genes": set()})
    workers = max(1, min(args.cpus, len(shard_tasks)))
    print(f"Processing {len(shard_tasks)} shards with {workers} workers", flush=True)

    with ThreadPoolExecutor(max_workers=workers) as pool:
        futures = [pool.submit(process_shard, task) for task in shard_tasks]
        for done, future in enumerate(as_completed(futures), start=1):
            result = future.result()
            print(
                f"[done {done}/{len(shard_tasks)}] shard {result['shard']['shard']}/{result['shard']['subshard']} "
                f"from {os.path.basename(result['vcf_path'])}",
                flush=True,
            )
            for gene_index, rows in result["rows_by_gene"].items():
                with open(gene_paths[gene_index], "a", newline="") as fh:
                    writer = csv.DictWriter(
                        fh,
                        delimiter="\t",
                        fieldnames=["gene_name", "gene_id", "rna_class", "chrom", "pos", "ref", "alt", "variant_id", "participant_id", "genotype"],
                    )
                    for row in rows:
                        writer.writerow(row)
            for gene_index, stats in result["gene_stats"].items():
                gene_agg[gene_index]["n_vcfs"] += stats["n_vcfs"]
                gene_agg[gene_index]["variant_ids"].update(result["gene_variant_ids"].get(gene_index, set()))
                gene_agg[gene_index]["participants"].update(result["gene_participants"].get(gene_index, set()))
            for gene_index, rows in result["rows_by_gene"].items():
                gene = gene_agg[gene_index]
                print(
                    f"  gene {gene['gene_name']} ({gene['gene_id']}): +{len(rows)} rows, "
                    f"{len(result['gene_variant_ids'].get(gene_index, set()))} variants in this shard",
                    flush=True,
                )
            for gene_index, rows in result["rows_by_gene"].items():
                gene_name = gene_agg[gene_index]["gene_name"]
                gene_id = gene_agg[gene_index]["gene_id"]
                for row in rows:
                    pid = row["participant_id"]
                    participants[pid]["variants"].add(row["variant_id"])
                    participants[pid]["genes"].add((gene_name, gene_id))

    write_gene_summary(args.out_prefix, gene_agg, gene_paths)
    write_participant_outputs(args.out_prefix, participants, participant_rows, args.group_col)


if __name__ == "__main__":
    main()
