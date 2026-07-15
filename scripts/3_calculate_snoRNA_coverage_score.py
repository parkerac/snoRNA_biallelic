#!/usr/bin/env python3
"""Estimate snoRNA coverage from AGGV3 site-QC VCFs using cohort MEDIAN_DP."""

import argparse
import csv
import gzip
import os
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, as_completed

try:
    from cyvcf2 import VCF as CyVCF
except Exception:  # pragma: no cover
    CyVCF = None


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def load_genes(path):
    genes = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"gene_index", "gene_name", "gene_id", "chrom", "start", "end"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            genes.append(
                {
                    "gene_index": int(row["gene_index"]),
                    "gene_name": row["gene_name"],
                    "gene_id": row["gene_id"],
                    "rna_class": row.get("rna_class", ""),
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                }
            )
    return sorted(genes, key=lambda g: (g["chrom"], g["start"], g["end"], g["gene_index"]))


def parse_shard_bed(path):
    shards = []
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) < 6:
                continue
            shards.append(
                {
                    "chrom": parts[0],
                    "start": int(parts[1]) + 1,
                    "end": int(parts[2]),
                    "shard": parts[4] if len(parts) > 4 else "",
                    "subshard": parts[5] if len(parts) > 5 else "",
                }
            )
    return sorted(shards, key=lambda s: (s["chrom"], s["start"], s["end"], s["shard"], s["subshard"]))


def gene_overlaps_interval(gene, start, end):
    return not (gene["end"] < start or gene["start"] > end)


def merge_gene_windows(genes):
    if not genes:
        return []
    windows = []
    current = {"chrom": genes[0]["chrom"], "start": genes[0]["start"], "end": genes[0]["end"], "genes": [genes[0]]}
    for gene in genes[1:]:
        if gene["start"] <= current["end"] + 1:
            current["end"] = max(current["end"], gene["end"])
            current["genes"].append(gene)
        else:
            windows.append(current)
            current = {"chrom": gene["chrom"], "start": gene["start"], "end": gene["end"], "genes": [gene]}
    windows.append(current)
    return windows


def has_index(vcf_path):
    return any(os.path.exists(vcf_path + suffix) for suffix in (".tbi", ".csi"))


def choose_region_mode(vcf_path, requested):
    if requested == "auto":
        if CyVCF is not None and has_index(vcf_path):
            return "cyvcf2"
        if shutil.which("tabix") and has_index(vcf_path):
            return "tabix"
        return "scan"
    if requested == "cyvcf2" and CyVCF is not None and has_index(vcf_path):
        return "cyvcf2"
    if requested == "tabix" and shutil.which("tabix") and has_index(vcf_path):
        return "tabix"
    return "scan"


def parse_info_text(info_text):
    info = {}
    for item in info_text.split(";"):
        if not item:
            continue
        if "=" in item:
            key, value = item.split("=", 1)
            info[key] = value
        else:
            info[item] = True
    return info


def scalar_float(value):
    if value is None or value == "":
        return None
    if hasattr(value, "__iter__") and not isinstance(value, (str, bytes)):
        for item in value:
            if item not in (None, ""):
                value = item
                break
        else:
            return None
    try:
        return float(value)
    except Exception:
        return None


def iter_records(vcf_path, regions, region_mode):
    if not regions:
        return

    mode = choose_region_mode(vcf_path, region_mode)
    regions = list(regions)

    if mode == "cyvcf2":
        vcf = CyVCF(vcf_path)
        for chrom, start, end in regions:
            for variant in vcf(f"{chrom}:{start}-{end}"):
                yield {
                    "chrom": variant.CHROM,
                    "pos": int(variant.POS),
                    "info": {"MEDIAN_DP": variant.INFO.get("MEDIAN_DP")},
                }
        return

    if mode == "tabix":
        cmd = ["tabix", "-h", vcf_path, *[f"{chrom}:{start}-{end}" for chrom, start, end in regions]]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        try:
            for line in proc.stdout:
                line = line.rstrip("\n")
                if not line or line.startswith("#"):
                    continue
                parts = line.split("\t")
                yield {
                    "chrom": parts[0],
                    "pos": int(parts[1]),
                    "info": parse_info_text(parts[7]),
                }
        finally:
            if proc.stdout:
                proc.stdout.close()
            rc = proc.wait()
            if rc != 0:
                raise RuntimeError(f"tabix failed for {vcf_path} with exit code {rc}")
        return

    with open_text(vcf_path) as fh:
        windows = list(regions)
        current = 0
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom = parts[0]
            pos = int(parts[1])
            while current < len(windows) and (chrom > windows[current][0] or (chrom == windows[current][0] and pos > windows[current][2])):
                current += 1
            if current >= len(windows):
                break
            wchrom, wstart, wend = windows[current]
            if chrom != wchrom or pos < wstart or pos > wend:
                continue
            yield {
                "chrom": chrom,
                "pos": pos,
                "info": parse_info_text(parts[7]),
            }


def build_shard_tasks(genes, shards, vcf_root, vcf_template):
    by_chrom = defaultdict(list)
    for gene in genes:
        by_chrom[gene["chrom"]].append(gene)

    tasks = []
    for shard_index, shard in enumerate(shards, start=1):
        shard_genes = [gene for gene in by_chrom.get(shard["chrom"], []) if gene_overlaps_interval(gene, shard["start"], shard["end"])]
        if not shard_genes:
            continue
        tasks.append(
            {
                "shard_index": shard_index,
                "shard": shard,
                "vcf_path": os.path.join(vcf_root, vcf_template.format(shard=shard["shard"], subshard=shard["subshard"])),
                "windows": merge_gene_windows(shard_genes),
                "gene_indices": [gene["gene_index"] for gene in shard_genes],
            }
        )
    return tasks


def process_shard(task):
    shard = task["shard"]
    region_mode = task["region_mode"]
    vcf_path = task["vcf_path"]
    windows = task["windows"]
    gene_stats = {gene["gene_index"]: {"n_sites": 0, "dp_sum": 0.0} for window in windows for gene in window["genes"]}
    records = 0

    if not os.path.exists(vcf_path):
        return {
            "shard_index": task["shard_index"],
            "total_shards": task["total_shards"],
            "shard": shard,
            "vcf_path": vcf_path,
            "gene_stats": gene_stats,
            "records": 0,
            "missing": True,
        }

    if region_mode == "cyvcf2":
        vcf = CyVCF(vcf_path)
        for window in windows:
            for variant in vcf(f"{window['chrom']}:{window['start']}-{window['end']}"):
                dp = scalar_float(variant.INFO.get("MEDIAN_DP"))
                if dp is None:
                    continue
                records += 1
                pos = int(variant.POS)
                for gene in window["genes"]:
                    if gene["start"] <= pos <= gene["end"]:
                        gene_stats[gene["gene_index"]]["n_sites"] += 1
                        gene_stats[gene["gene_index"]]["dp_sum"] += dp
        return {
            "shard_index": task["shard_index"],
            "total_shards": task["total_shards"],
            "shard": shard,
            "vcf_path": vcf_path,
            "gene_stats": gene_stats,
            "records": records,
            "missing": False,
        }

    regions = [(window["chrom"], window["start"], window["end"]) for window in windows]
    window_index = 0
    for record in iter_records(vcf_path, regions, region_mode):
        dp = scalar_float(record["info"].get("MEDIAN_DP"))
        if dp is None:
            continue
        records += 1
        pos = record["pos"]
        while window_index < len(windows) and pos > windows[window_index]["end"]:
            window_index += 1
        if window_index >= len(windows):
            break
        window = windows[window_index]
        if pos < window["start"] or pos > window["end"]:
            continue
        for gene in window["genes"]:
            if gene["start"] <= pos <= gene["end"]:
                gene_stats[gene["gene_index"]]["n_sites"] += 1
                gene_stats[gene["gene_index"]]["dp_sum"] += dp

    return {
        "shard_index": task["shard_index"],
        "total_shards": task["total_shards"],
        "shard": shard,
        "vcf_path": vcf_path,
        "gene_stats": gene_stats,
        "records": records,
        "missing": False,
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-summary", required=True, help="gene_summary.tsv from script 1")
    parser.add_argument("--shard-bed", required=True, help="biallelic_shards.bed")
    parser.add_argument("--site-qc-root", required=True, help="Directory containing shard/subshard site-QC VCFs")
    parser.add_argument(
        "--vcf-template",
        default="shard-{shard}/subshard-{subshard}/postproc/site_qc.vcf.gz",
        help="Relative site-QC VCF path template under --site-qc-root",
    )
    parser.add_argument("--out", required=True, help="Output TSV path")
    parser.add_argument("--cpus", type=int, default=16)
    parser.add_argument("--region-access", choices=["auto", "cyvcf2", "tabix", "scan"], default="auto")
    args = parser.parse_args()

    print(f"Loading genes from {args.gene_summary}", flush=True)
    genes = load_genes(args.gene_summary)
    print(f"Loaded {len(genes)} genes", flush=True)
    print(f"Loading shard BED from {args.shard_bed}", flush=True)
    shards = parse_shard_bed(args.shard_bed)

    tasks = build_shard_tasks(genes, shards, args.site_qc_root, args.vcf_template)
    total_tasks = len(tasks)
    for task in tasks:
        task["region_mode"] = args.region_access
        task["total_shards"] = total_tasks
        print(
            f"[queued task {task['shard_index']}/{total_tasks}] {task['shard']['shard']}/{task['shard']['subshard']}: "
            f"{len(task['windows'])} merged windows, {task['vcf_path']}",
            flush=True,
        )

    workers = max(1, min(args.cpus, total_tasks or 1))
    print(f"Processing {total_tasks} shards with {workers} workers", flush=True)

    gene_stats = {gene["gene_index"]: {"n_sites": 0, "dp_sum": 0.0} for gene in genes}
    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(process_shard, task) for task in tasks]
            for future in as_completed(futures):
                result = future.result()
                prefix = "missing" if result.get("missing") else "done"
                print(
                    f"[{prefix} task {result['shard_index']}/{result['total_shards']}] {result['shard']['shard']}/{result['shard']['subshard']} from {os.path.basename(result['vcf_path'])}: "
                    f"{result['records']} MEDIAN_DP records",
                    flush=True,
                )
                for gene_index, stats in result["gene_stats"].items():
                    gene_stats[gene_index]["n_sites"] += stats["n_sites"]
                    gene_stats[gene_index]["dp_sum"] += stats["dp_sum"]
    except (PermissionError, OSError) as exc:
        print(f"Process pool unavailable ({exc}); falling back to threads", flush=True)
        from concurrent.futures import ThreadPoolExecutor

        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(process_shard, task) for task in tasks]
            for future in as_completed(futures):
                result = future.result()
                prefix = "missing" if result.get("missing") else "done"
                print(
                    f"[{prefix} task {result['shard_index']}/{result['total_shards']}] {result['shard']['shard']}/{result['shard']['subshard']} from {os.path.basename(result['vcf_path'])}: "
                    f"{result['records']} MEDIAN_DP records",
                    flush=True,
                )
                for gene_index, stats in result["gene_stats"].items():
                    gene_stats[gene_index]["n_sites"] += stats["n_sites"]
                    gene_stats[gene_index]["dp_sum"] += stats["dp_sum"]

    print(f"Writing coverage score table to {args.out}", flush=True)
    with open(args.out, "w", newline="") as fh:
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
                "gene_length_bp",
                "n_site_qc_records",
                "coverage_score",
            ],
        )
        writer.writeheader()
        for gene in genes:
            stats = gene_stats[gene["gene_index"]]
            count = stats["n_sites"]
            writer.writerow(
                {
                    "gene_index": gene["gene_index"],
                    "gene_name": gene["gene_name"],
                    "gene_id": gene["gene_id"],
                    "rna_class": gene["rna_class"],
                    "chrom": gene["chrom"],
                    "start": gene["start"],
                    "end": gene["end"],
                    "gene_length_bp": gene["end"] - gene["start"] + 1,
                    "n_site_qc_records": count,
                    "coverage_score": stats["dp_sum"] / count if count else 0.0,
                }
            )


if __name__ == "__main__":
    main()
