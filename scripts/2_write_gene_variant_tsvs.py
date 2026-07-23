#!/usr/bin/env python3
"""Write rare snoRNA variant rows to one TSV per gene."""

import argparse
import csv
import gzip
import os
import re
import shutil
import subprocess
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor, ThreadPoolExecutor, as_completed

try:
    from cyvcf2 import VCF as CyVCF
except Exception:  # pragma: no cover
    CyVCF = None


RARE_AF_THRESHOLD = 0.005
DETAIL_FIELDS = [
    "gene_name",
    "gene_id",
    "rna_class",
    "chrom",
    "pos",
    "ref",
    "alt",
    "variant_id",
    "platekey",
    "genotype",
    "AF",
    "AC",
    "AN",
]


def open_text(path):
    return gzip.open(path, "rt") if str(path).endswith(".gz") else open(path)


def safe_name(value):
    value = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(value))
    return value.strip("._") or "unknown"


def load_coverage_summary(path):
    genes = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"gene_index", "gene_name", "gene_id", "rna_class", "chrom", "start", "end", "coverage_score"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            try:
                score = float(row["coverage_score"])
            except Exception:
                continue
            genes.append(
                {
                    "gene_index": int(row["gene_index"]),
                    "gene_name": row["gene_name"],
                    "gene_id": row["gene_id"],
                    "rna_class": row["rna_class"],
                    "chrom": row["chrom"],
                    "start": int(row["start"]),
                    "end": int(row["end"]),
                    "coverage_score": score,
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
            if len(parts) < 7:
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


def is_carrier(gt):
    return any(allele not in {".", "0"} for allele in gt.replace("|", "/").split("/"))


def format_gt_alleles(genotype):
    if not genotype:
        return "./."
    a1, a2 = genotype[:2]
    phased = len(genotype) > 2 and bool(genotype[2])
    sep = "|" if phased else "/"

    def fmt(allele):
        return "." if allele is None or allele < 0 else str(allele)

    return f"{fmt(a1)}{sep}{fmt(a2)}"


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


def info_to_string(value):
    if value is None:
        return ""
    if hasattr(value, "__iter__") and not isinstance(value, (str, bytes)):
        return ",".join("" if item is None else str(item) for item in value)
    return str(value)


def info_to_floats(value):
    if value is None or value == "":
        return []
    if hasattr(value, "__iter__") and not isinstance(value, (str, bytes)):
        values = value
    else:
        values = str(value).split(",")
    out = []
    for item in values:
        try:
            out.append(float(item))
        except Exception:
            pass
    return out


def extract_info_from_cyvcf2(variant):
    info = {}
    for key in ("AF", "AC", "AN"):
        try:
            info[key] = variant.INFO.get(key)
        except Exception:
            info[key] = None
    return info


def safe_region_strings(regions):
    return [f"{chrom}:{start}-{end}" for chrom, start, end in regions]


def iter_vcf_rows(vcf_path, regions, region_mode):
    mode = choose_region_mode(vcf_path, region_mode)
    regions = safe_region_strings(regions)

    if mode == "cyvcf2":
        vcf = CyVCF(vcf_path)
        samples = list(vcf.samples)
        for region in regions:
            for variant in vcf(region):
                yield {
                    "samples": samples,
                    "chrom": variant.CHROM,
                    "pos": int(variant.POS),
                    "ref": variant.REF,
                    "alt": ",".join(variant.ALT or []),
                    "genotypes": [format_gt_alleles(gt) for gt in variant.genotypes],
                    "info": extract_info_from_cyvcf2(variant),
                }
        return

    if mode == "tabix":
        cmd = ["tabix", "-h", vcf_path, *regions]
        proc = subprocess.Popen(cmd, stdout=subprocess.PIPE, text=True)
        try:
            samples = []
            for line in proc.stdout:
                line = line.rstrip("\n")
                if not line or line.startswith("##"):
                    continue
                if line.startswith("#CHROM"):
                    samples = line.split("\t")[9:]
                    continue
                parts = line.split("\t")
                yield {
                    "samples": samples,
                    "chrom": parts[0],
                    "pos": int(parts[1]),
                    "ref": parts[3],
                    "alt": parts[4],
                    "fmt": parts[8].split(":"),
                    "sample_fields": parts[9:],
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
        samples = []
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("##"):
                continue
            if line.startswith("#CHROM"):
                samples = line.split("\t")[9:]
                continue
            parts = line.split("\t")
            yield {
                "samples": samples,
                "chrom": parts[0],
                "pos": int(parts[1]),
                "ref": parts[3],
                "alt": parts[4],
                "fmt": parts[8].split(":"),
                "sample_fields": parts[9:],
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
    vcf_path = task["vcf_path"]
    region_mode = task["region_mode"]
    rows_by_gene = defaultdict(list)

    for window in task["windows"]:
        regions = [(window["chrom"], window["start"], window["end"])]
        for record in iter_vcf_rows(vcf_path, regions, region_mode):
            chrom = record["chrom"]
            pos = record["pos"]
            ref = record["ref"]
            alt = record["alt"]
            overlapping_genes = [gene for gene in window["genes"] if gene["start"] <= pos <= gene["end"]]
            if not overlapping_genes:
                continue

            af_values = info_to_floats(record["info"].get("AF"))
            if not af_values or min(af_values) >= RARE_AF_THRESHOLD:
                continue

            af_text = info_to_string(record["info"].get("AF"))
            ac_text = info_to_string(record["info"].get("AC"))
            an_text = info_to_string(record["info"].get("AN"))
            variant_id = f"{chrom}:{pos}:{ref}:{alt}"

            if "genotypes" in record:
                carrier_rows = [
                    (sample, gt)
                    for sample, gt in zip(record["samples"], record["genotypes"])
                    if gt not in {".", "./.", ".|."} and is_carrier(gt)
                ]
            else:
                fmt = record["fmt"]
                if "GT" not in fmt:
                    continue
                gt_idx = fmt.index("GT")
                carrier_rows = []
                for sample, sample_field in zip(record["samples"], record["sample_fields"]):
                    fields = sample_field.split(":")
                    if len(fields) <= gt_idx:
                        continue
                    gt = fields[gt_idx]
                    if gt in {".", "./.", ".|."} or not is_carrier(gt):
                        continue
                    carrier_rows.append((sample, gt))

            if not carrier_rows:
                continue

            for gene in overlapping_genes:
                gene_key = gene["gene_index"]
                for sample, gt in carrier_rows:
                    rows_by_gene[gene_key].append(
                        {
                            "gene_name": gene["gene_name"],
                            "gene_id": gene["gene_id"],
                            "rna_class": gene["rna_class"],
                            "chrom": chrom,
                            "pos": pos,
                            "ref": ref,
                            "alt": alt,
                            "variant_id": variant_id,
                            "platekey": sample,
                            "genotype": gt,
                            "AF": af_text,
                            "AC": ac_text,
                            "AN": an_text,
                        }
                    )

    return {
        "shard_index": task["shard_index"],
        "shard": shard,
        "vcf_path": vcf_path,
        "rows_by_gene": rows_by_gene,
        "gene_indices": task["gene_indices"],
    }


class GeneWriter:
    def __init__(self, detail_dir, gene):
        self.path = os.path.join(detail_dir, f"{safe_name(gene['gene_name'])}_{safe_name(gene['gene_id'])}.tsv")
        self.handle = open(self.path, "w", newline="")
        self.writer = csv.DictWriter(self.handle, delimiter="\t", fieldnames=DETAIL_FIELDS)
        self.writer.writeheader()

    def write_rows(self, rows):
        for row in rows:
            self.writer.writerow(row)
        self.handle.flush()

    def close(self):
        self.handle.close()


def run_tasks(tasks, workers, gene_writers):
    total_tasks = len(tasks)

    def consume_result(result):
        count = sum(len(rows) for rows in result["rows_by_gene"].values())
        print(
            f"[done task {result['shard_index']}/{total_tasks}] {result['shard']['shard']}/{result['shard']['subshard']} from {os.path.basename(result['vcf_path'])}: {count} rare variant rows",
            flush=True,
        )
        for gene_index, rows in result["rows_by_gene"].items():
            gene_writers[gene_index].write_rows(rows)

    try:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(process_shard, task) for task in tasks]
            for future in as_completed(futures):
                consume_result(future.result())
    except (PermissionError, OSError) as exc:
        print(f"Process pool unavailable ({exc}); falling back to threads", flush=True)
        with ThreadPoolExecutor(max_workers=workers) as pool:
            futures = [pool.submit(process_shard, task) for task in tasks]
            for future in as_completed(futures):
                consume_result(future.result())


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--coverage-summary", required=True, help="Coverage summary TSV from script 1")
    parser.add_argument("--min-coverage-score", type=float, default=20.0, help="Minimum coverage_score required to query a gene")
    parser.add_argument("--shard-bed", required=True)
    parser.add_argument("--vcf-root", required=True)
    parser.add_argument("--vcf-template", default="shard-{shard}/subshard-{subshard}/postproc/vcf/dragen.vcf.gz")
    parser.add_argument("--out-prefix", required=True)
    parser.add_argument("--cpus", type=int, default=16)
    parser.add_argument("--region-access", choices=["auto", "cyvcf2", "tabix", "scan"], default="auto")
    args = parser.parse_args()

    print(f"Loading genes from {args.coverage_summary}", flush=True)
    genes = load_coverage_summary(args.coverage_summary)
    print(f"Loaded {len(genes)} genes", flush=True)
    before = len(genes)
    genes = [gene for gene in genes if gene["coverage_score"] >= args.min_coverage_score]
    print(f"Coverage filter kept {len(genes)}/{before} genes at coverage_score >= {args.min_coverage_score}", flush=True)

    print(f"Loading shard BED from {args.shard_bed}", flush=True)
    shards = parse_shard_bed(args.shard_bed)

    out_prefix_dir = os.path.dirname(args.out_prefix)
    if out_prefix_dir:
        os.makedirs(out_prefix_dir, exist_ok=True)
    detail_dir = f"{args.out_prefix}.genes"
    os.makedirs(detail_dir, exist_ok=True)

    gene_writers = {}
    for gene in genes:
        gene_writers[gene["gene_index"]] = GeneWriter(detail_dir, gene)

    tasks = build_shard_tasks(genes, shards, args.vcf_root, args.vcf_template)
    total_tasks = len(tasks)
    for task in tasks:
        task["region_mode"] = args.region_access
        print(
            f"[queued task {task['shard_index']}/{total_tasks}] {task['shard']['shard']}/{task['shard']['subshard']}: "
            f"{len(task['windows'])} merged windows, {task['vcf_path']}",
            flush=True,
        )

    workers = max(1, min(args.cpus, total_tasks or 1))
    print(f"Processing {total_tasks} shards with {workers} workers", flush=True)
    run_tasks(tasks, workers, gene_writers)

    for writer in gene_writers.values():
        writer.close()


if __name__ == "__main__":
    main()
