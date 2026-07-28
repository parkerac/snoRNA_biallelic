#!/usr/bin/env python3
"""List participants with >=2 rare heterozygous variants in the same snoRNA gene."""

import argparse
import csv
import os
from collections import defaultdict


RARE_AF_THRESHOLD = 0.005


def is_heterozygous_alt(gt):
    if not gt:
        return False
    alleles = gt.replace("|", "/").split("/")
    return len(alleles) == 2 and alleles.count("0") == 1 and any(allele not in {".", "0"} for allele in alleles)


def parse_af(value):
    if value is None or value == "":
        return []
    out = []
    for item in str(value).split(","):
        try:
            out.append(float(item))
        except Exception:
            pass
    return out


def read_gene_tsv(path, threshold):
    hits = defaultdict(lambda: {"variants": {}, "gene_path": path})
    het_rows = []
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"platekey", "gene_name", "gene_id", "variant_id", "genotype", "AF"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            af_values = parse_af(row.get("AF"))
            if not af_values or min(af_values) >= threshold or not is_heterozygous_alt(row.get("genotype")):
                continue
            het_rows.append(
                {
                    "platekey": row["platekey"],
                    "gene_name": row["gene_name"],
                    "gene_id": row["gene_id"],
                    "variant_id": row["variant_id"],
                    "genotype": row["genotype"],
                    "AF": row.get("AF", ""),
                    "AC": row.get("AC", ""),
                    "AN": row.get("AN", ""),
                    "gene_tsv": path,
                }
            )
            key = (row["platekey"], row["gene_name"], row["gene_id"])
            hits[key]["variants"][row["variant_id"]] = row["genotype"]
    return hits, het_rows


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes-dir", required=True, help="Directory of per-gene TSVs from script 1")
    parser.add_argument("--out", required=True, help="Output TSV path")
    parser.add_argument("--all-het-out", help="Output TSV path for all rare heterozygous variants")
    parser.add_argument("--af-threshold", type=float, default=RARE_AF_THRESHOLD)
    args = parser.parse_args()

    aggregated = defaultdict(lambda: {"variants": {}, "gene_path": ""})
    all_het_rows = []
    for root, _, files in os.walk(args.genes_dir):
        for name in sorted(files):
            if not name.endswith(".tsv"):
                continue
            path = os.path.join(root, name)
            gene_hits, het_rows = read_gene_tsv(path, args.af_threshold)
            all_het_rows.extend(het_rows)
            for key, data in gene_hits.items():
                aggregated[key]["variants"].update(data["variants"])
                aggregated[key]["gene_path"] = path

    all_het_out = args.all_het_out or (
        args.out[:-4] + ".all_het.tsv" if args.out.endswith(".tsv") else f"{args.out}.all_het.tsv"
    )

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "platekey",
                "gene_name",
                "gene_id",
                "n_rare_heterozygous_variants",
                "variants",
                "genotypes",
                "gene_tsv",
            ],
        )
        writer.writeheader()
        for (platekey, gene_name, gene_id), data in sorted(aggregated.items()):
            if len(data["variants"]) < 2:
                continue
            items = sorted(data["variants"].items())
            writer.writerow(
                {
                    "platekey": platekey,
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "n_rare_heterozygous_variants": len(items),
                    "variants": ";".join(variant for variant, _ in items),
                    "genotypes": ";".join(genotype for _, genotype in items),
                    "gene_tsv": data["gene_path"],
                }
            )

    with open(all_het_out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "platekey",
                "gene_name",
                "gene_id",
                "variant_id",
                "genotype",
                "AF",
                "AC",
                "AN",
                "gene_tsv",
            ],
        )
        writer.writeheader()
        for row in sorted(all_het_rows, key=lambda r: (r["platekey"], r["gene_name"], r["gene_id"], r["variant_id"])):
            writer.writerow(row)


if __name__ == "__main__":
    main()
