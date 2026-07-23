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
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"participant_id", "gene_name", "gene_id", "variant_id", "genotype", "AF"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            af_values = parse_af(row.get("AF"))
            if not af_values or min(af_values) >= threshold or not is_heterozygous_alt(row.get("genotype")):
                continue
            key = (row["participant_id"], row["gene_name"], row["gene_id"])
            hits[key]["variants"][row["variant_id"]] = row["genotype"]
    return hits


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--genes-dir", required=True, help="Directory of per-gene TSVs from script 1")
    parser.add_argument("--out", required=True, help="Output TSV path")
    parser.add_argument("--af-threshold", type=float, default=RARE_AF_THRESHOLD)
    args = parser.parse_args()

    aggregated = defaultdict(lambda: {"variants": {}, "gene_path": ""})
    for root, _, files in os.walk(args.genes_dir):
        for name in sorted(files):
            if not name.endswith(".tsv"):
                continue
            path = os.path.join(root, name)
            for key, data in read_gene_tsv(path, args.af_threshold).items():
                aggregated[key]["variants"].update(data["variants"])
                aggregated[key]["gene_path"] = path

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "participant_id",
                "gene_name",
                "gene_id",
                "n_rare_heterozygous_variants",
                "variants",
                "genotypes",
                "gene_tsv",
            ],
        )
        writer.writeheader()
        for (participant_id, gene_name, gene_id), data in sorted(aggregated.items()):
            if len(data["variants"]) < 2:
                continue
            items = sorted(data["variants"].items())
            writer.writerow(
                {
                    "participant_id": participant_id,
                    "gene_name": gene_name,
                    "gene_id": gene_id,
                    "n_rare_heterozygous_variants": len(items),
                    "variants": ";".join(variant for variant, _ in items),
                    "genotypes": ";".join(genotype for _, genotype in items),
                    "gene_tsv": data["gene_path"],
                }
            )


if __name__ == "__main__":
    main()
