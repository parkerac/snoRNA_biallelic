#!/usr/bin/env python3
"""Calculate rare variant density per snoRNA from the gene summary table."""

import argparse
import csv


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--gene-summary", required=True, help="gene_summary.tsv from script 1")
    parser.add_argument("--out", required=True, help="Output TSV path")
    args = parser.parse_args()

    with open(args.gene_summary, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"gene_index", "gene_name", "gene_id", "chrom", "start", "end", "n_rare_variants"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{args.gene_summary} is missing columns: {', '.join(sorted(missing))}")
        rows = []
        for row in reader:
            start = int(row["start"])
            end = int(row["end"])
            length_bp = end - start + 1
            rare_variants = int(row["n_rare_variants"])
            rows.append(
                {
                    "gene_index": row["gene_index"],
                    "gene_name": row["gene_name"],
                    "gene_id": row["gene_id"],
                    "chrom": row["chrom"],
                    "start": start,
                    "end": end,
                    "gene_length_bp": length_bp,
                    "n_rare_variants": rare_variants,
                    "rare_variants_per_kb": rare_variants / (length_bp / 1000) if length_bp else 0.0,
                    "rare_variants_per_mb": rare_variants / (length_bp / 1_000_000) if length_bp else 0.0,
                    "rare_variants_per_bp": rare_variants / length_bp if length_bp else 0.0,
                }
            )

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(
            fh,
            delimiter="\t",
            fieldnames=[
                "gene_index",
                "gene_name",
                "gene_id",
                "chrom",
                "start",
                "end",
                "gene_length_bp",
                "n_rare_variants",
                "rare_variants_per_kb",
                "rare_variants_per_mb",
                "rare_variants_per_bp",
            ],
        )
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
