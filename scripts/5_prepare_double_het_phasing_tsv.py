#!/usr/bin/env python3
"""Prepare double-het non-coding RNA variant pairs for phase_nearby_variants.py."""

import argparse
import csv
import itertools


def parse_list(value):
    return [item for item in (value or "").split(";") if item]


def is_heterozygous_alt(gt):
    if not gt:
        return False
    alleles = gt.replace("|", "/").split("/")
    return len(alleles) == 2 and alleles.count("0") == 1 and any(allele not in {".", "0"} for allele in alleles)


def load_filepath_details(path):
    details = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {
            "sample",
            "bam",
            "vcf",
            "father_bam",
            "mother_bam",
            "father_vcf",
            "mother_vcf",
            "father_sample",
            "mother_sample",
        }
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{path} is missing columns: {', '.join(sorted(missing))}")
        for row in reader:
            details[row["sample"]] = row
    return details


def parse_variant(value):
    chrom, pos, ref, alt = value.split(":")
    return chrom, int(pos), ref, alt


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--double-het-tsv", required=True, help="TSV from script 4")
    parser.add_argument("--filepath-details", required=True, help="TSV with sample and family filepath columns")
    parser.add_argument("--reference-path", required=True, help="Reference FASTA path to add to every row")
    parser.add_argument("--out", required=True, help="Output TSV for phase_nearby_variants.py")
    args = parser.parse_args()

    filepath_details = load_filepath_details(args.filepath_details)
    rows = []
    missing_samples = set()

    with open(args.double_het_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        required = {"platekey", "gene_name", "gene_id", "variants", "genotypes"}
        missing = required - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"{args.double_het_tsv} is missing columns: {', '.join(sorted(missing))}")

        for row in reader:
            sample = row["platekey"]
            if sample not in filepath_details:
                missing_samples.add(sample)
                continue

            variants = parse_list(row.get("variants"))
            genotypes = parse_list(row.get("genotypes"))
            if len(variants) != len(genotypes):
                raise ValueError(f"{args.double_het_tsv} has mismatched variants/genotypes for platekey {sample}")

            het_variants = [(variant, genotype) for variant, genotype in zip(variants, genotypes) if is_heterozygous_alt(genotype)]
            if len(het_variants) < 2:
                continue

            for (variant1, _), (variant2, _) in itertools.combinations(het_variants, 2):
                chrom1, pos1, ref1, alt1 = parse_variant(variant1)
                chrom2, pos2, ref2, alt2 = parse_variant(variant2)
                if chrom1 != chrom2:
                    continue
                details = filepath_details[sample]
                rows.append(
                    {
                        "sample": details.get("sample", sample),
                        "platekey": sample,
                        "gene_name": row["gene_name"],
                        "gene_id": row["gene_id"],
                        "chrom": chrom1,
                        "pos1": pos1,
                        "ref1": ref1,
                        "alt1": alt1,
                        "pos2": pos2,
                        "ref2": ref2,
                        "alt2": alt2,
                        "reference": args.reference_path,
                        "bam": details.get("bam", ""),
                        "father_bam": details.get("father_bam", ""),
                        "mother_bam": details.get("mother_bam", ""),
                        "vcf": details.get("vcf", ""),
                        "father_vcf": details.get("father_vcf", ""),
                        "mother_vcf": details.get("mother_vcf", ""),
                        "father_sample": details.get("father_sample", ""),
                        "mother_sample": details.get("mother_sample", ""),
                    }
                )

    if missing_samples:
        print(f"Warning: {len(missing_samples)} samples were not found in {args.filepath_details}", flush=True)

    fieldnames = [
        "sample",
        "platekey",
        "gene_name",
        "gene_id",
        "chrom",
        "pos1",
        "ref1",
        "alt1",
        "pos2",
        "ref2",
        "alt2",
        "reference",
        "bam",
        "father_bam",
        "mother_bam",
        "vcf",
        "father_vcf",
        "mother_vcf",
        "father_sample",
        "mother_sample",
    ]
    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
