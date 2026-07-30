#!/usr/bin/env python3
"""Annotate one or more variants as inherited, de novo, or uncertain using family VCFs."""

import argparse
import csv
import os


MISSING_VALUES = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN"}
pysam = None


def parse_variant(value):
    chrom, pos, ref, alt = value.replace(",", ":").split(":")
    return chrom, int(pos), ref.upper(), alt.upper()


def split_variant_values(value):
    return [item.strip() for item in str(value).split(";") if item.strip()]


def contig_aliases(chrom):
    aliases = [chrom]
    aliases.append(chrom[3:] if chrom.startswith("chr") else f"chr{chrom}")
    if chrom in {"M", "MT", "chrM", "chrMT"}:
        aliases.extend(["M", "MT", "chrM", "chrMT"])
    return list(dict.fromkeys(aliases))


def resolve_contig(handle, chrom):
    references = set(getattr(handle, "references", ()) or getattr(handle.header, "references", ()) or [])
    if not references and hasattr(handle, "header") and hasattr(handle.header, "contigs"):
        references = set(handle.header.contigs)
    for alias in contig_aliases(chrom):
        if not references or alias in references:
            return alias
    raise ValueError(f"Contig {chrom} was not found; tried {', '.join(contig_aliases(chrom))}")


def real_alt_indices(rec):
    return {alt.upper(): i + 1 for i, alt in enumerate(rec.alts or []) if alt and not alt.startswith("<") and alt != "*"}


def gt_alleles(rec, sample):
    if not sample:
        sample = next(iter(rec.samples), None)
    if not sample or sample not in rec.samples:
        return None
    gt = rec.samples[sample].get("GT")
    if gt is None or any(a is None for a in gt):
        return None
    return set(gt)


def classify_record(rec, variant, sample):
    alleles = gt_alleles(rec, sample)
    if alleles is None:
        return "ambiguous"
    alt_indices = real_alt_indices(rec)
    if rec.pos == variant[1] and rec.ref.upper() == variant[2] and variant[3] in alt_indices:
        return "has_alt" if alt_indices[variant[3]] in alleles else "no_alt"
    if rec.start <= variant[1] - 1 < rec.stop and alleles == {0}:
        return "no_alt"
    return "ambiguous" if any(a > 0 for a in alleles) else "no_alt"


def status_for_vcf(vcf_path, sample, variant):
    if not vcf_path:
        return "missing"
    if not os.path.exists(vcf_path):
        return "missing_file"
    with pysam.VariantFile(vcf_path) as vcf:
        sample = sample or next(iter(vcf.header.samples), None)
        if not sample or sample not in vcf.header.samples:
            return "missing_sample"
        try:
            records = vcf.fetch(resolve_contig(vcf, variant[0]), max(0, variant[1] - 1), variant[1] + len(variant[2]))
        except Exception:
            return "no_record"
        saw_record = False
        for rec in records:
            if not (rec.start <= variant[1] - 1 < rec.stop):
                continue
            saw_record = True
            status = classify_record(rec, variant, sample)
            if status in {"has_alt", "no_alt"}:
                return status
        return "ambiguous" if saw_record else "no_record"


def classify_inheritance(child_status, mother_status, father_status):
    if child_status not in {"has_alt", "assumed_present"}:
        return "uncertain", "child_not_confirmed"
    if mother_status == "has_alt" and father_status == "no_alt":
        return "inherited", "maternal"
    if father_status == "has_alt" and mother_status == "no_alt":
        return "inherited", "paternal"
    if mother_status == "no_alt" and father_status == "no_alt":
        return "de_novo", "neither_parent_has_alt"
    return "uncertain", "parental_evidence_inconclusive"


def row_value(row, key):
    value = row.get(key)
    if value is None:
        return ""
    value = value.strip()
    return "" if value.upper() in MISSING_VALUES else value


def variant_from_row(row, variant_column):
    if variant_column and row_value(row, variant_column):
        return [parse_variant(value) for value in split_variant_values(row_value(row, variant_column))]
    required = ("chrom", "pos", "ref", "alt")
    if all(row_value(row, key) for key in required):
        return [(row_value(row, "chrom"), int(row_value(row, "pos")), row_value(row, "ref").upper(), row_value(row, "alt").upper())]
    raise ValueError(f"Could not find a variant in column {variant_column or 'variant_id'} or chrom/pos/ref/alt fields")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", required=True, help="TSV with one or more variants per row")
    parser.add_argument("--variant-column", default="variant_id", help="Variant column in chrom:pos:ref:alt format, optionally semicolon-separated")
    parser.add_argument("--sample-column", default="sample", help="Proband sample column in the proband VCF")
    parser.add_argument("--vcf-column", default="vcf", help="Proband VCF column used to confirm the variant is present")
    parser.add_argument("--mother-vcf-column", default="mother_vcf", help="Mother VCF column")
    parser.add_argument("--father-vcf-column", default="father_vcf", help="Father VCF column")
    parser.add_argument("--mother-sample-column", default="mother_sample", help="Mother sample column")
    parser.add_argument("--father-sample-column", default="father_sample", help="Father sample column")
    parser.add_argument("--out", required=True, help="Output TSV with inheritance annotation")
    parser.add_argument("--annotation-column", default="inheritance_annotation", help="Name of the output annotation column")
    parser.add_argument("--detail-column", default="inheritance_detail", help="Name of the output detail column")
    args = parser.parse_args()

    try:
        import pysam as pysam_module
    except ImportError as exc:
        raise SystemExit("This script requires pysam. Install with: pip install pysam") from exc

    global pysam
    pysam = pysam_module

    with open(args.input_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        if not reader.fieldnames:
            raise SystemExit(f"{args.input_tsv} has no header")
        rows = []
        for row in reader:
            variants = variant_from_row(row, args.variant_column)
            child_statuses = []
            mother_statuses = []
            father_statuses = []
            annotations = []
            details = []
            for variant in variants:
                child_status = status_for_vcf(row_value(row, args.vcf_column), row_value(row, args.sample_column), variant)
                if not row_value(row, args.vcf_column):
                    child_status = "assumed_present"
                mother_status = status_for_vcf(row_value(row, args.mother_vcf_column), row_value(row, args.mother_sample_column), variant)
                father_status = status_for_vcf(row_value(row, args.father_vcf_column), row_value(row, args.father_sample_column), variant)
                annotation, detail = classify_inheritance(child_status, mother_status, father_status)
                child_statuses.append(child_status)
                mother_statuses.append(mother_status)
                father_statuses.append(father_status)
                annotations.append(annotation)
                details.append(detail)
            row[args.annotation_column] = ";".join(annotations)
            row[args.detail_column] = ";".join(details)
            row["child_status"] = ";".join(child_statuses)
            row["mother_status"] = ";".join(mother_statuses)
            row["father_status"] = ";".join(father_statuses)
            rows.append(row)

    fieldnames = list(rows[0].keys()) if rows else list(reader.fieldnames or [])
    for extra in (args.annotation_column, args.detail_column, "child_status", "mother_status", "father_status"):
        if extra not in fieldnames:
            fieldnames.append(extra)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)


if __name__ == "__main__":
    main()
