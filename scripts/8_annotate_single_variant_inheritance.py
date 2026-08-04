#!/usr/bin/env python3
"""Annotate one or more variants as inherited, de novo, or uncertain using family VCFs."""

import argparse
import csv
import os
from concurrent.futures import ThreadPoolExecutor, as_completed


MISSING_VALUES = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN"}
DEFAULT_THREADS = 8
DEFAULT_ORIGIN_WINDOW = 5000
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


def variant_status_for_record(rec, variant, sample):
    alleles = gt_alleles(rec, sample)
    if alleles is None:
        return "ambiguous"
    alt_indices = real_alt_indices(rec)
    if rec.pos == variant[1] and rec.ref.upper() == variant[2] and variant[3] in alt_indices:
        return "has_alt" if alt_indices[variant[3]] in alleles else "no_alt"
    if rec.start <= variant[1] - 1 < rec.stop and alleles == {0}:
        return "no_alt"
    return "ambiguous" if any(a > 0 for a in alleles) else "no_alt"


def statuses_for_vcf(vcf_path, sample, variants):
    if not vcf_path:
        return ["missing"] * len(variants)
    if not os.path.exists(vcf_path):
        return ["missing_file"] * len(variants)
    with pysam.VariantFile(vcf_path) as vcf:
        sample = sample or next(iter(vcf.header.samples), None)
        if not sample or sample not in vcf.header.samples:
            return ["missing_sample"] * len(variants)
        statuses = ["no_record"] * len(variants)
        by_chrom = {}
        for i, variant in enumerate(variants):
            by_chrom.setdefault(variant[0], []).append((i, variant))
        for chrom, chrom_variants in by_chrom.items():
            try:
                records = list(vcf.fetch(resolve_contig(vcf, chrom), max(0, min(v[1] for _, v in chrom_variants) - 1), max(v[1] + len(v[2]) for _, v in chrom_variants)))
            except Exception:
                continue
            for i, variant in chrom_variants:
                saw_record = False
                status = "no_record"
                for rec in records:
                    if not (rec.start <= variant[1] - 1 < rec.stop):
                        continue
                    saw_record = True
                    status = variant_status_for_record(rec, variant, sample)
                    if status in {"has_alt", "no_alt"}:
                        break
                statuses[i] = "ambiguous" if saw_record and status not in {"has_alt", "no_alt"} else status
        return statuses


def classify_inheritance(mother_status, father_status):
    if mother_status == "has_alt" or father_status == "has_alt":
        return "inherited", "parent_carrier"
    if mother_status == "no_alt" and father_status == "no_alt":
        return "de_novo", "neither_parent_carries_variant"
    return "uncertain", "parental_evidence_inconclusive"


def classify_origin_hint(index, variants, mother_statuses, father_statuses, window):
    if mother_statuses[index] != "no_alt" or father_statuses[index] != "no_alt":
        return ""
    chrom, pos, _, _ = variants[index]
    maternal = paternal = 0
    for j, other in enumerate(variants):
        if j == index or other[0] != chrom or abs(other[1] - pos) > window:
            continue
        if mother_statuses[j] == "has_alt" and father_statuses[j] == "no_alt":
            maternal += 1
        elif father_statuses[j] == "has_alt" and mother_statuses[j] == "no_alt":
            paternal += 1
    if maternal and not paternal:
        return "maternal"
    if paternal and not maternal:
        return "paternal"
    return "uncertain"


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


def unique_variants(variants):
    seen = set()
    out = []
    for variant in variants:
        if variant not in seen:
            seen.add(variant)
            out.append(variant)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", required=True, help="TSV with one or more variants per row")
    parser.add_argument("--variant-column", default="variant_id", help="Variant column in chrom:pos:ref:alt format, optionally semicolon-separated")
    parser.add_argument("--sample-column", default="sample", help="Sample column to preserve in the output")
    parser.add_argument("--mother-vcf-column", default="mother_vcf", help="Mother VCF column")
    parser.add_argument("--father-vcf-column", default="father_vcf", help="Father VCF column")
    parser.add_argument("--mother-sample-column", default="mother_sample", help="Mother sample column")
    parser.add_argument("--father-sample-column", default="father_sample", help="Father sample column")
    parser.add_argument("--out", required=True, help="Output TSV with inheritance annotation")
    parser.add_argument("--annotation-column", default="inheritance_annotation", help="Name of the output annotation column")
    parser.add_argument("--detail-column", default="inheritance_detail", help="Name of the output detail column")
    parser.add_argument("--origin-column", default="parental_origin", help="Name of the output maternal/paternal origin hint column")
    parser.add_argument("--origin-window", type=int, default=DEFAULT_ORIGIN_WINDOW, help="Nearby base window used to infer de novo parent-of-origin from other informative variants")
    parser.add_argument("--threads", type=int, default=DEFAULT_THREADS, help="Number of worker threads")
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
        rows = list(reader)
        for row in rows:
            if args.sample_column in row and row_value(row, args.sample_column):
                row["sample"] = row_value(row, args.sample_column)
        total = len(rows)

    groups = {}
    order = []
    for row_index, row in enumerate(rows):
        key = (
            row_value(row, args.sample_column),
            row_value(row, args.mother_vcf_column),
            row_value(row, args.mother_sample_column),
            row_value(row, args.father_vcf_column),
            row_value(row, args.father_sample_column),
        )
        if key not in groups:
            groups[key] = []
            order.append(key)
        groups[key].append((row_index, row))

    def process_group(item):
        group_index, key = item
        grouped_rows = groups[key]
        row_variants = [(row_index, row, variant_from_row(row, args.variant_column)) for row_index, row in grouped_rows]
        all_variants = unique_variants([variant for _, _, variants in row_variants for variant in variants])
        mother_statuses = statuses_for_vcf(key[1], key[2], all_variants)
        father_statuses = statuses_for_vcf(key[3], key[4], all_variants)
        status_lookup = {variant: i for i, variant in enumerate(all_variants)}
        out = []
        for row_index, row, variants in row_variants:
            row_mother_statuses = [mother_statuses[status_lookup[variant]] for variant in variants]
            row_father_statuses = [father_statuses[status_lookup[variant]] for variant in variants]
            annotations = []
            details = []
            origin_hints = []
            for i, _variant in enumerate(variants):
                annotation, detail = classify_inheritance(row_mother_statuses[i], row_father_statuses[i])
                annotations.append(annotation)
                if annotation == "de_novo":
                    hint = classify_origin_hint(i, variants, row_mother_statuses, row_father_statuses, args.origin_window)
                    origin_hints.append(hint)
                    details.append(hint if hint else detail)
                else:
                    origin_hints.append("")
                    details.append(detail)
            row[args.annotation_column] = ";".join(annotations)
            row[args.detail_column] = ";".join(details)
            row[args.origin_column] = ";".join(origin_hints)
            row["mother_status"] = ";".join(row_mother_statuses)
            row["father_status"] = ";".join(row_father_statuses)
            out.append((row_index, row))
        return group_index, out

    rows_out = [None] * total
    with ThreadPoolExecutor(max_workers=max(1, args.threads)) as pool:
        futures = {pool.submit(process_group, (i, key)): i for i, key in enumerate(order, start=1)}
        done = 0
        for future in as_completed(futures):
            _, group_rows = future.result()
            for row_index, row in group_rows:
                rows_out[row_index] = row
            done += 1
            print(f"Finished sample group {done}/{len(order)}", flush=True)

    fieldnames = list(rows_out[0].keys()) if rows_out else list(reader.fieldnames or [])
    for extra in (args.annotation_column, args.detail_column, args.origin_column, "mother_status", "father_status"):
        if extra not in fieldnames:
            fieldnames.append(extra)

    with open(args.out, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows_out)


if __name__ == "__main__":
    main()
