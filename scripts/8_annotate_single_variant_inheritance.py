#!/usr/bin/env python3
"""Annotate one or more variants as inherited, de novo, or uncertain using family VCFs."""

import argparse
import csv
import os
from collections import defaultdict
from concurrent.futures import ThreadPoolExecutor, as_completed


MISSING_VALUES = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN"}
DEFAULT_THREADS = 8
DEFAULT_ORIGIN_WINDOW = 500
DEFAULT_MIN_MAPQ = 20
DEFAULT_MIN_BASEQ = 20
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


def is_cram(path):
    return path.lower().endswith(".cram")


def cached_reference(reference):
    return os.path.realpath(reference) if reference else None


def add_call(fragment, variant_index, allele, quality):
    old = fragment.get(variant_index)
    if old is None or (old[0] == allele and quality > old[1]):
        fragment[variant_index] = (allele, quality)
    elif old[0] != allele:
        fragment[variant_index] = None


def insertion_after(aligned_pairs, pair_index, read):
    bases = []
    for qpos, rpos in aligned_pairs[pair_index + 1 :]:
        if rpos is not None:
            break
        if qpos is not None:
            bases.append(read.query_sequence[qpos].upper())
    return "".join(bases)


def call_variant(read, variant, min_baseq):
    if read.is_unmapped or read.query_sequence is None:
        return None
    chrom, pos, ref, alt = variant
    start = pos - 1
    end = start + len(ref)
    pairs = read.get_aligned_pairs(matches_only=False)
    ref_to_qpos = {rpos: qpos for qpos, rpos in pairs if rpos is not None}
    qpos = [ref_to_qpos.get(rp) for rp in range(start, end)]
    quals = read.query_qualities or []
    observed = "".join("-" if q is None else read.query_sequence[q].upper() for q in qpos)
    baseq = min([quals[q] for q in qpos if q is not None] or [0])

    if len(ref) == len(alt):
        if None in qpos or baseq < min_baseq:
            return None
        if observed == ref:
            return 0, baseq
        if observed == alt:
            return 1, baseq
        return None

    if len(ref) > len(alt) and ref.startswith(alt):
        padded_alt = alt + "-" * (len(ref) - len(alt))
        if baseq >= min_baseq and observed == ref:
            return 0, baseq
        if baseq >= min_baseq and observed == padded_alt:
            return 1, baseq
        return None

    if len(alt) > len(ref) and alt.startswith(ref):
        if None in qpos or baseq < min_baseq or observed != ref:
            return None
        last_ref = end - 1
        pair_index = next((i for i, (q, r) in enumerate(pairs) if r == last_ref and q == qpos[-1]), None)
        inserted = insertion_after(pairs, pair_index, read) if pair_index is not None else ""
        if inserted == "":
            return 0, baseq
        if ref + inserted == alt:
            return 1, baseq
    return None


def read_fragments(bam_path, reference, variants, start, end, min_mapq, min_baseq, include_duplicates):
    fragments = defaultdict(dict)
    if not bam_path or not os.path.exists(bam_path):
        return fragments
    if is_cram(bam_path) and not reference:
        raise ValueError(f"{bam_path} is a CRAM; provide a reference")
    reference = cached_reference(reference)
    bam = pysam.AlignmentFile(bam_path, "rc" if is_cram(bam_path) else "rb", reference_filename=reference)
    try:
        fetch_chrom = resolve_contig(bam, variants[0][0])
        for read in bam.fetch(fetch_chrom, max(0, start - 1), end):
            if read.mapping_quality < min_mapq or read.is_secondary or read.is_supplementary or read.is_qcfail:
                continue
            if read.is_duplicate and not include_duplicates:
                continue
            key = read.query_name
            for i, variant in enumerate(variants):
                if read.reference_start > variant[1] - 1 or read.reference_end is None or read.reference_end < variant[1]:
                    continue
                call = call_variant(read, variant, min_baseq)
                if call:
                    add_call(fragments[key], i, call[0], min(read.mapping_quality, call[1]))
    finally:
        bam.close()
    return {name: {i: c for i, c in calls.items() if c is not None} for name, calls in fragments.items()}


def parent_origin_from_statuses(mother_status, father_status):
    if mother_status == "has_alt" and father_status == "no_alt":
        return "maternal"
    if father_status == "has_alt" and mother_status == "no_alt":
        return "paternal"
    return ""


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


def classify_origin_hint(index, variants, mother_statuses, father_statuses, window, fragments=None):
    if mother_statuses[index] != "no_alt" or father_statuses[index] != "no_alt":
        return ""
    chrom, pos, _, _ = variants[index]
    maternal = paternal = 0
    if fragments:
        origins = [parent_origin_from_statuses(m, f) for m, f in zip(mother_statuses, father_statuses)]
        for calls in fragments.values():
            if index not in calls:
                continue
            for j in calls:
                if j == index or variants[j][0] != chrom or abs(variants[j][1] - pos) > window:
                    continue
                origin = origins[j]
                if origin == "maternal":
                    maternal += 1
                elif origin == "paternal":
                    paternal += 1
    else:
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
    parser.add_argument("--bam-column", default="bam", help="Proband BAM/CRAM column")
    parser.add_argument("--reference-column", default="reference", help="Reference FASTA column")
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
    parser.add_argument("--min-mapq", type=int, default=DEFAULT_MIN_MAPQ, help="Minimum MAPQ for read-backed BAM calls")
    parser.add_argument("--min-baseq", type=int, default=DEFAULT_MIN_BASEQ, help="Minimum base quality for read-backed BAM calls")
    parser.add_argument("--include-duplicates", action="store_true", help="Include duplicate reads in BAM evidence")
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
            row_value(row, args.bam_column),
            row_value(row, args.reference_column),
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
        sample_bam, reference = key[1], key[2]
        mother_vcf, mother_sample = key[3], key[4]
        father_vcf, father_sample = key[5], key[6]
        mother_statuses = statuses_for_vcf(mother_vcf, mother_sample, all_variants)
        father_statuses = statuses_for_vcf(father_vcf, father_sample, all_variants)
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
                    nearby = [j for j, other in enumerate(variants) if j != i and other[0] == variants[i][0] and abs(other[1] - variants[i][1]) <= args.origin_window]
                    if nearby and sample_bam and os.path.exists(sample_bam):
                        local_indices = [i] + nearby
                        local_variants = [variants[j] for j in local_indices]
                        local_mother = [row_mother_statuses[j] for j in local_indices]
                        local_father = [row_father_statuses[j] for j in local_indices]
                        span_start = min(v[1] for v in local_variants) - args.origin_window
                        span_end = max(v[1] + len(v[2]) for v in local_variants) + args.origin_window
                        fragments = read_fragments(sample_bam, reference, local_variants, span_start, span_end, args.min_mapq, args.min_baseq, args.include_duplicates)
                        hint = classify_origin_hint(0, local_variants, local_mother, local_father, args.origin_window, fragments)
                    else:
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
