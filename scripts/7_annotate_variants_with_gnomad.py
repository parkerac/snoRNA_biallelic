#!/usr/bin/env python3
"""Annotate variant TSV rows with gnomAD frequency."""

import argparse
import csv
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

DEFAULT_BATCH_SIZE = 25
DEFAULT_WORKERS = 8
DEFAULT_FETCH_PADDING = 1
DEFAULT_GNOMAD_VCF_TEMPLATE = "gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"


def normalize_variant_id(value):
    value = str(value).strip()
    if not value:
        raise ValueError("empty variant id")
    parts = re.split(r"[:\-]", value)
    if len(parts) != 4:
        raise ValueError(f"{value!r} must look like chr:pos:ref:alt or chr-pos-ref-alt")
    chrom, pos, ref, alt = parts
    chrom = chrom[3:] if chrom.lower().startswith("chr") else chrom
    if chrom.upper() in {"MT", "M"}:
        chrom = "M"
    return f"{chrom}-{int(pos)}-{ref.upper()}-{alt.upper()}"


def split_variant_values(value):
    return [item.strip() for item in str(value).split(";") if item.strip()]


def contig_aliases(chrom):
    aliases = [chrom]
    aliases.append(chrom[3:] if chrom.startswith("chr") else f"chr{chrom}")
    if chrom in {"M", "MT", "chrM", "chrMT"}:
        aliases.extend(["M", "MT", "chrM", "chrMT"])
    return list(dict.fromkeys(aliases))


def gnomad_vcf_path_for_chrom(vcf_dir, template, chrom):
    return os.path.join(vcf_dir, template.format(chrom=chrom))


def open_gnomad_vcf(vcf_path):
    cyvcf2_error = None
    try:
        from cyvcf2 import VCF

        return "cyvcf2", VCF(vcf_path)
    except (ImportError, ValueError, OSError) as exc:
        cyvcf2_error = exc

    try:
        import pysam
    except ImportError as exc:
        if cyvcf2_error is not None:
            raise SystemExit(f"Failed to open gnomAD VCF with cyvcf2 ({cyvcf2_error}) and pysam is not installed") from exc
        raise SystemExit("This script requires either cyvcf2 or pysam to query a local indexed VCF") from exc

    try:
        return "pysam", pysam.VariantFile(vcf_path)
    except Exception as exc:
        if cyvcf2_error is not None:
            raise SystemExit(f"Failed to open gnomAD VCF with cyvcf2 ({cyvcf2_error}) and pysam ({exc})")
        raise


def ensure_indexed_vcf(vcf_path):
    if not os.path.exists(vcf_path):
        raise SystemExit(f"gnomAD VCF not found: {vcf_path}")
    if not (os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi")):
        raise SystemExit(f"Missing gnomAD VCF index for {vcf_path} (.tbi or .csi)")


def close_gnomad_vcf(kind, reader):
    close = getattr(reader, "close", None)
    if callable(close):
        close()


def reader_contigs(kind, reader):
    if kind == "cyvcf2":
        return set(reader.seqnames)
    return set(reader.header.contigs)


def fetch_region_records(kind, reader, chrom, start, end):
    if kind == "cyvcf2":
        region = f"{chrom}:{max(1, start)}-{end}"
        return list(reader(region))
    return list(reader.fetch(chrom, max(0, start - 1), end))


def record_contig(kind, record):
    return record.CHROM if kind == "cyvcf2" else record.chrom


def record_pos(kind, record):
    return int(record.POS if kind == "cyvcf2" else record.pos)


def record_ref(kind, record):
    return str(record.REF if kind == "cyvcf2" else record.ref).upper()


def record_alts(kind, record):
    values = record.ALT if kind == "cyvcf2" else record.alts
    return [str(value).upper() for value in (values or []) if value and value != "."]


def record_info(kind, record, key):
    if kind == "cyvcf2":
        info = record.INFO
        return info.get(key) if hasattr(info, "get") else getattr(info, key, None)
    info = record.info
    return info.get(key) if hasattr(info, "get") else getattr(info, key, None)


def coerce_info_value(value, alt_index=None):
    if value is None:
        return None
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        if alt_index is not None and alt_index < len(value):
            return value[alt_index]
        return value[0]
    if isinstance(value, str) and "," in value:
        parts = [part for part in value.split(",") if part != ""]
        if not parts:
            return None
        if alt_index is not None and alt_index < len(parts):
            value = parts[alt_index]
        else:
            value = parts[0]
    return value


def get_annotated_value(kind, record, keys, alt_index=None, numeric=float):
    for key in keys:
        value = coerce_info_value(record_info(kind, record, key), alt_index=alt_index)
        if value is None:
            continue
        try:
            return numeric(value)
        except Exception:
            continue
    return None


def parse_variant_result(kind, record, fallback_id, variant_id):
    if record is None:
        return {
            "gnomad_variant_id": fallback_id,
            "gnomad_ac": 0,
            "gnomad_an": 0,
            "gnomad_af": 0.0,
            "gnomad_nhomalt": 0,
            "gnomad_lookup_status": "not_found",
        }
    alts = record_alts(kind, record)
    if not alts or variant_id[3] not in alts:
        return {
            "gnomad_variant_id": fallback_id,
            "gnomad_ac": 0,
            "gnomad_an": 0,
            "gnomad_af": 0.0,
            "gnomad_nhomalt": 0,
            "gnomad_lookup_status": "not_found",
        }
    alt_index = alts.index(variant_id[3])
    ac = get_annotated_value(kind, record, ("AC_joint", "AC", "ac"), alt_index=alt_index, numeric=int) or 0
    an = get_annotated_value(kind, record, ("AN_joint", "AN", "an"), alt_index=None, numeric=int) or 0
    af = get_annotated_value(kind, record, ("AF_joint", "AF", "af"), alt_index=alt_index, numeric=float)
    if af is None:
        af = (ac / an) if an else 0.0
    nhomalt = get_annotated_value(kind, record, ("NHOMALT_joint", "nhomalt", "NHOMALT", "n_homalt", "HOMALT", "homozygote_count"), alt_index=alt_index, numeric=int) or 0
    return {
        "gnomad_variant_id": fallback_id,
        "gnomad_ac": ac,
        "gnomad_an": an,
        "gnomad_af": af,
        "gnomad_nhomalt": nhomalt,
        "gnomad_lookup_status": "found",
    }


def query_variant(kind, reader, variant_id):
    chrom, pos, ref, alt = variant_id.split("-")
    pos = int(pos)
    start = max(1, pos - DEFAULT_FETCH_PADDING)
    end = pos + len(ref) + DEFAULT_FETCH_PADDING
    for alias in contig_aliases(chrom):
        if alias not in reader_contigs(kind, reader):
            continue
        for record in fetch_region_records(kind, reader, alias, start, end):
            if record_contig(kind, record) not in contig_aliases(chrom):
                continue
            if record_pos(kind, record) != pos:
                continue
            if record_ref(kind, record) != ref:
                continue
            return parse_variant_result(kind, record, variant_id, (chrom, pos, ref, alt))
    return parse_variant_result(kind, None, variant_id, (chrom, pos, ref, alt))


def query_batch(vcf_path, batch):
    if not batch:
        return {}
    return _query_batch_single_vcf(vcf_path, batch)


def _query_batch_single_vcf(vcf_path, batch):
    ensure_indexed_vcf(vcf_path)
    kind, reader = open_gnomad_vcf(vcf_path)
    try:
        return {variant_id: query_variant(kind, reader, variant_id) for variant_id in batch}
    finally:
        close_gnomad_vcf(kind, reader)


def query_batch_by_chrom(vcf_dir, template, batch):
    if not batch:
        return {}
    grouped = {}
    for variant_id in batch:
        chrom = variant_id.split("-")[0]
        grouped.setdefault(gnomad_vcf_path_for_chrom(vcf_dir, template, chrom), []).append(variant_id)
    annotations = {}
    for vcf_path, variant_ids in grouped.items():
        annotations.update(_query_batch_single_vcf(vcf_path, variant_ids))
    return annotations


def chunked(values, size):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-tsv", required=True, help="Input TSV containing a variant column")
    parser.add_argument("--variant-column", default="variant_id", help="Column containing chr:pos:ref:alt variant IDs, optionally semicolon-separated")
    parser.add_argument("--out", required=True, help="Annotated output TSV path")
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--gnomad-vcf", help="Single local bgzipped gnomAD VCF (.vcf.gz)")
    group.add_argument("--gnomad-vcf-dir", help="Directory containing chromosome-specific gnomAD VCFs")
    parser.add_argument("--gnomad-vcf-template", default=DEFAULT_GNOMAD_VCF_TEMPLATE, help="Filename template under --gnomad-vcf-dir; use {chrom} for the chromosome name")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Number of unique variants to query per gnomAD request")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help="Number of query batches to run in parallel")
    args = parser.parse_args()

    if args.gnomad_vcf:
        ensure_indexed_vcf(args.gnomad_vcf)
    else:
        if not os.path.isdir(args.gnomad_vcf_dir):
            raise SystemExit(f"gnomAD VCF directory not found: {args.gnomad_vcf_dir}")

    with open(args.input_tsv, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        fieldnames = reader.fieldnames or []
        if args.variant_column not in fieldnames:
            raise ValueError(f"{args.input_tsv} is missing column {args.variant_column!r}")
        rows = list(reader)

    normalized_by_row = []
    unique_variants = []
    seen = set()
    for row in rows:
        variant_ids = [normalize_variant_id(value) for value in split_variant_values(row[args.variant_column])]
        if not variant_ids:
            raise ValueError(f"row is missing a variant in column {args.variant_column!r}")
        normalized_by_row.append(variant_ids)
        for variant_id in variant_ids:
            if variant_id not in seen:
                seen.add(variant_id)
                unique_variants.append(variant_id)

    print(f"Loaded {len(rows)} rows with {len(unique_variants)} unique variants", flush=True)

    batches = list(chunked(unique_variants, args.batch_size))
    workers = max(1, min(args.workers, len(batches) or 1))
    print(f"Querying local gnomAD VCF in {len(batches)} batches with {workers} workers", flush=True)

    annotations = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_to_batch = {}
        for i, batch in enumerate(batches, start=1):
            if args.gnomad_vcf:
                future = pool.submit(query_batch, args.gnomad_vcf, batch)
            else:
                future = pool.submit(query_batch_by_chrom, args.gnomad_vcf_dir, args.gnomad_vcf_template, batch)
            future_to_batch[future] = (i, batch)
        for future in as_completed(future_to_batch):
            batch_index, batch = future_to_batch[future]
            print(f"Finished gnomAD batch {batch_index}/{len(batches)}", flush=True)
            annotations.update(future.result())

    output_fields = fieldnames + ["gnomad_variant_id", "gnomad_ac", "gnomad_an", "gnomad_af", "gnomad_nhomalt", "gnomad_lookup_status"]

    with open(args.out, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, delimiter="\t", fieldnames=output_fields)
        writer.writeheader()
        for row, variant_ids in zip(rows, normalized_by_row):
            anns = [annotations[variant_id] for variant_id in variant_ids]
            output_row = dict(row)
            output_row["gnomad_variant_id"] = ";".join(ann["gnomad_variant_id"] for ann in anns)
            output_row["gnomad_ac"] = ";".join(str(ann["gnomad_ac"]) for ann in anns)
            output_row["gnomad_an"] = ";".join(str(ann["gnomad_an"]) for ann in anns)
            output_row["gnomad_af"] = ";".join(str(ann["gnomad_af"]) for ann in anns)
            output_row["gnomad_nhomalt"] = ";".join(str(ann["gnomad_nhomalt"]) for ann in anns)
            output_row["gnomad_lookup_status"] = ";".join(ann["gnomad_lookup_status"] for ann in anns)
            writer.writerow(output_row)

    print(f"Wrote annotated TSV to {args.out}", flush=True)


if __name__ == "__main__":
    main()
