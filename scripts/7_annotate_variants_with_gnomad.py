#!/usr/bin/env python3
"""Annotate variant TSV rows with gnomAD frequency."""

import argparse
import csv
import gzip
import os
import re
from concurrent.futures import ThreadPoolExecutor, as_completed

DEFAULT_BATCH_SIZE = 25
DEFAULT_WORKERS = 8
DEFAULT_FETCH_PADDING = 1
DEFAULT_GNOMAD_VCF_TEMPLATE = "gnomad.joint.v4.1.sites.chr{chrom}.vcf.bgz"

# Toggle detailed debug output when set by CLI
DEBUG = False


def debug_print(message):
    if DEBUG:
        print(f"DEBUG: {message}", flush=True)


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


def preview_vcf_header(vcf_path, max_lines=10):
    try:
        opener = gzip.open if vcf_path.endswith((".gz", ".bgz", ".bgzf")) else open
        with opener(vcf_path, "rt", errors="replace") as fh:
            lines = []
            for _ in range(max_lines):
                line = fh.readline()
                if not line:
                    break
                lines.append(line.rstrip("\n"))
            return lines
    except Exception as exc:
        return [f"<failed to read header preview: {exc}>"]


def open_gnomad_vcf(vcf_path):
    try:
        import pysam
    except ImportError as exc:
        raise SystemExit("This script requires pysam to query a local indexed VCF") from exc

    try:
        debug_print(f"opening VCF {vcf_path}")
        return pysam.VariantFile(vcf_path)
    except Exception as exc:
        header_preview = preview_vcf_header(vcf_path)
        message = [
            f"Failed to open gnomAD VCF with pysam: {exc}",
            f"VCF path: {vcf_path}",
            "Header preview:",
            *header_preview,
            "Common causes: bad/corrupt header, wrong file path, stale or missing index, or a file that is not a valid bgzipped VCF.",
        ]
        raise SystemExit("\n".join(message))


def ensure_indexed_vcf(vcf_path):
    if not os.path.exists(vcf_path):
        raise SystemExit(f"gnomAD VCF not found: {vcf_path}")
    if not (os.path.exists(vcf_path + ".tbi") or os.path.exists(vcf_path + ".csi")):
        raise SystemExit(f"Missing gnomAD VCF index for {vcf_path} (.tbi or .csi)")
    debug_print(f"found VCF and index for {vcf_path}")


def close_gnomad_vcf(reader):
    close = getattr(reader, "close", None)
    if callable(close):
        close()


def reader_contigs(reader):
    return set(reader.header.contigs)


def fetch_region_records(reader, chrom, start, end):
    return list(reader.fetch(chrom, max(0, start - 1), end))


def record_contig(record):
    return record.chrom


def record_pos(record):
    return int(record.pos)


def record_ref(record):
    return str(record.ref).upper()


def record_alts(record):
    values = record.alts
    # Normalize ALT values; handle bytes/bytearray and other container types
    alts = []
    for value in (values or []):
        if value is None or value == ".":
            continue
        try:
            if isinstance(value, (bytes, bytearray)):
                v = value.decode()
            else:
                v = str(value)
        except Exception:
            v = str(value)
        v = v.upper()
        if v:
            alts.append(v)
    return alts


def record_info(record, key):
    info = record.info
    val = info.get(key) if hasattr(info, "get") else getattr(info, key, None)
    if isinstance(val, (bytes, bytearray)):
        try:
            return val.decode()
        except Exception:
            return str(val)
    return val


def coerce_info_value(value, alt_index=None):
    if value is None:
        return None

    # Handle bytes/bytearray
    if isinstance(value, (bytes, bytearray)):
        try:
            value = value.decode()
        except Exception:
            value = str(value)

    # Convert numpy arrays or other array-like containers to list
    try:
        import numpy as _np

        if isinstance(value, _np.ndarray):
            value = value.tolist()
    except Exception:
        pass

    # Handle list/tuple-like values (per-allele annotations)
    if isinstance(value, (list, tuple)):
        if not value:
            return None
        # decode bytes inside containers
        converted = []
        for v in value:
            if isinstance(v, (bytes, bytearray)):
                try:
                    converted.append(v.decode())
                except Exception:
                    converted.append(str(v))
            else:
                converted.append(v)
        if alt_index is not None and alt_index < len(converted):
            return converted[alt_index]
        return converted[0]

    # Strings with comma-separated values (sometimes used for per-allele fields)
    if isinstance(value, str) and "," in value:
        parts = [part for part in value.split(",") if part != ""]
        if not parts:
            return None
        if alt_index is not None and alt_index < len(parts):
            value = parts[alt_index]
        else:
            value = parts[0]

    return value


def get_annotated_value(record, keys, alt_index=None, numeric=float):
    for key in keys:
        value = coerce_info_value(record_info(record, key), alt_index=alt_index)
        if value is None:
            continue
        try:
            return numeric(value)
        except Exception:
            try:
                return numeric(str(value))
            except Exception:
                continue
    return None


def parse_variant_result(record, fallback_id, variant_id):
    if record is None:
        return {
            "gnomad_variant_id": fallback_id,
            "gnomad_af": 0.0,
            "gnomad_nhomalt": 0,
            "gnomad_lookup_status": "not_found",
        }
    alts = record_alts(record)
    if not alts or variant_id[3] not in alts:
        return {
            "gnomad_variant_id": fallback_id,
            "gnomad_af": 0.0,
            "gnomad_nhomalt": 0,
            "gnomad_lookup_status": "not_found",
        }
    alt_index = alts.index(variant_id[3])
    af = get_annotated_value(record, ("AF_joint", "AF", "af"), alt_index=alt_index, numeric=float)
    nhomalt = get_annotated_value(record, ("nhomalt_joint", "NHOMALT_joint", "nhomalt", "NHOMALT", "n_homalt", "HOMALT", "homozygote_count"), alt_index=alt_index, numeric=int) or 0
    return {
        "gnomad_variant_id": fallback_id,
        "gnomad_af": af,
        "gnomad_nhomalt": nhomalt,
        "gnomad_lookup_status": "found",
    }


def query_variant(reader, variant_id):
    chrom, pos, ref, alt = variant_id.split("-")
    pos = int(pos)
    start = max(1, pos - DEFAULT_FETCH_PADDING)
    end = pos + len(ref) + DEFAULT_FETCH_PADDING
    for alias in contig_aliases(chrom):
        if alias not in reader_contigs(reader):
            continue
        for record in fetch_region_records(reader, alias, start, end):
            if DEBUG:
                try:
                    info_keys = list(record.info.keys()) if hasattr(record.info, 'keys') else []
                except Exception:
                    info_keys = []
                print(f"DEBUG: fetched record: CHROM={record.chrom} POS={record.pos} REF={record.ref} ALTS={record.alts} INFO_KEYS={info_keys}")
            if record_contig(record) not in contig_aliases(chrom):
                if DEBUG:
                    print(f"DEBUG: skipping record on contig mismatch: {record.chrom} not in {contig_aliases(chrom)}")
                continue
            if record_pos(record) != pos:
                if DEBUG:
                    print(f"DEBUG: skipping record at POS={record.pos} (want {pos})")
                continue
            if record_ref(record) != ref:
                if DEBUG:
                    print(f"DEBUG: REF mismatch: record REF={record_ref(record)} vs want {ref}")
                continue
            # At this point position and REF match; parse
            result = parse_variant_result(record, variant_id, (chrom, pos, ref, alt))
            if result["gnomad_lookup_status"] == "found":
                return result
            if DEBUG:
                print(f"DEBUG: ALT {alt} not found in record ALTS={record_alts(record)}")
    return parse_variant_result(None, variant_id, (chrom, pos, ref, alt))


def query_batch(vcf_path, batch):
    if not batch:
        return {}
    return _query_batch_single_vcf(vcf_path, batch)


def _query_batch_single_vcf(vcf_path, batch):
    ensure_indexed_vcf(vcf_path)
    reader = open_gnomad_vcf(vcf_path)
    try:
        return {variant_id: query_variant(reader, variant_id) for variant_id in batch}
    finally:
        close_gnomad_vcf(reader)


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
    parser.add_argument("--debug", action="store_true", help="Enable debug output to help diagnose matching issues")
    args = parser.parse_args()

    global DEBUG
    DEBUG = bool(getattr(args, "debug", False))

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
    debug_print(f"input TSV: {args.input_tsv}")
    debug_print(f"gnomAD source: {args.gnomad_vcf or args.gnomad_vcf_dir}")
    if args.gnomad_vcf:
        debug_print(f"gnomAD VCF header preview: {preview_vcf_header(args.gnomad_vcf)}")

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
            batch_annotations = future.result()
            debug_print(f"batch {batch_index}/{len(batches)} returned {len(batch_annotations)} annotations")
            annotations.update(batch_annotations)

    output_fields = fieldnames + ["gnomad_variant_id", "gnomad_af", "gnomad_nhomalt", "gnomad_lookup_status"]

    with open(args.out, "w", newline="") as out_fh:
        writer = csv.DictWriter(out_fh, delimiter="\t", fieldnames=output_fields)
        writer.writeheader()
        for row, variant_ids in zip(rows, normalized_by_row):
            anns = [annotations[variant_id] for variant_id in variant_ids]
            output_row = dict(row)
            output_row["gnomad_variant_id"] = ";".join(ann["gnomad_variant_id"] for ann in anns)
            output_row["gnomad_af"] = ";".join(str(ann["gnomad_af"]) for ann in anns)
            output_row["gnomad_nhomalt"] = ";".join(str(ann["gnomad_nhomalt"]) for ann in anns)
            output_row["gnomad_lookup_status"] = ";".join(ann["gnomad_lookup_status"] for ann in anns)
            writer.writerow(output_row)

    print(f"Wrote annotated TSV to {args.out}", flush=True)


if __name__ == "__main__":
    main()
