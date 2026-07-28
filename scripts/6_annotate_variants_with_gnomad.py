#!/usr/bin/env python3
"""Annotate variant TSV rows with gnomAD frequency."""

import argparse
import csv
import json
import re
import time
from concurrent.futures import ThreadPoolExecutor, as_completed
from urllib.error import HTTPError, URLError
from urllib.request import Request, urlopen


GNOMAD_API_URL = "https://gnomad.broadinstitute.org/api"
DEFAULT_DATASET = "gnomad_r4"
DEFAULT_BATCH_SIZE = 25
DEFAULT_WORKERS = 8
DEFAULT_SLEEP_SECONDS = 6
DEFAULT_RETRIES = 3
MAX_GRAPHQL_BATCH = 25


def scalar_int(value):
    try:
        return int(value)
    except Exception:
        return 0


def scalar_float(value):
    try:
        return float(value)
    except Exception:
        return None


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


def build_query(batch, dataset):
    lines = ["query VariantBatch($dataset: DatasetId!) {"]
    for i, variant_id in enumerate(batch):
        alias = f"v{i}"
        lines.append(
            f'  {alias}: variant(variantId: {json.dumps(variant_id)}, dataset: $dataset) '
            "{ variantId joint { ac an af nhomalt } }"
        )
    lines.append("}")
    return "\n".join(lines)


def post_query(api_url, query, dataset):
    payload = json.dumps({"query": query, "variables": {"dataset": dataset}}).encode()
    request = Request(
        api_url,
        data=payload,
        headers={
            "Content-Type": "application/json",
            "User-Agent": "Mozilla/5.0",
        },
    )
    with urlopen(request) as response:
        return json.loads(response.read().decode())


def query_batch(api_url, dataset, batch, retries, sleep_seconds):
    if not batch:
        return {}
    if len(batch) > MAX_GRAPHQL_BATCH:
        raise ValueError(f"batch size cannot exceed {MAX_GRAPHQL_BATCH}")

    query = build_query(batch, dataset)
    last_error = None
    for attempt in range(retries):
        try:
            result = post_query(api_url, query, dataset)
            if result.get("errors"):
                last_error = result["errors"]
                break
            data = result.get("data") or {}
            return {
                variant_id: parse_variant_result(data.get(f"v{i}"), variant_id)
                for i, variant_id in enumerate(batch)
            }
        except HTTPError as exc:
            last_error = f"HTTP {exc.code}: {exc.read().decode(errors='replace')}"
            if exc.code in {403, 429, 500, 502, 503, 504} and attempt < retries - 1:
                time.sleep(sleep_seconds * (attempt + 1))
                continue
            break
        except (URLError, TimeoutError, json.JSONDecodeError) as exc:
            last_error = str(exc)
            if attempt < retries - 1:
                time.sleep(sleep_seconds * (attempt + 1))
                continue
            break

    if len(batch) == 1:
        variant_id = batch[0]
        return {
            variant_id: {
                "gnomad_variant_id": variant_id,
                "gnomad_ac": 0,
                "gnomad_an": 0,
                "gnomad_af": 0.0,
                "gnomad_lookup_status": f"error:{last_error}" if last_error else "error",
            }
        }

    mid = len(batch) // 2
    left = query_batch(api_url, dataset, batch[:mid], retries, sleep_seconds)
    right = query_batch(api_url, dataset, batch[mid:], retries, sleep_seconds)
    left.update(right)
    return left


def chunked(values, size):
    for start in range(0, len(values), size):
        yield values[start : start + size]


def parse_variant_result(payload, fallback_id):
    if not payload or payload.get("variant") is None:
        return {
            "gnomad_variant_id": fallback_id,
            "gnomad_ac": 0,
            "gnomad_an": 0,
            "gnomad_af": 0.0,
            "gnomad_lookup_status": "not_found",
        }
    variant = payload["variant"]
    joint = variant.get("joint") or {}
    ac = scalar_int(joint.get("ac"))
    an = scalar_int(joint.get("an"))
    af = scalar_float(joint.get("af"))
    if af is None:
        af = (ac / an) if an else 0.0
    return {
        "gnomad_variant_id": variant.get("variantId") or fallback_id,
        "gnomad_ac": ac,
        "gnomad_an": an,
        "gnomad_af": af,
        "gnomad_nhomalt": scalar_int(joint.get("nhomalt")),
        "gnomad_lookup_status": "found",
    }


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--input-tsv", required=True, help="Input TSV containing a variant column")
    parser.add_argument("--variant-column", default="variant_id", help="Column containing chr:pos:ref:alt variant IDs, optionally semicolon-separated")
    parser.add_argument("--out", required=True, help="Annotated output TSV path")
    parser.add_argument("--dataset", default=DEFAULT_DATASET, help="gnomAD dataset id, for example gnomad_r4")
    parser.add_argument("--batch-size", type=int, default=DEFAULT_BATCH_SIZE, help="Number of unique variants to query per gnomAD request")
    parser.add_argument("--workers", type=int, default=DEFAULT_WORKERS, help="Number of query batches to run in parallel")
    parser.add_argument("--sleep-seconds", type=float, default=DEFAULT_SLEEP_SECONDS, help="Seconds to wait between failed retries")
    parser.add_argument("--retries", type=int, default=DEFAULT_RETRIES, help="Retry count for gnomAD requests")
    parser.add_argument("--api-url", default=GNOMAD_API_URL, help="gnomAD GraphQL endpoint")
    args = parser.parse_args()

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
    print(f"Querying gnomAD in {len(batches)} batches with {workers} workers", flush=True)

    annotations = {}
    with ThreadPoolExecutor(max_workers=workers) as pool:
        future_to_batch = {
            pool.submit(query_batch, args.api_url, args.dataset, batch, args.retries, args.sleep_seconds): (i, batch)
            for i, batch in enumerate(batches, start=1)
        }
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
