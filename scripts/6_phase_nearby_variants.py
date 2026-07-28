#!/usr/bin/env python3
"""Phase two nearby variants from local BAM/CRAM read evidence without downsampling."""

import argparse
import csv
import json
import math
import multiprocessing as mp
import os
from collections import defaultdict
from dataclasses import dataclass

pysam = None
MISSING_VALUES = {"", ".", "NA", "N/A", "NONE", "NULL", "NAN"}
REFERENCE_CACHE = {}
PAIR_COLUMNS = {
    "chrom",
    "pos1",
    "ref1",
    "alt1",
    "pos2",
    "ref2",
    "alt2",
    "bam",
    "reference",
    "vcf",
    "sample",
    "father_vcf",
    "mother_vcf",
    "father_sample",
    "mother_sample",
    "father_bam",
    "mother_bam",
}


@dataclass(frozen=True)
class Variant:
    chrom: str
    pos: int
    ref: str
    alt: str
    name: str
    target: bool = False


def parse_variant(value, name, target=True):
    parts = value.replace(",", ":").split(":")
    if len(parts) != 4:
        raise ValueError(f"{name} must look like chrom:pos:ref:alt")
    return Variant(parts[0], int(parts[1]), parts[2].upper(), parts[3].upper(), name, target)


def parse_gt(gt):
    if gt is None:
        return None
    alleles = [a for a in gt if a is not None]
    if len(alleles) != 2 or alleles[0] == alleles[1]:
        return None
    return set(alleles)


def same_variant(a, b):
    return b.chrom in contig_aliases(a.chrom) and a.pos == b.pos and a.ref == b.ref and a.alt == b.alt


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
    raise ValueError(f"Contig {chrom} was not found in input header; tried {', '.join(contig_aliases(chrom))}")


def load_bridge_variants(vcf_path, sample, chrom, start, end, targets, max_bridges):
    if not vcf_path:
        return []
    variants = []
    with pysam.VariantFile(vcf_path) as vcf:
        sample = sample or next(iter(vcf.header.samples), None)
        if not sample:
            raise ValueError("VCF has no samples; provide targets only or a genotyped VCF")
        for rec in vcf.fetch(resolve_contig(vcf, chrom), max(0, start - 1), end):
            real_alts = real_alt_indices(rec)
            if len(rec.ref) > 50 or len(real_alts) != 1:
                continue
            gt = parse_gt(rec.samples[sample].get("GT"))
            alt, alt_index = next(iter(real_alts.items()))
            if gt != {0, alt_index}:
                continue
            variant = Variant(rec.chrom, rec.pos, rec.ref.upper(), alt.upper(), rec.id or f"{rec.chrom}:{rec.pos}:{rec.ref}>{alt}")
            if any(same_variant(variant, target) for target in targets):
                continue
            variants.append(variant)
            if len(variants) >= max_bridges:
                break
    return variants


def infer_origin(mother_status, father_status):
    if mother_status == "low_quality" or father_status == "low_quality":
        return "low_quality"
    if mother_status == "has_alt" and father_status == "no_alt":
        return "maternal"
    if father_status == "has_alt" and mother_status == "no_alt":
        return "paternal"
    if mother_status == "no_alt" and father_status == "no_alt":
        return "conflicting"
    return "uninformative"


def gt_alleles(rec, sample):
    if not sample:
        sample = next(iter(rec.samples), None)
    if not sample or sample not in rec.samples:
        return None
    gt = rec.samples[sample].get("GT")
    if gt is None or any(a is None for a in gt):
        return None
    return set(gt)


def real_alt_indices(rec):
    return {alt.upper(): i + 1 for i, alt in enumerate(rec.alts or []) if alt and not alt.startswith("<") and alt != "*"}


def classify_parent_vcf_record(rec, variant, sample):
    alleles = gt_alleles(rec, sample)
    if alleles is None:
        return "ambiguous"
    if rec.pos == variant.pos and rec.ref.upper() == variant.ref and variant.alt in real_alt_indices(rec):
        return "has_alt" if real_alt_indices(rec)[variant.alt] in alleles else "no_alt"
    if rec.start <= variant.pos - 1 < rec.stop and alleles == {0}:
        return "no_alt"
    return "ambiguous" if any(a > 0 for a in alleles) else "no_alt"


def count_parent_vcf_variant(vcf, sample, variant):
    saw_record = False
    try:
        records = vcf.fetch(resolve_contig(vcf, variant.chrom), max(0, variant.pos - 1), variant.pos + len(variant.ref))
    except ValueError:
        return {"status": "no_record"}
    for rec in records:
        if not (rec.start <= variant.pos - 1 < rec.stop):
            continue
        saw_record = True
        status = classify_parent_vcf_record(rec, variant, sample)
        if status in {"has_alt", "no_alt"}:
            return {"status": status}
    return {"status": "ambiguous" if saw_record else "no_record"}


def parent_vcf_statuses_for_path(vcf_path, sample, variants):
    if not vcf_path:
        return {i: {"status": "missing"} for i in range(len(variants))}
    if not os.path.exists(vcf_path):
        return {i: {"status": "missing_file"} for i in range(len(variants))}
    with pysam.VariantFile(vcf_path) as vcf:
        sample = sample or next(iter(vcf.header.samples), None)
        if not sample or sample not in vcf.header.samples:
            return {i: {"status": "missing_sample"} for i in range(len(variants))}
        return {i: count_parent_vcf_variant(vcf, sample, variant) for i, variant in enumerate(variants)}


def parent_vcf_statuses(args, variants):
    return {
        "mother": parent_vcf_statuses_for_path(args.mother_vcf, args.mother_sample, variants),
        "father": parent_vcf_statuses_for_path(args.father_vcf, args.father_sample, variants),
    }


def parent_origins(statuses):
    return {i: infer_origin(statuses["mother"][i]["status"], statuses["father"][i]["status"]) for i in statuses["mother"]}


def phase_from_origins(origin1, origin2):
    if origin1 == "conflicting" or origin2 == "conflicting":
        return "conflicting"
    if origin1 not in {"maternal", "paternal"} or origin2 not in {"maternal", "paternal"}:
        return "ambiguous"
    return "cis" if origin1 == origin2 else "trans"


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
    start = variant.pos - 1
    end = start + len(variant.ref)
    pairs = read.get_aligned_pairs(matches_only=False)
    ref_to_qpos = {rpos: qpos for qpos, rpos in pairs if rpos is not None}
    qpos = [ref_to_qpos.get(pos) for pos in range(start, end)]
    quals = read.query_qualities or []
    observed = "".join("-" if q is None else read.query_sequence[q].upper() for q in qpos)
    baseq = min([quals[q] for q in qpos if q is not None] or [0])

    if len(variant.ref) == len(variant.alt):
        if None in qpos or baseq < min_baseq:
            return None
        if observed == variant.ref:
            return 0, baseq
        if observed == variant.alt:
            return 1, baseq
        return None

    if len(variant.ref) > len(variant.alt) and variant.ref.startswith(variant.alt):
        padded_alt = variant.alt + "-" * (len(variant.ref) - len(variant.alt))
        if baseq >= min_baseq and observed == variant.ref:
            return 0, baseq
        if baseq >= min_baseq and observed == padded_alt:
            return 1, baseq
        return None

    if len(variant.alt) > len(variant.ref) and variant.alt.startswith(variant.ref):
        if None in qpos or baseq < min_baseq or observed != variant.ref:
            return None
        last_ref = end - 1
        pair_index = next((i for i, (q, r) in enumerate(pairs) if r == last_ref and q == qpos[-1]), None)
        inserted = insertion_after(pairs, pair_index, read) if pair_index is not None else ""
        if inserted == "":
            return 0, baseq
        if variant.ref + inserted == variant.alt:
            return 1, baseq
    return None


def add_call(fragment, variant_index, allele, quality):
    old = fragment.get(variant_index)
    if old is None or old[0] == allele and quality > old[1]:
        fragment[variant_index] = (allele, quality)
    elif old[0] != allele:
        fragment[variant_index] = None


def is_cram(path):
    return path.lower().endswith(".cram")


def cached_reference(reference):
    if not reference:
        return None
    reference = os.path.realpath(reference)
    if reference not in REFERENCE_CACHE and pysam is not None:
        try:
            REFERENCE_CACHE[reference] = pysam.FastaFile(reference)
        except Exception:
            REFERENCE_CACHE[reference] = None
    return reference


def read_fragments(bam_path, reference, variants, start, end, min_mapq, min_baseq, include_duplicates):
    fragments = defaultdict(dict)
    if is_cram(bam_path) and not reference:
        raise ValueError(f"{bam_path} is a CRAM; provide the shared reference with --reference")
    reference = cached_reference(reference)
    bam = pysam.AlignmentFile(bam_path, "rc" if is_cram(bam_path) else "rb", reference_filename=reference)
    fetch_chrom = resolve_contig(bam, variants[0].chrom)
    for read in bam.fetch(fetch_chrom, max(0, start - 1), end):
        if read.mapping_quality < min_mapq or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.is_duplicate and not include_duplicates:
            continue
        key = read.query_name
        for i, variant in enumerate(variants):
            if read.reference_start > variant.pos - 1 or read.reference_end is None or read.reference_end < variant.pos:
                continue
            call = call_variant(read, variant, min_baseq)
            if call:
                add_call(fragments[key], i, call[0], min(read.mapping_quality, call[1]))
    bam.close()
    return {name: {i: c for i, c in calls.items() if c is not None} for name, calls in fragments.items()}


def count_parent_variant(bam_path, reference, variant, args):
    counts = {"ref": 0, "alt": 0}
    if not bam_path:
        counts["status"] = "missing"
        return counts
    if not os.path.exists(bam_path):
        counts["status"] = "missing_file"
        return counts
    if is_cram(bam_path) and not reference:
        counts["status"] = "missing_reference"
        return counts
    reference = cached_reference(reference)
    bam = pysam.AlignmentFile(bam_path, "rc" if is_cram(bam_path) else "rb", reference_filename=reference)
    fetch_chrom = resolve_contig(bam, variant.chrom)
    for read in bam.fetch(fetch_chrom, max(0, variant.pos - 2), variant.pos + len(variant.ref) + 1):
        if read.mapping_quality < args.min_mapq or read.is_secondary or read.is_supplementary or read.is_qcfail:
            continue
        if read.is_duplicate and not args.include_duplicates:
            continue
        call = call_variant(read, variant, args.min_baseq)
        if call:
            counts["alt" if call[0] else "ref"] += 1
    bam.close()
    depth = counts["ref"] + counts["alt"]
    alt_frac = counts["alt"] / depth if depth else 0.0
    counts["depth"] = depth
    counts["alt_frac"] = round(alt_frac, 4)
    if depth < args.min_parent_bam_depth:
        counts["status"] = "low_depth"
    elif counts["alt"] >= args.min_parent_bam_alt_depth and alt_frac >= args.min_parent_bam_alt_frac:
        counts["status"] = "has_alt"
    elif counts["alt"] <= args.max_parent_bam_alt_depth_for_ref and alt_frac <= args.max_parent_bam_alt_frac_for_ref:
        counts["status"] = "no_alt"
    else:
        counts["status"] = "ambiguous"
    return counts


def parent_bam_statuses(args, variants):
    out = {"mother": {}, "father": {}}
    for i, variant in enumerate(variants):
        out["mother"][i] = count_parent_variant(args.mother_bam, args.reference, variant, args)
        out["father"][i] = count_parent_variant(args.father_bam, args.reference, variant, args)
    return out


def parent_bam_origins(statuses):
    return parent_origins(statuses)


def edge_votes(fragments):
    votes = defaultdict(lambda: [0.0, 0.0, 0])
    informative = 0
    for calls in fragments.values():
        items = sorted(calls.items())
        if len(items) < 2:
            continue
        informative += 1
        for a in range(len(items)):
            i, (allele_i, qual_i) = items[a]
            for j, (allele_j, qual_j) in items[a + 1 :]:
                relation = allele_i ^ allele_j
                weight = max(1.0, min(60, qual_i, qual_j) / 10.0)
                votes[(i, j)][relation] += weight
                votes[(i, j)][2] += 1
    return votes, informative


def best_path(votes, source, target, skip_edge=None):
    graph = defaultdict(list)
    for (i, j), (cis, trans, _) in votes.items():
        if skip_edge and (i, j) == skip_edge:
            continue
        if cis == trans:
            continue
        relation = 0 if cis > trans else 1
        support = abs(cis - trans)
        graph[i].append((j, relation, support))
        graph[j].append((i, relation, support))
    best = {(source, 0): (math.inf, [source])}
    queue = [(source, 0)]
    while queue:
        node, parity = queue.pop(0)
        score, path = best[(node, parity)]
        for nxt, edge_parity, support in graph[node]:
            new_parity = parity ^ edge_parity
            new_score = min(score, support)
            state = (nxt, new_parity)
            if new_score > best.get(state, (-1, []))[0]:
                best[state] = (new_score, path + [nxt])
                queue.append(state)
    return best.get((target, 0), (0.0, [])), best.get((target, 1), (0.0, []))


def classify_pair(variants, fragments):
    read_votes, informative = edge_votes(fragments)
    direct_read = read_votes.get((0, 1), [0.0, 0.0, 0])
    read_cis_path, read_trans_path = best_path(read_votes, 0, 1, skip_edge=(0, 1))
    direct_delta = direct_read[0] - direct_read[1]
    path_delta = read_cis_path[0] - read_trans_path[0]
    score = direct_delta + path_delta
    if abs(score) < 1.0:
        phase = "ambiguous"
    else:
        phase = "cis" if score > 0 else "trans"
    return {
        "phase": phase,
        "score": round(score, 3),
        "direct_cis_weight": round(direct_read[0], 3),
        "direct_trans_weight": round(direct_read[1], 3),
        "direct_fragments": direct_read[2],
        "best_read_cis_path_score": round(read_cis_path[0], 3),
        "best_read_trans_path_score": round(read_trans_path[0], 3),
        "best_cis_path_score": round(read_cis_path[0], 3),
        "best_trans_path_score": round(read_trans_path[0], 3),
        "best_cis_path": ",".join(variants[i].name for i in read_cis_path[1]),
        "best_trans_path": ",".join(variants[i].name for i in read_trans_path[1]),
        "informative_fragments": informative,
        "called_fragments": sum(1 for calls in fragments.values() if calls),
    }


def write_fragment_evidence(path, variants, fragments):
    if not path:
        return
    with open(path, "w", newline="") as fh:
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=["fragment", "n_variants", "calls"])
        writer.writeheader()
        for name, calls in sorted(fragments.items()):
            if len(calls) < 2:
                continue
            call_text = ";".join(f"{variants[i].name}:{allele}:q{qual}" for i, (allele, qual) in sorted(calls.items()))
            writer.writerow({"fragment": name, "n_variants": len(calls), "calls": call_text})


def add_metadata(result, args, variant1, variant2, span_start, span_end, origins, bridge_count, method):
    variant1_origin = origins.get(0, "not_tested")
    variant2_origin = origins.get(1, "not_tested")
    result.update(
        {
            "variant1": f"{variant1.chrom}:{variant1.pos}:{variant1.ref}:{variant1.alt}",
            "variant2": f"{variant2.chrom}:{variant2.pos}:{variant2.ref}:{variant2.alt}",
            "region": f"{variant1.chrom}:{max(1, span_start)}-{span_end}",
            "bam": args.bam,
            "reference": args.reference or "",
            "vcf": args.vcf or "",
            "sample": args.sample or "",
            "variant1_origin": variant1_origin,
            "variant2_origin": variant2_origin,
            "bridge_variants": bridge_count,
            "method": method,
        }
    )
    return result


def empty_parent_bam_statuses():
    blank = {"status": "not_tested", "ref": 0, "alt": 0, "depth": 0, "alt_frac": 0.0}
    return {"mother": {0: dict(blank), 1: dict(blank)}, "father": {0: dict(blank), 1: dict(blank)}}


def empty_parent_vcf_statuses():
    blank = {"status": "not_tested"}
    return {"mother": {0: dict(blank), 1: dict(blank)}, "father": {0: dict(blank), 1: dict(blank)}}


def add_parent_vcf_metadata(result, args, statuses, phase):
    result.update(
        {
            "father_vcf": args.father_vcf or "",
            "mother_vcf": args.mother_vcf or "",
            "father_sample": args.father_sample or "",
            "mother_sample": args.mother_sample or "",
            "parent_vcf_phase": phase,
            "variant1_mother_vcf_status": statuses["mother"][0]["status"],
            "variant1_father_vcf_status": statuses["father"][0]["status"],
            "variant2_mother_vcf_status": statuses["mother"][1]["status"],
            "variant2_father_vcf_status": statuses["father"][1]["status"],
        }
    )
    return result


def add_parent_bam_metadata(result, args, statuses, phase):
    result.update(
        {
            "father_bam": args.father_bam or "",
            "mother_bam": args.mother_bam or "",
            "parent_bam_phase": phase,
            "variant1_mother_bam_status": statuses["mother"][0]["status"],
            "variant1_mother_bam_ref_depth": statuses["mother"][0].get("ref", 0),
            "variant1_mother_bam_alt_depth": statuses["mother"][0].get("alt", 0),
            "variant1_mother_bam_alt_frac": statuses["mother"][0].get("alt_frac", 0.0),
            "variant1_father_bam_status": statuses["father"][0]["status"],
            "variant1_father_bam_ref_depth": statuses["father"][0].get("ref", 0),
            "variant1_father_bam_alt_depth": statuses["father"][0].get("alt", 0),
            "variant1_father_bam_alt_frac": statuses["father"][0].get("alt_frac", 0.0),
            "variant2_mother_bam_status": statuses["mother"][1]["status"],
            "variant2_mother_bam_ref_depth": statuses["mother"][1].get("ref", 0),
            "variant2_mother_bam_alt_depth": statuses["mother"][1].get("alt", 0),
            "variant2_mother_bam_alt_frac": statuses["mother"][1].get("alt_frac", 0.0),
            "variant2_father_bam_status": statuses["father"][1]["status"],
            "variant2_father_bam_ref_depth": statuses["father"][1].get("ref", 0),
            "variant2_father_bam_alt_depth": statuses["father"][1].get("alt", 0),
            "variant2_father_bam_alt_frac": statuses["father"][1].get("alt_frac", 0.0),
        }
    )
    return result


def phase_one(args, variant1, variant2):
    if variant1.chrom != variant2.chrom:
        raise ValueError("Both variants must be on the same chromosome")
    if not args.bam:
        raise ValueError("Provide --bam or a bam column in --pairs-tsv for proband read-backed phasing")
    targets = [variant1, variant2]
    span_start = min(variant1.pos, variant2.pos) - args.window
    span_end = max(variant1.pos + len(variant1.ref), variant2.pos + len(variant2.ref)) + args.window
    bridges = load_bridge_variants(args.vcf, args.sample, variant1.chrom, span_start, span_end, targets, args.max_bridges)
    variants = targets + sorted(bridges, key=lambda v: v.pos)
    fragments = read_fragments(args.bam, args.reference, variants, span_start, span_end, args.min_mapq, args.min_baseq, args.include_duplicates)
    result = classify_pair(variants, fragments)
    origins = {i: "not_tested" for i in range(len(variants))}
    method = "read_backed"
    parent_vcf = empty_parent_vcf_statuses()
    parent_vcf_phase = "not_tested"
    parent_bam = empty_parent_bam_statuses()
    parent_bam_phase = "not_tested"
    if result["phase"] == "ambiguous":
        parent_vcf = parent_vcf_statuses(args, targets)
        vcf_origins = parent_origins(parent_vcf)
        origins.update(vcf_origins)
        parent_phase = phase_from_origins(origins.get(0, "not_tested"), origins.get(1, "not_tested"))
        parent_vcf_phase = parent_phase
        if parent_phase in {"cis", "trans"}:
            result["phase"] = parent_phase
            result["score"] = args.parent_vcf_weight if parent_phase == "cis" else -args.parent_vcf_weight
            method = "parent_vcf_rescue"
        else:
            parent_bam = parent_bam_statuses(args, targets)
            for i, origin in parent_bam_origins(parent_bam).items():
                if vcf_origins.get(i) not in {"maternal", "paternal"}:
                    origins[i] = origin
            parent_phase = phase_from_origins(origins.get(0, "not_tested"), origins.get(1, "not_tested"))
            parent_bam_phase = parent_phase
            if parent_phase in {"cis", "trans"}:
                result["phase"] = parent_phase
                result["score"] = args.parent_bam_weight if parent_phase == "cis" else -args.parent_bam_weight
                method = "parent_vcf_bam_rescue" if any(vcf_origins.get(i) in {"maternal", "paternal"} for i in (0, 1)) else "parent_bam_rescue"
            else:
                method = "read_backed_parent_uninformative"
    result = add_metadata(result, args, variant1, variant2, span_start, span_end, origins, len(bridges), method)
    result = add_parent_vcf_metadata(result, args, parent_vcf, parent_vcf_phase)
    result = add_parent_bam_metadata(result, args, parent_bam, parent_bam_phase)
    return result, variants, fragments


def add_input_metadata(result, metadata):
    for key, value in metadata.items():
        out_key = key
        if out_key in result:
            out_key = f"input_{key}"
        while out_key in result:
            out_key = f"input_{out_key}"
        result[out_key] = value
    return result


def phase_one_row(item):
    row_index, args, variant1, variant2, metadata = item
    result, _, _ = phase_one(args, variant1, variant2)
    result["row_index"] = row_index
    add_input_metadata(result, metadata)
    return result


def row_value(row, key, default=None, missing_overrides=False):
    if key not in row:
        return default
    value = row.get(key)
    if value is None:
        return None if missing_overrides else default
    value = value.strip()
    if value.upper() in MISSING_VALUES:
        return None if missing_overrides else default
    return value


def pair_rows(args):
    if args.pairs_tsv:
        with open(args.pairs_tsv, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            required = {"chrom", "pos1", "ref1", "alt1", "pos2", "ref2", "alt2"}
            missing = required - set(reader.fieldnames or [])
            if missing:
                raise ValueError(f"{args.pairs_tsv} is missing columns: {', '.join(sorted(missing))}")
            metadata_columns = [key for key in reader.fieldnames if key not in PAIR_COLUMNS]
            for row_index, row in enumerate(reader, start=1):
                row_args = argparse.Namespace(**vars(args))
                row_args.bam = row_value(row, "bam", args.bam)
                row_args.reference = row_value(row, "reference", args.reference)
                row_args.vcf = row_value(row, "vcf", args.vcf)
                row_args.sample = row_value(row, "sample", args.sample)
                row_args.father_vcf = row_value(row, "father_vcf", args.father_vcf, missing_overrides=True)
                row_args.mother_vcf = row_value(row, "mother_vcf", args.mother_vcf, missing_overrides=True)
                row_args.father_sample = row_value(row, "father_sample", args.father_sample)
                row_args.mother_sample = row_value(row, "mother_sample", args.mother_sample)
                row_args.father_bam = row_value(row, "father_bam", args.father_bam, missing_overrides=True)
                row_args.mother_bam = row_value(row, "mother_bam", args.mother_bam, missing_overrides=True)
                yield (
                    row_index,
                    row_args,
                    Variant(row["chrom"], int(row["pos1"]), row["ref1"].upper(), row["alt1"].upper(), "variant1", True),
                    Variant(row["chrom"], int(row["pos2"]), row["ref2"].upper(), row["alt2"].upper(), "variant2", True),
                    {key: row.get(key, "") for key in metadata_columns},
                )
    else:
        yield 1, args, parse_variant(args.variant1, "variant1"), parse_variant(args.variant2, "variant2"), {}


def init_worker():
    global pysam, REFERENCE_CACHE
    import pysam as pysam_module

    pysam = pysam_module
    REFERENCE_CACHE = {}


def job_reference_key(job):
    _, row_args, _, _, _ = job
    return os.path.realpath(row_args.reference) if row_args.reference else ""


def main():
    global pysam
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bam", help="Coordinate-sorted, indexed proband BAM or CRAM; required unless pairs TSV has bam")
    parser.add_argument("--reference", help="Reference FASTA; required for CRAM and recommended for indels; can be overridden per batch row")
    parser.add_argument("--variant1", help="First target as chrom:pos:ref:alt")
    parser.add_argument("--variant2", help="Second target as chrom:pos:ref:alt")
    parser.add_argument("--pairs-tsv", help="Batch TSV with target variants and optional sample/proband/parent input columns")
    parser.add_argument("--vcf", help="Optional indexed VCF/BCF of nearby heterozygous bridge variants; can be overridden per batch row")
    parser.add_argument("--sample", help="Sample name in VCF; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--father-vcf", help="Optional father VCF/gVCF for fast parent rescue; can be overridden per batch row")
    parser.add_argument("--mother-vcf", help="Optional mother VCF/gVCF for fast parent rescue; can be overridden per batch row")
    parser.add_argument("--father-sample", help="Sample name in father VCF/gVCF; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--mother-sample", help="Sample name in mother VCF/gVCF; can be overridden per batch row; defaults to first sample")
    parser.add_argument("--father-bam", help="Optional father BAM/CRAM for parent-BAM rescue; can be overridden per batch row")
    parser.add_argument("--mother-bam", help="Optional mother BAM/CRAM for parent-BAM rescue; can be overridden per batch row")
    parser.add_argument("--min-parent-bam-depth", type=int, default=8, help="Minimum parent BAM depth to call no_alt or has_alt")
    parser.add_argument("--min-parent-bam-alt-depth", type=int, default=3, help="Minimum ALT-supporting parent reads for has_alt")
    parser.add_argument("--min-parent-bam-alt-frac", type=float, default=0.2, help="Minimum parent ALT read fraction for has_alt")
    parser.add_argument("--max-parent-bam-alt-depth-for-ref", type=int, default=1, help="Maximum ALT-supporting parent reads for no_alt")
    parser.add_argument("--max-parent-bam-alt-frac-for-ref", type=float, default=0.05, help="Maximum parent ALT read fraction for no_alt")
    parser.add_argument("--parent-vcf-weight", type=float, default=10.0, help="Score assigned when parent VCF/gVCF rescue phases the target pair")
    parser.add_argument("--parent-bam-weight", type=float, default=10.0, help="Score assigned when parent BAM rescue phases the target pair")
    parser.add_argument("--window", type=int, default=1000, help="Bases to fetch around target pair")
    parser.add_argument("--max-bridges", type=int, default=200, help="Maximum nearby heterozygous variants to use")
    parser.add_argument("--min-mapq", type=int, default=20)
    parser.add_argument("--min-baseq", type=int, default=20)
    parser.add_argument("--include-duplicates", action="store_true")
    parser.add_argument("--threads", type=int, default=1, help="Worker processes for batch mode")
    parser.add_argument("--chunksize", type=int, default=1, help="Pairs assigned per worker task in batch mode")
    parser.add_argument("--out", required=True, help="Output result TSV")
    parser.add_argument("--evidence", help="Optional fragment evidence TSV; only valid for one pair")
    parser.add_argument("--json", help="Optional JSON output; only valid for one pair")
    args = parser.parse_args()

    try:
        import pysam as pysam_module
    except ImportError as exc:
        raise SystemExit("This script requires pysam. Install with: pip install pysam") from exc
    pysam = pysam_module

    if bool(args.pairs_tsv) == bool(args.variant1 or args.variant2):
        raise SystemExit("Provide either --pairs-tsv or both --variant1 and --variant2")
    if not args.pairs_tsv and not (args.variant1 and args.variant2):
        raise SystemExit("Provide both --variant1 and --variant2")
    if args.pairs_tsv and (args.evidence or args.json):
        raise SystemExit("--evidence and --json are only supported for single-pair runs")
    if args.threads < 1 or args.chunksize < 1:
        raise SystemExit("--threads and --chunksize must be at least 1")

    rows = []
    last = None
    jobs = list(pair_rows(args))
    if args.pairs_tsv:
        jobs = sorted(jobs, key=job_reference_key)
    if args.pairs_tsv and args.threads > 1:
        with mp.Pool(args.threads, initializer=init_worker) as pool:
            rows = list(pool.imap(phase_one_row, jobs, chunksize=args.chunksize))
    else:
        for row_index, row_args, variant1, variant2, metadata in jobs:
            result, variants, fragments = phase_one(row_args, variant1, variant2)
            result["row_index"] = row_index
            add_input_metadata(result, metadata)
            rows.append(result)
            last = result, variants, fragments
    if args.pairs_tsv:
        rows.sort(key=lambda row: row["row_index"])

    with open(args.out, "w", newline="") as fh:
        fieldnames = list(rows[0].keys()) if rows else []
        writer = csv.DictWriter(fh, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)

    if last and not args.pairs_tsv:
        result, variants, fragments = last
        write_fragment_evidence(args.evidence, variants, fragments)
        if args.json:
            with open(args.json, "w") as fh:
                json.dump(result, fh, indent=2)


if __name__ == "__main__":
    main()
