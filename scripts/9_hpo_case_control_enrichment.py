#!/usr/bin/env python3
"""Test HPO term enrichment between case and control participants."""

import argparse
import csv
import math
import os
import re
import sys
from collections import defaultdict


DEFAULT_CASE_VALUE = "case"
DEFAULT_CONTROL_VALUE = "control"
DEFAULT_MIN_DEPTH = 2
DEFAULT_MIN_CASE_COUNT = 2
DEFAULT_MIN_TOTAL_COUNT = 3


def require_columns(fieldnames, required, label):
    missing = sorted(set(required) - set(fieldnames or []))
    if missing:
        raise ValueError(f"{label} is missing required columns: {', '.join(missing)}")


def split_hpo_terms(value):
    if value is None:
        return []
    return [item.strip() for item in re.split(r"[;,]", str(value)) if item.strip()]


def normalize_term_id(hpo, value):
    if not value:
        return None
    try:
        term_id = hpo.get_term(value)
    except Exception:
        term_id = None
    return term_id


def term_id_value(node):
    identifier = getattr(node, "identifier", node)
    value = getattr(identifier, "value", None)
    return value if value else str(identifier)


def term_name(node):
    name = getattr(node, "name", None)
    return name if name else ""


def load_hpo(args):
    try:
        import hpotk
    except ImportError as exc:
        raise SystemExit("This script requires hpo-toolkit. Install with: pip install hpo-toolkit") from exc

    if args.hpo_json:
        if not os.path.exists(args.hpo_json):
            raise SystemExit(f"HPO JSON not found: {args.hpo_json}")
        return hpotk.load_minimal_ontology(args.hpo_json)

    store = hpotk.configure_ontology_store()
    return store.load_hpo(args.hpo_release)


def load_participants(path, participant_col, case_col, hpo_col, case_value, control_value):
    participants = {}
    groups = {}
    with open(path, newline="") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        require_columns(reader.fieldnames, {participant_col, case_col, hpo_col}, "Input TSV")
        for row in reader:
            participant_id = row.get(participant_col, "").strip()
            if not participant_id:
                continue
            group = row.get(case_col, "").strip().lower()
            if group not in {case_value, control_value}:
                raise ValueError(f"Unexpected group value {group!r} for participant {participant_id!r}")
            groups.setdefault(participant_id, group)
            if groups[participant_id] != group:
                raise ValueError(f"Participant {participant_id!r} has conflicting group values")
            participants.setdefault(participant_id, set()).update(split_hpo_terms(row.get(hpo_col)))
    return groups, participants


def expand_with_ancestors(hpo, terms):
    expanded = set()
    for raw_term in terms:
        term = normalize_term_id(hpo, raw_term)
        if term is None:
            continue
        for anc in hpo.graph.get_ancestors(term, include_source=True):
            expanded.add(term_id_value(anc))
    return expanded


def depth_to_root(hpo, term_id, root_id, cache):
    if term_id == root_id:
        return 0
    if term_id in cache:
        return cache[term_id]
    parents = list(hpo.graph.get_parents(term_id))
    if not parents:
        cache[term_id] = None
        return None
    depths = [depth_to_root(hpo, term_id_value(parent), root_id, cache) for parent in parents]
    depths = [depth for depth in depths if depth is not None]
    cache[term_id] = (1 + min(depths)) if depths else None
    return cache[term_id]


def log_comb(n, k):
    if k < 0 or k > n:
        return float("-inf")
    return math.lgamma(n + 1) - math.lgamma(k + 1) - math.lgamma(n - k + 1)


def fisher_exact_greater(a, b, c, d):
    row1 = a + b
    row2 = c + d
    col1 = a + c
    total = row1 + row2
    max_x = min(row1, col1)
    if a > max_x:
        return 0.0
    denom = log_comb(total, row1)

    def log_p(x):
        return log_comb(col1, x) + log_comb(total - col1, row1 - x) - denom

    logs = [log_p(x) for x in range(a, max_x + 1)]
    if not logs:
        return 1.0
    peak = max(logs)
    return sum(math.exp(value - peak) for value in logs) * math.exp(peak)


def odds_ratio(a, b, c, d, pseudocount=0.5):
    return ((a + pseudocount) * (d + pseudocount)) / ((b + pseudocount) * (c + pseudocount))


def benjamini_hochberg(pvalues):
    indexed = sorted(enumerate(pvalues), key=lambda item: item[1])
    qvalues = [1.0] * len(pvalues)
    prev = 1.0
    for rank in range(len(indexed), 0, -1):
        idx, pvalue = indexed[rank - 1]
        q = min(prev, (pvalue * len(pvalues)) / rank)
        qvalues[idx] = q
        prev = q
    return qvalues


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-tsv", required=True, help="TSV with participant_id, case, and hpo_ids columns")
    parser.add_argument("--participant-col", default="participant_id", help="Participant ID column")
    parser.add_argument("--case-col", default="case", help="Case/control label column")
    parser.add_argument("--hpo-col", default="hpo_ids", help="Comma- or semicolon-separated HPO IDs column")
    parser.add_argument("--case-value", default=DEFAULT_CASE_VALUE, help="Value indicating a case")
    parser.add_argument("--control-value", default=DEFAULT_CONTROL_VALUE, help="Value indicating a control")
    parser.add_argument("--out", required=True, help="Output TSV path")
    parser.add_argument("--hpo-json", help="Optional local HPO JSON file to load instead of downloading the ontology")
    parser.add_argument("--hpo-release", help="Optional HPO release tag to load with hpo-toolkit")
    parser.add_argument("--min-depth", type=int, default=DEFAULT_MIN_DEPTH, help="Minimum depth from the HPO root to test a term")
    parser.add_argument("--min-case-count", type=int, default=DEFAULT_MIN_CASE_COUNT, help="Minimum number of cases with a term")
    parser.add_argument("--min-total-count", type=int, default=DEFAULT_MIN_TOTAL_COUNT, help="Minimum total participants with a term")
    parser.add_argument("--pseudocount", type=float, default=0.5, help="Continuity correction for odds ratios")
    args = parser.parse_args()

    if not os.path.exists(args.input_tsv):
        raise SystemExit(f"Input TSV not found: {args.input_tsv}")

    hpo = load_hpo(args)
    root_id = term_id_value(hpo.graph.root)
    print(f"Loaded HPO ontology with root {root_id}", file=sys.stderr, flush=True)

    groups, participant_terms = load_participants(
        args.input_tsv,
        participant_col=args.participant_col,
        case_col=args.case_col,
        hpo_col=args.hpo_col,
        case_value=args.case_value.lower(),
        control_value=args.control_value.lower(),
    )
    cases = {pid for pid, group in groups.items() if group == args.case_value.lower()}
    controls = {pid for pid, group in groups.items() if group == args.control_value.lower()}
    if not cases or not controls:
        raise SystemExit("Need at least one case and one control participant")
    print(f"Loaded {len(cases)} cases and {len(controls)} controls", file=sys.stderr, flush=True)

    depth_cache = {}
    expanded = {}
    for participant_id, terms in participant_terms.items():
        expanded[participant_id] = expand_with_ancestors(hpo, terms)
    print(f"Expanded HPO annotations for {len(expanded)} participants", file=sys.stderr, flush=True)

    counts = defaultdict(lambda: {"case": 0, "control": 0})
    for participant_id in cases:
        for term in expanded.get(participant_id, set()):
            counts[term]["case"] += 1
    for participant_id in controls:
        for term in expanded.get(participant_id, set()):
            counts[term]["control"] += 1
    print(f"Collected counts for {len(counts)} HPO terms before filtering", file=sys.stderr, flush=True)

    rows = []
    for term_id, term_counts in counts.items():
        case_with = term_counts["case"]
        control_with = term_counts["control"]
        total_with = case_with + control_with
        if total_with < args.min_total_count:
            continue
        if case_with < args.min_case_count:
            continue
        if case_with <= control_with:
            continue
        depth = depth_to_root(hpo, term_id, root_id, depth_cache)
        if depth is None or depth < args.min_depth:
            continue
        if term_id == root_id:
            continue

        case_without = len(cases) - case_with
        control_without = len(controls) - control_with
        pvalue = fisher_exact_greater(case_with, case_without, control_with, control_without)
        rows.append({
            "hpo_id": term_id,
            "hpo_name": term_name(hpo.get_term(term_id)),
            "depth": depth,
            "cases_with_term": case_with,
            "cases_without_term": case_without,
            "controls_with_term": control_with,
            "controls_without_term": control_without,
            "case_fraction": case_with / len(cases),
            "control_fraction": control_with / len(controls),
            "fraction_difference": (case_with / len(cases)) - (control_with / len(controls)),
            "odds_ratio": odds_ratio(case_with, case_without, control_with, control_without, pseudocount=args.pseudocount),
            "p_value": pvalue,
        })

    print(f"{len(rows)} HPO terms passed the case-enrichment filters", file=sys.stderr, flush=True)

    if rows:
        qvalues = benjamini_hochberg([row["p_value"] for row in rows])
        for row, qvalue in zip(rows, qvalues):
            row["q_value"] = qvalue
        rows.sort(key=lambda row: (row["q_value"], row["p_value"], -row["fraction_difference"], row["hpo_id"]))
    else:
        rows = []

    fieldnames = [
        "hpo_id",
        "hpo_name",
        "depth",
        "cases_with_term",
        "cases_without_term",
        "controls_with_term",
        "controls_without_term",
        "case_fraction",
        "control_fraction",
        "fraction_difference",
        "odds_ratio",
        "p_value",
        "q_value",
    ]
    with open(args.out, "w", newline="") as out:
        writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames)
        writer.writeheader()
        for row in rows:
            writer.writerow({
                "hpo_id": row["hpo_id"],
                "hpo_name": row["hpo_name"],
                "depth": row["depth"],
                "cases_with_term": row["cases_with_term"],
                "cases_without_term": row["cases_without_term"],
                "controls_with_term": row["controls_with_term"],
                "controls_without_term": row["controls_without_term"],
                "case_fraction": f"{row['case_fraction']:.6g}",
                "control_fraction": f"{row['control_fraction']:.6g}",
                "fraction_difference": f"{row['fraction_difference']:.6g}",
                "odds_ratio": f"{row['odds_ratio']:.6g}",
                "p_value": f"{row['p_value']:.6g}",
                "q_value": f"{row['q_value']:.6g}",
            })
    print(f"Wrote {len(rows)} enriched HPO terms to {args.out}", file=sys.stderr, flush=True)


if __name__ == "__main__":
    main()
