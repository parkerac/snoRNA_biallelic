#!/usr/bin/env python3
"""Prepare hg38 tRNA domain intervals from GtRNAdb/tRNAscan-SE files."""

import argparse
import csv
import os

SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
PROJECTS_DIR = os.path.dirname(os.path.dirname(SCRIPT_DIR))
DEFAULT_TRNA_DIR = os.path.join(PROJECTS_DIR, "ncRNA_biallelic_project_files", "tRNA", "hg38-tRNAs")

DOMAINS = (
    ("acceptor_stem_5p", 1, 7),
    ("d_arm", 8, 26),
    ("anticodon_arm", 27, 43),
)


def read_name_map(path):
    with open(path) as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        return {row["tRNAscan-SE_id"]: row["GtRNAdb_id"] for row in reader}


def read_bed(path):
    out = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            chrom, start0, end, name, score, strand = parts[:6]
            out[name] = {
                "chrom": chrom,
                "start0": int(start0),
                "end": int(end),
                "score": score,
                "strand": strand,
                "mature_len": sum(int(x) for x in parts[10].rstrip(",").split(",") if x) if len(parts) >= 12 else int(end) - int(start0),
            }
    return out


def read_trnascan_out(path):
    records = {}
    with open(path) as fh:
        for line in fh:
            if not line.strip() or line.startswith(("-", "Sequence", "Name")):
                continue
            parts = line.rstrip("\n").split()
            if len(parts) < 15:
                continue
            intron_begin, intron_end = int(parts[6]), int(parts[7])
            begin, end = int(parts[2]), int(parts[3])
            strand = "+" if begin <= end else "-"
            intron = None
            if intron_begin and intron_end:
                if strand == "+":
                    intron = (intron_begin - begin + 1, intron_end - begin + 1)
                else:
                    intron = (begin - intron_begin + 1, begin - intron_end + 1)
                intron = tuple(sorted(intron))
            records[f"{parts[0]}.trna{parts[1]}"] = {
                "amino_acid": parts[4],
                "anticodon": parts[5],
                "score": parts[8],
                "origin": parts[14],
                "note": " ".join(parts[15:]) if len(parts) > 15 else "",
                "intron": intron,
            }
    return records


def mature_domains(mature_len):
    variable_end = mature_len - 25
    return list(DOMAINS) + [
        ("variable_region", 44, variable_end),
        ("t_arm", variable_end + 1, mature_len - 8),
        ("acceptor_stem_3p", mature_len - 7, mature_len),
        ("anticodon", 34, 36),
    ]


def mature_to_pre_segments(start, end, intron):
    positions = []
    intron_len = 0
    intron_start = None
    if intron:
        intron_start, intron_end = intron
        intron_len = intron_end - intron_start + 1
    for pos in range(start, end + 1):
        positions.append(pos + intron_len if intron_start and pos >= intron_start else pos)
    segments = []
    seg_start = prev = positions[0]
    for pos in positions[1:]:
        if pos == prev + 1:
            prev = pos
        else:
            segments.append((seg_start, prev))
            seg_start = prev = pos
    segments.append((seg_start, prev))
    return segments


def pre_to_bed(chrom_start0, chrom_end, strand, start, end):
    if strand == "+":
        return chrom_start0 + start - 1, chrom_start0 + end
    return chrom_end - end, chrom_end - start + 1


def wanted(record, keep_notes):
    return not keep_notes or any(note in record["note"] for note in keep_notes)


FIELDS = ["chrom", "start", "end", "name", "score", "strand", "trna_id", "trnascan_id", "amino_acid", "anticodon", "domain", "mature_start", "mature_end", "pretrna_start", "pretrna_end", "trnascan_score", "origin", "note"]


def write_domains(bed_path, name_map_path, trnascan_path, out_path, keep_notes, header):
    names = read_name_map(name_map_path)
    genes = read_bed(bed_path)
    scan = read_trnascan_out(trnascan_path)
    rows = []
    for scan_id, trna_id in names.items():
        if trna_id not in genes or scan_id not in scan or not wanted(scan[scan_id], keep_notes):
            continue
        gene = genes[trna_id]
        rec = scan[scan_id]
        for domain, start, end in mature_domains(gene["mature_len"]):
            if start > end:
                continue
            for part, (pre_start, pre_end) in enumerate(mature_to_pre_segments(start, end, rec["intron"]), start=1):
                start0, bed_end = pre_to_bed(gene["start0"], gene["end"], gene["strand"], pre_start, pre_end)
                rows.append((gene["chrom"], start0, bed_end, f"{trna_id}|{domain}|part{part}", gene["score"], gene["strand"], trna_id, scan_id, rec["amino_acid"], rec["anticodon"], domain, start, end, pre_start, pre_end, rec["score"], rec["origin"], rec["note"]))
        if rec["intron"]:
            start0, bed_end = pre_to_bed(gene["start0"], gene["end"], gene["strand"], *rec["intron"])
            rows.append((gene["chrom"], start0, bed_end, f"{trna_id}|intron|part1", gene["score"], gene["strand"], trna_id, scan_id, rec["amino_acid"], rec["anticodon"], "intron", "", "", rec["intron"][0], rec["intron"][1], rec["score"], rec["origin"], rec["note"]))

    rows.sort(key=lambda r: (r[0], r[1], r[2], r[3]))
    os.makedirs(os.path.dirname(out_path) or ".", exist_ok=True)
    with open(out_path, "w", newline="") as fh:
        writer = csv.writer(fh, delimiter="\t")
        if header:
            writer.writerow(FIELDS)
        writer.writerows(rows)
    return len(rows), len({row[6] for row in rows})


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--trna-dir", default=DEFAULT_TRNA_DIR)
    parser.add_argument("--out", default=os.path.join(DEFAULT_TRNA_DIR, "hg38-tRNA-domains.bed"))
    parser.add_argument("--include-all", action="store_true", help="Include secondary filtered and pseudo tRNAs as well as high-confidence tRNAs")
    parser.add_argument("--header", action="store_true", help="Write a header row; omit this for bedtools-compatible BED")
    args = parser.parse_args()
    keep_notes = [] if args.include_all else ["high confidence set"]
    n_rows, n_genes = write_domains(
        os.path.join(args.trna_dir, "hg38-tRNAs.bed"),
        os.path.join(args.trna_dir, "hg38-tRNAs_name_map.txt"),
        os.path.join(args.trna_dir, "hg38-tRNAs-detailed.out"),
        args.out,
        keep_notes,
        args.header,
    )
    print(f"Wrote {n_rows} domain intervals for {n_genes} tRNAs to {args.out}")


if __name__ == "__main__":
    main()
