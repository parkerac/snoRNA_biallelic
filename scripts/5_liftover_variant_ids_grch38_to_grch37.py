#!/usr/bin/env python3
"""Lift chr:pos:ref:alt variants from GRCh38 to GRCh37."""

import argparse
import csv
import os
import shutil
import subprocess
import tempfile


def open_text(path, mode="rt"):
    return open(path, mode) if str(path) == "-" else open(path, mode)


def revcomp(seq):
    return seq.translate(str.maketrans("ACGTacgt", "TGCAtgca"))[::-1]


def parse_variant(text):
    chrom, pos, ref, alt = text.rstrip("\n").split(":")
    return chrom, int(pos), ref, alt


def read_variants(path, column):
    if column:
        with open(path, newline="") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            if column not in (reader.fieldnames or []):
                raise ValueError(f"{path} is missing column: {column}")
            rows = list(reader)
            variants = [parse_variant(row[column]) if row.get(column) else None for row in rows]
            return rows, variants
    with open(path) as fh:
        variants = [parse_variant(line) for line in fh if line.strip() and not line.startswith("#")]
    return None, variants


class FastaIndex:
    def __init__(self, fasta_path):
        self.fasta_path = fasta_path
        self.records = {}
        with open(fasta_path + ".fai") as fh:
            for line in fh:
                name, length, offset, line_bases, line_width = line.rstrip("\n").split("\t")[:5]
                self.records[name] = (int(length), int(offset), int(line_bases), int(line_width))
        self.handle = open(fasta_path, "rb")

    def fetch(self, chrom, start, end):
        length, offset, line_bases, line_width = self.records[chrom]
        if start < 0 or end > length or start >= end:
            return ""
        first = start // line_bases
        last = (end - 1) // line_bases
        parts = []
        for line in range(first, last + 1):
            line_start = line * line_bases
            line_end = min(line_start + line_bases, end)
            chunk_start = max(start, line_start)
            chunk_len = line_end - chunk_start
            self.handle.seek(offset + line * line_width + (chunk_start - line_start))
            parts.append(self.handle.read(chunk_len).decode("ascii"))
        return "".join(parts)

    def close(self):
        self.handle.close()


def chrom_style(chrom, mode):
    if mode == "keep":
        return chrom
    if mode == "add_chr":
        return chrom if chrom.startswith("chr") else f"chr{chrom}"
    if mode == "strip_chr":
        return chrom[3:] if chrom.startswith("chr") else chrom
    raise ValueError(f"Unknown chrom mode: {mode}")


def liftover_variants(variants, chain_path, chrom_mode):
    if not shutil.which("liftOver"):
        raise RuntimeError("liftOver is not on PATH")

    with tempfile.TemporaryDirectory() as tmpdir:
        in_bed = os.path.join(tmpdir, "in.bed")
        out_bed = os.path.join(tmpdir, "out.bed")
        unmapped = os.path.join(tmpdir, "unmapped.bed")

        with open(in_bed, "w") as fh:
            for idx, (chrom, pos, ref, _) in enumerate(variants):
                chrom = chrom_style(chrom, chrom_mode)
                fh.write(f"{chrom}\t{pos - 1}\t{pos - 1 + len(ref)}\t{idx}\n")

        subprocess.run(["liftOver", in_bed, chain_path, out_bed, unmapped], check=True)

        mapped = {}
        with open(out_bed) as fh:
            for line in fh:
                chrom, start, end, idx = line.rstrip("\n").split("\t")[:4]
                mapped.setdefault(int(idx), []).append((chrom, int(start), int(end)))
        return mapped


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--variants", required=True, help="Text file with one chr:pos:ref:alt per line")
    parser.add_argument("--chain", required=True, help="GRCh38-to-GRCh37 UCSC over.chain file")
    parser.add_argument("--fasta", required=True, help="GRCh37 FASTA with .fai index")
    parser.add_argument("--out", required=True, help="Output file with lifted chr:pos:ref:alt variants")
    parser.add_argument("--column", help="TSV column containing chr:pos:ref:alt variant IDs")
    parser.add_argument("--lifted-column", default="lifted_variant_id", help="Output TSV column for lifted chr:pos:ref:alt values")
    parser.add_argument("--unmapped", help="Optional file for variants that do not lift cleanly")
    parser.add_argument("--chrom-mode", choices=["keep", "add_chr", "strip_chr"], default="keep")
    args = parser.parse_args()

    rows, variants = read_variants(args.variants, args.column)

    lifted = liftover_variants(variants, args.chain, args.chrom_mode)
    fasta = FastaIndex(args.fasta)
    unmapped_rows = []

    def lift_one(idx, chrom, pos, ref, alt):
        if chrom is None:
            return None, "missing_variant_id"
        hits = lifted.get(idx, [])
        if len(hits) != 1:
            return None, "no_unique_liftover"
        tchrom, tstart, tend = hits[0]
        target_ref = fasta.fetch(tchrom, tstart, tend).upper()
        ref_u = ref.upper()
        alt_u = alt.upper()
        if target_ref == ref_u:
            return f"{tchrom}:{tstart + 1}:{ref_u}:{alt_u}", ""
        if target_ref == revcomp(ref_u).upper():
            return f"{tchrom}:{tstart + 1}:{revcomp(ref_u)}:{revcomp(alt_u)}", ""
        return None, "reference_mismatch"

    if rows is None:
        with open(args.out, "w") as out:
            for idx, (chrom, pos, ref, alt) in enumerate(variants):
                lifted_variant, reason = lift_one(idx, chrom, pos, ref, alt)
                if lifted_variant is None:
                    unmapped_rows.append((chrom, pos, ref, alt, reason))
                    continue
                out.write(f"{lifted_variant}\n")
    else:
        fieldnames = list((rows[0].keys() if rows else []))
        if args.lifted_column not in fieldnames:
            fieldnames.append(args.lifted_column)
        with open(args.out, "w", newline="") as out:
            writer = csv.DictWriter(out, delimiter="\t", fieldnames=fieldnames, extrasaction="ignore")
            writer.writeheader()
            for idx, row in enumerate(rows):
                chrom, pos, ref, alt = variants[idx]
                lifted_variant, reason = lift_one(idx, chrom, pos, ref, alt)
                row = dict(row)
                row[args.lifted_column] = lifted_variant or ""
                writer.writerow(row)
                if lifted_variant is None:
                    unmapped_rows.append((chrom, pos, ref, alt, reason))

    fasta.close()

    if args.unmapped:
        with open(args.unmapped, "w") as fh:
            fh.write("chrom\tpos\tref\talt\treason\n")
            for row in unmapped_rows:
                fh.write("\t".join(map(str, row)) + "\n")


if __name__ == "__main__":
    main()
