#!/usr/bin/env python3
"""
GFF2BEDOrthoVenn.py — Convert a GFF3 file to the 5-column BED format
expected by OrthoVennPlus (https://orthovenn3.bioinfotoolkits.net/).

For every feature whose column-3 type matches --feature_type (default:
gene), writes one line:

    SeqID    GeneID    Start    End    Strand

GeneID is taken from the ID= attribute. Coordinates are copied as-is from
the GFF3 (1-based, inclusive) -- OrthoVennPlus expects this exact layout,
not 0-based BED. Rows are sorted by (SeqID, Start).

Usage
-----
    GFF2BEDOrthoVenn.py --gff annotation.gff3 --output annotation.bed
    GFF2BEDOrthoVenn.py --gff annotation.gff3 --dry_run
"""

import argparse
import re
import sys
from pathlib import Path

VERSION = "v0.0.1"


def parse_attributes(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        key, val = field.split("=", 1)
        attrs[key.strip()] = val.strip()
    return attrs


def natural_key(s: str) -> list:
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r"(\d+)", s)]


def extract_rows(path: Path, feature_type: str) -> tuple:
    """Return (rows, n_no_id). rows: list of (seqid, gene_id, start, end, strand)."""
    rows = []
    n_no_id = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            raw = line.rstrip("\n")
            if not raw.strip() or raw.startswith("#"):
                continue
            cols = raw.split("\t")
            if len(cols) < 9 or cols[2] != feature_type:
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                print(f"WARNING: skipping GFF3 line {lineno} "
                      f"(non-integer start/end)", file=sys.stderr)
                continue
            attrs = parse_attributes(cols[8])
            gene_id = attrs.get("ID")
            if not gene_id:
                n_no_id += 1
                continue
            rows.append((cols[0], gene_id, start, end, cols[6]))
    return rows, n_no_id


def write_bed(rows: list, out_fh) -> None:
    rows_sorted = sorted(rows, key=lambda r: (natural_key(r[0]), r[2]))
    for seqid, gene_id, start, end, strand in rows_sorted:
        out_fh.write(f"{seqid}\t{gene_id}\t{start}\t{end}\t{strand}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GFF2BEDOrthoVenn",
        description="Convert a GFF3 file to the 5-column BED format "
                     "expected by OrthoVennPlus.",
    )
    ap.add_argument("--gff", required=True, type=Path,
                    help="Input GFF3 file")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output BED file (default: stdout)")
    ap.add_argument("--feature_type", default="gene",
                    help="Column-3 feature type to extract (default: gene)")
    ap.add_argument("--dry_run", action="store_true",
                    help="Parse and report counts, then exit without "
                         "writing any output")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.gff.exists():
        print(f"ERROR: --gff file not found: {args.gff}", file=sys.stderr)
        sys.exit(1)

    rows, n_no_id = extract_rows(args.gff, args.feature_type)
    print(f"Loaded {len(rows)} '{args.feature_type}' feature(s) from "
          f"{args.gff.name}", file=sys.stderr)
    if n_no_id:
        print(f"WARNING: {n_no_id} '{args.feature_type}' feature(s) had no "
              f"ID attribute and were skipped", file=sys.stderr)
    if not rows:
        print(f"ERROR: no features of type '{args.feature_type}' with an "
              f"ID= attribute found. Nothing to write.", file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print(f"  Feature type : {args.feature_type}", file=sys.stderr)
        print(f"  Rows to write: {len(rows)}", file=sys.stderr)
        print("  Exiting (--dry_run). No output written.", file=sys.stderr)
        return

    if args.output:
        with open(args.output, "w") as out_fh:
            write_bed(rows, out_fh)
        print(f"Written: {args.output}", file=sys.stderr)
    else:
        write_bed(rows, sys.stdout)


if __name__ == "__main__":
    main()
