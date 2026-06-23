#!/usr/bin/env python3
"""
GFA2FASTA.py — Convert GFA format to FASTA.

Extracts segment sequences (S lines) from GFA version 1 or 2 files and
writes them to a FASTA file. The GFA version is auto-detected from the
header; if no version header is present, the format is inferred from
the structure of the first S line. Segments with no sequence ('*') are
skipped with a warning.

GFA1 S line:  S  <name>  <sequence>  [optional tags]
GFA2 S line:  S  <name>  <length>    <sequence>[+-]  [optional tags]

Usage
-----
    GFA2FASTA.py --input assembly.gfa --output assembly.fasta
    GFA2FASTA.py --input assembly.gfa --output assembly.fasta --summary summary.txt
"""

import argparse
import sys
from datetime import datetime
from pathlib import Path

VERSION = "v0.1.0"

FASTA_LINE_WIDTH = 60


# ── GFA parsing ────────────────────────────────────────────────────────────────

def detect_gfa_version(gfa_path: Path) -> str:
    """
    Return '1' or '2' by inspecting the H header line.
    Falls back to auto-detection from the first S line if no header is found.
    """
    with open(gfa_path) as fh:
        for line in fh:
            if line.startswith("H"):
                for field in line.rstrip("\n").split("\t")[1:]:
                    if field.startswith("VN:Z:"):
                        return "2" if field[5:].startswith("2") else "1"
            elif line.startswith("S"):
                parts = line.rstrip("\n").split("\t")
                # GFA2: field[2] is a numeric length, field[3] is the sequence
                if len(parts) >= 4 and parts[2].isdigit():
                    return "2"
                return "1"
    return "1"


def parse_gfa(gfa_path: Path, version: str) -> tuple:
    """
    Extract segment sequences from a GFA file.
    Returns (segments, no_seq_count) where segments is a list of
    (name, sequence) tuples and no_seq_count is the number of segments
    skipped due to missing sequence ('*').
    """
    segments     = []
    no_seq_count = 0

    with open(gfa_path) as fh:
        for lineno, line in enumerate(fh, 1):
            if not line.startswith("S"):
                continue
            parts = line.rstrip("\n").split("\t")

            if version == "2":
                if len(parts) < 4:
                    print(f"WARNING: malformed GFA2 S line at {lineno}, skipping",
                          file=sys.stderr)
                    continue
                name = parts[1]
                seq  = parts[3].rstrip("+-")  # strip orientation marker
            else:
                if len(parts) < 3:
                    print(f"WARNING: malformed GFA1 S line at {lineno}, skipping",
                          file=sys.stderr)
                    continue
                name = parts[1]
                seq  = parts[2]

            if seq == "*":
                print(f"WARNING: segment '{name}' has no sequence ('*'), skipping",
                      file=sys.stderr)
                no_seq_count += 1
                continue

            segments.append((name, seq))

    return segments, no_seq_count


# ── Writers ───────────────────────────────────────────────────────────────────

def write_fasta(segments: list, path: Path) -> None:
    with open(path, "w") as fh:
        for name, seq in segments:
            fh.write(f">{name}\n")
            for i in range(0, len(seq), FASTA_LINE_WIDTH):
                fh.write(seq[i:i + FASTA_LINE_WIDTH] + "\n")


def write_summary(segments: list, no_seq: int, version: str,
                  gfa_path: Path, fasta_path: Path, path: Path) -> None:
    total_len = sum(len(seq) for _, seq in segments)
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    col = 36
    with open(path, "w") as fh:
        fh.write(f"{'Date':<{col}}: {timestamp}\n")
        fh.write(f"{'Input GFA':<{col}}: {gfa_path}\n")
        fh.write(f"{'Output FASTA':<{col}}: {fasta_path}\n")
        fh.write(f"{'GFA version detected':<{col}}: {version}\n")
        fh.write(f"{'Segments with sequence (written)':<{col}}: {len(segments)}\n")
        fh.write(f"{'Segments without sequence (*)':<{col}}: {no_seq}\n")
        fh.write(f"{'Total length (bp)':<{col}}: {total_len}\n")


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GFA2FASTA",
        description="Convert GFA format (v1 or v2) to FASTA.",
    )
    ap.add_argument("--input",   required=True, type=Path,
                    help="Input GFA file")
    ap.add_argument("--output",  required=True, type=Path,
                    help="Output FASTA file")
    ap.add_argument("--summary", type=Path, default=None,
                    help="Optional summary report (plain text)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.input.exists():
        print(f"ERROR: --input file not found: {args.input}", file=sys.stderr)
        sys.exit(1)

    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"GFA2FASTA.py  {timestamp}", file=sys.stderr)

    version = detect_gfa_version(args.input)
    print(f"GFA version  : {version}", file=sys.stderr)

    segments, no_seq = parse_gfa(args.input, version)
    total_len = sum(len(seq) for _, seq in segments)

    print(f"Segments written  : {len(segments)}", file=sys.stderr)
    print(f"Segments skipped  : {no_seq}  (no sequence)", file=sys.stderr)
    print(f"Total length (bp) : {total_len}", file=sys.stderr)

    write_fasta(segments, args.output)
    print(f"Written: {args.output}", file=sys.stderr)

    if args.summary:
        write_summary(segments, no_seq, version, args.input, args.output, args.summary)
        print(f"Summary: {args.summary}", file=sys.stderr)


if __name__ == "__main__":
    main()
