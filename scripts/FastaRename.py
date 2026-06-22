#!/usr/bin/env python3
"""
FastaRename.py — Rename sequence IDs in a FASTA file using a TSV mapping.

The TSV file must have two tab-separated columns (no header required):
    old_seqid <TAB> new_seqid

Lines starting with '#' are treated as comments and ignored.
SeqIDs not found in the mapping are kept unchanged (with a warning).
Only the ID portion of the header (up to the first space) is renamed;
any description text is preserved.

Usage
-----
    FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv
    FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv --output renamed.fasta
"""

import argparse
import sys
from pathlib import Path

VERSION = "v0.0.1"

FASTA_LINE_WIDTH = 60


def load_mapping(tsv_path: Path) -> dict:
    mapping = {}
    with open(tsv_path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            parts = line.split("\t")
            if len(parts) < 2:
                print(f"WARNING: skipping malformed line {lineno}: {line!r}",
                      file=sys.stderr)
                continue
            old_id, new_id = parts[0].strip(), parts[1].strip()
            if old_id in mapping:
                print(f"WARNING: duplicate old ID '{old_id}' at line {lineno}, overwriting",
                      file=sys.stderr)
            mapping[old_id] = new_id
    return mapping


def _write_sequence(seq: list, out_fh) -> None:
    """Write buffered sequence lines wrapped at FASTA_LINE_WIDTH characters."""
    joined = "".join(seq)
    for i in range(0, len(joined), FASTA_LINE_WIDTH):
        out_fh.write(joined[i:i + FASTA_LINE_WIDTH] + "\n")


def rename_fasta(fasta_path: Path, mapping: dict, out_fh) -> tuple:
    """Write renamed FASTA to out_fh. Returns (renamed, unchanged) counts."""
    renamed   = 0
    unchanged = 0
    seq_buf   = []

    def flush(header_line):
        out_fh.write(header_line)
        _write_sequence(seq_buf, out_fh)
        seq_buf.clear()

    pending_header = None

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pending_header is not None:
                    flush(pending_header)
                header = line[1:].rstrip("\n")
                parts  = header.split(None, 1)
                seq_id = parts[0]
                description = (" " + parts[1]) if len(parts) > 1 else ""
                if seq_id in mapping:
                    pending_header = f">{mapping[seq_id]}{description}\n"
                    renamed += 1
                else:
                    print(f"WARNING: SeqID '{seq_id}' not in mapping, keeping original",
                          file=sys.stderr)
                    pending_header = f">{seq_id}{description}\n"
                    unchanged += 1
            else:
                seq_buf.append(line.rstrip("\n"))

    if pending_header is not None:
        flush(pending_header)

    return renamed, unchanged


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="FastaRename",
        description="Rename sequence IDs in a FASTA file using a TSV mapping.",
    )
    ap.add_argument("--fasta",  required=True, type=Path,
                    help="Input FASTA file")
    ap.add_argument("--tsv",    required=True, type=Path,
                    help="Two-column TSV mapping: old_seqid <TAB> new_seqid")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output FASTA file (default: stdout)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    for path, label in [(args.fasta, "--fasta"), (args.tsv, "--tsv")]:
        if not path.exists():
            print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    mapping = load_mapping(args.tsv)
    print(f"Loaded {len(mapping)} ID mappings from {args.tsv.name}", file=sys.stderr)

    if args.output:
        with open(args.output, "w") as out_fh:
            renamed, unchanged = rename_fasta(args.fasta, mapping, out_fh)
        print(f"Written: {args.output}", file=sys.stderr)
    else:
        renamed, unchanged = rename_fasta(args.fasta, mapping, sys.stdout)

    print(f"Renamed: {renamed}  |  Unchanged: {unchanged}", file=sys.stderr)


if __name__ == "__main__":
    main()
