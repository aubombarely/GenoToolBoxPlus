#!/usr/bin/env python3
"""
FastaRename.py — Rename sequence IDs in a FASTA file using a TSV mapping.

The TSV file must have two tab-separated columns (no header required):
    old_seqid <TAB> new_seqid

Lines starting with '#' are treated as comments and ignored.
SeqIDs not found in the mapping are kept unchanged (with a warning).
Only the ID portion of the header (up to the first space) is renamed;
any description text is preserved. Output sequences are wrapped at 60
characters per line.

Usage
-----
    FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv
    FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv --output renamed.fasta
    FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv --sort length_descending
"""

import argparse
import sys
from pathlib import Path

VERSION = "v0.0.1"

FASTA_LINE_WIDTH = 60

SORT_OPTIONS = [
    "length_ascending",
    "length_descending",
    "alphab_ascending",
    "alphab_descending",
]


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


def load_fasta(fasta_path: Path, mapping: dict) -> tuple:
    """
    Load all sequences from fasta_path, applying ID renaming.
    Returns (records, renamed, unchanged) where records is a list of
    (new_seq_id, description, sequence_str) tuples.
    """
    records   = []
    renamed   = 0
    unchanged = 0
    seq_buf   = []
    pending   = None  # (new_seq_id, description)

    def flush():
        seq_id, description = pending
        records.append((seq_id, description, "".join(seq_buf)))
        seq_buf.clear()

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if pending is not None:
                    flush()
                header = line[1:].rstrip("\n")
                parts  = header.split(None, 1)
                seq_id = parts[0]
                description = (" " + parts[1]) if len(parts) > 1 else ""
                if seq_id in mapping:
                    pending = (mapping[seq_id], description)
                    renamed += 1
                else:
                    print(f"WARNING: SeqID '{seq_id}' not in mapping, keeping original",
                          file=sys.stderr)
                    pending = (seq_id, description)
                    unchanged += 1
            else:
                seq_buf.append(line.rstrip("\n"))

    if pending is not None:
        flush()

    return records, renamed, unchanged


def sort_records(records: list, sort_order: str) -> list:
    if sort_order == "length_ascending":
        return sorted(records, key=lambda r: len(r[2]))
    elif sort_order == "length_descending":
        return sorted(records, key=lambda r: len(r[2]), reverse=True)
    elif sort_order == "alphab_ascending":
        return sorted(records, key=lambda r: r[0])
    elif sort_order == "alphab_descending":
        return sorted(records, key=lambda r: r[0], reverse=True)
    return records


def write_records(records: list, out_fh) -> None:
    for seq_id, description, seq in records:
        out_fh.write(f">{seq_id}{description}\n")
        for i in range(0, len(seq), FASTA_LINE_WIDTH):
            out_fh.write(seq[i:i + FASTA_LINE_WIDTH] + "\n")


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
    ap.add_argument("--sort",   choices=SORT_OPTIONS, default=None,
                    help="Sort output sequences (default: input order)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    for path, label in [(args.fasta, "--fasta"), (args.tsv, "--tsv")]:
        if not path.exists():
            print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    mapping = load_mapping(args.tsv)
    print(f"Loaded {len(mapping)} ID mappings from {args.tsv.name}", file=sys.stderr)

    records, renamed, unchanged = load_fasta(args.fasta, mapping)

    if args.sort:
        records = sort_records(records, args.sort)
        print(f"Sorted: {args.sort}", file=sys.stderr)

    if args.output:
        with open(args.output, "w") as out_fh:
            write_records(records, out_fh)
        print(f"Written: {args.output}", file=sys.stderr)
    else:
        write_records(records, sys.stdout)

    print(f"Renamed: {renamed}  |  Unchanged: {unchanged}", file=sys.stderr)


if __name__ == "__main__":
    main()
