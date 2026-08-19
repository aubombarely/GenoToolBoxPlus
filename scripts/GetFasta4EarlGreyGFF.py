#!/usr/bin/env python3
"""
GetFasta4EarlGreyGFF.py — Extract FASTA sequences for TE features from an
EarlGrey repeat-annotation GFF3.

Reads a genome FASTA and an EarlGrey GFF3 (column 3 = TE type, e.g.
LINE/L1, LTR/Copia; attributes include ID=<family ID>, e.g.
RND-1_FAMILY-789) and writes one FASTA record per GFF3 feature. Strand
'-' features are reverse-complemented. Extraction streams the genome
FASTA sequence-by-sequence (peak memory = one chromosome), so the GFF3
is grouped by seqid up front.

Output header format:
    >{ID}_{TYPE}_{SeqID}_{Start}

where TYPE has every '/' replaced with '_' (e.g. LTR/Copia -> LTR_Copia)
so the header/filename stays shell- and tool-safe.

Example
-------
GFF3 line:
    PhangAGP1C01  Earl_Grey  LINE/L1  1  2101  10800  +  .  TSTART=5686;TEND=7874;ID=RND-1_FAMILY-789;SHORTTE=F;KIMURA80=0.2841

produces:
    >RND-1_FAMILY-789_LINE_L1_PhangAGP1C01_1

Usage
-----
    GetFasta4EarlGreyGFF.py --fasta genome.fasta --gff repeats.gff3
    GetFasta4EarlGreyGFF.py --fasta genome.fasta --gff repeats.gff3 --output TEs.fasta
"""

import argparse
import re
import sys
from pathlib import Path

VERSION = "v0.0.1"

FASTA_LINE_WIDTH = 60
_COMPLEMENT = str.maketrans("ACGTNacgtn", "TGCANtgcan")


def revcomp(seq: str) -> str:
    return seq.translate(_COMPLEMENT)[::-1]


def sanitize_type(te_type: str) -> str:
    """Replace '/' with '_' so the TE type is safe in a FASTA header."""
    return te_type.replace("/", "_")


def parse_gff(gff_path: Path) -> dict:
    """
    Parse an EarlGrey GFF3, grouping features by seqid.
    Returns {seqid: [{"id", "type", "start", "end", "strand"}, ...]}.
    """
    by_seq: dict = {}
    n_no_id = 0
    with open(gff_path) as fh:
        for lineno, line in enumerate(fh, 1):
            if line.startswith("#") or not line.strip():
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 9:
                print(f"WARNING: skipping malformed GFF3 line {lineno} "
                      f"(fewer than 9 columns)", file=sys.stderr)
                continue
            seqid, _source, te_type, start, end, _score, strand = cols[:7]
            attrs = cols[8]
            try:
                start_i, end_i = int(start), int(end)
            except ValueError:
                print(f"WARNING: skipping GFF3 line {lineno} "
                      f"(non-integer start/end)", file=sys.stderr)
                continue

            id_m = re.search(r"ID=([^;]+)", attrs)
            if id_m:
                feat_id = id_m.group(1)
            else:
                n_no_id += 1
                feat_id = f"{seqid}_{start_i}_{end_i}"

            by_seq.setdefault(seqid, []).append({
                "id": feat_id, "type": te_type,
                "start": start_i, "end": end_i, "strand": strand,
            })

    if n_no_id:
        print(f"WARNING: {n_no_id} feature(s) had no ID= attribute; "
              f"used seqid_start_end as a fallback ID", file=sys.stderr)
    return by_seq


def extract_sequences(fasta_path: Path, by_seq: dict) -> tuple:
    """
    Stream the genome FASTA and extract feature sequences.
    Peak memory = one chromosome. Returns
    (records, n_extracted, n_skipped, fasta_seqids).
    records is a list of (header, sequence) tuples; fasta_seqids is the
    set of every seqid seen in the FASTA (used to report GFF3 seqids
    that have no match in the genome).
    """
    records: list = []
    n_extracted = 0
    n_skipped   = 0
    fasta_seqids: set = set()
    current_name: str | None  = None
    current_parts: list | None = None

    def flush(name, parts):
        nonlocal n_extracted, n_skipped
        if name not in by_seq or parts is None:
            return
        seq     = "".join(parts)
        seq_len = len(seq)
        for r in by_seq[name]:
            start, end = r["start"], r["end"]
            if start < 1 or end > seq_len or start > end:
                print(f"WARNING: feature {r['id']} ({name}:{start}-{end}) "
                      f"out of bounds (sequence length {seq_len}), skipping",
                      file=sys.stderr)
                n_skipped += 1
                continue
            sub = seq[start - 1:end]
            if r["strand"] == "-":
                sub = revcomp(sub)
            header = f"{r['id']}_{sanitize_type(r['type'])}_{name}_{start}"
            records.append((header, sub))
            n_extracted += 1

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                flush(current_name, current_parts)
                current_name  = line[1:].split()[0]
                fasta_seqids.add(current_name)
                current_parts = [] if current_name in by_seq else None
            elif current_parts is not None:
                current_parts.append(line)
    flush(current_name, current_parts)

    return records, n_extracted, n_skipped, fasta_seqids


def write_records(records: list, out_fh) -> None:
    for header, seq in records:
        out_fh.write(f">{header}\n")
        for i in range(0, len(seq), FASTA_LINE_WIDTH):
            out_fh.write(seq[i:i + FASTA_LINE_WIDTH] + "\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GetFasta4EarlGreyGFF",
        description="Extract FASTA sequences for TE features from an "
                     "EarlGrey repeat-annotation GFF3.",
    )
    ap.add_argument("--fasta",  required=True, type=Path,
                    help="Input genome FASTA")
    ap.add_argument("--gff",    required=True, type=Path,
                    help="EarlGrey GFF3 (column 3 = TE type; "
                         "ID= attribute = family ID)")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output FASTA file (default: stdout)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    for path, label in [(args.fasta, "--fasta"), (args.gff, "--gff")]:
        if not path.exists():
            print(f"ERROR: {label} file not found: {path}", file=sys.stderr)
            sys.exit(1)

    by_seq   = parse_gff(args.gff)
    n_feats  = sum(len(v) for v in by_seq.values())
    print(f"Loaded {n_feats} feature(s) across {len(by_seq)} sequence(s) "
          f"from {args.gff.name}", file=sys.stderr)

    records, n_extracted, n_skipped, fasta_seqids = extract_sequences(
        args.fasta, by_seq)

    n_missing_seqs = len(set(by_seq) - fasta_seqids)
    if n_missing_seqs:
        print(f"WARNING: {n_missing_seqs} GFF3 seqid(s) not found in "
              f"{args.fasta.name}; their features were skipped",
              file=sys.stderr)

    if args.output:
        with open(args.output, "w") as out_fh:
            write_records(records, out_fh)
        print(f"Written: {args.output}", file=sys.stderr)
    else:
        write_records(records, sys.stdout)

    print(f"Extracted: {n_extracted}  |  Skipped: {n_skipped}", file=sys.stderr)


if __name__ == "__main__":
    main()
