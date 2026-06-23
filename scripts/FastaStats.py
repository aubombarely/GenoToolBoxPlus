#!/usr/bin/env python3
"""
FastaStats.py — Compute common statistics for sequences in a FASTA file.

Assembly statistics:
    num_sequences, total_length, min_length, max_length, mean_length,
    median_length, n50, l50, n90, l90, gc_content, n_count, n_percent,
    seq_gt_1kb, seq_gt_10kb, seq_gt_100kb, seq_gt_1mb

Nucleotide composition:
    Count and percentage for every IUPAC base present in the sequences
    (A C G T U R Y M K S W H B V D N). Non-IUPAC characters are
    reported as 'other'.

Usage
-----
    FastaStats.py --fasta genome.fasta
    FastaStats.py --fasta genome.fasta --output genome_stats
    FastaStats.py --fasta genome.fasta --output genome_stats --format tsv,txt,json
"""

import argparse
import json
import sys
from datetime import datetime
from pathlib import Path
from statistics import median

VERSION = "v0.1.0"

FORMAT_EXT = {"tsv": ".tsv", "txt": ".txt", "json": ".json"}

ASSEMBLY_METRICS = [
    "num_sequences",
    "total_length",
    "min_length",
    "max_length",
    "mean_length",
    "median_length",
    "n50",
    "l50",
    "n90",
    "l90",
    "gc_content",
    "n_count",
    "n_percent",
    "seq_gt_1kb",
    "seq_gt_10kb",
    "seq_gt_100kb",
    "seq_gt_1mb",
]

IUPAC_BASES = list("ACGTURYMKSWHBVDN")


# ── Input ──────────────────────────────────────────────────────────────────────

def parse_formats(fmt_str: str) -> list:
    formats = []
    for f in fmt_str.split(","):
        f = f.strip().lower()
        if f not in FORMAT_EXT:
            print(f"ERROR: unknown format '{f}'. Valid: {', '.join(FORMAT_EXT)}",
                  file=sys.stderr)
            sys.exit(1)
        if f not in formats:
            formats.append(f)
    return formats


def load_fasta(fasta_path: Path) -> list:
    """Load sequences. Returns list of (seq_id, sequence_str) tuples."""
    sequences = []
    seq_id    = None
    seq_buf   = []

    with open(fasta_path) as fh:
        for line in fh:
            line = line.rstrip("\n")
            if line.startswith(">"):
                if seq_id is not None:
                    sequences.append((seq_id, "".join(seq_buf)))
                seq_id  = line[1:].split()[0]
                seq_buf = []
            else:
                seq_buf.append(line)
    if seq_id is not None:
        sequences.append((seq_id, "".join(seq_buf)))

    return sequences


# ── Statistics ─────────────────────────────────────────────────────────────────

def _nx_lx(lengths_desc: list, total: int, fraction: float) -> tuple:
    """Return (Nx, Lx) for a given fraction (0–1). lengths_desc sorted descending."""
    target = total * fraction
    cumsum = 0
    for i, length in enumerate(lengths_desc, 1):
        cumsum += length
        if cumsum >= target:
            return length, i
    return lengths_desc[-1], len(lengths_desc)


def compute_stats(sequences: list) -> tuple:
    """
    Return (assembly_stats dict, nt_composition dict).
    nt_composition maps base -> {"count": int, "percent": float}.
    """
    if not sequences:
        return ({m: 0 for m in ASSEMBLY_METRICS}, {})

    lengths      = [len(seq) for _, seq in sequences]
    lengths_desc = sorted(lengths, reverse=True)
    total        = sum(lengths)

    n50, l50 = _nx_lx(lengths_desc, total, 0.50)
    n90, l90 = _nx_lx(lengths_desc, total, 0.90)

    all_seq  = "".join(seq.upper() for _, seq in sequences)
    gc_count = all_seq.count("G") + all_seq.count("C")
    n_count  = all_seq.count("N")

    assembly = {
        "num_sequences": len(sequences),
        "total_length":  total,
        "min_length":    min(lengths),
        "max_length":    max(lengths),
        "mean_length":   round(total / len(lengths), 2),
        "median_length": median(lengths),
        "n50":           n50,
        "l50":           l50,
        "n90":           n90,
        "l90":           l90,
        "gc_content":    round(gc_count / total * 100, 4) if total else 0.0,
        "n_count":       n_count,
        "n_percent":     round(n_count  / total * 100, 4) if total else 0.0,
        "seq_gt_1kb":    sum(1 for l in lengths if l >      1_000),
        "seq_gt_10kb":   sum(1 for l in lengths if l >     10_000),
        "seq_gt_100kb":  sum(1 for l in lengths if l >    100_000),
        "seq_gt_1mb":    sum(1 for l in lengths if l >  1_000_000),
    }

    iupac_set = set(IUPAC_BASES)
    nt_comp   = {}
    for base in IUPAC_BASES:
        count = all_seq.count(base)
        if count:
            nt_comp[base] = {"count": count,
                             "percent": round(count / total * 100, 4)}
    other = sum(1 for c in all_seq if c not in iupac_set)
    if other:
        nt_comp["other"] = {"count": other,
                            "percent": round(other / total * 100, 4)}

    return assembly, nt_comp


# ── Writers ───────────────────────────────────────────────────────────────────

def write_tsv(assembly: dict, nt_comp: dict, path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("# Assembly statistics\n")
        fh.write("metric\tvalue\n")
        for m in ASSEMBLY_METRICS:
            fh.write(f"{m}\t{assembly[m]}\n")
        fh.write("\n# Nucleotide composition\n")
        fh.write("base\tcount\tpercent\n")
        for base, vals in nt_comp.items():
            fh.write(f"{base}\t{vals['count']}\t{vals['percent']}\n")


def write_txt(assembly: dict, nt_comp: dict, path: Path) -> None:
    col_w = max(len(m) for m in ASSEMBLY_METRICS)
    with open(path, "w") as fh:
        fh.write("=== Assembly Statistics ===\n")
        for m in ASSEMBLY_METRICS:
            fh.write(f"{m:<{col_w}}  {assembly[m]}\n")
        fh.write("\n=== Nucleotide Composition ===\n")
        fh.write(f"{'Base':<6}  {'Count':>12}  {'Percent (%)':>12}\n")
        for base, vals in nt_comp.items():
            fh.write(f"{base:<6}  {vals['count']:>12}  {vals['percent']:>12}\n")


def write_json(assembly: dict, nt_comp: dict, path: Path) -> None:
    data = {"assembly_stats": assembly, "nucleotide_composition": nt_comp}
    with open(path, "w") as fh:
        json.dump(data, fh, indent=2)
        fh.write("\n")


WRITERS = {"tsv": write_tsv, "txt": write_txt, "json": write_json}


# ── Stderr reporting ───────────────────────────────────────────────────────────

def _ascii_table(headers: list, rows: list) -> list:
    """Return list of strings forming an ASCII box table."""
    widths = [max(len(h), *(len(r[i]) for r in rows))
              for i, h in enumerate(headers)]
    sep = "+" + "+".join("-" * (w + 2) for w in widths) + "+"

    def fmt_row(cells):
        parts = []
        for i, cell in enumerate(cells):
            padded = cell.ljust(widths[i]) if i == 0 else cell.rjust(widths[i])
            parts.append(f" {padded} ")
        return "|" + "|".join(parts) + "|"

    lines = [sep, fmt_row(headers), sep]
    for r in rows:
        lines.append(fmt_row(r))
    lines.append(sep)
    return lines


def print_stderr(assembly: dict, nt_comp: dict, command: str) -> None:
    timestamp = datetime.now().strftime("%Y-%m-%d %H:%M:%S")
    print(f"FastaStats.py  {timestamp}", file=sys.stderr)
    print(f"Command: {command}", file=sys.stderr)
    print("", file=sys.stderr)

    # Assembly statistics table
    print("Assembly Statistics:", file=sys.stderr)
    rows = [[m, str(assembly[m])] for m in ASSEMBLY_METRICS]
    for line in _ascii_table(["Metric", "Value"], rows):
        print(line, file=sys.stderr)
    print("", file=sys.stderr)

    # Nucleotide composition table
    print("Nucleotide Composition:", file=sys.stderr)
    rows = [[b, str(v["count"]), str(v["percent"])]
            for b, v in nt_comp.items()]
    for line in _ascii_table(["Base", "Count", "Percent (%)"], rows):
        print(line, file=sys.stderr)


# ── Main ──────────────────────────────────────────────────────────────────────

def default_output(fasta_path: Path) -> str:
    stem = fasta_path.name
    for ext in (".fasta", ".fa", ".fna", ".ffn", ".faa", ".frn"):
        if stem.endswith(ext):
            return stem[: -len(ext)]
    return stem


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="FastaStats",
        description="Compute common statistics for sequences in a FASTA file.",
    )
    ap.add_argument("--fasta",  required=True, type=Path,
                    help="Input FASTA file")
    ap.add_argument("--output", default=None,
                    help="Output basename (default: derived from --fasta filename)")
    ap.add_argument("--format", default="tsv",
                    help="Comma-separated output formats: tsv, txt, json (default: tsv)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: --fasta file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    formats  = parse_formats(args.format)
    basename = args.output or default_output(args.fasta)
    command  = " ".join(sys.argv)

    print(f"Loading sequences from {args.fasta.name} ...", file=sys.stderr)
    sequences = load_fasta(args.fasta)
    print(f"  {len(sequences)} sequences loaded.", file=sys.stderr)

    assembly, nt_comp = compute_stats(sequences)

    for fmt in formats:
        out_path = Path(basename + FORMAT_EXT[fmt])
        WRITERS[fmt](assembly, nt_comp, out_path)
        print(f"Written: {out_path}", file=sys.stderr)

    print("", file=sys.stderr)
    print_stderr(assembly, nt_comp, command)


if __name__ == "__main__":
    main()
