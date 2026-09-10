#!/usr/bin/env python3
"""
SAM2SubgenomeBED.py — Tag a polyploid assembly by parental subgenome of
origin, from a SAM alignment against a concatenated parental reference.

Generic across any allopolyploid with two distinguishable parental
subgenomes: an assembly (e.g. Nicotiana tabacum, formed by hybridization
of N. tomentosiformis and N. sylvestris) mapped with minimap2 (e.g.
-a -x asm5/asm10) as the *query* against a single joint reference FASTA
containing both parental genomes, with sequence names distinguishable by
a prefix/pattern (e.g. Ntom_chr1 / Nsyl_chr1, or any other naming
scheme — set via --parent1_pattern/--parent2_pattern). The resulting SAM
is parsed here to classify each region of every query (assembly)
sequence by parent of origin, using whatever labels --parent1_name/
--parent2_name are given (default: "Parent1"/"Parent2").

For each query sequence, primary and supplementary alignment records are
used to recover which parent each part of the query best matches
(hard/soft clips are resolved back to original, un-reverse-complemented
query coordinates so intervals from + and - strand records are directly
comparable). Coverage is then tallied in fixed-size windows; each window
is labelled by whichever parent accounts for the larger share of
aligned bases, "ambiguous" if the two are too close to call, or
"unclassified" if too little of the window is covered by any alignment.
Adjacent windows with the same label are merged into single BED
intervals.

This does not resolve truly ambiguous regions (recent duplications,
shared repeats, rDNA) any better than the underlying alignment does —
those come out labelled "ambiguous" rather than a guess.

Usage
-----
    SAM2SubgenomeBED.py --sam tabacum_vs_joint.sam \\
        --parent1_name Ntom --parent1_pattern '^Ntom' \\
        --parent2_name Nsyl --parent2_pattern '^Nsyl' \\
        --output subgenomes.bed
    SAM2SubgenomeBED.py --sam wheat_vs_joint.sam \\
        --parent1_name A_genome --parent1_pattern '^chrA' \\
        --parent2_name B_genome --parent2_pattern '^chrB' \\
        --window_size 20000
    SAM2SubgenomeBED.py --sam tabacum_vs_joint.sam \\
        --parent1_pattern Ntom --parent2_pattern Nsyl --dry_run
"""

import argparse
import re
import sys
from pathlib import Path

VERSION = "v0.2.0"

DEFAULT_WINDOW_SIZE     = 10000
DEFAULT_MIN_MAPQ        = 5
DEFAULT_MIN_ALN_LEN     = 1000
DEFAULT_MIN_COVERAGE    = 0.5
DEFAULT_AMBIGUITY_MARGIN = 0.1

FLAG_UNMAPPED  = 0x4
FLAG_REVERSE   = 0x10
FLAG_SECONDARY = 0x100

CIGAR_OP_RE = re.compile(r"(\d+)([MIDNSHP=X])")


# ── SAM parsing ────────────────────────────────────────────────────────────

def classify_targets(sq_names: list, parent1_name: str, parent1_pattern: str,
                     parent2_name: str, parent2_pattern: str) -> dict:
    """Return {target_name: parent1_name|parent2_name|None} from @SQ SN values."""
    p1_re = re.compile(parent1_pattern)
    p2_re = re.compile(parent2_pattern)
    labels = {}
    for name in sq_names:
        is_p1 = p1_re.search(name) is not None
        is_p2 = p2_re.search(name) is not None
        if is_p1 and not is_p2:
            labels[name] = parent1_name
        elif is_p2 and not is_p1:
            labels[name] = parent2_name
        else:
            labels[name] = None  # unmatched or matches both patterns
    return labels


def parse_query_span(cigar: str, is_reverse: bool) -> tuple:
    """Return (q_start, q_end, q_total_len) in original, forward-strand
    query coordinates, resolving soft/hard clips on either strand."""
    ops = CIGAR_OP_RE.findall(cigar)
    i, j = 0, len(ops) - 1
    lead_clip = 0
    while i <= j and ops[i][1] in "SH":
        lead_clip += int(ops[i][0])
        i += 1
    trail_clip = 0
    while j >= i and ops[j][1] in "SH":
        trail_clip += int(ops[j][0])
        j -= 1
    aligned = sum(int(length) for length, op in ops[i:j + 1] if op in "MIX=")
    total_len = lead_clip + aligned + trail_clip
    if is_reverse:
        q_start = trail_clip
    else:
        q_start = lead_clip
    return q_start, q_start + aligned, total_len


def parse_sam(sam_path: Path, target_labels: dict, min_mapq: int,
              min_aln_len: int, include_secondary: bool) -> tuple:
    """Return (intervals, query_lengths, stats).
    intervals: {qname: [(start, end, label), ...]}
    query_lengths: {qname: total_len}
    """
    intervals = {}
    query_lengths = {}
    stats = {"records": 0, "unmapped": 0, "secondary_skipped": 0,
             "unclassified_target": 0, "filtered_mapq_or_len": 0, "used": 0}

    with open(sam_path) as fh:
        for line in fh:
            if line.startswith("@"):
                continue
            cols = line.rstrip("\n").split("\t")
            if len(cols) < 11:
                continue
            stats["records"] += 1
            qname, flag_s, rname, pos_s, mapq_s, cigar = cols[0:6]
            flag = int(flag_s)

            if flag & FLAG_UNMAPPED or rname == "*" or cigar == "*":
                stats["unmapped"] += 1
                continue
            if flag & FLAG_SECONDARY and not include_secondary:
                stats["secondary_skipped"] += 1
                continue

            label = target_labels.get(rname)
            if label is None:
                stats["unclassified_target"] += 1
                continue

            mapq = int(mapq_s)
            q_start, q_end, q_total = parse_query_span(cigar, bool(flag & FLAG_REVERSE))
            aln_len = q_end - q_start
            if mapq < min_mapq or aln_len < min_aln_len:
                stats["filtered_mapq_or_len"] += 1
                continue

            intervals.setdefault(qname, []).append((q_start, q_end, label))
            query_lengths[qname] = max(query_lengths.get(qname, 0), q_total)
            stats["used"] += 1

    return intervals, query_lengths, stats


# ── Windowing / classification ──────────────────────────────────────────────

def window_labels(query_len: int, aln_intervals: list, window_size: int,
                  min_coverage: float, ambiguity_margin: float,
                  parent1_name: str, parent2_name: str) -> list:
    """Return per-window (start, end, label, confidence) for one query."""
    n_windows = (query_len + window_size - 1) // window_size
    p1_bp = [0] * n_windows
    p2_bp = [0] * n_windows

    for start, end, label in aln_intervals:
        bucket = p1_bp if label == parent1_name else p2_bp
        w0 = start // window_size
        w1 = (end - 1) // window_size if end > start else w0
        for w in range(w0, min(w1, n_windows - 1) + 1):
            win_start = w * window_size
            win_end = min(win_start + window_size, query_len)
            overlap = min(end, win_end) - max(start, win_start)
            if overlap > 0:
                bucket[w] += overlap

    windows = []
    for w in range(n_windows):
        win_start = w * window_size
        win_end = min(win_start + window_size, query_len)
        win_len = win_end - win_start
        covered = min(p1_bp[w] + p2_bp[w], win_len)
        frac_covered = covered / win_len if win_len else 0.0

        if frac_covered < min_coverage or covered == 0:
            label, conf = "unclassified", 0.0
        else:
            total = p1_bp[w] + p2_bp[w]
            diff_frac = abs(p1_bp[w] - p2_bp[w]) / total
            if diff_frac < ambiguity_margin:
                label, conf = "ambiguous", 0.5 + diff_frac / 2
            else:
                label = parent1_name if p1_bp[w] > p2_bp[w] else parent2_name
                conf = max(p1_bp[w], p2_bp[w]) / total
        windows.append((win_start, win_end, label, conf))
    return windows


def merge_windows(windows: list) -> list:
    """Collapse adjacent windows sharing the same label into one interval,
    averaging confidence weighted by window length."""
    merged = []
    for start, end, label, conf in windows:
        if merged and merged[-1][2] == label:
            p_start, p_end, p_label, p_conf, p_len = merged[-1]
            new_len = p_len + (end - start)
            new_conf = (p_conf * p_len + conf * (end - start)) / new_len
            merged[-1] = (p_start, end, label, new_conf, new_len)
        else:
            merged.append((start, end, label, conf, end - start))
    return [(s, e, l, c) for s, e, l, c, _ in merged]


# ── Output ───────────────────────────────────────────────────────────────

def write_bed(fh, qname: str, merged: list) -> None:
    for start, end, label, conf in merged:
        score = round(conf * 1000)
        fh.write(f"{qname}\t{start}\t{end}\t{label}\t{score}\n")


# ── Main ─────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="SAM2SubgenomeBED",
        description="Classify regions of a polyploid assembly as "
                     "derived from one parental subgenome or the other, "
                     "from a SAM alignment against a concatenated "
                     "parental reference. Generic across any two-parent "
                     "allopolyploid via --parent1_*/--parent2_*.",
    )
    ap.add_argument("--sam", required=True, type=Path,
                    help="Input SAM file (assembly-vs-joint-reference "
                         "alignment, e.g. minimap2 -a -x asm5). BAM is not "
                         "supported (stdlib-only script) — convert first "
                         "with 'samtools view -h'.")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output BED file (default: stdout)")
    ap.add_argument("--parent1_name", default="Parent1",
                    help="Label written to the output BED for the first "
                         "parental subgenome (default: 'Parent1')")
    ap.add_argument("--parent1_pattern", required=True,
                    help="Regex matched against @SQ target names to "
                         "identify --parent1_name's subgenome (required)")
    ap.add_argument("--parent2_name", default="Parent2",
                    help="Label written to the output BED for the second "
                         "parental subgenome (default: 'Parent2')")
    ap.add_argument("--parent2_pattern", required=True,
                    help="Regex matched against @SQ target names to "
                         "identify --parent2_name's subgenome (required)")
    ap.add_argument("--window_size", type=int, default=DEFAULT_WINDOW_SIZE,
                    help=f"Window size in bp for majority-vote subgenome "
                         f"calls (default: {DEFAULT_WINDOW_SIZE})")
    ap.add_argument("--min_mapq", type=int, default=DEFAULT_MIN_MAPQ,
                    help=f"Minimum MAPQ for an alignment record to be "
                         f"used (default: {DEFAULT_MIN_MAPQ})")
    ap.add_argument("--min_aln_len", type=int, default=DEFAULT_MIN_ALN_LEN,
                    help=f"Minimum aligned query length in bp for an "
                         f"alignment record to be used "
                         f"(default: {DEFAULT_MIN_ALN_LEN})")
    ap.add_argument("--min_coverage", type=float, default=DEFAULT_MIN_COVERAGE,
                    help=f"Minimum fraction of a window that must be "
                         f"covered by used alignments to make a call, "
                         f"otherwise 'unclassified' "
                         f"(default: {DEFAULT_MIN_COVERAGE})")
    ap.add_argument("--ambiguity_margin", type=float,
                    default=DEFAULT_AMBIGUITY_MARGIN,
                    help=f"Minimum relative difference between "
                         f"--parent1_name and --parent2_name covered bp "
                         f"within a window to call a side; below this "
                         f"the window is 'ambiguous' "
                         f"(default: {DEFAULT_AMBIGUITY_MARGIN})")
    ap.add_argument("--include_secondary", action="store_true",
                    help="Also use secondary alignment records (FLAG "
                         "0x100), not just primary and supplementary "
                         "(default: primary + supplementary only)")
    ap.add_argument("--dry_run", action="store_true",
                    help="Parse the SAM header and report target "
                         "classification, then exit without processing "
                         "alignment records or writing output")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if args.parent1_name == args.parent2_name:
        print(f"ERROR: --parent1_name and --parent2_name must differ "
              f"(both were '{args.parent1_name}')", file=sys.stderr)
        sys.exit(1)

    if not args.sam.exists():
        print(f"ERROR: --sam file not found: {args.sam}", file=sys.stderr)
        sys.exit(1)

    sq_names = []
    with open(args.sam) as fh:
        for line in fh:
            if not line.startswith("@"):
                break
            if line.startswith("@SQ"):
                for field in line.rstrip("\n").split("\t"):
                    if field.startswith("SN:"):
                        sq_names.append(field[3:])

    if not sq_names:
        print(f"ERROR: no @SQ header lines found in {args.sam} — is this "
              f"a valid SAM file with a header (minimap2 -a)?",
              file=sys.stderr)
        sys.exit(1)

    target_labels = classify_targets(
        sq_names, args.parent1_name, args.parent1_pattern,
        args.parent2_name, args.parent2_pattern,
    )
    n_p1 = sum(1 for v in target_labels.values() if v == args.parent1_name)
    n_p2 = sum(1 for v in target_labels.values() if v == args.parent2_name)
    n_unclassified = sum(1 for v in target_labels.values() if v is None)

    print(f"Reference targets: {len(sq_names)} total, "
          f"{n_p1} {args.parent1_name}, {n_p2} {args.parent2_name}, "
          f"{n_unclassified} unclassified", file=sys.stderr)
    if n_unclassified:
        examples = [n for n, v in target_labels.items() if v is None][:3]
        print(f"  WARNING: {n_unclassified} target(s) matched neither/both "
              f"patterns and will be ignored, e.g.: {', '.join(examples)}",
              file=sys.stderr)

    if n_p1 == 0 or n_p2 == 0:
        print(f"\nERROR: --parent1_pattern ('{args.parent1_pattern}') "
              f"matched {n_p1} target(s) and --parent2_pattern "
              f"('{args.parent2_pattern}') matched {n_p2} target(s) — "
              f"both must be > 0.\n\n"
              f"  Example target names (first 5): "
              f"{', '.join(sq_names[:5])}\n\n"
              f"Likely cause: the patterns don't match this SAM's @SQ "
              f"naming scheme. Fix with "
              f"--parent1_pattern/--parent2_pattern.",
              file=sys.stderr)
        sys.exit(1)

    if args.dry_run:
        print(f"  Parent 1          : {args.parent1_name} "
              f"(pattern: {args.parent1_pattern})", file=sys.stderr)
        print(f"  Parent 2          : {args.parent2_name} "
              f"(pattern: {args.parent2_pattern})", file=sys.stderr)
        print(f"  Window size       : {args.window_size} bp", file=sys.stderr)
        print(f"  Min MAPQ          : {args.min_mapq}", file=sys.stderr)
        print(f"  Min aligned length: {args.min_aln_len} bp", file=sys.stderr)
        print(f"  Min window coverage: {args.min_coverage}", file=sys.stderr)
        print(f"  Ambiguity margin  : {args.ambiguity_margin}", file=sys.stderr)
        print(f"  Include secondary : {args.include_secondary}", file=sys.stderr)
        print(f"  Output            : {args.output or 'stdout'}", file=sys.stderr)
        print("  Exiting (--dry_run). No alignment records processed.",
              file=sys.stderr)
        return

    intervals, query_lengths, stats = parse_sam(
        args.sam, target_labels, args.min_mapq, args.min_aln_len,
        args.include_secondary,
    )
    print(f"Alignment records: {stats['records']} total, {stats['used']} used, "
          f"{stats['unmapped']} unmapped, "
          f"{stats['secondary_skipped']} secondary skipped, "
          f"{stats['unclassified_target']} on unclassified target, "
          f"{stats['filtered_mapq_or_len']} filtered (mapq/length)",
          file=sys.stderr)

    if not intervals:
        print("ERROR: no usable alignment records after filtering — "
              "nothing to classify. Check --min_mapq/--min_aln_len or the "
              "target patterns.", file=sys.stderr)
        sys.exit(1)

    out_fh = open(args.output, "w") if args.output else sys.stdout
    totals = {args.parent1_name: 0, args.parent2_name: 0,
              "ambiguous": 0, "unclassified": 0}
    try:
        for qname in sorted(intervals):
            windows = window_labels(
                query_lengths[qname], intervals[qname], args.window_size,
                args.min_coverage, args.ambiguity_margin,
                args.parent1_name, args.parent2_name,
            )
            merged = merge_windows(windows)
            write_bed(out_fh, qname, merged)
            for start, end, label, _ in merged:
                totals[label] += end - start
    finally:
        if args.output:
            out_fh.close()

    total_bp = sum(totals.values())
    print(f"\nSubgenome classification summary ({total_bp:,} bp total, "
          f"{len(intervals)} query sequence(s)):", file=sys.stderr)
    for label in (args.parent1_name, args.parent2_name, "ambiguous", "unclassified"):
        pct = 100 * totals[label] / total_bp if total_bp else 0.0
        print(f"  {label:<13}: {totals[label]:>14,} bp ({pct:5.1f}%)",
              file=sys.stderr)


if __name__ == "__main__":
    main()
