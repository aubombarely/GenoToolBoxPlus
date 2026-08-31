#!/usr/bin/env python3
"""
GWH4FastaRename.py — Rename sequence IDs in a FASTA file whose headers are
already in GWH (Genome Warehouse) format, using the SeqID/description
embedded in those headers. Does not download anything — operates on a
FASTA file you already have locally.

GWH headers look like (fields separated by whitespace, order/presence of
the optional category text and key=value fields varies by submission):

    >GWHBJVE00000001    OriSeqID=Chr13    Len=38709045
    >GWHDUBQ00000001    Chromosome 1    Complete=T    Circular=F    OriSeqID=Chr01    Len=45819973
    >GWHDUBQ00000014    OriSeqID=contig0001    Len=45917

OriSeqID is required on every record; Len and any other key=value fields
(Complete=, Circular=, ...) are optional and, if present, are carried
through to the equivalence table but not used for renaming. Any free-text
tokens (e.g. "Chromosome 1", "Mitochondrion") are used for classification
alongside OriSeqID.

Sequence ID renaming scheme (same categories as NCBI_DownloadGenome.py):
    CHR (chromosome) -> {prefix}C{NN}
        Number taken from "Chromosome N" text if present, else from
        OriSeqID when it is a strict Chr<digits> label (e.g. Chr01 -> 1).
        OriSeqID labels like ChrUN10 (unplaced) are NOT treated as CHR —
        "UN" makes the remainder non-numeric, so they fall through to
        SCF/CTG like any other non-chromosome sequence.
    MIT (mitochondrion) -> {prefix}MIT{NN}
    PLT (chloroplast/plastid) -> {prefix}PLT{NN}
    SCF (not chromosome/organelle, sequence contains any N) -> {prefix}SCF{NN}
    CTG (not chromosome/organelle, no N in the sequence) -> {prefix}CTG{NN}
Non-chromosome categories are numbered by order of appearance in the file,
not by any number in OriSeqID. Zero-padding width is the width of the
largest number in that category (minimum 2 digits).

Output FASTA headers are '>{new_id}' with the description stripped. An
equivalence table (gwh_id, oriseqid, category, new_seqid) is always
written alongside, or to stdout's sibling path when --equiv is omitted
and --output is given.

Usage
-----
    GWH4FastaRename.py --fasta genome.gwh.fasta --output renamed.fasta
    GWH4FastaRename.py --fasta genome.gwh.fasta --output renamed.fasta --prefix Sp
    GWH4FastaRename.py --fasta genome.gwh.fasta --output renamed.fasta --equiv equiv.tsv
    GWH4FastaRename.py --fasta genome.gwh.fasta --dry_run
"""

import argparse
import re
import sys
from pathlib import Path

VERSION = "v0.1.0"

FASTA_LINE_WIDTH = 60

ROMAN_RE = re.compile(r"^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$")
ROMAN_VALUES = [
    ("M", 1000), ("CM", 900), ("D", 500), ("CD", 400),
    ("C", 100), ("XC", 90), ("L", 50), ("XL", 40),
    ("X", 10), ("IX", 9), ("V", 5), ("IV", 4), ("I", 1),
]
SEX_CHROM_LETTERS = ("X", "Y", "W", "Z")

CHROM_TEXT_RE = re.compile(
    r"(?:chromosome|pseudomolecule|linkage[\s_]?group)[:\s]+([A-Za-z0-9]+)",
    re.IGNORECASE,
)
CHROM_ORISEQID_RE = re.compile(r"^chr0*([0-9]+[A-Za-z]?)$", re.IGNORECASE)
MITO_RE    = re.compile(r"mitochondri", re.IGNORECASE)
PLASTID_RE = re.compile(r"chloroplast|plastid", re.IGNORECASE)


def _roman_to_int(token: str) -> int:
    value = 0
    i = 0
    while i < len(token):
        pair = token[i:i + 2]
        if pair in dict(ROMAN_VALUES):
            value += dict(ROMAN_VALUES)[pair]
            i += 2
        else:
            value += dict(ROMAN_VALUES)[token[i]]
            i += 1
    return value


# ── GWH header parsing ──────────────────────────────────────────────────────

def parse_gwh_header(header: str) -> tuple:
    """Returns (gwh_id, kv_fields, category_text). kv_fields is a dict of
    key=value tokens (e.g. OriSeqID, Len, Complete, Circular); category_text
    is the remaining free-text tokens joined with a single space (e.g.
    "Chromosome 1")."""
    tokens = header.split()
    gwh_id = tokens[0] if tokens else ""
    kv, text_tokens = {}, []
    for tok in tokens[1:]:
        if "=" in tok:
            k, v = tok.split("=", 1)
            kv[k] = v
        else:
            text_tokens.append(tok)
    return gwh_id, kv, " ".join(text_tokens)


def _scan_gwh_fasta(fasta_path: Path) -> list:
    """Single streaming pass. Returns a list of dicts (gwh_id, oriseqid,
    kv, category_text, length, has_n) in file order. Never holds a full
    sequence in memory. Raises ValueError on a missing OriSeqID or a
    duplicate GWH ID."""
    records = []
    gwh_id, kv, category_text, length, has_n = None, {}, "", 0, False

    def flush():
        oriseqid = kv.get("OriSeqID")
        if not oriseqid:
            raise ValueError(f"record '{gwh_id}' has no OriSeqID= field in its header")
        records.append({"gwh_id": gwh_id, "oriseqid": oriseqid, "kv": kv,
                        "category_text": category_text, "length": length, "has_n": has_n})

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if gwh_id is not None:
                    flush()
                gwh_id, kv, category_text = parse_gwh_header(line[1:].rstrip("\n"))
                length, has_n = 0, False
            else:
                seq = line.rstrip("\n")
                length += len(seq)
                if not has_n and ("N" in seq or "n" in seq):
                    has_n = True
    if gwh_id is not None:
        flush()

    seen = {}
    duplicates = []
    for r in records:
        if r["gwh_id"] in seen:
            duplicates.append(r["gwh_id"])
        seen[r["gwh_id"]] = True
    if duplicates:
        examples = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(
            f"{fasta_path.name} has {len(set(duplicates))} duplicate GWH ID(s), "
            f"e.g. {examples} — every sequence in a FASTA must have a unique ID"
        )

    return records


# ── Classification / renaming ───────────────────────────────────────────────

def _classify_record(r: dict) -> tuple:
    """Returns (category, raw_chrom_token_or_None). category is one of
    CHR/MIT/PLT/SCF/CTG. Chromosome number is taken from free-text
    ("Chromosome N") first, then from a strict Chr<digits> OriSeqID label;
    OriSeqID labels like ChrUN10 don't match the strict pattern and fall
    through to SCF/CTG via has_n, same as any other non-chromosome sequence."""
    m = CHROM_TEXT_RE.search(r["category_text"])
    if m:
        return "CHR", m.group(1)
    m = CHROM_ORISEQID_RE.match(r["oriseqid"])
    if m:
        return "CHR", m.group(1)
    if MITO_RE.search(r["category_text"]) or MITO_RE.search(r["oriseqid"]):
        return "MIT", None
    if PLASTID_RE.search(r["category_text"]) or PLASTID_RE.search(r["oriseqid"]):
        return "PLT", None
    return ("SCF" if r["has_n"] else "CTG"), None


def _resolve_chromosome_numbers(classified: dict) -> dict:
    """For category CHR, resolve each raw token to (numeric_value_or_None,
    literal), same logic as NCBI_DownloadGenome.py: genome-wide Roman
    numeral usage is detected before resolving ambiguous single letters
    (X/Y/W/Z), and letter-suffixed polyploid labels (e.g. "1A") are kept
    literal since they aren't a plain number."""
    chrom_tokens = [tok for cat, tok in classified.values() if cat == "CHR"]
    genome_is_roman = any(
        len(t) > 1 and ROMAN_RE.match(t.upper()) and not t.isdigit()
        for t in chrom_tokens
    )

    resolved = {}
    for gwh_id, (cat, tok) in classified.items():
        if cat != "CHR":
            continue
        if tok.isdigit():
            resolved[gwh_id] = (int(tok), tok)
        elif len(tok) == 1 and tok.upper() in SEX_CHROM_LETTERS and not genome_is_roman:
            resolved[gwh_id] = (None, tok.upper())
        elif ROMAN_RE.match(tok.upper()) and (genome_is_roman or len(tok) > 1):
            resolved[gwh_id] = (_roman_to_int(tok.upper()), tok.upper())
        else:
            resolved[gwh_id] = (None, tok)
    return resolved


def build_rename_mapping(records: list, prefix: str) -> dict:
    """Returns {gwh_id: new_seq_id}, raising ValueError on a naming collision."""
    classified = {r["gwh_id"]: _classify_record(r) for r in records}
    chrom_numbers = _resolve_chromosome_numbers(classified)

    max_chrom_val = max([v for v, _ in chrom_numbers.values() if v is not None], default=0)
    chrom_width = max(2, len(str(max_chrom_val)))

    counts = {"SCF": 0, "CTG": 0, "MIT": 0, "PLT": 0}
    for _, (cat, _) in classified.items():
        if cat in counts:
            counts[cat] += 1
    widths = {cat: max(2, len(str(n))) for cat, n in counts.items()}

    mapping = {}
    counters = {"SCF": 0, "CTG": 0, "MIT": 0, "PLT": 0}
    for r in records:
        gwh_id = r["gwh_id"]
        cat, _ = classified[gwh_id]
        if cat == "CHR":
            value, literal = chrom_numbers[gwh_id]
            suffix = f"{value:0{chrom_width}d}" if value is not None else literal
            new_id = f"{prefix}C{suffix}"
        else:
            counters[cat] += 1
            new_id = f"{prefix}{cat}{counters[cat]:0{widths[cat]}d}"
        mapping[gwh_id] = new_id

    seen = {}
    for gwh_id, new_id in mapping.items():
        if new_id in seen:
            raise ValueError(
                f"Naming collision: both '{seen[new_id]}' and '{gwh_id}' would "
                f"rename to '{new_id}' — check the FASTA headers"
            )
        seen[new_id] = gwh_id
    return mapping, classified


def write_renamed_fasta(fasta_path: Path, mapping: dict, out_fh) -> None:
    with open(fasta_path) as in_fh:
        for line in in_fh:
            if line.startswith(">"):
                gwh_id = line[1:].split(None, 1)[0]
                out_fh.write(f">{mapping.get(gwh_id, gwh_id)}\n")
            else:
                out_fh.write(line)


def write_equiv_tsv(mapping: dict, classified: dict, records: list, out_path: Path) -> None:
    with open(out_path, "w") as fh:
        fh.write("gwh_id\toriseqid\tcategory\tnew_seqid\n")
        for r in records:
            gwh_id = r["gwh_id"]
            cat, _ = classified[gwh_id]
            fh.write(f"{gwh_id}\t{r['oriseqid']}\t{cat}\t{mapping[gwh_id]}\n")


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GWH4FastaRename",
        description="Rename sequence IDs in a FASTA file using GWH-format "
                    "headers already present in the file (no download).",
    )
    ap.add_argument("--fasta", required=True, type=Path,
                    help="Input FASTA file with GWH-format headers")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output renamed FASTA file (default: stdout)")
    ap.add_argument("--equiv", type=Path, default=None,
                    help="Equivalence TSV path (gwh_id, oriseqid, category, "
                         "new_seqid). Default: {output}.equiv_seqID.txt when "
                         "--output is given; skipped (with a note) when "
                         "writing to stdout and --equiv is not set")
    ap.add_argument("--prefix", default="Sp",
                    help="Prefix for renamed sequence IDs (default: Sp)")
    ap.add_argument("--dry_run", action="store_true",
                    help="Parse and classify the FASTA, print what would be "
                         "renamed, then exit without writing any output")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.fasta.exists():
        print(f"ERROR: --fasta file not found: {args.fasta}", file=sys.stderr)
        sys.exit(1)

    try:
        records = _scan_gwh_fasta(args.fasta)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)
    if not records:
        print(f"ERROR: no sequences found in {args.fasta}", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(records)} sequences from {args.fasta.name}", file=sys.stderr)

    try:
        mapping, classified = build_rename_mapping(records, args.prefix)
    except ValueError as e:
        print(f"ERROR: {e}", file=sys.stderr)
        sys.exit(1)

    cat_counts = {"CHR": 0, "SCF": 0, "CTG": 0, "MIT": 0, "PLT": 0}
    for cat, _ in classified.values():
        cat_counts[cat] += 1
    print(f"Classified: CHR={cat_counts['CHR']} SCF={cat_counts['SCF']} "
          f"CTG={cat_counts['CTG']} MIT={cat_counts['MIT']} PLT={cat_counts['PLT']}",
          file=sys.stderr)

    if args.dry_run:
        print("Dry run — no output will be written:", file=sys.stderr)
        for r in records:
            gwh_id = r["gwh_id"]
            cat, _ = classified[gwh_id]
            print(f"  {gwh_id} (OriSeqID={r['oriseqid']}, {cat}) -> {mapping[gwh_id]}",
                  file=sys.stderr)
        sys.exit(0)

    if args.output:
        with open(args.output, "w") as out_fh:
            write_renamed_fasta(args.fasta, mapping, out_fh)
        print(f"Written: {args.output}", file=sys.stderr)
        equiv_path = args.equiv or args.output.with_suffix(args.output.suffix + ".equiv_seqID.txt")
    else:
        write_renamed_fasta(args.fasta, mapping, sys.stdout)
        equiv_path = args.equiv

    if equiv_path:
        write_equiv_tsv(mapping, classified, records, equiv_path)
        print(f"Written: {equiv_path}", file=sys.stderr)
    else:
        print("NOTE: no --equiv path resolved (writing FASTA to stdout with no "
              "--output/--equiv given) — equivalence table not written",
              file=sys.stderr)

    print(f"Renamed: {len(mapping)} sequences", file=sys.stderr)


if __name__ == "__main__":
    main()
