#!/usr/bin/env python3
"""
NCBI_DownloadGenome.py — Download genome FASTA (and GFF3, if available) from
NCBI for a list of accessions, with optional systematic sequence ID renaming.

Accessions file format (same as annotseba's accessions.txt):
    species<TAB>accession<TAB>taxa_id[<TAB>prefix]

Lines starting with '#' are comments. taxa_id and prefix may be omitted or
set to "NA"; taxa_id is accepted but unused by this script (kept for file
compatibility). prefix falls back to --rename_prefix when omitted/NA.

For each accession, downloads via the NCBI Datasets REST API v2 and writes,
under --output/{species}_{accession}/:
    {species}_{accession}.fasta
    {species}_{accession}.gff3            (only if annotation exists)
    {species}_{accession}.equiv_seqID.txt (only with --rename_seqids)

Sequence ID renaming scheme (--rename_seqids):
    CHR (chromosome/pseudomolecule/linkage group/LG) -> {prefix}C{NN}
        Number taken from the description as-is: Arabic digits used
        directly, Roman numerals auto-converted when the genome uses them,
        single-letter sex chromosomes X/Y/W/Z kept literal (e.g. {prefix}CX),
        and letter-suffixed polyploid labels (e.g. "1A", "2B") kept literal
        too since they aren't a plain number.
    MIT (mitochondrion) -> {prefix}MIT{NN}
    PLT (chloroplast/plastid) -> {prefix}PLT{NN}
    SCF (not a chromosome/organelle, sequence contains any N) -> {prefix}SCF{NN}
    CTG (not a chromosome/organelle, no N in the sequence) -> {prefix}CTG{NN}
Non-chromosome categories are numbered by order of appearance in the file,
not by any number in their own description. Zero-padding width is the width
of the largest number in that category (minimum 2 digits).

A FASTA scan (used for renaming and/or --report_metrics) also enforces that
every SeqID is unique, raising a clear error otherwise, and flags any
non-standard/ambiguous IUPAC nucleotide codes found (anything besides
A/C/G/T/N) since those can break aligners or variant callers that assume a
plain 4-letter(+N) alphabet.

--report_metrics writes {output}/summary.tsv (one row per accession: seq_n,
assembly_size, avg_length, n50, l50, n90, l90, per-category counts
CHR/SCF/CTG/MIT/PLT, ambiguous-nt count+characters, annotation YES/NO) and
prints the same as an ASCII table to stderr.

Usage
-----
    NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/
    NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ \\
        --rename_seqids --rename_prefix Sp
    NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ --report_metrics
"""

import argparse
import csv
import re
import shutil
import ssl
import sys
import tempfile
import time
import urllib.error
import urllib.request
import zipfile
from pathlib import Path

VERSION = "v0.3.0"

FASTA_LINE_WIDTH = 60
DOWNLOAD_RETRIES = 3
RETRY_WAIT_S = 5
CHUNK_SIZE = 1024 * 1024  # 1 MB

NCBI_DOWNLOAD_URL = (
    "https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/{acc}/download"
    "?include_annotation_type=GENOME_FASTA,GENOME_GFF"
)

ROMAN_RE = re.compile(r"^M{0,4}(CM|CD|D?C{0,3})(XC|XL|L?X{0,3})(IX|IV|V?I{0,3})$")
ROMAN_VALUES = [
    ("M", 1000), ("CM", 900), ("D", 500), ("CD", 400),
    ("C", 100), ("XC", 90), ("L", 50), ("XL", 40),
    ("X", 10), ("IX", 9), ("V", 5), ("IV", 4), ("I", 1),
]
SEX_CHROM_LETTERS = ("X", "Y", "W", "Z")

CHROM_RE = re.compile(
    r"(?:chromosome|pseudomolecule|linkage[\s_]?group)[:\s]+([A-Za-z0-9]+)",
    re.IGNORECASE,
)
CHROM_LG_ABBR_RE = re.compile(r"\bLG[\s_-]?(\d+[A-Za-z]*)\b", re.IGNORECASE)
MITO_RE   = re.compile(r"mitochondri", re.IGNORECASE)
PLASTID_RE = re.compile(r"chloroplast|plastid", re.IGNORECASE)


def _ssl_context_with_certifi_fallback():
    """Some Python installs (notably python.org's macOS installer) ship
    without a populated root CA bundle, causing every HTTPS request to fail
    certificate verification. Falls back to the optional 'certifi' package's
    bundle if importable, purely as a nicety — not a hard dependency."""
    try:
        import certifi
        return ssl.create_default_context(cafile=certifi.where())
    except ImportError:
        return None


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


# ── Accessions file ────────────────────────────────────────────────────────

def load_accessions(path: Path, default_prefix: str) -> list:
    """Returns list of dicts: species, accession, taxa_id, prefix."""
    rows = []
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            fields = line.split("\t")
            if len(fields) < 2:
                print(f"WARNING: skipping malformed line {lineno}: {line!r}", file=sys.stderr)
                continue
            species, accession = fields[0].strip(), fields[1].strip()
            taxa_id = fields[2].strip() if len(fields) > 2 and fields[2].strip().upper() != "NA" else ""
            prefix = (fields[3].strip() if len(fields) > 3 and fields[3].strip().upper() != "NA"
                     else default_prefix)
            rows.append({"species": species, "accession": accession,
                        "taxa_id": taxa_id, "prefix": prefix})
    return rows


# ── Download ───────────────────────────────────────────────────────────────

def download_genome(accession: str, workdir: Path) -> tuple:
    """Download and extract genome+annotation zip for accession into workdir.
    Returns (fasta_path_or_None, gff_path_or_None)."""
    url = NCBI_DOWNLOAD_URL.format(acc=accession)
    zip_path = workdir / f"{accession}.zip"

    ssl_context = None  # None = urllib's normal default verification behavior
    last_err = None
    for attempt in range(1, DOWNLOAD_RETRIES + 1):
        try:
            with urllib.request.urlopen(url, timeout=120, context=ssl_context) as resp, \
                 open(zip_path, "wb") as out_fh:
                while True:
                    chunk = resp.read(CHUNK_SIZE)
                    if not chunk:
                        break
                    out_fh.write(chunk)
            last_err = None
            break
        except urllib.error.URLError as e:
            if isinstance(e.reason, ssl.SSLCertVerificationError) and ssl_context is None:
                fallback = _ssl_context_with_certifi_fallback()
                if fallback is not None:
                    print(f"  NOTE: this Python's default certificate store looks broken "
                          f"(common on python.org macOS installs) — retrying with the "
                          f"'certifi' package's bundle instead", file=sys.stderr)
                    ssl_context = fallback
                    continue
                last_err = e
                print(f"  ERROR: SSL certificate verification failed and 'certifi' is not "
                      f"installed to fall back on. Fix your Python installation's "
                      f"certificate store (on macOS python.org installs, run "
                      f"'/Applications/Python */Install Certificates.command'), or "
                      f"'pip install certifi' as a workaround.", file=sys.stderr)
                break
            last_err = e
            print(f"  WARNING: download attempt {attempt}/{DOWNLOAD_RETRIES} failed "
                  f"for {accession}: {e}", file=sys.stderr)
            if attempt < DOWNLOAD_RETRIES:
                time.sleep(RETRY_WAIT_S)
        except TimeoutError as e:
            last_err = e
            print(f"  WARNING: download attempt {attempt}/{DOWNLOAD_RETRIES} failed "
                  f"for {accession}: {e}", file=sys.stderr)
            if attempt < DOWNLOAD_RETRIES:
                time.sleep(RETRY_WAIT_S)

    if last_err is not None:
        print(f"  ERROR: could not download {accession} after {DOWNLOAD_RETRIES} attempts", file=sys.stderr)
        return None, None

    extract_dir = workdir / f"{accession}_extract"
    extract_dir.mkdir(exist_ok=True)
    try:
        with zipfile.ZipFile(zip_path) as zf:
            zf.extractall(extract_dir)
    except zipfile.BadZipFile:
        print(f"  ERROR: {accession} — response was not a valid zip "
              f"(accession may not exist)", file=sys.stderr)
        zip_path.unlink(missing_ok=True)
        shutil.rmtree(extract_dir, ignore_errors=True)
        return None, None

    fasta_matches = sorted(extract_dir.rglob("*.fna"))
    gff_matches   = sorted(extract_dir.rglob("*.gff"))

    fasta_path = fasta_matches[0] if fasta_matches else None
    gff_path   = gff_matches[0] if gff_matches else None

    if fasta_path is None:
        print(f"  WARNING: no genome FASTA found in {accession} download", file=sys.stderr)
    if gff_path is None:
        print(f"  NOTE: no GFF3 annotation available for {accession}", file=sys.stderr)

    return fasta_path, gff_path


# ── SeqID classification / renaming ───────────────────────────────────────

IUPAC_AMBIGUITY_CODES = "RYSWKMBDHV"  # standard = A/C/G/T/N only; anything
                                       # else (these + truly unexpected chars)
                                       # can break aligners/callers expecting
                                       # a plain 4-letter (+N) alphabet.


def _scan_fasta(fasta_path: Path) -> list:
    """Single streaming pass over the FASTA. Returns a list of dicts
    (seq_id, description, length, has_n, ambig_count, ambig_chars) in file
    order. Never holds a full sequence in memory — has_n is tracked
    incrementally and stops checking once an N is seen, so this stays a
    per-line scan even on long records. Raises ValueError if any seq_id is
    duplicated (a data-integrity problem that would otherwise silently
    corrupt classification/metrics/renaming, since results are keyed by
    seq_id)."""
    records = []
    seq_id, description, length, has_n = None, "", 0, False
    ambig_count, ambig_chars = 0, set()

    def flush():
        records.append({"seq_id": seq_id, "description": description,
                        "length": length, "has_n": has_n,
                        "ambig_count": ambig_count, "ambig_chars": ambig_chars})

    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                if seq_id is not None:
                    flush()
                header = line[1:].rstrip("\n")
                parts = header.split(None, 1)
                seq_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                length, has_n = 0, False
                ambig_count, ambig_chars = 0, set()
            else:
                seq = line.rstrip("\n")
                length += len(seq)
                if not has_n and ("N" in seq or "n" in seq):
                    has_n = True
                seq_upper = seq.upper()
                for c in IUPAC_AMBIGUITY_CODES:
                    n = seq_upper.count(c)
                    if n:
                        ambig_count += n
                        ambig_chars.add(c)
    if seq_id is not None:
        flush()

    seen_ids = {}
    duplicates = []
    for r in records:
        if r["seq_id"] in seen_ids:
            duplicates.append(r["seq_id"])
        seen_ids[r["seq_id"]] = True
    if duplicates:
        examples = ", ".join(sorted(set(duplicates))[:5])
        raise ValueError(
            f"{fasta_path.name} has {len(set(duplicates))} duplicate SeqID(s), "
            f"e.g. {examples} — every sequence in a FASTA must have a unique ID"
        )

    return records


def _classify_description(description: str) -> tuple:
    """Returns (category, raw_chrom_token_or_None) from description text alone
    (chromosome/pseudomolecule/linkage group/mitochondrion/chloroplast).
    Returns (None, None) when none of those match — caller resolves the
    remaining SCF-vs-CTG split from actual sequence content (has_n)."""
    m = CHROM_RE.search(description) or CHROM_LG_ABBR_RE.search(description)
    if m:
        return "CHR", m.group(1)
    if MITO_RE.search(description):
        return "MIT", None
    if PLASTID_RE.search(description):
        return "PLT", None
    return None, None


def classify_records(records: list) -> dict:
    """Returns {seq_id: (category, raw_chrom_token_or_None)}, category one of
    CHR/MIT/PLT/SCF/CTG. Non-chromosome/organelle sequences are SCF if they
    contain any N (gap-containing, typical of scaffolds), else CTG."""
    classified = {}
    for r in records:
        cat, tok = _classify_description(r["description"])
        if cat is None:
            cat = "SCF" if r["has_n"] else "CTG"
        classified[r["seq_id"]] = (cat, tok)
    return classified


def _resolve_chromosome_numbers(classified: dict) -> dict:
    """For category CHR, resolve each raw token to (numeric_value_or_None,
    literal). Detects genome-wide Roman-numeral usage before resolving
    ambiguous single letters (X/Y/W/Z), since a lone "X" could be the sex
    chromosome or Roman numeral 10 depending on how the rest of the genome
    is labeled. Letter-suffixed polyploid labels (e.g. "1A", "2B", common in
    wheat-like subgenomes) aren't digit-only, valid Roman numerals, or a
    lone sex-chromosome letter, so they fall through and are kept literal —
    exactly the "keep the original numbering" behavior wanted for those."""
    chrom_tokens = [tok for cat, tok in classified.values() if cat == "CHR"]
    genome_is_roman = any(
        len(t) > 1 and ROMAN_RE.match(t.upper()) and not t.isdigit()
        for t in chrom_tokens
    )

    resolved = {}
    for seq_id, (cat, tok) in classified.items():
        if cat != "CHR":
            continue
        if tok.isdigit():
            resolved[seq_id] = (int(tok), tok)
        elif len(tok) == 1 and tok.upper() in SEX_CHROM_LETTERS and not genome_is_roman:
            resolved[seq_id] = (None, tok.upper())
        elif ROMAN_RE.match(tok.upper()) and (genome_is_roman or len(tok) > 1):
            resolved[seq_id] = (_roman_to_int(tok.upper()), tok.upper())
        else:
            resolved[seq_id] = (None, tok)
    return resolved


def build_rename_mapping(fasta_path: Path, prefix: str) -> dict:
    """Returns {old_seq_id: new_seq_id}, raising ValueError on a naming collision."""
    records = _scan_fasta(fasta_path)
    classified = classify_records(records)
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
        seq_id = r["seq_id"]
        cat, _ = classified[seq_id]
        if cat == "CHR":
            value, literal = chrom_numbers[seq_id]
            suffix = f"{value:0{chrom_width}d}" if value is not None else literal
            new_id = f"{prefix}C{suffix}"
        else:
            counters[cat] += 1
            new_id = f"{prefix}{cat}{counters[cat]:0{widths[cat]}d}"
        mapping[seq_id] = new_id

    seen = {}
    for old_id, new_id in mapping.items():
        if new_id in seen:
            raise ValueError(
                f"Naming collision: both '{seen[new_id]}' and '{old_id}' would "
                f"rename to '{new_id}' — check the FASTA descriptions"
            )
        seen[new_id] = old_id
    return mapping


def write_renamed_fasta(fasta_path: Path, mapping: dict, out_path: Path,
                        strip_description: bool = False) -> None:
    with open(fasta_path) as in_fh, open(out_path, "w") as out_fh:
        for line in in_fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                parts = header.split(None, 1)
                seq_id = parts[0]
                if strip_description:
                    description = ""
                else:
                    description = (" " + parts[1]) if len(parts) > 1 else ""
                out_fh.write(f">{mapping.get(seq_id, seq_id)}{description}\n")
            else:
                out_fh.write(line)


def write_equiv_tsv(mapping: dict, records: list, out_path: Path) -> None:
    with open(out_path, "w") as fh:
        fh.write("# old_seqid\tnew_seqid\n")
        for r in records:
            fh.write(f"{r['seq_id']}\t{mapping[r['seq_id']]}\n")


def write_renamed_gff3(gff_path: Path, mapping: dict, out_path: Path) -> None:
    warned = set()
    with open(gff_path) as in_fh, open(out_path, "w") as out_fh:
        for line in in_fh:
            if line.startswith("##sequence-region"):
                parts = line.rstrip("\n").split(" ")
                if len(parts) > 1 and parts[1] in mapping:
                    parts[1] = mapping[parts[1]]
                out_fh.write(" ".join(parts) + "\n")
            elif line.startswith("#"):
                out_fh.write(line)
            else:
                fields = line.rstrip("\n").split("\t")
                if len(fields) >= 1:
                    if fields[0] in mapping:
                        fields[0] = mapping[fields[0]]
                    elif fields[0] not in warned:
                        print(f"  WARNING: GFF3 seqid '{fields[0]}' not in FASTA "
                              f"rename mapping, left unchanged", file=sys.stderr)
                        warned.add(fields[0])
                    out_fh.write("\t".join(fields) + "\n")
                else:
                    out_fh.write(line)


# ── Assembly metrics (--report_metrics) ───────────────────────────────────

METRIC_COLUMNS = ["species", "accession", "seq_n", "assembly_size", "avg_length",
                  "n50", "l50", "n90", "l90", "n_chr", "n_scf", "n_ctg", "n_mit",
                  "n_plt", "ambig_nt_count", "ambig_nt_chars", "annotation_gff"]


def _nx_lx(lengths_desc: list, total: int, fraction: float) -> tuple:
    """Return (Nx, Lx) for a given fraction (0-1). lengths_desc sorted descending."""
    target = total * fraction
    cumsum = 0
    for i, length in enumerate(lengths_desc, 1):
        cumsum += length
        if cumsum >= target:
            return length, i
    return lengths_desc[-1], len(lengths_desc)


def compute_simple_metrics(fasta_path: Path) -> dict:
    """Uses _scan_fasta (which also enforces unique SeqIDs) to compute
    seq_n/assembly_size/avg_length/n50/l50/n90/l90, per-category counts
    (CHR/SCF/CTG/MIT/PLT), and non-standard (ambiguous IUPAC) nt counts,
    all in the same pass. Raises ValueError on duplicate SeqIDs."""
    records = _scan_fasta(fasta_path)

    empty = {"seq_n": 0, "assembly_size": 0, "avg_length": 0, "n50": 0, "l50": 0,
             "n90": 0, "l90": 0, "n_chr": 0, "n_scf": 0, "n_ctg": 0, "n_mit": 0,
             "n_plt": 0, "ambig_nt_count": 0, "ambig_nt_chars": ""}
    if not records:
        return empty

    lengths = [r["length"] for r in records]
    lengths_desc = sorted(lengths, reverse=True)
    total = sum(lengths)
    n50, l50 = _nx_lx(lengths_desc, total, 0.50)
    n90, l90 = _nx_lx(lengths_desc, total, 0.90)

    classified = classify_records(records)
    cat_counts = {"CHR": 0, "SCF": 0, "CTG": 0, "MIT": 0, "PLT": 0}
    for cat, _ in classified.values():
        cat_counts[cat] += 1

    # A FASTA with no header descriptions at all (e.g. a previously
    # --strip_description'd output being re-scanned on a skip-existing run)
    # can't be classified into CHR/MIT/PLT from text — everything falls to
    # SCF/CTG purely by N-content, which would silently misreport those counts.
    if all(not r["description"] for r in records) and not (cat_counts["CHR"]
            or cat_counts["MIT"] or cat_counts["PLT"]):
        print(f"  WARNING: {fasta_path.name} has no header descriptions (likely "
              f"already renamed with --strip_description) — CHR/MIT/PLT counts "
              f"can't be recovered from this file; use --force to recompute "
              f"from a fresh download", file=sys.stderr)

    ambig_count = sum(r["ambig_count"] for r in records)
    ambig_chars = sorted(set().union(*(r["ambig_chars"] for r in records)))

    return {
        "seq_n": len(lengths),
        "assembly_size": total,
        "avg_length": round(total / len(lengths), 2),
        "n50": n50, "l50": l50,
        "n90": n90, "l90": l90,
        "n_chr": cat_counts["CHR"], "n_scf": cat_counts["SCF"],
        "n_ctg": cat_counts["CTG"], "n_mit": cat_counts["MIT"], "n_plt": cat_counts["PLT"],
        "ambig_nt_count": ambig_count,
        "ambig_nt_chars": ",".join(ambig_chars),
    }


def _ascii_table(headers: list, rows: list) -> list:
    """Return list of strings forming an ASCII box table."""
    widths = [max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
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


def write_metrics_tsv(rows: list, path: Path) -> None:
    with open(path, "w") as fh:
        fh.write("\t".join(METRIC_COLUMNS) + "\n")
        for row in rows:
            fh.write("\t".join(str(row[c]) for c in METRIC_COLUMNS) + "\n")


def print_metrics_report(rows: list) -> None:
    print("", file=sys.stderr)
    print("Assembly metrics summary:", file=sys.stderr)
    headers = ["Species", "Accession", "Seq_N", "Assembly_size", "Avg_length",
              "N50", "L50", "N90", "L90", "CHR", "SCF", "CTG", "MIT", "PLT",
              "Ambig_NT", "Ambig_chars", "Annotation (GFF)"]
    table_rows = [[str(row[c]) for c in METRIC_COLUMNS] for row in rows]
    for line in _ascii_table(headers, table_rows):
        print(line, file=sys.stderr)


# ── Main per-accession pipeline ───────────────────────────────────────────

def _metrics_row(species: str, accession: str, fasta_out: Path, gff_out: Path):
    """Returns a metrics dict, or None (with an error printed) if fasta_out
    fails the duplicate-SeqID sanity check."""
    try:
        metrics = compute_simple_metrics(fasta_out)
    except ValueError as e:
        print(f"  ERROR: {e}", file=sys.stderr)
        return None
    row = {"species": species, "accession": accession}
    row.update(metrics)
    row["annotation_gff"] = "YES" if gff_out.exists() and gff_out.stat().st_size > 0 else "NO"
    return row


def process_accession(row: dict, output_dir: Path, rename_seqids: bool, force: bool,
                      strip_description: bool = False, report_metrics: bool = False) -> tuple:
    """Returns (ok: bool, metrics_row: dict | None)."""
    species, accession, prefix = row["species"], row["accession"], row["prefix"]
    acc_dir = output_dir / f"{species}_{accession}"
    fasta_out = acc_dir / f"{species}_{accession}.fasta"
    gff_out   = acc_dir / f"{species}_{accession}.gff3"
    equiv_out = acc_dir / f"{species}_{accession}.equiv_seqID.txt"

    if not force and fasta_out.exists() and fasta_out.stat().st_size > 0:
        print(f"[{species} / {accession}] already downloaded, skipping (use --force to redo)",
              file=sys.stderr)
        metrics = _metrics_row(species, accession, fasta_out, gff_out) if report_metrics else None
        return True, metrics

    print(f"[{species} / {accession}] downloading...", file=sys.stderr)
    acc_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"ncbi_dl_{accession}_") as tmp:
        tmp_path = Path(tmp)
        fasta_src, gff_src = download_genome(accession, tmp_path)
        if fasta_src is None:
            return False, None

        if rename_seqids:
            try:
                mapping = build_rename_mapping(fasta_src, prefix)
            except ValueError as e:
                print(f"  ERROR: {e}", file=sys.stderr)
                return False, None
            records = _scan_fasta(fasta_src)
            write_renamed_fasta(fasta_src, mapping, fasta_out, strip_description)
            write_equiv_tsv(mapping, records, equiv_out)
            print(f"  Written: {fasta_out.name}, {equiv_out.name} "
                  f"({len(mapping)} sequences renamed)", file=sys.stderr)
            if gff_src is not None:
                write_renamed_gff3(gff_src, mapping, gff_out)
                print(f"  Written: {gff_out.name}", file=sys.stderr)
        else:
            shutil.copy(fasta_src, fasta_out)
            print(f"  Written: {fasta_out.name}", file=sys.stderr)
            if gff_src is not None:
                shutil.copy(gff_src, gff_out)
                print(f"  Written: {gff_out.name}", file=sys.stderr)

        # Classify/measure from fasta_src (pre-rename, always has full header
        # descriptions) rather than fasta_out, so --strip_description can't
        # blind CHR/MIT/PLT classification on a fresh download — only the
        # skip-existing path (no fasta_src available) falls back to fasta_out.
        metrics = _metrics_row(species, accession, fasta_src, gff_out) if report_metrics else None

    return True, metrics


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="NCBI_DownloadGenome",
        description="Download genome FASTA (and GFF3, if available) from NCBI "
                    "for a list of accessions, with optional SeqID renaming.",
    )
    ap.add_argument("--accessions", required=True, type=Path,
                    help="Accessions file: species<TAB>accession<TAB>taxa_id[<TAB>prefix] "
                         "(same format as annotseba's accessions.txt)")
    ap.add_argument("--output", required=True, type=Path,
                    help="Output directory (one subdirectory per accession)")
    ap.add_argument("--rename_seqids", action="store_true",
                    help="Rename sequence IDs to {prefix}C##/SCF##/CTG##/MIT##/PLT## "
                         "based on each sequence's description, and apply the same "
                         "renaming to the GFF3 if one was downloaded")
    ap.add_argument("--rename_prefix", default="Sp",
                    help="Default prefix for renamed sequence IDs when not set in the "
                         "accessions file (default: Sp)")
    ap.add_argument("--strip_description", action="store_true",
                    help="Drop the FASTA header description text when renaming "
                         "(requires --rename_seqids), leaving just '>{new_id}'")
    ap.add_argument("--report_metrics", action="store_true",
                    help="Compute simple assembly metrics (seq_n, assembly_size, "
                         "avg_length, N50/L50, N90/L90, CHR/SCF/CTG/MIT/PLT counts, "
                         "ambiguous-nt count, GFF annotation YES/NO) for each "
                         "processed accession; writes {output}/summary.tsv and "
                         "prints the same as an ASCII table")
    ap.add_argument("--force", action="store_true",
                    help="Re-download and re-process even if output already exists")
    ap.add_argument("--dry_run", action="store_true",
                    help="Validate the accessions file and print what would be "
                         "downloaded, then exit without any network calls")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.accessions.exists():
        print(f"ERROR: --accessions file not found: {args.accessions}", file=sys.stderr)
        sys.exit(1)

    if args.strip_description and not args.rename_seqids:
        print("WARNING: --strip_description has no effect without --rename_seqids "
              "(FASTA is copied through unmodified otherwise)", file=sys.stderr)

    rows = load_accessions(args.accessions, args.rename_prefix)
    if not rows:
        print(f"ERROR: no accessions parsed from {args.accessions}", file=sys.stderr)
        sys.exit(1)
    print(f"Loaded {len(rows)} accessions from {args.accessions.name}", file=sys.stderr)

    if args.dry_run:
        print("Dry run — no downloads will be executed:", file=sys.stderr)
        for row in rows:
            print(f"  [{row['species']} / {row['accession']}] prefix={row['prefix']!r} "
                  f"-> {args.output}/{row['species']}_{row['accession']}/", file=sys.stderr)
        sys.exit(0)

    args.output.mkdir(parents=True, exist_ok=True)

    n_ok, n_fail = 0, 0
    metrics_rows = []
    for row in rows:
        ok, metrics = process_accession(row, args.output, args.rename_seqids, args.force,
                                        args.strip_description, args.report_metrics)
        if ok:
            n_ok += 1
            if metrics is not None:
                metrics_rows.append(metrics)
        else:
            n_fail += 1

    print(f"Done: {n_ok} succeeded, {n_fail} failed (of {len(rows)} accessions)", file=sys.stderr)

    if args.report_metrics and metrics_rows:
        summary_path = args.output / "summary.tsv"
        write_metrics_tsv(metrics_rows, summary_path)
        print(f"Written: {summary_path}", file=sys.stderr)
        print_metrics_report(metrics_rows)

    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
