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

Sequence ID renaming scheme (--rename_seqids), inferred from each FASTA
header's description text:
    chromosome  -> {prefix}C{NN}   (number from the description; Arabic
                                     digits used directly, Roman numerals
                                     auto-converted when the genome uses
                                     them, single-letter sex chromosomes
                                     X/Y/W/Z kept literal e.g. {prefix}CX)
    mitochondrion -> {prefix}MIT{NN}
    chloroplast/plastid -> {prefix}PLT{NN}
    scaffold    -> {prefix}SCF{NN}
    (anything else, incl. unplaced/unlocalized contigs) -> {prefix}CTG{NN}
Non-chromosome categories are numbered by order of appearance in the file.
Zero-padding width is the width of the largest number in that category
(minimum 2 digits).

Usage
-----
    NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/
    NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ \\
        --rename_seqids --rename_prefix Sp
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

VERSION = "v0.1.0"

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

CHROM_RE  = re.compile(r"chromosome[:\s]+([A-Za-z0-9]+)", re.IGNORECASE)
MITO_RE   = re.compile(r"mitochondri", re.IGNORECASE)
PLASTID_RE = re.compile(r"chloroplast|plastid", re.IGNORECASE)
SCAFFOLD_RE = re.compile(r"scaffold", re.IGNORECASE)


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

def _parse_fasta_headers(fasta_path: Path) -> list:
    """Returns list of (seq_id, description) in file order, streaming."""
    headers = []
    with open(fasta_path) as fh:
        for line in fh:
            if line.startswith(">"):
                header = line[1:].rstrip("\n")
                parts = header.split(None, 1)
                seq_id = parts[0]
                description = parts[1] if len(parts) > 1 else ""
                headers.append((seq_id, description))
    return headers


def _classify(headers: list) -> dict:
    """Classify each seq_id into a category + raw chromosome token (if any).
    Returns {seq_id: (category, raw_token_or_None)}."""
    classified = {}
    for seq_id, description in headers:
        m = CHROM_RE.search(description)
        if m:
            classified[seq_id] = ("C", m.group(1))
        elif MITO_RE.search(description):
            classified[seq_id] = ("MIT", None)
        elif PLASTID_RE.search(description):
            classified[seq_id] = ("PLT", None)
        elif SCAFFOLD_RE.search(description):
            classified[seq_id] = ("SCF", None)
        else:
            classified[seq_id] = ("CTG", None)
    return classified


def _resolve_chromosome_numbers(classified: dict) -> dict:
    """For category C, resolve each raw token to (numeric_value_or_None, literal).
    Detects genome-wide Roman-numeral usage before resolving ambiguous single
    letters (X/Y/W/Z), since a lone "X" could be the sex chromosome or Roman
    numeral 10 depending on how the rest of the genome is labeled."""
    chrom_tokens = [tok for cat, tok in classified.values() if cat == "C"]
    genome_is_roman = any(
        len(t) > 1 and ROMAN_RE.match(t.upper()) and not t.isdigit()
        for t in chrom_tokens
    )

    resolved = {}
    for seq_id, (cat, tok) in classified.items():
        if cat != "C":
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
    headers = _parse_fasta_headers(fasta_path)
    classified = _classify(headers)
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
    for seq_id, description in headers:
        cat, _ = classified[seq_id]
        if cat == "C":
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


def write_equiv_tsv(mapping: dict, headers: list, out_path: Path) -> None:
    order = [seq_id for seq_id, _ in headers]
    with open(out_path, "w") as fh:
        fh.write("# old_seqid\tnew_seqid\n")
        for seq_id in order:
            fh.write(f"{seq_id}\t{mapping[seq_id]}\n")


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


# ── Main per-accession pipeline ───────────────────────────────────────────

def process_accession(row: dict, output_dir: Path, rename_seqids: bool, force: bool,
                      strip_description: bool = False) -> bool:
    species, accession, prefix = row["species"], row["accession"], row["prefix"]
    acc_dir = output_dir / f"{species}_{accession}"
    fasta_out = acc_dir / f"{species}_{accession}.fasta"
    gff_out   = acc_dir / f"{species}_{accession}.gff3"
    equiv_out = acc_dir / f"{species}_{accession}.equiv_seqID.txt"

    if not force and fasta_out.exists() and fasta_out.stat().st_size > 0:
        print(f"[{species} / {accession}] already downloaded, skipping (use --force to redo)",
              file=sys.stderr)
        return True

    print(f"[{species} / {accession}] downloading...", file=sys.stderr)
    acc_dir.mkdir(parents=True, exist_ok=True)

    with tempfile.TemporaryDirectory(prefix=f"ncbi_dl_{accession}_") as tmp:
        tmp_path = Path(tmp)
        fasta_src, gff_src = download_genome(accession, tmp_path)
        if fasta_src is None:
            return False

        if rename_seqids:
            try:
                mapping = build_rename_mapping(fasta_src, prefix)
            except ValueError as e:
                print(f"  ERROR: {e}", file=sys.stderr)
                return False
            headers = _parse_fasta_headers(fasta_src)
            write_renamed_fasta(fasta_src, mapping, fasta_out, strip_description)
            write_equiv_tsv(mapping, headers, equiv_out)
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

    return True


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
    for row in rows:
        ok = process_accession(row, args.output, args.rename_seqids, args.force,
                              args.strip_description)
        if ok:
            n_ok += 1
        else:
            n_fail += 1

    print(f"Done: {n_ok} succeeded, {n_fail} failed (of {len(rows)} accessions)", file=sys.stderr)
    if n_fail:
        sys.exit(1)


if __name__ == "__main__":
    main()
