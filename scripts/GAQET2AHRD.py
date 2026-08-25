#!/usr/bin/env python3
"""
GAQET2AHRD.py — Build an AHRD config from a GAQET run and (optionally) run
AHRD (https://github.com/groupschoof/AHRD) on it.

GAQET (https://github.com/victorgcb1987/GAQET2) writes a GAQET.log.txt inside
its run directory that records the exact diamond blastp commands used for
the TREMBL and SWISSPROT homology searches, e.g.:

    #TREMBL command used:
        diamond blastp --threads 48 --db .../uniprot_trembl_r2025_01.dmnd \\
            --query .../input_sequences/genome.proteins.fasta \\
            --out .../DIAMOND_run/genome.proteins.dmd.TREMBL.o6.txt

    #SWISSPROT command used:
        diamond blastp --threads 48 --db .../uniprot_sprot_r2025_01.dmnd \\
            --query .../input_sequences/genome.proteins.fasta \\
            --out .../DIAMOND_run/genome.proteins.dmd.SWISSPROT.o6.txt

This script parses those two commands out of GAQET.log.txt, derives the
FASTA-with-headers counterpart of each diamond .dmnd database (AHRD needs
the flat FASTA to read description lines, not the diamond index), writes an
AHRD YAML config from them, and by default invokes AHRD via
`java -jar $AHRD_JAR config.yml` (jar path from --ahrd_jar or the AHRD_JAR
environment variable). Use --skip_ahrd to only write the config.

If --ahrd_home (or the AHRD_HOME environment variable) is given, the
blacklist/filter/token_blacklist word-filtering files bundled with AHRD's
own distribution under {ahrd_home}/test/resources/ are wired into the
config automatically (filter is per-database: filter_descline_sprot.txt vs
filter_descline_trembl.txt; blacklist and token_blacklist are shared).
Any of --blacklist/--swissprot_filter/--trembl_filter/--token_blacklist can
override these individually.

GO term transfer is on by default: gene_ontology_result defaults to
goa_uniprot_all.gaf in the same directory as the SWISSPROT --db (the
standard layout alongside uniprot_sprot*.dmnd), override with
--gene_ontology_result or disable with --skip_go.

After a successful AHRD run, a <prefix>_AHRD.summary.txt is written next to
the AHRD output TSV (counts of proteins with a description/GO term(s) vs.
unknown, the AHRD-Quality-Code distribution, and the --top_n most abundant
descriptions), and the same summary is printed to stderr as ASCII tables.
Use --skip_summary to omit both.

--check_te_goterms cross-checks AHRD's GO terms against DETENGA's
dedicated TE calls, to catch
protein-coding gene models that are actually TE-derived. For each protein:
a GOTE_TAGGED of YES means every one of its GO terms is in the
TE-associated set (--te_goterms_file, default: a hardcoded list -- see
--print_te_associated_default_goterms); NO means it has GO terms but not
all are TE-associated; NA means it has no GO terms at all (no info).
DETENGA_TAGGED mirrors this using DETENGA's own combined call
(DeTEnGA_status containing "te" in either half, e.g. PteM0/P0Mte/PteMte,
means TE; PcpM0 means confirmed non-TE; a protein missing from DETENGA's
--detenga_csv entirely gets NA). Writes
<prefix>_TEGOterm_vs_DETENGA.tsv (ProteinID, AHRD_GO_TEs, DETENGA_TE,
GOTE_TAGGED, DETENGA_TAGGED) and prints an agreement summary.

Usage
-----
    GAQET2AHRD.py --gaqet_log GAQET.log.txt --ahrd_jar /opt/ahrd/ahrd.jar \\
        --ahrd_home /opt/ahrd
    AHRD_JAR=/opt/ahrd/ahrd.jar AHRD_HOME=/opt/ahrd \\
        GAQET2AHRD.py --gaqet_log GAQET.log.txt
    GAQET2AHRD.py --gaqet_log GAQET.log.txt --skip_ahrd
    GAQET2AHRD.py --gaqet_log GAQET.log.txt --dry_run
"""

import argparse
import glob
import os
import re
import subprocess
import sys
from collections import Counter
from pathlib import Path

VERSION = "v0.3.0"

DEFAULT_TOP_N = 10

DEFAULT_SWISSPROT_WEIGHT = 653
DEFAULT_TREMBL_WEIGHT = 904
DEFAULT_SWISSPROT_DESC_SCORE_BIT_SCORE_WEIGHT = 2.717061
DEFAULT_TREMBL_DESC_SCORE_BIT_SCORE_WEIGHT = 2.590211
DEFAULT_TOKEN_SCORE_BIT_SCORE_WEIGHT = 0.468
DEFAULT_TOKEN_SCORE_DATABASE_SCORE_WEIGHT = 0.2098
DEFAULT_TOKEN_SCORE_OVERLAP_SCORE_WEIGHT = 0.3221

# UniProt GOA GAF format: column 1 "UniProtKB", column 2 accession, column 4
# qualifier (excluded when it's a negation, "NOT|..."), column 5 GO term.
DEFAULT_REFERENCE_GO_REGEX = r"^UniProtKB\\t(?<shortAccession>[^\\t]+)\\t[^\\t]+\\t(?!NOT\\|)[^\\t]*\\t(?<goTerm>GO:\\d{7})"

# GO terms associated with transposable elements (transposase/reverse-
# transcriptase/integrase activity, transposition processes, retrotransposon
# nucleocapsid) -- used by --check_te_goterms, override with
# --te_goterms_file. (go_id, name, category)
DEFAULT_TE_GOTERMS = [
    ("GO:0004803", "transposase activity", "molecular_function"),
    ("GO:0003964", "RNA-directed DNA polymerase activity", "molecular_function"),
    ("GO:0004523", "RNA-DNA hybrid ribonuclease activity", "molecular_function"),
    ("GO:0004190", "aspartic-type endopeptidase activity", "molecular_function"),
    ("GO:0004521", "RNA endonuclease activity", "molecular_function"),
    ("GO:0006313", "transposition, DNA-mediated", "biological_process"),
    ("GO:0032197", "transposition, RNA-mediated", "biological_process"),
    ("GO:0015074", "DNA integration", "biological_process"),
    ("GO:0000943", "retrotransposon nucleocapsid", "cellular_component"),
    ("GO:0032068", "retrotransposon assembly", "cellular_component"),
]


def parse_gaqet_log(path: Path) -> dict:
    """Extract the TREMBL/SWISSPROT diamond commands from GAQET.log.txt.
    Returns a dict with proteins_fasta, swissprot_db, swissprot_out,
    trembl_db, trembl_out. Aborts with a clear error if either command
    (or any of its --db/--query/--out flags) can't be found."""
    text = path.read_text()

    parsed = {}
    proteins_fastas = set()
    for label, key in (("TREMBL", "trembl"), ("SWISSPROT", "swissprot")):
        cmd_match = re.search(
            rf"#{label} command used:\s*\n\s*(.+)", text, re.IGNORECASE
        )
        if not cmd_match:
            print(f"ERROR: could not find a '#{label} command used:' block "
                  f"in {path}", file=sys.stderr)
            sys.exit(1)
        cmd = cmd_match.group(1)

        db_match = re.search(r"--db\s+(\S+)", cmd)
        query_match = re.search(r"--query\s+(\S+)", cmd)
        out_match = re.search(r"--out\s+(\S+)", cmd)
        if not (db_match and query_match and out_match):
            print(f"ERROR: '#{label} command used:' line found but missing "
                  f"--db/--query/--out: {cmd}", file=sys.stderr)
            sys.exit(1)

        success_window = text[cmd_match.end():cmd_match.end() + 300]
        if "run successfully" not in success_window:
            print(f"WARNING: no success confirmation found after the "
                  f"{label} command in {path.name} — proceeding anyway",
                  file=sys.stderr)

        parsed[f"{key}_db"] = Path(db_match.group(1))
        parsed[f"{key}_out"] = Path(out_match.group(1))
        proteins_fastas.add(query_match.group(1))

    if len(proteins_fastas) != 1:
        print(f"ERROR: TREMBL and SWISSPROT commands use different --query "
              f"FASTA files ({sorted(proteins_fastas)}) — expected the same "
              f"protein FASTA for both", file=sys.stderr)
        sys.exit(1)
    parsed["proteins_fasta"] = Path(proteins_fastas.pop())
    return parsed


def derive_db_fasta(dmnd_path: Path) -> Path:
    """AHRD needs the flat FASTA (description headers) of each blast DB,
    not the diamond .dmnd index. Derived by swapping the extension."""
    return dmnd_path.with_suffix(".fasta")


def build_ahrd_config(
    proteins_fasta: Path, output_tsv: Path,
    swissprot_out: Path, swissprot_fasta: Path,
    trembl_out: Path, trembl_fasta: Path,
    swissprot_weight: int, trembl_weight: int,
    swissprot_desc_weight: float, trembl_desc_weight: float,
    token_bit_score_weight: float, token_database_score_weight: float,
    token_overlap_score_weight: float,
    blacklist: Path, swissprot_filter: Path, trembl_filter: Path,
    token_blacklist: Path, gene_ontology_result: Path,
) -> str:
    def db_block(name: str, weight: int, desc_weight: float,
                 out_file: Path, db_fasta: Path, filter_file: Path) -> str:
        lines = [
            f"  {name}:",
            f"    weight: {weight}",
            f"    description_score_bit_score_weight: {desc_weight}",
            f"    file: {out_file}",
            f"    database: {db_fasta}",
        ]
        if blacklist:
            lines.append(f"    blacklist: {blacklist}")
        if filter_file:
            lines.append(f"    filter: {filter_file}")
        if token_blacklist:
            lines.append(f"    token_blacklist: {token_blacklist}")
        return "\n".join(lines)

    go_block = ""
    if gene_ontology_result:
        go_block = (
            f"gene_ontology_result: {gene_ontology_result}\n"
            f'reference_go_regex: "{DEFAULT_REFERENCE_GO_REGEX}"\n'
            f"prefer_reference_with_go_annos: true\n"
        )

    return (
        f"proteins_fasta: {proteins_fasta}\n"
        f"token_score_bit_score_weight: {token_bit_score_weight}\n"
        f"token_score_database_score_weight: {token_database_score_weight}\n"
        f"token_score_overlap_score_weight: {token_overlap_score_weight}\n"
        f"{go_block}"
        f"output: {output_tsv}\n"
        f"blast_dbs:\n"
        f"{db_block('swissprot', swissprot_weight, swissprot_desc_weight, swissprot_out, swissprot_fasta, swissprot_filter)}\n"
        f"{db_block('trembl', trembl_weight, trembl_desc_weight, trembl_out, trembl_fasta, trembl_filter)}\n"
    )


def find_gaqet_prefix(gaqet_log: Path) -> str:
    """GAQET names its per-run stats file {prefix}_GAQET.stats.tsv in the
    same directory as GAQET.log.txt; used only to name AHRD's own output
    files. Falls back to 'AHRD' if no such file is found."""
    hits = glob.glob(str(gaqet_log.parent / "*_GAQET.stats.tsv"))
    if not hits:
        return "AHRD"
    return Path(hits[0]).name[: -len("_GAQET.stats.tsv")]


def _iter_ahrd_rows(tsv_path: Path):
    """Yield (protein_id, description, go_ids, quality_code) tuples from an
    AHRD output TSV. Column positions (besides protein_id, always column 0)
    are located by header keyword rather than fixed index, since the exact
    column set depends on the AHRD config (e.g. whether GO/InterPro were
    requested); AHRD's leading '# AHRD-Version ...' comment line is
    skipped. go_ids is a list of GO:####### tokens extracted from the GO
    column (handles both 'GO:XXXXXXX (name), GO:YYYYYYY (name)' and plain
    comma-separated forms)."""
    with open(tsv_path) as fh:
        idx = {}
        for line in fh:
            line = line.rstrip("\n")
            if not line or line.startswith("#"):
                continue
            cols = line.split("\t")
            if not idx:
                for i, col in enumerate(cols):
                    key = col.strip().lower()
                    if "quality-code" in key:
                        idx["quality"] = i
                    elif "gene-ontology" in key or "go-id" in key or key.startswith("go"):
                        idx["go"] = i
                    elif "description" in key and "interpro" not in key:
                        # AHRD's TSV also has an "Interpro-ID (Description)"
                        # column (empty unless InterPro was configured) --
                        # excluded so it can't clobber the real
                        # "Human-Readable-Description" column below.
                        idx["description"] = i
                continue

            protein_id = cols[0].strip() if cols else ""
            desc = cols[idx["description"]].strip() if "description" in idx and idx["description"] < len(cols) else ""
            go_raw = cols[idx["go"]].strip() if "go" in idx and idx["go"] < len(cols) else ""
            quality = cols[idx["quality"]].strip() if "quality" in idx and idx["quality"] < len(cols) else ""
            go_ids = re.findall(r"GO:\d{7}", go_raw)
            yield protein_id, desc, go_ids, quality


def summarize_ahrd_output(tsv_path: Path) -> dict:
    """Compute summary metrics over an AHRD output TSV. Returns a dict with
    n_total, n_with_desc, n_unknown, n_with_go, quality_counts (Counter),
    desc_counts (Counter)."""
    n_total = 0
    n_with_desc = 0
    n_with_go = 0
    quality_counts = Counter()
    desc_counts = Counter()

    for _protein_id, desc, go_ids, quality in _iter_ahrd_rows(tsv_path):
        n_total += 1
        if desc and not desc.lower().startswith("unknown protein"):
            n_with_desc += 1
            desc_counts[desc] += 1
        if go_ids:
            n_with_go += 1
        if quality:
            quality_counts[quality] += 1

    return {
        "n_total": n_total,
        "n_with_desc": n_with_desc,
        "n_unknown": n_total - n_with_desc,
        "n_with_go": n_with_go,
        "quality_counts": quality_counts,
        "desc_counts": desc_counts,
    }


def _ascii_table(headers: list, rows: list, left_align: set = None) -> list:
    """Return list of strings forming an ASCII box table. Column 0 is
    left-aligned by default (and any extra indices in left_align); the
    rest are right-aligned, matching the usual text-then-numbers layout."""
    left_align = left_align if left_align is not None else {0}
    widths = [max(len(h), *(len(r[i]) for r in rows)) if rows else len(h)
              for i, h in enumerate(headers)]
    sep = "+" + "+".join("-" * (w + 2) for w in widths) + "+"

    def fmt_row(cells):
        parts = []
        for i, cell in enumerate(cells):
            padded = cell.ljust(widths[i]) if i in left_align else cell.rjust(widths[i])
            parts.append(f" {padded} ")
        return "|" + "|".join(parts) + "|"

    lines = [sep, fmt_row(headers), sep]
    for r in rows:
        lines.append(fmt_row(r))
    lines.append(sep)
    return lines


def print_ahrd_summary_tables(summary: dict, top_n: int) -> None:
    def pct(n: int) -> str:
        return f"{100 * n / summary['n_total']:.1f}%" if summary["n_total"] else "0.0%"

    print("", file=sys.stderr)
    print("AHRD functional annotation summary:", file=sys.stderr)
    overview_rows = [
        ["Total proteins", str(summary["n_total"]), "100.0%"],
        ["With description", str(summary["n_with_desc"]), pct(summary["n_with_desc"])],
        ["Unknown / no description", str(summary["n_unknown"]), pct(summary["n_unknown"])],
        ["With GO term(s)", str(summary["n_with_go"]), pct(summary["n_with_go"])],
    ]
    for line in _ascii_table(["Metric", "Count", "Percent"], overview_rows):
        print(line, file=sys.stderr)

    print("", file=sys.stderr)
    print("AHRD-Quality-Code distribution:", file=sys.stderr)
    if summary["quality_counts"]:
        quality_rows = [[code, str(count), pct(count)]
                        for code, count in summary["quality_counts"].most_common()]
        for line in _ascii_table(["Quality Code", "Count", "Percent"], quality_rows):
            print(line, file=sys.stderr)
    else:
        print("  (no AHRD-Quality-Code column found in the output TSV)", file=sys.stderr)

    print("", file=sys.stderr)
    print(f"Top {top_n} most abundant descriptions:", file=sys.stderr)
    if summary["desc_counts"]:
        desc_rows = [[str(rank), str(count), desc] for rank, (desc, count)
                    in enumerate(summary["desc_counts"].most_common(top_n), start=1)]
        for line in _ascii_table(["Rank", "Count", "Description"], desc_rows, left_align={2}):
            print(line, file=sys.stderr)
    else:
        print("  (no annotated descriptions found)", file=sys.stderr)


def write_ahrd_summary(summary: dict, tsv_path: Path, out_path: Path, top_n: int) -> None:
    def pct(n: int) -> str:
        return f"{100 * n / summary['n_total']:.1f}%" if summary["n_total"] else "0.0%"

    lines = [
        "AHRD functional annotation summary",
        f"Source: {tsv_path}",
        "",
        f"Total proteins             : {summary['n_total']}",
        f"With description           : {summary['n_with_desc']} ({pct(summary['n_with_desc'])})",
        f"Unknown / no description   : {summary['n_unknown']} ({pct(summary['n_unknown'])})",
        f"With GO term(s)            : {summary['n_with_go']} ({pct(summary['n_with_go'])})",
        "",
        "AHRD-Quality-Code distribution:",
    ]
    if summary["quality_counts"]:
        for code, count in summary["quality_counts"].most_common():
            lines.append(f"  {code:<6} {count:>8}  ({pct(count)})")
    else:
        lines.append("  (no AHRD-Quality-Code column found in the output TSV)")

    lines += ["", f"Top {top_n} most abundant descriptions:"]
    if summary["desc_counts"]:
        for rank, (desc, count) in enumerate(summary["desc_counts"].most_common(top_n), start=1):
            lines.append(f"  {rank:>2}. {count:>6}  {desc}")
    else:
        lines.append("  (no annotated descriptions found)")

    out_path.write_text("\n".join(lines) + "\n")


def print_default_te_goterms() -> None:
    rows = [[go_id, category, name] for go_id, name, category in DEFAULT_TE_GOTERMS]
    for line in _ascii_table(["GO ID", "Category", "Name"], rows, left_align={0, 1, 2}):
        print(line, file=sys.stderr)


def load_te_goterms(path: Path) -> set:
    """Return the set of TE-associated GO IDs: DEFAULT_TE_GOTERMS unless
    --te_goterms_file is given, in which case each line's first
    GO:####### token is used (blank lines and #-comments skipped)."""
    if path is None:
        return {go_id for go_id, _, _ in DEFAULT_TE_GOTERMS}
    ids = set()
    for line in path.read_text().splitlines():
        line = line.strip()
        if not line or line.startswith("#"):
            continue
        match = re.match(r"GO:\d{7}", line)
        if match:
            ids.add(match.group(0))
    return ids


def load_detenga_status(csv_path: Path) -> dict:
    """Return {Transcript_ID: DeTEnGA_status} from DETENGA's
    <basename>_TE_summary.csv (semicolon-delimited). DeTEnGA_status
    encodes two independent TE calls as "P<interpro>M<tesort>", e.g.
    PcpM0 (Interpro: coding_sequence, TEsort: no hit), PteM0 (Interpro-only
    TE evidence), P0Mte (TEsort-only TE evidence), PteMte (both agree)."""
    status_by_id = {}
    with open(csv_path) as fh:
        header = None
        for line in fh:
            line = line.rstrip("\n")
            if not line:
                continue
            cols = line.split(";")
            if header is None:
                header = [c.strip().lower() for c in cols]
                continue
            row = dict(zip(header, cols))
            transcript_id = row.get("transcript_id", "").strip()
            if transcript_id:
                status_by_id[transcript_id] = row.get("detenga_status", "").strip()
    return status_by_id


def build_te_goterm_comparison(tsv_path: Path, te_goterms: set, detenga_status: dict) -> list:
    """Cross-check each AHRD-annotated protein's GO terms against the
    TE-associated GO set, and against DETENGA's own call. Returns a list of
    dicts (ProteinID, AHRD_GO_TEs, DETENGA_TE, GOTE_TAGGED, DETENGA_TAGGED)
    sorted by ProteinID. GOTE_TAGGED is YES only if EVERY GO term the
    protein has is TE-associated (a mix of TE and non-TE terms is NO, not
    YES); NA means the protein has no GO terms at all. DETENGA_TAGGED is
    YES/NO from DeTEnGA_status containing "te"; NA means the protein isn't
    in detenga_status at all."""
    rows = []
    for protein_id, _desc, go_ids, _quality in _iter_ahrd_rows(tsv_path):
        te_hits = [g for g in go_ids if g in te_goterms]
        if not go_ids:
            ahrd_go_tes, gote_tagged = "NA", "NA"
        elif len(te_hits) == len(go_ids):
            ahrd_go_tes, gote_tagged = ",".join(te_hits), "YES"
        else:
            ahrd_go_tes, gote_tagged = (",".join(te_hits) if te_hits else "-"), "NO"

        status = detenga_status.get(protein_id)
        if status is None:
            detenga_te, detenga_tagged = "NA", "NA"
        elif "te" in status.lower():
            detenga_te, detenga_tagged = status, "YES"
        else:
            detenga_te, detenga_tagged = "-", "NO"

        rows.append({
            "ProteinID": protein_id, "AHRD_GO_TEs": ahrd_go_tes,
            "DETENGA_TE": detenga_te, "GOTE_TAGGED": gote_tagged,
            "DETENGA_TAGGED": detenga_tagged,
        })
    rows.sort(key=lambda r: r["ProteinID"])
    return rows


def write_te_goterm_comparison(rows: list, out_path: Path) -> None:
    columns = ["ProteinID", "AHRD_GO_TEs", "DETENGA_TE", "GOTE_TAGGED", "DETENGA_TAGGED"]
    with open(out_path, "w") as fh:
        fh.write("\t".join(columns) + "\n")
        for row in rows:
            fh.write("\t".join(row[c] for c in columns) + "\n")


def print_te_goterm_agreement(rows: list) -> None:
    def is_na(row, key):
        return row[key] == "NA"

    both_te = sum(1 for r in rows if r["GOTE_TAGGED"] == "YES" and r["DETENGA_TAGGED"] == "YES")
    both_nonte = sum(1 for r in rows if r["GOTE_TAGGED"] == "NO" and r["DETENGA_TAGGED"] == "NO")
    go_only = sum(1 for r in rows if r["GOTE_TAGGED"] == "YES" and r["DETENGA_TAGGED"] == "NO")
    detenga_only = sum(1 for r in rows if r["GOTE_TAGGED"] == "NO" and r["DETENGA_TAGGED"] == "YES")
    incomplete = sum(1 for r in rows if is_na(r, "GOTE_TAGGED") or is_na(r, "DETENGA_TAGGED"))

    print("", file=sys.stderr)
    print("GOTE vs DETENGA agreement:", file=sys.stderr)
    agreement_rows = [
        ["Both TE (agree)", str(both_te)],
        ["Both non-TE (agree)", str(both_nonte)],
        ["GO-only TE (DETENGA: non-TE)", str(go_only)],
        ["DETENGA TE (GO: not TE-only)", str(detenga_only)],
        ["Incomplete info (either NA)", str(incomplete)],
    ]
    for line in _ascii_table(["Category", "Count"], agreement_rows):
        print(line, file=sys.stderr)
    print("Note: both methods can flag domesticated-TE-derived genes as "
          "false positives (real cellular function, but retained TE-like "
          "domains/GO terms) -- treat agreement/disagreement as a "
          "prioritization signal for manual review, not a final call.",
          file=sys.stderr)


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GAQET2AHRD",
        description="Build an AHRD config from a GAQET run's GAQET.log.txt "
                     "(TREMBL/SWISSPROT diamond commands) and, by default, "
                     "run AHRD on it.",
    )
    ap.add_argument("--gaqet_log", type=Path, default=None,
                    help="Path to a GAQET run's GAQET.log.txt (required "
                         "unless --print_te_associated_default_goterms)")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output directory for the AHRD config and results "
                         "(default: AHRD_run/ next to --gaqet_log)")
    ap.add_argument("--ahrd_jar", type=Path, default=None,
                    help="Path to ahrd.jar (default: $AHRD_JAR env var). "
                         "Required unless --skip_ahrd is set")
    ap.add_argument("--java_bin", default="java",
                    help="java executable to use (default: java)")
    ap.add_argument("--java_xmx", default=None,
                    help="JVM max heap size passed as -Xmx, e.g. 8g "
                         "(default: unset, use java's own default)")
    ap.add_argument("--swissprot_weight", type=int,
                    default=DEFAULT_SWISSPROT_WEIGHT,
                    help=f"AHRD blast_dbs.swissprot.weight "
                         f"(default: {DEFAULT_SWISSPROT_WEIGHT})")
    ap.add_argument("--trembl_weight", type=int,
                    default=DEFAULT_TREMBL_WEIGHT,
                    help=f"AHRD blast_dbs.trembl.weight "
                         f"(default: {DEFAULT_TREMBL_WEIGHT})")
    ap.add_argument("--swissprot_desc_weight", type=float,
                    default=DEFAULT_SWISSPROT_DESC_SCORE_BIT_SCORE_WEIGHT,
                    help="AHRD blast_dbs.swissprot.description_score_bit_"
                         f"score_weight (default: "
                         f"{DEFAULT_SWISSPROT_DESC_SCORE_BIT_SCORE_WEIGHT})")
    ap.add_argument("--trembl_desc_weight", type=float,
                    default=DEFAULT_TREMBL_DESC_SCORE_BIT_SCORE_WEIGHT,
                    help="AHRD blast_dbs.trembl.description_score_bit_"
                         f"score_weight (default: "
                         f"{DEFAULT_TREMBL_DESC_SCORE_BIT_SCORE_WEIGHT})")
    ap.add_argument("--token_bit_score_weight", type=float,
                    default=DEFAULT_TOKEN_SCORE_BIT_SCORE_WEIGHT,
                    help="AHRD token_score_bit_score_weight "
                         f"(default: {DEFAULT_TOKEN_SCORE_BIT_SCORE_WEIGHT})")
    ap.add_argument("--token_database_score_weight", type=float,
                    default=DEFAULT_TOKEN_SCORE_DATABASE_SCORE_WEIGHT,
                    help="AHRD token_score_database_score_weight (default: "
                         f"{DEFAULT_TOKEN_SCORE_DATABASE_SCORE_WEIGHT})")
    ap.add_argument("--token_overlap_score_weight", type=float,
                    default=DEFAULT_TOKEN_SCORE_OVERLAP_SCORE_WEIGHT,
                    help="AHRD token_score_overlap_score_weight (default: "
                         f"{DEFAULT_TOKEN_SCORE_OVERLAP_SCORE_WEIGHT})")
    ap.add_argument("--ahrd_home", type=Path, default=None,
                    help="Path to an AHRD source/distribution checkout "
                         "(default: $AHRD_HOME env var). When set, defaults "
                         "--blacklist/--swissprot_filter/--trembl_filter/"
                         "--token_blacklist to the word-filtering files "
                         "bundled under {ahrd_home}/test/resources/")
    ap.add_argument("--blacklist", type=Path, default=None,
                    help="AHRD blacklist_descline file, applied to both "
                         "swissprot and trembl (default: "
                         "{ahrd_home}/test/resources/blacklist_descline.txt "
                         "if --ahrd_home is set, else omitted)")
    ap.add_argument("--swissprot_filter", type=Path, default=None,
                    help="AHRD filter_descline file for swissprot (default: "
                         "{ahrd_home}/test/resources/filter_descline_sprot.txt "
                         "if --ahrd_home is set, else omitted)")
    ap.add_argument("--trembl_filter", type=Path, default=None,
                    help="AHRD filter_descline file for trembl (default: "
                         "{ahrd_home}/test/resources/filter_descline_trembl.txt "
                         "if --ahrd_home is set, else omitted)")
    ap.add_argument("--token_blacklist", type=Path, default=None,
                    help="AHRD blacklist_token file, applied to both "
                         "swissprot and trembl (default: "
                         "{ahrd_home}/test/resources/blacklist_token.txt "
                         "if --ahrd_home is set, else omitted)")
    ap.add_argument("--gene_ontology_result", type=Path, default=None,
                    help="GO Annotation (GAF) file for GO term transfer, "
                         "e.g. goa_uniprot_all.gaf (default: "
                         "goa_uniprot_all.gaf in the same directory as the "
                         "SWISSPROT --db from GAQET.log.txt)")
    ap.add_argument("--skip_go", action="store_true",
                    help="Do not transfer GO terms (omit gene_ontology_"
                         "result/reference_go_regex/prefer_reference_with_"
                         "go_annos from the config)")
    ap.add_argument("--top_n", type=int, default=DEFAULT_TOP_N,
                    help=f"Number of most abundant descriptions to list in "
                         f"the summary (default: {DEFAULT_TOP_N})")
    ap.add_argument("--skip_summary", action="store_true",
                    help="Do not write the <prefix>_AHRD.summary.txt "
                         "report after AHRD finishes")
    ap.add_argument("--check_te_goterms", action="store_true",
                    help="Cross-check AHRD's GO terms against DETENGA's TE "
                         "calls; writes <prefix>_TEGOterm_vs_DETENGA.tsv")
    ap.add_argument("--te_goterms_file", type=Path, default=None,
                    help="TE-associated GO term list, one GO:####### per "
                         "line (default: a hardcoded list -- see "
                         "--print_te_associated_default_goterms)")
    ap.add_argument("--detenga_csv", type=Path, default=None,
                    help="DETENGA's <basename>_TE_summary.csv (default: "
                         "DETENGA_run/{prefix}_TE_summary.csv next to "
                         "--gaqet_log)")
    ap.add_argument("--print_te_associated_default_goterms", action="store_true",
                    help="Print the hardcoded default TE-associated GO "
                         "term list as a table, then exit")
    ap.add_argument("--skip_ahrd", action="store_true",
                    help="Write the AHRD config only; do not invoke AHRD")
    ap.add_argument("--dry_run", action="store_true",
                    help="Parse GAQET.log.txt and print what would be "
                         "written/run, then exit without writing or "
                         "running anything")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if args.print_te_associated_default_goterms:
        print_default_te_goterms()
        return

    if args.gaqet_log is None:
        print("ERROR: --gaqet_log is required (unless "
              "--print_te_associated_default_goterms)", file=sys.stderr)
        sys.exit(1)
    if not args.gaqet_log.exists():
        print(f"ERROR: --gaqet_log file not found: {args.gaqet_log}",
              file=sys.stderr)
        sys.exit(1)

    te_goterms = None
    if args.check_te_goterms:
        if args.te_goterms_file is not None and not args.te_goterms_file.exists():
            print(f"ERROR: --te_goterms_file not found: "
                  f"{args.te_goterms_file}", file=sys.stderr)
            sys.exit(1)
        te_goterms = load_te_goterms(args.te_goterms_file)

    ahrd_jar = args.ahrd_jar or (Path(os.environ["AHRD_JAR"])
                                  if os.environ.get("AHRD_JAR") else None)
    if not args.skip_ahrd and ahrd_jar is None:
        print("ERROR: --ahrd_jar not given and $AHRD_JAR is not set. "
              "Either pass --ahrd_jar, export AHRD_JAR, or use --skip_ahrd "
              "to only write the config.", file=sys.stderr)
        sys.exit(1)

    ahrd_home = args.ahrd_home or (Path(os.environ["AHRD_HOME"])
                                    if os.environ.get("AHRD_HOME") else None)
    resources = ahrd_home / "test" / "resources" if ahrd_home else None
    blacklist = args.blacklist or (resources / "blacklist_descline.txt" if resources else None)
    swissprot_filter = args.swissprot_filter or (resources / "filter_descline_sprot.txt" if resources else None)
    trembl_filter = args.trembl_filter or (resources / "filter_descline_trembl.txt" if resources else None)
    token_blacklist = args.token_blacklist or (resources / "blacklist_token.txt" if resources else None)
    for label, path in (("blacklist", blacklist), ("swissprot filter", swissprot_filter),
                        ("trembl filter", trembl_filter), ("token_blacklist", token_blacklist)):
        if path is not None and not path.exists():
            print(f"WARNING: {label} file not found: {path}", file=sys.stderr)

    parsed = parse_gaqet_log(args.gaqet_log)
    print(f"Parsed GAQET.log.txt: proteins_fasta={parsed['proteins_fasta']}",
          file=sys.stderr)

    if not parsed["proteins_fasta"].exists():
        print(f"WARNING: proteins FASTA from GAQET.log.txt not found on "
              f"disk: {parsed['proteins_fasta']}", file=sys.stderr)
    for key in ("swissprot_out", "trembl_out"):
        if not parsed[key].exists():
            print(f"WARNING: diamond output from GAQET.log.txt not found "
                  f"on disk: {parsed[key]}", file=sys.stderr)

    swissprot_fasta = derive_db_fasta(parsed["swissprot_db"])
    trembl_fasta = derive_db_fasta(parsed["trembl_db"])
    for label, fasta in (("swissprot", swissprot_fasta), ("trembl", trembl_fasta)):
        if not fasta.exists():
            print(f"WARNING: derived {label} database FASTA not found: "
                  f"{fasta} (derived from {label} .dmnd path by swapping "
                  f"the extension; pass a matching file there or fix the "
                  f"path manually in the generated config)", file=sys.stderr)

    gene_ontology_result = None
    if not args.skip_go:
        gene_ontology_result = (args.gene_ontology_result
                                 or parsed["swissprot_db"].parent / "goa_uniprot_all.gaf")
        if not gene_ontology_result.exists():
            print(f"WARNING: gene_ontology_result file not found: "
                  f"{gene_ontology_result} (defaulted to goa_uniprot_all.gaf "
                  f"next to the SWISSPROT --db; pass --gene_ontology_result "
                  f"to point at the right file, or --skip_go to omit GO "
                  f"term transfer)", file=sys.stderr)

    prefix = find_gaqet_prefix(args.gaqet_log)
    output_dir = args.output or (args.gaqet_log.parent / "AHRD_run")
    config_path = output_dir / f"{prefix}_AHRD_config.yml"
    output_tsv = output_dir / f"{prefix}_AHRD.tsv"

    if args.dry_run:
        print(f"  Config      : {config_path}", file=sys.stderr)
        print(f"  AHRD output : {output_tsv}", file=sys.stderr)
        print(f"  SwissProt   : file={parsed['swissprot_out']} "
              f"database={swissprot_fasta}", file=sys.stderr)
        print(f"  TrEMBL      : file={parsed['trembl_out']} "
              f"database={trembl_fasta}", file=sys.stderr)
        print(f"  GO terms    : {gene_ontology_result if gene_ontology_result else 'skipped (--skip_go)'}",
              file=sys.stderr)
        if args.skip_ahrd:
            print("  AHRD run    : skipped (--skip_ahrd)", file=sys.stderr)
        else:
            print(f"  AHRD run    : {args.java_bin} "
                  f"{'-Xmx' + args.java_xmx + ' ' if args.java_xmx else ''}"
                  f"-jar {ahrd_jar} {config_path}", file=sys.stderr)
        if args.check_te_goterms:
            detenga_csv = args.detenga_csv or (
                args.gaqet_log.parent / "DETENGA_run" / f"{prefix}_TE_summary.csv")
            print(f"  TE GO check : {len(te_goterms)} TE GO term(s), "
                  f"DETENGA_csv={detenga_csv}"
                  f"{' (not found)' if not detenga_csv.exists() else ''}",
                  file=sys.stderr)
        print("  Exiting (--dry_run). Nothing written or run.",
              file=sys.stderr)
        return

    output_dir.mkdir(parents=True, exist_ok=True)
    config_text = build_ahrd_config(
        proteins_fasta=parsed["proteins_fasta"],
        output_tsv=output_tsv,
        swissprot_out=parsed["swissprot_out"], swissprot_fasta=swissprot_fasta,
        trembl_out=parsed["trembl_out"], trembl_fasta=trembl_fasta,
        swissprot_weight=args.swissprot_weight,
        trembl_weight=args.trembl_weight,
        swissprot_desc_weight=args.swissprot_desc_weight,
        trembl_desc_weight=args.trembl_desc_weight,
        token_bit_score_weight=args.token_bit_score_weight,
        token_database_score_weight=args.token_database_score_weight,
        token_overlap_score_weight=args.token_overlap_score_weight,
        blacklist=blacklist, swissprot_filter=swissprot_filter,
        trembl_filter=trembl_filter, token_blacklist=token_blacklist,
        gene_ontology_result=gene_ontology_result,
    )
    config_path.write_text(config_text)
    print(f"Written: {config_path}", file=sys.stderr)

    if args.skip_ahrd:
        print("AHRD run skipped (--skip_ahrd)", file=sys.stderr)
    else:
        cmd = [args.java_bin]
        if args.java_xmx:
            cmd.append(f"-Xmx{args.java_xmx}")
        cmd += ["-jar", str(ahrd_jar), str(config_path)]
        print(f"  $ {' '.join(cmd)}", file=sys.stderr)
        result = subprocess.run(cmd)
        if result.returncode != 0:
            print(f"ERROR: AHRD exited with code {result.returncode}",
                  file=sys.stderr)
            sys.exit(1)
        print(f"AHRD finished. Output: {output_tsv}", file=sys.stderr)

    # Summary and --check_te_goterms run against output_tsv whenever it
    # exists -- either just produced above, or from an earlier run when
    # --skip_ahrd is set (e.g. to only regenerate the summary/TE check).
    if not output_tsv.exists():
        if not args.skip_summary or args.check_te_goterms:
            print(f"WARNING: AHRD output TSV not found: {output_tsv} — "
                  f"skipping summary/TE check", file=sys.stderr)
        return

    if args.skip_summary:
        print("Summary skipped (--skip_summary)", file=sys.stderr)
    else:
        summary = summarize_ahrd_output(output_tsv)
        summary_path = output_tsv.with_suffix(".summary.txt")
        write_ahrd_summary(summary, output_tsv, summary_path, args.top_n)
        print(f"Written: {summary_path}", file=sys.stderr)
        print_ahrd_summary_tables(summary, args.top_n)

    if args.check_te_goterms:
        detenga_csv = args.detenga_csv or (
            args.gaqet_log.parent / "DETENGA_run" / f"{prefix}_TE_summary.csv")
        if detenga_csv.exists():
            detenga_status = load_detenga_status(detenga_csv)
        else:
            print(f"WARNING: --detenga_csv not found: {detenga_csv} — "
                  f"DETENGA_TE/DETENGA_TAGGED will be NA for every protein",
                  file=sys.stderr)
            detenga_status = {}

        te_rows = build_te_goterm_comparison(output_tsv, te_goterms, detenga_status)
        te_out_path = output_tsv.with_name(f"{prefix}_TEGOterm_vs_DETENGA.tsv")
        write_te_goterm_comparison(te_rows, te_out_path)
        print(f"Written: {te_out_path}", file=sys.stderr)
        print_te_goterm_agreement(te_rows)


if __name__ == "__main__":
    main()
