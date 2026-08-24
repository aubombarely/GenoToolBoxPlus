#!/usr/bin/env python3
"""
GAQET2AHRD.py — Build an AHRD config from a GAQET run and (optionally) run
AHRD (https://github.com/groupschoof/AHRD) on it.

GAQET (https://github.com/aubombarely/GAQET) writes a GAQET.log.txt inside
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
from pathlib import Path

VERSION = "v0.0.1"

DEFAULT_SWISSPROT_WEIGHT = 653
DEFAULT_TREMBL_WEIGHT = 904
DEFAULT_SWISSPROT_DESC_SCORE_BIT_SCORE_WEIGHT = 2.717061
DEFAULT_TREMBL_DESC_SCORE_BIT_SCORE_WEIGHT = 2.590211
DEFAULT_TOKEN_SCORE_BIT_SCORE_WEIGHT = 0.468
DEFAULT_TOKEN_SCORE_DATABASE_SCORE_WEIGHT = 0.2098
DEFAULT_TOKEN_SCORE_OVERLAP_SCORE_WEIGHT = 0.3221


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
    token_blacklist: Path,
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

    return (
        f"proteins_fasta: {proteins_fasta}\n"
        f"token_score_bit_score_weight: {token_bit_score_weight}\n"
        f"token_score_database_score_weight: {token_database_score_weight}\n"
        f"token_score_overlap_score_weight: {token_overlap_score_weight}\n"
        f"blast_dbs:\n"
        f"{db_block('swissprot', swissprot_weight, swissprot_desc_weight, swissprot_out, swissprot_fasta, swissprot_filter)}\n"
        f"{db_block('trembl', trembl_weight, trembl_desc_weight, trembl_out, trembl_fasta, trembl_filter)}\n"
        f"output: {output_tsv}\n"
    )


def find_gaqet_prefix(gaqet_log: Path) -> str:
    """GAQET names its per-run stats file {prefix}_GAQET.stats.tsv in the
    same directory as GAQET.log.txt; used only to name AHRD's own output
    files. Falls back to 'AHRD' if no such file is found."""
    hits = glob.glob(str(gaqet_log.parent / "*_GAQET.stats.tsv"))
    if not hits:
        return "AHRD"
    return Path(hits[0]).name[: -len("_GAQET.stats.tsv")]


def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GAQET2AHRD",
        description="Build an AHRD config from a GAQET run's GAQET.log.txt "
                     "(TREMBL/SWISSPROT diamond commands) and, by default, "
                     "run AHRD on it.",
    )
    ap.add_argument("--gaqet_log", required=True, type=Path,
                    help="Path to a GAQET run's GAQET.log.txt")
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
    ap.add_argument("--skip_ahrd", action="store_true",
                    help="Write the AHRD config only; do not invoke AHRD")
    ap.add_argument("--dry_run", action="store_true",
                    help="Parse GAQET.log.txt and print what would be "
                         "written/run, then exit without writing or "
                         "running anything")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.gaqet_log.exists():
        print(f"ERROR: --gaqet_log file not found: {args.gaqet_log}",
              file=sys.stderr)
        sys.exit(1)

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
        if args.skip_ahrd:
            print("  AHRD run    : skipped (--skip_ahrd)", file=sys.stderr)
        else:
            print(f"  AHRD run    : {args.java_bin} "
                  f"{'-Xmx' + args.java_xmx + ' ' if args.java_xmx else ''}"
                  f"-jar {ahrd_jar} {config_path}", file=sys.stderr)
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
    )
    config_path.write_text(config_text)
    print(f"Written: {config_path}", file=sys.stderr)

    if args.skip_ahrd:
        print("AHRD run skipped (--skip_ahrd)", file=sys.stderr)
        return

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


if __name__ == "__main__":
    main()
