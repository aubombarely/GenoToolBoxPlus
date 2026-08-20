#!/usr/bin/env python3
"""
GFF3RenameGenes.py — Systematically rename gene models in a GFF3 file.

Renumbers every gene, transcript, exon, CDS, and UTR feature using a
fixed, predictable scheme built from the sequence ID:

    Gene:       {SeqID}G{gene_number}              e.g. PhangC01G000010
    Transcript: {gene_ID}T{transcript_number}       e.g. PhangC01G000010T01
    Exon:       {transcript_ID}EXO{exon_number}     e.g. PhangC01G000010T01EXO01
    CDS:        {transcript_ID}CDS{cds_number}      e.g. PhangC01G000010T01CDS01
    UTR:        {transcript_ID}UTR{utr_number}      e.g. PhangC01G000010T01UTR01
                (five_prime_UTR and three_prime_UTR share one numbering,
                5' -> 3' along the transcript)

Gene numbers restart at 10 on every new SeqID and increment by 10
(10, 20, 30 ...), leaving gaps for later manual insertions. Transcript
numbers restart at 01 for every gene; exon/CDS/UTR numbers each restart
at 01 for every transcript (independent counters per feature type) and
follow transcript direction (ascending start on '+', descending on '-').

--prefix_geneid replaces the literal SeqID with one fixed prefix across
the whole file (e.g. all chromosomes get "Phan" instead of their own
PhangC01/PhangC02/...); since that prefix is then no longer unique per
SeqID, gene numbering switches to one continuous count across the whole
file instead of restarting per SeqID. --after_seqid_tag inserts an extra
tag between the SeqID (or --prefix_geneid) and the 'G', e.g.
PhangC01ANN2G000010 with --after_seqid_tag ANN2 -- useful when the same
assembly has more than one annotation and gene IDs need to stay
distinguishable between them; it does not affect numbering.

Parent hierarchy is rebuilt from the *new* IDs. The feature's original
ID (or, if it had none, "{seqid}:{start}-{end}") is preserved as a
trailing OldFeatID= attribute; every other original attribute (Note=,
product=, Dbxref=, ...) is kept as-is.

A "gene" is any column-3 type matching --gene_type (default: gene).
A transcript is any feature whose Parent points to a gene (regardless
of its own column-3 type: mRNA, tRNA, ncRNA, ...). A subfeature is any
feature whose Parent points to a transcript. Everything else (pragma
lines, comments, and features outside this hierarchy, e.g. standalone
region/chromosome lines) passes through unchanged.

Before renaming, a structural sanity check runs over the whole file:
duplicate IDs, dangling Parent references, missing IDs on genes/
transcripts, invalid coordinates, seqid mismatches between a feature and
its Parent (all error-level, abort by default), plus genes/transcripts
with no children, out-of-range child coordinates, strand mismatches, and
likely exact-duplicate features from merged annotation runs (all
warning-level, reported but non-fatal). This is deliberately generic
rather than tool-specific, since output from different annotation tools
(BRAKER, MAKER, HELIXER, EVIANN, EGAPX, ANNEVO, TransDecoder, ...) fails
in different ways. Use --check_only to just run the check and exit,
--force to rename anyway despite errors, or --skip_sanity_check to
bypass the check entirely.

Usage
-----
    GFF3RenameGenes.py --gff annotation.gff3
    GFF3RenameGenes.py --gff annotation.gff3 --output renamed.gff3
    GFF3RenameGenes.py --gff annotation.gff3 --check_only
    GFF3RenameGenes.py --gff annotation.gff3 --dry_run
"""

import argparse
import re
import sys
from pathlib import Path

VERSION = "v0.0.1"

DEFAULT_GENE_PAD       = 6
DEFAULT_GENE_STEP      = 10
DEFAULT_GENE_START     = 10
DEFAULT_TRANSCRIPT_PAD = 2
DEFAULT_FEATURE_PAD    = 2

UTR_TYPES = {"five_prime_utr", "three_prime_utr", "5'utr", "3'utr", "utr"}
SUFFIX_MAP = {"exon": "EXO", "cds": "CDS"}  # UTR_TYPES handled separately
# Both matched case-insensitively (lower-cased ftype) -- annotation tools
# vary in casing (BRAKER/MAKER/HELIXER/EVIANN/EGAPX/ANNEVO/TransDecoder),
# and the fallback (first 3 letters of the type, upper-cased) already
# self-corrects for exon/CDS casing variants by coincidence, but not for
# five_prime_UTR/three_prime_UTR, so casing is normalized explicitly here.


# ── GFF3 parsing ──────────────────────────────────────────────────────────────

def parse_attributes(attr_str: str) -> dict:
    attrs = {}
    for field in attr_str.strip().rstrip(";").split(";"):
        field = field.strip()
        if not field or "=" not in field:
            continue
        key, val = field.split("=", 1)
        attrs[key.strip()] = val.strip()
    return attrs


def format_attributes(attrs: dict) -> str:
    return ";".join(f"{k}={v}" for k, v in attrs.items())


def natural_key(s: str) -> list:
    return [int(t) if t.isdigit() else t.lower()
            for t in re.split(r"(\d+)", s)]


def read_gff(path: Path) -> tuple:
    """Return (pragma_lines, features). features is a list of dicts:
    seqid, source, ftype, start, end, score, strand, phase, attrs,
    old_id (str or None), old_parents (list of str), lineno."""
    pragma_lines = []
    features = []
    n_skipped = 0
    with open(path) as fh:
        for lineno, line in enumerate(fh, 1):
            raw = line.rstrip("\n")
            if not raw.strip():
                continue
            if raw.startswith("#"):
                pragma_lines.append(raw)
                continue
            cols = raw.split("\t")
            if len(cols) < 9:
                print(f"WARNING: skipping malformed GFF3 line {lineno} "
                      f"(fewer than 9 columns)", file=sys.stderr)
                n_skipped += 1
                continue
            try:
                start, end = int(cols[3]), int(cols[4])
            except ValueError:
                print(f"WARNING: skipping GFF3 line {lineno} "
                      f"(non-integer start/end)", file=sys.stderr)
                n_skipped += 1
                continue
            attrs = parse_attributes(cols[8])
            old_id = attrs.get("ID")
            old_parents = [p.strip() for p in attrs.get("Parent", "").split(",")
                           if p.strip()]
            features.append({
                "seqid": cols[0], "source": cols[1], "ftype": cols[2],
                "start": start, "end": end, "score": cols[5],
                "strand": cols[6], "phase": cols[7], "attrs": attrs,
                "old_id": old_id, "old_parents": old_parents,
                "lineno": lineno,
            })
    if n_skipped:
        print(f"WARNING: {n_skipped} malformed line(s) skipped entirely",
              file=sys.stderr)
    return pragma_lines, features


# ── Classification (shared by sanity-check and renaming) ──────────────────────

def classify_features(features: list, gene_type: str) -> dict:
    """
    Split features into the gene/transcript/subfeature/passthrough
    hierarchy, purely from ID/Parent linkage (independent of the actual
    column-3 type, aside from --gene_type itself). Returns a dict with
    genes, transcripts, subfeatures, passthrough (lists) and by_id (dict
    old_id -> feature, for features that had an ID).
    """
    genes = [f for f in features if f["ftype"] == gene_type]
    gene_ids = {f["old_id"] for f in genes if f["old_id"]}

    transcripts = [f for f in features
                   if f["ftype"] != gene_type
                   and any(p in gene_ids for p in f["old_parents"])]
    transcript_ids = {f["old_id"] for f in transcripts if f["old_id"]}
    transcript_obj_ids = {id(f) for f in transcripts}

    subfeatures = [f for f in features
                   if f["old_id"] not in gene_ids
                   and id(f) not in transcript_obj_ids
                   and any(p in transcript_ids for p in f["old_parents"])]

    handled = set(id(f) for f in genes) | transcript_obj_ids \
        | set(id(f) for f in subfeatures)
    passthrough = [f for f in features if id(f) not in handled]

    by_id = {f["old_id"]: f for f in features if f["old_id"]}

    return {
        "genes": genes, "gene_ids": gene_ids,
        "transcripts": transcripts, "transcript_ids": transcript_ids,
        "subfeatures": subfeatures, "passthrough": passthrough,
        "by_id": by_id,
    }


# ── Sanity check ──────────────────────────────────────────────────────────────

def validate_gff(features: list, classified: dict) -> tuple:
    """
    Structural/referential integrity checks, run before renaming so a
    malformed input (duplicate IDs, dangling Parent references, bad
    coordinates ...) is caught with a clear diagnostic instead of
    silently producing a corrupted rename. Different annotation tools
    (BRAKER, MAKER, HELIXER, EVIANN, EGAPX, ANNEVO, TransDecoder ...)
    have different failure modes, so this is deliberately generic rather
    than tool-specific. Returns (errors, warnings) -- errors indicate the
    renaming would be unreliable; warnings are informational.
    """
    errors: list = []
    warnings: list = []

    genes        = classified["genes"]
    gene_ids     = classified["gene_ids"]
    transcripts  = classified["transcripts"]
    transcript_ids = classified["transcript_ids"]
    subfeatures  = classified["subfeatures"]
    by_id        = classified["by_id"]

    # 1. Duplicate IDs anywhere in the file
    id_seen: dict = {}
    for f in features:
        if f["old_id"]:
            id_seen.setdefault(f["old_id"], []).append(f)
    dup_examples = []
    n_dup = 0
    for old_id, flist in id_seen.items():
        if len(flist) > 1:
            n_dup += 1
            if len(dup_examples) < 5:
                lines = ",".join(str(x["lineno"]) for x in flist)
                dup_examples.append(f"'{old_id}' on lines {lines}")
    if n_dup:
        errors.append(f"{n_dup} ID(s) used on more than one line, e.g.: "
                      + "; ".join(dup_examples))

    # 2. Genes/transcripts missing an ID (subfeatures may lack one --
    #    handled via the OldFeatID seqid:start-end fallback)
    n_gene_no_id = sum(1 for g in genes if not g["old_id"])
    if n_gene_no_id:
        errors.append(f"{n_gene_no_id} '{genes[0]['ftype'] if genes else '?'}' "
                      f"feature(s) have no ID attribute")
    n_transcript_no_id = sum(1 for t in transcripts if not t["old_id"])
    if n_transcript_no_id:
        errors.append(f"{n_transcript_no_id} transcript-level feature(s) "
                      f"(Parent points to a gene) have no ID attribute")

    # 3. Dangling Parent references
    n_dangling = 0
    dangling_examples = []
    for f in features:
        for p in f["old_parents"]:
            if p not in id_seen:
                n_dangling += 1
                if len(dangling_examples) < 5:
                    dangling_examples.append(f"line {f['lineno']} Parent={p}")
    if n_dangling:
        errors.append(f"{n_dangling} Parent reference(s) point to an ID that "
                      f"doesn't exist anywhere in the file, e.g.: "
                      + "; ".join(dangling_examples))

    # 4. Invalid coordinates
    n_bad_coord = 0
    for f in features:
        if f["start"] < 1 or f["start"] > f["end"]:
            n_bad_coord += 1
    if n_bad_coord:
        errors.append(f"{n_bad_coord} feature(s) have invalid coordinates "
                      f"(start < 1 or start > end)")

    # 5. seqid mismatch between a feature and its Parent
    n_seqid_mismatch = 0
    for f in features:
        for p in f["old_parents"]:
            parent = by_id.get(p)
            if parent is not None and parent["seqid"] != f["seqid"]:
                n_seqid_mismatch += 1
    if n_seqid_mismatch:
        errors.append(f"{n_seqid_mismatch} feature(s) are on a different "
                      f"SeqID than their Parent (corrupted or merged GFF3)")

    # ── Warnings (non-fatal) ────────────────────────────────────────────────

    # 6. Genes with no transcript children
    genes_with_children = {p for t in transcripts for p in t["old_parents"]}
    n_gene_no_transcript = sum(1 for g in genes
                               if g["old_id"] and g["old_id"] not in genes_with_children)
    if n_gene_no_transcript:
        warnings.append(f"{n_gene_no_transcript} gene(s) have no transcript "
                        f"children")

    # 7. Transcripts with no exon/CDS/UTR children
    transcripts_with_children = {p for s in subfeatures for p in s["old_parents"]}
    n_transcript_no_children = sum(1 for t in transcripts
                                   if t["old_id"] and t["old_id"] not in transcripts_with_children)
    if n_transcript_no_children:
        warnings.append(f"{n_transcript_no_children} transcript(s) have no "
                        f"exon/CDS/UTR children")

    # 8. Child coordinate range not contained within its Parent's range
    n_out_of_bounds = 0
    for f in transcripts + subfeatures:
        for p in f["old_parents"]:
            parent = by_id.get(p)
            if parent is not None and (f["start"] < parent["start"] or f["end"] > parent["end"]):
                n_out_of_bounds += 1
                break
    if n_out_of_bounds:
        warnings.append(f"{n_out_of_bounds} feature(s) extend beyond their "
                        f"Parent's coordinate range")

    # 9. Strand mismatch between a feature and its Parent
    n_strand_mismatch = 0
    for f in transcripts + subfeatures:
        for p in f["old_parents"]:
            parent = by_id.get(p)
            if (parent is not None and parent["strand"] in ("+", "-")
                    and f["strand"] in ("+", "-") and parent["strand"] != f["strand"]):
                n_strand_mismatch += 1
                break
    if n_strand_mismatch:
        warnings.append(f"{n_strand_mismatch} feature(s) have a different "
                        f"strand than their Parent")

    # 10. Exact-duplicate features (possible artifact of merging multiple
    #     annotation runs, e.g. combining several BRAKER predictions)
    seen_exact: dict = {}
    for f in features:
        key = (f["seqid"], f["start"], f["end"], f["strand"], f["ftype"],
              tuple(sorted(f["old_parents"])))
        seen_exact[key] = seen_exact.get(key, 0) + 1
    n_exact_dup = sum(c - 1 for c in seen_exact.values() if c > 1)
    if n_exact_dup:
        warnings.append(f"{n_exact_dup} feature(s) appear to be exact "
                        f"duplicates (same seqid/coordinates/strand/type/"
                        f"parent) -- possibly from merging multiple "
                        f"annotation runs")

    return errors, warnings


# ── Renaming ──────────────────────────────────────────────────────────────────

def rename_features(classified: dict, gene_type: str,
                     gene_pad: int, gene_step: int, gene_start: int,
                     transcript_pad: int, feature_pad: int,
                     prefix_geneid: str | None = None,
                     after_seqid_tag: str | None = None) -> dict:
    """
    Compute new IDs and rebuild Parent links in place (feature['attrs']).
    Returns a stats dict for the run summary.

    prefix_geneid, if given, replaces the literal SeqID in the gene ID
    with a single fixed prefix for the whole file -- since that prefix is
    no longer unique per SeqID, gene numbering switches from restarting
    at gene_start on every new SeqID to one continuous count across the
    whole file (still SeqID-ordered, then by start), so IDs stay unique.
    after_seqid_tag, if given, is inserted between the SeqID (or
    prefix_geneid) and the 'G' -- e.g. for tagging genes from a specific
    annotation round when the same assembly has multiple annotations.
    Does not affect numbering (still restarts per real SeqID unless
    prefix_geneid is also set).
    """
    genes       = classified["genes"]
    transcripts = classified["transcripts"]
    subfeatures = classified["subfeatures"]
    passthrough = classified["passthrough"]

    id_map: dict = {}
    tag = after_seqid_tag or ""

    # ── Genes: numbered per SeqID (or continuously, if prefix_geneid is
    #    set), ordered by start ──────────────────────────────────────────
    genes_by_seq: dict = {}
    for g in genes:
        genes_by_seq.setdefault(g["seqid"], []).append(g)

    n_genes_no_id = 0
    n = gene_start  # only used continuously when prefix_geneid is set
    for seqid in sorted(genes_by_seq, key=natural_key):
        seq_genes = sorted(genes_by_seq[seqid],
                           key=lambda g: (g["start"], g["old_id"] or ""))
        if prefix_geneid is None:
            n = gene_start  # restart per SeqID
        base = prefix_geneid if prefix_geneid is not None else seqid
        for g in seq_genes:
            new_id = f"{base}{tag}G{n:0{gene_pad}d}"
            g["new_id"] = new_id
            if g["old_id"]:
                id_map[g["old_id"]] = new_id
            else:
                n_genes_no_id += 1
            n += gene_step

    # ── Transcripts: number per gene, ordered by start ─────────────────────
    transcripts_by_gene: dict = {}
    for t in transcripts:
        primary_parent = t["old_parents"][0] if t["old_parents"] else None
        transcripts_by_gene.setdefault(primary_parent, []).append(t)

    n_transcripts_no_id = 0
    for g in genes:
        t_list = sorted(transcripts_by_gene.get(g["old_id"], []),
                        key=lambda t: (t["start"], t["old_id"] or ""))
        for i, t in enumerate(t_list, start=1):
            new_id = f"{g['new_id']}T{i:0{transcript_pad}d}"
            t["new_id"] = new_id
            if t["old_id"]:
                id_map[t["old_id"]] = new_id
            else:
                n_transcripts_no_id += 1

    # ── Subfeatures: number per (transcript, suffix), strand-aware order ───
    sub_by_transcript: dict = {}
    for s in subfeatures:
        primary_parent = s["old_parents"][0] if s["old_parents"] else None
        sub_by_transcript.setdefault(primary_parent, []).append(s)

    n_other_type = 0
    n_sub_no_id  = 0
    for t in transcripts:
        s_list = sub_by_transcript.get(t["old_id"], [])
        groups: dict = {}
        for s in s_list:
            ftype_lc = s["ftype"].lower()
            if ftype_lc in UTR_TYPES:
                suffix = "UTR"
            elif ftype_lc in SUFFIX_MAP:
                suffix = SUFFIX_MAP[ftype_lc]
            else:
                suffix = s["ftype"].upper()[:3]
                n_other_type += 1
            groups.setdefault(suffix, []).append(s)

        reverse = (t["strand"] == "-")
        for suffix, group in groups.items():
            ordered = sorted(group, key=lambda s: s["start"], reverse=reverse)
            for i, s in enumerate(ordered, start=1):
                new_id = f"{t['new_id']}{suffix}{i:0{feature_pad}d}"
                s["new_id"] = new_id
                if s["old_id"]:
                    id_map[s["old_id"]] = new_id
                else:
                    n_sub_no_id += 1

    # ── Rewrite attributes ──────────────────────────────────────────────────
    n_unresolved_parent = 0
    for group, is_top in ((genes, True), (transcripts, False), (subfeatures, False)):
        for f in group:
            old_attrs = f["attrs"]
            new_attrs = {}
            new_attrs["ID"] = f["new_id"]
            if not is_top:
                new_parents = []
                for p in f["old_parents"]:
                    if p in id_map:
                        new_parents.append(id_map[p])
                    else:
                        new_parents.append(p)
                        n_unresolved_parent += 1
                new_attrs["Parent"] = ",".join(new_parents)
            for k, v in old_attrs.items():
                if k in ("ID", "Parent"):
                    continue
                new_attrs[k] = v
            old_feat_id = f["old_id"] or f"{f['seqid']}:{f['start']}-{f['end']}"
            new_attrs["OldFeatID"] = old_feat_id
            f["attrs"] = new_attrs

    stats = {
        "n_genes": len(genes),
        "n_transcripts": len(transcripts),
        "n_subfeatures": len(subfeatures),
        "n_passthrough": len(passthrough),
        "n_genes_no_id": n_genes_no_id,
        "n_transcripts_no_id": n_transcripts_no_id,
        "n_sub_no_id": n_sub_no_id,
        "n_other_type": n_other_type,
        "n_unresolved_parent": n_unresolved_parent,
        "n_seqids": len(genes_by_seq),
    }
    return stats


# ── Output ────────────────────────────────────────────────────────────────────

def write_gff(pragma_lines: list, features: list, gene_type: str, out_fh) -> None:
    out_fh.write("##gff-version 3\n")
    for line in pragma_lines:
        if line.strip() != "##gff-version 3":
            out_fh.write(line + "\n")

    genes = [f for f in features if f["ftype"] == gene_type]

    # Build parent -> children index using ORIGINAL old_id linkage
    # (already validated during renaming).
    children_of_old_id: dict = {}
    for f in features:
        for p in f["old_parents"]:
            children_of_old_id.setdefault(p, []).append(f)

    passthrough_seen = set()

    def write_feature(f):
        cols = [f["seqid"], f["source"], f["ftype"], str(f["start"]),
                str(f["end"]), f["score"], f["strand"], f["phase"],
                format_attributes(f["attrs"])]
        out_fh.write("\t".join(cols) + "\n")

    genes_by_seq: dict = {}
    for g in genes:
        genes_by_seq.setdefault(g["seqid"], []).append(g)

    for seqid in sorted(genes_by_seq, key=natural_key):
        seq_genes = sorted(genes_by_seq[seqid], key=lambda g: g["start"])
        for g in seq_genes:
            write_feature(g)
            passthrough_seen.add(id(g))
            transcripts = sorted(children_of_old_id.get(g["old_id"], []),
                                 key=lambda t: t["start"])
            for t in transcripts:
                if id(t) in passthrough_seen:
                    continue  # already written under an earlier multi-parent gene
                write_feature(t)
                passthrough_seen.add(id(t))
                subs = sorted(children_of_old_id.get(t["old_id"], []),
                              key=lambda s: s["start"])
                for s in subs:
                    if id(s) in passthrough_seen:
                        continue  # shared feature (multi-parent Parent=), write once
                    write_feature(s)
                    passthrough_seen.add(id(s))

    # Anything not written above (passthrough features not part of the
    # gene hierarchy) -- emitted afterward, in original file order.
    for f in features:
        if id(f) not in passthrough_seen:
            write_feature(f)


# ── Main ──────────────────────────────────────────────────────────────────────

def main(argv=None):
    ap = argparse.ArgumentParser(
        prog="GFF3RenameGenes",
        description="Systematically rename gene models in a GFF3 file "
                     "using a fixed SeqID-based numbering scheme.",
    )
    ap.add_argument("--gff", required=True, type=Path,
                    help="Input GFF3 file")
    ap.add_argument("--output", type=Path, default=None,
                    help="Output GFF3 file (default: stdout)")
    ap.add_argument("--gene_type", default="gene",
                    help="Column-3 feature type treated as a gene "
                         "(default: gene)")
    ap.add_argument("--gene_pad", type=int, default=DEFAULT_GENE_PAD,
                    help=f"Zero-padding width for gene numbers "
                         f"(default: {DEFAULT_GENE_PAD})")
    ap.add_argument("--gene_step", type=int, default=DEFAULT_GENE_STEP,
                    help=f"Increment between consecutive gene numbers "
                         f"(default: {DEFAULT_GENE_STEP})")
    ap.add_argument("--gene_start", type=int, default=DEFAULT_GENE_START,
                    help=f"First gene number on each SeqID "
                         f"(default: {DEFAULT_GENE_START})")
    ap.add_argument("--transcript_pad", type=int, default=DEFAULT_TRANSCRIPT_PAD,
                    help=f"Zero-padding width for transcript numbers "
                         f"(default: {DEFAULT_TRANSCRIPT_PAD})")
    ap.add_argument("--feature_pad", type=int, default=DEFAULT_FEATURE_PAD,
                    help=f"Zero-padding width for exon/CDS/UTR numbers "
                         f"(default: {DEFAULT_FEATURE_PAD})")
    ap.add_argument("--prefix_geneid", default=None,
                    help="Replace the literal SeqID with this fixed prefix "
                         "in every gene ID (default: use each feature's own "
                         "SeqID). Since the prefix is then no longer unique "
                         "per SeqID, gene numbering switches from "
                         "restarting at --gene_start on every new SeqID to "
                         "one continuous count across the whole file, so "
                         "IDs stay unique.")
    ap.add_argument("--after_seqid_tag", default=None,
                    help="Insert this tag between the SeqID (or "
                         "--prefix_geneid) and the 'G' in gene IDs, e.g. "
                         "PhangC01ANN2G000010 with --after_seqid_tag ANN2 -- "
                         "useful when the same assembly has multiple "
                         "annotations to keep their gene IDs distinguishable "
                         "(default: no tag). Does not affect numbering.")
    ap.add_argument("--dry_run", action="store_true",
                    help="Parse and report counts, then exit without "
                         "writing any output")
    ap.add_argument("--skip_sanity_check", action="store_true",
                    help="Skip the pre-renaming GFF3 structural sanity "
                         "check (duplicate IDs, dangling Parents, bad "
                         "coordinates, etc.)")
    ap.add_argument("--check_only", action="store_true",
                    help="Run the sanity check and print the report, "
                         "then exit without renaming or writing output")
    ap.add_argument("--force", action="store_true",
                    help="Proceed with renaming even if the sanity check "
                         "finds error-level problems (renaming may be "
                         "unreliable for the affected features)")
    ap.add_argument("--version", action="version", version=f"%(prog)s {VERSION}")
    args = ap.parse_args(argv)

    if not args.gff.exists():
        print(f"ERROR: --gff file not found: {args.gff}", file=sys.stderr)
        sys.exit(1)

    pragma_lines, features = read_gff(args.gff)
    print(f"Loaded {len(features)} feature(s) from {args.gff.name}",
          file=sys.stderr)

    n_genes_in = sum(1 for f in features if f["ftype"] == args.gene_type)
    if n_genes_in == 0:
        print(f"ERROR: no features of type '{args.gene_type}' found "
              f"(--gene_type). Nothing to rename.", file=sys.stderr)
        sys.exit(1)

    classified = classify_features(features, args.gene_type)

    if not args.skip_sanity_check:
        print(f"Sanity check: {len(features)} feature(s), "
              f"{len(classified['genes'])} gene(s), "
              f"{len(classified['transcripts'])} transcript(s), "
              f"{len(classified['subfeatures'])} exon/CDS/UTR feature(s), "
              f"{len(classified['passthrough'])} passthrough feature(s)",
              file=sys.stderr)
        errors, warnings = validate_gff(features, classified)
        for w in warnings:
            print(f"  SANITY WARNING: {w}", file=sys.stderr)
        for e in errors:
            print(f"  SANITY ERROR: {e}", file=sys.stderr)
        if errors and not args.force and not args.check_only:
            print(f"\nERROR: {len(errors)} structural problem(s) found -- "
                  f"renaming would likely be unreliable. Fix the input GFF3, "
                  f"or rerun with --force to proceed anyway (affected "
                  f"features may end up with an incorrect Parent/OldFeatID), "
                  f"or --skip_sanity_check to bypass this check entirely.",
                  file=sys.stderr)
            sys.exit(1)
        elif errors and args.force:
            print(f"  --force set: proceeding despite {len(errors)} "
                  f"error(s)", file=sys.stderr)
        elif not errors and not warnings:
            print("  No issues found.", file=sys.stderr)
    else:
        print("Sanity check skipped (--skip_sanity_check)", file=sys.stderr)

    if args.check_only:
        print("Exiting (--check_only). No output written.", file=sys.stderr)
        return

    if args.dry_run:
        print(f"  Genes found ({args.gene_type}): {n_genes_in}", file=sys.stderr)
        if args.prefix_geneid is not None:
            print(f"  Gene numbering: continuous (--prefix_geneid set) "
                  f"start={args.gene_start} step={args.gene_step} "
                  f"pad={args.gene_pad}", file=sys.stderr)
            print(f"  Gene ID prefix: '{args.prefix_geneid}' "
                  f"(replaces SeqID)", file=sys.stderr)
        else:
            print(f"  Gene numbering: restarts per SeqID "
                  f"start={args.gene_start} step={args.gene_step} "
                  f"pad={args.gene_pad}", file=sys.stderr)
        if args.after_seqid_tag:
            print(f"  After-SeqID tag: '{args.after_seqid_tag}'", file=sys.stderr)
        print(f"  Transcript numbering: pad={args.transcript_pad}", file=sys.stderr)
        print(f"  Exon/CDS/UTR numbering: pad={args.feature_pad}", file=sys.stderr)
        print("  Exiting (--dry_run). No output written.", file=sys.stderr)
        return

    stats = rename_features(
        classified,
        gene_type=args.gene_type,
        gene_pad=args.gene_pad,
        gene_step=args.gene_step,
        gene_start=args.gene_start,
        transcript_pad=args.transcript_pad,
        feature_pad=args.feature_pad,
        prefix_geneid=args.prefix_geneid,
        after_seqid_tag=args.after_seqid_tag,
    )

    if args.output:
        with open(args.output, "w") as out_fh:
            write_gff(pragma_lines, features, args.gene_type, out_fh)
        print(f"Written: {args.output}", file=sys.stderr)
    else:
        write_gff(pragma_lines, features, args.gene_type, sys.stdout)

    print(f"Genes: {stats['n_genes']}  (across {stats['n_seqids']} SeqID(s))  |  "
          f"Transcripts: {stats['n_transcripts']}  |  "
          f"Subfeatures: {stats['n_subfeatures']}  |  "
          f"Passthrough (unrenamed): {stats['n_passthrough']}", file=sys.stderr)
    if stats["n_genes_no_id"]:
        print(f"WARNING: {stats['n_genes_no_id']} gene(s) had no original ID "
              f"(OldFeatID fell back to seqid:start-end)", file=sys.stderr)
    if stats["n_transcripts_no_id"]:
        print(f"WARNING: {stats['n_transcripts_no_id']} transcript(s) had no "
              f"original ID", file=sys.stderr)
    if stats["n_sub_no_id"]:
        print(f"WARNING: {stats['n_sub_no_id']} exon/CDS/UTR feature(s) had "
              f"no original ID", file=sys.stderr)
    if stats["n_other_type"]:
        print(f"WARNING: {stats['n_other_type']} subfeature(s) had a type "
              f"other than exon/CDS/UTR; suffixed with the first 3 letters "
              f"of their own type instead", file=sys.stderr)
    if stats["n_unresolved_parent"]:
        print(f"WARNING: {stats['n_unresolved_parent']} Parent reference(s) "
              f"pointed to an ID outside the renamed gene/transcript/"
              f"subfeature hierarchy; left unchanged", file=sys.stderr)


if __name__ == "__main__":
    main()
