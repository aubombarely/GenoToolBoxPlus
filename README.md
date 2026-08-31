# GenoToolBoxPlus

<img src="https://img.shields.io/badge/version-v0.3.0-teal"/> <img src="https://img.shields.io/badge/python-3.9%2B-blue"/> <img src="https://img.shields.io/badge/platform-Linux%20%7C%20macOS-lightgrey"/> [Changelog](CHANGELOG.md)

A collection of general-purpose command-line scripts for genomics and genome annotation tasks. See [`CITATION.cff`](CITATION.cff) for how to cite this collection.

## Index

**FASTA_Utilities**
- [FastaRename.py](#fastarenamepy) — Rename sequence IDs in a FASTA file using a two-column TSV mapping
- [FastaStats.py](#fastastatspy) — Compute per-assembly and nucleotide-composition statistics for a FASTA file
- [GFA2FASTA.py](#gfa2fastapy) — Convert GFA (v1 or v2) assembly graph segments to FASTA

**GenomicData_Download**
- [NCBI_DownloadGenome.py](#ncbi_downloadgenomepy) — Download genome FASTA/GFF3 from NCBI with optional SeqID renaming
- [GWH4FastaRename.py](#gwh4fastarenamepy) — Rename SeqIDs in a local FASTA using GWH headers already present in the file, no download

**GFF_Utilities**
- [GFF3RenameGenes.py](#gff3renamegenespy) — Systematically rename gene models in a GFF3 file

**ThirdPartyTool_Utilities**
- [GetFasta4EarlGreyGFF.py](#getfasta4earlgreygffpy) — Extract FASTA sequences for TE features from an EarlGrey GFF3
- [GFF2BEDOrthoVenn.py](#gff2bedorthovennpy) — Convert a GFF3 file to the 5-column BED format expected by OrthoVennPlus
- [GAQET2AHRD.py](#gaqet2ahrdpy) — Build an AHRD config from a GAQET run and run AHRD

## Requirements

- Python 3.9+
- No external dependencies (standard library only)

## Scripts

<br>

---

---

## <img src="https://img.shields.io/badge/-FASTA__Utilities-4C9BE8?style=for-the-badge" height="42" alt="FASTA_Utilities">

- [FastaRename.py](#fastarenamepy) — Rename sequence IDs in a FASTA file using a two-column TSV mapping
- [FastaStats.py](#fastastatspy) — Compute per-assembly and nucleotide-composition statistics for a FASTA file
- [GFA2FASTA.py](#gfa2fastapy) — Convert GFA (v1 or v2) assembly graph segments to FASTA

<br>

### FastaRename.py

Rename sequence IDs in a FASTA file using a two-column TSV mapping.

**Usage**

```bash
FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv
FastaRename.py --fasta genome.fasta --tsv id_mapping.tsv --output renamed.fasta
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--fasta` | Yes | Input FASTA file |
| `--tsv` | Yes | Two-column TSV mapping: `old_seqid <TAB> new_seqid` |
| `--output` | No | Output FASTA file (default: stdout) |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**TSV format**

```
# Lines starting with # are ignored
chr1_old    chr1
chr2_old    chr2
scaffold_1  chr3
```

**Notes**
- Only the ID portion of the header (up to the first space) is renamed; description text is preserved.
- SeqIDs not found in the mapping are kept unchanged with a warning.
- Output sequences are wrapped at 60 characters per line.
- Stats (renamed/unchanged counts) are printed to stderr.

<br>

---

### FastaStats.py

Compute per-assembly and nucleotide-composition statistics for a FASTA file.

**Usage**

```bash
FastaStats.py --fasta genome.fasta
FastaStats.py --fasta genome.fasta --output genome_stats
FastaStats.py --fasta genome.fasta --output genome_stats --format tsv,txt,json
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--fasta` | Yes | Input FASTA file |
| `--output` | No | Output basename (default: derived from `--fasta` filename, stripping `.fasta`/`.fa`/`.fna`/`.ffn`/`.faa`/`.frn`) |
| `--format` | No | Comma-separated output formats: `tsv`, `txt`, `json` (default: `tsv`) |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Assembly statistics reported**

`num_sequences`, `total_length`, `min_length`, `max_length`, `mean_length`,
`median_length`, `n50`, `l50`, `n90`, `l90`, `gc_content`, `n_count`,
`n_percent`, `seq_gt_1kb`, `seq_gt_10kb`, `seq_gt_100kb`, `seq_gt_1mb`.

**Nucleotide composition reported**

Count and percentage for every IUPAC base present (`A C G T U R Y M K S W H
B V D N`); any other character is reported under `other`.

**Notes**
- One `--output` file is written per requested `--format` (e.g. `--format
  tsv,json` writes both `<basename>.tsv` and `<basename>.json`).
- A human-readable summary (both tables) is always printed to stderr,
  regardless of `--format`.
- GC content, N content, and per-base percentages are computed over total
  assembly length.

<br>

---

### GFA2FASTA.py

Extract segment sequences from a GFA (v1 or v2) assembly graph and write them
to a FASTA file.

**Usage**

```bash
GFA2FASTA.py --input assembly.gfa --output assembly.fasta
GFA2FASTA.py --input assembly.gfa --output assembly.fasta --summary summary.txt
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--input` | Yes | Input GFA file |
| `--output` | Yes | Output FASTA file |
| `--summary` | No | Optional plain-text summary report |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Notes**
- GFA version is auto-detected from the `H` header line's `VN:Z:` tag; if no
  header is present, it is inferred from the shape of the first `S` line
  (GFA2 `S` lines carry a numeric length field before the sequence).
- Segments with no sequence (`*`) are skipped with a warning and excluded
  from both the FASTA output and the total-length count.
- GFA2 sequence orientation markers (`+`/`-`) are stripped before writing.
- Output sequences are wrapped at 60 characters per line.
- Stats (segments written/skipped, total length) are printed to stderr.


<br>

---

---

## <img src="https://img.shields.io/badge/-GenomicData__Download-F5A623?style=for-the-badge" height="42" alt="GenomicData_Download">

- [NCBI_DownloadGenome.py](#ncbi_downloadgenomepy) — Download genome FASTA/GFF3 from NCBI with optional SeqID renaming
- [GWH4FastaRename.py](#gwh4fastarenamepy) — Rename SeqIDs in a local FASTA using GWH headers already present in the file, no download

<br>

### NCBI_DownloadGenome.py

Download genome FASTA (and GFF3, if available) from NCBI for a list of
accessions, with optional systematic sequence ID renaming.

**Usage**

```bash
NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/
NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ --rename_seqids
NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ --rename_seqids --strip_description
NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ --dry_run
NCBI_DownloadGenome.py --accessions accessions.txt --output genomes/ --rename_seqids --report_metrics
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--accessions` | Yes | Accessions file: `species<TAB>accession<TAB>taxa_id[<TAB>prefix]` (same format as [annotseba](../annotseba)'s `accessions.txt`) |
| `--output` | Yes | Output directory (one subdirectory per accession) |
| `--rename_seqids` | No | Rename sequence IDs (see scheme below) and apply the same renaming to the GFF3, if one was downloaded |
| `--rename_prefix` | No | Default prefix for renamed sequence IDs when not set in the accessions file (default: `Sp`) |
| `--strip_description` | No | Drop the FASTA header description text when renaming (requires `--rename_seqids`), leaving just `>{new_id}` |
| `--report_metrics` | No | Compute simple assembly metrics for each processed accession; writes `{output}/summary.tsv` and prints an ASCII table (see below) |
| `--force` | No | Re-download and re-process even if output already exists |
| `--dry_run` | No | Validate the accessions file and print what would be downloaded, then exit without any network calls |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Accessions file format**

```
# Format: species<TAB>accession<TAB>taxa_id[<TAB>prefix]
Saccharomyces_cerevisiae	GCF_000146045.2	4932	Sc
Arabidopsis_thaliana	GCA_000001735.4	3702	At
```
`taxa_id` is accepted for compatibility with annotseba's file but is unused
by this script. `taxa_id`/`prefix` may be `NA` or omitted; `prefix` then
falls back to `--rename_prefix`.

**Output layout**

```
{output}/
├── {species}_{accession}/
│   ├── {species}_{accession}.fasta
│   ├── {species}_{accession}.gff3             (only if annotation exists)
│   └── {species}_{accession}.equiv_seqID.txt  (only with --rename_seqids)
└── summary.tsv                                (only with --report_metrics)
```

**SeqID renaming scheme (`--rename_seqids`)**

| Category | New ID | How it's detected | Number source |
|---|---|---|---|
| Chromosome | `{prefix}C##` | Description says chromosome / pseudomolecule / linkage group / `LG` | Parsed from the description as-is |
| Mitochondrion | `{prefix}MIT##` | Description says mitochondrion | Order of appearance |
| Chloroplast / plastid | `{prefix}PLT##` | Description says chloroplast / plastid | Order of appearance |
| Scaffold | `{prefix}SCF##` | Not chromosome/organelle, and the **sequence itself contains an N** | Order of appearance |
| Contig | `{prefix}CTG##` | Not chromosome/organelle, and the sequence has **no N** | Order of appearance |

Chromosome numbers accept Arabic digits (`chromosome 1`), Roman numerals
(`chromosome I`, common in yeast/fungal genomes — auto-detected per genome,
not per sequence, and converted to Arabic for zero-padding), single-letter
sex chromosomes (`X`/`Y`/`W`/`Z`, kept literal e.g. `{prefix}CX` unless the
rest of the genome is clearly Roman-numbered), and letter-suffixed polyploid
labels (`1A`, `2B`, `3D`, common in wheat-like subgenomes — kept literal,
original numbering preserved, since these aren't a plain integer). Chromosome
zero-padding width is the width of the largest plain number in that genome
(minimum 2 digits); SCF/CTG/MIT/PLT are numbered by order of appearance in
the file, not by any number in their own description.

The SCF-vs-CTG split is decided from actual sequence content, not the
header text — a scaffold-like keyword in the description is not required
(and isn't checked); what matters is whether the sequence contains any `N`.

**Sanity checks**

Every FASTA scan (for renaming and/or `--report_metrics`) verifies all
SeqIDs are unique, aborting that accession with a clear error otherwise —
duplicate IDs would otherwise silently corrupt classification, metrics, and
renaming, since results are keyed by SeqID. It also flags any non-standard
IUPAC ambiguity codes found (anything besides `A`/`C`/`G`/`T`/`N`, e.g. `Y`,
`R`, `W`), which can break aligners or variant callers that assume a plain
4-letter(+N) alphabet.

**Assembly metrics (`--report_metrics`)**

For each processed accession (whether just downloaded or already present
from a previous run): `seq_n`, `assembly_size`, `avg_length`, `n50`, `l50`,
`n90`, `l90`, per-category counts `n_chr`/`n_scf`/`n_ctg`/`n_mit`/`n_plt`,
`ambig_nt_count`/`ambig_nt_chars` (non-standard IUPAC codes found, if any),
and `annotation_gff` (`YES`/`NO`, whether a GFF3 was downloaded). Written as
`{output}/summary.tsv` (one row per accession) and printed as an ASCII table
to stderr, e.g.:

```
+--------------------------+-----------------+-------+---------------+------------+--------+-----+--------+-----+-----+-----+-----+-----+-----+----------+-------------+------------------+
| Species                  |       Accession | Seq_N | Assembly_size | Avg_length |    N50 | L50 |    N90 | L90 | CHR | SCF | CTG | MIT | PLT | Ambig_NT | Ambig_chars | Annotation (GFF) |
+--------------------------+-----------------+-------+---------------+------------+--------+-----+--------+-----+-----+-----+-----+-----+-----+----------+-------------+------------------+
| Saccharomyces_cerevisiae | GCF_000146045.2 |    17 |      12157105 |  715123.82 | 924431 |   6 | 439888 |  13 |  16 |   0 |   0 |   1 |   0 |        0 |             |              YES |
+--------------------------+-----------------+-------+---------------+------------+--------+-----+--------+-----+-----+-----+-----+-----+-----+----------+-------------+------------------+
```

On a fresh download, metrics are computed from the original downloaded
FASTA (before any renaming), so `--rename_seqids --strip_description`
together can't blind CHR/MIT/PLT classification by removing the header
descriptions it depends on. On an accession skipped as already-downloaded
(no fresh copy available), metrics fall back to the file in `--output`; if
that file was previously written with `--strip_description`, CHR/MIT/PLT
counts can't be recovered from it and a warning is printed — use `--force`
to recompute from a fresh download in that case.

**Notes**
- Downloads via the [NCBI Datasets REST API v2](https://www.ncbi.nlm.nih.gov/datasets/docs/v2/reference-docs/rest-api/) directly (`urllib`/`zipfile`, both standard library) — no `datasets` CLI tool required.
- If an accession has no GFF3 annotation on NCBI, the FASTA is still downloaded and a note (not an error) is printed to stderr.
- Already-downloaded accessions (non-empty `.fasta` already present) are skipped unless `--force` is given.
- Retries each download up to 3 times on transient network errors.
- If your Python installation's default certificate store is broken (a
  known issue with python.org's macOS installer, which ships without a
  populated root CA bundle) downloads will fail with an SSL verification
  error; the script automatically falls back to the `certifi` package's
  bundle if it happens to be installed, and otherwise prints the exact fix
  (`Install Certificates.command`, or `pip install certifi`).
- Stats (sequences renamed, files written) are printed to stderr.

<br>

### GWH4FastaRename.py

Rename sequence IDs in a FASTA file whose headers are already in GWH
(Genome Warehouse) format. Unlike `NCBI_DownloadGenome.py`, this script
does not download anything — it operates on a FASTA file you already have
locally, and derives the new SeqID/description from the GWH header fields
already present in that file.

**Usage**

```bash
GWH4FastaRename.py --fasta genome.gwh.fasta --output renamed.fasta
GWH4FastaRename.py --fasta genome.gwh.fasta --output renamed.fasta --prefix Sp
GWH4FastaRename.py --fasta genome.gwh.fasta --output renamed.fasta --equiv equiv.tsv
GWH4FastaRename.py --fasta genome.gwh.fasta --dry_run
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--fasta` | Yes | Input FASTA file with GWH-format headers |
| `--output` | No | Output renamed FASTA file (default: stdout) |
| `--equiv` | No | Equivalence TSV path (`gwh_id`, `oriseqid`, `category`, `new_seqid`). Default: `{output}.equiv_seqID.txt` when `--output` is given; skipped (with a note) when writing to stdout and `--equiv` isn't set |
| `--prefix` | No | Prefix for renamed sequence IDs (default: `Sp`) |
| `--dry_run` | No | Parse and classify the FASTA, print what would be renamed, then exit without writing any output |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**GWH header format**

Fields are whitespace-separated; the optional free-text category and
`Complete=`/`Circular=` fields vary by submission, but `OriSeqID=` is
required on every record:

```
>GWHBJVE00000001    OriSeqID=Chr13    Len=38709045
>GWHDUBQ00000001    Chromosome 1    Complete=T    Circular=F    OriSeqID=Chr01    Len=45819973
>GWHDUBQ00000014    OriSeqID=contig0001    Len=45917
```

**SeqID renaming scheme**

Same CHR/SCF/CTG/MIT/PLT categories, numbering, and zero-padding rules as
`NCBI_DownloadGenome.py` (see above), with the chromosome number sourced
from GWH fields instead of an NCBI-style description:

| Category | New ID | How it's detected | Number source |
|---|---|---|---|
| Chromosome | `{prefix}C##` | Free text says "Chromosome N", **or** `OriSeqID` is a strict `Chr<digits>` label (e.g. `Chr01`) | Parsed from whichever matched |
| Mitochondrion | `{prefix}MIT##` | Free text or `OriSeqID` says mitochondrion | Order of appearance |
| Chloroplast / plastid | `{prefix}PLT##` | Free text or `OriSeqID` says chloroplast / plastid | Order of appearance |
| Scaffold | `{prefix}SCF##` | Not chromosome/organelle, and the **sequence itself contains an N** | Order of appearance |
| Contig | `{prefix}CTG##` | Not chromosome/organelle, and the sequence has **no N** | Order of appearance |

`OriSeqID` labels like `ChrUN10` (unplaced/unlocalized chromosome-scale
sequences, a common GWH convention) are **not** treated as chromosomes —
the letters after `Chr` make the label fail the strict `Chr<digits>`
pattern, so these fall through to the SCF/CTG split by sequence content
like any other non-chromosome sequence.

**Output**

Renamed FASTA headers are `>{new_id}` with the original description
stripped. The equivalence table has one row per sequence:

```
gwh_id	oriseqid	category	new_seqid
GWHBJVE00000001	Chr13	CHR	SpC13
GWHBJVE00000002	Chr4	CHR	SpC04
GWHBJVE00000014	ChrUN10	CTG	SpCTG01
```

**Sanity checks**

Every GWH ID must be unique and every record must have an `OriSeqID=`
field; both are validated in a single streaming pass over the file before
any renaming happens, aborting with a clear error otherwise.

**Notes**
- Standard library only — no network access, no external dependencies.
- Streams the FASTA line-by-line; never holds a full sequence in memory.

<br>

---

---

## <img src="https://img.shields.io/badge/-GFF__Utilities-888888?style=for-the-badge" height="42" alt="GFF_Utilities">

- [GFF3RenameGenes.py](#gff3renamegenespy) — Systematically rename gene models in a GFF3 file

<br>

### GFF3RenameGenes.py

Systematically rename every gene, transcript, exon, CDS, and UTR ID in a
GFF3 file using a fixed, predictable scheme built from the sequence ID.

**Usage**

```bash
GFF3RenameGenes.py --gff annotation.gff3
GFF3RenameGenes.py --gff annotation.gff3 --output renamed.gff3
GFF3RenameGenes.py --gff annotation.gff3 --check_only
GFF3RenameGenes.py --gff annotation.gff3 --dry_run
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--gff` | Yes | Input GFF3 file |
| `--output` | No | Output GFF3 file (default: stdout) |
| `--gene_type` | No | Column-3 feature type treated as a gene (default: `gene`) |
| `--gene_pad` | No | Zero-padding width for gene numbers (default: 6) |
| `--gene_step` | No | Increment between consecutive gene numbers (default: 10) |
| `--gene_start` | No | First gene number on each SeqID (default: 10) |
| `--transcript_pad` | No | Zero-padding width for transcript numbers (default: 2) |
| `--feature_pad` | No | Zero-padding width for exon/CDS/UTR numbers (default: 2) |
| `--prefix_geneid` | No | Replace the literal SeqID with this fixed prefix in every gene ID (default: use each feature's own SeqID); switches gene numbering to one continuous count across the whole file, since the prefix is then no longer unique per SeqID |
| `--after_seqid_tag` | No | Insert this tag between the SeqID (or `--prefix_geneid`) and the `G`, e.g. `PhangC01ANN2G000010` with `--after_seqid_tag ANN2` — useful when the same assembly has multiple annotations (default: no tag); does not affect numbering |
| `--skip_sanity_check` | No | Skip the pre-renaming structural sanity check |
| `--check_only` | No | Run the sanity check and print the report, then exit without renaming |
| `--force` | No | Proceed with renaming even if the sanity check finds error-level problems |
| `--dry_run` | No | Parse and report counts, then exit without writing output |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Sanity check**

Runs automatically before renaming (deliberately generic, not tuned to
any one annotation tool, since BRAKER/MAKER/HELIXER/EVIANN/EGAPX/ANNEVO/
TransDecoder output all fail in different ways):

*Error-level* (abort by default, since renaming would be unreliable):
duplicate IDs, genes/transcripts with no ID, `Parent` referencing an ID
that doesn't exist anywhere in the file, invalid coordinates
(`start < 1` or `start > end`), and a feature on a different SeqID than
its `Parent` (a strong sign of a corrupted or badly merged GFF3).

*Warning-level* (reported, non-fatal): genes with no transcript
children, transcripts with no exon/CDS/UTR children, a feature whose
coordinates extend beyond its `Parent`'s range, a feature on a different
strand than its `Parent`, and features that look like exact duplicates
(same seqid/coordinates/strand/type/parent — a common artifact of
merging multiple annotation runs, e.g. combining several BRAKER
predictions; can also false-positive on two genuinely distinct
transcripts that happen to share the same outer span, so treat it as a
prompt to check, not a guaranteed problem).

Use `--check_only` to run just the check (e.g. as a pre-flight QC step
on a fresh tool output before deciding whether to rename at all),
`--force` to rename anyway despite errors (affected features may end up
with an incorrect `Parent`/`OldFeatID`), or `--skip_sanity_check` to
bypass the check entirely.

**Naming scheme**

```
Gene:       {SeqID}{tag}G{gene_number}          e.g. PhangC01G000010
Transcript: {gene_ID}T{transcript_number}       e.g. PhangC01G000010T01
Exon:       {transcript_ID}EXO{exon_number}     e.g. PhangC01G000010T01EXO01
CDS:        {transcript_ID}CDS{cds_number}      e.g. PhangC01G000010T01CDS01
UTR:        {transcript_ID}UTR{utr_number}      e.g. PhangC01G000010T01UTR01
```

- Gene numbers restart at 10 on every new SeqID and increment by 10 (10, 20, 30 ...), leaving gaps for later manual insertions — unless `--prefix_geneid` is set, in which case numbering becomes one continuous count across the whole file instead, so IDs built from a shared, non-SeqID-derived prefix stay unique.
- `--after_seqid_tag` inserts a fixed tag between the SeqID (or `--prefix_geneid`) and the `G`; it doesn't affect numbering either way.
- Transcript numbers restart at 01 for every gene, ordered by start coordinate.
- Exon, CDS, and UTR numbers are independent counters that each restart at 01 for every transcript, and follow transcript direction (ascending start on `+` strand, descending on `-` strand) — not raw genomic order.
- `five_prime_UTR` and `three_prime_UTR` share one continuous UTR numbering, 5' → 3' along the transcript.
- A "gene" is any column-3 type matching `--gene_type`. A transcript is any feature whose `Parent` points to a gene, regardless of its own column-3 type (`mRNA`, `tRNA`, `ncRNA`, ...). A subfeature is any feature whose `Parent` points to a transcript; an unrecognized subfeature type is suffixed with the first 3 letters of its own type instead of EXO/CDS/UTR, with a warning.
- Every other original attribute (`Note=`, `product=`, `Dbxref=`, ...) is preserved; only `ID`/`Parent` are replaced, and the original ID (or `{seqid}:{start}-{end}` if the feature had none) is appended as a trailing `OldFeatID=`.
- Features shared by multiple parents (`Parent=mRNA1,mRNA2`) are written once, with both renamed parents listed.
- Pragma lines, comments, and any feature outside the gene/transcript/subfeature hierarchy (e.g. a standalone `region` line) pass through unchanged, with no `OldFeatID` added.
- Output is grouped by SeqID (natural sort order) then by gene → transcript → subfeature, not raw input file order.
- Stats (gene/transcript/subfeature/passthrough counts, and warnings for missing IDs or unresolved Parent references) are printed to stderr.


<br>

---

---

## <img src="https://img.shields.io/badge/-ThirdPartyTool__Utilities-E8604C?style=for-the-badge" height="42" alt="ThirdPartyTool_Utilities">

- [GetFasta4EarlGreyGFF.py](#getfasta4earlgreygffpy) — Extract FASTA sequences for TE features from an EarlGrey GFF3
- [GFF2BEDOrthoVenn.py](#gff2bedorthovennpy) — Convert a GFF3 file to the 5-column BED format expected by OrthoVennPlus
- [GAQET2AHRD.py](#gaqet2ahrdpy) — Build an AHRD config from a GAQET run and run AHRD

<br>

### GetFasta4EarlGreyGFF.py

Extract FASTA sequences for TE features from an
[EarlGrey](https://github.com/TobyBaril/EarlGrey) repeat-annotation GFF3,
where column 3 is the TE type (e.g. `LINE/L1`, `LTR/Copia`) and the
attributes carry an `ID=` (the repeat family ID, e.g. `RND-1_FAMILY-789`).

**Rationale**

Earl Grey annotates TE locations in a GFF3 but doesn't itself export the
matching nucleotide sequences per family — that's a separate manual step
most pipelines need for downstream work (building/curating a repeat
library, running the family through a classifier, BLASTing a specific
insertion). This script closes that gap: it reads Earl Grey's GFF3 and
genome FASTA directly, so no intermediate `bedtools getfasta` command with
its own coordinate/strand bookkeeping is needed, and the output headers are
already tagged with family ID, TE type, and genomic location for
traceability back to the annotation.

**Usage**

```bash
GetFasta4EarlGreyGFF.py --fasta genome.fasta --gff repeats.gff3
GetFasta4EarlGreyGFF.py --fasta genome.fasta --gff repeats.gff3 --output TEs.fasta
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--fasta` | Yes | Input genome FASTA |
| `--gff` | Yes | EarlGrey GFF3 (column 3 = TE type; `ID=` attribute = family ID) |
| `--output` | No | Output FASTA file (default: stdout) |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Header format**

```
>{ID}_{TYPE}_{SeqID}_{Start}
```

`TYPE` has every `/` replaced with `_` (e.g. `LTR/Copia` → `LTR_Copia`) so
the header stays shell- and tool-safe. Example — this GFF3 line:

```
PhangAGP1C01  Earl_Grey  LINE/L1  1  2101  10800  +  .  TSTART=5686;TEND=7874;ID=RND-1_FAMILY-789;SHORTTE=F;KIMURA80=0.2841
```

produces:

```
>RND-1_FAMILY-789_LINE_L1_PhangAGP1C01_1
```

**Notes**
- Strand-aware: `-` strand features are reverse-complemented; `+`/`.`/unset are extracted forward.
- Extraction streams the genome FASTA sequence-by-sequence (peak memory = one chromosome), never loading the whole genome at once.
- Features with no `ID=` attribute fall back to `{seqid}_{start}_{end}` as the ID, with a warning.
- Features whose coordinates fall outside their sequence's length are skipped with a warning (not silently dropped or truncated).
- GFF3 seqids with no match in the genome FASTA are reported once as a count, not per-feature.
- Output sequences are wrapped at 60 characters per line.
- Extraction stats (extracted/skipped counts) are printed to stderr.

<br>

---

### GFF2BEDOrthoVenn.py

Convert a GFF3 file to the 5-column BED format expected by
[OrthoVennPlus](https://orthovenn3.bioinfotoolkits.net/):
`SeqID  GeneID  Start  End  Strand`.

**Rationale**

OrthoVennPlus needs gene coordinates in its own 5-column layout, not
standard GFF3 or 0-based BED, and every genome annotation pipeline in this
workspace outputs GFF3. Rather than hand-writing an `awk`/`cut` one-liner
per project (and re-deriving the right column order and 1-based
coordinates each time), this script does the conversion directly from the
GFF3 `gene` features, so a GAQET/annotation run's output can be fed into
an OrthoVennPlus comparison without an intermediate manual reformatting
step.

**Usage**

```bash
GFF2BEDOrthoVenn.py --gff annotation.gff3 --output annotation.bed
GFF2BEDOrthoVenn.py --gff annotation.gff3 --dry_run
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--gff` | Yes | Input GFF3 file |
| `--output` | No | Output BED file (default: stdout) |
| `--feature_type` | No | Column-3 feature type to extract (default: `gene`) |
| `--dry_run` | No | Parse and report counts, then exit without writing output |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Notes**
- `GeneID` is taken from the `ID=` attribute; features of `--feature_type` with no `ID=` are skipped with a warning.
- Coordinates are copied as-is from the GFF3 (1-based, inclusive) — this matches OrthoVennPlus's expected layout, not 0-based BED.
- Rows are sorted by `(SeqID, Start)`, SeqID in natural sort order.
- Stats (rows written, skipped-for-missing-ID count) are printed to stderr.

<br>

---

### GAQET2AHRD.py

Parse a [GAQET](https://github.com/victorgcb1987/GAQET2) run's `GAQET.log.txt`
for the exact TREMBL/SWISSPROT `diamond blastp` commands it used, build an
[AHRD](https://github.com/groupschoof/AHRD) YAML config from them, and (by
default) run AHRD via `java -jar $AHRD_JAR config.yml`.

**Rationale**

GAQET already runs the two diamond homology searches AHRD needs (TREMBL
and SWISSPROT) as part of its QC pipeline, but doesn't run AHRD itself or
generate its config — assembling an AHRD YAML by hand means re-typing the
same diamond output paths already sitting in `GAQET.log.txt`, guessing at
weight/scoring parameters, and remembering to wire up the GOA file for GO
term transfer, all error-prone and easy to get subtly wrong (e.g. omitting
`gene_ontology_result` silently produces a config that runs fine but
transfers zero GO terms). This script closes that integration gap: it
reuses the search results GAQET already computed instead of re-running
diamond, and encodes a validated set of AHRD parameters as defaults, so
functional annotation becomes a single command that follows directly from
a GAQET run rather than a separate, manually-configured step.

**Usage**

```bash
GAQET2AHRD.py --gaqet_log GAQET.log.txt --ahrd_jar /opt/ahrd/ahrd.jar --ahrd_home /opt/ahrd
AHRD_JAR=/opt/ahrd/ahrd.jar AHRD_HOME=/opt/ahrd GAQET2AHRD.py --gaqet_log GAQET.log.txt
GAQET2AHRD.py --gaqet_log GAQET.log.txt --skip_ahrd
GAQET2AHRD.py --gaqet_log GAQET.log.txt --dry_run
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--gaqet_log` | Yes | Path to a GAQET run's `GAQET.log.txt` |
| `--output` | No | Output directory for the config and AHRD results (default: `AHRD_run/` next to `--gaqet_log`) |
| `--ahrd_jar` | No | Path to `ahrd.jar` (default: `$AHRD_JAR` env var); required unless `--skip_ahrd` |
| `--ahrd_home` | No | Path to an AHRD checkout (default: `$AHRD_HOME` env var); when set, defaults the filter files below to `{ahrd_home}/test/resources/` |
| `--java_bin` | No | `java` executable to use (default: `java`) |
| `--java_xmx` | No | JVM max heap, e.g. `8g` (default: unset) |
| `--swissprot_weight` | No | `blast_dbs.swissprot.weight` (default: `653`) |
| `--trembl_weight` | No | `blast_dbs.trembl.weight` (default: `904`) |
| `--swissprot_desc_weight` | No | `blast_dbs.swissprot.description_score_bit_score_weight` (default: `2.717061`) |
| `--trembl_desc_weight` | No | `blast_dbs.trembl.description_score_bit_score_weight` (default: `2.590211`) |
| `--token_bit_score_weight` | No | `token_score_bit_score_weight` (default: `0.468`) |
| `--token_database_score_weight` | No | `token_score_database_score_weight` (default: `0.2098`) |
| `--token_overlap_score_weight` | No | `token_score_overlap_score_weight` (default: `0.3221`) |
| `--blacklist` | No | AHRD `blacklist_descline` file, shared by both DBs (default: derived from `--ahrd_home`, else omitted) |
| `--swissprot_filter` | No | AHRD `filter_descline` file for swissprot (default: derived from `--ahrd_home`, else omitted) |
| `--trembl_filter` | No | AHRD `filter_descline` file for trembl (default: derived from `--ahrd_home`, else omitted) |
| `--token_blacklist` | No | AHRD `blacklist_token` file, shared by both DBs (default: derived from `--ahrd_home`, else omitted) |
| `--gene_ontology_result` | No | GO Annotation (GAF) file for GO term transfer (default: `goa_uniprot_all.gaf` in the same directory as the SWISSPROT `--db` from `GAQET.log.txt`) |
| `--skip_go` | No | Do not transfer GO terms (omit `gene_ontology_result`/`reference_go_regex`/`prefer_reference_with_go_annos` from the config) |
| `--top_n` | No | Number of most abundant descriptions to list in the summary (default: `10`) |
| `--skip_summary` | No | Do not write `<prefix>_AHRD.summary.txt` or print the summary tables after AHRD finishes |
| `--check_te_goterms` | No | Cross-check AHRD's GO terms against DETENGA's TE calls; writes `<prefix>_TEGOterm_vs_DETENGA.tsv` |
| `--te_goterms_file` | No | TE-associated GO term list, one `GO:#######` per line (default: hardcoded list — see `--print_te_associated_default_goterms`) |
| `--detenga_csv` | No | DETENGA's `<basename>_TE_summary.csv` (default: `DETENGA_run/{prefix}_TE_summary.csv` next to `--gaqet_log`) |
| `--print_te_associated_default_goterms` | No | Print the hardcoded default TE-associated GO term list as a table, then exit (does not require `--gaqet_log`) |
| `--skip_ahrd` | No | Write the config only; do not invoke AHRD |
| `--dry_run` | No | Parse the log and print what would be written/run, then exit |
| `--version` | No | Show version and exit |
| `--help` | No | Show help and exit |

**Notes**
- AHRD needs the flat FASTA (with description headers) of each blast DB, not the diamond `.dmnd` index used for the search — derived by swapping the `--db` path's extension to `.fasta`; a warning is printed if that file isn't found, since the path may need manual correction in the generated config.
- The AHRD output filename prefix is taken from GAQET's own `{prefix}_GAQET.stats.tsv` file, found by globbing `--gaqet_log`'s directory; falls back to `AHRD` if not found.
- A warning (not an error) is printed if the log's success marker (`run successfully`) isn't found after either command, or if any referenced file (proteins FASTA, diamond output, derived DB FASTA, blacklist/filter files, GO GAF file) doesn't exist on disk — the config is still written so paths can be corrected by hand if needed.
- GO term transfer is on by default (`gene_ontology_result`/`reference_go_regex`/`prefer_reference_with_go_annos`); use `--skip_go` to omit it entirely.
- Default weights match a validated real-world AHRD config tuned for plant genome annotation; override any of them per-run as needed.

**Summary report**

After a successful AHRD run, `<prefix>_AHRD.summary.txt` is written next
to the AHRD output TSV, and the same summary is printed to stderr as ASCII
tables (skip both with `--skip_summary`):

```
AHRD functional annotation summary:
+--------------------------+-------+---------+
| Metric                   | Count | Percent |
+--------------------------+-------+---------+
| Total proteins           | 34521 |  100.0% |
| With description         | 29876 |   86.5% |
| Unknown / no description |  4645 |   13.5% |
| With GO term(s)          | 24103 |   69.8% |
+--------------------------+-------+---------+

AHRD-Quality-Code distribution:
+--------------+-------+---------+
| Quality Code | Count | Percent |
+--------------+-------+---------+
| ***          | 18204 |   52.7% |
| **-          |  6532 |   18.9% |
| ...          |   ... |     ... |
+--------------+-------+---------+

Top 10 most abundant descriptions:
+------+-------+--------------------------------------+
| Rank | Count | Description                          |
+------+-------+--------------------------------------+
|    1 |   412 | Putative disease resistance protein  |
|    2 |   298 | Reverse transcriptase                |
|  ... |   ... | ...                                   |
+------+-------+--------------------------------------+
```

Column positions (description, GO term, AHRD-Quality-Code) are located by
header keyword rather than a fixed index, since the exact column set
depends on the AHRD config (e.g. whether GO/InterPro were requested). A
protein counts as "with description" unless its description is empty or
starts with `Unknown protein` (AHRD's own placeholder for no hit).

**TE GO-term vs DETENGA cross-check (`--check_te_goterms`)**

Cross-checks AHRD's transferred GO terms against DETENGA's own dedicated
TE calls, to catch protein-coding gene models that are actually
TE-derived. Writes `<prefix>_TEGOterm_vs_DETENGA.tsv`:

| Column | Meaning |
|---|---|
| `ProteinID` | AHRD's `Protein-Accession` |
| `AHRD_GO_TEs` | Comma-separated TE-associated GO terms the protein has (`-` if it has GO terms but none are TE-associated; `NA` if it has no GO terms at all) |
| `DETENGA_TE` | DETENGA's `DeTEnGA_status` for this protein when it indicates a TE (`-` if DETENGA explicitly called it non-TE; `NA` if the protein isn't in `--detenga_csv` at all) |
| `GOTE_TAGGED` | `YES` if **at least one** GO term the protein has is TE-associated; `NO` if it has GO terms but none are TE-associated; `NA` if it has no GO terms |
| `DETENGA_TAGGED` | `YES`/`NO` from whether `DeTEnGA_status` contains `"te"` (e.g. `PteM0`, `P0Mte`, `PteMte` vs `PcpM0`); `NA` if the protein isn't in `--detenga_csv` |

`NA` always means "no info available" (protein missing from the relevant
source), never "checked and not a TE" — that case is `-`/`NO`.

An agreement summary is printed to stderr:

```
GOTE vs DETENGA agreement:
+-------------------------------+-------+
| Category                      | Count |
+-------------------------------+-------+
| Both TE (agree)               |  1834 |
| Both non-TE (agree)           | 28442 |
| GO-only TE (DETENGA: non-TE)  |   112 |
| DETENGA TE (GO: not TE-only)  |   890 |
| Incomplete info (either NA)   |  3138 |
+-------------------------------+-------+
```

Use `--print_te_associated_default_goterms` to see the hardcoded default
TE GO term list (transposase/reverse-transcriptase/integrase molecular
functions, transposition biological processes, retrotransposon
nucleocapsid/assembly components); override with `--te_goterms_file`
(one `GO:#######` per line).

**Notes**
- `GOTE_TAGGED` uses an "at least one" rule, not "only TE terms": real TE
  proteins routinely also carry generic companion GO terms (e.g. DNA
  binding, zinc ion binding) alongside the TE-specific ones, so requiring
  *every* GO term to be TE-associated essentially never matches in
  practice — a genuine `Transposase`-annotated protein can easily have
  `GO:0003677` (DNA binding) and `GO:0008270` (zinc ion binding) alongside
  `GO:0004803` (transposase activity).
- Both methods can flag **domesticated TE-derived genes** as false
  positives — genes that now have a real cellular function but retained
  TE-like domains or GO terms from their evolutionary origin. Treat
  agreement/disagreement as a prioritization signal for manual review,
  not a final call.
- DETENGA can catch TEs that AHRD's GO transfer misses entirely (e.g. a
  protein with no GO terms at all, where `GOTE_TAGGED=NA` but
  `DETENGA_TAGGED=YES`) — this is expected, since the two methods use
  independent evidence (reference GOA GO terms vs. TEsort/InterPro
  domain calls).

## Third-party tools and citations

These scripts don't bundle or depend on the tools below at import time (no
external Python dependencies, per the design principles above) — they
generate inputs/configs for them or parse their output. Cite the
corresponding tool if you use it via one of these scripts:

| Tool | Used by | Citation |
|---|---|---|
| [AHRD](https://github.com/groupschoof/AHRD) | `GAQET2AHRD.py` | Hallab A. et al. *AHRD — Automated Assignment of Human Readable Descriptions.* github.com/groupschoof/AHRD |
| [DIAMOND](https://github.com/bbuchfink/diamond) | `GAQET2AHRD.py` (consumes its output) | Buchfink B, Reuter K, Drost HG. Sensitive protein alignments at tree-of-life scale using DIAMOND. *Nat Methods.* 2021;18:366–368. doi:[10.1038/s41592-021-01101-x](https://doi.org/10.1038/s41592-021-01101-x) |
| [Earl Grey](https://github.com/TobyBaril/EarlGrey) | `GetFasta4EarlGreyGFF.py` | Baril T, Galbraith J, Hayward A. Earl Grey: A Fully Automated User-Friendly Transposable Element Annotation and Analysis Pipeline. *Mol Biol Evol.* 2024;41(4):msae068. doi:[10.1093/molbev/msae068](https://doi.org/10.1093/molbev/msae068) |
| [OrthoVenn3](https://orthovenn3.bioinfotoolkits.net/) | `GFF2BEDOrthoVenn.py` | Sun J. et al. OrthoVenn3: an integrated platform for exploring and visualizing orthologous data across genomes. *Nucleic Acids Res.* 2023;51(W1):W397–W403. doi:[10.1093/nar/gkad313](https://doi.org/10.1093/nar/gkad313) |
| [NCBI Datasets](https://www.ncbi.nlm.nih.gov/datasets/) | `NCBI_DownloadGenome.py` | O'Leary NA. et al. Exploring and retrieving sequence and metadata for species across the tree of life with NCBI Datasets. *Sci Data.* 2024;11:732. doi:[10.1038/s41597-024-03571-y](https://doi.org/10.1038/s41597-024-03571-y) |
| [GAQET2](https://github.com/victorgcb1987/GAQET2) | `GAQET2AHRD.py` (input format) | victorgcb1987. *GAQET2.* github.com/victorgcb1987/GAQET2 |
