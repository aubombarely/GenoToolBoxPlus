# GenoToolBoxPlus

A collection of general-purpose command-line scripts for genomics and genome annotation tasks.

## Index

- [FastaRename.py](#fastarenamepy) — Rename sequence IDs in a FASTA file using a two-column TSV mapping
- [FastaStats.py](#fastastatspy) — Compute per-assembly and nucleotide-composition statistics for a FASTA file
- [GFA2FASTA.py](#gfa2fastapy) — Convert GFA (v1 or v2) assembly graph segments to FASTA
- [NCBI_DownloadGenome.py](#ncbi_downloadgenomepy) — Download genome FASTA/GFF3 from NCBI with optional SeqID renaming
- [GetFasta4EarlGreyGFF.py](#getfasta4earlgreygffpy) — Extract FASTA sequences for TE features from an EarlGrey GFF3
- [GFF3RenameGenes.py](#gff3renamegenespy) — Systematically rename gene models in a GFF3 file

## Requirements

- Python 3.9+
- No external dependencies (standard library only)

## Scripts

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

### GetFasta4EarlGreyGFF.py

Extract FASTA sequences for TE features from an
[EarlGrey](https://github.com/TobyBaril/EarlGrey) repeat-annotation GFF3,
where column 3 is the TE type (e.g. `LINE/L1`, `LTR/Copia`) and the
attributes carry an `ID=` (the repeat family ID, e.g. `RND-1_FAMILY-789`).

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
