# GenoToolBoxPlus

A collection of general-purpose command-line scripts for genomics and genome annotation tasks.

## Index

- [FastaRename.py](#fastarenamepy) — Rename sequence IDs in a FASTA file using a two-column TSV mapping
- [FastaStats.py](#fastastatspy) — Compute per-assembly and nucleotide-composition statistics for a FASTA file
- [GFA2FASTA.py](#gfa2fastapy) — Convert GFA (v1 or v2) assembly graph segments to FASTA
- [NCBI_DownloadGenome.py](#ncbi_downloadgenomepy) — Download genome FASTA/GFF3 from NCBI with optional SeqID renaming
- [GetFasta4EarlGreyGFF.py](#getfasta4earlgreygffpy) — Extract FASTA sequences for TE features from an EarlGrey GFF3

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
