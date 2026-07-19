# GenoToolBoxPlus

A collection of general-purpose command-line scripts for genomics and genome annotation tasks.

## Index

- [FastaRename.py](#fastarenamepy) — Rename sequence IDs in a FASTA file using a two-column TSV mapping
- [FastaStats.py](#fastastatspy) — Compute per-assembly and nucleotide-composition statistics for a FASTA file
- [GFA2FASTA.py](#gfa2fastapy) — Convert GFA (v1 or v2) assembly graph segments to FASTA
- [NCBI_DownloadGenome.py](#ncbi_downloadgenomepy) — Download genome FASTA/GFF3 from NCBI with optional SeqID renaming

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
```

**Arguments**

| Argument | Required | Description |
|----------|----------|-------------|
| `--accessions` | Yes | Accessions file: `species<TAB>accession<TAB>taxa_id[<TAB>prefix]` (same format as [annotseba](../annotseba)'s `accessions.txt`) |
| `--output` | Yes | Output directory (one subdirectory per accession) |
| `--rename_seqids` | No | Rename sequence IDs (see scheme below) and apply the same renaming to the GFF3, if one was downloaded |
| `--rename_prefix` | No | Default prefix for renamed sequence IDs when not set in the accessions file (default: `Sp`) |
| `--strip_description` | No | Drop the FASTA header description text when renaming (requires `--rename_seqids`), leaving just `>{new_id}` |
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
{output}/{species}_{accession}/
├── {species}_{accession}.fasta
├── {species}_{accession}.gff3             (only if annotation exists)
└── {species}_{accession}.equiv_seqID.txt  (only with --rename_seqids)
```

**SeqID renaming scheme (`--rename_seqids`)**

Inferred from each FASTA header's description text:

| Category | New ID | Number source |
|---|---|---|
| Chromosome | `{prefix}C##` | Parsed from the description (`chromosome 1`, `chromosome I`, `chromosome X`) |
| Mitochondrion | `{prefix}MIT##` | Order of appearance |
| Chloroplast / plastid | `{prefix}PLT##` | Order of appearance |
| Scaffold | `{prefix}SCF##` | Order of appearance |
| Anything else (unplaced/unlocalized contigs, etc.) | `{prefix}CTG##` | Order of appearance |

Chromosome numbers accept Arabic digits (`chromosome 1`) or Roman numerals
(`chromosome I`, common in yeast/fungal genomes) — Roman numerals are
auto-detected per genome (not per sequence) and converted to Arabic for
zero-padding. Single-letter sex chromosomes (`X`/`Y`/`W`/`Z`) are kept
literal (e.g. `{prefix}CX`), unless the rest of the genome is clearly
Roman-numbered, in which case a lone `X` is treated as Roman numeral 10.
Zero-padding width is the width of the largest number in that category
(minimum 2 digits); non-chromosome categories are numbered by order of
appearance in the file, not by any number in their own description.

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
