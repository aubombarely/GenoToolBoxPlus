# GenoToolBoxPlus — Project Notes

Collection of general-purpose, zero-dependency CLI scripts for genomics and
genome annotation tasks. Each script is self-contained and uses the Python
standard library only.

**Current version:** v0.2.1 (overall repo/citation version, tracked in
`CITATION.cff` and `CHANGELOG.md`; each script also carries its own
independent `VERSION` string for `--version`/troubleshooting)

This project follows the shared coding blueprint at `../CLAUDE.md`.
Apply those standards to any new scripts added here.

---

## Scripts

**FASTA_Utilities**

| Script | Purpose |
|---|---|
| `scripts/FastaRename.py` | Rename sequence IDs in a FASTA file using a two-column TSV mapping |
| `scripts/FastaStats.py` | Compute per-sequence and summary statistics for a FASTA file |
| `scripts/GFA2FASTA.py` | Convert GFA assembly graph format to FASTA |

**GenomicData_Download**

| Script | Purpose |
|---|---|
| `scripts/NCBI_DownloadGenome.py` | Download genome FASTA/GFF3 from NCBI (accessions.txt, same format as annotseba), optional SeqID renaming |

**GFF_Utilities**

| Script | Purpose |
|---|---|
| `scripts/GFF3RenameGenes.py` | Systematically rename gene/transcript/exon/CDS/UTR IDs in a GFF3 using a fixed SeqID-based numbering scheme, preserving old IDs as OldFeatID= |

**ThirdPartyTool_Utilities**

| Script | Purpose |
|---|---|
| `scripts/GetFasta4EarlGreyGFF.py` | Extract FASTA sequences for TE features from an EarlGrey repeat-annotation GFF3, strand-aware, sanitized headers |
| `scripts/GFF2BEDOrthoVenn.py` | Convert a GFF3 file to the 5-column BED format (SeqID, GeneID, Start, End, Strand) expected by OrthoVennPlus |
| `scripts/GAQET2AHRD.py` | Parse a GAQET run's GAQET.log.txt for its TREMBL/SWISSPROT diamond commands, build an AHRD YAML config, and run AHRD |

## Design principles for this collection

- **No external dependencies** — standard library only; each script must run
  in any Python 3.10+ environment without conda.
- **Single-file scripts** — each script is fully self-contained; no shared
  utility modules.
- **Streaming I/O** — always process FASTA/GFA line-by-line; never load an
  entire genome into memory.
- **Structured output not required** — these are small utilities; the
  `results/workdir/logs/` layout is overkill. A single `--output` file
  argument (defaulting to stdout) is sufficient.
- **Run log not required** — too heavyweight for one-liner utilities.
  Add `VERSION`, `--version`, and clear `--help` instead.

## Blueprint compliance status

- [x] `VERSION` string per script
- [x] `--version` argument
- [ ] `--dry_run` flag — add to new scripts
- [x] `CHANGELOG.md` — added at v0.2.1
- [x] README version badge — added at v0.2.1

## FAIR compliance status

- [x] `LICENSE` — MIT, added 2026-06-24
- [x] `CITATION.cff` — author, ORCID, version, keywords, repository URL
- [ ] Zenodo DOI — mint after first public release; add `doi:` field to `CITATION.cff`
- [ ] bio.tools registration — register after first public release
