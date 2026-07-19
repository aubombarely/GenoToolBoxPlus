# GenoToolBoxPlus — Project Notes

Collection of general-purpose, zero-dependency CLI scripts for genomics and
genome annotation tasks. Each script is self-contained and uses the Python
standard library only.

**Current version:** v0.0.1 (per script)

This project follows the shared coding blueprint at `../CLAUDE.md`.
Apply those standards to any new scripts added here.

---

## Scripts

| Script | Purpose |
|---|---|
| `scripts/FastaRename.py` | Rename sequence IDs in a FASTA file using a two-column TSV mapping |
| `scripts/FastaStats.py` | Compute per-sequence and summary statistics for a FASTA file |
| `scripts/GFA2FASTA.py` | Convert GFA assembly graph format to FASTA |
| `scripts/NCBI_DownloadGenome.py` | Download genome FASTA/GFF3 from NCBI (accessions.txt, same format as annotseba), optional SeqID renaming |

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
- [ ] `CHANGELOG.md` — add when scripts reach v0.1.0
- [ ] README version badge — add when CHANGELOG is in place

## FAIR compliance status

- [x] `LICENSE` — MIT, added 2026-06-24
- [x] `CITATION.cff` — author, ORCID, version, keywords, repository URL
- [ ] Zenodo DOI — mint after first public release; add `doi:` field to `CITATION.cff`
- [ ] bio.tools registration — register after first public release
