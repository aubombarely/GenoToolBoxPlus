# Changelog

All notable changes to GenoToolBoxPlus are documented here. See `CITATION.cff`
for the version to cite.

## [v0.2.5] — 2026-08-25

### Added
- `GAQET2AHRD.py`: `--check_te_goterms` cross-checks AHRD's transferred GO
  terms against DETENGA's own TE calls (`DeTEnGA_status`), to catch
  protein-coding gene models that are actually TE-derived. Writes
  `<prefix>_TEGOterm_vs_DETENGA.tsv` (ProteinID, AHRD_GO_TEs, DETENGA_TE,
  GOTE_TAGGED, DETENGA_TAGGED) and prints an agreement summary.
  `--te_goterms_file` overrides the hardcoded default TE-associated GO
  term list; `--print_te_associated_default_goterms` prints it.
  `--detenga_csv` overrides the default `DETENGA_run/{prefix}_TE_summary.csv`
  path. `GOTE_TAGGED=YES` requires every one of a protein's GO terms to be
  TE-associated (strict "only" rule); `NA` in either tag column means "no
  info available", not "checked and not a TE".
- `GAQET2AHRD.py`: `--skip_ahrd` no longer always skips the summary/TE
  check — if the AHRD output TSV already exists from an earlier run, the
  summary and `--check_te_goterms` still run against it (useful for
  regenerating just those without re-running AHRD).

## [v0.2.4] — 2026-08-25

### Fixed
- `GAQET2AHRD.py`: the summary always reported 0% of proteins with a
  description. AHRD's TSV also has an `Interpro-ID (Description)` column
  (present, and empty, even when InterPro isn't configured); the
  description-column detector matched both it and the real
  `Human-Readable-Description` column and kept the later (empty) match,
  so every description read as blank. Fixed by excluding any column whose
  header contains "interpro". Also skip `#`-prefixed comment lines (e.g.
  AHRD's leading `# AHRD-Version 3.11` line) explicitly.

## [v0.2.3] — 2026-08-25

### Added
- `GAQET2AHRD.py`: the functional-annotation summary is now also printed
  to stderr as ASCII tables at the end of a run (same metrics as
  `<prefix>_AHRD.summary.txt`), not just written to file.

## [v0.2.2] — 2026-08-25

### Added
- `GAQET2AHRD.py`: `<prefix>_AHRD.summary.txt` written after a successful
  AHRD run — counts of proteins with a description/GO term(s) vs. unknown,
  the AHRD-Quality-Code distribution, and the `--top_n` (default 10) most
  abundant descriptions. Disable with `--skip_summary`.

## [v0.2.1] — 2026-08-25

### Added
- `scripts/GFF3RenameGenes.py` — systematic gene/transcript/exon/CDS/UTR
  renaming from a fixed SeqID-based numbering scheme, with a pre-renaming
  structural sanity check.
- `scripts/GetFasta4EarlGreyGFF.py` — extract FASTA sequences for TE
  features from an EarlGrey repeat-annotation GFF3.
- `scripts/GFF2BEDOrthoVenn.py` — convert a GFF3 file to the 5-column BED
  format expected by OrthoVennPlus.
- `scripts/GAQET2AHRD.py` — parse a GAQET run's `GAQET.log.txt` for its
  TREMBL/SWISSPROT diamond commands, build an AHRD YAML config, and run
  AHRD.
- README: script groups (`FASTA_Utilities`, `GenomicData_Download`,
  `GFF_Utilities`, `ThirdPartyTool_Utilities`), a "Third-party tools and
  citations" section, and a Rationale subsection for each
  `ThirdPartyTool_Utilities` script.
- `NCBI_DownloadGenome.py`: `--strip_description` and `--report_metrics`
  (assembly metrics summary table), content-based SCF/CTG classification,
  CHR/organelle counts, ambiguous-nucleotide and duplicate-ID checks.
- `FastaRename.py`: `--sort` option.

### Fixed
- `GFF3RenameGenes.py`: handle multi-segment features sharing one ID (e.g.
  a CDS split across several exons) as a single logical feature during
  renaming, instead of assigning a separate ID per physical line.
- `GAQET2AHRD.py`: generated AHRD config was missing `gene_ontology_result`
  entirely, so AHRD transferred zero GO terms; now defaults it to
  `goa_uniprot_all.gaf` next to the SWISSPROT `--db` (override with
  `--gene_ontology_result`, disable with `--skip_go`).
- Citation for GAQET corrected to `github.com/victorgcb1987/GAQET2` (was
  incorrectly pointing at a different repository).

### Changed
- `CITATION.cff` version and abstract updated to reflect the full current
  script set.

## [v0.0.1] — 2026-06-24

### Added
- Initial release: `FastaRename.py`, `FastaStats.py`, `GFA2FASTA.py`,
  `NCBI_DownloadGenome.py`.
- `LICENSE` (MIT) and `CITATION.cff`.
