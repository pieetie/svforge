# Changelog

## [1.0.1] — 2026-04-28

### Fixed

- `gen-pair` now produces a two-sample somatic VCF (#7).
- VCF output is now tabix-indexable; previously inter-chromosomal BND records broke sort order (#6).
- README DOI and license badges (#8, #9).

### Changed

- **BREAKING**: `gen-pair` replaces `--out-tumor` / `--out-normal` with a single `--out` flag.

## [1.0.0] — 2026-04-25

Initial release.

- `svforge gen` and `svforge gen-pair` to generate synthetic VCFs (single-sample or coherent tumor/normal pair)
- Manta (VCFv4.1) and DELLY (VCFv4.2) writers
- Injection via `--gnomad-fraction` and `--blacklist-fraction`
- `INFO/SVFORGE_SOURCE` tag on every record (`bank` / `gnomad` / `blacklist`) for self-verifying pipelines
- `svforge validate` to confirm injected records match the bundled catalogs exactly
- Chromosome filtering via `--chromosomes`
- Custom header support via `--header-template PATH`
- Deterministic sampling via `--seed`; effective seed always logged in the output VCF
- BCF / VCF.gz / VCF output formats
- Default hg38 SV bank with weighted templates across DEL, DUP, INV, INS, BND
- Writer plugin system: third-party callers can register via `svforge.writers` entry point
- Non-configurable `##svforgeWarning=SYNTHETIC_DATA_DO_NOT_USE_FOR_CLINICAL_DIAGNOSIS` injected in every VCF header
- `sanitize_command()` strips absolute paths from the logged command line so user home directories and cluster paths never leak into generated VCFs
- Python 3.10+, hg38 only
