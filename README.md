<h1>
<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/pieetie/svforge/main/docs/assets/vector-svg/logos-no-background/logo-red-and-white.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/pieetie/svforge/main/docs/assets/vector-svg/logos-no-background/logo-red-and-black.svg">
  <img src="https://raw.githubusercontent.com/pieetie/svforge/main/docs/assets/vector-svg/logos-no-background/logo-red-and-black.svg" alt="svForge logo" width="450">
</picture>
</h1>

### Generate synthetic SV VCFs to stress-test your pipelines with confidence

[![PyPI version](https://img.shields.io/pypi/v/svforge.svg)](https://pypi.org/project/svforge/)
[![License](https://img.shields.io/pypi/l/svforge)](https://pypi.org/project/svforge/)
[![DOI](https://zenodo.org/badge/1218168425.svg)](https://doi.org/10.5281/zenodo.19762333)

---

**svForge** produces caller-specific VCFs (Manta, DELLY) in VCF / VCF.gz / BCF format with fine-grained control over variability (HOMLEN, SVLEN, VAF) and realistic artefact injection (SVs in ENCODE blacklist regions, gnomAD germline SVs).

Designed to be modular, it is easy to adapt to your own use case. 
You can tune generation parameters, plug in new callers, and customize the workflow without reworking the whole tool.

## Installation

```bash
pip install svforge
```

Or from source:

```bash
git clone https://github.com/pieetie/svforge
cd svforge
pip install -e ".[dev,test]"
```

## Quick start

For ready-to-run command lines (single sample, tumor/normal pair, validation, banks, and dev checks), see [`docs/ready-to-use.md`](docs/ready-to-use.md).

## Typical use cases

- Validate downstream filters (for example, `SVFORGE_SOURCE=gnomad` records should disappear after your gnomAD filtering step).
- Validate ENCODE blacklist annotation logic (for example, `SVFORGE_SOURCE=blacklist` records should receive your expected poor-mappability label).
- Run reproducible CI regression tests with fixed seeds, without committing generated VCF files.
- Smoke-test deployments and pipeline changes in seconds instead of rerunning full variant callers on BAM files.
- Reproduce specific scenarios and edge cases (cross-chrom BNDs, contig-edge events, controlled VAF/HOMLEN ranges) for debugging and QA.
- Demo or onboard safely with realistic SV VCFs and no patient data.

## CLI

```
svforge gen          # one VCF for one sample
svforge gen-pair     # tumor + normal VCFs for somatic pipelines
svforge validate     # self-consistency check of injected SVs
svforge bank list    # list built-in banks
svforge bank show    # dump a bank as YAML
svforge callers      # list registered writers
```

Run `svforge <cmd> --help` for the full flag list.

## Credits

- **Logo and visual identity** - [Elisa Perrin](https://www.linkedin.com/in/elisaperrin/)
- **Claude** (Anthropic) - assisted with tests, documentation, refactoring, and release tooling (Ruff linting/formatting, CI cleanup)

## License

<picture>
  <source media="(prefers-color-scheme: dark)" srcset="https://raw.githubusercontent.com/pieetie/svforge/main/docs/assets/vector-svg/icons-on-background/icon-nb-app.svg">
  <source media="(prefers-color-scheme: light)" srcset="https://raw.githubusercontent.com/pieetie/svforge/main/docs/assets/vector-svg/icons-on-background/icon-bn-app.svg">
  <img src="https://raw.githubusercontent.com/pieetie/svforge/main/docs/assets/vector-svg/icons-on-background/icon-bn-app.svg" alt="svForge icon" width="14" style="vertical-align: -1px;">
</picture> Distributed under the <a href="./LICENSE">MIT License</a>.
