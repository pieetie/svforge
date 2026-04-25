# Contributing to svForge

## Development setup

```bash
python -m venv .venv
source .venv/bin/activate
pip install -e ".[dev,test]"
```

## Contribution workflow

1. Open an issue first and describe the change (new caller, new template, fix, etc.). Wait for feedback before coding.
2. Create a dedicated branch from `main` named `{issue-number}-short-description`. Never work directly on `main`.
3. Use commit messages that reference the issue number: `feat: XXX #12`, `fix: XXX (#15)`.
4. Update `CHANGELOG.md`.
5. Run the local quality gate before pushing. If it is not all green, the PR is not reviewed:

```bash
ruff check .
ruff format --check .
mypy src
pytest -q
```

6. Push and open a PR to `main`. In the PR description, add `Closes #<issue-number>` for automatic issue closing at merge.
7. Use squash merge. The squash commit message must be clean and final (not a draft branch history dump).
8. Delete the branch after merge.

## Common contribution cases

- Add a new caller (GRIDSS, SvABA, Lumpy, Smoove, ...). This is the main contribution path and a core modularity goal of the project.
- Add a new header template for an uncovered genome or scenario (for example `hg19`, germline mode, tumor-only mode).
- Enrich the bundled mini-catalogs (`gnomad_hg38_mini.tsv`, `blacklist_hg38_mini.tsv`) with more entries or more diversity.
- Improve the default bank (`default_hg38.yaml`) with more representative SV templates.
- Report a bug with a minimal command that reproduces it.
- Improve documentation.
- Add tests for an uncovered edge case (for example negative `SVLEN`, BND to alt contig, etc.).
