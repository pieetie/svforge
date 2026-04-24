# svForge - ready-to-use commands

Test **outputs** in the examples below are written under `data_local/gen-test/` 

```bash
mkdir -p data_local/gen-test
```

---

## 1. Environment and install

```bash
python -m venv .venv
source .venv/bin/activate          # Windows: .venv\Scripts\activate
pip install -U pip
pip install -e ".[dev,test]"       # app + ruff, mypy, pytest
```

```bash
svforge --version
```

**More logging** (optional):

```bash
svforge -v gen --help
svforge -vv gen --help
```

---

## 2. Discover built-in options

| Command | What it does |
|--------|----------------|
| `svforge callers` | Prints every **writer** name for `--caller` (e.g. `manta`, `delly`). |
| `svforge bank list` | Lists **built-in bank** names. |
| `svforge bank show` | Dumps the default bank (`default_hg38`) as YAML. |
| `svforge bank show <name>` | Dumps a specific built-in bank. |

```bash
svforge callers
svforge bank list
svforge bank show
svforge bank show default_hg38
```

---

## 3. Generate one VCF (`gen`)

**Minimal example** (Manta-style VCF, 50 records, one sample):

```bash
svforge gen --caller manta --out data_local/gen-test/synthetic.vcf.gz --n 50 --sample-name TUMOR01
```

**Same result every time** (fixed seed):

```bash
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 100 --sample-name S1 --seed 42
```

**DELLY** (output can be `.bcf` if you set the path):

```bash
svforge gen --caller delly --out data_local/gen-test/synthetic.bcf --n 30 --sample-name S1 --seed 1
```

**Only some SV types** (comma-separated, upper case) :

```bash
svforge gen --caller manta --out data_local/gen-test/only_del_dup.vcf.gz --n 20 --sample-name S1 \
  --svtypes DEL,DUP
```

**Custom bank (YAML on disk):**

```bash
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 20 --sample-name S1 --bank /path/to/bank.yaml
```

**Custom VCF header template:**

```bash
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 10 --sample-name S1 \
  --header-template /path/to/header.template
```

**Reference build** (v1 is mainly **hg38**) :

```bash
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 10 --sample-name S1 --genome hg38
```

**Sampling ranges** (examples):

```bash
# SV length in bp, inclusive range
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 15 --sample-name S1 --svlen-range 100,5000

# Microhomology length
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 15 --sample-name S1 --homlen-range 5,30

# Variant allele frequency range (0–1)
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 15 --sample-name S1 --vaf-range 0.2,0.9
```

**Only some chromosomes** (with or without `chr`) :

```bash
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 20 --sample-name S1 --chromosomes chr1,chr7
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 20 --sample-name S1 --chromosomes 1,7
```

**Bundled mini-catalogs** (fractions in `[0, 1]`; the rest is synthetic) :

```bash
# gnomAD-style mini TSV
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 50 --sample-name S1 --gnomad-fraction 0.2

# ENCODE blacklist-style mini TSV
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 50 --sample-name S1 --blacklist-fraction 0.1

# Both
svforge gen --caller manta --out data_local/gen-test/out.vcf.gz --n 100 --sample-name S1 --gnomad-fraction 0.15 --blacklist-fraction 0.1 --seed 7
```

---

## 4. Tumor + normal pair (`gen-pair`)

**Example:**

```bash
svforge gen-pair \
  --caller manta \
  --out-tumor data_local/gen-test/tumor.vcf.gz \
  --out-normal data_local/gen-test/normal.vcf.gz \
  --n-somatic 30 \
  --n-germline 10 \
  --tumor-sample-name TUMOR_01 \
  --normal-sample-name NORMAL_01 \
  --seed 99
```

Same **sampling flags** as `gen` work here (`--bank`, `--gnomad-fraction`, `--blacklist-fraction`, `--svtypes`, ranges, etc.).

---

## 5. Check a VCF you generated (`validate`)

Exit code `0` = pass; non-zero = fail.

**Print to the terminal:**

```bash
svforge validate --vcf data_local/gen-test/synthetic.vcf.gz
```

**Write a TSV report:**

```bash
svforge validate --vcf data_local/gen-test/synthetic.vcf.gz --report-tsv data_local/gen-test/validate_report.tsv
```

---

## 6. Developer commands (from repository root)

| Step | Command |
|------|--------|
| Lint | `ruff check .` |
| Auto-fix (safe) | `ruff check . --fix` |
| Format | `ruff format .` |
| Types | `mypy src` |
| Tests | `pytest -q` |
| Tests + coverage | `pytest --cov=svforge --cov-report=term-missing` |
| One file | `pytest tests/test_cli.py` |

**All in one go:**

```bash
ruff check . && ruff format --check . && mypy src && pytest -q
```

---

## 7. Optional: **bcftools** (sanity on the file)

```bash
bcftools view -h data_local/gen-test/synthetic.vcf.gz
bcftools view -H data_local/gen-test/synthetic.vcf.gz | head
```

---

## 8. Help in the terminal

```bash
svforge --help
svforge gen --help
svforge gen-pair --help
svforge validate --help
svforge bank --help
```
