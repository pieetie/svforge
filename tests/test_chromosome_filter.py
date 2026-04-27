"""
Chromosome-filter tests for ``--chromosomes`` on gen and gen-pair

Requirements from the V1 robustness pass:

1. ``--chromosomes chr1`` only produces records on ``chr1`` and no
   cross-chromosome BND
2. ``--chromosomes chr8,chr14`` permits chr8<->chr14 BND events
3. ``--chromosomes 1,7`` is normalised to ``chr1,chr7``
4. ``--chromosomes chrZZ`` raises a validation error with a non-zero exit code
5. A filter that empties the pool raises an explicit error
"""

from __future__ import annotations

import re
from pathlib import Path

import pysam

from svforge.cli import main
from svforge.core.genome import normalize_chromosomes

FIXTURES = Path(__file__).parent / "fixtures"
BANK = FIXTURES / "mini_bank.yaml"

def _records(vcf: Path) -> list[pysam.VariantRecord]:
    with pysam.VariantFile(str(vcf)) as vf:
        return list(vf)

def test_single_chromosome_filters_all_records(tmp_path: Path) -> None:
    out = tmp_path / "out.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(out),
            "--n",
            "20",
            "--sample-name",
            "S",
            "--chromosomes",
            "chr1",
            "--seed",
            "7",
            "--svtypes",
            "DEL,DUP,INV,INS",
        ]
    )
    assert rc == 0
    recs = _records(out)
    assert recs
    for rec in recs:
        assert rec.chrom == "chr1"

def test_two_chromosomes_permit_bnd_between_them(tmp_path: Path) -> None:
    out = tmp_path / "out.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "delly",
            "--bank",
            str(BANK),
            "--out",
            str(out),
            "--n",
            "30",
            "--sample-name",
            "S",
            "--chromosomes",
            "chr8,chr14",
            "--seed",
            "11",
        ]
    )
    assert rc == 0
    recs = _records(out)
    assert recs
    allowed = {"chr8", "chr14"}
    for rec in recs:
        assert rec.chrom in allowed
        if rec.info.get("SVTYPE") == "BND":
            alt = str(rec.alts[0]) if rec.alts else ""
            mate_match = re.search(r"[\[\]]([^:\[\]]+):\d+[\[\]]", alt)
            if mate_match:
                assert mate_match.group(1) in allowed

def test_alias_digits_normalize_to_chr_prefix() -> None:
    assert normalize_chromosomes(["1", "7"], "hg38") == ["chr1", "chr7"]
    assert normalize_chromosomes(["X", "MT"], "hg38") == ["chrX", "chrM"]
    assert normalize_chromosomes(["chr3"], "hg38") == ["chr3"]

def test_cli_alias_digits_normalize(tmp_path: Path) -> None:
    out = tmp_path / "out.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(out),
            "--n",
            "5",
            "--sample-name",
            "S",
            "--chromosomes",
            "1",
            "--seed",
            "0",
            "--svtypes",
            "DEL",
        ]
    )
    assert rc == 0
    for rec in _records(out):
        assert rec.chrom == "chr1"

def test_unknown_chromosome_exits_non_zero(tmp_path: Path) -> None:
    out = tmp_path / "out.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(out),
            "--n",
            "1",
            "--sample-name",
            "S",
            "--chromosomes",
            "chrZZ",
            "--seed",
            "0",
        ]
    )
    assert rc != 0

def test_empty_pool_raises_explicit_error(tmp_path: Path) -> None:
    """
    Filter requests a chromosome that no template in the bank covers
    """
    constrained = tmp_path / "bank.yaml"
    constrained.write_text(
        "name: bounded\n"
        "genome: hg38\n"
        "templates:\n"
        "  - svtype: DEL\n"
        "    svlen: [100, 1000]\n"
        "    homlen: [0, 5]\n"
        "    weight: 1.0\n"
        "    chroms: [chr1, chr2]\n",
        encoding="utf-8",
    )
    out = tmp_path / "out.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(constrained),
            "--out",
            str(out),
            "--n",
            "5",
            "--sample-name",
            "S",
            "--chromosomes",
            "chr21,chr22",
            "--seed",
            "0",
        ]
    )
    assert rc != 0

def test_gen_pair_accepts_chromosome_filter(tmp_path: Path) -> None:
    paired = tmp_path / "paired.vcf"
    rc = main(
        [
            "gen-pair",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(paired),
            "--n-somatic",
            "5",
            "--n-germline",
            "5",
            "--tumor-sample-name",
            "T",
            "--normal-sample-name",
            "N",
            "--chromosomes",
            "chr1,chr2",
            "--seed",
            "3",
            "--svtypes",
            "DEL,DUP",
        ]
    )
    assert rc == 0
    for rec in _records(paired):
        assert rec.chrom in {"chr1", "chr2"}
