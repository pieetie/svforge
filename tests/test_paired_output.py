"""
Tests for paired (somatic) VCF output structure.

Guard rail for issue #7 regression.
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from svforge.cli import main as cli_main


def _run_svforge(args: list[str]) -> None:
    rc = cli_main(args)
    if rc != 0:
        raise RuntimeError(f"svforge exited with {rc}")


@pytest.fixture
def paired_vcf(tmp_path: Path) -> Path:
    out = tmp_path / "somaticSV.vcf.gz"
    _run_svforge(
        [
            "gen-pair",
            "--caller",
            "manta",
            "--n-somatic",
            "10",
            "--n-germline",
            "5",
            "--tumor-sample-name",
            "MY_TUMOR",
            "--normal-sample-name",
            "MY_NORMAL",
            "--seed",
            "42",
            "--out",
            str(out),
        ]
    )
    return out


def test_paired_vcf_has_two_sample_columns(paired_vcf: Path) -> None:
    """Paired VCF must have NORMAL and TUMOR columns in Manta order."""
    with pysam.VariantFile(str(paired_vcf)) as vf:
        samples = list(vf.header.samples)
    assert samples == ["MY_NORMAL", "MY_TUMOR"], (
        f"Expected [NORMAL, TUMOR] columns, got {samples}"
    )


def test_paired_vcf_contains_both_somatic_and_germline(paired_vcf: Path) -> None:
    """Paired VCF must contain both somatic and germline records."""
    somatic_count = 0
    germline_count = 0
    with pysam.VariantFile(str(paired_vcf)) as vf:
        for rec in vf:
            if rec.info.get("SOMATIC"):
                somatic_count += 1
            else:
                germline_count += 1
    assert somatic_count > 0, "No somatic records found"
    assert germline_count > 0, "No germline records found"


def test_somatic_records_have_ref_only_normal_column(paired_vcf: Path) -> None:
    """For SOMATIC records, NORMAL column should have only ref support."""
    with pysam.VariantFile(str(paired_vcf)) as vf:
        for rec in vf:
            if rec.info.get("SOMATIC"):
                pr = rec.samples["MY_NORMAL"]["PR"]
                assert pr[1] == 0, (
                    f"Record {rec.id}: NORMAL PR alt should be 0 for "
                    f"somatic variant, got {pr}"
                )


def test_germline_records_have_alt_support_in_both_columns(paired_vcf: Path) -> None:
    """Germline records must have alt support in both NORMAL and TUMOR."""
    with pysam.VariantFile(str(paired_vcf)) as vf:
        for rec in vf:
            if not rec.info.get("SOMATIC"):
                normal_pr = rec.samples["MY_NORMAL"]["PR"]
                tumor_pr = rec.samples["MY_TUMOR"]["PR"]
                assert normal_pr[1] > 0, (
                    f"Record {rec.id}: germline should have alt support "
                    f"in NORMAL, got {normal_pr}"
                )
                assert tumor_pr[1] > 0, (
                    f"Record {rec.id}: germline should have alt support "
                    f"in TUMOR, got {tumor_pr}"
                )


def test_paired_vcf_is_tabix_indexable(paired_vcf: Path) -> None:
    """Paired VCF must be tabix-indexable (no #6 regression)."""
    pysam.tabix_index(str(paired_vcf), preset="vcf", force=True)
    assert (paired_vcf.parent / (paired_vcf.name + ".tbi")).exists()


@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_paired_works_for_both_callers(tmp_path: Path, caller: str) -> None:
    """Both Manta and DELLY support paired output with 2 sample columns."""
    out = tmp_path / f"{caller}_somatic.vcf"
    _run_svforge(
        [
            "gen-pair",
            "--caller",
            caller,
            "--n-somatic",
            "5",
            "--n-germline",
            "3",
            "--tumor-sample-name",
            "T",
            "--normal-sample-name",
            "N",
            "--seed",
            "42",
            "--out",
            str(out),
        ]
    )
    with pysam.VariantFile(str(out)) as vf:
        samples = list(vf.header.samples)
    assert len(samples) == 2


def test_caller_specific_column_order(tmp_path: Path) -> None:
    """Manta uses NORMAL,TUMOR order; DELLY uses TUMOR,NORMAL."""
    manta_out = tmp_path / "manta.vcf"
    _run_svforge(
        [
            "gen-pair",
            "--caller",
            "manta",
            "--n-somatic",
            "3",
            "--n-germline",
            "2",
            "--tumor-sample-name",
            "T",
            "--normal-sample-name",
            "N",
            "--seed",
            "42",
            "--out",
            str(manta_out),
        ]
    )
    with pysam.VariantFile(str(manta_out)) as vf:
        assert list(vf.header.samples) == ["N", "T"]

    delly_out = tmp_path / "delly.vcf"
    _run_svforge(
        [
            "gen-pair",
            "--caller",
            "delly",
            "--n-somatic",
            "3",
            "--n-germline",
            "2",
            "--tumor-sample-name",
            "T",
            "--normal-sample-name",
            "N",
            "--seed",
            "42",
            "--out",
            str(delly_out),
        ]
    )
    with pysam.VariantFile(str(delly_out)) as vf:
        assert list(vf.header.samples) == ["T", "N"]
