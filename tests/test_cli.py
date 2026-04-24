"""
CLI integration tests -- each subcommand end-to-end
"""

from __future__ import annotations

from pathlib import Path

import pysam

from svforge.cli import main

def test_cli_callers(capsys) -> None:  # type: ignore[no-untyped-def]
    rc = main(["callers"])
    assert rc == 0
    out = capsys.readouterr().out.splitlines()
    assert "manta" in out
    assert "delly" in out

def test_cli_bank_list(capsys) -> None:  # type: ignore[no-untyped-def]
    rc = main(["bank", "list"])
    assert rc == 0
    assert "default_hg38" in capsys.readouterr().out

def test_cli_bank_show(capsys) -> None:  # type: ignore[no-untyped-def]
    rc = main(["bank", "show", "default_hg38"])
    assert rc == 0
    out = capsys.readouterr().out
    assert "name: default_hg38" in out
    assert "templates:" in out

def test_cli_gen_manta_vcf_gz(tmp_path: Path, mini_bank_path: Path) -> None:
    out = tmp_path / "out.vcf.gz"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--out",
            str(out),
            "--n",
            "15",
            "--sample-name",
            "T01",
            "--bank",
            str(mini_bank_path),
            "--seed",
            "42",
        ]
    )
    assert rc == 0
    assert out.exists()
    with pysam.VariantFile(str(out)) as vf:
        assert sum(1 for _ in vf) > 0

def test_cli_gen_pair(tmp_path: Path, mini_bank_path: Path) -> None:
    tumor = tmp_path / "tumor.vcf.gz"
    normal = tmp_path / "normal.vcf.gz"
    rc = main(
        [
            "gen-pair",
            "--caller",
            "delly",
            "--out-tumor",
            str(tumor),
            "--out-normal",
            str(normal),
            "--n-somatic",
            "5",
            "--n-germline",
            "6",
            "--tumor-sample-name",
            "T01",
            "--normal-sample-name",
            "N01",
            "--bank",
            str(mini_bank_path),
            "--seed",
            "99",
            "--svtypes",
            "DEL,DUP,INV",
        ]
    )
    assert rc == 0
    with pysam.VariantFile(str(tumor)) as vf:
        tumor_ids = {r.id for r in vf if r.id}
    with pysam.VariantFile(str(normal)) as vf:
        normal_ids = {r.id for r in vf if r.id}
    assert normal_ids.issubset(tumor_ids)
    assert len(tumor_ids - normal_ids) >= 5

def test_cli_validate_pass(
    tmp_path: Path,
    mini_bank_path: Path,
) -> None:
    gen_out = tmp_path / "gen.vcf.gz"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--out",
            str(gen_out),
            "--n",
            "40",
            "--sample-name",
            "T01",
            "--bank",
            str(mini_bank_path),
            "--seed",
            "7",
            "--svtypes",
            "DEL,DUP,INV",
            "--blacklist-fraction",
            "0.10",
            "--gnomad-fraction",
            "0.20",
        ]
    )
    assert rc == 0

    report = tmp_path / "report.tsv"
    rc = main(
        [
            "validate",
            "--vcf",
            str(gen_out),
            "--report-tsv",
            str(report),
        ]
    )
    assert rc == 0
    content = report.read_text()
    assert "metric" in content
    assert "gnomad_self_consistency\t" in content
    assert "blacklist_self_consistency\t" in content
    assert content.count("PASS") == 2
    assert "FAIL" not in content
