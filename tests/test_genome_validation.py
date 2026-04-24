"""
Validate the ``--genome`` CLI option

V1 supports hg38 only. Any other value must exit non-zero with a clear
message that points the user at ``--header-template`` for self-service
override
"""

from __future__ import annotations

from pathlib import Path

import pytest

from svforge.cli import main
from svforge.core.genome import validate_genome

BANK = Path(__file__).parent / "fixtures" / "mini_bank.yaml"

def test_validate_genome_accepts_hg38() -> None:
    assert validate_genome("hg38") == "hg38"

@pytest.mark.parametrize("bad", ["hg19", "t2t", "chm13", "invalid"])
def test_validate_genome_rejects_everything_else(bad: str) -> None:
    with pytest.raises(ValueError, match="No bundled header template"):
        validate_genome(bad)

def test_validate_genome_error_points_to_header_template() -> None:
    with pytest.raises(ValueError, match="--header-template PATH"):
        validate_genome("hg19")

def test_cli_gen_with_hg38_succeeds(tmp_path: Path) -> None:
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(tmp_path / "out.vcf"),
            "--n",
            "2",
            "--sample-name",
            "S",
            "--genome",
            "hg38",
            "--seed",
            "1",
            "--svtypes",
            "DEL",
        ]
    )
    assert rc == 0

def test_cli_gen_with_hg19_exits_non_zero(tmp_path: Path) -> None:
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(tmp_path / "out.vcf"),
            "--n",
            "1",
            "--sample-name",
            "S",
            "--genome",
            "hg19",
            "--seed",
            "0",
        ]
    )
    assert rc != 0

def test_cli_gen_with_invalid_genome_exits_non_zero(tmp_path: Path) -> None:
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(BANK),
            "--out",
            str(tmp_path / "out.vcf"),
            "--n",
            "1",
            "--sample-name",
            "S",
            "--genome",
            "does_not_exist",
            "--seed",
            "0",
        ]
    )
    assert rc != 0
