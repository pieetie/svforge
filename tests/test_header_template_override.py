"""
User-supplied header template override: ``--header-template PATH``

1. A valid template with all required placeholders is used as-is
2. A template missing a required placeholder raises an explicit error
3. Omitting ``--header-template`` keeps the bundled template active
"""

from __future__ import annotations

from pathlib import Path

import pytest

from svforge.cli import main
from svforge.writers import get_writer

FIXTURES = Path(__file__).parent / "fixtures"
BANK = FIXTURES / "mini_bank.yaml"
CUSTOM = FIXTURES / "custom_manta.template"
BROKEN = FIXTURES / "broken_manta.template"

def test_override_template_is_used(tmp_path: Path) -> None:
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
            "--seed",
            "0",
            "--svtypes",
            "DEL",
            "--chromosomes",
            "chr1",
            "--header-template",
            str(CUSTOM),
        ]
    )
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    assert "CUSTOM_SVFORGE_TAG" in text, "custom template marker not propagated"

def test_broken_template_raises_missing_placeholder(tmp_path: Path) -> None:
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
            "--seed",
            "0",
            "--header-template",
            str(BROKEN),
        ]
    )
    assert rc != 0

def test_missing_placeholder_raised_via_writer_api() -> None:
    writer = get_writer("manta")
    with pytest.raises(ValueError, match=r"missing required placeholders.*\{REFERENCE_FASTA\}"):
        writer.header_lines("S", template_override=BROKEN)

def test_bundled_template_used_when_override_absent(tmp_path: Path) -> None:
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
            "--seed",
            "0",
            "--svtypes",
            "DEL",
            "--chromosomes",
            "chr1",
        ]
    )
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    assert "CUSTOM_SVFORGE_TAG" not in text
    # the bundled Manta template advertises GenerateSVCandidates
    assert "GenerateSVCandidates" in text

def test_nonexistent_template_raises(tmp_path: Path) -> None:
    writer = get_writer("manta")
    with pytest.raises(FileNotFoundError):
        writer.header_lines("S", template_override=tmp_path / "nope.template")
