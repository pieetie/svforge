"""
Provenance tests: every VCF produced by svforge must carry the six
``##svforge*`` header tags, in the spec order, right after ``##fileformat``.
The warning tag must be the exact invariant string
``SYNTHETIC_DATA_DO_NOT_USE_FOR_CLINICAL_DIAGNOSIS`` and must not be
removable via CLI or otherwise. ``sanitize_command`` must strip absolute
paths to their basename on POSIX, macOS and Windows-style inputs
"""

from __future__ import annotations

from pathlib import Path

import pytest

from svforge import __version__
from svforge.cli import main
from svforge.core.provenance import (
    WARNING_VALUE,
    build_svforge_tags,
    sanitize_command,
)
from svforge.io.vcf_writer import write_vcf
from svforge.writers import get_writer

_EXPECTED_KEYS = (
    "##svforgeVersion",
    "##svforgeCommand",
    "##svforgeSeed",
    "##svforgeCaller",
    "##svforgeWarning",
    "##svforgeDocumentation",
)

def test_build_tags_returns_six_tags_in_order() -> None:
    tags = build_svforge_tags(caller="manta", seed=42, argv=["svforge", "gen"])
    assert len(tags) == 6
    keys = [t.split("=", 1)[0] for t in tags]
    assert tuple(keys) == _EXPECTED_KEYS

def test_warning_tag_value_is_the_invariant_string() -> None:
    tags = build_svforge_tags(caller="delly", seed=0, argv=["svforge"])
    warning = next(t for t in tags if t.startswith("##svforgeWarning="))
    assert warning == f"##svforgeWarning={WARNING_VALUE}"
    assert WARNING_VALUE == "SYNTHETIC_DATA_DO_NOT_USE_FOR_CLINICAL_DIAGNOSIS"

def test_version_tag_matches_package_version() -> None:
    tags = build_svforge_tags(caller="manta", seed=1, argv=["svforge"])
    assert f"##svforgeVersion={__version__}" in tags

def test_caller_tag_reflects_caller_argument() -> None:
    assert "##svforgeCaller=manta" in build_svforge_tags(caller="manta", seed=0, argv=["x"])
    assert "##svforgeCaller=delly" in build_svforge_tags(caller="delly", seed=0, argv=["x"])

@pytest.mark.parametrize(
    ("argv", "expected"),
    [
        (
            ["svforge", "gen", "--bank", "/mnt/beegfs/user/secret/bank.yaml"],
            "svforge gen --bank bank.yaml",
        ),
        (
            ["svforge", "gen", "--bank=/home/alice/banks/my_bank.yaml"],
            "svforge gen --bank=my_bank.yaml",
        ),
        (
            ["svforge", "--out", "C:\\Users\\bob\\run\\out.vcf"],
            "svforge --out out.vcf",
        ),
        (
            ["svforge", "--seed", "42"],
            "svforge --seed 42",
        ),
        (
            ["/Users/name with space/.venv/bin/svforge", "gen"],
            "svforge gen",
        ),
        (
            ["svforge", "--bank=/home/alice/my banks/bank.yaml"],
            "svforge --bank=bank.yaml",
        ),
        (
            ["svforge", "--out", "C:\\Users\\Name With Space\\svforge.exe"],
            "svforge --out svforge.exe",
        ),
        (
            ["svforge", "--out", "C:/Users/Name With Space/out.vcf"],
            "svforge --out out.vcf",
        ),
    ],
)
def test_sanitize_command_strips_absolute_paths(argv: list[str], expected: str) -> None:
    assert sanitize_command(argv) == expected

def test_sanitize_command_preserves_relative_paths() -> None:
    argv = ["svforge", "gen", "--bank", "banks/local.yaml", "--out", "./out.vcf"]
    assert sanitize_command(argv) == "svforge gen --bank banks/local.yaml --out ./out.vcf"

@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_header_lines_inject_tags_in_order_after_fileformat(caller: str) -> None:
    writer = get_writer(caller)
    tags = build_svforge_tags(caller=caller, seed=123, argv=["svforge", "gen"])
    header = writer.header_lines("S1", provenance_tags=tags)

    assert header[0].startswith("##fileformat=")
    injected = header[1 : 1 + len(_EXPECTED_KEYS)]
    assert [line.split("=", 1)[0] for line in injected] == list(_EXPECTED_KEYS)

@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_generated_vcf_contains_all_six_tags(tmp_path: Path, caller: str) -> None:
    writer = get_writer(caller)
    tags = build_svforge_tags(caller=caller, seed=123, argv=["svforge", "gen"])
    header = writer.header_lines("S1", provenance_tags=tags)

    out = tmp_path / f"{caller}.vcf"
    write_vcf(out, header, [])

    lines = out.read_text(encoding="utf-8").splitlines()
    assert lines[0].startswith("##fileformat=")
    injected = lines[1 : 1 + len(_EXPECTED_KEYS)]
    assert [line.split("=", 1)[0] for line in injected] == list(_EXPECTED_KEYS), (
        f"six ##svforge* tags must appear on lines 2-7 of the generated {caller} VCF; "
        f"found {injected!r}"
    )

def test_cli_gen_injects_provenance_and_warning(tmp_path: Path) -> None:
    """
    End-to-end: ``svforge gen`` writes a VCF with all six tags
    """
    bank = Path(__file__).parent / "fixtures" / "mini_bank.yaml"
    out = tmp_path / "out.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--bank",
            str(bank),
            "--out",
            str(out),
            "--n",
            "3",
            "--sample-name",
            "S1",
            "--seed",
            "42",
        ]
    )
    assert rc == 0
    text = out.read_text(encoding="utf-8")
    for key in _EXPECTED_KEYS:
        assert key in text, f"missing {key}"
    assert WARNING_VALUE in text
    assert "##svforgeSeed=42" in text
