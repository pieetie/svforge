"""
Tests for bank loading and the built-in default bank
"""

from __future__ import annotations

from pathlib import Path

import pytest

from svforge.core.bank import list_builtin_banks, load_bank

def test_load_builtin_default_bank() -> None:
    bank = load_bank("default_hg38")
    assert bank.name == "default_hg38"
    assert bank.genome == "hg38"
    assert len(bank.templates) >= 15
    assert bank.svtypes() == {"DEL", "DUP", "INV", "INS", "BND"}

def test_list_builtin_banks_contains_default() -> None:
    assert "default_hg38" in list_builtin_banks()

def test_load_mini_bank_from_path(mini_bank_path: Path) -> None:
    bank = load_bank(mini_bank_path)
    assert bank.name == "mini_bank"
    assert len(bank.templates) == 5

def test_load_bank_missing() -> None:
    with pytest.raises(FileNotFoundError):
        load_bank("does_not_exist_anywhere")

def test_invalid_bank_missing_templates(tmp_path: Path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: bad\ngenome: hg38\ntemplates: []\n")
    with pytest.raises(ValueError, match="templates"):
        load_bank(bad)

def test_invalid_genome(tmp_path: Path) -> None:
    bad = tmp_path / "bad.yaml"
    bad.write_text("name: bad\ngenome: t2t\ntemplates:\n  - {svtype: DEL, svlen: [100, 200]}\n")
    with pytest.raises(ValueError, match="genome"):
        load_bank(bad)
