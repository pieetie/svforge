"""
Shared fixtures for the svforge test suite
"""

from __future__ import annotations

from pathlib import Path

import pytest

from svforge.core.bank import Bank, load_bank
from svforge.core.models import SV
from svforge.core.sampler import SamplerConfig, sample

FIXTURES_DIR = Path(__file__).parent / "fixtures"

@pytest.fixture
def fixtures_dir() -> Path:
    return FIXTURES_DIR

@pytest.fixture
def mini_bank_path(fixtures_dir: Path) -> Path:
    return fixtures_dir / "mini_bank.yaml"

@pytest.fixture
def mini_blacklist_path(fixtures_dir: Path) -> Path:
    return fixtures_dir / "mini_blacklist.bed"

@pytest.fixture
def mini_bank(mini_bank_path: Path) -> Bank:
    return load_bank(mini_bank_path)

@pytest.fixture
def sample_svs(mini_bank: Bank) -> list[SV]:
    cfg = SamplerConfig(n=10, seed=42)
    return sample(mini_bank, cfg)
