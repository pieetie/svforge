"""
Tests for the sampler
"""

from __future__ import annotations

from svforge.core.bank import Bank
from svforge.core.sampler import SamplerConfig, sample, sample_pair

def test_sample_deterministic_with_seed(mini_bank: Bank) -> None:
    cfg = SamplerConfig(n=30, seed=42)
    a = sample(mini_bank, cfg)
    b = sample(mini_bank, cfg)
    assert [sv.id for sv in a] == [sv.id for sv in b]
    assert [(sv.chrom, sv.pos, sv.svtype, sv.svlen) for sv in a] == [
        (sv.chrom, sv.pos, sv.svtype, sv.svlen) for sv in b
    ]

def test_sample_different_seeds_differ(mini_bank: Bank) -> None:
    a = sample(mini_bank, SamplerConfig(n=20, seed=1))
    b = sample(mini_bank, SamplerConfig(n=20, seed=2))
    assert [sv.id for sv in a] != [sv.id for sv in b]

def test_sample_count(mini_bank: Bank) -> None:
    svs = sample(mini_bank, SamplerConfig(n=50, seed=0))
    assert len(svs) == 50

def test_sample_restricted_svtypes(mini_bank: Bank) -> None:
    svs = sample(mini_bank, SamplerConfig(n=40, seed=0, svtypes=frozenset({"DEL"})))
    assert {sv.svtype for sv in svs} == {"DEL"}

def test_sample_svlen_range_clamp(mini_bank: Bank) -> None:
    svs = sample(
        mini_bank,
        SamplerConfig(n=40, seed=0, svtypes=frozenset({"DEL"}), svlen_range=(200, 300)),
    )
    for sv in svs:
        assert 200 <= sv.svlen <= 300

def test_sample_vaf_range(mini_bank: Bank) -> None:
    svs = sample(mini_bank, SamplerConfig(n=30, seed=0, vaf_range=(0.2, 0.8)))
    for sv in svs:
        assert 0.2 <= sv.vaf <= 0.8

def test_sample_blacklist_injection(mini_bank: Bank) -> None:
    cfg = SamplerConfig(
        n=20,
        seed=0,
        svtypes=frozenset({"DEL", "DUP", "INV"}),
        blacklist_fraction=0.4,
    )
    svs = sample(mini_bank, cfg)
    assert sum(1 for sv in svs if sv.source == "blacklist") == 8

def test_sample_pair_germline_shared(mini_bank: Bank) -> None:
    pair = sample_pair(
        mini_bank,
        n_somatic=10,
        n_germline=15,
        cfg=SamplerConfig(n=0, seed=7, svtypes=frozenset({"DEL", "DUP", "INV"})),
    )
    germline_ids = pair.germline_ids
    somatic_ids = pair.somatic_ids
    assert len(germline_ids) == 15
    assert len(somatic_ids) == 10
    assert germline_ids.isdisjoint(somatic_ids)

    normal_coords = {(sv.id, sv.chrom, sv.pos, sv.end) for sv in pair.normal}
    tumor_germ = {(sv.id, sv.chrom, sv.pos, sv.end) for sv in pair.tumor if sv.id in germline_ids}
    assert normal_coords == tumor_germ
