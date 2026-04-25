"""
Tests for the bundled mini injection catalogs
"""

from __future__ import annotations

from svforge.core.genome import HG38_CONTIGS
from svforge.core.injection_catalogs import load_blacklist_catalog, load_gnomad_catalog


def test_gnomad_catalog_parses_with_min_entries() -> None:
    catalog = load_gnomad_catalog()
    assert len(catalog) >= 50

def test_blacklist_catalog_parses_with_min_entries() -> None:
    catalog = load_blacklist_catalog()
    assert len(catalog) >= 50

def test_gnomad_chroms_and_positions_within_contigs() -> None:
    for entry in load_gnomad_catalog():
        assert entry.chrom in HG38_CONTIGS, entry
        assert entry.end_chrom in HG38_CONTIGS, entry
        assert 1 <= entry.pos <= HG38_CONTIGS[entry.chrom]
        assert 1 <= entry.end <= HG38_CONTIGS[entry.end_chrom]

def test_blacklist_chroms_and_positions_within_contigs() -> None:
    for entry in load_blacklist_catalog():
        assert entry.chrom in HG38_CONTIGS, entry
        assert 1 <= entry.pos <= HG38_CONTIGS[entry.chrom]
        assert 1 <= entry.end <= HG38_CONTIGS[entry.chrom]
        assert entry.end >= entry.pos

def test_gnomad_has_four_or_more_svtypes() -> None:
    svtypes = {entry.svtype for entry in load_gnomad_catalog()}
    assert len(svtypes) >= 4, svtypes

def test_catalogs_cover_at_least_10_chroms_combined() -> None:
    chroms = {e.chrom for e in load_gnomad_catalog()} | {e.chrom for e in load_blacklist_catalog()}
    assert len(chroms) >= 10, chroms

def test_gnomad_source_ids_unique() -> None:
    ids = [e.source_id for e in load_gnomad_catalog()]
    assert len(ids) == len(set(ids))

def test_blacklist_source_ids_unique() -> None:
    ids = [e.source_id for e in load_blacklist_catalog()]
    assert len(ids) == len(set(ids))
