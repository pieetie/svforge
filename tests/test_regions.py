"""
Tests for BED region overlap index
"""

from __future__ import annotations

import gzip
from pathlib import Path

from svforge.core.regions import RegionSet, load_bed

def test_overlap_point() -> None:
    rs = RegionSet([("chr1", 100, 200)])
    assert rs.overlaps_point("chr1", 150)
    assert not rs.overlaps_point("chr1", 250)
    assert not rs.overlaps_point("chr2", 150)

def test_overlap_vcf_coords() -> None:
    rs = RegionSet([("chr1", 1000, 2000)])
    assert rs.overlaps_vcf("chr1", 1500, 1600)
    assert rs.overlaps_vcf("chr1", 500, 1500)
    assert not rs.overlaps_vcf("chr1", 2001, 3000)

def test_overlap_bed_half_open() -> None:
    rs = RegionSet([("chr1", 100, 200)])
    assert rs.overlaps_bed("chr1", 100, 101)
    assert rs.overlaps_bed("chr1", 199, 200)
    assert not rs.overlaps_bed("chr1", 200, 300)

def test_empty_regionset() -> None:
    rs = RegionSet()
    assert len(rs) == 0
    assert not rs.overlaps_point("chr1", 1)

def test_load_bed_plain(mini_blacklist_path: Path) -> None:
    rs = load_bed(mini_blacklist_path)
    assert len(rs) >= 8
    assert "chr1" in rs.chromosomes
    assert rs.overlaps_bed("chr1", 1_050_000, 1_060_000)

def test_load_bed_gzipped(tmp_path: Path, mini_blacklist_path: Path) -> None:
    src = mini_blacklist_path.read_bytes()
    gz = tmp_path / "bl.bed.gz"
    gz.write_bytes(gzip.compress(src))
    rs = load_bed(gz)
    assert len(rs) >= 8
    assert rs.overlaps_bed("chr2", 2_010_000, 2_020_000)
