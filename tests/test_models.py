"""
Tests for core domain models
"""

from __future__ import annotations

import pytest

from svforge.core.models import SV, Breakpoint, SVPair


class TestBreakpoint:
    def test_valid(self) -> None:
        bp = Breakpoint("chr1", 100, "+")
        assert bp.chrom == "chr1"
        assert bp.pos == 100
        assert bp.strand == "+"

    def test_invalid_pos(self) -> None:
        with pytest.raises(ValueError, match="pos"):
            Breakpoint("chr1", 0)

    def test_invalid_strand(self) -> None:
        with pytest.raises(ValueError, match="strand"):
            Breakpoint("chr1", 100, "x")  #type : ignore[arg-type]

class TestSV:
    def test_valid_deletion(self) -> None:
        sv = SV(id="x", svtype="DEL", chrom="chr1", pos=100, end=200, svlen=100)
        assert sv.svlen == 100
        assert sv.breakpoint1 == Breakpoint("chr1", 100, "+")
        assert sv.breakpoint2 == Breakpoint("chr1", 200, "-")

    def test_valid_bnd(self) -> None:
        sv = SV(
            id="x",
            svtype="BND",
            chrom="chr1",
            pos=100,
            end=100,
            svlen=1,
            mate_chrom="chr2",
            mate_pos=200,
            strands="+-",
        )
        assert sv.breakpoint2.chrom == "chr2"
        assert sv.breakpoint2.pos == 200

    def test_unknown_svtype(self) -> None:
        with pytest.raises(ValueError, match="svtype"):
            SV(id="x", svtype="WAT", chrom="chr1", pos=1, end=2, svlen=1)

    def test_bnd_requires_mate(self) -> None:
        with pytest.raises(ValueError, match="mate"):
            SV(id="x", svtype="BND", chrom="chr1", pos=1, end=1, svlen=1)

    def test_end_before_pos(self) -> None:
        with pytest.raises(ValueError, match="end"):
            SV(id="x", svtype="DEL", chrom="chr1", pos=200, end=100, svlen=1)

    def test_invalid_vaf(self) -> None:
        with pytest.raises(ValueError, match="VAF"):
            SV(id="x", svtype="DEL", chrom="chr1", pos=1, end=2, svlen=1, vaf=1.5)

    def test_invalid_strands(self) -> None:
        with pytest.raises(ValueError, match="strands"):
            SV(id="x", svtype="DEL", chrom="chr1", pos=1, end=2, svlen=1, strands="xy")

class TestSVPair:
    def _mk(self, sv_id: str, pos: int) -> SV:
        return SV(id=sv_id, svtype="DEL", chrom="chr1", pos=pos, end=pos + 50, svlen=50)

    def test_somatic_and_germline_partition(self) -> None:
        germ = self._mk("g1", 1_000)
        som = self._mk("s1", 2_000)
        pair = SVPair(tumor=[germ, som], normal=[germ])
        assert pair.germline_ids == {"g1"}
        assert pair.somatic_ids == {"s1"}

    def test_empty(self) -> None:
        pair = SVPair()
        assert pair.somatic_ids == set()
        assert pair.germline_ids == set()
