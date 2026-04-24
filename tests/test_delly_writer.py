"""
Tests for the DELLY writer with pysam round-trip validation
"""

from __future__ import annotations

from pathlib import Path

import pysam

from svforge.core.models import SV
from svforge.io.vcf_writer import write_vcf
from svforge.writers import get_writer

def _svs() -> list[SV]:
    return [
        SV(id="del1", svtype="DEL", chrom="chr1", pos=100_000, end=101_000, svlen=1_000),
        SV(id="dup1", svtype="DUP", chrom="chr2", pos=200_000, end=205_000, svlen=5_000),
        SV(
            id="inv1",
            svtype="INV",
            chrom="chr3",
            pos=300_000,
            end=310_000,
            svlen=10_000,
            strands="++",
        ),
        SV(id="ins1", svtype="INS", chrom="chr4", pos=400_000, end=400_000, svlen=200),
        SV(
            id="bnd1",
            svtype="BND",
            chrom="chr5",
            pos=500_000,
            end=500_000,
            svlen=1,
            mate_chrom="chr7",
            mate_pos=700_000,
            strands="+-",
        ),
    ]

def test_delly_one_record_per_event(tmp_path: Path) -> None:
    writer = get_writer("delly")
    svs = _svs()
    header = writer.header_lines("SAMP01")
    records = writer.format_records(svs, "SAMP01")
    assert len(records) == len(svs)

    out = tmp_path / "out.vcf"
    write_vcf(out, header, records)
    with pysam.VariantFile(str(out)) as vf:
        got = list(vf)
    assert len(got) == len(svs)

def test_delly_bnd_has_chr2_pos2(tmp_path: Path) -> None:
    writer = get_writer("delly")
    sv = SV(
        id="bnd1",
        svtype="BND",
        chrom="chr5",
        pos=500_000,
        end=500_000,
        svlen=1,
        mate_chrom="chr7",
        mate_pos=700_000,
        strands="+-",
    )
    out = tmp_path / "d.vcf"
    write_vcf(out, writer.header_lines("S"), writer.format_records([sv], "S"))
    with pysam.VariantFile(str(out)) as vf:
        rec = next(iter(vf))
    assert rec.info["CHR2"] == "chr7"
    assert int(rec.info["POS2"]) == 700_000
    assert "CT" in rec.info

def test_delly_roundtrip_bcf(tmp_path: Path) -> None:
    writer = get_writer("delly")
    svs = _svs()
    out = tmp_path / "out.bcf"
    write_vcf(out, writer.header_lines("S"), writer.format_records(svs, "S"))
    with pysam.VariantFile(str(out)) as vf:
        assert sum(1 for _ in vf) == len(svs)
