"""
Tests for the Manta writer with pysam round-trip validation
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

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
        SV(
            id="ins1",
            svtype="INS",
            chrom="chr4",
            pos=400_000,
            end=400_000,
            svlen=200,
            ins_seq="ACGT",
        ),
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

def test_manta_header_and_records_parse(tmp_path: Path) -> None:
    writer = get_writer("manta")
    svs = _svs()
    header = writer.header_lines("TUMOR01")
    records = writer.format_records(svs, "TUMOR01")

    out = tmp_path / "out.vcf"
    write_vcf(out, header, records)

    with pysam.VariantFile(str(out)) as vf:
        got = list(vf)
    # BND emits two records (mate pair); other svtypes emit one each
    assert len(got) == len(svs) + 1

    for rec in got:
        assert "SVTYPE" in rec.info

def test_manta_bnd_mates_have_mateid() -> None:
    writer = get_writer("manta")
    svs = [
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
        )
    ]
    records = writer.format_records(svs, "S")
    assert len(records) == 2
    assert "MATEID=bnd1_2" in records[0]
    assert "MATEID=bnd1_1" in records[1]

def test_manta_writes_vcf_gz_and_bcf(tmp_path: Path) -> None:
    writer = get_writer("manta")
    svs = _svs()
    header = writer.header_lines("S")
    records = writer.format_records(svs, "S")

    for suffix in (".vcf", ".vcf.gz", ".bcf"):
        out = tmp_path / f"out{suffix}"
        write_vcf(out, header, records)
        with pysam.VariantFile(str(out)) as vf:
            assert sum(1 for _ in vf) == len(svs) + 1

def test_manta_rejects_unknown_svtype() -> None:
    writer = get_writer("manta")
    bad = [
        SV(
            id="b1",
            svtype="DEL",
            chrom="chr1",
            pos=1_000,
            end=2_000,
            svlen=1_000,
        )
    ]
    # sanity: valid svtypes are accepted; CallerWriter validates itself
    writer.format_records(bad, "S")
    with pytest.raises(ValueError, match="svtype"):
        SV(id="b2", svtype="XYZ", chrom="chr1", pos=1, end=2, svlen=1)
