"""
Integration tests for :func:`~svforge.core.sampler.sample_pair`: writing tumor vs.
normal-only SV lists as separate single-sample VCFs (distinct from CLI ``gen-pair``).
"""

from __future__ import annotations

from pathlib import Path

import pysam

from svforge.core.bank import Bank
from svforge.core.sampler import SamplerConfig, sample_pair
from svforge.io.vcf_writer import write_vcf
from svforge.writers import get_writer


def test_sample_pair_writes_two_consistent_vcfs(tmp_path: Path, mini_bank: Bank) -> None:
    pair = sample_pair(
        mini_bank,
        n_somatic=8,
        n_germline=10,
        cfg=SamplerConfig(
            n=0,
            seed=13,
            svtypes=frozenset({"DEL", "DUP", "INV"}),
            vaf_range=(0.3, 0.9),
        ),
    )

    writer = get_writer("manta")
    tumor_out = tmp_path / "tumor.vcf.gz"
    normal_out = tmp_path / "normal.vcf.gz"

    tumor_hdr = writer.header_lines("TUMOR01")
    normal_hdr = writer.header_lines("NORMAL01")
    write_vcf(
        tumor_out,
        tumor_hdr,
        writer.format_records_sorted(pair.tumor, "TUMOR01", tumor_hdr),
    )
    write_vcf(
        normal_out,
        normal_hdr,
        writer.format_records_sorted(pair.normal, "NORMAL01", normal_hdr),
    )

    with pysam.VariantFile(str(tumor_out)) as vf:
        tumor_ids = {rec.id for rec in vf if rec.id}
    with pysam.VariantFile(str(normal_out)) as vf:
        normal_ids = {rec.id for rec in vf if rec.id}

    assert normal_ids.issubset(tumor_ids)
    somatic_ids = tumor_ids - normal_ids
    assert len(somatic_ids) >= 8
    assert len(normal_ids) >= 10
