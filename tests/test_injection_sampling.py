"""
Injection sampling: SVFORGE_SOURCE counts + exact catalog match
"""

from __future__ import annotations

from pathlib import Path

import pysam
import pytest

from svforge.cli import main
from svforge.core.bank import Bank
from svforge.core.injection_catalogs import load_blacklist_catalog, load_gnomad_catalog
from svforge.core.sampler import SamplerConfig, sample

def _record_coords(rec: pysam.VariantRecord) -> tuple[str, int, int]:
    svtype = rec.info.get("SVTYPE")
    if isinstance(svtype, tuple):
        svtype = svtype[0] if svtype else None
    pos = int(rec.pos)
    if str(svtype) == "INS":
        svlen = rec.info.get("SVLEN")
        if isinstance(svlen, tuple):
            svlen = svlen[0] if svlen else None
        if svlen is not None:
            return (str(rec.chrom), pos, pos + abs(int(svlen)))
    return (str(rec.chrom), pos, int(rec.stop))

def _count_sources(vcf_path: Path) -> dict[str, int]:
    counts: dict[str, int] = {"bank": 0, "gnomad": 0, "blacklist": 0}
    seen_events: set[str] = set()
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            event = rec.info.get("EVENT")
            if isinstance(event, tuple):
                event = event[0] if event else None
            key = str(event) if event else f"{rec.chrom}:{rec.pos}:{rec.id}"
            if event and event in seen_events:
                continue
            seen_events.add(key)
            src = rec.info.get("SVFORGE_SOURCE")
            if isinstance(src, tuple):
                src = src[0] if src else None
            if src:
                counts[str(src)] = counts.get(str(src), 0) + 1
    return counts

def test_gnomad_fraction_produces_matching_records(tmp_path: Path, mini_bank_path: Path) -> None:
    out = tmp_path / "g.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--out",
            str(out),
            "--n",
            "50",
            "--sample-name",
            "T",
            "--bank",
            str(mini_bank_path),
            "--seed",
            "42",
            "--svtypes",
            "DEL,DUP,INV,INS",
            "--gnomad-fraction",
            "0.20",
        ]
    )
    assert rc == 0
    counts = _count_sources(out)
    assert counts["gnomad"] == 10

    gnomad_coords = {(e.chrom, e.pos, e.end) for e in load_gnomad_catalog()}
    with pysam.VariantFile(str(out)) as vf:
        for rec in vf:
            src = rec.info.get("SVFORGE_SOURCE")
            if isinstance(src, tuple):
                src = src[0] if src else None
            if str(src) != "gnomad":
                continue
            coord = _record_coords(rec)
            assert coord in gnomad_coords, coord

def test_blacklist_fraction_produces_matching_records(tmp_path: Path, mini_bank_path: Path) -> None:
    out = tmp_path / "b.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--out",
            str(out),
            "--n",
            "50",
            "--sample-name",
            "T",
            "--bank",
            str(mini_bank_path),
            "--seed",
            "99",
            "--svtypes",
            "DEL,DUP,INV",
            "--blacklist-fraction",
            "0.20",
        ]
    )
    assert rc == 0
    counts = _count_sources(out)
    assert counts["blacklist"] == 10

    bl_coords = {(e.chrom, e.pos, e.end) for e in load_blacklist_catalog()}
    with pysam.VariantFile(str(out)) as vf:
        for rec in vf:
            src = rec.info.get("SVFORGE_SOURCE")
            if isinstance(src, tuple):
                src = src[0] if src else None
            if str(src) != "blacklist":
                continue
            coord = _record_coords(rec)
            assert coord in bl_coords, coord

def test_combined_fractions(tmp_path: Path, mini_bank_path: Path) -> None:
    out = tmp_path / "c.vcf"
    rc = main(
        [
            "gen",
            "--caller",
            "manta",
            "--out",
            str(out),
            "--n",
            "100",
            "--sample-name",
            "T",
            "--bank",
            str(mini_bank_path),
            "--seed",
            "1",
            "--svtypes",
            "DEL,DUP,INV",
            "--gnomad-fraction",
            "0.20",
            "--blacklist-fraction",
            "0.10",
        ]
    )
    assert rc == 0
    counts = _count_sources(out)
    assert counts["gnomad"] == 20
    assert counts["blacklist"] == 10
    assert counts["bank"] >= 70

def test_sum_of_fractions_over_one_is_rejected(mini_bank: Bank) -> None:
    with pytest.raises(ValueError, match=r"<= 1\.0"):
        sample(
            mini_bank,
            SamplerConfig(n=10, seed=0, gnomad_fraction=0.6, blacklist_fraction=0.5),
        )

def test_sampling_without_replacement_when_possible(mini_bank: Bank) -> None:
    svs = sample(
        mini_bank,
        SamplerConfig(
            n=50,
            seed=3,
            svtypes=frozenset({"DEL", "DUP", "INV", "INS"}),
            gnomad_fraction=0.20,
        ),
    )
    gnomad_ids = [sv.info_extra.get("SVFORGE_SOURCE_ID") for sv in svs if sv.source == "gnomad"]
    assert len(gnomad_ids) == 10
    assert len(set(gnomad_ids)) == 10
