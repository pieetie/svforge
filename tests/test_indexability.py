"""
Tests for VCF indexability and sort order.

These tests guard against regressions of issue #6 (BND inter-chrom sort).
"""

from __future__ import annotations

import re
from collections import defaultdict
from pathlib import Path

import pysam
import pytest

from svforge.cli import main as cli_main

FIXTURES_DIR = Path(__file__).parent / "fixtures"
BND_ONLY_BANK = FIXTURES_DIR / "bnd_only_hg38.yaml"


def _run_svforge(args: list[str]) -> None:
    rc = cli_main(args)
    if rc != 0:
        raise RuntimeError(f"svforge exited with {rc}")


def _extract_contig_order(vcf_path: Path) -> list[str]:
    contigs: list[str] = []
    with pysam.VariantFile(str(vcf_path)) as vf:
        for line in str(vf.header).splitlines():
            m = re.match(r"^##contig=<ID=([^,>]+)", line)
            if m:
                contigs.append(m.group(1))
    return contigs


def _records_in_order(vcf_path: Path) -> bool:
    contig_order = _extract_contig_order(vcf_path)
    contig_idx = {c: i for i, c in enumerate(contig_order)}
    last_key = (-1, -1)
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            key = (contig_idx.get(rec.chrom, 10**9), rec.pos)
            if key < last_key:
                return False
            last_key = key
    return True


def _info_first(rec: pysam.VariantRecord, key: str) -> str | int | None:
    v = rec.info.get(key)
    if v is None:
        return None
    return v[0] if isinstance(v, (tuple, list)) else v


def _has_inter_chrom_bnd(caller: str, vcf_path: Path) -> bool:
    """True if at least one BND with primary and mate on different contigs."""
    with pysam.VariantFile(str(vcf_path)) as vf:
        if caller == "delly":
            for rec in vf:
                if _info_first(rec, "SVTYPE") != "BND":
                    continue
                chr2 = _info_first(rec, "CHR2")
                if chr2 is not None and rec.chrom != str(chr2):
                    return True
            return False
        # manta: two lines per EVENT; inter-chrom iff the two breakends differ in CHROM
        by_event: dict[str, set[str]] = defaultdict(set)
        for rec in vf:
            if _info_first(rec, "SVTYPE") != "BND":
                continue
            event = _info_first(rec, "EVENT")
            if event is None:
                continue
            by_event[str(event)].add(rec.chrom)
        return any(len(chroms) > 1 for chroms in by_event.values())


@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_tabix_indexable_simple(tmp_path: Path, caller: str) -> None:
    """Generated VCFs without BND should be tabix-indexable."""
    out = tmp_path / "test.vcf.gz"
    _run_svforge(
        [
            "gen",
            "--caller",
            caller,
            "--n",
            "20",
            "--sample-name",
            "TEST",
            "--seed",
            "42",
            "--svtypes",
            "DEL,DUP,INV,INS",
            "--out",
            str(out),
        ]
    )
    pysam.tabix_index(str(out), preset="vcf", force=True)
    assert (tmp_path / "test.vcf.gz.tbi").exists()


@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_tabix_indexable_with_bnd(tmp_path: Path, caller: str) -> None:
    """Generated VCFs with BND inter-chrom must remain tabix-indexable (issue #6)."""
    out = tmp_path / "test.vcf.gz"
    _run_svforge(
        [
            "gen",
            "--caller",
            caller,
            "--n",
            "30",
            "--sample-name",
            "TEST",
            "--seed",
            "42",
            "--bank",
            str(BND_ONLY_BANK),
            "--svtypes",
            "BND",
            "--out",
            str(out),
        ]
    )
    with pysam.VariantFile(str(out)) as vf:
        bnd_count = sum(1 for r in vf if _info_first(r, "SVTYPE") == "BND")
    assert bnd_count > 0, (
        "Test setup error: no BND records generated; issue #6 regression is not exercised."
    )
    assert _has_inter_chrom_bnd(caller, out), (
        "Test setup error: no inter-chromosomal BND; issue #6 regression is not exercised."
    )
    pysam.tabix_index(str(out), preset="vcf", force=True)
    assert (tmp_path / "test.vcf.gz.tbi").exists()


@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_records_sorted_by_header_contig_order(tmp_path: Path, caller: str) -> None:
    """All records must be sorted by (contig_index, pos), header order."""
    out = tmp_path / "test.vcf"
    _run_svforge(
        [
            "gen",
            "--caller",
            caller,
            "--n",
            "50",
            "--sample-name",
            "TEST",
            "--seed",
            "42",
            "--gnomad-fraction",
            "0.3",
            "--out",
            str(out),
        ]
    )
    assert _records_in_order(out), (
        f"Records in {out} are not sorted by (contig_index, pos)"
    )


def test_region_query_chr1_after_index(tmp_path: Path) -> None:
    """After indexing, pysam fetch for chr1 returns only chr1 records."""
    out = tmp_path / "test.vcf.gz"
    _run_svforge(
        [
            "gen",
            "--caller",
            "manta",
            "--n",
            "100",
            "--sample-name",
            "TEST",
            "--seed",
            "42",
            "--out",
            str(out),
        ]
    )
    pysam.tabix_index(str(out), preset="vcf", force=True)
    with pysam.VariantFile(str(out)) as vf:
        chr1_records = list(vf.fetch("chr1"))
    assert all(r.chrom == "chr1" for r in chr1_records)


def test_sort_preserves_record_content(tmp_path: Path) -> None:
    """Re-sorting data lines by header (contig index, pos) is a no-op — order matches sort key."""
    out = tmp_path / "test.vcf"
    _run_svforge(
        [
            "gen",
            "--caller",
            "manta",
            "--n",
            "30",
            "--sample-name",
            "TEST",
            "--seed",
            "42",
            "--gnomad-fraction",
            "0.3",
            "--out",
            str(out),
        ]
    )
    with out.open(encoding="utf-8") as fh:
        records = [line.rstrip("\n") for line in fh if not line.startswith("#")]
    contig_order = _extract_contig_order(out)
    contig_idx = {c: i for i, c in enumerate(contig_order)}

    def _line_sort_key(line: str) -> tuple[int, int]:
        parts = line.split("\t", 2)
        chrom, pos_s = parts[0], parts[1]
        return (contig_idx.get(chrom, 10**9), int(pos_s))

    sorted_records = sorted(records, key=_line_sort_key)
    assert records == sorted_records
