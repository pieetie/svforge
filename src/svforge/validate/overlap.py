"""
gnomAD SV overlap with 2 bp breakpoint tolerance

Semantics: an SV matches gnomAD if there exists a gnomAD SV with both
breakpoints within :data:`DEFAULT_TOLERANCE_BP` base-pairs on the same
chromosomes. This mirrors how clinical SV pipelines actually filter
germline events -- the svtype is not part of the match condition since
callers often disagree on DEL vs CNV vs DUP for the same locus
"""

from __future__ import annotations

from collections import defaultdict
from collections.abc import Iterator
from dataclasses import dataclass
from pathlib import Path

import pysam

from svforge.core.models import SV

DEFAULT_TOLERANCE_BP = 2


@dataclass(slots=True, frozen=True)
class _IndexedSV:
    svtype: str
    chrom: str
    pos: int
    chrom2: str
    pos2: int


class GnomadIndex:
    """
    Chromosome-keyed index of gnomAD SV sites for tolerant breakpoint lookup
    """

    __slots__ = ("_by_chrom",)

    def __init__(self) -> None:
        self._by_chrom: dict[str, list[_IndexedSV]] = defaultdict(list)

    def add(self, entry: _IndexedSV) -> None:
        self._by_chrom[entry.chrom].append(entry)

    def overlaps(self, sv: SV, tolerance: int = DEFAULT_TOLERANCE_BP) -> bool:
        """
        True if ``sv`` has a gnomAD match with both breakpoints within ``tolerance``
        """
        sv_chrom2 = sv.mate_chrom if sv.svtype == "BND" else sv.chrom
        sv_pos2 = sv.mate_pos if sv.svtype == "BND" else sv.end
        if sv_chrom2 is None or sv_pos2 is None:
            return False

        for entry in self._by_chrom.get(sv.chrom, ()):
            if abs(entry.pos - sv.pos) > tolerance:
                continue
            if entry.chrom2 != sv_chrom2:
                continue
            if abs(entry.pos2 - sv_pos2) > tolerance:
                continue
            return True
        return False

    def __len__(self) -> int:
        return sum(len(v) for v in self._by_chrom.values())


def load_gnomad_index(vcf_path: str | Path) -> GnomadIndex:
    """
    Parse a gnomAD SV VCF into a :class:`GnomadIndex`
    """
    index = GnomadIndex()
    for entry in _iter_gnomad(Path(vcf_path)):
        index.add(entry)
    return index


def _iter_gnomad(path: Path) -> Iterator[_IndexedSV]:
    with pysam.VariantFile(str(path)) as vf:
        for rec in vf:
            svtype = _extract_svtype(rec)
            if svtype is None:
                continue
            chrom = str(rec.chrom)
            pos = int(rec.pos)
            chrom2, pos2 = _extract_mate_coords(rec, svtype, chrom, pos)
            yield _IndexedSV(svtype=svtype, chrom=chrom, pos=pos, chrom2=chrom2, pos2=pos2)


def _extract_svtype(rec: pysam.VariantRecord) -> str | None:
    svtype = _safe_info(rec, "SVTYPE")
    if isinstance(svtype, tuple):
        svtype = svtype[0] if svtype else None
    if svtype is None:
        return None
    return str(svtype).upper()


def _extract_mate_coords(
    rec: pysam.VariantRecord,
    svtype: str,
    chrom: str,
    pos: int,
) -> tuple[str, int]:
    chrom2 = _safe_info(rec, "CHR2") or _safe_info(rec, "CHROM2")
    pos2 = _safe_info(rec, "POS2") or _safe_info(rec, "END2")
    end = _safe_info(rec, "END")
    if end is None:
        end = rec.stop

    if isinstance(chrom2, tuple):
        chrom2 = chrom2[0] if chrom2 else None
    if isinstance(pos2, tuple):
        pos2 = pos2[0] if pos2 else None

    if svtype == "BND" and chrom2 is not None and pos2 is not None:
        return str(chrom2), _to_int(pos2)
    if end is None:
        return chrom, pos
    return chrom, _to_int(end)


def _to_int(value: object) -> int:
    """
    Safely coerce a pysam INFO scalar (already proven non-None) to ``int``
    """
    if isinstance(value, (int, str, float)):
        return int(value)
    return int(str(value))


def _safe_info(rec: pysam.VariantRecord, key: str) -> object:
    """
    pysam raises ValueError if ``key`` is not declared in the VCF header

    We want missing-key to mean "absent" when polling optional INFO tags
    across heterogeneous gnomAD / caller dialects
    """
    if key not in rec.header.info:
        return None
    try:
        return rec.info.get(key)
    except (KeyError, ValueError):
        return None
