"""
BED region index with pure-Python interval trees

Designed for blacklist- and gnomAD-region overlap queries. Supports
plain and gzip-compressed BED files via the standard library only; no
system dependency on tabix or bgzip is introduced here
"""

from __future__ import annotations

import gzip
from collections import defaultdict
from collections.abc import Iterable, Iterator
from pathlib import Path

from intervaltree import IntervalTree


class RegionSet:
    """
    In-memory set of half-open [start, end) intervals per chromosome

    BED uses 0-based half-open coordinates. We keep that convention
    internally and convert at the boundaries (callers usually work in
    1-based closed VCF coordinates)
    """

    __slots__ = ("_trees",)

    def __init__(self, intervals: Iterable[tuple[str, int, int]] | None = None) -> None:
        self._trees: dict[str, IntervalTree] = defaultdict(IntervalTree)
        if intervals is not None:
            for chrom, start, end in intervals:
                self.add(chrom, start, end)

    def add(self, chrom: str, start: int, end: int) -> None:
        """
        Add one half-open BED interval
        """
        if end <= start:
            return
        self._trees[chrom].addi(start, end)

    def overlaps_bed(self, chrom: str, start: int, end: int) -> bool:
        """
        True if any indexed region overlaps BED interval [start, end)
        """
        tree = self._trees.get(chrom)
        if tree is None or end <= start:
            return False
        return bool(tree.overlaps(start, end))

    def overlaps_vcf(self, chrom: str, pos: int, end: int) -> bool:
        """
        True if any indexed region overlaps a 1-based closed [pos, end] SV
        """
        return self.overlaps_bed(chrom, pos - 1, end)

    def overlaps_point(self, chrom: str, pos: int) -> bool:
        """
        True if any indexed region contains the 1-based VCF position ``pos``
        """
        tree = self._trees.get(chrom)
        if tree is None:
            return False
        return bool(tree.at(pos - 1))

    def __len__(self) -> int:
        return sum(len(t) for t in self._trees.values())

    @property
    def chromosomes(self) -> set[str]:
        return set(self._trees.keys())

    def intervals_on(self, chrom: str) -> list[tuple[int, int]]:
        """
        Return ``[(start, end)]`` BED intervals indexed on ``chrom``
        """
        tree = self._trees.get(chrom)
        if tree is None:
            return []
        return [(iv.begin, iv.end) for iv in tree]

    def all_intervals(self) -> list[tuple[str, int, int]]:
        """
        Return ``[(chrom, start, end)]`` for every indexed interval
        """
        return [(chrom, iv.begin, iv.end) for chrom, tree in self._trees.items() for iv in tree]


def load_bed(path: str | Path) -> RegionSet:
    """
    Load a (optionally gzip-compressed) BED file into a :class:`RegionSet`

    Only the first three columns are read. Lines starting with ``#`` or
    the ``browser``/``track`` header lines are skipped
    """
    rs = RegionSet()
    for chrom, start, end in _iter_bed(Path(path)):
        rs.add(chrom, start, end)
    return rs


def _iter_bed(path: Path) -> Iterator[tuple[str, int, int]]:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "rt", encoding="utf-8") as fh:
        for raw in fh:
            line = raw.strip()
            if not line or line.startswith(("#", "browser", "track")):
                continue
            parts = line.split("\t")
            if len(parts) < 3:
                continue
            try:
                start = int(parts[1])
                end = int(parts[2])
            except ValueError:
                continue
            yield parts[0], start, end
