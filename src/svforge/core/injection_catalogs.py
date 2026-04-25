"""
Load and expose the bundled real-catalog injection entries

Two TSV files ship in :mod:`svforge.data.injections`:

- ``gnomad_hg38_mini.tsv`` -- curated real gnomAD v4.1 SV sites
- ``blacklist_hg38_mini.tsv`` -- curated ENCODE hg38 blacklist v2 regions

When ``svforge gen`` is invoked with ``--gnomad-fraction`` or
``--blacklist-fraction`` the sampler pulls the requested number of
entries from these catalogs verbatim (no synthesis, no spatial bias).
Downstream pipelines recognise the injected SVs directly because the
CHROM/POS/END coincide with real catalog rows, giving ``svforge
validate`` a zero-ambiguity self-consistency check to verify that a
generated VCF was actually produced by this build
"""

from __future__ import annotations

import csv
from collections.abc import Iterator
from dataclasses import dataclass
from functools import cache
from importlib import resources

GNOMAD_TSV = "gnomad_hg38_mini.tsv"
BLACKLIST_TSV = "blacklist_hg38_mini.tsv"

@dataclass(frozen=True, slots=True)
class GnomadEntry:
    """
    One curated gnomAD SV site
    """

    chrom: str
    pos: int
    end: int
    end_chrom: str
    svtype: str
    source_id: str

@dataclass(frozen=True, slots=True)
class BlacklistEntry:
    """
    One curated ENCODE blacklist region
    """

    chrom: str
    pos: int
    end: int
    region_type: str
    source_id: str


@cache
def load_gnomad_catalog() -> tuple[GnomadEntry, ...]:
    """
    Return the bundled gnomAD mini catalog, cached for the process life
    """
    return tuple(_iter_gnomad(_read_tsv(GNOMAD_TSV)))


@cache
def load_blacklist_catalog() -> tuple[BlacklistEntry, ...]:
    """
    Return the bundled blacklist mini catalog, cached for the process life
    """
    return tuple(_iter_blacklist(_read_tsv(BLACKLIST_TSV)))

def _read_tsv(name: str) -> list[dict[str, str]]:
    resource = resources.files("svforge.data.injections").joinpath(name)
    if not resource.is_file():
        raise RuntimeError(f"Injection catalog {name!r} is missing from the package")
    text = resource.read_text(encoding="utf-8")
    reader = csv.DictReader(text.splitlines(), delimiter="\t")
    rows = list(reader)
    if not rows:
        raise RuntimeError(f"Injection catalog {name!r} is empty")
    return rows



def _iter_gnomad(rows: list[dict[str, str]]) -> Iterator[GnomadEntry]:
    for row in rows:
        yield GnomadEntry(
            chrom=row["chrom"],
            pos=int(row["pos"]),
            end=int(row["end"]),
            end_chrom=row.get("end_chrom") or row["chrom"],
            svtype=row["svtype"].upper(),
            source_id=row["source_id"],
        )


def _iter_blacklist(rows: list[dict[str, str]]) -> Iterator[BlacklistEntry]:
    for row in rows:
        yield BlacklistEntry(
            chrom=row["chrom"],
            pos=int(row["pos"]),
            end=int(row["end"]),
            region_type=row["region_type"],
            source_id=row["source_id"],
        )
