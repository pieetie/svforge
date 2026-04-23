"""
Blacklist annotation: tag SVs that overlap ENCODE blacklist regions

A record is flagged ``poor_mappability=True`` in its INFO if any of its
breakpoints falls within a blacklist interval. For intra-chromosomal SVs,
we also check the spanning interval [POS, END]; for BND we check both
breakends.

Counts returned to the validate report are **event-based** (a Manta BND
mate pair is one SV event, not two records), see :func:`event_id`
"""

from __future__ import annotations

from pathlib import Path

import pysam

from svforge.core.regions import RegionSet

POOR_MAPPABILITY_FLAG = "poor_mappability"
POOR_MAPPABILITY_HEADER = (
    f"##INFO=<ID={POOR_MAPPABILITY_FLAG},Number=0,Type=Flag,"
    f'Description="SV overlaps an ENCODE blacklist region">'
)


def event_id(rec: pysam.VariantRecord) -> tuple[str, str]:
    """
    Return a canonical event identifier for ``rec``

    An SV event may span multiple VCF records (Manta emits a BND as a
    mate pair). This helper returns a hashable key that is identical
    across the records belonging to the same event so that counts are
    taken at the event level, not the record level

    Resolution order:

    1. ``INFO/EVENT`` (Manta groups mates under a shared event ID)
    2. ``INFO/MATEID`` canonicalised as a sorted pair with the record ID
    3. ``ID`` fallback (record-unique; one event per record)
    """
    event = _first_info(rec, "EVENT")
    if event is not None:
        return ("event", str(event))
    mate = _first_info(rec, "MATEID")
    if mate is not None and rec.id:
        pair = tuple(sorted([str(rec.id), str(mate)]))
        return ("mate", f"{pair[0]}|{pair[1]}")
    return ("id", str(rec.id or f"_anon_{rec.chrom}_{rec.pos}"))


def annotate_vcf(
    vcf_in: str | Path,
    vcf_out: str | Path,
    blacklist: RegionSet,
) -> tuple[int, int]:
    """
    Copy ``vcf_in`` to ``vcf_out`` tagging blacklist-overlapping SVs

    Returns ``(n_total_events, n_flagged_events)`` -- mate pairs count as
    a single event (see :func:`event_id`). Per-record flagging is
    preserved in the output VCF so downstream consumers still see the
    ``poor_mappability`` tag on every record whose own breakend falls
    inside a blacklisted interval
    """
    in_path = Path(vcf_in)
    out_path = Path(vcf_out)
    out_path.parent.mkdir(parents=True, exist_ok=True)

    total_events: set[tuple[str, str]] = set()
    flagged_events: set[tuple[str, str]] = set()

    with pysam.VariantFile(str(in_path)) as vin:
        if POOR_MAPPABILITY_FLAG not in vin.header.info:
            vin.header.add_line(POOR_MAPPABILITY_HEADER)

        mode = _write_mode_for(out_path)
        with pysam.VariantFile(str(out_path), mode, header=vin.header) as vout:  # type: ignore[arg-type]
            for rec in vin:
                eid = event_id(rec)
                total_events.add(eid)
                if _record_overlaps(rec, blacklist):
                    rec.info[POOR_MAPPABILITY_FLAG] = True
                    flagged_events.add(eid)
                vout.write(rec)

    return len(total_events), len(flagged_events)


def _write_mode_for(path: Path) -> str:
    name = path.name.lower()
    if name.endswith(".vcf.gz"):
        return "wz"
    if name.endswith(".bcf"):
        return "wb"
    return "w"


def _record_overlaps(rec: pysam.VariantRecord, blacklist: RegionSet) -> bool:
    chrom = str(rec.chrom)
    pos = int(rec.pos)
    svtype = _first_info(rec, "SVTYPE")

    if svtype == "BND":
        if blacklist.overlaps_point(chrom, pos):
            return True
        chr2 = _first_info(rec, "CHR2")
        pos2 = _first_info(rec, "POS2")
        if chr2 is not None and pos2 is not None:
            return blacklist.overlaps_point(str(chr2), _to_int(pos2))
        return False

    end = int(rec.stop)
    return blacklist.overlaps_vcf(chrom, pos, end)


def _to_int(value: object) -> int:
    """
    Safely coerce a pysam INFO scalar (already proven non-None) to ``int``
    """
    if isinstance(value, (int, str, float)):
        return int(value)
    return int(str(value))


def _first_info(rec: pysam.VariantRecord, key: str) -> object:
    """
    Return the first value for INFO key ``key`` or ``None`` when absent

    pysam raises ValueError if ``key`` is not declared in the header, so we
    guard against both that and normal absence
    """
    if key not in rec.header.info:
        return None
    try:
        value = rec.info.get(key)
    except (KeyError, ValueError):
        return None
    if isinstance(value, tuple):
        return value[0] if value else None
    return value
