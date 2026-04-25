"""
Manta VCF writer

Reproduces the VCF 4.1 dialect emitted by Illumina Manta 1.6.0 for
somatic SV calls. The header is loaded verbatim from the reference
template shipped under :mod:`svforge.data.headers` so that downstream
parsers cannot distinguish svforge output from a real Manta run, except
for the ``##svforge*`` provenance tags

BND events are emitted as two mate records with ``MATEID`` cross-reference,
matching real Manta output. DEL/DUP/INV/INS events use symbolic alleles
"""

from __future__ import annotations

from typing import ClassVar

from svforge import __version__
from svforge.core.models import SV
from svforge.writers.base import CallerWriter


class MantaWriter(CallerWriter):
    """
    Manta 1.6.x VCF 4.1 writer
    """

    name: ClassVar[str] = "manta"
    supported_svtypes: ClassVar[frozenset[str]] = frozenset({"DEL", "DUP", "INV", "INS", "BND"})

    FILEFORMAT: ClassVar[str] = "VCFv4.1"
    SOURCE_STRING: ClassVar[str] = f"svforge-{__version__} (mimics GenerateSVCandidates 1.6.0)"
    SAMPLE_COLUMN_ORDER: ClassVar[tuple[str, ...]] = ("NORMAL", "TUMOR")
    TEMPLATE_BASENAME: ClassVar[str] = "manta_somatic"
    REQUIRED_PLACEHOLDERS: ClassVar[tuple[str, ...]] = (
        "{FILE_DATE}",
        "{REFERENCE_FASTA}",
        "{CMDLINE}",
        "{TUMOR_SAMPLE}",
        "{NORMAL_SAMPLE}",
    )

    def format_record(self, sv: SV, sample_name: str) -> list[str]:
        if sv.svtype == "BND":
            return _bnd_mate_records(sv)
        return [_symbolic_record(sv)]


def _sample_depth(vaf: float) -> tuple[int, int]:
    """
    Return (ref_count, alt_count) consistent with ``vaf`` at 30x total
    """
    total = 30
    alt = max(1, round(total * vaf))
    ref = max(0, total - alt)
    return ref, alt


def _sample_column(sv: SV) -> tuple[str, str]:
    """
    Return (FORMAT, sample) columns for a Manta record

    Only the FORMAT fields declared in the reference template (``PR``, ``SR``)
    are emitted so that the output parses cleanly through :mod:`pysam`
    """
    ref, alt = _sample_depth(sv.vaf)
    split_alt = max(1, alt // 2)
    split_ref = max(0, ref // 2)
    fmt = "PR:SR"
    sample = f"{ref},{alt}:{split_ref},{split_alt}"
    return fmt, sample


def _base_info(sv: SV) -> list[str]:
    info: list[str] = [f"SVTYPE={sv.svtype}"]
    if sv.svtype in {"DEL", "DUP", "INV"}:
        signed_len = -sv.svlen if sv.svtype == "DEL" else sv.svlen
        info.append(f"END={sv.end}")
        info.append(f"SVLEN={signed_len}")
    elif sv.svtype == "INS":
        info.append(f"END={sv.pos}")
        info.append(f"SVLEN={sv.svlen}")
        info.append(f"SVINSLEN={sv.svlen}")
        if sv.ins_seq:
            info.append(f"SVINSSEQ={sv.ins_seq}")
    if sv.homlen:
        info.append(f"HOMLEN={sv.homlen}")
        if sv.homseq:
            info.append(f"HOMSEQ={sv.homseq}")
    if sv.origin == "somatic" and sv.svtype != "BND":
        info.append("SOMATIC")
        info.append("SOMATICSCORE=60")
    info.append(f"SVFORGE_SOURCE={sv.source}")
    for key, value in sv.info_extra.items():
        info.append(f"{key}={value}" if value else key)
    return info


def _symbolic_record(sv: SV) -> str:
    alt = f"<{sv.svtype}>"
    info = ";".join(_base_info(sv))
    fmt, sample = _sample_column(sv)
    return "\t".join(
        [
            sv.chrom,
            str(sv.pos),
            sv.id,
            sv.ref_base,
            alt,
            ".",
            sv.filter,
            info,
            fmt,
            sample,
        ]
    )


def _bnd_alt(ref_base: str, mate_chrom: str, mate_pos: int, strands: str) -> str:
    """
    Return the Manta-style BND ALT string for a breakend

    VCF 4.2 mate syntax:

    - ``++`` (``t[p[``): REF joins region right of ``p``
    - ``+-`` (``t]p]``): REF joins region left of ``p``
    - ``-+`` (``[p[t``): region right of ``p`` joins REF
    - ``--`` (``]p]t``): region left of ``p`` joins REF
    """
    loc = f"{mate_chrom}:{mate_pos}"
    s1, s2 = strands[0], strands[1]
    if s1 == "+" and s2 == "+":
        return f"{ref_base}[{loc}["
    if s1 == "+" and s2 == "-":
        return f"{ref_base}]{loc}]"
    if s1 == "-" and s2 == "+":
        return f"[{loc}[{ref_base}"
    return f"]{loc}]{ref_base}"


def _mate_strands(strands: str) -> str:
    """
    Flip strands for the mate record (breakend pair)
    """
    s1, s2 = strands[0], strands[1]
    return f"{s2}{s1}"


def _bnd_mate_records(sv: SV) -> list[str]:
    if sv.mate_chrom is None or sv.mate_pos is None:
        raise ValueError(f"BND {sv.id!r} missing mate coordinates")
    id1 = f"{sv.id}_1"
    id2 = f"{sv.id}_2"

    alt1 = _bnd_alt(sv.ref_base, sv.mate_chrom, sv.mate_pos, sv.strands)
    alt2 = _bnd_alt(sv.ref_base, sv.chrom, sv.pos, _mate_strands(sv.strands))

    info1 = ["SVTYPE=BND", f"MATEID={id2}", f"EVENT={sv.id}"]
    info2 = ["SVTYPE=BND", f"MATEID={id1}", f"EVENT={sv.id}"]
    if sv.homlen:
        info1.append(f"HOMLEN={sv.homlen}")
        info2.append(f"HOMLEN={sv.homlen}")
    if sv.origin == "somatic":
        for info in (info1, info2):
            info.append("SOMATIC")
            info.append("SOMATICSCORE=60")
    info1.append(f"SVFORGE_SOURCE={sv.source}")
    info2.append(f"SVFORGE_SOURCE={sv.source}")

    fmt, sample = _sample_column(sv)
    rec1 = "\t".join(
        [
            sv.chrom,
            str(sv.pos),
            id1,
            sv.ref_base,
            alt1,
            ".",
            sv.filter,
            ";".join(info1),
            fmt,
            sample,
        ]
    )
    rec2 = "\t".join(
        [
            sv.mate_chrom,
            str(sv.mate_pos),
            id2,
            sv.ref_base,
            alt2,
            ".",
            sv.filter,
            ";".join(info2),
            fmt,
            sample,
        ]
    )
    return [rec1, rec2]


from svforge.writers._registry import register_writer  # noqa: E402

register_writer("manta")(MantaWriter)
