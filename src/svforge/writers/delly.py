"""
DELLY VCF writer

Reproduces the VCF 4.2 dialect emitted by DELLY 1.x for somatic SV calls.
The header is loaded verbatim from the reference template shipped under
:mod:`svforge.data.headers`; only ``##svforge*`` provenance tags, the
``##source`` line and dynamic fields (``##fileDate``) are rewritten.
Unlike Manta, DELLY emits a single record per BND event and encodes mate
coordinates via ``CHR2``/``POS2`` INFO fields plus a connection type
``CT`` (e.g. ``3to5`` for DEL-like, ``5to3`` for DUP-like)
"""

from __future__ import annotations

from typing import ClassVar

from svforge import __version__
from svforge.core.models import SV
from svforge.writers.base import CallerWriter


class DellyWriter(CallerWriter):
    """
    DELLY 1.x VCF 4.2 writer
    """

    name: ClassVar[str] = "delly"
    supported_svtypes: ClassVar[frozenset[str]] = frozenset({"DEL", "DUP", "INV", "INS", "BND"})

    FILEFORMAT: ClassVar[str] = "VCFv4.2"
    SOURCE_STRING: ClassVar[str] = f"svforge-{__version__} (mimics DELLY)"
    SAMPLE_COLUMN_ORDER: ClassVar[tuple[str, ...]] = ("TUMOR", "NORMAL")
    TEMPLATE_BASENAME: ClassVar[str] = "delly_somatic"
    REQUIRED_PLACEHOLDERS: ClassVar[tuple[str, ...]] = (
        "{FILE_DATE}",
        "{REFERENCE_FASTA}",
        "{TUMOR_SAMPLE}",
        "{NORMAL_SAMPLE}",
    )

    def format_record(self, sv: SV, sample_name: str) -> list[str]:
        return [_delly_record(sv)]


def _ct_for(sv: SV) -> str:
    """
    Return the DELLY connection-type string for an SV
    """
    match sv.svtype:
        case "DEL":
            return "3to5"
        case "DUP":
            return "5to3"
        case "INV":
            return "3to3" if sv.strands == "++" else "5to5"
        case "INS":
            return "3to5"
        case "BND":
            s1, s2 = sv.strands[0], sv.strands[1]
            left = "3" if s1 == "+" else "5"
            right = "5" if s2 == "-" else "3"
            return f"{left}to{right}"
        case _:  # pragma: no cover
            return "NtoN"


def _sample_column(sv: SV) -> tuple[str, str]:
    total = 30
    alt = max(1, round(total * sv.vaf))
    ref = max(0, total - alt)
    junction_alt = max(1, alt // 2)
    junction_ref = max(0, ref // 2)
    gq = 60
    gl = (
        "-50,-5,0"
        if sv.genotype == "1/1"
        else ("0,-5,-50" if sv.genotype == "0/0" else "-10,-1,-10")
    )
    fmt = "GT:GL:GQ:FT:RC:RCL:RCR:RDCN:DR:DV:RR:RV"
    sample = f"{sv.genotype}:{gl}:{gq}:PASS:40:20:20:2:{ref}:{alt}:{junction_ref}:{junction_alt}"
    return fmt, sample


def _base_info(sv: SV) -> list[str]:
    ct = _ct_for(sv)
    end = sv.end if sv.svtype != "BND" else sv.pos
    chr2 = sv.mate_chrom if sv.svtype == "BND" else sv.chrom
    pos2 = sv.mate_pos if sv.svtype == "BND" else sv.end
    info = [
        "PRECISE",
        f"SVTYPE={sv.svtype}",
        "SVMETHOD=EMBL.DELLYv1.2.6-svforge",
        f"CHR2={chr2}",
        f"END={end}",
        f"POS2={pos2}",
        "PE=20",
        "MAPQ=60",
        f"CT={ct}",
        "CIPOS=-10,10",
        "CIEND=-10,10",
        "SR=5",
        "SRQ=1",
        "CE=1.9",
    ]
    if sv.svtype == "INS":
        info.append(f"SVLEN={sv.svlen}")
    if sv.homlen:
        info.append(f"HOMLEN={sv.homlen}")
    if sv.origin == "somatic":
        info.append("SOMATIC")
    info.append(f"SVFORGE_SOURCE={sv.source}")
    for key, value in sv.info_extra.items():
        info.append(f"{key}={value}" if value else key)
    return info


def _delly_record(sv: SV) -> str:
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


from svforge.writers._registry import register_writer  # noqa: E402

register_writer("delly")(DellyWriter)
