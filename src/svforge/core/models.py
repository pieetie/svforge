"""
Domain models shared by sampler, writers and validator

These dataclasses are the contract that every other module depends on.
A writer only ever receives :class:`SV` instances and must never assume
anything beyond what is declared here
"""

from __future__ import annotations

from dataclasses import dataclass, field
from typing import Literal

SVType = Literal["DEL", "DUP", "INV", "INS", "BND"]
Genotype = Literal["0/0", "0/1", "1/1"]
Origin = Literal["somatic", "germline"]
SVSource = Literal["bank", "gnomad", "blacklist"]

VALID_SVTYPES: frozenset[str] = frozenset({"DEL", "DUP", "INV", "INS", "BND"})
VALID_GENOTYPES: frozenset[str] = frozenset({"0/0", "0/1", "1/1"})
VALID_SOURCES: frozenset[str] = frozenset({"bank", "gnomad", "blacklist"})


@dataclass(slots=True)
class Breakpoint:
    """
    A single breakpoint on the reference

    ``strand`` encodes the orientation at the breakpoint: ``+`` means the
    joined sequence extends to the right of ``pos``, ``-`` means to the left
    """

    chrom: str
    pos: int
    strand: Literal["+", "-"] = "+"

    def __post_init__(self) -> None:
        if self.pos < 1:
            raise ValueError(f"Breakpoint pos must be 1-based positive, got {self.pos}")
        if self.strand not in {"+", "-"}:
            raise ValueError(f"Breakpoint strand must be '+' or '-', got {self.strand!r}")


@dataclass(slots=True)
class SV:
    """
    Caller-agnostic structural-variant record

    Writers format this into their own VCF dialect. The sampler produces it
    with deterministic fields so two runs with the same seed yield identical
    records
    """

    id: str
    svtype: str
    chrom: str
    pos: int
    end: int
    svlen: int
    mate_chrom: str | None = None
    mate_pos: int | None = None
    strands: str = "+-"
    homlen: int = 0
    homseq: str = ""
    vaf: float = 0.5
    genotype: str = "0/1"
    ref_base: str = "N"
    ins_seq: str = ""
    filter: str = "PASS"
    origin: Origin = "somatic"
    source: SVSource = "bank"
    info_extra: dict[str, str] = field(default_factory=dict)

    def __post_init__(self) -> None:
        if self.svtype not in VALID_SVTYPES:
            raise ValueError(f"Unknown svtype {self.svtype!r} (valid: {sorted(VALID_SVTYPES)})")
        if self.genotype not in VALID_GENOTYPES:
            raise ValueError(
                f"Unknown genotype {self.genotype!r} (valid: {sorted(VALID_GENOTYPES)})"
            )
        if self.pos < 1:
            raise ValueError(f"SV pos must be 1-based positive, got {self.pos}")
        if self.end < self.pos and self.svtype != "BND":
            raise ValueError(f"SV end ({self.end}) must be >= pos ({self.pos}) for {self.svtype}")
        if self.homlen < 0:
            raise ValueError(f"HOMLEN must be non-negative, got {self.homlen}")
        if not 0.0 <= self.vaf <= 1.0:
            raise ValueError(f"VAF must be in [0, 1], got {self.vaf}")
        if self.svtype == "BND" and (self.mate_chrom is None or self.mate_pos is None):
            raise ValueError("BND SV requires mate_chrom and mate_pos")
        if len(self.strands) != 2 or any(c not in "+-" for c in self.strands):
            raise ValueError(f"strands must be two chars in '+-', got {self.strands!r}")
        if self.source not in VALID_SOURCES:
            raise ValueError(f"Unknown source {self.source!r} (valid: {sorted(VALID_SOURCES)})")

    @property
    def breakpoint1(self) -> Breakpoint:
        """
        Primary breakpoint (chrom, pos, strand1)
        """
        return Breakpoint(self.chrom, self.pos, self.strands[0])  # type: ignore[arg-type]

    @property
    def breakpoint2(self) -> Breakpoint:
        """
        Secondary breakpoint (mate for BND, end for intra-chrom SVs)
        """
        chrom = self.mate_chrom if self.svtype == "BND" else self.chrom
        pos = self.mate_pos if self.svtype == "BND" else self.end
        assert chrom is not None and pos is not None
        return Breakpoint(chrom, pos, self.strands[1])  # type: ignore[arg-type]


@dataclass(slots=True)
class SVPair:
    """
    Tumor + normal SV sets for paired somatic pipelines

    ``tumor`` contains all SVs present in the tumor VCF (somatic + germline)
    and ``normal`` contains only the germline subset. Germline SVs appear in
    both lists with the same id and coordinates; only VAF/genotype may
    differ between tumor and normal for a given germline SV
    """

    tumor: list[SV] = field(default_factory=list)
    normal: list[SV] = field(default_factory=list)

    @property
    def somatic_ids(self) -> set[str]:
        normal_ids = {sv.id for sv in self.normal}
        return {sv.id for sv in self.tumor if sv.id not in normal_ids}

    @property
    def germline_ids(self) -> set[str]:
        normal_ids = {sv.id for sv in self.normal}
        return {sv.id for sv in self.tumor if sv.id in normal_ids}
