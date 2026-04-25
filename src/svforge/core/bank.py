"""
SV banks: YAML-backed catalogue of weighted SV templates

A bank is a set of :class:`SVTemplate` entries. Each template describes a
family of SVs (an svtype plus coordinate / length ranges) from which the
sampler draws concrete SVs. Banks are pure metadata: no genomic sequence
ever ships with the package
"""

from __future__ import annotations

from dataclasses import dataclass, field
from importlib import resources
from pathlib import Path
from typing import Any

import yaml

from svforge.core.genome import GenomeBuild, validate_genome

BUILTIN_BANKS_PACKAGE = "svforge.data.banks"
BUILTIN_BANK_SUFFIX = ".yaml"


@dataclass(slots=True, frozen=True)
class SVTemplate:
    """
    A weighted template describing how to draw one SV
    """

    svtype: str
    svlen_min: int
    svlen_max: int
    homlen_min: int = 0
    homlen_max: int = 0
    weight: float = 1.0
    chroms: tuple[str, ...] = ()

    def __post_init__(self) -> None:
        if self.svlen_min < 1:
            raise ValueError(f"svlen_min must be >= 1, got {self.svlen_min}")
        if self.svlen_max < self.svlen_min:
            raise ValueError(
                f"svlen_max ({self.svlen_max}) must be >= svlen_min ({self.svlen_min})"
            )
        if self.homlen_min < 0 or self.homlen_max < self.homlen_min:
            raise ValueError(f"invalid homlen range [{self.homlen_min}, {self.homlen_max}]")
        if self.weight <= 0:
            raise ValueError(f"weight must be > 0, got {self.weight}")


@dataclass(slots=True)
class Bank:
    """
    A named collection of SV templates bound to a genome build
    """

    name: str
    genome: GenomeBuild
    templates: list[SVTemplate] = field(default_factory=list)

    def __post_init__(self) -> None:
        if not self.templates:
            raise ValueError(f"Bank {self.name!r} has no templates")

    def by_svtype(self, svtype: str) -> list[SVTemplate]:
        """
        Return all templates that produce ``svtype`` SVs
        """
        return [t for t in self.templates if t.svtype == svtype]

    def svtypes(self) -> set[str]:
        """
        Distinct svtypes covered by this bank
        """
        return {t.svtype for t in self.templates}


def _template_from_dict(raw: dict[str, Any]) -> SVTemplate:
    try:
        svtype = str(raw["svtype"])
        svlen = raw["svlen"]
        svlen_min, svlen_max = int(svlen[0]), int(svlen[1])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"Invalid template entry: {raw!r} ({exc})") from exc

    homlen = raw.get("homlen", [0, 0])
    homlen_min, homlen_max = int(homlen[0]), int(homlen[1])
    weight = float(raw.get("weight", 1.0))
    chroms = tuple(str(c) for c in raw.get("chroms", []))

    return SVTemplate(
        svtype=svtype,
        svlen_min=svlen_min,
        svlen_max=svlen_max,
        homlen_min=homlen_min,
        homlen_max=homlen_max,
        weight=weight,
        chroms=chroms,
    )


def load_bank(source: str | Path) -> Bank:
    """
    Load a bank by name (built-in) or by filesystem path

    If ``source`` names an existing file, it is parsed as YAML. Otherwise
    ``source`` is looked up among the banks packaged under
    ``svforge.data.banks``
    """
    text = _read_bank_text(source)
    data = yaml.safe_load(text)
    if not isinstance(data, dict):
        raise ValueError(f"Bank source {source!r} must be a YAML mapping at top level")

    name = str(data.get("name", str(source)))
    genome_raw = str(data.get("genome", "hg38"))
    genome: GenomeBuild = validate_genome(genome_raw)

    raw_templates = data.get("templates")
    if not isinstance(raw_templates, list) or not raw_templates:
        raise ValueError(f"Bank {name!r} must define a non-empty 'templates' list")

    templates = [_template_from_dict(t) for t in raw_templates]
    return Bank(name=name, genome=genome, templates=templates)


def _read_bank_text(source: str | Path) -> str:
    path = Path(source)
    if path.exists():
        return path.read_text(encoding="utf-8")

    resource = resources.files(BUILTIN_BANKS_PACKAGE).joinpath(f"{source}{BUILTIN_BANK_SUFFIX}")
    if not resource.is_file():
        raise FileNotFoundError(
            f"Bank {source!r} not found: neither a file path nor a built-in bank. "
            f"Available: {sorted(list_builtin_banks())}"
        )
    return resource.read_text(encoding="utf-8")


def list_builtin_banks() -> list[str]:
    """
    Return the names of banks shipped with the package
    """
    root = resources.files(BUILTIN_BANKS_PACKAGE)
    names: list[str] = []
    for entry in root.iterdir():
        entry_name = entry.name
        if entry_name.endswith(BUILTIN_BANK_SUFFIX):
            names.append(entry_name.removesuffix(BUILTIN_BANK_SUFFIX))
    return sorted(names)
