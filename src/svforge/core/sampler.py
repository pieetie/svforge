"""
Sample SVs from a :class:`Bank` with reproducible RNG

Knobs:

- ``svtypes``: restrict to a subset of svtypes (e.g. ``{"DEL", "DUP"}``)
- ``svlen_range`` / ``homlen_range``: hard CLI overrides, clamp template ranges
- ``vaf_range``: uniform VAF draw
- ``gnomad_fraction`` / ``blacklist_fraction``: fraction of the requested SV
  count drawn directly from the bundled mini catalogs (Phase C). The
  catalog rows (CHROM/POS/END) are reproduced verbatim so downstream
  pipelines match them exactly. Sampling is without replacement when
  possible; with replacement otherwise

``sample_pair`` composes ``sample`` twice to produce a consistent tumor/normal
SV set where germline events are duplicated across both samples with the same
coordinates but possibly different VAFs
"""

from __future__ import annotations

import logging
import random
import uuid
from collections.abc import Sequence
from dataclasses import dataclass
from typing import Literal, TypeVar

from svforge.core.bank import Bank, SVTemplate
from svforge.core.genome import GenomeBuild, get_contigs
from svforge.core.injection_catalogs import (
    BlacklistEntry,
    GnomadEntry,
    load_blacklist_catalog,
    load_gnomad_catalog,
)
from svforge.core.models import SV, SVPair

_CHROM_SAMPLE_MARGIN = 10_000
_BLACKLIST_SVTYPES = ("DEL", "DUP", "INV")

_T = TypeVar("_T")

log = logging.getLogger(__name__)


@dataclass(slots=True)
class SamplerConfig:
    """
    Knobs that shape a sampling run
    """

    n: int
    genome: GenomeBuild = "hg38"
    svtypes: frozenset[str] | None = None
    svlen_range: tuple[int, int] | None = None
    homlen_range: tuple[int, int] | None = None
    vaf_range: tuple[float, float] = (0.3, 1.0)
    blacklist_fraction: float = 0.0
    gnomad_fraction: float = 0.0
    seed: int | None = None
    origin: Literal["somatic", "germline"] = "somatic"
    id_prefix: str = "svforge"
    chroms: frozenset[str] | None = None

    def __post_init__(self) -> None:
        if self.n < 0:
            raise ValueError(f"n must be >= 0, got {self.n}")
        if not 0.0 <= self.blacklist_fraction <= 1.0:
            raise ValueError("blacklist_fraction must be in [0, 1]")
        if not 0.0 <= self.gnomad_fraction <= 1.0:
            raise ValueError("gnomad_fraction must be in [0, 1]")
        total = self.blacklist_fraction + self.gnomad_fraction
        if total > 1.0 + 1e-9:
            raise ValueError(
                f"blacklist_fraction + gnomad_fraction must be <= 1.0, got {total:.3f}"
            )
        if self.vaf_range[0] > self.vaf_range[1] or not (
            0.0 <= self.vaf_range[0] <= 1.0 and 0.0 <= self.vaf_range[1] <= 1.0
        ):
            raise ValueError(f"invalid vaf_range: {self.vaf_range}")


def sample(bank: Bank, cfg: SamplerConfig) -> list[SV]:
    """
    Draw ``cfg.n`` SVs from ``bank`` respecting ``cfg`` constraints
    """
    rng = random.Random(cfg.seed)
    contigs = get_contigs(cfg.genome)
    if cfg.chroms is not None:
        contigs = {c: n for c, n in contigs.items() if c in cfg.chroms}

    templates = _filter_templates(bank.templates, cfg.svtypes, cfg.chroms)
    if not templates:
        if cfg.chroms is not None:
            raise ValueError(
                f"No SV in bank matches the requested chromosomes: "
                f"{sorted(cfg.chroms)}. Enrich the bank or broaden the filter"
            )
        raise ValueError(f"No template in bank {bank.name!r} matches svtypes={cfg.svtypes}")

    weights = [t.weight for t in templates]

    n_blacklist = round(cfg.n * cfg.blacklist_fraction)
    n_gnomad = round(cfg.n * cfg.gnomad_fraction)
    if n_blacklist + n_gnomad > cfg.n:
        n_gnomad = max(0, cfg.n - n_blacklist)
    n_bank = max(0, cfg.n - n_blacklist - n_gnomad)

    svs: list[SV] = []
    svs.extend(_draw_bank(n_bank, templates, weights, contigs, cfg, rng))
    if n_gnomad:
        svs.extend(_draw_gnomad(n_gnomad, cfg, rng))
    if n_blacklist:
        svs.extend(_draw_blacklist(n_blacklist, cfg, rng))

    rng.shuffle(svs)
    svs.sort(key=lambda sv: (_chrom_key(sv.chrom), sv.pos))
    return svs


def sample_pair(
    bank: Bank,
    n_somatic: int,
    n_germline: int,
    cfg: SamplerConfig,
    germline_vaf_range: tuple[float, float] = (0.45, 1.0),
) -> SVPair:
    """
    Draw a coherent tumor + normal SV pair

    Germline SVs appear in both VCFs with the same id/coords. Their VAF in
    the normal is drawn from ``germline_vaf_range`` (~0.5 het / ~1.0 hom)
    while their VAF in the tumor stays within ``cfg.vaf_range``
    """
    if n_somatic < 0 or n_germline < 0:
        raise ValueError("n_somatic and n_germline must be >= 0")

    base_seed = cfg.seed

    somatic_cfg = _replace(cfg, n=n_somatic, origin="somatic", seed=base_seed)
    somatic = sample(bank, somatic_cfg) if n_somatic else []

    germline_seed = None if base_seed is None else base_seed + 1
    germline_cfg = _replace(cfg, n=n_germline, origin="germline", seed=germline_seed)
    germline_tumor = sample(bank, germline_cfg) if n_germline else []

    rng_normal = random.Random(None if base_seed is None else base_seed + 2)
    germline_normal: list[SV] = []
    for sv in germline_tumor:
        vaf = rng_normal.uniform(*germline_vaf_range)
        genotype = "1/1" if vaf >= 0.9 else "0/1"
        germline_normal.append(
            SV(
                id=sv.id,
                svtype=sv.svtype,
                chrom=sv.chrom,
                pos=sv.pos,
                end=sv.end,
                svlen=sv.svlen,
                mate_chrom=sv.mate_chrom,
                mate_pos=sv.mate_pos,
                strands=sv.strands,
                homlen=sv.homlen,
                homseq=sv.homseq,
                vaf=vaf,
                genotype=genotype,
                ref_base=sv.ref_base,
                ins_seq=sv.ins_seq,
                filter=sv.filter,
                origin="germline",
                source=sv.source,
                info_extra=dict(sv.info_extra),
            )
        )

    tumor = sorted(somatic + germline_tumor, key=lambda sv: (_chrom_key(sv.chrom), sv.pos))
    normal = sorted(germline_normal, key=lambda sv: (_chrom_key(sv.chrom), sv.pos))
    return SVPair(tumor=tumor, normal=normal)


def _filter_templates(
    templates: Sequence[SVTemplate],
    svtypes: frozenset[str] | None,
    chroms: frozenset[str] | None,
) -> list[SVTemplate]:
    out: list[SVTemplate] = []
    for t in templates:
        if svtypes is not None and t.svtype not in svtypes:
            continue
        if chroms is not None and t.chroms and not any(c in chroms for c in t.chroms):
            continue
        out.append(t)
    return out


def _draw_bank(
    count: int,
    templates: Sequence[SVTemplate],
    weights: Sequence[float],
    contigs: dict[str, int],
    cfg: SamplerConfig,
    rng: random.Random,
) -> list[SV]:
    out: list[SV] = []
    for _ in range(count):
        template = rng.choices(list(templates), weights=list(weights), k=1)[0]
        sv = _materialize(template, contigs, cfg, rng)
        out.append(sv)
    return out


def _draw_gnomad(count: int, cfg: SamplerConfig, rng: random.Random) -> list[SV]:
    catalog = list(load_gnomad_catalog())
    eligible = _filter_gnomad_catalog(catalog, cfg)
    if not eligible:
        raise ValueError(
            "No gnomAD catalog entry matches the active filter "
            f"(svtypes={cfg.svtypes}, chroms={sorted(cfg.chroms) if cfg.chroms else None})"
        )
    picks = _sample_without_replacement(eligible, count, rng)
    return [_gnomad_entry_to_sv(entry, cfg, rng) for entry in picks]


def _draw_blacklist(count: int, cfg: SamplerConfig, rng: random.Random) -> list[SV]:
    catalog = list(load_blacklist_catalog())
    eligible = _filter_blacklist_catalog(catalog, cfg)
    if not eligible:
        raise ValueError(
            "No blacklist catalog entry matches the active filter "
            f"(chroms={sorted(cfg.chroms) if cfg.chroms else None})"
        )
    picks = _sample_without_replacement(eligible, count, rng)
    return [_blacklist_entry_to_sv(entry, cfg, rng) for entry in picks]


def _filter_gnomad_catalog(catalog: list[GnomadEntry], cfg: SamplerConfig) -> list[GnomadEntry]:
    out: list[GnomadEntry] = []
    for entry in catalog:
        if cfg.svtypes is not None and entry.svtype not in cfg.svtypes:
            continue
        if cfg.chroms is not None and entry.chrom not in cfg.chroms:
            continue
        if cfg.chroms is not None and entry.end_chrom not in cfg.chroms:
            continue
        out.append(entry)
    return out


def _filter_blacklist_catalog(
    catalog: list[BlacklistEntry], cfg: SamplerConfig
) -> list[BlacklistEntry]:
    out: list[BlacklistEntry] = []
    eligible_svtypes = tuple(cfg.svtypes) if cfg.svtypes is not None else _BLACKLIST_SVTYPES
    if not any(t in _BLACKLIST_SVTYPES for t in eligible_svtypes):
        raise ValueError(
            "Blacklist injection requires at least one of DEL/DUP/INV in svtypes; "
            f"got {sorted(eligible_svtypes)}"
        )
    for entry in catalog:
        if cfg.chroms is not None and entry.chrom not in cfg.chroms:
            continue
        out.append(entry)
    return out


def _sample_without_replacement(pool: Sequence[_T], count: int, rng: random.Random) -> list[_T]:
    if count <= len(pool):
        return rng.sample(list(pool), count)
    log.warning(
        "Requested %d injections from catalog of size %d; falling back to sampling "
        "with replacement",
        count,
        len(pool),
    )
    return [rng.choice(list(pool)) for _ in range(count)]


def _gnomad_entry_to_sv(entry: GnomadEntry, cfg: SamplerConfig, rng: random.Random) -> SV:
    vaf = rng.uniform(*cfg.vaf_range)
    genotype = "1/1" if vaf >= 0.9 else "0/1"
    sv_id = f"{cfg.id_prefix}_{uuid.UUID(int=rng.getrandbits(128)).hex[:10]}"

    if entry.svtype == "BND":
        return SV(
            id=sv_id,
            svtype="BND",
            chrom=entry.chrom,
            pos=entry.pos,
            end=entry.pos,
            svlen=1,
            mate_chrom=entry.end_chrom,
            mate_pos=entry.end,
            strands="+-",
            vaf=vaf,
            genotype=genotype,
            origin=cfg.origin,
            source="gnomad",
            info_extra={"SVFORGE_SOURCE_ID": entry.source_id},
        )

    svlen = max(1, entry.end - entry.pos)
    return SV(
        id=sv_id,
        svtype=entry.svtype,
        chrom=entry.chrom,
        pos=entry.pos,
        end=entry.end,
        svlen=svlen,
        strands=_default_strands(entry.svtype),
        vaf=vaf,
        genotype=genotype,
        origin=cfg.origin,
        source="gnomad",
        info_extra={"SVFORGE_SOURCE_ID": entry.source_id},
    )


def _blacklist_entry_to_sv(entry: BlacklistEntry, cfg: SamplerConfig, rng: random.Random) -> SV:
    vaf = rng.uniform(*cfg.vaf_range)
    genotype = "1/1" if vaf >= 0.9 else "0/1"
    sv_id = f"{cfg.id_prefix}_{uuid.UUID(int=rng.getrandbits(128)).hex[:10]}"

    svtype_pool = (
        [t for t in _BLACKLIST_SVTYPES if cfg.svtypes is None or t in cfg.svtypes]
        if cfg.svtypes is not None
        else list(_BLACKLIST_SVTYPES)
    )
    svtype = rng.choice(svtype_pool)
    svlen = max(1, entry.end - entry.pos)
    return SV(
        id=sv_id,
        svtype=svtype,
        chrom=entry.chrom,
        pos=entry.pos,
        end=entry.end,
        svlen=svlen,
        strands=_default_strands(svtype),
        vaf=vaf,
        genotype=genotype,
        origin=cfg.origin,
        source="blacklist",
        info_extra={"SVFORGE_SOURCE_ID": entry.source_id},
    )


def _materialize(
    template: SVTemplate,
    contigs: dict[str, int],
    cfg: SamplerConfig,
    rng: random.Random,
) -> SV:
    allowed_chroms = _template_chrom_pool(template, contigs, cfg)
    if not allowed_chroms:
        raise ValueError(
            f"Template {template.svtype!r} has no chromosome compatible with the "
            f"current filter (chroms={sorted(cfg.chroms) if cfg.chroms else None})"
        )
    chrom = rng.choice(allowed_chroms)
    chrom_len = contigs[chrom]

    svlen_min, svlen_max = _clamp_range((template.svlen_min, template.svlen_max), cfg.svlen_range)
    homlen_min, homlen_max = _clamp_range(
        (template.homlen_min, template.homlen_max), cfg.homlen_range
    )

    svlen = 1 if template.svtype == "BND" else rng.randint(svlen_min, svlen_max)

    homlen = rng.randint(homlen_min, homlen_max)
    vaf = rng.uniform(*cfg.vaf_range)
    genotype = "1/1" if vaf >= 0.9 else "0/1"
    sv_id = f"{cfg.id_prefix}_{uuid.UUID(int=rng.getrandbits(128)).hex[:10]}"

    if template.svtype == "BND":
        mate_candidates = _bnd_mate_candidates(chrom, contigs, cfg)
        mate_chrom = rng.choice(mate_candidates)
        pos = rng.randint(
            _CHROM_SAMPLE_MARGIN, max(_CHROM_SAMPLE_MARGIN + 1, chrom_len - _CHROM_SAMPLE_MARGIN)
        )
        mate_pos = rng.randint(
            _CHROM_SAMPLE_MARGIN,
            max(_CHROM_SAMPLE_MARGIN + 1, contigs[mate_chrom] - _CHROM_SAMPLE_MARGIN),
        )
        strands = rng.choice(("+-", "-+", "++", "--"))
        return SV(
            id=sv_id,
            svtype="BND",
            chrom=chrom,
            pos=pos,
            end=pos,
            svlen=svlen,
            mate_chrom=mate_chrom,
            mate_pos=mate_pos,
            strands=strands,
            homlen=homlen,
            vaf=vaf,
            genotype=genotype,
            origin=cfg.origin,
            source="bank",
        )

    max_pos = chrom_len - svlen - _CHROM_SAMPLE_MARGIN
    if max_pos <= _CHROM_SAMPLE_MARGIN:
        pos = _CHROM_SAMPLE_MARGIN
    else:
        pos = rng.randint(_CHROM_SAMPLE_MARGIN, max_pos)
    end = pos + svlen

    strands = _default_strands(template.svtype)
    return SV(
        id=sv_id,
        svtype=template.svtype,
        chrom=chrom,
        pos=pos,
        end=end,
        svlen=svlen,
        strands=strands,
        homlen=homlen,
        vaf=vaf,
        genotype=genotype,
        origin=cfg.origin,
        source="bank",
    )


def _template_chrom_pool(
    template: SVTemplate, contigs: dict[str, int], cfg: SamplerConfig
) -> list[str]:
    base = template.chroms if template.chroms else tuple(contigs.keys())
    pool = [c for c in base if c in contigs]
    if cfg.chroms is not None:
        pool = [c for c in pool if c in cfg.chroms]
    return pool


def _bnd_mate_candidates(chrom: str, contigs: dict[str, int], cfg: SamplerConfig) -> list[str]:
    """
    Return chroms eligible as a BND mate: different from ``chrom`` when
    possible, restricted to the user filter when one is active
    """
    universe = list(contigs.keys())
    if cfg.chroms is not None:
        universe = [c for c in universe if c in cfg.chroms]
    others = [c for c in universe if c != chrom]
    if others:
        return others
    return universe or [chrom]


def _default_strands(svtype: str) -> str:
    match svtype:
        case "DEL":
            return "+-"
        case "DUP":
            return "-+"
        case "INV":
            return "++"
        case "INS":
            return "+-"
        case _:
            return "+-"


def _clamp_range(
    template_range: tuple[int, int], override: tuple[int, int] | None
) -> tuple[int, int]:
    lo, hi = template_range
    if override is None:
        return lo, hi
    o_lo, o_hi = override
    new_lo = max(lo, o_lo)
    new_hi = min(hi, o_hi)
    if new_hi < new_lo:
        new_lo = new_hi = max(o_lo, min(o_hi, lo))
    return new_lo, new_hi


def _replace(cfg: SamplerConfig, **kw: object) -> SamplerConfig:
    data: dict[str, object] = {
        "n": cfg.n,
        "genome": cfg.genome,
        "svtypes": cfg.svtypes,
        "svlen_range": cfg.svlen_range,
        "homlen_range": cfg.homlen_range,
        "vaf_range": cfg.vaf_range,
        "blacklist_fraction": cfg.blacklist_fraction,
        "gnomad_fraction": cfg.gnomad_fraction,
        "seed": cfg.seed,
        "origin": cfg.origin,
        "id_prefix": cfg.id_prefix,
        "chroms": cfg.chroms,
    }
    for k, v in kw.items():
        data[k] = v
    return SamplerConfig(**data)  # type: ignore[arg-type]


def _chrom_key(chrom: str) -> tuple[int, str]:
    stripped = chrom.removeprefix("chr")
    if stripped.isdigit():
        return (int(stripped), "")
    order = {"X": 23, "Y": 24, "M": 25, "MT": 25}.get(stripped, 99)
    return (order, stripped)
