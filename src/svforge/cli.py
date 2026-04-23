"""
svforge command-line interface

Subcommands:
- ``gen``        single-sample synthetic VCF
- ``gen-pair``   tumor + normal paired VCFs
- ``validate``   end-to-end sanity check
- ``bank``       list / show built-in banks
- ``callers``    list every registered writer
"""

from __future__ import annotations

import argparse
import contextlib
import logging
import secrets
import sys
from collections.abc import Sequence
from pathlib import Path

import pysam
import yaml

from svforge import __version__
from svforge.core.bank import Bank, list_builtin_banks, load_bank
from svforge.core.genome import GenomeBuild, normalize_chromosomes, validate_genome
from svforge.core.models import SV
from svforge.core.provenance import build_svforge_tags
from svforge.core.regions import RegionSet, load_bed
from svforge.core.sampler import SamplerConfig, sample, sample_pair
from svforge.io.vcf_writer import write_vcf
from svforge.validate.annotate import POOR_MAPPABILITY_FLAG, annotate_vcf, event_id
from svforge.validate.overlap import DEFAULT_TOLERANCE_BP, GnomadIndex, load_gnomad_index
from svforge.validate.report import build_report, write_tsv, write_tsv_path
from svforge.writers import CallerWriter, available_writers, get_writer

_SEED_MAX = 2**31 - 1

log = logging.getLogger("svforge")


def main(argv: Sequence[str] | None = None) -> int:
    """
    Entry point invoked by the ``svforge`` console script
    """
    parser = _build_parser()
    args = parser.parse_args(argv)
    _setup_logging(args.verbose if hasattr(args, "verbose") else 0)

    handler = _DISPATCH.get(args.command)
    if handler is None:
        parser.print_help()
        return 2
    try:
        return handler(args)
    except (ValueError, FileNotFoundError, KeyError) as exc:
        log.error("%s: %s", type(exc).__name__, exc)
        return 1


def _build_parser() -> argparse.ArgumentParser:
    p = argparse.ArgumentParser(
        prog="svforge",
        description="Synthetic VCF generator for structural variants",
    )
    p.add_argument("--version", action="version", version=f"svforge {__version__}")
    p.add_argument("-v", "--verbose", action="count", default=0, help="Increase log verbosity")

    sub = p.add_subparsers(dest="command", required=True)

    _add_gen_parser(sub)
    _add_gen_pair_parser(sub)
    _add_validate_parser(sub)
    _add_bank_parser(sub)
    _add_callers_parser(sub)

    return p


def _svtype_list(value: str) -> frozenset[str]:
    items = [v.strip().upper() for v in value.split(",") if v.strip()]
    if not items:
        raise argparse.ArgumentTypeError("svtypes must not be empty")
    return frozenset(items)


def _int_pair(value: str) -> tuple[int, int]:
    try:
        lo, hi = (int(v) for v in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected MIN,MAX ints, got {value!r}") from exc
    if hi < lo:
        raise argparse.ArgumentTypeError(f"MAX < MIN in {value!r}")
    return lo, hi


def _float_pair(value: str) -> tuple[float, float]:
    try:
        lo, hi = (float(v) for v in value.split(","))
    except ValueError as exc:
        raise argparse.ArgumentTypeError(f"expected MIN,MAX floats, got {value!r}") from exc
    if hi < lo:
        raise argparse.ArgumentTypeError(f"MAX < MIN in {value!r}")
    return lo, hi


def _chromosome_list(value: str) -> list[str]:
    items = [v.strip() for v in value.split(",") if v.strip()]
    if not items:
        raise argparse.ArgumentTypeError("--chromosomes must not be empty")
    return items


def _add_common_sampling_args(p: argparse.ArgumentParser) -> None:
    p.add_argument("--caller", required=True, choices=sorted(available_writers()))
    p.add_argument("--bank", default="default_hg38", help="Built-in bank name or path to YAML bank")
    p.add_argument("--seed", type=int, default=None)
    p.add_argument(
        "--genome",
        default="hg38",
        help="Reference build for contig lengths (V1: hg38 only)",
    )
    p.add_argument(
        "--header-template",
        type=Path,
        default=None,
        metavar="PATH",
        help="Override the bundled VCF header template with a user-supplied one",
    )
    p.add_argument("--svtypes", type=_svtype_list, default=None)
    p.add_argument("--blacklist-bed", type=Path, default=None)
    p.add_argument("--blacklist-fraction", type=float, default=0.0)
    p.add_argument("--gnomad-bed", type=Path, default=None)
    p.add_argument("--gnomad-fraction", type=float, default=0.0)
    p.add_argument("--homlen-range", type=_int_pair, default=None, metavar="MIN,MAX")
    p.add_argument("--svlen-range", type=_int_pair, default=None, metavar="MIN,MAX")
    p.add_argument("--vaf-range", type=_float_pair, default=(0.3, 1.0), metavar="MIN,MAX")
    p.add_argument(
        "--chromosomes",
        type=_chromosome_list,
        default=None,
        metavar="CHR[,CHR...]",
        help="Restrict SV generation to this subset of chromosomes (e.g. chr1,chr7 or 1,7)",
    )


def _add_gen_parser(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    p = sub.add_parser("gen", help="Generate one synthetic VCF for one sample")
    p.add_argument("--out", type=Path, required=True, help="Output .vcf / .vcf.gz / .bcf")
    p.add_argument("--n", type=int, required=True)
    p.add_argument("--sample-name", required=True)
    _add_common_sampling_args(p)


def _add_gen_pair_parser(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    p = sub.add_parser("gen-pair", help="Generate paired tumor/normal VCFs")
    p.add_argument("--out-tumor", type=Path, required=True)
    p.add_argument("--out-normal", type=Path, required=True)
    p.add_argument("--n-somatic", type=int, required=True)
    p.add_argument("--n-germline", type=int, required=True)
    p.add_argument("--tumor-sample-name", required=True)
    p.add_argument("--normal-sample-name", required=True)
    _add_common_sampling_args(p)


def _add_validate_parser(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    p = sub.add_parser("validate", help="Sanity-check a svforge-generated VCF")
    p.add_argument("--vcf", type=Path, required=True)
    p.add_argument("--blacklist-bed", type=Path, required=True)
    p.add_argument("--gnomad-vcf", type=Path, required=True)
    p.add_argument("--expected-blacklist-fraction", type=float, required=True)
    p.add_argument("--expected-gnomad-fraction", type=float, required=True)
    p.add_argument("--tolerance", type=float, default=0.20)
    p.add_argument("--report-tsv", type=Path, default=None)
    p.add_argument(
        "--annotated-vcf",
        type=Path,
        default=None,
        help="If set, write the blacklist-annotated VCF to this path",
    )


def _add_bank_parser(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    p = sub.add_parser("bank", help="Inspect built-in banks")
    bsub = p.add_subparsers(dest="bank_cmd", required=True)
    bsub.add_parser("list", help="List built-in banks")
    show = bsub.add_parser("show", help="Dump a bank as YAML")
    show.add_argument("name", nargs="?", default="default_hg38")


def _add_callers_parser(sub: argparse._SubParsersAction[argparse.ArgumentParser]) -> None:
    sub.add_parser("callers", help="List every registered writer")


def _cmd_gen(args: argparse.Namespace) -> int:
    bank = load_bank(args.bank)
    effective_seed = _effective_seed(args.seed)
    args.seed = effective_seed
    cfg = _make_sampler_config(args, bank)
    svs = sample(bank, cfg)
    writer = get_writer(args.caller)
    provenance = build_svforge_tags(caller=args.caller, seed=effective_seed, argv=sys.argv)
    _write_sample_vcf(
        writer,
        svs,
        args.sample_name,
        args.out,
        args.genome,
        provenance,
        header_template_override=args.header_template,
    )
    log.info("Wrote %d SVs to %s (seed=%d)", len(svs), args.out, effective_seed)
    return 0


def _cmd_gen_pair(args: argparse.Namespace) -> int:
    bank = load_bank(args.bank)
    effective_seed = _effective_seed(args.seed)
    args.seed = effective_seed
    cfg = _make_sampler_config(args, bank, n_override=0)
    pair = sample_pair(bank, args.n_somatic, args.n_germline, cfg)
    writer = get_writer(args.caller)
    provenance = build_svforge_tags(caller=args.caller, seed=effective_seed, argv=sys.argv)
    _write_sample_vcf(
        writer,
        pair.tumor,
        args.tumor_sample_name,
        args.out_tumor,
        args.genome,
        provenance,
        header_template_override=args.header_template,
    )
    _write_sample_vcf(
        writer,
        pair.normal,
        args.normal_sample_name,
        args.out_normal,
        args.genome,
        provenance,
        header_template_override=args.header_template,
    )
    log.info(
        "Tumor: %d SVs (%d somatic + %d germline), Normal: %d SVs -> %s, %s (seed=%d)",
        len(pair.tumor),
        len(pair.somatic_ids),
        len(pair.germline_ids),
        len(pair.normal),
        args.out_tumor,
        args.out_normal,
        effective_seed,
    )
    return 0


def _effective_seed(seed: int | None) -> int:
    """
    Pick the effective seed: the user's if supplied, otherwise a fresh one

    The result is always logged in the ``##svforgeSeed`` provenance tag so a
    run is always reproducible from its output VCF
    """
    if seed is not None:
        return int(seed)
    return secrets.randbelow(_SEED_MAX)


def _cmd_validate(args: argparse.Namespace) -> int:
    blacklist = load_bed(args.blacklist_bed)
    gnomad = load_gnomad_index(args.gnomad_vcf)

    annotated_path = args.annotated_vcf
    if annotated_path is None:
        annotated_path = args.vcf.with_suffix(args.vcf.suffix + ".annotated.vcf")
    total, flagged = annotate_vcf(args.vcf, annotated_path, blacklist)
    gnomad_hits = _count_gnomad_hits(annotated_path, gnomad)

    checks = build_report(
        total_svs=total,
        n_gnomad_overlaps=gnomad_hits,
        n_blacklist_flags=flagged,
        expected_gnomad_fraction=args.expected_gnomad_fraction,
        expected_blacklist_fraction=args.expected_blacklist_fraction,
        tolerance=args.tolerance,
    )

    if args.report_tsv:
        write_tsv_path(checks, args.report_tsv)
    else:
        write_tsv(checks, sys.stdout)

    if args.annotated_vcf is None:
        _maybe_unlink(annotated_path)

    return 0 if all(c.passed for c in checks) else 1


def _cmd_bank(args: argparse.Namespace) -> int:
    if args.bank_cmd == "list":
        for name in list_builtin_banks():
            print(name)
        return 0
    if args.bank_cmd == "show":
        bank = load_bank(args.name)
        yaml.safe_dump(_bank_to_dict(bank), sys.stdout, sort_keys=False)
        return 0
    return 2  # pragma: no cover -- argparse enforces required subcommand


def _cmd_callers(_args: argparse.Namespace) -> int:
    for name in available_writers():
        print(name)
    return 0


_DISPATCH = {
    "gen": _cmd_gen,
    "gen-pair": _cmd_gen_pair,
    "validate": _cmd_validate,
    "bank": _cmd_bank,
    "callers": _cmd_callers,
}


def _make_sampler_config(
    args: argparse.Namespace,
    bank: Bank,
    n_override: int | None = None,
) -> SamplerConfig:
    genome: GenomeBuild = validate_genome(args.genome)
    blacklist = load_bed(args.blacklist_bed) if args.blacklist_bed else None
    gnomad = _maybe_load_gnomad(args.gnomad_bed)

    chroms: frozenset[str] | None = None
    raw_chroms = getattr(args, "chromosomes", None)
    if raw_chroms:
        chroms = frozenset(normalize_chromosomes(raw_chroms, genome))

    return SamplerConfig(
        n=n_override if n_override is not None else args.n,
        genome=genome,
        svtypes=args.svtypes,
        svlen_range=args.svlen_range,
        homlen_range=args.homlen_range,
        vaf_range=args.vaf_range,
        blacklist=blacklist,
        blacklist_fraction=args.blacklist_fraction,
        gnomad=gnomad,
        gnomad_fraction=args.gnomad_fraction,
        seed=args.seed,
        chroms=chroms,
    )


def _maybe_load_gnomad(path: Path | None) -> RegionSet | None:
    """
    Accept either a BED (preferred) or a VCF and build a :class:`RegionSet`
    """
    if path is None:
        return None
    lower = path.name.lower()
    if lower.endswith(".bed") or lower.endswith(".bed.gz"):
        return load_bed(path)
    if lower.endswith((".vcf", ".vcf.gz", ".bcf")):
        return _region_set_from_vcf(path)
    return load_bed(path)


def _region_set_from_vcf(path: Path) -> RegionSet:
    rs = RegionSet()
    with pysam.VariantFile(str(path)) as vf:
        for rec in vf:
            chrom = str(rec.chrom)
            start = int(rec.pos) - 1
            end = int(rec.stop)
            if end <= start:
                end = start + 1
            rs.add(chrom, start, end)
    return rs


def _write_sample_vcf(
    writer: CallerWriter,
    svs: list[SV],
    sample_name: str,
    out_path: Path,
    genome: GenomeBuild,
    provenance_tags: Sequence[str],
    *,
    header_template_override: Path | None = None,
) -> None:
    header = writer.header_lines(
        sample_name,
        genome=genome,
        provenance_tags=provenance_tags,
        template_override=header_template_override,
    )
    records = writer.format_records(svs, sample_name)
    write_vcf(out_path, header, records)


def _count_gnomad_hits(vcf_path: Path, gnomad_index: GnomadIndex) -> int:
    """
    Event-level count of gnomAD overlaps (mate pairs count once)
    """
    hit_events: set[tuple[str, str]] = set()
    with pysam.VariantFile(str(vcf_path)) as vf:
        for rec in vf:
            sv = _pysam_to_sv(rec)
            if sv is None:
                continue
            if gnomad_index.overlaps(sv, tolerance=DEFAULT_TOLERANCE_BP):
                hit_events.add(event_id(rec))
    return len(hit_events)


def _pysam_to_sv(rec: pysam.VariantRecord) -> SV | None:
    svtype_raw = _safe_info(rec, "SVTYPE")
    if isinstance(svtype_raw, tuple):
        svtype_raw = svtype_raw[0] if svtype_raw else None
    if svtype_raw is None:
        return None
    svtype = str(svtype_raw).upper()

    chrom = str(rec.chrom)
    pos = int(rec.pos)
    end = int(rec.stop)

    chr2 = _safe_info(rec, "CHR2")
    if isinstance(chr2, tuple):
        chr2 = chr2[0] if chr2 else None
    pos2 = _safe_info(rec, "POS2")
    if isinstance(pos2, tuple):
        pos2 = pos2[0] if pos2 else None

    mate_chrom: str | None = None
    mate_pos: int | None = None
    if svtype == "BND":
        if chr2 is not None and pos2 is not None:
            mate_chrom = str(chr2)
            mate_pos = _to_int(pos2)
        else:
            alt = str(rec.alts[0]) if rec.alts else ""
            parsed = _parse_bnd_alt(alt)
            if parsed is None:
                return None
            mate_chrom, mate_pos = parsed
        end = pos

    svlen_raw = _safe_info(rec, "SVLEN")
    if isinstance(svlen_raw, tuple):
        svlen_raw = svlen_raw[0] if svlen_raw else None
    svlen = abs(_to_int(svlen_raw)) if svlen_raw is not None else max(1, end - pos)
    if svtype == "BND":
        svlen = 1

    flag_raw = _safe_info(rec, POOR_MAPPABILITY_FLAG)
    info_extra: dict[str, str] = {}
    if bool(flag_raw):
        info_extra[POOR_MAPPABILITY_FLAG] = ""

    return SV(
        id=str(rec.id) if rec.id else f"rec_{pos}",
        svtype=svtype,
        chrom=chrom,
        pos=pos,
        end=end,
        svlen=svlen,
        mate_chrom=mate_chrom,
        mate_pos=mate_pos,
        info_extra=info_extra,
    )


def _safe_info(rec: pysam.VariantRecord, key: str) -> object:
    """
    Return INFO value for ``key`` or None if the header lacks the declaration
    """
    if key not in rec.header.info:
        return None
    try:
        return rec.info.get(key)
    except (KeyError, ValueError):
        return None


def _to_int(value: object) -> int:
    """
    Coerce a pysam INFO scalar (already proven non-None) to ``int``
    """
    if isinstance(value, (int, str, float)):
        return int(value)
    return int(str(value))


def _parse_bnd_alt(alt: str) -> tuple[str, int] | None:
    for opener, closer in (("[", "["), ("]", "]")):
        if opener in alt and closer in alt:
            try:
                inside = alt.split(opener, 1)[1].split(closer, 1)[0]
                chrom, pos_str = inside.split(":", 1)
                return chrom, int(pos_str)
            except (IndexError, ValueError):
                continue
    return None


def _bank_to_dict(bank: Bank) -> dict[str, object]:
    return {
        "name": bank.name,
        "genome": bank.genome,
        "templates": [
            {
                "svtype": t.svtype,
                "svlen": [t.svlen_min, t.svlen_max],
                "homlen": [t.homlen_min, t.homlen_max],
                "weight": t.weight,
                "chroms": list(t.chroms),
            }
            for t in bank.templates
        ],
    }


def _maybe_unlink(path: Path) -> None:
    with contextlib.suppress(FileNotFoundError):
        path.unlink()


def _setup_logging(verbose: int) -> None:
    level = logging.WARNING
    if verbose == 1:
        level = logging.INFO
    elif verbose >= 2:
        level = logging.DEBUG
    logging.basicConfig(
        level=level,
        format="%(asctime)s %(levelname)s %(name)s: %(message)s",
    )
