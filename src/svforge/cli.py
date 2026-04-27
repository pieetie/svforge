"""
svforge command-line interface

Subcommands:

- ``gen``        single-sample synthetic VCF
- ``gen-pair``   tumor + normal paired VCFs
- ``validate``   self-consistency check against bundled injection catalogs
- ``bank``       list / show built-in banks
- ``callers``    list every registered writer
"""

from __future__ import annotations

import argparse
import logging
import secrets
import sys
from collections.abc import Sequence
from pathlib import Path

import yaml

from svforge import __version__
from svforge.core.bank import Bank, list_builtin_banks, load_bank
from svforge.core.genome import GenomeBuild, normalize_chromosomes, validate_genome
from svforge.core.models import SV
from svforge.core.provenance import build_svforge_tags
from svforge.core.sampler import SamplerConfig, sample, sample_pair
from svforge.io.vcf_writer import write_vcf
from svforge.validate.annotate import ValidationResult, validate_vcf
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
    p.add_argument(
        "--blacklist-fraction",
        type=float,
        default=0.0,
        help="Fraction of SVs drawn from the bundled ENCODE blacklist mini catalog",
    )
    p.add_argument(
        "--gnomad-fraction",
        type=float,
        default=0.0,
        help="Fraction of SVs drawn from the bundled gnomAD mini catalog",
    )
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
    p = sub.add_parser(
        "validate",
        help="Self-consistency check: every SVFORGE_SOURCE-tagged record "
        "matches the bundled mini catalog",
    )
    p.add_argument("--vcf", type=Path, required=True)
    p.add_argument(
        "--report-tsv",
        type=Path,
        default=None,
        help="Write the divergence report to this TSV (default: stdout)",
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
    result = validate_vcf(args.vcf)
    text = _format_validation_report(result)
    if args.report_tsv:
        args.report_tsv.parent.mkdir(parents=True, exist_ok=True)
        args.report_tsv.write_text(text, encoding="utf-8")
    else:
        sys.stdout.write(text)
    return 0 if result.passed else 1


def _format_validation_report(result: ValidationResult) -> str:
    lines = ["metric\ttotal\tmatched\tstatus"]
    for metric, total, matched in (
        ("gnomad_self_consistency", result.gnomad_total, result.gnomad_matched),
        ("blacklist_self_consistency", result.blacklist_total, result.blacklist_matched),
    ):
        status = "PASS" if total == matched else "FAIL"
        lines.append(f"{metric}\t{total}\t{matched}\t{status}")
    if result.divergences:
        lines.append("")
        lines.append("# divergences (record_id, source, chrom, pos, end)")
        for div in result.divergences:
            lines.append(f"{div.record_id}\t{div.source}\t{div.chrom}\t{div.pos}\t{div.end}")
    return "\n".join(lines) + "\n"


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
        blacklist_fraction=args.blacklist_fraction,
        gnomad_fraction=args.gnomad_fraction,
        seed=args.seed,
        chroms=chroms,
    )


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
    records = writer.format_records_sorted(svs, sample_name, header)
    write_vcf(out_path, header, records)


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
