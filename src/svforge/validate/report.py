"""
Fraction report: compare observed vs expected blacklist/gnomAD overlap rates
"""

from __future__ import annotations

from dataclasses import dataclass
from pathlib import Path
from typing import TextIO


@dataclass(slots=True, frozen=True)
class FractionCheck:
    """
    One row of the validation report
    """

    metric: str
    expected: float
    observed: float
    tolerance: float
    total: int
    count: int

    @property
    def passed(self) -> bool:
        if self.expected == 0.0:
            return self.observed == 0.0
        return abs(self.observed - self.expected) / self.expected <= self.tolerance

    def as_row(self) -> list[str]:
        return [
            self.metric,
            _fmt(self.expected),
            _fmt(self.observed),
            _fmt(self.tolerance),
            str(self.total),
            str(self.count),
            "PASS" if self.passed else "FAIL",
        ]


HEADER = ["metric", "expected", "observed", "tolerance", "total", "count", "status"]


def build_report(
    total_svs: int,
    n_gnomad_overlaps: int,
    n_blacklist_flags: int,
    expected_gnomad_fraction: float,
    expected_blacklist_fraction: float,
    tolerance: float,
) -> list[FractionCheck]:
    """
    Compose the full set of checks for the validate subcommand
    """
    obs_gnomad = 0.0 if total_svs == 0 else n_gnomad_overlaps / total_svs
    obs_blacklist = 0.0 if total_svs == 0 else n_blacklist_flags / total_svs
    return [
        FractionCheck(
            metric="gnomad_overlap_fraction",
            expected=expected_gnomad_fraction,
            observed=obs_gnomad,
            tolerance=tolerance,
            total=total_svs,
            count=n_gnomad_overlaps,
        ),
        FractionCheck(
            metric="blacklist_flag_fraction",
            expected=expected_blacklist_fraction,
            observed=obs_blacklist,
            tolerance=tolerance,
            total=total_svs,
            count=n_blacklist_flags,
        ),
    ]


def write_tsv(checks: list[FractionCheck], out: TextIO) -> None:
    """
    Write a TSV report with a header row and one row per check
    """
    out.write("\t".join(HEADER))
    out.write("\n")
    for chk in checks:
        out.write("\t".join(chk.as_row()))
        out.write("\n")


def write_tsv_path(checks: list[FractionCheck], path: str | Path) -> Path:
    """
    Write the TSV report to ``path`` (parent dir created if needed)
    """
    p = Path(path)
    p.parent.mkdir(parents=True, exist_ok=True)
    with p.open("w", encoding="utf-8") as fh:
        write_tsv(checks, fh)
    return p


def _fmt(x: float) -> str:
    return f"{x:.6f}"
