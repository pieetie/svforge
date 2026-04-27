"""
Unified VCF / VCF.gz / BCF writer

Writers produce text VCF lines. Plain ``.vcf`` and ``.vcf.gz`` outputs are
emitted as raw bytes so that the header the writer built (``##fileformat``
immediately followed by the six ``##svforge*`` tags, then the rest of the
template) reaches disk verbatim. pysam is only used for ``.bcf``, where
BCF encoding is required; in that mode pysam always declares
``##FILTER=<ID=PASS>`` itself, which is an accepted drift for BCF and is
documented in the release notes
"""

from __future__ import annotations

import tempfile
from collections.abc import Iterable
from pathlib import Path

import pysam

from svforge.writers.base import VCFRecord


def detect_mode(path: Path) -> str:
    """
    Map a path suffix to the pysam write-mode string

    Returns one of ``"w"`` (plain VCF), ``"wz"`` (bgzipped VCF) or
    ``"wb"`` (BCF). Raises :class:`ValueError` for unknown suffixes
    """
    name = path.name.lower()
    if name.endswith(".vcf.gz"):
        return "wz"
    if name.endswith(".bcf"):
        return "wb"
    if name.endswith(".vcf"):
        return "w"
    raise ValueError(f"Cannot infer output format from {path!r}: expected .vcf, .vcf.gz or .bcf")


def write_vcf(
    out_path: str | Path,
    header_lines: Iterable[str],
    record_lines: Iterable[VCFRecord],
) -> Path:
    """
    Materialise a text VCF to a file in the format implied by ``out_path``

    Plain ``.vcf`` is written as raw text and the input header bytes hit
    disk verbatim. ``.vcf.gz`` is bgzipped from the same text buffer.
    ``.bcf`` is encoded via pysam (which introduces a pysam-managed
    ``##FILTER=<ID=PASS>`` header line; accepted drift, see module docstring)
    """
    out = Path(out_path)
    out.parent.mkdir(parents=True, exist_ok=True)
    mode = detect_mode(out)

    header_text = _lines_to_text(header_lines)
    record_text = _records_to_text(record_lines)

    if mode == "w":
        out.write_text(header_text + record_text, encoding="utf-8")
        return out

    if mode == "wz":
        with tempfile.TemporaryDirectory(prefix="svforge_") as tmpdir:
            staging = Path(tmpdir) / "staging.vcf"
            staging.write_text(header_text + record_text, encoding="utf-8")
            if out.exists():
                out.unlink()
            pysam.tabix_compress(str(staging), str(out), force=True)
        return out

    with tempfile.TemporaryDirectory(prefix="svforge_") as tmpdir:
        staging = Path(tmpdir) / "staging.vcf"
        staging.write_text(header_text + record_text, encoding="utf-8")
        with (
            pysam.VariantFile(str(staging)) as vin,
            pysam.VariantFile(str(out), mode, header=vin.header) as vout,  # type: ignore[arg-type]
        ):
            for rec in vin:
                vout.write(rec)

    return out


def _lines_to_text(lines: Iterable[str]) -> str:
    parts: list[str] = []
    for line in lines:
        parts.append(line if line.endswith("\n") else line + "\n")
    return "".join(parts)


def _records_to_text(records: Iterable[VCFRecord]) -> str:
    parts: list[str] = []
    for rec in records:
        line = rec.line
        parts.append(line if line.endswith("\n") else line + "\n")
    return "".join(parts)
