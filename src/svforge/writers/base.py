"""
Abstract writer: every caller formatter inherits from :class:`CallerWriter`

Writers know nothing of the sampler, the bank or the CLI. They only know
the :class:`SV` contract and produce VCF text lines (header + records)
which :mod:`svforge.io.vcf_writer` later materialises as VCF / VCF.gz / BCF

Headers are built from real-world reference templates shipped under
:mod:`svforge.data.headers` so that the output is indistinguishable from
the target caller except for the ``##svforge*`` provenance tags injected
right after ``##fileformat``. Users who need a different build or a custom
caller variant can supply their own template via ``--header-template PATH``
as long as it carries the placeholders listed in :attr:`CallerWriter.REQUIRED_PLACEHOLDERS`
"""

from __future__ import annotations

import datetime as _dt
import logging
import re
from abc import ABC, abstractmethod
from collections.abc import Iterable, Sequence
from dataclasses import dataclass
from functools import cache
from importlib import resources
from pathlib import Path
from typing import ClassVar

from svforge.core.genome import GenomeBuild, validate_genome
from svforge.core.models import SV

SYNTHETIC_REFERENCE = "svforge_synthetic_reference.fa"

_LOG = logging.getLogger(__name__)

_CONTIG_RE = re.compile(r"^##contig=<ID=([^,>]+)")


@dataclass(frozen=True, slots=True)
class VCFRecord:
    """
    A rendered VCF data line, tagged with sort key for ordering before write.

    chrom and pos are duplicated from the line for sort efficiency; line itself
    is the fully-formatted tab-separated VCF record without trailing newline.
    """

    chrom: str
    pos: int
    line: str


def extract_contig_order(header_lines: Sequence[str]) -> list[str]:
    """Return contig IDs in the order they appear in the header."""
    contigs: list[str] = []
    for line in header_lines:
        m = _CONTIG_RE.match(line)
        if m:
            contigs.append(m.group(1))
    return contigs


class CallerWriter(ABC):
    """
    Caller-specific VCF dialect writer

    Attributes:
        name: registry key, set by :func:`register_writer`
        supported_svtypes: subset of ``{"DEL", "DUP", "INV", "INS", "BND"}``
        FILEFORMAT: VCF version string emitted by the caller (``VCFv4.1`` or
            ``VCFv4.2``)
        SOURCE_STRING: value of the ``##source`` header line; svforge always
            acknowledges that it is mimicking the upstream caller
        SAMPLE_COLUMN_ORDER: order of the sample columns in paired
            (somatic) headers, e.g. ``("NORMAL", "TUMOR")`` for Manta
        TEMPLATE_BASENAME: stem used to locate the reference header
            template in :mod:`svforge.data.headers` (``f"{stem}_{genome}.vcf.template"``)
        REQUIRED_PLACEHOLDERS: placeholders a user-supplied template must
            contain to be considered valid for this writer
    """

    name: ClassVar[str] = ""
    supported_svtypes: ClassVar[frozenset[str]] = frozenset({"DEL", "DUP", "INV", "INS", "BND"})

    FILEFORMAT: ClassVar[str] = "VCFv4.2"
    SOURCE_STRING: ClassVar[str] = "svforge"
    SAMPLE_COLUMN_ORDER: ClassVar[tuple[str, ...]] = ("TUMOR", "NORMAL")
    TEMPLATE_BASENAME: ClassVar[str] = ""
    REQUIRED_PLACEHOLDERS: ClassVar[tuple[str, ...]] = (
        "{FILE_DATE}",
        "{REFERENCE_FASTA}",
        "{TUMOR_SAMPLE}",
        "{NORMAL_SAMPLE}",
    )

    def header_lines(
        self,
        sample_name: str,
        genome: GenomeBuild = "hg38",
        *,
        provenance_tags: Sequence[str] = (),
        template_override: Path | None = None,
    ) -> list[str]:
        """
        Return the full VCF header for a single-sample VCF

        Built from the caller's reference template with the ``NORMAL`` column
        dropped and the remaining sample slot filled with ``sample_name``
        """
        template = self._load_template(genome, template_override)
        return _render_template(
            template,
            tumor_sample=sample_name,
            normal_sample=None,
            sample_order=self.SAMPLE_COLUMN_ORDER,
            source_string=self.SOURCE_STRING,
            provenance_tags=provenance_tags,
        )

    def header_lines_paired(
        self,
        tumor_sample: str,
        normal_sample: str,
        genome: GenomeBuild = "hg38",
        *,
        provenance_tags: Sequence[str] = (),
        template_override: Path | None = None,
    ) -> list[str]:
        """
        Return the full VCF header for a paired tumor/normal (somatic) VCF

        Column order follows :attr:`SAMPLE_COLUMN_ORDER`. The rest of the
        header is preserved verbatim from the reference template
        """
        template = self._load_template(genome, template_override)
        return _render_template(
            template,
            tumor_sample=tumor_sample,
            normal_sample=normal_sample,
            sample_order=self.SAMPLE_COLUMN_ORDER,
            source_string=self.SOURCE_STRING,
            provenance_tags=provenance_tags,
        )

    def _load_template(
        self,
        genome: GenomeBuild,
        override: Path | None,
    ) -> tuple[str, ...]:
        if not self.TEMPLATE_BASENAME:
            raise RuntimeError(f"{type(self).__name__} did not set TEMPLATE_BASENAME")
        if override is not None:
            return _load_override_template(override, self.name, self.REQUIRED_PLACEHOLDERS)
        validate_genome(genome)
        name = f"{self.TEMPLATE_BASENAME}_{genome}.vcf.template"
        return _load_bundled_template(name)

    @abstractmethod
    def format_record(self, sv: SV, sample_name: str) -> list[VCFRecord]:
        """
        Return one or more VCF records for ``sv``

        Manta returns two records for a BND event (one per mate). DELLY and
        simple symbolic records (DEL/DUP/INV/INS) return a single-item list
        """

    def format_records(self, svs: Iterable[SV], sample_name: str) -> list[VCFRecord]:
        """
        Format a whole SV collection to a list of VCF records (not sorted)
        """
        records: list[VCFRecord] = []
        for sv in svs:
            if sv.svtype not in self.supported_svtypes:
                raise ValueError(f"Writer {self.name!r} does not support svtype {sv.svtype!r}")
            records.extend(self.format_record(sv, sample_name))
        return records

    def format_records_sorted(
        self,
        svs: Iterable[SV],
        sample_name: str,
        header_lines: Sequence[str],
    ) -> list[VCFRecord]:
        """
        Like :meth:`format_records` but ordered by ``(contig_index, pos)``.

        Contig order is taken from ``##contig=`` lines in the rendered header, not
        from a hardcoded genome list, so it stays aligned with
        user-supplied ``--header-template`` and other template overrides.
        """
        raw = self.format_records(svs, sample_name)
        return self.sort_records(raw, extract_contig_order(header_lines))

    def sort_records(
        self,
        records: Sequence[VCFRecord],
        contig_order: Sequence[str],
    ) -> list[VCFRecord]:
        """
        Sort VCF records by (contig_index, pos), using header contig order.
        """
        contig_index = {c: i for i, c in enumerate(contig_order)}
        warned: set[str] = set()

        def key(rec: VCFRecord) -> tuple[int, int]:
            idx = contig_index.get(rec.chrom)
            if idx is None:
                if rec.chrom not in warned:
                    _LOG.warning("Records on undeclared contig %s", rec.chrom)
                    warned.add(rec.chrom)
                idx = 10**9
            return (idx, rec.pos)

        return sorted(records, key=key)


@cache
def _load_bundled_template(name: str) -> tuple[str, ...]:
    try:
        resource = resources.files("svforge.data.headers").joinpath(name)
    except ModuleNotFoundError as exc:  # pragma: no cover -- packaging-time guard
        raise RuntimeError("svforge.data.headers not installed") from exc
    if not resource.is_file():
        raise ValueError(
            f"No bundled header template {name!r}. Provide your own via --header-template PATH."
        )
    text = resource.read_text(encoding="utf-8")
    return tuple(text.splitlines())


def _load_override_template(
    path: Path,
    caller: str,
    required: Sequence[str],
) -> tuple[str, ...]:
    if not path.is_file():
        raise FileNotFoundError(f"Header template {path} does not exist")
    text = path.read_text(encoding="utf-8")
    missing = [ph for ph in required if ph not in text]
    if missing:
        raise ValueError(
            f"Header template {path} is missing required placeholders for "
            f"caller {caller!r}: {', '.join(missing)}"
        )
    return tuple(text.splitlines())


def _render_template(
    template: Sequence[str],
    *,
    tumor_sample: str,
    normal_sample: str | None,
    sample_order: Sequence[str],
    source_string: str,
    provenance_tags: Sequence[str],
) -> list[str]:
    """
    Apply placeholder substitutions and inject svforge + source lines
    """
    lines = list(template)
    if not lines or not lines[0].startswith("##fileformat="):
        raise RuntimeError("template must start with ##fileformat=")

    file_date = _dt.datetime.now(_dt.timezone.utc).strftime("%Y%m%d")
    cmdline = _cmdline_from_tags(provenance_tags)
    out: list[str] = [lines[0]]
    out.extend(provenance_tags)

    saw_source = False
    for raw in lines[1:]:
        line = (
            raw.replace("{FILE_DATE}", file_date)
            .replace("{CMDLINE}", cmdline)
            .replace("{REFERENCE_FASTA}", SYNTHETIC_REFERENCE)
        )
        if line.startswith("#CHROM"):
            line = _render_chrom_line(line, tumor_sample, normal_sample, sample_order)
            if not saw_source:
                out.append(f"##source={source_string}")
                saw_source = True
            out.append(line)
            continue
        if line.startswith("##source"):
            out.append(f"##source={source_string}")
            saw_source = True
            continue
        out.append(line)
    return out


def _render_chrom_line(
    line: str,
    tumor_sample: str,
    normal_sample: str | None,
    sample_order: Sequence[str],
) -> str:
    fields = line.split("\t")
    if len(fields) < 10:
        raise RuntimeError(f"malformed #CHROM line in template: {line!r}")
    header = fields[:9]
    if normal_sample is None:
        return "\t".join([*header, tumor_sample])
    mapping = {"TUMOR": tumor_sample, "NORMAL": normal_sample}
    ordered = [mapping[role] for role in sample_order]
    return "\t".join([*header, *ordered])


def _cmdline_from_tags(provenance_tags: Sequence[str]) -> str:
    for tag in provenance_tags:
        if tag.startswith("##svforgeCommand="):
            return tag.split("=", 1)[1]
    return ""
