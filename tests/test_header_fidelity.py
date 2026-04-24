"""
Header-fidelity test: svforge output must match the real-caller template

The test generates a paired tumor/normal header for each writer with zero
records and compares it line-by-line to the reference template. Only five
categories of lines are ignored (they are either dynamic or deliberately
different):

- ``##svforge*``  (not present in the reference)
- ``##source``    (delibarately different from real output)
- ``##fileDate``  (dynamic)
- ``##cmdline``   (dynamic)
- ``##reference`` (dynamic: substituted with a synthetic placeholder)

Any other divergence means the writer drifted from the real caller and
must be fixed before release
"""

from __future__ import annotations

from importlib import resources
from pathlib import Path

import pytest

from svforge.core.provenance import build_svforge_tags
from svforge.io.vcf_writer import write_vcf
from svforge.writers import get_writer
from svforge.writers.base import CallerWriter

def _load_template(name: str) -> list[str]:
    return (
        resources.files("svforge.data.headers")
        .joinpath(name)
        .read_text(encoding="utf-8")
        .splitlines()
    )

def _filter_ignored(lines: list[str]) -> list[str]:
    return [
        line
        for line in lines
        if not line.startswith("##svforge")
        and not line.startswith("##source")
        and not line.startswith("##fileDate")
        and not line.startswith("##cmdline")
        and not line.startswith("##reference")
    ]

def _substitute_samples(lines: list[str], tumor: str, normal: str) -> list[str]:
    out: list[str] = []
    for line in lines:
        if line.startswith("#CHROM"):
            out.append(line.replace("{TUMOR_SAMPLE}", tumor).replace("{NORMAL_SAMPLE}", normal))
        else:
            out.append(line)
    return out

@pytest.mark.parametrize(
    ("caller", "template_name"),
    [
        ("manta", "manta_somatic_hg38.vcf.template"),
        ("delly", "delly_somatic_hg38.vcf.template"),
    ],
)
def test_paired_header_matches_template(caller: str, template_name: str) -> None:
    writer: CallerWriter = get_writer(caller)
    tumor, normal = "TUMOR01", "NORMAL01"
    tags = build_svforge_tags(caller=caller, seed=42, argv=["svforge", "gen-pair"])

    generated = writer.header_lines_paired(tumor, normal, provenance_tags=tags)
    template = _substitute_samples(_load_template(template_name), tumor, normal)

    assert _filter_ignored(generated) == _filter_ignored(template)

@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_svforge_tags_are_injected_right_after_fileformat(caller: str) -> None:
    writer: CallerWriter = get_writer(caller)
    tags = build_svforge_tags(caller=caller, seed=0, argv=["svforge"])
    lines = writer.header_lines_paired("T", "N", provenance_tags=tags)

    assert lines[0].startswith("##fileformat=")
    for idx, tag in enumerate(tags, start=1):
        assert lines[idx] == tag

def test_templates_contain_no_sensitive_data() -> None:
    """
    Guard the committed templates against leaks of internal run metadata
    """
    for name in (
        "manta_somatic_hg38.vcf.template",
        "delly_somatic_hg38.vcf.template",
    ):
        text = resources.files("svforge.data.headers").joinpath(name).read_text(encoding="utf-8")
        for bad in ("P69", "pmichonn", "beegfs", "/mnt/"):
            assert bad not in text, f"{bad!r} leaked into {name}"

@pytest.mark.parametrize("caller", ["manta", "delly"])
def test_single_sample_header_keeps_one_sample_column(caller: str) -> None:
    writer: CallerWriter = get_writer(caller)
    lines = writer.header_lines("ONLY")
    chrom = next(line for line in lines if line.startswith("#CHROM"))
    assert chrom.split("\t")[-1] == "ONLY"
    assert chrom.count("\t") == 9  # 9 fixed + 1 sample = 10 fields

def test_hg19_has_no_bundled_template() -> None:
    writer = get_writer("manta")
    with pytest.raises(ValueError, match="No bundled header template"):
        writer.header_lines("S", genome="hg19")

def test_no_pass_filter_in_manta_output(tmp_path: Path) -> None:
    """
    Manta template does not declare ``##FILTER=<ID=PASS>`` so the
    generated plain VCF must not either (the real Manta somatic header
    doesn't carry it; pysam was previously injecting it, see B1)
    """
    writer = get_writer("manta")
    tags = build_svforge_tags(caller="manta", seed=0, argv=["svforge", "gen"])
    header = writer.header_lines("S", provenance_tags=tags)

    out = tmp_path / "manta_header_only.vcf"
    write_vcf(out, header, [])

    text = out.read_text(encoding="utf-8")
    assert "##FILTER=<ID=PASS" not in text, (
        "Manta output must not carry ##FILTER=<ID=PASS>; Manta's somatic "
        "template doesn't declare it"
    )
