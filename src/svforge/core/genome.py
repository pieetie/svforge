"""
Reference genome contig sizes

svForge V1 supports hg38 only: the bundled header templates and contig
lengths are hg38-specific. Attempting to use any other build raises a
clear error pointing to ``--header-template`` for self-service override
"""

from __future__ import annotations

from typing import Literal

GenomeBuild = Literal["hg38"]

HG38_CONTIGS: dict[str, int] = {
    "chr1": 248_956_422,
    "chr2": 242_193_529,
    "chr3": 198_295_559,
    "chr4": 190_214_555,
    "chr5": 181_538_259,
    "chr6": 170_805_979,
    "chr7": 159_345_973,
    "chr8": 145_138_636,
    "chr9": 138_394_717,
    "chr10": 133_797_422,
    "chr11": 135_086_622,
    "chr12": 133_275_309,
    "chr13": 114_364_328,
    "chr14": 107_043_718,
    "chr15": 101_991_189,
    "chr16": 90_338_345,
    "chr17": 83_257_441,
    "chr18": 80_373_285,
    "chr19": 58_617_616,
    "chr20": 64_444_167,
    "chr21": 46_709_983,
    "chr22": 50_818_468,
    "chrX": 156_040_895,
    "chrY": 57_227_415,
    "chrM": 16_569,
}

SUPPORTED_GENOMES: tuple[str, ...] = ("hg38",)


def validate_genome(genome: str) -> GenomeBuild:
    """
    Return ``genome`` if supported, else raise ``ValueError`` with a
    message pointing the user at the ``--header-template`` escape hatch
    """
    if genome in SUPPORTED_GENOMES:
        return "hg38"
    raise ValueError(
        f"No bundled header template for genome {genome!r}.\n"
        f"Provide your own via --header-template PATH, or use "
        f"--genome {SUPPORTED_GENOMES[0]}."
    )


def get_contigs(genome: GenomeBuild) -> dict[str, int]:
    """
    Return the contig -> length mapping for the given build
    """
    validate_genome(genome)
    return dict(HG38_CONTIGS)


def contig_length(genome: GenomeBuild, chrom: str) -> int:
    """
    Length of ``chrom`` on ``genome``, raises KeyError if unknown
    """
    return get_contigs(genome)[chrom]


def normalize_chromosomes(
    chroms: list[str] | tuple[str, ...],
    genome: GenomeBuild,
) -> list[str]:
    """
    Normalise a list of chromosome names against ``genome``

    Accepts short aliases (``"1"``, ``"X"``, ``"M"``, ``"MT"``) and adds the
    ``chr`` prefix as needed. Unknown chromosomes raise :class:`ValueError`
    with the list of valid names
    """
    valid = get_contigs(genome)
    out: list[str] = []
    for raw in chroms:
        name = raw.strip()
        if not name:
            continue
        candidate = name if name.startswith("chr") else f"chr{name}"
        if candidate == "chrMT":
            candidate = "chrM"
        if candidate not in valid:
            raise ValueError(
                f"Unknown chromosome {raw!r} for genome {genome}. Valid: {', '.join(valid.keys())}"
            )
        if candidate not in out:
            out.append(candidate)
    return out
