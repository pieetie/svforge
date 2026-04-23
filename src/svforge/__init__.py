"""
svForge -- synthetic structural-variant VCF generator
"""

from __future__ import annotations

from importlib.metadata import PackageNotFoundError
from importlib.metadata import version as _pkg_version

try:
    __version__ = _pkg_version("svforge")
except PackageNotFoundError:  # pragma: no cover -- running from a source checkout
    __version__ = "0.0.0+local"

__all__ = ["__version__"]
