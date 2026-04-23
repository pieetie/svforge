"""
Writer registry: in-tree @register_writer + external entry-points

External packages add a caller by exposing a subclass of
:class:`~svforge.writers.base.CallerWriter` under the
``svforge.writers`` entry-point group in their ``pyproject.toml``::

    [project.entry-points."svforge.writers"]
    mycaller = "mypkg.writers:MyCallerWriter"
"""

from __future__ import annotations

from svforge.writers import delly, manta  # noqa: F401  ensure in-tree registration
from svforge.writers._registry import (
    available_writers,
    get_writer,
    register_writer,
)
from svforge.writers.base import CallerWriter

__all__ = [
    "CallerWriter",
    "available_writers",
    "get_writer",
    "register_writer",
]
