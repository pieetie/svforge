"""
Writer registry -- kept in its own module to avoid circular imports
"""

from __future__ import annotations

from collections.abc import Iterable
from functools import cache
from importlib.metadata import entry_points

from svforge.writers.base import CallerWriter

_ENTRY_POINT_GROUP = "svforge.writers"
_REGISTRY: dict[str, type[CallerWriter]] = {}


def register_writer(name: str):  # type: ignore[no-untyped-def]
    """
    Decorator that registers a :class:`CallerWriter` subclass under ``name``
    """

    def decorator(cls: type[CallerWriter]) -> type[CallerWriter]:
        if not issubclass(cls, CallerWriter):
            raise TypeError(f"{cls!r} is not a CallerWriter subclass")
        _REGISTRY[name] = cls
        cls.name = name
        return cls

    return decorator


def get_writer(name: str) -> CallerWriter:
    """
    Instantiate the writer registered under ``name``
    """
    cls = _lookup(name)
    if cls is None:
        raise KeyError(f"No writer registered under {name!r}. Known: {sorted(available_writers())}")
    return cls()


def available_writers() -> Iterable[str]:
    """
    Names of every known writer (in-tree + entry-point)
    """
    names = set(_REGISTRY.keys()) | set(_load_entry_point_classes().keys())
    return sorted(names)


def _lookup(name: str) -> type[CallerWriter] | None:
    if name in _REGISTRY:
        return _REGISTRY[name]
    return _load_entry_point_classes().get(name)


@cache
def _load_entry_point_classes() -> dict[str, type[CallerWriter]]:
    found: dict[str, type[CallerWriter]] = {}
    for ep in entry_points(group=_ENTRY_POINT_GROUP):
        try:
            cls = ep.load()
        except Exception:  # pragma: no cover -- defensive
            continue
        if isinstance(cls, type) and issubclass(cls, CallerWriter):
            found[ep.name] = cls
    return found
