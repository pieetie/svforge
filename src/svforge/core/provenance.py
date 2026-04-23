"""
Synthetic-data provenance helpers

Every VCF produced by svforge carries six ``##svforge*`` tags immediately
after ``##fileformat`` so that a downstream consumer can always tell the
data is synthetic and trace it back to the exact command line and seed
that generated it

The warning tag :data:`WARNING_VALUE` is deliberately not configurable: no
CLI flag, no environment variable, no hidden parameter can remove it.
That invariant is enforced by :func:`build_svforge_tags`
"""

from __future__ import annotations

from collections.abc import Sequence
from pathlib import PurePosixPath, PureWindowsPath

from svforge import __version__

WARNING_VALUE = "SYNTHETIC_DATA_DO_NOT_USE_FOR_CLINICAL_DIAGNOSIS"
DOCUMENTATION_URL = "https://github.com/pieetie/svforge"

_SVFORGE_TAG_ORDER: tuple[str, ...] = (
    "##svforgeVersion",
    "##svforgeCommand",
    "##svforgeSeed",
    "##svforgeCaller",
    "##svforgeWarning",
    "##svforgeDocumentation",
)


def build_svforge_tags(
    *,
    caller: str,
    seed: int,
    argv: Sequence[str],
    version: str = __version__,
) -> list[str]:
    """
    Return the six provenance header lines, in spec order

    ``seed`` is the effective seed, even when the user did not pass one
    explicitly. ``argv`` is typically :data:`sys.argv`; it is piped through
    :func:`sanitize_command` before being embedded in the header
    """
    command = sanitize_command(argv)
    tags = [
        f"##svforgeVersion={version}",
        f"##svforgeCommand={command}",
        f"##svforgeSeed={seed}",
        f"##svforgeCaller={caller}",
        f"##svforgeWarning={WARNING_VALUE}",
        f"##svforgeDocumentation={DOCUMENTATION_URL}",
    ]
    _assert_ordered(tags)
    return tags


def sanitize_command(argv: Sequence[str]) -> str:
    """
    Return a shell-safe command line with absolute paths stripped

    Each argv token is inspected individually: if the whole token (or the
    value portion of a ``--flag=VALUE`` form) is an absolute path -- POSIX
    ``/foo/bar``, macOS ``/Users/alice/...`` (which may contain whitespace
    on the dev workstation), or Windows-style ``C:\\path\\bar`` -- it is
    replaced by its basename so that cluster-internal or user-home
    directories never leak into a VCF header
    """
    parts: list[str] = []
    for token in argv:
        if "=" in token and token.startswith("--"):
            flag, _, value = token.partition("=")
            if _is_absolute_path(value):
                parts.append(f"{flag}={_basename(value)}")
                continue
        if _is_absolute_path(token):
            parts.append(_basename(token))
            continue
        parts.append(token)
    return " ".join(parts)


def _is_absolute_path(s: str) -> bool:
    if not s:
        return False
    if s.startswith("/"):
        return True
    return len(s) >= 3 and s[0].isalpha() and s[1:3] in (":\\", ":/")


def _basename(path: str) -> str:
    if len(path) >= 3 and path[0].isalpha() and path[1:3] in (":\\", ":/"):
        name = PureWindowsPath(path).name
    else:
        name = PurePosixPath(path).name
    return name or path


def _assert_ordered(tags: Sequence[str]) -> None:
    keys = [t.split("=", 1)[0] for t in tags]
    if tuple(keys) != _SVFORGE_TAG_ORDER:
        raise AssertionError(f"svforge tag order drift: {keys} != {_SVFORGE_TAG_ORDER}")
