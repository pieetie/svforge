"""
Allow ``python -m svforge`` to invoke the CLI
"""

from __future__ import annotations

from svforge.cli import main

if __name__ == "__main__":
    raise SystemExit(main())
