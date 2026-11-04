#!/usr/bin/env python3
"""Fetches the reference data and runs every validation in order, then rewrites docs/VALIDATION.md."""
from __future__ import annotations

import pathlib
import subprocess
import sys

HERE = pathlib.Path(__file__).resolve().parent
SCRIPTS = sorted(p for p in HERE.glob("v[0-9][0-9]_*.py"))


def main() -> int:
    subprocess.run([sys.executable, str(HERE.parent / "scripts" / "fetch_data.py")], check=True)
    failed = []
    for script in SCRIPTS:
        print(f"== {script.name}", flush=True)
        if subprocess.run([sys.executable, str(script)], cwd=HERE).returncode != 0:
            failed.append(script.name)
    subprocess.run([sys.executable, str(HERE / "report.py")], cwd=HERE, check=True)
    print("failed:", failed or "none")
    return 1 if failed else 0


if __name__ == "__main__":
    sys.exit(main())
