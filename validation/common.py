"""Shared helpers for the validation scripts."""
from __future__ import annotations

import json
import os
import pathlib
import subprocess
import tempfile

REPO = pathlib.Path(__file__).resolve().parents[1]
RESULTS = REPO / "validation" / "results"
DATA = pathlib.Path(os.environ.get("RASHEVSKY_DATA", pathlib.Path.home() / ".cache/rashevsky-work/data"))
CORE = REPO / "models" / "e_coli_core" / "e_coli_core.xml.gz"
BIN = os.environ.get(
    "RASHEVSKY_BIN",
    str(pathlib.Path(os.environ.get("CARGO_TARGET_DIR", REPO / "target")) / "release" / "rashevsky"),
)


def engine(*args: str):
    """Runs the engine's CLI and returns its parsed JSON output."""
    with tempfile.NamedTemporaryFile(suffix=".json", dir=os.environ.get("RASHEVSKY_TMP")) as handle:
        subprocess.run([BIN, "--out", handle.name, *map(str, args)], check=True)
        return json.load(open(handle.name))


def write_result(name: str, payload: dict) -> None:
    """Records a validation result; `payload['pass']` is the verdict."""
    RESULTS.mkdir(parents=True, exist_ok=True)
    (RESULTS / f"{name}.json").write_text(json.dumps(payload, indent=2, sort_keys=True, default=lambda o: o.item() if hasattr(o, "item") else str(o)) + "\n")
    print(f"{'PASS' if payload['pass'] else 'FAIL'} {name}")


def monk_table(name: str = "Table 1"):
    """Rows of Supplementary Data 11 of Monk et al. 2017, cached as JSON next to the workbook."""
    import openpyxl
    book = pathlib.Path(os.environ.get("RASHEVSKY_MONK", DATA.parent / "papers/monk/M40.xlsx"))
    cache = book.with_suffix("." + name.replace(" ", "_") + ".json")
    if cache.exists():
        return json.load(open(cache))
    wb = openpyxl.load_workbook(book, read_only=True, data_only=True)
    rows = list(wb[name].iter_rows(values_only=True))
    head = [str(c) for c in rows[0]]
    out = [dict(zip(head, r)) for r in rows[1:]]
    json.dump(out, open(cache, "w"))
    return out
