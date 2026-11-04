#!/usr/bin/env python3
"""Goes through every workspace of the real window with real input, on a private display.

Reports errors raised in the interface and the map's frame rate, and writes screenshots. A nonzero exit
means the interface raised an error. Set CARGO_TARGET_DIR to where the release binary is, and
RASHEVSKY_DATA to a folder holding Human-GEM.xml to include the Atlas with a human model.

Usage: qa.py [screenshot folder]
"""
from __future__ import annotations

import os
import pathlib
import sys
import time

sys.path.insert(0, str(pathlib.Path(__file__).resolve().parent))
import drive

OUT = pathlib.Path(sys.argv[1]) if len(sys.argv) > 1 else pathlib.Path(__file__).resolve().parents[2] / "docs" / "screenshots"
RUN_BUTTON = (1486, 52)
# The navigation rail, from the top. Keys would work too, unless a slider has the focus.
RAIL = {"map": 100, "screen": 165, "ranges": 228, "samples": 292, "kinetics": 355, "atlas": 419, "models": 482, "checks": 546}


def go(session, view, settle=0.8):
    session.click(38, RAIL[view], settle=settle)


def problems(session):
    return [line for line in session.read_log().splitlines() if line.startswith("QA ") and ("error" in line or "rejection" in line)]


def main() -> int:
    OUT.mkdir(parents=True, exist_ok=True)
    found = []
    with drive.Session() as s:
        time.sleep(3)
        s.shot(OUT / "map.png")
        fps = [line for line in s.read_log().splitlines() if "fps" in line]
        print("map frame rate:", fps[-1] if fps else "no frames reported")
        s.click(598, 198, settle=1.0)  # PGI, on the spine between G6P and F6P
        s.key("k", settle=1.4)
        s.shot(OUT / "knockout.png")
        s.key("Escape", settle=0.5)
        s.click(967, 923, settle=1.0)  # Restore all, so that what follows starts from the whole network
        s.click(107, 863, settle=1.0)  # Oxygen uptake to zero: succinate couples to growth only without it
        go(s, "screen", settle=0.6)
        s.click(1151, 53, settle=0.5)  # At once: every pair
        s.key("Down", settle=0.2)
        s.key("Return", settle=0.5)
        s.click(*RUN_BUTTON, settle=20)
        s.shot(OUT / "screen.png")
        go(s, "map", settle=0.6)
        s.click(290, 863, settle=1.0)  # Oxygen back to unlimited
        go(s, "ranges", settle=2.5)
        s.shot(OUT / "ranges.png")
        go(s, "samples", settle=0.6)
        s.click(*RUN_BUTTON, settle=6)
        s.shot(OUT / "samples.png")
        go(s, "kinetics", settle=2.5)
        s.shot(OUT / "kinetics.png")
        s.click(636, 52, settle=3)  # the steady state and control tab
        s.shot(OUT / "control.png")
        go(s, "models", settle=0.8)
        s.shot(OUT / "models.png")
        go(s, "checks", settle=1.2)
        s.shot(OUT / "checks.png")
        go(s, "map", settle=1)
        found += problems(s)
    human = pathlib.Path(os.environ.get("RASHEVSKY_DATA", "")) / "Human-GEM.xml"
    if human.exists():
        with drive.Session(args=[str(human)]) as s:
            time.sleep(10)
            s.shot(OUT / "human-map.png")
            go(s, "atlas", settle=7)
            s.shot(OUT / "atlas.png")
            found += problems(s)
    print("errors in the interface:", found or "none")
    return 1 if found else 0


if __name__ == "__main__":
    sys.exit(main())
