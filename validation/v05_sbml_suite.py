#!/usr/bin/env python3
"""SBML Test Suite 3.5.0 semantic cases, with both integrators, and the stochastic DSMTS cases.

Applicable semantic cases: all except those that use the comp package, which this engine does not
implement. Stochastic cases: the 39 StochasticTimeCourse cases of the DSMTS; the others test the
distrib package. The suite's criterion is Z in meanRange and Y in sdRange per time point, and it
states that a correct simulator can fail 0-1 times per test, 2-3 times warranting a check.

Pass: at least 95 % of applicable semantic cases with Rodas4 and with BDF, and no stochastic case
with more than 3 failing points per output.
"""
from __future__ import annotations

import collections
import csv
import math
import re
import sys
import time

from common import DATA, engine, write_result

SEMANTIC = DATA / "semantic" / "semantic"
STOCHASTIC = DATA / "stochastic" / "stochastic"
RUNS = 10000


def semantic(method: str):
    t0 = time.time()
    results = engine("sbml-suite", SEMANTIC, "--method", method)
    elapsed = time.time() - t0
    applicable = [c for c in results if c["applicable"]]
    passed = [c for c in applicable if c["passed"]]
    tags = collections.Counter()
    for c in applicable:
        if not c["passed"]:
            tags.update(c["test_tags"])
    return {
        "method": method,
        "total_cases": len(results),
        "applicable": len(applicable),
        "passed": len(passed),
        "percent": 100.0 * len(passed) / len(applicable),
        "failed": [{"case": c["case"], "detail": c["detail"]} for c in applicable if not c["passed"]],
        "fbc_cases_passed": sum(1 for c in passed if c["test_type"].startswith("FluxBalance")),
        "fast_reaction_cases_passed": sum(1 for c in passed if "FastReaction" in c["test_tags"]),
        "failed_tags": dict(tags),
        "seconds": elapsed,
    }


def stochastic():
    out = {}
    for case_dir in sorted(p for p in STOCHASTIC.iterdir() if p.is_dir()):
        case = case_dir.name
        text = (case_dir / f"{case}-model.m").read_text(errors="ignore")
        if "StochasticTimeCourse" not in text:
            continue
        settings = dict(
            line.split(":", 1) for line in (case_dir / f"{case}-settings.txt").read_text().splitlines() if ":" in line
        )
        start, duration, steps = float(settings["start"]), float(settings["duration"]), int(settings["steps"])
        rows = list(csv.DictReader(open(case_dir / f"{case}-results.csv")))
        variables = [v.strip() for v in settings["variables"].split(",")]
        mean_range = tuple(float(x) for x in re.findall(r"-?\d+\.?\d*", settings["meanRange"]))
        sd_range = tuple(float(x) for x in re.findall(r"-?\d+\.?\d*", settings["sdRange"]))
        sbml = sorted(case_dir.glob("*-sbml-l3v2.xml")) or sorted(case_dir.glob("*-sbml-l*.xml"))
        amounts = ",".join(f"{v}#a" for v in variables)
        res = engine("ssa", sbml[-1], "--end", duration, "--steps", steps, "--runs", RUNS, "--seed", 20261104, "--symbols", amounts)
        z_fail = y_fail = 0
        for v in variables:
            mean = res["mean"][f"{v}#a"]
            sd = res["sd"][f"{v}#a"]
            for k, row in enumerate(rows):
                mu, sigma = float(row[f"{v}-mean"]), float(row[f"{v}-sd"])
                if sigma > 0:
                    z = math.sqrt(RUNS) * (mean[k] - mu) / sigma
                    y = math.sqrt(RUNS / 2) * (sd[k] ** 2 / sigma**2 - 1)
                    if not mean_range[0] <= z <= mean_range[1]:
                        z_fail += 1
                    if not sd_range[0] <= y <= sd_range[1]:
                        y_fail += 1
                else:
                    if abs(mean[k] - mu) > 1e-9:
                        z_fail += 1
        out[case] = {"z_failures": z_fail, "y_failures": y_fail, "points": len(rows) * len(variables)}
        print(f"  stochastic {case}: Z fails {z_fail}, Y fails {y_fail} of {len(rows) * len(variables)}", flush=True)
    return out


def main() -> int:
    rodas = semantic("rodas4")
    print(f"rodas4 {rodas['passed']}/{rodas['applicable']} = {rodas['percent']:.2f} % in {rodas['seconds']:.0f} s", flush=True)
    bdf = semantic("bdf")
    print(f"bdf {bdf['passed']}/{bdf['applicable']} = {bdf['percent']:.2f} % in {bdf['seconds']:.0f} s", flush=True)
    stoch = stochastic()
    worst_points = max((c["z_failures"] + c["y_failures"]) for c in stoch.values())
    passed = rodas["percent"] >= 95 and bdf["percent"] >= 95 and worst_points <= 3
    write_result("v05_sbml_suite", {
        "suite": "SBML Test Suite 3.5.0",
        "rodas4": rodas,
        "bdf": bdf,
        "stochastic_runs_per_case": RUNS,
        "stochastic": stoch,
        "stochastic_cases": len(stoch),
        "stochastic_failed_points_total": sum(c["z_failures"] + c["y_failures"] for c in stoch.values()),
        "pass": passed,
    })
    print(f"stochastic {len(stoch)} cases, worst failing points in a case {worst_points}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
