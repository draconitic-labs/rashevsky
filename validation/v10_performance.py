#!/usr/bin/env python3
"""Measured speed on this machine, and FVA on a genome-scale model against COBRApy.

Times are wall clock on the machine named in the results, with whatever else it was running, so they
are indicative. COBRApy uses GLPK with 4 processes; the engine uses HiGHS with all logical cores.

Pass: FVA on iML1515 within 1e-5 of COBRApy over reactions with finite ranges, at 100 % of the optimum.
"""
from __future__ import annotations

import os
import platform
import statistics
import subprocess
import sys
import time

import cobra
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion

from common import BIN, CORE, DATA, ICH360, engine, write_result


def timed(fn, repeats=3):
    runs = []
    result = None
    for _ in range(repeats):
        started = time.perf_counter()
        result = fn()
        runs.append(time.perf_counter() - started)
    return statistics.median(runs), result


def cpu_name() -> str:
    for line in open("/proc/cpuinfo"):
        if line.startswith("model name"):
            return line.split(":", 1)[1].strip()
    return platform.processor()


def main() -> int:
    table = {}
    iml = DATA / "iML1515.xml.gz"
    human = DATA / "Human-GEM.xml"

    for name, path in (("iCH360", ICH360), ("e_coli_core", CORE), ("iML1515", iml), ("Human-GEM", human)):
        seconds, _ = timed(lambda p=path: engine("info", p))
        fba, _ = timed(lambda p=path: engine("fba", p))
        ref_load, model = timed(lambda p=path: cobra.io.read_sbml_model(str(p)), repeats=1)
        ref_fba, _ = timed(lambda m=model: m.slim_optimize())
        table[name] = {"rashevsky_read_seconds": seconds, "rashevsky_read_and_fba_seconds": fba, "cobrapy_read_seconds": ref_load, "cobrapy_fba_seconds": ref_fba}
        print(name, table[name], flush=True)

    fva_s, ours = timed(lambda: engine("fva", iml), repeats=1)
    model = cobra.io.read_sbml_model(str(iml))
    ref_s, ref = timed(lambda: flux_variability_analysis(model, processes=4), repeats=1)
    worst = 0.0
    worst_id = None
    for row in ours:
        d = max(abs(row["minimum"] - ref.loc[row["id"], "minimum"]), abs(row["maximum"] - ref.loc[row["id"], "maximum"]))
        if d > worst:
            worst, worst_id = d, row["id"]
    table["iML1515_fva"] = {"reactions": len(ours), "rashevsky_seconds": fva_s, "cobrapy_glpk_4_processes_seconds": ref_s, "max_abs_difference": worst, "worst_reaction": worst_id}
    print(table["iML1515_fva"], flush=True)

    del_s, deletions = timed(lambda: engine("deletion", iml, "--kind", "gene"), repeats=1)
    model = cobra.io.read_sbml_model(str(iml))
    ref_del, _ = timed(lambda: single_gene_deletion(model, processes=4), repeats=1)
    table["iML1515_single_gene_deletion"] = {"genes": len(deletions), "rashevsky_seconds": del_s, "cobrapy_glpk_4_processes_seconds": ref_del}
    print(table["iML1515_single_gene_deletion"], flush=True)

    tasks_s, _ = timed(lambda: engine("tasks", human, DATA / "metabolicTasks_Essential.txt"), repeats=1)
    table["human_gem_essential_tasks_seconds"] = tasks_s

    sim_s, _ = timed(lambda: engine("simulate", DATA / "BIOMD0000000064.xml", "--end", 100, "--steps", 1000, "--method", "bdf", "--symbols", "G6P"))
    table["teusink_time_course_100_min_1000_points_seconds"] = sim_s

    sample_s = timed(lambda: subprocess.run([BIN, "sample", str(CORE), "--samples", "2500", "--thinning", "100", "--fraction", "0.9"], check=True, capture_output=True), repeats=1)[0]
    table["core_chrr_10000_samples_seconds"] = sample_s

    passed = worst < 1e-5
    write_result("v10_performance", {"cpu": cpu_name(), "logical_cores": os.cpu_count(), "benchmarks": table, "pass": passed})
    print("worst FVA difference", worst, worst_id)
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
