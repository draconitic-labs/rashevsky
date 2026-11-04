#!/usr/bin/env python3
"""SBML written by this engine, read back by other tools.

Each model goes through `rashevsky convert`. libSBML must find no errors of fatal or error severity
in the core structure, COBRApy must read the same reaction counts and growth rate, and libRoadRunner
must reach the same steady state on the rewritten kinetic model.

Pass: growth within 1e-6 of the original, steady-state species within 1e-6 relative.
"""
from __future__ import annotations

import pathlib
import subprocess
import sys
import tempfile

import cobra
import libsbml
import roadrunner

from common import BIN, CORE, DATA, ICH360, REPO, write_result


def convert(source, target) -> None:
    subprocess.run([BIN, "convert", str(source), str(target)], check=True)


def main() -> int:
    checks = {}
    with tempfile.TemporaryDirectory() as scratch:
        for name, source in (("iCH360", ICH360), ("e_coli_core", CORE), ("iML1515", DATA / "iML1515.xml.gz"), ("Human-GEM", DATA / "Human-GEM.xml")):
            out = pathlib.Path(scratch) / f"{name}.xml"
            convert(source, out)
            document = libsbml.readSBML(str(out))
            severe = sum(
                1 for i in range(document.getNumErrors())
                if document.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR
            )
            original = cobra.io.read_sbml_model(str(source))
            again = cobra.io.read_sbml_model(str(out))
            a, b = original.slim_optimize(), again.slim_optimize()
            checks[name] = {
                "libsbml_errors": severe,
                "reactions": [len(original.reactions), len(again.reactions)],
                "metabolites": [len(original.metabolites), len(again.metabolites)],
                "genes": [len(original.genes), len(again.genes)],
                "objective_original": a,
                "objective_rewritten": b,
                "pass": severe == 0 and len(original.reactions) == len(again.reactions) and len(original.genes) == len(again.genes) and abs(a - b) < 1e-6,
            }
            print(name, checks[name], flush=True)
        teusink = REPO / "models" / "teusink" / "BIOMD0000000064.xml"
        out = pathlib.Path(scratch) / "teusink.xml"
        convert(teusink, out)
        states = []
        for path in (teusink, out):
            rr = roadrunner.RoadRunner(str(path))
            rr.conservedMoietyAnalysis = True
            rr.simulate(0, 200, 2)
            rr.steadyState()
            states.append({s: rr[f"[{s}]"] for s in rr.model.getFloatingSpeciesIds()})
        worst = max(abs(states[0][k] - states[1][k]) / max(abs(states[0][k]), 1e-9) for k in states[0])
        checks["teusink_kinetic"] = {"species": len(states[0]), "max_rel_difference_roadrunner": worst, "pass": worst < 1e-6}
        print(checks["teusink_kinetic"])
    passed = all(c["pass"] for c in checks.values())
    write_result("v09_sbml_export", {"checks": checks, "libsbml_version": libsbml.getLibSBMLDottedVersion(), "pass": passed})
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
