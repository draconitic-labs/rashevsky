#!/usr/bin/env python3
"""Recon3D (Brunk et al. 2018) against COBRApy: loading, growth, flux variability, deletions, export.

Pass: the same reactions, metabolites and genes; biomass objective within 1e-6 relative; FVA at 100 % of the
optimum within 1e-5 of COBRApy on 300 reactions chosen by a fixed seed; the objective after each of the 2,248 gene
deletions within 1e-6 relative; the rewritten file gives the same objective in COBRApy.
"""
from __future__ import annotations

import pathlib
import random
import subprocess
import sys
import tempfile

import cobra
import libsbml
from cobra.flux_analysis import flux_variability_analysis, single_gene_deletion

from common import BIN, DATA, engine, write_result

MODEL = DATA / "Recon3D.xml.gz"


def main() -> int:
    cobra.Configuration().solver = "glpk"
    model = cobra.io.read_sbml_model(str(MODEL))
    info = engine("info", MODEL)
    shape = {
        "reactions": [info["reactions"], len(model.reactions)],
        "metabolites": [info["metabolites"], len(model.metabolites)],
        "genes": [info["genes"], len(model.genes)],
    }
    print(shape, flush=True)

    reference = model.slim_optimize()
    ours = engine("fba", MODEL)
    growth_error = abs(ours["objective"] - reference) / abs(reference)
    print("objective", ours["objective"], reference, flush=True)

    ids = sorted(r.id for r in model.reactions)
    chosen = random.Random(20261104).sample(ids, 300)
    ref_fva = flux_variability_analysis(model, reaction_list=chosen, fraction_of_optimum=1.0, processes=2)
    rows = engine("fva", MODEL, "--fraction", "1.0", "--reactions", ",".join(chosen))
    worst = 0.0
    worst_id = None
    for row in rows:
        d = max(abs(row["minimum"] - ref_fva.loc[row["id"], "minimum"]), abs(row["maximum"] - ref_fva.loc[row["id"], "maximum"]))
        if d > worst:
            worst, worst_id = d, row["id"]
    print("fva", worst, worst_id, flush=True)

    ref_deletions = single_gene_deletion(model, processes=2)
    deletion_worst = 0.0
    deletion_id = None
    changed = 0
    lethal = 0
    our_deletions = engine("deletion", MODEL, "--kind", "gene", "--order", "1")
    for row in our_deletions:
        theirs = ref_deletions.knockout[row["ids"][0]]["growth"].iloc[0]
        theirs = 0.0 if theirs != theirs else theirs
        d = abs(row["growth"] - theirs) / abs(reference)
        if d > deletion_worst:
            deletion_worst, deletion_id = d, row["ids"][0]
        changed += abs(theirs - reference) > 1e-6 * abs(reference)
        lethal += abs(theirs) < 1e-6 * abs(reference)
    print("deletions", len(our_deletions), deletion_worst, deletion_id, "changed", changed, "lethal", lethal, flush=True)

    with tempfile.TemporaryDirectory() as scratch:
        out = pathlib.Path(scratch) / "Recon3D.xml"
        subprocess.run([BIN, "convert", str(MODEL), str(out)], check=True)
        document = libsbml.readSBML(str(out))
        severe = sum(1 for i in range(document.getNumErrors()) if document.getError(i).getSeverity() >= libsbml.LIBSBML_SEV_ERROR)
        again = cobra.io.read_sbml_model(str(out))
        rewritten = again.slim_optimize()
    print("export", severe, rewritten, flush=True)

    checks = {
        "shape": {**shape, "pass": all(a == b for a, b in shape.values())},
        "objective": {"cobrapy": reference, "rashevsky": ours["objective"], "relative_difference": growth_error, "tolerance": 1e-6, "pass": growth_error < 1e-6},
        "fva": {"reactions": len(chosen), "max_abs_difference": worst, "worst_reaction": worst_id, "tolerance": 1e-5, "pass": worst < 1e-5},
        "gene_deletions": {"genes": len(our_deletions), "growth_changed": changed, "lethal": lethal, "max_relative_difference": deletion_worst, "worst_gene": deletion_id, "tolerance": 1e-6, "pass": len(our_deletions) == len(model.genes) and deletion_worst < 1e-6},
        "export": {"libsbml_errors": severe, "objective_rewritten": rewritten, "reaction_count": len(again.reactions), "pass": severe == 0 and len(again.reactions) == len(model.reactions) and abs(rewritten - reference) / abs(reference) < 1e-6},
    }
    passed = all(c["pass"] for c in checks.values())
    write_result("v11_recon3d", {"pass": passed, "cobrapy_version": cobra.__version__, "checks": checks})
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
