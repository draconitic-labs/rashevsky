#!/usr/bin/env python3
"""Teusink et al. 2000 yeast glycolysis, BIOMD0000000064, and metabolic control analysis.

Reference values: Table 4 of Teusink et al., Eur. J. Biochem. 267 (2000) 5313, as reproduced by the
BioModels curators in the model's notes, with their curated model column. Independent simulators on
the same SBML file: libRoadRunner for the steady state, COPASI for the scaled control coefficients.
(libRoadRunner's control matrices raise an exception on this model.)

Pass: fluxes and concentrations within one unit in the last digit printed in the curators' table;
steady state within 1e-5 relative of libRoadRunner; control coefficients within 1e-4 absolute of
COPASI; control theorems hold to 1e-8 on every model listed.
"""
from __future__ import annotations

import sys

import numpy as np
import roadrunner

from common import DATA, engine, write_result

MODEL = DATA / "BIOMD0000000064.xml"
# (reaction, multiplier, published, curated model value, unit of last printed digit)
FLUXES = [
    ("vGLK", 1, 88, 88, 1),
    ("vADH", 1, 129, 129, 1),
    ("vGLYCO", 1, 6, 6, 1),
    ("vTreha", 2, 4.8, 4.8, 0.1),
    ("vG3PDH", 1, 18.2, 18.2, 0.1),
    ("vSUC", 1, 3.6, 3.6, 0.1),
]
# (species, published, curated model value)
CONCENTRATIONS = [
    ("G6P", 1.07, 1.03),
    ("F6P", 0.11, 0.11),
    ("F16P", 0.6, 0.6),
    ("P3G", 0.36, 0.36),
    ("P2G", 0.04, 0.04),
    ("PEP", 0.07, 0.07),
    ("PYR", 8.52, 8.52),
    ("ACE", 0.17, 0.17),
    ("ATP", 2.51, 2.51),
    ("ADP", 1.29, 1.29),
    ("AMP", 0.3, 0.3),
    ("NAD", 1.55, 1.55),
    ("NADH", 0.04, 0.04),
]
THEOREM_MODELS = {
    "BIOMD0000000064": "Teusink 2000 yeast glycolysis",
    "BIOMD0000000176": "Conant and Wolfe 2007 yeast glycolysis",
    "BIOMD0000000172": "Pritchard and Kell 2002 yeast glycolysis",
    "BIOMD0000000071": "Bakker 2001 Trypanosoma glycolysis",
    "BIOMD0000000061": "Hynne 2001 yeast glycolysis",
}
THEOREM_TOL = 1e-8


def copasi_control(mca: dict) -> tuple[float, float, str]:
    """Largest absolute differences from COPASI's scaled flux and concentration control matrices."""
    import basico
    import COPASI

    dm = basico.load_model(str(MODEL))
    basico.run_steadystate()
    task = dm.getTask("Metabolic Control Analysis")
    task.getProblem().setSteadyStateRequested(True)
    task.process(True)
    method = task.getMethod()
    reactions = basico.get_reactions()
    rid = {name: row.sbml_id for name, row in reactions.iterrows()}
    species = basico.get_species()
    sid = {name: row.sbml_id for name, row in species.iterrows()}

    def table(annotated):
        rows = [str(x) for x in annotated.getAnnotationsString(0)]
        cols = [str(x) for x in annotated.getAnnotationsString(1)]
        array = annotated.getArray() if hasattr(annotated, "getArray") else annotated
        matrix = {}
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                index = COPASI.SizeTStdVector()
                index.append(i)
                index.append(j)
                matrix[(r, c)] = array.get(index) if hasattr(array, "get") else annotated.get(i, j)
        return matrix

    cj = table(method.getScaledFluxCCAnn())
    cs = table(method.getScaledConcentrationCCAnn())
    clean = lambda name: name.strip("()").split("{")[0]
    ours_reac = mca["reactions"]
    ours_states = mca["states"]
    cj_ours = np.array(mca["flux_control"])
    cs_ours = np.array(mca["concentration_control"])
    worst_j = 0.0
    for (r, c), value in cj.items():
        # COPASI's tables carry an extra "Summation Error" row and column.
        if clean(r) not in rid or clean(c) not in rid:
            continue
        i = ours_reac.index(rid[clean(r)])
        j = ours_reac.index(rid[clean(c)])
        worst_j = max(worst_j, abs(value - cj_ours[i, j]))
    worst_s = 0.0
    for (r, c), value in cs.items():
        name = sid.get(clean(r))
        state = f"{name}#a"
        if state not in ours_states or clean(c) not in rid:
            continue
        i = ours_states.index(state)
        j = ours_reac.index(rid[clean(c)])
        worst_s = max(worst_s, abs(value - cs_ours[i, j]))
    return worst_j, worst_s, COPASI.CVersion.VERSION.getVersion()


def main() -> int:
    mca = engine("mca", MODEL, "--integrate", 200)
    state = mca["steady_state"]
    table = []
    ok = True
    for rid, mult, published, curated, unit in FLUXES:
        value = state["fluxes"][rid] * mult
        good = abs(value - curated) <= unit
        ok &= good
        table.append({"quantity": rid, "rashevsky": value, "published": published, "curated_model": curated, "tolerance": unit, "pass": good})
    for sid, published, curated in CONCENTRATIONS:
        value = state["species"][sid]
        good = abs(value - curated) <= 0.01
        ok &= good
        table.append({"quantity": sid, "rashevsky": value, "published": published, "curated_model": curated, "tolerance": 0.01, "pass": good})

    # Independent simulator.
    rr = roadrunner.RoadRunner(str(MODEL))
    rr.conservedMoietyAnalysis = True
    rr.simulate(0, 200, 2)
    rr.steadyState()
    worst = 0.0
    worst_id = None
    for sid, ours in state["species"].items():
        try:
            ref = rr[f"[{sid}]"]
        except RuntimeError:
            continue
        d = abs(ours - ref) / max(abs(ref), 1e-9)
        if d > worst:
            worst, worst_id = d, sid
    flux_worst = 0.0
    for rid, ours in state["fluxes"].items():
        ref = rr[rid]
        flux_worst = max(flux_worst, abs(ours - ref) / max(abs(ref), 1e-9))
    # Control coefficients from COPASI's metabolic control analysis task on the same file.
    cj_worst, cs_worst, copasi_version = copasi_control(mca)
    agree = worst < 1e-5 and flux_worst < 1e-5 and cj_worst < 1e-4 and cs_worst < 1e-4

    # Control theorems on all models.
    theorems = {}
    for m, title in THEOREM_MODELS.items():
        try:
            res = engine("mca", DATA / f"{m}.xml", "--integrate", 500)
        except Exception as exc:  # noqa: BLE001
            theorems[m] = {"title": title, "error": str(exc), "pass": False}
            continue
        # The theorems are statements about scaled coefficients; unscaled residuals carry units.
        scaled = {k: v for k, v in res["theorems"].items() if not k.startswith("unscaled")}
        worst_t = max(scaled.values())
        theorems[m] = {
            "title": title,
            "reactions": len(res["reactions"]),
            "states": len(res["states"]),
            "conserved_moieties": res["conserved_moieties"],
            "residuals": res["theorems"],
            "worst": worst_t,
            "pass": worst_t < THEOREM_TOL,
        }
    theorem_ok = all(t["pass"] for t in theorems.values())
    passed = ok and agree and theorem_ok
    write_result("v04_teusink_mca", {
        "table": table,
        "roadrunner_version": roadrunner.__version__,
        "steady_state_max_rel_difference_species": worst,
        "worst_species": worst_id,
        "steady_state_max_rel_difference_fluxes": flux_worst,
        "copasi_version": copasi_version,
        "flux_control_max_abs_difference_vs_copasi": cj_worst,
        "concentration_control_max_abs_difference_vs_copasi": cs_worst,
        "theorems": theorems,
        "pass": passed,
    })
    print(f"table ok {ok}; roadrunner species {worst:.2e} ({worst_id}) fluxes {flux_worst:.2e} COPASI Cj {cj_worst:.2e} Cs {cs_worst:.2e}")
    for row in table:
        print(f"  {'ok ' if row['pass'] else 'BAD'} {row['quantity']:7s} {row['rashevsky']:.4f} published {row['published']} curated {row['curated_model']}")
    for m, t in theorems.items():
        print(f"  {m} {t.get('title')}: {t.get('worst', t.get('error'))}")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
