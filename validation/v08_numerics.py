#!/usr/bin/env python3
"""Numerical studies: integrator convergence, dynamic FBA convergence and TFA consistency.

Integrators: on the Teusink glycolysis model, the error at t = 2 against a tight reference falls with
the tolerance for Rodas4 and BDF, and the two agree. The fixed-step order of Rodas4 is established in
the Rust test `rodas_converges_with_order_four`, and Robertson's problem is checked there too.

Dynamic FBA: e_coli_core on glucose through depletion. The scheme is exact while uptake stays fixed and
first-order once the substrate runs out, so the final biomass converges with order 1 in the step.

TFA: with eQuilibrator energies, growth stays at the FBA optimum, no flux runs against its DrG', and
each DrG' agrees with DrG'0 + RT sum n ln c within two standard uncertainties.
"""
from __future__ import annotations

import json
import math
import sys
import tempfile

from common import CORE, DATA, REPO, engine, write_result

TEUSINK = DATA / "BIOMD0000000064.xml"
THERMO = REPO / "data" / "thermo" / "e_coli_core.json"


def integrators() -> dict:
    tight = engine("simulate", TEUSINK, "--end", 2, "--steps", 4, "--method", "rodas4", "--rtol", 1e-12, "--atol", 1e-14, "--symbols", "G6P,F16P,ATP,NADH,PYR")
    rows = {}
    for method in ("rodas4", "bdf"):
        for rtol in (1e-4, 1e-6, 1e-8, 1e-10):
            res = engine("simulate", TEUSINK, "--end", 2, "--steps", 4, "--method", method, "--rtol", rtol, "--atol", rtol * 1e-3, "--symbols", "G6P,F16P,ATP,NADH,PYR")
            err = max(abs(res["series"][k][i] - tight["series"][k][i]) / max(abs(tight["series"][k][i]), 1e-12) for k in tight["series"] for i in range(5))
            rows[f"{method}_{rtol:g}"] = {"max_rel_error": err, "steps": res["stats"]["steps"], "rejected": res["stats"]["rejected"], "rhs_evals": res["stats"]["rhs_evals"]}
    ok = all(rows[f"{m}_1e-10"]["max_rel_error"] < rows[f"{m}_1e-06"]["max_rel_error"] < rows[f"{m}_0.0001"]["max_rel_error"] for m in ("rodas4", "bdf"))
    ok &= rows["rodas4_1e-10"]["max_rel_error"] < 1e-7 and rows["bdf_1e-10"]["max_rel_error"] < 1e-6
    return {"study": rows, "pass": ok}


def dfba() -> dict:
    settings = {
        "biomass": 0.01,
        "t_end": 10.0,
        "dt": 0.0,
        "substrates": [{"exchange": "EX_glc__D_e", "initial": 20.0, "vmax": 10.0, "km": 0.5}],
        "tracked": [["EX_ac_e", 0.0]],
    }
    out = {}
    finals = []
    for dt in (0.04, 0.02, 0.01, 0.005):
        settings["dt"] = dt
        with tempfile.NamedTemporaryFile("w", suffix=".json", delete=False) as handle:
            json.dump(settings, handle)
        res = engine("dfba", CORE, handle.name)
        biomass = res["biomass"]
        glucose = [c[0] for c in res["concentrations"]]
        finals.append(biomass[-1])
        out[f"dt_{dt:g}"] = {"biomass_end": biomass[-1], "glucose_end": glucose[-1]}
    # Successive differences halve for a first-order scheme.
    d1, d2, d3 = abs(finals[1] - finals[0]), abs(finals[2] - finals[1]), abs(finals[3] - finals[2])
    orders = [math.log2(d1 / d2), math.log2(d2 / d3)]
    out["observed_order"] = orders
    out["pass"] = all(0.7 < o < 1.3 for o in orders)
    return out


def tfa() -> dict:
    fba = engine("fba", CORE)["objective"]
    res = engine("tfa", CORE, THERMO)
    data = json.load(open(THERMO))
    gibbs = res["gibbs"]
    flux = res["fluxes"]["fluxes"]
    against = {k: (gibbs[k], flux[k]) for k in gibbs if flux[k] * gibbs[k] > 1e-7}
    # Recompute each DrG' from the concentrations reported, with an independent stoichiometry parser.
    from v02_helpers import reaction_table

    rt = 8.314462618e-3 * data["temperature"]
    table = reaction_table(CORE)
    worst = 0.0
    for rid, energy in data["reactions"].items():
        term = sum(c * math.log(res["concentrations"][m]) for m, c in table[rid].items() if m not in data["excluded_metabolites"])
        expected = energy["dg0"] + rt * term
        worst = max(worst, abs(gibbs[rid] - expected) - 2.0 * energy["sigma"])
    return {
        "fba_objective": fba,
        "tfa_objective": res["fluxes"]["objective"],
        "reactions_with_energies": len(data["reactions"]),
        "source": data["source"],
        "fluxes_against_gibbs": against,
        "largest_excess_beyond_two_sigma": worst,
        "pass": not against and worst < 1e-6 and abs(fba - res["fluxes"]["objective"]) < 1e-6,
    }


def main() -> int:
    results = {"integrators": integrators(), "dfba": dfba(), "tfa": tfa()}
    results["pass"] = all(v["pass"] for v in results.values())
    write_result("v08_numerics", results)
    print(json.dumps({k: v["pass"] for k, v in results.items() if isinstance(v, dict)}))
    print(json.dumps(results["integrators"]["study"], indent=0)[:900])
    print(results["dfba"]["observed_order"], {k: v for k, v in results["tfa"].items() if k != "source"})
    return 0 if results["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
