#!/usr/bin/env python3
"""Flux sampling on e_coli_core: CHRR and OptGP against COBRApy's OptGP.

The biomass reaction is held at 90 % of its optimum or more, which makes the polytope small enough
for all samplers to mix. Means are compared in units of the pooled standard deviation of the sampled
flux, since the absolute scale differs across reactions.

Pass: split R-hat at most 1.05 for CHRR; over reactions that vary, the mean of each sampler within
0.15 of a flux standard deviation of the mean of the others.
"""
from __future__ import annotations

import json
import subprocess
import sys

import cobra
import numpy as np
from cobra.sampling import OptGPSampler

from common import BIN, CORE, write_result

RHAT_LIMIT = 1.05
MEAN_LIMIT = 0.15


def run_ours(method: str, thinning: int) -> dict:
    out = subprocess.run(
        [BIN, "sample", CORE, "--method", method, "--chains", "4", "--samples", "2500", "--thinning", str(thinning),
         "--burn-in", "20000", "--seed", "11", "--fraction", "0.9"],
        check=True, capture_output=True, text=True,
    )
    return json.loads(out.stdout)


def main() -> int:
    # OptGP mixes more slowly per move, so it gets more thinning to reach a comparable sample size.
    chrr = run_ours("chrr", 100)
    optgp = run_ours("optgp", 1000)
    model = cobra.io.read_sbml_model(str(CORE))
    optimum = model.slim_optimize()
    model.reactions.BIOMASS_Ecoli_core_w_GAM.lower_bound = 0.9 * optimum
    sampler = OptGPSampler(model, processes=4, thinning=100, seed=11)
    ref = sampler.sample(10000)
    valid = sampler.validate(ref.to_numpy())
    ref = ref[valid == "v"]

    reactions = chrr["reactions"]
    worst = {"chrr_vs_cobrapy": (0.0, None), "optgp_vs_cobrapy": (0.0, None), "chrr_vs_optgp": (0.0, None)}
    compared = 0
    for k, rid in enumerate(reactions):
        sd = float(ref[rid].std())
        if sd < 1e-6:
            continue
        compared += 1
        means = {"chrr": chrr["diagnostics"]["mean"][k], "optgp": optgp["diagnostics"]["mean"][k], "cobrapy": float(ref[rid].mean())}
        for name, (a, b) in {"chrr_vs_cobrapy": ("chrr", "cobrapy"), "optgp_vs_cobrapy": ("optgp", "cobrapy"), "chrr_vs_optgp": ("chrr", "optgp")}.items():
            d = abs(means[a] - means[b]) / sd
            if d > worst[name][0]:
                worst[name] = (d, rid)
    passed = chrr["diagnostics"]["max_split_rhat"] <= RHAT_LIMIT and all(v[0] <= MEAN_LIMIT for v in worst.values())
    write_result("v06_sampling", {
        "polytope_dimension": chrr["dimension"],
        "reactions_compared": compared,
        "cobrapy_valid_samples": int(len(ref)),
        "chrr_diagnostics": {k: v for k, v in chrr["diagnostics"].items() if k.startswith(("max", "min"))},
        "optgp_diagnostics": {k: v for k, v in optgp["diagnostics"].items() if k.startswith(("max", "min"))},
        "chrr_seconds": chrr.get("seconds"),
        "optgp_seconds": optgp.get("seconds"),
        "max_mean_difference_in_flux_sd": {k: {"value": v[0], "reaction": v[1]} for k, v in worst.items()},
        "limit": MEAN_LIMIT,
        "pass": passed,
    })
    print(chrr["diagnostics"]["max_split_rhat"], chrr["diagnostics"]["min_ess"], optgp["diagnostics"]["max_split_rhat"], optgp["diagnostics"]["min_ess"], worst, "dim", chrr["dimension"])
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
