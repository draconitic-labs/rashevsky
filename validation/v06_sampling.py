#!/usr/bin/env python3
"""Flux sampling on e_coli_core and iCH360: CHRR and OptGP against COBRApy's OptGP.

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

from common import BIN, CORE, ICH360, write_result


class AllWarmupPoints(OptGPSampler):
    """COBRApy drops warmup points it finds redundant, and its chains can only move along the points kept.
    On iCH360 that loses every point with flux through ASNS1, whose split with ASNS2 is then frozen, so the
    reference for that model keeps all the points. `stock_cobrapy_on_ich360` measures the difference."""

    def _is_redundant(self, matrix, cutoff=None):
        return np.zeros(matrix.shape[0], dtype=bool)


RHAT_LIMIT = 1.05
MEAN_LIMIT = 0.15


def run_ours(path, method: str, thinning: int) -> dict:
    out = subprocess.run(
        [BIN, "sample", path, "--method", method, "--chains", "4", "--samples", "2500", "--thinning", str(thinning),
         "--burn-in", "20000", "--seed", "11", "--fraction", "0.9"],
        check=True, capture_output=True, text=True,
    )
    return json.loads(out.stdout)


def at_biomass_floor(path, biomass: str):
    model = cobra.io.read_sbml_model(str(path))
    optimum = model.slim_optimize()
    model.reactions.get_by_id(biomass).lower_bound = 0.9 * optimum
    return model


def stock_cobrapy_on_ich360(chrr: dict) -> dict:
    """What COBRApy's own warmup filter does to ASNS1 and ASNS2, against CHRR's means."""
    sampler = OptGPSampler(at_biomass_floor(ICH360, "Biomass"), processes=4, thinning=200, seed=11)
    samples = sampler.sample(1500)
    samples = samples[sampler.validate(samples.to_numpy()) == "v"]
    means = dict(zip(chrr["reactions"], chrr["diagnostics"]["mean"]))
    return {
        "valid_samples": int(len(samples)),
        "asns1_largest_flux": float(samples["ASNS1"].abs().max()),
        "asns1_mean_chrr": means["ASNS1"],
        "asns2_distance_from_chrr_in_flux_sd": abs(means["ASNS2"] - float(samples["ASNS2"].mean())) / float(samples["ASNS2"].std()),
    }


def compare(path, biomass: str, chrr_moves: int, optgp_moves: int, reference_moves: int, sampler_class=OptGPSampler) -> dict:
    # OptGP mixes more slowly per move, so it gets more thinning to reach a comparable sample size.
    chrr = run_ours(path, "chrr", chrr_moves)
    optgp = run_ours(path, "optgp", optgp_moves)
    sampler = sampler_class(at_biomass_floor(path, biomass), processes=4, thinning=reference_moves, seed=11)
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
    stock = stock_cobrapy_on_ich360(chrr) if sampler_class is AllWarmupPoints else None
    passed = chrr["diagnostics"]["max_split_rhat"] <= RHAT_LIMIT and all(v[0] <= MEAN_LIMIT for v in worst.values())
    print(chrr["diagnostics"]["max_split_rhat"], chrr["diagnostics"]["min_ess"], optgp["diagnostics"]["max_split_rhat"], optgp["diagnostics"]["min_ess"], worst)
    return {
        **({"stock_cobrapy": stock} if stock else {}),
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
    }


def main() -> int:
    result = compare(CORE, "BIOMASS_Ecoli_core_w_GAM", 100, 1000, 100)
    # iCH360 has about twice the free dimensions, so every sampler gets more moves between samples.
    result["ich360"] = compare(ICH360, "Biomass", 500, 20000, 2000, AllWarmupPoints)
    result["pass"] = result["pass"] and result["ich360"]["pass"]
    write_result("v06_sampling", result)
    return 0 if result["pass"] else 1


if __name__ == "__main__":
    sys.exit(main())
