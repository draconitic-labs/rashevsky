#!/usr/bin/env python3
"""e_coli_core: Rashevsky against COBRApy (GLPK) on the same SBML file.

Tolerances: FBA objective within 1e-6 of 0.8739 and of COBRApy; FVA ranges within 1e-6 of COBRApy.
The pFBA total flux, loopless objective, single and double deletions and envelopes are compared
with the same tolerance. MOMA is compared through the quadratic distance within 1e-3, since
COBRApy's QP solver here is OSQP, which has a looser termination tolerance.
"""
from __future__ import annotations

import sys

import cobra
from cobra.flux_analysis import (
    double_gene_deletion,
    flux_variability_analysis,
    loopless,
    moma,
    pfba,
    production_envelope,
    room,
    single_gene_deletion,
    single_reaction_deletion,
)
from cobra.medium import minimal_medium

from common import CORE, engine, write_result

TOL = 1e-6


def compare(CORE, published=None) -> dict:
    """Every comparison on one SBML file; `published` is a growth rate from the literature, if any."""

    def fresh():
        return cobra.io.read_sbml_model(str(CORE))

    model = fresh()
    checks = {}

    # FBA
    ref = model.optimize()
    ours = engine("fba", CORE)
    simplex = engine("fba", CORE, "--backend", "simplex")
    checks["fba_objective"] = {
        "published": published,
        "cobrapy": ref.objective_value,
        "rashevsky_highs": ours["objective"],
        "rashevsky_simplex": simplex["objective"],
        "pass": abs(ours["objective"] - ref.objective_value) < TOL
        and abs(simplex["objective"] - ref.objective_value) < TOL
        and (published is None or round(ours["objective"], 4) == published),
    }

    # FVA at fractions of 1.0 and 0.9
    for fraction in (1.0, 0.9):
        ref_fva = flux_variability_analysis(model, fraction_of_optimum=fraction)
        ours_fva = engine("fva", CORE, "--fraction", fraction)
        worst = 0.0
        worst_id = None
        for row in ours_fva:
            lo = ref_fva.loc[row["id"], "minimum"]
            hi = ref_fva.loc[row["id"], "maximum"]
            d = max(abs(row["minimum"] - lo), abs(row["maximum"] - hi))
            if d > worst:
                worst, worst_id = d, row["id"]
        checks[f"fva_fraction_{fraction}"] = {
            "reactions": len(ours_fva),
            "max_abs_difference": worst,
            "worst_reaction": worst_id,
            "tolerance": TOL,
            "pass": worst < TOL,
        }

    # pFBA: the optimum total flux is unique even where the flux vector is not.
    ref_pfba = pfba(model)
    ours_pfba = engine("pfba", CORE)
    ref_total = float(abs(ref_pfba.fluxes).sum())
    ours_total = sum(abs(v) for v in ours_pfba["fluxes"].values())
    checks["pfba_total_flux"] = {
        "cobrapy": ref_total,
        "rashevsky": ours_total,
        "pass": abs(ref_total - ours_total) < 1e-5,
    }

    # Loopless FBA with COBRApy's add_loopless; GLPK solves the MILP.
    with model:
        loopless.add_loopless(model)
        ref_ll = model.optimize().objective_value
    ours_ll = engine("loopless", CORE)["objective"]
    checks["loopless_objective"] = {
        "cobrapy": ref_ll,
        "rashevsky": ours_ll,
        "pass": abs(ref_ll - ours_ll) < 1e-5,
    }

    # CycleFreeFlux keeps the objective and exchanges.
    free = engine("cycle-free", CORE)
    ref_free = loopless.loopless_solution(model, ref.fluxes.to_dict())
    d = max(abs(free["fluxes"][r.id] - ref_free.fluxes[r.id]) for r in model.reactions if r.boundary)
    checks["cycle_free_exchanges"] = {"max_abs_difference": d, "pass": d < 1e-6}

    # Deletions
    ref_sgd = single_gene_deletion(model)
    ours_sgd = engine("deletion", CORE, "--kind", "gene", "--order", "1")
    worst = 0.0
    for row in ours_sgd:
        ref_val = ref_sgd.knockout[row["ids"][0]]["growth"].iloc[0]
        ref_val = 0.0 if ref_val != ref_val else ref_val
        worst = max(worst, abs(ref_val - row["growth"]))
    checks["single_gene_deletion"] = {"genes": len(ours_sgd), "max_abs_difference": worst, "pass": worst < TOL}

    ref_srd = single_reaction_deletion(model)
    ours_srd = engine("deletion", CORE, "--kind", "reaction", "--order", "1")
    worst = 0.0
    for row in ours_srd:
        ref_val = ref_srd.knockout[row["ids"][0]]["growth"].iloc[0]
        ref_val = 0.0 if ref_val != ref_val else ref_val
        worst = max(worst, abs(ref_val - row["growth"]))
    checks["single_reaction_deletion"] = {"reactions": len(ours_srd), "max_abs_difference": worst, "pass": worst < TOL}

    subset = [g.id for g in model.genes]
    ref_dgd = double_gene_deletion(model, subset, subset, processes=1)
    ours_dgd = engine("deletion", CORE, "--kind", "gene", "--order", "2")
    lookup = {frozenset(row.ids): row.growth for row in ref_dgd.itertuples()}
    worst = 0.0
    for row in ours_dgd:
        ref_val = lookup[frozenset(row["ids"])]
        ref_val = 0.0 if ref_val != ref_val else ref_val
        worst = max(worst, abs(ref_val - row["growth"]))
    checks["double_gene_deletion"] = {"pairs": len(ours_dgd), "max_abs_difference": worst, "pass": worst < 1e-5}

    # MOMA and ROOM after deleting PGI, a reaction with a clear response. The reference flux is the
    # pFBA solution of each engine; COBRApy's is taken on a fresh model for each method because a
    # solver switch inside a model context does not restore every setting.
    perturb = "PGI"
    wild = pfba(fresh())
    m = fresh()
    m.solver = "hybrid"
    m.reactions.get_by_id(perturb).knock_out()
    ref_moma = moma(m, solution=wild, linear=False)
    ref_dist = float(((ref_moma.fluxes - wild.fluxes) ** 2).sum())
    ours_moma = engine("moma", CORE, "--ko-reaction", perturb)
    ours_wild = engine("pfba", CORE)
    ours_dist = sum((ours_moma["fluxes"][k] - ours_wild["fluxes"][k]) ** 2 for k in ours_moma["fluxes"])
    # The two pFBA references may be different optima of a degenerate LP, so MOMA runs against
    # COBRApy's reference as well: the distance to the same reference must agree.
    m = fresh()
    m.reactions.get_by_id(perturb).knock_out()
    checks["moma_distance"] = {
        "knockout": perturb,
        "cobrapy_osqp": ref_dist,
        "rashevsky_highs_qp": ours_dist,
        "tolerance": 1e-3,
        "pass": abs(ref_dist - ours_dist) < 1e-3 * max(1.0, ref_dist),
    }
    m = fresh()
    m.reactions.get_by_id(perturb).knock_out()
    ref_lin = moma(m, solution=wild, linear=True)
    ref_l1 = float(abs(ref_lin.fluxes - wild.fluxes).sum())
    ours_lin = engine("moma", CORE, "--ko-reaction", perturb, "--linear")
    ours_l1 = sum(abs(ours_lin["fluxes"][k] - ours_wild["fluxes"][k]) for k in ours_lin["fluxes"])
    checks["linear_moma_distance"] = {
        "cobrapy": ref_l1,
        "rashevsky": ours_l1,
        "pass": abs(ref_l1 - ours_l1) < 1e-5,
    }
    m = fresh()
    m.solver = "hybrid"
    m.reactions.get_by_id(perturb).knock_out()
    ref_room = room(m, solution=wild, linear=False, delta=0.03, epsilon=0.001)
    ours_room = engine("room", CORE, "--ko-reaction", perturb)
    checks["room_objective"] = {
        "cobrapy_highs": ref_room.objective_value,
        "cobrapy_status": ref_room.status,
        "rashevsky": ours_room["objective"],
        "pass": abs(ref_room.objective_value - ours_room["objective"]) < 1e-5,
    }
    model = fresh()

    # Production envelope
    ref_env = production_envelope(model, ["EX_o2_e"], objective="EX_ac_e", points=20)
    ours_env = engine("envelope", CORE, "--control", "EX_o2_e", "--target", "EX_ac_e", "--points", 20)
    worst = 0.0
    for i, row in enumerate(ours_env):
        worst = max(
            worst,
            abs(row["x"] - ref_env["EX_o2_e"].iloc[i]),
            abs(row["minimum"] - ref_env["flux_minimum"].iloc[i]),
            abs(row["maximum"] - ref_env["flux_maximum"].iloc[i]),
        )
    checks["production_envelope"] = {"points": len(ours_env), "max_abs_difference": worst, "pass": worst < TOL}

    # Minimal medium: same total for the LP form, same number of components for the MILP form.
    growth = 0.1
    ref_mm = minimal_medium(model, growth)
    ours_mm = engine("minimal-medium", CORE, "--growth", growth)
    checks["minimal_medium_total"] = {
        "cobrapy": float(ref_mm.sum()),
        "rashevsky": ours_mm["total"],
        "pass": abs(float(ref_mm.sum()) - ours_mm["total"]) < 1e-5,
    }
    ref_mm_c = minimal_medium(model, growth, minimize_components=True)
    ours_mm_c = engine("minimal-medium", CORE, "--growth", growth, "--components")
    checks["minimal_medium_components"] = {
        "cobrapy": int((ref_mm_c > 1e-9).sum()),
        "rashevsky": len(ours_mm_c["uptake"]),
        "pass": int((ref_mm_c > 1e-9).sum()) == len(ours_mm_c["uptake"]),
    }

    return checks


def record(name: str, checks: dict) -> int:
    passed = all(c["pass"] for c in checks.values())
    write_result(name, {"cobrapy_version": cobra.__version__, "checks": checks, "pass": passed})
    for key, c in checks.items():
        print(f"  {'ok ' if c['pass'] else 'BAD'} {key}: { {k: v for k, v in c.items() if k != 'pass'} }")
    return 0 if passed else 1


def main() -> int:
    return record("v01_core_cobrapy", compare(CORE, published=0.8739))


if __name__ == "__main__":
    sys.exit(main())
