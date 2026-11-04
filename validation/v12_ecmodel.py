#!/usr/bin/env python3
"""Enzyme constraints against the published GECKO model eciML1515 (ecModels, batch version).

The same optimum is reached four ways on seven conditions, all from the numbers in the published file:

1. The engine's enzyme module on the plain network: `_REV` and isozyme reactions folded back, with the turnover
   numbers, molecular weights and protein pool given as data.
2. The engine reading the GECKO file directly, where enzymes are metabolites.
3. HiGHS on the LP with the enzyme rows eliminated by algebra, leaving one protein-mass row with well scaled
   coefficients. Its dual vector gives an upper bound on the true optimum that holds for any dual vector, because
   every variable has finite bounds, so this value is the reference once the bound closes on it.
4. COBRApy with GLPK on the GECKO file.

COBRApy cannot read the file as published: it decodes `R_protein__32__pseudoreaction` to an identifier with a
space and optlang rejects it, so COBRApy gets a copy with that one encoding replaced.

Pass: routes 1, 2 and 4 within 1e-6 relative of the certified optimum, and the certificate closed to 1e-8 relative.
"""
from __future__ import annotations

import json
import pathlib
import re
import sys
import tempfile

import cobra
import highspy
import numpy as np
import scipy.sparse as sp
from cobra.util import create_stoichiometric_matrix

from common import DATA, engine, write_result

EC = DATA / "ecmodels" / "eciML1515_batch.xml"
NAME = re.compile(r"^(?P<base>.+?)(?P<rev>_REV)?(?:No\d+)?$")
GLUCOSE = "EX_glc__D_e"
OXYGEN = "EX_o2_e"
POOL = "prot_pool_exchange"
DROP_FROM_PLAIN = ("prot_", "pmet_")


def is_protein(metabolite_id: str) -> bool:
    return metabolite_id.startswith("prot_") and not metabolite_id.startswith("prot_pool")


def stoichiometry(reaction: cobra.Reaction) -> dict:
    return {m: c for m, c in reaction.metabolites.items() if not m.id.startswith(DROP_FROM_PLAIN)}


def folded_network(ec: cobra.Model) -> tuple[cobra.Model, list, list, dict]:
    """The published network without its enzyme machinery, and the enzyme data that machinery encoded."""
    weights = {}
    for r in ec.reactions:
        if r.id.startswith("draw_prot_"):
            target = next(m for m in r.metabolites if is_protein(m.id))
            pool = next(m for m in r.metabolites if m.id.startswith("prot_pool"))
            weights[target.id] = -r.metabolites[pool]

    forms: dict[str, dict] = {}
    catalysis = []
    for r in ec.reactions:
        if r.id.startswith("draw_prot_") or r.id == POOL:
            continue
        if r.id.startswith("arm_"):
            forms.setdefault(r.id[4:], {})["arm"] = r
            continue
        m = NAME.match(r.id)
        slot = forms.setdefault(m["base"], {}).setdefault("reverse" if m["rev"] else "forward", [])
        slot.append(r)
        subunits = {x.id: -c for x, c in r.metabolites.items() if is_protein(x.id)}
        if subunits:
            # The file's coefficient on a subunit is copies / kcat. The largest sets the turnover number.
            turnover = 1.0 / max(subunits.values())
            forward = m["rev"] is None
            catalysis.append({
                "reaction": m["base"],
                "subunits": [[k, c * turnover] for k, c in sorted(subunits.items())],
                "kcat_forward": turnover if forward else 0.0,
                "kcat_reverse": 0.0 if forward else turnover,
            })

    def direction(base: str, key: str, name: str) -> tuple[dict, float, float] | None:
        form = forms.get(base, {})
        arm = forms.get(name, {}).get("arm")
        members = form.get(key)
        if not members and arm is None:
            return None
        s = stoichiometry(members[0]) if members else {}
        low, high = (members[0].lower_bound, members[0].upper_bound) if members else (0.0, 0.0)
        if arm is not None:
            s = {**stoichiometry(arm), **s}
            low, high = arm.lower_bound, arm.upper_bound
        return s, low, high

    plain = cobra.Model("eciML1515_folded")
    report = {"reverse_not_the_opposite_of_forward": [], "one_direction_costed": []}
    for base in sorted(k for k in forms if not k.endswith("_REV")):
        forward = direction(base, "forward", base)
        reverse = direction(base, "reverse", f"{base}_REV")
        if forward is None and reverse is None:
            continue
        reaction = cobra.Reaction(base)
        if forward is not None:
            s, low, high = forward
        else:
            s = {m: -c for m, c in reverse[0].items()}
            low, high = 0.0, 0.0
        if reverse is not None:
            low = -reverse[2]
            if forward is not None and any(abs(s.get(m, 0.0) + reverse[0].get(m, 0.0)) > 1e-9 for m in set(s) | set(reverse[0])):
                report["reverse_not_the_opposite_of_forward"].append(base)
        reaction.add_metabolites({cobra.Metabolite(m.id, compartment=m.compartment): c for m, c in s.items()})
        reaction.bounds = (low, high)
        plain.add_reactions([reaction])
        costed = {("forward" if c["kcat_forward"] else "reverse") for c in catalysis if c["reaction"] == base}
        if len(costed) == 1 and reverse is not None and forward is not None:
            report["one_direction_costed"].append(base)
    objective = next(r for r in ec.reactions if r.objective_coefficient)
    plain.objective = objective.id
    enzymes = [{"id": k, "molecular_weight": w} for k, w in sorted(weights.items())]
    return plain, enzymes, catalysis, report


def certified_optimum(ec: cobra.Model, edits: dict) -> tuple[float, float, float]:
    """Optimum with the enzyme rows eliminated, its dual upper bound, and the protein-mass row's activity."""
    S = create_stoichiometric_matrix(ec, array_type="lil").tocoo()
    reactions = [r.id for r in ec.reactions]
    metabolites = [m.id for m in ec.metabolites]
    low = np.array([r.lower_bound for r in ec.reactions])
    high = np.array([r.upper_bound for r in ec.reactions])
    for name, value in edits.items():
        high[reactions.index(name)] = value
    objective = np.array([r.objective_coefficient for r in ec.reactions])
    pool_row = next(i for i, m in enumerate(metabolites) if m.startswith("prot_pool"))
    weight = {}
    for i, j, v in zip(S.row, S.col, S.data):
        if reactions[j].startswith("draw_prot_") and i != pool_row:
            column = S.col == j
            weight[i] = -S.data[column & (S.row == pool_row)][0] / v
    cost = np.zeros(len(reactions))
    for i, j, v in zip(S.row, S.col, S.data):
        if i in weight and not reactions[j].startswith("draw_prot_"):
            cost[j] -= v * weight[i]
    columns = [j for j, n in enumerate(reactions) if not n.startswith("draw_prot_") and n != POOL]
    rows = [i for i in range(len(metabolites)) if i not in weight and i != pool_row]
    pool = high[reactions.index(POOL)]
    A = sp.vstack([S.tocsr()[rows, :].tocsc()[:, columns], sp.csr_matrix(cost[columns] / pool)]).tocsc()
    row_high = np.r_[np.zeros(len(rows)), 1.0]
    row_low = np.r_[np.zeros(len(rows)), -np.inf]
    solver = highspy.Highs()
    solver.setOptionValue("output_flag", False)
    solver.setOptionValue("small_matrix_value", 1e-12)
    solver.setOptionValue("primal_feasibility_tolerance", 1e-10)
    solver.setOptionValue("dual_feasibility_tolerance", 1e-10)
    lp = highspy.HighsLp()
    lp.num_col_, lp.num_row_ = A.shape[1], A.shape[0]
    lp.col_cost_, lp.col_lower_, lp.col_upper_ = -objective[columns], low[columns], high[columns]
    lp.row_lower_, lp.row_upper_ = row_low, row_high
    lp.a_matrix_.format_ = highspy.MatrixFormat.kColwise
    lp.a_matrix_.start_, lp.a_matrix_.index_, lp.a_matrix_.value_ = A.indptr, A.indices, A.data
    solver.passModel(lp)
    solver.run()
    solution = solver.getSolution()
    optimum = -solver.getInfo().objective_function_value
    # For y >= 0 on the mass row, max c'x <= y_pool + sum_j max over x_j of (c_j - A_j'y) x_j, for any other y.
    best = np.inf
    for sign in (1.0, -1.0):
        y = sign * np.array(solution.row_dual)
        if y[-1] < 0:
            y = y.copy()
            y[-1] = 0.0
        reduced = objective[columns] - A.T @ y
        best = min(best, float(np.sum(np.maximum(reduced * low[columns], reduced * high[columns]))) + y[-1])
    activity = float((A @ np.array(solution.col_value))[-1]) * pool
    return optimum, best, activity


def main() -> int:
    with tempfile.TemporaryDirectory() as scratch:
        scratch = pathlib.Path(scratch)
        clean = scratch / "eciML1515_batch_clean.xml"
        clean.write_text(EC.read_text().replace("__32__", "_"))
        ec = cobra.io.read_sbml_model(str(clean))
        plain, enzymes, catalysis, report = folded_network(ec)
        pool = ec.reactions.get_by_id(POOL).upper_bound
        plain_path = scratch / "plain.xml"
        cobra.io.write_sbml_model(plain, str(plain_path))
        table = scratch / "enzymes.json"
        table.write_text(json.dumps({"enzymes": enzymes, "catalysis": catalysis, "protein_pool": pool}))
        shared = sum(1 for e in enzymes if sum(1 for c in catalysis if any(s[0] == e["id"] for s in c["subunits"])) > 1)
        complexes = sum(1 for c in catalysis if len(c["subunits"]) > 1)
        print(len(plain.reactions), "folded reactions,", len(enzymes), "enzymes,", shared, "shared,", len(catalysis), "catalyses,", complexes, "complexes, pool", pool, report, flush=True)

        conditions = [
            ("batch, unlimited glucose", 1.0, None, True),
            ("glucose uptake 10", 1.0, 10.0, True),
            ("glucose uptake 5", 1.0, 5.0, True),
            ("glucose uptake 1", 1.0, 1.0, True),
            ("half the protein pool", 0.5, None, True),
            ("twice the protein pool", 2.0, None, True),
            ("no oxygen", 1.0, None, False),
        ]
        rows = []
        for label, factor, glucose, oxygen in conditions:
            edits = {POOL: pool * factor}
            edits_plain, edits_gecko = [], ["--bound", f"{POOL}=0:{pool * factor}"]
            if glucose is not None:
                edits[f"{GLUCOSE}_REV"] = glucose
                edits_plain += ["--bound", f"{GLUCOSE}=-{glucose}:{plain.reactions.get_by_id(GLUCOSE).upper_bound}"]
                edits_gecko += ["--bound", f"{GLUCOSE}_REV=0:{glucose}"]
            if not oxygen:
                edits[f"{OXYGEN}_REV"] = 0.0
                edits_plain += ["--bound", f"{OXYGEN}=0:{plain.reactions.get_by_id(OXYGEN).upper_bound}"]
                edits_gecko += ["--bound", f"{OXYGEN}_REV=0:0"]
            with ec:
                for name, value in edits.items():
                    ec.reactions.get_by_id(name).upper_bound = value
                glpk = ec.slim_optimize()
            optimum, bound, mass = certified_optimum(ec, edits)
            module = engine("ecfba", plain_path, table, "--pool", pool * factor, *edits_plain)
            direct = engine("fba", EC, *edits_gecko)
            module_objective = module["fluxes"]["objective"]
            row = {
                "condition": label,
                "certified_optimum": optimum,
                "certificate_upper_bound": bound,
                "certificate_gap": (bound - optimum) / optimum,
                "engine_enzyme_module": module_objective,
                "engine_gecko_file": direct["objective"],
                "cobrapy_glpk_gecko_file": glpk,
                "module_relative_difference": abs(module_objective - optimum) / optimum,
                "direct_relative_difference": abs(direct["objective"] - optimum) / optimum,
                "cobrapy_relative_difference": abs(glpk - optimum) / optimum,
                "protein_mass_used": module["pool_used"],
                "protein_mass_limit": pool * factor,
            }
            row["pass"] = (
                row["module_relative_difference"] < 1e-6
                and row["direct_relative_difference"] < 1e-6
                and abs(row["certificate_gap"]) < 1e-8
                and row["cobrapy_relative_difference"] < 1e-6
            )
            rows.append(row)
            print(f"{label}: optimum {optimum:.10f} bound {bound:.10f} module {module_objective:.10f} ({row['module_relative_difference']:.1e}) direct {direct['objective']:.10f} ({row['direct_relative_difference']:.1e}) glpk {glpk:.10f} ({row['cobrapy_relative_difference']:.1e})", flush=True)
    passed = all(r["pass"] for r in rows) and not report["reverse_not_the_opposite_of_forward"] and not report["one_direction_costed"]
    write_result("v12_ecmodel", {
        "pass": passed,
        "cobrapy_version": cobra.__version__,
        "source_commit": "8c50008decb93a3dbf31949836ff61eec57d572b",
        "folded_reactions": len(plain.reactions),
        "enzymes": len(enzymes),
        "shared_enzymes": shared,
        "catalyses": len(catalysis),
        "complexes": complexes,
        "structure_report": report,
        "protein_pool": pool,
        "conditions": rows,
    })
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
