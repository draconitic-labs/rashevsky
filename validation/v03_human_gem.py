#!/usr/bin/env python3
"""Human-GEM: load, element and charge balance, and the 57 essential metabolic tasks.

References: Robinson et al., Sci. Signal. 13 (2020) eaaz1482, which defines the essential task list
and requires the model to pass all of it; the list and the model are read from the Human-GEM
repository. Cross-checks: the Python port of RAVEN's `check_tasks` on the same SBML file and task
list, and COBRApy's `check_mass_balance` per reaction.

Pass: every essential task passes in both engines with identical verdicts, the verification list
agrees with RAVEN, and the set of imbalanced reactions equals COBRApy's.
"""
from __future__ import annotations

import sys
import time

import cobra
from raven_toolbox.tasks import check_tasks

from common import DATA, engine, write_result

MODEL = DATA / "Human-GEM.xml"
ESSENTIAL = DATA / "metabolicTasks_Essential.txt"
VERIFY = DATA / "metabolicTasks_VerifyModel.txt"


def main() -> int:
    t0 = time.time()
    info = engine("info", MODEL)
    load_seconds = time.time() - t0

    t0 = time.time()
    ours_essential = engine("tasks", MODEL, ESSENTIAL)
    ours_seconds = time.time() - t0
    ours_verify = engine("tasks", MODEL, VERIFY)

    model = cobra.io.read_sbml_model(str(MODEL))
    model.solver = "glpk"
    t0 = time.time()
    ref_essential = check_tasks(model, str(ESSENTIAL), close_boundaries=True)
    ref_seconds = time.time() - t0
    ref_verify = check_tasks(model, str(VERIFY), close_boundaries=True)

    def verdicts(results):
        return [(r.id, r.passed) for r in results]

    ours_v = [(r["id"], r["passed"]) for r in ours_essential]
    ours_v2 = [(r["id"], r["passed"]) for r in ours_verify]
    essential_failed = [r["id"] for r in ours_essential if not r["passed"]]

    # Balance
    balance = engine("balance", MODEL)
    ref_imbalanced = {}
    for r in model.reactions:
        if len(r.metabolites) == 1:
            continue
        try:
            mb = r.check_mass_balance()
        except ValueError:
            continue
        mb = {k: v for k, v in mb.items() if abs(v) > 1e-6}
        if mb:
            ref_imbalanced[r.id] = mb
    only_ours = sorted(set(balance["imbalanced"]) - set(ref_imbalanced))
    # COBRApy counts a metabolite without a formula as empty and reports the reaction as imbalanced;
    # this engine reports such a reaction as undetermined. Those are the only differences allowed.
    undetermined = set(balance["formula_missing"])
    only_ref_all = set(ref_imbalanced) - set(balance["imbalanced"])
    only_ref = sorted(only_ref_all - undetermined)

    internal = info["reactions"] - info["exchanges"]
    passed = (
        len(ours_essential) == 57
        and not essential_failed
        and ours_v == verdicts(ref_essential)
        and ours_v2 == verdicts(ref_verify)
        and not only_ours
        and not only_ref
    )
    write_result("v03_human_gem", {
        "model": info,
        "load_seconds": load_seconds,
        "essential_tasks": len(ours_essential),
        "essential_passed": sum(r["passed"] for r in ours_essential),
        "essential_failed": essential_failed,
        "essential_agree_with_raven": ours_v == verdicts(ref_essential),
        "raven_essential_passed": sum(r.passed for r in ref_essential),
        "verification_tasks": len(ours_verify),
        "verification_passed": sum(r["passed"] for r in ours_verify),
        "verification_agree_with_raven": ours_v2 == verdicts(ref_verify),
        "raven_verification_passed": sum(r.passed for r in ref_verify),
        "seconds_essential_rashevsky": ours_seconds,
        "seconds_essential_raven_glpk": ref_seconds,
        "internal_reactions": internal,
        "balanced_reactions": balance["balanced"],
        "formula_missing": len(undetermined),
        "formula_missing_reported_imbalanced_by_cobrapy": len(only_ref_all & undetermined),
        "imbalanced_reactions": len(balance["imbalanced"]),
        "balanced_percent_of_internal": 100.0 * balance["balanced"] / internal,
        "imbalanced_only_in_rashevsky": only_ours[:30],
        "imbalanced_only_in_cobrapy": only_ref[:30],
        "pass": passed,
    })
    print(
        f"essential {sum(r['passed'] for r in ours_essential)}/{len(ours_essential)} "
        f"(RAVEN {sum(r.passed for r in ref_essential)}), verification "
        f"{sum(r['passed'] for r in ours_verify)}/{len(ours_verify)} (RAVEN {sum(r.passed for r in ref_verify)}); "
        f"balance {balance['balanced']}/{internal}, imbalanced {len(balance['imbalanced'])}, "
        f"undetermined {len(undetermined)} ({len(only_ref_all & undetermined)} of them imbalanced in COBRApy), other differences {len(only_ours)}+{len(only_ref)}; "
        f"seconds ours {ours_seconds:.1f} RAVEN {ref_seconds:.1f}"
    )
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
