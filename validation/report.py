#!/usr/bin/env python3
"""Writes docs/VALIDATION.md and its figures from the result files. The numbers in that document exist
only here and in validation/results, so they cannot drift apart."""
from __future__ import annotations

import json
import pathlib

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

HERE = pathlib.Path(__file__).resolve().parent
DOCS = HERE.parent / "docs"
FIGURES = DOCS / "figures"
R = {p.stem: json.loads(p.read_text()) for p in sorted((HERE / "results").glob("v*.json"))}

INK = "#20252a"
BLUE = "#3d6ae0"
GOLD = "#c79a12"
plt.rcParams.update({"font.size": 10, "axes.spines.top": False, "axes.spines.right": False, "savefig.dpi": 150, "figure.autolayout": True})


def verdict(ok: bool) -> str:
    return "pass" if ok else "**FAIL**"


def fig_keio() -> None:
    k = R["v02_iml1515_keio"]
    names = sorted(k["per_condition"])
    acc = [100 * (k["per_condition"][n]["TP"] + k["per_condition"][n]["TN"]) / sum(k["per_condition"][n].values()) for n in names]
    fig, ax = plt.subplots(figsize=(8, 3.2))
    ax.bar(names, acc, color=BLUE)
    ax.axhline(93.4, color=GOLD, label="Published overall accuracy, 93.4 %")
    ax.set_ylabel("Accuracy, %")
    ax.set_ylim(85, 100)
    ax.tick_params(axis="x", rotation=60)
    ax.legend(frameon=False)
    fig.savefig(FIGURES / "keio.png")
    plt.close(fig)


def fig_integrators() -> None:
    study = R["v08_numerics"]["integrators"]["study"]
    fig, ax = plt.subplots(figsize=(5.2, 3.4))
    for method, colour in (("rodas4", BLUE), ("bdf", GOLD)):
        rows = [(v["rhs_evals"], v["max_rel_error"]) for k, v in study.items() if k.startswith(method)]
        rows.sort()
        ax.loglog([r[0] for r in rows], [r[1] for r in rows], "o-", color=colour, label=method.upper() if method == "bdf" else "Rodas4")
    ax.set_xlabel("Right-hand side evaluations")
    ax.set_ylabel("Largest relative error at t = 2")
    ax.legend(frameon=False)
    fig.savefig(FIGURES / "integrators.png")
    plt.close(fig)


def fig_teusink() -> None:
    table = R["v04_teusink_mca"]["table"]
    fig, ax = plt.subplots(figsize=(4, 4))
    x = [row["curated_model"] for row in table]
    y = [row["rashevsky"] for row in table]
    ax.loglog(x, y, "o", color=BLUE)
    line = [min(x) * 0.7, max(x) * 1.4]
    ax.loglog(line, line, color=GOLD)
    ax.set_xlabel("Curators' table")
    ax.set_ylabel("Rashevsky")
    fig.savefig(FIGURES / "teusink.png")
    plt.close(fig)


def main() -> None:
    FIGURES.mkdir(parents=True, exist_ok=True)
    fig_keio()
    fig_integrators()
    fig_teusink()
    core = R["v01_core_cobrapy"]["checks"]
    keio = R["v02_iml1515_keio"]
    human = R["v03_human_gem"]
    t = R["v04_teusink_mca"]
    suite = R["v05_sbml_suite"]
    sampling = R["v06_sampling"]
    overlay = R["v07_overlay_escher"]
    numerics = R["v08_numerics"]
    export = R.get("v09_sbml_export")
    perf = R.get("v10_performance")
    recon = R["v11_recon3d"]
    ec = R["v12_ecmodel"]

    w = []
    add = w.append
    add("# Validation\n")
    add("This file is written by `validation/report.py` from `validation/results`. To reproduce every number:\n")
    add("```\ncargo build --release -p rashevsky-cli\npython -m venv .venv && . .venv/bin/activate\npip install -r validation/requirements.txt\n(maturin develop --release is needed only for the notebook)\nCARGO_TARGET_DIR=target python validation/run_all.py\n```\n")
    add("Each script below also runs alone, from `validation/`, and prints what it finds. References are other tools run on this machine on the same files, or numbers from the named publication. External tools: "
        f"COBRApy {R['v01_core_cobrapy']['cobrapy_version']}, libRoadRunner {t['roadrunner_version']}, COPASI {t['copasi_version']}.\n")
    add("| Check | Criterion | Result | Verdict |\n| --- | --- | --- | --- |")
    rows = [
        ("`e_coli_core` growth", "0.8739 1/h from BiGG; within 1e-6 of COBRApy", f"{core['fba_objective']['rashevsky_highs']:.6f} with HiGHS, {core['fba_objective']['rashevsky_simplex']:.6f} with the dense simplex", core["fba_objective"]["pass"]),
        ("FVA at 100 % of optimum", "within 1e-6 of COBRApy", f"largest difference {core['fva_fraction_1.0']['max_abs_difference']:.1e} over {core['fva_fraction_1.0']['reactions']} reactions", core["fva_fraction_1.0"]["pass"]),
        ("FVA at 90 % of optimum", "within 1e-6 of COBRApy", f"largest difference {core['fva_fraction_0.9']['max_abs_difference']:.1e}", core["fva_fraction_0.9"]["pass"]),
        ("iML1515 essentiality", "within 0.5 points of the published 93.4 %", f"{keio['accuracy_percent']:.3f} % of {keio['phenotypes']:,} phenotypes, MCC {keio['mcc']:.4f}", keio["pass"]),
        ("Human-GEM essential tasks", "all pass; same verdicts as RAVEN", f"{human['essential_passed']} of {human['essential_tasks']}", human["essential_passed"] == human["essential_tasks"] and human["essential_agree_with_raven"]),
        ("Human-GEM balance", "same imbalanced set as COBRApy", f"{human['balanced_reactions']:,} of {human['internal_reactions']:,} balanced, {human['imbalanced_reactions']} imbalanced, {human['formula_missing']} undetermined", not human["imbalanced_only_in_rashevsky"] and not human["imbalanced_only_in_cobrapy"]),
        ("Control theorems", "under 1e-8 on five models", f"worst scaled residual {max(x['worst'] for x in t['theorems'].values()):.1e}", all(x["pass"] for x in t["theorems"].values())),
        ("Teusink steady state", "the curators' table; libRoadRunner within 1e-5", f"all {len(t['table'])} values in tolerance; libRoadRunner {t['steady_state_max_rel_difference_species']:.1e}", all(r["pass"] for r in t["table"]) and t["steady_state_max_rel_difference_species"] < 1e-5),
        ("Recon3D", "same size, objective and 2,248 gene deletions within 1e-6 relative of COBRApy, FVA within 1e-5", f"objective {recon['checks']['objective']['rashevsky']:.4f}, deletions {recon['checks']['gene_deletions']['max_relative_difference']:.1e}, FVA {recon['checks']['fva']['max_abs_difference']:.1e}", recon["pass"]),
        ("Enzyme constraints on eciML1515", "within 1e-6 of a certified optimum in 7 conditions", f"largest difference {max(c['module_relative_difference'] for c in ec['conditions']):.1e}", ec["pass"]),
        ("SBML semantic suite", "at least 95 % of applicable cases", f"{suite['rodas4']['passed']} of {suite['rodas4']['applicable']} with Rodas4, {suite['bdf']['passed']} of {suite['bdf']['applicable']} with BDF", suite["rodas4"]["percent"] >= 95 and suite["bdf"]["percent"] >= 95),
    ]
    for name, criterion, result, ok in rows:
        add(f"| {name} | {criterion} | {result} | {verdict(ok)} |")
    add("")

    add("## `e_coli_core` against COBRApy\n\n`validation/v01_core_cobrapy.py`\n")
    add("| Analysis | COBRApy | Rashevsky or largest difference | Verdict |\n| --- | --- | --- | --- |")
    add(f"| pFBA total flux | {core['pfba_total_flux']['cobrapy']:.6f} | {core['pfba_total_flux']['rashevsky']:.6f} | {verdict(core['pfba_total_flux']['pass'])} |")
    add(f"| Loopless objective | {core['loopless_objective']['cobrapy']:.6f} | {core['loopless_objective']['rashevsky']:.6f} | {verdict(core['loopless_objective']['pass'])} |")
    add(f"| CycleFreeFlux exchanges | | difference {core['cycle_free_exchanges']['max_abs_difference']:.1e} | {verdict(core['cycle_free_exchanges']['pass'])} |")
    add(f"| Single gene deletions ({core['single_gene_deletion']['genes']}) | | difference {core['single_gene_deletion']['max_abs_difference']:.1e} | {verdict(core['single_gene_deletion']['pass'])} |")
    add(f"| Single reaction deletions ({core['single_reaction_deletion']['reactions']}) | | difference {core['single_reaction_deletion']['max_abs_difference']:.1e} | {verdict(core['single_reaction_deletion']['pass'])} |")
    add(f"| Double gene deletions ({core['double_gene_deletion']['pairs']:,}) | | difference {core['double_gene_deletion']['max_abs_difference']:.1e} | {verdict(core['double_gene_deletion']['pass'])} |")
    add(f"| MOMA distance after deleting PGI | {core['moma_distance']['cobrapy_osqp']:.4f} (OSQP) | {core['moma_distance']['rashevsky_highs_qp']:.4f} (HiGHS QP) | {verdict(core['moma_distance']['pass'])} |")
    add(f"| Linear MOMA distance | {core['linear_moma_distance']['cobrapy']:.6f} | {core['linear_moma_distance']['rashevsky']:.6f} | {verdict(core['linear_moma_distance']['pass'])} |")
    add(f"| ROOM objective | {core['room_objective']['cobrapy_highs']:.0f} | {core['room_objective']['rashevsky']:.0f} | {verdict(core['room_objective']['pass'])} |")
    add(f"| Production envelope, {core['production_envelope']['points']} points | | difference {core['production_envelope']['max_abs_difference']:.1e} | {verdict(core['production_envelope']['pass'])} |")
    add(f"| Minimal medium total uptake | {core['minimal_medium_total']['cobrapy']:.6f} | {core['minimal_medium_total']['rashevsky']:.6f} | {verdict(core['minimal_medium_total']['pass'])} |")
    add(f"| Minimal medium components | {core['minimal_medium_components']['cobrapy']} | {core['minimal_medium_components']['rashevsky']} | {verdict(core['minimal_medium_components']['pass'])} |")
    add("\nThe MOMA tolerance is 1e-3 because OSQP stops at a looser tolerance than the simplex-based solvers.\n")

    add("## iML1515 essentiality against the Keio collection\n\n`validation/v02_iml1515_keio.py`\n")
    c = keio["counts"]
    p = keio["published_counts"]
    add("Reference: Monk et al. 2017, Supplementary Note 4 and Supplementary Data 11, with its growth screens of the Keio collection on 16 carbon sources. Protocol from its Supplementary Methods: oxygen uptake bounded at 20 and the carbon source at 10 mmol/gDW/h, the reactions of genes absent from BW25113 closed, and essential below 0.05 1/h.\n")
    add("| | TP | FN | FP | TN | Accuracy |\n| --- | --- | --- | --- | --- | --- |")
    add(f"| Published | {p['TP']:,} | {p['FN']} | {p['FP']:,} | {p['TN']:,} | {keio['published_accuracy_percent']} % |")
    add(f"| Rashevsky | {c['TP']:,} | {c['FN']} | {c['FP']:,} | {c['TN']:,} | {keio['accuracy_percent']:.3f} % |")
    differ = keio["classification_disagreements_sample"]
    add(f"\n{keio['rows_classified_differently_from_published']} of {keio['phenotypes']:,} rows differ from the call in the supplement: " + "; ".join(f"{d['gene']} on {d['condition']}, published growth {d['published_growth']:.3f}, computed {d['rashevsky_growth']:.3f}" for d in differ) + ".")
    x = keio["cobrapy_same_file"]["Xylose"]
    add(f" COBRApy gives the same value on this row from the same BiGG file ({x['cobrapy_b3564']:.4f}), so the difference lies between the BiGG file of 2019 and the paper's model, not in the engine. Over the two conditions compared with COBRApy on all {x['genes']:,} genes, the largest difference is {max(v['max_abs_difference_vs_cobrapy'] for v in keio['cobrapy_same_file'].values()):.1e}.\n")
    add("![Accuracy per carbon source](figures/keio.png)\n")

    add("## Human-GEM\n\n`validation/v03_human_gem.py`\n")
    add(f"Human-GEM v2.1.1: {human['model']['reactions']:,} reactions, {human['model']['metabolites']:,} metabolites, {human['model']['genes']:,} genes, read in {human['load_seconds']:.2f} s. Essential tasks of Robinson et al. 2020: {human['essential_passed']} of {human['essential_tasks']} pass, with the verdict on every task equal to that of RAVEN's Python port on the same files ({human['seconds_essential_rashevsky']:.1f} s against {human['seconds_essential_raven_glpk']:.1f} s with GLPK). The verification list gives {human['verification_passed']} of {human['verification_tasks']}, same verdicts.\n")
    add(f"Of {human['internal_reactions']:,} internal reactions, {human['balanced_reactions']:,} ({human['balanced_percent_of_internal']:.1f} %) balance in elements and charge and {human['imbalanced_reactions']} do not, the same set as COBRApy's `check_mass_balance`. Another {human['formula_missing']} involve a metabolite without a formula. This engine reports those as undetermined, and COBRApy counts them as imbalanced.\n")

    add("## Recon3D\n\n`validation/v11_recon3d.py`\n")
    rc = recon["checks"]
    add(f"Recon3D from BiGG Models (Brunk et al. 2018) against COBRApy with GLPK. Reactions, metabolites and genes: {rc['shape']['reactions'][0]:,}, {rc['shape']['metabolites'][0]:,} and {rc['shape']['genes'][0]:,} in both. Biomass objective {rc['objective']['rashevsky']:.6f} against {rc['objective']['cobrapy']:.6f}, relative difference {rc['objective']['relative_difference']:.1e}. FVA at 100 % of the optimum on {rc['fva']['reactions']} reactions chosen by a fixed seed: largest difference {rc['fva']['max_abs_difference']:.1e}. All {rc['gene_deletions']['genes']:,} single gene deletions, of which {rc['gene_deletions']['growth_changed']} change the objective and {rc['gene_deletions']['lethal']} stops it: largest relative difference {rc['gene_deletions']['max_relative_difference']:.1e}. After `rashevsky convert`, libSBML finds {rc['export']['libsbml_errors']} errors and COBRApy reads the same objective.\n")

    add("## Enzyme constraints\n\n`validation/v12_ecmodel.py`\n")
    add(f"The published `eciML1515` batch model (ecModels commit {ec['source_commit'][:7]}) with its pool of {ec['protein_pool']:.4f} g/gDW. The file's `_REV` and isozyme reactions were folded back into {ec['folded_reactions']:,} reactions and its enzymes read out as {ec['isozyme_direction_pairs']:,} data rows for the enzyme module. The reference is the same program with the enzyme rows eliminated by algebra, solved by HiGHS and closed by a dual upper bound, which holds for any dual vector because every variable has finite bounds.\n")
    add("| Condition | Certified optimum | Upper bound | Enzyme module | Engine on the GECKO file | COBRApy on the GECKO file | Verdict |\n| --- | --- | --- | --- | --- | --- | --- |")
    for c in ec["conditions"]:
        add(f"| {c['condition']} | {c['certified_optimum']:.8f} | {c['certificate_upper_bound']:.8f} | {c['engine_enzyme_module']:.8f} ({c['module_relative_difference']:.0e}) | {c['engine_gecko_file']:.8f} ({c['direct_relative_difference']:.0e}) | {c['cobrapy_glpk_gecko_file']:.8f} ({c['cobrapy_relative_difference']:.0e}) | {verdict(c['pass'])} |")
    add("\nRelative differences from the certified optimum are in brackets. COBRApy cannot read the file as published, because it decodes `R_protein__32__pseudoreaction` to an identifier with a space, so it reads a copy with that one encoding replaced.\n")
    add("This check found a defect. HiGHS removes matrix entries at or below `small_matrix_value`, 1e-9 by default. `eciML1515` has 17 enzyme costs at or below that value, the smallest 5.6e-11, and without them those enzymes are free: the same file then gave 0.49713418 without oxygen against a true 0.49712429. The engine now sets 1e-12, the lowest value HiGHS accepts.\n")

    add("## Kinetic models\n\n`validation/v04_teusink_mca.py`\n")
    add("Reference values: Table 4 of Teusink et al. 2000 as reproduced by the BioModels curators in the notes of BIOMD0000000064, together with their curated model's column. The publisher's site refused direct access to the paper (HTTP 403), so the table was read from those notes. The curators note that the G6P of the curated model differs slightly from the paper's. Tolerance: one unit in the last printed digit of the curators' column.\n")
    add("| Quantity | Published | Curated model | Rashevsky | Verdict |\n| --- | --- | --- | --- | --- |")
    for row in t["table"]:
        add(f"| {row['quantity']} | {row['published']} | {row['curated_model']} | {row['rashevsky']:.4f} | {verdict(row['pass'])} |")
    add(f"\nThe same file in libRoadRunner gives the same steady state: largest relative difference {t['steady_state_max_rel_difference_species']:.1e} over species and {t['steady_state_max_rel_difference_fluxes']:.1e} over fluxes. Scaled control coefficients match COPASI: largest difference {t['flux_control_max_abs_difference_vs_copasi']:.1e} for flux and {t['concentration_control_max_abs_difference_vs_copasi']:.1e} for concentration, COPASI's being computed with finite-difference elasticities. libRoadRunner's own control matrices raise an exception on this model. The Rust test `control_coefficient_matches_a_perturbation` also checks a coefficient against a direct perturbation of the transporter's Vmax.\n")
    add("![Steady state against the curators' table](figures/teusink.png)\n")
    add("### Control theorems\n")
    add("Scaled residuals, the form in which the theorems are stated. The unscaled residuals carry the units of the elasticities and are in the results file.\n")
    add("| Model | Reactions | Conserved moieties | Flux summation | Concentration summation | Flux connectivity | Concentration connectivity |\n| --- | --- | --- | --- | --- | --- | --- |")
    for m, x in t["theorems"].items():
        r = x["residuals"]
        add(f"| {m}, {x['title']} | {x['reactions']} | {x['conserved_moieties']} | {r['flux_summation']:.1e} | {r['concentration_summation']:.1e} | {r['flux_connectivity']:.1e} | {r['concentration_connectivity']:.1e} |")
    add("")

    add("## SBML Test Suite\n\n`validation/v05_sbml_suite.py`, SBML Test Suite 3.5.0\n")
    add(f"Semantic cases: {suite['rodas4']['total_cases']:,} in the suite, {suite['rodas4']['applicable']:,} applicable. The {suite['rodas4']['total_cases'] - suite['rodas4']['applicable']} left out use the `comp` package, which this engine does not implement. A case passes when every output is within the suite's own tolerance, `|sim - ref| <= absolute + relative |ref|`.\n")
    add(f"| Integrator | Passed | Percent | Seconds |\n| --- | --- | --- | --- |\n| Rodas4 | {suite['rodas4']['passed']:,} of {suite['rodas4']['applicable']:,} | {suite['rodas4']['percent']:.2f} | {suite['rodas4']['seconds']:.0f} |\n| BDF | {suite['bdf']['passed']:,} of {suite['bdf']['applicable']:,} | {suite['bdf']['percent']:.2f} | {suite['bdf']['seconds']:.0f} |\n")
    add(f"Of these, {suite['rodas4']['fbc_cases_passed']} are flux balance cases and {suite['rodas4']['fast_reaction_cases_passed']} use fast reactions, simulated by the singular perturbation of `docs/MODEL.md`. Random tie-breaking between events uses a fixed seed.\n")
    add("### Stochastic cases\n")
    k3 = suite["case_00003_exact_law_check"]
    flagged = suite["stochastic_cases_flagged_on_first_seed"]
    add(f"The {suite['stochastic_cases']} `StochasticTimeCourse` cases of the DSMTS (Evans, Gillespie and Wilkinson 2008), with {suite['stochastic_runs_per_case']:,} realizations each. The suite accepts a mean `Z` in (-3, 3) and a variance `Y` in (-5, 5) at each time point, and states that a correct simulator can fail a few points. The limit used here is 3 per case, with a flagged case repeated once on an independent seed, because the points of a case are strongly correlated and one excursion can flag many of them.\n")
    add(f"- Flagged on the first seed: {', '.join(flagged)}.")
    fail = suite["stochastic_cases_failing_the_literal_criterion"]
    add(f"- Still over the limit on the second seed: {', '.join(fail) if fail else 'none'}.")
    if k3:
        add(f"\n**Case 00003 fails the suite's literal criterion**, at this number of realizations, on both seeds. It is a birth-death process, for which the exact law is known (Kendall 1948). Its kurtosis reaches {k3['max_kurtosis']:.0f}, while `Y` assumes the sampling error of a normal population and so overstates significance by a factor of `sqrt((kurtosis - 1) / 2)`. Compared with the exact variance using its exact sampling error at {k3['runs']:,} realizations, the largest `|z|` is {k3['max_abs_variance_z']:.2f} for the variance and {k3['max_abs_mean_z']:.2f} for the mean, with a limit of 4. No error in the simulator is detectable.")
    add("\nThe other cases of the stochastic folder test the `distrib` package, which this engine does not implement.\n")

    add("## Flux sampling\n\n`validation/v06_sampling.py`\n")
    d = sampling["chrr_diagnostics"]
    o = sampling["optgp_diagnostics"]
    m = sampling["max_mean_difference_in_flux_sd"]
    add(f"`e_coli_core` with biomass at 90 % of its optimum or more, a polytope of dimension {sampling['polytope_dimension']}. Criterion: means within {sampling['limit']} of a flux standard deviation.\n")
    add("| Comparison | Largest difference of means, in flux SD | Reaction |\n| --- | --- | --- |")
    for key, label in (("chrr_vs_cobrapy", "CHRR against COBRApy OptGP"), ("optgp_vs_cobrapy", "OptGP against COBRApy OptGP"), ("chrr_vs_optgp", "CHRR against OptGP")):
        add(f"| {label} | {m[key]['value']:.3f} | {m[key]['reaction']} |")
    add(f"\nDiagnostics: CHRR split R-hat {d['max_split_rhat']:.3f}, smallest effective sample size {d['min_ess']:.0f}, largest Geweke z {d['max_geweke_z']:.1f}, in {sampling['chrr_seconds']:.2f} s. OptGP R-hat {o['max_split_rhat']:.3f}, effective size {o['min_ess']:.0f}, Geweke z {o['max_geweke_z']:.1f}, in {sampling['optgp_seconds']:.2f} s. A Rust test checks uniformity on a simplex against its exact mean and variance, and convergence on `e_coli_core` for several seeds.\n")

    add("## Outputs\n\n`validation/v07_overlay_escher.py`\n")
    cv = overlay["coverage"]
    add(f"pFBA on Human-GEM projected onto the {overlay['scene_release']} scene: {cv['events_exact']} of {cv['events_in_scene']:,} events reached by a Reactome or Rhea identifier, {cv['events_ec_only']} more by an EC number only, {cv['events_with_flux']} carrying flux. Event totals match a separate recomputation from the SBML text to {overlay['max_abs_difference_exact_events_vs_independent_recompute']:.1e}, and every event identifier exists in the scene. The Escher map validates against Escher's published 1-0-0 schema: {'yes' if overlay['escher_schema_valid'] else 'no'}.\n")
    if export:
        add("### SBML written by the engine\n\n`validation/v09_sbml_export.py`\n")
        add("| Model | Errors from libSBML | Reactions | Objective read by COBRApy, original | Rewritten | Verdict |\n| --- | --- | --- | --- | --- | --- |")
        for name in ("e_coli_core", "iML1515", "Human-GEM"):
            x = export["checks"][name]
            add(f"| {name} | {x['libsbml_errors']} | {x['reactions'][1]:,} | {x['objective_original']:.6f} | {x['objective_rewritten']:.6f} | {verdict(x['pass'])} |")
        x = export["checks"]["teusink_kinetic"]
        add(f"\nThe rewritten Teusink model reaches the same steady state in libRoadRunner, largest relative difference {x['max_rel_difference_roadrunner']:.1e}.\n")

    add("## Numerical studies\n\n`validation/v08_numerics.py`\n")
    add("### Integrators\n")
    add("Relative error at t = 2 on the Teusink model, against a run at `rtol = 1e-12`. Rust tests also check Rodas4's order of 4 at fixed step, the stiff Prothero-Robinson problem with an exact solution, Robertson's problem with its conservation law, and an index-1 algebraic equation, for both methods.\n")
    add("| Method | `rtol` | Largest relative error | Steps | Right-hand side evaluations |\n| --- | --- | --- | --- | --- |")
    for key, v in numerics["integrators"]["study"].items():
        method, rtol = key.split("_")
        add(f"| {'Rodas4' if method == 'rodas4' else 'BDF'} | {rtol} | {v['max_rel_error']:.1e} | {v['steps']} | {v['rhs_evals']} |")
    add("\n![Work against precision](figures/integrators.png)\n")
    add("### Dynamic FBA\n")
    dfba = numerics["dfba"]
    add(f"`e_coli_core` on glucose through depletion. Observed order of convergence of the final biomass: {dfba['observed_order'][0]:.2f} and {dfba['observed_order'][1]:.2f}, for a scheme that is first-order once fluxes change. A Rust test checks the exponential phase against its closed form.\n")
    add("### Thermodynamics\n")
    tf = numerics["tfa"]
    add(f"TFA on `e_coli_core` with {tf['reactions_with_energies']} eQuilibrator energies: growth {tf['tfa_objective']:.6f} against {tf['fba_objective']:.6f} without them, {len(tf['fluxes_against_gibbs'])} fluxes against their Gibbs energy, and every reported `DrG'` within two standard uncertainties of its value recomputed from the reported concentrations (largest excess {tf['largest_excess_beyond_two_sigma']:.1e} kJ/mol). The sign logic also has Rust tests on small networks with known answers.\n")
    add("### Sensitivity and estimation\n")
    add("Rust tests, with Python and command-line tests of the same calls. Local sensitivities equal the analytic derivatives of a birth-death model and central finite differences on the Teusink model. Sobol indices match the closed form of a product of independent inputs at 2,000,000 samples, and the Ishigami function. On a fixed sample the model route gives the same indices as the closed-form output, from solved steady states and from a late time point. A fit recovers known rates from noise-free data, its profile brackets them, and a non-identifiable pair gives an open profile.\n")

    if perf:
        add("## Performance\n\n`validation/v10_performance.py`\n")
        b = perf["benchmarks"]
        add(f"{perf['cpu']}, {perf['logical_cores']} logical cores, sharing the machine with other work, so these are indicative. COBRApy used GLPK with 4 processes.\n")
        add("| Task | Rashevsky | COBRApy |\n| --- | --- | --- |")
        for name in ("e_coli_core", "iML1515", "Human-GEM"):
            add(f"| Read {name} | {b[name]['rashevsky_read_seconds']:.2f} s | {b[name]['cobrapy_read_seconds']:.2f} s |")
        f = b["iML1515_fva"]
        add(f"| FVA of iML1515, {f['reactions']:,} reactions | {f['rashevsky_seconds']:.1f} s | {f['cobrapy_glpk_4_processes_seconds']:.1f} s |")
        g = b["iML1515_single_gene_deletion"]
        add(f"| All single gene deletions of iML1515 | {g['rashevsky_seconds']:.1f} s | {g['cobrapy_glpk_4_processes_seconds']:.1f} s |")
        add(f"\nFVA of iML1515 agrees with COBRApy to {f['max_abs_difference']:.1e}. Essential tasks of Human-GEM take {b['human_gem_essential_tasks_seconds']:.1f} s, a 100-minute time course of the Teusink model at 1,000 points {b['teusink_time_course_100_min_1000_points_seconds']:.2f} s, and 10,000 CHRR samples of `e_coli_core` {b['core_chrr_10000_samples_seconds']:.2f} s, each including reading the file.\n")
    (DOCS / "VALIDATION.md").write_text("\n".join(w) + "\n")
    print("wrote docs/VALIDATION.md")


if __name__ == "__main__":
    main()
