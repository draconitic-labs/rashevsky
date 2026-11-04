#!/usr/bin/env python3
"""iML1515 single-gene essentiality against the Keio screens of Monk et al. 2017.

Reference: Monk et al., Nat. Biotechnol. 35 (2017) 904, Supplementary Note 4 and Supplementary
Data 11, Table 1: 23,617 gene/condition phenotypes on 16 carbon sources, reported accuracy 93.4 %.
Protocol from the paper's Supplementary Methods: oxygen uptake bounded at 20 and the carbon source
at 10 mmol/gDW/h; FHL, and the reactions of genes absent from BW25113 (ARAI, RBK_L1, RMPA, LYXI,
RMI, RMK, LACZ), closed; a gene is computationally essential below 0.05 1/h.

Pass: accuracy within 0.5 percentage points of the published 93.4 %.
"""
from __future__ import annotations

import collections
import statistics
import sys
import time

from common import DATA, engine, monk_table, write_result

MODEL = DATA / "iML1515.xml.gz"
PUBLISHED = 93.4
CARBON = {
    "Acetate": "EX_ac_e",
    "Galactose": "EX_gal_e",
    "Gluconate": "EX_glcn_e",
    "Glucose": "EX_glc__D_e",
    "Glycerol": "EX_glyc_e",
    "L-Alanine": "EX_ala__L_e",
    "Lactate": "EX_lac__L_e",
    "Maltose": "EX_malt_e",
    "Mannitol": "EX_mnl_e",
    "N-acetyl-glucosamine": "EX_acgam_e",
    "Oxoglutarate": "EX_akg_e",
    "Pyruvate": "EX_pyr_e",
    "Ribose": "EX_rib__D_e",
    "Sorbitol": "EX_sbt__D_e",
    "Succinate": "EX_succ_e",
    "Xylose": "EX_xyl__D_e",
}
CLOSED = ["FHL", "ARAI", "RBK_L1", "RMPA", "LYXI", "RMI", "RMK", "LACZ"]
THRESHOLD = 0.05


def condition_bounds(condition: str) -> list[str]:
    args = ["--bound", "EX_o2_e=-20:1000"]
    if condition != "Glucose":
        args += ["--bound", "EX_glc__D_e=0:1000"]
    args += ["--bound", f"{CARBON[condition]}=-10:1000"]
    for rid in CLOSED:
        args += ["--bound", f"{rid}=0:0"]
    return args


def main() -> int:
    rows = monk_table()
    by_condition = collections.defaultdict(list)
    for r in rows:
        by_condition[r["Condition"]].append(r)

    totals = collections.Counter()
    per_condition = {}
    wild_type = {}
    growth_diffs = []
    disagreements = []
    missing_genes = set()
    started = time.time()
    for condition, items in sorted(by_condition.items()):
        bounds = condition_bounds(condition)
        wt = engine("fba", MODEL, *bounds)["objective"]
        deletions = engine("deletion", MODEL, "--kind", "gene", *bounds)
        growth = {d["ids"][0]: d["growth"] for d in deletions}
        counts = collections.Counter()
        for r in items:
            gene = r["Gene"]
            if gene not in growth:
                missing_genes.add(gene)
                continue
            exp_growth = float(r["Measured_growth"]) >= float(r["Growth_cutoff"])
            mod_growth = growth[gene] >= THRESHOLD
            label = ("T" if exp_growth == mod_growth else "F") + ("P" if mod_growth else "N")
            counts[label] += 1
            growth_diffs.append(abs(growth[gene] - float(r["Computational_growth"])))
            if label != r["Prediction"]:
                disagreements.append({
                    "gene": gene, "condition": condition, "published": r["Prediction"], "rashevsky": label,
                    "published_growth": float(r["Computational_growth"]), "rashevsky_growth": growth[gene],
                })
        wild_type[condition] = wt
        totals.update(counts)
        per_condition[condition] = dict(counts)
        print(f"  {condition:22s} wt {wt:.5f} {dict(counts)}", flush=True)
    elapsed = time.time() - started

    # The same deletions in COBRApy on the same SBML file, to separate engine from model-file effects.
    cobra_check = {}
    import cobra
    from cobra.flux_analysis import single_gene_deletion
    for condition in ("Glucose", "Xylose"):
        model = cobra.io.read_sbml_model(str(MODEL))
        for arg in zip(condition_bounds(condition)[1::2]):
            rid, rng = arg[0].split("=")
            lo, hi = rng.split(":")
            model.reactions.get_by_id(rid).bounds = (float(lo), float(hi))
        t0 = time.time()
        ref = single_gene_deletion(model, processes=4)
        ref_time = time.time() - t0
        ours = {d["ids"][0]: d["growth"] for d in engine("deletion", MODEL, "--kind", "gene", *condition_bounds(condition))}
        worst, worst_gene = 0.0, None
        for gene, value in ours.items():
            rv = ref.knockout[gene]["growth"].iloc[0]
            rv = 0.0 if rv != rv else float(rv)
            if abs(rv - value) > worst:
                worst, worst_gene = abs(rv - value), gene
        cobra_check[condition] = {
            "genes": len(ours),
            "max_abs_difference_vs_cobrapy": worst,
            "worst_gene": worst_gene,
            "cobrapy_b3564": float(ref.knockout["b3564"]["growth"].iloc[0]),
            "rashevsky_b3564": ours["b3564"],
            "cobrapy_seconds_4_processes": ref_time,
        }
        print(f"  cobrapy {condition}: {cobra_check[condition]}", flush=True)

    n = sum(totals.values())
    accuracy = 100.0 * (totals["TP"] + totals["TN"]) / n
    tp, tn, fp, fn = totals["TP"], totals["TN"], totals["FP"], totals["FN"]
    mcc = (tp * tn - fp * fn) / ((tp + fp) * (tp + fn) * (tn + fp) * (tn + fn)) ** 0.5
    passed = abs(accuracy - PUBLISHED) <= 0.5 and all(
        c["max_abs_difference_vs_cobrapy"] < 1e-6 for c in cobra_check.values()
    )
    write_result("v02_iml1515_keio", {
        "reference": "Monk et al. 2017, Supplementary Data 11, Table 1",
        "published_accuracy_percent": PUBLISHED,
        "published_counts": {"TP": 19233, "FN": 295, "FP": 1245, "TN": 2844},
        "accuracy_percent": accuracy,
        "phenotypes": n,
        "counts": dict(totals),
        "mcc": mcc,
        "per_condition": per_condition,
        "wild_type_growth": wild_type,
        "rows_classified_differently_from_published": len(disagreements),
        "classification_disagreements_sample": disagreements[:40],
        "growth_abs_diff_vs_published": {
            "median": statistics.median(growth_diffs),
            "p99": sorted(growth_diffs)[int(0.99 * len(growth_diffs))],
            "max": max(growth_diffs),
        },
        "genes_in_table_not_in_model": sorted(missing_genes),
        "elapsed_seconds": elapsed,
        "cobrapy_same_file": cobra_check,
        "pass": passed,
    })
    print(f"accuracy {accuracy:.3f} % over {n} phenotypes; counts {dict(totals)}; MCC {mcc:.4f}; "
          f"{len(disagreements)} rows differ from the published call; {elapsed:.1f} s")
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
