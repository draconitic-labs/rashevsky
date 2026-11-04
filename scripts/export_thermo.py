#!/usr/bin/env python3
"""Exports standard transformed reaction energies for a BiGG-namespace model from eQuilibrator.

The energies come from the component contribution method (Noor et al. 2013; Beber et al. 2022) through
the `equilibrator-api` package, whose compound database is CC BY 4.0 (Zenodo 10.5281/zenodo.4128543).
Reactions with metabolites eQuilibrator cannot resolve, transport across a membrane, and pseudo
reactions such as biomass and boundary reactions are left out; the engine leaves those unconstrained.

Usage: export_thermo.py MODEL.xml[.gz] OUT.json [--ph 7.5] [--pmg 3] [--ionic 0.25] [--temperature 298.15]
"""
from __future__ import annotations

import argparse
import gzip
import json
import re


def read_species_and_reactions(path: str):
    text = gzip.open(path, "rt").read() if path.endswith(".gz") else open(path).read()
    species = dict(re.findall(r'<species[^>]*?id="M_([^"]+)"[^>]*?compartment="([^"]+)"', text))
    reactions = {}
    for block in re.findall(r"<reaction .*?</reaction>", text, flags=re.S):
        rid = re.search(r'id="R_([^"]+)"', block).group(1)
        stoich = {}
        for tag, sign in (("listOfReactants", -1.0), ("listOfProducts", 1.0)):
            section = re.search(rf"<{tag}>(.*?)</{tag}>", block, flags=re.S)
            if not section:
                continue
            for ref in re.findall(r"<speciesReference [^>]*>", section.group(1)):
                met = re.search(r'species="M_([^"]+)"', ref).group(1)
                coef = re.search(r'stoichiometry="([^"]+)"', ref)
                stoich[met] = stoich.get(met, 0.0) + sign * (float(coef.group(1)) if coef else 1.0)
        reactions[rid] = stoich
    return species, reactions


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("model")
    parser.add_argument("out")
    parser.add_argument("--ph", type=float, default=7.5)
    parser.add_argument("--pmg", type=float, default=3.0)
    parser.add_argument("--ionic", type=float, default=0.25)
    parser.add_argument("--temperature", type=float, default=298.15)
    args = parser.parse_args()

    from equilibrator_api import ComponentContribution, Q_, Reaction

    cc = ComponentContribution()
    cc.p_h = Q_(args.ph)
    cc.p_mg = Q_(args.pmg)
    cc.ionic_strength = Q_(f"{args.ionic}M")
    cc.temperature = Q_(f"{args.temperature}K")
    species, reactions = read_species_and_reactions(args.model)
    excluded = sorted(m for m in species if m.rsplit("_", 1)[0] in ("h", "h2o"))

    cache: dict[str, object] = {}

    def compound(met: str):
        base = met.rsplit("_", 1)[0]
        if base not in cache:
            cache[base] = cc.get_compound(f"bigg.metabolite:{base}")
        return cache[base]

    out = {}
    skipped = {}
    for rid, stoich in reactions.items():
        if len(stoich) < 2 or rid.startswith(("BIOMASS", "ATPM", "EX_")):
            skipped[rid] = "boundary or pseudo reaction"
            continue
        compartments = {species[m] for m in stoich}
        if len(compartments) > 1:
            skipped[rid] = "crosses a membrane"
            continue
        try:
            sparse = {}
            for met, coef in stoich.items():
                c = compound(met)
                if c is None:
                    raise KeyError(met)
                sparse[c] = sparse.get(c, 0.0) + coef
            reaction = Reaction(sparse)
            if not reaction.is_balanced():
                skipped[rid] = "unbalanced in eQuilibrator's protonation states"
                continue
            estimate = cc.standard_dg_prime(reaction)
            value = estimate.value.m_as("kJ/mol")
            sigma = estimate.error.m_as("kJ/mol")
        except Exception as exc:  # noqa: BLE001
            skipped[rid] = f"{type(exc).__name__}: {exc}"[:120]
            continue
        out[rid] = {"dg0": float(value), "sigma": float(sigma)}
    result = {
        "temperature": args.temperature,
        "reactions": out,
        "excluded_metabolites": excluded,
        "concentration_bounds": {},
        "source": (
            f"eQuilibrator component contribution, equilibrator-api, pH {args.ph}, pMg {args.pmg}, "
            f"I {args.ionic} M, {args.temperature} K; data CC BY 4.0, doi:10.5281/zenodo.4128543"
        ),
        "skipped": skipped,
    }
    with open(args.out, "w") as handle:
        json.dump(result, handle, indent=1, sort_keys=True)
    print(f"{len(out)} reactions with energies, {len(skipped)} skipped")


if __name__ == "__main__":
    main()
