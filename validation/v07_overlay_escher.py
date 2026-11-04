#!/usr/bin/env python3
"""Flux overlay on the Biochemistry site's map, and Escher map output.

Overlay: pFBA on Human-GEM with its own bounds is projected onto the site's scene. Pass requires every
event identifier to exist in the scene, event fluxes to equal the sums of their reactions in an
independent recomputation from the model's annotations, and at least 300 events reached by an exact
identifier. Escher: the map validates against Escher's published 1-0-0 JSON schema, and reaction data
covers every reaction.
"""
from __future__ import annotations

import gzip
import json
import re
import sys
import urllib.request

import jsonschema

from common import DATA, ICH360, REPO, engine, write_result

SCENE = REPO / "data" / "biochemistry" / "scene.json"
RHEA = REPO / "data" / "biochemistry" / "rhea-events.json"
MODEL = DATA / "Human-GEM.xml"
RESULTS = REPO / "validation" / "results"
SCHEMA_URL = "https://raw.githubusercontent.com/opencobra/escher/master/py/escher/static/jsonschema/1-0-0"


def main() -> int:
    svg = RESULTS / "human_gem_overlay.svg"
    overlay = engine("overlay", MODEL, "--scene", SCENE, "--rhea", RHEA, "--svg", svg)
    scene = json.load(open(SCENE))
    known = {e for r in scene["reactions"] for e in r["e"]} | {e for r in scene["processes"] for e in r["e"]}
    unknown = [e for e in overlay["events"] if e not in known]

    # Recompute exact-match event totals straight from the SBML text with a separate parser.
    fluxes = engine("pfba", MODEL)["fluxes"]
    rhea = json.load(open(RHEA))
    text = open(MODEL, errors="ignore").read()
    expected: dict[str, float] = {}
    for block in re.findall(r"<reaction .*?</reaction>", text, flags=re.S):
        rid = re.search(r'id="R_([^"]+)"', block).group(1)
        hits = set(e for e in re.findall(r"identifiers\.org/reactome/(R-HSA-\d+)", block) if e in known)
        for code in re.findall(r"identifiers\.org/rhea/(\d+)", block):
            hits |= {e for e in rhea.get(code, []) if e in known}
        for event in hits:
            expected[event] = expected.get(event, 0.0) + abs(fluxes[rid])
    # pFBA in two runs of the same deterministic solve must agree; compare event by event.
    worst = 0.0
    for event, value in overlay["events"].items():
        if value["match"] != "ec":
            worst = max(worst, abs(value["flux"] - expected.get(event, 0.0)))
    coverage = overlay["coverage"]

    # Escher.
    map_path = RESULTS / "ich360_escher_map.json"
    data_path = RESULTS / "ich360_reaction_data.json"
    engine("escher", ICH360, "--map", map_path, "--data", data_path)
    doc = json.load(open(map_path))
    data = json.load(open(data_path))
    try:
        schema = json.load(urllib.request.urlopen(SCHEMA_URL, timeout=30))
        jsonschema.validate(doc, schema)
        schema_ok = True
        schema_note = SCHEMA_URL
    except jsonschema.ValidationError as exc:
        schema_ok = False
        schema_note = str(exc)[:300]
    n_rxn = engine("info", ICH360)["reactions"]
    passed = not unknown and worst < 1e-6 and coverage["events_exact"] >= 300 and schema_ok and len(data) == n_rxn
    write_result("v07_overlay_escher", {
        "scene_release": overlay["scene_release"],
        "objective_human_gem": overlay["objective"],
        "coverage": coverage,
        "unknown_event_ids": unknown[:10],
        "max_abs_difference_exact_events_vs_independent_recompute": worst,
        "overlay_svg": str(svg.relative_to(REPO)),
        "escher_schema": schema_note,
        "escher_schema_valid": schema_ok,
        "escher_reactions_in_data": len(data),
        "pass": passed,
    })
    print(coverage, "objective", overlay["objective"], "worst", worst, "schema", schema_ok, schema_note[:200])
    return 0 if passed else 1


if __name__ == "__main__":
    sys.exit(main())
