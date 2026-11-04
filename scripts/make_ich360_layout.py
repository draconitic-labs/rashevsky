#!/usr/bin/env python3
"""Writes models/ich360/layout.json from the authors' Escher map of iCH360.

Positions, routes and pathway names come from the map in the iCH360 repository (MIT), which
scripts/fetch_data.py downloads at a pinned commit. Each reaction's substrates and products come from
the model's own stoichiometry, so the drawing cannot disagree with the model about direction. As on the
Escher map, cofactors stay off the pathways: only metabolites marked primary are drawn.
"""
from __future__ import annotations

import gzip
import json
import math
import pathlib
import re
import sys
from collections import deque

ROOT = pathlib.Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "validation"))
from common import DATA, ICH360  # noqa: E402
from v02_helpers import reaction_table  # noqa: E402

# Map units to layout units: the renderer's line and text sizes suit a world a few thousand wide.
SCALE = 0.4
MARGIN = 60
# A route further than this from a straight line keeps its curve.
CURVE_TOLERANCE = 5.0
EXCHANGE_LENGTH = 70
LONGEST_NAME = 24
REGIONS = {
    "Glycolysis": "Glycolysis",
    "TCA cycle": "Citric acid cycle",
    "Pentose phosphate pathway": "Pentose phosphate",
    "glucose uptake": "Glucose uptake",
    "Pyrimidines de novo biosynthesis": "Pyrimidines",
    "Purines de novo biosynthesis": "Purines",
    "Fatty acids: initiation": "Fatty acid initiation",
    "Fatty acids: elongation": "Fatty acid elongation",
    "Biosynthesis of unsaturated fatty acids": "Unsaturated fatty acids",
    "Biosynthesis of aromatic amino acids": "Aromatic amino acids",
    "Histidine biosynthesis": "Histidine",
    "Sulfate uptake and metabolism": "Sulfate",
    "C1 Pool": "C1 pool",
    "Electron transport chain": "Electron transport",
    "Energy and maintenance": "Energy and maintenance",
    "reduction of redox proteins": "Redox proteins",
    "uptake of water, phosphorus, and nitrogen": "Water, phosphate and nitrogen uptake",
}
# Where a name would land on a line at the map's own position: x, y and which end of the text is there.
REGION_AT = {
    "Glucose uptake": (901, 136, "end"),
    "Glycolysis": (901, 975, "end"),
    "Pentose phosphate": (963, 364, "start"),
    "Sulfate": (123, 214, "start"),
    "Energy and maintenance": (2583, 557, "start"),
    "Citric acid cycle": (1787, 1255, "middle"),
    "C1 pool": (239, 1172, "start"),
    "Fatty acid initiation": (773, 2032, "end"),
    "Fatty acid elongation": (1863, 2351, "start"),
    "Water, phosphate and nitrogen uptake": (212, 2783, "start"),
}
# Common names, where the model's would not fit on the map.
SHORT = {
    "g6p_c": "G6P", "f6p_c": "F6P", "fdp_c": "FBP", "dhap_c": "DHAP", "g3p_c": "GAP", "13dpg_c": "1,3-BPG",
    "3pg_c": "3PG", "2pg_c": "2PG", "pep_c": "PEP", "pyr_c": "Pyruvate", "6pgl_c": "6PGL", "6pgc_c": "6PG",
    "ru5p__D_c": "Ru5P", "r5p_c": "R5P", "xu5p__D_c": "Xu5P", "s7p_c": "S7P", "e4p_c": "E4P", "cit_c": "Citrate",
    "acon_C_c": "cis-Aconitate", "icit_c": "Isocitrate", "akg_c": "2-Oxoglutarate", "succoa_c": "Succinyl-CoA",
    "succ_c": "Succinate", "fum_c": "Fumarate", "mal__L_c": "Malate", "oaa_c": "Oxaloacetate", "glx_c": "Glyoxylate",
    "accoa_c": "Acetyl-CoA",
}


def species_names(path: pathlib.Path) -> dict[str, str]:
    text = gzip.open(path, "rt").read()
    names = {}
    for tag in re.findall(r"<species [^>]*>", text):
        names[re.search(r'id="M_([^"]+)"', tag).group(1)] = (re.search(r'name="([^"]*)"', tag) or [None, ""])[1]
    return names


def label_for(met: str, names: dict[str, str]) -> str:
    if met in SHORT:
        return SHORT[met]
    name = names.get(met, "")
    return name if name and len(name) <= LONGEST_NAME else met.rsplit("_", 1)[0]


def route(reaction: dict, nodes: dict, start: str, goal: str) -> list[tuple[float, float]]:
    """Points along the drawn path between two nodes of a reaction, curves included."""
    links: dict[str, list] = {}
    for seg in reaction["segments"].values():
        a, b = seg["from_node_id"], seg["to_node_id"]
        links.setdefault(a, []).append((b, seg["b1"], seg["b2"]))
        links.setdefault(b, []).append((a, seg["b2"], seg["b1"]))
    came: dict[str, tuple] = {start: None}
    queue = deque([start])
    while queue and goal not in came:
        at = queue.popleft()
        for to, b1, b2 in links.get(at, []):
            if to not in came:
                came[to] = (at, b1, b2)
                queue.append(to)
    if goal not in came:
        return []
    chain = []
    at = goal
    while came[at] is not None:
        chain.append((came[at][0], at, came[at][1], came[at][2]))
        at = came[at][0]
    points: list[tuple[float, float]] = []
    for a, b, b1, b2 in reversed(chain):
        p0, p3 = nodes[a], nodes[b]
        c1 = b1 or p0
        c2 = b2 or p3
        for i in range(12):
            t = i / 12
            u = 1 - t
            points.append((
                u**3 * p0["x"] + 3 * u * u * t * c1["x"] + 3 * u * t * t * c2["x"] + t**3 * p3["x"],
                u**3 * p0["y"] + 3 * u * u * t * c1["y"] + 3 * u * t * t * c2["y"] + t**3 * p3["y"],
            ))
    points.append((nodes[goal]["x"], nodes[goal]["y"]))
    return points


def fit_control(points: list[tuple[float, float]]) -> tuple[tuple[float, float], float]:
    """The control point of the quadratic through the two ends closest to the route, and how far the
    route strays from a straight line."""
    lengths = [0.0]
    for p, q in zip(points, points[1:]):
        lengths.append(lengths[-1] + math.dist(p, q))
    total = lengths[-1] or 1.0
    a, d = points[0], points[-1]
    numerator = [0.0, 0.0]
    denominator = 0.0
    stray = 0.0
    chord = math.dist(a, d) or 1.0
    for p, s in zip(points, lengths):
        t = s / total
        w = 2 * t * (1 - t)
        denominator += w * w
        for axis in (0, 1):
            numerator[axis] += w * (p[axis] - (1 - t) ** 2 * a[axis] - t * t * d[axis])
        stray = max(stray, abs((d[0] - a[0]) * (a[1] - p[1]) - (a[0] - p[0]) * (d[1] - a[1])) / chord)
    if denominator == 0:
        return ((a[0] + d[0]) / 2, (a[1] + d[1]) / 2), 0.0
    return (numerator[0] / denominator, numerator[1] / denominator), stray


def label_side(x: float, y: float, towards: list[tuple[float, float]]) -> str:
    """The side of a node farthest from the lines leaving it; above unless another is clearly better."""
    sides = {"t": -math.pi / 2, "r": 0.0, "b": math.pi / 2, "l": math.pi}
    angles = [math.atan2(b - y, a - x) for a, b in towards if (a, b) != (x, y)]
    if not angles:
        return "t"

    def clearance(side: float) -> float:
        return min(abs((side - angle + math.pi) % (2 * math.pi) - math.pi) for angle in angles)

    best = max(("r", "l", "t", "b"), key=lambda s: clearance(sides[s]))
    return "t" if clearance(sides["t"]) >= clearance(sides[best]) - 0.2 else best


def main() -> None:
    source = json.load(open(DATA / "iCH360_full_map.json"))[1]
    table = reaction_table(ICH360)
    names = species_names(ICH360)
    nodes = source["nodes"]
    drawn = {r["bigg_id"]: r for r in source["reactions"].values() if r["bigg_id"] in table}

    primary = {}
    for reaction in drawn.values():
        for seg in reaction["segments"].values():
            for end in (seg["from_node_id"], seg["to_node_id"]):
                node = nodes[end]
                if node["node_type"] == "metabolite" and node["node_is_primary"]:
                    primary[end] = node
    x0 = min(n["x"] for n in primary.values())
    y0 = min(n["y"] for n in primary.values())
    fx = lambda x: round((x - x0) * SCALE + MARGIN, 1)  # noqa: E731
    fy = lambda y: round((y - y0) * SCALE + MARGIN, 1)  # noqa: E731

    keys: dict[str, str] = {}
    used: dict[str, int] = {}
    reactions: dict[str, dict] = {}
    towards: dict[str, list] = {}
    for rid in sorted(drawn):
        reaction = drawn[rid]
        stoich = table[rid]
        mid = next(
            end
            for seg in reaction["segments"].values()
            for end in (seg["from_node_id"], seg["to_node_id"])
            if nodes[end]["node_type"] == "midmarker"
        )
        ends = {
            end
            for seg in reaction["segments"].values()
            for end in (seg["from_node_id"], seg["to_node_id"])
            if nodes[end]["node_type"] == "metabolite" and nodes[end]["bigg_id"] in stoich
        }
        chosen = sorted(e for e in ends if nodes[e]["node_is_primary"])
        if not chosen:
            raise SystemExit(f"{rid} has no primary metabolite on the map")
        entry: dict = {"s": [], "p": [], "m": [fx(nodes[mid]["x"]), fy(nodes[mid]["y"])]}
        curves = {}
        for end in chosen:
            met = nodes[end]["bigg_id"]
            if end not in keys:
                used[met] = used.get(met, 0) + 1
                keys[end] = met if used[met] == 1 else f"{met}@{used[met]}"
            key = keys[end]
            (entry["s"] if stoich[met] < 0 else entry["p"]).append(key)
            towards.setdefault(key, []).append((fx(nodes[mid]["x"]), fy(nodes[mid]["y"])))
            points = route(reaction, nodes, end, mid)
            if len(points) > 2:
                control, stray = fit_control([((x - x0) * SCALE + MARGIN, (y - y0) * SCALE + MARGIN) for x, y in points])
                if stray > CURVE_TOLERANCE:
                    curves[key] = [round(control[0], 1), round(control[1], 1)]
        if curves:
            entry["c"] = curves
        if len(stoich) == 1:
            entry["boundary"] = True
            n = nodes[chosen[0]]
            dx, dy = nodes[mid]["x"] - n["x"], nodes[mid]["y"] - n["y"]
            size = math.hypot(dx, dy) or 1.0
            entry["end"] = [round(fx(n["x"]) + EXCHANGE_LENGTH * dx / size, 1), round(fy(n["y"]) + EXCHANGE_LENGTH * dy / size, 1)]
            del entry["m"]
            entry.pop("c", None)
        reactions[rid] = entry

    out_nodes = {}
    for end, key in keys.items():
        node = nodes[end]
        met = node["bigg_id"]
        x, y = fx(node["x"]), fy(node["y"])
        item = {"x": x, "y": y, "label": label_for(met, names), "external": met.endswith("_e")}
        if key != met:
            item["metabolite"] = met
        side = label_side(x, y, towards.get(key, []))
        if side != "t":
            item["side"] = side
        out_nodes[key] = item

    width = max(n["x"] for n in out_nodes.values()) + MARGIN
    height = max(n["y"] for n in out_nodes.values()) + MARGIN
    regions = []
    for text in source["text_labels"].values():
        if text["text"] in REGIONS:
            label = REGIONS[text["text"]]
            x, y, anchor = REGION_AT.get(label, (fx(text["x"]), fy(text["y"]), "start"))
            regions.append({"label": label, "x": x, "y": y, "anchor": anchor})
    missing = sorted(set(table) - set(reactions) - {"Biomass"})
    if missing:
        raise SystemExit(f"not on the map: {missing}")
    reactions["Biomass"] = {"biomass": True, "s": [], "p": [], "m": BIOMASS_AT}
    layout = {
        "world": [round(width), round(height)],
        "nodes": out_nodes,
        "reactions": reactions,
        "regions": regions,
        "source": "the authors' Escher map at commit 08dd438, scripts/make_ich360_layout.py",
    }
    out = ROOT / "models" / "ich360" / "layout.json"
    out.write_text(json.dumps(layout, indent=0, sort_keys=True) + "\n")
    curved = sum(1 for r in reactions.values() if "c" in r)
    print(f"{len(reactions)} reactions, {curved} curved, {len(out_nodes)} nodes, world {layout['world']}, {out.stat().st_size} bytes")


# An empty part of the map, for the ring that shows growth.
BIOMASS_AT = [2559, 1111]

if __name__ == "__main__":
    main()
