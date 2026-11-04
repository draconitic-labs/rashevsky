#!/usr/bin/env python3
"""Builds the Rhea index of the Biochemistry site's map from its expanded map file.

Usage: build_scene_index.py map-expanded.json OUT.json
The result maps a Rhea identifier to the Reactome events of the map that carry it.
"""
import json
import sys

events = json.load(open(sys.argv[1]))["events"]
index: dict[str, list[str]] = {}
for event in events:
    for rhea in event.get("rhea_ids", []):
        index.setdefault(rhea.split(":")[-1], []).append(event["id"])
json.dump({k: sorted(set(v)) for k, v in sorted(index.items())}, open(sys.argv[2], "w"), indent=0, sort_keys=True)
print(len(index), "Rhea identifiers")
