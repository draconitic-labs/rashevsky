#!/usr/bin/env python3
"""Writes THIRD-PARTY-NOTICES.md: the license text of everything linked into a release.

MIT, BSD and Apache terms ask for their copyright lines and text in every copy, so the list of names is
not enough. Texts come from the crates' own files, as published. `--check` only compares, for CI.
"""
from __future__ import annotations

import json
import os
import pathlib
import re
import subprocess
import sys

ROOT = pathlib.Path(__file__).resolve().parents[1]
OUT = ROOT / "THIRD-PARTY-NOTICES.md"
WORKSPACES = [
    (ROOT / "Cargo.toml", "x86_64-unknown-linux-gnu"),
    (ROOT / "Cargo.toml", "wasm32-unknown-unknown"),
    (ROOT / "desktop" / "src-tauri" / "Cargo.toml", "x86_64-unknown-linux-gnu"),
]
TEXT_NAME = re.compile(r"^(licen[cs]e|copying|copyright|notice|unlicense|mit|apache)", re.I)
SKIP_DIRS = {"tests", "test", "examples", "benches", "docs", "doc", "target", ".git", ".github"}
# Which side of a dual license this project takes, first found wins.
PREFERRED = ["MIT", "Apache-2.0", "BSD-3-Clause", "BSD-2-Clause", "ISC", "Zlib", "Unlicense", "0BSD", "MPL-2.0"]
HAS_THE_OTHER_SIDE = {
    "MIT": re.compile(r"apache", re.I),
    "Apache-2.0": re.compile(r"(^|[-_.])mit([-_.]|$)", re.I),
}
STOCK_MIT = """Permission is hereby granted, free of charge, to any person obtaining a copy of this software and
associated documentation files (the "Software"), to deal in the Software without restriction, including
without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the
following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial
portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED
TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL
THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF
CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
DEALINGS IN THE SOFTWARE."""


def linked_packages(manifest: pathlib.Path, platform: str) -> dict[str, dict]:
    """Packages reachable from the workspace through normal dependencies, which are the ones linked."""
    meta = json.loads(subprocess.run(
        ["cargo", "metadata", "--format-version", "1", "--locked", "--manifest-path", str(manifest), "--filter-platform", platform],
        check=True, capture_output=True, text=True,
    ).stdout)
    nodes = {n["id"]: n for n in meta["resolve"]["nodes"]}
    packages = {p["id"]: p for p in meta["packages"]}
    seen: set[str] = set()
    todo = list(meta["workspace_members"])
    while todo:
        at = todo.pop()
        if at in seen:
            continue
        seen.add(at)
        for dep in nodes[at]["deps"]:
            if any(kind["kind"] is None for kind in dep["dep_kinds"]):
                todo.append(dep["pkg"])
    members = set(meta["workspace_members"])
    # A package with no source is a path dependency: this project's own engine, under its own license.
    return {i: packages[i] for i in seen - members if packages[i].get("source")}


def chosen_license(expression: str | None) -> str:
    if not expression:
        return ""
    names = re.split(r"\s+(?:OR|AND)\s+|/", expression.replace("(", " ").replace(")", " ").strip())
    names = [n.strip() for n in names if n.strip()]
    if " AND " in expression:
        return ""
    for name in PREFERRED:
        if name in names:
            return name
    return names[0] if names else ""


def license_files(package: dict, license_name: str) -> list[pathlib.Path]:
    root = pathlib.Path(package["manifest_path"]).parent
    found = []
    other = HAS_THE_OTHER_SIDE.get(license_name)
    for folder, dirs, files in os.walk(root):
        depth = len(pathlib.Path(folder).relative_to(root).parts)
        dirs[:] = sorted(d for d in dirs if d.lower() not in SKIP_DIRS and depth < 4)
        for name in sorted(files):
            if not TEXT_NAME.match(name) or name.lower().endswith((".rs", ".py", ".h", ".c", ".cpp", ".toml", ".json", ".yml")):
                continue
            if depth == 0 and other and other.search(name):
                continue
            found.append(pathlib.Path(folder) / name)
    return found


def main() -> int:
    packages: dict[str, dict] = {}
    for manifest, platform in WORKSPACES:
        packages.update(linked_packages(manifest, platform))
    rows = []
    groups: dict[str, list[str]] = {}
    without = []
    for package in sorted(packages.values(), key=lambda p: (p["name"], p["version"])):
        label = f"{package['name']} {package['version']}"
        rows.append(f"| {package['name']} | {package['version']} | {package.get('license') or 'see its text below'} |")
        name = chosen_license(package.get("license"))
        files = license_files(package, name)
        if not files:
            without.append(label)
            authors = ", ".join(package.get("authors") or []) or f"the authors of {package['name']}"
            if name == "MIT":
                groups.setdefault(f"Copyright (c) {authors}\n\n{STOCK_MIT}", []).append(label)
            else:
                groups.setdefault(f"Licensed under {package.get('license')}, by {authors}. The published crate carries no license file; the text is at https://spdx.org/licenses/{name}.html and the source at {package.get('repository') or 'https://crates.io/crates/' + package['name']}", []).append(label)
            continue
        for path in files:
            text = path.read_text(errors="replace").strip()
            relative = path.relative_to(pathlib.Path(package["manifest_path"]).parent)
            key = text
            groups.setdefault(key, []).append(label if len(relative.parts) == 1 else f"{label} ({relative})")
    lines = [
        "# Third-party software",
        "",
        "Rashevsky links the HiGHS solver, built from source by the `highs-sys` crate, and the Rust crates below.",
        "Fonts and data are listed in `docs/DATA.md`. `scripts/make_third_party_notices.py` writes this file from `cargo metadata` and the texts the crates publish.",
        "Where a crate offers a choice of licenses, this project takes the first of MIT and Apache-2.0 that it offers.",
        "",
        "| Crate | Version | License |",
        "| --- | --- | --- |",
        *rows,
        "",
        "## License texts",
        "",
    ]
    for text, labels in sorted(groups.items(), key=lambda item: item[1][0]):
        lines += [f"### {', '.join(labels)}", "", "````text", text, "````", ""]
    body = "\n".join(lines)
    if "--check" in sys.argv:
        if not OUT.exists() or OUT.read_text() != body:
            print("THIRD-PARTY-NOTICES.md is out of date: run scripts/make_third_party_notices.py", file=sys.stderr)
            return 1
        print("up to date")
        return 0
    OUT.write_text(body)
    print(f"{len(rows)} crates, {len(groups)} distinct texts, {OUT.stat().st_size:,} bytes; no file in: {without}")
    return 0


if __name__ == "__main__":
    sys.exit(main())
