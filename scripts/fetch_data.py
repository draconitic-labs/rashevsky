#!/usr/bin/env python3
"""Downloads the reference data the validation scripts compare against, and checks its digests.

None of this data is in the repository: each file keeps its own license, listed in docs/DATA.md.
Usage: fetch_data.py [--dir PATH]. The default folder is $RASHEVSKY_DATA, or `.data` in the repository.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import pathlib
import urllib.request
import zipfile

ROOT = pathlib.Path(__file__).resolve().parents[1]
MANIFEST = ROOT / "validation" / "data-manifest.json"


def digest(path: pathlib.Path) -> str:
    h = hashlib.sha256()
    with open(path, "rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    default = os.environ.get("RASHEVSKY_DATA", str(ROOT / ".data"))
    parser.add_argument("--dir", default=default)
    args = parser.parse_args()
    folder = pathlib.Path(args.dir)
    folder.mkdir(parents=True, exist_ok=True)
    for item in json.loads(MANIFEST.read_text()):
        path = folder / item["file"]
        if not path.exists():
            print(f"downloading {item['file']}", flush=True)
            request = urllib.request.Request(item["url"], headers={"User-Agent": "rashevsky-validation"})
            with urllib.request.urlopen(request, timeout=600) as reply, open(path, "wb") as out:
                out.write(reply.read())
        found = digest(path)
        if found != item["sha256"]:
            # A changed upstream file changes the results, so it must not pass silently.
            print(f"DIGEST DIFFERS for {item['file']}: found {found}, pinned {item['sha256']}")
        else:
            print(f"ok {item['file']}")
        if item.get("unzip"):
            target = folder / item["unzip"]
            if not target.exists():
                with zipfile.ZipFile(path) as archive:
                    archive.extractall(target)


if __name__ == "__main__":
    main()
