"""Keep bundled scientific references pinned and usable without network access."""
from __future__ import annotations

import hashlib
import importlib.util
import json
import pathlib
import tempfile
import unittest
import zipfile
from unittest.mock import patch


ROOT = pathlib.Path(__file__).resolve().parents[2]
SPEC = importlib.util.spec_from_file_location("fetch_data", ROOT / "scripts/fetch_data.py")
fetch_data = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fetch_data)


class ReferenceDataTests(unittest.TestCase):
    def test_bundled_models_match_the_original_pins(self):
        items = json.loads(fetch_data.MANIFEST.read_text())
        bundled = [item for item in items if item.get("bundled")]
        self.assertEqual(len(bundled), 5)
        for item in bundled:
            with self.subTest(file=item["file"]):
                data = (ROOT / item["bundled"]).read_bytes()
                self.assertEqual(hashlib.sha256(data).hexdigest(), item["sha256"])

    def test_bundled_references_work_with_network_unavailable(self):
        items = [item for item in json.loads(fetch_data.MANIFEST.read_text()) if item.get("bundled")]
        with tempfile.TemporaryDirectory() as folder:
            folder = pathlib.Path(folder)
            manifest = folder / "manifest.json"
            manifest.write_text(json.dumps(items))
            cache = folder / "cache"
            with (
                patch.object(fetch_data, "MANIFEST", manifest),
                patch("sys.argv", ["fetch_data.py", "--dir", str(cache)]),
                patch.object(
                    fetch_data.urllib.request,
                    "urlopen",
                    side_effect=AssertionError("network must not be needed"),
                ),
            ):
                fetch_data.main()
            for item in items:
                self.assertEqual(
                    (cache / item["file"]).read_bytes(),
                    (ROOT / item["bundled"]).read_bytes(),
                )

    def test_a_digest_mismatch_never_extracts_an_archive(self):
        with tempfile.TemporaryDirectory() as folder:
            folder = pathlib.Path(folder)
            cache = folder / "cache"
            cache.mkdir()
            with zipfile.ZipFile(cache / "reference.zip", "w") as archive:
                archive.writestr("model.xml", "unverified model")
            manifest = folder / "manifest.json"
            manifest.write_text(json.dumps([
                {"file": "reference.zip", "sha256": "0" * 64, "unzip": "extracted"}
            ]))
            with (
                patch.object(fetch_data, "MANIFEST", manifest),
                patch("sys.argv", ["fetch_data.py", "--dir", str(cache)]),
                self.assertRaises(SystemExit),
            ):
                fetch_data.main()
            self.assertFalse((cache / "extracted").exists())


if __name__ == "__main__":
    unittest.main()
