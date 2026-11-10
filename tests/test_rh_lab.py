"""No third-party dependencies: run with python -m unittest discover -s tests."""

import json
import tempfile
import unittest
from fractions import Fraction
from pathlib import Path

from rh_lab.catalog import ROOT, load, validate
from rh_lab.controls import on_line, zero_residual
from rh_lab.runner import run


class LabTests(unittest.TestCase):
    def test_offline_quartet_exact(self):
        self.assertEqual(zero_residual(), (0, 0))
        a, b = Fraction(1, 4), Fraction(2)
        for t in [Fraction(i, 10) for i in range(-250, 251)]:
            self.assertGreaterEqual(on_line(t, a, b), 4 * a * a * b * b)

    def test_traversal_rejected(self):
        good = {"id": "bad", "path": "../../outside.py", "method": "exact_rational",
                "evidence": "numerical", "claim": "bad", "assumptions": ["x"],
                "falsifiers": ["x"], "limitations": ["x"],
                "timeout_s": 2, "default_run": True}
        with self.assertRaisesRegex(ValueError, "escapes repository"):
            validate(good, ROOT)

    def test_replay_manifest(self):
        with tempfile.TemporaryDirectory() as temp:
            home = Path(temp)
            folder = home / "runs"
            item = {"id": "control", "path": "rh_lab/controls.py", "method": "exact_rational",
                    "evidence": "paper_proved", "claim": "fixture", "assumptions": ["x"],
                    "falsifiers": ["x"], "limitations": ["x"],
                    "timeout_s": 15, "default_run": True}
            manifest = home / "catalog.json"
            manifest.write_text(json.dumps({"version": 1, "experiments": [item]}), encoding="utf-8")
            record = run("control", root=ROOT, catalog=manifest, output_dir=folder, replay=True)
            self.assertTrue(record["passed"])
            self.assertTrue(record["deterministic"])
            self.assertEqual(len(list(folder.glob("*.json"))), 1)
            self.assertEqual(json.loads(Path(record["artifact_path"]).read_text())["code_sha256"], record["code_sha256"])

    def test_declared_stdin_pipeline(self):
        with tempfile.TemporaryDirectory() as temp:
            home = Path(temp)
            base = {"path": "rh_lab/controls.py", "method": "exact_rational",
                    "evidence": "numerical", "claim": "test", "assumptions": ["test"],
                    "falsifiers": ["test"], "limitations": ["test"],
                    "timeout_s": 15, "default_run": False}
            manifest = home / "catalog.json"
            manifest.write_text(json.dumps({"version": 1, "experiments": [
                dict(base, id="producer"), dict(base, id="consumer", stdin_from="producer")]}),
                                encoding="utf-8")
            result = run("consumer", root=ROOT, catalog=manifest,
                         output_dir=home / "runs", replay=True)
            self.assertTrue(result["passed"])
            self.assertIsNotNone(result["first"]["input_sha256"])
            self.assertIsNotNone(result["stdin_source_sha256"])

    def test_existing_registry_when_checkout_available(self):
        if not (ROOT / "rh_lab" / "controls.py").exists():
            self.skipTest("full repo checkout required for catalog validation")
        registry = load()
        self.assertIn("fake_quartet_control", registry)
        self.assertEqual(registry["fake_quartet_control"]["path"], "rh_lab/controls.py")
        self.assertTrue(registry["fake_quartet_control"]["default_run"])
        for key in ("s9_symbolic_identities", "dh_phi_identities",
                    "modular_thermality_moments", "levy_archimedean_kernel",
                    "pole_renorm_prime_measure"):
            self.assertIn(key, registry)
            self.assertFalse(registry[key]["default_run"])
            self.assertTrue((ROOT / registry[key]["path"]).is_file())


if __name__ == "__main__":
    unittest.main()
