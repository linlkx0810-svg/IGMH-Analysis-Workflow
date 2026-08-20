import csv
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run_workflow


class WorkflowTests(unittest.TestCase):
    def test_fragment_to_vmd_indices(self):
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1-72"), "index 0 to 71")
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1-3,5"), "index 0 to 2 4")
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1\u20133,5"), "index 0 to 2 4")
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1\u20143,5"), "index 0 to 2 4")

    def test_parse_systems(self):
        systems = run_workflow.parse_systems(["TS3a=TS3a_SP_PCM.fchk", "TS3b=TS3b_SP_PCM.fchk"])
        self.assertEqual(systems["TS3a"], "TS3a_SP_PCM.fchk")
        self.assertEqual(systems["TS3b"], "TS3b_SP_PCM.fchk")

    def test_read_fragment_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "fragments.csv"
            with path.open("w", newline="", encoding="utf-8") as handle:
                writer = csv.writer(handle)
                writer.writerow(["system", "fragment1", "fragment2", "center_selection"])
                writer.writerow(["TS3a", "1-72", "73-94", "index 0 to 93"])
            configs = run_workflow.read_fragment_file(path)
        self.assertEqual(configs["TS3a"].fragment1, "1-72")
        self.assertEqual(configs["TS3a"].fragment2, "73-94")
        self.assertEqual(configs["TS3a"].center_selection, "index 0 to 93")

    def test_combine_panels(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            left = tmp_path / "left.png"
            right = tmp_path / "right.png"
            out_clean = tmp_path / "combined_clean.png"
            out_labeled = tmp_path / "combined_labeled.png"
            Image.new("RGB", (100, 120), "white").save(left)
            Image.new("RGB", (120, 100), "white").save(right)

            run_workflow.combine_panels(left, right, out_clean, "A", "B", add_labels=False, dry_run=False)
            run_workflow.combine_panels(left, right, out_labeled, "A", "B", add_labels=True, dry_run=False)

            self.assertTrue(out_clean.exists())
            self.assertTrue(out_labeled.exists())
            with Image.open(out_clean) as img:
                self.assertEqual(img.size, (224, 100))
            with Image.open(out_labeled) as img:
                self.assertGreater(img.height, 100)

    def test_help_command(self):
        result = subprocess.run(
            [sys.executable, str(ROOT / "run_workflow.py"), "--help"],
            text=True,
            stdout=subprocess.PIPE,
            stderr=subprocess.PIPE,
            check=False,
        )
        self.assertEqual(result.returncode, 0, result.stderr)
        self.assertIn("Multiwfn IGMH", result.stdout)


if __name__ == "__main__":
    unittest.main()