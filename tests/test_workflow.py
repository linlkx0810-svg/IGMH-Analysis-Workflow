import csv
import json
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path

from PIL import Image, ImageDraw

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

import run_workflow
from igmh_alignment import apply_transform, kabsch_align
from igmh_comparison import build_comparison_plan, load_comparison_config
from igmh_fragments import FragmentValidationError, fragment_to_vmd_indices, validate_fragments
from igmh_pair_parser import parse_pair_value_file, read_pair_outputs, rank_pairs
from igmh_validation import final_status, inspect_png_framing, pair_integration_status, validate_pair_outputs, write_manifest


def write_cube(path: Path, coords):
    lines = [
        "Synthetic cube",
        "For tests",
        f" {len(coords)} 0.0 0.0 0.0",
        " 1 1.0 0.0 0.0 0.0",
        " 1 0.0 1.0 0.0 0.0",
        " 1 0.0 0.0 1.0 0.0",
    ]
    for atomic_number, xyz in coords:
        lines.append(f" {atomic_number} 0.0 {xyz[0]} {xyz[1]} {xyz[2]}")
    lines.append(" 0.0")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


class WorkflowTests(unittest.TestCase):
    def test_fragment_to_vmd_indices(self):
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1-72"), "index 0 to 71")
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1-3,5"), "index 0 to 2 4")
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1\u20133,5"), "index 0 to 2 4")
        self.assertEqual(run_workflow.fragment_to_vmd_indices("1\u20143,5"), "index 0 to 2 4")
        self.assertEqual(fragment_to_vmd_indices("1-3,5-6,8"), "index 0 to 2 4 to 5 7")

    def test_fragment_validation_non_contiguous_and_uncovered(self):
        result = validate_fragments("1-3,5", "6,8-9", 10)
        self.assertEqual(result.fragment1_atoms, (1, 2, 3, 5))
        self.assertEqual(result.fragment2_atoms, (6, 8, 9))
        self.assertEqual(result.uncovered_atoms, (4, 7, 10))

    def test_fragment_overlap_detection(self):
        with self.assertRaises(FragmentValidationError):
            validate_fragments("1-5", "5-8", 10)

    def test_fragment_out_of_range_detection(self):
        with self.assertRaises(FragmentValidationError):
            validate_fragments("1-5", "11", 10)

    def test_fragment_duplicate_detection(self):
        with self.assertRaises(FragmentValidationError):
            validate_fragments("1-3,3", "4", 10)

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

    def test_pair_integration_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            log = Path(tmp) / "out.log"
            log.write_text("Integrating grid... 85.2 %\n", encoding="utf-8")
            self.assertEqual(pair_integration_status(log).status, "FAIL")
            log.write_text("Integrating grid... 100.0 %\n", encoding="utf-8")
            self.assertEqual(pair_integration_status(log).status, "PASS")

    def test_missing_pair_outputs_fail(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            log = tmp_path / "out.log"
            log.write_text("100.0 %\n", encoding="utf-8")
            statuses = validate_pair_outputs(tmp_path, log, parse_pair_value_file)
        self.assertEqual([status.status for status in statuses], ["PASS", "FAIL", "FAIL"])

    def test_pair_parsers_and_ranking(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            (tmp_path / "atmdg.txt").write_text("Atomic pair delta-g indices\n 79 119 : 0.0123\n 5 116 : 0.0200\n", encoding="utf-8")
            (tmp_path / "IBSIW.txt").write_text("IBSIW index\n 79 119 : 0.51\n 5 116 : 0.80\n", encoding="utf-8")
            pairs = parse_pair_value_file(tmp_path / "atmdg.txt")
            self.assertIn((79, 119), pairs)
            records = read_pair_outputs(tmp_path)
            ranked = rank_pairs(records, top_n=1)
        self.assertEqual(ranked[0].key, (5, 116))
        self.assertAlmostEqual(ranked[0].deltaGpair, 0.0200)
        self.assertAlmostEqual(ranked[0].IBSIW, 0.80)

    def test_kabsch_alignment_and_rmsd(self):
        reference = [[0, 0, 0], [1, 0, 0], [0, 1, 0], [0, 0, 1]]
        target = [[1, 2, 3], [1, 3, 3], [0, 2, 3], [1, 2, 4]]
        result = kabsch_align(reference, target)
        aligned = apply_transform(target, result.rotation, result.translation)
        self.assertLess(result.rmsd, 1e-10)
        self.assertAlmostEqual(float(aligned[1][0]), 1.0, places=8)

    def test_comparison_config_alignment_camera_scale_metadata(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            ref_dir = tmp_path / "Ref_IGMH_files"
            tar_dir = tmp_path / "Tar_IGMH_files"
            ref_dir.mkdir()
            tar_dir.mkdir()
            ref_coords = [(6, (0, 0, 0)), (1, (1, 0, 0)), (1, (0, 1, 0)), (1, (0, 0, 1))]
            tar_coords = [(6, (1, 2, 3)), (1, (2, 2, 3)), (1, (1, 3, 3)), (1, (1, 2, 4))]
            write_cube(ref_dir / "sl2r.cub", ref_coords)
            write_cube(tar_dir / "sl2r.cub", tar_coords)
            config = tmp_path / "comparison.json"
            config.write_text(json.dumps({
                "reference_system": "Ref",
                "panel_order": ["Ref", "Tar"],
                "alignment_mappings": {"Tar": {"reference_atoms": [1, 2, 3, 4], "target_atoms": [1, 2, 3, 4]}},
                "view_matrix": [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]],
            }), encoding="utf-8")
            self.assertEqual(load_comparison_config(config)["reference_system"], "Ref")
            plan = build_comparison_plan(config, {"Ref": ref_dir, "Tar": tar_dir}, tmp_path)
            self.assertTrue(plan.transform_tcl.exists())
        self.assertEqual(plan.panel_order, ["Ref", "Tar"])
        self.assertLess(plan.alignment_rmsd_A["Tar"], 1e-8)
        self.assertGreater(plan.scale, 0)
        self.assertIn("VMD molecule global_matrix", plan.volume_registration)

    def test_framing_validation_thresholds(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "framed.png"
            img = Image.new("RGB", (100, 100), "white")
            draw = ImageDraw.Draw(img)
            draw.rectangle((5, 5, 94, 94), fill="black")
            img.save(path)
            statuses, summary = inspect_png_framing(path, expected_size=(100, 100))
        self.assertEqual(statuses[-1].status, "PASS")
        self.assertGreaterEqual(summary["min_margin_percent"], 3.0)

    def test_combine_panels_and_order(self):
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
            labels = run_workflow.ordered_labels({"B": right, "A": left}, ["A", "B"])
            self.assertEqual(labels, ["A", "B"])
            with Image.open(out_clean) as img:
                self.assertEqual(img.size, (224, 100))
            with Image.open(out_labeled) as img:
                self.assertGreater(img.height, 100)

    def test_manifest_generation_and_strict_status(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "manifest.json"
            write_manifest(path, {"final_status": "FAIL", "comparison": {"camera_shared": True}})
            data = json.loads(path.read_text(encoding="utf-8"))
        self.assertEqual(data["final_status"], "FAIL")
        self.assertTrue(data["comparison"]["camera_shared"])
        self.assertEqual(final_status([]), "PASS")

    def test_strict_failure_behavior(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp_path = Path(tmp)
            rc = run_workflow.main([
                "--validate-only", "--strict", "--input-dir", str(tmp_path), "--output-dir", str(tmp_path),
                "--figures-dir", str(tmp_path), "--manifest", str(tmp_path / "manifest.json"),
                "--systems", "Missing=missing.fchk",
            ])
        self.assertEqual(rc, 1)

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
        self.assertIn("--panel-order", result.stdout)
        self.assertIn("--pair-analysis", result.stdout)
        self.assertIn("--comparison-config", result.stdout)


if __name__ == "__main__":
    unittest.main()
