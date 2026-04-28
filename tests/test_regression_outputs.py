import filecmp
import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from tests.helpers.compare_newick_semantics import compare_tree_files


RUN_TREESHRINK = REPO_ROOT / "run_treeshrink.py"
INPUT_ROOT = REPO_ROOT / "tests" / "fixtures" / "inputs"
GOLDEN_ROOT = REPO_ROOT / "tests" / "fixtures" / "golden"
QUANTILES = "0.01 0.05 0.1"


class GoldenOutputRegressionTests(unittest.TestCase):
    def run_treeshrink(self, dataset, mode, outdir):
        return subprocess.run(
            [
                sys.executable,
                str(RUN_TREESHRINK),
                "-t",
                str(INPUT_ROOT / ("%s.trees" % dataset)),
                "-o",
                str(outdir),
                "-O",
                "output",
                "-m",
                mode,
                "-q",
                QUANTILES,
                "-f",
            ],
            cwd=str(REPO_ROOT),
            text=True,
            capture_output=True,
        )

    def assert_golden_case(self, dataset, mode, golden_name):
        with tempfile.TemporaryDirectory() as tmpdir:
            outdir = Path(tmpdir) / "out"
            result = self.run_treeshrink(dataset, mode, outdir)

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Output files written to %s with prefix output." % outdir, result.stdout)
            self.compare_output_dir(GOLDEN_ROOT / golden_name, outdir)

    def test_mm_allgenes_matches_v139(self):
        self.assert_golden_case("mm", "all-genes", "mm_allgenes_v139")

    def test_mm_pergene_matches_v139(self):
        self.assert_golden_case("mm", "per-gene", "mm_pergene_v139")

    def test_mm_perspecies_matches_v139(self):
        self.assert_golden_case("mm", "per-species", "mm_perspecies_v139")

    def test_kp_allgenes_matches_v139(self):
        self.assert_golden_case("kp", "all-genes", "kp_allgenes_v139")

    def test_kp_pergene_matches_v139(self):
        self.assert_golden_case("kp", "per-gene", "kp_pergene_v139")

    def test_kp_perspecies_matches_v139(self):
        self.assert_golden_case("kp", "per-species", "kp_perspecies_v139")

    def test_frogs_allgenes_matches_v139(self):
        self.assert_golden_case("frogs", "all-genes", "frogs_allgenes_v139")

    def test_frogs_pergene_matches_v139(self):
        self.assert_golden_case("frogs", "per-gene", "frogs_pergene_v139")

    def test_frogs_perspecies_matches_v139(self):
        self.assert_golden_case("frogs", "per-species", "frogs_perspecies_v139")

    def compare_output_dir(self, golden_dir, actual_dir):
        golden_files = self.output_files(golden_dir)
        actual_files = self.output_files(actual_dir)

        self.assertEqual(
            sorted(str(path) for path in golden_files),
            sorted(str(path) for path in actual_files),
            "generated output file set differs from golden fixture",
        )

        for rel_path in golden_files:
            golden_path = golden_dir / rel_path
            actual_path = actual_dir / rel_path
            if rel_path.suffix == ".trees":
                mismatches = compare_tree_files(golden_path, actual_path, length_tol=1e-9)
                if mismatches:
                    self.fail("%s Newick semantic mismatch:\n%s" % (rel_path, "\n".join(mismatches[:10])))
            else:
                self.assertTrue(
                    filecmp.cmp(golden_path, actual_path, shallow=False),
                    "%s differs from golden fixture" % rel_path,
                )

    def output_files(self, directory):
        files = set()
        for root, _, filenames in os.walk(directory):
            root_path = Path(root)
            for filename in filenames:
                path = root_path / filename
                rel_path = path.relative_to(directory)
                if rel_path.suffix == ".log":
                    continue
                files.add(rel_path)
        return files


if __name__ == "__main__":
    unittest.main()
