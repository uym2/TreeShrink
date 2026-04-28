import os
import subprocess
import sys
import tempfile
import unittest
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
RUN_TREESHRINK = REPO_ROOT / "run_treeshrink.py"


def write_tree_fixture(path):
    with open(path, "w") as tree_file:
        for i in range(1, 11):
            tree_file.write(
                "((A:0.%d,B:0.2):0.1,"
                "(C:0.3,(D:0.4,(E:0.5,(F:0.6,(G:0.7,H:0.8):0.1):0.1):0.1):0.1):0.1);\n"
                % i
            )


class RuntimeLoggingTests(unittest.TestCase):
    def run_treeshrink(self, *args, cwd=None):
        return subprocess.run(
            [sys.executable, str(RUN_TREESHRINK), *args],
            cwd=str(cwd or REPO_ROOT),
            text=True,
            capture_output=True,
        )

    def test_normal_run_writes_log_and_preserves_terminal_output(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            tree_file = tmpdir / "input.trees"
            outdir = tmpdir / "out"
            write_tree_fixture(tree_file)

            result = self.run_treeshrink(
                "-t", str(tree_file),
                "-o", str(outdir),
                "-O", "output",
                "-f",
            )

            self.assertEqual(result.returncode, 0, result.stderr)
            self.assertIn("Launching TREESHRINK version", result.stdout)
            self.assertIn("There are only 10 gene trees in the dataset.", result.stdout)
            self.assertIn("Writing output", result.stdout)
            self.assertIn("Output files written to %s with prefix output." % outdir, result.stdout)

            log_text = (outdir / "output.log").read_text()
            self.assertIn("Launching TREESHRINK version", log_text)
            self.assertIn("TREESHRINK was called as follow", log_text)
            self.assertIn("There are only 10 gene trees in the dataset.", log_text)
            self.assertIn("Writing output", log_text)
            self.assertIn("Output files written to %s with prefix output." % outdir, log_text)

    def test_repeated_run_without_force_uses_incremented_log_prefix(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            tree_file = tmpdir / "input.trees"
            outdir = tmpdir / "out"
            write_tree_fixture(tree_file)

            first = self.run_treeshrink("-t", str(tree_file), "-o", str(outdir), "-O", "output", "-f")
            self.assertEqual(first.returncode, 0, first.stderr)

            second = self.run_treeshrink("-t", str(tree_file), "-o", str(outdir), "-O", "output")
            self.assertEqual(second.returncode, 0, second.stderr)

            self.assertTrue((outdir / "output.log").is_file())
            self.assertTrue((outdir / "output1.log").is_file())
            self.assertIn("Automatically changes prefix to 'output1'", second.stdout)
            self.assertIn("Automatically changes prefix to 'output1'", (outdir / "output1.log").read_text())
            self.assertIn("Output files written to %s with prefix output1." % outdir, second.stdout)

    def test_force_uses_requested_log_prefix_and_overwrites_existing_log(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            tree_file = tmpdir / "input.trees"
            outdir = tmpdir / "out"
            write_tree_fixture(tree_file)

            first = self.run_treeshrink("-t", str(tree_file), "-o", str(outdir), "-O", "output", "-f")
            self.assertEqual(first.returncode, 0, first.stderr)

            log_path = outdir / "output.log"
            log_path.write_text("stale log contents\n")

            second = self.run_treeshrink("-t", str(tree_file), "-o", str(outdir), "-O", "output", "-f")
            self.assertEqual(second.returncode, 0, second.stderr)

            log_text = log_path.read_text()
            self.assertNotIn("stale log contents", log_text)
            self.assertIn("With --force, all existing files with prefix 'output' will be overrided", log_text)
            self.assertIn("Output files written to %s with prefix output." % outdir, log_text)

    def test_version_and_help_do_not_create_logs(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)

            version = self.run_treeshrink("-v", cwd=tmpdir)
            self.assertEqual(version.returncode, 0, version.stderr)
            self.assertEqual(version.stdout.strip(), "1.4.0")
            self.assertEqual(list(tmpdir.iterdir()), [])

            help_result = self.run_treeshrink(cwd=tmpdir)
            self.assertEqual(help_result.returncode, 0, help_result.stderr)
            self.assertIn("usage:", help_result.stdout)
            self.assertEqual(list(tmpdir.iterdir()), [])

    def test_exception_after_logging_starts_is_written_to_log_and_stderr(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            tmpdir = Path(tmpdir)
            outdir = tmpdir / "out"
            missing_tree = tmpdir / "missing.trees"

            result = self.run_treeshrink(
                "-t", str(missing_tree),
                "-o", str(outdir),
                "-O", "output",
                "-f",
            )

            self.assertNotEqual(result.returncode, 0)
            self.assertIn("FileNotFoundError", result.stderr)

            log_text = (outdir / "output.log").read_text()
            self.assertIn("Launching TREESHRINK version", log_text)
            self.assertIn("TREESHRINK was called as follow", log_text)
            self.assertIn("FileNotFoundError", log_text)


if __name__ == "__main__":
    unittest.main()
