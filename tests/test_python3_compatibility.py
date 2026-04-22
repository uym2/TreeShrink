import sys
import tempfile
import unittest
import inspect
from pathlib import Path


REPO_ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(REPO_ROOT))

from treeshrink.alignment import MultiLocusDataset
from treeshrink.sequence_lib import hash_taxon_seq
from treeshrink._vendor import dendropy


class Python3CompatibilityTests(unittest.TestCase):
    def write_fasta(self, directory, filename="input.fasta"):
        path = Path(directory) / filename
        with open(path, "w") as fasta:
            fasta.write(">taxon_one\n")
            fasta.write("AC-GT\n")
            fasta.write(">taxon_two\n")
            fasta.write("A--GT\n")
        return path

    def test_hash_taxon_seq_uses_python3_file_iteration(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self.write_fasta(tmpdir)

            self.assertEqual(
                hash_taxon_seq(str(fasta_path)),
                {
                    "taxon_one": "ACGT",
                    "taxon_two": "AGT",
                },
            )

    def test_multilocus_dataset_reads_fasta_with_text_mode(self):
        with tempfile.TemporaryDirectory() as tmpdir:
            fasta_path = self.write_fasta(tmpdir)
            dataset = MultiLocusDataset()

            dataset.read_files([str(fasta_path)], datatype="DNA")

            self.assertEqual(len(dataset), 1)
            self.assertEqual(dataset[0].datatype, "DNA")
            labels = dataset[0].dataset.taxon_namespaces[0].labels()
            self.assertEqual(set(labels), {"taxon_one", "taxon_two"})

    def test_dendropy_is_private_vendor_package(self):
        dendropy_path = Path(inspect.getfile(dendropy)).resolve()

        self.assertIn("treeshrink", dendropy_path.parts)
        self.assertIn("_vendor", dendropy_path.parts)
        self.assertEqual(dendropy.__name__, "treeshrink._vendor.dendropy")


if __name__ == "__main__":
    unittest.main()
