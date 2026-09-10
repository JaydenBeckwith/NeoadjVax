import csv
import importlib.util
import json
from pathlib import Path
import re
import tempfile
import unittest


ROOT = Path(__file__).resolve().parent / "dna_somatic"
SPEC = importlib.util.spec_from_file_location("somatic_fixture", ROOT / "make_fixture.py")
fixture = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fixture)
CHECK_SPEC = importlib.util.spec_from_file_location("somatic_check", ROOT / "check_results.py")
checker = importlib.util.module_from_spec(CHECK_SPEC)
CHECK_SPEC.loader.exec_module(checker)


class SomaticFixture(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def test_reproducible_reference_and_longitudinal_samples(self):
        first = fixture.generate(self.root / "one", sam_only=True)
        second = fixture.generate(self.root / "two", sam_only=True)
        self.assertEqual((self.root / "one/genome.fa").read_bytes(), (self.root / "two/genome.fa").read_bytes())
        self.assertEqual(len(first["bams"]), 8)
        self.assertEqual({(row["patient_id"], row["timepoint"]) for row in first["rna_samples"]},
                         {("P1", "PRE"), ("P1", "DAY42"), ("P2", "PRE"), ("P2", "DAY42")})
        self.assertFalse(second["bam_generation_complete"])

    def test_sam_coordinates_pairs_sequences_and_spliced_reads(self):
        manifest = fixture.generate(self.root / "data", sam_only=True)
        for bam in manifest["bams"]:
            sam = Path(bam["path"]).with_suffix(".sam").read_text()
            self.assertIn("@HD\tVN:1.6\tSO:coordinate", sam)
            self.assertIn(f"SM:{bam['sample']}", sam)
            rows = [line.split("\t") for line in sam.splitlines() if not line.startswith("@")]
            self.assertEqual(len(rows), 48)
            positions = [int(row[3]) for row in rows]
            self.assertEqual(positions, sorted(positions))
            groups = {}
            for row in rows:
                groups.setdefault(row[0], []).append(row)
                consumed = sum(int(n) for n, operation in re.findall(r"(\d+)([MIDNSHP=X])", row[5]) if operation in "MIS=X")
                self.assertEqual(consumed, len(row[9]))
                self.assertEqual(len(row[9]), len(row[10]))
            self.assertEqual(len(groups), 24)
            for pair in groups.values():
                self.assertEqual({row[1] for row in pair}, {"99", "147"})
                self.assertEqual(pair[0][3], pair[1][7])
                self.assertEqual(pair[1][3], pair[0][7])
                self.assertEqual(int(pair[0][8]), -int(pair[1][8]))
            self.assertEqual(any("N" in row[5] for row in rows), bam["rna"])

    def test_will_not_overwrite_existing_data(self):
        marker = self.root / "important.txt"
        marker.write_text("preserve")
        with self.assertRaisesRegex(ValueError, "empty"):
            fixture.generate(self.root, sam_only=True)
        self.assertEqual(marker.read_text(), "preserve")

    def test_checker_rejects_incomplete_bam_generation(self):
        fixture.generate(self.root / "data", sam_only=True)
        with self.assertRaisesRegex(ValueError, "generated/validated"):
            checker.check(self.root / "data/manifest.json", self.root / "missing_trace", self.root)

    def test_checker_detects_silent_dropped_tasks(self):
        manifest = fixture.generate(self.root / "data", sam_only=True)
        manifest["bam_generation_complete"] = True
        path = self.root / "data/manifest.json"
        path.write_text(json.dumps(manifest))
        trace = self.root / "trace.tsv"
        trace.write_text("process\tstatus\texit\nDNA_VARIANT_CALLING:MUTECT2\tCOMPLETED\t0\n")
        with self.assertRaisesRegex(ValueError, "expected 4 tasks"):
            checker.check(path, trace, self.root)


if __name__ == "__main__":
    unittest.main()
