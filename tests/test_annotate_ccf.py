import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("annotate_ccf", Path(__file__).resolve().parents[1] / "bin" / "annotate_ccf.py")
ccf = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ccf)

SEGMENTS = '"chromosome"\t"start.pos"\t"end.pos"\t"Bf"\t"CNt"\t"A"\t"B"\n' \
           '"chr1"\t1\t1000\t0.5\t2\t1\t1\n' \
           '"chr1"\t2000\t5000\t0.33\t4\t3\t1\n'
REPORT = "ID\tGene\tDNA VAF\n" \
         "chr1-499-500-A-T\tG1\t0.25\n" \
         "chr1-2999-3000-C-G\tG2\t0.6\n" \
         "chr1-1499-1500-G-A\tG3\t0.1\n" \
         "chr2-99-100-T-C\tG4\tNA\n"


class AnnotateCcf(unittest.TestCase):
    def run_ccf(self, *extra):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "report.tsv").write_text(REPORT)
            (tmp / "segments.txt").write_text(SEGMENTS)
            args = [str(tmp / "report.tsv"), "--purity", "0.5", "-o", str(tmp / "out.tsv")]
            args += [a.replace("SEGMENTS", str(tmp / "segments.txt")) for a in extra]
            ccf.main(args)
            with (tmp / "out.tsv").open() as handle:
                return {row["Gene"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_clonal_heterozygous_on_diploid_segment(self):
        row = self.run_ccf("--segments", "SEGMENTS")["G1"]
        # VAF 0.25 at purity 0.5, CN 2: one mutated copy in every tumour cell
        self.assertEqual(row["CN Source"], "sequenza_segment")
        self.assertEqual(row["Mutated Copies"], "1.000")
        self.assertEqual(row["CCF"], "1.000")

    def test_amplified_segment_raises_multiplicity(self):
        row = self.run_ccf("--segments", "SEGMENTS")["G2"]
        # u = 0.6 * (0.5*4 + 2*0.5) / 0.5 = 3.6 -> m capped at major CN 3 -> CCF 1
        self.assertEqual(row["Tumor CN"], "4")
        self.assertEqual(row["Multiplicity"], "3")
        self.assertEqual(row["CCF"], "1.000")

    def test_falls_back_to_ploidy_between_segments(self):
        row = self.run_ccf("--segments", "SEGMENTS", "--ploidy", "3")["G3"]
        # u = 0.1 * (0.5*3 + 1) / 0.5 = 0.5 -> subclonal
        self.assertEqual(row["CN Source"], "ploidy")
        self.assertEqual(row["CCF"], "0.500")

    def test_without_ploidy_or_segments_assumes_normal_copy_number(self):
        rows = self.run_ccf()
        self.assertEqual(rows["G1"]["CN Source"], "assumed_normal")
        self.assertEqual(rows["G4"]["CCF"], "NA")

    def test_male_x_uses_single_normal_copy(self):
        self.assertEqual(ccf.ccf_columns(0.5, "chrX", 10, 1.0, None, "male", None)[4], "1")
        self.assertEqual(ccf.ccf_columns(0.5, "chrX", 10, 1.0, None, "female", None)[4], "2")


if __name__ == "__main__":
    unittest.main()
