import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("parse_hlahd_result", Path(__file__).resolve().parents[1] / "bin" / "parse_hlahd_result.py")
hlahd = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(hlahd)

RESULT = """A\tHLA-A*24:02:01:01\tHLA-A*02:07:01
B\tHLA-B*40:06:01:01\tHLA-B*40:06:01:02
C\tHLA-C*08:01:01\tHLA-C*15:02:01
E\tHLA-E*01:01:01\tHLA-E*01:03:02
DRB1\tHLA-DRB1*09:01:02\tHLA-DRB1*15:02:01
DRB3\tNot typed
DRB4\tHLA-DRB4*01:03:01\t-
DQA1\tHLA-DQA1*01:03:01\tHLA-DQA1*03:02
DQB1\tHLA-DQB1*06:01:01\tHLA-DQB1*03:03:02
DPA1\tHLA-DPA1*02:02:02\t-
DPB1\tHLA-DPB1*05:01:01\tHLA-DPB1*09:01:01
"""


class ParseHlaHd(unittest.TestCase):
    def test_two_field_names_in_pvacseq_format(self):
        alleles = hlahd.pvacseq_alleles(hlahd.parse_result(RESULT))
        self.assertEqual(alleles[:5], ["HLA-A*24:02", "HLA-A*02:07", "HLA-B*40:06", "HLA-C*08:01", "HLA-C*15:02"])
        self.assertIn("DRB1*09:01", alleles)
        self.assertIn("DRB4*01:03", alleles)
        self.assertNotIn("HLA-E*01:01", alleles)

    def test_homozygous_collapses_to_one_allele(self):
        calls = hlahd.parse_result(RESULT)
        self.assertEqual(calls["B"], ["B*40:06"])
        self.assertEqual(calls["DRB3"], [])

    def test_class_ii_alpha_beta_pairs_include_trans(self):
        alleles = hlahd.pvacseq_alleles(hlahd.parse_result(RESULT))
        dq = [a for a in alleles if a.startswith("DQA1")]
        self.assertEqual(len(dq), 4)
        self.assertIn("DQA1*01:03-DQB1*03:03", dq)
        self.assertEqual([a for a in alleles if a.startswith("DPA1")], ["DPA1*02:02-DPB1*05:01", "DPA1*02:02-DPB1*09:01"])

    def test_null_alleles_dropped(self):
        calls = hlahd.parse_result("A\tHLA-A*02:01:01\tHLA-A*24:09N\n")
        self.assertEqual(calls["A"], ["A*02:01"])

    def test_missing_class_i_is_an_error(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = Path(tmp) / "x_final.result.txt"
            result.write_text("A\tHLA-A*02:01:01\t-\nB\tNot typed\nC\tHLA-C*07:02\t-\n")
            with self.assertRaises(SystemExit):
                hlahd.main([str(result), "--pvacseq-out", str(Path(tmp) / "a.txt"), "--table-out", str(Path(tmp) / "a.tsv")])

    def test_writes_comma_joined_allele_file(self):
        with tempfile.TemporaryDirectory() as tmp:
            result = Path(tmp) / "x_final.result.txt"
            result.write_text(RESULT)
            out = Path(tmp) / "alleles.txt"
            hlahd.main([str(result), "--pvacseq-out", str(out), "--table-out", str(Path(tmp) / "a.tsv")])
            self.assertTrue(out.read_text().startswith("HLA-A*24:02,HLA-A*02:07,"))


if __name__ == "__main__":
    unittest.main()
