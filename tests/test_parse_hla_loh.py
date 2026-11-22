"""parse_hla_loh.py: SpecHLA merge.hla.copy.txt -> tidy per-locus calls."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "parse_hla_loh.py"
spec = importlib.util.spec_from_file_location("parse_hla_loh", BIN)
parse_hla_loh = importlib.util.module_from_spec(spec)
spec.loader.exec_module(parse_hla_loh)

HEADER = ["Sample", "HLA", "Allele1", "Allele2", "copyratio", "KeptHLA", "LossHLA",
          "Freq1", "Freq2", "Purity", "Het_num", "LOH", "SampleDir", "Homozygous", "Class"]

# Rows shaped like real cal.hla.copy.pl output.
ROWS = [
    # heterozygous, both alleles present -> retained
    ["P1", "A", "A*11:01:01:01", "A*02:01:01:01", "1:1", "A*02:01:01:01", "A*11:01:01:01",
     "0.42", "0.58", "0.5", "99", "N", "P1", "N", "I"],
    # one allele at zero copies -> loss of C*04:01
    ["P1", "C", "C*01:02:01:01", "C*04:01:01:01", "2:0", "C*01:02:01:01", "C*04:01:01:01",
     "0.622", "0.378", "0.5", "36", "Y", "P1", "N", "I"],
    # imbalance only (no zero side), SpecHLA still says Y
    ["P1", "B", "B*07:02:01:01", "B*35:01:01:01", "3:1", "B*07:02:01:01", "B*35:01:01:01",
     "0.7", "0.3", "0.5", "58", "Y", "P1", "N", "I"],
    # homozygous at two fields -> nothing to lose, not evaluable
    ["P1", "DQB1", "DQB1*06:02:01:01", "DQB1*06:02:01:01", "1:3", "DQB1*06:02:01:01",
     "DQB1*06:02:01:01", "0.216", "0.784", "0.5", "30", "N", "P1", "Y", "II"],
    # too few heterozygous sites -> not evaluable
    ["P1", "DRB1", "DRB1*15:01:01:01", "DRB1*03:01:01:01", "2:0", "DRB1*15:01:01:01",
     "DRB1*03:01:01:01", "0.9", "0.1", "0.5", "3", "Y", "P1", "N", "II"],
]


def write_loh(rows, header=HEADER, delimiter="\t"):
    path = Path(tempfile.mkdtemp()) / "merge.hla.copy.txt"
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter=delimiter, lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def run(rows, **kwargs):
    path = write_loh(rows)
    out = path.parent / "out.tsv"
    argv = [str(path), "--patient", "P1", "-o", str(out)]
    for key, value in kwargs.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    parse_hla_loh.main(argv)
    with out.open(newline="") as handle:
        return {row["locus"]: row for row in csv.DictReader(handle, delimiter="\t")}


class TestNaming(unittest.TestCase):
    def test_two_field_truncation(self):
        self.assertEqual(parse_hla_loh.two_field("A*02:01:01:01"), "A*02:01")
        self.assertEqual(parse_hla_loh.two_field("DQA1*01:02"), "DQA1*01:02")

    def test_class_i_gets_the_hla_prefix(self):
        self.assertEqual(parse_hla_loh.pvacseq_name("A*02:01:01:01"), "HLA-A*02:01")
        self.assertEqual(parse_hla_loh.pvacseq_name("C*04:01:01:01"), "HLA-C*04:01")

    def test_class_ii_does_not(self):
        self.assertEqual(parse_hla_loh.pvacseq_name("DRB1*15:01:01:01"), "DRB1*15:01")
        self.assertEqual(parse_hla_loh.pvacseq_name("DQB1*06:02:01:01"), "DQB1*06:02")

    def test_unparseable_allele(self):
        self.assertEqual(parse_hla_loh.pvacseq_name(""), "")
        self.assertEqual(parse_hla_loh.pvacseq_name("no_match"), "")

    def test_already_prefixed_allele_is_not_doubled(self):
        self.assertEqual(parse_hla_loh.pvacseq_name("HLA-A*02:01"), "HLA-A*02:01")


class TestCalls(unittest.TestCase):
    def setUp(self):
        self.table = run(ROWS)

    def test_retained_locus(self):
        self.assertEqual(self.table["A"]["call"], "retained")
        self.assertEqual(self.table["A"]["evaluable"], "TRUE")
        self.assertEqual(self.table["A"]["lost_allele_pvacseq"], "")

    def test_copy_loss_names_the_lost_allele(self):
        row = self.table["C"]
        self.assertEqual(row["call"], "loss")
        self.assertEqual(row["lost_allele_pvacseq"], "HLA-C*04:01")
        self.assertEqual(row["kept_allele_pvacseq"], "HLA-C*01:02")
        self.assertEqual((row["cn1"], row["cn2"]), ("2", "0"))

    def test_imbalance_is_not_a_loss_by_default(self):
        self.assertEqual(self.table["B"]["call"], "imbalance")
        self.assertEqual(self.table["B"]["lost_allele_pvacseq"], "")

    def test_imbalance_becomes_a_loss_under_the_spechla_definition(self):
        table = run(ROWS, definition="spechla")
        self.assertEqual(table["B"]["call"], "loss")
        self.assertEqual(table["B"]["lost_allele_pvacseq"], "HLA-B*35:01")

    def test_homozygous_locus_is_not_evaluable(self):
        row = self.table["DQB1"]
        self.assertEqual(row["call"], "not_evaluable")
        self.assertEqual(row["heterozygous_2field"], "FALSE")

    def test_too_few_heterozygous_sites_is_not_evaluable(self):
        self.assertEqual(self.table["DRB1"]["call"], "not_evaluable")

    def test_lowering_min_het_makes_that_locus_callable(self):
        table = run(ROWS, min_het=1)
        self.assertEqual(table["DRB1"]["call"], "loss")
        self.assertEqual(table["DRB1"]["lost_allele_pvacseq"], "DRB1*03:01")

    def test_hla_class_column(self):
        self.assertEqual(self.table["A"]["hla_class"], "I")
        self.assertEqual(self.table["DRB1"]["hla_class"], "II")

    def test_patient_column(self):
        self.assertEqual(self.table["A"]["patient"], "P1")


class TestPurityFilter(unittest.TestCase):
    def low_purity_rows(self):
        rows = [list(row) for row in ROWS]
        for row in rows:
            row[9] = "0.1"
        return rows

    def test_low_purity_is_not_called(self):
        table = run(self.low_purity_rows())
        self.assertEqual(table["C"]["call"], "low_purity")
        self.assertEqual(table["C"]["evaluable"], "FALSE")

    def test_disabling_the_purity_filter(self):
        table = run(self.low_purity_rows(), min_purity=0)
        self.assertEqual(table["C"]["call"], "loss")


class TestInputShapes(unittest.TestCase):
    def test_long_collated_csv(self):
        header = ["cohort", "sample_id", "locus", "allele1", "allele2", "copyratio",
                  "kept", "lost", "freq1", "freq2", "purity", "het_num", "loh"]
        rows = [["NeoTrio", "P1", "C", "C*01:02:01:01", "C*04:01:01:01", "2:0",
                 "C*01:02:01:01", "C*04:01:01:01", "0.6", "0.4", "0.5", "36", "Y"]]
        path = write_loh(rows, header=header, delimiter=",")
        out = path.parent / "out.tsv"
        parse_hla_loh.main([str(path), "--patient", "P1", "-o", str(out)])
        with out.open(newline="") as handle:
            table = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(table[0]["call"], "loss")
        self.assertEqual(table[0]["lost_allele_pvacseq"], "HLA-C*04:01")

    def test_locus_column_is_required(self):
        path = write_loh([["P1", "x"]], header=["Sample", "Something"])
        with self.assertRaises(SystemExit):
            parse_hla_loh.main([str(path), "-o", str(path.parent / "out.tsv")])

    def test_a_differently_named_single_sample_is_kept(self):
        rows = [list(ROWS[1])]
        rows[0][0] = "SpecHLA_run_7"
        table = run(rows)
        self.assertEqual(table["C"]["call"], "loss")
        self.assertEqual(table["C"]["patient"], "P1")

    def test_several_samples_with_no_match_is_an_error(self):
        rows = [list(ROWS[0]), list(ROWS[1])]
        rows[0][0], rows[1][0] = "other_a", "other_b"
        path = write_loh(rows)
        with self.assertRaises(SystemExit):
            parse_hla_loh.main([str(path), "--patient", "P1", "-o", str(path.parent / "out.tsv")])

    def test_empty_file(self):
        path = Path(tempfile.mkdtemp()) / "empty.txt"
        path.write_text("")
        out = path.parent / "out.tsv"
        parse_hla_loh.main([str(path), "-o", str(out)])
        self.assertEqual(out.read_text().strip().split("\t")[0], "patient")


if __name__ == "__main__":
    unittest.main()
