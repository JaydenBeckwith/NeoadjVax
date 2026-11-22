"""agretopicity.py and pvactools_output.py: the standalone calculation.

Both are meant to be run by hand against a pVACseq or pVACfuse output
directory, so the tests drive the command line the same way rather than
only the functions behind it.
"""
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(BIN))


def load(name):
    spec = importlib.util.spec_from_file_location(name, BIN / f"{name}.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


agr = load("agretopicity")
pvo = load("pvactools_output")

SEQ_COLUMNS = ["ID", "Gene", "AA Change", "Allele", "Best Peptide",
               "IC50 MT", "IC50 WT", "%ile MT", "%ile WT", "DNA VAF", "Tier"]
FUSE_COLUMNS = ["ID", "Gene", "Best Peptide", "Allele", "IC50 MT", "%ile MT",
                "Expr", "Read Support", "Tier"]


def seq_row(id_, mt, wt, gene="BRAF"):
    return {**{c: "NA" for c in SEQ_COLUMNS},
            "ID": id_, "Gene": gene, "IC50 MT": mt, "IC50 WT": wt, "DNA VAF": "0.4"}


def write_report(rows, columns, path):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


def read(path):
    with open(path, newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def pvacseq_tree(root, rows=None):
    """A pVACseq output directory with all three class reports."""
    rows = rows if rows is not None else [seq_row("chr1-100-100-A-T", "12", "900")]
    for sub, label in (("MHC_Class_I", "MHC_I"), ("MHC_Class_II", "MHC_II"), ("combined", "Combined")):
        write_report(rows, SEQ_COLUMNS, root / "pvacseq" / sub / f"tumor.{label}.all_epitopes.aggregated.tsv")
    return root / "pvacseq"


def pvacfuse_tree(root):
    rows = [{**{c: "NA" for c in FUSE_COLUMNS}, "ID": "EML4--ALK_chr2:42295516",
             "Gene": "EML4--ALK", "IC50 MT": "35", "Read Support": "30", "Expr": "2.4"}]
    write_report(rows, FUSE_COLUMNS, root / "pvacfuse" / "MHC_Class_I" / "P1.MHC_I.all_epitopes.aggregated.tsv")
    return root / "pvacfuse"


class TestCalculation(unittest.TestCase):
    def test_a_foreign_mutant(self):
        ratio, dai = agr.agretopicity(12.0, 900.0)
        self.assertAlmostEqual(ratio, 12 / 900)
        self.assertAlmostEqual(dai, 888)

    def test_a_mutant_no_better_than_the_wild_type(self):
        ratio, _ = agr.agretopicity(420.0, 450.0)
        self.assertGreater(ratio, 0.5)

    def test_a_mutant_that_binds_worse_than_the_wild_type(self):
        ratio, dai = agr.agretopicity(1000.0, 100.0)
        self.assertEqual(ratio, 10.0)
        self.assertEqual(dai, -900)

    def test_a_missing_wild_type_prediction(self):
        self.assertEqual(agr.agretopicity(12.0, None), (None, None))

    def test_a_missing_mutant_prediction(self):
        self.assertEqual(agr.agretopicity(None, 900.0), (None, None))

    def test_a_zero_wild_type_ic50_is_degenerate_not_infinite(self):
        ratio, dai = agr.agretopicity(12.0, 0.0)
        self.assertIsNone(ratio)
        self.assertEqual(dai, -12)

    def test_dai_separates_what_the_ratio_cannot(self):
        # Both ratios are close to 1; the differences are 100 and 1.
        weak = agr.agretopicity(900.0, 1000.0)
        strong = agr.agretopicity(9.0, 10.0)
        self.assertAlmostEqual(weak[0], strong[0])
        self.assertNotAlmostEqual(weak[1], strong[1])


class TestClassify(unittest.TestCase):
    def test_bands(self):
        self.assertEqual(agr.classify(0.01, 0.1, 0.5), "foreign")
        self.assertEqual(agr.classify(0.3, 0.1, 0.5), "intermediate")
        self.assertEqual(agr.classify(0.9, 0.1, 0.5), "self_like")
        self.assertEqual(agr.classify(None, 0.1, 0.5), "NA")

    def test_boundaries_are_inclusive(self):
        self.assertEqual(agr.classify(0.1, 0.1, 0.5), "foreign")
        self.assertEqual(agr.classify(0.5, 0.1, 0.5), "self_like")


class TestFindReports(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())

    def test_a_pvacseq_directory_gives_all_three_classes(self):
        reports = pvo.find_reports(pvacseq_tree(self.root))
        self.assertEqual(sorted(label for label, _ in reports), ["Combined", "MHC_I", "MHC_II"])

    def test_a_single_file_is_taken_as_given(self):
        path = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, self.root / "x.MHC_I.tsv")
        self.assertEqual(pvo.find_reports(path), [("MHC_I", path)])

    def test_a_directory_with_nothing_in_it(self):
        self.assertEqual(pvo.find_reports(self.root), [])

    def test_a_path_that_does_not_exist(self):
        self.assertEqual(pvo.find_reports(self.root / "nope"), [])

    def test_a_flat_folder_of_reports(self):
        write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS,
                     self.root / "flat" / "tumor.MHC_I.all_epitopes.aggregated.tsv")
        self.assertEqual(len(pvo.find_reports(self.root / "flat")), 1)


class TestKind(unittest.TestCase):
    def test_a_somatic_report(self):
        self.assertEqual(pvo.kind(SEQ_COLUMNS), "somatic")

    def test_a_fusion_report(self):
        self.assertEqual(pvo.kind(FUSE_COLUMNS), "fusion")

    def test_an_empty_header(self):
        self.assertEqual(pvo.kind([]), "unknown")

    def test_something_else_entirely(self):
        self.assertEqual(pvo.kind(["a", "b"]), "unknown")

    def test_wild_type_columns_decide_it(self):
        self.assertEqual(pvo.kind(["ID", "IC50 MT", "IC50 WT", "Expr"]), "somatic")


class TestCommandLine(unittest.TestCase):
    def setUp(self):
        self.root = Path(tempfile.mkdtemp())

    def test_one_report_in_one_out(self):
        path = write_report([seq_row("v1", "12", "900"), seq_row("v2", "420", "450")],
                            SEQ_COLUMNS, self.root / "in.Combined.tsv")
        out = self.root / "out.tsv"
        agr.main([str(path), "-o", str(out)])
        rows = read(out)
        self.assertEqual(rows[0]["Agretopicity Class"], "foreign")
        self.assertEqual(rows[1]["Agretopicity Class"], "self_like")

    def test_the_original_columns_survive(self):
        path = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, self.root / "in.Combined.tsv")
        out = self.root / "out.tsv"
        agr.main([str(path), "-o", str(out)])
        self.assertEqual(read(out)[0]["Gene"], "BRAF")
        self.assertEqual(read(out)[0]["IC50 MT"], "12")

    def test_a_whole_pvacseq_directory(self):
        agr.main([str(pvacseq_tree(self.root)), "--outdir", str(self.root / "out")])
        written = sorted(p.name for p in (self.root / "out").glob("*.agretopicity.tsv"))
        self.assertEqual(len(written), 3)

    def test_a_directory_needs_outdir_not_o(self):
        with self.assertRaises(SystemExit):
            agr.main([str(pvacseq_tree(self.root)), "-o", str(self.root / "out.tsv")])

    def test_no_output_argument_at_all(self):
        path = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, self.root / "in.Combined.tsv")
        with self.assertRaises(SystemExit):
            agr.main([str(path)])

    def test_an_input_that_is_not_a_report(self):
        with self.assertRaises(SystemExit):
            agr.main([str(self.root / "missing"), "-o", str(self.root / "out.tsv")])

    def test_thresholds_are_arguments(self):
        path = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, self.root / "in.Combined.tsv")
        out = self.root / "out.tsv"
        agr.main([str(path), "-o", str(out), "--foreign", "0.001"])
        self.assertEqual(read(out)[0]["Agretopicity Class"], "intermediate")

    def test_foreign_above_self_like_is_rejected(self):
        path = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, self.root / "in.Combined.tsv")
        with self.assertRaises(SystemExit):
            agr.main([str(path), "-o", str(self.root / "o.tsv"), "--foreign", "0.9", "--self-like", "0.1"])


class TestFusionInput(unittest.TestCase):
    """A fusion peptide has no wild-type counterpart, so there is nothing
    to divide by. That is different from a missing value."""

    def setUp(self):
        self.root = Path(tempfile.mkdtemp())
        agr.main([str(pvacfuse_tree(self.root)), "--outdir", str(self.root / "out")])
        self.rows = read(next((self.root / "out").glob("*.agretopicity.tsv")))

    def test_it_runs_rather_than_failing(self):
        self.assertEqual(len(self.rows), 1)

    def test_the_class_says_it_does_not_apply(self):
        self.assertEqual(self.rows[0]["Agretopicity Class"], "not_applicable")

    def test_the_numbers_are_na(self):
        self.assertEqual(self.rows[0]["Agretopicity"], "NA")
        self.assertEqual(self.rows[0]["DAI"], "NA")

    def test_the_fusion_columns_survive(self):
        self.assertEqual(self.rows[0]["Read Support"], "30")


class TestTimepoints(unittest.TestCase):
    VCF_HEADER = ("##fileformat=VCFv4.2\n"
                  "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnormal\ttumor\n")

    def vcf(self, path, depth, vaf):
        line = "\t".join(["chr1", "100", ".", "A", "T", ".", "PASS", ".", "GT:RDP:RAF",
                          "0/0:0:0", f"0/1:{depth}:{vaf}"])
        path.write_text(self.VCF_HEADER + line + "\n")
        return path

    def test_both_timepoints_are_joined_on(self):
        root = Path(tempfile.mkdtemp())
        report = write_report([seq_row("chr1-100-100-A-T", "12", "900")],
                              SEQ_COLUMNS, root / "in.Combined.tsv")
        baseline = self.vcf(root / "baseline.vcf", 80, "0.40")
        week6 = self.vcf(root / "week6.vcf", 55, "0.22")
        out = root / "out.tsv"
        agr.main([str(report), "-o", str(out),
                  "--timepoint", f"Baseline={baseline}", "--timepoint", f"Week 6={week6}"])
        row = read(out)[0]
        self.assertEqual(row["RNA Depth (Baseline)"], "80")
        self.assertEqual(row["RNA VAF (Week 6)"], "0.22")
        # and the agretopicity is still there alongside them
        self.assertEqual(row["Agretopicity Class"], "foreign")

    def test_a_malformed_timepoint(self):
        root = Path(tempfile.mkdtemp())
        report = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, root / "in.Combined.tsv")
        with self.assertRaises(SystemExit):
            agr.main([str(report), "-o", str(root / "o.tsv"), "--timepoint", "no-equals-sign"])

    def test_a_timepoint_vcf_that_does_not_exist(self):
        root = Path(tempfile.mkdtemp())
        report = write_report([seq_row("v1", "12", "900")], SEQ_COLUMNS, root / "in.Combined.tsv")
        with self.assertRaises(SystemExit):
            agr.main([str(report), "-o", str(root / "o.tsv"), "--timepoint", f"Week 6={root}/nope.vcf"])


if __name__ == "__main__":
    unittest.main()
