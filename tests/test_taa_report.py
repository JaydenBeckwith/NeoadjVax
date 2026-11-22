"""taa_report.py: shared tumour antigens as an expression question."""
import csv
import importlib.util
import sys
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin"
sys.path.insert(0, str(BIN))
spec = importlib.util.spec_from_file_location("taa_report", BIN / "taa_report.py")
taa = importlib.util.module_from_spec(spec)
spec.loader.exec_module(taa)

PANEL = """# a comment line
MAGEA3\tcancer_testis\tMAGE-A3

TYR\tmelanocyte_differentiation\tTyrosinase
MUC1\toverexpressed\tMUC1
"""

# build_expression_db.py long format.
DB_HEADER = ["ENSEMBL_ID", "GENE_NAME", "GROUP", "QUANTILE", "TPM"]
GTEX_ROWS = [
    # MAGEA3: testis only
    ["ENSG1", "MAGEA3", "Testis", "0.5", "30"],
    ["ENSG1", "MAGEA3", "Skin - Sun Exposed (Lower leg)", "0.5", "0.01"],
    ["ENSG1", "MAGEA3", "Lung", "0.5", "0.0"],
    # TYR: skin
    ["ENSG2", "TYR", "Skin - Sun Exposed (Lower leg)", "0.5", "12"],
    ["ENSG2", "TYR", "Testis", "0.5", "0.0"],
    ["ENSG2", "TYR", "Lung", "0.5", "0.0"],
    # MUC1: everywhere
    ["ENSG3", "MUC1", "Skin - Sun Exposed (Lower leg)", "0.5", "40"],
    ["ENSG3", "MUC1", "Lung", "0.5", "120"],
    ["ENSG3", "MUC1", "Testis", "0.5", "35"],
    ["ENSG3", "MUC1", "Pancreas", "0.5", "88"],
    ["ENSG3", "MUC1", "Thyroid", "0.5", "70"],
]
TCGA_ROWS = [
    ["ENSG1", "MAGEA3", "TCGA-SKCM", "0.5", "8"],
    ["ENSG2", "TYR", "TCGA-SKCM", "0.5", "150"],
    ["ENSG3", "MUC1", "TCGA-SKCM", "0.5", "30"],
]


def write(directory, name, header, rows):
    path = Path(directory) / name
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(header)
        writer.writerows(rows)
    return path


def run(sample_tpm=None, with_gtex=True, with_tcga=True, **kwargs):
    directory = Path(tempfile.mkdtemp())
    panel = directory / "panel.tsv"
    panel.write_text(PANEL)
    out = directory / "taa.tsv"
    argv = ["--panel", str(panel), "-o", str(out)]
    if sample_tpm is not None:
        rows = [[f"ENSG{i}", gene, str(tpm)] for i, (gene, tpm) in enumerate(sample_tpm.items(), start=1)]
        argv += ["--sample-tpm", str(write(directory, "tpm.tsv", ["ENSEMBL_ID", "GENE_NAME", "TPM"], rows))]
    if with_gtex:
        argv += ["--gtex", str(write(directory, "gtex.tsv", DB_HEADER, GTEX_ROWS)),
                 "--gtex-tissue", "Skin - Sun Exposed (Lower leg)"]
    if with_tcga:
        argv += ["--tcga", str(write(directory, "tcga.tsv", DB_HEADER, TCGA_ROWS)),
                 "--tcga-cancer-type", "TCGA-SKCM"]
    for key, value in kwargs.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    taa.main(argv)
    with out.open(newline="") as handle:
        return {row["Gene"]: row for row in csv.DictReader(handle, delimiter="\t")}


class TestPanelParsing(unittest.TestCase):
    def test_comments_and_blank_lines_are_skipped(self):
        path = Path(tempfile.mkdtemp()) / "p.tsv"
        path.write_text(PANEL)
        rows = taa.read_panel(path)
        self.assertEqual([r[0] for r in rows], ["MAGEA3", "TYR", "MUC1"])

    def test_class_and_note_are_kept(self):
        path = Path(tempfile.mkdtemp()) / "p.tsv"
        path.write_text(PANEL)
        self.assertEqual(taa.read_panel(path)[0], ("MAGEA3", "cancer_testis", "MAGE-A3"))

    def test_gene_only_rows_are_allowed(self):
        path = Path(tempfile.mkdtemp()) / "p.tsv"
        path.write_text("MAGEA1\n")
        self.assertEqual(taa.read_panel(path), [("MAGEA1", "", "")])

    def test_an_empty_panel_is_an_error(self):
        path = Path(tempfile.mkdtemp()) / "p.tsv"
        path.write_text("# nothing\n")
        with self.assertRaises(SystemExit):
            taa.read_panel(path)

    def test_the_shipped_panel_parses(self):
        rows = taa.read_panel(BIN.parent / "assets" / "taa_panel.tsv")
        genes = [r[0] for r in rows]
        self.assertIn("CTAG1B", genes)
        self.assertIn("PMEL", genes)
        self.assertEqual(len(genes), len(set(genes)), "a gene is listed twice")
        self.assertTrue(all(r[1] for r in rows), "every gene needs a class")


class TestSpecificity(unittest.TestCase):
    def test_testis_only_is_privileged(self):
        self.assertEqual(taa.specificity({"Testis": 30, "Lung": 0.0}, 1.0), "testis_restricted")

    def test_nothing_above_the_threshold(self):
        self.assertEqual(taa.specificity({"Lung": 0.2}, 1.0), "not_expressed")

    def test_a_few_tissues(self):
        self.assertEqual(taa.specificity({"Skin": 10, "Lung": 5}, 1.0), "tissue_restricted")

    def test_many_tissues(self):
        values = {f"T{i}": 10 for i in range(6)}
        self.assertEqual(taa.specificity(values, 1.0), "broadly_expressed")

    def test_no_reference_data(self):
        self.assertEqual(taa.specificity({}, 1.0), "NA")


class TestCalls(unittest.TestCase):
    def test_expressed_and_testis_restricted_is_a_good_target(self):
        row = run(sample_tpm={"MAGEA3": 25, "TYR": 200, "MUC1": 50})["MAGEA3"]
        self.assertEqual(row["Call"], "good_target")
        self.assertEqual(row["GTEx Specificity"], "testis_restricted")
        self.assertEqual(row["Evidence"], "sample_rna")

    def test_expressed_in_normal_skin_carries_risk(self):
        row = run(sample_tpm={"MAGEA3": 25, "TYR": 200, "MUC1": 50})["TYR"]
        self.assertEqual(row["Call"], "target_with_risk")

    def test_broadly_expressed_carries_risk(self):
        row = run(sample_tpm={"MAGEA3": 25, "TYR": 200, "MUC1": 50})["MUC1"]
        self.assertEqual(row["GTEx Specificity"], "broadly_expressed")
        self.assertEqual(row["Call"], "target_with_risk")

    def test_not_expressed_in_this_tumour(self):
        row = run(sample_tpm={"MAGEA3": 0.0, "TYR": 200, "MUC1": 50})["MAGEA3"]
        self.assertEqual(row["Call"], "not_expressed")
        self.assertEqual(row["Sample Expressed"], "no")

    def test_falls_back_to_the_cohort_without_patient_rna(self):
        row = run()["MAGEA3"]
        self.assertEqual(row["Evidence"], "tcga_cohort")
        self.assertEqual(row["Call"], "good_target")

    def test_no_rna_and_no_cohort_is_no_data(self):
        row = run(with_tcga=False)["MAGEA3"]
        self.assertEqual(row["Call"], "no_data")

    def test_no_gtex_means_risk_is_unknown_not_absent(self):
        row = run(sample_tpm={"MAGEA3": 25}, with_gtex=False)["MAGEA3"]
        self.assertEqual(row["GTEx Specificity"], "NA")
        self.assertEqual(row["Call"], "target_with_risk")

    def test_thresholds_are_arguments(self):
        row = run(sample_tpm={"MAGEA3": 5}, min_tpm=10)["MAGEA3"]
        self.assertEqual(row["Call"], "not_expressed")


class TestOutput(unittest.TestCase):
    def test_columns(self):
        row = run(sample_tpm={"MAGEA3": 25})["MAGEA3"]
        for column in taa.COLUMNS:
            self.assertIn(column, row)

    def test_reported_tissue_and_max_tissue(self):
        row = run(sample_tpm={"MAGEA3": 25})["MAGEA3"]
        self.assertEqual(row["GTEx Max Tissue"], "Testis")
        self.assertEqual(float(row["GTEx Tissue TPM"]), 0.01)

    def test_good_targets_sort_first(self):
        directory = Path(tempfile.mkdtemp())
        panel = directory / "panel.tsv"
        panel.write_text(PANEL)
        out = directory / "taa.tsv"
        tpm = write(directory, "tpm.tsv", ["ENSEMBL_ID", "GENE_NAME", "TPM"],
                    [["ENSG1", "MAGEA3", "25"], ["ENSG2", "TYR", "200"], ["ENSG3", "MUC1", "50"]])
        taa.main(["--panel", str(panel), "--sample-tpm", str(tpm), "-o", str(out),
                  "--gtex", str(write(directory, "g.tsv", DB_HEADER, GTEX_ROWS)),
                  "--gtex-tissue", "Skin - Sun Exposed (Lower leg)"])
        with out.open(newline="") as handle:
            genes = [row["Gene"] for row in csv.DictReader(handle, delimiter="\t")]
        self.assertEqual(genes[0], "MAGEA3")

    def test_an_unknown_tissue_name_is_an_error(self):
        with self.assertRaises(SystemExit):
            run(gtex_tissue="Not A Tissue")

    def test_gtex_without_a_tissue_is_an_error(self):
        directory = Path(tempfile.mkdtemp())
        panel = directory / "panel.tsv"
        panel.write_text(PANEL)
        with self.assertRaises(SystemExit):
            taa.main(["--panel", str(panel), "-o", str(directory / "o.tsv"),
                      "--gtex", str(write(directory, "g.tsv", DB_HEADER, GTEX_ROWS))])


if __name__ == "__main__":
    unittest.main()
