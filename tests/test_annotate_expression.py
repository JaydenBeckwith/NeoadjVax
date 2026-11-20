import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "annotate_expression", Path(__file__).resolve().parents[1] / "bin" / "annotate_expression.py")
annotate = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(annotate)

SKIN = "Skin - Sun Exposed (Lower leg)"
DB_HEADER = "ENSEMBL_ID\tGENE_NAME\tGROUP\tQUANTILE\tTPM\n"
# TP53 broadly expressed; CTAG1B testis only; A1BG nowhere; KRT80 in one tissue
GTEX = DB_HEADER + "".join(
    f"{gene}\t{name}\t{tissue}\t0.5\t{tpm:.4f}\n"
    for gene, name, values in [
        ("ENSG00000141510", "TP53", {SKIN: 12.5, "Testis": 30.1, "Whole Blood": 8.0,
                                     "Liver": 9.0, "Lung": 7.0, "Brain - Cortex": 6.0,
                                     "Pancreas": 5.0, "Thyroid": 4.0, "Ovary": 3.0, "Spleen": 2.0}),
        ("ENSG00000213088", "CTAG1B", {SKIN: 0.0, "Testis": 25.0, "Whole Blood": 0.0,
                                       "Liver": 0.0, "Lung": 0.0, "Brain - Cortex": 0.0,
                                       "Pancreas": 0.0, "Thyroid": 0.0, "Ovary": 0.0, "Spleen": 0.0}),
        ("ENSG00000121410", "A1BG", {SKIN: 0.1, "Testis": 0.2, "Whole Blood": 0.05,
                                     "Liver": 0.0, "Lung": 0.0, "Brain - Cortex": 0.0,
                                     "Pancreas": 0.0, "Thyroid": 0.0, "Ovary": 0.0, "Spleen": 0.0}),
        ("ENSG00000158988", "KRT80", {SKIN: 20.0, "Testis": 0.0, "Whole Blood": 0.0,
                                      "Liver": 0.0, "Lung": 0.0, "Brain - Cortex": 0.0,
                                      "Pancreas": 0.0, "Thyroid": 0.0, "Ovary": 0.0, "Spleen": 0.0}),
    ]
    for tissue, tpm in values.items())

TCGA = DB_HEADER \
    + "ENSG00000141510\tTP53\tTCGA-SKCM\t0.5\t15.0000\n" \
    + "ENSG00000141510\tTP53\tTCGA-LUAD\t0.5\t3.0000\n" \
    + "ENSG00000121410\tA1BG\tTCGA-SKCM\t0.5\t0.2000\n"

SAMPLE_TPM = "ENSEMBL_ID\tGENE_NAME\tGENE_BIOTYPE\tCOUNTS\tLENGTH\tTPM\tEXPRESSION_QUARTILE\n" \
             "ENSG00000141510\tTP53\tprotein_coding\t400\t2000\t120.0000\t4\n" \
             "ENSG00000121410\tA1BG\tprotein_coding\t0\t500\t0.0000\t0\n" \
             "ENSG00000213088\tCTAG1B\tprotein_coding\t2\t900\t0.5000\t1\n"

REPORT = "ID\tGene\tDNA VAF\n" \
         "chr17-7676154-7676154-G-A\tTP53\t0.42\n" \
         "chr6-1-1-C-T\tCTAG1B\t0.30\n" \
         "chr19-1-1-C-T\tA1BG\t0.10\n" \
         "chr12-1-1-C-T\tKRT80\t0.20\n" \
         "chr1-1-1-C-T\tNOSUCHGENE\t0.10\n"


class AnnotateExpression(unittest.TestCase):
    def run_annotate(self, *extra, report=REPORT):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "report.tsv").write_text(report)
            (tmp / "gtex.tsv").write_text(GTEX)
            (tmp / "tcga.tsv").write_text(TCGA)
            (tmp / "sample.tsv").write_text(SAMPLE_TPM)
            args = [str(tmp / "report.tsv")] + [a.replace("TMP", str(tmp)) for a in extra]
            args += ["-o", str(tmp / "out.tsv")]
            annotate.main(args)
            with (tmp / "out.tsv").open() as handle:
                return {row["Gene"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def with_rna(self, *extra):
        return self.run_annotate("--gtex", "TMP/gtex.tsv", "--gtex-tissue", SKIN,
                                 "--tcga", "TMP/tcga.tsv", "--tcga-cancer-type", "TCGA-SKCM",
                                 "--sample-tpm", "TMP/sample.tsv", *extra)

    def without_rna(self, *extra):
        return self.run_annotate("--gtex", "TMP/gtex.tsv", "--gtex-tissue", SKIN,
                                 "--tcga", "TMP/tcga.tsv", "--tcga-cancer-type", "TCGA-SKCM", *extra)

    def test_reports_the_matched_normal_tissue(self):
        self.assertAlmostEqual(float(self.with_rna()["TP53"][f"GTEx {SKIN} TPM"]), 12.5, places=3)

    def test_summarises_expression_across_all_tissues(self):
        row = self.with_rna()["TP53"]
        self.assertAlmostEqual(float(row["GTEx Max TPM"]), 30.1, places=3)
        self.assertEqual(row["GTEx Tissues Expressed"], "10")

    def test_calls_a_testis_only_gene_testis_restricted(self):
        # the class of target a melanoma vaccine wants: silent in every
        # somatic tissue, so the testis columns must not be averaged away
        self.assertEqual(self.with_rna()["CTAG1B"]["GTEx Specificity"], "testis_restricted")

    def test_calls_a_silent_gene_not_expressed(self):
        self.assertEqual(self.with_rna()["A1BG"]["GTEx Specificity"], "not_expressed")

    def test_calls_a_widely_expressed_gene_broadly_expressed(self):
        self.assertEqual(self.with_rna()["TP53"]["GTEx Specificity"], "broadly_expressed")

    def test_calls_a_single_tissue_gene_tissue_restricted(self):
        self.assertEqual(self.with_rna()["KRT80"]["GTEx Specificity"], "tissue_restricted")

    def test_reports_the_matched_tcga_cohort_only(self):
        self.assertAlmostEqual(float(self.with_rna()["TP53"]["TCGA TCGA-SKCM TPM"]), 15.0, places=3)

    def test_sample_rna_is_the_expression_evidence_when_present(self):
        row = self.with_rna()["TP53"]
        self.assertEqual(row["Expression Evidence"], "sample_rna")
        self.assertAlmostEqual(float(row["Sample TPM"]), 120.0, places=3)
        self.assertEqual(row["Expression Flag"], "OK")

    def test_zero_sample_tpm_flags_no_expression(self):
        self.assertEqual(self.with_rna()["A1BG"]["Expression Flag"], "NO_EXPRESSION")

    def test_low_sample_tpm_flags_low_expression(self):
        self.assertEqual(self.with_rna()["CTAG1B"]["Expression Flag"], "LOW_EXPRESSION")

    def test_min_tpm_threshold_moves_the_flag(self):
        self.assertEqual(self.with_rna("--min-tpm", "200")["TP53"]["Expression Flag"], "LOW_EXPRESSION")

    def test_falls_back_to_the_cohort_without_rna(self):
        row = self.without_rna()["TP53"]
        self.assertEqual(row["Sample TPM"], "NA")
        self.assertEqual(row["Expression Evidence"], "tcga_cohort")
        self.assertEqual(row["Expression Flag"], "OK")

    def test_cohort_fallback_flags_a_gene_the_cohort_barely_expresses(self):
        self.assertEqual(self.without_rna()["A1BG"]["Expression Flag"], "LOW_EXPRESSION")

    def test_unknown_gene_gets_no_data_rather_than_a_zero(self):
        row = self.with_rna()["NOSUCHGENE"]
        self.assertEqual(row[f"GTEx {SKIN} TPM"], "NA")
        self.assertEqual(row["Expression Flag"], "NO_DATA")

    def test_original_columns_are_preserved(self):
        self.assertEqual(self.with_rna()["TP53"]["DNA VAF"], "0.42")

    def test_candidate_with_several_gene_symbols_resolves_the_known_one(self):
        report = "ID\tGene\tDNA VAF\nchr1-1-1-C-T\tNOSUCHGENE,TP53\t0.3\n"
        rows = self.run_annotate("--gtex", "TMP/gtex.tsv", "--gtex-tissue", SKIN, report=report)
        self.assertAlmostEqual(float(list(rows.values())[0][f"GTEx {SKIN} TPM"]), 12.5, places=3)

    def test_matches_on_ensembl_id_when_that_is_what_the_report_carries(self):
        report = "ID\tGene\tDNA VAF\nchr17-1-1-G-A\tENSG00000141510\t0.3\n"
        rows = self.run_annotate("--gtex", "TMP/gtex.tsv", "--gtex-tissue", SKIN, report=report)
        self.assertAlmostEqual(float(rows["ENSG00000141510"][f"GTEx {SKIN} TPM"]), 12.5, places=3)

    def test_unknown_tissue_name_fails_loudly(self):
        with self.assertRaises(SystemExit):
            self.run_annotate("--gtex", "TMP/gtex.tsv", "--gtex-tissue", "Pineal Gland")

    def test_unknown_cohort_name_fails_loudly(self):
        with self.assertRaises(SystemExit):
            self.run_annotate("--tcga", "TMP/tcga.tsv", "--tcga-cancer-type", "TCGA-NOPE")

    def test_table_needs_its_group_name(self):
        with self.assertRaises(SystemExit):
            self.run_annotate("--gtex", "TMP/gtex.tsv")

    def test_runs_with_neither_table(self):
        row = self.run_annotate("--sample-tpm", "TMP/sample.tsv")["TP53"]
        self.assertEqual(row["Expression Flag"], "OK")
        self.assertNotIn("GTEx Median TPM", row)


if __name__ == "__main__":
    unittest.main()
