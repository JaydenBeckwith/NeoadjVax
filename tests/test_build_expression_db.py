import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "build_expression_db", Path(__file__).resolve().parents[1] / "bin" / "build_expression_db.py")
build = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(build)

GCT = "#1.2\n3\t3\n" \
      "Name\tDescription\tSkin - Sun Exposed (Lower leg)\tTestis\tWhole Blood\n" \
      "ENSG00000141510.16\tTP53\t12.5\t30.1\t8.0\n" \
      "ENSG00000213088.9\tCTAG1B\t0.0\t25.0\t0.0\n" \
      "ENSG00000121410.11\tA1BG\t0.1\tNA\t0.05\n"

# log2(TPM+1) values: 4 -> 15, 5 -> 31, 2 -> 3
MATRIX = "sample\tTCGA-A1\tTCGA-A2\tTCGA-B1\tGTEX-1\n" \
         "ENSG00000141510.16\t4.0\t5.0\t2.0\t9.9\n"
GROUPS = "TCGA-A1\tTCGA-SKCM\nTCGA-A2\tTCGA-SKCM\nTCGA-B1\tTCGA-LUAD\n"
GTF = 'chr17\tx\tgene\t1\t100\t.\t-\t.\tgene_id "ENSG00000141510.16"; gene_name "TP53";\n'


class BuildExpressionDb(unittest.TestCase):
    def build(self, *extra, files=None):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            for name, text in (files or {}).items():
                (tmp / name).write_text(text)
            args = [a.replace("TMP", str(tmp)) for a in extra] + ["-o", str(tmp / "out.tsv")]
            build.main(args)
            with (tmp / "out.tsv").open() as handle:
                return list(csv.DictReader(handle, delimiter="\t"))

    def gct_rows(self):
        return self.build("--from-gct", "TMP/gtex.gct", files={"gtex.gct": GCT})

    def test_gct_gives_one_row_per_gene_and_tissue(self):
        rows = self.gct_rows()
        # 3 genes x 3 tissues, less the one NA cell
        self.assertEqual(len(rows), 8)

    def test_gct_keeps_tissue_values_and_symbols(self):
        row = next(r for r in self.gct_rows() if r["ENSEMBL_ID"] == "ENSG00000141510" and r["GROUP"] == "Testis")
        self.assertEqual(row["GENE_NAME"], "TP53")
        self.assertEqual(float(row["TPM"]), 30.1)
        self.assertEqual(row["QUANTILE"], "0.5")

    def test_gct_strips_gene_version(self):
        self.assertTrue(all(not r["ENSEMBL_ID"].count(".") for r in self.gct_rows()))

    def test_rejects_a_file_that_is_not_a_gct(self):
        with self.assertRaises(SystemExit):
            self.build("--from-gct", "TMP/x.gct", files={"x.gct": "gene\tvalue\nTP53\t1\n"})

    def matrix_rows(self, *extra):
        return self.build("--from-matrix", "TMP/m.tsv", "--groups", "TMP/groups.tsv", *extra,
                          files={"m.tsv": MATRIX, "groups.tsv": GROUPS, "genes.gtf": GTF})

    def test_matrix_groups_samples_by_cohort(self):
        groups = {row["GROUP"] for row in self.matrix_rows("--expm1")}
        self.assertEqual(groups, {"TCGA-SKCM", "TCGA-LUAD"})

    def test_samples_absent_from_the_groups_table_are_ignored(self):
        # GTEX-1 is the highest value in the matrix; it must not reach any group
        self.assertTrue(all(float(row["TPM"]) < 900 for row in self.matrix_rows("--expm1")))

    def test_expm1_undoes_the_log_transform(self):
        row = next(r for r in self.matrix_rows("--expm1") if r["GROUP"] == "TCGA-LUAD")
        self.assertAlmostEqual(float(row["TPM"]), 3.0, places=3)

    def test_without_expm1_values_pass_through(self):
        row = next(r for r in self.matrix_rows() if r["GROUP"] == "TCGA-LUAD")
        self.assertAlmostEqual(float(row["TPM"]), 2.0, places=3)

    def test_requested_quantiles_are_all_stored(self):
        rows = [r for r in self.matrix_rows("--expm1", "--quantiles", "0.25", "0.5", "0.75")
                if r["GROUP"] == "TCGA-SKCM"]
        self.assertEqual({r["QUANTILE"] for r in rows}, {"0.25", "0.5", "0.75"})

    def test_gene_names_fills_symbols_a_matrix_lacks(self):
        self.assertEqual(self.matrix_rows("--expm1")[0]["GENE_NAME"], "NA")
        named = self.matrix_rows("--expm1", "--gene-names", "TMP/genes.gtf")
        self.assertEqual(named[0]["GENE_NAME"], "TP53")

    def test_rejects_groups_that_match_no_matrix_column(self):
        with self.assertRaises(SystemExit):
            self.build("--from-matrix", "TMP/m.tsv", "--groups", "TMP/groups.tsv",
                       files={"m.tsv": MATRIX, "groups.tsv": "OTHER-1\tCOHORT\n"})

    def test_matrix_needs_a_groups_table(self):
        with self.assertRaises(SystemExit):
            self.build("--from-matrix", "TMP/m.tsv", files={"m.tsv": MATRIX})


if __name__ == "__main__":
    unittest.main()
