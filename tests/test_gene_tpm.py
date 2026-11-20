import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("gene_tpm", Path(__file__).resolve().parents[1] / "bin" / "gene_tpm.py")
gene_tpm = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(gene_tpm)

# Lengths chosen so the expected TPMs are exact: RPKs are 200, 100, 0, 10.
COUNTS = "# Program:featureCounts\n" \
         "Geneid\tChr\tStart\tEnd\tStrand\tLength\tsample.bam\n" \
         "ENSG00000141510.16\tchr17\t1\t100\t-\t2000\t400\n" \
         "ENSG00000157764.13\tchr7\t1\t100\t-\t1000\t100\n" \
         "ENSG00000121410.11\tchr19\t1\t100\t+\t500\t0\n" \
         "ENSG00000012048.22\tchr17\t1\t100\t-\t4000\t40\n"

GTF = 'chr17\tx\tgene\t1\t100\t.\t-\t.\tgene_id "ENSG00000141510.16"; gene_type "protein_coding"; gene_name "TP53";\n' \
      'chr7\tx\tgene\t1\t100\t.\t-\t.\tgene_id "ENSG00000157764.13"; gene_type "lncRNA"; gene_name "BRAF";\n'


class GeneTpm(unittest.TestCase):
    def run_tpm(self, counts=COUNTS, gtf=None):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "counts.txt").write_text(counts)
            args = [str(tmp / "counts.txt"), "-o", str(tmp / "out.tsv")]
            if gtf is not None:
                (tmp / "genes.gtf").write_text(gtf)
                args += ["--gtf", str(tmp / "genes.gtf")]
            gene_tpm.main(args)
            with (tmp / "out.tsv").open() as handle:
                return {row["ENSEMBL_ID"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_tpm_sums_to_a_million(self):
        rows = self.run_tpm()
        self.assertAlmostEqual(sum(float(row["TPM"]) for row in rows.values()), 1e6, places=1)

    def test_tpm_is_length_normalised(self):
        rows = self.run_tpm()
        # 400 reads over 2 kb against 100 reads over 1 kb: twice the RPK
        self.assertAlmostEqual(float(rows["ENSG00000141510"]["TPM"]),
                               2 * float(rows["ENSG00000157764"]["TPM"]), places=1)

    def test_version_suffix_stripped_to_match_reference_tables(self):
        self.assertIn("ENSG00000012048", self.run_tpm())

    def test_unexpressed_gene_is_quartile_zero(self):
        self.assertEqual(self.run_tpm()["ENSG00000121410"]["EXPRESSION_QUARTILE"], "0")

    def test_quartiles_rank_expressed_genes(self):
        rows = self.run_tpm()
        self.assertGreater(int(rows["ENSG00000141510"]["EXPRESSION_QUARTILE"]),
                           int(rows["ENSG00000012048"]["EXPRESSION_QUARTILE"]))

    def test_gtf_supplies_names_and_biotypes(self):
        rows = self.run_tpm(gtf=GTF)
        self.assertEqual(rows["ENSG00000141510"]["GENE_NAME"], "TP53")
        self.assertEqual(rows["ENSG00000157764"]["GENE_BIOTYPE"], "lncRNA")
        # a gene the GTF does not mention keeps its ID as the name
        self.assertEqual(rows["ENSG00000012048"]["GENE_NAME"], "ENSG00000012048")

    def test_several_bam_columns_are_summed(self):
        counts = "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\tb.bam\n" \
                 "ENSG00000141510.16\tchr17\t1\t100\t-\t1000\t300\t100\n"
        self.assertEqual(self.run_tpm(counts)["ENSG00000141510"]["COUNTS"], "400")

    def test_rejects_a_table_that_is_not_featurecounts(self):
        with self.assertRaises(SystemExit):
            self.run_tpm("gene\tcount\nTP53\t5\n")

    def test_rejects_an_all_zero_table(self):
        counts = "Geneid\tChr\tStart\tEnd\tStrand\tLength\ta.bam\n" \
                 "ENSG00000141510.16\tchr17\t1\t100\t-\t1000\t0\n"
        with self.assertRaises(SystemExit):
            self.run_tpm(counts)


if __name__ == "__main__":
    unittest.main()
