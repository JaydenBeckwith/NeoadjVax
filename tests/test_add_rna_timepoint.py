import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("add_rna_timepoint", Path(__file__).resolve().parents[1] / "bin" / "add_rna_timepoint.py")
rna = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rna)

VCF = "##fileformat=VCFv4.2\n" \
      "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnormal\ttumor\n" \
      "chr1\t500\t.\tA\tT\t.\tPASS\t.\tGT:AD:RAD:RDP:RAF\t0/0:30,0\t0/1:20,10:40,12:52:0.231\n" \
      "chr1\t900\t.\tCAG\tC\t.\tPASS\t.\tGT:AD:RAD:RDP:RAF\t0/0:30,0\t0/1:20,10:8,0:8:0\n" \
      "chr1\t1200\t.\tG\tGTT\t.\tPASS\t.\tGT:AD\t0/0:30,0\t0/1:20,10\n"
REPORT = "ID\tGene\n" \
         "chr1-499-500-A-T\tSNV\n" \
         "chr1-900-902-CAG-C\tDEL\n" \
         "chr1-1200-1200-G-GTT\tINS_NO_RNA\n" \
         "chr2-99-100-T-C\tABSENT\n"


class AddRnaTimepoint(unittest.TestCase):
    def run_script(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "report.tsv").write_text(REPORT)
            (tmp / "week6.vcf").write_text(VCF)
            rna.main([str(tmp / "report.tsv"), "--vcf", str(tmp / "week6.vcf"), "--label", "week 6", "-o", str(tmp / "out.tsv")])
            with (tmp / "out.tsv").open() as handle:
                return {row["Gene"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_snv_counts_joined_by_pvacseq_id(self):
        row = self.run_script()["SNV"]
        self.assertEqual(row["RNA Depth (week 6)"], "52")
        self.assertEqual(row["RNA VAF (week 6)"], "0.231")

    def test_deletion_starts_at_vcf_pos(self):
        row = self.run_script()["DEL"]
        self.assertEqual((row["RNA Depth (week 6)"], row["RNA VAF (week 6)"]), ("8", "0"))

    def test_missing_counts_and_variants_are_na(self):
        rows = self.run_script()
        for gene in ("INS_NO_RNA", "ABSENT"):
            self.assertEqual((rows[gene]["RNA Depth (week 6)"], rows[gene]["RNA VAF (week 6)"]), ("NA", "NA"))


if __name__ == "__main__":
    unittest.main()
