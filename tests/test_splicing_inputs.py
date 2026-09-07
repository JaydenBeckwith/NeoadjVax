import csv
import gzip
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("splicing", Path(__file__).resolve().parents[1] / "bin/prepare_splicing_inputs.py")
splicing = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(splicing)
HEADER = ('##fileformat=VCFv4.2\n##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI">\n'
          '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n')


class SplicingInputs(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def vcf(self, records, threshold=0.5, header=HEADER, compressed=False):
        source = self.root / ("input.vcf.gz" if compressed else "input.vcf")
        with (gzip.open(source, "wt") if compressed else source.open("w")) as out:
            out.write(header + records)
        output = self.root / "out.tsv"
        qc = splicing.spliceai_table(source, output, self.root / "qc.json", threshold)
        with output.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t")), qc

    def counts(self, text, sample="", min_reads=3):
        source = self.root / "counts.txt"
        source.write_text(text)
        evidence = self.root / "evidence.tsv"
        splicing.leafcutter_sample(source, sample, self.root / "selected.gz", evidence, min_reads)
        with evidence.open() as handle:
            return list(csv.DictReader(handle, delimiter="\t"))

    def test_standard_info_and_allele_matching(self):
        rows, qc = self.vcf("chr1\t100\t.\tA\tC,G\t.\tPASS\tSpliceAI=C|GENE|0.1|0|0|0|4|0|0|0,G|GENE|0.9|0.1|0|0|5|-1|0|0\n", compressed=True)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["ALT"], "G")
        self.assertEqual(rows[0]["DS_AL"], "0.1")  # R filters each effect, not the entire row.
        self.assertEqual(qc["below_threshold"], 1)

    def test_pass_only(self):
        rows, qc = self.vcf("chr1\t100\t.\tA\tG\t.\t.\tSpliceAI=G|GENE|0.9|0|0|0|5|0|0|0\n")
        self.assertFalse(rows)
        self.assertEqual(qc["non_pass"], 1)

    def test_missing_annotations_counted(self):
        rows, qc = self.vcf("chr1\t100\t.\tA\tG\t.\tPASS\t.\n")
        self.assertFalse(rows)
        self.assertEqual(qc["missing_annotation"], 1)

    def test_header_only_vcf(self):
        self.assertEqual(self.vcf("")[0], [])

    def test_no_spliceai_header_is_error(self):
        with self.assertRaisesRegex(ValueError, "headers"):
            self.vcf("", header=HEADER.replace('##INFO=<ID=SpliceAI,Number=.,Type=String,Description="SpliceAI">\n', ""))

    def test_invalid_score_or_offset(self):
        for score, offset in [("nan", "0"), ("1.1", "0"), ("-0.1", "0"), ("0.8", ".")]:
            with self.subTest(score=score, offset=offset), self.assertRaises(ValueError):
                self.vcf(f"chr1\t100\t.\tA\tG\t.\tPASS\tSpliceAI=G|GENE|{score}|0|0|0|{offset}|0|0|0\n")

    def test_missing_scores(self):
        rows, qc = self.vcf("chr1\t100\t.\tA\tG\t.\tPASS\tSpliceAI=G|GENE|.|.|.|.|.|.|.|.\n")
        self.assertFalse(rows)
        self.assertEqual(qc["missing_scores"], 1)

    def test_malformed_annotation_or_wrong_allele(self):
        for annotation in ["G|GENE|0.9", "T|GENE|0.9|0|0|0|5|0|0|0"]:
            with self.subTest(annotation=annotation), self.assertRaises(ValueError):
                self.vcf(f"chr1\t100\t.\tA\tG\t.\tPASS\tSpliceAI={annotation}\n")

    def test_threshold_zero_not_replaced_by_default(self):
        rows, _ = self.vcf("chr1\t100\t.\tA\tG\t.\tPASS\tSpliceAI=G|GENE|0.01|0|0|0|5|0|0|0\n", threshold=0)
        self.assertEqual(len(rows), 1)

    def test_cohort_selects_only_this_timepoint(self):
        data = "chrom P1_day0 P1_day99\nchr1:160:295:clu_1_+ 0/50 5/50\nchr1:160:301:clu_1_+ 50/50 45/50\n"
        self.assertEqual(len(self.counts(data, "P1_day0")), 1)
        rows = self.counts(data, "P1_day99")
        self.assertEqual(rows[0]["junc_id"], "chr1:160-295:+")
        self.assertEqual(rows[0]["rna_junction_reads"], "5")

    def test_negative_strand_no_extra_coordinate_offset(self):
        row, = self.counts("chrom S\nchr1:966:1101:clu_2_- 3/50\n")
        self.assertEqual(row["junc_id"], "chr1:966-1101:-")

    def test_ambiguous_or_missing_sample(self):
        for sample in ["", "missing"]:
            with self.subTest(sample=sample), self.assertRaises(ValueError):
                self.counts("chrom A B\nchr1:10:20:clu_1_+ 3/5 2/5\n", sample)

    def test_zero_low_and_empty_counts(self):
        self.assertFalse(self.counts("chrom A\nchr1:10:20:clu_1_+ 0/0\nchr1:20:30:clu_1_+ 2/5\n"))
        self.assertFalse(self.counts("chrom A\n"))

    def test_legacy_or_numerator_only_counts_rejected(self):
        for text in ["chrom S\nchr1:10:20:clu_1 3/4\n", "chrom S\nchr1:10:20:clu_1_+ 3\n"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.counts(text)

    def test_bad_counts_and_duplicates(self):
        for text in ["chrom S\nchr1:10:20:clu_1_+ 5/4\n",
                     "chrom S\nchr1:10:20:clu_1_+ 3/5\nchr1:10:20:clu_2_+ 3/5\n",
                     "chrom S\nchr1:20:10:clu_1_+ 3/5\n"]:
            with self.subTest(text=text), self.assertRaises(ValueError):
                self.counts(text)


if __name__ == "__main__":
    unittest.main()
