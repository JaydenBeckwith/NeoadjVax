"""Tests for evidence parsing and BAM/reference compatibility."""
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "run_ervcaller", Path(__file__).resolve().parents[1] / "bin" / "run_ervcaller.py")
runner = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(runner)
HEADER = "##fileformat=VCFv4.2\n#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tP1.tumor\n"
RECORD = (
    "chr1\t5617379\t.\tT\t<INS_MEI:HERV>\t.\t.\t"
    "TSD=NULL,NULL;INFOR=HERVK,1,7831,7831,+,4;CR=64;SR=3;GTF=YES;GR=1.000\t"
    "GT:GQ:GL:DPN:DPI\t1/1:40:0,0,1:0:67\n"
)


class EvidenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.path = Path(self.temp.name) / "sample.vcf"

    def parse(self, content, **kwargs):
        self.path.write_text(content)
        return runner.read_vcf(self.path, "P1", "tumor", **kwargs)

    def test_keeps_evidence_and_raw_coordinate(self):
        row, = self.parse(HEADER + RECORD)
        self.assertEqual(row["caller_pos"], "5617379")
        self.assertEqual(row["te_name"], "HERVK")
        self.assertEqual(row["chimeric_reads"], "64")
        self.assertEqual(row["DPI"], "67")
        self.assertEqual(row["DPN"], "0")
        self.assertEqual(set(row), set(runner.FIELDS))
        self.assertNotIn("somatic", row)

    def test_format_order_is_not_assumed(self):
        row, = self.parse(HEADER + RECORD.replace(
            "GT:GQ:GL:DPN:DPI\t1/1:40:0,0,1:0:67",
            "DPI:GL:DPN:GQ:GT\t67:0,0,1:0:40:1/1"))
        self.assertEqual((row["GT"], row["DPI"], row["DPN"]), ("1/1", "67", "0"))

    def test_info_order_is_not_assumed(self):
        row, = self.parse(HEADER + RECORD.replace(
            "CR=64;SR=3;GTF=YES;GR=1.000", "GR=1.000;GTF=YES;SR=3;CR=64"))
        self.assertEqual(row["split_reads"], "3")

    def test_keeps_low_quality_calls(self):
        row, = self.parse(HEADER + RECORD.replace("1/1:40:", "0/1:0:"))
        self.assertEqual(row["GQ"], "0")

    def test_ungenotyped_calls_are_not_dropped(self):
        row, = self.parse(HEADER + RECORD.replace("GTF=YES;GR=1.000", "GTF=NO;GR=NULL")
                         .replace("1/1:40:0,0,1:0:67", "./.:.:.,.,.:.:67"))
        self.assertEqual((row["GT"], row["genotyped"]), ("./.", "NO"))

    def test_header_only_is_valid_zero_calls(self):
        self.assertEqual(self.parse(HEADER), [])

    def test_blank_output_is_not_silently_zero_calls(self):
        with self.assertRaisesRegex(ValueError, "Missing VCF header"):
            self.parse("")

    def test_missing_column_header(self):
        with self.assertRaises(ValueError):
            self.parse("##fileformat=VCFv4.2\n" + RECORD)

    def test_tumour_normal_mismatch_rejected(self):
        with self.assertRaisesRegex(ValueError, "sample mismatch"):
            self.parse(HEADER + RECORD, expected_sample="P1.normal")

    def test_multiple_samples_rejected(self):
        with self.assertRaisesRegex(ValueError, "single-sample"):
            self.parse(HEADER.replace("P1.tumor\n", "P1.tumor\tP1.normal\n"))

    def test_invalid_record_fields(self):
        for content in (
            RECORD.replace("SR=3", "SR=NA"),
            RECORD.replace("INFOR=HERVK,1,7831,7831,+,4", "INFOR=HERVK"),
            RECORD.replace("5617379", "-1"),
            RECORD.replace("1/1:40:0,0,1:0:67", "1/1:40"),
        ):
            with self.subTest(content=content), self.assertRaises(ValueError):
                self.parse(HEADER + content)

    def test_source_zero_coordinate_is_retained_not_normalized(self):
        row, = self.parse(HEADER + RECORD.replace("5617379", "0"))
        self.assertEqual(row["caller_pos"], "0")


class ReferenceTests(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.fai = Path(self.temp.name) / "genome.fa.fai"
        self.fai.write_text("chr1\t1000\t6\t60\t61\nchr2\t2000\t0\t60\t61\n")

    def test_matching_reference(self):
        self.assertEqual(runner.check_reference_bam(
            "@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:1000\n", self.fai), "chr1")

    def test_wrong_build_or_contig_name_rejected(self):
        for sq in ("@SQ\tSN:chr1\tLN:999", "@SQ\tSN:1\tLN:1000"):
            with self.subTest(sq=sq), self.assertRaisesRegex(ValueError, "mismatch"):
                runner.check_reference_bam("@HD\tSO:coordinate\n" + sq, self.fai)

    def test_unsorted_bam_rejected(self):
        with self.assertRaisesRegex(ValueError, "coordinate-sorted"):
            runner.check_reference_bam("@HD\tSO:queryname\n@SQ\tSN:chr1\tLN:1000\n", self.fai)

    def test_missing_dictionary_rejected(self):
        with self.assertRaises(ValueError):
            runner.check_reference_bam("@HD\tSO:coordinate\n", self.fai)


if __name__ == "__main__":
    unittest.main()
