import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location("annotate_ccf", Path(__file__).resolve().parents[1] / "bin" / "annotate_ccf.py")
ccf = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(ccf)

SEGMENTS = '"chromosome"\t"start.pos"\t"end.pos"\t"Bf"\t"CNt"\t"A"\t"B"\n' \
           '"chr1"\t1\t1000\t0.5\t2\t1\t1\n' \
           '"chr1"\t2000\t5000\t0.33\t4\t3\t1\n'
REPORT = "ID\tGene\tDNA VAF\n" \
         "chr1-499-500-A-T\tG1\t0.25\n" \
         "chr1-2999-3000-C-G\tG2\t0.6\n" \
         "chr1-1499-1500-G-A\tG3\t0.1\n" \
         "chr2-99-100-T-C\tG4\tNA\n"


class AnnotateCcf(unittest.TestCase):
    def run_ccf(self, *extra):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "report.tsv").write_text(REPORT)
            (tmp / "segments.txt").write_text(SEGMENTS)
            args = [str(tmp / "report.tsv"), "--purity", "0.5", "-o", str(tmp / "out.tsv")]
            args += [a.replace("SEGMENTS", str(tmp / "segments.txt")) for a in extra]
            ccf.main(args)
            with (tmp / "out.tsv").open() as handle:
                return {row["Gene"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def test_clonal_heterozygous_on_diploid_segment(self):
        row = self.run_ccf("--segments", "SEGMENTS")["G1"]
        # VAF 0.25 at purity 0.5, CN 2: one mutated copy in every tumour cell
        self.assertEqual(row["CN Source"], "sequenza_segment")
        self.assertEqual(row["Mutated Copies"], "1.000")
        self.assertEqual(row["CCF"], "1.000")

    def test_amplified_segment_raises_multiplicity(self):
        row = self.run_ccf("--segments", "SEGMENTS")["G2"]
        # u = 0.6 * (0.5*4 + 2*0.5) / 0.5 = 3.6 -> m capped at major CN 3 -> CCF 1
        self.assertEqual(row["Tumor CN"], "4")
        self.assertEqual(row["Multiplicity"], "3")
        self.assertEqual(row["CCF"], "1.000")

    def test_falls_back_to_ploidy_between_segments(self):
        row = self.run_ccf("--segments", "SEGMENTS", "--ploidy", "3")["G3"]
        # u = 0.1 * (0.5*3 + 1) / 0.5 = 0.5 -> subclonal
        self.assertEqual(row["CN Source"], "ploidy")
        self.assertEqual(row["CCF"], "0.500")

    def test_without_ploidy_or_segments_assumes_normal_copy_number(self):
        rows = self.run_ccf()
        self.assertEqual(rows["G1"]["CN Source"], "assumed_normal")
        self.assertEqual(rows["G4"]["CCF"], "NA")

    def test_male_x_uses_single_normal_copy(self):
        self.assertEqual(ccf.ccf_columns(0.5, "chrX", 10, 1.0, None, "male", None)[4], "1")
        self.assertEqual(ccf.ccf_columns(0.5, "chrX", 10, 1.0, None, "female", None)[4], "2")


def tsv(*rows):
    return "".join("\t".join(str(cell) for cell in row) + "\n" for row in rows)


SEQ_REPORT = tsv(["ID", "Gene", "IC50 MT", "IC50 WT", "DNA VAF"],
                 ["chr1-99-100-A-T", "BRAF", 12, 900, 0.4])
FUSE_REPORT = tsv(["ID", "Gene", "IC50 MT", "Expr", "Read Support"],
                  ["EML4--ALK_chr2:42295516", "EML4--ALK", 35, 2.4, 30])

VCF_HEADER = ("##fileformat=VCFv4.2\n"
              + tsv(["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO",
                     "FORMAT", "normal", "tumor"]))


def rna_vcf(path, depth, vaf):
    """A readcount VCF with one variant, as bam-readcount annotation leaves it."""
    path.write_text(VCF_HEADER + tsv(["chr1", 100, ".", "A", "T", ".", "PASS", ".",
                                      "GT:RDP:RAF", "0/0:0:0", f"0/1:{depth}:{vaf}"]))
    return path


def first_row(path):
    with open(path, newline="") as handle:
        return next(csv.DictReader(handle, delimiter="\t"))


class VariantPosition(unittest.TestCase):
    def test_a_pvacseq_id(self):
        self.assertEqual(ccf.variant_position("chr1-499-500-A-T"), ("chr1", 500))

    def test_a_fusion_id_has_no_locus(self):
        self.assertIsNone(ccf.variant_position("EML4--ALK_chr2:42295516-chr2:29223528"))

    def test_a_short_id(self):
        self.assertIsNone(ccf.variant_position("something"))


class DirectoryInput(unittest.TestCase):
    """Pointed straight at a pVACseq or pVACfuse output directory."""

    def tree(self, text, tool):
        root = Path(tempfile.mkdtemp())
        for sub, label in (("MHC_Class_I", "MHC_I"), ("combined", "Combined")):
            directory = root / tool / sub
            directory.mkdir(parents=True)
            (directory / f"tumor.{label}.all_epitopes.aggregated.tsv").write_text(text)
        return root / tool

    def annotate(self, text, tool):
        root = self.tree(text, tool)
        out = root.parent / "out"
        ccf.main([str(root), "--purity", "0.5", "--outdir", str(out)])
        return out

    def test_every_class_report_is_annotated(self):
        out = self.annotate(SEQ_REPORT, "pvacseq")
        self.assertEqual(len(list(out.glob("*.ccf.tsv"))), 2)

    def test_the_ccf_column_is_filled(self):
        out = self.annotate(SEQ_REPORT, "pvacseq")
        self.assertNotEqual(first_row(next(out.glob("*.ccf.tsv")))["CCF"], "NA")

    def test_a_directory_with_o_instead_of_outdir(self):
        root = self.tree(SEQ_REPORT, "pvacseq")
        with self.assertRaises(SystemExit):
            ccf.main([str(root), "--purity", "0.5", "-o", str(root / "x.tsv")])

    def test_fusion_candidates_get_no_ccf(self):
        row = first_row(next(self.annotate(FUSE_REPORT, "pvacfuse").glob("*.ccf.tsv")))
        self.assertEqual(row["CN Source"], "not_applicable")
        self.assertEqual(row["CCF"], "NA")

    def test_fusion_rows_keep_their_own_columns(self):
        row = first_row(next(self.annotate(FUSE_REPORT, "pvacfuse").glob("*.ccf.tsv")))
        self.assertEqual(row["Read Support"], "30")

    def test_an_input_that_is_not_a_report(self):
        root = Path(tempfile.mkdtemp())
        with self.assertRaises(SystemExit):
            ccf.main([str(root / "nope"), "--purity", "0.5", "-o", str(root / "o.tsv")])


class Timepoints(unittest.TestCase):
    def test_both_timepoints_land_beside_the_ccf(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "report.tsv").write_text(SEQ_REPORT)
        baseline = rna_vcf(tmp / "baseline.vcf", 80, "0.40")
        week6 = rna_vcf(tmp / "week6.vcf", 55, "0.22")
        out = tmp / "out.tsv"
        ccf.main([str(tmp / "report.tsv"), "--purity", "0.5", "-o", str(out),
                  "--timepoint", f"Baseline={baseline}", "--timepoint", f"Week 6={week6}"])
        row = first_row(out)
        self.assertEqual(row["RNA Depth (Baseline)"], "80")
        self.assertEqual(row["RNA VAF (Week 6)"], "0.22")
        self.assertNotEqual(row["CCF"], "NA")

    def test_a_missing_timepoint_vcf(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "report.tsv").write_text(SEQ_REPORT)
        with self.assertRaises(SystemExit):
            ccf.main([str(tmp / "report.tsv"), "--purity", "0.5", "-o", str(tmp / "o.tsv"),
                      "--timepoint", f"Week 6={tmp / 'nope.vcf'}"])

    def test_a_malformed_timepoint(self):
        tmp = Path(tempfile.mkdtemp())
        (tmp / "report.tsv").write_text(SEQ_REPORT)
        with self.assertRaises(SystemExit):
            ccf.main([str(tmp / "report.tsv"), "--purity", "0.5", "-o", str(tmp / "o.tsv"),
                      "--timepoint", "no-equals-sign"])


if __name__ == "__main__":
    unittest.main()
