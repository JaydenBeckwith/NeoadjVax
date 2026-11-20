import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
SPEC = importlib.util.spec_from_file_location("patient_report", BIN / "patient_report.py")
pr = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(pr)

CSQ_HEADER = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
              'Format: Allele|Consequence|SYMBOL|HGVSp|Protein_position|Amino_acids|CANONICAL|MANE_SELECT">\n')


def record(chrom, pos, ref, alt, consequence, gene, hgvsp="", dp=100, ad=40):
    info = f"CSQ={alt}|{consequence}|{gene}|{hgvsp}|||YES|"
    tumour = f"0/1:{dp - ad},{ad}:{dp}:{ad / dp:.3f}:20,10:30:0.333"
    return f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:AD:DP:AF:RAD:RDP:RAF\t0/0:30,0:30:0\t{tumour}\n"


HEADER = "##fileformat=VCFv4.2\n" + CSQ_HEADER + "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnormal\ttumor\n"
# chr1: ...CCAGTT... so pos 3 is C with C 5' (dipyrimidine); pos 8 is G with A 3' on the plus strand
GENOME = ">chr1\nACCAGTTGAAGGCTTACG\n"


def vcf(*records):
    return HEADER + "".join(records)


class Subtype(unittest.TestCase):
    def drivers(self, *records):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "x.vcf"
            path.write_text(vcf(*records))
            variants = pr.rs.read_vcf(path)
        genes = {"BRAF", "NRAS", "KRAS", "HRAS", "NF1", "CDKN2A"}
        return pr.driver_mutations(variants, genes, {})

    def test_braf_v600_is_class_i_and_wins(self):
        drivers = self.drivers(record("chr1", 3, "C", "T", "missense_variant", "BRAF", "p.Val600Glu"),
                               record("chr1", 8, "G", "A", "stop_gained", "NF1", "p.Arg100Ter"))
        subtype, braf = pr.melanoma_subtype(drivers)
        self.assertEqual(subtype, "BRAF hotspot")
        self.assertIn("V600E", braf)
        self.assertIn("class I", braf)

    def test_nras_q61_is_ras_subtype(self):
        subtype, braf = pr.melanoma_subtype(self.drivers(record("chr1", 3, "C", "T", "missense_variant", "NRAS", "p.Gln61Arg")))
        self.assertTrue(subtype.startswith("RAS hotspot"))
        self.assertEqual(braf, "BRAF wild-type")

    def test_non_hotspot_braf_reported_as_other(self):
        subtype, braf = pr.melanoma_subtype(self.drivers(record("chr1", 3, "C", "T", "missense_variant", "BRAF", "p.Ala100Thr")))
        self.assertEqual(subtype, "Triple wild-type")
        self.assertEqual(braf, "BRAF other (A100T)")

    def test_nf1_then_triple_wild_type(self):
        self.assertEqual(pr.melanoma_subtype(self.drivers(record("chr1", 8, "G", "A", "stop_gained", "NF1", "p.Arg100Ter")))[0], "NF1 mutant")
        self.assertEqual(pr.melanoma_subtype([])[0], "Triple wild-type")


class Helpers(unittest.TestCase):
    def test_tex_escaping(self):
        self.assertEqual(pr.tex("a_b & 5% #1"), r"a\_b \& 5\% \#1")

    def test_hla_split_and_zygosity(self):
        class_i, class_ii = pr.split_hla("HLA-A*02:01,HLA-A*02:01,HLA-B*07:02,HLA-B*08:01,DRB1*15:01")
        self.assertEqual(class_ii, ["DRB1*15:01"])
        zygosity = {row[0]: row[2] for row in pr.class_i_zygosity(class_i)}
        self.assertEqual(zygosity["HLA-A"], "homozygous")
        self.assertEqual(zygosity["HLA-B"], "heterozygous")
        self.assertEqual(zygosity["HLA-C"], "not typed")

    def test_bed_size_merges_overlaps(self):
        with tempfile.TemporaryDirectory() as tmp:
            bed = Path(tmp) / "t.bed"
            bed.write_text("chr1\t0\t1000000\nchr1\t500000\t1500000\nchr2\t0\t500000\n")
            self.assertAlmostEqual(pr.bed_size_mb(bed), 2.0)

    def test_uv_signature_uses_reference_context(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "g.fa").write_text(GENOME)
            (tmp / "g.fa.fai").write_text("chr1\t18\t6\t18\t19\n")
            (tmp / "x.vcf").write_text(vcf(record("chr1", 3, "C", "T", "missense_variant", "G1"),
                                           record("chr1", 8, "G", "A", "missense_variant", "G2"),
                                           record("chr1", 10, "A", "G", "missense_variant", "G3")))
            fasta = pr.Fasta(tmp / "g.fa")
            self.assertEqual(fasta.base("chr1", 3), "C")
            uv = pr.uv_signature(pr.rs.read_vcf(tmp / "x.vcf"), fasta)
            fasta.close()
        self.assertEqual(uv["ct"], 2)
        self.assertAlmostEqual(uv["dipyrimidine_fraction"], 1.0)


FUSION_HEADER = "ID\tGene\tBest Peptide\tBest Transcript\tAllele\tIC50 MT\t%ile MT\tExpr\tRead Support\tTier\n"
FUSION_ROWS = (FUSION_HEADER
               + "1\tEML4_ALK\tKLMNPQRST\tENST1\tHLA-A*02:01\t45\t0.3\t2.5\t30\tPass\n"
               + "2\tAAA_BBB\tPQRSTUVWY\tENST2\tHLA-A*02:01\t80\t0.6\t0.05\t2\tLowReadSupport\n"
               + "3\tCCC_DDD\tWYACDEFGH\tENST3\tHLA-A*02:01\t9000\t40\t5.0\t50\tPoor\n")


class Fusions(unittest.TestCase):
    def test_status_from_the_directory(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            self.assertEqual(pr.read_pvacfuse(None)[0], "not_run")
            self.assertEqual(pr.read_pvacfuse(tmp / "missing")[0], "not_run")
            (tmp / "pvacfuse").mkdir()
            (tmp / "pvacfuse" / "no_fusion_proteins.txt").write_text("nothing to predict\n")
            self.assertEqual(pr.read_pvacfuse(tmp / "pvacfuse")[0], "no_proteins")

    def test_reads_per_class_reports(self):
        with tempfile.TemporaryDirectory() as tmp:
            d = Path(tmp) / "pvacfuse"
            (d / "MHC_Class_I").mkdir(parents=True)
            (d / "MHC_Class_II").mkdir(parents=True)
            (d / "MHC_Class_I" / "P1.MHC_I.all_epitopes.aggregated.tsv").write_text(FUSION_ROWS)
            (d / "MHC_Class_II" / "P1.MHC_II.all_epitopes.aggregated.tsv").write_text(FUSION_HEADER)
            status, reports = pr.read_pvacfuse(d)
        self.assertEqual(status, "ok")
        self.assertEqual(len(reports["I"]), 3)
        self.assertEqual(reports["II"], [])

    def test_ranking_prefers_supported_fusions(self):
        rows = [r for r in csv.DictReader(FUSION_ROWS.splitlines(), delimiter="\t")]
        ranked = pr.rank_fusions(rows, 500, 2.0, 5.0, 0.1)
        self.assertEqual([r["Gene"] for r in ranked], ["EML4_ALK", "AAA_BBB"])   # CCC_DDD is not a binder

    def test_melanoma_partner_note(self):
        self.assertIn("Spitz", pr.fusion_note("EML4_ALK"))
        self.assertIn("resistance", pr.fusion_note("AGK-BRAF"))
        self.assertEqual(pr.fusion_note("AAA_BBB"), "")

    def test_counts_unique_fusion_names(self):
        with tempfile.TemporaryDirectory() as tmp:
            path = Path(tmp) / "sf.tsv"
            path.write_text("#FusionName\tJunctionReadCount\nA--B\t10\nA--B\t4\nC--D\t7\n")
            self.assertEqual(pr.count_fusions(path), 2)
        self.assertIsNone(pr.count_fusions(None))


class Report(unittest.TestCase):
    def test_fusion_section_renders_when_the_branch_ran(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "x.vcf").write_text(vcf(record("chr1", 3, "C", "T", "missense_variant", "BRAF", "p.Val600Glu")))
            d = tmp / "pvacfuse"
            (d / "MHC_Class_I").mkdir(parents=True)
            (d / "MHC_Class_I" / "P1.MHC_I.all_epitopes.aggregated.tsv").write_text(FUSION_ROWS)
            (tmp / "sf.tsv").write_text("#FusionName\tJunctionReadCount\nEML4--ALK\t30\n")
            pr.main(["--patient", "P1", "--vcf", str(tmp / "x.vcf"), "--purity", "0.6", "--callable-mb", "10",
                     "--pvacfuse-dir", str(d), "--fusion-predictions", str(tmp / "sf.tsv"),
                     "-o", str(tmp / "r.tex")])
            text = (tmp / "r.tex").read_text()
        self.assertIn("STAR-Fusion called 1 fusion for this patient", text)
        self.assertIn("EML4--ALK", text)
        self.assertNotIn("CCC\\_DDD", text)        # not a binder
        self.assertIn("Spitz", text)

    def test_fusion_section_when_no_coding_fusion_protein(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "x.vcf").write_text(vcf(record("chr1", 3, "C", "T", "missense_variant", "BRAF", "p.Val600Glu")))
            d = tmp / "pvacfuse"
            d.mkdir()
            (d / "no_fusion_proteins.txt").write_text("nothing\n")
            pr.main(["--patient", "P1", "--vcf", str(tmp / "x.vcf"), "--purity", "0.6", "--callable-mb", "10",
                     "--pvacfuse-dir", str(d), "-o", str(tmp / "r.tex")])
            text = (tmp / "r.tex").read_text()
        self.assertIn("no coding fusion protein", text)
        self.assertIn("not a failure", text)

    def test_minimal_report_renders_with_missing_sections_marked(self):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            (tmp / "x.vcf").write_text(vcf(record("chr1", 3, "C", "T", "missense_variant", "BRAF", "p.Val600Glu")))
            pr.main(["--patient", "P_1", "--vcf", str(tmp / "x.vcf"), "--purity", "0.6",
                     "--callable-mb", "0.1", "-o", str(tmp / "r.tex")])
            text = (tmp / "r.tex").read_text()
        self.assertIn(r"\begin{document}", text)
        self.assertIn("P\\_1", text)
        self.assertIn("BRAF hotspot", text)
        self.assertIn("TMB-high", text)                   # 1 variant / 0.1 Mb
        self.assertIn("RNA support was not assessed", text)
        self.assertIn("fusion branch did not run", text)
        self.assertIn("Allele-specific HLA LOH was not assessed", text)
        self.assertIn("Research use only", text)
        self.assertTrue(text.rstrip().endswith(r"\end{document}"))


if __name__ == "__main__":
    unittest.main()
