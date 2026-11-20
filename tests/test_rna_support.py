import csv
import importlib.util
from pathlib import Path
import sys
import tempfile
import unittest

BIN = Path(__file__).resolve().parents[1] / "bin"
sys.path.insert(0, str(BIN))
SPEC = importlib.util.spec_from_file_location("rna_support", BIN / "rna_support.py")
rs = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(rs)

CSQ_HEADER = ('##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. '
              'Format: Allele|Consequence|SYMBOL|HGVSp|Protein_position|Amino_acids|CANONICAL|MANE_SELECT">\n')


def csq(consequence, gene, hgvsp="", canonical="YES"):
    other = f"T|non_coding_transcript_exon_variant|{gene}|||||"
    return f"CSQ={other},T|{consequence}|{gene}|{hgvsp}|||{canonical}|"


def vcf(records):
    """records: (chrom, pos, ref, alt, info, dna (dp, alt), rna (dp, alt))"""
    lines = ["##fileformat=VCFv4.2\n", CSQ_HEADER,
             "#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\tFORMAT\tnormal\ttumor\n"]
    for chrom, pos, ref, alt, info, (dp, ad), (rdp, rad) in records:
        tumour = f"0/1:{dp - ad},{ad}:{dp}:{ad / dp:.3f}:{rdp - rad},{rad}:{rdp}:{(rad / rdp if rdp else 0):.3f}"
        lines.append(f"{chrom}\t{pos}\t.\t{ref}\t{alt}\t.\tPASS\t{info}\tGT:AD:DP:AF:RAD:RDP:RAF\t0/0:30,0:30:0\t{tumour}\n")
    return "".join(lines)


# purity 1, no ploidy -> CN 2: VAF 0.5 is clonal, VAF 0.1 is CCF 0.2
SITES = [
    ("chr1", 100, "C", "T", csq("missense_variant", "G1", "ENSP1.1:p.Val600Glu"), (100, 50), (40, 10), (40, 8)),   # both
    ("chr1", 200, "C", "T", csq("missense_variant", "G2"), (100, 50), (40, 10), (40, 0)),    # baseline only
    ("chr1", 300, "G", "A", csq("stop_gained", "G3"), (100, 50), (40, 0), (40, 6)),          # week 6 only
    ("chr1", 400, "G", "A", csq("missense_variant", "G4"), (100, 10), (0, 0), (0, 0)),       # DNA only, subclonal
    ("chr1", 500, "C", "T", csq("synonymous_variant", "G5"), (100, 50), (40, 20), (40, 20)), # not counted
    ("chr1", 600, "AT", "A", csq("frameshift_variant", "G6"), (100, 50), (40, 1), (40, 1)),  # 1 alt read: DNA only
]

CLASS_I = ("ID\tGene\tBest Peptide\tAllele\tIC50 MT\t%ile MT\n"
           "chr1-99-100-C-T\tG1\tPEPTIDEAA\tHLA-A*02:01\t45\t0.3\n"
           "chr1-299-300-G-A\tG3\tPEPTIDEBB\tHLA-A*02:01\t4000\t1.5\n"
           "chr1-399-400-G-A\tG4\tPEPTIDECC\tHLA-B*07:02\t120\t0.9\n"
           "chr1-199-200-C-T\tG2\tPEPTIDEDD\tHLA-B*07:02\t9000\t30\n"
           "chr1-600-601-AT-A\tG6\tPEPTIDEEE\tHLA-B*07:02\t80\t0.5\n")
CLASS_II = ("ID\tGene\tBest Peptide\tAllele\tIC50 MT\t%ile MT\n"
            "chr1-199-200-C-T\tG2\tLONGPEPTIDEAAAAA\tDRB1*15:01\t3000\t12\n")


class RnaSupport(unittest.TestCase):
    def run_support(self, labels=("Baseline", "Week 6")):
        tmp = Path(self.tmp.name)
        (tmp / "base.vcf").write_text(vcf([s[:6] + (s[6],) for s in SITES]))
        (tmp / "w6.vcf").write_text(vcf([s[:6] + (s[7],) for s in SITES]))
        (tmp / "i.tsv").write_text(CLASS_I)
        (tmp / "ii.tsv").write_text(CLASS_II)
        files = {"Baseline": tmp / "base.vcf", "Week 6": tmp / "w6.vcf"}
        args = ["--patient", "P1", "--vcf", str(tmp / "base.vcf"), "--purity", "1.0",
                "--class-i", str(tmp / "i.tsv"), "--class-ii", str(tmp / "ii.tsv"),
                "--no-plot", "-o", str(tmp / "P1")]
        for label in labels:
            args += ["--timepoint", f"{label}={files[label]}"]
        rs.main(args)
        with open(tmp / "P1.rna_support.variants.tsv") as handle:
            variants = {r["Gene"]: r for r in csv.DictReader(handle, delimiter="\t")}
        with open(tmp / "P1.rna_support.summary.tsv") as handle:
            summary = {r["Set"]: r for r in csv.DictReader(handle, delimiter="\t")}
        return variants, summary

    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()

    def tearDown(self):
        self.tmp.cleanup()

    def test_categories_for_two_timepoints(self):
        variants, _ = self.run_support()
        self.assertEqual(variants["G1"]["RNA Support"], "Expressed both")
        self.assertEqual(variants["G2"]["RNA Support"], "Baseline only")
        self.assertEqual(variants["G3"]["RNA Support"], "Week 6 only")
        self.assertEqual(variants["G4"]["RNA Support"], "DNA only")
        self.assertEqual(variants["G6"]["RNA Support"], "DNA only")

    def test_only_non_synonymous_counted(self):
        variants, summary = self.run_support()
        self.assertNotIn("G5", variants)
        self.assertEqual(summary["non_synonymous"]["Total"], "5")

    def test_no_coverage_kept_apart_from_no_alt_reads(self):
        variants, _ = self.run_support()
        self.assertEqual(variants["G4"]["RNA Status (Baseline)"], "no_coverage")
        self.assertEqual(variants["G6"]["RNA Status (Baseline)"], "covered_no_alt")

    def test_clonality_and_neoantigen_sets(self):
        variants, summary = self.run_support()
        self.assertEqual(variants["G4"]["Clonal"], "no")
        self.assertEqual(variants["G1"]["Clonal"], "yes")
        # class I binders: G1 (IC50), G3 (percentile), G4 (subclonal), G6 (indel ID with start = POS)
        self.assertEqual(summary["neoag_class_i"]["Total"], "4")
        self.assertEqual(summary["neoag_class_i_clonal"]["Total"], "3")
        self.assertEqual(summary["neoag_class_i_clonal"]["Expressed both"], "1")
        # G2's class II row misses both thresholds
        self.assertEqual(summary["neoag_class_ii"]["Total"], "0")

    def test_single_timepoint(self):
        variants, summary = self.run_support(labels=("Baseline",))
        self.assertEqual(variants["G3"]["RNA Support"], "DNA only")
        self.assertEqual(variants["G2"]["RNA Support"], "Expressed (Baseline)")
        self.assertIn("Expressed (Baseline)", summary["non_synonymous"])

    def test_protein_change_shortened(self):
        self.assertEqual(rs.short_protein_change("ENSP1.1:p.Val600Glu"), "V600E")
        self.assertEqual(rs.short_protein_change("ENSP1.1:p.Arg80Ter"), "R80*")

    def test_categories_generalise_beyond_two(self):
        cats = rs.categories(["T0", "T1", "T2"])
        self.assertEqual(cats[0], "Expressed all")
        self.assertEqual(cats[-1], "DNA only")
        self.assertEqual(rs.category_for({"T0", "T2"}, ["T0", "T1", "T2"]), "Expressed T0 + T2")


if __name__ == "__main__":
    unittest.main()
