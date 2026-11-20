import csv
import importlib.util
from pathlib import Path
import tempfile
import unittest

SPEC = importlib.util.spec_from_file_location(
    "gene_panel_report", Path(__file__).resolve().parents[1] / "bin" / "gene_panel_report.py")
panel = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(panel)

SKIN = "Skin - Sun Exposed (Lower leg)"

GTF = "".join(
    f'{chrom}\tx\tgene\t{start}\t{end}\t.\t+\t.\tgene_id "{gid}"; gene_name "{name}";\n'
    for chrom, start, end, gid, name in [
        ("chr15", 44711487, 44718877, "ENSG00000166710.21", "B2M"),
        ("chr7", 140719327, 140924929, "ENSG00000157764.13", "BRAF"),
        ("chr1", 64833229, 64966411, "ENSG00000162434.11", "JAK1"),
        ("chr17", 7668402, 7687550, "ENSG00000141510.16", "TP53"),
    ])

# B2M on a 1-copy LOH segment, BRAF amplified, JAK1 copy-neutral
SEGMENTS = '"chromosome"\t"start.pos"\t"end.pos"\t"CNt"\t"A"\t"B"\n' \
           '"chr15"\t40000000\t50000000\t1\t1\t0\n' \
           '"chr7"\t140000000\t141000000\t6\t4\t2\n' \
           '"chr1"\t60000000\t70000000\t2\t1\t1\n'

GTEX = "ENSEMBL_ID\tGENE_NAME\tGROUP\tQUANTILE\tTPM\n" + "".join(
    f"{gid}\t{name}\t{SKIN}\t{q}\t{tpm}\n"
    for gid, name, low, high in [("ENSG00000166710", "B2M", 300, 800),
                                 ("ENSG00000157764", "BRAF", 5, 15),
                                 ("ENSG00000162434", "JAK1", 20, 40)]
    for q, tpm in (("0.25", low), ("0.75", high)))

GTEX_MEDIAN_ONLY = "ENSEMBL_ID\tGENE_NAME\tGROUP\tQUANTILE\tTPM\n" \
                   f"ENSG00000166710\tB2M\t{SKIN}\t0.5\t500.0000\n"

SAMPLE = "ENSEMBL_ID\tGENE_NAME\tTPM\tEXPRESSION_QUARTILE\n" \
         "ENSG00000166710\tB2M\t12.0000\t1\n" \
         "ENSG00000157764\tBRAF\t60.0000\t4\n" \
         "ENSG00000162434\tJAK1\t25.0000\t3\n"

VEP = '##fileformat=VCFv4.2\n' \
      '##INFO=<ID=CSQ,Number=.,Type=String,Description="Consequence annotations from Ensembl VEP. ' \
      'Format: Allele|Consequence|SYMBOL|HGVSp|Amino_acids|Protein_position">\n' \
      '#CHROM\tPOS\tID\tREF\tALT\tQUAL\tFILTER\tINFO\n' \
      'chr7\t140753336\t.\tA\tT\t.\tPASS\tCSQ=T|missense_variant|BRAF|ENSP00000288602.6:p.Val600Glu|V/E|600\n' \
      'chr1\t64844000\t.\tC\tCA\t.\tPASS\tCSQ=CA|frameshift_variant|JAK1|ENSP00000343204.3:p.Ser500fs|S/X|500\n' \
      'chr9\t1\t.\tC\tT\t.\tPASS\tCSQ=T|missense_variant|NOTINPANEL|ENSP1:p.Ala1Val|A/V|1\n'


class GenePanelReport(unittest.TestCase):
    def run_panel(self, *extra, files=None):
        with tempfile.TemporaryDirectory() as tmp:
            tmp = Path(tmp)
            default = {"gtf": GTF, "seg.txt": SEGMENTS, "gtex.tsv": GTEX,
                       "sample.tsv": SAMPLE, "vep.vcf": VEP, "median.tsv": GTEX_MEDIAN_ONLY}
            default.update(files or {})
            for name, text in default.items():
                (tmp / name).write_text(text)
            args = [a.replace("TMP", str(tmp)) for a in extra] + ["-o", str(tmp / "out.tsv")]
            panel.main(args)
            with (tmp / "out.tsv").open() as handle:
                return {row["Gene"]: row for row in csv.DictReader(handle, delimiter="\t")}

    def full(self, *extra):
        return self.run_panel("--sample-tpm", "TMP/sample.tsv",
                              "--gtex", "TMP/gtex.tsv", "--gtex-tissue", SKIN,
                              "--vep-vcf", "TMP/vep.vcf",
                              "--segments", "TMP/seg.txt", "--gtf", "TMP/gtf",
                              "--ploidy", "2", *extra)

    def test_panel_covers_presentation_and_drivers(self):
        rows = self.full()
        for gene in ("HLA-A", "B2M", "TAP1", "PSMB9", "CIITA", "NLRC5", "JAK1",
                     "BRAF", "NRAS", "NF1", "TP53", "CDKN2A", "PTEN"):
            self.assertIn(gene, rows, gene)

    def test_pathways_group_the_panel(self):
        self.assertEqual(self.full()["B2M"]["Pathway"], "apm_class_i")
        self.assertEqual(self.full()["BRAF"]["Pathway"], "driver_mapk")

    def test_collapsed_expression_and_loh_flag_loss_of_function(self):
        # the classic immune-escape finding: B2M down and single-copy with LOH
        row = self.full()["B2M"]
        self.assertEqual(row["GTEx Call"], "DOWN")
        self.assertEqual(row["CN Call"], "LOSS,LOH")
        self.assertEqual(row["Panel Flag"], "LOSS_OF_FUNCTION")

    def test_truncating_mutation_alone_flags_loss_of_function(self):
        # JAK1 reads normally and is copy-neutral; the frameshift is the finding
        row = self.full()["JAK1"]
        self.assertEqual(row["GTEx Call"], "NORMAL")
        self.assertEqual(row["CN Call"], "NEUTRAL")
        self.assertEqual(row["Truncating"], "yes")
        self.assertEqual(row["Panel Flag"], "LOSS_OF_FUNCTION")

    def test_activating_driver_is_altered_not_loss_of_function(self):
        row = self.full()["BRAF"]
        self.assertIn("p.Val600Glu", row["Mutations"])
        self.assertEqual(row["Truncating"], "no")
        self.assertEqual(row["CN Call"], "AMP")
        self.assertEqual(row["Panel Flag"], "ALTERED")

    def test_intact_gene_is_ok(self):
        # JAK1 without the VCF: expressed normally, copy-neutral, unmutated
        rows = self.run_panel("--sample-tpm", "TMP/sample.tsv",
                              "--gtex", "TMP/gtex.tsv", "--gtex-tissue", SKIN,
                              "--segments", "TMP/seg.txt", "--gtf", "TMP/gtf", "--ploidy", "2")
        row = rows["JAK1"]
        self.assertEqual((row["GTEx Call"], row["CN Call"], row["Mutations"]), ("NORMAL", "NEUTRAL", "-"))
        self.assertEqual(row["Panel Flag"], "OK")

    def test_gene_with_no_data_at_all_is_na(self):
        self.assertEqual(self.full()["HLA-DRB1"]["Panel Flag"], "NA")

    def test_mutations_outside_the_panel_are_ignored(self):
        self.assertNotIn("NOTINPANEL", self.full())

    def test_copy_number_is_called_against_sample_ploidy(self):
        # at ploidy 4, two copies is a relative loss
        rows = self.run_panel("--segments", "TMP/seg.txt", "--gtf", "TMP/gtf", "--ploidy", "4")
        self.assertIn("LOSS", rows["JAK1"]["CN Call"])
        self.assertEqual(self.full()["JAK1"]["CN Call"], "NEUTRAL")

    def test_median_only_table_falls_back_to_fold_change(self):
        rows = self.run_panel("--sample-tpm", "TMP/sample.tsv",
                              "--gtex", "TMP/median.tsv", "--gtex-tissue", SKIN, "--fold", "4")
        row = rows["B2M"]
        self.assertEqual(row["GTEx Basis"], "median_x4")
        self.assertEqual(row["GTEx Range"], "125.00-2000.00")
        self.assertEqual(row["GTEx Call"], "DOWN")

    def test_quantile_range_is_used_when_the_table_has_one(self):
        self.assertEqual(self.full()["B2M"]["GTEx Basis"], "q0.25-q0.75")

    def test_runs_with_no_inputs_but_the_panel(self):
        rows = self.run_panel()
        self.assertTrue(all(row["Panel Flag"] == "NA" for row in rows.values()))

    def test_tolerates_an_empty_vcf(self):
        rows = self.run_panel("--vep-vcf", "TMP/empty.vcf", files={"empty.vcf": ""})
        self.assertEqual(rows["BRAF"]["Mutations"], "-")

    def test_custom_panel_replaces_the_built_in_one(self):
        rows = self.run_panel("--panel", "TMP/p.tsv", "--vep-vcf", "TMP/vep.vcf",
                              files={"p.tsv": "gene\tpathway\nBRAF\tmy_genes\n"})
        self.assertEqual(list(rows), ["BRAF"])
        self.assertEqual(rows["BRAF"]["Pathway"], "my_genes")

    def test_segments_without_a_gtf_is_rejected(self):
        with self.assertRaises(SystemExit):
            self.run_panel("--segments", "TMP/seg.txt")

    def test_unknown_tissue_name_fails_loudly(self):
        with self.assertRaises(SystemExit):
            self.run_panel("--gtex", "TMP/gtex.tsv", "--gtex-tissue", "Pineal Gland")

    def test_gene_straddling_two_segments_reports_lowest_minor_copy(self):
        segments = '"chromosome"\t"start.pos"\t"end.pos"\t"CNt"\t"A"\t"B"\n' \
                   '"chr17"\t7000000\t7670000\t2\t1\t1\n' \
                   '"chr17"\t7670001\t8000000\t2\t2\t0\n'
        rows = self.run_panel("--segments", "TMP/s2.txt", "--gtf", "TMP/gtf", "--ploidy", "2",
                              files={"s2.txt": segments})
        self.assertIn("LOH", rows["TP53"]["CN Call"])


if __name__ == "__main__":
    unittest.main()
