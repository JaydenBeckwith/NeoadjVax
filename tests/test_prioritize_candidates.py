"""prioritize_candidates.py: flags, reject reasons and ranking."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "prioritize_candidates.py"
spec = importlib.util.spec_from_file_location("prioritize_candidates", BIN)
prioritize = importlib.util.module_from_spec(spec)
spec.loader.exec_module(prioritize)

COLUMNS = ["ID", "Gene", "Allele", "Best Peptide", "IC50 MT", "IC50 WT", "%ile MT",
           "RNA Expr", "RNA VAF", "RNA Depth", "DNA VAF", "DNA Depth", "CCF",
           "Expression Flag", "GTEx Specificity", "Sample TPM"]

# A clean candidate: strong binder, foreign, clonal, expressed, RNA-supported.
GOOD = {
    "ID": "chr1-100-100-A-T", "Gene": "GENE1", "Allele": "HLA-A*02:01",
    "Best Peptide": "SLYNTVATL", "IC50 MT": "12", "IC50 WT": "900", "%ile MT": "0.1",
    "RNA Expr": "40", "RNA VAF": "0.4", "RNA Depth": "80", "DNA VAF": "0.45",
    "DNA Depth": "120", "CCF": "0.98", "Expression Flag": "OK",
    "GTEx Specificity": "tissue_restricted", "Sample TPM": "40",
}

LOH_COLUMNS = ["patient", "locus", "hla_class", "allele1_pvacseq", "allele2_pvacseq",
               "call", "lost_allele_pvacseq"]


def write_report(rows):
    path = Path(tempfile.mkdtemp()) / "report.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({**{c: "NA" for c in COLUMNS}, **row})
    return path


def write_loh(rows, directory):
    path = Path(directory) / "loh.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=LOH_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({**{c: "" for c in LOH_COLUMNS}, **row})
    return path


def run(rows, loh_rows=None, **kwargs):
    report = write_report(rows)
    out = report.parent / "out.tsv"
    argv = [str(report), "-o", str(out)]
    if loh_rows is not None:
        argv += ["--hla-loh", str(write_loh(loh_rows, report.parent))]
    for key, value in kwargs.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    prioritize.main(argv)
    with out.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def flags_of(row):
    return set(filter(None, row["Flags"].split(";")))


class TestCleanCandidate(unittest.TestCase):
    def test_no_flags(self):
        row = run([GOOD])[0]
        self.assertEqual(flags_of(row), set())
        self.assertEqual(row["Selection"], "Top ranking")
        self.assertEqual(row["Reject Reason"], "")
        self.assertEqual(row["Priority Rank"], "1")

    def test_agretopicity_and_dai(self):
        row = run([GOOD])[0]
        self.assertAlmostEqual(float(row["Agretopicity"]), 12 / 900, places=4)
        self.assertAlmostEqual(float(row["DAI"]), 888, places=1)

    def test_original_columns_are_kept(self):
        row = run([GOOD])[0]
        self.assertEqual(row["Gene"], "GENE1")
        self.assertEqual(row["Best Peptide"], "SLYNTVATL")


class TestFlags(unittest.TestCase):
    def test_weak_binder_is_rejected(self):
        row = run([{**GOOD, "IC50 MT": "4000", "%ile MT": "8"}])[0]
        self.assertIn("LOW_BINDING", flags_of(row))
        self.assertEqual(row["Selection"], "Rejected")
        self.assertEqual(row["Reject Reason"], "predicted binding too weak")

    def test_weak_on_one_measure_only_is_not_rejected(self):
        # Strong IC50, poor percentile: the two predictors disagree, so it stays.
        row = run([{**GOOD, "%ile MT": "8"}])[0]
        self.assertNotIn("LOW_BINDING", flags_of(row))

    def test_self_like_is_a_caution(self):
        row = run([{**GOOD, "IC50 MT": "400", "IC50 WT": "420"}])[0]
        self.assertIn("SELF_LIKE", flags_of(row))
        self.assertEqual(row["Selection"], "Caution")

    def test_subclonal(self):
        row = run([{**GOOD, "CCF": "0.2"}])[0]
        self.assertIn("SUBCLONAL", flags_of(row))

    def test_no_expression_is_rejected(self):
        row = run([{**GOOD, "Expression Flag": "NO_EXPRESSION"}])[0]
        self.assertIn("NO_EXPRESSION", flags_of(row))
        self.assertEqual(row["Selection"], "Rejected")

    def test_missing_expression_data_is_only_a_caution(self):
        row = run([{**GOOD, "Expression Flag": "NO_DATA"}])[0]
        self.assertIn("NO_EXPRESSION_DATA", flags_of(row))
        self.assertEqual(row["Selection"], "Caution")

    def test_rna_covered_but_unsupported_is_rejected(self):
        row = run([{**GOOD, "RNA Depth": "60", "RNA VAF": "0"}])[0]
        self.assertIn("NO_RNA_SUPPORT", flags_of(row))

    def test_no_rna_coverage_is_only_a_caution(self):
        row = run([{**GOOD, "RNA Depth": "NA", "RNA VAF": "NA"}])[0]
        self.assertIn("NO_RNA_DATA", flags_of(row))
        self.assertNotIn("NO_RNA_SUPPORT", flags_of(row))
        self.assertEqual(row["Selection"], "Caution")

    def test_broadly_expressed_normal_tissue(self):
        row = run([{**GOOD, "GTEx Specificity": "broadly_expressed"}])[0]
        self.assertIn("NORMAL_TISSUE", flags_of(row))

    def test_low_vaf_and_depth(self):
        row = run([{**GOOD, "DNA VAF": "0.01", "DNA Depth": "8"}])[0]
        self.assertIn("LOW_DNA_VAF", flags_of(row))
        self.assertIn("LOW_DEPTH", flags_of(row))

    def test_thresholds_are_arguments(self):
        row = run([{**GOOD, "CCF": "0.2"}], min_ccf=0.1)[0]
        self.assertNotIn("SUBCLONAL", flags_of(row))


class TestHlaLoh(unittest.TestCase):
    def test_candidate_on_a_lost_allele_is_rejected(self):
        loh = [{"locus": "A", "allele1_pvacseq": "HLA-A*02:01", "allele2_pvacseq": "HLA-A*01:01",
                "call": "loss", "lost_allele_pvacseq": "HLA-A*02:01"}]
        row = run([GOOD], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "lost")
        self.assertIn("HLA_LOST", flags_of(row))
        self.assertEqual(row["Reject Reason"], "restricting HLA allele lost in the tumour")

    def test_the_surviving_allele_of_a_lost_locus_is_retained(self):
        loh = [{"locus": "A", "allele1_pvacseq": "HLA-A*02:01", "allele2_pvacseq": "HLA-A*01:01",
                "call": "loss", "lost_allele_pvacseq": "HLA-A*01:01"}]
        row = run([GOOD], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "retained")
        self.assertEqual(flags_of(row), set())

    def test_not_evaluable_locus_is_a_caution(self):
        loh = [{"locus": "A", "allele1_pvacseq": "HLA-A*02:01", "allele2_pvacseq": "HLA-A*02:01",
                "call": "not_evaluable", "lost_allele_pvacseq": ""}]
        row = run([GOOD], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "not_evaluable")
        self.assertIn("HLA_NOT_EVALUABLE", flags_of(row))
        self.assertEqual(row["Selection"], "Caution")

    def test_allele_absent_from_the_loh_table(self):
        loh = [{"locus": "B", "allele1_pvacseq": "HLA-B*07:02", "allele2_pvacseq": "HLA-B*08:01",
                "call": "retained", "lost_allele_pvacseq": ""}]
        row = run([GOOD], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "NA")
        self.assertEqual(flags_of(row), set())

    def test_no_loh_file_leaves_the_column_na(self):
        row = run([GOOD])[0]
        self.assertEqual(row["HLA LOH"], "NA")

    def test_class_ii_pair_is_lost_if_either_chain_is(self):
        loh = [
            {"locus": "DQA1", "allele1_pvacseq": "DQA1*01:02", "allele2_pvacseq": "DQA1*05:01",
             "call": "retained", "lost_allele_pvacseq": ""},
            {"locus": "DQB1", "allele1_pvacseq": "DQB1*06:02", "allele2_pvacseq": "DQB1*02:01",
             "call": "loss", "lost_allele_pvacseq": "DQB1*06:02"},
        ]
        row = run([{**GOOD, "Allele": "DQA1*01:02-DQB1*06:02"}], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "lost")

    def test_class_ii_pair_retained_when_both_chains_are(self):
        loh = [
            {"locus": "DQA1", "allele1_pvacseq": "DQA1*01:02", "allele2_pvacseq": "DQA1*05:01",
             "call": "retained", "lost_allele_pvacseq": ""},
            {"locus": "DQB1", "allele1_pvacseq": "DQB1*06:02", "allele2_pvacseq": "DQB1*02:01",
             "call": "retained", "lost_allele_pvacseq": ""},
        ]
        row = run([{**GOOD, "Allele": "DQA1*01:02-DQB1*06:02"}], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "retained")

    def test_drb1_is_not_split_on_its_hyphen(self):
        loh = [{"locus": "DRB1", "allele1_pvacseq": "DRB1*15:01", "allele2_pvacseq": "DRB1*03:01",
                "call": "loss", "lost_allele_pvacseq": "DRB1*15:01"}]
        row = run([{**GOOD, "Allele": "DRB1*15:01"}], loh_rows=loh)[0]
        self.assertEqual(row["HLA LOH"], "lost")


class TestRanking(unittest.TestCase):
    def test_rejected_candidates_rank_last(self):
        rows = run([
            {**GOOD, "ID": "bad", "Expression Flag": "NO_EXPRESSION"},
            {**GOOD, "ID": "good"},
        ])
        self.assertEqual([r["ID"] for r in rows], ["good", "bad"])
        self.assertEqual([r["Priority Rank"] for r in rows], ["1", "2"])

    def test_clonal_outranks_subclonal(self):
        rows = run([
            {**GOOD, "ID": "sub", "CCF": "0.3"},
            {**GOOD, "ID": "clonal", "CCF": "0.95"},
        ])
        self.assertEqual([r["ID"] for r in rows], ["clonal", "sub"])

    def test_more_foreign_outranks_self_like_at_equal_clonality(self):
        rows = run([
            {**GOOD, "ID": "selfish", "IC50 MT": "100", "IC50 WT": "250"},
            {**GOOD, "ID": "foreign", "IC50 MT": "100", "IC50 WT": "9000"},
        ])
        self.assertEqual([r["ID"] for r in rows], ["foreign", "selfish"])

    def test_fewer_cautions_outrank_more(self):
        rows = run([
            {**GOOD, "ID": "two", "CCF": "0.2", "DNA VAF": "0.01"},
            {**GOOD, "ID": "one", "CCF": "0.2"},
        ])
        self.assertEqual([r["ID"] for r in rows], ["one", "two"])

    def test_stable_for_identical_rows(self):
        rows = run([{**GOOD, "ID": "b"}, {**GOOD, "ID": "a"}])
        self.assertEqual([r["ID"] for r in rows], ["a", "b"])


class TestEdges(unittest.TestCase):
    def test_empty_report_writes_a_header(self):
        rows = run([])
        self.assertEqual(rows, [])

    def test_all_values_missing(self):
        row = run([{c: "NA" for c in COLUMNS}])[0]
        self.assertEqual(row["Agretopicity"], "NA")
        self.assertEqual(row["DAI"], "NA")
        # No binding data at all is not a binding rejection.
        self.assertNotIn("LOW_BINDING", flags_of(row))

    def test_zero_wild_type_ic50_does_not_divide_by_zero(self):
        row = run([{**GOOD, "IC50 WT": "0"}])[0]
        self.assertEqual(row["Agretopicity"], "NA")


if __name__ == "__main__":
    unittest.main()
