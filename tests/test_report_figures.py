"""report_figures.py: the data behind each figure, and the files it writes.

The data preparation is pure and is tested unconditionally. The drawing
needs matplotlib, which the pVACtools container has but a bare checkout
may not, so those tests skip rather than fail when it is missing.
"""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "report_figures.py"
spec = importlib.util.spec_from_file_location("report_figures", BIN)
figures = importlib.util.module_from_spec(spec)
spec.loader.exec_module(figures)

HAS_MATPLOTLIB = figures.pyplot() is not None

CANDIDATE_COLUMNS = ["ID", "Gene", "CCF", "Agretopicity", "Sample TPM", "RNA Expr",
                     "Selection", "Flags", "Priority Rank"]
LOH_COLUMNS = ["patient", "locus", "hla_class", "allele1_pvacseq", "allele2_pvacseq",
               "call", "lost_allele_pvacseq"]


def candidate(id_, **kwargs):
    row = {c: "NA" for c in CANDIDATE_COLUMNS}
    row.update({"ID": id_, "Gene": "GENE1", "CCF": "0.9", "Agretopicity": "0.05",
                "Sample TPM": "20", "Selection": "Top ranking", "Flags": "", "Priority Rank": "1"})
    row.update(kwargs)
    return row


def write(rows, columns, name="t.tsv", directory=None):
    path = Path(directory or tempfile.mkdtemp()) / name
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    return path


class TestNumber(unittest.TestCase):
    def test_values_and_blanks(self):
        self.assertEqual(figures.number("0.5"), 0.5)
        self.assertIsNone(figures.number("NA"))
        self.assertIsNone(figures.number(""))
        self.assertIsNone(figures.number(None))
        self.assertIsNone(figures.number("not a number"))

    def test_infinities_are_not_values(self):
        self.assertIsNone(figures.number("inf"))
        self.assertIsNone(figures.number("nan"))


class TestOneRowPerVariant(unittest.TestCase):
    def test_a_variant_in_two_classes_counts_once(self):
        rows = [candidate("v1"), candidate("v1"), candidate("v2")]
        self.assertEqual(len(figures.one_row_per_variant(rows)), 2)

    def test_rows_without_an_id_are_all_kept(self):
        rows = [{"ID": ""}, {"ID": ""}]
        self.assertEqual(len(figures.one_row_per_variant(rows)), 2)


class TestLandscapePoints(unittest.TestCase):
    def test_a_complete_candidate_becomes_a_point(self):
        points, dropped = figures.landscape_points([candidate("v1")])
        self.assertEqual(dropped, 0)
        self.assertEqual(points[0]["ccf"], 0.9)
        self.assertEqual(points[0]["agretopicity"], 0.05)
        self.assertEqual(points[0]["expression"], 20.0)

    def test_missing_either_axis_is_counted_not_silently_dropped(self):
        points, dropped = figures.landscape_points([
            candidate("v1", CCF="NA"),
            candidate("v2", Agretopicity="NA"),
            candidate("v3"),
        ])
        self.assertEqual(len(points), 1)
        self.assertEqual(dropped, 2)

    def test_zero_agretopicity_cannot_go_on_a_log_axis(self):
        points, dropped = figures.landscape_points([candidate("v1", Agretopicity="0")])
        self.assertEqual((len(points), dropped), (0, 1))

    def test_expression_falls_back_to_the_pvacseq_column(self):
        points, _ = figures.landscape_points([candidate("v1", **{"Sample TPM": "NA", "RNA Expr": "7"})])
        self.assertEqual(points[0]["expression"], 7.0)

    def test_expression_may_be_missing_entirely(self):
        points, _ = figures.landscape_points([candidate("v1", **{"Sample TPM": "NA", "RNA Expr": "NA"})])
        self.assertIsNone(points[0]["expression"])

    def test_selection_defaults_to_rejected_when_blank(self):
        points, _ = figures.landscape_points([candidate("v1", Selection="")])
        self.assertEqual(points[0]["selection"], "Rejected")

    def test_each_variant_appears_once(self):
        points, _ = figures.landscape_points([candidate("v1"), candidate("v1")])
        self.assertEqual(len(points), 1)


class TestFlagCounts(unittest.TestCase):
    def test_counts_and_order(self):
        rows = [candidate("v1", Flags="SUBCLONAL;SELF_LIKE"),
                candidate("v2", Flags="SUBCLONAL"),
                candidate("v3", Flags="")]
        counts = figures.flag_counts(rows)
        self.assertEqual(counts[0][:2], ("SUBCLONAL", 2))
        self.assertEqual(counts[1][:2], ("SELF_LIKE", 1))

    def test_rejections_are_marked_apart_from_cautions(self):
        counts = dict((flag, rejects) for flag, _, rejects in
                      figures.flag_counts([candidate("v1", Flags="HLA_LOST;SUBCLONAL")]))
        self.assertTrue(counts["HLA_LOST"])
        self.assertFalse(counts["SUBCLONAL"])

    def test_every_reject_flag_is_recognised(self):
        flags = ";".join(figures.REJECT_FLAGS)
        counts = figures.flag_counts([candidate("v1", Flags=flags)])
        self.assertTrue(all(rejects for _, _, rejects in counts))

    def test_no_flags_at_all(self):
        self.assertEqual(figures.flag_counts([candidate("v1")]), [])

    def test_a_variant_in_two_classes_is_counted_once(self):
        rows = [candidate("v1", Flags="SUBCLONAL"), candidate("v1", Flags="SUBCLONAL")]
        self.assertEqual(figures.flag_counts(rows)[0][1], 1)


class TestLohEntries(unittest.TestCase):
    def rows(self):
        return [
            {"locus": "DRB1", "call": "not_evaluable", "allele1_pvacseq": "DRB1*15:01",
             "allele2_pvacseq": "DRB1*15:01", "lost_allele_pvacseq": ""},
            {"locus": "B", "call": "retained", "allele1_pvacseq": "HLA-B*07:02",
             "allele2_pvacseq": "HLA-B*08:01", "lost_allele_pvacseq": ""},
            {"locus": "A", "call": "loss", "allele1_pvacseq": "HLA-A*02:01",
             "allele2_pvacseq": "HLA-A*01:01", "lost_allele_pvacseq": "HLA-A*02:01"},
        ]

    def test_class_i_comes_first(self):
        entries = figures.loh_entries(self.rows())
        self.assertEqual([locus for locus, _, _ in entries], ["A", "B", "DRB1"])

    def test_a_loss_names_the_lost_allele(self):
        entries = figures.loh_entries(self.rows())
        self.assertIn("lost HLA-A*02:01", entries[0][2])

    def test_a_retained_locus_does_not(self):
        entries = figures.loh_entries(self.rows())
        self.assertNotIn("lost", entries[1][2])

    def test_rows_without_a_locus_are_skipped(self):
        self.assertEqual(figures.loh_entries([{"locus": "", "call": "loss"}]), [])

    def test_a_blank_call_is_treated_as_not_evaluable(self):
        entries = figures.loh_entries([{"locus": "A", "call": "", "allele1_pvacseq": "HLA-A*02:01"}])
        self.assertEqual(entries[0][1], "not_evaluable")


class TestReadTsv(unittest.TestCase):
    def test_a_missing_file_is_empty_not_an_error(self):
        self.assertEqual(figures.read_tsv(Path(tempfile.mkdtemp()) / "nope.tsv"), [])

    def test_no_path_at_all(self):
        self.assertEqual(figures.read_tsv(None), [])

    def test_an_empty_file(self):
        path = Path(tempfile.mkdtemp()) / "empty.tsv"
        path.write_text("")
        self.assertEqual(figures.read_tsv(path), [])

    def test_csv_is_read_too(self):
        path = Path(tempfile.mkdtemp()) / "t.csv"
        path.write_text("locus,call\nA,loss\n")
        self.assertEqual(figures.read_tsv(path)[0]["call"], "loss")


@unittest.skipUnless(HAS_MATPLOTLIB, "matplotlib is not installed")
class TestRendering(unittest.TestCase):
    def setUp(self):
        self.dir = Path(tempfile.mkdtemp())
        self.prefix = str(self.dir / "P1")
        self.candidates = write([candidate("v1"), candidate("v2", Selection="Rejected",
                                                            Flags="LOW_BINDING;SUBCLONAL")],
                                CANDIDATE_COLUMNS, "cand.tsv", self.dir)
        self.loh = write([{"patient": "P1", "locus": "A", "hla_class": "I",
                           "allele1_pvacseq": "HLA-A*02:01", "allele2_pvacseq": "HLA-A*01:01",
                           "call": "loss", "lost_allele_pvacseq": "HLA-A*02:01"}],
                         LOH_COLUMNS, "loh.tsv", self.dir)

    def run_main(self, *extra):
        return figures.main(["--prefix", self.prefix, *extra])

    def test_all_three_figures_are_written(self):
        self.assertEqual(self.run_main("--prioritized", str(self.candidates),
                                       "--hla-loh", str(self.loh)), 0)
        for name in figures.FIGURES:
            for suffix in ("png", "pdf"):
                path = Path(f"{self.prefix}.{name}.{suffix}")
                self.assertTrue(path.is_file(), path)
                self.assertGreater(path.stat().st_size, 1000, f"{path} looks empty")

    def test_only_selects_a_subset(self):
        self.run_main("--prioritized", str(self.candidates), "--only", "flags")
        self.assertTrue(Path(f"{self.prefix}.flags.png").is_file())
        self.assertFalse(Path(f"{self.prefix}.candidates.png").is_file())

    def test_format_selects_the_file_types(self):
        self.run_main("--prioritized", str(self.candidates), "--only", "flags", "--format", "png")
        self.assertTrue(Path(f"{self.prefix}.flags.png").is_file())
        self.assertFalse(Path(f"{self.prefix}.flags.pdf").is_file())

    def test_an_unknown_figure_name_is_an_error(self):
        with self.assertRaises(SystemExit):
            self.run_main("--prioritized", str(self.candidates), "--only", "nope")

    def test_no_loh_table_means_no_loh_figure(self):
        self.run_main("--prioritized", str(self.candidates))
        self.assertFalse(Path(f"{self.prefix}.hla_loh.png").is_file())

    def test_an_empty_candidate_table_still_draws_something(self):
        empty = write([], CANDIDATE_COLUMNS, "empty.tsv", self.dir)
        self.assertEqual(self.run_main("--prioritized", str(empty)), 0)
        self.assertTrue(Path(f"{self.prefix}.candidates.png").is_file())
        self.assertTrue(Path(f"{self.prefix}.flags.png").is_file())

    def test_candidates_with_no_plottable_axes_still_draw(self):
        unplottable = write([candidate("v1", CCF="NA", Agretopicity="NA")],
                            CANDIDATE_COLUMNS, "bad.tsv", self.dir)
        self.assertEqual(self.run_main("--prioritized", str(unplottable), "--only", "candidates"), 0)
        self.assertTrue(Path(f"{self.prefix}.candidates.png").is_file())

    def test_nothing_requested_writes_nothing_without_failing(self):
        self.assertEqual(self.run_main("--only", "candidates"), 0)
        self.assertFalse(Path(f"{self.prefix}.candidates.png").is_file())


if __name__ == "__main__":
    unittest.main()
