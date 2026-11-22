"""clonal_structure.py: CCF clusters and week 6 persistence."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "clonal_structure.py"
spec = importlib.util.spec_from_file_location("clonal_structure", BIN)
clonal = importlib.util.module_from_spec(spec)
spec.loader.exec_module(clonal)

COLUMNS = ["ID", "Gene", "AA Change", "CCF", "RNA Depth", "RNA VAF",
           "RNA Depth (week 6)", "RNA VAF (week 6)", "Selection", "Priority Rank"]


def variant(id_, ccf, **kwargs):
    row = {c: "NA" for c in COLUMNS}
    row.update({"ID": id_, "Gene": f"G{id_}", "CCF": str(ccf), "Selection": "Top ranking"})
    row.update(kwargs)
    return row


def run(rows, **kwargs):
    directory = Path(tempfile.mkdtemp())
    report = directory / "report.tsv"
    with report.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)
    prefix = str(directory / "P1")
    argv = [str(report), "--prefix", prefix, "--no-plot"]
    for key, value in kwargs.items():
        argv += [f"--{key.replace('_', '-')}", str(value)]
    clonal.main(argv)

    def read(name):
        with open(f"{prefix}.{name}.tsv", newline="") as handle:
            return list(csv.DictReader(handle, delimiter="\t"))
    return read("clonal_structure"), read("clonal_variants")


class TestClustering(unittest.TestCase):
    def test_one_tight_group_is_one_cluster(self):
        self.assertEqual(len(clonal.cluster_ccfs([0.9, 0.92, 0.95, 0.88], 0.15, 3)), 1)

    def test_a_wide_gap_splits(self):
        clusters = clonal.cluster_ccfs([0.9, 0.92, 0.95, 0.3, 0.32, 0.28], 0.15, 3)
        self.assertEqual(len(clusters), 2)

    def test_highest_ccf_cluster_comes_first(self):
        clusters = clonal.cluster_ccfs([0.3, 0.32, 0.28, 0.9, 0.92, 0.95], 0.15, 3)
        self.assertGreater(sum(clusters[0]) / len(clusters[0]), sum(clusters[1]) / len(clusters[1]))

    def test_a_stray_variant_is_merged_into_its_nearest_neighbour(self):
        clusters = clonal.cluster_ccfs([0.9, 0.92, 0.95, 0.5], 0.15, 3)
        self.assertEqual(len(clusters), 1)
        self.assertIn(0.5, clusters[0])

    def test_min_cluster_of_one_keeps_everything(self):
        clusters = clonal.cluster_ccfs([0.9, 0.92, 0.5], 0.15, 1)
        self.assertEqual(len(clusters), 2)

    def test_empty_input(self):
        self.assertEqual(clonal.cluster_ccfs([], 0.15, 3), [])

    def test_gap_is_an_argument(self):
        self.assertEqual(len(clonal.cluster_ccfs([0.9, 0.7], 0.1, 1)), 2)
        self.assertEqual(len(clonal.cluster_ccfs([0.9, 0.7], 0.3, 1)), 1)


class TestRnaStatus(unittest.TestCase):
    def test_detected(self):
        self.assertEqual(clonal.rna_status(50, 0.3, 5), "detected")

    def test_covered_but_absent(self):
        self.assertEqual(clonal.rna_status(50, 0.0, 5), "not_detected")

    def test_too_shallow_is_no_coverage(self):
        self.assertEqual(clonal.rna_status(2, 0.0, 5), "no_coverage")

    def test_no_data_is_no_coverage(self):
        self.assertEqual(clonal.rna_status(None, None, 5), "no_coverage")


class TestChange(unittest.TestCase):
    def test_persisted(self):
        self.assertEqual(clonal.change("detected", "detected"), "persisted")

    def test_lost(self):
        self.assertEqual(clonal.change("detected", "not_detected"), "lost")

    def test_gained(self):
        self.assertEqual(clonal.change("not_detected", "detected"), "gained")

    def test_absent_throughout(self):
        self.assertEqual(clonal.change("not_detected", "not_detected"), "absent_throughout")

    def test_missing_coverage_is_never_a_change(self):
        self.assertEqual(clonal.change("detected", "no_coverage"), "not_evaluable")
        self.assertEqual(clonal.change("no_coverage", "detected"), "not_evaluable")


class TestTables(unittest.TestCase):
    def setUp(self):
        self.rows = [
            variant("v1", 0.95, **{"RNA Depth": "40", "RNA VAF": "0.3",
                                   "RNA Depth (week 6)": "40", "RNA VAF (week 6)": "0.25"}),
            variant("v2", 0.92, **{"RNA Depth": "40", "RNA VAF": "0.3",
                                   "RNA Depth (week 6)": "40", "RNA VAF (week 6)": "0"}),
            variant("v3", 0.90, **{"RNA Depth": "40", "RNA VAF": "0.3",
                                   "RNA Depth (week 6)": "40", "RNA VAF (week 6)": "0.2"}),
            variant("v4", 0.30, Selection="Rejected"),
            variant("v5", 0.28),
            variant("v6", 0.32),
        ]
        self.clusters, self.variants = run(self.rows)

    def test_two_clusters_found(self):
        self.assertEqual(len(self.clusters), 2)

    def test_clonal_and_subclonal_labels(self):
        self.assertEqual(self.clusters[0]["Label"], "clonal")
        self.assertEqual(self.clusters[1]["Label"], "subclonal")

    def test_variant_counts(self):
        self.assertEqual([int(c["Variants"]) for c in self.clusters], [3, 3])

    def test_rejected_candidates_are_counted_separately(self):
        self.assertEqual(int(self.clusters[1]["Candidates Kept"]), 2)

    def test_week6_persistence_is_a_fraction_of_evaluable_variants(self):
        self.assertEqual(self.clusters[0]["Week 6 Detected"], "2")
        self.assertEqual(self.clusters[0]["Week 6 Not Detected"], "1")
        self.assertAlmostEqual(float(self.clusters[0]["Week 6 Persistence"]), 2 / 3, places=2)

    def test_clusters_without_rna_report_na_persistence(self):
        self.assertEqual(self.clusters[1]["Week 6 Persistence"], "NA")

    def test_per_variant_table_has_a_row_each(self):
        self.assertEqual(len(self.variants), 6)

    def test_per_variant_change_column(self):
        changes = {row["ID"]: row["Change"] for row in self.variants}
        self.assertEqual(changes["v1"], "persisted")
        self.assertEqual(changes["v2"], "lost")
        self.assertEqual(changes["v4"], "not_evaluable")

    def test_variants_carry_their_cluster(self):
        clusters = {row["ID"]: row["Cluster"] for row in self.variants}
        self.assertEqual(clusters["v1"], "1")
        self.assertEqual(clusters["v5"], "2")

    def test_columns(self):
        for column in clonal.CLUSTER_COLUMNS:
            self.assertIn(column, self.clusters[0])
        for column in clonal.VARIANT_COLUMNS:
            self.assertIn(column, self.variants[0])


class TestEdges(unittest.TestCase):
    def test_one_variant_per_id_even_with_several_epitopes(self):
        rows = [variant("v1", 0.9), variant("v1", 0.9), variant("v2", 0.91)]
        clusters, variants = run(rows)
        self.assertEqual(len(variants), 2)
        self.assertEqual(int(clusters[0]["Variants"]), 2)

    def test_variants_without_a_ccf_are_kept_but_unclustered(self):
        rows = [variant("v1", 0.9), variant("v2", 0.91), variant("v3", 0.92),
                {**variant("v4", 0.0), "CCF": "NA"}]
        clusters, variants = run(rows)
        unclustered = [v for v in variants if v["Cluster"] == "NA"]
        self.assertEqual(len(unclustered), 1)
        self.assertEqual(unclustered[0]["Cluster Label"], "unclustered")

    def test_an_empty_report_produces_empty_tables(self):
        clusters, variants = run([])
        self.assertEqual(clusters, [])
        self.assertEqual(variants, [])

    def test_a_report_with_no_ccf_column_at_all(self):
        directory = Path(tempfile.mkdtemp())
        report = directory / "report.tsv"
        report.write_text("ID\tGene\nv1\tG1\n")
        prefix = str(directory / "P1")
        clonal.main([str(report), "--prefix", prefix, "--no-plot"])
        with open(f"{prefix}.clonal_variants.tsv", newline="") as handle:
            rows = list(csv.DictReader(handle, delimiter="\t"))
        self.assertEqual(rows[0]["CCF"], "NA")

    def test_clonal_threshold_is_an_argument(self):
        rows = [variant("v1", 0.9), variant("v2", 0.91), variant("v3", 0.92)]
        clusters, _ = run(rows, clonal_ccf=0.99)
        self.assertEqual(clusters[0]["Label"], "subclonal")


if __name__ == "__main__":
    unittest.main()
