"""design_long_peptides.py: ranked candidates -> peptides you can order."""
import csv
import importlib.util
import tempfile
import unittest
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "design_long_peptides.py"
spec = importlib.util.spec_from_file_location("design_long_peptides", BIN)
design = importlib.util.module_from_spec(spec)
spec.loader.exec_module(design)

REPORT_COLUMNS = ["ID", "Gene", "AA Change", "Allele", "Best Peptide", "CCF",
                  "Agretopicity", "HLA LOH", "Priority Rank", "Selection", "Flags"]

# A 60-residue protein with a substitution at position 31 (0-based 30).
WT_PROTEIN = "AAAAAGGGGGSSSSSTTTTTVVVVVLLLLL" + "KKKKKDDDDDEEEEENNNNNQQQQQPPPPP"
MT_PROTEIN = WT_PROTEIN[:30] + "W" + WT_PROTEIN[31:]


def write_fasta(entries):
    path = Path(tempfile.mkdtemp()) / "protein.fasta"
    with path.open("w", newline="") as handle:
        for header, sequence in entries:
            handle.write(f">{header}\n{sequence}\n")
    return path


def write_report(rows, directory):
    path = Path(directory) / "report.tsv"
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=REPORT_COLUMNS, delimiter="\t",
                                lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for row in rows:
            writer.writerow({**{c: "" for c in REPORT_COLUMNS}, **row})
    return path


def run(rows, fasta_entries, **kwargs):
    fasta = write_fasta(fasta_entries)
    report = write_report(rows, fasta.parent)
    out = fasta.parent / "peptides.tsv"
    argv = [str(report), "--protein-fasta", str(fasta), "--patient", "P1", "-o", str(out)]
    for key, value in kwargs.items():
        flag = f"--{key.replace('_', '-')}"
        argv += [flag] if value is True else [flag, str(value)]
    design.main(argv)
    with out.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


class TestMutatedRegion(unittest.TestCase):
    def test_substitution(self):
        self.assertEqual(design.mutated_region("ABCDE", "ABXDE"), (2, 3, "substitution"))

    def test_two_adjacent_changes(self):
        self.assertEqual(design.mutated_region("ABCDE", "ABXYE"), (2, 4, "multi_residue"))

    def test_frameshift_runs_to_the_end(self):
        start, end, kind = design.mutated_region("ABCWXYZ", "ABCDE")
        self.assertEqual((start, end, kind), (3, 7, "frameshift_or_indel"))

    def test_identical_sequences(self):
        self.assertEqual(design.mutated_region("ABC", "ABC")[2], "none")

    def test_no_wild_type_partner(self):
        self.assertEqual(design.mutated_region("ABCDE", ""), (0, 5, "unknown"))


class TestWindow(unittest.TestCase):
    def test_centred_on_the_mutation(self):
        start, end = design.window("X" * 60, 30, 31, 25, False)
        self.assertEqual(end - start, 25)
        self.assertTrue(start <= 30 < end)
        self.assertEqual(start, 18)

    def test_pulled_inside_at_the_n_terminus(self):
        start, end = design.window("X" * 60, 1, 2, 25, False)
        self.assertEqual((start, end), (0, 25))

    def test_pulled_inside_at_the_c_terminus(self):
        start, end = design.window("X" * 60, 58, 59, 25, False)
        self.assertEqual((start, end), (35, 60))

    def test_short_protein_gives_a_short_window(self):
        start, end = design.window("X" * 12, 5, 6, 25, False)
        self.assertEqual((start, end), (0, 12))

    def test_frameshift_is_anchored_at_the_junction(self):
        start, end = design.window("X" * 60, 30, 60, 25, True)
        self.assertEqual((start, end), (30, 55))


class TestDesign(unittest.TestCase):
    def setUp(self):
        self.candidate = {
            "ID": "chr1-100-100-A-T", "Gene": "GENE1", "AA Change": "L31W",
            "Allele": "HLA-A*02:01", "Best Peptide": MT_PROTEIN[26:35],
            "CCF": "0.9", "Priority Rank": "1", "Selection": "Top ranking",
        }
        self.fasta = [("MT.1.GENE1.ENST001.missense", MT_PROTEIN),
                      ("WT.1.GENE1.ENST001.missense", WT_PROTEIN)]

    def test_one_peptide_per_candidate(self):
        rows = run([self.candidate], self.fasta)
        self.assertEqual(len(rows), 1)

    def test_peptide_length_and_mutation_position(self):
        row = run([self.candidate], self.fasta)[0]
        self.assertEqual(int(row["Length"]), 25)
        position = int(row["Mutation Position"])
        self.assertEqual(row["Long Peptide"][position - 1], "W")

    def test_peptide_is_a_slice_of_the_mutant_protein(self):
        row = run([self.candidate], self.fasta)[0]
        self.assertIn(row["Long Peptide"], MT_PROTEIN)

    def test_it_contains_the_predicted_epitope(self):
        row = run([self.candidate], self.fasta)[0]
        self.assertEqual(row["Contains Best Peptide"], "yes")
        self.assertIn(self.candidate["Best Peptide"], row["Long Peptide"])

    def test_max_length_is_respected(self):
        row = run([self.candidate], self.fasta, max_length=20)[0]
        self.assertEqual(int(row["Length"]), 20)

    def test_peptide_id_and_patient(self):
        row = run([self.candidate], self.fasta)[0]
        self.assertEqual(row["Patient"], "P1")
        self.assertEqual(row["Peptide ID"], "P1_LP01_GENE1")

    def test_candidate_columns_are_carried_through(self):
        row = run([self.candidate], self.fasta)[0]
        self.assertEqual(row["Allele"], "HLA-A*02:01")
        self.assertEqual(row["CCF"], "0.9")
        self.assertEqual(row["Priority Rank"], "1")

    def test_rejected_candidates_are_skipped(self):
        rejected = {**self.candidate, "Selection": "Rejected"}
        self.assertEqual(run([rejected], self.fasta), [])

    def test_rejected_can_be_included(self):
        rejected = {**self.candidate, "Selection": "Rejected"}
        self.assertEqual(len(run([rejected], self.fasta, include_rejected=True)), 1)

    def test_top_limits_the_output(self):
        second = {**self.candidate, "ID": "chr2-1-1-A-T", "AA Change": "L31W_b"}
        rows = run([self.candidate, second], self.fasta, top=1)
        self.assertEqual(len(rows), 1)

    def test_candidate_with_no_protein_sequence_is_skipped(self):
        orphan = {**self.candidate, "Gene": "ABSENT", "Best Peptide": "QQQQQQQQQ"}
        self.assertEqual(run([orphan], [("MT.9.OTHER.x.y", "MMMMMMMMMMMM")]), [])

    def test_frameshift_peptide_starts_at_the_junction(self):
        # A frameshift changes the length as well as the residues.
        mt = WT_PROTEIN[:30] + "W" * 35
        candidate = {**self.candidate, "AA Change": "L31fs", "Best Peptide": mt[30:39]}
        row = run([candidate], [("MT.2.GENE1.x.fs", mt), ("WT.2.GENE1.x.fs", WT_PROTEIN)])[0]
        self.assertEqual(row["Mutation Type"], "frameshift_or_indel")
        self.assertEqual(int(row["Mutation Position"]), 1)
        self.assertTrue(row["Long Peptide"].startswith("W"))


class TestMerging(unittest.TestCase):
    def test_two_mutations_in_one_window_give_one_peptide(self):
        mt = WT_PROTEIN[:30] + "W" + WT_PROTEIN[31:33] + "W" + WT_PROTEIN[34:]
        fasta = [("MT.1.GENE1.x.y", mt), ("WT.1.GENE1.x.y", WT_PROTEIN)]
        base = {"Gene": "GENE1", "Allele": "HLA-A*02:01", "Selection": "Top ranking"}
        rows = run([
            {**base, "AA Change": "L31W", "Best Peptide": mt[26:35], "Priority Rank": "1"},
            {**base, "AA Change": "T34W", "Best Peptide": mt[30:39], "Priority Rank": "2"},
        ], fasta)
        self.assertEqual(len(rows), 1)
        self.assertEqual(rows[0]["Covers"], "L31W;T34W")

    def test_distant_mutations_give_two_peptides(self):
        mt = "W" + WT_PROTEIN[1:59] + "W"
        fasta = [("MT.1.GENE1.x.y", mt), ("WT.1.GENE1.x.y", WT_PROTEIN)]
        # Two separate entries, since one fasta pair can only carry one
        # difference pattern; this is the usual case of two transcripts.
        fasta2 = [("MT.1.GENE1.x.y", "W" + WT_PROTEIN[1:]),
                  ("WT.1.GENE1.x.y", WT_PROTEIN),
                  ("MT.2.GENE1.x.y", WT_PROTEIN[:59] + "W"),
                  ("WT.2.GENE1.x.y", WT_PROTEIN)]
        base = {"Gene": "GENE1", "Selection": "Top ranking"}
        rows = run([
            {**base, "AA Change": "A1W", "Best Peptide": ("W" + WT_PROTEIN[1:])[0:9]},
            {**base, "AA Change": "P60W", "Best Peptide": (WT_PROTEIN[:59] + "W")[51:60]},
        ], fasta2)
        self.assertEqual(len(rows), 2)


class TestTrimming(unittest.TestCase):
    def test_bad_nterm_residue_is_trimmed(self):
        # Mutation in the middle, cysteine at the window's N terminus.
        wt = "AAAAAAAAAACAAAAAAAAAAAAAAAAAAAAAAAAAAAAA"
        mt = wt[:22] + "W" + wt[23:]
        fasta = [("MT.1.GENE1.x.y", mt), ("WT.1.GENE1.x.y", wt)]
        candidate = {"Gene": "GENE1", "Best Peptide": mt[18:27], "Selection": "Top ranking",
                     "AA Change": "A23W"}
        row = run([candidate], fasta, max_length=13, min_length=11, bad_nterm="C")[0]
        self.assertFalse(row["Long Peptide"].startswith("C"))

    def test_trimming_never_cuts_past_the_mutation(self):
        # Every residue is "bad", so trimming must stop at the min flank.
        wt = "C" * 40
        mt = wt[:20] + "W" + wt[21:]
        fasta = [("MT.1.GENE1.x.y", mt), ("WT.1.GENE1.x.y", wt)]
        candidate = {"Gene": "GENE1", "Best Peptide": mt[16:25], "Selection": "Top ranking"}
        row = run([candidate], fasta, bad_nterm="C", bad_cterm="C", min_length=15, min_flank=3)[0]
        self.assertGreaterEqual(int(row["Length"]), 15)
        self.assertIn("W", row["Long Peptide"])


class TestNotes(unittest.TestCase):
    def test_free_cysteine_is_warned_about(self):
        _, warnings = design.peptide_notes("AAACAAA", 1.0, 4)
        self.assertTrue(any("free cysteine" in w for w in warnings))

    def test_two_cysteines_warn_about_disulfides(self):
        _, warnings = design.peptide_notes("ACAACAA", 1.0, 4)
        self.assertTrue(any("disulfide" in w for w in warnings))

    def test_hydrophobic_peptide_is_warned_about(self):
        metrics, warnings = design.peptide_notes("IIIIVVVVLLLL", 1.0, 99)
        self.assertGreater(float(metrics["GRAVY"]), 1.0)
        self.assertTrue(any("solubility" in w for w in warnings))

    def test_net_charge(self):
        metrics, _ = design.peptide_notes("KKRDE", 1.0, 99)
        self.assertEqual(metrics["Net Charge"], 1)

    def test_longest_run(self):
        metrics, warnings = design.peptide_notes("AGGGGSA", 1.0, 4)
        self.assertEqual(metrics["Longest Run"], 4)
        self.assertTrue(any("run" in w for w in warnings))

    def test_clean_peptide_has_no_warnings(self):
        metrics, warnings = design.peptide_notes("SLYNTVATLGSAK", 1.0, 4)
        self.assertEqual(warnings, [])
        self.assertEqual(metrics["Cys"], 0)


class TestFastaOutput(unittest.TestCase):
    def test_fasta_is_written(self):
        fasta = write_fasta([("MT.1.GENE1.x.y", MT_PROTEIN), ("WT.1.GENE1.x.y", WT_PROTEIN)])
        report = write_report([{"Gene": "GENE1", "AA Change": "L31W",
                                "Best Peptide": MT_PROTEIN[26:35], "Selection": "Top ranking"}],
                              fasta.parent)
        out = fasta.parent / "peptides.tsv"
        fasta_out = fasta.parent / "peptides.fasta"
        design.main([str(report), "--protein-fasta", str(fasta), "--patient", "P1",
                     "-o", str(out), "--fasta-out", str(fasta_out)])
        entries = design.read_fasta(fasta_out)
        self.assertEqual(len(entries), 1)
        self.assertEqual(entries[0][0], "P1_LP01_GENE1")
        self.assertEqual(len(entries[0][1]), 25)


class TestFastaParsing(unittest.TestCase):
    def test_wrapped_sequences_are_joined(self):
        path = Path(tempfile.mkdtemp()) / "wrapped.fasta"
        path.write_text(">MT.1.G.x.y\nAAAA\nBBBB\n>WT.1.G.x.y\nAAAA\nCBBB\n")
        pairs = design.pair_entries(design.read_fasta(path))
        self.assertEqual(pairs[0]["mt"], "AAAABBBB")
        self.assertEqual(pairs[0]["wt"], "AAAACBBB")

    def test_entries_are_paired_by_index(self):
        entries = [("MT.7.GENE1.t.c", "AAAW"), ("WT.7.GENE1.t.c", "AAAA")]
        pairs = design.pair_entries(entries)
        self.assertEqual(len(pairs), 1)
        self.assertEqual(pairs[0]["gene"], "GENE1")

    def test_a_plain_fasta_is_treated_as_mutant_only(self):
        pairs = design.pair_entries([("GENE1_peptide", "AAAW")])
        self.assertEqual(pairs[0]["mt"], "AAAW")
        self.assertEqual(pairs[0]["wt"], "")


if __name__ == "__main__":
    unittest.main()
