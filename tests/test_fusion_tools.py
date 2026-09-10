import argparse
from contextlib import closing
import csv
import importlib.util
import json
from pathlib import Path
import sqlite3
import tempfile
import unittest
from unittest.mock import patch


SPEC = importlib.util.spec_from_file_location("fusion", Path(__file__).resolve().parents[1] / "bin/fusion_tools.py")
fusion = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(fusion)

STAR = {"#FusionName": "GENE1--GENE2", "LeftGene": "GENE1^ENSG00000000001.1",
        "RightGene": "GENE2^ENSG00000000002.2", "LeftBreakpoint": "chr1:123:+",
        "RightBreakpoint": "chr2:456:-", "JunctionReadCount": "3", "SpanningFragCount": "2", "FFPM": "0.1"}
ARRIBA = {"#gene1": "GENE1", "gene2": "GENE2", "gene_id1": "ENSG00000000001.1",
          "gene_id2": "ENSG00000000002.2", "breakpoint1": "1:123", "breakpoint2": "2:456",
          "split_reads1": "2", "split_reads2": "1", "discordant_mates": "2"}


class FusionTools(unittest.TestCase):
    def setUp(self):
        self.temp = tempfile.TemporaryDirectory()
        self.addCleanup(self.temp.cleanup)
        self.root = Path(self.temp.name)

    def calls(self, rows, template=STAR):
        path = self.root / "caller.tsv"
        with path.open("w", newline="") as handle:
            writer = csv.DictWriter(handle, delimiter="\t", fieldnames=list(template))
            writer.writeheader()
            writer.writerows(rows)
        return path

    def database(self, release=111):
        path = self.root / f"agfusion.homo_sapiens.{release}.db"
        with closing(sqlite3.connect(path)) as db:
            db.execute(f"CREATE TABLE homo_sapiens_{release} (stable_id TEXT)")
            db.execute(f"CREATE TABLE homo_sapiens_{release}_transcript (transcript_stable_id TEXT)")
            db.commit()
        return path

    def annotate_args(self, rows):
        cache = self.root / "cache"
        cache.mkdir()
        (cache / "reference.db").touch()
        return argparse.Namespace(input=self.calls(rows), output=self.root / "agfusion_out",
                                  caller="starfusion", min_reads=5, min_ffpm=0.1,
                                  database=self.database(), cache=cache, noncanonical=False)

    def protein(self, directory, seq="MAAAAA*VVVVVVVVVVVV"):
        event = directory / "GENE1-123_GENE2-456"
        event.mkdir(parents=True)
        (event / "GENE1-GENE2_protein.fa").write_text(">test\n" + seq + "\n")
        (event / "GENE1-GENE2.exons.txt").write_text("exon_chr\n1\n")

    def predict_args(self, caller="starfusion"):
        models = self.root / "models"
        models.mkdir(exist_ok=True)
        (models / "manifest.csv").touch()
        return argparse.Namespace(input=self.root / "agfusion_out", output=self.root / "pvacfuse_output",
                                  sample="P1_D42_" + caller, caller=caller,
                                  alleles="HLA-A*02:01|HLA-B*07:02", algorithms="MHCflurry",
                                  threads=4, lengths_i="8,9,10,11", lengths_ii="15",
                                  min_reads=5, min_ffpm=0.1, binding_threshold=500, models=models, iedb=None)

    def test_native_starfusion_counts_and_columns(self):
        fields, rows, audit = fusion.read_calls(self.calls([STAR]), "starfusion", 5, 0.1)
        self.assertEqual(rows, [STAR])
        self.assertIn("#FusionName", fields)
        self.assertEqual(audit[0]["read_support"], 5)
        self.assertEqual(audit[0]["decision"], "retained")

    def test_arriba_counts_are_not_fabricated_ffpm(self):
        _, rows, audit = fusion.read_calls(self.calls([ARRIBA], ARRIBA), "arriba", 5, 0.1)
        self.assertEqual(rows, [ARRIBA])
        self.assertEqual(audit[0]["read_support"], 5)
        self.assertEqual(audit[0]["ffpm"], "NA")

    def test_filter_keeps_explicit_reasons(self):
        row = dict(STAR, JunctionReadCount="0", FFPM="0.01")
        _, rows, audit = fusion.read_calls(self.calls([row]), "starfusion", 5, 0.1)
        self.assertFalse(rows)
        self.assertEqual(audit[0]["decision"], "low_read_support;low_ffpm")

    def test_header_only_is_valid_empty_result(self):
        _, rows, audit = fusion.read_calls(self.calls([]), "starfusion", 5, 0.1)
        self.assertEqual((rows, audit), ([], []))

    def test_zero_byte_file_is_not_valid_empty_result(self):
        source = self.root / "empty.tsv"
        source.touch()
        with self.assertRaises(ValueError):
            fusion.read_calls(source, "starfusion", 5, 0.1)

    def test_malformed_rows_and_missing_columns_fail(self):
        for content in ("#FusionName\tLeftGene\nx\ty\n", "\t".join(STAR) + "\nshort\n"):
            source = self.root / "malformed.tsv"
            source.write_text(content)
            with self.assertRaises(ValueError):
                fusion.read_calls(source, "starfusion", 5, 0.1)

    def test_bad_breakpoints_gene_ids_and_counts_fail(self):
        for field, value in (("LeftBreakpoint", "chr1:0:+"), ("LeftBreakpoint", "chr1:10"),
                             ("LeftGene", "GENE1"), ("LeftGene", "GENE1^ENSG1,ENSG2"),
                             ("JunctionReadCount", "-1"), ("JunctionReadCount", "1.5"),
                             ("FFPM", "nan"), ("FFPM", "inf")):
            with self.subTest(field=field, value=value), self.assertRaises(ValueError):
                fusion.read_calls(self.calls([dict(STAR, **{field: value})]), "starfusion", 5, 0.1)

    def test_empty_annotation_does_not_launch_agfusion(self):
        args = self.annotate_args([])
        with patch.object(fusion, "run_logged") as runner:
            qc = fusion.annotate(args)
        runner.assert_not_called()
        self.assertEqual(qc["status"], "no_fusions")
        self.assertEqual(qc["protein_records"], 0)
        self.assertTrue((args.output / "input_audit.tsv").exists())

    def test_annotation_command_and_marker(self):
        args = self.annotate_args([STAR])
        args.noncanonical = True
        def run(command, log, env):
            log.write_text("INFO annotated\n")
            self.protein(args.output)
            self.assertEqual(command[:2], ["agfusion", "batch"])
            self.assertEqual(command[command.index("-a") + 1], "starfusion")
            self.assertIn("--middlestar", command)
            self.assertIn("--noncanonical", command)
            self.assertEqual(env["PYENSEMBL_CACHE_DIR"], str(args.cache.resolve()))
        with patch.object(fusion, "run_logged", side_effect=run):
            self.assertEqual(fusion.annotate(args)["protein_records"], 1)

    def test_agfusion_zero_exit_with_logged_error_is_failure(self):
        args = self.annotate_args([STAR])
        with patch.object(fusion, "run_logged", side_effect=lambda cmd, log, env: log.write_text("ERROR gene not found")):
            with self.assertRaisesRegex(RuntimeError, "annotation errors"):
                fusion.annotate(args)

    def test_missing_or_multiple_fusion_markers_fail(self):
        for seq in ("AAAAAVVVVV", "AAAA*VVV*"):
            with self.subTest(sequence=seq):
                directory = self.root / str(len(seq))
                self.protein(directory, seq)
                with self.assertRaisesRegex(ValueError, "middlestar"):
                    fusion.protein_records(directory)

    def test_missing_exons_fail(self):
        directory = self.root / "annotations"
        self.protein(directory)
        next(directory.glob("*/*.exons.txt")).unlink()
        with self.assertRaisesRegex(ValueError, "exon table"):
            fusion.protein_records(directory)

    def test_database_filename_and_internal_release_agree(self):
        db = self.database(111)
        self.assertEqual(fusion.check_database(db), 111)
        wrong = db.with_name("agfusion.homo_sapiens.110.db")
        db.rename(wrong)
        with self.assertRaisesRegex(ValueError, "release 110"):
            fusion.check_database(wrong)

    def test_hla_cannot_be_null_or_shell_content(self):
        self.assertEqual(fusion.normalize_alleles("HLA-A*02:01|HLA-A*02:01,HLA-B*07:02"), ["HLA-A*02:01", "HLA-B*07:02"])
        for value in ("", "  ", "HLA-A*02:01'; bad"):
            with self.assertRaises(ValueError):
                fusion.normalize_alleles(value)

    def test_pvacfuse_native_starfusion_evidence_is_forwarded(self):
        args = self.predict_args()
        cmd, _, _ = fusion.prediction_command(args)
        self.assertEqual(cmd[4], "HLA-A*02:01,HLA-B*07:02")
        self.assertIn("--starfusion-file", cmd)
        self.assertEqual(cmd[cmd.index("--starfusion-file") + 1], str(args.input / "filtered_fusions.tsv"))
        self.assertEqual(cmd[cmd.index("-t") + 1], "4")

    def test_arriba_never_masquerades_as_starfusion_evidence(self):
        cmd, _, _ = fusion.prediction_command(self.predict_args("arriba"))
        self.assertNotIn("--starfusion-file", cmd)

    def test_multiple_algorithms_are_separate_arguments_and_offline(self):
        args = self.predict_args()
        args.algorithms = "MHCflurry,NetMHCIIpan"
        with self.assertRaisesRegex(ValueError, "local IEDB"):
            fusion.prediction_command(args)
        args.iedb = self.root / "iedb"
        for relative in ("mhc_i/src/predict_binding.py", "mhc_ii/mhc_II_binding.py"):
            executable = args.iedb / relative
            executable.parent.mkdir(parents=True, exist_ok=True)
            executable.touch()
        command, _, _ = fusion.prediction_command(args)
        self.assertEqual(command[5:7], ["MHCflurry", "NetMHCIIpan"])
        self.assertNotIn("MHCflurry,NetMHCIIpan", command)

    def test_bad_predictors_or_missing_models_fail(self):
        args = self.predict_args()
        for value in ("", "unknown", "MHCflurry,MHCflurry"):
            args.algorithms = value
            with self.assertRaises(ValueError):
                fusion.prediction_command(args)
        args.algorithms = "MHCflurry"
        args.models = None
        with self.assertRaisesRegex(ValueError, "manifest.csv"):
            fusion.prediction_command(args)

    def test_empty_annotation_writes_status_not_fake_binding_report(self):
        args = self.predict_args()
        args.input.mkdir()
        (args.input / "annotation_qc.json").write_text('{"status":"no_fusions"}')
        with patch.object(fusion, "validate_prediction_pairs", return_value=([["MHCflurry", "HLA-A*02:01"]], [])), patch.object(fusion, "run_logged") as runner:
            result = fusion.predict(args)
        runner.assert_not_called()
        self.assertEqual(result["status"], "no_predictable_fusion_proteins")
        self.assertEqual(list(args.output.glob("*.tsv")), [])
        self.assertTrue((args.output / "run_status.json").is_file())

    def test_nonempty_proteins_without_report_fail(self):
        args = self.predict_args()
        self.protein(args.input)
        (args.input / "annotation_qc.json").write_text('{"status":"annotated"}')
        with patch.object(fusion, "validate_prediction_pairs", return_value=([["MHCflurry", "HLA-A*02:01"]], [])), patch.object(fusion, "run_logged"):
            with self.assertRaisesRegex(RuntimeError, "no all_epitopes report"):
                fusion.predict(args)


if __name__ == "__main__":
    unittest.main()
