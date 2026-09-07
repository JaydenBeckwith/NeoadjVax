#!/usr/bin/env python3
"""Exercise the pinned upstream clusterer and our adapter on real BED12 junctions."""
import argparse
import csv
import importlib.util
from pathlib import Path
import subprocess
import sys
import tempfile

parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("cluster_script")
args = parser.parse_args()
clusterer = Path(args.cluster_script).resolve()
spec = importlib.util.spec_from_file_location("adapter", Path(__file__).resolve().parents[2] / "bin/prepare_splicing_inputs.py")
adapter = importlib.util.module_from_spec(spec)
spec.loader.exec_module(adapter)
with tempfile.TemporaryDirectory(prefix="leafcutter-test-") as scratch:
    root = Path(scratch)
    bed = root / "sample.junc"
    lines = []
    for left, right, reads, strand in [(160, 295, 5, "+"), (160, 298, 5, "+"),
                                        (160, 301, 50, "+"), (966, 1101, 5, "-"), (960, 1101, 45, "-")]:
        start, end = left - 10, right - 1 + 10
        lines.append(f"chr1\t{start}\t{end}\tJUNC\t{reads}\t{strand}\t{start}\t{end}\t255,0,0\t2\t10,10\t0,{end-start-10}")
    bed.write_text("\n".join(lines) + "\n")
    (root / "junctions.txt").write_text("sample.junc\n")
    subprocess.run([sys.executable, str(clusterer), "-j", str(root / "junctions.txt"),
                    "-o", "sample", "-r", root.as_posix(), "-m", "3", "-p", "0.001", "-l", "500000"], check=True, cwd=root)
    adapter.leafcutter_sample(root / "sample_perind.counts.gz", "", root / "selected.gz", root / "evidence.tsv")
    with (root / "evidence.tsv").open() as handle:
        observed = {row["junc_id"]: int(row["rna_junction_reads"]) for row in csv.DictReader(handle, delimiter="\t")}
    assert observed["chr1:160-295:+"] == 5
    assert observed["chr1:966-1101:-"] == 5
    assert len(observed) == 5
    # This upstream version applies -m to EACH junction, not only cluster totals.
    subprocess.run([sys.executable, str(clusterer), "-j", str(root / "junctions.txt"),
                    "-o", "strict", "-r", root.as_posix(), "-m", "50"], check=True, cwd=root)
    adapter.leafcutter_sample(root / "strict_perind.counts.gz", "", root / "strict.gz", root / "strict.tsv")
    with (root / "strict.tsv").open() as handle:
        assert list(csv.DictReader(handle, delimiter="\t")) == []
print("PASS: upstream LeafCutter BED12 clustering, strand/coordinate conversion and count adapter")
