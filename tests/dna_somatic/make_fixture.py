#!/usr/bin/env python3
"""Generate deterministic synthetic tumour/normal DNA inputs; no patient data.

P1 enters as already-aligned BAMs with HLA alleles in the samplesheet and a
Sequenza segments file; P2 enters as FASTQ pairs with no HLA alleles, so it
goes through alignment and HLA-HD. Reference/resource files are empty
placeholders: only valid with nextflow -stub-run.
"""
import argparse
import csv
import gzip
import json
from pathlib import Path
import random
import shutil
import subprocess

PATIENTS = ("P1", "P2")
PAIRS = 24
SCATTER = ("chr1",)


def write_csv(path, fields, rows):
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def read_pairs(sample):
    """(name, start, mate) for each synthetic read pair."""
    return [(f"{sample}_pair{n:03d}", 200 + n * 7, 350 + n * 7) for n in range(PAIRS)]


def write_sam(path, reference, sample):
    """SAM SEQ is in alignment orientation, including reverse-strand mates."""
    records = []
    for name, start, mate in read_pairs(sample):
        span = mate + 100 - start
        records += [
            (start, f"{name}\t99\tchr1\t{start + 1}\t60\t100M\t=\t{mate + 1}\t{span}\t{reference[start:start + 100]}\t{'I' * 100}\tRG:Z:{sample}"),
            (mate, f"{name}\t147\tchr1\t{mate + 1}\t60\t100M\t=\t{start + 1}\t{-span}\t{reference[mate:mate + 100]}\t{'I' * 100}\tRG:Z:{sample}"),
        ]
    header = (f"@HD\tVN:1.6\tSO:coordinate\n@SQ\tSN:chr1\tLN:{len(reference)}\n"
              f"@RG\tID:{sample}\tSM:{sample}\tLB:synthetic\tPL:ILLUMINA\tPU:CI\n")
    path.write_text(header + "\n".join(record for _, record in sorted(records)) + "\n", encoding="ascii")


def write_fastq(r1, r2, reference, sample):
    complement = str.maketrans("ACGT", "TGCA")
    with gzip.open(r1, "wt") as one, gzip.open(r2, "wt") as two:
        for name, start, mate in read_pairs(sample):
            one.write(f"@{name}/1\n{reference[start:start + 100]}\n+\n{'I' * 100}\n")
            two.write(f"@{name}/2\n{reference[mate:mate + 100].translate(complement)[::-1]}\n+\n{'I' * 100}\n")


def generate(root, samtools="samtools", sam_only=False):
    root = Path(root).resolve()
    root.mkdir(parents=True, exist_ok=True)
    if any(root.iterdir()):
        raise ValueError(f"Fixture output directory must be empty: {root}")
    if not sam_only and not shutil.which(samtools):
        raise ValueError("Install samtools or pass --samtools /path/to/samtools")
    rng = random.Random(104729)
    reference = "".join(rng.choice("ACGT") for _ in range(10000))
    (root / "genome.fa").write_text(">chr1\n" + "\n".join(reference[i:i + 80] for i in range(0, len(reference), 80)) + "\n")
    for name in ("known_dbsnp", "known_mills", "known_1000g", "germline_resource", "panel_of_normals", "contamination_sites"):
        (root / f"{name}.vcf.gz").write_bytes(b"")
        (root / f"{name}.vcf.gz.tbi").write_bytes(b"")
    for name in ("vep_cache", "vep_plugins"):
        (root / name).mkdir()
        (root / name / "STUB_ONLY.txt").write_text("Not a reference cache. Used only with nextflow -stub-run.\n")

    rows, inputs = [], []
    for patient in PATIENTS:
        row = {"patient_id": patient}
        for sample_type in ("tumor", "normal"):
            sample = f"{patient}_{sample_type}"
            if patient == "P1":
                sam, bam = root / f"{sample}.sam", root / f"{sample}.bam"
                write_sam(sam, reference, sample)
                if not sam_only:
                    subprocess.run([samtools, "sort", "-o", str(bam), str(sam)], check=True)
                    subprocess.run([samtools, "index", str(bam)], check=True)
                    subprocess.run([samtools, "quickcheck", "-v", str(bam)], check=True)
                    count = int(subprocess.check_output([samtools, "view", "-c", str(bam)], text=True))
                    if count != PAIRS * 2:
                        raise ValueError(f"{bam}: expected {PAIRS * 2} reads, got {count}")
                row[f"{sample_type}_bam"] = str(bam)
                inputs.append({"path": str(bam), "sample": sample, "kind": "bam", "reads": PAIRS * 2})
            else:
                r1, r2 = root / f"{sample}_R1.fastq.gz", root / f"{sample}_R2.fastq.gz"
                write_fastq(r1, r2, reference, sample)
                row[f"{sample_type}_r1"], row[f"{sample_type}_r2"] = str(r1), str(r2)
                inputs.append({"path": str(r1), "sample": sample, "kind": "fastq", "reads": PAIRS})
        if patient == "P1":
            segments = root / "P1_segments.txt"
            segments.write_text('"chromosome"\t"start.pos"\t"end.pos"\t"CNt"\t"A"\t"B"\n"chr1"\t1\t10000\t3\t2\t1\n')
            row.update(purity="0.62", ploidy="3.1", sequenza_segments=str(segments), sex="female",
                       hla_alleles="HLA-A*02:01|HLA-B*07:02|DRB1*15:01")
        else:
            row.update(purity="0.45", ploidy="", sequenza_segments="", sex="", hla_alleles="")
        rows.append(row)

    fields = ["patient_id", "tumor_r1", "tumor_r2", "normal_r1", "normal_r2", "tumor_bam", "normal_bam",
              "purity", "ploidy", "sequenza_segments", "sex", "hla_alleles"]
    write_csv(root / "samplesheet.csv", fields, rows)
    manifest = {"synthetic": True, "test_type": "wiring_only", "patients": list(PATIENTS),
                "fastq_patients": ["P2"], "hla_typed_patients": ["P2"], "scatter": list(SCATTER),
                "inputs": inputs, "bam_generation_complete": not sam_only}
    (root / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n")
    return manifest


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--outdir", type=Path, required=True)
    parser.add_argument("--samtools", default="samtools")
    parser.add_argument("--sam-only", action="store_true", help="Developer inspection only; CI always builds/validates BAMs")
    args = parser.parse_args()
    result = generate(args.outdir, args.samtools, args.sam_only)
    print(f"Generated {len(result['inputs'])} synthetic read sets in {args.outdir}; BAMs complete: {result['bam_generation_complete']}")
