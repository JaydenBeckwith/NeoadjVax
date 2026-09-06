#!/usr/bin/env python3
"""Run ERVcaller v1.4 on one DNA BAM and retain unfiltered insertion evidence.

No tumour/normal subtraction or RNA-locus matching is performed. Positions
are preserved verbatim as caller_pos: upstream v1.4 subtracts one in its VCF
writer, so downstream coordinate conversion requires explicit validation.
"""
import argparse
import csv
import hashlib
import json
from pathlib import Path
import re
import shutil
import subprocess
import sys


FIELDS = [
    "patient_id", "sample_type", "sample_id", "chrom", "caller_pos", "ref", "alt",
    "filter", "te_name", "te_start", "te_end", "te_length", "orientation",
    "breakpoint_status", "chimeric_reads", "split_reads", "genotyped",
    "insertion_read_ratio", "GT", "GQ", "GL", "DPN", "DPI",
]
SAFE_PATH = re.compile(r"[A-Za-z0-9_./-]+")


def read_vcf(path, patient, role, expected_sample=None):
    """Validate a single-sample ERVcaller VCF; retain all rows and raw positions."""
    rows = []
    sample = None
    fileformat = False
    with Path(path).open(encoding="utf-8") as handle:
        for number, line in enumerate(handle, 1):
            if line.startswith("##fileformat=VCFv4."):
                fileformat = True
            elif line.startswith("##"):
                continue
            elif line.startswith("#CHROM\t"):
                if sample is not None:
                    raise ValueError("Repeated VCF column header")
                columns = line.rstrip("\r\n").split("\t")
                if columns[:9] != ["#CHROM", "POS", "ID", "REF", "ALT", "QUAL", "FILTER", "INFO", "FORMAT"] or len(columns) != 10:
                    raise ValueError("Expected a single-sample ERVcaller VCF")
                sample = columns[9]
                if expected_sample and sample != expected_sample:
                    raise ValueError(f"VCF sample mismatch: {sample!r} != {expected_sample!r}")
            elif not line.strip():
                continue
            else:
                columns = line.rstrip("\r\n").split("\t")
                if sample is None or len(columns) != 10:
                    raise ValueError(f"Malformed VCF row {number}")
                if not columns[1].isdigit():
                    raise ValueError(f"Invalid caller position on row {number}")
                info = dict(item.split("=", 1) for item in columns[7].split(";") if "=" in item)
                te = info.get("INFOR", "").split(",")
                if len(te) != 6 or not {"CR", "SR", "GTF", "GR"} <= info.keys():
                    raise ValueError(f"Missing ERVcaller evidence fields on row {number}")
                for field in ("CR", "SR"):
                    if not info[field].isdigit():
                        raise ValueError(f"Invalid {field} on row {number}")
                keys = columns[8].split(":")
                values = columns[9].split(":")
                if len(keys) != len(values) or len(set(keys)) != len(keys):
                    raise ValueError(f"Malformed FORMAT data on row {number}")
                fmt = dict(zip(keys, values))
                if not {"GT", "GQ", "GL", "DPN", "DPI"} <= fmt.keys():
                    raise ValueError(f"Missing genotype fields on row {number}")
                rows.append(dict(zip(FIELDS, [
                    patient, role, sample, columns[0], columns[1], columns[3],
                    columns[4], columns[6], *te, info["CR"], info["SR"],
                    info["GTF"], info["GR"], *(fmt[k] for k in ("GT", "GQ", "GL", "DPN", "DPI")),
                ])))
    if not fileformat or sample is None:
        raise ValueError(
            "Missing VCF header. ERVcaller can produce an empty file on zero calls "
            "or an internal error; inspect the log. No successful zero-call result is assumed."
        )
    return rows


def check_reference_bam(header, fai):
    lengths = {}
    for line in Path(fai).read_text().splitlines():
        fields = line.split("\t")
        if len(fields) < 2 or fields[0] in lengths:
            raise ValueError("Malformed or duplicate FASTA index record")
        lengths[fields[0]] = int(fields[1])
    sq = []
    coordinate_sorted = False
    for line in header.splitlines():
        tags = dict(part.split(":", 1) for part in line.split("\t")[1:] if ":" in part)
        if line.startswith("@HD\t"):
            coordinate_sorted = tags.get("SO") == "coordinate"
        if line.startswith("@SQ\t"):
            name = tags.get("SN")
            if name not in lengths or lengths[name] != int(tags.get("LN", -1)):
                raise ValueError(f"BAM/reference contig mismatch: {name}")
            sq.append(name)
    if not coordinate_sorted or not sq:
        raise ValueError("A coordinate-sorted BAM with an @SQ dictionary is required")
    return sq[0]


def sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def command(args):
    return subprocess.run(args, check=True, text=True, capture_output=True)


def run(args):
    if not re.fullmatch(r"[A-Za-z0-9][A-Za-z0-9_.-]*", args.patient):
        raise ValueError("Unsafe patient identifier")
    if min(args.read_length, args.min_reads, args.min_split, args.threads) < 1:
        raise ValueError("Lengths, counts and threads must be positive")
    if args.min_split > args.read_length:
        raise ValueError("Minimum split length exceeds read length")
    sample = f"{args.patient}.{args.role}"
    # Keep aliases in the Nextflow task path. Upstream interpolates paths into
    # shell commands internally, so only shell-safe task paths can be accepted.
    refs = Path(args.references).absolute()
    caller = Path(args.caller).absolute()
    out = Path(args.outdir).absolute()
    work = Path("ervcaller_work").absolute()
    for path in (refs, caller, out, work):
        if not SAFE_PATH.fullmatch(str(path)):
            raise ValueError(f"ERVcaller requires Linux paths without whitespace/shell characters: {path}")
    for tool in ("perl", "bwa", "samtools", "Rscript", "extractSoftclipped"):
        if shutil.which(tool) is None:
            raise ValueError(f"Missing dependency in the ERVcaller container: {tool}")
    if not caller.is_file() or not (caller.parent / "Scripts").is_dir():
        raise ValueError("The ERVcaller installation must include its adjacent Scripts directory")
    for ref in ("genome.fa", "te.fa"):
        for suffix in ("", ".amb", ".ann", ".bwt", ".pac", ".sa", ".fai"):
            path = refs / (ref + suffix)
            if not path.is_file() or path.stat().st_size == 0:
                raise ValueError(f"Missing or empty prepared reference component: {path}")
    bam, bai = Path(args.bam).resolve(), Path(args.bai).resolve()
    if not bam.is_file() or not bai.is_file():
        raise ValueError("BAM and BAI are both required")
    command(["samtools", "quickcheck", "-v", str(bam)])
    header = command(["samtools", "view", "-H", str(bam)]).stdout
    first_contig = check_reference_bam(header, refs / "genome.fa.fai")
    work.mkdir(exist_ok=False)
    (work / f"{sample}.bam").symlink_to(bam)
    (work / f"{sample}.bam.bai").symlink_to(bai)
    # A region query forces index access (idxstats can fall back to scanning).
    command(["samtools", "view", "-c", str(work / f"{sample}.bam"), f"{first_contig}:1-1"])
    out.mkdir(parents=True, exist_ok=True)
    argv = [
        "perl", str(caller), "-i", sample, "-f", ".bam",
        "-H", str(refs / "genome.fa"), "-T", str(refs / "te.fa"),
        "-I", str(work) + "/", "-O", str(work) + "/",
        "-d", "WGS", "-s", "paired-end", "-r", str(args.read_length),
        "-n", str(args.min_reads), "-S", str(args.min_split), "-t", str(args.threads),
    ]
    if args.bwa_mem:
        argv.append("-B")
    if args.genotype:
        argv.append("-G")
    log_path = out / f"{sample}.log"
    with log_path.open("w") as log:
        subprocess.run(argv, stdout=log, stderr=subprocess.STDOUT, check=True, text=True)
    log_text = log_path.read_text(errors="replace")
    if re.search(r"command not found|Can't exec|Execution halted|No such file or directory", log_text, re.I):
        raise ValueError(f"ERVcaller reported an internal failure; inspect {log_path}")
    raw = work / f"{sample}.vcf"
    rows = read_vcf(raw, args.patient, args.role, expected_sample=sample)
    shutil.copyfile(raw, out / raw.name)
    with (out / f"{sample}.evidence.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, FIELDS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    status = {
        "patient_id": args.patient, "sample_type": args.role, "sample_id": sample,
        "insertion_records": len(rows), "genotyping_requested": args.genotype,
        "caller_sha256": sha256(caller), "genome_fai_sha256": sha256(refs / "genome.fa.fai"),
        "te_fasta_sha256": sha256(refs / "te.fa"),
        "command": argv, "samtools_version": command(["samtools", "--version"]).stdout.splitlines()[0],
        "interpretation": "Unfiltered non-reference insertion evidence; no somatic or peptide classification.",
        "coordinates": "caller_pos retained verbatim; validate v1.4 coordinate convention before genomic joins.",
    }
    manifest = refs / "references.sha256"
    if manifest.exists():
        status["reference_checksums"] = manifest.read_text()
    (out / f"{sample}.status.json").write_text(json.dumps(status, indent=2) + "\n")


def parser():
    p = argparse.ArgumentParser(description=__doc__)
    for name in ("patient", "bam", "bai", "references", "caller"):
        p.add_argument("--" + name, required=True)
    p.add_argument("--role", choices=("tumor", "normal"), required=True)
    p.add_argument("--read-length", type=int, required=True)
    p.add_argument("--min-reads", type=int, default=3)
    p.add_argument("--min-split", type=int, default=20)
    p.add_argument("--threads", type=int, default=2)
    p.add_argument("--bwa-mem", action="store_true")
    p.add_argument("--genotype", action="store_true")
    p.add_argument("--outdir", default="results")
    return p


if __name__ == "__main__":
    try:
        run(parser().parse_args())
    except (ValueError, OSError, subprocess.CalledProcessError) as exc:
        print(f"ERVcaller wrapper failed: {exc}", file=sys.stderr)
        if isinstance(exc, subprocess.CalledProcessError) and exc.stderr:
            print(exc.stderr, file=sys.stderr)
        sys.exit(1)
