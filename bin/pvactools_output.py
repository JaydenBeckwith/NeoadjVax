#!/usr/bin/env python3
"""
Find and classify pVACtools aggregated reports, so a script can be pointed
straight at a pvacseq/ or pvacfuse/ output directory.

Both tools lay their results out the same way:

    pvacseq/MHC_Class_I/tumor.MHC_I.all_epitopes.aggregated.tsv
    pvacseq/MHC_Class_II/tumor.MHC_II.all_epitopes.aggregated.tsv
    pvacseq/combined/tumor.Combined.all_epitopes.aggregated.tsv

so one glob finds either, and the MHC class comes out of the file name.

The two are not interchangeable in what they carry, and that difference
decides what can be computed from them:

    somatic (pVACseq)  has IC50 WT / %ile WT and DNA VAF, so agretopicity
                       and cancer cell fraction are both computable.
    fusion (pVACfuse)  has neither. A fusion junction peptide has no
                       wild-type counterpart to compare against, so
                       agretopicity is undefined rather than missing; and
                       fusions are called from RNA, so there is no DNA VAF
                       to turn into a CCF.

That is why `kind()` exists: a script should say which of the two it is
looking at and decline the calculation that does not apply, rather than
filling a column with NA that reads like a failed run.
"""
import csv
from pathlib import Path

AGGREGATED_GLOB = "*.all_epitopes.aggregated.tsv"

# A somatic report is identified by the wild-type columns pVACseq writes
# and pVACfuse has no equivalent of.
SOMATIC_MARKERS = ("IC50 WT", "%ile WT", "DNA VAF")


def mhc_class(path):
    """'MHC_I', 'MHC_II', 'Combined', or the parent directory's name.

    pVACtools names the file <sample>.<class>.all_epitopes.aggregated.tsv,
    so the class is the second dot-separated field. A renamed file falls
    back to the directory it sits in, which carries the same information.
    """
    fields = Path(path).name.split(".")
    if len(fields) > 1 and fields[1]:
        candidate = fields[1]
        if candidate.lower().startswith("mhc") or candidate.lower() == "combined":
            return candidate
    parent = Path(path).parent.name
    return parent if parent else "unknown"


def find_reports(path):
    """A file or a pvacseq/pvacfuse directory -> [(mhc class, Path)].

    A file is taken as given. A directory is searched one level down
    (MHC_Class_I/, MHC_Class_II/, combined/) and then at the top level, so
    both a tool's own output tree and a flat folder of reports work.
    """
    path = Path(path)
    if path.is_file():
        return [(mhc_class(path), path)]
    if not path.is_dir():
        return []
    found = sorted(path.glob(f"*/{AGGREGATED_GLOB}")) or sorted(path.glob(AGGREGATED_GLOB))
    return [(mhc_class(p), p) for p in found]


def header(path):
    """The report's column names, or [] for an empty file."""
    with open(path, newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        try:
            return next(reader)
        except StopIteration:
            return []


def kind(columns):
    """'somatic', 'fusion' or 'unknown' from a report's columns."""
    if not columns:
        return "unknown"
    if any(marker in columns for marker in SOMATIC_MARKERS):
        return "somatic"
    # pVACfuse's own columns. Checked second so a somatic report that
    # happens to carry an Expr column is still called somatic.
    if any(marker in columns for marker in ("Read Support", "Expr")):
        return "fusion"
    return "unknown"


def describe(path):
    """-> (mhc class, kind) for one report."""
    return mhc_class(path), kind(header(path))
