#!/usr/bin/env python3
"""
Collect a patient's result tables into one .xlsx for whoever orders the
peptides.

The pipeline's own outputs are TSVs and a LaTeX report: good for a
pipeline, wrong for the person deciding what to synthesise, who wants to
sort, filter, hide columns and write in a comment. This writes one sheet
per table, with the header frozen and filterable, rejected rows in red and
cautioned rows in amber, and a legend sheet that spells out every flag so
the colours do not have to be explained again.

The file is written with the standard library only - no openpyxl, no
pandas - because this runs inside the pVACtools container and adding a
dependency there to format a spreadsheet is a bad trade. The format is
plain SpreadsheetML with inline strings: cells that parse as numbers are
written as numbers so sorting works, everything else as text.

Sheet names come from --sheet NAME=path, in the order given. Names are
truncated to Excel's 31 characters and its forbidden characters are
replaced, so a long table name cannot produce a file Excel refuses.

Usage:
    write_workbook.py -o P1.xlsx \\
        --sheet 'Class I=P1.MHC_I.prioritized.tsv' \\
        --sheet 'Long peptides=P1.Combined.long_peptides.tsv' \\
        --sheet 'HLA LOH=P1.hla_loh.tsv'
"""
import argparse
import csv
import re
import sys
import zipfile
from pathlib import Path
from xml.sax.saxutils import escape

# Cell styles, in the order they appear in styles.xml below.
STYLE_NORMAL, STYLE_HEADER, STYLE_REJECTED, STYLE_CAUTION, STYLE_WRAP = 0, 1, 2, 3, 4

FLAG_LEGEND = [
    ("Selection: Top ranking", "No criterion failed. These are the candidates to make."),
    ("Selection: Caution", "Usable, but something is weaker than its neighbours. Read the Flags column."),
    ("Selection: Rejected", "Should not go into a vaccine. The Reject Reason column says why."),
    ("HLA_LOST", "The restricting HLA allele is deleted in the tumour, so this epitope cannot be presented at all."),
    ("LOW_BINDING", "Predicted binding is weak on both percentile and IC50."),
    ("NO_EXPRESSION", "The gene is not expressed in this tumour, so the protein is not there to be processed."),
    ("NO_RNA_SUPPORT", "RNA reads cover the variant site but none carry the variant."),
    ("SELF_LIKE", "The mutant peptide binds no better than its wild-type counterpart, which normal cells also present."),
    ("LOW_EXPRESSION", "The gene is expressed only weakly."),
    ("NORMAL_TISSUE", "The gene is broadly expressed in healthy tissue, so an on-target off-tumour response is a risk."),
    ("SUBCLONAL", "Present in a minority of tumour cells (CCF below the threshold), so it cannot cover the whole tumour."),
    ("LOW_DNA_VAF", "Low variant allele fraction in tumour DNA."),
    ("LOW_DEPTH", "Low DNA coverage at the variant, so every other number for it is noisy."),
    ("NO_EXPRESSION_DATA", "No expression evidence either way. Not a mark against the candidate."),
    ("NO_RNA_DATA", "No RNA coverage at the variant. Not a mark against the candidate."),
    ("HLA_NOT_EVALUABLE", "HLA-LOH could not be called for this allele, usually a homozygous locus or too few heterozygous sites."),
    ("Agretopicity", "Mutant IC50 divided by wild-type IC50. Below 1 the mutation is what creates the epitope; near 1 the wild type is presented too."),
    ("CCF", "Cancer cell fraction: the share of tumour cells carrying the variant, corrected for purity and local copy number."),
    ("Warnings (long peptides)", "Synthesis notes only: cysteine and methionine counts, hydrophobicity (GRAVY), residue runs. Nothing was dropped for these."),
]

ILLEGAL = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f]")
NUMERIC = re.compile(r"^-?(\d+\.?\d*|\.\d+)([eE][-+]?\d+)?$")


def column_name(index):
    """1 -> A, 27 -> AA."""
    name = ""
    while index > 0:
        index, remainder = divmod(index - 1, 26)
        name = chr(ord("A") + remainder) + name
    return name


def sheet_name(raw, taken):
    """An Excel-legal, unique sheet name."""
    name = re.sub(r"[\[\]:*?/\\]", "-", str(raw)).strip("'") or "Sheet"
    name = name[:31]
    if name.lower() in taken:
        for suffix in range(2, 100):
            candidate = f"{name[:31 - len(str(suffix)) - 1]}_{suffix}"
            if candidate.lower() not in taken:
                name = candidate
                break
    taken.add(name.lower())
    return name


def row_style(row, header):
    """Rejected rows red, cautioned rows amber, from the Selection column."""
    if "Selection" not in header:
        return STYLE_NORMAL
    value = (row[header.index("Selection")] if header.index("Selection") < len(row) else "").strip().lower()
    if value.startswith("reject"):
        return STYLE_REJECTED
    if value.startswith("caution"):
        return STYLE_CAUTION
    return STYLE_NORMAL


def cell_xml(reference, value, style):
    """One <c>. Numbers stay numbers so the column sorts as one."""
    text = ILLEGAL.sub("", "" if value is None else str(value))
    attributes = f' s="{style}"' if style else ""
    if text == "":
        return f'<c r="{reference}"{attributes}/>'
    if NUMERIC.match(text):
        # Keep long digit strings (identifiers) as text: Excel would
        # otherwise round them or show them in scientific notation.
        if not (len(text) > 11 and text.isdigit()):
            return f'<c r="{reference}"{attributes}><v>{text}</v></c>'
    return f'<c r="{reference}"{attributes} t="inlineStr"><is><t xml:space="preserve">{escape(text)}</t></is></c>'


def sheet_xml(rows):
    """Rows of (values, style) -> a worksheet part."""
    header = [str(v) for v in rows[0][0]] if rows else []
    width = max((len(values) for values, _ in rows), default=1)
    body = []
    for number, (values, style) in enumerate(rows, start=1):
        cells = "".join(
            cell_xml(f"{column_name(index + 1)}{number}", value, style)
            for index, value in enumerate(values)
        )
        body.append(f'<row r="{number}">{cells}</row>')
    dimension = f"A1:{column_name(max(width, 1))}{max(len(rows), 1)}"
    # Freeze the header and turn on the filter dropdowns, which is the
    # whole point of handing this out as a spreadsheet.
    panes = (
        '<sheetViews><sheetView workbookViewId="0">'
        '<pane ySplit="1" topLeftCell="A2" activePane="bottomLeft" state="frozen"/>'
        "</sheetView></sheetViews>"
    )
    columns = "".join(
        f'<col min="{i + 1}" max="{i + 1}" width="{min(60, max(10, len(str(name)) + 4))}" customWidth="1"/>'
        for i, name in enumerate(header)
    )
    autofilter = f'<autoFilter ref="{dimension}"/>' if len(rows) > 1 else ""
    return (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<worksheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
        f'<dimension ref="{dimension}"/>{panes}'
        f'{f"<cols>{columns}</cols>" if columns else ""}'
        f"<sheetData>{''.join(body)}</sheetData>{autofilter}</worksheet>"
    )


STYLES_XML = (
    '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
    '<styleSheet xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main">'
    '<fonts count="4">'
    '<font><sz val="11"/><name val="Calibri"/></font>'
    '<font><b/><sz val="11"/><color rgb="FFFFFFFF"/><name val="Calibri"/></font>'
    '<font><sz val="11"/><color rgb="FF9C0006"/><name val="Calibri"/></font>'
    '<font><sz val="11"/><color rgb="FF9C6500"/><name val="Calibri"/></font>'
    "</fonts>"
    '<fills count="5">'
    '<fill><patternFill patternType="none"/></fill>'
    '<fill><patternFill patternType="gray125"/></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FF31557A"/><bgColor indexed="64"/></patternFill></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FFFFC7CE"/><bgColor indexed="64"/></patternFill></fill>'
    '<fill><patternFill patternType="solid"><fgColor rgb="FFFFEB9C"/><bgColor indexed="64"/></patternFill></fill>'
    "</fills>"
    '<borders count="1"><border/></borders>'
    '<cellStyleXfs count="1"><xf numFmtId="0" fontId="0" fillId="0" borderId="0"/></cellStyleXfs>'
    '<cellXfs count="5">'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0"/>'
    '<xf numFmtId="0" fontId="1" fillId="2" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="2" fillId="3" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="3" fillId="4" borderId="0" xfId="0" applyFont="1" applyFill="1"/>'
    '<xf numFmtId="0" fontId="0" fillId="0" borderId="0" xfId="0" applyAlignment="1">'
    '<alignment vertical="top" wrapText="1"/></xf>'
    "</cellXfs>"
    '<cellStyles count="1"><cellStyle name="Normal" xfId="0" builtinId="0"/></cellStyles>'
    "</styleSheet>"
)


def write_xlsx(path, sheets):
    """sheets: [(name, [(values, style), ...])] -> an .xlsx at path."""
    names = [name for name, _ in sheets]
    sheet_entries = "".join(
        f'<sheet name="{escape(name)}" sheetId="{i + 1}" r:id="rId{i + 1}"/>'
        for i, name in enumerate(names)
    )
    workbook = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<workbook xmlns="http://schemas.openxmlformats.org/spreadsheetml/2006/main"'
        ' xmlns:r="http://schemas.openxmlformats.org/officeDocument/2006/relationships">'
        f"<sheets>{sheet_entries}</sheets></workbook>"
    )
    relationships = "".join(
        f'<Relationship Id="rId{i + 1}" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/worksheet" Target="worksheets/sheet{i + 1}.xml"/>'
        for i in range(len(names))
    )
    workbook_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        f'{relationships}<Relationship Id="rIdStyles" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/styles" Target="styles.xml"/>'
        "</Relationships>"
    )
    overrides = "".join(
        f'<Override PartName="/xl/worksheets/sheet{i + 1}.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.worksheet+xml"/>'
        for i in range(len(names))
    )
    content_types = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Types xmlns="http://schemas.openxmlformats.org/package/2006/content-types">'
        '<Default Extension="rels" ContentType="application/vnd.openxmlformats-package.relationships+xml"/>'
        '<Default Extension="xml" ContentType="application/xml"/>'
        '<Override PartName="/xl/workbook.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet.main+xml"/>'
        '<Override PartName="/xl/styles.xml" ContentType="application/vnd.openxmlformats-officedocument.spreadsheetml.styles+xml"/>'
        f"{overrides}</Types>"
    )
    root_rels = (
        '<?xml version="1.0" encoding="UTF-8" standalone="yes"?>'
        '<Relationships xmlns="http://schemas.openxmlformats.org/package/2006/relationships">'
        '<Relationship Id="rId1" Type="http://schemas.openxmlformats.org/officeDocument/2006/relationships/officeDocument" Target="xl/workbook.xml"/>'
        "</Relationships>"
    )

    with zipfile.ZipFile(path, "w", zipfile.ZIP_DEFLATED) as archive:
        archive.writestr("[Content_Types].xml", content_types)
        archive.writestr("_rels/.rels", root_rels)
        archive.writestr("xl/workbook.xml", workbook)
        archive.writestr("xl/_rels/workbook.xml.rels", workbook_rels)
        archive.writestr("xl/styles.xml", STYLES_XML)
        for index, (_, rows) in enumerate(sheets, start=1):
            archive.writestr(f"xl/worksheets/sheet{index}.xml", sheet_xml(rows))


def read_table(path, max_rows=None):
    """A TSV or CSV -> [(values, style)], header styled and rows coloured."""
    text = Path(path).read_text()
    if not text.strip():
        return []
    delimiter = "\t" if "\t" in text.splitlines()[0] else ","
    rows = list(csv.reader(text.splitlines(), delimiter=delimiter))
    if not rows:
        return []
    header = [str(value) for value in rows[0]]
    body = rows[1:max_rows + 1] if max_rows else rows[1:]
    return [(header, STYLE_HEADER)] + [(row, row_style(row, header)) for row in body]


def legend_rows():
    return [(["Column or flag", "What it means"], STYLE_HEADER)] + [
        ([name, text], STYLE_WRAP) for name, text in FLAG_LEGEND
    ]


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter)
    parser.add_argument("--sheet", action="append", default=[], metavar="NAME=PATH",
                        help="a table to add, in order; repeat for each sheet")
    parser.add_argument("--max-rows", type=int, help="keep only the first N data rows of each sheet")
    parser.add_argument("--no-legend", action="store_true", help="leave out the legend sheet")
    parser.add_argument("-o", "--output", type=Path, required=True)
    args = parser.parse_args(argv)

    sheets, taken, missing = [], set(), []
    for item in args.sheet:
        name, separator, path = item.partition("=")
        if not separator:
            parser.error(f"--sheet needs NAME=PATH, got: {item}")
        if not Path(path).exists():
            missing.append(path)
            continue
        rows = read_table(path, args.max_rows)
        if not rows:
            missing.append(path)
            continue
        sheets.append((sheet_name(name, taken), rows))

    if not args.no_legend:
        sheets.append((sheet_name("Legend", taken), legend_rows()))
    if not sheets:
        parser.error("no readable tables to write")

    write_xlsx(args.output, sheets)
    print(f"{args.output.name}: {len(sheets)} sheets ({', '.join(name for name, _ in sheets)})", file=sys.stderr)
    for path in missing:
        print(f"note: skipped {path} (missing or empty)", file=sys.stderr)


if __name__ == "__main__":
    main()
