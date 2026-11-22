"""write_workbook.py: result tables -> one .xlsx, written with the stdlib."""
import csv
import importlib.util
import tempfile
import unittest
import xml.etree.ElementTree as ET
import zipfile
from pathlib import Path

BIN = Path(__file__).resolve().parent.parent / "bin" / "write_workbook.py"
spec = importlib.util.spec_from_file_location("write_workbook", BIN)
workbook = importlib.util.module_from_spec(spec)
spec.loader.exec_module(workbook)

NS = {"m": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}


def write_tsv(rows, directory=None, name="table.tsv"):
    path = Path(directory or tempfile.mkdtemp()) / name
    with path.open("w", newline="") as handle:
        csv.writer(handle, delimiter="\t", lineterminator="\n").writerows(rows)
    return path


def build(sheets, **kwargs):
    """sheets: {name: rows} -> the opened zipfile of the written workbook."""
    directory = Path(tempfile.mkdtemp())
    argv = []
    for index, (name, rows) in enumerate(sheets.items()):
        path = write_tsv(rows, directory, f"table{index}.tsv")
        argv += ["--sheet", f"{name}={path}"]
    out = directory / "book.xlsx"
    argv += ["-o", str(out)]
    for key, value in kwargs.items():
        flag = f"--{key.replace('_', '-')}"
        argv += [flag] if value is True else [flag, str(value)]
    workbook.main(argv)
    return zipfile.ZipFile(out), out


def cells(archive, sheet_index=1):
    """-> {'A1': (text, style)} for one worksheet."""
    root = ET.fromstring(archive.read(f"xl/worksheets/sheet{sheet_index}.xml"))
    found = {}
    for cell in root.iter(f"{{{NS['m']}}}c"):
        reference = cell.get("r")
        value = cell.find("m:v", NS)
        inline = cell.find("m:is/m:t", NS)
        text = value.text if value is not None else (inline.text if inline is not None else "")
        found[reference] = (text, cell.get("s") or "0")
    return found


SIMPLE = [["Gene", "CCF", "Selection"],
          ["GENE1", "0.9", "Top ranking"],
          ["GENE2", "0.2", "Caution"],
          ["GENE3", "0.5", "Rejected"]]


class TestColumnNames(unittest.TestCase):
    def test_single_letters(self):
        self.assertEqual(workbook.column_name(1), "A")
        self.assertEqual(workbook.column_name(26), "Z")

    def test_two_letters(self):
        self.assertEqual(workbook.column_name(27), "AA")
        self.assertEqual(workbook.column_name(52), "AZ")
        self.assertEqual(workbook.column_name(53), "BA")


class TestSheetNames(unittest.TestCase):
    def test_illegal_characters_are_replaced(self):
        self.assertEqual(workbook.sheet_name("Class I/II: [draft]", set()), "Class I-II- -draft-")

    def test_long_names_are_truncated(self):
        name = workbook.sheet_name("x" * 60, set())
        self.assertEqual(len(name), 31)

    def test_duplicates_are_made_unique(self):
        taken = set()
        first = workbook.sheet_name("Candidates", taken)
        second = workbook.sheet_name("Candidates", taken)
        self.assertEqual(first, "Candidates")
        self.assertNotEqual(first, second)

    def test_empty_name(self):
        self.assertEqual(workbook.sheet_name("", set()), "Sheet")


class TestWorkbookStructure(unittest.TestCase):
    def setUp(self):
        self.archive, self.path = build({"Candidates": SIMPLE})

    def test_required_parts_exist(self):
        names = self.archive.namelist()
        for part in ("[Content_Types].xml", "_rels/.rels", "xl/workbook.xml",
                     "xl/_rels/workbook.xml.rels", "xl/styles.xml", "xl/worksheets/sheet1.xml"):
            self.assertIn(part, names)

    def test_every_part_is_well_formed_xml(self):
        for name in self.archive.namelist():
            ET.fromstring(self.archive.read(name))

    def test_the_zip_is_not_corrupt(self):
        self.assertIsNone(self.archive.testzip())

    def test_sheet_is_named_and_registered(self):
        root = ET.fromstring(self.archive.read("xl/workbook.xml"))
        names = [s.get("name") for s in root.iter(f"{{{NS['m']}}}sheet")]
        self.assertEqual(names, ["Candidates", "Legend"])

    def test_a_legend_sheet_is_added(self):
        self.assertIn("xl/worksheets/sheet2.xml", self.archive.namelist())
        text = self.archive.read("xl/worksheets/sheet2.xml").decode()
        self.assertIn("HLA_LOST", text)

    def test_legend_can_be_turned_off(self):
        archive, _ = build({"Candidates": SIMPLE}, no_legend=True)
        root = ET.fromstring(archive.read("xl/workbook.xml"))
        self.assertEqual([s.get("name") for s in root.iter(f"{{{NS['m']}}}sheet")], ["Candidates"])

    def test_header_is_frozen_and_filterable(self):
        text = self.archive.read("xl/worksheets/sheet1.xml").decode()
        self.assertIn('state="frozen"', text)
        self.assertIn("<autoFilter", text)


class TestCells(unittest.TestCase):
    def setUp(self):
        self.archive, _ = build({"Candidates": SIMPLE})
        self.cells = cells(self.archive)

    def test_header_text_and_style(self):
        self.assertEqual(self.cells["A1"], ("Gene", str(workbook.STYLE_HEADER)))
        self.assertEqual(self.cells["C1"][0], "Selection")

    def test_text_cell(self):
        self.assertEqual(self.cells["A2"][0], "GENE1")

    def test_numbers_are_written_as_numbers(self):
        root = ET.fromstring(self.archive.read("xl/worksheets/sheet1.xml"))
        cell = [c for c in root.iter(f"{{{NS['m']}}}c") if c.get("r") == "B2"][0]
        self.assertIsNone(cell.get("t"))
        self.assertEqual(cell.find("m:v", NS).text, "0.9")

    def test_rejected_rows_are_styled(self):
        self.assertEqual(self.cells["A4"][1], str(workbook.STYLE_REJECTED))

    def test_cautioned_rows_are_styled(self):
        self.assertEqual(self.cells["A3"][1], str(workbook.STYLE_CAUTION))

    def test_top_ranking_rows_are_plain(self):
        self.assertEqual(self.cells["A2"][1], str(workbook.STYLE_NORMAL))

    def test_long_identifiers_stay_text(self):
        archive, _ = build({"S": [["ID"], ["123456789012345"]]})
        root = ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        cell = [c for c in root.iter(f"{{{NS['m']}}}c") if c.get("r") == "A2"][0]
        self.assertEqual(cell.get("t"), "inlineStr")

    def test_xml_special_characters_are_escaped(self):
        archive, _ = build({"S": [["Note"], ["a < b & c > d"]]})
        self.assertEqual(cells(archive)["A2"][0], "a < b & c > d")

    def test_empty_cells_are_written_without_a_value(self):
        archive, _ = build({"S": [["A", "B"], ["x", ""]]})
        text = archive.read("xl/worksheets/sheet1.xml").decode()
        self.assertIn('<c r="B2"/>', text)

    def test_peptide_sequences_are_not_mangled(self):
        archive, _ = build({"S": [["Long Peptide"], ["SLYNTVATLGSAKQEM"]]})
        self.assertEqual(cells(archive)["A2"][0], "SLYNTVATLGSAKQEM")


class TestInputHandling(unittest.TestCase):
    def test_csv_is_read_too(self):
        directory = Path(tempfile.mkdtemp())
        path = directory / "t.csv"
        path.write_text("Gene,CCF\nGENE1,0.9\n")
        out = directory / "book.xlsx"
        workbook.main(["--sheet", f"Candidates={path}", "-o", str(out)])
        self.assertEqual(cells(zipfile.ZipFile(out))["A2"][0], "GENE1")

    def test_max_rows_truncates(self):
        archive, _ = build({"Candidates": SIMPLE}, max_rows=1)
        self.assertNotIn("A3", cells(archive))

    def test_a_missing_table_is_skipped_not_fatal(self):
        directory = Path(tempfile.mkdtemp())
        present = write_tsv(SIMPLE, directory)
        out = directory / "book.xlsx"
        workbook.main(["--sheet", f"Candidates={present}",
                       "--sheet", f"Gone={directory / 'nope.tsv'}", "-o", str(out)])
        root = ET.fromstring(zipfile.ZipFile(out).read("xl/workbook.xml"))
        self.assertEqual([s.get("name") for s in root.iter(f"{{{NS['m']}}}sheet")],
                         ["Candidates", "Legend"])

    def test_an_empty_table_is_skipped(self):
        directory = Path(tempfile.mkdtemp())
        empty = directory / "empty.tsv"
        empty.write_text("")
        out = directory / "book.xlsx"
        workbook.main(["--sheet", f"Empty={empty}", "-o", str(out)])
        root = ET.fromstring(zipfile.ZipFile(out).read("xl/workbook.xml"))
        self.assertEqual([s.get("name") for s in root.iter(f"{{{NS['m']}}}sheet")], ["Legend"])

    def test_sheet_without_a_name_is_an_error(self):
        directory = Path(tempfile.mkdtemp())
        path = write_tsv(SIMPLE, directory)
        with self.assertRaises(SystemExit):
            workbook.main(["--sheet", str(path), "-o", str(directory / "b.xlsx")])

    def test_nothing_readable_at_all_is_an_error(self):
        directory = Path(tempfile.mkdtemp())
        with self.assertRaises(SystemExit):
            workbook.main(["--sheet", f"Gone={directory / 'nope.tsv'}",
                           "--no-legend", "-o", str(directory / "b.xlsx")])

    def test_sheet_order_follows_the_arguments(self):
        archive, _ = build({"First": SIMPLE, "Second": SIMPLE})
        root = ET.fromstring(archive.read("xl/workbook.xml"))
        self.assertEqual([s.get("name") for s in root.iter(f"{{{NS['m']}}}sheet")],
                         ["First", "Second", "Legend"])

    def test_ragged_rows_do_not_break_the_sheet(self):
        archive, _ = build({"S": [["A", "B", "C"], ["x"], ["y", "z"]]})
        ET.fromstring(archive.read("xl/worksheets/sheet1.xml"))
        self.assertEqual(cells(archive)["A3"][0], "y")


if __name__ == "__main__":
    unittest.main()
