from datetime import datetime
from typing import Optional, List, Dict, Any
from zipfile import ZipFile
import xml.etree.ElementTree as ET

import pandas as pd


class ExcelsheetParser:
    KEY_COLUMN = "NCATSDPI_Variable_Name"
    VALUE_COLUMN = "Submitter_Value"
    MAPPED_VALUE_COLUMN = "Submitter_Variable_Name"

    file_path: str
    sheet_dfs: Dict[str, pd.DataFrame]
    _parameter_map_cache: Dict[str, Dict[str, str]]
    _sheet_xml_paths: Dict[str, str]

    def __init__(self, file_path: str):
        self.file_path = file_path
        self.sheet_dfs = self._read_all_sheets()
        self._parameter_map_cache = {}
        self._sheet_xml_paths = self._load_sheet_xml_paths()

    @classmethod
    def from_sheets(cls, file_name: str, sheets: Dict[str, pd.DataFrame]):
        """Construct a parser from sheet data supplied by another transport."""
        parser = cls.__new__(cls)
        parser.file_path = file_name
        parser.sheet_dfs = {}
        for name, frame in sheets.items():
            frame = frame.copy()
            frame.columns = [c.strip() if isinstance(c, str) else c for c in frame.columns]
            frame = frame.apply(lambda col: col.map(
                lambda value: value.replace('\xa0', ' ').strip()
                if isinstance(value, str) else value
            ))
            parser.sheet_dfs[name] = frame
        parser._parameter_map_cache = {}
        parser._sheet_xml_paths = {}
        return parser

    def _read_all_sheets(self) -> Dict[str, pd.DataFrame]:
        xls = pd.ExcelFile(self.file_path)
        sheets = {}
        for sheet in xls.sheet_names:
            df = xls.parse(sheet, keep_default_na=False)
            df.columns = [c.strip() if isinstance(c, str) else c for c in df.columns]
            df = df.apply(lambda col: col.map(
                lambda x: x.replace('\xa0', ' ').strip() if isinstance(x, str) else x
            ))
            sheets[sheet] = df
        return sheets

    def _get_one_value(self, sheet_df: pd.DataFrame, data_key: str) -> Any:
        matches = sheet_df[sheet_df[self.KEY_COLUMN] == data_key][self.VALUE_COLUMN]
        if len(matches) == 0:
            raise LookupError(f"Could not find row for '{data_key}'", sheet_df[self.KEY_COLUMN].values)
        if len(matches) != 1:
            raise LookupError(f"Found multiple rows for '{data_key}' There must be one and only one", matches.values)
        value = matches.values[0]
        if pd.isna(value) or value == '':
            return None
        return value

    def get_one_string(self, sheet_name: str, data_key: str) -> Optional[str]:
        sheet_df = self.sheet_dfs[sheet_name]
        value = self._get_one_value(sheet_df, data_key)
        if value is None:
            return None
        return str(value).strip()

    def get_one_string_list(self, sheet_name: str, data_key: str, delimiter: str = '|') -> List[str]:
        sheet_df = self.sheet_dfs[sheet_name]
        value_list = self._get_one_value(sheet_df, data_key)
        if value_list is None:
            return []
        return [val.strip() for val in str(value_list).split(delimiter)]

    def get_one_date(self, sheet_name: str, data_key: str) -> datetime:
        sheet_df = self.sheet_dfs[sheet_name]
        date_string = self._get_one_value(sheet_df, data_key)
        if isinstance(date_string, datetime):
            return date_string
        if date_string is None:
            raw_value = self._get_raw_submitter_value(sheet_name, data_key)
            if raw_value not in (None, ""):
                date_string = raw_value
        if isinstance(date_string, (int, float)) and not isinstance(date_string, bool):
            date_string = str(int(date_string))
        if isinstance(date_string, str):
            stripped = date_string.strip()
            if len(stripped) == 8 and stripped.isdigit():
                return datetime.strptime(stripped, "%Y%m%d")
            try:
                return datetime.fromisoformat(stripped)
            except ValueError:
                pass
        return datetime.strptime(str(date_string), "%Y%m%d")

    def get_parameter_map(self, sheet_name: str) -> Dict[str, str]:
        if sheet_name in self._parameter_map_cache:
            return self._parameter_map_cache[sheet_name]

        sheet_df = self.sheet_dfs[sheet_name]
        data_frame = sheet_df.dropna(subset=[self.KEY_COLUMN, self.MAPPED_VALUE_COLUMN])
        series = pd.Series(data_frame[self.MAPPED_VALUE_COLUMN].values, index=data_frame[self.KEY_COLUMN])
        result = series.to_dict()

        self._parameter_map_cache[sheet_name] = result
        return result

    def safe_get_string(self, sheet_name: str, data_key: str) -> Optional[str]:
        """Like get_one_string, but returns None instead of raising if the key is missing."""
        try:
            return self.get_one_string(sheet_name, data_key)
        except (LookupError, KeyError):
            return None

    def get_mapped_value(self, sheet_name: str, key: str) -> str:
        param_map = self.get_parameter_map(sheet_name=sheet_name)
        return param_map[key]

    def _load_sheet_xml_paths(self) -> Dict[str, str]:
        ns = {
            "main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main",
            "rel": "http://schemas.openxmlformats.org/officeDocument/2006/relationships",
            "pkg": "http://schemas.openxmlformats.org/package/2006/relationships",
        }
        with ZipFile(self.file_path) as zf:
            workbook_root = ET.fromstring(zf.read("xl/workbook.xml"))
            rels_root = ET.fromstring(zf.read("xl/_rels/workbook.xml.rels"))
        rel_map = {
            rel.attrib["Id"]: rel.attrib["Target"]
            for rel in rels_root.findall("pkg:Relationship", ns)
        }
        sheet_paths = {}
        for sheet in workbook_root.findall("main:sheets/main:sheet", ns):
            name = sheet.attrib["name"]
            rel_id = sheet.attrib[f"{{{ns['rel']}}}id"]
            target = rel_map[rel_id]
            sheet_paths[name] = target if target.startswith("xl/") else f"xl/{target}"
        return sheet_paths

    def _get_raw_submitter_value(self, sheet_name: str, data_key: str) -> Optional[str]:
        sheet_df = self.sheet_dfs[sheet_name]
        matches = sheet_df.index[sheet_df[self.KEY_COLUMN] == data_key].tolist()
        if len(matches) != 1:
            return None

        sheet_xml_path = self._sheet_xml_paths.get(sheet_name)
        if not sheet_xml_path:
            return None

        excel_row_number = matches[0] + 2  # header row + 1-based indexing
        value_column_number = sheet_df.columns.get_loc(self.VALUE_COLUMN) + 1
        cell_ref = f"{self._column_number_to_letters(value_column_number)}{excel_row_number}"

        ns = {"main": "http://schemas.openxmlformats.org/spreadsheetml/2006/main"}
        with ZipFile(self.file_path) as zf:
            sheet_root = ET.fromstring(zf.read(sheet_xml_path))

        for cell in sheet_root.findall(f".//main:c[@r='{cell_ref}']", ns):
            raw_value = cell.findtext("main:v", default=None, namespaces=ns)
            if raw_value is not None:
                return raw_value.strip()
        return None

    @staticmethod
    def _column_number_to_letters(column_number: int) -> str:
        letters = []
        n = column_number
        while n > 0:
            n, remainder = divmod(n - 1, 26)
            letters.append(chr(65 + remainder))
        return "".join(reversed(letters))
