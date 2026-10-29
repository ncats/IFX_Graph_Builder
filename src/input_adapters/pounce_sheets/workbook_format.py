"""Detect and normalize POUNCE workbook templates at the input boundary.

The original workbook remains untouched.  Both validation and graph building
consume the same in-memory sheet names, map headers, and field keys.
"""

from dataclasses import dataclass
import re
from typing import Literal

from src.input_adapters.excel_sheet_adapter import ExcelsheetParser


WorkbookRole = Literal["project", "experiment", "stats"]

_ROLE_SIGNATURES = {
    "project": {
        "anchor": "ProjectMeta",
        "maps": {"BioSampleMap", "BiospecimenMap"},
        "required_maps": {"BioSampleMap", "BiospecimenMap"},
    },
    "experiment": {
        "anchor": "ExperimentMeta",
        "maps": {"RunBioSampleMap", "GeneMap", "MetabMap", "ProteinMap"},
        "required_maps": {"RunBioSampleMap"},
    },
    "stats": {
        "anchor": "StatsResultsMeta",
        "maps": {"EffectSize_Map"},
        "required_maps": {"EffectSize_Map"},
    },
}
_V2_MAP_HEADER = ExcelsheetParser.MAPPED_VALUE_COLUMN
_V3_MAP_HEADER = "Submitter_Column_Name"


class UnsupportedPounceFormat(ValueError):
    """A workbook does not have one unambiguous supported template format."""


def detect_pounce_version(parser: ExcelsheetParser, role: WorkbookRole) -> int:
    """Return 2 or 3 from workbook contents, independent of its filename.

    Every present map sheet must agree.  A missing role-specific anchor or a
    mixture of map headers is an error rather than an arbitrary guess.
    """
    if role not in _ROLE_SIGNATURES:
        raise ValueError(f"Unknown POUNCE workbook role: {role}")
    signature = _ROLE_SIGNATURES[role]
    map_names = signature["maps"]
    anchor = signature["anchor"]
    required_maps = signature["required_maps"]
    actual = set(parser.sheet_dfs)
    if anchor not in actual:
        raise UnsupportedPounceFormat(
            f"{parser.file_path}: expected {role} workbook with sheet '{anchor}'"
        )
    if not actual & required_maps:
        raise UnsupportedPounceFormat(
            f"{parser.file_path}: missing {role} mapping sheet; expected one of {sorted(required_maps)}"
        )
    present = sorted(actual & map_names)
    versions = set()
    for name in present:
        columns = set(parser.sheet_dfs[name].columns)
        has_v2 = _V2_MAP_HEADER in columns
        has_v3 = _V3_MAP_HEADER in columns
        if has_v2 == has_v3:
            raise UnsupportedPounceFormat(
                f"{parser.file_path} / {name}: expected exactly one of "
                f"'{_V2_MAP_HEADER}' or '{_V3_MAP_HEADER}'"
            )
        versions.add(2 if has_v2 else 3)
    if len(versions) != 1:
        raise UnsupportedPounceFormat(
            f"{parser.file_path}: mapping sheets disagree about template version"
        )
    version = versions.pop()
    if role == "project" and "BiospecimenMap" in actual and version != 3:
        raise UnsupportedPounceFormat(
            f"{parser.file_path}: BiospecimenMap requires a v3 map header"
        )
    return version


@dataclass(frozen=True)
class PounceWorkbook:
    parser: ExcelsheetParser
    role: WorkbookRole
    version: int
    original_sheets: dict[str, str]

    def original_sheet(self, normalized_name: str) -> str:
        return self.original_sheets.get(normalized_name, normalized_name)


_V3_ORIGINAL_KEYS = {
    "ProjectMeta": {
        "lab_groups": "ncats_lab_groups",
    },
    "RunBioSampleMap": {
        "run_biosample_id": "runbiosample_id",
        "batch": "runbiosample_batch",
    },
    "MetabMap": {
        "metab_id": "metabfeature_id",
    },
    "ProteinMap": {
        "gene_name": "gene_symbol",
    },
    "StatsResultsMeta": {
        "experiment_name": "experiment_linkage",
        "ES": "ES_1",
        "ESPval": "ES_1_Pval",
        "ESadjPval": "ES_1_adjPval",
        "ES2": "ES_2",
    },
}


def restore_pounce_issue(issue, original_sheets: dict[str, str], version: int) -> None:
    """Translate normalized validation locations back to submitted workbook labels."""
    if version != 3:
        return
    canonical_sheet = issue.sheet
    if canonical_sheet:
        replacements = dict(_V3_ORIGINAL_KEYS.get(canonical_sheet, {}))
        if canonical_sheet == "BioSampleMap":
            for key in (issue.field, issue.column):
                if isinstance(key, str):
                    match = re.fullmatch(r"exposure(\d+)_(time|time_unit)", key)
                    if match:
                        suffix = "duration" if match.group(2) == "time" else "duration_unit"
                        replacements[key] = f"exposure{match.group(1)}_{suffix}"
        if issue.field in replacements:
            issue.field = replacements[issue.field]
        if issue.column in replacements:
            issue.column = replacements[issue.column]
        map_prefix = f"'{canonical_sheet}' maps '"
        mapped_key_end = issue.message.find("' → '", len(map_prefix))
        if issue.message.startswith(map_prefix) and mapped_key_end != -1:
            # Only the standard key is normalized. The submitter's column name
            # after the arrow must remain exactly as it appeared in the file.
            key_start = len(map_prefix)
            key = issue.message[key_start:mapped_key_end]
            issue.message = (
                issue.message[:key_start]
                + replacements.get(key, key)
                + issue.message[mapped_key_end:]
            )
        else:
            for normalized, original in replacements.items():
                issue.message = re.sub(
                    rf"(?<!\w){re.escape(normalized)}(?!\w)", original, issue.message
                )
        issue.sheet = original_sheets.get(canonical_sheet, canonical_sheet)
    for normalized, original in original_sheets.items():
        if normalized != original:
            issue.message = issue.message.replace(normalized, original)


def load_pounce_workbook(file_path: str, role: WorkbookRole) -> PounceWorkbook:
    """Load a workbook and normalize v3 layout to the shared parser contract."""
    parser = ExcelsheetParser(file_path=file_path)
    return normalize_pounce_workbook(parser, role)


def normalize_pounce_workbook(parser: ExcelsheetParser, role: WorkbookRole) -> PounceWorkbook:
    """Normalize an already loaded workbook without changing its source file."""
    file_path = parser.file_path
    version = detect_pounce_version(parser, role)
    parser.pounce_version = version
    original_sheets = {name: name for name in parser.sheet_dfs}
    if version == 2:
        return PounceWorkbook(parser, role, version, original_sheets)

    # The v3 templates use a new label for the submitter's mapped column.
    for name, frame in parser.sheet_dfs.items():
        if name.endswith("Map") or name == "EffectSize_Map":
            if _V3_MAP_HEADER in frame.columns:
                frame.rename(columns={_V3_MAP_HEADER: _V2_MAP_HEADER}, inplace=True)

    if role == "project":
        if "BiospecimenMap" in parser.sheet_dfs:
            if "BioSampleMap" in parser.sheet_dfs or "BioSampleMeta" in parser.sheet_dfs:
                raise UnsupportedPounceFormat(
                    f"{file_path}: both Biospecimen and BioSample sheets are present"
                )
            for old, new in (("BiospecimenMap", "BioSampleMap"),
                             ("BiospecimenMeta", "BioSampleMeta")):
                if old not in parser.sheet_dfs:
                    raise UnsupportedPounceFormat(f"{file_path}: missing required sheet '{old}'")
                parser.sheet_dfs[new] = parser.sheet_dfs.pop(old)
                original_sheets[new] = old
        _rename_keys(parser, "ProjectMeta", {
            "ncats_lab_groups": "lab_groups",
        })
        frame = parser.sheet_dfs.get("BioSampleMap")
        if frame is not None:
            key_column = ExcelsheetParser.KEY_COLUMN
            frame[key_column] = frame[key_column].replace({
                r"^exposure(\d+)_duration$": r"exposure\1_time",
                r"^exposure(\d+)_duration_unit$": r"exposure\1_time_unit",
            }, regex=True)
    elif role == "experiment":
        _rename_keys(parser, "RunBioSampleMap", {
            "runbiosample_id": "run_biosample_id",
            "runbiosample_batch": "batch",
        })
        if "ProteinMap" in parser.sheet_dfs and "PeakData" in parser.sheet_dfs:
            if "MetabMap" in parser.sheet_dfs:
                raise UnsupportedPounceFormat(
                    f"{file_path}: PeakData has both protein and metabolite mappings"
                )
            if "ProteinData" in parser.sheet_dfs or "ProteinDataMeta" in parser.sheet_dfs:
                raise UnsupportedPounceFormat(
                    f"{file_path}: both PeakData and ProteinData sheets are present"
                )
            for old, new in (("PeakDataMeta", "ProteinDataMeta"),
                             ("PeakData", "ProteinData")):
                if old in parser.sheet_dfs:
                    parser.sheet_dfs[new] = parser.sheet_dfs.pop(old)
                    original_sheets[new] = old
            _rename_keys(parser, "ProteinDataMeta", {"peakdata_tag": "proteindata_tag"})
        if "MetabMap" in parser.sheet_dfs:
            _rename_keys(parser, "MetabMap", {"metabfeature_id": "metab_id"})
        if "ProteinMap" in parser.sheet_dfs:
            _rename_keys(parser, "ProteinMap", {"gene_symbol": "gene_name"})
    elif role == "stats":
        _rename_keys(parser, "StatsResultsMeta", {
            "experiment_linkage": "experiment_name",
            "ES_1": "ES",
            "ES_1_Pval": "ESPval",
            "ES_1_adjPval": "ESadjPval",
            "ES_2": "ES2",
        })
    return PounceWorkbook(parser, role, version, original_sheets)


def _rename_keys(parser: ExcelsheetParser, sheet: str, replacements: dict[str, str]) -> None:
    if sheet not in parser.sheet_dfs:
        return
    frame = parser.sheet_dfs[sheet]
    key_column = ExcelsheetParser.KEY_COLUMN
    if key_column in frame.columns:
        keys = set(frame[key_column].dropna())
        for old, new in replacements.items():
            if old in keys and new in keys:
                raise UnsupportedPounceFormat(
                    f"{parser.file_path} / {sheet}: both '{old}' and '{new}' are present"
                )
        frame[key_column] = frame[key_column].replace(replacements)
