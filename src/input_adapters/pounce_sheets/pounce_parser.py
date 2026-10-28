"""Parse POUNCE Excel workbooks into typed parsed objects.

``PounceParser`` uses ``sheet_field`` metadata on the parsed dataclasses to
drive Excel extraction.  It is shared by both the validation workflow (metadata
only, no DB) and the ETL pipeline (which additionally processes data matrices).
"""

import re
from collections import defaultdict
import pandas as pd
from typing import Dict, List, Optional, Set, Type

from src.core.validator import ValidationError
from src.input_adapters.excel_sheet_adapter import ExcelsheetParser
from src.input_adapters.pounce_sheets.parsed_classes import (
    ParsedBiosample,
    ParsedBiospecimen,
    ParsedDemographics,
    ParsedExperiment,
    ParsedExposure,
    ParsedGene,
    ParsedMetab,
    ParsedPeakDataMeta,
    ParsedPerson,
    ParsedProtein,
    ParsedProteinDataMeta,
    ParsedProject,
    ParsedRawDataMeta,
    ParsedRunBiosample,
    ParsedStatsResultsMeta,
)
from src.input_adapters.pounce_sheets.constants import ExperimentWorkbook, ProjectWorkbook, StatsResultsWorkbook
from src.input_adapters.pounce_sheets.parsed_pounce_data import ParsedPounceData
from src.input_adapters.pounce_sheets.sheet_field import get_sheet_fields
from src.input_adapters.pounce_sheets.workbook_format import load_pounce_workbook, restore_pounce_issue

# Sheet sets for structural validation, derived from constants.
_PROJECT_REQUIRED_SHEETS: Set[str] = {
    ProjectWorkbook.ProjectSheet.name,
    ProjectWorkbook.BiosampleMapSheet.name,
    ProjectWorkbook.BiosampleMetaSheet.name,
}

_EXPERIMENT_REQUIRED_SHEETS: Set[str] = {
    ExperimentWorkbook.ExperimentSheet.name,
    ExperimentWorkbook.RunBioSampleMapSheet.name,
    ExperimentWorkbook.RunBioSampleMetaSheet.name,
}
_EXPERIMENT_RECOGNIZED_SHEETS: Set[str] = _EXPERIMENT_REQUIRED_SHEETS | {
    ExperimentWorkbook.GeneMapSheet.name,
    ExperimentWorkbook.GeneMetaSheet.name,
    ExperimentWorkbook.ProteinMapSheet.name,
    ExperimentWorkbook.ProteinMetaSheet.name,
    ExperimentWorkbook.ProteinDataMetaSheet.name,
    ExperimentWorkbook.ProteinDataSheet.name,
    ExperimentWorkbook.MetabMapSheet.name,
    ExperimentWorkbook.MetabMetaSheet.name,
    ExperimentWorkbook.RawDataMetaSheet.name,
    ExperimentWorkbook.RawDataSheet.name,
    ExperimentWorkbook.PeakDataMetaSheet.name,
    ExperimentWorkbook.PeakDataSheet.name,
}

_STATS_REQUIRED_SHEETS: Set[str] = {
    StatsResultsWorkbook.StatsResultsMetaSheet.name,
}
_STATS_RECOGNIZED_SHEETS: Set[str] = _STATS_REQUIRED_SHEETS | {
    StatsResultsWorkbook.StatsReadyDataSheet.name,
    StatsResultsWorkbook.EffectSizeMapSheet.name,
    StatsResultsWorkbook.EffectSizeSheet.name,
}

# Build a mapping of sheet_name -> set of known NCATSDPI keys (including indexed templates).
# Used to warn when a submitted file contains unrecognized variable names.
_KNOWN_KEYS_BY_SHEET: Dict[str, Set[str]] = defaultdict(set)
for _cls in [
    ParsedProject, ParsedBiosample, ParsedBiospecimen, ParsedDemographics,
    ParsedExposure, ParsedExperiment, ParsedGene, ParsedMetab,
    ParsedProtein, ParsedPeakDataMeta, ParsedProteinDataMeta, ParsedRawDataMeta,
    ParsedRunBiosample, ParsedStatsResultsMeta,
]:
    for _f, _meta in get_sheet_fields(_cls):
        if _meta["sheet"]:
            _KNOWN_KEYS_BY_SHEET[_meta["sheet"]].add(_meta["key"])
            for _alias in _meta.get("aliases", []):
                _KNOWN_KEYS_BY_SHEET[_meta["sheet"]].add(_alias)

# EffectSize_Map has no parsed dataclass (parsing not yet implemented); seed manually.
_KNOWN_KEYS_BY_SHEET[StatsResultsWorkbook.EffectSizeMapSheet.name] = {
    StatsResultsWorkbook.EffectSizeMapSheet.Key.gene_id,
    StatsResultsWorkbook.EffectSizeMapSheet.Key.metabolite_id,
    StatsResultsWorkbook.EffectSizeMapSheet.Key.protein_id,
    "ES_result{}", "Pval_result{}", "AdjPval_result{}",
}

# Placeholder key used in EffectSize_Map rows that don't yet have an NCATSDPI name.
# Treated as NA: silently skipped in unknown-key checks and column-existence checks.
_PLACEHOLDER_KEY = "NA (submitter creates new values)"


def _restore_issue_sheets(issues, workbook) -> None:
    """Show submitted sheet names after checking the normalized workbook."""
    for issue in issues:
        restore_pounce_issue(issue, workbook.original_sheets, workbook.version)


class PounceParser:
    """Reads POUNCE Excel workbooks into :class:`ParsedPounceData`."""

    # ------------------------------------------------------------------
    # Meta-sheet parsing (single-row key-value sheets like ProjectMeta)
    # ------------------------------------------------------------------

    @staticmethod
    def _parse_meta_field(parser: ExcelsheetParser, sheet: str, key: str, parse_type: str):
        """Read a single value from a meta sheet using the appropriate parse strategy."""
        if parse_type == "string":
            return parser.safe_get_string(sheet, key)
        if parse_type == "string_list":
            try:
                return parser.get_one_string_list(sheet, key)
            except (LookupError, KeyError):
                return []
        if parse_type == "date":
            try:
                return parser.get_one_date(sheet, key).date()
            except (LookupError, KeyError, ValueError):
                return None
        if parse_type == "int":
            raw = parser.safe_get_string(sheet, key)
            if raw is None:
                return None
            try:
                return int(float(raw))
            except (ValueError, TypeError):
                return None
        if parse_type == "float":
            raw = parser.safe_get_string(sheet, key)
            if raw is None:
                return None
            try:
                return float(raw)
            except (ValueError, TypeError):
                return None
        if parse_type == "bool":
            raw = parser.safe_get_string(sheet, key)
            if raw is None:
                return None
            return raw.lower() in ("true", "yes", "1")
        if parse_type == "category":
            return parser.safe_get_string(sheet, key)
        return parser.safe_get_string(sheet, key)

    @classmethod
    def _parse_meta_field_with_aliases(
        cls,
        parser: ExcelsheetParser,
        sheet: str,
        key: str,
        aliases: list[str],
        parse_type: str
    ):
        candidate_keys = [key, *aliases]
        for candidate in candidate_keys:
            value = cls._parse_meta_field(parser, sheet, candidate, parse_type)
            if value is not None and value != "" and value != []:
                return value
        return None

    @classmethod
    def parse_meta_sheet(cls, parser: ExcelsheetParser, dataclass_type: Type):
        """Parse a meta sheet into an instance of *dataclass_type*.

        Iterates over every ``sheet_field``-annotated field, reads the value
        from the corresponding sheet/key, and constructs the object.
        """
        kwargs = {}
        for f, meta in get_sheet_fields(dataclass_type):
            sheet = meta["sheet"]
            key = meta["key"]
            aliases = meta.get("aliases", [])
            parse_type = meta["parse"]
            indexed = meta["indexed"]

            if indexed:
                # Indexed fields like attached_file_{}: collect 1, 2, 3, ...
                values = []
                idx = 1
                while True:
                    actual_key = key.format(idx)
                    alias_keys = [alias.format(idx) for alias in aliases]
                    val = cls._parse_meta_field_with_aliases(parser, sheet, actual_key, alias_keys, parse_type)
                    if val is None or val == "" or val == []:
                        break
                    values.append(val)
                    idx += 1
                kwargs[f.name] = values
            else:
                kwargs[f.name] = cls._parse_meta_field_with_aliases(parser, sheet, key, aliases, parse_type)

        return dataclass_type(**kwargs)

    # ------------------------------------------------------------------
    # Mapped-sheet parsing (multi-row sheets like BioSampleMeta)
    # ------------------------------------------------------------------

    @staticmethod
    def _get_column_name(param_map: dict, key: str, index: int = None, aliases: Optional[List[str]] = None) -> Optional[str]:
        """Resolve a column name from the parameter map, optionally with an index."""
        candidate_keys = [key, *(aliases or [])]
        for candidate in candidate_keys:
            actual_key = candidate.format(index) if index is not None else candidate
            col = param_map.get(actual_key)
            if col not in ("", "NA", "N/A", None):
                return col
        return None

    @staticmethod
    def _get_row_value(row, column_name: Optional[str], is_list: bool = False):
        """Read a cell value from a DataFrame row."""
        if column_name is None or column_name == "":
            return None
        val = row.get(column_name)
        if pd.isna(val) or val == "" or val == "NA" or val == "N/A":
            return None
        if is_list:
            return [v.strip() for v in str(val).split("|") if v.strip()]
        return str(val).strip()

    @classmethod
    def _parse_mapped_row(cls, row, param_map: dict, dataclass_type: Type, index: int = None):
        """Parse one data row into an instance of *dataclass_type*.

        For indexed classes (e.g. ParsedExposure), *index* is the slot number.
        """
        kwargs = {}
        for f, meta in get_sheet_fields(dataclass_type):
            key = meta["key"]
            aliases = meta.get("aliases", [])
            parse_type = meta["parse"]
            indexed = meta["indexed"]

            if indexed:
                if index is not None:
                    col = cls._get_column_name(param_map, key, index, aliases)
                else:
                    col = cls._get_column_name(param_map, key, aliases=aliases)
            else:
                col = cls._get_column_name(param_map, key, aliases=aliases)

            is_list = parse_type == "string_list"
            raw = cls._get_row_value(row, col, is_list=is_list)

            if parse_type == "int" and raw is not None:
                try:
                    raw = int(float(raw))
                except (ValueError, TypeError):
                    raw = None
            elif parse_type == "float" and raw is not None:
                try:
                    raw = float(raw)
                except (ValueError, TypeError):
                    raw = None

            kwargs[f.name] = raw

        return dataclass_type(**kwargs)

    # ------------------------------------------------------------------
    # High-level workbook parsers
    # ------------------------------------------------------------------

    def parse_project(self, parser: ExcelsheetParser):
        """Parse the project workbook (metadata only).

        Returns ``(ParsedPounceData, List[ValidationError])``.
        """
        data = ParsedPounceData()
        issues: List[ValidationError] = []

        issues.extend(self._check_sheets(parser, _PROJECT_REQUIRED_SHEETS, _PROJECT_REQUIRED_SHEETS))

        # --- ProjectMeta (meta sheet) ---
        issues.extend(self._check_unknown_keys(parser, "ProjectMeta", "ProjectMeta"))
        data.project = self.parse_meta_sheet(parser, ParsedProject)

        # --- People (from project meta sheet) ---
        project = data.project
        owners = self._build_persons(
            project.owner_names or [], project.owner_emails or [], "Owner")
        if len(owners) > 0:
            owners[0].roles.append("Contact")

        collaborators = self._build_persons(
            project.collaborator_names or [], project.collaborator_emails or [], "Collaborator")
        data.people = owners + collaborators

        # --- BioSampleMeta via BioSampleMap (mapped sheet) ---
        if "BioSampleMap" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "BioSampleMap", "BioSampleMeta"))
        if "BioSampleMap" in parser.sheet_dfs and "BioSampleMeta" in parser.sheet_dfs:
            biosample_map = parser.get_parameter_map("BioSampleMap")
            issues.extend(self._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", biosample_map))
            sheet_df = parser.sheet_dfs["BioSampleMeta"]

            # Determine how many exposure slots exist
            exposure_indices = self._detect_exposure_indices(biosample_map)

            for _, row in sheet_df.iterrows():
                data.biosamples.append(
                    self._parse_mapped_row(row, biosample_map, ParsedBiosample))
                data.biospecimens.append(
                    self._parse_mapped_row(row, biosample_map, ParsedBiospecimen))
                data.demographics.append(
                    self._parse_demographics_row(row, biosample_map))

                for exp_idx in exposure_indices:
                    parsed_exp = self._parse_mapped_row(
                        row, biosample_map, ParsedExposure, index=exp_idx)
                    parsed_exp.exposure_index = exp_idx
                    data.exposures.append(parsed_exp)

            data.param_maps[ProjectWorkbook.BiosampleMapSheet.name] = biosample_map

        return data, issues

    def parse_experiment(self, parser: ExcelsheetParser,
                          valid_biosample_ids: Optional[Set[str]] = None):
        """Parse an experiment workbook (metadata only, skips data matrices).

        ``valid_biosample_ids`` is the set of biosample IDs from the project's
        BioSampleMeta.  When provided, every biosample_id in RunBioSampleMeta
        is checked against that set.

        Returns ``(ParsedPounceData, List[ValidationError])``.
        """
        data = ParsedPounceData()
        issues: List[ValidationError] = []

        issues.extend(self._check_sheets(parser, _EXPERIMENT_REQUIRED_SHEETS, _EXPERIMENT_RECOGNIZED_SHEETS))

        # --- ExperimentMeta (meta sheet) ---
        issues.extend(self._check_unknown_keys(parser, "ExperimentMeta", "ExperimentMeta"))
        data.experiments.append(self.parse_meta_sheet(parser, ParsedExperiment))

        # --- RunBioSampleMeta via RunBioSampleMap (mapped sheet) ---
        if "RunBioSampleMap" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "RunBioSampleMap", "RunBioSampleMeta"))
            run_sample_map = parser.get_parameter_map("RunBioSampleMap")
            data.param_maps[ExperimentWorkbook.RunBioSampleMapSheet.name] = run_sample_map
        if "RunBioSampleMap" in parser.sheet_dfs and "RunBioSampleMeta" in parser.sheet_dfs:
            issues.extend(self._check_mapped_columns(parser, "RunBioSampleMap", "RunBioSampleMeta", run_sample_map))
            for _, row in parser.sheet_dfs["RunBioSampleMeta"].iterrows():
                data.run_biosamples.append(
                    self._parse_mapped_row(row, run_sample_map, ParsedRunBiosample))

        # Cross-reference: every biosample_id in RunBioSampleMeta must exist in BioSampleMeta.
        if valid_biosample_ids is not None and "RunBioSampleMap" in parser.sheet_dfs:
            biosample_id_col = run_sample_map.get(ProjectWorkbook.BiosampleMapSheet.Key.biosample_id)
            if biosample_id_col and biosample_id_col not in ("NA", "N/A", ""):
                issues.extend(self._check_data_matrix_references(
                    parser, ExperimentWorkbook.RunBioSampleMetaSheet.name,
                    biosample_id_col, valid_biosample_ids,
                    "biosample_id '{value}' in RunBioSampleMeta does not exist in BioSampleMeta",
                ))

        # --- GeneMeta via GeneMap (mapped sheet) ---
        if "GeneMap" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "GeneMap", "GeneMeta"))
            gene_map = parser.get_parameter_map("GeneMap")
            data.param_maps[ExperimentWorkbook.GeneMapSheet.name] = gene_map
        if "GeneMap" in parser.sheet_dfs and "GeneMeta" in parser.sheet_dfs:
            issues.extend(self._check_mapped_columns(parser, "GeneMap", "GeneMeta", gene_map))
            for _, row in parser.sheet_dfs["GeneMeta"].iterrows():
                data.genes.append(
                    self._parse_mapped_row(row, gene_map, ParsedGene))

        # --- ProteinMeta via ProteinMap (mapped sheet) ---
        if "ProteinMap" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "ProteinMap", "ProteinMeta"))
            protein_map = parser.get_parameter_map("ProteinMap")
            data.param_maps[ExperimentWorkbook.ProteinMapSheet.name] = protein_map
        if "ProteinMap" in parser.sheet_dfs and "ProteinMeta" in parser.sheet_dfs:
            issues.extend(self._check_mapped_columns(parser, "ProteinMap", "ProteinMeta", protein_map))
            for _, row in parser.sheet_dfs["ProteinMeta"].iterrows():
                data.proteins.append(
                    self._parse_mapped_row(row, protein_map, ParsedProtein))

        # --- MetabMeta via MetabMap (mapped sheet) ---
        if "MetabMap" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "MetabMap", "MetabMeta"))
            metab_map = parser.get_parameter_map("MetabMap")
            data.param_maps[ExperimentWorkbook.MetabMapSheet.name] = metab_map
        if "MetabMap" in parser.sheet_dfs and "MetabMeta" in parser.sheet_dfs:
            issues.extend(self._check_mapped_columns(parser, "MetabMap", "MetabMeta", metab_map))
            for _, row in parser.sheet_dfs["MetabMeta"].iterrows():
                metabolite = self._parse_mapped_row(row, metab_map, ParsedMetab)
                metabolite.metab_id_is_feature = getattr(parser, "pounce_version", 2) == 3
                metabolite.source_file = getattr(parser, "file_path", None)
                data.metabolites.append(metabolite)

        # --- Data matrix ID cross-references (per-experiment, before merging) ---
        # Check that every gene/metabolite ID in the data matrix exists in the meta sheet.
        if "GeneMap" in parser.sheet_dfs and data.genes:
            gene_id_col = gene_map.get(ExperimentWorkbook.GeneMapSheet.Key.gene_id)
            if gene_id_col and gene_id_col not in ("NA", "N/A", ""):
                valid_gene_ids = {g.gene_id for g in data.genes if g.gene_id}
                issues.extend(self._check_data_matrix_references(
                    parser, ExperimentWorkbook.RawDataSheet.name,
                    gene_id_col, valid_gene_ids,
                    f"gene_id '{{value}}' in RawData does not exist in GeneMeta",
                ))

        if "MetabMap" in parser.sheet_dfs and data.metabolites:
            metab_id_col = metab_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_id)
            if metab_id_col and metab_id_col not in ("NA", "N/A", ""):
                valid_metab_ids = {m.metab_id for m in data.metabolites if m.metab_id}
                issues.extend(self._check_data_matrix_references(
                    parser, ExperimentWorkbook.PeakDataSheet.name,
                    metab_id_col, valid_metab_ids,
                    f"metab_id '{{value}}' in PeakData does not exist in MetabMeta",
                ))

        if "ProteinMap" in parser.sheet_dfs and data.proteins:
            protein_id_col = protein_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_id)
            if protein_id_col and protein_id_col not in ("NA", "N/A", ""):
                valid_protein_ids = {p.protein_id for p in data.proteins if p.protein_id}
                issues.extend(self._check_data_matrix_references(
                    parser, ExperimentWorkbook.ProteinDataSheet.name,
                    protein_id_col, valid_protein_ids,
                    f"protein_id '{{value}}' in ProteinData does not exist in ProteinMeta",
                ))

        # Check that data matrix sample column names match run_biosample_ids in RunBioSampleMeta.
        valid_run_ids = {rb.run_biosample_id for rb in data.run_biosamples if rb.run_biosample_id}
        if valid_run_ids:
            if "GeneMap" in parser.sheet_dfs:
                gene_id_col = gene_map.get(ExperimentWorkbook.GeneMapSheet.Key.gene_id)
                if gene_id_col and gene_id_col not in ("NA", "N/A", ""):
                    issues.extend(self._check_matrix_sample_columns(
                        parser, ExperimentWorkbook.RawDataSheet.name,
                        exclude_columns={gene_id_col},
                        valid_sample_ids=valid_run_ids,
                        message_template="column '{value}' in RawData does not match any run_biosample_id in RunBioSampleMeta",
                    ))
            if "MetabMap" in parser.sheet_dfs:
                metab_id_col = metab_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_id)
                if metab_id_col and metab_id_col not in ("NA", "N/A", ""):
                    issues.extend(self._check_matrix_sample_columns(
                        parser, ExperimentWorkbook.PeakDataSheet.name,
                        exclude_columns={metab_id_col},
                        valid_sample_ids=valid_run_ids,
                        message_template="column '{value}' in PeakData does not match any run_biosample_id in RunBioSampleMeta",
                    ))
            if "ProteinMap" in parser.sheet_dfs:
                protein_id_col = protein_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_id)
                if protein_id_col and protein_id_col not in ("NA", "N/A", ""):
                    issues.extend(self._check_matrix_sample_columns(
                        parser, ExperimentWorkbook.ProteinDataSheet.name,
                        exclude_columns={protein_id_col},
                        valid_sample_ids=valid_run_ids,
                        message_template="column '{value}' in ProteinData does not match any run_biosample_id in RunBioSampleMeta",
                    ))

        # --- PeakDataMeta (meta sheet) ---
        if "PeakDataMeta" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "PeakDataMeta", "PeakDataMeta"))
            data.peak_data_meta.append(self.parse_meta_sheet(parser, ParsedPeakDataMeta))

        # --- RawDataMeta (meta sheet) ---
        if "RawDataMeta" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "RawDataMeta", "RawDataMeta"))
            data.raw_data_meta.append(self.parse_meta_sheet(parser, ParsedRawDataMeta))

        if "ProteinDataMeta" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "ProteinDataMeta", "ProteinDataMeta"))
            data.protein_data_meta.append(self.parse_meta_sheet(parser, ParsedProteinDataMeta))

        return data, issues

    def parse_stats_results(self, parser: ExcelsheetParser,
                             valid_gene_ids: Optional[Set[str]] = None,
                             valid_metab_ids: Optional[Set[str]] = None,
                             valid_protein_ids: Optional[Set[str]] = None,
                             valid_run_biosample_ids: Optional[Set[str]] = None,
                             expected_analyte: Optional[str] = None):
        """Parse a stats results workbook (metadata only).

        ``valid_gene_ids``, ``valid_metab_ids``, ``valid_protein_ids``, and ``valid_run_biosample_ids``
        are ID sets from the paired experiment.  When provided, ID values in
        StatsReadyData and EffectSize are checked against them, and StatsReadyData
        sample column names are checked against run_biosample_ids.

        Returns ``(ParsedPounceData, List[ValidationError])``.
        """
        data = ParsedPounceData()
        issues: List[ValidationError] = []

        issues.extend(self._check_sheets(parser, _STATS_REQUIRED_SHEETS, _STATS_RECOGNIZED_SHEETS))

        if "StatsResultsMeta" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "StatsResultsMeta", "StatsResultsMeta"))
            data.stats_results.append(self.parse_meta_sheet(parser, ParsedStatsResultsMeta))
            if getattr(parser, "pounce_version", 2) == 2:
                required_stats_keys = ("ES", "ESPval", "ESadjPval")
            else:
                required_stats_keys = ("ES", "ES_1_prefix")
            for key in required_stats_keys:
                if not parser.safe_get_string("StatsResultsMeta", key):
                    issues.append(ValidationError(
                        severity="error", entity="stats_results", field=key,
                        sheet="StatsResultsMeta", source_file=parser.file_path,
                        message=f"StatsResultsMeta is missing required value: {key}",
                    ))

        # --- EffectSize_Map (map sheet) ---
        if "EffectSize_Map" in parser.sheet_dfs:
            issues.extend(self._check_unknown_keys(parser, "EffectSize_Map", "EffectSize_Map"))
            effect_size_map = parser.get_parameter_map("EffectSize_Map")
            data.param_maps[StatsResultsWorkbook.EffectSizeMapSheet.name] = effect_size_map
            if expected_analyte:
                analyte_keys = {
                    "gene": "gene_id", "metabolite": "metabolite_id", "protein": "protein_id",
                }
                active_keys = {
                    family for family, key in analyte_keys.items()
                    if effect_size_map.get(key) not in (None, "", "NA", "N/A")
                }
                if expected_analyte not in active_keys:
                    issues.append(ValidationError(
                        severity="error", entity="parse", field="EffectSize_Map",
                        sheet="EffectSize_Map", source_file=parser.file_path,
                        message=(f"StatsResults analyte mapping does not match its paired "
                                 f"{expected_analyte} experiment"),
                    ))
            # Validate that every mapped column name exists in the EffectSize sheet.
            # EffectSize_Map values may be pipe-delimited ("col_name | property");
            # only the first segment is the actual column name.
            # We read the raw DataFrame (not param_map) because duplicate placeholder
            # keys collapse in to_dict(), losing most rows.
            issues.extend(self._check_effect_size_columns(parser))

            # Cross-reference IDs in StatsReadyData and EffectSize against the
            # paired experiment's meta sheets.  The analyte ID column names come
            # from EffectSize_Map (same map applies to both stats sheets).
            _stats_data_sheets = [
                StatsResultsWorkbook.StatsReadyDataSheet.name,
                StatsResultsWorkbook.EffectSizeSheet.name,
            ]
            if valid_gene_ids:
                gene_id_col = effect_size_map.get(StatsResultsWorkbook.EffectSizeMapSheet.Key.gene_id)
                if gene_id_col and gene_id_col not in ("NA", "N/A", ""):
                    for sheet in _stats_data_sheets:
                        issues.extend(self._check_data_matrix_references(
                            parser, sheet, gene_id_col, valid_gene_ids,
                            f"gene_id '{{value}}' in {sheet} does not exist in paired experiment's GeneMeta",
                        ))

            if valid_metab_ids:
                metab_id_col = effect_size_map.get(StatsResultsWorkbook.EffectSizeMapSheet.Key.metabolite_id)
                if metab_id_col and metab_id_col not in ("NA", "N/A", ""):
                    for sheet in _stats_data_sheets:
                        issues.extend(self._check_data_matrix_references(
                            parser, sheet, metab_id_col, valid_metab_ids,
                            f"metab_id '{{value}}' in {sheet} does not exist in paired experiment's MetabMeta",
                        ))

            if valid_protein_ids:
                protein_id_col = effect_size_map.get(StatsResultsWorkbook.EffectSizeMapSheet.Key.protein_id)
                if protein_id_col and protein_id_col not in ("NA", "N/A", ""):
                    for sheet in _stats_data_sheets:
                        issues.extend(self._check_data_matrix_references(
                            parser, sheet, protein_id_col, valid_protein_ids,
                            f"protein_id '{{value}}' in {sheet} does not exist in paired experiment's ProteinMeta",
                        ))

            # Check StatsReadyData sample column names against paired run_biosample_ids.
            # Exclude whichever analyte ID column(s) are configured in EffectSize_Map.
            if valid_run_biosample_ids:
                analyte_cols = set()
                _g = effect_size_map.get(StatsResultsWorkbook.EffectSizeMapSheet.Key.gene_id)
                if _g and _g not in ("NA", "N/A", ""):
                    analyte_cols.add(_g)
                _m = effect_size_map.get(StatsResultsWorkbook.EffectSizeMapSheet.Key.metabolite_id)
                if _m and _m not in ("NA", "N/A", ""):
                    analyte_cols.add(_m)
                _p = effect_size_map.get(StatsResultsWorkbook.EffectSizeMapSheet.Key.protein_id)
                if _p and _p not in ("NA", "N/A", ""):
                    analyte_cols.add(_p)
                issues.extend(self._check_matrix_sample_columns(
                    parser, StatsResultsWorkbook.StatsReadyDataSheet.name,
                    exclude_columns=analyte_cols,
                    valid_sample_ids=valid_run_biosample_ids,
                    message_template="column '{value}' in StatsReadyData does not match any run_biosample_id in RunBioSampleMeta",
                ))

        return data, issues

    def parse_all(self, project_file: str,
                  experiment_files: List[str] = None,
                  stats_files: List[str] = None,
                  workbook_loader=None):
        """Parse all workbooks into a single :class:`ParsedPounceData`.

        Returns ``(ParsedPounceData, List[ValidationError])``.
        """
        experiment_files = experiment_files or []
        stats_files = stats_files or []
        workbook_loader = workbook_loader or load_pounce_workbook

        project_workbook = workbook_loader(project_file, "project")
        combined, issues = self.parse_project(project_workbook.parser)
        combined.workbook_versions[project_file] = project_workbook.version
        combined.original_sheets[project_file] = project_workbook.original_sheets
        _restore_issue_sheets(issues, project_workbook)

        # Build the valid biosample ID set from the project (available for all experiments).
        valid_biosample_ids = {str(b.biosample_id).strip() for b in combined.biosamples if b.biosample_id} or None

        for i, exp_file in enumerate(experiment_files):
            exp_workbook = workbook_loader(exp_file, "experiment")
            exp_data, exp_issues = self.parse_experiment(
                exp_workbook.parser, valid_biosample_ids=valid_biosample_ids)
            combined.workbook_versions[exp_file] = exp_workbook.version
            combined.original_sheets[exp_file] = exp_workbook.original_sheets
            _restore_issue_sheets(exp_issues, exp_workbook)
            combined.experiments.extend(exp_data.experiments)
            combined.genes.extend(exp_data.genes)
            combined.proteins.extend(exp_data.proteins)
            combined.metabolites.extend(exp_data.metabolites)
            combined.protein_data_meta.extend(exp_data.protein_data_meta)
            combined.peak_data_meta.extend(exp_data.peak_data_meta)
            combined.raw_data_meta.extend(exp_data.raw_data_meta)
            combined.run_biosamples.extend(exp_data.run_biosamples)
            combined.param_maps.update(exp_data.param_maps)
            issues.extend(exp_issues)

            # Process the paired stats file (if any) immediately after its experiment
            # so we can pass the experiment's gene/metabolite ID sets for cross-referencing.
            if i < len(stats_files):
                valid_gene_ids = {str(g.gene_id).strip() for g in exp_data.genes if g.gene_id} or None
                valid_metab_ids = {str(m.metab_id).strip() for m in exp_data.metabolites if m.metab_id} or None
                valid_protein_ids = {str(p.protein_id).strip() for p in exp_data.proteins if p.protein_id} or None
                valid_run_biosample_ids = {str(rb.run_biosample_id).strip() for rb in exp_data.run_biosamples if rb.run_biosample_id} or None
                stats_workbook = workbook_loader(stats_files[i], "stats")
                stats_parser = stats_workbook.parser
                stats_data, stats_issues = self.parse_stats_results(
                    stats_parser,
                    valid_gene_ids=valid_gene_ids,
                    valid_metab_ids=valid_metab_ids,
                    valid_protein_ids=valid_protein_ids,
                    valid_run_biosample_ids=valid_run_biosample_ids,
                    expected_analyte=(
                        "gene" if "GeneMap" in exp_workbook.parser.sheet_dfs else
                        "protein" if "ProteinMap" in exp_workbook.parser.sheet_dfs else
                        "metabolite" if "MetabMap" in exp_workbook.parser.sheet_dfs else None
                    ),
                )
                combined.stats_results.extend(stats_data.stats_results)
                combined.param_maps.update(stats_data.param_maps)
                combined.workbook_versions[stats_files[i]] = stats_workbook.version
                combined.original_sheets[stats_files[i]] = stats_workbook.original_sheets
                _restore_issue_sheets(stats_issues, stats_workbook)
                issues.extend(stats_issues)

        # Any stats files beyond the number of experiment files (unusual) get
        # parsed without ID cross-reference checks.
        for stats_file in stats_files[len(experiment_files):]:
            stats_workbook = workbook_loader(stats_file, "stats")
            stats_parser = stats_workbook.parser
            stats_data, stats_issues = self.parse_stats_results(stats_parser)
            combined.stats_results.extend(stats_data.stats_results)
            combined.param_maps.update(stats_data.param_maps)
            combined.workbook_versions[stats_file] = stats_workbook.version
            combined.original_sheets[stats_file] = stats_workbook.original_sheets
            _restore_issue_sheets(stats_issues, stats_workbook)
            issues.extend(stats_issues)

        return combined, issues

    # ------------------------------------------------------------------
    # Helpers
    # ------------------------------------------------------------------

    @staticmethod
    def _check_sheets(parser: ExcelsheetParser, required: Set[str],
                      recognized: Set[str]) -> List[ValidationError]:
        """Return errors for missing required sheets and warnings for unrecognized extra sheets."""
        issues: List[ValidationError] = []
        actual = set(parser.sheet_dfs.keys())
        for sheet in sorted(required - actual):
            issues.append(ValidationError(
                severity="error",
                entity="parse",
                field=sheet,
                message=f"Required sheet '{sheet}' is missing from the workbook",
                source_file=parser.file_path,
            ))
        for sheet in sorted(actual - recognized):
            issues.append(ValidationError(
                severity="warning",
                entity="parse",
                field=sheet,
                message=f"Unrecognized sheet '{sheet}' — it will not be imported",
                sheet=sheet,
                source_file=parser.file_path,
            ))
        return issues

    @staticmethod
    def _check_unknown_keys(parser: ExcelsheetParser, sheet_name: str,
                             data_sheet_name: str) -> List[ValidationError]:
        """Return warnings for unrecognized NCATSDPI keys present in a sheet.

        *sheet_name* is the sheet to inspect (e.g. ``"BioSampleMap"``).
        *data_sheet_name* is the sheet whose known keys apply (e.g. ``"BioSampleMeta"``).
        """
        issues: List[ValidationError] = []
        if sheet_name not in parser.sheet_dfs:
            return issues
        df = parser.sheet_dfs[sheet_name]
        key_col = ExcelsheetParser.KEY_COLUMN
        if key_col not in df.columns:
            return issues
        known = _KNOWN_KEYS_BY_SHEET.get(data_sheet_name, set())
        for raw_key in df[key_col].dropna():
            key = str(raw_key).strip()
            if not key or key == _PLACEHOLDER_KEY:
                continue
            normalized = re.sub(r'\d+', '{}', key)
            if key not in known and normalized not in known:
                issues.append(ValidationError(
                    severity="warning",
                    entity="parse",
                    field=key,
                    message=f"Unrecognized NCATSDPI_Variable_Name: '{key}' — this field will not be used",
                    sheet=sheet_name,
                    column=key,
                    source_file=parser.file_path,
                ))
        return issues

    @staticmethod
    def _check_mapped_columns(
        parser: ExcelsheetParser,
        map_sheet: str,
        meta_sheet: str,
        param_map: dict,
    ) -> List[ValidationError]:
        """Return errors for column names configured in a map sheet that do not
        exist as headers in the corresponding meta sheet.

        A missing mapped column means every row silently yields ``None`` for
        that field, which is almost always a submitter typo.
        """
        issues: List[ValidationError] = []
        if meta_sheet not in parser.sheet_dfs:
            return issues
        meta_columns = set(parser.sheet_dfs[meta_sheet].columns)
        for ncatsdpi_key, column_name in param_map.items():
            if not column_name or column_name in ("NA", "N/A", ""):
                continue
            if column_name not in meta_columns:
                issues.append(ValidationError(
                    severity="error",
                    entity="parse",
                    field=ncatsdpi_key,
                    message=(
                        f"'{map_sheet}' maps '{ncatsdpi_key}' → '{column_name}', "
                        f"but column '{column_name}' does not exist in '{meta_sheet}'"
                    ),
                    sheet=map_sheet,
                    column=ncatsdpi_key,
                    source_file=parser.file_path,
                ))
        return issues

    @staticmethod
    def _check_matrix_sample_columns(
        parser: ExcelsheetParser,
        data_sheet: str,
        exclude_columns: Set[str],
        valid_sample_ids: Set[str],
        message_template: str,
    ) -> List[ValidationError]:
        """Check that every column header in a data matrix (excluding analyte ID columns)
        exists as a run_biosample_id in the paired experiment's RunBioSampleMeta.

        ``exclude_columns`` is the set of analyte ID column names to skip.
        ``message_template`` may contain ``{value}`` for the offending column name.
        """
        if data_sheet not in parser.sheet_dfs:
            return []
        issues = []
        for col in parser.sheet_dfs[data_sheet].columns:
            str_col = str(col).strip()
            if not str_col or str_col in exclude_columns or str_col in ("NA", "N/A"):
                continue
            if str_col not in valid_sample_ids:
                issues.append(ValidationError(
                    severity="error",
                    entity="parse",
                    field=str_col,
                    message=message_template.format(value=str_col),
                    sheet=data_sheet,
                    column=str_col,
                    source_file=parser.file_path,
                ))
        return issues

    @staticmethod
    def _check_data_matrix_references(
        parser: ExcelsheetParser,
        data_sheet: str,
        id_column: str,
        valid_ids: set,
        message_template: str,
    ) -> List[ValidationError]:
        """Check that every ID value in a data matrix column exists in a known-good set.

        Run per-experiment (before merging) so each file's data matrix is checked
        against only its own meta sheet — not the union of all experiments.
        ``id_column`` is the resolved column name (from the map sheet param_map).
        ``message_template`` may contain ``{value}`` for the offending value.
        """
        if data_sheet not in parser.sheet_dfs:
            return []
        df = parser.sheet_dfs[data_sheet]
        if id_column not in df.columns:
            return []  # _check_mapped_columns already reported this
        issues = []
        for i, val in enumerate(df[id_column]):
            if pd.isna(val):
                continue
            str_val = str(val).strip()
            if not str_val or str_val in ("NA", "N/A"):
                continue
            if str_val not in valid_ids:
                issues.append(ValidationError(
                    severity="error",
                    entity="parse",
                    field=id_column,
                    message=message_template.format(value=str_val),
                    sheet=data_sheet,
                    column=id_column,
                    row=i,
                    source_file=parser.file_path,
                ))
        return issues

    @staticmethod
    def _check_effect_size_columns(parser: ExcelsheetParser) -> List[ValidationError]:
        """Check that every column name configured in EffectSize_Map exists in EffectSize.

        Reads the raw DataFrame (not the param_map dict) because multiple rows
        share the placeholder key and would collapse to one entry in to_dict().
        Values may be pipe-delimited ("col_name | property"); only the first
        segment is the actual column name.
        """
        if "EffectSize_Map" not in parser.sheet_dfs or "EffectSize" not in parser.sheet_dfs:
            return []
        issues: List[ValidationError] = []
        es_columns = set(parser.sheet_dfs["EffectSize"].columns)
        es_map_df = parser.sheet_dfs["EffectSize_Map"]
        key_col = ExcelsheetParser.KEY_COLUMN
        val_col = ExcelsheetParser.MAPPED_VALUE_COLUMN
        if val_col not in es_map_df.columns:
            return []
        for _, row in es_map_df.iterrows():
            raw_key = row.get(key_col)
            raw_val = row.get(val_col)
            if pd.isna(raw_val):
                continue
            col_name = str(raw_val).split("|")[0].strip()
            if not col_name or col_name in ("NA", "N/A", ""):
                continue
            if col_name not in es_columns:
                ncatsdpi_key = str(raw_key).strip() if raw_key and not pd.isna(raw_key) else _PLACEHOLDER_KEY
                issues.append(ValidationError(
                    severity="error",
                    entity="parse",
                    field=ncatsdpi_key,
                    message=(
                        f"'EffectSize_Map' maps '{ncatsdpi_key}' → '{col_name}', "
                        f"but column '{col_name}' does not exist in 'EffectSize'"
                    ),
                    sheet="EffectSize_Map",
                    column=ncatsdpi_key,
                    source_file=parser.file_path,
                ))
        return issues

    @staticmethod
    def _build_persons(names: List[str], emails: List[str], role: str) -> List[ParsedPerson]:
        if not names:
            return []
        if emails and len(names) == len(emails):
            return [ParsedPerson(name=n, email=e, roles=[role]) for n, e in zip(names, emails)]
        return [ParsedPerson(name=n, roles=[role]) for n in names]

    @classmethod
    def _detect_exposure_indices(cls, param_map: dict) -> List[int]:
        """Return the list of valid exposure slot numbers (1, 2, ...) from the map."""
        indices = []
        idx = 1
        while True:
            names_col = cls._get_column_name(param_map, "exposure{}_names", idx)
            type_col = cls._get_column_name(param_map, "exposure{}_type", idx)
            cat_col = cls._get_column_name(param_map, "exposure{}_category", idx)
            if all(c is None for c in [names_col, type_col, cat_col]):
                break
            indices.append(idx)
            idx += 1
        return indices

    @classmethod
    def _parse_demographics_row(cls, row, param_map: dict) -> ParsedDemographics:
        """Parse demographics fields from a biosample row.

        Demographics uses indexed fields for category columns, so we handle
        those specially here.
        """
        kwargs = {}
        for f, meta in get_sheet_fields(ParsedDemographics):
            key = meta["key"]
            parse_type = meta["parse"]
            indexed = meta["indexed"]

            if indexed:
                # Collect category1, category2, ... columns
                values = []
                idx = 1
                while True:
                    col = cls._get_column_name(param_map, key, idx)
                    if col is None:
                        break
                    val = cls._get_row_value(row, col)
                    if val is not None:
                        values.append(f"{col}:{val}")
                    idx += 1
                kwargs[f.name] = values if values else []
            else:
                col = cls._get_column_name(param_map, key)
                kwargs[f.name] = cls._get_row_value(row, col)

        return ParsedDemographics(**kwargs)
