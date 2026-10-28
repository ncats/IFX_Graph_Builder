"""Template detection and normalization at the POUNCE workbook boundary."""

from types import SimpleNamespace

import pandas as pd
import pytest

from src.input_adapters.pounce_sheets.workbook_format import (
    UnsupportedPounceFormat,
    detect_pounce_version,
    normalize_pounce_workbook,
    restore_pounce_issue,
)
from src.core.validator import ValidationError
from src.input_adapters.pounce_sheets.mapping_coverage import check_metabolite_coverage
from src.input_adapters.pounce_sheets.parsed_classes import ParsedMetab
from src.input_adapters.pounce_sheets.pounce_node_builder import PounceNodeBuilder
from src.models.pounce.experiment import Experiment


def _map(keys, header):
    return pd.DataFrame({"NCATSDPI_Variable_Name": keys, header: [f"col{i}" for i in range(len(keys))]})


def _book(sheets, filename="renamed.xlsx"):
    return SimpleNamespace(file_path=filename, sheet_dfs=sheets)


def test_version_detection_uses_content_not_filename():
    v2 = _book({"ProjectMeta": pd.DataFrame(), "BioSampleMap": _map(["biosample_id"], "Submitter_Variable_Name")}, "project_v3.xlsx")
    v3 = _book({"ProjectMeta": pd.DataFrame(), "BioSampleMap": _map(["biosample_id"], "Submitter_Column_Name")}, "project_v2.xlsx")
    assert detect_pounce_version(v2, "project") == 2
    assert detect_pounce_version(v3, "project") == 3


def test_experiment_rejects_conflicting_map_headers():
    workbook = _book({
        "ExperimentMeta": pd.DataFrame(),
        "RunBioSampleMap": _map(["runbiosample_id"], "Submitter_Column_Name"),
        "GeneMap": _map(["gene_id"], "Submitter_Variable_Name"),
    })
    with pytest.raises(UnsupportedPounceFormat, match="disagree"):
        detect_pounce_version(workbook, "experiment")


def test_v3_proteomics_layout_normalizes_to_shared_sheet_contract():
    workbook = _book({
        "ProjectMeta": pd.DataFrame({"NCATSDPI_Variable_Name": ["ncats_lab_groups"], "Submitter_Value": ["ACC"]}),
        "BiospecimenMap": _map(["biosample_id", "exposure2_duration"], "Submitter_Column_Name"),
        "BiospecimenMeta": pd.DataFrame({"col0": ["S1"], "col1": ["2 hours"]}),
    })
    normalized = normalize_pounce_workbook(workbook, "project")
    assert normalized.version == 3
    assert normalized.original_sheet("BioSampleMeta") == "BiospecimenMeta"
    assert "BioSampleMap" in workbook.sheet_dfs
    assert "Submitter_Variable_Name" in workbook.sheet_dfs["BioSampleMap"].columns
    assert workbook.sheet_dfs["BioSampleMap"]["NCATSDPI_Variable_Name"].tolist() == [
        "biosample_id", "exposure2_time",
    ]
    assert workbook.sheet_dfs["ProjectMeta"]["NCATSDPI_Variable_Name"].tolist() == ["lab_groups"]


def test_v3_metabolite_feature_and_standard_id_remain_distinct():
    workbook = _book({
        "ExperimentMeta": pd.DataFrame(),
        "RunBioSampleMap": _map(["runbiosample_id"], "Submitter_Column_Name"),
        "MetabMap": _map(["metabfeature_id", "standardized_metab_id"], "Submitter_Column_Name"),
    })
    normalize_pounce_workbook(workbook, "experiment")
    assert workbook.sheet_dfs["MetabMap"]["NCATSDPI_Variable_Name"].tolist() == [
        "metab_id", "standardized_metab_id",
    ]


def test_v3_stats_metadata_normalizes_primary_effect_size_descriptions():
    workbook = _book({
        "StatsResultsMeta": pd.DataFrame({
            "NCATSDPI_Variable_Name": ["ES_1", "ES_1_Pval", "ES_1_prefix"],
            "Submitter_Value": ["fold change", "test", "ES"],
        }),
        "EffectSize_Map": _map(["gene_id", "ES_result1"], "Submitter_Column_Name"),
    })
    normalize_pounce_workbook(workbook, "stats")
    assert workbook.sheet_dfs["StatsResultsMeta"]["NCATSDPI_Variable_Name"].tolist() == [
        "ES", "ESPval", "ES_1_prefix",
    ]


def test_v3_feature_id_is_not_treated_as_a_canonical_metabolite_id():
    class Resolver:
        def resolve_internal(self, stubs):
            return {stub.id: [stub] for stub in stubs if stub.id == "CHEBI:123"}

    metabolites = [
        ParsedMetab(metab_id="feature-1", metab_id_is_feature=True),
        ParsedMetab(metab_id="feature-2", standardized_metab_id="CHEBI:123", metab_id_is_feature=True),
        ParsedMetab(metab_id="CHEBI:123"),
    ]
    coverage = check_metabolite_coverage(metabolites, Resolver())
    assert coverage.total == 3
    assert coverage.mapped == 2
    assert coverage.unmapped_ids == ["feature-1"]


def test_v3_coverage_keeps_same_feature_id_from_different_experiments():
    class Resolver:
        def resolve_internal(self, stubs):
            return {stub.id: [stub] for stub in stubs if stub.id == "CHEBI:123"}

    metabolites = [
        ParsedMetab(metab_id="feature-1", standardized_metab_id="CHEBI:123", metab_id_is_feature=True, source_file="exp1.xlsx"),
        ParsedMetab(metab_id="feature-1", metab_id_is_feature=True, source_file="exp2.xlsx"),
    ]
    coverage = check_metabolite_coverage(metabolites, Resolver())
    assert (coverage.total, coverage.mapped, coverage.unmapped_ids) == (2, 1, ["feature-1"])


def test_v3_metabolite_graph_edge_uses_standard_id_only():
    parser = _book({"MetabMeta": pd.DataFrame({
        "feature": ["feature-1", "feature-2", "feature-3"],
        "standard": ["CHEBI:123", "NA", "N/A"],
    })})
    parser.pounce_version = 3
    parser.get_parameter_map = lambda sheet: {
        "metab_id": "feature", "standardized_metab_id": "standard",
    }
    nodes, edges = PounceNodeBuilder()._parse_metabolites(Experiment(id="experiment-1"), parser)
    assert [node.id for node in nodes] == ["experiment-1-feature-1", "experiment-1-feature-2", "experiment-1-feature-3"]
    assert [edge.end_node.id for edge in edges] == ["CHEBI:123"]


def test_v3_effect_size_map_key_defines_stat_type_for_custom_prefix():
    parser = _book({"EffectSize_Map": pd.DataFrame({
        "NCATSDPI_Variable_Name": ["gene_id", "ES_result1", "Pval_result1", "AdjPval_result1"],
        "Submitter_Variable_Name": ["gene", "logFC_case_vs_ctrl", "prob_case_vs_ctrl", "q_case_vs_ctrl"],
    })})
    parser.pounce_version = 3
    parser.safe_get_string = lambda sheet, key: {
        "ES_1_prefix": "logFC",
        "ES_1_Pval_prefix": "prob",
        "ES_1_adjPval_prefix": "q",
    }.get(key)
    analyte_id, columns = PounceNodeBuilder._parse_effect_size_map(parser)
    assert analyte_id == "gene"
    assert [columns[name].stat_type for name in ["logFC_case_vs_ctrl", "prob_case_vs_ctrl", "q_case_vs_ctrl"]] == [
        "effect_size", "pvalue", "adjusted_pvalue",
    ]
    assert columns["logFC_case_vs_ctrl"].comparison == "case_VS_ctrl"


def test_v3_issue_uses_original_sheet_and_key_names():
    issue = ValidationError(
        severity="error", entity="parse", field="run_biosample_id",
        sheet="RunBioSampleMap", column="run_biosample_id",
        message="'BioSampleMap' maps 'run_biosample_id' but column is missing in 'BioSampleMeta'",
    )
    restore_pounce_issue(issue, {"BioSampleMap": "BiospecimenMap", "BioSampleMeta": "BiospecimenMeta"}, 3)
    assert issue.field == "runbiosample_id"
    assert issue.column == "runbiosample_id"
    assert issue.message == "'BiospecimenMap' maps 'runbiosample_id' but column is missing in 'BiospecimenMeta'"
