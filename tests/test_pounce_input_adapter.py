import pandas as pd
import pytest
from unittest.mock import patch, MagicMock
from datetime import datetime

from src.input_adapters.pounce_sheets.parsed_classes import ParsedPerson, ParsedBiosample, ParsedBiospecimen, ParsedExperiment, ParsedRunBiosample, ParsedProject
from src.input_adapters.pounce_sheets.parsed_pounce_data import ParsedPounceData
from src.input_adapters.pounce_sheets.pounce_input_adapter import PounceInputAdapter
from src.input_adapters.pounce_sheets.pounce_node_builder import PounceNodeBuilder, _iter_pairs
from src.input_adapters.excel_sheet_adapter import ExcelsheetParser
from src.models.pounce.project import Project, AccessLevel
from src.models.pounce.experiment import Experiment, ExperimentBiospecimenEdge, RunBiosample
from src.models.protein import MeasuredProtein
from src.models.pounce.dataset import Dataset, DatasetProteinEdge
from src.models.pounce.workbook_artifact import WorkbookArtifact
from src.constants import DataSourceName


def test_experiment_stats_pairing_keeps_missing_middle_slot():
    assert list(_iter_pairs(["exp1", "exp2", "exp3"], ["stats1", None, "stats3"])) == [
        ("exp1", "stats1"), ("exp2", None), ("exp3", "stats3"),
    ]
    with pytest.raises(ValueError, match="use None"):
        list(_iter_pairs(["exp1", "exp2"], ["stats2"]))


def test_missing_stats_experiment_does_not_require_its_analyte_in_other_stats_map():
    from pathlib import Path

    root = Path(__file__).resolve().parents[1] / "input_files/manual/pounce/OLD/templates"
    adapter = PounceInputAdapter(
        project_file=str(root / "POUNCE_Project_Metabo_Template_v2.xlsx"),
        experiment_files=[
            str(root / "POUNCE_Experiment_Metabo_Template_v2.xlsx"),
            str(root / "POUNCE_Experiment_Proteomics_Template_v2.xlsx"),
        ],
        stats_results_files=[None, str(root / "POUNCE_StatsResults_Proteomics_Template_v2.xlsx")],
        validators_config="src/use_cases/pounce/pounce_validators.yaml",
    )
    data = adapter.get_validation_data()
    content_issues = [issue for validator in adapter.get_validators() for issue in validator.validate(data)]
    assert not any(
        issue.sheet == "EffectSize_Map" and issue.field == "metabolite_id"
        for issue in content_issues
    )


def create_mock_parser(project_data: dict, biosample_map_data: dict = None, biosample_meta_df: pd.DataFrame = None):
    """Helper to create a PounceInputAdapter with mocked parser data."""
    mock_parser = MagicMock(spec=ExcelsheetParser)

    # Mock get_one_string
    def mock_get_one_string(sheet_name, data_key):
        if sheet_name == "ProjectMeta":
            return project_data.get(data_key)
        return None
    mock_parser.get_one_string = MagicMock(side_effect=mock_get_one_string)

    # Mock get_one_string_list
    def mock_get_one_string_list(sheet_name, data_key):
        if sheet_name == "ProjectMeta":
            val = project_data.get(data_key)
            if val is None:
                return []
            if isinstance(val, list):
                return val
            return [val]
        return []
    mock_parser.get_one_string_list = MagicMock(side_effect=mock_get_one_string_list)

    # Mock get_one_date
    def mock_get_one_date(sheet_name, data_key):
        if sheet_name == "ProjectMeta" and data_key == "date":
            return datetime(2024, 3, 15)
        return None
    mock_parser.get_one_date = MagicMock(side_effect=mock_get_one_date)

    # Mock get_parameter_map
    def mock_get_parameter_map(sheet_name):
        if sheet_name == "BioSampleMap" and biosample_map_data:
            return biosample_map_data
        return {}
    mock_parser.get_parameter_map = MagicMock(side_effect=mock_get_parameter_map)

    # Mock sheet_dfs
    mock_parser.sheet_dfs = {}
    if biosample_meta_df is not None:
        mock_parser.sheet_dfs["BioSampleMeta"] = biosample_meta_df

    return mock_parser


def create_adapter_with_mock(project_data: dict, biosample_map_data: dict = None, biosample_meta_df: pd.DataFrame = None):
    """Helper to create a PounceInputAdapter with mocked parser."""
    mock_parser = create_mock_parser(project_data, biosample_map_data, biosample_meta_df)

    with patch.object(PounceInputAdapter, '__init__', lambda self, **kwargs: None):
        adapter = PounceInputAdapter()
        adapter.project_file = "dummy.xlsx"
        adapter.project_parser = mock_parser
        adapter.experiment_files = []
        adapter.stats_results_files = []

    return adapter


# --- _person_nodes_and_edges tests ---

def test_person_nodes_and_edges_with_names_and_emails():
    proj = Project(id="PROJ001", name="Test")
    people = [
        ParsedPerson(name="Alice", email="alice@example.com", roles=["Owner"]),
        ParsedPerson(name="Bob", email="bob@example.com", roles=["Owner"]),
    ]
    nodes, edges = PounceNodeBuilder._person_nodes_and_edges(proj, people)
    assert len(nodes) == 2
    assert nodes[0].id == "alice@example.com"
    assert nodes[0].email == "alice@example.com"
    assert nodes[1].id == "bob@example.com"
    assert nodes[1].email == "bob@example.com"


def test_person_nodes_and_edges_with_no_email():
    proj = Project(id="PROJ001", name="Test")
    people = [ParsedPerson(name="Alice", roles=["Owner"])]
    nodes, edges = PounceNodeBuilder._person_nodes_and_edges(proj, people)
    assert len(nodes) == 1
    assert nodes[0].id == "alice"
    assert nodes[0].email is None


def test_person_nodes_and_edges_with_empty_people():
    proj = Project(id="PROJ001", name="Test")
    nodes, edges = PounceNodeBuilder._person_nodes_and_edges(proj, [])
    assert nodes == []
    assert edges == []


def test_person_nodes_and_edges_creates_correct_edges():
    proj = Project(id="PROJ001", name="Test")
    people = [
        ParsedPerson(name="Alice", email="alice@example.com", roles=["Owner"]),
        ParsedPerson(name="Bob", email="bob@example.com", roles=["Collaborator"]),
    ]
    nodes, edges = PounceNodeBuilder._person_nodes_and_edges(proj, people)
    assert len(edges) == 2
    assert edges[0].roles == ["Owner"]
    assert edges[1].roles == ["Collaborator"]


# --- get_datasource_name tests ---

def test_get_datasource_name_returns_ncats_pounce():
    # Arrange
    adapter = create_adapter_with_mock({})

    # Act
    result = adapter.get_datasource_name()

    # Assert
    assert result == DataSourceName.NCATSPounce


def test_project_id_builder_uses_hyphens_not_colons():
    parsed_project = ParsedProject(
        project_name="Colon Test Project",
        date=datetime(2024, 3, 15),
        owner_emails=["owner@example.com"],
        collaborator_emails=["collab@example.com"],
    )

    project_id = PounceNodeBuilder._project_id(parsed_project)

    assert project_id.startswith("pounce_")
    assert ":" not in project_id


def test_experiment_id_builder_uses_hyphens_not_colons():
    parsed_experiment = ParsedExperiment(
        experiment_name="Experiment: 1",
        date=datetime(2024, 3, 16),
        platform_type="bulk_rnaseq",
        platform_name="NovaSeq",
        platform_output_type="counts",
    )

    experiment_id = PounceNodeBuilder._experiment_id("pounce-proj-1234abcd", parsed_experiment)

    assert experiment_id.startswith("pounce-proj-1234abcd_")
    assert ":" not in experiment_id


def test_experiment_nodes_emit_deduped_experiment_biospecimen_edges():
    builder = PounceNodeBuilder()
    project = Project(id="PROJ001", name="Test")
    project_data = ParsedPounceData(
        biosamples=[
            ParsedBiosample(biosample_id="BS1", biosample_type="sample"),
            ParsedBiosample(biosample_id="BS2", biosample_type="sample"),
        ],
        biospecimens=[
            ParsedBiospecimen(biospecimen_id="SPEC1", biospecimen_type="specimen"),
            ParsedBiospecimen(biospecimen_id="SPEC1", biospecimen_type="specimen"),
        ],
    )
    list(builder._biosample_nodes(project, project_data))

    exp_data = ParsedPounceData(
        experiments=[ParsedExperiment(experiment_name="Experiment 1")],
        run_biosamples=[
            ParsedRunBiosample(run_biosample_id="RB1", biosample_id="BS1"),
            ParsedRunBiosample(run_biosample_id="RB2", biosample_id="BS2"),
        ],
    )
    exp_parser = MagicMock(spec=ExcelsheetParser)
    exp_parser.file_path = "experiment.xlsx"
    exp_parser.sheet_dfs = {}

    with patch.object(PounceNodeBuilder, "_compute_experiment_counts", return_value=(0, 0)):
        batches = list(builder._experiment_nodes(project, exp_data, exp_parser))

    experiment_biospecimen_edges = [
        edge
        for batch in batches
        for edge in batch
        if isinstance(edge, ExperimentBiospecimenEdge)
    ]

    assert len(experiment_biospecimen_edges) == 1
    assert experiment_biospecimen_edges[0].start_node.id.startswith("PROJ001_experiment-1-")
    assert experiment_biospecimen_edges[0].end_node.id == "PROJ001-SPEC1"


def test_experiment_nodes_propagate_run_biosample_batch():
    builder = PounceNodeBuilder()
    project = Project(id="PROJ001", name="Test")
    project_data = ParsedPounceData(
        biosamples=[ParsedBiosample(biosample_id="BS1", biosample_type="sample")],
        biospecimens=[ParsedBiospecimen(biospecimen_id="SPEC1", biospecimen_type="specimen")],
    )
    list(builder._biosample_nodes(project, project_data))

    exp_data = ParsedPounceData(
        experiments=[ParsedExperiment(experiment_name="Experiment 1")],
        run_biosamples=[ParsedRunBiosample(run_biosample_id="RB1", biosample_id="BS1", batch="3")],
    )
    exp_parser = MagicMock(spec=ExcelsheetParser)
    exp_parser.file_path = "experiment.xlsx"
    exp_parser.sheet_dfs = {}

    with patch.object(PounceNodeBuilder, "_compute_experiment_counts", return_value=(0, 0)):
        batches = list(builder._experiment_nodes(project, exp_data, exp_parser))

    run_biosamples = [
        node
        for batch in batches
        for node in batch
        if isinstance(node, RunBiosample)
    ]

    assert len(run_biosamples) == 1
    assert run_biosamples[0].batch == "3"


def test_experiment_nodes_inherit_project_access():
    builder = PounceNodeBuilder()
    project = Project(id="PROJ001", name="Test", access=AccessLevel.ncats)
    exp_data = ParsedPounceData(
        experiments=[ParsedExperiment(experiment_name="Experiment 1")],
    )
    exp_parser = MagicMock(spec=ExcelsheetParser)
    exp_parser.file_path = "experiment.xlsx"
    exp_parser.sheet_dfs = {}

    with patch.object(PounceNodeBuilder, "_compute_experiment_counts", return_value=(0, 0)):
        first_batch = next(builder._experiment_nodes(project, exp_data, exp_parser))

    experiment = next(node for node in first_batch if isinstance(node, Experiment))

    assert experiment.access == AccessLevel.ncats


def test_parse_data_matrix_inherits_experiment_access_for_dataset():
    builder = PounceNodeBuilder()
    experiment = Experiment(id="EXP001", name="Experiment 1", access=AccessLevel.private)
    parser = MagicMock(spec=ExcelsheetParser)
    parser.sheet_dfs = {"ProteinData": pd.DataFrame()}
    parser.safe_get_string.return_value = None

    first_batch = next(builder._parse_data_matrix(
        experiment,
        meta_sheet="ProteinDataMeta",
        data_sheet="ProteinData",
        analyte_id_col="protein_id",
        data_type="protein data",
        parser=parser,
    ))
    dataset = next(node for node in first_batch if isinstance(node, Dataset))

    assert dataset.access == AccessLevel.private


def test_project_node_carries_workbook_artifact():
    parsed_project = ParsedProject(project_name="Test Project")

    project = PounceNodeBuilder._project_node(parsed_project, workbook_local_path="/tmp/test_project.xlsx")

    assert isinstance(project.workbook, WorkbookArtifact)
    assert project.workbook.original_filename == "test_project.xlsx"
    assert project.workbook._local_path == "/tmp/test_project.xlsx"


def test_build_analyte_edges_includes_proteins():
    builder = PounceNodeBuilder()
    parent = MagicMock()
    df = pd.DataFrame(index=["Q09666"])
    measured = MeasuredProtein(id="UniProtKB:Q09666", name="AHNAK")
    builder._ctx.protein_by_raw_id["Q09666"] = measured

    edges = builder._build_analyte_edges(parent, df, MagicMock, MagicMock, DatasetProteinEdge)

    assert len(edges) == 1
    assert isinstance(edges[0], DatasetProteinEdge)
    assert edges[0].end_node.id == "UniProtKB:Q09666"
