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
from src.models.pounce.experiment import Experiment, ExperimentBiospecimenEdge, ExperimentMeasuredGeneEdge, PlatformType, RunBiosample
from src.models.protein import MeasuredProtein
from src.models.gene import MeasuredGene, MeasuredGeneEdge
from src.models.pounce.dataset import Dataset, DatasetGeneEdge, DatasetMetaboliteEdge, DatasetProteinEdge
from src.models.pounce.raw_data_files import RawDataFiles, ExperimentRawDataFilesEdge, RunBiosampleRawDataFilesEdge
from src.models.pounce.workbook_artifact import WorkbookArtifact
from src.constants import DataSourceName


def test_experiment_stats_pairing_keeps_missing_middle_slot():
    assert list(_iter_pairs(["exp1", "exp2", "exp3"], ["stats1", None, "stats3"])) == [
        ("exp1", "stats1"), ("exp2", None), ("exp3", "stats3"),
    ]
    with pytest.raises(ValueError, match="use None"):
        list(_iter_pairs(["exp1", "exp2"], ["stats2"]))


def test_gene_matrix_without_id_mapping_fails_with_actionable_error():
    parser = MagicMock(spec=ExcelsheetParser)
    parser.file_path = "sc.xlsx"
    parser.sheet_dfs = {"RawData": pd.DataFrame([{"Symbol": "TP53", "Sample1": 2}])}
    with pytest.raises(ValueError, match="sc.xlsx / RawData.*analyte ID column"):
        PounceNodeBuilder._prepare_data_frame(parser, "RawData", None)
    parser.sheet_dfs["RawData"] = pd.DataFrame()
    assert PounceNodeBuilder._prepare_data_frame(parser, "RawData", None) is None


def test_populated_gene_meta_with_na_mapping_fails_before_silent_skip():
    parser = MagicMock(spec=ExcelsheetParser)
    parser.file_path = "sc.xlsx"
    parser.sheet_dfs = {"GeneMeta": pd.DataFrame([{"Symbol": "TP53"}])}
    parser.get_parameter_map.return_value = {"gene_id": "NA"}
    with pytest.raises(ValueError, match="GeneMeta.*GeneMap gene_id mapping"):
        list(PounceNodeBuilder()._parse_genes(None, parser))


def test_mixed_single_cell_features_keep_source_ids_without_false_gene_links():
    parser = MagicMock(spec=ExcelsheetParser)
    parser.sheet_dfs = {"GeneMeta": pd.DataFrame([
        {"GeneID": "ENSG00000141510.17", "Symbol": "TP53"},
        {"GeneID": "ENSG00000141510.17_PAR_Y", "Symbol": "TP53_PAR_Y"},
        {"GeneID": " gene-S ", "Symbol": "CoV2_S"},
    ])}
    parser.get_parameter_map.return_value = {
        "gene_id": "GeneID", "hgnc_gene_symbol": "Symbol",
    }
    builder = PounceNodeBuilder()
    batches = list(builder._parse_genes(Experiment(id="E1", platform_type=PlatformType.SC_RNASEQ), parser))
    nodes = [item for item in batches[0] if isinstance(item, MeasuredGene)]
    edges = [item for item in batches[1] if isinstance(item, MeasuredGeneEdge)]
    assert {node.id for node in nodes} == {
        "ENSEMBL:ENSG00000141510.17",
        "E1:measured_gene:ENSG00000141510.17_PAR_Y",
        "E1:measured_gene:gene-S",
    }
    assert len(edges) == 1
    assert edges[0].end_node.id == "ENSEMBL:ENSG00000141510.17"
    experiment_edges = [item for item in batches[2] if isinstance(item, ExperimentMeasuredGeneEdge)]
    assert len(experiment_edges) == 3
    assert experiment_edges[2].details[0].source_feature_id == " gene-S "
    assert experiment_edges[2].details[0].source_symbol == "CoV2_S"
    from src.output_adapters.arango_output_adapter import ArangoOutputAdapter
    output = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    assert output.clean_dict(experiment_edges[2], True)["details"] == [
        {
            "source_feature_id": " gene-S ", "source_symbol": "CoV2_S",
            "source_biotype": None, "run_biosample_ids": [],
        }
    ]
    assert builder._ctx.gene_by_raw_id[" gene-S "] is builder._ctx.gene_by_raw_id["gene-S"]
    matrix_edges = builder._build_analyte_edges(
        Dataset(id="D1"), pd.DataFrame(index=[" gene-S "]),
        DatasetGeneEdge, DatasetMetaboliteEdge, DatasetProteinEdge,
    )
    assert len(matrix_edges) == 1
    assert matrix_edges[0].end_node.id == "E1:measured_gene:gene-S"


def test_sc_gene_membership_rejects_run_from_another_experiment():
    parser = MagicMock(spec=ExcelsheetParser)
    parser.file_path = "second_experiment.xlsx"
    parser.sheet_dfs = {"GeneMeta": pd.DataFrame([{
        "GeneID": "ENSG00000141510.17", "Runs": "R1",
    }])}
    parser.get_parameter_map.return_value = {"gene_id": "GeneID", "runbiosample": "Runs"}
    builder = PounceNodeBuilder()
    builder._ctx.run_biosample_by_original_id["R1"] = RunBiosample(id="first-run")
    with pytest.raises(ValueError, match="unknown runbiosample IDs.*R1"):
        list(builder._parse_genes(
            Experiment(id="E2", platform_type=PlatformType.SC_RNASEQ),
            parser, valid_run_ids={"R2"},
        ))


def test_all_null_matrix_id_column_has_actionable_error():
    parser = MagicMock(spec=ExcelsheetParser)
    parser.file_path = "sc.xlsx"
    parser.sheet_dfs = {"RawData": pd.DataFrame([{"GeneID": None, "Sample1": 2}])}
    with pytest.raises(ValueError, match="mapped analyte ID column 'GeneID' has no values"):
        PounceNodeBuilder._prepare_data_frame(parser, "RawData", "GeneID")


def test_direct_build_rejects_sc_stats_without_gene_id(monkeypatch):
    from types import SimpleNamespace

    project_data = ParsedPounceData(project=ParsedProject(project_name="P1"))
    experiment_data = ParsedPounceData(
        experiments=[ParsedExperiment(platform_type="sc_rnaseq")],
    )
    mock_parser = MagicMock()
    mock_parser.parse_project.return_value = (project_data, [])
    mock_parser.parse_experiment.return_value = (experiment_data, [])
    mock_parser.parse_stats_results.return_value = (ParsedPounceData(), [])
    stats_sheet_parser = SimpleNamespace(
        file_path="stats.xlsx",
        sheet_dfs={"EffectSize": pd.DataFrame([{"Symbol": "TP53"}])},
    )
    workbooks = {
        "project.xlsx": SimpleNamespace(parser=SimpleNamespace()),
        "sc.xlsx": SimpleNamespace(parser=SimpleNamespace()),
        "stats.xlsx": SimpleNamespace(parser=stats_sheet_parser),
    }
    monkeypatch.setattr(
        "src.input_adapters.pounce_sheets.pounce_node_builder.load_pounce_workbook",
        lambda path, role: workbooks[path],
    )
    monkeypatch.setattr(PounceNodeBuilder, "_project_node", lambda *args, **kwargs: object())
    monkeypatch.setattr(PounceNodeBuilder, "_person_nodes_and_edges", lambda *args: ([], []))
    monkeypatch.setattr(PounceNodeBuilder, "_biosample_nodes", lambda *args: iter(()))
    with pytest.raises(ValueError, match="stats.xlsx.*sc_rnaseq gene-level.*GeneMeta"):
        list(PounceNodeBuilder().build(
            "project.xlsx", ["sc.xlsx"], ["stats.xlsx"], mock_parser,
        ))


def test_direct_build_rejects_invalid_sc_file_manifest(monkeypatch):
    from types import SimpleNamespace
    from src.core.validator import ValidationError

    project_data = ParsedPounceData(project=ParsedProject(project_name="P1"))
    experiment_data = ParsedPounceData(
        experiments=[ParsedExperiment(platform_type="sc_rnaseq")],
    )
    parser = MagicMock()
    parser.parse_project.return_value = (project_data, [])
    parser.parse_experiment.return_value = (experiment_data, [ValidationError(
        severity="error", entity="raw_data_files", field="seurat_file",
        sheet="RawData", message="provide a complete trio or seurat_file",
    )])
    monkeypatch.setattr(
        "src.input_adapters.pounce_sheets.pounce_node_builder.load_pounce_workbook",
        lambda path, role: SimpleNamespace(parser=SimpleNamespace(file_path=path)),
    )
    monkeypatch.setattr(PounceNodeBuilder, "_project_node", lambda *args, **kwargs: object())
    monkeypatch.setattr(PounceNodeBuilder, "_person_nodes_and_edges", lambda *args: ([], []))
    monkeypatch.setattr(PounceNodeBuilder, "_biosample_nodes", lambda *args: iter(()))
    with pytest.raises(ValueError, match="sc.xlsx / RawData: provide a complete trio"):
        list(PounceNodeBuilder().build("project.xlsx", ["sc.xlsx"], [None], parser))


def test_direct_build_rejects_sc_stats_gene_missing_from_gene_meta(monkeypatch):
    from types import SimpleNamespace
    from src.core.validator import ValidationError
    from src.input_adapters.pounce_sheets.parsed_classes import ParsedGene

    parser = MagicMock()
    parser.parse_project.return_value = (
        ParsedPounceData(project=ParsedProject(project_name="P1")), [],
    )
    parser.parse_experiment.return_value = (
        ParsedPounceData(
            experiments=[ParsedExperiment(platform_type="sc_rnaseq")],
            genes=[ParsedGene(gene_id="ENSG1")],
            param_maps={"GeneMap": {"gene_id": "GeneID"}},
        ), [],
    )
    parser.parse_stats_results.return_value = (ParsedPounceData(), [ValidationError(
        severity="error", entity="stats_results", field="gene_id", sheet="EffectSize",
        message="gene_id 'ENSG2' in EffectSize does not exist in paired experiment's GeneMeta",
    )])
    workbooks = {
        "project.xlsx": SimpleNamespace(parser=SimpleNamespace(file_path="project.xlsx")),
        "sc.xlsx": SimpleNamespace(parser=SimpleNamespace(file_path="sc.xlsx")),
        "stats.xlsx": SimpleNamespace(parser=SimpleNamespace(file_path="stats.xlsx")),
    }
    monkeypatch.setattr(
        "src.input_adapters.pounce_sheets.pounce_node_builder.load_pounce_workbook",
        lambda path, role: workbooks[path],
    )
    monkeypatch.setattr(PounceNodeBuilder, "_project_node", lambda *args, **kwargs: object())
    monkeypatch.setattr(PounceNodeBuilder, "_person_nodes_and_edges", lambda *args: ([], []))
    monkeypatch.setattr(PounceNodeBuilder, "_biosample_nodes", lambda *args: iter(()))
    with pytest.raises(ValueError, match="stats.xlsx / EffectSize: gene_id 'ENSG2'"):
        list(PounceNodeBuilder().build("project.xlsx", ["sc.xlsx"], ["stats.xlsx"], parser))
    assert parser.parse_stats_results.call_args.kwargs["valid_gene_ids"] == {"ENSG1"}


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


def test_biospecimen_ids_with_repeated_whitespace_share_graph_identity():
    builder = PounceNodeBuilder()
    project = Project(id="PROJ001", name="Test")
    project_data = ParsedPounceData(
        biosamples=[
            ParsedBiosample(biosample_id="BS1", biosample_type="sample"),
            ParsedBiosample(biosample_id="BS2", biosample_type="sample"),
        ],
        biospecimens=[
            ParsedBiospecimen(biospecimen_id="CTBPRO2,  GFP-HUVECs", biospecimen_type="specimen"),
            ParsedBiospecimen(biospecimen_id="CTBPRO2, GFP-HUVECs", biospecimen_type="specimen"),
        ],
    )

    batches = list(builder._biosample_nodes(project, project_data))
    biospecimens = batches[0]
    biosamples = batches[1]
    assert len(biospecimens) == 1
    assert biospecimens[0].id == "PROJ001-CTBPRO2, GFP-HUVECs"
    assert biospecimens[0].original_id == "CTBPRO2,  GFP-HUVECs"
    assert len({sample.id for sample in biosamples}) == 2
    assert len(batches[-1]) == 2
    assert all(edge.end_node.id == biospecimens[0].id for edge in batches[-1])


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


def test_geo_metadata_reaches_project_experiment_and_run_biosample():
    builder = PounceNodeBuilder()
    parsed_project = ParsedProject(
        project_name="GEO study", source_record="https://example.org/GSE1",
        geo_title="Deposited title", geo_submission_date="Jun 11 2026",
        curation_notes="Mapped from GEO",
    )
    project = builder._project_node(parsed_project)
    assert project.source_record == "https://example.org/GSE1"
    assert project.geo_title == "Deposited title"
    assert project.curation_notes == "Mapped from GEO"

    exp_data = ParsedPounceData(
        experiments=[ParsedExperiment(
            experiment_name="GEO experiment", geo_submission_date="Jun 11 2026",
            curation_notes="Experiment curation",
        )],
        run_biosamples=[ParsedRunBiosample(
            run_biosample_id="R1", biosample_id="B1", runbiosample_file="GSM1",
            geo_sample_accession="GSM1", geo_sample_url="https://example.org/SAMN1",
            geo_sra_url="https://example.org/SRX1", replicate_label="rep 1",
            library_text="GRCh38",
        )],
    )
    exp_parser = MagicMock(spec=ExcelsheetParser)
    exp_parser.file_path = "experiment.xlsx"
    exp_parser.sheet_dfs = {}
    with patch.object(PounceNodeBuilder, "_compute_experiment_counts", return_value=(0, 0)):
        batches = list(builder._experiment_nodes(project, exp_data, exp_parser))
    nodes = [node for batch in batches for node in batch]
    experiment = next(node for node in nodes if isinstance(node, Experiment))
    run = next(node for node in nodes if isinstance(node, RunBiosample))
    assert experiment.geo_submission_date == "Jun 11 2026"
    assert experiment.curation_notes == "Experiment curation"
    assert run.geo_sample_accession == "GSM1"
    assert run.geo_sra_url == "https://example.org/SRX1"
    assert run.source_file_name == "GSM1"
    assert run.library_text == "GRCh38"


def test_geo_biosample_characteristics_and_raw_data_notes_are_preserved():
    builder = PounceNodeBuilder()
    project = Project(id="PROJ001", name="GEO study")
    project_data = ParsedPounceData(
        biosamples=[ParsedBiosample(
            biosample_id="B1", biosample_type="sample",
            geo_characteristics='{"tissue": "brain"}',
        )],
        biospecimens=[ParsedBiospecimen(
            biospecimen_id="S1", biospecimen_type="specimen",
        )],
    )
    biosample_batches = list(builder._biosample_nodes(project, project_data))
    biosample = next(
        node for batch in biosample_batches for node in batch
        if getattr(node, "original_id", None) == "B1"
    )
    assert biosample.geo_characteristics == '{"tissue": "brain"}'

    parser = MagicMock(spec=ExcelsheetParser)
    parser.sheet_dfs = {"RawData": pd.DataFrame()}
    parser.safe_get_string.side_effect = lambda _sheet, key: {
        "reported_data_processing": "Deposited process",
        "data_representation": "Count matrix",
        "source_feature_type_note": "Source feature note",
    }.get(key)
    batches = list(builder._parse_data_matrix(
        Experiment(id="E1"), "RawDataMeta", "RawData", parser=parser,
    ))
    dataset = next(node for batch in batches for node in batch if isinstance(node, Dataset))
    assert dataset.reported_data_processing == "Deposited process"
    assert dataset.data_representation == "Count matrix"
    assert dataset.source_feature_type_note == "Source feature note"


def test_sc_raw_data_emits_file_references_without_matrix_dataset():
    builder = PounceNodeBuilder()
    run = RunBiosample(id="P1-R1")
    builder._ctx.run_biosample_by_original_id["R1"] = run
    parser = MagicMock(spec=ExcelsheetParser)
    parser.sheet_dfs = {"RawData": pd.DataFrame([{
        "runbiosample_id": "R1", "seurat_file": "https://example.org/object.rds",
        "GEO_sample": "GSM1", "feature_count": 100,
    }])}
    batches = list(builder._parse_sc_raw_data_files(Experiment(id="E1"), parser))
    records = [record for batch in batches for record in batch]
    file_set = next(record for record in records if isinstance(record, RawDataFiles))
    assert file_set.id.startswith("E1:raw_data_files:")
    assert file_set.seurat_file == "https://example.org/object.rds"
    assert file_set.geo_sample == "GSM1"
    assert file_set.feature_count == 100
    assert not any(isinstance(record, Dataset) for record in records)
    assert any(isinstance(record, ExperimentRawDataFilesEdge) for record in records)
    assert any(isinstance(record, RunBiosampleRawDataFilesEdge) for record in records)


def test_sc_experiment_routes_raw_data_to_file_manifest():
    builder = PounceNodeBuilder()
    project = Project(id="P1", name="Single cell")
    exp_data = ParsedPounceData(
        experiments=[ParsedExperiment(experiment_name="SC", platform_type="sc_rnaseq")],
        run_biosamples=[ParsedRunBiosample(run_biosample_id="R1", biosample_id="B1")],
    )
    parser = MagicMock(spec=ExcelsheetParser)
    parser.file_path = "sc.xlsx"
    parser.sheet_dfs = {"RawData": pd.DataFrame([{
        "runbiosample_id": "R1", "seurat_file": "https://example.org/object.rds",
    }])}
    with patch.object(PounceNodeBuilder, "_compute_experiment_counts", return_value=(0, 1)):
        batches = list(builder._experiment_nodes(project, exp_data, parser))
    records = [record for batch in batches for record in batch]
    assert any(isinstance(record, RawDataFiles) for record in records)
    assert not any(isinstance(record, Dataset) for record in records)


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
