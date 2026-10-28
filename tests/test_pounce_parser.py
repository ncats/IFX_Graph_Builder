"""Tests for PounceParser."""
import pytest
import pandas as pd
from datetime import date
from unittest.mock import MagicMock

from src.input_adapters.pounce_sheets.validator_loader import load_validators
from src.input_adapters.excel_sheet_adapter import ExcelsheetParser
from src.input_adapters.pounce_sheets.pounce_parser import PounceParser
from src.input_adapters.pounce_sheets.parsed_classes import ParsedPerson
from src.models.pounce.project import ProjectType


# ---------------------------------------------------------------------------
# Helpers
# ---------------------------------------------------------------------------

def make_mock_parser(meta_data: dict = None, sheet_dfs: dict = None, param_maps: dict = None):
    """Build a mock ExcelsheetParser.

    Args:
        meta_data: {sheet_name: {ncatsdpi_key: value}} used by safe_get_string /
                   get_one_string_list / get_one_date.
        sheet_dfs: {sheet_name: DataFrame} returned by parser.sheet_dfs.
        param_maps: {sheet_name: {key: col}} returned by get_parameter_map.
    """
    meta_data = meta_data or {}
    sheet_dfs = sheet_dfs or {}
    param_maps = param_maps or {}

    mock = MagicMock(spec=ExcelsheetParser)
    mock.file_path = "test_file.xlsx"
    mock.sheet_dfs = sheet_dfs

    def safe_get_string(sheet, key):
        val = meta_data.get(sheet, {}).get(key)
        return str(val).strip() if isinstance(val, str) else (str(val) if val is not None else None)

    def get_one_string_list(sheet, key):
        val = meta_data.get(sheet, {}).get(key)
        if val is None:
            return []
        return val if isinstance(val, list) else [str(val)]

    def get_one_date(sheet, key):
        val = meta_data.get(sheet, {}).get(key)
        if val is None:
            raise KeyError(key)
        return val  # caller does .date() on it

    def get_parameter_map(sheet):
        return param_maps.get(sheet, {})

    mock.safe_get_string = MagicMock(side_effect=safe_get_string)
    mock.get_one_string_list = MagicMock(side_effect=get_one_string_list)
    mock.get_one_date = MagicMock(side_effect=get_one_date)
    mock.get_parameter_map = MagicMock(side_effect=get_parameter_map)

    return mock


# ---------------------------------------------------------------------------
# _parse_meta_field
# ---------------------------------------------------------------------------

class TestParseMetaField:

    def test_string_returns_value(self):
        parser = make_mock_parser({"S": {"k": "hello"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "string") == "hello"

    def test_string_returns_none_for_missing_key(self):
        parser = make_mock_parser()
        assert PounceParser._parse_meta_field(parser, "S", "missing", "string") is None

    def test_string_list_returns_list(self):
        parser = make_mock_parser({"S": {"k": ["a", "b"]}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "string_list") == ["a", "b"]

    def test_bool_yes_is_true(self):
        parser = make_mock_parser({"S": {"k": "yes"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is True

    def test_bool_true_is_true(self):
        parser = make_mock_parser({"S": {"k": "true"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is True

    def test_bool_one_is_true(self):
        parser = make_mock_parser({"S": {"k": "1"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is True

    def test_bool_no_is_false(self):
        parser = make_mock_parser({"S": {"k": "no"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is False

    def test_bool_false_string_is_false(self):
        parser = make_mock_parser({"S": {"k": "false"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is False

    def test_bool_missing_is_none(self):
        parser = make_mock_parser()
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is None

    def test_bool_is_case_insensitive(self):
        parser = make_mock_parser({"S": {"k": "YES"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "bool") is True

    def test_int_parses_integer_string(self):
        parser = make_mock_parser({"S": {"k": "42"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "int") == 42

    def test_int_truncates_float_string(self):
        parser = make_mock_parser({"S": {"k": "3.9"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "int") == 3

    def test_int_returns_none_for_invalid(self):
        parser = make_mock_parser({"S": {"k": "not_a_number"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "int") is None

    def test_int_returns_none_for_missing(self):
        parser = make_mock_parser()
        assert PounceParser._parse_meta_field(parser, "S", "k", "int") is None

    def test_float_parses_value(self):
        parser = make_mock_parser({"S": {"k": "3.14"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "float") == pytest.approx(3.14)

    def test_float_returns_none_for_invalid(self):
        parser = make_mock_parser({"S": {"k": "not_a_number"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "float") is None

    def test_float_returns_none_for_missing(self):
        parser = make_mock_parser()
        assert PounceParser._parse_meta_field(parser, "S", "k", "float") is None

    def test_date_returns_date_object(self):
        mock_dt = MagicMock()
        mock_dt.date.return_value = date(2024, 3, 15)
        parser = make_mock_parser({"S": {"k": mock_dt}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "date") == date(2024, 3, 15)

    def test_date_returns_none_when_key_missing(self):
        parser = make_mock_parser()
        assert PounceParser._parse_meta_field(parser, "S", "missing", "date") is None

    def test_category_returns_string(self):
        parser = make_mock_parser({"S": {"k": "Homo sapiens"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "category") == "Homo sapiens"

    def test_unknown_parse_type_falls_back_to_safe_get_string(self):
        parser = make_mock_parser({"S": {"k": "value"}})
        assert PounceParser._parse_meta_field(parser, "S", "k", "custom") == "value"


# ---------------------------------------------------------------------------
# _get_column_name
# ---------------------------------------------------------------------------

class TestGetColumnName:

    def test_returns_mapped_column(self):
        assert PounceParser._get_column_name({"biosample_id": "SampleID"}, "biosample_id") == "SampleID"

    def test_returns_none_for_missing_key(self):
        assert PounceParser._get_column_name({}, "missing") is None

    def test_returns_none_for_empty_string(self):
        assert PounceParser._get_column_name({"k": ""}, "k") is None

    def test_returns_none_for_na(self):
        assert PounceParser._get_column_name({"k": "NA"}, "k") is None

    def test_returns_none_for_n_slash_a(self):
        assert PounceParser._get_column_name({"k": "N/A"}, "k") is None

    def test_returns_none_for_none(self):
        assert PounceParser._get_column_name({"k": None}, "k") is None

    def test_indexed_key_formats_with_index(self):
        result = PounceParser._get_column_name({"exposure1_names": "DrugName"}, "exposure{}_names", index=1)
        assert result == "DrugName"

    def test_indexed_key_returns_none_for_wrong_index(self):
        result = PounceParser._get_column_name({"exposure1_names": "DrugName"}, "exposure{}_names", index=2)
        assert result is None


# ---------------------------------------------------------------------------
# _get_row_value
# ---------------------------------------------------------------------------

class TestGetRowValue:

    def test_returns_string_value(self):
        assert PounceParser._get_row_value({"Col": "S1"}, "Col") == "S1"

    def test_strips_whitespace(self):
        assert PounceParser._get_row_value({"Col": "  S1  "}, "Col") == "S1"

    def test_returns_numeric_value_unchanged(self):
        assert PounceParser._get_row_value({"Col": 42}, "Col") == "42"

    def test_returns_none_for_none_column_name(self):
        assert PounceParser._get_row_value({"Col": "S1"}, None) is None

    def test_returns_none_for_empty_column_name(self):
        assert PounceParser._get_row_value({"Col": "S1"}, "") is None

    def test_returns_none_for_none_cell_value(self):
        assert PounceParser._get_row_value({"Col": None}, "Col") is None

    def test_returns_none_for_empty_string_value(self):
        assert PounceParser._get_row_value({"Col": ""}, "Col") is None

    def test_returns_none_for_na_value(self):
        assert PounceParser._get_row_value({"Col": "NA"}, "Col") is None

    def test_returns_none_for_n_slash_a_value(self):
        assert PounceParser._get_row_value({"Col": "N/A"}, "Col") is None

    def test_list_splits_by_pipe(self):
        result = PounceParser._get_row_value({"Col": "Cancer|Diabetes"}, "Col", is_list=True)
        assert result == ["Cancer", "Diabetes"]

    def test_list_trims_whitespace(self):
        result = PounceParser._get_row_value({"Col": " Cancer | Diabetes "}, "Col", is_list=True)
        assert result == ["Cancer", "Diabetes"]

    def test_list_single_value(self):
        result = PounceParser._get_row_value({"Col": "Cancer"}, "Col", is_list=True)
        assert result == ["Cancer"]

    def test_list_returns_none_for_na(self):
        assert PounceParser._get_row_value({"Col": "NA"}, "Col", is_list=True) is None


# ---------------------------------------------------------------------------
# _build_persons
# ---------------------------------------------------------------------------

class TestBuildPersons:

    def test_returns_empty_for_no_names(self):
        assert PounceParser._build_persons([], [], "Owner") == []

    def test_builds_persons_with_matching_emails(self):
        result = PounceParser._build_persons(["Alice", "Bob"], ["a@x.com", "b@x.com"], "Owner")
        assert len(result) == 2
        assert result[0] == ParsedPerson(name="Alice", email="a@x.com", roles=["Owner"])
        assert result[1] == ParsedPerson(name="Bob", email="b@x.com", roles=["Owner"])

    def test_drops_emails_when_counts_differ(self):
        result = PounceParser._build_persons(["Alice", "Bob"], ["a@x.com"], "Owner")
        assert all(p.email is None for p in result)

    def test_no_emails_list(self):
        result = PounceParser._build_persons(["Alice"], [], "Collaborator")
        assert result == [ParsedPerson(name="Alice", email=None, roles=["Collaborator"])]

    def test_assigns_correct_role(self):
        result = PounceParser._build_persons(["Alice"], ["a@x.com"], "Collaborator")
        assert result[0].roles == ["Collaborator"]


def test_parse_experiment_parses_protein_metadata_and_protein_data_meta():
    parser = make_mock_parser(
        meta_data={
            "ExperimentMeta": {"experiment_name": "Proteomics Experiment"},
            "ProteinDataMeta": {
                "pre_processing_description": "search",
                "peri_processing_description": "normalize",
                "proteindata_tag": "LFQ_Intensity",
            },
        },
        sheet_dfs={
            "ExperimentMeta": pd.DataFrame(),
            "RunBioSampleMap": pd.DataFrame(columns=["NCATSDPI_Variable_Name"]),
            "RunBioSampleMeta": pd.DataFrame(columns=["run_biosample_id", "biosample_id"]),
            "ProteinMap": pd.DataFrame(columns=["NCATSDPI_Variable_Name"]),
            "ProteinMeta": pd.DataFrame([
                {"Uniprot_Acc": "Q09666", "Protein_Description": "AHNAK", "Gene_Symbol": "AHNAK"}
            ]),
            "ProteinDataMeta": pd.DataFrame(columns=["NCATSDPI_Variable_Name"]),
            "ProteinData": pd.DataFrame(columns=["Uniprot_Acc", "S1"]),
        },
        param_maps={
            "RunBioSampleMap": {
                "run_biosample_id": "run_biosample_id",
                "biosample_id": "biosample_id",
            },
            "ProteinMap": {
                "protein_id": "Uniprot_Acc",
                "protein_name": "Protein_Description",
                "gene_name": "Gene_Symbol",
            },
        },
    )

    data, issues = PounceParser().parse_experiment(parser)

    assert not [issue for issue in issues if issue.severity == "error"]
    assert len(data.proteins) == 1
    assert data.proteins[0].protein_id == "Q09666"
    assert data.proteins[0].protein_name == "AHNAK"
    assert data.proteins[0].gene_name == "AHNAK"
    assert len(data.protein_data_meta) == 1
    assert data.protein_data_meta[0].proteindata_tag == "LFQ_Intensity"


# ---------------------------------------------------------------------------
# _detect_exposure_indices
# ---------------------------------------------------------------------------

class TestDetectExposureIndices:

    def test_no_exposures_returns_empty(self):
        assert PounceParser._detect_exposure_indices({}) == []

    def test_single_exposure_by_names_key(self):
        assert PounceParser._detect_exposure_indices({"exposure1_names": "DrugName"}) == [1]

    def test_single_exposure_by_type_key(self):
        assert PounceParser._detect_exposure_indices({"exposure1_type": "Chemical"}) == [1]

    def test_multiple_contiguous_exposures(self):
        param_map = {
            "exposure1_names": "Drug1",
            "exposure2_names": "Drug2",
            "exposure3_type": "Chemical",
        }
        assert PounceParser._detect_exposure_indices(param_map) == [1, 2, 3]

    def test_stops_at_gap(self):
        param_map = {
            "exposure1_names": "Drug1",
            # no exposure2
            "exposure3_names": "Drug3",
        }
        assert PounceParser._detect_exposure_indices(param_map) == [1]

    def test_all_na_values_count_as_absent(self):
        param_map = {
            "exposure1_names": "NA",
            "exposure1_type": "NA",
            "exposure1_category": "NA",
        }
        assert PounceParser._detect_exposure_indices(param_map) == []


# ---------------------------------------------------------------------------
# parse_project
# ---------------------------------------------------------------------------

def _make_project_parser(meta: dict = None, biosample_map: dict = None, biosample_rows: list = None):
    """Convenience builder for project workbook mocks."""
    meta = meta or {}
    sheet_dfs = {}
    param_maps = {}

    if biosample_map is not None and biosample_rows is not None:
        sheet_dfs["BioSampleMap"] = pd.DataFrame()
        sheet_dfs["BioSampleMeta"] = pd.DataFrame(biosample_rows)
        param_maps["BioSampleMap"] = biosample_map

    return make_mock_parser(
        meta_data={"ProjectMeta": meta},
        sheet_dfs=sheet_dfs,
        param_maps=param_maps,
    )


BASIC_BIOSAMPLE_MAP = {
    "biosample_id": "SampleID",
    "biosample_type": "SampleType",
    "biospecimen_id": "SpecimenID",
    "biospecimen_type": "SpecimenType",
    "biospecimen_description": "Desc",
    "organism_names": "Organism",
    "disease_names": "Disease",
}


class TestParseProject:

    def test_parses_project_name(self):
        parser = _make_project_parser({"project_name": "My Project"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.project_name == "My Project"

    def test_parses_description(self):
        parser = _make_project_parser({"description": "A great project"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.description == "A great project"

    def test_parses_privacy_type(self):
        parser = _make_project_parser({"privacy_type": "public"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.privacy_type == "public"

    def test_parses_keywords_list(self):
        parser = _make_project_parser({"keywords": ["cancer", "genomics"]})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.keywords == ["cancer", "genomics"]

    def test_parses_lab_groups(self):
        parser = _make_project_parser({"lab_groups": ["Lab A", "Lab B"]})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.lab_groups == ["Lab A", "Lab B"]

    def test_parses_project_type(self):
        parser = _make_project_parser({"project_type": ["Transcriptomics"]})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.project_type == ["Transcriptomics"]

    def test_project_type_validator_allows_new_exposure_characterization_value(self):
        parser = _make_project_parser({"project_type": ["Drug or other Exposure Characterization"]})
        data, _ = PounceParser().parse_project(parser)
        validators = load_validators("src/use_cases/pounce/pounce_validators.yaml")

        project_type_validator = next(
            validator
            for validator in validators
            if validator.sheet == "ProjectMeta" and validator.column == "project_type"
        )

        assert project_type_validator.validate(data) == []
        assert ProjectType.parse("Drug or other Exposure Characterization") is (
            ProjectType.DRUG_OR_OTHER_EXPOSURE_CHARACTERIZATION
        )

    def test_rd_tag_yes_is_true(self):
        parser = _make_project_parser({"RD_tag": "yes"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.rd_tag is True

    def test_rd_tag_true_is_true(self):
        parser = _make_project_parser({"RD_tag": "true"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.rd_tag is True

    def test_rd_tag_no_is_false(self):
        parser = _make_project_parser({"RD_tag": "no"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.rd_tag is False

    def test_rd_tag_missing_is_none(self):
        parser = _make_project_parser({})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.rd_tag is None

    def test_biosample_preparation(self):
        parser = _make_project_parser({"biosample_preparation": "Standard prep"})
        data, _ = PounceParser().parse_project(parser)
        assert data.project.biosample_preparation == "Standard prep"

    def test_builds_owners_with_emails(self):
        parser = _make_project_parser({
            "owner_name": ["Alice", "Bob"],
            "owner_email": ["alice@x.com", "bob@x.com"],
        })
        data, _ = PounceParser().parse_project(parser)
        assert len(data.people) == 2
        assert data.people[0] == ParsedPerson(name="Alice", email="alice@x.com", roles=["Owner", "Contact"])
        assert data.people[1] == ParsedPerson(name="Bob", email="bob@x.com", roles=["Owner"])

    def test_builds_collaborators(self):
        parser = _make_project_parser({
            "collaborator_name": ["Carol"],
            "collaborator_email": ["carol@x.com"],
        })
        data, _ = PounceParser().parse_project(parser)
        collabs = [p for p in data.people if p.roles == ["Collaborator"]]
        assert len(collabs) == 1
        assert collabs[0].name == "Carol"
        assert collabs[0].email == "carol@x.com"

    def test_owners_and_collaborators_combined(self):
        parser = _make_project_parser({
            "owner_name": ["Alice"],
            "owner_email": ["alice@x.com"],
            "collaborator_name": ["Bob"],
            "collaborator_email": ["bob@x.com"],
        })
        data, _ = PounceParser().parse_project(parser)
        assert len(data.people) == 2
        roles = {tuple(p.roles) for p in data.people}
        assert roles == {("Owner", "Contact"), ("Collaborator",)}

    def test_no_people_when_names_absent(self):
        parser = _make_project_parser({})
        data, _ = PounceParser().parse_project(parser)
        assert data.people == []

    def test_parses_two_biosamples(self):
        rows = [
            {"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
             "SpecimenType": "human", "Desc": "d", "Organism": "Homo sapiens", "Disease": ""},
            {"SampleID": "S2", "SampleType": "tissue", "SpecimenID": "SP1",
             "SpecimenType": "human", "Desc": "d", "Organism": "Homo sapiens", "Disease": ""},
        ]
        parser = _make_project_parser({}, BASIC_BIOSAMPLE_MAP, rows)
        data, _ = PounceParser().parse_project(parser)
        assert len(data.biosamples) == 2
        assert data.biosamples[0].biosample_id == "S1"
        assert data.biosamples[1].biosample_id == "S2"

    def test_parses_biospecimen_fields(self):
        rows = [{"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
                 "SpecimenType": "human", "Desc": "tissue desc", "Organism": "Homo sapiens", "Disease": ""}]
        parser = _make_project_parser({}, BASIC_BIOSAMPLE_MAP, rows)
        data, _ = PounceParser().parse_project(parser)
        spec = data.biospecimens[0]
        assert spec.biospecimen_id == "SP1"
        assert spec.biospecimen_type == "human"
        assert spec.organism_names == "Homo sapiens"

    def test_disease_names_parsed_as_list(self):
        rows = [{"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
                 "SpecimenType": "human", "Desc": "d", "Organism": "H.s.",
                 "Disease": "Cancer|Diabetes"}]
        parser = _make_project_parser({}, BASIC_BIOSAMPLE_MAP, rows)
        data, _ = PounceParser().parse_project(parser)
        assert data.biospecimens[0].disease_names == ["Cancer", "Diabetes"]

    def test_no_biosamples_when_sheets_absent(self):
        parser = _make_project_parser({"project_id": "P1"})
        data, _ = PounceParser().parse_project(parser)
        assert data.biosamples == []
        assert data.biospecimens == []
        assert data.exposures == []

    def test_parses_single_exposure_slot(self):
        biosample_map = {**BASIC_BIOSAMPLE_MAP,
                         "exposure1_names": "DrugName",
                         "exposure1_type": "DrugType",
                         "exposure1_category": "DrugCat"}
        rows = [{"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
                 "SpecimenType": "human", "Desc": "d", "Organism": "H.s.", "Disease": "",
                 "DrugName": "Aspirin", "DrugType": "NSAID", "DrugCat": "drug"}]
        parser = _make_project_parser({}, biosample_map, rows)
        data, _ = PounceParser().parse_project(parser)
        assert len(data.exposures) == 1
        assert data.exposures[0].names == ["Aspirin"]
        assert data.exposures[0].type == "NSAID"
        assert data.exposures[0].exposure_index == 1

    def test_parses_multiple_exposure_slots_per_row(self):
        biosample_map = {**BASIC_BIOSAMPLE_MAP,
                         "exposure1_names": "Drug1Name",
                         "exposure1_type": "Drug1Type",
                         "exposure1_category": "Drug1Cat",
                         "exposure2_names": "Drug2Name",
                         "exposure2_type": "Drug2Type",
                         "exposure2_category": "Drug2Cat"}
        rows = [{"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
                 "SpecimenType": "human", "Desc": "d", "Organism": "H.s.", "Disease": "",
                 "Drug1Name": "Aspirin", "Drug1Type": "NSAID", "Drug1Cat": "drug",
                 "Drug2Name": "Caffeine", "Drug2Type": "Stimulant", "Drug2Cat": "drug"}]
        parser = _make_project_parser({}, biosample_map, rows)
        data, _ = PounceParser().parse_project(parser)
        assert len(data.exposures) == 2
        assert data.exposures[0].exposure_index == 1
        assert data.exposures[1].exposure_index == 2

    def test_exposure_count_multiplies_by_row_count(self):
        biosample_map = {**BASIC_BIOSAMPLE_MAP,
                         "exposure1_names": "DrugName",
                         "exposure1_type": "DrugType",
                         "exposure1_category": "DrugCat"}
        rows = [
            {"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
             "SpecimenType": "human", "Desc": "d", "Organism": "H.s.", "Disease": "",
             "DrugName": "Aspirin", "DrugType": "NSAID", "DrugCat": "drug"},
            {"SampleID": "S2", "SampleType": "tissue", "SpecimenID": "SP1",
             "SpecimenType": "human", "Desc": "d", "Organism": "H.s.", "Disease": "",
             "DrugName": "Ibuprofen", "DrugType": "NSAID", "DrugCat": "drug"},
        ]
        parser = _make_project_parser({}, biosample_map, rows)
        data, _ = PounceParser().parse_project(parser)
        # 2 rows × 1 exposure slot = 2 exposures
        assert len(data.exposures) == 2

    def test_biosample_param_map_is_stored(self):
        rows = [{"SampleID": "S1", "SampleType": "blood", "SpecimenID": "SP1",
                 "SpecimenType": "human", "Desc": "d", "Organism": "H.s.", "Disease": ""}]
        parser = _make_project_parser({}, BASIC_BIOSAMPLE_MAP, rows)
        data, _ = PounceParser().parse_project(parser)
        assert data.param_maps.get("BioSampleMap") == BASIC_BIOSAMPLE_MAP


# ---------------------------------------------------------------------------
# parse_experiment
# ---------------------------------------------------------------------------

def _make_experiment_parser(meta: dict = None, run_map: dict = None, run_rows: list = None):
    meta = meta or {}
    sheet_dfs = {}
    param_maps = {}

    if run_map is not None and run_rows is not None:
        sheet_dfs["RunBioSampleMap"] = pd.DataFrame()
        sheet_dfs["RunBioSampleMeta"] = pd.DataFrame(run_rows)
        param_maps["RunBioSampleMap"] = run_map

    return make_mock_parser(
        meta_data={"ExperimentMeta": meta},
        sheet_dfs=sheet_dfs,
        param_maps=param_maps,
    )


BASIC_RUN_MAP = {
    "run_biosample_id": "RunID",
    "biosample_id": "SampleID",
    "biological_replicate_number": "BioRep",
    "technical_replicate_number": "TechRep",
    "biosample_run_order": "RunOrder",
    "batch": "Batch",
}


class TestParseExperiment:

    def test_parses_experiment_name(self):
        parser = _make_experiment_parser({"experiment_name": "RNA-seq run 1"})
        data, _ = PounceParser().parse_experiment(parser)
        assert data.experiments[0].experiment_name == "RNA-seq run 1"

    def test_parses_experiment_type(self):
        parser = _make_experiment_parser({"experiment_type": "Transcriptomics"})
        data, _ = PounceParser().parse_experiment(parser)
        assert data.experiments[0].experiment_type == "Transcriptomics"

    def test_parses_platform_fields(self):
        parser = _make_experiment_parser({
            "platform_type": "RNA-seq",
            "platform_name": "Illumina NovaSeq",
            "platform_provider": "Illumina",
        })
        data, _ = PounceParser().parse_experiment(parser)
        exp = data.experiments[0]
        assert exp.platform_type == "RNA-seq"
        assert exp.platform_name == "Illumina NovaSeq"
        assert exp.platform_provider == "Illumina"

    def test_parses_run_biosamples(self):
        rows = [
            {"RunID": "R1", "SampleID": "S1", "BioRep": "1", "TechRep": "1", "RunOrder": "1"},
            {"RunID": "R2", "SampleID": "S2", "BioRep": "1", "TechRep": "1", "RunOrder": "2"},
        ]
        parser = _make_experiment_parser({}, BASIC_RUN_MAP, rows)
        data, _ = PounceParser().parse_experiment(parser)
        assert len(data.run_biosamples) == 2
        assert data.run_biosamples[0].run_biosample_id == "R1"
        assert data.run_biosamples[1].run_biosample_id == "R2"

    def test_run_biosample_links_to_biosample_id(self):
        rows = [{"RunID": "R1", "SampleID": "S1", "BioRep": "1", "TechRep": "1", "RunOrder": "1"}]
        parser = _make_experiment_parser({}, BASIC_RUN_MAP, rows)
        data, _ = PounceParser().parse_experiment(parser)
        assert data.run_biosamples[0].biosample_id == "S1"

    def test_parses_run_biosample_batch(self):
        rows = [{"RunID": "R1", "SampleID": "S1", "BioRep": "1", "TechRep": "1", "RunOrder": "1", "Batch": "2"}]
        parser = _make_experiment_parser({}, BASIC_RUN_MAP, rows)
        data, _ = PounceParser().parse_experiment(parser)
        assert data.run_biosamples[0].batch == "2"

    def test_no_run_biosamples_when_sheets_absent(self):
        parser = _make_experiment_parser({"experiment_name": "EXP001"})
        data, _ = PounceParser().parse_experiment(parser)
        assert data.run_biosamples == []

    def test_always_produces_one_experiment(self):
        parser = _make_experiment_parser({})
        data, _ = PounceParser().parse_experiment(parser)
        assert len(data.experiments) == 1


# ---------------------------------------------------------------------------
# parse_stats_results
# ---------------------------------------------------------------------------

def _make_stats_parser(meta: dict = None):
    meta = meta or {}
    return make_mock_parser(
        meta_data={"StatsResultsMeta": meta},
        sheet_dfs={"StatsResultsMeta": pd.DataFrame()},
    )


class TestParseStatsResults:

    def test_parses_stats_results_name(self):
        parser = _make_stats_parser({"statsresults_name": "My Stats"})
        data, _ = PounceParser().parse_stats_results(parser)
        assert len(data.stats_results) == 1
        assert data.stats_results[0].statsresults_name == "My Stats"

    def test_parses_stats_description(self):
        parser = _make_stats_parser({"stats_description": "DESeq2 analysis"})
        data, _ = PounceParser().parse_stats_results(parser)
        assert data.stats_results[0].stats_description == "DESeq2 analysis"

    def test_parses_lead_informatician(self):
        parser = _make_stats_parser({"lead_informatician": "Dr Smith"})
        data, _ = PounceParser().parse_stats_results(parser)
        assert data.stats_results[0].lead_informatician == "Dr Smith"

    def test_no_stats_results_when_sheet_absent(self):
        mock = make_mock_parser(sheet_dfs={})
        data, _ = PounceParser().parse_stats_results(mock)
        assert data.stats_results == []


# ---------------------------------------------------------------------------
# _parse_demographics_row
# ---------------------------------------------------------------------------

class TestParseDemographicsRow:

    def test_parses_basic_fields(self):
        row = {"Age": "45", "Race": "White", "Ethnicity": "Non-Hispanic", "Sex": "F"}
        param_map = {"age": "Age", "race": "Race", "ethnicity": "Ethnicity", "sex": "Sex"}
        result = PounceParser._parse_demographics_row(row, param_map)
        assert result.age == "45"
        assert result.race == "White"
        assert result.ethnicity == "Non-Hispanic"
        assert result.sex == "F"

    def test_parses_indexed_demographic_categories(self):
        row = {"DemoCat1": "GroupA", "DemoCat2": "GroupB"}
        param_map = {
            "demographic1_category": "DemoCat1",
            "demographic2_category": "DemoCat2",
        }
        result = PounceParser._parse_demographics_row(row, param_map)
        assert "DemoCat1:GroupA" in result.demographic_categories
        assert "DemoCat2:GroupB" in result.demographic_categories

    def test_parses_indexed_phenotype_categories(self):
        row = {"Pheno1": "Hypertension"}
        param_map = {"phenotype1_category": "Pheno1"}
        result = PounceParser._parse_demographics_row(row, param_map)
        assert "Pheno1:Hypertension" in result.phenotype_categories

    def test_skips_null_category_values(self):
        row = {"DemoCat1": "GroupA", "DemoCat2": None}
        param_map = {
            "demographic1_category": "DemoCat1",
            "demographic2_category": "DemoCat2",
        }
        result = PounceParser._parse_demographics_row(row, param_map)
        assert len(result.demographic_categories) == 1
        assert "DemoCat1:GroupA" in result.demographic_categories

    def test_empty_demographics_when_no_columns_mapped(self):
        result = PounceParser._parse_demographics_row({}, {})
        assert result.age is None
        assert result.race is None
        assert result.demographic_categories == []
        assert result.phenotype_categories == []

    def test_demographic_categories_stop_when_column_absent(self):
        # Only demographic1 mapped, no demographic2 → stops after 1
        row = {"DemoCat1": "GroupA"}
        param_map = {"demographic1_category": "DemoCat1"}
        result = PounceParser._parse_demographics_row(row, param_map)
        assert len(result.demographic_categories) == 1


# ---------------------------------------------------------------------------
# parse_all (integration)
# ---------------------------------------------------------------------------

# ---------------------------------------------------------------------------
# _check_mapped_columns
# ---------------------------------------------------------------------------

class TestCheckMappedColumns:

    def _make_parser(self, meta_columns, param_map):
        mock = MagicMock(spec=ExcelsheetParser)
        mock.file_path = "test.xlsx"
        mock.sheet_dfs = {
            "BioSampleMeta": pd.DataFrame(columns=meta_columns)
        }
        mock.get_parameter_map = MagicMock(return_value=param_map)
        return mock

    def test_all_columns_present_returns_no_errors(self):
        parser = self._make_parser(["SampleID", "SampleType"], {"biosample_id": "SampleID", "biosample_type": "SampleType"})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_id": "SampleID", "biosample_type": "SampleType"})
        assert issues == []

    def test_missing_column_returns_error(self):
        parser = self._make_parser(["SampleID"], {"biosample_id": "SampleID", "biosample_type": "WrongColumn"})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_id": "SampleID", "biosample_type": "WrongColumn"})
        assert len(issues) == 1
        assert "WrongColumn" in issues[0].message
        assert issues[0].severity == "error"

    def test_error_references_ncatsdpi_key_and_column(self):
        parser = self._make_parser([], {"biosample_id": "Typo"})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_id": "Typo"})
        assert issues[0].field == "biosample_id"
        assert issues[0].column == "biosample_id"
        assert "Typo" in issues[0].message

    def test_na_value_is_skipped(self):
        parser = self._make_parser([], {"biosample_id": "NA"})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_id": "NA"})
        assert issues == []

    def test_empty_value_is_skipped(self):
        parser = self._make_parser([], {"biosample_id": ""})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_id": ""})
        assert issues == []

    def test_meta_sheet_absent_returns_no_errors(self):
        mock = MagicMock(spec=ExcelsheetParser)
        mock.file_path = "test.xlsx"
        mock.sheet_dfs = {}  # meta sheet not present
        issues = PounceParser._check_mapped_columns(mock, "BioSampleMap", "BioSampleMeta", {"biosample_id": "SampleID"})
        assert issues == []

    def test_multiple_missing_columns_all_reported(self):
        parser = self._make_parser([], {"biosample_id": "Col1", "biosample_type": "Col2"})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_id": "Col1", "biosample_type": "Col2"})
        assert len(issues) == 2

    def test_error_includes_both_sheet_names_in_message(self):
        parser = self._make_parser(["SampleID"], {"biosample_type": "BadCol"})
        issues = PounceParser._check_mapped_columns(parser, "BioSampleMap", "BioSampleMeta", {"biosample_type": "BadCol"})
        assert "BioSampleMap" in issues[0].message
        assert "BioSampleMeta" in issues[0].message


class TestParseAll:

    def test_parse_all_no_experiment_or_stats_files(self, tmp_path, monkeypatch):
        """parse_all with empty experiment/stats lists delegates to parse_project."""
        captured = {}

        def fake_parse_project(self, parser):
            from src.input_adapters.pounce_sheets.parsed_pounce_data import ParsedPounceData
            from src.input_adapters.pounce_sheets.parsed_classes import ParsedProject
            d = ParsedPounceData()
            d.project = ParsedProject(project_name="P1")
            captured["called"] = True
            return d, []

        monkeypatch.setattr(PounceParser, "parse_project", fake_parse_project)
        from types import SimpleNamespace
        monkeypatch.setattr(
            "src.input_adapters.pounce_sheets.pounce_parser.load_pounce_workbook",
            lambda file_path, role: SimpleNamespace(
                parser=MagicMock(spec=ExcelsheetParser, sheet_dfs={}),
                version=2,
                original_sheets={},
                original_sheet=lambda sheet: sheet,
            ),
        )

        result, _ = PounceParser().parse_all("fake_project.xlsx")
        assert captured.get("called") is True
        assert result.project.project_name == "P1"
        assert result.experiments == []
        assert result.stats_results == []
