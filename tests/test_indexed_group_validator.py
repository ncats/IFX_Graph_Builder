"""Tests for map-sheet validators and the indexed_group YAML loader."""

import pytest
from src.core.validator import IndexedGroupValidator, RequiredMapKeyValidator
from src.input_adapters.pounce_sheets.validator_loader import load_validators
from src.input_adapters.pounce_sheets.parsed_pounce_data import ParsedPounceData


_REQUIRED = ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
_SHEET = "BioSampleMap"


def _make_data(param_map: dict) -> ParsedPounceData:
    data = ParsedPounceData()
    data.param_maps = {_SHEET: param_map}
    return data


def _make_validator(required=None):
    return IndexedGroupValidator(sheet=_SHEET, required_templates=required or _REQUIRED)


# ---------------------------------------------------------------------------
# IndexedGroupValidator unit tests
# ---------------------------------------------------------------------------

class TestIndexedGroupValidatorDirect:

    def test_no_exposure_keys_passes(self):
        data = _make_data({"biosample_id": "SampleID"})
        assert _make_validator().validate(data) == []

    def test_configured_exposure_category_is_optional_in_pounce_rules(self):
        validator = next(
            v for v in load_validators("src/use_cases/pounce/pounce_validators.yaml")
            if isinstance(v, IndexedGroupValidator) and v.sheet == "BioSampleMap"
        )
        assert validator.validate(_make_data({
            "exposure1_names": "DrugName", "exposure1_type": "DrugType",
        })) == []
        errors = validator.validate(_make_data({"exposure1_category": "DrugCat"}))
        assert {error.column for error in errors} == {
            "exposure1_names", "exposure1_type",
        }

    def test_complete_exposure1_passes(self):
        data = _make_data({
            "exposure1_names": "DrugName",
            "exposure1_type": "DrugType",
            "exposure1_category": "DrugCat",
        })
        assert _make_validator().validate(data) == []

    def test_optional_only_triggers_check(self):
        # exposure1_concentration is present but the three required ones are not
        data = _make_data({"exposure1_concentration": "ConcentrationCol"})
        errors = _make_validator().validate(data)
        assert len(errors) == 3  # all three required keys missing

    def test_missing_type_produces_error(self):
        data = _make_data({
            "exposure1_names": "DrugName",
            "exposure1_category": "DrugCat",
            # exposure1_type missing
        })
        errors = _make_validator().validate(data)
        assert len(errors) == 1
        assert "exposure1_type" in errors[0].message
        assert errors[0].column == "exposure1_type"

    def test_missing_names_and_category_produces_two_errors(self):
        data = _make_data({"exposure1_type": "DrugType"})
        errors = _make_validator().validate(data)
        assert len(errors) == 2
        missing_cols = {e.column for e in errors}
        assert missing_cols == {"exposure1_names", "exposure1_category"}

    def test_na_value_counts_as_missing(self):
        data = _make_data({
            "exposure1_names": "DrugName",
            "exposure1_type": "NA",       # NA = not configured
            "exposure1_category": "DrugCat",
        })
        errors = _make_validator().validate(data)
        assert len(errors) == 1
        assert errors[0].column == "exposure1_type"

    def test_empty_value_counts_as_missing(self):
        data = _make_data({
            "exposure1_names": "DrugName",
            "exposure1_type": "",
            "exposure1_category": "DrugCat",
        })
        errors = _make_validator().validate(data)
        assert len(errors) == 1

    def test_two_complete_exposure_groups_pass(self):
        data = _make_data({
            "exposure1_names": "Drug1Name",
            "exposure1_type": "Drug1Type",
            "exposure1_category": "Drug1Cat",
            "exposure2_names": "Drug2Name",
            "exposure2_type": "Drug2Type",
            "exposure2_category": "Drug2Cat",
        })
        assert _make_validator().validate(data) == []

    def test_second_group_incomplete_errors_on_correct_index(self):
        data = _make_data({
            "exposure1_names": "Drug1Name",
            "exposure1_type": "Drug1Type",
            "exposure1_category": "Drug1Cat",
            "exposure2_names": "Drug2Name",
            # exposure2_type and exposure2_category missing
        })
        errors = _make_validator().validate(data)
        assert len(errors) == 2
        assert all("2" in e.column for e in errors)

    def test_error_has_correct_entity_and_sheet(self):
        data = _make_data({"exposure1_names": "DrugName"})
        errors = _make_validator().validate(data)
        assert errors[0].entity == "param_maps"
        assert errors[0].sheet == _SHEET
        assert errors[0].severity == "error"

    def test_error_field_is_template_not_concrete_key(self):
        # field should be the template ("exposure{}_type"), not the concrete key
        data = _make_data({"exposure1_names": "DrugName"})
        errors = _make_validator().validate(data)
        fields = {e.field for e in errors}
        assert "exposure{}_type" in fields
        assert "exposure{}_category" in fields

    def test_empty_param_maps_passes(self):
        data = ParsedPounceData()  # param_maps is empty dict
        assert _make_validator().validate(data) == []

    def test_sheet_absent_from_param_maps_passes(self):
        data = _make_data({})
        data.param_maps = {}  # sheet not present at all
        assert _make_validator().validate(data) == []


# ---------------------------------------------------------------------------
# load_validators with indexed_group
# ---------------------------------------------------------------------------

class TestLoadValidatorsIndexedGroup:

    def test_loads_indexed_group_validator(self):
        config = {
            "BioSampleMap": {
                "indexed_group": [
                    ["exposure{}_names", "exposure{}_type", "exposure{}_category"],
                ]
            }
        }
        validators = load_validators(config)
        assert len(validators) == 1
        assert isinstance(validators[0], IndexedGroupValidator)

    def test_loaded_validator_sheet(self):
        config = {
            "BioSampleMap": {
                "indexed_group": [
                    ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
                ]
            }
        }
        v = load_validators(config)[0]
        assert v.sheet == "BioSampleMap"

    def test_loaded_validator_required_templates(self):
        templates = ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
        config = {"BioSampleMap": {"indexed_group": [templates]}}
        v = load_validators(config)[0]
        assert v._required_templates == templates

    def test_multiple_groups_produce_multiple_validators(self):
        config = {
            "BioSampleMap": {
                "indexed_group": [
                    ["exposure{}_names", "exposure{}_type", "exposure{}_category"],
                    ["condition{}_category"],
                ]
            }
        }
        validators = load_validators(config)
        assert len(validators) == 2
        assert all(isinstance(v, IndexedGroupValidator) for v in validators)

    def test_mismatched_prefixes_raises(self):
        config = {
            "BioSampleMap": {
                "indexed_group": [
                    ["exposure{}_names", "condition{}_category"],  # different prefixes
                ]
            }
        }
        with pytest.raises(ValueError, match="prefixes"):
            load_validators(config)

    def test_indexed_group_alongside_required(self):
        config = {
            "BioSampleMap": {
                "required": ["biosample_id"],
                "indexed_group": [
                    ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
                ],
            }
        }
        validators = load_validators(config)
        types = [type(v).__name__ for v in validators]
        assert "RequiredMapKeyValidator" in types
        assert "IndexedGroupValidator" in types

    def test_loaded_validator_catches_incomplete_group(self):
        config = {
            "BioSampleMap": {
                "indexed_group": [
                    ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
                ]
            }
        }
        v = load_validators(config)[0]
        data = _make_data({"exposure1_names": "DrugName"})
        errors = v.validate(data)
        assert len(errors) == 2  # type and category missing

    def test_loaded_validator_passes_complete_group(self):
        config = {
            "BioSampleMap": {
                "indexed_group": [
                    ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
                ]
            }
        }
        v = load_validators(config)[0]
        data = _make_data({
            "exposure1_names": "DrugName",
            "exposure1_type": "DrugType",
            "exposure1_category": "DrugCat",
        })
        assert v.validate(data) == []


class TestLoadValidatorsMapSheets:

    def test_protein_map_required_uses_map_key_validator(self):
        config = {
            "ProteinMap": {
                "required": ["protein_id"],
            }
        }
        validators = load_validators(config)
        assert len(validators) == 1
        validator = validators[0]
        assert isinstance(validator, RequiredMapKeyValidator)
        assert validator.sheet == "ProteinMap"
        assert validator.field == "protein_id"

    def test_full_pounce_validator_config_loads(self):
        validators = load_validators("src/use_cases/pounce/pounce_validators.yaml")
        assert validators
