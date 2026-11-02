"""Tests for ParallelListsValidator and the parallel_lists YAML loader."""

import pytest
from src.core.validator import ParallelListsValidator
from src.input_adapters.pounce_sheets.validator_loader import load_validators
from src.input_adapters.pounce_sheets.parsed_pounce_data import ParsedPounceData
from src.input_adapters.pounce_sheets.parsed_classes import ParsedProject


# ---------------------------------------------------------------------------
# ParallelListsValidator unit tests
# ---------------------------------------------------------------------------

def _make_data(owner_names=None, owner_emails=None,
               collaborator_names=None, collaborator_emails=None):
    """Build a minimal ParsedPounceData with a project."""
    proj = ParsedProject(
        owner_names=owner_names,
        owner_emails=owner_emails,
        collaborator_names=collaborator_names,
        collaborator_emails=collaborator_emails,
    )
    data = ParsedPounceData()
    data.project = proj
    return data


def _make_validator(fields=None, columns=None):
    fields = fields or ["owner_names", "owner_emails"]
    columns = columns or ["owner_name", "owner_email"]
    return ParallelListsValidator(
        entity="project",
        fields=fields,
        columns=columns,
        sheet="ProjectMeta",
    )


class TestParallelListsValidatorDirect:

    def test_equal_length_lists_pass(self):
        data = _make_data(owner_names=["Alice", "Bob"], owner_emails=["a@x.com", "b@x.com"])
        errors = _make_validator().validate(data)
        assert errors == []

    def test_single_entry_each_passes(self):
        data = _make_data(owner_names=["Alice"], owner_emails=["a@x.com"])
        errors = _make_validator().validate(data)
        assert errors == []

    def test_both_empty_passes(self):
        data = _make_data(owner_names=[], owner_emails=[])
        errors = _make_validator().validate(data)
        assert errors == []

    def test_both_none_passes(self):
        data = _make_data(owner_names=None, owner_emails=None)
        errors = _make_validator().validate(data)
        assert errors == []

    def test_length_mismatch_raises_error(self):
        data = _make_data(owner_names=["Alice", "Bob"], owner_emails=["a@x.com"])
        errors = _make_validator().validate(data)
        assert len(errors) == 1
        assert "owner_name=2" in errors[0].message
        assert "owner_email=1" in errors[0].message

    def test_error_has_correct_entity_and_sheet(self):
        data = _make_data(owner_names=["Alice", "Bob"], owner_emails=["a@x.com"])
        errors = _make_validator().validate(data)
        assert errors[0].entity == "project"
        assert errors[0].sheet == "ProjectMeta"
        assert errors[0].severity == "error"

    def test_one_none_one_populated_is_mismatch(self):
        # None treated as length 0; ["Alice"] is length 1
        data = _make_data(owner_names=["Alice"], owner_emails=None)
        errors = _make_validator().validate(data)
        assert len(errors) == 1

    def test_three_fields_all_equal_passes(self):
        proj = ParsedProject(
            owner_names=["Alice", "Bob"],
            owner_emails=["a@x.com", "b@x.com"],
            collaborator_names=["Carol", "Dave"],
        )
        data = ParsedPounceData()
        data.project = proj
        validator = ParallelListsValidator(
            entity="project",
            fields=["owner_names", "owner_emails", "collaborator_names"],
            columns=["owner_name", "owner_email", "collaborator_name"],
            sheet="ProjectMeta",
        )
        assert validator.validate(data) == []

    def test_three_fields_one_mismatch_raises_error(self):
        proj = ParsedProject(
            owner_names=["Alice", "Bob"],
            owner_emails=["a@x.com", "b@x.com"],
            collaborator_names=["Carol"],
        )
        data = ParsedPounceData()
        data.project = proj
        validator = ParallelListsValidator(
            entity="project",
            fields=["owner_names", "owner_emails", "collaborator_names"],
            columns=["owner_name", "owner_email", "collaborator_name"],
            sheet="ProjectMeta",
        )
        errors = validator.validate(data)
        assert len(errors) == 1

    def test_missing_entity_returns_no_errors(self):
        data = ParsedPounceData()  # data.project is None
        errors = _make_validator().validate(data)
        assert errors == []


# ---------------------------------------------------------------------------
# load_validators with parallel_lists
# ---------------------------------------------------------------------------

class TestLoadValidatorsParallelLists:

    def test_collaborator_email_is_optional_but_supplied_lists_must_match(self):
        validators = load_validators("src/use_cases/pounce/pounce_validators.yaml")
        required_email = [
            v for v in validators
            if v.sheet == "ProjectMeta" and v.column == "collaborator_email"
        ]
        assert required_email == []
        collaborator_list = next(
            v for v in validators
            if isinstance(v, ParallelListsValidator)
            and v.columns == ["collaborator_name", "collaborator_email"]
        )
        assert collaborator_list.validate(_make_data(
            collaborator_names=["Alice", "Bob"], collaborator_emails=[]
        )) == []
        assert len(collaborator_list.validate(_make_data(
            collaborator_names=["Alice", "Bob"], collaborator_emails=["alice@example.org"]
        ))) == 1

    def test_loads_parallel_lists_validator(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [
                    ["owner_name", "owner_email"],
                ]
            }
        }
        validators = load_validators(config)
        assert len(validators) == 1
        v = validators[0]
        assert isinstance(v, ParallelListsValidator)

    def test_loaded_validator_entity_and_sheet(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [["owner_name", "owner_email"]]
            }
        }
        v = load_validators(config)[0]
        assert v.entity == "project"
        assert v.sheet == "ProjectMeta"

    def test_loaded_validator_fields_resolved(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [["owner_name", "owner_email"]]
            }
        }
        v = load_validators(config)[0]
        # Python field names (not YAML keys)
        assert v.fields == ["owner_names", "owner_emails"]
        assert v.columns == ["owner_name", "owner_email"]

    def test_two_groups_produce_two_validators(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [
                    ["owner_name", "owner_email"],
                    ["collaborator_name", "collaborator_email"],
                ]
            }
        }
        validators = load_validators(config)
        assert len(validators) == 2
        assert all(isinstance(v, ParallelListsValidator) for v in validators)

    def test_parallel_lists_alongside_required(self):
        config = {
            "ProjectMeta": {
                "required": ["project_name"],
                "parallel_lists": [["owner_name", "owner_email"]],
            }
        }
        validators = load_validators(config)
        types = [type(v).__name__ for v in validators]
        assert "RequiredValidator" in types
        assert "ParallelListsValidator" in types

    def test_unknown_key_raises(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [["owner_name", "nonexistent_key"]]
            }
        }
        with pytest.raises(KeyError, match="nonexistent_key"):
            load_validators(config)

    def test_loaded_validator_catches_mismatch(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [["owner_name", "owner_email"]]
            }
        }
        v = load_validators(config)[0]
        data = _make_data(owner_names=["Alice", "Bob"], owner_emails=["a@x.com"])
        errors = v.validate(data)
        assert len(errors) == 1

    def test_loaded_validator_passes_on_match(self):
        config = {
            "ProjectMeta": {
                "parallel_lists": [["owner_name", "owner_email"]]
            }
        }
        v = load_validators(config)[0]
        data = _make_data(owner_names=["Alice", "Bob"], owner_emails=["a@x.com", "b@x.com"])
        assert v.validate(data) == []
