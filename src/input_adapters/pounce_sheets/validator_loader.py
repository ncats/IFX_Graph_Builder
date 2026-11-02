import importlib
from typing import List

import yaml

from src.core.object_factory import create_object_from_config
from src.core.validator import AllowedValuesValidator, ConditionalRequiredMapKeyValidator, ConditionalRequiredValidator, IndexedGroupValidator, ParallelListsValidator, RequiredMapKeyValidator, RequiredValidator, Validator
from src.input_adapters.pounce_sheets.constants import ExperimentWorkbook, ProjectWorkbook, StatsResultsWorkbook
from src.input_adapters.pounce_sheets.parsed_classes import (
    ParsedBiosample, ParsedBiospecimen, ParsedExperiment, ParsedExposure,
    ParsedGene, ParsedMetab, ParsedPeakDataMeta, ParsedProject, ParsedProtein,
    ParsedProteinDataMeta, ParsedRawDataMeta, ParsedRunBiosample, ParsedStatsResultsMeta,
)
from src.input_adapters.pounce_sheets.sheet_field import get_sheet_fields

# Map (sheet_name, ncatsdpi_key) -> (entity_name, python_field_name)
# Derived from sheet_field metadata on all parsed dataclasses.
_ENTITY_CLASSES = [
    ("project", ParsedProject),
    ("biosamples", ParsedBiosample),
    ("biospecimens", ParsedBiospecimen),
    ("exposures", ParsedExposure),
    ("experiments", ParsedExperiment),
    ("genes", ParsedGene),
    ("proteins", ParsedProtein),
    ("metabolites", ParsedMetab),
    ("peak_data_meta", ParsedPeakDataMeta),
    ("protein_data_meta", ParsedProteinDataMeta),
    ("raw_data_meta", ParsedRawDataMeta),
    ("run_biosamples", ParsedRunBiosample),
    ("stats_results", ParsedStatsResultsMeta),
]
_SHEET_KEY_MAP = {
    (meta["sheet"], meta["key"]): (entity_name, f.name)
    for entity_name, cls in _ENTITY_CLASSES
    for f, meta in get_sheet_fields(cls)
}

# Map sheets whose YAML entries represent configuration checks (key present in
# param_maps[sheet]) rather than content checks (row values non-null).
# Derived from constants so it stays in sync automatically.
_MAP_SHEETS: set = {
    ProjectWorkbook.BiosampleMapSheet.name,
    ExperimentWorkbook.RunBioSampleMapSheet.name,
    ExperimentWorkbook.MetabMapSheet.name,
    ExperimentWorkbook.GeneMapSheet.name,
    ExperimentWorkbook.ProteinMapSheet.name,
    StatsResultsWorkbook.EffectSizeMapSheet.name,
}

_DEFAULT_REQUIRED_MESSAGE = "Sheet {sheet} is missing required field: {field}"
_DEFAULT_MAP_REQUIRED_MESSAGE = "Map sheet {sheet} is missing configuration for: {field}"


def _resolve_allowed_values(values) -> List:
    """If values is a dict with an 'enum' key, import that class and return its display labels.

    YAML form::

        project_type:
          enum: src.models.pounce.project.ProjectType

    Works for both LabeledIntEnum (uses .label) and SimpleEnum (uses .value).
    """
    if isinstance(values, dict) and "enum" in values:
        module_path, class_name = values["enum"].rsplit(".", 1)
        module = importlib.import_module(module_path)
        cls = getattr(module, class_name)
        return [m.label if hasattr(m, "label") else m.value for m in cls]
    return values


def load_validators(yaml_file_or_dict) -> List[Validator]:
    if not isinstance(yaml_file_or_dict, dict):
        with open(yaml_file_or_dict) as f:
            yaml_file_or_dict = yaml.safe_load(f)
    config = yaml_file_or_dict

    required_message = config.get("required_message", _DEFAULT_REQUIRED_MESSAGE)

    validators = []
    for top_key, top_value in config.items():
        if top_key == "required_message":
            continue

        if top_key == "cross_entity":
            for entry in top_value:
                if "import" in entry:
                    validators.append(create_object_from_config(entry))
            continue

        # Sheet-based format: top_key is a sheet name
        if isinstance(top_value, dict) and (
            "required" in top_value or "allowed_values" in top_value
            or "parallel_lists" in top_value or "indexed_group" in top_value
            or "required_if" in top_value
        ):
            sheet = top_key

            for ncatsdpi_key in top_value.get("required", []):
                if sheet in _MAP_SHEETS:
                    # Map sheet: check that the key is configured (column mapped)
                    validators.append(RequiredMapKeyValidator(
                        entity="param_maps",
                        field=ncatsdpi_key,
                        message=_DEFAULT_MAP_REQUIRED_MESSAGE.format(sheet=sheet, field=ncatsdpi_key),
                        sheet=sheet,
                        column=ncatsdpi_key,
                    ))
                else:
                    if (sheet, ncatsdpi_key) not in _SHEET_KEY_MAP:
                        known_sheets = sorted({s for s, _ in _SHEET_KEY_MAP})
                        raise KeyError(
                            f"No parsed field for sheet='{sheet}', key='{ncatsdpi_key}'. "
                            f"Known data sheets: {known_sheets}. "
                            f"Note: use the *Meta sheet name (e.g. BioSampleMeta), not the Map sheet."
                        )
                    entity, field = _SHEET_KEY_MAP[(sheet, ncatsdpi_key)]
                    validators.append(RequiredValidator(
                        entity=entity,
                        field=field,
                        message=required_message.format(sheet=sheet, field=ncatsdpi_key),
                        sheet=sheet,
                        column=ncatsdpi_key,
                    ))

            for ncatsdpi_key, values in top_value.get("allowed_values", {}).items():
                values = _resolve_allowed_values(values)
                if (sheet, ncatsdpi_key) not in _SHEET_KEY_MAP:
                    known_sheets = sorted({s for s, _ in _SHEET_KEY_MAP})
                    raise KeyError(
                        f"No parsed field for sheet='{sheet}', key='{ncatsdpi_key}'. "
                        f"Known data sheets: {known_sheets}. "
                        f"Note: use the *Meta sheet name (e.g. BioSampleMeta), not the Map sheet."
                    )
                entity, field = _SHEET_KEY_MAP[(sheet, ncatsdpi_key)]
                validators.append(AllowedValuesValidator(
                    entity=entity,
                    field=field,
                    message=f"Field '{ncatsdpi_key}' in {sheet} must be one of: {', '.join(str(v) for v in values)}",
                    values=values,
                    sheet=sheet,
                    column=ncatsdpi_key,
                ))

            for group in top_value.get("parallel_lists", []):
                # A group may mark a field optional when its entire list is absent.
                if isinstance(group, dict):
                    keys = group["fields"]
                    optional_keys = set(group.get("optional_if_empty", []))
                    if not optional_keys.issubset(keys):
                        raise ValueError(
                            f"optional_if_empty must name fields in parallel_lists: {group}"
                        )
                else:
                    keys = group
                    optional_keys = set()
                resolved = []
                for ncatsdpi_key in keys:
                    if (sheet, ncatsdpi_key) not in _SHEET_KEY_MAP:
                        known_sheets = sorted({s for s, _ in _SHEET_KEY_MAP})
                        raise KeyError(
                            f"No parsed field for sheet='{sheet}', key='{ncatsdpi_key}'. "
                            f"Known data sheets: {known_sheets}. "
                            f"Note: use the *Meta sheet name (e.g. BioSampleMeta), not the Map sheet."
                        )
                    entity, python_field = _SHEET_KEY_MAP[(sheet, ncatsdpi_key)]
                    resolved.append((entity, python_field, ncatsdpi_key))

                entities = {r[0] for r in resolved}
                if len(entities) != 1:
                    raise ValueError(
                        f"All fields in a parallel_lists group must belong to the same entity. "
                        f"Sheet '{sheet}', group {group} resolved to entities: {entities}"
                    )

                validators.append(ParallelListsValidator(
                    entity=resolved[0][0],
                    fields=[r[1] for r in resolved],
                    columns=[r[2] for r in resolved],
                    sheet=sheet,
                    optional_if_empty=[field for _, field, key in resolved if key in optional_keys],
                ))

            for group in top_value.get("indexed_group", []):
                # group is a list of key templates with {}, e.g.
                # ["exposure{}_names", "exposure{}_type", "exposure{}_category"]
                prefixes = {t.split("{}")[0] for t in group}
                if len(prefixes) != 1:
                    raise ValueError(
                        f"All templates in an indexed_group must share the same prefix before '{{}}'. "
                        f"Sheet '{sheet}', group {group} has prefixes: {prefixes}"
                    )
                validators.append(IndexedGroupValidator(
                    sheet=sheet,
                    required_templates=group,
                ))

            for ncatsdpi_key, condition in top_value.get("required_if", {}).items():
                when_field_key = condition["when_field"]
                when_values = condition["when_values"]

                if sheet in _MAP_SHEETS:
                    # Map sheet: check the param_map config, condition on a different entity.
                    # Search all parsed fields for when_field_key to find its entity.
                    when_matches = [
                        (e, f) for (s, k), (e, f) in _SHEET_KEY_MAP.items()
                        if k == when_field_key
                    ]
                    if not when_matches:
                        raise KeyError(
                            f"No parsed field found for when_field='{when_field_key}' "
                            f"(used in required_if on map sheet '{sheet}')."
                        )
                    when_entity, when_python_field = when_matches[0]
                    validators.append(ConditionalRequiredMapKeyValidator(
                        sheet=sheet,
                        field=ncatsdpi_key,
                        when_entity=when_entity,
                        when_field=when_python_field,
                        when_values=when_values,
                        message=_DEFAULT_MAP_REQUIRED_MESSAGE.format(sheet=sheet, field=ncatsdpi_key),
                        column=ncatsdpi_key,
                    ))
                else:
                    if (sheet, ncatsdpi_key) not in _SHEET_KEY_MAP:
                        known_sheets = sorted({s for s, _ in _SHEET_KEY_MAP})
                        raise KeyError(
                            f"No parsed field for sheet='{sheet}', key='{ncatsdpi_key}'. "
                            f"Known data sheets: {known_sheets}. "
                            f"Note: use the *Meta sheet name (e.g. BioSampleMeta), not the Map sheet."
                        )
                    entity, field = _SHEET_KEY_MAP[(sheet, ncatsdpi_key)]

                    if (sheet, when_field_key) not in _SHEET_KEY_MAP:
                        known_sheets = sorted({s for s, _ in _SHEET_KEY_MAP})
                        raise KeyError(
                            f"No parsed field for sheet='{sheet}', key='{when_field_key}' (used in required_if.when_field). "
                            f"Known data sheets: {known_sheets}."
                        )
                    _, when_python_field = _SHEET_KEY_MAP[(sheet, when_field_key)]

                    validators.append(ConditionalRequiredValidator(
                        entity=entity,
                        field=field,
                        when_field=when_python_field,
                        when_values=when_values,
                        message=required_message.format(sheet=sheet, field=ncatsdpi_key),
                        sheet=sheet,
                        column=ncatsdpi_key,
                    ))

    return validators
