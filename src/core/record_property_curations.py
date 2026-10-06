"""Projection helpers for graph-neutral node property curations."""

from __future__ import annotations

import copy
import json
from typing import Any

from src.shared.metabolite_structure_chemistry import (
    calculate_metabolite_chem_props_derivatives,
    calculate_smiles_chemistry,
)


CURATION_ORIGINAL_FIELD = "_curation_original"
MISSING_ORIGINAL_MARKER = {"_odin_field_was_missing": True}

PROTECTED_FIELD_NAMES = frozenset({
    "id", "xref", "provenance", "sources", "start_node", "end_node",
    "_id", "_key", "_rev", "_from", "_to", "start_id", "end_id",
    "creation", "updates", "resolved_ids", "entity_resolution",
    CURATION_ORIGINAL_FIELD,
    "source", "source_id", "chem_data_source", "chem_source_id",
    "prefix", "structure_components",
})
METABOLITE_STRUCTURE_INPUT_FIELDS = frozenset({
    "iso_smiles", "isomeric_smiles", "canonical_smiles",
})


def is_protected_curation_field(field_name: str) -> bool:
    name = str(field_name or "")
    return (
        not name
        or name.startswith("_")
        or name in PROTECTED_FIELD_NAMES
        or name.startswith("calculated_")
        or name.startswith("derived_")
        or name.startswith("structure_calculation_")
    )


def canonical_path(path: list[Any]) -> str:
    return json.dumps(path, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def display_path(path: list[Any]) -> str:
    parts = []
    for token in path:
        if isinstance(token, str):
            parts.append(token)
        else:
            match = token.get("match") or {}
            selector = ", ".join(f"{key}={value}" for key, value in sorted(match.items()))
            parts[-1] = f"{parts[-1]}[{selector}]"
    return ".".join(parts)


def resolve_parent_and_field(document: dict, path: list[Any]) -> tuple[dict, str]:
    current: Any = document
    for index, token in enumerate(path):
        if isinstance(token, str):
            if is_protected_curation_field(token):
                raise ValueError(f"Field {token!r} is managed by the graph build")
            if index == len(path) - 1:
                if not isinstance(current, dict):
                    raise ValueError(f"Cannot address {display_path(path)} on a non-object value")
                return current, token
            if not isinstance(current, dict) or token not in current:
                raise ValueError(f"Curation path does not exist: {display_path(path)}")
            current = current[token]
            continue
        match = token.get("match") if isinstance(token, dict) else None
        if not isinstance(current, list) or not isinstance(match, dict) or not match:
            raise ValueError(f"Invalid list selector in curation path: {display_path(path)}")
        matches = [
            item for item in current
            if isinstance(item, dict)
            and all(item.get(key) == value for key, value in match.items())
        ]
        if len(matches) != 1:
            raise ValueError(
                f"Curation path selector matched {len(matches)} records: {display_path(path)}"
            )
        current = matches[0]
    raise ValueError("Curation path must end in a field")


def schema_for_path(schema_fields: dict, path: list[Any]) -> Any:
    current_fields = schema_fields or {}
    descriptor: Any = None
    for token in path:
        if isinstance(token, str):
            if is_protected_curation_field(token) or token not in current_fields:
                raise ValueError(f"Field is not declared curatable: {display_path(path)}")
            descriptor = current_fields[token]
            if isinstance(descriptor, dict) and descriptor.get("type") == "object":
                current_fields = descriptor.get("fields") or {}
            elif isinstance(descriptor, dict) and descriptor.get("type") == "list":
                current_fields = {}
            else:
                current_fields = {}
            continue
        if not isinstance(descriptor, dict) or descriptor.get("type") != "list":
            raise ValueError(f"List selector does not follow a declared list: {display_path(path)}")
        if descriptor.get("item_type") != "object":
            raise ValueError(f"List selector requires object items: {display_path(path)}")
        current_fields = descriptor.get("fields") or {}
        descriptor = {"type": "object", "fields": current_fields}
    if _aggregate_schema_contains_protected_fields(descriptor):
        raise ValueError(
            f"Aggregate edit could replace protected nested fields: {display_path(path)}"
        )
    return descriptor


def _aggregate_schema_contains_protected_fields(descriptor: Any) -> bool:
    if not isinstance(descriptor, dict):
        return False
    schema_type = descriptor.get("type")
    # Arango nested objects are schemaless. Replacing an aggregate object could
    # inject undeclared identity/provenance keys even when its declared schema
    # looks safe, so object values are only editable through declared leaf paths.
    if schema_type in {"dict", "object"}:
        return True
    if schema_type != "list":
        return False
    item_type = descriptor.get("item_type")
    return not (
        isinstance(item_type, str)
        and item_type in {"str", "int", "float", "bool", "date", "datetime"}
    )


def validate_value_for_schema(value: Any, descriptor: Any, path: list[Any]) -> None:
    if value is None:
        return
    schema_type = descriptor.get("type") if isinstance(descriptor, dict) else descriptor
    valid = True
    if schema_type == "bool":
        valid = type(value) is bool
    elif schema_type == "int":
        valid = type(value) is int
    elif schema_type == "float":
        valid = type(value) in {int, float}
    elif schema_type in {"str", "date", "datetime"}:
        valid = isinstance(value, str)
    elif schema_type == "list":
        valid = isinstance(value, list)
    elif schema_type in {"dict", "object"}:
        valid = isinstance(value, dict)
    if not valid:
        raise ValueError(
            f"Curated value for {display_path(path)} does not match schema type {schema_type}"
        )
    if schema_type == "list":
        item_descriptor = descriptor.get("item_type", "str")
        for item in value:
            validate_value_for_schema(item, item_descriptor, path)


def apply_record_property_decision(
    document: dict,
    decision,
) -> tuple[dict, dict]:
    """Return a projected document and one application report."""
    projected = copy.deepcopy(document)
    parent, field_name = resolve_parent_and_field(projected, decision.path)
    originals = parent.get(CURATION_ORIGINAL_FIELD)
    if originals is not None and not isinstance(originals, dict):
        raise ValueError(
            f"{CURATION_ORIGINAL_FIELD} must be an object at {display_path(decision.path)}"
        )
    originals = originals or {}
    stored_original = originals.get(field_name)
    original_was_missing = stored_original == MISSING_ORIGINAL_MARKER
    baseline_exists = (
        False
        if original_was_missing
        else field_name in originals or field_name in parent
    )
    baseline = (
        None
        if original_was_missing
        else stored_original if field_name in originals else parent.get(field_name)
    )
    report = {
        "curation_type": decision.curation_type,
        "batch_id": decision.batch_id,
        "operation_id": decision.source_operation.get("operation_id"),
        "action": decision.mode,
        "target": decision.target,
        "path": decision.path,
        "path_label": display_path(decision.path),
        "previous": parent.get(field_name),
        "result": decision.value if decision.mode == "set" else baseline,
    }
    if decision.mode == "remove_override":
        if field_name in originals:
            if original_was_missing:
                parent.pop(field_name, None)
            else:
                parent[field_name] = copy.deepcopy(stored_original)
            originals.pop(field_name, None)
            if originals:
                parent[CURATION_ORIGINAL_FIELD] = originals
            else:
                parent.pop(CURATION_ORIGINAL_FIELD, None)
        report["status"] = "restored"
        return projected, report
    observed_exists = getattr(decision, "observed_exists", True)
    if baseline_exists != observed_exists or (
        baseline_exists and baseline != decision.observed_value
    ):
        if baseline == decision.value:
            report["status"] = "redundant"
            return projected, report
        raise ValueError(
            f"Stale curation for {decision.target['model_type']}:{decision.target['id']} "
            f"{display_path(decision.path)}: observed {decision.observed_value!r}, "
            f"loaded {baseline!r}"
        )
    if field_name not in originals:
        originals[field_name] = (
            copy.deepcopy(baseline)
            if baseline_exists
            else copy.deepcopy(MISSING_ORIGINAL_MARKER)
        )
    parent[CURATION_ORIGINAL_FIELD] = originals
    parent[field_name] = copy.deepcopy(decision.value)
    update = format_curation_update(decision, baseline, decision.value)
    updates = list(projected.get("updates") or [])
    if update not in updates:
        updates.append(update)
    projected["updates"] = updates
    published_date = str(decision.published_at or "")[:10] or "None"
    curation_source = (
        f"Manual Curation\t{decision.batch_id}\t{published_date}\tNone"
    )
    sources = list(projected.get("sources") or [])
    if curation_source not in sources:
        sources.append(curation_source)
    projected["sources"] = sources
    report["status"] = "applied"
    return projected, report


def recalculate_curated_structure_derivatives(
    document: dict,
    *,
    model_type: str,
    decisions: list,
) -> dict:
    """Refresh protected structure derivatives after explicit source-field edits."""
    if model_type == "ChemicalEntity":
        recalculated = calculate_smiles_chemistry(
            document.get("smiles"), "smiles"
        )
        for field_name in (
            "calculated_mw",
            "calculated_monoisotopic_mass",
            "structure_calculation_input_field",
            "structure_calculation_method",
            "structure_calculation_method_version",
            "structure_calculation_error",
        ):
            document[field_name] = recalculated.get(field_name)
        document["structure_components"] = recalculated.get(
            "structure_components", []
        )
        return document

    if model_type != "MetaboliteIdentifier":
        return document

    touched_properties: dict[int, dict] = {}
    for decision in decisions:
        path = (
            decision.get("path")
            if isinstance(decision, dict)
            else getattr(decision, "path", None)
        )
        if not (
            isinstance(path, list)
            and len(path) == 3
            and path[0] == "chem_props"
            and isinstance(path[1], dict)
            and path[2] in METABOLITE_STRUCTURE_INPUT_FIELDS
        ):
            continue
        properties, _field_name = resolve_parent_and_field(document, path)
        touched_properties[id(properties)] = properties

    for properties in touched_properties.values():
        properties.update(calculate_metabolite_chem_props_derivatives(properties))
    return document


def format_curation_update(decision, previous: Any, current: Any) -> str:
    publisher = decision.published_by or {}
    curator = publisher.get("name") or publisher.get("id") or "unknown curator"
    note = str(decision.source_operation.get("note") or "").replace("\t", " ").replace("\n", " ")
    provenance = (
        f"curation batch={decision.batch_id}; operation="
        f"{decision.source_operation.get('operation_id') or ''}; curator={curator}; "
        f"published={decision.published_at}; rationale={note}"
    )
    return "\t".join([
        display_path(decision.path),
        json.dumps(previous, sort_keys=True, default=str),
        json.dumps(current, sort_keys=True, default=str),
        provenance,
        "CurationOverlay",
    ])


def project_record_decisions(document: dict, decisions: list, schema_fields: dict,
                             *, model_type: str) -> tuple[dict, list[dict]]:
    """Apply validated decisions to one effective record without graph mutation."""
    reports = []
    for decision in decisions:
        descriptor = schema_for_path(schema_fields, decision.path)
        if decision.mode == "set":
            validate_value_for_schema(decision.value, descriptor, decision.path)
        document, report = apply_record_property_decision(document, decision)
        reports.append(report)
    if decisions:
        document = recalculate_curated_structure_derivatives(
            document, model_type=model_type, decisions=decisions,
        )
    return document, reports
