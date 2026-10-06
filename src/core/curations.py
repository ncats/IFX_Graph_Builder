"""Typed, reusable curation streams shared by ODIN builds and applications."""

from __future__ import annotations

import hashlib
import json
import math
import re
import types
from dataclasses import dataclass, field, fields
from typing import Any, Iterable, Optional, Union, get_args, get_origin, get_type_hints


FORMAT_VERSION = 2
PUBLISHED_PREFIX = "curations/v2"

METABOLITE_EQUIVALENCE_EDGES = "metabolite_equivalence_edges"
METABOLITE_EXPECTED_CLIQUES = "metabolite_expected_cliques"
METABOLITE_RECORD_PROPERTIES = "metabolite_record_properties"
CHEBI_RECORD_PROPERTIES = "chebi_record_properties"
METABOLITE_RECORD_SUPPRESSIONS = "metabolite_record_suppressions"
METABOLITE_MW_ADJUDICATIONS = "metabolite_mw_adjudications"
MW_REVIEW_EVIDENCE_VERSION = "mw-review-evidence-v1"


@dataclass(frozen=True)
class CuratablePropertyDefinition:
    model_type: str
    property_name: str
    value_type: type
    nullable: bool

    @property
    def value_type_name(self) -> str:
        return {
            bool: "boolean",
            str: "string",
            int: "integer",
            float: "number",
        }[self.value_type]

    def metadata(self) -> dict:
        return {
            "model_type": self.model_type,
            "name": self.property_name,
            "label": self.property_name.replace("_", " ").strip().title(),
            "value_type": self.value_type_name,
            "nullable": self.nullable,
        }


@dataclass(frozen=True)
class CurationTypeDefinition:
    id: str
    actions: frozenset[str]
    model_types: tuple[str, ...] = ()
    edge_types: tuple[str, ...] = ()
    operation_contract: Optional[str] = None
    fixed_curation_set: Optional[str] = None


_FRAMEWORK_PROPERTY_DENYLIST = frozenset({
    "id", "xref", "provenance", "sources", "start_node", "end_node",
    "_id", "_key", "_rev", "_from", "_to", "start_id", "end_id",
    "creation", "updates", "resolved_ids", "entity_resolution",
    "_curation_original",
})
_MODEL_PROPERTY_DENYLIST = {
    "MetaboliteIdentifier": frozenset({"prefix", "hmdb_status"}),
}
_MODEL_PROPERTY_ALLOWLIST = {
    "ChemicalEntity": frozenset({
        "charge", "formula", "inchi", "inchi_key", "mass",
        "monoisotopic_mass", "smiles", "wurcs",
    }),
}
_SUPPORTED_PROPERTY_TYPES = frozenset({bool, str, int, float})


def _curatable_model_class(model_type: str):
    if model_type == "MetaboliteIdentifier":
        from src.models.metabolite_harmonization import MetaboliteIdentifier
        return MetaboliteIdentifier
    if model_type == "ChemicalEntity":
        from src.models.chebi import ChemicalEntity
        return ChemicalEntity
    return None


def _unwrap_optional(type_hint) -> tuple[Any, bool]:
    origin = get_origin(type_hint)
    if origin in (Union, types.UnionType):
        arguments = get_args(type_hint)
        non_none = tuple(argument for argument in arguments if argument is not type(None))
        if len(non_none) == 1 and len(non_none) != len(arguments):
            return non_none[0], True
    return type_hint, False


def curatable_property_definitions(model_type: str) -> tuple[CuratablePropertyDefinition, ...]:
    model_class = _curatable_model_class(model_type)
    if model_class is None:
        raise ValueError(f"Unsupported curation target model: {model_type!r}")
    type_hints = get_type_hints(model_class)
    denied = _FRAMEWORK_PROPERTY_DENYLIST | _MODEL_PROPERTY_DENYLIST.get(model_type, frozenset())
    definitions = []
    for model_field in fields(model_class):
        if model_field.name.startswith("_") or model_field.name in denied:
            continue
        allowlist = _MODEL_PROPERTY_ALLOWLIST.get(model_type)
        if allowlist is not None and model_field.name not in allowlist:
            continue
        value_type, nullable = _unwrap_optional(type_hints.get(model_field.name, model_field.type))
        if value_type not in _SUPPORTED_PROPERTY_TYPES:
            continue
        definitions.append(CuratablePropertyDefinition(
            model_type=model_type,
            property_name=model_field.name,
            value_type=value_type,
            nullable=nullable,
        ))
    return tuple(definitions)


def curatable_property_definition(
    model_type: str,
    property_name: str,
) -> Optional[CuratablePropertyDefinition]:
    return next((
        definition
        for definition in curatable_property_definitions(model_type)
        if definition.property_name == property_name
    ), None)


CURATION_TYPES = {
    CHEBI_RECORD_PROPERTIES: CurationTypeDefinition(
        id=CHEBI_RECORD_PROPERTIES,
        actions=frozenset({"set_properties"}),
        model_types=("ChemicalEntity",),
        operation_contract="record_properties",
        fixed_curation_set="chebi",
    ),
    METABOLITE_RECORD_PROPERTIES: CurationTypeDefinition(
        id=METABOLITE_RECORD_PROPERTIES,
        actions=frozenset({"set_properties"}),
        model_types=("MetaboliteIdentifier",),
        operation_contract="record_properties",
    ),
    METABOLITE_EQUIVALENCE_EDGES: CurationTypeDefinition(
        id=METABOLITE_EQUIVALENCE_EDGES,
        actions=frozenset({"remove_edge", "retain_edge"}),
        edge_types=("MetaboliteIdentifierMappingEdge",),
    ),
    METABOLITE_EXPECTED_CLIQUES: CurationTypeDefinition(
        id=METABOLITE_EXPECTED_CLIQUES,
        actions=frozenset({"assert_same_clique", "retire_assertion"}),
    ),
    METABOLITE_RECORD_SUPPRESSIONS: CurationTypeDefinition(
        id=METABOLITE_RECORD_SUPPRESSIONS,
        actions=frozenset({"suppress_record", "restore_record"}),
        model_types=("MetaboliteIdentifier",),
    ),
    METABOLITE_MW_ADJUDICATIONS: CurationTypeDefinition(
        id=METABOLITE_MW_ADJUDICATIONS,
        actions=frozenset({"accept_mw_discrepancy", "reopen_mw_discrepancy"}),
    ),
}


def is_record_property_type(curation_type: str) -> bool:
    return validate_curation_type(curation_type).operation_contract == "record_properties"


def record_property_type_for_model(model_type: str) -> Optional[str]:
    matches = [
        definition.id
        for definition in CURATION_TYPES.values()
        if definition.operation_contract == "record_properties"
        and model_type in definition.model_types
    ]
    if len(matches) > 1:
        raise ValueError(
            f"Multiple record-property curation types own model {model_type!r}: "
            + ", ".join(sorted(matches))
        )
    return matches[0] if matches else None


def record_property_curation_set_for_model(
    model_type: str,
    default_curation_set: str,
) -> str:
    curation_type = record_property_type_for_model(model_type)
    if not curation_type:
        return default_curation_set
    return CURATION_TYPES[curation_type].fixed_curation_set or default_curation_set


def validate_curation_type(curation_type: str) -> CurationTypeDefinition:
    definition = CURATION_TYPES.get(str(curation_type or "").strip())
    if definition is None:
        raise ValueError(f"Unsupported curation type: {curation_type!r}")
    return definition


def validate_operation(curation_type: str, operation: dict) -> dict:
    definition = validate_curation_type(curation_type)
    if not isinstance(operation, dict):
        raise ValueError("Curation operation must be an object")
    action = str(operation.get("action") or "").strip()
    if action not in definition.actions:
        raise ValueError(f"Action {action!r} is not supported by {curation_type}")

    if is_record_property_type(curation_type):
        _validate_record_property_operation(operation, definition)
    if action in {"remove_edge", "retain_edge"}:
        edge_type = str(operation.get("edge_type") or "").strip()
        left = str(operation.get("start_id") or "").strip()
        right = str(operation.get("end_id") or "").strip()
        if edge_type not in definition.edge_types:
            raise ValueError(f"Edge type {edge_type!r} is not registered for {curation_type}")
        if not left or not right or left == right:
            raise ValueError("Edge curation requires two different endpoint ids")
        if operation.get("symmetric") is not True:
            raise ValueError(f"{edge_type} curation requires symmetric=true")
    if action in {"suppress_record", "restore_record"}:
        target = operation.get("target")
        if not isinstance(target, dict) or target.get("kind") != "node":
            raise ValueError("Record suppression curations require a node target")
        model_type = str(target.get("model_type") or "").strip()
        target_id = str(target.get("id") or "").strip()
        if model_type not in definition.model_types:
            raise ValueError(f"Model type {model_type!r} is not registered for {curation_type}")
        if not target_id:
            raise ValueError("Record suppression target id is required")
        note = str(operation.get("note") or "").strip()
        if action == "suppress_record" and not note:
            raise ValueError("Suppressing a record requires a rationale")
    if action in {"accept_mw_discrepancy", "reopen_mw_discrepancy"}:
        target = operation.get("target")
        if not isinstance(target, dict) or target.get("kind") != "validation_finding":
            raise ValueError("MW adjudications require a validation-finding target")
        required_target = {
            "curation_set": "metabolite_harmonization",
            "check": "mw_spread",
        }
        for key, expected in required_target.items():
            if target.get(key) != expected:
                raise ValueError(f"MW adjudication target {key} must be {expected!r}")
        finding_id = str(target.get("finding_id") or "").strip()
        anchor_id = str(target.get("anchor_id") or "").strip()
        if not re.fullmatch(r"mw-[0-9a-f]{24}", finding_id):
            raise ValueError("MW adjudication requires a valid finding id")
        if not anchor_id or ":" not in anchor_id:
            raise ValueError("MW adjudication requires a stable anchor identifier")
        fingerprint = str(operation.get("observed_evidence_fingerprint") or "").strip()
        if not re.fullmatch(r"[0-9a-f]{64}", fingerprint):
            raise ValueError("MW adjudication requires the reviewed evidence fingerprint")
        member_ids = operation.get("observed_member_ids")
        if not isinstance(member_ids, list) or anchor_id not in member_ids:
            raise ValueError("MW adjudication requires the complete reviewed member list")
        if len(member_ids) != len(set(member_ids)) or any(
            not isinstance(identifier, str) or ":" not in identifier
            for identifier in member_ids
        ):
            raise ValueError("MW adjudication member identifiers must be unique prefixed strings")
        note = str(operation.get("note") or "").strip()
        if not note:
            raise ValueError("MW adjudication requires an explanatory note")
        if action == "accept_mw_discrepancy":
            if not str(operation.get("reason") or "").strip():
                raise ValueError("Accepting an MW discrepancy requires a reason")
            supporting_ids = operation.get("supporting_ids") or []
            if not isinstance(supporting_ids, list) or any(
                identifier not in member_ids for identifier in supporting_ids
            ):
                raise ValueError("Supporting identifiers must belong to the reviewed finding")
        if operation.get("observed_evidence_snapshot") is not None:
            _validate_mw_review_evidence_snapshot(
                operation["observed_evidence_snapshot"], member_ids
            )
    return operation


def _validate_mw_review_evidence_snapshot(snapshot: object, member_ids: list[str]) -> None:
    if not isinstance(snapshot, dict):
        raise ValueError("MW reviewed evidence snapshot must be an object")
    if snapshot.get("version") != MW_REVIEW_EVIDENCE_VERSION:
        raise ValueError("MW reviewed evidence snapshot has an unsupported version")
    validator_version = snapshot.get("validator_version")
    if validator_version is not None and not isinstance(validator_version, str):
        raise ValueError("MW reviewed evidence validator version must be text or null")
    threshold = snapshot.get("threshold")
    decimal_pattern = r"-?(?:0|[1-9][0-9]*)(?:\.[0-9]+)?"
    if threshold is not None and (
        not isinstance(threshold, str)
        or not re.fullmatch(decimal_pattern, threshold)
    ):
        raise ValueError("MW reviewed evidence threshold must be a decimal string or null")

    mass_observations = snapshot.get("mass_observations")
    if not isinstance(mass_observations, list):
        raise ValueError("MW reviewed evidence mass observations must be a list")
    normalized_masses = []
    for observation in mass_observations:
        if not isinstance(observation, dict):
            raise ValueError("MW reviewed mass observations must be objects")
        member_id = observation.get("member_id")
        channel = observation.get("channel")
        value = observation.get("value")
        if member_id not in member_ids:
            raise ValueError("MW reviewed mass observation member must belong to the finding")
        if channel not in {"average", "monoisotopic", "unspecified"}:
            raise ValueError("MW reviewed mass observation has an unsupported channel")
        if not isinstance(value, str) or not re.fullmatch(decimal_pattern, value):
            raise ValueError("MW reviewed mass observation value must be a decimal string")
        normalized_masses.append((member_id, channel, value))
    if normalized_masses != sorted(set(normalized_masses)):
        raise ValueError("MW reviewed mass observations must be unique and sorted")

    component_matches = snapshot.get("component_matches")
    if not isinstance(component_matches, list):
        raise ValueError("MW reviewed component matches must be a list")
    normalized_components = []
    for match in component_matches:
        if not isinstance(match, dict):
            raise ValueError("MW reviewed component matches must be objects")
        channel = match.get("channel")
        whole_member_id = match.get("whole_member_id")
        component_member_id = match.get("component_member_id")
        whole_mass = match.get("whole_mass")
        component_mass = match.get("component_mass")
        if channel not in {"average", "monoisotopic", "unspecified"}:
            raise ValueError("MW reviewed component match has an unsupported channel")
        if whole_member_id not in member_ids or component_member_id not in member_ids:
            raise ValueError("MW reviewed component match members must belong to the finding")
        if any(
            not isinstance(value, str) or not re.fullmatch(decimal_pattern, value)
            for value in (whole_mass, component_mass)
        ):
            raise ValueError("MW reviewed component masses must be decimal strings")
        formula = match.get("component_formula")
        if formula is not None and not isinstance(formula, str):
            raise ValueError("MW reviewed component formula must be text or null")
        normalized_components.append((
            channel, whole_member_id, whole_mass,
            component_member_id, component_mass, formula or "",
        ))
    if normalized_components != sorted(set(normalized_components)):
        raise ValueError("MW reviewed component matches must be unique and sorted")


def _validate_record_property_operation(
    operation: dict,
    definition: CurationTypeDefinition,
) -> None:
    target = operation.get("target")
    if not isinstance(target, dict) or target.get("kind") != "node":
        raise ValueError("Record property curations require a node target")
    for field_name, label in (
        ("curation_set", "curation set"),
        ("model_type", "target model type"),
        ("id", "target id"),
    ):
        if not str(target.get(field_name) or "").strip():
            raise ValueError(f"Record property curation {label} is required")
    if target["model_type"] not in definition.model_types:
        raise ValueError(
            f"Model type {target['model_type']!r} is not registered for {definition.id}"
        )
    if (
        definition.fixed_curation_set is not None
        and target["curation_set"] != definition.fixed_curation_set
    ):
        raise ValueError(
            f"{definition.id} curation set must be {definition.fixed_curation_set!r}"
        )
    if definition.id == CHEBI_RECORD_PROPERTIES and not re.fullmatch(
        r"CHEBI:[0-9]+", str(target.get("id") or "")
    ):
        raise ValueError("ChEBI record-property curations require a CHEBI identifier")
    note = str(operation.get("note") or "").strip()
    if not note:
        raise ValueError("Record property curation requires a rationale")
    decisions = operation.get("decisions")
    if not isinstance(decisions, list) or not decisions:
        raise ValueError("Record property curation requires at least one field decision")
    canonical_paths = []
    for decision in decisions:
        if not isinstance(decision, dict):
            raise ValueError("Record property decisions must be objects")
        path = decision.get("path")
        _validate_record_property_path(path)
        allowlist = _MODEL_PROPERTY_ALLOWLIST.get(target["model_type"])
        if allowlist is not None and path[0] not in allowlist:
            raise ValueError(
                f"Field {path[0]!r} is not curatable for {target['model_type']}"
            )
        canonical_path = canonical_json(path)
        if canonical_path in canonical_paths:
            raise ValueError("Record property decisions must not repeat a field path")
        canonical_paths.append(canonical_path)
        mode = str(decision.get("mode") or "").strip()
        if mode not in {"set", "remove_override"}:
            raise ValueError("Record property decision mode must be set or remove_override")
        if mode == "set" and "value" not in decision:
            raise ValueError("Setting a record property requires a value")
        if mode == "remove_override" and "value" in decision:
            raise ValueError("Removing a record property override must not include a value")
        if mode == "set" and "observed_value" not in decision:
            raise ValueError("Setting a record property requires the observed loaded value")
        if "observed_exists" in decision and type(decision["observed_exists"]) is not bool:
            raise ValueError("Record property observed_exists must be true or false")


def _validate_record_property_path(path: Any) -> None:
    from src.core.record_property_curations import is_protected_curation_field

    if not isinstance(path, list) or not path:
        raise ValueError("Record property decision path must be a non-empty list")
    expect_field = True
    for token in path:
        if expect_field:
            if not isinstance(token, str) or not token.strip():
                raise ValueError("Record property field path tokens must be non-empty strings")
            field_name = token.strip()
            if is_protected_curation_field(field_name):
                raise ValueError(f"Field {field_name!r} is managed by the graph build")
            expect_field = False
            continue
        if isinstance(token, str):
            field_name = token.strip()
            if not field_name:
                raise ValueError("Record property field path tokens must be non-empty strings")
            if is_protected_curation_field(field_name):
                raise ValueError(f"Field {field_name!r} is managed by the graph build")
            continue
        if not isinstance(token, dict) or set(token) != {"match"}:
            raise ValueError("Nested list paths require an exact-match selector")
        match = token.get("match")
        if not isinstance(match, dict) or not match:
            raise ValueError("Nested list selectors require at least one match field")
        if any(
            not isinstance(key, str) or not key.strip() or key.startswith("_")
            for key in match
        ):
            raise ValueError("Nested list selector fields must be named public fields")
        expect_field = True
    if expect_field:
        raise ValueError("Record property path cannot end with a list selector")


def _validate_property_value(definition: CuratablePropertyDefinition, value: Any) -> None:
    if value is None:
        if not definition.nullable:
            raise ValueError(
                f"Property {definition.model_type}.{definition.property_name} does not allow null"
            )
        return
    expected = definition.value_type
    if expected is float:
        valid = type(value) in {int, float} and math.isfinite(float(value))
    else:
        valid = type(value) is expected
    if not valid:
        raise ValueError(
            f"Property {definition.model_type}.{definition.property_name} requires "
            f"{definition.value_type_name}{' or null' if definition.nullable else ''}"
        )


def operation_subject(curation_type: str, operation: dict) -> tuple[str, ...]:
    """Return the stable subject whose latest operation wins."""
    validate_operation(curation_type, operation)
    action = operation["action"]
    if action == "set_properties":
        raise ValueError("set_properties has multiple subjects; use operation_subjects")
    if action in {"remove_edge", "retain_edge"}:
        edge_type = operation["edge_type"]
        left, right = sorted((operation["start_id"], operation["end_id"]))
        return ("edge", edge_type, left, right)
    if action in {"suppress_record", "restore_record"}:
        target = operation["target"]
        return ("record_harmonization", target["model_type"], target["id"])
    if action in {"accept_mw_discrepancy", "reopen_mw_discrepancy"}:
        target = operation["target"]
        return (
            "validation_finding",
            target["curation_set"],
            target["check"],
            target["anchor_id"],
        )
    assertion_id = str(operation.get("assertion_id") or "").strip()
    if not assertion_id:
        raise ValueError("Assertion curation requires assertion_id")
    return ("assertion", assertion_id)


def operation_subjects(curation_type: str, operation: dict) -> tuple[tuple[str, ...], ...]:
    validate_operation(curation_type, operation)
    if is_record_property_type(curation_type):
        target = operation["target"]
        return tuple(
            (
                "record_property",
                target["curation_set"],
                target["kind"],
                target["model_type"],
                target["id"],
                canonical_json(decision["path"]),
            )
            for decision in operation["decisions"]
        )
    return (operation_subject(curation_type, operation),)


@dataclass(frozen=True)
class ResolvedRecordPropertyDecision:
    curation_type: str
    target: dict
    path: list[Any]
    mode: str
    value: Any
    observed_value: Any
    observed_exists: bool
    batch_id: str
    published_at: str
    published_by: Optional[dict]
    source_operation: dict

    @property
    def subject(self) -> tuple[str, ...]:
        return (
            "record_property",
            self.target["curation_set"],
            self.target["kind"],
            self.target["model_type"],
            self.target["id"],
            canonical_json(self.path),
        )


def record_property_decisions_from_operation(
    curation_type: str,
    operation: dict,
    *,
    batch_id: str,
    published_at: str,
    published_by: Optional[dict],
) -> tuple[ResolvedRecordPropertyDecision, ...]:
    validate_operation(curation_type, operation)
    if not is_record_property_type(curation_type):
        return ()
    return tuple(
        ResolvedRecordPropertyDecision(
            curation_type=curation_type,
            target=dict(operation["target"]),
            path=list(decision["path"]),
            mode=decision["mode"],
            value=decision.get("value"),
            observed_value=decision.get("observed_value"),
            observed_exists=decision.get("observed_exists", True),
            batch_id=batch_id,
            published_at=published_at,
            published_by=published_by,
            source_operation=dict(operation),
        )
        for decision in operation["decisions"]
    )


def canonical_json(payload: Any) -> str:
    return json.dumps(payload, sort_keys=True, separators=(",", ":"), ensure_ascii=False)


def payload_sha256(payload: Any) -> str:
    return hashlib.sha256(canonical_json(payload).encode("utf-8")).hexdigest()


def batch_key(curation_type: str, batch_id: str) -> str:
    validate_curation_type(curation_type)
    if not re.fullmatch(r"[A-Za-z0-9._-]+", str(batch_id or "")):
        raise ValueError("Invalid curation batch id")
    return f"{PUBLISHED_PREFIX}/{curation_type}/batches/{batch_id}.json"


def manifest_key(curation_type: str) -> str:
    validate_curation_type(curation_type)
    return f"{PUBLISHED_PREFIX}/{curation_type}/manifest.json"


@dataclass(frozen=True)
class ResolvedCurationOperation:
    curation_type: str
    operation: dict
    batch_id: str
    published_at: str
    published_by: Optional[dict]


@dataclass(frozen=True)
class ResolvedCurationBatch:
    """Portable publication and origin metadata for one immutable batch."""

    batch_id: str
    name: str
    description: str
    created_at: str
    published_at: str
    created_by: Optional[dict]
    source: Optional[dict]
    object_key: str
    sha256: str


@dataclass
class CurationSnapshot:
    curation_type: str
    manifest_revision: int
    manifest_hash: str
    batch_ids: list[str]
    batch_hashes: list[str]
    operations: list[ResolvedCurationOperation]
    fingerprint: str
    source_uri: str
    batches_by_id: dict[str, ResolvedCurationBatch] = field(default_factory=dict)
    _active_by_subject: dict[tuple[str, ...], ResolvedCurationOperation] = field(
        default_factory=dict,
        repr=False,
    )
    _active_record_property_decisions: dict[
        tuple[str, ...], ResolvedRecordPropertyDecision
    ] = field(default_factory=dict, repr=False)

    @property
    def active_operations(self) -> list[ResolvedCurationOperation]:
        return list(self._active_by_subject.values())

    @property
    def active_record_property_decisions(self) -> list[ResolvedRecordPropertyDecision]:
        return list(self._active_record_property_decisions.values())

    def record_property_decisions_for_set(
        self, curation_set: str
    ) -> list[ResolvedRecordPropertyDecision]:
        return [
            decision
            for decision in self.active_record_property_decisions
            if decision.target.get("curation_set") == curation_set
        ]

    def metadata(self) -> dict:
        return {
            "curation_type": self.curation_type,
            "manifest_revision": self.manifest_revision,
            "manifest_hash": self.manifest_hash,
            "batch_ids": self.batch_ids,
            "batch_hashes": self.batch_hashes,
            "resolved_operation_fingerprint": self.fingerprint,
            "active_operation_count": (
                len(self._active_by_subject)
                + len(self._active_record_property_decisions)
            ),
            "source_uri": self.source_uri,
        }

    def batch(self, batch_id: str) -> Optional[ResolvedCurationBatch]:
        return self.batches_by_id.get(batch_id)


def _read_json(storage, key: str) -> tuple[dict, str]:
    raw = storage.read_text(key)
    try:
        payload = json.loads(raw)
    except json.JSONDecodeError as exc:
        raise ValueError(f"Invalid curation JSON at s3://{storage.bucket}/{key}: {exc}") from exc
    if not isinstance(payload, dict):
        raise ValueError(f"Curation JSON must be an object at s3://{storage.bucket}/{key}")
    return payload, raw


def _is_missing_object_error(exc: Exception) -> bool:
    error_code = str(getattr(exc, "response", {}).get("Error", {}).get("Code", ""))
    return isinstance(exc, KeyError) or error_code in {"404", "NoSuchKey", "NotFound"}


def resolve_curation_type(storage, curation_type: str, *, allow_missing: bool = False) -> CurationSnapshot:
    validate_curation_type(curation_type)
    key = manifest_key(curation_type)
    try:
        manifest, _raw_manifest = _read_json(storage, key)
    except Exception as exc:
        if not _is_missing_object_error(exc):
            raise
        if not allow_missing:
            raise ValueError(f"Missing curation manifest: s3://{storage.bucket}/{key}")
        manifest = {
            "format_version": FORMAT_VERSION,
            "curation_type": curation_type,
            "revision": 0,
            "batches": [],
        }
    if manifest.get("format_version") != FORMAT_VERSION:
        raise ValueError(f"Unsupported curation manifest version in {key}")
    if manifest.get("curation_type") != curation_type:
        raise ValueError(f"Curation manifest type mismatch in {key}")

    return _resolve_curation_manifest(storage, curation_type, manifest)


def replay_curation_snapshot(storage, metadata: dict) -> CurationSnapshot:
    """Resolve exact ordered batches recorded by a stage, never the latest manifest."""
    curation_type = metadata["curation_type"]
    validate_curation_type(curation_type)
    ids, hashes = metadata["batch_ids"], metadata["batch_hashes"]
    if len(ids) != len(hashes) or len(set(ids)) != len(ids):
        raise ValueError("Invalid recorded curation batch list")
    expected_uri = f"s3://{storage.bucket}/{manifest_key(curation_type)}"
    if metadata.get("source_uri") != expected_uri:
        raise ValueError(f"Recorded curation storage differs from {expected_uri}")
    manifest = {
        "format_version": FORMAT_VERSION,
        "curation_type": curation_type,
        "revision": metadata["manifest_revision"],
        "batches": [{"batch_id": bid, "sha256": sha} for bid, sha in zip(ids, hashes)],
    }
    snapshot = _resolve_curation_manifest(storage, curation_type, manifest)
    if snapshot.fingerprint != metadata["resolved_operation_fingerprint"]:
        raise ValueError(f"Recorded curation fingerprint mismatch: {curation_type}")
    snapshot.manifest_hash = metadata["manifest_hash"]
    return snapshot


def _resolve_curation_manifest(storage, curation_type: str, manifest: dict) -> CurationSnapshot:
    key = manifest_key(curation_type)
    operations: list[ResolvedCurationOperation] = []
    batch_ids: list[str] = []
    batch_hashes: list[str] = []
    active_by_subject: dict[tuple[str, ...], ResolvedCurationOperation] = {}
    active_record_property_decisions: dict[
        tuple[str, ...], ResolvedRecordPropertyDecision
    ] = {}
    batches_by_id: dict[str, ResolvedCurationBatch] = {}
    for entry in manifest.get("batches") or []:
        batch_id = str(entry.get("batch_id") or "").strip()
        object_key = str(entry.get("object_key") or batch_key(curation_type, batch_id))
        batch, _raw_batch = _read_json(storage, object_key)
        actual_hash = payload_sha256(batch)
        expected_hash = str(entry.get("sha256") or "")
        if not expected_hash or actual_hash != expected_hash:
            raise ValueError(f"Curation batch hash mismatch: s3://{storage.bucket}/{object_key}")
        if batch.get("format_version") != FORMAT_VERSION or batch.get("curation_type") != curation_type:
            raise ValueError(f"Curation batch contract mismatch: s3://{storage.bucket}/{object_key}")
        if batch.get("curation_batch_id") != batch_id:
            raise ValueError(f"Curation batch id mismatch: s3://{storage.bucket}/{object_key}")
        batch_ids.append(batch_id)
        batch_hashes.append(actual_hash)
        batches_by_id[batch_id] = ResolvedCurationBatch(
            batch_id=batch_id,
            name=str(batch.get("name") or ""),
            description=str(batch.get("description") or ""),
            created_at=str(batch.get("created_at") or ""),
            published_at=str(batch.get("published_at") or batch.get("created_at") or ""),
            created_by=(dict(batch["created_by"]) if isinstance(batch.get("created_by"), dict) else None),
            source=(dict(batch["source"]) if isinstance(batch.get("source"), dict) else None),
            object_key=object_key,
            sha256=actual_hash,
        )
        for operation in batch.get("operations") or []:
            validate_operation(curation_type, operation)
            resolved = ResolvedCurationOperation(
                curation_type=curation_type,
                operation=dict(operation),
                batch_id=batch_id,
                published_at=str(batch.get("published_at") or batch.get("created_at") or ""),
                published_by=batch.get("created_by"),
            )
            operations.append(resolved)
            if is_record_property_type(curation_type):
                for decision in record_property_decisions_from_operation(
                    curation_type,
                    operation,
                    batch_id=batch_id,
                    published_at=resolved.published_at,
                    published_by=resolved.published_by,
                ):
                    active_record_property_decisions[decision.subject] = decision
            else:
                active_by_subject[operation_subject(curation_type, operation)] = resolved

    fingerprint_payload = [
        {
            "subject": list(subject),
            "batch_id": resolved.batch_id,
            "operation": resolved.operation,
        }
        for subject, resolved in sorted(active_by_subject.items())
    ]
    fingerprint_payload.extend({
        "subject": list(subject),
        "batch_id": decision.batch_id,
        "decision": {
            "mode": decision.mode,
            "value": decision.value,
            "observed_value": decision.observed_value,
            "observed_exists": decision.observed_exists,
        },
    } for subject, decision in sorted(active_record_property_decisions.items()))
    fingerprint_payload.sort(key=lambda item: item["subject"])
    return CurationSnapshot(
        curation_type=curation_type,
        manifest_revision=int(manifest.get("revision") or 0),
        manifest_hash=payload_sha256(manifest),
        batch_ids=batch_ids,
        batch_hashes=batch_hashes,
        operations=operations,
        fingerprint=payload_sha256(fingerprint_payload),
        source_uri=f"s3://{storage.bucket}/{key}",
        batches_by_id=batches_by_id,
        _active_by_subject=active_by_subject,
        _active_record_property_decisions=active_record_property_decisions,
    )


def resolve_curation_types(
    storage,
    curation_types: Iterable[str],
    *,
    allow_missing: bool = False,
) -> dict[str, CurationSnapshot]:
    return {
        curation_type: resolve_curation_type(storage, curation_type, allow_missing=allow_missing)
        for curation_type in dict.fromkeys(curation_types)
    }
