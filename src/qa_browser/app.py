import argparse
import copy
import csv
import gzip
import hashlib
import importlib.util
import inspect
import io
import json
import logging
import math
import os
import re
import shutil
import sys
import threading
import time
import uuid
from concurrent.futures import ThreadPoolExecutor
from pathlib import Path
from datetime import datetime, timedelta, timezone
from typing import Any, Optional, Dict, Iterable, List
from urllib.parse import quote as url_quote, urlencode

import urllib3
import yaml
from arango import ArangoClient
from fastapi import FastAPI, Request, Form, HTTPException
from fastapi.responses import HTMLResponse, StreamingResponse, RedirectResponse, FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.templating import Jinja2Templates
from sqlalchemy import create_engine, text, inspect as sa_inspect
from sqlalchemy.engine import Engine
from starlette.concurrency import run_in_threadpool

import uvicorn

logger = logging.getLogger(__name__)

from src.core.curations import (
    CHEBI_RECORD_PROPERTIES,
    METABOLITE_EQUIVALENCE_EDGES,
    METABOLITE_EXPECTED_CLIQUES,
    METABOLITE_MW_ADJUDICATIONS,
    METABOLITE_RECORD_PROPERTIES,
    METABOLITE_RECORD_SUPPRESSIONS,
    curatable_property_definitions,
    payload_sha256,
    record_property_curation_set_for_model,
    record_property_type_for_model,
    resolve_curation_type,
    validate_operation,
)
from src.core.graph_build_identity import source_build_fingerprint
from src.core.record_property_curations import (
    CURATION_ORIGINAL_FIELD,
    project_record_decisions,
    MISSING_ORIGINAL_MARKER,
    apply_record_property_decision,
    canonical_path as _canonical_record_property_path,
    display_path as _display_record_property_path,
    is_protected_curation_field,
    recalculate_curated_structure_derivatives,
    resolve_parent_and_field,
    schema_for_path,
    validate_value_for_schema,
)
from src.shared.metabolite_generic_structure import (
    GENERIC_STRUCTURE_CLASSIFIER_VERSION,
    generic_structure_evidence as _generic_structure_evidence,
    structure_values as _structure_values_for_generic_check,
    value_has_generic_structure_token as _value_has_generic_structure_token,
)
from src.shared.metabolite_mass_validation import (
    MW_VALIDATION_VERSION,
    assess_mass_profiles,
    cluster_mass_observations,
    compare_mw_review_evidence,
    mass_validation_finding_metadata,
    mw_review_evidence_snapshot,
)
from src.qa_browser.build_provenance import extract_build_inputs
from src.qa_browser.harmonization_stage_maintenance import (
    delete_one_orphan_stage,
    delete_stage_artifacts,
)
from src.qa_browser.metabolite_curation_review import (
    ReviewFilters,
    build_active_curation_review,
)
from src.qa_browser.disease_id_graph import (
    DOWNLOADABLE_FILES,
    REVIEW_DECISION_OPTIONS,
    REVIEW_INTAKE_COLUMNS,
    build_review_queue,
    build_clinical_descendant_graph_payload,
    build_disease_graph_payload,
    build_hierarchy_graph_payload,
    build_provenance_chain,
    bulk_resolve_names,
    compute_dashboard_stats,
    compute_download_preview_counts,
    compute_source_agreement_matrix,
    compute_source_flow_data,
    compute_version_diff,
    export_filtered_download,
    export_filtered_xlsx,
    export_flagged_tsv,
    export_review_intake_template,
    export_search_tsv,
    export_sssom,
    iter_sssom_bytes,
    find_concept_neighbors,
    load_baseline_data,
    load_disease_graph_data,
    load_version_data,
    load_flagged_concepts,
    _normalize_download_mode,
    parse_disease_ids,
    resolve_concept,
    resolve_name_candidates,
    search_concepts,
    stream_filtered_download,
)
from src.qa_browser.cure_entity_resolver import build_cure_entity_resolver_index
from src.qa_browser.curation_cart import (
    add_cart_operation,
    load_cart,
    publish_cart,
    remove_cart_operation,
)
try:
    from src.qa_browser.metabolite_anomer_rule import (
        CHEBI_POLYMER_ROOT_IDS,
        FREE_ANOMER_ALGORITHM_VERSION,
        FREE_ANOMER_RULE_ID,
        evaluate_free_anomer_candidate,
        free_anomeric_oh_atoms,
        normalized_free_anomer_structure,
        summarize_free_anomer_decisions,
    )
    _metabolite_free_anomer_available = True
except ModuleNotFoundError as exc:
    if exc.name != "rdkit":
        raise
    CHEBI_POLYMER_ROOT_IDS = ("CHEBI:60027", "CHEBI:18154")
    FREE_ANOMER_ALGORITHM_VERSION = "free-anomer-v2"
    FREE_ANOMER_RULE_ID = "merge_free_anomeric_forms"
    _metabolite_free_anomer_available = False

    def _missing_free_anomer_dependency(*_args, **_kwargs):
        raise RuntimeError("The free-anomer metabolite rule requires RDKit, but RDKit is not installed.")

    free_anomeric_oh_atoms = _missing_free_anomer_dependency
    normalized_free_anomer_structure = _missing_free_anomer_dependency
    evaluate_free_anomer_candidate = _missing_free_anomer_dependency

    def summarize_free_anomer_decisions(_decisions):
        return {
            "free_anomer_candidate_pair_count": 0,
            "free_anomer_accepted_pair_count": 0,
            "free_anomer_rejected_pair_count": 0,
            "free_anomer_rejection_reason_counts": {"rdkit_unavailable": 1},
        }
from src.qa_browser.ramp_id_graph import set_ramp_diagnosis_file
from src.qa_browser.target_id_graph import (
    TARGET_REVIEW_INTAKE_COLUMNS as TARGET_REVIEW_INTAKE_COLUMNS,
    append_target_review_rows,
    build_batch_review_payload,
    build_pharos_ppi_payload,
    build_target_graph_payload,
    build_target_review_queue,
    bulk_resolve_targets,
    compute_target_stats,
    compute_target_version_diff,
    export_divergences,
    export_entity_ids,
    export_targets,
    fetch_string_ppi,
    invalidate_pharos_tdl_cache,
    load_pharos_tdl_data,
    load_target_graph_data,
    lookup_tdl_by_symbols,
    load_target_version_data,
    mark_rows_resolved_by_registry_ids,
    mark_rows_resolved_by_triage,
    match_targets_to_tdl,
    resolve_target_query,
    search_targets,
    validate_and_build_target_review_rows,
)
from src.qa_browser.variant_id_graph import (
    VARIANT_REVIEW_DECISION_OPTIONS,
    VARIANT_REVIEW_INTAKE_COLUMNS,
    build_variant_graph_payload,
    build_variant_review_queue,
    compute_variant_stats,
    compute_variant_version_diff,
    export_variant_review_intake_template,
    export_variants,
    load_variant_graph_data,
    load_variant_version_data,
    search_variants,
)
from src.qa_browser.drug_id_graph import (
    DRUG_REVIEW_DECISION_OPTIONS,
    DRUG_REVIEW_INTAKE_COLUMNS,
    build_drug_graph_payload,
    build_drug_review_queue,
    compute_drug_stats,
    compute_drug_version_diff,
    export_drug_sssom,
    export_drugs,
    export_drug_review_intake_template,
    find_drug_neighbors,
    iter_drug_sssom_bytes,
    load_drug_graph_data,
    load_drug_version_data,
    resolve_drug,
    search_drugs,
    summarize_drug_source_versions,
    _validate_drug_resolver_index,
)
from src.qa_browser.pathway_id_graph import (
    PATHWAY_REVIEW_DECISION_OPTIONS,
    PATHWAY_REVIEW_INTAKE_COLUMNS,
    build_batch_pathway_review_payload,
    load_pathway_graph_data,
    load_pathway_categories,
    mark_pathway_rows_resolved,
    resolve_pathway_queries,
    export_pathway_resolver_results,
    search_pathways,
    compute_pathway_stats,
    build_pathway_graph_payload,
    build_cross_entity_summary,
    build_pathway_review_queue,
    export_pathways,
    export_pathway_review_intake_template,
    load_pathway_version_data,
    compute_pathway_version_diff,
    DEFAULT_PATHWAY_COLUMNS,
)
from src.qa_browser.drug_resolver import get_ncats_property_catalog, resolve_and_enrich
from src.qa_browser.variant_resolver import (
    resolve_and_enrich as variant_resolve_and_enrich,
    get_myvariant_field_catalog,
)

urllib3.disable_warnings(urllib3.exceptions.InsecureRequestWarning)

BASE_DIR = Path(__file__).parent
TEMPLATES_DIR = BASE_DIR / "templates"
STATIC_DIR = BASE_DIR / "static"
RAMP_ID_QA_ABOUT_DIR = STATIC_DIR / "ramp_id_qa_about"
DISEASE_APP_GRAPH_BUNDLED_DIR = BASE_DIR / "data" / "disease_app_graph"
TARGET_APP_GRAPH_BUNDLED_DIR = BASE_DIR / "data" / "target_app_graph"
VARIANT_APP_GRAPH_BUNDLED_DIR = BASE_DIR / "data" / "variant_app_graph"
DRUG_APP_GRAPH_BUNDLED_DIR = BASE_DIR / "data" / "drug_app_graph"
PATHWAY_APP_GRAPH_BUNDLED_DIR = BASE_DIR / "data" / "pathway_app_graph"


app = FastAPI(title="QA Browser")
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")
app.mount(
    "/ramp-id-qa/about",
    StaticFiles(directory=str(RAMP_ID_QA_ABOUT_DIR), html=True),
    name="ramp_id_qa_about",
)
templates = Jinja2Templates(directory=str(TEMPLATES_DIR))
templates.env.globals["root_path"] = ""
templates.env.globals["style_version"] = hashlib.sha1(
    (STATIC_DIR / "style.css").read_bytes()
).hexdigest()[:12]
templates.env.globals["metabolite_visuals_version"] = hashlib.sha1(
    (STATIC_DIR / "metabolite_harmonization_visuals.js").read_bytes()
).hexdigest()[:12]

# -- Global state set at startup --
_client: Optional[ArangoClient] = None
_credentials: dict = {}
_mysql_credentials: dict = {}
_mysql_sources: dict = {}
_mysql_db_engines: dict = {}
_mysql_inspector_cache: dict = {}   # db_name -> CachableInspector data
_object_storage_credentials: dict = {}
_parquet_storage_credentials: dict = {}
_disease_graph_dir: str = ""
_baseline_graph_dir: str = ""
_target_graph_dir: str = ""
_variant_graph_dir: str = ""
_variant_review_file: str = ""
_drug_graph_dir: str = ""
_drug_resolver_index: str = ""
_verify_drug_resolver_checksum: bool = False
_drug_review_file: str = ""
_pathway_graph_dir: str = ""
_pathway_categories_file: str = ""
_pathway_review_file: str = ""
_pathway_review_lock = threading.Lock()
DRUG_RESOLVER_MAX_QUERIES = 2000
_DRUG_RESOLVER_JOB_TTL_SECONDS = 3600
_drug_resolver_jobs: dict[str, dict[str, Any]] = {}
_drug_resolver_jobs_lock = threading.Lock()
_disease_review_file: str = ""
_disease_review_lock = threading.Lock()
_drug_review_lock = threading.Lock()
_variant_review_lock = threading.Lock()
_cure_entity_resolver_cache: dict = {
    "graph_dir": None,
    "index": None,
    "loaded_at": None,
}
_cure_entity_resolver_lock = threading.Lock()
_HARMONIZATION_PIPELINE_COLLECTION = "HarmonizationPipeline"
_HARMONIZATION_PIPELINE_RUN_COLLECTION = "HarmonizationPipelineRun"
_HARMONIZATION_RUNNER_BOOT_ID = uuid.uuid4().hex
_HARMONIZATION_RECONCILIATION_RETRY_SECONDS = 30
_harmonization_reconciliation_lock = threading.Lock()
_harmonization_reconciliation_last_attempt = 0.0
_harmonization_pipeline_mutation_lock = threading.Lock()
_HARMONIZATION_CURATION_CACHE_TTL_SECONDS = 30
_harmonization_curation_cache_lock = threading.Lock()
_harmonization_curation_refresh_lock = threading.Lock()
_harmonization_curation_cache: dict[str, dict] = {}
_harmonization_curation_cache_generation = 0
_HARMONIZATION_STAGE_COLLECTION = "HarmonizationStage"
_HARMONIZED_METABOLITE_COLLECTION = "HarmonizedMetabolite"
_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION = "HarmonizedMetaboliteMemberEdge"
_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION = "HarmonizationStageEvidenceEdge"
_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION = "HarmonizationStageActiveIdentifierChunk"
_HARMONIZATION_ENGINE_VERSION = "staged-pipeline-v9"
_HARMONIZATION_AQL_BATCH_SIZE = 1000
# Pipeline mutation and stage retirement share one process-wide lock in the
# single-worker QA deployment, so an unreferenced complete stage is immediately
# safe to claim. Delaying it would require a poller or delayed durable job to
# guarantee that a newly orphaned stage is revisited.
_HARMONIZATION_ORPHAN_GRACE_PERIOD = timedelta(0)
_HARMONIZATION_CLEANUP_MAX_RETRIES = 2
_HARMONIZATION_SOURCE_COLLECTIONS = (
    "MetaboliteIdentifier",
    "MetaboliteIdentifierMappingEdge",
    "ChemicalEntity",
    "IsAEdge",
    "MetabolitePathwayEdge",
    "HmdbMetaboliteProteinAssociationEdge",
    "RheaMetaboliteReactionEdge",
    "MetaboliteClassificationEdge",
    "HmdbMetaboliteOntologyEdge",
)
_HMDB_IGNORED_PREFIX_DEFAULTS = [
    "BiGG",
    "BioCyc",
    "DRUGBANK",
    "DRUGBANK.METABOLITE",
    "FoodDB",
    "KNApSAcK",
    "METLIN",
    "MetaGene",
    "NugoWiki",
    "PDB.HET",
    "PhenolExplorer.COMPOUND",
    "PhenolExplorer.METABOLITE",
    "Wikipedia",
]
_WIKIPATHWAYS_IGNORED_PREFIX_DEFAULTS = [
    "PID.PATHWAY",
    "Reactome",
    "LipidBank",
    "PUBCHEM.SUBSTANCE",
    "DRUGBANK",
    "ChEMBL.COMPOUND",
    "KEGG.DRUG",
    "TTD.DRUG",
    "InChIKey",
    "UniProtKB",
    "PharmGKB.DRUG",
]
_WIKIPATHWAYS_IGNORED_SOURCE_FIELD_DEFAULTS = [
    "bdbKeggCompound",
]
_LIPIDMAPS_IGNORED_PREFIX_DEFAULTS = []
_METABOLITE_HARMONIZATION_RULES = [
    {
        "id": "apply_curations",
        "group": "Curation",
        "label": "Apply curations",
        "description": "Apply selected reusable curation types as an overlay without changing the evidence graph.",
        "parameters": [
            {
                "id": "curation_types",
                "label": "Curation types",
                "type": "textarea",
                "default": "\n".join([
                    CHEBI_RECORD_PROPERTIES,
                    METABOLITE_EQUIVALENCE_EDGES,
                    METABOLITE_RECORD_SUPPRESSIONS,
                    METABOLITE_RECORD_PROPERTIES,
                ]),
                "placeholder": "chebi_record_properties\nmetabolite_equivalence_edges\nmetabolite_record_suppressions\nmetabolite_record_properties",
            },
        ],
    },
    {
        "id": "ignore_generic_structure_mismatch",
        "group": "Pruning",
        "label": "Ignore generic/non-generic edges",
        "description": "Drop equivalence edges between explicitly generic and explicitly specific structures; keep unclassified endpoints for validation.",
    },
    {
        "id": "ignore_hmdb_prefixes",
        "group": "Pruning",
        "label": "Ignore HMDB prefixes",
        "description": "Drop HMDB-reported equivalence edges to configured identifier prefixes.",
        "parameters": [
            {
                "id": "prefixes",
                "label": "Prefixes",
                "type": "textarea",
                "default": "\n".join(_HMDB_IGNORED_PREFIX_DEFAULTS),
                "placeholder": "FoodDB\nKNApSAcK\nDRUGBANK",
            },
        ],
    },
    {
        "id": "ignore_wikipathways_prefixes",
        "group": "Pruning",
        "label": "Ignore WikiPathways prefixes",
        "description": "Drop WikiPathways-reported equivalence edges to configured identifier prefixes or BridgeDb source fields.",
        "parameters": [
            {
                "id": "prefixes",
                "label": "Prefixes",
                "type": "textarea",
                "default": "\n".join(_WIKIPATHWAYS_IGNORED_PREFIX_DEFAULTS),
                "placeholder": "PUBCHEM.SUBSTANCE\nKEGG.DRUG\nReactome",
            },
            {
                "id": "source_fields",
                "label": "Source fields",
                "type": "textarea",
                "default": "\n".join(_WIKIPATHWAYS_IGNORED_SOURCE_FIELD_DEFAULTS),
                "placeholder": "bdbKeggCompound",
            },
        ],
    },
    {
        "id": "ignore_lipidmaps_prefixes",
        "group": "Pruning",
        "label": "Ignore LipidMaps prefixes",
        "description": "Drop LipidMaps-reported equivalence edges to configured identifier prefixes.",
        "parameters": [
            {
                "id": "prefixes",
                "label": "Prefixes",
                "type": "textarea",
                "default": "\n".join(_LIPIDMAPS_IGNORED_PREFIX_DEFAULTS),
                "placeholder": "Optional, for experiments only",
            },
        ],
    },
    {
        "id": "force_expected_clique_assertions",
        "group": "Cleanup",
        "label": "Give up: Make Assertions True",
        "description": "Add explicitly synthetic, curation-provenanced edges so every active expected-clique assertion merges. Use only as a visible fallback before cleanup.",
    },
    {
        "id": "remove_enrichment_only_metabolites",
        "group": "Cleanup",
        "label": "Cleanup: Remove Secondary-Source-Only Metabolites",
        "description": "Remove harmonized metabolites whose only sourcing is from secondary databases (RefMet, PubChem, ChEBI, LipidMaps) and have no meaningful biological connections. Use this last.",
    },
    {
        "id": "merge_shared_inchikey_prefix",
        "group": "Merging",
        "label": "Merge by InChIKey prefix",
        "description": "Use the source-reported InChIKey when available; otherwise use one calculated from the stored structure. Merge identifiers with the same first block.",
    },
    {
        "id": "merge_shared_inchikey_duplex",
        "group": "Merging",
        "label": "Merge by InChIKey duplex",
        "description": "Use the source-reported InChIKey when available; otherwise use one calculated from the stored structure. Merge identifiers with the same first two blocks.",
    },
    {
        "id": "merge_inchikey_by_mw_cutoff",
        "group": "Merging",
        "label": "Merge by InChIKey with MW cutoff",
        "description": "Use the source-reported InChIKey when available; otherwise use one calculated from the stored structure. Use the first two blocks below the cutoff and the first block at or above it.",
        "parameters": [
            {
                "id": "mw_cutoff",
                "label": "MW cutoff",
                "type": "number",
                "default": 500,
                "min": 0,
                "step": 1,
                "unit": "Da",
            },
        ],
    },
    {
        "id": FREE_ANOMER_RULE_ID,
        "group": "Chemistry-based merging",
        "label": "Merge Free Anomeric Forms",
        "description": "Merge direct ChEBI child/parent forms that differ only at one free hemiacetal or hemiketal hydroxyl; preserve all other stereochemistry and exclude polymers.",
    },
]
_metabolite_snapshot_jobs: dict = {}
_metabolite_snapshot_jobs_lock = threading.Lock()
_METABOLITE_MW_SPREAD_WARNING_THRESHOLD = 0.10
_METABOLITE_MW_SPREAD_WARNING_LIMIT = 50
_METABOLITE_DENYLIST_STILL_MERGED_LIMIT = 50
_METABOLITE_GENERIC_STRUCTURE_WARNING_LIMIT = 50
_METABOLITE_COMPARE_MAX_IDS = 500
_METABOLITE_INTERACTIVE_CLIQUE_NODE_LIMIT = 100
_METABOLITE_INTERACTIVE_CLIQUE_EDGE_LIMIT = 300
_demo_queries_enabled = os.getenv("QA_BROWSER_ENABLE_POUNCE_DEMOS", "").lower() in {
    "1", "true", "yes", "on"
}


def _root_path(request: Optional[Request] = None) -> str:
    if request is not None:
        root_path = request.scope.get("root_path")
        if root_path:
            return str(root_path).rstrip("/")
    return str(templates.env.globals.get("root_path") or "").rstrip("/")


def _app_path(path: str, request: Optional[Request] = None) -> str:
    normalized_path = path if path.startswith("/") else f"/{path}"
    return f"{_root_path(request)}{normalized_path}"


def _redirect_to(path: str, *, request: Optional[Request] = None, status_code: int = 303) -> RedirectResponse:
    return RedirectResponse(url=_app_path(path, request), status_code=status_code)


def _is_hidden_qa_browser_collection(collection_name: str) -> bool:
    return collection_name.startswith("_")


def _is_qa_visible_edge_collection(collection_name: str) -> bool:
    return collection_name in {
        _HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION,
        _HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION,
    }


def _harmonization_code_fingerprint() -> str:
    """Short hash of the functions that most affect harmonization results.

    Stored in each stage doc so mismatches across a pipeline's stages are
    immediately detectable (different fingerprint = different server code).
    """
    code = "\n".join(
        inspect.getsource(fn)
        for fn in [
            _value_has_generic_structure_token,
            _load_generic_structure_classifications,
            _load_generic_structure_ids,
            _wikipathways_xref_only_ids_from_rows,
            _filter_identifier_support_for_rules,
            _active_metabolite_identifier_mapping_edges_for_rules,
            _effective_inchi_key_matches,
            _iter_metabolite_identifier_inchi_key_matches,
            _build_harmonized_groups,
            free_anomeric_oh_atoms,
            normalized_free_anomer_structure,
            evaluate_free_anomer_candidate,
            _free_anomer_rule_edges,
        ]
    )
    return hashlib.sha1(code.encode()).hexdigest()[:12]


def _metabolite_generic_structure_classification(node: dict, curation_state: Optional[dict] = None) -> dict:
    structures = []
    structure_records = (
        node.get("generic_structure_evidence")
        if node.get("generic_structure_evidence") is not None
        else (node.get("chem_props") or [])
    )
    for prop in structure_records:
        structures.append({
            "source": prop.get("source"),
            "source_id": prop.get("source_id"),
            "smiles": prop.get("smiles"),
            "formula": prop.get("formula"),
            "inchi": prop.get("inchi"),
        })
    chemical_entity = node.get("chemical_entity") or {}
    if any(chemical_entity.get(field) for field in ("smiles", "formula", "inchi")):
        structures.append({
            "source": "ChEBI",
            "source_id": chemical_entity.get("id"),
            "smiles": chemical_entity.get("smiles"),
            "formula": chemical_entity.get("formula"),
            "inchi": chemical_entity.get("inchi"),
        })
    evidence = _generic_structure_evidence(structures)
    graph_values = node.get("editable_properties") or {}
    detected = (
        graph_values.get("is_generic_structure")
        if "is_generic_structure" in graph_values
        else node.get("is_generic_structure")
    )
    state_available = (curation_state or {}).get(
        "record_property_state_available", True
    )
    state_error = (curation_state or {}).get("record_property_state_error")
    decision = (
        ((curation_state or {}).get("record_property_decisions") or {})
        .get(node.get("id"), {})
        .get(_canonical_record_property_path(["is_generic_structure"]))
        or {}
    )
    has_override = decision.get("mode") == "set"
    published_override = decision.get("value") if has_override else None
    effective = published_override if has_override else detected
    provenance = decision or None
    return {
        "detected": detected,
        "published_override": published_override,
        "has_published_override": has_override,
        "effective": effective,
        "evidence": evidence,
        "reason": (
            f"Generic marker in {evidence[0]['field']} from {evidence[0].get('source') or 'structure evidence'}"
            if evidence
            else (
                "Persisted graph classification is specific"
                if detected is False
                else "No persisted generic-structure classification"
            )
        ),
        "provenance": provenance,
        "state_available": state_available,
        "state_error": state_error,
    }


def _metabolite_harmonization_participation(
    identifier: str,
    curation_state: Optional[dict] = None,
) -> dict:
    action = ((curation_state or {}).get("record_states") or {}).get(identifier)
    decision = ((curation_state or {}).get("record_decisions") or {}).get(identifier)
    return {
        "status": "suppressed" if action == "suppress_record" else "active",
        "is_suppressed": action == "suppress_record",
        "published_action": action,
        "decision": decision,
        "state_available": (curation_state or {}).get("suppression_state_available", True),
        "state_error": (curation_state or {}).get("suppression_state_error"),
    }


_INCHI_KEY_PREFIX_RE = re.compile(r"^[A-Z]{14}$")
_INCHI_KEY_DUPLEX_RE = re.compile(r"^[A-Z]{14}-[A-Z]{10}$")
_INCHI_KEY_FULL_RE = re.compile(r"^([A-Z]{14})-([A-Z]{10})-[A-Z]$")


def _normalize_inchi_key_prefix(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    prefix = str(value).strip().upper()
    if not prefix:
        return None
    prefix = prefix.split("-", 1)[0]
    if not _INCHI_KEY_PREFIX_RE.match(prefix):
        return None
    return prefix


def _normalize_inchi_key_duplex(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip().upper()
    if not text:
        return None
    match = _INCHI_KEY_FULL_RE.match(text)
    if match:
        return f"{match.group(1)}-{match.group(2)}"
    parts = text.split("-")
    if len(parts) >= 2:
        duplex = "-".join(parts[:2])
        if _INCHI_KEY_DUPLEX_RE.match(duplex):
            return duplex
    return None


def _parse_optional_float(value) -> Optional[float]:
    if value is None:
        return None
    try:
        return float(str(value).strip())
    except (TypeError, ValueError):
        return None


def _identifier_prefix(identifier: Optional[str]) -> Optional[str]:
    if not identifier or ":" not in str(identifier):
        return None
    prefix = str(identifier).split(":", 1)[0].strip()
    return prefix or None


def _normalize_prefix_list(value) -> List[str]:
    if value is None:
        return []
    if isinstance(value, (list, tuple, set)):
        raw_items = [str(item) for item in value]
    else:
        raw_items = re.split(r"[\s,;]+", str(value))
    seen = set()
    prefixes = []
    for item in raw_items:
        prefix = item.strip()
        if not prefix:
            continue
        prefix_key = prefix.lower()
        if prefix_key in seen:
            continue
        seen.add(prefix_key)
        prefixes.append(prefix)
    return prefixes


def _normalize_ramp_denylist_identifier(value: Optional[str]) -> Optional[str]:
    if value is None:
        return None
    text = str(value).strip()
    if not text or ":" not in text:
        return None
    prefix, local_id = text.split(":", 1)
    prefix_key = prefix.strip().lower()
    local_id = local_id.strip()
    if not local_id:
        return None
    prefix_map = {
        "hmdb": "HMDB",
        "chebi": "CHEBI",
        "kegg": "KEGG.COMPOUND",
        "kegg.compound": "KEGG.COMPOUND",
        "pubchem": "PUBCHEM.COMPOUND",
        "pubchem.compound": "PUBCHEM.COMPOUND",
        "chemspider": "ChemSpider",
        "refmet": "REFMET",
        "lipidmaps": "LIPIDMAPS",
        "cas": "CAS",
        "wikidata": "Wikidata",
    }
    normalized_prefix = prefix_map.get(prefix_key, prefix.strip())
    return f"{normalized_prefix}:{local_id}"


def _metabolite_edge_decisions_from_curation_batch(
    batch: dict,
) -> List[tuple[tuple[str, str], str]]:
    if batch.get("format_version") != 1:
        raise ValueError(f"unsupported curation format_version: {batch.get('format_version')!r}")
    if batch.get("graph") != "metabolite_harmonization":
        return []

    decisions = []
    for operation in batch.get("operations") or []:
        if (
            operation.get("action") not in {"remove_edge", "retain_edge"}
            or operation.get("edge_type") != "MetaboliteIdentifierMappingEdge"
        ):
            continue
        if operation.get("symmetric") is not True:
            continue
        left = _normalize_ramp_denylist_identifier(operation.get("start_id"))
        right = _normalize_ramp_denylist_identifier(operation.get("end_id"))
        if left and right and left != right:
            decisions.append((tuple(sorted((left, right))), operation["action"]))
    return decisions


def _metabolite_edge_removal_pairs_from_curation_batch(batch: dict) -> set[tuple[str, str]]:
    pair_states = {}
    for pair, action in _metabolite_edge_decisions_from_curation_batch(batch):
        pair_states[pair] = action
    return {pair for pair, action in pair_states.items() if action == "remove_edge"}


def _normalize_expected_clique_member_ids(value: Any) -> List[str]:
    raw_ids = _parse_metabolite_identifier_query(
        " ".join(str(item) for item in value) if isinstance(value, (list, tuple, set)) else str(value or "")
    )
    normalized = []
    seen = set()
    for raw_id in raw_ids:
        identifier = _normalize_ramp_denylist_identifier(raw_id)
        if not identifier:
            raise HTTPException(status_code=400, detail=f"Invalid prefixed metabolite identifier: {raw_id}")
        if identifier not in seen:
            seen.add(identifier)
            normalized.append(identifier)
    normalized.sort()
    if len(normalized) < 2:
        raise HTTPException(status_code=400, detail="An expected-clique assertion requires at least two identifiers.")
    if len(normalized) > _METABOLITE_COMPARE_MAX_IDS:
        raise HTTPException(
            status_code=400,
            detail=f"An expected-clique assertion may contain at most {_METABOLITE_COMPARE_MAX_IDS} identifiers.",
        )
    return normalized


def _metabolite_expected_clique_assertion_operation(
    member_ids: Any,
    name: str,
    rationale: str = "",
    missing_member_ids_at_add_time: Optional[List[str]] = None,
) -> dict:
    normalized_ids = _normalize_expected_clique_member_ids(member_ids)
    clean_name = str(name or "").strip()
    if not clean_name:
        raise HTTPException(status_code=400, detail="Assertion name is required.")
    assertion_id = "same-clique-" + hashlib.sha256(
        json.dumps(normalized_ids, separators=(",", ":")).encode("utf-8")
    ).hexdigest()[:24]
    operation = {
        "action": "assert_same_clique",
        "assertion_id": assertion_id,
        "assertion_type": "expected_same_clique",
        "member_ids": normalized_ids,
        "name": clean_name,
    }
    clean_rationale = str(rationale or "").strip()
    if clean_rationale:
        operation["rationale"] = clean_rationale
    missing = sorted(set(missing_member_ids_at_add_time or []))
    if missing:
        operation["missing_member_ids_at_add_time"] = missing
    return operation


def _metabolite_retire_assertion_operation(assertion_id: str, note: str = "") -> dict:
    clean_id = str(assertion_id or "").strip()
    if not re.fullmatch(r"same-clique-[0-9a-f]{24}", clean_id):
        raise HTTPException(status_code=400, detail="A valid expected-clique assertion ID is required.")
    operation = {"action": "retire_assertion", "assertion_id": clean_id}
    clean_note = str(note or "").strip()
    if clean_note:
        operation["note"] = clean_note
    return operation


_MW_ADJUDICATION_REASONS = (
    ("protonation_or_charge_state", "Protonation or charge state"),
    ("salt_or_counterion", "Salt or counterion"),
    ("hydrate_or_solvate", "Hydrate or solvate"),
    ("multicomponent_form", "Multicomponent form"),
    ("source_mass_convention", "Source mass convention"),
    ("polymer_or_repeat_unit", "Polymer or repeat unit"),
    ("other_reviewed_explanation", "Other reviewed explanation"),
)


def _metabolite_mw_adjudication_operation(
    action: str,
    stage_key: str,
    finding_id: str,
    reason: str,
    supporting_ids: Any,
    note: str,
) -> dict:
    if action not in {"accept_mw_discrepancy", "reopen_mw_discrepancy"}:
        raise HTTPException(status_code=400, detail="Unsupported MW adjudication action.")
    stage = _get_harmonization_stage(str(stage_key or "").strip())
    validation = _with_mw_finding_metadata(
        _prioritize_kegg_validation_samples(
            (stage.get("validation") or {}).get("mw_spread") or {}
        )
    )
    finding = next((
        item for item in validation.get("warnings") or []
        if item.get("finding_id") == str(finding_id or "").strip()
    ), None)
    if finding is None:
        raise HTTPException(
            status_code=409,
            detail="This MW finding is not present in the requested saved stage. Reload it before curating.",
        )
    member_ids = finding.get("member_ids") or []
    anchor_id = finding["anchor_id"]
    parsed_supporting_ids = _parse_metabolite_identifier_query(
        " ".join(str(item) for item in supporting_ids)
        if isinstance(supporting_ids, list)
        else str(supporting_ids or "")
    )
    unknown_supporting = sorted(set(parsed_supporting_ids) - set(member_ids))
    if unknown_supporting:
        raise HTTPException(
            status_code=400,
            detail="Supporting IDs are not members of this finding: " + ", ".join(unknown_supporting),
        )
    operation = {
        "action": action,
        "target": {
            "kind": "validation_finding",
            "curation_set": "metabolite_harmonization",
            "check": "mw_spread",
            "finding_id": finding["finding_id"],
            "anchor_id": anchor_id,
        },
        "observed_evidence_fingerprint": finding["evidence_fingerprint"],
        "observed_member_ids": member_ids,
        "reviewed_stage_key": str(stage_key or "").strip(),
        "observed_evidence_snapshot": mw_review_evidence_snapshot(
            finding,
            validator_version=validation.get("version"),
            threshold=validation.get("threshold"),
        ),
        "note": str(note or "").strip(),
    }
    if action == "accept_mw_discrepancy":
        valid_reasons = {value for value, _label in _MW_ADJUDICATION_REASONS}
        clean_reason = str(reason or "").strip()
        if clean_reason not in valid_reasons:
            raise HTTPException(status_code=400, detail="Choose a valid MW acceptance reason.")
        operation["reason"] = clean_reason
        operation["supporting_ids"] = sorted(set(parsed_supporting_ids))
    try:
        validate_operation(METABOLITE_MW_ADJUDICATIONS, operation)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return operation


def _missing_metabolite_identifier_ids(member_ids: List[str]) -> List[str]:
    db = get_db("metabolite_harmonization")
    found_ids = set(db.aql.execute(
        """
        FOR metabolite IN MetaboliteIdentifier
          FILTER metabolite.id IN @member_ids
          RETURN metabolite.id
        """,
        bind_vars={"member_ids": member_ids},
        max_runtime=120,
    ))
    return sorted(set(member_ids) - found_ids)


def _metabolite_assertion_decisions_from_curation_batch(batch: dict) -> List[tuple[str, Optional[dict]]]:
    if batch.get("format_version") != 1:
        raise ValueError(f"unsupported curation format_version: {batch.get('format_version')!r}")
    if batch.get("graph") != "metabolite_harmonization":
        return []
    decisions = []
    for operation in batch.get("operations") or []:
        action = operation.get("action")
        assertion_id = str(operation.get("assertion_id") or "").strip()
        if action == "retire_assertion" and assertion_id:
            decisions.append((assertion_id, None))
            continue
        if action != "assert_same_clique" or not assertion_id:
            continue
        member_ids = _normalize_expected_clique_member_ids(operation.get("member_ids") or [])
        expected_id = _metabolite_expected_clique_assertion_operation(
            member_ids,
            operation.get("name"),
            operation.get("rationale") or "",
            operation.get("missing_member_ids_at_add_time") or [],
        )["assertion_id"]
        if assertion_id != expected_id:
            raise ValueError(f"Expected-clique assertion ID does not match its members: {assertion_id}")
        decisions.append((assertion_id, dict(operation)))
    return decisions


def _metabolite_curation_publication_time(batch: dict, key: str) -> datetime:
    value = str(batch.get("published_at") or batch.get("created_at") or "").strip()
    if not value:
        return datetime.min.replace(tzinfo=timezone.utc)
    try:
        parsed = datetime.fromisoformat(value.replace("Z", "+00:00"))
    except ValueError as exc:
        raise ValueError(f"Invalid curation publication timestamp in {key}: {value!r}") from exc
    if parsed.tzinfo is None:
        parsed = parsed.replace(tzinfo=timezone.utc)
    return parsed.astimezone(timezone.utc)


def _record_property_snapshot_state(
    snapshot,
    *,
    curation_set: str,
    model_type: str,
) -> dict:
    decisions = sorted([
        decision
        for decision in snapshot.record_property_decisions_for_set(
            curation_set
        )
        if decision.target.get("model_type") == model_type
    ], key=lambda decision: (
        decision.target.get("id") or "",
        _canonical_record_property_path(decision.path),
    ))
    fingerprint = payload_sha256([
        {
            "target": decision.target,
            "path": decision.path,
            "mode": decision.mode,
            "value": decision.value,
            "observed_value": decision.observed_value,
            "observed_exists": decision.observed_exists,
            "batch_id": decision.batch_id,
        }
        for decision in decisions
    ])
    selected_batch_ids = {decision.batch_id for decision in decisions}
    metadata = {
        **snapshot.metadata(),
        "selected_curation_set": curation_set,
        "selected_model_type": model_type,
        "selected_batch_ids": [
            batch_id for batch_id in snapshot.batch_ids
            if batch_id in selected_batch_ids
        ],
        "selected_active_decision_count": len(decisions),
        "selected_resolved_operation_fingerprint": fingerprint,
    }
    return {
        "decisions": decisions,
        "batch_ids": metadata["selected_batch_ids"],
        "fingerprint": fingerprint,
        "metadata": metadata,
    }


def _metabolite_record_property_snapshot_state(snapshot) -> dict:
    return _record_property_snapshot_state(
        snapshot,
        curation_set="metabolite_harmonization",
        model_type="MetaboliteIdentifier",
    )


def _invalidate_harmonization_curation_cache() -> None:
    global _harmonization_curation_cache_generation
    with _harmonization_curation_cache_lock:
        _harmonization_curation_cache.clear()
        _harmonization_curation_cache_generation += 1


def _resolve_metabolite_curations(
    storage=None,
    prefix: Optional[str] = None,
) -> dict:
    """Resolve the typed metabolite curation streams.

    The returned state includes graph-specific decisions, general record-property
    overlays, and the immutable metadata for every snapshot. ``prefix`` is
    rejected so runtime code cannot accidentally dual-read v1.
    """
    if prefix is not None:
        raise ValueError("Custom graph-scoped curation prefixes are no longer supported; use v2 manifests")
    if storage is None:
        if not _object_storage_credentials:
            raise RuntimeError("Object-storage credentials are required to load metabolite curations")
        storage = _storage_from_credentials(_object_storage_credentials, use_internal_url=False)
    curation_types = (
        METABOLITE_EQUIVALENCE_EDGES,
        METABOLITE_EXPECTED_CLIQUES,
        METABOLITE_RECORD_SUPPRESSIONS,
        METABOLITE_RECORD_PROPERTIES,
        CHEBI_RECORD_PROPERTIES,
    )
    # Each immutable curation stream has its own manifest and batch chain. Read
    # the independent streams concurrently so a cold workbench request pays one
    # object-storage round trip sequence rather than five in series.
    with ThreadPoolExecutor(max_workers=len(curation_types)) as executor:
        futures = {
            curation_type: executor.submit(
                resolve_curation_type,
                storage,
                curation_type,
                allow_missing=True,
            )
            for curation_type in curation_types
        }
        snapshots = {
            curation_type: future.result()
            for curation_type, future in futures.items()
        }
    edge_snapshot = snapshots[METABOLITE_EQUIVALENCE_EDGES]
    assertion_snapshot = snapshots[METABOLITE_EXPECTED_CLIQUES]
    suppression_snapshot = snapshots[METABOLITE_RECORD_SUPPRESSIONS]
    record_property_snapshot = snapshots[METABOLITE_RECORD_PROPERTIES]
    chebi_record_property_snapshot = snapshots[CHEBI_RECORD_PROPERTIES]
    record_property_state = _metabolite_record_property_snapshot_state(
        record_property_snapshot
    )
    chebi_record_property_state = _record_property_snapshot_state(
        chebi_record_property_snapshot,
        curation_set="chebi",
        model_type="ChemicalEntity",
    )
    record_property_decisions = {}
    for decision in record_property_state["decisions"]:
        record_property_decisions.setdefault(decision.target["id"], {})[
            _canonical_record_property_path(decision.path)
        ] = {
            "mode": decision.mode,
            "value": decision.value,
            "curation_batch_id": decision.batch_id,
            "published_at": decision.published_at,
            "published_by": decision.published_by,
        }

    pair_states = {}
    pair_decisions = {}
    for resolved in edge_snapshot.active_operations:
        operation = resolved.operation
        left = _normalize_ramp_denylist_identifier(operation.get("start_id"))
        right = _normalize_ramp_denylist_identifier(operation.get("end_id"))
        if left and right and left != right:
            pair = tuple(sorted((left, right)))
            pair_states[pair] = operation["action"]
            pair_decisions[pair] = {
                "action": operation["action"],
                "batch_id": getattr(resolved, "batch_id", None),
                "published_at": getattr(resolved, "published_at", None),
                "published_by": getattr(resolved, "published_by", None),
                "note": operation.get("note"),
            }
    pairs = {pair for pair, action in pair_states.items() if action == "remove_edge"}

    assertions = []
    for resolved in assertion_snapshot.active_operations:
        operation = resolved.operation
        if operation["action"] != "assert_same_clique":
            continue
        assertions.append({
            **operation,
            "curation_batch_id": resolved.batch_id,
            "published_at": resolved.published_at,
            "published_by": resolved.published_by,
        })
    assertions.sort(key=lambda item: (item.get("name", "").casefold(), item["assertion_id"]))

    record_states = {}
    record_decisions = {}
    for resolved in suppression_snapshot.active_operations:
        operation = resolved.operation
        target_id = _normalize_ramp_denylist_identifier(
            (operation.get("target") or {}).get("id")
        )
        if not target_id:
            continue
        record_states[target_id] = operation["action"]
        record_decisions[target_id] = {
            "action": operation["action"],
            "batch_id": resolved.batch_id,
            "published_at": resolved.published_at,
            "published_by": resolved.published_by,
            "note": operation.get("note"),
        }
    suppressed_identifier_ids = {
        identifier
        for identifier, action in record_states.items()
        if action == "suppress_record"
    }

    apply_fingerprint = payload_sha256({
        METABOLITE_EQUIVALENCE_EDGES: edge_snapshot.fingerprint,
        METABOLITE_RECORD_SUPPRESSIONS: suppression_snapshot.fingerprint,
        METABOLITE_RECORD_PROPERTIES: record_property_state["fingerprint"],
        CHEBI_RECORD_PROPERTIES: chebi_record_property_state["fingerprint"],
    })

    result = {
        "pairs": pairs,
        "pair_states": pair_states,
        "pair_decisions": pair_decisions,
        "batch_ids": edge_snapshot.batch_ids,
        "fingerprint": apply_fingerprint,
        "edge_fingerprint": edge_snapshot.fingerprint,
        "assertions": assertions,
        "assertion_batch_ids": assertion_snapshot.batch_ids,
        "assertion_fingerprint": assertion_snapshot.fingerprint,
        "suppressed_identifier_ids": suppressed_identifier_ids,
        "record_states": record_states,
        "record_decisions": record_decisions,
        "suppression_batch_ids": suppression_snapshot.batch_ids,
        "suppression_fingerprint": suppression_snapshot.fingerprint,
        "record_property_snapshot": record_property_snapshot,
        "record_property_decisions": record_property_decisions,
        "record_property_batch_ids": record_property_state["batch_ids"],
        "record_property_fingerprint": record_property_state["fingerprint"],
        "record_property_state_available": True,
        "record_property_state_error": None,
        "chebi_record_property_snapshot": chebi_record_property_snapshot,
        "chebi_record_property_batch_ids": chebi_record_property_state["batch_ids"],
        "chebi_record_property_fingerprint": chebi_record_property_state["fingerprint"],
        "snapshots": {
            METABOLITE_EQUIVALENCE_EDGES: edge_snapshot.metadata(),
            METABOLITE_EXPECTED_CLIQUES: assertion_snapshot.metadata(),
            METABOLITE_RECORD_SUPPRESSIONS: suppression_snapshot.metadata(),
            METABOLITE_RECORD_PROPERTIES: record_property_state["metadata"],
            CHEBI_RECORD_PROPERTIES: chebi_record_property_state["metadata"],
        },
        "suppression_state_available": True,
        "suppression_state_error": None,
        "prefix": f"s3://{storage.bucket}/curations/v2/",
    }
    return result


def _load_metabolite_curations(
    storage=None,
    prefix: Optional[str] = None,
    force_refresh: bool = False,
) -> dict:
    if storage is not None:
        return _resolve_metabolite_curations(storage=storage, prefix=prefix)
    with _harmonization_curation_cache_lock:
        cache_generation = _harmonization_curation_cache_generation
    if not force_refresh:
        with _harmonization_curation_cache_lock:
            cached = _harmonization_curation_cache.get("metabolite_curations")
            if cached and cached["expires_at"] > time.monotonic():
                return cached["value"]
        if cached:
            # A workbench read may use the last successfully resolved immutable
            # snapshot while one background thread refreshes it. This prevents
            # the polling UI from stalling every time the short TTL expires.
            if _harmonization_curation_refresh_lock.acquire(blocking=False):
                def refresh_cache() -> None:
                    try:
                        result = _resolve_metabolite_curations(prefix=prefix)
                        with _harmonization_curation_cache_lock:
                            if (
                                cache_generation
                                == _harmonization_curation_cache_generation
                            ):
                                _harmonization_curation_cache["metabolite_curations"] = {
                                    "value": result,
                                    "expires_at": (
                                        time.monotonic()
                                        + _HARMONIZATION_CURATION_CACHE_TTL_SECONDS
                                    ),
                                }
                    except Exception as exc:
                        print(f"Warning: could not refresh metabolite curations: {exc}")
                    finally:
                        _harmonization_curation_refresh_lock.release()

                threading.Thread(
                    target=refresh_cache,
                    name="metabolite-curation-cache-refresh",
                    daemon=True,
                ).start()
            return cached["value"]
    # Only one request refreshes the remote immutable snapshots at a time.
    with _harmonization_curation_refresh_lock:
        if not force_refresh:
            with _harmonization_curation_cache_lock:
                cached = _harmonization_curation_cache.get("metabolite_curations")
                if cached and cached["expires_at"] > time.monotonic():
                    return cached["value"]
        result = _resolve_metabolite_curations(prefix=prefix)
        with _harmonization_curation_cache_lock:
            if cache_generation == _harmonization_curation_cache_generation:
                _harmonization_curation_cache["metabolite_curations"] = {
                    "value": result,
                    "expires_at": (
                        time.monotonic()
                        + _HARMONIZATION_CURATION_CACHE_TTL_SECONDS
                    ),
                }
        return result


def _load_metabolite_edge_removal_curations(
    storage=None,
    prefix: Optional[str] = None,
    force_refresh: bool = False,
) -> dict:
    """Compatibility wrapper for callers using the historical private name."""
    return _load_metabolite_curations(
        storage=storage,
        prefix=prefix,
        force_refresh=force_refresh,
    )


def _curation_cart_storage():
    if not _object_storage_credentials:
        raise HTTPException(
            status_code=503,
            detail="Object-storage credentials are required to use curation carts.",
        )
    return _storage_from_credentials(_object_storage_credentials, use_internal_url=False)


def _load_record_property_snapshot(model_type: str):
    curation_type = record_property_type_for_model(model_type)
    if not curation_type:
        raise ValueError(
            f"No record-property curation stream is registered for {model_type}"
        )
    return resolve_curation_type(
        _curation_cart_storage(), curation_type, allow_missing=True
    )


def _record_property_observation(document: dict, path: list) -> tuple[bool, Any]:
    parent, field_name = resolve_parent_and_field(document, path)
    originals = parent.get(CURATION_ORIGINAL_FIELD) or {}
    if field_name in originals:
        if originals[field_name] == MISSING_ORIGINAL_MARKER:
            return False, None
        return True, originals[field_name]
    return field_name in parent, parent.get(field_name)


def _record_property_baseline(document: dict, path: list) -> Any:
    return _record_property_observation(document, path)[1]


def _record_property_fields(document: dict, schema_entry: dict, snapshot) -> list[dict]:
    """Build graph-neutral editors from the persisted collection schema."""
    schema_fields = schema_entry.get("fields") or {}
    published = {
        _canonical_record_property_path(decision.path): decision
        for decision in snapshot.record_property_decisions_for_set(document["_curation_set"])
        if decision.target.get("model_type") == document["_curation_model_type"]
        and decision.target.get("id") == document.get("id")
    }
    fields = []
    model_type = document["_curation_model_type"]
    top_level_allowlist = (
        {
            definition.property_name
            for definition in curatable_property_definitions(model_type)
        }
        if model_type == "ChemicalEntity" else None
    )

    def add_field(path, descriptor):
        baseline_present, baseline = _record_property_observation(document, path)
        decision = published.get(_canonical_record_property_path(path))
        effective = (
            decision.value
            if decision is not None and decision.mode == "set"
            else baseline
        )
        schema_type = descriptor.get("type") if isinstance(descriptor, dict) else descriptor
        editor_type = schema_type if schema_type in {"bool", "int", "float", "str"} else "json"
        fields.append({
            "name": _display_record_property_path(path),
            "label": (
                "Generic structure"
                if path == ["is_generic_structure"]
                else str(path[-1]).replace("_", " ").title()
            ),
            "path": path,
            "path_json": json.dumps(path, separators=(",", ":")),
            "editor_type": editor_type,
            "baseline": baseline,
            "baseline_present": baseline_present,
            "baseline_json": json.dumps(baseline, indent=2, default=str),
            "published": decision.value if decision is not None and decision.mode == "set" else None,
            "has_published_override": decision is not None and decision.mode == "set",
            "effective": effective,
            "effective_json": json.dumps(effective, indent=2, default=str),
            "batch_id": decision.batch_id if decision else None,
            "published_at": decision.published_at if decision else None,
            "published_by": decision.published_by if decision else None,
            "rationale": (decision.source_operation.get("note") or "") if decision else "",
        })

    def unique_selector(items: list[dict]) -> Optional[tuple[str, ...]]:
        candidates = (
            ("source", "source_id"),
            ("chem_data_source", "chem_source_id"),
            ("id",),
            ("value", "source"),
        )
        for keys in candidates:
            if not all(all(key in item for key in keys) for item in items):
                continue
            values = [tuple(item.get(key) for key in keys) for item in items]
            if len(values) == len(set(values)):
                return keys
        return None

    def visit(value, descriptor, path):
        schema_type = descriptor.get("type") if isinstance(descriptor, dict) else descriptor
        if schema_type == "object" and isinstance(value, dict):
            for child_name, child_descriptor in (descriptor.get("fields") or {}).items():
                if is_protected_curation_field(child_name) or child_name not in value:
                    continue
                visit(value[child_name], child_descriptor, [*path, child_name])
            return
        if (
            schema_type == "list"
            and isinstance(value, list)
            and value
            and descriptor.get("item_type") == "object"
            and all(isinstance(item, dict) for item in value)
        ):
            selector_keys = unique_selector(value)
            if selector_keys:
                for item in value:
                    selector = {key: item.get(key) for key in selector_keys}
                    selected_path = [*path, {"match": selector}]
                    for child_name, child_descriptor in (descriptor.get("fields") or {}).items():
                        if is_protected_curation_field(child_name) or child_name not in item:
                            continue
                        visit(item[child_name], child_descriptor, [*selected_path, child_name])
                return
        try:
            schema_for_path(schema_fields, path)
        except ValueError:
            return
        add_field(path, descriptor)

    for field_name, descriptor in schema_fields.items():
        if is_protected_curation_field(field_name):
            continue
        if top_level_allowlist is not None and field_name not in top_level_allowlist:
            continue
        visit(document.get(field_name), descriptor, [field_name])
    return fields


def _build_record_property_operation(
    document: dict,
    *,
    curation_set: str,
    model_type: str,
    schema_fields: dict,
    raw_decisions: list[dict],
    note: str,
) -> tuple[str, dict]:
    curation_type = record_property_type_for_model(model_type)
    if not curation_type:
        raise ValueError(
            f"No record-property curation stream is registered for {model_type}"
        )
    decisions = []
    for raw_decision in raw_decisions:
        path = raw_decision.get("path") if isinstance(raw_decision, dict) else None
        mode = str((raw_decision or {}).get("mode") or "set")
        descriptor = schema_for_path(schema_fields, path)
        observed_exists, baseline = _record_property_observation(document, path)
        decision = {"path": path, "mode": mode}
        if mode == "set":
            value = raw_decision.get("value")
            validate_value_for_schema(value, descriptor, path)
            decision.update({
                "value": value,
                "observed_value": baseline,
                "observed_exists": observed_exists,
            })
        decisions.append(decision)
    operation = {
        "action": "set_properties",
        "target": {
            "kind": "node",
            "curation_set": curation_set,
            "model_type": model_type,
            "id": document["id"],
        },
        "decisions": decisions,
        "note": str(note or "").strip(),
    }
    validate_operation(curation_type, operation)
    return curation_type, operation


_INLINE_CURATION_MARKER = "_qa_inline_curation"


def _inline_curation_display_value(value: Any, path: Optional[list] = None) -> str:
    if path == ["is_generic_structure"]:
        if value is True:
            return "Generic"
        if value is False:
            return "Specific / non-generic"
        return "Unknown / unclassified"
    if value is None:
        return "null"
    if value is True:
        return "true"
    if value is False:
        return "false"
    if isinstance(value, str):
        return value
    return json.dumps(value, indent=2, default=str)


def _document_with_inline_curations(document: dict, published_fields: list[dict]) -> dict:
    """Annotate a display-only copy with published values at their exact paths."""
    display_document = copy.deepcopy(document)
    for field in published_fields:
        try:
            parent, field_name = resolve_parent_and_field(
                display_document, field.get("path") or []
            )
        except (TypeError, ValueError):
            continue
        parent[field_name] = {
            _INLINE_CURATION_MARKER: True,
            "curated_display": _inline_curation_display_value(
                field.get("effective"), field.get("path")
            ),
            "original_display": _inline_curation_display_value(
                field.get("baseline"), field.get("path")
            ),
            "original_missing_display": (
                "not classified"
                if field.get("path") == ["is_generic_structure"]
                else "not present"
            ),
            "original_present": field.get("baseline_present", True),
        }
    return display_document


def _template_supports_inline_curations(template_name: str) -> bool:
    return template_name == "document.html"


def _record_curation_summary(published_fields: list[dict]) -> dict:
    """Summarize active field curations without repeating per-field metadata."""
    curator_names = []
    missing_curator_count = 0
    publication_dates = []
    missing_date_count = 0
    for field in published_fields:
        published_by = field.get("published_by") or {}
        curator_name = str(
            published_by.get("name") or published_by.get("id") or ""
        ).strip()
        if curator_name:
            if curator_name not in curator_names:
                curator_names.append(curator_name)
        else:
            missing_curator_count += 1

        published_at = str(field.get("published_at") or "").strip()
        try:
            published_date = datetime.fromisoformat(
                published_at.replace("Z", "+00:00")
            ).date()
        except ValueError:
            missing_date_count += 1
        else:
            publication_dates.append(published_date)

    visible_curators = curator_names[:2]
    additional_curator_count = max(0, len(curator_names) - len(visible_curators))
    if visible_curators:
        curator_text = "Curated by " + ", ".join(visible_curators)
        if additional_curator_count:
            curator_text += f" + {additional_curator_count} other"
        if missing_curator_count:
            curator_text += (
                f" · {missing_curator_count} curation"
                f"{'s' if missing_curator_count != 1 else ''} with no curator recorded"
            )
    else:
        curator_text = "Curator not recorded"

    if publication_dates:
        first_date = min(publication_dates)
        last_date = max(publication_dates)
        first_label = f"{first_date.strftime('%b')} {first_date.day}, {first_date.year}"
        last_label = f"{last_date.strftime('%b')} {last_date.day}, {last_date.year}"
        date_text = (
            f"Published {first_label}"
            if first_date == last_date
            else f"Published {first_label}–{last_label}"
        )
        if missing_date_count:
            date_text += (
                f" · {missing_date_count} curation"
                f"{'s' if missing_date_count != 1 else ''} with no publication date recorded"
            )
    else:
        date_text = "Publication date not recorded"

    return {
        "field_count": len(published_fields),
        "curator_text": curator_text,
        "date_text": date_text,
    }


def _record_property_cart_view(cart: dict) -> dict:
    def decision_row(decision: dict) -> dict:
        path = decision.get("path") or []
        generic_structure = path == ["is_generic_structure"]
        value = decision.get("value")
        if generic_structure:
            value_label = (
                "Generic"
                if value is True
                else "Specific / non-generic"
                if value is False
                else "Unknown / unclassified"
            )
        else:
            value_label = json.dumps(value, default=str)
        return {
            **decision,
            "path_label": (
                "Generic structure" if generic_structure
                else _display_record_property_path(path)
            ),
            "value_json": json.dumps(value, default=str),
            "value_label": value_label,
            "restore_label": (
                "use evidence graph value"
                if generic_structure else "restore loaded value"
            ),
        }

    operations = []
    for operation in cart.get("operations") or []:
        target = operation.get("target") or {}
        operations.append({
            **operation,
            "target_label": str(target.get("id") or ""),
            "decision_rows": [
                decision_row(decision)
                for decision in operation.get("decisions") or []
            ],
        })
    return {**cart, "operations": operations}


def _curator_identity(request: Request, payload: Optional[dict] = None) -> tuple[str, str]:
    payload = payload or {}
    header_identity = next((
        request.headers.get(header)
        for header in (
            "x-auth-request-email",
            "x-auth-request-user",
            "x-forwarded-email",
            "x-forwarded-user",
            "remote-user",
        )
        if request.headers.get(header)
    ), None)
    supplied_identity = str(payload.get("curator") or "").strip()
    curator_id = str(header_identity or supplied_identity).strip()
    curator_name = str(payload.get("curator_name") or curator_id).strip()
    if not curator_id:
        raise HTTPException(status_code=400, detail="Enter your curator name before using the curation cart.")
    return curator_id, curator_name


def _metabolite_edge_decision_operation(
    action: str,
    start_id: str,
    end_id: str,
    note: str = "",
) -> dict:
    if action not in {"remove_edge", "retain_edge"}:
        raise HTTPException(status_code=400, detail="Edge curation action must be remove_edge or retain_edge.")
    left = _normalize_ramp_denylist_identifier(start_id)
    right = _normalize_ramp_denylist_identifier(end_id)
    if not left or not right:
        raise HTTPException(status_code=400, detail="Both edge endpoints must be valid prefixed identifiers.")
    if left == right:
        raise HTTPException(status_code=400, detail="A curation cannot target a self-edge.")
    left, right = sorted((left, right))
    operation = {
        "action": action,
        "edge_type": "MetaboliteIdentifierMappingEdge",
        "start_id": left,
        "end_id": right,
        "symmetric": True,
    }
    clean_note = str(note or "").strip()
    if clean_note:
        operation["note"] = clean_note
    return operation


def _metabolite_edge_removal_operation(start_id: str, end_id: str, note: str = "") -> dict:
    return _metabolite_edge_decision_operation("remove_edge", start_id, end_id, note)


def _metabolite_property_operation(
    identifier: str,
    values: dict,
    remove_overrides: List[str],
    note: str = "",
) -> dict:
    normalized_id = _normalize_ramp_denylist_identifier(identifier)
    if not normalized_id:
        raise HTTPException(status_code=400, detail="A valid prefixed metabolite identifier is required.")
    try:
        db = get_db("metabolite_harmonization")
        documents = list(db.aql.execute(
            """
            FOR d IN MetaboliteIdentifier
              FILTER d.id == @target_id
              LIMIT 2
              RETURN d
            """,
            bind_vars={"target_id": normalized_id},
            max_runtime=30,
        ))
        if len(documents) != 1:
            raise ValueError(
                f"MetaboliteIdentifier:{normalized_id} matched {len(documents)} records"
            )
        schema_fields = (
            _get_collection_schema_entry(db, "MetaboliteIdentifier").get("fields")
            or {}
        )
        raw_decisions = [
            {"path": [property_name], "mode": "set", "value": value}
            for property_name, value in values.items()
        ]
        raw_decisions.extend(
            {"path": [property_name], "mode": "remove_override"}
            for property_name in remove_overrides
        )
        _curation_type, operation = _build_record_property_operation(
            documents[0],
            curation_set="metabolite_harmonization",
            model_type="MetaboliteIdentifier",
            schema_fields=schema_fields,
            raw_decisions=raw_decisions,
            note=note,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return operation


def _metabolite_generic_structure_operation(identifier: str, value: Any, note: str = "") -> dict:
    if value == "detected":
        return _metabolite_property_operation(identifier, {}, ["is_generic_structure"], note)
    if value is None or type(value) is bool:
        return _metabolite_property_operation(
            identifier, {"is_generic_structure": value}, [], note
        )
    raise HTTPException(status_code=400, detail="Generic-structure value must be true, false, null, or detected.")


def _metabolite_record_suppression_operation(
    action: str,
    identifier: str,
    note: str = "",
) -> dict:
    if action not in {"suppress_record", "restore_record"}:
        raise HTTPException(status_code=400, detail="Record action must suppress or restore harmonization.")
    normalized_id = _normalize_ramp_denylist_identifier(identifier)
    if not normalized_id:
        raise HTTPException(status_code=400, detail="A valid prefixed metabolite identifier is required.")
    operation = {
        "action": action,
        "target": {
            "kind": "node",
            "model_type": "MetaboliteIdentifier",
            "id": normalized_id,
        },
    }
    clean_note = str(note or "").strip()
    if clean_note:
        operation["note"] = clean_note
    try:
        validate_operation(METABOLITE_RECORD_SUPPRESSIONS, operation)
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    return operation


_QA_CURATION_TYPES = (
    CHEBI_RECORD_PROPERTIES,
    METABOLITE_EQUIVALENCE_EDGES,
    METABOLITE_EXPECTED_CLIQUES,
    METABOLITE_MW_ADJUDICATIONS,
    METABOLITE_RECORD_PROPERTIES,
    METABOLITE_RECORD_SUPPRESSIONS,
)


def _combined_curation_cart(storage, curator_id: str, curator_name: str) -> dict:
    carts = {
        curation_type: load_cart(storage, curation_type, curator_id, curator_name)
        for curation_type in _QA_CURATION_TYPES
    }
    operations = []
    for curation_type, cart in carts.items():
        for operation in cart.get("operations") or []:
            presented = {**operation, "curation_type": curation_type}
            if curation_type in {CHEBI_RECORD_PROPERTIES, METABOLITE_RECORD_PROPERTIES}:
                presented = _record_property_cart_view(
                    {"operations": [presented]}
                )["operations"][0]
            operations.append(presented)
    return {
        "format_version": 2,
        "curator": {"id": curator_id, "name": curator_name},
        "operations": operations,
        "operation_count": len(operations),
        "carts": carts,
    }


def _normalize_metabolite_rule_ids(rule_ids: Optional[List[str]]) -> List[str]:
    allowed = {rule["id"] for rule in _METABOLITE_HARMONIZATION_RULES}
    normalized = []
    for rule_id in rule_ids or []:
        if rule_id == "ignore_ramp_mapping_denylist":
            rule_id = "apply_curations"
        if rule_id in allowed and rule_id not in normalized:
            normalized.append(rule_id)
    return normalized


def _default_metabolite_rule_parameters(rule_id: str) -> dict:
    rule = next((item for item in _METABOLITE_HARMONIZATION_RULES if item["id"] == rule_id), None)
    if not rule:
        return {}
    return {
        parameter["id"]: parameter.get("default")
        for parameter in rule.get("parameters", [])
    }


def _normalize_metabolite_rule_parameters(rule_ids: List[str], rule_parameters: Optional[dict]) -> dict:
    normalized = {}
    submitted = rule_parameters or {}
    for rule_id in rule_ids:
        rule = next((item for item in _METABOLITE_HARMONIZATION_RULES if item["id"] == rule_id), None)
        defaults = _default_metabolite_rule_parameters(rule_id)
        if not defaults:
            continue
        parameter_types = {
            parameter["id"]: parameter.get("type")
            for parameter in (rule or {}).get("parameters", [])
        }
        normalized[rule_id] = {}
        for parameter_id, default_value in defaults.items():
            raw_value = submitted.get(rule_id, {}).get(parameter_id, default_value)
            if parameter_id == "mw_cutoff":
                parsed_value = _parse_optional_float(raw_value)
                normalized[rule_id][parameter_id] = parsed_value if parsed_value is not None else default_value
            elif parameter_types.get(parameter_id) == "textarea":
                normalized[rule_id][parameter_id] = _normalize_prefix_list(raw_value)
            elif parameter_types.get(parameter_id) == "boolean":
                if isinstance(raw_value, bool):
                    normalized[rule_id][parameter_id] = raw_value
                else:
                    normalized[rule_id][parameter_id] = str(raw_value).strip().lower() not in {
                        "0", "false", "no", "off", "disabled",
                    }
            else:
                normalized[rule_id][parameter_id] = raw_value
    return normalized


def _metabolite_rule_metadata(rule_ids: List[str], rule_parameters: Optional[dict] = None) -> List[dict]:
    rule_by_id = {rule["id"]: rule for rule in _METABOLITE_HARMONIZATION_RULES}
    parameter_by_rule_id = rule_parameters or {}
    metadata = []
    for rule_id in rule_ids:
        if rule_id not in rule_by_id:
            continue
        item = {
            "id": rule_id,
            "label": rule_by_id[rule_id]["label"],
            "description": rule_by_id[rule_id]["description"],
        }
        if parameter_by_rule_id.get(rule_id):
            item["parameters"] = parameter_by_rule_id[rule_id]
        metadata.append(item)
    return metadata


def _metabolite_harmonization_evidence_display(edge: dict) -> dict:
    """Return the graph label and class for a materialized evidence edge."""
    sources = edge.get("sources") or []
    details = edge.get("details") or []
    rule_id = edge.get("rule_id") or next(
        (detail.get("rule_id") for detail in details if detail.get("rule_id")),
        None,
    )
    is_ifx_rule = bool(rule_id) or "IFX Harmonization Rule" in sources
    if not is_ifx_rule:
        return {
            "label": "equiv" if not sources else ", ".join(sources[:3]),
            "class_name": "metabolite-equivalence-edge",
            "rule_id": None,
            "rule_label": None,
        }

    rule_metadata = _metabolite_rule_metadata([rule_id]) if rule_id else []
    rule_label = rule_metadata[0]["label"] if rule_metadata else rule_id
    return {
        "label": f"IFX: {rule_label}" if rule_label else "IFX Harmonization Rule",
        "class_name": "metabolite-equivalence-edge ifx-harmonization-rule-edge",
        "rule_id": rule_id,
        "rule_label": rule_label,
    }


def _chunked_records(records: List[dict], chunk_size: int = 10000):
    for index in range(0, len(records), chunk_size):
        yield records[index:index + chunk_size]


def _load_generic_structure_ids(db, curated_overrides: Optional[dict[str, bool]] = None) -> set:
    classifications = _load_generic_structure_classifications(db)
    classifications.update(curated_overrides or {})
    return {
        identifier
        for identifier, value in classifications.items()
        if value is True
    }


def _load_generic_structure_classifications(db) -> dict[str, Optional[bool]]:
    return {
        row["id"]: row["is_generic_structure"]
        for row in db.aql.execute(
        """
        FOR d IN MetaboliteIdentifier
          FILTER d.is_generic_structure == true OR d.is_generic_structure == false
          RETURN {id: d.id, is_generic_structure: d.is_generic_structure}
        """,
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
        )
    }


def _load_record_overlays(
    db,
    snapshot,
    *,
    curation_set: str,
    model_type: str,
) -> tuple[dict[str, dict], list[dict]]:
    """Project published field corrections without mutating the evidence graph."""
    decisions = [
        decision
        for decision in snapshot.record_property_decisions_for_set(
            curation_set
        )
        if decision.target.get("model_type") == model_type
    ]
    target_ids = sorted({decision.target.get("id") for decision in decisions})
    if not target_ids:
        return {}, []
    documents = {
        document["id"]: document
        for document in db.aql.execute(
            """
            FOR d IN @@collection
              FILTER d.id IN @target_ids
              RETURN d
            """,
            bind_vars={"@collection": model_type, "target_ids": target_ids},
            max_runtime=60,
        )
    }
    missing = sorted(set(target_ids) - set(documents))
    if missing:
        raise ValueError(
            "Record-property curation targets are missing from the evidence graph: "
            + ", ".join(missing[:10])
        )
    schema_fields = (_get_collection_schema_entry(
        db, model_type
    ).get("fields") or {})
    if not schema_fields:
        raise ValueError(
            f"{model_type} does not publish a schema for record-property curations"
        )
    overlays = dict(documents)
    reports = []
    decisions_by_target = {}
    for decision in decisions:
        decisions_by_target.setdefault(decision.target["id"], []).append(decision)
    for target_id, target_decisions in decisions_by_target.items():
        overlays[target_id], target_reports = project_record_decisions(
            overlays[target_id], target_decisions, schema_fields, model_type=model_type,
        )
        reports.extend(target_reports)
    return overlays, reports


def _load_metabolite_identifier_record_overlays(db, snapshot):
    return _load_record_overlays(
        db,
        snapshot,
        curation_set="metabolite_harmonization",
        model_type="MetaboliteIdentifier",
    )


def _effective_inchi_key_matches(
    chem_props: Iterable[dict],
    mode: str,
    key_kind: str = "effective",
) -> tuple[List[str], bool]:
    """Select InChIKey matches per chemistry record with reported-first fallback."""
    if key_kind not in {"effective", "reported", "derived"}:
        raise ValueError(f"Unsupported InChIKey kind: {key_kind}")
    matches = set()
    used_derived_fallback = False
    for prop in chem_props or []:
        reported_full = _normalize_inchi_key_duplex(prop.get("inchi_key"))
        reported_prefixes = {
            value
            for value in (
                _normalize_inchi_key_prefix(prop.get("inchi_key")),
                _normalize_inchi_key_prefix(prop.get("inchi_key_prefix")),
            )
            if value
        }
        derived_full = _normalize_inchi_key_duplex(prop.get("derived_inchi_key"))
        derived_prefixes = {
            value
            for value in (
                _normalize_inchi_key_prefix(prop.get("derived_inchi_key")),
                _normalize_inchi_key_prefix(prop.get("derived_inchi_key_prefix")),
            )
            if value
        }
        if mode == "prefix":
            reported_matches = reported_prefixes
            derived_matches = derived_prefixes
        elif mode == "duplex":
            reported_matches = {reported_full} if reported_full else set()
            derived_matches = {derived_full} if derived_full else set()
        else:
            raise ValueError(f"Unsupported InChIKey match mode: {mode}")

        if key_kind == "reported":
            selected = reported_matches
        elif key_kind == "derived":
            selected = derived_matches
        elif reported_matches:
            selected = reported_matches
        else:
            selected = derived_matches
            used_derived_fallback = used_derived_fallback or bool(selected)
        matches.update(selected)
    return sorted(matches), used_derived_fallback


def _iter_metabolite_identifier_inchi_key_matches(
    db,
    mode: str,
    mw_cutoff: Optional[float] = None,
    key_kind: str = "effective",
    record_overlays: Optional[dict[str, dict]] = None,
):
    if key_kind not in {"effective", "reported", "derived"}:
        raise ValueError(f"Unsupported InChIKey kind: {key_kind}")
    cursor = db.aql.execute(
        """
        FOR d IN MetaboliteIdentifier
          LET chem_props = (
            FOR prop IN d.chem_props || []
              FILTER prop.inchi_key != null
                OR prop.inchi_key_prefix != null
                OR prop.derived_inchi_key != null
                OR prop.derived_inchi_key_prefix != null
              RETURN KEEP(
                prop,
                "inchi_key",
                "inchi_key_prefix",
                "derived_inchi_key",
                "derived_inchi_key_prefix"
              )
          )
          LET masses = UNIQUE(
            FLATTEN(
              FOR prop IN d.chem_props || []
                RETURN [prop.mw, prop.monoisotopic_mass]
            )
          )
          FILTER LENGTH(chem_props) > 0 OR d.id IN @overlay_ids
          RETURN {id: d.id, chem_props: chem_props, masses: masses}
        """,
        bind_vars={"overlay_ids": sorted((record_overlays or {}).keys())},
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
    )
    for row in cursor:
        overlay = (record_overlays or {}).get(row["id"])
        if overlay is not None:
            overlay_chem_props = overlay.get("chem_props") or []
            key_fields = (
                "inchi_key",
                "inchi_key_prefix",
                "derived_inchi_key",
                "derived_inchi_key_prefix",
            )
            row = {
                **row,
                "chem_props": [
                    {key: prop.get(key) for key in key_fields}
                    for prop in overlay_chem_props
                    if any(prop.get(key) is not None for key in key_fields)
                ],
                "masses": [
                    value
                    for prop in overlay_chem_props
                    for value in (prop.get("mw"), prop.get("monoisotopic_mass"))
                ],
            }
        numeric_masses = [
            parsed
            for parsed in (_parse_optional_float(value) for value in row.get("masses", []))
            if parsed is not None
        ]
        mass = min(numeric_masses) if numeric_masses else None
        effective_mode = mode
        if mode == "mw_cutoff":
            effective_mode = "prefix" if mass is not None and mass >= (mw_cutoff or 500) else "duplex"
        matches, used_derived_fallback = _effective_inchi_key_matches(
            row.get("chem_props") or [],
            effective_mode,
            key_kind,
        )
        if matches:
            yield row["id"], matches, effective_mode, mass, used_derived_fallback


def _ensure_harmonization_pipeline_collections(db) -> None:
    for collection_name in (
        _HARMONIZATION_PIPELINE_COLLECTION,
        _HARMONIZATION_PIPELINE_RUN_COLLECTION,
        _HARMONIZATION_STAGE_COLLECTION,
        _HARMONIZED_METABOLITE_COLLECTION,
        _HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION,
    ):
        if not db.has_collection(collection_name):
            db.create_collection(collection_name)
    for collection_name in (
        _HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION,
        _HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION,
    ):
        if not db.has_collection(collection_name):
            db.create_collection(collection_name, edge=True)
    graph = db.graph("graph") if db.has_graph("graph") else db.create_graph("graph")
    definitions = {
        definition["edge_collection"]: definition
        for definition in graph.edge_definitions()
    }

    def upsert_edge_definition(collection_name: str, from_collections: List[str], to_collections: List[str]) -> None:
        existing = definitions.get(collection_name)
        if existing is None:
            graph.create_edge_definition(collection_name, from_collections, to_collections)
            return
        updated_from = sorted(set(existing.get("from_vertex_collections", []) + from_collections))
        updated_to = sorted(set(existing.get("to_vertex_collections", []) + to_collections))
        if (
            updated_from != sorted(existing.get("from_vertex_collections", []))
            or updated_to != sorted(existing.get("to_vertex_collections", []))
        ):
            graph.replace_edge_definition(collection_name, updated_from, updated_to)

    upsert_edge_definition(
        _HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION,
        [_HARMONIZED_METABOLITE_COLLECTION],
        ["MetaboliteIdentifier"],
    )
    upsert_edge_definition(
        _HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION,
        ["MetaboliteIdentifier"],
        ["MetaboliteIdentifier"],
    )


def _canonical_json_digest(payload: dict, length: int = 16) -> str:
    return hashlib.sha1(json.dumps(payload, sort_keys=True, separators=(",", ":")).encode("utf-8")).hexdigest()[:length]


def _harmonization_source_build_fingerprint(db) -> Optional[str]:
    """Return the authoritative ETL build identity without scanning graph data."""
    if not db.has_collection("metadata_store"):
        return None
    metadata_doc = db.collection("metadata_store").get("etl_metadata") or {}
    metadata = metadata_doc.get("value") or {}
    return source_build_fingerprint(metadata)


def _harmonization_graph_fingerprint(db) -> dict:
    row = list(db.aql.execute(
        """
        RETURN {
          metabolite_identifier_count: LENGTH(MetaboliteIdentifier),
          equivalence_edge_count: LENGTH(MetaboliteIdentifierMappingEdge)
        }
        """,
        max_runtime=120,
    ))[0]
    source_collection_revisions = {
        collection_name: db.collection(collection_name).revision()
        for collection_name in _HARMONIZATION_SOURCE_COLLECTIONS
        if db.has_collection(collection_name)
    }
    source_graph_fingerprint = _harmonization_source_build_fingerprint(db)
    source_graph_content_fingerprint = payload_sha256({
        "source_graph_fingerprint": source_graph_fingerprint,
        "source_collection_revisions": source_collection_revisions,
    })
    return {
        "database": "metabolite_harmonization",
        "engine_version": _HARMONIZATION_ENGINE_VERSION,
        "source_graph_fingerprint": source_graph_fingerprint,
        "source_graph_content_fingerprint": source_graph_content_fingerprint,
        "source_collection_revisions": source_collection_revisions,
        # Counts are useful diagnostics but are not graph identity: two rebuilds
        # can change records without changing collection cardinality. Arango's
        # collection revision tokens make deterministic stage reuse sensitive
        # to inserts, deletes, and updates without scanning every document.
        **row,
    }


def _harmonization_stage_key(rule_ids: List[str], rule_parameters: dict, graph_fingerprint: dict) -> str:
    digest = _canonical_json_digest({
        "kind": "harmonization_stage",
        "graph_fingerprint": graph_fingerprint,
        "rule_ids": rule_ids,
        "rule_parameters": rule_parameters,
    })
    if not rule_ids:
        return f"baseline-{digest}"
    return f"stage-{len(rule_ids):02d}-{digest}"


def _cumulative_rule_parameters(rule_ids: List[str], rule_parameters: dict) -> dict:
    return {
        rule_id: rule_parameters.get(rule_id, {})
        for rule_id in rule_ids
        if rule_parameters.get(rule_id)
    }


def _delete_harmonization_stage_artifacts(db, stage_key: str) -> None:
    delete_stage_artifacts(db, stage_key)


def _stage_is_referenced_by_any_run(db, stage_key: str) -> bool:
    if not db.has_collection(_HARMONIZATION_PIPELINE_RUN_COLLECTION):
        return False
    rows = list(db.aql.execute(
        f"""
        FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
          FILTER @stage_key IN (r.stage_keys || [])
          LIMIT 1
          RETURN 1
        """,
        bind_vars={"stage_key": stage_key},
        max_runtime=120,
    ))
    return bool(rows)


def _delete_unreferenced_harmonization_stages(db, orphaned_before: Optional[str] = None) -> List[str]:
    if (
        not db.has_collection(_HARMONIZATION_STAGE_COLLECTION)
        or not db.has_collection(_HARMONIZATION_PIPELINE_RUN_COLLECTION)
    ):
        return []
    cutoff = orphaned_before or (
        datetime.now(timezone.utc) - _HARMONIZATION_ORPHAN_GRACE_PERIOD
    ).isoformat()
    candidate_stage_keys = list(db.aql.execute(
        f"""
        FOR s IN {_HARMONIZATION_STAGE_COLLECTION}
          FILTER s.status == "complete"
          LET stage_time = s.updated_at != null ? s.updated_at : s.created_at
          FILTER stage_time != null AND stage_time < @orphaned_before
          LET referenced = FIRST(
            FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
              FILTER s._key IN (r.stage_keys || [])
              LIMIT 1
              RETURN true
          )
          FILTER referenced == null
          SORT s._key
          RETURN s._key
        """,
        bind_vars={"orphaned_before": cutoff},
        max_runtime=120,
    ))

    deleted_stage_keys = []
    stage_collection = db.collection(_HARMONIZATION_STAGE_COLLECTION)
    for stage_key in candidate_stage_keys:
        # Recheck immediately before deletion in case another pipeline run reused
        # this deterministic stage while the orphan scan was in progress.
        if _stage_is_referenced_by_any_run(db, stage_key):
            continue
        _delete_harmonization_stage_artifacts(db, stage_key)
        if stage_collection.has(stage_key):
            stage_collection.delete(stage_key)
        deleted_stage_keys.append(stage_key)
    return deleted_stage_keys


def _delete_previous_harmonization_pipeline_runs(db, pipeline_key: str, current_run_key: str) -> dict:
    previous_runs = list(db.aql.execute(
        f"""
        FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
          FILTER r.pipeline_key == @pipeline_key
          FILTER r._key != @current_run_key
          RETURN KEEP(r, "_key", "stage_keys")
        """,
        bind_vars={"pipeline_key": pipeline_key, "current_run_key": current_run_key},
        max_runtime=120,
    ))
    previous_run_keys = [run["_key"] for run in previous_runs]
    if previous_run_keys:
        db.aql.execute(
            f"""
            FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
              FILTER r._key IN @run_keys
              REMOVE r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
            """,
            bind_vars={"run_keys": previous_run_keys},
            max_runtime=120,
        )

    # Stages are deterministic, shared, immutable cache entries. Deleting their
    # large artifact sets inline creates a check-then-delete race with another
    # pipeline that is concurrently reusing the same stage. Reclaim them only
    # through a future exclusive maintenance operation.
    return {
        "deleted_run_keys": previous_run_keys,
        "deleted_stage_keys": [],
        "deleted_superseded_stage_keys": [],
        "deleted_orphan_stage_keys": [],
    }


def _delete_harmonization_pipeline(db, pipeline_key: str) -> dict:
    _ensure_harmonization_pipeline_collections(db)
    pipeline = _get_harmonization_pipeline(db, pipeline_key)
    db.aql.execute(
        f"""
        FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
          FILTER r.pipeline_key == @pipeline_key
          REMOVE r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
        """,
        bind_vars={"pipeline_key": pipeline_key},
        max_runtime=120,
    )
    pipeline_collection = db.collection(_HARMONIZATION_PIPELINE_COLLECTION)
    if pipeline_collection.has(pipeline_key):
        pipeline_collection.delete(pipeline_key)
    return {
        "pipeline": pipeline,
        "deleted_stage_keys": [],
    }


def _load_metabolite_identifier_source_support(db) -> Dict[str, set]:
    support_by_id: Dict[str, set] = {}
    cursor = db.aql.execute(
        """
        FOR d IN MetaboliteIdentifier
          LET support = UNIQUE(FLATTEN([
            (FOR item IN d.sources || [] RETURN IS_OBJECT(item) ? (item.name != null ? item.name : (item.source != null ? item.source : item.id)) : item),
            (FOR item IN d.names || [] RETURN item.source),
            (FOR item IN d.synonyms || [] RETURN item.source),
            (FOR item IN d.chem_props || [] RETURN item.source)
          ]))
          RETURN {id: d.id, support: support}
        """,
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
    )
    for row in cursor:
        sources = {
            _metabolite_support_source_name(source)
            for source in row.get("support") or []
            if source is not None and str(source).strip()
        }
        support_by_id[row["id"]] = sources or {"__unsourced_node__"}
    return support_by_id


def _metabolite_support_source_name(source) -> str:
    return str(source).strip().split("\t", 1)[0]


def _wikipathways_xref_only_ids_from_rows(
    mapping_rows: Iterable[dict],
    pathway_primary_ids: Iterable[str],
    ignored_source_fields: Iterable[str],
) -> set:
    ignored_fields = {
        str(source_field).strip().lower()
        for source_field in ignored_source_fields
        if source_field is not None and str(source_field).strip()
    }
    if not ignored_fields:
        return set()

    primary_ids = {
        identifier
        for identifier in pathway_primary_ids
        if identifier
    }
    source_fields_by_xref_id: Dict[str, set] = {}
    for edge in mapping_rows:
        start_id = edge.get("start_id")
        end_id = edge.get("end_id")
        for detail in edge.get("details") or []:
            if str(detail.get("source") or "").strip().lower() != "wikipathways":
                continue
            source_id = detail.get("source_id")
            source_field = str(detail.get("source_field") or "").strip().lower()
            if not source_id or not source_field:
                continue
            primary_ids.add(source_id)
            if source_id == start_id:
                xref_id = end_id
            elif source_id == end_id:
                xref_id = start_id
            else:
                continue
            if xref_id and xref_id != source_id:
                source_fields_by_xref_id.setdefault(xref_id, set()).add(source_field)

    return {
        identifier
        for identifier, source_fields in source_fields_by_xref_id.items()
        if identifier not in primary_ids
        and source_fields
        and source_fields.issubset(ignored_fields)
    }


def _load_wikipathways_xref_only_ids_for_ignored_source_fields(
    db,
    ignored_source_fields: Iterable[str],
) -> set:
    mapping_rows = db.aql.execute(
        """
        FOR e IN MetaboliteIdentifierMappingEdge
          FILTER LENGTH(
            FOR detail IN e.details || []
              FILTER LOWER(detail.source || "") == "wikipathways"
              RETURN 1
          ) > 0
          RETURN {
            start_id: e.start_id,
            end_id: e.end_id,
            details: e.details || []
          }
        """,
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
    )
    pathway_primary_ids = db.aql.execute(
        """
        FOR e IN MetabolitePathwayEdge
          FILTER LENGTH(
            FOR detail IN e.details || []
              FILTER LOWER(detail.source || "") == "wikipathways"
              RETURN 1
          ) > 0
          RETURN DISTINCT e.start_id
        """,
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
    )
    return _wikipathways_xref_only_ids_from_rows(
        mapping_rows,
        pathway_primary_ids,
        ignored_source_fields,
    )


def _filter_identifier_support_for_rules(
    support_by_id: Dict[str, set],
    rule_ids: List[str],
    rule_parameters: dict,
    wikipathways_ignored_source_field_ids: Optional[set] = None,
) -> Dict[str, set]:
    filtered = {identifier: set(sources) for identifier, sources in support_by_id.items()}
    source_prefix_rules = [
        ("ignore_hmdb_prefixes", "hmdb"),
        ("ignore_wikipathways_prefixes", "wikipathways"),
        ("ignore_lipidmaps_prefixes", "lipidmaps"),
    ]
    for rule_id, source_name in source_prefix_rules:
        if rule_id not in rule_ids:
            continue
        prefixes = {
            prefix.lower()
            for prefix in rule_parameters.get(rule_id, {}).get("prefixes", [])
        }
        for identifier, sources in filtered.items():
            if (_identifier_prefix(identifier) or "").lower() in prefixes:
                filtered[identifier] = {
                    source for source in sources
                    if str(source).lower() != source_name
                }
    if "ignore_wikipathways_prefixes" in rule_ids:
        for identifier in wikipathways_ignored_source_field_ids or set():
            if identifier not in filtered:
                continue
            filtered[identifier] = {
                source for source in filtered[identifier]
                if str(source).lower() != "wikipathways"
            }
    return filtered


def _active_metabolite_identifier_mapping_edges_for_rules(
    db,
    active_ids: set,
    rule_ids: List[str],
    rule_parameters: dict,
    ramp_mapping_denylist_pairs: Optional[set[tuple[str, str]]] = None,
    generic_structure_classifications: Optional[dict[str, Optional[bool]]] = None,
) -> tuple[List[dict], dict]:
    generic_structure_ids = (
        {
            identifier
            for identifier, value in (generic_structure_classifications or {}).items()
            if value is True
        }
        if "ignore_generic_structure_mismatch" in rule_ids
        else set()
    )
    ramp_mapping_denylist_pairs = ramp_mapping_denylist_pairs or set()
    hmdb_ignored_prefixes = set()
    if "ignore_hmdb_prefixes" in rule_ids:
        hmdb_ignored_prefixes = {
            prefix.lower()
            for prefix in rule_parameters.get("ignore_hmdb_prefixes", {}).get("prefixes", [])
        }
    wikipathways_ignored_prefixes = set()
    wikipathways_ignored_source_fields = set()
    if "ignore_wikipathways_prefixes" in rule_ids:
        wikipathways_ignored_prefixes = {
            prefix.lower()
            for prefix in rule_parameters.get("ignore_wikipathways_prefixes", {}).get("prefixes", [])
        }
        wikipathways_ignored_source_fields = {
            field
            for field in rule_parameters.get("ignore_wikipathways_prefixes", {}).get("source_fields", [])
        }
    lipidmaps_ignored_prefixes = set()
    if "ignore_lipidmaps_prefixes" in rule_ids:
        lipidmaps_ignored_prefixes = {
            prefix.lower()
            for prefix in rule_parameters.get("ignore_lipidmaps_prefixes", {}).get("prefixes", [])
        }

    active_edges = []
    summary = {
        "generic_structure_identifier_count": len(generic_structure_ids),
        "ramp_denylist_pair_count": len(ramp_mapping_denylist_pairs),
        "hmdb_ignored_prefixes": sorted(hmdb_ignored_prefixes),
        "hmdb_ignored_prefix_count": len(hmdb_ignored_prefixes),
        "wikipathways_ignored_prefixes": sorted(wikipathways_ignored_prefixes),
        "wikipathways_ignored_prefix_count": len(wikipathways_ignored_prefixes),
        "wikipathways_ignored_source_fields": sorted(wikipathways_ignored_source_fields),
        "wikipathways_ignored_source_field_count": len(wikipathways_ignored_source_fields),
        "lipidmaps_ignored_prefixes": sorted(lipidmaps_ignored_prefixes),
        "lipidmaps_ignored_prefix_count": len(lipidmaps_ignored_prefixes),
        "ignored_edge_count": 0,
        "generic_structure_ignored_edge_count": 0,
        "ramp_denylist_ignored_edge_count": 0,
        "hmdb_prefix_ignored_detail_count": 0,
        "hmdb_prefix_ignored_edge_count": 0,
        "wikipathways_prefix_ignored_detail_count": 0,
        "wikipathways_prefix_ignored_edge_count": 0,
        "lipidmaps_prefix_ignored_detail_count": 0,
        "lipidmaps_prefix_ignored_edge_count": 0,
    }
    cursor = db.aql.execute(
        """
        FOR e IN MetaboliteIdentifierMappingEdge
          RETURN {
            key: e._key,
            id: e.id || e._key,
            start_id: e.start_id,
            end_id: e.end_id,
            details: e.details || []
          }
        """,
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
    )
    for edge in cursor:
        start_id = edge["start_id"]
        end_id = edge["end_id"]
        if start_id not in active_ids or end_id not in active_ids:
            summary["ignored_edge_count"] += 1
            continue
        start_generic_classification = (generic_structure_classifications or {}).get(start_id)
        end_generic_classification = (generic_structure_classifications or {}).get(end_id)
        if (
            "ignore_generic_structure_mismatch" in rule_ids
            and {start_generic_classification, end_generic_classification} == {True, False}
        ):
            summary["ignored_edge_count"] += 1
            summary["generic_structure_ignored_edge_count"] += 1
            continue
        if "apply_curations" in rule_ids and tuple(sorted((start_id, end_id))) in ramp_mapping_denylist_pairs:
            summary["ignored_edge_count"] += 1
            summary["ramp_denylist_ignored_edge_count"] += 1
            continue
        details = edge.get("details") or []
        active_details = details
        source_prefix_rules = [
            ("hmdb", hmdb_ignored_prefixes, "hmdb_prefix"),
            ("wikipathways", wikipathways_ignored_prefixes, "wikipathways_prefix"),
            ("lipidmaps", lipidmaps_ignored_prefixes, "lipidmaps_prefix"),
        ]
        for source_name, ignored_prefixes, summary_prefix in source_prefix_rules:
            if not ignored_prefixes:
                continue
            endpoint_has_ignored_prefix = (
                (_identifier_prefix(start_id) or "").lower() in ignored_prefixes
                or (_identifier_prefix(end_id) or "").lower() in ignored_prefixes
            )
            if endpoint_has_ignored_prefix:
                before_count = len(active_details)
                active_details = [
                    detail for detail in active_details
                    if str(detail.get("source") or "").lower() != source_name
                ]
                removed_count = before_count - len(active_details)
                summary[f"{summary_prefix}_ignored_detail_count"] += removed_count
                if removed_count and not active_details:
                    summary[f"{summary_prefix}_ignored_edge_count"] += 1
        if wikipathways_ignored_source_fields:
            before_count = len(active_details)
            active_details = [
                detail for detail in active_details
                if not (
                    str(detail.get("source") or "").lower() == "wikipathways"
                    and detail.get("source_field") in wikipathways_ignored_source_fields
                )
            ]
            removed_count = before_count - len(active_details)
            if removed_count:
                summary["wikipathways_prefix_ignored_detail_count"] += removed_count
                if not active_details:
                    summary["wikipathways_prefix_ignored_edge_count"] += 1
        if details and not active_details:
            summary["ignored_edge_count"] += 1
            continue
        sources = sorted({
            str(detail.get("source")).strip()
            for detail in active_details
            if detail.get("source") is not None and str(detail.get("source")).strip()
        })
        active_edges.append({
            **edge,
            "details": active_details,
            "sources": sources,
        })
    summary["active_edge_count"] = len(active_edges)
    return active_edges, summary


def _expected_clique_assertion_edges(
    active_ids: set,
    assertions: List[dict],
) -> tuple[List[dict], dict]:
    synthetic_edges = []
    skipped_assertions = []
    applied_assertion_count = 0
    for assertion in assertions:
        member_ids = assertion.get("member_ids") or []
        missing_member_ids = sorted(set(member_ids) - active_ids)
        if missing_member_ids:
            skipped_assertions.append({
                "assertion_id": assertion.get("assertion_id"),
                "name": assertion.get("name"),
                "missing_member_ids": missing_member_ids,
            })
            continue
        if len(member_ids) < 2:
            continue
        applied_assertion_count += 1
        anchor_id = member_ids[0]
        for end_id in member_ids[1:]:
            assertion_id = assertion["assertion_id"]
            edge_key = f"assertion-{assertion_id}-{_canonical_json_digest({'start_id': anchor_id, 'end_id': end_id})}"
            detail = {
                "source": "Expected Clique Assertion",
                "synthetic": True,
                "fallback": True,
                "assertion_id": assertion_id,
                "assertion_name": assertion.get("name"),
                "rationale": assertion.get("rationale"),
                "curation_batch_id": assertion.get("curation_batch_id"),
                "published_at": assertion.get("published_at"),
                "published_by": assertion.get("published_by"),
            }
            synthetic_edges.append({
                "key": edge_key,
                "id": f"SyntheticExpectedCliqueEdge:{edge_key}",
                "start_id": anchor_id,
                "end_id": end_id,
                "sources": ["Expected Clique Assertion"],
                "details": [detail],
                "synthetic": True,
                "fallback": True,
                "assertion_id": assertion_id,
            })
    return synthetic_edges, {
        "expected_clique_assertion_count": len(assertions),
        "expected_clique_assertion_applied_count": applied_assertion_count,
        "expected_clique_assertion_skipped_count": len(skipped_assertions),
        "expected_clique_synthetic_edge_count": len(synthetic_edges),
        "expected_clique_assertion_skipped_samples": skipped_assertions[:20],
    }


def _load_chebi_descendant_ids(db, root_ids: Iterable[str]) -> set:
    descendants = set()
    for root_id in root_ids:
        rows = db.aql.execute(
            """
            LET root = DOCUMENT("ChemicalEntity", @root_id)
            FOR vertex IN 0..100 INBOUND root IsAEdge
              OPTIONS {order: "bfs", uniqueVertices: "global"}
              RETURN DISTINCT vertex.id
            """,
            bind_vars={"root_id": root_id},
            batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
            max_runtime=600,
        )
        descendants.update(identifier for identifier in rows if identifier)
    return descendants


def _free_anomer_rule_edges(
    db,
    active_ids: set,
) -> tuple[List[dict], dict]:
    """Build auditable synthetic edges for qualifying ChEBI free-anomer pairs."""
    polymer_ids = _load_chebi_descendant_ids(db, CHEBI_POLYMER_ROOT_IDS)
    candidate_rows = list(db.aql.execute(
        """
        FOR edge IN IsAEdge
          FILTER edge.start_id NOT IN @polymer_ids
          FILTER edge.end_id NOT IN @polymer_ids
          LET child = DOCUMENT("ChemicalEntity", edge.start_id)
          LET parent = DOCUMENT("ChemicalEntity", edge.end_id)
          LET child_identifier = DOCUMENT("MetaboliteIdentifier", edge.start_id)
          LET parent_identifier = DOCUMENT("MetaboliteIdentifier", edge.end_id)
          FILTER child != null AND parent != null
          FILTER child_identifier != null AND parent_identifier != null
          RETURN {
            relationship_key: edge._key,
            relationship_id: edge.id,
            relationship_predicate: edge.source_predicate,
            child: KEEP(child, "id", "name", "formula", "charge", "smiles"),
            parent: KEEP(parent, "id", "name", "formula", "charge", "smiles")
          }
        """,
        bind_vars={
            "polymer_ids": sorted(polymer_ids),
        },
        batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        max_runtime=600,
    ))
    candidate_rows = [
        row
        for row in candidate_rows
        if row.get("child", {}).get("id") in active_ids
        and row.get("parent", {}).get("id") in active_ids
    ]
    decisions = [evaluate_free_anomer_candidate(row) for row in candidate_rows]
    edges = []
    for decision in decisions:
        if not decision.get("accepted"):
            continue
        child_id = decision["child_id"]
        parent_id = decision["parent_id"]
        edge_key = f"free-anomer-{_canonical_json_digest({'child_id': child_id, 'parent_id': parent_id, 'algorithm_version': FREE_ANOMER_ALGORITHM_VERSION})}"
        detail = {
            "source": "IFX Harmonization Rule",
            "source_field": "free_anomeric_center",
            "source_id": decision.get("relationship_id") or decision.get("relationship_key") or child_id,
            "rule_id": FREE_ANOMER_RULE_ID,
            "algorithm_version": FREE_ANOMER_ALGORITHM_VERSION,
            "decision": decision["reason"],
            "chebi_relationship": decision["relationship_predicate"],
            "child_name": decision.get("child_name"),
            "parent_name": decision.get("parent_name"),
            "normalized_inchi_key": decision["normalized_inchi_key"],
            "child_source_smiles": decision["child_structure"]["source_smiles"],
            "parent_source_smiles": decision["parent_structure"]["source_smiles"],
            "child_source_inchi_key": decision["child_structure"]["source_inchi_key"],
            "parent_source_inchi_key": decision["parent_structure"]["source_inchi_key"],
            "child_anomeric_atom_index": decision["child_structure"]["anomeric_atom_index"],
            "parent_anomeric_atom_index": decision["parent_structure"]["anomeric_atom_index"],
        }
        edges.append({
            "key": edge_key,
            "id": f"SyntheticFreeAnomerEdge:{edge_key}",
            "start_id": child_id,
            "end_id": parent_id,
            "sources": ["IFX Harmonization Rule"],
            "details": [detail],
            "synthetic": True,
            "fallback": False,
            "rule_id": FREE_ANOMER_RULE_ID,
            "algorithm_version": FREE_ANOMER_ALGORITHM_VERSION,
        })
    summary = {
        **summarize_free_anomer_decisions(decisions),
        "free_anomer_algorithm_version": FREE_ANOMER_ALGORITHM_VERSION,
        "free_anomer_polymer_excluded_identifier_count": len(polymer_ids),
        "free_anomer_candidate_identifier_count": len({
            identifier
            for decision in decisions
            for identifier in (decision.get("child_id"), decision.get("parent_id"))
            if identifier
        }),
        "free_anomer_synthetic_edge_count": len(edges),
    }
    return edges, summary


def _build_harmonized_groups(
    db,
    active_ids: set,
    active_edges: List[dict],
    rule_ids: List[str],
    rule_parameters: dict,
    record_overlays: Optional[dict[str, dict]] = None,
) -> tuple[List[List[str]], dict]:
    parent: Dict[str, str] = {}
    size: Dict[str, int] = {}

    def find(identifier: str) -> str:
        while parent[identifier] != identifier:
            parent[identifier] = parent[parent[identifier]]
            identifier = parent[identifier]
        return identifier

    def add(identifier: str) -> None:
        if identifier not in parent:
            parent[identifier] = identifier
            size[identifier] = 1

    def union(left: str, right: str) -> bool:
        add(left)
        add(right)
        left_root = find(left)
        right_root = find(right)
        if left_root == right_root:
            return False
        if size[left_root] < size[right_root]:
            left_root, right_root = right_root, left_root
        parent[right_root] = left_root
        size[left_root] += size[right_root]
        return True

    for identifier in active_ids:
        add(identifier)
    for edge in active_edges:
        union(edge["start_id"], edge["end_id"])

    merge_summary = {
        "inchi_key_prefix_count": 0,
        "inchi_key_prefix_identifier_count": 0,
        "inchi_key_prefix_merge_count": 0,
        "inchi_key_duplex_count": 0,
        "inchi_key_duplex_identifier_count": 0,
        "inchi_key_duplex_merge_count": 0,
        "inchi_key_mw_cutoff_prefix_identifier_count": 0,
        "inchi_key_mw_cutoff_duplex_identifier_count": 0,
        "inchi_key_mw_cutoff_identifier_count": 0,
        "inchi_key_mw_cutoff_merge_count": 0,
        "inchi_key_prefix_derived_fallback_identifier_count": 0,
        "inchi_key_duplex_derived_fallback_identifier_count": 0,
        "inchi_key_mw_cutoff_derived_fallback_identifier_count": 0,
        "derived_inchi_key_mw_cutoff_prefix_identifier_count": 0,
        "derived_inchi_key_mw_cutoff_duplex_identifier_count": 0,
        "derived_inchi_key_mw_cutoff_identifier_count": 0,
        "derived_inchi_key_mw_cutoff_merge_count": 0,
    }

    def merge_by_inchi_key_matches(
        mode: str,
        mw_cutoff: Optional[float] = None,
        key_kind: str = "effective",
        use_record_overlays: bool = False,
    ) -> dict:
        first_identifier_by_match: Dict[str, str] = {}
        identifier_count = 0
        merge_count = 0
        prefix_identifier_count = 0
        duplex_identifier_count = 0
        derived_fallback_identifier_count = 0
        match_rows = (
            _iter_metabolite_identifier_inchi_key_matches(
                db, mode, mw_cutoff, key_kind, record_overlays
            )
            if use_record_overlays
            else _iter_metabolite_identifier_inchi_key_matches(
                db, mode, mw_cutoff, key_kind
            )
        )
        for identifier, matches, effective_mode, _mass, used_derived_fallback in match_rows:
            if identifier not in active_ids:
                continue
            identifier_count += 1
            if effective_mode == "prefix":
                prefix_identifier_count += 1
            if effective_mode == "duplex":
                duplex_identifier_count += 1
            if used_derived_fallback:
                derived_fallback_identifier_count += 1
            for match in matches:
                match_key = f"{effective_mode}:{match}"
                if match_key not in first_identifier_by_match:
                    first_identifier_by_match[match_key] = identifier
                    continue
                if union(first_identifier_by_match[match_key], identifier):
                    merge_count += 1
        return {
            "match_count": len(first_identifier_by_match),
            "identifier_count": identifier_count,
            "merge_count": merge_count,
            "prefix_identifier_count": prefix_identifier_count,
            "duplex_identifier_count": duplex_identifier_count,
            "derived_fallback_identifier_count": derived_fallback_identifier_count,
        }

    def curations_precede(rule_id: str) -> bool:
        return (
            "apply_curations" in rule_ids
            and rule_id in rule_ids
            and rule_ids.index("apply_curations") < rule_ids.index(rule_id)
        )

    if "merge_shared_inchikey_prefix" in rule_ids:
        stats = merge_by_inchi_key_matches(
            "prefix",
            use_record_overlays=curations_precede("merge_shared_inchikey_prefix"),
        )
        merge_summary["inchi_key_prefix_count"] = stats["match_count"]
        merge_summary["inchi_key_prefix_identifier_count"] = stats["identifier_count"]
        merge_summary["inchi_key_prefix_merge_count"] = stats["merge_count"]
        merge_summary["inchi_key_prefix_derived_fallback_identifier_count"] = stats[
            "derived_fallback_identifier_count"
        ]
    if "merge_shared_inchikey_duplex" in rule_ids:
        stats = merge_by_inchi_key_matches(
            "duplex",
            use_record_overlays=curations_precede("merge_shared_inchikey_duplex"),
        )
        merge_summary["inchi_key_duplex_count"] = stats["match_count"]
        merge_summary["inchi_key_duplex_identifier_count"] = stats["identifier_count"]
        merge_summary["inchi_key_duplex_merge_count"] = stats["merge_count"]
        merge_summary["inchi_key_duplex_derived_fallback_identifier_count"] = stats[
            "derived_fallback_identifier_count"
        ]
    if "merge_inchikey_by_mw_cutoff" in rule_ids:
        mw_cutoff = rule_parameters.get("merge_inchikey_by_mw_cutoff", {}).get("mw_cutoff", 500)
        stats = merge_by_inchi_key_matches(
            "mw_cutoff",
            mw_cutoff,
            use_record_overlays=curations_precede("merge_inchikey_by_mw_cutoff"),
        )
        merge_summary["inchi_key_mw_cutoff_identifier_count"] = stats["identifier_count"]
        merge_summary["inchi_key_mw_cutoff_merge_count"] = stats["merge_count"]
        merge_summary["inchi_key_mw_cutoff_prefix_identifier_count"] = stats["prefix_identifier_count"]
        merge_summary["inchi_key_mw_cutoff_duplex_identifier_count"] = stats["duplex_identifier_count"]
        merge_summary["inchi_key_mw_cutoff_derived_fallback_identifier_count"] = stats[
            "derived_fallback_identifier_count"
        ]
    if "merge_derived_inchikey_by_mw_cutoff" in rule_ids:
        mw_cutoff = rule_parameters.get("merge_derived_inchikey_by_mw_cutoff", {}).get("mw_cutoff", 500)
        stats = merge_by_inchi_key_matches(
            "mw_cutoff",
            mw_cutoff,
            "derived",
            use_record_overlays=curations_precede(
                "merge_derived_inchikey_by_mw_cutoff"
            ),
        )
        merge_summary["derived_inchi_key_mw_cutoff_identifier_count"] = stats["identifier_count"]
        merge_summary["derived_inchi_key_mw_cutoff_merge_count"] = stats["merge_count"]
        merge_summary["derived_inchi_key_mw_cutoff_prefix_identifier_count"] = stats["prefix_identifier_count"]
        merge_summary["derived_inchi_key_mw_cutoff_duplex_identifier_count"] = stats["duplex_identifier_count"]

    members_by_root: Dict[str, List[str]] = {}
    for identifier in active_ids:
        members_by_root.setdefault(find(identifier), []).append(identifier)
    non_singleton_groups = [
        sorted(members)
        for members in members_by_root.values()
        if len(members) > 1
    ]
    non_singleton_groups.sort(key=lambda members: (-len(members), members[0]))
    return non_singleton_groups, merge_summary


_SECONDARY_SOURCES = {"refmet", "pubchem", "chebi", "lipidmaps"}

_BIO_CONTEXT_EDGE_COLLECTIONS = [
    "MetabolitePathwayEdge",
    "HmdbMetaboliteProteinAssociationEdge",
    "RheaMetaboliteReactionEdge",
    "MetaboliteClassificationEdge",
    "HmdbMetaboliteOntologyEdge",
]


def _load_bio_connected_ids(db) -> set:
    connected = set()
    for coll in _BIO_CONTEXT_EDGE_COLLECTIONS:
        for doc_id in db.aql.execute(
            f"FOR e IN {coll} RETURN DISTINCT e._from",
            batch_size=_HARMONIZATION_AQL_BATCH_SIZE,
        ):
            connected.add(doc_id.split("/", 1)[-1])
    return connected


def _remove_enrichment_only_metabolites(
    active_ids: set,
    active_edges: List[dict],
    groups: List[List[str]],
    support_by_id: Dict[str, set],
    bio_connected_ids: set,
) -> tuple[set, List[dict], List[List[str]], dict]:
    removed_ids = set()

    def is_secondary_only(identifiers: Iterable[str]) -> bool:
        sources = {
            _metabolite_support_source_name(source).lower()
            for identifier in identifiers
            for source in support_by_id.get(identifier, set())
            if source is not None and str(source).strip()
        }
        if not sources or not (sources <= _SECONDARY_SOURCES):
            return False
        return not any(identifier in bio_connected_ids for identifier in identifiers)

    kept_groups = []
    removed_group_count = 0
    for members in groups:
        if is_secondary_only(members):
            removed_group_count += 1
            removed_ids.update(members)
        else:
            kept_groups.append(members)

    grouped_ids = {member_id for members in groups for member_id in members}
    enrichment_singleton_ids = {
        identifier
        for identifier in active_ids
        if identifier not in grouped_ids and is_secondary_only([identifier])
    }
    removed_ids.update(enrichment_singleton_ids)

    filtered_active_ids = set(active_ids) - removed_ids
    filtered_active_edges = [
        edge for edge in active_edges
        if edge["start_id"] in filtered_active_ids and edge["end_id"] in filtered_active_ids
    ]
    return filtered_active_ids, filtered_active_edges, kept_groups, {
        "enrichment_only_removed_metabolite_count": removed_group_count + len(enrichment_singleton_ids),
        "enrichment_only_removed_group_count": removed_group_count,
        "enrichment_only_removed_singleton_count": len(enrichment_singleton_ids),
        "enrichment_only_removed_identifier_count": len(removed_ids),
        "enrichment_only_removed_edge_count": len(active_edges) - len(filtered_active_edges),
    }


def _load_metabolite_identifier_handles(db, identifier_ids: set) -> Dict[str, str]:
    handles_by_id: Dict[str, str] = {}
    for identifiers in _chunked_records(sorted(identifier_ids), 5000):
        rows = db.aql.execute(
            """
            FOR d IN MetaboliteIdentifier
              FILTER d.id IN @ids
              RETURN {id: d.id, handle: d._id}
            """,
            bind_vars={"ids": identifiers},
            max_runtime=120,
        )
        for row in rows:
            handles_by_id[row["id"]] = row["handle"]
    missing = sorted(set(identifier_ids) - set(handles_by_id))
    if missing:
        sample = ", ".join(missing[:10])
        if len(missing) > 10:
            sample += f", ... +{len(missing) - 10} more"
        raise ValueError(f"Cannot materialize harmonization stage; missing MetaboliteIdentifier handles for {sample}")
    return handles_by_id


def _materialize_harmonization_stage(
    db,
    stage_key: str,
    stage_doc: dict,
    active_ids: set,
    active_edges: List[dict],
    groups: List[List[str]],
    generic_structure_classifications: Optional[dict[str, Optional[bool]]] = None,
) -> None:
    _delete_harmonization_stage_artifacts(db, stage_key)
    stage_collection = db.collection(_HARMONIZATION_STAGE_COLLECTION)
    materialization_started_at = datetime.now(timezone.utc).isoformat()
    stage_doc.update({
        "status": "materializing",
        "materialization_started_at": materialization_started_at,
        "materialization_completed_at": None,
        "materialization_failed_at": None,
        "materialization_error": None,
    })
    stage_collection.insert(stage_doc, overwrite=True)
    try:
        identifier_handles = _load_metabolite_identifier_handles(db, active_ids)

        active_chunk_collection = db.collection(_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION)
        active_ids_sorted = sorted(active_ids)
        for chunk_index, identifiers in enumerate(_chunked_records(active_ids_sorted, 20000)):
            active_chunk_collection.insert({
                "_key": f"{stage_key}-{chunk_index:05d}",
                "id": f"HarmonizationStageActiveIdentifierChunk:{stage_key}-{chunk_index:05d}",
                "stage_key": stage_key,
                "chunk_index": chunk_index,
                "identifier_ids": identifiers,
            }, overwrite=True)

        evidence_collection = db.collection(_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION)
        evidence_batch = []
        for edge in active_edges:
            edge_key = f"{stage_key}-{_canonical_json_digest({'edge_key': edge.get('key'), 'start_id': edge.get('start_id'), 'end_id': edge.get('end_id')})}"
            evidence_batch.append({
                "_key": edge_key,
                "_from": identifier_handles[edge["start_id"]],
                "_to": identifier_handles[edge["end_id"]],
                "id": f"HarmonizationStageEvidenceEdge:{edge_key}",
                "stage_key": stage_key,
                "raw_edge_key": edge.get("key"),
                "raw_edge_id": edge.get("id"),
                "start_id": edge["start_id"],
                "end_id": edge["end_id"],
                "sources": edge.get("sources") or [],
                "source_count": len(edge.get("sources") or []),
                "detail_count": len(edge.get("details") or []),
                "details": edge.get("details") or [],
                "synthetic": bool(edge.get("synthetic")),
                "fallback": bool(edge.get("fallback")),
                "assertion_id": edge.get("assertion_id"),
                "rule_id": edge.get("rule_id"),
                "algorithm_version": edge.get("algorithm_version"),
            })
            if len(evidence_batch) >= 10000:
                evidence_collection.insert_many(evidence_batch, overwrite=True, silent=True)
                evidence_batch = []
        if evidence_batch:
            evidence_collection.insert_many(evidence_batch, overwrite=True, silent=True)

        metabolite_collection = db.collection(_HARMONIZED_METABOLITE_COLLECTION)
        member_collection = db.collection(_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION)
        metabolite_batch = []
        member_batch = []
        for rank, members in enumerate(groups, start=1):
            generic_summary = _generic_structure_clique_summary(
                members,
                generic_structure_classifications or {},
            )
            member_digest = _canonical_json_digest({"members": members}, 16)
            harmonized_key = f"{stage_key}-{rank:06d}-{member_digest}"
            harmonized_id = f"HarmonizedMetabolite:{harmonized_key}"
            metabolite_batch.append({
                "_key": harmonized_key,
                "id": harmonized_id,
                "stage_key": stage_key,
                "stage_id": stage_doc["id"],
                "name": f"{stage_doc['name']} group {rank}",
                "size": len(members),
                "rank_by_size": rank,
                "representative_id": members[0],
                "sample_member_ids": members[:25],
                "member_hash": member_digest,
                "generic_structure_status": generic_summary["status"],
                "generic_member_count": generic_summary["generic_member_count"],
                "specific_member_count": generic_summary["specific_member_count"],
                "unknown_member_count": generic_summary["unknown_member_count"],
                "generic_member_samples": generic_summary["generic_member_samples"],
                "specific_member_samples": generic_summary["specific_member_samples"],
                "unknown_member_samples": generic_summary["unknown_member_samples"],
            })
            for member_index, member_id in enumerate(members):
                member_edge_key = f"{harmonized_key}-{_canonical_json_digest({'member': member_id})}"
                member_batch.append({
                    "_key": member_edge_key,
                    "_from": f"{_HARMONIZED_METABOLITE_COLLECTION}/{harmonized_key}",
                    "_to": identifier_handles[member_id],
                    "id": f"HarmonizedMetaboliteMemberEdge:{member_edge_key}",
                    "stage_key": stage_key,
                    "stage_id": stage_doc["id"],
                    "harmonized_metabolite_id": harmonized_id,
                    "member_id": member_id,
                    "member_index": member_index,
                })
            if len(metabolite_batch) >= 5000:
                metabolite_collection.insert_many(metabolite_batch, overwrite=True, silent=True)
                metabolite_batch = []
            if len(member_batch) >= 10000:
                member_collection.insert_many(member_batch, overwrite=True, silent=True)
                member_batch = []
        if metabolite_batch:
            metabolite_collection.insert_many(metabolite_batch, overwrite=True, silent=True)
        for batch in _chunked_records(member_batch):
            member_collection.insert_many(batch, overwrite=True, silent=True)
    except Exception as exc:
        failed_at = datetime.now(timezone.utc).isoformat()
        failure_update = {
            "_key": stage_key,
            "status": "failed",
            "updated_at": failed_at,
            "materialization_failed_at": failed_at,
            "materialization_error": str(exc),
        }
        try:
            stage_collection.update(failure_update)
        except Exception:
            pass
        stage_doc.update(failure_update)
        raise

    completed_at = datetime.now(timezone.utc).isoformat()
    completed_update = {
        "_key": stage_key,
        "status": "complete",
        "updated_at": completed_at,
        "materialization_completed_at": completed_at,
    }
    stage_collection.update(completed_update)
    stage_doc.update(completed_update)


def _record_properties_apply_to_generic_rule(rule_ids: List[str]) -> bool:
    """Curated record properties affect the generic rule only when applied earlier."""
    if "apply_curations" not in rule_ids or "ignore_generic_structure_mismatch" not in rule_ids:
        return False
    return rule_ids.index("apply_curations") < rule_ids.index("ignore_generic_structure_mismatch")


def _generic_structure_classification_fingerprint(
    classifications: dict[str, Optional[bool]],
) -> str:
    digest = hashlib.sha256()
    for identifier, value in sorted(classifications.items()):
        digest.update(identifier.encode("utf-8"))
        digest.update(b"\0")
        digest.update(b"true" if value is True else b"false" if value is False else b"null")
        digest.update(b"\n")
    return digest.hexdigest()


def _ensure_harmonization_stage(
    db,
    rule_ids: List[str],
    rule_parameters: dict,
    graph_fingerprint: dict,
    display_name: str,
    stage_index: int,
    mass_values_provider=None,
    generic_structure_classifications_provider=None,
    curation_state=None,
    record_overlay_mass_values_provider=None,
) -> dict:
    apply_curations_enabled = "apply_curations" in rule_ids
    selected_curation_types = set(
        rule_parameters.get("apply_curations", {}).get("curation_types", [])
    ) if apply_curations_enabled else set()
    unsupported_curation_types = selected_curation_types - {
        CHEBI_RECORD_PROPERTIES,
        METABOLITE_EQUIVALENCE_EDGES,
        METABOLITE_RECORD_PROPERTIES,
        METABOLITE_RECORD_SUPPRESSIONS,
    }
    if unsupported_curation_types:
        raise ValueError(
            "Apply curations does not support: " + ", ".join(sorted(unsupported_curation_types))
        )
    apply_edge_curations = METABOLITE_EQUIVALENCE_EDGES in selected_curation_types
    apply_record_suppressions = METABOLITE_RECORD_SUPPRESSIONS in selected_curation_types
    apply_record_properties = METABOLITE_RECORD_PROPERTIES in selected_curation_types
    apply_chebi_record_properties = CHEBI_RECORD_PROPERTIES in selected_curation_types
    assertion_fallback_enabled = "force_expected_clique_assertions" in rule_ids
    if apply_curations_enabled or assertion_fallback_enabled:
        if curation_state is None:
            curation_state = _load_metabolite_edge_removal_curations()
    else:
        curation_state = {
            "pairs": set(),
            "batch_ids": [],
            "fingerprint": None,
            "assertions": [],
            "assertion_batch_ids": [],
            "assertion_fingerprint": None,
            "suppressed_identifier_ids": set(),
            "record_states": {},
            "record_decisions": {},
            "suppression_batch_ids": [],
            "suppression_fingerprint": None,
            "record_property_snapshot": None,
            "record_property_batch_ids": [],
            "record_property_fingerprint": None,
            "chebi_record_property_snapshot": None,
            "chebi_record_property_batch_ids": [],
            "chebi_record_property_fingerprint": None,
            "snapshots": {},
            "prefix": None,
        }
    stage_fingerprint = dict(graph_fingerprint)
    if "ignore_generic_structure_mismatch" in rule_ids:
        stage_fingerprint["generic_structure_classifier_version"] = (
            GENERIC_STRUCTURE_CLASSIFIER_VERSION
        )
    if apply_curations_enabled:
        stage_fingerprint["curation_fingerprint"] = payload_sha256({
            curation_type: (
                curation_state.get("edge_fingerprint")
                if curation_type == METABOLITE_EQUIVALENCE_EDGES
                else (
                    curation_state.get("suppression_fingerprint")
                    if curation_type == METABOLITE_RECORD_SUPPRESSIONS
                    else (
                        curation_state.get("record_property_fingerprint")
                        if curation_type == METABOLITE_RECORD_PROPERTIES
                        else curation_state.get("chebi_record_property_fingerprint")
                    )
                )
            )
            for curation_type in sorted(selected_curation_types)
        })
    if assertion_fallback_enabled:
        stage_fingerprint["assertion_fingerprint"] = curation_state["assertion_fingerprint"]
    baseline_generic_classifications = dict(
        generic_structure_classifications_provider()
        if generic_structure_classifications_provider is not None
        else _load_generic_structure_classifications(db)
    )
    generic_structure_classification_fingerprint = _generic_structure_classification_fingerprint(
        baseline_generic_classifications
    )
    stage_fingerprint["generic_structure_classification_fingerprint"] = (
        generic_structure_classification_fingerprint
    )
    stage_key = _harmonization_stage_key(rule_ids, rule_parameters, stage_fingerprint)
    stage_collection = db.collection(_HARMONIZATION_STAGE_COLLECTION)
    existing = stage_collection.get(stage_key)
    if existing and existing.get("status") == "complete":
        return existing
    if existing and existing.get("status") == "deleting":
        # Cleanup deletes the stage document last. A pipeline asking for the
        # same deterministic stage resumes that interrupted retirement before
        # rebuilding it from scratch.
        _delete_harmonization_stage_artifacts(db, stage_key)
        current = stage_collection.get(stage_key)
        if current and current.get("status") == "deleting":
            stage_collection.delete(stage_key)
        existing = None

    record_overlays = {}
    record_property_reports = []
    chebi_record_overlays = {}
    chebi_record_property_reports = []
    if apply_record_properties:
        record_overlays, record_property_reports = (
            _load_metabolite_identifier_record_overlays(
                db, curation_state["record_property_snapshot"]
            )
        )
    if apply_chebi_record_properties:
        chebi_record_overlays, chebi_record_property_reports = _load_record_overlays(
            db,
            curation_state["chebi_record_property_snapshot"],
            curation_set="chebi",
            model_type="ChemicalEntity",
        )
    record_generic_override_ids = {
        report["target"]["id"]
        for report in record_property_reports
        if report.get("path") == ["is_generic_structure"]
        and report.get("action") == "set"
    }
    record_generic_overrides = {
        identifier: record_overlays[identifier].get("is_generic_structure")
        for identifier in record_generic_override_ids
    }
    created_at = datetime.now(timezone.utc).isoformat()
    support_by_id = _load_metabolite_identifier_source_support(db)
    wikipathways_ignored_source_fields = (
        rule_parameters.get("ignore_wikipathways_prefixes", {}).get("source_fields", [])
        if "ignore_wikipathways_prefixes" in rule_ids
        else []
    )
    wikipathways_ignored_source_field_ids = (
        _load_wikipathways_xref_only_ids_for_ignored_source_fields(
            db,
            wikipathways_ignored_source_fields,
        )
        if wikipathways_ignored_source_fields
        else set()
    )
    filtered_support_by_id = _filter_identifier_support_for_rules(
        support_by_id,
        rule_ids,
        rule_parameters,
        wikipathways_ignored_source_field_ids,
    )
    active_ids = {
        identifier
        for identifier, sources in filtered_support_by_id.items()
        if sources
    }
    suppressed_identifier_ids = (
        set(curation_state.get("suppressed_identifier_ids") or set())
        if apply_record_suppressions
        else set()
    )
    applied_suppressed_identifier_ids = active_ids & suppressed_identifier_ids
    active_ids -= suppressed_identifier_ids
    effective_generic_classifications = dict(baseline_generic_classifications)
    if apply_record_properties:
        effective_generic_classifications.update(record_generic_overrides)
    pruning_generic_classifications = dict(baseline_generic_classifications)
    record_properties_apply_to_generic_rule = (
        apply_record_properties
        and _record_properties_apply_to_generic_rule(rule_ids)
    )
    if record_properties_apply_to_generic_rule:
        pruning_generic_classifications.update(record_generic_overrides)
    active_edges, edge_summary = _active_metabolite_identifier_mapping_edges_for_rules(
        db,
        active_ids,
        rule_ids,
        rule_parameters,
        curation_state["pairs"] if apply_edge_curations else set(),
        pruning_generic_classifications,
    )
    for edge in active_edges:
        active_ids.add(edge["start_id"])
        active_ids.add(edge["end_id"])
    assertion_edge_summary = {
        "expected_clique_assertion_count": 0,
        "expected_clique_assertion_applied_count": 0,
        "expected_clique_assertion_skipped_count": 0,
        "expected_clique_synthetic_edge_count": 0,
        "expected_clique_assertion_skipped_samples": [],
    }
    free_anomer_summary = {
        "free_anomer_candidate_pair_count": 0,
        "free_anomer_accepted_pair_count": 0,
        "free_anomer_rejected_pair_count": 0,
        "free_anomer_rejection_reason_counts": {},
        "free_anomer_candidate_identifier_count": 0,
        "free_anomer_synthetic_edge_count": 0,
    }
    if FREE_ANOMER_RULE_ID in rule_ids:
        free_anomer_edges, free_anomer_summary = _free_anomer_rule_edges(
            db,
            active_ids,
        )
        active_edges.extend(free_anomer_edges)
    if assertion_fallback_enabled:
        assertion_edges, assertion_edge_summary = _expected_clique_assertion_edges(
            active_ids,
            curation_state.get("assertions") or [],
        )
        active_edges.extend(assertion_edges)
    groups, merge_summary = _build_harmonized_groups(
        db,
        active_ids,
        active_edges,
        rule_ids,
        rule_parameters,
        record_overlays if apply_record_properties else None,
    )
    enrichment_only_summary = {
        "enrichment_only_removed_metabolite_count": 0,
        "enrichment_only_removed_group_count": 0,
        "enrichment_only_removed_singleton_count": 0,
        "enrichment_only_removed_identifier_count": 0,
        "enrichment_only_removed_edge_count": 0,
    }
    if "remove_enrichment_only_metabolites" in rule_ids:
        bio_connected_ids = _load_bio_connected_ids(db)
        active_ids, active_edges, groups, enrichment_only_summary = _remove_enrichment_only_metabolites(
            active_ids,
            active_edges,
            groups,
            filtered_support_by_id,
            bio_connected_ids,
        )
    non_singleton_member_count = sum(len(members) for members in groups)
    singleton_count = max(len(active_ids) - non_singleton_member_count, 0)
    summary = {
        **graph_fingerprint,
        **edge_summary,
        **merge_summary,
        **assertion_edge_summary,
        **free_anomer_summary,
        **enrichment_only_summary,
        "wikipathways_ignored_source_field_identifier_count": len(
            wikipathways_ignored_source_field_ids
        ),
        "active_identifier_count": len(active_ids),
        "suppressed_identifier_count": len(applied_suppressed_identifier_ids),
        "active_edge_count": len(active_edges),
        "harmonized_metabolite_count": len(groups),
        "non_singleton_clique_count": len(groups),
        "singleton_identifier_count": singleton_count,
        "clique_count": len(groups) + singleton_count,
        "largest_clique_sizes": [len(members) for members in groups[:10]],
    }
    if apply_record_properties or apply_chebi_record_properties:
        mass_values_by_id = (
            record_overlay_mass_values_provider(
                record_overlays,
                chebi_record_overlays,
            )
            if record_overlay_mass_values_provider is not None
            else _load_metabolite_identifier_mass_values(
                db,
                record_overlays,
                chebi_record_overlays,
            )
        )
    else:
        mass_values_by_id = mass_values_provider() if mass_values_provider is not None else _load_metabolite_identifier_mass_values(db)
    mw_validation = _build_harmonization_stage_mw_validation(groups, mass_values_by_id)
    denylist_validation = _build_harmonization_stage_denylist_validation(
        groups,
        curation_state["pairs"],
        apply_edge_curations,
        suppressed_identifier_ids=applied_suppressed_identifier_ids,
    )
    generic_structure_validation = _build_harmonization_stage_generic_structure_validation(
        groups,
        active_edges,
        effective_generic_classifications,
        "ignore_generic_structure_mismatch" in rule_ids,
    )
    summary.update({
        "generic_structure_inconsistent_clique_count": (
            generic_structure_validation["inconsistent_clique_count"]
        ),
        "generic_structure_unknown_clique_count": (
            generic_structure_validation["generic_unknown_clique_count"]
        ),
    })
    stage_doc = {
        "_key": stage_key,
        "id": f"HarmonizationStage:{stage_key}",
        "name": display_name,
        "created_at": created_at,
        "updated_at": created_at,
        "status": "complete",
        "graph_database": "metabolite_harmonization",
        "engine_version": _HARMONIZATION_ENGINE_VERSION,
        "generic_structure_classifier_version": (
            GENERIC_STRUCTURE_CLASSIFIER_VERSION
            if "ignore_generic_structure_mismatch" in rule_ids
            else None
        ),
        "code_fingerprint": _harmonization_code_fingerprint(),
        "generic_structure_classification_fingerprint": generic_structure_classification_fingerprint,
        "curation_types": sorted(selected_curation_types),
        "curation_fingerprint": stage_fingerprint.get("curation_fingerprint"),
        "curation_batch_ids": curation_state["batch_ids"] if apply_edge_curations else [],
        "suppression_curation_batch_ids": (
            curation_state.get("suppression_batch_ids") or []
            if apply_record_suppressions else []
        ),
        "suppression_curation_fingerprint": (
            curation_state.get("suppression_fingerprint") if apply_record_suppressions else None
        ),
        "record_property_curation_batch_ids": (
            curation_state.get("record_property_batch_ids") or []
            if apply_record_properties else []
        ),
        "record_property_curation_fingerprint": (
            curation_state.get("record_property_fingerprint")
            if apply_record_properties else None
        ),
        "record_property_curation_count": (
            sum(
                report.get("action") == "set"
                for report in record_property_reports
            ) if apply_record_properties else 0
        ),
        "chebi_record_property_curation_batch_ids": (
            curation_state.get("chebi_record_property_batch_ids") or []
            if apply_chebi_record_properties else []
        ),
        "chebi_record_property_curation_fingerprint": (
            curation_state.get("chebi_record_property_fingerprint")
            if apply_chebi_record_properties else None
        ),
        "chebi_record_property_curation_count": (
            sum(
                report.get("action") == "set"
                for report in chebi_record_property_reports
            ) if apply_chebi_record_properties else 0
        ),
        "curation_snapshots": {
            curation_type: curation_state.get("snapshots", {}).get(curation_type)
            for curation_type in sorted(selected_curation_types)
        },
        "assertion_fingerprint": curation_state.get("assertion_fingerprint"),
        "assertion_batch_ids": curation_state.get("assertion_batch_ids") or [],
        "curation_prefix": curation_state["prefix"],
        "stage_index": stage_index,
        "rule_ids": rule_ids,
        "rule_parameters": rule_parameters,
        "rules": _metabolite_rule_metadata(rule_ids, rule_parameters),
        "summary": summary,
        "validation": {
            "mw_spread": mw_validation,
            "denylist_still_merged": denylist_validation,
            "generic_structure_consistency": generic_structure_validation,
        },
        "materialization": {
            "active_identifier_chunks": _HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION,
            "evidence_edges": _HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION,
            "harmonized_metabolites": _HARMONIZED_METABOLITE_COLLECTION,
            "member_edges": _HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION,
            "singleton_note": "Singleton active identifiers are stored in active identifier chunks, not as HarmonizedMetabolite nodes.",
        },
    }
    _materialize_harmonization_stage(
        db,
        stage_key,
        stage_doc,
        active_ids,
        active_edges,
        groups,
        effective_generic_classifications,
    )
    return stage_doc


def _stage_display_label(stage: dict, rule_parameters: Optional[dict] = None) -> str:
    rule_ids = stage.get("rule_ids") or []
    if not rule_ids:
        return "Baseline"
    rule_metadata = _metabolite_rule_metadata([rule_ids[-1]], rule_parameters or {})
    if rule_metadata:
        return rule_metadata[0]["label"]
    return rule_ids[-1]


def _list_harmonization_pipelines(
    limit: int = 25,
    pipeline_key: Optional[str] = None,
) -> List[dict]:
    db = get_db("metabolite_harmonization")
    if not db.has_collection(_HARMONIZATION_PIPELINE_COLLECTION):
        return []
    pipeline_filter = "FILTER p._key == @pipeline_key" if pipeline_key else ""
    bind_vars = {"limit": limit}
    if pipeline_key:
        bind_vars["pipeline_key"] = pipeline_key
    pipelines = list(db.aql.execute(
        f"""
        FOR p IN {_HARMONIZATION_PIPELINE_COLLECTION}
          {pipeline_filter}
          SORT p.updated_at DESC
          LIMIT @limit
          LET runs = (
            FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
              FILTER r.pipeline_key == p._key
              SORT r.created_at DESC
              LIMIT 5
              RETURN r
          )
          RETURN MERGE(p, {{runs: runs}})
        """,
        bind_vars=bind_vars,
        max_runtime=120,
    ))
    for pipeline in pipelines:
        rule_parameters = pipeline.get("rule_parameters") or {}
        for run in pipeline.get("runs") or []:
            stages = run.get("stages") or []
            stages.sort(key=_stage_sort_key)
            for stage in stages:
                stage["display_label"] = _stage_display_label(stage, rule_parameters)
    return pipelines


def _build_harmonization_pipeline_tree(
    pipelines: List[dict],
    pipeline_stage_stats: Optional[List[dict]] = None,
) -> dict:
    """Build a trie of the latest run's exact stage sequence for each pipeline."""
    root = {"children": {}}
    represented_pipeline_count = 0
    stage_stats_by_key = {
        stage["_key"]: stage.get("overview_stats") or {}
        for pipeline_stats in pipeline_stage_stats or []
        for stage in pipeline_stats.get("stages") or []
        if stage.get("_key")
    }

    for pipeline in pipelines:
        latest_run = (pipeline.get("runs") or [None])[0]
        stages = sorted((latest_run or {}).get("stages") or [], key=_stage_sort_key)
        if not stages:
            continue
        represented_pipeline_count += 1
        pipeline_key = pipeline.get("_key") or pipeline.get("name") or "pipeline"
        current = root
        for stage in stages:
            stage_key = stage.get("_key")
            if not stage_key:
                continue
            current = current["children"].setdefault(stage_key, {
                "stage_key": stage_key,
                "label": stage.get("display_label") or _stage_display_label(
                    stage,
                    pipeline.get("rule_parameters") or {},
                ),
                "stage_index": stage.get("stage_index"),
                "stats": stage_stats_by_key.get(stage_key),
                "pipeline_keys": set(),
                "terminals": [],
                "children": {},
            })
            current["pipeline_keys"].add(pipeline_key)
        if current is not root:
            current["terminals"].append({
                "pipeline_key": pipeline_key,
                "pipeline_name": pipeline.get("name") or pipeline_key,
                "run_key": latest_run.get("_key"),
                "status": latest_run.get("status"),
            })

    def render_node(node: dict) -> dict:
        children = sorted(
            (render_node(child) for child in node["children"].values()),
            key=lambda child: (
                child["stage_index"] if child["stage_index"] is not None else 10**9,
                child["label"],
                child["stage_key"],
            ),
        )
        pipeline_keys = sorted(node["pipeline_keys"])
        return {
            "stage_key": node["stage_key"],
            "label": node["label"],
            "stage_index": node["stage_index"],
            "stats": node.get("stats"),
            "is_missing": node.get("stats") is None,
            "pipeline_keys": pipeline_keys,
            "pipeline_count": len(pipeline_keys),
            "is_shared": len(pipeline_keys) > 1,
            "terminals": sorted(
                node["terminals"],
                key=lambda terminal: (terminal["pipeline_name"], terminal["pipeline_key"]),
            ),
            "children": children,
        }

    roots = sorted(
        (render_node(node) for node in root["children"].values()),
        key=lambda node: (
            node["stage_index"] if node["stage_index"] is not None else 10**9,
            node["label"],
            node["stage_key"],
        ),
    )
    return {
        "roots": roots,
        "pipeline_count": len(pipelines),
        "represented_pipeline_count": represented_pipeline_count,
        "not_run_count": len(pipelines) - represented_pipeline_count,
    }


def _pipeline_selected_curation_types(pipeline: dict) -> set[str]:
    rule_ids = pipeline.get("rule_ids") or []
    if "apply_curations" in rule_ids:
        configured = (
            (pipeline.get("rule_parameters") or {})
            .get("apply_curations", {})
            .get("curation_types", [])
        )
        return set(configured or [])
    if "ignore_ramp_mapping_denylist" in rule_ids:
        return {METABOLITE_EQUIVALENCE_EDGES}
    return set()


def _pipeline_curation_fingerprint(pipeline: dict, curation_state: dict) -> Optional[str]:
    selected_types = _pipeline_selected_curation_types(pipeline)
    if not selected_types:
        return None
    fingerprints = {
        CHEBI_RECORD_PROPERTIES: curation_state.get("chebi_record_property_fingerprint"),
        METABOLITE_EQUIVALENCE_EDGES: curation_state.get("edge_fingerprint"),
        METABOLITE_RECORD_PROPERTIES: curation_state.get("record_property_fingerprint"),
        METABOLITE_RECORD_SUPPRESSIONS: curation_state.get("suppression_fingerprint"),
    }
    return payload_sha256({
        curation_type: fingerprints.get(curation_type)
        for curation_type in sorted(selected_types)
    })


def _annotate_harmonization_pipeline_curation_status(
    pipelines: List[dict],
    current_curation_fingerprint: Any,
    current_assertion_fingerprint: Optional[str] = None,
    available_stage_keys: Optional[set[str]] = None,
    current_source_graph_fingerprint: Optional[str] = None,
) -> None:
    for pipeline in pipelines:
        rule_ids = pipeline.get("rule_ids") or []
        uses_edge_curations = bool(
            {"apply_curations", "ignore_ramp_mapping_denylist"} & set(rule_ids)
        )
        uses_assertion_curations = "force_expected_clique_assertions" in rule_ids
        uses_curations = uses_edge_curations or uses_assertion_curations
        latest_complete_run = next(
            (run for run in pipeline.get("runs") or [] if run.get("status") == "complete"),
            None,
        )
        latest_run = (pipeline.get("runs") or [None])[0]
        expected_stage_keys = list(dict.fromkeys(
            (latest_complete_run.get("stage_keys") or [
                stage.get("_key")
                for stage in latest_complete_run.get("stages") or []
                if stage.get("_key")
            ]) if latest_complete_run else []
        ))
        missing_stage_keys = (
            [
                stage_key for stage_key in expected_stage_keys
                if stage_key not in available_stage_keys
            ]
            if available_stage_keys is not None
            and latest_run is latest_complete_run
            else []
        )
        completed_source_graph_fingerprint = (
            latest_complete_run.get("source_graph_content_fingerprint")
            if latest_complete_run else None
        )
        if latest_complete_run and not completed_source_graph_fingerprint:
            completed_stages = latest_complete_run.get("stages") or []
            completed_source_graph_fingerprint = next((
                (stage.get("summary") or {}).get("source_graph_content_fingerprint")
                for stage in completed_stages
                if (stage.get("summary") or {}).get("source_graph_content_fingerprint")
            ), None)
        source_graph_changed = bool(
            latest_complete_run
            and current_source_graph_fingerprint
            and completed_source_graph_fingerprint != current_source_graph_fingerprint
        )
        completed_engine_version = (
            latest_complete_run.get("engine_version")
            if latest_complete_run
            else None
        ) or pipeline.get("engine_version")
        engine_version_changed = bool(
            latest_complete_run
            and completed_engine_version != _HARMONIZATION_ENGINE_VERSION
        )
        pipeline_curation_fingerprint = (
            _pipeline_curation_fingerprint(pipeline, current_curation_fingerprint)
            if isinstance(current_curation_fingerprint, dict)
            else current_curation_fingerprint
        )
        edge_curations_changed = bool(
            uses_edge_curations
            and latest_complete_run
            and latest_complete_run.get("curation_fingerprint") != pipeline_curation_fingerprint
        )
        assertion_curations_changed = bool(
            uses_assertion_curations
            and latest_complete_run
            and latest_complete_run.get("assertion_fingerprint") != current_assertion_fingerprint
        )
        needs_sync = edge_curations_changed or assertion_curations_changed
        changed_rule_ids = set()
        if edge_curations_changed:
            changed_rule_ids.add(
                "apply_curations"
                if "apply_curations" in rule_ids
                else "ignore_ramp_mapping_denylist"
            )
        if assertion_curations_changed:
            changed_rule_ids.add("force_expected_clique_assertions")
        changed_rule_indexes = [
            index
            for index, rule_id in enumerate(rule_ids, start=1)
            if rule_id in changed_rule_ids
        ]
        pipeline["uses_curations"] = uses_curations
        pipeline["uses_edge_curations"] = uses_edge_curations
        pipeline["uses_assertion_curations"] = uses_assertion_curations
        pipeline["edge_curations_changed"] = edge_curations_changed
        pipeline["assertion_curations_changed"] = assertion_curations_changed
        pipeline["needs_curation_sync"] = needs_sync
        pipeline["missing_stage_keys"] = missing_stage_keys
        pipeline["missing_stage_count"] = len(missing_stage_keys)
        pipeline["has_stage_integrity_failure"] = bool(missing_stage_keys)
        pipeline["source_graph_changed"] = source_graph_changed
        pipeline["completed_source_graph_fingerprint"] = completed_source_graph_fingerprint
        pipeline["current_source_graph_fingerprint"] = current_source_graph_fingerprint
        pipeline["engine_version_changed"] = engine_version_changed
        pipeline["completed_engine_version"] = completed_engine_version
        pipeline["current_engine_version"] = _HARMONIZATION_ENGINE_VERSION
        pipeline["sync_from_stage_index"] = min(changed_rule_indexes) if changed_rule_indexes else None
        pipeline["current_curation_fingerprint"] = pipeline_curation_fingerprint if uses_edge_curations else None
        pipeline["current_assertion_fingerprint"] = current_assertion_fingerprint if uses_assertion_curations else None
        pipeline["latest_complete_run_key"] = latest_complete_run.get("_key") if latest_complete_run else None
        pipeline["run_action_enabled"] = True
        if latest_run and latest_run.get("status") == "failed":
            pipeline["run_button_label"] = "Retry failed pipeline"
            pipeline["run_action_kind"] = "retry"
        elif latest_run and latest_run.get("status") in {"queued", "running", "cleaning_up"}:
            pipeline["run_button_label"] = "Running..."
            pipeline["run_action_kind"] = "running"
            pipeline["run_action_enabled"] = False
        elif latest_complete_run is None:
            pipeline["run_button_label"] = "Run pipeline"
            pipeline["run_action_kind"] = "run"
        elif missing_stage_keys:
            pipeline["run_button_label"] = "Repair pipeline"
            pipeline["run_action_kind"] = "repair"
        elif source_graph_changed:
            pipeline["run_button_label"] = "Re-run pipeline"
            pipeline["run_action_kind"] = "rebuild"
        elif engine_version_changed:
            pipeline["run_button_label"] = "Re-run pipeline"
            pipeline["run_action_kind"] = "rebuild"
        elif needs_sync:
            pipeline["run_button_label"] = "Sync curations"
            pipeline["run_action_kind"] = "sync"
        else:
            pipeline["run_button_label"] = "Up to date"
            pipeline["run_action_kind"] = "current"
            pipeline["run_action_enabled"] = False


def _list_harmonization_stages(limit: int = 50) -> List[dict]:
    db = get_db("metabolite_harmonization")
    if not db.has_collection(_HARMONIZATION_STAGE_COLLECTION):
        return []
    return list(db.aql.execute(
        f"""
        LET recent_stages = (
          FOR s IN {_HARMONIZATION_STAGE_COLLECTION}
            FILTER s.status == "complete"
            SORT s.created_at DESC
            LIMIT @limit
            RETURN s
        )
        FOR s IN recent_stages
          SORT s.created_at ASC
          RETURN s
        """,
        bind_vars={"limit": limit},
        max_runtime=120,
    ))


def _harmonization_jobs_by_pipeline_key(jobs: List[dict]) -> Dict[str, List[dict]]:
    jobs_by_pipeline_key: Dict[str, List[dict]] = {}
    for job in jobs:
        pipeline_key = job.get("pipeline_key")
        if not pipeline_key:
            continue
        jobs_by_pipeline_key.setdefault(pipeline_key, []).append(job)
    for pipeline_key, pipeline_jobs in jobs_by_pipeline_key.items():
        pipeline_jobs.sort(key=lambda job: job.get("created_at") or "", reverse=True)
        successful_actions = set()
        visible_jobs = []
        for job in pipeline_jobs:
            action = job.get("action")
            if job.get("status") == "failed" and action in successful_actions:
                continue
            visible_jobs.append(job)
            if job.get("status") == "complete":
                successful_actions.add(action)
        jobs_by_pipeline_key[pipeline_key] = visible_jobs
    return jobs_by_pipeline_key


def _annotate_harmonization_pipeline_run_progress(
    pipelines: List[dict],
    jobs: List[dict],
) -> None:
    pipelines_by_key = {pipeline.get("_key"): pipeline for pipeline in pipelines}
    for pipeline in pipelines:
        pipeline["active_stage_index"] = None
        pipeline["active_stage_text"] = None
        active_run = next(
            (
                run
                for run in pipeline.get("runs") or []
                if run.get("status") in {"running", "cleaning_up"}
            ),
            None,
        )
        if active_run is None:
            continue
        if active_run.get("status") == "cleaning_up":
            pipeline["active_stage_text"] = "Cleaning up old snapshots"
            continue
        completed_stage_count = len(active_run.get("stage_keys") or [])
        rules = pipeline.get("rules") or []
        if completed_stage_count == 0:
            pipeline["active_stage_index"] = 0
            pipeline["active_stage_text"] = "Baseline"
        elif completed_stage_count <= len(rules):
            active_rule = rules[completed_stage_count - 1]
            active_label = active_rule.get("label") or active_rule.get("id")
            pipeline["active_stage_index"] = completed_stage_count
            pipeline["active_stage_text"] = f"Step {completed_stage_count}: {active_label}"
        else:
            pipeline["active_stage_text"] = "Cleaning up old snapshots"

    for job in jobs:
        if job.get("action") != "run_pipeline" or job.get("status") != "running":
            continue
        pipeline = pipelines_by_key.get(job.get("pipeline_key"))
        if pipeline:
            job["active_stage_text"] = pipeline.get("active_stage_text")


def _list_harmonization_stage_overview_stats(
    limit: int = 100,
    stage_keys: Optional[List[str]] = None,
    include_distribution: bool = True,
) -> List[dict]:
    db = get_db("metabolite_harmonization")
    if not db.has_collection(_HARMONIZATION_STAGE_COLLECTION):
        return []
    unique_stage_keys = list(dict.fromkeys(stage_keys or []))
    if unique_stage_keys:
        stage_filter = 'FILTER s._key IN @stage_keys AND s.status == "complete"'
        bind_vars = {"limit": max(limit, len(unique_stage_keys)), "stage_keys": unique_stage_keys}
    else:
        stage_filter = "FILTER s.status == \"complete\""
        bind_vars = {"limit": limit, "stage_keys": []}
    stages = list(db.aql.execute(
        f"""
        FOR s IN {_HARMONIZATION_STAGE_COLLECTION}
          {stage_filter}
          SORT s.created_at DESC
          LIMIT @limit
          LET mw = s.validation.mw_spread
          LET overview_mw = (
            IS_OBJECT(mw) AND IS_ARRAY(mw.finding_index)
            ? MERGE(
                KEEP(
                  mw,
                  "computed", "warning_count", "error_count",
                  "review_warning_count", "accepted_count",
                  "stale_acceptance_count"
                ),
                {{
                  finding_index: (
                    FOR finding IN mw.finding_index
                      RETURN KEEP(
                        finding,
                        "finding_id", "member_ids", "severity",
                        "evidence_fingerprint", "legacy_evidence_fingerprint"
                      )
                  )
                }}
              )
            : mw
          )
          RETURN MERGE(
            KEEP(s, "_key", "name", "created_at", "stage_index", "rule_ids", "summary"),
            {{
              validation: {{
                mw_spread: overview_mw,
                denylist_still_merged: UNSET(
                  s.validation.denylist_still_merged || {{}}, "warnings"
                ),
                generic_structure_consistency: UNSET(
                  s.validation.generic_structure_consistency || {{}}, "warnings"
                )
              }}
            }}
          )
        """,
        bind_vars=bind_vars,
        max_runtime=120,
    ))
    if not stages:
        return []
    stage_keys = [stage["_key"] for stage in stages]
    distribution_by_stage = {}
    if include_distribution and db.has_collection(_HARMONIZED_METABOLITE_COLLECTION):
        distribution_rows = list(db.aql.execute(
            f"""
            FOR c IN {_HARMONIZED_METABOLITE_COLLECTION}
              FILTER c.stage_key IN @stage_keys
              COLLECT stage_key = c.stage_key
              AGGREGATE
                non_singleton_groups = COUNT(),
                non_singleton_members = SUM(c.size),
                max_size = MAX(c.size),
                average_size = AVG(c.size),
                size_2 = SUM(c.size == 2 ? 1 : 0),
                size_3_5 = SUM(c.size >= 3 AND c.size <= 5 ? 1 : 0),
                size_6_10 = SUM(c.size >= 6 AND c.size <= 10 ? 1 : 0),
                size_11_50 = SUM(c.size >= 11 AND c.size <= 50 ? 1 : 0),
                size_51_100 = SUM(c.size >= 51 AND c.size <= 100 ? 1 : 0),
                size_101_plus = SUM(c.size >= 101 ? 1 : 0)
              RETURN {{
                stage_key: stage_key,
                non_singleton_groups: non_singleton_groups,
                non_singleton_members: non_singleton_members,
                max_size: max_size,
                average_size: average_size,
                bins: {{
                  size_2: size_2,
                  size_3_5: size_3_5,
                  size_6_10: size_6_10,
                  size_11_50: size_11_50,
                  size_51_100: size_51_100,
                  size_101_plus: size_101_plus
                }}
              }}
            """,
            bind_vars={"stage_keys": stage_keys},
            max_runtime=120,
        ))
        distribution_by_stage = {
            row["stage_key"]: row
            for row in distribution_rows
        }
    try:
        mw_adjudication_state = _load_metabolite_mw_adjudications()
    except Exception:
        mw_adjudication_state = None
    for stage in stages:
        summary = stage.get("summary") or {}
        distribution = distribution_by_stage.get(stage["_key"], {})
        mw_validation = _mw_validation_overview_from_doc(
            stage, mw_adjudication_state
        )
        denylist_validation = _harmonization_stage_denylist_validation_from_doc(stage)
        generic_structure_validation = _harmonization_stage_generic_structure_validation_from_doc(stage)
        stage["overview_stats"] = {
            "active_identifier_count": summary.get("active_identifier_count", 0),
            "clique_count": summary.get("clique_count", 0),
            "non_singleton_clique_count": summary.get(
                "non_singleton_clique_count",
                distribution.get("non_singleton_groups", 0),
            ),
            "singleton_identifier_count": summary.get("singleton_identifier_count", 0),
            "active_edge_count": summary.get("active_edge_count", 0),
            "ignored_edge_count": summary.get("ignored_edge_count", 0),
            "assertion_fallback_enabled": "force_expected_clique_assertions" in (stage.get("rule_ids") or []),
            "assertion_fallback_applied_count": summary.get("expected_clique_assertion_applied_count", 0),
            "assertion_fallback_skipped_count": summary.get("expected_clique_assertion_skipped_count", 0),
            "assertion_synthetic_edge_count": summary.get("expected_clique_synthetic_edge_count", 0),
            "mw_warning_count": mw_validation.get("warning_count", 0),
            "mw_error_count": mw_validation.get(
                "error_count", mw_validation.get("warning_count", 0)
            ),
            "mw_review_warning_count": mw_validation.get("review_warning_count", 0),
            "mw_accepted_count": mw_validation.get("accepted_count", 0),
            "mw_stale_acceptance_count": mw_validation.get("stale_acceptance_count", 0),
            "mw_validation_computed": mw_validation.get("computed", False),
            "denylist_warning_count": denylist_validation.get("warning_count", 0),
            "denylist_validation_computed": denylist_validation.get("computed", False),
            "denylist_rule_enabled": denylist_validation.get("rule_enabled", False),
            "generic_structure_validation_computed": generic_structure_validation.get("computed", False),
            "generic_structure_rule_enabled": generic_structure_validation.get("rule_enabled", False),
            "generic_structure_failure_count": generic_structure_validation.get(
                "inconsistent_clique_count", 0
            ),
            "generic_structure_gap_count": generic_structure_validation.get(
                "generic_unknown_clique_count", 0
            ),
            "max_size": distribution.get("max_size") or (
                (summary.get("largest_clique_sizes") or [0])[0]
            ),
            "average_non_singleton_size": distribution.get("average_size") or 0,
            "largest_clique_sizes": summary.get("largest_clique_sizes") or [],
            "bins": distribution.get("bins") or {
                "size_2": 0,
                "size_3_5": 0,
                "size_6_10": 0,
                "size_11_50": 0,
                "size_51_100": 0,
                "size_101_plus": 0,
            },
        }
    return stages


def _stage_sort_key(stage: dict) -> tuple:
    stage_index = stage.get("stage_index")
    if stage_index is None:
        stage_index = len(stage.get("rule_ids") or [])
    return (
        stage_index,
        len(stage.get("rule_ids") or []),
        stage.get("_key") or "",
    )


def _evaluate_expected_clique_assertions(
    assertions: List[dict],
    stage_keys: List[str],
    active_memberships: Iterable[dict],
    clique_memberships: Iterable[dict],
) -> dict:
    ordered_stage_keys = list(dict.fromkeys(stage_keys))
    active_by_stage: Dict[str, set[str]] = {stage_key: set() for stage_key in ordered_stage_keys}
    for row in active_memberships:
        if row.get("stage_key") in active_by_stage and row.get("member_id"):
            active_by_stage[row["stage_key"]].add(row["member_id"])

    clique_by_stage: Dict[str, Dict[str, dict]] = {stage_key: {} for stage_key in ordered_stage_keys}
    for row in clique_memberships:
        stage_key = row.get("stage_key")
        member_id = row.get("member_id")
        if stage_key in clique_by_stage and member_id:
            clique_by_stage[stage_key][member_id] = row

    by_stage = {}
    matrix_rows = [{**assertion, "stage_results": []} for assertion in assertions]
    matrix_by_id = {row["assertion_id"]: row for row in matrix_rows}
    for stage_key in ordered_stage_keys:
        results = []
        for assertion in assertions:
            member_ids = assertion.get("member_ids") or []
            missing_ids = [member_id for member_id in member_ids if member_id not in active_by_stage[stage_key]]
            partitions_by_token: Dict[str, dict] = {}
            for member_id in member_ids:
                if member_id in missing_ids:
                    continue
                membership = clique_by_stage[stage_key].get(member_id)
                if membership:
                    token = membership.get("clique_id") or f"clique:{member_id}"
                    partition = partitions_by_token.setdefault(token, {
                        "clique_id": membership.get("clique_id"),
                        "rank_by_size": membership.get("rank_by_size"),
                        "clique_size": membership.get("clique_size"),
                        "representative_id": membership.get("representative_id"),
                        "member_ids": [],
                    })
                else:
                    token = f"singleton:{member_id}"
                    partition = partitions_by_token.setdefault(token, {
                        "clique_id": None,
                        "rank_by_size": None,
                        "clique_size": 1,
                        "representative_id": member_id,
                        "member_ids": [],
                    })
                partition["member_ids"].append(member_id)
            partitions = sorted(
                partitions_by_token.values(),
                key=lambda item: (item.get("representative_id") or "", item.get("clique_id") or ""),
            )
            if missing_ids:
                status = "missing"
            elif len(partitions) == 1:
                status = "pass"
            else:
                status = "split"
            result = {
                "assertion_id": assertion["assertion_id"],
                "name": assertion.get("name") or assertion["assertion_id"],
                "rationale": assertion.get("rationale"),
                "member_ids": member_ids,
                "status": status,
                "missing_member_ids": missing_ids,
                "partitions": partitions,
                "inspection_query": urlencode({"id": " ".join(member_ids), "stages": stage_key}),
                "published_at": assertion.get("published_at"),
                "published_by": assertion.get("published_by"),
                "curation_batch_id": assertion.get("curation_batch_id"),
            }
            results.append(result)
            matrix_by_id[assertion["assertion_id"]]["stage_results"].append({
                "stage_key": stage_key,
                "status": status,
                "missing_member_ids": missing_ids,
                "partitions": partitions,
            })
        status_counts = {
            status: sum(result["status"] == status for result in results)
            for status in ("pass", "split", "missing")
        }
        by_stage[stage_key] = {
            "stage_key": stage_key,
            "assertion_count": len(results),
            "pass_count": status_counts["pass"],
            "split_count": status_counts["split"],
            "missing_count": status_counts["missing"],
            "failure_count": status_counts["split"] + status_counts["missing"],
            "results": results,
        }
    return {"by_stage": by_stage, "matrix_rows": matrix_rows}


def _load_expected_clique_assertion_results(
    db,
    assertions: List[dict],
    stage_keys: List[str],
) -> dict:
    ordered_stage_keys = list(dict.fromkeys(stage_keys))
    member_ids = sorted({
        member_id
        for assertion in assertions
        for member_id in assertion.get("member_ids") or []
    })
    if not assertions or not ordered_stage_keys or not member_ids:
        return _evaluate_expected_clique_assertions(assertions, ordered_stage_keys, [], [])

    active_rows = list(db.aql.execute(
        f"""
        FOR chunk IN {_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION}
          FILTER chunk.stage_key IN @stage_keys
          FOR member_id IN @member_ids
            FILTER member_id IN (chunk.identifier_ids || [])
            RETURN DISTINCT {{stage_key: chunk.stage_key, member_id: member_id}}
        """,
        bind_vars={"stage_keys": ordered_stage_keys, "member_ids": member_ids},
        max_runtime=120,
    ))
    clique_rows = list(db.aql.execute(
        f"""
        FOR member IN MetaboliteIdentifier
          FILTER member.id IN @member_ids
          FOR clique, edge IN INBOUND member {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
            FILTER edge.stage_key IN @stage_keys
            RETURN {{
              stage_key: edge.stage_key,
              member_id: member.id,
              clique_id: clique._id,
              rank_by_size: clique.rank_by_size,
              clique_size: clique.size,
              representative_id: clique.representative_id
            }}
        """,
        bind_vars={"stage_keys": ordered_stage_keys, "member_ids": member_ids},
        max_runtime=120,
    ))
    return _evaluate_expected_clique_assertions(
        assertions,
        ordered_stage_keys,
        active_rows,
        clique_rows,
    )


def _expected_clique_assertion_matrix(assertion_results: dict, stages: List[dict]) -> dict:
    stage_by_key = {stage["_key"]: stage for stage in stages}
    stage_keys = set(stage_by_key)
    rows = []
    for assertion in assertion_results.get("matrix_rows") or []:
        stage_results = [
            result
            for result in assertion.get("stage_results") or []
            if result.get("stage_key") in stage_keys
        ]
        first_failure = next((result for result in stage_results if result["status"] != "pass"), None)
        final_result = stage_results[-1] if stage_results else None
        rows.append({
            **assertion,
            "stage_results": stage_results,
            "first_failure": ({
                **first_failure,
                "display_label": stage_by_key.get(first_failure["stage_key"], {}).get("display_label"),
            } if first_failure else None),
            "final_status": final_result.get("status") if final_result else None,
        })
    return {
        "assertion_count": len(rows),
        "stages": stages,
        "rows": rows,
    }


def _list_harmonization_pipeline_stage_overview_stats(
    pipelines: List[dict],
    curation_state: Optional[dict] = None,
    include_distribution: bool = True,
    include_assertion_matrix: bool = True,
) -> List[dict]:
    if not pipelines:
        return []
    stage_keys = []
    for pipeline in pipelines:
        latest_run = (pipeline.get("runs") or [None])[0]
        for stage in (latest_run or {}).get("stages") or []:
            stage_key = stage.get("_key")
            if stage_key and stage_key not in stage_keys:
                stage_keys.append(stage_key)
    if not stage_keys:
        return []
    stats_by_key = {
        stage["_key"]: stage
        for stage in _list_harmonization_stage_overview_stats(
            limit=max(len(stage_keys), 100),
            stage_keys=stage_keys,
            include_distribution=include_distribution,
        )
    }
    pipeline_stats = []
    assertions = (curation_state or {}).get("assertions") or []
    assertion_results = _load_expected_clique_assertion_results(
        get_db("metabolite_harmonization"),
        assertions,
        stage_keys,
    )
    for stage_key, stage_assertions in assertion_results["by_stage"].items():
        stage = stats_by_key.get(stage_key)
        if stage:
            stage["overview_stats"].update({
                "assertion_count": stage_assertions["assertion_count"],
                "assertion_pass_count": stage_assertions["pass_count"],
                "assertion_failure_count": stage_assertions["failure_count"],
                "assertion_split_count": stage_assertions["split_count"],
                "assertion_missing_count": stage_assertions["missing_count"],
            })
    for pipeline in pipelines:
        latest_run = (pipeline.get("runs") or [None])[0]
        stages = [
            stats_by_key[stage["_key"]]
            for stage in (latest_run or {}).get("stages") or []
            if stage.get("_key") in stats_by_key
        ]
        for stage in stages:
            stage["display_label"] = _stage_display_label(stage, pipeline.get("rule_parameters") or {})
        stages.sort(key=_stage_sort_key)
        if not stages:
            continue
        fingerprints = {s.get("code_fingerprint") for s in stages if s.get("code_fingerprint")}
        pipeline_stats.append({
            "pipeline_key": pipeline.get("_key"),
            "pipeline_name": pipeline.get("name"),
            "run": latest_run,
            "stages": stages,
            "code_fingerprint_mismatch": len(fingerprints) > 1,
            "code_fingerprints": sorted(fingerprints),
            "expected_clique_assertions": (
                _expected_clique_assertion_matrix(assertion_results, stages)
                if include_assertion_matrix
                else {"assertion_count": len(assertions), "stages": [], "rows": []}
            ),
        })
    return pipeline_stats


def _get_harmonization_pipeline(db, pipeline_key: str) -> dict:
    if not db.has_collection(_HARMONIZATION_PIPELINE_COLLECTION):
        raise ValueError("Harmonization pipeline collection does not exist.")
    pipeline = db.collection(_HARMONIZATION_PIPELINE_COLLECTION).get(pipeline_key)
    if not pipeline:
        raise ValueError(f"Harmonization pipeline {pipeline_key} does not exist.")
    return pipeline


def _rename_harmonization_pipeline(pipeline_key: str, pipeline_name: str) -> dict:
    db = get_db("metabolite_harmonization")
    _ensure_harmonization_pipeline_collections(db)
    pipeline = _get_harmonization_pipeline(db, pipeline_key)
    display_name = pipeline_name.strip()
    if not display_name:
        raise ValueError("Pipeline name cannot be empty.")
    updated_at = datetime.now(timezone.utc).isoformat()
    updated = {
        "_key": pipeline_key,
        "name": display_name,
        "updated_at": updated_at,
    }
    db.collection(_HARMONIZATION_PIPELINE_COLLECTION).update(updated)
    if db.has_collection(_HARMONIZATION_PIPELINE_RUN_COLLECTION):
        db.aql.execute(
            f"""
            FOR r IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
              FILTER r.pipeline_key == @pipeline_key
              UPDATE r WITH {{pipeline_name: @pipeline_name}} IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
            """,
            bind_vars={"pipeline_key": pipeline_key, "pipeline_name": display_name},
            max_runtime=120,
        )
    return {**pipeline, **updated}


def _upsert_harmonization_pipeline(
    pipeline_key: str,
    pipeline_name: str,
    rule_ids: List[str],
    rule_parameters: dict,
) -> dict:
    db = get_db("metabolite_harmonization")
    _ensure_harmonization_pipeline_collections(db)
    normalized_rule_ids = _normalize_metabolite_rule_ids(rule_ids)
    normalized_rule_parameters = _normalize_metabolite_rule_parameters(normalized_rule_ids, rule_parameters)
    now = datetime.now(timezone.utc).isoformat()
    display_name = pipeline_name.strip() or "Untitled harmonization pipeline"
    if pipeline_key:
        existing = db.collection(_HARMONIZATION_PIPELINE_COLLECTION).get(pipeline_key)
        if not existing:
            raise ValueError(f"Harmonization pipeline {pipeline_key} does not exist.")
        created_at = existing.get("created_at") or now
    else:
        pipeline_key = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
        created_at = now
    steps = []
    for index, rule_id in enumerate(normalized_rule_ids, start=1):
        steps.append({
            "index": index,
            "rule_id": rule_id,
            "parameters": normalized_rule_parameters.get(rule_id, {}),
            "rule": _metabolite_rule_metadata([rule_id], normalized_rule_parameters)[0],
        })
    pipeline_doc = {
        "_key": pipeline_key,
        "id": f"HarmonizationPipeline:{pipeline_key}",
        "name": display_name,
        "created_at": created_at,
        "updated_at": now,
        "engine_version": _HARMONIZATION_ENGINE_VERSION,
        "rule_ids": normalized_rule_ids,
        "rule_parameters": normalized_rule_parameters,
        "rules": _metabolite_rule_metadata(normalized_rule_ids, normalized_rule_parameters),
        "steps": steps,
        "status": "ready",
    }
    db.collection(_HARMONIZATION_PIPELINE_COLLECTION).insert(pipeline_doc, overwrite=True)
    return pipeline_doc


def _harmonization_runner_instance_id() -> str:
    configured_id = os.getenv("QA_BROWSER_RUNNER_INSTANCE_ID", "").strip()
    if configured_id:
        return configured_id
    return f"qa-browser-local-{os.getenv('USER', 'unknown')}"


def _reconcile_interrupted_harmonization_runs(
    db,
    runner_instance_id: Optional[str] = None,
    runner_boot_id: Optional[str] = None,
    interrupted_at: Optional[str] = None,
) -> List[str]:
    instance_id = runner_instance_id or _harmonization_runner_instance_id()
    boot_id = runner_boot_id or _HARMONIZATION_RUNNER_BOOT_ID
    timestamp = interrupted_at or datetime.now(timezone.utc).isoformat()
    interrupted_runs = list(db.aql.execute(
        f"""
        FOR run IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
          FILTER run.status IN ["running", "cleaning_up"]
          FILTER run.runner_instance_id == @runner_instance_id
          FILTER run.runner_boot_id != @runner_boot_id
          UPDATE run WITH {{
            status: "failed",
            failure_kind: "interrupted",
            interruption_reason: "qa_browser_restarted",
            completed_at: @interrupted_at,
            interrupted_at: @interrupted_at,
            error: "QA Browser restarted before this pipeline run completed."
          }} IN {_HARMONIZATION_PIPELINE_RUN_COLLECTION}
          RETURN {{
            run_key: NEW._key,
            pipeline_key: NEW.pipeline_key
          }}
        """,
        bind_vars={
            "runner_instance_id": instance_id,
            "runner_boot_id": boot_id,
            "interrupted_at": timestamp,
        },
        max_runtime=120,
    ))
    if interrupted_runs:
        list(db.aql.execute(
            f"""
            FOR interrupted_run IN @interrupted_runs
              FOR pipeline IN {_HARMONIZATION_PIPELINE_COLLECTION}
                FILTER pipeline._key == interrupted_run.pipeline_key
                FILTER pipeline.latest_run_key == interrupted_run.run_key
                FILTER pipeline.status IN ["running", "cleaning_up"]
                UPDATE pipeline WITH {{
                  status: "failed",
                  updated_at: @interrupted_at
                }} IN {_HARMONIZATION_PIPELINE_COLLECTION}
                RETURN NEW._key
            """,
            bind_vars={
                "interrupted_runs": interrupted_runs,
                "interrupted_at": timestamp,
            },
            max_runtime=120,
        ))
    return [
        interrupted_run["run_key"]
        for interrupted_run in interrupted_runs
        if interrupted_run.get("run_key")
    ]


def _try_reconcile_interrupted_harmonization_runs(force: bool = False) -> List[str]:
    global _harmonization_reconciliation_last_attempt
    with _harmonization_reconciliation_lock:
        now_monotonic = time.monotonic()
        if (
            not force
            and now_monotonic - _harmonization_reconciliation_last_attempt
            < _HARMONIZATION_RECONCILIATION_RETRY_SECONDS
        ):
            return []
        _harmonization_reconciliation_last_attempt = now_monotonic
        try:
            return _reconcile_interrupted_harmonization_runs(
                get_db("metabolite_harmonization")
            )
        except Exception as exc:
            print(f"Warning: could not reconcile interrupted harmonization runs: {exc}")
            return []


def _run_harmonization_pipeline(pipeline_key: str) -> dict:
    db = get_db("metabolite_harmonization")
    _ensure_harmonization_pipeline_collections(db)
    pipeline = _get_harmonization_pipeline(db, pipeline_key)
    started = time.time()
    created_at = datetime.now(timezone.utc).isoformat()
    run_key = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    run_doc = {
        "_key": run_key,
        "id": f"HarmonizationPipelineRun:{run_key}",
        "pipeline_key": pipeline_key,
        "pipeline_id": pipeline["id"],
        "pipeline_name": pipeline.get("name"),
        "engine_version": _HARMONIZATION_ENGINE_VERSION,
        "runner_instance_id": _harmonization_runner_instance_id(),
        "runner_boot_id": _HARMONIZATION_RUNNER_BOOT_ID,
        "created_at": created_at,
        "started_at": created_at,
        "completed_at": None,
        "status": "running",
        "stage_keys": [],
        "stages": [],
        "error": None,
    }
    run_collection = db.collection(_HARMONIZATION_PIPELINE_RUN_COLLECTION)
    run_collection.insert(run_doc, overwrite=True)
    try:
        db.collection(_HARMONIZATION_PIPELINE_COLLECTION).update({
            "_key": pipeline_key,
            "latest_run_key": run_key,
            "updated_at": created_at,
            "status": "running",
        })
        graph_fingerprint = _harmonization_graph_fingerprint(db)
        run_collection.update({
            "_key": run_key,
            "source_graph_fingerprint": graph_fingerprint.get("source_graph_fingerprint"),
            "source_graph_content_fingerprint": graph_fingerprint.get(
                "source_graph_content_fingerprint"
            ),
        })
        stage_docs = []
        mass_values_cache = {}
        def mass_values_provider():
            if "values" not in mass_values_cache:
                mass_values_cache["values"] = _load_metabolite_identifier_mass_values(db)
            return mass_values_cache["values"]

        generic_structure_classifications_cache = {}
        def generic_structure_classifications_provider():
            if "values" not in generic_structure_classifications_cache:
                generic_structure_classifications_cache["values"] = (
                    _load_generic_structure_classifications(db)
                )
            return generic_structure_classifications_cache["values"]

        normalized_rule_ids = pipeline.get("rule_ids") or []
        normalized_rule_parameters = pipeline.get("rule_parameters") or {}
        curation_state = (
            _load_metabolite_edge_removal_curations(force_refresh=True)
            if (
                "apply_curations" in normalized_rule_ids
                or "force_expected_clique_assertions" in normalized_rule_ids
            )
            else None
        )
        record_overlay_mass_values_provider = (
            _cached_record_overlay_mass_values_provider(db)
        )
        if curation_state is not None:
            selected_curation_fingerprint = _pipeline_curation_fingerprint(
                pipeline, curation_state
            )
            run_collection.update({
                "_key": run_key,
                "curation_fingerprint": selected_curation_fingerprint,
                "curation_batch_ids": curation_state["batch_ids"],
                "assertion_fingerprint": curation_state["assertion_fingerprint"],
                "assertion_batch_ids": curation_state["assertion_batch_ids"],
                "curation_prefix": curation_state["prefix"],
            })
        cumulative_rule_ids: List[str] = []
        cumulative_rule_parameters = _cumulative_rule_parameters(cumulative_rule_ids, normalized_rule_parameters)
        baseline_stage = _ensure_harmonization_stage(
            db,
            cumulative_rule_ids,
            cumulative_rule_parameters,
            graph_fingerprint,
            f"{pipeline.get('name') or 'Pipeline'}: baseline",
            0,
            mass_values_provider,
            generic_structure_classifications_provider,
            curation_state,
            record_overlay_mass_values_provider,
        )
        stage_docs.append(baseline_stage)
        run_collection.update({
            "_key": run_key,
            "stage_keys": [stage["_key"] for stage in stage_docs],
            "stages": [
                {
                    "_key": stage["_key"],
                    "name": stage.get("name"),
                    "stage_index": stage.get("stage_index"),
                    "rule_ids": stage.get("rule_ids") or [],
                    "summary": stage.get("summary") or {},
                }
                for stage in stage_docs
            ],
        })
        for stage_index, rule_id in enumerate(normalized_rule_ids, start=1):
            cumulative_rule_ids = normalized_rule_ids[:stage_index]
            cumulative_rule_parameters = _cumulative_rule_parameters(cumulative_rule_ids, normalized_rule_parameters)
            rule_metadata = _metabolite_rule_metadata([rule_id], normalized_rule_parameters)
            rule_label = rule_metadata[0]["label"] if rule_metadata else rule_id
            stage_docs.append(_ensure_harmonization_stage(
                db,
                cumulative_rule_ids,
                cumulative_rule_parameters,
                graph_fingerprint,
                f"{pipeline.get('name') or 'Pipeline'}: {stage_index}. {rule_label}",
                stage_index,
                mass_values_provider,
                generic_structure_classifications_provider,
                curation_state,
                record_overlay_mass_values_provider,
            ))
            run_collection.update({
                "_key": run_key,
                "stage_keys": [stage["_key"] for stage in stage_docs],
                "stages": [
                    {
                        "_key": stage["_key"],
                        "name": stage.get("name"),
                        "stage_index": stage.get("stage_index"),
                        "rule_ids": stage.get("rule_ids") or [],
                        "summary": stage.get("summary") or {},
                    }
                    for stage in stage_docs
                ],
            })
        stages_completed_at = datetime.now(timezone.utc).isoformat()
        cleaning_update = {
            "_key": run_key,
            "status": "cleaning_up",
            "stages_completed_at": stages_completed_at,
            "stage_elapsed_seconds": round(time.time() - started, 1),
            "stage_keys": [stage["_key"] for stage in stage_docs],
            "stages": [
                {
                    "_key": stage["_key"],
                    "name": stage.get("name"),
                    "stage_index": stage.get("stage_index"),
                    "rule_ids": stage.get("rule_ids") or [],
                    "summary": stage.get("summary") or {},
                }
                for stage in stage_docs
            ],
        }
        run_collection.update(cleaning_update)
        db.collection(_HARMONIZATION_PIPELINE_COLLECTION).update({
            "_key": pipeline_key,
            "latest_run_key": run_key,
            "latest_stage_key": stage_docs[-1]["_key"] if stage_docs else None,
            "updated_at": stages_completed_at,
            "status": "cleaning_up",
        })
        cleanup = {
            "deleted_run_keys": [],
            "deleted_stage_keys": [],
            "deleted_superseded_stage_keys": [],
            "deleted_orphan_stage_keys": [],
        }
        cleanup_error = None
        try:
            cleanup = _delete_previous_harmonization_pipeline_runs(db, pipeline_key, run_key)
        except Exception as cleanup_exc:
            cleanup_error = str(cleanup_exc)
        completed_at = datetime.now(timezone.utc).isoformat()
        completed_update = {
            "_key": run_key,
            "status": "complete",
            "completed_at": completed_at,
            "elapsed_seconds": round(time.time() - started, 1),
            "replaced_run_keys": cleanup["deleted_run_keys"],
            "deleted_superseded_stage_keys": cleanup["deleted_superseded_stage_keys"],
            "deleted_orphan_stage_keys": cleanup["deleted_orphan_stage_keys"],
            "cleanup_error": cleanup_error,
        }
        run_collection.update(completed_update)
        db.collection(_HARMONIZATION_PIPELINE_COLLECTION).update({
            "_key": pipeline_key,
            "latest_run_key": run_key,
            "latest_stage_key": stage_docs[-1]["_key"] if stage_docs else None,
            "engine_version": _HARMONIZATION_ENGINE_VERSION,
            "updated_at": completed_at,
            "status": "complete",
        })
        return {**run_doc, **cleaning_update, **completed_update}
    except Exception as exc:
        failed_at = datetime.now(timezone.utc).isoformat()
        run_collection.update({
            "_key": run_key,
            "status": "failed",
            "completed_at": failed_at,
            "error": str(exc),
        })
        db.collection(_HARMONIZATION_PIPELINE_COLLECTION).update({
            "_key": pipeline_key,
            "updated_at": failed_at,
            "status": "failed",
        })
        raise


def _load_harmonization_pipeline_workbench() -> dict:
    _try_reconcile_interrupted_harmonization_runs()
    jobs = _list_metabolite_snapshot_jobs()
    active_pipeline_keys = _active_harmonization_pipeline_job_keys(jobs)
    deleting_pipeline_keys = _active_harmonization_delete_pipeline_keys(jobs)
    active_jobs = [
        job for job in jobs
        if job.get("status") in {"queued", "running"}
    ]
    try:
        db = get_db("metabolite_harmonization")
        _ensure_harmonization_pipeline_collections(db)
        pipelines = _list_harmonization_pipelines()
        curation_state = _load_metabolite_edge_removal_curations()
        pipeline_stage_stats = _list_harmonization_pipeline_stage_overview_stats(
            pipelines,
            curation_state,
            include_distribution=False,
            include_assertion_matrix=False,
        )
        available_stage_keys = {
            stage.get("_key")
            for pipeline_stats in pipeline_stage_stats
            for stage in pipeline_stats.get("stages") or []
            if stage.get("_key")
        }
        current_graph_fingerprint = (
            None if active_jobs else _harmonization_graph_fingerprint(db)
        )
        _annotate_harmonization_pipeline_curation_status(
            pipelines,
            curation_state,
            curation_state["assertion_fingerprint"],
            available_stage_keys,
            (current_graph_fingerprint or {}).get(
                "source_graph_content_fingerprint"
            ),
        )
        _annotate_harmonization_pipeline_run_progress(pipelines, jobs)
        return {
            "available_rules": _METABOLITE_HARMONIZATION_RULES,
            "pipelines": pipelines,
            "pipeline_tree": _build_harmonization_pipeline_tree(pipelines, pipeline_stage_stats),
            "pipeline_stage_stats": pipeline_stage_stats,
            "jobs": jobs,
            "active_jobs": active_jobs,
            "jobs_by_pipeline_key": _harmonization_jobs_by_pipeline_key(jobs),
            "active_job_count": len(active_jobs),
            "active_pipeline_keys": sorted(active_pipeline_keys),
            "deleting_pipeline_keys": sorted(deleting_pipeline_keys),
            "error": None,
        }
    except Exception as exc:
        return {
            "available_rules": _METABOLITE_HARMONIZATION_RULES,
            "pipelines": [],
            "pipeline_tree": {
                "roots": [],
                "pipeline_count": 0,
                "represented_pipeline_count": 0,
                "not_run_count": 0,
            },
            "pipeline_stage_stats": [],
            "jobs": jobs,
            "active_jobs": active_jobs,
            "jobs_by_pipeline_key": _harmonization_jobs_by_pipeline_key(jobs),
            "active_job_count": len(active_jobs),
            "active_pipeline_keys": sorted(active_pipeline_keys),
            "deleting_pipeline_keys": sorted(deleting_pipeline_keys),
            "error": str(exc),
        }


def _load_harmonization_pipeline_page(pipeline_key: str) -> dict:
    db = get_db("metabolite_harmonization")
    _ensure_harmonization_pipeline_collections(db)
    pipelines = _list_harmonization_pipelines(limit=1, pipeline_key=pipeline_key)
    if not pipelines:
        raise ValueError(f"Harmonization pipeline {pipeline_key} does not exist.")
    pipeline = pipelines[0]
    curation_state = _load_metabolite_edge_removal_curations()
    stage_stats = _list_harmonization_pipeline_stage_overview_stats(pipelines, curation_state)
    return {
        "pipeline": pipeline,
        "overview": {
            "pipeline_stage_stats": stage_stats,
        },
    }


def _get_harmonization_stage(stage_key: str) -> dict:
    db = get_db("metabolite_harmonization")
    if not db.has_collection(_HARMONIZATION_STAGE_COLLECTION):
        raise ValueError("Harmonization stage collection does not exist.")
    stage = db.collection(_HARMONIZATION_STAGE_COLLECTION).get(stage_key)
    if not stage:
        raise ValueError(f"Harmonization stage {stage_key} does not exist.")
    if stage.get("status") == "deleting":
        raise ValueError(
            f"Harmonization stage {stage_key} is being retired; reload after cleanup completes."
        )
    return stage


def _legacy_mw_finding_metadata(validation: dict, warning: dict) -> dict:
    """Give pre-fingerprint stages a deterministic, reviewable finding identity."""
    member_ids = sorted(set(
        warning.get("member_ids")
        or str(warning.get("comparison_ids") or "").split()
        or warning.get("sample_member_ids")
        or []
    ))
    anchor_candidates = _kegg_first(member_ids)
    identity = {
        "check": "mw_spread",
        "validator_version": validation.get("version"),
        "member_ids": member_ids,
    }
    evidence = {
        **identity,
        "fingerprint_version": "legacy-saved-finding-v1",
        "threshold": validation.get("threshold"),
        "mass_summaries": warning.get("mass_summaries"),
        "channel_results": warning.get("channel_results"),
        "component_matches": warning.get("component_matches") or [],
        "member_mass_examples": warning.get("member_mass_examples") or [],
    }
    return {
        "finding_id": "mw-" + hashlib.sha256(
            json.dumps(identity, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest()[:24],
        "anchor_id": anchor_candidates[0] if anchor_candidates else None,
        "member_ids": member_ids,
        "evidence_fingerprint": hashlib.sha256(
            json.dumps(evidence, sort_keys=True, separators=(",", ":")).encode("utf-8")
        ).hexdigest(),
        "evidence_fingerprint_version": "legacy-saved-finding-v1",
    }


def _compact_mw_finding_examples(warning: dict, limit: int = 3) -> List[dict]:
    examples = warning.get("member_mass_examples") or []
    if len(examples) <= limit:
        return examples

    def values(item: dict) -> List[float]:
        return [
            float(value)
            for value in [
                *(item.get("average_masses") or []),
                *(item.get("monoisotopic_masses") or []),
                *(item.get("masses") or []),
            ]
            if value is not None
        ]

    with_values = [item for item in examples if values(item)]
    selected = []
    if with_values:
        selected.extend([
            min(with_values, key=lambda item: min(values(item))),
            max(with_values, key=lambda item: max(values(item))),
        ])
    explanatory_ids = [
        identifier
        for match in warning.get("component_matches") or []
        for identifier in (match.get("component_member_id"), match.get("whole_member_id"))
        if identifier
    ]
    selected.extend(
        item for identifier in explanatory_ids
        for item in examples
        if item.get("member_id") == identifier
    )
    selected.extend(examples)
    return list({item.get("member_id"): item for item in selected}.values())[:limit]


def _mw_finding_mass_value_rows(finding: dict) -> List[dict]:
    """Group exact validator masses by channel and contributing identifier."""
    rows = []
    for channel, label in (
        ("average", "Average"),
        ("monoisotopic", "Monoisotopic"),
        ("unspecified", "Unspecified (legacy)"),
    ):
        identifiers_by_mass: dict[float, set[str]] = {}
        for example in finding.get("member_mass_examples") or []:
            identifier = str(example.get("member_id") or "").strip()
            if not identifier:
                continue
            values = (
                example.get("masses")
                if channel == "unspecified"
                else example.get(f"{channel}_masses")
            ) or []
            for raw_value in values:
                parsed = _parse_metabolite_mass(raw_value)
                if parsed is None:
                    continue
                identifiers_by_mass.setdefault(parsed, set()).add(identifier)
        rows.extend({
            "channel": channel,
            "channel_label": label,
            "mass": mass,
            "mass_label": format(mass, ".12g"),
            "identifiers": _kegg_first(sorted(identifiers)),
        } for mass, identifiers in sorted(
            identifiers_by_mass.items(), key=lambda item: item[0], reverse=True
        ))
    return rows


def _mw_finding_mass_cluster_rows(
    finding: dict,
    mass_value_rows: Optional[List[dict]] = None,
) -> List[dict]:
    """Summarize exact saved evidence with the validator's clustering semantics."""
    exact_rows = mass_value_rows or _mw_finding_mass_value_rows(finding)
    channel_results = finding.get("channel_results") or {}
    legacy_sampled_evidence = (
        finding.get("evidence_fingerprint_version") == "legacy-saved-finding-v1"
    )
    cluster_rows = []
    for channel, label in (
        ("average", "Average"),
        ("monoisotopic", "Monoisotopic"),
        ("unspecified", "Unspecified (legacy)"),
    ):
        channel_exact_rows = [row for row in exact_rows if row["channel"] == channel]
        result = channel_results.get(channel) or {}
        summary = result.get("summary") or {}
        if not channel_exact_rows and not summary:
            continue
        observations = [
            {"member_id": identifier, "value": row["mass"]}
            for row in channel_exact_rows
            for identifier in row["identifiers"]
        ]
        clusters = cluster_mass_observations(observations) if observations else []
        reported_cluster_count = result.get("mass_cluster_count")
        evidence_complete = not legacy_sampled_evidence and (
            reported_cluster_count is None
            or reported_cluster_count == len(clusters)
        )
        display_clusters = []
        if evidence_complete:
            for cluster in sorted(
                clusters,
                key=lambda items: max(item["value"] for item in items),
                reverse=True,
            ):
                values = [item["value"] for item in cluster]
                identifiers = _kegg_first(sorted({
                    item["member_id"] for item in cluster
                }))
                minimum = min(values)
                maximum = max(values)
                minimum_label = format(minimum, ".6g")
                maximum_label = format(maximum, ".6g")
                display_clusters.append({
                    "mass_label": (
                        f"≈{minimum_label} Da"
                        if minimum_label == maximum_label
                        else f"{minimum_label}–{maximum_label} Da"
                    ),
                    "identifier_count": len(identifiers),
                })
        fallback_range_label = None
        if not evidence_complete and summary.get("min") is not None:
            minimum_label = format(float(summary["min"]), ".6g")
            maximum_label = format(float(summary.get("max", summary["min"])), ".6g")
            fallback_range_label = (
                f"≈{minimum_label} Da"
                if minimum_label == maximum_label
                else f"{minimum_label}–{maximum_label} Da"
            )
        cluster_rows.append({
            "channel": channel,
            "channel_label": label,
            "cluster_count": (
                reported_cluster_count
                if reported_cluster_count is not None
                else len(clusters)
            ),
            "explained_cluster_count": result.get("explained_mass_cluster_count"),
            "clusters": display_clusters,
            "fallback_range_label": fallback_range_label,
            "evidence_complete": evidence_complete,
        })
    return cluster_rows


def _with_mw_finding_metadata(validation: dict) -> dict:
    normalized = dict(validation)
    normalized["warnings"] = []
    for warning in validation.get("warnings") or []:
        normalized_warning = {
            **warning,
            **(
                {
                    key: warning.get(key)
                    for key in (
                        "finding_id", "anchor_id", "member_ids",
                        "evidence_fingerprint", "evidence_fingerprint_version",
                    )
                }
                if all(warning.get(key) for key in (
                    "finding_id", "anchor_id", "member_ids", "evidence_fingerprint"
                ))
                else _legacy_mw_finding_metadata(validation, warning)
            ),
        }
        normalized_warning["compact_member_mass_examples"] = (
            _compact_mw_finding_examples(normalized_warning)
        )
        normalized["warnings"].append(normalized_warning)
    if isinstance(validation.get("finding_index"), list):
        normalized["finding_index"] = [dict(finding) for finding in validation["finding_index"]]
    return normalized


def _load_metabolite_mw_adjudications(storage=None) -> dict:
    use_cache = storage is None
    if use_cache:
        with _harmonization_curation_cache_lock:
            cache_generation = _harmonization_curation_cache_generation
            cached = _harmonization_curation_cache.get("mw_adjudications")
            if cached and cached["expires_at"] > time.monotonic():
                return cached["value"]
    if storage is None:
        storage = _curation_cart_storage()
    snapshot = resolve_curation_type(
        storage, METABOLITE_MW_ADJUDICATIONS, allow_missing=True
    )
    decisions = []
    # MW review can branch when a clique splits. Read immutable history newest-first
    # rather than collapsing decisions by anchor, so each descendant keeps context.
    for resolved in reversed(snapshot.operations):
        operation = resolved.operation
        decisions.append({
            **operation,
            "curation_batch_id": resolved.batch_id,
            "published_at": resolved.published_at,
            "published_by": resolved.published_by,
        })
    result = {"snapshot": snapshot.metadata(), "decisions": decisions}
    if use_cache:
        with _harmonization_curation_cache_lock:
            if cache_generation == _harmonization_curation_cache_generation:
                _harmonization_curation_cache["mw_adjudications"] = {
                    "value": result,
                    "expires_at": (
                        time.monotonic()
                        + _HARMONIZATION_CURATION_CACHE_TTL_SECONDS
                    ),
                }
    return result


def _mw_adjudication_matches_members(decision: dict, member_ids: Iterable[str]) -> bool:
    """Match a review lineage even when its original anchor left the clique."""
    members = set(member_ids)
    target = decision.get("target") or {}
    if target.get("anchor_id") in members:
        return True
    return bool(members.intersection(decision.get("observed_member_ids") or []))


def _mw_decision_matches_evidence(
    validation: dict,
    decision: dict,
    finding: dict,
) -> bool:
    observed = decision.get("observed_evidence_fingerprint")
    if observed == finding.get("evidence_fingerprint"):
        return True
    if observed == finding.get("legacy_evidence_fingerprint"):
        return True
    if not any(
        key in finding
        for key in (
            "mass_summaries", "channel_results", "component_matches",
            "member_mass_examples",
        )
    ):
        return False
    legacy = _legacy_mw_finding_metadata(validation, finding)
    return observed == legacy.get("evidence_fingerprint")


def _mw_finding_with_adjudication(
    validation: dict,
    finding: dict,
    decisions: List[dict],
) -> dict:
    members = set(finding.get("member_ids") or [])
    decision = next((
        decision for decision in decisions
        if _mw_adjudication_matches_members(decision, members)
    ), None)
    review_status = "unresolved"
    if decision and decision.get("action") == "accept_mw_discrepancy":
        review_status = (
            "accepted"
            if _mw_decision_matches_evidence(validation, decision, finding)
            else "stale"
        )
    return {
        **finding,
        "detected_severity": finding.get("severity", "error"),
        "effective_severity": (
            "warning" if review_status == "accepted" else finding.get("severity", "error")
        ),
        "review_status": review_status,
        "adjudication": decision,
    }


def _mw_validation_overview_from_doc(
    stage: dict,
    adjudication_state: Optional[dict],
) -> dict:
    """Compute card counts without formatting full finding evidence."""
    validation = (stage.get("validation") or {}).get("mw_spread")
    if not isinstance(validation, dict):
        return {"computed": False}
    finding_index = validation.get("finding_index")
    if not isinstance(finding_index, list):
        # Legacy stages did not persist a complete finding index. Their scalar
        # totals remain the only complete overview available.
        return {
            **validation,
            "accepted_count": validation.get("accepted_count", 0),
            "stale_acceptance_count": validation.get("stale_acceptance_count", 0),
        }
    decisions = (adjudication_state or {}).get("decisions") or []
    findings = [
        _mw_finding_with_adjudication(validation, finding, decisions)
        for finding in finding_index
    ]
    detected_error_count = sum(
        finding.get("detected_severity") == "error" for finding in findings
    )
    accepted_error_count = sum(
        finding.get("detected_severity") == "error"
        and finding.get("review_status") == "accepted"
        for finding in findings
    )
    return {
        "computed": validation.get("computed", False),
        "warning_count": validation.get("warning_count", len(findings)),
        "error_count": max(0, detected_error_count - accepted_error_count),
        "review_warning_count": sum(
            finding.get("detected_severity") == "warning" for finding in findings
        ),
        "accepted_count": sum(
            finding.get("review_status") == "accepted" for finding in findings
        ),
        "stale_acceptance_count": sum(
            finding.get("review_status") == "stale" for finding in findings
        ),
    }


def _mw_review_sort_key(finding: dict) -> tuple:
    review_status = finding.get("review_status")
    severity = finding.get("detected_severity") or finding.get("severity", "error")
    if review_status == "accepted":
        priority = 3
    elif review_status == "stale":
        priority = 1
    elif severity == "warning":
        priority = 2
    else:
        priority = 0
    return (
        priority,
        -(finding.get("spread_percent") or 0),
        -(finding.get("size") or 0),
        finding.get("finding_id") or "",
    )


def _with_mw_adjudications(validation: dict, adjudication_state: Optional[dict]) -> dict:
    normalized = _with_mw_finding_metadata(validation)
    decisions = (adjudication_state or {}).get("decisions") or []

    warnings = [
        _mw_finding_with_adjudication(validation, warning, decisions)
        for warning in normalized.get("warnings") or []
    ]
    warnings.sort(key=_mw_review_sort_key)
    normalized["warnings"] = warnings
    has_full_finding_index = isinstance(normalized.get("finding_index"), list)
    warnings_by_finding_id = {
        finding.get("finding_id"): finding
        for finding in warnings
        if finding.get("finding_id")
    }
    counted_findings = (
        [
            _mw_finding_with_adjudication(
                validation,
                warnings_by_finding_id.get(finding.get("finding_id"), finding),
                decisions,
            )
            for finding in normalized["finding_index"]
        ]
        if has_full_finding_index else warnings
    )
    if has_full_finding_index:
        normalized["finding_index"] = counted_findings

    # New stages carry a lightweight index for every finding. Older stages retain
    # only a display sample, so preserve their full scalar totals as a fallback.
    detected_error_count = normalized.get("error_count")
    if has_full_finding_index or detected_error_count is None:
        detected_error_count = sum(
            finding.get("detected_severity") == "error" for finding in counted_findings
        )
    detected_review_warning_count = normalized.get("review_warning_count")
    if has_full_finding_index or detected_review_warning_count is None:
        detected_review_warning_count = sum(
            finding.get("detected_severity") == "warning" for finding in counted_findings
        )
    accepted_errors = sum(
        finding.get("detected_severity") == "error"
        and finding.get("review_status") == "accepted"
        for finding in counted_findings
    )
    normalized["detected_error_count"] = detected_error_count
    normalized["error_count"] = max(0, detected_error_count - accepted_errors)
    normalized["review_warning_count"] = detected_review_warning_count
    normalized["accepted_count"] = sum(
        finding.get("review_status") == "accepted" for finding in counted_findings
    )
    normalized["stale_acceptance_count"] = sum(
        finding.get("review_status") == "stale" for finding in counted_findings
    )
    return normalized


def _build_mw_validation_review(stage_key: str, finding_id: str) -> dict:
    stage = _get_harmonization_stage(stage_key)
    try:
        adjudications = _load_metabolite_mw_adjudications()
    except Exception:
        adjudications = None
    validation = _harmonization_stage_mw_validation_from_doc(stage, adjudications)
    finding = next((
        item for item in validation.get("warnings") or []
        if item.get("finding_id") == finding_id
    ), None)
    if finding is None:
        raise ValueError(
            f"MW finding {finding_id} is not present in stage {stage_key}."
        )
    mass_value_rows = _mw_finding_mass_value_rows(finding)
    finding = {
        **finding,
        "mass_value_rows": mass_value_rows,
        "mass_cluster_rows": _mw_finding_mass_cluster_rows(
            finding, mass_value_rows
        ),
        "mass_evidence_id_count": len({
            identifier
            for row in mass_value_rows
            for identifier in row["identifiers"]
        }),
    }
    adjudication = finding.get("adjudication") or {}
    current_snapshot = mw_review_evidence_snapshot(
        finding,
        validator_version=validation.get("version"),
        threshold=validation.get("threshold"),
    )
    if finding.get("review_status") == "stale" and adjudication:
        finding["change_summary"] = compare_mw_review_evidence(
            current_member_ids=finding.get("member_ids") or [],
            current_snapshot=current_snapshot,
            current_fingerprint=finding.get("evidence_fingerprint") or "",
            previous_member_ids=adjudication.get("observed_member_ids") or [],
            previous_snapshot=adjudication.get("observed_evidence_snapshot"),
            previous_fingerprint=(
                adjudication.get("observed_evidence_fingerprint") or ""
            ),
        )
        prior_supporting_ids = adjudication.get("supporting_ids") or []
        current_members = set(finding.get("member_ids") or [])
        finding["reusable_supporting_ids"] = [
            identifier for identifier in prior_supporting_ids
            if identifier in current_members
        ]
        finding["removed_supporting_ids"] = [
            identifier for identifier in prior_supporting_ids
            if identifier not in current_members
        ]
    return {
        "stage_key": stage_key,
        "stage_name": stage.get("name"),
        "validator_version": validation.get("version"),
        "threshold_percent": validation.get("threshold_percent"),
        "finding": finding,
        "reason_options": [
            {"value": value, "label": label}
            for value, label in _MW_ADJUDICATION_REASONS
        ],
    }


def _find_mw_review_stage(
    finding_id: str,
    anchor_id: str,
    member_ids: List[str],
) -> Optional[dict]:
    """Find the newest exact finding, or its closest current descendant."""
    db = get_db("metabolite_harmonization")
    if not db.has_collection(_HARMONIZATION_STAGE_COLLECTION):
        return None
    candidates = list(db.aql.execute(
        f"""
        FOR stage IN {_HARMONIZATION_STAGE_COLLECTION}
          FILTER stage.status == "complete"
          LET mw = stage.validation.mw_spread
          // The detail page needs the full warning payload. The lightweight
          // finding_index also contains findings beyond the saved display
          // limit, but those cannot be opened or curated from this stage.
          LET findings = IS_ARRAY(mw.warnings) ? mw.warnings : []
          FOR finding IN findings
            LET members = IS_ARRAY(finding.member_ids) ? finding.member_ids : []
            LET exact = finding.finding_id == @finding_id
            LET overlap = LENGTH(INTERSECTION(members, @member_ids))
            FILTER exact OR overlap >= 2 OR @anchor_id IN members
            SORT exact DESC, overlap DESC, stage.created_at DESC
            LIMIT 1
            RETURN {{
              stage_key: stage._key,
              finding_id: finding.finding_id,
              exact: exact,
              overlap: overlap
            }}
        """,
        bind_vars={
            "finding_id": finding_id,
            "anchor_id": anchor_id,
            "member_ids": member_ids,
        },
        max_runtime=120,
    ))
    return candidates[0] if candidates else None


def _harmonization_stage_mw_validation_from_doc(
    stage: dict,
    adjudication_state: Optional[dict] = None,
    threshold: float = _METABOLITE_MW_SPREAD_WARNING_THRESHOLD,
    limit: int = _METABOLITE_MW_SPREAD_WARNING_LIMIT,
) -> dict:
    existing = (stage.get("validation") or {}).get("mw_spread")
    if isinstance(existing, dict):
        return _with_mw_adjudications(
            _prioritize_kegg_validation_samples(existing), adjudication_state
        )
    return {
        "computed": False,
        "version": None,
        "threshold": threshold,
        "threshold_percent": threshold * 100,
        "warning_count": 0,
        "issue_count": 0,
        "error_count": 0,
        "review_warning_count": 0,
        "display_limit": limit,
        "warnings": [],
    }


def _harmonization_stage_generic_structure_validation_from_doc(stage: dict) -> dict:
    existing = (stage.get("validation") or {}).get("generic_structure_consistency")
    if isinstance(existing, dict):
        return _prioritize_kegg_validation_samples(existing)
    return {
        "computed": False,
        "rule_enabled": "ignore_generic_structure_mismatch" in (stage.get("rule_ids") or []),
        "inconsistent_clique_count": 0,
        "generic_unknown_clique_count": 0,
        "actionable_raw_edge_count": 0,
        "warning_count": 0,
        "display_limit": _METABOLITE_GENERIC_STRUCTURE_WARNING_LIMIT,
        "warnings": [],
    }


def _load_harmonization_stage_stats(stage_key: str) -> dict:
    db = get_db("metabolite_harmonization")
    stage = _get_harmonization_stage(stage_key)
    counts = {
        "harmonized_metabolite_count": 0,
        "member_edge_count": 0,
        "evidence_edge_count": 0,
        "active_identifier_count_from_chunks": 0,
    }
    if db.has_collection(_HARMONIZED_METABOLITE_COLLECTION):
        counts["harmonized_metabolite_count"] = list(db.aql.execute(
            f"""
            RETURN LENGTH(
              FOR d IN {_HARMONIZED_METABOLITE_COLLECTION}
                FILTER d.stage_key == @stage_key
                RETURN 1
            )
            """,
            bind_vars={"stage_key": stage_key},
            max_runtime=120,
        ))[0]
    if db.has_collection(_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION):
        counts["member_edge_count"] = list(db.aql.execute(
            f"""
            RETURN LENGTH(
              FOR d IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
                FILTER d.stage_key == @stage_key
                RETURN 1
            )
            """,
            bind_vars={"stage_key": stage_key},
            max_runtime=120,
        ))[0]
    if db.has_collection(_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION):
        counts["evidence_edge_count"] = list(db.aql.execute(
            f"""
            RETURN LENGTH(
              FOR d IN {_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION}
                FILTER d.stage_key == @stage_key
                RETURN 1
            )
            """,
            bind_vars={"stage_key": stage_key},
            max_runtime=120,
        ))[0]
    if db.has_collection(_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION):
        counts["active_identifier_count_from_chunks"] = list(db.aql.execute(
            f"""
            RETURN SUM(
              FOR d IN {_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION}
                FILTER d.stage_key == @stage_key
                RETURN LENGTH(d.identifier_ids || [])
            )
            """,
            bind_vars={"stage_key": stage_key},
            max_runtime=120,
        ))[0] or 0
    largest_groups = []
    if db.has_collection(_HARMONIZED_METABOLITE_COLLECTION):
        largest_groups = list(db.aql.execute(
            f"""
            FOR d IN {_HARMONIZED_METABOLITE_COLLECTION}
              FILTER d.stage_key == @stage_key
              SORT d.size DESC, d._key ASC
              LIMIT 20
              RETURN KEEP(d, "_key", "id", "name", "size", "rank_by_size", "representative_id", "sample_member_ids")
            """,
            bind_vars={"stage_key": stage_key},
            max_runtime=120,
        ))
    source_counts = []
    if db.has_collection(_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION):
        source_counts = list(db.aql.execute(
            f"""
            FOR e IN {_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION}
              FILTER e.stage_key == @stage_key
              FOR source IN e.sources || []
                COLLECT source_name = source WITH COUNT INTO edge_count
                SORT edge_count DESC, source_name ASC
                LIMIT 20
                RETURN {{source: source_name, edge_count: edge_count}}
            """,
            bind_vars={"stage_key": stage_key},
            max_runtime=120,
        ))
    curation_state = _load_metabolite_edge_removal_curations()
    expected_clique_assertions = _load_expected_clique_assertion_results(
        db,
        curation_state.get("assertions") or [],
        [stage_key],
    )["by_stage"].get(stage_key, {
        "stage_key": stage_key,
        "assertion_count": 0,
        "pass_count": 0,
        "split_count": 0,
        "missing_count": 0,
        "failure_count": 0,
        "results": [],
    })
    expected_clique_assertions["fingerprint"] = curation_state.get("assertion_fingerprint")
    try:
        mw_adjudication_state = _load_metabolite_mw_adjudications()
    except Exception:
        mw_adjudication_state = None
    mw_validation = _harmonization_stage_mw_validation_from_doc(
        stage, mw_adjudication_state
    )
    denylist_validation = _harmonization_stage_denylist_validation_from_doc(stage)
    generic_structure_validation = _harmonization_stage_generic_structure_validation_from_doc(stage)
    generic_structure_validation = {
        **generic_structure_validation,
        "warning_count": generic_structure_validation.get("inconsistent_clique_count", 0),
        "warnings": [
            warning
            for warning in generic_structure_validation.get("warnings") or []
            if warning.get("status") == "inconsistent"
        ],
    }
    review_ranks = {
        int(warning["rank_by_size"])
        for warning in generic_structure_validation.get("warnings") or []
        if warning.get("rank_by_size") is not None
    }
    review_id_by_rank = _load_stage_validation_review_ids(db, stage_key, review_ranks)
    generic_structure_validation = _with_validation_review_ids(
        generic_structure_validation,
        review_id_by_rank,
    )
    return {
        "stage": stage,
        "summary": stage.get("summary") or {},
        "counts": counts,
        "largest_groups": largest_groups,
        "source_counts": source_counts,
        "mw_validation": mw_validation,
        "denylist_validation": denylist_validation,
        "generic_structure_validation": generic_structure_validation,
        "expected_clique_assertions": expected_clique_assertions,
    }


def _build_harmonization_stage_cart_flags(
    operations: List[dict],
    member_rank_by_id: Dict[str, int],
    active_edge_pairs: set[tuple[str, str]],
    warning_ranks: set[int],
    denylist_warning_pairs: Optional[set[tuple[str, str]]] = None,
) -> List[dict]:
    denylist_warning_pairs = denylist_warning_pairs or set()
    decisions_by_rank: Dict[int, set[tuple[str, str, str]]] = {}
    updates_by_rank: Dict[int, Dict[str, List[dict]]] = {}
    validation_decisions_by_rank: Dict[int, List[dict]] = {}
    for operation in operations or []:
        if operation.get("action") in {"accept_mw_discrepancy", "reopen_mw_discrepancy"}:
            target = operation.get("target") or {}
            anchor_id = target.get("anchor_id")
            rank = member_rank_by_id.get(anchor_id)
            if rank in warning_ranks:
                validation_decisions_by_rank.setdefault(rank, []).append({
                    "action": operation["action"],
                    "finding_id": target.get("finding_id"),
                    "note": operation.get("note"),
                })
            continue
        if operation.get("action") == "set_properties" and operation.get("decisions"):
            target = operation.get("target") or {}
            target_id = _normalize_ramp_denylist_identifier(target.get("id"))
            rank = member_rank_by_id.get(target_id)
            if target_id and rank in warning_ranks:
                changes = [{
                    "property": _display_record_property_path(decision.get("path") or []),
                    "mode": decision.get("mode"),
                    "value": decision.get("value"),
                } for decision in operation.get("decisions") or []]
                updates_by_rank.setdefault(rank, {})[target_id] = changes
            continue
        if (
            operation.get("action") not in {"remove_edge", "retain_edge"}
            or operation.get("edge_type") != "MetaboliteIdentifierMappingEdge"
        ):
            continue
        left = _normalize_ramp_denylist_identifier(operation.get("start_id"))
        right = _normalize_ramp_denylist_identifier(operation.get("end_id"))
        if not left or not right or left == right:
            continue
        pair = tuple(sorted((left, right)))
        left_rank = member_rank_by_id.get(left)
        right_rank = member_rank_by_id.get(right)
        if (
            left_rank is None
            or left_rank != right_rank
            or left_rank not in warning_ranks
        ):
            continue
        action = operation["action"]
        if action == "remove_edge" and pair not in active_edge_pairs:
            continue
        if action == "retain_edge" and pair not in denylist_warning_pairs:
            continue
        decisions_by_rank.setdefault(left_rank, set()).add((action, *pair))
    return [
        {
            "rank_by_size": rank,
            "queued_edge_count": len(decisions),
            "updated_identifier_count": len(updates_by_rank.get(rank, {})),
            "edges": [
                {"action": action, "start_id": left, "end_id": right}
                for action, left, right in sorted(decisions)
            ],
            "identifier_updates": [
                {"target_id": target_id, "changes": changes}
                for target_id, changes in sorted(updates_by_rank.get(rank, {}).items())
            ],
            **({
                "validation_decision_count": len(validation_decisions_by_rank[rank]),
                "validation_decisions": validation_decisions_by_rank[rank],
            } if rank in validation_decisions_by_rank else {}),
        }
        for rank in sorted(
            set(decisions_by_rank) | set(updates_by_rank) | set(validation_decisions_by_rank)
        )
        for decisions in [decisions_by_rank.get(rank, set())]
    ]


def _harmonization_stage_cart_warning_ranks(stage: dict) -> set[int]:
    warnings = [
        *(_harmonization_stage_mw_validation_from_doc(stage).get("warnings") or []),
        *(_harmonization_stage_denylist_validation_from_doc(stage).get("warnings") or []),
        *(_harmonization_stage_generic_structure_validation_from_doc(stage).get("warnings") or []),
    ]
    return {
        int(warning["rank_by_size"])
        for warning in warnings
        if warning.get("rank_by_size") is not None
    }


def _load_harmonization_stage_cart_flags(stage_key: str, operations: List[dict]) -> dict:
    stage = _get_harmonization_stage(stage_key)
    warning_ranks = _harmonization_stage_cart_warning_ranks(stage)
    denylist_warning_pairs = {
        tuple(sorted((pair["left_id"], pair["right_id"])))
        for warning in (_harmonization_stage_denylist_validation_from_doc(stage).get("warnings") or [])
        for pair in (warning.get("pairs") or [])
        if pair.get("left_id") and pair.get("right_id")
    }
    edge_endpoint_ids = sorted({
        normalized
        for operation in operations or []
        if (
            operation.get("action") in {"remove_edge", "retain_edge"}
            and operation.get("edge_type") == "MetaboliteIdentifierMappingEdge"
        )
        for normalized in (
            _normalize_ramp_denylist_identifier(operation.get("start_id")),
            _normalize_ramp_denylist_identifier(operation.get("end_id")),
        )
        if normalized
    })
    property_target_ids = sorted({
        normalized
        for operation in operations or []
        if operation.get("action") == "set_properties"
        for normalized in [
            _normalize_ramp_denylist_identifier((operation.get("target") or {}).get("id"))
        ]
        if normalized
    })
    validation_anchor_ids = sorted({
        normalized
        for operation in operations or []
        if operation.get("action") in {"accept_mw_discrepancy", "reopen_mw_discrepancy"}
        for normalized in [
            _normalize_ramp_denylist_identifier((operation.get("target") or {}).get("anchor_id"))
        ]
        if normalized
    })
    subject_ids = sorted(
        set(edge_endpoint_ids) | set(property_target_ids) | set(validation_anchor_ids)
    )
    if not subject_ids or not warning_ranks:
        return {"stage_key": stage_key, "flags": [], "flagged_warning_count": 0}

    db = get_db("metabolite_harmonization")
    member_rows = list(db.aql.execute(
        f"""
        FOR e IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
          FILTER e.stage_key == @stage_key
          FILTER e.member_id IN @subject_ids
          LET clique = DOCUMENT(e._from)
          FILTER clique != null
          RETURN {{member_id: e.member_id, rank_by_size: clique.rank_by_size}}
        """,
        bind_vars={"stage_key": stage_key, "subject_ids": subject_ids},
        max_runtime=120,
    ))
    member_rank_by_id = {
        row["member_id"]: int(row["rank_by_size"])
        for row in member_rows
        if row.get("member_id") and row.get("rank_by_size") is not None
    }
    edge_rows = list(db.aql.execute(
        f"""
        FOR e IN {_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION}
          FILTER e.stage_key == @stage_key
          FILTER e.start_id IN @endpoint_ids AND e.end_id IN @endpoint_ids
          RETURN {{start_id: e.start_id, end_id: e.end_id}}
        """,
        bind_vars={"stage_key": stage_key, "endpoint_ids": edge_endpoint_ids},
        max_runtime=120,
    ))
    active_edge_pairs = {
        tuple(sorted((row["start_id"], row["end_id"])))
        for row in edge_rows
        if row.get("start_id") and row.get("end_id")
    }
    flags = _build_harmonization_stage_cart_flags(
        operations,
        member_rank_by_id,
        active_edge_pairs,
        warning_ranks,
        denylist_warning_pairs,
    )
    return {
        "stage_key": stage_key,
        "flags": flags,
        "flagged_warning_count": len(flags),
    }


def _load_stage_harmonized_member_sets(stage_key: str) -> List[dict]:
    db = get_db("metabolite_harmonization")
    if (
        not db.has_collection(_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION)
        or not db.has_collection(_HARMONIZED_METABOLITE_COLLECTION)
    ):
        return []
    return list(db.aql.execute(
        f"""
        FOR c IN {_HARMONIZED_METABOLITE_COLLECTION}
          FILTER c.stage_key == @stage_key
          LET member_ids = (
            FOR e IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
              FILTER e._from == c._id
              SORT e.member_id
              RETURN e.member_id
          )
          FILTER LENGTH(member_ids) > 1
          SORT c.rank_by_size ASC
          RETURN {{
            clique_key: c._key,
            clique_id: c.id,
            size: LENGTH(member_ids),
            rank_by_size: c.rank_by_size,
            representative_id: c.representative_id,
            generic_structure_status: c.generic_structure_status,
            generic_member_count: c.generic_member_count || 0,
            specific_member_count: c.specific_member_count || 0,
            unknown_member_count: c.unknown_member_count || 0,
            generic_member_samples: c.generic_member_samples || [],
            specific_member_samples: c.specific_member_samples || [],
            unknown_member_samples: c.unknown_member_samples || [],
            member_ids: member_ids,
            signature: CONCAT_SEPARATOR("\\n", member_ids)
          }}
        """,
        bind_vars={"stage_key": stage_key},
        batch_size=1000,
        max_runtime=600,
    ))


def _load_metabolite_identifier_display_map(member_ids: List[str]) -> Dict[str, dict]:
    if not member_ids:
        return {}
    db = get_db("metabolite_harmonization")
    rows = list(db.aql.execute(
        """
        FOR node IN MetaboliteIdentifier
          FILTER node.id IN @member_ids
          LET chemical_entity = STARTS_WITH(node.id, "CHEBI:") ? DOCUMENT("ChemicalEntity", node.id) : null
          LET node_names = UNIQUE(
            FOR name IN node.names || []
              FILTER name.value != null AND name.value != ""
              RETURN name.value
          )
          LET chemical_names = UNIQUE(APPEND(
            chemical_entity != null && chemical_entity.name != null && chemical_entity.name != "" ? [chemical_entity.name] : [],
            (
              FOR synonym IN chemical_entity != null ? (chemical_entity.synonyms || []) : []
                FILTER synonym.value != null AND synonym.value != ""
                RETURN synonym.value
            )
          ))
          RETURN {
            id: node.id,
            label: LENGTH(node_names) > 0 ? node_names[0] : (
              LENGTH(chemical_names) > 0 ? chemical_names[0] : node.id
            ),
            names: SLICE(UNIQUE(APPEND(node_names, chemical_names)), 0, 6),
            prefix: node.prefix
          }
        """,
        bind_vars={"member_ids": member_ids},
        max_runtime=120,
    ))
    return {row["id"]: row for row in rows}


def _sample_member_ids(member_ids: set, limit: int = 20) -> List[str]:
    return sorted(member_ids)[:limit]


def _kegg_first(items: Iterable, identifier_getter=lambda item: item) -> list:
    """Move KEGG identifiers first while preserving all other relative order."""
    return sorted(
        items,
        key=lambda item: not str(identifier_getter(item) or "").upper().startswith("KEGG"),
    )


def _prioritize_kegg_validation_samples(validation: dict) -> dict:
    """Order the identifier samples shown by stage validation tables."""
    normalized = dict(validation)
    normalized_warnings = []
    for warning in validation.get("warnings") or []:
        normalized_warning = dict(warning)
        for field in (
            "sample_member_ids",
            "generic_member_samples",
            "specific_member_samples",
            "unknown_member_samples",
        ):
            if isinstance(warning.get(field), list):
                normalized_warning[field] = _kegg_first(warning[field])
        if isinstance(warning.get("member_mass_examples"), list):
            normalized_warning["member_mass_examples"] = _kegg_first(
                warning["member_mass_examples"],
                lambda item: item.get("member_id"),
            )
        if isinstance(warning.get("comparison_ids"), str):
            normalized_warning["comparison_ids"] = " ".join(
                _kegg_first(warning["comparison_ids"].split())
            )
        normalized_warnings.append(normalized_warning)
    normalized["warnings"] = normalized_warnings
    return normalized


def _with_validation_review_ids(
    validation: dict,
    review_id_by_rank: Dict[int, str],
) -> dict:
    """Add a display-only review ID without changing clique representatives."""
    normalized = dict(validation)
    normalized["warnings"] = [
        {
            **warning,
            "review_id": review_id_by_rank.get(
                int(warning["rank_by_size"]),
                warning.get("representative_id"),
            ) if warning.get("rank_by_size") is not None else warning.get("representative_id"),
        }
        for warning in validation.get("warnings") or []
    ]
    return normalized


def _load_stage_validation_review_ids(
    db,
    stage_key: str,
    ranks: set[int],
) -> Dict[int, str]:
    if (
        not ranks
        or not db.has_collection(_HARMONIZED_METABOLITE_COLLECTION)
        or not db.has_collection(_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION)
    ):
        return {}
    rows = db.aql.execute(
        f"""
        FOR clique IN {_HARMONIZED_METABOLITE_COLLECTION}
          FILTER clique.stage_key == @stage_key
          FILTER clique.rank_by_size IN @ranks
          LET kegg_id = FIRST(
            FOR member IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
              FILTER member._from == clique._id
              FILTER STARTS_WITH(UPPER(member.member_id), "KEGG")
              SORT member.member_id
              LIMIT 1
              RETURN member.member_id
          )
          RETURN {{
            rank_by_size: clique.rank_by_size,
            review_id: kegg_id || clique.representative_id
          }}
        """,
        bind_vars={"stage_key": stage_key, "ranks": sorted(ranks)},
        max_runtime=120,
    )
    return {
        int(row["rank_by_size"]): row["review_id"]
        for row in rows
        if row.get("rank_by_size") is not None and row.get("review_id")
    }


def _format_compare_member_samples(member_ids: set, display_by_id: Dict[str, dict], limit: int = 20) -> List[dict]:
    samples = []
    for member_id in _sample_member_ids(member_ids, limit):
        display = display_by_id.get(member_id) or {}
        samples.append({
            "id": member_id,
            "label": display.get("label") or member_id,
            "names": display.get("names") or [],
            "prefix": display.get("prefix") or _identifier_prefix(member_id),
        })
    return samples


def _snapshot_compare_involved_preview(
    member_ids: set,
    display_by_id: Dict[str, dict],
    limit: int = 24,
) -> dict:
    sampled_ids = _sample_member_ids(member_ids, limit)
    names = []
    seen_names = set()
    for member_id in sampled_ids:
        label = (display_by_id.get(member_id) or {}).get("label")
        if not label or label == member_id:
            continue
        normalized_label = label.casefold()
        if normalized_label in seen_names:
            continue
        seen_names.add(normalized_label)
        names.append(label)
    return {
        "count": len(member_ids),
        "ids": sampled_ids,
        "names": names,
        "remaining_id_count": max(len(member_ids) - len(sampled_ids), 0),
    }


def _build_snapshot_comparison_review_fields(
    component: dict,
    left_key: str,
    right_key: str,
) -> None:
    review_ids = []
    for clique in component.get("right_cliques") or []:
        if clique.get("generic_structure_status") not in {"inconsistent", "generic_unknown"}:
            continue
        for member_id in [
            *(clique.get("generic_member_samples") or []),
            *(clique.get("specific_member_samples") or []),
            *(clique.get("unknown_member_samples") or []),
        ]:
            if member_id not in review_ids:
                review_ids.append(member_id)
            if len(review_ids) >= 12:
                break
        if len(review_ids) >= 12:
            break
    for assertion in component.get("newly_failing_assertions") or []:
        for member_id in assertion.get("member_ids") or []:
            if member_id not in review_ids:
                review_ids.append(member_id)
            if len(review_ids) >= 12:
                break
        if len(review_ids) >= 12:
            break
    for bucket in (
        component.get("_added_members", set()),
        component.get("_removed_members", set()),
        component.get("_retained_members", set()),
    ):
        for member_id in _sample_member_ids(bucket, 12):
            if member_id not in review_ids:
                review_ids.append(member_id)
            if len(review_ids) >= 12:
                break
        if len(review_ids) >= 12:
            break
    component["review_ids"] = review_ids
    query_params = {
        "id": " ".join(review_ids),
        "stages": ",".join([left_key, right_key]),
    }
    component["review_query_string"] = urlencode(query_params)


def _newly_failing_expected_clique_assertions(
    assertion_results: dict,
    left_key: str,
    right_key: str,
) -> List[dict]:
    left_results = {
        result["assertion_id"]: result
        for result in assertion_results.get("by_stage", {}).get(left_key, {}).get("results", [])
    }
    right_results = {
        result["assertion_id"]: result
        for result in assertion_results.get("by_stage", {}).get(right_key, {}).get("results", [])
    }
    newly_failing = []
    for assertion_id, left_result in left_results.items():
        right_result = right_results.get(assertion_id)
        if left_result.get("status") != "pass" or not right_result or right_result.get("status") == "pass":
            continue
        newly_failing.append({
            **right_result,
            "left_status": left_result.get("status"),
            "right_status": right_result.get("status"),
            "inspection_query": urlencode({
                "id": " ".join(right_result.get("member_ids") or []),
                "stages": ",".join([left_key, right_key]),
            }),
        })
    return sorted(
        newly_failing,
        key=lambda assertion: (
            (assertion.get("name") or "").casefold(),
            assertion.get("assertion_id") or "",
        ),
    )


def _snapshot_compare_clique_summary(clique: dict, sample_members: Optional[set] = None) -> dict:
    return {
        "clique_key": clique.get("clique_key"),
        "clique_id": clique.get("clique_id"),
        "size": clique.get("size"),
        "rank_by_size": clique.get("rank_by_size"),
        "representative_id": clique.get("representative_id"),
        "generic_structure_status": clique.get("generic_structure_status"),
        "generic_member_count": clique.get("generic_member_count", 0),
        "specific_member_count": clique.get("specific_member_count", 0),
        "unknown_member_count": clique.get("unknown_member_count", 0),
        "generic_member_samples": clique.get("generic_member_samples") or [],
        "specific_member_samples": clique.get("specific_member_samples") or [],
        "unknown_member_samples": clique.get("unknown_member_samples") or [],
        "member_ids": _sample_member_ids(sample_members if sample_members is not None else clique.get("member_set", set()), 12),
    }


def _snapshot_compare_nontrivial_transitions(
    source_items: List[dict],
    target_items: List[dict],
    min_size: int,
    target_label: str,
) -> List[dict]:
    transitions = []
    for source_item in source_items:
        fragments = []
        source_members = source_item["member_set"]
        for target_item in target_items:
            overlap_members = source_members & target_item["member_set"]
            if len(overlap_members) < min_size:
                continue
            fragments.append({
                "target_clique": _snapshot_compare_clique_summary(target_item, overlap_members),
                "overlap_count": len(overlap_members),
            })
        if len(fragments) < 2:
            continue
        fragments.sort(key=lambda item: (-item["overlap_count"], item["target_clique"]["clique_key"]))
        transitions.append({
            "source_clique": _snapshot_compare_clique_summary(source_item),
            "target_label": target_label,
            "fragment_count": len(fragments),
            "member_count": sum(fragment["overlap_count"] for fragment in fragments),
            "fragments": fragments,
        })
    transitions.sort(key=lambda item: (-item["fragment_count"], -item["member_count"], item["source_clique"]["clique_key"]))
    return transitions


def _snapshot_comparison_pipeline_context(
    left_snapshot: dict,
    right_snapshot: dict,
    pipelines: List[dict],
) -> dict:
    left_key = left_snapshot.get("_key")
    right_key = right_snapshot.get("_key")
    run_matches = []
    for pipeline in pipelines:
        for run in pipeline.get("runs") or []:
            stages = sorted(run.get("stages") or [], key=_stage_sort_key)
            stage_keys = [stage.get("_key") for stage in stages]
            if left_key in stage_keys or right_key in stage_keys:
                run_matches.append({
                    "pipeline": pipeline,
                    "run": run,
                    "stages": stages,
                    "stage_keys": stage_keys,
                })

    shared_run = next(
        (match for match in run_matches if left_key in match["stage_keys"] and right_key in match["stage_keys"]),
        None,
    )
    selected_matches = []
    if shared_run:
        selected_matches = [shared_run]
    else:
        left_match = next((match for match in run_matches if left_key in match["stage_keys"]), None)
        right_match = next((match for match in run_matches if right_key in match["stage_keys"]), None)
        selected_matches = [match for match in (left_match, right_match) if match]

    if not selected_matches:
        selected_matches = [
            {
                "pipeline": {"_key": "unassigned", "name": "Saved stages"},
                "run": {},
                "stages": [left_snapshot, right_snapshot],
                "stage_keys": [left_key, right_key],
            }
        ]

    lane_stage_key_sets = [set(match["stage_keys"]) for match in selected_matches]
    shared_stage_keys = (
        set.intersection(*lane_stage_key_sets)
        if len(lane_stage_key_sets) > 1
        else set(selected_matches[0]["stage_keys"])
    )
    common_prefix_keys = []
    if len(selected_matches) > 1:
        for stage_keys in zip(*(match["stage_keys"] for match in selected_matches)):
            if len(set(stage_keys)) != 1:
                break
            common_prefix_keys.append(stage_keys[0])

    lanes = []
    for match in selected_matches:
        lane_stages = []
        for fallback_index, stage in enumerate(match["stages"]):
            stage_key = stage.get("_key")
            stage_index = stage.get("stage_index")
            if stage_index is None:
                stage_index = fallback_index
            lane_stages.append({
                "key": stage_key,
                "index": stage_index,
                "number": stage_index,
                "label": stage.get("display_label") or _stage_display_label(
                    stage,
                    match["pipeline"].get("rule_parameters") or {},
                ),
                "is_left": stage_key == left_key,
                "is_right": stage_key == right_key,
                "is_shared": stage_key in shared_stage_keys,
            })
        lanes.append({
            "pipeline_key": match["pipeline"].get("_key"),
            "pipeline_name": match["pipeline"].get("name") or match["pipeline"].get("_key"),
            "run_key": match["run"].get("_key"),
            "stages": lane_stages,
        })

    stage_by_key = {
        stage["key"]: stage
        for lane in lanes
        for stage in lane["stages"]
    }
    max_stage_index = max(
        (stage["index"] for lane in lanes for stage in lane["stages"]),
        default=0,
    )
    trunk_stage_keys = (
        set(common_prefix_keys)
        if len(lanes) > 1
        else {stage["key"] for stage in lanes[0]["stages"]}
    )
    trunk_stages = [
        stage
        for stage in lanes[0]["stages"]
        if stage["key"] in trunk_stage_keys
    ]
    branches = [
        {
            **{key: value for key, value in lane.items() if key != "stages"},
            "stages": [stage for stage in lane["stages"] if stage["key"] not in trunk_stage_keys],
        }
        for lane in lanes
    ] if len(lanes) > 1 else []
    branch_stage_indexes = [
        stage["index"]
        for branch in branches
        for stage in branch["stages"]
    ]
    return {
        "lanes": lanes,
        "trunk_stages": trunk_stages,
        "branches": branches,
        "same_pipeline": len(lanes) == 1,
        "shared_stage_count": len(shared_stage_keys),
        "common_prefix_count": len(common_prefix_keys),
        "divergence_after_label": (
            stage_by_key.get(common_prefix_keys[-1], {}).get("label")
            if common_prefix_keys
            else None
        ),
        "branch_start_column": min(branch_stage_indexes, default=max_stage_index + 1),
        "column_count": max_stage_index + 1,
    }


def _load_metabolite_snapshot_comparison(
    left_key: str,
    right_key: str,
    limit: int = 100,
) -> dict:
    if not left_key or not right_key:
        raise ValueError("Choose two stages to compare.")
    if left_key == right_key:
        raise ValueError("Choose two different stages.")

    left_snapshot = _get_harmonization_stage(left_key)
    right_snapshot = _get_harmonization_stage(right_key)
    left_cliques = _load_stage_harmonized_member_sets(left_key)
    right_cliques = _load_stage_harmonized_member_sets(right_key)
    left_generic_validation = _harmonization_stage_generic_structure_validation_from_doc(left_snapshot)
    right_generic_validation = _harmonization_stage_generic_structure_validation_from_doc(right_snapshot)
    left_generic_failures = {
        clique["signature"] for clique in left_cliques
        if clique.get("generic_structure_status") == "inconsistent"
    }
    right_generic_failures = {
        clique["signature"] for clique in right_cliques
        if clique.get("generic_structure_status") == "inconsistent"
    }
    left_generic_gaps = {
        clique["signature"] for clique in left_cliques
        if clique.get("generic_structure_status") == "generic_unknown"
    }
    right_generic_gaps = {
        clique["signature"] for clique in right_cliques
        if clique.get("generic_structure_status") == "generic_unknown"
    }
    generic_validation_comparison = {
        "computed": bool(left_generic_validation.get("computed") and right_generic_validation.get("computed")),
        "left_computed": bool(left_generic_validation.get("computed")),
        "right_computed": bool(right_generic_validation.get("computed")),
        "left": left_generic_validation,
        "right": right_generic_validation,
        "new_failure_count": len(right_generic_failures - left_generic_failures),
        "persistent_failure_count": len(right_generic_failures & left_generic_failures),
        "resolved_failure_count": len(left_generic_failures - right_generic_failures),
        "new_gap_count": len(right_generic_gaps - left_generic_gaps),
        "persistent_gap_count": len(right_generic_gaps & left_generic_gaps),
        "resolved_gap_count": len(left_generic_gaps - right_generic_gaps),
        # These saved warning samples also cover unchanged cliques, which are not
        # represented in the changed-component cards below.
        "destination_warnings": right_generic_validation.get("warnings") or [],
    }
    curations = _load_metabolite_edge_removal_curations()
    assertion_results = _load_expected_clique_assertion_results(
        get_db("metabolite_harmonization"),
        curations.get("assertions") or [],
        [left_key, right_key],
    )
    newly_failing_assertions = _newly_failing_expected_clique_assertions(
        assertion_results,
        left_key,
        right_key,
    )
    right_signatures = {clique["signature"] for clique in right_cliques}
    left_signatures = {clique["signature"] for clique in left_cliques}
    changed_left = [
        {**clique, "side": "left", "member_set": set(clique["member_ids"])}
        for clique in left_cliques
        if clique["signature"] not in right_signatures
    ]
    changed_right = [
        {**clique, "side": "right", "member_set": set(clique["member_ids"])}
        for clique in right_cliques
        if clique["signature"] not in left_signatures
    ]
    node_parent = {}
    node_payload = {}
    members_by_change_node = {}

    def node_id(side: str, clique_key: str) -> str:
        return f"{side}:{clique_key}"

    def add_node(identifier: str, payload: dict) -> None:
        node_parent.setdefault(identifier, identifier)
        node_payload[identifier] = payload

    def find(identifier: str) -> str:
        while node_parent[identifier] != identifier:
            node_parent[identifier] = node_parent[node_parent[identifier]]
            identifier = node_parent[identifier]
        return identifier

    def union(left: str, right: str) -> None:
        left_root = find(left)
        right_root = find(right)
        if left_root != right_root:
            node_parent[right_root] = left_root

    member_to_change_nodes: Dict[str, List[str]] = {}
    for clique in [*changed_left, *changed_right]:
        identifier = node_id(clique["side"], clique["clique_key"])
        add_node(identifier, clique)
        members_by_change_node[identifier] = clique["member_set"]
        for member_id in clique["member_set"]:
            member_to_change_nodes.setdefault(member_id, []).append(identifier)

    for change_nodes in member_to_change_nodes.values():
        if len(change_nodes) < 2:
            continue
        first = change_nodes[0]
        for other in change_nodes[1:]:
            union(first, other)

    nodes_by_component: Dict[str, List[str]] = {}
    for identifier in node_parent:
        nodes_by_component.setdefault(find(identifier), []).append(identifier)

    components = []
    nontrivial_split_candidates = []
    nontrivial_merge_candidates = []
    recombination_candidates = []
    nontrivial_split_min_size = 2
    for component_nodes in nodes_by_component.values():
        left_items = [node_payload[node] for node in component_nodes if node_payload[node]["side"] == "left"]
        right_items = [node_payload[node] for node in component_nodes if node_payload[node]["side"] == "right"]
        left_members = set().union(*(item["member_set"] for item in left_items)) if left_items else set()
        right_members = set().union(*(item["member_set"] for item in right_items)) if right_items else set()
        retained_members = left_members & right_members
        added_members = right_members - left_members
        removed_members = left_members - right_members
        if not left_items and not right_items:
            continue
        change_score = len(added_members) + len(removed_members) + max(len(left_items), len(right_items))
        change_type = "reconfigured"
        if not left_items:
            change_type = "new clique"
        elif not right_items:
            change_type = "removed clique"
        elif len(left_items) == 1 and len(right_items) > 1:
            change_type = "split"
        elif len(left_items) > 1 and len(right_items) == 1:
            change_type = "merged"
        elif len(added_members) or len(removed_members):
            change_type = "membership changed"
        split_transitions = _snapshot_compare_nontrivial_transitions(
            left_items,
            right_items,
            nontrivial_split_min_size,
            "To",
        )
        merge_transitions = _snapshot_compare_nontrivial_transitions(
            right_items,
            left_items,
            nontrivial_split_min_size,
            "From",
        )
        split_fragment_count = sum(transition["fragment_count"] for transition in split_transitions)
        split_member_count = sum(transition["member_count"] for transition in split_transitions)
        merge_fragment_count = sum(transition["fragment_count"] for transition in merge_transitions)
        merge_member_count = sum(transition["member_count"] for transition in merge_transitions)
        is_recombination = bool(split_transitions and merge_transitions)
        component_members = left_members | right_members
        component_assertions = [
            assertion
            for assertion in newly_failing_assertions
            if set(assertion.get("member_ids") or []).issubset(component_members)
        ]
        component_left_failure_signatures = {
            item["signature"] for item in left_items
            if item.get("generic_structure_status") == "inconsistent"
        }
        component_right_failure_signatures = {
            item["signature"] for item in right_items
            if item.get("generic_structure_status") == "inconsistent"
        }
        component_left_gap_signatures = {
            item["signature"] for item in left_items
            if item.get("generic_structure_status") == "generic_unknown"
        }
        component_right_gap_signatures = {
            item["signature"] for item in right_items
            if item.get("generic_structure_status") == "generic_unknown"
        }
        components.append({
            "change_type": change_type,
            "change_score": change_score,
            "nontrivial_split_count": len(split_transitions),
            "nontrivial_split_fragment_count": split_fragment_count,
            "nontrivial_split_member_count": split_member_count,
            "nontrivial_merge_count": len(merge_transitions),
            "nontrivial_merge_fragment_count": merge_fragment_count,
            "nontrivial_merge_member_count": merge_member_count,
            "is_recombination": is_recombination,
            "nontrivial_split_min_size": nontrivial_split_min_size,
            "split_transitions": split_transitions,
            "merge_transitions": merge_transitions,
            "left_clique_count": len(left_items),
            "right_clique_count": len(right_items),
            "left_member_count": len(left_members),
            "right_member_count": len(right_members),
            "retained_count": len(retained_members),
            "added_count": len(added_members),
            "removed_count": len(removed_members),
            "newly_failing_assertion_count": len(component_assertions),
            "newly_failing_assertions": component_assertions,
            "generic_structure_new_failure_count": len(
                component_right_failure_signatures - component_left_failure_signatures
            ),
            "generic_structure_destination_failure_count": len(component_right_failure_signatures),
            "generic_structure_destination_gap_count": len(component_right_gap_signatures),
            "left_cliques": sorted(left_items, key=lambda item: (item.get("rank_by_size") or 10**12, item["clique_key"])),
            "right_cliques": sorted(right_items, key=lambda item: (item.get("rank_by_size") or 10**12, item["clique_key"])),
            "_retained_members": retained_members,
            "_added_members": added_members,
            "_removed_members": removed_members,
        })
        if split_transitions:
            nontrivial_split_candidates.append(components[-1])
        if merge_transitions:
            nontrivial_merge_candidates.append(components[-1])
        if is_recombination:
            recombination_candidates.append(components[-1])

    components.sort(
        key=lambda item: (
            -item["generic_structure_new_failure_count"],
            -item["generic_structure_destination_failure_count"],
            -item["generic_structure_destination_gap_count"],
            -item["newly_failing_assertion_count"],
            -int(item["is_recombination"]),
            -int(bool(item["nontrivial_split_count"] or item["nontrivial_merge_count"])),
            -item["change_score"],
            -max(item["left_member_count"], item["right_member_count"]),
            item["change_type"],
        )
    )
    limited_components = components[:limit]
    display_components = limited_components
    sample_ids = sorted({
        member_id
        for component in display_components
        for bucket in ("_added_members", "_removed_members", "_retained_members")
        for member_id in _sample_member_ids(component[bucket], 24)
    })
    display_by_id = _load_metabolite_identifier_display_map(sample_ids)
    for component in display_components:
        _build_snapshot_comparison_review_fields(component, left_key, right_key)
        component["visualization_ids"] = list(dict.fromkeys(
            identifier
            for clique in [*component["right_cliques"], *component["left_cliques"]]
            for identifier in [
                clique.get("representative_id"),
                *(clique.get("generic_member_samples") or []),
                *(clique.get("specific_member_samples") or []),
                *(clique.get("unknown_member_samples") or []),
            ]
            if identifier
        ))[:_METABOLITE_COMPARE_MAX_IDS]
        added_members = component.pop("_added_members")
        removed_members = component.pop("_removed_members")
        retained_members = component.pop("_retained_members")
        component["involved_preview"] = _snapshot_compare_involved_preview(
            added_members | removed_members | retained_members,
            display_by_id,
        )
        component["added_samples"] = _format_compare_member_samples(added_members, display_by_id)
        component["removed_samples"] = _format_compare_member_samples(removed_members, display_by_id)
        component["retained_samples"] = _format_compare_member_samples(retained_members, display_by_id, 12)
        for clique in [*component["left_cliques"], *component["right_cliques"]]:
            clique["member_ids"] = clique["member_ids"][:12]
            clique.pop("member_set", None)
            clique.pop("signature", None)

    pipelines = _list_harmonization_pipelines(limit=100)
    return {
        "left_snapshot": left_snapshot,
        "right_snapshot": right_snapshot,
        "pipeline_context": _snapshot_comparison_pipeline_context(
            left_snapshot,
            right_snapshot,
            pipelines,
        ),
        "left_materialized_clique_count": len(left_cliques),
        "right_materialized_clique_count": len(right_cliques),
        "unchanged_clique_count": len(left_cliques) - len(changed_left),
        "changed_left_clique_count": len(changed_left),
        "changed_right_clique_count": len(changed_right),
        "changed_component_count": len(components),
        "newly_failing_assertion_count": len(newly_failing_assertions),
        "newly_failing_assertions": newly_failing_assertions,
        "generic_structure_validation": generic_validation_comparison,
        "assertion_failure_component_count": sum(
            component["newly_failing_assertion_count"] > 0
            for component in components
        ),
        "nontrivial_split_count": sum(
            not component["is_recombination"]
            for component in nontrivial_split_candidates
        ),
        "nontrivial_merge_count": sum(
            not component["is_recombination"]
            for component in nontrivial_merge_candidates
        ),
        "recombination_count": len(recombination_candidates),
        "nontrivial_split_min_size": nontrivial_split_min_size,
        "display_limit": limit,
        "components": limited_components,
        "snapshots": _list_harmonization_stages(),
        "comparison_kind": "stage",
    }


def _load_metabolite_snapshot_comparison_visualization(
    identifier_ids: List[str],
    left_key: str,
    right_key: str,
) -> dict:
    if not identifier_ids:
        raise ValueError("No changed-clique identifiers were provided.")
    if len(identifier_ids) > _METABOLITE_COMPARE_MAX_IDS:
        raise ValueError(f"Visualize {_METABOLITE_COMPARE_MAX_IDS} or fewer identifiers at a time.")
    if not left_key or not right_key or left_key == right_key:
        raise ValueError("Choose two different stages to visualize.")
    snapshot_union = _load_metabolite_snapshot_union(identifier_ids, [left_key, right_key])
    snapshot_graphs = snapshot_union.get("snapshot_graphs") or []
    return {
        "query_ids": identifier_ids,
        "snapshot_graphs": snapshot_graphs,
        "sankey": _build_metabolite_snapshot_sankey(
            snapshot_graphs,
            [left_key, right_key],
        ),
    }


def _list_metabolite_snapshot_jobs() -> List[dict]:
    with _metabolite_snapshot_jobs_lock:
        jobs = [
            dict(job)
            for job in _metabolite_snapshot_jobs.values()
            if not job.get("housekeeping")
        ]
    return sorted(jobs, key=lambda job: job.get("created_at") or "", reverse=True)[:10]


def _active_harmonization_pipeline_job_keys(jobs: Optional[List[dict]] = None) -> set:
    job_rows = jobs if jobs is not None else _list_metabolite_snapshot_jobs()
    return {
        job.get("pipeline_key")
        for job in job_rows
        if job.get("pipeline_key") and job.get("status") in {"queued", "running"}
        and job.get("action") == "run_pipeline"
    }


def _active_harmonization_delete_pipeline_keys(jobs: Optional[List[dict]] = None) -> set:
    job_rows = jobs if jobs is not None else _list_metabolite_snapshot_jobs()
    return {
        job.get("pipeline_key")
        for job in job_rows
        if job.get("pipeline_key")
        and job.get("action") == "delete_pipeline"
        and job.get("status") in {"queued", "running"}
    }


def _has_active_harmonization_pipeline_mutation(exclude_job_id: Optional[str] = None) -> bool:
    with _metabolite_snapshot_jobs_lock:
        return any(
            job.get("id") != exclude_job_id
            and job.get("action") in {"run_pipeline", "delete_pipeline"}
            and job.get("status") in {"queued", "running"}
            for job in _metabolite_snapshot_jobs.values()
        )


def _enqueue_harmonization_orphan_cleanup(retry_count: int = 0) -> str:
    """Queue one low-priority cleanup pass when no pipeline mutation is waiting."""
    with _metabolite_snapshot_jobs_lock:
        if any(
            job.get("status") in {"queued", "running"}
            and job.get("action") in {"run_pipeline", "delete_pipeline", "cleanup_orphan_stages"}
            for job in _metabolite_snapshot_jobs.values()
        ):
            return ""
        completed_cleanup_job_ids = [
            existing_job_id
            for existing_job_id, job in _metabolite_snapshot_jobs.items()
            if job.get("action") == "cleanup_orphan_stages"
            and job.get("status") in {"complete", "failed"}
        ]
        for completed_job_id in completed_cleanup_job_ids:
            _metabolite_snapshot_jobs.pop(completed_job_id, None)
        job_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
        _metabolite_snapshot_jobs[job_id] = {
            "id": job_id,
            "action": "cleanup_orphan_stages",
            "label": "Retire unused harmonization stage",
            "housekeeping": True,
            "snapshot_key": None,
            "pipeline_key": None,
            "status": "queued",
            "created_at": datetime.now(timezone.utc).isoformat(),
            "started_at": None,
            "completed_at": None,
            "message": None,
            "error": None,
            "retry_count": retry_count,
        }
    thread = threading.Thread(
        target=_run_metabolite_snapshot_job,
        args=(job_id, "cleanup_orphan_stages", {"retry_count": retry_count}),
        daemon=True,
    )
    thread.start()
    return job_id


def _schedule_next_harmonization_orphan_cleanup(
    retry_count: int = 0,
    delay_seconds: float = 0.25,
) -> None:
    """Yield briefly, then enqueue another one-shot pass if the app is idle."""
    timer = threading.Timer(
        delay_seconds,
        _enqueue_harmonization_orphan_cleanup,
        kwargs={"retry_count": retry_count},
    )
    timer.daemon = True
    timer.start()


def _update_metabolite_snapshot_job(job_id: str, **updates) -> None:
    with _metabolite_snapshot_jobs_lock:
        job = _metabolite_snapshot_jobs.get(job_id)
        if job is None:
            return
        job.update(updates)


def _run_metabolite_snapshot_job(job_id: str, action: str, payload: dict) -> None:
    _update_metabolite_snapshot_job(job_id, status="running", started_at=datetime.now(timezone.utc).isoformat())
    cleanup_deleted_stage = False
    cleanup_retry_count = None
    trigger_cleanup = action in {"run_pipeline", "delete_pipeline"}
    try:
        # Stage snapshots are shared across pipelines. Keep all mutations in
        # one process-wide critical section so two jobs cannot materialize or
        # retire the same deterministic stage concurrently.
        with _harmonization_pipeline_mutation_lock:
            if action == "run_pipeline":
                run = _run_harmonization_pipeline(payload["pipeline_key"])
                _update_metabolite_snapshot_job(
                    job_id,
                    status="complete",
                    completed_at=datetime.now(timezone.utc).isoformat(),
                    pipeline_key=payload["pipeline_key"],
                    run_key=run.get("_key"),
                    stage_keys=run.get("stage_keys") or [],
                    message=f"Completed pipeline {run.get('pipeline_name') or payload['pipeline_key']}",
                )
            elif action == "delete_pipeline":
                result = _delete_harmonization_pipeline(get_db("metabolite_harmonization"), payload["pipeline_key"])
                _update_metabolite_snapshot_job(
                    job_id,
                    status="complete",
                    completed_at=datetime.now(timezone.utc).isoformat(),
                    pipeline_key=payload["pipeline_key"],
                    stage_keys=result.get("deleted_stage_keys") or [],
                    message=f"Deleted pipeline {result.get('pipeline', {}).get('name') or payload['pipeline_key']}",
                )
            elif action == "cleanup_orphan_stages":
                if _has_active_harmonization_pipeline_mutation(exclude_job_id=job_id):
                    _update_metabolite_snapshot_job(
                        job_id,
                        status="complete",
                        completed_at=datetime.now(timezone.utc).isoformat(),
                        message="Deferred stage cleanup while pipeline work is queued",
                    )
                else:
                    result = delete_one_orphan_stage(
                        get_db("metabolite_harmonization"),
                        grace_period=_HARMONIZATION_ORPHAN_GRACE_PERIOD,
                    )
                    cleanup_deleted_stage = result is not None
                    _update_metabolite_snapshot_job(
                        job_id,
                        status="complete",
                        completed_at=datetime.now(timezone.utc).isoformat(),
                        stage_keys=[result["stage_key"]] if result else [],
                        message=(
                            f"Retired unused stage {result['stage_key']}"
                            if result else "No unused harmonization stages were ready for cleanup"
                        ),
                    )
            else:
                raise ValueError(f"Unknown harmonization job action: {action}")
    except Exception as exc:
        if action == "cleanup_orphan_stages":
            prior_retry_count = int(payload.get("retry_count") or 0)
            logger.exception(
                "Harmonization orphan-stage cleanup failed (attempt %s of %s)",
                prior_retry_count + 1,
                _HARMONIZATION_CLEANUP_MAX_RETRIES + 1,
            )
            if prior_retry_count < _HARMONIZATION_CLEANUP_MAX_RETRIES:
                cleanup_retry_count = prior_retry_count + 1
        _update_metabolite_snapshot_job(
            job_id,
            status="failed",
            completed_at=datetime.now(timezone.utc).isoformat(),
            error=str(exc),
        )
    finally:
        if cleanup_deleted_stage:
            _schedule_next_harmonization_orphan_cleanup()
        elif cleanup_retry_count is not None:
            _schedule_next_harmonization_orphan_cleanup(
                retry_count=cleanup_retry_count,
                delay_seconds=2.0,
            )
        elif trigger_cleanup:
            _enqueue_harmonization_orphan_cleanup()


def _enqueue_metabolite_snapshot_job(action: str, label: str, payload: dict) -> str:
    pipeline_key = payload.get("pipeline_key")
    active_pipeline_keys = _active_harmonization_pipeline_job_keys()
    deleting_pipeline_keys = _active_harmonization_delete_pipeline_keys()
    if action == "run_pipeline" and pipeline_key in (active_pipeline_keys | deleting_pipeline_keys):
        return ""
    if action == "delete_pipeline" and pipeline_key in (active_pipeline_keys | deleting_pipeline_keys):
        return ""
    job_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    job = {
        "id": job_id,
        "action": action,
        "label": label,
        "snapshot_key": payload.get("snapshot_key"),
        "pipeline_key": payload.get("pipeline_key"),
        "status": "queued",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "started_at": None,
        "completed_at": None,
        "message": None,
        "error": None,
    }
    with _metabolite_snapshot_jobs_lock:
        _metabolite_snapshot_jobs[job_id] = job
    thread = threading.Thread(
        target=_run_metabolite_snapshot_job,
        args=(job_id, action, payload),
        daemon=True,
    )
    thread.start()
    return job_id


def _load_metabolite_snapshot_memberships(
    identifier_ids: List[str],
    snapshot_keys: Optional[List[str]] = None,
) -> List[dict]:
    if not identifier_ids:
        return []
    db = get_db("metabolite_harmonization")
    if (
        not db.has_collection(_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION)
        or not db.has_collection(_HARMONIZED_METABOLITE_COLLECTION)
        or not db.has_collection(_HARMONIZATION_STAGE_COLLECTION)
    ):
        return []
    selected_snapshot_keys = snapshot_keys or []
    stage_filter = "FILTER e.stage_key IN @snapshot_keys" if selected_snapshot_keys else ""
    return list(db.aql.execute(
        f"""
        FOR node IN MetaboliteIdentifier
          FILTER node.id IN @identifier_ids
          FOR e IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
            FILTER e._to == node._id
            {stage_filter}
            LET clique = DOCUMENT(e._from)
            LET snapshot = DOCUMENT("{_HARMONIZATION_STAGE_COLLECTION}", e.stage_key)
            FILTER clique != null AND snapshot != null
            SORT snapshot.created_at DESC, clique.rank_by_size ASC, e.member_id
            RETURN {{
              member_id: e.member_id,
              snapshot_key: snapshot._key,
              snapshot_name: snapshot.name,
              snapshot_created_at: snapshot.created_at,
              rule_ids: snapshot.rule_ids,
              rules: snapshot.rules,
              clique_id: clique.id,
              clique_key: clique._key,
              clique_size: clique.size,
              clique_rank_by_size: clique.rank_by_size,
              representative_id: clique.representative_id,
              sample_member_ids: clique.sample_member_ids
            }}
        """,
        bind_vars={
            "identifier_ids": identifier_ids,
            **({"snapshot_keys": selected_snapshot_keys} if selected_snapshot_keys else {}),
        },
        max_runtime=120,
    ))


def _parse_metabolite_mass(value) -> Optional[float]:
    if value is None:
        return None
    try:
        parsed = float(str(value).strip())
    except (TypeError, ValueError):
        return None
    if math.isnan(parsed):
        return None
    return parsed


def _metabolite_mass_summary(values: List[float]) -> Optional[dict]:
    sorted_values = sorted(values)
    if not sorted_values:
        return None
    midpoint = len(sorted_values) // 2
    if len(sorted_values) % 2:
        median = sorted_values[midpoint]
    else:
        median = (sorted_values[midpoint - 1] + sorted_values[midpoint]) / 2
    return {
        "count": len(sorted_values),
        "min": sorted_values[0],
        "median": median,
        "max": sorted_values[-1],
    }


def _metabolite_member_mass_summary(member_ids: List[str], node_by_member_id: Dict[str, dict]) -> Optional[dict]:
    masses = []
    for member_id in member_ids:
        node = node_by_member_id.get(member_id) or {}
        masses.extend(node.get("masses") or [])
    return _metabolite_mass_summary(masses)


def _metabolite_mass_summary_label(summary: Optional[dict]) -> Optional[str]:
    if not summary:
        return None
    min_mass = summary.get("min")
    median_mass = summary.get("median")
    max_mass = summary.get("max")
    count = summary.get("count")
    if min_mass is None or median_mass is None or max_mass is None:
        return None
    if min_mass == max_mass:
        return f"MW {min_mass:.4g} (n={count})"
    return f"MW {min_mass:.4g}-{max_mass:.4g}; med {median_mass:.4g} (n={count})"


def _metabolite_compact_mass_label(summary: Optional[dict]) -> Optional[str]:
    if not summary:
        return None
    median_mass = summary.get("median")
    if median_mass is None:
        return None
    return f"MW {median_mass:.4g}"


def _metabolite_mw_spread_percent(summary: Optional[dict]) -> Optional[float]:
    if not summary:
        return None
    min_mass = summary.get("min")
    max_mass = summary.get("max")
    if min_mass is None or max_mass is None or min_mass <= 0:
        return None
    return (max_mass - min_mass) / min_mass


def _metabolite_member_mass_values(member_rows: List[dict]) -> List[float]:
    masses = []
    for member in member_rows or []:
        for value in member.get("raw_masses") or []:
            parsed = _parse_metabolite_mass(value)
            if parsed is not None:
                masses.append(parsed)
    return masses


def _metabolite_member_mass_examples(member_rows: List[dict], limit: int = 8) -> List[dict]:
    examples = []
    for member in _kegg_first(member_rows or [], lambda item: item.get("member_id")):
        member_id = member.get("member_id")
        if not member_id:
            continue
        parsed_masses = []
        for value in member.get("raw_masses") or []:
            parsed = _parse_metabolite_mass(value)
            if parsed is not None:
                parsed_masses.append(parsed)
        if not parsed_masses:
            continue
        examples.append({
            "member_id": member_id,
            "masses": sorted(set(parsed_masses)),
        })
        if len(examples) >= limit:
            break
    return examples


def _metabolite_member_mass_profile_examples(
    member_ids: List[str],
    profiles_by_id: Dict[str, dict],
    limit: int = 8,
) -> List[dict]:
    examples = []
    for member_id in _kegg_first(member_ids):
        whole = (profiles_by_id.get(member_id) or {}).get("whole") or {}
        average_masses = sorted(set(whole.get("average") or []))
        monoisotopic_masses = sorted(set(whole.get("monoisotopic") or []))
        if not average_masses and not monoisotopic_masses:
            continue
        examples.append({
            "member_id": member_id,
            "average_masses": average_masses,
            "monoisotopic_masses": monoisotopic_masses,
        })
        if len(examples) >= limit:
            break
    return examples


def _metabolite_mass_channel_labels(channel_results: dict) -> List[str]:
    labels = []
    for channel, label in (
        ("average", "Average MW"),
        ("monoisotopic", "Monoisotopic mass"),
    ):
        summary = (channel_results.get(channel) or {}).get("summary")
        if not summary:
            continue
        if summary["min"] == summary["max"]:
            labels.append(f"{label} {summary['min']:.4g} (n={summary['count']})")
        else:
            labels.append(
                f"{label} {summary['min']:.4g}-{summary['max']:.4g}; "
                f"med {summary['median']:.4g} (n={summary['count']})"
            )
    return labels


def _load_metabolite_identifier_mass_values(
    db,
    record_overlays: Optional[dict[str, dict]] = None,
    chemical_entity_overlays: Optional[dict[str, dict]] = None,
) -> Dict[str, dict]:
    mass_values_by_id: Dict[str, dict] = {}

    def profile(identifier: str) -> dict:
        return mass_values_by_id.setdefault(identifier, {
            "whole": {"average": [], "monoisotopic": []},
            "whole_observations": [],
            "components": [],
        })

    def add_values(target: list, values: Iterable[Any]) -> None:
        for value in values:
            parsed = _parse_metabolite_mass(value)
            if parsed is not None and parsed > 0:
                target.append(parsed)

    if db.has_collection("MetaboliteIdentifier"):
        for row in db.aql.execute(
            """
            FOR d IN MetaboliteIdentifier
              LET chemistry = (
                FOR prop IN d.chem_props || []
                  RETURN {
                    source: prop.source,
                    source_id: prop.source_id,
                    molecular_formula: prop.molecular_formula,
                    average: [prop.mw, prop.calculated_mw],
                    monoisotopic: [prop.monoisotopic_mass, prop.calculated_monoisotopic_mass],
                    components: prop.structure_components || []
                  }
              )
              FILTER LENGTH(chemistry) > 0
              RETURN {id: d.id, chemistry: chemistry}
            """,
            max_runtime=300,
        ):
            target = profile(row["id"])
            overlay = (record_overlays or {}).get(row["id"])
            chemistry_records = (
                overlay.get("chem_props") or []
                if overlay is not None
                else row.get("chemistry") or []
            )
            for chemistry in chemistry_records:
                if overlay is not None:
                    chemistry = {
                        "source": chemistry.get("source"),
                        "source_id": chemistry.get("source_id"),
                        "molecular_formula": chemistry.get("molecular_formula"),
                        "average": [
                            chemistry.get("mw"),
                            chemistry.get("calculated_mw"),
                        ],
                        "monoisotopic": [
                            chemistry.get("monoisotopic_mass"),
                            chemistry.get("calculated_monoisotopic_mass"),
                        ],
                        "components": chemistry.get("structure_components") or [],
                    }
                add_values(target["whole"]["average"], chemistry.get("average") or [])
                add_values(target["whole"]["monoisotopic"], chemistry.get("monoisotopic") or [])
                for channel in ("average", "monoisotopic"):
                    for raw_value in chemistry.get(channel) or []:
                        parsed = _parse_metabolite_mass(raw_value)
                        if parsed is not None and parsed > 0:
                            target["whole_observations"].append({
                                "model_type": "MetaboliteIdentifier",
                                "source": chemistry.get("source"),
                                "source_id": chemistry.get("source_id"),
                                "channel": channel,
                                "value": parsed,
                                "molecular_formula": chemistry.get("molecular_formula"),
                            })
                for component_index, component in enumerate(chemistry.get("components") or []):
                    average = _parse_metabolite_mass(component.get("mw"))
                    monoisotopic = _parse_metabolite_mass(component.get("monoisotopic_mass"))
                    if average is None and monoisotopic is None:
                        continue
                    target["components"].append({
                        "source": chemistry.get("source"),
                        "source_id": chemistry.get("source_id"),
                        "component_index": component_index,
                        "average": average,
                        "monoisotopic": monoisotopic,
                        "molecular_formula": component.get("molecular_formula"),
                        "smiles": component.get("smiles"),
                    })
    if db.has_collection("ChemicalEntity"):
        for row in db.aql.execute(
            """
            FOR d IN ChemicalEntity
              LET chemistry = {
                source: "ChEBI",
                source_id: d.id,
                molecular_formula: d.formula,
                average: [d.mass, d.calculated_mw],
                monoisotopic: [d.monoisotopic_mass, d.calculated_monoisotopic_mass],
                components: d.structure_components || []
              }
              FILTER chemistry.average[0] != null
                OR chemistry.average[1] != null
                OR chemistry.monoisotopic[0] != null
                OR chemistry.monoisotopic[1] != null
                OR LENGTH(chemistry.components) > 0
                OR d.id IN @overlay_ids
              RETURN {id: d.id, chemistry: chemistry}
            """,
            bind_vars={
                "overlay_ids": sorted((chemical_entity_overlays or {}).keys()),
            },
            max_runtime=300,
        ):
            target = profile(row["id"])
            chemical_overlay = (chemical_entity_overlays or {}).get(row["id"])
            chemistry = (
                {
                    "source": "ChEBI",
                    "source_id": chemical_overlay.get("id"),
                    "molecular_formula": chemical_overlay.get("formula"),
                    "average": [
                        chemical_overlay.get("mass"),
                        chemical_overlay.get("calculated_mw"),
                    ],
                    "monoisotopic": [
                        chemical_overlay.get("monoisotopic_mass"),
                        chemical_overlay.get("calculated_monoisotopic_mass"),
                    ],
                    "components": chemical_overlay.get("structure_components") or [],
                }
                if chemical_overlay is not None
                else row.get("chemistry") or {}
            )
            add_values(target["whole"]["average"], chemistry.get("average") or [])
            add_values(
                target["whole"]["monoisotopic"],
                chemistry.get("monoisotopic") or [],
            )
            for channel in ("average", "monoisotopic"):
                for raw_value in chemistry.get(channel) or []:
                    parsed = _parse_metabolite_mass(raw_value)
                    if parsed is not None and parsed > 0:
                        target["whole_observations"].append({
                            "model_type": "ChemicalEntity",
                            "source": chemistry.get("source"),
                            "source_id": chemistry.get("source_id"),
                            "channel": channel,
                            "value": parsed,
                            "molecular_formula": chemistry.get("molecular_formula"),
                        })
            for component_index, component in enumerate(chemistry.get("components") or []):
                average = _parse_metabolite_mass(component.get("mw"))
                monoisotopic = _parse_metabolite_mass(component.get("monoisotopic_mass"))
                if average is None and monoisotopic is None:
                    continue
                target["components"].append({
                    "source": "ChEBI",
                    "source_id": row["id"],
                    "component_index": component_index,
                    "average": average,
                    "monoisotopic": monoisotopic,
                    "molecular_formula": component.get("molecular_formula"),
                    "smiles": component.get("smiles"),
                })
    for value in mass_values_by_id.values():
        for channel in value["whole"]:
            value["whole"][channel] = sorted(set(value["whole"][channel]))
        value["whole_observations"] = sorted(
            {
                json.dumps(item, sort_keys=True, separators=(",", ":")): item
                for item in value["whole_observations"]
            }.values(),
            key=lambda item: (
                item.get("source") or "", item.get("source_id") or "",
                item.get("channel") or "", item.get("value") or 0,
            ),
        )
    return mass_values_by_id


def _cached_record_overlay_mass_values_provider(db):
    """Cache one immutable pipeline run's overlay-aware mass profile."""
    cache = {}

    def provider(record_overlays, chemical_entity_overlays=None):
        if "values" not in cache:
            cache["values"] = _load_metabolite_identifier_mass_values(
                db,
                record_overlays,
                chemical_entity_overlays,
            )
        return cache["values"]

    return provider


def _generic_structure_clique_summary(
    members: List[str],
    classifications: dict[str, Optional[bool]],
    sample_limit: int = 8,
) -> dict:
    members_by_state = {"generic": [], "specific": [], "unknown": []}
    for identifier in members:
        value = classifications.get(identifier)
        state = "generic" if value is True else "specific" if value is False else "unknown"
        members_by_state[state].append(identifier)
    has_generic = bool(members_by_state["generic"])
    has_specific = bool(members_by_state["specific"])
    has_unknown = bool(members_by_state["unknown"])
    if has_generic and has_specific:
        status = "inconsistent"
    elif has_generic and has_unknown:
        status = "generic_unknown"
    else:
        status = "consistent"
    return {
        "status": status,
        "generic_member_count": len(members_by_state["generic"]),
        "specific_member_count": len(members_by_state["specific"]),
        "unknown_member_count": len(members_by_state["unknown"]),
        "generic_member_samples": _kegg_first(members_by_state["generic"])[:sample_limit],
        "specific_member_samples": _kegg_first(members_by_state["specific"])[:sample_limit],
        "unknown_member_samples": _kegg_first(members_by_state["unknown"])[:sample_limit],
    }


def _build_harmonization_stage_generic_structure_validation(
    groups: List[List[str]],
    active_edges: List[dict],
    classifications: dict[str, Optional[bool]],
    rule_enabled: bool,
    limit: int = _METABOLITE_GENERIC_STRUCTURE_WARNING_LIMIT,
) -> dict:
    group_rank_by_id = {
        identifier: rank
        for rank, members in enumerate(groups, start=1)
        for identifier in members
    }
    boundary_edges_by_rank: dict[int, list[dict]] = {}
    actionable_raw_edges = set()
    for edge in active_edges:
        start_id = edge["start_id"]
        end_id = edge["end_id"]
        rank = group_rank_by_id.get(start_id)
        if rank is None or rank != group_rank_by_id.get(end_id):
            continue
        start_value = classifications.get(start_id)
        end_value = classifications.get(end_id)
        if True not in {start_value, end_value} or start_value == end_value:
            continue
        other_value = end_value if start_value is True else start_value
        boundary_kind = "generic_specific" if other_value is False else "generic_unknown"
        boundary = {
            "start_id": start_id,
            "end_id": end_id,
            "start_classification": (
                "generic" if start_value is True else "specific" if start_value is False else "unknown"
            ),
            "end_classification": (
                "generic" if end_value is True else "specific" if end_value is False else "unknown"
            ),
            "boundary_kind": boundary_kind,
            "sources": edge.get("sources") or [],
            "synthetic": bool(edge.get("synthetic")),
            "rule_id": edge.get("rule_id"),
            "assertion_id": edge.get("assertion_id"),
        }
        boundary_edges_by_rank.setdefault(rank, []).append(boundary)
        if boundary_kind == "generic_specific" and not boundary["synthetic"]:
            actionable_raw_edges.add(tuple(sorted((start_id, end_id))))

    warnings = []
    inconsistent_clique_count = 0
    generic_unknown_clique_count = 0
    for rank, members in enumerate(groups, start=1):
        clique_summary = _generic_structure_clique_summary(members, classifications)
        if clique_summary["status"] == "consistent":
            continue
        if clique_summary["status"] == "inconsistent":
            inconsistent_clique_count += 1
        else:
            generic_unknown_clique_count += 1
        boundary_edges = boundary_edges_by_rank.get(rank) or []
        warnings.append({
            "rank_by_size": rank,
            "representative_id": members[0] if members else None,
            "size": len(members),
            **clique_summary,
            "boundary_edge_count": len(boundary_edges),
            "boundary_edges": boundary_edges[:12],
            "comparison_ids": " ".join(
                list(dict.fromkeys([
                    *clique_summary["generic_member_samples"],
                    *clique_summary["specific_member_samples"],
                    *clique_summary["unknown_member_samples"],
                ]))[:_METABOLITE_COMPARE_MAX_IDS]
            ),
        })
    warnings.sort(key=lambda item: (
        0 if item["status"] == "inconsistent" else 1,
        -item["boundary_edge_count"],
        -item["size"],
        item["rank_by_size"],
    ))
    return {
        "computed": True,
        "rule_enabled": rule_enabled,
        "inconsistent_clique_count": inconsistent_clique_count,
        "generic_unknown_clique_count": generic_unknown_clique_count,
        "actionable_raw_edge_count": len(actionable_raw_edges),
        "warning_count": inconsistent_clique_count + generic_unknown_clique_count,
        "display_limit": limit,
        "warnings": warnings[:limit],
    }


def _build_harmonization_stage_mw_validation(
    groups: List[List[str]],
    mass_values_by_id: Dict[str, Any],
    threshold: float = _METABOLITE_MW_SPREAD_WARNING_THRESHOLD,
    limit: int = _METABOLITE_MW_SPREAD_WARNING_LIMIT,
) -> dict:
    warnings = []
    for rank, members in enumerate(groups, start=1):
        display_members = _kegg_first(members)
        profiles_by_id = {
            member_id: (
                mass_values_by_id.get(member_id)
                if isinstance(mass_values_by_id.get(member_id), dict)
                else {
                    "whole": {
                        "average": [],
                        "monoisotopic": mass_values_by_id.get(member_id) or [],
                    },
                    "components": [],
                }
            )
            for member_id in members
        }
        assessment = assess_mass_profiles(
            members,
            profiles_by_id,
            spread_threshold=threshold,
        )
        if assessment is None:
            continue
        summary = assessment["mass_summary"]
        spread = max(
            result.get("spread") or 0
            for result in assessment["channel_results"].values()
        )
        finding_metadata = mass_validation_finding_metadata(
            members,
            profiles_by_id,
            assessment,
            spread_threshold=threshold,
        )
        warnings.append({
            **finding_metadata,
            "rank_by_size": rank,
            "representative_id": members[0] if members else None,
            "size": len(members),
            "sample_member_ids": display_members[:12],
            "mass_summary": summary,
            "mass_summaries": assessment["mass_summaries"],
            "mass_label": " · ".join(
                _metabolite_mass_channel_labels(assessment["channel_results"])
            ),
            "mass_channel_labels": _metabolite_mass_channel_labels(
                assessment["channel_results"]
            ),
            "spread": spread,
            "spread_percent": spread * 100,
            "severity": assessment["severity"],
            "reason": assessment["reason"],
            "channel_results": assessment["channel_results"],
            "component_matches": assessment["component_matches"],
            "member_mass_examples": _metabolite_member_mass_profile_examples(
                display_members, profiles_by_id, limit=len(display_members),
            ),
            "source_mass_observations": [
                {"member_id": member_id, **observation}
                for member_id in display_members
                for observation in (
                    profiles_by_id.get(member_id, {}).get("whole_observations") or []
                )
            ],
            "comparison_ids": " ".join(display_members[:_METABOLITE_COMPARE_MAX_IDS]),
        })
    warnings.sort(key=lambda item: (
        0 if item["severity"] == "error" else 1,
        -(item["spread"] or 0),
        -(item["size"] or 0),
        item.get("representative_id") or "",
    ))
    error_count = sum(item["severity"] == "error" for item in warnings)
    review_warning_count = sum(item["severity"] == "warning" for item in warnings)
    finding_index = [
        {
            **{
                key: item.get(key)
                for key in (
                "finding_id",
                "anchor_id",
                "member_ids",
                "evidence_fingerprint",
                "evidence_fingerprint_version",
                "severity",
                )
            },
            "legacy_evidence_fingerprint": _legacy_mw_finding_metadata(
                {
                    "version": MW_VALIDATION_VERSION,
                    "threshold": threshold,
                },
                item,
            )["evidence_fingerprint"],
        }
        for item in warnings
    ]
    return {
        "computed": True,
        "version": MW_VALIDATION_VERSION,
        "threshold": threshold,
        "threshold_percent": threshold * 100,
        "warning_count": len(warnings),
        "issue_count": len(warnings),
        "error_count": error_count,
        "review_warning_count": review_warning_count,
        "display_limit": limit,
        "finding_index": finding_index,
        "warnings": warnings[:limit],
    }


def _build_harmonization_stage_denylist_validation(
    groups: List[List[str]],
    denylist_pairs: set[tuple[str, str]],
    rule_enabled: bool,
    limit: int = _METABOLITE_DENYLIST_STILL_MERGED_LIMIT,
    suppressed_identifier_ids: Optional[set[str]] = None,
) -> dict:
    suppressed_identifier_ids = suppressed_identifier_ids or set()
    group_index_by_id: Dict[str, int] = {}
    for rank, members in enumerate(groups, start=1):
        for member_id in members:
            group_index_by_id[member_id] = rank

    pairs_by_rank: Dict[int, List[tuple[str, str]]] = {}
    suppressed_incident_pairs = {
        pair
        for pair in denylist_pairs
        if pair[0] in suppressed_identifier_ids or pair[1] in suppressed_identifier_ids
    }
    for left, right in sorted(denylist_pairs - suppressed_incident_pairs):
        left_rank = group_index_by_id.get(left)
        right_rank = group_index_by_id.get(right)
        if left_rank is None or right_rank is None or left_rank != right_rank:
            continue
        pairs_by_rank.setdefault(left_rank, []).append((left, right))

    total_pair_count = sum(len(pairs) for pairs in pairs_by_rank.values())

    warnings = []
    for rank, pairs in pairs_by_rank.items():
        members = groups[rank - 1]
        flagged_ids = sorted({identifier for pair in pairs for identifier in pair})
        remaining_members = [member_id for member_id in members if member_id not in flagged_ids]
        ordered_member_ids = _kegg_first(flagged_ids + remaining_members)
        warnings.append({
            "rank_by_size": rank,
            "size": len(members),
            "pair_count": len(pairs),
            "pairs": [{"left_id": left, "right_id": right} for left, right in pairs],
            "sample_member_ids": ordered_member_ids[:12],
            "comparison_ids": " ".join(ordered_member_ids[:_METABOLITE_COMPARE_MAX_IDS]),
        })
    warnings.sort(key=lambda item: (-(item["pair_count"] or 0), -(item["size"] or 0), item["rank_by_size"]))
    return {
        "computed": True,
        "rule_enabled": rule_enabled,
        "denylist_pair_count": len(denylist_pairs),
        "evaluated_denylist_pair_count": len(denylist_pairs) - len(suppressed_incident_pairs),
        "suppressed_incident_pair_count": len(suppressed_incident_pairs),
        "warning_count": total_pair_count,
        "affected_clique_count": len(pairs_by_rank),
        "display_limit": limit,
        "warnings": warnings[:limit],
    }


def _harmonization_stage_denylist_validation_from_doc(
    stage: dict,
    limit: int = _METABOLITE_DENYLIST_STILL_MERGED_LIMIT,
) -> dict:
    existing = (stage.get("validation") or {}).get("denylist_still_merged")
    if isinstance(existing, dict):
        return _prioritize_kegg_validation_samples(existing)
    return {
        "computed": False,
        "rule_enabled": False,
        "denylist_pair_count": 0,
        "evaluated_denylist_pair_count": 0,
        "suppressed_incident_pair_count": 0,
        "warning_count": 0,
        "affected_clique_count": 0,
        "display_limit": limit,
        "warnings": [],
    }


def _harmonization_denylist_review_from_stage(stage: dict, rank_by_size: int) -> Optional[dict]:
    validation = _harmonization_stage_denylist_validation_from_doc(stage)
    warning = next(
        (
            item
            for item in validation.get("warnings") or []
            if item.get("rank_by_size") == rank_by_size
        ),
        None,
    )
    if warning is None:
        return None
    return {
        "stage_key": stage.get("_key"),
        "rank_by_size": rank_by_size,
        "pairs": warning.get("pairs") or [],
    }


def _normalize_harmonization_denylist_pair(value: str) -> tuple[str, str]:
    identifiers = _parse_metabolite_identifier_query(value)
    if len(identifiers) != 2:
        raise ValueError("denylist_pair must contain exactly two prefixed metabolite identifiers.")
    normalized = [_normalize_ramp_denylist_identifier(identifier) for identifier in identifiers]
    if any(identifier is None for identifier in normalized) or normalized[0] == normalized[1]:
        raise ValueError("denylist_pair must contain two distinct prefixed metabolite identifiers.")
    return tuple(sorted(normalized))


def _build_harmonization_denylist_pair_review(
    pair: tuple[str, str],
    result: dict,
    curation_state: dict,
) -> dict:
    stage_states = []
    focus = None
    for section in result.get("snapshot_graphs") or []:
        memberships = {}
        for clique in section.get("cliques") or []:
            for identifier in pair:
                if identifier in (clique.get("member_ids") or []):
                    memberships[identifier] = clique
        together = len(memberships) == 2 and len({
            clique.get("clique_key") for clique in memberships.values()
        }) == 1
        stage_state = {
            "stage_key": section.get("snapshot_key"),
            "stage_name": section.get("snapshot_name"),
            "together": together,
            "found_count": len(memberships),
            "memberships": [
                {
                    "identifier": identifier,
                    "clique_key": memberships.get(identifier, {}).get("clique_key"),
                    "clique_size": memberships.get(identifier, {}).get("clique_size"),
                    "is_singleton": memberships.get(identifier, {}).get("is_singleton_bucket", False),
                }
                for identifier in pair
            ],
        }
        stage_states.append(stage_state)
        if focus is None and together:
            shared = memberships[pair[0]]
            focus = {
                "stage_key": section.get("snapshot_key"),
                "clique_key": shared.get("clique_key"),
            }

    suppressed_endpoint_ids = sorted(
        set(pair) & set(curation_state.get("suppressed_identifier_ids") or set())
    )
    decision = (curation_state.get("pair_decisions") or {}).get(pair)
    active_action = (curation_state.get("pair_states") or {}).get(pair)
    if focus is None:
        for stage_state in stage_states:
            membership = next(
                (item for item in stage_state["memberships"] if item.get("clique_key")),
                None,
            )
            if membership:
                focus = {
                    "stage_key": stage_state["stage_key"],
                    "clique_key": membership["clique_key"],
                }
                break
    return {
        "pair": {"left_id": pair[0], "right_id": pair[1]},
        "active_action": active_action,
        "is_active_removal": active_action == "remove_edge",
        "retain_allowed": active_action == "remove_edge" and not suppressed_endpoint_ids,
        "suppressed_endpoint_ids": suppressed_endpoint_ids,
        "decision": decision,
        "stage_states": stage_states,
        "focus": focus,
    }


def _metabolite_query_id_labels(query_ids: List[str], max_ids: int = 4) -> List[str]:
    if not query_ids:
        return []
    labels = [f"ID {query_id}" for query_id in query_ids[:max_ids]]
    if len(query_ids) > max_ids:
        labels.append(f"+{len(query_ids) - max_ids} more query IDs")
    return labels


def _metabolite_identifier_source_linkout(identifier: str) -> Optional[dict]:
    prefix, separator, local_id = str(identifier or "").partition(":")
    if not separator or not local_id:
        return None
    linkouts = {
        "CHEBI": (
            "ChEBI",
            f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={url_quote(identifier, safe='')}",
        ),
        "HMDB": ("HMDB", f"https://hmdb.ca/metabolites/{url_quote(local_id, safe='')}"),
        "KEGG.COMPOUND": ("KEGG Compound", f"https://www.kegg.jp/entry/{url_quote(local_id, safe='')}"),
        "LIPIDMAPS": ("LIPID MAPS", f"https://www.lipidmaps.org/databases/lmsd/{url_quote(local_id, safe='')}"),
        "REFMET": (
            "RefMet",
            f"https://www.metabolomicsworkbench.org/databases/refmet/refmet_details.php?REFMET_ID={url_quote(local_id, safe='')}",
        ),
        "PUBCHEM.COMPOUND": ("PubChem", f"https://pubchem.ncbi.nlm.nih.gov/compound/{url_quote(local_id, safe='')}"),
        "ChemSpider": ("ChemSpider", f"https://www.chemspider.com/Chemical-Structure.{url_quote(local_id, safe='')}.html"),
        "SwissLipids": ("SwissLipids", f"https://www.swisslipids.org/#/entity/{url_quote(local_id, safe=':')}"),
        "Wikidata": ("Wikidata", f"https://www.wikidata.org/wiki/{url_quote(local_id, safe='')}"),
        "CAS": ("CAS Common Chemistry", f"https://commonchemistry.cas.org/detail?cas_rn={url_quote(local_id, safe='')}"),
        "DRUGBANK": ("DrugBank", f"https://go.drugbank.com/drugs/{url_quote(local_id, safe='')}"),
        "InChIKey": (
            "PubChem InChIKey lookup",
            f"https://pubchem.ncbi.nlm.nih.gov/compound/{url_quote(local_id, safe='-')}",
        ),
    }
    linkout = linkouts.get(prefix)
    if linkout is None:
        return None
    label, url = linkout
    return {"label": label, "url": url}


def _parse_metabolite_snapshot_key_filter(value: str) -> List[str]:
    snapshot_keys = []
    seen = set()
    for token in re.split(r"[\s,;|]+", value or ""):
        snapshot_key = token.strip()
        if not snapshot_key or snapshot_key in seen:
            continue
        seen.add(snapshot_key)
        snapshot_keys.append(snapshot_key)
    return snapshot_keys


def _default_metabolite_snapshot_key_filter() -> List[str]:
    stages = _list_harmonization_stages()
    if not stages:
        return []
    latest_stage_key = stages[-1].get("_key")
    return [latest_stage_key] if latest_stage_key else []


def _default_metabolite_denylist_pair_stage_filter() -> List[str]:
    """Show a current graph baseline and final stage for pair-level evidence."""
    stages = _list_harmonization_stages()
    if not stages:
        return []
    latest = next((
        stage
        for stage in reversed(stages)
        if {"apply_curations", "ignore_ramp_mapping_denylist"}
        & set(stage.get("rule_ids") or [])
    ), stages[-1])

    def graph_fingerprint(stage: dict) -> tuple:
        summary = stage.get("summary") or {}
        return (
            stage.get("engine_version"),
            summary.get("database"),
            summary.get("metabolite_identifier_count"),
            summary.get("equivalence_edge_count"),
        )

    latest_graph_fingerprint = graph_fingerprint(latest)
    matching_baselines = [
        stage
        for stage in stages
        if int(stage.get("stage_index") or 0) == 0
        and not (stage.get("rule_ids") or [])
        and graph_fingerprint(stage) == latest_graph_fingerprint
    ]
    selected = [
        stage.get("_key")
        for stage in [matching_baselines[-1] if matching_baselines else None, latest]
        if stage and stage.get("_key")
    ]
    return list(dict.fromkeys(selected))


def _metabolite_display_member_ids_by_clique(
    rows: List[dict],
    identifier_ids: List[str],
    direct_neighbor_rows: List[dict],
    limit: int = _METABOLITE_INTERACTIVE_CLIQUE_NODE_LIMIT,
) -> Dict[tuple, List[str]]:
    direct_neighbor_ids_by_stage: Dict[str, List[str]] = {}
    for neighbor in direct_neighbor_rows:
        direct_neighbor_ids_by_stage.setdefault(neighbor["stage_key"], []).append(neighbor["id"])
    member_ids_by_clique: Dict[tuple, List[str]] = {}
    for row in rows:
        member_ids_by_clique.setdefault(
            (row["snapshot_key"], row["clique_key"]), []
        ).append(row["member_id"])
    display_member_ids_by_clique: Dict[tuple, List[str]] = {}
    for clique_ref, clique_member_ids in member_ids_by_clique.items():
        snapshot_key, _clique_key = clique_ref
        clique_member_id_set = set(clique_member_ids)
        prioritized_ids = [
            *[member_id for member_id in identifier_ids if member_id in clique_member_id_set],
            *[
                member_id
                for member_id in direct_neighbor_ids_by_stage.get(snapshot_key, [])
                if member_id in clique_member_id_set
            ],
            *clique_member_ids,
        ]
        display_member_ids_by_clique[clique_ref] = list(dict.fromkeys(prioritized_ids))[:limit]
    return display_member_ids_by_clique


def _load_metabolite_snapshot_union(
    identifier_ids: List[str],
    snapshot_keys: Optional[List[str]] = None,
) -> dict:
    selected_snapshot_keys = snapshot_keys or []
    memberships = _load_metabolite_snapshot_memberships(
        identifier_ids,
        selected_snapshot_keys,
    )
    db = get_db("metabolite_harmonization")
    existing_query_ids = sorted(list(db.aql.execute(
        """
        FOR node IN MetaboliteIdentifier
          FILTER node.id IN @ids
          COLLECT existing_id = node.id
          SORT existing_id
          RETURN existing_id
        """,
        bind_vars={"ids": identifier_ids},
        max_runtime=120,
    )))
    if not memberships and not existing_query_ids:
        return {
            "memberships": [],
            "member_ids": [],
            "snapshot_graphs": [],
            "sankey": {"nodes": [], "links": []},
        }
    clique_keys = []
    for membership in memberships:
        clique_key = membership.get("clique_key")
        if clique_key and clique_key not in clique_keys:
            clique_keys.append(clique_key)
    clique_member_ids = []
    clique_vertex_ids = []
    if clique_keys:
        clique_vertex_ids = [
            f"{_HARMONIZED_METABOLITE_COLLECTION}/{clique_key}"
            for clique_key in clique_keys
        ]
        clique_member_ids = list(db.aql.execute(
            f"""
            FOR e IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
              FILTER e._from IN @clique_vertex_ids
              FILTER DOCUMENT(e._to) != null
              COLLECT member_id = e.member_id
              SORT member_id
              RETURN member_id
            """,
            bind_vars={"clique_vertex_ids": clique_vertex_ids},
            max_runtime=120,
        ))
    union_member_ids = sorted(set(clique_member_ids) | set(existing_query_ids))
    if not union_member_ids:
        return {
            "memberships": memberships,
            "member_ids": [],
            "snapshot_graphs": [],
            "sankey": {"nodes": [], "links": []},
        }
    rows = list(db.aql.execute(
        f"""
        FOR e IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
          FILTER e._from IN @clique_vertex_ids
          LET node = DOCUMENT(e._to)
          LET clique = DOCUMENT(e._from)
          LET snapshot = DOCUMENT("{_HARMONIZATION_STAGE_COLLECTION}", e.stage_key)
          FILTER node != null AND clique != null AND snapshot != null
          LET names = node.names || []
          SORT snapshot.created_at, clique.rank_by_size, e.member_index
          RETURN {{
            member_id: e.member_id,
            member_label: LENGTH(names) > 0 ? names[0].value : e.member_id,
            member_prefix: node.prefix,
            name_count: LENGTH(node.names || []),
            synonym_count: LENGTH(node.synonyms || []),
            chem_prop_count: LENGTH(node.chem_props || []),
            snapshot_key: snapshot._key,
            snapshot_name: snapshot.name,
            snapshot_created_at: snapshot.created_at,
            rules: snapshot.rules,
            clique_id: clique.id,
            clique_key: clique._key,
            clique_size: clique.size,
            clique_rank_by_size: clique.rank_by_size
          }}
        """,
        bind_vars={"clique_vertex_ids": clique_vertex_ids},
        max_runtime=120,
    ))
    detail_stage_keys = selected_snapshot_keys or sorted({
        row["snapshot_key"] for row in rows
    })
    direct_neighbor_rows = list(db.aql.execute(
        f"""
        LET query_nodes = (
          FOR node IN MetaboliteIdentifier
            FILTER node.id IN @identifier_ids
            RETURN node._id
        )
        FOR stage_key IN @snapshot_keys
          LET outgoing = (
            FOR handle IN query_nodes
              FOR e IN {_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION}
                FILTER e._from == handle
                FILTER e.stage_key == stage_key
                LIMIT @neighbor_limit
                RETURN e._to
          )
          LET incoming = (
            FOR handle IN query_nodes
              FOR e IN {_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION}
                FILTER e._to == handle
                FILTER e.stage_key == stage_key
                LIMIT @neighbor_limit
                RETURN e._from
          )
          FOR handle IN UNION_DISTINCT(outgoing, incoming)
            LET neighbor = DOCUMENT(handle)
            FILTER neighbor != null
            SORT neighbor.id
            RETURN {{id: neighbor.id, stage_key}}
        """,
        bind_vars={
            "identifier_ids": identifier_ids,
            "neighbor_limit": _METABOLITE_INTERACTIVE_CLIQUE_NODE_LIMIT,
            "snapshot_keys": detail_stage_keys,
        },
        max_runtime=120,
    )) if db.has_collection(_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION) else []
    display_member_ids_by_clique = _metabolite_display_member_ids_by_clique(
        rows,
        identifier_ids,
        direct_neighbor_rows,
    )
    detail_member_ids = list(existing_query_ids)
    detail_member_id_set = set(detail_member_ids)
    for displayed_ids in display_member_ids_by_clique.values():
        for member_id in displayed_ids:
            if member_id not in detail_member_id_set:
                detail_member_id_set.add(member_id)
                detail_member_ids.append(member_id)
    curatable_property_names = [
        definition.property_name
        for definition in curatable_property_definitions("MetaboliteIdentifier")
    ]
    node_rows = list(db.aql.execute(
        """
        FOR node IN MetaboliteIdentifier
          FILTER node.id IN @member_ids
          LET chemical_entity = STARTS_WITH(node.id, "CHEBI:") ? DOCUMENT("ChemicalEntity", node.id) : null
          LET names = node.names || []
          LET name_values = UNIQUE(
            FOR name IN names
              FILTER name.value != null AND name.value != ""
              RETURN name.value
          )
          LET chemical_entity_names = UNIQUE(APPEND(
            chemical_entity != null && chemical_entity.name != null && chemical_entity.name != "" ? [chemical_entity.name] : [],
            (
              FOR synonym IN chemical_entity != null ? (chemical_entity.synonyms || []) : []
                FILTER synonym.value != null AND synonym.value != ""
                RETURN synonym.value
            )
          ))
          LET synonym_values = UNIQUE(
            FOR synonym IN node.synonyms || []
              FILTER synonym.value != null AND synonym.value != ""
              RETURN synonym.value
          )
          LET display_names = UNIQUE(APPEND(
            LENGTH(name_values) > 0 ? name_values : chemical_entity_names,
            LENGTH(name_values) > 0 ? chemical_entity_names : synonym_values
          ))
          LET raw_masses = UNIQUE(FLATTEN(
            APPEND(
              (FOR prop IN node.chem_props || [] RETURN [
                prop.mw,
                prop.monoisotopic_mass,
                prop.calculated_mw,
                prop.calculated_monoisotopic_mass
              ]),
              chemical_entity == null ? [] : [[chemical_entity.mass, chemical_entity.monoisotopic_mass]]
            )
          ))
          LET chem_prop_summaries = (
            FOR prop IN node.chem_props || []
              RETURN {
                source: prop.source,
                source_id: prop.source_id,
                name: prop.common_name || prop.iupac_name,
                formula: prop.molecular_formula,
                mw: prop.mw,
                monoisotopic_mass: prop.monoisotopic_mass,
                calculated_mw: prop.calculated_mw,
                calculated_monoisotopic_mass: prop.calculated_monoisotopic_mass,
                structure_components: prop.structure_components || [],
                structure_calculation_input_field: prop.structure_calculation_input_field,
                structure_calculation_method: prop.structure_calculation_method,
                structure_calculation_method_version: prop.structure_calculation_method_version,
                structure_calculation_error: prop.structure_calculation_error,
                inchi_key_prefix: prop.inchi_key_prefix,
                inchi_key: prop.inchi_key,
                derived_inchi_key_prefix: prop.derived_inchi_key_prefix,
                derived_inchi_key: prop.derived_inchi_key,
                derived_inchi_key_input_field: prop.derived_inchi_key_input_field,
                derived_inchi_key_method: prop.derived_inchi_key_method,
                derived_inchi_key_method_version: prop.derived_inchi_key_method_version,
                derived_inchi_key_error: prop.derived_inchi_key_error,
                smiles: prop.iso_smiles || prop.isomeric_smiles || prop.canonical_smiles,
                inchi: prop.inchi
              }
          )
          RETURN {
            id: node.id,
            editable_properties: KEEP(node, @curatable_property_names),
            label: LENGTH(display_names) > 0 ? display_names[0] : node.id,
            names: SLICE(display_names, 0, 12),
            prefix: node.prefix,
            name_count: LENGTH(node.names || []),
            synonym_count: LENGTH(node.synonyms || []),
            chem_prop_count: LENGTH(node.chem_props || []),
            raw_masses: raw_masses,
            chem_props: SLICE(chem_prop_summaries, 0, 8),
            generic_structure_evidence: (
              FOR prop IN chem_prop_summaries
                RETURN KEEP(prop, "source", "source_id", "smiles", "formula", "inchi")
            ),
            formulas: UNIQUE(
              FOR prop IN chem_prop_summaries
                FILTER prop.formula != null AND prop.formula != ""
                RETURN prop.formula
            ),
            smiles: UNIQUE(
              FOR prop IN chem_prop_summaries
                FILTER prop.smiles != null AND prop.smiles != ""
                RETURN prop.smiles
            ),
            inchi_keys: UNIQUE(
              FOR prop IN chem_prop_summaries
                FILTER prop.inchi_key != null AND prop.inchi_key != ""
                RETURN prop.inchi_key
            ),
            derived_inchi_keys: UNIQUE(
              FOR prop IN chem_prop_summaries
                FILTER prop.derived_inchi_key != null AND prop.derived_inchi_key != ""
                RETURN prop.derived_inchi_key
            ),
            chemical_entity: chemical_entity == null ? null : {
              id: chemical_entity.id,
              name: chemical_entity.name,
              formula: chemical_entity.formula,
              mass: chemical_entity.mass,
              monoisotopic_mass: chemical_entity.monoisotopic_mass,
              inchi_key: chemical_entity.inchi_key,
              smiles: chemical_entity.smiles,
              inchi: chemical_entity.inchi
            }
          }
        """,
        bind_vars={
            "member_ids": detail_member_ids,
            "curatable_property_names": curatable_property_names,
        },
        max_runtime=120,
    ))
    mass_rows = list(db.aql.execute(
        """
        FOR node IN MetaboliteIdentifier
          FILTER node.id IN @member_ids
          LET chemical_entity = STARTS_WITH(node.id, "CHEBI:") ? DOCUMENT("ChemicalEntity", node.id) : null
          RETURN {
            id: node.id,
            raw_masses: UNIQUE(FLATTEN(
              APPEND(
                (FOR prop IN node.chem_props || [] RETURN [
                  prop.mw,
                  prop.monoisotopic_mass,
                  prop.calculated_mw,
                  prop.calculated_monoisotopic_mass
                ]),
                chemical_entity == null ? [] : [[chemical_entity.mass, chemical_entity.monoisotopic_mass]]
              )
            ))
          }
        """,
        bind_vars={"member_ids": union_member_ids},
        max_runtime=120,
    ))
    mass_values_by_id = {
        row["id"]: [
            mass
            for mass in (_parse_metabolite_mass(value) for value in row.get("raw_masses") or [])
            if mass is not None
        ]
        for row in mass_rows
    }
    active_member_ids_by_stage: Dict[str, set[str]] = {}
    if db.has_collection(_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION):
        active_member_ids_by_stage = {
            stage_key: set() for stage_key in detail_stage_keys
        }
        active_rows = db.aql.execute(
            f"""
            FOR chunk IN {_HARMONIZATION_STAGE_ACTIVE_IDENTIFIER_CHUNK_COLLECTION}
              FILTER chunk.stage_key IN @stage_keys
              FOR member_id IN chunk.identifier_ids || []
                FILTER member_id IN @member_ids
                RETURN DISTINCT {{stage_key: chunk.stage_key, member_id: member_id}}
            """,
            bind_vars={
                "stage_keys": detail_stage_keys,
                "member_ids": union_member_ids,
            },
            max_runtime=120,
        )
        for row in active_rows:
            active_member_ids_by_stage.setdefault(row["stage_key"], set()).add(row["member_id"])
    display_selections = [
        {
            "stage_key": snapshot_key,
            "clique_vertex_id": f"{_HARMONIZED_METABOLITE_COLLECTION}/{clique_key}",
            "member_ids": member_ids,
        }
        for (snapshot_key, clique_key), member_ids in display_member_ids_by_clique.items()
    ]
    mapping_edge_rows = list(db.aql.execute(
        f"""
        FOR selection IN @display_selections
          LET member_nodes = (
            FOR membership IN {_HARMONIZED_METABOLITE_MEMBER_EDGE_COLLECTION}
              FILTER membership._from == selection.clique_vertex_id
              FILTER membership.member_id IN selection.member_ids
              COLLECT handle = membership._to, id = membership.member_id
              RETURN {{handle, id}}
          )
          LET member_id_by_handle = ZIP(member_nodes[*].handle, member_nodes[*].id)
          FOR node IN member_nodes
            FOR e IN {_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION}
              FILTER e._from == node.handle
              FILTER e.stage_key == selection.stage_key
              FILTER HAS(member_id_by_handle, e._to)
              RETURN {{
                id: e.id || e._key,
                key: e._key,
                stage_key: e.stage_key,
                start_id: e.start_id,
                end_id: e.end_id,
                sources: e.sources || [],
                details: e.details || [],
                raw_edge_key: e.raw_edge_key,
                raw_edge_id: e.raw_edge_id,
                synthetic: e.synthetic || false,
                rule_id: e.rule_id,
                algorithm_version: e.algorithm_version
              }}
        """,
        bind_vars={"display_selections": display_selections},
        max_runtime=120,
    )) if db.has_collection(_HARMONIZATION_STAGE_EVIDENCE_EDGE_COLLECTION) else []
    highlighted_id_set = set(identifier_ids)
    sections_by_snapshot: Dict[str, dict] = {}
    node_by_member_id: Dict[str, dict] = {
        row["id"]: row
        for row in node_rows
    }
    try:
        display_curation_state = _load_metabolite_edge_removal_curations()
    except Exception:
        display_curation_state = {
            "record_property_decisions": {},
            "record_property_state_available": False,
            "record_property_state_error": (
                "Published record-property curation state is unavailable. "
                "Classification editing is disabled."
            ),
            "suppression_state_available": False,
            "suppression_state_error": (
                "Published suppression state is unavailable. Record participation editing is disabled."
            ),
        }
    for node in node_by_member_id.values():
        node["record_url"] = (
            f"{_root_path()}/db/metabolite_harmonization/collection/"
            f"MetaboliteIdentifier/doc/{url_quote(str(node.get('id') or ''), safe='')}"
        )
        chemical_entity = node.get("chemical_entity") or {}
        if chemical_entity.get("id"):
            chemical_entity["record_url"] = (
                f"{_root_path()}/db/metabolite_harmonization/collection/"
                f"ChemicalEntity/doc/"
                f"{url_quote(str(chemical_entity['id']), safe='')}"
            )
        node["generic_structure"] = _metabolite_generic_structure_classification(
            node,
            display_curation_state,
        )
        node["harmonization_participation"] = _metabolite_harmonization_participation(
            node.get("id"),
            display_curation_state,
        )
        source_linkout = _metabolite_identifier_source_linkout(node.get("id"))
        if source_linkout:
            node["source_url"] = source_linkout["url"]
            node["source_label"] = source_linkout["label"]
        node["masses"] = [
            mass
            for mass in (_parse_metabolite_mass(value) for value in node.get("raw_masses") or [])
            if mass is not None
        ]
        mass_summary = _metabolite_mass_summary(node["masses"])
        mass_label = _metabolite_mass_summary_label(mass_summary)
        compact_mass_label = _metabolite_compact_mass_label(mass_summary)
        node["mass_summary"] = mass_summary
        node["node_label"] = "\n".join([
            part
            for part in (node.get("id"), compact_mass_label)
            if part
        ])
        node["node_detail_label"] = "\n".join([part for part in [node.get("label") or node.get("id"), mass_label] if part])
    for snapshot in _list_harmonization_stages():
        if selected_snapshot_keys and snapshot["_key"] not in selected_snapshot_keys:
            continue
        sections_by_snapshot.setdefault(snapshot["_key"], {
            "snapshot_key": snapshot["_key"],
            "snapshot_name": snapshot.get("name"),
            "snapshot_created_at": snapshot.get("created_at"),
            "rules": snapshot.get("rules") or [],
            "cliques_by_key": {},
        })
    for row in rows:
        member_node = node_by_member_id.setdefault(row["member_id"], {
            "id": row["member_id"],
            "label": row["member_label"],
            "node_label": row["member_label"],
            "names": [row["member_label"]] if row.get("member_label") else [],
            "prefix": row["member_prefix"],
            "name_count": row["name_count"],
            "synonym_count": row["synonym_count"],
            "chem_prop_count": row["chem_prop_count"],
        })
        member_node.setdefault("masses", mass_values_by_id.get(row["member_id"], []))
        snapshot_key = row["snapshot_key"]
        section = sections_by_snapshot.setdefault(snapshot_key, {
            "snapshot_key": snapshot_key,
            "snapshot_name": row["snapshot_name"],
            "snapshot_created_at": row["snapshot_created_at"],
            "rules": row.get("rules") or [],
            "cliques_by_key": {},
        })
        clique = section["cliques_by_key"].setdefault(row["clique_key"], {
            "snapshot_key": snapshot_key,
            "snapshot_name": row["snapshot_name"],
            "snapshot_created_at": row["snapshot_created_at"],
            "rules": row.get("rules") or [],
            "clique_id": row["clique_id"],
            "clique_key": row["clique_key"],
            "clique_size": row["clique_size"],
            "clique_rank_by_size": row["clique_rank_by_size"],
            "member_ids": [],
        })
        if row["member_id"] not in clique["member_ids"]:
            clique["member_ids"].append(row["member_id"])
    union_member_id_set = set(union_member_ids)
    for section in sections_by_snapshot.values():
        represented_member_ids = {
            member_id
            for clique in section["cliques_by_key"].values()
            for member_id in clique["member_ids"]
        }
        missing_member_ids = union_member_id_set - represented_member_ids
        active_ids_for_stage = active_member_ids_by_stage.get(section["snapshot_key"])
        inactive_member_ids = sorted(
            missing_member_ids - active_ids_for_stage
            if active_ids_for_stage is not None
            else set()
        )
        singleton_member_ids = sorted(missing_member_ids - set(inactive_member_ids))
        if singleton_member_ids:
            singleton_bucket_key = f"{section['snapshot_key']}__singleton_bucket"
            section["cliques_by_key"][singleton_bucket_key] = {
                "snapshot_key": section["snapshot_key"],
                "snapshot_name": section["snapshot_name"],
                "snapshot_created_at": section["snapshot_created_at"],
                "rules": section["rules"],
                "clique_id": f"MetaboliteHarmonizationClique:{singleton_bucket_key}",
                "clique_key": singleton_bucket_key,
                "clique_size": len(singleton_member_ids),
                "clique_rank_by_size": None,
                "member_ids": singleton_member_ids,
                "is_singleton_bucket": True,
            }
        if inactive_member_ids:
            inactive_bucket_key = f"{section['snapshot_key']}__inactive_bucket"
            section["cliques_by_key"][inactive_bucket_key] = {
                "snapshot_key": section["snapshot_key"],
                "snapshot_name": section["snapshot_name"],
                "snapshot_created_at": section["snapshot_created_at"],
                "rules": section["rules"],
                "clique_id": f"MetaboliteHarmonizationInactive:{inactive_bucket_key}",
                "clique_key": inactive_bucket_key,
                "clique_size": len(inactive_member_ids),
                "clique_rank_by_size": None,
                "member_ids": inactive_member_ids,
                "is_inactive_bucket": True,
            }

    snapshot_sections = []
    for section in sections_by_snapshot.values():
        cliques = []
        for clique in sorted(
            section["cliques_by_key"].values(),
            key=lambda item: (
                item.get("is_inactive_bucket", False),
                item.get("is_singleton_bucket", False),
                item["clique_rank_by_size"] if item["clique_rank_by_size"] is not None else 10**12,
                item["clique_key"],
            ),
        ):
            member_name_samples = []
            member_id_samples = []
            member_name_by_id = {}
            member_mass_by_id = {}
            seen_member_names = set()
            for member_id in clique["member_ids"]:
                node = node_by_member_id.get(member_id, {"id": member_id, "label": member_id, "node_label": member_id, "names": []})
                member_mass_by_id[member_id] = node.get("masses") or []
                display_names = node.get("names") or []
                if display_names:
                    member_name_by_id[member_id] = display_names[0]
                for display_name in display_names:
                    if display_name and display_name not in seen_member_names:
                        seen_member_names.add(display_name)
                        if len(member_name_samples) < 25:
                            member_name_samples.append(display_name)
                if len(member_id_samples) < 25:
                    member_id_samples.append(member_id)
            selected_query_member_ids = [
                query_id
                for query_id in identifier_ids
                if query_id in clique["member_ids"]
            ]
            clique_member_id_set = set(clique["member_ids"])
            clique_mapping_edges = [
                edge
                for edge in mapping_edge_rows
                if (
                    (not edge.get("stage_key") or edge.get("stage_key") == section["snapshot_key"])
                    and edge.get("start_id") in clique_member_id_set
                    and edge.get("end_id") in clique_member_id_set
                )
            ]
            selected_query_member_id_set = set(selected_query_member_ids)
            elements = []
            displayed_member_ids = display_member_ids_by_clique.get(
                (section["snapshot_key"], clique["clique_key"]),
                clique["member_ids"][:_METABOLITE_INTERACTIVE_CLIQUE_NODE_LIMIT],
            )
            displayed_member_id_set = set(displayed_member_ids)
            for member_id in displayed_member_ids:
                node = node_by_member_id.get(member_id, {"id": member_id, "label": member_id, "node_label": member_id, "prefix": "unknown"})
                classes = f"metabolite-id-node prefix-{str(node.get('prefix') or 'unknown').lower().replace('.', '-')}"
                effective_generic = (node.get("generic_structure") or {}).get("effective")
                if effective_generic is True:
                    classes += " generic-structure-generic"
                    node["generic_structure_state"] = "generic"
                elif effective_generic is False:
                    classes += " generic-structure-specific"
                    node["generic_structure_state"] = "specific / non-generic"
                else:
                    classes += " generic-structure-unknown"
                    node["generic_structure_state"] = "unknown"
                if member_id in highlighted_id_set:
                    classes += " selected-query"
                if clique.get("is_inactive_bucket"):
                    classes += " harmonization-inactive"
                elements.append({
                    "data": {
                        **node,
                        "selected": member_id in highlighted_id_set,
                        "is_active_in_stage": not clique.get("is_inactive_bucket", False),
                        "stage_key": section["snapshot_key"],
                    },
                    "classes": classes,
                })
            displayed_mapping_edges = [
                edge
                for edge in clique_mapping_edges
                if (
                    edge.get("start_id") in displayed_member_id_set
                    and edge.get("end_id") in displayed_member_id_set
                )
            ]
            displayed_mapping_edges.sort(key=lambda edge: (
                not (
                    edge.get("start_id") in selected_query_member_id_set
                    or edge.get("end_id") in selected_query_member_id_set
                ),
                edge.get("key") or edge.get("id") or "",
            ))
            rendered_mapping_edges = displayed_mapping_edges[
                :_METABOLITE_INTERACTIVE_CLIQUE_EDGE_LIMIT
            ]
            for edge in rendered_mapping_edges:
                start_id = edge.get("start_id")
                end_id = edge.get("end_id")
                sources = edge.get("sources") or []
                evidence_display = _metabolite_harmonization_evidence_display(edge)
                elements.append({
                    "data": {
                        "id": f"metabolite-equivalence::{clique['clique_key']}::{edge.get('key') or edge.get('id')}",
                        "source": start_id,
                        "target": end_id,
                        "label": evidence_display["label"],
                        "kind": "MetaboliteIdentifierMappingEdge",
                        "start_id": start_id,
                        "end_id": end_id,
                        "sources": ", ".join(sources),
                        "source_count": len(sources),
                        "detail_count": len(edge.get("details") or []),
                        "raw_edge_key": edge.get("raw_edge_key"),
                        "raw_edge_id": edge.get("raw_edge_id"),
                        "synthetic": bool(edge.get("synthetic")),
                        "rule_id": evidence_display["rule_id"],
                        "rule_label": evidence_display["rule_label"],
                        "algorithm_version": edge.get("algorithm_version"),
                        "snapshot": section["snapshot_name"],
                        "snapshot_key": section["snapshot_key"],
                    },
                    "classes": evidence_display["class_name"],
                })
            cliques.append({
                **clique,
                "member_name_samples": member_name_samples,
                "member_id_samples": member_id_samples,
                "member_name_by_id": member_name_by_id,
                "member_mass_by_id": member_mass_by_id,
                "selected_query_member_ids": selected_query_member_ids,
                "mass_summary": _metabolite_member_mass_summary(clique["member_ids"], node_by_member_id),
                "display_member_count": len(displayed_member_ids),
                "display_edge_count": len(rendered_mapping_edges),
                "display_center_ids": selected_query_member_ids,
                "display_truncated": len(clique["member_ids"]) > len(displayed_member_ids),
                "edge_display_truncated": len(displayed_mapping_edges) > len(rendered_mapping_edges),
                "elements": elements,
            })
        snapshot_sections.append({
            **{key: value for key, value in section.items() if key != "cliques_by_key"},
            "cliques": cliques,
            "materialized_clique_count": sum(1 for clique in cliques if not clique.get("is_singleton_bucket")),
            "singleton_bucket_count": sum(1 for clique in cliques if clique.get("is_singleton_bucket")),
            "singleton_bucket_member_count": sum(
                len(clique.get("member_ids", []))
                for clique in cliques
                if clique.get("is_singleton_bucket")
            ),
            "inactive_bucket_count": sum(1 for clique in cliques if clique.get("is_inactive_bucket")),
            "inactive_bucket_member_count": sum(
                len(clique.get("member_ids", []))
                for clique in cliques
                if clique.get("is_inactive_bucket")
            ),
        })

    sankeys = _build_metabolite_snapshot_sankeys(
        snapshot_sections,
        selected_snapshot_keys,
    )
    return {
        "memberships": memberships,
        "member_ids": union_member_ids,
        "snapshot_graphs": snapshot_sections,
        "sankeys": sankeys,
    }


def _build_metabolite_snapshot_sankeys(
    snapshot_sections: List[dict],
    selected_snapshot_keys: Optional[List[str]] = None,
) -> List[dict]:
    """Build stage-flow Sankeys, preserving explicit comparisons across pipelines."""
    sections_by_key = {
        section["snapshot_key"]: section
        for section in snapshot_sections
    }
    ordered_selected_keys = [
        stage_key
        for stage_key in selected_snapshot_keys or []
        if stage_key in sections_by_key
    ]
    if len(ordered_selected_keys) > 1:
        selected_sections = [sections_by_key[stage_key] for stage_key in ordered_selected_keys]
        return [{
            "pipeline_key": "selected-stage-comparison",
            "pipeline_name": "Selected stage comparison",
            "stage_count": len(selected_sections),
            "sankey": _build_metabolite_snapshot_sankey(
                selected_sections,
                ordered_selected_keys,
            ),
        }]

    pipelines = _list_harmonization_pipelines(limit=100)
    sankeys = []
    assigned_stage_keys = set()
    for pipeline in pipelines:
        latest_run = (pipeline.get("runs") or [None])[0] or {}
        ordered_stage_keys = [
            stage.get("_key")
            for stage in sorted(latest_run.get("stages") or [], key=_stage_sort_key)
            if stage.get("_key") in sections_by_key
        ]
        if not ordered_stage_keys:
            continue
        assigned_stage_keys.update(ordered_stage_keys)
        sections = [sections_by_key[key] for key in ordered_stage_keys]
        sankeys.append({
            "pipeline_key": pipeline.get("_key"),
            "pipeline_name": pipeline.get("name") or pipeline.get("_key"),
            "stage_count": len(ordered_stage_keys),
            "sankey": _build_metabolite_snapshot_sankey(sections, ordered_stage_keys),
        })

    unassigned_sections = [
        section for section in snapshot_sections
        if section.get("snapshot_key") not in assigned_stage_keys
    ]
    if unassigned_sections:
        unassigned_sections.sort(key=lambda section: section.get("snapshot_created_at") or "")
        sankeys.append({
            "pipeline_key": "unassigned",
            "pipeline_name": "Stages not assigned to a pipeline",
            "stage_count": len(unassigned_sections),
            "sankey": _build_metabolite_snapshot_sankey(
                unassigned_sections,
                [section["snapshot_key"] for section in unassigned_sections],
            ),
        })
    return sankeys


def _build_metabolite_snapshot_sankey(
    snapshot_sections: List[dict],
    ordered_snapshot_keys: Optional[List[str]] = None,
) -> dict:
    nodes = []
    links_by_pair: Dict[tuple[str, str], dict] = {}
    node_ids = set()
    if ordered_snapshot_keys:
        section_by_key = {section["snapshot_key"]: section for section in snapshot_sections}
        ordered_sections = [
            section_by_key[key]
            for key in ordered_snapshot_keys
            if key in section_by_key
        ]
    else:
        ordered_sections = sorted(snapshot_sections, key=lambda section: section.get("snapshot_created_at") or "")
    for stage_index, section in enumerate(ordered_sections):
        for clique in section.get("cliques", []):
            if clique.get("is_inactive_bucket"):
                continue
            node_id = f"{section['snapshot_key']}::{clique['clique_key']}"
            if node_id not in node_ids:
                node_ids.add(node_id)
                mass_label = _metabolite_mass_summary_label(clique.get("mass_summary"))
                detail_label = mass_label or ("Single ID" if clique.get("is_singleton_bucket") else None)
                node_name_parts = [section.get("snapshot_name") or section["snapshot_key"]]
                if detail_label:
                    node_name_parts.append(detail_label)
                node_name_parts.extend(_metabolite_query_id_labels(clique.get("selected_query_member_ids") or []))
                nodes.append({
                    "id": node_id,
                    "name": "\\n".join(node_name_parts),
                    "snapshot": section.get("snapshot_name"),
                    "snapshot_key": section.get("snapshot_key"),
                    "clique_key": clique.get("clique_key"),
                    "clique_size": clique.get("clique_size"),
                    "is_singleton_bucket": clique.get("is_singleton_bucket", False),
                    "name_samples": clique.get("member_name_samples", [])[:25],
                    "member_id_samples": clique.get("member_id_samples", [])[:25],
                    "selected_query_member_ids": clique.get("selected_query_member_ids", []),
                    "mass_summary": clique.get("mass_summary"),
                    "stageIndex": stage_index,
                })
    clique_lookup_by_section = {
        section["snapshot_key"]: {
            f"{section['snapshot_key']}::{clique['clique_key']}": clique
            for clique in section.get("cliques", [])
        }
        for section in ordered_sections
    }
    for left, right in zip(ordered_sections, ordered_sections[1:]):
        left_by_member = {}
        right_by_member = {}
        for clique in left.get("cliques", []):
            if clique.get("is_inactive_bucket"):
                continue
            for member_id in clique.get("member_ids", []):
                left_by_member[member_id] = f"{left['snapshot_key']}::{clique['clique_key']}"
        for clique in right.get("cliques", []):
            if clique.get("is_inactive_bucket"):
                continue
            for member_id in clique.get("member_ids", []):
                right_by_member[member_id] = f"{right['snapshot_key']}::{clique['clique_key']}"
        for member_id in sorted(set(left_by_member) & set(right_by_member)):
            pair = (left_by_member[member_id], right_by_member[member_id])
            link = links_by_pair.setdefault(pair, {
                "source": pair[0],
                "target": pair[1],
                "value": 0,
                "member_ids": [],
                "member_names": [],
                "masses": [],
            })
            link["value"] += 1
            source_clique = clique_lookup_by_section.get(left["snapshot_key"], {}).get(pair[0], {})
            link["masses"].extend(source_clique.get("member_mass_by_id", {}).get(member_id, []))
            if len(link["member_ids"]) < 25:
                link["member_ids"].append(member_id)
                member_name = source_clique.get("member_name_by_id", {}).get(member_id, member_id)
                if member_name != member_id:
                    link["member_names"].append(member_name)
    links = []
    for link in links_by_pair.values():
        link["mass_summary"] = _metabolite_mass_summary(link.pop("masses", []))
        links.append(link)
    return {"nodes": nodes, "links": links}


def _load_metabolite_identifier_qa(
    identifier: str,
    snapshot_keys: Optional[List[str]] = None,
    include_snapshot_union: bool = True,
) -> dict:
    if not identifier:
        raise HTTPException(status_code=400, detail="Provide a metabolite identifier.")
    db = get_db("metabolite_harmonization")
    query = """
    LET metabolite_doc = FIRST(
      FOR d IN MetaboliteIdentifier
        FILTER d.id == @identifier
        LIMIT 1
        RETURN d
    )
    LET metabolite = metabolite_doc == null ? null : UNSET(metabolite_doc, "_id", "_key", "_rev")
    LET chemical_entity_doc = FIRST(
      FOR d IN ChemicalEntity
        FILTER d.id == @identifier
        LIMIT 1
        RETURN d
    )
    LET chemical_entity = chemical_entity_doc == null ? null : UNSET(chemical_entity_doc, "_id", "_key", "_rev")
    LET mapping_edges = (
      FOR e IN MetaboliteIdentifierMappingEdge
        FILTER e.start_id == @identifier OR e.end_id == @identifier
        SORT e.start_id, e.end_id
        LIMIT 500
        LET neighbor_id = e.start_id == @identifier ? e.end_id : e.start_id
        LET neighbor = FIRST(
          FOR n IN MetaboliteIdentifier
            FILTER n.id == neighbor_id
            LIMIT 1
            RETURN KEEP(n, "id", "prefix", "names", "synonyms", "sources")
        )
        RETURN {
          start_id: e.start_id,
          end_id: e.end_id,
          neighbor_id: neighbor_id,
          neighbor: neighbor,
          details: e.details,
          sources: e.sources
        }
    )
    LET mapping_edge_count = LENGTH(
      FOR e IN MetaboliteIdentifierMappingEdge
        FILTER e.start_id == @identifier OR e.end_id == @identifier
        RETURN 1
    )
    LET chebi_bridge_edges = (
      FOR e IN ChebiChemicalEntityMetaboliteIdentifierEdge
        FILTER e.start_id == @identifier OR e.end_id == @identifier
        LIMIT 50
        RETURN UNSET(e, "_id", "_key", "_rev", "_from", "_to")
    )
    RETURN {
      query_id: @identifier,
      found: metabolite != null,
      metabolite: metabolite,
      chemical_entity: chemical_entity,
      mapping_edges: mapping_edges,
      mapping_edge_count: mapping_edge_count,
      mapping_edges_truncated: mapping_edge_count > LENGTH(mapping_edges),
      chebi_bridge_edges: chebi_bridge_edges
    }
    """
    cursor = db.aql.execute(query, bind_vars={"identifier": identifier}, max_runtime=120)
    rows = list(cursor)
    if not rows:
        raise HTTPException(status_code=404, detail=f"No result for {identifier}")
    result = rows[0]
    if include_snapshot_union:
        snapshot_union = _load_metabolite_snapshot_union([identifier], snapshot_keys)
        result["snapshot_memberships"] = snapshot_union["memberships"]
        result["snapshot_union_ids"] = snapshot_union["member_ids"]
        result["snapshot_graphs"] = snapshot_union["snapshot_graphs"]
        result["snapshot_sankeys"] = snapshot_union["sankeys"]
    return result


def _parse_metabolite_identifier_query(value: str) -> List[str]:
    tokens = (
        value.replace(",", " ")
        .replace("|", " ")
        .replace(";", " ")
        .replace("\n", " ")
        .split()
    )
    ids = []
    seen = set()
    for token in tokens:
        identifier = token.strip()
        if identifier and identifier not in seen:
            ids.append(identifier)
            seen.add(identifier)
    return ids


def _load_metabolite_identifier_summaries(identifiers: List[str]) -> List[dict]:
    """Load only the identifier fields rendered by the HTML query page, in one AQL call."""
    if not identifiers:
        return []
    db = get_db("metabolite_harmonization")
    return list(db.aql.execute(
        """
        FOR identifier IN @identifiers
          LET metabolite_doc = FIRST(
            FOR d IN MetaboliteIdentifier
              FILTER d.id == identifier
              LIMIT 1
              RETURN d
          )
          LET chemical_entity_doc = FIRST(
            FOR d IN ChemicalEntity
              FILTER d.id == identifier
              LIMIT 1
              RETURN d
          )
          RETURN {
            query_id: identifier,
            found: metabolite_doc != null,
            metabolite: metabolite_doc == null ? null : {
              id: metabolite_doc.id,
              prefix: metabolite_doc.prefix,
              names: SLICE(metabolite_doc.names || [], 0, 25)
            },
            chemical_entity: chemical_entity_doc == null ? null : {
              id: chemical_entity_doc.id,
              name: chemical_entity_doc.name
            }
          }
        """,
        bind_vars={"identifiers": identifiers},
        max_runtime=120,
    ))


def _load_metabolite_identifier_qa_many(
    value: str,
    snapshot_keys: Optional[List[str]] = None,
    include_identifier_details: bool = True,
) -> dict:
    query_ids = _parse_metabolite_identifier_query(value)
    if not query_ids:
        raise HTTPException(status_code=400, detail="Provide at least one metabolite identifier.")
    if len(query_ids) > _METABOLITE_COMPARE_MAX_IDS:
        raise HTTPException(
            status_code=400,
            detail=f"Compare {_METABOLITE_COMPARE_MAX_IDS} or fewer identifiers at a time.",
        )

    if include_identifier_details:
        results = [
            _load_metabolite_identifier_qa(identifier, snapshot_keys, include_snapshot_union=False)
            for identifier in query_ids
        ]
    else:
        results = _load_metabolite_identifier_summaries(query_ids)
    found_results = [result for result in results if result.get("found")]
    missing_results = [result for result in results if not result.get("found")]
    preview_ids = {
        (result.get("metabolite") or {}).get("id")
        for result in found_results
        if (result.get("metabolite") or {}).get("id")
    }
    preview_display_by_id = {}
    for result in found_results:
        metabolite = result.get("metabolite") or {}
        identifier = metabolite.get("id")
        if not identifier:
            continue
        names = metabolite.get("names") or []
        preferred_name = next(
            (name.get("value") for name in names if name.get("value")),
            None,
        )
        chemical_entity = result.get("chemical_entity") or {}
        preview_display_by_id[identifier] = {
            "label": preferred_name or chemical_entity.get("name") or identifier,
        }
    combined_result = {
        "query_id": " ".join(query_ids),
        "query_ids": query_ids,
        "found": bool(found_results),
        "results": results,
        "found_results": found_results,
        "missing_results": missing_results,
        "query_preview": _snapshot_compare_involved_preview(
            preview_ids,
            preview_display_by_id,
        ),
    }
    snapshot_union = _load_metabolite_snapshot_union(query_ids, snapshot_keys)
    combined_result["snapshot_memberships"] = snapshot_union["memberships"]
    combined_result["snapshot_union_ids"] = snapshot_union["member_ids"]
    combined_result["snapshot_graphs"] = snapshot_union["snapshot_graphs"]
    combined_result["snapshot_sankeys"] = snapshot_union["sankeys"]
    combined_result["selected_snapshot_keys"] = snapshot_keys or []
    return combined_result


def _slugify_mysql_source(value: str) -> str:
    slug = "".join(ch.lower() if ch.isalnum() else "-" for ch in value).strip("-")
    while "--" in slug:
        slug = slug.replace("--", "-")
    return slug or "mysql"


def _register_mysql_source(source_id: str, creds: dict, label: Optional[str] = None):
    _mysql_sources[source_id] = {
        "credentials": creds or {},
        "label": label or source_id,
    }


def _get_mysql_source(source_id: str) -> dict:
    source = _mysql_sources.get(source_id)
    if source:
        return source
    if source_id == "default" and _mysql_credentials:
        return {"credentials": _mysql_credentials, "label": "default"}
    raise KeyError(f"Unknown MySQL source: {source_id}")


def _get_parquet_buffer(file_ref: str):
    """Fetch a parquet file from object storage and return (BytesIO, size_bytes).

    file_ref must be an s3:// URI produced by the ETL pipeline.
    Returns (None, None) if the file cannot be fetched.
    """
    if not file_ref.startswith("s3://"):
        return None, None
    try:
        import io

        without_prefix = file_ref[len("s3://"):]
        bucket, key = without_prefix.split("/", 1)
        credentials_options = []
        if _parquet_storage_credentials:
            credentials_options.append(("parquet storage", _parquet_storage_credentials))
        if _object_storage_credentials:
            credentials_options.append(("object storage", _object_storage_credentials))
        if not credentials_options:
            return None, None

        errors = []
        for label, credentials in credentials_options:
            try:
                storage = _storage_from_credentials(credentials, use_internal_url=True)
                response = storage.client().get_object(Bucket=bucket, Key=key)
                data = response["Body"].read()
                return io.BytesIO(data), response["ContentLength"]
            except Exception as exc:
                errors.append(f"{label}: {exc}")
        raise RuntimeError("; ".join(errors))
    except Exception as e:
        raise RuntimeError(f"Failed to fetch {file_ref} from object storage: {e}") from e


def get_client() -> ArangoClient:
    global _client
    if _client is None:
        url = _credentials.get("internal_url") or _credentials.get("url", "http://localhost:8529")
        _client = ArangoClient(hosts=url, request_timeout=120, verify_override=False)
    return _client


def get_db(name: str):
    client = get_client()
    return client.db(name, username=_credentials.get("user", "root"),
                     password=_credentials.get("password", "password"))


def get_sys_db():
    client = get_client()
    return client.db("_system", username=_credentials.get("user", "root"),
                     password=_credentials.get("password", "password"))


def get_mysql_engine(source_id: str = "default") -> Optional[Engine]:
    """Get a MySQL engine (no specific database) for listing databases."""
    try:
        source = _get_mysql_source(source_id)
    except KeyError:
        return None

    credentials = source["credentials"]
    if not credentials:
        return None

    cache_key = f"{source_id}::_root"
    if cache_key not in _mysql_db_engines:
        host = credentials.get("url", "localhost")
        port = credentials.get("port", 3306)
        user = credentials.get("user", "root")
        password = url_quote(credentials.get("password", ""), safe="")
        _mysql_db_engines[cache_key] = create_engine(
            f"mysql+pymysql://{user}:{password}@{host}:{port}",
            pool_pre_ping=True, pool_size=2,
        )
    return _mysql_db_engines[cache_key]


def get_mysql_db_engine(db_name: str, source_id: str = "default") -> Engine:
    """Get a MySQL engine scoped to a specific database."""
    source = _get_mysql_source(source_id)
    credentials = source["credentials"]
    cache_key = f"{source_id}::{db_name}"
    if cache_key not in _mysql_db_engines:
        host = credentials.get("url", "localhost")
        port = credentials.get("port", 3306)
        user = credentials.get("user", "root")
        password = url_quote(credentials.get("password", ""), safe="")
        _mysql_db_engines[cache_key] = create_engine(
            f"mysql+pymysql://{user}:{password}@{host}:{port}/{db_name}",
            pool_pre_ping=True, pool_size=2,
        )
    return _mysql_db_engines[cache_key]


def get_mysql_inspector(db_name: str, source_id: str = "default"):
    """Return cached schema metadata for a database.

    Cached as a plain dict so we don't hold live Inspector objects.
    Call invalidate_mysql_inspector(db_name) to force a refresh.
    """
    cache_key = f"{source_id}::{db_name}"
    if cache_key not in _mysql_inspector_cache:
        engine = get_mysql_db_engine(db_name, source_id=source_id)
        insp = sa_inspect(engine)
        table_names = insp.get_table_names()
        meta = {}
        for tbl in table_names:
            meta[tbl] = {
                "columns": insp.get_columns(tbl),
                "pk": insp.get_pk_constraint(tbl).get("constrained_columns", []),
                "fks": insp.get_foreign_keys(tbl),
            }
        _mysql_inspector_cache[cache_key] = meta
    return _mysql_inspector_cache[cache_key]


def invalidate_mysql_inspector(db_name: str, source_id: str = "default"):
    """Drop the cached schema for a database so the next request re-fetches it."""
    _mysql_inspector_cache.pop(f"{source_id}::{db_name}", None)


# ── Jinja2 filters ──────────────────────────────────────────────────────────

def json_pretty(value):
    try:
        return json.dumps(value, indent=2, default=str)
    except Exception:
        return str(value)


templates.env.filters["json_pretty"] = json_pretty


def _get_graph_views(db) -> dict:
    if not db.has_collection("metadata_store"):
        return {}
    try:
        store = db.collection("metadata_store")
        doc = store.get("graph_views")
        if not doc:
            return {}
        return doc.get("value", {}).get("views", {}) or {}
    except Exception:
        return {}


def _parse_iso_datetime(value: Optional[str]) -> Optional[datetime]:
    if not value or not isinstance(value, str):
        return None
    try:
        dt = datetime.fromisoformat(value.replace("Z", "+00:00"))
        if dt.tzinfo is not None:
            dt = dt.astimezone(timezone.utc).replace(tzinfo=None)
        return dt
    except Exception:
        return None


def _format_duration_seconds(total_seconds: float | int | None) -> str:
    if total_seconds is None:
        return ""
    total_seconds = int(total_seconds)
    hours, remainder = divmod(total_seconds, 3600)
    minutes, seconds = divmod(remainder, 60)
    if hours:
        return f"{hours}h {minutes}m {seconds}s"
    if minutes:
        return f"{minutes}m {seconds}s"
    return f"{seconds}s"


def _format_elapsed(started_at: Optional[str], finished_at: Optional[str] = None) -> str:
    start_dt = _parse_iso_datetime(started_at)
    if start_dt is None:
        return ""
    end_dt = _parse_iso_datetime(finished_at) or datetime.utcnow()
    return _format_duration_seconds(max((end_dt - start_dt).total_seconds(), 0))


def _summarize_checkpoint(checkpoint: dict) -> dict:
    adapters = checkpoint.get("adapters") or {}
    rows = []
    completed = 0
    running = 0
    failed = 0
    total = 0
    total_records = 0

    for adapter_name, metadata in sorted(
        adapters.items(),
        key=lambda item: (item[1].get("adapter_position") or 10**9, item[0])
    ):
        status = metadata.get("status", "unknown")
        if status == "completed":
            completed += 1
        elif status == "running":
            running += 1
        elif status == "failed":
            failed += 1
        total = max(total, metadata.get("adapter_total") or 0)
        total_records += metadata.get("records_written") or 0
        rows.append({
            "name": adapter_name,
            "status": status,
            "adapter_position": metadata.get("adapter_position"),
            "adapter_total": metadata.get("adapter_total"),
            "records_written": metadata.get("records_written"),
            "started_at": metadata.get("started_at"),
            "completed_at": metadata.get("completed_at"),
            "failed_at": metadata.get("failed_at"),
            "error_message": metadata.get("error_message"),
            "elapsed": _format_elapsed(
                metadata.get("started_at"),
                metadata.get("completed_at") or metadata.get("failed_at"),
            ),
        })

    if not total:
        total = len(rows)

    run_id = checkpoint.get("run_id")
    run_name = (run_id or checkpoint.get("_key") or "Unknown run").split("/")[-1]
    run_name = run_name.replace(".yaml", "").replace("_", " ")

    return {
        "run_id": run_id,
        "run_label": run_name,
        "last_updated": checkpoint.get("last_updated"),
        "adapters": rows,
        "completed": completed,
        "running": running,
        "failed": failed,
        "known_total": total,
        "seen_total": len(rows),
        "total_records": total_records,
        "is_active": running > 0,
        "status": (
            "failed" if failed > 0 else
            "running" if running > 0 else
            "completed" if completed > 0 and completed >= total else
            "partial" if completed > 0 else
            "unknown"
        ),
    }


def _get_build_status(db) -> Optional[dict]:
    if not db.has_collection("metadata_store"):
        return None
    try:
        cursor = db.aql.execute("""
            FOR d IN metadata_store
                FILTER d.type == "etl_checkpoint"
                SORT d.last_updated DESC
                RETURN d
        """)
        checkpoints = list(cursor)
    except Exception:
        return None

    if not checkpoints:
        return None

    run_summaries = [_summarize_checkpoint(checkpoint) for checkpoint in checkpoints]

    etl_meta = None
    try:
        cursor = db.aql.execute('RETURN DOCUMENT("metadata_store", "etl_metadata").value')
        results = list(cursor)
        if results and results[0]:
            etl_meta = results[0]
    except Exception:
        etl_meta = None

    latest_checkpoint = checkpoints[0]
    checkpoint_updated_dt = _parse_iso_datetime(latest_checkpoint.get("last_updated"))
    etl_run_dt = _parse_iso_datetime((etl_meta or {}).get("run_date"))
    any_failed = any(summary["failed"] > 0 for summary in run_summaries)
    any_running = any(summary["running"] > 0 for summary in run_summaries)
    any_partial = any(
        summary["completed"] > 0 and summary["completed"] < summary["known_total"]
        for summary in run_summaries
    )
    all_complete = bool(run_summaries) and all(
        summary["completed"] >= summary["known_total"] and summary["known_total"] > 0
        for summary in run_summaries
    )

    if any_failed:
        overall_status = "failed"
        overall_message = "One or more ETL runs have failed adapters."
    elif any_running:
        overall_status = "running"
        overall_message = "One or more ETL runs are still running."
    elif all_complete:
        if etl_run_dt and (checkpoint_updated_dt is None or etl_run_dt >= checkpoint_updated_dt):
            overall_status = "completed"
            overall_message = "ETL and post-processing completed."
        else:
            overall_status = "post_processing"
            overall_message = "Adapters completed. Post-processing or cleanup is still running."
    elif any_partial:
        overall_status = "partial"
        overall_message = "Some ETL runs are only partially complete."
    else:
        overall_status = "unknown"
        overall_message = "Build state is not yet clear from metadata."

    return {
        "run_id": latest_checkpoint.get("run_id"),
        "last_updated": latest_checkpoint.get("last_updated"),
        "etl_run_date": (etl_meta or {}).get("run_date"),
        "runs": run_summaries,
        "completed": sum(summary["completed"] for summary in run_summaries),
        "running": sum(summary["running"] for summary in run_summaries),
        "failed": sum(summary["failed"] for summary in run_summaries),
        "known_total": sum(summary["known_total"] for summary in run_summaries),
        "seen_total": sum(summary["seen_total"] for summary in run_summaries),
        "total_records": sum(summary["total_records"] for summary in run_summaries),
        "is_active": any_running,
        "overall_status": overall_status,
        "overall_message": overall_message,
    }


def _get_collection_schema_entry(db, coll_name: str) -> dict:
    if not db.has_collection("metadata_store"):
        return {}
    try:
        store = db.collection("metadata_store")
        doc = store.get("collection_schemas")
        if not doc:
            return {}
        return (doc.get("collections") or {}).get(coll_name, {}) or {}
    except Exception:
        return {}


def _edge_endpoint_pairs(db, edge_collection: str) -> List[dict]:
    if not db.has_collection(edge_collection):
        return []
    query = """
    FOR e IN @@edge_collection
      COLLECT from_collection = SPLIT(e._from, "/")[0],
              to_collection = SPLIT(e._to, "/")[0]
      WITH COUNT INTO count
      SORT from_collection, to_collection
      RETURN {
        from_collection: from_collection,
        to_collection: to_collection,
        count: count
      }
    """
    try:
        return list(db.aql.execute(query, bind_vars={"@edge_collection": edge_collection}))
    except Exception:
        return []


def _edge_definitions_with_endpoint_pairs(db) -> List[dict]:
    if not db.has_graph("graph"):
        return []
    graph = db.graph("graph")
    edge_defs = graph.edge_definitions()
    enriched_defs = []
    for edge_def in edge_defs:
        enriched_def = dict(edge_def)
        enriched_def["endpoint_pairs"] = _edge_endpoint_pairs(db, edge_def["edge_collection"])
        enriched_defs.append(enriched_def)
    return enriched_defs


def _get_collection_facet_metadata(db, coll_name: str) -> dict:
    schema_entry = _get_collection_schema_entry(db, coll_name)
    facet_metadata = schema_entry.get("facet_metadata") or {}
    category_fields = list(facet_metadata.get("category_fields") or [])
    if not category_fields:
        category_fields = _infer_collection_category_fields(db, coll_name)
    return {
        "category_fields": category_fields,
        "numeric_fields": list(facet_metadata.get("numeric_fields") or []),
    }


def _get_collection_search_metadata(db, coll_name: str) -> dict:
    schema_entry = _get_collection_schema_entry(db, coll_name)
    search_metadata = schema_entry.get("search_metadata") or {}
    text_fields = list(search_metadata.get("text_fields") or [])
    if not text_fields:
        text_fields = _infer_collection_search_fields(db, coll_name, schema_entry=schema_entry)
    return {
        "text_fields": text_fields,
    }


def _infer_collection_category_fields(db, coll_name: str) -> List[str]:
    try:
        coll = db.collection(coll_name)
        indexes = coll.indexes()
    except Exception:
        return []

    # Older graphs only expose index intent indirectly. Restrict the fallback to
    # simple single-field hash indexes and suppress obviously noisy/internal fields.
    suppressed_fields = {
        "_key",
        "_id",
        "_rev",
        "_from",
        "_to",
        "id",
        "xref",
        "resolved_ids",
        "provenance",
        "creation",
        "updates",
    }
    inferred_fields = []
    for index in indexes:
        if index.get("type") != "hash":
            continue
        fields = index.get("fields") or []
        if len(fields) != 1:
            continue
        field = fields[0]
        if not field or field.startswith("_") or field in suppressed_fields:
            continue
        inferred_fields.append(field)
    return sorted(set(inferred_fields))


def _infer_collection_search_fields(db, coll_name: str, schema_entry: dict = None) -> List[str]:
    schema_entry = schema_entry or _get_collection_schema_entry(db, coll_name)
    schema_fields = set((schema_entry.get("fields") or {}).keys())
    preferred_fields = [
        "id",
        "name",
        "symbol",
        "preferred_symbol",
        "description",
        "label",
        "gene_name",
        "type",
        "xref",
        "uniprot_id",
        "ensembl_id",
        "refseq_id",
        "ncbi_id",
    ]
    return [field for field in preferred_fields if field in schema_fields]


def _parse_collection_facet_filters(request: Request, category_fields: List[str]) -> Dict[str, List[str]]:
    allowed_fields = set(category_fields)
    active_filters = {}
    for field in category_fields:
        key = f"facet_{field}"
        values = [value for value in request.query_params.getlist(key) if value != ""]
        if values:
            active_filters[field] = values
    for key in request.query_params.keys():
        if not key.startswith("facet_"):
            continue
        field = key[len("facet_"):]
        if field not in allowed_fields:
            continue
        values = [value for value in request.query_params.getlist(key) if value != ""]
        if values:
            active_filters[field] = values
    return active_filters


def _parse_collection_search_term(request: Request) -> str:
    return request.query_params.get("q", "").strip()


def _build_collection_query_params(page: int, page_size: int, facet_filters: Dict[str, List[str]],
                                   search_term: str = "", overrides: dict = None) -> List[tuple]:
    params = [("page", page), ("page_size", page_size)]
    if search_term:
        params.append(("q", search_term))
    merged_filters = {field: list(values) for field, values in facet_filters.items()}
    overrides = overrides or {}
    for field, values in overrides.items():
        merged_filters[field] = list(values)
    for field in sorted(merged_filters):
        for value in merged_filters[field]:
            params.append((f"facet_{field}", value))
    return params


def _build_collection_url(db_name: str, coll_name: str, page: int, page_size: int,
                          facet_filters: Dict[str, List[str]], search_term: str = "",
                          overrides: dict = None) -> str:
    params = _build_collection_query_params(page, page_size, facet_filters, search_term=search_term, overrides=overrides)
    query_string = urlencode(params, doseq=True)
    return f"{_app_path(f'/db/{db_name}/collection/{coll_name}')}?{query_string}"


def _build_collection_stats_url(db_name: str, coll_name: str, facet_filters: Dict[str, List[str]],
                                search_term: str = "") -> str:
    params = []
    if search_term:
        params.append(("q", search_term))
    for field in sorted(facet_filters):
        for value in facet_filters[field]:
            params.append((f"facet_{field}", value))
    query_string = urlencode(params, doseq=True)
    base_url = _app_path(f"/db/{db_name}/collection/{coll_name}/stats")
    return f"{base_url}?{query_string}" if query_string else base_url


def _build_collection_facets_url(db_name: str, coll_name: str, facet_filters: Dict[str, List[str]],
                                 search_term: str = "", page_size: int = 25) -> str:
    params = [("page_size", page_size)]
    if search_term:
        params.append(("q", search_term))
    for field in sorted(facet_filters):
        for value in facet_filters[field]:
            params.append((f"facet_{field}", value))
    query_string = urlencode(params, doseq=True)
    base_url = _app_path(f"/db/{db_name}/collection/{coll_name}/facets")
    return f"{base_url}?{query_string}" if query_string else base_url


def _build_collection_facet_url(db_name: str, coll_name: str, field: str,
                                facet_filters: Dict[str, List[str]],
                                search_term: str = "", page_size: int = 25) -> str:
    params = [("page_size", page_size)]
    if search_term:
        params.append(("q", search_term))
    for filter_field in sorted(facet_filters):
        for value in facet_filters[filter_field]:
            params.append((f"facet_{filter_field}", value))
    query_string = urlencode(params, doseq=True)
    base_url = _app_path(f"/db/{db_name}/collection/{coll_name}/facet/{url_quote(field, safe='')}")
    return f"{base_url}?{query_string}" if query_string else base_url


def _build_collection_download_url(db_name: str, coll_name: str, page: int, page_size: int,
                                   facet_filters: Dict[str, List[str]], search_term: str = "") -> str:
    params = _build_collection_query_params(page, page_size, facet_filters, search_term=search_term)
    query_string = urlencode(params, doseq=True)
    base_url = _app_path(f"/db/{db_name}/collection/{coll_name}/download.csv")
    return f"{base_url}?{query_string}" if query_string else base_url


def _get_search_constraint_clause(search_fields: List[str], search_term: str, bind_vars: dict,
                                  variable: str = "doc") -> str:
    if not search_term or not search_fields:
        return ""

    bind_vars["search_term"] = search_term.lower()
    field_bind_names = []
    for idx, field in enumerate(search_fields):
        bind_name = f"search_field_{idx}"
        bind_vars[bind_name] = field
        field_bind_names.append(bind_name)

    field_clauses = []
    for bind_name in field_bind_names:
        field_clauses.append(
            "("
            f"HAS({variable}, @{bind_name}) && {variable}[@{bind_name}] != null && "
            "("
            f"IS_ARRAY({variable}[@{bind_name}]) "
            f"? LENGTH(FOR item IN {variable}[@{bind_name}] FILTER CONTAINS(LOWER(TO_STRING(item)), @search_term) RETURN 1) > 0 "
            f": CONTAINS(LOWER(TO_STRING({variable}[@{bind_name}])), @search_term)"
            ")"
            ")"
        )
    return " OR ".join(field_clauses)


def _get_filter_constraint_clause(filter_settings: Dict[str, List[str]], bind_vars: dict, variable: str = "doc") -> str:
    clauses = []
    for idx, (field, values) in enumerate(filter_settings.items()):
        field_bind = f"facet_field_{idx}"
        values_bind = f"facet_values_{idx}"
        bind_vars[field_bind] = field
        bind_vars[values_bind] = [_coerce_facet_filter_value(value) for value in values]
        clauses.append(
            "("
            f"IS_ARRAY({variable}[@{field_bind}]) "
            f"? LENGTH(INTERSECTION({variable}[@{field_bind}], @{values_bind})) > 0 "
            f": {variable}[@{field_bind}] IN @{values_bind}"
            ")"
        )
    return " AND ".join(clauses)


def _get_facet_clause(field_bind_name: str, top_bind_name: str, variable: str = "doc") -> str:
    return f"""
        LET values = (
            !HAS({variable}, @{field_bind_name}) || {variable}[@{field_bind_name}] == null
                ? [null]
                : (IS_ARRAY({variable}[@{field_bind_name}]) ? UNIQUE({variable}[@{field_bind_name}]) : [{variable}[@{field_bind_name}]])
        )
            FOR item IN values
                COLLECT value = item WITH COUNT INTO count
                SORT count DESC, value
                LIMIT @{top_bind_name}
                RETURN {{ value, count }}"""


def _normalize_facet_value(value):
    return "null" if value is None else str(value)


def _format_facet_value(value):
    return "missing" if value is None else str(value)


def _coerce_facet_filter_value(value: str):
    lowered = value.lower()
    if lowered == "true":
        return True
    if lowered == "false":
        return False
    if lowered == "null":
        return None
    return value


def _build_collection_constraints(active_filters: Dict[str, List[str]], search_fields: List[str],
                                  search_term: str, bind_vars: dict, variable: str = "doc") -> str:
    clauses = []
    filter_clause = _get_filter_constraint_clause(active_filters, bind_vars, variable=variable)
    if filter_clause:
        clauses.append(f"({filter_clause})")
    search_clause = _get_search_constraint_clause(search_fields, search_term, bind_vars, variable=variable)
    if search_clause:
        clauses.append(f"({search_clause})")
    return " AND ".join(clauses)


def _build_collection_facet_panels(db, db_name: str, coll_name: str, page_size: int,
                                   category_fields: List[str], active_filters: Dict[str, List[str]],
                                   search_fields: List[str], search_term: str,
                                   top: int = 20) -> List[dict]:
    ordered_fields = sorted(
        category_fields,
        key=lambda field: (0 if active_filters.get(field) else 1, field),
    )
    return [
        _build_collection_facet_panel(
            db=db,
            db_name=db_name,
            coll_name=coll_name,
            page_size=page_size,
            field=field,
            active_filters=active_filters,
            search_fields=search_fields,
            search_term=search_term,
            top=top,
        )
        for field in ordered_fields
    ]


def _build_collection_facet_panel(db, db_name: str, coll_name: str, page_size: int,
                                  field: str, active_filters: Dict[str, List[str]],
                                  search_fields: List[str], search_term: str,
                                  top: int = 20) -> dict:
    other_filters = {k: v for k, v in active_filters.items() if k != field}
    bind_vars = {
        "facet_field": field,
        "facet_top": top,
    }
    filter_clause = _build_collection_constraints(
        active_filters=other_filters,
        search_fields=search_fields,
        search_term=search_term,
        bind_vars=bind_vars,
    )
    query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            {_get_facet_clause('facet_field', 'facet_top')}
    """
    rows = list(db.aql.execute(query, bind_vars=bind_vars))
    selected_values = set(active_filters.get(field, []))
    facet_values = []
    for row in rows:
        normalized_value = _normalize_facet_value(row.get("value"))
        updated_values = sorted(selected_values ^ {normalized_value})
        overrides = dict(active_filters)
        if updated_values:
            overrides[field] = updated_values
        else:
            overrides[field] = []
        facet_values.append({
            "value": normalized_value,
            "label": _format_facet_value(row.get("value")),
            "count": row.get("count", 0),
            "selected": normalized_value in selected_values,
            "href": _build_collection_url(
                db_name=db_name,
                coll_name=coll_name,
                page=1,
                page_size=page_size,
                facet_filters=active_filters,
                search_term=search_term,
                overrides=overrides,
            ),
        })
    return {
        "field": field,
        "title": field.replace("_", " "),
        "values": facet_values,
        "selected_values": sorted(selected_values),
        "clear_href": _build_collection_url(
            db_name=db_name,
            coll_name=coll_name,
            page=1,
            page_size=page_size,
            facet_filters=active_filters,
            search_term=search_term,
            overrides={field: []},
        ),
    }


def _build_active_filter_summary(db_name: str, coll_name: str, page_size: int,
                                 active_filters: Dict[str, List[str]], search_term: str = "") -> List[dict]:
    summaries = []
    for field in sorted(active_filters):
        selected_values = list(active_filters[field])
        values = []
        for value in selected_values:
            remaining_values = [item for item in selected_values if item != value]
            overrides = dict(active_filters)
            if remaining_values:
                overrides[field] = remaining_values
            else:
                overrides[field] = []
            values.append({
                "label": _format_facet_value(None if value == "null" else value),
                "remove_href": _build_collection_url(
                    db_name=db_name,
                    coll_name=coll_name,
                    page=1,
                    page_size=page_size,
                    facet_filters=active_filters,
                    search_term=search_term,
                    overrides=overrides,
                ),
            })
        summaries.append({
            "field": field,
            "values": values,
            "clear_href": _build_collection_url(
                db_name=db_name,
                coll_name=coll_name,
                page=1,
                page_size=page_size,
                facet_filters=active_filters,
                search_term=search_term,
                overrides={field: []},
            ),
        })
    return summaries


def _get_collection_preview_fields(is_edge: bool, facet_metadata: dict, search_metadata: dict) -> List[str]:
    """Fields to show in collection tables without loading every large document field."""
    fields = ["_key"]
    if is_edge:
        fields.extend(["_from", "_to"])
    priority = [
        "id",
        "name",
        "symbol",
        "preferred_symbol",
        "type",
        "description",
        "label",
        "gene_name",
        "uniprot_id",
        "ensembl_id",
        "refseq_id",
        "ncbi_id",
    ]
    fields.extend(priority)
    fields.extend(sorted(facet_metadata.get("category_fields") or []))
    fields.extend(sorted(facet_metadata.get("numeric_fields") or []))
    fields.extend(search_metadata.get("text_fields") or [])
    return list(dict.fromkeys(fields))


def _build_collection_facet_loaders(db_name: str, coll_name: str, page_size: int,
                                    category_fields: List[str],
                                    active_filters: Dict[str, List[str]],
                                    search_term: str = "") -> List[dict]:
    ordered_fields = sorted(
        category_fields,
        key=lambda field: (0 if active_filters.get(field) else 1, field),
    )
    return [
        {
            "field": field,
            "title": field.replace("_", " "),
            "url": _build_collection_facet_url(
                db_name=db_name,
                coll_name=coll_name,
                field=field,
                facet_filters=active_filters,
                search_term=search_term,
                page_size=page_size,
            ),
        }
        for field in ordered_fields
    ]


def _get_dashboard_collection_summaries(db) -> List[dict]:
    collections = []
    for coll in db.collections():
        if _is_hidden_qa_browser_collection(coll["name"]):
            continue
        coll_type = coll.get("type")
        is_edge = coll_type in ("edge", 3)
        count_cursor = db.aql.execute(
            f"FOR doc IN `{coll['name']}` COLLECT WITH COUNT INTO c RETURN c"
        )
        count = list(count_cursor)[0]
        collections.append({
            "name": coll["name"],
            "type": "edge" if is_edge else "document",
            "count": count,
            "qa_visible": _is_qa_visible_edge_collection(coll["name"]),
        })
    collections.sort(key=lambda c: c["name"])
    return collections


def _get_dashboard_collection_shell(db) -> List[dict]:
    collections = []
    for coll in db.collections():
        if _is_hidden_qa_browser_collection(coll["name"]):
            continue
        coll_type = coll.get("type")
        is_edge = coll_type in ("edge", 3)
        collections.append({
            "name": coll["name"],
            "type": "edge" if is_edge else "document",
            "count": None,
            "qa_visible": _is_qa_visible_edge_collection(coll["name"]),
        })
    collections.sort(key=lambda c: c["name"])
    return collections


def _build_browser_home_context(request: Request) -> dict:
    # Arango databases
    arango_databases = []
    arango_url = _credentials.get("url", "")
    if _credentials:
        try:
            sys_db = get_sys_db()
            arango_databases = [db for db in sys_db.databases() if not db.startswith("_")]
        except Exception:
            pass

    # MySQL databases
    mysql_sources = []
    mysql_databases = []
    system_dbs = {"information_schema", "mysql", "performance_schema", "sys"}
    for source_id, source in _mysql_sources.items():
        credentials = source["credentials"]
        engine = get_mysql_engine(source_id)
        if not engine:
            continue
        host = credentials.get("url", "localhost")
        port = credentials.get("port", 3306)
        mysql_url = f"{host}:{port}"
        try:
            with engine.connect() as conn:
                result = conn.execute(text("SHOW DATABASES"))
                dbs = [row[0] for row in result if row[0] not in system_dbs]
                print(f"MySQL databases for {source_id}: {dbs}")
                mysql_databases.extend(dbs)
                mysql_sources.append({
                    "id": source_id,
                    "label": source["label"],
                    "url": mysql_url,
                    "databases": dbs,
                })
        except Exception:
            print(f"Failed to connect to MySQL source {source_id}:", sys.exc_info()[1])
            pass

    return {
        "request": request,
        "databases": arango_databases,
        "arango_url": arango_url,
        "mysql_databases": mysql_databases,
        "mysql_sources": mysql_sources,
        "demo_queries_enabled": _demo_queries_enabled,
    }


def _storage_from_credentials(credentials_config: dict, *, use_internal_url: bool):
    from src.infrastructure.object_storage import AwsAssumeRoleCredentials
    from src.infrastructure.object_storage import AwsAssumeRoleStorage
    from src.infrastructure.object_storage import S3CompatibleStorage
    from src.shared.db_credentials import DBCredentials

    if "role_arn" in credentials_config or credentials_config.get("type") == "aws_assume_role":
        return AwsAssumeRoleStorage(AwsAssumeRoleCredentials.from_yaml(credentials_config))
    return S3CompatibleStorage(DBCredentials.from_yaml(credentials_config), use_internal_url=use_internal_url)


@app.get("/", response_class=HTMLResponse)
async def home(request: Request):
    return templates.TemplateResponse(request, "home.html", {
        "request": request,
    })


@app.get("/qa-browser", response_class=HTMLResponse)
async def qa_browser_home(request: Request):
    return templates.TemplateResponse(
        request,
        "qa_browser_home.html",
        _build_browser_home_context(request),
    )


@app.get("/ramp-id-qa", response_class=HTMLResponse)
def ramp_id_qa(
    request: Request,
    id: str = "",
    ids: str = "",
    stages: str = "",
    denylist_rank: int = 0,
    denylist_pair: str = "",
    mw_finding: str = "",
):
    query_id = (id or ids or "").strip()
    pair = None
    pair_error = None
    if denylist_pair.strip():
        try:
            pair = _normalize_harmonization_denylist_pair(denylist_pair)
            query_ids = _parse_metabolite_identifier_query(query_id)
            query_id = " ".join(dict.fromkeys([*query_ids, *pair]))
        except ValueError as exc:
            pair_error = str(exc)
    selected_stage_keys = _parse_metabolite_snapshot_key_filter(stages)
    if query_id and not selected_stage_keys:
        selected_stage_keys = (
            _default_metabolite_denylist_pair_stage_filter()
            if pair
            else _default_metabolite_snapshot_key_filter()
        )
    result = None
    error = None
    overview = None
    denylist_review = None
    denylist_pair_review = None
    mw_validation_review = None
    if pair_error:
        error = pair_error
    elif query_id:
        if mw_finding:
            try:
                if len(selected_stage_keys) != 1:
                    raise ValueError("An MW finding link requires exactly one saved stage.")
                mw_validation_review = _build_mw_validation_review(
                    selected_stage_keys[0], mw_finding
                )
            except Exception as exc:
                error = str(exc)
        try:
            result = _load_metabolite_identifier_qa_many(
                query_id,
                selected_stage_keys,
                include_identifier_details=False,
            )
            if denylist_rank and len(selected_stage_keys) == 1:
                stage = _get_harmonization_stage(selected_stage_keys[0])
                denylist_review = _harmonization_denylist_review_from_stage(stage, denylist_rank)
            if pair:
                denylist_pair_review = _build_harmonization_denylist_pair_review(
                    pair,
                    result,
                    _load_metabolite_edge_removal_curations(),
                )
        except Exception as exc:
            graph_error = str(exc)
            error = (
                f"{error} Graph review also failed: {graph_error}"
                if error else graph_error
            )
    else:
        overview = _load_harmonization_pipeline_workbench()
    return templates.TemplateResponse(request, "ramp_id_qa.html", {
        "request": request,
        "query_id": query_id,
        "selected_snapshot_keys": selected_stage_keys,
        "selected_stage_keys": selected_stage_keys,
        "result": result,
        "denylist_review": denylist_review,
        "denylist_pair_review": denylist_pair_review,
        "denylist_pair": denylist_pair,
        "mw_finding": mw_finding,
        "mw_validation_review": mw_validation_review,
        "error": error,
        "overview": overview,
    })


@app.get("/ramp-id-qa/curations", response_class=HTMLResponse)
def ramp_id_qa_curations(request: Request):
    snapshots = {}
    stream_errors = {}
    error = None
    try:
        storage = _curation_cart_storage()
        with ThreadPoolExecutor(max_workers=len(_QA_CURATION_TYPES)) as executor:
            futures = {
                curation_type: executor.submit(
                    resolve_curation_type,
                    storage,
                    curation_type,
                    allow_missing=True,
                )
                for curation_type in _QA_CURATION_TYPES
            }
            for curation_type, future in futures.items():
                try:
                    snapshots[curation_type] = future.result()
                except Exception as exc:
                    stream_errors[curation_type] = str(exc)
        query_params = request.query_params
        try:
            page = max(1, int(query_params.get("page") or 1))
        except (TypeError, ValueError):
            page = 1
        index = build_active_curation_review(
            snapshots,
            ReviewFilters(
                query=str(query_params.get("q") or "").strip(),
                curation_types=tuple(query_params.getlist("curation_type")),
                sources=tuple(value.upper() for value in query_params.getlist("source")),
                batches=tuple(query_params.getlist("batch")),
                curators=tuple(query_params.getlist("curator")),
                page=page,
            ),
            stream_errors=stream_errors,
        )
        index["loaded_at"] = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    except Exception as exc:
        error = str(exc)
        index = build_active_curation_review({}, ReviewFilters())
        index["loaded_at"] = datetime.now(timezone.utc).isoformat().replace("+00:00", "Z")
    return templates.TemplateResponse(request, "ramp_id_curations.html", {
        "request": request,
        "index": index,
        "error": error,
    })


@app.get("/ramp-id-qa/curations/mw-review")
def ramp_id_qa_curation_mw_review(
    finding_id: str,
    anchor_id: str,
    id: str = "",
):
    member_ids = _parse_metabolite_identifier_query(id)
    candidate = _find_mw_review_stage(finding_id, anchor_id, member_ids)
    if candidate:
        query = urlencode({
            "id": " ".join(member_ids),
            "stages": candidate["stage_key"],
            "mw_finding": candidate["finding_id"],
        })
    else:
        query = urlencode({"id": " ".join(member_ids or [anchor_id])})
    return _redirect_to(f"/ramp-id-qa?{query}")


@app.post("/ramp-id-qa/pipelines")
async def ramp_id_qa_save_pipeline(
    request: Request,
    pipeline_key: str = Form(""),
    pipeline_name: str = Form(""),
    enabled_rules: List[str] = Form([]),
    run_now: str = Form(""),
):
    form = await request.form()
    rule_parameters = {}
    for key, value in form.multi_items():
        if not key.startswith("rule_param__"):
            continue
        _prefix, rule_id, parameter_id = key.split("__", 2)
        rule_parameters.setdefault(rule_id, {})[parameter_id] = value
    pipeline = _upsert_harmonization_pipeline(pipeline_key, pipeline_name, enabled_rules, rule_parameters)
    if run_now:
        _enqueue_metabolite_snapshot_job(
            "run_pipeline",
            f"Run {pipeline.get('name')}",
            {"pipeline_key": pipeline["_key"]},
        )
    return _redirect_to("/ramp-id-qa")


@app.get("/ramp-id-qa/pipelines/{pipeline_key}", response_class=HTMLResponse)
def ramp_id_qa_pipeline_page(request: Request, pipeline_key: str):
    pipeline = None
    overview = {"pipeline_stage_stats": []}
    error = None
    try:
        page = _load_harmonization_pipeline_page(pipeline_key)
        pipeline = page["pipeline"]
        overview = page["overview"]
    except Exception as exc:
        error = str(exc)
    return templates.TemplateResponse(request, "ramp_id_pipeline.html", {
        "request": request,
        "pipeline_key": pipeline_key,
        "pipeline": pipeline,
        "overview": overview,
        "error": error,
    })


@app.post("/ramp-id-qa/pipelines/{pipeline_key}/run")
def ramp_id_qa_run_pipeline(pipeline_key: str):
    _enqueue_metabolite_snapshot_job(
        "run_pipeline",
        f"Run {pipeline_key}",
        {"pipeline_key": pipeline_key},
    )
    return _redirect_to("/ramp-id-qa")


@app.post("/ramp-id-qa/pipelines/{pipeline_key}/rename")
def ramp_id_qa_rename_pipeline(pipeline_key: str, pipeline_name: str = Form("")):
    _rename_harmonization_pipeline(pipeline_key, pipeline_name)
    return _redirect_to("/ramp-id-qa")


@app.post("/ramp-id-qa/pipelines/{pipeline_key}/delete")
def ramp_id_qa_delete_pipeline(pipeline_key: str):
    try:
        pipeline = _get_harmonization_pipeline(get_db("metabolite_harmonization"), pipeline_key)
        pipeline_label = pipeline.get("name") or pipeline_key
    except Exception:
        pipeline_label = pipeline_key
    _enqueue_metabolite_snapshot_job(
        "delete_pipeline",
        f"Delete {pipeline_label}",
        {"pipeline_key": pipeline_key},
    )
    return _redirect_to("/ramp-id-qa")


@app.get("/ramp-id-qa/stage-comparison", response_class=HTMLResponse)
def ramp_id_qa_stage_comparison(
    request: Request,
    left_stage: str = "",
    right_stage: str = "",
    limit: int = 100,
):
    comparison = None
    error = None
    normalized_limit = max(10, min(limit, 250))
    try:
        comparison = _load_metabolite_snapshot_comparison(left_stage, right_stage, normalized_limit)
    except Exception as exc:
        error = str(exc)
    if comparison is None:
        comparison = {
            "left_snapshot": None,
            "right_snapshot": None,
            "snapshots": _list_harmonization_stages(),
            "components": [],
            "display_limit": normalized_limit,
            "comparison_kind": "stage",
        }
    return templates.TemplateResponse(request, "ramp_id_snapshot_compare.html", {
        "request": request,
        "left_snapshot_key": left_stage,
        "right_snapshot_key": right_stage,
        "limit": normalized_limit,
        "comparison": comparison,
        "error": error,
    })


@app.get("/ramp-id-qa/api/stage-comparison-visualization")
def ramp_id_qa_stage_comparison_visualization(
    ids: str = "",
    left_stage: str = "",
    right_stage: str = "",
):
    return _load_metabolite_snapshot_comparison_visualization(
        _parse_metabolite_identifier_query(ids),
        left_stage,
        right_stage,
    )


@app.get("/ramp-id-qa/stages/{stage_key}", response_class=HTMLResponse)
def ramp_id_qa_stage_stats(request: Request, stage_key: str):
    stats = None
    error = None
    try:
        stats = _load_harmonization_stage_stats(stage_key)
    except Exception as exc:
        error = str(exc)
    return templates.TemplateResponse(request, "ramp_id_stage_stats.html", {
        "request": request,
        "stage_key": stage_key,
        "stats": stats,
        "error": error,
    })


@app.get("/ramp-id-qa/harmonization-jobs", response_class=HTMLResponse)
def ramp_id_qa_harmonization_jobs(request: Request):
    overview = _load_harmonization_pipeline_workbench()
    return templates.TemplateResponse(request, "ramp_id_pipeline_table.html", {
        "request": request,
        "overview": overview,
    })


@app.get("/ramp-id-qa/api/metabolite")
def ramp_id_qa_metabolite(id: str = "", ids: str = "", stages: str = ""):
    query_id = (id or ids or "").strip()
    query_ids = _parse_metabolite_identifier_query(query_id)
    selected_snapshot_keys = _parse_metabolite_snapshot_key_filter(stages)
    if query_ids and not selected_snapshot_keys:
        selected_snapshot_keys = _default_metabolite_snapshot_key_filter()
    if len(query_ids) == 1:
        return _load_metabolite_identifier_qa(query_ids[0], selected_snapshot_keys)
    return _load_metabolite_identifier_qa_many(query_id, selected_snapshot_keys)


@app.get("/api/curation-cart")
@app.get("/ramp-id-qa/api/curation-cart")
def ramp_id_qa_curation_cart(request: Request, curator: str = "", curator_name: str = ""):
    curator_id, display_name = _curator_identity(request, {
        "curator": curator,
        "curator_name": curator_name,
    })
    try:
        return _combined_curation_cart(
            _curation_cart_storage(), curator_id, display_name
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/ramp-id-qa/api/curation-cart/items")
async def ramp_id_qa_add_curation_cart_item(request: Request):
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    curator_id, curator_name = _curator_identity(request, payload)
    action = str(payload.get("action") or "remove_edge")
    if action == "assert_same_clique":
        curation_type = METABOLITE_EXPECTED_CLIQUES
        member_ids = _normalize_expected_clique_member_ids(payload.get("member_ids") or "")
        operation = _metabolite_expected_clique_assertion_operation(
            member_ids,
            str(payload.get("name") or ""),
            str(payload.get("rationale") or ""),
            await run_in_threadpool(_missing_metabolite_identifier_ids, member_ids),
        )
    elif action == "retire_assertion":
        curation_type = METABOLITE_EXPECTED_CLIQUES
        operation = _metabolite_retire_assertion_operation(
            str(payload.get("assertion_id") or ""),
            str(payload.get("note") or ""),
        )
    elif action in {"accept_mw_discrepancy", "reopen_mw_discrepancy"}:
        curation_type = METABOLITE_MW_ADJUDICATIONS
        operation = _metabolite_mw_adjudication_operation(
            action,
            str(payload.get("stage_key") or ""),
            str(payload.get("finding_id") or ""),
            str(payload.get("reason") or ""),
            payload.get("supporting_ids") or "",
            str(payload.get("note") or ""),
        )
    elif action == "set_properties":
        curation_type = METABOLITE_RECORD_PROPERTIES
        values = payload.get("values", {})
        remove_overrides = payload.get("remove_overrides", [])
        if not isinstance(values, dict) or not isinstance(remove_overrides, list):
            raise HTTPException(
                status_code=400,
                detail="Property values must be an object and remove_overrides must be a list.",
            )
        operation = _metabolite_property_operation(
            str(payload.get("target_id") or ""),
            values,
            remove_overrides,
            str(payload.get("note") or ""),
        )
    elif action in {"suppress_record", "restore_record"}:
        curation_type = METABOLITE_RECORD_SUPPRESSIONS
        operation = _metabolite_record_suppression_operation(
            action,
            str(payload.get("target_id") or ""),
            str(payload.get("note") or ""),
        )
        published_state = await run_in_threadpool(_load_metabolite_edge_removal_curations)
        published_action = (published_state.get("record_states") or {}).get(
            operation["target"]["id"]
        )
        if action == "restore_record" and published_action != "suppress_record":
            raise HTTPException(
                status_code=409,
                detail="This identifier is no longer actively suppressed. Reload the record before restoring it.",
            )
    else:
        curation_type = METABOLITE_EQUIVALENCE_EDGES
        operation = _metabolite_edge_decision_operation(
            action,
            str(payload.get("start_id") or ""),
            str(payload.get("end_id") or ""),
            str(payload.get("note") or ""),
        )
        if operation["action"] == "retain_edge":
            pair = tuple(sorted((operation["start_id"], operation["end_id"])))
            published_state = await run_in_threadpool(_load_metabolite_edge_removal_curations)
            if (published_state.get("pair_states") or {}).get(pair) != "remove_edge":
                raise HTTPException(
                    status_code=409,
                    detail="This pair is no longer an active published removal. Reload the review before retaining it.",
                )
    try:
        await run_in_threadpool(
            add_cart_operation,
            _curation_cart_storage(),
            curation_type,
            curator_id,
            curator_name,
            operation,
            bool(payload.get("replace_target", False)),
        )
        return await run_in_threadpool(
            _combined_curation_cart,
            _curation_cart_storage(),
            curator_id,
            curator_name,
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.get("/ramp-id-qa/api/expected-clique-assertions")
def ramp_id_qa_expected_clique_assertions(stages: str = ""):
    stage_keys = _parse_metabolite_snapshot_key_filter(stages)
    curation_state = _load_metabolite_edge_removal_curations()
    assertions = curation_state.get("assertions") or []
    evaluated = _load_expected_clique_assertion_results(
        get_db("metabolite_harmonization"),
        assertions,
        stage_keys,
    )
    return {
        "assertion_fingerprint": curation_state.get("assertion_fingerprint"),
        "assertion_batch_ids": curation_state.get("assertion_batch_ids") or [],
        "assertions": assertions,
        "stages": stage_keys,
        **evaluated,
    }


@app.get("/ramp-id-qa/api/stages/{stage_key}/curation-cart-flags")
def ramp_id_qa_stage_curation_cart_flags(
    stage_key: str,
    request: Request,
    curator: str = "",
    curator_name: str = "",
):
    curator_id, display_name = _curator_identity(request, {
        "curator": curator,
        "curator_name": curator_name,
    })
    try:
        cart = _combined_curation_cart(
            _curation_cart_storage(), curator_id, display_name
        )
        return _load_harmonization_stage_cart_flags(stage_key, cart.get("operations") or [])
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.delete("/api/curation-cart/items/{operation_id}")
@app.delete("/ramp-id-qa/api/curation-cart/items/{operation_id}")
async def ramp_id_qa_remove_curation_cart_item(operation_id: str, request: Request):
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    curator_id, curator_name = _curator_identity(request, payload)
    try:
        storage = _curation_cart_storage()
        requested_type = str(payload.get("curation_type") or "").strip()
        candidate_types = [requested_type] if requested_type in _QA_CURATION_TYPES else list(_QA_CURATION_TYPES)
        for curation_type in candidate_types:
            cart = await run_in_threadpool(load_cart, storage, curation_type, curator_id, curator_name)
            if any(item.get("operation_id") == operation_id for item in cart.get("operations") or []):
                await run_in_threadpool(
                    remove_cart_operation,
                    storage,
                    curation_type,
                    curator_id,
                    curator_name,
                    operation_id,
                )
                break
        return await run_in_threadpool(
            _combined_curation_cart, storage, curator_id, curator_name
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/curation-cart/publish")
@app.post("/ramp-id-qa/api/curation-cart/publish")
async def ramp_id_qa_publish_curation_cart(request: Request):
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    curator_id, curator_name = _curator_identity(request, payload)
    try:
        storage = _curation_cart_storage()
        batch_name = str(payload.get("batch_name") or "").strip()
        if not batch_name:
            raise ValueError("batch name is required")
        carts = {}
        for curation_type in _QA_CURATION_TYPES:
            cart = await run_in_threadpool(
                load_cart, storage, curation_type, curator_id, curator_name
            )
            if cart.get("operations"):
                for operation in cart["operations"]:
                    validate_operation(curation_type, operation)
                if curation_type in {CHEBI_RECORD_PROPERTIES, METABOLITE_RECORD_PROPERTIES}:
                    curation_sets = {
                        (operation.get("target") or {}).get("curation_set")
                        for operation in cart["operations"]
                    }
                    curation_sets.discard(None)
                    if len(curation_sets) != 1:
                        raise ValueError(
                            "Record-property curations must target exactly one graph per batch"
                        )
                carts[curation_type] = cart
        if not carts:
            raise ValueError("The curation review is empty")
        published = []
        for curation_type in carts:
            try:
                result = await run_in_threadpool(
                    publish_cart,
                    storage,
                    curation_type,
                    curator_id,
                    curator_name,
                    batch_name,
                    str(payload.get("description") or ""),
                )
                _invalidate_harmonization_curation_cache()
            except Exception as exc:
                if not published:
                    raise
                remaining = await run_in_threadpool(
                    _combined_curation_cart,
                    storage,
                    curator_id,
                    curator_name,
                )
                remaining_types = sorted({
                    operation.get("curation_type")
                    for operation in remaining.get("operations") or []
                    if operation.get("curation_type")
                })
                return JSONResponse(status_code=207, content={
                    "partial": True,
                    "published": published,
                    "failed_type": curation_type,
                    "error": str(exc),
                    "remaining_cart": remaining,
                    "remaining_types": remaining_types,
                    "operation_count": sum(item["operation_count"] for item in published),
                    "batch_count": len(published),
                })
            published.append({"curation_type": curation_type, **result})
        return {
            "published": published,
            "operation_count": sum(item["operation_count"] for item in published),
            "batch_count": len(published),
        }
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc
    except HTTPException:
        raise
    except Exception as exc:
        raise HTTPException(status_code=502, detail=f"Curation publication failed: {exc}") from exc


def _load_cure_entity_resolver():
    graph_dir = _disease_graph_dir or ""
    target_graph_dir = _target_graph_dir or ""
    cache_key = (graph_dir, target_graph_dir)
    with _cure_entity_resolver_lock:
        if (
            _cure_entity_resolver_cache.get("index") is not None
            and _cure_entity_resolver_cache.get("graph_dir") == cache_key
        ):
            return _cure_entity_resolver_cache["index"]
        index = build_cure_entity_resolver_index(graph_dir, target_graph_dir)
        _cure_entity_resolver_cache.update({
            "graph_dir": cache_key,
            "index": index,
            "loaded_at": datetime.now(timezone.utc).isoformat(),
        })
        return index


def _parse_cure_resolver_queries(value) -> list[str]:
    if isinstance(value, list):
        raw_values = value
    else:
        raw_values = re.split(r"[\r\n]+", str(value or ""))
    queries = []
    for raw in raw_values:
        text = str(raw or "").strip()
        if text:
            queries.append(text)
    return queries


@app.get("/cure-entity-resolver", response_class=HTMLResponse)
def cure_entity_resolver(request: Request):
    return templates.TemplateResponse(request, "cure_entity_resolver.html", {
        "request": request,
    })


@app.get("/cure-entity-resolver/api/stats")
def cure_entity_resolver_stats():
    index = _load_cure_entity_resolver()
    stats = index.stats()
    stats["loaded_at"] = _cure_entity_resolver_cache.get("loaded_at")
    stats["disease_graph_dir"] = _disease_graph_dir
    stats["target_graph_dir"] = _target_graph_dir
    return stats


@app.post("/cure-entity-resolver/api/resolve")
async def cure_entity_resolver_resolve(request: Request):
    payload = await request.json()
    queries = _parse_cure_resolver_queries(
        payload.get("queries") if "queries" in payload else payload.get("text", "")
    )
    if not queries:
        raise HTTPException(status_code=400, detail="Provide at least one query.")
    if len(queries) > 1000:
        raise HTTPException(status_code=400, detail="At most 1000 queries are allowed per request.")
    entity_type = str(payload.get("entity_type") or "auto").strip().lower()
    allowed_entity_types = {
        "auto",
        "all",
        "clinical",
        "drug",
        "disease",
        "phenotype",
        "adverse_event",
        "gene",
        "protein",
        "sequence_variant",
    }
    if entity_type not in allowed_entity_types:
        allowed = ", ".join(sorted(allowed_entity_types))
        raise HTTPException(status_code=400, detail=f"entity_type must be one of: {allowed}.")
    try:
        top_k = int(payload.get("top_k", 5))
    except Exception:
        top_k = 5
    top_k = max(1, min(top_k, 10))

    index = _load_cure_entity_resolver()
    results = await run_in_threadpool(
        index.resolve_many,
        queries,
        entity_type=entity_type,
        top_k=top_k,
    )
    return {
        "query_count": len(queries),
        "entity_type": entity_type,
        "top_k": top_k,
        "index": index.stats(),
        "results": results,
    }


# ---------------------------------------------------------------------------
# Target Harmonizer Explorer
# ---------------------------------------------------------------------------


_target_qc_dir: str = ""


def _load_target_graph():
    if not _target_graph_dir:
        raise HTTPException(status_code=500, detail="No --target-graph-dir configured.")
    return load_target_graph_data(_target_graph_dir, qc_dir=_target_qc_dir)


def _target_manifest_snapshot() -> dict[str, Any]:
    if not _target_graph_dir:
        return {}
    manifest_path = Path(_target_graph_dir) / "manifest.json"
    if not manifest_path.exists():
        return {}
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    return manifest if isinstance(manifest, dict) else {}


def _target_version_display(manifest: dict[str, Any]) -> str:
    raw = (
        manifest.get("target_release_version")
        or manifest.get("release_version")
        or manifest.get("pipeline_version")
        or "—"
    )
    text = str(raw).strip()
    if not text or text == "—":
        return "—"
    if text.startswith("v."):
        return text
    if text.startswith("v"):
        return f"v.{text[1:]}"
    return f"v.{text}"


def _target_updated_last(manifest: dict[str, Any]) -> str:
    raw = manifest.get("updated_last") or manifest.get("generated_at") or "—"
    text = str(raw).strip()
    return text[:10] if text and text != "—" else "—"


def _load_variant_graph():
    if not _variant_graph_dir:
        raise HTTPException(status_code=500, detail="No --variant-graph-dir configured.")
    return load_variant_graph_data(_variant_graph_dir)


def _variant_manifest_snapshot() -> dict[str, Any]:
    if not _variant_graph_dir:
        return {}
    manifest_path = Path(_variant_graph_dir) / "manifest.json"
    if not manifest_path.exists():
        return {}
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    return manifest if isinstance(manifest, dict) else {}


def _variant_version_display(manifest: dict[str, Any]) -> str:
    raw = manifest.get("version") or manifest.get("variant_release_version") or "—"
    text = str(raw).strip()
    if not text or text == "—":
        return "—"
    if text.startswith("v."):
        return text
    if text.startswith("v"):
        return f"v.{text[1:]}"
    return f"v.{text}"


def _variant_updated_last(manifest: dict[str, Any]) -> str:
    raw = manifest.get("updated_last") or manifest.get("createdAt") or manifest.get("generated_at") or "—"
    text = str(raw).strip()
    return text[:10] if text and text != "—" else "—"


def _load_drug_graph():
    if not _drug_graph_dir:
        raise HTTPException(status_code=500, detail="No --drug-graph-dir configured.")
    return load_drug_graph_data(
        _drug_graph_dir,
        resolver_index_path=_drug_resolver_index or None,
        require_resolver_index=bool(_drug_resolver_index),
        verify_resolver_checksum=_verify_drug_resolver_checksum,
    )


def _drug_manifest_snapshot() -> dict[str, Any]:
    if not _drug_graph_dir:
        return {}
    manifest_path = Path(_drug_graph_dir) / "manifest.json"
    if not manifest_path.exists():
        return {}
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    return manifest if isinstance(manifest, dict) else {}


def _drug_version_display(manifest: dict[str, Any]) -> str:
    raw = manifest.get("version") or manifest.get("drug_release_version") or "—"
    text = str(raw).strip()
    if not text or text == "—":
        return "—"
    if text.startswith("v."):
        return text
    if text.startswith("v"):
        return f"v.{text[1:]}"
    return f"v.{text}"


def _drug_updated_last(manifest: dict[str, Any]) -> str:
    raw = manifest.get("updated_last") or manifest.get("createdAt") or manifest.get("generated_at") or "—"
    text = str(raw).strip()
    return text[:10] if text and text != "—" else "—"


@app.get("/target-id-qa", response_class=HTMLResponse)
def target_id_qa(request: Request, ids: str = "", tab: str = ""):
    selected_ids = " | ".join([part for part in re.split(r"[\s,|]+", ids or "") if part])
    default_tab = "graph" if selected_ids else "dashboard"
    manifest = _target_manifest_snapshot()
    return templates.TemplateResponse(request, "target_id_qa.html", {
        "request": request,
        "selected_ids": selected_ids,
        "active_tab": tab or default_tab,
        "manifest": manifest,
        "target_version_display": _target_version_display(manifest),
        "target_updated_last": _target_updated_last(manifest),
    })


@app.get("/target-id-qa/api/stats")
def target_id_qa_stats():
    data = _load_target_graph()
    return compute_target_stats(data)


@app.get("/target-id-qa/api/search")
def target_id_qa_search(
    q: str = "",
    target_type: str = "",
    namespace: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_target_graph()
    return search_targets(data, q=q, target_type=target_type, namespace=namespace, page=page, per_page=per_page)


@app.get("/target-id-qa/api/resolve-quick")
def target_id_qa_resolve_quick(q: str = "", limit: int = 8, target_type: str = ""):
    if not q.strip():
        raise HTTPException(status_code=400, detail="q parameter required.")
    queries = [part.strip() for part in re.split(r"[\n|]+", q) if part.strip()]
    data = _load_target_graph()
    if len(queries) == 1:
        return resolve_target_query(data, queries[0], limit=max(1, min(limit, 50)), target_type=target_type)
    return bulk_resolve_targets(data, queries[:500], limit=max(1, min(limit, 50)), target_type=target_type)


@app.post("/target-id-qa/api/resolve")
async def target_id_qa_resolve(request: Request, limit: int = 8, target_type: str = ""):
    payload = await request.json()
    if isinstance(payload, list):
        queries = payload
    elif isinstance(payload, dict):
        queries = payload.get("queries") or payload.get("targets") or payload.get("names") or []
        limit = payload.get("limit", limit)
        target_type = payload.get("target_type", target_type)
    else:
        queries = []
    if isinstance(queries, str):
        queries = [line.strip() for line in queries.splitlines() if line.strip()]
    if not isinstance(queries, list) or not queries:
        raise HTTPException(status_code=400, detail="Request body must include a queries list.")
    data = _load_target_graph()
    return bulk_resolve_targets(
        data,
        [str(q).strip() for q in queries if str(q).strip()][:500],
        limit=max(1, min(int(limit), 50)),
        target_type=str(target_type or ""),
    )


@app.get("/target-id-qa/api/graph")
def target_id_qa_graph(ids: str = ""):
    data = _load_target_graph()
    return build_target_graph_payload(data, ids)


@app.get("/target-id-qa/download-filtered")
def target_id_qa_download_filtered(
    q: str = "",
    target_type: str = "",
    namespace: str = "",
    columns: str = "",
    format: str = "tsv",
):
    data = _load_target_graph()
    fmt = "csv" if format == "csv" else "tsv"
    selected_columns = [col.strip() for col in columns.split(",") if col.strip()]
    content = export_targets(data, q=q, target_type=target_type, namespace=namespace, fmt=fmt, columns=selected_columns)
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    release = data.manifest.get("target_release_version") or data.manifest.get("release_version") or "current"
    release = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(release)).strip("_") or "current"
    filename = f"ODIN_{release}_targets.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/target-id-qa/download-divergences")
def target_id_qa_download_divergences(
    entity_type: str = "",
    triage_category: str = "",
    review_group: str = "",
    status: str = "",
    format: str = "tsv",
):
    data = _load_target_graph()
    if not data.review_can_write:
        raise HTTPException(status_code=403, detail="Target QC export is disabled in public/read-only mode.")
    fmt = "csv" if format == "csv" else "tsv"
    content = export_divergences(
        data, entity_type=entity_type, triage_category=triage_category,
        review_group=review_group, status=status, fmt=fmt,
    )
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    release = data.manifest.get("target_release_version") or data.manifest.get("release_version") or "current"
    release = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(release)).strip("_") or "current"
    filename = f"ODIN_{release}_divergences.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/target-id-qa/api/entity-ids/{entity_type}")
def target_id_qa_entity_ids(entity_type: str, format: str = "tsv"):
    """Download gene/protein/transcript IDs in pipeline-native column format.

    This endpoint provides the harmonized ID files in the format expected by
    the ``target_graph`` ArangoDB build pipeline (Keith's ingestion).

    Examples::

        GET /target-id-qa/api/entity-ids/gene
        GET /target-id-qa/api/entity-ids/protein?format=csv
        GET /target-id-qa/api/entity-ids/transcript
    """
    entity_type = entity_type.lower().strip()
    if entity_type not in ("gene", "protein", "transcript"):
        raise HTTPException(status_code=400, detail=f"Unknown entity type: {entity_type!r}. Use gene, protein, or transcript.")
    data = _load_target_graph()
    fmt = "csv" if format == "csv" else "tsv"
    content = export_entity_ids(data, entity_type=entity_type, fmt=fmt)
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    release = data.manifest.get("target_release_version") or data.manifest.get("release_version") or "current"
    release = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(release)).strip("_") or "current"
    filename = f"{entity_type}_ids_{release}.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/target-id-qa/api/review-queue")
def target_id_qa_review_queue(
    entity_type: str = "",
    divergence_type: str = "",
    source: str = "",
    status: str = "open",
    auto_decision: str = "",
    triage_category: str = "",
    review_group: str = "",
    q: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_target_graph()
    return build_target_review_queue(
        data,
        entity_type=entity_type,
        divergence_type=divergence_type,
        source=source,
        status=status,
        auto_decision=auto_decision,
        triage_category=triage_category,
        review_group=review_group,
        q=q,
        page=page,
        per_page=per_page,
    )


def _default_target_review_file() -> str:
    # Use the loaded data's qc_dir if available
    try:
        data = _load_target_graph()
        if data.qc_dir:
            return str(Path(data.qc_dir) / "review" / "app_completed" / "target_app_review_decisions.tsv")
    except Exception:
        pass
    if not _target_graph_dir:
        return ""
    graph_dir = Path(_target_graph_dir)
    if graph_dir.name == "app_graph":
        return str(graph_dir.parent / "qc" / "review" / "app_completed" / "target_app_review_decisions.tsv")
    return str(graph_dir / "qc" / "review" / "app_completed" / "target_app_review_decisions.tsv")


@app.post("/target-id-qa/api/review-decisions")
async def target_id_qa_save_review_decisions(request: Request):
    """Append app-entered target review decisions to an IFX Harmonizers intake TSV."""
    data = _load_target_graph()
    if not data.review_can_write:
        raise HTTPException(status_code=403, detail="Target review decisions are disabled in public/read-only mode.")
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    rows = validate_and_build_target_review_rows(payload)
    review_path = _default_target_review_file()
    if not review_path:
        raise HTTPException(status_code=500, detail="No target graph dir configured.")
    path = append_target_review_rows(rows, review_path)
    # Mark resolved in memory so they disappear from the "open" view
    resolved_ids = {r["Registry ID"] for r in rows if r.get("Registry ID")}
    mark_rows_resolved_by_registry_ids(data, resolved_ids)
    return {
        "saved": len(rows),
        "review_file": path,
        "message": "Saved review decision(s). Import with IFX Harmonizers target_review_intake.",
    }


@app.post("/target-id-qa/api/batch-review")
async def target_id_qa_batch_review(request: Request):
    """Batch-apply a review decision to all items matching a triage category."""
    data = _load_target_graph()
    if not data.review_can_write:
        raise HTTPException(status_code=403, detail="Target batch review is disabled in public/read-only mode.")
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    triage_cat = str(payload.get("triage_category") or "").strip()
    review_group = str(payload.get("review_group") or "").strip()
    decision = str(payload.get("review_decision") or "").strip()
    if not triage_cat and not review_group:
        raise HTTPException(status_code=400, detail="triage_category or review_group is required.")
    if not decision:
        raise HTTPException(status_code=400, detail="review_decision is required.")

    rows = build_batch_review_payload(
        data,
        triage_category=triage_cat,
        review_decision=decision,
        entity_type=str(payload.get("entity_type") or ""),
        source=str(payload.get("source") or ""),
        status=str(payload.get("status") or "open"),
        review_group=review_group,
        reviewed_by=str(payload.get("reviewed_by") or ""),
    )
    if not rows:
        return {"saved": 0, "message": "No matching items found."}
    review_path = _default_target_review_file()
    if not review_path:
        raise HTTPException(status_code=500, detail="No target graph dir configured.")
    path = append_target_review_rows(rows, review_path)
    # Mark resolved in memory so they disappear from "open" view
    mark_rows_resolved_by_triage(
        data,
        triage_category=triage_cat,
        entity_type=str(payload.get("entity_type") or ""),
        source=str(payload.get("source") or ""),
        status=str(payload.get("status") or "open"),
        review_group=review_group,
    )
    return {
        "saved": len(rows),
        "review_file": path,
        "message": f"Batch-applied '{decision}' to {len(rows):,} {review_group or triage_cat} items.",
    }


@app.get("/target-id-qa/api/review-template")
def target_id_qa_review_template():
    """Return a TSV with just the header row of target review intake columns."""
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=TARGET_REVIEW_INTAKE_COLUMNS, delimiter="\t")
    writer.writeheader()
    return StreamingResponse(
        io.BytesIO(buf.getvalue().encode("utf-8")),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="target_review_intake_template.tsv"'},
    )


# ---------------------------------------------------------------------------
# Target Harmonizer — Version Comparison
# ---------------------------------------------------------------------------

def _candidate_target_version_roots() -> list[Path]:
    roots: list[Path] = []
    if _target_graph_dir:
        graph_path = Path(_target_graph_dir)
        if graph_path.name == "app_graph" and graph_path.parent.name.startswith("v"):
            roots.append(graph_path.parent.parent)
        elif graph_path.name.startswith("v") or graph_path.parent.name == "target_app_graph":
            roots.append(graph_path.parent)
        else:
            roots.append(graph_path)
    roots.append(TARGET_APP_GRAPH_BUNDLED_DIR)

    seen: set[str] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            key = str(root.resolve())
        except OSError:
            key = str(root)
        if key not in seen and root.is_dir():
            seen.add(key)
            unique.append(root)
    return unique


def _normalize_target_version(version: str) -> str:
    text = (version or "").strip()
    if text.startswith("v."):
        text = text[2:]
    elif text.startswith("v"):
        text = text[1:]
    return text


def _target_version_sort_key(version: str) -> tuple[int, ...]:
    parts = [int(part) for part in re.findall(r"\d+", _normalize_target_version(version))]
    return tuple(parts or [0])


def _discover_target_versions() -> list[dict]:
    version_by_key: dict[str, dict] = {}
    current_path = Path(_target_graph_dir).resolve() if _target_graph_dir else None

    for root in _candidate_target_version_roots():
        if not root.is_dir():
            continue
        for child in root.iterdir():
            graph_child = child
            directory_name = child.name
            if child.is_dir() and not (child / "manifest.json").exists() and (child / "app_graph" / "manifest.json").exists():
                graph_child = child / "app_graph"
            if not graph_child.is_dir() or not (graph_child / "manifest.json").exists():
                continue
            if not (graph_child / "target_nodes.tsv").exists():
                continue
            with open(graph_child / "manifest.json", encoding="utf-8") as fh:
                try:
                    manifest = json.load(fh)
                except (json.JSONDecodeError, ValueError) as exc:
                    logger.warning("Skipping %s: bad manifest.json: %s", graph_child, exc)
                    continue
            version = _normalize_target_version(
                manifest.get("target_release_version", "") or manifest.get("pipeline_version", "") or manifest.get("version", "") or directory_name
            )
            if not version:
                continue
            resolved = graph_child.resolve()
            entry = {
                "version": version,
                "directory_name": directory_name,
                "path": str(graph_child),
                "updated_at": manifest.get("updated_last") or manifest.get("generated_at", ""),
                "node_rows": manifest.get("files", {}).get("target_nodes.tsv", {}).get("rows"),
                "edge_rows": manifest.get("files", {}).get("target_edges.tsv", {}).get("rows"),
                "current": bool(current_path and resolved == current_path),
            }
            existing = version_by_key.get(version)
            if existing is None or entry["current"]:
                version_by_key[version] = entry

    # Fallback: match by directory name
    if _target_graph_dir and not any(v.get("current") for v in version_by_key.values()):
        graph_dir_name = Path(_target_graph_dir).name
        for entry in version_by_key.values():
            if entry["directory_name"] == graph_dir_name:
                entry["current"] = True
                break

    return sorted(
        version_by_key.values(),
        key=lambda row: _target_version_sort_key(row["version"]),
    )


def _default_target_version_pair(versions: list[dict]) -> tuple[str, str]:
    if not versions:
        return "", ""
    current = next((v["version"] for v in versions if v.get("current")), versions[-1]["version"])
    baseline = ""
    if len(versions) >= 2:
        current_index = next(
            (idx for idx, v in enumerate(versions) if v["version"] == current),
            len(versions) - 1,
        )
        baseline_index = max(0, current_index - 1)
        baseline = versions[baseline_index]["version"]
    return baseline, current


def _target_version_dir(version: str) -> Path | None:
    wanted = _normalize_target_version(version)
    for entry in _discover_target_versions():
        if entry["version"] == wanted:
            return Path(entry["path"])
    return None


def _target_precomputed_diff_path(from_version: str, to_version: str, to_dir: Path) -> Path | None:
    from_token = f"v{_normalize_target_version(from_version)}"
    to_token = f"v{_normalize_target_version(to_version)}"
    candidates = [
        to_dir / f"target_version_diff_{from_token}_to_{to_token}.json",
        to_dir.parent / to_token / f"target_version_diff_{from_token}_to_{to_token}.json",
    ]
    for candidate in candidates:
        if candidate.exists():
            return candidate
    return None


@app.get("/target-id-qa/api/versions")
def target_id_qa_versions():
    """List bundled/versioned target app_graph datasets available for diffs."""
    versions = _discover_target_versions()
    default_from, default_to = _default_target_version_pair(versions)
    return {
        "versions": versions,
        "default_from_version": default_from,
        "default_to_version": default_to,
    }


@app.get("/target-id-qa/api/version-diff")
def target_id_qa_version_diff(
    from_version: Optional[str] = None,
    to_version: Optional[str] = None,
):
    """Compute delta between two versioned target app_graph datasets."""
    if not _target_graph_dir:
        raise HTTPException(status_code=500, detail="No --target-graph-dir configured.")
    versions = _discover_target_versions()
    default_from, default_to = _default_target_version_pair(versions)
    from_version = _normalize_target_version(from_version or default_from)
    to_version = _normalize_target_version(to_version or default_to)
    if not from_version or not to_version:
        raise HTTPException(status_code=404, detail="No version pair available for comparison.")
    if from_version == to_version:
        raise HTTPException(status_code=400, detail="Choose two different target graph versions.")

    from_dir = _target_version_dir(from_version)
    to_dir = _target_version_dir(to_version)
    if from_dir is None:
        raise HTTPException(status_code=404, detail=f"Target graph version not found: {from_version}")
    if to_dir is None:
        raise HTTPException(status_code=404, detail=f"Target graph version not found: {to_version}")

    precomputed = _target_precomputed_diff_path(from_version, to_version, to_dir)
    if precomputed:
        with open(precomputed, encoding="utf-8") as fh:
            return json.load(fh)

    baseline = load_target_version_data(from_dir)
    current = _load_target_graph() if to_dir.resolve() == Path(_target_graph_dir).resolve() else load_target_version_data(to_dir)
    return compute_target_version_diff(current, baseline)


# ---------------------------------------------------------------------------
# Pharos TDL + STRING PPI (Target Harmonizer)
# ---------------------------------------------------------------------------


def _get_target_graph_db():
    """Get ArangoDB connection to target_graph, raising a clear HTTP error on failure."""
    try:
        return get_db("target_graph")
    except Exception as exc:
        raise HTTPException(
            status_code=503,
            detail=f"Cannot connect to target_graph ArangoDB: {exc}. "
                   "Ensure ArangoDB credentials are configured (--credentials).",
        ) from exc


@app.get("/target-id-qa/api/pharos-tdl")
async def target_pharos_tdl():
    """Return Pharos TDL distribution for all protein targets in the loaded graph."""
    data = _load_target_graph()
    db = _get_target_graph_db()
    pharos = await run_in_threadpool(load_pharos_tdl_data, db)
    return match_targets_to_tdl(data, pharos)


@app.post("/target-id-qa/api/pharos-tdl/refresh")
async def target_pharos_tdl_refresh():
    """Clear the cached Pharos TDL data and reload from ArangoDB."""
    invalidate_pharos_tdl_cache()
    data = _load_target_graph()
    db = _get_target_graph_db()
    pharos = await run_in_threadpool(load_pharos_tdl_data, db)
    result = match_targets_to_tdl(data, pharos)
    result["refreshed"] = True
    return result


@app.post("/target-id-qa/api/pharos-ppi")
async def target_pharos_ppi(request: Request):
    """Build a STRING PPI network colored by Pharos TDL for a gene list."""
    body = await request.json()
    genes = body.get("genes") or []
    if isinstance(genes, str):
        genes = [g.strip() for g in genes.split("\n") if g.strip()]
    genes = sorted(set(genes))
    if not genes:
        raise HTTPException(status_code=400, detail="No genes provided.")
    if len(genes) > 500:
        raise HTTPException(status_code=400, detail="Too many genes (max 500).")
    required_score = int(body.get("required_score", 700))
    db = _get_target_graph_db()
    pharos = await run_in_threadpool(load_pharos_tdl_data, db)
    ppi_edges = await run_in_threadpool(fetch_string_ppi, genes, required_score)
    return build_pharos_ppi_payload(pharos, genes, ppi_edges)


@app.post("/target-id-qa/api/pharos-tdl/lookup")
async def target_pharos_tdl_lookup(request: Request):
    """Look up Pharos TDL classification for a list of gene symbols."""
    body = await request.json()
    genes = body.get("genes") or []
    if isinstance(genes, str):
        genes = [g.strip() for g in genes.split("\n") if g.strip()]
    genes = [g.strip() for g in genes if g.strip()]
    if not genes:
        raise HTTPException(status_code=400, detail="No genes provided.")
    if len(genes) > 2000:
        raise HTTPException(status_code=400, detail="Too many genes (max 2000).")
    db = _get_target_graph_db()
    pharos = await run_in_threadpool(load_pharos_tdl_data, db)
    return lookup_tdl_by_symbols(pharos, genes)


# ---------------------------------------------------------------------------
# Variant Harmonizer Explorer
# ---------------------------------------------------------------------------


@app.get("/variant-id-qa", response_class=HTMLResponse)
def variant_id_qa(request: Request, ids: str = "", tab: str = ""):
    selected_ids = " | ".join([part for part in re.split(r"[\s,|]+", ids or "") if part])
    default_tab = "graph" if selected_ids else "dashboard"
    manifest = _variant_manifest_snapshot()
    return templates.TemplateResponse(request, "variant_id_qa.html", {
        "request": request,
        "selected_ids": selected_ids,
        "active_tab": tab or default_tab,
        "manifest": manifest,
        "variant_version_display": _variant_version_display(manifest),
        "variant_updated_last": _variant_updated_last(manifest),
    })


@app.get("/variant-id-qa/api/stats")
def variant_id_qa_stats():
    data = _load_variant_graph()
    return compute_variant_stats(data)


# ---------------------------------------------------------------------------
# Variant Version Diff
# ---------------------------------------------------------------------------

VARIANT_APP_GRAPH_BUNDLED_DIR = Path(__file__).resolve().parent / "data" / "variant_app_graph"


def _candidate_variant_version_roots() -> list[Path]:
    roots: list[Path] = []
    if _variant_graph_dir:
        graph_path = Path(_variant_graph_dir)
        if graph_path.name == "app_graph" and graph_path.parent.name.startswith("v"):
            roots.append(graph_path.parent.parent)
        elif graph_path.name.startswith("v") or graph_path.parent.name in {"variant_app_graph", "variant_data"}:
            roots.append(graph_path.parent)
        else:
            roots.append(graph_path)
    roots.append(VARIANT_APP_GRAPH_BUNDLED_DIR)
    seen: set[str] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            key = str(root.resolve())
        except OSError:
            continue
        if key not in seen:
            seen.add(key)
            unique.append(root)
    return unique


def _normalize_variant_version(version: str) -> str:
    text = (version or "").strip()
    if text.startswith("v."):
        text = text[2:]
    elif text.startswith("v"):
        text = text[1:]
    return text


def _variant_version_sort_key(version: str) -> tuple[int, ...]:
    parts = [int(part) for part in re.findall(r"\d+", _normalize_variant_version(version))]
    return tuple(parts or [0])


def _discover_variant_versions() -> list[dict]:
    version_by_key: dict[str, dict] = {}
    current_path = Path(_variant_graph_dir).resolve() if _variant_graph_dir else None

    for root in _candidate_variant_version_roots():
        if not root.is_dir():
            continue
        for child in root.iterdir():
            graph_child = child
            directory_name = child.name
            if child.is_dir() and not (child / "manifest.json").exists() and (child / "app_graph" / "manifest.json").exists():
                graph_child = child / "app_graph"
            if not graph_child.is_dir() or not (graph_child / "manifest.json").exists():
                continue
            if not (graph_child / "variant_nodes.tsv").exists():
                continue
            with open(graph_child / "manifest.json", encoding="utf-8") as fh:
                try:
                    manifest = json.load(fh)
                except (json.JSONDecodeError, ValueError):
                    continue
            version = _normalize_variant_version(
                manifest.get("variant_harmonizer_version", "") or manifest.get("pipeline_version", "") or manifest.get("version", "") or directory_name
            )
            if not version:
                continue
            resolved = graph_child.resolve()
            entry = {
                "version": version,
                "directory_name": directory_name,
                "path": str(graph_child),
                "updated_at": manifest.get("updated_last") or manifest.get("generated_at", ""),
                "node_rows": manifest.get("files", {}).get("variant_nodes.tsv", {}).get("rows"),
                "current": bool(current_path and resolved == current_path),
            }
            existing = version_by_key.get(version)
            if existing is None or entry["current"]:
                version_by_key[version] = entry

    if _variant_graph_dir and not any(v.get("current") for v in version_by_key.values()):
        graph_dir_name = Path(_variant_graph_dir).name
        for entry in version_by_key.values():
            if entry["directory_name"] == graph_dir_name:
                entry["current"] = True
                break

    return sorted(version_by_key.values(), key=lambda row: _variant_version_sort_key(row["version"]))


def _default_variant_version_pair(versions: list[dict]) -> tuple[str, str]:
    if not versions:
        return "", ""
    current = next((v["version"] for v in versions if v.get("current")), versions[-1]["version"])
    baseline = ""
    if len(versions) >= 2:
        current_index = next(
            (idx for idx, v in enumerate(versions) if v["version"] == current),
            len(versions) - 1,
        )
        baseline_index = max(0, current_index - 1)
        baseline = versions[baseline_index]["version"]
    return baseline, current


def _variant_version_dir(version: str) -> Path | None:
    wanted = _normalize_variant_version(version)
    for entry in _discover_variant_versions():
        if entry["version"] == wanted:
            return Path(entry["path"])
    return None


@app.get("/variant-id-qa/api/versions")
def variant_id_qa_versions():
    """List bundled/versioned variant app_graph datasets available for diffs."""
    versions = _discover_variant_versions()
    default_from, default_to = _default_variant_version_pair(versions)
    return {
        "versions": versions,
        "default_from_version": default_from,
        "default_to_version": default_to,
    }


@app.get("/variant-id-qa/api/version-diff")
def variant_id_qa_version_diff(
    from_version: Optional[str] = None,
    to_version: Optional[str] = None,
):
    """Compute delta between two versioned variant app_graph datasets."""
    if not _variant_graph_dir:
        raise HTTPException(status_code=500, detail="No --variant-graph-dir configured.")
    versions = _discover_variant_versions()
    default_from, default_to = _default_variant_version_pair(versions)
    from_version = _normalize_variant_version(from_version or default_from)
    to_version = _normalize_variant_version(to_version or default_to)
    if not from_version or not to_version:
        raise HTTPException(status_code=404, detail="No version pair available for comparison.")
    if from_version == to_version:
        raise HTTPException(status_code=400, detail="Choose two different variant graph versions.")

    from_dir = _variant_version_dir(from_version)
    to_dir = _variant_version_dir(to_version)
    if from_dir is None:
        raise HTTPException(status_code=404, detail=f"Variant graph version not found: {from_version}")
    if to_dir is None:
        raise HTTPException(status_code=404, detail=f"Variant graph version not found: {to_version}")

    # Check for precomputed diff JSON
    for candidate in [
        to_dir / f"version_diff_from_v{from_version}.json",
        to_dir.parent / f"v{to_version}" / f"version_diff_from_v{from_version}.json",
    ]:
        if candidate.exists():
            with open(candidate, encoding="utf-8") as fh:
                return json.load(fh)

    baseline = load_variant_version_data(from_dir)
    current = _load_variant_graph() if to_dir.resolve() == Path(_variant_graph_dir).resolve() else load_variant_version_data(to_dir)
    return compute_variant_version_diff(current, baseline)


def _default_variant_review_file() -> str:
    if not _variant_graph_dir:
        return ""
    graph_dir = Path(_variant_graph_dir)
    if graph_dir.name == "app_graph":
        return str(graph_dir.parent / "review_intake" / "variant_app_review_decisions.tsv")
    return str(graph_dir / "review_intake" / "variant_app_review_decisions.tsv")


def _resolved_variant_review_file() -> str:
    return _variant_review_file or _default_variant_review_file()


def _variant_review_payload_rows(payload: dict) -> list[dict[str, str]]:
    raw_rows = payload.get("decisions") or payload.get("rows") or []
    if isinstance(raw_rows, dict):
        raw_rows = [raw_rows]
    if not isinstance(raw_rows, list):
        raise HTTPException(status_code=400, detail="decisions must be a list.")

    reviewed_by = str(payload.get("reviewed_by") or os.getenv("USER") or "app_review").strip()
    reviewed_at = datetime.now(timezone.utc).isoformat()
    allowed = {opt["value"] for opt in VARIANT_REVIEW_DECISION_OPTIONS}
    rows: list[dict[str, str]] = []
    for raw in raw_rows:
        if not isinstance(raw, dict):
            continue
        registry_id = str(raw.get("registry_id") or raw.get("Registry ID") or "").strip()
        decision = str(raw.get("review_decision") or raw.get("Human decision") or "").strip()
        if not registry_id:
            raise HTTPException(status_code=400, detail="Each variant review decision needs registry_id.")
        if not decision:
            raise HTTPException(status_code=400, detail="Each variant review decision needs review_decision.")
        if decision not in allowed:
            raise HTTPException(status_code=400, detail=f"Unsupported variant review_decision: {decision}")

        rows.append({
            "App review ID": str(raw.get("app_review_id") or raw.get("App review ID") or uuid.uuid4()).strip(),
            "Registry ID": registry_id,
            "Variant ID": str(raw.get("variant_id") or raw.get("Variant ID") or "").strip(),
            "Variant Key": str(raw.get("variant_key") or raw.get("Variant Key") or "").strip(),
            "Primary ID": str(raw.get("primary_id") or raw.get("Primary ID") or "").strip(),
            "Name": str(raw.get("name") or raw.get("Name") or "").strip(),
            "Gene": str(raw.get("gene_symbol") or raw.get("gene_curie") or raw.get("Gene") or "").strip(),
            "Issue type": str(raw.get("issue_type") or raw.get("Issue type") or "").strip(),
            "Issue label": str(raw.get("issue_label") or raw.get("Issue label") or "").strip(),
            "Review group": str(raw.get("review_group") or raw.get("Review group") or "").strip(),
            "Severity": str(raw.get("severity") or raw.get("Severity") or "").strip(),
            "Source": str(raw.get("source") or raw.get("Source") or "").strip(),
            "Source value": str(raw.get("source_value") or raw.get("Source value") or "").strip(),
            "Target ID": str(raw.get("target_id") or raw.get("Target ID") or "").strip(),
            "Target label": str(raw.get("target_label") or raw.get("Target label") or "").strip(),
            "Target category": str(raw.get("target_category") or raw.get("Target category") or "").strip(),
            "Target mapping status": str(raw.get("target_mapping_status") or raw.get("Target mapping status") or "").strip(),
            "Target mapping note": str(raw.get("target_mapping_note") or raw.get("Target mapping note") or "").strip(),
            "Affected variants": str(raw.get("affected_variant_count") or raw.get("Affected variants") or "").strip(),
            "Affected edges": str(raw.get("affected_edge_count") or raw.get("Affected edges") or "").strip(),
            "Example variant IDs": str(raw.get("example_variant_ids") or raw.get("Example variant IDs") or "").strip(),
            "Example variant keys": str(raw.get("example_variant_keys") or raw.get("Example variant keys") or "").strip(),
            "Recommended action": str(raw.get("recommended_action") or raw.get("Recommended action") or "").strip(),
            "Human decision": decision,
            "Resolution": str(raw.get("resolution") or raw.get("Resolution") or "").strip(),
            "Reviewer notes": str(raw.get("notes") or raw.get("Reviewer notes") or "").strip(),
            "Reviewed by": reviewed_by,
            "Reviewed at": reviewed_at,
        })
    if not rows:
        raise HTTPException(status_code=400, detail="No variant review decisions provided.")
    return rows


def _append_variant_review_rows(rows: list[dict[str, str]]) -> str:
    review_file = _resolved_variant_review_file()
    if not review_file:
        raise HTTPException(status_code=500, detail="No variant review file configured.")
    path = Path(review_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with _variant_review_lock:
        exists = path.exists() and path.stat().st_size > 0
        with open(path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=VARIANT_REVIEW_INTAKE_COLUMNS,
                delimiter="\t",
                extrasaction="ignore",
            )
            if not exists:
                writer.writeheader()
            writer.writerows(rows)
    return str(path)


@app.get("/variant-id-qa/api/review-queue")
def variant_id_qa_review_queue(
    issue_type: str = "",
    review_group: str = "",
    source: str = "",
    severity: str = "",
    status: str = "open",
    q: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_variant_graph()
    return build_variant_review_queue(
        data,
        issue_type=issue_type,
        review_group=review_group,
        source=source,
        severity=severity,
        status=status,
        q=q,
        page=page,
        per_page=per_page,
    )


@app.get("/variant-id-qa/download-review-template")
def variant_id_qa_download_review_template(
    issue_type: str = "",
    review_group: str = "",
    source: str = "",
    severity: str = "",
    status: str = "open",
    q: str = "",
):
    data = _load_variant_graph()
    tsv_content = export_variant_review_intake_template(
        data,
        issue_type=issue_type,
        review_group=review_group,
        source=source,
        severity=severity,
        status=status,
        q=q,
    )
    return StreamingResponse(
        io.StringIO(tsv_content),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="variant_review_intake_template.tsv"'},
    )


@app.post("/variant-id-qa/api/review-decisions")
async def variant_id_qa_review_decisions(payload: dict):
    rows = _variant_review_payload_rows(payload)
    path = _append_variant_review_rows(rows)
    return {
        "ok": True,
        "saved_rows": len(rows),
        "review_file": path,
        "message": "Saved review decision(s). Import with IFX Harmonizers variant_qc.",
    }


@app.get("/variant-id-qa/api/search")
def variant_id_qa_search(
    q: str = "",
    source: str = "",
    gene: str = "",
    significance: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_variant_graph()
    return search_variants(
        data, q=q, source=source, gene=gene,
        significance=significance, page=page, per_page=per_page,
    )


@app.get("/variant-id-qa/api/graph")
def variant_id_qa_graph(ids: str = "", include_overlap: bool = True, include_disagreements: bool = True):
    data = _load_variant_graph()
    return build_variant_graph_payload(
        data,
        ids,
        include_overlap=include_overlap,
        include_disagreements=include_disagreements,
    )


@app.get("/variant-id-qa/download-filtered")
def variant_id_qa_download_filtered(
    q: str = "",
    source: str = "",
    gene: str = "",
    significance: str = "",
    columns: str = "",
    format: str = "tsv",
):
    data = _load_variant_graph()
    fmt = "csv" if format == "csv" else "tsv"
    selected_columns = [col.strip() for col in columns.split(",") if col.strip()]
    content = export_variants(
        data, q=q, source=source, gene=gene, significance=significance,
        fmt=fmt, columns=selected_columns,
    )
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    release = (data.manifest.get("version") or data.manifest.get("variant_release_version") or "current")
    release = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(release)).strip("_") or "current"
    filename = f"ODIN_{release}_variants.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )

@app.get("/variant-id-qa/api/myvariant-fields")
def variant_id_qa_myvariant_fields():
    """Return the full MyVariant.info field catalog for the frontend selector."""
    return get_myvariant_field_catalog()


@app.post("/variant-id-qa/api/resolve")
async def variant_id_qa_resolve(body: dict):
    """Batch resolve + enrich variant queries against local graph and MyVariant.info."""
    queries = body.get("queries", [])
    if not queries:
        raise HTTPException(status_code=400, detail="No queries provided")
    if len(queries) > 50:
        raise HTTPException(status_code=400, detail="Maximum 50 queries per request")
    data = _load_variant_graph()
    myvariant_fields = body.get("myvariant_fields")
    if myvariant_fields is not None and not isinstance(myvariant_fields, list):
        myvariant_fields = None
    try:
        workers = min(int(body.get("workers", 4)), 8)
    except (ValueError, TypeError):
        workers = 4
    enable_myvariant = body.get("enable_myvariant", True)
    if isinstance(enable_myvariant, str):
        enable_myvariant = enable_myvariant.lower() not in ("false", "0", "no", "")
    return variant_resolve_and_enrich(
        queries=queries,
        data=data,
        enable_myvariant=enable_myvariant,
        myvariant_fields=myvariant_fields,
        workers=workers,
    )


@app.get("/variant-id-qa/api/resolve-quick")
async def variant_id_qa_resolve_quick(q: str = ""):
    """Quick local-only resolution (no MyVariant enrichment, instant response)."""
    if not q.strip():
        return {"results": [], "stats": {"total": 0}}
    queries = [s.strip() for s in q.split("|") if s.strip()][:50]
    data = _load_variant_graph()
    return variant_resolve_and_enrich(
        queries=queries,
        data=data,
        enable_myvariant=False,
    )


# ---------------------------------------------------------------------------
# Drug Harmonizer Explorer
# ---------------------------------------------------------------------------


@app.get("/drug-id-qa", response_class=HTMLResponse)
def drug_id_qa(request: Request, ids: str = "", tab: str = ""):
    selected_ids = " | ".join([part for part in re.split(r"[\s,|]+", ids or "") if part])
    default_tab = "graph" if selected_ids else "dashboard"
    manifest = _drug_manifest_snapshot()
    resolver_index_available = False
    try:
        drug_data = _load_drug_graph()
        drug_total = f"{len(drug_data.nodes):,}"
        resolver_index_available = bool(
            drug_data.resolver_index_path and drug_data.resolver_index_path.exists()
        )
    except Exception:
        if _drug_resolver_index:
            raise
        drug_total = "412,000+"
    return templates.TemplateResponse(request, "drug_id_qa.html", {
        "request": request,
        "selected_ids": selected_ids,
        "active_tab": tab or default_tab,
        "manifest": manifest,
        "drug_version_display": _drug_version_display(manifest),
        "drug_updated_last": _drug_updated_last(manifest),
        "drug_total_drugs": drug_total,
        "drug_resolver_index_available": resolver_index_available,
    })


@app.get("/drug-id-qa/api/stats")
def drug_id_qa_stats():
    data = _load_drug_graph()
    return compute_drug_stats(data)


@app.get("/drug-id-qa/api/source-versions")
def drug_id_qa_source_versions():
    data = _load_drug_graph()
    return {"sources": summarize_drug_source_versions(data)}


# ---------------------------------------------------------------------------
# Drug Version Diff
# ---------------------------------------------------------------------------

DRUG_APP_GRAPH_BUNDLED_DIR = Path(__file__).resolve().parent / "data" / "drug_app_graph"


def _candidate_drug_version_roots() -> list[Path]:
    roots: list[Path] = []
    if _drug_graph_dir:
        graph_path = Path(_drug_graph_dir)
        if graph_path.name == "app_graph" and graph_path.parent.name.startswith("v"):
            roots.append(graph_path.parent.parent)
        elif graph_path.name.startswith("v") or graph_path.parent.name in {"drug_app_graph", "drug_data"}:
            roots.append(graph_path.parent)
        else:
            roots.append(graph_path)
    roots.append(DRUG_APP_GRAPH_BUNDLED_DIR)
    seen: set[str] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            key = str(root.resolve())
        except OSError:
            continue
        if key not in seen:
            seen.add(key)
            unique.append(root)
    return unique


def _normalize_drug_version(version: str) -> str:
    text = (version or "").strip()
    if text.startswith("v."):
        text = text[2:]
    elif text.startswith("v"):
        text = text[1:]
    return text


def _drug_version_sort_key(version: str) -> tuple[int, ...]:
    parts = [int(part) for part in re.findall(r"\d+", _normalize_drug_version(version))]
    return tuple(parts or [0])


def _discover_drug_versions() -> list[dict]:
    version_by_key: dict[str, dict] = {}
    current_path = Path(_drug_graph_dir).resolve() if _drug_graph_dir else None

    for root in _candidate_drug_version_roots():
        if not root.is_dir():
            continue
        for child in root.iterdir():
            graph_child = child
            directory_name = child.name
            if child.is_dir() and not (child / "manifest.json").exists() and (child / "app_graph" / "manifest.json").exists():
                graph_child = child / "app_graph"
            if not graph_child.is_dir() or not (graph_child / "manifest.json").exists():
                continue
            if not (graph_child / "drug_nodes.tsv").exists():
                continue
            with open(graph_child / "manifest.json", encoding="utf-8") as fh:
                try:
                    manifest = json.load(fh)
                except (json.JSONDecodeError, ValueError):
                    continue
            version = _normalize_drug_version(
                manifest.get("drug_harmonizer_version", "") or manifest.get("pipeline_version", "") or manifest.get("version", "") or directory_name
            )
            if not version:
                continue
            resolved = graph_child.resolve()
            entry = {
                "version": version,
                "directory_name": directory_name,
                "path": str(graph_child),
                "updated_at": manifest.get("updated_last") or manifest.get("generated_at", ""),
                "node_rows": manifest.get("files", {}).get("drug_nodes.tsv", {}).get("rows"),
                "current": bool(current_path and resolved == current_path),
            }
            existing = version_by_key.get(version)
            if existing is None or existing.get("diff_only") or entry["current"]:
                version_by_key[version] = entry

            # A compact audited release diff can preserve the prior release
            # comparison without bundling another large preview graph.
            for diff_path in graph_child.glob("version_diff_from_v*.json"):
                try:
                    diff = json.loads(diff_path.read_text(encoding="utf-8"))
                except (OSError, json.JSONDecodeError, ValueError):
                    continue
                baseline_version = _normalize_drug_version(str(diff.get("baseline_version") or ""))
                if not baseline_version or baseline_version in version_by_key:
                    continue
                summary = diff.get("summary") or {}
                version_by_key[baseline_version] = {
                    "version": baseline_version,
                    "directory_name": "audited diff baseline",
                    "path": "",
                    "updated_at": diff.get("generated_at", ""),
                    "node_rows": summary.get("drugs_old"),
                    "current": False,
                    "diff_only": True,
                }

    if _drug_graph_dir and not any(v.get("current") for v in version_by_key.values()):
        graph_dir_name = Path(_drug_graph_dir).name
        for entry in version_by_key.values():
            if entry["directory_name"] == graph_dir_name:
                entry["current"] = True
                break

    return sorted(version_by_key.values(), key=lambda row: _drug_version_sort_key(row["version"]))


def _default_drug_version_pair(versions: list[dict]) -> tuple[str, str]:
    if not versions:
        return "", ""
    current = next((v["version"] for v in versions if v.get("current")), versions[-1]["version"])
    baseline = ""
    if len(versions) >= 2:
        current_index = next(
            (idx for idx, v in enumerate(versions) if v["version"] == current),
            len(versions) - 1,
        )
        baseline_index = max(0, current_index - 1)
        baseline = versions[baseline_index]["version"]
    return baseline, current


def _drug_version_dir(version: str) -> Path | None:
    wanted = _normalize_drug_version(version)
    for entry in _discover_drug_versions():
        if entry["version"] == wanted:
            return Path(entry["path"]) if entry.get("path") else None
    return None


@app.get("/drug-id-qa/api/versions")
def drug_id_qa_versions():
    """List bundled/versioned drug app_graph datasets available for diffs."""
    versions = _discover_drug_versions()
    default_from, default_to = _default_drug_version_pair(versions)
    return {
        "versions": versions,
        "default_from_version": default_from,
        "default_to_version": default_to,
    }


@app.get("/drug-id-qa/api/version-diff")
def drug_id_qa_version_diff(
    from_version: Optional[str] = None,
    to_version: Optional[str] = None,
):
    """Compute delta between two versioned drug app_graph datasets."""
    if not _drug_graph_dir:
        raise HTTPException(status_code=500, detail="No --drug-graph-dir configured.")
    versions = _discover_drug_versions()
    default_from, default_to = _default_drug_version_pair(versions)
    from_version = _normalize_drug_version(from_version or default_from)
    to_version = _normalize_drug_version(to_version or default_to)
    if not from_version or not to_version:
        raise HTTPException(status_code=404, detail="No version pair available for comparison.")
    if from_version == to_version:
        raise HTTPException(status_code=400, detail="Choose two different drug graph versions.")

    to_dir = _drug_version_dir(to_version)
    if to_dir is None:
        raise HTTPException(status_code=404, detail=f"Drug graph version not found: {to_version}")

    # Check for precomputed diff JSON
    for candidate in [
        to_dir / f"version_diff_from_v{from_version}.json",
        to_dir.parent / f"v{to_version}" / f"version_diff_from_v{from_version}.json",
    ]:
        if candidate.exists():
            with open(candidate, encoding="utf-8") as fh:
                return json.load(fh)

    from_dir = _drug_version_dir(from_version)
    if from_dir is None:
        raise HTTPException(status_code=404, detail=f"Drug graph version not found: {from_version}")
    baseline = load_drug_version_data(from_dir)
    current = _load_drug_graph() if to_dir.resolve() == Path(_drug_graph_dir).resolve() else load_drug_version_data(to_dir)
    return compute_drug_version_diff(current, baseline)


@app.get("/drug-id-qa/api/search")
def drug_id_qa_search(
    q: str = "",
    source: str = "",
    tier: str = "",
    key_type: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_drug_graph()
    return search_drugs(data, q=q, source=source, tier=tier, key_type=key_type, page=page, per_page=per_page)


@app.get("/drug-id-qa/api/graph")
def drug_id_qa_graph(ids: str = "", include_targets: bool = True):
    data = _load_drug_graph()
    return build_drug_graph_payload(data, ids, include_targets=include_targets)


def _default_drug_review_file() -> str:
    if not _drug_graph_dir:
        return ""
    graph_dir = Path(_drug_graph_dir)
    if graph_dir.name == "app_graph":
        return str(graph_dir.parent / "review_intake" / "drug_app_review_decisions.tsv")
    return str(graph_dir / "review_intake" / "drug_app_review_decisions.tsv")


def _resolved_drug_review_file() -> str:
    return _drug_review_file or _default_drug_review_file()


def _drug_review_payload_rows(payload: dict) -> list[dict[str, str]]:
    raw_rows = payload.get("decisions") or payload.get("rows") or []
    if isinstance(raw_rows, dict):
        raw_rows = [raw_rows]
    if not isinstance(raw_rows, list):
        raise HTTPException(status_code=400, detail="decisions must be a list.")

    reviewed_by = str(payload.get("reviewed_by") or os.getenv("USER") or "app_review").strip()
    reviewed_at = datetime.now(timezone.utc).isoformat()
    allowed = {opt["value"] for opt in DRUG_REVIEW_DECISION_OPTIONS}
    rows: list[dict[str, str]] = []
    for raw in raw_rows:
        if not isinstance(raw, dict):
            continue
        registry_id = str(raw.get("registry_id") or raw.get("Registry ID") or "").strip()
        decision = str(raw.get("review_decision") or raw.get("Human decision") or "").strip()
        if not registry_id:
            raise HTTPException(status_code=400, detail="Each drug review decision needs registry_id.")
        if not decision:
            raise HTTPException(status_code=400, detail="Each drug review decision needs review_decision.")
        if decision not in allowed:
            raise HTTPException(status_code=400, detail=f"Unsupported drug review_decision: {decision}")

        rows.append({
            "App review ID": str(raw.get("app_review_id") or raw.get("App review ID") or uuid.uuid4()).strip(),
            "Registry ID": registry_id,
            "Drug ID": str(raw.get("drug_id") or raw.get("Drug ID") or "").strip(),
            "Drug Name": str(raw.get("standard_name") or raw.get("Drug Name") or "").strip(),
            "Primary ID": str(raw.get("primary_id") or raw.get("Primary ID") or "").strip(),
            "Source-standard primary ID": str(raw.get("source_standard_primary_id") or raw.get("Source-standard primary ID") or "").strip(),
            "Issue type": str(raw.get("issue_type") or raw.get("Issue type") or "").strip(),
            "Issue label": str(raw.get("issue_label") or raw.get("Issue label") or "").strip(),
            "Severity": str(raw.get("severity") or raw.get("Severity") or "").strip(),
            "Source": str(raw.get("source") or raw.get("Source") or "").strip(),
            "Source value": str(raw.get("source_value") or raw.get("Source value") or "").strip(),
            "Consensus value": str(raw.get("consensus_value") or raw.get("Consensus value") or "").strip(),
            "NodeNorm canonical CURIE": str(raw.get("nodenorm_canonical_curie") or raw.get("NodeNorm canonical CURIE") or "").strip(),
            "NodeNorm canonical label": str(raw.get("nodenorm_canonical_label") or raw.get("NodeNorm canonical label") or "").strip(),
            "NodeNorm validation status": str(raw.get("nodenorm_validation_status") or raw.get("NodeNorm validation status") or "").strip(),
            "Recommended action": str(raw.get("recommended_action") or raw.get("Recommended action") or "").strip(),
            "Human decision": decision,
            "Resolution": str(raw.get("resolution") or raw.get("Resolution") or "").strip(),
            "Reviewer notes": str(raw.get("notes") or raw.get("Reviewer notes") or "").strip(),
            "Reviewed by": reviewed_by,
            "Reviewed at": reviewed_at,
        })
    if not rows:
        raise HTTPException(status_code=400, detail="No drug review decisions provided.")
    return rows


def _append_drug_review_rows(rows: list[dict[str, str]]) -> str:
    review_file = _resolved_drug_review_file()
    if not review_file:
        raise HTTPException(status_code=500, detail="No drug review file configured.")
    path = Path(review_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with _drug_review_lock:
        exists = path.exists() and path.stat().st_size > 0
        with open(path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=DRUG_REVIEW_INTAKE_COLUMNS,
                delimiter="\t",
                extrasaction="ignore",
            )
            if not exists:
                writer.writeheader()
            writer.writerows(rows)
    return str(path)


@app.get("/drug-id-qa/api/review-queue")
def drug_id_qa_review_queue(
    issue_type: str = "",
    source: str = "",
    severity: str = "",
    status: str = "open",
    q: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_drug_graph()
    return build_drug_review_queue(
        data,
        issue_type=issue_type,
        source=source,
        severity=severity,
        status=status,
        q=q,
        page=page,
        per_page=per_page,
    )


@app.get("/drug-id-qa/download-review-template")
def drug_id_qa_download_review_template(
    issue_type: str = "",
    source: str = "",
    severity: str = "",
    status: str = "open",
    q: str = "",
):
    data = _load_drug_graph()
    tsv_content = export_drug_review_intake_template(
        data,
        issue_type=issue_type,
        source=source,
        severity=severity,
        status=status,
        q=q,
    )
    return StreamingResponse(
        io.BytesIO(tsv_content.encode("utf-8")),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="drug_review_intake_template.tsv"'},
    )


@app.post("/drug-id-qa/api/review-decisions")
async def drug_id_qa_save_review_decisions(request: Request):
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    rows = _drug_review_payload_rows(payload)
    path = _append_drug_review_rows(rows)
    return {
        "saved": len(rows),
        "review_file": path,
        "message": "Saved review decision(s). Import on the next IFX Harmonizers DRUGS drug_qc run.",
    }


@app.get("/drug-id-qa/download-filtered")
def drug_id_qa_download_filtered(
    q: str = "",
    source: str = "",
    tier: str = "",
    key_type: str = "",
    columns: str = "",
    format: str = "tsv",
):
    data = _load_drug_graph()
    fmt = "csv" if format == "csv" else "tsv"
    selected_columns = [col.strip() for col in columns.split(",") if col.strip()]
    content = export_drugs(data, q=q, source=source, tier=tier, key_type=key_type, fmt=fmt, columns=selected_columns)
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    release = data.manifest.get("version") or data.manifest.get("drug_release_version") or "current"
    release = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(release)).strip("_") or "current"
    filename = f"ODIN_{release}_drugs.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


@app.get("/drug-id-qa/download/{filename}")
def drug_id_qa_download_file(filename: str):
    if not _drug_graph_dir:
        raise HTTPException(status_code=500, detail="No drug graph dir configured.")
    safe_name = re.sub(r"[^A-Za-z0-9_.-]", "", filename)
    path = Path(_drug_graph_dir) / safe_name
    if not path.exists() or not path.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {safe_name}")
    return FileResponse(str(path), filename=safe_name)


@app.get("/drug-id-qa/download-sssom")
def drug_id_qa_download_sssom(sources: str = ""):
    data = _load_drug_graph()
    include_sources = [s.strip() for s in sources.split(",") if s.strip()] or None
    content = iter_drug_sssom_bytes(data, include_sources=include_sources)
    return StreamingResponse(
        io.BytesIO(content),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="drug_sssom.tsv"'},
    )


@app.get("/drug-id-qa/api/neighbors/{drug_id:path}")
def drug_id_qa_neighbors(drug_id: str):
    data = _load_drug_graph()
    return find_drug_neighbors(data, drug_id)


@app.get("/api/v1/drug/resolve/{query_id:path}")
def api_drug_resolve(query_id: str):
    data = _load_drug_graph()
    return resolve_drug(data, query_id)


@app.get("/drug-id-qa/api/ncats-properties")
def drug_id_qa_ncats_properties():
    """Return the full NCATS property catalog for the frontend selector."""
    return get_ncats_property_catalog()


def _drug_resolver_payload_bool(body: dict, key: str, default: bool = True) -> bool:
    value = body.get(key, default)
    if value is None:
        return default
    if isinstance(value, str):
        return value.strip().lower() not in ("false", "0", "no", "")
    return bool(value)


def _drug_resolver_request_payload(body: dict) -> dict[str, Any]:
    queries = body.get("queries", [])
    if not isinstance(queries, list):
        queries = [queries]
    queries = [str(query).strip() for query in queries if str(query or "").strip()]
    if not queries:
        raise HTTPException(status_code=400, detail="No queries provided")
    if len(queries) > DRUG_RESOLVER_MAX_QUERIES:
        raise HTTPException(status_code=400, detail=f"Maximum {DRUG_RESOLVER_MAX_QUERIES} queries per request")

    ncats_props = body.get("ncats_props")  # None -> use defaults
    if ncats_props is not None and not isinstance(ncats_props, list):
        ncats_props = None
    try:
        workers = min(max(int(body.get("workers", 4)), 1), 8)
    except (TypeError, ValueError):
        workers = 4

    return {
        "queries": queries,
        "enable_ncats": _drug_resolver_payload_bool(body, "enable_ncats", True),
        "enable_pubchem": _drug_resolver_payload_bool(body, "enable_pubchem", True),
        "enable_pharos": _drug_resolver_payload_bool(body, "enable_pharos", False),
        "enable_inxight": _drug_resolver_payload_bool(body, "enable_inxight", False),
        "enable_openfda": _drug_resolver_payload_bool(body, "enable_openfda", False),
        "enable_chebi": _drug_resolver_payload_bool(body, "enable_chebi", False),
        "workers": workers,
        "delay": 0.15,
        "ncats_props": ncats_props,
    }


def _drug_resolver_source_labels(payload: dict[str, Any]) -> list[str]:
    labels: list[str] = []
    if payload.get("enable_ncats"):
        labels.append("NCATS Resolver")
    if payload.get("enable_pubchem"):
        labels.append("PubChem (local misses)")
    if payload.get("enable_pharos"):
        labels.append("Pharos")
    if payload.get("enable_inxight"):
        labels.append("Inxight/GSRS")
    if payload.get("enable_openfda"):
        labels.append("openFDA")
    if payload.get("enable_chebi"):
        labels.append("ChEBI")
    return labels


def _cleanup_drug_resolver_jobs() -> None:
    cutoff = time.time() - _DRUG_RESOLVER_JOB_TTL_SECONDS
    with _drug_resolver_jobs_lock:
        old_ids = [
            job_id for job_id, job in _drug_resolver_jobs.items()
            if float(job.get("updated_epoch") or job.get("created_epoch") or 0) < cutoff
        ]
        for job_id in old_ids:
            _drug_resolver_jobs.pop(job_id, None)


def _update_drug_resolver_job(job_id: str, **updates: Any) -> None:
    with _drug_resolver_jobs_lock:
        job = _drug_resolver_jobs.get(job_id)
        if job is None:
            return
        job.update(updates)
        job["updated_at"] = datetime.now(timezone.utc).isoformat()
        job["updated_epoch"] = time.time()


def _run_drug_resolver_job(job_id: str, payload: dict[str, Any]) -> None:
    source_labels = _drug_resolver_source_labels(payload)
    source_text = ", ".join(source_labels) if source_labels else "local graph only"
    _update_drug_resolver_job(
        job_id,
        status="running",
        started_at=datetime.now(timezone.utc).isoformat(),
        message=f"Started local harmonizer lookup plus {source_text}.",
    )

    def progress(event: dict[str, Any]) -> None:
        total = int(event.get("total") or 0)
        completed = int(event.get("completed") or 0)
        stage = str(event.get("stage") or "running")
        if stage == "local_complete":
            progress_pct = 5 if source_labels else 100
        elif total:
            progress_pct = min(99, 5 + int((completed / total) * 94))
        else:
            progress_pct = 0
        _update_drug_resolver_job(
            job_id,
            stage=stage,
            completed=completed,
            total=total,
            progress=progress_pct,
            local_resolved=event.get("local_resolved"),
            sources_found=event.get("sources_found") or {},
            sources_failed=event.get("sources_failed") or {},
            message=event.get("message") or "",
        )

    try:
        data = _load_drug_graph()
        result = resolve_and_enrich(data, progress_callback=progress, **payload)
        _update_drug_resolver_job(
            job_id,
            status="complete",
            stage="complete",
            completed=len(payload.get("queries") or []),
            total=len(payload.get("queries") or []),
            progress=100,
            completed_at=datetime.now(timezone.utc).isoformat(),
            result=result,
            message="Drug resolver enrichment complete.",
        )
    except Exception as exc:
        _update_drug_resolver_job(
            job_id,
            status="failed",
            stage="failed",
            completed_at=datetime.now(timezone.utc).isoformat(),
            error=str(exc),
            message=f"Drug resolver job failed: {exc}",
        )


def _enqueue_drug_resolver_job(payload: dict[str, Any]) -> dict[str, Any]:
    _cleanup_drug_resolver_jobs()
    job_id = f"{datetime.now(timezone.utc).strftime('%Y%m%dT%H%M%SZ')}-{uuid.uuid4().hex[:8]}"
    source_labels = _drug_resolver_source_labels(payload)
    job = {
        "id": job_id,
        "status": "queued",
        "stage": "queued",
        "created_at": datetime.now(timezone.utc).isoformat(),
        "created_epoch": time.time(),
        "updated_at": None,
        "updated_epoch": time.time(),
        "started_at": None,
        "completed_at": None,
        "progress": 0,
        "completed": 0,
        "total": len(payload.get("queries") or []),
        "local_resolved": None,
        "sources_requested": source_labels,
        "sources_found": {},
        "sources_failed": {},
        "message": "Queued drug resolver job.",
        "error": None,
        "result": None,
    }
    with _drug_resolver_jobs_lock:
        _drug_resolver_jobs[job_id] = job
    thread = threading.Thread(target=_run_drug_resolver_job, args=(job_id, payload), daemon=True)
    thread.start()
    return {"job_id": job_id, "status": "queued", "total": job["total"], "sources_requested": source_labels}


@app.post("/drug-id-qa/api/resolve")
async def drug_id_qa_resolve(body: dict):
    """Batch resolve + enrich drug queries against local graph and live APIs."""
    payload = _drug_resolver_request_payload(body)
    data = _load_drug_graph()
    return await run_in_threadpool(
        resolve_and_enrich,
        data,
        **payload,
    )


@app.post("/drug-id-qa/api/resolve-job")
async def drug_id_qa_resolve_job(body: dict):
    """Run a larger drug resolver batch in the background with pollable progress."""
    payload = _drug_resolver_request_payload(body)
    return _enqueue_drug_resolver_job(payload)


@app.get("/drug-id-qa/api/resolve-job/{job_id}")
def drug_id_qa_resolve_job_status(job_id: str):
    _cleanup_drug_resolver_jobs()
    with _drug_resolver_jobs_lock:
        job = _drug_resolver_jobs.get(job_id)
        if not job:
            raise HTTPException(status_code=404, detail="Resolver job not found")
        return dict(job)


@app.get("/drug-id-qa/api/resolve-quick")
async def drug_id_qa_resolve_quick(q: str = ""):
    """Quick local-only resolution (no enrichment, instant response)."""
    if not q.strip():
        return {"results": [], "stats": {"total": 0}}
    queries = [s.strip() for s in q.split("|") if s.strip()][:DRUG_RESOLVER_MAX_QUERIES]
    data = _load_drug_graph()
    return await run_in_threadpool(
        resolve_and_enrich,
        data,
        queries=queries,
        enable_ncats=False,
        enable_pubchem=False,
        enable_pharos=False,
        enable_inxight=False,
        enable_openfda=False,
        enable_chebi=False,
    )


# ---------------------------------------------------------------------------
# Pathway Harmonizer Explorer
# ---------------------------------------------------------------------------

PATHWAY_APP_GRAPH_BUNDLED_DIR = Path(__file__).resolve().parent / "data" / "pathway_app_graph"


def _load_pathway_graph():
    if not _pathway_graph_dir:
        raise HTTPException(status_code=500, detail="No --pathway-graph-dir configured.")
    return load_pathway_graph_data(_pathway_graph_dir)


def _pathway_manifest_snapshot() -> dict[str, Any]:
    if not _pathway_graph_dir:
        return {}
    manifest_path = Path(_pathway_graph_dir) / "manifest.json"
    if not manifest_path.exists():
        return {}
    try:
        with open(manifest_path, encoding="utf-8") as fh:
            manifest = json.load(fh)
    except (OSError, json.JSONDecodeError):
        return {}
    return manifest if isinstance(manifest, dict) else {}


def _pathway_version_display(manifest: dict[str, Any]) -> str:
    raw = manifest.get("version") or "—"
    text = str(raw).strip()
    if not text or text == "—":
        return "—"
    if text.startswith("v"):
        return text
    return f"v{text}"


def _pathway_updated_last(manifest: dict[str, Any]) -> str:
    raw = manifest.get("updated_last") or manifest.get("generated_at") or "—"
    text = str(raw).strip()
    return text[:10] if text and text != "—" else "—"


@app.get("/pathway-id-qa", response_class=HTMLResponse)
def pathway_id_qa(request: Request, ids: str = "", tab: str = ""):
    selected_ids = " | ".join([part for part in re.split(r"[\s,|]+", ids or "") if part])
    default_tab = "graph" if selected_ids else "dashboard"
    manifest = _pathway_manifest_snapshot()
    return templates.TemplateResponse(request, "pathway_id_qa.html", {
        "request": request,
        "selected_ids": selected_ids,
        "active_tab": tab or default_tab,
        "manifest": manifest,
        "pathway_version_display": _pathway_version_display(manifest),
        "pathway_updated_last": _pathway_updated_last(manifest),
        "pathway_columns": DEFAULT_PATHWAY_COLUMNS,
    })


@app.get("/pathway-id-qa/api/stats")
def pathway_id_qa_stats():
    data = _load_pathway_graph()
    return compute_pathway_stats(data)


@app.get("/pathway-id-qa/api/search")
def pathway_id_qa_search(
    q: str = "",
    source: str = "",
    category: str = "",
    min_sources: int = 0,
    page: int = 1,
    per_page: int = 50,
):
    data = _load_pathway_graph()
    return search_pathways(data, q=q, source=source, min_sources=min_sources, category=category, page=page, per_page=per_page)


@app.get("/pathway-id-qa/api/graph")
def pathway_id_qa_graph(ids: str = ""):
    data = _load_pathway_graph()
    return build_pathway_graph_payload(data, ids)


@app.get("/pathway-id-qa/download-filtered")
def pathway_id_qa_download_filtered(
    q: str = "",
    source: str = "",
    min_sources: int = 0,
    columns: str = "",
    format: str = "tsv",
):
    data = _load_pathway_graph()
    fmt = "csv" if format == "csv" else "tsv"
    selected_columns = [col.strip() for col in columns.split(",") if col.strip()]
    content = export_pathways(data, q=q, source=source, min_sources=min_sources, fmt=fmt, columns=selected_columns)
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    release = data.manifest.get("version") or "current"
    release = re.sub(r"[^A-Za-z0-9_.-]+", "_", str(release)).strip("_") or "current"
    filename = f"ODIN_{release}_pathways.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ---------------------------------------------------------------------------
# Pathway Resolver
# ---------------------------------------------------------------------------


@app.post("/pathway-id-qa/api/resolve")
async def pathway_id_qa_resolve(request: Request):
    """Resolve pathway names/IDs to IFXPathway records."""
    data = _load_pathway_graph()
    payload = await request.json()
    queries = payload.get("queries", [])
    if not isinstance(queries, list):
        raise HTTPException(status_code=400, detail="queries must be a list of strings")
    results = resolve_pathway_queries(data, queries[:500])
    return {"results": results, "total": len(results)}


@app.post("/pathway-id-qa/download-resolver-results")
async def pathway_id_qa_download_resolver_results(
    request: Request,
    format: str = "tsv",
):
    """Download resolver results as TSV/CSV."""
    data = _load_pathway_graph()
    form = await request.form()
    raw_queries = str(form.get("queries", ""))
    queries = [q.strip() for q in raw_queries.split("\n") if q.strip()]
    results = resolve_pathway_queries(data, queries[:500])
    fmt = "csv" if format == "csv" else "tsv"
    content = export_pathway_resolver_results(results, fmt=fmt)
    media_type = "text/csv" if fmt == "csv" else "text/tab-separated-values"
    filename = f"pathway_resolver_results.{fmt}"
    return StreamingResponse(
        io.StringIO(content),
        media_type=media_type,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ---------------------------------------------------------------------------
# Pathway Cross-Entity Links
# ---------------------------------------------------------------------------


@app.get("/pathway-id-qa/api/cross-links")
def pathway_id_qa_cross_links(pathway_id: str = ""):
    """Return disease and drug connections for a pathway via its gene members."""
    data = _load_pathway_graph()
    # Try to get other entity data if loaded
    disease_data = None
    drug_data = None
    try:
        if _disease_graph_dir:
            disease_data = load_disease_graph_data(_disease_graph_dir)
    except Exception:
        pass
    try:
        if _drug_graph_dir:
            from src.qa_browser.drug_id_graph import load_drug_graph_data
            drug_data = load_drug_graph_data(_drug_graph_dir)
    except Exception:
        pass
    return build_cross_entity_summary(
        data, pathway_id,
        disease_data=disease_data,
        drug_data=drug_data,
    )


# ---------------------------------------------------------------------------
# Pathway Version Diff
# ---------------------------------------------------------------------------


def _candidate_pathway_version_roots() -> list[Path]:
    roots: list[Path] = []
    if _pathway_graph_dir:
        graph_path = Path(_pathway_graph_dir)
        if graph_path.name == "app_graph" and graph_path.parent.name.startswith("v"):
            roots.append(graph_path.parent.parent)
        elif graph_path.name.startswith("v") or graph_path.parent.name in {"pathway_app_graph", "pathway_data"}:
            roots.append(graph_path.parent)
        else:
            roots.append(graph_path)
    roots.append(PATHWAY_APP_GRAPH_BUNDLED_DIR)
    seen: set[str] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            key = str(root.resolve())
        except OSError:
            continue
        if key not in seen:
            seen.add(key)
            unique.append(root)
    return unique


def _normalize_pathway_version(version: str) -> str:
    text = (version or "").strip()
    if text.startswith("v."):
        text = text[2:]
    elif text.startswith("v"):
        text = text[1:]
    return text


def _pathway_version_sort_key(version: str) -> tuple[int, ...]:
    parts = [int(part) for part in re.findall(r"\d+", _normalize_pathway_version(version))]
    return tuple(parts or [0])


def _discover_pathway_versions() -> list[dict]:
    version_by_key: dict[str, dict] = {}
    current_path = Path(_pathway_graph_dir).resolve() if _pathway_graph_dir else None

    for root in _candidate_pathway_version_roots():
        if not root.is_dir():
            continue
        for child in root.iterdir():
            graph_child = child
            directory_name = child.name
            if child.is_dir() and not (child / "manifest.json").exists() and (child / "app_graph" / "manifest.json").exists():
                graph_child = child / "app_graph"
            if not graph_child.is_dir() or not (graph_child / "manifest.json").exists():
                continue
            if not (graph_child / "pathway_nodes.tsv").exists():
                continue
            with open(graph_child / "manifest.json", encoding="utf-8") as fh:
                try:
                    manifest = json.load(fh)
                except (json.JSONDecodeError, ValueError):
                    continue
            # The enclosing vX.Y.Z directory is the release identity. Older
            # pathway packages incorrectly wrote a fixed app-schema version
            # (0.1.0) into every manifest, so the directory must take priority.
            version = _normalize_pathway_version(directory_name)
            if not re.fullmatch(r"\d+\.\d+\.\d+", version):
                version = _normalize_pathway_version(manifest.get("version", ""))
            if not version:
                continue
            resolved = graph_child.resolve()
            entry = {
                "version": version,
                "directory_name": directory_name,
                "path": str(graph_child),
                "updated_at": manifest.get("updated_last") or manifest.get("generated_at", ""),
                "node_rows": manifest.get("nodes", {}).get("count") if isinstance(manifest.get("nodes"), dict) else None,
                "current": bool(current_path and resolved == current_path),
            }
            existing = version_by_key.get(version)
            if existing is None or entry["current"]:
                version_by_key[version] = entry

    if _pathway_graph_dir and not any(v.get("current") for v in version_by_key.values()):
        graph_dir_name = Path(_pathway_graph_dir).name
        for entry in version_by_key.values():
            if entry["directory_name"] == graph_dir_name:
                entry["current"] = True
                break

    return sorted(version_by_key.values(), key=lambda row: _pathway_version_sort_key(row["version"]))


def _default_pathway_version_pair(versions: list[dict]) -> tuple[str, str]:
    if not versions:
        return "", ""
    current = next((v["version"] for v in versions if v.get("current")), versions[-1]["version"])
    baseline = ""
    if len(versions) >= 2:
        current_index = next(
            (idx for idx, v in enumerate(versions) if v["version"] == current),
            len(versions) - 1,
        )
        baseline_index = max(0, current_index - 1)
        baseline = versions[baseline_index]["version"]
    return baseline, current


def _pathway_version_dir(version: str) -> Path | None:
    wanted = _normalize_pathway_version(version)
    for entry in _discover_pathway_versions():
        if entry["version"] == wanted:
            return Path(entry["path"])
    return None


@app.get("/pathway-id-qa/api/versions")
def pathway_id_qa_versions():
    """List bundled/versioned pathway app_graph datasets available for diffs."""
    versions = _discover_pathway_versions()
    default_from, default_to = _default_pathway_version_pair(versions)
    return {
        "versions": versions,
        "default_from_version": default_from,
        "default_to_version": default_to,
    }


@app.get("/pathway-id-qa/api/version-diff")
def pathway_id_qa_version_diff(
    from_version: Optional[str] = None,
    to_version: Optional[str] = None,
):
    """Compute delta between two versioned pathway app_graph datasets."""
    if not _pathway_graph_dir:
        raise HTTPException(status_code=500, detail="No --pathway-graph-dir configured.")
    versions = _discover_pathway_versions()
    default_from, default_to = _default_pathway_version_pair(versions)
    from_version = _normalize_pathway_version(from_version or default_from)
    to_version = _normalize_pathway_version(to_version or default_to)
    if not from_version or not to_version:
        raise HTTPException(status_code=404, detail="No version pair available for comparison.")
    if from_version == to_version:
        raise HTTPException(status_code=400, detail="Choose two different pathway graph versions.")

    from_dir = _pathway_version_dir(from_version)
    to_dir = _pathway_version_dir(to_version)
    if from_dir is None:
        raise HTTPException(status_code=404, detail=f"Pathway graph version not found: {from_version}")
    if to_dir is None:
        raise HTTPException(status_code=404, detail=f"Pathway graph version not found: {to_version}")

    # Prefer the versioned release artifact owned by the harmonizer. Dynamic
    # computation remains a compatibility fallback for older packages.
    precomputed = to_dir / f"version_diff_from_v{from_version}.json"
    if precomputed.is_file():
        try:
            payload = json.loads(precomputed.read_text(encoding="utf-8"))
            if payload.get("diff_schema_version") == "1.0" and payload.get("summary"):
                return payload
        except (OSError, json.JSONDecodeError):
            pass
    baseline = load_pathway_version_data(from_dir)
    current = _load_pathway_graph() if to_dir.resolve() == Path(_pathway_graph_dir).resolve() else load_pathway_version_data(to_dir)
    diff = compute_pathway_version_diff(current, baseline)
    diff["computed_fallback"] = True
    return diff


# ---------------------------------------------------------------------------
# Pathway review helpers & routes
# ---------------------------------------------------------------------------

def _default_pathway_review_file() -> str:
    if not _pathway_graph_dir:
        return ""
    graph_dir = Path(_pathway_graph_dir)
    if graph_dir.name == "app_graph":
        return str(graph_dir.parent / "review_intake" / "pathway_app_review_decisions.tsv")
    return str(graph_dir / "review_intake" / "pathway_app_review_decisions.tsv")


def _resolved_pathway_review_file() -> str:
    return _pathway_review_file or _default_pathway_review_file()


def _pathway_review_payload_rows(payload: dict) -> list[dict[str, str]]:
    raw_rows = payload.get("decisions") or payload.get("rows") or []
    if isinstance(raw_rows, dict):
        raw_rows = [raw_rows]
    if not isinstance(raw_rows, list):
        raise HTTPException(status_code=400, detail="decisions must be a list.")

    reviewed_by = str(payload.get("reviewed_by") or os.getenv("USER") or "app_review").strip()
    reviewed_at = datetime.now(timezone.utc).isoformat()
    allowed = {opt["value"] for opt in PATHWAY_REVIEW_DECISION_OPTIONS}
    rows: list[dict[str, str]] = []
    for raw in raw_rows:
        if not isinstance(raw, dict):
            continue
        pathway_id = str(raw.get("ncats_pathway_id") or raw.get("Pathway ID") or "").strip()
        decision = str(raw.get("review_decision") or raw.get("Human decision") or "").strip()
        if not pathway_id:
            raise HTTPException(status_code=400, detail="Each pathway review decision needs ncats_pathway_id.")
        if not decision:
            raise HTTPException(status_code=400, detail="Each pathway review decision needs review_decision.")
        if decision not in allowed:
            raise HTTPException(status_code=400, detail=f"Unsupported pathway review_decision: {decision}")
        rows.append({
            "App review ID": str(raw.get("app_review_id") or raw.get("App review ID") or uuid.uuid4()).strip(),
            "Pathway ID": pathway_id,
            "Pathway Name": str(raw.get("pathway_name") or raw.get("Pathway Name") or "").strip(),
            "Scenario": str(raw.get("scenario") or raw.get("Scenario") or "").strip(),
            "Severity": str(raw.get("severity") or raw.get("Severity") or "").strip(),
            "Detail": str(raw.get("detail") or raw.get("Detail") or "").strip(),
            "Status": str(raw.get("status") or raw.get("Status") or "").strip(),
            "Human decision": decision,
            "Resolution": str(raw.get("resolution") or raw.get("Resolution") or "").strip(),
            "Reviewer notes": str(raw.get("notes") or raw.get("Reviewer notes") or "").strip(),
            "Reviewed by": reviewed_by,
            "Reviewed at": reviewed_at,
        })
    if not rows:
        raise HTTPException(status_code=400, detail="No pathway review decisions provided.")
    return rows


def _append_pathway_review_rows(rows: list[dict[str, str]]) -> str:
    review_file = _resolved_pathway_review_file()
    if not review_file:
        raise HTTPException(status_code=500, detail="No pathway review file configured.")
    path = Path(review_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with _pathway_review_lock:
        exists = path.exists() and path.stat().st_size > 0
        with open(path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=PATHWAY_REVIEW_INTAKE_COLUMNS,
                delimiter="\t",
                extrasaction="ignore",
            )
            if not exists:
                writer.writeheader()
            writer.writerows(rows)
    return str(path)


@app.get("/pathway-id-qa/api/review-queue")
def pathway_id_qa_review_queue(
    scenario: str = "",
    severity: str = "",
    status: str = "open",
    q: str = "",
    page: int = 1,
    per_page: int = 50,
):
    data = _load_pathway_graph()
    return build_pathway_review_queue(
        data,
        scenario=scenario,
        severity=severity,
        status=status,
        q=q,
        page=page,
        per_page=per_page,
    )


@app.get("/pathway-id-qa/download-review-template")
def pathway_id_qa_download_review_template(
    scenario: str = "",
    severity: str = "",
    status: str = "open",
    q: str = "",
):
    data = _load_pathway_graph()
    tsv_content = export_pathway_review_intake_template(
        data,
        scenario=scenario,
        severity=severity,
        status=status,
        q=q,
    )
    return StreamingResponse(
        io.BytesIO(tsv_content.encode("utf-8")),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="pathway_review_intake_template.tsv"'},
    )


@app.post("/pathway-id-qa/api/review-decisions")
async def pathway_id_qa_save_review_decisions(request: Request):
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    rows = _pathway_review_payload_rows(payload)
    path = _append_pathway_review_rows(rows)
    return {
        "saved": len(rows),
        "review_file": path,
        "message": "Saved review decision(s). Import on the next IFX Harmonizers PATHWAYS pathway_qc run.",
    }


@app.post("/pathway-id-qa/api/batch-review")
async def pathway_id_qa_batch_review(request: Request):
    """Batch-apply a review decision to all items matching a scenario."""
    data = _load_pathway_graph()
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    scenario = str(payload.get("scenario") or "").strip()
    decision = str(payload.get("review_decision") or "").strip()
    if not scenario:
        raise HTTPException(status_code=400, detail="scenario is required.")
    if not decision:
        raise HTTPException(status_code=400, detail="review_decision is required.")
    rows = build_batch_pathway_review_payload(
        data,
        scenario=scenario,
        review_decision=decision,
        severity=str(payload.get("severity") or ""),
        status=str(payload.get("status") or "open"),
        reviewed_by=str(payload.get("reviewed_by") or ""),
    )
    if not rows:
        return {"saved": 0, "message": "No matching items found."}
    path = _append_pathway_review_rows(rows)
    mark_pathway_rows_resolved(
        data,
        scenario=scenario,
        severity=str(payload.get("severity") or ""),
        status=str(payload.get("status") or "open"),
    )
    return {
        "saved": len(rows),
        "review_file": path,
        "message": f"Batch-applied '{decision}' to {len(rows):,} {scenario} items.",
    }


@app.get("/disease-id-qa", response_class=HTMLResponse)
def disease_id_qa(request: Request, ids: str = "", tab: str = ""):
    selected_ids = " | ".join(parse_disease_ids(ids))
    error_message = None
    manifest: dict = {}
    if _disease_graph_dir:
        try:
            data = load_disease_graph_data(_disease_graph_dir)
            manifest = data.manifest
        except Exception as exc:
            error_message = str(exc)
    else:
        error_message = "No --disease-graph-dir configured."
    # Default to graph tab when IDs are provided, dashboard otherwise
    active_tab = tab if tab else ("graph" if ids else "dashboard")
    return templates.TemplateResponse(request, "disease_id_qa.html", {
        "request": request,
        "ids": ids,
        "selected_ids": selected_ids,
        "manifest": manifest,
        "flagged_rows": [],
        "error_message": error_message,
        "active_tab": active_tab,
    })


@app.get("/disease-id-qa/api/graph")
def disease_id_qa_graph(ids: str = ""):
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    target_data = None
    if _target_graph_dir:
        try:
            target_data = load_target_graph_data(_target_graph_dir)
        except Exception:
            pass
    return build_disease_graph_payload(data, parse_disease_ids(ids), target_data=target_data)


@app.get("/disease-id-qa/download/{filename}")
def disease_id_qa_download(filename: str):
    """Serve a file from the disease app_graph directory."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if filename not in DOWNLOADABLE_FILES:
        raise HTTPException(status_code=404, detail=f"File not available: {filename}")
    file_path = Path(_disease_graph_dir) / filename
    if not file_path.is_file():
        raise HTTPException(status_code=404, detail=f"File not found: {filename}")
    media_type = "application/json" if filename.endswith(".json") else "text/tab-separated-values"
    return FileResponse(
        path=file_path,
        media_type=media_type,
        filename=filename,
    )


@app.get("/disease-id-qa/download-flagged")
def disease_id_qa_download_flagged():
    """Download concepts that still need review/action as a TSV file."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    tsv_content = export_flagged_tsv(data)
    return StreamingResponse(
        io.BytesIO(tsv_content.encode("utf-8")),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="needs_action_concepts.tsv"'},
    )


def _default_disease_review_file() -> str:
    if not _disease_graph_dir:
        return ""
    graph_dir = Path(_disease_graph_dir)
    if graph_dir.name == "app_graph":
        return str(graph_dir.parent / "qc" / "review" / "app_completed" / "disease_app_review_decisions.tsv")
    return str(graph_dir / "review_intake" / "disease_app_review_decisions.tsv")


def _resolved_disease_review_file() -> str:
    return _disease_review_file or _default_disease_review_file()


def _review_payload_rows(payload: dict) -> list[dict[str, str]]:
    raw_rows = payload.get("decisions") or payload.get("rows") or []
    if isinstance(raw_rows, dict):
        raw_rows = [raw_rows]
    if not isinstance(raw_rows, list):
        raise HTTPException(status_code=400, detail="decisions must be a list.")

    reviewed_by = str(payload.get("reviewed_by") or os.getenv("USER") or "app_review").strip()
    reviewed_at = datetime.now(timezone.utc).isoformat()
    rows: list[dict[str, str]] = []
    for raw in raw_rows:
        if not isinstance(raw, dict):
            continue
        registry_id = str(raw.get("registry_id") or raw.get("Registry ID") or "").strip()
        decision = str(raw.get("review_decision") or raw.get("Review decision") or "").strip()
        if not registry_id:
            raise HTTPException(status_code=400, detail="Each review decision needs registry_id.")
        if not decision:
            raise HTTPException(status_code=400, detail="Each review decision needs review_decision.")
        if decision not in {opt["value"] for opt in REVIEW_DECISION_OPTIONS}:
            raise HTTPException(status_code=400, detail=f"Unsupported review_decision: {decision}")

        rows.append({
            "Registry ID": registry_id,
            "Review decision": decision,
            "Corrected xref / replacement": str(raw.get("replacement") or raw.get("Corrected xref / replacement") or "").strip(),
            "Reviewer notes": str(raw.get("notes") or raw.get("Reviewer notes") or "").strip(),
            "Primary Xref": str(raw.get("primary_xref") or raw.get("Primary Xref") or "").strip(),
            "Disease Name": str(raw.get("standard_name") or raw.get("Disease Name") or "").strip(),
            "NCATS Disease ID": str(raw.get("ncats_disease_id") or raw.get("NCATS Disease ID") or "").strip(),
            "Review category": str(raw.get("category") or raw.get("Review category") or "").strip(),
            "Xref ID": str(raw.get("xref_id") or raw.get("Xref ID") or "").strip(),
            "Xref Namespace": str(raw.get("xref_namespace") or raw.get("Xref Namespace") or "").strip(),
            "Decision type": str(raw.get("decision_type") or raw.get("Decision type") or "").strip(),
            "Scenario ID": str(raw.get("scenario_id") or raw.get("Scenario ID") or "").strip(),
            "Source value": str(raw.get("source_value") or raw.get("Source value") or "").strip(),
            "Consensus value": str(raw.get("consensus_value") or raw.get("Consensus value") or "").strip(),
            "Evidence note": str(raw.get("evidence_note") or raw.get("Evidence note") or "").strip(),
            "Reviewed by": reviewed_by,
            "Reviewed at": reviewed_at,
            "App review ID": str(raw.get("app_review_id") or raw.get("App review ID") or uuid.uuid4()).strip(),
        })
    if not rows:
        raise HTTPException(status_code=400, detail="No review decisions provided.")
    return rows


def _append_disease_review_rows(rows: list[dict[str, str]]) -> str:
    review_file = _resolved_disease_review_file()
    if not review_file:
        raise HTTPException(status_code=500, detail="No disease review file configured.")
    path = Path(review_file)
    path.parent.mkdir(parents=True, exist_ok=True)
    with _disease_review_lock:
        exists = path.exists() and path.stat().st_size > 0
        with open(path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=REVIEW_INTAKE_COLUMNS,
                delimiter="\t",
                extrasaction="ignore",
            )
            if not exists:
                writer.writeheader()
            writer.writerows(rows)
    return str(path)


@app.get("/disease-id-qa/api/review-queue")
def disease_id_qa_review_queue(
    category: str = "",
    source: str = "",
    status: str = "open",
    q: str = "",
    page: int = 1,
    per_page: int = 50,
):
    """Return clustered disease review items for the Review tab."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    return build_review_queue(
        data,
        category=category,
        source=source,
        status=status,
        q=q,
        page=page,
        per_page=per_page,
    )


@app.get("/disease-id-qa/download-review-template")
def disease_id_qa_download_review_template(
    category: str = "",
    source: str = "",
    status: str = "open",
    q: str = "",
):
    """Download an IFX Harmonizers review intake template for the current queue."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    tsv_content = export_review_intake_template(
        data,
        category=category,
        source=source,
        status=status,
        q=q,
    )
    return StreamingResponse(
        io.BytesIO(tsv_content.encode("utf-8")),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="disease_review_intake_template.tsv"'},
    )


@app.post("/disease-id-qa/api/review-decisions")
async def disease_id_qa_save_review_decisions(request: Request):
    """Append app-entered review decisions to an IFX Harmonizers intake TSV."""
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    rows = _review_payload_rows(payload)
    path = _append_disease_review_rows(rows)
    return {
        "saved": len(rows),
        "review_file": path,
        "message": "Saved review decision(s). Import with IFX Harmonizers disease_review_intake.",
    }


@app.get("/disease-id-qa/api/search")
def disease_id_qa_search(
    q: str = "",
    filter: str = "all",
    page: int = 1,
    per_page: int = 50,
    confidence_tier: str = "",
    is_rare: str = "",
    source: str = "",
    quality: str = "",
    disease_type: str = "",
    sort: str = "",
    sort_dir: str = "asc",
):
    """Paginated search across all disease concepts."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    per_page = max(1, min(per_page, 200))
    return search_concepts(
        data, q=q, filter_mode=filter, page=page, per_page=per_page,
        confidence_tier=confidence_tier, is_rare=is_rare, source=source, quality=quality,
        disease_type=disease_type, sort=sort, sort_dir=sort_dir,
    )


@app.get("/disease-id-qa/download-search")
def disease_id_qa_download_search(q: str = "", filter: str = "all"):
    """Download current search results as a TSV file."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    tsv_content = export_search_tsv(data, q=q, filter_mode=filter)
    return StreamingResponse(
        io.BytesIO(tsv_content.encode("utf-8")),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="disease_search_results.tsv"'},
    )


@app.get("/disease-id-qa/api/dashboard")
def disease_id_qa_dashboard():
    """Return dashboard aggregate stats."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    return compute_dashboard_stats(data)


@app.get("/disease-id-qa/api/download-preview")
def disease_id_qa_download_preview(
    q: str = "",
    filter: str = "all",
    confidence_tier: str = "",
    is_rare: str = "",
    source: str = "",
    quality: str = "",
    disease_type: str = "",
    mode: str = "concepts",
    include_sources: str = "",
    exclude_obsolete: str = "",
):
    """Return exact row counts for the selected download table mode."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    if include_sources == "__none__":
        src_set = set()
    else:
        src_set = (
            {s.strip().upper() for s in include_sources.split(",") if s.strip()}
            if include_sources else None
        )
    return compute_download_preview_counts(
        data,
        q=q,
        filter_mode=filter,
        confidence_tier=confidence_tier,
        is_rare=is_rare,
        source=source,
        quality=quality,
        disease_type=disease_type,
        mode=mode,
        include_sources=src_set,
        exclude_obsolete=exclude_obsolete.lower() == "true",
    )


@app.get("/disease-id-qa/api/resolve-quick")
def disease_id_qa_resolve_quick(q: str = "", limit: int = 8, include_obsolete: bool = False):
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if not q.strip():
        raise HTTPException(status_code=400, detail="q parameter required.")
    queries = [part.strip() for part in re.split(r"[\n|]+", q) if part.strip()]
    data = load_disease_graph_data(_disease_graph_dir)
    if len(queries) == 1:
        return resolve_name_candidates(
            data,
            queries[0],
            limit=max(1, min(limit, 50)),
            include_obsolete=include_obsolete,
        )
    return bulk_resolve_names(
        data,
        queries[:500],
        limit=max(1, min(limit, 50)),
        include_obsolete=include_obsolete,
    )


@app.post("/disease-id-qa/api/resolve")
async def disease_id_qa_resolve(request: Request, limit: int = 8, include_obsolete: bool = False):
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    payload = await request.json()
    if isinstance(payload, list):
        queries = payload
    elif isinstance(payload, dict):
        queries = payload.get("queries") or payload.get("names") or []
        limit = payload.get("limit", limit)
        raw_include_obsolete = payload.get("include_obsolete", include_obsolete)
        include_obsolete = str(raw_include_obsolete).lower() in {"1", "true", "yes", "on"}
    else:
        queries = []
    if isinstance(queries, str):
        queries = [line.strip() for line in queries.splitlines() if line.strip()]
    if not isinstance(queries, list) or not queries:
        raise HTTPException(status_code=400, detail="Request body must include a queries list.")
    data = load_disease_graph_data(_disease_graph_dir)
    return bulk_resolve_names(
        data,
        [str(q).strip() for q in queries if str(q).strip()][:500],
        limit=max(1, min(int(limit), 50)),
        include_obsolete=include_obsolete,
    )


@app.get("/disease-id-qa/api/suggest")
def disease_id_qa_suggest(q: str = "", limit: int = 8):
    """Return autocomplete suggestions for the search input."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if len(q.strip()) < 2:
        return {"suggestions": []}
    data = load_disease_graph_data(_disease_graph_dir)
    result = resolve_name_candidates(data, q, limit=max(1, min(limit, 20)))
    suggestions = []
    seen_pxrefs: set = set()
    for cand in result.get("candidates", []):
        pxref = cand.get("primary_xref") or cand.get("pxref", "")
        if pxref in seen_pxrefs:
            continue
        seen_pxrefs.add(pxref)
        concept = data.concepts_by_pxref.get(pxref, {})
        suggestions.append({
            "primary_xref": pxref,
            "standard_name": concept.get("standard_name", ""),
            "confidence_tier": concept.get("confidence_tier", ""),
            "n_sources": concept.get("n_sources", ""),
            "disease_type": concept.get("disease_type", ""),
        })
        if len(suggestions) >= limit:
            break
    return {"suggestions": suggestions}


@app.get("/disease-id-qa/api/review-history")
def disease_id_qa_review_history(primary_xref: str = ""):
    """Return review decision history for a concept from the intake TSV."""
    review_file = _resolved_disease_review_file()
    if not review_file or not Path(review_file).is_file():
        return {"decisions": []}
    decisions = []
    with open(review_file, newline="", encoding="utf-8") as fh:
        reader = csv.DictReader(fh, delimiter="\t")
        for row in reader:
            if primary_xref and row.get("Primary Xref", "") != primary_xref:
                continue
            decisions.append({
                "review_decision": row.get("Review decision", ""),
                "reviewed_by": row.get("Reviewed by", ""),
                "reviewed_at": row.get("Reviewed at", ""),
                "notes": row.get("Reviewer notes", ""),
                "replacement": row.get("Corrected xref / replacement", ""),
                "category": row.get("Review category", ""),
                "xref_id": row.get("Xref ID", ""),
                "xref_namespace": row.get("Xref Namespace", ""),
                "decision_type": row.get("Decision type", ""),
                "evidence_note": row.get("Evidence note", ""),
                "disease_name": row.get("Disease Name", ""),
            })
    decisions.sort(key=lambda d: d.get("reviewed_at", ""), reverse=True)
    return {"decisions": decisions}


@app.get("/disease-id-qa/api/compare")
def disease_id_qa_compare(a: str = "", b: str = ""):
    """Compare two disease concepts side by side."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if not a or not b:
        raise HTTPException(status_code=400, detail="Both 'a' and 'b' parameters are required.")
    data = load_disease_graph_data(_disease_graph_dir)
    concept_a = resolve_concept(data, a)
    concept_b = resolve_concept(data, b)
    if not concept_a:
        raise HTTPException(status_code=404, detail=f"Concept not found: {a}")
    if not concept_b:
        raise HTTPException(status_code=404, detail=f"Concept not found: {b}")
    xrefs_a = {x["xref_id"] for x in concept_a.get("xrefs", []) if x.get("xref_id")}
    xrefs_b = {x["xref_id"] for x in concept_b.get("xrefs", []) if x.get("xref_id")}
    return {
        "concept_a": concept_a,
        "concept_b": concept_b,
        "shared_xrefs": sorted(xrefs_a & xrefs_b),
        "unique_a": sorted(xrefs_a - xrefs_b),
        "unique_b": sorted(xrefs_b - xrefs_a),
    }


@app.get("/disease-id-qa/download-filtered")
def disease_id_qa_download_filtered(
    q: str = "",
    filter: str = "all",
    confidence_tier: str = "",
    is_rare: str = "",
    source: str = "",
    quality: str = "",
    disease_type: str = "",
    columns: str = "",
    mode: str = "concepts",
    format: str = "tsv",
    include_sources: str = "",
    exclude_obsolete: str = "",
):
    """Download filtered disease harmonizer tables with selectable columns."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    col_list = [c.strip() for c in columns.split(",") if c.strip()] if columns else None
    if include_sources == "__none__":
        src_set = set()
    else:
        src_set = (
            {s.strip().upper() for s in include_sources.split(",") if s.strip()}
            if include_sources else None
        )
    excl_obs = exclude_obsolete.lower() == "true"
    mode_key = _normalize_download_mode(mode)
    version = str(data.manifest.get("pipeline_version", "") or "unknown").strip()
    version_tag = f"ODINv{version}" if version and not version.startswith("v") else f"ODIN{version}"
    table_name = {
        "xref_edges": "disease_xref_edges",
        "clinical_descendants": "disease_clinical_descendants",
    }.get(mode_key, "diseases")
    row_filters_active = any([
        q.strip(),
        filter not in ("", "all"),
        confidence_tier.strip(),
        is_rare.strip(),
        source.strip(),
        quality.strip(),
        disease_type.strip(),
        mode_key != "concepts" and bool(src_set),
        mode_key != "concepts" and excl_obs,
    ])
    suffix = "_filtered" if row_filters_active else ""

    if format == "xlsx":
        content = export_filtered_xlsx(
            data, q=q, filter_mode=filter,
            confidence_tier=confidence_tier, is_rare=is_rare, source=source, quality=quality,
            disease_type=disease_type,
            columns=col_list, mode=mode_key, include_sources=src_set,
            exclude_obsolete=excl_obs,
        )
        filename = f"{version_tag}_{table_name}{suffix}.xlsx"
        return StreamingResponse(
            io.BytesIO(content),
            media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            headers={"Content-Disposition": f'attachment; filename="{filename}"'},
        )

    def _stream():
        for line in stream_filtered_download(
            data, q=q, filter_mode=filter,
            confidence_tier=confidence_tier, is_rare=is_rare, source=source, quality=quality,
            disease_type=disease_type,
            columns=col_list, mode=mode_key, fmt=format, include_sources=src_set,
            exclude_obsolete=excl_obs,
        ):
            yield line.encode("utf-8")

    if format == "csv":
        media = "text/csv"
        ext = "csv"
    else:
        media = "text/tab-separated-values"
        ext = "tsv"
    filename = f"{version_tag}_{table_name}{suffix}.{ext}"
    return StreamingResponse(
        _stream(),
        media_type=media,
        headers={"Content-Disposition": f'attachment; filename="{filename}"'},
    )


# ---------------------------------------------------------------------------
# Disease Explorer — SSSOM download
# ---------------------------------------------------------------------------

@app.get("/disease-id-qa/download-sssom")
def disease_id_qa_download_sssom(include_sources: str = ""):
    """Download xref edges in SSSOM TSV format."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    src_set = (
        {s.strip() for s in include_sources.split(",") if s.strip()}
        if include_sources else None
    )
    return StreamingResponse(
        iter_sssom_bytes(data, include_sources=src_set),
        media_type="text/tab-separated-values",
        headers={"Content-Disposition": 'attachment; filename="disease_mappings.sssom.tsv"'},
    )


# ---------------------------------------------------------------------------
# Target Explorer — Public REST API
# ---------------------------------------------------------------------------

@app.get("/api/v1/target/resolve-name")
def api_resolve_target_name(q: str = "", limit: int = 10, target_type: str = ""):
    """Resolve a free-text target name, symbol, ID, or ambiguous target-like label."""
    if not q.strip():
        raise HTTPException(status_code=400, detail="q parameter required.")
    data = _load_target_graph()
    return resolve_target_query(
        data,
        q,
        limit=max(1, min(limit, 50)),
        target_type=target_type,
    )


@app.post("/api/v1/target/bulk-resolve")
async def api_bulk_resolve_target_names(request: Request, limit: int = 8, target_type: str = ""):
    """Resolve a list of target names, symbols, IDs, or ambiguous target-like labels."""
    payload = await request.json()
    if isinstance(payload, list):
        queries = payload
    elif isinstance(payload, dict):
        queries = payload.get("queries") or payload.get("targets") or payload.get("names") or []
        limit = payload.get("limit", limit)
        target_type = payload.get("target_type", target_type)
    else:
        queries = []
    if isinstance(queries, str):
        queries = [line.strip() for line in queries.splitlines() if line.strip()]
    if not isinstance(queries, list) or not queries:
        raise HTTPException(status_code=400, detail="Request body must include a queries list.")
    data = _load_target_graph()
    return bulk_resolve_targets(
        data,
        [str(q).strip() for q in queries if str(q).strip()][:500],
        limit=max(1, min(int(limit), 50)),
        target_type=str(target_type or ""),
    )


# ---------------------------------------------------------------------------
# Disease Explorer — Public REST API
# ---------------------------------------------------------------------------

@app.get("/api/v1/disease/resolve/{curie:path}")
def api_resolve_disease(curie: str):
    """Resolve any disease CURIE to its harmonized concept with all xrefs."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    result = resolve_concept(data, curie)
    if result is None:
        raise HTTPException(status_code=404, detail=f"CURIE not found: {curie}")
    return result


@app.get("/api/v1/disease/search")
def api_search_diseases(q: str = "", limit: int = 20, offset: int = 0,
                        confidence_tier: str = "", source: str = ""):
    """Search disease concepts by name or ID."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    limit = max(1, min(limit, 100))
    offset = max(0, offset)
    page = (offset // limit) + 1
    result = search_concepts(
        data, q=q, page=page, per_page=limit,
        confidence_tier=confidence_tier, source=source,
    )
    return {
        "results": result["rows"],
        "total": result["total"],
        "limit": limit,
        "offset": offset,
    }


@app.get("/api/v1/disease/resolve-name")
def api_resolve_disease_name(q: str = "", limit: int = 10, include_obsolete: bool = False):
    """Resolve a free-text disease name, synonym, xref label, or CURIE."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if not q.strip():
        raise HTTPException(status_code=400, detail="q parameter required.")
    data = load_disease_graph_data(_disease_graph_dir)
    limit = max(1, min(limit, 50))
    return resolve_name_candidates(
        data,
        q,
        limit=limit,
        include_obsolete=include_obsolete,
    )


@app.post("/api/v1/disease/bulk-resolve")
async def api_bulk_resolve_disease_names(
    request: Request,
    limit: int = 5,
    include_obsolete: bool = False,
):
    """Resolve a list of disease names, synonyms, xref labels, or CURIEs."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    payload = await request.json()
    if isinstance(payload, list):
        queries = payload
    elif isinstance(payload, dict):
        queries = payload.get("queries") or payload.get("names") or []
    else:
        queries = []
    if isinstance(queries, str):
        queries = [line.strip() for line in queries.splitlines() if line.strip()]
    if not isinstance(queries, list) or not queries:
        raise HTTPException(status_code=400, detail="Request body must include a queries list.")
    queries = [str(q).strip() for q in queries if str(q).strip()][:500]
    data = load_disease_graph_data(_disease_graph_dir)
    limit = max(1, min(limit, 20))
    return bulk_resolve_names(
        data,
        queries,
        limit=limit,
        include_obsolete=include_obsolete,
    )


@app.get("/api/v1/disease/stats")
def api_disease_stats():
    """Return aggregate harmonization statistics."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    stats = compute_dashboard_stats(data)
    return {
        "total_concepts": stats["total_concepts"],
        "total_edges": stats["total_edges"],
        "total_hierarchy_edges": stats.get("total_hierarchy_edges", 0),
        "total_clinical_descendant_edges": stats.get("total_clinical_descendant_edges", 0),
        "sources_covered": len(stats["source_coverage"]),
        "avg_xrefs_per_concept": stats["avg_xrefs_per_concept"],
        "multi_source_count": stats.get("multi_source_count", 0),
        "multi_source_percent": stats.get("multi_source_percent", 0),
        "source_only_count": stats.get("source_only_count", 0),
        "concepts_with_clinical_codes": stats.get("concepts_with_clinical_codes", 0),
        "concepts_with_clinical_descendants": stats.get("concepts_with_clinical_descendants", 0),
        "clinical_descendant_summary": stats.get("clinical_descendant_summary", {}),
        "rare_count": stats["rare_count"],
        "flagged_count": stats["flagged_count"],
        "needs_action_count": stats.get("needs_action_count", stats["flagged_count"]),
        "needs_review_count": stats.get("needs_review_count", 0),
        "qc_decision_count": stats.get("qc_decision_count", 0),
        "auto_cleared_count": stats.get("auto_cleared_count", 0),
        "confidence_distribution": stats["confidence_distribution"],
        "disease_type_distribution": stats.get("disease_type_distribution", {}),
        "pipeline_version": data.manifest.get("pipeline_version", ""),
        "generated_at": data.manifest.get("generated_at", ""),
    }


@app.get("/api/v1/disease/sources")
def api_disease_sources():
    """Return available source namespaces with coverage counts."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    stats = compute_dashboard_stats(data)
    return {
        "sources": [
            {"namespace": ns, "concept_count": count}
            for ns, count in sorted(stats["source_coverage"].items(), key=lambda x: -x[1])
        ],
    }


# ---------------------------------------------------------------------------
# Disease Explorer — Source Versions (from pipeline metadata)
# ---------------------------------------------------------------------------


@app.get("/disease-id-qa/api/source-versions")
def disease_id_qa_source_versions():
    """Return source database versions from the pipeline metadata catalog."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    graph_dir = Path(_disease_graph_dir)
    # metadata lives at ../metadata/ relative to app_graph/
    if graph_dir.name == "app_graph":
        metadata_dir = graph_dir.parent / "metadata"
    else:
        metadata_dir = graph_dir / "metadata"
    catalog_candidates = [
        graph_dir / "disease_source_catalog.tsv",
        metadata_dir / "disease_source_catalog.tsv",
    ]
    catalog_path = next((path for path in catalog_candidates if path.exists()), None)
    if catalog_path is None:
        return {"sources": [], "catalog_found": False}
    sources = []

    def _clean(value: object) -> str:
        text = str(value or "").strip()
        return "" if text.lower() == "unknown" else text

    def _nested(dct: dict, *keys: str) -> str:
        if not isinstance(dct, dict):
            return ""
        for key in keys:
            value = _clean(dct.get(key))
            if value:
                return value
        return ""

    def _load_row_metadata(row: dict) -> dict:
        meta_name = row.get("metadata_file", "").strip()
        if not meta_name:
            return {}
        meta_path = Path(meta_name)
        if not meta_path.is_absolute():
            meta_path = metadata_dir / "source" / meta_name
        if not meta_path.exists():
            return {}
        try:
            with open(meta_path, "r", encoding="utf-8") as fh:
                return json.load(fh)
        except Exception as exc:
            logger.warning("Failed to read source metadata %s: %s", meta_path, exc)
            return {}

    def _resolve_metadata_path(path_text: str) -> Optional[Path]:
        path_text = _clean(path_text)
        if not path_text:
            return None
        raw_path = Path(path_text)
        if raw_path.is_absolute() and raw_path.exists():
            return raw_path
        for base in [Path.cwd(), *metadata_dir.parents]:
            candidate = base / raw_path
            if candidate.exists():
                return candidate
        return None

    def _raw_obo_version(raw_file: str) -> str:
        raw_path = _resolve_metadata_path(raw_file)
        if not raw_path:
            return ""
        try:
            with open(raw_path, "r", encoding="utf-8", errors="replace") as fh:
                for line_no, line in enumerate(fh, start=1):
                    if line_no > 200 or line.startswith("[Term]"):
                        break
                    line = line.strip()
                    if line.startswith("property_value: owl:versionInfo"):
                        match = re.search(r'"([^"]+)"', line)
                        if match:
                            return match.group(1)
                    if line.startswith("data-version:"):
                        match = re.search(r"\d{4}-\d{2}-\d{2}", line)
                        return match.group(0) if match else line.split(":", 1)[1].strip()
                    if line.startswith("date:"):
                        match = re.search(r"(\d{2}):(\d{2}):(\d{4})", line)
                        if match:
                            day, month, year = match.groups()
                            return f"{year}-{month}-{day}"
        except Exception as exc:
            logger.warning("Failed to read source OBO version from %s: %s", raw_path, exc)
        return ""

    def _freshest_file_date(meta: dict) -> str:
        files = meta.get("files", []) if isinstance(meta, dict) else []
        if not isinstance(files, list):
            return ""
        dates = []
        for file_meta in files:
            if not isinstance(file_meta, dict):
                continue
            date_value = _nested(file_meta, "last_modified_iso", "last_modified")
            if date_value:
                dates.append(date_value)
        return max(dates) if dates else ""

    def _date_only(value: str) -> str:
        value = _clean(value)
        return value[:10] if value else ""

    try:
        import csv
        with open(catalog_path, "r") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                name = row.get("source_name", "").strip()
                if not name:
                    continue
                meta = _load_row_metadata(row)
                source_meta = meta.get("source", {}) if isinstance(meta.get("source"), dict) else {}
                release_info = (
                    meta.get("release_info", {})
                    if isinstance(meta.get("release_info"), dict)
                    else {}
                )
                file_release_date = _freshest_file_date(meta)
                raw_obo_version = _raw_obo_version(meta.get("raw_file", ""))
                version = (
                    _clean(row.get("source_version"))
                    or _nested(meta, "source_version")
                    or _nested(source_meta, "version")
                    or _nested(release_info, "releaseVersion", "version", "releaseDate")
                    or raw_obo_version
                    or _date_only(file_release_date)
                    or "unknown"
                )
                release_date = (
                    _clean(row.get("release_date"))
                    or _nested(meta, "release_date", "last_modified_iso")
                    or _nested(source_meta, "last_modified_iso", "last_modified")
                    or _nested(release_info, "releaseDate", "last_modified_iso", "last_modified")
                    or raw_obo_version
                    or file_release_date
                )
                download_start = (
                    _clean(row.get("download_start"))
                    or _nested(meta, "download_start", "timestamp_start")
                )
                download_date = _date_only(download_start)
                tf_output_rows = _clean(row.get("tf_output_rows"))
                tf_input_rows = _clean(row.get("tf_input_rows"))
                status = _clean(row.get("status"))
                changed = _clean(row.get("changed"))
                url = (
                    _clean(row.get("url"))
                    or _nested(meta, "url", "download_url")
                )
                source_release_url = _nested(meta, "source_release_url", "html_url")
                note = ""
                if name.upper() == "SNOMEDCT":
                    note = "UMLS SNOMEDCT subset; source version is the UMLS release."
                elif name.upper() == "UMLS":
                    note = "UMLS-derived labels/xrefs; access depends on local UMLS credentials."
                elif version == "unknown":
                    note = "Source-native version was not captured by the download metadata."
                sources.append({
                    "name": name,
                    "version": version,
                    "release_date": release_date,
                    "source_release_date": release_date,
                    "download_date": download_date,
                    "odin_download_date": download_date,
                    "transform_input_rows": tf_input_rows,
                    "transform_output_rows": tf_output_rows,
                    "status": status,
                    "changed": changed,
                    "note": note,
                    "url": url,
                    "source_release_url": source_release_url,
                })
    except Exception as exc:
        logger.warning("Failed to read source catalog %s: %s", catalog_path, exc)
        return {"sources": [], "catalog_found": False, "error": str(exc)}
    return {"sources": sources, "catalog_found": True}


# ---------------------------------------------------------------------------
# Disease Explorer — Agreement Matrix & Source Flows
# ---------------------------------------------------------------------------

@app.get("/disease-id-qa/api/agreement-matrix")
def disease_id_qa_agreement_matrix():
    """Return pairwise source agreement matrix for the heatmap."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    return compute_source_agreement_matrix(data)


@app.get("/disease-id-qa/api/source-flows")
def disease_id_qa_source_flows():
    """Return source namespace to confidence tier flow data for the Sankey diagram."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    return compute_source_flow_data(data)


# ---------------------------------------------------------------------------
# Disease Explorer — Version Diff
# ---------------------------------------------------------------------------


def _normalize_disease_version(version: str) -> str:
    text = (version or "").strip()
    if text.startswith("v."):
        text = text[2:]
    elif text.startswith("v"):
        text = text[1:]
    return text


def _disease_version_sort_key(version: str) -> tuple[int, ...]:
    parts = [int(part) for part in re.findall(r"\d+", _normalize_disease_version(version))]
    return tuple(parts or [0])


def _versioned_app_graph_dirs(
    root: Path,
    required_file: str,
    release_fields: tuple[str, ...] = ("version",),
) -> list[Path]:
    """Return valid app-graph directories ordered by manifest semantic version.

    Release directories may be named ``vX.Y.Z`` or ``current``. The manifest
    is authoritative so a newer ``current`` release is not ignored merely
    because older version-named directories are also bundled.
    """
    if not root.is_dir():
        return []

    releases_by_version: dict[tuple[int, ...], Path] = {}
    for path in root.iterdir():
        if not path.is_dir() or not (path / required_file).exists():
            continue
        try:
            manifest = json.loads((path / "manifest.json").read_text(encoding="utf-8"))
        except (OSError, json.JSONDecodeError, ValueError):
            continue
        raw_version = next(
            (manifest.get(key) for key in release_fields if manifest.get(key)),
            path.name if path.name.startswith("v") else "",
        )
        if not raw_version:
            continue
        version_key = _disease_version_sort_key(str(raw_version))
        existing = releases_by_version.get(version_key)
        if existing is None or path.name == "current":
            releases_by_version[version_key] = path

    return [releases_by_version[key] for key in sorted(releases_by_version)]


def _candidate_disease_version_roots() -> list[Path]:
    roots: list[Path] = []
    if _disease_graph_dir:
        graph_path = Path(_disease_graph_dir)
        if graph_path.name == "app_graph" and graph_path.parent.name.startswith("v"):
            roots.append(graph_path.parent.parent)
        elif graph_path.name.startswith("v") or graph_path.parent.name == "disease_app_graph":
            roots.append(graph_path.parent)
        else:
            roots.append(graph_path)
    roots.append(DISEASE_APP_GRAPH_BUNDLED_DIR)

    seen: set[str] = set()
    unique: list[Path] = []
    for root in roots:
        try:
            key = str(root.resolve())
        except OSError:
            key = str(root)
        if key not in seen and root.is_dir():
            seen.add(key)
            unique.append(root)
    return unique


def _discover_disease_versions() -> list[dict]:
    version_by_key: dict[str, dict] = {}
    current_path = Path(_disease_graph_dir).resolve() if _disease_graph_dir else None
    baseline_path = Path(_baseline_graph_dir).resolve() if _baseline_graph_dir else None

    for root in _candidate_disease_version_roots():
        for child in root.iterdir():
            graph_child = child
            directory_name = child.name
            if child.is_dir() and not (child / "manifest.json").exists() and (child / "app_graph" / "manifest.json").exists():
                graph_child = child / "app_graph"
            if not graph_child.is_dir() or not (graph_child / "manifest.json").exists():
                continue
            with open(graph_child / "manifest.json", encoding="utf-8") as fh:
                try:
                    manifest = json.load(fh)
                except (json.JSONDecodeError, ValueError) as exc:
                    logger.warning("Skipping %s: bad manifest.json: %s", graph_child, exc)
                    continue
            version = _normalize_disease_version(
                manifest.get("pipeline_version", "") or manifest.get("version", "") or directory_name
            )
            if not version:
                continue
            resolved = graph_child.resolve()
            entry = {
                "version": version,
                "directory_name": directory_name,
                "path": str(graph_child),
                "updated_at": manifest.get("generated_at", ""),
                "contract_version": manifest.get("contract_version", ""),
                "concept_rows": manifest.get("files", {}).get("disease_concepts.tsv", {}).get("rows"),
                "xref_edge_rows": manifest.get("files", {}).get("disease_xref_edges.tsv", {}).get("rows"),
                "current": bool(current_path and resolved == current_path),
                "baseline": bool(baseline_path and resolved == baseline_path),
            }
            existing = version_by_key.get(version)
            if existing is None or entry["current"] or entry["baseline"]:
                version_by_key[version] = entry

    # Fallback: if nothing matched by resolved path, match by directory name
    if _disease_graph_dir and not any(v.get("current") for v in version_by_key.values()):
        graph_dir_name = Path(_disease_graph_dir).name
        for entry in version_by_key.values():
            if entry["directory_name"] == graph_dir_name:
                entry["current"] = True
                break
    if _baseline_graph_dir and not any(v.get("baseline") for v in version_by_key.values()):
        baseline_dir_name = Path(_baseline_graph_dir).name
        for entry in version_by_key.values():
            if entry["directory_name"] == baseline_dir_name:
                entry["baseline"] = True
                break

    return sorted(
        version_by_key.values(),
        key=lambda row: _disease_version_sort_key(row["version"]),
    )


def _default_disease_version_pair(versions: list[dict]) -> tuple[str, str]:
    if not versions:
        return "", ""
    current = next((v["version"] for v in versions if v.get("current")), versions[-1]["version"])
    baseline = next((v["version"] for v in versions if v.get("baseline")), "")
    if not baseline and len(versions) >= 2:
        current_index = next(
            (idx for idx, v in enumerate(versions) if v["version"] == current),
            len(versions) - 1,
        )
        baseline_index = max(0, current_index - 1)
        baseline = versions[baseline_index]["version"]
    return baseline, current


def _disease_version_dir(version: str) -> Path | None:
    wanted = _normalize_disease_version(version)
    for entry in _discover_disease_versions():
        if entry["version"] == wanted:
            return Path(entry["path"])
    return None


@app.get("/disease-id-qa/api/versions")
def disease_id_qa_versions():
    """List bundled/versioned disease app_graph datasets available for diffs."""
    versions = _discover_disease_versions()
    default_from, default_to = _default_disease_version_pair(versions)
    return {
        "versions": versions,
        "default_from_version": default_from,
        "default_to_version": default_to,
    }


@app.get("/disease-id-qa/api/version-diff")
def disease_id_qa_version_diff(
    from_version: Optional[str] = None,
    to_version: Optional[str] = None,
):
    """Compute delta between two versioned disease app_graph datasets."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    versions = _discover_disease_versions()
    default_from, default_to = _default_disease_version_pair(versions)
    from_version = _normalize_disease_version(from_version or default_from)
    to_version = _normalize_disease_version(to_version or default_to)
    if not from_version or not to_version:
        raise HTTPException(status_code=404, detail="No version pair available for comparison.")
    if from_version == to_version:
        raise HTTPException(status_code=400, detail="Choose two different disease graph versions.")

    from_dir = _disease_version_dir(from_version)
    to_dir = _disease_version_dir(to_version)
    if from_dir is None:
        raise HTTPException(status_code=404, detail=f"Disease graph version not found: {from_version}")
    if to_dir is None:
        raise HTTPException(status_code=404, detail=f"Disease graph version not found: {to_version}")

    baseline = load_version_data(from_dir)
    current = load_disease_graph_data(_disease_graph_dir) if to_dir.resolve() == Path(_disease_graph_dir).resolve() else load_version_data(to_dir)
    return compute_version_diff(current, baseline)


# ---------------------------------------------------------------------------
# Disease Explorer — Graph Neighborhood & Provenance
# ---------------------------------------------------------------------------

@app.get("/disease-id-qa/api/graph/neighbors")
def disease_id_qa_graph_neighbors(pxref: str = ""):
    """Return Cytoscape elements for concepts sharing xrefs with the given concept."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if not pxref:
        raise HTTPException(status_code=400, detail="pxref parameter required.")
    data = load_disease_graph_data(_disease_graph_dir)
    neighbors = find_concept_neighbors(data, pxref, max_neighbors=10)
    if not neighbors:
        return {"elements": []}
    # Build graph payload for the neighbors (include the original concept too)
    all_ids = [pxref] + neighbors
    target_data = None
    if _target_graph_dir:
        try:
            target_data = load_target_graph_data(_target_graph_dir)
        except Exception:
            pass
    return build_disease_graph_payload(data, all_ids, target_data=target_data)


@app.get("/disease-id-qa/api/graph/hierarchy")
def disease_id_qa_graph_hierarchy(pxref: str = ""):
    """Return Cytoscape elements for immediate MONDO parent/child context."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if not pxref:
        raise HTTPException(status_code=400, detail="pxref parameter required.")
    data = load_disease_graph_data(_disease_graph_dir)
    return build_hierarchy_graph_payload(data, pxref)


@app.get("/disease-id-qa/api/graph/clinical-descendants")
def disease_id_qa_graph_clinical_descendants(pxref: str = "", limit: int = 80):
    """Return Cytoscape elements for SNOMED/ICD descendant candidate mappings."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    if not pxref:
        raise HTTPException(status_code=400, detail="pxref parameter required.")
    limit = max(1, min(limit, 500))
    data = load_disease_graph_data(_disease_graph_dir)
    return build_clinical_descendant_graph_payload(data, pxref, limit=limit)


@app.get("/disease-id-qa/api/provenance/{pxref:path}")
def disease_id_qa_provenance(pxref: str):
    """Return the provenance chain for a concept."""
    if not _disease_graph_dir:
        raise HTTPException(status_code=500, detail="No --disease-graph-dir configured.")
    data = load_disease_graph_data(_disease_graph_dir)
    chain = build_provenance_chain(data, pxref)
    if not chain:
        raise HTTPException(status_code=404, detail=f"Concept not found: {pxref}")
    return {"primary_xref": pxref, "steps": chain}


@app.get("/db/{db_name}", response_class=HTMLResponse)
def dashboard(request: Request, db_name: str):
    db = get_db(db_name)
    collections = _get_dashboard_collection_shell(db)

    # Edge definitions for schema summary
    edge_defs = _edge_definitions_with_endpoint_pairs(db)

    # ETL metadata
    etl_meta = None
    if db.has_collection("metadata_store"):
        try:
            cursor = db.aql.execute('RETURN DOCUMENT("metadata_store", "etl_metadata").value')
            results = list(cursor)
            if results and results[0]:
                etl_meta = results[0]
        except Exception:
            pass

    graph_views = _get_graph_views(db)
    build_inputs = extract_build_inputs(etl_meta)
    return templates.TemplateResponse(request, "dashboard.html", {
        "request": request,
        "db_name": db_name,
        "collections": collections,
        "edge_defs": edge_defs,
        "etl_meta": etl_meta,
        "build_inputs": build_inputs,
        "graph_views": graph_views,
        "doc_count": None,
        "edge_count": None,
        "counts_url": _app_path(f"/db/{db_name}/collections", request),
    })


@app.get("/db/{db_name}/collections", response_class=HTMLResponse)
def dashboard_collections(request: Request, db_name: str):
    db = get_db(db_name)
    collections = _get_dashboard_collection_summaries(db)
    return templates.TemplateResponse(request, "dashboard_collections.html", {
        "request": request,
        "db_name": db_name,
        "collections": collections,
        "doc_count": sum(c["count"] for c in collections if c["type"] == "document"),
        "edge_count": sum(c["count"] for c in collections if c["type"] == "edge"),
    })


@app.get("/db/{db_name}/build-status", response_class=HTMLResponse)
async def build_status_page(request: Request, db_name: str):
    db = get_db(db_name)
    build_status = _get_build_status(db)
    return templates.TemplateResponse(request, "build_status.html", {
        "request": request,
        "db_name": db_name,
        "build_status": build_status,
    })


@app.get("/db/{db_name}/view/{view_id}/preview", response_class=HTMLResponse)
async def preview_graph_view(request: Request, db_name: str, view_id: str, limit: int = 50):
    db = get_db(db_name)
    graph_views = _get_graph_views(db)
    graph_view = graph_views.get(view_id)

    if not graph_view:
        return HTMLResponse(f"Graph view '{view_id}' not found.", status_code=404)

    if graph_view.get("query_language") != "aql":
        return HTMLResponse("Only AQL graph views are supported.", status_code=400)

    if graph_view.get("output_format") != "jsonl":
        return HTMLResponse("Only JSONL graph views can be previewed.", status_code=400)

    query = graph_view.get("query")
    if not query:
        return HTMLResponse("Graph view is missing query metadata.", status_code=400)

    preview_limit = max(1, min(limit, 200))
    preview_query = f"""
    LET graph_view_rows = (
      {query}
    )
    RETURN {{
      total_count: LENGTH(graph_view_rows),
      rows: (
        FOR row IN graph_view_rows
          LIMIT @limit
          RETURN row
      )
    }}
    """

    try:
        result = next(iter(db.aql.execute(preview_query, bind_vars={"limit": preview_limit}, max_runtime=60)), None)
    except Exception as exc:
        return templates.TemplateResponse(request, "graph_view_preview.html", {
            "request": request,
            "db_name": db_name,
            "view_id": view_id,
            "graph_view": graph_view,
            "rows": [],
            "total_count": None,
            "preview_limit": preview_limit,
            "error": str(exc),
        })

    return templates.TemplateResponse(request, "graph_view_preview.html", {
        "request": request,
        "db_name": db_name,
        "view_id": view_id,
        "graph_view": graph_view,
        "rows": (result or {}).get("rows") or [],
        "total_count": (result or {}).get("total_count"),
        "preview_limit": preview_limit,
        "error": None,
    })


@app.get("/db/{db_name}/view/{view_id}")
async def execute_graph_view(db_name: str, view_id: str):
    db = get_db(db_name)
    graph_views = _get_graph_views(db)
    graph_view = graph_views.get(view_id)

    if not graph_view:
        return HTMLResponse(f"Graph view '{view_id}' not found.", status_code=404)

    if graph_view.get("query_language") != "aql":
        return HTMLResponse("Only AQL graph views are supported.", status_code=400)

    output_format = graph_view.get("output_format")
    if output_format not in {"csv", "jsonl"}:
        return HTMLResponse("Only CSV and JSONL graph views are supported.", status_code=400)

    query = graph_view.get("query")
    columns = graph_view.get("columns") or []
    if not query:
        return HTMLResponse("Graph view is missing query metadata.", status_code=400)
    if output_format == "csv" and not columns:
        return HTMLResponse("CSV graph view is missing columns metadata.", status_code=400)

    try:
        rows = list(db.aql.execute(query, batch_size=5000, max_runtime=600))
    except Exception as exc:
        return HTMLResponse(f"Failed to execute graph view '{view_id}': {exc}", status_code=500)

    if output_format == "jsonl":
        def generate_jsonl():
            for row in rows:
                yield json.dumps(row, default=str) + "\n"

        return StreamingResponse(
            generate_jsonl(),
            media_type="application/x-ndjson",
            headers={"Content-Disposition": f"attachment; filename={view_id}.jsonl"},
        )

    output = io.StringIO()
    writer = csv.DictWriter(output, fieldnames=columns, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        normalized_row = {}
        for column in columns:
            value = row.get(column) if isinstance(row, dict) else None
            if isinstance(value, (dict, list)):
                normalized_row[column] = json.dumps(value, default=str)
            elif value is None:
                normalized_row[column] = ""
            else:
                normalized_row[column] = value
        writer.writerow(normalized_row)
    output.seek(0)

    return StreamingResponse(
        output,
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={view_id}.csv"},
    )


@app.get("/db/{db_name}/collection/{coll_name}", response_class=HTMLResponse)
def collection_browser(request: Request, db_name: str, coll_name: str,
                       page: int = 1, page_size: int = 25):
    db = get_db(db_name)
    skip = (page - 1) * page_size
    facet_metadata = _get_collection_facet_metadata(db, coll_name)
    search_metadata = _get_collection_search_metadata(db, coll_name)
    category_fields = sorted(facet_metadata.get("category_fields") or [])
    search_fields = list(search_metadata.get("text_fields") or [])
    active_filters = _parse_collection_facet_filters(request, category_fields)
    search_term = _parse_collection_search_term(request)
    filter_bind_vars = {}
    filter_clause = _build_collection_constraints(
        active_filters=active_filters,
        search_fields=search_fields,
        search_term=search_term,
        bind_vars=filter_bind_vars,
    )
    htmx = request.headers.get("HX-Request") == "true"

    if not htmx:
        coll = db.collection(coll_name)
        is_edge = coll.properties().get("type") in ("edge", 3)
        active_filter_summary = _build_active_filter_summary(
            db_name=db_name,
            coll_name=coll_name,
            page_size=page_size,
            active_filters=active_filters,
            search_term=search_term,
        )
        stats_url = _build_collection_stats_url(db_name, coll_name, active_filters, search_term=search_term)
        facets_url = _build_collection_facets_url(
            db_name=db_name,
            coll_name=coll_name,
            facet_filters=active_filters,
            search_term=search_term,
            page_size=page_size,
        )
        facet_loaders = _build_collection_facet_loaders(
            db_name=db_name,
            coll_name=coll_name,
            page_size=page_size,
            category_fields=category_fields,
            active_filters=active_filters,
            search_term=search_term,
        )
        rows_url = _build_collection_url(
            db_name=db_name,
            coll_name=coll_name,
            page=page,
            page_size=page_size,
            facet_filters=active_filters,
            search_term=search_term,
        )
        download_url = _build_collection_download_url(
            db_name=db_name,
            coll_name=coll_name,
            page=page,
            page_size=page_size,
            facet_filters=active_filters,
            search_term=search_term,
        )
        clear_search_url = _build_collection_url(
            db_name=db_name,
            coll_name=coll_name,
            page=1,
            page_size=page_size,
            facet_filters=active_filters,
            search_term="",
        )
        return templates.TemplateResponse(request, "collection.html", {
            "request": request,
            "db_name": db_name,
            "coll_name": coll_name,
            "docs": [],
            "columns": [],
            "is_edge": is_edge,
            "total": None,
            "page": page,
            "page_size": page_size,
            "total_pages": 1,
            "has_facets": bool(category_fields),
            "facet_panels": [],
            "facet_loaders": facet_loaders,
            "active_filters": active_filters,
            "active_filter_summary": active_filter_summary,
            "search_term": search_term,
            "search_fields": search_fields,
            "stats_url": stats_url,
            "facets_url": facets_url,
            "download_url": download_url,
            "clear_search_url": clear_search_url,
            "rows_url": rows_url,
            "loading_table": True,
            "prev_url": None,
            "next_url": None,
        })

    # Use AQL for both count and list so they always agree
    count_query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            COLLECT WITH COUNT INTO c
            RETURN c
    """
    count_cursor = db.aql.execute(count_query, bind_vars=filter_bind_vars)
    total = list(count_cursor)[0]
    total_pages = max(1, (total + page_size - 1) // page_size)
    page = min(page, total_pages)
    skip = (page - 1) * page_size
    list_bind_vars = {**filter_bind_vars, "skip": skip, "top": page_size}

    coll = db.collection(coll_name)
    is_edge = coll.properties().get("type") in ("edge", 3)
    preview_fields = _get_collection_preview_fields(
        is_edge=is_edge,
        facet_metadata=facet_metadata,
        search_metadata=search_metadata,
    )
    list_bind_vars["return_fields"] = preview_fields

    # Fetch documents
    query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            SORT doc._key ASC
            LIMIT @skip, @top
            RETURN KEEP(doc, @return_fields)
    """
    cursor = db.aql.execute(query, bind_vars=list_bind_vars)
    docs = list(cursor)
    # Discover columns from this page of results
    columns = _discover_columns(docs, is_edge)
    active_filter_summary = _build_active_filter_summary(
        db_name=db_name,
        coll_name=coll_name,
        page_size=page_size,
        active_filters=active_filters,
        search_term=search_term,
    )
    stats_url = _build_collection_stats_url(db_name, coll_name, active_filters, search_term=search_term)
    facets_url = _build_collection_facets_url(
        db_name=db_name,
        coll_name=coll_name,
        facet_filters=active_filters,
        search_term=search_term,
        page_size=page_size,
    )
    facet_loaders = _build_collection_facet_loaders(
        db_name=db_name,
        coll_name=coll_name,
        page_size=page_size,
        category_fields=category_fields,
        active_filters=active_filters,
        search_term=search_term,
    )
    download_url = _build_collection_download_url(db_name, coll_name, page, page_size, active_filters, search_term=search_term)
    clear_search_url = _build_collection_url(
        db_name=db_name,
        coll_name=coll_name,
        page=1,
        page_size=page_size,
        facet_filters=active_filters,
        search_term="",
    )
    prev_url = None
    next_url = None
    if page > 1:
        prev_url = _build_collection_url(db_name, coll_name, page - 1, page_size, active_filters, search_term=search_term)
    if page < total_pages:
        next_url = _build_collection_url(db_name, coll_name, page + 1, page_size, active_filters, search_term=search_term)

    return templates.TemplateResponse(request, "collection_rows.html", {
        "request": request,
        "db_name": db_name,
        "coll_name": coll_name,
        "docs": docs,
        "columns": columns,
        "is_edge": is_edge,
        "total": total,
        "page": page,
        "page_size": page_size,
        "total_pages": total_pages,
        "has_facets": bool(category_fields),
        "facet_panels": [],
        "facet_loaders": facet_loaders,
        "active_filters": active_filters,
        "active_filter_summary": active_filter_summary,
        "search_term": search_term,
        "search_fields": search_fields,
        "stats_url": stats_url,
        "facets_url": facets_url,
        "download_url": download_url,
        "clear_search_url": clear_search_url,
        "rows_url": "",
        "loading_table": False,
        "prev_url": prev_url,
        "next_url": next_url,
    })


@app.get("/db/{db_name}/collection/{coll_name}/download.csv")
async def collection_download(request: Request, db_name: str, coll_name: str,
                              page: int = 1, page_size: int = 25):
    db = get_db(db_name)
    facet_metadata = _get_collection_facet_metadata(db, coll_name)
    search_metadata = _get_collection_search_metadata(db, coll_name)
    category_fields = sorted(facet_metadata.get("category_fields") or [])
    search_fields = list(search_metadata.get("text_fields") or [])
    active_filters = _parse_collection_facet_filters(request, category_fields)
    search_term = _parse_collection_search_term(request)
    filter_bind_vars = {}
    filter_clause = _build_collection_constraints(
        active_filters=active_filters,
        search_fields=search_fields,
        search_term=search_term,
        bind_vars=filter_bind_vars,
    )

    coll = db.collection(coll_name)
    is_edge = coll.properties().get("type") in ("edge", 3)
    preview_fields = _get_collection_preview_fields(
        is_edge=is_edge,
        facet_metadata=facet_metadata,
        search_metadata=search_metadata,
    )

    # Match the export columns to the current list-page view by rediscovering
    # columns from the currently visible page, then export all filtered rows.
    page = max(page, 1)
    page_size = max(page_size, 1)
    skip = (page - 1) * page_size
    preview_bind_vars = {**filter_bind_vars, "skip": skip, "top": page_size, "return_fields": preview_fields}
    preview_query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            SORT doc._key ASC
            LIMIT @skip, @top
            RETURN KEEP(doc, @return_fields)
    """
    preview_docs = list(db.aql.execute(preview_query, bind_vars=preview_bind_vars))
    columns = _discover_columns(preview_docs, is_edge)
    if not columns:
        columns = ["_key"]

    export_query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            SORT doc._key ASC
            RETURN doc
    """
    cursor = db.aql.execute(export_query, bind_vars=filter_bind_vars)

    def generate_csv():
        buffer = io.StringIO()
        writer = csv.DictWriter(buffer, fieldnames=columns, extrasaction="ignore")
        writer.writeheader()
        yield buffer.getvalue()
        buffer.seek(0)
        buffer.truncate(0)

        for doc in cursor:
            row = {column: _normalize_csv_value(doc.get(column)) for column in columns}
            writer.writerow(row)
            yield buffer.getvalue()
            buffer.seek(0)
            buffer.truncate(0)

    filename = f"{db_name}_{coll_name}.csv"
    return StreamingResponse(
        generate_csv(),
        media_type="text/csv",
        headers={"Content-Disposition": f"attachment; filename={filename}"},
    )


@app.get("/db/{db_name}/collection/{coll_name}/facets", response_class=HTMLResponse)
def collection_facets(request: Request, db_name: str, coll_name: str, page_size: int = 25):
    """Facet panels for a collection (loaded after the main table)."""
    db = get_db(db_name)
    facet_metadata = _get_collection_facet_metadata(db, coll_name)
    search_metadata = _get_collection_search_metadata(db, coll_name)
    category_fields = sorted(facet_metadata.get("category_fields") or [])
    search_fields = list(search_metadata.get("text_fields") or [])
    active_filters = _parse_collection_facet_filters(request, category_fields)
    search_term = _parse_collection_search_term(request)
    facet_panels = _build_collection_facet_panels(
        db=db,
        db_name=db_name,
        coll_name=coll_name,
        page_size=page_size,
        category_fields=category_fields,
        active_filters=active_filters,
        search_fields=search_fields,
        search_term=search_term,
    )
    return templates.TemplateResponse(request, "collection_facets.html", {
        "request": request,
        "facet_panels": facet_panels,
    })


@app.get("/db/{db_name}/collection/{coll_name}/facet/{field}", response_class=HTMLResponse)
def collection_facet(request: Request, db_name: str, coll_name: str, field: str, page_size: int = 25):
    """Single facet panel for a collection (loaded independently via HTMX)."""
    db = get_db(db_name)
    facet_metadata = _get_collection_facet_metadata(db, coll_name)
    search_metadata = _get_collection_search_metadata(db, coll_name)
    category_fields = sorted(facet_metadata.get("category_fields") or [])
    if field not in category_fields:
        return templates.TemplateResponse(request, "collection_facet_panel.html", {
            "request": request,
            "panel": {
                "field": field,
                "title": field.replace("_", " "),
                "values": [],
                "selected_values": [],
                "clear_href": "#",
            },
        })
    search_fields = list(search_metadata.get("text_fields") or [])
    active_filters = _parse_collection_facet_filters(request, category_fields)
    search_term = _parse_collection_search_term(request)
    panel = _build_collection_facet_panel(
        db=db,
        db_name=db_name,
        coll_name=coll_name,
        page_size=page_size,
        field=field,
        active_filters=active_filters,
        search_fields=search_fields,
        search_term=search_term,
    )
    return templates.TemplateResponse(request, "collection_facet_panel.html", {
        "request": request,
        "panel": panel,
    })


@app.get("/db/{db_name}/collection/{coll_name}/stats", response_class=HTMLResponse)
def collection_stats(request: Request, db_name: str, coll_name: str):
    """Field coverage stats for a collection (loaded via HTMX)."""
    db = get_db(db_name)
    facet_metadata = _get_collection_facet_metadata(db, coll_name)
    search_metadata = _get_collection_search_metadata(db, coll_name)
    category_fields = sorted(facet_metadata.get("category_fields") or [])
    search_fields = list(search_metadata.get("text_fields") or [])
    active_filters = _parse_collection_facet_filters(request, category_fields)
    search_term = _parse_collection_search_term(request)
    filter_bind_vars = {}
    filter_clause = _build_collection_constraints(
        active_filters=active_filters,
        search_fields=search_fields,
        search_term=search_term,
        bind_vars=filter_bind_vars,
    )

    count_query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            COLLECT WITH COUNT INTO c
            RETURN c
    """
    count_results = list(db.aql.execute(count_query, bind_vars=filter_bind_vars))
    total = count_results[0] if count_results else 0
    sample_size = min(total, 500)
    sample_bind_vars = {**filter_bind_vars, "sample_size": sample_size}

    # Sample documents to discover fields
    query = f"""
        FOR doc IN `{coll_name}`
            {f"FILTER {filter_clause}" if filter_clause else ""}
            SORT RAND()
            LIMIT @sample_size
            RETURN ATTRIBUTES(doc)
    """
    cursor = db.aql.execute(query, bind_vars=sample_bind_vars)
    all_attrs = list(cursor)

    field_counts = {}
    for attrs in all_attrs:
        for attr in attrs:
            if attr.startswith("_"):
                continue
            field_counts[attr] = field_counts.get(attr, 0) + 1

    stats = []
    for field, count in sorted(field_counts.items(), key=lambda x: -x[1]):
        pct = round(100 * count / sample_size, 1) if sample_size > 0 else 0
        stats.append({"field": field, "count": count, "pct": pct, "sample_size": sample_size})

    return templates.TemplateResponse(request, "stats_partial.html", {
        "request": request,
        "stats": stats,
        "sample_size": sample_size,
        "total": total,
    })


@app.get("/db/{db_name}/collection/{coll_name}/doc/{doc_key:path}/parquet-stats", response_class=HTMLResponse)
async def parquet_stats(request: Request, db_name: str, coll_name: str, doc_key: str):
    """Load parquet file stats for a Dataset document (called via HTMX)."""
    db = get_db(db_name)
    coll = db.collection(coll_name)
    error = None
    stats = None

    try:
        doc = coll.get(doc_key)
        file_ref = doc.get("file_reference") if doc else None
        if not file_ref:
            error = "No file_reference found on this document."
        else:
            import pyarrow.parquet as pq
            buf, content_length = _get_parquet_buffer(file_ref)
            if buf is None:
                error = f"Could not fetch parquet file: {file_ref}"
            else:
                pf = pq.ParquetFile(buf)
                metadata = pf.metadata

                # Read into pandas for descriptive stats
                table = pf.read()
                df = table.to_pandas()

                col_stats = []
                for col_name in df.columns:
                    col = df[col_name]
                    info = {"name": col_name, "dtype": str(col.dtype), "non_null": int(col.count())}
                    if col.dtype.kind in ("f", "i", "u"):  # numeric
                        info["min"] = f"{col.min():.4g}"
                        info["max"] = f"{col.max():.4g}"
                        info["mean"] = f"{col.mean():.4g}"
                        info["std"] = f"{col.std():.4g}"
                    else:
                        info["unique"] = int(col.nunique())
                    col_stats.append(info)

                stats = {
                    "file_path": file_ref,
                    "file_size_mb": round(content_length / (1024 * 1024), 2),
                    "num_rows": metadata.num_rows,
                    "num_columns": metadata.num_columns,
                    "num_row_groups": metadata.num_row_groups,
                    "index_name": df.index.name,
                    "index_count": len(df.index),
                    "columns": col_stats,
                    "head": df.head(5).to_html(classes="parquet-table", border=0),
                }
    except Exception as e:
        error = str(e)

    return templates.TemplateResponse(request, "parquet_stats.html", {
        "request": request,
        "stats": stats,
        "error": error,
    })


@app.get("/api/record-curation-cart")
def record_curation_cart(request: Request, curator: str = "", curator_name: str = ""):
    curator_id, display_name = _curator_identity(request, {
        "curator": curator,
        "curator_name": curator_name,
    })
    try:
        return _combined_curation_cart(
            _curation_cart_storage(), curator_id, display_name
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.post("/api/record-curation-cart/items")
async def add_record_curation_cart_item(request: Request):
    payload = await request.json()
    if not isinstance(payload, dict):
        raise HTTPException(status_code=400, detail="Request body must be a JSON object.")
    curator_id, curator_name = _curator_identity(request, payload)
    db_name = str(
        payload.get("database_name") or payload.get("curation_set") or ""
    ).strip()
    coll_name = str(payload.get("model_type") or "").strip()
    doc_key = str(payload.get("doc_key") or "").strip()
    try:
        db = get_db(db_name)
        if not db.has_collection(coll_name):
            raise ValueError(f"Collection does not exist: {coll_name}")
        coll = db.collection(coll_name)
        if coll.properties().get("type") in ("edge", 3):
            raise ValueError("Edge property curations are not supported yet")
        document = coll.get(doc_key)
        if not document or not document.get("id"):
            raise ValueError("The curation target record no longer exists")
        target_curation_set = record_property_curation_set_for_model(
            coll_name, db_name
        )
        schema_fields = (_get_collection_schema_entry(db, coll_name).get("fields") or {})
        if not schema_fields:
            raise ValueError("This collection does not publish a curatable schema")
        raw_decisions = payload.get("decisions")
        if not isinstance(raw_decisions, list) or not raw_decisions:
            raise ValueError("Choose at least one property to curate")
        curation_type, operation = _build_record_property_operation(
            document,
            curation_set=target_curation_set,
            model_type=coll_name,
            schema_fields=schema_fields,
            raw_decisions=raw_decisions,
            note=str(payload.get("note") or ""),
        )
        storage = _curation_cart_storage()
        existing_cart = await run_in_threadpool(
            load_cart, storage, curation_type, curator_id, curator_name
        )
        existing_sets = {
            (item.get("target") or {}).get("curation_set")
            for item in existing_cart.get("operations") or []
        }
        existing_sets.discard(None)
        if existing_sets and existing_sets != {target_curation_set}:
            raise ValueError(
                "Publish or clear the current review before curating a different graph"
            )
        await run_in_threadpool(
            add_cart_operation,
            storage,
            curation_type,
            curator_id,
            curator_name,
            operation,
            bool(payload.get("replace_target", False)),
        )
        return await run_in_threadpool(
            _combined_curation_cart, storage, curator_id, curator_name
        )
    except ValueError as exc:
        raise HTTPException(status_code=400, detail=str(exc)) from exc


@app.delete("/api/record-curation-cart/items/{operation_id}")
async def remove_record_curation_cart_item(operation_id: str, request: Request):
    return await ramp_id_qa_remove_curation_cart_item(operation_id, request)


@app.post("/api/record-curation-cart/publish")
async def publish_record_curation_cart(request: Request):
    return await ramp_id_qa_publish_curation_cart(request)


@app.get("/db/{db_name}/collection/{coll_name}/doc/{doc_key:path}", response_class=HTMLResponse)
async def document_detail(request: Request, db_name: str, coll_name: str, doc_key: str):
    db = get_db(db_name)
    coll = db.collection(coll_name)

    try:
        doc = coll.get(doc_key)
    except Exception:
        doc = None

    is_edge = coll.properties().get("type") in ("edge", 3) if doc else False

    # Parse provenance/sources for nicer display
    sources = []
    if doc and "sources" in doc and doc["sources"]:
        for src in doc["sources"]:
            if isinstance(src, str) and "\t" in src:
                parts = src.split("\t")
                sources.append({
                    "name": parts[0] if len(parts) > 0 else "",
                    "version": parts[1] if len(parts) > 1 else "",
                    "version_date": parts[2] if len(parts) > 2 else "",
                    "download_date": parts[3] if len(parts) > 3 else "",
                })
            else:
                sources.append({"name": str(src), "version": "", "version_date": "", "download_date": ""})

    # Categorize fields for display
    scalar_fields, list_fields, nested_fields = _categorize_document_fields(doc)

    # Find connected nodes via graph traversal (vertex nodes only)
    outgoing_linked_groups = []
    incoming_linked_groups = []
    linked_aql = ""
    if doc and not is_edge and db.has_graph("graph"):
        try:
            _SAMPLE_LIMIT = 20
            linked_aql = (
                f"LET outgoing_counts = (\n"
                f"    FOR v, e IN 1..1 OUTBOUND '{doc['_id']}' GRAPH 'graph'\n"
                f"    COLLECT coll = SPLIT(v._id, '/')[0], edge_coll = SPLIT(e._id, '/')[0] WITH COUNT INTO cnt\n"
                f"    RETURN {{coll, edge_coll, cnt}}\n"
                f")\n"
                f"LET incoming_counts = (\n"
                f"    FOR v, e IN 1..1 INBOUND '{doc['_id']}' GRAPH 'graph'\n"
                f"    COLLECT coll = SPLIT(v._id, '/')[0], edge_coll = SPLIT(e._id, '/')[0] WITH COUNT INTO cnt\n"
                f"    RETURN {{coll, edge_coll, cnt}}\n"
                f")\n"
                f"RETURN {{\n"
                f"    outgoing: (\n"
                f"        FOR c IN outgoing_counts\n"
                f"            SORT c.cnt DESC, c.edge_coll\n"
                f"            LET samples = (\n"
                f"                FOR v2, e2 IN 1..1 OUTBOUND '{doc['_id']}' GRAPH 'graph'\n"
                f"                FILTER SPLIT(v2._id, '/')[0] == c.coll\n"
                f"                FILTER SPLIT(e2._id, '/')[0] == c.edge_coll\n"
                f"                LIMIT {_SAMPLE_LIMIT}\n"
                f"                RETURN {{key: v2._key, label: v2.name || v2.symbol || v2._key}}\n"
                f"            )\n"
                f"            RETURN {{collection: c.coll, edge_type: c.edge_coll, count: c.cnt, nodes: samples}}\n"
                f"    ),\n"
                f"    incoming: (\n"
                f"        FOR c IN incoming_counts\n"
                f"            SORT c.cnt DESC, c.edge_coll\n"
                f"            LET samples = (\n"
                f"                FOR v2, e2 IN 1..1 INBOUND '{doc['_id']}' GRAPH 'graph'\n"
                f"                FILTER SPLIT(v2._id, '/')[0] == c.coll\n"
                f"                FILTER SPLIT(e2._id, '/')[0] == c.edge_coll\n"
                f"                LIMIT {_SAMPLE_LIMIT}\n"
                f"                RETURN {{key: v2._key, label: v2.name || v2.symbol || v2._key}}\n"
                f"            )\n"
                f"            RETURN {{collection: c.coll, edge_type: c.edge_coll, count: c.cnt, nodes: samples}}\n"
                f"    )\n"
                f"}}"
            )
            cursor = db.aql.execute(
                """
                LET outgoing_counts = (
                    FOR v, e IN 1..1 OUTBOUND @node_id GRAPH 'graph'
                    COLLECT coll = SPLIT(v._id, '/')[0], edge_coll = SPLIT(e._id, '/')[0] WITH COUNT INTO cnt
                    RETURN {coll, edge_coll, cnt}
                )
                LET incoming_counts = (
                    FOR v, e IN 1..1 INBOUND @node_id GRAPH 'graph'
                    COLLECT coll = SPLIT(v._id, '/')[0], edge_coll = SPLIT(e._id, '/')[0] WITH COUNT INTO cnt
                    RETURN {coll, edge_coll, cnt}
                )
                RETURN {
                    outgoing: (
                        FOR c IN outgoing_counts
                            SORT c.cnt DESC, c.edge_coll
                            LET samples = (
                                FOR v2, e2 IN 1..1 OUTBOUND @node_id GRAPH 'graph'
                                FILTER SPLIT(v2._id, '/')[0] == c.coll
                                FILTER SPLIT(e2._id, '/')[0] == c.edge_coll
                                LIMIT @sample_limit
                                RETURN {key: v2._key, label: v2.name || v2.symbol || v2._key}
                            )
                            RETURN {collection: c.coll, edge_type: c.edge_coll, count: c.cnt, nodes: samples}
                    ),
                    incoming: (
                        FOR c IN incoming_counts
                            SORT c.cnt DESC, c.edge_coll
                            LET samples = (
                                FOR v2, e2 IN 1..1 INBOUND @node_id GRAPH 'graph'
                                FILTER SPLIT(v2._id, '/')[0] == c.coll
                                FILTER SPLIT(e2._id, '/')[0] == c.edge_coll
                                LIMIT @sample_limit
                                RETURN {key: v2._key, label: v2.name || v2.symbol || v2._key}
                            )
                            RETURN {collection: c.coll, edge_type: c.edge_coll, count: c.cnt, nodes: samples}
                    )
                }
                """,
                bind_vars={"node_id": doc["_id"], "sample_limit": _SAMPLE_LIMIT},
                max_runtime=15,
            )
            linked_result = list(cursor)
            if linked_result:
                outgoing_linked_groups = linked_result[0].get("outgoing", [])
                incoming_linked_groups = linked_result[0].get("incoming", [])
        except Exception:
            pass

    template_name = _get_document_template(db_name, coll_name)
    record_curation = None
    if doc:
        target_curation_set = record_property_curation_set_for_model(
            coll_name, db_name
        )
        record_curation = {
            "available": False,
            "curation_set": target_curation_set,
            "database_name": db_name,
            "model_type": coll_name,
            "target_id": doc.get("id"),
            "doc_key": doc_key,
            "fields": [],
            "application_mode": (
                "evidence_overlay"
                if db_name == "metabolite_harmonization"
                else "materialized"
            ),
            "published_fields": [],
            "published_summary": None,
            "curation_type": record_property_type_for_model(coll_name),
        }
        if is_edge:
            record_curation["error"] = "Edge property curation is not supported yet."
        elif not doc.get("id"):
            record_curation["error"] = "This record has no stable ID to curate."
        elif not record_curation["curation_type"]:
            record_curation["error"] = (
                f"No record-property curation stream is registered for {coll_name}."
            )
        else:
            schema_entry = _get_collection_schema_entry(db, coll_name)
            if not schema_entry.get("fields"):
                record_curation["error"] = "No editable schema fields are published for this collection."
            else:
                try:
                    snapshot = await run_in_threadpool(
                        _load_record_property_snapshot, coll_name
                    )
                    curation_document = dict(doc)
                    curation_document["_curation_set"] = target_curation_set
                    curation_document["_curation_model_type"] = coll_name
                    fields = _record_property_fields(
                        curation_document, schema_entry, snapshot
                    )
                    record_curation["fields"] = fields
                    record_curation["published_fields"] = [
                        field for field in fields
                        if field.get("has_published_override")
                    ]
                    if record_curation["published_fields"]:
                        record_curation["published_summary"] = (
                            _record_curation_summary(
                                record_curation["published_fields"]
                            )
                        )
                    record_curation["available"] = bool(fields)
                    if not fields:
                        record_curation["error"] = "This record has no editable field values."
                except Exception as exc:
                    record_curation["error"] = f"Curation is temporarily unavailable: {exc}"

    if (
        _template_supports_inline_curations(template_name)
        and record_curation
        and record_curation["published_fields"]
    ):
        display_doc = _document_with_inline_curations(
            doc, record_curation["published_fields"]
        )
        scalar_fields, list_fields, nested_fields = _categorize_document_fields(
            display_doc
        )

    context = {
        "request": request,
        "db_name": db_name,
        "coll_name": coll_name,
        "doc": doc,
        "doc_key": doc_key,
        "is_edge": is_edge,
        "sources": sources,
        "scalar_fields": scalar_fields,
        "list_fields": list_fields,
        "nested_fields": nested_fields,
        "outgoing_linked_groups": outgoing_linked_groups,
        "incoming_linked_groups": incoming_linked_groups,
        "linked_aql": linked_aql,
        "this_node_label": f"This {coll_name}",
        "record_curation": record_curation,
    }
    if template_name == "cure_case_report_document.html":
        context.update(_get_adjacent_collection_docs(db, coll_name, doc))
        context.update(_build_cure_case_report_context(db_name, db, doc))
    return templates.TemplateResponse(request, template_name, context)


@app.get("/db/{db_name}/schema", response_class=HTMLResponse)
async def schema_view(request: Request, db_name: str):
    db = get_db(db_name)
    edge_defs = _edge_definitions_with_endpoint_pairs(db)

    # Build mermaid diagram
    mermaid_lines = ["graph LR"]
    nodes_seen = set()
    for ed in edge_defs:
        edge_name = ed["edge_collection"]
        for endpoint_pair in ed.get("endpoint_pairs") or []:
            frm = endpoint_pair["from_collection"]
            to = endpoint_pair["to_collection"]
            safe_from = frm.replace(":", "_").replace(" ", "_")
            safe_to = to.replace(":", "_").replace(" ", "_")
            if safe_from not in nodes_seen:
                mermaid_lines.append(f'    {safe_from}["{frm}"]')
                nodes_seen.add(safe_from)
            if safe_to not in nodes_seen:
                mermaid_lines.append(f'    {safe_to}["{to}"]')
                nodes_seen.add(safe_to)
            mermaid_lines.append(f"    {safe_from} -->|{edge_name}| {safe_to}")

    mermaid_text = "\n".join(mermaid_lines)

    return templates.TemplateResponse(request, "schema.html", {
        "request": request,
        "db_name": db_name,
        "edge_defs": edge_defs,
        "mermaid_text": mermaid_text,
    })


@app.get("/db/{db_name}/aql", response_class=HTMLResponse)
async def aql_page(request: Request, db_name: str):
    return templates.TemplateResponse(request, "aql.html", {
        "request": request,
        "db_name": db_name,
        "results": None,
        "query": "",
        "error": None,
        "columns": [],
    })


@app.post("/db/{db_name}/aql", response_class=HTMLResponse)
async def aql_execute(request: Request, db_name: str, query: str = Form(...)):
    results = None
    error = None
    columns = []
    try:
        db = get_db(db_name)
        cursor = db.aql.execute(query, max_runtime=60)
        results = list(cursor)
        # Auto-detect columns from results
        if results and isinstance(results[0], dict):
            col_set = set()
            for row in results[:50]:
                col_set.update(row.keys())
            columns = sorted(col_set)
    except Exception as e:
        error = str(e)

    htmx = request.headers.get("HX-Request") == "true"
    template = "aql_results.html" if htmx else "aql.html"

    return templates.TemplateResponse(request, template, {
        "request": request,
        "db_name": db_name,
        "results": results,
        "query": query,
        "error": error,
        "columns": columns,
    })


# ── MySQL Routes ─────────────────────────────────────────────────────────────

def _mysql_template_context(source_id: str, db_name: str) -> dict:
    source = _get_mysql_source(source_id)
    return {
        "mysql_source_id": source_id,
        "mysql_source_label": source["label"],
        "db_name": db_name,
    }


@app.get("/mysql/{source_id}/{db_name}", response_class=HTMLResponse)
async def mysql_dashboard(request: Request, source_id: str, db_name: str):
    engine = get_mysql_db_engine(db_name, source_id=source_id)
    schema_meta = get_mysql_inspector(db_name, source_id=source_id)

    # SHOW TABLE STATUS returns approximate row counts — much faster than COUNT(*) for InnoDB
    with engine.connect() as conn:
        status_rows = conn.execute(text("SHOW TABLE STATUS")).mappings().all()
    row_counts = {r["Name"]: (r["Rows"] or 0) for r in status_rows}

    tables = []
    for table_name, meta in schema_meta.items():
        tables.append({
            "name": table_name,
            "count": row_counts.get(table_name, 0),
            "fk_count": len(meta["fks"]),
        })
    tables.sort(key=lambda t: t["name"])

    fk_defs = []
    for table_name, meta in schema_meta.items():
        for fk in meta["fks"]:
            fk_defs.append({
                "from_table": table_name,
                "from_columns": ", ".join(fk["constrained_columns"]),
                "to_table": fk["referred_table"],
                "to_columns": ", ".join(fk["referred_columns"]),
            })

    return templates.TemplateResponse(request, "mysql_dashboard.html", {
        "request": request,
        "tables": tables,
        "fk_defs": fk_defs,
        "table_count": len(tables),
        "total_rows": sum(t["count"] for t in tables),
        **_mysql_template_context(source_id, db_name),
    })


@app.get("/mysql/{source_id}/{db_name}/table/{table_name}", response_class=HTMLResponse)
async def mysql_table_browser(request: Request, source_id: str, db_name: str, table_name: str,
                               page: int = 1, page_size: int = 25):
    engine = get_mysql_db_engine(db_name, source_id=source_id)
    meta = get_mysql_inspector(db_name, source_id=source_id)[table_name]

    columns_info = meta["columns"]
    pk_cols = meta["pk"]
    column_names = [c["name"] for c in columns_info]

    # Column ordering: PK first, then priority names, then FK columns, then alpha
    fks = meta["fks"]
    fk_col_names = set()
    for fk in fks:
        fk_col_names.update(fk["constrained_columns"])

    priority = list(pk_cols)
    for name in ["name", "symbol", "type", "description"]:
        if name in column_names and name not in priority:
            priority.append(name)
    fk_ordered = [c for c in column_names if c in fk_col_names and c not in priority]
    remaining = sorted(set(column_names) - set(priority) - set(fk_ordered))
    ordered_columns = priority + fk_ordered + remaining

    offset = (page - 1) * page_size
    with engine.connect() as conn:
        total = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).scalar()
        rows = conn.execute(
            text(f"SELECT * FROM `{table_name}` LIMIT :limit OFFSET :offset"),
            {"limit": page_size, "offset": offset}
        ).mappings().all()

    total_pages = max(1, (total + page_size - 1) // page_size)

    htmx = request.headers.get("HX-Request") == "true"
    template = "mysql_table_rows.html" if htmx else "mysql_table.html"

    return templates.TemplateResponse(request, template, {
        "request": request,
        "table_name": table_name,
        "rows": [dict(r) for r in rows],
        "columns": ordered_columns,
        "pk_cols": pk_cols,
        "total": total,
        "page": page,
        "page_size": page_size,
        "total_pages": total_pages,
        **_mysql_template_context(source_id, db_name),
    })


@app.get("/mysql/{source_id}/{db_name}/table/{table_name}/stats", response_class=HTMLResponse)
async def mysql_table_stats(request: Request, source_id: str, db_name: str, table_name: str):
    """Column coverage stats for a MySQL table (loaded via HTMX)."""
    engine = get_mysql_db_engine(db_name, source_id=source_id)
    columns = get_mysql_inspector(db_name, source_id=source_id)[table_name]["columns"]

    with engine.connect() as conn:
        total = conn.execute(text(f"SELECT COUNT(*) FROM `{table_name}`")).scalar()

        stats = []
        for col in columns:
            col_name = col["name"]
            non_null = conn.execute(
                text(f"SELECT COUNT(`{col_name}`) FROM `{table_name}`")
            ).scalar()
            pct = round(100 * non_null / total, 1) if total > 0 else 0
            stats.append({"field": col_name, "count": non_null, "pct": pct})

    return templates.TemplateResponse(request, "stats_partial.html", {
        "request": request,
        "stats": stats,
        "sample_size": total,
        "total": total,
    })


@app.get("/mysql/{source_id}/{db_name}/table/{table_name}/row/{pk_value:path}", response_class=HTMLResponse)
async def mysql_row_detail(request: Request, source_id: str, db_name: str, table_name: str, pk_value: str):
    engine = get_mysql_db_engine(db_name, source_id=source_id)
    meta = get_mysql_inspector(db_name, source_id=source_id)[table_name]
    pk_cols = meta["pk"]

    # Parse pk_value: "123" for single PK, "col1=val1/col2=val2" for composite
    where_clauses = []
    bind_params = {}
    if "=" in pk_value:
        parts = pk_value.split("/")
        for part in parts:
            col, val = part.split("=", 1)
            where_clauses.append(f"`{col}` = :pk_{col}")
            bind_params[f"pk_{col}"] = val
    elif pk_cols:
        where_clauses.append(f"`{pk_cols[0]}` = :pk_val")
        bind_params["pk_val"] = pk_value

    where_sql = " AND ".join(where_clauses) if where_clauses else "1=0"

    with engine.connect() as conn:
        result = conn.execute(text(f"SELECT * FROM `{table_name}` WHERE {where_sql}"), bind_params)
        row = result.mappings().first()

    # Build FK links for navigation
    fk_links = {}
    if row:
        for fk in meta["fks"]:
            for local_col, ref_col in zip(fk["constrained_columns"], fk["referred_columns"]):
                val = row.get(local_col)
                if val is not None:
                    fk_links[local_col] = {
                        "table": fk["referred_table"],
                        "column": ref_col,
                        "value": val,
                        "url": _app_path(f"/mysql/{source_id}/{db_name}/table/{fk['referred_table']}/row/{val}", request),
                    }

    # Get column metadata
    columns_info = {c["name"]: str(c["type"]) for c in meta["columns"]}

    return templates.TemplateResponse(request, "mysql_row.html", {
        "request": request,
        "table_name": table_name,
        "row": dict(row) if row else None,
        "pk_value": pk_value,
        "pk_cols": pk_cols,
        "fk_links": fk_links,
        "columns_info": columns_info,
        **_mysql_template_context(source_id, db_name),
    })


@app.get("/mysql/{source_id}/{db_name}/schema", response_class=HTMLResponse)
async def mysql_schema(request: Request, source_id: str, db_name: str):
    schema_meta = get_mysql_inspector(db_name, source_id=source_id)

    fk_defs = []
    mermaid_lines = ["erDiagram"]

    for table_name, meta in schema_meta.items():
        for fk in meta["fks"]:
            ref_table = fk["referred_table"]
            label = ", ".join(fk["constrained_columns"])
            safe_from = table_name.replace(" ", "_")
            safe_to = ref_table.replace(" ", "_")
            mermaid_lines.append(f'    {safe_from} }}|--|| {safe_to} : "{label}"')
            fk_defs.append({
                "from_table": table_name,
                "from_columns": ", ".join(fk["constrained_columns"]),
                "to_table": ref_table,
                "to_columns": ", ".join(fk["referred_columns"]),
            })

    mermaid_text = "\n".join(mermaid_lines)

    return templates.TemplateResponse(request, "mysql_schema.html", {
        "request": request,
        "fk_defs": fk_defs,
        "mermaid_text": mermaid_text,
        **_mysql_template_context(source_id, db_name),
    })


@app.post("/mysql/{source_id}/{db_name}/refresh-schema", response_class=HTMLResponse)
async def mysql_refresh_schema(request: Request, source_id: str, db_name: str):
    """Bust the schema cache for a database and redirect to dashboard."""
    invalidate_mysql_inspector(db_name, source_id=source_id)
    return _redirect_to(f"/mysql/{source_id}/{db_name}", request=request, status_code=303)


@app.get("/mysql/{source_id}/{db_name}/sql", response_class=HTMLResponse)
async def sql_page(request: Request, source_id: str, db_name: str):
    return templates.TemplateResponse(request, "mysql_sql.html", {
        "request": request,
        "results": None,
        "query": "",
        "error": None,
        "columns": [],
        **_mysql_template_context(source_id, db_name),
    })


@app.post("/mysql/{source_id}/{db_name}/sql", response_class=HTMLResponse)
async def sql_execute(request: Request, source_id: str, db_name: str, query: str = Form(...)):
    results = None
    error = None
    columns = []

    # Read-only guard
    query_upper = query.strip().upper()
    allowed = ("SELECT", "SHOW", "DESCRIBE", "DESC", "EXPLAIN")
    if not any(query_upper.startswith(kw) for kw in allowed):
        error = "Only SELECT, SHOW, DESCRIBE, and EXPLAIN queries are allowed."
    else:
        try:
            engine = get_mysql_db_engine(db_name, source_id=source_id)
            with engine.connect() as conn:
                result = conn.execute(text(query))
                if result.returns_rows:
                    columns = list(result.keys())
                    results = [dict(row) for row in result.mappings().all()]
                else:
                    results = []
        except Exception as e:
            error = str(e)

    htmx = request.headers.get("HX-Request") == "true"
    template = "mysql_sql_results.html" if htmx else "mysql_sql.html"

    return templates.TemplateResponse(request, template, {
        "request": request,
        "results": results,
        "query": query,
        "error": error,
        "columns": columns,
        **_mysql_template_context(source_id, db_name),
    })


@app.get("/mysql/{db_name}", response_class=HTMLResponse)
async def mysql_dashboard_default(request: Request, db_name: str):
    return _redirect_to(f"/mysql/default/{db_name}", request=request, status_code=307)


@app.get("/mysql/{db_name}/table/{table_name}", response_class=HTMLResponse)
async def mysql_table_browser_default(request: Request, db_name: str, table_name: str, page: int = 1, page_size: int = 25):
    return _redirect_to(
        f"/mysql/default/{db_name}/table/{table_name}?page={page}&page_size={page_size}",
        request=request,
        status_code=307,
    )


@app.get("/mysql/{db_name}/table/{table_name}/stats", response_class=HTMLResponse)
async def mysql_table_stats_default(request: Request, db_name: str, table_name: str):
    return _redirect_to(f"/mysql/default/{db_name}/table/{table_name}/stats", request=request, status_code=307)


@app.get("/mysql/{db_name}/table/{table_name}/row/{pk_value:path}", response_class=HTMLResponse)
async def mysql_row_detail_default(request: Request, db_name: str, table_name: str, pk_value: str):
    return _redirect_to(f"/mysql/default/{db_name}/table/{table_name}/row/{pk_value}", request=request, status_code=307)


@app.get("/mysql/{db_name}/schema", response_class=HTMLResponse)
async def mysql_schema_default(request: Request, db_name: str):
    return _redirect_to(f"/mysql/default/{db_name}/schema", request=request, status_code=307)


@app.post("/mysql/{db_name}/refresh-schema", response_class=HTMLResponse)
async def mysql_refresh_schema_default(request: Request, db_name: str):
    return _redirect_to(f"/mysql/default/{db_name}/refresh-schema", request=request, status_code=307)


@app.get("/mysql/{db_name}/sql", response_class=HTMLResponse)
async def sql_page_default(request: Request, db_name: str):
    return _redirect_to(f"/mysql/default/{db_name}/sql", request=request, status_code=307)


@app.post("/mysql/{db_name}/sql", response_class=HTMLResponse)
async def sql_execute_default(request: Request, db_name: str):
    return _redirect_to(f"/mysql/default/{db_name}/sql", request=request, status_code=307)


# ── Helpers ──────────────────────────────────────────────────────────────────

def _discover_columns(docs: list, is_edge: bool) -> list:
    """Pick the most useful columns from a set of documents."""
    if not docs:
        return []
    col_set = set()
    for doc in docs:
        col_set.update(doc.keys())

    # Remove internal arango keys except _key
    internal = {"_id", "_rev"}
    if not is_edge:
        internal.update({"_from", "_to"})
    col_set -= internal

    # Order: _key first, then _from/_to for edges, then id, then alpha
    priority = ["_key", "_from", "_to", "id", "name", "symbol", "type", "description"]
    ordered = [c for c in priority if c in col_set]
    remaining = sorted(col_set - set(ordered))
    return ordered + remaining


def _truncate(value, max_len=80):
    """Truncate a value for table display."""
    s = str(value)
    if len(s) > max_len:
        return s[:max_len] + "..."
    return s


def _normalize_csv_value(value):
    if isinstance(value, (dict, list)):
        return json.dumps(value, default=str)
    if value is None:
        return ""
    return value


def _categorize_document_fields(doc: dict | None):
    skip_keys = {"_key", "_id", "_rev", "_from", "_to", "sources", "provenance", "xref", "labels", "creation"}
    scalar_fields = []
    list_fields = []
    nested_fields = []

    if doc:
        for key, val in doc.items():
            if key in skip_keys or key.startswith("_"):
                continue
            if (
                val is None
                or isinstance(val, (str, int, float, bool))
                or (isinstance(val, dict) and val.get(_INLINE_CURATION_MARKER))
            ):
                scalar_fields.append((key, val))
            elif isinstance(val, dict):
                nested_fields.append((key, val))
            elif isinstance(val, list):
                if val and isinstance(val[0], dict):
                    nested_fields.append((key, val))
                elif val:
                    list_fields.append((key, val))
                else:
                    scalar_fields.append((key, val))
            else:
                scalar_fields.append((key, val))

    return scalar_fields, list_fields, nested_fields


def _get_adjacent_collection_docs(db, coll_name: str, doc: dict | None) -> dict:
    if not doc:
        return {"previous_doc": None, "next_doc": None, "adjacent_error": None}

    current_key = doc.get("_key")
    if not current_key:
        return {"previous_doc": None, "next_doc": None, "adjacent_error": "missing _key"}

    try:
        cursor = db.aql.execute(
            f"""
            LET previous_doc = FIRST(
                FOR d IN `{coll_name}`
                FILTER d._key < @current_key
                SORT d._key DESC
                LIMIT 1
                RETURN KEEP(d, "_key", "id", "name")
            )
            LET next_doc = FIRST(
                FOR d IN `{coll_name}`
                FILTER d._key > @current_key
                SORT d._key ASC
                LIMIT 1
                RETURN KEEP(d, "_key", "id", "name")
            )
            RETURN {{
                previous_doc: previous_doc,
                next_doc: next_doc
            }}
            """,
            bind_vars={"current_key": current_key},
            max_runtime=10,
        )
        result = next(iter(cursor), {"previous_doc": None, "next_doc": None})
        result["adjacent_error"] = None
        return result
    except Exception as exc:
        return {"previous_doc": None, "next_doc": None, "adjacent_error": str(exc)}


def _build_cure_case_report_context(db_name: str, db, doc: dict | None) -> dict:
    if not doc or not db.has_graph("graph"):
        return {"cure_case_url": _build_cure_case_url(db_name, doc)}

    try:
        cursor = db.aql.execute(
            """
            FOR v, e IN 1..1 OUTBOUND @node_id GRAPH 'graph'
            FILTER SPLIT(v._id, '/')[0] == 'Patient'
            FILTER SPLIT(e._id, '/')[0] == 'CaseReportPatientEdge'
            LIMIT 1
            RETURN v
            """,
            bind_vars={"node_id": doc["_id"]},
            max_runtime=10,
        )
        patient_doc = next(iter(cursor), None)
    except Exception:
        patient_doc = None

    try:
        reporter_cursor = db.aql.execute(
            """
            FOR v, e IN 1..1 OUTBOUND @node_id GRAPH 'graph'
            FILTER SPLIT(v._id, '/')[0] == 'Reporter'
            FILTER SPLIT(e._id, '/')[0] == 'CaseReportReporterEdge'
            LIMIT 1
            RETURN v
            """,
            bind_vars={"node_id": doc["_id"]},
            max_runtime=10,
        )
        reporter_doc = next(iter(reporter_cursor), None)
    except Exception:
        reporter_doc = None

    try:
        clinical_context_cursor = db.aql.execute(
            """
            LET via_patient = (
              FOR patient, report_edge IN 1..1 OUTBOUND @node_id GRAPH 'graph'
                FILTER SPLIT(patient._id, '/')[0] == 'Patient'
                FILTER SPLIT(report_edge._id, '/')[0] == 'CaseReportPatientEdge'
                FOR clinical_context, clinical_context_edge IN 1..1 OUTBOUND patient._id GRAPH 'graph'
                  FILTER SPLIT(clinical_context._id, '/')[0] == 'ClinicalContext'
                  FILTER SPLIT(clinical_context_edge._id, '/')[0] == 'PatientClinicalContextEdge'
                  LIMIT 1
                  RETURN clinical_context
            )
            LET via_report = (
              FOR clinical_context, clinical_context_edge IN 1..1 OUTBOUND @node_id GRAPH 'graph'
                FILTER SPLIT(clinical_context._id, '/')[0] == 'ClinicalContext'
                FILTER SPLIT(clinical_context_edge._id, '/')[0] == 'CaseReportClinicalContextEdge'
                LIMIT 1
                RETURN clinical_context
            )
            RETURN FIRST(APPEND(via_patient, via_report))
            """,
            bind_vars={"node_id": doc["_id"]},
            max_runtime=10,
        )
        clinical_context_doc = next(iter(clinical_context_cursor), None)
    except Exception:
        clinical_context_doc = None

    try:
        primary_episode_cursor = db.aql.execute(
            """
            FOR patient, report_edge IN 1..1 OUTBOUND @node_id GRAPH 'graph'
            FILTER SPLIT(patient._id, '/')[0] == 'Patient'
            FILTER SPLIT(report_edge._id, '/')[0] == 'CaseReportPatientEdge'
            FOR episode, patient_edge IN 1..1 OUTBOUND patient._id GRAPH 'graph'
            FILTER SPLIT(episode._id, '/')[0] == 'Episode'
            FILTER SPLIT(patient_edge._id, '/')[0] == 'PersonEpisodeEdge'
            FILTER episode.role == 'primary'
            LIMIT 1
            RETURN episode
            """,
            bind_vars={"node_id": doc["_id"]},
            max_runtime=10,
        )
        primary_episode_doc = next(iter(primary_episode_cursor), None)
    except Exception:
        primary_episode_doc = None

    try:
        acute_episode_cursor = db.aql.execute(
            """
            FOR patient, report_edge IN 1..1 OUTBOUND @node_id GRAPH 'graph'
            FILTER SPLIT(patient._id, '/')[0] == 'Patient'
            FILTER SPLIT(report_edge._id, '/')[0] == 'CaseReportPatientEdge'
            FOR episode, patient_edge IN 1..1 OUTBOUND patient._id GRAPH 'graph'
            FILTER SPLIT(episode._id, '/')[0] == 'Episode'
            FILTER SPLIT(patient_edge._id, '/')[0] == 'PersonEpisodeEdge'
            FILTER episode.role == 'contextual'
            FILTER LENGTH(
                FOR condition, condition_edge IN 1..1 OUTBOUND episode._id GRAPH 'graph'
                  FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                  FILTER SPLIT(condition_edge._id, '/')[0] == 'EpisodeConditionEdge'
                  FILTER condition.slug == 'acute-covid-19'
                  LIMIT 1
                  RETURN 1
            ) > 0
            LIMIT 1
            RETURN episode
            """,
            bind_vars={"node_id": doc["_id"]},
            max_runtime=10,
        )
        acute_episode_doc = next(iter(acute_episode_cursor), None)
    except Exception:
        acute_episode_doc = None

    try:
        pregnancy_episode_cursor = db.aql.execute(
            """
            FOR patient, report_edge IN 1..1 OUTBOUND @node_id GRAPH 'graph'
            FILTER SPLIT(patient._id, '/')[0] == 'Patient'
            FILTER SPLIT(report_edge._id, '/')[0] == 'CaseReportPatientEdge'
            FOR episode, patient_edge IN 1..1 OUTBOUND patient._id GRAPH 'graph'
            FILTER SPLIT(episode._id, '/')[0] == 'Episode'
            FILTER SPLIT(patient_edge._id, '/')[0] == 'PersonEpisodeEdge'
            FILTER episode.role == 'contextual'
            FILTER LENGTH(
                FOR condition, condition_edge IN 1..1 OUTBOUND episode._id GRAPH 'graph'
                  FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                  FILTER SPLIT(condition_edge._id, '/')[0] == 'EpisodeConditionEdge'
                  FILTER condition.slug == 'pregnancy'
                  LIMIT 1
                  RETURN 1
            ) > 0
            LIMIT 1
            RETURN episode
            """,
            bind_vars={"node_id": doc["_id"]},
            max_runtime=10,
        )
        pregnancy_episode_doc = next(iter(pregnancy_episode_cursor), None)
    except Exception:
        pregnancy_episode_doc = None

    patient_scalar_fields, patient_list_fields, patient_nested_fields = _categorize_document_fields(patient_doc)
    reporter_scalar_fields, reporter_list_fields, reporter_nested_fields = _categorize_document_fields(reporter_doc)
    clinical_context_scalar_fields, clinical_context_list_fields, clinical_context_nested_fields = _categorize_document_fields(clinical_context_doc)
    primary_episode_scalar_fields, primary_episode_list_fields, primary_episode_nested_fields = _categorize_document_fields(primary_episode_doc)
    acute_episode_scalar_fields, acute_episode_list_fields, acute_episode_nested_fields = _categorize_document_fields(acute_episode_doc)
    pregnancy_episode_scalar_fields, pregnancy_episode_list_fields, pregnancy_episode_nested_fields = _categorize_document_fields(pregnancy_episode_doc)

    episode_relationship_cards = []
    clinical_context_condition_cards = []
    clinical_context_phenotype_cards = []
    background_context_doc = None
    background_context_scalar_fields = []
    background_context_list_fields = []
    background_context_nested_fields = []
    perinatal_context_doc = None
    perinatal_context_scalar_fields = []
    perinatal_context_list_fields = []
    perinatal_context_nested_fields = []
    perinatal_context_phenotype_cards = []
    rasopathies_diagnosis_cards = []
    rasopathies_drug_treatment_cards = []
    patient_prior_condition_cards = []
    background_regular_medicine_cards = []
    background_immunosuppressant_cards = []
    primary_episode_post_covid_condition_cards = []
    exposure_cards = []
    acute_exposure_cards = []
    pregnancy_exposure_cards = []
    acute_complication_cards = []
    acute_vaccination_card = None
    phenotype_cards = []
    therapy_cards = []
    treatment_cards = []
    outcome_cards = []

    if clinical_context_doc:
        try:
            clinical_context_condition_cursor = db.aql.execute(
                """
                FOR condition, condition_edge IN 1..1 OUTBOUND @clinical_context_id GRAPH 'graph'
                FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                FILTER SPLIT(condition_edge._id, '/')[0] == 'ClinicalContextConditionEdge'
                SORT condition.name, condition._key
                RETURN {
                    condition: condition,
                    condition_edge: condition_edge
                }
                """,
                bind_vars={"clinical_context_id": clinical_context_doc["_id"]},
                max_runtime=15,
            )
            for row in clinical_context_condition_cursor:
                condition_doc = row.get("condition")
                condition_edge_doc = row.get("condition_edge")
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(condition_doc)
                clinical_context_condition_cards.append({
                    "condition_doc": condition_doc,
                    "condition_edge_doc": condition_edge_doc,
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                })
        except Exception:
            clinical_context_condition_cards = []

        try:
            clinical_context_phenotype_cursor = db.aql.execute(
                """
                LET via_finding = (
                  FOR finding, finding_edge IN 1..1 OUTBOUND @clinical_context_id GRAPH 'graph'
                    FILTER SPLIT(finding._id, '/')[0] == 'Finding'
                    FILTER SPLIT(finding_edge._id, '/')[0] == 'ClinicalContextFindingEdge'
                    FOR phenotype, phenotype_edge IN 1..1 OUTBOUND finding._id GRAPH 'graph'
                      FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                      FILTER SPLIT(phenotype_edge._id, '/')[0] == 'FindingPhenotypeEdge'
                      RETURN {
                        phenotype: phenotype,
                        finding: finding,
                        finding_edge: finding_edge,
                        phenotype_edge: phenotype_edge,
                        sort_group: finding.group,
                        sort_name: phenotype.name
                      }
                )
                LET direct = (
                  FOR phenotype, phenotype_edge IN 1..1 OUTBOUND @clinical_context_id GRAPH 'graph'
                    FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                    FILTER SPLIT(phenotype_edge._id, '/')[0] == 'ClinicalContextPhenotypeEdge'
                    RETURN {
                      phenotype: phenotype,
                      finding: null,
                      finding_edge: null,
                      phenotype_edge: phenotype_edge,
                      sort_group: phenotype_edge.group,
                      sort_name: phenotype.name
                    }
                )
                FOR row IN (LENGTH(via_finding) > 0 ? via_finding : direct)
                  SORT row.sort_group, row.sort_name, row.phenotype._key
                  RETURN row
                """,
                bind_vars={"clinical_context_id": clinical_context_doc["_id"]},
                max_runtime=15,
            )
            for row in clinical_context_phenotype_cursor:
                phenotype_doc = row.get("phenotype")
                finding_doc = row.get("finding")
                finding_edge_doc = row.get("finding_edge")
                phenotype_edge_doc = row.get("phenotype_edge")
                display_edge_doc = finding_doc or phenotype_edge_doc
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(phenotype_doc)
                edge_scalar_fields, edge_list_fields, edge_nested_fields = _categorize_document_fields(display_edge_doc)
                clinical_context_phenotype_cards.append({
                    "phenotype_doc": phenotype_doc,
                    "finding_doc": finding_doc,
                    "finding_edge_doc": finding_edge_doc,
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                    "edge_doc": display_edge_doc,
                    "phenotype_edge_doc": phenotype_edge_doc,
                    "edge_scalar_fields": edge_scalar_fields,
                    "edge_list_fields": edge_list_fields,
                    "edge_nested_fields": edge_nested_fields,
                })
        except Exception:
            clinical_context_phenotype_cards = []

        try:
            diagnosis_cursor = db.aql.execute(
                """
                FOR diagnosis, diagnosis_edge IN 1..1 OUTBOUND @clinical_context_id GRAPH 'graph'
                FILTER SPLIT(diagnosis._id, '/')[0] == 'Diagnosis'
                FILTER SPLIT(diagnosis_edge._id, '/')[0] == 'ClinicalContextDiagnosisEdge'
                LET conditions = (
                  FOR condition, condition_edge IN 1..1 OUTBOUND diagnosis._id GRAPH 'graph'
                    FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                    FILTER SPLIT(condition_edge._id, '/')[0] == 'DiagnosisConditionEdge'
                    SORT condition.name, condition._key
                    RETURN {
                      condition: condition,
                      condition_edge: condition_edge
                    }
                )
                LET genes = (
                  FOR gene, gene_edge IN 1..1 OUTBOUND diagnosis._id GRAPH 'graph'
                    FILTER SPLIT(gene._id, '/')[0] == 'Gene'
                    FILTER SPLIT(gene_edge._id, '/')[0] == 'DiagnosisGeneEdge'
                    SORT gene.symbol, gene.name, gene._key
                    RETURN {
                      gene: gene,
                      gene_edge: gene_edge
                    }
                )
                LET variants = (
                  FOR variant, variant_edge IN 1..1 OUTBOUND diagnosis._id GRAPH 'graph'
                    FILTER SPLIT(variant._id, '/')[0] == 'GeneVariant'
                    FILTER SPLIT(variant_edge._id, '/')[0] == 'DiagnosisGeneVariantEdge'
                    LET linked_genes = (
                      FOR gene, gene_variant_edge IN 1..1 INBOUND variant._id GRAPH 'graph'
                        FILTER SPLIT(gene._id, '/')[0] == 'Gene'
                        FILTER SPLIT(gene_variant_edge._id, '/')[0] == 'GeneGeneVariantEdge'
                        SORT gene.symbol, gene.name, gene._key
                        RETURN {
                          gene: gene,
                          gene_variant_edge: gene_variant_edge
                        }
                    )
                    SORT variant.source_gene_symbol, variant.variant_label, variant._key
                    RETURN {
                      variant: variant,
                      variant_edge: variant_edge,
                      linked_genes: linked_genes
                    }
                )
                SORT diagnosis._key
                RETURN {
                  diagnosis: diagnosis,
                  diagnosis_edge: diagnosis_edge,
                  conditions: conditions,
                  genes: genes,
                  variants: variants
                }
                """,
                bind_vars={"clinical_context_id": clinical_context_doc["_id"]},
                max_runtime=15,
            )
            for row in diagnosis_cursor:
                diagnosis_doc = row.get("diagnosis")
                diagnosis_scalar_fields, diagnosis_list_fields, diagnosis_nested_fields = _categorize_document_fields(diagnosis_doc)
                condition_cards = []
                for condition_row in row.get("conditions") or []:
                    condition_doc = condition_row.get("condition")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(condition_doc)
                    condition_cards.append({
                        "condition_doc": condition_doc,
                        "condition_edge_doc": condition_row.get("condition_edge"),
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
                gene_cards = []
                for gene_row in row.get("genes") or []:
                    gene_doc = gene_row.get("gene")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(gene_doc)
                    gene_cards.append({
                        "gene_doc": gene_doc,
                        "gene_edge_doc": gene_row.get("gene_edge"),
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
                variant_cards = []
                for variant_row in row.get("variants") or []:
                    variant_doc = variant_row.get("variant")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(variant_doc)
                    linked_gene_cards = []
                    for linked_gene_row in variant_row.get("linked_genes") or []:
                        gene_doc = linked_gene_row.get("gene")
                        linked_gene_cards.append({
                            "gene_doc": gene_doc,
                            "gene_variant_edge_doc": linked_gene_row.get("gene_variant_edge"),
                        })
                    variant_cards.append({
                        "variant_doc": variant_doc,
                        "variant_edge_doc": variant_row.get("variant_edge"),
                        "linked_gene_cards": linked_gene_cards,
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
                rasopathies_diagnosis_cards.append({
                    "diagnosis_doc": diagnosis_doc,
                    "diagnosis_edge_doc": row.get("diagnosis_edge"),
                    "diagnosis_scalar_fields": diagnosis_scalar_fields,
                    "diagnosis_list_fields": diagnosis_list_fields,
                    "diagnosis_nested_fields": diagnosis_nested_fields,
                    "condition_cards": condition_cards,
                    "gene_cards": gene_cards,
                    "variant_cards": variant_cards,
                })
        except Exception:
            rasopathies_diagnosis_cards = []

    finding_group_order = [
        "Cardiac",
        "Endocrine/Growth",
        "Gastrointestinal",
        "Hematologic/Oncologic",
        "Lymphatic/Immunologic",
        "Neurologic/Audiologic",
        "Opthalmalogic",
        "Diagnoses not listed above",
    ]
    finding_group_rank = {name: index for index, name in enumerate(finding_group_order)}
    grouped_clinical_context_phenotype_cards = []
    if clinical_context_phenotype_cards:
        grouped = {}
        for card in clinical_context_phenotype_cards:
            group_name = ((card.get("edge_doc") or {}).get("group") or "Ungrouped").strip()
            grouped.setdefault(group_name, []).append(card)
        for group_name, cards in grouped.items():
            cards.sort(
                key=lambda card: (
                    ((card.get("phenotype_doc") or {}).get("name") or "").lower(),
                    ((card.get("phenotype_doc") or {}).get("_key") or ""),
                )
            )
        grouped_clinical_context_phenotype_cards = [
            {
                "group_name": group_name,
                "cards": grouped[group_name],
            }
            for group_name in sorted(
                grouped.keys(),
                key=lambda name: (
                    finding_group_rank.get(name, 10_000),
                    name == "Ungrouped",
                    name.lower(),
                ),
            )
        ]

    if primary_episode_doc:
        if acute_episode_doc:
            try:
                episode_relationship_cursor = db.aql.execute(
                    """
                    FOR episode, episode_edge IN 1..1 OUTBOUND @acute_episode_id GRAPH 'graph'
                    FILTER SPLIT(episode._id, '/')[0] == 'Episode'
                    FILTER SPLIT(episode_edge._id, '/')[0] == 'EpisodeEpisodeEdge'
                    FILTER episode_edge.relationship_type == 'precedes'
                    FILTER episode._id == @primary_episode_id
                    LIMIT 1
                    RETURN episode_edge
                    """,
                    bind_vars={
                        "acute_episode_id": acute_episode_doc["_id"],
                        "primary_episode_id": primary_episode_doc["_id"],
                    },
                    max_runtime=10,
                )
                relationship_doc = next(iter(episode_relationship_cursor), None)
                if relationship_doc:
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(relationship_doc)
                    episode_relationship_cards.append({
                        "relationship_doc": relationship_doc,
                        "label": "precedes",
                        "left_tab": "acute-covid",
                        "right_tab": "long-covid",
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
            except Exception:
                pass
        if pregnancy_episode_doc:
            try:
                overlap_relationship_cursor = db.aql.execute(
                    """
                    FOR episode, episode_edge IN 1..1 OUTBOUND @pregnancy_episode_id GRAPH 'graph'
                    FILTER SPLIT(episode._id, '/')[0] == 'Episode'
                    FILTER SPLIT(episode_edge._id, '/')[0] == 'EpisodeEpisodeEdge'
                    FILTER episode_edge.relationship_type == 'overlaps'
                    FILTER episode._id == @primary_episode_id
                    LIMIT 1
                    RETURN episode_edge
                    """,
                    bind_vars={
                        "pregnancy_episode_id": pregnancy_episode_doc["_id"],
                        "primary_episode_id": primary_episode_doc["_id"],
                    },
                    max_runtime=10,
                )
                relationship_doc = next(iter(overlap_relationship_cursor), None)
                if relationship_doc:
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(relationship_doc)
                    episode_relationship_cards.append({
                        "relationship_doc": relationship_doc,
                        "label": "overlaps",
                        "left_tab": "long-covid",
                        "right_tab": "pregnancy",
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
            except Exception:
                pass

        if patient_doc:
            try:
                background_context_cursor = db.aql.execute(
                    """
                    FOR background_context, context_edge IN 1..1 OUTBOUND @person_id GRAPH 'graph'
                    FILTER SPLIT(background_context._id, '/')[0] == 'BackgroundContext'
                    FILTER SPLIT(context_edge._id, '/')[0] == 'PersonBackgroundContextEdge'
                    LIMIT 1
                    RETURN background_context
                    """,
                    bind_vars={"person_id": patient_doc["_id"]},
                    max_runtime=15,
                )
                background_context_doc = next(iter(background_context_cursor), None)
            except Exception:
                background_context_doc = None

            background_context_scalar_fields, background_context_list_fields, background_context_nested_fields = _categorize_document_fields(background_context_doc)

            if background_context_doc:
                try:
                    patient_condition_cursor = db.aql.execute(
                        """
                        FOR condition, condition_edge IN 1..1 OUTBOUND @background_context_id GRAPH 'graph'
                        FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                        FILTER SPLIT(condition_edge._id, '/')[0] == 'BackgroundContextConditionEdge'
                        FILTER condition_edge.relationship_type == 'prior_comorbidity'
                        SORT condition.name, condition._key
                        RETURN condition
                        """,
                        bind_vars={"background_context_id": background_context_doc["_id"]},
                        max_runtime=15,
                    )
                    for condition_doc in patient_condition_cursor:
                        scalar_fields, list_fields, nested_fields = _categorize_document_fields(condition_doc)
                        patient_prior_condition_cards.append({
                            "condition_doc": condition_doc,
                            "scalar_fields": scalar_fields,
                            "list_fields": list_fields,
                            "nested_fields": nested_fields,
                        })
                except Exception:
                    patient_prior_condition_cards = []

                try:
                    background_exposure_cursor = db.aql.execute(
                        """
                        FOR exposure, context_edge IN 1..1 OUTBOUND @background_context_id GRAPH 'graph'
                        FILTER SPLIT(exposure._id, '/')[0] == 'Exposure'
                        FILTER SPLIT(context_edge._id, '/')[0] == 'BackgroundContextExposureEdge'
                        LET drug_doc = FIRST(
                            FOR drug, drug_edge IN 1..1 OUTBOUND exposure._id GRAPH 'graph'
                            FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                            FILTER SPLIT(drug_edge._id, '/')[0] == 'ExposureDrugEdge'
                            LIMIT 1
                            RETURN drug
                        )
                        SORT context_edge.relationship_type, drug_doc.name, exposure.long_drug_name, exposure._key
                        RETURN {
                            exposure_doc: exposure,
                            drug_doc: drug_doc,
                            relationship_type: context_edge.relationship_type
                        }
                        """,
                        bind_vars={"background_context_id": background_context_doc["_id"]},
                        max_runtime=15,
                    )
                    for row in background_exposure_cursor:
                        exposure_doc = row.get("exposure_doc")
                        scalar_fields, list_fields, nested_fields = _categorize_document_fields(exposure_doc)
                        card = {
                            "exposure_doc": exposure_doc,
                            "drug_doc": row.get("drug_doc"),
                            "scalar_fields": scalar_fields,
                            "list_fields": list_fields,
                            "nested_fields": nested_fields,
                        }
                        if row.get("relationship_type") == "immunosuppressant":
                            background_immunosuppressant_cards.append(card)
                        else:
                            background_regular_medicine_cards.append(card)
                except Exception:
                    background_regular_medicine_cards = []
                    background_immunosuppressant_cards = []

    if patient_doc:
        try:
            perinatal_context_cursor = db.aql.execute(
                """
                FOR perinatal_context, context_edge IN 1..1 OUTBOUND @person_id GRAPH 'graph'
                FILTER SPLIT(perinatal_context._id, '/')[0] == 'PerinatalContext'
                FILTER SPLIT(context_edge._id, '/')[0] == 'PatientPerinatalContextEdge'
                LIMIT 1
                RETURN perinatal_context
                """,
                bind_vars={"person_id": patient_doc["_id"]},
                max_runtime=15,
            )
            perinatal_context_doc = next(iter(perinatal_context_cursor), None)
        except Exception:
            perinatal_context_doc = None

        perinatal_context_scalar_fields, perinatal_context_list_fields, perinatal_context_nested_fields = _categorize_document_fields(perinatal_context_doc)
        if perinatal_context_doc:
            try:
                perinatal_phenotype_cursor = db.aql.execute(
                    """
                    LET via_finding = (
                      FOR finding, finding_edge IN 1..1 OUTBOUND @perinatal_context_id GRAPH 'graph'
                        FILTER SPLIT(finding._id, '/')[0] == 'Finding'
                        FILTER SPLIT(finding_edge._id, '/')[0] == 'PerinatalContextFindingEdge'
                        FOR phenotype, phenotype_edge IN 1..1 OUTBOUND finding._id GRAPH 'graph'
                          FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                          FILTER SPLIT(phenotype_edge._id, '/')[0] == 'FindingPhenotypeEdge'
                          RETURN {
                            phenotype: phenotype,
                            finding: finding,
                            finding_edge: finding_edge,
                            phenotype_edge: phenotype_edge,
                            sort_name: phenotype.name
                          }
                    )
                    LET direct = (
                      FOR phenotype, phenotype_edge IN 1..1 OUTBOUND @perinatal_context_id GRAPH 'graph'
                        FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                        FILTER SPLIT(phenotype_edge._id, '/')[0] == 'PerinatalContextPhenotypeEdge'
                        RETURN {
                          phenotype: phenotype,
                          finding: null,
                          finding_edge: null,
                          phenotype_edge: phenotype_edge,
                          sort_name: phenotype.name
                        }
                    )
                    FOR row IN (LENGTH(via_finding) > 0 ? via_finding : direct)
                      SORT row.sort_name, row.phenotype._key
                      RETURN row
                    """,
                    bind_vars={"perinatal_context_id": perinatal_context_doc["_id"]},
                    max_runtime=15,
                )
                for row in perinatal_phenotype_cursor:
                    phenotype_doc = row.get("phenotype")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(phenotype_doc)
                    perinatal_context_phenotype_cards.append({
                        "phenotype_doc": phenotype_doc,
                        "finding_doc": row.get("finding"),
                        "finding_edge_doc": row.get("finding_edge"),
                        "edge_doc": row.get("finding") or row.get("phenotype_edge"),
                        "phenotype_edge_doc": row.get("phenotype_edge"),
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
            except Exception:
                perinatal_context_phenotype_cards = []

        try:
            rasopathies_treatment_cursor = db.aql.execute(
                """
                FOR treatment, treatment_edge IN 1..1 OUTBOUND @clinical_context_id GRAPH 'graph'
                FILTER SPLIT(treatment._id, '/')[0] == 'DrugTreatment'
                FILTER SPLIT(treatment_edge._id, '/')[0] == 'ClinicalContextDrugTreatmentEdge'
                LET drug_doc = FIRST(
                  FOR drug, drug_edge IN 1..1 OUTBOUND treatment._id GRAPH 'graph'
                    FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                    FILTER SPLIT(drug_edge._id, '/')[0] == 'DrugTreatmentDrugEdge'
                    LIMIT 1
                    RETURN drug
                )
                LET responses = (
                  FOR response, response_edge IN 1..1 OUTBOUND treatment._id GRAPH 'graph'
                    FILTER SPLIT(response._id, '/')[0] == 'TreatmentResponse'
                    FILTER SPLIT(response_edge._id, '/')[0] == 'DrugTreatmentResponseEdge'
                    LET finding_doc = FIRST(
                      FOR finding, finding_edge IN 1..1 OUTBOUND response._id GRAPH 'graph'
                        FILTER SPLIT(finding._id, '/')[0] == 'Finding'
                        FILTER SPLIT(finding_edge._id, '/')[0] == 'TreatmentResponseFindingEdge'
                        LIMIT 1
                        RETURN finding
                    )
                    LET phenotype_docs = (
                      FOR phenotype, phenotype_edge IN 1..1 OUTBOUND finding_doc._id GRAPH 'graph'
                        FILTER finding_doc != null
                        FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                        FILTER SPLIT(phenotype_edge._id, '/')[0] == 'FindingPhenotypeEdge'
                        SORT phenotype.name, phenotype._key
                        RETURN {
                          phenotype: phenotype,
                          phenotype_edge: phenotype_edge
                        }
                    )
                    SORT response.source_target_index, response._key
                    RETURN {
                      response_doc: response,
                      response_edge_doc: response_edge,
                      finding_doc: finding_doc,
                      phenotype_docs: phenotype_docs
                    }
                )
                LET adverse_events = (
                  FOR phenotype, adverse_event_edge IN 1..1 OUTBOUND treatment._id GRAPH 'graph'
                    FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                    FILTER SPLIT(adverse_event_edge._id, '/')[0] == 'DrugTreatmentAdverseEventEdge'
                    SORT adverse_event_edge.source_adverse_event_index, phenotype.name, phenotype._key
                    RETURN {
                      phenotype: phenotype,
                      adverse_event_edge: adverse_event_edge
                    }
                )
                SORT treatment.source_treatment_index, drug_doc.name, treatment._key
                RETURN {
                  treatment_doc: treatment,
                  treatment_edge_doc: treatment_edge,
                  drug_doc: drug_doc,
                  responses: responses,
                  adverse_events: adverse_events
                }
                """,
                bind_vars={"clinical_context_id": clinical_context_doc["_id"]},
                max_runtime=20,
            )
            for row in rasopathies_treatment_cursor:
                treatment_doc = row.get("treatment_doc")
                treatment_scalar_fields, treatment_list_fields, treatment_nested_fields = _categorize_document_fields(treatment_doc)
                drug_doc = row.get("drug_doc")
                drug_scalar_fields, drug_list_fields, drug_nested_fields = _categorize_document_fields(drug_doc)
                response_cards = []
                for response_row in row.get("responses") or []:
                    response_doc = response_row.get("response_doc")
                    response_scalar_fields, response_list_fields, response_nested_fields = _categorize_document_fields(response_doc)
                    finding_doc = response_row.get("finding_doc")
                    finding_scalar_fields, finding_list_fields, finding_nested_fields = _categorize_document_fields(finding_doc)
                    phenotype_cards = []
                    for phenotype_row in response_row.get("phenotype_docs") or []:
                        phenotype_doc = phenotype_row.get("phenotype")
                        phenotype_scalar_fields, phenotype_list_fields, phenotype_nested_fields = _categorize_document_fields(phenotype_doc)
                        phenotype_cards.append({
                            "phenotype_doc": phenotype_doc,
                            "phenotype_edge_doc": phenotype_row.get("phenotype_edge"),
                            "scalar_fields": phenotype_scalar_fields,
                            "list_fields": phenotype_list_fields,
                            "nested_fields": phenotype_nested_fields,
                        })
                    response_cards.append({
                        "response_doc": response_doc,
                        "response_edge_doc": response_row.get("response_edge_doc"),
                        "scalar_fields": response_scalar_fields,
                        "list_fields": response_list_fields,
                        "nested_fields": response_nested_fields,
                        "finding_doc": finding_doc,
                        "finding_scalar_fields": finding_scalar_fields,
                        "finding_list_fields": finding_list_fields,
                        "finding_nested_fields": finding_nested_fields,
                        "phenotype_cards": phenotype_cards,
                    })
                adverse_event_cards = []
                for adverse_event_row in row.get("adverse_events") or []:
                    adverse_event_doc = adverse_event_row.get("phenotype")
                    adverse_event_edge_doc = adverse_event_row.get("adverse_event_edge")
                    adverse_event_scalar_fields, adverse_event_list_fields, adverse_event_nested_fields = _categorize_document_fields(adverse_event_doc)
                    adverse_event_edge_scalar_fields, adverse_event_edge_list_fields, adverse_event_edge_nested_fields = _categorize_document_fields(adverse_event_edge_doc)
                    adverse_event_cards.append({
                        "phenotype_doc": adverse_event_doc,
                        "adverse_event_edge_doc": adverse_event_edge_doc,
                        "scalar_fields": adverse_event_scalar_fields,
                        "list_fields": adverse_event_list_fields,
                        "nested_fields": adverse_event_nested_fields,
                        "edge_scalar_fields": adverse_event_edge_scalar_fields,
                        "edge_list_fields": adverse_event_edge_list_fields,
                        "edge_nested_fields": adverse_event_edge_nested_fields,
                    })
                rasopathies_drug_treatment_cards.append({
                    "treatment_doc": treatment_doc,
                    "treatment_edge_doc": row.get("treatment_edge_doc"),
                    "drug_doc": drug_doc,
                    "treatment_scalar_fields": treatment_scalar_fields,
                    "treatment_list_fields": treatment_list_fields,
                    "treatment_nested_fields": treatment_nested_fields,
                    "drug_scalar_fields": drug_scalar_fields,
                    "drug_list_fields": drug_list_fields,
                    "drug_nested_fields": drug_nested_fields,
                    "response_cards": response_cards,
                    "adverse_event_cards": adverse_event_cards,
                })
        except Exception:
            rasopathies_drug_treatment_cards = []

        try:
            primary_post_covid_condition_cursor = db.aql.execute(
                """
                FOR condition, condition_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                FILTER SPLIT(condition_edge._id, '/')[0] == 'EpisodeConditionEdge'
                FILTER condition_edge.relationship_type == 'comorbidity'
                SORT condition.name, condition._key
                RETURN condition
                """,
                bind_vars={"episode_id": primary_episode_doc["_id"]},
                max_runtime=15,
            )
            for condition_doc in primary_post_covid_condition_cursor:
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(condition_doc)
                primary_episode_post_covid_condition_cards.append({
                    "condition_doc": condition_doc,
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                })
        except Exception:
            primary_episode_post_covid_condition_cards = []

        try:
            exposure_cursor = db.aql.execute(
                """
                FOR exposure, episode_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                FILTER SPLIT(exposure._id, '/')[0] == 'Exposure'
                FILTER SPLIT(episode_edge._id, '/')[0] == 'EpisodeExposureEdge'
                LET drug_doc = FIRST(
                    FOR drug, drug_edge IN 1..1 OUTBOUND exposure._id GRAPH 'graph'
                    FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                    FILTER SPLIT(drug_edge._id, '/')[0] == 'ExposureDrugEdge'
                    LIMIT 1
                    RETURN drug
                )
                LET adverse_events = (
                    FOR ae, ae_edge IN 1..1 OUTBOUND exposure._id GRAPH 'graph'
                    FILTER SPLIT(ae._id, '/')[0] IN ['Phenotype', 'AdverseEvent']
                    FILTER SPLIT(ae_edge._id, '/')[0] == 'ExposureAdverseEventEdge'
                    SORT ae.name
                    RETURN {
                        id: ae._id,
                        key: ae._key,
                        collection: SPLIT(ae._id, '/')[0],
                        name: ae.name,
                        outcomes: ae_edge.outcomes || []
                    }
                )
                SORT exposure.long_drug_name, exposure._key
                RETURN {
                    exposure_doc: exposure,
                    drug_doc: drug_doc,
                    adverse_events: adverse_events
                }
                """,
                bind_vars={"episode_id": primary_episode_doc["_id"]},
                max_runtime=15,
            )
            for row in exposure_cursor:
                exposure_doc = row.get("exposure_doc")
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(exposure_doc)
                exposure_cards.append({
                    "exposure_doc": exposure_doc,
                    "drug_doc": row.get("drug_doc"),
                    "adverse_events": row.get("adverse_events", []),
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                })
        except Exception:
            exposure_cards = []

        if acute_episode_doc:
            try:
                acute_exposure_cursor = db.aql.execute(
                    """
                    FOR exposure, episode_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                    FILTER SPLIT(exposure._id, '/')[0] == 'Exposure'
                    FILTER SPLIT(episode_edge._id, '/')[0] == 'EpisodeExposureEdge'
                    LET drug_doc = FIRST(
                        FOR drug, drug_edge IN 1..1 OUTBOUND exposure._id GRAPH 'graph'
                        FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                        FILTER SPLIT(drug_edge._id, '/')[0] == 'ExposureDrugEdge'
                        LIMIT 1
                        RETURN drug
                    )
                    SORT exposure.long_drug_name, exposure._key
                    RETURN {
                        exposure_doc: exposure,
                        drug_doc: drug_doc,
                        adverse_events: []
                    }
                    """,
                    bind_vars={"episode_id": acute_episode_doc["_id"]},
                    max_runtime=15,
                )
                for row in acute_exposure_cursor:
                    exposure_doc = row.get("exposure_doc")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(exposure_doc)
                    acute_exposure_cards.append({
                        "exposure_doc": exposure_doc,
                        "drug_doc": row.get("drug_doc"),
                        "adverse_events": [],
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
            except Exception:
                acute_exposure_cards = []

            try:
                acute_complication_cursor = db.aql.execute(
                    """
                    FOR condition, complication_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                    FILTER SPLIT(condition._id, '/')[0] == 'Condition'
                    FILTER SPLIT(complication_edge._id, '/')[0] == 'EpisodeConditionEdge'
                    FILTER complication_edge.relationship_type == 'complication'
                    SORT condition.name, condition._key
                    RETURN condition
                    """,
                    bind_vars={"episode_id": acute_episode_doc["_id"]},
                    max_runtime=15,
                )
                for condition_doc in acute_complication_cursor:
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(condition_doc)
                    acute_complication_cards.append({
                        "condition_doc": condition_doc,
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
            except Exception:
                acute_complication_cards = []

            try:
                acute_vaccination_cursor = db.aql.execute(
                    """
                    FOR vaccination_event, event_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                    FILTER SPLIT(vaccination_event._id, '/')[0] == 'VaccinationEvent'
                    FILTER SPLIT(event_edge._id, '/')[0] == 'VaccinationEventEpisodeEdge'
                    LET vaccines = (
                        FOR vaccine, vaccine_edge IN 1..1 OUTBOUND vaccination_event._id GRAPH 'graph'
                        FILTER SPLIT(vaccine._id, '/')[0] == 'Vaccine'
                        FILTER SPLIT(vaccine_edge._id, '/')[0] == 'VaccinationEventVaccineEdge'
                        SORT vaccine.name, vaccine._key
                        RETURN vaccine
                    )
                    LIMIT 1
                    RETURN {
                        vaccination_event_doc: vaccination_event,
                        event_edge: event_edge,
                        vaccines: vaccines
                    }
                    """,
                    bind_vars={"episode_id": acute_episode_doc["_id"]},
                    max_runtime=15,
                )
                acute_vaccination_row = next(iter(acute_vaccination_cursor), None)
                if acute_vaccination_row:
                    vaccination_event_doc = acute_vaccination_row.get("vaccination_event_doc")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(vaccination_event_doc)
                    acute_vaccination_card = {
                        "vaccination_event_doc": vaccination_event_doc,
                        "event_edge": acute_vaccination_row.get("event_edge") or {},
                        "vaccines": acute_vaccination_row.get("vaccines", []),
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    }
            except Exception:
                acute_vaccination_card = None

        if pregnancy_episode_doc:
            try:
                pregnancy_exposure_cursor = db.aql.execute(
                    """
                    FOR exposure, episode_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                    FILTER SPLIT(exposure._id, '/')[0] == 'Exposure'
                    FILTER SPLIT(episode_edge._id, '/')[0] == 'EpisodeExposureEdge'
                    LET drug_doc = FIRST(
                        FOR drug, drug_edge IN 1..1 OUTBOUND exposure._id GRAPH 'graph'
                        FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                        FILTER SPLIT(drug_edge._id, '/')[0] == 'ExposureDrugEdge'
                        LIMIT 1
                        RETURN drug
                    )
                    SORT exposure.long_drug_name, exposure._key
                    RETURN {
                        exposure_doc: exposure,
                        drug_doc: drug_doc,
                        adverse_events: []
                    }
                    """,
                    bind_vars={"episode_id": pregnancy_episode_doc["_id"]},
                    max_runtime=15,
                )
                for row in pregnancy_exposure_cursor:
                    exposure_doc = row.get("exposure_doc")
                    scalar_fields, list_fields, nested_fields = _categorize_document_fields(exposure_doc)
                    pregnancy_exposure_cards.append({
                        "exposure_doc": exposure_doc,
                        "drug_doc": row.get("drug_doc"),
                        "adverse_events": [],
                        "scalar_fields": scalar_fields,
                        "list_fields": list_fields,
                        "nested_fields": nested_fields,
                    })
            except Exception:
                pregnancy_exposure_cards = []

        try:
            phenotype_cursor = db.aql.execute(
                """
                FOR phenotype, phenotype_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                FILTER SPLIT(phenotype_edge._id, '/')[0] == 'EpisodePhenotypeEdge'
                SORT phenotype.name, phenotype._key
                RETURN {
                    phenotype_doc: phenotype,
                    severity: phenotype_edge.severity
                }
                """,
                bind_vars={"episode_id": primary_episode_doc["_id"]},
                max_runtime=15,
            )
            for row in phenotype_cursor:
                phenotype_doc = row.get("phenotype_doc")
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(phenotype_doc)
                phenotype_cards.append({
                    "phenotype_doc": phenotype_doc,
                    "severity": row.get("severity"),
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                })
        except Exception:
            phenotype_cards = []

        try:
            therapy_cursor = db.aql.execute(
                """
                FOR therapy, therapy_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                FILTER SPLIT(therapy._id, '/')[0] == 'Therapy'
                FILTER SPLIT(therapy_edge._id, '/')[0] == 'EpisodeTherapyEdge'
                SORT therapy.name, therapy._key
                RETURN therapy
                """,
                bind_vars={"episode_id": primary_episode_doc["_id"]},
                max_runtime=15,
            )
            for therapy_doc in therapy_cursor:
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(therapy_doc)
                therapy_cards.append({
                    "therapy_doc": therapy_doc,
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                })
        except Exception:
            therapy_cards = []

        try:
            treatment_cursor = db.aql.execute(
                """
                FOR exposure, episode_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                FILTER SPLIT(exposure._id, '/')[0] == 'Exposure'
                FILTER SPLIT(episode_edge._id, '/')[0] == 'EpisodeExposureEdge'
                FOR treatment, treatment_edge IN 1..1 OUTBOUND exposure._id GRAPH 'graph'
                FILTER SPLIT(treatment._id, '/')[0] == 'Treatment'
                FILTER SPLIT(treatment_edge._id, '/')[0] == 'TreatmentExposureEdge'
                COLLECT treatment_id = treatment._id INTO grouped = {
                    treatment: treatment,
                    exposure: exposure
                }
                LET treatment_doc = FIRST(grouped[*].treatment)
                LET mapped_exposures = (
                    FOR row IN grouped
                    LET drug_doc = FIRST(
                        FOR drug, drug_edge IN 1..1 OUTBOUND row.exposure._id GRAPH 'graph'
                        FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                        FILTER SPLIT(drug_edge._id, '/')[0] == 'ExposureDrugEdge'
                        LIMIT 1
                        RETURN drug
                    )
                    SORT drug_doc.name, row.exposure.long_drug_name, row.exposure._key
                    RETURN {
                        exposure_id: row.exposure._id,
                        drug_name: drug_doc.name || row.exposure.long_drug_name || row.exposure._key
                    }
                )
                SORT treatment_doc._key
                RETURN {
                    treatment_doc: treatment_doc,
                    mapped_exposures: mapped_exposures
                }
                """,
                bind_vars={"episode_id": primary_episode_doc["_id"]},
                max_runtime=15,
            )
            for row in treatment_cursor:
                treatment_doc = row.get("treatment_doc")
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(treatment_doc)
                treatment_cards.append({
                    "treatment_doc": treatment_doc,
                    "mapped_exposures": row.get("mapped_exposures", []),
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                })
        except Exception:
            treatment_cards = []

        try:
            outcome_cursor = db.aql.execute(
                """
                FOR outcome, episode_outcome_edge IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                FILTER SPLIT(outcome._id, '/')[0] == 'Outcome'
                FILTER SPLIT(episode_outcome_edge._id, '/')[0] == 'EpisodeOutcomeEdge'
                LET treatment_doc = FIRST(
                    FOR treatment, outcome_edge IN 1..1 INBOUND outcome._id GRAPH 'graph'
                    FILTER SPLIT(treatment._id, '/')[0] == 'Treatment'
                    FILTER SPLIT(outcome_edge._id, '/')[0] == 'TreatmentOutcomeEdge'
                    LIMIT 1
                    RETURN treatment
                )
                LET phenotype_doc = FIRST(
                    FOR phenotype, phenotype_edge IN 1..1 OUTBOUND outcome._id GRAPH 'graph'
                    FILTER SPLIT(phenotype._id, '/')[0] == 'Phenotype'
                    FILTER SPLIT(phenotype_edge._id, '/')[0] == 'OutcomePhenotypeEdge'
                    LIMIT 1
                    RETURN phenotype
                )
                LET phenotype_severity = FIRST(
                    FOR phenotype2, phenotype_edge2 IN 1..1 OUTBOUND @episode_id GRAPH 'graph'
                    FILTER SPLIT(phenotype2._id, '/')[0] == 'Phenotype'
                    FILTER SPLIT(phenotype_edge2._id, '/')[0] == 'EpisodePhenotypeEdge'
                    FILTER phenotype_doc != null AND phenotype2._id == phenotype_doc._id
                    LIMIT 1
                    RETURN phenotype_edge2.severity
                )
                LET mapped_exposures = (
                    FOR exposure2, exposure_edge2 IN 1..1 INBOUND treatment_doc._id GRAPH 'graph'
                    FILTER SPLIT(exposure2._id, '/')[0] == 'Exposure'
                    FILTER SPLIT(exposure_edge2._id, '/')[0] == 'TreatmentExposureEdge'
                    LET drug_doc = FIRST(
                        FOR drug, drug_edge IN 1..1 OUTBOUND exposure2._id GRAPH 'graph'
                        FILTER SPLIT(drug._id, '/')[0] == 'Drug'
                        FILTER SPLIT(drug_edge._id, '/')[0] == 'ExposureDrugEdge'
                        LIMIT 1
                        RETURN drug
                    )
                    SORT drug_doc.name, exposure2.long_drug_name, exposure2._key
                    RETURN {
                        exposure_id: exposure2._id,
                        exposure_doc: exposure2,
                        drug_doc: drug_doc,
                        drug_name: drug_doc.name || exposure2.long_drug_name || exposure2._key
                    }
                )
                SORT outcome._key
                RETURN {
                    outcome_doc: outcome,
                    treatment_doc: treatment_doc,
                    phenotype_doc: phenotype_doc,
                    phenotype_severity: phenotype_severity,
                    mapped_exposures: mapped_exposures
                }
                """,
                bind_vars={"episode_id": primary_episode_doc["_id"]},
                max_runtime=20,
            )
            for row in outcome_cursor:
                outcome_doc = row.get("outcome_doc")
                effect_display = _get_outcome_effect_display((outcome_doc or {}).get("effect"))
                scalar_fields, list_fields, nested_fields = _categorize_document_fields(outcome_doc)
                phenotype_doc = row.get("phenotype_doc")
                phenotype_scalar_fields, phenotype_list_fields, phenotype_nested_fields = _categorize_document_fields(phenotype_doc)
                treatment_doc = row.get("treatment_doc")
                treatment_scalar_fields, treatment_list_fields, treatment_nested_fields = _categorize_document_fields(treatment_doc)
                outcome_cards.append({
                    "outcome_doc": outcome_doc,
                    "scalar_fields": scalar_fields,
                    "list_fields": list_fields,
                    "nested_fields": nested_fields,
                    "treatment_doc": treatment_doc,
                    "treatment_scalar_fields": treatment_scalar_fields,
                    "treatment_list_fields": treatment_list_fields,
                    "treatment_nested_fields": treatment_nested_fields,
                    "mapped_exposures": row.get("mapped_exposures", []),
                    "phenotype_doc": phenotype_doc,
                    "phenotype_severity": row.get("phenotype_severity"),
                    "phenotype_scalar_fields": phenotype_scalar_fields,
                    "phenotype_list_fields": phenotype_list_fields,
                    "phenotype_nested_fields": phenotype_nested_fields,
                    "effect_rank": effect_display["rank"],
                    "effect_tone": effect_display["tone"],
                    "effect_pct": effect_display["pct"],
                })
            outcome_cards.sort(
                key=lambda card: (
                    card["effect_rank"],
                    (card["outcome_doc"] or {}).get("effect") or "",
                    (card["outcome_doc"] or {}).get("raw_symptom_name") or "",
                    (card["outcome_doc"] or {}).get("_key") or "",
                )
            )
        except Exception:
            outcome_cards = []

    return {
        "cure_case_url": _build_cure_case_url(db_name, doc),
        "patient_doc": patient_doc,
        "patient_scalar_fields": patient_scalar_fields,
        "patient_list_fields": patient_list_fields,
        "patient_nested_fields": patient_nested_fields,
        "reporter_doc": reporter_doc,
        "reporter_scalar_fields": reporter_scalar_fields,
        "reporter_list_fields": reporter_list_fields,
        "reporter_nested_fields": reporter_nested_fields,
        "clinical_context_doc": clinical_context_doc,
        "clinical_context_scalar_fields": clinical_context_scalar_fields,
        "clinical_context_list_fields": clinical_context_list_fields,
        "clinical_context_nested_fields": clinical_context_nested_fields,
        "clinical_context_condition_cards": clinical_context_condition_cards,
        "clinical_context_phenotype_cards": clinical_context_phenotype_cards,
        "grouped_clinical_context_phenotype_cards": grouped_clinical_context_phenotype_cards,
        "background_context_doc": background_context_doc,
        "background_context_scalar_fields": background_context_scalar_fields,
        "background_context_list_fields": background_context_list_fields,
        "background_context_nested_fields": background_context_nested_fields,
        "perinatal_context_doc": perinatal_context_doc,
        "perinatal_context_scalar_fields": perinatal_context_scalar_fields,
        "perinatal_context_list_fields": perinatal_context_list_fields,
        "perinatal_context_nested_fields": perinatal_context_nested_fields,
        "perinatal_context_phenotype_cards": perinatal_context_phenotype_cards,
        "rasopathies_diagnosis_cards": rasopathies_diagnosis_cards,
        "rasopathies_drug_treatment_cards": rasopathies_drug_treatment_cards,
        "episode_relationship_cards": episode_relationship_cards,
        "patient_prior_condition_cards": patient_prior_condition_cards,
        "background_regular_medicine_cards": background_regular_medicine_cards,
        "background_immunosuppressant_cards": background_immunosuppressant_cards,
        "primary_episode_doc": primary_episode_doc,
        "primary_episode_scalar_fields": primary_episode_scalar_fields,
        "primary_episode_list_fields": primary_episode_list_fields,
        "primary_episode_nested_fields": primary_episode_nested_fields,
        "primary_episode_post_covid_condition_cards": primary_episode_post_covid_condition_cards,
        "acute_episode_doc": acute_episode_doc,
        "acute_episode_scalar_fields": acute_episode_scalar_fields,
        "acute_episode_list_fields": acute_episode_list_fields,
        "acute_episode_nested_fields": acute_episode_nested_fields,
        "pregnancy_episode_doc": pregnancy_episode_doc,
        "pregnancy_episode_scalar_fields": pregnancy_episode_scalar_fields,
        "pregnancy_episode_list_fields": pregnancy_episode_list_fields,
        "pregnancy_episode_nested_fields": pregnancy_episode_nested_fields,
        "pregnancy_episode_exposure_cards": pregnancy_exposure_cards,
        "acute_episode_vaccination_card": acute_vaccination_card,
        "acute_episode_exposure_cards": acute_exposure_cards,
        "acute_episode_complication_cards": acute_complication_cards,
        "primary_episode_exposure_cards": exposure_cards,
        "primary_episode_phenotype_cards": phenotype_cards,
        "primary_episode_therapy_cards": therapy_cards,
        "primary_episode_treatment_cards": treatment_cards,
        "primary_episode_outcome_cards": outcome_cards,
    }


def _build_cure_case_url(db_name: str, doc: dict | None) -> str | None:
    if not doc:
        return None
    if doc.get("case_report_url"):
        return doc.get("case_report_url")
    report_id = doc.get("id") or doc.get("_key")
    if not report_id:
        return None

    route_slug = _get_cure_route_slug(db_name, doc.get("form_type"))
    if not route_slug:
        return None
    return f"https://cure.ncats.io/explore/{route_slug}/case-reports/case-details/{report_id}"


def _is_cure_db(db_name: str) -> bool:
    normalized = (db_name or "").strip().lower()
    return normalized == "cure" or normalized.startswith("cure_")


def _get_cure_route_slug(db_name: str, form_type: str | None = None) -> str | None:
    normalized_db_name = (db_name or "").strip().lower()
    if normalized_db_name == "cure" or normalized_db_name.startswith("cure_pasc"):
        return "long-covid"
    if normalized_db_name.startswith("cure_rasopathies"):
        return "rasopathies"

    return {
        "pasc": "long-covid",
        "rasopathies": "rasopathies",
    }.get((form_type or "").strip().lower())


def _get_outcome_effect_display(effect: str | None) -> dict:
    normalized = (effect or "").strip().lower()
    mapping = {
        "complete symptom resolution": {"rank": 0, "tone": "positive", "pct": 100},
        "significant symptom improvement": {"rank": 1, "tone": "positive", "pct": 82},
        "moderate symptom improvement": {"rank": 2, "tone": "positive", "pct": 64},
        "mild symptom improvement": {"rank": 3, "tone": "positive", "pct": 46},
        "symptom was unchanged": {"rank": 4, "tone": "neutral", "pct": 0},
        "unknown": {"rank": 5, "tone": "unknown", "pct": 0},
        "mild symptom worsening": {"rank": 6, "tone": "negative", "pct": 40},
        "moderate symptom worsening": {"rank": 7, "tone": "negative", "pct": 62},
        "significant symptom worsening": {"rank": 8, "tone": "negative", "pct": 84},
    }
    default = {"rank": 9, "tone": "unknown", "pct": 0}
    return mapping.get(normalized, default)


templates.env.filters["truncate_val"] = _truncate


# ── Document template dispatch ───────────────────────────────────────────────

def _get_document_template(db_name: str, coll_name: str) -> str:
    if coll_name == "CaseReport" and _is_cure_db(db_name):
        return "cure_case_report_document.html"
    return "document.html"


# ── Demo routes ──────────────────────────────────────────────────────────────

if _demo_queries_enabled:
    import src.qa_browser.demo_routes as _demo_module  # noqa: E402
    app.include_router(_demo_module.router)
    _demo_module.set_templates(templates)

# ── POUNCE validation routes ─────────────────────────────────────────────────

import src.qa_browser.pounce_routes as _pounce_module  # noqa: E402
app.include_router(_pounce_module.router)
_pounce_module.set_templates(templates)
_pounce_module.set_mysql_engine_getter(get_mysql_db_engine)
# set_pounce_config is called in main() after args are parsed

# ── Feedback routes ───────────────────────────────────────────────────────────

import src.qa_browser.feedback_routes as _feedback_module  # noqa: E402
app.include_router(_feedback_module.router)
_feedback_module.set_templates(templates)
# set_feedback_file is called in main() after args are parsed


# ── Entrypoint ───────────────────────────────────────────────────────────────

def main():
    parser = argparse.ArgumentParser(description="QA Browser")
    parser.add_argument("--credentials", "-c",
                        default="./src/use_cases/secrets/local_arangodb.yaml",
                        help="Path to ArangoDB credentials YAML file")
    parser.add_argument("--mysql-credentials", "-m",
                        action="append",
                        default=[],
                        help="Path to a MySQL credentials YAML file; repeat to load multiple MySQL servers")
    parser.add_argument("--object-storage-credentials", "-s",
                        default=None,
                        help="Path to object-storage credentials YAML file")
    parser.add_argument("--parquet-storage-credentials",
                        default=None,
                        help="Path to object-storage credentials YAML for existing Dataset parquet files")
    parser.add_argument("--port", "-p", type=int, default=8050)
    parser.add_argument("--host", default="127.0.0.1")
    parser.add_argument("--root-path", default="", help="ASGI root path for running behind a sub-path proxy (e.g. /odin-qa)")
    parser.add_argument("--pounce-config", "-P", default="./src/use_cases/pounce/pounce.yaml",
                        help="Path to pounce.yaml — used to load resolvers for mapping coverage checks")
    parser.add_argument("--smtp-credentials", "-e",
                        default=None,
                        help="Path to SMTP credentials YAML (host, port, user, password, from_address, to_address, use_tls)")
    parser.add_argument("--feedback-file", "-f",
                        default=None,
                        help="Path to JSON file for storing feedback comments (created if missing)")
    parser.add_argument("--ramp-diagnosis-file",
                        default=os.getenv("QA_BROWSER_RAMP_DIAGNOSIS_FILE", "input_files/ramp_id_diagnoses.json"),
                        help="Path to JSON file for storing legacy RaMP merge diagnoses")
    parser.add_argument("--pounce-project-base-url",
                        default="",
                        help="Public base URL for project detail links, e.g. https://pounce-ci.ncats.nih.gov/project")
    parser.add_argument("--disease-graph-dir",
                        default="",
                        help="Path to disease app_graph/ directory (TSV files + manifest.json)")
    parser.add_argument("--disease-review-file",
                        default="",
                        help="Path to IFX Harmonizers-compatible disease review intake TSV written by the Review tab")
    parser.add_argument("--baseline-graph-dir",
                        default="",
                        help="Path to baseline disease app_graph/ directory for version diff")
    parser.add_argument("--target-graph-dir",
                        default="",
                        help="Path to target app_graph/ directory (target_nodes.tsv + manifest.json)")
    parser.add_argument("--target-qc-dir",
                        default="",
                        help="Path to target pipeline qc/ directory (divergence registries + provenance CSVs)")
    parser.add_argument("--variant-graph-dir",
                        default="",
                        help="Path to variant app_graph/ directory (variant_nodes.tsv + manifest.json)")
    parser.add_argument("--variant-review-file",
                        default="",
                        help="Path to IFX Harmonizers-compatible variant review intake TSV written by the Review tab")
    parser.add_argument("--drug-graph-dir",
                        default="",
                        help="Path to drug app_graph/ directory (drug_nodes.tsv + manifest.json)")
    parser.add_argument("--drug-resolver-index",
                        default="",
                        help="Path to the version-matched complete drug_resolver_index.sqlite sidecar")
    parser.add_argument("--verify-drug-resolver-checksum",
                        action="store_true",
                        help="Verify the complete drug resolver sidecar SHA-256 at startup (recommended for staging/release checks)")
    parser.add_argument("--drug-review-file",
                        default="",
                        help="Path to IFX Harmonizers-compatible drug review intake TSV written by the Review tab")
    parser.add_argument("--pathway-graph-dir",
                        default="",
                        help="Path to pathway app_graph/ directory (pathway_nodes.tsv + manifest.json)")
    parser.add_argument("--pathway-categories-file",
                        default="",
                        help="Path to bioplanet_categories.csv for functional category annotations")
    parser.add_argument("--pathway-review-file",
                        default="",
                        help="Path to IFX Harmonizers-compatible pathway review intake TSV written by the Review tab")
    args = parser.parse_args()

    global _credentials, _mysql_credentials, _mysql_sources, _minio_credentials, _object_storage_credentials, _parquet_storage_credentials, _disease_graph_dir, _disease_review_file, _baseline_graph_dir, _target_graph_dir, _target_qc_dir, _variant_graph_dir, _variant_review_file, _drug_graph_dir, _drug_resolver_index, _verify_drug_resolver_checksum, _drug_review_file, _pathway_graph_dir, _pathway_categories_file, _pathway_review_file
    templates.env.globals["root_path"] = args.root_path.rstrip("/")
    cred_path = Path(args.credentials)
    if cred_path.exists():
        with open(cred_path) as f:
            _credentials = yaml.safe_load(f)
        print(f"Loaded ArangoDB credentials from {cred_path}")
    else:
        print(f"Warning: {cred_path} not found, using defaults")
        _credentials = {"url": "http://localhost:8529", "user": "root", "password": "password"}

    for index, mysql_cred_path in enumerate(args.mysql_credentials, start=1):
        mysql_path = Path(mysql_cred_path)
        if not mysql_path.exists():
            print(f"Warning: MySQL credentials file {mysql_path} not found")
            continue

        with open(mysql_path) as f:
            creds = yaml.safe_load(f) or {}

        source_id = "default" if not _mysql_sources else _slugify_mysql_source(mysql_path.stem)
        if source_id in _mysql_sources:
            suffix = 2
            while f"{source_id}-{suffix}" in _mysql_sources:
                suffix += 1
            source_id = f"{source_id}-{suffix}"

        label = "default" if source_id == "default" else mysql_path.stem
        _register_mysql_source(source_id, creds, label=label)
        if index == 1:
            _mysql_credentials = creds
        print(f"Loaded MySQL credentials from {mysql_path} as source '{source_id}'")

    if _demo_queries_enabled:
        _demo_module.set_mysql_credentials(_mysql_credentials)

    if args.object_storage_credentials:
        storage_path = Path(args.object_storage_credentials)
        if storage_path.exists():
            with open(storage_path) as f:
                _object_storage_credentials = yaml.safe_load(f)
            print(f"Loaded object-storage credentials from {storage_path}")
        else:
            print(f"Warning: object-storage credentials file {storage_path} not found")

    if args.parquet_storage_credentials:
        parquet_storage_path = Path(args.parquet_storage_credentials)
        if parquet_storage_path.exists():
            with open(parquet_storage_path) as f:
                _parquet_storage_credentials = yaml.safe_load(f)
            print(f"Loaded parquet storage credentials from {parquet_storage_path}")
        else:
            print(f"Warning: parquet storage credentials file {parquet_storage_path} not found")

    _pounce_module.set_pounce_config(args.pounce_config)
    _pounce_module.set_smtp_config(args.smtp_credentials)
    _pounce_module.set_public_project_base_url(args.pounce_project_base_url)
    _feedback_module.set_feedback_file(args.feedback_file)
    set_ramp_diagnosis_file(args.ramp_diagnosis_file)

    _disease_graph_dir = args.disease_graph_dir
    if not _disease_graph_dir:
        bundled_dir = DISEASE_APP_GRAPH_BUNDLED_DIR
        versions = _versioned_app_graph_dirs(
            bundled_dir,
            "disease_concepts.tsv",
            ("disease_harmonizer_version", "disease_release_version", "pipeline_version", "version"),
        )
        if versions:
            _disease_graph_dir = str(versions[-1])
            print(f"Auto-detected bundled disease data: {_disease_graph_dir}")
    if _disease_graph_dir:
        print(f"Disease graph dir: {_disease_graph_dir}")
    _disease_review_file = args.disease_review_file
    if _resolved_disease_review_file():
        print(f"Disease review intake file: {_resolved_disease_review_file()}")

    _baseline_graph_dir = args.baseline_graph_dir
    if not _baseline_graph_dir:
        bundled_dir = DISEASE_APP_GRAPH_BUNDLED_DIR
        versions = _versioned_app_graph_dirs(
            bundled_dir,
            "disease_concepts.tsv",
            ("disease_harmonizer_version", "disease_release_version", "pipeline_version", "version"),
        )
        # Use the second-highest version as baseline (highest is current)
        if len(versions) >= 2:
            _baseline_graph_dir = str(versions[-2])
            print(f"Auto-detected baseline disease data: {_baseline_graph_dir}")
    if _baseline_graph_dir:
        print(f"Baseline graph dir: {_baseline_graph_dir}")

    _target_graph_dir = args.target_graph_dir
    if not _target_graph_dir:
        bundled_dir = TARGET_APP_GRAPH_BUNDLED_DIR
        versions = _versioned_app_graph_dirs(
            bundled_dir,
            "target_nodes.tsv",
            ("target_release_version", "target_harmonizer_version", "version"),
        )
        if versions:
            _target_graph_dir = str(versions[-1])
            print(f"Auto-detected bundled target data: {_target_graph_dir}")
    if _target_graph_dir:
        print(f"Target graph dir: {_target_graph_dir}")
    _target_qc_dir = args.target_qc_dir
    if _target_qc_dir:
        print(f"Target QC dir: {_target_qc_dir}")

    _variant_graph_dir = args.variant_graph_dir
    if not _variant_graph_dir:
        bundled_dir = VARIANT_APP_GRAPH_BUNDLED_DIR
        versions = _versioned_app_graph_dirs(
            bundled_dir,
            "variant_nodes.tsv",
            ("variant_harmonizer_version", "variant_release_version", "version"),
        )
        if versions:
            _variant_graph_dir = str(versions[-1])
            print(f"Auto-detected bundled variant data: {_variant_graph_dir}")
    if _variant_graph_dir:
        print(f"Variant graph dir: {_variant_graph_dir}")
    _variant_review_file = args.variant_review_file
    if _resolved_variant_review_file():
        print(f"Variant review intake file: {_resolved_variant_review_file()}")

    _drug_graph_dir = args.drug_graph_dir
    _drug_resolver_index = args.drug_resolver_index
    _verify_drug_resolver_checksum = args.verify_drug_resolver_checksum
    if not _drug_graph_dir:
        bundled_dir = DRUG_APP_GRAPH_BUNDLED_DIR
        versions = _versioned_app_graph_dirs(
            bundled_dir,
            "drug_nodes.tsv",
            ("drug_harmonizer_version", "drug_release_version", "version"),
        )
        if versions:
            _drug_graph_dir = str(versions[-1])
            print(f"Auto-detected bundled drug data: {_drug_graph_dir}")
    if _drug_graph_dir:
        print(f"Drug graph dir: {_drug_graph_dir}")
    if _drug_resolver_index:
        print(f"Drug resolver index: {_drug_resolver_index}")
        resolver_metadata = _validate_drug_resolver_index(
            Path(_drug_resolver_index),
            _drug_manifest_snapshot(),
            verify_checksum=_verify_drug_resolver_checksum,
        )
        print(
            "Validated drug resolver index: "
            f"v{str(resolver_metadata.get('harmonizer_version', '')).lstrip('v')}, "
            f"{resolver_metadata.get('node_count', 'unknown')} nodes, "
            f"scope={resolver_metadata.get('scope', 'unknown')}"
        )
    _drug_review_file = args.drug_review_file
    if _resolved_drug_review_file():
        print(f"Drug review intake file: {_resolved_drug_review_file()}")

    _pathway_graph_dir = args.pathway_graph_dir
    if not _pathway_graph_dir:
        bundled_dir = PATHWAY_APP_GRAPH_BUNDLED_DIR
        versions = _versioned_app_graph_dirs(
            bundled_dir,
            "pathway_nodes.tsv",
            ("pathway_harmonizer_version", "version", "pipeline_version"),
        )
        if versions:
            _pathway_graph_dir = str(versions[-1])
            print(f"Auto-detected bundled pathway data: {_pathway_graph_dir}")
    if _pathway_graph_dir:
        print(f"Pathway graph dir: {_pathway_graph_dir}")
    _pathway_categories_file = args.pathway_categories_file
    if _pathway_graph_dir and _pathway_categories_file:
        try:
            data = load_pathway_graph_data(_pathway_graph_dir)
            load_pathway_categories(_pathway_categories_file, data)
            print(f"Loaded pathway categories: {len(data.categories_by_pathway)} pathways mapped")
        except Exception as exc:
            print(f"Warning: could not load pathway categories: {exc}")
    _pathway_review_file = args.pathway_review_file
    if _resolved_pathway_review_file():
        print(f"Pathway review intake file: {_resolved_pathway_review_file()}")

    interrupted_run_keys = _try_reconcile_interrupted_harmonization_runs(force=True)
    if interrupted_run_keys:
        print(
            "Marked interrupted harmonization runs as failed: "
            + ", ".join(interrupted_run_keys)
        )

    print(f"Starting QA Browser at http://{args.host}:{args.port}")
    uvicorn.run(app, host=args.host, port=args.port, root_path=args.root_path)


if __name__ == "__main__":
    main()
