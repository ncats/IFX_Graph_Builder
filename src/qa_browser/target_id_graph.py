"""Target harmonizer app graph loader and explorer helpers."""
from __future__ import annotations

import csv
import io
import json
import logging
import os
import re
import threading
import uuid
from collections import Counter, OrderedDict, defaultdict
from datetime import datetime, timezone
from difflib import SequenceMatcher
from pathlib import Path
from typing import Any

from fastapi import HTTPException


_VALIDATION_SOURCES = {"nodenorm", "error"}


class TargetGraphData:
    __slots__ = (
        "nodes",
        "nodes_by_id",
        "ids_to_targets",
        "relation_edges",
        "relation_edges_by_target",
        "relation_edge_index",
        "type_index",
        "exact_terms",
        "token_index",
        "manifest",
        "_stats",
        "divergence_gene",
        "divergence_protein",
        "divergence_transcript",
        "_divergence_stats",
        "gene_biotype_counts",
        "qc_dir",
        "public_qc_summary",
        "review_mode",
        "review_can_write",
        "divergence_by_target_id",
    )

    def __init__(self) -> None:
        self.nodes: list[dict[str, str]] = []
        self.nodes_by_id: dict[str, dict[str, str]] = {}
        self.ids_to_targets: dict[str, list[str]] = defaultdict(list)
        self.relation_edges: list[dict[str, str]] = []
        self.relation_edges_by_target: dict[str, list[dict[str, str]]] = defaultdict(list)
        self.relation_edge_index: dict[tuple[str, str, str], dict[str, str]] = {}
        self.type_index: dict[str, list[str]] = defaultdict(list)
        self.exact_terms: dict[str, set[str]] = defaultdict(set)
        self.token_index: dict[str, set[str]] = defaultdict(set)
        self.manifest: dict[str, Any] = {}
        self._stats: dict[str, Any] | None = None
        self.divergence_gene: list[dict[str, str]] = []
        self.divergence_protein: list[dict[str, str]] = []
        self.divergence_transcript: list[dict[str, str]] = []
        self._divergence_stats: dict[str, Any] | None = None
        self.gene_biotype_counts: dict[str, int] = {}
        self.qc_dir: str = ""
        self.public_qc_summary: dict[str, Any] = {}
        self.review_mode: str = "none"
        self.review_can_write: bool = False
        self.divergence_by_target_id: dict[str, list[dict[str, str]]] = defaultdict(list)


_singleton: TargetGraphData | None = None
_singleton_lock = threading.Lock()
_target_review_lock = threading.Lock()

TARGET_REVIEW_DECISION_OPTIONS = [
    {"value": "accept_consensus", "label": "Accept consensus"},
    {"value": "accept_sole_source", "label": "Accept sole source"},
    {"value": "accept_low_score", "label": "Accept low score"},
    {"value": "trust_source", "label": "Trust source value"},
    {"value": "trust_refseq", "label": "Trust RefSeq"},
    {"value": "trust_uniprot", "label": "Trust UniProt"},
    {"value": "trust_ensembl", "label": "Trust Ensembl"},
    {"value": "trust_ncbi", "label": "Trust NCBI"},
    {"value": "trust_hgnc", "label": "Trust HGNC"},
    {"value": "cosmetic_ignore", "label": "Cosmetic / ignore"},
    {"value": "use_replacement", "label": "Use replacement"},
    {"value": "needs_expert_review", "label": "Needs expert review"},
]

# ── Triage categories ──
# Classify each divergence row into a triage bucket so users can batch-decide
# the easy ones first.  Categories ordered from easiest to hardest.
TRIAGE_CATEGORIES = [
    {"key": "coverage_gap", "label": "Coverage Gap",
     "description": "Only one source provides this field — expected, not an error."},
    {"key": "low_mapping", "label": "Low Mapping Support",
     "description": "Mapping support ratio is below threshold but may still be correct."},
    {"key": "cosmetic", "label": "Cosmetic Difference",
     "description": "Minor formatting or notation difference between sources."},
    {"key": "source_conflict", "label": "Source Conflict",
     "description": "Sources disagree on a substantive value — needs human review."},
]
_TRIAGE_ORDER = {cat["key"]: i for i, cat in enumerate(TRIAGE_CATEGORIES)}

_SOURCE_TO_DECISION: dict[str, str] = {
    "refseq": "trust_refseq",
    "uniprot": "trust_uniprot",
    "ensembl": "trust_ensembl",
    "ncbi": "trust_ncbi",
    "hgnc": "trust_hgnc",
}

_TRIAGE_SUGGESTED_DECISION: dict[str, str] = {
    "coverage_gap": "accept_sole_source",
    "low_mapping": "accept_low_score",
    "cosmetic": "cosmetic_ignore",
    "source_conflict": "",  # default — overridden per-row by xref suggestion
}


REVIEW_GROUPS = [
    {
        "key": "protein_refseq_authority",
        "label": "Protein RefSeq Authority",
        "description": "RefSeq protein-accession disagreements where the QC policy already prefers RefSeq for RefSeq IDs.",
        "recommended_decision": "trust_refseq",
    },
    {
        "key": "protein_uniprot_authority",
        "label": "Protein UniProt Authority",
        "description": "UniProt accession disagreements where the QC policy already prefers UniProt for UniProt IDs.",
        "recommended_decision": "trust_uniprot",
    },
    {
        "key": "protein_ensembl_authority",
        "label": "Protein Ensembl Authority",
        "description": "Ensembl protein-ID disagreements where the QC policy already prefers Ensembl for Ensembl IDs.",
        "recommended_decision": "trust_ensembl",
    },
    {
        "key": "protein_isoform_context",
        "label": "Protein Isoform Context",
        "description": "Protein identifier conflicts needing isoform/canonical-product context before accepting one source.",
        "recommended_decision": "needs_expert_review",
    },
    {
        "key": "gene_ncbi_authority",
        "label": "Gene NCBI Authority",
        "description": "NCBI Gene identifier conflicts where the QC policy already prefers NCBI for NCBI Gene IDs.",
        "recommended_decision": "trust_ncbi",
    },
    {
        "key": "gene_hgnc_authority",
        "label": "Gene HGNC Authority",
        "description": "HGNC identifier or symbol conflicts where the QC policy already prefers HGNC.",
        "recommended_decision": "trust_hgnc",
    },
    {
        "key": "gene_symbol_conflict",
        "label": "Gene Symbol Conflict",
        "description": "Gene symbol disagreements, including aliases and changed symbols, that need symbol-history review.",
        "recommended_decision": "needs_expert_review",
    },
    {
        "key": "gene_hgnc_identifier_conflict",
        "label": "Gene HGNC ID Conflict",
        "description": "HGNC identifier disagreements without a current high-confidence authority decision.",
        "recommended_decision": "needs_expert_review",
    },
    {
        "key": "gene_omim_identifier_conflict",
        "label": "Gene OMIM ID Conflict",
        "description": "OMIM identifier disagreements that should be checked against gene/disease inheritance context.",
        "recommended_decision": "needs_expert_review",
    },
    {
        "key": "gene_ensembl_identifier_conflict",
        "label": "Gene Ensembl ID Conflict",
        "description": "Ensembl gene identifier disagreements without a current high-confidence authority decision.",
        "recommended_decision": "needs_expert_review",
    },
    {
        "key": "coverage_gap",
        "label": "Coverage Gap",
        "description": "Only one source provides the field. This is usually expected coverage difference.",
        "recommended_decision": "accept_sole_source",
    },
    {
        "key": "low_mapping_support",
        "label": "Low Mapping Support",
        "description": "Identifier mapping support is below threshold but the row can often be accepted as low-confidence evidence.",
        "recommended_decision": "accept_low_score",
    },
    {
        "key": "cosmetic_metadata",
        "label": "Cosmetic Metadata",
        "description": "Formatting, punctuation, or metadata-only differences that do not change the identifier mapping.",
        "recommended_decision": "cosmetic_ignore",
    },
    {
        "key": "other_conflict",
        "label": "Other Conflict",
        "description": "Conflict pattern not yet assigned to a more specific review group.",
        "recommended_decision": "needs_expert_review",
    },
]
_REVIEW_GROUP_BY_KEY = {group["key"]: group for group in REVIEW_GROUPS}
_REVIEW_GROUP_ORDER = {group["key"]: i for i, group in enumerate(REVIEW_GROUPS)}


def classify_triage_category(row: dict[str, str]) -> str:
    """Assign a triage category to a divergence row.

    Returns the cached ``_triage_category`` if already computed (set once
    at load time), otherwise computes it fresh.

    Categories:
        coverage_gap   — sole_source with no scenario (expected gaps)
        low_mapping    — low_score divergence type
        cosmetic       — already auto-decided as cosmetic_ignore, or
                         description-only conflicts (scenario G06)
        source_conflict — everything else (real conflicts needing review)
    """
    cached = row.get("_triage_category")
    if cached:
        return cached

    div_type = (row.get("divergence_type", "") or "").lower()
    auto_dec = (row.get("auto_decision", "") or "").lower()
    scenario = (row.get("scenario_id", "") or "").upper()

    if div_type == "sole_source":
        return "coverage_gap"
    if div_type == "low_score":
        return "low_mapping"
    if auto_dec == "cosmetic_ignore" or scenario in ("G06", "G07"):
        return "cosmetic"
    return "source_conflict"


def classify_review_group(entity_type: str, row: dict[str, str]) -> str:
    """Assign a specific review group within the broader triage category.

    The group is deliberately policy-facing: it tells a reviewer what kind of
    decision can be made in bulk and which conflicts still need scientific
    inspection.
    """
    cached = row.get("_review_group")
    if cached:
        return cached

    etype = (entity_type or "").strip().lower()
    div_type = (row.get("divergence_type", "") or "").lower()
    auto_dec = (row.get("auto_decision", "") or "").lower()
    scenario = (row.get("scenario_id", "") or "").upper()
    namespace = row.get("namespace", "") or ""

    if div_type == "sole_source":
        return "coverage_gap"
    if div_type == "low_score":
        return "low_mapping_support"
    if auto_dec == "cosmetic_ignore" or scenario in ("G06", "G07"):
        return "cosmetic_metadata"

    if etype == "protein":
        if auto_dec == "trust_refseq":
            return "protein_refseq_authority"
        if auto_dec == "trust_uniprot":
            return "protein_uniprot_authority"
        if auto_dec == "trust_ensembl":
            return "protein_ensembl_authority"
        if namespace in {"Ensembl_Protein_ID", "RefSeq_Protein_ID", "UniProt_ID"}:
            return "protein_isoform_context"

    if etype == "gene":
        if auto_dec == "trust_ncbi":
            return "gene_ncbi_authority"
        if auto_dec == "trust_hgnc":
            return "gene_hgnc_authority"
        if namespace == "Symbol":
            return "gene_symbol_conflict"
        if namespace == "HGNC_ID":
            return "gene_hgnc_identifier_conflict"
        if namespace == "OMIM_ID":
            return "gene_omim_identifier_conflict"
        if namespace == "Ensembl_ID":
            return "gene_ensembl_identifier_conflict"

    return "other_conflict"


def _compute_review_group_summary(
    data: TargetGraphData,
    review_group: str,
    entity_type: str = "",
    status: str = "open",
) -> dict[str, Any]:
    """Compute pattern summary for a review group."""
    auto_decisions: Counter[str] = Counter()
    confidences: list[float] = []
    scenarios: Counter[str] = Counter()
    namespaces: Counter[str] = Counter()
    sample_rationales: list[str] = []
    total = 0

    for etype, rows in [
        ("gene", data.divergence_gene),
        ("protein", data.divergence_protein),
        ("transcript", data.divergence_transcript),
    ]:
        if entity_type and etype != entity_type:
            continue
        for row in rows:
            if classify_review_group(etype, row) != review_group:
                continue
            if status and (row.get("status", "") or "open").lower() != status:
                continue
            total += 1
            ad = row.get("auto_decision", "") or ""
            if ad:
                auto_decisions[ad] += 1
            try:
                confidences.append(float(row.get("auto_confidence", "")))
            except (ValueError, TypeError):
                pass
            sc = row.get("scenario_id", "") or ""
            if sc:
                scenarios[sc] += 1
            ns = row.get("namespace", "") or ""
            if ns:
                namespaces[ns] += 1
            rat = row.get("auto_rationale", "") or ""
            if rat and len(sample_rationales) < 3:
                sample_rationales.append(rat)

    if not total:
        return {}

    top_decision, top_count = auto_decisions.most_common(1)[0] if auto_decisions else ("", 0)
    homogeneous = top_count == total

    return {
        "total": total,
        "auto_decision_breakdown": dict(auto_decisions.most_common()),
        "top_decision": top_decision,
        "top_decision_pct": round(100 * top_count / total, 1) if total else 0,
        "homogeneous": homogeneous,
        "avg_confidence": round(sum(confidences) / len(confidences), 3) if confidences else None,
        "scenarios": dict(scenarios.most_common()),
        "namespaces": dict(namespaces.most_common()),
        "sample_rationales": sample_rationales,
    }


# ── Cross-reference enrichment ──
# Map divergence-registry namespace names to the ID prefix used in the
# target graph ``ids`` pipe-list.
_NAMESPACE_TO_ID_PREFIX: dict[str, str] = {
    "Ensembl_Protein_ID": "ENSEMBL:",
    "Ensembl_ID": "ENSEMBL:",
    "Ensembl_Transcript_ID": "ENSEMBL:",
    "RefSeq_Protein_ID": "RefSeq:",
    "RefSeq_Transcript_ID": "RefSeq:",
    "UniProt_ID": "UniProtKB:",
    "NCBI_ID": "NCBIGene:",
    "HGNC_ID": "HGNC:",
    "OMIM_ID": "MIM:",
}


def _lookup_id_in_graph(
    data: TargetGraphData, raw_value: str, namespace: str,
) -> dict[str, Any] | None:
    """Look up a single ID in the target graph.

    Returns ``None`` if the namespace is not mappable (e.g. Description,
    Location, Symbol).  Otherwise returns a dict with ``exists`` bool plus
    target details when found.

    Handles edge cases:
    - Values that already include the prefix (e.g. ``HGNC:6943``)
    - Pipe-delimited compound values (e.g. ``MIM:312095|465000``)
    - Versioned IDs (e.g. ``ENSP00000223321.4``)
    """
    prefix = _NAMESPACE_TO_ID_PREFIX.get(namespace)
    if prefix is None:
        return None
    val = (raw_value or "").strip()
    if not val:
        return None

    # Handle pipe-delimited compound values — split and look up the first
    if "|" in val:
        parts = [p.strip() for p in val.split("|") if p.strip()]
        if len(parts) > 1:
            first_lu = _lookup_id_in_graph(data, parts[0], namespace)
            return first_lu

    def _try_lookup(v: str) -> list[str]:
        """Try prefixed, already-prefixed, and direct lookups."""
        # If value already starts with the expected prefix, use as-is
        if v.startswith(prefix):
            hits = data.ids_to_targets.get(v, [])
            if hits:
                return hits
        # Try adding prefix
        key = f"{prefix}{v}"
        hits = data.ids_to_targets.get(key, [])
        if hits:
            return hits
        # Try direct (value might be indexed without prefix)
        return data.ids_to_targets.get(v, [])

    targets = _try_lookup(val)

    if not targets:
        # Try version-stripped (e.g. ENSP00000223321.4 → ENSP00000223321)
        val_no_ver = re.sub(r"\.\d+$", "", val)
        if val_no_ver != val:
            targets = _try_lookup(val_no_ver)

    if not targets:
        return {"exists": False, "value": raw_value}

    results = []
    for tid in targets[:3]:
        node = data.nodes_by_id.get(tid)
        if node:
            # Parse pipe-delimited ids into a dict: {"UniProtKB": "K7EP69", ...}
            id_map: dict[str, str] = {}
            for xref in _split_pipe(node.get("ids", "")):
                if ":" in xref:
                    pfx, val = xref.split(":", 1)
                    id_map[pfx] = val
            results.append({
                "target_id": tid,
                "primary_id": node.get("primary_id", ""),
                "symbol": node.get("symbol", ""),
                "name": node.get("name", ""),
                "canonical_status": node.get("canonical_status", ""),
                "target_type": node.get("target_type", ""),
                "identifiers": id_map,
            })
    return {"exists": True, "value": raw_value, "targets": results}


def _enrich_conflict_xrefs(
    data: TargetGraphData, row: dict[str, str],
) -> dict[str, Any]:
    """Enrich a conflict row with cross-reference lookups from the graph.

    Called only for the paginated output rows (typically 50) so cost is low.
    Returns a dict keyed by source name with lookup results plus parent-gene
    information when multiple IDs resolve to the same gene.
    """
    namespace = row.get("namespace", "")
    source_value = row.get("source_value", "")
    source_name = (row.get("source", "") or "").lower()
    consensus_value = row.get("consensus_value", "")
    other_sources_raw = row.get("other_sources", "")

    lookups: dict[str, Any] = {}

    # Look up the primary source_value
    if source_value:
        lu = _lookup_id_in_graph(data, source_value, namespace)
        if lu is not None:
            lookups[source_name] = lu

    # Parse other_sources (protein format: "uniprot: ENSP00000455744, ensembl: ...")
    if other_sources_raw:
        for part in re.split(r"[,;]\s*", other_sources_raw):
            match = re.match(r"^\s*(\w+)\s*:\s*(.+?)\s*$", part)
            if match:
                other_src = match.group(1).lower()
                other_val = match.group(2)
                lu = _lookup_id_in_graph(data, other_val, namespace)
                if lu is not None:
                    lookups[other_src] = lu

    # For gene conflicts (no other_sources), look up consensus_value too
    if consensus_value and not other_sources_raw:
        lu = _lookup_id_in_graph(data, consensus_value, namespace)
        if lu is not None:
            lookups["consensus"] = lu

    # Look up the current target in the graph by row_id to show its
    # harmonised identifier set ("what does our graph say about this target?")
    row_id = row.get("row_id", "")
    self_target: dict[str, Any] | None = None
    if row_id:
        # row_id for protein divergence rows is typically a UniProt accession
        for try_ns in ("UniProt_ID", "RefSeq_Protein_ID", "Ensembl_Protein_ID",
                        "NCBI_ID", "HGNC_ID", "Ensembl_ID"):
            lu = _lookup_id_in_graph(data, row_id, try_ns)
            if lu and lu.get("exists") and lu.get("targets"):
                self_target = lu["targets"][0]
                break

    if not lookups and not self_target:
        return {"lookups": lookups}

    # Determine if all found IDs map to the same parent gene
    all_target_ids: set[str] = set()
    for lu in lookups.values():
        if lu.get("exists"):
            for t in lu.get("targets", []):
                all_target_ids.add(t["target_id"])

    result: dict[str, Any] = {"lookups": lookups}
    if self_target:
        result["self_target"] = self_target

    if len(all_target_ids) > 1:
        parent_genes: dict[str, dict[str, str]] = {}
        for tid in all_target_ids:
            for edge in data.relation_edges_by_target.get(tid, []):
                if edge.get("target_id") == tid:
                    parent = edge.get("source_id", "")
                    parent_node = data.nodes_by_id.get(parent)
                    if parent_node and parent_node.get("target_type") == "gene":
                        parent_genes[tid] = {
                            "gene_id": parent,
                            "symbol": parent_node.get("symbol", ""),
                        }
        if parent_genes:
            result["parent_genes"] = parent_genes
            unique_genes = {pg["gene_id"] for pg in parent_genes.values()}
            result["same_gene"] = len(unique_genes) == 1
            if len(unique_genes) == 1:
                gene = next(iter(parent_genes.values()))
                result["shared_gene_symbol"] = gene["symbol"]

    return result


def _compute_xref_suggestion(
    data: TargetGraphData,
    row: dict[str, str],
    xref_context: dict[str, Any],
) -> dict[str, str]:
    """Compute an automated resolution suggestion from xref enrichment.

    Heuristics:
    - One ID obsolete (not in graph), other exists → trust valid source (high)
    - Both exist, same gene, one maps to target matching row_id → trust that source (medium)
    - Both exist, different genes → needs_expert_review (high)
    - Both exist, same gene, can't pick → needs_expert_review (low)
    - No lookups → no suggestion
    """
    lookups = xref_context.get("lookups", {})
    if not lookups:
        return {}

    # Partition into existing vs missing
    existing: dict[str, dict[str, Any]] = {}
    missing: dict[str, dict[str, Any]] = {}
    for src_name, lu in lookups.items():
        if lu.get("exists"):
            existing[src_name] = lu
        else:
            missing[src_name] = lu

    # Scenario: one obsolete, one exists
    if missing and existing and len(existing) == 1:
        valid_src = next(iter(existing))
        decision = _SOURCE_TO_DECISION.get(valid_src, "trust_source")
        obsolete_srcs = ", ".join(missing.keys())
        return {
            "suggested_decision": decision,
            "suggestion_rationale": f"{obsolete_srcs} ID not in graph (likely obsolete); {valid_src} ID is valid",
            "suggestion_confidence": "high",
        }

    # All exist — check gene relationships
    if len(existing) >= 2:
        same_gene = xref_context.get("same_gene")
        row_id = row.get("row_id", "")

        if same_gene is False:
            # Different genes — needs expert
            return {
                "suggested_decision": "needs_expert_review",
                "suggestion_rationale": "IDs map to different genes — requires expert judgment",
                "suggestion_confidence": "high",
            }

        if same_gene is True:
            # Same gene — try to identify which source matches the row target
            for src_name, lu in existing.items():
                for t in lu.get("targets", []):
                    primary_id = t.get("primary_id", "")
                    target_id = t.get("target_id", "")
                    # Protein rows: row_id is UniProt accession, primary_id is UniProtKB:accession
                    if row_id and (row_id in primary_id or row_id in target_id):
                        decision = _SOURCE_TO_DECISION.get(src_name, "trust_source")
                        gene_sym = xref_context.get("shared_gene_symbol", "same gene")
                        return {
                            "suggested_decision": decision,
                            "suggestion_rationale": f"Both map to {gene_sym}; {src_name} matches target {row_id}",
                            "suggestion_confidence": "medium",
                        }

            # Gene rows: check if row_id maps via ENSEMBL: prefix
            if row_id:
                gene_key = f"ENSEMBL:{row_id}"
                gene_targets = data.ids_to_targets.get(gene_key, [])
                if gene_targets:
                    for src_name, lu in existing.items():
                        for t in lu.get("targets", []):
                            if t.get("target_id") in gene_targets:
                                decision = _SOURCE_TO_DECISION.get(src_name, "trust_source")
                                gene_sym = xref_context.get("shared_gene_symbol", "same gene")
                                return {
                                    "suggested_decision": decision,
                                    "suggestion_rationale": f"Both map to {gene_sym}; {src_name} maps to row target",
                                    "suggestion_confidence": "medium",
                                }

            # Same gene but can't determine correct source
            return {
                "suggested_decision": "needs_expert_review",
                "suggestion_rationale": "Both IDs map to same gene but cannot determine correct mapping",
                "suggestion_confidence": "low",
            }

    return {}


TARGET_REVIEW_INTAKE_COLUMNS = [
    "Registry ID", "Entity Type", "Review decision", "Replacement value",
    "Reviewer notes", "Standard Name", "Row ID", "Source", "Namespace",
    "Divergence Type", "Source Value", "Consensus Value", "Scenario ID",
    "Auto Decision", "Auto Confidence", "Auto Rationale",
    "Reviewed by", "Reviewed at", "App review ID",
]


def _read_tsv(path: Path) -> list[dict[str, str]]:
    with open(path, newline="", encoding="utf-8") as fh:
        return list(csv.DictReader(fh, delimiter="\t"))


_ROW_ID_PREFIXES = ("UniProtKB:", "NCBIGene:", "HGNC:", "ENSEMBL:", "RefSeq:", "MIM:")


def _build_divergence_target_index(data: TargetGraphData) -> None:
    """Build mapping from target_id to open divergence rows for conflict flagging."""
    data.divergence_by_target_id = defaultdict(list)
    for rows in (data.divergence_gene, data.divergence_protein, data.divergence_transcript):
        for row in rows:
            if (row.get("status", "") or "open").lower() == "resolved":
                continue
            row_id = row.get("row_id", "")
            if not row_id:
                continue
            matched = False
            for pfx in _ROW_ID_PREFIXES:
                for tid in data.ids_to_targets.get(f"{pfx}{row_id}", []):
                    data.divergence_by_target_id[tid].append(row)
                    matched = True
            if not matched:
                for tid in data.ids_to_targets.get(row_id, []):
                    data.divergence_by_target_id[tid].append(row)


def _load_public_target_qc_bundle(data: TargetGraphData, data_dir: Path) -> None:
    """Load deployable, read-only target QC artifacts.

    Public deployments should not ship the full local review registry or allow
    write-back decisions.  They can still expose reviewed/open counts and a
    compact open-review queue for transparency.
    """
    summary_path = data_dir / "target_qc_public_summary.json"
    rows_path = data_dir / "target_review_public.tsv"

    if summary_path.exists():
        try:
            with open(summary_path, encoding="utf-8") as fh:
                parsed = json.load(fh)
            if isinstance(parsed, dict):
                data.public_qc_summary = parsed
        except (OSError, json.JSONDecodeError) as exc:
            print(f"  Warning: could not read public target QC summary: {exc}")

    if rows_path.exists():
        try:
            rows = _read_tsv(rows_path)
        except OSError as exc:
            print(f"  Warning: could not read public target review rows: {exc}")
            rows = []
        for row in rows:
            etype = (row.get("entity_type") or "").strip().lower()
            if etype == "gene":
                data.divergence_gene.append(row)
            elif etype == "protein":
                data.divergence_protein.append(row)
            elif etype == "transcript":
                data.divergence_transcript.append(row)

    if data.public_qc_summary or data.divergence_gene or data.divergence_protein or data.divergence_transcript:
        data.review_mode = "public"
        data.review_can_write = False
        for etype, rows in [
            ("gene", data.divergence_gene),
            ("protein", data.divergence_protein),
            ("transcript", data.divergence_transcript),
        ]:
            for row in rows:
                if not row.get("status"):
                    row["status"] = "open"
                row["_triage_category"] = classify_triage_category(row)
                row["_review_group"] = classify_review_group(etype, row)
        _build_divergence_target_index(data)
        print(
            "  Loaded public target QC bundle: "
            f"{len(data.divergence_gene) + len(data.divergence_protein) + len(data.divergence_transcript):,} open rows"
        )


def _norm(value: Any) -> str:
    text = "" if value is None else str(value)
    text = text.lower()
    text = re.sub(r"['`’]", "", text)
    text = re.sub(r"[^a-z0-9]+", " ", text)
    return re.sub(r"\s+", " ", text).strip()


def _tokens(value: Any) -> set[str]:
    return {token for token in _norm(value).split() if len(token) > 1}


def _split_pipe(value: Any) -> list[str]:
    text = "" if value is None else str(value).strip()
    if not text or text.lower() in {"nan", "none", "null"}:
        return []
    return [part.strip() for part in text.split("|") if part.strip()]


def _index_exact(data: TargetGraphData, term: str, target_id: str) -> None:
    norm = _norm(term)
    if norm:
        data.exact_terms[norm].add(target_id)


def _index_tokens(data: TargetGraphData, term: str, target_id: str) -> None:
    for token in _tokens(term):
        data.token_index[token].add(target_id)


def _index_row(data: TargetGraphData, row: dict[str, str]) -> None:
    target_id = row.get("target_id", "")
    if not target_id:
        return
    data.nodes.append(row)
    data.nodes_by_id[target_id] = row
    data.type_index[row.get("target_type", "unknown")].append(target_id)

    for term in (
        target_id,
        row.get("primary_id", ""),
        row.get("symbol", ""),
        row.get("name", ""),
    ):
        _index_exact(data, term, target_id)
    for term in (row.get("symbol", ""), row.get("name", "")):
        _index_tokens(data, term, target_id)
    for alias in _split_pipe(row.get("description", "")):
        if len(alias) <= 90:
            _index_tokens(data, alias, target_id)
    for xref in _split_pipe(row.get("ids", "")):
        data.ids_to_targets[xref].append(target_id)
        _index_exact(data, xref, target_id)
        # Also index the bare value (without prefix) so users can search
        # for "Q6N017" instead of needing "UniProtKB:Q6N017"
        if ":" in xref:
            bare = xref.split(":", 1)[1]
            _index_exact(data, bare, target_id)


def _index_edge(data: TargetGraphData, row: dict[str, str]) -> None:
    source_id = row.get("source_id", "")
    target_id = row.get("target_id", "")
    predicate = row.get("predicate", "")
    if not source_id or not target_id:
        return
    key = (source_id, predicate, target_id)
    existing = data.relation_edge_index.get(key)
    if existing is not None:
        _merge_relation_edge(existing, row)
        return
    data.relation_edge_index[key] = row
    data.relation_edges.append(row)
    data.relation_edges_by_target[source_id].append(row)
    data.relation_edges_by_target[target_id].append(row)


def _merge_relation_edge(existing: dict[str, str], incoming: dict[str, str]) -> None:
    for field in ("relation_kind", "evidence_identifier", "evidence_namespace", "evidence_source", "support_tier"):
        existing[field] = _pipe_union(existing.get(field, ""), incoming.get(field, ""))
    existing["support_ratio"] = _max_numeric_text(existing.get("support_ratio", ""), incoming.get("support_ratio", ""))
    existing["support_score"] = _max_numeric_text(existing.get("support_score", ""), incoming.get("support_score", ""))


def _pipe_union(*values: Any) -> str:
    output: list[str] = []
    for value in values:
        for part in re.split(r"[|;]", "" if value is None else str(value)):
            part = part.strip()
            if part and part not in output:
                output.append(part)
    return "|".join(output)


def _max_numeric_text(existing: str, incoming: str) -> str:
    existing = str(existing or "").strip()
    incoming = str(incoming or "").strip()
    if not existing:
        return incoming
    if not incoming:
        return existing
    try:
        return existing if float(existing) >= float(incoming) else incoming
    except ValueError:
        return _pipe_union(existing, incoming)


def load_target_graph_data(
    data_dir: str | Path,
    qc_dir: str | Path = "",
) -> TargetGraphData:
    global _singleton
    data_dir = Path(data_dir)
    if _singleton is not None:
        return _singleton
    with _singleton_lock:
        if _singleton is not None:
            return _singleton
        if not data_dir.is_dir():
            raise HTTPException(status_code=500, detail=f"Target app graph dir not found: {data_dir}")
        nodes_path = data_dir / "target_nodes.tsv"
        if not nodes_path.exists():
            raise HTTPException(status_code=500, detail=f"Target nodes file not found: {nodes_path}")

        data = TargetGraphData()
        manifest_path = data_dir / "manifest.json"
        if manifest_path.exists():
            with open(manifest_path, encoding="utf-8") as fh:
                data.manifest = json.load(fh)
        with open(nodes_path, newline="", encoding="utf-8") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                _index_row(data, row)
        edges_path = data_dir / "target_edges.tsv"
        if edges_path.exists():
            with open(edges_path, newline="", encoding="utf-8") as fh:
                for row in csv.DictReader(fh, delimiter="\t"):
                    _index_edge(data, row)

        divergence_files = [
            ("gene_divergence_registry.tsv", "divergence_gene"),
            ("protein_divergence_registry.tsv", "divergence_protein"),
            ("transcript_divergence_registry.tsv", "divergence_transcript"),
        ]

        def _usable_qc_dir(path: Path | None) -> bool:
            return bool(path and path.is_dir() and any((path / filename).exists() for filename, _ in divergence_files))

        # Resolve QC directory: explicit > manifest-inferred > relative to data_dir.
        # Use None as the empty sentinel; Path() points at cwd and would falsely
        # enable internal review mode in deployed bundles.
        resolved_qc: Path | None = Path(qc_dir) if qc_dir else None
        if not _usable_qc_dir(resolved_qc):
            # Try to infer from manifest source_files (paths relative to pipeline root)
            source_files = data.manifest.get("source_files", {})
            for rel_path in source_files.values():
                if "target_data" in rel_path:
                    # Walk up from data_dir to find pipeline root
                    for ancestor in [data_dir] + list(data_dir.parents):
                        candidate = ancestor / rel_path
                        if candidate.exists():
                            # Found pipeline root — qc is at src/data/publicdata/target_data/qc
                            td_idx = rel_path.find("target_data/")
                            if td_idx >= 0:
                                resolved_qc = ancestor / rel_path[:td_idx + len("target_data")] / "qc"
                                if _usable_qc_dir(resolved_qc):
                                    break
                    if _usable_qc_dir(resolved_qc):
                        break
        if not _usable_qc_dir(resolved_qc):
            # Fallback: relative to data_dir
            for candidate in [
                data_dir.parent / "qc",
                data_dir / "qc",
            ]:
                if _usable_qc_dir(candidate):
                    resolved_qc = candidate
                    break
        data.qc_dir = str(resolved_qc) if _usable_qc_dir(resolved_qc) else ""

        # Load divergence registries from qc/ directory
        if data.qc_dir:
            data.review_mode = "internal"
            data.review_can_write = True
            qc_path = Path(data.qc_dir)
            for filename, attr in divergence_files:
                path = qc_path / filename
                if path.exists():
                    rows = _read_tsv(path)
                    setattr(data, attr, rows)
                    if rows:
                        print(f"  Loaded {len(rows):,} {attr.replace('divergence_', '')} divergence rows")

            # Pre-compute triage and review groups on every divergence row (avoids
            # recomputing on each API request across 847K rows).
            for etype, rows in [
                ("gene", data.divergence_gene),
                ("protein", data.divergence_protein),
                ("transcript", data.divergence_transcript),
            ]:
                for row in rows:
                    row["_triage_category"] = classify_triage_category(row)
                    row["_review_group"] = classify_review_group(etype, row)

            # Reconcile previously-saved review decisions so resolved rows
            # don't reappear as "open" after a server restart.
            _reconcile_prior_decisions(data, qc_path)

            # Load gene biotype from provenance CSV if available
            _load_gene_biotype(data, qc_path)

            # Build target → divergence index for conflict edge flagging
            _build_divergence_target_index(data)
        else:
            _load_public_target_qc_bundle(data, data_dir)
            _build_divergence_target_index(data)

        _singleton = data
        print(f"Target graph loaded: {len(data.nodes):,} targets, {len(data.relation_edges):,} relation edges")
        if data.qc_dir:
            print(f"  QC dir: {data.qc_dir}")
        if data.review_mode != "none":
            print(f"  Target review mode: {data.review_mode} (write={data.review_can_write})")
        return data


def _load_gene_biotype(data: TargetGraphData, qc_dir: Path) -> None:
    """Load gene biotype distribution from gene_mapping_provenance.csv."""
    # Look for provenance CSV relative to qc_dir (sibling cleaned/sources/)
    candidates = [
        qc_dir.parent / "cleaned" / "sources" / "gene_mapping_provenance.csv",
    ]
    for path in candidates:
        if not path.exists():
            continue
        biotype: Counter[str] = Counter()
        with open(path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh)
            for row in reader:
                # Prefer Ensembl biotype (most granular), fall back to NCBI, then HGNC
                bt = (
                    row.get("ensembl_gene_type", "").strip()
                    or row.get("ncbi_gene_type", "").strip()
                    or row.get("hgnc_gene_type", "").strip()
                )
                if bt:
                    biotype[bt] += 1
        data.gene_biotype_counts = dict(biotype.most_common())
        if biotype:
            print(f"  Loaded gene biotype distribution: {len(biotype)} types, {sum(biotype.values()):,} genes")
        return


def _reconcile_prior_decisions(data: TargetGraphData, qc_dir: Path) -> None:
    """Read the review decisions TSV and mark matching rows as resolved.

    This ensures previously-saved decisions survive a server restart.
    Only reads the 'Registry ID' column to build a lookup set, so the
    cost is proportional to the decisions file size, not divergence rows.
    """
    decisions_path = qc_dir / "review" / "app_completed" / "target_app_review_decisions.tsv"
    if not decisions_path.exists():
        return
    resolved_ids: set[str] = set()
    try:
        with open(decisions_path, newline="", encoding="utf-8") as fh:
            reader = csv.DictReader(fh, delimiter="\t")
            for row in reader:
                rid = row.get("Registry ID", "").strip()
                if rid:
                    resolved_ids.add(rid)
    except Exception as exc:
        print(f"  Warning: could not read prior decisions: {exc}")
        return
    if not resolved_ids:
        return
    count = 0
    for rows in (data.divergence_gene, data.divergence_protein, data.divergence_transcript):
        for row in rows:
            if row.get("registry_id", "") in resolved_ids:
                row["status"] = "resolved"
                count += 1
    print(f"  Reconciled {count:,} previously-resolved rows from {len(resolved_ids):,} decisions")


def compute_divergence_stats(data: TargetGraphData) -> dict[str, Any]:
    """Single-pass divergence statistics across gene + protein + transcript."""
    if data._divergence_stats is not None:
        return data._divergence_stats
    if data.review_mode == "public" and isinstance(data.public_qc_summary.get("divergence"), dict):
        result = dict(data.public_qc_summary["divergence"])
        result.setdefault("review_groups", REVIEW_GROUPS)
        result.setdefault("review_mode", data.review_mode)
        result.setdefault("review_can_write", data.review_can_write)
        data._divergence_stats = result
        return result

    entity_counts: dict[str, int] = {}
    type_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    auto_decision_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    namespace_counts: Counter[str] = Counter()
    triage_counts: Counter[str] = Counter()
    review_group_counts: Counter[str] = Counter()
    review_group_open_counts: Counter[str] = Counter()
    review_group_resolved_counts: Counter[str] = Counter()
    confidence_values: list[float] = []
    source_breakdown: dict[str, dict[str, Any]] = {}
    open_conflicts = 0
    actionable_count = 0
    total = 0
    # Track open conflict auto-decision breakdown for priority guidance
    open_conflict_auto_decisions: Counter[str] = Counter()

    for entity_type, rows in [
        ("gene", data.divergence_gene),
        ("protein", data.divergence_protein),
        ("transcript", data.divergence_transcript),
    ]:
        entity_counts[entity_type] = len(rows)
        for row in rows:
            total += 1
            div_type = row.get("divergence_type", "") or ""
            status = row.get("status", "") or "open"
            auto_dec = row.get("auto_decision", "") or ""
            source = row.get("source", "") or ""
            ns = row.get("namespace", "") or ""

            type_counts[div_type] += 1
            status_counts[status] += 1
            if auto_dec:
                auto_decision_counts[auto_dec] += 1
            if source:
                source_counts[source] += 1
            if ns:
                namespace_counts[ns] += 1

            conf_str = row.get("auto_confidence", "")
            try:
                conf = float(conf_str)
                confidence_values.append(conf)
            except (ValueError, TypeError):
                pass

            triage_counts[classify_triage_category(row)] += 1
            review_group = classify_review_group(entity_type, row)
            review_group_counts[review_group] += 1
            if status == "resolved":
                review_group_resolved_counts[review_group] += 1
            else:
                review_group_open_counts[review_group] += 1

            if div_type == "conflict" and status == "open":
                open_conflicts += 1
                open_conflict_auto_decisions[auto_dec or "no_decision"] += 1
            if auto_dec == "needs_expert_review":
                actionable_count += 1

            # Per-source breakdown
            if source:
                bd = source_breakdown.setdefault(source, {
                    "total": 0, "conflict_count": 0, "sole_source_count": 0,
                    "low_score_count": 0, "confidence_sum": 0.0,
                    "confidence_n": 0, "needs_expert_count": 0,
                })
                bd["total"] += 1
                if div_type == "conflict":
                    bd["conflict_count"] += 1
                elif div_type == "sole_source":
                    bd["sole_source_count"] += 1
                elif div_type == "low_score":
                    bd["low_score_count"] += 1
                try:
                    bd["confidence_sum"] += float(conf_str)
                    bd["confidence_n"] += 1
                except (ValueError, TypeError):
                    pass
                if auto_dec == "needs_expert_review":
                    bd["needs_expert_count"] += 1

    # Build confidence histogram (0.0–1.0 in 10 buckets)
    confidence_distribution: dict[str, int] = {}
    if confidence_values:
        for v in confidence_values:
            bucket = min(int(v * 10), 9)
            label = f"{bucket * 10}-{(bucket + 1) * 10}%"
            confidence_distribution[label] = confidence_distribution.get(label, 0) + 1

    # Finalize source breakdown with avg_confidence
    for bd in source_breakdown.values():
        n = bd.pop("confidence_n")
        s = bd.pop("confidence_sum")
        bd["avg_confidence"] = round(s / n, 3) if n else None

    result: dict[str, Any] = {
        "total_divergences": total,
        "divergences_by_entity": entity_counts,
        "divergence_type_counts": dict(type_counts.most_common()),
        "status_counts": dict(status_counts.most_common()),
        "auto_decision_counts": dict(auto_decision_counts.most_common()),
        "source_divergence_counts": dict(source_counts.most_common()),
        "namespace_divergence_counts": dict(namespace_counts.most_common()),
        "confidence_distribution": confidence_distribution,
        "source_divergence_breakdown": dict(
            sorted(source_breakdown.items(), key=lambda kv: -kv[1]["total"])
        ),
        "open_conflicts": open_conflicts,
        "actionable_count": actionable_count,
        "open_conflict_auto_decisions": dict(open_conflict_auto_decisions.most_common()),
        "triage_category_counts": dict(
            sorted(triage_counts.items(), key=lambda kv: _TRIAGE_ORDER.get(kv[0], 9))
        ),
        "review_group_counts": dict(
            sorted(review_group_counts.items(), key=lambda kv: _REVIEW_GROUP_ORDER.get(kv[0], 99))
        ),
        "review_group_open_counts": dict(
            sorted(review_group_open_counts.items(), key=lambda kv: _REVIEW_GROUP_ORDER.get(kv[0], 99))
        ),
        "review_group_resolved_counts": dict(
            sorted(review_group_resolved_counts.items(), key=lambda kv: _REVIEW_GROUP_ORDER.get(kv[0], 99))
        ),
        "review_groups": REVIEW_GROUPS,
        "review_mode": data.review_mode,
        "review_can_write": data.review_can_write,
    }
    data._divergence_stats = result
    return result


def compute_target_stats(data: TargetGraphData) -> dict[str, Any]:
    if data._stats is not None:
        # Divergence stats may have been invalidated independently
        if "divergence" not in data._stats:
            data._stats["divergence"] = compute_divergence_stats(data)
        return data._stats
    type_counts = Counter(row.get("target_type", "unknown") or "unknown" for row in data.nodes)
    namespace_counts: Counter[str] = Counter()
    canonical_counts: Counter[str] = Counter()
    # Sources per entity type (excluding validation-only sources like nodenorm)
    gene_sources: Counter[str] = Counter()
    protein_sources: Counter[str] = Counter()
    transcript_sources: Counter[str] = Counter()
    canonical_count = 0
    annotation_counts: Counter[str] = Counter(data.manifest.get("annotation_counts") or {})
    has_manifest_annotation_counts = bool(annotation_counts)
    for row in data.nodes:
        target_type = row.get("target_type", "")
        for ns in _split_pipe(row.get("id_namespaces", "")):
            namespace_counts[ns] += 1
        # Canonical status only tracked for proteins
        if target_type == "protein":
            status = row.get("canonical_status", "") or "not_applicable"
            canonical_counts[status] += 1
        if _canonical_rank(row) == 0:
            canonical_count += 1
        for src in _split_pipe(row.get("source_namespaces", "").replace(",", "|")):
            if src.lower() in _VALIDATION_SOURCES:
                continue
            if target_type == "gene":
                gene_sources[src] += 1
            elif target_type == "protein":
                protein_sources[src] += 1
            elif target_type == "transcript":
                transcript_sources[src] += 1
        if not has_manifest_annotation_counts:
            annotations = _parse_annotations(row)
            annotation_counts.update(annotations.keys())
    edge_counts = Counter(data.manifest.get("edge_predicate_counts") or {})
    if not edge_counts:
        edge_counts = Counter(row.get("predicate", "unknown") or "unknown" for row in data.relation_edges)
    # All real sources combined (for "sources covered" stat)
    all_sources: set[str] = set(gene_sources) | set(protein_sources) | set(transcript_sources)
    stats = {
        "total_targets": len(data.nodes),
        "type_counts": dict(type_counts.most_common()),
        "identifier_namespace_counts": dict(namespace_counts.most_common()),
        "canonical_status_counts": dict(canonical_counts.most_common(20)),
        "gene_sources": dict(gene_sources.most_common()),
        "protein_sources": dict(protein_sources.most_common()),
        "transcript_sources": dict(transcript_sources.most_common()),
        "edge_predicate_counts": dict(edge_counts.most_common()),
        "annotation_counts": dict(annotation_counts.most_common()),
        "download_column_groups": _download_column_groups(annotation_counts),
        "manifest": data.manifest,
        "canonical_count": canonical_count,
        "sources_covered": len(all_sources),
        "gene_biotype_counts": data.gene_biotype_counts,
        "divergence": compute_divergence_stats(data),
        "review_mode": data.review_mode,
        "review_can_write": data.review_can_write,
    }
    data._stats = stats
    return stats


def _resolve_target_id(data: TargetGraphData, value: str) -> str | None:
    value = (value or "").strip()
    if value in data.nodes_by_id:
        return value
    candidates = data.ids_to_targets.get(value)
    if candidates:
        return sorted(candidates, key=lambda tid: _target_id_rank(data, tid))[0]
    norm = _norm(value)
    exact = data.exact_terms.get(norm)
    if exact:
        return sorted(exact, key=lambda tid: _target_id_rank(data, tid))[0]
    return None


def _target_score(query_norm: str, row: dict[str, str]) -> float:
    values = [row.get("symbol", ""), row.get("name", ""), row.get("primary_id", ""), row.get("target_id", "")]
    norms = [_norm(v) for v in values if v]
    if query_norm in norms:
        return 1.0
    best = max((SequenceMatcher(None, query_norm, n).ratio() for n in norms), default=0.0)
    if any(n.startswith(query_norm) for n in norms if query_norm):
        best = max(best, 0.92)
    return best


def _type_rank(target_type: str) -> int:
    return {"gene": 0, "protein": 1, "transcript": 2}.get(target_type, 9)


def _canonical_rank(row: dict[str, str]) -> int:
    status = (row.get("canonical_status", "") or "").lower()
    category = (row.get("category", "") or "").lower()
    quality = (row.get("quality_note", "") or "").lower()
    if "true|canonical" in status or category == "canonical" or quality == "reviewed_representative":
        return 0
    if "canonical" in status:
        return 1
    if "alternate" in category:
        return 3
    return 2


def _target_id_rank(data: TargetGraphData, target_id: str) -> tuple[int, int, str]:
    row = data.nodes_by_id.get(target_id, {})
    return (
        _type_rank(row.get("target_type", "")),
        _canonical_rank(row),
        row.get("symbol", "") or row.get("name", "") or target_id,
    )


def _parse_annotations(row: dict[str, str]) -> dict[str, Any]:
    raw = row.get("source_annotations", "")
    if not raw:
        return {}
    try:
        parsed = json.loads(raw)
    except json.JSONDecodeError:
        return {}
    return parsed if isinstance(parsed, dict) else {}


def _target_to_row(row: dict[str, str]) -> dict[str, Any]:
    return {
        "target_id": row.get("target_id", ""),
        "target_type": row.get("target_type", ""),
        "primary_id": row.get("primary_id", ""),
        "symbol": row.get("symbol", ""),
        "name": row.get("name", ""),
        "description": row.get("description", ""),
        "category": row.get("category", ""),
        "ids": row.get("ids", ""),
        "id_namespaces": row.get("id_namespaces", ""),
        "source_namespaces": row.get("source_namespaces", ""),
        "mapping_ratio": row.get("mapping_ratio", ""),
        "quality_note": row.get("quality_note", ""),
        "canonical_status": row.get("canonical_status", ""),
        "canonical_ifx_id": row.get("canonical_ifx_id", ""),
        "source_annotations": row.get("source_annotations", ""),
        "createdAt": row.get("createdAt", ""),
        "updatedAt": row.get("updatedAt", ""),
    }


def _edge_classes(predicate: str) -> str:
    if predicate == "biolink:transcribed_to":
        return "target-relation-edge transcription-edge"
    if predicate == "biolink:translates_to":
        return "target-relation-edge translation-edge"
    if predicate == "biolink:has_gene_product":
        return "target-relation-edge gene-product-edge"
    return "target-relation-edge related-edge"


def _edge_relationship_metadata(predicate: str) -> dict[str, str]:
    if predicate == "biolink:transcribed_to":
        return {
            "display_label": "transcribed_to",
            "relationship_label": "Gene -> Transcript",
            "relationship_description": "Gene is transcribed to an RNA transcript. Identifier evidence usually comes from Ensembl or RefSeq transcript mappings.",
            "relationship_category": "Transcript",
        }
    if predicate == "biolink:translates_to":
        return {
            "display_label": "translates_to",
            "relationship_label": "Transcript -> Protein",
            "relationship_description": "Transcript is translated into a protein or protein isoform.",
            "relationship_category": "Translates to",
        }
    if predicate == "biolink:has_gene_product":
        return {
            "display_label": "has_gene_product",
            "relationship_label": "Gene -> Protein",
            "relationship_description": "Gene has a protein product. Protein identifiers are UniProt-centered, with Ensembl/RefSeq/NodeNorm evidence retained as provenance.",
            "relationship_category": "Protein",
        }
    return {
        "display_label": predicate.replace("biolink:", "") or "related_to",
        "relationship_label": "Related target",
        "relationship_description": "Target relationship retained from the IFX Harmonizers target workflow.",
        "relationship_category": "Related",
    }


def _build_edge_element(
    source_id: str,
    target_id: str,
    edge_row: dict[str, str],
) -> dict[str, Any]:
    predicate = edge_row.get("predicate", "") or "biolink:related_to"
    evidence_key = edge_row.get("evidence_identifier", "")
    edge_id = f"{source_id}->{target_id}:{predicate}:{evidence_key}"
    relationship = _edge_relationship_metadata(predicate)
    return {
        "data": {
            "id": edge_id,
            "source": source_id,
            "target": target_id,
            "label": relationship["display_label"],
            "kind": "target_relation_edge",
            "predicate": predicate,
            "relationship_label": relationship["relationship_label"],
            "relationship_description": relationship["relationship_description"],
            "relationship_category": relationship["relationship_category"],
            "relation_kind": edge_row.get("relation_kind", ""),
            "evidence_identifier": evidence_key,
            "evidence_namespace": edge_row.get("evidence_namespace", ""),
            "evidence_source": edge_row.get("evidence_source", ""),
            "support_tier": edge_row.get("support_tier", ""),
            "support_ratio": edge_row.get("support_ratio", ""),
            "support_score": edge_row.get("support_score", ""),
        },
        "classes": _edge_classes(predicate),
    }


def _annotate_conflict_edge(
    data: TargetGraphData,
    element: dict[str, Any],
) -> None:
    """Annotate an edge element with open-review info from either endpoint."""
    payload = element.get("data", {})
    endpoint_ids = [payload.get("source", ""), payload.get("target", "")]
    div_rows: list[dict[str, str]] = []
    seen_registry_ids: set[str] = set()
    for node_id in endpoint_ids:
        for row in data.divergence_by_target_id.get(node_id, []):
            registry_id = row.get("registry_id", "")
            if registry_id and registry_id in seen_registry_ids:
                continue
            if registry_id:
                seen_registry_ids.add(registry_id)
            div_rows.append(row)
    if not div_rows:
        return
    element["classes"] += " has-conflict unresolved-qc-edge"
    element["data"]["has_conflict"] = True
    element["data"]["review_status"] = "unresolved"
    element["data"]["qc_indicator"] = "Unresolved QC finding"
    element["data"]["conflict_count"] = len(div_rows)
    decisions = Counter(r.get("auto_decision", "") for r in div_rows if r.get("auto_decision"))
    element["data"]["conflict_decisions"] = ", ".join(f"{d} ({c})" for d, c in decisions.most_common(3))
    namespaces = Counter(r.get("namespace", "") for r in div_rows if r.get("namespace"))
    sources = Counter(r.get("source", "") for r in div_rows if r.get("source"))
    groups = Counter((r.get("review_group") or r.get("_review_group") or "") for r in div_rows if (r.get("review_group") or r.get("_review_group")))
    element["data"]["conflict_namespaces"] = ", ".join(f"{k} ({v})" for k, v in namespaces.most_common(4))
    element["data"]["conflict_sources"] = ", ".join(f"{k} ({v})" for k, v in sources.most_common(4))
    element["data"]["conflict_review_groups"] = ", ".join(f"{k} ({v})" for k, v in groups.most_common(4))
    sample = next((r.get("auto_rationale", "") for r in div_rows if r.get("auto_rationale")), "")
    element["data"]["conflict_rationale"] = sample


def _annotation_group(field: str) -> str:
    lower = field.lower()
    if lower.startswith("hgnc_") or lower.startswith("ncbi_") or lower in {"consolidated_location", "ensembl_strand"}:
        return "Gene annotations"
    if lower.startswith("uniprot") or lower.startswith("protein_") or lower.startswith("mapping_"):
        return "Protein annotations"
    if lower.startswith("ensembl_trans") or lower.startswith("refseq_"):
        return "Transcript annotations"
    if lower.endswith("_provenance") or lower.endswith("_mapping_score") or lower.startswith("total_mapping"):
        return "Provenance and support"
    return "Source annotations"


def _download_column_groups(annotation_counts: Counter[str]) -> list[dict[str, Any]]:
    groups: dict[str, list[str]] = {
        "Core": [
            "target_id", "target_type", "primary_id", "symbol", "name", "description", "category",
        ],
        "Identifiers": [
            "ids", "id_namespaces", "source_namespaces",
        ],
        "Quality and status": [
            "mapping_ratio", "quality_note", "canonical_status", "canonical_ifx_id", "createdAt", "updatedAt",
        ],
    }
    for field, _count in annotation_counts.most_common():
        groups.setdefault(_annotation_group(field), []).append(field)
    return [{"label": label, "columns": columns} for label, columns in groups.items()]


def _filter_targets(data: TargetGraphData, q: str = "", target_type: str = "", namespace: str = "") -> list[dict[str, str]]:
    q_norm = _norm(q)
    target_type = (target_type or "").strip().lower()
    namespace = (namespace or "").strip()

    if q_norm:
        candidate_ids: set[str] = set(data.exact_terms.get(q_norm, set()))
        for token in _tokens(q):
            ids = data.token_index.get(token, set())
            if len(ids) <= 5000:
                candidate_ids.update(ids)
        candidates = [data.nodes_by_id[tid] for tid in candidate_ids if tid in data.nodes_by_id]
    else:
        candidates = data.nodes

    filtered: list[dict[str, str]] = []
    for row in candidates:
        if target_type and row.get("target_type", "").lower() != target_type:
            continue
        if namespace and namespace not in _split_pipe(row.get("id_namespaces", "")):
            continue
        filtered.append(row)

    sort_key = (
        (lambda row: (
            -_target_score(q_norm, row),
            _type_rank(row.get("target_type", "")),
            _canonical_rank(row),
            row.get("symbol", "") or row.get("name", ""),
        ))
        if q_norm else
        (lambda row: (
            _type_rank(row.get("target_type", "")),
            _canonical_rank(row),
            row.get("symbol", "") or row.get("name", ""),
        ))
    )
    filtered.sort(key=sort_key)
    return filtered


def search_targets(
    data: TargetGraphData,
    q: str = "",
    target_type: str = "",
    namespace: str = "",
    page: int = 1,
    per_page: int = 50,
) -> dict[str, Any]:
    q_norm = _norm(q)
    target_type = (target_type or "").strip().lower()
    namespace = (namespace or "").strip()
    per_page = max(1, min(per_page, 200))

    filtered = _filter_targets(data, q=q, target_type=target_type, namespace=namespace)

    total = len(filtered)
    total_pages = max(1, (total + per_page - 1) // per_page)
    page = max(1, min(page, total_pages))
    start = (page - 1) * per_page
    rows = [_target_to_row(row) for row in filtered[start:start + per_page]]
    return {"rows": rows, "total": total, "page": page, "per_page": per_page, "total_pages": total_pages}


def _resolver_norm(value: Any) -> str:
    text = "" if value is None else str(value)
    text = text.replace("α", " alpha ").replace("β", " beta ").replace("γ", " gamma ").replace("δ", " delta ")
    text = text.replace("Α", " alpha ").replace("Β", " beta ").replace("Γ", " gamma ").replace("Δ", " delta ")
    return _norm(text)


def _resolver_tokens(value: Any) -> set[str]:
    return {token for token in _resolver_norm(value).split() if len(token) > 1}


def _target_query_type(query: str) -> dict[str, Any]:
    """Classify free text before matching so non-target labels are not forced."""
    norm = _resolver_norm(query)
    tokens = set(norm.split())
    warnings: list[str] = []
    query_type = "gene_or_protein_name"
    authoritative = True

    id_patterns = [
        r"ifx(gene|protein|transcript):[a-z0-9]+",
        r"uniprotkb:[a-z0-9]+",
        r"hgnc:\d+",
        r"ncbigene:\d+",
        r"ensembl:(ensg|ensp|enst)\d+",
        r"refseq:(nm|nr|np|xp|xm)_\d+(\.\d+)?",
        r"(ensg|ensp|enst)\d+",
        r"[opq][0-9][a-z0-9]{3}[0-9](-\d+)?",
        r"[a-nr-z][0-9]([a-z][a-z0-9]{2}[0-9]){1,2}(-\d+)?",
        r"(nm|nr|np|xp|xm)_\d+(\.\d+)?",
    ]
    if any(re.fullmatch(pattern, norm) for pattern in id_patterns):
        return {
            "query_type": "identifier",
            "authoritative_mapping_allowed": True,
            "warnings": [],
        }

    if "/" in str(query) or ";" in str(query) or " fusion " in f" {norm} ":
        query_type = "protein_complex_or_multi_target"
        authoritative = False
        warnings.append("multi_target_or_fusion_label")

    process_terms = {
        "apoptosis", "differentiation", "biosynthetic", "process",
        "pathway", "response", "stress", "alkylation", "signaling",
    }
    if norm == "cell growth" or tokens & process_terms:
        query_type = "biological_process_or_pathway"
        authoritative = False
        warnings.append("not_a_single_target_label")

    if "cell line" in norm or norm in {"hepg2", "ht 1080", "a673"} or tokens & {"cells", "leukemia"}:
        query_type = "cell_line_or_cell_context"
        authoritative = False
        warnings.append("cell_line_or_context_not_target")

    if tokens & {"trypanosoma", "brucei"}:
        query_type = "organism_or_pathogen"
        authoritative = False
        warnings.append("organism_not_human_target")

    if norm in {"dna", "rna"}:
        query_type = "molecule_or_macromolecule"
        authoritative = False
        warnings.append("macromolecule_not_specific_gene_or_protein")

    broad_family_labels = {
        "tubulin",
        "d2 like dopamine receptor",
        "retinoic acid receptor",
        "serine threonine protein kinase raf",
        "neurotrophic tyrosine kinase receptor",
        "dual specificity phosphatases",
        "dual specific phosphatases",
    }
    if authoritative and (
        norm in broad_family_labels
        or "phosphatases" in tokens
        or ("like" in tokens and "receptor" in tokens)
        or norm.endswith(" kinase raf")
    ):
        query_type = "protein_family_or_multi_target"
        authoritative = False
        warnings.append("family_label_may_have_multiple_targets")

    return {
        "query_type": query_type,
        "authoritative_mapping_allowed": authoritative,
        "warnings": sorted(set(warnings)),
    }


def _target_term_priority(field: str, row: dict[str, str]) -> float:
    base = {
        "target_id": 1.0,
        "primary_id": 0.99,
        "external_id": 0.98,
        "symbol": 0.97,
        "name": 0.92,
        "description": 0.74,
    }.get(field, 0.82)
    target_type = (row.get("target_type") or "").lower()
    if target_type == "gene" and field == "symbol":
        base += 0.01
    if target_type == "protein" and field == "name":
        base += 0.01
    return min(base, 1.0)


def _add_target_resolver_term(
    terms: list[dict[str, Any]],
    seen: set[tuple[str, str, str]],
    row: dict[str, str],
    term: str,
    field: str,
) -> None:
    term = str(term or "").strip()
    if not term or term.lower() in {"nan", "none", "null"}:
        return
    norm = _resolver_norm(term)
    if not norm:
        return
    key = (row.get("target_id", ""), norm, field)
    if key in seen:
        return
    seen.add(key)
    terms.append({
        "target_id": row.get("target_id", ""),
        "term": term,
        "norm": norm,
        "field": field,
    })


def _target_row_resolver_terms(row: dict[str, str]) -> list[dict[str, Any]]:
    terms: list[dict[str, Any]] = []
    seen: set[tuple[str, str, str]] = set()
    for field in ("target_id", "primary_id", "symbol", "name"):
        _add_target_resolver_term(terms, seen, row, row.get(field, ""), field)
    for alias in _split_pipe(row.get("description", "")):
        if len(alias) <= 90:
            _add_target_resolver_term(terms, seen, row, alias, "description")
    for xref in _split_pipe(row.get("ids", "")):
        _add_target_resolver_term(terms, seen, row, xref, "external_id")
        if ":" in xref:
            _add_target_resolver_term(terms, seen, row, xref.split(":", 1)[1], "external_id")
    return terms


def _target_lexical_score(query_norm: str, term_norm: str) -> float:
    if not query_norm or not term_norm:
        return 0.0
    if query_norm == term_norm:
        return 1.0
    query_tokens = _resolver_tokens(query_norm)
    term_tokens = _resolver_tokens(term_norm)
    if not query_tokens or not term_tokens:
        return 0.0
    overlap = len(query_tokens & term_tokens)
    containment = overlap / max(1, min(len(query_tokens), len(term_tokens)))
    jaccard = overlap / max(1, len(query_tokens | term_tokens))
    sequence = SequenceMatcher(None, query_norm, term_norm).ratio()
    substring_score = 0.0
    if query_norm in term_norm or term_norm in query_norm:
        shorter = min(len(query_norm), len(term_norm))
        longer = max(len(query_norm), len(term_norm))
        substring_score = 0.78 + 0.18 * (shorter / max(1, longer))
    return round(max(sequence, containment * 0.95, jaccard * 0.9, substring_score), 4)


def _target_resolver_status(score: float, lexical: float, allowed: bool) -> str:
    if not allowed:
        if lexical >= 0.9:
            return "context_match_not_authoritative"
        return "classified_non_target"
    if lexical >= 0.98 and score >= 0.9:
        return "exact"
    if score >= 0.82:
        return "probable"
    if score >= 0.68:
        return "possible"
    return "weak"


def resolve_target_query(
    data: TargetGraphData,
    query: str,
    limit: int = 8,
    target_type: str = "",
) -> dict[str, Any]:
    """Resolve a free-text target label while preserving non-target classifications."""
    query = str(query or "").strip()
    norm = _resolver_norm(query)
    classification = _target_query_type(query)
    if not norm:
        return {
            "query": query,
            "normalized_query": "",
            **classification,
            "best": None,
            "candidates": [],
        }

    target_type = (target_type or "").strip().lower()
    allowed = bool(classification.get("authoritative_mapping_allowed"))
    if not allowed:
        return {
            "query": query,
            "normalized_query": norm,
            **classification,
            "best": None,
            "candidates": [],
            "not_mapped_reason": "classified_as_non_single_target_or_non_target",
        }

    candidate_terms: dict[str, tuple[dict[str, Any], float]] = {}

    direct_target_id = _resolve_target_id(data, query)
    if direct_target_id and (not target_type or data.nodes_by_id.get(direct_target_id, {}).get("target_type", "").lower() == target_type):
        candidate_terms[direct_target_id] = ({
            "target_id": direct_target_id,
            "term": query,
            "norm": norm,
            "field": "target_id" if direct_target_id == query else "external_id",
        }, 1.0)

    token_candidates: set[str] = set()
    for token in _resolver_tokens(query):
        ids = data.token_index.get(token, set())
        if len(ids) <= 10000:
            token_candidates.update(ids)
    exact = data.exact_terms.get(norm, set())
    token_candidates.update(exact)

    candidate_ids = set(token_candidates) | set(candidate_terms.keys())
    if not candidate_ids and norm:
        # Final bounded fallback for short/noisy labels: reuse the explorer search
        # ranking but cap work aggressively so resolver calls stay interactive.
        candidate_ids.update(row.get("target_id", "") for row in _filter_targets(data, q=query)[:250])

    for term_tid in candidate_ids:
        row = data.nodes_by_id.get(term_tid, {})
        if not row:
            continue
        if target_type and row.get("target_type", "").lower() != target_type:
            continue
        for term in _target_row_resolver_terms(row):
            lexical = _target_lexical_score(norm, term.get("norm", ""))
            if lexical < 0.58:
                continue
            current = candidate_terms.get(term_tid)
            if current is None or lexical > current[1]:
                candidate_terms[term_tid] = (term, lexical)

    candidates: list[dict[str, Any]] = []
    for target_id, (term, lexical) in candidate_terms.items():
        row = data.nodes_by_id.get(target_id, {})
        if not row:
            continue
        if target_type and row.get("target_type", "").lower() != target_type:
            continue
        quality = _canonical_rank(row)
        canonical_score = {0: 1.0, 1: 0.92, 2: 0.84, 3: 0.78}.get(quality, 0.75)
        field_prior = _target_term_priority(term.get("field", ""), row)
        type_prior = 1.0 - min(_type_rank(row.get("target_type", "")), 3) * 0.025
        score = 0.84 * lexical + 0.08 * field_prior + 0.05 * canonical_score + 0.03 * type_prior
        if not allowed:
            score = min(score, 0.74)
        score = max(0.0, min(1.0, score))
        warnings = list(classification.get("warnings", []))
        if row.get("target_type", "") == "protein" and row.get("canonical_status", "") and _canonical_rank(row) > 0:
            warnings.append("protein_isoform_or_alternate_product")
        if score < 0.68:
            warnings.append("low_resolver_score")
        candidates.append({
            **_target_to_row(row),
            "resolver_score": round(score, 4),
            "lexical_score": round(lexical, 4),
            "match_status": _target_resolver_status(score, lexical, allowed),
            "query_type": classification.get("query_type", ""),
            "authoritative_mapping_allowed": allowed,
            "matched_term": term.get("term", ""),
            "matched_field": term.get("field", ""),
            "warnings": sorted(set(warnings)),
            "href": f"/target-id-qa?ids={row.get('target_id', '')}&tab=graph",
        })

    candidates.sort(key=lambda row: (
        -row["resolver_score"],
        -row["lexical_score"],
        _type_rank(row.get("target_type", "")),
        _canonical_rank(row),
        row.get("symbol", "") or row.get("name", ""),
    ))
    candidates = candidates[: max(1, min(limit, 50))]
    if len(candidates) > 1 and candidates[0]["resolver_score"] - candidates[1]["resolver_score"] <= 0.03:
        candidates[0]["warnings"] = sorted(set(candidates[0].get("warnings", []) + ["ambiguous_close_candidate"]))

    return {
        "query": query,
        "normalized_query": norm,
        **classification,
        "best": candidates[0] if candidates else None,
        "candidates": candidates,
    }


def bulk_resolve_targets(
    data: TargetGraphData,
    queries: list[str],
    limit: int = 8,
    target_type: str = "",
) -> dict[str, Any]:
    clean_queries = [str(q).strip() for q in queries if str(q).strip()]
    limit = max(1, min(limit, 50))
    return {
        "count": len(clean_queries),
        "limit": limit,
        "target_type": target_type,
        "results": [
            resolve_target_query(data, query, limit=limit, target_type=target_type)
            for query in clean_queries
        ],
    }


def _is_canonical(row: dict[str, str]) -> bool:
    return _canonical_rank(row) == 0


def _short_primary_id(primary_id: str) -> str:
    if not primary_id:
        return ""
    if ":" in primary_id:
        return primary_id.split(":", 1)[1]
    return primary_id


def _target_node_label(row: dict[str, str]) -> str:
    """Symbol-first graph label with identifier context kept secondary."""
    target_type = (row.get("target_type") or "").lower()
    symbol = row.get("symbol", "") or ""
    name = row.get("name", "") or ""
    primary_id = row.get("primary_id", "") or ""
    short_id = _short_primary_id(primary_id)
    annotations = _parse_annotations(row)
    isoform = str(annotations.get("uniprot_isoform") or "").strip()

    if target_type == "gene":
        label = symbol or name or short_id or row.get("target_id", "")
    elif target_type == "protein":
        label = symbol or name or "Protein"
        secondary = isoform or short_id
        if secondary and secondary not in label:
            label = f"{label}\n{secondary}"
    elif target_type == "transcript":
        label = symbol or name or "Transcript"
        if short_id and short_id not in label:
            label = f"{label}\n{short_id}"
    else:
        label = symbol or name or short_id or row.get("target_id", "")

    if _is_canonical(row):
        label = f"{label}\n\u2605 canonical"
    return label


def _build_product_element(
    data: TargetGraphData,
    product_id: str,
    parent_id: str,
    edge_row: dict[str, str],
) -> tuple[dict[str, Any], dict[str, Any]] | None:
    """Build a (node, edge) pair for one relation product."""
    product = data.nodes_by_id.get(product_id)
    if not product:
        return None
    product_type = product.get("target_type", "")
    canonical = _is_canonical(product)
    classes = f"target {product_type}"
    if canonical:
        classes += " canonical"
    if _parse_annotations(product).get("uniprot_isoform"):
        classes += " isoform"
    div_rows = data.divergence_by_target_id.get(product_id, [])
    if div_rows:
        classes += " has-conflict"
    node = {
        "data": {
            "id": product_id,
            "label": _target_node_label(product),
            **_target_to_row(product),
        },
        "classes": classes,
    }
    edge = _build_edge_element(parent_id, product_id, edge_row)
    _annotate_conflict_edge(data, edge)
    return node, edge


def _collect_products_by_type(
    data: TargetGraphData,
    parent_id: str,
    max_expandable: int = 25,
) -> tuple[list[dict], dict[str, list[dict]], dict[str, int]]:
    """Collect products split into initial (canonical) and expandable.

    Shows the proper hierarchy: gene → transcript → protein.
    The initial view connects the canonical protein to its translating
    transcript via ``biolink:translates_to`` (not gene→protein).

    Returns (initial_elements, expandable_by_type, total_counts).
    """
    # Step 1: Collect gene→transcript and gene→protein edges
    transcript_info: list[tuple[str, dict[str, str], dict[str, str]]] = []
    protein_info: list[tuple[str, dict[str, str], dict[str, str]]] = []
    seen: set[str] = set()

    for edge_row in data.relation_edges_by_target.get(parent_id, []):
        predicate = edge_row.get("predicate", "")
        if predicate not in ("biolink:transcribed_to", "biolink:has_gene_product"):
            continue
        source_id = edge_row.get("source_id", "")
        target_id = edge_row.get("target_id", "")
        product_id = target_id if source_id == parent_id else source_id
        if product_id == parent_id or product_id in seen:
            continue
        seen.add(product_id)
        product = data.nodes_by_id.get(product_id)
        if not product:
            continue
        product_type = product.get("target_type", "")
        if product_type == "transcript":
            transcript_info.append((product_id, product, edge_row))
        elif product_type == "protein":
            protein_info.append((product_id, product, edge_row))

    # Sort: canonical first
    transcript_info.sort(key=lambda x: _canonical_rank(x[1]))
    protein_info.sort(key=lambda x: _canonical_rank(x[1]))

    total_counts: dict[str, int] = {}
    if transcript_info:
        total_counts["transcript"] = len(transcript_info)
    if protein_info:
        total_counts["protein"] = len(protein_info)
    if not transcript_info and not protein_info:
        return [], {}, total_counts

    # Step 2: Build translates_to index (protein -> (transcript, edge))
    transcript_id_set = {t[0] for t in transcript_info}
    transcript_edge_by_id = {tid: edge for tid, _row, edge in transcript_info}
    translates_to: dict[str, tuple[str, dict[str, str]]] = {}
    for pid, _prow, _pedge in protein_info:
        for edge_row in data.relation_edges_by_target.get(pid, []):
            if edge_row.get("predicate") != "biolink:translates_to":
                continue
            src = edge_row.get("source_id", "")
            tgt = edge_row.get("target_id", "")
            tid = src if tgt == pid else tgt
            if tid in transcript_id_set:
                translates_to[pid] = (tid, edge_row)
                break

    # Step 3: Build initial elements — gene → transcript → protein
    initial: list[dict] = []
    initial_transcript_id: str | None = None
    initial_protein_id: str | None = None

    if protein_info:
        initial_protein_id = protein_info[0][0]

        if initial_protein_id in translates_to:
            # Use the transcript that translates to the canonical protein
            initial_transcript_id, tt_edge = translates_to[initial_protein_id]
            gene_product_edge = protein_info[0][2]
            for tid, _trow, tedge in transcript_info:
                if tid == initial_transcript_id:
                    result = _build_product_element(data, tid, parent_id, tedge)
                    if result:
                        initial.extend(result)
                    break
            # Connect protein to transcript via translates_to
            result = _build_product_element(data, initial_protein_id, initial_transcript_id, tt_edge)
            if result:
                initial.extend(result)
            gp_edge = _build_edge_element(parent_id, initial_protein_id, gene_product_edge)
            _annotate_conflict_edge(data, gp_edge)
            initial.append(gp_edge)
        else:
            # Fallback: no translates_to found
            if transcript_info:
                initial_transcript_id = transcript_info[0][0]
                result = _build_product_element(data, initial_transcript_id, parent_id, transcript_info[0][2])
                if result:
                    initial.extend(result)
            result = _build_product_element(data, initial_protein_id, parent_id, protein_info[0][2])
            if result:
                initial.extend(result)
    elif transcript_info:
        initial_transcript_id = transcript_info[0][0]
        result = _build_product_element(data, initial_transcript_id, parent_id, transcript_info[0][2])
        if result:
            initial.extend(result)

    # Step 4: Expandable products
    expandable_by_type: dict[str, list[dict]] = {}

    remaining_transcripts = [
        (tid, trow, tedge) for tid, trow, tedge in transcript_info
        if tid != initial_transcript_id
    ][:max_expandable]
    if remaining_transcripts:
        exp_t: list[dict] = []
        for tid, _trow, tedge in remaining_transcripts:
            result = _build_product_element(data, tid, parent_id, tedge)
            if result:
                exp_t.extend(result)
        if exp_t:
            expandable_by_type["transcript"] = exp_t

    remaining_proteins = [
        (pid, prow, pedge) for pid, prow, pedge in protein_info
        if pid != initial_protein_id
    ][:max_expandable]
    if remaining_proteins:
        exp_p: list[dict] = []
        for pid, _prow, pedge in remaining_proteins:
            used_translation_edge = False
            if pid in translates_to:
                tt_tid, tt_edge = translates_to[pid]
                if tt_tid != initial_transcript_id:
                    transcript_gene_edge = transcript_edge_by_id.get(tt_tid)
                    if transcript_gene_edge:
                        transcript_result = _build_product_element(data, tt_tid, parent_id, transcript_gene_edge)
                        if transcript_result:
                            exp_p.extend(transcript_result)
                result = _build_product_element(data, pid, tt_tid, tt_edge)
                used_translation_edge = True
            else:
                result = _build_product_element(data, pid, parent_id, pedge)
            if result:
                exp_p.extend(result)
                if used_translation_edge:
                    extra_edge = _build_edge_element(parent_id, pid, pedge)
                    _annotate_conflict_edge(data, extra_edge)
                    exp_p.append(extra_edge)
        if exp_p:
            expandable_by_type["protein"] = exp_p

    return initial, expandable_by_type, total_counts


def build_target_graph_payload(
    data: TargetGraphData,
    ids: str,
) -> dict[str, Any]:
    selected: list[str] = []
    for token in re.split(r"[\s,|]+", ids or ""):
        target_id = _resolve_target_id(data, token)
        if target_id and target_id not in selected:
            selected.append(target_id)
    if not selected:
        return {"elements": [], "selected": []}

    elements: list[dict[str, Any]] = []
    seen_nodes: set[str] = set()
    seen_edges: set[str] = set()

    def add_node(node_id: str, label: str, classes: str, payload: dict[str, Any] | None = None) -> None:
        if node_id in seen_nodes:
            return
        seen_nodes.add(node_id)
        elements.append({"data": {"id": node_id, "label": label, **(payload or {})}, "classes": classes})

    def add_edge(source: str, target: str, label: str, classes: str, payload: dict[str, Any] | None = None) -> None:
        payload = payload or {}
        evidence_key = payload.get("evidence_identifier", "")
        edge_id = f"{source}->{target}:{label}:{evidence_key}"
        if edge_id in seen_edges:
            return
        seen_edges.add(edge_id)
        elements.append({
            "data": {"id": edge_id, "source": source, "target": target, "label": label, **payload},
            "classes": classes,
        })

    def add_element(element: dict[str, Any]) -> None:
        payload = element.get("data", {})
        element_id = payload.get("id", "")
        if payload.get("source") and payload.get("target"):
            if element_id in seen_edges:
                return
            seen_edges.add(element_id)
            elements.append(element)
            return
        if element_id in seen_nodes:
            return
        seen_nodes.add(element_id)
        elements.append(element)

    def add_target_node(target_id: str) -> None:
        row = data.nodes_by_id[target_id]
        target_type = row.get("target_type", "")
        canonical = _is_canonical(row)
        classes = f"target {target_type}"
        if canonical:
            classes += " canonical"
        if _parse_annotations(row).get("uniprot_isoform"):
            classes += " isoform"
        if data.divergence_by_target_id.get(target_id):
            classes += " has-conflict"
        add_node(target_id, _target_node_label(row), classes, _target_to_row(row))

    # Build selected nodes and keep gene anchors for product expansion.
    context_genes: list[str] = []

    def add_context_gene(gene_id: str) -> None:
        if gene_id and gene_id not in context_genes:
            context_genes.append(gene_id)

    for target_id in selected:
        add_target_node(target_id)
        row = data.nodes_by_id[target_id]
        if row.get("target_type") == "gene":
            add_context_gene(target_id)
            continue
        for edge_row in data.relation_edges_by_target.get(target_id, []):
            predicate = edge_row.get("predicate", "")
            if predicate not in {
                "biolink:transcribed_to",
                "biolink:translates_to",
                "biolink:has_gene_product",
            }:
                continue
            source_id = edge_row.get("source_id", "")
            edge_target_id = edge_row.get("target_id", "")
            if not source_id or not edge_target_id:
                continue
            for node_id in (source_id, edge_target_id):
                if node_id in data.nodes_by_id:
                    add_target_node(node_id)
                    if data.nodes_by_id[node_id].get("target_type") == "gene":
                        add_context_gene(node_id)
            edge_element = _build_edge_element(source_id, edge_target_id, edge_row)
            _annotate_conflict_edge(data, edge_element)
            add_element(edge_element)

    # Add initial canonical products and collect expandable products from gene anchors.
    expandable_transcripts: dict[str, list[dict]] = {}
    expandable_proteins: dict[str, list[dict]] = {}
    product_counts: dict[str, dict[str, int]] = {}
    for target_id in context_genes:
        initial, expandable_by_type, counts = _collect_products_by_type(data, target_id)
        # Add canonical transcript + protein to initial elements
        for element in initial:
            add_element(element)
        if expandable_by_type.get("transcript"):
            expandable_transcripts[target_id] = expandable_by_type["transcript"]
        if expandable_by_type.get("protein"):
            expandable_proteins[target_id] = expandable_by_type["protein"]
        if counts:
            product_counts[target_id] = counts
        # Annotate gene node with total and expandable counts
        for el in elements:
            if el.get("data", {}).get("id") == target_id:
                for ptype, count in counts.items():
                    el["data"][f"total_{ptype}_count"] = str(count)
                    el["data"][f"expandable_{ptype}_count"] = str(max(0, count - 1))
                break

    result: dict[str, Any] = {"elements": elements, "selected": selected}
    if expandable_transcripts:
        result["expandableTranscripts"] = expandable_transcripts
    if expandable_proteins:
        result["expandableProteins"] = expandable_proteins
    if product_counts:
        result["productCounts"] = product_counts
    return result


def export_targets(
    data: TargetGraphData,
    q: str = "",
    target_type: str = "",
    namespace: str = "",
    fmt: str = "tsv",
    columns: list[str] | None = None,
) -> str:
    rows = _filter_targets(data, q=q, target_type=target_type, namespace=namespace)
    default_columns = [
        "target_id", "target_type", "primary_id", "symbol", "name", "description", "category",
        "ids", "id_namespaces", "source_namespaces", "mapping_ratio", "quality_note",
        "canonical_status", "canonical_ifx_id",
    ]
    columns = [col for col in (columns or default_columns) if col]
    if not columns:
        columns = default_columns
    delimiter = "," if fmt == "csv" else "\t"
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
    writer.writeheader()
    for row in rows:
        base = _target_to_row(row)
        annotations = _parse_annotations(row)
        writer.writerow({col: base.get(col, annotations.get(col, "")) for col in columns})
    return buf.getvalue()


# ---------------------------------------------------------------------------
# Pipeline-format entity ID exports (for Keith's target_graph ingestion)
# ---------------------------------------------------------------------------

def _extract_id(ids_field: str, prefix: str) -> str:
    """Extract ID value(s) for a given namespace prefix from the pipe-delimited ids field."""
    vals = []
    for raw in ids_field.split("|"):
        raw = raw.strip()
        if raw.startswith(prefix + ":"):
            vals.append(raw[len(prefix) + 1:])
    return "|".join(vals)


# Column definitions for each entity type — maps output column name to
# either a direct field from _target_to_row, a namespace prefix to extract
# from the ``ids`` field, or a key from ``source_annotations`` JSON.

_GENE_ID_COLUMNS = [
    ("ncats_gene_id", "field", "target_id"),
    ("consolidated_NCBI_id", "id", "NCBIGene"),
    ("consolidated_hgnc_id", "id", "HGNC"),
    ("consolidated_symbol", "field", "symbol"),
    ("consolidated_gene_id", "id", "ENSEMBL"),
    ("consolidated_mim_id", "id", "OMIM"),
    ("consolidated_description", "field", "description"),
    ("consolidated_gene_type", "field", "category"),
    ("ids", "field", "ids"),
    ("id_namespaces", "field", "id_namespaces"),
    ("source_namespaces", "field", "source_namespaces"),
    ("Total_Mapping_Ratio", "field", "mapping_ratio"),
    ("quality_note", "field", "quality_note"),
    ("Symbol_Provenance", "annotation", "Symbol_Provenance"),
    ("Description_Provenance", "annotation", "Description_Provenance"),
    ("NCBI_ID_Provenance", "annotation", "NCBI_ID_Provenance"),
    ("HGNC_ID_Provenance", "annotation", "HGNC_ID_Provenance"),
    ("Ensembl_ID_Provenance", "annotation", "Ensembl_ID_Provenance"),
    ("Location_Provenance", "annotation", "Location_Provenance"),
    ("Mapping_Support_Tier", "annotation", "Mapping_Support_Tier"),
    ("Total_Mapping_Score", "annotation", "Total_Mapping_Score"),
    ("createdAt", "field", "createdAt"),
    ("updatedAt", "field", "updatedAt"),
]

_PROTEIN_ID_COLUMNS = [
    ("ncats_protein_id", "field", "target_id"),
    ("uniprot_id", "id", "UniProtKB"),
    ("consolidated_ensembl_protein_id", "id", "ENSEMBL"),
    ("consolidated_ensembl_transcript_id", "annotation", "consolidated_ensembl_transcript_id"),
    ("consolidated_refseq_protein", "id", "RefSeq"),
    ("consolidated_symbol", "field", "symbol"),
    ("combined_protein_name", "field", "name"),
    ("uniprot_entryType", "field", "description"),
    ("is_canonical", "annotation", "is_canonical"),
    ("canonical_isoform_status", "field", "canonical_status"),
    ("canonical_isoform_status_source", "annotation", "canonical_isoform_status_source"),
    ("uniprot_computationally_mapped_isoform", "annotation", "uniprot_computationally_mapped_isoform"),
    ("protein_group_anchor", "annotation", "protein_group_anchor"),
    ("canonical_ifx_id", "field", "canonical_ifx_id"),
    ("protein_grouping_reason", "annotation", "protein_grouping_reason"),
    ("ensembl_canonical", "annotation", "ensembl_canonical"),
    ("ensembl_canonical_token", "annotation", "ensembl_canonical_token"),
    ("uniprot_isoform", "annotation", "uniprot_isoform"),
    ("uniref100_cluster_id", "annotation", "uniref100_cluster_id"),
    ("parent_uniref100_cluster_id", "annotation", "parent_uniref100_cluster_id"),
    ("same_uniref100_as_parent", "annotation", "same_uniref100_as_parent"),
    ("uniprot_annotationScore", "annotation", "uniprot_annotationScore"),
    ("ensembl_transcript_tsl", "annotation", "ensembl_transcript_tsl"),
    ("refseq_method", "annotation", "refseq_method"),
    ("protein_name_score", "annotation", "protein_name_score"),
    ("protein_name_method", "annotation", "protein_name_method"),
    ("uniprot_secondaryAccessions", "annotation", "uniprot_secondaryAccessions"),
    ("ids", "field", "ids"),
    ("id_namespaces", "field", "id_namespaces"),
    ("source_namespaces", "field", "source_namespaces"),
    ("Total_Mapping_Ratio", "field", "mapping_ratio"),
    ("quality_note", "field", "quality_note"),
    ("Mapping_Support_Tier", "annotation", "Mapping_Support_Tier"),
    ("Total_Mapping_Score", "annotation", "Total_Mapping_Score"),
    ("UniProt_ID_Provenance", "annotation", "UniProt_ID_Provenance"),
    ("Ensembl_ID_Provenance", "annotation", "Ensembl_ID_Provenance"),
    ("RefSeq_ID_Provenance", "annotation", "RefSeq_ID_Provenance"),
    ("createdAt", "field", "createdAt"),
    ("updatedAt", "field", "updatedAt"),
]

_TRANSCRIPT_ID_COLUMNS = [
    ("ncats_transcript_id", "field", "target_id"),
    ("ensembl_transcript_id", "id", "ENSEMBL"),
    ("consolidated_symbol", "field", "symbol"),
    ("consolidated_refseq_rna", "id", "RefSeq"),
    ("ids", "field", "ids"),
    ("id_namespaces", "field", "id_namespaces"),
    ("source_namespaces", "field", "source_namespaces"),
    ("Total_Mapping_Ratio", "field", "mapping_ratio"),
    ("quality_note", "field", "quality_note"),
    ("canonical_status", "field", "canonical_status"),
    ("Mapping_Support_Tier", "annotation", "Mapping_Support_Tier"),
    ("Total_Mapping_Score", "annotation", "Total_Mapping_Score"),
    ("Ensembl_Transcript_ID_Provenance", "annotation", "Ensembl_Transcript_ID_Provenance"),
    ("RefSeq_Provenance", "annotation", "RefSeq_Provenance"),
    ("createdAt", "field", "createdAt"),
    ("updatedAt", "field", "updatedAt"),
]

_ENTITY_COLUMN_MAP = {
    "gene": _GENE_ID_COLUMNS,
    "protein": _PROTEIN_ID_COLUMNS,
    "transcript": _TRANSCRIPT_ID_COLUMNS,
}


def export_entity_ids(
    data: TargetGraphData,
    entity_type: str,
    fmt: str = "tsv",
) -> str:
    """Export entity IDs in the pipeline-native column format for ingestion.

    Reconstructs the original column names (``ncats_gene_id``,
    ``uniprot_id``, etc.) from the app-graph's normalised fields so Keith's
    ``target_graph`` build pipeline can consume the output directly.
    """
    entity_type = entity_type.lower().strip()
    col_spec = _ENTITY_COLUMN_MAP.get(entity_type)
    if not col_spec:
        raise ValueError(f"Unknown entity type: {entity_type!r}")

    columns = [c[0] for c in col_spec]
    delimiter = "," if fmt == "csv" else "\t"
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
    writer.writeheader()

    for node in data.nodes:
        if (node.get("target_type") or "").lower() != entity_type:
            continue

        base = _target_to_row(node)
        annotations = _parse_annotations(node)
        ids_field = base.get("ids", "")

        out: dict[str, str] = {}
        for col_name, source_type, source_key in col_spec:
            if source_type == "field":
                out[col_name] = base.get(source_key, "")
            elif source_type == "id":
                out[col_name] = _extract_id(ids_field, source_key)
            elif source_type == "annotation":
                out[col_name] = str(annotations.get(source_key, ""))

        writer.writerow(out)

    return buf.getvalue()


_DIVERGENCE_EXPORT_COLUMNS = [
    "registry_id", "entity_type", "standard_name", "row_id", "source",
    "namespace", "divergence_type", "source_value", "consensus_value",
    "other_sources", "n_sources_agree", "n_sources_total",
    "auto_decision", "auto_confidence", "auto_rationale",
    "scenario_id", "status", "finding_summary",
]


def export_divergences(
    data: TargetGraphData,
    entity_type: str = "",
    triage_category: str = "",
    review_group: str = "",
    status: str = "",
    fmt: str = "tsv",
) -> str:
    """Export divergence/QC rows as TSV or CSV for offline analysis."""
    entity_type = (entity_type or "").strip().lower()
    triage_filter = (triage_category or "").strip().lower()
    review_group_filter = (review_group or "").strip().lower()
    status_filter = (status or "").strip().lower()

    delimiter = "," if fmt == "csv" else "\t"
    columns = _DIVERGENCE_EXPORT_COLUMNS
    buf = io.StringIO()
    writer = csv.DictWriter(buf, fieldnames=columns, delimiter=delimiter, extrasaction="ignore")
    writer.writeheader()

    for etype, rows in [
        ("gene", data.divergence_gene),
        ("protein", data.divergence_protein),
        ("transcript", data.divergence_transcript),
    ]:
        if entity_type and etype != entity_type:
            continue
        for row in rows:
            row_status = (row.get("status", "") or "open").lower()
            if status_filter and row_status != status_filter:
                continue
            if triage_filter and classify_triage_category(row) != triage_filter:
                continue
            if review_group_filter and classify_review_group(etype, row) != review_group_filter:
                continue
            out = {col: row.get(col, "") for col in columns}
            out["entity_type"] = etype
            writer.writerow(out)

    return buf.getvalue()


def build_target_review_queue(
    data: TargetGraphData,
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
) -> dict[str, Any]:
    """Build a paginated, filterable review queue from divergence registries."""
    entity_type = (entity_type or "").strip().lower()
    divergence_type = (divergence_type or "").strip().lower()
    source_filter = (source or "").strip().lower()
    status_filter = (status or "").strip().lower()
    auto_decision_filter = (auto_decision or "").strip().lower()
    triage_filter = (triage_category or "").strip().lower()
    review_group_filter = (review_group or "").strip().lower()
    q_lower = (q or "").strip().lower()
    summary_divergence = (
        data.public_qc_summary.get("divergence", {})
        if data.review_mode == "public" and isinstance(data.public_qc_summary.get("divergence"), dict)
        else {}
    )

    # ── Single pass: total counts + build combined list ──
    triage_category_totals: Counter[str] = Counter(summary_divergence.get("triage_category_counts") or {})
    review_group_totals: Counter[str] = Counter(summary_divergence.get("review_group_counts") or {})
    total_all = 0
    total_resolved = 0
    if summary_divergence:
        total_all = int(summary_divergence.get("total_divergences") or 0)
        total_resolved = int((summary_divergence.get("status_counts") or {}).get("resolved") or 0)
    combined: list[tuple[str, dict[str, str], str, str]] = []
    for etype, rows in [
        ("gene", data.divergence_gene),
        ("protein", data.divergence_protein),
        ("transcript", data.divergence_transcript),
    ]:
        for row in rows:
            if not summary_divergence:
                total_all += 1
            row_triage = classify_triage_category(row)
            row_group = classify_review_group(etype, row)
            if not summary_divergence:
                triage_category_totals[row_triage] += 1
                review_group_totals[row_group] += 1
            if (row.get("status", "") or "open").lower() == "resolved":
                if not summary_divergence:
                    total_resolved += 1
            combined.append((etype, row, row_triage, row_group))

    # ── Shared-filter pass: apply all filters EXCEPT entity_type, triage, and review group ──
    def _passes_shared(etype: str, row: dict[str, str], row_triage: str, row_group: str) -> bool:
        row_div_type = (row.get("divergence_type", "") or "").lower()
        row_status = (row.get("status", "") or "open").lower()
        row_source = (row.get("source", "") or "").lower()
        row_auto_dec = (row.get("auto_decision", "") or "").lower()
        if divergence_type and row_div_type != divergence_type:
            return False
        if status_filter and row_status != status_filter:
            return False
        if source_filter and row_source != source_filter:
            return False
        if auto_decision_filter and row_auto_dec != auto_decision_filter:
            return False
        if q_lower:
            haystack = " ".join([
                row.get("standard_name", ""),
                row.get("row_id", ""),
                row.get("source_value", ""),
                row.get("consensus_value", ""),
                row.get("auto_rationale", ""),
                row.get("registry_id", ""),
            ]).lower()
            if q_lower not in haystack:
                return False
        return True

    # Two-pass counting for cross-filter stability:
    # - entity_type_counts: from base_pool filtered by triage (but NOT entity_type)
    # - triage_counts: from base_pool filtered by entity_type (but NOT triage)
    entity_type_counts: Counter[str] = Counter()
    triage_counts: Counter[str] = Counter()
    review_group_counts: Counter[str] = Counter()
    divergence_type_counts: Counter[str] = Counter()
    source_counts: Counter[str] = Counter()
    status_counts: Counter[str] = Counter()
    auto_decision_counts: Counter[str] = Counter()
    filtered: list[tuple[str, dict[str, str], str, str]] = []

    for etype, row, row_triage, row_group in combined:
        if not _passes_shared(etype, row, row_triage, row_group):
            continue
        row_div_type = (row.get("divergence_type", "") or "").lower()
        row_status = (row.get("status", "") or "open").lower()
        row_source = (row.get("source", "") or "").lower()
        row_auto_dec = (row.get("auto_decision", "") or "").lower()

        # Entity counts: apply triage + review group filters but NOT entity_type
        if (not triage_filter or row_triage == triage_filter) and (
            not review_group_filter or row_group == review_group_filter
        ):
            entity_type_counts[etype] += 1

        # Triage counts: apply entity_type + review group filters but NOT triage
        if (not entity_type or etype == entity_type) and (
            not review_group_filter or row_group == review_group_filter
        ):
            triage_counts[row_triage] += 1

        # Review group counts: apply entity_type + triage filters but NOT review group
        if (not entity_type or etype == entity_type) and (
            not triage_filter or row_triage == triage_filter
        ):
            review_group_counts[row_group] += 1

        # All filters for final list and remaining facets
        if entity_type and etype != entity_type:
            continue
        if triage_filter and row_triage != triage_filter:
            continue
        if review_group_filter and row_group != review_group_filter:
            continue

        divergence_type_counts[row_div_type] += 1
        if row_source:
            source_counts[row_source] += 1
        status_counts[row_status] += 1
        if row_auto_dec:
            auto_decision_counts[row_auto_dec] += 1
        filtered.append((etype, row, row_triage, row_group))

    # Sort: source_conflict first, then cosmetic, low_mapping, coverage_gap.
    # Within each category, lowest confidence first.
    def _sort_key(item: tuple[str, dict[str, str], str, str]) -> tuple:
        etype, row, triage, group = item
        try:
            conf = float(row.get("auto_confidence", "1.0"))
        except (ValueError, TypeError):
            conf = 1.0
        return (
            _TRIAGE_ORDER.get(triage, 9),
            _REVIEW_GROUP_ORDER.get(group, 99),
            conf,
            row.get("standard_name", "").lower(),
        )

    filtered.sort(key=_sort_key)

    # Paginate
    total = len(filtered)
    per_page = max(1, min(per_page, 200))
    total_pages = max(1, (total + per_page - 1) // per_page)
    page = max(1, min(page, total_pages))
    start = (page - 1) * per_page
    rows_out: list[dict[str, Any]] = []
    for etype, row, triage, group in filtered[start:start + per_page]:
        group_info = _REVIEW_GROUP_BY_KEY.get(group, {})
        out = dict(row)
        out.pop("_triage_category", None)
        out.pop("_review_group", None)
        out["entity_type"] = etype
        out["triage_category"] = triage
        out["review_group"] = group
        out["review_group_label"] = group_info.get("label", group)
        out["review_group_description"] = group_info.get("description", "")
        out["review_group_recommended_decision"] = group_info.get("recommended_decision", "")
        out["suggested_decision"] = (
            group_info.get("recommended_decision")
            or _TRIAGE_SUGGESTED_DECISION.get(triage, "")
        )
        # Enrich conflict rows with cross-reference lookups (only for page)
        if triage == "source_conflict":
            xref_ctx = _enrich_conflict_xrefs(data, row)
            out["xref_context"] = xref_ctx
            # Compute automated suggestion from xref enrichment
            suggestion = _compute_xref_suggestion(data, row, xref_ctx)
            if suggestion:
                out["suggested_decision"] = suggestion.get("suggested_decision", "")
                out["suggestion_rationale"] = suggestion.get("suggestion_rationale", "")
                out["suggestion_confidence"] = suggestion.get("suggestion_confidence", "")
        rows_out.append(out)

    # Compute review group pattern summary when a review group filter is active
    review_group_summary: dict[str, Any] = {}
    if review_group_filter:
        review_group_summary = _compute_review_group_summary(
            data, review_group_filter, entity_type=entity_type, status=status_filter,
        )

    return {
        "rows": rows_out,
        "total": total,
        "page": page,
        "per_page": per_page,
        "total_pages": total_pages,
        "entity_type_counts": dict(entity_type_counts.most_common()),
        "divergence_type_counts": dict(divergence_type_counts.most_common()),
        "source_counts": dict(source_counts.most_common()),
        "status_counts": dict(status_counts.most_common()),
        "auto_decision_counts": dict(auto_decision_counts.most_common()),
        "triage_category_counts": dict(
            sorted(triage_counts.items(), key=lambda kv: _TRIAGE_ORDER.get(kv[0], 9))
        ),
        "review_group_counts": dict(
            sorted(review_group_counts.items(), key=lambda kv: _REVIEW_GROUP_ORDER.get(kv[0], 99))
        ),
        "triage_categories": TRIAGE_CATEGORIES,
        "review_groups": REVIEW_GROUPS,
        "triage_category_totals": dict(triage_category_totals.most_common()),
        "review_group_totals": dict(
            sorted(review_group_totals.items(), key=lambda kv: _REVIEW_GROUP_ORDER.get(kv[0], 99))
        ),
        "total_all": total_all,
        "total_resolved": total_resolved,
        "decision_options": TARGET_REVIEW_DECISION_OPTIONS,
        "review_mode": data.review_mode,
        "can_write_review": data.review_can_write,
        "public_read_only": not data.review_can_write,
        "review_group_summary": review_group_summary,
    }


def build_batch_review_payload(
    data: TargetGraphData,
    triage_category: str,
    review_decision: str,
    entity_type: str = "",
    source: str = "",
    status: str = "open",
    review_group: str = "",
    reviewed_by: str = "",
) -> list[dict[str, str]]:
    """Build review rows for all items matching a triage category.

    Used for batch-applying a decision to an entire category at once.
    """
    entity_type = (entity_type or "").strip().lower()
    source_filter = (source or "").strip().lower()
    status_filter = (status or "").strip().lower()
    review_group_filter = (review_group or "").strip().lower()
    reviewed_by = (reviewed_by or os.getenv("USER") or "app_batch_review").strip()
    reviewed_at = datetime.now(timezone.utc).isoformat()

    allowed = {opt["value"] for opt in TARGET_REVIEW_DECISION_OPTIONS}
    if review_decision not in allowed:
        raise HTTPException(status_code=400, detail=f"Unsupported review_decision: {review_decision}")

    rows: list[dict[str, str]] = []
    for etype, divergence_rows in [
        ("gene", data.divergence_gene),
        ("protein", data.divergence_protein),
        ("transcript", data.divergence_transcript),
    ]:
        if entity_type and etype != entity_type:
            continue
        for row in divergence_rows:
            row_status = (row.get("status", "") or "open").lower()
            row_source = (row.get("source", "") or "").lower()
            if status_filter and row_status != status_filter:
                continue
            if source_filter and row_source != source_filter:
                continue
            if triage_category and classify_triage_category(row) != triage_category:
                continue
            if review_group_filter and classify_review_group(etype, row) != review_group_filter:
                continue
            rows.append({
                "Registry ID": row.get("registry_id", ""),
                "Entity Type": etype,
                "Review decision": review_decision,
                "Replacement value": "",
                "Reviewer notes": f"Batch: {review_group_filter or triage_category}",
                "Standard Name": row.get("standard_name", ""),
                "Row ID": row.get("row_id", ""),
                "Source": row.get("source", ""),
                "Namespace": row.get("namespace", ""),
                "Divergence Type": row.get("divergence_type", ""),
                "Source Value": row.get("source_value", ""),
                "Consensus Value": row.get("consensus_value", ""),
                "Scenario ID": row.get("scenario_id", ""),
                "Auto Decision": row.get("auto_decision", ""),
                "Auto Confidence": row.get("auto_confidence", ""),
                "Auto Rationale": row.get("auto_rationale", ""),
                "Reviewed by": reviewed_by,
                "Reviewed at": reviewed_at,
                "App review ID": str(uuid.uuid4()),
            })
    return rows


def _invalidate_divergence_caches(data: TargetGraphData) -> None:
    """Clear cached divergence stats so they are recomputed after status changes.

    Only invalidates ``_divergence_stats`` and the nested ``divergence``
    key inside ``_stats``.  The expensive node-level target stats (type
    counts, namespace coverage, etc.) never change from review decisions
    and are preserved.
    """
    data._divergence_stats = None
    if data._stats is not None:
        data._stats.pop("divergence", None)


def mark_rows_resolved_by_registry_ids(
    data: TargetGraphData, registry_ids: set[str],
) -> int:
    """Mark matching divergence rows as resolved in memory."""
    count = 0
    for rows in (data.divergence_gene, data.divergence_protein, data.divergence_transcript):
        for row in rows:
            if row.get("registry_id", "") in registry_ids and row.get("status", "open") != "resolved":
                row["status"] = "resolved"
                count += 1
    if count:
        _invalidate_divergence_caches(data)
    return count


def mark_rows_resolved_by_triage(
    data: TargetGraphData,
    triage_category: str,
    entity_type: str = "",
    source: str = "",
    status: str = "open",
    review_group: str = "",
) -> int:
    """Mark all divergence rows matching a triage category as resolved in memory."""
    entity_type = (entity_type or "").strip().lower()
    source_filter = (source or "").strip().lower()
    status_filter = (status or "").strip().lower()
    review_group_filter = (review_group or "").strip().lower()
    count = 0
    for etype, rows in [
        ("gene", data.divergence_gene),
        ("protein", data.divergence_protein),
        ("transcript", data.divergence_transcript),
    ]:
        if entity_type and etype != entity_type:
            continue
        for row in rows:
            row_status = (row.get("status", "") or "open").lower()
            row_source = (row.get("source", "") or "").lower()
            if status_filter and row_status != status_filter:
                continue
            if source_filter and row_source != source_filter:
                continue
            if triage_category and classify_triage_category(row) != triage_category:
                continue
            if review_group_filter and classify_review_group(etype, row) != review_group_filter:
                continue
            row["status"] = "resolved"
            count += 1
    if count:
        _invalidate_divergence_caches(data)
    return count


def validate_and_build_target_review_rows(payload: dict) -> list[dict[str, str]]:
    """Validate review payload and return rows ready for TSV append."""
    raw_rows = payload.get("decisions") or payload.get("rows") or []
    if isinstance(raw_rows, dict):
        raw_rows = [raw_rows]
    if not isinstance(raw_rows, list):
        raise HTTPException(status_code=400, detail="decisions must be a list.")

    reviewed_by = str(payload.get("reviewed_by") or os.getenv("USER") or "app_review").strip()
    reviewed_at = datetime.now(timezone.utc).isoformat()
    allowed = {opt["value"] for opt in TARGET_REVIEW_DECISION_OPTIONS}
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
        if decision not in allowed:
            raise HTTPException(status_code=400, detail=f"Unsupported review_decision: {decision}")

        rows.append({
            "Registry ID": registry_id,
            "Entity Type": str(raw.get("entity_type") or raw.get("Entity Type") or "").strip(),
            "Review decision": decision,
            "Replacement value": str(raw.get("replacement") or raw.get("Replacement value") or "").strip(),
            "Reviewer notes": str(raw.get("notes") or raw.get("Reviewer notes") or "").strip(),
            "Standard Name": str(raw.get("standard_name") or raw.get("Standard Name") or "").strip(),
            "Row ID": str(raw.get("row_id") or raw.get("Row ID") or "").strip(),
            "Source": str(raw.get("source") or raw.get("Source") or "").strip(),
            "Namespace": str(raw.get("namespace") or raw.get("Namespace") or "").strip(),
            "Divergence Type": str(raw.get("divergence_type") or raw.get("Divergence Type") or "").strip(),
            "Source Value": str(raw.get("source_value") or raw.get("Source Value") or "").strip(),
            "Consensus Value": str(raw.get("consensus_value") or raw.get("Consensus Value") or "").strip(),
            "Scenario ID": str(raw.get("scenario_id") or raw.get("Scenario ID") or "").strip(),
            "Auto Decision": str(raw.get("auto_decision") or raw.get("Auto Decision") or "").strip(),
            "Auto Confidence": str(raw.get("auto_confidence") or raw.get("Auto Confidence") or "").strip(),
            "Auto Rationale": str(raw.get("auto_rationale") or raw.get("Auto Rationale") or "").strip(),
            "Reviewed by": reviewed_by,
            "Reviewed at": reviewed_at,
            "App review ID": str(raw.get("app_review_id") or uuid.uuid4()),
        })
    if not rows:
        raise HTTPException(status_code=400, detail="No review decisions provided.")
    return rows


def append_target_review_rows(rows: list[dict[str, str]], review_file_path: str | Path) -> str:
    """Thread-safe append of review rows to TSV file."""
    path = Path(review_file_path)
    path.parent.mkdir(parents=True, exist_ok=True)
    with _target_review_lock:
        exists = path.exists() and path.stat().st_size > 0
        with open(path, "a", newline="", encoding="utf-8") as fh:
            writer = csv.DictWriter(
                fh,
                fieldnames=TARGET_REVIEW_INTAKE_COLUMNS,
                delimiter="\t",
                extrasaction="ignore",
            )
            if not exists:
                writer.writeheader()
            writer.writerows(rows)
    return str(path)


# ---------------------------------------------------------------------------
# Version Comparison
# ---------------------------------------------------------------------------

_TARGET_VERSION_CACHE_MAX = 4
_target_version_cache: OrderedDict[str, TargetGraphData] = OrderedDict()


def load_target_version_data(data_dir: str | Path) -> TargetGraphData:
    """Load any versioned target app_graph directory, cached by absolute path."""
    data_dir = Path(data_dir)
    cache_key = str(data_dir.resolve())
    if cache_key in _target_version_cache:
        _target_version_cache.move_to_end(cache_key)
        return _target_version_cache[cache_key]
    data = TargetGraphData()
    _load_target_files(data, data_dir)
    _target_version_cache[cache_key] = data
    while len(_target_version_cache) > _TARGET_VERSION_CACHE_MAX:
        _target_version_cache.popitem(last=False)
    print(f"  Loaded target version {data_dir.name}: {len(data.nodes):,} nodes, {len(data.relation_edges):,} edges")
    return data


def _load_target_files(data: TargetGraphData, graph_dir: Path) -> None:
    """Load just the TSV files from an app_graph directory (no QC/divergence)."""
    nodes_file = graph_dir / "target_nodes.tsv"
    edges_file = graph_dir / "target_edges.tsv"
    manifest_file = graph_dir / "manifest.json"

    if manifest_file.exists():
        try:
            with open(manifest_file, encoding="utf-8") as fh:
                data.manifest = json.load(fh)
        except (OSError, json.JSONDecodeError):
            pass

    if nodes_file.exists():
        for row in _read_tsv(nodes_file):
            data.nodes.append(row)
            tid = row.get("target_id", "")
            if tid:
                data.nodes_by_id[tid] = row
                ttype = row.get("target_type", "")
                if ttype:
                    data.type_index[ttype].append(tid)
                for raw_id in row.get("ids", "").split("|"):
                    raw_id = raw_id.strip()
                    if raw_id:
                        data.ids_to_targets[raw_id].append(tid)

    if edges_file.exists():
        for row in _read_tsv(edges_file):
            _index_edge(data, row)


def _target_type_distribution(data: TargetGraphData) -> Counter:
    return Counter(row.get("target_type", "unknown") for row in data.nodes)


def _target_namespace_distribution(data: TargetGraphData) -> Counter:
    counts: Counter = Counter()
    for row in data.nodes:
        for ns in row.get("id_namespaces", "").split("|"):
            ns = ns.strip()
            if ns:
                counts[ns] += 1
    return counts


def _target_edge_key(edge: dict) -> tuple:
    return (
        edge.get("source_id", ""),
        edge.get("predicate", ""),
        edge.get("target_id", ""),
    )


def _distribution_delta(
    baseline: Counter,
    current: Counter,
    key_name: str,
) -> list[dict[str, Any]]:
    rows = []
    for key in sorted(set(baseline) | set(current)):
        old_count = int(baseline.get(key, 0))
        new_count = int(current.get(key, 0))
        rows.append({
            key_name: key,
            "old_count": old_count,
            "new_count": new_count,
            "delta": new_count - old_count,
        })
    return sorted(rows, key=lambda r: (-abs(r["delta"]), r[key_name]))


def _target_source_version_map(data: TargetGraphData) -> dict[str, dict[str, str]]:
    rows = data.manifest.get("source_versions") or []
    out: dict[str, dict[str, str]] = {}
    if not isinstance(rows, list):
        return out
    for raw in rows:
        if not isinstance(raw, dict):
            continue
        source = str(raw.get("source") or raw.get("source_name") or "").strip()
        if not source:
            continue
        out[source] = {
            "source": source,
            "version": str(raw.get("version") or raw.get("source_version") or "").strip(),
            "download_start": str(raw.get("download_start") or "").strip(),
            "download_end": str(raw.get("download_end") or "").strip(),
            "metadata_file": str(raw.get("metadata_file") or "").strip(),
        }
    return out


def _source_version_changes(
    baseline: TargetGraphData,
    current: TargetGraphData,
) -> list[dict[str, Any]]:
    baseline_versions = _target_source_version_map(baseline)
    current_versions = _target_source_version_map(current)
    ordered_sources = []
    for source in list(current_versions) + list(baseline_versions):
        if source not in ordered_sources:
            ordered_sources.append(source)

    rows: list[dict[str, Any]] = []
    for source in ordered_sources:
        old = baseline_versions.get(source, {})
        new = current_versions.get(source, {})
        old_version = old.get("version", "")
        new_version = new.get("version", "")
        if old and not new:
            status = "removed"
        elif new and not old:
            status = "added"
        elif old_version != new_version:
            status = "version_changed"
        else:
            status = "unchanged"
        rows.append({
            "source": source,
            "old_version": old_version,
            "new_version": new_version,
            "status": status,
            "old_download_date": old.get("download_start", "")[:10],
            "new_download_date": new.get("download_start", "")[:10],
        })
    return rows


def _format_signed_delta(value: int, label: str) -> str:
    sign = "+" if value > 0 else ""
    return f"{label} {sign}{value:,}"


_CHANGE_OWNER_EXCLUDED_NAMESPACES = {"NodeNorm", "UMLS", "UniRef100"}


def _target_identifier_values(row: dict[str, str], *, identity_only: bool = False) -> list[str]:
    values: list[str] = []
    seen: set[str] = set()
    for field in ("primary_id", "ids"):
        for value in _split_pipe(row.get(field, "")):
            namespace = value.split(":", 1)[0] if ":" in value else ""
            if identity_only and namespace in _CHANGE_OWNER_EXCLUDED_NAMESPACES:
                continue
            if value not in seen:
                seen.add(value)
                values.append(value)
    return values


def _target_identifier_owner_index(data: TargetGraphData) -> dict[str, set[str]]:
    owners: dict[str, set[str]] = defaultdict(set)
    for target_id, row in data.nodes_by_id.items():
        for value in _target_identifier_values(row, identity_only=True):
            owners[value].add(target_id)
    return owners


def _target_change_where(row: dict[str, str]) -> str:
    sources = _split_pipe(str(row.get("source_namespaces", "")).replace(",", "|"))
    namespaces = _split_pipe(row.get("id_namespaces", ""))
    parts: list[str] = []
    if sources:
        parts.append("source streams: " + ", ".join(sources))
    if namespaces:
        parts.append("identifier namespaces: " + ", ".join(namespaces))
    return "; ".join(parts) or "not captured in app graph"


def _target_related_owner_text(
    row: dict[str, str],
    owner_index: dict[str, set[str]],
    owner_rows: dict[str, dict[str, str]],
    current_target_id: str,
) -> str:
    target_type = row.get("target_type", "")
    related_same_type: list[str] = []
    related_other_type: list[str] = []
    seen: set[str] = set()
    for value in _target_identifier_values(row, identity_only=True):
        for owner in sorted(owner_index.get(value, set())):
            if owner == current_target_id or owner in seen:
                continue
            seen.add(owner)
            owner_type = owner_rows.get(owner, {}).get("target_type", "")
            if owner_type == target_type:
                related_same_type.append(owner)
            else:
                related_other_type.append(owner)
    related = related_same_type or related_other_type
    if not related:
        return ""
    shown = related[:3]
    suffix = "" if len(related) <= 3 else f" +{len(related) - 3} more"
    return ", ".join(shown) + suffix


def _target_added_reason(
    row: dict[str, str],
    baseline_owner_index: dict[str, set[str]],
    baseline_nodes_by_id: dict[str, dict[str, str]],
) -> str:
    target_id = row.get("target_id", "")
    previous_owner_text = _target_related_owner_text(row, baseline_owner_index, baseline_nodes_by_id, target_id)
    if previous_owner_text:
        return (
            "Identifier evidence existed in the baseline under "
            f"{previous_owner_text}; current harmonization re-keyed or split the target node."
        )

    target_type = row.get("target_type", "")
    sources = {s.lower() for s in _split_pipe(str(row.get("source_namespaces", "")).replace(",", "|"))}
    namespaces = set(_split_pipe(row.get("id_namespaces", "")))
    if target_type == "transcript":
        return "New transcript node from the current Ensembl/RefSeq transcript snapshot or transcript-to-protein expansion."
    if target_type == "protein":
        if "UniRef100" in namespaces:
            return "New protein/isoform node from current UniProt plus UniRef100 grouping evidence."
        return "New protein node from current UniProt, RefSeq, Ensembl, or NodeNorm-backed protein harmonization."
    if target_type == "gene":
        if {"hgnc", "ncbi", "ensembl"} & sources:
            return "New gene node from current HGNC/NCBI/Ensembl source harmonization."
        return "New gene node in the current target harmonizer output."
    return "New target node in the current harmonized app graph."


def _target_removed_reason(
    row: dict[str, str],
    current_owner_index: dict[str, set[str]],
    current_nodes_by_id: dict[str, dict[str, str]],
) -> str:
    target_id = row.get("target_id", "")
    current_owner_text = _target_related_owner_text(row, current_owner_index, current_nodes_by_id, target_id)
    if current_owner_text:
        return (
            "Identifier evidence is now represented under current target node(s) "
            f"{current_owner_text}; likely re-keyed, merged, or regrouped during harmonization."
        )

    target_type = row.get("target_type", "")
    sources = _split_pipe(str(row.get("source_namespaces", "")).replace(",", "|"))
    source_text = ", ".join(sources) if sources else "previous source"
    if target_type == "transcript":
        return f"Baseline transcript was not reproduced by the current Ensembl/RefSeq snapshot ({source_text})."
    if target_type == "protein":
        return f"Baseline protein/isoform was not reproduced by the current UniProt/RefSeq/Ensembl harmonization ({source_text})."
    if target_type == "gene":
        return f"Baseline gene was not reproduced by the current HGNC/NCBI/Ensembl harmonization ({source_text})."
    return f"Baseline target was not reproduced by the current source refresh ({source_text})."


def _target_change_drivers(
    source_version_rows: list[dict[str, Any]],
    type_delta: list[dict[str, Any]],
    namespace_delta: list[dict[str, Any]],
    predicate_delta: list[dict[str, Any]],
) -> list[dict[str, Any]]:
    source_rows = {row["source"]: row for row in source_version_rows}
    type_counts = {row["target_type"]: int(row.get("delta") or 0) for row in type_delta}
    ns_counts = {row["namespace"]: int(row.get("delta") or 0) for row in namespace_delta}
    pred_counts = {row["predicate"]: int(row.get("delta") or 0) for row in predicate_delta}

    def version_text(source: str) -> str:
        row = source_rows.get(source, {})
        old = row.get("old_version") or "not present"
        new = row.get("new_version") or "not present"
        return f"{old} -> {new}"

    nodenorm_changed = any(
        source_rows.get(source, {}).get("status") != "unchanged"
        for source in ("NodeNorm Gene", "NodeNorm Protein")
    )
    nodenorm_interpretation = (
        "NodeNorm advanced in this release; the NodeNorm namespace delta reflects expanded validator coverage, not a new primary source authority."
        if nodenorm_changed
        else "NodeNorm stayed version-stable here, so it is validator context rather than the main driver of the large added/removed counts."
    )

    drivers = [
        {
            "source": "Ensembl BioMart",
            "version_change": version_text("Ensembl BioMart"),
            "observed_delta": "; ".join([
                _format_signed_delta(ns_counts.get("ENSEMBL", 0), "ENSEMBL identifiers"),
                _format_signed_delta(type_counts.get("transcript", 0), "transcript nodes"),
                _format_signed_delta(pred_counts.get("biolink:transcribed_to", 0), "gene->transcript edges"),
            ]),
            "interpretation": "Likely largest contributor to the transcript expansion and ENSEMBL identifier growth after moving to the newer Ensembl snapshot.",
        },
        {
            "source": "UniProt",
            "version_change": version_text("UniProt"),
            "observed_delta": "; ".join([
                _format_signed_delta(ns_counts.get("UniProtKB", 0), "UniProtKB identifiers"),
                _format_signed_delta(ns_counts.get("UniRef100", 0), "UniRef100 identifiers"),
                _format_signed_delta(type_counts.get("protein", 0), "protein nodes"),
                _format_signed_delta(pred_counts.get("biolink:translates_to", 0), "transcript->protein edges"),
            ]),
            "interpretation": "Likely largest contributor to protein and isoform growth after the UniProt release advanced and UniRef100 enrichment completed.",
        },
        {
            "source": "RefSeq (NCBI)",
            "version_change": version_text("RefSeq (NCBI)"),
            "observed_delta": _format_signed_delta(ns_counts.get("RefSeq", 0), "RefSeq identifiers"),
            "interpretation": "Contributes transcript/protein identifier support and can change mapping support even when target identity is driven by Ensembl or UniProt.",
        },
        {
            "source": "NCBI Gene Info / HGNC",
            "version_change": f"NCBI {version_text('NCBI Gene Info')}; HGNC {version_text('HGNC')}",
            "observed_delta": "; ".join([
                _format_signed_delta(ns_counts.get("NCBIGene", 0), "NCBI Gene identifiers"),
                _format_signed_delta(ns_counts.get("HGNC", 0), "HGNC identifiers"),
                _format_signed_delta(type_counts.get("gene", 0), "gene nodes"),
            ]),
            "interpretation": "Gene-space changed only slightly; most of the release movement is not from new genes.",
        },
        {
            "source": "NodeNorm",
            "version_change": f"Gene {version_text('NodeNorm Gene')}; Protein {version_text('NodeNorm Protein')}",
            "observed_delta": "; ".join([
                _format_signed_delta(ns_counts.get("NodeNorm", 0), "NodeNorm-backed identifiers"),
                _format_signed_delta(ns_counts.get("UMLS", 0), "UMLS identifiers"),
            ]),
            "interpretation": nodenorm_interpretation,
        },
    ]
    return drivers


def compute_target_version_diff(
    current: TargetGraphData,
    baseline: TargetGraphData,
) -> dict[str, Any]:
    """Compute delta between two versioned target datasets."""
    current_ids = set(current.nodes_by_id.keys())
    baseline_ids = set(baseline.nodes_by_id.keys())

    added_ids = current_ids - baseline_ids
    removed_ids = baseline_ids - current_ids
    shared_ids = current_ids & baseline_ids

    # Classify added/removed by type
    added_by_type = Counter(
        current.nodes_by_id[tid].get("target_type", "unknown") for tid in added_ids
    )
    removed_by_type = Counter(
        baseline.nodes_by_id[tid].get("target_type", "unknown") for tid in removed_ids
    )

    # Detect field changes on shared targets
    symbol_changes: list[dict[str, str]] = []
    name_changes: list[dict[str, str]] = []
    id_changes: list[dict[str, str]] = []
    canonical_changes: list[dict[str, str]] = []

    for tid in shared_ids:
        cur = current.nodes_by_id[tid]
        base = baseline.nodes_by_id[tid]

        cur_sym = cur.get("symbol", "")
        base_sym = base.get("symbol", "")
        if cur_sym != base_sym:
            symbol_changes.append({
                "target_id": tid,
                "target_type": cur.get("target_type", ""),
                "old_symbol": base_sym,
                "new_symbol": cur_sym,
            })

        cur_name = cur.get("name", "")
        base_name = base.get("name", "")
        if cur_name != base_name:
            name_changes.append({
                "target_id": tid,
                "target_type": cur.get("target_type", ""),
                "symbol": cur_sym or base_sym,
                "old_name": base_name,
                "new_name": cur_name,
            })

        cur_ids = cur.get("ids", "")
        base_ids = base.get("ids", "")
        if cur_ids != base_ids:
            id_changes.append({
                "target_id": tid,
                "target_type": cur.get("target_type", ""),
                "symbol": cur_sym or base_sym,
                "old_ids": base_ids,
                "new_ids": cur_ids,
            })

        cur_canon = cur.get("canonical_status", "")
        base_canon = base.get("canonical_status", "")
        if cur_canon != base_canon:
            canonical_changes.append({
                "target_id": tid,
                "target_type": cur.get("target_type", ""),
                "symbol": cur_sym or base_sym,
                "old_canonical_status": base_canon,
                "new_canonical_status": cur_canon,
            })

    # Edge diffs
    current_edge_keys = {_target_edge_key(e) for e in current.relation_edges}
    baseline_edge_keys = {_target_edge_key(e) for e in baseline.relation_edges}
    edges_added = len(current_edge_keys - baseline_edge_keys)
    edges_removed = len(baseline_edge_keys - current_edge_keys)

    # Distributions
    baseline_type_dist = _target_type_distribution(baseline)
    current_type_dist = _target_type_distribution(current)
    baseline_ns_dist = _target_namespace_distribution(baseline)
    current_ns_dist = _target_namespace_distribution(current)

    # Predicate distribution
    baseline_pred_dist = Counter(e.get("predicate", "unknown") for e in baseline.relation_edges)
    current_pred_dist = Counter(e.get("predicate", "unknown") for e in current.relation_edges)

    # Data layer summary
    data_layer_changes = [
        {
            "layer": "Target nodes",
            "old_count": len(baseline.nodes),
            "new_count": len(current.nodes),
        },
        {
            "layer": "Genes",
            "old_count": len(baseline.type_index.get("gene", [])),
            "new_count": len(current.type_index.get("gene", [])),
        },
        {
            "layer": "Proteins",
            "old_count": len(baseline.type_index.get("protein", [])),
            "new_count": len(current.type_index.get("protein", [])),
        },
        {
            "layer": "Transcripts",
            "old_count": len(baseline.type_index.get("transcript", [])),
            "new_count": len(current.type_index.get("transcript", [])),
        },
        {
            "layer": "Relation edges",
            "old_count": len(baseline.relation_edges),
            "new_count": len(current.relation_edges),
        },
    ]
    for row in data_layer_changes:
        row["delta"] = row["new_count"] - row["old_count"]

    # Added/removed target listings (truncated)
    baseline_owner_index = _target_identifier_owner_index(baseline)
    current_owner_index = _target_identifier_owner_index(current)
    added_targets: list[dict[str, str]] = []
    for tid in added_ids:
        row = current.nodes_by_id[tid]
        added_targets.append({
            "target_id": tid,
            "target_type": row.get("target_type", ""),
            "symbol": row.get("symbol", ""),
            "name": row.get("name", ""),
            "where": _target_change_where(row),
            "reason": _target_added_reason(row, baseline_owner_index, baseline.nodes_by_id),
        })
    added_targets = sorted(added_targets, key=lambda r: r["target_id"])

    removed_targets: list[dict[str, str]] = []
    for tid in removed_ids:
        row = baseline.nodes_by_id[tid]
        removed_targets.append({
            "target_id": tid,
            "target_type": row.get("target_type", ""),
            "symbol": row.get("symbol", ""),
            "name": row.get("name", ""),
            "where": _target_change_where(row),
            "reason": _target_removed_reason(row, current_owner_index, current.nodes_by_id),
        })
    removed_targets = sorted(removed_targets, key=lambda r: r["target_id"])

    baseline_version = baseline.manifest.get("target_release_version", "") or baseline.manifest.get("pipeline_version", "unknown")
    current_version = current.manifest.get("target_release_version", "") or current.manifest.get("pipeline_version", "unknown")

    source_version_changes = _source_version_changes(baseline, current)
    type_delta = _distribution_delta(baseline_type_dist, current_type_dist, "target_type")
    namespace_delta = _distribution_delta(baseline_ns_dist, current_ns_dist, "namespace")
    predicate_delta = _distribution_delta(baseline_pred_dist, current_pred_dist, "predicate")

    return {
        "summary": {
            "targets_added": len(added_ids),
            "targets_removed": len(removed_ids),
            "symbol_changes": len(symbol_changes),
            "name_changes": len(name_changes),
            "identifier_changes": len(id_changes),
            "canonical_changes": len(canonical_changes),
            "edges_added": edges_added,
            "edges_removed": edges_removed,
        },
        "data_layer_changes": data_layer_changes,
        "added_targets": added_targets[:500],
        "removed_targets": removed_targets[:500],
        "symbol_changes": symbol_changes[:500],
        "name_changes": name_changes[:500],
        "identifier_changes": id_changes[:500],
        "canonical_changes": canonical_changes[:500],
        "truncation": {
            "added_targets": {"total": len(added_targets), "shown": min(len(added_targets), 500)},
            "removed_targets": {"total": len(removed_targets), "shown": min(len(removed_targets), 500)},
            "symbol_changes": {"total": len(symbol_changes), "shown": min(len(symbol_changes), 500)},
            "name_changes": {"total": len(name_changes), "shown": min(len(name_changes), 500)},
            "identifier_changes": {"total": len(id_changes), "shown": min(len(id_changes), 500)},
            "canonical_changes": {"total": len(canonical_changes), "shown": min(len(canonical_changes), 500)},
        },
        "target_type_distribution": {
            "delta": type_delta,
        },
        "namespace_distribution": {
            "delta": namespace_delta,
        },
        "predicate_distribution": {
            "delta": predicate_delta,
        },
        "source_version_changes": source_version_changes,
        "change_drivers": _target_change_drivers(
            source_version_changes,
            type_delta,
            namespace_delta,
            predicate_delta,
        ),
        "change_driver_note": "Source attribution is inferred from source-version changes plus namespace, node-type, and predicate deltas in the final app graph. Treat it as release-diff context, not a formal upstream changelog.",
        "added_by_type": dict(added_by_type.most_common()),
        "removed_by_type": dict(removed_by_type.most_common()),
        "baseline_version": baseline_version,
        "current_version": current_version,
    }


# ---------------------------------------------------------------------------
# Pharos TDL (from live target_graph ArangoDB) + STRING PPI
# ---------------------------------------------------------------------------

_log = logging.getLogger(__name__)

TDL_COLOR_MAP = {
    "Tclin": "#1f77b4",
    "Tchem": "#2ca02c",
    "Tbio": "#ff7f0e",
    "Tdark": "#d62728",
    "Unknown": "#999999",
}

# AQL mirrors the current_tdls graph view in target_graph.yaml
_CURRENT_TDLS_AQL = """\
FOR pro IN `Protein`
  RETURN {
    id: pro.id,
    uniprot_id: pro.uniprot_id,
    tdl: pro.tdl,
    symbol: pro.symbol,
    name: pro.name,
    idg_family: pro.idg_family,
    canonical_isoform_status: pro.canonical_isoform_status,
    tdl_ligand_count: pro.tdl_meta.tdl_ligand_count,
    tdl_drug_count: pro.tdl_meta.tdl_drug_count,
    tdl_go_term_count: pro.tdl_meta.tdl_go_term_count,
    tdl_generif_count: pro.tdl_meta.tdl_generif_count,
    tdl_pm_score: pro.tdl_meta.tdl_pm_score,
    tdl_antibody_count: pro.tdl_meta.tdl_antibody_count
  }
"""


class PharosTDLData:
    """In-memory index of Pharos TDL data from the live target_graph ArangoDB."""

    __slots__ = (
        "by_ifx_id",
        "by_uniprot",
        "by_symbol",
        "canonical_by_symbol",
        "tdl_counts",
        "loaded",
        "source",
    )

    def __init__(self) -> None:
        self.by_ifx_id: dict[str, dict] = {}
        self.by_uniprot: dict[str, dict] = {}
        self.by_symbol: dict[str, list[dict]] = defaultdict(list)
        self.canonical_by_symbol: dict[str, dict] = {}
        self.tdl_counts: Counter = Counter()
        self.loaded: bool = False
        self.source: str = ""


_pharos_tdl: PharosTDLData | None = None
_pharos_tdl_lock = threading.Lock()


def load_pharos_tdl_data(db) -> PharosTDLData:
    """Load Pharos TDL data from the live target_graph ArangoDB.

    Runs the ``current_tdls`` AQL query against the Protein collection and
    builds in-memory indexes by IFX ID, UniProt accession, and gene symbol.

    Parameters
    ----------
    db : arango.database.StandardDatabase
        An ArangoDB connection to the ``target_graph`` database,
        typically obtained via ``get_db("target_graph")``.
    """
    global _pharos_tdl
    if _pharos_tdl is not None and _pharos_tdl.loaded:
        return _pharos_tdl
    with _pharos_tdl_lock:
        if _pharos_tdl is not None and _pharos_tdl.loaded:
            return _pharos_tdl

        db_name = getattr(db, "name", "target_graph")
        _log.info("Loading Pharos TDL from ArangoDB: %s", db_name)

        data = PharosTDLData()
        data.source = f"arangodb:{db_name}/Protein"

        try:
            cursor = db.aql.execute(_CURRENT_TDLS_AQL, batch_size=5000)
        except Exception as exc:
            raise HTTPException(
                status_code=503,
                detail=f"Cannot query target_graph ArangoDB ({db_name}): {exc}. "
                       "Ensure ArangoDB is reachable and credentials are configured.",
            ) from exc
        for row in cursor:
            ifx_id = (row.get("id") or "").strip()
            uniprot = (row.get("uniprot_id") or "").strip()
            symbol = (row.get("symbol") or "").strip()
            tdl = (row.get("tdl") or "").strip() or "Unknown"

            entry = {
                "id": ifx_id,
                "uniprot_id": uniprot,
                "symbol": symbol,
                "tdl": tdl,
                "name": (row.get("name") or "").strip(),
                "canonical_isoform_status": (row.get("canonical_isoform_status") or "").strip(),
                "idg_family": (row.get("idg_family") or "").strip(),
                "tdl_drug_count": row.get("tdl_drug_count") or 0,
                "tdl_ligand_count": row.get("tdl_ligand_count") or 0,
                "tdl_go_term_count": row.get("tdl_go_term_count") or 0,
                "tdl_generif_count": row.get("tdl_generif_count") or 0,
                "tdl_pm_score": row.get("tdl_pm_score") or 0,
                "tdl_antibody_count": row.get("tdl_antibody_count") or 0,
            }
            if ifx_id:
                data.by_ifx_id[ifx_id] = entry
            if uniprot:
                data.by_uniprot[uniprot] = entry
            if symbol:
                data.by_symbol[symbol].append(entry)
            if tdl:
                data.tdl_counts[tdl] += 1

        # Build canonical_by_symbol: prefer canonical isoform, then first entry
        for sym, entries in data.by_symbol.items():
            canonical = None
            for e in entries:
                status = e.get("canonical_isoform_status", "").lower()
                if "canonical" in status:
                    canonical = e
                    break
            data.canonical_by_symbol[sym] = canonical or entries[0]

        data.loaded = True
        _pharos_tdl = data
        _log.info(
            "Pharos TDL loaded: %d proteins, %d UniProt, %d symbols — %s",
            len(data.by_ifx_id),
            len(data.by_uniprot),
            len(data.by_symbol),
            dict(data.tdl_counts.most_common()),
        )
        return data


def invalidate_pharos_tdl_cache() -> None:
    """Clear the cached Pharos TDL data so the next call re-queries ArangoDB."""
    global _pharos_tdl
    with _pharos_tdl_lock:
        _pharos_tdl = None


def match_targets_to_tdl(
    data: TargetGraphData, pharos: PharosTDLData,
) -> dict[str, Any]:
    """Match protein nodes in the target graph to Pharos TDL classifications.

    Matching strategy (in order):
    1. IFXProtein target_id → ``pharos.by_ifx_id`` (direct pipeline ID match)
    2. UniProt ID from node ``ids`` field → ``pharos.by_uniprot``
    3. Gene symbol fallback → ``pharos.canonical_by_symbol``
    """
    tdl_dist: Counter = Counter()
    matched = 0
    unmatched = 0
    targets: list[dict[str, Any]] = []

    for node in data.nodes:
        if (node.get("target_type") or "").lower() != "protein":
            continue

        symbol = (node.get("symbol") or "").strip()
        target_id = node.get("target_id", "")
        ids_field = node.get("ids") or ""

        tdl_entry: dict | None = None

        # Strategy 1: direct IFXProtein ID match
        if target_id:
            tdl_entry = pharos.by_ifx_id.get(target_id)

        # Strategy 2: match by UniProt ID
        if not tdl_entry:
            for raw_id in ids_field.split("|"):
                raw_id = raw_id.strip()
                if raw_id.startswith("UniProtKB:"):
                    up_acc = raw_id[len("UniProtKB:"):]
                    tdl_entry = pharos.by_uniprot.get(up_acc)
                    if tdl_entry:
                        break

        # Strategy 3: match by symbol
        if not tdl_entry and symbol:
            tdl_entry = pharos.canonical_by_symbol.get(symbol)

        if tdl_entry:
            tdl = tdl_entry.get("tdl", "Unknown") or "Unknown"
            matched += 1
        else:
            tdl = "Unknown"
            unmatched += 1

        tdl_dist[tdl] += 1
        targets.append({
            "target_id": target_id,
            "symbol": symbol,
            "tdl": tdl,
            "name": node.get("name", ""),
            "canonical_status": node.get("canonical_status", ""),
            "pharos_name": tdl_entry.get("name", "") if tdl_entry else "",
            "idg_family": tdl_entry.get("idg_family", "") if tdl_entry else "",
        })

    # Sort TDL distribution in standard order
    tdl_order = ["Tclin", "Tchem", "Tbio", "Tdark", "Unknown"]
    ordered_dist = {k: tdl_dist.get(k, 0) for k in tdl_order if tdl_dist.get(k, 0) > 0}

    return {
        "tdl_distribution": ordered_dist,
        "matched": matched,
        "unmatched": unmatched,
        "total_proteins": matched + unmatched,
        "source": pharos.source,
        "targets": targets,
    }


def lookup_tdl_by_symbols(
    pharos: PharosTDLData, symbols: list[str],
) -> dict[str, Any]:
    """Look up Pharos TDL classification for a list of gene symbols.

    Each symbol is matched against ``pharos.canonical_by_symbol``.  Returns
    per-gene results with TDL metrics plus an aggregate TDL distribution.
    """
    results: list[dict[str, Any]] = []
    not_found_symbols: list[str] = []
    tdl_dist: Counter = Counter()

    for symbol in symbols:
        entry = pharos.canonical_by_symbol.get(symbol)
        if entry:
            tdl = entry.get("tdl", "Unknown") or "Unknown"
            tdl_dist[tdl] += 1
            results.append({
                "symbol": symbol,
                "tdl": tdl,
                "name": entry.get("name", ""),
                "uniprot_id": entry.get("uniprot_id", ""),
                "idg_family": entry.get("idg_family", ""),
                "tdl_drug_count": entry.get("tdl_drug_count", 0),
                "tdl_ligand_count": entry.get("tdl_ligand_count", 0),
                "tdl_go_term_count": entry.get("tdl_go_term_count", 0),
                "tdl_generif_count": entry.get("tdl_generif_count", 0),
                "tdl_pm_score": entry.get("tdl_pm_score", 0),
                "tdl_antibody_count": entry.get("tdl_antibody_count", 0),
            })
        else:
            not_found_symbols.append(symbol)

    tdl_order = ["Tclin", "Tchem", "Tbio", "Tdark"]
    ordered_dist = {k: tdl_dist.get(k, 0) for k in tdl_order if tdl_dist.get(k, 0) > 0}

    return {
        "results": results,
        "tdl_distribution": ordered_dist,
        "found": len(results),
        "not_found": len(not_found_symbols),
        "not_found_symbols": not_found_symbols,
    }


# ── STRING PPI ──

STRING_NETWORK_URL = "https://string-db.org/api/tsv/network"


def fetch_string_ppi(
    genes: list[str],
    required_score: int = 700,
    species: int = 9606,
) -> list[dict[str, Any]]:
    """Query STRING for protein-protein interactions among a gene list.

    Returns a list of edge dicts with ``preferredName_A``, ``preferredName_B``,
    and ``score``.
    """
    import requests
    from requests.adapters import HTTPAdapter, Retry

    if not genes:
        return []

    session = requests.Session()
    retries = Retry(
        total=4,
        backoff_factor=0.8,
        status_forcelist=[429, 500, 502, 503, 504],
        allowed_methods=["GET"],
    )
    session.mount("https://", HTTPAdapter(max_retries=retries))

    identifiers = "%0d".join(genes)
    params = {
        "identifiers": identifiers,
        "species": species,
        "required_score": required_score,
        "caller_identity": "odin_qa_browser",
    }
    resp = session.get(STRING_NETWORK_URL, params=params, timeout=60)
    resp.raise_for_status()

    edges: list[dict[str, Any]] = []
    reader = csv.DictReader(io.StringIO(resp.text), delimiter="\t")
    for row in reader:
        score_raw = row.get("score", "0")
        try:
            score = float(score_raw)
        except (ValueError, TypeError):
            score = 0.0
        # STRING sometimes returns 0–1 scale; normalize to 0–1000
        if score <= 1.0:
            score = score * 1000.0
        edges.append({
            "preferredName_A": row.get("preferredName_A", ""),
            "preferredName_B": row.get("preferredName_B", ""),
            "score": round(score),
            "stringId_A": row.get("stringId_A", ""),
            "stringId_B": row.get("stringId_B", ""),
        })

    return edges


# ── Cytoscape payload builder ──


def build_pharos_ppi_payload(
    pharos: PharosTDLData,
    genes: list[str],
    ppi_edges: list[dict[str, Any]],
) -> dict[str, Any]:
    """Build a Cytoscape.js-compatible payload for the PPI network.

    Nodes are colored by TDL class; edges carry STRING scores.
    """
    gene_set = set(genes)
    tdl_map: dict[str, str] = {}
    nodes: list[dict] = []

    for g in genes:
        entry = pharos.canonical_by_symbol.get(g)
        tdl = (entry.get("tdl") if entry else None) or "Unknown"
        tdl_map[g] = tdl
        nodes.append({
            "data": {
                "id": g,
                "label": g,
                "tdl": tdl,
                "color": TDL_COLOR_MAP.get(tdl, TDL_COLOR_MAP["Unknown"]),
                "pharos_name": (entry.get("name", "") if entry else ""),
                "idg_family": (entry.get("idg_family", "") if entry else ""),
            },
        })

    edges: list[dict] = []
    for e in ppi_edges:
        a = e.get("preferredName_A", "")
        b = e.get("preferredName_B", "")
        if a in gene_set and b in gene_set:
            score = e.get("score", 0)
            edges.append({
                "data": {
                    "id": f"{a}--{b}",
                    "source": a,
                    "target": b,
                    "score": score,
                    "width": 1 + (score / 1000) * 6,
                },
            })

    tdl_dist: Counter = Counter(tdl_map.values())
    tdl_order = ["Tclin", "Tchem", "Tbio", "Tdark", "Unknown"]
    ordered_dist = {k: tdl_dist.get(k, 0) for k in tdl_order if tdl_dist.get(k, 0) > 0}

    return {
        "elements": nodes + edges,
        "tdl_distribution": ordered_dist,
        "gene_count": len(nodes),
        "edge_count": len(edges),
        "tdl_map": tdl_map,
    }
