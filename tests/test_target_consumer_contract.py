import csv
from pathlib import Path

import pytest

from src.id_resolvers import sqlite_cache_resolver
from src.id_resolvers.target_consumer_resolver import (
    TargetConsumerGeneResolver,
    TargetConsumerProteinResolver,
    mapping_curie,
)
from src.id_resolvers.target_consumer_pharos_resolver import TargetConsumerPharosCanonicalResolver
from src.interfaces.id_resolver import IdMatch, IdResolver, NoMatchBehavior, MultiMatchBehavior
from src.input_adapters.target_graph.target_consumer import (
    TargetConsumerGeneAdapter,
    TargetConsumerProteinAdapter,
    TargetConsumerRelationshipAdapter,
)
from src.models.gene import Gene
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.transcript import GeneTranscriptEdge
from src.output_adapters.arango_output_adapter import ArangoOutputAdapter


def test_pharos_gene_route_retains_only_unmatched_ifx_genes():
    resolver = object.__new__(TargetConsumerPharosCanonicalResolver)
    resolver.allow_unmatched_prefixes = frozenset({"IFXGene"})
    resolver.no_match_behavior = NoMatchBehavior.Allow
    resolver.multi_match_behavior = MultiMatchBehavior.All
    entity_map = {
        IdResolver.MatchKeys.matched: {},
        IdResolver.MatchKeys.newborns: {},
        IdResolver.MatchKeys.unmatched: {
            "IFXGene:GENE": Gene(id="IFXGene:GENE"),
            "NCBIGene:123": Gene(id="NCBIGene:123"),
        },
    }
    assert list(resolver.parse_entity_map(entity_map)) == ["IFXGene:GENE"]


class _Release:
    snapshot_id = "ifx_harmonizers:targets:2.5.1"
    version = "2.5.1"

    def __init__(self, root: Path):
        self.root = root
        self.manifest = {
            "files": [{"path": "app_graph/target_mappings.tsv", "sha256": "test-sha"}]
        }

    def file(self, name):
        path = self.root / name
        if not path.is_file():
            raise FileNotFoundError(path)
        return path

    def version_info(self):
        return DatasourceVersionInfo(version="2.5.1", version_date=None, download_date=None)


class _Release252(_Release):
    snapshot_id = "ifx_harmonizers:targets:2.5.2"
    version = "2.5.2"

    def version_info(self):
        return DatasourceVersionInfo(version="2.5.2", version_date=None, download_date=None)


def _write_tsv(path, columns, rows):
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def _mapping(ifx_id, entity_type, exact, base, namespace):
    return {
        "ifx_id": ifx_id,
        "entity_type": entity_type,
        "equivalent_id": exact,
        "equivalent_id_base": base,
        "namespace": namespace,
        "source": "uniprot" if namespace == "UniProtKB" else "ensembl",
        "source_version": "116",
        "confidence": "high",
    }


def test_252_match_types_keep_all_lookups_grouped_on_xref(tmp_path, monkeypatch):
    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    rows = [
        {**_mapping("IFXGene:EGFR", "gene", "ENSG00000146648", "ENSG00000146648", "ENSEMBL"), "match_type": "exact"},
        {**_mapping("IFXGene:EGFR", "gene", "EGFR", "EGFR", "Symbol"), "match_type": "symbol"},
        {**_mapping("IFXGene:EGFR", "gene", "ERBB", "ERBB", "Symbol"), "match_type": "previous_symbol"},
        {**_mapping("IFXGene:EGFR", "gene", "123", "123", "Orphanet"), "match_type": "exact"},
    ]
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(rows[0]), rows)
    resolver = TargetConsumerGeneResolver(_Release252(tmp_path), types=["Gene"])

    match = resolver.resolve_internal([Gene(id="Symbol:EGFR")])["Symbol:EGFR"][0]
    assert match.match == "IFXGene:EGFR"
    assert match.context == ["symbol"]
    result = resolver.get_merged_map([Gene(id="Symbol:EGFR")], {"Symbol:EGFR": [match]})
    node = result[IdResolver.MatchKeys.matched]["Symbol:EGFR"]
    assert {xref.id_str() for xref in node.xref} == set().union(*node.xref_match_types.values())
    assert node.xref_match_types == {
        "exact": ["ENSEMBL:ENSG00000146648", "IFXGene:EGFR", "orphanet:123"],
        "previous_symbol": ["Symbol:ERBB"],
        "symbol": ["Symbol:EGFR"],
    }
    assert ArangoOutputAdapter.clean_dict(object.__new__(ArangoOutputAdapter), node, convert_dates=True)["xref_match_types"] == node.xref_match_types


def test_252_requires_known_match_type(tmp_path, monkeypatch):
    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    row = {**_mapping("IFXGene:EGFR", "gene", "EGFR", "EGFR", "Symbol"), "match_type": "guess"}
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(row), [row])
    with pytest.raises(ValueError, match="unknown match_type"):
        TargetConsumerGeneResolver(_Release252(tmp_path), types=["Gene"])


def test_252_rejects_conflicting_types_for_one_lookup(tmp_path, monkeypatch):
    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    row = {**_mapping("IFXGene:EGFR", "gene", "EGFR", "EGFR", "Symbol"), "match_type": "symbol"}
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(row), [
        row, {**row, "match_type": "previous_symbol"},
    ])
    with pytest.raises(ValueError, match="conflicting match types"):
        TargetConsumerGeneResolver(_Release252(tmp_path), types=["Gene"])


def test_252_pharos_projection_groups_routed_and_source_lookup_types(tmp_path, monkeypatch):
    from src.models.protein import Protein

    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    source = _Release252(tmp_path)
    source.manifest["files"].extend([
        {"path": "app_graph/target_relationships.tsv", "sha256": "relationship-sha"},
        {"path": "app_graph/target_protein_payload.tsv", "sha256": "payload-sha"},
    ])
    rows = [
        {**_mapping("IFXProtein:CANON", "protein", "P00533-1", "P00533-1", "UniProtKB"), "match_type": "exact"},
        {**_mapping("IFXProtein:CANON", "protein", "P00533", "P00533", "UniProtKB"), "match_type": "uniprot_entry"},
        {**_mapping("IFXProtein:CANON", "protein", "QOLD", "QOLD", "UniProtKB"), "match_type": "secondary_accession"},
        {**_mapping("IFXGene:GENE", "gene", "EGFR", "EGFR", "Symbol"), "match_type": "symbol"},
    ]
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(rows[0]), rows)
    _write_tsv(tmp_path / "app_graph/target_protein_payload.tsv", ["ifx_id", "is_canonical"], [
        {"ifx_id": "IFXProtein:CANON", "is_canonical": "True"},
    ])
    relations = [
        {"subject_id": "IFXGene:GENE", "object_id": "IFXProtein:CANON", "relation_kind": "gene_to_protein_resolved_parent", "status": "selected"},
        {"subject_id": "IFXGene:GENE", "object_id": "IFXTranscript:TX", "relation_kind": "gene_to_transcript", "status": "selected"},
        {"subject_id": "IFXTranscript:TX", "object_id": "IFXProtein:CANON", "relation_kind": "transcript_to_protein", "status": "selected"},
        {"subject_id": "Symbol:SHARED", "object_id": "IFXGene:GENE", "relation_kind": "shared_identifier_claim", "status": "source_claim"},
    ]
    _write_tsv(tmp_path / "app_graph/target_relationships.tsv", list(relations[0]), relations)

    resolver = TargetConsumerPharosCanonicalResolver(source, types=["Protein"], multi_match_behavior="All")
    match = resolver.resolve_internal([Protein(id="Symbol:EGFR")])["Symbol:EGFR"][0]
    assert match.match == "IFXProtein:CANON"
    assert match.context == ["symbol"]
    node = resolver.parse_flat_node_list_from_map(
        resolver.get_merged_map([Protein(id="Symbol:EGFR")], {"Symbol:EGFR": [match]})
    )[0]
    assert node.xref_match_types["symbol"] == ["Symbol:EGFR"]
    assert node.xref_match_types["secondary_accession"] == ["UniProtKB:QOLD"]
    assert node.xref_match_types["uniprot_entry"] == ["UniProtKB:P00533"]
    assert "IFXGene:GENE" in node.xref_match_types["routed"]
    assert resolver.resolve_internal([Protein(id="Symbol:SHARED")]) == {}


def test_252_payload_requires_matching_release_version(tmp_path):
    source = _Release252(tmp_path)
    row = {
        "ifx_id": "IFXProtein:ONE", "release_version": "v2.5.1", "aliases": "",
        "is_canonical": "True", "createdAt": "2026-04-21T10:28:24",
        "updatedAt": "2026-10-19T10:28:24",
    }
    _write_tsv(tmp_path / "app_graph/target_protein_payload.tsv", list(row), [row])
    with pytest.raises(ValueError, match="expected 'v2.5.2'"):
        list(TargetConsumerProteinAdapter(source).get_all())
    row["release_version"] = "v2.5.2"
    _write_tsv(tmp_path / "app_graph/target_protein_payload.tsv", list(row), [row])
    nodes = list(TargetConsumerProteinAdapter(source).get_all())
    assert nodes[0][0].target_payload["release_version"] == "v2.5.2"


def test_contract_resolver_uses_only_mapping_rows_and_version_free_lookup(tmp_path, monkeypatch):
    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    path = tmp_path / "app_graph/target_mappings.tsv"
    rows = [
        _mapping("IFXGene:EGFR", "gene", "ENSG00000146648", "ENSG00000146648", "ENSEMBL"),
        _mapping("IFXProtein:CANON", "protein", "P00533-1", "P00533-1", "UniProtKB"),
        _mapping("IFXProtein:CANON", "protein", "ENSP00000275493.2", "ENSP00000275493", "ENSEMBL"),
        _mapping("IFXProtein:ALT", "protein", "P00533-2", "P00533-2", "UniProtKB"),
    ]
    _write_tsv(path, list(rows[0]), rows)
    source = _Release(tmp_path)
    genes = TargetConsumerGeneResolver(source, types=["Gene"])
    proteins = TargetConsumerProteinResolver(source, types=["Protein"])

    assert genes.resolve_internal([Gene(id="ENSEMBL:ENSG00000146648")])["ENSEMBL:ENSG00000146648"][0].match == "IFXGene:EGFR"
    from src.models.protein import Protein
    assert proteins.resolve_internal([Protein(id="ENSEMBL:ENSP00000275493")])["ENSEMBL:ENSP00000275493"][0].match == "IFXProtein:CANON"
    assert proteins.resolve_internal([Protein(id="UniProtKB:P00533-2")])["UniProtKB:P00533-2"][0].match == "IFXProtein:ALT"
    assert proteins.resolve_internal([Protein(id="UniProtKB:P00533")]) == {}
    assert genes.resolve_internal([Gene(id="Symbol:EGFR")]) == {}


def test_contract_resolver_rejects_ambiguous_base(tmp_path, monkeypatch):
    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    rows = [
        _mapping("IFXProtein:ONE", "protein", "ENSP123.1", "ENSP123", "ENSEMBL"),
        _mapping("IFXProtein:TWO", "protein", "ENSP123.2", "ENSP123", "ENSEMBL"),
    ]
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(rows[0]), rows)
    with pytest.raises(ValueError, match="mapping conflict"):
        TargetConsumerProteinResolver(_Release(tmp_path), types=["Protein"])


def test_contract_payload_and_selected_edge_keep_source_rows(tmp_path):
    gene_row = {
        "ifx_id": "IFXGene:EGFR",
        "release_version": "v2.5.1",
        "aliases": "EGFR=approved_symbol|ERBB=previous_symbol",
        "consolidated_symbol": "EGFR",
        "consolidated_description": "epidermal growth factor receptor",
        "createdAt": "2026-04-21T10:28:24",
        "updatedAt": "2026-10-19T10:28:24",
        "Total_Mapping_Ratio": "0.9",
    }
    _write_tsv(tmp_path / "app_graph/target_gene_payload.tsv", list(gene_row), [gene_row])
    node = next(TargetConsumerGeneAdapter(_Release(tmp_path)).get_all())[0]
    assert node.id == "IFXGene:EGFR"
    assert node.symbol == "EGFR"
    assert node.target_payload["aliases"] == "EGFR=approved_symbol|ERBB=previous_symbol"

    selected = {
        "subject_id": "IFXGene:EGFR", "predicate": "biolink:transcribed_to",
        "biolink_predicate": "biolink:transcribed_to", "object_id": "IFXTranscript:ONE",
        "relation_kind": "gene_to_transcript", "status": "selected",
        "evidence_source": "transcript_ids.ensembl_gene_id", "evidence": "gene_to_transcript",
        "confidence": "high", "source_version": "116",
    }
    candidate = {**selected, "object_id": "IFXTranscript:TWO", "status": "candidate"}
    deferred = {**selected, "object_id": "IFXProtein:ONE", "relation_kind": "gene_to_canonical_protein", "predicate": "ifx:gene_has_canonical_product"}
    _write_tsv(tmp_path / "app_graph/target_relationships.tsv", list(selected), [selected, candidate, deferred])
    edges = [edge for batch in TargetConsumerRelationshipAdapter(_Release(tmp_path)).get_all() for edge in batch]
    assert len(edges) == 1
    assert isinstance(edges[0], GeneTranscriptEdge)
    assert {key: edges[0].details[0].to_dict()[key] for key in selected} == selected
    serializer = object.__new__(ArangoOutputAdapter)
    converted_node = serializer.clean_dict(node, convert_dates=True)
    converted_edge = serializer.clean_dict(edges[0], convert_dates=True)
    assert converted_node["target_payload"]["aliases"] == gene_row["aliases"]
    assert {key: converted_edge["details"][0][key] for key in selected} == selected
    schema = ArangoOutputAdapter._introspect_dataclass(GeneTranscriptEdge)
    detail_fields = schema["details"]["fields"]
    assert {"predicate", "relation_kind", "status", "evidence_source", "evidence", "source_version"} <= set(detail_fields)


def test_mapping_curie_preserves_identifier_family():
    assert mapping_curie("HGNC", "HGNC:3236") == "HGNC:3236"
    assert mapping_curie("OMIM", "MIM:131550") == "OMIM:131550"
    assert mapping_curie("IMGT", "IMGT/GENE-DB:TRAC") == "IMGT:TRAC"
    assert mapping_curie("Symbol", "HHC2:066522") == "Symbol:HHC2:066522"
    with pytest.raises(ValueError, match="Unexpected prefix"):
        mapping_curie("UniProtKB", "HGNC:3236")


def test_contract_protein_canonical_filter_keeps_complete_payload(tmp_path):
    rows = [
        {
            "ifx_id": "IFXProtein:CANON", "release_version": "v2.5.1", "aliases": "EGFR=protein_symbol",
            "is_canonical": "True", "uniprot_id": "P00533", "uniprot_isoform": "P00533-1",
            "combined_protein_name": "EGFR", "createdAt": "2026-04-21T10:28:24",
            "updatedAt": "2026-10-19T10:28:24",
        },
        {
            "ifx_id": "IFXProtein:ISOFORM", "release_version": "v2.5.1", "aliases": "EGFR=protein_symbol",
            "is_canonical": "False", "uniprot_id": "P00533", "uniprot_isoform": "P00533-2",
            "combined_protein_name": "EGFR isoform", "createdAt": "2026-04-21T10:28:24",
            "updatedAt": "2026-10-19T10:28:24",
        },
    ]
    _write_tsv(tmp_path / "app_graph/target_protein_payload.tsv", list(rows[0]), rows)
    source = _Release(tmp_path)

    all_nodes = [node for batch in TargetConsumerProteinAdapter(source).get_all() for node in batch]
    canonical_nodes = [node for batch in TargetConsumerProteinAdapter(source, canonical_only=True).get_all() for node in batch]

    assert {node.id for node in all_nodes} == {"IFXProtein:CANON", "IFXProtein:ISOFORM"}
    assert [node.id for node in canonical_nodes] == ["IFXProtein:CANON"]
    assert canonical_nodes[0].target_payload == rows[0]


def test_pharos_canonical_projection_includes_both_gene_paths_and_isoforms(tmp_path, monkeypatch):
    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    source = _Release(tmp_path)
    source.manifest["files"].extend([
        {"path": "app_graph/target_relationships.tsv", "sha256": "relationship-sha"},
        {"path": "app_graph/target_protein_payload.tsv", "sha256": "payload-sha"},
    ])
    mappings = [
        _mapping("IFXProtein:CANON", "protein", "P00533-1", "P00533-1", "UniProtKB"),
        _mapping("IFXProtein:ISO", "protein", "P00533-2", "P00533-2", "UniProtKB"),
        _mapping("IFXProtein:CANON", "protein", "NP_123.2", "NP_123", "RefSeq"),
        _mapping("IFXTranscript:TX", "transcript", "ENST00001.2", "ENST00001", "ENSEMBL"),
        _mapping("IFXGene:GENE", "gene", "ENSG00001.3", "ENSG00001", "ENSEMBL"),
    ]
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(mappings[0]), mappings)
    _write_tsv(tmp_path / "app_graph/target_protein_payload.tsv", ["ifx_id", "is_canonical"], [
        {"ifx_id": "IFXProtein:CANON", "is_canonical": "True"},
        {"ifx_id": "IFXProtein:ISO", "is_canonical": "False"},
    ])
    relations = [
        ("IFXProtein:ISO", "IFXProtein:CANON", "isoform_of_canonical_protein"),
        ("IFXTranscript:TX", "IFXProtein:ISO", "transcript_to_protein"),
        ("IFXGene:GENE", "IFXProtein:ISO", "gene_to_protein_resolved_parent"),
        ("IFXGene:GENE", "IFXTranscript:TX", "gene_to_transcript"),
        ("IFXGene:GENE", "IFXProtein:CANON", "gene_to_canonical_protein"),
        ("P00533", "IFXProtein:CANON", "uniprot_entry_to_displayed_protein"),
    ]
    rows = [dict(subject_id=a, object_id=b, relation_kind=k, status="selected") for a, b, k in relations]
    rows.append({**rows[2], "subject_id": "IFXGene:IGNORED", "status": "candidate"})
    _write_tsv(tmp_path / "app_graph/target_relationships.tsv", list(rows[0]), rows)

    from src.models.protein import Protein
    resolver = TargetConsumerPharosCanonicalResolver(source, types=["Protein"])
    for lookup in ("IFXGene:GENE", "ENSEMBL:ENSG00001", "IFXTranscript:TX",
                   "ENSEMBL:ENST00001", "IFXProtein:ISO", "UniProtKB:P00533-2",
                   "UniProtKB:P00533"):
        match = resolver.resolve_internal([Protein(id=lookup)])[lookup][0]
        assert match.match == "IFXProtein:CANON"
        assert "ENSEMBL:ENSG00001" in match.equivalent_ids
        assert "UniProtKB:P00533-2" in match.equivalent_ids
    assert resolver.resolve_internal([Protein(id="ENSEMBL:ENST00001.2")])["ENSEMBL:ENST00001.2"][0].match == "IFXProtein:CANON"
    fallback = resolver.resolve_internal([Protein(id="ENSEMBL:ENST00001.9")])["ENSEMBL:ENST00001.9"][0]
    assert fallback.match == "IFXProtein:CANON"
    assert "version_free_fallback" in fallback.context
    assert resolver.resolve_internal([Protein(id="RefSeq:NP_123.9")])["RefSeq:NP_123.9"][0].match == "IFXProtein:CANON"
    canonical = Protein(id="IFXProtein:CANON")
    resolved = resolver.get_merged_map([canonical], resolver.resolve_internal([canonical]))
    xrefs = {xref.id_str() for xref in resolved[IdResolver.MatchKeys.matched][canonical.id].xref}
    assert "ENSEMBL:ENST00001" in xrefs
    assert "ENSEMBL:ENST00001.2" not in xrefs
    assert "RefSeq:NP_123" in xrefs
    assert "RefSeq:NP_123.2" not in xrefs
    assert "UniProtKB:P00533-2" in xrefs
    assert resolver.resolve_internal([Protein(id="IFXGene:IGNORED")]) == {}


def test_pharos_fanout_keeps_each_canonical_targets_xrefs_separate():
    from src.models.protein import Protein

    resolver = object.__new__(TargetConsumerPharosCanonicalResolver)
    IdResolver.__init__(resolver, types=["Protein"], multi_match_behavior="All")
    shared = "IFXGene:SHARED"
    matches = [
        IdMatch(input=shared, match="IFXProtein:ONE",
                equivalent_ids=[shared, "UniProtKB:ONE"]),
        IdMatch(input=shared, match="IFXProtein:TWO",
                equivalent_ids=[shared, "UniProtKB:TWO"]),
    ]
    result = resolver.get_merged_map([Protein(id=shared)], {shared: matches})
    nodes = resolver.parse_flat_node_list_from_map(result)
    by_id = {node.id: {xref.id_str() for xref in node.xref} for node in nodes}
    assert by_id == {
        "IFXProtein:ONE": {shared, "UniProtKB:ONE"},
        "IFXProtein:TWO": {shared, "UniProtKB:TWO"},
    }


def test_pharos_exact_version_precedes_ambiguous_base_fallback(tmp_path, monkeypatch):
    from src.models.protein import Protein

    monkeypatch.setattr(sqlite_cache_resolver, "DEFAULT_SQLITE_RESOLVER_CACHE_DIR", tmp_path / "cache")
    source = _Release(tmp_path)
    source.manifest["files"].extend([
        {"path": "app_graph/target_relationships.tsv", "sha256": "relationship-sha"},
        {"path": "app_graph/target_protein_payload.tsv", "sha256": "payload-sha"},
    ])
    mappings = [
        _mapping("IFXProtein:ONE", "protein", "ENSPX.1", "ENSPX", "ENSEMBL"),
        _mapping("IFXProtein:TWO", "protein", "ENSPX.2", "ENSPX", "ENSEMBL"),
    ]
    _write_tsv(tmp_path / "app_graph/target_mappings.tsv", list(mappings[0]), mappings)
    _write_tsv(tmp_path / "app_graph/target_protein_payload.tsv", ["ifx_id", "is_canonical"], [
        {"ifx_id": "IFXProtein:ONE", "is_canonical": "True"},
        {"ifx_id": "IFXProtein:TWO", "is_canonical": "True"},
    ])
    _write_tsv(tmp_path / "app_graph/target_relationships.tsv",
               ["subject_id", "object_id", "relation_kind", "status"], [])

    resolver = TargetConsumerPharosCanonicalResolver(source, types=["Protein"], multi_match_behavior="All")
    exact = resolver.resolve_internal([Protein(id="ENSEMBL:ENSPX.1")])["ENSEMBL:ENSPX.1"]
    fallback = resolver.resolve_internal([Protein(id="ENSEMBL:ENSPX.9")])["ENSEMBL:ENSPX.9"]
    assert {match.match for match in exact} == {"IFXProtein:ONE"}
    assert {match.match for match in fallback} == {"IFXProtein:ONE", "IFXProtein:TWO"}
    for match in fallback:
        assert "version_free_fallback" in match.context
        xrefs = {xref.id_str() for xref in resolver._xrefs_for_match(match, [])}
        assert xrefs == {match.match, "ENSEMBL:ENSPX"}
