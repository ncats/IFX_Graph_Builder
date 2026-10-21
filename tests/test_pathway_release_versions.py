import json

from src.qa_browser import app as qa_app
from types import SimpleNamespace

from src.qa_browser.pathway_id_graph import (
    PathwayGraphData,
    _release_version_from_path,
    build_cross_entity_summary,
)


def _write_graph(root, version):
    graph = root / f"v{version}" / "app_graph"
    graph.mkdir(parents=True)
    # Reproduce the historical producer bug: both manifests claimed 0.1.0.
    (graph / "manifest.json").write_text(json.dumps({"version": "0.1.0"}))
    (graph / "pathway_nodes.tsv").write_text("ncats_pathway_id\n")
    (graph / "pathway_edges.tsv").write_text("source_id\ttarget_id\n")
    return graph


def test_release_version_is_derived_from_version_directory(tmp_path):
    graph = _write_graph(tmp_path, "1.2.3")
    assert _release_version_from_path(graph) == "1.2.3"


def test_bundled_graph_selection_uses_latest_manifest_version(tmp_path):
    old = tmp_path / "v2.1.1"
    old.mkdir()
    (old / "manifest.json").write_text(json.dumps({"target_release_version": "v2.1.1"}))
    (old / "target_nodes.tsv").write_text("target_id\n")

    current = tmp_path / "current"
    current.mkdir()
    (current / "manifest.json").write_text(json.dumps({"target_release_version": "v2.2.0"}))
    (current / "target_nodes.tsv").write_text("target_id\n")

    releases = qa_app._versioned_app_graph_dirs(
        tmp_path, "target_nodes.tsv", ("target_release_version", "version")
    )

    assert releases == [old, current]


def test_bundled_graph_selection_prefers_current_for_same_version(tmp_path):
    versioned = tmp_path / "v2.2.0"
    versioned.mkdir()
    (versioned / "manifest.json").write_text(json.dumps({"target_release_version": "v2.2.0"}))
    (versioned / "target_nodes.tsv").write_text("target_id\n")

    current = tmp_path / "current"
    current.mkdir()
    (current / "manifest.json").write_text(json.dumps({"target_release_version": "v2.2.0"}))
    (current / "target_nodes.tsv").write_text("target_id\n")

    releases = qa_app._versioned_app_graph_dirs(
        tmp_path, "target_nodes.tsv", ("target_release_version", "version")
    )

    assert releases == [current]


def test_bundled_graph_selection_skips_invalid_manifests(tmp_path):
    invalid = tmp_path / "v9.9.9"
    invalid.mkdir()
    (invalid / "manifest.json").write_text("not json")
    (invalid / "target_nodes.tsv").write_text("target_id\n")

    assert qa_app._versioned_app_graph_dirs(
        tmp_path, "target_nodes.tsv", ("target_release_version", "version")
    ) == []


def test_bundled_graph_selection_uses_entity_release_field_precedence(tmp_path):
    older = tmp_path / "v1.0.0"
    older.mkdir()
    (older / "manifest.json").write_text(json.dumps({
        "version": "99.0.0",
        "variant_harmonizer_version": "1.0.0",
    }))
    (older / "variant_nodes.tsv").write_text("variant_id\n")

    newer = tmp_path / "current"
    newer.mkdir()
    (newer / "manifest.json").write_text(json.dumps({
        "version": "1.0.0",
        "variant_harmonizer_version": "1.1.0",
    }))
    (newer / "variant_nodes.tsv").write_text("variant_id\n")

    releases = qa_app._versioned_app_graph_dirs(
        tmp_path,
        "variant_nodes.tsv",
        ("variant_harmonizer_version", "variant_release_version", "version"),
    )

    assert releases == [older, newer]


def test_cross_entity_summary_uses_association_and_edge_indexes():
    pathway = PathwayGraphData()
    pathway.nodes_by_id["IFXPathway:1"] = {"ncats_pathway_id": "IFXPathway:1"}
    pathway.edges_by_pathway["IFXPathway:1"] = [{
        "target_id": "IFXGene:1", "target_label": "GENE1"
    }]
    disease = SimpleNamespace(
        nodes=[{"ncats_disease_id": "IFXDisease:1", "consolidated_disease_name": "Disease"}],
        associations_by_ncats_id={"IFXDisease:1": [{"ncats_gene_id": "IFXGene:1"}]},
    )
    drug = SimpleNamespace(
        nodes=[{"drug_id": "IFXDrug:1", "standard_name": "Drug"}],
        edges_by_drug={"IFXDrug:1": [{
            "target_id": "IFXProtein:1", "target_label": "Gene one protein",
            "target_gene_id": "IFXGene:1", "target_symbol": "GENE1",
            "target_category": "biolink:Protein",
            "relation_kind": "drug_target",
        }]},
    )

    result = build_cross_entity_summary(
        pathway, "IFXPathway:1", disease_data=disease, drug_data=drug
    )

    assert result["disease_links"][0]["disease_id"] == "IFXDisease:1"
    assert result["drug_links"][0]["drug_id"] == "IFXDrug:1"
