"""Keep the DISEASES projection distinct from the source ontology edges."""

import csv
import hashlib
import json

import pytest

from src.input_adapters.disease_harmonizer.hierarchy import (
    DiseaseConsumerHierarchyAdapter, HIERARCHY_MEMBER,
)
from src.input_adapters.mondo.mondo_adapter import MondoDiseaseParentEdgeAdapter
from src.models.disease import DiseaseParentEdge, MondoDiseaseParentEdge
from src.models.datasource_version_info import DatasourceVersionInfo


COLUMNS = [
    "parent_primary_xref", "parent_ncats_disease_id", "parent_label",
    "child_primary_xref", "child_ncats_disease_id", "child_label",
    "relationship_type", "predicate", "source",
]


class Source:
    snapshot_id = "ifx_harmonizers:diseases:2.3.0"

    def __init__(self, path):
        self.path = path
        self.manifest = {"files": [{
            "path": HIERARCHY_MEMBER,
            "sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
        }]}

    def file(self, name):
        assert name == HIERARCHY_MEMBER
        return self.path

    def version_info(self):
        return DatasourceVersionInfo(version="2.3.0", version_date=None, download_date=None)


def _row(**changes):
    row = {
        "parent_primary_xref": "MONDO:1",
        "parent_ncats_disease_id": "IFXDisease:PARENT",
        "parent_label": "Parent",
        "child_primary_xref": "MONDO:2",
        "child_ncats_disease_id": "IFXDisease:CHILD",
        "child_label": "Child",
        "relationship_type": "parent_child",
        "predicate": "biolink:subclass_of",
        "source": "MONDO",
    }
    row.update(changes)
    return row


def _source(tmp_path, rows):
    path = tmp_path / "disease_hierarchy_edges.tsv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return Source(path)


def test_ifx_hierarchy_emits_disease_parent_edges_with_original_row(tmp_path):
    row = _row()
    edge = next(DiseaseConsumerHierarchyAdapter(_source(tmp_path, [row])).get_all())[0]
    assert type(edge) is DiseaseParentEdge
    assert (edge.start_node.id, edge.end_node.id) == ("IFXDisease:CHILD", "IFXDisease:PARENT")
    assert edge.source == "MONDO"
    assert edge.hierarchy_payload == row
    assert DiseaseParentEdge.from_dict({
        "start_id": edge.start_node.id, "end_id": edge.end_node.id,
        "source": edge.source, "hierarchy_payload": edge.hierarchy_payload,
    }).hierarchy_payload == row


def test_ifx_hierarchy_rejects_corrupt_or_malformed_release(tmp_path):
    source = _source(tmp_path, [_row()])
    source.manifest["files"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checksum mismatch"):
        DiseaseConsumerHierarchyAdapter(source)
    source = _source(tmp_path, [_row(source="DO")])
    with pytest.raises(ValueError, match="unexpected hierarchy semantics"):
        list(DiseaseConsumerHierarchyAdapter(source).get_all())


def test_raw_mondo_uses_distinct_parent_edge_class(tmp_path):
    path = tmp_path / "mondo.json"
    path.write_text(json.dumps({"graphs": [{
        "nodes": [
            {"id": "http://purl.obolibrary.org/obo/MONDO_1", "type": "CLASS"},
            {"id": "http://purl.obolibrary.org/obo/MONDO_2", "type": "CLASS"},
        ],
        "edges": [{
            "sub": "http://purl.obolibrary.org/obo/MONDO_2",
            "obj": "http://purl.obolibrary.org/obo/MONDO_1",
            "pred": "is_a",
        }],
    }]}), encoding="utf-8")

    class MondoSource:
        def file(self, _name=None):
            return path

        def version_info(self):
            return DatasourceVersionInfo(version="test", version_date=None, download_date=None)

    edge = next(MondoDiseaseParentEdgeAdapter(MondoSource()).get_all())[0]
    assert type(edge) is MondoDiseaseParentEdge
    assert (edge.start_node.id, edge.end_node.id, edge.source) == ("MONDO:2", "MONDO:1", "MONDO")
