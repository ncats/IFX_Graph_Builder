"""Disease Ontology's OBO adapter preserves the graph and table contracts."""

from src.input_adapters.disease_ontology.do_adapter import (
    DODiseaseAdapter,
    DODiseaseParentEdgeAdapter,
)
from src.input_adapters.pharos_source_tcrd.ontology_tables import (
    DOTableAdapter,
    DOTableParentEdgeAdapter,
)


class Source:
    def __init__(self, path):
        self.path = path

    def file(self, name=None):
        assert name == "doid.obo"
        return self.path

    def version_info(self):
        return None


def test_obo_disease_and_pharos_table_adapters(tmp_path):
    path = tmp_path / "doid.obo"
    path.write_text(
        "format-version: 1.2\n"
        "data-version: releases/2026-09-30/doid.obo\n\n"
        "[Term]\n"
        "id: DOID:1\n"
        "name: Parent disease\n"
        "def: \"A parent.\" []\n\n"
        "[Term]\n"
        "id: DOID:2\n"
        "name: Child disease\n"
        "def: \"A child with \\\"quotes\\\".\" []\n"
        "synonym: \"Other alias\" EXACT []\n"
        "synonym: \"Child alias\" EXACT []\n"
        "synonym: \"A related alias\" RELATED []\n"
        "comment: Reviewed term\n"
        "is_a: DOID:1 ! Parent disease\n"
        "relationship: part_of DOID:1 ! Non-parent link\n\n"
        "[Term]\n"
        "id: DOID:3\n"
        "name: Obsolete disease\n"
        "is_obsolete: true\n",
        encoding="utf-8",
    )
    source = Source(path)

    diseases = next(DODiseaseAdapter(source).get_all())
    assert {term.id for term in diseases} == {"DOID:1", "DOID:2"}
    child = next(term for term in diseases if term.id == "DOID:2")
    assert child.type == "CLASS"
    assert child.do_description == 'A child with "quotes".'
    assert child.synonyms == ["Child alias", "Other alias", "A related alias"]
    assert child.comments == ["Reviewed term"]

    edges = next(DODiseaseParentEdgeAdapter(source).get_all())
    assert [(edge.start_node.id, edge.end_node.id) for edge in edges] == [
        ("DOID:2", "DOID:1")
    ]
    assert [(term.id, term.do_description) for term in next(DOTableAdapter(source).get_all())] == [
        ("DOID:1", "A parent."),
        ("DOID:2", 'A child with "quotes".'),
    ]
    table_edges = next(DOTableParentEdgeAdapter(source).get_all())
    assert [(edge.start_node.id, edge.end_node.id) for edge in table_edges] == [
        ("DOID:2", "DOID:1")
    ]
