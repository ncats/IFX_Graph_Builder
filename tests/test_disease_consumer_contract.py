import csv
import hashlib
from collections import namedtuple
from types import SimpleNamespace

import pytest

from src.id_resolvers.disease_consumer_resolver import DiseaseConsumerResolver
from src.input_adapters.disease_harmonizer.consumer import DiseaseConsumerPayloadAdapter
from src.interfaces.id_resolver import IdResolver
from src.models.disease import Disease
from src.models.datasource_version_info import DatasourceVersionInfo
from src.output_adapters.arango_output_adapter import ArangoOutputAdapter
from src.output_adapters.sql_converters.tcrd import TCRDOutputConverter
from src.output_adapters.mysql_output_adapter import TCRDOutputAdapter


COLUMNS = [
    "ncats_disease_id", "primary_xref", "MONDO_xref", "GARD_xref",
    "DOID_xref", "MEDGEN_xref", "MESH_xref", "non_exact_xrefs", "standard_id", "nn_curie",
    "harmonized_name", "is_rare", "rare_basis", "review_note",
]


class _Release:
    snapshot_id = "ifx_harmonizers:diseases:2.3.0"

    def __init__(self, path, edges_path, labels_path):
        self.path = path
        self.edges_path = edges_path
        self.labels_path = labels_path
        self.manifest = {"files": [
            {"path": "disease_ids.tsv", "sha256": hashlib.sha256(path.read_bytes()).hexdigest()},
            {"path": "app_graph/disease_xref_edges.tsv", "sha256": hashlib.sha256(edges_path.read_bytes()).hexdigest()},
            {"path": "app_graph/xref_labels.tsv", "sha256": hashlib.sha256(labels_path.read_bytes()).hexdigest()},
        ]}

    def file(self, name):
        return {
            "disease_ids.tsv": self.path,
            "app_graph/disease_xref_edges.tsv": self.edges_path,
            "app_graph/xref_labels.tsv": self.labels_path,
        }[name]

    def version_info(self):
        return DatasourceVersionInfo(version="2.3.0", version_date=None, download_date=None)


def _release(tmp_path, rows, edges=(), labels=None):
    path = tmp_path / "disease_ids.tsv"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=COLUMNS, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    edges_path = tmp_path / "disease_xref_edges.tsv"
    with edges_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(
            handle, fieldnames=["ncats_disease_id", "xref_id", "match_type", "xref_is_obsolete"], delimiter="\t"
        )
        writer.writeheader()
        writer.writerows({**edge, "xref_is_obsolete": edge.get("xref_is_obsolete", "false")} for edge in edges)
    labels_path = tmp_path / "xref_labels.tsv"
    with labels_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["xref_id", "is_obsolete"], delimiter="\t")
        writer.writeheader()
        writer.writerows(labels if labels is not None else [
            {"xref_id": "MONDO:0000001", "is_obsolete": "false"}
        ])
    return _Release(path, edges_path, labels_path)


def _row(**changes):
    row = {
        "ncats_disease_id": "IFXDisease:ONE",
        "primary_xref": "MONDO:0000001",
        "MONDO_xref": "MONDO:0000001",
        "GARD_xref": "GARD:0000001",
        "DOID_xref": "DOID:1",
        "MEDGEN_xref": "MEDGEN:123",
        "MESH_xref": "",
        "non_exact_xrefs": "MESH:D1=related",
        "standard_id": "MONDO:9999999",
        "nn_curie": "MONDO:8888888",
        "harmonized_name": "Example disease",
        "is_rare": "true",
        "rare_basis": "GARD:0000001",
        "review_note": "keep this field",
    }
    row.update(changes)
    return row


def test_resolver_uses_only_exact_aliases_and_allows_source_keyed_misses(tmp_path):
    release = _release(tmp_path, [_row()])
    resolver = DiseaseConsumerResolver(
        release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"],
        no_match_behavior="Allow", multi_match_behavior="Error",
    )
    matches = resolver.resolve_internal([
        Disease(id="GARD:0000001"), Disease(id="MESH:D1"),
        Disease(id="MONDO:9999999"), Disease(id="MONDO:8888888"),
    ])
    assert matches["GARD:0000001"][0].match == "IFXDisease:ONE"
    assert set(matches) == {"GARD:0000001"}
    resolved = resolver.resolve_nodes([Disease(id="GARD:0000001"), Disease(id="MESH:D1")])
    nodes = resolver.parse_flat_node_list_from_map(resolved)
    assert {node.id for node in nodes} == {"IFXDisease:ONE", "MESH:D1"}
    assert next(node for node in nodes if node.id == "MESH:D1").resolver_miss is True
    assert next(node for node in nodes if node.id == "MESH:D1").xref_match_types == {}
    canonical = resolved[IdResolver.MatchKeys.matched]["GARD:0000001"]
    assert "GARD:0000001" in {xref.id_str() for xref in canonical.xref}
    assert canonical.xref_match_types == {
        "exact": sorted(xref.id_str() for xref in canonical.xref)
    }
    assert "medgen:123" in canonical.xref_match_types["exact"]
    assert "MESH:D1" not in canonical.xref_match_types["exact"]
    assert "MONDO:9999999" not in canonical.xref_match_types["exact"]
    assert ArangoOutputAdapter.clean_dict(
        object.__new__(ArangoOutputAdapter), canonical, convert_dates=True
    )["xref_match_types"] == canonical.xref_match_types


def test_duplicate_exact_alias_and_corrupt_snapshot_fail(tmp_path):
    release = _release(tmp_path, [
        _row(),
        _row(ncats_disease_id="IFXDisease:TWO", primary_xref="DOID:2",
             MONDO_xref="", GARD_xref="GARD:0000001", DOID_xref="DOID:2"),
    ])
    with pytest.raises(ValueError, match="belongs to both"):
        DiseaseConsumerResolver(release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"])
    release.manifest["files"][0]["sha256"] = "0" * 64
    with pytest.raises(ValueError, match="checksum mismatch"):
        DiseaseConsumerPayloadAdapter(release)


def test_resolver_indexes_exact_edges_and_keeps_ambiguity_explicit(tmp_path):
    release = _release(tmp_path, [
        _row(MESH_xref="MESH:D031901"),
        _row(ncats_disease_id="IFXDisease:TWO", primary_xref="MONDO:0000002",
             MONDO_xref="MONDO:0000002", GARD_xref="", DOID_xref="", MEDGEN_xref="",
             MESH_xref=""),
    ], edges=[
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "MESH:D031901", "match_type": "broad"},
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "MESH:D2", "match_type": "exact"},
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "NCIT:C3396", "match_type": "exact"},
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "NCIt:C3396", "match_type": "broad"},
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "ICD11:2B50", "match_type": "related"},
        {"ncats_disease_id": "IFXDisease:TWO", "xref_id": "MESH:D2", "match_type": "exact"},
    ])
    resolver = DiseaseConsumerResolver(
        release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"],
        no_match_behavior="Allow", multi_match_behavior="Error",
    )
    matches = resolver.resolve_internal([Disease(id="MESH:D2"), Disease(id="MESH:D031901")])
    assert [match.match for match in matches["MESH:D2"]] == [
        "IFXDisease:ONE", "IFXDisease:TWO"
    ]
    assert [match.match for match in matches["MESH:D031901"]] == ["IFXDisease:ONE"]
    ambiguous = resolver.resolve_nodes([Disease(id="MESH:D2")])
    with pytest.raises(ValueError, match="resolves to multiple targets"):
        resolver.parse_flat_node_list_from_map(ambiguous)
    resolved = resolver.resolve_nodes([Disease(id="IFXDisease:ONE"), Disease(id="IFXDisease:TWO")])
    one = resolved[IdResolver.MatchKeys.matched]["IFXDisease:ONE"]
    two = resolved[IdResolver.MatchKeys.matched]["IFXDisease:TWO"]
    assert "MESH:D031901" in one.xref_match_types["broad"]
    assert "MESH:D031901" not in one.xref_match_types["exact"]
    assert "MESH:D2" in one.xref_match_types["exact"]
    assert "MESH:D2" in two.xref_match_types["exact"]
    assert "NCIT:C3396" in one.xref_match_types["exact"]
    assert "NCIt:C3396" in one.xref_match_types["broad"]
    assert "ICD11:2B50" in one.xref_match_types["related"]
    assert {xref.id_str() for xref in one.xref} >= {"NCIT:C3396", "NCIt:C3396"}
    converter = TCRDOutputConverter()
    converter.id_mapping["ncats_disease"] = {"IFXDisease:ONE": 1}
    graph_node = ArangoOutputAdapter.clean_dict(
        object.__new__(ArangoOutputAdapter), one, convert_dates=True
    )
    mysql_rows = converter.disease_entity_xref_converter(graph_node)
    assert {(row.source_id, row.match_type) for row in mysql_rows} >= {
        ("NCIT:C3396", "exact"), ("NCIt:C3396", "broad"),
        ("ICD11:2B50", "related"),
    }


def test_obsolete_labels_exclude_aliases_and_edges_from_resolution_and_xrefs(tmp_path):
    release = _release(tmp_path, [_row(MESH_xref="MESH:OLD")], edges=[
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "MESH:OLD", "match_type": "exact"},
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "MESH:NEW", "match_type": "exact"},
        {"ncats_disease_id": "IFXDisease:ONE", "xref_id": "MESH:EDGE_OLD",
         "match_type": "exact", "xref_is_obsolete": "true"},
    ], labels=[
        {"xref_id": "MESH:OLD", "is_obsolete": "true"},
        {"xref_id": "MESH:NEW", "is_obsolete": "false"},
    ])
    resolver = DiseaseConsumerResolver(
        release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"],
    )
    assert resolver.resolve_internal([Disease(id="MESH:NEW")])["MESH:NEW"][0].match == "IFXDisease:ONE"
    assert resolver.resolve_internal([Disease(id="MESH:EDGE_OLD")]) == {}
    assert resolver.resolve_internal([Disease(id="MESH:OLD")]) == {}
    obsolete = resolver.resolve_nodes([Disease(id="MESH:OLD")])
    assert not obsolete[IdResolver.MatchKeys.matched]
    unmatched = resolver.parse_flat_node_list_from_map(obsolete)
    assert len(unmatched) == 1
    assert unmatched[0].id == "MESH:OLD"
    assert unmatched[0].resolver_miss is True
    assert unmatched[0].xref_match_types == {}
    node = resolver.resolve_nodes([Disease(id="IFXDisease:ONE")])[IdResolver.MatchKeys.matched]["IFXDisease:ONE"]
    all_xrefs = {xref.id_str() for xref in node.xref}
    assert "MESH:NEW" in all_xrefs
    assert "MESH:OLD" not in all_xrefs
    assert "MESH:EDGE_OLD" not in all_xrefs



def test_obsolete_gard_source_stays_unmatched_through_mysql_conversion(tmp_path):
    old_id = "GARD:0000963"
    release = _release(tmp_path, [_row(GARD_xref=f"GARD:0000001|{old_id}")], labels=[
        {"xref_id": old_id, "is_obsolete": "true"},
        {"xref_id": "MONDO:0000001", "is_obsolete": "false"},
    ])
    resolver = DiseaseConsumerResolver(
        release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"],
        no_match_behavior="Allow", multi_match_behavior="Error",
    )
    resolved = resolver.resolve_nodes([
        Disease(id=old_id, name="RDAS disease", rare_disease=True, provenance="RDAS")
    ])
    assert not resolved[IdResolver.MatchKeys.matched]
    node = resolver.parse_flat_node_list_from_map(resolved)[0]
    assert node.id == old_id
    assert node.resolver_miss is True
    assert node.xref_match_types == {}
    graph_node = ArangoOutputAdapter.clean_dict(
        object.__new__(ArangoOutputAdapter), node, convert_dates=True
    )
    converter = TCRDOutputConverter()
    mysql_node = converter.ncats_disease_converter(graph_node)
    assert mysql_node.ifx_id is None
    assert mysql_node.gard_rare is True
    assert converter.disease_entity_xref_converter(graph_node) == []
    canonical = resolver.resolve_nodes([Disease(id="IFXDisease:ONE")])[IdResolver.MatchKeys.matched]["IFXDisease:ONE"]
    assert old_id not in {xref.id_str() for xref in canonical.xref}


def test_obsolete_primary_and_duplicate_labels_fail(tmp_path):
    release = _release(tmp_path, [_row()], labels=[
        {"xref_id": "MONDO:0000001", "is_obsolete": "true"}
    ])
    with pytest.raises(ValueError, match="obsolete primary_xref"):
        DiseaseConsumerResolver(release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"])

    release = _release(tmp_path, [_row()], labels=[
        {"xref_id": "MESH:D1", "is_obsolete": "true"},
        {"xref_id": "MESH:D1", "is_obsolete": "false"},
    ])
    with pytest.raises(ValueError, match="duplicates MESH:D1"):
        DiseaseConsumerResolver(release, cache_path=str(tmp_path / "disease.sqlite"), types=["Disease"])


def test_label_checksum_change_rebuilds_resolver_cache(tmp_path):
    release = _release(tmp_path, [_row(MESH_xref="MESH:OLD")], labels=[
        {"xref_id": "MESH:OLD", "is_obsolete": "false"}
    ])
    cache_path = str(tmp_path / "disease.sqlite")
    original = DiseaseConsumerResolver(release, cache_path=cache_path, types=["Disease"])
    assert original.resolve_internal([Disease(id="MESH:OLD")])
    original.connection.close()
    original.connection = None

    release.labels_path.write_text(
        "xref_id\tis_obsolete\nMESH:OLD\ttrue\n", encoding="utf-8"
    )
    release.manifest["files"][2]["sha256"] = hashlib.sha256(
        release.labels_path.read_bytes()
    ).hexdigest()
    updated = DiseaseConsumerResolver(release, cache_path=cache_path, types=["Disease"])
    assert updated.connection.execute(
        "SELECT COUNT(*) FROM matches WHERE match = 'MESH:OLD'"
    ).fetchone()[0] == 0


def test_payload_keeps_every_column_and_round_trips_through_disease_model(tmp_path):
    row = _row()
    adapter = DiseaseConsumerPayloadAdapter(_release(tmp_path, [row]))
    node = next(adapter.get_all())[0]
    assert node.id == "IFXDisease:ONE"
    assert node.name == "Example disease"
    assert node.rare_disease is True
    assert node.disease_harmonizer_payload == row
    assert Disease.from_dict(node.__dict__).disease_harmonizer_payload == row


def test_pharos_uses_exact_mondo_alias_and_release_rarity_for_ifx_disease(tmp_path):
    converter = TCRDOutputConverter()
    converter._known_mondo_ids.add("MONDO:0000001")
    payload = _row()
    node = converter.ncats_disease_converter({
        "id": "IFXDisease:ONE", "name": "Example disease", "rare_disease": False,
        "disease_harmonizer_payload": payload, "provenance": "DISEASES\t2.3.0\tNone\tNone",
    })
    assert node.mondoid == "MONDO:0000001"
    assert node.gard_rare is True
    converter.id_mapping["protein"] = {"IFXProtein:ONE": 1}
    association = converter.disease_converter({
        "start_id": "IFXProtein:ONE", "end_id": "IFXDisease:ONE",
        "details": [{"source": "example", "source_id": "DOID:1"}],
        "provenance": "example\t1\tNone\tNone",
    })
    assert association[0].mondoid == "MONDO:0000001"

def test_pharos_keeps_disease_when_exact_mondo_alias_is_absent_from_table():
    converter = TCRDOutputConverter()
    # MONDO:0100332 occurs in DISEASES 2.3.0 but not the pinned MONDO table.
    node = converter.ncats_disease_converter({
        "id": "IFXDisease:ONE", "name": "Example disease",
        "disease_harmonizer_payload": _row(
            primary_xref="MONDO:0100332", MONDO_xref="MONDO:0100332"
        ),
        "provenance": "DISEASES\t2.3.0",
    })
    assert node.mondoid is None
    assert node.gard_rare is True


def test_pharos_resume_restores_ifx_disease_mondoid_from_persisted_row(monkeypatch):
    persisted = namedtuple("PersistedDisease", "id name uniprot_description do_description mondo_description mondoid")(
        42, "Example disease", None, None, None, "MONDO:0000001"
    )
    duplicate = persisted._replace(id=43, mondoid=None)

    class Query:
        def __init__(self, rows):
            self.rows = rows

        def filter(self, *_args):
            return self

        def all(self):
            return self.rows

        def __iter__(self):
            return iter(self.rows)

    class Session:
        def query(self, *columns):
            names = tuple(column.name for column in columns)
            if names == ("id", "mondoid"):
                return Query([(persisted.id, persisted.mondoid), (duplicate.id, duplicate.mondoid)])
            if names == ("id", "ifx_id"):
                return Query([(42, "IFXDisease:ONE"), (43, "IFXDisease:TWO")])
            if names == ("id", "name", "uniprot_description", "do_description", "mondo_description", "mondoid"):
                return Query([persisted, duplicate])
            if names in {("id", "name"), ("id", "pcid")}:
                return Query([])
            raise AssertionError(names)

    class Graph:
        def has_collection(self, name):
            return name == "Disease"

    class Arango:
        def __init__(self, *_args):
            pass

        def get_db(self):
            return Graph()

        def runQuery(self, _query):
            return [
                {"id": "IFXDisease:ONE", "name": "Example disease"},
                {"id": "IFXDisease:TWO", "name": "Example disease"},
            ]

    monkeypatch.setattr("src.output_adapters.mysql_output_adapter.ArangoAdapter", Arango)
    adapter = TCRDOutputAdapter.__new__(TCRDOutputAdapter)
    adapter.output_converter = TCRDOutputConverter()
    adapter.output_converter.id_mapping = {}
    adapter.source_graph_credentials = SimpleNamespace()
    adapter.source_graph_database = "working_diseases_230"
    adapter._preload_graph_backed_id_mappings(Session())
    converter = adapter.output_converter
    assert converter.id_mapping["ncats_disease"]["IFXDisease:ONE"] == 42
    assert converter.id_mapping["ncats_disease"]["IFXDisease:TWO"] == 43
    assert converter._disease_mondoid_by_id["IFXDisease:ONE"] == "MONDO:0000001"
    assert converter._disease_mondoid_by_id["IFXDisease:TWO"] is None
    assert converter.resolve_id("ncats_disease", "IFXDisease:THREE") == 44
    converter.id_mapping["protein"] = {"IFXProtein:ONE": 1}
    association = converter.disease_converter({
        "start_id": "IFXProtein:ONE", "end_id": "IFXDisease:ONE",
        "details": [{"source": "example", "source_id": "DOID:1"}],
        "provenance": "example\t1\tNone\tNone",
    })
    assert association[0].mondoid == "MONDO:0000001"

    class MissingIdentitySession(Session):
        def query(self, *columns):
            if tuple(column.name for column in columns) == ("id", "ifx_id"):
                return Query([])
            return super().query(*columns)

    missing_identity_adapter = TCRDOutputAdapter.__new__(TCRDOutputAdapter)
    missing_identity_adapter.output_converter = TCRDOutputConverter()
    missing_identity_adapter.output_converter.id_mapping = {}
    missing_identity_adapter.source_graph_credentials = SimpleNamespace()
    missing_identity_adapter.source_graph_database = "working_diseases_230"
    with pytest.raises(RuntimeError, match="verified exact ID backfill"):
        missing_identity_adapter._preload_graph_backed_id_mappings(MissingIdentitySession())


def test_pharos_requires_unique_ifx_id_schema_before_conversion(monkeypatch):
    class Inspector:
        def __init__(self, columns, unique, nullable=True):
            self.columns = columns
            self.unique = unique
            self.nullable = nullable

        def get_columns(self, _table):
            return [{"name": name, "nullable": self.nullable} for name in self.columns]

        def get_unique_constraints(self, _table):
            return [{"column_names": ["ifx_id"]}] if self.unique else []

        def get_indexes(self, _table):
            return []

    adapter = TCRDOutputAdapter.__new__(TCRDOutputAdapter)
    adapter.get_engine = lambda: object()
    for inspector in (
        Inspector({"id"}, False),
        Inspector({"id", "ifx_id"}, False),
        Inspector({"id", "ifx_id"}, True, nullable=False),
    ):
        monkeypatch.setattr("src.output_adapters.mysql_output_adapter.sqlalchemy_inspect", lambda _engine: inspector)
        with pytest.raises(RuntimeError, match="nullable, unique ifx_id"):
            adapter._validate_disease_identity_schema()
    monkeypatch.setattr(
        "src.output_adapters.mysql_output_adapter.sqlalchemy_inspect",
        lambda _engine: Inspector({"id", "ifx_id"}, True),
    )
    adapter._validate_disease_identity_schema()
