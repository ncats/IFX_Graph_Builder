import json
from datetime import date
from pathlib import Path
from types import SimpleNamespace

import pytest

from src.core.curations import (
    CHEBI_RECORD_PROPERTIES,
    METABOLITE_EQUIVALENCE_EDGES,
    METABOLITE_RECORD_PROPERTIES,
    record_property_decisions_from_operation,
)
from src.models.protein import Protein
from src.models.test_models import TestEdge, TestNode
from src.interfaces.resolver_metadata import resolver_fingerprints_by_type
from src.output_adapters.arango_output_adapter import ArangoOutputAdapter
from arango.exceptions import DocumentUpdateError
from src.shared.record_merger import FieldConflictBehavior
from src.models.registry_dataset import RegistryDataset, RegistryDatasetKind


class FakeDocumentUpdateError(DocumentUpdateError):
    def __init__(self):
        Exception.__init__(self, "fake 413")


class FakeCollection:
    def __init__(self):
        self.insert_calls = []
        self.update_calls = []

    def insert_many(self, docs, overwrite=False):
        self.insert_calls.append({
            "docs": docs,
            "overwrite": overwrite,
        })
        return []

    def update_many(self, docs, merge=True, keep_none=False, check_rev=False):
        self.update_calls.append({
            "docs": docs,
            "merge": merge,
            "keep_none": keep_none,
            "check_rev": check_rev,
        })
        return []


class FailingUpdateCollection(FakeCollection):
    def update_many(self, docs, merge=True, keep_none=False, check_rev=False):
        self.update_calls.append({
            "docs": docs,
            "merge": merge,
            "keep_none": keep_none,
            "check_rev": check_rev,
        })
        if len(docs) > 1:
            raise FakeDocumentUpdateError()
        return []


class FakeGraph:
    def __init__(self, edge_collection=None):
        self.edge_collection = edge_collection or FakeCollection()

    def has_edge_collection(self, label):
        return False

    def create_edge_definition(self, label, from_vertex_collections, to_vertex_collections):
        return self.edge_collection


class FakeCursor(list):
    pass


class FakeAql:
    def execute(self, query):
        if "RETURN COUNT" in query:
            return FakeCursor([4])
        if "COLLECT value = item" in query:
            return FakeCursor([
                {"value": "hmdb\t5.0\t2026-01-01\t2026-01-02", "count": 2},
                {"value": "", "count": 1},
            ])
        if "COLLECT combo = key" in query:
            return FakeCursor([
                {"combination": "hmdb\t5.0\t2026-01-01\t2026-01-02", "count": 2},
                {"combination": "", "count": 1},
            ])
        return FakeCursor([])


class FakeMetadataDb:
    aql = FakeAql()

    def collections(self):
        return [
            {"system": False, "name": "TestCollection"},
            {"system": False, "name": "MetaboliteHarmonizationClique"},
        ]


class FakeMalformedAql:
    def execute(self, query):
        if "RETURN COUNT" in query:
            return FakeCursor([1])
        if "COLLECT value = item" in query:
            return FakeCursor([{"value": "bad-source-fragment", "count": 1}])
        return FakeCursor([])


class FakeMalformedMetadataDb:
    aql = FakeMalformedAql()

    def collections(self):
        return [{"system": False, "name": "TestCollection"}]


def make_protein(protein_id: str, name: str, entity_resolution: str, provenance: str) -> Protein:
    protein = Protein(id=protein_id, name=name)
    protein.entity_resolution = entity_resolution
    protein.provenance = provenance
    return protein


def build_adapter(existing_nodes, collection):
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter._collection_schemas = {}
    adapter._graph_views = []
    adapter._graph_view_source_yaml = None
    adapter.object_storage = None
    adapter.database_name = "test_db"
    adapter._handle_dataset_nodes = lambda objects: None
    adapter._handle_pounce_workbook_nodes = lambda objects: None
    adapter.get_db = lambda: object()
    adapter.get_graph = lambda: FakeGraph()
    adapter.create_indexes = lambda obj_cls, coll: None
    adapter.get_existing_nodes = lambda db, label, obj_list, skip_merge=False: (collection, existing_nodes)
    return adapter


def test_get_existing_nodes_uses_independent_bounded_queries():
    collection = FakeCollection()
    calls = []

    def execute(query, *, bind_vars, batch_size, max_runtime):
        calls.append((bind_vars, batch_size))
        assert len(bind_vars["keys"]) < batch_size
        return [
            {"_key": key, "id": key}
            for key in bind_vars["keys"]
            if key.endswith("0")
        ]

    db = SimpleNamespace(
        has_collection=lambda label: True,
        collection=lambda label: collection,
        aql=SimpleNamespace(execute=execute),
    )
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    objects = [{"id": f"P{i}", "name": "protein"} for i in range(2001)]

    returned_collection, existing = adapter.get_existing_nodes(db, "Protein", objects)

    assert returned_collection is collection
    assert [len(call[0]["keys"]) for call in calls] == [1000, 1000, 1]
    assert len(existing) == 201
    assert all(call[0]["fields"] == calls[0][0]["fields"] for call in calls)


def test_arango_output_adapter_etl_metadata_includes_readable_resolver_metadata():
    data_source = RegistryDataset(
        kind=RegistryDatasetKind.SOURCE,
        source="target_graph",
        dataset="disease_ids",
        version="deps-test",
        version_date=None,
        download_date=date(2026, 6, 12),
        snapshot_id="target_graph:disease_ids:deps-test",
        manifest_uri="s3://ifx-registry/sources/target_graph/disease_ids/deps-test/manifest.yaml",
        manifest={"kind": "source_snapshot", "files": []},
        local_dir=Path("/tmp/ifx-registry-cache/target_graph/disease_ids/deps-test"),
    )
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    resolver_metadata = resolver_fingerprints_by_type([{
        "label": "disease_ids",
        "import": "./src/id_resolvers/disease_resolver.py",
        "class": "DiseaseIdResolver",
        "kwargs": {
            "data_source": data_source,
            "types": ["Disease"],
            "multi_match_behavior": "All",
        },
    }])
    adapter.set_resolver_metadata(
        resolver_fingerprints_by_type=resolver_metadata,
        source_yaml="./src/use_cases/pharos/pharos.yaml",
    )

    metadata = adapter.get_etl_metadata()

    disease_metadata = metadata["resolver_metadata"]["by_type"]["Disease"]
    assert metadata["resolver_metadata"]["source_yaml"] == "./src/use_cases/pharos/pharos.yaml"
    assert disease_metadata["class"] == "DiseaseIdResolver"
    assert disease_metadata["dataset_inputs"][0]["snapshot_id"] == "target_graph:disease_ids:deps-test"
    assert disease_metadata["kwargs"]["data_source"]["snapshot_id"] == "target_graph:disease_ids:deps-test"
    assert "local_dir" not in disease_metadata["kwargs"]["data_source"]
    assert disease_metadata["fingerprint"]


def test_arango_output_adapter_get_metadata_skips_empty_source_combinations_and_qa_artifacts():
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.metadata_store_label = "metadata_store"
    adapter.get_db = lambda: FakeMetadataDb()

    metadata = adapter.get_metadata()

    collection = metadata.collections[0]
    assert collection.name == "TestCollection"
    assert collection.total_count == 4
    assert [source.name for source in collection.sources] == ["hmdb"]
    assert collection.marginal_source_counts == {"hmdb": 2}
    assert collection.joint_source_counts == {"hmdb": 2}
    assert [collection.name for collection in metadata.collections] == ["TestCollection"]


def test_arango_output_adapter_get_metadata_keeps_source_name_without_optional_metadata():
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.metadata_store_label = "metadata_store"
    adapter.get_db = lambda: FakeMalformedMetadataDb()

    metadata = adapter.get_metadata()

    source = metadata.collections[0].sources[0]
    assert source.name == "bad-source-fragment"
    assert source.version is None


def test_arango_output_adapter_resolver_metadata_serializes_registry_dataset(tmp_path):
    data_source = RegistryDataset(
        kind=RegistryDatasetKind.SOURCE,
        source="cure",
        dataset="cure_id_labels",
        version="deps-test",
        version_date=None,
        download_date=date(2026, 6, 12),
        snapshot_id="cure:cure_id_labels:deps-test",
        manifest_uri="s3://ifx-registry/sources/cure/cure_id_labels/deps-test/manifest.yaml",
        manifest={"kind": "source_snapshot", "files": []},
        local_dir=tmp_path,
    )
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    resolver_metadata = resolver_fingerprints_by_type([{
        "label": "cure_id_labels",
        "import": "./src/id_resolvers/cure_id_label_resolver.py",
        "class": "CureIdLabelResolver",
        "kwargs": {
            "data_source": data_source,
            "types": ["Gene"],
        },
    }])
    adapter.set_resolver_metadata(
        resolver_fingerprints_by_type=resolver_metadata,
        source_yaml="./src/use_cases/cure/cure_rasopathies.yaml",
    )

    metadata = adapter.get_etl_metadata()

    json.dumps(metadata)
    snapshot_metadata = metadata["resolver_metadata"]["by_type"]["Gene"]["dataset_inputs"][0]
    assert snapshot_metadata["snapshot_id"] == "cure:cure_id_labels:deps-test"
    assert snapshot_metadata["kind"] == "source_snapshot"


def test_arango_output_adapter_etl_metadata_includes_registry_datasets():
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.set_registry_dataset_metadata([{
        "source": "cure",
        "dataset": "case_reports",
        "version": "reports_20260612T182139Z",
        "snapshot_id": "cure:case_reports:reports_20260612T182139Z",
        "usages": ["adapter:CUREAdapter"],
    }])

    metadata = adapter.get_etl_metadata()

    assert metadata["registry_datasets"] == [{
        "source": "cure",
        "dataset": "case_reports",
        "version": "reports_20260612T182139Z",
        "snapshot_id": "cure:case_reports:reports_20260612T182139Z",
        "usages": ["adapter:CUREAdapter"],
    }]


def test_arango_output_adapter_merges_existing_resolver_metadata_for_post_processing():
    existing = resolver_fingerprints_by_type([{
        "label": "disease_ids",
        "import": "./src/id_resolvers/disease_resolver.py",
        "class": "DiseaseIdResolver",
        "kwargs": {
            "file_path": "./input_files/manual/target_graph/disease_ids.tsv",
            "multi_match_behavior": "All",
            "types": ["Disease"],
        },
    }])
    current = resolver_fingerprints_by_type([{
        "label": "tcrd_targets",
        "import": "./src/id_resolvers/target_graph_resolver.py",
        "class": "TCRDTargetResolver",
        "kwargs": {
            "canonical_type": "Protein",
            "collapse_to_canonical": True,
            "types": ["Protein", "Gene", "Transcript"],
        },
    }])

    merged = ArangoOutputAdapter._merge_resolver_metadata(
        {
            "source_yaml": "./src/use_cases/pharos/pharos.yaml",
            "by_type": existing,
        },
        {
            "source_yaml": "./src/use_cases/pharos/pharos_aql_post.yaml",
            "by_type": current,
        },
    )

    assert sorted(merged["by_type"]) == ["Disease", "Gene", "Protein", "Transcript"]
    assert merged["source_yamls"] == [
        "./src/use_cases/pharos/pharos.yaml",
        "./src/use_cases/pharos/pharos_aql_post.yaml",
    ]
    assert merged["summary"]["Disease"]["class"] == "DiseaseIdResolver"
    assert merged["summary"]["Protein"]["class"] == "TCRDTargetResolver"


def test_arango_output_adapter_merges_existing_registry_datasets_for_post_processing():
    merged = ArangoOutputAdapter._merge_registry_datasets(
        [{
            "source": "cure",
            "dataset": "case_reports",
            "version": "reports_20260612T182139Z",
            "snapshot_id": "cure:case_reports:reports_20260612T182139Z",
            "usages": ["adapter:CUREAdapter"],
        }],
        [{
            "source": "cure",
            "dataset": "case_reports",
            "version": "reports_20260612T182139Z",
            "snapshot_id": "cure:case_reports:reports_20260612T182139Z",
            "manifest_uri": "s3://ifx-registry/sources/cure/case_reports/reports_20260612T182139Z/manifest.yaml",
            "usages": ["adapter:RasopathiesAdapter"],
        }, {
            "source": "cure",
            "dataset": "curated_concepts",
            "version": "2026-05-14",
            "snapshot_id": "cure:curated_concepts:2026-05-14",
            "usages": ["resolver:cure_id_labels"],
        }],
    )

    assert [dataset["snapshot_id"] for dataset in merged] == [
        "cure:case_reports:reports_20260612T182139Z",
        "cure:curated_concepts:2026-05-14",
    ]
    assert merged[0]["manifest_uri"] == "s3://ifx-registry/sources/cure/case_reports/reports_20260612T182139Z/manifest.yaml"
    assert merged[0]["usages"] == ["adapter:CUREAdapter", "adapter:RasopathiesAdapter"]


def test_arango_output_adapter_keeps_same_snapshot_id_in_separate_kinds():
    snapshot_id = "example:records:1"

    merged = ArangoOutputAdapter._merge_registry_datasets(
        [{"kind": "source_snapshot", "snapshot_id": snapshot_id}],
        [{"kind": "derived_snapshot", "snapshot_id": snapshot_id}],
    )

    assert [(entry["kind"], entry["snapshot_id"]) for entry in merged] == [
        ("derived_snapshot", snapshot_id),
        ("source_snapshot", snapshot_id),
    ]


def test_get_node_merge_fetch_fields_keeps_only_merge_relevant_fields():
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    obj_list = [{
        "id": "IFXProtein:P1",
        "name": "Alpha",
        "pm_score": [1.0],
        "entity_resolution": "resolver-1",
        "provenance": "source-1",
        "_internal": "ignored",
    }]

    fields = adapter.get_node_merge_fetch_fields(obj_list)

    assert fields == [
        "_key",
        "creation",
        "id",
        "name",
        "pm_score",
        "resolved_ids",
        "updates",
    ]


def test_document_handle_uses_same_safe_key_as_node_writes():
    assert (
        ArangoOutputAdapter.document_handle("MetaboliteIdentifier", "CAS:100-09-4")
        == "MetaboliteIdentifier/CAS:100_minus_09_minus_4"
    )


def test_store_writes_edge_handles_with_safe_document_keys():
    collection = FakeCollection()
    adapter = build_adapter([], collection)
    adapter.get_graph = lambda: FakeGraph(collection)
    edge = TestEdge(
        start_node=TestNode(id="CAS:100-09-4"),
        end_node=TestNode(id="PUBCHEM.COMPOUND:62698"),
        provenance="test-source",
    )
    edge.entity_resolution = "test-resolver"

    adapter.store([edge], single_source=True)

    edge_doc = collection.insert_calls[0]["docs"][0]
    assert edge_doc["_from"] == "TestNode/CAS:100_minus_09_minus_4"
    assert edge_doc["_to"] == "TestNode/PUBCHEM.COMPOUND:62698"
    assert edge_doc["start_id"] == "CAS:100-09-4"
    assert edge_doc["end_id"] == "PUBCHEM.COMPOUND:62698"


def test_store_uses_update_many_for_existing_nodes_and_insert_many_for_new_nodes():
    existing_nodes = [{
        "_key": "IFXProtein:P1",
        "id": "IFXProtein:P1",
        "name": "Old name",
        "creation": "creation-source",
        "updates": ["existing-update"],
        "resolved_ids": ["resolver-old"],
    }]
    collection = FakeCollection()
    adapter = build_adapter(existing_nodes, collection)

    existing = make_protein("IFXProtein:P1", "New name", "resolver-new", "source-new")
    new = make_protein("IFXProtein:P2", "Brand new", "resolver-brand-new", "source-brand-new")

    adapter.store([existing, new], single_source=False)

    assert len(collection.update_calls) == 1
    assert len(collection.insert_calls) == 1

    updated_docs = collection.update_calls[0]["docs"]
    inserted_docs = collection.insert_calls[0]["docs"]

    assert len(updated_docs) == 1
    assert updated_docs[0]["_key"] == "IFXProtein:P1"
    assert updated_docs[0]["creation"] == "creation-source"
    assert sorted(updated_docs[0]["resolved_ids"]) == ["resolver-new", "resolver-old"]
    assert "existing-update" in updated_docs[0]["updates"]
    assert any("name\tOld name\tNew name\tsource-new\tKeepFirst" == update for update in updated_docs[0]["updates"])
    assert updated_docs[0]["name"] == "Old name"

    assert len(inserted_docs) == 1
    assert inserted_docs[0]["_key"] == "IFXProtein:P2"
    assert inserted_docs[0]["creation"] == "source-brand-new"
    assert inserted_docs[0]["resolved_ids"] == ["resolver-brand-new"]
    assert collection.insert_calls[0]["overwrite"] is False
    assert collection.update_calls[0]["merge"] is True
    assert collection.update_calls[0]["keep_none"] is False
    assert collection.update_calls[0]["check_rev"] is False


def test_store_uses_requested_field_conflict_behavior_for_existing_nodes():
    existing_nodes = [{
        "_key": "IFXProtein:P1",
        "id": "IFXProtein:P1",
        "name": "Old name",
        "creation": "creation-source",
        "updates": [],
        "resolved_ids": ["resolver-old"],
    }]
    collection = FakeCollection()
    adapter = build_adapter(existing_nodes, collection)

    incoming = make_protein("IFXProtein:P1", "New name", "resolver-new", "source-new")

    adapter.store(
        [incoming],
        single_source=False,
        field_conflict_behavior=FieldConflictBehavior.KeepLast,
    )

    updated_docs = collection.update_calls[0]["docs"]

    assert updated_docs[0]["name"] == "New name"
    assert any("name\tOld name\tNew name\tsource-new\tKeepLast" == update for update in updated_docs[0]["updates"])


def test_update_many_with_backoff_splits_on_document_update_error():
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    collection = FailingUpdateCollection()
    records = [
        {"_key": "a", "id": "a"},
        {"_key": "b", "id": "b"},
        {"_key": "c", "id": "c"},
        {"_key": "d", "id": "d"},
    ]

    adapter.update_many_with_backoff(collection, records, label="Protein", kind="node")

    assert [len(call["docs"]) for call in collection.update_calls] == [4, 2, 1, 1, 2, 1, 1]


class CurationAql:
    def __init__(self, lookup_rows, edge_rows=None, fail_on_update_number=None):
        self.lookup_rows = lookup_rows
        self.edge_rows = ["edge-1"] if edge_rows is None else edge_rows
        self.fail_on_update_number = fail_on_update_number
        self.update_count = 0
        self.calls = []

    def execute(self, query, bind_vars=None, **kwargs):
        self.calls.append({"query": query, "bind_vars": bind_vars or {}})
        if "FILTER d.id" in query:
            return FakeCursor(self.lookup_rows)
        if "RETURN e._key" in query:
            return FakeCursor(self.edge_rows)
        if "UPDATE @target_key" in query or "REPLACE @target_key" in query:
            self.update_count += 1
            if self.update_count == self.fail_on_update_number:
                raise RuntimeError("simulated mutation failure")
        return FakeCursor(["updated"])


class CurationCollection:
    def __init__(self):
        self.deleted = []

    def delete_many(self, keys):
        self.deleted.append(keys)


class CurationDb:
    def __init__(self, lookup_rows, edge_rows=None, fail_on_update_number=None):
        self.aql = CurationAql(
            lookup_rows,
            edge_rows=edge_rows,
            fail_on_update_number=fail_on_update_number,
        )
        self.edge_collection = CurationCollection()
        self.committed = False
        self.aborted = False

    def has_collection(self, _name):
        return True

    def collection(self, _name):
        return self.edge_collection

    def begin_transaction(self, write):
        self.transaction_collections = write
        return self

    def commit_transaction(self):
        self.committed = True

    def abort_transaction(self):
        self.aborted = True


def curation_snapshot(curation_type, operation):
    resolved = SimpleNamespace(operation=operation, batch_id="batch-1")
    record_property_decisions = record_property_decisions_from_operation(
        curation_type,
        operation,
        batch_id="batch-1",
        published_at="2026-09-24T12:00:00Z",
        published_by={"id": "keith"},
    )
    return SimpleNamespace(
        active_operations=[] if record_property_decisions else [resolved],
        active_record_property_decisions=list(record_property_decisions),
        metadata=lambda: {"curation_type": curation_type, "manifest_revision": 1},
    )


def test_apply_generic_record_property_curation_writes_effective_and_original_values():
    db = CurationDb([{
        "_key": "REFMET:1",
        "id": "REFMET:1",
        "formula": "C6H1005",
        "sources": ["RefMet"],
    }])
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.database_name = "metabolite_harmonization"
    adapter._collection_schemas = {
        "MetaboliteIdentifier": {"fields": {"formula": "str"}},
    }
    adapter.get_db = lambda: db
    operation = {
        "action": "set_properties",
        "target": {
            "kind": "node",
            "curation_set": "metabolite_harmonization",
            "model_type": "MetaboliteIdentifier",
            "id": "REFMET:1",
        },
        "decisions": [{
            "path": ["formula"], "mode": "set",
            "value": "C6H10O5", "observed_value": "C6H1005",
        }],
        "note": "Correct a source typo",
    }

    result = adapter.apply_curation_snapshots({
        METABOLITE_RECORD_PROPERTIES: curation_snapshot(
            METABOLITE_RECORD_PROPERTIES, operation
        ),
    })

    patch = db.aql.calls[1]["bind_vars"]["patch"]
    assert "REPLACE @target_key" in db.aql.calls[1]["query"]
    assert patch["formula"] == "C6H10O5"
    assert patch["_curation_original"] == {"formula": "C6H1005"}
    assert patch["sources"][-1].startswith("Manual Curation\tbatch-1")
    assert result["applied"] == 1


def test_apply_chebi_record_curation_is_portable_and_recalculates_structure(monkeypatch):
    db = CurationDb([{
        "_key": "CHEBI:137735",
        "id": "CHEBI:137735",
        "smiles": "bad-smiles",
        "mass": "336.40",
        "calculated_mw": "336.40",
        "structure_components": [{"mw": "336.40"}],
    }])
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.database_name = "standalone_chebi_graph"
    adapter._collection_schemas = {
        "ChemicalEntity": {"fields": {"smiles": "str", "mass": "str"}},
    }
    adapter.get_db = lambda: db
    monkeypatch.setattr(
        "src.core.record_property_curations.calculate_smiles_chemistry",
        lambda smiles, input_field: {
            "calculated_mw": "284.31",
            "calculated_monoisotopic_mass": "284.265337",
            "structure_components": [],
            "structure_calculation_input_field": input_field,
        },
    )
    operation = {
        "action": "set_properties",
        "target": {
            "kind": "node",
            "curation_set": "chebi",
            "model_type": "ChemicalEntity",
            "id": "CHEBI:137735",
        },
        "decisions": [
            {
                "path": ["smiles"],
                "mode": "set",
                "value": "corrected-smiles",
                "observed_value": "bad-smiles",
            },
            {
                "path": ["mass"],
                "mode": "set",
                "value": "284.31",
                "observed_value": "336.40",
            },
        ],
        "note": "Correct deuterium atoms misread as methyl groups",
    }

    adapter.apply_curation_snapshots({
        CHEBI_RECORD_PROPERTIES: curation_snapshot(
            CHEBI_RECORD_PROPERTIES, operation
        ),
    })

    patch = db.aql.calls[1]["bind_vars"]["patch"]
    assert patch["smiles"] == "corrected-smiles"
    assert patch["calculated_mw"] == "284.31"
    assert patch["calculated_monoisotopic_mass"] == "284.265337"
    assert patch["structure_components"] == []


def test_apply_metabolite_structure_curation_recalculates_selected_chem_props():
    corrected_smiles = "[2H]O[2H]"
    db = CurationDb([{
        "_key": "CHEBI:139244",
        "id": "CHEBI:139244",
        "chem_props": [
            {
                "source": "ChEBI",
                "source_id": "CHEBI:139244",
                "iso_smiles": "CC",
                "calculated_mw": "30.07",
                "calculated_monoisotopic_mass": "30.046950192",
                "derived_inchi_key": "OTMSDBZUPAUEDD-UHFFFAOYSA-N",
                "structure_components": [{"smiles": "CC", "mw": "30.07"}],
            },
            {
                "source": "HMDB",
                "source_id": "HMDB0000001",
                "iso_smiles": "C",
                "calculated_mw": "16.043",
            },
        ],
    }])
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.database_name = "metabolite_harmonization"
    adapter._collection_schemas = {
        "MetaboliteIdentifier": {"fields": {"chem_props": {
            "type": "list",
            "item_type": "object",
            "fields": {
                "source": "str",
                "source_id": "str",
                "iso_smiles": "str",
                "calculated_mw": "str",
                "calculated_monoisotopic_mass": "str",
                "derived_inchi_key": "str",
                "structure_components": {"type": "list", "item_type": "object"},
            },
        }}},
    }
    adapter.get_db = lambda: db
    operation = {
        "action": "set_properties",
        "target": {
            "kind": "node",
            "curation_set": "metabolite_harmonization",
            "model_type": "MetaboliteIdentifier",
            "id": "CHEBI:139244",
        },
        "decisions": [{
            "path": [
                "chem_props",
                {"match": {"source": "ChEBI", "source_id": "CHEBI:139244"}},
                "iso_smiles",
            ],
            "mode": "set",
            "value": corrected_smiles,
            "observed_value": "CC",
        }],
        "note": "Restore deuterium atoms",
    }

    adapter.apply_curation_snapshots({
        METABOLITE_RECORD_PROPERTIES: curation_snapshot(
            METABOLITE_RECORD_PROPERTIES, operation
        ),
    })

    patch = db.aql.calls[1]["bind_vars"]["patch"]
    chebi = patch["chem_props"][0]
    assert chebi["iso_smiles"] == corrected_smiles
    assert chebi["calculated_mw"] == "20.027203556"
    assert chebi["calculated_monoisotopic_mass"] == "20.023118176"
    assert chebi["derived_inchi_key"] == "XLYOFNOQVPJJNP-ZSJDYOACSA-N"
    assert chebi["structure_components"][0]["molecular_formula"] == "D2O"
    assert patch["chem_props"][1]["calculated_mw"] == "16.043"


def test_record_property_restore_replaces_persisted_document_without_override_fields():
    db = CurationDb([{
        "_key": "node-1",
        "id": "REFMET:1",
        "formula": "C6H10O5",
        "_curation_original": {"formula": "C6H1005"},
    }])
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.database_name = "metabolite_harmonization"
    adapter._collection_schemas = {
        "MetaboliteIdentifier": {"fields": {"formula": "str"}},
    }
    adapter.get_db = lambda: db
    operation = {
        "action": "set_properties",
        "target": {
            "kind": "node",
            "curation_set": "metabolite_harmonization",
            "model_type": "MetaboliteIdentifier",
            "id": "REFMET:1",
        },
        "decisions": [{"path": ["formula"], "mode": "remove_override"}],
        "note": "Restore the loaded formula",
    }

    adapter.apply_curation_snapshots({
        METABOLITE_RECORD_PROPERTIES: curation_snapshot(
            METABOLITE_RECORD_PROPERTIES, operation
        ),
    })

    mutation = db.aql.calls[1]
    assert "REPLACE @target_key" in mutation["query"]
    assert mutation["bind_vars"]["patch"]["formula"] == "C6H1005"
    assert "_curation_original" not in mutation["bind_vars"]["patch"]


def test_apply_edge_removal_deletes_the_resolved_edge_key():
    db = CurationDb([])
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.get_db = lambda: db
    operation = {
        "action": "remove_edge",
        "edge_type": "MetaboliteIdentifierMappingEdge",
        "start_id": "CHEBI:1",
        "end_id": "HMDB:1",
        "symmetric": True,
    }

    adapter.apply_curation_snapshots({
        METABOLITE_EQUIVALENCE_EDGES: curation_snapshot(
            METABOLITE_EQUIVALENCE_EDGES,
            operation,
        ),
    })

    assert db.edge_collection.deleted == [["edge-1"]]


def test_apply_curations_preflights_all_targets_before_any_mutation():
    db = CurationDb(
        [{"_key": "node-1", "id": "KEGG.COMPOUND:C00001"}],
        edge_rows=[],
    )
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.database_name = "metabolite_harmonization"
    adapter._collection_schemas = {
        "MetaboliteIdentifier": {"fields": {"is_generic_structure": "bool"}},
    }
    adapter.get_db = lambda: db
    property_operation = {
        "action": "set_properties",
        "target": {
            "kind": "node",
            "curation_set": "metabolite_harmonization",
            "model_type": "MetaboliteIdentifier",
            "id": "KEGG.COMPOUND:C00001",
        },
        "decisions": [{
            "path": ["is_generic_structure"], "mode": "set", "value": True,
            "observed_value": None, "observed_exists": False,
        }],
        "note": "Classify generic structure",
    }
    edge_operation = {
        "action": "remove_edge",
        "edge_type": "MetaboliteIdentifierMappingEdge",
        "start_id": "CHEBI:1",
        "end_id": "HMDB:1",
        "symmetric": True,
    }

    with pytest.raises(RuntimeError, match="is missing"):
        adapter.apply_curation_snapshots({
            METABOLITE_RECORD_PROPERTIES: curation_snapshot(
                METABOLITE_RECORD_PROPERTIES, property_operation
            ),
            METABOLITE_EQUIVALENCE_EDGES: curation_snapshot(
                METABOLITE_EQUIVALENCE_EDGES, edge_operation
            ),
        })

    assert not any("REPLACE @target_key" in call["query"] for call in db.aql.calls)
    assert db.edge_collection.deleted == []


def test_apply_curations_aborts_transaction_when_a_later_mutation_fails():
    db = CurationDb(
        [{"_key": "node-1", "id": "KEGG.COMPOUND:C00001"}],
        fail_on_update_number=2,
    )
    adapter = ArangoOutputAdapter.__new__(ArangoOutputAdapter)
    adapter.database_name = "metabolite_harmonization"
    adapter._collection_schemas = {
        "MetaboliteIdentifier": {"fields": {"is_generic_structure": "bool"}},
    }
    adapter.get_db = lambda: db
    operations = [
        {
                "action": "set_properties",
                "target": {
                    "kind": "node",
                    "curation_set": "metabolite_harmonization",
                    "model_type": "MetaboliteIdentifier",
                    "id": identifier,
                },
                "decisions": [{
                    "path": ["is_generic_structure"], "mode": "set", "value": value,
                    "observed_value": None, "observed_exists": False,
                }],
                "note": "Classify generic structure",
            }
        for identifier, value in [
            ("KEGG.COMPOUND:C00001", True),
            ("KEGG.COMPOUND:C00002", False),
        ]
    ]
    snapshot = SimpleNamespace(
        active_operations=[],
        active_record_property_decisions=[
            record_property_decisions_from_operation(
                METABOLITE_RECORD_PROPERTIES,
                operation,
                batch_id="batch-1",
                published_at="2026-09-24T12:00:00Z",
                published_by=None,
            )[0]
            for operation in operations
        ],
        metadata=lambda: {"curation_type": METABOLITE_RECORD_PROPERTIES},
    )

    with pytest.raises(RuntimeError, match="simulated mutation failure"):
        adapter.apply_curation_snapshots({METABOLITE_RECORD_PROPERTIES: snapshot})

    assert db.aborted is True
    assert db.committed is False
