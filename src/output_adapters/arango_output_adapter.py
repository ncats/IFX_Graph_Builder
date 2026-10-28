import dataclasses
import hashlib
import json
import os
import platform
import socket
import time
from types import UnionType
from datetime import datetime, date, timezone
from enum import Enum
from typing import Type, List, get_origin, get_args, Union

from arango.exceptions import DocumentInsertError, DocumentUpdateError
from src.core.decorators import collect_facets, collect_indexed_fields, collect_search_fields
from src.core.curations import (
    METABOLITE_EQUIVALENCE_EDGES,
    is_record_property_type,
    validate_curation_type,
)
from src.core.record_property_curations import (
    apply_record_property_decision,
    recalculate_curated_structure_derivatives,
    schema_for_path,
    validate_value_for_schema,
)
from src.interfaces.metadata import DatabaseMetadata, CollectionMetadata, get_git_metadata
from src.interfaces.output_adapter import OutputAdapter
from src.interfaces.resolver_metadata import resolver_fingerprint_summary
from src.models.datasource_version_info import DataSourceDetails
from src.shared.arango_adapter import ArangoAdapter
from src.shared.record_merger import RecordMerger, FieldConflictBehavior

from src.shared.db_credentials import DBCredentials
from src.infrastructure.object_storage import (
    AwsAssumeRoleCredentials,
    AwsAssumeRoleStorage,
    S3CompatibleStorage,
)

class ArangoOutputAdapter(OutputAdapter, ArangoAdapter):
    NODE_MERGE_METADATA_FIELDS = ("_key", "id", "creation", "updates", "resolved_ids")
    BUILD_DIAGNOSTICS_COLLECTION = "build_diagnostics"

    def __init__(self, credentials, database_name, object_storage_credentials=None):
        self._collection_schemas = {}
        self._graph_views = []
        self._graph_view_source_yaml = None
        self._resolver_fingerprints_by_type = {}
        self._resolver_source_yaml = None
        self._registry_datasets = []
        self._curation_metadata = []
        self._curation_application_reports = []
        self.object_storage = self._object_storage_from_credentials(object_storage_credentials)
        super().__init__(credentials=credentials, database_name=database_name)

    @staticmethod
    def _object_storage_from_credentials(credentials: dict | None):
        if not credentials:
            return None
        if "role_arn" in credentials or credentials.get("type") == "aws_assume_role":
            return AwsAssumeRoleStorage(AwsAssumeRoleCredentials.from_yaml(credentials))
        return S3CompatibleStorage(DBCredentials(**credentials))

    def set_graph_views_metadata(self, graph_views=None, source_yaml=None):
        self._graph_views = graph_views or []
        self._graph_view_source_yaml = source_yaml

    def set_resolver_metadata(self, resolver_fingerprints_by_type=None, source_yaml=None):
        self._resolver_fingerprints_by_type = resolver_fingerprints_by_type or {}
        self._resolver_source_yaml = source_yaml

    def set_registry_dataset_metadata(self, registry_datasets=None):
        self._registry_datasets = registry_datasets or []

    def supports_curations(self) -> bool:
        return True

    def apply_curation_snapshots(self, snapshots: dict, phase: str = "final") -> dict:
        if phase not in {"pre_post", "final"}:
            raise ValueError(f"Unsupported curation application phase: {phase}")
        db = self.get_db()
        plans = []
        reports = []
        record_decisions_by_target = {}
        # Resolve every target before performing the first write. A missing or
        # ambiguous target therefore fails without leaving an earlier curation applied.
        for curation_type, snapshot in snapshots.items():
            if (
                curation_type != METABOLITE_EQUIVALENCE_EDGES
                and not is_record_property_type(curation_type)
            ):
                raise RuntimeError(
                    f"Arango graph builds do not support applying {curation_type!r}"
                )
            if is_record_property_type(curation_type):
                definition = validate_curation_type(curation_type)
                for decision in snapshot.active_record_property_decisions:
                    target = decision.target
                    if (
                        definition.fixed_curation_set is None
                        and target.get("curation_set") != self.database_name
                    ):
                        continue
                    target_key = (target["model_type"], target["id"])
                    record_decisions_by_target.setdefault(target_key, []).append(decision)
                continue
            if phase == "pre_post":
                continue
            for resolved in snapshot.active_operations:
                operation = resolved.operation
                action = operation["action"]
                edge_type = operation["edge_type"]
                if not db.has_collection(edge_type):
                    raise RuntimeError(f"Curation target collection does not exist: {edge_type}")
                matches = list(db.aql.execute(
                    f"""
                    FOR e IN `{edge_type}`
                      FILTER (e.start_id == @left AND e.end_id == @right)
                         OR (@symmetric AND e.start_id == @right AND e.end_id == @left)
                      RETURN e._key
                    """,
                    bind_vars={
                        "left": operation["start_id"],
                        "right": operation["end_id"],
                        "symmetric": operation.get("symmetric") is True,
                    },
                    max_runtime=120,
                ))
                if not matches:
                    raise RuntimeError(
                        f"Curation target edge {edge_type}:{operation['start_id']}→{operation['end_id']} is missing"
                    )
                plans.append({
                    "kind": "edge",
                    "collection_name": edge_type,
                    "keys": matches,
                    "delete": action == "remove_edge",
                    "report": {
                        "curation_type": curation_type,
                        "batch_id": resolved.batch_id,
                        "operation_id": operation.get("operation_id"),
                        "action": action,
                        "target_count": len(matches),
                        "status": "applied" if action == "remove_edge" else "retained",
                    },
                })

        for (collection_name, target_id), decisions in record_decisions_by_target.items():
            if not db.has_collection(collection_name):
                raise RuntimeError(
                    f"Curation target collection does not exist: {collection_name}"
                )
            rows = list(db.aql.execute(
                f"""
                FOR d IN `{collection_name}`
                  FILTER d.id == @target_id
                  LIMIT 2
                  RETURN d
                """,
                bind_vars={"target_id": target_id},
                max_runtime=120,
            ))
            if len(rows) != 1:
                raise RuntimeError(
                    f"Curation target {collection_name}:{target_id} matched {len(rows)} documents"
                )
            schema_fields = (
                (self._collection_schemas.get(collection_name) or {}).get("fields") or {}
            )
            projected = rows[0]
            target_reports = []
            for decision in decisions:
                descriptor = schema_for_path(schema_fields, decision.path)
                if decision.mode == "set":
                    validate_value_for_schema(decision.value, descriptor, decision.path)
                projected, report = apply_record_property_decision(projected, decision)
                report["phase"] = phase
                target_reports.append(report)
            projected = recalculate_curated_structure_derivatives(
                projected,
                model_type=collection_name,
                decisions=decisions,
            )
            patch = {
                key: value
                for key, value in projected.items()
                if not key.startswith("_") or key == "_curation_original"
            }
            plans.append({
                "kind": "record_properties",
                "collection_name": collection_name,
                "target_key": rows[0]["_key"],
                "patch": patch,
                "reports": target_reports,
            })

        transaction = db.begin_transaction(
            write=sorted({plan["collection_name"] for plan in plans})
        ) if plans else None
        write_db = transaction or db
        try:
            for plan in plans:
                if plan["kind"] == "record_properties":
                    write_verb = "REPLACE"
                    write_options = (
                        "" if plan["kind"] == "record_properties"
                        else "OPTIONS {keepNull: true}"
                    )
                    updated = list(write_db.aql.execute(
                        f"""
                        {write_verb} @target_key WITH @patch IN `{plan['collection_name']}`
                          {write_options}
                          RETURN NEW._key
                        """,
                        bind_vars={
                            "target_key": plan["target_key"],
                            "patch": plan["patch"],
                        },
                        max_runtime=120,
                    ))
                    if len(updated) != 1:
                        raise RuntimeError(
                            f"Curation target disappeared during application: "
                            f"{plan['collection_name']}:{plan['target_key']}"
                        )
                elif plan["delete"]:
                    write_db.collection(plan["collection_name"]).delete_many(plan["keys"])
                reports.extend(plan.get("reports") or [plan["report"]])
            if transaction is not None:
                transaction.commit_transaction()
        except Exception:
            if transaction is not None:
                transaction.abort_transaction()
            raise
        metadata_by_type = {
            item.get("curation_type"): item
            for item in getattr(self, "_curation_metadata", [])
        }
        metadata_by_type.update({
            snapshot.metadata().get("curation_type"): snapshot.metadata()
            for snapshot in snapshots.values()
        })
        self._curation_metadata = list(metadata_by_type.values())
        application_reports = getattr(self, "_curation_application_reports", [])
        application_reports.extend(reports)
        self._curation_application_reports = application_reports
        return {
            "applied": sum(report["status"] == "applied" for report in reports),
            "reports": reports,
        }

    @staticmethod
    def _introspect_dataclass(cls) -> dict:
        """Extract a schema dict from a dataclass by inspecting its fields and type hints."""
        skip_fields = {'labels', 'start_node', 'end_node'}
        result = {}

        for f in dataclasses.fields(cls):
            if f.name.startswith('_') or f.name in skip_fields:
                continue
            result[f.name] = ArangoOutputAdapter._type_hint_to_schema(f.type)

        return result

    @staticmethod
    def _type_hint_to_schema(type_hint) -> Union[str, dict]:
        """Convert a Python type hint to a simple schema descriptor."""
        # Unwrap Optional[X] → X
        origin = get_origin(type_hint)
        if origin is Union:
            args = [a for a in get_args(type_hint) if a is not type(None)]
            if len(args) == 1:
                type_hint = args[0]
                origin = get_origin(type_hint)

        # Handle string annotations (forward references)
        if isinstance(type_hint, str):
            return "str"

        # Primitives
        if type_hint is str:
            return "str"
        if type_hint is int:
            return "int"
        if type_hint is float:
            return "float"
        if type_hint is bool:
            return "bool"
        if type_hint is date:
            return "date"
        if type_hint is datetime:
            return "datetime"

        # Enums → "str" (LabeledIntEnum stores its .label, a string, in Arango)
        if isinstance(type_hint, type) and issubclass(type_hint, Enum):
            return "str"

        # List[X]
        if origin is list:
            item_args = get_args(type_hint)
            if item_args:
                item_type = item_args[0]
                # List of dataclasses → nested object
                if isinstance(item_type, type) and dataclasses.is_dataclass(item_type):
                    return {
                        "type": "list",
                        "item_type": "object",
                        "fields": ArangoOutputAdapter._introspect_dataclass(item_type)
                    }
                # List of LabeledIntEnum → labels stored as strings, but record enum class
                # path so the MySQL converter can build a proper integer lookup table
                if isinstance(item_type, type) and issubclass(item_type, int) and issubclass(item_type, Enum):
                    return {
                        "type": "list",
                        "item_type": "str",
                        "enum": f"{item_type.__module__}.{item_type.__qualname__}",
                    }
                inner = ArangoOutputAdapter._type_hint_to_schema(item_type)
                return {"type": "list", "item_type": inner}
            return {"type": "list", "item_type": "str"}

        # Dict[K, V]
        if origin is dict:
            args = get_args(type_hint)
            key_type = ArangoOutputAdapter._type_hint_to_schema(args[0]) if args else "str"
            val_type = ArangoOutputAdapter._type_hint_to_schema(args[1]) if len(args) > 1 else "str"
            return {"type": "dict", "key_type": key_type, "value_type": val_type}

        # Dataclass objects
        if isinstance(type_hint, type) and dataclasses.is_dataclass(type_hint):
            return {
                "type": "object",
                "fields": ArangoOutputAdapter._introspect_dataclass(type_hint)
            }

        return "str"

    def _ensure_bucket(self, s3_client, bucket: str):
        from botocore.exceptions import ClientError
        try:
            s3_client.head_bucket(Bucket=bucket)
        except ClientError as e:
            if e.response['Error']['Code'] in ('404', 'NoSuchBucket', '403'):
                s3_client.create_bucket(Bucket=bucket)

    def _handle_dataset_nodes(self, objects):
        """Upload DataFrame from any Dataset or StatsResult node to object storage as Parquet."""
        from src.models.pounce.dataset import Dataset
        from src.models.pounce.stats_result import StatsResult

        if not self.object_storage:
            return

        import io
        import pyarrow as pa
        import pyarrow.parquet as pq

        storage = self.object_storage
        s3 = storage.client()
        storage.ensure_bucket()

        for obj in objects:
            if isinstance(obj, (Dataset, StatsResult)) and obj._data_frame is not None:
                safe_id = self.safe_key(obj.id).replace(':', '_')
                key = f"{self.database_name}/parquet/{safe_id}.parquet"

                buf = io.BytesIO()
                table = pa.Table.from_pandas(obj._data_frame)
                pq.write_table(table, buf)
                buf.seek(0)

                s3.put_object(Bucket=storage.bucket, Key=key, Body=buf)
                print(f"Uploaded Parquet to object storage: s3://{storage.bucket}/{key} ({obj.row_count} rows x {obj.column_count} cols)")

                obj.file_reference = f"s3://{storage.bucket}/{key}"
                obj._data_frame = None

    def _project_id_for_workbook_owner(self, obj) -> str:
        from src.models.pounce.project import Project
        from src.models.pounce.experiment import Experiment
        from src.models.pounce.stats_result import StatsResult

        if isinstance(obj, Project):
            return obj.id
        if isinstance(obj, Experiment):
            return obj.id.rsplit("_", 1)[0]
        if isinstance(obj, StatsResult):
            experiment_id = obj.id.split(":", 1)[0]
            return experiment_id.rsplit("_", 1)[0]
        raise TypeError(f"Unsupported workbook owner type: {type(obj).__name__}")

    def _handle_pounce_workbook_nodes(self, objects):
        from src.models.pounce.project import Project
        from src.models.pounce.experiment import Experiment
        from src.models.pounce.stats_result import StatsResult
        from src.models.pounce.workbook_artifact import WorkbookArtifact

        if not self.object_storage:
            return

        storage = self.object_storage
        s3 = storage.client()
        storage.ensure_bucket()

        for obj in objects:
            if not isinstance(obj, (Project, Experiment, StatsResult)):
                continue
            workbook = getattr(obj, "workbook", None)
            if not isinstance(workbook, WorkbookArtifact):
                continue
            local_path = workbook._local_path
            if workbook.file_reference and str(workbook.file_reference).startswith("s3://"):
                continue
            if not local_path:
                continue
            if not os.path.exists(local_path):
                raise FileNotFoundError(f"Workbook file not found for {type(obj).__name__} {obj.id}: {local_path}")

            project_id = self._project_id_for_workbook_owner(obj)
            original_name = workbook.original_filename or os.path.basename(local_path)
            key = f"{self.database_name}/workbooks/{self.safe_key(project_id)}/{original_name}"
            with open(local_path, "rb") as handle:
                s3.put_object(Bucket=storage.bucket, Key=key, Body=handle)
            workbook.original_filename = original_name
            workbook.file_reference = f"s3://{storage.bucket}/{key}"
            workbook._local_path = None

    def create_indexes(self, cls: Type, collection):
        indexed_fields = collect_indexed_fields(cls)
        categories, numerics = collect_facets(cls)
        # Dataclass fields retain the declared types without evaluating unrelated
        # forward references (for example ProteinPathwayEdge.start_node).
        field_types = {field.name: field.type for field in dataclasses.fields(cls)}

        def index_path(field: str) -> str:
            if "[*]" in field:
                return field
            hint = field_types.get(field)
            if get_origin(hint) in (Union, UnionType):
                members = [member for member in get_args(hint) if member is not type(None)]
                if len(members) == 1:
                    hint = members[0]
            return f"{field}[*]" if get_origin(hint) is list else field

        existing_indexes = collection.indexes()
        existing_fields = {tuple(index['fields']) for index in existing_indexes}

        # Index scalar values directly and individual members of marked list fields.
        for field in sorted({index_path(field) for field in indexed_fields | categories}):
            field_tuple = (field,)
            if field_tuple not in existing_fields:
                if "[*]" in field:
                    print(f"Creating PERSISTENT array index on: {field}")
                    collection.add_persistent_index(fields=[field], sparse=True)
                elif hasattr(collection, "add_hash_index"):
                    print(f"Creating HASH index on: {field}")
                    collection.add_hash_index(fields=[field], sparse=True)
                else:
                    print(f"Creating PERSISTENT index on: {field}")
                    collection.add_persistent_index(fields=[field], sparse=True)
                existing_fields.add(field_tuple)

        # Create persistent indexes for numeric fields
        for field in sorted({index_path(field) for field in numerics}):
            field_tuple = (field,)
            if field_tuple not in existing_fields:
                print(f"Creating PERSISTENT index on: {field}")
                collection.add_persistent_index(fields=[field], sparse=True)
                existing_fields.add(field_tuple)


    @staticmethod
    def document_handle(collection_name: str, semantic_id: str) -> str:
        return f"{collection_name}/{ArangoOutputAdapter.safe_key(semantic_id)}"


    def store(self, objects, single_source = False,
              field_conflict_behavior: FieldConflictBehavior = FieldConflictBehavior.KeepFirst) -> bool:

        def generate_edge_key(from_node, to_node, edge_type):
            edge_hash = hashlib.sha256(f"{from_node}|{to_node}".encode("utf-8")).hexdigest()
            return f"{self.safe_key(edge_type)}__{edge_hash}"

        def get_insert_failures(insert_result):
            failed = []
            if not insert_result:
                return failed
            for result in insert_result:
                if isinstance(result, Exception):
                    failed.append(result)
                elif isinstance(result, dict) and result.get("error"):
                    failed.append(result)
            return failed

        merger = RecordMerger(field_conflict_behavior=field_conflict_behavior)

        if not isinstance(objects, list):
            objects = [objects]

        self._handle_dataset_nodes(objects)
        self._handle_pounce_workbook_nodes(objects)

        db = self.get_db()
        graph = self.get_graph()
        object_groups = self.sort_and_convert_objects(objects, convert_dates=True)

        # Collect schema info from each object group
        for obj_list, labels, is_relationship, start_labels, end_labels, obj_cls in object_groups.values():
            label = labels[0]
            indexed_fields = collect_indexed_fields(obj_cls)
            categories, numerics = collect_facets(obj_cls)
            text_fields = collect_search_fields(obj_cls)
            if label not in self._collection_schemas:
                schema_entry = {
                    "fields": self._introspect_dataclass(obj_cls),
                    "index_metadata": {
                        "fields": sorted(indexed_fields),
                    },
                    "facet_metadata": {
                        "category_fields": sorted(categories),
                        "numeric_fields": sorted(numerics),
                    },
                    "search_metadata": {
                        "text_fields": sorted(text_fields),
                    },
                }
                if is_relationship:
                    schema_entry["type"] = "edge"
                    schema_entry["from_collections"] = start_labels
                    schema_entry["to_collections"] = end_labels
                else:
                    schema_entry["type"] = "document"
                self._collection_schemas[label] = schema_entry
            else:
                existing = self._collection_schemas[label]
                existing_index_metadata = existing.setdefault("index_metadata", {
                    "fields": [],
                })
                existing_facet_metadata = existing.setdefault("facet_metadata", {
                    "category_fields": [],
                    "numeric_fields": [],
                })
                existing_search_metadata = existing.setdefault("search_metadata", {
                    "text_fields": [],
                })
                existing_index_metadata["fields"] = sorted(
                    set(existing_index_metadata.get("fields", [])) | indexed_fields
                )
                existing_facet_metadata["category_fields"] = sorted(
                    set(existing_facet_metadata.get("category_fields", [])) | categories
                )
                existing_facet_metadata["numeric_fields"] = sorted(
                    set(existing_facet_metadata.get("numeric_fields", [])) | numerics
                )
                existing_search_metadata["text_fields"] = sorted(
                    set(existing_search_metadata.get("text_fields", [])) | text_fields
                )
            if is_relationship and label in self._collection_schemas:
                # Merge from/to collections for edges seen from multiple sources
                existing = self._collection_schemas[label]
                existing["from_collections"] = list(set(existing.get("from_collections", []) + start_labels))
                existing["to_collections"] = list(set(existing.get("to_collections", []) + end_labels))

        for obj_list, labels, is_relationship, start_labels, end_labels, obj_cls in object_groups.values():
            label = labels[0]
            if is_relationship:
                if not graph.has_edge_collection(label):
                    edge_collection = graph.create_edge_definition(label, start_labels, end_labels)
                else:
                    edge_definition = [definition for definition in graph.edge_definitions() if definition['edge_collection'] == label][0]
                    updated_from = list(set(edge_definition['from_vertex_collections'] + start_labels))
                    updated_to = list(set(edge_definition['to_vertex_collections'] + end_labels))
                    edge_collection = graph.replace_edge_definition(label, updated_from, updated_to)

                self.create_indexes(obj_cls, edge_collection)
                keys = [generate_edge_key(obj['start_id'], obj['end_id'], label) for obj in obj_list]

                existing_record_map = {}
                if not single_source:
                    existing_edges = edge_collection.get_many(keys)
                    existing_record_map = {
                        (record['start_id'], record['end_id']): record for record in existing_edges
                    }
                merged_records = merger.merge_records(obj_list, existing_record_map, nodes_or_edges='edges')

                edges = []
                for obj in merged_records:
                    edge = {
                        **obj,
                        "_from": self.document_handle(start_labels[0], obj["start_id"]),
                        "_to": self.document_handle(end_labels[0], obj["end_id"]),
                        "_key": generate_edge_key(obj["start_id"], obj["end_id"], label)
                    }
                    edges.append(edge)

                self.insert_many_with_backoff(
                    edge_collection,
                    edges,
                    overwrite=True,
                    label=label,
                    kind="edge"
                )
            else:
                collection, existing_nodes = self.get_existing_nodes(db, label, obj_list, skip_merge=single_source)

                self.create_indexes(obj_cls, collection)
                existing_keys = {record['_key'] for record in existing_nodes}
                existing_record_map = {record['id']: record for record in existing_nodes}
                merged_nodes = merger.merge_records(obj_list, existing_record_map, nodes_or_edges='nodes')

                node_payloads = [{**obj, "_key": self.safe_key(obj["id"])} for obj in merged_nodes]

                if single_source:
                    self.insert_many_with_backoff(
                        collection,
                        node_payloads,
                        overwrite=True,
                        label=label,
                        kind="node"
                    )
                    continue

                inserts = [obj for obj in node_payloads if obj["_key"] not in existing_keys]
                updates = [obj for obj in node_payloads if obj["_key"] in existing_keys]

                if inserts:
                    self.insert_many_with_backoff(
                        collection,
                        inserts,
                        overwrite=False,
                        label=label,
                        kind="node"
                    )

                if updates:
                    self.update_many_with_backoff(
                        collection,
                        updates,
                        label=label,
                        kind="node"
                    )

        return True

    def get_existing_nodes(self, db, label, obj_list, skip_merge = False):
        if not db.has_collection(label):
            collection = db.create_collection(label)
        else:
            collection = db.collection(label)
        if skip_merge:
            return collection, []
        keys = [self.safe_key(obj['id']) for obj in obj_list]
        fetch_fields = self.get_node_merge_fetch_fields(obj_list)
        existing_nodes = []
        # Keep each result within one response so a long merge does not depend
        # on a server cursor surviving across fetch requests.
        for offset in range(0, len(keys), 1000):
            key_batch = keys[offset:offset + 1000]
            existing_nodes.extend(db.aql.execute(
                """
                FOR key IN @keys
                  LET doc = DOCUMENT(@collection, key)
                  FILTER doc != null
                  RETURN KEEP(doc, @fields)
                """,
                bind_vars={
                    "collection": label,
                    "keys": key_batch,
                    "fields": fetch_fields,
                },
                batch_size=len(key_batch) + 1,
                max_runtime=600,
            ))
        return collection, existing_nodes

    def get_node_merge_fetch_fields(self, obj_list):
        fields = set(self.NODE_MERGE_METADATA_FIELDS)
        for obj in obj_list:
            for key in obj.keys():
                if key.startswith('_'):
                    continue
                if key in {"provenance", "entity_resolution"}:
                    continue
                fields.add(key)
        return sorted(fields)

    def insert_many_with_backoff(self, collection, records, overwrite=False, label="", kind="node"):
        self._insert_many_with_backoff(
            collection,
            records,
            overwrite=overwrite,
            label=label,
            kind=kind
        )

    def _insert_many_with_backoff(self, collection, records, overwrite=False, label="", kind="node"):
        try:
            insert_result = collection.insert_many(records, overwrite=overwrite)
        except DocumentInsertError:
            if len(records) == 1:
                raise
            midpoint = len(records) // 2
            print(
                f"{kind} insert batch failed for {label}; "
                f"splitting {len(records)} records into {midpoint} and {len(records) - midpoint}"
            )
            self._insert_many_with_backoff(collection, records[:midpoint], overwrite=overwrite, label=label, kind=kind)
            self._insert_many_with_backoff(collection, records[midpoint:], overwrite=overwrite, label=label, kind=kind)
            return

        failed = []
        if insert_result:
            for result in insert_result:
                if isinstance(result, DocumentInsertError):
                    failed.append(result)
                elif isinstance(result, dict) and result.get("error"):
                    failed.append(result)
        if failed:
            print(f"{kind} insert failures for {label}:")
            for failure in failed[:10]:
                print(failure)
            raise Exception(f"failed to insert {len(failed)} {kind} records into {label}")

    def update_many_with_backoff(self, collection, records, label="", kind="node"):
        self._update_many_with_backoff(collection, records, label=label, kind=kind)

    def _update_many_with_backoff(self, collection, records, label="", kind="node"):
        try:
            update_result = collection.update_many(
                records,
                merge=True,
                keep_none=False,
                check_rev=False
            )
        except DocumentUpdateError:
            if len(records) == 1:
                raise
            midpoint = len(records) // 2
            print(
                f"{kind} update batch failed for {label}; "
                f"splitting {len(records)} records into {midpoint} and {len(records) - midpoint}"
            )
            self._update_many_with_backoff(collection, records[:midpoint], label=label, kind=kind)
            self._update_many_with_backoff(collection, records[midpoint:], label=label, kind=kind)
            return

        failed = []
        if update_result:
            for result in update_result:
                if isinstance(result, DocumentUpdateError):
                    failed.append(result)
                elif isinstance(result, dict) and result.get("error"):
                    failed.append(result)
        if failed:
            print(f"{kind} update failures for {label}:")
            for failure in failed[:10]:
                print(failure)
            raise Exception(f"failed to update {len(failed)} {kind} records into {label}")

    def _delete_object_storage_prefix(self, prefix: str) -> None:
        if not self.object_storage:
            return

        from botocore.exceptions import ClientError

        storage = self.object_storage
        s3 = storage.client()

        try:
            s3.head_bucket(Bucket=storage.bucket)
        except ClientError:
            return

        paginator = s3.get_paginator("list_objects_v2")
        total_deleted = 0
        for page in paginator.paginate(Bucket=storage.bucket, Prefix=prefix):
            contents = page.get("Contents", [])
            if not contents:
                continue
            objects = [{"Key": item["Key"]} for item in contents]
            for i in range(0, len(objects), 1000):
                chunk = objects[i:i + 1000]
                s3.delete_objects(Bucket=storage.bucket, Delete={"Objects": chunk})
                total_deleted += len(chunk)
        if total_deleted:
            print(f"Deleted {total_deleted} object storage objects under s3://{storage.bucket}/{prefix}")

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        sys_db = self.client.db('_system', username=self.credentials.user,
                            password=self.credentials.password)

        effective_truncate = True if truncate_tables is None else truncate_tables

        if sys_db.has_database(self.database_name):
            if effective_truncate:
                sys_db.delete_database(self.database_name)
                sys_db.create_database(self.database_name)
                self._delete_object_storage_prefix(f"{self.database_name}/")
        else:
            sys_db.create_database(self.database_name)

        return True

    def _checkpoint_doc_key(self, run_id: str) -> str:
        return f"etl_checkpoint__{self.safe_key(run_id)}"

    def _run_lock_key(self, run_id: str) -> str:
        return f"etl_run_lock__{self.safe_key(run_id)}"

    @staticmethod
    def _owner_is_live(owner: dict) -> bool:
        if owner.get("hostname") != socket.gethostname() or not isinstance(owner.get("pid"), int):
            return True
        try:
            os.kill(owner["pid"], 0)
        except ProcessLookupError:
            return False
        except PermissionError:
            return True
        return True

    def acquire_run_lock(self, run_id: str, attempt_id: str) -> None:
        store = self.get_metadata_store(truncate=False)
        key = self._run_lock_key(run_id)
        payload = {
            "_key": key, "type": "etl_run_lock", "run_id": run_id,
            "attempt_id": attempt_id, "hostname": socket.gethostname(), "pid": os.getpid(),
            "started_at": datetime.now(timezone.utc).isoformat(),
        }
        try:
            store.insert(payload)
            return
        except DocumentInsertError:
            owner = store.get(key)
            if not owner or self._owner_is_live(owner):
                raise RuntimeError(f"Build {run_id} already has an active or unverifiable run lock")
            # Revision checking allows only one contender to recover a dead process's lock.
            self.get_db().aql.execute(
                "REMOVE { _key: @key, _rev: @rev } IN metadata_store OPTIONS { ignoreRevs: false }",
                bind_vars={"key": key, "rev": owner["_rev"]},
            )
            store.insert(payload)

    def release_run_lock(self, run_id: str, attempt_id: str) -> None:
        key = self._run_lock_key(run_id)
        store = self.get_metadata_store(truncate=False)
        owner = store.get(key)
        if owner and owner.get("attempt_id") == attempt_id:
            self.get_db().aql.execute(
                "REMOVE { _key: @key, _rev: @rev } IN metadata_store OPTIONS { ignoreRevs: false }",
                bind_vars={"key": key, "rev": owner["_rev"]},
            )

    def _read_checkpoint_doc(self, run_id: str) -> dict:
        store = self.get_metadata_store(truncate=False)
        return store.get(self._checkpoint_doc_key(run_id)) or {}

    def _write_checkpoint_doc(self, run_id: str, doc: dict) -> None:
        store = self.get_metadata_store(truncate=False)
        previous = store.get(self._checkpoint_doc_key(run_id)) or {}
        payload = {
            **{key: value for key, value in previous.items() if not key.startswith("_")},
            "_key": self._checkpoint_doc_key(run_id),
            "type": "etl_checkpoint",
            "run_id": run_id,
            **doc,
            "last_updated": datetime.now().isoformat(),
        }
        store.insert(payload, overwrite=True)

    def begin_run_attempt(self, run_id: str, resume: bool, attempt_id: str) -> None:
        doc = self._read_checkpoint_doc(run_id)
        attempts = list(doc.get("attempts") or [])
        now = datetime.now(timezone.utc).isoformat()
        if attempts and attempts[-1].get("status") == "running":
            if not resume:
                raise RuntimeError(f"Build {run_id} already has an active attempt")
            owner = attempts[-1]
            if self._owner_is_live(owner):
                raise RuntimeError(f"Build {run_id} is still owned by process {owner.get('pid', 'unknown')}")
            attempts[-1] = {**attempts[-1], "status": "interrupted"}
        attempts.append({
            "id": attempt_id, "resume": resume, "started_at": now,
            "status": "running", "hostname": socket.gethostname(), "pid": os.getpid(),
        })
        self._write_checkpoint_doc(run_id, {"attempts": attempts})

    def finish_run_attempt(self, run_id: str, attempt_id: str, outcome: str, elapsed_seconds: float,
                           error_message: str | None = None) -> None:
        doc = self._read_checkpoint_doc(run_id)
        attempts = list(doc.get("attempts") or [])
        if not attempts or attempts[-1].get("id") != attempt_id or attempts[-1].get("status") != "running":
            raise RuntimeError(f"Build attempt {attempt_id} is no longer active")
        attempts[-1] = {
            **attempts[-1], "status": outcome, "ended_at": datetime.now(timezone.utc).isoformat(),
            "elapsed_seconds": elapsed_seconds, "error_message": error_message,
        }
        self._write_checkpoint_doc(run_id, {"attempts": attempts})

    def _upsert_collection_schemas_doc(self) -> None:
        if not self._collection_schemas:
            return
        store = self.get_metadata_store(truncate=False)
        existing_doc = store.get("collection_schemas") or {}
        merged_collections = {**existing_doc.get("collections", {}), **self._collection_schemas}
        store.insert({
            "_key": "collection_schemas",
            "collections": merged_collections,
        }, overwrite=True)

    def get_completed_adapter_names(self, run_id: str) -> set[str]:
        doc = self._read_checkpoint_doc(run_id)
        adapters = doc.get("adapters", {}) or {}
        return {
            adapter_name for adapter_name, metadata in adapters.items()
            if metadata.get("status") == "completed"
        }

    def reset_run_state(self, run_id: str) -> None:
        store = self.get_metadata_store(truncate=False)
        existing = store.get(self._checkpoint_doc_key(run_id))
        if existing:
            store.delete(existing["_key"])
        if self.get_db().has_collection(self.BUILD_DIAGNOSTICS_COLLECTION):
            self.get_db().aql.execute(
                "FOR d IN build_diagnostics FILTER d.run_id == @run_id REMOVE d IN build_diagnostics",
                bind_vars={"run_id": run_id},
                max_runtime=600,
            )

    def _build_diagnostics_collection(self):
        db = self.get_db()
        name = self.BUILD_DIAGNOSTICS_COLLECTION
        if not db.has_collection(name):
            collection = db.create_collection(name)
            collection.add_persistent_index(fields=["run_id", "kind", "adapter"])
            return collection
        return db.collection(name)

    @staticmethod
    def _diagnostic_key(parts: list[str]) -> str:
        return hashlib.sha256(json.dumps(parts, ensure_ascii=False).encode("utf-8")).hexdigest()

    def set_build_run_id(self, run_id: str) -> None:
        self._build_run_id = run_id
        self._build_diagnostics_collection()
        self.get_metadata_store(truncate=False).insert(
            {
                "_key": f"build_diagnostics__{self.safe_key(run_id)}",
                "type": "build_diagnostics_status",
                "run_id": run_id,
                "status": "partial",
                "updated_at": datetime.now(timezone.utc).isoformat(),
            },
            overwrite=True,
        )

    def record_unmatched_ids(self, run_id: str, adapter_name: str, phase: str, issues: list[dict]) -> None:
        if not issues:
            return
        self._build_diagnostics_collection()
        now = datetime.now(timezone.utc).isoformat()
        records = []
        for issue in issues:
            records.append({
                "_key": self._diagnostic_key([
                    "unmatched_id", run_id, adapter_name, phase,
                    issue["node_type"], issue["source_id"], issue["endpoint_role"],
                ]),
                "kind": "unmatched_id",
                "run_id": run_id,
                "phase": phase,
                "adapter": adapter_name,
                "node_type": issue["node_type"],
                "source_id": issue["source_id"],
                "endpoint_role": issue["endpoint_role"],
                "occurrences": issue["occurrences"],
                "last_seen_at": now,
            })
        for start in range(0, len(records), 1000):
            self.get_db().aql.execute(
                """
                FOR issue IN @issues
                    UPSERT { _key: issue._key }
                    INSERT issue
                    UPDATE {
                        occurrences: OLD.occurrences + issue.occurrences,
                        last_seen_at: issue.last_seen_at
                    } IN build_diagnostics
                """,
                bind_vars={"issues": records[start:start + 1000]},
                max_runtime=600,
            )

    def _record_dangling_edges(self, run_id: str, collection_name: str, edges: list[dict], status: str) -> None:
        if not edges:
            return
        self._build_diagnostics_collection()
        records = [{
            "_key": self._diagnostic_key(["dangling_edge", run_id, collection_name, edge["_key"]]),
            "kind": "dangling_edge",
            "run_id": run_id,
            "phase": "cleanup",
            "edge_collection": collection_name,
            "edge_key": edge["_key"],
            "from_id": edge["_from"],
            "to_id": edge["_to"],
            "missing_from": edge["missing_from"],
            "missing_to": edge["missing_to"],
            "status": status,
        } for edge in edges]
        for start in range(0, len(records), 1000):
            self.get_db().aql.execute(
                """
                FOR issue IN @issues
                    UPSERT { _key: issue._key }
                    INSERT issue
                    UPDATE { status: issue.status } IN build_diagnostics
                """,
                bind_vars={"issues": records[start:start + 1000]},
                max_runtime=600,
            )

    def finish_build_diagnostics(self, run_id: str) -> None:
        self._reconcile_observed_dangling_edges(run_id)
        counts = self.get_db().aql.execute(
            """
            FOR d IN build_diagnostics
                FILTER d.run_id == @run_id
                COLLECT kind = d.kind, status = d.status WITH COUNT INTO count
                RETURN {kind, status, count}
            """,
            bind_vars={"run_id": run_id},
            max_runtime=600,
        )
        summary = {f"{row['kind']}:{row.get('status') or 'recorded'}": row['count'] for row in counts}
        if summary.get("dangling_edge:observed", 0):
            raise RuntimeError("Build diagnostics still contain dangling edges whose deletion was not confirmed")
        self.get_metadata_store(truncate=False).insert(
            {
                "_key": f"build_diagnostics__{self.safe_key(run_id)}",
                "type": "build_diagnostics_status",
                "run_id": run_id,
                "status": "complete",
                "summary": summary,
                "updated_at": datetime.now(timezone.utc).isoformat(),
            },
            overwrite=True,
        )
        print(f"Build diagnostics: {self.BUILD_DIAGNOSTICS_COLLECTION} {summary}")

    def _reconcile_observed_dangling_edges(self, run_id: str) -> None:
        db = self.get_db()
        observed = db.aql.execute(
            "FOR d IN build_diagnostics FILTER d.run_id == @run_id AND d.kind == 'dangling_edge' AND d.status == 'observed' RETURN d",
            bind_vars={"run_id": run_id},
            max_runtime=600,
        )
        missing_keys = []
        for issue in observed:
            collection_name = issue["edge_collection"]
            if not db.has_collection(collection_name) or db.collection(collection_name).get(issue["edge_key"]) is None:
                missing_keys.append(issue["_key"])
        for start in range(0, len(missing_keys), 1000):
            db.aql.execute(
                "FOR key IN @keys UPDATE key WITH { status: 'deleted' } IN build_diagnostics",
                bind_vars={"keys": missing_keys[start:start + 1000]},
                max_runtime=600,
            )

    def mark_adapter_running(self, run_id: str, adapter_name: str, adapter_position: int | None = None,
                             adapter_total: int | None = None) -> None:
        self._build_diagnostics_collection()
        self.get_db().aql.execute(
            """
            FOR d IN build_diagnostics
                FILTER d.run_id == @run_id AND d.kind == "unmatched_id" AND d.adapter == @adapter
                REMOVE d IN build_diagnostics
            """,
            bind_vars={"run_id": run_id, "adapter": adapter_name},
            max_runtime=600,
        )
        doc = self._read_checkpoint_doc(run_id)
        adapters = dict(doc.get("adapters", {}) or {})
        previous = dict(adapters.get(adapter_name, {}) or {})
        adapters[adapter_name] = {
            **previous,
            "status": "running",
            "adapter_position": adapter_position,
            "adapter_total": adapter_total,
            "started_at": datetime.now(timezone.utc).isoformat(),
        }
        self._write_checkpoint_doc(run_id, {"adapters": adapters})

    def mark_adapter_completed(self, run_id: str, adapter_name: str, records_written: int = 0,
                               adapter_position: int | None = None, adapter_total: int | None = None) -> None:
        doc = self._read_checkpoint_doc(run_id)
        adapters = dict(doc.get("adapters", {}) or {})
        previous = dict(adapters.get(adapter_name, {}) or {})
        adapters[adapter_name] = {
            **previous,
            "status": "completed",
            "adapter_position": adapter_position,
            "adapter_total": adapter_total,
            "records_written": records_written,
            "completed_at": datetime.now(timezone.utc).isoformat(),
        }
        self._write_checkpoint_doc(run_id, {"adapters": adapters})

    def mark_adapter_failed(self, run_id: str, adapter_name: str, error_message: str | None = None,
                            adapter_position: int | None = None, adapter_total: int | None = None) -> None:
        doc = self._read_checkpoint_doc(run_id)
        adapters = dict(doc.get("adapters", {}) or {})
        previous = dict(adapters.get(adapter_name, {}) or {})
        adapters[adapter_name] = {
            **previous,
            "status": "failed",
            "adapter_position": adapter_position,
            "adapter_total": adapter_total,
            "error_message": error_message,
            "failed_at": datetime.now(timezone.utc).isoformat(),
        }
        self._write_checkpoint_doc(run_id, {"adapters": adapters})

    def flush_incremental_metadata(self) -> None:
        self._upsert_collection_schemas_doc()


    def get_metadata(self) -> DatabaseMetadata:
        collections: List[CollectionMetadata] = []
        ignore_list = [self.metadata_store_label, self.BUILD_DIAGNOSTICS_COLLECTION]
        ignore_prefixes = ("MetaboliteHarmonizationClique",)

        db = self.get_db()

        for collection in db.collections():
            if collection['system']:
                continue
            name = collection['name']
            if name in ignore_list:
                continue
            if any(name.startswith(prefix) for prefix in ignore_prefixes):
                continue
            count_query = f"""
                    RETURN COUNT(
                        FOR doc IN `{name}`
                        RETURN 1)
                    """

            cursor = db.aql.execute(count_query)
            coll_obj = CollectionMetadata(name=name, total_count=cursor.pop())

            source_count_query = f"""
                    FOR doc IN `{collection['name']}`
                        LET values = UNIQUE(doc.sources || [])
                        FOR item IN values
                            COLLECT value = item WITH COUNT INTO count
                        SORT count DESC
                        RETURN {{ value, count }}
                    """

            cursor = db.aql.execute(source_count_query)
            res = list(cursor)

            for row in res:
                count = row['count']
                source_tsv = row['value']
                if not source_tsv:
                    continue
                dsd = DataSourceDetails.parse_tsv(source_tsv)
                coll_obj.sources.append(dsd)
                coll_obj.marginal_source_counts[dsd.name] = count

            upset_count_query = f"""FOR doc IN `{collection['name']}`
                LET agg = doc.sources
                LET sortedAgg = SORTED(agg)
                LET key = CONCAT_SEPARATOR("|", sortedAgg)
                COLLECT combo = key WITH COUNT INTO count
                SORT count DESC
                RETURN {{
                    combination: combo,
                    count: count
            }}"""

            cursor = db.aql.execute(upset_count_query)
            res = list(cursor)

            for row in res:
                count = row['count']
                combined_source_tsv = row['combination']
                sources: List[str] = []
                for source_tsv in combined_source_tsv.split('|'):
                    if not source_tsv:
                        continue
                    sources.append(DataSourceDetails.parse_tsv(source_tsv).name)
                if not sources:
                    continue
                coll_obj.joint_source_counts['|'.join(sources)] = count

            collections.append(coll_obj)

        return DatabaseMetadata(collections=collections)

    def do_post_processing(self, clean_edges: bool = True) -> None:
        if clean_edges:
            self.clean_up_dangling_edges()

        existing_store = self.get_metadata_store(truncate=False)
        existing_graph_views = {}
        existing_resolver_metadata = {}
        existing_registry_datasets = []
        try:
            existing_doc = existing_store.get("collection_schemas")
            if existing_doc:
                self._collection_schemas = {**existing_doc.get("collections", {}), **self._collection_schemas}
        except Exception:
            pass
        try:
            existing_doc = existing_store.get("graph_views")
            if existing_doc:
                existing_graph_views = existing_doc.get("value", {}).get("views", {}) or {}
        except Exception:
            pass
        try:
            existing_doc = existing_store.get("etl_metadata")
            if existing_doc:
                existing_metadata = existing_doc.get("value", {}) or {}
                existing_resolver_metadata = existing_metadata.get("resolver_metadata", {}) or {}
                existing_registry_datasets = existing_metadata.get("registry_datasets", []) or []
        except Exception:
            pass

        metadata_store = self.get_metadata_store(truncate=False)

        db_metadata = self.get_metadata()
        metadata_store.insert({
            "_key": "database_metadata",
            "value": db_metadata.to_dict()
        }, overwrite=True)

        etl_metadata = self.get_etl_metadata()
        etl_metadata["resolver_metadata"] = self._merge_resolver_metadata(
            existing_resolver_metadata,
            etl_metadata.get("resolver_metadata", {}),
        )
        etl_metadata["registry_datasets"] = self._merge_registry_datasets(
            existing_registry_datasets,
            etl_metadata.get("registry_datasets", []),
        )
        metadata_store.insert({
            "_key": "etl_metadata",
            "value": etl_metadata
        }, overwrite=True)

        if self._collection_schemas:
            metadata_store.insert({
                "_key": "collection_schemas",
                "collections": self._collection_schemas
            }, overwrite=True)
            print(f"Wrote collection schemas for {len(self._collection_schemas)} collections")

        graph_views = self.get_graph_views_metadata(existing_graph_views)
        if graph_views:
            metadata_store.insert({
                "_key": "graph_views",
                "value": {
                    "views": graph_views
                }
            }, overwrite=True)
            print(f"Wrote graph views metadata for {len(graph_views)} views")

    def get_etl_metadata(self):
        git_info = get_git_metadata()
        resolver_source_yaml = getattr(self, "_resolver_source_yaml", None)
        resolver_fingerprints_by_type = getattr(self, "_resolver_fingerprints_by_type", {})
        etl_metadata = {
            "_key": f"etl_run_{datetime.now().isoformat()}",
            "type": "etl_run",
            "run_date": datetime.now().isoformat(),
            "source_yaml": resolver_source_yaml,
            "resolver_metadata": {
                "source_yaml": resolver_source_yaml,
                "by_type": resolver_fingerprints_by_type,
                "summary": resolver_fingerprint_summary(resolver_fingerprints_by_type),
            },
            "registry_datasets": getattr(self, "_registry_datasets", []),
            "curations": getattr(self, "_curation_metadata", []),
            "curation_application_reports": getattr(self, "_curation_application_reports", []),
            "runner": os.getenv("USER", "unknown"),
            "git_info": git_info,
            "hostname": socket.gethostname(),
            "platform": platform.system(),
            "platform_version": platform.version(),
            "python_version": platform.python_version(),
        }
        return etl_metadata

    @staticmethod
    def _merge_resolver_metadata(existing_metadata=None, current_metadata=None):
        existing_metadata = existing_metadata or {}
        current_metadata = current_metadata or {}

        merged_by_type = {
            **(existing_metadata.get("by_type") or {}),
            **(current_metadata.get("by_type") or {}),
        }
        source_yamls = []
        for metadata in (existing_metadata, current_metadata):
            source_yaml = metadata.get("source_yaml")
            if source_yaml and source_yaml not in source_yamls:
                source_yamls.append(source_yaml)

        return {
            "source_yaml": current_metadata.get("source_yaml") or existing_metadata.get("source_yaml"),
            "source_yamls": source_yamls,
            "by_type": merged_by_type,
            "summary": resolver_fingerprint_summary(merged_by_type),
        }

    @staticmethod
    def _merge_registry_datasets(existing_datasets=None, current_datasets=None):
        datasets_by_identity = {}
        for dataset in (existing_datasets or []) + (current_datasets or []):
            if not isinstance(dataset, dict):
                continue
            snapshot_id = dataset.get("snapshot_id")
            kind = dataset.get("kind") or "source_snapshot"
            if not snapshot_id:
                continue
            identity = (kind, snapshot_id)
            merged = datasets_by_identity.setdefault(
                identity,
                {**dataset, "kind": kind, "usages": []},
            )
            usages = set(merged.get("usages") or [])
            usages.update(dataset.get("usages") or [])
            merged.update({key: value for key, value in dataset.items() if value is not None})
            merged["usages"] = sorted(usages)
        return sorted(
            datasets_by_identity.values(),
            key=lambda item: (item["kind"], item["snapshot_id"]),
        )

    def get_graph_views_metadata(self, existing_graph_views=None):
        existing_graph_views = existing_graph_views or {}
        merged_graph_views = dict(existing_graph_views)

        for graph_view in self._graph_views:
            if not graph_view or "id" not in graph_view:
                continue
            view_id = graph_view["id"]
            normalized_graph_view = dict(graph_view)
            if self._graph_view_source_yaml and "defined_in_yaml" not in normalized_graph_view:
                normalized_graph_view["defined_in_yaml"] = self._graph_view_source_yaml
            merged_graph_views[view_id] = normalized_graph_view

        return merged_graph_views


    def _get_document_existence_map(self, db, doc_ids, check_batch_size: int = 10000):
        existence_map = {}
        unique_doc_ids = list(dict.fromkeys(doc_ids))

        for start in range(0, len(unique_doc_ids), check_batch_size):
            batch_doc_ids = unique_doc_ids[start:start + check_batch_size]
            batch_results = list(db.aql.execute(
                """
                FOR doc_id IN @doc_ids
                    LET exists = DOCUMENT(doc_id) != null
                    RETURN {id: doc_id, exists: exists}
                """,
                bind_vars={"doc_ids": batch_doc_ids},
                batch_size=max(100, min(len(batch_doc_ids), 1000)),
                max_runtime=600,
            ))
            existence_map.update({item['id']: item['exists'] for item in batch_results})

        return existence_map

    def clean_up_dangling_edges(self, batch_size: int = 250000, check_batch_size: int = 10000):
        db = self.get_db()
        graph = self.get_graph()
        for edge_collection in graph.edge_definitions():
            collection_name = edge_collection['edge_collection']
            print(f'cleaning up {collection_name}')

            total_deleted = 0
            last_key = ''

            while True:
                start_time = time.time()

                # Get a batch of edges first, then check them
                key_filter = f"FILTER e._key > '{last_key}'" if last_key else ""

                # Step 1: Get edge batch
                cursor = db.aql.execute(
                    f"""
                    FOR e IN `{collection_name}`
                        {key_filter}
                        SORT e._key
                        LIMIT {batch_size}
                        RETURN {{_key: e._key, _from: e._from, _to: e._to}}
                    """,
                    batch_size=max(100, min(batch_size, 1000)),
                    max_runtime=600,
                )

                edges = list(cursor)
                if not edges:
                    break

                last_key = edges[-1]['_key']

                # Step 2: Check which ones are dangling (batch the DOCUMENT calls)
                from_exists = self._get_document_existence_map(
                    db,
                    [e['_from'] for e in edges],
                    check_batch_size=check_batch_size,
                )
                to_exists = self._get_document_existence_map(
                    db,
                    [e['_to'] for e in edges],
                    check_batch_size=check_batch_size,
                )

                # Find dangling edges
                dangling_keys = []
                dangling_edges = []
                for edge in edges:
                    missing_from = not from_exists.get(edge['_from'], True)
                    missing_to = not to_exists.get(edge['_to'], True)
                    if missing_from or missing_to:
                        dangling_keys.append(edge['_key'])
                        dangling_edges.append({**edge, "missing_from": missing_from, "missing_to": missing_to})

                # Step 3: Delete dangling edges by key
                if dangling_keys:
                    run_id = getattr(self, "_build_run_id", None)
                    if run_id:
                        self._record_dangling_edges(run_id, collection_name, dangling_edges, "observed")
                    db.aql.execute(
                        f"""
                        FOR key IN @keys
                            REMOVE key IN `{collection_name}`
                        """,
                        bind_vars={"keys": dangling_keys},
                        batch_size=max(100, min(len(dangling_keys), 1000)),
                        max_runtime=600,
                    )
                    if run_id:
                        self._record_dangling_edges(run_id, collection_name, dangling_edges, "deleted")

                deleted_count = len(dangling_keys)
                total_deleted += deleted_count

                print(f"Processed {len(edges)} edges, deleted {deleted_count} dangling, "
                      f"total deleted: {total_deleted}, time: {time.time() - start_time:.1f}s")

                if len(edges) < batch_size:
                    break

            print(f"Completed {collection_name}: {total_deleted} total edges deleted")
