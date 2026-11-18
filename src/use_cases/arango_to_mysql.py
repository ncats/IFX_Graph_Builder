"""
Generic Arango → MySQL converter.

Reads collection_schemas metadata from an ArangoDB metadata_store,
dynamically generates MySQL tables via SQLAlchemy, and copies all data over.
Works with any Arango DB that was built by ArangoOutputAdapter.
"""
import importlib
import json
import re
from itertools import product

from sqlalchemy import MetaData, Table, Column, String, Text, Integer, Float, Boolean, ForeignKey, Index, LargeBinary
from sqlalchemy.dialects.mysql import LONGBLOB

from src.input_adapters.sql_adapter import MySqlAdapter
from src.infrastructure.object_storage import ObjectStorageCredentials, object_storage_from_credentials
from src.shared.arango_adapter import ArangoAdapter
from src.shared.db_credentials import DBCredentials


# Map schema type strings to SQLAlchemy column types
_SQL_TYPE_MAP = {
    "str": Text,
    "int": Integer,
    "float": Float,
    "bool": Boolean,
    "date": Text,       # stored as ISO strings in Arango
    "datetime": Text,   # stored as ISO strings in Arango
}

# Fields from the base Node class that don't need MySQL tables
_SKIP_FIELDS = {"sources", "xref", "provenance"}


def _camel_to_snake(name: str) -> str:
    """Convert CamelCase to snake_case."""
    s = re.sub(r'([A-Z]+)([A-Z][a-z])', r'\1_\2', name)
    s = re.sub(r'([a-z0-9])([A-Z])', r'\1_\2', s)
    return s.lower()


def _edge_table_name(schema: dict) -> str:
    """Derive a MySQL table name for an edge collection: 'from_to_to'."""
    from_name = _camel_to_snake(schema["from_collections"][0])
    to_name = _camel_to_snake(schema["to_collections"][0])
    return f"{from_name}_to_{to_name}"


def _is_workbook_artifact_schema(field_schema: dict) -> bool:
    if not isinstance(field_schema, dict) or field_schema.get("type") != "object":
        return False
    fields = field_schema.get("fields", {})
    return {"file_reference", "original_filename", "media_type"}.issubset(fields.keys())


class ArangoToMySqlConverter(ArangoAdapter):
    """Reads from ArangoDB, writes to MySQL using schema metadata."""

    def __init__(self, arango_credentials: DBCredentials, arango_db_name: str,
                 mysql_credentials: DBCredentials, mysql_db_name: str,
                 object_storage_credentials: ObjectStorageCredentials = None):
        super().__init__(credentials=arango_credentials, database_name=arango_db_name)
        self.mysql = MySqlAdapter(mysql_credentials)
        self.mysql_db_name = mysql_db_name
        self.object_storage = (
            object_storage_from_credentials(object_storage_credentials)
            if object_storage_credentials
            else None
        )
        self.sa_metadata = MetaData()
        # enum_class_path -> (Table, enum_class) for LabeledIntEnum lookup tables
        self._enum_lookup_tables: dict = {}

    def convert(self, batch_size: int = 10000):
        """Run the full conversion from Arango to MySQL."""
        schemas = self._read_schemas()
        if not schemas:
            raise RuntimeError("No collection_schemas found in metadata_store")

        # Plan data tables — identifies which edge collections are replaced by data table FKs
        data_implicit_edges, data_table_configs = self._plan_data_tables(schemas)

        # Pass 1: create all tables
        document_collections = {}
        edge_collections = {}

        for collection_name, schema in schemas.items():
            if schema["type"] == "document":
                table, child_tables, object_tables, dict_tables = self._create_document_table(collection_name, schema["fields"])
                document_collections[collection_name] = (table, child_tables, object_tables, dict_tables, schema)
            elif collection_name not in data_implicit_edges:
                table = self._create_edge_table(collection_name, schema)
                edge_collections[collection_name] = (table, schema)

        data_tables = {}
        for config in data_table_configs:
            table = self._create_data_table(config)
            data_tables[config["table_name"]] = (table, config)

        self.mysql.recreate_mysql_db(self.mysql_db_name)
        engine = self.mysql.get_engine()
        self.sa_metadata.create_all(engine)
        self._populate_enum_lookup_tables(engine)

        # Pass 2: copy documents and edges
        for collection_name, (table, child_tables, object_tables, dict_tables, schema) in document_collections.items():
            self._copy_document_collection(engine, collection_name, table, child_tables, object_tables, dict_tables, schema, batch_size)

        for collection_name, (table, schema) in edge_collections.items():
            self._copy_edge_collection(engine, collection_name, table, schema, batch_size)

        # Pass 3: melt parquet data into typed data tables
        if data_tables:
            self._melt_parquet_data(engine, data_tables)

        print("Conversion complete.")

    def _read_schemas(self) -> dict:
        """Read collection_schemas from Arango metadata_store."""
        db = self.get_db()
        if not db.has_collection(self.metadata_store_label):
            return {}
        store = db.collection(self.metadata_store_label)
        doc = store.get("collection_schemas")
        if doc is None:
            return {}
        return doc.get("collections", {})

    # --- Data table planning ---

    def _plan_data_tables(self, schemas: dict):
        """Identify data-implicit edges and plan typed data tables.

        For collections with file_reference (Dataset, StatsResult), outgoing edges
        are replaced by foreign keys in data tables. Analyte targets (Gene, Metabolite)
        become row FKs. Sample targets (RunBiosample) become column FKs.

        Returns (data_implicit_edges, data_table_configs).
        """
        edge_schemas = {k: v for k, v in schemas.items() if v["type"] == "edge"}
        file_ref_collections = {
            k for k, v in schemas.items()
            if v["type"] == "document" and "file_reference" in v.get("fields", {})
        }

        data_implicit_edges = set()
        data_table_configs = []

        for parent_coll in file_ref_collections:
            parent_table = _camel_to_snake(parent_coll)
            analyte_targets = []
            sample_target = None

            for edge_name, edge_schema in edge_schemas.items():
                if parent_coll not in edge_schema.get("from_collections", []):
                    continue

                target_coll = edge_schema["to_collections"][0]
                target_table = _camel_to_snake(target_coll)

                non_file_ref_count = self._non_file_ref_source_count(target_coll, edge_schemas, file_ref_collections)
                if non_file_ref_count > 1:
                    # Metadata entity (e.g. Person) — keep as a regular edge table
                    continue

                data_implicit_edges.add(edge_name)

                if non_file_ref_count == 1:
                    # Sample dimension (e.g. RunBiosample) — becomes column FK in data tables
                    sample_target = {"collection": target_coll, "table": target_table, "edge": edge_name}
                else:
                    # Pure analyte (e.g. Gene, Metabolite) — becomes row FK in data tables
                    analyte_targets.append({"collection": target_coll, "table": target_table, "edge": edge_name})

            for analyte in analyte_targets:
                table_name = f"{analyte['table']}_{parent_table}__data"
                data_table_configs.append({
                    "table_name": table_name,
                    "parent_collection": parent_coll,
                    "parent_table": parent_table,
                    "analyte": analyte,
                    "sample": sample_target,
                })

        return data_implicit_edges, data_table_configs

    @staticmethod
    def _non_file_ref_source_count(target_coll: str, edge_schemas: dict, file_ref_collections: set) -> int:
        """Count distinct non-file-reference collections that have edges pointing to target_coll.

        0 → pure analyte (Gene, Metabolite): only ever targeted by data collections.
        1 → sample dimension (RunBiosample): linked from exactly one non-data collection (Biosample).
        2+ → metadata entity (Person): referenced by multiple high-level collections (Project, Experiment).
             These should remain as regular edge tables, not be folded into data tables.
        """
        sources = set()
        for edge_schema in edge_schemas.values():
            if target_coll not in edge_schema.get("to_collections", []):
                continue
            from_coll = edge_schema["from_collections"][0]
            if from_coll not in file_ref_collections:
                sources.add(from_coll)
        return len(sources)

    # --- Enum lookup tables ---

    def _get_or_create_enum_lookup_table(self, enum_path: str) -> tuple:
        """Return (Table, enum_class) for a LabeledIntEnum, creating the SA Table if needed."""
        if enum_path in self._enum_lookup_tables:
            return self._enum_lookup_tables[enum_path]
        module_path, class_name = enum_path.rsplit(".", 1)
        cls = getattr(importlib.import_module(module_path), class_name)
        table_name = _camel_to_snake(class_name)
        table = Table(
            table_name, self.sa_metadata,
            Column("id", Integer, primary_key=True, autoincrement=False),
            Column("label", String(255), nullable=False),
        )
        self._enum_lookup_tables[enum_path] = (table, cls)
        return table, cls

    def _populate_enum_lookup_tables(self, engine):
        """Insert enum members into every lookup table created during schema planning."""
        for table, cls in self._enum_lookup_tables.values():
            rows = [{"id": m.value, "label": m.label} for m in cls]
            with engine.connect() as conn:
                conn.execute(table.insert(), rows)
                conn.commit()
            print(f"  {table.name}: {len(rows)} enum values")

    # --- Table creation ---

    def _create_document_table(self, name: str, fields: dict) -> tuple:
        """Create a SQLAlchemy Table for a document collection.

        Returns (table, child_tables, object_tables, dict_tables) where:
          child_tables:  {field_name: (Table, {sub_field_name: Table})}  — for list fields
          object_tables: {field_name: (Table, {sub_field_name: Table})}  — for object fields
          dict_tables:   {field_name: Table}                             — for dict fields
        """
        table_name = _camel_to_snake(name)
        columns = []
        child_tables = {}
        object_tables = {}
        dict_tables = {}

        for field_name, field_schema in fields.items():
            if field_name in _SKIP_FIELDS:
                continue
            if isinstance(field_schema, dict):
                field_type = field_schema.get("type")
                if field_type == "list":
                    child_table, child_gc_tables = self._create_child_table(table_name, field_name, field_schema)
                    child_tables[field_name] = (child_table, child_gc_tables)
                    continue
                elif field_type == "object":
                    obj_table, grandchild_tables = self._create_object_table(table_name, field_name, field_schema)
                    object_tables[field_name] = (obj_table, grandchild_tables)
                    continue
                elif field_type == "dict":
                    dict_table = self._create_dict_table(table_name, field_name, field_schema)
                    dict_tables[field_name] = dict_table
                    continue
            col_type = _SQL_TYPE_MAP.get(field_schema, Text)
            if field_name == "id":
                columns.append(Column("id", String(255), primary_key=True))
            else:
                columns.append(Column(field_name, col_type, nullable=True))

        table = Table(table_name, self.sa_metadata, *columns)
        return table, child_tables, object_tables, dict_tables

    def _create_object_table(self, parent_table_name: str, field_name: str, field_schema: dict) -> tuple:
        """Create a linked table for an embedded object field.

        Named {parent}__{field} (e.g. 'exposure__category', 'biosample__demographics').
        The parent FK is used as the primary key since this is always a 1:1 relationship,
        avoiding duplicate key errors from object ids that are only locally unique.
        Grandchild list tables reference back via that same FK column.

        Returns (table, grandchild_tables) where grandchild_tables is
        {sub_field_name: Table} for any list sub-fields.
        """
        obj_table_name = f"{parent_table_name}__{_camel_to_snake(field_name)}"
        parent_fk_col = f"{parent_table_name}_id"
        sub_fields = field_schema.get("fields", {})

        # Parent FK doubles as PK — uniqueness is guaranteed by the 1:1 relationship
        columns = [
            Column(parent_fk_col, String(255), ForeignKey(f"{parent_table_name}.id"), primary_key=True),
        ]

        grandchild_tables = {}

        for sub_name, sub_schema in sub_fields.items():
            if sub_name in _SKIP_FIELDS:
                continue
            if isinstance(sub_schema, dict):
                sub_type = sub_schema.get("type")
                if sub_type == "list":
                    # Grandchild table references this table via parent_fk_col (its PK)
                    gc_table, _ = self._create_child_table(obj_table_name, sub_name, sub_schema,
                                                           parent_pk_col=parent_fk_col)
                    grandchild_tables[sub_name] = gc_table
                else:
                    # Nested object-within-object: store as JSON blob
                    columns.append(Column(sub_name, Text, nullable=True))
            else:
                col_type = _SQL_TYPE_MAP.get(sub_schema, Text)
                columns.append(Column(sub_name, col_type, nullable=True))

        if _is_workbook_artifact_schema(field_schema):
            columns.append(Column("content_blob", LargeBinary().with_variant(LONGBLOB(), "mysql"), nullable=True))

        table = Table(obj_table_name, self.sa_metadata, *columns)
        return table, grandchild_tables

    @staticmethod
    def _list_item_object_fields(field_schema: dict) -> dict | None:
        """Return the object fields dict if this list schema describes a list of objects, else None.

        Handles two formats produced by ArangoOutputAdapter._type_hint_to_schema:
          - item_type == "object" with fields at the top level  (List[dataclass] case)
          - item_type == {"type": "object", "fields": {...}}      (defensive)
        """
        item_type = field_schema.get("item_type")
        if item_type == "object":
            return field_schema.get("fields", {})
        if isinstance(item_type, dict) and item_type.get("type") == "object":
            return item_type.get("fields", {})
        return None

    def _create_child_table(self, parent_name: str, field_name: str, field_schema: dict,
                            parent_pk_col: str = "id", parent_pk_type=None) -> tuple:
        """Create a child table for a list field.

        Returns (table, grandchild_tables) where grandchild_tables is
        {sub_field_name: Table} for any nested list sub-fields (list-of-object case).
        When grandchild tables exist, the child table gets an auto-increment 'id' PK
        so grandchild rows can reference it.
        """
        table_name = f"{parent_name}__{field_name}"
        fk = ForeignKey(f"{parent_name}.{parent_pk_col}")
        pk_col_type = parent_pk_type if parent_pk_type is not None else String(255)

        obj_fields = self._list_item_object_fields(field_schema)
        grandchild_tables = {}

        if obj_fields is not None:
            # Detect nested list sub-fields that need their own grandchild tables
            nested_list_subs = {
                k for k, v in obj_fields.items()
                if k not in _SKIP_FIELDS and isinstance(v, dict) and v.get("type") == "list"
            }

            columns = []
            if nested_list_subs:
                # Auto-increment PK so grandchild rows can reference this table
                columns.append(Column("id", Integer, primary_key=True))
            columns.append(Column("parent_id", pk_col_type, fk, nullable=False))

            for sub_name, sub_schema in obj_fields.items():
                if sub_name in _SKIP_FIELDS:
                    continue
                if sub_name in nested_list_subs:
                    gc_table, _ = self._create_child_table(table_name, sub_name, sub_schema,
                                                           parent_pk_col="id", parent_pk_type=Integer)
                    grandchild_tables[sub_name] = gc_table
                    continue
                col_type = _SQL_TYPE_MAP.get(sub_schema, Text) if isinstance(sub_schema, str) else Text
                columns.append(Column(sub_name, col_type, nullable=True))
            columns.append(Index(f"ix_{table_name}_parent", "parent_id"))
        else:
            item_type = field_schema.get("item_type")
            enum_path = field_schema.get("enum")
            if enum_path:
                # LabeledIntEnum: child FK column named {lookup_table}_id
                lookup_table, _ = self._get_or_create_enum_lookup_table(enum_path)
                value_col = Column(f"{lookup_table.name}_id", Integer, ForeignKey(f"{lookup_table.name}.id"), nullable=True)
            else:
                col_type = _SQL_TYPE_MAP.get(item_type, Text) if isinstance(item_type, str) else Text
                value_col = Column("value", col_type, nullable=True)
            columns = [
                Column("parent_id", pk_col_type, fk, nullable=False),
                value_col,
                Index(f"ix_{table_name}_parent", "parent_id"),
            ]

        return Table(table_name, self.sa_metadata, *columns), grandchild_tables

    def _create_dict_table(self, parent_name: str, field_name: str, field_schema: dict) -> Table:
        """Create a key-value child table for a Dict field.

        Named {parent}__{field} with columns (parent_id FK, key, value).
        """
        table_name = f"{parent_name}__{field_name}"
        value_schema = field_schema.get("value_type")
        val_type = _SQL_TYPE_MAP.get(value_schema, Text) if isinstance(value_schema, str) else Text
        return Table(
            table_name, self.sa_metadata,
            Column("parent_id", String(255), ForeignKey(f"{parent_name}.id"), nullable=False),
            Column("key", String(255), nullable=False),
            Column("value", val_type, nullable=True),
            Index(f"ix_{table_name}_parent", "parent_id"),
        )

    def _create_edge_table(self, name: str, schema: dict) -> Table:
        """Create a SQLAlchemy Table for an edge collection."""
        table_name = _edge_table_name(schema)
        from_table = _camel_to_snake(schema["from_collections"][0])
        to_table = _camel_to_snake(schema["to_collections"][0])
        columns = [
            Column("from_id", String(255), ForeignKey(f"{from_table}.id"), nullable=False),
            Column("to_id", String(255), ForeignKey(f"{to_table}.id"), nullable=False),
        ]

        for field_name, field_schema in schema.get("fields", {}).items():
            if field_name in _SKIP_FIELDS:
                continue
            if isinstance(field_schema, dict):
                if (field_schema.get("type") == "list"
                        and isinstance(field_schema.get("item_type"), str)
                        and field_schema.get("item_type") != "object"):
                    col_name = field_name[:-1] if field_name.endswith("s") else field_name
                    col_type = _SQL_TYPE_MAP.get(field_schema["item_type"], Text)
                    columns.append(Column(col_name, col_type, nullable=True))
                continue
            col_type = _SQL_TYPE_MAP.get(field_schema, Text)
            columns.append(Column(field_name, col_type, nullable=True))

        table = Table(table_name, self.sa_metadata,
                      *columns,
                      Index(f"ix_{table_name}_from", "from_id"),
                      Index(f"ix_{table_name}_to", "to_id"))
        return table

    def _create_data_table(self, config: dict) -> Table:
        """Create a typed data table for melted parquet data with proper FKs."""
        table_name = config["table_name"]
        parent_table = config["parent_table"]
        analyte = config["analyte"]
        sample = config["sample"]

        columns = [
            Column(f"{parent_table}_id", String(255), ForeignKey(f"{parent_table}.id"), nullable=False),
            Column(f"{analyte['table']}_id", String(255), ForeignKey(f"{analyte['table']}.id"), nullable=False),
        ]

        if sample:
            columns.append(
                Column(f"{sample['table']}_id", String(255), ForeignKey(f"{sample['table']}.id"), nullable=False))
        else:
            columns.append(Column("column_name", String(255), nullable=False))

        columns.append(Column("value", Float, nullable=True))

        indexes = [
            Index(f"ix_{table_name}_parent", f"{parent_table}_id"),
            Index(f"ix_{table_name}_analyte", f"{analyte['table']}_id"),
        ]
        if sample:
            indexes.append(Index(f"ix_{table_name}_sample", f"{sample['table']}_id"))

        return Table(table_name, self.sa_metadata, *columns, *indexes)

    # --- Data copying ---

    def _read_collection_paginated(self, collection_name: str, batch_size: int):
        """Read documents from an Arango collection in batches using cursor pagination."""
        db = self.get_db()
        offset = 0
        while True:
            cursor = db.aql.execute(f"""
                FOR doc IN `{collection_name}`
                    SORT doc._key
                    LIMIT {offset}, {batch_size}
                    RETURN doc
            """)
            docs = list(cursor)
            if not docs:
                break
            yield docs
            offset += batch_size
            if len(docs) < batch_size:
                break

    def _copy_document_collection(self, engine, collection_name: str, table: Table,
                                  child_tables: dict, object_tables: dict, dict_tables: dict,
                                  schema: dict, batch_size: int):
        """Copy a document collection from Arango to MySQL."""
        fields = {k: v for k, v in schema["fields"].items() if k not in _SKIP_FIELDS}
        list_fields = {k for k, v in fields.items() if isinstance(v, dict) and v.get("type") == "list"}
        object_fields = set(object_tables.keys())
        dict_fields = set(dict_tables.keys())
        scalar_fields = set(fields.keys()) - list_fields - object_fields - dict_fields
        total = 0

        # Separate list fields: those whose child table has nested grandchild tables vs flat
        child_fields_with_gc = {f for f in list_fields if f in child_tables and child_tables[f][1]}
        child_fields_flat = list_fields - child_fields_with_gc

        for docs in self._read_collection_paginated(collection_name, batch_size):
            rows = []
            child_rows_flat = {field: [] for field in child_fields_flat}
            # Each entry is (child_row_dict, {sub_name: [scalar_values]})
            child_rows_with_gc = {field: [] for field in child_fields_with_gc}
            object_rows = {field: [] for field in object_fields}
            dict_rows = {field: [] for field in dict_fields}
            grandchild_rows = {
                field: {sub: [] for sub in grandchild_tables}
                for field, (_, grandchild_tables) in object_tables.items()
            }

            parent_table_name = _camel_to_snake(collection_name)

            for doc in docs:
                row = {}
                for field in scalar_fields:
                    val = doc.get(field)
                    if isinstance(val, (dict, list)):
                        val = json.dumps(val)
                    row[field] = val
                rows.append(row)

                doc_id = doc.get("id")

                # List fields
                for list_field in list_fields:
                    values = doc.get(list_field)
                    if not values or not isinstance(values, list):
                        continue
                    obj_fields_map = self._list_item_object_fields(fields[list_field])
                    if obj_fields_map is not None:
                        if list_field in child_fields_with_gc:
                            _, child_gc_tables = child_tables[list_field]
                            for item in values:
                                if isinstance(item, dict):
                                    child_row = {"parent_id": doc_id}
                                    gc_data = {sub_name: [] for sub_name in child_gc_tables}
                                    for k, v in item.items():
                                        if k in _SKIP_FIELDS:
                                            continue
                                        if k in child_gc_tables:
                                            if isinstance(v, list):
                                                gc_data[k] = v
                                            # else: null/empty — no column for this in child_table
                                        elif isinstance(v, (dict, list)):
                                            child_row[k] = json.dumps(v)
                                        else:
                                            child_row[k] = v
                                    child_rows_with_gc[list_field].append((child_row, gc_data))
                        else:
                            for item in values:
                                if isinstance(item, dict):
                                    child_rows_flat[list_field].append({"parent_id": doc_id, **{
                                        k: json.dumps(v) if isinstance(v, (dict, list)) else v
                                        for k, v in item.items() if k not in _SKIP_FIELDS
                                    }})
                    else:
                        enum_path = fields[list_field].get("enum") if isinstance(fields[list_field], dict) else None
                        if enum_path and enum_path in self._enum_lookup_tables:
                            lookup_table, enum_cls = self._enum_lookup_tables[enum_path]
                            col_name = f"{lookup_table.name}_id"
                            label_to_id = {m.label: m.value for m in enum_cls}
                            for item in values:
                                int_id = label_to_id.get(item)
                                if int_id is not None:
                                    child_rows_flat[list_field].append({"parent_id": doc_id, col_name: int_id})
                        else:
                            for item in values:
                                child_rows_flat[list_field].append({"parent_id": doc_id, "value": item})

                # Object fields: extract nested dicts into linked tables
                for obj_field in object_fields:
                    nested = doc.get(obj_field)
                    if not nested or not isinstance(nested, dict):
                        continue

                    obj_schema = fields[obj_field]
                    parent_fk_col = f"{parent_table_name}_id"

                    # PK of the object table is parent_fk_col = doc_id
                    obj_row = {parent_fk_col: doc_id}

                    for sub_name, sub_schema in obj_schema.get("fields", {}).items():
                        if sub_name in _SKIP_FIELDS:
                            continue
                        val = nested.get(sub_name)
                        if isinstance(sub_schema, dict):
                            sub_type = sub_schema.get("type")
                            if sub_type == "list":
                                # Populate grandchild rows; parent_id references doc_id
                                # (the object table's PK is the parent FK, not an object id)
                                list_vals = val if isinstance(val, list) else []
                                item_obj_fields = self._list_item_object_fields(sub_schema)
                                for item in list_vals:
                                    if item_obj_fields is not None:
                                        if isinstance(item, dict):
                                            gc_row = {"parent_id": doc_id}
                                            for k, v in item.items():
                                                if k not in _SKIP_FIELDS:
                                                    gc_row[k] = json.dumps(v) if isinstance(v, (dict, list)) else v
                                            grandchild_rows[obj_field][sub_name].append(gc_row)
                                    else:
                                        grandchild_rows[obj_field][sub_name].append(
                                            {"parent_id": doc_id, "value": item}
                                        )
                                continue
                            else:
                                val = json.dumps(val) if val is not None else None
                        elif isinstance(val, (dict, list)):
                            val = json.dumps(val)
                        obj_row[sub_name] = val

                    if _is_workbook_artifact_schema(obj_schema):
                        file_ref = nested.get("file_reference")
                        if file_ref:
                            obj_row["content_blob"] = self._get_s3_buffer(file_ref).read()

                    object_rows[obj_field].append(obj_row)

                # Dict fields: extract key-value pairs
                for dict_field in dict_fields:
                    val = doc.get(dict_field)
                    if val and isinstance(val, dict):
                        for k, v in val.items():
                            if isinstance(v, (dict, list)):
                                v = json.dumps(v)
                            dict_rows[dict_field].append({"parent_id": doc_id, "key": k, "value": v})

            with engine.connect() as conn:
                if rows:
                    conn.execute(table.insert(), rows)

                # Bulk insert flat child rows (no nested grandchild tables)
                for field, c_rows in child_rows_flat.items():
                    if c_rows and field in child_tables:
                        child_table, _ = child_tables[field]
                        conn.execute(child_table.insert(), c_rows)

                # Row-by-row insert for child rows that have nested grandchild tables,
                # capturing the auto-increment PK to use as parent_id for grandchild rows
                for field, rows_with_gc in child_rows_with_gc.items():
                    child_table, child_gc_tables = child_tables[field]
                    for child_row, gc_data in rows_with_gc:
                        result = conn.execute(child_table.insert().values(**child_row))
                        child_id = result.inserted_primary_key[0]
                        for sub_name, gc_values in gc_data.items():
                            if gc_values and sub_name in child_gc_tables:
                                gc_table = child_gc_tables[sub_name]
                                gc_rows = [{"parent_id": child_id, "value": v} for v in gc_values]
                                conn.execute(gc_table.insert(), gc_rows)

                for obj_field, o_rows in object_rows.items():
                    if o_rows:
                        obj_table, grandchild_tables = object_tables[obj_field]
                        conn.execute(obj_table.insert(), o_rows)
                    for sub_name, gc_rows in grandchild_rows[obj_field].items():
                        _, grandchild_tables = object_tables[obj_field]
                        if gc_rows and sub_name in grandchild_tables:
                            conn.execute(grandchild_tables[sub_name].insert(), gc_rows)

                for dict_field, d_rows in dict_rows.items():
                    if d_rows:
                        conn.execute(dict_tables[dict_field].insert(), d_rows)

                conn.commit()

            total += len(docs)

        print(f"  {collection_name}: {total} rows")

    def _copy_edge_collection(self, engine, collection_name: str, table: Table,
                              schema: dict, batch_size: int):
        """Copy an edge collection from Arango to MySQL."""
        fields = schema.get("fields", {})
        scalar_fields = {k for k, v in fields.items()
                         if isinstance(v, str) and k not in _SKIP_FIELDS}
        list_expand_fields = {
            k: k[:-1] if k.endswith("s") else k
            for k, v in fields.items()
            if isinstance(v, dict) and v.get("type") == "list"
            and isinstance(v.get("item_type"), str) and v.get("item_type") != "object"
            and k not in _SKIP_FIELDS
        }
        total = 0

        for docs in self._read_collection_paginated(collection_name, batch_size):
            rows = []
            for doc in docs:
                base = {
                    "from_id": doc.get("start_id"),
                    "to_id": doc.get("end_id"),
                }
                for field in scalar_fields:
                    base[field] = doc.get(field)

                if not list_expand_fields:
                    rows.append(base)
                    continue

                lists = [doc.get(field) or [None] for field in list_expand_fields]
                col_names = list(list_expand_fields.values())
                for combo in product(*lists):
                    rows.append({**base, **dict(zip(col_names, combo))})

            with engine.connect() as conn:
                if rows:
                    try:
                        conn.execute(table.insert(), rows)
                    except Exception as e:
                        conn.rollback()
                        print(f"  Batch insert failed for {collection_name}: {e}")
                        print(f"  Retrying row-by-row to identify bad rows...")
                        bad_rows = []
                        for row in rows:
                            try:
                                conn.execute(table.insert(), [row])
                            except Exception:
                                bad_rows.append(row)
                                print(bad_rows)
                                print('that row was bad, stopping here so you can figure it out')
                                raise Exception(bad_rows)
                        if bad_rows:
                            print(f"  {len(bad_rows)} bad rows in {collection_name}:")
                            for r in bad_rows[:20]:
                                print(f"    from_id={r.get('from_id')}  to_id={r.get('to_id')}")
                            if len(bad_rows) > 20:
                                print(f"    ... and {len(bad_rows) - 20} more")
                conn.commit()

            total += len(docs)

        print(f"  {collection_name}: {total} edges")

    # --- Object-storage / parquet helpers ---

    def _get_s3_buffer(self, file_ref: str):
        """Fetch an S3-backed object and return a BytesIO buffer."""
        import io

        if not file_ref.startswith("s3://"):
            raise ValueError(f"Expected s3:// URI, got: {file_ref}")
        if not self.object_storage:
            raise RuntimeError("object_storage_credentials are required to read parquet files")

        without_prefix = file_ref[len("s3://"):]
        bucket, key = without_prefix.split("/", 1)
        s3 = self.object_storage.client()
        response = s3.get_object(Bucket=bucket, Key=key)
        return io.BytesIO(response["Body"].read())

    def _get_parquet_buffer(self, file_ref: str):
        """Fetch a parquet file from object storage and return a BytesIO buffer."""
        return self._get_s3_buffer(file_ref)

    # --- Parquet data melting ---

    def _melt_parquet_data(self, engine, data_tables: dict):
        """Melt parquet files into typed data tables with proper FKs."""
        import pyarrow.parquet as pq
        import pandas as pd

        db = self.get_db()

        # Group configs by parent collection
        configs_by_parent = {}
        for table_name, (table, config) in data_tables.items():
            parent = config["parent_collection"]
            configs_by_parent.setdefault(parent, []).append((table, config))

        for parent_coll, configs in configs_by_parent.items():
            # Pre-fetch which documents belong to which analyte type
            # by querying each analyte edge collection for distinct start_ids
            doc_to_config = {}
            for table, config in configs:
                edge_name = config["analyte"]["edge"]
                if not db.has_collection(edge_name):
                    continue
                cursor = db.aql.execute(f"FOR e IN `{edge_name}` RETURN DISTINCT e.start_id")
                for doc_id in cursor:
                    doc_to_config[doc_id] = (table, config)

            # Read all documents with file_reference
            cursor = db.aql.execute(f"""
                FOR doc IN `{parent_coll}`
                    FILTER doc.file_reference != null
                    RETURN {{id: doc.id, file_reference: doc.file_reference}}
            """)

            for doc in cursor:
                doc_id = doc["id"]
                file_ref = doc["file_reference"]

                if doc_id not in doc_to_config:
                    print(f"  Warning: no analyte edges found for {doc_id}, skipping")
                    continue

                table, config = doc_to_config[doc_id]
                parent_table = config["parent_table"]
                analyte_table = config["analyte"]["table"]
                sample = config["sample"]

                try:
                    buf = self._get_parquet_buffer(file_ref)
                except Exception as e:
                    print(f"  Warning: could not fetch parquet for {doc_id}: {e}")
                    continue

                df = pq.read_table(buf).to_pandas()
                index_name = df.index.name or "index"

                melted = df.reset_index().melt(
                    id_vars=[index_name],
                    var_name="col_id",
                    value_name="value"
                )

                rows = []
                for _, r in melted.iterrows():
                    row = {
                        f"{parent_table}_id": doc_id,
                        f"{analyte_table}_id": str(r[index_name]),
                        "value": float(r["value"]) if pd.notna(r["value"]) else None,
                    }
                    if sample:
                        row[f"{sample['table']}_id"] = str(r["col_id"])
                    else:
                        row["column_name"] = str(r["col_id"])
                    rows.append(row)

                if rows:
                    with engine.connect() as conn:
                        for i in range(0, len(rows), 10000):
                            conn.execute(table.insert(), rows[i:i + 10000])
                        conn.commit()
                    print(f"  {config['table_name']}: {len(rows)} values from {doc_id}")
