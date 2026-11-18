import io
import json

from sqlalchemy import MetaData, create_engine, select

from src.use_cases.arango_to_mysql import ArangoToMySqlConverter


WORKBOOK_SCHEMA = {
    "type": "object",
    "fields": {
        "file_reference": "str",
        "original_filename": "str",
        "media_type": "str",
    },
}


def _make_converter():
    converter = object.__new__(ArangoToMySqlConverter)
    converter.sa_metadata = MetaData()
    converter.object_storage = object()
    return converter


def test_create_object_table_adds_content_blob_for_workbook_artifact():
    converter = _make_converter()

    table, _ = converter._create_object_table("project", "workbook", WORKBOOK_SCHEMA)

    assert "content_blob" in table.c


def test_copy_document_collection_materializes_workbook_blob():
    converter = _make_converter()
    schema = {
        "fields": {
            "id": "str",
            "name": "str",
            "workbook": WORKBOOK_SCHEMA,
        }
    }
    table, child_tables, object_tables, dict_tables = converter._create_document_table("Project", schema["fields"])
    workbook_table, _ = object_tables["workbook"]

    engine = create_engine("sqlite:///:memory:")
    converter.sa_metadata.create_all(engine)
    converter._read_collection_paginated = lambda collection_name, batch_size: [[
        {
            "id": "PROJ001",
            "name": "Project 1",
            "workbook": {
                "file_reference": "s3://aws-ifx-registry/test_pounce/workbooks/PROJ001/test.xlsx",
                "original_filename": "test.xlsx",
                "media_type": "application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            },
        }
    ]]
    converter._get_s3_buffer = lambda file_ref: io.BytesIO(b"fake workbook bytes")

    converter._copy_document_collection(
        engine,
        "Project",
        table,
        child_tables,
        object_tables,
        dict_tables,
        schema,
        batch_size=1000,
    )

    with engine.connect() as conn:
        row = conn.execute(select(workbook_table)).mappings().one()

    assert row["project_id"] == "PROJ001"
    assert row["original_filename"] == "test.xlsx"
    assert row["content_blob"] == b"fake workbook bytes"


def test_copy_document_collection_preserves_structured_dict_values_as_json():
    converter = _make_converter()
    schema = {
        "fields": {
            "id": "str",
            "xref_match_types": {
                "type": "dict",
                "key_type": "str",
                "value_type": {"type": "list", "item_type": "str"},
            },
            "calculated_properties": {
                "type": "dict",
                "key_type": "str",
                "value_type": "float",
            },
        }
    }
    table, child_tables, object_tables, dict_tables = converter._create_document_table("Gene", schema["fields"])
    engine = create_engine("sqlite:///:memory:")
    converter.sa_metadata.create_all(engine)
    converter._read_collection_paginated = lambda collection_name, batch_size: [[{
        "id": "GENE:1",
        "xref_match_types": {"exact": ["Ensembl:ENSG1", "NCBIGene:1"]},
        "calculated_properties": {"score": 0.75},
    }]]

    converter._copy_document_collection(
        engine, "Gene", table, child_tables, object_tables, dict_tables, schema, batch_size=1000,
    )

    with engine.connect() as conn:
        match_row = conn.execute(select(dict_tables["xref_match_types"])).mappings().one()
        property_row = conn.execute(select(dict_tables["calculated_properties"])).mappings().one()

    assert match_row["parent_id"] == "GENE:1"
    assert match_row["key"] == "exact"
    assert json.loads(match_row["value"]) == ["Ensembl:ENSG1", "NCBIGene:1"]
    assert property_row["value"] == 0.75
