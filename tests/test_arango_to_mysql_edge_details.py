from sqlalchemy import Column, MetaData, String, Table, create_engine, select
import pytest

from src.use_cases.arango_to_mysql import ArangoToMySqlConverter


EDGE_SCHEMA = {
    "type": "edge",
    "from_collections": ["Experiment"],
    "to_collections": ["MeasuredGene"],
    "fields": {
        "details": {
            "type": "list", "item_type": "object",
            "fields": {
                "source_feature_id": "str",
                "source_symbol": "str",
                "run_biosample_ids": {"type": "list", "item_type": "str"},
            },
        },
    },
}


def test_edge_object_details_survive_mysql_conversion():
    converter = object.__new__(ArangoToMySqlConverter)
    converter.sa_metadata = MetaData()
    Table("experiment", converter.sa_metadata, Column("id", String(255), primary_key=True))
    Table("measured_gene", converter.sa_metadata, Column("id", String(255), primary_key=True))
    edge_table = converter._create_edge_table("ExperimentMeasuredGeneEdge", EDGE_SCHEMA)
    detail_table, grandchild_tables = converter._edge_detail_tables["ExperimentMeasuredGeneEdge"]["details"]
    runs_table = grandchild_tables["run_biosample_ids"]
    assert "details" not in edge_table.c
    assert all(len(index.name) <= 64 for index in runs_table.indexes)
    engine = create_engine("sqlite:///:memory:")
    converter.sa_metadata.create_all(engine)

    details = [{
        "source_feature_id": "gene-S", "source_symbol": "CoV2_S",
        "run_biosample_ids": ["R1", "R2"],
    }]
    converter._read_collection_paginated = lambda name, size: [[{
        "_key": "edge-1", "start_id": "E1", "end_id": "F1", "details": details,
    }]]
    converter._copy_edge_collection(
        engine, "ExperimentMeasuredGeneEdge", edge_table, EDGE_SCHEMA, batch_size=1000,
    )

    with engine.connect() as conn:
        row = conn.execute(select(edge_table)).mappings().one()
        detail_row = conn.execute(select(detail_table)).mappings().one()
        run_rows = conn.execute(select(runs_table).order_by(runs_table.c.ordinal)).mappings().all()
    assert row["from_id"] == "E1"
    assert row["to_id"] == "F1"
    assert row["id"] == "edge-1"
    assert detail_row["parent_id"] == "edge-1"
    assert detail_row["source_feature_id"] == "gene-S"
    assert detail_row["source_symbol"] == "CoV2_S"
    assert detail_row["ordinal"] == 0
    assert [run["value"] for run in run_rows] == ["R1", "R2"]
    assert all(run["parent_id"] == detail_row["id"] for run in run_rows)


def test_edge_object_details_reject_malformed_payload():
    converter = object.__new__(ArangoToMySqlConverter)
    converter.sa_metadata = MetaData()
    Table("experiment", converter.sa_metadata, Column("id", String(255), primary_key=True))
    Table("measured_gene", converter.sa_metadata, Column("id", String(255), primary_key=True))
    edge_table = converter._create_edge_table("ExperimentMeasuredGeneEdge", EDGE_SCHEMA)
    engine = create_engine("sqlite:///:memory:")
    converter.sa_metadata.create_all(engine)
    converter._read_collection_paginated = lambda name, size: [[{
        "_key": "edge-1", "start_id": "E1", "end_id": "F1", "details": "not a list",
    }]]

    with pytest.raises(ValueError, match="edge-1.*details.*list of objects"):
        converter._copy_edge_collection(
            engine, "ExperimentMeasuredGeneEdge", edge_table, EDGE_SCHEMA, batch_size=1000,
        )
