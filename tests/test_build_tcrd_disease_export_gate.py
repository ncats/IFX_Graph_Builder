import pytest

from src.output_adapters import mysql_output_adapter


def test_shared_tcrd_output_gate_stops_truncation_and_resume(monkeypatch):
    output = object.__new__(mysql_output_adapter.TCRDOutputAdapter)
    output.source_graph_credentials = object()
    output.source_graph_database = "pharos"
    output.database_name = "pharos400"
    graph = type("Graph", (), {"has_collection": lambda self, name: True})()
    queries = []

    class FakeArangoAdapter:
        def __init__(self, credentials, database_name):
            assert credentials is output.source_graph_credentials
            assert database_name == "pharos"

        def get_db(self):
            return graph

        def runQuery(self, query):
            queries.append(query)
            return [True]

    monkeypatch.setattr(mysql_output_adapter, "ArangoAdapter", FakeArangoAdapter)
    monkeypatch.setattr(
        mysql_output_adapter.MySQLOutputAdapter,
        "create_or_truncate_datastore",
        lambda *_args, **_kwargs: pytest.fail("MySQL truncation was reached"),
    )

    with pytest.raises(RuntimeError, match="No MySQL tables were truncated"):
        output.create_or_truncate_datastore()
    with pytest.raises(RuntimeError, match="gene-projected"):
        output.do_pre_processing()
    assert len(queries) == 2
    assert all("gene_assessments" in query for query in queries)
