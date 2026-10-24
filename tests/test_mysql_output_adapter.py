from sqlalchemy.exc import IntegrityError, OperationalError
from datetime import date, datetime
from pathlib import Path

import pytest

from src.interfaces.metadata import CollectionMetadata, DatabaseMetadata
from src.interfaces.resolver_metadata import resolver_fingerprints_by_type
from src.models.datasource_version_info import DataSourceDetails
from src.output_adapters.mysql_output_adapter import MySQLOutputAdapter, TCRDOutputAdapter
from src.models.registry_dataset import RegistryDataset, RegistryDatasetKind
from src.shared.db_credentials import DBCredentials
from src.shared.sqlalchemy_tables.pharos_tables_new import TDL_info
from src.shared.sqlalchemy_tables.test_tables import AutoIncNode


def test_serialize_rows_fills_missing_optional_columns_with_none():
    rows = MySQLOutputAdapter._serialize_rows([
        TDL_info(itype="Ab Count", protein_id=1, integer_value=75),
        TDL_info(itype="JensenLab PubMed Score", protein_id=1, number_value=25.269722),
    ])

    assert rows == [
        {
            "itype": "Ab Count",
            "protein_id": 1,
            "integer_value": 75,
            "number_value": None,
        },
        {
            "itype": "JensenLab PubMed Score",
            "protein_id": 1,
            "integer_value": None,
            "number_value": 25.269722,
        },
    ]


def test_mysql_output_adapter_selects_database_on_init():
    credentials = DBCredentials(
        url="localhost",
        user="tester",
        password="secret",
        schema=None,
    )

    adapter = MySQLOutputAdapter(credentials=credentials, database_name="pharos400", truncate_tables=False)

    assert adapter.credentials.schema == "pharos400"
    assert adapter.connection_string.endswith("/pharos400")


def test_mysql_output_adapter_store_does_not_use_insert_ignore():
    credentials = DBCredentials(
        url="localhost",
        user="tester",
        password="secret",
        schema=None,
    )
    adapter = MySQLOutputAdapter(credentials=credentials, database_name="pharos400", truncate_tables=False)

    class FakeConverter:
        def get_object_converters(self, _obj_cls):
            return lambda _obj: AutoIncNode(identifier="pathway-1", value="demo")

    class FakeSession:
        def __init__(self):
            self.executed = []

        def execute(self, stmt, rows):
            self.executed.append((stmt, rows))

        def commit(self):
            pass

        def rollback(self):
            pass

        def close(self):
            pass

    fake_session = FakeSession()
    adapter.output_converter = FakeConverter()
    adapter.get_session = lambda: fake_session
    adapter.sort_and_convert_objects = lambda _objects, keep_nested_objects=True: {
        "AutoIncNode": ([{"dummy": "value"}], ["AutoIncNode"], False, None, None, object)
    }

    adapter.store(["unused"])

    assert len(fake_session.executed) == 1
    stmt, rows = fake_session.executed[0]
    assert stmt._prefixes == ()
    assert rows == [{"identifier": "pathway-1", "value": "demo"}]


def test_mysql_output_adapter_store_inserts_in_chunks():
    credentials = DBCredentials(
        url="localhost",
        user="tester",
        password="secret",
        schema=None,
    )
    adapter = MySQLOutputAdapter(credentials=credentials, database_name="pharos400", truncate_tables=False)
    adapter.conversion_chunk_size = 2
    adapter.insert_batch_size = 3

    class FakeConverter:
        def get_object_converters(self, _obj_cls):
            return lambda obj: [
                AutoIncNode(identifier=f"{obj['id']}-a", value="demo"),
                AutoIncNode(identifier=f"{obj['id']}-b", value="demo"),
            ]

    class FakeSession:
        def __init__(self):
            self.executed = []
            self.commit_calls = 0

        def execute(self, stmt, rows):
            self.executed.append((stmt, rows))

        def commit(self):
            self.commit_calls += 1

        def rollback(self):
            pass

        def close(self):
            pass

    fake_session = FakeSession()
    adapter.output_converter = FakeConverter()
    adapter.get_session = lambda: fake_session
    adapter.sort_and_convert_objects = lambda _objects, keep_nested_objects=True: {
        "AutoIncNode": (
            [{"id": "one"}, {"id": "two"}, {"id": "three"}],
            ["AutoIncNode"],
            False,
            None,
            None,
            object,
        )
    }

    adapter.store(["unused"])

    assert [len(rows) for _, rows in fake_session.executed] == [3, 1, 2]
    assert fake_session.commit_calls == 3


def test_mysql_output_adapter_retries_disconnect_with_smaller_chunks():
    credentials = DBCredentials(
        url="localhost",
        user="tester",
        password="secret",
        schema=None,
    )
    adapter = MySQLOutputAdapter(
        credentials=credentials,
        database_name="pharos400",
        truncate_tables=False,
    )
    adapter.conversion_chunk_size = 10
    adapter.insert_batch_size = 10
    adapter.min_insert_batch_size = 1

    class FakeConverter:
        def get_object_converters(self, _obj_cls):
            return lambda obj: AutoIncNode(identifier=obj["id"], value="demo")

    class FakeSession:
        def __init__(self, controller):
            self.controller = controller

        def execute(self, stmt, rows):
            self.controller["executed"].append((stmt, rows))
            if len(rows) > 2:
                raise OperationalError(
                    "stmt",
                    rows,
                    Exception(2013, "Lost connection to MySQL server during query ([Errno 110] Connection timed out)"),
                )

        def commit(self):
            self.controller["commit_calls"] += 1

        def rollback(self):
            self.controller["rollback_calls"] += 1

        def close(self):
            self.controller["close_calls"] += 1

    controller = {
        "executed": [],
        "commit_calls": 0,
        "rollback_calls": 0,
        "close_calls": 0,
    }

    adapter.output_converter = FakeConverter()
    adapter.get_session = lambda: FakeSession(controller)
    adapter.sort_and_convert_objects = lambda _objects, keep_nested_objects=True: {
        "AutoIncNode": (
            [{"id": "one"}, {"id": "two"}, {"id": "three"}, {"id": "four"}, {"id": "five"}],
            ["AutoIncNode"],
            False,
            None,
            None,
            object,
        )
    }

    adapter.store(["unused"])

    assert [len(rows) for _, rows in controller["executed"]] == [5, 2, 2, 1]
    assert controller["commit_calls"] == 3
    assert controller["rollback_calls"] == 1


def test_tcrd_output_adapter_preloads_mappings_in_pre_processing():
    credentials = DBCredentials(
        url="localhost",
        user="tester",
        password="secret",
        schema=None,
    )
    adapter = TCRDOutputAdapter(credentials=credentials, database_name="pharos400", truncate_tables=False)

    class FakeConverter:
        def __init__(self):
            self.sessions = []

        def preload_id_mappings(self, session):
            self.sessions.append(session)

    class FakeSession:
        def __init__(self):
            self.closed = False

        def close(self):
            self.closed = True

    fake_converter = FakeConverter()
    fake_session = FakeSession()
    adapter.output_converter = fake_converter
    adapter.get_session = lambda: fake_session

    adapter.do_pre_processing()

    assert fake_converter.sessions == [fake_session]
    assert fake_session.closed is True


def _resolver_dataset(version="deps-test"):
    return RegistryDataset(
        kind=RegistryDatasetKind.SOURCE,
        source="target_graph",
        dataset="disease_ids",
        version=version,
        version_date=None,
        download_date=date(2026, 6, 12),
        snapshot_id=f"target_graph:disease_ids:{version}",
        manifest_uri=f"s3://ifx-registry/sources/target_graph/disease_ids/{version}/manifest.yaml",
        manifest={"kind": "source_snapshot", "files": []},
        local_dir=Path(f"/tmp/ifx-registry-cache/sources/target_graph/disease_ids/{version}"),
    )


def _resolver_config(class_name="DiseaseIdResolver", version="deps-test", types=None):
    return [{
        "label": "disease_ids",
        "import": "./src/id_resolvers/disease_resolver.py",
        "class": class_name,
        "kwargs": {
            "data_source": _resolver_dataset(version=version),
            "types": types or ["Disease"],
            "multi_match_behavior": "All",
        },
    }]


def _adapter_with_resolver_metadata(expected_metadata, graph_metadata):
    adapter = TCRDOutputAdapter.__new__(TCRDOutputAdapter)
    adapter.source_graph_credentials = object()
    adapter.source_graph_database = "pharos"
    adapter._resolver_fingerprints_by_type = expected_metadata
    adapter._resolver_source_yaml = "./src/use_cases/pharos/tcrd.yaml"
    adapter._load_graph_etl_metadata = lambda: graph_metadata
    return adapter


def test_tcrd_output_adapter_accepts_matching_source_graph_resolver_metadata():
    metadata = resolver_fingerprints_by_type(_resolver_config())
    adapter = _adapter_with_resolver_metadata(metadata, {
        "source_yaml": "./src/use_cases/pharos/pharos.yaml",
        "resolver_metadata": {
            "source_yaml": "./src/use_cases/pharos/pharos.yaml",
            "by_type": metadata,
        },
    })

    adapter._validate_source_graph_resolver_metadata()


def test_tcrd_output_adapter_requires_source_graph_resolver_metadata():
    metadata = resolver_fingerprints_by_type(_resolver_config())
    adapter = _adapter_with_resolver_metadata(metadata, {
        "source_yaml": "./src/use_cases/pharos/pharos.yaml",
    })

    with pytest.raises(RuntimeError, match="resolver metadata is missing"):
        adapter._validate_source_graph_resolver_metadata()


def test_tcrd_output_adapter_rejects_mismatched_source_graph_resolver_metadata():
    expected = resolver_fingerprints_by_type(_resolver_config())
    actual = resolver_fingerprints_by_type(_resolver_config(version="deps-other"))
    adapter = _adapter_with_resolver_metadata(expected, {
        "source_yaml": "./src/use_cases/pharos/pharos.yaml",
        "resolver_metadata": {
            "source_yaml": "./src/use_cases/pharos/pharos.yaml",
            "by_type": actual,
        },
    })

    with pytest.raises(RuntimeError, match="mismatched_types=\\['Disease'\\]"):
        adapter._validate_source_graph_resolver_metadata()


def test_tcrd_output_adapter_rejects_fingerprint_drift():
    expected = resolver_fingerprints_by_type(_resolver_config())
    actual = resolver_fingerprints_by_type(_resolver_config())
    actual["Disease"] = {
        **actual["Disease"],
        "fingerprint": "different-expanded-dependency-fingerprint",
    }
    adapter = _adapter_with_resolver_metadata(expected, {
        "source_yaml": "./src/use_cases/pharos/pharos.yaml",
        "resolver_metadata": {
            "source_yaml": "./src/use_cases/pharos/pharos.yaml",
            "by_type": actual,
        },
    })

    with pytest.raises(RuntimeError, match="mismatched_types=\\['Disease'\\]"):
        adapter._validate_source_graph_resolver_metadata()


def test_tcrd_output_adapter_rejects_graph_resolver_with_broader_type_scope():
    expected = resolver_fingerprints_by_type(_resolver_config())
    actual = resolver_fingerprints_by_type(_resolver_config(types=["Ligand", "Disease"]))
    adapter = _adapter_with_resolver_metadata(expected, {
        "source_yaml": "./src/use_cases/pharos/pharos.yaml",
        "resolver_metadata": {
            "source_yaml": "./src/use_cases/pharos/pharos.yaml",
            "by_type": actual,
        },
    })

    with pytest.raises(RuntimeError, match="mismatched_types=\\['Disease'\\]"):
        adapter._validate_source_graph_resolver_metadata()


def test_tcrd_resolver_validation_precedes_mysql_recreation(monkeypatch):
    adapter = TCRDOutputAdapter.__new__(TCRDOutputAdapter)
    calls = []

    def reject_source_graph():
        calls.append("validate")
        raise RuntimeError("source graph resolver mismatch")

    def recreate_mysql(self, truncate_tables=None):
        calls.append("recreate")
        return True

    adapter._validate_source_graph_resolver_metadata = reject_source_graph
    monkeypatch.setattr(MySQLOutputAdapter, "create_or_truncate_datastore", recreate_mysql)

    with pytest.raises(RuntimeError, match="source graph resolver mismatch"):
        adapter.create_or_truncate_datastore()
    assert calls == ["validate"]


def test_tcrd_output_adapter_builds_transitive_closure_with_self_rows():
    closure = TCRDOutputAdapter._build_transitive_closure(
        node_ids=["A", "B", "C", "D"],
        direct_edges=[("B", "A"), ("C", "B")],
    )

    assert closure == [
        ("A", "A"),
        ("B", "A"),
        ("B", "B"),
        ("C", "A"),
        ("C", "B"),
        ("C", "C"),
        ("D", "D"),
    ]


def test_tcrd_output_adapter_builds_deduped_data_source_version_rows():
    metadata = DatabaseMetadata(collections=[
        CollectionMetadata(
            name="Protein",
            total_count=1,
            sources=[
                DataSourceDetails(
                    name="UniProt",
                    version="2026_01",
                    version_date=None,
                    download_date=None,
                ),
                DataSourceDetails(
                    name="GOA (GO)",
                    version=None,
                    version_date=None,
                    download_date=None,
                ),
            ],
        ),
        CollectionMetadata(
            name="ProteinGoTermEdge",
            total_count=1,
            sources=[
                DataSourceDetails(
                    name="GOA (GO)",
                    version=None,
                    version_date=None,
                    download_date=None,
                ),
                DataSourceDetails(
                    name="GOA (UniProt)",
                    version=None,
                    version_date=None,
                    download_date=None,
                ),
            ],
        ),
    ])

    rows = TCRDOutputAdapter._build_data_source_version_rows(metadata)

    assert rows == [
        {
            "data_source": "GOA (GO)",
            "version": None,
            "version_date": None,
            "download_date": None,
            "collections": "Protein, ProteinGoTermEdge",
        },
        {
            "data_source": "GOA (UniProt)",
            "version": None,
            "version_date": None,
            "download_date": None,
            "collections": "ProteinGoTermEdge",
        },
        {
            "data_source": "UniProt",
            "version": "2026_01",
            "version_date": None,
            "download_date": None,
            "collections": "Protein",
        },
    ]


class FakeExecuteSession:
    def __init__(self):
        self.statements = []
        self.commit_calls = 0
        self.rollback_calls = 0
        self.close_calls = 0

    def execute(self, stmt):
        self.statements.append(str(stmt))

    def commit(self):
        self.commit_calls += 1

    def rollback(self):
        self.rollback_calls += 1

    def close(self):
        self.close_calls += 1


def test_tcrd_output_adapter_populates_disease_summary_fields_with_idempotent_sql():
    session = FakeExecuteSession()

    TCRDOutputAdapter._populate_ncats_disease_summary_fields(session)

    assert len(session.statements) == 9
    assert "SET\n              `target_count` = 0" in session.statements[0]
    assert "`direct_target_count` = 0" in session.statements[0]
    assert "`maxTDL` = NULL" in session.statements[0]
    assert "COUNT(DISTINCT d.`protein_id`) AS `direct_target_count`" in session.statements[1]
    assert "DROP TEMPORARY TABLE IF EXISTS `ncats_disease_descendant_stage`" in session.statements[2]
    assert "CREATE TEMPORARY TABLE `ncats_disease_descendant_stage`" in session.statements[3]
    assert "FROM `ncats_disease`" in session.statements[4]
    assert "JOIN `ancestry_mondo` ancestry" in session.statements[5]
    assert "COUNT(DISTINCT d.`protein_id`) AS `target_count`" in session.statements[6]
    assert "DROP TEMPORARY TABLE IF EXISTS `ncats_disease_descendant_stage`" in session.statements[7]
    assert "SET n.`maxTDL` = CASE" in session.statements[8]


def test_tcrd_output_adapter_populates_disease_summary_fields_max_tdl_from_direct_targets():
    session = FakeExecuteSession()

    TCRDOutputAdapter._populate_ncats_disease_summary_fields(session)

    max_tdl_sql = session.statements[-1]
    assert "JOIN `ncats_d2da` x" in max_tdl_sql
    assert "JOIN `disease` d" in max_tdl_sql
    assert "JOIN `t2tc`" in max_tdl_sql
    assert "JOIN `target` tgt" in max_tdl_sql
    assert "WHERE x.`direct` = 1" in max_tdl_sql
    assert "WHEN tgt.`tdl` = 'Tclin' THEN 4" in max_tdl_sql
    assert "WHEN tgt.`tdl` = 'Tchem' THEN 3" in max_tdl_sql
    assert "WHEN tgt.`tdl` = 'Tbio' THEN 2" in max_tdl_sql


def test_tcrd_output_adapter_populates_ligand_summary_fields_with_idempotent_sql():
    session = FakeExecuteSession()

    TCRDOutputAdapter._populate_ncats_ligand_summary_fields(session)

    assert len(session.statements) == 2
    assert "SET\n              `actCnt` = 0" in session.statements[0]
    assert "`targetCount` = 0" in session.statements[0]
    assert "FROM `ncats_ligand_activity` a" in session.statements[1]
    assert "COUNT(a.`id`) AS `actCnt`" in session.statements[1]
    assert "COUNT(DISTINCT a.`target_id`) AS `targetCount`" in session.statements[1]


def test_tcrd_output_adapter_populates_typeahead_index_with_idempotent_sql():
    session = FakeExecuteSession()

    TCRDOutputAdapter._populate_typeahead_index(session)

    assert session.statements[0] == "DELETE FROM `ncats_typeahead_index`"
    assert len(session.statements) == 13
    assert all("`category`, `value`, `reference_id`" in stmt for stmt in session.statements[1:])
    assert any("'Targets'" in stmt for stmt in session.statements)
    assert any("'Diseases'" in stmt for stmt in session.statements)


def test_tcrd_output_adapter_builds_etl_run_rows():
    graph_etl = {
        "run_date": "2026-04-04T11:49:01.558897",
        "runner": "kelleherkj",
        "hostname": "ncatslnifxdvv01",
        "platform": "Linux",
        "platform_version": "#181-Ubuntu SMP",
        "python_version": "3.11.11",
        "git_info": {
            "git_commit": "abc123",
            "git_branch": "main",
            "git_tag": None,
        },
    }

    rows = TCRDOutputAdapter._build_etl_run_rows(graph_etl, "pharos400")

    assert len(rows) == 2
    assert rows[0] == {
        "stage": "graph_build",
        "database_name": "pharos400",
        "run_date": datetime.fromisoformat("2026-04-04T11:49:01.558897"),
        "git_commit": "abc123",
        "git_branch": "main",
        "git_tag": None,
        "runner": "kelleherkj",
        "hostname": "ncatslnifxdvv01",
        "platform": "Linux",
        "platform_version": "#181-Ubuntu SMP",
        "python_version": "3.11.11",
    }
    assert rows[1]["stage"] == "graph_to_mysql"
    assert rows[1]["database_name"] == "pharos400"
    assert isinstance(rows[1]["run_date"], datetime)


def test_mysql_output_adapter_detects_fk_integrity_error():
    exc = IntegrityError(
        "stmt",
        {},
        Exception("Cannot add or update a child row: a foreign key constraint fails"),
    )

    assert MySQLOutputAdapter._is_fk_integrity_error(exc) is True


def test_mysql_output_adapter_store_invokes_fk_diagnosis():
    credentials = DBCredentials(
        url="localhost",
        user="tester",
        password="secret",
        schema=None,
    )
    adapter = MySQLOutputAdapter(credentials=credentials, database_name="pharos400", truncate_tables=False)

    class FakeConverter:
        def get_object_converters(self, _obj_cls):
            return lambda _obj: AutoIncNode(identifier="pathway-1", value="demo")

    class FakeSession:
        def execute(self, stmt, rows):
            raise IntegrityError(
                "stmt",
                rows,
                Exception("Cannot add or update a child row: a foreign key constraint fails"),
            )

        def commit(self):
            pass

        def rollback(self):
            pass

        def close(self):
            pass

    calls = []

    def fake_diagnose(table_class, rows):
        calls.append((table_class, rows))
        raise RuntimeError("isolated bad row")

    adapter.output_converter = FakeConverter()
    adapter.get_session = lambda: FakeSession()
    adapter._diagnose_fk_batch_failure = fake_diagnose
    adapter.sort_and_convert_objects = lambda _objects, keep_nested_objects=True: {
        "AutoIncNode": ([{"dummy": "value"}], ["AutoIncNode"], False, None, None, object)
    }

    try:
        adapter.store(["unused"])
    except RuntimeError as exc:
        assert str(exc) == "isolated bad row"
    else:
        raise AssertionError("Expected store() to raise the diagnostic RuntimeError")

    assert len(calls) == 1
    table_class, rows = calls[0]
    assert table_class is AutoIncNode
    assert rows == [{"identifier": "pathway-1", "value": "demo"}]
