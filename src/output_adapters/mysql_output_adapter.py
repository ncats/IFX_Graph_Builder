from abc import ABC
from datetime import datetime
import os
import platform
import socket
from src.interfaces.metadata import DatabaseMetadata, get_git_metadata
from sqlalchemy import case, func, text
from sqlalchemy import inspect as sqlalchemy_inspect
from sqlalchemy.dialects.mysql import insert as mysql_insert
from sqlalchemy.exc import IntegrityError, OperationalError
from src.input_adapters.sql_adapter import MySqlAdapter
from src.interfaces.output_adapter import OutputAdapter
from src.interfaces.resolver_metadata import resolver_fingerprint_summary
from src.output_adapters.sql_converters.output_converter_base import SQLOutputConverter
from src.output_adapters.sql_converters.tcrd import TCRDOutputConverter
from src.output_adapters.sql_converters.test import TestSQLOutputConverter
from src.shared.arango_adapter import ArangoAdapter
from src.shared.db_credentials import DBCredentials
from src.shared.sqlalchemy_tables.pharos_tables_new import (
    AncestryDO,
    AncestryDTO,
    AncestryMONDO,
    AncestryUBERON,
    DataSourceVersion,
    DO,
    DOParent,
    DTO,
    DTOParent,
    ETLAdapterRun,
    ETLRun,
    Mondo,
    MondoParent,
    NcatsDisease,
    NcatsTypeaheadIndex,
    PantherClass,
    Tissue,
    TinxImportance,
    Uberon,
    UberonParent,
    WordCount,
)
from src.shared.sqlalchemy_tables.test_tables import Base as TestBase


TYPEAHEAD_INSERTS = [
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT DISTINCT 'Diseases', LEFT(`ncats_name`, 255), LEFT(`ncats_name`, 255)
    FROM `disease`
    WHERE `ncats_name` IS NOT NULL
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT DISTINCT 'Drugs', LEFT(`name`, 255), LEFT(`name`, 255)
    FROM `ncats_ligands`
    WHERE `isDrug` = 1
      AND `name` IS NOT NULL
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT
      'Family',
      CASE `fam`
        WHEN 'IC' THEN 'Ion Channel'
        WHEN 'TF; Epigenetic' THEN 'TF-Epigenetic'
        WHEN 'TF' THEN 'Transcription Factor'
        WHEN 'NR' THEN 'Nuclear Receptor'
        ELSE LEFT(`fam`, 255)
      END,
      NULL
    FROM `target`
    WHERE `fam` IS NOT NULL
    GROUP BY `fam`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT
      'Genes',
      `sym`,
      CASE WHEN COUNT(*) = 1 THEN MAX(`sym`) ELSE MAX(`uniprot`) END
    FROM `protein`
    WHERE `sym` IS NOT NULL
    GROUP BY `sym`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT DISTINCT CONCAT('GO ', `go_type`), LEFT(`go_term_text`, 255), NULL
    FROM `goa`
    WHERE `go_type` IN ('Function', 'Process')
      AND `go_term_text` IS NOT NULL
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT 'GWAS', `trait`, NULL
    FROM `tiga`
    WHERE `trait` IS NOT NULL
    GROUP BY `trait`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT 'IMPC Phenotype', `term_name`, NULL
    FROM `phenotype`
    WHERE `ptype` = 'IMPC'
      AND `term_name` IS NOT NULL
    GROUP BY `term_name`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT 'Interacting Virus', v.`name`, NULL
    FROM `viral_ppi` p
    JOIN `viral_protein` vp
      ON p.`viral_protein_id` = vp.`id`
    JOIN `virus` v
      ON vp.`virus_id` = v.`virusTaxid`
    WHERE p.`finalLR` >= 500
      AND v.`name` IS NOT NULL
    GROUP BY v.`virusTaxid`, v.`name`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT 'JAX/MGI Phenotype', `term_name`, NULL
    FROM `phenotype`
    WHERE `ptype` = 'JAX/MGI Human Ortholog Phenotype'
      AND `term_name` IS NOT NULL
    GROUP BY `term_name`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT DISTINCT CONCAT(`pwtype`, ' Pathway'), LEFT(`name`, 255), NULL
    FROM `pathway`
    WHERE `pwtype` IN ('WikiPathways', 'KEGG', 'Reactome')
      AND `name` IS NOT NULL
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT 'UniProt Keyword', `xtra`, NULL
    FROM `xref`
    WHERE `xtype` = 'UniProt Keyword'
      AND `xtra` IS NOT NULL
    GROUP BY `xtra`
    """,
    """
    INSERT INTO `ncats_typeahead_index` (`category`, `value`, `reference_id`)
    SELECT
      'Targets',
      t.`name`,
      COALESCE(MAX(p.`sym`), MAX(p.`uniprot`))
    FROM `target` t
    JOIN `t2tc` x
      ON x.`target_id` = t.`id`
    JOIN `protein` p
      ON x.`protein_id` = p.`id`
    WHERE t.`name` IS NOT NULL
    GROUP BY t.`name`
    """,
]


class MySQLOutputAdapter(OutputAdapter, MySqlAdapter, ABC):
    database_name: str
    truncate_tables: bool
    output_converter: SQLOutputConverter
    conversion_chunk_size: int = 500
    insert_batch_size: int = 100_000
    min_insert_batch_size: int = 100
    adapter_run_model = None

    def __init__(
        self,
        credentials: DBCredentials,
        database_name: str,
        truncate_tables: bool = True,
    ):
        self.database_name = database_name
        self.truncate_tables = truncate_tables
        self._current_run_id = None
        self._current_adapter_name = None
        self._current_adapter_stats = None
        OutputAdapter.__init__(self)
        MySqlAdapter.__init__(self, credentials)
        self.update_database(database_name)

    @staticmethod
    def _serialize_rows(converted_objects):
        raw_rows = []
        for obj in converted_objects:
            mapper = sqlalchemy_inspect(type(obj))
            attr_to_column = {
                attr.key: attr.columns[0].name
                for attr in mapper.column_attrs
                if len(attr.columns) == 1
            }
            raw_rows.append({
                attr_to_column.get(k, k): v
                for k, v in obj.__dict__.items()
                if k != '_sa_instance_state'
            })
        all_keys = set().union(*(row.keys() for row in raw_rows))
        return [{key: row.get(key) for key in all_keys} for row in raw_rows]

    @staticmethod
    def _is_fk_integrity_error(exc: Exception) -> bool:
        message = str(exc).lower()
        return isinstance(exc, IntegrityError) and (
            "foreign key constraint fails" in message
            or "cannot add or update a child row" in message
        )

    def _diagnose_fk_batch_failure(self, table_class, rows):
        stmt = mysql_insert(table_class.__table__)
        print(f"Batch insert failed for {table_class.__name__}; retrying row-by-row to isolate FK issue")
        bad_rows = []

        for index, row in enumerate(rows, start=1):
            temp_session = self.get_session()
            try:
                temp_session.execute(stmt, [row])
                temp_session.flush()
            except Exception as row_exc:
                temp_session.rollback()
                bad_rows.append((index, row, row_exc))
                print(f"Bad row {index} for {table_class.__name__}: {row}")
                print(f"Row error: {row_exc}")
                break
            finally:
                temp_session.rollback()
                temp_session.close()

        if bad_rows:
            index, row, row_exc = bad_rows[0]
            raise RuntimeError(
                f"Foreign key insert failed for {table_class.__name__} at row {index}: {row}"
            ) from row_exc
        raise RuntimeError(
            f"Foreign key insert failed for {table_class.__name__}, but row-by-row replay did not isolate a bad row"
        )

    @staticmethod
    def _chunked(items, chunk_size: int):
        for index in range(0, len(items), chunk_size):
            yield items[index:index + chunk_size]

    @staticmethod
    def _is_retryable_disconnect_error(exc: Exception) -> bool:
        if not isinstance(exc, OperationalError):
            return False

        message = str(exc).lower()
        if "lost connection to mysql server during query" in message:
            return True
        if "server has gone away" in message:
            return True
        if "connection timed out" in message:
            return True

        orig = getattr(exc, "orig", None)
        args = getattr(orig, "args", ())
        code = args[0] if args else None
        return code in {2006, 2013}

    def _execute_insert_chunk(self, stmt, table_class, rows, stats, retry_depth: int = 0):
        session = self.get_session()
        try:
            insert_start = datetime.now()
            session.execute(stmt, rows)
            session.commit()
            if stats is not None:
                stats["insert_seconds"] += (datetime.now() - insert_start).total_seconds()
                stats["inserted_row_count"] += len(rows)
            return
        except Exception as exc:
            session.rollback()
            if self._is_fk_integrity_error(exc):
                self._diagnose_fk_batch_failure(table_class, rows)
            if self._is_retryable_disconnect_error(exc) and len(rows) > self.min_insert_batch_size:
                next_chunk_size = max(self.min_insert_batch_size, len(rows) // 2)
                print(
                    f"Retryable MySQL insert failure for {table_class.__name__} "
                    f"with {len(rows)} rows; retrying in chunks of {next_chunk_size}"
                )
                for row_chunk in self._chunked(rows, next_chunk_size):
                    self._execute_insert_chunk(
                        stmt,
                        table_class,
                        row_chunk,
                        stats,
                        retry_depth=retry_depth + 1,
                    )
                return
            raise
        finally:
            session.close()

    @staticmethod
    def _new_adapter_stats() -> dict:
        return {
            "store_calls": 0,
            "input_object_count": 0,
            "converted_object_count": 0,
            "inserted_row_count": 0,
            "conversion_seconds": 0.0,
            "serialization_seconds": 0.0,
            "insert_seconds": 0.0,
            "total_store_seconds": 0.0,
        }

    def _supports_adapter_run_metadata(self) -> bool:
        return self.adapter_run_model is not None

    def _get_or_create_adapter_run_row(self, session, run_id: str, adapter_name: str):
        return (
            session.query(self.adapter_run_model)
            .filter(
                self.adapter_run_model.run_id == run_id,
                self.adapter_run_model.adapter_name == adapter_name,
            )
            .one_or_none()
        ) or self.adapter_run_model(
            run_id=run_id,
            database_name=self.database_name,
            adapter_name=adapter_name,
            status="pending",
        )

    def _upsert_adapter_run_metadata(
        self,
        run_id: str,
        adapter_name: str,
        *,
        status: str | None = None,
        adapter_position: int | None = None,
        adapter_total: int | None = None,
        records_written: int | None = None,
        started_at: datetime | None = None,
        completed_at: datetime | None = None,
        failed_at: datetime | None = None,
        error_message: str | None = None,
    ) -> None:
        if not self._supports_adapter_run_metadata():
            return

        session = self.get_session()
        try:
            row = self._get_or_create_adapter_run_row(session, run_id, adapter_name)
            if row.id is None:
                session.add(row)

            row.database_name = self.database_name
            if status is not None:
                row.status = status
            if adapter_position is not None:
                row.adapter_position = adapter_position
            if adapter_total is not None:
                row.adapter_total = adapter_total
            if records_written is not None:
                row.records_written = records_written
            if started_at is not None:
                row.started_at = started_at
            if completed_at is not None:
                row.completed_at = completed_at
            if failed_at is not None:
                row.failed_at = failed_at
            if error_message is not None:
                row.error_message = error_message

            stats = self._current_adapter_stats if (
                run_id == self._current_run_id and adapter_name == self._current_adapter_name
            ) else None
            if stats:
                row.store_calls = stats["store_calls"]
                row.input_object_count = stats["input_object_count"]
                row.converted_object_count = stats["converted_object_count"]
                row.inserted_row_count = stats["inserted_row_count"]
                row.conversion_seconds = stats["conversion_seconds"]
                row.serialization_seconds = stats["serialization_seconds"]
                row.insert_seconds = stats["insert_seconds"]
                row.total_store_seconds = stats["total_store_seconds"]

            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()

    def store(self, objects, single_source=False, field_conflict_behavior=None) -> bool:
        if not isinstance(objects, list):
            objects = [objects]

        store_start = datetime.now()
        stats = self._current_adapter_stats
        if stats is not None:
            stats["store_calls"] += 1
            stats["input_object_count"] += len(objects)

        object_groups = self.sort_and_convert_objects(objects, keep_nested_objects=True)
        try:
            for obj_list, labels, is_relationship, start_labels, end_labels, obj_cls in object_groups.values():
                converters = self.output_converter.get_object_converters(obj_cls)
                if converters is None:
                    continue

                if not isinstance(converters, list):
                    converters = [converters]

                for converter in converters:
                    start_time = datetime.now()
                    table_class = None
                    stmt = None
                    inserted_count = 0

                    for obj_chunk in self._chunked(obj_list, self.conversion_chunk_size):
                        conversion_start = datetime.now()
                        converted_objects = []
                        for obj in obj_chunk:
                            result = converter(obj)
                            if isinstance(result, list):
                                converted_objects.extend(result)
                            elif result is not None:
                                converted_objects.append(result)
                        if stats is not None:
                            stats["conversion_seconds"] += (datetime.now() - conversion_start).total_seconds()
                            stats["converted_object_count"] += len(converted_objects)

                        if not converted_objects:
                            continue

                        if table_class is None:
                            table_class = converted_objects[0].__class__
                            stmt = mysql_insert(table_class.__table__)
                            if table_class is TinxImportance:
                                stmt = stmt.on_duplicate_key_update(
                                    score=func.greatest(table_class.score, stmt.inserted.score),
                                    doid=case(
                                        (stmt.inserted.score > table_class.score, stmt.inserted.doid),
                                        else_=table_class.doid,
                                    ),
                                )
                            if table_class is WordCount:
                                stmt = stmt.on_duplicate_key_update(count=stmt.inserted.count)
                            print(f"Inserting objects of type {table_class.__name__}")

                        serialization_start = datetime.now()
                        rows = self._serialize_rows(converted_objects)
                        if stats is not None:
                            stats["serialization_seconds"] += (datetime.now() - serialization_start).total_seconds()
                        inserted_count += len(rows)

                        for row_chunk in self._chunked(rows, self.insert_batch_size):
                            self._execute_insert_chunk(stmt, table_class, row_chunk, stats)

                    if table_class is None:
                        continue

                    print(f"Inserted {inserted_count} objects of type {table_class.__name__}")
                    duration = (datetime.now() - start_time).total_seconds()
                    print(f"Processed {inserted_count} objects in {duration:.2f} seconds.")

            return True

        except Exception as e:
            print("Error during insert:", e)
            raise

        finally:
            if stats is not None:
                stats["total_store_seconds"] += (datetime.now() - store_start).total_seconds()

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        effective_truncate = self.truncate_tables if truncate_tables is None else truncate_tables
        self.recreate_mysql_db(self.database_name, effective_truncate)
        return True

    def get_completed_adapter_names(self, run_id: str) -> set[str]:
        if not self._supports_adapter_run_metadata():
            return set()

        session = self.get_session()
        try:
            rows = (
                session.query(self.adapter_run_model.adapter_name)
                .filter(
                    self.adapter_run_model.run_id == run_id,
                    self.adapter_run_model.status == "completed",
                )
                .all()
            )
            return {row[0] for row in rows}
        finally:
            session.close()

    def reset_run_state(self, run_id: str) -> None:
        if not self._supports_adapter_run_metadata():
            return

        session = self.get_session()
        try:
            (
                session.query(self.adapter_run_model)
                .filter(self.adapter_run_model.run_id == run_id)
                .delete(synchronize_session=False)
            )
            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()

        if self._current_run_id == run_id:
            self._current_run_id = None
            self._current_adapter_name = None
            self._current_adapter_stats = None

    def mark_adapter_running(self, run_id: str, adapter_name: str, adapter_position: int | None = None,
                             adapter_total: int | None = None) -> None:
        self._current_run_id = run_id
        self._current_adapter_name = adapter_name
        self._current_adapter_stats = self._new_adapter_stats()
        self._upsert_adapter_run_metadata(
            run_id,
            adapter_name,
            status="running",
            adapter_position=adapter_position,
            adapter_total=adapter_total,
            started_at=datetime.utcnow(),
        )

    def mark_adapter_completed(self, run_id: str, adapter_name: str, records_written: int = 0,
                               adapter_position: int | None = None, adapter_total: int | None = None) -> None:
        self._upsert_adapter_run_metadata(
            run_id,
            adapter_name,
            status="completed",
            adapter_position=adapter_position,
            adapter_total=adapter_total,
            records_written=records_written,
            completed_at=datetime.utcnow(),
        )
        if self._current_run_id == run_id and self._current_adapter_name == adapter_name:
            self._current_run_id = None
            self._current_adapter_name = None
            self._current_adapter_stats = None

    def mark_adapter_failed(self, run_id: str, adapter_name: str, error_message: str | None = None,
                            adapter_position: int | None = None, adapter_total: int | None = None) -> None:
        self._upsert_adapter_run_metadata(
            run_id,
            adapter_name,
            status="failed",
            adapter_position=adapter_position,
            adapter_total=adapter_total,
            failed_at=datetime.utcnow(),
            error_message=error_message,
        )
        if self._current_run_id == run_id and self._current_adapter_name == adapter_name:
            self._current_run_id = None
            self._current_adapter_name = None
            self._current_adapter_stats = None

    def flush_incremental_metadata(self) -> None:
        if self._current_run_id and self._current_adapter_name:
            self._upsert_adapter_run_metadata(
                self._current_run_id,
                self._current_adapter_name,
                status="running",
            )


class TestOutputAdapter(MySQLOutputAdapter):
    output_converter: TestSQLOutputConverter

    def __init__(self, credentials: DBCredentials, database_name: str):
        MySQLOutputAdapter.__init__(self, credentials, database_name)
        self.output_converter = TestSQLOutputConverter(sql_base=TestBase)

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        super().create_or_truncate_datastore(truncate_tables=truncate_tables)
        self.output_converter.sql_base.metadata.create_all(self.get_engine())
        return True


class TCRDOutputAdapter(MySQLOutputAdapter):
    output_converter = TCRDOutputConverter
    adapter_run_model = ETLAdapterRun

    def __init__(
        self,
        credentials: DBCredentials,
        database_name: str,
        truncate_tables: bool,
        source_graph_credentials: DBCredentials | dict | None = None,
        source_graph_database: str | None = None,
    ):
        MySQLOutputAdapter.__init__(
            self,
            credentials,
            database_name,
            truncate_tables=truncate_tables,
        )
        self.output_converter = TCRDOutputConverter()
        self.source_graph_credentials = self._coerce_db_credentials(source_graph_credentials)
        self.source_graph_database = source_graph_database
        self._resolver_fingerprints_by_type = {}
        self._resolver_source_yaml = None

    def set_resolver_metadata(self, resolver_fingerprints_by_type=None, source_yaml=None):
        self._resolver_fingerprints_by_type = resolver_fingerprints_by_type or {}
        self._resolver_source_yaml = source_yaml

    def do_pre_processing(self) -> None:
        self._validate_source_graph_resolver_metadata()
        session = self.get_session()
        try:
            self.output_converter.preload_id_mappings(session)
            self._preload_graph_backed_id_mappings(session)
        finally:
            session.close()

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        self._validate_source_graph_resolver_metadata()
        super().create_or_truncate_datastore(truncate_tables=truncate_tables)
        self.output_converter.sql_base.metadata.create_all(self.get_engine())
        session = self.get_session()
        try:
            self.output_converter.preload_id_mappings(session)
            self._preload_graph_backed_id_mappings(session)
        finally:
            session.close()
        return True

    def _validate_source_graph_resolver_metadata(self) -> None:
        if self.source_graph_credentials is None or not self.source_graph_database:
            return
        if not self._resolver_fingerprints_by_type:
            return

        graph_etl_metadata = self._load_graph_etl_metadata()
        graph_resolver_metadata = graph_etl_metadata.get("resolver_metadata") or {}
        graph_fingerprints = graph_resolver_metadata.get("by_type")
        if not graph_fingerprints:
            raise RuntimeError(
                "Source graph resolver metadata is missing. Rebuild the source graph with resolver "
                "fingerprint metadata before running TCRD MySQL conversion, or disable resume and "
                "use a graph whose resolver provenance is known."
            )

        mismatches = []
        for node_type, expected in self._resolver_fingerprints_by_type.items():
            actual = graph_fingerprints.get(node_type)
            if actual is None:
                mismatches.append({
                    "type": node_type,
                    "expected": expected,
                    "actual": None,
                })
                continue
            if actual.get("fingerprint") != expected.get("fingerprint"):
                mismatches.append({
                    "type": node_type,
                    "expected": expected,
                    "actual": actual,
                })

        if mismatches:
            expected_summary = resolver_fingerprint_summary(self._resolver_fingerprints_by_type)
            actual_summary = resolver_fingerprint_summary(graph_fingerprints)
            raise RuntimeError(
                "Source graph resolver metadata does not match this MySQL conversion config. "
                f"source_graph={self.source_graph_database}; "
                f"source_graph_yaml={graph_resolver_metadata.get('source_yaml') or graph_etl_metadata.get('source_yaml')}; "
                f"source_graph_resolver_yamls={graph_resolver_metadata.get('source_yamls') or []}; "
                f"mysql_yaml={self._resolver_source_yaml}; "
                f"mismatched_types={[row['type'] for row in mismatches]}; "
                f"graph_resolvers={actual_summary}; mysql_resolvers={expected_summary}"
            )

    def _preload_graph_backed_id_mappings(self, session) -> None:
        if not isinstance(self.output_converter, TCRDOutputConverter):
            return

        converter = self.output_converter
        ncats_mapping = converter.id_mapping.setdefault("ncats_disease", {})

        # Direct MONDO-backed mapping is stable and should always win when available.
        for disease_id, mondoid in session.query(NcatsDisease.id, NcatsDisease.mondoid).filter(
            NcatsDisease.mondoid.isnot(None)
        ):
            if mondoid:
                ncats_mapping[mondoid] = disease_id

        # Resume-safe mappings for other generated IDs can come straight from MySQL.
        converter.id_mapping["tissue"] = {
            name: tissue_id
            for tissue_id, name in session.query(Tissue.id, Tissue.name).all()
            if name
        }
        converter.id_mapping["panther_class"] = {
            pcid: class_id
            for class_id, pcid in session.query(PantherClass.id, PantherClass.pcid).all()
            if pcid
        }

        if self.source_graph_credentials is None or not self.source_graph_database:
            return

        existing_rows = session.query(
            NcatsDisease.id,
            NcatsDisease.name,
            NcatsDisease.uniprot_description,
            NcatsDisease.do_description,
            NcatsDisease.mondo_description,
            NcatsDisease.mondoid,
        ).all()

        signature_to_id = {}
        ambiguous_signatures = set()
        for row in existing_rows:
            signature = (
                row.name,
                row.uniprot_description,
                row.do_description,
                row.mondo_description,
            )
            if signature in ambiguous_signatures:
                continue
            if signature in signature_to_id and signature_to_id[signature] != row.id:
                ambiguous_signatures.add(signature)
                del signature_to_id[signature]
                continue
            signature_to_id[signature] = row.id

        adapter = ArangoAdapter(self.source_graph_credentials, self.source_graph_database)
        if adapter.get_db().has_collection("VirusViralProteinEdge"):
            converter.id_mapping["viral_protein_to_virus"] = {
                row["start_id"]: row["end_id"].split(":", 1)[1]
                for row in adapter.runQuery("""
                    FOR rel IN `VirusViralProteinEdge`
                        RETURN KEEP(rel, "start_id", "end_id")
                """)
                if row.get("start_id") and row.get("end_id")
            }

        if not adapter.get_db().has_collection("Disease"):
            return
        disease_rows = adapter.runQuery("""
            FOR d IN `Disease`
                RETURN KEEP(d, "id", "name", "uniprot_description", "do_description", "mondo_description")
        """)

        mapped = 0
        for row in disease_rows:
            disease_id = row.get("id")
            if not disease_id or disease_id in ncats_mapping:
                continue

            signature = (
                row.get("name") or disease_id,
                row.get("uniprot_description"),
                row.get("do_description"),
                row.get("mondo_description"),
            )
            existing_id = signature_to_id.get(signature)
            if existing_id is None:
                continue
            ncats_mapping[disease_id] = existing_id
            mapped += 1

        if mapped:
            print(f"Preloaded {mapped} graph-backed ncats_disease ids for resume safety")

    @staticmethod
    def _build_transitive_closure(node_ids, direct_edges):
        parent_map = {node_id: set() for node_id in node_ids}
        for child_id, parent_id in direct_edges:
            if child_id in parent_map and parent_id in parent_map:
                parent_map[child_id].add(parent_id)

        closure = set()
        for node_id in node_ids:
            closure.add((node_id, node_id))
            stack = list(parent_map.get(node_id, ()))
            visited = set()
            while stack:
                ancestor_id = stack.pop()
                if ancestor_id in visited:
                    continue
                visited.add(ancestor_id)
                closure.add((node_id, ancestor_id))
                stack.extend(parent_map.get(ancestor_id, ()))
        return sorted(closure)

    @staticmethod
    def _coerce_db_credentials(raw_credentials):
        if raw_credentials is None:
            return None
        if isinstance(raw_credentials, DBCredentials):
            return raw_credentials
        if isinstance(raw_credentials, dict):
            return DBCredentials.from_yaml(raw_credentials)
        raise TypeError(f"Unsupported credential type: {type(raw_credentials)}")

    @staticmethod
    def _build_data_source_version_rows(database_metadata: DatabaseMetadata):
        row_map = {}
        for collection in database_metadata.collections:
            for source in collection.sources:
                key = (
                    source.name,
                    source.version,
                    source.version_date,
                    source.download_date,
                )
                if key not in row_map:
                    row_map[key] = {
                        "data_source": source.name,
                        "version": source.version,
                        "version_date": source.version_date,
                        "download_date": source.download_date,
                        "collections": set(),
                    }
                row_map[key]["collections"].add(collection.name)

        rows = []
        for key in sorted(
            row_map,
            key=lambda item: (
                item[0],
                item[1] or "",
                item[2].isoformat() if item[2] else "",
                item[3].isoformat() if item[3] else "",
            ),
        ):
            row = row_map[key]
            rows.append({
                "data_source": row["data_source"],
                "version": row["version"],
                "version_date": row["version_date"],
                "download_date": row["download_date"],
                "collections": ", ".join(sorted(row["collections"])),
            })
        return rows

    def _load_graph_database_metadata(self) -> DatabaseMetadata:
        if self.source_graph_credentials is None or not self.source_graph_database:
            return DatabaseMetadata(collections=[])

        adapter = ArangoAdapter(self.source_graph_credentials, self.source_graph_database)
        store = adapter.get_metadata_store(False)
        metadata_doc = store.get("database_metadata") or {}
        return DatabaseMetadata.from_dict(metadata_doc.get("value", []))

    def _load_graph_etl_metadata(self) -> dict:
        if self.source_graph_credentials is None or not self.source_graph_database:
            return {}

        adapter = ArangoAdapter(self.source_graph_credentials, self.source_graph_database)
        store = adapter.get_metadata_store(False)
        metadata_doc = store.get("etl_metadata") or {}
        return metadata_doc.get("value", {}) or {}

    @staticmethod
    def _build_etl_run_rows(graph_etl_metadata: dict, mysql_database_name: str):
        rows = []

        if graph_etl_metadata and graph_etl_metadata.get("run_date"):
            graph_git = graph_etl_metadata.get("git_info", {}) or {}
            rows.append({
                "stage": "graph_build",
                "database_name": mysql_database_name,
                "run_date": datetime.fromisoformat(graph_etl_metadata["run_date"]),
                "git_commit": graph_git.get("git_commit"),
                "git_branch": graph_git.get("git_branch"),
                "git_tag": graph_git.get("git_tag"),
                "runner": graph_etl_metadata.get("runner"),
                "hostname": graph_etl_metadata.get("hostname"),
                "platform": graph_etl_metadata.get("platform"),
                "platform_version": graph_etl_metadata.get("platform_version"),
                "python_version": graph_etl_metadata.get("python_version"),
            })

        mysql_git = get_git_metadata()
        rows.append({
            "stage": "graph_to_mysql",
            "database_name": mysql_database_name,
            "run_date": datetime.now(),
            "git_commit": mysql_git.get("git_commit"),
            "git_branch": mysql_git.get("git_branch"),
            "git_tag": mysql_git.get("git_tag"),
            "runner": os.getenv("USER", "unknown"),
            "hostname": socket.gethostname(),
            "platform": platform.system(),
            "platform_version": platform.version(),
            "python_version": platform.python_version(),
        })

        return rows

    def _populate_data_source_version_table(self, session):
        database_metadata = self._load_graph_database_metadata()
        rows = self._build_data_source_version_rows(database_metadata)
        session.query(DataSourceVersion).delete()
        if rows:
            session.bulk_insert_mappings(DataSourceVersion, rows)

    def _populate_etl_run_table(self, session):
        graph_etl_metadata = self._load_graph_etl_metadata()
        rows = self._build_etl_run_rows(graph_etl_metadata, self.database_name)
        session.query(ETLRun).delete()
        if rows:
            session.bulk_insert_mappings(ETLRun, rows)

    def _populate_ancestry_table(self, session, node_cls, parent_cls, ancestry_cls, node_key, parent_key):
        node_ids = [row[0] for row in session.query(node_cls).with_entities(getattr(node_cls, node_key)).all() if row[0]]
        direct_edges = session.query(parent_cls).with_entities(
            getattr(parent_cls, node_key),
            getattr(parent_cls, parent_key),
        ).all()
        closure = self._build_transitive_closure(node_ids, direct_edges)
        session.query(ancestry_cls).delete()
        if closure:
            session.bulk_insert_mappings(
                ancestry_cls,
                [{"oid": oid, "ancestor_id": ancestor_id} for oid, ancestor_id in closure],
            )

    @staticmethod
    def _populate_dto_parent_column(session):
        session.query(DTO).update({DTO.parent_id: None})
        direct_edges = session.query(DTOParent.dtoid, DTOParent.parent_id).all()
        if direct_edges:
            session.bulk_update_mappings(
                DTO,
                [{"dtoid": dtoid, "parent_id": parent_id} for dtoid, parent_id in direct_edges],
            )

    @staticmethod
    def _populate_ncats_disease_summary_fields(session):
        session.execute(text("""
            UPDATE `ncats_disease`
            SET
              `target_count` = 0,
              `direct_target_count` = 0,
              `maxTDL` = NULL
        """))

        session.execute(text("""
            UPDATE `ncats_disease` n
            JOIN (
              SELECT
                x.`ncats_disease_id`,
                COUNT(DISTINCT d.`protein_id`) AS `direct_target_count`
              FROM `ncats_d2da` x
              JOIN `disease` d
                ON d.`id` = x.`disease_assoc_id`
              WHERE x.`direct` = 1
                AND d.`protein_id` IS NOT NULL
              GROUP BY x.`ncats_disease_id`
            ) q
              ON q.`ncats_disease_id` = n.`id`
            SET n.`direct_target_count` = q.`direct_target_count`
        """))

        session.execute(text("DROP TEMPORARY TABLE IF EXISTS `ncats_disease_descendant_stage`"))
        session.execute(text("""
            CREATE TEMPORARY TABLE `ncats_disease_descendant_stage` (
              `ancestor_id` INT NOT NULL,
              `descendant_id` INT NOT NULL,
              PRIMARY KEY (`ancestor_id`, `descendant_id`)
            )
        """))
        session.execute(text("""
            INSERT IGNORE INTO `ncats_disease_descendant_stage`
              (`ancestor_id`, `descendant_id`)
            SELECT `id`, `id`
            FROM `ncats_disease`
        """))
        session.execute(text("""
            INSERT IGNORE INTO `ncats_disease_descendant_stage`
              (`ancestor_id`, `descendant_id`)
            SELECT
              parent.`id` AS `ancestor_id`,
              child.`id` AS `descendant_id`
            FROM `ncats_disease` parent
            JOIN `ancestry_mondo` ancestry
              ON ancestry.`ancestor_id` = parent.`mondoid`
            JOIN `ncats_disease` child
              ON child.`mondoid` = ancestry.`oid`
            WHERE parent.`mondoid` IS NOT NULL
        """))
        session.execute(text("""
            UPDATE `ncats_disease` n
            JOIN (
              SELECT
                stage.`ancestor_id` AS `ncats_disease_id`,
                COUNT(DISTINCT d.`protein_id`) AS `target_count`
              FROM `ncats_disease_descendant_stage` stage
              JOIN `ncats_d2da` x
                ON x.`ncats_disease_id` = stage.`descendant_id`
              JOIN `disease` d
                ON d.`id` = x.`disease_assoc_id`
              WHERE d.`protein_id` IS NOT NULL
              GROUP BY stage.`ancestor_id`
            ) q
              ON q.`ncats_disease_id` = n.`id`
            SET n.`target_count` = q.`target_count`
        """))
        session.execute(text("DROP TEMPORARY TABLE IF EXISTS `ncats_disease_descendant_stage`"))

        session.execute(text("""
            UPDATE `ncats_disease` n
            JOIN (
              SELECT
                n2.`id`,
                MAX(CASE
                  WHEN tgt.`tdl` = 'Tclin' THEN 4
                  WHEN tgt.`tdl` = 'Tchem' THEN 3
                  WHEN tgt.`tdl` = 'Tbio' THEN 2
                  ELSE 1
                END) AS `tempTDL`
              FROM `ncats_disease` n2
              JOIN `ncats_d2da` x
                ON n2.`id` = x.`ncats_disease_id`
              JOIN `disease` d
                ON x.`disease_assoc_id` = d.`id`
              JOIN `t2tc`
                ON d.`protein_id` = `t2tc`.`protein_id`
              JOIN `target` tgt
                ON `t2tc`.`target_id` = tgt.`id`
              WHERE x.`direct` = 1
              GROUP BY n2.`id`
            ) q
              ON q.`id` = n.`id`
            SET n.`maxTDL` = CASE
              WHEN q.`tempTDL` = 4 THEN 'Tclin'
              WHEN q.`tempTDL` = 3 THEN 'Tchem'
              WHEN q.`tempTDL` = 2 THEN 'Tbio'
              ELSE 'Tdark'
            END
        """))

    @staticmethod
    def _populate_ncats_ligand_summary_fields(session):
        session.execute(text("""
            UPDATE `ncats_ligands`
            SET
              `actCnt` = 0,
              `targetCount` = 0
        """))
        session.execute(text("""
            UPDATE `ncats_ligands` l
            JOIN (
              SELECT
                a.`ncats_ligand_id`,
                COUNT(a.`id`) AS `actCnt`,
                COUNT(DISTINCT a.`target_id`) AS `targetCount`
              FROM `ncats_ligand_activity` a
              GROUP BY a.`ncats_ligand_id`
            ) q
              ON q.`ncats_ligand_id` = l.`id`
            SET
              l.`actCnt` = q.`actCnt`,
              l.`targetCount` = q.`targetCount`
        """))

    @staticmethod
    def _populate_typeahead_index(session):
        session.execute(text("DELETE FROM `ncats_typeahead_index`"))
        for insert_sql in TYPEAHEAD_INSERTS:
            session.execute(text(insert_sql))

    def do_post_processing(self, clean_edges: bool = True) -> None:
        NcatsTypeaheadIndex.__table__.create(self.get_engine(), checkfirst=True)
        session = self.get_session()
        try:
            self._populate_data_source_version_table(session)
            self._populate_etl_run_table(session)
            self._populate_dto_parent_column(session)
            self._populate_ancestry_table(
                session=session,
                node_cls=DO,
                parent_cls=DOParent,
                ancestry_cls=AncestryDO,
                node_key="doid",
                parent_key="parent_id",
            )
            self._populate_ancestry_table(
                session=session,
                node_cls=DTO,
                parent_cls=DTOParent,
                ancestry_cls=AncestryDTO,
                node_key="dtoid",
                parent_key="parent_id",
            )
            self._populate_ancestry_table(
                session=session,
                node_cls=Mondo,
                parent_cls=MondoParent,
                ancestry_cls=AncestryMONDO,
                node_key="mondoid",
                parent_key="parentid",
            )
            self._populate_ancestry_table(
                session=session,
                node_cls=Uberon,
                parent_cls=UberonParent,
                ancestry_cls=AncestryUBERON,
                node_key="uid",
                parent_key="parent_id",
            )
            self._populate_ncats_disease_summary_fields(session)
            self._populate_ncats_ligand_summary_fields(session)
            self._populate_typeahead_index(session)
            session.commit()
        except Exception:
            session.rollback()
            raise
        finally:
            session.close()
