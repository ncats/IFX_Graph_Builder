"""A stage-based SQLite build must finish both post-processing passes before publication."""

import json
import sqlite3
from types import SimpleNamespace

import pytest

from src.use_cases.ramp import build_sqlite


def test_stage_command_defaults_to_complete_sqlite():
    args = build_sqlite.parser().parse_args([
        "--stage-id", "stage-07-test", "--output", "ramp.sqlite", "--overwrite"])
    assert not args.diagnostic
    assert not args.base_only
    assert args.overwrite


def test_complete_postprocessing_runs_human_and_pathway_passes(monkeypatch, tmp_path):
    from src.use_cases.ramp import postprocess_sqlite

    calls = []
    monkeypatch.setattr(postprocess_sqlite, "main", lambda argv: calls.append(argv) or 0)
    args = build_sqlite.parser().parse_args([
        "--stage-id", "stage-test", "--output", str(tmp_path / "ramp.sqlite")])
    build_sqlite._postprocess_complete(tmp_path / "building.sqlite", args, calls.append)
    assert calls[0].startswith("Calculating entity status")
    assert "--with-human-flags" in calls[1]
    assert calls[2].startswith("Calculating pathway similarity")
    assert "--pathway-only" in calls[3]


def test_failed_postprocessing_preserves_existing_output(monkeypatch, tmp_path):
    output = tmp_path / "ramp.sqlite"
    output.write_bytes(b"previous database")
    temporary = tmp_path / "building.sqlite"
    temporary.write_bytes(b"new base tables")
    args = SimpleNamespace(stage_id="stage-test", output=output, overwrite=True)
    monkeypatch.setattr(build_sqlite, "_postprocess_complete",
                        lambda *_: (_ for _ in ()).throw(RuntimeError("pathway pass failed")))
    timer = SimpleNamespace(mark=lambda _: None)
    with pytest.raises(RuntimeError, match="pathway pass failed"):
        build_sqlite._finish_complete(temporary, args, timer)
    assert output.read_bytes() == b"previous database"


def test_complete_manifest_and_atomic_publish(tmp_path):
    output = tmp_path / "ramp.sqlite"
    output.write_bytes(b"previous database")
    temporary = tmp_path / "building.sqlite"
    with sqlite3.connect(temporary) as db:
        db.execute("CREATE TABLE ramp_export_metadata (key TEXT, value TEXT)")
        db.execute("INSERT INTO ramp_export_metadata VALUES (?,?)", (
            "manifest", json.dumps({
                "stage_id": "HarmonizationStage:stage-test", "scope": "base_tables",
                "status": "complete",
                "pending_tables": [], "pending_fields": [],
                "post_processing": {
                    "pathway_tables": {"status": "computed"},
                    "human_reaction_flags": {"policy": "computed"},
                },
            })))
    assert build_sqlite._validate_complete(temporary, "stage-test")["scope"] == "base_tables"
    build_sqlite._publish_complete(temporary, output, overwrite=True)
    assert not temporary.exists()
    assert output.read_bytes() != b"previous database"


def test_incomplete_manifest_is_not_published(tmp_path):
    temporary = tmp_path / "building.sqlite"
    with sqlite3.connect(temporary) as db:
        db.execute("CREATE TABLE ramp_export_metadata (key TEXT, value TEXT)")
        db.execute("INSERT INTO ramp_export_metadata VALUES (?,?)", (
            "manifest", json.dumps({
                "stage_id": "HarmonizationStage:stage-test", "scope": "base_tables",
                "pending_tables": ["pathway_similarity"], "pending_fields": [],
            })))
    with pytest.raises(ValueError, match="pending"):
        build_sqlite._validate_complete(temporary, "stage-test")


def test_main_runs_postprocessing_and_report_after_replacing_output(tmp_path, monkeypatch):
    from src.core.registry_integration import RegistryIntegration
    from src.infrastructure import object_storage
    from src.shared import arango_adapter
    from src.use_cases.ramp import sqlite_gene_identity, sqlite_protein_annotations, sqlite_stage

    credentials = tmp_path / "graph.yaml"
    credentials.write_text("url: http://unused\nuser: unused\npassword: unused\n")
    output = tmp_path / "ramp.sqlite"
    output.write_bytes(b"previous database")
    previous_report = tmp_path / "comparison.json"
    previous_report.write_text("{}")
    events = []
    monkeypatch.setattr(arango_adapter, "ArangoAdapter",
                        lambda *_args: SimpleNamespace(get_db=lambda: object()))
    monkeypatch.setattr(object_storage, "load_object_storage_credentials", lambda _: object())
    monkeypatch.setattr(object_storage, "object_storage_from_credentials", lambda _: object())
    monkeypatch.setattr(sqlite_stage, "StageReader", lambda *_args: SimpleNamespace(metadata={}))
    monkeypatch.setattr(RegistryIntegration, "connect", lambda *_args: object())
    identity = SimpleNamespace(input_file=tmp_path / "unused", provenance={})
    monkeypatch.setattr(sqlite_gene_identity.GeneIdentity, "from_stage",
                        lambda *_args, **_kw: identity)
    monkeypatch.setattr(sqlite_protein_annotations.ProteinAnnotations, "from_file",
                        lambda *_args: object())

    def export(_reader, path, **kwargs):
        assert kwargs["source_only"] is False
        events.append("export")
        path.write_bytes(b"base")
        return {"row_counts": {}, "pathway_source_ids_at_cutoff": 0,
                "pathway_assertions_excluded": 0,
                "lipidmaps_class_level4_assertions_excluded": 0}

    def postprocess(path, _args, _progress):
        assert output.read_bytes() == b"previous database"
        assert path.read_bytes() == b"base"
        events.append("postprocess")
        path.write_bytes(b"complete")

    monkeypatch.setattr(build_sqlite, "export_sqlite", export)
    monkeypatch.setattr(build_sqlite, "_postprocess_complete", postprocess)
    monkeypatch.setattr(build_sqlite, "_validate_complete",
                        lambda _path, _stage: {"row_counts": {"entity_status_info": 3}})
    def refresh(_previous, sqlite_path, _html):
        assert _previous == previous_report
        assert sqlite_path == output
        assert output.read_bytes() == b"complete"
        events.append("report")
        return [tmp_path / "comparison.html", previous_report]
    monkeypatch.setattr(build_sqlite, "_refresh_report_from", refresh)
    assert build_sqlite.main([
        "--stage-id", "stage-test", "--output", str(output),
        "--overwrite", "--graph-credentials", str(credentials),
        "--refresh-report-from", str(previous_report)]) == 0
    assert events == ["export", "postprocess", "report"]
    assert output.read_bytes() == b"complete"
    output.write_bytes(b"previous database")
    monkeypatch.setattr(build_sqlite, "_refresh_report_from",
                        lambda *_args: (_ for _ in ()).throw(FileNotFoundError("missing baseline")))
    assert build_sqlite.main([
        "--stage-id", "stage-test", "--output", str(output),
        "--overwrite", "--graph-credentials", str(credentials),
        "--refresh-report-from", str(previous_report)]) == 2
    assert output.read_bytes() == b"complete"


def test_report_refresh_reuses_ordered_historical_inputs(tmp_path, monkeypatch):
    from src.use_cases.ramp import compare_source

    old_a = tmp_path / "a.sqlite"
    old_b = tmp_path / "b.sqlite"
    latest = tmp_path / "new.sqlite"
    for path in (old_a, old_b, latest):
        path.write_bytes(b"fixture")
    previous = tmp_path / "comparison.json"
    previous.write_text(json.dumps({"databases": [
        {"label": "Released A", "path": str(old_a)},
        {"label": "Released B", "path": str(old_b)},
        {"label": "New harmonized build", "path": "/old/location.sqlite"},
    ]}))
    seen = []
    monkeypatch.setattr(compare_source, "compare",
                        lambda databases: seen.extend(databases) or {"databases": []})
    monkeypatch.setattr(compare_source, "write",
                        lambda _report, output: (output, output.with_suffix(".json")))
    html = tmp_path / "comparison.html"
    assert build_sqlite._refresh_report_from(previous, latest, html) == (html, previous)
    assert [(label, path) for label, path in seen] == [
        ("Released A", old_a), ("Released B", old_b), ("New harmonized build", latest)]


def test_report_refresh_fails_loudly_for_missing_baseline(tmp_path):
    previous = tmp_path / "comparison.json"
    previous.write_text(json.dumps({"databases": [
        {"label": "Released A", "path": str(tmp_path / "missing.sqlite")},
        {"label": "New harmonized build", "path": "/old/location.sqlite"},
    ]}))
    with pytest.raises(FileNotFoundError, match="Historical comparison SQLite"):
        build_sqlite._refresh_report_from(
            previous, tmp_path / "new.sqlite", tmp_path / "comparison.html")


def test_last_pass_marks_manifest_and_db_version_complete(tmp_path):
    from tests.test_ramp_pathway_similarity import make_database
    from src.use_cases.ramp.postprocess_sqlite import postprocess_pathway_similarity

    path = tmp_path / "ramp.sqlite"
    make_database(path)
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE db_version (version_notes TEXT)")
        db.execute("INSERT INTO db_version VALUES ('Base tables only; post-processing pending')")
        manifest = json.loads(db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
        manifest.update(scope="base_tables", status="post_processing_pending",
                        pending_fields=[], post_processing={
                            "human_reaction_flags": {"policy": "computed"}})
        db.execute("UPDATE ramp_export_metadata SET value=? WHERE key='manifest'",
                   (json.dumps(manifest),))
    postprocess_pathway_similarity(path)
    with sqlite3.connect(path) as db:
        manifest = json.loads(db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
        assert manifest["status"] == "complete"
        assert db.execute("SELECT version_notes FROM db_version").fetchone() == (
            "Post-processing complete. See ramp_export_metadata.",)
