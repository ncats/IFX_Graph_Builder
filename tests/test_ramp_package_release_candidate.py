"""A candidate must preserve its build while labeling all consumer-facing metadata."""

import gzip
import json
import sqlite3

import pytest

from src.use_cases.ramp.package_release_candidate import package_candidate


def _database(path, *, complete=True):
    manifest = {
        "stage_id": "HarmonizationStage:stage-test", "scope": "base_tables",
        "status": "complete" if complete else "post_processing_pending",
        "pending_tables": [] if complete else ["pathway_similarity"],
        "pending_fields": [], "release_version": "unreleased",
        "release_ready": False, "created_at": "2026-10-13T12:00:00+00:00",
        "post_processing": {
            "pathway_tables": {"status": "computed"},
            "human_reaction_flags": {"policy": "computed"},
        },
    }
    with sqlite3.connect(path) as db:
        db.execute("CREATE TABLE ramp_export_metadata (key TEXT PRIMARY KEY, value TEXT)")
        db.execute("INSERT INTO ramp_export_metadata VALUES ('manifest', ?)",
                   (json.dumps(manifest),))
        db.execute("CREATE TABLE db_version (ramp_version TEXT)")
        db.execute("INSERT INTO db_version VALUES ('unreleased')")
        db.execute("CREATE TABLE version_info (ramp_db_version TEXT)")
        db.executemany("INSERT INTO version_info VALUES (?)",
                       [("unreleased",), ("unreleased",)])


def test_package_labels_copy_without_changing_build(tmp_path):
    source = tmp_path / "base.sqlite"
    output = tmp_path / "RaMP_SQLite_v3.0.13.sqlite.gz"
    _database(source)
    original = source.read_bytes()

    receipt = package_candidate(source, output, "3.0.13")
    candidate = tmp_path / "candidate.sqlite"
    candidate.write_bytes(gzip.decompress(output.read_bytes()))
    with sqlite3.connect(candidate) as db:
        manifest = json.loads(db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
        assert db.execute("SELECT ramp_version FROM db_version").fetchone()[0] == "3.0.13"
        assert db.execute("SELECT DISTINCT ramp_db_version FROM version_info").fetchall() == [("3.0.13",)]
        assert db.execute("PRAGMA integrity_check").fetchone() == ("ok",)
    assert manifest["release_version"] == receipt["version"] == "3.0.13"
    assert manifest["created_at"] == "2026-10-13T12:00:00+00:00"
    assert manifest["release_ready"] is True
    assert "candidate_packaging" not in manifest
    assert receipt["packaged_at"]
    assert source.read_bytes() == original
    with pytest.raises(FileExistsError):
        package_candidate(source, output, "3.0.13")


def test_incomplete_database_cannot_be_packaged(tmp_path):
    source = tmp_path / "base.sqlite"
    output = tmp_path / "RaMP_SQLite_v3.0.13.sqlite.gz"
    _database(source, complete=False)
    with pytest.raises(ValueError, match="pending"):
        package_candidate(source, output, "3.0.13")
    assert not output.exists()


def test_github_artifact_requires_final_version_label(tmp_path):
    source = tmp_path / "base.sqlite"
    _database(source)
    with pytest.raises(ValueError, match="final X.Y.Z"):
        package_candidate(source, tmp_path / "RaMP_SQLite_v3.0.13-rc1.sqlite.gz",
                          "3.0.13-rc1")
