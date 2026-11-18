import gzip
import hashlib
import sqlite3
from dataclasses import replace
from pathlib import Path

import pytest

from src.input_adapters.sqlite_ramp import registry_sqlite
from src.input_adapters.sqlite_ramp.metabolite_adapter import MetaboliteAdapter
from src.models.registry_dataset import RegistryDataset, RegistryDatasetKind


def _dataset(tmp_path: Path, *, version: str = "3.0.12-ramp4.0", invalid: bool = False):
    tmp_path.mkdir(parents=True, exist_ok=True)
    sqlite_path = tmp_path / "source.sqlite"
    with sqlite3.connect(sqlite_path) as db:
        db.executescript("""
            CREATE TABLE db_version (ramp_version TEXT, load_timestamp TEXT, version_notes TEXT);
            INSERT INTO db_version VALUES ('3.0.12', '2026-01-01', 'fixture');
            CREATE TABLE version_info (data_source_id TEXT, data_source_name TEXT,
                data_source_url TEXT, data_source_version TEXT, status TEXT);
            CREATE TABLE source (sourceId TEXT, rampId TEXT, geneOrCompound TEXT);
            CREATE TABLE analyte (rampId TEXT, common_name TEXT, type TEXT);
            CREATE TABLE pathway (pathwayRampId TEXT);
            CREATE TABLE analytehaspathway (rampId TEXT, pathwayRampId TEXT);
            INSERT INTO analyte VALUES ('RAMP_C1', 'test metabolite', 'compound');
        """)
    artifact = tmp_path / "RaMP_SQLite_v3.0.12.sqlite.gz"
    with gzip.open(artifact, "wb") as handle:
        handle.write(b"not SQLite" if invalid else sqlite_path.read_bytes())
    digest = hashlib.sha256(artifact.read_bytes()).hexdigest()
    return RegistryDataset(
        kind=RegistryDatasetKind.SOURCE,
        source="ramp",
        dataset="sqlite_database",
        version=version,
        version_date=None,
        download_date=None,
        snapshot_id=f"ramp:sqlite_database:{version}",
        manifest_uri="fixture",
        manifest={"files": [{"path": artifact.name, "sha256": digest}]},
        local_dir=tmp_path,
    )


def test_ramp_registry_adapter_uses_snapshot_version_and_shared_sqlite(tmp_path, monkeypatch):
    monkeypatch.setattr(registry_sqlite.tempfile, "gettempdir", lambda: str(tmp_path / "cache"))
    dataset = _dataset(tmp_path)
    first = MetaboliteAdapter(data_source=dataset)
    second = MetaboliteAdapter(data_source=dataset)
    assert first.connection_string == second.connection_string
    assert first.get_version().version == "3.0.12-ramp4.0"
    assert first.get_database_version().id == "3.0.12"
    assert [m.id for batch in first.get_all() for m in batch] == ["RAMP_C1"]


def test_ramp_registry_rejects_invalid_sqlite_and_version(tmp_path, monkeypatch):
    monkeypatch.setattr(registry_sqlite.tempfile, "gettempdir", lambda: str(tmp_path / "cache"))
    with pytest.raises(ValueError, match="not a SQLite database"):
        registry_sqlite.ramp_sqlite_file(_dataset(tmp_path / "invalid", invalid=True))
    with pytest.raises(ValueError, match="disagrees with Registry snapshot"):
        registry_sqlite.ramp_sqlite_file(_dataset(tmp_path / "newer", version="3.0.13"))


def test_ramp_registry_rejects_checksum_mismatch(tmp_path, monkeypatch):
    monkeypatch.setattr(registry_sqlite.tempfile, "gettempdir", lambda: str(tmp_path / "cache"))
    dataset = _dataset(tmp_path)
    dataset = replace(dataset, manifest={"files": [{"path": dataset.file().name, "sha256": "0" * 64}]})
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        registry_sqlite.ramp_sqlite_file(dataset)


def test_ramp_registry_reports_corrupt_gzip(tmp_path, monkeypatch):
    monkeypatch.setattr(registry_sqlite.tempfile, "gettempdir", lambda: str(tmp_path / "cache"))
    dataset = _dataset(tmp_path)
    dataset.file().write_bytes(b"not gzip")
    with pytest.raises(ValueError, match="Cannot materialize RaMP Registry SQLite"):
        registry_sqlite.ramp_sqlite_file(dataset)
