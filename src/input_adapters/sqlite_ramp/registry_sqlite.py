"""Materialize a pinned RaMP Registry gzip as a reusable SQLite file."""

import fcntl
import gzip
import hashlib
import os
from pathlib import Path
import shutil
import sqlite3
import tempfile

from src.models.registry_dataset import RegistryDataset, RegistryDatasetKind


def ramp_sqlite_file(data_source: RegistryDataset) -> Path:
    if data_source.kind is not RegistryDatasetKind.SOURCE or (
        data_source.source, data_source.dataset
    ) != ("ramp", "sqlite_database"):
        raise ValueError("RaMP SQLite requires a ramp:sqlite_database source snapshot")

    artifact = data_source.file()
    files = data_source.manifest.get("files") or []
    if not artifact.name.endswith(".sqlite.gz"):
        raise ValueError(f"RaMP Registry artifact is not a gzipped SQLite file: {artifact}")
    expected_sha = files[0].get("sha256") if len(files) == 1 and isinstance(files[0], dict) else None
    if not isinstance(expected_sha, str) or len(expected_sha) != 64 or any(c not in "0123456789abcdef" for c in expected_sha.lower()):
        raise ValueError(f"RaMP Registry snapshot {data_source.snapshot_id} lacks a valid file SHA-256")
    expected_sha = expected_sha.lower()

    cache_dir = Path(tempfile.gettempdir()) / "ifx-odin-ramp-sqlite"
    cache_dir.mkdir(parents=True, exist_ok=True)
    target = cache_dir / f"{expected_sha}.sqlite"
    with (cache_dir / f"{expected_sha}.lock").open("a+b") as lock:
        fcntl.flock(lock, fcntl.LOCK_EX)
        if target.is_file():
            _validate_sqlite(target, data_source)
            return target

        digest = hashlib.sha256()
        temporary = None
        try:
            with tempfile.NamedTemporaryFile(dir=cache_dir, suffix=".sqlite", delete=False) as output:
                temporary = Path(output.name)
                try:
                    with artifact.open("rb") as source:
                        reader = _HashingReader(source, digest)
                        with gzip.GzipFile(fileobj=reader) as decoded:
                            shutil.copyfileobj(decoded, output, length=4 * 1024 * 1024)
                except (OSError, EOFError) as error:
                    raise ValueError(f"Cannot materialize RaMP Registry SQLite from {artifact}") from error
                output.flush()
                os.fsync(output.fileno())
            if digest.hexdigest() != expected_sha:
                raise ValueError(f"RaMP Registry artifact SHA-256 mismatch: {artifact}")
            _validate_sqlite(temporary, data_source)
            os.replace(temporary, target)
            return target
        finally:
            if temporary is not None:
                temporary.unlink(missing_ok=True)


class _HashingReader:
    def __init__(self, source, digest):
        self.source = source
        self.digest = digest

    def read(self, size=-1):
        chunk = self.source.read(size)
        self.digest.update(chunk)
        return chunk


def _validate_sqlite(path: Path, data_source: RegistryDataset) -> None:
    with path.open("rb") as handle:
        header = handle.read(16)
    if path.stat().st_size < 100 or header != b"SQLite format 3\x00":
        raise ValueError(f"RaMP Registry artifact is not a SQLite database: {path}")
    try:
        with sqlite3.connect(f"{path.resolve().as_uri()}?mode=ro", uri=True) as db:
            versions = db.execute("SELECT ramp_version FROM db_version").fetchall()
            if len(versions) != 1 or versions[0][0] != data_source.version.split("-", 1)[0]:
                raise ValueError(
                    f"RaMP SQLite db_version {versions!r} disagrees with Registry snapshot "
                    f"{data_source.snapshot_id}"
                )
            for table in ("source", "analyte", "pathway", "analytehaspathway"):
                if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone() is None:
                    raise ValueError(f"RaMP SQLite is missing required table {table}: {path}")
    except sqlite3.DatabaseError as error:
        raise ValueError(f"RaMP Registry SQLite is invalid: {path}") from error
