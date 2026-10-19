"""Label and package a completed RaMP SQLite as a local release candidate."""

import argparse
from datetime import datetime, timezone
import gzip
import hashlib
import json
import os
from pathlib import Path
import re
import shutil
import sqlite3
import tempfile

from src.use_cases.ramp.build_sqlite import _validate_complete


VERSION_PATTERN = re.compile(r"\d+\.\d+\.\d+\Z")
CHUNK_SIZE = 4 * 1024 * 1024


def _sha256(path):
    digest = hashlib.sha256()
    with Path(path).open("rb") as handle:
        for chunk in iter(lambda: handle.read(CHUNK_SIZE), b""):
            digest.update(chunk)
    return digest.hexdigest()


def _read_version_state(path):
    with sqlite3.connect(f"{Path(path).resolve().as_uri()}?mode=ro", uri=True) as db:
        row = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
        if row is None:
            raise ValueError("SQLite export manifest is missing")
        manifest = json.loads(row[0])
        versions = db.execute("SELECT ramp_version FROM db_version").fetchall()
        if len(versions) != 1:
            raise ValueError(f"Expected one db_version row, found {len(versions)}")
        source_versions = db.execute(
            "SELECT DISTINCT ramp_db_version FROM version_info").fetchall()
        if len(source_versions) != 1 or source_versions[0][0] != versions[0][0]:
            raise ValueError("version_info and db_version labels disagree")
        if manifest.get("release_version") != versions[0][0]:
            raise ValueError("Export manifest and db_version labels disagree")
        return manifest


def _validate(path, expected_version=None):
    manifest = _read_version_state(path)
    _validate_complete(path, manifest["stage_id"])
    if expected_version and manifest["release_version"] != expected_version:
        raise ValueError("Candidate version does not match the requested version")
    with sqlite3.connect(f"{Path(path).resolve().as_uri()}?mode=ro", uri=True) as db:
        if db.execute("PRAGMA foreign_key_check").fetchone():
            raise ValueError("SQLite foreign-key check failed")
    return manifest


def package_candidate(source, output, version):
    source, output = Path(source).resolve(), Path(output)
    if not VERSION_PATTERN.fullmatch(version):
        raise ValueError("Version must be a final X.Y.Z label")
    expected_name = f"RaMP_SQLite_v{version}.sqlite.gz"
    if output.name != expected_name:
        raise ValueError(f"Output filename must be {expected_name}")
    if not source.is_file() or source.is_symlink():
        raise ValueError(f"Source must be an ordinary SQLite file: {source}")
    if output.exists() or output.is_symlink():
        raise FileExistsError(f"Candidate already exists: {output}")
    for suffix in ("-wal", "-shm", "-journal"):
        if Path(str(source) + suffix).exists():
            raise ValueError(f"Source has an active SQLite sidecar: {source}{suffix}")
    _validate(source)
    source_sha256 = _sha256(source)
    output.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix="ramp-candidate-", dir=output.parent) as folder:
        candidate = Path(folder) / expected_name.removesuffix(".gz")
        compressed = Path(folder) / expected_name
        shutil.copyfile(source, candidate)
        if _sha256(candidate) != source_sha256:
            raise ValueError("Source SQLite changed while the candidate copy was made")
        original = _validate(candidate)
        with sqlite3.connect(candidate) as db:
            db.execute("BEGIN IMMEDIATE")
            db.execute("UPDATE db_version SET ramp_version=?", (version,))
            db.execute("UPDATE version_info SET ramp_db_version=?", (version,))
            updated = dict(original)
            updated["release_version"] = version
            updated["release_ready"] = True
            updated.pop("candidate_packaging", None)
            db.execute("UPDATE ramp_export_metadata SET value=? WHERE key='manifest'",
                       (json.dumps(updated, sort_keys=True),))
            db.commit()
        checked = _validate(candidate, version)
        candidate_sha256 = _sha256(candidate)
        with candidate.open("rb") as raw, compressed.open("xb") as target:
            with gzip.GzipFile(fileobj=target, mode="wb", filename="", mtime=0) as gz:
                shutil.copyfileobj(raw, gz, CHUNK_SIZE)
            target.flush()
            os.fsync(target.fileno())
        decoded_sha256 = hashlib.sha256()
        with gzip.open(compressed, "rb") as decoded:
            for chunk in iter(lambda: decoded.read(CHUNK_SIZE), b""):
                decoded_sha256.update(chunk)
        if decoded_sha256.hexdigest() != candidate_sha256:
            raise ValueError("Compressed candidate does not reproduce the validated SQLite")
        compressed_sha256 = _sha256(compressed)
        os.link(compressed, output)  # Atomic no-clobber publication.
    return {
        "version": version,
        "packaged_at": datetime.now(timezone.utc).isoformat(),
        "source": str(source),
        "source_sha256": source_sha256,
        "sqlite_sha256": candidate_sha256,
        "gzip_sha256": compressed_sha256,
        "gzip_bytes": output.stat().st_size,
        "output": str(output.resolve()),
        "stage_id": checked["stage_id"],
        "created_at": checked["created_at"],
        "status": checked["status"],
        "release_ready": checked["release_ready"],
    }


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sqlite", required=True, type=Path,
                        help="Completed build_sqlite output; never modified")
    parser.add_argument("--version", required=True, help="Candidate version, e.g. 3.0.13")
    parser.add_argument("--output", required=True, type=Path,
                        help="Versioned .sqlite.gz output; existing files are refused")
    args = parser.parse_args(argv)
    print(json.dumps(package_candidate(args.sqlite, args.output, args.version), indent=2))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
