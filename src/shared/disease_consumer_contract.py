"""Validated access to the published DISEASES consumer node file."""

from __future__ import annotations

import csv
import hashlib
from pathlib import Path


REQUIRED_COLUMNS = {
    "ncats_disease_id",
    "primary_xref",
    "harmonized_name",
    "is_rare",
    "rare_basis",
    "MONDO_xref",
    "GARD_xref",
    "non_exact_xrefs",
}


def disease_ids_file(data_source) -> tuple[Path, str]:
    """Return the materialized file only when its pinned checksum is intact."""
    path, expected_sha256 = validated_disease_snapshot_file(data_source, "disease_ids.tsv")
    with path.open(newline="", encoding="utf-8") as handle:
        columns = next(csv.reader(handle, delimiter="\t"), [])
    missing = REQUIRED_COLUMNS - set(columns)
    if missing:
        raise ValueError(f"{data_source.snapshot_id} disease_ids.tsv lacks columns {sorted(missing)}")
    if len(columns) != len(set(columns)):
        raise ValueError(f"{data_source.snapshot_id} disease_ids.tsv has duplicate column names")
    return path, expected_sha256


def validated_disease_snapshot_file(data_source, member: str) -> tuple[Path, str]:
    """Check a DISEASES member against its Registry manifest checksum."""
    snapshot_id = data_source.snapshot_id
    if not snapshot_id.startswith("ifx_harmonizers:diseases:"):
        raise ValueError(f"Expected a DISEASES derived snapshot, got {snapshot_id!r}")
    entries = data_source.manifest.get("files") or []
    entry = next(
        (item for item in entries if isinstance(item, dict) and item.get("path") == member), None
    )
    if not entry or not entry.get("sha256"):
        raise ValueError(f"{snapshot_id} has no verified {member} checksum")
    path = Path(data_source.file(member))
    with path.open("rb") as handle:
        actual_sha256 = hashlib.file_digest(handle, "sha256").hexdigest()
    expected_sha256 = entry["sha256"]
    if actual_sha256 != expected_sha256:
        raise ValueError(
            f"{snapshot_id} {member} checksum mismatch: "
            f"expected {expected_sha256}, got {actual_sha256}"
        )
    return path, expected_sha256


def iter_disease_rows(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        xref_columns = [
            column for column in reader.fieldnames or []
            if column.endswith("_xref") and column != "primary_xref"
        ]
        for row_number, row in enumerate(reader, start=2):
            if None in row or any(value is None for value in row.values()):
                raise ValueError(f"{path}:{row_number} has malformed TSV fields")
            canonical_id = row["ncats_disease_id"].strip()
            if not canonical_id.startswith("IFXDisease:") or not canonical_id.split(":", 1)[1]:
                raise ValueError(f"{path}:{row_number} has invalid IFXDisease ID {canonical_id!r}")
            primary_xref = row["primary_xref"].strip()
            if not _is_curie(primary_xref):
                raise ValueError(f"{path}:{row_number} has invalid primary_xref {primary_xref!r}")
            exact_aliases = set()
            for column in xref_columns:
                for token in row[column].split("|"):
                    alias = token.strip()
                    if not alias:
                        continue
                    if not _is_curie(alias):
                        raise ValueError(f"{path}:{row_number} has invalid {column} alias {alias!r}")
                    exact_aliases.add(alias)
            if primary_xref not in exact_aliases:
                raise ValueError(
                    f"{path}:{row_number} primary_xref {primary_xref!r} is not an exact alias"
                )
            yield row_number, row, exact_aliases


def _is_curie(value: str) -> bool:
    if ":" not in value:
        return False
    prefix, identifier = value.split(":", 1)
    return bool(prefix.strip() and identifier.strip())
