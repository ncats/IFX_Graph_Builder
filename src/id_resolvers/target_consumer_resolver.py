"""Identity-only TARGETS consumer-contract resolvers."""

from __future__ import annotations

import csv
from collections import defaultdict
from pathlib import Path

from src.id_resolvers.sqlite_cache_resolver import MatchingPair, SqliteCacheResolver


_NAMESPACE_PREFIX = {
    "ENSEMBL": "ENSEMBL",
    "NCBIGene": "NCBIGene",
    "HGNC": "HGNC",
    "OMIM": "OMIM",
    "RefSeq": "RefSeq",
    "UniProtKB": "UniProtKB",
    "miRBase": "miRBase",
    "IMGT": "IMGT",
    "VEGA": "Vega",
    "Orphanet": "orphanet",
    "Symbol": "Symbol",
}
_MATCH_TYPES = {"exact", "symbol", "previous_symbol", "secondary_accession", "uniprot_entry"}
_REQUIRED_COLUMNS = {
    "ifx_id", "entity_type", "equivalent_id", "equivalent_id_base",
    "namespace", "source", "source_version", "confidence",
}
_ENTITY_PREFIX = {"gene": "IFXGene:", "transcript": "IFXTranscript:", "protein": "IFXProtein:"}


def mapping_curie(namespace: str, identifier: str) -> str:
    """Express one published lookup ID in ODIN's existing prefix syntax."""
    try:
        prefix = _NAMESPACE_PREFIX[namespace]
    except KeyError as exc:
        raise ValueError(f"Unknown TARGETS mapping namespace {namespace!r}") from exc
    value = identifier.strip()
    if not value:
        raise ValueError(f"Empty TARGETS {namespace} equivalent ID")
    if value.startswith(f"{prefix}:"):
        return value
    if namespace == "OMIM" and value.startswith("MIM:"):
        return f"OMIM:{value[4:]}"
    if namespace == "IMGT" and value.startswith("IMGT/GENE-DB:"):
        return f"IMGT:{value[len('IMGT/GENE-DB:'):]}"
    if namespace == "Symbol":
        # Source-defined symbols can themselves contain colons.
        return f"Symbol:{value}"
    if ":" in value:
        raise ValueError(f"Unexpected prefix in TARGETS {namespace} ID {value!r}")
    return f"{prefix}:{value}"


def annotate_xref_match_types(resolver, entity_map):
    """Group the persisted lookup surface by the cache's published match types."""
    matched = entity_map[resolver.MatchKeys.matched]
    newborns = entity_map[resolver.MatchKeys.newborns]
    nodes = [*matched.values(), *(node for group in newborns.values() for node in group)]
    nodes_by_id = defaultdict(list)
    for node in nodes:
        nodes_by_id[node.id].append(node)
    ids = sorted(nodes_by_id)
    with resolver._sqlite_lock:
        cursor = resolver.connection.cursor()
        try:
            for start in range(0, len(ids), 250):
                chunk = ids[start:start + 250]
                by_id = defaultdict(lambda: defaultdict(set))
                cursor.execute(
                    f"SELECT id, match, type FROM matches WHERE id IN ({','.join('?' * len(chunk))})",
                    chunk,
                )
                for canonical_id, lookup_id, match_type in cursor:
                    by_id[canonical_id][match_type].add(lookup_id)
                for canonical_id in chunk:
                    for node in nodes_by_id[canonical_id]:
                        xrefs = {xref.id_str() for xref in node.xref or []}
                        groups = {
                            match_type: sorted(values & xrefs)
                            for match_type, values in by_id[canonical_id].items()
                            if values & xrefs
                        }
                        if set().union(*groups.values()) != xrefs:
                            raise ValueError(f"TARGETS lookup metadata does not cover all xrefs for {node.id}")
                        node.xref_match_types = dict(sorted(groups.items()))
        finally:
            cursor.close()
    return entity_map


class TargetConsumerResolver(SqliteCacheResolver):
    entity_type: str
    name = "TARGETS consumer mappings"

    def __init__(self, data_source, **kwargs):
        self.mapping_path = Path(data_source.file("app_graph/target_mappings.tsv"))
        self.snapshot_id = data_source.snapshot_id
        if not self.mapping_path.is_file():
            raise FileNotFoundError(self.mapping_path)
        with self.mapping_path.open(newline="", encoding="utf-8") as handle:
            columns = set(next(csv.reader(handle, delimiter="\t")))
        missing = _REQUIRED_COLUMNS - columns
        if missing:
            raise ValueError(f"TARGETS mapping file lacks columns: {sorted(missing)}")
        self.has_match_type = "match_type" in columns
        if tuple(int(part) for part in data_source.version.split(".")) >= (2, 5, 2) and not self.has_match_type:
            raise ValueError("TARGETS 2.5.2+ mapping file lacks required match_type column")
        manifest_files = data_source.manifest.get("files") or []
        match = next((entry for entry in manifest_files if entry.get("path") == "app_graph/target_mappings.tsv"), None)
        if not match or not match.get("sha256"):
            raise ValueError("TARGETS snapshot has no verified mapping checksum")
        self.mapping_sha256 = match["sha256"]
        super().__init__(**kwargs)

    def get_version_info(self) -> str:
        return f"{self.snapshot_id}:{self.mapping_sha256}:{self.entity_type}:typed-lookups-v2"

    def matching_ids(self):
        with self.mapping_path.open(newline="", encoding="utf-8") as handle:
            for row_number, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
                if row["entity_type"] != self.entity_type:
                    continue
                if not row["ifx_id"].startswith(_ENTITY_PREFIX[self.entity_type]):
                    raise ValueError(f"TARGETS mapping row {row_number} has wrong IFX ID: {row['ifx_id']!r}")
                if row["confidence"] not in {"high", "medium", "low"}:
                    raise ValueError(f"TARGETS mapping row {row_number} has invalid confidence {row['confidence']!r}")
                match_type = row["match_type"] if self.has_match_type else row["namespace"]
                if self.has_match_type and match_type not in _MATCH_TYPES:
                    raise ValueError(f"TARGETS mapping row {row_number} has unknown match_type {match_type!r}")
                target = row["ifx_id"]
                yield MatchingPair(id=target, match=target, type="exact")
                for column in ("equivalent_id", "equivalent_id_base"):
                    yield MatchingPair(
                        id=target,
                        match=mapping_curie(row["namespace"], row[column]),
                        type=match_type,
                    )

    def get_merged_map(self, entries, id_map, allow_retype=False):
        result = super().get_merged_map(entries, id_map, allow_retype=allow_retype)
        return annotate_xref_match_types(self, result) if self.has_match_type else result

    def populate_lookup_db(self):
        """Stream the large mapping table into SQLite and reject ambiguous IDs."""
        self.connection.execute("CREATE UNIQUE INDEX IF NOT EXISTS mapping_pair_type ON matches (id, match, type)")
        batch = []
        for pair in self.matching_ids():
            batch.append((pair.id, pair.match, pair.type))
            if len(batch) >= 10000:
                self.connection.executemany("INSERT OR IGNORE INTO matches VALUES (?, ?, ?)", batch)
                batch.clear()
        if batch:
            self.connection.executemany("INSERT OR IGNORE INTO matches VALUES (?, ?, ?)", batch)
        type_conflict = self.connection.execute(
            "SELECT id, match, MIN(type), MAX(type) FROM matches "
            "GROUP BY id, match HAVING MIN(type) != MAX(type) LIMIT 1"
        ).fetchone()
        if type_conflict:
            raise ValueError(
                f"TARGETS {self.entity_type} lookup {type_conflict[1]} on {type_conflict[0]} "
                f"has conflicting match types {type_conflict[2]} and {type_conflict[3]}"
            )
        conflict = self.connection.execute(
            "SELECT match, MIN(id), MAX(id) FROM matches GROUP BY match "
            "HAVING MIN(id) != MAX(id) LIMIT 1"
        ).fetchone()
        if conflict:
            raise ValueError(
                f"TARGETS {self.entity_type} mapping conflict for {conflict[0]}: "
                f"{conflict[1]} versus {conflict[2]}"
            )
        self.connection.commit()
        self.store_file_metadata()


class TargetConsumerGeneResolver(TargetConsumerResolver):
    entity_type = "gene"


class TargetConsumerTranscriptResolver(TargetConsumerResolver):
    entity_type = "transcript"


class TargetConsumerProteinResolver(TargetConsumerResolver):
    entity_type = "protein"
