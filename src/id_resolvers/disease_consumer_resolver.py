"""Disease identity and typed xrefs from the pinned DISEASES release."""

from __future__ import annotations

import csv
from collections import defaultdict
from itertools import chain

from src.id_resolvers.sqlite_cache_resolver import (
    DEFAULT_SQLITE_RESOLVER_CACHE_DIR,
    MatchingPair,
    SqliteCacheResolver,
)
from src.interfaces.id_resolver import MultiMatchBehavior
from src.models.node import EquivalentId
from src.shared.disease_consumer_contract import (
    disease_ids_file,
    iter_disease_rows,
    validated_disease_snapshot_file,
)


_MATCH_TYPES = {"exact", "broad", "narrow", "related", "unclassified"}


class DiseaseConsumerResolver(SqliteCacheResolver):
    name = "DISEASES consumer exact aliases"
    cache_schema_version = "v5"

    def __init__(self, data_source, cache_path: str | None = None, **kwargs):
        self.path, checksum = disease_ids_file(data_source)
        self.edges_path, edges_checksum = validated_disease_snapshot_file(
            data_source, "app_graph/disease_xref_edges.tsv"
        )
        self.labels_path, labels_checksum = validated_disease_snapshot_file(
            data_source, "app_graph/xref_labels.tsv"
        )
        self.obsolete_ids = self._read_obsolete_ids()
        self.cache_path = cache_path
        self.version_key = (
            f"{data_source.snapshot_id}:{checksum}:{edges_checksum}:"
            f"{labels_checksum}:{self.cache_schema_version}"
        )
        super().__init__(**kwargs)

    def _read_obsolete_ids(self):
        obsolete = set()
        seen = set()
        with self.labels_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"xref_id", "is_obsolete"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(
                    f"{self.labels_path} lacks {sorted(required - set(reader.fieldnames or []))}"
                )
            for row_number, row in enumerate(reader, start=2):
                source_id = (row["xref_id"] or "").strip()
                status = (row["is_obsolete"] or "").strip().lower()
                if not source_id or ":" not in source_id or status not in {"true", "false"}:
                    raise ValueError(f"{self.labels_path}:{row_number} has invalid obsolete label")
                if source_id in seen:
                    raise ValueError(f"{self.labels_path}:{row_number} duplicates {source_id}")
                seen.add(source_id)
                if status == "true":
                    obsolete.add(source_id)
        if not seen:
            raise ValueError(f"{self.labels_path} has no xref labels")
        return obsolete

    def cache_location(self):
        return self.cache_path or str(
            DEFAULT_SQLITE_RESOLVER_CACHE_DIR / f"{self.__class__.__name__}.sqlite"
        )

    def get_version_info(self) -> str:
        return self.version_key

    def db_is_corrupt(self):
        return super().db_is_corrupt() or not self.table_exists("xref_annotations")

    def matching_ids(self):
        for row_number, row, exact_aliases in iter_disease_rows(self.path):
            canonical_id = row["ncats_disease_id"].strip()
            if row["primary_xref"].strip() in self.obsolete_ids:
                raise ValueError(f"{self.path}:{row_number} has obsolete primary_xref")
            yield MatchingPair(id=canonical_id, match=canonical_id, type="canonical")
            for alias in sorted(exact_aliases):
                if alias not in self.obsolete_ids:
                    yield MatchingPair(id=canonical_id, match=alias, type="exact")
        with self.edges_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if (row["match_type"] == "exact"
                    and row["xref_id"] not in self.obsolete_ids
                    and row["xref_is_obsolete"].strip().lower() == "false"):
                    yield MatchingPair(
                        id=row["ncats_disease_id"], match=row["xref_id"], type="exact"
                    )

    def resolve_internal(self, input_nodes):
        matches = super().resolve_internal(input_nodes)
        for input_id, candidates in matches.items():
            candidates.sort(key=lambda candidate: (candidate.match != input_id, candidate.match))
        return matches

    def parse_flat_node_list_from_map(self, entity_map):
        ambiguous = entity_map[self.MatchKeys.newborns]
        if self.multi_match_behavior == MultiMatchBehavior.Error and ambiguous:
            source_id = min(ambiguous)
            targets = [entity_map[self.MatchKeys.matched][source_id].id]
            targets.extend(node.id for node in ambiguous[source_id])
            raise ValueError(
                f"DISEASES identifier {source_id!r} resolves to multiple targets "
                f"{sorted(targets)}"
            )
        return super().parse_flat_node_list_from_map(entity_map)

    def get_merged_map(self, entries, id_map, allow_retype=False):
        entity_map = super().get_merged_map(entries, id_map, allow_retype=allow_retype)
        matched = entity_map[self.MatchKeys.matched].values()
        newborns = (
            node
            for group in entity_map[self.MatchKeys.newborns].values()
            for node in group
        )
        nodes = list(chain(matched, newborns))
        by_id = defaultdict(list)
        for node in nodes:
            by_id[node.id].append(node)
        if not by_id:
            return entity_map
        canonical_ids = sorted(by_id)
        rows = []
        with self._sqlite_lock:
            for start in range(0, len(canonical_ids), 900):
                chunk = canonical_ids[start:start + 900]
                placeholders = ",".join("?" for _ in chunk)
                rows.extend(self.connection.execute(
                    f"SELECT id, match, type FROM xref_annotations WHERE id IN ({placeholders})",
                    chunk,
                ).fetchall())
        annotations = defaultdict(lambda: defaultdict(set))
        for canonical_id, source_id, match_type in rows:
            annotations[canonical_id][match_type].add(source_id)
        for canonical_id, group in by_id.items():
            typed = annotations[canonical_id]
            if not typed.get("exact"):
                raise ValueError(f"DISEASES resolved {canonical_id} without exact xrefs")
            xrefs = sorted({source_id for values in typed.values() for source_id in values})
            for node in group:
                node.xref = [self._raw_xref(source_id) for source_id in xrefs]
                node.xref_match_types = {
                    match_type: sorted(values)
                    for match_type, values in sorted(typed.items())
                }
        return entity_map

    @staticmethod
    def _raw_xref(source_id):
        prefix, local_id = source_id.split(":", 1)
        return EquivalentId(id=local_id, type=prefix)

    def populate_lookup_db(self):
        self.connection.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS disease_pair ON matches (id, match)"
        )
        self.connection.execute(
            "CREATE TABLE xref_annotations (id TEXT, match TEXT, type TEXT, "
            "PRIMARY KEY (id, match))"
        )
        batch = []
        annotation_batch = []
        seen_ids = set()
        for row_number, row, exact_aliases in iter_disease_rows(self.path):
            canonical_id = row["ncats_disease_id"].strip()
            if row["primary_xref"].strip() in self.obsolete_ids:
                raise ValueError(
                    f"{self.path}:{row_number} has obsolete primary_xref {row['primary_xref']!r}"
                )
            if canonical_id in seen_ids:
                raise ValueError(f"{self.path}:{row_number} duplicates {canonical_id}")
            seen_ids.add(canonical_id)
            batch.append((canonical_id, canonical_id, "canonical"))
            current_aliases = exact_aliases - self.obsolete_ids
            batch.extend((canonical_id, alias, "exact") for alias in current_aliases)
            annotation_batch.append((canonical_id, canonical_id, "exact"))
            annotation_batch.extend(
                (canonical_id, EquivalentId.parse(alias).id_str(), "exact")
                for alias in current_aliases
            )
            if len(batch) >= 10000:
                self.connection.executemany(
                    "INSERT OR IGNORE INTO matches VALUES (?, ?, ?)", batch
                )
                batch.clear()
            if len(annotation_batch) >= 10000:
                self.connection.executemany(
                    "INSERT OR REPLACE INTO xref_annotations VALUES (?, ?, ?)", annotation_batch
                )
                annotation_batch.clear()
        if not seen_ids:
            raise ValueError(f"{self.path} has no DISEASES rows")
        if batch:
            self.connection.executemany(
                "INSERT OR IGNORE INTO matches VALUES (?, ?, ?)", batch
            )
        if annotation_batch:
            self.connection.executemany(
                "INSERT OR REPLACE INTO xref_annotations VALUES (?, ?, ?)", annotation_batch
            )
        conflict = self.connection.execute(
            "SELECT match, MIN(id), MAX(id) FROM matches GROUP BY match "
            "HAVING MIN(id) != MAX(id) LIMIT 1"
        ).fetchone()
        if conflict:
            raise ValueError(
                f"DISEASES exact alias {conflict[0]!r} belongs to both "
                f"{conflict[1]} and {conflict[2]}"
            )
        with self.edges_path.open(newline="", encoding="utf-8") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"ncats_disease_id", "xref_id", "match_type", "xref_is_obsolete"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"{self.edges_path} lacks {sorted(required - set(reader.fieldnames or []))}")
            edge_batch = []
            for row_number, row in enumerate(reader, start=2):
                canonical_id = row["ncats_disease_id"].strip()
                source_id = row["xref_id"].strip()
                match_type = row["match_type"].strip()
                edge_obsolete = row["xref_is_obsolete"].strip().lower()
                if canonical_id not in seen_ids or ":" not in source_id or match_type not in _MATCH_TYPES:
                    raise ValueError(f"{self.edges_path}:{row_number} has invalid disease xref edge")
                if edge_obsolete not in {"true", "false"}:
                    raise ValueError(f"{self.edges_path}:{row_number} has invalid xref_is_obsolete")
                if source_id in self.obsolete_ids or edge_obsolete == "true":
                    continue
                edge_batch.append((canonical_id, source_id, match_type))
                if len(edge_batch) >= 10000:
                    self._store_edge_annotations(edge_batch)
                    edge_batch.clear()
            if edge_batch:
                self._store_edge_annotations(edge_batch)
        self.store_file_metadata()

    def _store_edge_annotations(self, rows):
        self.connection.executemany(
            "INSERT INTO xref_annotations VALUES (?, ?, ?) "
            "ON CONFLICT(id, match) DO UPDATE SET type = excluded.type",
            rows,
        )
        self.connection.executemany(
            "INSERT OR IGNORE INTO matches VALUES (?, ?, ?)",
            ((canonical_id, source_id, "exact")
             for canonical_id, source_id, match_type in rows if match_type == "exact"),
        )
