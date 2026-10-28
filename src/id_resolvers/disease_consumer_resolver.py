"""Exact disease identity from the pinned DISEASES consumer release."""

from __future__ import annotations

from itertools import chain

from src.id_resolvers.sqlite_cache_resolver import (
    DEFAULT_SQLITE_RESOLVER_CACHE_DIR,
    MatchingPair,
    SqliteCacheResolver,
)
from src.shared.disease_consumer_contract import disease_ids_file, iter_disease_rows


class DiseaseConsumerResolver(SqliteCacheResolver):
    name = "DISEASES consumer exact aliases"
    cache_schema_version = "v1"

    def __init__(self, data_source, cache_path: str | None = None, **kwargs):
        self.path, checksum = disease_ids_file(data_source)
        self.cache_path = cache_path
        self.version_key = f"{data_source.snapshot_id}:{checksum}:{self.cache_schema_version}"
        super().__init__(**kwargs)

    def cache_location(self):
        return self.cache_path or str(
            DEFAULT_SQLITE_RESOLVER_CACHE_DIR / f"{self.__class__.__name__}.sqlite"
        )

    def get_version_info(self) -> str:
        return self.version_key

    def matching_ids(self):
        for _row_number, row, exact_aliases in iter_disease_rows(self.path):
            canonical_id = row["ncats_disease_id"].strip()
            yield MatchingPair(id=canonical_id, match=canonical_id, type="canonical")
            for alias in sorted(exact_aliases):
                yield MatchingPair(id=canonical_id, match=alias, type="exact")

    def get_merged_map(self, entries, id_map, allow_retype=False):
        entity_map = super().get_merged_map(entries, id_map, allow_retype=allow_retype)
        matched = entity_map[self.MatchKeys.matched].values()
        newborns = (
            node
            for group in entity_map[self.MatchKeys.newborns].values()
            for node in group
        )
        for node in chain(matched, newborns):
            # The released *_xref columns are exact aliases. Use the resolver's
            # normalized xrefs so prefix spelling matches the stored graph.
            exact_ids = sorted({xref.id_str() for xref in node.xref or []})
            if not exact_ids:
                raise ValueError(f"DISEASES resolved {node.id} without exact xrefs")
            node.xref_match_types = {"exact": exact_ids}
        return entity_map

    def populate_lookup_db(self):
        self.connection.execute(
            "CREATE UNIQUE INDEX IF NOT EXISTS disease_pair ON matches (id, match)"
        )
        batch = []
        seen_ids = set()
        for row_number, row, exact_aliases in iter_disease_rows(self.path):
            canonical_id = row["ncats_disease_id"].strip()
            if canonical_id in seen_ids:
                raise ValueError(f"{self.path}:{row_number} duplicates {canonical_id}")
            seen_ids.add(canonical_id)
            batch.append((canonical_id, canonical_id, "canonical"))
            batch.extend((canonical_id, alias, "exact") for alias in exact_aliases)
            if len(batch) >= 10000:
                self.connection.executemany(
                    "INSERT OR IGNORE INTO matches VALUES (?, ?, ?)", batch
                )
                batch.clear()
        if not seen_ids:
            raise ValueError(f"{self.path} has no DISEASES rows")
        if batch:
            self.connection.executemany(
                "INSERT OR IGNORE INTO matches VALUES (?, ?, ?)", batch
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
        self.connection.commit()
        self.store_file_metadata()
