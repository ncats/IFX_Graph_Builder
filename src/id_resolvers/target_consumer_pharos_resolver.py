"""Pharos lookup projection over the TARGETS consumer contract.

Pharos uses Protein.xref as a lookup surface.  Its routed gene, transcript,
and isoform IDs therefore join the exact mapping IDs on canonical proteins.
The general Target Graph mapping resolvers remain identity-only.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from src.id_resolvers.sqlite_cache_resolver import SqliteCacheResolver
from src.id_resolvers.target_consumer_resolver import (
    _MATCH_TYPES,
    _REQUIRED_COLUMNS,
    annotate_xref_match_types,
    mapping_curie,
)
from src.interfaces.id_resolver import IdMatch
from src.models.node import EquivalentId
from src.models.protein import Protein


_RELATION_KINDS = {
    "isoform_of_canonical_protein",
    "transcript_to_protein",
    "gene_to_protein_resolved_parent",
    "gene_to_transcript",
    "gene_to_canonical_protein",
    "uniprot_entry_to_displayed_protein",
}
_RELATION_COLUMNS = {"subject_id", "object_id", "relation_kind", "status"}
_PROJECTION_VERSION = "pharos-canonical-selected-paths-v5-typed-lookups"
_VERSIONED_LOOKUP = re.compile(r"^((?:ENSEMBL|RefSeq):[^:]+)\.([0-9]+)$")


class TargetConsumerPharosCanonicalResolver(SqliteCacheResolver):
    """Map all selected target routes and their aliases onto canonical proteins."""

    name = "TARGETS Pharos canonical projection"

    def __init__(self, data_source, allow_unmatched_prefixes=None, **kwargs):
        self.allow_unmatched_prefixes = frozenset(allow_unmatched_prefixes or ())
        self.snapshot_id = data_source.snapshot_id
        self.mapping_path = Path(data_source.file("app_graph/target_mappings.tsv"))
        self.relationship_path = Path(data_source.file("app_graph/target_relationships.tsv"))
        self.protein_path = Path(data_source.file("app_graph/target_protein_payload.tsv"))
        manifest_files = {entry["path"]: entry for entry in data_source.manifest.get("files", [])}
        self.checksums = []
        for name in (
            "app_graph/target_mappings.tsv",
            "app_graph/target_relationships.tsv",
            "app_graph/target_protein_payload.tsv",
        ):
            checksum = manifest_files.get(name, {}).get("sha256")
            if not checksum:
                raise ValueError(f"TARGETS snapshot has no verified checksum for {name}")
            self.checksums.append(checksum)
        for path, required in (
            (self.mapping_path, _REQUIRED_COLUMNS),
            (self.relationship_path, _RELATION_COLUMNS),
            (self.protein_path, {"ifx_id", "is_canonical"}),
        ):
            with path.open(newline="", encoding="utf-8") as handle:
                columns = set(next(csv.reader(handle, delimiter="\t")))
            missing = required - columns
            if missing:
                raise ValueError(f"TARGETS file {path} lacks columns: {sorted(missing)}")
        with self.mapping_path.open(newline="", encoding="utf-8") as handle:
            self.has_match_type = "match_type" in next(csv.reader(handle, delimiter="\t"))
        if tuple(int(part) for part in data_source.version.split(".")) >= (2, 5, 2) and not self.has_match_type:
            raise ValueError("TARGETS 2.5.2+ mapping file lacks required match_type column")
        super().__init__(canonical_class=Protein, **kwargs)

    def _eligible_unmatched(self, entity_map):
        prefixes = getattr(self, "allow_unmatched_prefixes", frozenset())
        if not prefixes:
            return entity_map
        filtered = dict(entity_map)
        filtered[self.MatchKeys.unmatched] = {
            source_id: node
            for source_id, node in entity_map[self.MatchKeys.unmatched].items()
            if source_id.split(":", 1)[0] in prefixes
        }
        return filtered

    def parse_entity_map(self, entity_map):
        return super().parse_entity_map(self._eligible_unmatched(entity_map))

    def parse_flat_node_list_from_map(self, entity_map):
        return super().parse_flat_node_list_from_map(self._eligible_unmatched(entity_map))

    def get_version_info(self) -> str:
        return ":".join((self.snapshot_id, _PROJECTION_VERSION, *self.checksums))

    def matching_ids(self):
        raise NotImplementedError("Pharos routes are built from joined contract tables")

    def _xrefs_for_match(self, match, merged_xrefs):
        # A lookup ID may fan out to several canonical targets. Each receives
        # only its own routed IDs, never the other target's lookup surface.
        identifiers = set(match.equivalent_ids)
        return list({
            EquivalentId.parse(value)
            for value in identifiers
            if not ((version := _VERSIONED_LOOKUP.fullmatch(value)) and version.group(1) in identifiers)
        })

    def get_merged_map(self, entries, id_map, allow_retype=False):
        result = super().get_merged_map(entries, id_map, allow_retype=allow_retype)
        return annotate_xref_match_types(self, result) if getattr(self, "has_match_type", False) else result

    def resolve_internal(self, input_nodes):
        # Published exact versions retain priority. A different source version
        # can still resolve through the published version-free lookup key.
        matches = super().resolve_internal(input_nodes)
        fallback = {}
        for node in input_nodes:
            if node.id in matches:
                continue
            version = _VERSIONED_LOOKUP.fullmatch(node.id)
            if version:
                fallback[node.id] = version.group(1)
        if fallback:
            base_matches = super().resolve_internal([Protein(id=base) for base in set(fallback.values())])
            for input_id, base in fallback.items():
                matches_for_base = base_matches.get(base, [])
                if matches_for_base:
                    matches[input_id] = [
                        IdMatch(
                            input=input_id,
                            match=result.match,
                            equivalent_ids=result.equivalent_ids,
                            context=[*(result.context or []), "version_free_fallback"],
                        )
                        for result in matches_for_base
                    ]
        return matches

    def populate_lookup_db(self):
        db = self.connection
        db.execute("CREATE TEMP TABLE canonical (id TEXT PRIMARY KEY)")
        canonical = []
        with self.protein_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["is_canonical"] == "True":
                    canonical.append((row["ifx_id"],))
        db.executemany("INSERT INTO canonical VALUES (?)", canonical)
        if not canonical:
            raise ValueError("TARGETS protein payload contains no canonical proteins")

        db.execute("CREATE TEMP TABLE selected_relation (subject TEXT, object TEXT, kind TEXT)")
        batch = []
        with self.relationship_path.open(newline="", encoding="utf-8") as handle:
            for row in csv.DictReader(handle, delimiter="\t"):
                if row["status"] != "selected" or row["relation_kind"] not in _RELATION_KINDS:
                    continue
                batch.append((row["subject_id"], row["object_id"], row["relation_kind"]))
                if len(batch) >= 10000:
                    db.executemany("INSERT INTO selected_relation VALUES (?, ?, ?)", batch)
                    batch.clear()
        if batch:
            db.executemany("INSERT INTO selected_relation VALUES (?, ?, ?)", batch)
        db.execute("CREATE INDEX selected_relation_kind_subject ON selected_relation (kind, subject)")
        db.execute("CREATE INDEX selected_relation_kind_object ON selected_relation (kind, object)")

        db.execute("CREATE TEMP TABLE route (entity TEXT, canonical TEXT, path TEXT)")
        db.execute("CREATE UNIQUE INDEX route_pair ON route (entity, canonical, path)")
        db.execute("INSERT INTO route SELECT id, id, 'canonical' FROM canonical")
        db.execute("""
            INSERT OR IGNORE INTO route
            SELECT r.subject, c.id, 'isoform'
            FROM selected_relation r JOIN canonical c ON c.id = r.object
            WHERE r.kind = 'isoform_of_canonical_protein'
        """)
        db.execute("CREATE INDEX route_entity ON route (entity)")
        db.execute("""
            INSERT OR IGNORE INTO route
            SELECT r.subject, p.canonical, 'transcript_to_protein'
            FROM selected_relation r JOIN route p ON p.entity = r.object
            WHERE r.kind = 'transcript_to_protein'
        """)
        db.execute("""
            INSERT OR IGNORE INTO route
            SELECT r.subject, p.canonical, 'gene_to_protein_resolved_parent'
            FROM selected_relation r JOIN route p ON p.entity = r.object
            WHERE r.kind = 'gene_to_protein_resolved_parent'
              AND p.path IN ('canonical', 'isoform')
        """)
        db.execute("""
            INSERT OR IGNORE INTO route
            SELECT r.subject, t.canonical, 'gene_via_transcript'
            FROM selected_relation r JOIN route t ON t.entity = r.object
            WHERE r.kind = 'gene_to_transcript' AND t.path = 'transcript_to_protein'
        """)
        db.execute("""
            INSERT OR IGNORE INTO route
            SELECT r.subject, c.id, 'gene_to_canonical_protein'
            FROM selected_relation r JOIN canonical c ON c.id = r.object
            WHERE r.kind = 'gene_to_canonical_protein'
        """)
        db.execute("""
            INSERT OR IGNORE INTO route
            SELECT 'UniProtKB:' || r.subject, p.canonical, 'uniprot_entry_display'
            FROM selected_relation r JOIN route p ON p.entity = r.object
            WHERE r.kind = 'uniprot_entry_to_displayed_protein'
              AND p.path IN ('canonical', 'isoform')
        """)

        db.execute("CREATE TEMP TABLE mapping_alias (entity TEXT, alias TEXT, namespace TEXT)")
        batch = []
        with self.mapping_path.open(newline="", encoding="utf-8") as handle:
            for number, row in enumerate(csv.DictReader(handle, delimiter="\t"), start=2):
                if row["entity_type"] not in {"gene", "transcript", "protein"}:
                    raise ValueError(f"TARGETS mapping row {number} has unknown entity type {row['entity_type']!r}")
                if row["confidence"] not in {"high", "medium", "low"}:
                    raise ValueError(f"TARGETS mapping row {number} has invalid confidence {row['confidence']!r}")
                match_type = row["match_type"] if self.has_match_type else row["namespace"]
                if self.has_match_type and match_type not in _MATCH_TYPES:
                    raise ValueError(f"TARGETS mapping row {number} has unknown match_type {match_type!r}")
                for column in ("equivalent_id", "equivalent_id_base"):
                    batch.append((row["ifx_id"], mapping_curie(row["namespace"], row[column]), match_type))
                if len(batch) >= 10000:
                    db.executemany("INSERT INTO mapping_alias VALUES (?, ?, ?)", batch)
                    batch.clear()
        if batch:
            db.executemany("INSERT INTO mapping_alias VALUES (?, ?, ?)", batch)
        db.execute("CREATE INDEX mapping_alias_entity ON mapping_alias (entity)")
        db.execute("CREATE UNIQUE INDEX match_pair_type ON matches (id, match, type)")
        # Several selected paths may reach the same canonical protein. The
        # lookup metadata records that an IFX ID was routed; relationship rows
        # remain the authority for which paths support that route.
        db.execute("""
            INSERT OR IGNORE INTO matches
            SELECT canonical, entity, CASE WHEN path = 'canonical' THEN 'canonical' ELSE 'routed' END
            FROM route
        """)
        db.execute("""
            INSERT OR IGNORE INTO matches
            SELECT route.canonical, mapping_alias.alias, mapping_alias.namespace
            FROM route JOIN mapping_alias ON mapping_alias.entity = route.entity
        """)
        db.commit()
        self.store_file_metadata()
