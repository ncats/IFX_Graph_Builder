"""Read a completed harmonization stage and its effective source records."""
from collections import defaultdict

from src.core.curations import is_record_property_type, replay_curation_snapshot
from src.core.graph_build_identity import source_build_fingerprint
from src.core.record_property_curations import project_record_decisions


COLLECTIONS = (
    "MetaboliteIdentifier", "GeneIdentifier", "ProteinIdentifier", "PathwayIdentifier",
    "MetabolitePathwayEdge", "GenePathwayEdge", "ProteinPathwayEdge",
    "HmdbMetaboliteProteinAssociationEdge", "MetaboliteClassificationTerm",
    "MetaboliteClassificationEdge", "MetaboliteClassificationParentEdge",
    "HmdbOntologyTerm", "HmdbOntologyParentEdge", "HmdbMetaboliteOntologyEdge",
    "RheaReaction", "RheaMetaboliteReactionEdge", "RheaProteinReactionEdge",
    "RheaReactionClass", "RheaReactionClassParentEdge", "RheaReactionReactionClassEdge",
    "ChemicalEntity", "ChebiChemicalEntityMetaboliteIdentifierEdge",
    "HarmonizationStageActiveIdentifierChunk", "HarmonizedMetabolite",
    "HarmonizedMetaboliteMemberEdge", "metadata_store",
)


class StageReader:
    def __init__(self, db, stage_id, storage):
        self.db = db
        self.stage_id = stage_id.removeprefix("HarmonizationStage:")
        self.stage = db.collection("HarmonizationStage").get(self.stage_id)
        if not self.stage or self.stage.get("status") != "complete":
            raise ValueError(f"Stage {self.stage_id} is missing or incomplete")
        expected = self.stage.get("summary", {}).get("source_collection_revisions")
        if not expected:
            raise ValueError("Stage has no source revision record; rebuild the stage before export")
        self.revisions = {c: db.collection(c).revision() for c in set(COLLECTIONS) | set(expected)}
        changed = [c for c, revision in expected.items() if self.revisions[c] != revision]
        if changed:
            raise ValueError(f"Stage source evidence has changed: {', '.join(changed)}. Rebuild the stage.")
        self.metadata = db.collection("metadata_store").get("etl_metadata")["value"]
        recorded_build = self.stage["summary"].get("source_graph_fingerprint")
        if not recorded_build or source_build_fingerprint(self.metadata) != recorded_build:
            raise ValueError("Graph build metadata differs from the selected stage; rebuild the stage before export")
        if not self.metadata.get("registry_datasets"):
            raise ValueError("Graph build has no Registry input metadata")
        self.schemas = db.collection("metadata_store").get("collection_schemas")["collections"]
        missing_rhea_source = list(db.aql.execute(
            'FOR e IN RheaProteinReactionEdge FILTER !IS_STRING(e.source_id) OR e.source_id == "" '
            'LIMIT 1 RETURN e._key', max_runtime=120,
        ))
        if missing_rhea_source:
            raise ValueError('Rhea protein-reaction edges lack original source_id; refresh or repair the Rhea ingest before export')
        self.decisions = defaultdict(lambda: defaultdict(list))
        recorded = self.stage.get("curation_snapshots") or {}
        for kind in self.stage.get("curation_types") or []:
            if not recorded.get(kind):
                raise ValueError(f"Stage is missing recorded curation snapshot: {kind}")
        for kind, metadata in recorded.items():
            if not metadata:
                continue
            snapshot = replay_curation_snapshot(storage, metadata)
            if is_record_property_type(kind):
                selected_set = metadata.get("selected_curation_set")
                selected_model = metadata.get("selected_model_type")
                for decision in snapshot.active_record_property_decisions:
                    if selected_set and decision.target.get("curation_set") != selected_set:
                        continue
                    model, identifier = decision.target["model_type"], decision.target["id"]
                    if selected_model and model != selected_model:
                        continue
                    self.decisions[model][identifier].append(decision)
        # Validate all recorded corrections, including suppressed records and ChEBI
        # annotations used upstream to establish the persisted metabolite groups.
        for model, by_id in self.decisions.items():
            found = set()
            for doc in self.db.aql.execute(
                "FOR d IN @@c FILTER d.id IN @ids RETURN d",
                bind_vars={"@c": model, "ids": list(by_id)}, max_runtime=120,
            ):
                self.effective(model, doc)
                found.add(doc["id"])
            if found != set(by_id):
                raise ValueError(f"Missing curation targets in {model}: {sorted(set(by_id) - found)[:5]}")

    def effective(self, collection, doc):
        decisions = self.decisions[collection].get(doc.get("id"), [])
        if not decisions:
            return doc
        return project_record_decisions(
            doc, decisions, self.schemas[collection]["fields"], model_type=collection,
        )[0]

    def records(self, collection):
        cursor = self.db.aql.execute(
            "FOR d IN @@c SORT d._key RETURN d", bind_vars={"@c": collection},
            batch_size=5000, stream=True, max_runtime=3600,
        )
        try:
            for doc in cursor:
                yield self.effective(collection, doc)
        finally:
            cursor.close(ignore_missing=True)

    def groups(self):
        active = set()
        for chunk in self.db.aql.execute(
            "FOR d IN HarmonizationStageActiveIdentifierChunk FILTER d.stage_key == @s RETURN d.identifier_ids",
            bind_vars={"s": self.stage_id}, max_runtime=120,
        ):
            if active.intersection(chunk) or len(set(chunk)) != len(chunk):
                raise ValueError("Duplicate active stage identifiers")
            active.update(chunk)
        groups, seen = defaultdict(list), set()
        for row in self.db.aql.execute(
            "FOR d IN HarmonizedMetaboliteMemberEdge FILTER d.stage_key == @s "
            "RETURN KEEP(d, 'harmonized_metabolite_id', 'member_id')",
            bind_vars={"s": self.stage_id}, batch_size=10000, max_runtime=120,
        ):
            member = row["member_id"]
            if member not in active or member in seen:
                raise ValueError(f"Invalid stage group membership: {member}")
            groups[row["harmonized_metabolite_id"]].append(member)
            seen.add(member)
        expected_groups = {d["id"]: d["size"] for d in self.db.aql.execute(
            "FOR d IN HarmonizedMetabolite FILTER d.stage_key == @s RETURN KEEP(d, 'id', 'size')",
            bind_vars={"s": self.stage_id}, max_runtime=120,
        )}
        if {k: len(v) for k, v in groups.items()} != expected_groups:
            raise ValueError("Stage group documents and full membership edges disagree")
        singletons = active - seen
        summary = self.stage["summary"]
        observed = {"active_identifier_count": len(active), "non_singleton_clique_count": len(groups),
                    "singleton_identifier_count": len(singletons), "clique_count": len(groups) + len(singletons)}
        if any(summary.get(k) != v for k, v in observed.items()):
            raise ValueError(f"Stage membership totals disagree with its summary: {observed}")
        return sorted([tuple(sorted(v)) for v in groups.values()] + [(i,) for i in singletons])

    def verify_unchanged(self):
        current = self.db.collection("HarmonizationStage").get(self.stage_id)
        if not current or current["_rev"] != self.stage["_rev"]:
            raise ValueError("Stage changed during export; rerun against a stable stage")
        changed = [c for c, revision in self.revisions.items() if self.db.collection(c).revision() != revision]
        if changed:
            raise ValueError(f"Graph changed during export: {', '.join(changed)}")
