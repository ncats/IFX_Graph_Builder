"""Export RaMP lookup tables from a completed harmonization stage."""
import argparse
from datetime import datetime, timezone
from pathlib import Path
import sys
import time


def _duration(seconds):
    tenths = round(seconds * 10)
    hours, remainder = divmod(tenths, 36_000)
    minutes, remainder = divmod(remainder, 600)
    return f"{hours:02d}:{minutes:02d}:{remainder / 10:04.1f}"


class BuildTimer:
    """Report completed phase and whole-build wall time using a monotonic clock."""

    def __init__(self):
        self.started = time.perf_counter()
        self.phase_started = self.started
        self.phase = None

    def mark(self, message):
        now = time.perf_counter()
        if self.phase is not None:
            print(f"Finished {self.phase}: {_duration(now - self.phase_started)}", flush=True)
        print(f"{message} (elapsed {_duration(now - self.started)})", flush=True)
        self.phase = message
        self.phase_started = now

    def finish(self, *, failed=False):
        now = time.perf_counter()
        if self.phase is not None:
            phase_label = "Stopped during" if failed else "Finished"
            print(f"{phase_label} {self.phase}: {_duration(now - self.phase_started)}", flush=True)
        label = "Elapsed before failure" if failed else "Total build time"
        print(f"{label}: {_duration(now - self.started)}", flush=True)


class SourceOnlyWriter:
    """Run the established projection while retaining reviewed diagnostic tables."""
    def __init__(self, writer):
        self.writer = writer
        self.db = writer.db
        self.columns = writer.columns

    def add(self, table, **row):
        if table in ("analyte", "source", "analytesynonym", "pathway", "ontology",
                     "analytehaspathway", "analytehasontology"):
            self.writer.add(table, **row)

    def flush(self, table):
        self.writer.flush(table)


def export_sqlite(reader, output, *, gene_identity, protein_annotations=None, release_version="unreleased", overwrite=False, source_only=False, progress=print):
    import yaml
    from src.use_cases.ramp.sqlite_projection import Projection, write_versions
    from src.use_cases.ramp.sqlite_writer import SQLiteWriter
    from src.use_cases.ramp.sqlite_gene_identity import POLICY

    gene_identity.validate(reader.metadata)
    policy = yaml.safe_load(Path(__file__).with_name("sqlite_ontology_policy.yaml").read_text())
    writer = SQLiteWriter(output, overwrite=overwrite)
    timestamp = datetime.now(timezone.utc).isoformat()
    try:
        writer.initialize()
        projection_writer = SourceOnlyWriter(writer) if source_only else writer
        projection = Projection(reader, projection_writer, policy["ontology"]["denylist"], gene_identity, progress,
                                protein_annotations=protein_annotations)
        projection.run(source_only=source_only)
        if source_only:
            writer.retain_tables({"analyte", "source", "analytesynonym", "pathway", "ontology",
                                  "analytehaspathway", "analytehasontology"})
        else:
            write_versions(writer, reader, release_version, timestamp,
                           kegg_via_hmdb='kegg' in projection.pathway_sources.values())
        progress("Validating SQLite and checking source revisions")
        manifest = {
            "format_version": 1, "scope": "lookup_table_diagnostic" if source_only else "base_tables",
            "status": "lookup_table_diagnostic" if source_only else "post_processing_pending",
            "release_ready": False, "release_version": release_version, "created_at": timestamp,
            "stage_id": reader.stage["id"], "stage_revision": reader.stage.get("_rev"),
            "stage_summary": reader.stage["summary"],
            "curation_snapshots": reader.stage.get("curation_snapshots", {}),
            "graph_build": reader.metadata, "source_revisions": reader.revisions,
            "id_policy": "fresh deterministic IDs; no cross-release continuity",
            "gene_policy": POLICY,
            "gene_resolution": gene_identity.manifest(),
            "protein_annotations": projection.protein_annotations.manifest(),
            "association_source_id_policy": "Preserve edge source_id / detail source_id, gene_id, protein_id; register against endpoint RAMP IDs without regrouping; missing evidence fails export",
            "hmdb_status": projection.hmdb_status_manifest,
            "chebi_chemistry_fallback": projection.chebi_chemistry_fallback,
            "pathway_attribution_policy": "HMDB-supplied KEGG pathways use kegg in pathway.type and analytehaspathway.pathwaySource; input provenance remains HMDB",
            "omitted_tables": ["reaction_protein2met"],
            "pending_tables": ["entity_status_info", "pathway_similarity", "pathway_duplicates"],
            "pending_fields": ["source.pathwayCount", "ontology.metCount", "db_version.*intersects*",
                               "reaction.has_human_prot", "reaction.only_human_mets",
                               "reaction2met.is_cofactor (when absent from source)",
                               "reaction2protein.is_reviewed (when absent from source)"],
            "pending_value_policy": "NULL where nullable, -1 where legacy schema requires an integer",
            "excluded_metabolite_edges": dict(projection.skipped),
        }
        if source_only:
            manifest["included_tables"] = ["analyte", "source", "analytesynonym", "pathway",
                                           "ontology", "analytehaspathway", "analytehasontology"]
            manifest["pending_tables"] = []
            manifest["pending_fields"] = ["source.pathwayCount", "ontology.metCount"]
        return writer.finish(manifest, reader.verify_unchanged)
    finally:
        writer.close()


def parser():
    result = argparse.ArgumentParser(
        description="Export RaMP lookup tables from a harmonization stage (full base export optional).",
        epilog="Example: python -m src.use_cases.ramp.build_sqlite --stage-id stage-07-9cfef0c333530e83 --output ramp-base.sqlite",
    )
    result.add_argument("--stage-id", required=True, help="Exact graph stage key or HarmonizationStage ID")
    result.add_argument("--output", required=True, type=Path, help="SQLite output path")
    result.add_argument('--overwrite', action='store_true', help='Replace an existing output only after the new database passes validation')
    result.add_argument('--full-base', action='store_true', help='Build the full base-table export instead of the default source-table diagnostic')
    result.add_argument("--release-version", default="unreleased", help="Output version label (default: unreleased)")
    result.add_argument("--database", default="metabolite_harmonization")
    result.add_argument("--graph-credentials", type=Path, default=Path("src/use_cases/secrets/ifxdev_arangodb.yaml"), help="Arango credential YAML file")
    result.add_argument("--curation-credentials", type=Path, default=Path("src/use_cases/secrets/aws_ifx_registry.yaml"), help="Object-storage credential YAML file for recorded curation batches")
    result.add_argument("--registry-credentials", type=Path, default=Path("src/use_cases/secrets/aws_ifx_registry.yaml"), help="Registry credentials for the graph's pinned UniProt input")
    result.add_argument("--registry-cache-dir", type=Path, default=Path('/var/tmp/ifx-registry-cache'))
    result.add_argument('--uniprot-file', choices=['uniprot-human-reviewed.json.gz', 'uniprot-human.json.gz'],
                        default='uniprot-human.json.gz', help='Resolver input within the recorded snapshot (default: all human, including unreviewed)')
    return result


def main(argv=None):
    args = parser().parse_args(argv)
    from src.use_cases.ramp.sqlite_writer import SQLiteWriter
    try:
        SQLiteWriter.validate_output(args.output, overwrite=args.overwrite)
    except (OSError, ValueError) as exc:
        parser().error(str(exc))
    if not args.stage_id.strip() or not args.release_version.strip():
        parser().error("Stage ID and release version must not be blank")
    timer = BuildTimer()
    try:
        import yaml
        from src.infrastructure.object_storage import load_object_storage_credentials, object_storage_from_credentials
        from src.shared.arango_adapter import ArangoAdapter
        from src.shared.db_credentials import DBCredentials
        from src.use_cases.ramp.sqlite_stage import StageReader
        from src.use_cases.ramp.sqlite_gene_identity import GeneIdentity
        from src.core.registry_integration import RegistryIntegration

        scope = "base tables; post-processing pending" if args.full_base else "lookup-table diagnostic"
        print(f"Stage: {args.stage_id}\nGraph: {args.database}\nOutput: {args.output}\nScope: {scope}", flush=True)
        timer.mark("Connecting to graph and storage")
        credentials = DBCredentials.from_yaml(yaml.safe_load(args.graph_credentials.read_text()))
        db = ArangoAdapter(credentials, args.database).get_db()
        storage = object_storage_from_credentials(load_object_storage_credentials(args.curation_credentials))
        timer.mark("Reading stage and replaying recorded corrections")
        reader = StageReader(db, args.stage_id, storage)
        timer.mark('Building gene/protein resolver from the recorded UniProt snapshot')
        registry = RegistryIntegration.connect({'credentials': args.registry_credentials,
                                               'cache_dir': args.registry_cache_dir})
        identity = GeneIdentity.from_stage(reader.metadata, registry, file_name=args.uniprot_file)
        from src.use_cases.ramp.sqlite_protein_annotations import ProteinAnnotations
        annotations = ProteinAnnotations.from_file(identity.input_file, identity.provenance)
        result = export_sqlite(reader, args.output, gene_identity=identity, protein_annotations=annotations, release_version=args.release_version,
                               overwrite=args.overwrite, source_only=not args.full_base, progress=timer.mark)
    except Exception as exc:
        timer.finish(failed=True)
        print(f"RaMP SQLite export failed: {exc}", file=sys.stderr)
        return 1
    timer.finish()
    print(f"Created {args.output}\n" + ("Base tables complete; post-processing pending." if args.full_base
                                   else "Lookup-table diagnostic complete; not a release database."))
    if args.full_base:
        print("Manifest: SELECT value FROM ramp_export_metadata WHERE key = 'manifest';")
    for table, count in sorted(result["row_counts"].items()):
        print(f"  {table}: {count:,}")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
