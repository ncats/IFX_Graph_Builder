"""Build a complete RaMP SQLite from a validated harmonization stage."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import sqlite3
import sys
import tempfile
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
                     "analytehaspathway", "analytehasontology", "catalyzed", "metabolite_class",
                     "reaction", "reaction2met", "reaction2protein", "reaction_ec_class",
                     "chem_props", "version_info"):
            self.writer.add(table, **row)

    def flush(self, table):
        self.writer.flush(table)


def export_sqlite(reader, output, *, gene_identity, protein_annotations=None, release_version="unreleased", overwrite=False,
                  source_only=False, pathway_association_cutoff=25000,
                  include_lipidmaps_class_level4=False, source_ontology_policy='legacy', progress=print):
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
        projection = Projection(reader, projection_writer, policy["ontology"], gene_identity, progress,
                                protein_annotations=protein_annotations,
                                pathway_association_cutoff=pathway_association_cutoff,
                                include_lipidmaps_class_level4=include_lipidmaps_class_level4,
                                source_ontology_policy=source_ontology_policy)
        projection.run(source_only=source_only)
        if source_only:
            writer.retain_tables({"analyte", "source", "analytesynonym", "pathway", "ontology",
                                  "analytehaspathway", "analytehasontology", "catalyzed", "metabolite_class",
                                  "reaction", "reaction2met", "reaction2protein", "reaction_ec_class",
                                  "chem_props", "version_info", "ramp_export_metadata"})
        write_versions(writer, reader, release_version, timestamp,
                       kegg_via_hmdb='kegg' in projection.pathway_sources.values(),
                       include_db_version=not source_only)
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
            "reaction_cofactor_policy": "ChEBI CHEBI:23357 cofactor role, role-bearing chemicals, and is_a chemical descendants; classify reported ChEBI participant IDs, not RaMP groups",
            "cofactor_role_count": projection.cofactor_role_count,
            "cofactor_chemical_id_count": projection.cofactor_chemical_id_count,
            "cofactor_reaction_assertions": projection.cofactor_reaction_assertions,
            "reaction_protein_direction_policy": "UN only; directional protein evidence remains in source lookups",
            "non_un_protein_assertions_excluded": projection.non_un_protein_assertions_excluded,
            "hmdb_status": projection.hmdb_status_manifest,
            "chebi_chemistry_fallback": projection.chebi_chemistry_fallback,
            "pathway_attribution_policy": "HMDB-supplied KEGG pathways use kegg in pathway.type and analytehaspathway.pathwaySource; input provenance remains HMDB",
            "omitted_tables": ["reaction_protein2met"],
            "pending_tables": ["entity_status_info", "pathway_similarity", "pathway_duplicates"],
            "pending_fields": ["source.pathwayCount", "ontology.metCount", "db_version.*intersects*",
                               "reaction.has_human_prot", "reaction.only_human_mets"],
            "protein_review_status_policy": "Use source review status when available; -1 means the source did not report it",
            "pending_value_policy": "NULL where nullable, -1 where legacy schema requires an integer",
            "excluded_metabolite_edges": dict(projection.skipped),
            "pathway_association_cutoff": pathway_association_cutoff,
            "pathway_source_ids_at_cutoff": projection.pathway_source_ids_at_cutoff,
            "pathway_assertions_excluded": projection.pathway_assertions_excluded,
            "include_lipidmaps_class_level4": include_lipidmaps_class_level4,
            "lipidmaps_class_level4_assertions_excluded": projection.lipidmaps_class_level4_assertions_excluded,
            "source_ontology_policy": source_ontology_policy,
        }
        if source_only:
            manifest["included_tables"] = ["analyte", "source", "analytesynonym", "pathway",
                                           "ontology", "analytehaspathway", "analytehasontology",
                                           "catalyzed", "metabolite_class", "reaction", "reaction2met",
                                           "reaction2protein", "reaction_ec_class", "chem_props",
                                           "version_info", "ramp_export_metadata"]
            manifest["pending_tables"] = []
            manifest["pending_fields"] = ["source.pathwayCount", "ontology.metCount"]
        return writer.finish(manifest, reader.verify_unchanged)
    finally:
        writer.close()


def parser():
    result = argparse.ArgumentParser(
        description="Build and post-process a complete RaMP SQLite from a validated harmonization stage.",
        epilog="Example: python -m src.use_cases.ramp.build_sqlite --stage-id stage-07-9cfef0c333530e83 --output ramp-base.sqlite",
    )
    result.add_argument("--stage-id", required=True, help="Exact graph stage key or HarmonizationStage ID")
    result.add_argument("--output", required=True, type=Path, help="SQLite output path")
    result.add_argument('--overwrite', action='store_true', help='Replace an existing output only after the new database passes validation')
    scope = result.add_mutually_exclusive_group()
    scope.add_argument('--diagnostic', action='store_true',
                       help='Export lookup/association/reaction diagnostic tables only; incomplete for the R package')
    scope.add_argument('--base-only', action='store_true',
                       help='Export all base tables, leaving post-processing pending')
    scope.add_argument('--full-base', action='store_true',
                       help='Deprecated alias for the default complete build')
    result.add_argument('--refresh-report-from', type=Path, metavar='COMPARISON_JSON',
                        help='After a complete build, refresh an existing comparison using its recorded historical inputs')
    result.add_argument('--report-output', type=Path, metavar='HTML',
                        help='Comparison HTML path (default: COMPARISON_JSON with .html suffix)')
    result.add_argument('--pathway-association-cutoff', type=nonnegative_int, default=25000,
                        help='Omit all metabolite/gene pathway links for a reported source ID with this many or more links, per source; 0 disables (default: 25000)')
    result.add_argument('--include-lipidmaps-class-level4', action='store_true',
                        help='Include LipidMaps CLASS_LEVEL4 in metabolite_class (excluded by default to match legacy RaMP)')
    result.add_argument('--source-ontology-policy', choices=('legacy', 'expanded'), default='legacy',
                        help='Select HMDB Source ontology terms: legacy named terms (default) or all non-denied terms')
    result.add_argument("--release-version", default="unreleased", help="Output version label (default: unreleased)")
    result.add_argument("--database", default="metabolite_harmonization")
    result.add_argument("--graph-credentials", type=Path, default=Path("src/use_cases/secrets/ifxdev_arangodb.yaml"), help="Arango credential YAML file")
    result.add_argument("--curation-credentials", type=Path, default=Path("src/use_cases/secrets/aws_ifx_registry.yaml"), help="Object-storage credential YAML file for recorded curation batches")
    result.add_argument("--registry-credentials", type=Path, default=Path("src/use_cases/secrets/aws_ifx_registry.yaml"), help="Registry credentials for the graph's pinned UniProt input")
    result.add_argument("--registry-cache-dir", type=Path, default=Path('/var/tmp/ifx-registry-cache'))
    result.add_argument('--uniprot-file', choices=['uniprot-human-reviewed.json.gz', 'uniprot-human.json.gz'],
                        default='uniprot-human.json.gz', help='Resolver input within the recorded snapshot (default: all human, including unreviewed)')
    return result


def nonnegative_int(value):
    try:
        result = int(value)
    except ValueError as exc:
        raise argparse.ArgumentTypeError('expected a nonnegative integer') from exc
    if result < 0:
        raise argparse.ArgumentTypeError('expected a nonnegative integer')
    return result


def _validate_complete(path, stage_id):
    with sqlite3.connect(f"{Path(path).resolve().as_uri()}?mode=ro", uri=True) as db:
        row = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
        if not row:
            raise ValueError("SQLite export manifest is missing")
        manifest = json.loads(row[0])
        if manifest.get("stage_id", "").split(":")[-1] != stage_id.split(":")[-1]:
            raise ValueError("SQLite manifest stage does not match the requested stage")
        if (manifest.get("scope") != "base_tables"
                or manifest.get("pending_tables") or manifest.get("pending_fields")):
            raise ValueError("SQLite export still has pending tables or fields")
        if manifest.get("status") != "complete":
            raise ValueError("SQLite export is not marked complete")
        post = manifest.get("post_processing", {})
        if post.get("pathway_tables", {}).get("status") != "computed":
            raise ValueError("Pathway similarity post-processing is incomplete")
        if not post.get("human_reaction_flags"):
            raise ValueError("Human reaction flags were not calculated")
        if db.execute("PRAGMA integrity_check").fetchone() != ("ok",):
            raise ValueError("SQLite integrity check failed")
    return manifest


def _publish_complete(temporary, output, *, overwrite):
    from src.use_cases.ramp.sqlite_writer import SQLiteWriter

    SQLiteWriter.validate_output(output, overwrite=overwrite)
    with Path(temporary).open("rb") as handle:
        os.fsync(handle.fileno())
    if overwrite:
        os.replace(temporary, output)
    else:
        # A second build may have created the output since the initial check.
        os.link(temporary, output)
        Path(temporary).unlink()


def _postprocess_complete(path, args, progress):
    from src.use_cases.ramp.postprocess_sqlite import main as postprocess_main

    progress("Calculating entity status, lookup counts, and human reaction flags")
    normal = ["--sqlite", str(path), "--with-human-flags",
              "--database", args.database,
              "--graph-credentials", str(args.graph_credentials),
              "--registry-credentials", str(args.registry_credentials),
              "--registry-cache-dir", str(args.registry_cache_dir)]
    if postprocess_main(normal) != 0:
        raise RuntimeError("RaMP entity and human-scope post-processing failed")
    progress("Calculating pathway similarity and duplicates")
    if postprocess_main(["--sqlite", str(path), "--pathway-only"]) != 0:
        raise RuntimeError("RaMP pathway similarity post-processing failed")


def _finish_complete(temporary, args, timer):
    _postprocess_complete(temporary, args, timer.mark)
    timer.mark("Validating and publishing complete SQLite")
    manifest = _validate_complete(temporary, args.stage_id)
    _publish_complete(temporary, args.output, overwrite=args.overwrite)
    return manifest


def _refresh_report_from(previous_report, sqlite_output, html_output):
    from src.use_cases.ramp.compare_source import compare, write

    previous = json.loads(Path(previous_report).read_text())
    entries = previous.get('databases')
    if not isinstance(entries, list) or len(entries) < 2:
        raise ValueError('Previous comparison must contain historical databases and a new-build entry')
    final = entries[-1]
    if not isinstance(final, dict) or final.get('label') != 'New harmonized build':
        raise ValueError('Previous comparison must end with the New harmonized build entry')
    databases = []
    for entry in entries[:-1]:
        if not isinstance(entry, dict) or not entry.get('label') or not entry.get('path'):
            raise ValueError('Previous comparison has an invalid historical database entry')
        path = Path(entry['path'])
        if not path.is_file():
            raise FileNotFoundError(f'Historical comparison SQLite is missing: {path}')
        databases.append((entry['label'], path))
    databases.append((final['label'], Path(sqlite_output)))
    return write(compare(databases), html_output)


def main(argv=None):
    args = parser().parse_args(argv)
    if args.report_output and not args.refresh_report_from:
        parser().error('--report-output requires --refresh-report-from')
    if args.refresh_report_from and (args.diagnostic or args.base_only):
        parser().error('Report refresh requires a complete SQLite build')
    from src.use_cases.ramp.sqlite_writer import SQLiteWriter
    try:
        SQLiteWriter.validate_output(args.output, overwrite=args.overwrite)
    except (OSError, ValueError) as exc:
        parser().error(str(exc))
    if not args.stage_id.strip() or not args.release_version.strip():
        parser().error("Stage ID and release version must not be blank")
    timer = BuildTimer()
    complete = not args.diagnostic and not args.base_only
    temporary = None
    try:
        import yaml
        from src.infrastructure.object_storage import load_object_storage_credentials, object_storage_from_credentials
        from src.shared.arango_adapter import ArangoAdapter
        from src.shared.db_credentials import DBCredentials
        from src.use_cases.ramp.sqlite_stage import StageReader
        from src.use_cases.ramp.sqlite_gene_identity import GeneIdentity
        from src.core.registry_integration import RegistryIntegration

        scope = ("complete SQLite" if complete else
                 "base tables; post-processing pending" if args.base_only else
                 "lookup-table diagnostic")
        print(f"Stage: {args.stage_id}\nGraph: {args.database}\nOutput: {args.output}\nScope: {scope}", flush=True)
        output = args.output
        if complete:
            output.parent.mkdir(parents=True, exist_ok=True)
            fd, name = tempfile.mkstemp(prefix=f".{output.name}.building-", suffix=".sqlite",
                                         dir=output.parent)
            os.close(fd)
            temporary = Path(name)
            temporary.unlink()
            output = temporary
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
        result = export_sqlite(reader, output, gene_identity=identity, protein_annotations=annotations, release_version=args.release_version,
                               overwrite=args.overwrite if not complete else False, source_only=args.diagnostic,
                               pathway_association_cutoff=args.pathway_association_cutoff,
                               include_lipidmaps_class_level4=args.include_lipidmaps_class_level4,
                               source_ontology_policy=args.source_ontology_policy,
                               progress=timer.mark)
        if complete:
            completed_manifest = _finish_complete(output, args, timer)
    except Exception as exc:
        timer.finish(failed=True)
        print(f"RaMP SQLite export failed: {exc}", file=sys.stderr)
        return 1
    finally:
        if temporary is not None:
            temporary.unlink(missing_ok=True)
    timer.finish()
    print(f"Created {args.output}\n" + ("Complete SQLite." if complete else
          "Base tables complete; post-processing pending." if args.base_only else
          "Lookup-table diagnostic complete; not a release database."))
    print("Pathway association cutoff: " + ("disabled" if args.pathway_association_cutoff == 0
           else f"{args.pathway_association_cutoff:,} links per reported source ID")
          + f"; excluded {result['pathway_source_ids_at_cutoff']:,} source IDs and "
            f"{result['pathway_assertions_excluded']:,} reported links")
    print("LipidMaps class level 4: " + ("included" if args.include_lipidmaps_class_level4 else "excluded")
          + f"; excluded {result['lipidmaps_class_level4_assertions_excluded']:,} reported assertions")
    print(f"HMDB Source ontology policy: {args.source_ontology_policy}")
    if not args.diagnostic:
        print("Manifest: SELECT value FROM ramp_export_metadata WHERE key = 'manifest';")
    row_counts = completed_manifest["row_counts"] if complete else result["row_counts"]
    for table, count in sorted(row_counts.items()):
        print(f"  {table}: {count:,}")
    if args.refresh_report_from:
        report_output = args.report_output or args.refresh_report_from.with_suffix('.html')
        try:
            for path in _refresh_report_from(args.refresh_report_from, args.output, report_output):
                print(f"Comparison report: {path}")
        except (OSError, ValueError, sqlite3.Error) as exc:
            print(f"SQLite complete at {args.output}; comparison report failed: {exc}", file=sys.stderr)
            return 2
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
