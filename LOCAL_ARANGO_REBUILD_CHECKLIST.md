# Local Arango Rebuild Checklist After the IFX Registry Cutover

Created 2026-09-15 for the Registry cutover in `a981dc3`. This is a local
working checklist, not a claim that any shared database is ready to replace.
The list is limited to Arango databases whose build configuration changed as
part of the cutover.

Status checked against ifxdev on 2026-09-17. A database is checked off only
when its completed build timestamp is later than the Registry cutover commit
at 2026-09-15 10:34 EDT. None currently meet that criterion.

## Primary graphs

- [ ] `target_graph`
  - Last completed ETL metadata: 2026-07-10 14:35 EDT (pre-cutover).
  - Build: `src/use_cases/pharos/target_graph.yaml`
  - Post adapters and graph views are in the build YAML.
  - Registry changes: source-backed target, HCOP, and tissue resolvers;
    external and derived dataset references.
- [ ] `pharos`
  - Last completed ETL metadata: 2026-07-08 10:13 EDT (pre-cutover).
  - Build: `src/use_cases/pharos/pharos.yaml`
  - Post adapters and graph views are in the build YAML.
  - Registry changes: source-backed TCRD, HCOP, and tissue resolvers;
    external and derived dataset references.
  - Follow-up outside this checklist: rebuild/validate `pharos400` from this
    graph only after the Arango graph passes validation.
- [ ] `metabolite_harmonization`
  - Last completed ETL metadata: 2026-08-26 14:10 EDT (pre-cutover).
  - Build: `src/use_cases/ramp/ramp.yaml`
  - Do not also run `src/use_cases/working.yaml` against the same database;
    choose the intended configuration.
  - Registry change: explicit `derived_snapshot` input for PubChem molecular
    information.
- [ ] `cure_rasopathies`
  - Last completed ETL metadata: 2026-09-15 10:13 EDT (pre-cutover by
    approximately 21 minutes).
  - Build: `src/use_cases/cure/cure_rasopathies.yaml`
  - Registry change: the CURE label resolver now consumes the pinned curated
    concepts source dataset directly.
- [ ] `frdb`
  - Status on ifxdev: database not found.
  - Build: `src/use_cases/frdb.yaml`
  - Registry change: removed the legacy Translator resolver snapshot.
- [ ] `ccle`
  - Status on ifxdev: database not found.
  - Build: `src/use_cases/pounce/ccle.yaml`
  - Registry change: the gene resolver now consumes its pinned source dataset
    directly.
  - Note: `src/use_cases/pounce/build_ccle.py` currently points at the wrong
    YAML path, so do not use that entry point without correcting it.

## Temporary and comparison graphs

Rebuild these only if they are still being used.

- [ ] `impatient_target_graph`
  - Last completed ETL metadata: 2026-06-23 01:12 EDT (pre-cutover).
  - Shared build: `src/use_cases/pharos/impatient_target_graph.yaml`
  - Shared post-process: `src/use_cases/pharos/impatient_target_graph_aql_post.yaml`
  - Registry changes: source-backed target resolvers plus explicit derived and
    external dataset kinds.
  - [ ] Local reconstruction configuration migrated before rebuilding.
    `src/use_cases/pharos/impatient_target_graph_local.yaml` and its local
    post-process YAML still use the removed `resolver_snapshot` contract. Read
    `src/use_cases/pharos/impatient_target_graph_local_refresh.md` before
    changing them; keep that work local until Jess approves promotion.
- [ ] `impatient_pharos`
  - Last completed ETL metadata: 2026-06-18 05:46 EDT (pre-cutover).
  - Build: `src/use_cases/pharos/impatient_pharos.yaml`
  - Post adapters and graph views are in the build YAML.
  - Registry changes: source-backed TCRD resolver plus explicit external
    dataset kinds.
- [ ] `pharos_current_tdls_reconstruction`
  - Last completed ETL metadata: 2026-06-30 21:35 EDT (pre-cutover).
  - Build: `src/use_cases/pharos/pharos_current_tdls_reconstruction.yaml`
  - Post-process:
    `src/use_cases/pharos/pharos_current_tdls_reconstruction_aql_post.yaml`
  - Registry changes: source-backed TCRD resolver plus explicit external
    dataset kinds.

## Per-database completion checks

- [ ] Run a fresh build, not `--resume`, so old Registry metadata and resolver
  fingerprints cannot survive from the pre-cutover graph.
- [ ] Confirm configured post adapters completed; active Pharos and Target Graph
  build scripts run them from the primary YAML. Older paired builds still run
  their companion post-processing YAML.
- [ ] Confirm `metadata_store/etl_metadata` records Registry datasets with both
  `kind` and `snapshot_id`.
- [ ] Confirm resolver fingerprints contain `dataset_inputs` and no legacy
  `resolver_snapshot` dependency.
- [ ] Review adapter failures/skips and representative collection and edge
  counts before marking the database complete.
- [ ] Validate any downstream MySQL or RDF export only after its source Arango
  database is complete.

## Deliberately not included

- `pounce`, `cure_pasc`, `chebi`, and `gramp`: their graph YAML did not change
  in the Registry cutover commit.
- `pharos400` and `pharos400_working`: MySQL databases, not Arango databases.
- Snakemake, ETL, and conversion execution: left for the user to run.
