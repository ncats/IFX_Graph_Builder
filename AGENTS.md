## Specialist Agent Workflow

- Project-scoped specialist definitions live in `.codex/agents/`; the workflow
  and role boundaries are summarized in `designs/ai_specialist_agents.md`.
- Before implementing changes to a user interface, interaction, navigation,
  information hierarchy, or user-visible workflow, delegate a read-only
  consultation to `ux_designer`. Bring material product choices and its
  recommended default back to the user before implementation.
- Before implementation, delegate a read-only consultation to
  `software_architect` when a change materially alters an ownership boundary,
  public API contract, persistence model, dependency direction, or introduces
  a new architectural pattern. Do not require this consultation for routine
  extensions that follow an established design.
- When both consultations apply, they may run in parallel. Wait for both and
  reconcile their recommendations before editing.
- After implementing and testing a feature or substantial refactor, delegate a
  read-only review to `code_reviewer` before the final handoff.
  Address Blocking and Important findings. If that causes substantive changes,
  ask the same reviewer for one focused re-check.
- Do not use specialist passes for documentation-only changes, trivial
  mechanical edits, or when the user explicitly asks not to delegate. The user
  may explicitly request any specialist at any time.
- The primary agent remains responsible for decisions, edits, validation, and
  the final answer. Specialists advise and return concise evidence; they do not
  make parallel code changes.
- Treat user feedback during development as possible evidence that a specialist
  role should improve. When a comment reflects a reusable preference, decision
  lens, or workflow that belongs to one of the roles, seriously evaluate it and
  update that agent's instructions when appropriate. Keep feedback that is
  specific to one feature in the feature's design or implementation instead of
  overfitting a general-purpose agent. Tell the user when feedback changes an
  agent's lasting instructions.

## Current Local TDL Handoff

- For Jess's local `impatient_target_graph` reconstruction, read
  `src/use_cases/pharos/impatient_target_graph_local_refresh.md` before editing
  Pharos graph YAML, target resolvers, TDL override handling, or local registry
  cache code.
- The local reconstruction is scratch/test work until Jess decides whether to
  ask Keith to review it. Do not promote local YAML or registry-cache changes to
  production Pharos configs without explicit approval.

## Ingest Workflow Preferences

- For any ingest task, start by reading `playbooks/ingest_playbook.md` for new sources or `playbooks/data_update_playbook.md` for refreshes.
- For new ingest work, start with source discovery: identify the files, inspect real payload shape, and validate assumptions before implementing adapter/model changes.
- For new ingest work, pause after discovery and provide a short implementation plan; get user confirmation before making code changes.
- For source refreshes, profile the changed payloads before changing adapters, resolvers, or downstream conversion logic.
- For Snakemake runs and ETL executions, let the user run those steps unless they explicitly ask otherwise.

## Ingest Standards

- Start new ingest development in `src/use_cases/working.yaml`; only promote to `src/use_cases/pharos/target_graph.yaml` after validation.
- Keep source-specific mapping/coverage decisions documented in a design doc under `designs/`.
- Prefer deriving and persisting datasource metadata (`version`, `version_date`, `download_date`) during download/prep, then have adapters consume it.
- Keep first-pass ingest scope intentionally minimal, then expand in follow-up iterations.
- Use stable IDs and consistent prefixes across nodes/edges.
- Validate assumptions against real payloads (field presence, cardinality, identifier shape) before finalizing model changes.
- Avoid speculative parsing when source text is ambiguous; preserve source text when parsing would be lossy.
- In adapters, dedupe repeated entities in-memory using deterministic IDs.
- Emit ingest output in type-grouped batches (for example: primary nodes, related nodes, then edges).
- When validating a new source, compare against any prior TCRD or legacy ETL only after confirming what the current raw payload actually contains.

## Adapter vs Resolver Boundary

- Use the plain-language test: an ID resolver answers **"who is this?"**; an adapter answers **"what does this source say about it?"**
- Adapters parse source payloads and emit graph structure using the identifiers the source actually provides.
- Adapters may normalize identifier syntax only within the same identifier family when the transformation is lossless and source-defined, for example trimming whitespace or converting `DOID_1234` to `DOID:1234`.
- Adapters must not perform cross-identifier reconciliation, for example symbol-to-UniProt, UniProt-isoform-to-canonical-accession, Entrez-to-Ensembl, disease-xref-to-MONDO, or source-specific IDs to IFX canonical IDs.
- Cross-identifier reconciliation belongs in `IdResolver` implementations and resolver configuration. If a source ID needs to map to an existing `Protein`, `Gene`, `Disease`, `Tissue`, or other canonical node, emit the source ID and add or configure the resolver path.
- Adapters should not invent `xref` values to make downstream joins work. Canonical/equivalent `xref` content should come from resolver matches unless the source explicitly provides those exact equivalent IDs.
- Filtering is different from resolving. Adapters should still apply source-scope filters that decide whether a row belongs in the ingest at all, such as human-only filters, evidence thresholds, or source-declared subsets.
- If implementation feels like it needs a lookup table from one identifier family to another, stop and decide whether that lookup belongs in a resolver, a converter, or a documented post-processing step rather than inside the adapter.

## Harmonizer vs Graph Builder Boundary

- A harmonized or Registry-derived artifact is a versioned input, not permission to bypass the adapter/resolver boundary. Classify the artifact by what it represents before wiring it into a graph build.
- Prefer using a harmonizer output in two explicit ways when it contains both identity mappings and harmonizer-owned annotations:
  1. Build an `IdResolver` from its source-ID-to-preferred-ID mappings.
  2. Use a separate adapter to ingest only the harmonizer's own scores, methods, confidence, review decisions, or other derived assertions.
- Keep primary-source payloads independently ingestible through their source adapters. Graph configuration must remain able to include, exclude, version, filter, and model sources independently when their semantics differ.
- Do not use a combined application, presentation, QA, or explorer graph as a drop-in replacement for primary-source adapters unless an explicit design decision transfers source selection, payload ownership, scientific policy, and provenance responsibility to that upstream product.
- Treat these as warning signs that a proposed handoff has combined **"who"** and **"what"**:
  - the artifact already uses preferred IFX IDs for records originating from several providers;
  - provider rows or edges have been collapsed before the graph builder sees them;
  - the adapter reports the harmonizer as its datasource even though its fields or evidence originated from multiple primary sources;
  - adding a primary-source field requires changing and rebuilding the harmonizer;
  - one adapter implicitly enables or disables several providers together.
- Source names retained in a combined string or nested evidence blob do not replace field-level provenance or independent source selection. Information omitted upstream cannot be reconstructed by the graph builder.
- Treat evidence-bearing relationships more strictly than preferred scalar annotations. A harmonizer may publish a preferred name, symbol, location, or similar consensus value when it also publishes the selection method, score, and evidence trail. It must not invisibly filter, sample, rank, cap, or collapse primary-source evidence edges intended for graph building. Those decisions belong in visible graph configuration or source adapters unless an explicitly approved derived-evidence contract preserves provider identity, raw endpoints, row-level evidence, and lossless source selection.
- Separate ownership of Target Graph, Pharos, or another graph-building use case does not waive this boundary. Each use case may choose different sources and policies while still keeping identity resolution separate from source payload modeling.
- Before implementing an adapter over a pre-resolved or multi-provider artifact, read `designs/harmonizer_graph_builder_boundary.md`, document the handoff contract, and pause for explicit architectural approval.

## Profiling And Planning Tools

- Examine the old Pharos loader code when relevant, usually under `https://github.com/unmtransinfo/TCRD/tree/master/loaders`.
- Examine current Pharos MySQL contents when relevant using `src/use_cases/secrets/pharos_credentials.yaml` against the `pharos319` schema.
- Examine new Pharos MySQL contents when relevant using `src/use_cases/secrets/pharos_credentials.yaml` against the `pharos400` schema.
- Examine the graph staging database on `ifxdev` when relevant, especially for ingests that may already land in the graph built by `build_pharos.py`.
- Examine the raw input files directly and verify whether their data lands correctly in the graph database and in MySQL outputs when a MySQL path exists.

## Fail Loudly on Broken Requirements

- Diagnose the root cause when required inputs, assumptions, or validation fail.
  Do not invent a fallback, substitute another source or older data, skip a
  required step, return an empty/default result, or weaken a check just to make
  the run succeed.
- Stop the affected operation with an actionable error that identifies what
  failed and preserves the original cause. Do not report incomplete work or
  degraded data as a successful build, capture, or release.
- Retries and optional behavior are appropriate only when the existing contract
  allows them and they preserve correctness. Make any degraded state visible;
  a warning is not a substitute for failing a required operation.
- Fix the underlying issue. If continuing requires changing the intended
  behavior or data contract, explain the tradeoff and get an explicit decision
  rather than silently introducing a workaround. Keep the fix minimal.

## Lessons Learned

- Keep adapters focused on source parsing and structural graph emission; move cross-ontology ID normalization to resolvers.
- When canonical/equivalent `xref` content should come from identifier reconciliation, populate it via `IdResolver` matches rather than inventing it in the adapter. The ETL framework applies resolvers in `InputAdapter.get_resolved_and_provenanced_list()`, and `IdResolver.get_merged_map()` writes `IdMatch.equivalent_ids` onto resolved nodes as `xref`.
- Raw source adapters should not populate `sources` or `provenance`; let the framework stamp canonical datasource/version metadata during ETL.
- For ontology xrefs, maintain an explicit allowlist and perform case-insensitive prefix checks.
- When adding new datasource version handling, use named parameters for `DatasourceVersionInfo` to avoid argument-order regressions.
- When an edge can be emitted by multiple sources and later merged, keep source-specific payload in a `details` list instead of top-level edge fields.
- Put disease descriptions on `Disease` nodes (for example `uniprot_description`), not on disease association edge details.
- Prefer modern ontology-backed evidence codes from the source (for example ECO) over recreating legacy source-specific evidence encodings from older Pharos/TCRD tables.
- If a build shows inexplicable dataclass/constructor signature errors for shared graph models, suspect import-order-sensitive model initialization. Prefer adding the affected modules to the centralized model preload in `src/use_cases/build_from_yaml.py` rather than scattering one-off imports through adapters or ETL code.

## Workflow References

- Use `playbooks/ingest_playbook.md` for new-source ingest work.
- Use `playbooks/data_update_playbook.md` for source refreshes and payload-drift investigations.
- Use `playbooks/pharos_table_update_playbook.md` when updating the Pipeline Status Table in `src/use_cases/pharos/TCRD_TODO.md`.

## Documentation Model

- `AGENTS.md` provides repo-level instructions and routing.
- `playbooks/` contains repeatable repo workflows for specific classes of work.
- `designs/` contains source-specific investigation notes, design decisions, and validation outcomes.
- For ingest work, `AGENTS.md` should point to the relevant playbook, and the playbook should point to the expected design-doc artifact.

## Project Overview

IFX_ODIN (Ontology and Data Integration Network) is a modular ETL platform developed at NCATS for biomedical data integration. It powers translational science platforms including Pharos, RaMP-DB, GSRS, and CURE ID.

**Primary contacts:** Jessica Maine, Keith Kelleher

## Tech Stack

- **Language:** Python 3.11+
- **Workflow:** Snakemake
- **ORM:** SQLAlchemy
- **Data models:** dataclasses, LinkML
- **Databases:** MySQL, ArangoDB, Neo4j/Memgraph
- **APIs:** FastAPI, Strawberry GraphQL

## Project Structure

```
src/
├── constants.py           # DataSourceName, Prefix enums (200+)
├── core/                  # ETL engine, config loading
├── input_adapters/        # Data source adapters (40+ sources)
├── output_adapters/       # Database writers (MySQL, ArangoDB, Neo4j, JSONL)
├── interfaces/            # Abstract base classes
├── models/                # Dataclass models (Node, Relationship, etc.)
├── id_resolvers/          # Entity normalization
├── shared/                # Utilities, SQLAlchemy tables
└── use_cases/             # Build scripts, YAML configs, secrets
workflows/                 # Snakemake files for data download
tests/                     # pytest tests
*-container/               # Docker compose for databases
```

## Build & Run Commands

```bash
# Environment setup
python -m venv .venv
source .venv/bin/activate
export PYTHONPATH=/path/to/IFX_ODIN:$PYTHONPATH
pip install -r requirements.txt -c constraints.txt

# Run tests
pytest

# Download data files
cd workflows && snakemake -j 4

# Example ETL entrypoint a user may run for validation
python -c "
from src.use_cases.build_from_yaml import BuildGraphFromYaml
builder = BuildGraphFromYaml('src/use_cases/pounce_v2.yaml')
builder.truncate_datastore()
builder.do_etl(do_post_processing=True)
"

# Start containers
cd mysql-container && docker compose up -d
cd memgraph-container && docker compose up -d
```

## Architecture

**Adapter Pattern + ETL Pipeline:**

1. **InputAdapters** - Extract from data sources via `get_all()` generator
2. **OutputAdapters** - Write to databases via `store()`
3. **ETL Engine** - Orchestrates pipeline with resolver_map
4. **Id Resolvers** - Entity normalization with configurable match behavior

```python
# Creating an input adapter
class CustomAdapter(InputAdapter):
    def get_all(self) -> Generator[List[Union[Node, Relationship]], None, None]:
        yield batch  # Batches of 25,000 by default

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.XXX
```

## Key Patterns

- **Generators** for memory-efficient batch processing
- **YAML configs** drive adapter selection and parameters
- **Provenance tracking** on all entities
- **@facets decorator** marks searchable fields on models
- **RecordMerger** handles field conflicts during data integration

## Database Access

Credentials in `src/use_cases/secrets/local_credentials.yaml`:
- MySQL: localhost:3306
- Memgraph Lab: http://localhost:3000
- Memgraph query: localhost:7687

## Workflow Conventions

- **Always load the latest harmonizer release.** For every versioned entity
  explorer (disease, target, drug, variant, and pathway), discover all valid
  bundled releases and select the highest semantic version when no explicit
  graph directory is configured. Explicit CLI/deployment configuration retains
  precedence. A
  legacy `current/` directory must not override a newer versioned release, and
  a newer `current/` manifest must not be ignored because older `v*` folders
  exist.
- **Complete the release handoff.** After a successful harmonizer release,
  compare its app-ready package with the IFX_ODIN bundle and update the bundle
  or deployment configuration in the same body of work. Verify manifest
  version, required files, row counts, release diffs, and source catalog.
- **Never relabel stale data.** The version shown in the UI must match the data
  actually loaded. If a version-matched large resolver/index must be mounted
  outside Git, fail or clearly report the missing deployment artifact rather
  than silently presenting bounded preview coverage as complete.

- **Do not assume Git staging is the user's unit of work.** The user organizes changes in IntelliJ changelists. Prefer leaving commit grouping to the user unless they explicitly ask for staging or committing help, and be careful not to mix unrelated changelists.
- **Do not keep one-time scripts in the repository by default.** Run migration,
  audit, repair, and release-compilation helpers from a temporary location when
  they are genuinely single-use. Preserve the resulting reviewed artifacts,
  reports, or durable product code where appropriate, but add a script to the
  repository only when it defines a repeatable supported workflow.
- **Keep changes in collaborator-owned repositories minimal.** When working in
  IFX_Harmonizers or another repository primarily owned by a collaborator,
  make the smallest boundary-preserving change needed for the approved task.
  Preserve recognizable entry points, class names, configuration shape, file
  paths, and workflow structure unless a change is necessary. Keep unrelated
  refactoring or cleanup as a separate proposal.
- **For new ingest sources**, read `playbooks/ingest_playbook.md` before doing discovery or code changes.
- **For source refreshes**, read `playbooks/data_update_playbook.md` before investigating payload drift or changing code.
- **For Pipeline Status Table updates**, read `playbooks/pharos_table_update_playbook.md` before changing `src/use_cases/pharos/TCRD_TODO.md`.
- **Do not skip discovery.** Inspect the real files first, record mapping decisions in `designs/`, then propose the implementation plan before editing code.
- **Do not run Snakemake or ETL on the user's behalf by default.** Prepare the changes and tell the user exactly what to run for validation.

## Key Files

| File | Purpose |
|------|---------|
| `src/core/etl.py` | Main ETL orchestrator |
| `src/core/config.py` | YAML config parsing |
| `src/interfaces/input_adapter.py` | Base adapter class |
| `src/interfaces/output_adapter.py` | Output writer base |
| `src/models/node.py` | Core Node/Relationship classes |
| `src/use_cases/build_from_yaml.py` | Generic YAML-driven builder |

## Current Work

- POUNCE v2 database with UniProt protein annotations
- Schema migration: v1 models in `models/pounce_v1/`, v2 in `models/pounce/`
- SQLAlchemy table refactoring (Sample -> Biosample)

### POUNCE Input Adapter (`src/input_adapters/pounce_sheets/`)

Parses three Excel workbooks per project (all configured via `pounce_v2.yaml`):

| Workbook | Sheets | Purpose |
|----------|--------|---------|
| **Project** | ProjectMeta, BioSampleMap, BioSampleMeta | Project metadata, biosamples, biospecimens, exposures |
| **Experiment** | ExperimentMeta, RunBioSampleMap/Meta, GeneMap/Meta or MetabMap/Meta, RawData/Meta or PeakData/Meta | Experiment metadata, run-level samples, analytes, raw data matrix |
| **StatsResults** | StatsResultsMeta, StatsReadyData, EffectSize_Map, EffectSize | Statistical analysis results (normalized data, effect sizes) |

Sheet constants are in `constants.py` (`ProjectWorkbook`, `ExperimentWorkbook`, `StatsResultsWorkbook`). Map sheets define column name mappings from NCATS standard names to submitter names. Meta sheets contain the actual data rows.

**Data model graph:**
```
Project -> Experiment -> Dataset -> RunBiosample -> Biosample <- Project
                                -> Gene/Metabolite
Project -> Biosample -> Biospecimen
                     -> Exposure
```

- `Dataset` nodes carry a `_data_frame` field written to parquet by the ArangoOutputAdapter
- Dataset edges to RunBiosample and Gene/Metabolite make the data matrix relationships explicit in the graph
- `_parse_data_matrix()` is generic: handles RawData, PeakData, and StatsReadyData via `parser`, `analyte_id_col`, and `default_data_type` parameters
- Empty data sheets (e.g. PeakData for metabo) are skipped gracefully

**Supported templates:** Transcriptomics (genes, RawData) and Metabolomics (metabolites, PeakData). EffectSize parsing is not yet implemented.
