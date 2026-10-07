# RaMP builds

`build_ramp.py` builds the source-evidence graph. `build_sqlite.py` currently
exports a seven-table lookup and association diagnostic from a completed
harmonization stage.

## Export the lookup tables

Run from the IFX_ODIN repository root with the project environment:

```bash
.venv/bin/python -m src.use_cases.ramp.build_sqlite \
  --stage-id stage-07-b477c68b941a5659 \
  --output output_files/ramp/ramp-base.sqlite
```

The only required arguments are stage ID and output path. Optional arguments:

- `--release-version`: output label, default `unreleased`. Does not make the
  intermediate artifact release-ready.
- `--overwrite`: replace an existing output after the new export passes validation.
- `--full-base`: opt in to the earlier full base-table export. By default, the
  diagnostic contains `analyte`, `source`, `analytesynonym`, `pathway`,
  `ontology`, `analytehaspathway`, and `analytehasontology`;
  it is not an R-package database.
- `--database`: defaults to `metabolite_harmonization`.
- `--graph-credentials`: Arango credential YAML file, defaults to
  `src/use_cases/secrets/ifxdev_arangodb.yaml`.
- `--curation-credentials`: object-storage credential YAML file, defaults to
  `src/use_cases/secrets/aws_ifx_registry.yaml`.
- `--registry-credentials`: Registry credential YAML, with the same default.
- `--registry-cache-dir`: local Registry cache, default `/var/tmp/ifx-registry-cache`.
- `--uniprot-file`: `uniprot-human.json.gz` (default, reviewed and unreviewed
  human entries, matching the old RaMP scope) or `uniprot-human-reviewed.json.gz`.
  Both are read from the exact UniProt human snapshot
  recorded in the graph build, currently `uniprot:human:2026_03`.

The command prints elapsed wall time for each completed phase and the total
build (`HH:MM:SS.s`). If a phase fails, it prints the time spent in that phase
and the elapsed time before failure.

The command reads the existing graph and immutable curation batches. It does not
rebuild the graph, perform another harmonization pass, or publish the stage to
the Registry. Avoid graph builds and stage changes while exporting. The command
requires a graph schema with `MetaboliteIdentifier.hmdb_status`; older graphs
must be rebuilt and harmonized before export. The command
rejects source revisions inconsistent with the selected stage and verifies the
stage and all read collections again before publishing the SQLite file.
The stage predates this exporter and records historical revisions for only the
collections used in harmonization. Direct manual edits to other collections
before export cannot be detected retrospectively; export from the unchanged
source build. The build fingerprint detects a different ETL run or input versions.

Existing files are refused by default. To rebuild at the same path:

```bash
.venv/bin/python -m src.use_cases.ramp.build_sqlite \
  --stage-id stage-07-b477c68b941a5659 \
  --output output_files/ramp/ramp-base.sqlite --overwrite
```

Rows are written into a temporary database; only a validated artifact is
published. With `--overwrite`, publication atomically replaces the output and
failure preserves the prior file. Allow space for both the old database and the
new build, including raw rows, indexes, and compaction. Use an ordinary output
file, not a symlink, and close any process writing to it; SQLite journal/WAL
sidecars prevent replacement.

## Current scope

The default export focuses on lookup data. `source` contains all IDs each input used to report
retained RaMP data, including associations and chemistry. A source row maps a
reported ID to a RaMP analyte and attributes it to the reporting input;
`dataSource` and the ID namespace `IDtype` can differ. Direct chemistry
assertions register their reported source ID. A ChEBI bridge adds a ChEBI row
only when its ChemicalEntity provides retained fallback chemistry. UniProt
lookup aliases remain attributed to UniProt. The diagnostic is not release-ready.
`analytesynonym` keeps each source-owned metabolite name or synonym with its
attributed source. HMDB protein names, gene symbols, and synonyms stay under
HMDB; matched UniProt protein names and gene symbols/aliases from the same pinned
resolver input appear under `uniprot`. Names from merged protein nodes are not
copied to Rhea or Reactome simply because those sources supplied an association.
WikiPathways HGNC-symbol links from the pinned RDF are stored as source-attributed
names on the corresponding gene or protein identifier and exported as `wiki`
synonyms. The current graph was backfilled in place; the audit and backup are in
`output_files/ramp/wikipathways-symbol-repair/`.
Rows are deduplicated by all four stored columns, retaining the same spelling
under different sources or RaMP IDs.

The diagnostic also carries pathway and HMDB ontology associations with their
target lookup tables. `analytehaspathway.pathwaySource` is the stored source
attribution, including `kegg` for HMDB-supplied KEGG pathways. The ontology
association table has no source column; its HMDB attribution follows from its
HMDB ontology terms. `ontology.metCount` remains pending.

The earlier full base-table export remains available with `--full-base`. It
retains the legacy table/column/index contract with these explicit differences:

- `reaction_protein2met` is omitted, as agreed.
- `version_info.data_source_snapshot_ids` lists the exact Registry input IDs.
- `ramp_export_metadata` stores one JSON manifest under key `manifest`.
- `entity_status_info`, `pathway_similarity`, and `pathway_duplicates` are empty
  until post-processing is implemented.
- Pending counts, derived flags, and intersections are `NULL` where permitted,
  or `-1` where the legacy schema requires an integer. The manifest enumerates
  these fields; they must not be interpreted as completed calculations.
- Gene/protein identifiers present in the RaMP inputs connect all their returned
  canonical UniProt matches into groups. Overlapping groups merge transitively
  and receive one sequential `RAMP_G_*` ID. Unused UniProt aliases cannot connect
  groups; unmatched identifiers retain separate typed identities.
  Original IDs, source attributions, and relationship
  evidence are preserved; canonical UniProt IDs are also queryable in `source`.
  When merged protein records have different types, `catalyzed.proteinType`
  retains sorted distinct values separated by `; ` (missing values are `Unknown`).
  This input-driven collapsing policy preserves the legacy RaMP grouping approach.
  The source graph remains unchanged.

The scientific grouping decision is isolated in `sqlite_gene_grouping.py` as
`collapse_input_matches` (policy `input_identifier_connected_components`, version
2). `sqlite_gene_identity.py` handles the pinned resolver and collects its matches;
`sqlite_projection.py` assigns RaMP IDs and writes rows. The grouping policy can
therefore be changed independently of resolver construction or database writing.

Resolution uses the existing `UniProtResolver`, including its secondary-accession
and isoform-to-canonical mapping, explicit NCBI Gene and Ensembl cross-references,
and source identifiers explicitly typed as symbols. It does not use arbitrary
display names as identity evidence. The export manifest records the resolver,
policy name/version, input snapshot, filename, checksum, and resolved/ambiguous/unmatched
input counts, plus multi-accession group counts and the largest group size.
An ambiguous input now connects its returned matches; that count does not mean
it remains standalone. A missing input or mismatch with the graph's recorded checksum fails
before publication. The Registry pin already appears in `version_info`.

Fresh deterministic RaMP IDs are assigned for this input. They are not stable
across changed inputs or historical RaMP releases. Singleton metabolites are
included. Suppressed and cleanup-excluded metabolites do not reappear through
their pathway, reaction, classification, or protein-association edges; skipped
edge counts are included in the manifest.

Source identifier spelling is converted to the legacy namespace conventions
(for example, `CHEBI:15377` becomes `chebi:15377`); original source IDs are retained.
Unknown namespaces are preserved. Metabolite `analyte.common_name` is the most
frequent nonblank `source.commonName` within its RaMP ID, counted case-insensitively
across the finalized source rows as in the legacy build. Tied names prefer the
shortest spelling, then lexical order;
groups without named rows fall back to their most frequent source ID, then their
RaMP ID. Individual source rows retain
source-specific names. HMDB, RefMet, and WikiPathways names on reported
cross-reference IDs come from the same provider record; Rhea names come from
its reaction evidence. For LipidMaps `source.commonName`, a LipidMaps
`ABBREVIATION` takes precedence, followed by a `SYNONYMS` value, then `NAME`,
`COMMON_NAME`, and `SYSTEMATIC_NAME` from the same identifier. Names are never
copied to another LipidMaps-reported ID in a RaMP group.
Names are not borrowed from another provider to fill missing source rows.
For gene/protein `analyte.common_name`, the exporter prefers a source-reported
`gene_name` symbol, then the pinned UniProt protein name, then another available
name. An identifier with its namespace prefix is used only when the group has
no readable name; this fallback remains searchable in `source`.
The current PubChem pin includes compound Titles, which are used as
PubChem-attributed source names when present. CIDs without a Title retain an
empty PubChem source name.
Chemistry retains source-reported keys and effective
curated values; derived InChIKeys do not replace reported keys.
If an active `MetaboliteIdentifier` has no `chem_props`, a same-CURIE link to a
`ChemicalEntity` supplies a ChEBI `chem_props` row from its reported structure,
formula, mass, InChI, and InChIKey. The corresponding source lookup row is
attributed to ChEBI only when fallback chemistry is retained. The
export manifest counts bridged identifiers, fallback chemistry rows, and groups
gaining their first chemistry row.

After gene/protein groups are assigned, the exporter adds lookup identifiers
from their matched records in the same verified UniProt input. These include
primary/secondary accessions, gene names and synonyms, Entrez Gene IDs, HGNC IDs
(`hgnc:<number>`), and Ensembl transcript/gene/protein IDs with version suffixes
removed. Added rows use `dataSource='uniprot'`; source-reported rows keep their
provider attribution. Shared aliases can point to several existing RAMP_Gs and
never cause additional group merging. Protein display names and arbitrary
cross-references are not added as identifier aliases. The annotation manifest
records the selection policy and accession/alias counts across the input file.

Class and ontology ancestors are expanded explicitly. Ontology exports apply
the legacy denylist in `sqlite_ontology_policy.yaml`. Identical output rows are
deduplicated; conflicting rows sharing a legacy primary key cause failure.

## Inspect and validate the output

The command prints per-table counts and scope. The full base export also stores
stage identity, source dataset metadata, curation snapshots, and graph revisions
inside its SQLite metadata table. Both modes can be checked with:

```sql
PRAGMA integrity_check;
PRAGMA foreign_key_check;
```

`ramp_export_metadata` and `version_info` are available only with `--full-base`.

Tests (no network or graph required):

```bash
.venv/bin/python -m pytest tests/test_ramp_sqlite_export.py \
  tests/test_record_property_curations.py tests/test_curations.py \
  tests/test_metabolite_harmonization_workbench_rules.py --no-cov -q
```

The full export and R-package validation are run by the maintainer. Post-processing,
curation audit columns, and optional consumer changes are separate subsequent
phases. See `designs/ramp_sqlite_export_design.md` for discovery and decisions.

## Compare SQLite releases

Generate a standalone HTML report and matching JSON from two or more local
SQLite files. Every comparison table uses the preceding displayed database for
absolute and percentage changes. Materialize historical `ramp:sqlite_database:<version>` artifacts from
the Registry and decompress them before running:

```bash
.venv/bin/python -m src.use_cases.ramp.compare_source \
  --database 'Released 2.3.2=/path/to/RaMP_SQLite_v2.3.2.sqlite' \
  --database 'Released 2.4.0=/path/to/RaMP_SQLite_v2.4.0.sqlite' \
  --database 'Released 3.0.7=/path/to/RaMP_SQLite_v3.0.7.sqlite' \
  --database 'Unreleased 3.0.12=/path/to/RaMP_SQLite_v3.0.12.sqlite' \
  --database 'New harmonized build=output_files/ramp/ramp-base.sqlite' \
  --output output_files/ramp/comparison.html
```

The command opens inputs read-only and needs no graph connection. It first
compares distinct metabolite and gene/protein RaMP IDs in `analyte`, then shows
the `analyte` cells for D-glucose and EGFR matched across builds by their
source IDs. Two subsequent tables count distinct RaMP IDs represented by each
`source.dataSource`, separately for metabolites and genes/proteins. Expandable
examples show a matched metabolite and gene/protein `source` record per input
where available, including all eight source columns. Synonym sections then compare row counts and distinct
RaMP IDs by attributed source and analyte type, with one real stored example
per populated database/source/type. Example rows are paired by analyte type
within each input; an unsupported type is labeled explicitly. A missing synonym table is shown as
"Table absent" rather than zero. Separate tabs compare HMDB ontology links by
ontology type and pathway links by stored source and analyte type. They show
association rows, distinct RaMP IDs and distinct terms/pathways, plus examples
matched by source identities across builds. Missing association tables are
shown as "Table absent". The lookup ambiguity and integrity checks follow.
Deltas compare each build with the preceding one. RaMP IDs are counted within
a build, never matched across releases. Count increases/decreases are shaded
green/red even when the previous count was zero; example field-coverage shading
uses percentage-point change in populated rows, with its signed change shown
as text.

```bash
.venv/bin/python -m pytest tests/test_ramp_source_comparison.py --no-cov -q
```
