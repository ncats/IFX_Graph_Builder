# RaMP SQLite export: discovery and phased implementation

Status: phase-one exporter implemented, 2026-10-05; full graph export remains maintainer-run.

The one-time audit and graph-repair paths mentioned below are members of
`output_files/ramp/one-time-audits-20261008.tar.gz`; their original loose files
were removed during output cleanup.

## Accepted user decisions

- Omit `reaction_protein2met` going forward.
- Assign fresh RaMP IDs. Cross-release continuity was never promised; historical
  ID mapping is deferred unless the team requests it. Make assignment repeatable
  for identical export inputs without implying stability across changed inputs.
- Initial input: `stage-07-9cfef0c333530e83` in `metabolite_harmonization`.
  Read it directly; do not require Registry registration/publication of the stage.
- Include pinned input version IDs in output version metadata. Keep source
  snapshot IDs, curation batch identities, graph stage identity, and output
  release identity distinct.
- The user identifies `https://github.com/ncats/ifx-frontend-library` as the
  website repository and reports that it accesses data through the R package.
  No frontend change is planned for phase one.

## Scope and sequence

The output is a RaMP application database projected from the selected metabolite
harmonization stage and its associated source evidence. Use a dedicated exporter
with explicit legacy table mappings. Do not derive the application schema from
arbitrary graph collection metadata or repeat identity harmonization.

The user requested this order:

1. Basic table structure and graph-to-row mappings.
2. Post-processing.
3. Audit fields documenting source-data curations.
4. Optional improvements, including version-compatible RaMP-DB data-access changes.

## Deferred compatibility improvements

- `source.pathwayCount` is an analyte-level count of distinct non-HMDB
  pathways, repeated on every source ID row for that RaMP ID. The legacy
  schema and R package expect it there, so the current export preserves it.
  In a later version-compatible schema/access-layer update, make `analyte`
  the authoritative home for this value and retire the repeated copy only
  after existing R-package queries continue to work.
- For future SQLite builds, calculate the RaMP-specific reaction human-scope
  flags during projection, while the pinned ChEBI and Rhea evidence is already
  available. The existing SQLite requires one graph-backed post-processing pass
  to recover the historical directional-protein flag; later count refreshes
  can use SQLite alone. Do not turn the RaMP ChEBI-root policy into a generic
  `ChemicalEntity` assertion.

Preserve existing application behavior initially. Additional chemistry providers
are in scope; their rows should not be excluded just because older releases lack
them. Audit-field design is deferred, but using the effective curated values is
part of correct base-data export.

## Verified comparison artifacts

Both files were materialized through the Registry client into a temporary cache,
decompressed, and inspected using SQLite read-only connections.

| Role | Registry snapshot | Compressed artifact SHA-256 |
| --- | --- | --- |
| Released baseline | `ramp:sqlite_database:3.0.7` | `bc1f7b04fb66f1017a542a5d7a4116c32338f26120beb22121c8cfae0ceb1e1e` |
| Unreleased comparison | `ramp:sqlite_database:3.0.12-ramp4.0` | `969f619dfb2fd9a7b827f1bf0efa5f2cfb46fb9c65ee43806fdca25abbd63adf` |

The second manifest explicitly records `release_status: unreleased`, branch
`ramp4.0`, and internal database version `3.0.12`. Keep these identities distinct.
Registry capture/download dates must not be presented as original release dates.

Both databases have **20 tables and 59 indexes**. The complete ordered
`sqlite_master` rows (`type`, `name`, `tbl_name`, `sql`) are identical. This gives
phase one a shared concrete schema baseline, including constraints and indexes.
Other historical Registry releases have not yet been profiled.

## Table inventory and candidate input mappings

These mappings are planning hypotheses from the graph models and legacy loader;
they still require validation against the selected live stage and real records.

| Table | Candidate input / responsibility |
| --- | --- |
| `analyte` | Selected metabolite groups; gene/protein identity policy still to trace |
| `source` | Source identifiers linked to exported analytes; defer `pathwayCount` calculation |
| `analytesynonym` | Source-attributed names and synonyms linked to exported analytes |
| `chem_props` | Effective curated chemistry records with provider and source ID retained |
| `metabolite_class` | Source classification terms and membership/parent relationships |
| `ontology` | Exportable HMDB ontology terms; defer `metCount` calculation |
| `analytehasontology` | Membership plus explicit legacy ontology traversal/filter rules |
| `pathway` | Pathway identifiers and source-specific annotations |
| `analytehaspathway` | Metabolite, protein, and gene pathway relationships |
| `catalyzed` | Metabolite/protein association evidence; legacy row semantics to verify |
| `reaction` | Rhea reaction records; identify fields needing derivation |
| `reaction2met` | Reaction participants mapped to exported metabolites |
| `reaction2protein` | Reaction/protein relationships mapped to exported gene/analyte IDs |
| `reaction_ec_class` | Reaction-to-EC classifications and hierarchy |
| `reaction_protein2met` | Omit, explicitly approved by the user |
| `db_version` | Export artifact version and four exclusive source-intersection JSON fields, filled by repeatable SQLite post-processing |
| `version_info` | Actual contributing source versions, with gaps reported explicitly |
| `entity_status_info` | Rerunnable post-processing summaries from stored source-bearing tables; providers discovered dynamically, legacy KEGG aliases combined |
| `pathway_similarity` | Separate rerunnable SQLite-only post-processing pass writes scope-specific zlib-compressed sparse Jaccard rows |
| `pathway_duplicates` | Same pass writes exact combined-analyte membership pairs for eligible pathways |

The backend's current `schema/RaMP_SQLite_BASE.sqlite` has only 19 tables and
adds `version_info.data_source_snapshot_ids`. It is not an exact schema template
for either comparison artifact: `reaction_protein2met` is absent. That table
contains 207,246 rows in 3.0.7 and 237,866 in the unreleased artifact. The user
identifies it as a redundant cross-product of `reaction2met` and
`reaction2protein`. No table-name references were found in the current local
RaMP-DB R code or the local `origin/ramp3.0` and `origin/ramp3.0-refmetfix` R trees.
The user approved omitting this redundant materialization and confirmed the
website's R-package data-access route. Schema parity checks should record this as an
intentional omission rather than recreating an unused table automatically.

## Chemistry and source metadata findings

| `chem_data_source` | 3.0.7 rows | Unreleased 3.0.12 rows |
| --- | ---: | ---: |
| `chebi` | 24,089 | 21,648 |
| `hmdb` | 217,776 | 217,776 |
| `lipidmaps` | 47,889 | 48,506 |
| `pubchem` | 0 | 149,017 |

Both artifacts have nine current `version_info` sources: HMDB, Reactome,
WikiPathways, KEGG, ChEBI, LipidMaps, Rhea, PFOCR, and RefMet. Neither contains a
PubChem version row, despite PubChem data in the unreleased artifact. This is a
confirmed metadata omission. Completeness for supporting resources such as
UniProt and ExPASy needs a separate dependency audit.

The existing `chem_props` columns are `ramp_id`, `chem_data_source`,
`chem_source_id`, `iso_smiles`, `inchi_key_prefix`, `inchi_key`, `inchi`, `mw`,
`monoisotop_mass`, `common_name`, and `mol_formula`. Preserve these names and
meanings initially. Existing design guidance requires source-reported InChIKeys;
do not silently substitute structure-derived keys.

RaMP-DB's `R/dataAccess.R` discovers chemical-property columns with
`PRAGMA table_info(chem_props)`. `R/rampChemPropQueries.R` treats all columns other
than three identifier/provider fields as selectable properties, and its default
query uses `SELECT *`. Later audit-column additions therefore need an explicit
data-access review; additive SQL columns are not automatically invisible to R.
The R layer already checks for `analyte.common_name`, `pathway_similarity`, and
`pathway_duplicates`, providing an existing capability-detection pattern.

## Effective stage data

In `src/qa_browser/app.py`, stage persistence stores harmonized groups and member
edges. Group records contain stage identity, a member hash, representative ID,
and a bounded `sample_member_ids` preview. Export must use complete membership
edges, never the preview list.

`_load_record_overlays()` applies published property decisions without mutating
the evidence graph. Reading raw `MetaboliteIdentifier.chem_props` alone could
therefore lose accepted corrections. Before implementation, trace how the
selected stage pins its curation snapshot and exposes effective records. Do not
combine historical group membership with an unrelated latest curation state.

Stage-generated group IDs include stage/rank information. Generate fresh RaMP
IDs deterministically from full membership rather than depending on query order.
Gene/protein, pathway, ontology, and reaction ID mappings still need inspection.

### Selected-stage observations

Read-only inspection confirmed that `stage-07-9cfef0c333530e83` is complete,
named "RaMP-ish - after s3 updates: 7. Cleanup: Remove Secondary-Source-Only
Metabolites", and was materialized on 2026-10-01.

- 147,526 non-singleton groups were counted directly.
- 540,006 membership edges were counted directly.
- Its summary reports 643,970 active identifiers and 103,964 singleton identifiers,
  yielding 251,490 metabolite groups including singletons.
- Singletons are stored implicitly in `HarmonizationStageActiveIdentifierChunk`;
  they are not individual `HarmonizedMetabolite` documents. Reconstruct them from
  active identifiers minus complete group membership and verify these totals.
- The stage records 98 metabolite property corrections and 16 ChEBI property
  corrections, with batch IDs, hashes, manifest revisions, and fingerprints in
  `curation_snapshots`. Effective export values must use that recorded state.
- `metadata_store.etl_metadata.value.registry_datasets` contains explicit source
  and derived snapshot IDs, manifest URIs, file hashes, usages, and derived-input
  references. This is a concrete input for output version metadata, including
  PubChem, without registering the stage itself.

The current stage's recorded source collection revisions can be used to check
whether its underlying graph evidence has changed. Version metadata must describe
the graph actually read, not simply today's YAML pins. Full source consistency
checks and stage-specific curation replay remain implementation prerequisites.

## Phase-one completion criteria and next work

Read-only software-architect consultation recommends four cohesive components:
a stage reader, pure table-row projection, an explicit SQLite writer, and a thin
CLI under the RaMP use case. No new ETL adapter or second resolver pass is needed.
Reuse core curation operation reduction and property projection; add narrow
recorded-snapshot replay rather than resolving the latest manifest. Extract any
shared effective-record projection from the browser instead of importing its
FastAPI application into the exporter. Retain stage/source/curation identities in
an export manifest from the first increment, independently of later SQLite audit
fields. Write to a temporary artifact and promote only after validation.

1. Inspect the selected stage, full memberships, effective curated records, and
   source metadata; record mapping coverage and gaps for every base table.
2. Trace legacy IDs, row cardinalities, null/default handling, classification
   flattening, and `reaction_protein2met` population against consumer queries.
3. Architecture consultation completed; the user explicitly authorized implementing
   the use case on 2026-10-05.
4. Implement explicit schema creation and base-table projection, then verify
   schema/index parity, referential consistency, and representative row behavior.

The intermediate artifact is not a release-ready replacement until phase-two
post-processing passes consumer checks. Defer expensive derived calculations,
new audit fields, and optional redesigns while retaining their required legacy
table/column structure. Any phase-one placeholder policy must be explicit.

No full graph export, graph mutation, or RaMP-DB code change was run. Read-only
stage-reader validation on 2026-10-05 replayed all 98 metabolite and 16 ChEBI
property decisions and reconstructed 251,490 groups / 643,970 identifiers.

## Evidence locations

- `src/use_cases/ramp/ramp.yaml`
- `src/models/metabolite_harmonization.py`
- `src/qa_browser/app.py`: stage persistence and record-property overlays
- `designs/metabolite_harmonization_design.md`: chemistry and ontology export rules
- `../ramp-backend-ncats/config/db_load_resource_config.txt`
- `../ramp-backend-ncats/main/mainSqliteDBLoad.py`
- `../ramp-backend-ncats/src/util/SQLiteDBBulkLoader.py`
- `../RaMP-DB/R/dataAccess.R`
- `../RaMP-DB/R/rampChemPropQueries.R`


## Implemented phase-one contract

Entry point: `python -m src.use_cases.ramp.build_sqlite --stage-id ID --output FILE`.
See `src/use_cases/ramp/README.md` for arguments and limitations. The user request
to implement supersedes the earlier discovery-only status.

The exporter uses `sqlite_stage.py` for immutable-batch replay and graph reads,
`sqlite_projection.py` for explicit mappings, and `sqlite_writer.py` for batched
writes, exact-row deduplication, integrity validation, and no-clobber publication.
It shares property projection with the browser through core helpers and does not
import the browser application.

The copied 3.0.7 DDL omits the approved redundant table and adds the requested
snapshot-ID column plus an embedded JSON manifest table. Base-row projection does
not calculate post-processing statistics. Pending values use NULL where allowed
and -1 for required integers; this is an explicitly incomplete artifact.

Gene and protein nodes initially retain distinct typed identities, following the
architect's conservative recommendation. Their shared symbols do not imply
identity. The user was notified of this provisional default while implementation
continued; consolidation requires an explicit identity input before asserting
legacy-equivalent gene counts.

Source data annotations remain attached to source identities. Primary name
selection is deterministic by documented provider preference. Class/ontology
ancestors follow the graph parent edges; ontology denylist policy is copied into
this use case. Scientific policy refinements remain visible follow-up work.

The 2026-10-09 product review exposed a missing part of that policy: the old
HMDB parser first selected 16 named `Source` terms, and the later denylist
removed broad categories. The current graph contains `Plant` and `Microbe` as
parent terms, while the first SQLite projection selected only Source leaves;
this dropped both query terms and admitted hundreds of specific Source leaves.
The export now defaults to `--source-ontology-policy legacy`, applying the old
Source allowlist before the denylist and expanding descendant associations to
retained parents. `--source-ontology-policy expanded` selects all non-denied
Source terms, including parents, and the manifest records the choice. Other
ontology types keep their established leaf/Health-condition-parent rule. This is an export-only
correction and does not require a graph or harmonization-stage rebuild.

Validation on 2026-10-05:

- 178 targeted tests passed across exporter, curation, property projection, and
  harmonization workbench suites.
- All 19 retained tables match the released baseline's column definitions,
  apart from the intentional `version_info.data_source_snapshot_ids` addition.
- Live read-only checks confirmed matching build fingerprint and source revisions,
  recorded curation replay, group totals, and streaming of 24,631 gene records.
- All 1,726 classification terms were checked for conflicting same-level
  ancestors; none were found. Chemistry source IDs matched their owning records.
- Code review's build-metadata consistency finding was fixed and regression
  tested. Focused re-review found no outstanding Blocking or Important issues.
- Full-stage export and downstream R-package validation remain user-run.

The browser and exporter share the source-build fingerprint calculation. The
reader rejects different ETL builds and checks revisions again before publication.
Because historical stages did not record every export collection's revision,
manual edits to those additional collections before export cannot be detected
retrospectively; this limitation is documented in the use-case README.

## UniProt-centered RaMP gene/protein identities (2026-10-05)

The initial typed-singleton export policy is superseded for new exports. The user
approved resolving genes and proteins with the same UniProt pin as the graph and
assigning sequential RAMP_G IDs to resulting accessions. This is a downstream
application projection; primary-source adapters and graph identities do not change.

`sqlite_gene_identity.py` reuses `UniProtResolver`. Exactly one canonical match
forms a shared group; unmatched or ambiguous identifiers remain typed singletons.
Ambiguous genes never bridge proteins, and associations are not expanded to all
candidate proteins. Isoform-to-canonical resolution follows the existing resolver.
Canonical accessions are queryable as UniProt-attributed source rows. Original
source IDs and provider annotations survive, and every relationship uses the same
mapping. Preferred names and sequential IDs are deterministic.

The shared UniProt parser needed explicit GeneID aliases (previously absent).
Explicit Symbol/HGNC.SYMBOL identifiers use only symbol/synonym contexts, and
Ensembl namespace spelling/version normalization happens in the resolver. RefSeq
parsing is outside this change because the inspected RaMP input has no RefSeq IDs.

The exact graph-recorded human snapshot is required, currently
`uniprot:human:2026_03`. RaMP defaults to the all-human file, including reviewed and unreviewed entries,
matching the historical backend. The generic resolver keeps its reviewed-only
default; RaMP passes its file choice explicitly. Reviewed-only remains an explicit
CLI option within the same pin. Snapshot and file checksum
must match graph metadata. Export resolution provenance is separate from the
historical graph-build metadata; no latest-version lookup is used.

Read-only discovery of the pinned reviewed file found 20,431 entries, 19,154 with
GeneID cross-references. Coverage over distinct original source IDs in the existing
base SQLite (not a new export):

| Family | Unique match | Ambiguous | Unmatched |
| --- | ---: | ---: | ---: |
| Entrez | 14,403 | 35 | 1,689 |
| Ensembl | 6,874 | 30 | 507 |
| Symbol | 984 | 5 | 33 |
| UniProt | 13,372 | 10 | 30,515 |
| Wikidata | 0 | 0 | 71 |

These matches reach 16,557 distinct canonical accessions. The survey counts
source identifiers rather than typed graph nodes, so it is not a predicted final
analyte count. Large unmatched UniProt coverage reflects the reviewed-only scope;
all-human coverage has not been inferred from this result.

Validation covers canonical/secondary/isoform and GeneID grouping, ambiguous and
unmatched preservation, deterministic order, all relationship endpoint paths,
source attribution, and snapshot/file mismatch failures. No full graph export was
run; existing SQLite reports still describe the prior build until regenerated.

Historical scope was checked in ramp-backend-ncats commit c5a61b9: the UniProt
parser loaded both human Swiss-Prot and human TrEMBL DAT files. Current backend
RegistryConfig uses uniprot-human.json.gz. Human-only and reviewed-only are
separate restrictions. The RaMP exporter therefore defaults to all-human.
Merged catalyzed pairs retain sorted distinct protein types separated by `; `;
missing source types remain Unknown. This preserves the legacy one-row primary
key and avoids choosing one conflicting annotation.

All-human audit against the same pin subsequently found:

| Family | Unique match | Ambiguous | Unmatched |
| --- | ---: | ---: | ---: |
| Entrez | 6,967 | 7,498 | 1,662 |
| Ensembl | 914 | 6,031 | 466 |
| Symbol | 79 | 912 | 31 |
| UniProt | 41,024 | 10 | 2,863 |
| Wikidata | 0 | 0 | 71 |

41,034 / 43,897 distinct UniProt IDs have human-snapshot matches. All 31,566
Rhea-attributed UniProt IDs match. HMDB has 6,074 matched / 8,292 total;
Reactome 11,687 / 12,155; WikiPathways 2,660 / 2,869. Provider counts overlap.
Unmatched does not establish non-human identity: deprecated accessions and
cross-source version differences need separate investigation. Multiple human
proteins per gene substantially increase ambiguity under the all-human scope;
the implementation retains those ambiguous genes instead of arbitrarily picking
an accession. Full output counts require a fresh stage export by the maintainer.

User-authorized live UniProt taxonomy follow-up on the 2,863 unmatched IDs:
2,792 matched primary-accession syntax and were queried. Of the full 2,863,
2,694 are confirmed non-human by returned taxonomy; 169 remain unresolved
(including 71 non-primary-format values, and inactive/missing-taxonomy records).
No queried unmatched primary accessions were returned with human taxonomy.
The 2,694 represent 6.14% of the original 43,897 UniProt source IDs. Attribution
counts overlap: HMDB 2,165, Reactome 455, WikiPathways 105. Rhea had no unmatched
IDs in the human-snapshot audit. Live taxonomy is supplementary evidence, not a
replacement of the pinned input. No species-based filtering was introduced.
Detailed results are in output_files/ramp/unmatched-uniprot-taxonomy.json and CSV.

The exporter now accepts --overwrite (API overwrite=True). It validates a fresh
temporary database and verifies the source before atomic replacement; failures
preserve the old output. Default publication still refuses existing files.
Symlinks and SQLite sidecars are rejected; close writers before replacement.

## Current gene/protein grouping policy: input-driven collapse (version 2)

The user explicitly chose legacy-compatible collapsing for the first SQLite
export, superseding the unique-match-only policy above. This is a RaMP application
identity decision, not a change to graph identifiers or primary-source adapters.

The dedicated `sqlite_gene_grouping.py` module owns `collapse_input_matches` and
policy name/version. It receives only typed input identifier -> returned canonical
UniProt matches. Each used identifier connects all its matches; overlapping sets
merge transitively. It never traverses the resolver's alias catalog. A shared
alias absent from the input cannot merge proteins, even if it appears in their
display names. An unmatched typed identifier remains a singleton. Canonical
matches need not have their own input nodes to participate in a used mapping.

`sqlite_gene_identity.py` owns pinned resolver construction and match collection.
It continues using the resolver's preferred-context matches, including primary
accession precedence; lower-priority raw aliases are not added to those matches.
`sqlite_projection.py` gathers input records, applies the completed group mapping,
assigns sequential RAMP_G IDs, preserves all original IDs and every returned
canonical accession, and remaps associations consistently. Association fan-out
is not used. The grouping function can be replaced without changing the resolver
or writer. No plugin/configuration framework was introduced.

Manifest policy version 2 identifies this decision. Ambiguous counts describe
multi-match inputs that now connect groups, not standalone analytes. Additional
counts expose multi-accession components and largest component size. IDs can be
renumbered when groups change, consistent with the accepted fresh-ID policy.

Validation includes unused aliases, used multi-match aliases, transitive closure,
typed unmatched IDs, canonical matches without input protein nodes, deterministic
ordering, source-row retention, all relationship endpoint paths, and resolver
primary-accession precedence. Full stage export remains a maintainer-run step.

### Display metadata, protein annotations, and HMDB status (2026-10-05)

The column audit found infrastructure strings in legacy display metadata,
missing protein names, and uniformly absent HMDB status. Read-only inspections
covered releases 2.3.2, 2.4.0, 3.0.7, unreleased 3.0.12, and the local export.

- `sqlite_version_metadata.py` owns readable source names, public homepages,
  release labels and recorded release dates. Distinct dataset releases remain
  labeled; identical labels collapse. Hash/dependency snapshots are described
  as snapshots/derived data with explicitly labeled download dates, never as
  upstream release numbers. Exact pins and manifest URLs remain in provenance.
- `sqlite_protein_annotations.py` reads full names independently of the resolver
  from the same SHA-verified pinned UniProt JSON. Real payload profiling found
  95,262 records with recommended names and 115,444 with submitted names, with
  no missing names. Recommended, submitted, then alternative full names are
  preferred. Annotation lookup uses primary accessions; identity reconciliation
  remains in the resolver. Canonical UniProt source rows get accession-specific
  names even when multiple accessions share a RAMP_G ID. Existing provider names
  are not relabeled as UniProt annotations. Reaction protein display names prefer
  full names over gene symbols; group display names still prefer gene symbols.
- HMDB raw XML has `<status>` (sample: HMDB0000001, quantified), but the adapter
  ignored it and the model had no slot; the exporter then hard-coded absence.
  The adapter/model now preserve `hmdb_status` on primary HMDB identifiers.
  `sqlite_hmdb_status.py` explicitly applies legacy priority
  quantified > detected > expected > predicted over active members only, then
  propagates the group's selected status to all its metabolite source rows.
  Unexpected nonempty status values fail instead of silently losing information.
- The user chose an upstream graph/stage rebuild rather than an export-time HMDB
  supplement. No graph/stage mutation or ETL was run. Re-exporting an old stage
  warns about missing status fields; the manifest records missing-field and
  status coverage counts. It cannot reconstruct statuses that were never ingested.

Validation covers metadata formatting, file-integrity rejection, recommended /
submitted / alternative protein names, accession-specific names in collapsed
protein groups, HMDB parser preservation, priority, active-member filtering and
cross-provider group propagation. The user must rebuild the graph and curated
stage, then export the new stage ID with `--overwrite`; existing SQLite/report
artifacts still describe the earlier export until regenerated.

#### Authorized existing-graph annotation repair

The user subsequently requested an in-place HMDB status backfill instead of
waiting for a rebuild. The local archive matched the recorded
`hmdb:metabolites_xml:5.0` SHA-256
`260ffcfa1e284de1390866aaaf5694b13a6dd68f7b77a6587669c6fffb59340b`.
All 217,920 primary records matched graph IDs and were patched: 3,385 quantified,
20,924 detected, 98,256 expected, and 95,355 predicted. Secondary HMDB IDs did not
receive invented source statuses. The exporter uses schema support to distinguish
an older ingest from legitimate secondary records lacking a source status.

The temporary repair verified unchanged non-status content for all 988,243
metabolite identifiers. It preserved the ETL execution identity and all stage
artifacts, added the model-derived schema field, and rebased source revision /
content fingerprint metadata for 14 previously current completed stages and two
completed runs (including embedded stage summaries). Original stage keys remain
as an explicitly recorded historical annotation-repair exception, not a claim
that the pipeline ran again. Future pipeline execution derives new stage keys.
Repair provenance is stored in `metadata_store/hmdb_status_backfill_20261005`
and stage summaries; the latter carry it into the SQLite export manifest.
Backups, input hashes, update journal, and verification are under
`output_files/ramp/hmdb-status-repair/`. No SQLite export was executed.

The separate KEGG-pathway question was checked against the existing SQLite:
all 363 KEGG pathways remain identifiable as `pathway.type = 'hmdb'` and
`pathway.pathwayCategory = 'kegg'`, with original KEGG `mapNNNNN` source IDs.
They can be reported as KEGG via HMDB without changing source provenance.

The user subsequently chose legacy KEGG attribution in the SQLite presentation:
HMDB pathway records with category `kegg` export as `pathway.type = 'kegg'` and
all corresponding metabolite/gene/protein links use
`analytehaspathway.pathwaySource = 'kegg'`. Internal graph joins still use the
actual HMDB provider. A KEGG version row explicitly says `From HMDB (...)` and
carries the HMDB metabolite snapshot ID. Other HMDB pathways remain HMDB.
This changes the next SQLite export; it does not rewrite graph provenance.

### Association source IDs and source-table contract (2026-10-05)

The user confirmed that association provenance is independent of clique
membership. The graph audit found all seven existing source-ID contracts fully
populated and equal to their identifier endpoints; protein/Rhea edges alone
lacked an explicit source identifier. The adapter previously canonicalized
secondary accessions, but reconstruction from this build's pinned files confirmed
that every retained Rhea accession already matched its primary endpoint.

The Rhea adapter now uses the pinned human UniProt primary/secondary accession
set solely for inclusion filtering. It emits the accession reported in the TSV
`ID` column as both the ProteinIdentifier ID and RheaProteinReactionEdge.source_id.
Secondary-to-primary identity reconciliation remains with GeneIdentity's resolver.
Display names for secondary input accessions use their unique resolver match,
not an arbitrary name from the collapsed group.

`sqlite_source_evidence.py` defines the eight existing edge contracts:
metabolite details use source_id; gene/protein pathway details retain gene_id /
protein_id; both Rhea edge types use source_id. Every assertion registers its
source ID and actual reporting provider against the endpoint's already assigned
RAMP ID. Evidence does not create new identity groups or additional union rules.
HMDB protein accessions retained in node properties and relationship details
also populate HMDB source rows. Resolver-added canonical rows remain attributed
to UniProt. This does not add unused resolver aliases or ingest WikiPathways
cross-references missing upstream.

Source rows are accumulated on disk with one row per (sourceId, rampId,
dataSource). Existing nonempty provider labels win over empty evidence labels;
equally ranked differing labels are selected lexically for deterministic output.
HMDB status propagates from the established metabolite group and KEGG identifier
provider attribution remains consistent with the existing exporter policy.
Known active metabolite identifiers pointing at a different curated group fail
validation; source evidence cannot reconnect a curated split.

Existing class_source_id, met_source_id, and reaction2protein.uniprot columns now
read the edge evidence rather than reconstructing from endpoints. Pathway,
ontology and catalyzed tables retain their existing schemas and RAMP endpoints;
their source identifiers enter source. This preserves lookup coverage but does
not add a new per-association audit schema.

Missing evidence fields fail clearly, with no endpoint fallback. StageReader
checks missing Rhea protein source_id before expensive projection/resolver work.
The user subsequently authorized a provenance-preserving graph repair from the
pinned inputs, described below. No SQLite export was run for this feature. The
earlier HMDB-status-only repair remains independently valid.

Validation includes raw secondary-accession retention with human filtering,
all eight source-ID contracts, multi-provider evidence, differing evidence IDs
versus endpoints, source-row deduplication/name preservation, HMDB protein IDs,
missing-field rejection, curated-split rejection, unchanged group counts, and
unchanged legacy relationship schemas.

### Existing Rhea graph repair (2026-10-05)

Reconstructed all 116,540 protein/reaction associations from the SHA-verified
`rhea:reaction_bundle:142` mapping files and `uniprot:human:2026_03` human scope.
Each stored edge matched exactly one reported accession, equal to its existing
protein endpoint. Added that accession as `source_id` on all 116,540 edges.
Verified every other edge field unchanged, with no new nodes or associations and
no changes to group membership or other source collections.

Updated the stored edge schema and recorded `rhea_source_id_backfill_20261005`
in graph metadata and compatible completed stage/run annotation-repair records.
Rhea protein edges are outside the harmonization content-fingerprint collection
set, so stage fingerprints and tracked source revisions remain unchanged.
Backups, pinned-input reconstruction hashes, and verification artifacts are in
`output_files/ramp/rhea-source-id-repair/`. This was an explicit annotation repair,
not a graph rebuild or SQLite export.
Post-repair `StageReader` validation accepted `stage-07-9cfef0c333530e83`,
including its pinned curations and unchanged-source checks. Fourteen compatible
stages and two completed runs carry the repair audit record.

### UniProt lookup enrichment (2026-10-06)

The user approved consistent UniProt enrichment across all source providers in
the SQLite extraction. Preserve original provider evidence, but attribute added
lookup aliases to UniProt rather than the legacy Reactome/Rhea attribution.
`ProteinAnnotations` reads aliases during its existing SHA-verified stream of
the resolver's exact pinned input. Only records in a group's already resolved
canonical-accession set contribute lookup rows, after group assignment.

Policy v1 includes primary and secondary UniProt accessions; gene names and
synonyms; GeneID cross-references; numeric HGNC IDs (`hgnc:<number>`, prefix once);
and Ensembl transcript IDs plus explicit GeneId/ProteinId properties. Ensembl
versions are stripped to match existing lookup conventions. Protein full/short
names, entry names, locus/ORF names, RefSeq, STRING, and other cross-references
are excluded from this first policy. They are not all equivalent identifier
types, and the resolver's broad display-name alias parser is intentionally not
used for export enrichment.

Shared aliases remain multiple lookup targets; they never feed identity collapse
or create associations. The existing source-row accumulator deduplicates within
each RAMP_G/provider. Enrichment does not populate analytesynonym. The annotation
manifest records policy version, exact input provenance, and per-namespace
accession/alias pair counts across the input file (not output-row counts).
Tests compare all non-metadata/non-source tables before and after enrichment,
verify unchanged provider rows, ambiguous lookup targets, collapsed-group
deduplication, unrepresented-record exclusion, and namespace/field allowlists.
Validation: 52 focused tests passed. A full SHA-verified parse of
`uniprot:human:2026_03` succeeded for 210,706 records, yielding 283,071 UniProt,
196,584 symbol, 39,878 GeneID, 647,000 Ensembl, and 146,604 HGNC accession/alias
pairs before restricting to represented groups. No production export was run.

### Synonym lookup attribution (2026-10-07)

The diagnostic export now includes `analytesynonym` alongside `analyte` and
`source`. Metabolite names and synonyms carry their `MetaboliteName.source`
attribution; records without that source fail export rather than borrowing every
provider on the merged identifier. HMDB `ProteinIdentifier.name`, `gene_name`,
and string synonyms originate in the HMDB protein adapter and are attributed
only to HMDB. For represented gene/protein groups, pinned UniProt recommended
full protein names and explicit gene names/synonyms from the verified resolver
file are attributed to UniProt. They never affect group construction. Identical
four-column rows collapse; names shared by providers retain one row per provider.

Legacy Rhea gene synonyms were obtained from a human UniProt download, while
Reactome gene symbols were filled from UniProt, then both were labeled with the
association provider in legacy `analytesynonym`. The new source labels describe
the actual name provider. WikiPathways' old gene symbol came from its RDF
HGNC-symbol cross-reference. The comparison report separates synonym row volume
from distinct RaMP-ID coverage by attributed source. Each example follows one
source ID across builds and shows names attributed to the selected provider;
missing synonym tables are not zero-filled. ChEBI three-star SDF names are
placed on their own `MetaboliteIdentifier` nodes and flow through this same
source-attributed synonym path after rebuilding the graph and stage.

### Pathway and ontology diagnostic extension (2026-10-07)

The default diagnostic now retains `pathway`, `ontology`, `analytehaspathway`,
and `analytehasontology` in addition to the three lookup tables. The projection
uses the same pathway and ontology logic as the full base export. It skips
unrelated catalysis and classification work. The report compares association
rows, distinct analytes, and distinct targets by stored pathway source and
analyte type; HMDB ontology links are grouped by ontology type. The ontology
link has no source column, so its HMDB attribution is inferred from the target
vocabulary. Examples match analytes by `source.sourceId` and targets by source
pathway identity or ontology type/name, with ambiguous matches labeled.
Source, synonym, and pathway example panels pair metabolite and gene/protein
records under each input where both types exist. Missing types are labeled;
ontology's gene side is structurally not applicable.

Old WikiPathways gene synonyms are not UniProt enrichment: the legacy loader
reads WikiPathways RDF `bdbHgncSymbol` into `geneInfoDictionary.common_name`
and writes it with `wiki` attribution. The adapter now retains that predicate
as `source_names` on the `GeneIdentifier` or `ProteinIdentifier` represented by
the same RDF subject. The exporter writes those explicitly attributed symbols
to `analytesynonym` with `source=wiki`; source-owned names also supply the
identifier's `source.commonName`. This is distinct from corrected UniProt
attribution for old Rhea and Reactome gene synonyms.

The pinned `wikipathways:rdf_wp:2026-09-10` archive has 2,208 Turtle files;
1,131 human pathway files yielded symbols on 12,795 distinct gene identifiers
and 2,621 distinct protein identifiers. RDF types and identifier families can
cross: a `wp:Protein` subject may carry a gene ID, and a `wp:GeneProduct`
subject may carry a protein ID. The first audited patch covered 12,324 gene
and 2,131 protein nodes; review caught the crossed cases, so a second audited
patch covered the remaining 471 gene and 490 protein nodes. Every expected ID
matched an existing graph node with WikiPathways provenance. Backups and plans
are under `output_files/ramp/wikipathways-symbol-repair/` and
`output_files/ramp/wikipathways-symbol-cross-type-repair/`. The stored collection schemas were updated
for the new field. GeneIdentifier and ProteinIdentifier are not pinned in the
harmonization stage's source revisions, so stage summaries and fingerprints
were left intact; `StageReader` validated the latest completed stage afterward.

#### Authorized current-graph ChEBI name repair

The user requested an in-place repair instead of another graph rebuild. The
current graph already held the exact pinned `chebi:three_star_sdf:2026-09-09`
name in same-ID `chem_props.common_name`; 52,956 ChEBI identifiers had such a
name, none had a ChEBI-attributed `names` entry, and all had ChEBI provenance.
The repair appended one `MetaboliteName` with source `ChEBI` and source field
`ChEBI NAME` to each of those identifiers, leaving existing names and chemistry
intact. All 52,956 patched names were verified against their same-ID chemistry.

The current pipeline rules do not use names to form metabolite groups. Eight
completed stage summaries and their completed pipeline run were rebased to the
new `MetaboliteIdentifier` revision and source-content fingerprint, with an
explicit `chebi_name_backfill` annotation. The latest stage
`stage-07-b477c68b941a5659` passed the exporter's `StageReader` validation.
Original stage keys remain as a documented annotation-repair exception; the
pipeline was not rerun. Reversible field backups, stage/run backups, and the
repair plan are under `output_files/ramp/chebi-name-repair/`. A new SQLite export
is still needed to populate ChEBI `analytesynonym` rows.

For LipidMaps `source.commonName`, prefer the identifier's `ABBREVIATION`, then
the last value in its LipidMaps `SYNONYMS` field, then `NAME`, `COMMON_NAME`, or
`SYSTEMATIC_NAME`. `LIPIDMAPS:LMFA00000001` demonstrates the distinction:
the graph has `Acetylenic acids` as a synonym and `FA 40:7;O3` as its
abbreviation; legacy RaMP stored the former, while the new preference selects
the more informative abbreviation. This rule uses only names attached
to the same identifier. Legacy RaMP also copied a name among a metabolite's
LipidMaps-reported IDs, but the current export deliberately leaves an alias ID
unnamed when it has no source-owned name, as requested. The field choice cannot
restore those aliases' name coverage.

### ChEBI ChemicalEntity chemistry fallback (2026-10-06)

The user's requested fallback reads the ChEBI `ChemicalEntity` only through the
explicit same-CURIE `ChebiChemicalEntityMetaboliteIdentifierEdge`. It writes a
ChEBI `chem_props` row only when the linked, stage-active
`MetaboliteIdentifier.chem_props` list is empty and the chemical has at least
one source chemistry field (SMILES, InChI, InChIKey, formula, average mass, or
monoisotopic mass). Existing identifier chemistry takes precedence. The ChEBI
source lookup is attributed to ChEBI; no identity groups change. The bridge
collection is added to the export's revision watch set, and the manifest records
fallback rows and groups gaining their first chemistry row.

Read-only graph audit against `stage-07-9cfef0c333530e83` found 82,733 ChEBI
links, 19,450 active identifiers with empty `chem_props`, and 18,595 of those
with usable linked chemical fields. Projecting these would add 18,595 rows to
the existing 464,654 `chem_props` rows (including 27,832 existing ChEBI rows).
They span 18,465 RAMP_C groups; 3,616 groups would gain their first chemistry
row, taking coverage from 244,254 to 247,870 distinct groups. The other 14,979
rows add ChEBI provenance to groups already having chemistry from another ID.
The separate 855 active empty identifiers whose linked chemicals have no source
chemistry remain empty. Counts are measured against the local SQLite and the
selected live graph stage; the next export will record its own actual counts.

Follow-up source-table audit found all 47,282 active ChEBI-bridged identifiers
already had a `source` row in the pre-fallback SQLite, but only 27,832 had one
attributed to ChEBI. The other 19,450 were present under other providers because
the exporter had not projected the bridge itself into `source`; 18,595 of them
also have usable linked chemistry and 855 do not. The exporter now registers a
ChEBI source row for every active same-CURIE bridge, independent of chemistry.
This changes provider attribution in the next SQLite build; it does not add new
metabolite identities or merge groups. The older `ramp-base.sqlite` report still
reflects the pre-bridge-attribution build.
