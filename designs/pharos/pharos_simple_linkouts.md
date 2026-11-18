# Pharos Simple Linkouts

## Goal

Add Pharos target-level external linkouts to the graph, then export them to the
Pharos MySQL `affiliate` and `extlink` tables.

This is Pharos presentation data. Wire it into the Pharos build path only, not
the shared target graph build.

## Graph Shape

Use one provider node and one protein-to-provider edge per concrete linkout:

```text
Protein --ProteinExternalLinkEdge { url, source_id, source_id_type }--> ExternalLinkProvider
```

Proposed models:

- `ExternalLinkProvider`
  - `id`: stable key, for example `external_link_provider:PubChem`
  - `source`: legacy/source key, for example `PubChem`
  - `display_name`
  - `description`
- `ProteinExternalLinkEdge`
  - `start_node`: resolved `Protein`
  - `end_node`: `ExternalLinkProvider`
  - `url`: concrete outbound URL
  - `source_id`: identifier used by the provider URL or list
  - `source_id_type`: `uniprot`, `symbol`, `ensembl_gene`, etc.

Downstream MySQL mapping:

- `ExternalLinkProvider` -> `affiliate`
- `ProteinExternalLinkEdge` -> `extlink`
- `affiliate.link_count` stores the number of linkout edges for that provider so
  the API can sort providers without relying on hidden `affiliate.id` ordering.

When a graph edge has multiple `details` entries, the MySQL export writes one
`extlink` row using the detail with the lexicographically smallest `source_id`.
The graph keeps all details; MySQL gets one deterministic URL per
protein/provider edge.

## MySQL Schema Decision

Legacy Pharos stores provider metadata in `affiliate` and per-protein URLs in
`extlink`. The old schema used an enum for `extlink.source`; this caused PubChem
rows to be stored with a blank source because `PubChem` was missing from the
enum.

Current IFX_ODIN should use an explicit relationship instead:

- `affiliate.source` is unique.
- `extlink.source` is `VARCHAR(255) NOT NULL`.
- `extlink.source` has a foreign key to `affiliate.source`.

Legacy `pharos319` correction already applied on 2026-05-19:

- Added `PubChem` to the live `extlink.source` enum.
- Updated `20,412` blank-source PubChem rows to `source = 'PubChem'`.

## Adapter vs Resolver Boundary

Adapters and linkout generators should preserve the source identifier they
actually have. They should not perform cross-identifier reconciliation.

Examples:

- GlyGen emits `uniprot_canonical_ac` as provided, including isoform suffixes.
- Dark Kinome emits the source gene symbol from the page.
- RESOLUTE emits source gene symbols and any returned NextProt/Ensembl protein
  identifiers.
- TIGA reuses the Ensembl gene ID already present in TIGA association details.

Mapping those identifiers to canonical Pharos `Protein` nodes belongs in
resolver configuration, conversion logic, or a documented post-resolution
materialization step.

## Source Plan

| Source | Status | How We Get Coverage | Source ID | URL Rule | Graph Materialization |
| --- | --- | --- | --- | --- | --- |
| PubChem | Implement | Generate from every Pharos `Protein` with `uniprot_id` | UniProt accession | `https://pubchem.ncbi.nlm.nih.gov/protein/<uniprot_id>` | Create provider node and one edge per protein with `uniprot_id` |
| ARCHS4 | Implement | Generate from every Pharos `Protein` with `sym` | Gene symbol | `https://archs4.org/gene/<symbol>` | Create provider node and one edge per protein with `sym`; no first-pass probing |
| GlyGen | Implement | Fetch current GlyGen API list: search human proteins, then CSV download via `data/list_download` | GlyGen `uniprot_canonical_ac` | `https://glygen.org/protein/<base UniProt accession>` | Create edges for GlyGen list entries that resolve to Pharos proteins |
| Dark Kinome | Implement | Fetch and parse `https://darkkinome.org/data` Kinase List | Gene symbol | `https://darkkinome.org/kinase/<symbol>` | Create edges for source symbols that resolve to Pharos proteins |
| RESOLUTE | Implement | Query `https://re-solute.eu/api/graphql` for `genes(condition: {isSlc: true})` | Gene symbol; optional NextProt/Ensembl protein IDs | `https://re-solute.eu/knowledgebase/gene/<symbol>` | Create edges for RESOLUTE genes that resolve to Pharos proteins |
| TIGA | Implement | Derive from existing TIGA ingest; no separate download | Ensembl gene ID | `https://unmtid-shinyapps.net/shiny/tiga/?gene=<ensg>` | Create edges after TIGA protein endpoints are resolved/materialized |
| LinkedOmicsKB | Implement | Fetch full current list from `https://kb.linkedomics.org/data/list/gene` | Source gene symbol, preserved exactly | `https://kb.linkedomics.org/gene/<source symbol>` | Create edges for source symbols that resolve to Pharos proteins |
| ProKinO | Defer | No reliable current endpoint/artifact found | Unknown/currently unavailable | Legacy links appear broken | Do not emit |
| GENEVA | Defer | Current service appears unavailable | Gene symbol if service returns | `https://genevatool.org/gene_table?gene_name=<symbol>` | Do not emit while service is unavailable |
| Reactome | Defer | Stable target/entity pages need Reactome stable IDs not currently on protein nodes | Reactome stable ID | `https://reactome.org/content/detail/<Reactome stable ID>` | Keep pathway-level Reactome links only |
| ClinGen | Defer as simple linkout | Has list/report datasets, not protein/gene landing pages | HGNC gene ID/symbol in reports | Report/list URLs vary by curation type | Treat as future clinical curation ingest, not an `extlink` provider |

## Acquisition Notes

### Static From Protein Fields

These providers do not need upstream coverage files for the first pass:

- PubChem from `Protein.uniprot_id`.
- ARCHS4 from `Protein.sym`.

These should be generated after protein nodes are available in the Pharos graph.

### Source Lists

These should be captured as Registry snapshots, then read
by a Pharos-only linkout adapter:

- GlyGen human protein CSV.
- Dark Kinome kinase list page.
- RESOLUTE SLC GraphQL result.
- LinkedOmicsKB full gene JSON list.

Each source-list adapter should emit source identifiers as-is and rely on the
resolver/materialization path to connect to Pharos proteins.

### Existing Ingest-Derived

TIGA should not be fetched again for linkouts. It should be derived from the
existing TIGA ingest after `ProteinGwasTraitEdge` endpoints are resolved.

## Provider Details

Downloaded source-list counts from the Pharos Snakemake workflow on 2026-05-19:

| Source | Useful Current Endpoint / Artifact | Downloaded Count | Notes |
| --- | --- | ---: | --- |
| GlyGen | `POST https://api.glygen.org/protein/search_simple/` with `{"term_category": "organism", "term": "human"}`, then `POST /data/list_download/` CSV | `20,659` rows | matches discovery count confirmed by user |
| Dark Kinome | `https://darkkinome.org/data` | `162` rows | legacy had `161` rows |
| RESOLUTE | `https://re-solute.eu/api/graphql` | `522` rows | current API reports `522` SLC genes; legacy had `451` rows |
| LinkedOmicsKB | `https://kb.linkedomics.org/data/list/gene` | `19,701` rows | legacy had `19,183` rows |

Generated or existing-ingest provider counts:

| Source | Artifact / Rule | Count Reference |
| --- | --- | ---: |
| TIGA | `input_files/auto/tiga/tiga_gene-trait_stats.tsv` and provenance file | legacy had `19,816` rows over `18,005` proteins |
| ARCHS4 | generated from Pharos symbols | legacy had `20,238` rows, matching proteins with populated symbols |
| PubChem | generated from Pharos UniProt IDs | legacy had `20,412` rows, matching proteins with UniProt IDs |

## Deferred Sources

| Source | Reason |
| --- | --- |
| ProKinO | Site/downloads appear unreliable and old links appear broken. Revisit only with a stable current artifact or endpoint. |
| GENEVA | `genevatool.org` did not return useful content; HTTP returned `503` in discovery. |
| Reactome | Target-level stable pages require Reactome stable IDs. Current graph has pathway-level Reactome links, which should remain the Reactome representation for now. |
| ClinGen | Current data is list/report oriented: Gene-Disease Validity, Dosage Sensitivity, Actionability, Variant Pathogenicity. It is a good future clinical curation source, but not a simple protein linkout provider. |

## Implementation Direction

Build a Pharos-only graph linkout adapter or post-resolution materialization step
that emits:

1. `ExternalLinkProvider` nodes for implemented providers.
2. `ProteinExternalLinkEdge` edges for resolved proteins and concrete URLs.
3. Enough `source_id` detail to explain how each URL was generated.

The first pass should keep linkout logic minimal:

- Generate PubChem and ARCHS4 from existing protein fields.
- Use current source lists for GlyGen, Dark Kinome, RESOLUTE, and LinkedOmicsKB.
- Derive TIGA from existing TIGA graph data.
- Skip deferred providers entirely.

After validation, export provider nodes and linkout edges to `affiliate` and
`extlink`.
