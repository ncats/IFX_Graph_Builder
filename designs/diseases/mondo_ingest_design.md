# Mondo Ingest Design

## Purpose
Ingest Mondo disease ontology terms and parent-child hierarchy into the target graph ETL in a scoped, incremental way.

## Source Data
- Source page: https://mondo.monarchinitiative.org/pages/download/
- File used: `mondo.json`
- Direct URL: `https://purl.obolibrary.org/obo/mondo.json`
- Download target path: `input_files/auto/mondo/mondo.json`

## Current Workflow Wiring
- Historical download rule (removed): `workflows/pharos.Snakefile` (`download_mondo`)
- Working ETL config: `src/use_cases/working.yaml`
  - `MondoDiseaseAdapter`
  - `MondoDiseaseParentEdgeAdapter`

## Scope (v1)
### Included
- Disease nodes with canonical MONDO IDs (`MONDO:NNNNNNN`) where node `type == "CLASS"`
- Disease hierarchy edges where relationship is `is_a` and both endpoints are MONDO IDs
- Disease fields:
  - `id`
  - `name`
  - `type`
  - `definition`
  - `subsets`
  - `synonyms`
  - `comments`

### Excluded (for now)
- Non-MONDO node IDs as first-class Disease nodes
- Non-`CLASS` MONDO nodes (e.g., ontology/property terms)
- Deprecated MONDO nodes (`meta.deprecated == true`)
- Non-`is_a` ontology relationships
- MONDO xrefs mapped directly to `Disease.xref` (handled by IdResolver flow)

## Data Model Mapping

### Disease node
Source: `graphs[0].nodes[*]`

- `Disease.id`:
  - from node `id`
  - normalize `http://purl.obolibrary.org/obo/MONDO_0000001` -> `MONDO:0000001`
- `Disease.name`: node `lbl`
- `Disease.type`: node `type`
- `Disease.definition`: `meta.definition.val`
- `Disease.subsets`: readable labels resolved from subset declarations
  - source subset URI from `meta.subsets[*].val`
  - resolve via declaration node (`type == "PROPERTY"`) with matching `id`
  - use first declaration comment as label when available
  - fallback to URI token (fragment or trailing path segment) when declaration label is missing
- `Disease.synonyms`: values from `meta.synonyms[*].val`
- `Disease.comments`: values from `meta.comments[*].val` (or raw string entries)

### Disease parent edge
Source: `graphs[0].edges[*]`

- Include only edges where:
  - `pred == "is_a"`
  - normalized `sub` and `obj` are both MONDO IDs
- Map to `DiseaseParentEdge`:
  - `start_node`: child Disease (`sub`)
  - `end_node`: parent Disease (`obj`)
  - `source`: `"MONDO"`

## Implemented Components
- `src/input_adapters/mondo/mondo_adapter.py`
  - `MondoDiseaseAdapter`
  - `MondoDiseaseParentEdgeAdapter`
- `src/models/disease.py`
  - expanded `Disease` fields
  - added `DiseaseParentEdge`
- `src/constants.py`
  - added `DataSourceName.Mondo`

## Version Strategy
Current adapter behavior:
- Reads graph metadata from `mondo.json` (`graphs[0].meta.basicPropertyValues`)
- Reads `graphs[0].meta.version` first (release IRI) and parses release date when available
- Uses `basicPropertyValues` only as fallback (`versionInfo` and `oboInOwl#date`)
- Emits compact version tag `vYYYY-MM-DD` when a release date is found
- Falls back to graph id if no version/date signal is available
- `download_date` from local file mtime

Observed from current file:
- version metadata may be sparse in JSON export.

Recommended follow-up:
- Preserve deterministic MONDO version and date metadata in its Registry snapshot.

## Validation Checklist
For a current snapshot, expected order-of-magnitude checks:
- MONDO disease nodes: ~30k
- MONDO->MONDO `is_a` edges: ~40k
- ID format: all nodes `MONDO:NNNNNNN`
- Edge endpoint filter: no non-MONDO IDs in `DiseaseParentEdge`
- Field presence sanity:
  - `name` on nearly all nodes
  - `synonyms`, `subsets`, `definition` on substantial subsets
  - `deprecated` present on a minority of nodes

## Open Decisions (Post-v1)
- Whether deprecated nodes should remain in graph or be filtered
- Whether to keep comments verbatim vs transform
- Whether to ingest additional relationship predicates beyond `is_a`
- Whether to persist xref mappings as explicit edges in this ingest path
