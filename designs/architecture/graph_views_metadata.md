# Graph Views In Metadata Store

## Goal

Add graph-owned live views/exports that `qa_browser` can discover and execute
generically.


The first target is `current_tdls`. The active Pharos builds now declare this
view beside `post_adapters` in their primary YAML, for example
[`pharos.yaml`](/Users/kelleherkj/IdeaProjects/IFX_ODIN/src/use_cases/pharos/pharos.yaml).

## Key Decisions

- `qa_browser` stays generic and must not hardcode graph-specific views.
- View definitions are authored in YAML and persisted into the graph's
  `metadata_store`.
- YAML uses a `graph_views` list because that is readable in config.
- Persisted metadata uses an id-keyed dictionary because ETL phases may rerun
  and need deterministic overwrite behavior.
- Aggregation happens during metadata write-back, not in a higher-level build
  orchestrator.

## Why This Design

The older GraphQL-era API exposed graph-specific endpoints such as
`current_tdls`, but that logic lived outside the graph and outside the QA
browser.

The QA browser is more maintainable precisely because it is generic. The graph
should therefore advertise its own named views, and the browser should only:

1. read advertised views from metadata
2. list them
3. execute them according to declared metadata

## YAML Authoring Shape

Add a top-level `graph_views` block to whichever YAML phase owns the view:

- main ETL YAML when the view depends on main-ETL fields
- post-processing YAML when the view depends on post-ETL fields

Example:

```yaml
graph_views:
  - id: current_tdls
    label: Current TDLs
    description: Export current protein TDL assignments
    type: export
    output_format: csv
    query_language: aql
    columns:
      - id
      - uniprot_id
      - tdl
      - uniprot_reviewed
      - uniprot_annotationScore
      - name
      - xref
      - idg_family
      - uniprot_function
      - symbol
      - ncbi_id
      - ensembl_id
      - canonical_isoform_status
      - uniprot_isoform
      - tdl_ligand_count
      - tdl_drug_count
      - tdl_go_term_count
      - tdl_generif_count
      - tdl_pm_score
      - tdl_antibody_count
    query: |
      FOR pro IN `biolink:Protein`
        RETURN {
          id: pro.id,
          uniprot_id: pro.uniprot_id,
          tdl: pro.tdl,
          uniprot_reviewed: pro.uniprot_reviewed,
          uniprot_annotationScore: pro.uniprot_annotationScore,
          name: pro.name,
          xref: pro.xref,
          idg_family: pro.idg_family,
          uniprot_function: pro.uniprot_function,
          symbol: pro.symbol,
          ncbi_id: pro.ncbi_id,
          ensembl_id: pro.ensembl_id,
          canonical_isoform_status: pro.canonical_isoform_status,
          uniprot_isoform: pro.uniprot_isoform,
          tdl_ligand_count: pro.tdl_ligand_count,
          tdl_drug_count: pro.tdl_drug_count,
          tdl_go_term_count: pro.tdl_go_term_count,
          tdl_generif_count: pro.tdl_generif_count,
          tdl_pm_score: pro.tdl_pm_score,
          tdl_antibody_count: pro.tdl_antibody_count
        }
```

First-pass scope:

- live AQL exports only
- no parameters yet
- explicit columns for deterministic CSV output

## Persisted Metadata Shape

Persist a single `metadata_store` document with key `graph_views`.

Suggested shape:

```json
{
  "_key": "graph_views",
  "value": {
    "views": {
      "current_tdls": {
        "id": "current_tdls",
        "label": "Current TDLs",
        "description": "Export current protein TDL assignments",
        "type": "export",
        "output_format": "csv",
        "query_language": "aql",
        "columns": ["id", "uniprot_id", "tdl"],
        "query": "FOR pro IN `biolink:Protein` RETURN {...}",
        "defined_in_yaml": "./src/use_cases/pharos/pharos.yaml"
      }
    }
  }
}
```

Why a dict in persisted metadata:

- post-processing can be rerun after a crash
- a rerun should overwrite the same view cleanly
- multiple ETL phases can contribute views without duplicate list entries

## Merge Semantics

Build scripts such as
[`build_pharos.py`](/Users/kelleherkj/IdeaProjects/IFX_ODIN/src/use_cases/pharos/build_pharos.py)
run primary ETL and post-processing ETL as separate `BuildGraphFromYaml(...)`
invocations. There is no current orchestration layer that combines their config
up front.

Because of that, each ETL run should:

1. read the current `metadata_store/graph_views` document
2. convert its YAML `graph_views` list into an id-keyed map
3. merge by `id`
4. overwrite any ids defined in the current run
5. preserve ids not mentioned in the current run
6. write the merged canonical document back

This makes the metadata update idempotent across repeated ETL phases.

## QA Browser Contract

The QA browser should remain generic and only assume:

1. a graph may advertise `metadata_store/graph_views`
2. each view declares enough metadata to render or download it
3. CSV-style views provide explicit columns and an AQL query

The browser should not:

- infer views from collection names
- special-case `current_tdls`
- contain Pharos-specific logic

## Relation To Provenance

This is adjacent to, but separate from, graph build provenance.

It is reasonable for `metadata_store` to eventually also expose:

- source YAML files
- contributing adapters
- build timestamps

That should remain a separate metadata concern from `graph_views`, even if both
are shown in the QA browser.

## Recommended First Implementation

1. Add `graph_views` to the graph build YAML (currently
   [`pharos.yaml`](/Users/kelleherkj/IdeaProjects/IFX_ODIN/src/use_cases/pharos/pharos.yaml)
   and [`target_graph.yaml`](/Users/kelleherkj/IdeaProjects/IFX_ODIN/src/use_cases/pharos/target_graph.yaml)).
2. Extend metadata writing so each ETL run merges its declared views into
   `metadata_store/graph_views`.
3. Update `qa_browser` to list and execute advertised views generically.

This keeps the browser generic while letting graphs expose recurring workflow
exports such as `current_tdls`.
