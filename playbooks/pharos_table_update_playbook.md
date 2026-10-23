# Pharos Pipeline Status Table Playbook

## Goal
Keep the Pipeline Status Table in `src/use_cases/pharos/TCRD_TODO.md` aligned with the current Pharos build and TCRD export path.

This table is intentionally **protein-facing**. It is not a literal dump of raw adapter-emitted relationship classes.

## What To Read

1. `src/use_cases/pharos/pharos.yaml`
2. `src/use_cases/pharos/tcrd.yaml`
3. The adapter source files wired in those YAMLs
4. `src/output_adapters/sql_converters/tcrd.py`
5. When needed, the SQLAlchemy schema files under `src/shared/sqlalchemy_tables/`

## Core Interpretation Rule

- Treat the table as a **Pharos/TCRD semantic status table**, not a raw edge-collection inventory.
- If a source lands as gene-based edges in the graph but is side-lifted through the `TCRDTargetResolver` into the protein-oriented Pharos/TCRD view, keep that source on the protein-facing concept row.
- Typical examples:
  - CTD belongs on `ProteinDiseaseEdge` for this table, even though the source adapter emits `GeneDiseaseEdge`.
  - WikiPathways and PathwayCommons belong on `ProteinPathwayEdge` for this table, even though their source adapters emit `GenePathwayEdge`.

## Table Semantics

Columns: `Concept | Data Sources (→ graph) | Arango Type | MySQL Tables (graph → TCRD)`

- **Concept**
  - Use the protein-facing Pharos/TCRD concept, not the adapter class name.
- **Data Sources**
  - Use one `[x]` entry per configured source that contributes to the concept.
  - Include a short note when the source is side-lifted rather than emitted directly as the protein-facing edge type.
- **Arango Type**
  - Use the primary graph class or the protein-facing edge type used by the TCRD export path.
  - Do not add parallel gene-edge rows just because the source adapter emits a gene-based edge first, unless the table is intentionally being redesigned away from protein-facing semantics.
- **MySQL Tables**
  - Mark only what the current converter and schema actually populate.
  - If there is no current converter output for that concept, use `[ ] TBD`.
  - Do not carry forward table names that only existed in older schemas.

## Update Rules

- One row per concept.
- Multiple sources can share the same row.
- Post-processing adapters from `pharos.yaml`'s `post_adapters` section stay below the separator row.
- Do not modify the Code Style TODOs section.
- Do not modify the Planned Data Sources section while updating the table.

## MySQL Verification Rules

- Verify converter output in `src/output_adapters/sql_converters/tcrd.py`.
- Verify the current schema in `src/shared/sqlalchemy_tables/pharos_tables_new.py`.
- If an old table name appears only in `pharos_tables_old.py`, do not mark it as current.
- Example:
  - `pathway_type` should not be listed as a current MySQL table because the current schema stores pathway type inline as `pathway.pwtype`.

## Recommended Workflow

1. Read the two YAML files to see what is wired.
2. Read each relevant adapter `get_all()` to see what nodes and edges it contributes.
3. Apply the protein-facing side-lift rule where appropriate.
4. Read the TCRD converter to see which concepts have downstream output.
5. Confirm any ambiguous table names against `pharos_tables_new.py`.
6. Update only the Pipeline Status Table unless the task explicitly includes adjacent checklist items.
