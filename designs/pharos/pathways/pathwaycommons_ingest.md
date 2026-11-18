# PathwayCommons Ingest Design

## Background

PathwayCommons (PC) is an aggregator of curated biological pathway data from multiple
upstream sources (PID, INOH, HumanCyc, NetPath, Panther, PathBank, Biofactoid). The
old TCRD loader (`loaders/load-PathwayCommons.py`, Steve Mathias) downloaded a UniProt
GMT file from PathwayCommons, filtered out sources already covered elsewhere (Reactome,
KEGG, WikiPathways), and inserted rows into the `pathway` table with
`pwtype = 'PathwayCommons: <source>'`. This design replaces that approach using the
IFX_ODIN ETL pattern.

## Inputs

Downloaded via `workflows/pharos.Snakefile` into `input_files/auto/pathwaycommons/`:

- `pc-hgnc.gmt.gz` — HGNC-symbol GMT from `https://download.baderlab.org/PathwayCommons/PC2/v14/pc-hgnc.gmt.gz`
- `pathwaycommons_version.tsv` — derived from `original_datasources.txt` header line:
  - `version`: PC major version number (e.g. `14`)
  - `version_date`: release date parsed from the `#CPATH2: PC version N DD Mon YYYY` header (e.g. `2024-05-21`)

For a new complete release, capture and pin a new Registry snapshot from
`https://download.baderlab.org/PathwayCommons/PC2/`.

## File Format (surveyed 2026-03-16)

GMT format: one pathway per line, tab-separated:

```
pathway_id \t name: ...; original_datasource: ...; organism: 9606; idtype: hgnc.symbol \t GENE1 \t GENE2 \t ...
```

- Column 0: pathway URI (format varies by source — see below)
- Column 1: metadata string with semicolon-separated key-value pairs
- Columns 2+: HGNC gene symbols (all human, organism is always 9606)

### Pathway ID formats by source

| Source     | Example ID                                               |
|------------|----------------------------------------------------------|
| biofactoid | `biofactoid:00e3ac12-8c96-45df-a59c-e212724962b6`        |
| humancyc   | `humancyc:Pathway4378797`                                |
| inoh       | `inoh:id1001348414_Calcineurin_activation_signaling`     |
| netpath    | `netpath:Pathway_Alpha6Beta4Integrin`                    |
| panther    | `https://identifiers.org/panther.pathway:P00001`         |
| pathbank   | `http://bioregistry.io/pathbank:SMP0000001`              |
| pid        | `pid:pid_10168`                                          |
| kegg       | `http://bioregistry.io/kegg.pathway:hsa00010`            |
| reactome   | `http://bioregistry.io/reactome:R-HSA-1059683`           |

### Observed shape (v14, 2024-05-21)

| Source     | Pathways | Notes                              |
|------------|----------|------------------------------------|
| reactome   | 2160     | **excluded** — loaded separately   |
| pathbank   | 834      | included (new in v14)              |
| inoh       | 576      | included                           |
| humancyc   | 352      | included                           |
| biofactoid | 263      | included (new in v14)              |
| pid        | 223      | included                           |
| panther    | 150      | included                           |
| kegg       | 83       | included — KEGG not loaded elsewhere in Pharos |
| netpath    | 27       | included                           |
| **total kept** | **2508 pathways** | after deduplication              |
| **edges**  | **44,703** | protein→pathway edges            |

Note: pathway count after dedup (2508) is slightly higher than unique pathway IDs (2425)
because a small number of pathway IDs appear with identical names across sources in the
metadata. Deduplication is by pathway_id.

## Scope

- Pathway nodes for all non-excluded sources
- ProteinPathwayRelationship edges linking proteins to pathways via HGNC symbol resolution
- **No pathway hierarchy** — `pathways.txt.gz` contains hierarchy data but it is not loaded
  in this design (hierarchy would require `PathwayParentEdge` support for PC pathway IDs)

## Mapping to the data model

### Nodes

- **Pathway** (`src/models/pathway.py`)
  - `id`: compact `prefix:localid` form — URL prefixes (`http://bioregistry.io/`, `https://identifiers.org/`) are stripped (e.g. `pathbank:SMP0000001`, `kegg.pathway:hsa00010`, `pid:pid_10168`)
  - `source_id`: local identifier only — everything after the last `:` (e.g. `SMP0000001`, `hsa00010`, `pid_10168`)
  - `type`: `"PathwayCommons"` — consistent with `"Reactome"`, `"WikiPathways"`
  - `original_datasource`: sub-source name (e.g. `"pid"`, `"inoh"`, `"kegg"`) — new field added to `Pathway` model
  - `name`: value of `name:` key in metadata column
  - `url`: pathway URI when it starts with `http`/`https`; `None` otherwise

### Edges

- **ProteinPathwayRelationship** (`src/models/pathway.py`)
  - `start_node`: `Protein(id="Symbol:<symbol>")` — resolved to canonical UniProt protein by `tcrd_targets` resolver
  - `end_node`: `Pathway(id=pathway_uri)`
  - `source`: `"PathwayCommons"`

## Adapters (`src/input_adapters/pathwaycommons/pathwaycommons_pathways.py`)

- **`PathwayCommonsBaseAdapter`** — extends `FlatFileAdapter`
  - `get_original_datasource_name()` → `DataSourceName.PathwayCommons`
  - `get_version()` → reads `pathwaycommons_version.tsv` for `version` and `version_date`
  - `_iter_parsed_lines()` — yields `(clean_id, source_id, name, original_datasource, url, genes)` tuples;
    excludes `reactome` and `wikipathways` only (kegg is included)

- **`PathwayCommonsPathwayAdapter`**
  - Input: `pc-hgnc.gmt.gz`
  - Output: `Pathway` nodes (deduplicated by clean_id)

- **`PathwayCommonsProteinPathwayEdgeAdapter`**
  - Input: `pc-hgnc.gmt.gz`
  - Output: `ProteinPathwayRelationship` edges
  - Emits `Protein(id="Symbol:<symbol>")` as `start_node`; resolver maps to canonical protein

## Model changes

- `src/models/pathway.py`: added `original_datasource: Optional[str] = None` to `Pathway`

## Constants

Add `PathwayCommons = "PathwayCommons"` to `DataSourceName` in `src/constants.py`.

## Version strategy

PathwayCommons releases are versioned numerically (v12, v14, etc.). The download URL is
hardcoded to the current release. The Snakefile rule parses the version number and release
date from the `#CPATH2:` comment line in `original_datasources.txt` and writes them to
`pathwaycommons_version.tsv`. For future releases, capture this metadata with the Registry snapshot.

## YAML wiring

Add to both `src/use_cases/pharos/pharos.yaml` and `src/use_cases/pharos/target_graph.yaml`:

```yaml
- import: ./src/input_adapters/pathwaycommons/pathwaycommons_pathways.py
  class: PathwayCommonsPathwayAdapter
  kwargs:
    file_path: ./input_files/auto/pathwaycommons/pc-hgnc.gmt.gz
    version_file_path: ./input_files/auto/pathwaycommons/pathwaycommons_version.tsv

- import: ./src/input_adapters/pathwaycommons/pathwaycommons_pathways.py
  class: PathwayCommonsProteinPathwayEdgeAdapter
  kwargs:
    file_path: ./input_files/auto/pathwaycommons/pc-hgnc.gmt.gz
    version_file_path: ./input_files/auto/pathwaycommons/pathwaycommons_version.tsv
```

## Comparison with old TCRD loader

| Aspect | Old loader | This design |
|--------|------------|-------------|
| Gene ID source | UniProt accessions (`pc.UniProt.hgnc.gmt`) | HGNC symbols (`pc-hgnc.gmt.gz`) |
| Gene→Protein mapping | Direct MySQL lookup by UniProt | `tcrd_targets` resolver via `HGNC:` prefix |
| Sources kept | humancyc, inoh, netpath, panther, pid | Same + kegg, pathbank, biofactoid (new in v14) |
| Sources excluded | reactome, kegg, wikipathways | reactome only |
| Sub-source stored as | `pwtype = "PathwayCommons: pid"` (compound) | `type = "PathwayCommons"` + `original_datasource = "pid"` (separate fields) |
| PC version | pc11 (2015-era) | v14 (2024-05-21) |
| Storage | MySQL `pathway` table | ArangoDB graph → TCRD via existing pathway converter |

## Notes

- The file contains only human entries (organism: 9606); no species filtering needed.
- HGNC symbols are used as-is; wrap as `EquivalentId(id=symbol, type=Prefix.Symbol)` to
  produce `Symbol:TP53` etc. for the resolver. Note: `Prefix.HGNC` is for numeric HGNC IDs
  (e.g. `HGNC:6023`), not gene symbols.
- The `pathways.txt.gz` hierarchy file (pathway → sub-pathways) is not currently loaded.
  If pathway hierarchy support is added for PathwayCommons, it would emit `PathwayParentEdge`
  records from the `DIRECT_SUB_PATHWAY_URIS` column.
- The TCRD MySQL pathway converter (`src/input_adapters/pharos_arango/tcrd/pathway.py`)
  already handles `ProteinPathwayRelationship`; PathwayCommons edges will flow through it
  automatically once in the graph.
