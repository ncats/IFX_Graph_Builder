# Pharos MySQL Summary Fields

## Goal

Reconstruct the legacy Pharos summary fields that `pharos400` needs for API/UI
compatibility, without materializing these derived counts onto graph nodes.

The graph should keep the normalized relationships:

- `ProteinDiseaseEdge`
- `DiseaseParentEdge`
- `ProteinLigandEdge`
- `ProteinLigandEdge.meets_idg_cutoff`

The summary values are easy to calculate from those relationships or from the
exported MySQL tables, and the MySQL side is the authoritative compatibility
surface for these legacy fields.

Downstream MySQL persistence fields:

- `ncats_disease.target_count`
- `ncats_disease.direct_target_count`
- `ncats_disease.maxTDL`
- `ncats_ligands.targetCount`
- `ncats_ligands.actCnt`

Do not add graph fields for:

- disease associated target count
- disease highest TDL
- ligand associated target count
- ligand cutoff-qualified associated target count
- ligand activity count
- ligand cutoff-qualified activity count

If a graph consumer needs these values, calculate them at query time from graph
edges.

## Existing Processing Pattern

Pharos already has graph post-processing for fields that affect graph semantics,
such as ligand cutoff flags and protein TDL:

Relevant existing adapters:

- `set_ligand_activity_flag.py`
  - sets `ProteinLigandEdge.meets_idg_cutoff`
  - cutoff is family-specific:
    - kinase: `7.52288`
    - ion channel: `5`
    - GPCR: `7`
    - default: `6`
- `tdl_input_adapter.py`
  - computes `Protein.tdl`
  - also computes `Protein.tdl_meta`
  - depends on `ProteinLigandEdge.meets_idg_cutoff`

The summary fields in this design should not be implemented as graph
post-processing input adapters. They should be populated during MySQL
materialization/post-processing after the graph has been exported to the TCRD
tables.

## Legacy Schema Findings

Read-only inspection of `pharos319` on 2026-05-26 showed:

`ncats_disease` columns relevant to this work:

- `maxTDL varchar(6)`
- `target_count int(11)`
- `direct_target_count int(11)`

All were populated for all `13,953` `ncats_disease` rows:

```sql
SELECT
  COUNT(*) AS ncats_disease_rows,
  SUM(target_count IS NOT NULL) AS disease_target_count_populated,
  SUM(direct_target_count IS NOT NULL) AS disease_direct_target_count_populated,
  SUM(maxTDL IS NOT NULL) AS disease_maxTDL_populated
FROM ncats_disease;
```

Result:

- rows: `13,953`
- `target_count` populated: `13,953`
- `direct_target_count` populated: `13,953`
- `maxTDL` populated: `13,953`

`ncats_ligands` columns relevant to this work:

- `actCnt int(11)`
- `targetCount int(11)`

Both were populated for all `355,932` `ncats_ligands` rows:

```sql
SELECT
  COUNT(*) AS ligand_rows,
  SUM(actCnt IS NOT NULL) AS ligand_actCnt_populated,
  SUM(targetCount IS NOT NULL) AS ligand_targetCount_populated
FROM ncats_ligands;
```

Result:

- rows: `355,932`
- `actCnt` populated: `355,932`
- `targetCount` populated: `355,932`

Live `pharos400` also already has these columns, plus `provenance` and
`ncats_disease.novelty`. No schema migration is needed for the legacy-compatible
summary fields currently present in `pharos400`.

## Disease Count Semantics

The old migration source at
`https://github.com/ncats/pharos-ETL/blob/main/tcrd2pharos/migrations/20210713101200_createTargetCountColumn.ts`
was not accessible from this environment; GitHub returned `404` for the repo
and raw file. The local `pharos-graphql-server` checkout still documents and
uses the legacy semantics:

- `src/config/diseaseModelConfig.ts`
  - `target_count`: "Count of proteins associated with the disease, or child diseases"
  - `direct_target_count`: "Count of proteins documented to be directly associated with the disease, not a child disease"
- `src/models/target/targetList.ts`
  - disease target-list queries include the selected disease MONDO ID or any
    `ancestry_mondo.oid` row where `ancestor_id` is the selected disease MONDO ID

`ncats_d2da.direct` describes disease hierarchy: `1` links a target association
to the named disease, while the broader target count rolls up associations from
descendant diseases. It does not describe whether the source assertion named a
gene, transcript, or protein. Pharos deliberately resolves gene and transcript
evidence to canonical protein targets, as it does for CTD disease associations
and TIGA traits. DISEASES gene assessments therefore do not by themselves
prevent TCRD export. The graph retains DISEASES source subject IDs in edge
details; the current MySQL `disease.did` stores the source disease ID, not the
original gene or protein subject ID.
Negated DISEASES assertions remain in graph details but are omitted from TCRD
`disease` and `ncats_d2da` rows, whose schema has no assertion-polarity field.

`pharos319.ncats_disease.direct_target_count` matches the distinct directly
associated target count through `ncats_d2da.direct = 1`:

```sql
SELECT
  SUM(calc_direct = direct_target_count) AS direct_matches,
  COUNT(*) AS rows_checked
FROM (
  SELECT
    n.id,
    n.direct_target_count,
    COUNT(DISTINCT d.protein_id) AS calc_direct
  FROM ncats_disease n
  LEFT JOIN ncats_d2da x
    ON x.ncats_disease_id = n.id
   AND x.direct = 1
  LEFT JOIN disease d
    ON d.id = x.disease_assoc_id
  GROUP BY n.id
) q;
```

Result:

- `direct_matches`: `13,953 / 13,953`

`target_count` is broader than direct disease-association rows. Comparing it to
direct/all `ncats_d2da` rows only matched `12,696 / 13,953` diseases, and top
legacy diseases show large rollups:

- `neoplasm`: `target_count = 16,382`, `direct_target_count = 24`
- `cancer`: `target_count = 16,039`, `direct_target_count = 901`
- `nervous system disorder`: `target_count = 14,453`, `direct_target_count = 53`

This confirms `target_count` is an ancestry-aware rollup over descendant MONDO
disease associations, not just direct `ProteinDiseaseEdge` rows. The MySQL
post-processing step should reproduce the `ancestry_mondo` behavior: selected
disease plus all descendants in the MONDO disease hierarchy.

### Proposed Disease MySQL Fields

Populate these fields on `ncats_disease`:

- `target_count`
  - ancestry-aware distinct associated target count
  - selected disease plus descendants from `ancestry_mondo`
- `direct_target_count`
  - distinct directly associated target count
  - direct `ncats_d2da -> disease.protein_id` only
- `maxTDL`
  - highest target development level among direct associated targets

TDL ordering should be:

```text
Tclin > Tchem > Tbio > Tdark
```

`maxTDL` should use direct associated targets. The legacy migration supplied by
Keith confirms it joins `ncats_disease -> ncats_d2da -> disease -> t2tc ->
target` and takes the max TDL rank. There is no MONDO ancestry join in that
calculation. The ranking intentionally mirrors the legacy migration; live
`pharos400.target.tdl` currently uses only `Tclin`, `Tchem`, `Tbio`, and
`Tdark`.

### Disease SQL Shape

Direct count:

```sql
UPDATE ncats_disease n
JOIN (
  SELECT
    x.ncats_disease_id,
    COUNT(DISTINCT d.protein_id) AS direct_target_count
  FROM ncats_d2da x
  JOIN disease d
    ON d.id = x.disease_assoc_id
  WHERE x.direct = 1
    AND d.protein_id IS NOT NULL
  GROUP BY x.ncats_disease_id
) q
  ON q.ncats_disease_id = n.id
SET n.direct_target_count = q.direct_target_count;
```

Ancestry-aware count should include direct targets for the disease and direct
targets for all descendant diseases. It should use the MySQL `ancestry_mondo`
closure after MySQL ontology post-processing populates it:

```sql
UPDATE ncats_disease n
JOIN (
  SELECT
    parent.id AS ncats_disease_id,
    COUNT(DISTINCT d.protein_id) AS target_count
  FROM ncats_disease parent
  LEFT JOIN ncats_disease child
    ON child.mondoid = parent.mondoid
    OR child.mondoid IN (
      SELECT a.oid
      FROM ancestry_mondo a
      WHERE a.ancestor_id = parent.mondoid
    )
  LEFT JOIN ncats_d2da x
    ON x.ncats_disease_id = child.id
  LEFT JOIN disease d
    ON d.id = x.disease_assoc_id
   AND d.protein_id IS NOT NULL
  GROUP BY parent.id
) q
  ON q.ncats_disease_id = n.id
SET n.target_count = q.target_count;
```

This query is intentionally shown as semantics rather than final tuned SQL. The
implementation should avoid correlated subqueries if performance is poor; a
temporary/staging table for `(ancestor_ncats_disease_id, descendant_ncats_disease_id)`
may be preferable.

`maxTDL`:

```sql
UPDATE ncats_disease n
JOIN (
  SELECT
    n2.id,
    MAX(CASE
      WHEN target.tdl = 'Tclin' THEN 4
      WHEN target.tdl = 'Tchem' THEN 3
      WHEN target.tdl = 'Tbio' THEN 2
      ELSE 1
    END) AS tempTDL
  FROM ncats_disease n2
  JOIN ncats_d2da x
    ON n2.id = x.ncats_disease_id
  JOIN disease d
    ON x.disease_assoc_id = d.id
  JOIN t2tc
    ON d.protein_id = t2tc.protein_id
  JOIN target
    ON t2tc.target_id = target.id
  GROUP BY n2.id
) q
  ON q.id = n.id
SET n.maxTDL = CASE
  WHEN q.tempTDL = 4 THEN 'Tclin'
  WHEN q.tempTDL = 3 THEN 'Tchem'
  WHEN q.tempTDL = 2 THEN 'Tbio'
  ELSE 'Tdark'
END;
```

## Ligand Count Semantics

Read-only `pharos319` checks show:

- `ncats_ligands.actCnt` equals `COUNT(ncats_ligand_activity.id)`
- `ncats_ligands.targetCount` equals `COUNT(DISTINCT ncats_ligand_activity.target_id)`

Sample/top rows:

- acetazolamide: `actCnt = 1,886`, `targetCount = 13`
- staurosporine: `actCnt = 1,196`, `targetCount = 265`
- vorinostat: `actCnt = 985`, `targetCount = 15`

A 10,000-row grouped check matched both fields exactly:

```sql
SELECT
  SUM(calc_act = actCnt) AS act_matches,
  SUM(calc_target = targetCount) AS target_matches,
  COUNT(*) AS rows_checked
FROM (
  SELECT
    l.id,
    l.actCnt,
    l.targetCount,
    COUNT(a.id) AS calc_act,
    COUNT(DISTINCT a.target_id) AS calc_target
  FROM ncats_ligands l
  LEFT JOIN ncats_ligand_activity a
    ON a.ncats_ligand_id = l.id
  GROUP BY l.id
  LIMIT 10000
) q;
```

Result:

- `act_matches`: `10,000 / 10,000`
- `target_matches`: `10,000 / 10,000`

### Proposed Ligand MySQL Fields

Populate these fields on `ncats_ligands`:

- `actCnt`
  - count of all activity details emitted to `ProteinLigandEdge.details`
  - exported rows in `ncats_ligand_activity`
- `targetCount`
  - count of distinct targets connected by `ProteinLigandEdge`
  - exported rows in `ncats_ligand_activity`

Do not persist cutoff-qualified ligand counts in `pharos400` unless new columns
are explicitly added. If the graph/API needs cutoff-qualified counts, calculate
them from `ProteinLigandEdge.meets_idg_cutoff` at query time.

### Ligand SQL Shape

```sql
UPDATE ncats_ligands l
JOIN (
  SELECT
    a.ncats_ligand_id,
    COUNT(a.id) AS actCnt,
    COUNT(DISTINCT a.target_id) AS targetCount
  FROM ncats_ligand_activity a
  GROUP BY a.ncats_ligand_id
) q
  ON q.ncats_ligand_id = l.id
SET
  l.actCnt = q.actCnt,
  l.targetCount = q.targetCount;
```

## MySQL Persistence Plan

Do not update the graph models or SQL converters to copy summary values from
graph nodes.

Instead, add MySQL-side post-processing after table export and ontology ancestry
post-processing:

- populate `ancestry_mondo` from `mondo_parent`
- populate `ncats_disease.target_count`
- populate `ncats_disease.direct_target_count`
- populate `ncats_disease.maxTDL`
- populate `ncats_ligands.actCnt`
- populate `ncats_ligands.targetCount`

The natural location is `MysqlOutputAdapter.do_post_processing()` or helper
methods called from there, because it already builds ancestry tables.

### Graph vs MySQL Calculation Boundary

Summary fields are a downstream materialization concern, not source graph
content.

Direct graph queries remain simple when needed:

- direct disease target count: count inbound `ProteinDiseaseEdge`
- direct disease highest TDL: max over inbound `ProteinDiseaseEdge -> Protein.tdl`
- ligand activity/target counts: count inbound `ProteinLigandEdge` and details

But Pharos compatibility should be reconstructed from MySQL tables because:

- the API/UI consumes the MySQL fields
- `target_count` depends on `ancestry_mondo`
- `pharos319` compatibility is defined by MySQL migrations and SQL query
  semantics, not by graph node denormalization

Current inspection notes:

- `pharos319` uses an older MONDO snapshot; for `MONDO:0005070` / neoplasm,
  `ancestry_mondo` has `792` descendant `oid` values.
- the current `ifxdev` `pharos` graph has `39,765` MONDO `DiseaseParentEdge`
  rows and graph traversal for `MONDO:0005070` finds `4,397` descendants.
- live `pharos400` has `26,660` `mondo` rows and `39,765` `mondo_parent` rows,
  but `ancestry_mondo` is currently empty, so MySQL post-processing must run
  before MySQL-side ancestry parity can be checked.

The safest implementation is therefore to avoid graph summary persistence and
populate only the MySQL legacy columns from MySQL tables after export.

## Validation Plan

Graph validation:

- no new graph summary fields should be present
- existing graph relationships remain queryable:
  - `ProteinDiseaseEdge`
  - `DiseaseParentEdge`
  - `ProteinLigandEdge`
  - `ProteinLigandEdge.meets_idg_cutoff`

MySQL validation after export to `pharos400_working`:

```sql
SELECT
  COUNT(*) AS disease_rows,
  SUM(target_count IS NOT NULL) AS target_count_populated,
  SUM(direct_target_count IS NOT NULL) AS direct_target_count_populated,
  SUM(maxTDL IS NOT NULL) AS maxTDL_populated
FROM ncats_disease;
```

```sql
SELECT
  COUNT(*) AS ligand_rows,
  SUM(actCnt IS NOT NULL) AS actCnt_populated,
  SUM(targetCount IS NOT NULL) AS targetCount_populated
FROM ncats_ligands;
```

Spot-check disease direct counts:

```sql
SELECT
  SUM(calc_direct = direct_target_count) AS direct_matches,
  COUNT(*) AS rows_checked
FROM (
  SELECT
    n.id,
    n.direct_target_count,
    COUNT(DISTINCT d.protein_id) AS calc_direct
  FROM ncats_disease n
  LEFT JOIN ncats_d2da x
    ON x.ncats_disease_id = n.id
   AND x.direct = 1
  LEFT JOIN disease d
    ON d.id = x.disease_assoc_id
  GROUP BY n.id
) q;
```

Spot-check ligand counts:

```sql
SELECT
  SUM(calc_act = actCnt) AS act_matches,
  SUM(calc_target = targetCount) AS target_matches,
  COUNT(*) AS rows_checked
FROM (
  SELECT
    l.id,
    l.actCnt,
    l.targetCount,
    COUNT(a.id) AS calc_act,
    COUNT(DISTINCT a.target_id) AS calc_target
  FROM ncats_ligands l
  LEFT JOIN ncats_ligand_activity a
    ON a.ncats_ligand_id = l.id
  GROUP BY l.id
) q;
```

## Open Questions

- Decide whether cutoff-qualified ligand counts are needed by a MySQL-backed
  API/UI surface. If yes, add explicit `pharos400` columns and a separate
  migration/design note. Do not overload `actCnt` or `targetCount`.
- Tune the ancestry-expanded `target_count` SQL on `pharos400_working`; a
  staging table may be needed for acceptable runtime.
