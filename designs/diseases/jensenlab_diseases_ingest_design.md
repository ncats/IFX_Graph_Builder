# JensenLab DISEASES Ingest Design

## Goal
Add a modern JensenLab DISEASES ingest for protein-disease associations using the current live Jensen download surface, without recreating the old TCRD preprocessing pipeline.

## Source Files

Primary ingest inputs:
- `https://download.jensenlab.org/human_disease_knowledge_filtered.tsv`
- `https://download.jensenlab.org/human_disease_experiments_filtered.tsv`
- `https://download.jensenlab.org/human_disease_textmining_filtered.tsv`

Useful secondary files:
- `https://download.jensenlab.org/human_disease_benchmark.tsv`
- `https://download.jensenlab.org/human_disease_integrated_full.tsv`
- `https://download.jensenlab.org/human_disease_knowledge_full.tsv`
- `https://download.jensenlab.org/human_disease_experiments_full.tsv`
- `https://download.jensenlab.org/human_disease_textmining_full.tsv`

Out of scope for first pass:
- `disease_textmining_mentions.tsv`
  This appears to be mention-level support data rather than the right first-pass association file.

## Version Strategy

- Download the three filtered files into `input_files/auto/jensenlab/`.
- Record a shared `diseases_version.tsv`.
- Use the latest `Last-Modified` date across the three filtered URLs as `version_date`.
- Leave `version` empty unless Jensen exposes a stable release string later.

## Observed Payload Shape

All three filtered files are headerless TSVs with 7 columns.

### `human_disease_knowledge_filtered.tsv`
Columns:
1. protein identifier
2. gene symbol / label
3. disease identifier
4. disease name
5. source label
6. evidence text
7. confidence

Sample:
`ENSP00000001146  CYP26B1  DOID:2340  Craniosynostosis  UniProtKB-KW  CURATED  4`

### `human_disease_experiments_filtered.tsv`
Columns:
1. protein identifier
2. gene symbol / label
3. disease identifier
4. disease name
5. source label
6. evidence text
7. confidence

Sample:
`ENSP00000000233  ARF5  DOID:11714  Gestational diabetes  TIGA  MeanRankScore = 36  1.004`

### `human_disease_textmining_filtered.tsv`
Columns:
1. protein-like identifier
2. gene symbol / label
3. disease identifier
4. disease name
5. z-score
6. confidence
7. source URL

Sample:
`18S_rRNA  18S_rRNA  DOID:9643  Babesiosis  7.244  3.622  https://diseases.jensenlab.org/Entity?...`

Notes:
- Current files include stable disease IDs such as `DOID:*`, and also other systems such as `ICD10:*`.
- Target identifiers are mostly ENSP, but some text-mining rows contain non-protein-like identifiers.
- We will preserve these source IDs and let the existing resolvers determine what maps cleanly.

## Legacy Comparison

Old TCRD used `load-JensenLab-DISEASES.py` and loaded:
- `human_disease_knowledge_filtered.tsv`
- `human_disease_experiments_filtered.tsv`
- `human_disease_textmining_filtered.tsv`

Historical `pharos319` `disease.dtype` values include:
- `JensenLab Knowledge UniProtKB-KW`
- `JensenLab Knowledge MedlinePlus`
- `JensenLab Knowledge AmyCo`
- `JensenLab Experiment TIGA`
- `JensenLab Text Mining`

The implemented adapter preserves this `dtype` fragmentation by storing the old-style `dtype` string in each detail `source` field.

## Implemented Mapping

- Emit `Disease` nodes keyed by the source disease ID.
- Emit `ProteinDiseaseEdge` edges keyed by source protein ID to source disease ID.
- Use one detail object per source line.
- Preserve disease IDs as provided.
- Preserve Jensen channel fragmentation in `detail.source`.
- Preserve Jensen evidence payload on the detail object:
  - `evidence_terms`
  - `confidence`
  - `zscore`
  - `url`

Detail `source` values:
- knowledge: `JensenLab Knowledge {source_label}`
- experiments: `JensenLab Experiment {source_label}`
- text mining: `JensenLab Text Mining`

## Filtering

Implemented filtering:
- skip rows missing protein ID, disease ID, or disease name
- do not pre-filter non-ENSP target IDs
- rely on the target resolver to prune unresolved targets
- do not remap disease IDs in the adapter
- apply a configurable text-mining z-score cutoff when desired
  - current working and promoted full-config setting: `textmining_min_zscore: 6.0`
  - rationale: this closely matches the historical `pharos319` Jensen text-mining score floor and row volume

Observed text-mining score drift:
- `pharos319` Jensen text mining had `min(zscore) ~= 6.0`, `avg(zscore) ~= 6.66`, and `7,404` rows
- the current Jensen `human_disease_textmining_filtered.tsv` extends down to `zscore = 3.0`
- without a local cutoff, the modern ingest produced `258,220` text-mining rows
- with `zscore >= 6`, the modern score distribution aligns much more closely with legacy Pharos behavior

## Implemented Wiring

- Add a new adapter under `src/input_adapters/jensenlab/`.
- Pin an exact JensenLab source snapshot in the Registry.
- Wire the adapter into `src/use_cases/working.yaml`.
- Promote the validated adapter into:
  - `src/use_cases/pharos/pharos.yaml`
  - `src/use_cases/pharos/target_graph.yaml`

## pharos400 Translation Expectations

When `ProteinDiseaseEdge` details are converted into `pharos400.disease` rows:
- `detail.source` should become `disease.dtype`
- `detail.evidence_terms` should become `disease.evidence`
- `detail.confidence` should become `disease.conf`
- `detail.zscore` should become `disease.zscore`
- `detail.url` may be carried into `disease.reference`
- `detail.source_id` should populate `disease.did`
- resolved disease IDs should continue to key the association identity and `ncats_d2da` links
- `disease.mondoid` should remain a best-effort FK-backed resolved MONDO mapping

This keeps the modern graph ingest aligned with the useful Jensen columns that were populated in `pharos319`.

## Validation Outcomes

Working graph checks confirmed:
- `Disease` nodes and `ProteinDiseaseEdge` edges are present in `test_pharos`
- expected detail sources are present:
  - `JensenLab Knowledge UniProtKB-KW`
  - `JensenLab Knowledge MedlinePlus`
  - `JensenLab Knowledge AmyCo`
  - `JensenLab Experiment TIGA`
  - `JensenLab Text Mining`
- text-mining detail payload carries `zscore` and `url`
- experiment and knowledge details carry `confidence` and `evidence_terms`

Working MySQL checks confirmed:
- `disease.dtype` fragmentation matches the intended Jensen source breakdown
- `disease.did` now preserves source IDs
- `disease.conf` populates for knowledge, experiment, and text-mining rows
- `disease.zscore` and `disease.reference` populate for text-mining rows
- `disease.mondoid` populates for most rows when the resolved MONDO ID exists in the loaded ontology tables

Remaining interpretation:
- count differences vs `pharos319` are now primarily source-version drift, not converter loss
- the main intentional behavior choice is the `textmining_min_zscore: 6.0` cutoff to stay closer to historical Pharos text-mining scope
