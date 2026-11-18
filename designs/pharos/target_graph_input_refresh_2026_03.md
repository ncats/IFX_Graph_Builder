# Target Graph Input Refresh, March 2026

## Goal
Integrate the refreshed target graph identifier exports shared on March 24, 2026 while preserving the prior manual inputs for reference and limiting the first-pass code changes to `working.yaml`.

## Source Files

New package received from:
- `/Users/kelleherkj/Downloads/260315Targets (1)`

Files supplied:
- `gene_ids.tsv`
- `protein_ids.tsv`
- `transcript_ids.tsv`
- `uniprotkb_mapping.csv`

Current active manual inputs before refresh:
- `input_files/manual/target_graph/gene_ids.csv`
- `input_files/manual/target_graph/protein_ids.csv`
- `input_files/manual/target_graph/transcript_ids.csv`
- `input_files/manual/target_graph/uniprotkb_mapping.csv`

Subsequent replacement drops received during March 2026 evaluation:
- `/Users/kelleherkj/Downloads/gene_ids.tsv`
- `/Users/kelleherkj/Downloads/protein_ids.tsv`
- `/Users/kelleherkj/Downloads/transcript_ids.tsv`
- `/Users/kelleherkj/Downloads/260325Targets/`

## Discovery Summary

### File format changes
- `gene_ids`, `protein_ids`, and `transcript_ids` are now TSV rather than CSV.
- `uniprotkb_mapping.csv` remains comma-delimited.
- Existing shared parsing in `src/shared/csv_parser.py` assumed comma-delimited input, so the new TSV exports are not drop-in compatible without parser changes.
- Parsing now keys off the file extension only: `.tsv` uses tab delimiters and other flat files continue to use commas.

### Schema comparison

#### Gene
- Column count stayed at 28.
- Column set is unchanged.
- Differences are limited to column ordering and refreshed row content.

#### Protein
- Column count changed from 43 to 36.
- New file drops:
  - `uniprot_secondaryAccessions`
  - `uniprot_uniProtkbId`
  - `uniprot_sequence`
  - `uniprot_references`
  - `uniprot_feature_Domain`
  - `uniprot_feature_Region`
  - `uniprot_feature_Coiled_coil`
  - `uniprot_FUNCTION`
- New file adds:
  - `ensembl_uniprot_transcript_mismatch`

#### Transcript
- Column count changed from 21 to 17.
- New file drops:
  - `ensembl_transcript_name`
  - `ensembl_trans_bp_start`
  - `ensembl_trans_bp_end`
  - `ensembl_refseq_NM`

### Data coverage impact

#### Accepted protein-field loss
- The dropped target graph protein annotations were populated in the old file.
- For rebuild purposes, losing target graph copies of sequence and function is acceptable because richer protein annotations already come from UniProt in the Pharos build.
- The target graph protein file still provides the ID linkage needed by the target graph adapters and resolvers.

#### Remaining transcript impact
- Transcript nodes will no longer receive genomic start/end coordinates from target graph.
- Transcript equivalent IDs will no longer include the `ensembl_refseq_NM` source.
- `ensembl_transcript_name` aliases will no longer be emitted from this source.
- Remaining transcript identifiers and core fields are still present, so transcript linkage remains usable.

### Content refresh notes
- External identifier overlap is high but not exact.
- Internal IFX IDs are regenerated in the new exports for matching external records, so downstream comparisons should be anchored on external identifiers rather than previous IFX IDs.

## Subsequent Drops

### March 24 top-level Downloads drop

This follow-up set improved the transcript side but was not adopted as the active package.

Improvements:
- `transcript_ids.tsv` used only `IFXTranscript` primary IDs.
- `Total_Mapping_Ratio` in `gene_ids.tsv` was no longer pipe-delimited.

Remaining issues:
- `gene_ids.tsv` still carried multi-valued `consolidated_NCBI_id` and `ensembl_strand`.
- `transcript_ids.tsv` timestamps were consistently ISO seconds-only.
- `protein_ids.tsv` had a serious reviewed-protein duplication problem:
  - `42,562` reviewed rows but still only `20,431` unique reviewed UniProt IDs
  - `10,587` reviewed UniProt accessions appeared under multiple distinct `IFXProtein` IDs
  - unreviewed rows did not show the same duplication pattern

Because that duplication changed the effective reviewed protein node count without adding new unique reviewed proteins, the top-level Downloads protein file was rejected.

### March 25 `260325Targets`

This later drop became the preferred active set.

Observed simplifications and improvements:
- `transcript_ids.tsv`
  - all primary IDs are `IFXTranscript`
  - dates are consistent seconds-only ISO timestamps
- `protein_ids.tsv`
  - `237,969` total rows
  - `20,431` reviewed rows
  - no reviewed UniProt duplication
  - restores richer UniProt-derived columns such as:
    - `uniprot_secondaryAccessions`
    - `uniprot_uniProtkbId`
    - `uniprot_sequence`
    - `uniprot_references`
    - `uniprot_FUNCTION`
  - also retains `ensembl_uniprot_transcript_mismatch`
- `gene_ids.tsv`
  - `Total_Mapping_Ratio` remains singleton-valued
  - `consolidated_NCBI_id` and `ensembl_strand` are still multi-valued in a minority of rows and still require defensive parser logic

## Current Active Set

The currently preferred target graph manual inputs are:
- `input_files/manual/target_graph/gene_ids.tsv`
- `input_files/manual/target_graph/protein_ids.tsv`
- `input_files/manual/target_graph/transcript_ids.tsv`
- `input_files/manual/target_graph/uniprotkb_mapping_20260507.csv`

The active TSVs came from:
- `/Users/kelleherkj/Downloads/260325Targets/`

Archived intermediate sets were preserved under:
- `input_files/manual/target_graph/archive/2026_03_24_jess_refresh/`
- `input_files/manual/target_graph/archive/2026_03_26_260325Targets/`

## Implementation Plan

1. Update shared flat-file parsing so `.tsv` files are read as tab-delimited.
2. Keep the target graph parser behavior tolerant of missing optional fields.
3. Preserve the pre-refresh manual files under an archive/reference location before swapping active files.
4. Copy the refreshed files into `input_files/manual/target_graph/` without changing the currently referenced Pharos YAML inputs yet.
5. Use the new TSV files plus a date-stamped UniProt mapping file from `src/use_cases/working.yaml` first.
6. Leave promotion into `src/use_cases/pharos/target_graph.yaml` and `src/use_cases/pharos/pharos.yaml` until validation is complete.

## File Preservation Plan

Preserve the prior manual inputs under a dated reference directory:
- `input_files/manual/target_graph/archive/2025_05_snapshot/`

Files to preserve:
- `gene_ids.csv`
- `protein_ids.csv`
- `transcript_ids.csv`
- `uniprotkb_mapping.csv`

This keeps the prior package available for schema checks, regression comparisons, and resolver debugging.

For the first-pass working config:
- keep the prior CSV inputs in place for existing YAMLs
- copy in `gene_ids.tsv`, `protein_ids.tsv`, and `transcript_ids.tsv` alongside them
- copy the refreshed mapping file as `uniprotkb_mapping_20260507.csv`

## Validation Notes

- First-pass validation should use `src/use_cases/working.yaml`.
- The user will acquire required source snapshots and run ETL steps.
- After the user runs the build, compare target graph node and edge counts plus resolver coverage before promoting config changes into the Pharos YAMLs.

## May 2026 Protein Refresh

The May 2026 target graph protein refresh changed the Pharos-selection contract.

Verified ingest rules:
- `pharos.yaml` should keep only rows where `is_canonical == True`
- `target_graph.yaml` should ingest all protein rows
- `canonical_isoform_status` should be carried onto `Protein` nodes
- `current_tdls.csv` should expose `canonical_isoform_status`
- `uniprot_secondaryAccessions` remains part of `protein_ids.tsv` and should continue to feed resolver aliases

As a result:
- the older `collapse_reviewed_targets` heuristic is retired
- Pharos protein node selection is now a direct `canonical_only` adapter filter
- full target graph ingest preserves the upstream `canonical`, `isoform`, and `alternate_product` labeling without local collapse logic

### May 19 protein_ids.tsv replacement

Source file:
- `/Users/kelleherkj/Downloads/protein_ids_2026-05-19.tsv`

Active replacement:
- `input_files/manual/target_graph/protein_ids.tsv`

Archived prior active file:
- `input_files/manual/target_graph/archive/2026_05_19_pre_protein_refresh/protein_ids.tsv`

Payload profile:
- Row count stayed stable at `166,833` data rows.
- Primary `ncats_protein_id` key set stayed stable: no added IDs, removed IDs, or duplicate primary IDs.
- The refreshed file adds `uniprot_uniProtkbId` and `uniprot_sequence_length`.
- `uniprot_uniProtkbId` is already consumed by `TargetGraphProteinParser.get_equivalent_ids()`, so no parser change is needed for that restored alias field.
- `uniprot_sequence_length` is currently extra payload and is not consumed by target graph parsing.
- Drift in shared columns is mostly timestamp refresh. Resolver-relevant drift observed during profiling:
  - `92` rows changed `is_canonical`
  - `370` rows changed `canonical_ifx_id`
  - `345` rows changed `protein_grouping_reason`

Validation should start with the existing narrow target graph working path before broader Pharos rebuilds.

### Pharos canonical target alias handling

For Pharos/TCRD target resolution, `collapse_to_canonical` keeps canonical proteins as the only output target IDs but allows non-canonical rows to contribute aliases to their mapped canonical target. If a non-canonical protein row has a `canonical_ifx_id` that points to a protein row surviving the canonical collapse filter, the resolver maps that row's protein aliases plus linked transcript and gene aliases to the canonical `IFXProtein` ID. Non-canonical rows without a valid canonical target are still dropped.
