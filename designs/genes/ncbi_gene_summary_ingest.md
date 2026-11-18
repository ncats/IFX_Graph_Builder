# NCBI Gene Summary Ingest

## Source

- URL: `https://ftp.ncbi.nlm.nih.gov/gene/DATA/gene_summary.gz`
- Format: gzip-compressed TSV
- Columns observed on 2026-05-20:
  - `#tax_id`
  - `GeneID`
  - `Source`
  - `Summary`

## Payload Profile

Profiled from a live download to `/private/tmp/ncbi_gene_summary.gz`.

- Total rows: 2,949,218
- Human rows (`#tax_id == 9606`): 151,658
- Distinct human `GeneID` values: 151,658
- Duplicate human `GeneID` values: 0
- Empty human summaries: 0
- Human summary length:
  - minimum: 69 characters
  - median: 383 characters
  - maximum: 3,446 characters
- Human source breakdown:
  - `RefSeq`: 143,352
  - `Alliance of Genome Resources`: 7,478
  - `OMIM`: 821
  - `other`: 7

Example human row:

```text
9606	1	RefSeq	The protein encoded by this gene is a plasma glycoprotein of unknown function. ...
```

## Identifier Coverage

The source provides NCBI Gene IDs only. The adapter emits `Gene` nodes with source IDs shaped as `NCBIGene:<GeneID>`.

Coverage against current target resolver inputs:

- Human summary GeneIDs: 151,658
- Target graph GeneIDs with NCBI IDs: 193,802
- Human summary IDs matching target graph GeneIDs: 151,657
- Pharos canonical protein file NCBI Gene IDs: 19,295
- Human summary IDs matching Pharos canonical protein file GeneIDs: 18,710

Expected target graph landing scope is the matched `Gene` subset. Expected Pharos landing scope is the subset of emitted `Gene` summaries that the existing `tcrd_targets` resolver can retype onto canonical `Protein` nodes.

## Provisional Mapping

First-pass scope:

- Filter to human rows only (`#tax_id == 9606`).
- Emit `Gene(id="NCBIGene:<GeneID>")`.
- Store the source text as a new source-specific field, `ncbi_gene_summary`, on both `Gene` and `Protein`.
- In target graph builds, the value can remain on `Gene`; in Pharos builds, the existing `tcrd_targets` resolver can retype compatible `Gene` payloads onto canonical `Protein`.
- Preserve `Source` as a companion field only if needed for downstream display or auditing; otherwise rely on normal datasource provenance.
- Do not overwrite `description` or `uniprot_function`, which currently contain UniProt function text.

## Version Strategy

No discrete release number was identified during discovery. For a future refresh:

- Acquire and pin an exact NCBI publication snapshot in the Registry.
- Capture `Last-Modified` as `version_date`.
- Capture run date as `download_date`.
- Preserve the version and dates in Registry snapshot metadata for the adapter.

## Downstream MySQL Mapping

Legacy `pharos319.tdl_info` stores NCBI gene summaries as:

- `itype`: `NCBI Gene Summary`
- `protein_id`: populated
- `string_value`: summary text
- numeric/date/boolean value columns: null

Legacy profile on 2026-05-20:

- `NCBI Gene Summary`: 12,900 rows
- `NCBI Gene PubMed Count`: 20,153 rows

The current TCRD converter should export `Protein.ncbi_gene_summary` to `tdl_info` using the same legacy `itype` and `string_value` pattern.

Working MySQL validation after `working_mysql.yaml` build:

- `pharos400_working.tdl_info` rows with `itype = 'NCBI Gene Summary'`: 19,110
- all 19,110 rows have `protein_id`
- all 19,110 rows have `string_value`
- `target_id`, `integer_value`, and `number_value` are null for these rows
- string length range: 71 to 3,446 characters

## Open Questions

- Whether to store the source label from the TSV (`RefSeq`, `OMIM`, etc.) as a separate field such as `ncbi_gene_summary_source`.
