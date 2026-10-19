# DrugCentral release 56 dump cutover

## Source profile

The Registry source snapshot `drugcentral:drug_exports:56-export1` represents
DrugCentral database version 56 (2026-09-02), captured 2026-10-05. The
published export contract declares UTF-8 TSV with CSV quoting and empty fields
for SQL `NULL`. It contains eight members, including the source SQL and
contract. The three ODIN graph adapters need these members:

| File | Rows | Adapter use |
| --- | ---: | --- |
| `drugcentral_structures.smiles.tsv` | 5,116 | Drug ligand ID, name, SMILES; indication drug name |
| `drugcentral_drug_target_interaction.tsv.gz` | 16,926 | Human protein-ligand activity and indication protein join |
| `drugcentral_indications.tsv` | 12,231 | Drug-level disease indications |

Structure IDs and names are complete and unique; 736 SMILES values are blank.
All interaction rows are `Homo sapiens` with populated structure IDs and
accessions; 642 accession cells contain `|`-separated IDs. `ACT_VALUE` is
blank on 1,302 rows, activity PMIDs on 16,240, and MOA PMIDs on 16,563.
Numeric values that are present parse as floats or integers as expected. All
interaction and indication structure IDs exist in the structure export.

The interaction export has 16,908 distinct rows under the old SQL projection
used by `ProteinDrugEdgeAdapter`. The indication file has 2,845 distinct drug
IDs, 9,031 rows whose drug has a human interaction, and 63,524 distinct rows
under the former SQL indication join and projection. This matches the
published export contract. Indication rows are drug level; no protein ID is
present on them.

## Mapping

- The adapters use the same exact Registry source snapshot and read named
  files through `RegistryDataset.file()`. The framework supplies the snapshot's
  version and dates as provenance. PostgreSQL credentials and the version-54
  external pin are removed from the four Pharos YAML consumers.
- Structure `ID` remains the DrugCentral ligand ID. Interaction `ACCESSION`
  remains a UniProt ID for resolver handling. The protein-ligand adapter
  preserves activity details, numeric types, multi-accession expansion, SQL
  projection deduplication, and grouping by protein/ligand endpoints.
- The indication adapter joins `STRUCT_ID` to the human interaction export and
  structure name, then keeps the existing UMLS/local disease ID choice and
  protein-disease edge detail deduplication. It does not infer a different
  disease identifier from DOID or SNOMEDCT.
- Missing files, required columns, or malformed rows fail. The adapters have
  no live database fallback.

## Validation and handoff

Direct adapter parsing of the materialized Registry files emitted 5,116
DrugCentral ligand nodes, 18,228 protein-ligand edges with 18,731 activity
details, and 2,101 disease nodes plus 45,187 protein-disease edges with 66,627
indication details. These counts reflect multi-accession expansion and graph
deduplication, so they differ from raw export row counts. The focused fixture
tests cover the SQL `DISTINCT` projection, multi-accession expansion, numeric
activity fields, indication join, missing required columns, and all four YAML
pins. They passed without ETL.

The user runs ETL builds. Start with the narrow `impatient_pharos` build,
check DrugCentral ligand, protein-ligand, disease, and protein-disease counts
and sample details, then run full `target_graph` and `pharos` builds. Compare
DrugCentral-only graph counts and edge endpoint integrity with this source
profile and the existing `pharos319`/`pharos400` outputs where appropriate.
