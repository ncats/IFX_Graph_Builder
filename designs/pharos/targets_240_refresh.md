# TARGETS 2.4.0 Registry refresh

## Published input

`ifx_harmonizers:targets:2.4.0` was published on 2026-10-08 from Harmonizers
revision `a6d596d`. The ODIN Registry client materialized and checksum-verified
the release manifest, three ID files, `uniprot_mapping.csv`, and
`protein_identifier_links.tsv`. The published validation reports source health
and graph coherence passed. Its 18 members include 17 data files plus
`release_manifest.json`.

| Member | Rows | Size |
| --- | ---: | ---: |
| `gene_ids.tsv` | 228,209 | 101,193,554 bytes |
| `transcript_ids.tsv` | 858,959 | 372,228,914 bytes |
| `protein_ids.tsv` | 365,985 | 285,950,626 bytes |
| `uniprot_mapping.csv` | 147,520 | 86,478,765 bytes |
| `protein_identifier_links.tsv` | 1,180,651 | 119,795,574 bytes |

The prior release had 228,373 genes, 858,959 transcripts, and 232,837
proteins. All 858,959 transcript rows have `parent_ifx_gene_id`; 858,497 are
high confidence and 462 medium. All those parents exist in `gene_ids.tsv`.
Protein rows include 232,837 with a UniProt accession and 133,148 Ensembl-only
rows without one. Protein parents are present on 302,397 rows (54,590 high,
247,807 medium) and absent on 63,588. All populated parents exist in the gene
file. `uniprot_NCBI_id` has pipe-separated values on 163 rows; it cannot be
used as a single gene edge endpoint. The protein file has 20,739 canonical,
146,861 isoform, 65,237 alternate-product, and 133,148 Ensembl-only rows.

`uniprot_mapping.csv` has 140 columns. Of the ID families ODIN currently
resolves, `uniprot_ccds_id` is pipe-separated on 6,976 rows and
`uniprot_xref_ChEMBL` on 38. The link table has 1,010,445 `primary`, 126,742
`secondary`, 6,379 `shared`, and 37,085 `not_retained` rows. The last status
was absent from the preliminary design-chat description. All its protein IDs
exist in `protein_ids.tsv`. The active consolidated identifier fields in the
protein file remain the resolver input; the link table preserves weaker and
discarded claims for investigation, not automatic equivalence.

## Mapping decisions

- All related ID files and the UniProt mapping now use one explicit derived
  Registry snapshot reference. The protein resolvers read
  `uniprot_mapping.csv` from the same dataset that supplies `protein_ids.tsv`;
  no second mapping snapshot argument is accepted. Protein adapters read only
  the protein file because they do not use mapping IDs. Source adapters such
  as UniProt remain independently configured.
- Gene-transcript and gene-protein edges use the published IFX parent ID.
  `parent_resolution_source` and `parent_confidence` are retained as edge
  details. Arango stores each as a `details` entry with
  `resolution_source` and `confidence` keys. Blank protein parents produce no
  gene edge. Existing pre-2.4.0 local input files retain their older parser
  path.
- The protein resolver reads active consolidated IDs and individual mapping
  values. The TCRD resolver excludes Ensembl-only protein rows because the
  Pharos protein set is UniProt scoped. The broader Target Graph retains them.
  Pharos protein-node configurations retain `canonical_only`; all 20,739
  published canonical rows have a UniProt accession. The adapter has no
  `uniprot_only` mode, so it remains free to emit Ensembl-only proteins in
  Target Graph builds.
  It expands gene aliases only from the published protein parent, not from
  every transcript cross-reference. On 241 isoform rows the published gene
  parent differs from the canonical protein's parent, so collapse does not
  transfer an isoform's gene assertion to its canonical protein.
- On the published 2.4.0 file, Pharos's collapse retains 20,739 canonical
  protein IDs and maps 146,861 isoform rows to those IDs. The 65,237
  alternate-product rows have no canonical IFX parent and do not collapse;
  the 133,148 Ensembl-only rows are outside Pharos's UniProt scope. Among
  retained rows, no UniProt accession maps to more than one canonical IFX
  protein.
- No confidence cutoff is applied. All populated parent links in 2.4.0 are
  high or medium, so a high-plus-medium filter would remove no parent edges.
  The labels are preserved on Target Graph edges for later impact analysis.
- `protein_identifier_links.tsv` is not ingested as graph evidence. It is a
  harmonizer identity-claim artifact; secondary, shared, and not-retained
  claims must not silently become resolver equivalences or primary-source
  evidence edges.

## Build handoff

Start with `src/use_cases/pharos/impatient_target_graph.yaml`, then inspect
Gene, Transcript, Protein, GeneTranscriptEdge, GeneProteinEdge, and
TranscriptProteinEdge counts and endpoint integrity. Compare the graph's
protein scope and confidence distribution with the figures above. Continue
with `impatient_pharos`, `target_graph`, `pharos`, and the MySQL conversion only
after that narrow graph is clean. The user runs Snakemake and ETL builds.

From the repository root, the build entrypoints are:

```bash
.venv/bin/python -m src.use_cases.pharos.build_impatient_target_graph
.venv/bin/python -m src.use_cases.pharos.build_impatient_pharos
.venv/bin/python -m src.use_cases.pharos.build_target_graph
.venv/bin/python -m src.use_cases.pharos.build_pharos
.venv/bin/python -m src.use_cases.pharos.build_tcrd
```

These commands confirm before truncating. Use `--yes` only when the chosen
database is ready to be replaced. `build_tcrd` targets `pharos400`.

The local scratch `impatient_target_graph_local.yaml` and its resolver cache
were not promoted or repinned.

## Registry YAML audit (2026-10-09)

Run `.venv/bin/python -m src.use_cases.audit_registry_freshness <yaml>` before
each build. The audit now prints every source provider shared by a direct YAML
pin and any transitive derived input. It compares versions only for the same
dataset. When the derived release has no pin for that dataset but uses other
datasets from the same source, it marks the case for review and lists those
inputs so dataset naming and version alignment can be assessed. An exact-dataset
pin mismatch or unavailable derived lineage causes a nonzero exit. Review cases
alone do not fail the audit.

The published Targets 2.4.0 direct source inputs are registered and current
under the available automated checks. The full `target_graph.yaml` audit
reports that the derived release has no `uniprot:human` pin corresponding to
the graph's `uniprot:human:2026_02`; its other UniProt inputs carry three
`2026_03` pins and one `2026_03-export1` pin. The derived release likewise has
no `ncbi:gene_summary` or `ncbi:publications` pin corresponding to the graph's
pins; its other NCBI datasets are reported separately. These are not
exact-dataset mismatches.

The existing freshness audit also reports roughly two dozen known updates and one dependent
SureChEMBL rebuild for `target_graph.yaml`; `pharos.yaml` reports 27 known
updates and one dependent rebuild. The narrow `impatient_target_graph.yaml`
has two unregistered derived drug graph pins, and its post-processing YAML has
an unregistered `target_graph:tdl_updates:2026-09-04` pin. These are hard
build blockers for those exact Registry references. The other post-processing
YAMLs have no known automatic update but retain manual freshness checks for
the Targets 2.3.0 baseline and TDL updates. No ETL was run.
