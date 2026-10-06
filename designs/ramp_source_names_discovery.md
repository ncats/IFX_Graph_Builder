# RaMP `source.commonName` discovery (2026-10-06)

The source table records IDs each input uses to report retained RaMP data. Its
`commonName` is provider-attributed input to the legacy most-common-name
selection; a name borrowed from another provider would count as false support.
Do not fill missing provider names from the harmonized metabolite or another
provider. `priorityHMDBStatus` is a metabolite-group attribute repeated on
source rows, not a status claimed by PubChem or PFOCR.

## Observations from the pinned stage and inputs

- The new diagnostic SQLite has 1,416 WikiPathways-KEGG source IDs and 1,360
  RaMP IDs, versus 268 IDs and 262 RaMP IDs in unreleased 3.0.12. All 268 old
  IDs remain; 1,148 IDs are added, 1,147 of them unnamed. In the graph, 1,518
  KEGG compound nodes carry WikiPathways provenance, but only 269 carry a
  WikiPathways name. The new equivalence adapter registers KEGG cross-references
  from human WikiPathways RDF as WikiPathways source IDs, whereas the old source
  table mostly reflected the KEGG IDs with pathway context. This is a coverage
  change in ID evidence, not a new KEGG data source.
- HMDB XML gives each metabolite record a `name` and cross-reference IDs such as
  CAS and KEGG. The adapter labels the primary HMDB node but emits unnamed
  cross-reference stubs. The export consequently has names on 217,920 of
  462,812 HMDB rows and none on 5,900 HMDB-KEGG rows. The XML record's own name
  can be attached to its cross-reference nodes with HMDB provenance; this is
  within-record attribution, not a cross-provider fallback.
- RefMet's pinned CSV has 208,170 populated `refmet_name` cells. Its adapter
  labels the primary RefMet ID but leaves its five cross-reference families
  unnamed. The export has names on 15,933 of 67,186 RefMet rows.
- The pinned WikiPathways RDF labels primary metabolite subjects, while some
  bdb cross-references become unnamed nodes. Carrying the subject label to
  cross-references would be source-native within-record attribution. It may be
  ambiguous when multiple labeled subjects report the same cross-reference;
  retain all labels as evidence and use a deterministic selection rule in the
  export.
- Every one of 344,152 graph `RheaMetaboliteReactionEdge` records has a `name`
  from the Rhea compound RDF. The exporter currently ignores it, leaving all
  14,986 Rhea source rows unnamed. This can be fixed in the exporter without a
  graph rebuild.
- The pinned LipidMaps SDF has 50,572 entries: `NAME` on 45,528,
  `SYSTEMATIC_NAME` on 46,257, `ABBREVIATION` on 34,645, and `SYNONYMS` on
  26,103; `COMMON_NAME` is absent. For `LMFA00000001`, `NAME` is
  `2-methoxy-12-methyloctadec-17-en-5-ynoyl anhydride`, `ABBREVIATION` is
  `FA 40:7;O3`, and `SYNONYMS` contains `Acetylenic acids`. The old loader
  collected all names and overwrote the per-ID common name as it read them,
  explaining why the old row used the broad last synonym. The new exporter
  selects the long source `NAME`; a better display-name policy needs an explicit
  preference among source fields, without claiming the synonym is more specific.
- PubChem's pinned `cid_molecular_info.tsv` has 172,451 rows and 172,316
  populated `iupac_name` cells, but no common name or compound title field. The
  current graph faithfully retains IUPAC in chemistry, not in names. Calling
  it a common name would mislabel the value. PubChem's official PUG REST
  compound-property API documents `Title` as the compound summary-page title:
  https://pubchem.ncbi.nlm.nih.gov/pcfe/docs/markdown/pug-rest.md . A separately
  pinned CID-to-Title snapshot is the source-native path to common names; the
  current molecular-info snapshot cannot provide them.
  A bounded live response for CIDs 1 and 2244 had the documented shape
  `PropertyTable.Properties[{CID, Title}]` and titles
  `Acetyl-DL-carnitine` and `Aspirin`, respectively.
  A later bounded live check found that CID 25 is a valid compound but its
  Title response is `PropertyTable.Properties[{CID: 25}]`, with the `Title`
  field absent. The Registry recipe accepts this complete CID row, preserves
  the compound, and projects a blank title; it still rejects a response that
  omits the requested CID row entirely.
- Reactome's pinned `ChEBI2Reactome_All_Levels.txt` has ChEBI ID, pathway ID,
  URL, pathway name, evidence code, and species, but no metabolite name. The
  old Reactome parser looked up names via ChEBI; its source attribution was not
  Reactome-native. The unreleased 3.0.12 SQLite already has no Reactome names.
  PFOCR likewise provides association IDs rather than metabolite labels. Leave
  these provider names empty unless a provider-owned name field is found.
- There is no KEGG compound-name input in the current RaMP graph. `hmdb_kegg`
  and `wikipathways_kegg` are KEGG IDs supplied by HMDB and WikiPathways;
  names on those rows should be attributed to those reporting providers.

## Implemented name policy

1. Preserve the source-row no-cross-provider-fallback rule. The exporter uses
   Rhea edge `name` directly when registering its source ID.
2. HMDB, RefMet, and WikiPathways mapping adapters now attach the reporting
   record's name to its cross-reference ID stubs with the same provider and
   source field. Repeated assertions for the same ID/name are deduplicated;
   differing source-owned labels remain available for graph merging. The user
   will rebuild the graph and validated harmonization stage, then export the
   two-table SQLite and compare field coverage.
3. LipidMaps source rows prefer `ABBREVIATION`, then `NAME`, then
   `SYSTEMATIC_NAME`; its chemical `common_name` remains available at the same
   priority as `NAME`. The old "last synonym wins" behavior is not reproduced.
   Reactome and PFOCR common names stay empty without provider-owned labels.
4. IFX_Registry's new `compound_records` recipe pins exact CID-to-Title API
   results beside the compound records; `cid_molecular_info` projects a `title`
   column only from that pinned evidence. The ODIN PubChem adapter uses `title`
   as a PubChem name when present, but does not treat IUPAC as a common name.
   RaMP now pins `pubchem:cid_molecular_info:deps-74c850f19f6c`, derived from
   `pubchem:compound_records:deps-076242a4ae17`. Registry validation reports
   172,451 molecular-info rows and 171,776 populated Titles. A new graph build
   and harmonization stage are required before these names appear in SQLite.
   Direct TSV profiling confirmed the `title` column and those counts. Examples
   include CID 1 `Acetyl-DL-carnitine`, CID 962 `Water`, and CID 2244 `Aspirin`;
   CID 25 is present with a blank Title. A bounded adapter read emits the
   populated values as PubChem-attributed `MetaboliteName` records.

The user confirmed the source-owned propagation policy and LipidMaps ranking
after discovery. The user runs graph/ETL rebuilds.
