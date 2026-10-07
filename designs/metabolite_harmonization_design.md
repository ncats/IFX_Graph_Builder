# Metabolite Harmonization Graph Design

## Purpose

The `metabolite_harmonization` graph is the evidence layer for the next
generation of RaMP metabolite harmonization.

The graph should preserve what source systems say about metabolite identifiers,
names, chemistry, pathways, proteins, and ontology context. It should not decide
the final harmonized metabolite groups during ingest. Those decisions belong in
explicit harmonization rules and saved snapshots that can be inspected in the
Harmonization Workbench.

## High-Level Model

The central node is `MetaboliteIdentifier`. Examples include `HMDB:...`,
`CHEBI:...`, `PUBCHEM.COMPOUND:...`, `LIPIDMAPS:...`, `REFMET:...`,
`InChIKey:...`, and other source-provided IDs.

Important supporting nodes include:

- `PathwayIdentifier` for source pathway IDs.
- `ProteinIdentifier` and `GeneIdentifier` for source pathway/protein context.
- Rhea-specific reaction nodes and EC reaction class nodes.
- ChEBI ontology nodes such as `ChemicalEntity` and role nodes.
- Source-specific metabolite classification terms from HMDB/ClassyFire and
  LipidMaps.
- HMDB ontology terms for source-reported metabolite context such as biofluid,
  tissue, source, health condition, subcellular location, and application.
- `MetaboliteHarmonizationClique*` nodes created by saved workbench snapshots.

Important edges include:

- source-reported metabolite ID equivalence edges;
- metabolite/protein/gene pathway membership edges;
- HMDB metabolite-protein association edges;
- Rhea metabolite-reaction, protein-reaction, reaction-class, reaction-direction,
  and EC class-parent edges;
- Reactome pathway hierarchy edges;
- ChEBI ontology edges;
- ChEBI `ChemicalEntity` to `MetaboliteIdentifier` bridge edges;
- metabolite-to-most-specific-class and class-parent edges for HMDB/ClassyFire
  and LipidMaps classification context;
- HMDB ontology parent edges and metabolite-to-HMDB-ontology membership edges;
- snapshot clique membership edges.

## What We Keep

Adapters should be source-preserving by default. If a source provides an
identifier or association and we can normalize it into a stable CURIE-like ID,
we generally keep it.

Current retained evidence includes:

- metabolite identifier mappings from HMDB, WikiPathways, LipidMaps, and RefMet;
- source-attributed metabolite names and synonyms;
- chemical properties from HMDB, ChEBI, LipidMaps, and PubChem-derived files;
- pathway context from HMDB/SMPDB, WikiPathways, Reactome, and PFOCR;
- HMDB metabolite-protein associations;
- HMDB ClassyFire `super_class`, `class`, and `sub_class` taxonomy levels;
- LipidMaps `CATEGORY`, `MAIN_CLASS`, `SUB_CLASS`, and sparse `CLASS_LEVEL4`
  classification levels;
- HMDB ontology terms and metabolite memberships from the source XML ontology
  tree, limited to the known HMDB ontology category branches used by RaMP:
  `Biofluid and excreta`, `Health condition`, `Industrial application`,
  `Organ and components`, `Source`, `Subcellular`, and
  `Tissue and substructures`;
- active Rhea reactions, using Rhea's native reaction/direction/participant
  model;
- Rhea participant IDs as `MetaboliteIdentifier` nodes, including ChEBI-backed
  participants, Rhea generic compounds, and Rhea polymer compounds;
- human-filtered Rhea UniProt protein mappings;
- ExPASy enzyme classes used as Rhea reaction class metadata;
- Reactome pathway hierarchy;
- the ChEBI full ontology graph;
- bridge edges linking matching ChEBI `ChemicalEntity` and
  `MetaboliteIdentifier` nodes;
- saved harmonization snapshots and their clique memberships.

Names, synonyms, and chemical properties are stored on the relevant identifier
node as source-attributed nested data. We do not copy a name from one identifier
to another just because a source reports an equivalence.

Biological associations are likewise stored only on the identifier the source
used for that assertion. An adapter must not copy a pathway, reaction,
classification, ontology, or protein association onto every equivalent
identifier. The equivalence edges already connect those identifiers, and
downstream clique materialization can traverse or project the association when
needed. Association evidence uses `source_id` to preserve the source-native
metabolite identifier even if a downstream representation later rewires the
edge to a harmonized metabolite.

The adapter audit found that HMDB, Reactome, PFOCR, Rhea, and LipidMaps already
follow this invariant. WikiPathways pathway context was the exception: it
expanded each RDF subject through all of its BridgeDb `bdb*` fields. The
pathway adapter now emits metabolite, gene, and protein pathway edges only for
the normalized RDF subject. The separate WikiPathways metabolite-equivalence
adapter continues to emit the BridgeDb identifiers and mapping edges.

InChIKey can appear in two ways:

- as source-reported identifier evidence, for example `InChIKey:...`;
- as chemical-property evidence inside `chem_props`.

Rules must be explicit about which InChIKey evidence they use.

Chemical-property records preserve source-reported `inchi_key` and
`inchi_key_prefix` fields unchanged. During chemistry ingest, a standard
InChIKey is also calculated from the first available SMILES field in this
order: `iso_smiles`, `isomeric_smiles`, then `canonical_smiles`. The calculated
value is stored separately as `derived_inchi_key` and
`derived_inchi_key_prefix`, together with its input field, calculation method,
toolkit version, and any calculation error. A derived key never overwrites or
masquerades as a source-reported key.

LipidMaps records may carry a molecular drawing in the SDF mol block while
omitting the `SMILES`, `INCHI`, `INCHI_KEY`, and formula tags. When a reported
SMILES tag is absent, the LipidMaps chemistry adapter converts a nonempty mol
block to canonical isomeric SMILES and stores it in `canonical_smiles`. This is
normal adapter parsing of the versioned LipidMaps payload, so its datasource
provenance is supplied by the framework. Fully specified mol blocks also
produce a `derived_inchi_key` with `sdf_mol_block` recorded as its input.
Generic/query atoms such as `R`, `R#`, and `R1` are preserved as SMILES
wildcards, make the generic-structure classifier return true, and intentionally
do not produce an InChIKey. Zero-atom mol blocks remain unclassified.

The current ChEBI three-star SDF uses uppercase property tags including
`INCHIKEY`, `INCHI`, `MONOISOTOPIC_MASS`, `ChEBI NAME`, and `FORMULA`. Older
three-star payloads used the mixed-case tags `InChIKey`, `InChI`,
`Monoisotopic Mass`, `ChEBI Name`, and `Formulae`. The ChEBI chemistry adapter
accepts both explicit tag sets, preferring the current tags. This compatibility
is intentional: a case-sensitive legacy parser silently retained only ChEBI
ID, SMILES, and molecular weight after the upstream tag change, which removed
reported keys and other chemistry fields without failing the build.
The SDF's ChEBI name is also stored as a source-attributed name on the same
`MetaboliteIdentifier`, with the actual SDF tag as `source_field`. It remains
in `chem_props.common_name` for chemistry consumers. This lets downstream
RaMP exports include the ChEBI-reported name in `analytesynonym` without
assigning a name from another identifier in the metabolite group.

The reported-key harmonization rules continue to consume only `inchi_key` and
`inchi_key_prefix`. A separate optional derived-key rule can include calculated
keys in a pipeline. Downstream compatibility exports such as RaMP SQLite should
continue to export source-reported keys unless the export contract is
explicitly expanded to include derived chemistry.

Rhea participant names, HTML names, and formula strings are stored as
reaction-participant edge context. They are not stored as `chem_props`, because
Rhea does not provide structure-derived SMILES, InChI, or InChIKey values in the
reaction bundle.

HMDB and LipidMaps classifications are source-specific evidence, not canonical
chemical ontology reconciliation. Term IDs therefore use IFX source-stable
prefixes such as `HMDB.CLASS:...`, `LipidMaps.CLASS:...`, and
`HMDB.ONTOLOGY:...` rather than pretending these labels are ChEBI, UBERON,
MONDO, or other canonical ontology identifiers.

Metabolite membership edges link only to the most specific source-provided class
or HMDB ontology child term. Broader class and ontology membership should be
derived by traversing parent edges, rather than storing redundant direct
metabolite edges at every level.

## What We Skip

The graph is broad, but it is not “everything from everywhere.”

Current intentional skips:

- non-human WikiPathways RDF pathway files;
- non-human PFOCR pathway files;
- pathway sources that are not part of this RaMP-focused graph yet, such as the
  current Pathway Commons usage;
- cross-identifier reconciliation inside adapters;
- final metabolite granularity decisions during ingest;
- implicit RaMP legacy prefix allowlists or denylists inside adapters.
- Rhea metabolite equivalence edges. Rhea generic and polymer accessions become
  `MetaboliteIdentifier` nodes and participate in reactions, but Rhea does not
  make them equivalent to ChEBI or other metabolite IDs.
- HMDB ontology branches outside the RaMP ontology categories. The HMDB XML
  contains other nested context such as pathway labels; these are not emitted as
  `HmdbOntologyTerm` nodes because they duplicate pathway-specific graph
  adapters and can otherwise become misleading `ontology_type` values.

RefMet-only orphan handling is not an ingest concern. The harmonization
workbench handles it through the broader secondary-source-only cleanup rule,
which also covers PubChem-, ChEBI-, and LipidMaps-only results without meaningful
biological connections.

RaMP's legacy SQLite build applies additional ontology export logic after HMDB
parsing. In `ramp-backend-ncats`, `hmdbData.py` extracts only specific HMDB
ontology branches and `EntityBuilder.recordOntology()` applies
`config/options.yml` category-specific denylist terms before writing
`ontology` and `analytehasontology`. The graph ingest intentionally does not
apply that full SQLite denylist yet; the RaMP-compatible SQLite export should
apply the denylist and any exportable-term rules when flattening graph evidence
into legacy tables.

Rhea `has_human_protein`, `only_human_metabolites`, and `is_cofactor` style
flags are treated as derived annotations rather than source fields for the
first-pass ingest.

## Harmonization Rules

The graph stores evidence. Harmonization rules decide how to group IDs.

The workbench can create saved snapshots from rule pipelines such as:

- use all source-reported equivalence edges;
- apply the RaMP mapping deny list;
- ignore generic/non-generic structure equivalences;
- merge by InChIKey prefix;
- merge by InChIKey first two components;
- switch between InChIKey strategies by molecular-weight cutoff.

Rules should be independently configurable, ordered, and inspectable. A saved
snapshot should record which rules ran, their order, parameters, and summary
statistics.

### Molecular-Weight Cutoff Policy

The molecular-weight-sensitive InChIKey rule uses chemical properties from the
current harmonization graph. With the default cutoff of 500, it uses the
stereochemistry-sensitive InChIKey duplex below the cutoff and the
connectivity-only InChIKey prefix at or above the cutoff.

This policy intentionally does not reproduce decisions caused by missing
molecular-weight values in a legacy RaMP SQLite build. When the current graph
has a supported molecular weight that legacy SQLite lacked, the graph value
controls the strategy. Resulting merges are expected harmonization outcomes,
not compatibility defects merely because SQLite kept the identifiers separate.

The WikiPathways prefix rule may also ignore configured BridgeDb source fields.
If an identifier is introduced only as a WikiPathways xref through ignored
fields, the rule removes its WikiPathways node support as well as its mapping
detail. It remains active when it is also a WikiPathways source subject, is
introduced through an allowed WikiPathways field, or has support from another
datasource.

## Known Granularity Problems

The graph intentionally does not solve metabolite granularity at ingest time.
These are rule and curation questions:

- generic structures with R groups;
- generic IDs linked to specific structures;
- salts and parent compounds;
- open-chain and cyclic carbohydrate forms;
- stereochemistry that matters in some contexts and not others;
- lipids where analytical and biological granularity differ;
- source cliques that share InChIKey evidence but should not always merge.

The goal is to keep enough evidence for experts to inspect these cases, compare
rules, and decide which behavior should become RaMP release behavior.

## Current Configuration

The promoted RaMP metabolite harmonization graph lives in:

- `src/use_cases/ramp/ramp.yaml`
- `src/use_cases/ramp/build_ramp.py`

`src/use_cases/working.yaml` remains a scratch configuration for experiments
before changes are promoted.

For structural model changes, prefer a full truncate/rebuild. Use `--resume`
only when the graph shape has not changed.

## Deferred Decisions

- How to export a backward-compatible RaMP SQLite database from this graph.
- Which rule profile becomes the released RaMP harmonization policy.
- How expert curation decisions should be stored and replayed. See
  `designs/metabolite_harmonization_curation_design.md`.
- Whether additional class systems from RefMet or ChEBI should become
  first-class graph context.
- Whether broader pathway sources belong in this graph later.
