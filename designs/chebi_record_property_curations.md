# Reusable ChEBI record-property curations

## Decision

ChEBI ontology-record corrections are published in the dedicated
`chebi_record_properties` curation stream. They are not added to
`metabolite_record_properties`.

The portable target contract is:

```json
{
  "kind": "node",
  "curation_set": "chebi",
  "model_type": "ChemicalEntity",
  "id": "CHEBI:137735"
}
```

The fixed `chebi` curation set identifies the source record rather than the
database in which it happens to be loaded. Any graph build that contains the
target `ChemicalEntity` collection can explicitly opt into this bundle.

## Representation boundary

RaMP metabolite harmonization currently exposes two distinct ChEBI-backed
records, and they remain separately curatable:

- `ChemicalEntity` is the ChEBI ontology record and uses
  `chebi_record_properties`.
- `MetaboliteIdentifier.chem_props` is the ChEBI SDF-derived metabolite record
  and continues to use `metabolite_record_properties` with the
  `metabolite_harmonization` curation set.

A correction to one representation does not silently propagate to the other.
When both source records are wrong, both are curated explicitly. The MW finding
detail links to each applicable record editor so reviewers can see which
representation supplied the evidence.

## Editable and derived fields

The ChEBI bundle allows source-property corrections to `charge`, `formula`,
`inchi`, `inchi_key`, `mass`, `monoisotopic_mass`, `smiles`, and `wurcs`.
Calculated structure fields are protected from direct editing.

After applying a ChEBI record curation, structure-derived chemistry is
recomputed from the effective curated SMILES. This replaces calculated masses,
component chemistry, method metadata, and calculation errors together, so a
portable curation cannot leave stale calculated evidence in either a
harmonization stage or an opted-in graph build.

The same replacement rule applies independently to a curated structure field
inside `MetaboliteIdentifier.chem_props`. A decision on `iso_smiles`,
`isomeric_smiles`, or `canonical_smiles` recalculates the selected chemistry
record after all its explicit property decisions have been projected. The
calculation uses that established precedence order and replaces calculated
masses, component chemistry, derived InChIKey fields, method metadata, and
errors as one unit. Reported masses, formulas, and InChIKeys remain source
evidence and change only through explicit decisions.

## Application and provenance

Published batches retain the ordinary immutable curation-batch provenance,
observed values, rationale, curator, and publication metadata. Harmonization
stages resolve and fingerprint the ChEBI bundle independently from the
metabolite-record bundle. Existing pipeline definitions remain unchanged;
new pipeline definitions include the ChEBI stream by default.

There is no migration of existing metabolite-record curations. If an existing
curation actually describes the ChEBI ontology record, it should be reviewed
and republished explicitly in the ChEBI bundle rather than reinterpreted in
place.
