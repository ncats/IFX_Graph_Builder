# ChEBI data-quality findings

## Purpose

This report records internally inconsistent ChEBI records found while reviewing
molecular-weight validation warnings in the metabolite harmonization graph. It
is intended both as a source-feedback stub for ChEBI and as input to defensive
chemistry handling in IFX_ODIN.

The findings below were confirmed against ChEBI release 255, loaded on
2026-09-09, and the corresponding LIPID MAPS records. Structure-derived values
were independently recalculated with RDKit 2026.03.6.

## Deuterium represented as methyl substituents

The first four confirmed errors share the same signature:

- The ChEBI name, synonym, isotopic InChI, InChIKey, and any LIPID MAPS
  cross-reference identify the intended deuterated compound.
- The ChEBI SMILES has a carbon atom where the corresponding LIPID MAPS SMILES
  has an explicit deuterium (`[2H]`).
- Relative to the ordinary-hydrogen molecular formula, every affected isotope
  site therefore adds `CH2`: `D4` produces an erroneous `+C4H8`, `D5` produces
  `+C5H10`, and so on.
- ChEBI's formula and reported masses are internally consistent with that
  corrupted SMILES, but inconsistent with the preserved isotopic InChI and
  InChIKey.

| ChEBI ID | Intended entity / corroborating record | Correct formula | Correct monoisotopic mass | ChEBI reported formula / mass | Structure-key evidence | Assessment |
|---|---|---:|---:|---:|---|---|
| CHEBI:137828 | Arachidonic acid-d8 / LIPIDMAPS:LMFA01030003 | C20H24D8O2 | 312.290444 | C28H48O2 / 416.36543 | Preserved key `YZXBAPSDXZZRGB-FBFLGLDDSA-N`; corrupted SMILES derives `RLKSKSKUURKFND-FKJQSWGXSA-N` | Eight deuterium sites became eight methyl substituents |
| CHEBI:137735 | Linoleic acid-d4 / LIPIDMAPS:LMFA01030809 | C18H28D4O2 | 284.265337 | C22H40O2 / 336.30283 | Preserved key `OYHQOLUKZRVURQ-GPTDOFLBSA-N`; corrupted SMILES derives `NWLFFMSMSQPLAF-AUYXYSRISA-N` | Four deuterium sites became four methyl substituents |
| CHEBI:137823 | 9,14,19,19,19-pentadeuterio-1alpha,25-dihydroxyprevitamin D3 / LIPIDMAPS:LMST03020097 | C27H39D5O3 | 421.360429 | C32H54O3 / 486.4073 | Preserved key `DOIZGAFWGREMOD-SGWDTOPTSA-N`; corrupted SMILES derives `NAQZUBHPYOOMOP-WIQRVGQPSA-N` | Five deuterium sites became five methyl substituents |
| CHEBI:139244 | gemini-0097 (hexadeutero vitamin-D analogue) | C32H38D6F6O4 | 612.352040 | C38H56F6O4 / 690.40828 | Preserved key `VSOWXEHVBCDXAY-CVRMBSQLSA-N`; corrupted SMILES derives `UQBNMLKSOFEAQX-DAPZISOPSA-N` | Six deuterium sites became six methyl substituents |

The first three ChEBI records are 2-star entries and reach the harmonization
graph through the ChEBI FULL ontology `ChemicalEntity` path. `CHEBI:139244` is
a 3-star entry and demonstrates that the same defect also reaches the ChEBI
SDF chemistry path.

## Harmonization assessment

The LIPID MAPS-to-ChEBI equivalences above should not be denied. For the mapped
pairs, ChEBI's name, isotopic InChI, complete InChIKey, and explicit LIPID MAPS
xref all corroborate the identity. The invalid values are ChEBI's SMILES,
formula, average mass, and monoisotopic mass.

These cases also show that agreement between a reported formula/mass and a
reported SMILES is not sufficient validation: all three can originate from the
same corrupted structure representation. A disagreement between the reported
complete InChIKey and the InChIKey derived from SMILES is the stronger warning.

## Recommended follow-up

1. Profile all ChEBI records for disagreement between the reported InChIKey and
   the key independently derived from SMILES, with special attention to
   isotopically modified compounds.
2. When the reported key and SMILES-derived key disagree, quarantine the
   SMILES-derived formula and masses instead of using them for merging or MW
   validation.
3. Evaluate deriving chemistry from a valid isotopic InChI as a fallback. Keep
   the source fields and the reason for selecting or rejecting each value.
4. Send the affected IDs, preserved InChIs, corrupted SMILES, and corrected
   calculated values to the ChEBI maintainers.
