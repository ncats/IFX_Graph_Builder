# HMDB XML versus SDF chemistry profile

## Purpose

Investigate why `HMDB:HMDB0009911` and other HMDB identifiers directly linked
to generic `KEGG.COMPOUND:C00626` have no `chem_props` after a clean metabolite
harmonization rebuild.

## Inputs

Both artifacts are registered as HMDB 5.0 with version date 2021-11-17 and
download date 2026-06-30.

- `hmdb:metabolites_xml:5.0`: `hmdb_metabolites.zip` (953,913,010 bytes)
- `hmdb:structures_sdf:5.0`: `structures.zip` (96,479,839 bytes)

The comparison used primary HMDB accessions and the following corresponding
fields:

| XML | SDF |
|---|---|
| `chemical_formula` | `FORMULA` |
| `average_molecular_weight` | `MOLECULAR_WEIGHT` |
| `monisotopic_molecular_weight` | `EXACT_MASS` |
| `smiles` | `SMILES` |
| `inchi` | `INCHI_IDENTIFIER` |
| `inchikey` | `INCHI_KEY` |

## Accession coverage

| Population | Records |
|---|---:|
| XML | 217,920 |
| SDF | 217,776 |
| Present in both | 217,776 |
| XML only | 144 |
| SDF only | 0 |
| Duplicate primary accessions | 0 in either artifact |

The SDF accession population is therefore a strict subset of the XML
population. Before the XML fallback was added, all 144 XML-only identifiers
existed in the rebuilt graph, but none had `chem_props`; all consequently had
an unknown generic-structure classification.

Of the 144 XML-only records:

- 123 have formula, average mass, monoisotopic mass, and InChIKey;
- 120 also have SMILES;
- 119 also have InChI; and
- 21 have none of the six chemistry fields.

## Agreement on the 217,776 shared accessions

| Field | Exact/normalized agreement | Rate |
|---|---:|---:|
| Formula | 217,737 / 217,776 | 99.9821% |
| InChI | 217,737 / 217,776 | 99.9821% |
| InChIKey | 217,737 / 217,776 | 99.9821% |
| Average mass | 217,717 / 217,756 numeric pairs | 99.9821% |
| Monoisotopic mass | 217,717 / 217,756 numeric pairs | 99.9821% |
| SMILES text | 217,014 / 217,776 | 99.6501% |

The 762 SMILES text differences are mostly alternate serializations: 723 have
the same reported InChIKey. The remaining 39 coincide with the formula, InChI,
InChIKey, and mass disagreement population.

Of those 39 structural-property disagreements:

- 26 retain the same first two InChIKey blocks. They are mostly protonation
  differences: the XML adds a proton layer and the final InChIKey block changes
  from `N` to `O` while connectivity and stereochemistry remain the same.
- 13 have different InChIKey connectivity blocks and represent genuinely
  different structures. Ten of these are concentrated in
  `HMDB0241915`–`HMDB0241924`, with the other three at `HMDB0013132`,
  `HMDB0034051`, and `HMDB0242435`. Several SDF formulas are visibly
  incompatible with their HMDB names, so these should be treated as source
  inconsistencies rather than reconciled automatically.

The maximum shared-record mass disagreement is about 275 Da, but only the same
39 records disagree at all; the other 217,717 numeric pairs are exactly equal.

## Internal validation of XML-only chemistry

The 123 XML-only records with chemistry are internally consistent:

- all 123 formulas were parseable;
- every reported monoisotopic mass was within 0.00001 Da of the value
  calculated from its formula;
- median absolute formula/mass difference was 0.000000092 Da;
- all 120 SMILES values parsed successfully with RDKit; and
- every RDKit-derived InChIKey exactly matched the XML-reported InChIKey.

All 23 chemistry-free graph records directly linked to
`KEGG.COMPOUND:C00626`, including `HMDB0009911`, belong to this XML-only set.
All 23 have complete core XML chemistry, formula-derived masses within 0.001
Da, and exact reported-versus-derived InChIKey agreement.

For `HMDB0009911`, XML reports:

- formula `C49H83O13P`;
- average molecular weight `911.1493`;
- monoisotopic molecular weight `910.557129254`;
- InChIKey `OMBBFDFDCMFWTB-YMBWUVAOSA-N`; and
- complete SMILES and InChI values.

## Recommended ingestion rule

Keep the SDF chemistry for every accession present in the SDF, preserving the
existing behavior and avoiding any automatic choice among the 39 conflicting
shared records. For primary HMDB accessions absent from the SDF, use the XML
chemistry fields as a fallback. This would add trustworthy chemistry to 123
currently unclassified graph identifiers, including all 23 affected
`C00626` neighbors, while leaving the 21 XML-only records with no source
chemistry unknown.

Implementation should ensure the fallback is selected by accession coverage,
not by individual missing fields on overlapping records, so a single HMDB
identifier still receives one coherent source chemistry record.
