# Metabolite deny-list multicomponent review

Status: read-only pre-release audit; no curation batches or manifests were
changed.

## Question

Some equivalence removals may have been curated in response to molecular-weight
spread even though the difference is explained by a counterion, hydrate,
repeated molecule, or another multicomponent representation. Those differences
should remain visible as MW warnings without automatically implying that the
equivalence edge is wrong.

## Validation semantics

The current harmonization workbench already treats molecular-weight spread as
a non-blocking warning. `_build_harmonization_stage_mw_validation` records
`warning_count` and warning examples, while a successfully materialized stage
is stored with `status: complete`. No MW validation result automatically
removes an edge or fails a pipeline.

Therefore, no validation-severity code change is required to allow these
cliques. The cleanup target is the published equivalence curation stream.

## Audit basis

- Curation manifest: `metabolite_equivalence_edges`, revision 21
- Active resolved decisions: 1,460
- Active `remove_edge` decisions: 1,449
- Active `retain_edge` decisions: 11
- Pre-curation stage: `baseline-eb6c4290c65248a3`
- Immediately post-curation stage: `stage-01-53ce1914ac6a1220`

The pre-curation stage was used as the source of clique membership. Endpoint
chemistry was then checked against the entire small pre-curation clique and the
two components produced by applying curations.

### Effect of the active removals

| Result | Pair count |
|---|---:|
| Endpoints in the same pre-curation clique | 1,377 |
| Endpoints separated after curations | 592 |
| Endpoints still connected after curations through another path | 785 |
| Endpoints not together in the pre-curation materialization | 72 |

The 592 effective splits arose from 411 distinct pre-curation cliques. Of
those, 273 cliques had a molecular-weight spread greater than 10%, accounting
for 439 separated deny-list pairs. This establishes the review population but
does not prove that MW spread was the curator's reason: the consolidated
pre-release baseline operations do not contain notes.

## High-signal multicomponent candidates

All eight candidates below are active removals from
`pre-release-equivalence-denylist-baseline-f98ae4a78a6766a9`. Each pair was in
one clique before curations, is split by the curation step, and had a baseline
MW spread above 10%.

| Denied pair | Pre-curation chemistry | MW spread | Assessment |
|---|---|---:|---|
| CHEBI:46659 ↔ HMDB:HMDB0014476 | Ipratropium bromide, 412.368 Da ↔ ipratropium cation, 332.222 Da | 24.1% | Strong retain candidate: bromide counterion is the difference; the clique also contains ChEBI, PubChem, and RefMet records for the cation |
| CHEBI:72449 ↔ REFMET:RM0199390 | Malachite green chloride, 364.171 Da ↔ malachite green cation, 329.202 Da | 10.8% | Strong retain candidate: chloride counterion is the difference; PubChem agrees with the RefMet cation |
| HMDB:HMDB0030524 ↔ PUBCHEM.COMPOUND:21831736 | One C7H5NS2 unit, 166.986 Da ↔ two identical units, 333.973 Da | 100.3% | Strong retain candidate under a stoichiometry-insensitive identity policy; PubChem 697993 and the rest of the clique represent the single unit |
| HMDB:HMDB0031543 ↔ PUBCHEM.COMPOUND:16211403 | One C6H8O2 unit, 112.052 Da ↔ two C6H8O2 components represented in keto/enol forms, 224.105 Da | 100.1% | Retain candidate if duplicated tautomeric components are intentionally collapsed |
| HMDB:HMDB0034885 ↔ PUBCHEM.COMPOUND:18772438 | One C15H18N2O3 unit, 274.132 Da ↔ two units, 548.263 Da | 100.1% | Strong retain candidate under a stoichiometry-insensitive identity policy |
| HMDB:HMDB0034913 ↔ PUBCHEM.COMPOUND:21873022 | One C16H20N2O3 unit, 288.147 Da ↔ two units, 576.295 Da | 100.1% | Strong retain candidate under a stoichiometry-insensitive identity policy; PubChem 54744 represents the single unit |
| HMDB:HMDB0041886 ↔ PUBCHEM.COMPOUND:6917719 | Enalaprilat anhydrous, 348.169 Da ↔ enalaprilat dihydrate, 384.190 Da | 10.4% | Strong retain candidate: two water components are the difference; ChEBI, PubChem 5462501, and RefMet agree on the anhydrous parent |
| HMDB:HMDB0041934 ↔ PUBCHEM.COMPOUND:49776908 | C9H13N3O6 compound, 259.080 Da ↔ that component plus a distinct C5H5N5O component, 410.130 Da | 58.4% | Manual review: this is genuinely multicomponent, but the added component is not a counterion, water, or duplicate, so equivalence depends on the intended policy for complexes/mixtures |

The first seven are good candidates for explicit `retain_edge` decisions after
curator confirmation. The eighth should not be automatically restored merely
because its structure is multicomponent.

## Detection limits and next pass

This high-confidence screen required at least one of:

- the same standardized parent structure after removing small fragments;
- a repeated structure component;
- one formula appearing as a dot-separated component of the other; or
- an exact two- or three-fold formula relationship.

It intentionally favors precision. It can miss deny-list pairs whose endpoints
have no structure or formula even when other members of their pre-curation
components provide enough indirect evidence. A broader candidate report should
compare the chemistry consensus of the two post-curation components within
each of the 273 affected pre-curation MW-warning cliques, while keeping
salt/hydrate/duplicate cases separate from complexes and unrelated source-data
errors.
