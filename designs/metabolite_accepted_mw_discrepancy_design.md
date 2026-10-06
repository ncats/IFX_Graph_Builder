# Accepted Metabolite MW Discrepancies

## Status

Implemented in the QA Browser curation workflow. MW adjudications are stored
in the dedicated `metabolite_mw_adjudications` v2 curation stream and are
overlaid on immutable saved-stage validation results.

## Problem

Molecular-weight validation can correctly detect a large spread while still
being unable to determine that the spread is chemically expected. A common
case is a clique containing both a parent compound and a salt, hydrate,
solvate, or other multicomponent form.

The component-aware validator can explain such a difference when at least one
source provides a parseable multicomponent structure. Some source records,
however, provide only a name, formula, identifier, and whole-record mass. The
mass is valid for the entity represented by that source, but the graph does not
contain enough structure on that record to decompose it automatically.

Removing an equivalence edge or correcting the source mass would be the wrong
response when expert review concludes that the clique and source assertions
are valid. Leaving the finding as an unresolved error also misstates the
reviewed quality of the harmonization result.

## Motivating example

The reviewed clique containing miconazole and miconazole nitrate includes, in
part:

- `PUBCHEM.COMPOUND:4189`, a structured miconazole record;
- `PUBCHEM.COMPOUND:68553`, a structured miconazole nitrate record whose
  components expose miconazole and nitric acid;
- `REFMET:RM0200855`, miconazole;
- `REFMET:RM0137167`, miconazole nitrate with a nitrate-inclusive mass but no
  structure from which ODIN can derive components.

The larger RefMet mass is not evidence that its mapping is wrong. The source
record represents the nitrate form, and the structured PubChem neighbor
provides the explanation that the RefMet record cannot express on its own.
Differences between exact-mass and monoisotopic-mass conventions may remain
visible but are not the identity decision being adjudicated here.

## Semantic decision

Support a durable **accepted MW discrepancy** adjudication. It means:

> A reviewer has examined the current molecular-weight evidence and accepts
> the observed spread as chemically explained for this harmonized group.

Acceptance is validation metadata only. It must not:

- remove, retain, or create identifier-mapping edges;
- replace or suppress source chemical properties;
- imply that every member has the same molecular formula or whole-record mass;
- make the discrepancy disappear from the audit trail.

An accepted discrepancy remains visible as a reviewed warning, is removed from
the unresolved-error count, and is counted separately from both errors and
automatic component-match warnings. A curator can reopen it through the same
cart and publication workflow.

## Curation shape

This belongs in a dedicated validation-adjudication stream rather than
`metabolite_record_properties` or `metabolite_equivalence_edges`:

```json
{
  "action": "accept_mw_discrepancy",
  "target": {
    "kind": "validation_finding",
    "curation_set": "metabolite_harmonization",
    "check": "mw_spread",
    "finding_id": "mw-<finding hash>",
    "anchor_id": "PUBCHEM.COMPOUND:4189"
  },
  "reason": "salt_or_counterion",
  "supporting_ids": ["PUBCHEM.COMPOUND:68553"],
  "observed_member_ids": [
    "PUBCHEM.COMPOUND:4189",
    "PUBCHEM.COMPOUND:68553",
    "REFMET:RM0137167",
    "REFMET:RM0200855"
  ],
  "observed_evidence_fingerprint": "<content hash>",
  "note": "RefMet reports the nitrate-inclusive mass; PubChem exposes the parent and nitric-acid components."
}
```

The target does not depend on clique rank, stage ID, or another rebuild-specific
identifier. Decisions are evaluated from the immutable publication history,
newest first. The newest decision whose observed members overlap the current
finding controls that lineage: an exact-fingerprint acceptance is accepted, a
changed-fingerprint acceptance is stale, and a reopen is unresolved. Keeping
history rather than collapsing by anchor preserves both descendants when an
accepted clique later splits.

Candidate reason codes include:

- `salt_or_counterion`;
- `hydrate_or_solvate`;
- `multicomponent_form`;
- `source_mass_convention`;
- `protonation_or_charge_state`;
- `polymer_or_repeat_unit`;
- `other_reviewed_explanation`.

A later decision must be able to retire or reopen an acceptance without
deleting its history.

## Staleness and safety

Acceptance must be bounded to the evidence the reviewer saw. A bare clique ID
or member list is too broad: clique membership and chemical properties can
change between source releases and pipeline stages.

The adjudication should therefore carry a deterministic evidence fingerprint
covering the relevant members and their mass-bearing evidence, including as
available:

- source and source record ID;
- reported average, exact, or monoisotopic mass and its declared channel;
- molecular formula;
- structure-derived whole and component masses;
- the structured record used to explain the discrepancy.

If membership or evidence changes, the acceptance becomes stale and the
validation finding reopens for review. Membership is intentionally strict:
any identifier joining or leaving the clique invalidates the reviewed
fingerprint, even when that identifier carries no MW evidence.

## Presentation expectations

The QA Browser should distinguish at least:

- unresolved MW error;
- automatically explained component-match warning;
- reviewed and accepted MW discrepancy;
- stale acceptance requiring renewed review.

The reviewer should be able to see the original mass clusters, component
evidence, rationale, supporting records, curator, publication time, and the
stage against which the decision is currently being evaluated.

Saved stages keep full scalar counts plus a lightweight identity/fingerprint
index for every MW finding; detailed display evidence remains capped. This lets
accepted and stale totals remain exact even when a finding is outside the
displayed sample. Older stages without the index retain their stored totals and
can still use their persisted display findings.

Potential future work is limited to promoting recurring human explanations
into automatic validator rules and deciding whether exact mass warrants its own
source channel. Those changes do not alter the accepted-discrepancy contract.
