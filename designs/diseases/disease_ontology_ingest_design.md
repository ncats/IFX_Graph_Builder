# Disease Ontology ingest design

## Purpose and current source

Ingest Disease Ontology (DO) terms and direct parent relationships into the Pharos and Target Graph builds. The current input is the Registry snapshot `disease_ontology:ontology:2026-09-30`, containing `doid.obo` (SHA-256 `51e717b6b9f5d391f2793f9e16faf9c8629a661ab8f21c3a024458129f182ad2`). The adapter requests `doid.obo` by name and takes version metadata from the Registry snapshot.

The upstream `doid.obo` and the Harmonizers handoff previously named `HumanDO.obo` have the same bytes. ODIN uses the upstream filename because Disease Ontology describes its separate historical `HumanDO.obo` product as equivalent to `doid-non-classified.obo`.

## Prior ETL reference

The old Pharos319 TCRD pipeline had 9,233 `do` rows, 10,880 `do_parent` rows, and 85,072 `ancestry_do` rows. Those counts reflect an older DO release. `ancestry_do` is the transitive closure computed from direct `do_parent` links.

## Mapping

The OBO reader uses `obonet`. It emits active `DOID:` terms as `Disease` nodes, preserving the term ID, name, definition (`do_description`), scoped synonyms, comments, and `CLASS` type. Synonyms are sorted by OBO scope (broad, exact, narrow, related) and value. The DOID resolver may map a term to a canonical MONDO `Disease`; DO-only terms retain their DOID IDs.

`DODiseaseParentEdgeAdapter` emits direct `is_a` relationships when both endpoints are active DOID terms. The edge uses `DODiseaseParentEdge`, distinct from the MONDO hierarchy edge, with `source: DO`. Transitive ancestry is calculated downstream. The graph adapters do not emit cross-ontology xrefs; identifier reconciliation belongs in the resolver.

The same OBO reader feeds `DOTableAdapter` and `DOTableParentEdgeAdapter` for the Pharos `do` and `do_parent` tables. The table path is configured in `src/use_cases/pharos/tcrd.yaml`. The graph adapters are configured in `src/use_cases/pharos/pharos.yaml` and `target_graph.yaml`.

## 2026-09-30 OBO cutover validation

The cutover compared the Registry's former `doid.json` with the exact OBO file. Both yielded 12,335 active DOID terms and 17,479 direct parent links. IDs, names, definitions, comments, parent pairs, and synonym sets matched. Scope-and-value sorting reproduced the former JSON synonym order for every term in this release.

No graph ETL or MySQL conversion was run as part of the OBO code change. For graph validation, run a narrow build with `DODiseaseAdapter` and `DODiseaseParentEdgeAdapter`, then inspect node and edge counts and parent endpoints. For the TCRD path, include `DOTableAdapter` and `DOTableParentEdgeAdapter` and verify `do`, `do_parent`, and `ancestry_do`.

## Open decisions

- Whether to ingest DO cross-ontology xrefs (MESH, NCI, OMIM, and others) as explicit references.
- Whether to link DOID and MONDO through shared xrefs in post-processing.
