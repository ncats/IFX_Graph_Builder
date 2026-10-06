"""RaMP's explicit, input-driven gene/protein collapsing policy.

Only matches for identifiers present in the RaMP inputs create connections.
This module deliberately has no access to the resolver's wider alias catalog.
"""
from dataclasses import dataclass
from typing import Mapping, Iterable

TypedIdentifier = tuple[str, str]
POLICY_NAME = 'input_identifier_connected_components'
POLICY_VERSION = 2
POLICY = ('Collapse canonical UniProt matches connected by identifiers present in RaMP inputs; '
          'merge transitively, ignore unused aliases, and retain unmatched typed singletons')


@dataclass(frozen=True)
class GeneGroups:
    input_to_group: dict[TypedIdentifier, TypedIdentifier]
    canonical_accessions: dict[TypedIdentifier, tuple[str, ...]]


def collapse_input_matches(matches: Mapping[TypedIdentifier, Iterable[str]]) -> GeneGroups:
    """Connect all returned matches per input, never all aliases per protein.

    Canonical accessions need not themselves occur as input nodes. They can
    connect two used identifiers, but unused aliases cannot create connections.
    Unmatched input keys preserve their original collection/type.
    """
    parents = {}

    def root(accession):
        parents.setdefault(accession, accession)
        current = accession
        while parents[current] != current:
            current = parents[current]
        while accession != current:
            parent = parents[accession]
            parents[accession] = current
            accession = parent
        return current

    normalized = {key: tuple(sorted(set(values))) for key, values in matches.items()}
    for values in normalized.values():
        if not values:
            continue
        for value in values:
            left, right = sorted((root(values[0]), root(value)))
            parents[right] = left

    accessions = {}
    for accession in sorted(parents):
        accessions.setdefault(('UniProtKB', root(accession)), []).append(accession)
    inputs = {key: ('UniProtKB', root(values[0])) if values else key
              for key, values in sorted(normalized.items())}
    for key in inputs.values():
        accessions.setdefault(key, [])
    return GeneGroups(inputs, {key: tuple(values) for key, values in sorted(accessions.items())})
