"""Derive legacy RaMP cofactor flags from the pinned ChEBI role graph."""

from collections import defaultdict


COFACTOR_ROLE = "CHEBI:23357"


def cofactor_chemical_ids(reader):
    """Return chemicals below the cofactor role, including chemical subclasses.

    ChEBI ``IsAEdge`` points from child to parent; ``HasBiologicalRoleEdge``
    points from a chemical entity to its role. The caller tests the reported
    ChEBI participant ID, never the harmonized RaMP metabolite group.
    """
    roles = {doc["id"] for doc in reader.records("BiologicalRole")}
    if COFACTOR_ROLE not in roles:
        raise ValueError(f"ChEBI cofactor role {COFACTOR_ROLE} is missing from the graph")

    children = defaultdict(set)
    for edge in reader.records("IsAEdge"):
        child, parent = edge["start_id"], edge["end_id"]
        children[parent].add(child)

    descendants = {COFACTOR_ROLE}
    pending = [COFACTOR_ROLE]
    while pending:
        for child in children[pending.pop()]:
            if child in roles and child not in descendants:
                descendants.add(child)
                pending.append(child)

    chemicals = set()
    for edge in reader.records("HasBiologicalRoleEdge"):
        if edge["end_id"] in descendants:
            chemicals.add(edge["start_id"])
    pending = list(chemicals)
    while pending:
        for child in children[pending.pop()]:
            if child not in chemicals:
                chemicals.add(child)
                pending.append(child)
    return chemicals, descendants
