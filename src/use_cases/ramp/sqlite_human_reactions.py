"""Legacy RaMP human-scope flags from original graph evidence."""

from collections import defaultdict
import gzip
import hashlib
import ijson


HUMAN_CHEBI_ROOTS = frozenset({
    'CHEBI:77746', 'CHEBI:85234', 'CHEBI:84087', 'CHEBI:76967', 'CHEBI:49103',
})


def _records(graph, collection, fields):
    select = ','.join(repr(field) for field in fields)
    cursor = graph.aql.execute(
        f'FOR d IN @@collection RETURN KEEP(d,{select})',
        bind_vars={'@collection': collection}, batch_size=5000, stream=True,
        max_runtime=3600,
    )
    try:
        yield from cursor
    finally:
        cursor.close(ignore_missing=True)


def human_chebi_ids(graph):
    """Follow the old OBO mixed is_a/has_biological_role closure."""
    present = {doc['id'] for doc in _records(graph, 'BiologicalRole', ('id',))}
    missing = HUMAN_CHEBI_ROOTS - present
    if missing:
        raise ValueError(f'Human ChEBI roots are missing from graph: {sorted(missing)}')
    children = defaultdict(set)
    for collection in ('IsAEdge', 'HasBiologicalRoleEdge'):
        for edge in _records(graph, collection, ('start_id', 'end_id')):
            children[edge['end_id']].add(edge['start_id'])
    descendants = set(HUMAN_CHEBI_ROOTS)
    pending = list(HUMAN_CHEBI_ROOTS)
    while pending:
        for child in children[pending.pop()]:
            if child not in descendants:
                descendants.add(child)
                pending.append(child)
    return descendants


def human_primary_accessions(path, expected_sha256):
    with path.open('rb') as handle:
        digest = hashlib.file_digest(handle, 'sha256').hexdigest()
    if digest != expected_sha256:
        raise ValueError('UniProt file differs from the exact resolver input recorded in SQLite')
    accessions = set()
    with gzip.open(path, 'rb') as handle:
        for record in ijson.items(handle, 'results.item'):
            accession = record.get('primaryAccession')
            if not isinstance(accession, str) or not accession:
                raise ValueError('Pinned human UniProt entry lacks primaryAccession')
            accessions.add('UniProtKB:' + accession)
    if not accessions:
        raise ValueError('Pinned human UniProt input has no primary accessions')
    return accessions


def human_reaction_flags(graph, human_accessions):
    human_chemicals = human_chebi_ids(graph)
    reactions = {doc['id'] for doc in _records(graph, 'RheaReaction', ('id',))}
    if not reactions:
        raise ValueError('RheaReaction is empty in graph')
    flags = {reaction: [0, 1] for reaction in reactions}
    for edge in _records(graph, 'RheaMetaboliteReactionEdge', ('end_id', 'source_id')):
        reaction = edge['end_id']
        if reaction not in flags:
            raise ValueError(f'Rhea metabolite edge targets missing reaction {reaction}')
        raw = edge.get('source_id')
        if not raw:
            raise ValueError(f'Rhea metabolite edge lacks reported ChEBI ID for {reaction}')
        if raw not in human_chemicals:
            flags[reaction][1] = 0
    for edge in _records(graph, 'RheaProteinReactionEdge', ('end_id', 'source_id')):
        reaction = edge['end_id']
        if reaction not in flags:
            raise ValueError(f'Rhea protein edge targets missing reaction {reaction}')
        raw = edge.get('source_id')
        if not raw:
            raise ValueError(f'Rhea protein edge lacks reported UniProt ID for {reaction}')
        if raw in human_accessions:
            flags[reaction][0] = 1
    return {reaction: tuple(values) for reaction, values in flags.items()}, {
        'human_chebi_ids': len(human_chemicals),
        'human_primary_uniprot_accessions': len(human_accessions),
        'reactions': len(reactions),
        'has_human_prot': sum(values[0] for values in flags.values()),
        'only_human_mets': sum(values[1] for values in flags.values()),
    }
