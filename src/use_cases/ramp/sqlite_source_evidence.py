"""Source assertion extraction and disk-backed source row registration.

Evidence maps onto an already assigned RAMP identity; it never changes grouping.
"""
from collections import Counter

DETAIL_FIELDS = {
    'MetabolitePathwayEdge': 'source_id',
    'GenePathwayEdge': 'gene_id',
    'ProteinPathwayEdge': 'protein_id',
    'HmdbMetaboliteOntologyEdge': 'source_id',
    'MetaboliteClassificationEdge': 'source_id',
    'HmdbMetaboliteProteinAssociationEdge': 'source_id',
}


def evidence(collection, edge):
    """Return (provider, raw identifier, detail) assertions; never infer from endpoints."""
    field = DETAIL_FIELDS.get(collection)
    details = edge.get('details') if field else [dict(edge, source='rhea')]
    field = field or 'source_id'
    if not details:
        raise ValueError(f'{collection} {edge.get("start_id")} lacks source evidence; refresh the ingest before export')
    for detail in details:
        identifier, provider = detail.get(field), detail.get('source')
        if not isinstance(identifier, str) or not identifier.strip() or ':' not in identifier or not provider:
            raise ValueError(f'{collection} {edge.get("start_id")} lacks {field}/provider; refresh or repair the ingest before export')
        yield provider, identifier, detail


class SourceRows:
    """One row per source ID / RAMP ID / provider, independent of evidence order."""
    def __init__(self, writer):
        self.writer = writer
        self.pending = []
        writer.db.execute('''CREATE TEMP TABLE ramp_source_rows (
            sourceId TEXT, rampId TEXT, IDtype TEXT, geneOrCompound TEXT,
            commonName TEXT, priorityHMDBStatus TEXT, dataSource TEXT, pathwayCount INTEGER,
            name_rank INTEGER, PRIMARY KEY(sourceId,rampId,dataSource))''')

    def add(self, identifier, rid, kind, provider, name, status, name_rank=0):
        self.pending.append((identifier, rid, identifier.split(':', 1)[0], kind,
                             name, status, provider, -1, name_rank if name else 999))
        if len(self.pending) >= 5000:
            self.flush()

    def flush(self):
        self.writer.db.executemany('''INSERT INTO ramp_source_rows VALUES (?,?,?,?,?,?,?,?,?)
            ON CONFLICT(sourceId,rampId,dataSource) DO UPDATE SET
            commonName = CASE WHEN excluded.name_rank < name_rank OR
                (excluded.name_rank = name_rank AND excluded.commonName < commonName)
                THEN excluded.commonName ELSE commonName END,
            name_rank = MIN(name_rank,excluded.name_rank)''', self.pending)
        self.pending.clear()

    def finish(self):
        self.flush()
        for row in self.writer.db.execute('SELECT sourceId,rampId,IDtype,geneOrCompound,commonName,priorityHMDBStatus,dataSource,pathwayCount FROM ramp_source_rows ORDER BY sourceId,rampId,dataSource'):
            self.writer.add('source', **dict(zip(self.writer.columns['source'], row)))
        self.writer.db.execute('DROP TABLE ramp_source_rows')


def representative_compound_names(db):
    """Choose each compound's most frequent source-row name, as legacy RaMP did."""
    result = {}
    current = None
    name_counts = Counter()
    source_id_counts = Counter()
    spellings = {}
    source_spellings = {}

    def selected():
        counts = name_counts if name_counts else source_id_counts
        if not counts:
            return current
        highest = max(counts.values())
        winners = (value for value, count in counts.items() if count == highest)
        if name_counts:
            key = min(winners, key=lambda value: (len(spellings[value]), value))
            return spellings[key]
        key = min(winners)
        return source_spellings[key]

    for ramp_id, source_id, name in db.execute(
        'SELECT rampId, sourceId, commonName FROM "_raw_source" '
        "WHERE geneOrCompound='compound' ORDER BY rampId, sourceId, dataSource"
    ):
        if ramp_id != current:
            if current is not None:
                result[current] = selected()
            current = ramp_id
            name_counts.clear()
            source_id_counts.clear()
            spellings.clear()
            source_spellings.clear()
        normalized_id = source_id.lower()
        source_id_counts[normalized_id] += 1
        source_spellings.setdefault(normalized_id, source_id)
        if name not in (None, '', 'NA', 'None'):
            normalized = name.lower()
            name_counts[normalized] += 1
            spellings.setdefault(normalized, name)
    if current is not None:
        result[current] = selected()
    return result
