"""Legacy-compatible sparse pathway overlap data from an existing RaMP SQLite."""

from bisect import bisect_right
from collections import defaultdict
import time
import zlib


# Each scope has its own ordered pathway universe and therefore its own indices.
SCOPES = (
    ('metabolite', 'metabolite_blob', 'metabolite_count', 5, "a.rampId LIKE 'RAMP_C%'"),
    ('gene', 'gene_blob', 'gene_count', 5, "a.rampId LIKE 'RAMP_G%'"),
    ('analyte', 'analyte_blob', None, 10, '1=1'),
)


def load_scope(db, *, cutoff, analyte_filter):
    """Load distinct pathway memberships in the exact legacy sort order."""
    query = f'''
        SELECT a.pathwayRampId,a.rampId
        FROM analytehaspathway a JOIN pathway p ON p.pathwayRampId=a.pathwayRampId
        WHERE p.type!='hmdb' AND {analyte_filter}
        GROUP BY a.pathwayRampId,a.rampId
        ORDER BY a.pathwayRampId,a.rampId'''
    pathways = []
    members = []
    current = None
    ids = []
    for pathway_id, analyte_id in db.execute(query):
        if pathway_id != current:
            if current is not None and len(ids) >= cutoff:
                pathways.append(current)
                members.append(frozenset(ids))
            current, ids = pathway_id, []
        ids.append(analyte_id)
    if current is not None and len(ids) >= cutoff:
        pathways.append(current)
        members.append(frozenset(ids))
    return pathways, members


def duplicate_pairs(pathways, members):
    """Return one pair per identical combined-analyte pathway membership."""
    seen = defaultdict(list)
    for pathway_id, analytes in zip(pathways, members, strict=True):
        seen[analytes].append(pathway_id)
    return [(left, right) for group in seen.values() if len(group) > 1
            for i, left in enumerate(group) for right in group[i + 1:]]


def sparse_rows(pathways, members, *, label, progress=None, progress_every=5000):
    """Yield compressed rows without allocating a dense pathway matrix."""
    postings = defaultdict(list)
    for index, analytes in enumerate(members):
        for analyte in analytes:
            postings[analyte].append(index)
    sizes = [len(analytes) for analytes in members]
    started = time.monotonic()
    last_logged = started
    positive_pairs = 0
    for index, (pathway_id, analytes) in enumerate(zip(pathways, members, strict=True)):
        overlaps = defaultdict(int)
        for analyte in analytes:
            posting = postings[analyte]
            for partner in posting[bisect_right(posting, index):]:
                overlaps[partner] += 1
        pieces = [f'{index},-']
        last_index = index
        for partner in sorted(overlaps):
            shared = overlaps[partner]
            score = round(1000 * shared / (sizes[index] + sizes[partner] - shared))
            if score > 0:
                pieces.append(f'{partner - last_index},{score}')
                last_index = partner
                positive_pairs += 1
        blob = zlib.compress('|'.join(pieces).encode('ascii'))
        yield pathway_id, sizes[index], blob
        done = index + 1
        now = time.monotonic()
        if progress and (done % progress_every == 0 or done == len(pathways)
                         or now - last_logged >= 30):
            progress(f'{label}: {done:,}/{len(pathways):,} pathways '
                     f'({done / len(pathways):.0%}); {positive_pairs:,} nonzero pairs; '
                     f'{now - started:.1f}s')
            last_logged = now


def populate_temp_tables(db, *, progress=None):
    """Build scope blobs and duplicate pairs in temporary tables for atomic replacement."""
    db.execute('''CREATE TEMP TABLE _ramp_similarity_new (
        pathwayRampId TEXT PRIMARY KEY, analyte_blob BLOB,
        metabolite_blob BLOB, gene_blob BLOB,
        metabolite_count INTEGER, gene_count INTEGER)''')
    db.execute('''CREATE TEMP TABLE _ramp_duplicates_new (
        pathwayRampId1 TEXT NOT NULL, pathwayRampId2 TEXT NOT NULL)''')
    summaries = {}
    for label, blob_column, count_column, cutoff, analyte_filter in SCOPES:
        started = time.monotonic()
        pathways, members = load_scope(db, cutoff=cutoff, analyte_filter=analyte_filter)
        if progress:
            progress(f'{label}: {len(pathways):,} eligible pathways (minimum {cutoff} '
                     f'distinct {label} IDs); calculating sparse overlaps')
        if label == 'analyte':
            pairs = duplicate_pairs(pathways, members)
            db.executemany('INSERT INTO _ramp_duplicates_new VALUES (?,?)', pairs)
        sql = (f'''INSERT INTO _ramp_similarity_new (pathwayRampId,{blob_column}'''
               + (f',{count_column}' if count_column else '') + ') VALUES ('
               + ','.join('?' for _ in range(3 if count_column else 2)) + ') '
               + f'ON CONFLICT(pathwayRampId) DO UPDATE SET {blob_column}=excluded.{blob_column}'
               + (f',{count_column}=excluded.{count_column}' if count_column else ''))
        batch = []
        for pathway_id, count, blob in sparse_rows(
                pathways, members, label=label, progress=progress):
            batch.append((pathway_id, blob, count) if count_column else (pathway_id, blob))
            if len(batch) == 500:
                db.executemany(sql, batch)
                batch.clear()
        if batch:
            db.executemany(sql, batch)
        summaries[label] = {'eligible_pathways': len(pathways),
                            'cutoff': cutoff, 'elapsed_seconds': round(time.monotonic() - started, 1)}
        if label == 'analyte':
            summaries['duplicate_pairs'] = len(pairs)
        if progress:
            progress(f'{label}: finished in {summaries[label]["elapsed_seconds"]:.1f}s')
    orphan = db.execute('''SELECT pathwayRampId FROM _ramp_similarity_new
        WHERE metabolite_count IS NULL AND gene_count IS NULL LIMIT 1''').fetchone()
    if orphan:
        raise ValueError(f'Combined-only pathway is absent from the legacy row union: {orphan[0]}')
    summaries['rows'] = db.execute('SELECT count(*) FROM _ramp_similarity_new').fetchone()[0]
    return summaries
