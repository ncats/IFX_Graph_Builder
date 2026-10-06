"""Legacy RaMP HMDB status priority, applied only to active group members."""
from collections import Counter

PRIORITY = {'predicted': 1, 'expected': 2, 'detected': 3, 'quantified': 4}


def group_statuses(records, membership):
    statuses = {}
    counts = Counter(records_missing_field=0, records_with_status=0, records_without_status=0)
    for doc in records:
        identifier = doc['id']
        if identifier not in membership or not identifier.lower().startswith('hmdb:'):
            continue
        if 'hmdb_status' not in doc:
            counts['records_missing_field'] += 1
        value = doc.get('hmdb_status')
        if not value:
            counts['records_without_status'] += 1
            continue
        if value not in PRIORITY:
            raise ValueError(f'Unknown HMDB status {value!r} on {identifier}')
        counts['records_with_status'] += 1
        rid = membership[identifier]
        if PRIORITY[value] > PRIORITY.get(statuses.get(rid), 0):
            statuses[rid] = value
    return statuses, {'policy': 'quantified > detected > expected > predicted; propagate to every source row in group',
                      **counts, 'groups_with_status': len(statuses)}
