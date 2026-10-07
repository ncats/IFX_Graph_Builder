"""Read-only comparison of RaMP lookup and reviewed association tables."""
import argparse
from datetime import datetime, timezone
import html
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile


EXAMPLES = (
    ('metabolite', 'Metabolite example: D-glucose', 'hmdb:HMDB0000122', 'compound'),
    ('gene', 'Gene/protein example: EGFR', 'uniprot:P00533', 'gene'),
)
EXAMPLE_FIELDS = ('rampId', 'type', 'common_name')
SOURCE_FIELDS = ('sourceId', 'rampId', 'IDtype', 'geneOrCompound', 'commonName',
                 'priorityHMDBStatus', 'dataSource', 'pathwayCount')


def populated_expression(field):
    column = '"' + field.replace('"', '""') + '"'
    value = f'lower(trim(CAST({column} AS TEXT)))'
    excluded = ("'na','n/a'" if field in ('commonName', 'common_name') else
                "'no_hmdb_status'" if field == 'priorityHMDBStatus' else
                "'-1'" if field == 'pathwayCount' else None)
    return f'{column} IS NOT NULL AND trim(CAST({column} AS TEXT)) != \'\'' + (
        f' AND {value} NOT IN ({excluded})' if excluded else '')


def field_coverage(db, table, fields, columns, group_fields):
    """Count meaningful values for each displayed field within each row group."""
    select = ', '.join(group_fields)
    counts = ', '.join(f'SUM(CASE WHEN {populated_expression(field)} THEN 1 ELSE 0 END)'
                       if field in columns else 'NULL' for field in fields)
    query = f'SELECT {select}, count(*), {counts} FROM {table} GROUP BY {select}'
    return {tuple(row[:len(group_fields)]): {
        'total': row[len(group_fields)],
        'fields': {field: row[len(group_fields) + 1 + i] for i, field in enumerate(fields)},
    } for row in db.execute(query)}


def coverage_cell(coverage, field):
    if not coverage or not coverage['total']:
        return 'No rows'
    filled = coverage['fields'][field]
    if filled is None:
        return 'Column absent'
    total = coverage['total']
    return f'{filled:,} / {total:,} ({filled / total:.1%})'


def coverage_percent(coverage, field):
    if not coverage or not coverage['total'] or coverage['fields'][field] is None:
        return None
    return 100 * coverage['fields'][field] / coverage['total']


def change_style(delta, magnitude):
    if delta is None or delta == 0:
        return ''
    endpoint = (173, 222, 188) if delta > 0 else (246, 185, 185)
    strength = max(0.08, min(abs(magnitude), 100) / 100)
    rgb = ','.join(str(round(255 + (value - 255) * strength)) for value in endpoint)
    return f' style="background-color:rgb({rgb})"'


def association_coverage(db, label, table, target, required, query):
    columns = {row[1] for row in db.execute(f'PRAGMA table_info({table})')}
    if not columns:
        return False, {}
    if not required <= columns:
        raise ValueError(f'{label}: {table} is missing {sorted(required - columns)}')
    if not db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (target,)).fetchone():
        raise ValueError(f'{label}: {table} requires missing {target} table')
    return True, {key: {'rows': rows, 'analytes': analytes, 'targets': targets}
                  for key, rows, analytes, targets in db.execute(query)}


def profile(label, path):
    path = Path(path).resolve(strict=True)
    db = sqlite3.connect(path.as_uri() + '?mode=ro', uri=True)
    try:
        db.execute('PRAGMA query_only=ON')
        db.execute('BEGIN')
        columns = {row[1] for row in db.execute('PRAGMA table_info(source)')}
        required = {'sourceId', 'rampId', 'IDtype', 'geneOrCompound', 'dataSource'}
        if not required <= columns:
            raise ValueError(f'{label}: source is missing {sorted(required - columns)}')
        analyte_columns = [row[1] for row in db.execute('PRAGMA table_info(analyte)')]
        if not {'rampId', 'type'} <= set(analyte_columns):
            raise ValueError(f'{label}: analyte is missing rampId or type')
        analyte_coverage = field_coverage(db, 'analyte', EXAMPLE_FIELDS,
                                          set(analyte_columns), ('type',))
        source_field_counts = field_coverage(db, 'source', SOURCE_FIELDS,
                                             columns, ('trim(dataSource)', 'geneOrCompound'))
        ramp_ids = {
            'Metabolite RaMP IDs': db.execute("SELECT count(DISTINCT rampId) FROM analyte WHERE type='compound'").fetchone()[0],
            'Gene/protein RaMP IDs': db.execute("SELECT count(DISTINCT rampId) FROM analyte WHERE type='gene'").fetchone()[0],
        }
        source_coverage = {'metabolite_sources': {}, 'gene_sources': {}}
        for provider, kind, count in db.execute('''
                SELECT coalesce(nullif(trim(dataSource),''),'(missing provider)'),
                       geneOrCompound,count(DISTINCT rampId)
                FROM source GROUP BY coalesce(nullif(trim(dataSource),''),'(missing provider)'),
                                     geneOrCompound'''):
            section = 'metabolite_sources' if kind == 'compound' else 'gene_sources' if kind == 'gene' else None
            if section:
                source_coverage[section][provider] = count
        examples = {}
        for key, title, source_id, kind in EXAMPLES:
            matched = [row[0] for row in db.execute(
                'SELECT DISTINCT rampId FROM source WHERE sourceId=? AND geneOrCompound=?',
                (source_id, kind))]
            if len(matched) != 1:
                raise ValueError(f'{label}: example {source_id} maps to {len(matched)} {kind} RaMP IDs; expected one')
            row = db.execute('SELECT * FROM analyte WHERE rampId=?', (matched[0],)).fetchone()
            if row is None:
                raise ValueError(f'{label}: example {source_id} references missing analyte {matched[0]}')
            stored = dict(zip(analyte_columns, row))
            examples[key] = {
                'title': title, 'match_source_id': source_id,
                'cells': {field: {'state': 'column_absent', 'value': None} if field not in stored
                          else {'state': 'null' if stored[field] is None else 'available',
                                'value': stored[field]}
                          for field in EXAMPLE_FIELDS},
            }
        lookup = {}
        for kind in ('compound', 'gene'):
            lookup[f'{kind} IDs with multiple RaMP IDs'] = db.execute('''
                SELECT count(*) FROM (SELECT sourceId FROM source WHERE geneOrCompound=?
                GROUP BY sourceId HAVING count(DISTINCT rampId)>1)''', (kind,)).fetchone()[0]
        lookup['IDs used for both entity types'] = db.execute('''
            SELECT count(*) FROM (SELECT sourceId FROM source
            GROUP BY sourceId HAVING count(DISTINCT geneOrCompound)>1)''').fetchone()[0]
        lookup['Rows missing ID, RaMP ID, type, or provider'] = db.execute('''
            SELECT count(*) FROM source WHERE sourceId IS NULL OR trim(sourceId)=''
            OR rampId IS NULL OR trim(rampId)='' OR IDtype IS NULL OR trim(IDtype)=''
            OR geneOrCompound IS NULL
            OR trim(geneOrCompound)='' OR dataSource IS NULL OR trim(dataSource)='' ''').fetchone()[0]
        lookup['Duplicate ID–RaMP ID–provider rows'] = db.execute('''
            SELECT coalesce(sum(n-1),0) FROM (SELECT count(*) n FROM source
            GROUP BY sourceId,rampId,dataSource HAVING n>1)''').fetchone()[0]
        has_analyte = db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='analyte'").fetchone()
        lookup['Rows with missing analyte'] = (db.execute('''
            SELECT count(*) FROM source s WHERE NOT EXISTS
            (SELECT 1 FROM analyte a WHERE a.rampId=s.rampId)''').fetchone()[0]
            if has_analyte else None)
        synonym_columns = {row[1] for row in db.execute('PRAGMA table_info(analytesynonym)')}
        synonym_sources = {'compound': {}, 'gene': {}}
        if synonym_columns:
            required_synonym = {'Synonym', 'rampId', 'geneOrCompound', 'source'}
            if not required_synonym <= synonym_columns:
                raise ValueError(f'{label}: analytesynonym is missing {sorted(required_synonym - synonym_columns)}')
            for provider, kind, count, ramp_ids_with_synonyms in db.execute('''
                    SELECT coalesce(nullif(trim(source),''),'(missing source)'), geneOrCompound,
                           count(*),count(DISTINCT rampId)
                    FROM analytesynonym
                    GROUP BY coalesce(nullif(trim(source),''),'(missing source)'), geneOrCompound'''):
                if kind not in synonym_sources:
                    continue
                synonym_sources[kind][provider] = {'rows': count, 'ramp_ids': ramp_ids_with_synonyms}
        ontology_present, ontology_associations = association_coverage(
            db, label, 'analytehasontology', 'ontology',
            {'rampCompoundId', 'rampOntologyId'}, '''
                SELECT 'HMDB ontology / ' || coalesce(nullif(trim(o.HMDBOntologyType),''),'(missing type)'),
                       count(*),count(DISTINCT a.rampCompoundId),count(DISTINCT a.rampOntologyId)
                FROM analytehasontology a JOIN ontology o USING(rampOntologyId)
                GROUP BY coalesce(nullif(trim(o.HMDBOntologyType),''),'(missing type)')''')
        if ontology_present:
            if db.execute('''SELECT 1 FROM analytehasontology a WHERE NOT EXISTS
                (SELECT 1 FROM analyte x WHERE x.rampId=a.rampCompoundId AND x.type='compound')
                LIMIT 1''').fetchone():
                raise ValueError(f'{label}: ontology association has a missing or non-metabolite analyte')
            rows, analytes, targets = db.execute('''
                SELECT count(*),count(DISTINCT rampCompoundId),count(DISTINCT rampOntologyId)
                FROM analytehasontology''').fetchone()
            if sum(value['rows'] for value in ontology_associations.values()) != rows:
                raise ValueError(f'{label}: ontology association has a missing ontology target')
            ontology_associations['HMDB ontology / All types'] = {
                'rows': rows, 'analytes': analytes, 'targets': targets}
        pathway_present, pathway_associations = association_coverage(
            db, label, 'analytehaspathway', 'pathway',
            {'rampId', 'pathwayRampId', 'pathwaySource'}, '''
                SELECT CASE WHEN x.type='compound' THEN 'Metabolites' ELSE 'Genes/proteins' END
                       || ' / ' || coalesce(nullif(trim(a.pathwaySource),''),'(missing source)'),
                       count(*),count(DISTINCT a.rampId),count(DISTINCT a.pathwayRampId)
                FROM analytehaspathway a JOIN analyte x ON x.rampId=a.rampId
                GROUP BY x.type,coalesce(nullif(trim(a.pathwaySource),''),'(missing source)')''')
        if pathway_present:
            if db.execute('''SELECT 1 FROM analytehaspathway a WHERE NOT EXISTS
                (SELECT 1 FROM pathway p WHERE p.pathwayRampId=a.pathwayRampId)
                LIMIT 1''').fetchone():
                raise ValueError(f'{label}: pathway association has a missing pathway target')
            for kind, title in (('compound', 'Metabolites'), ('gene', 'Genes/proteins')):
                rows, analytes, targets = db.execute('''
                    SELECT count(*),count(DISTINCT a.rampId),count(DISTINCT a.pathwayRampId)
                    FROM analytehaspathway a JOIN analyte x ON x.rampId=a.rampId
                    WHERE x.type=?''', (kind,)).fetchone()
                pathway_associations[f'{title} / All sources'] = {
                    'rows': rows, 'analytes': analytes, 'targets': targets}
            if sum(value['rows'] for key, value in pathway_associations.items()
                   if not key.endswith('/ All sources')) != sum(
                       pathway_associations[f'{title} / All sources']['rows']
                       for title in ('Metabolites', 'Genes/proteins')):
                raise ValueError(f'{label}: pathway association has an unsupported analyte type or missing analyte')
        manifest = {}
        if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='ramp_export_metadata'").fetchone():
            row = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
            if row:
                manifest = json.loads(row[0])
        table_names = {name for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")}
        scope = manifest.get('scope') or ('source_table_diagnostic' if table_names == {'analyte', 'source'}
                                          else 'lookup_table_diagnostic' if table_names in (
                                              {'analyte', 'source', 'analytesynonym'},
                                              {'analyte', 'source', 'analytesynonym', 'pathway', 'ontology',
                                               'analytehaspathway', 'analytehasontology'})
                                          else 'historical full database')
        return {'label': label, 'path': str(path), 'scope': scope,
                'ramp_ids': ramp_ids, 'examples': examples, 'lookup_quality': lookup,
                'analyte_field_coverage': {kind: values for (kind,), values in analyte_coverage.items()},
                'source_field_coverage': {
                    kind: {provider or '(missing provider)': values
                           for (provider, row_kind), values in source_field_counts.items()
                           if row_kind == kind}
                    for kind in ('compound', 'gene')},
                'synonym_table_present': bool(synonym_columns),
                'synonym_sources': synonym_sources,
                'ontology_table_present': ontology_present,
                'pathway_association_table_present': pathway_present,
                'ontology_associations': ontology_associations,
                'pathway_associations': pathway_associations,
                **source_coverage,
                'source_rows': db.execute('SELECT count(*) FROM source').fetchone()[0]}
    finally:
        db.close()


def compare(databases):
    report = {'format_version': 1, 'generated_at': datetime.now(timezone.utc).isoformat(),
            'definition': 'source records identifiers each input used to report retained RaMP data; analytesynonym records names by attributed source',
            'databases': [profile(label, path) for label, path in databases]}
    report['source_examples'] = source_examples(report['databases'])
    report['synonym_examples'] = matched_synonym_examples(report['databases'])
    report['ontology_example'] = matched_ontology_example(report['databases'])
    report['pathway_examples'] = matched_pathway_examples(report['databases'])
    return report


def _association_databases(profiles):
    dbs = [sqlite3.connect(Path(profile['path']).as_uri() + '?mode=ro', uri=True)
           for profile in profiles]
    for db in dbs:
        db.execute('PRAGMA query_only=ON')
    return dbs


def _identity_rows(db, source_id, kind):
    return [row[0] for row in db.execute('''
        SELECT DISTINCT rampId FROM source WHERE sourceId=? AND geneOrCompound=? LIMIT 2''',
        (source_id, kind))]


def matched_ontology_example(profiles):
    available = [i for i, profile in enumerate(profiles) if profile['ontology_table_present']
                 and profile['ontology_associations']['HMDB ontology / All types']['rows']]
    if not available:
        return None
    dbs = _association_databases(profiles)
    try:
        latest = available[-1]
        preferred = EXAMPLES[0][2]
        candidates = [(preferred, row[0], row[1]) for row in dbs[latest].execute('''
            SELECT DISTINCT o.HMDBOntologyType,o.commonName
            FROM source s JOIN analytehasontology a ON a.rampCompoundId=s.rampId
            JOIN ontology o ON o.rampOntologyId=a.rampOntologyId
            WHERE s.sourceId=? AND s.geneOrCompound='compound'
              AND o.HMDBOntologyType IS NOT NULL AND o.commonName IS NOT NULL
            ORDER BY o.HMDBOntologyType,o.commonName LIMIT 30''', (preferred,))]
        if not candidates:
            for rid, term_type, name in dbs[latest].execute('''
                SELECT a.rampCompoundId,o.HMDBOntologyType,o.commonName
                FROM analytehasontology a JOIN ontology o USING(rampOntologyId)
                WHERE o.HMDBOntologyType IS NOT NULL AND o.commonName IS NOT NULL
                LIMIT 30'''):
                source_ids = [row[0] for row in dbs[latest].execute('''
                    SELECT DISTINCT sourceId FROM source
                    WHERE rampId=? AND geneOrCompound='compound' AND sourceId IS NOT NULL
                    ORDER BY CASE WHEN IDtype='hmdb' THEN 0 ELSE 1 END,sourceId LIMIT 3''', (rid,))]
                candidates.extend((source_id, term_type, name) for source_id in source_ids)
        if not candidates:
            raise ValueError('Ontology associations exist but no stable example identity is available')

        def row(index, source_id, term_type, name):
            if not profiles[index]['ontology_table_present']:
                return {'state': 'table_absent'}
            rids = _identity_rows(dbs[index], source_id, 'compound')
            if len(rids) != 1:
                return {'state': 'source_id_absent' if not rids else 'ambiguous_source_id'}
            terms = [value[0] for value in dbs[index].execute('''
                SELECT rampOntologyId FROM ontology
                WHERE HMDBOntologyType=? AND commonName=? LIMIT 2''', (term_type, name))]
            if len(terms) != 1:
                return {'state': 'term_absent' if not terms else 'ambiguous_term', 'ramp_id': rids[0]}
            found = dbs[index].execute('''SELECT 1 FROM analytehasontology
                WHERE rampCompoundId=? AND rampOntologyId=? LIMIT 1''', (rids[0], terms[0])).fetchone()
            return {'state': 'available' if found else 'association_absent',
                    'ramp_id': rids[0], 'target_id': terms[0]}

        selected = max(candidates, key=lambda candidate: (
            sum(row(i, *candidate)['state'] == 'available' for i in available),
            candidate[0] == preferred))
        return {'source_id': selected[0], 'term_type': selected[1], 'term_name': selected[2],
                'rows': [row(i, *selected) for i in range(len(profiles))]}
    finally:
        for db in dbs:
            db.close()


def matched_pathway_examples(profiles):
    groups = sorted({key for profile in profiles for key in profile['pathway_associations']
                     if not key.endswith('/ All sources')})
    dbs = _association_databases(profiles)
    try:
        examples = {}
        for group in groups:
            kind = 'compound' if group.startswith('Metabolites / ') else 'gene'
            provider = group.split(' / ', 1)[1]
            available = [i for i, profile in enumerate(profiles)
                         if profile['pathway_association_table_present']
                         and profile['pathway_associations'].get(group, {}).get('rows', 0)]
            latest = available[-1]
            preferred = EXAMPLES[0 if kind == 'compound' else 1][2]
            candidates = [(preferred, row[0], row[1]) for row in dbs[latest].execute('''
                SELECT DISTINCT p.type,p.sourceId FROM source s
                JOIN analytehaspathway a ON a.rampId=s.rampId
                JOIN pathway p ON p.pathwayRampId=a.pathwayRampId
                WHERE s.sourceId=? AND s.geneOrCompound=? AND a.pathwaySource=?
                  AND p.type IS NOT NULL AND p.sourceId IS NOT NULL
                ORDER BY p.sourceId LIMIT 20''', (preferred, kind, provider))]
            if not candidates:
                for rid, pathway_type, pathway_id in dbs[latest].execute('''
                    SELECT a.rampId,p.type,p.sourceId FROM analytehaspathway a
                    JOIN pathway p ON p.pathwayRampId=a.pathwayRampId
                    WHERE a.pathwaySource=? AND a.rampId LIKE ?
                      AND p.type IS NOT NULL AND p.sourceId IS NOT NULL
                    LIMIT 20''', (provider, 'RAMP_C_%' if kind == 'compound' else 'RAMP_G_%')):
                    source_ids = [row[0] for row in dbs[latest].execute('''
                        SELECT DISTINCT sourceId FROM source WHERE rampId=? AND geneOrCompound=?
                        AND sourceId IS NOT NULL
                        ORDER BY CASE WHEN IDtype=? THEN 0 ELSE 1 END,sourceId LIMIT 3''',
                        (rid, kind, 'hmdb' if kind == 'compound' else 'uniprot'))]
                    candidates.extend((source_id, pathway_type, pathway_id)
                                      for source_id in source_ids)
            if not candidates:
                raise ValueError(f'Pathway associations exist but no example identity is available: {group}')

            def row(index, source_id, pathway_type, pathway_id):
                if not profiles[index]['pathway_association_table_present']:
                    return {'state': 'table_absent'}
                rids = _identity_rows(dbs[index], source_id, kind)
                if len(rids) != 1:
                    return {'state': 'source_id_absent' if not rids else 'ambiguous_source_id'}
                paths = list(dbs[index].execute('''
                    SELECT pathwayRampId,pathwayName FROM pathway
                    WHERE type=? AND sourceId=? LIMIT 2''', (pathway_type, pathway_id)))
                if len(paths) != 1:
                    return {'state': 'pathway_absent' if not paths else 'ambiguous_pathway',
                            'ramp_id': rids[0]}
                found = dbs[index].execute('''SELECT 1 FROM analytehaspathway
                    WHERE rampId=? AND pathwayRampId=? AND pathwaySource=? LIMIT 1''',
                    (rids[0], paths[0][0], provider)).fetchone()
                return {'state': 'available' if found else 'association_absent',
                        'ramp_id': rids[0], 'target_id': paths[0][0], 'target_name': paths[0][1],
                        'pathway_source': provider}

            selected = max(candidates, key=lambda candidate: (
                sum(row(i, *candidate)['state'] == 'available' for i in available),
                candidate[0] == preferred))
            examples[group] = {'source_id': selected[0], 'pathway_type': selected[1],
                               'pathway_source_id': selected[2],
                               'rows': [row(i, *selected) for i in range(len(profiles))]}
        return examples
    finally:
        for db in dbs:
            db.close()


def matched_synonym_examples(profiles):
    """Compare source-attributed names for one shared source ID per group."""
    groups = sorted({(kind, provider) for profile in profiles
                     for kind, entries in profile['synonym_sources'].items()
                     for provider in entries})
    dbs = [sqlite3.connect(Path(profile['path']).as_uri() + '?mode=ro', uri=True)
           for profile in profiles]
    try:
        for db in dbs:
            db.execute('PRAGMA query_only=ON')
        matches = {}
        def ramp_ids(index, kind, source_id):
            key = (index, kind, source_id)
            if key not in matches:
                matches[key] = [row[0] for row in dbs[index].execute('''
                    SELECT DISTINCT rampId FROM source
                    WHERE sourceId=? AND geneOrCompound=? LIMIT 2''', (source_id, kind))]
            return matches[key]
        examples = {}
        for kind, provider in groups:
            latest = max(i for i, profile in enumerate(profiles)
                         if provider in profile['synonym_sources'][kind])
            candidate_ramp_ids = [row[0] for row in dbs[latest].execute('''
                SELECT DISTINCT rampId FROM analytesynonym
                WHERE source=? AND geneOrCompound=? AND rampId IS NOT NULL LIMIT 20''',
                (provider, kind))]
            candidates = set()
            for ramp_id in candidate_ramp_ids:
                candidates.update(row[0] for row in dbs[latest].execute('''
                    SELECT DISTINCT sourceId FROM source
                    WHERE rampId=? AND geneOrCompound=? AND sourceId IS NOT NULL
                    ORDER BY CASE WHEN IDtype=? THEN 0 ELSE 1 END, sourceId LIMIT 6''',
                    (ramp_id, kind, 'hmdb' if kind == 'compound' else 'uniprot')))
            preferred_id = EXAMPLES[0 if kind == 'compound' else 1][2]
            preferred_ramp_ids = ramp_ids(latest, kind, preferred_id)
            if len(preferred_ramp_ids) == 1 and dbs[latest].execute('''
                    SELECT 1 FROM analytesynonym
                    WHERE rampId=? AND geneOrCompound=? AND source=? LIMIT 1''',
                    (preferred_ramp_ids[0], kind, provider)).fetchone():
                candidates.add(preferred_id)
            if not candidates:
                raise ValueError(f'No source ID can anchor {kind} / {provider} synonym examples')
            def score(source_id):
                mapped = [ramp_ids(i, kind, source_id) for i in range(len(dbs))]
                return (sum(len(ids) == 1 for ids in mapped),
                        -sum(len(ids) > 1 for ids in mapped),
                        source_id == preferred_id,
                        source_id.startswith('hmdb:' if kind == 'compound' else 'uniprot:'))
            source_id = max(sorted(candidates, reverse=True), key=score)
            rows = []
            for i, profile in enumerate(profiles):
                ids = ramp_ids(i, kind, source_id)
                if not ids:
                    rows.append({'state': 'source_id_absent'})
                elif len(ids) > 1:
                    rows.append({'state': 'ambiguous_source_id', 'matching_ids': len(ids)})
                elif not profile['synonym_table_present']:
                    rows.append({'state': 'table_absent', 'ramp_id': ids[0]})
                else:
                    names = [row[0] for row in dbs[i].execute('''
                        SELECT Synonym FROM analytesynonym
                        WHERE rampId=? AND geneOrCompound=? AND source=?''',
                        (ids[0], kind, provider))]
                    names.sort(key=lambda value: (value is None, (value or '').casefold(), value or ''))
                    rows.append({'state': 'available' if names else 'no_attributed_synonyms',
                                 'ramp_id': ids[0], 'count': len(names), 'synonyms': names[:3]})
            examples[(kind, provider)] = {'source_id': source_id, 'rows': rows}
        return {kind: {provider: example for (row_kind, provider), example in examples.items()
                       if row_kind == kind} for kind in ('compound', 'gene')}
    finally:
        for db in dbs:
            db.close()


def source_examples(profiles):
    """Choose one uniquely mapped ID per provider and analyte type."""
    providers = sorted({provider for profile in profiles
                        for section in ('metabolite_sources', 'gene_sources')
                        for provider in profile[section] if provider != '(missing provider)'})
    dbs = [sqlite3.connect(Path(profile['path']).as_uri() + '?mode=ro', uri=True)
           for profile in profiles]
    try:
        for db in dbs:
            db.row_factory = sqlite3.Row
            db.execute('PRAGMA query_only=ON')
        examples = {}
        for provider in providers:
            examples[provider] = {}
            for kind, section in (('compound', 'metabolite_sources'), ('gene', 'gene_sources')):
                available_indices = [i for i, profile in enumerate(profiles)
                                     if profile[section].get(provider, 0) > 0]
                if not available_indices:
                    examples[provider][kind] = None
                    continue
                latest = available_indices[-1]
                available = len(available_indices)
                candidates = dbs[latest].execute('''
                SELECT DISTINCT sourceId FROM source WHERE dataSource=? AND geneOrCompound=?
                AND sourceId IS NOT NULL ORDER BY sourceId LIMIT 1000''', (provider, kind)).fetchall()
                if not candidates:
                    # Preserve trimmed-provider grouping for atypical spelling.
                    candidates = dbs[latest].execute('''
                    SELECT DISTINCT sourceId FROM source WHERE trim(dataSource)=?
                    AND geneOrCompound=? AND sourceId IS NOT NULL
                    ORDER BY sourceId LIMIT 1000''', (provider, kind)).fetchall()
                selected = None
                first_candidate = None
                best_score = -1
                for candidate_number, (source_id,) in enumerate(candidates, 1):
                    rows = [db.execute('''SELECT * FROM source WHERE sourceId=? AND trim(dataSource)=?
                         AND geneOrCompound=? ORDER BY rampId,sourceId''',
                                   (source_id, provider, kind)).fetchall()
                            for db in dbs]
                    if first_candidate is None:
                        first_candidate = (source_id, rows)
                    if any(len(found) > 1 for found in rows):
                        continue
                    score = sum(bool(found) for found in rows)
                    if score > best_score:
                        selected, best_score = (source_id, rows), score
                    if score == available:
                        break
                    if candidate_number >= 250 and selected is not None:
                        break
                if selected is None:
                    selected = first_candidate
                if selected is None:
                    raise ValueError(f'No source example for {provider} / {kind}')
                source_id, rows = selected
                compared = []
                for profile, found in zip(profiles, rows):
                    if found:
                        compared.append({'state': 'multiple_mappings' if len(found) > 1 else 'available',
                                     'matching_rows': len(found),
                                     'cells': {field: found[0][field] if field in found[0].keys()
                                               else {'state': 'column_absent'}
                                               for field in SOURCE_FIELDS}})
                    elif not profile[section].get(provider):
                        state = ('provider_absent' if not (profile['metabolite_sources'].get(provider)
                                                       or profile['gene_sources'].get(provider))
                             else 'entity_type_absent')
                        compared.append({'state': state, 'cells': None})
                    else:
                        compared.append({'state': 'example_id_absent', 'cells': None})
                examples[provider][kind] = {'source_id': source_id, 'entity_type': kind,
                                            'selected_from': profiles[latest]['label'], 'rows': compared}
        return examples
    finally:
        for db in dbs:
            db.close()


def cell(value, previous):
    if value is None:
        return {'value': None, 'delta': None, 'percent_change': None}
    delta = value - previous if previous is not None else None
    percent = 100 * delta / previous if delta is not None and previous else None
    return {'value': value, 'delta': delta, 'percent_change': percent}


def tables(report):
    profiles = report['databases']
    result = {}
    for section, fields in [('ramp_ids', (None,)), ('metabolite_sources', (None,)),
                            ('gene_sources', (None,)),
                            ('compound_synonyms', ('rows', 'ramp_ids')),
                            ('gene_synonyms', ('rows', 'ramp_ids')),
                            ('ontology_associations', ('rows', 'analytes', 'targets')),
                            ('pathway_associations', ('rows', 'analytes', 'targets')),
                            ('lookup_quality', (None,))]:
        synonym_kind = {'compound_synonyms': 'compound', 'gene_synonyms': 'gene'}.get(section)
        association_table = {'ontology_associations': 'ontology_table_present',
                             'pathway_associations': 'pathway_association_table_present'}.get(section)
        keys = sorted({key for p in profiles for key in (p['synonym_sources'][synonym_kind]
                       if synonym_kind else p[section])})
        if section == 'ramp_ids':
            keys = ['Metabolite RaMP IDs', 'Gene/protein RaMP IDs']
        if association_table:
            keys.sort(key=lambda key: (not key.endswith(('All types', 'All sources')), key))
        result[section] = []
        for key in keys:
            for field in fields:
                values = [(None if not p['synonym_table_present'] else
                           p['synonym_sources'][synonym_kind].get(key, {}).get(field, 0))
                          if synonym_kind else
                          (None if not p[association_table] else
                           p[section].get(key, {}).get(field, 0)) if association_table else
                          p[section].get(key, 0)
                          for p in profiles]
                metric_field = ('RaMP IDs' if field in ('ramp_ids', 'analytes') else
                                'ontology terms' if section == 'ontology_associations' and field == 'targets' else
                                'pathways' if section == 'pathway_associations' and field == 'targets' else field)
                result[section].append({'metric': f'{key} / {metric_field}' if field else key,
                                        'cells': [cell(value, values[i-1] if i else None)
                                                  for i, value in enumerate(values)]})
    return result


def render(report):
    esc = lambda value: html.escape(str(value))
    profiles = report['databases']
    sections = tables(report)
    names = [('ramp_ids', "RaMP IDs represented in 'analyte'"),
             ('metabolite_sources', 'Metabolites by input'),
             ('gene_sources', 'Genes/proteins by input'),
             ('source_examples', 'Example source rows'),
             ('lookup_quality', 'Lookup quality'),
             ('compound_synonyms', 'Metabolite synonyms by attributed source'),
             ('gene_synonyms', 'Gene/protein synonyms by attributed source'),
             ('synonym_examples', 'Example synonym rows'),
             ('ontology_associations', 'Metabolite ontology associations by HMDB type'),
             ('ontology_examples', 'Example ontology association'),
             ('pathway_associations', 'Pathway associations by source and analyte type'),
             ('pathway_examples', 'Example pathway associations')]
    panel_for = {'ramp_ids': 'analyte',
                 'metabolite_sources': 'source', 'gene_sources': 'source',
                 'source_examples': 'source', 'lookup_quality': 'source',
                 'compound_synonyms': 'analytesynonym', 'gene_synonyms': 'analytesynonym',
                 'synonym_examples': 'analytesynonym',
                 'ontology_associations': 'analytehasontology',
                 'ontology_examples': 'analytehasontology',
                 'pathway_associations': 'analytehaspathway',
                 'pathway_examples': 'analytehaspathway'}
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
             '<title>RaMP lookup and association comparison</title>',
             '<style>body{font:16px/1.5 system-ui,sans-serif;margin:0;background:#f5f7fa;color:#182635}main{max-width:1500px;margin:auto;padding:24px}p{max-width:100ch}.scroll{overflow-x:auto;background:white;border:1px solid #d2dbe3}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px 12px;border-bottom:1px solid #dde4eb;text-align:right;min-width:140px}th:first-child,td:first-child{text-align:left;min-width:300px}th{background:#eaf0f6}tbody tr:nth-child(even){background:#f8fafc}tbody tr.coverage{background:#eef3f7;color:#42566a;font-size:12px}tbody tr.coverage td{padding-top:4px;padding-bottom:8px}small{display:block;color:#35475a}article{background:white;border-left:4px solid #52738c;padding:10px 16px;margin:10px 0}.tabs{display:flex;flex-wrap:wrap;gap:6px;border-bottom:2px solid #c5d2de;margin:24px 0 18px}.tabs button{font:inherit;color:#23527c;background:#eaf0f6;border:1px solid #c5d2de;border-bottom:0;border-radius:6px 6px 0 0;padding:9px 18px;cursor:pointer}.tabs button[aria-selected=true]{background:white;color:#182635;font-weight:700;border-color:#52738c}.tabs button:focus-visible{outline:3px solid #23527c;outline-offset:2px}.tab-panel[hidden]{display:none}.example th,.example td{text-align:left}.example th:first-child,.example td:first-child{min-width:220px}.synonym-example td{overflow-wrap:anywhere;white-space:normal;max-width:650px}@media(max-width:600px){main{padding:12px}.tabs button{flex:1 1 auto}}@media print{.tabs{display:none}.tab-panel[hidden]{display:block}}</style><main>',
             '<h1>RaMP lookup and association comparison</h1>']
    def coverage_td(current, previous, field):
        current_percent = coverage_percent(current, field)
        previous_percent = coverage_percent(previous, field)
        points = (current_percent - previous_percent
                  if current_percent is not None and previous_percent is not None else None)
        result = f'<td{change_style(points, points)}>{esc(coverage_cell(current, field))}'
        if points is not None:
            result += f'<small>Δ {points:+.1f} pp</small>'
        return result + '</td>'
    tab_labels = [('analyte', 'Analyte'), ('source', 'Source'),
                  ('analytesynonym', 'Synonyms'), ('analytehasontology', 'Ontology'),
                  ('analytehaspathway', 'Pathways'), ('about', 'About this report')]
    parts.append('<div class="tabs" role="tablist" aria-label="SQLite tables">' +
                 ''.join(f'<button type="button" role="tab" id="tab-{panel}" aria-controls="{panel}" aria-selected="{str(index == 0).lower()}" tabindex="{0 if index == 0 else -1}">{label}</button>'
                         for index, (panel, label) in enumerate(tab_labels)) + '</div>')
    current_panel = None
    for i, (key, title) in enumerate(names):
        panel = panel_for[key]
        if panel != current_panel:
            if current_panel is not None:
                parts.append('</section>')
            parts.append(f'<section class="tab-panel" id="{panel}" role="tabpanel" aria-labelledby="tab-{panel}">')
            current_panel = panel
        parts.append(f'<h2 id="s{i}">{esc(title)}</h2>')
        if key == 'ramp_ids':
            parts.append('<p>Metabolites and genes/proteins are counted separately. These are within-build analyte counts, not stable IDs across releases.</p>')
        if key in ('metabolite_sources', 'gene_sources'):
            parts.append('<p>Each cell counts distinct RaMP IDs with at least one ID attributed to that input in <code>source</code>. An analyte may occur under several inputs, so the rows do not add up to the analyte total. KEGG IDs supplied by HMDB or WikiPathways retain their separate legacy source labels. Changes describe coverage, not quality.</p>')
        if key == 'source_examples':
            parts.append('<p>Each input has separate metabolite and gene/protein examples when those rows exist. Each example follows one stored source ID across builds. Field coverage counts populated values across all rows for that input and analyte type; NULL, blanks, and known placeholders are unpopulated.</p>')
            for provider, by_kind in report['source_examples'].items():
                parts.append(f'<details><summary>{esc(provider)}</summary>')
                for kind, label in (('compound', 'Metabolite'), ('gene', 'Gene/protein')):
                    example = by_kind[kind]
                    if example is None:
                        parts.append(f'<p><strong>{label}:</strong> No {label.lower()} rows for this input in any displayed build.</p>')
                        continue
                    parts.append(f'<details><summary>{label} · <code>{esc(example["source_id"])}</code></summary>')
                    parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in SOURCE_FIELDS) + '</tr></thead><tbody>')
                    for profile_index, (profile, row) in enumerate(zip(profiles, example['rows'])):
                        parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                        if row['state'] in ('available', 'multiple_mappings'):
                            for field in SOURCE_FIELDS:
                                value = row['cells'][field]
                                display = ('Column absent' if isinstance(value, dict)
                                       and value.get('state') == 'column_absent' else
                                       esc('NULL' if value is None else value))
                                if field == 'rampId' and row['state'] == 'multiple_mappings':
                                    display += f'<small>Showing 1 of {row["matching_rows"]} mappings</small>'
                                parts.append(f'<td>{display}</td>')
                        else:
                            status = {'provider_absent': 'Input absent',
                                  'entity_type_absent': 'No rows for this analyte type',
                                  'example_id_absent': 'Example ID absent'}[row['state']]
                            parts.append(f'<td colspan="{len(SOURCE_FIELDS)}">{esc(status)}</td>')
                        parts.append('</tr>')
                        coverage = profile['source_field_coverage'].get(kind, {}).get(provider)
                        previous = (profiles[profile_index - 1]['source_field_coverage']
                                    .get(kind, {}).get(provider) if profile_index else None)
                        parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                        parts.extend(coverage_td(coverage, previous, field) for field in SOURCE_FIELDS)
                        parts.append('</tr>')
                    parts.append('</tbody></table></div></details>')
                parts.append('</details>')
            continue
        if key in ('compound_synonyms', 'gene_synonyms'):
            parts.append('<p>Rows count stored source-attributed names; RaMP IDs count analytes with at least one such name. A name can occur under multiple sources. Old Rhea and Reactome gene names include UniProt enrichment; movement to <code>uniprot</code> can reflect corrected attribution. Old WikiPathways gene symbols came from WikiPathways itself. Earlier new builds omitted those source-owned symbols; the current exporter reads them from the patched graph.</p>')
        if key == 'synonym_examples':
            parts.append('<p>Each input has separate metabolite and gene/protein examples where attributed synonyms exist. Each example follows the same source ID across builds. The RaMP ID is resolved separately in each database; up to three names are shown per build.</p>')
            providers = sorted({provider for entries in report['synonym_examples'].values()
                                for provider in entries})
            for provider in providers:
                parts.append(f'<details><summary>{esc(provider)}</summary>')
                for kind, label in (('compound', 'Metabolite'), ('gene', 'Gene/protein')):
                    example = report['synonym_examples'][kind].get(provider)
                    if example is None:
                        parts.append(f'<p><strong>{label}:</strong> No {label.lower()} synonyms attributed to this input in any displayed build.</p>')
                        continue
                    parts.append(f'<details><summary>{label} · <code>{esc(example["source_id"])}</code></summary>')
                    parts.append(f'<p>Identity anchor: <code>{esc(example["source_id"])}</code>; attributed synonym source: <code>{esc(provider)}</code>.</p>')
                    parts.append('<div class="scroll"><table class="example synonym-example"><thead><tr><th>Database</th><th>RaMP ID</th><th>Synonym count</th><th>Synonyms</th></tr></thead><tbody>')
                    for profile, row in zip(profiles, example['rows']):
                        parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                        if row['state'] in ('available', 'no_attributed_synonyms'):
                            parts.append(f'<td>{esc(row["ramp_id"])}</td><td>{row["count"]:,}</td>')
                            parts.append('<td>' + ('<br>'.join(esc('NULL' if name is None else name)
                                               for name in row['synonyms'])
                                               if row['synonyms'] else 'No synonyms attributed to this source') + '</td>')
                        else:
                            status = {'source_id_absent': 'Source ID absent',
                                  'ambiguous_source_id': 'Source ID maps to multiple RaMP IDs',
                                  'table_absent': 'Synonym table absent'}[row['state']]
                            parts.append(f'<td colspan="3">{esc(status)}</td>')
                        parts.append('</tr>')
                    parts.append('</tbody></table></div></details>')
                parts.append('</details>')
            continue
        if key == 'ontology_associations':
            parts.append('<p>Rows count analyte–ontology links; RaMP IDs and ontology terms are distinct within each HMDB ontology type. This table has no source column, so HMDB attribution is inferred from its HMDB ontology terms. Types can overlap on an analyte.</p>')
        if key == 'pathway_associations':
            parts.append('<p>Rows count analyte–pathway links by the stored <code>pathwaySource</code>. RaMP IDs and pathways are distinct within each source; an analyte or pathway may occur under several sources.</p>')
        if key in ('ontology_examples', 'pathway_examples'):
            if key == 'ontology_examples':
                groups = [('HMDB ontology', [('Metabolite', report['ontology_example']),
                                            ('Gene/protein', None)])]
            else:
                providers = sorted({group.split(' / ', 1)[1] for group in report['pathway_examples']})
                groups = [(provider, [
                    ('Metabolite', report['pathway_examples'].get(f'Metabolites / {provider}')),
                    ('Gene/protein', report['pathway_examples'].get(f'Genes/proteins / {provider}'))])
                    for provider in providers]
            if not groups:
                parts.append('<p>No associations are available for an example.</p>')
            for provider, pairs in groups:
                parts.append(f'<details><summary>{esc(provider)}</summary>')
                for label, example in pairs:
                    if example is None:
                        status = ('Not applicable: ontology links only metabolites.'
                                  if key == 'ontology_examples' else
                                  f'No {label.lower()} pathway associations for this source in any displayed build.')
                        parts.append(f'<p><strong>{label}:</strong> {esc(status)}</p>')
                        continue
                    target = (f'{example["term_type"]} / {example["term_name"]}'
                              if key == 'ontology_examples' else
                              f'{example["pathway_type"]} / {example["pathway_source_id"]}')
                    parts.append(f'<details><summary>{label} · <code>{esc(example["source_id"])}</code> → {esc(target)}</summary>')
                    parts.append('<p>Matched across builds by the analyte source ID and the target’s source identity. RaMP IDs are shown from each build.</p>')
                    parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th><th>Analyte RaMP ID</th><th>Target RaMP ID</th><th>Target name</th><th>Association</th></tr></thead><tbody>')
                    for profile, row in zip(profiles, example['rows']):
                        status = {'available': 'Present', 'association_absent': 'Absent',
                              'source_id_absent': 'Analyte source ID absent',
                              'ambiguous_source_id': 'Analyte source ID maps to multiple RaMP IDs',
                              'term_absent': 'Ontology term absent',
                              'ambiguous_term': 'Ontology term is ambiguous',
                              'pathway_absent': 'Pathway absent',
                              'ambiguous_pathway': 'Pathway identity is ambiguous',
                              'table_absent': 'Table absent'}[row['state']]
                        target_name = example['term_name'] if key == 'ontology_examples' else row.get('target_name')
                        parts.append(f'<tr><td>{esc(profile["label"])}</td><td>{esc(row.get("ramp_id", "—"))}</td><td>{esc(row.get("target_id", "—"))}</td><td>{esc(target_name or "—")}</td><td>{esc(status)}</td></tr>')
                    parts.append('</tbody></table></div></details>')
                parts.append('</details>')
            continue
        if key == 'lookup_quality':
            parts.append('<p>A reported ID can map to several gene/protein RaMP IDs. These are candidate lookups, not guaranteed unique resolutions.</p>')
        parts.append('<div class="scroll"><table><thead><tr><th>Metric</th>' + ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
        for row in sections[key]:
            parts.append(f'<tr><td>{esc(row["metric"])}</td>')
            for profile, datum in zip(profiles, row['cells']):
                value, delta, percent = datum['value'], datum['delta'], datum['percent_change']
                style = change_style(delta, percent if percent is not None and math.isfinite(percent) else 100)
                missing = ('Table absent' if (
                    key in ('compound_synonyms', 'gene_synonyms') and not profile['synonym_table_present']
                    or key == 'ontology_associations' and not profile['ontology_table_present']
                    or key == 'pathway_associations' and not profile['pathway_association_table_present'])
                           else 'Unavailable')
                parts.append(f'<td{style}>{value:,}' if value is not None else f'<td{style}>{missing}')
                if delta is not None:
                    parts.append(f'<small>Δ {delta:+,} ({percent:+.1f}%)</small>' if percent is not None
                                 else f'<small>Δ {delta:+,} ({"from zero" if delta > 0 else "—"})</small>')
                parts.append('</td>')
            parts.append('</tr>')
        parts.append('</tbody></table></div>')
        if key == 'ramp_ids':
            for example_key, title, source_id, _ in EXAMPLES:
                parts.append(f'<h3>{esc(title)}</h3><p>Matched by <code>{esc(source_id)}</code> in each database’s <code>source</code> table. Each example comes from one <code>analyte</code> row; field coverage counts populated values across all analytes of that type. NULL, blanks, and known placeholders are unpopulated.</p>')
                parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in EXAMPLE_FIELDS) + '</tr></thead><tbody>')
                for profile_index, profile in enumerate(profiles):
                    parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                    for field in EXAMPLE_FIELDS:
                        datum = profile['examples'][example_key]['cells'][field]
                        display = ('Column absent' if datum['state'] == 'column_absent' else
                                   'NULL' if datum['state'] == 'null' else str(datum['value']))
                        parts.append(f'<td>{esc(display)}</td>')
                    parts.append('</tr>')
                    coverage = profile['analyte_field_coverage'].get('compound' if example_key == 'metabolite' else 'gene')
                    previous = (profiles[profile_index - 1]['analyte_field_coverage']
                                .get('compound' if example_key == 'metabolite' else 'gene')
                                if profile_index else None)
                    parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                    parts.extend(coverage_td(coverage, previous, field) for field in EXAMPLE_FIELDS)
                    parts.append('</tr>')
                parts.append('</tbody></table></div>')
    parts.append('</section>')
    parts.append('<section class="tab-panel" id="about" role="tabpanel" aria-labelledby="tab-about">'
                 '<h2>About this report</h2>'
                 '<p>The analyte table counts distinct RaMP IDs. Its examples match records across builds by source ID, since RaMP IDs can change.</p>'
                 '<p>Each change compares with the preceding build. Fresh RaMP IDs are counted within a build and are never matched across builds.</p>'
                 '<p>Green = increased; red = decreased versus the preceding build. Color shows direction, not quality. Field coverage compares percentage populated.</p>'
                 f'<p>Generated {esc(report["generated_at"])}.</p>')
    for p in profiles:
        parts.append(f'<article><strong>{esc(p["label"])}</strong> · {p["source_rows"]:,} source rows · {esc(p["scope"])}<br><code>{esc(p["path"])}</code></article>')
    parts.append('</section>')
    parts.append("""</main><script>
(() => {
  const tabs = Array.from(document.querySelectorAll('[role="tab"]'));
  const panels = new Map(Array.from(document.querySelectorAll('[role="tabpanel"]'),
    panel => [panel.id, panel]));
  function selectedFromHash() {
    const id = location.hash.slice(1);
    return panels.has(id) ? id : 'analyte';
  }
  function activate(id, updateHash = false, focus = false) {
    if (!panels.has(id)) id = 'analyte';
    for (const tab of tabs) {
      const selected = tab.getAttribute('aria-controls') === id;
      tab.setAttribute('aria-selected', String(selected));
      tab.tabIndex = selected ? 0 : -1;
      if (selected && focus) tab.focus();
    }
    for (const [panelId, panel] of panels) panel.hidden = panelId !== id;
    if (updateHash && location.hash !== '#' + id) history.pushState(null, '', '#' + id);
  }
  tabs.forEach((tab, index) => {
    tab.addEventListener('click', () => activate(tab.getAttribute('aria-controls'), true));
    tab.addEventListener('keydown', event => {
      let next;
      if (event.key === 'ArrowRight' || event.key === 'ArrowDown') next = (index + 1) % tabs.length;
      else if (event.key === 'ArrowLeft' || event.key === 'ArrowUp') next = (index - 1 + tabs.length) % tabs.length;
      else if (event.key === 'Home') next = 0;
      else if (event.key === 'End') next = tabs.length - 1;
      else return;
      event.preventDefault();
      activate(tabs[next].getAttribute('aria-controls'), true, true);
    });
  });
  window.addEventListener('hashchange', () => activate(selectedFromHash()));
  window.addEventListener('popstate', () => activate(selectedFromHash()));
  activate(selectedFromHash());
})();
</script></html>""")
    return '\n'.join(parts)


def write(report, output):
    output = Path(output).resolve()
    if output.suffix.lower() != '.html':
        raise ValueError('Report output must end in .html')
    targets = (output, output.with_suffix('.json'))
    inputs = [Path(p['path']) for p in report['databases']]
    if any(target in inputs or any(target.exists() and target.samefile(p) for p in inputs) for target in targets):
        raise ValueError('Report output must not overwrite an input database')
    output.parent.mkdir(parents=True, exist_ok=True)
    report = {**report, 'comparison_tables': tables(report)}
    for target, content in zip(targets, (render(report), json.dumps(report, indent=2))):
        fd, temporary = tempfile.mkstemp(dir=output.parent, suffix='.tmp')
        try:
            with os.fdopen(fd, 'w') as handle:
                handle.write(content)
            os.replace(temporary, target)
        finally:
            Path(temporary).unlink(missing_ok=True)
    return targets


def main(argv=None):
    parser = argparse.ArgumentParser(description='Compare RaMP lookup and association tables across SQLite builds')
    parser.add_argument('--database', action='append', required=True, metavar='LABEL=PATH')
    parser.add_argument('--output', required=True, type=Path)
    args = parser.parse_args(argv)
    databases = []
    for item in args.database:
        label, sep, path = item.partition('=')
        if not sep or not label.strip() or not path:
            parser.error('Each --database must be LABEL=PATH')
        databases.append((label.strip(), path))
    if len(databases) < 2 or len({label for label, _ in databases}) != len(databases):
        parser.error('Provide at least two databases with unique labels')
    try:
        report = compare(databases)
        for path in write(report, args.output):
            print(path)
    except (OSError, sqlite3.Error, ValueError) as exc:
        parser.exit(1, f'Source comparison failed: {exc}\n')


if __name__ == '__main__':
    main()
