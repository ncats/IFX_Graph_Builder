"""Read-only comparison of RaMP lookup, association, chemistry, and version tables."""
import argparse
from datetime import datetime, timezone
import html
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile
import zlib

from src.use_cases.ramp.pathway_similarity import SCOPES as PATHWAY_SCOPES


EXAMPLES = (
    ('metabolite', 'Metabolite example: D-glucose', 'hmdb:HMDB0000122', 'compound'),
    ('gene', 'Gene/protein example: EGFR', 'uniprot:P00533', 'gene'),
)
EXAMPLE_FIELDS = ('rampId', 'type', 'common_name')
SOURCE_FIELDS = ('sourceId', 'rampId', 'IDtype', 'geneOrCompound', 'commonName',
                 'priorityHMDBStatus', 'dataSource', 'pathwayCount')
CATALYZED_FIELDS = ('rampCompoundId', 'rampGeneId', 'proteinType')
CLASS_FIELDS = ('ramp_id', 'class_source_id', 'class_level_name', 'class_name', 'source')
REACTION_FIELDS = ('ramp_rxn_id', 'rxn_source_id', 'status', 'is_transport', 'direction',
                   'label', 'equation', 'html_equation', 'ec_num', 'has_human_prot', 'only_human_mets')
REACTION_MET_FIELDS = ('ramp_rxn_id', 'rxn_source_id', 'ramp_cmpd_id', 'substrate_product',
                       'met_source_id', 'met_name', 'is_cofactor')
REACTION_PROTEIN_FIELDS = ('ramp_rxn_id', 'rxn_source_id', 'ramp_gene_id', 'uniprot',
                           'protein_name', 'is_reviewed')
REACTION_EC_FIELDS = ('ramp_rxn_id', 'rxn_source_id', 'rxn_class_ec', 'ec_level',
                      'rxn_class', 'rxn_class_hierarchy')
CHEM_FIELDS = ('ramp_id', 'chem_data_source', 'chem_source_id', 'iso_smiles',
               'inchi_key_prefix', 'inchi_key', 'inchi', 'mw', 'monoisotop_mass',
               'common_name', 'mol_formula')
DB_VERSION_FIELDS = ('met_intersects_json', 'gene_intersects_json',
                     'met_intersects_json_pw_mapped', 'gene_intersects_json_pw_mapped')
DB_VERSION_SCOPE_LABELS = {
    'met_intersects_json': 'Metabolites / all source-bearing analytes',
    'gene_intersects_json': 'Genes/proteins / all source-bearing analytes',
    'met_intersects_json_pw_mapped': 'Metabolites / mapped to non-SMPDB pathway',
    'gene_intersects_json_pw_mapped': 'Genes/proteins / mapped to non-SMPDB pathway',
}
PATHWAY_SIMILARITY_EXAMPLE = 'R-HSA-109581'
PATHWAY_DUPLICATE_EXAMPLE = ('R-HSA-110329', 'R-HSA-73928')
PATHWAY_OVERLAP_SCOPES = tuple((blob, cutoff, condition)
                               for blob in ('analyte_blob', 'metabolite_blob', 'gene_blob')
                               for _, scope_blob, _, cutoff, condition in PATHWAY_SCOPES
                               if scope_blob == blob)
REACTION_TABLES = {
    'reaction': REACTION_FIELDS,
    'reaction2met': REACTION_MET_FIELDS,
    'reaction2protein': REACTION_PROTEIN_FIELDS,
    'reaction_ec_class': REACTION_EC_FIELDS,
}
REACTION_ANCHORS = {
    'reaction': ('rxn_source_id',),
    'reaction2met': ('rxn_source_id', 'met_source_id', 'substrate_product'),
    'reaction2protein': ('rxn_source_id', 'uniprot'),
    'reaction_ec_class': ('rxn_source_id', 'rxn_class_ec', 'ec_level'),
}


def populated_expression(field):
    column = '"' + field.replace('"', '""') + '"'
    value = f'lower(trim(CAST({column} AS TEXT)))'
    excluded = ("'na','n/a'" if field in ('commonName', 'common_name') else
                "'no_hmdb_status'" if field == 'priorityHMDBStatus' else
                "'unknown','na','n/a'" if field == 'proteinType' else
                "'-1'" if field in ('pathwayCount', 'is_reviewed', 'is_cofactor',
                                      'has_human_prot', 'only_human_mets') else None)
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


def reaction_profiles(db, label):
    """Profile reaction tables without assuming older releases have every column."""
    result = {}
    for table, fields in REACTION_TABLES.items():
        columns = {row[1] for row in db.execute(f'PRAGMA table_info({table})')}
        required = set(REACTION_ANCHORS[table]) | {'ramp_rxn_id'}
        if table == 'reaction2met':
            required.add('ramp_cmpd_id')
        if table == 'reaction2protein':
            required.add('ramp_gene_id')
        if columns and not required <= columns:
            raise ValueError(f'{label}: {table} is missing {sorted(required - columns)}')
        counts = {}
        coverage = None
        if columns:
            participants = (',count(DISTINCT ramp_cmpd_id)' if table == 'reaction2met' else
                            ',count(DISTINCT ramp_gene_id)' if table == 'reaction2protein' else '')
            def record(key, where='', params=()):
                values = db.execute(f'''SELECT count(*),count(DISTINCT rxn_source_id){participants}
                    FROM {table} {where}''', params).fetchone()
                counts[key] = {'rows': values[0], 'reactions': values[1]}
                if table in ('reaction2met', 'reaction2protein'):
                    counts[key]['analytes'] = values[2]
            record('All rows')
            if table != 'reaction' and db.execute(
                    "SELECT 1 FROM sqlite_master WHERE type='table' AND name='reaction'").fetchone():
                record('Missing reaction target', '''WHERE NOT EXISTS
                    (SELECT 1 FROM reaction r WHERE r.ramp_rxn_id=''' + table + '''.ramp_rxn_id)''')
            if table == 'reaction':
                for field, title in (('status', 'Status'), ('direction', 'Direction')):
                    for (value,) in db.execute(f'SELECT DISTINCT {field} FROM reaction ORDER BY {field}'):
                        record(f'{title} / {value if value is not None else "NULL"}',
                               f'WHERE {field} IS ?', (value,))
            elif table == 'reaction2met':
                for (side,) in db.execute('SELECT DISTINCT substrate_product FROM reaction2met ORDER BY 1'):
                    record(f'Side / {side if side is not None else "NULL"}',
                           'WHERE substrate_product IS ?', (side,))
                if 'is_cofactor' in columns:
                    for (flag,) in db.execute('SELECT DISTINCT is_cofactor FROM reaction2met ORDER BY 1'):
                        record(f'Cofactor / {flag if flag is not None else "NULL"}',
                               'WHERE is_cofactor IS ?', (flag,))
            elif table == 'reaction2protein' and 'is_reviewed' in columns:
                for (flag,) in db.execute('SELECT DISTINCT is_reviewed FROM reaction2protein ORDER BY 1'):
                    record(f'Reviewed / {flag if flag is not None else "NULL"}',
                           'WHERE is_reviewed IS ?', (flag,))
            elif table == 'reaction_ec_class':
                for (level,) in db.execute('SELECT DISTINCT ec_level FROM reaction_ec_class ORDER BY 1'):
                    record(f'EC level / {level if level is not None else "NULL"}',
                           'WHERE ec_level IS ?', (level,))
            coverage = field_coverage(db, table, fields, columns, ('1',)).get((1,))
        result[table] = {'present': bool(columns), 'columns': sorted(columns),
                         'counts': counts, 'field_coverage': coverage}
    return result


def chemistry_profile(db, label):
    columns = {row[1] for row in db.execute('PRAGMA table_info(chem_props)')}
    if not columns:
        return {'present': False, 'counts': {}, 'field_coverage': {}}
    required = {'ramp_id', 'chem_data_source', 'chem_source_id'}
    if not required <= columns:
        raise ValueError(f'{label}: chem_props is missing {sorted(required - columns)}')
    if db.execute('''SELECT 1 FROM chem_props c WHERE NOT EXISTS
        (SELECT 1 FROM analyte a WHERE a.rampId=c.ramp_id AND a.type='compound')
        LIMIT 1''').fetchone():
        raise ValueError(f'{label}: chem_props has a missing or non-metabolite analyte')
    source = "coalesce(nullif(lower(trim(chem_data_source)),''),'(missing source)')"
    counts = {name: {'rows': rows, 'source_ids': ids, 'analytes': analytes}
              for name, rows, ids, analytes in db.execute(f'''
                  SELECT {source},count(*),count(DISTINCT chem_source_id),
                         count(DISTINCT ramp_id) FROM chem_props GROUP BY 1''')}
    counts['All sources'] = {'rows': db.execute('SELECT count(*) FROM chem_props').fetchone()[0],
                             'source_ids': db.execute('SELECT count(DISTINCT chem_source_id) FROM chem_props').fetchone()[0],
                             'analytes': db.execute('SELECT count(DISTINCT ramp_id) FROM chem_props').fetchone()[0]}
    coverage = field_coverage(db, 'chem_props', CHEM_FIELDS, columns, (source,))
    return {'present': True, 'counts': counts,
            'field_coverage': {name: values for (name,), values in coverage.items()}}


def source_versions_profile(db, label):
    columns = {row[1] for row in db.execute('PRAGMA table_info(version_info)')}
    if not columns:
        return {'present': False, 'current': {}, 'rows': [], 'columns': []}
    required = {'data_source_id', 'data_source_version', 'status'}
    if not required <= columns:
        raise ValueError(f'{label}: version_info is missing {sorted(required - columns)}')
    names = [row[1] for row in db.execute('PRAGMA table_info(version_info)')]
    current = {}
    rows = []
    order = ','.join(field for field in ('db_mod_date', 'ramp_db_version', 'data_source_id')
                     if field in columns)
    for raw in db.execute('SELECT * FROM version_info' + (f' ORDER BY {order}' if order else '')):
        record = dict(zip(names, raw))
        rows.append(record)
        if str(record['status']).lower() == 'current':
            source = str(record['data_source_id']).strip().lower()
            if not source or source in current:
                raise ValueError(f'{label}: version_info has an empty or duplicate current source: {source}')
            current[source] = record
    return {'present': True, 'current': current, 'rows': rows, 'columns': names}


def entity_status_profile(db, label):
    columns = [row[1] for row in db.execute('PRAGMA table_info(entity_status_info)')]
    if not columns:
        return {'present': False, 'rows': [], 'by_key': {}}
    required = {'status_category', 'entity_source_id', 'entity_source_name', 'entity_count'}
    if not required <= set(columns):
        raise ValueError(f'{label}: entity_status_info is missing {sorted(required - set(columns))}')
    rows = [dict(zip(columns, row)) for row in db.execute(
        'SELECT * FROM entity_status_info ORDER BY status_category,entity_source_id')]
    by_key = {}
    for row in rows:
        key = row['status_category'] + ' / ' + row['entity_source_id']
        if key in by_key:
            raise ValueError(f'{label}: duplicate entity status row {key}')
        by_key[key] = row
    return {'present': True, 'rows': rows, 'by_key': by_key}


def db_version_profile(db, label):
    columns = [row[1] for row in db.execute('PRAGMA table_info(db_version)')]
    if not columns:
        return {'present': False, 'columns': [], 'rows': [], 'latest': None, 'intersections': {}}
    required = {'ramp_version', 'load_timestamp'} | set(DB_VERSION_FIELDS)
    if not required <= set(columns):
        raise ValueError(f'{label}: db_version is missing {sorted(required - set(columns))}')
    rows = [dict(zip(columns, row)) for row in db.execute(
        'SELECT * FROM db_version ORDER BY load_timestamp DESC')]
    latest = rows[0] if rows else None
    intersections = {}
    for field in DB_VERSION_FIELDS:
        value = latest[field] if latest else None
        if value is None:
            intersections[field] = {'state': 'row_absent' if latest is None else 'null'}
            continue
        try:
            parsed = json.loads(value)
            if not isinstance(parsed, list):
                raise ValueError('expected a JSON array')
            by_sets = {}
            for item in parsed:
                if not isinstance(item, dict) or set(item) != {'id', 'sets', 'size'}:
                    raise ValueError('expected {id,sets,size} records')
                if (not isinstance(item['id'], str) or not isinstance(item['sets'], list)
                        or not all(isinstance(source, str) and source for source in item['sets'])
                        or len(set(item['sets'])) != len(item['sets'])
                        or not isinstance(item['size'], int) or isinstance(item['size'], bool)
                        or item['size'] < 0):
                    raise ValueError('invalid intersection ID, sources, or count')
                key = tuple(sorted(item['sets']))
                by_sets[key] = by_sets.get(key, 0) + item['size']
            intersections[field] = {
                'state': 'empty' if not parsed else 'valid',
                'total': sum(by_sets.values()), 'combinations': len(by_sets),
                'records': len(parsed), 'duplicate_combination_records': len(parsed) - len(by_sets),
                'zero_size_records': sum(item['size'] == 0 for item in parsed),
                'sources': sorted({source for sets in by_sets for source in sets}),
                'by_sets': {' & '.join(key): size for key, size in by_sets.items()},
                'combination_data': [{'sets': list(key), 'size': size}
                                     for key, size in by_sets.items()],
            }
        except (json.JSONDecodeError, ValueError, TypeError) as exc:
            intersections[field] = {'state': 'malformed', 'error': str(exc)}
    return {'present': True, 'columns': columns, 'rows': rows,
            'latest': latest, 'intersections': intersections}


def _pathway_scope_ids(db, cutoff, analyte_filter):
    return [row[0] for row in db.execute(f'''
        SELECT a.pathwayRampId FROM analytehaspathway a
        JOIN pathway p ON p.pathwayRampId=a.pathwayRampId
        WHERE p.type!='hmdb' AND {analyte_filter}
        GROUP BY a.pathwayRampId HAVING count(DISTINCT a.rampId)>={cutoff}
        ORDER BY a.pathwayRampId''')]


def _blob_preview(db, blob, anchor_id, cutoff, analyte_filter):
    if blob is None:
        return {'state': 'null', 'bytes': None, 'preview': []}
    try:
        ids = _pathway_scope_ids(db, cutoff, analyte_filter)
        decoded = zlib.decompress(blob).decode('ascii')
        tokens = decoded.split('|')
        index_str, sentinel = tokens[0].split(',')
        index = int(index_str)
        if sentinel != '-' or not 0 <= index < len(ids) or ids[index] != anchor_id:
            raise ValueError('anchor index does not match the eligible pathway order')
        partner_index = index
        preview = []
        for token in tokens[1:4]:
            delta_str, score_str = token.split(',')
            delta, score = int(delta_str), int(score_str)
            partner_index += delta
            if delta <= 0 or not index < partner_index < len(ids) or not 0 < score <= 1000:
                raise ValueError('partner delta or score is invalid')
            partner_id = ids[partner_index]
            row = db.execute('SELECT sourceId,pathwayName FROM pathway WHERE pathwayRampId=?',
                             (partner_id,)).fetchone()
            if row is None:
                raise ValueError(f'partner pathway is missing: {partner_id}')
            preview.append({'pathwayRampId': partner_id, 'sourceId': row[0],
                            'pathwayName': row[1], 'score': score / 1000})
        return {'state': 'valid', 'bytes': len(blob), 'anchor_index': index,
                'preview': preview, 'more_partners': max(0, len(tokens) - 1 - len(preview))}
    except (ValueError, TypeError, UnicodeError, zlib.error) as exc:
        return {'state': 'malformed', 'bytes': len(blob), 'error': str(exc), 'preview': []}


def pathway_overlap_profile(db, label):
    similarity_columns = tuple(row[1] for row in db.execute('PRAGMA table_info(pathway_similarity)'))
    duplicate_columns = tuple(row[1] for row in db.execute('PRAGMA table_info(pathway_duplicates)'))
    result = {'similarity_present': bool(similarity_columns),
              'duplicates_present': bool(duplicate_columns),
              'similarity_rows': None, 'duplicate_rows': None,
              'field_coverage': {}, 'blob_bytes': {}, 'by_type': {},
              'similarity_example': {'state': 'table_absent'},
              'duplicate_example': {'state': 'table_absent'}}
    if similarity_columns:
        required = {'pathwayRampId', 'analyte_blob', 'metabolite_blob', 'gene_blob',
                    'metabolite_count', 'gene_count'}
        if not required <= set(similarity_columns):
            raise ValueError(f'{label}: pathway_similarity is missing {sorted(required - set(similarity_columns))}')
        fields = ('pathwayRampId', 'analyte_blob', 'metabolite_blob', 'gene_blob',
                  'metabolite_count', 'gene_count')
        row = db.execute('SELECT count(*),' + ','.join(f'count({field})' for field in fields) +
                         ',' + ','.join(f'coalesce(sum(length({field})),0)' for field in
                                        ('analyte_blob', 'metabolite_blob', 'gene_blob')) +
                         ' FROM pathway_similarity').fetchone()
        result['similarity_rows'] = row[0]
        result['field_coverage'] = dict(zip(fields, row[1:1 + len(fields)]))
        result['blob_bytes'] = dict(zip(('analyte_blob', 'metabolite_blob', 'gene_blob'),
                                        row[1 + len(fields):]))
        result['by_type'] = {kind: count for kind, count in db.execute('''
            SELECT p.type,count(*) FROM pathway_similarity s
            JOIN pathway p ON p.pathwayRampId=s.pathwayRampId GROUP BY p.type''')}
        path = db.execute('''SELECT pathwayRampId,sourceId,pathwayName,type
            FROM pathway WHERE sourceId=?''', (PATHWAY_SIMILARITY_EXAMPLE,)).fetchall()
        if result['similarity_rows'] == 0:
            result['similarity_example'] = {'state': 'table_empty'}
        elif len(path) > 1:
            raise ValueError(f'{label}: example pathway source ID is not unique')
        elif not path:
            result['similarity_example'] = {'state': 'source_absent'}
        else:
            pathway_id, source_id, name, kind = path[0]
            stored = db.execute('''SELECT analyte_blob,metabolite_blob,gene_blob,
                metabolite_count,gene_count FROM pathway_similarity
                WHERE pathwayRampId=?''', (pathway_id,)).fetchone()
            if stored is None:
                result['similarity_example'] = {'state': 'not_eligible',
                                                'pathwayRampId': pathway_id,
                                                'sourceId': source_id}
            else:
                blobs = stored[:3]
                result['similarity_example'] = {
                    'state': 'available', 'pathwayRampId': pathway_id,
                    'sourceId': source_id, 'pathwayName': name, 'type': kind,
                    'metabolite_count': stored[3], 'gene_count': stored[4],
                    'blobs': {column: _blob_preview(db, blob, pathway_id, cutoff, condition)
                              for (column, cutoff, condition), blob in
                              zip(PATHWAY_OVERLAP_SCOPES, blobs, strict=True)},
                }
    if duplicate_columns:
        required = {'pathwayRampId1', 'pathwayRampId2'}
        if not required <= set(duplicate_columns):
            raise ValueError(f'{label}: pathway_duplicates is missing {sorted(required - set(duplicate_columns))}')
        result['duplicate_rows'] = db.execute('SELECT count(*) FROM pathway_duplicates').fetchone()[0]
        paths = {}
        for source_id in PATHWAY_DUPLICATE_EXAMPLE:
            found = db.execute('''SELECT pathwayRampId,sourceId,pathwayName
                FROM pathway WHERE sourceId=?''', (source_id,)).fetchall()
            if len(found) > 1:
                raise ValueError(f'{label}: duplicate example source ID {source_id} is not unique')
            if found:
                paths[source_id] = found[0]
        if result['duplicate_rows'] == 0:
            result['duplicate_example'] = {'state': 'table_empty'}
        elif len(paths) != 2:
            result['duplicate_example'] = {'state': 'source_absent'}
        else:
            left, right = (paths[source] for source in PATHWAY_DUPLICATE_EXAMPLE)
            found = db.execute('''SELECT pathwayRampId1,pathwayRampId2
                FROM pathway_duplicates WHERE (pathwayRampId1=? AND pathwayRampId2=?)
                   OR (pathwayRampId1=? AND pathwayRampId2=?)''',
                (left[0], right[0], right[0], left[0])).fetchall()
            result['duplicate_example'] = {'state': 'available' if found else 'not_stored',
                                           'left': left, 'right': right}
    return result


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
        catalyzed_columns = {row[1] for row in db.execute('PRAGMA table_info(catalyzed)')}
        catalyzed_present = bool(catalyzed_columns)
        catalyzed_counts = {}
        catalyzed_field_coverage = None
        if catalyzed_present:
            missing = set(CATALYZED_FIELDS) - catalyzed_columns
            if missing:
                raise ValueError(f'{label}: catalyzed is missing {sorted(missing)}')
            if db.execute('''SELECT 1 FROM catalyzed c
                LEFT JOIN analyte m ON m.rampId=c.rampCompoundId AND m.type='compound'
                LEFT JOIN analyte g ON g.rampId=c.rampGeneId AND g.type='gene'
                WHERE m.rampId IS NULL OR g.rampId IS NULL LIMIT 1''').fetchone():
                raise ValueError(f'{label}: catalyzed has a missing or wrong-type analyte')
            rows, metabolites, genes = db.execute('''SELECT count(*),
                count(DISTINCT rampCompoundId),count(DISTINCT rampGeneId) FROM catalyzed''').fetchone()
            catalyzed_counts = {'All associations': {
                'rows': rows, 'metabolites': metabolites, 'genes': genes}}
            catalyzed_field_coverage = field_coverage(
                db, 'catalyzed', CATALYZED_FIELDS, catalyzed_columns, ('1',)
            ).get((1,))
        class_columns = {row[1] for row in db.execute('PRAGMA table_info(metabolite_class)')}
        class_present = bool(class_columns)
        class_counts = {}
        class_field_coverage = {}
        if class_present:
            missing = set(CLASS_FIELDS) - class_columns
            if missing:
                raise ValueError(f'{label}: metabolite_class is missing {sorted(missing)}')
            if db.execute('''SELECT 1 FROM metabolite_class c WHERE NOT EXISTS
                (SELECT 1 FROM analyte a WHERE a.rampId=c.ramp_id AND a.type='compound')
                LIMIT 1''').fetchone():
                raise ValueError(f'{label}: metabolite_class has a missing or non-metabolite analyte')
            for source, level, rows, analytes, source_ids in db.execute('''
                SELECT coalesce(nullif(trim(source),''),'(missing source)'),
                       coalesce(nullif(trim(class_level_name),''),'(missing level)'),
                       count(*),count(DISTINCT ramp_id),count(DISTINCT class_source_id)
                FROM metabolite_class GROUP BY 1,2'''):
                class_counts[f'{source} / {level}'] = {
                    'rows': rows, 'analytes': analytes, 'source_ids': source_ids}
            for source, rows, analytes, source_ids in db.execute('''
                SELECT coalesce(nullif(trim(source),''),'(missing source)'),
                       count(*),count(DISTINCT ramp_id),count(DISTINCT class_source_id)
                FROM metabolite_class GROUP BY 1'''):
                class_counts[f'{source} / All levels'] = {
                    'rows': rows, 'analytes': analytes, 'source_ids': source_ids}
            class_field_coverage = field_coverage(
                db, 'metabolite_class', CLASS_FIELDS, class_columns,
                ('source', 'class_level_name'))
        reactions = reaction_profiles(db, label)
        chemistry = chemistry_profile(db, label)
        versions = source_versions_profile(db, label)
        entity_status = entity_status_profile(db, label)
        db_version = db_version_profile(db, label)
        pathway_overlap = pathway_overlap_profile(db, label)
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
                                               'analytehaspathway', 'analytehasontology'},
                                              {'analyte', 'source', 'analytesynonym', 'pathway', 'ontology',
                                               'analytehaspathway', 'analytehasontology',
                                               'catalyzed', 'metabolite_class'})
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
                'catalyzed_table_present': catalyzed_present,
                'catalyzed_counts': catalyzed_counts,
                'catalyzed_field_coverage': catalyzed_field_coverage,
                'class_table_present': class_present,
                'class_counts': class_counts,
                'class_field_coverage': {source: {level: values for (row_source, level), values
                                                 in class_field_coverage.items() if row_source == source}
                                         for source, _ in class_field_coverage},
                'reaction_tables': reactions,
                'chemistry': chemistry,
                'version_info': versions,
                'entity_status_info': entity_status,
                'db_version': db_version,
                'pathway_overlap': pathway_overlap,
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
    report['catalyzed_example'] = matched_catalyzed_example(report['databases'])
    report['class_examples'] = matched_class_examples(report['databases'])
    report['reaction_examples'] = matched_reaction_examples(report['databases'])
    report['chemistry_examples'] = matched_chemistry_examples(report['databases'])
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


def matched_chemistry_examples(profiles):
    """Choose one reported chemistry identity per provider, favoring cross-build matches."""
    providers = sorted({source for profile in profiles
                        for source in profile['chemistry']['counts'] if source != 'All sources'})
    if not providers:
        return {}
    dbs = _association_databases(profiles)
    try:
        candidates = set()
        for db, profile in zip(dbs, profiles):
            if not profile['chemistry']['present']:
                continue
            for provider_name in providers:
                if provider_name not in profile['chemistry']['counts']:
                    continue
                candidates.update((provider_name, source_id) for (source_id,) in db.execute('''
                    SELECT DISTINCT chem_source_id FROM chem_props
                    WHERE lower(trim(chem_data_source))=? AND chem_source_id IS NOT NULL
                    LIMIT 30''', (provider_name,)))
        found = []
        for db, profile in zip(dbs, profiles):
            matches = {}
            if profile['chemistry']['present']:
                for rowid, source, source_id in db.execute('''
                    SELECT rowid,chem_data_source,chem_source_id FROM chem_props'''):
                    key = ((source or '').strip().lower(), source_id)
                    if key in candidates:
                        matches.setdefault(key, []).append(rowid)
            found.append(matches)
        examples = {}
        for provider_name in providers:
            options = [key for key in candidates if key[0] == provider_name]
            if not options:
                raise ValueError(f'chem_props provider {provider_name} has rows but no chemistry source ID')
            selected = min(options, key=lambda key: (
                -sum(key in rows for rows in found),
                -(key in found[-1]),
                -(key in found[-2]) if len(found) > 1 else 0,
                sum(len(rows.get(key, ())) > 1 for rows in found),
                key[1]))
            rows = []
            for db, profile, matches in zip(dbs, profiles, found):
                if not profile['chemistry']['present']:
                    rows.append({'state': 'table_absent'})
                elif selected not in matches:
                    rows.append({'state': 'source_id_absent'})
                else:
                    columns = [row[1] for row in db.execute('PRAGMA table_info(chem_props)')]
                    values = db.execute('SELECT * FROM chem_props WHERE rowid=?',
                                        (matches[selected][0],)).fetchone()
                    rows.append({'state': 'available', 'matched_rows': len(matches[selected]),
                                 'cells': dict(zip(columns, values))})
            examples[provider_name] = {'chem_source_id': selected[1], 'rows': rows}
        return examples
    finally:
        for db in dbs:
            db.close()


def matched_reaction_examples(profiles):
    """Match rows by source identities, never by release-local RaMP IDs."""
    dbs = _association_databases(profiles)
    try:
        examples = {}
        for table, fields in REACTION_TABLES.items():
            available = [i for i, p in enumerate(profiles)
                         if p['reaction_tables'][table]['present']
                         and p['reaction_tables'][table]['counts']['All rows']['rows']]
            if not available:
                examples[table] = None
                continue
            anchor_fields = REACTION_ANCHORS[table]
            latest = available[-1]
            select = ','.join(anchor_fields)
            candidates = list(dbs[latest].execute(
                f'SELECT DISTINCT {select} FROM {table} WHERE rxn_source_id=? LIMIT 30',
                ('rhea:10000',)))
            if not candidates:
                candidates = list(dbs[latest].execute(
                    f'SELECT DISTINCT {select} FROM {table} LIMIT 100'))
            if not candidates:
                raise ValueError(f'{table} has rows but no source-matched example')

            def row(index, anchor):
                if not profiles[index]['reaction_tables'][table]['present']:
                    return {'state': 'table_absent'}
                clauses = ' AND '.join(f'{field} IS ?' for field in anchor_fields)
                found = list(dbs[index].execute(
                    f'SELECT * FROM {table} WHERE {clauses} LIMIT 2', anchor))
                if len(found) > 1:
                    return {'state': 'ambiguous_row'}
                if not found:
                    reaction_exists = dbs[index].execute(
                        'SELECT 1 FROM reaction WHERE rxn_source_id=? LIMIT 1',
                        (anchor[0],)).fetchone() if profiles[index]['reaction_tables']['reaction']['present'] else None
                    return {'state': 'association_absent' if reaction_exists else 'reaction_absent'}
                columns = [c[1] for c in dbs[index].execute(f'PRAGMA table_info({table})')]
                stored = dict(zip(columns, found[0]))
                return {'state': 'available',
                        'cells': {field: stored[field] if field in stored else {'state': 'column_absent'}
                                  for field in fields}}

            selected = max(candidates, key=lambda anchor: (
                sum(row(i, anchor)['state'] == 'available' for i in available),
                anchor[0] == 'rhea:10000'))
            examples[table] = {'anchor': dict(zip(anchor_fields, selected)),
                               'rows': [row(i, selected) for i in range(len(profiles))]}
        return examples
    finally:
        for db in dbs:
            db.close()


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


def matched_catalyzed_example(profiles):
    available = [i for i, p in enumerate(profiles) if p['catalyzed_table_present']
                 and p['catalyzed_counts']['All associations']['rows']]
    if not available:
        return None
    dbs = _association_databases(profiles)
    try:
        latest = available[-1]
        db = dbs[latest]
        preferred = EXAMPLES[0][2]
        candidates = []

        def add_candidates(compound_id, gene_ids):
            for gene_id in gene_ids:
                protein_ids = [row[0] for row in db.execute('''
                    SELECT DISTINCT sourceId FROM source
                    WHERE rampId=? AND geneOrCompound='gene' AND sourceId LIKE 'uniprot:%'
                    ORDER BY sourceId LIMIT 3''', (gene_id,))]
                for protein_id in protein_ids:
                    candidates.append((compound_id, protein_id))

        preferred_rids = _identity_rows(db, preferred, 'compound')
        if len(preferred_rids) == 1:
            add_candidates(preferred, [row[0] for row in db.execute('''
                SELECT rampGeneId FROM catalyzed WHERE rampCompoundId=? LIMIT 30''',
                (preferred_rids[0],))])
        if not candidates:
            for compound_rid, gene_rid in db.execute('''
                    SELECT rampCompoundId,rampGeneId FROM catalyzed LIMIT 30'''):
                compound_ids = [row[0] for row in db.execute('''
                    SELECT DISTINCT sourceId FROM source
                    WHERE rampId=? AND geneOrCompound='compound' AND sourceId LIKE 'hmdb:%'
                    ORDER BY sourceId LIMIT 2''', (compound_rid,))]
                for compound_id in compound_ids:
                    add_candidates(compound_id, [gene_rid])
        if not candidates:
            raise ValueError('Catalyzed associations exist but no stable HMDB/UniProt example pair is available')

        def row(index, compound_id, protein_id):
            if not profiles[index]['catalyzed_table_present']:
                return {'state': 'table_absent'}
            compounds = _identity_rows(dbs[index], compound_id, 'compound')
            genes = _identity_rows(dbs[index], protein_id, 'gene')
            if len(compounds) != 1 or len(genes) != 1:
                return {'state': 'identity_absent' if not compounds or not genes else 'ambiguous_identity'}
            match = dbs[index].execute('''SELECT rampCompoundId,rampGeneId,proteinType
                FROM catalyzed WHERE rampCompoundId=? AND rampGeneId=? LIMIT 1''',
                (compounds[0], genes[0])).fetchone()
            return {'state': 'available' if match else 'association_absent',
                    'cells': dict(zip(CATALYZED_FIELDS, match)) if match else
                    {'rampCompoundId': compounds[0], 'rampGeneId': genes[0], 'proteinType': None}}

        selected = max(dict.fromkeys(candidates), key=lambda pair: (
            sum(row(i, *pair)['state'] == 'available' for i in available),
            pair[0] == preferred))
        return {'compound_source_id': selected[0], 'protein_source_id': selected[1],
                'rows': [row(i, *selected) for i in range(len(profiles))]}
    finally:
        for db in dbs:
            db.close()


def matched_class_examples(profiles):
    providers = sorted({key.split(' / ', 1)[0] for p in profiles
                        for key in p['class_counts']})
    dbs = _association_databases(profiles)
    try:
        examples = {}
        for provider in providers:
            available = [i for i, p in enumerate(profiles)
                         if p['class_table_present']
                         and p['class_counts'].get(f'{provider} / All levels', {}).get('rows')]
            latest = available[-1]
            db = dbs[latest]
            preferred = EXAMPLES[0][2]
            candidates = list(db.execute('''
                SELECT DISTINCT class_source_id,class_level_name FROM metabolite_class
                WHERE source=? AND class_source_id=? ORDER BY class_level_name LIMIT 20''',
                (provider, preferred)))
            if not candidates:
                candidates = list(db.execute('''
                    SELECT DISTINCT class_source_id,class_level_name FROM metabolite_class
                    WHERE source=? ORDER BY class_source_id,class_level_name LIMIT 20''',
                    (provider,)))
            if not candidates:
                raise ValueError(f'Class rows exist but no example is available: {provider}')

            def row(index, source_id, level):
                if not profiles[index]['class_table_present']:
                    return {'state': 'table_absent'}
                rids = _identity_rows(dbs[index], source_id, 'compound')
                if len(rids) != 1:
                    return {'state': 'identity_absent' if not rids else 'ambiguous_identity'}
                matches = list(dbs[index].execute('''
                    SELECT ramp_id,class_source_id,class_level_name,class_name,source
                    FROM metabolite_class WHERE ramp_id=? AND class_source_id=?
                      AND class_level_name=? AND source=? LIMIT 2''',
                    (rids[0], source_id, level, provider)))
                if len(matches) > 1:
                    return {'state': 'ambiguous_class_row'}
                return {'state': 'available' if matches else 'association_absent',
                        'cells': dict(zip(CLASS_FIELDS, matches[0])) if matches else
                        {'ramp_id': rids[0], 'class_source_id': source_id,
                         'class_level_name': level, 'class_name': None, 'source': provider}}

            selected = max(candidates, key=lambda pair: (
                sum(row(i, *pair)['state'] == 'available' for i in available),
                pair[0] == preferred))
            examples[provider] = {'source_id': selected[0], 'class_level_name': selected[1],
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
                            ('catalyzed_counts', ('rows', 'metabolites', 'genes')),
                            ('class_counts', ('rows', 'analytes', 'source_ids')),
                            ('chemistry_counts', ('rows', 'source_ids', 'analytes')),
                            ('reaction', ('rows', 'reactions')),
                            ('reaction2met', ('rows', 'reactions', 'analytes')),
                            ('reaction2protein', ('rows', 'reactions', 'analytes')),
                            ('reaction_ec_class', ('rows', 'reactions')),
                            ('lookup_quality', (None,))]:
        reaction_table = section if section in REACTION_TABLES else None
        chemistry_table = section == 'chemistry_counts'
        synonym_kind = {'compound_synonyms': 'compound', 'gene_synonyms': 'gene'}.get(section)
        association_table = {'ontology_associations': 'ontology_table_present',
                             'pathway_associations': 'pathway_association_table_present',
                             'catalyzed_counts': 'catalyzed_table_present',
                             'class_counts': 'class_table_present'}.get(section)
        keys = sorted({key for p in profiles for key in (p['reaction_tables'][reaction_table]['counts']
                       if reaction_table else p['synonym_sources'][synonym_kind]
                       if synonym_kind else p['chemistry']['counts']
                       if chemistry_table else p[section])})
        if section == 'ramp_ids':
            keys = ['Metabolite RaMP IDs', 'Gene/protein RaMP IDs']
        if association_table:
            keys.sort(key=lambda key: (not key.endswith(('All types', 'All sources')), key))
        if section == 'class_counts':
            keys.sort(key=lambda key: (key.split(' / ', 1)[0], not key.endswith('All levels'), key))
        if reaction_table:
            keys.sort(key=lambda key: (key != 'All rows', key))
        if chemistry_table:
            keys.sort(key=lambda key: (key != 'All sources', key))
        result[section] = []
        for key in keys:
            for field in fields:
                values = [(None if not p['reaction_tables'][reaction_table]['present']
                           or reaction_table == 'reaction2protein' and key.startswith('Reviewed / ')
                           and 'is_reviewed' not in p['reaction_tables'][reaction_table]['columns'] else
                           p['reaction_tables'][reaction_table]['counts'].get(key, {}).get(field, 0))
                          if reaction_table else
                          (None if not p['synonym_table_present'] else
                           p['synonym_sources'][synonym_kind].get(key, {}).get(field, 0))
                          if synonym_kind else
                          (None if not p['chemistry']['present'] else
                           p['chemistry']['counts'].get(key, {}).get(field, 0))
                          if chemistry_table else
                          (None if not p[association_table] else
                           p[section].get(key, {}).get(field, 0)) if association_table else
                          p[section].get(key, 0)
                          for p in profiles]
                metric_field = ('RaMP IDs' if field in ('ramp_ids', 'analytes') else
                                'metabolite RaMP IDs' if field == 'metabolites' else
                                'gene/protein RaMP IDs' if field == 'genes' else
                                'source IDs' if field == 'source_ids' else
                                'metabolite RaMP IDs' if section == 'reaction2met' and field == 'analytes' else
                                'gene/protein RaMP IDs' if section == 'reaction2protein' and field == 'analytes' else
                                'Rhea reactions' if field == 'reactions' else
                                'ontology terms' if section == 'ontology_associations' and field == 'targets' else
                                'pathways' if section == 'pathway_associations' and field == 'targets' else field)
                result[section].append({'metric': f'{key} / {metric_field}' if field else key,
                                        'cells': [cell(value, values[i-1] if i else None)
                                                  for i, value in enumerate(values)]})
    return result


def render_upset_svg(item, *, title, shared_max, limit=20):
    """Render an accessible, dependency-free exact-intersection UpSet plot."""
    combinations = sorted(item['combination_data'],
                          key=lambda row: (-row['size'], tuple(row['sets'])))[:limit]
    source_totals = {source: sum(row['size'] for row in item['combination_data']
                                 if source in row['sets']) for source in item['sources']}
    sources = sorted(item['sources'], key=lambda source: (-source_totals[source], source))
    left, step, bar_bottom, bar_height, row_step = 215, 47, 188, 135, 25
    matrix_top = 218
    width = max(620, left + len(combinations) * step + 30)
    height = matrix_top + len(sources) * row_step + 24
    esc = lambda value: html.escape(str(value), quote=True)
    plotted = sum(row['size'] for row in combinations)
    other = item['total'] - plotted
    pieces = [f'<svg class="upset-svg" viewBox="0 0 {width} {height}" '
              f'width="{width}" height="{height}" role="img" '
              f'aria-label="{esc(title)}: top {len(combinations)} of {item["combinations"]} '
              f'source combinations, representing {plotted:,} of {item["total"]:,} analytes">',
              f'<title>{esc(title)} source intersections</title>',
              f'<desc>Bars show exact combination sizes on a shared linear scale up to '
              f'{shared_max:,}. Filled dots show sources in each combination. '
              f'{other:,} analytes are in combinations outside this plot.</desc>',
              f'<text x="{left}" y="24" class="upset-axis-label">'
              f'Exact intersection size · shared scale 0–{shared_max:,}</text>',
              f'<line x1="{left}" y1="{bar_bottom}" x2="{width-20}" '
              f'y2="{bar_bottom}" class="upset-axis"/>']
    for index, row in enumerate(combinations):
        x = left + index * step + 7
        bar = round(bar_height * row['size'] / shared_max) if shared_max else 0
        if row['size'] and bar < 2:
            bar = 2
        combination_label = ', '.join(row['sets']) or '(no sources)'
        pieces.append(f'<rect x="{x}" y="{bar_bottom-bar}" width="27" '
                      f'height="{bar}" class="upset-bar"><title>'
                      f'{esc(combination_label)}: {row["size"]:,} analytes</title></rect>')
        pieces.append(f'<text x="{x+13}" y="{bar_bottom-bar-6}" '
                      f'text-anchor="middle" class="upset-count">{row["size"]:,}</text>')
        selected_rows = [i for i, source in enumerate(sources) if source in row['sets']]
        if len(selected_rows) > 1:
            y1 = matrix_top + selected_rows[0] * row_step
            y2 = matrix_top + selected_rows[-1] * row_step
            pieces.append(f'<line x1="{x+13}" y1="{y1}" x2="{x+13}" '
                          f'y2="{y2}" class="upset-connector"/>')
    for source_index, source in enumerate(sources):
        y = matrix_top + source_index * row_step
        pieces.append(f'<text x="{left-15}" y="{y+4}" text-anchor="end" '
                      f'class="upset-source">{esc(source)} · {source_totals[source]:,}</text>')
        pieces.append(f'<line x1="{left}" y1="{y+row_step/2}" x2="{width-20}" '
                      f'y2="{y+row_step/2}" class="upset-row"/>')
        for index, row in enumerate(combinations):
            x = left + index * step + 20
            selected = source in row['sets']
            pieces.append(f'<circle cx="{x}" cy="{y}" r="5" '
                          f'class="upset-dot{" selected" if selected else ""}"/>')
    pieces.append('</svg>')
    return ''.join(pieces)


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
             ('pathway_examples', 'Example pathway associations'),
             ('catalyzed_counts', 'Catalyzed metabolite–gene associations'),
             ('catalyzed_example', 'Example catalyzed association'),
             ('class_counts', 'Metabolite classes by source and level'),
             ('class_examples', 'Example metabolite class rows'),
             ('chemistry_counts', 'Chemical properties by source'),
             ('chemistry_examples', 'Example chemical property rows'),
             ('reaction', 'Reactions'), ('reaction_example', 'Example reaction'),
             ('reaction2met', 'Reaction–metabolite participants'),
             ('reaction2met_example', 'Example reaction–metabolite row'),
             ('reaction2protein', 'Reaction–protein participants'),
             ('reaction2protein_example', 'Example reaction–protein row'),
             ('reaction_ec_class', 'Reaction EC classes'),
             ('reaction_ec_class_example', 'Example reaction EC class row'),
             ('version_info', 'Recorded source versions')]
    panel_for = {'ramp_ids': 'analyte',
                 'metabolite_sources': 'source', 'gene_sources': 'source',
                 'source_examples': 'source', 'lookup_quality': 'source',
                 'compound_synonyms': 'analytesynonym', 'gene_synonyms': 'analytesynonym',
                 'synonym_examples': 'analytesynonym',
                 'ontology_associations': 'analytehasontology',
                 'ontology_examples': 'analytehasontology',
                 'pathway_associations': 'analytehaspathway',
                 'pathway_examples': 'analytehaspathway',
                 'catalyzed_counts': 'catalyzed', 'catalyzed_example': 'catalyzed',
                 'class_counts': 'metabolite_class', 'class_examples': 'metabolite_class',
                 'chemistry_counts': 'chem_props', 'chemistry_examples': 'chem_props',
                 'version_info': 'version_info',
                 **{table: table for table in REACTION_TABLES},
                 **{f'{table}_example': table for table in REACTION_TABLES}}
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
             '<title>RaMP lookup and association comparison</title>',
             '<style>body{font:16px/1.5 system-ui,sans-serif;margin:0;background:#f5f7fa;color:#182635}main{max-width:1500px;margin:auto;padding:24px}p{max-width:100ch}.scroll{overflow-x:auto;background:white;border:1px solid #d2dbe3}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px 12px;border-bottom:1px solid #dde4eb;text-align:right;min-width:140px}th:first-child,td:first-child{text-align:left;min-width:300px}th{background:#eaf0f6}tbody tr:nth-child(even){background:#f8fafc}tbody tr.coverage{background:#eef3f7;color:#42566a;font-size:12px}tbody tr.coverage td{padding-top:4px;padding-bottom:8px}small{display:block;color:#35475a}article{background:white;border-left:4px solid #52738c;padding:10px 16px;margin:10px 0}.tabs{display:flex;flex-wrap:wrap;gap:6px;border-bottom:2px solid #c5d2de;margin:24px 0 18px}.tabs button{font:inherit;color:#23527c;background:#eaf0f6;border:1px solid #c5d2de;border-bottom:0;border-radius:6px 6px 0 0;padding:9px 18px;cursor:pointer}.tabs button[aria-selected=true]{background:white;color:#182635;font-weight:700;border-color:#52738c}.tabs button:focus-visible{outline:3px solid #23527c;outline-offset:2px}.tab-panel[hidden]{display:none}.example th,.example td{text-align:left}.example td{overflow-wrap:anywhere;white-space:normal}.example th:first-child,.example td:first-child{min-width:220px}.synonym-example td{max-width:650px}@media(max-width:600px){main{padding:12px}.tabs button{flex:1 1 auto}}@media print{.tabs{display:none}.tab-panel[hidden]{display:block}}</style><main>',
             '<style>.scope-choices{display:flex;flex-wrap:wrap;gap:8px;margin:14px 0}'
             '.scope-choices button{font:inherit;padding:7px 12px;background:#eaf0f6;color:#23527c;'
             'border:1px solid #b8c9d8;border-radius:5px;cursor:pointer}'
             '.scope-choices button[aria-pressed=true]{background:#23527c;color:white}'
             '.upset-card{background:white;border:1px solid #d2dbe3;border-radius:5px;padding:12px;margin:14px 0}'
             '.upset-scroll{overflow-x:auto}.upset-svg{max-width:none;height:auto}'
             '.upset-bar,.upset-dot.selected{fill:#276895}.upset-dot{fill:white;stroke:#8da4b6;stroke-width:1.5}'
             '.upset-connector{stroke:#276895;stroke-width:2}.upset-axis{stroke:#597080;stroke-width:1.5}'
             '.upset-row{stroke:#e5ebf0;stroke-width:1}.upset-count{font:11px system-ui,sans-serif;fill:#193b57}'
             '.upset-source{font:12px system-ui,sans-serif;fill:#233c50}'
             '.upset-axis-label{font:12px system-ui,sans-serif;fill:#42566a}'
             '.upset-card p{margin:6px 0 12px}.db-version-raw pre{white-space:pre-wrap;overflow-wrap:anywhere}'
             '</style>',
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
                  ('db_version', 'Version & UpSet'),
                  ('analytesynonym', 'Synonyms'), ('analytehasontology', 'Ontology'),
                  ('analytehaspathway', 'Pathways'),
                  ('pathway_similarity', 'Pathway similarity'),
                  ('catalyzed', 'Catalyzed'),
                  ('metabolite_class', 'Metabolite class'),
                  ('chem_props', 'Chemical properties'),
                  ('reaction', 'Reaction'), ('reaction2met', 'Reaction–metabolite'),
                  ('reaction2protein', 'Reaction–protein'),
                  ('reaction_ec_class', 'Reaction EC class'),
                  ('version_info', 'Source versions'),
                  ('entity_status_info', 'Entity counts'),
                  ('about', 'About this report')]
    parts.append('<div class="tabs" role="tablist" aria-label="SQLite tables">' +
                 ''.join(f'<button type="button" role="tab" id="tab-{panel}" aria-controls="{panel}" aria-selected="{str(index == 0).lower()}" tabindex="{0 if index == 0 else -1}">{esc(label)}</button>'
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
        if key == 'catalyzed_example':
            example = report['catalyzed_example']
            if example is None:
                parts.append('<p>No catalyzed associations are available for an example.</p>')
                continue
            parts.append(f'<p>Matched across builds by metabolite <code>{esc(example["compound_source_id"])}</code> and protein <code>{esc(example["protein_source_id"])}</code>. The table has no source column. “Unknown” protein type is counted as unpopulated in field coverage.</p>')
            parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                         ''.join(f'<th>{esc(field)}</th>' for field in CATALYZED_FIELDS) +
                         '<th>Association</th></tr></thead><tbody>')
            for index, (profile, row) in enumerate(zip(profiles, example['rows'])):
                status = {'available': 'Present', 'association_absent': 'Absent',
                          'identity_absent': 'Endpoint source ID absent',
                          'ambiguous_identity': 'Endpoint source ID maps to multiple RaMP IDs',
                          'table_absent': 'Table absent'}[row['state']]
                parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                for field in CATALYZED_FIELDS:
                    value = row.get('cells', {}).get(field)
                    parts.append(f'<td>{esc(value if value is not None else "—")}</td>')
                parts.append(f'<td>{esc(status)}</td></tr>')
                parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                previous = profiles[index - 1]['catalyzed_field_coverage'] if index else None
                parts.extend(coverage_td(profile['catalyzed_field_coverage'], previous, field)
                             for field in CATALYZED_FIELDS)
                parts.append('<td>—</td></tr>')
            parts.append('</tbody></table></div>')
            continue
        if key == 'class_examples':
            parts.append('<p>Each input’s example follows a reported metabolite ID and class level across builds. All stored columns are shown. Field coverage applies to the same source and level across the whole table.</p>')
            for provider, example in report['class_examples'].items():
                parts.append(f'<details><summary>{esc(provider)} · <code>{esc(example["source_id"])}</code> · {esc(example["class_level_name"])}</summary>')
                parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in CLASS_FIELDS) +
                             '<th>Association</th></tr></thead><tbody>')
                for index, (profile, row) in enumerate(zip(profiles, example['rows'])):
                    status = {'available': 'Present', 'association_absent': 'Absent',
                              'identity_absent': 'Metabolite source ID absent',
                              'ambiguous_identity': 'Source ID maps to multiple RaMP IDs',
                              'ambiguous_class_row': 'Multiple class rows',
                              'table_absent': 'Table absent'}[row['state']]
                    parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                    for field in CLASS_FIELDS:
                        value = row.get('cells', {}).get(field)
                        parts.append(f'<td>{esc(value if value is not None else "—")}</td>')
                    parts.append(f'<td>{esc(status)}</td></tr>')
                    coverage = profile['class_field_coverage'].get(provider, {}).get(example['class_level_name'])
                    previous = (profiles[index - 1]['class_field_coverage'].get(provider, {})
                                .get(example['class_level_name']) if index else None)
                    parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                    parts.extend(coverage_td(coverage, previous, field) for field in CLASS_FIELDS)
                    parts.append('<td>—</td></tr>')
                parts.append('</tbody></table></div></details>')
            continue
        if key == 'chemistry_examples':
            parts.append('<p>Each example follows one reported chemistry ID within a source across builds. '
                         'All stored columns are shown. The coverage row measures populated values across '
                         'all chemistry rows from that source; RaMP IDs can change between builds.</p>')
            for provider_name, example in report['chemistry_examples'].items():
                parts.append(f'<details><summary>{esc(provider_name)} · '
                             f'<code>{esc(example["chem_source_id"])}</code></summary>')
                parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in CHEM_FIELDS) +
                             '<th>Row</th></tr></thead><tbody>')
                for index, (profile, row) in enumerate(zip(profiles, example['rows'])):
                    parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                    for field in CHEM_FIELDS:
                        value = row.get('cells', {}).get(field)
                        display = ('—' if row['state'] != 'available' else
                                   'Column absent' if field not in row['cells'] else
                                   'NULL' if value is None else value)
                        parts.append(f'<td>{esc(display)}</td>')
                    status = {'available': 'Present', 'source_id_absent': 'Source ID absent',
                              'table_absent': 'Table absent'}[row['state']]
                    if row.get('matched_rows', 0) > 1:
                        status += f' ({row["matched_rows"]} matching rows; first shown)'
                    parts.append(f'<td>{esc(status)}</td></tr>')
                    coverage = profile['chemistry']['field_coverage'].get(provider_name)
                    previous = (profiles[index - 1]['chemistry']['field_coverage'].get(provider_name)
                                if index else None)
                    parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                    parts.extend(coverage_td(coverage, previous, field) for field in CHEM_FIELDS)
                    parts.append('<td>—</td></tr>')
                parts.append('</tbody></table></div></details>')
            continue
        if key == 'version_info':
            parts.append('<p>The first table compares current recorded versions by source. '
                         'Every stored row, including archived versions, appears below with all columns. '
                         'Registry snapshot IDs are shown in the full rows when the column exists.</p>')
            sources = sorted({source for profile in profiles
                              for source in profile['version_info']['current']})
            parts.append('<div class="scroll"><table><thead><tr><th>Source</th>' +
                         ''.join(f'<th>{esc(profile["label"])}</th>' for profile in profiles) +
                         '</tr></thead><tbody>')
            for source in sources:
                parts.append(f'<tr><td>{esc(source)}</td>')
                for profile in profiles:
                    version = profile['version_info']['current'].get(source)
                    display = (str(version['data_source_version']) if version else
                               'No current row' if profile['version_info']['present'] else 'Table absent')
                    parts.append(f'<td>{esc(display)}</td>')
                parts.append('</tr>')
            parts.append('</tbody></table></div>')
            columns = list(dict.fromkeys(column for profile in profiles
                                         for column in profile['version_info']['columns']))
            for index, profile in enumerate(profiles):
                info = profile['version_info']
                parts.append(f'<details{" open" if index >= len(profiles) - 2 else ""}>'
                             f'<summary>{esc(profile["label"])} · {len(info["rows"]):,} stored rows</summary>')
                if not info['present']:
                    parts.append('<p>Table absent.</p></details>')
                    continue
                parts.append('<div class="scroll"><table class="example"><thead><tr>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in columns) +
                             '</tr></thead><tbody>')
                for row in info['rows']:
                    parts.append('<tr>' + ''.join(
                        f'<td>{esc("Column absent" if field not in row else "NULL" if row[field] is None else row[field])}</td>'
                        for field in columns) + '</tr>')
                parts.append('</tbody></table></div></details>')
            continue
        if key.endswith('_example') and key[:-8] in REACTION_TABLES:
            table = key[:-8]
            example = report['reaction_examples'][table]
            if example is None:
                parts.append('<p>No reaction rows are available for an example.</p>')
                continue
            anchor = ', '.join(f'{field}={value}' for field, value in example['anchor'].items())
            parts.append(f'<p>Matched across builds by source identity: <code>{esc(anchor)}</code>. '
                         'RaMP IDs are assigned within each build. Field coverage counts populated values '
                         'across the whole table; −1 means unknown and is treated as unpopulated.</p>')
            parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                         ''.join(f'<th>{esc(field)}</th>' for field in REACTION_TABLES[table]) +
                         '<th>Row</th></tr></thead><tbody>')
            for index, (profile, row) in enumerate(zip(profiles, example['rows'])):
                status = {'available': 'Present', 'association_absent': 'Row absent',
                          'reaction_absent': 'Rhea reaction absent',
                          'ambiguous_row': 'Source identity matches multiple rows',
                          'table_absent': 'Table absent'}[row['state']]
                parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                for field in REACTION_TABLES[table]:
                    value = row.get('cells', {}).get(field)
                    display = ('Column absent' if isinstance(value, dict) and
                               value.get('state') == 'column_absent' else
                               '—' if value is None else value)
                    parts.append(f'<td>{esc(display)}</td>')
                parts.append(f'<td>{esc(status)}</td></tr>')
                coverage = profile['reaction_tables'][table]['field_coverage']
                previous = (profiles[index - 1]['reaction_tables'][table]['field_coverage']
                            if index else None)
                parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                parts.extend(coverage_td(coverage, previous, field)
                             for field in REACTION_TABLES[table])
                parts.append('<td>—</td></tr>')
            parts.append('</tbody></table></div>')
            continue
        if key == 'catalyzed_counts':
            parts.append('<p>Rows are distinct metabolite–gene pairs. The table has no source column; the counts show all catalyzed associations and the number of distinct metabolite and gene/protein RaMP IDs they cover.</p>')
        if key == 'class_counts':
            parts.append('<p>Rows are metabolite classifications grouped by their stored source and class level. An analyte may have several levels, so distinct RaMP ID counts across levels do not add together.</p>')
        if key == 'chemistry_counts':
            parts.append('<p>Rows count stored property records. Source IDs count distinct reported chemistry identifiers; RaMP IDs count distinct metabolites with chemistry from each source. A metabolite may have chemistry from several sources.</p>')
        if key in REACTION_TABLES:
            parts.append('<p>Rows and distinct Rhea reactions are counted within this table. '
                         'Participant tables also count distinct RaMP analytes; those IDs can change across builds. '
                         '“Missing reaction target” counts links whose stored reaction RaMP ID has no reaction row. '
                         'Status: 1 approved, 0 preliminary, −1 obsolete. '
                         'Substrate/product: 0 left, 1 right. For cofactor, review, and human-scope flags, '
                         '−1 means unknown. Group counts can overlap on a reaction.</p>')
        if key == 'lookup_quality':
            parts.append('<p>A reported ID can map to several gene/protein RaMP IDs. These are candidate lookups, not guaranteed unique resolutions.</p>')
        parts.append('<div class="scroll"><table><thead><tr><th>Metric</th>' + ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
        for row in sections[key]:
            parts.append(f'<tr><td>{esc(row["metric"])}</td>')
            for profile, datum in zip(profiles, row['cells']):
                value, delta, percent = datum['value'], datum['delta'], datum['percent_change']
                style = change_style(delta, percent if percent is not None and math.isfinite(percent) else 100)
                missing = ('Column absent' if key == 'reaction2protein'
                           and row['metric'].startswith('Reviewed / ')
                           and profile['reaction_tables'][key]['present']
                           and 'is_reviewed' not in profile['reaction_tables'][key]['columns']
                           else 'Table absent' if (
                    key in ('compound_synonyms', 'gene_synonyms') and not profile['synonym_table_present']
                    or key == 'ontology_associations' and not profile['ontology_table_present']
                    or key == 'pathway_associations' and not profile['pathway_association_table_present']
                    or key == 'catalyzed_counts' and not profile['catalyzed_table_present']
                    or key == 'class_counts' and not profile['class_table_present']
                    or key == 'chemistry_counts' and not profile['chemistry']['present']
                    or key in REACTION_TABLES and not profile['reaction_tables'][key]['present'])
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
    parts.append('<section class="tab-panel" id="entity_status_info" role="tabpanel" aria-labelledby="tab-entity_status_info">'
                 '<h2>Entity counts by category and source</h2>'
                 '<p>These are the stored rows in <code>entity_status_info</code>, used by the R package count table. '
                 'The analyte-overlap UpSet plot uses intersection JSON shown in the Database version and overlaps tab. '
                 'A missing row differs from a stored zero; changes compare with the preceding build only when both rows exist.</p>')
    status_keys = sorted({(row['status_category'], row['entity_source_id'])
                          for p in profiles for row in p['entity_status_info']['rows']})
    if status_keys:
        parts.append('<div class="scroll"><table><thead><tr><th>Category / source ID</th>' +
                     ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
        for category, source in status_keys:
            key = category + ' / ' + source
            parts.append(f'<tr><td>{esc(category)} / <code>{esc(source)}</code></td>')
            for index, p in enumerate(profiles):
                status = p['entity_status_info']
                row = status['by_key'].get(key)
                previous = profiles[index - 1]['entity_status_info']['by_key'].get(key) if index else None
                if row is None:
                    parts.append('<td>' + ('Table absent' if not status['present'] else 'Row absent') + '</td>')
                    continue
                count = row['entity_count']
                delta = count - previous['entity_count'] if previous is not None else None
                pct = (100 * delta / previous['entity_count'] if delta is not None and previous['entity_count']
                       else 100 if delta is not None and count else None)
                parts.append(f'<td{change_style(delta, pct or 0)}>{count:,}<small>{esc(row["entity_source_name"])}' +
                             (f' · Δ {delta:+,}' if delta is not None else '') + '</small></td>')
            parts.append('</tr>')
        parts.append('</tbody></table></div>')
    else:
        parts.append('<p>No entity-count rows are stored in the compared databases.</p>')
    parts.append('<h3>All stored rows</h3>')
    for p in profiles:
        status = p['entity_status_info']
        parts.append(f'<details><summary>{esc(p["label"])} · ' +
                     ('table absent' if not status['present'] else f'{len(status["rows"]):,} rows') + '</summary>')
        if status['present']:
            parts.append('<div class="scroll"><table><thead><tr><th>status_category</th><th>entity_source_id</th>'
                         '<th>entity_source_name</th><th>entity_count</th></tr></thead><tbody>')
            for row in status['rows']:
                parts.append('<tr>' + ''.join(f'<td>{esc(row[field])}</td>' for field in
                             ('status_category', 'entity_source_id', 'entity_source_name', 'entity_count')) + '</tr>')
            parts.append('</tbody></table></div>')
        parts.append('</details>')
    parts.append('</section>')
    parts.append('<section class="tab-panel" id="pathway_similarity" role="tabpanel" '
                 'aria-labelledby="tab-pathway_similarity"><h2>Pathway similarity and exact duplicates</h2>'
                 '<p>Only non-SMPDB pathways qualify. Combined-analyte rows require at least 10 distinct '
                 'RaMP IDs; metabolite and gene rows each require at least 5. Each blob stores positive '
                 'Jaccard similarities as compressed pathway-index deltas and integer scores divided by 1,000. '
                 'NULL means this pathway did not qualify for that scope. An empty table has no computed results.</p>')
    overlap_metrics = (
        ('Similarity rows', lambda x: x['similarity_rows'], 'similarity_present'),
        ('Combined blobs', lambda x: x['field_coverage'].get('analyte_blob'), 'similarity_present'),
        ('Metabolite blobs', lambda x: x['field_coverage'].get('metabolite_blob'), 'similarity_present'),
        ('Gene blobs', lambda x: x['field_coverage'].get('gene_blob'), 'similarity_present'),
        ('Duplicate pairs', lambda x: x['duplicate_rows'], 'duplicates_present'),
        ('Compressed BLOB bytes', lambda x: sum(x['blob_bytes'].values()), 'similarity_present'),
    )
    parts.append('<div class="scroll"><table><thead><tr><th>Stored data</th>' +
                 ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
    for title, value_of, present_key in overlap_metrics:
        parts.append(f'<tr><td>{esc(title)}</td>')
        for index, p in enumerate(profiles):
            data = p['pathway_overlap']
            if not data[present_key]:
                parts.append('<td>Table absent</td>')
                continue
            value = value_of(data)
            previous_data = profiles[index - 1]['pathway_overlap'] if index else None
            previous = value_of(previous_data) if previous_data and previous_data[present_key] else None
            delta = value - previous if previous is not None else None
            pct = (100 * delta / previous if delta is not None and previous else
                   100 if delta is not None and value else 0)
            suffix = ' bytes' if title == 'Compressed BLOB bytes' else ''
            parts.append(f'<td{change_style(delta, pct)}>{value:,}{suffix}' +
                         (f'<small>Δ {delta:+,}</small>' if delta is not None else '') +
                         ('<small>Present but empty</small>' if value == 0 and title in
                          ('Similarity rows', 'Duplicate pairs') else '') + '</td>')
        parts.append('</tr>')
    parts.append('</tbody></table></div>')
    pathway_types = sorted({kind for p in profiles
                            for kind in p['pathway_overlap']['by_type']})
    if pathway_types:
        parts.append('<h3>Similarity rows by pathway type</h3><div class="scroll"><table><thead>'
                     '<tr><th>Pathway type</th>' +
                     ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) +
                     '</tr></thead><tbody>')
        for kind in pathway_types:
            parts.append(f'<tr><td>{esc(kind)}</td>')
            for index, p in enumerate(profiles):
                data = p['pathway_overlap']
                if not data['similarity_present']:
                    parts.append('<td>Table absent</td>')
                    continue
                value = data['by_type'].get(kind, 0)
                previous_data = profiles[index - 1]['pathway_overlap'] if index else None
                previous = previous_data['by_type'].get(kind, 0) if previous_data and previous_data['similarity_present'] else None
                delta = value - previous if previous is not None else None
                pct = (100 * delta / previous if delta is not None and previous else
                       100 if delta is not None and value else 0)
                parts.append(f'<td{change_style(delta, pct)}>{value:,}' +
                             (f'<small>Δ {delta:+,}</small>' if delta is not None else '') + '</td>')
            parts.append('</tr>')
        parts.append('</tbody></table></div>')
    parts.append('<h3>Column coverage</h3><div class="scroll"><table><thead><tr><th>pathway_similarity column</th>' +
                 ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
    for field in ('pathwayRampId', 'analyte_blob', 'metabolite_blob', 'gene_blob',
                  'metabolite_count', 'gene_count'):
        parts.append(f'<tr><td><code>{esc(field)}</code></td>')
        for index, p in enumerate(profiles):
            data = p['pathway_overlap']
            if not data['similarity_present']:
                parts.append('<td>Table absent</td>')
                continue
            total = data['similarity_rows']
            filled = data['field_coverage'][field]
            previous_data = profiles[index - 1]['pathway_overlap'] if index else None
            previous_pct = (100 * previous_data['field_coverage'][field] /
                            previous_data['similarity_rows'] if previous_data and
                            previous_data['similarity_present'] and previous_data['similarity_rows'] else None)
            pct = 100 * filled / total if total else None
            points = pct - previous_pct if pct is not None and previous_pct is not None else None
            parts.append(f'<td{change_style(points, points or 0)}>' +
                         (f'{filled:,} / {total:,} ({pct:.1f}%)' if pct is not None else '0 rows') +
                         (f'<small>Δ {points:+.1f} pp</small>' if points is not None else '') + '</td>')
        parts.append('</tr>')
    parts.append('</tbody></table></div>')
    parts.append(f'<h3>Matched pathway example · {esc(PATHWAY_SIMILARITY_EXAMPLE)}</h3>')
    for p in profiles:
        example = p['pathway_overlap']['similarity_example']
        parts.append(f'<details><summary>{esc(p["label"])} · {esc(example["state"].replace("_", " "))}</summary>')
        if example['state'] == 'available':
            parts.append(f'<p><code>{esc(example["pathwayRampId"])}</code> · '
                         f'{esc(example["pathwayName"])} · {esc(example["type"])}; '
                         f'metabolite_count={esc(example["metabolite_count"])}; '
                         f'gene_count={esc(example["gene_count"])}</p>')
            for field in ('analyte_blob', 'metabolite_blob', 'gene_blob'):
                blob = example['blobs'][field]
                parts.append(f'<p><code>{field}</code>: {esc(blob["state"])}' +
                             (f' · {blob["bytes"]:,} compressed bytes' if blob['bytes'] is not None else '') +
                             (f' · {esc(blob["error"])}' if blob['state'] == 'malformed' else '') + '</p>')
                if blob['state'] == 'valid':
                    if blob['preview']:
                        parts.append('<div class="scroll"><table><thead><tr><th>Partner source ID</th>'
                                     '<th>Partner RaMP ID</th><th>Similarity</th></tr></thead><tbody>')
                        for partner in blob['preview']:
                            parts.append(f'<tr><td>{esc(partner["sourceId"])}</td>'
                                         f'<td>{esc(partner["pathwayRampId"])}</td>'
                                         f'<td>{partner["score"]:.3f}</td></tr>')
                        parts.append('</tbody></table></div>')
                    else:
                        parts.append('<p>No positive similarity partners are stored for this scope.</p>')
                    if blob['more_partners']:
                        parts.append(f'<p>{blob["more_partners"]:,} more stored partners.</p>')
        parts.append('</details>')
    parts.append('<h3>Matched exact-duplicate example · ' +
                 ' / '.join(esc(source) for source in PATHWAY_DUPLICATE_EXAMPLE) + '</h3>')
    for p in profiles:
        example = p['pathway_overlap']['duplicate_example']
        parts.append(f'<details><summary>{esc(p["label"])} · {esc(example["state"].replace("_", " "))}</summary>')
        if 'left' in example:
            left, right = example['left'], example['right']
            parts.append(f'<p>{esc(left[1])} · <code>{esc(left[0])}</code> · {esc(left[2])}<br>'
                         f'{esc(right[1])} · <code>{esc(right[0])}</code> · {esc(right[2])}</p>')
        parts.append('</details>')
    parts.append('</section>')
    parts.append('<section class="tab-panel" id="db_version" role="tabpanel" aria-labelledby="tab-db_version">'
                 '<h2>Database version and source overlaps</h2>'
                 '<p>The R package reads the four intersection JSON fields from the latest <code>db_version</code> row. '
                 'Each count is an exact, mutually exclusive combination of reporting sources. Pathway-mapped scopes '
                 'include only analytes linked to a non-SMPDB pathway. Changes compare totals with the preceding build; '
                 'generated intersection IDs are not matched across versions.</p>')
    parts.append('<div class="scroll"><table><thead><tr><th>Stored version field</th>' +
                 ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
    for field in ('ramp_version', 'load_timestamp', 'version_notes', 'db_sql_url'):
        parts.append(f'<tr><td><code>{field}</code></td>')
        for p in profiles:
            version = p['db_version']
            value = ('Table absent' if not version['present'] else
                     'Row absent' if version['latest'] is None else
                     'Column absent' if field not in version['columns'] else
                     'NULL' if version['latest'][field] is None else
                     version['latest'][field])
            parts.append(f'<td>{esc(value)}</td>')
        parts.append('</tr>')
    parts.append('</tbody></table></div>')
    parts.append('<h3>Exact source intersections</h3><div class="scroll"><table><thead><tr><th>Scope</th>' +
                 ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
    for field in DB_VERSION_FIELDS:
        parts.append(f'<tr><td>{esc(DB_VERSION_SCOPE_LABELS[field])}</td>')
        for index, p in enumerate(profiles):
            version = p['db_version']
            item = version['intersections'].get(field)
            if not version['present']:
                parts.append('<td>Table absent</td>')
                continue
            if item is None:
                parts.append('<td>Column absent</td>')
                continue
            state = item['state']
            if state not in ('valid', 'empty'):
                value = {'row_absent': 'Row absent', 'null': 'NULL cell',
                         'malformed': 'Malformed JSON: ' + item.get('error', '')}[state]
                parts.append(f'<td>{esc(value)}</td>')
                continue
            previous = profiles[index - 1]['db_version']['intersections'].get(field) if index else None
            delta = item['total'] - previous['total'] if previous and previous['state'] in ('valid', 'empty') else None
            pct = (100 * delta / previous['total'] if delta is not None and previous['total']
                   else 100 if delta is not None and item['total'] else None)
            label = 'Empty array · ' if state == 'empty' else ''
            unit = 'analyte' if item['total'] == 1 else 'analytes'
            parts.append(f'<td{change_style(delta, pct or 0)}>{label}{item["total"]:,} {unit}'
                         f'<small>{item["combinations"]:,} combinations · {len(item["sources"])} sources'
                         + (f' · {item["duplicate_combination_records"]} repeated-combination rows'
                            if item['duplicate_combination_records'] else '')
                         + (f' · {item["zero_size_records"]} zero-size rows' if item['zero_size_records'] else '')
                         + (f' · Δ {delta:+,}' if delta is not None else '')
                         + '</small><small>' + esc(', '.join(item['sources']) or 'No sources') + '</small></td>')
        parts.append('</tr>')
    parts.append('</tbody></table></div>')
    parts.append('<h3>UpSet plots by build</h3><p>Choose an analyte and pathway scope. '
                 'Bars show exact, mutually exclusive source combinations; filled dots show the sources '
                 'in each combination. The largest 20 combinations are plotted on one shared linear scale '
                 'within each scope. The numerical tables below give exact counts.</p>')
    parts.append('<div class="scope-choices" role="group" aria-label="UpSet plot scope">' +
                 ''.join(f'<button type="button" class="upset-scope-choice" '
                         f'data-scope="scope-{i}" aria-pressed="{str(i == 0).lower()}">'
                         f'{esc(DB_VERSION_SCOPE_LABELS[field])}</button>'
                         for i, field in enumerate(DB_VERSION_FIELDS)) + '</div>')
    for scope_index, field in enumerate(DB_VERSION_FIELDS):
        valid = [p['db_version']['intersections'].get(field, {}) for p in profiles]
        shared_max = max((row['size'] for item in valid if item.get('state') == 'valid'
                          for row in item['combination_data']), default=0)
        parts.append(f'<div class="upset-scope-panel" data-scope="scope-{scope_index}"'
                     + (' hidden' if scope_index else '') + '>')
        for p in profiles:
            item = p['db_version']['intersections'].get(field)
            parts.append(f'<article class="upset-card"><h4>{esc(p["label"])}</h4>')
            if not p['db_version']['present']:
                parts.append('<p>Table absent; no plot available.</p>')
            elif item is None:
                parts.append('<p>Column absent; no plot available.</p>')
            elif item['state'] == 'row_absent':
                parts.append('<p>Version row absent; no plot available.</p>')
            elif item['state'] == 'null':
                parts.append('<p>NULL intersection field; no plot available.</p>')
            elif item['state'] == 'malformed':
                parts.append(f'<p>Malformed intersection JSON: {esc(item["error"])}</p>')
            elif item['state'] == 'empty':
                parts.append('<p>Empty intersection array; no source combinations recorded.</p>')
            else:
                plotted_rows = min(20, item['combinations'])
                plotted = sum(row['size'] for row in sorted(
                    item['combination_data'], key=lambda row: (-row['size'], tuple(row['sets'])))[:20])
                other = item['total'] - plotted
                share = 100 * plotted / item['total'] if item['total'] else 0
                parts.append(f'<p>Showing {plotted_rows:,} of {item["combinations"]:,} combinations '
                             f'({plotted:,} of {item["total"]:,} analytes; {share:.1f}%). '
                             f'Other combinations: {item["combinations"] - plotted_rows:,} '
                             f'containing {other:,} analytes.</p>')
                parts.append('<div class="upset-scroll">' + render_upset_svg(
                    item, title=p['label'] + ' / ' + DB_VERSION_SCOPE_LABELS[field],
                    shared_max=shared_max) + '</div>')
            parts.append('</article>')
        parts.append('</div>')
    for field in DB_VERSION_FIELDS:
        combinations = sorted({name for p in profiles
                               for name in p['db_version']['intersections'].get(field, {}).get('by_sets', {})},
                              key=lambda name: (-max((p['db_version']['intersections'].get(field, {})
                                                      .get('by_sets', {}).get(name, 0) for p in profiles)), name))[:30]
        parts.append(f'<details><summary>{esc(DB_VERSION_SCOPE_LABELS[field])} · top 30 source combinations</summary>')
        if combinations:
            parts.append('<div class="scroll"><table><thead><tr><th>Source combination</th>' +
                         ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
            for name in combinations:
                parts.append(f'<tr><td>{esc(name)}</td>')
                for index, p in enumerate(profiles):
                    intersection = p['db_version']['intersections'].get(field, {})
                    current = intersection.get('by_sets', {}).get(name)
                    previous = (profiles[index - 1]['db_version']['intersections'].get(field, {})
                                .get('by_sets', {}).get(name) if index else None)
                    delta = current - previous if current is not None and previous is not None else None
                    parts.append(f'<td{change_style(delta, 100 * delta / previous if delta is not None and previous else 100 if delta else 0)}>' +
                                 (f'{current:,}' if current is not None else 'Absent') +
                                 (f'<small>Δ {delta:+,}</small>' if delta is not None else '') + '</td>')
                parts.append('</tr>')
            parts.append('</tbody></table></div>')
        else:
            parts.append('<p>No stored source combinations.</p>')
        parts.append('</details>')
    parts.append('<h3>All stored db_version rows</h3>')
    for p in profiles:
        version = p['db_version']
        parts.append(f'<details><summary>{esc(p["label"])} · ' +
                     ('table absent' if not version['present'] else f'{len(version["rows"]):,} rows') + '</summary>')
        for row in version['rows']:
            parts.append(f'<details><summary>{esc(row.get("ramp_version"))} · {esc(row.get("load_timestamp"))}</summary>')
            for field in version['columns']:
                value = row[field]
                if field in DB_VERSION_FIELDS and value is not None:
                    parts.append(f'<details class="db-version-raw"><summary><code>{esc(field)}</code> '
                                 f'· {len(value):,} characters</summary><pre>{esc(value)}</pre></details>')
                else:
                    parts.append(f'<p><code>{esc(field)}</code>: {esc("NULL" if value is None else value)}</p>')
            parts.append('</details>')
        parts.append('</details>')
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
  const scopeButtons = Array.from(document.querySelectorAll('.upset-scope-choice'));
  const scopePanels = Array.from(document.querySelectorAll('.upset-scope-panel'));
  scopeButtons.forEach(button => button.addEventListener('click', () => {
    const selected = button.dataset.scope;
    scopeButtons.forEach(item => item.setAttribute('aria-pressed', String(item === button)));
    scopePanels.forEach(panel => { panel.hidden = panel.dataset.scope !== selected; });
  }));
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
    parser = argparse.ArgumentParser(description='Compare RaMP SQLite tables across builds')
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
