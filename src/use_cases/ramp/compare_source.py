"""Read-only comparison of the RaMP source lookup table across SQLite builds."""
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
        manifest = {}
        if db.execute("SELECT 1 FROM sqlite_master WHERE type='table' AND name='ramp_export_metadata'").fetchone():
            row = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
            if row:
                manifest = json.loads(row[0])
        table_names = {name for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")}
        scope = manifest.get('scope') or ('source_table_diagnostic' if table_names == {'analyte', 'source'}
                                          else 'historical full database')
        return {'label': label, 'path': str(path), 'scope': scope,
                'ramp_ids': ramp_ids, 'examples': examples, 'lookup_quality': lookup,
                'analyte_field_coverage': {kind: values for (kind,), values in analyte_coverage.items()},
                'source_field_coverage': {
                    kind: {provider or '(missing provider)': values
                           for (provider, row_kind), values in source_field_counts.items()
                           if row_kind == kind}
                    for kind in ('compound', 'gene')},
                **source_coverage,
                'source_rows': db.execute('SELECT count(*) FROM source').fetchone()[0]}
    finally:
        db.close()


def compare(databases):
    report = {'format_version': 1, 'generated_at': datetime.now(timezone.utc).isoformat(),
            'definition': 'source records identifiers each input used to report retained RaMP data, including associations and chemistry',
            'databases': [profile(label, path) for label, path in databases]}
    report['source_examples'] = source_examples(report['databases'])
    return report


def source_examples(profiles):
    """Choose one uniquely mapped ID per provider and compare its stored rows."""
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
            latest = max(i for i, profile in enumerate(profiles)
                         if provider in profile['metabolite_sources'] or provider in profile['gene_sources'])
            section, kind = (('metabolite_sources', 'compound')
                             if profiles[latest]['metabolite_sources'].get(provider) else
                             ('gene_sources', 'gene'))
            available = sum(profile[section].get(provider, 0) > 0 for profile in profiles)
            candidates = dbs[latest].execute('''
                SELECT DISTINCT sourceId FROM source WHERE dataSource=? AND geneOrCompound=?
                AND sourceId IS NOT NULL ORDER BY sourceId LIMIT 1000''', (provider, kind)).fetchall()
            if not candidates:
                # Existing databases use exact labels; this fallback preserves
                # the trimmed provider grouping for atypical stored spelling.
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
                # Search farther only if the first page contains no unique row.
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
            examples[provider] = {'source_id': source_id, 'entity_type': kind,
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
                            ('gene_sources', (None,)), ('lookup_quality', (None,))]:
        keys = sorted({key for p in profiles for key in p[section]})
        if section == 'ramp_ids':
            keys = ['Metabolite RaMP IDs', 'Gene/protein RaMP IDs']
        result[section] = []
        for key in keys:
            for field in fields:
                values = [p[section].get(key, {}).get(field, 0) if field else p[section].get(key)
                          for p in profiles]
                result[section].append({'metric': f'{key} / {field.replace("_", " ")}' if field else key,
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
             ('lookup_quality', 'Lookup quality')]
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
             '<title>RaMP source-table comparison</title>',
             '<style>body{font:16px/1.5 system-ui,sans-serif;margin:0;background:#f5f7fa;color:#182635}main{max-width:1500px;margin:auto;padding:24px}p{max-width:100ch}.scroll{overflow-x:auto;background:white;border:1px solid #d2dbe3}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px 12px;border-bottom:1px solid #dde4eb;text-align:right;min-width:140px}th:first-child,td:first-child{text-align:left;min-width:300px}th{background:#eaf0f6}tbody tr:nth-child(even){background:#f8fafc}tbody tr.coverage{background:#eef3f7;color:#42566a;font-size:12px}tbody tr.coverage td{padding-top:4px;padding-bottom:8px}small{display:block;color:#35475a}article{background:white;border-left:4px solid #52738c;padding:10px 16px;margin:10px 0}nav a{margin-right:20px;color:#23527c}.example th,.example td{text-align:left}.example th:first-child,.example td:first-child{min-width:220px}@media(max-width:600px){main{padding:12px}}</style><main>',
             '<h1>RaMP source-table comparison</h1>',
             '<p>The first table counts distinct RaMP IDs in <code>analyte</code>. The examples below match the same records across builds by a source ID, since RaMP IDs can change.</p>',
             '<p>Each change compares with the preceding build. Fresh RaMP IDs are counted within a build and are never matched across builds.</p>',
             f'<p>Generated {esc(report["generated_at"])}.</p>']
    for p in profiles:
        parts.append(f'<article><strong>{esc(p["label"])}</strong> · {p["source_rows"]:,} source rows · {esc(p["scope"])}<br><code>{esc(p["path"])}</code></article>')
    parts.append('<nav>' + ''.join(f'<a href="#s{i}">{esc(title)}</a>' for i, (_, title) in enumerate(names)) + '</nav>')
    for i, (key, title) in enumerate(names):
        parts.append(f'<h2 id="s{i}">{esc(title)}</h2>')
        if key == 'ramp_ids':
            parts.append('<p>Metabolites and genes/proteins are counted separately. These are within-build analyte counts, not stable IDs across releases.</p>')
        if key in ('metabolite_sources', 'gene_sources'):
            parts.append('<p>Each cell counts distinct RaMP IDs with at least one ID attributed to that input in <code>source</code>. An analyte may occur under several inputs, so the rows do not add up to the analyte total. KEGG IDs supplied by HMDB or WikiPathways retain their separate legacy source labels. Changes describe coverage, not quality.</p>')
        if key == 'source_examples':
            parts.append('<p>Each example follows one stored source ID and analyte type across builds. Field coverage below each example counts populated values across all rows for that input and analyte type; NULL, blanks, and known placeholders are unpopulated. Expand an input to compare all eight <code>source</code> columns.</p>')
            for provider, example in report['source_examples'].items():
                parts.append(f'<details><summary>{esc(provider)} · {esc(example["entity_type"])} · <code>{esc(example["source_id"])}</code></summary>')
                parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in SOURCE_FIELDS) + '</tr></thead><tbody>')
                for profile, row in zip(profiles, example['rows']):
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
                    coverage = profile['source_field_coverage'].get(example['entity_type'], {}).get(provider)
                    parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                    parts.extend(f'<td>{esc(coverage_cell(coverage, field))}</td>' for field in SOURCE_FIELDS)
                    parts.append('</tr>')
                parts.append('</tbody></table></div></details>')
            continue
        if key == 'lookup_quality':
            parts.append('<p>A reported ID can map to several gene/protein RaMP IDs. These are candidate lookups, not guaranteed unique resolutions.</p>')
        parts.append('<div class="scroll"><table><thead><tr><th>Metric</th>' + ''.join(f'<th>{esc(p["label"])}</th>' for p in profiles) + '</tr></thead><tbody>')
        for row in sections[key]:
            parts.append(f'<tr><td>{esc(row["metric"])}</td>')
            for datum in row['cells']:
                value, delta, percent = datum['value'], datum['delta'], datum['percent_change']
                style = ''
                if percent is not None and math.isfinite(percent) and percent:
                    endpoint = (173, 222, 188) if percent > 0 else (246, 185, 185)
                    strength = min(abs(percent), 100) / 100
                    rgb = ','.join(str(round(255 + (x-255)*strength)) for x in endpoint)
                    style = f' style="background-color:rgb({rgb})"'
                parts.append(f'<td{style}>{value:,}' if value is not None else f'<td{style}>Unavailable')
                if delta is not None:
                    parts.append(f'<small>Δ {delta:+,} ({percent:+.1f}%)</small>' if percent is not None else f'<small>Δ {delta:+,} (—)</small>')
                parts.append('</td>')
            parts.append('</tr>')
        parts.append('</tbody></table></div>')
        if key == 'ramp_ids':
            for example_key, title, source_id, _ in EXAMPLES:
                parts.append(f'<h3>{esc(title)}</h3><p>Matched by <code>{esc(source_id)}</code> in each database’s <code>source</code> table. Each example comes from one <code>analyte</code> row; field coverage counts populated values across all analytes of that type. NULL, blanks, and known placeholders are unpopulated.</p>')
                parts.append('<div class="scroll"><table class="example"><thead><tr><th>Database</th>' +
                             ''.join(f'<th>{esc(field)}</th>' for field in EXAMPLE_FIELDS) + '</tr></thead><tbody>')
                for profile in profiles:
                    parts.append(f'<tr><td>{esc(profile["label"])}</td>')
                    for field in EXAMPLE_FIELDS:
                        datum = profile['examples'][example_key]['cells'][field]
                        display = ('Column absent' if datum['state'] == 'column_absent' else
                                   'NULL' if datum['state'] == 'null' else str(datum['value']))
                        parts.append(f'<td>{esc(display)}</td>')
                    parts.append('</tr>')
                    coverage = profile['analyte_field_coverage'].get('compound' if example_key == 'metabolite' else 'gene')
                    parts.append(f'<tr class="coverage"><td>{esc(profile["label"])} · Field coverage</td>')
                    parts.extend(f'<td>{esc(coverage_cell(coverage, field))}</td>' for field in EXAMPLE_FIELDS)
                    parts.append('</tr>')
                parts.append('</tbody></table></div>')
    parts.append('</main></html>')
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
    parser = argparse.ArgumentParser(description='Compare only RaMP source lookup tables across SQLite builds')
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
