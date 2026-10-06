import json
import sqlite3

from src.use_cases.ramp.compare_source import compare, tables, write


def make_database(path, rows, *, names=False):
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE analyte(rampId TEXT PRIMARY KEY,type TEXT);
            CREATE TABLE source(sourceId TEXT,rampId TEXT,IDtype TEXT,
                                geneOrCompound TEXT,commonName TEXT,
                                priorityHMDBStatus TEXT,dataSource TEXT,pathwayCount INTEGER);
        ''')
        if names:
            db.execute('ALTER TABLE analyte ADD COLUMN common_name TEXT')
        for rid in {row[1] for row in rows}:
            kind = 'compound' if rid.startswith('C') else 'gene'
            if names:
                db.execute('INSERT INTO analyte VALUES (?,?,?)',
                           (rid, kind, 'D-Glucose' if rid == 'C1' else 'EGFR' if rid == 'G1' else None))
            else:
                db.execute('INSERT INTO analyte VALUES (?,?)', (rid, kind))
        db.executemany('''INSERT INTO source
            (sourceId,rampId,IDtype,geneOrCompound,dataSource) VALUES (?,?,?,?,?)''', rows)


def test_focused_source_report_counts_and_previous_build_deltas(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    make_database(old, [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
                        ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')])
    make_database(new, [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
                        ('chebi:2', 'C2', 'chebi', 'compound', 'pfocr'),
                        ('chebi:2', 'C2', 'chebi', 'compound', 'chebi'),
                        ('chebi:3', 'C2', 'chebi', 'compound', 'chebi'),
                        ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')], names=True)
    with sqlite3.connect(new) as db:
        db.execute("INSERT INTO analyte VALUES ('C3','compound',NULL)")
    report = compare([('Old', old), ('New', new)])
    assert report['databases'][1]['scope'] == 'source_table_diagnostic'
    assert report['databases'][1]['ramp_ids'] == {'Metabolite RaMP IDs': 3, 'Gene/protein RaMP IDs': 1}
    assert report['databases'][1]['metabolite_sources']['chebi'] == 1
    assert report['databases'][1]['metabolite_sources']['pfocr'] == 1
    assert report['databases'][1]['gene_sources']['reactome'] == 1
    assert report['source_examples']['chebi']['rows'][0]['state'] == 'provider_absent'
    assert report['source_examples']['chebi']['rows'][1]['cells']['rampId'] == 'C2'
    assert report['source_examples']['hmdb']['rows'][0]['state'] == 'available'
    assert report['databases'][0]['analyte_field_coverage']['compound']['fields']['common_name'] is None
    assert report['databases'][1]['analyte_field_coverage']['compound']['fields']['common_name'] == 1
    assert report['databases'][1]['analyte_field_coverage']['compound']['total'] == 3
    assert report['databases'][1]['source_field_coverage']['compound']['chebi']['fields']['commonName'] == 0
    rows = tables(report)['ramp_ids']
    assert [row['metric'] for row in rows] == ['Metabolite RaMP IDs', 'Gene/protein RaMP IDs']
    row = rows[0]
    assert [c['value'] for c in row['cells']] == [1, 3]
    assert row['cells'][1]['delta'] == 2
    html_path, json_path = write(report, tmp_path / 'source.html')
    page = html_path.read_text()
    assert 'RaMP source-table comparison' in page
    assert "RaMP IDs represented in &#x27;analyte&#x27;" in page
    assert 'hmdb:HMDB0000122' in page
    assert 'uniprot:P00533' in page
    assert 'Column absent' in page
    assert report['databases'][0]['examples']['metabolite']['cells']['common_name']['state'] == 'column_absent'
    assert report['databases'][1]['examples']['gene']['cells']['common_name']['value'] == 'EGFR'
    assert 'Reported IDs by input and analyte type' not in page
    assert 'Reported ID namespaces' not in page
    assert 'Metabolites by input' in page
    assert 'Genes/proteins by input' in page
    assert 'Example source rows' in page
    assert 'priorityHMDBStatus' in page
    assert 'Field coverage' in page
    assert '1 / 3 (33.3%)' in page
    assert '0 / 2 (0.0%)' in page
    assert 'Input absent' in page
    assert 'Chemistry completeness' not in page
    assert set(json.loads(json_path.read_text())['comparison_tables']) == {
        'ramp_ids', 'metabolite_sources', 'gene_sources', 'lookup_quality'}


def test_missing_provider_is_counted_without_becoming_an_example(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    base = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')]
    make_database(old, base)
    make_database(new, base + [('chebi:2', 'C2', 'chebi', 'compound', 'chebi')])
    with sqlite3.connect(new) as db:
        db.executemany('''INSERT INTO source
            (sourceId,rampId,IDtype,geneOrCompound,dataSource) VALUES (?,?,?,?,?)''', [
            ('x:1', 'C1', 'x', 'compound', None),
            ('x:2', 'C2', 'x', 'compound', ''),
        ])
    report = compare([('Old', old), ('New', new)])
    assert report['databases'][1]['metabolite_sources']['(missing provider)'] == 2
    assert '(missing provider)' not in report['source_examples']


def test_example_selection_scans_past_ambiguous_first_page(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    base = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')]
    make_database(old, base)
    make_database(new, base + [('chebi:2', 'C2', 'chebi', 'compound', 'chebi')])
    ambiguous = [(f'rhea:{i:03d}', rid, 'rhea', 'compound', 'rhea')
                 for i in range(250) for rid in ('C1', 'C2')]
    with sqlite3.connect(new) as db:
        db.executemany('''INSERT INTO source
            (sourceId,rampId,IDtype,geneOrCompound,dataSource) VALUES (?,?,?,?,?)''',
                       ambiguous + [('rhea:999', 'C1', 'rhea', 'compound', 'rhea')])
    report = compare([('Old', old), ('New', new)])
    assert report['source_examples']['rhea']['source_id'] == 'rhea:999'


def test_all_ambiguous_provider_still_has_explicit_example(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    base = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')]
    make_database(old, base)
    make_database(new, base + [('chebi:2', 'C2', 'chebi', 'compound', 'chebi'),
                               ('rhea:1', 'C1', 'rhea', 'compound', 'rhea'),
                               ('rhea:1', 'C2', 'rhea', 'compound', 'rhea')])
    report = compare([('Old', old), ('New', new)])
    example = report['source_examples']['rhea']['rows'][1]
    assert example['state'] == 'multiple_mappings'
    assert example['matching_rows'] == 2
    assert 'Showing 1 of 2 mappings' in write(report, tmp_path / 'report.html')[0].read_text()


def test_field_coverage_excludes_placeholders_and_counts_zero(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('hmdb:HMDB0000123', 'C2', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')]
    make_database(old, rows, names=True)
    make_database(new, rows, names=True)
    with sqlite3.connect(new) as db:
        db.execute("UPDATE analyte SET common_name='NA' WHERE rampId='C2'")
        db.execute("UPDATE source SET commonName='NA', priorityHMDBStatus='no_HMDB_status', pathwayCount=-1 WHERE rampId='C1'")
        db.execute("UPDATE source SET commonName='Water', priorityHMDBStatus='quantified', pathwayCount=0 WHERE rampId='C2'")
    report = compare([('Old', old), ('New', new)])
    analyte = report['databases'][1]['analyte_field_coverage']['compound']
    source = report['databases'][1]['source_field_coverage']['compound']['hmdb']
    assert (analyte['total'], analyte['fields']['common_name']) == (2, 1)
    assert source['total'] == 2
    assert source['fields']['commonName'] == 1
    assert source['fields']['priorityHMDBStatus'] == 1
    assert source['fields']['pathwayCount'] == 1


def test_optional_source_columns_can_be_absent(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')]
    make_database(old, rows, names=True)
    make_database(new, rows, names=True)
    with sqlite3.connect(old) as db:
        db.executescript('''
            CREATE TABLE source_minimal(sourceId TEXT,rampId TEXT,IDtype TEXT,
                                        geneOrCompound TEXT,dataSource TEXT);
            INSERT INTO source_minimal SELECT sourceId,rampId,IDtype,geneOrCompound,dataSource FROM source;
            DROP TABLE source;
            ALTER TABLE source_minimal RENAME TO source;
        ''')
    report = compare([('Old', old), ('New', new)])
    assert report['databases'][0]['source_field_coverage']['compound']['hmdb']['fields']['commonName'] is None
    assert report['source_examples']['hmdb']['rows'][0]['cells']['commonName'] == {'state': 'column_absent'}
    assert 'Column absent' in write(report, tmp_path / 'source.html')[0].read_text()
