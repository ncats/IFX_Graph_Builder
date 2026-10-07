import json
import sqlite3

from src.use_cases.ramp.compare_source import compare, tables, write


def make_database(path, rows, *, names=False, synonyms=None):
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
        if synonyms is not None:
            db.execute('''CREATE TABLE analytesynonym
                          (Synonym TEXT COLLATE NOCASE,rampId TEXT,geneOrCompound TEXT,source TEXT)''')
            db.executemany('INSERT INTO analytesynonym VALUES (?,?,?,?)', synonyms)


def add_associations(path, *, compound_id, gene_id, include_gene_pathway=True):
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE ontology(rampOntologyId TEXT,commonName TEXT,HMDBOntologyType TEXT,metCount INTEGER);
            CREATE TABLE analytehasontology(rampCompoundId TEXT,rampOntologyId TEXT);
            CREATE TABLE pathway(pathwayRampId TEXT,sourceId TEXT,type TEXT,pathwayCategory TEXT,pathwayName TEXT);
            CREATE TABLE analytehaspathway(rampId TEXT,pathwayRampId TEXT,pathwaySource TEXT);
        ''')
        db.execute("INSERT INTO ontology VALUES ('OL1','Blood','Location',1)")
        db.execute('INSERT INTO analytehasontology VALUES (?,?)', (compound_id, 'OL1'))
        db.execute("INSERT INTO pathway VALUES ('P1','WP1','wiki',NULL,'Example pathway')")
        db.execute('INSERT INTO analytehaspathway VALUES (?,?,?)', (compound_id, 'P1', 'wiki'))
        if include_gene_pathway:
            db.execute('INSERT INTO analytehaspathway VALUES (?,?,?)', (gene_id, 'P1', 'wiki'))


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
    assert report['source_examples']['chebi']['compound']['rows'][0]['state'] == 'provider_absent'
    assert report['source_examples']['chebi']['compound']['rows'][1]['cells']['rampId'] == 'C2'
    assert report['source_examples']['chebi']['gene'] is None
    assert report['source_examples']['hmdb']['compound']['rows'][0]['state'] == 'available'
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
    assert 'RaMP lookup and association comparison' in page
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
    assert 'No gene/protein rows for this input' in page
    assert 'Δ +1 (from zero)' in page
    assert 'background-color:rgb(173,222,188)' in page
    assert 'Chemistry completeness' not in page
    assert page.count('<button type="button" role="tab"') == 6
    assert page.count('<section class="tab-panel"') == 6
    assert page.index('id="analyte" role="tabpanel"') < page.index('id="source" role="tabpanel"')
    assert page.index('id="source" role="tabpanel"') < page.index('id="analytesynonym" role="tabpanel"')
    assert page.index('id="analytesynonym" role="tabpanel"') < page.index('id="analytehasontology" role="tabpanel"')
    assert page.index('id="analytehasontology" role="tabpanel"') < page.index('id="analytehaspathway" role="tabpanel"')
    assert page.index('id="analytehaspathway" role="tabpanel"') < page.index('id="about" role="tabpanel"')
    assert page.index('Example source rows') < page.index('id="analytesynonym" role="tabpanel"')
    assert page.index('Generated ') > page.index('id="about" role="tabpanel"')
    assert "history.pushState(null, '', '#' + id)" in page
    assert '.tab-panel[hidden]{display:block}' in page
    assert set(json.loads(json_path.read_text())['comparison_tables']) == {
        'ramp_ids', 'metabolite_sources', 'gene_sources',
        'compound_synonyms', 'gene_synonyms', 'ontology_associations',
        'pathway_associations', 'lookup_quality'}


def test_synonym_report_counts_sources_and_shows_real_examples(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    source_rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
                   ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')]
    make_database(old, source_rows, synonyms=[
        ('Glucose', 'C1', 'compound', 'hmdb'),
        ('Dextrose', 'C1', 'compound', 'hmdb'),
        ('EGFR', 'G1', 'gene', 'reactome'),
    ])
    make_database(new, source_rows, synonyms=[
        ('Dextrose', 'C1', 'compound', 'hmdb'),
        ('EGFR', 'G1', 'gene', 'uniprot'),
    ])
    report = compare([('Old', old), ('New', new)])
    compound = next(row for row in tables(report)['compound_synonyms']
                    if row['metric'] == 'hmdb / rows')
    assert [datum['value'] for datum in compound['cells']] == [2, 1]
    assert compound['cells'][1]['delta'] == -1
    assert report['databases'][0]['synonym_sources']['gene']['reactome'] == {'rows': 1, 'ramp_ids': 1}
    assert report['databases'][1]['synonym_sources']['gene']['uniprot'] == {'rows': 1, 'ramp_ids': 1}
    page = write(report, tmp_path / 'synonyms.html')[0].read_text()
    assert 'Metabolite synonyms by attributed source' in page
    assert 'Gene/protein synonyms by attributed source' in page
    assert 'Example synonym rows' in page
    assert 'Glucose' in page and 'Dextrose' in page and 'EGFR' in page
    assert 'Identity anchor:' in page
    assert 'No synonyms attributed to this source' in page
    assert '<summary>hmdb</summary>' in page
    assert 'No gene/protein synonyms attributed to this input' in page


def test_source_examples_pair_metabolite_and_gene_for_same_input(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    for path, compound, gene in ((old, 'C1', 'G1'), (new, 'C8', 'G8')):
        make_database(path, [
            ('hmdb:HMDB0000122', compound, 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', gene, 'uniprot', 'gene', 'hmdb')])
    report = compare([('Old', old), ('New', new)])
    assert report['source_examples']['hmdb']['compound']['source_id'] == 'hmdb:HMDB0000122'
    assert report['source_examples']['hmdb']['gene']['source_id'] == 'uniprot:P00533'
    assert [row['cells']['rampId'] for row in report['source_examples']['hmdb']['compound']['rows']] == ['C1', 'C8']
    assert [row['cells']['rampId'] for row in report['source_examples']['hmdb']['gene']['rows']] == ['G1', 'G8']
    page = write(report, tmp_path / 'both.html')[0].read_text()
    assert '<summary>hmdb</summary>' in page
    assert '<summary>Metabolite · <code>hmdb:HMDB0000122</code></summary>' in page
    assert '<summary>Gene/protein · <code>uniprot:P00533</code></summary>' in page


def test_association_tabs_compare_sources_and_match_examples_across_ramp_ids(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    for path, compound_id, gene_id, gene_pathway in (
            (old, 'C1', 'G1', False), (new, 'C8', 'G8', True)):
        make_database(path, [
            ('hmdb:HMDB0000122', compound_id, 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', gene_id, 'uniprot', 'gene', 'uniprot')])
        add_associations(path, compound_id=compound_id, gene_id=gene_id,
                         include_gene_pathway=gene_pathway)
    report = compare([('Old', old), ('New', new)])
    assert report['databases'][1]['ontology_associations']['HMDB ontology / Location'] == {
        'rows': 1, 'analytes': 1, 'targets': 1}
    assert report['databases'][0]['pathway_associations']['Genes/proteins / All sources']['rows'] == 0
    assert report['databases'][1]['pathway_associations']['Genes/proteins / wiki']['rows'] == 1
    gene_rows = next(row['cells'] for row in tables(report)['pathway_associations']
                     if row['metric'] == 'Genes/proteins / wiki / rows')
    assert [cell['value'] for cell in gene_rows] == [0, 1]
    assert gene_rows[1]['delta'] == 1
    assert [row['ramp_id'] for row in report['ontology_example']['rows']] == ['C1', 'C8']
    gene_example = report['pathway_examples']['Genes/proteins / wiki']
    assert [row['state'] for row in gene_example['rows']] == ['association_absent', 'available']
    assert [row['ramp_id'] for row in gene_example['rows']] == ['G1', 'G8']
    page = write(report, tmp_path / 'associations.html')[0].read_text()
    assert 'Metabolite ontology associations by HMDB type' in page
    assert 'Pathway associations by source and analyte type' in page
    assert 'Example ontology association' in page
    assert 'Example pathway associations' in page
    assert '<summary>wiki</summary>' in page
    assert 'Gene/protein · <code>uniprot:P00533</code>' in page
    assert 'Not applicable: ontology links only metabolites' in page
    assert 'Old WikiPathways gene symbols came from WikiPathways itself' in page
    assert 'current exporter reads them from the patched graph' in page
    assert 'Δ +1 (from zero)' in page


def test_synonym_report_distinguishes_missing_table_from_zero(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')]
    make_database(old, rows)
    make_database(new, rows, synonyms=[('Glucose', 'C1', 'compound', 'hmdb')])
    report = compare([('Old', old), ('New', new)])
    values = next(row for row in tables(report)['compound_synonyms']
                  if row['metric'] == 'hmdb / rows')['cells']
    assert [datum['value'] for datum in values] == [None, 1]
    assert values[1]['delta'] is None
    page = write(report, tmp_path / 'missing.html')[0].read_text()
    assert 'Table absent' in page


def test_synonym_examples_follow_same_source_id_across_changed_ramp_ids(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    make_database(old, [('hmdb:1', 'C1', 'hmdb', 'compound', 'hmdb'),
                        ('hmdb:2', 'C2', 'hmdb', 'compound', 'hmdb'),
                        ('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
                        ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')],
                  synonyms=[('unrelated old', 'C1', 'compound', 'hmdb'),
                            ('same entity old', 'C2', 'compound', 'hmdb')])
    make_database(new, [('hmdb:1', 'C8', 'hmdb', 'compound', 'hmdb'),
                        ('hmdb:2', 'C9', 'hmdb', 'compound', 'hmdb'),
                        ('hmdb:HMDB0000122', 'C8', 'hmdb', 'compound', 'hmdb'),
                        ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')],
                  synonyms=[('same entity new', 'C9', 'compound', 'hmdb'),
                            (None, 'C8', 'compound', 'hmdb'),
                            ('unrelated new', 'C8', 'compound', 'hmdb')])
    example = compare([('Old', old), ('New', new)])['synonym_examples']['compound']['hmdb']
    assert example['source_id'] == 'hmdb:HMDB0000122'
    assert [row['ramp_id'] for row in example['rows']] == ['C1', 'C8']
    assert [row['synonyms'] for row in example['rows']] == [['unrelated old'], ['unrelated new', None]]


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
    assert report['source_examples']['rhea']['compound']['source_id'] == 'rhea:999'


def test_all_ambiguous_provider_still_has_explicit_example(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    base = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'reactome')]
    make_database(old, base)
    make_database(new, base + [('chebi:2', 'C2', 'chebi', 'compound', 'chebi'),
                               ('rhea:1', 'C1', 'rhea', 'compound', 'rhea'),
                               ('rhea:1', 'C2', 'rhea', 'compound', 'rhea')])
    report = compare([('Old', old), ('New', new)])
    example = report['source_examples']['rhea']['compound']['rows'][1]
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
    with sqlite3.connect(old) as db:
        db.execute("UPDATE source SET commonName='Old source name' WHERE geneOrCompound='compound'")
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
    page = write(report, tmp_path / 'coverage.html')[0].read_text()
    assert 'Δ -50.0 pp' in page
    assert 'Δ +50.0 pp' in page
    assert 'Green = increased; red = decreased' in page


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
    assert report['source_examples']['hmdb']['compound']['rows'][0]['cells']['commonName'] == {'state': 'column_absent'}
    assert 'Column absent' in write(report, tmp_path / 'source.html')[0].read_text()
