import json
import sqlite3

from src.use_cases.ramp.compare_source import (
    compare, db_version_profile, render_upset_svg, tables, write,
)


def test_upset_svg_limits_combinations_and_labels_each_bar():
    rows = [{'sets': [f'S{i}'], 'size': 25 - i} for i in range(25)]
    plot = render_upset_svg({
        'combination_data': rows, 'sources': [f'S{i}' for i in range(25)],
        'total': sum(row['size'] for row in rows), 'combinations': 25,
    }, title='Fixture / metabolites', shared_max=25)
    assert plot.count('class="upset-bar"') == 20
    assert '<title>S0: 25 analytes</title>' in plot
    assert '15 analytes are in combinations outside this plot' in plot


def test_db_version_profile_keeps_all_rows_and_reports_intersection_states():
    with sqlite3.connect(':memory:') as db:
        db.executescript('''
            CREATE TABLE db_version(
                ramp_version TEXT, load_timestamp TEXT, version_notes TEXT,
                met_intersects_json TEXT, gene_intersects_json TEXT,
                met_intersects_json_pw_mapped TEXT, gene_intersects_json_pw_mapped TEXT,
                db_sql_url TEXT);
            INSERT INTO db_version VALUES ('old','2025-01-01',NULL,NULL,NULL,NULL,NULL,NULL);
        ''')
        db.execute('INSERT INTO db_version VALUES (?,?,?,?,?,?,?,?)', (
            'new', '2026-01-01', 'test',
            json.dumps([{'id': 'x', 'sets': ['HMDB', 'KEGG'], 'size': 3},
                        {'id': 'y', 'sets': ['HMDB'], 'size': 0},
                        {'id': 'z', 'sets': ['KEGG', 'HMDB'], 'size': 2}]),
            '[]', '{bad json', None, None))
        result = db_version_profile(db, 'Fixture')
    assert len(result['rows']) == 2
    assert result['latest']['ramp_version'] == 'new'
    mets = result['intersections']['met_intersects_json']
    assert (mets['total'], mets['combinations'], mets['zero_size_records']) == (5, 2, 1)
    assert (mets['records'], mets['duplicate_combination_records']) == (3, 1)
    assert mets['by_sets']['HMDB & KEGG'] == 5
    assert mets['sources'] == ['HMDB', 'KEGG']
    assert result['intersections']['gene_intersects_json']['state'] == 'empty'
    assert result['intersections']['met_intersects_json_pw_mapped']['state'] == 'malformed'
    assert result['intersections']['gene_intersects_json_pw_mapped']['state'] == 'null'


def test_db_version_report_shows_populated_overlap_and_all_raw_rows(tmp_path):
    path = tmp_path / 'current.sqlite'
    make_database(path, [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
                         ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')])
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE db_version(
                ramp_version TEXT, load_timestamp TEXT, version_notes TEXT,
                met_intersects_json TEXT, gene_intersects_json TEXT,
                met_intersects_json_pw_mapped TEXT, gene_intersects_json_pw_mapped TEXT,
                db_sql_url TEXT);
            INSERT INTO db_version VALUES
                ('older','2025-01-01',NULL,NULL,NULL,NULL,NULL,NULL);
        ''')
        db.execute('INSERT INTO db_version VALUES (?,?,?,?,?,?,?,?)', (
            'new', '2026-01-01', 'test',
            json.dumps([{'id': 'cmpd_src_set_1', 'sets': ['HMDB'], 'size': 1}]),
            json.dumps([{'id': 'gene_src_set_1', 'sets': ['UniProt'], 'size': 1}]),
            '[]', '[]', None))
    report = compare([('Current', path)])
    html_path, _ = write(report, tmp_path / 'comparison.html')
    page = html_path.read_text()
    assert 'id="db_version" role="tabpanel"' in page
    assert 'Version &amp; UpSet' in page
    assert page.count('class="upset-scope-choice"') == 4
    assert page.count('class="upset-svg"') == 2
    assert 'filled dots show the sources' in page
    assert '1 analyte' in page
    assert 'HMDB' in page and 'UniProt' in page
    assert 'older · 2025-01-01' in page
    assert 'cmpd_src_set_1' in page


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


def add_catalyzed_classes(path, *, compound_id, gene_id, protein_type, class_name):
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE catalyzed(rampCompoundId TEXT,rampGeneId TEXT,proteinType TEXT);
            CREATE TABLE metabolite_class(ramp_id TEXT,class_source_id TEXT,
                class_level_name TEXT,class_name TEXT,source TEXT);
        ''')
        db.execute('INSERT INTO catalyzed VALUES (?,?,?)',
                   (compound_id, gene_id, protein_type))
        db.execute('INSERT INTO metabolite_class VALUES (?,?,?,?,?)',
                   (compound_id, 'hmdb:HMDB0000122', 'ClassyFire_class', class_name, 'hmdb'))
        db.execute('INSERT INTO metabolite_class VALUES (?,?,?,?,?)',
                   (compound_id, 'lipidmaps:L1', 'LipidMaps_category', 'Fatty Acyls', 'lipidmaps'))


def add_reactions(path, *, reaction_id, compound_id, gene_id, modern):
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE reaction(ramp_rxn_id TEXT,rxn_source_id TEXT,status INTEGER,
                is_transport INTEGER,direction TEXT,label TEXT,equation TEXT,
                html_equation TEXT,ec_num TEXT,has_human_prot INTEGER,only_human_mets INTEGER);
            CREATE TABLE reaction2met(ramp_rxn_id TEXT,rxn_source_id TEXT,ramp_cmpd_id TEXT,
                substrate_product INTEGER,met_source_id TEXT,met_name TEXT,is_cofactor INTEGER);
        ''')
        db.execute('''CREATE TABLE reaction2protein(ramp_rxn_id TEXT,rxn_source_id TEXT,
            ramp_gene_id TEXT,uniprot TEXT,protein_name TEXT%s)''' %
            (',is_reviewed INTEGER' if modern else ''))
        db.execute('INSERT INTO reaction VALUES (?,?,?,?,?,?,?,?,?,?,?)',
                   (reaction_id, 'rhea:10000', 1, 0, 'UN', 'Water reaction',
                    'A + H2O = B', '<b>A</b> + H2O = B', '1.1.1.1', -1, -1))
        db.execute('INSERT INTO reaction2met VALUES (?,?,?,?,?,?,?)',
                   (reaction_id, 'rhea:10000', compound_id, 0, 'chebi:1', 'Water', -1))
        db.execute('INSERT INTO reaction2protein VALUES (%s)' %
                   ','.join('?' for _ in range(6 if modern else 5)),
                   (reaction_id, 'rhea:10000', gene_id, 'uniprot:P00533', 'EGFR') +
                   ((1,) if modern else ()))
        if modern:
            db.execute('''CREATE TABLE reaction_ec_class(ramp_rxn_id TEXT,rxn_source_id TEXT,
                rxn_class_ec TEXT,ec_level INTEGER,rxn_class TEXT,rxn_class_hierarchy TEXT)''')
            db.execute('INSERT INTO reaction_ec_class VALUES (?,?,?,?,?,?)',
                       (reaction_id, 'rhea:10000', '1.1.1.1', 4, 'Oxidoreductases',
                        'Oxidoreductases | EC 1.1.1.1'))


def test_entity_status_report_shows_stored_zero_absent_row_and_source_names(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')]
    for path in (old, new):
        make_database(path, rows)
        with sqlite3.connect(path) as db:
            db.execute('''CREATE TABLE entity_status_info(status_category TEXT,
                entity_source_id TEXT,entity_source_name TEXT,entity_count INTEGER)''')
    with sqlite3.connect(old) as db:
        db.execute("INSERT INTO entity_status_info VALUES ('Metabolites','hmdb','HMDB',2)")
    with sqlite3.connect(new) as db:
        db.execute("INSERT INTO entity_status_info VALUES ('Metabolites','hmdb','HMDB',0)")
        db.execute("INSERT INTO entity_status_info VALUES ('Metabolites','novel','Novel Input',1)")
    report = compare([('Old', old), ('New', new)])
    assert report['databases'][1]['entity_status_info']['by_key']['Metabolites / hmdb']['entity_count'] == 0
    page = write(report, tmp_path / 'comparison.html')[0].read_text()
    assert 'Entity counts by category and source' in page
    assert 'Metabolites / <code>hmdb</code>' in page
    assert 'Novel Input' in page
    assert 'Row absent' in page
    assert 'Δ -2' in page


def test_chemistry_and_all_source_version_rows_are_reported(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    source_rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
                   ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'hmdb')]
    for path in (old, new):
        make_database(path, source_rows)
    with sqlite3.connect(old) as db:
        db.executescript('''
            CREATE TABLE chem_props(ramp_id TEXT,chem_data_source TEXT,chem_source_id TEXT,
                                    iso_smiles TEXT,common_name TEXT,mw REAL);
            CREATE TABLE version_info(ramp_db_version TEXT,db_mod_date TEXT,status TEXT,
                                      data_source_id TEXT,data_source_name TEXT,
                                      data_source_url TEXT,data_source_version TEXT);
        ''')
        db.execute("INSERT INTO chem_props VALUES ('C1','hmdb','hmdb:HMDB0000122',NULL,'Glucose',180.0)")
        db.executemany('INSERT INTO version_info VALUES (?,?,?,?,?,?,?)', [
            ('2.3', '2023-01-01', 'archive', 'hmdb', 'HMDB', 'https://hmdb.ca/', 'v4'),
            ('3.0', '2025-01-01', 'current', 'hmdb', 'HMDB', 'https://hmdb.ca/', 'v5'),
        ])
    with sqlite3.connect(new) as db:
        db.executescript('''
            CREATE TABLE chem_props(ramp_id TEXT,chem_data_source TEXT,chem_source_id TEXT,
                                    iso_smiles TEXT,common_name TEXT,mw REAL);
            CREATE TABLE version_info(ramp_db_version TEXT,db_mod_date TEXT,status TEXT,
                                      data_source_id TEXT,data_source_name TEXT,
                                      data_source_url TEXT,data_source_version TEXT,
                                      data_source_snapshot_ids TEXT);
        ''')
        db.executemany('INSERT INTO chem_props VALUES (?,?,?,?,?,?)', [
            ('C1', 'hmdb', 'hmdb:HMDB0000122', 'C(C1C(C(C(C(O1)O)O)O)O)O', 'Glucose', 180.0),
            ('C1', 'pubchem', 'pubchem:5793', None, 'D-Glucose', 180.0),
        ])
        db.executemany('INSERT INTO version_info VALUES (?,?,?,?,?,?,?,?)', [
            ('4.0', '2026-01-01', 'current', 'hmdb', 'HMDB', 'https://hmdb.ca/', 'v5', '["hmdb:5"]'),
            ('4.0', '2026-01-01', 'current', 'pubchem', 'PubChem', 'https://pubchem.ncbi.nlm.nih.gov/',
             'Derived dataset', '["pubchem:1"]'),
        ])
    report = compare([('Old', old), ('New', new)])
    counts = tables(report)['chemistry_counts']
    hmdb_rows = next(row for row in counts if row['metric'] == 'hmdb / rows')
    pubchem_rows = next(row for row in counts if row['metric'] == 'pubchem / rows')
    assert [cell['value'] for cell in hmdb_rows['cells']] == [1, 1]
    assert [cell['value'] for cell in pubchem_rows['cells']] == [0, 1]
    assert report['chemistry_examples']['hmdb']['chem_source_id'] == 'hmdb:HMDB0000122'
    assert [row['state'] for row in report['chemistry_examples']['hmdb']['rows']] == ['available', 'available']
    assert report['chemistry_examples']['pubchem']['rows'][0]['state'] == 'source_id_absent'
    assert len(report['databases'][0]['version_info']['rows']) == 2
    assert report['databases'][0]['version_info']['rows'][0]['status'] == 'archive'
    assert len(report['databases'][1]['version_info']['rows']) == 2
    page = write(report, tmp_path / 'chemistry.html')[0].read_text()
    assert 'Chemical properties by source' in page
    assert 'pubchem:5793' in page
    assert '1 / 1 (100.0%)' in page
    assert 'v4' in page and 'v5' in page
    assert 'hmdb:5' in page
    assert 'Column absent' in page
    assert '2 stored rows' in page


def test_reaction_tabs_compare_stable_source_rows_and_missing_legacy_columns(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    for path, reaction_id, compound_id, gene_id, modern in (
            (old, 'RAMP_R_000000001', 'C1', 'G1', False),
            (new, 'RAMP_R_000000999', 'C8', 'G8', True)):
        make_database(path, [
            ('hmdb:HMDB0000122', compound_id, 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', gene_id, 'uniprot', 'gene', 'uniprot')])
        add_reactions(path, reaction_id=reaction_id, compound_id=compound_id,
                      gene_id=gene_id, modern=modern)
    report = compare([('Old', old), ('New', new)])
    for table in ('reaction', 'reaction2met', 'reaction2protein'):
        rows = report['reaction_examples'][table]['rows']
        assert [row['state'] for row in rows] == ['available', 'available']
        assert [row['cells']['ramp_rxn_id'] for row in rows] == [
            'RAMP_R_000000001', 'RAMP_R_000000999']
        assert [cell['value'] for cell in tables(report)[table][0]['cells']] == [1, 1]
    assert report['reaction_examples']['reaction2met']['anchor'] == {
        'rxn_source_id': 'rhea:10000', 'met_source_id': 'chebi:1',
        'substrate_product': 0}
    assert report['reaction_examples']['reaction2protein']['rows'][0]['cells']['is_reviewed'] == {
        'state': 'column_absent'}
    assert report['reaction_examples']['reaction_ec_class']['rows'][0]['state'] == 'table_absent'
    assert report['reaction_examples']['reaction_ec_class']['rows'][1]['state'] == 'available'
    assert report['databases'][0]['reaction_tables']['reaction2met']['field_coverage']['fields']['is_cofactor'] == 0
    assert [cell['value'] for cell in tables(report)['reaction_ec_class'][0]['cells']] == [None, 1]
    page = write(report, tmp_path / 'reactions.html')[0].read_text()
    assert 'Reaction–metabolite participants' in page
    assert 'Reaction–protein participants' in page
    assert 'Reaction EC classes' in page
    assert 'rhea:10000' in page and 'chebi:1' in page and 'uniprot:P00533' in page
    assert 'Table absent' in page and 'Column absent' in page


def test_reaction_report_surfaces_links_to_missing_reaction_rows(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    for path in (old, new):
        make_database(path, [
            ('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'uniprot')])
        add_reactions(path, reaction_id='RAMP_R_1', compound_id='C1', gene_id='G1', modern=True)
    with sqlite3.connect(old) as db:
        db.execute('INSERT INTO reaction2met VALUES (?,?,?,?,?,?,?)',
                   ('RAMP_R_missing', 'rhea:99999', 'C1', 0, 'chebi:1', 'Water', 0))
    report = compare([('Old', old), ('New', new)])
    metric = next(row for row in tables(report)['reaction2met']
                  if row['metric'] == 'Missing reaction target / rows')
    assert [cell['value'] for cell in metric['cells']] == [1, 0]
    assert metric['cells'][1]['delta'] == -1
    assert 'Missing reaction target' in write(report, tmp_path / 'orphan.html')[0].read_text()


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
    assert page.count('<button type="button" role="tab"') == 17
    assert page.count('<section class="tab-panel"') == 17
    assert 'id="pathway_similarity" role="tabpanel"' in page
    assert 'Database version and source overlaps' in page
    assert 'id="db_version" role="tabpanel"' in page
    assert page.index('id="analyte" role="tabpanel"') < page.index('id="source" role="tabpanel"')
    assert page.index('id="source" role="tabpanel"') < page.index('id="analytesynonym" role="tabpanel"')
    assert page.index('id="analytesynonym" role="tabpanel"') < page.index('id="analytehasontology" role="tabpanel"')
    assert page.index('id="analytehasontology" role="tabpanel"') < page.index('id="analytehaspathway" role="tabpanel"')
    assert page.index('id="analytehaspathway" role="tabpanel"') < page.index('id="catalyzed" role="tabpanel"')
    assert page.index('id="catalyzed" role="tabpanel"') < page.index('id="metabolite_class" role="tabpanel"')
    assert page.index('id="metabolite_class" role="tabpanel"') < page.index('id="about" role="tabpanel"')
    assert page.index('id="metabolite_class" role="tabpanel"') < page.index('id="reaction" role="tabpanel"')
    assert page.index('id="reaction" role="tabpanel"') < page.index('id="reaction2met" role="tabpanel"')
    assert page.index('id="reaction2met" role="tabpanel"') < page.index('id="reaction2protein" role="tabpanel"')
    assert page.index('id="reaction2protein" role="tabpanel"') < page.index('id="reaction_ec_class" role="tabpanel"')
    assert page.index('id="reaction_ec_class" role="tabpanel"') < page.index('id="about" role="tabpanel"')
    assert page.index('Example source rows') < page.index('id="analytesynonym" role="tabpanel"')
    assert page.index('Generated ') > page.index('id="about" role="tabpanel"')
    assert "history.pushState(null, '', '#' + id)" in page
    assert '.tab-panel[hidden]{display:block}' in page
    assert set(json.loads(json_path.read_text())['comparison_tables']) == {
        'ramp_ids', 'metabolite_sources', 'gene_sources',
        'compound_synonyms', 'gene_synonyms', 'ontology_associations',
        'pathway_associations', 'catalyzed_counts', 'class_counts', 'chemistry_counts', 'lookup_quality',
        'reaction', 'reaction2met', 'reaction2protein', 'reaction_ec_class'}


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


def test_catalyzed_and_class_tabs_match_stable_source_ids_and_show_coverage(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    for path, compound, gene, protein_type, class_name in (
            (old, 'C1', 'G1', 'Enzyme', 'Carbohydrates'),
            (new, 'C8', 'G8', 'Unknown', 'Organic compounds')):
        make_database(path, [
            ('hmdb:HMDB0000122', compound, 'hmdb', 'compound', 'hmdb'),
            ('lipidmaps:L1', compound, 'lipidmaps', 'compound', 'lipidmaps'),
            ('uniprot:P00533', gene, 'uniprot', 'gene', 'hmdb')])
        add_catalyzed_classes(path, compound_id=compound, gene_id=gene,
                              protein_type=protein_type, class_name=class_name)
    report = compare([('Old', old), ('New', new)])
    catalyzed = report['catalyzed_example']
    assert (catalyzed['compound_source_id'], catalyzed['protein_source_id']) == (
        'hmdb:HMDB0000122', 'uniprot:P00533')
    assert [row['cells']['rampCompoundId'] for row in catalyzed['rows']] == ['C1', 'C8']
    assert [row['cells']['proteinType'] for row in catalyzed['rows']] == ['Enzyme', 'Unknown']
    assert report['databases'][1]['catalyzed_field_coverage']['fields']['proteinType'] == 0
    assert [row['cells']['class_name'] for row in report['class_examples']['hmdb']['rows']] == [
        'Carbohydrates', 'Organic compounds']
    assert [row['cells']['ramp_id'] for row in report['class_examples']['lipidmaps']['rows']] == [
        'C1', 'C8']
    assert [cell['value'] for cell in next(row for row in tables(report)['class_counts']
           if row['metric'] == 'hmdb / ClassyFire_class / rows')['cells']] == [1, 1]
    page = write(report, tmp_path / 'classes.html')[0].read_text()
    assert 'Catalyzed metabolite–gene associations' in page
    assert 'Metabolite classes by source and level' in page
    assert 'Example catalyzed association' in page
    assert 'Example metabolite class rows' in page
    assert 'uniprot:P00533' in page and 'lipidmaps:L1' in page
    assert '0 / 1 (0.0%)' in page


def test_new_association_tables_distinguish_absent_from_zero(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    rows = [('hmdb:HMDB0000122', 'C1', 'hmdb', 'compound', 'hmdb'),
            ('uniprot:P00533', 'G1', 'uniprot', 'gene', 'hmdb')]
    make_database(old, rows)
    make_database(new, rows)
    add_catalyzed_classes(new, compound_id='C1', gene_id='G1',
                          protein_type='Enzyme', class_name='Carbohydrates')
    with sqlite3.connect(new) as db:
        db.execute("INSERT INTO source(sourceId,rampId,IDtype,geneOrCompound,dataSource) "
                   "VALUES ('lipidmaps:L1','C1','lipidmaps','compound','lipidmaps')")
    report = compare([('Old', old), ('New', new)])
    catalyzed = next(row for row in tables(report)['catalyzed_counts']
                     if row['metric'] == 'All associations / rows')
    assert [cell['value'] for cell in catalyzed['cells']] == [None, 1]
    assert report['catalyzed_example']['rows'][0]['state'] == 'table_absent'
    assert report['class_examples']['hmdb']['rows'][0]['state'] == 'table_absent'
    page = write(report, tmp_path / 'absent.html')[0].read_text()
    assert 'Table absent' in page


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
