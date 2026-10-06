import hashlib
import copy
import json
import os
import sqlite3

import pytest

from src.use_cases.ramp.compare_sqlite import compare, get_metric, main, prepare_tables, render_html, write_report, metabolite_group_metrics, GROUP_METRICS


def database(path, *, new=False):
    with sqlite3.connect(path) as db:
        db.executescript("""
            CREATE TABLE analyte(rampId TEXT PRIMARY KEY,type TEXT);
            CREATE TABLE source(sourceId TEXT,rampId TEXT,IDtype TEXT,geneOrCompound TEXT,dataSource TEXT);
            CREATE TABLE chem_props(ramp_id TEXT,chem_data_source TEXT,chem_source_id TEXT,mw REAL,inchi_key TEXT);
            CREATE TABLE db_version(ramp_version TEXT);
            CREATE TABLE version_info(data_source_id TEXT,data_source_version TEXT,status TEXT);
            CREATE TABLE entity_status_info(n INTEGER);
            CREATE TABLE ramp_export_metadata(key TEXT,value TEXT);
        """)
        rid = "fresh" if new else "old"
        db.execute("INSERT INTO analyte VALUES (?, 'compound')", (rid,))
        db.executemany("INSERT INTO source VALUES (?,?,?,?,?)", [
            ("hmdb:HMDB1" if new else "HMDB:HMDB1", rid, "hmdb", "compound", "hmdb"),
            ("chebi:2" if new else "chebi:1", rid, "chebi", "compound", "hmdb"),
            ("hmdb:HMDB1" if new else "HMDB:HMDB1", rid, "hmdb", "compound", "chebi"),
        ])
        db.executemany("INSERT INTO chem_props VALUES (?, 'pubchem', ?, ?, ?)", [
            (rid, "pubchem:1", None, ""), (rid, "pubchem:1", 1.0, "KEY"), (rid, "pubchem:2", None, None),
        ])
        db.execute("INSERT INTO version_info VALUES ('hmdb','5.0','current')")
        db.execute("INSERT INTO db_version VALUES (?)", ("unreleased" if new else "3.0.7",))
        if new:
            db.execute("INSERT INTO ramp_export_metadata VALUES ('manifest', ?)", (json.dumps({
                "status": "post_processing_pending", "pending_tables": ["entity_status_info"],
                "omitted_tables": ["reaction_protein2met"], "gene_policy": "Separate gene/protein identities",
            }),))
        else:
            db.execute("INSERT INTO entity_status_info VALUES (4)")
            db.execute("CREATE TABLE reaction_protein2met(n INTEGER)")
            db.execute("INSERT INTO reaction_protein2met VALUES (1)")


def test_comparison_uses_source_ids_and_distinct_chemistry_denominators(tmp_path):
    old, new = tmp_path / "old.sqlite", tmp_path / "new.sqlite"
    database(old)
    database(new, new=True)
    hashes = [hashlib.sha256(p.read_bytes()).hexdigest() for p in (old, new)]
    report = prepare_tables(compare([("Old", old), ("New", new)], progress=lambda _: None))
    profile = report["databases"][1]
    assert profile["source_id_changes"]["hmdb"] == {"retained": 1, "added": 0, "removed": 0}
    assert profile["source_id_changes"]["chebi"] == {"retained": 0, "added": 1, "removed": 1}
    assert profile["sections"]["Source attribution"]["hmdb / compound / source IDs"]["value"] == 2
    assert profile["sections"]["Chemistry completeness"]["pubchem / mw"] == {"value": 1, "denominator": 2, "state": "available"}
    assert profile["sections"]["Chemistry completeness"]["pubchem / iso_smiles"]["state"] == "unavailable"
    assert "pubchem" in " ".join(profile["notes"])
    assert get_metric(profile, "Table rows", "reaction_protein2met")["state"] == "omitted"
    pending = next(r for r in report["comparison_tables"]["Table rows"] if r["metric"] == "entity_status_info")["cells"][1]
    assert pending["text"] == "Pending (0 rows)"
    assert pending["delta"] is None
    assert hashes == [hashlib.sha256(p.read_bytes()).hexdigest() for p in (old, new)]
    html_path, json_path = write_report(report, tmp_path / "report.html")
    assert json.loads(json_path.read_text()) == report
    assert "1 / 2 (50.0%)" in html_path.read_text()


def test_group_size_bins_cover_boundaries_without_counting_duplicate_evidence():
    with sqlite3.connect(':memory:') as db:
        db.executescript('CREATE TABLE analyte(rampId TEXT,type TEXT); CREATE TABLE source(sourceId TEXT,rampId TEXT);')
        sizes = [0,1,2,5,6,10,11,20,21,50,51,100,101,500,501,1000,1001]
        for i, size in enumerate(sizes):
            rid = str(i)
            db.execute("INSERT INTO analyte VALUES (?, 'compound')", (rid,))
            db.executemany('INSERT INTO source VALUES (?,?)', [(f'x:{j}',rid) for j in range(size)] * 2)
        db.execute("INSERT INTO analyte VALUES ('gene','gene')")
        db.executemany("INSERT INTO source VALUES (?,'gene')", [(str(i),) for i in range(2000)])
        metrics = metabolite_group_metrics(db)
        assert [metrics[k]['value'] for k in GROUP_METRICS] == [1,1,2,2,2,2,2,2,2,1,1001]
        assert sum(metrics[k]['value'] for k in GROUP_METRICS[:-1]) == len(sizes)


def test_multiple_ramp_mappings_are_counted_within_entity_type(tmp_path):
    path = tmp_path / 'mappings.sqlite'
    database(path)
    with sqlite3.connect(path) as db:
        db.executemany('INSERT INTO source VALUES (?,?,?,?,?)', [
            ('shared-across-types', 'C1', 'test', 'compound', 'hmdb'),
            ('shared-across-types', 'G1', 'test', 'gene', 'hmdb'),
            ('shared-across-types', 'G1', 'test', 'gene', 'uniprot'),
            ('multi-metabolite', 'C1', 'test', 'compound', 'hmdb'),
            ('multi-metabolite', 'C2', 'test', 'compound', 'hmdb'),
            ('multi-metabolite', 'C2', 'test', 'compound', 'chebi'),
            ('multi-gene', 'G1', 'test', 'gene', 'reactome'),
            ('multi-gene', 'G2', 'test', 'gene', 'rhea'),
            ('multi-gene', 'G2', 'test', 'gene', 'uniprot'),
        ])
    report = prepare_tables(compare([('Build', path)], progress=lambda _: None))
    core = report['databases'][0]['sections']['Core coverage']
    assert core['Source IDs linked to multiple metabolite RaMP IDs']['value'] == 1
    assert core['Source IDs linked to multiple gene/protein RaMP IDs']['value'] == 1
    assert core['Source IDs linked to both metabolites and genes/proteins']['value'] == 1
    assert 'Source IDs linked to multiple RaMP IDs' not in core
    rows = [row['metric'] for row in report['comparison_tables']['Core coverage']]
    index = rows.index('Unique source identifiers')
    assert rows[index + 1:index + 4] == [
        'Source IDs linked to multiple metabolite RaMP IDs',
        'Source IDs linked to multiple gene/protein RaMP IDs',
        'Source IDs linked to both metabolites and genes/proteins']
    page = render_html(report)
    assert 'Genes and proteins are counted together' in page
    assert 'These counts can overlap' in page
    assert 'Source IDs linked to multiple RaMP IDs' not in page


def test_multiple_ramp_mappings_require_entity_type_column(tmp_path):
    path = tmp_path / 'missing_type.sqlite'
    database(path)
    with sqlite3.connect(path) as db:
        db.execute('ALTER TABLE source RENAME TO source_with_type')
        db.execute('CREATE TABLE source AS SELECT sourceId,rampId,IDtype,dataSource FROM source_with_type')
        db.execute('DROP TABLE source_with_type')
    core = compare([('Old', path)], progress=lambda _: None)['databases'][0]['sections']['Core coverage']
    assert all(core[name]['state'] == 'unavailable' for name in (
        'Source IDs linked to multiple metabolite RaMP IDs',
        'Source IDs linked to multiple gene/protein RaMP IDs',
        'Source IDs linked to both metabolites and genes/proteins'))


def test_metabolites_lead_with_total_and_distinct_provider_coverage(tmp_path):
    old, new = tmp_path / 'old.sqlite', tmp_path / 'new.sqlite'
    database(old)
    database(new, new=True)
    with sqlite3.connect(new) as db:
        db.execute("INSERT INTO analyte VALUES ('C2','compound')")
        db.execute("INSERT INTO analyte VALUES ('G1','gene')")
        db.executemany('INSERT INTO source VALUES (?,?,?,?,?)', [
            ('hmdb:other', 'C2', 'hmdb', 'compound', 'hmdb'),
            ('hmdb:other', 'C2', 'hmdb', 'compound', 'hmdb'),
            ('pubchem:42', 'C2', 'pubchem', 'compound', 'pubchem'),
            ('uniprot:P1', 'G1', 'uniprot', 'gene', 'uniprot'),
        ])
    report = prepare_tables(compare([('Old', old), ('New', new)], progress=lambda _: None))
    assert next(iter(report['comparison_tables'])) == 'Metabolites'
    rows = {row['metric']: row['cells'] for row in report['comparison_tables']['Metabolites']}
    assert next(iter(rows)) == 'Total metabolites'
    assert [c['value'] for c in rows['Total metabolites']] == [1, 2]
    assert rows['Total metabolites'][1]['delta'] == 1
    assert [c['value'] for c in rows['With data from HMDB']] == [1, 2]
    assert rows['With data from HMDB'][1]['delta'] == 1
    assert [c['value'] for c in rows['With data from ChEBI']] == [1, 1]
    assert [c['value'] for c in rows['With data from PubChem']] == [0, 1]
    assert not any('UniProt' in key for key in rows)
    assert report['databases'][1]['sections']['Core coverage']['Metabolites']['value'] == 2
    assert report['databases'][1]['sections']['Source attribution']['hmdb / compound / rows']['value'] == 4
    page = render_html(report)
    assert page.index('<h2 id=\'s0\'>Metabolites</h2>') < page.index('<h2 id=\'s1\'>Core coverage</h2>')
    assert '<strong>Total metabolites</strong>' in page
    assert 'source counts do not add up to the total' in page
    assert 'Chemistry recorded only in chem_props' in page


def test_metabolite_summary_uses_core_exclusion_and_reports_missing_source_schema(tmp_path):
    old, hidden, missing = (tmp_path / f'{name}.sqlite' for name in ('old', 'hidden', 'missing'))
    for path in (old, hidden, missing):
        database(path)
    with sqlite3.connect(missing) as db:
        db.execute('ALTER TABLE source RENAME TO source_with_type')
        db.execute('CREATE TABLE source AS SELECT sourceId,rampId,IDtype,dataSource FROM source_with_type')
        db.execute('DROP TABLE source_with_type')
    report = compare([('Old', old), ('Hidden', hidden), ('Missing', missing)], progress=lambda _: None)
    report['core_coverage_excluded_labels'] = ['Hidden']
    prepare_tables(report)
    assert [d['label'] for d in report['table_databases']['Metabolites']] == ['Old', 'Missing']
    assert [d['label'] for d in report['table_databases']['Core coverage']] == ['Old', 'Missing']
    assert [d['label'] for d in report['table_databases']['Source attribution']] == ['Old', 'Hidden', 'Missing']
    rows = {row['metric']: row['cells'] for row in report['comparison_tables']['Metabolites']}
    assert [c['value'] for c in rows['Total metabolites']] == [1, 1]
    assert rows['With data from HMDB'][1]['state'] == 'unavailable'


def test_metabolite_provider_case_variants_use_distinct_union(tmp_path):
    path = tmp_path / 'providers.sqlite'
    database(path)
    with sqlite3.connect(path) as db:
        db.executemany('INSERT INTO source VALUES (?,?,?,?,?)', [
            ('chebi:2', 'C1', 'chebi', 'compound', 'ChEBI'),
            ('chebi:3', 'C2', 'chebi', 'compound', 'chebi'),
            ('chebi:4', 'C2', 'chebi', 'compound', 'ChEBI'),
        ])
    profile = compare([('Build', path)], progress=lambda _: None)['databases'][0]
    assert profile['sections']['Metabolites']['With data from ChEBI']['value'] == 3
    # The existing raw-attribution rows retain their original provider spelling.
    assert profile['sections']['Source attribution']['ChEBI / compound / analytes']['value'] == 2
    assert profile['sections']['Source attribution']['chebi / compound / analytes']['value'] == 2


def test_single_version_matrix_and_simplified_html(tmp_path):
    old, new = tmp_path/'old.sqlite', tmp_path/'new.sqlite'
    database(old); database(new)
    with sqlite3.connect(new) as db:
        db.execute("INSERT INTO version_info VALUES ('ChEBI','255','current')")
        db.execute("INSERT INTO version_info VALUES ('chebi','254','current')")
    report = prepare_tables(compare([('Old',old),('New',new)],progress=lambda _:None))
    assert report['source_version_table'] == [
        {'source':'chebi','versions':['Not recorded','254; 255']},
        {'source':'hmdb','versions':['5.0','5.0']}]
    page = render_html(report)
    version_section = page.split('<h2>Recorded source versions</h2>')[1].split('<h2>Methodology</h2>')[0]
    assert version_section.count('<table>') == 1
    assert '<th>Old</th><th>New</th>' in version_section
    assert 'Pinned snapshot IDs' not in version_section
    assert 'Retained, added, and removed' not in page
    assert [r['metric'] for r in report['comparison_tables']['Metabolite group sizes']] == GROUP_METRICS


def test_all_tables_use_previous_visible_database(tmp_path):
    path = tmp_path / 'input.sqlite'
    database(path)
    report = compare([('First',path)],progress=lambda _:None)
    template = report['databases'][0]
    report['databases'] = []
    for label,value in [('First',10),('Second',20),('Hidden',90),('New',35)]:
        profile = copy.deepcopy(template)
        profile['label'] = label
        profile['sections']['Core coverage']['Metabolites']['value'] = value
        profile['sections']['Table rows']['source']['value'] = value
        report['databases'].append(profile)
    report['core_coverage_excluded_labels'] = ['Hidden']
    prepare_tables(report)
    assert report['table_databases']['Core coverage'] == [
        {'label':'First','compared_with':None}, {'label':'Second','compared_with':'First'},
        {'label':'New','compared_with':'Second'}]
    row = next(r for r in report['comparison_tables']['Core coverage'] if r['metric']=='Metabolites')
    assert [c['delta'] for c in row['cells']] == [None,10,15]
    assert row['cells'][-1]['percent_change'] == 75
    other = next(r for r in report['comparison_tables']['Table rows'] if r['metric']=='source')
    assert [c['delta'] for c in other['cells']] == [None,10,70,-55]
    assert other['cells'][-1]['percent_change'] == pytest.approx(-100*55/90)
    assert report['table_databases']['Table rows'][-1]['compared_with'] == 'Hidden'
    page = render_html(report)
    metabolites = page.split("<h2 id='s0'>")[1].split("<h2 id='s1'>")[0]
    core = page.split("<h2 id='s1'>")[1].split("<h2 id='s2'>")[0]
    assert '<th>Hidden</th>' not in metabolites
    assert '<th>Hidden</th>' not in core
    assert '<th>Hidden</th>' in page.split("<h2 id='s2'>")[1]
    report['databases'][1]['sections']['Core coverage']['Metabolites'] = {'value':None,'state':'unavailable','denominator':None}
    prepare_tables(report)
    row = next(r for r in report['comparison_tables']['Core coverage'] if r['metric']=='Metabolites')
    assert row['cells'][-1]['delta'] is None


def test_unavailable_is_not_zero_and_orphans_include_nulls(tmp_path):
    old, new = tmp_path / "old.sqlite", tmp_path / "new.sqlite"
    database(old)
    database(new)
    with sqlite3.connect(new) as db:
        db.execute("DELETE FROM source")
        db.executemany("INSERT INTO source(sourceId,rampId) VALUES (?,?)", [("x:1", None), ("x:2", "missing")])
        db.execute("DROP TABLE chem_props")
    report = prepare_tables(compare([("Old", old), ("New", new)], progress=lambda _: None))
    p = report["databases"][1]
    assert p["sections"]["Relationship integrity"]["source.rampId → analyte: orphan rows"]["value"] == 2
    assert get_metric(p, "Chemistry completeness", "pubchem / mw")["state"] == "unavailable"
    assert get_metric(p, "Source identifiers by namespace", "hmdb")["value"] == 0
    missing = next(r for r in report["comparison_tables"]["Core coverage"] if r["metric"] == "Metabolites with chemistry")["cells"][1]
    assert missing["delta"] is None


def test_html_escapes_labels_and_metadata(tmp_path):
    path = tmp_path / "a.sqlite"
    database(path)
    report = prepare_tables(compare([("<script>alert(1)</script>", path), ("Other", path)], progress=lambda _: None))
    rendered = render_html(report)
    assert "<script>" not in rendered
    assert "&lt;script&gt;" in rendered


def test_absent_provider_and_missing_column_is_unavailable_not_loss(tmp_path):
    old, new = tmp_path / "old.sqlite", tmp_path / "new.sqlite"
    database(old)
    database(new)
    with sqlite3.connect(new) as db:
        db.execute("DROP TABLE chem_props")
        db.execute("CREATE TABLE chem_props(ramp_id TEXT,chem_data_source TEXT,chem_source_id TEXT,inchi_key TEXT)")
    report = prepare_tables(compare([("Old", old), ("New", new)], progress=lambda _: None))
    row = next(r for r in report["comparison_tables"]["Chemistry completeness"] if r["metric"] == "pubchem / mw")
    assert row["cells"][1]["state"] == "unavailable"
    assert row["cells"][1]["delta"] is None
    assert get_metric(report["databases"][1], "Chemistry completeness", "pubchem / inchi_key")["value"] == 0


def test_report_refuses_to_replace_input_or_hard_link(tmp_path):
    path = tmp_path / "input.sqlite"
    database(path)
    report = prepare_tables(compare([("A", path), ("B", path)], progress=lambda _: None))
    output = tmp_path / "report.html"
    os.link(path, output)
    with pytest.raises(ValueError, match="overwrite"):
        write_report(report, output)


def test_missing_input_does_not_create_database(tmp_path):
    missing = tmp_path / "missing.sqlite"
    with pytest.raises(FileNotFoundError):
        compare([("missing", missing)], progress=lambda _: None)
    assert not missing.exists()


@pytest.mark.parametrize("values", [["same=a", "same=b"], ["only=a"], ["missing_equals", "b=b"]])
def test_cli_rejects_invalid_inputs(values, tmp_path):
    argv = ["--output", str(tmp_path / "report.html")]
    for value in values:
        argv.extend(["--database", value])
    with pytest.raises(SystemExit):
        main(argv)
