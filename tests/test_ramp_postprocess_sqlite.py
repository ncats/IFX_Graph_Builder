"""Derived RaMP counts can be refreshed without returning to the graph."""
import json
import sqlite3

import pytest

from src.use_cases.ramp import postprocess_sqlite
from src.use_cases.ramp.postprocess_sqlite import postprocess_entity_status
from src.use_cases.ramp import sqlite_human_reactions
from tests.test_ramp_sqlite_export import FixtureReader, export_sqlite


def test_postprocess_counts_sources_and_is_repeatable(tmp_path):
    path = tmp_path / 'ramp.sqlite'
    export_sqlite(FixtureReader(), path, source_only=True, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        compound = db.execute("SELECT rampId FROM analyte WHERE type='compound' ORDER BY rampId LIMIT 1").fetchone()[0]
        # A new input must be discovered without changing a hard-coded status source list.
        db.execute('''INSERT INTO source
            (sourceId,rampId,IDtype,geneOrCompound,commonName,priorityHMDBStatus,dataSource,pathwayCount)
            VALUES (?,?,?,?,?,?,?,?)''',
            ('novel:X1', compound, 'novel', 'compound', None, None, 'novel', -1))
        db.execute('''INSERT INTO version_info
            (ramp_db_version,db_mod_date,status,data_source_id,data_source_name,
             data_source_url,data_source_version)
            VALUES ('test','2026-10-08','current','novel','Novel Input','','test')''')
        db.execute("INSERT INTO pathway VALUES ('RAMP_P_EXTRA','WPX','wiki',NULL,'Extra path')")
        db.execute('INSERT INTO analytehaspathway VALUES (?,?,?)',
                   (compound, 'RAMP_P_EXTRA', 'wiki'))
        db.commit()

    first = postprocess_entity_status(path)
    with sqlite3.connect(path) as db:
        rows = db.execute('SELECT * FROM entity_status_info ORDER BY 1,2').fetchall()
        manifest = json.loads(db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
        assert ('Metabolites', 'novel', 'Novel Input', 1) in rows
        assert {row[0] for row in rows} == {
            'Metabolites', 'Genes', 'Pathways', 'Metabolite-Pathway Associations',
            'Gene-Pathway Associations', 'Metabolite-Reaction Associations',
            'Gene-Reaction Associations', 'Metabolite-Gene Associations',
            'Chemical Property Records',
        }
        assert db.execute('SELECT DISTINCT pathwayCount FROM source WHERE rampId=?',
                          (compound,)).fetchall() == [(1,)]
        assert db.execute('SELECT metCount FROM ontology').fetchall() == [(1,)]
        assert manifest['row_counts']['entity_status_info'] == len(rows)
        assert manifest['row_counts']['db_version'] == 1
        assert manifest['row_counts']['pathway_duplicates'] == 0
        assert manifest['row_counts']['pathway_similarity'] == 0
        assert manifest['post_processing']['pathway_tables']['status'] == 'schema_only'
        for table, (_, expected) in postprocess_sqlite.EMPTY_PATHWAY_TABLES.items():
            assert tuple(row[1] for row in db.execute(f'PRAGMA table_info({table})')) == expected
            assert db.execute(f'SELECT count(*) FROM {table}').fetchone()[0] == 0
        version = db.execute('SELECT ramp_version,load_timestamp,met_intersects_json FROM db_version').fetchone()
        assert version[:2] == (manifest['release_version'], manifest['created_at'])
        assert sum(item['size'] for item in json.loads(version[2])) == db.execute(
            "SELECT count(DISTINCT rampId) FROM source WHERE geneOrCompound='compound'").fetchone()[0]
        assert 'source.pathwayCount' not in manifest['pending_fields']
        assert 'ontology.metCount' not in manifest['pending_fields']
        assert 'entity_status_info' in manifest['included_tables']
    assert first['rows'] == len(rows)
    postprocess_entity_status(path)
    with sqlite3.connect(path) as db:
        assert db.execute('SELECT * FROM entity_status_info ORDER BY 1,2').fetchall() == rows
        assert db.execute('SELECT ramp_version,load_timestamp FROM db_version').fetchone() == version[:2]


def test_intersections_partition_sources_and_exclude_smpdb_pathways():
    with sqlite3.connect(':memory:') as db:
        db.executescript('''
            CREATE TABLE source(rampId TEXT,geneOrCompound TEXT,dataSource TEXT);
            CREATE TABLE pathway(pathwayRampId TEXT,pathwayCategory TEXT);
            CREATE TABLE analytehaspathway(rampId TEXT,pathwayRampId TEXT);
            CREATE TABLE version_info(data_source_id TEXT,data_source_name TEXT,status TEXT);
        ''')
        db.executemany('INSERT INTO source VALUES (?,?,?)', [
            ('C1','compound','hmdb'), ('C1','compound','hmdb_kegg'),
            ('C1','compound','wikipathways_kegg'), ('C1','compound','wiki'),
            ('C1','compound','chebi'), ('C2','compound','hmdb'),
            ('G1','gene','uniprot'),
        ])
        db.executemany('INSERT INTO pathway VALUES (?,?)',
                       [('P1',None), ('P2','smpdb3'), ('P3','hmdb:smpdb2')])
        db.executemany('INSERT INTO analytehaspathway VALUES (?,?)',
                       [('C1','P1'), ('C2','P2'), ('G1','P3')])
        result, summary = postprocess_sqlite.db_version_intersections(db)
    global_mets = json.loads(result['met_intersects_json'])
    mapped_mets = json.loads(result['met_intersects_json_pw_mapped'])
    assert [(row['sets'],row['size']) for row in global_mets] == [
        (['ChEBI','HMDB','KEGG','WikiPathways'],1), (['HMDB'],1)]
    assert [(row['sets'],row['size']) for row in mapped_mets] == [
        (['ChEBI','HMDB','KEGG','WikiPathways'],1)]
    assert json.loads(result['gene_intersects_json_pw_mapped']) == []
    assert summary['met_intersects_json']['analytes'] == 2
    assert summary['met_intersects_json_pw_mapped']['analytes'] == 1


def test_intersections_reject_colliding_source_display_names():
    with sqlite3.connect(':memory:') as db:
        db.executescript('''
            CREATE TABLE source(rampId TEXT,geneOrCompound TEXT,dataSource TEXT);
            CREATE TABLE pathway(pathwayRampId TEXT,pathwayCategory TEXT);
            CREATE TABLE analytehaspathway(rampId TEXT,pathwayRampId TEXT);
            CREATE TABLE version_info(data_source_id TEXT,data_source_name TEXT,status TEXT);
            INSERT INTO source VALUES ('C1','compound','novel_a');
            INSERT INTO source VALUES ('C2','compound','novel_b');
            INSERT INTO version_info VALUES ('novel_a','Shared','current');
            INSERT INTO version_info VALUES ('novel_b','Shared','current');
        ''')
        with pytest.raises(ValueError, match='ambiguous'):
            postprocess_sqlite.db_version_intersections(db)


def test_postprocess_rolls_back_on_invalid_source(tmp_path):
    path = tmp_path / 'ramp.sqlite'
    export_sqlite(FixtureReader(), path, source_only=True, progress=lambda _: None)
    postprocess_entity_status(path)
    with sqlite3.connect(path) as db:
        prior_rows = db.execute('SELECT * FROM entity_status_info ORDER BY 1,2').fetchall()
        prior_manifest = db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0]
        db.execute("UPDATE source SET geneOrCompound='ambiguous' WHERE rowid=(SELECT min(rowid) FROM source)")
        db.commit()
    with pytest.raises(ValueError, match='Unsupported source analyte type'):
        postprocess_entity_status(path)
    with sqlite3.connect(path) as db:
        assert db.execute('SELECT * FROM entity_status_info ORDER BY 1,2').fetchall() == prior_rows
        assert db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0] == prior_manifest


def test_human_reaction_flags_use_mixed_chebi_closure_and_all_raw_edges(monkeypatch):
    records = {
        'BiologicalRole': [{'id': root} for root in sqlite_human_reactions.HUMAN_CHEBI_ROOTS],
        'IsAEdge': [{'start_id': 'CHEBI:child-role', 'end_id': 'CHEBI:77746'},
                    {'start_id': 'CHEBI:child-metabolite', 'end_id': 'CHEBI:role-bearing'}],
        'HasBiologicalRoleEdge': [
            {'start_id': 'CHEBI:role-bearing', 'end_id': 'CHEBI:child-role'}],
        'RheaReaction': [{'id': 'RHEA:1'}, {'id': 'RHEA:2'}],
        'RheaMetaboliteReactionEdge': [
            {'end_id': 'RHEA:1', 'source_id': 'CHEBI:child-metabolite'},
            {'end_id': 'RHEA:2', 'source_id': 'CHEBI:outside'}],
        'RheaProteinReactionEdge': [
            {'end_id': 'RHEA:1', 'source_id': 'UniProtKB:HUMAN'},
            {'end_id': 'RHEA:2', 'source_id': 'UniProtKB:NONHUMAN'}],
    }
    monkeypatch.setattr(sqlite_human_reactions, '_records',
                        lambda graph, collection, fields: iter(records[collection]))
    flags, summary = sqlite_human_reactions.human_reaction_flags(None, {'UniProtKB:HUMAN'})
    assert flags == {'RHEA:1': (1, 1), 'RHEA:2': (0, 0)}
    assert summary['human_chebi_ids'] == len(sqlite_human_reactions.HUMAN_CHEBI_ROOTS) + 3


def test_postprocess_writes_supplied_human_flags(tmp_path):
    path = tmp_path / 'ramp.sqlite'
    export_sqlite(FixtureReader(), path, source_only=True, progress=lambda _: None)
    postprocess_entity_status(path, human_flags={'RHEA:1': (1, 0)},
                              human_manifest={'policy': 'fixture'})
    with sqlite3.connect(path) as db:
        assert db.execute('SELECT has_human_prot,only_human_mets FROM reaction').fetchone() == (1, 0)
        manifest = json.loads(db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
        assert 'reaction.has_human_prot' not in manifest['pending_fields']
        assert manifest['post_processing']['human_reaction_flags'] == {'policy': 'fixture'}


def test_postprocess_cli_refreshes_existing_sqlite_without_graph_on_rerun(
        tmp_path, monkeypatch, capsys):
    path = tmp_path / 'ramp.sqlite'
    export_sqlite(FixtureReader(), path, source_only=True, progress=lambda _: None)
    calls = []

    def graph_flags(manifest, **kwargs):
        calls.append(manifest['stage_id'])
        return {'RHEA:1': (1, 1)}, {'policy': 'fixture'}

    monkeypatch.setattr(postprocess_sqlite, 'graph_human_flags', graph_flags)
    assert postprocess_sqlite.main(['--sqlite', str(path), '--with-human-flags']) == 0
    assert calls == ['HarmonizationStage:test']
    assert 'Updated reaction human-scope flags: 1 reactions' in capsys.readouterr().out
    assert postprocess_sqlite.main(['--sqlite', str(path)]) == 0
    assert calls == ['HarmonizationStage:test']
    with sqlite3.connect(path) as db:
        assert db.execute('SELECT has_human_prot,only_human_mets FROM reaction').fetchone() == (1, 1)
