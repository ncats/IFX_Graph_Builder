"""The compressed pathway rows must match the legacy R decoding contract."""

import json
import sqlite3
import zlib

import pytest

from src.use_cases.ramp.pathway_similarity import sparse_rows
from src.use_cases.ramp.postprocess_sqlite import postprocess_pathway_similarity
from src.use_cases.ramp.compare_source import pathway_overlap_profile


def make_database(path):
    with sqlite3.connect(path) as db:
        db.executescript('''
            CREATE TABLE pathway(pathwayRampId TEXT PRIMARY KEY,type TEXT,
                                 sourceId TEXT,pathwayName TEXT);
            CREATE TABLE analytehaspathway(rampId TEXT,pathwayRampId TEXT);
            CREATE TABLE ramp_export_metadata(key TEXT PRIMARY KEY,value TEXT);
        ''')
        db.execute('INSERT INTO ramp_export_metadata VALUES (?,?)', (
            'manifest', json.dumps({'scope': 'lookup_table_diagnostic',
                                    'included_tables': [], 'pending_tables': [
                                        'pathway_similarity', 'pathway_duplicates']})))
        memberships = {
            'P1': ({f'RAMP_C_{i}' for i in range(1, 6)} |
                   {f'RAMP_G_{i}' for i in range(1, 6)}),
            'P2': ({f'RAMP_C_{i}' for i in range(1, 6)} |
                   {f'RAMP_G_{i}' for i in range(1, 6)}),
            'P3': ({f'RAMP_C_{i}' for i in range(1, 5)} |
                   {f'RAMP_G_{i}' for i in range(1, 7)}),
            'P4': ({f'RAMP_C_{i}' for i in range(1, 6)} |
                   {f'RAMP_G_{i}' for i in range(6, 11)}),
            'P5': ({f'RAMP_C_{i}' for i in range(1, 6)} |
                   {f'RAMP_G_{i}' for i in range(1, 6)}),
            'P6': {f'RAMP_C_{i}' for i in range(20, 25)},
            'P7': {f'RAMP_G_{i}' for i in range(20, 25)},
        }
        source_ids = {'P1': 'R-HSA-110329', 'P2': 'R-HSA-73928',
                      'P3': 'R-HSA-109581'}
        db.executemany('INSERT INTO pathway VALUES (?,?,?,?)',
                       [(p, 'hmdb' if p == 'P5' else 'wiki',
                         source_ids.get(p, p), f'Example {p}') for p in memberships])
        db.executemany('INSERT INTO analytehaspathway VALUES (?,?)',
                       [(rid, p) for p, ids in memberships.items() for rid in ids])


def test_similarity_rows_use_scope_specific_indices_counts_and_duplicates(tmp_path):
    path = tmp_path / 'ramp.sqlite'
    make_database(path)
    progress = []
    first = postprocess_pathway_similarity(path, progress=progress.append)
    assert first['rows'] == 6
    assert first['duplicate_pairs'] == 1
    assert first['analyte']['eligible_pathways'] == 4
    assert first['metabolite']['eligible_pathways'] == 4
    assert first['gene']['eligible_pathways'] == 5
    assert any('eligible pathways' in message for message in progress)
    with sqlite3.connect(path) as db:
        stored = db.execute('''SELECT pathwayRampId,analyte_blob,metabolite_blob,
            gene_blob,metabolite_count,gene_count
            FROM pathway_similarity ORDER BY pathwayRampId''').fetchall()
        assert [row[0] for row in stored] == ['P1', 'P2', 'P3', 'P4', 'P6', 'P7']
        assert zlib.decompress(stored[0][1]).decode() == '0,-|1,1000|1,818|1,333'
        assert zlib.decompress(stored[0][2]).decode() == '0,-|1,1000|1,1000'
        assert zlib.decompress(stored[0][3]).decode().startswith('0,-|1,1000|1,833')
        assert stored[2][2] is None and stored[2][4] is None
        assert stored[4][1] is None and stored[4][3] is None
        assert (stored[4][4], stored[4][5]) == (5, None)
        assert stored[5][1] is None and stored[5][2] is None
        assert (stored[5][4], stored[5][5]) == (None, 5)
        assert db.execute('SELECT * FROM pathway_duplicates').fetchall() == [('P1', 'P2')]
        manifest = json.loads(db.execute(
            "SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
        assert manifest['post_processing']['pathway_tables']['status'] == 'computed'
        assert manifest['pending_tables'] == []
        report = pathway_overlap_profile(db, 'Fixture')
        assert report['similarity_rows'] == 6
        assert report['duplicate_rows'] == 1
        assert report['similarity_example']['state'] == 'available'
        assert report['similarity_example']['blobs']['analyte_blob']['state'] == 'valid'
        assert report['duplicate_example']['state'] == 'available'
    second = postprocess_pathway_similarity(path)
    assert second['rows'] == first['rows']
    with sqlite3.connect(path) as db:
        assert db.execute('''SELECT pathwayRampId,analyte_blob,metabolite_blob,
            gene_blob,metabolite_count,gene_count
            FROM pathway_similarity ORDER BY pathwayRampId''').fetchall() == stored


def test_sparse_rows_omit_shared_pairs_that_round_to_zero():
    first = frozenset(range(1201))
    second = frozenset(range(1200, 2401))
    rows = list(sparse_rows(['P1', 'P2'], [first, second], label='fixture'))
    assert zlib.decompress(rows[0][2]).decode() == '0,-'
    assert zlib.decompress(rows[1][2]).decode() == '1,-'


def test_failed_rebuild_preserves_existing_pathway_tables(tmp_path, monkeypatch):
    path = tmp_path / 'ramp.sqlite'
    make_database(path)
    postprocess_pathway_similarity(path)
    with sqlite3.connect(path) as db:
        previous = db.execute('SELECT pathwayRampId,analyte_blob FROM pathway_similarity ORDER BY 1').fetchall()
        manifest = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0]

    def fail(db, *, progress):
        db.execute('DELETE FROM pathway_similarity')
        raise ValueError('fixture failure')

    monkeypatch.setattr('src.use_cases.ramp.pathway_similarity.populate_temp_tables', fail)
    with pytest.raises(ValueError, match='fixture failure'):
        postprocess_pathway_similarity(path)
    with sqlite3.connect(path) as db:
        assert db.execute('SELECT pathwayRampId,analyte_blob FROM pathway_similarity ORDER BY 1').fetchall() == previous
        assert db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0] == manifest
