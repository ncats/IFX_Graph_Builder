import sqlite3

from src.use_cases.ramp.build_sqlite import export_sqlite
from tests.test_ramp_sqlite_export import FixtureReader, resolved_fixture


def test_source_only_artifact_keeps_lookup_and_required_analytes(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['RheaReaction'][0]['status'] = 'UNSUPPORTED_FOR_FULL_EXPORT'
    reader.data['HmdbMetaboliteOntologyEdge'][1]['details'][0]['source_id'] = 'HMDB:EXCLUDED_ONLY'
    path = tmp_path / 'source.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity, source_only=True,
                             progress=lambda _: None)
    with sqlite3.connect(path) as db:
        tables = {name for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        assert tables == {'analyte', 'source'}
        assert db.execute('SELECT count(*) FROM source').fetchone()[0] > 0
        assert db.execute("SELECT count(*) FROM source WHERE sourceId='hmdb:EXCLUDED_ONLY'").fetchone() == (0,)
        assert db.execute('PRAGMA foreign_key_check').fetchall() == []
    assert manifest['scope'] == 'source_table_diagnostic'
    assert manifest['release_ready'] is False
