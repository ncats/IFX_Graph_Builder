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
        assert tables == {'analyte', 'source', 'analytesynonym', 'pathway', 'ontology',
                          'analytehaspathway', 'analytehasontology'}
        assert db.execute('SELECT count(*) FROM source').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM analytesynonym').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM analytehaspathway').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM analytehasontology').fetchone()[0] > 0
        assert db.execute("SELECT count(*) FROM source WHERE sourceId='hmdb:EXCLUDED_ONLY'").fetchone() == (0,)
        assert db.execute('PRAGMA foreign_key_check').fetchall() == []
    assert manifest['scope'] == 'lookup_table_diagnostic'
    assert manifest['included_tables'] == ['analyte', 'source', 'analytesynonym', 'pathway',
                                          'ontology', 'analytehaspathway', 'analytehasontology']
    assert manifest['release_ready'] is False


def test_wikipathways_symbols_keep_own_attribution_on_gene_and_protein(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    gene = reader.data['GeneIdentifier'][0]
    protein = reader.data['ProteinIdentifier'][0]
    for doc, symbol in ((gene, 'GENE_FROM_WIKI'), (protein, 'PROTEIN_FROM_WIKI')):
        doc['sources'].append('WikiPathways')
        doc['source_names'] = [{
            'value': symbol, 'source': 'WikiPathways', 'source_field': 'wp:bdbHgncSymbol'}]
    path = tmp_path / 'wiki-names.sqlite'
    export_sqlite(reader, path, gene_identity=identity, source_only=True,
                  progress=lambda _: None)
    with sqlite3.connect(path) as db:
        rows = db.execute("SELECT Synonym,source FROM analytesynonym WHERE source='wiki' ORDER BY Synonym").fetchall()
        assert rows == [('GENE_FROM_WIKI', 'wiki'), ('PROTEIN_FROM_WIKI', 'wiki')]
        assert db.execute("SELECT commonName FROM source WHERE sourceId='entrez:1' AND dataSource='wiki'").fetchone() == ('GENE_FROM_WIKI',)
        assert db.execute("SELECT commonName FROM source WHERE sourceId='uniprot:P1' AND dataSource='wiki'").fetchone() == ('PROTEIN_FROM_WIKI',)
        assert db.execute('PRAGMA foreign_key_check').fetchall() == []


def test_wikipathways_source_name_does_not_borrow_unattributed_gene_name(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['GeneIdentifier'][0]['sources'].append('WikiPathways')
    reader.data['ProteinIdentifier'][0]['sources'].append('WikiPathways')
    path = tmp_path / 'unnamed-wiki.sqlite'
    export_sqlite(reader, path, gene_identity=identity, source_only=True,
                  progress=lambda _: None)
    with sqlite3.connect(path) as db:
        rows = db.execute("SELECT sourceId,commonName FROM source WHERE dataSource='wiki' ORDER BY sourceId").fetchall()
        assert rows == [('entrez:1', None), ('uniprot:P1', None)]
        assert db.execute("SELECT count(*) FROM analytesynonym WHERE source='wiki'").fetchone() == (0,)
