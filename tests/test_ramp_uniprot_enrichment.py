import sqlite3

from src.use_cases.ramp.build_sqlite import export_sqlite
from src.use_cases.ramp.sqlite_protein_annotations import ProteinAnnotations
from tests.test_ramp_sqlite_export import resolved_fixture


def test_lookup_allowlist_and_namespace_spelling():
    record = {
        'primaryAccession': 'P1', 'secondaryAccessions': ['OLD1', 'OLD1'],
        'uniProtkbId': 'ENTRY_HUMAN',
        'genes': [{'geneName': {'value': 'GENE'}, 'synonyms': [{'value': 'ALIAS'}],
                   'orderedLocusNames': [{'value': 'LOCUS'}]}],
        'uniProtKBCrossReferences': [
            {'database': 'GeneID', 'id': '123'},
            {'database': 'HGNC', 'id': 'HGNC:456'},
            {'database': 'HGNC', 'id': '456'},
            {'database': 'Ensembl', 'id': 'ENST1.2', 'properties': [
                {'key': 'GeneId', 'value': 'ENSG1.3'},
                {'key': 'ProteinId', 'value': 'ENSP1.4'},
                {'key': 'Description', 'value': 'Not an identifier'}]},
            {'database': 'RefSeq', 'id': 'NP_123.1'},
        ],
    }
    assert set(ProteinAnnotations.lookup_identifiers(record)) == {
        'UniProtKB:P1', 'UniProtKB:OLD1', 'HGNC.SYMBOL:GENE', 'HGNC.SYMBOL:ALIAS',
        'NCBIGene:123', 'hgnc:456', 'Ensembl:ENST1', 'Ensembl:ENSG1', 'Ensembl:ENSP1',
    }
    assert ProteinAnnotations.lookup_identifiers({}) == ()


def snapshot(path):
    with sqlite3.connect(path) as db:
        return {name: sorted(db.execute('SELECT * FROM "' + name + '"').fetchall(), key=repr)
                for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table'")
                if name not in ('source', 'analytesynonym', 'ramp_export_metadata', 'version_info', 'db_version')}


def test_enrichment_shared_aliases_do_not_merge_groups_or_add_associations(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    baseline, enriched = tmp_path / 'baseline.sqlite', tmp_path / 'enriched.sqlite'
    export_sqlite(reader, baseline, gene_identity=identity, progress=lambda _: None)
    annotations = ProteinAnnotations.from_file(identity.input_file, identity.provenance)
    manifest = export_sqlite(reader, enriched, gene_identity=identity, protein_annotations=annotations,
                            progress=lambda _: None)
    assert snapshot(baseline) == snapshot(enriched)
    with sqlite3.connect(baseline) as old, sqlite3.connect(enriched) as new:
        # Every provider assertion is retained exactly; new aliases belong only to UniProt.
        query = "select * from source where dataSource != 'uniprot' order by sourceId,rampId,dataSource"
        assert old.execute(query).fetchall() == new.execute(query).fetchall()
        for alias in ('gene_symbol:SAME', 'entrez:9'):
            assert new.execute('select count(distinct rampId),group_concat(distinct dataSource) from source where sourceId=?',
                               (alias,)).fetchone() == (2, 'uniprot')
        assert new.execute("select count(*) from source where sourceId='uniprot:OLD1' and dataSource='uniprot'").fetchone() == (1,)
        assert new.execute("select count(*) from analytesynonym where source='uniprot' and Synonym='SAME'").fetchone() == (2,)
        assert old.execute("select count(*) from analytesynonym where source='uniprot'").fetchone() == (0,)
    assert manifest['gene_resolution'] == identity.manifest()
    assert manifest['protein_annotations']['lookup_alias_policy_version'] == 1
    assert manifest['protein_annotations']['file_accession_alias_counts_by_namespace']['HGNC.SYMBOL'] == 2


def test_collapsed_groups_dedupe_aliases_and_skip_unmatched_records(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['GeneIdentifier'].append({'id': 'NCBIGene:9', 'sources': ['HMDB']})
    annotations = ProteinAnnotations.from_file(identity.input_file, identity.provenance)
    annotations.aliases['UniProtKB:NOT_IN_GRAPH'] = ('HGNC.SYMBOL:DO_NOT_EXPORT',)
    path = tmp_path / 'collapsed.sqlite'
    export_sqlite(reader, path, gene_identity=identity, protein_annotations=annotations, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select count(*) from source where sourceId='gene_symbol:SAME' and dataSource='uniprot'").fetchone() == (1,)
        assert db.execute("select count(*) from source where sourceId='gene_symbol:DO_NOT_EXPORT'").fetchone() == (0,)
