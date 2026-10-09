import sqlite3

import pytest

from src.use_cases.ramp.build_sqlite import export_sqlite
from tests.test_ramp_sqlite_export import FixtureReader, edge, node, resolved_fixture


def test_source_only_artifact_keeps_lookup_and_required_analytes(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['HmdbMetaboliteOntologyEdge'][1]['details'][0]['source_id'] = 'HMDB:EXCLUDED_ONLY'
    path = tmp_path / 'source.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity, source_only=True,
                             progress=lambda _: None)
    with sqlite3.connect(path) as db:
        tables = {name for (name,) in db.execute("SELECT name FROM sqlite_master WHERE type='table'")}
        assert tables == {'analyte', 'source', 'analytesynonym', 'pathway', 'ontology',
                          'analytehaspathway', 'analytehasontology', 'catalyzed',
                          'metabolite_class', 'reaction', 'reaction2met',
                          'reaction2protein', 'reaction_ec_class', 'chem_props',
                          'version_info', 'ramp_export_metadata'}
        assert db.execute('SELECT count(*) FROM source').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM analytesynonym').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM analytehaspathway').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM analytehasontology').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM catalyzed').fetchone()[0] > 0
        assert db.execute('SELECT count(*) FROM metabolite_class').fetchone()[0] > 0
        for table in ('reaction', 'reaction2met', 'reaction2protein', 'reaction_ec_class'):
            assert db.execute(f'SELECT count(*) FROM {table}').fetchone()[0] > 0
        assert db.execute("SELECT chem_data_source,chem_source_id FROM chem_props").fetchall() == [
            ('hmdb', 'hmdb:HMDB1')]
        assert db.execute("SELECT data_source_id,data_source_snapshot_ids FROM version_info WHERE data_source_id='hmdb'").fetchone() == (
            'hmdb', '["hmdb:metabolites_xml:5.0"]')
        assert db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0]
        assert db.execute("SELECT count(*) FROM source WHERE sourceId='hmdb:EXCLUDED_ONLY'").fetchone() == (0,)
        assert db.execute('PRAGMA foreign_key_check').fetchall() == []
    assert manifest['scope'] == 'lookup_table_diagnostic'
    assert manifest['included_tables'] == ['analyte', 'source', 'analytesynonym', 'pathway',
                                          'ontology', 'analytehaspathway', 'analytehasontology',
                                          'catalyzed', 'metabolite_class', 'reaction', 'reaction2met',
                                          'reaction2protein', 'reaction_ec_class', 'chem_props',
                                          'version_info', 'ramp_export_metadata']
    assert manifest['release_ready'] is False


@pytest.mark.parametrize('source_only', [True, False])
@pytest.mark.parametrize('include_level4', [False, True])
def test_lipidmaps_level4_is_optional_without_losing_broader_classes_or_source(
        tmp_path, source_only, include_level4):
    reader, identity = resolved_fixture(tmp_path)
    reader.groups = lambda: [('CHEBI:1', 'HMDB:HMDB1', 'LIPIDMAPS:LM1'), ('HMDB:HMDB2',)]
    reader.data['MetaboliteIdentifier'].append(node('LIPIDMAPS:LM1'))
    levels = [('category', 'LipidMaps_category'),
              ('main', 'LipidMaps_main_class'),
              ('sub', 'LipidMaps_sub_class'),
              ('level4', 'LipidMaps_class_level4')]
    reader.data['MetaboliteClassificationTerm'] += [
        node(f'LMC:{name}', name=name, level_name=level, source='LipidMaps')
        for name, level in levels]
    reader.data['MetaboliteClassificationParentEdge'] += [
        edge(f'LMC:{parent}', f'LMC:{child}')
        for parent, child in zip(('category', 'main', 'sub'), ('main', 'sub', 'level4'))]
    reader.data['MetaboliteClassificationEdge'].append(
        edge('LIPIDMAPS:LM1', 'LMC:level4',
             details=[{'source': 'LipidMaps', 'source_id': 'LIPIDMAPS:LM1'}]))
    path = tmp_path / f'lipidmaps-{source_only}-{include_level4}.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity,
                             source_only=source_only,
                             include_lipidmaps_class_level4=include_level4,
                             progress=lambda _: None)
    with sqlite3.connect(path) as db:
        classes = dict(db.execute(
            "SELECT class_level_name,class_name FROM metabolite_class "
            "WHERE class_source_id='LIPIDMAPS:LM1'"))
        assert {level for _, level in levels[:3]} <= classes.keys()
        assert ('LipidMaps_class_level4' in classes) is include_level4
        assert db.execute("SELECT count(*) FROM source WHERE sourceId='LIPIDMAPS:LM1' ").fetchone()[0] > 0
        if source_only:
            from json import loads
            saved = loads(db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0])
            assert saved['include_lipidmaps_class_level4'] is include_level4
            assert saved['lipidmaps_class_level4_assertions_excluded'] == (0 if include_level4 else 1)
    assert manifest['include_lipidmaps_class_level4'] is include_level4
    assert manifest['lipidmaps_class_level4_assertions_excluded'] == (0 if include_level4 else 1)


def test_lipidmaps_level4_cli_defaults_off_and_can_be_enabled():
    from src.use_cases.ramp.build_sqlite import parser
    required = ['--stage-id', 's', '--output', 'x.sqlite']
    assert parser().parse_args(required).include_lipidmaps_class_level4 is False
    assert parser().parse_args(required + ['--include-lipidmaps-class-level4']).include_lipidmaps_class_level4 is True


@pytest.mark.parametrize('source_only', [True, False])
def test_hmdb_intermediate_health_terms_match_legacy_ontology(tmp_path, source_only):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['HmdbOntologyTerm'] += [
        node('ONT:health', name='Health condition', term_type='parent',
             ontology_type='Health condition'),
        node('ONT:cancer', name='Cancer', term_type='parent',
             ontology_type='Health condition'),
        node('ONT:diagnosis', name='Specific diagnosis', term_type='child',
             ontology_type='Health condition'),
        node('ONT:denied', name='Meningitis', term_type='child',
             ontology_type='Health condition'),
    ]
    reader.data['HmdbOntologyParentEdge'] += [
        edge('ONT:health', 'ONT:cancer'),
        edge('ONT:cancer', 'ONT:diagnosis'),
        edge('ONT:health', 'ONT:denied'),
    ]
    reader.data['HmdbMetaboliteOntologyEdge'] += [
        edge('HMDB:HMDB1', 'ONT:diagnosis'),
        edge('HMDB:HMDB1', 'ONT:denied'),
    ]
    path = tmp_path / 'health.sqlite'
    export_sqlite(reader, path, gene_identity=identity, source_only=source_only,
                  progress=lambda _: None)
    with sqlite3.connect(path) as db:
        names = {name for (name,) in db.execute(
            "SELECT commonName FROM ontology WHERE HMDBOntologyType='Health condition'")}
        assert names == {'Cancer', 'Specific diagnosis'}
        linked = {name for (name,) in db.execute(
            "SELECT o.commonName FROM analytehasontology a JOIN ontology o USING (rampOntologyId) "
            "WHERE a.rampCompoundId='RAMP_C_000000001' "
            "AND o.HMDBOntologyType='Health condition'")}
        assert linked == names


@pytest.mark.parametrize('source_only', [True, False])
@pytest.mark.parametrize('source_policy', ['legacy', 'expanded'])
def test_hmdb_source_terms_follow_selected_policy(tmp_path, source_only, source_policy):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['HmdbOntologyTerm'] += [
        node('ONT:source', name='Source', term_type='parent', ontology_type='Source'),
        node('ONT:plant', name='Plant', term_type='parent', ontology_type='Source'),
        node('ONT:microbe', name='Microbe', term_type='parent', ontology_type='Source'),
        node('ONT:fabaceae', name='Fabaceae', term_type='child', ontology_type='Source'),
        node('ONT:ecoli', name='Escherichia coli', term_type='child', ontology_type='Source'),
        node('ONT:food', name='Food', term_type='child', ontology_type='Source'),
    ]
    reader.data['HmdbOntologyParentEdge'] += [
        edge('ONT:source', 'ONT:plant'), edge('ONT:source', 'ONT:microbe'),
        edge('ONT:source', 'ONT:food'), edge('ONT:plant', 'ONT:fabaceae'),
        edge('ONT:microbe', 'ONT:ecoli'),
    ]
    reader.data['HmdbMetaboliteOntologyEdge'] += [
        edge('HMDB:HMDB1', 'ONT:fabaceae'),
        edge('HMDB:HMDB2', 'ONT:ecoli'),
        edge('HMDB:HMDB1', 'ONT:food'),
    ]
    path = tmp_path / 'source-ontology.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity, source_only=source_only,
                             source_ontology_policy=source_policy, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert {name for (name,) in db.execute(
            "SELECT commonName FROM ontology WHERE HMDBOntologyType='Source'")} == {
                'Plant', 'Microbe'} | ({'Fabaceae', 'Escherichia coli'}
                                     if source_policy == 'expanded' else set())
        assert set(db.execute(
            "SELECT a.rampCompoundId,o.commonName FROM analytehasontology a "
            "JOIN ontology o USING (rampOntologyId) WHERE o.HMDBOntologyType='Source'")) == {
                ('RAMP_C_000000001', 'Plant'), ('RAMP_C_000000002', 'Microbe')} | (
                    {('RAMP_C_000000001', 'Fabaceae'),
                     ('RAMP_C_000000002', 'Escherichia coli')}
                    if source_policy == 'expanded' else set())
    assert manifest['source_ontology_policy'] == source_policy


def test_source_ontology_policy_cli_defaults_to_legacy():
    from src.use_cases.ramp.build_sqlite import parser
    required = ['--stage-id', 's', '--output', 'x.sqlite']
    assert parser().parse_args(required).source_ontology_policy == 'legacy'
    assert parser().parse_args(required + ['--source-ontology-policy', 'expanded']).source_ontology_policy == 'expanded'


@pytest.mark.parametrize('source_only', [True, False])
@pytest.mark.parametrize('cutoff,hmdb_links,excluded', [(2, 0, 2), (3, 2, 0), (0, 2, 0)])
def test_pathway_cutoff_is_per_reported_id_and_source(
        tmp_path, source_only, cutoff, hmdb_links, excluded):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['PathwayIdentifier'].append(
        node('SMPDB:PW2', category='smpdb3', names=[{'value': 'Second path', 'source': 'HMDB'}]))
    reader.data['PathwayIdentifier'].append({
        'id': 'WP:1', 'sources': ['WikiPathways'], 'category': 'wiki',
        'names': [{'value': 'Wiki path', 'source': 'WikiPathways'}]})
    reader.data['MetabolitePathwayEdge'] += [
        edge('HMDB:HMDB2', 'SMPDB:PW1'),
        edge('HMDB:HMDB2', 'SMPDB:PW2'),
        edge('HMDB:HMDB2', 'WP:1', details=[{'source': 'WikiPathways'}]),
    ]
    path = tmp_path / f'cutoff-{cutoff}.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity, source_only=source_only,
                             pathway_association_cutoff=cutoff, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        rows = dict(db.execute(
            "SELECT pathwaySource,count(*) FROM analytehaspathway "
            "WHERE rampId='RAMP_C_000000002' GROUP BY pathwaySource"))
        assert rows.get('hmdb', 0) == hmdb_links
        assert rows.get('wiki', 0) == 1
        assert db.execute("SELECT count(*) FROM analytehaspathway WHERE rampId LIKE 'RAMP_G_%'").fetchone() == (1,)
        assert db.execute("SELECT count(*) FROM source WHERE sourceId='hmdb:HMDB2'").fetchone()[0] > 0
    assert manifest['pathway_association_cutoff'] == cutoff
    assert manifest['pathway_assertions_excluded'] == excluded
    assert manifest['pathway_source_ids_at_cutoff'] == (1 if cutoff == 2 else 0)


def test_pathway_association_cutoff_cli_validation():
    from src.use_cases.ramp.build_sqlite import parser
    args = parser().parse_args(['--stage-id', 's', '--output', 'x.sqlite'])
    assert args.pathway_association_cutoff == 25000
    with pytest.raises(SystemExit):
        parser().parse_args(['--stage-id', 's', '--output', 'x.sqlite',
                             '--pathway-association-cutoff', '-1'])


def test_pathway_cutoff_also_applies_to_wikipathways(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    for index in (1, 2):
        pathway_id = f'WP:{index}'
        reader.data['PathwayIdentifier'].append({
            'id': pathway_id, 'sources': ['WikiPathways'], 'category': 'wiki',
            'names': [{'value': f'Wiki path {index}', 'source': 'WikiPathways'}]})
        reader.data['MetabolitePathwayEdge'].append(
            edge('HMDB:HMDB2', pathway_id, details=[{'source': 'WikiPathways'}]))
    path = tmp_path / 'wiki-cutoff.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity, source_only=True,
                             pathway_association_cutoff=2, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute(
            "SELECT count(*) FROM analytehaspathway "
            "WHERE rampId='RAMP_C_000000002' AND pathwaySource='wiki'").fetchone() == (0,)
        assert db.execute(
            "SELECT count(*) FROM analytehaspathway "
            "WHERE rampId='RAMP_C_000000001' AND pathwaySource='hmdb'").fetchone() == (1,)
    assert manifest['pathway_source_ids_at_cutoff'] == 1
    assert manifest['pathway_assertions_excluded'] == 2


def test_pathway_cutoff_also_applies_to_gene_and_protein_source_ids(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['PathwayIdentifier'].append(
        node('SMPDB:PW2', category='smpdb3', names=[{'value': 'Second path', 'source': 'HMDB'}]))
    reader.data['GenePathwayEdge'].append(edge('NCBIGene:1', 'SMPDB:PW2'))
    reader.data['ProteinPathwayEdge'].append(edge('UniProtKB:P1', 'SMPDB:PW2'))
    path = tmp_path / 'gene-cutoff.sqlite'
    manifest = export_sqlite(reader, path, gene_identity=identity, source_only=True,
                             pathway_association_cutoff=2, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute(
            "SELECT count(*) FROM analytehaspathway WHERE rampId LIKE 'RAMP_G_%'").fetchone() == (0,)
        assert db.execute(
            "SELECT count(*) FROM source WHERE sourceId IN ('entrez:1','uniprot:P1')"
        ).fetchone()[0] > 0
    assert manifest['pathway_source_ids_at_cutoff'] == 2
    assert manifest['pathway_assertions_excluded'] == 4


def test_cutoff_does_not_hide_missing_pathway_nodes(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['MetabolitePathwayEdge'] += [
        edge('HMDB:HMDB2', 'SMPDB:PW1'), edge('HMDB:HMDB2', 'SMPDB:MISSING')]
    with pytest.raises(ValueError, match='Missing provider pathway'):
        export_sqlite(reader, tmp_path / 'invalid.sqlite', gene_identity=identity,
                      source_only=True, pathway_association_cutoff=2,
                      progress=lambda _: None)


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
