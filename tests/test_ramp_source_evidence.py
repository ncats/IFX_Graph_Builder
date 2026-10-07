import sqlite3
import pytest
from tests.test_ramp_sqlite_export import FixtureReader, export_sqlite, resolved_fixture
from src.use_cases.ramp.build_sqlite import export_sqlite as run_export
from src.use_cases.ramp.sqlite_protein_annotations import ProteinAnnotations
from src.use_cases.ramp.sqlite_projection import name_candidates
from src.use_cases.ramp.sqlite_source_evidence import representative_compound_names


DETAIL_TYPES = {
    'MetabolitePathwayEdge': 'source_id', 'HmdbMetaboliteOntologyEdge': 'source_id',
    'HmdbMetaboliteProteinAssociationEdge': 'source_id', 'MetaboliteClassificationEdge': 'source_id',
    'GenePathwayEdge': 'gene_id', 'ProteinPathwayEdge': 'protein_id',
}


def test_metabolite_name_falls_back_to_most_common_source_id_when_names_missing():
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE _raw_source (rampId TEXT, sourceId TEXT, commonName TEXT, dataSource TEXT, geneOrCompound TEXT)')
        db.executemany('INSERT INTO _raw_source VALUES (?,?,?,?,?)', [
            ('C1', 'CAS:123', None, 'hmdb', 'compound'),
            ('C1', 'CAS:123', 'NA', 'wiki', 'compound'),
            ('C1', 'HMDB:X', 'None', 'hmdb', 'compound'),
            ('G1', 'UniProtKB:P1', 'EGFR', 'uniprot', 'gene'),
        ])
        assert representative_compound_names(db) == {'C1': 'CAS:123'}


def test_metabolite_name_ties_prefer_shorter_name_then_lexical_order():
    with sqlite3.connect(':memory:') as db:
        db.execute('CREATE TABLE _raw_source (rampId TEXT, sourceId TEXT, commonName TEXT, dataSource TEXT, geneOrCompound TEXT)')
        db.executemany('INSERT INTO _raw_source VALUES (?,?,?,?,?)', [
            ('C1', 'HMDB:1', 'Long chemical name', 'hmdb', 'compound'),
            ('C1', 'HMDB:2', 'Glucose', 'hmdb', 'compound'),
            ('C2', 'HMDB:3', 'Water', 'hmdb', 'compound'),
            ('C2', 'HMDB:4', 'Aqua!', 'hmdb', 'compound'),
        ])
        assert representative_compound_names(db) == {'C1': 'Glucose', 'C2': 'Aqua!'}


def test_all_evidence_ids_reach_source_without_changing_ramp_groups(tmp_path):
    reader = FixtureReader()
    expected = {}
    for i, (collection, field) in enumerate(DETAIL_TYPES.items()):
        raw = f'HMDB:reported{i}' if field == 'source_id' else f'UniProtKB:reported{i}'
        reader.data[collection][0]['details'][0][field] = raw
        expected[raw.replace('HMDB:', 'hmdb:').replace('UniProtKB:', 'uniprot:')] = 'hmdb'
    reader.data['RheaMetaboliteReactionEdge'][0]['source_id'] = 'CHEBI:reported7'
    reader.data['RheaProteinReactionEdge'][0]['source_id'] = 'UniProtKB:reported8'
    expected.update({'chebi:reported7': 'rhea', 'uniprot:reported8': 'rhea'})
    path = tmp_path / 'evidence.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        for identifier, src in expected.items():
            assert db.execute('select count(*) from source where sourceId=? and dataSource=?', (identifier, src)).fetchone() == (1,)
        assert db.execute('select distinct class_source_id from metabolite_class').fetchall() == [('hmdb:reported3',)]
        assert db.execute('select met_source_id from reaction2met').fetchone() == ('chebi:reported7',)
        assert db.execute('select uniprot from reaction2protein').fetchone() == ('uniprot:reported8',)
        assert db.execute("select count(*) from analyte where type='gene'").fetchone() == (3,)
        assert db.execute("select count(*) from analyte where type='compound'").fetchone() == (2,)
        assert db.execute("select count(*) from source where sourceId='hmdb:HMDBPP1' and geneOrCompound='gene'").fetchone() == (1,)
        for table in ('analytehaspathway','analytehasontology','catalyzed'):
            assert 'source_id' not in {r[1] for r in db.execute(f'pragma table_info({table})')}


@pytest.mark.parametrize('collection,field', [*DETAIL_TYPES.items(), ('RheaMetaboliteReactionEdge','source_id'), ('RheaProteinReactionEdge','source_id')])
def test_missing_source_evidence_is_not_reconstructed_from_endpoint(tmp_path, collection, field):
    reader = FixtureReader()
    target = reader.data[collection][0]
    if collection in DETAIL_TYPES:
        target['details'][0].pop(field)
    else:
        target.pop(field)
    with pytest.raises(ValueError, match='refresh or repair'):
        export_sqlite(reader, tmp_path / 'invalid.sqlite', progress=lambda _: None)
    assert not (tmp_path / 'invalid.sqlite').exists()


def test_source_evidence_cannot_reconnect_curated_metabolite_split(tmp_path):
    reader = FixtureReader()
    reader.data['MetabolitePathwayEdge'][0]['details'][0]['source_id'] = 'HMDB:HMDB2'
    with pytest.raises(ValueError, match='contradicts metabolite groups'):
        export_sqlite(reader, tmp_path / 'split.sqlite', progress=lambda _: None)


def test_node_and_repeated_edge_evidence_yield_one_named_source_row(tmp_path):
    reader = FixtureReader()
    reader.data['MetabolitePathwayEdge'].append(reader.data['MetabolitePathwayEdge'][0])
    path = tmp_path / 'dedup.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select commonName from source where sourceId='hmdb:HMDB1' and dataSource='hmdb'").fetchall() == [('Water',)]


def test_ontology_edge_excluded_from_sqlite_does_not_create_source_lookup(tmp_path):
    reader = FixtureReader()
    reader.data['HmdbMetaboliteOntologyEdge'][1]['details'][0]['source_id'] = 'HMDB:EXCLUDED_ONLY'
    path = tmp_path / 'ontology.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select count(*) from source where sourceId='hmdb:EXCLUDED_ONLY'").fetchone() == (0,)


def test_secondary_accession_name_uses_its_resolver_match_not_whole_group(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['ProteinIdentifier'].append({'id': 'UniProtKB:OLD1', 'sources': ['Rhea']})
    reader.data['RheaProteinReactionEdge'][0].update(start_id='UniProtKB:OLD1', source_id='UniProtKB:OLD1')
    path = tmp_path / 'secondary.sqlite'
    run_export(reader, path, gene_identity=identity,
               protein_annotations=ProteinAnnotations({'UniProtKB:P1': 'Protein one'}), progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute('select uniprot,protein_name from reaction2protein').fetchone() == ('uniprot:OLD1', 'Protein one')
        assert db.execute("select count(distinct rampId) from source where sourceId in ('uniprot:OLD1','uniprot:P1')").fetchone() == (1,)


def test_merged_pathway_edge_preserves_each_reporting_provider(tmp_path):
    reader = FixtureReader()
    reader.data['PathwayIdentifier'][0]['sources'].append('Reactome')
    reader.data['PathwayIdentifier'][0]['names'].append({'value': 'Reactome path', 'source': 'Reactome'})
    reader.data['MetabolitePathwayEdge'][0]['details'].append({'source': 'Reactome', 'source_id': 'CHEBI:1'})
    path = tmp_path / 'providers.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select count(*) from source where sourceId='chebi:1' and dataSource='reactome'").fetchone() == (1,)
        assert db.execute("select distinct pathwaySource from analytehaspathway where rampId='RAMP_C_000000001' order by pathwaySource").fetchall() == [('hmdb',), ('reactome',)]


def test_rhea_edge_uses_its_own_name_and_lipidmaps_prefers_abbreviation(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][1]['names'].extend([
        'Unattributed name',
        {'value': 'long systematic lipid name', 'source': 'LipidMaps', 'source_field': 'NAME'},
        {'value': 'FA 40:7;O3', 'source': 'LipidMaps', 'source_field': 'ABBREVIATION'},
    ])
    reader.data['MetaboliteIdentifier'][1]['sources'].append('LipidMaps')
    assert min(n for n in name_candidates(reader.data['MetaboliteIdentifier'][1])
               if n[1] == 'lipidmaps')[2] == 'FA 40:7;O3'
    reader.data['MetaboliteIdentifier'][1]['names'].remove('Unattributed name')
    path = tmp_path / 'source-names.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select commonName from source where sourceId='chebi:1' and dataSource='rhea'").fetchone() == ('Water',)
        assert db.execute("select commonName from source where sourceId='chebi:1' and dataSource='lipidmaps'").fetchone() == ('FA 40:7;O3',)
        assert db.execute("select commonName from source where sourceId='chebi:1' and dataSource='hmdb'").fetchone() == (None,)


def test_lipidmaps_source_name_prefers_abbreviation_then_own_synonym(tmp_path):
    reader = FixtureReader()
    lipid = reader.data['MetaboliteIdentifier'][1]
    lipid['sources'].append('LipidMaps')
    lipid['names'].extend([
        {'value': 'Long chemical name', 'source': 'LipidMaps', 'source_field': 'NAME'},
        {'value': 'FA 40:7;O3', 'source': 'LipidMaps', 'source_field': 'ABBREVIATION'},
    ])
    lipid['synonyms'] = [{'value': 'Acetylenic acids', 'source': 'LipidMaps',
                          'source_field': 'SYNONYMS'}]
    path = tmp_path / 'lipidmaps-name.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT commonName FROM source WHERE sourceId='chebi:1' AND dataSource='lipidmaps'").fetchone() == ('FA 40:7;O3',)
        assert db.execute("SELECT Synonym FROM analytesynonym WHERE source='lipidmaps' ORDER BY Synonym").fetchall() == [
            ('Acetylenic acids',), ('FA 40:7;O3',), ('Long chemical name',)]
    lipid['names'] = [name for name in lipid['names']
                      if not isinstance(name, dict) or name.get('source_field') != 'ABBREVIATION']
    fallback_path = tmp_path / 'lipidmaps-synonym-fallback.sqlite'
    export_sqlite(reader, fallback_path, progress=lambda _: None)
    with sqlite3.connect(fallback_path) as db:
        assert db.execute("SELECT commonName FROM source WHERE sourceId='chebi:1' AND dataSource='lipidmaps'").fetchone() == ('Acetylenic acids',)


def test_metabolite_synonym_without_source_fails_export(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][0]['synonyms'].append({'value': 'Unattributed'})
    with pytest.raises(ValueError, match='Metabolite name lacks source attribution: HMDB:HMDB1'):
        export_sqlite(reader, tmp_path / 'missing-provenance.sqlite', progress=lambda _: None)
