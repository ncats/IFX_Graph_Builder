import gzip
import hashlib
import json
import sqlite3
from pathlib import Path

import pytest

from src.use_cases.ramp.sqlite_version_metadata import display_metadata
from src.use_cases.ramp.sqlite_protein_annotations import ProteinAnnotations
from src.use_cases.ramp.sqlite_hmdb_status import group_statuses
from tests.test_ramp_sqlite_export import FixtureReader, export_sqlite


def test_version_display_preserves_multiple_releases_and_snapshot_meaning():
    assert display_metadata('rhea', [{'dataset': 'reaction_bundle', 'version': '142', 'version_date': '2026-09-02'}]) == {
        'data_source_name': 'Rhea', 'data_source_url': 'https://www.rhea-db.org/',
        'data_source_version': 'Release 142 (2026-09-02)'}
    assert display_metadata('hmdb', [{'dataset': d, 'version': '5.0'} for d in ['proteins_xml', 'metabolites_xml']])['data_source_version'] == 'v5.0'
    value = display_metadata('chebi', [{'dataset': 'ontology_full', 'version': '255'}, {'dataset': 'three_star_sdf', 'version': '2026-09-09'}])['data_source_version']
    assert value == 'Ontology: Release 255; Structures: 2026-09-09'
    assert display_metadata('refmet', [{'dataset': 'metabolites_csv', 'version': 'sha256-abc', 'download_date': '2026-09-01'}])['data_source_version'] == 'downloaded 2026-09-01'
    assert display_metadata('pubchem', [{'dataset': 'cid_molecular_info', 'version': 'deps-abc', 'download_date': '2026-10-06'}])['data_source_version'] == 'downloaded 2026-10-06'
    assert display_metadata('pubchem', [{'dataset': 'cid_molecular_info', 'version': 'deps-abc'}])['data_source_version'] == 'Derived dataset; release not recorded'
    assert display_metadata('hmdb', [{'dataset': 'metabolites_xml'}])['data_source_version'] == 'Release not recorded'


def test_annotation_file_names_and_integrity(tmp_path):
    path = tmp_path / 'names.json.gz'
    with gzip.open(path, 'wt') as f:
        json.dump({'results': [
            {'primaryAccession': 'P1', 'secondaryAccessions': ['OLD'], 'proteinDescription': {'recommendedName': {'fullName': {'value': 'Protein one'}}}},
            {'primaryAccession': 'P2', 'proteinDescription': {'submissionNames': [{'fullName': {'value': 'Protein two'}}]}},
            {'primaryAccession': 'P3', 'proteinDescription': {'alternativeNames': [{'fullName': {'value': 'Protein three'}}]}},
        ]}, f)
    provenance = {'sha256': hashlib.sha256(path.read_bytes()).hexdigest(), 'snapshot_id': 'uniprot:human:test', 'file': path.name}
    annotations = ProteinAnnotations.from_file(path, provenance)
    assert annotations.names == {f'UniProtKB:P{i}': f'Protein {word}' for i, word in enumerate(['one', 'two', 'three'], 1)}
    assert annotations.name('UniProtKB:OLD') is None  # Reconciliation belongs to the resolver.
    with pytest.raises(ValueError, match='differs'):
        ProteinAnnotations.from_file(path, dict(provenance, sha256='wrong'))


def test_hmdb_priority_and_missing_status_values():
    mapping = {'HMDB:1': 'C1', 'HMDB:2': 'C1', 'HMDB:3': 'C2'}
    result, stats = group_statuses([{'id': 'HMDB:1', 'hmdb_status': 'expected'},
                                   {'id': 'HMDB:2', 'hmdb_status': 'quantified'},
                                   {'id': 'HMDB:3'}, {'id': 'HMDB:excluded', 'hmdb_status': 'quantified'}], mapping)
    assert result == {'C1': 'quantified'}
    assert stats['records_missing_field'] == 1
    with pytest.raises(ValueError, match='Unknown HMDB status'):
        group_statuses([{'id': 'HMDB:1', 'hmdb_status': 'unexpected'}], mapping)


def test_projection_status_propagation_and_protein_full_name(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][0]['hmdb_status'] = 'quantified'
    annotations = ProteinAnnotations({'UniProtKB:P1': 'Full UniProt protein name'})
    path = tmp_path / 'ramp.sqlite'
    export_sqlite(reader, path, protein_annotations=annotations, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select distinct priorityHMDBStatus from source where sourceId in ('hmdb:HMDB1','chebi:1')").fetchall() == [('quantified',)]
        assert db.execute('select protein_name from reaction2protein').fetchone() == ('Full UniProt protein name',)


def test_collapsed_proteins_keep_accession_specific_names(tmp_path):
    from tests.test_ramp_sqlite_export import resolved_fixture
    from src.use_cases.ramp.build_sqlite import export_sqlite as run_export
    reader, identity = resolved_fixture(tmp_path)
    reader.data['GeneIdentifier'].append({'id': 'NCBIGene:9', 'sources': ['HMDB']})
    annotations = ProteinAnnotations({'UniProtKB:P1': 'First protein', 'UniProtKB:P2': 'Second protein'})
    path = tmp_path / 'collapsed.sqlite'
    run_export(reader, path, gene_identity=identity, protein_annotations=annotations, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        rows = db.execute("select sourceId,rampId,commonName from source where dataSource='uniprot' order by sourceId").fetchall()
        assert [(r[0], r[2]) for r in rows] == [('uniprot:P1', 'First protein'), ('uniprot:P2', 'Second protein')]
        assert rows[0][1] == rows[1][1]


@pytest.mark.parametrize('source_only', [True, False])
def test_gene_display_name_prefers_available_protein_name_over_identifier(tmp_path, source_only):
    from tests.test_ramp_sqlite_export import resolved_fixture
    from src.use_cases.ramp.build_sqlite import export_sqlite as run_export
    reader, identity = resolved_fixture(tmp_path)
    reader.data['ProteinIdentifier'][0].pop('gene_name')
    path = tmp_path / 'gene-name.sqlite'
    run_export(reader, path, gene_identity=identity,
               protein_annotations=ProteinAnnotations({'UniProtKB:P1': 'Protein one'}),
               source_only=source_only, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute('''SELECT a.common_name FROM analyte a JOIN source s ON s.rampId=a.rampId
                             WHERE s.sourceId='entrez:1' AND a.type='gene' LIMIT 1''').fetchone() == ('Protein one',)


def test_gene_only_group_uses_name_from_its_canonical_uniprot_match(tmp_path):
    from tests.test_ramp_sqlite_export import resolved_fixture
    from src.use_cases.ramp.build_sqlite import export_sqlite as run_export
    reader, identity = resolved_fixture(tmp_path)
    reader.data['ProteinIdentifier'].clear()
    reader.data['ProteinPathwayEdge'].clear()
    reader.data['HmdbMetaboliteProteinAssociationEdge'].clear()
    reader.data['RheaProteinReactionEdge'].clear()
    path = tmp_path / 'gene-only.sqlite'
    run_export(reader, path, gene_identity=identity,
               protein_annotations=ProteinAnnotations({'UniProtKB:P1': 'Protein one'}),
               source_only=True, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute('''SELECT a.common_name FROM analyte a JOIN source s ON s.rampId=a.rampId
                             WHERE s.sourceId='entrez:1' AND a.type='gene' LIMIT 1''').fetchone() == ('Protein one',)


@pytest.mark.parametrize('source_only', [True, False])
def test_synonym_rows_preserve_annotation_owner_after_merging(tmp_path, source_only):
    from tests.test_ramp_sqlite_export import resolved_fixture
    from src.use_cases.ramp.build_sqlite import export_sqlite as run_export
    reader, identity = resolved_fixture(tmp_path)
    metabolite = reader.data['MetaboliteIdentifier'][0]
    metabolite['names'].append({'value': 'Water', 'source': 'LipidMaps'})
    metabolite['synonyms'].append({'value': 'Aqua', 'source': 'HMDB'})
    protein = reader.data['ProteinIdentifier'][0]
    protein['sources'].extend(['Reactome', 'RHEA'])
    protein['synonyms'] = ['HMDB protein alias']
    annotations = ProteinAnnotations(
        {'UniProtKB:P1': 'Pinned protein'},
        aliases={'UniProtKB:P1': ['HGNC.SYMBOL:SAME', 'HGNC.SYMBOL:ALIAS']},
    )
    path = tmp_path / 'synonyms.sqlite'
    run_export(reader, path, gene_identity=identity, protein_annotations=annotations,
               source_only=source_only, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        rows = db.execute('''SELECT Synonym,geneOrCompound,source,COUNT(*)
                             FROM analytesynonym GROUP BY Synonym COLLATE BINARY,geneOrCompound,source
                             ORDER BY geneOrCompound,source,Synonym''').fetchall()
        assert ('Water', 'compound', 'hmdb', 1) in rows
        assert ('Water', 'compound', 'lipidmaps', 1) in rows
        assert ('Aqua', 'compound', 'hmdb', 1) in rows
        assert ('SAME', 'gene', 'hmdb', 2) in rows  # Same name, two RaMP_G groups.
        assert ('Protein one', 'gene', 'hmdb', 1) in rows
        assert ('HMDB protein alias', 'gene', 'hmdb', 1) in rows
        assert ('Pinned protein', 'gene', 'uniprot', 1) in rows
        assert ('ALIAS', 'gene', 'uniprot', 1) in rows
        assert ('SAME', 'gene', 'uniprot', 1) in rows
        assert not any(source in ('reactome', 'rhea') and kind == 'gene'
                       for _, kind, source, _ in rows)


def test_hmdb_status_survives_source_parser_and_primary_node():
    import xml.etree.ElementTree as ET
    from src.input_adapters.metabolite_harmonization.hmdb import HmdbMetaboliteEquivalenceAdapter as Adapter
    record = Adapter._parse_metabolite(ET.fromstring('<metabolite xmlns="http://www.hmdb.ca"><accession>HMDB1</accession><status>quantified</status><name>Water</name></metabolite>'))
    assert Adapter._primary_node('HMDB:HMDB1', record).hmdb_status == 'quantified'


def test_status_schema_allows_secondary_identifiers_without_source_status(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][0]['hmdb_status'] = 'quantified'
    messages = []
    manifest = export_sqlite(reader, tmp_path / 'status.sqlite', progress=messages.append)
    assert manifest['hmdb_status']['records_missing_field'] > 0
    assert manifest['hmdb_status']['graph_schema_supports_status']
    assert not any('rebuild the graph' in m for m in messages)


def test_export_rejects_graph_without_hmdb_status_schema(tmp_path):
    reader = FixtureReader()
    reader.schemas = {'MetaboliteIdentifier': {'fields': {}}}
    with pytest.raises(ValueError, match='lacks hmdb_status; rebuild the graph'):
        export_sqlite(reader, tmp_path / 'old-graph.sqlite', progress=lambda _: None)


def test_kegg_pathways_use_legacy_attribution_with_hmdb_provenance(tmp_path):
    reader = FixtureReader()
    reader.data['PathwayIdentifier'].append({
        'id': 'KEGG.PATHWAY:map00010', 'sources': ['HMDB'], 'category': 'kegg',
        'names': [{'value': 'Glycolysis / Gluconeogenesis', 'source': 'HMDB'}]})
    for collection, start in [('MetabolitePathwayEdge', 'HMDB:HMDB1'),
                              ('GenePathwayEdge', 'NCBIGene:1'),
                              ('ProteinPathwayEdge', 'UniProtKB:P1')]:
        reader.data[collection].append({'start_id': start, 'end_id': 'KEGG.PATHWAY:map00010',
                                       'details': [{'source': 'HMDB', 'source_id': start, 'gene_id': start, 'protein_id': start}]})
    path = tmp_path / 'kegg.sqlite'
    manifest = export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select type,pathwayCategory,pathwayName from pathway where sourceId='map00010'").fetchone() == ('kegg', 'kegg', 'Glycolysis / Gluconeogenesis')
        assert db.execute("select ap.pathwaySource,count(*) from analytehaspathway ap join pathway p using(pathwayRampId) where p.sourceId='map00010' group by ap.pathwaySource").fetchall() == [('kegg', 3)]
        assert db.execute("select type from pathway where sourceId='PW1'").fetchone() == ('hmdb',)
        assert db.execute("select data_source_version,data_source_snapshot_ids from version_info where data_source_id='kegg'").fetchone() == ('From HMDB (v5.0)', '["hmdb:metabolites_xml:5.0"]')
    assert manifest['graph_build']['registry_datasets'][0]['source'] == 'hmdb'
