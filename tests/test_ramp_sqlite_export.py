import copy
import gzip
import hashlib
import json
import sqlite3
from types import SimpleNamespace

import pytest

from src.core.curations import (
    METABOLITE_RECORD_PROPERTIES, batch_key, manifest_key, payload_sha256,
    replay_curation_snapshot, resolve_curation_type,
)
from src.use_cases.ramp.build_sqlite import export_sqlite as run_export, main
from src.use_cases.ramp.sqlite_gene_identity import GeneIdentity, FILE_NAME
from src.use_cases.ramp.sqlite_stage import StageReader
from src.core.graph_build_identity import source_build_fingerprint


def node(identifier, **kwargs):
    return {"id": identifier, "sources": ["HMDB\t5.0\t2021-11-17\t2026-06-30"], **kwargs}


def edge(start, end, **kwargs):
    kwargs.setdefault('source_id', start)
    details = kwargs.setdefault('details', [{'source': 'HMDB'}])
    for detail in details:
        detail.setdefault('source_id', start)
        detail.setdefault('gene_id', start)
        detail.setdefault('protein_id', start)
        if start.startswith('UniProtKB:') or end.startswith('UniProtKB:'):
            detail.setdefault('hmdb_protein_accession', 'HMDBP' + (start if start.startswith('UniProtKB:') else end).split(':')[1])
    return {"start_id": start, "end_id": end, **kwargs}


def export_sqlite(reader, output, **kwargs):
    """Existing projection fixtures deliberately have no UniProt matches."""
    identity = GeneIdentity(SimpleNamespace(resolve_internal=lambda nodes: {n.id: [] for n in nodes}),
                            {'snapshot_id': 'uniprot:human:test', 'file': FILE_NAME, 'sha256': 'fixture'})
    return run_export(reader, output, gene_identity=identity, **kwargs)


class FixtureReader:
    stage = {"id": "HarmonizationStage:test", "summary": {"clique_count": 2}}
    metadata = {"registry_datasets": [{"source": "hmdb", "dataset": "metabolites_xml",
                "version": "5.0", "snapshot_id": "hmdb:metabolites_xml:5.0",
                "manifest_uri": "s3://registry/hmdb/manifest.yaml"},
                {'source': 'uniprot', 'dataset': 'human', 'version': 'test',
                 'snapshot_id': 'uniprot:human:test', 'files': [{'path': FILE_NAME, 'sha256': 'fixture'}]}]}
    revisions = {"MetaboliteIdentifier": "revision-1"}

    def __init__(self):
        self.schemas = {'MetaboliteIdentifier': {'fields': {'hmdb_status': 'str'}}}
        self.data = {
            "MetaboliteIdentifier": [
                node("HMDB:HMDB1", names=[{"value": "Water", "source": "HMDB"}],
                     synonyms=[{"value": "Aqua", "source": "HMDB"}],
                     chem_props=[{"source": "HMDB", "source_id": "HMDB:HMDB1", "mw": "18.015",
                                  "molecular_formula": "H2O", "inchi_key": "REPORTED-KEY",
                                  "derived_inchi_key": "DO-NOT-EXPORT"}]),
                node("CHEBI:1", names=[{"value": "Oxidane", "source": "ChEBI"}]),
                node("HMDB:HMDB2", names=[{"value": "Singleton", "source": "HMDB"}]),
                node("HMDB:REMOVED", names=[{"value": "Excluded", "source": "HMDB"}]),
            ],
            "GeneIdentifier": [node("NCBIGene:1", names=["SAME"])],
            "ProteinIdentifier": [node("UniProtKB:P1", gene_name="SAME", name="Protein one", is_reviewed=True),
                                  node("UniProtKB:P2", gene_name="SAME", name="Protein two")],
            "PathwayIdentifier": [node("SMPDB:PW1", category="smpdb3", names=[{"value": "Path", "source": "HMDB"}])],
            "MetabolitePathwayEdge": [edge("HMDB:HMDB1", "SMPDB:PW1", details=[{"source": "HMDB"}]),
                                      edge("CHEBI:1", "SMPDB:PW1", details=[{"source": "HMDB"}]),
                                      edge("HMDB:REMOVED", "SMPDB:PW1", details=[{"source": "HMDB"}])],
            "GenePathwayEdge": [edge("NCBIGene:1", "SMPDB:PW1", details=[{"source": "HMDB"}])],
            "ProteinPathwayEdge": [edge("UniProtKB:P1", "SMPDB:PW1", details=[{"source": "HMDB"}])],
            "HmdbMetaboliteProteinAssociationEdge": [edge("HMDB:HMDB1", "UniProtKB:P1")],
            "MetaboliteClassificationTerm": [node("CLASS:root", name="Root", level_name="superclass", source="HMDB"),
                                              node("CLASS:leaf", name="Leaf", level_name="class", source="HMDB")],
            "MetaboliteClassificationParentEdge": [edge("CLASS:root", "CLASS:leaf")],
            "MetaboliteClassificationEdge": [edge("HMDB:HMDB1", "CLASS:leaf")],
            "HmdbOntologyTerm": [node("ONT:root", name="Biofluid and excreta", term_type="parent", ontology_type="Biofluid and excreta"),
                                 node("ONT:saliva", name="Saliva", term_type="child", ontology_type="Biofluid and excreta"),
                                 node("ONT:plasma", name="Plasma", term_type="child", ontology_type="Biofluid and excreta")],
            "HmdbOntologyParentEdge": [edge("ONT:root", "ONT:saliva"), edge("ONT:root", "ONT:plasma")],
            "HmdbMetaboliteOntologyEdge": [edge("HMDB:HMDB1", "ONT:saliva"), edge("HMDB:HMDB1", "ONT:plasma")],
            "RheaReaction": [node("RHEA:1", status="Approved", is_transport=False, direction="UN", label="water reaction")],
            "BiologicalRole": [node("CHEBI:23357", name="cofactor")],
            "RheaReactionClass": [node("EC:1.-.-.-", name="Root EC", ec_level=1), node("EC:1.1.1.1", name="Leaf EC", ec_level=4)],
            "RheaReactionClassParentEdge": [edge("EC:1.1.1.1", "EC:1.-.-.-")],
            "RheaReactionReactionClassEdge": [edge("RHEA:1", "EC:1.1.1.1")],
            "RheaMetaboliteReactionEdge": [edge("CHEBI:1", "RHEA:1", side="left", name="Water")],
            "RheaProteinReactionEdge": [edge("UniProtKB:P1", "RHEA:1")],
        }
        self.checked = False

    def records(self, collection):
        # StageReader provides stable ordering; emulate its _key ordering.
        return sorted(copy.deepcopy(self.data.get(collection, [])), key=lambda d: d.get("id", d.get("start_id", "")))

    def groups(self):
        return [("CHEBI:1", "HMDB:HMDB1"), ("HMDB:HMDB2",)]

    def verify_unchanged(self):
        self.checked = True


def test_base_export_preserves_consumer_contract_and_does_not_invent_identity(tmp_path):
    reader = FixtureReader()
    path = tmp_path / "ramp.sqlite"
    manifest = export_sqlite(reader, path, progress=lambda _: None)
    assert reader.checked
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT count(*) FROM analyte WHERE type='compound'").fetchone() == (2,)
        assert db.execute("SELECT count(*) FROM analyte WHERE type='gene'").fetchone() == (3,)
        assert db.execute("SELECT common_name FROM analyte WHERE rampId='RAMP_C_000000001'").fetchone() == ("Water",)
        assert db.execute("SELECT count(*) FROM analytehaspathway WHERE rampId='RAMP_C_000000001'").fetchone() == (1,)
        assert db.execute("SELECT inchi_key,mw FROM chem_props WHERE chem_source_id='hmdb:HMDB1'").fetchone() == ("REPORTED-KEY", 18.015)
        assert db.execute("SELECT class_name FROM metabolite_class ORDER BY class_name").fetchall() == [("Leaf",), ("Root",)]
        assert db.execute("SELECT commonName FROM ontology").fetchall() == [("Saliva",)]
        assert db.execute("SELECT count(*) FROM analytehasontology").fetchone() == (1,)
        assert db.execute("SELECT uniprot,is_reviewed FROM reaction2protein").fetchone() == ("uniprot:P1", 1)
        assert db.execute("SELECT substrate_product,is_cofactor FROM reaction2met").fetchone() == (0, 0)
        assert db.execute("SELECT has_human_prot,only_human_mets FROM reaction").fetchone() == (-1, -1)
        assert db.execute("PRAGMA foreign_key_check").fetchall() == []
        assert db.execute("SELECT name FROM sqlite_master WHERE name='reaction_protein2met'").fetchall() == []
        assert db.execute("SELECT name FROM sqlite_master WHERE name LIKE '_raw_%'").fetchall() == []
        assert json.loads(db.execute("SELECT data_source_snapshot_ids FROM version_info").fetchone()[0]) == ["hmdb:metabolites_xml:5.0"]
        assert json.loads(db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()[0]) == manifest
        # Representative R data-access join, using existing column names.
        assert db.execute("SELECT a.common_name, cp.iso_smiles FROM reaction2met rm JOIN analyte a ON a.rampId=rm.ramp_cmpd_id LEFT JOIN chem_props cp ON cp.chem_source_id=rm.met_source_id").fetchone()[0] == "Water"
    assert manifest["release_ready"] is False
    assert manifest["excluded_metabolite_edges"] == {"MetabolitePathwayEdge": 1}


@pytest.mark.parametrize('source_only', [True, False])
def test_reaction_cofactor_uses_reported_chebi_id_and_proteins_use_un_only(tmp_path, source_only):
    reader = FixtureReader()
    reader.data['BiologicalRole'].append(node('CHEBI:23354', name='coenzyme'))
    reader.data['IsAEdge'] = [edge('CHEBI:23354', 'CHEBI:23357'),
                              edge('CHEBI:3', 'CHEBI:1')]
    reader.data['HasBiologicalRoleEdge'] = [edge('CHEBI:1', 'CHEBI:23354')]
    reader.data['MetaboliteIdentifier'] += [node('CHEBI:2'), node('CHEBI:3')]
    reader.groups = lambda: [('CHEBI:1', 'CHEBI:2', 'CHEBI:3', 'HMDB:HMDB1'), ('HMDB:HMDB2',)]
    reader.data['RheaReaction'].append(
        node('RHEA:2', status='Approved', is_transport=False, direction='LR', label='directional'))
    reader.data['RheaMetaboliteReactionEdge'] += [
        edge('CHEBI:2', 'RHEA:1', side='right'),
        edge('CHEBI:3', 'RHEA:1', side='right'),
        edge('CHEBI:1', 'RHEA:2', side='left'),
    ]
    reader.data['RheaProteinReactionEdge'].append(edge('UniProtKB:P1', 'RHEA:2'))
    path = tmp_path / 'reaction.sqlite'
    manifest = export_sqlite(reader, path, source_only=source_only, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute(
            'SELECT met_source_id,is_cofactor FROM reaction2met ORDER BY rxn_source_id,met_source_id'
        ).fetchall() == [('chebi:1', 1), ('chebi:2', 0), ('chebi:3', 1), ('chebi:1', 1)]
        assert db.execute('SELECT rxn_source_id FROM reaction2protein').fetchall() == [('rhea:1',)]
        assert db.execute('SELECT rxn_source_id FROM reaction ORDER BY rxn_source_id').fetchall() == [
            ('rhea:1',), ('rhea:2',)]
        assert db.execute('SELECT count(*) FROM source WHERE sourceId=? AND dataSource=?',
                          ('uniprot:P1', 'rhea')).fetchone()[0] > 0
    assert manifest['cofactor_role_count'] == 2
    assert manifest['cofactor_chemical_id_count'] == 2
    assert manifest['cofactor_reaction_assertions'] == 3
    assert manifest['non_un_protein_assertions_excluded'] == 1


def test_reaction_cofactor_role_is_required(tmp_path):
    reader = FixtureReader()
    reader.data['BiologicalRole'] = []
    with pytest.raises(ValueError, match='CHEBI:23357.*missing'):
        export_sqlite(reader, tmp_path / 'missing-role.sqlite', progress=lambda _: None)


@pytest.mark.parametrize('source_only', [True, False])
def test_metabolite_analyte_name_uses_most_common_source_name(tmp_path, source_only):
    reader = FixtureReader()
    for i, spelling in enumerate(('D-Glucose', 'd-glucose', 'D-GLUCOSE'), 1):
        reader.data['RheaMetaboliteReactionEdge'].append(
            edge('CHEBI:1', 'RHEA:1', source_id=f'CHEBI:vote{i}', name=spelling, side='left')
        )
    path = tmp_path / 'majority.sqlite'
    export_sqlite(reader, path, source_only=source_only, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT common_name FROM analyte WHERE rampId='RAMP_C_000000001'").fetchone() == ('D-Glucose',)
        assert db.execute("SELECT count(*) FROM source WHERE rampId='RAMP_C_000000001' AND lower(commonName)='d-glucose'").fetchone() == (3,)


def test_shuffled_input_keeps_ids_and_rows(tmp_path):
    left, right = FixtureReader(), FixtureReader()
    for values in right.data.values():
        values.reverse()
    export_sqlite(left, tmp_path / "a.sqlite", progress=lambda _: None)
    export_sqlite(right, tmp_path / "b.sqlite", progress=lambda _: None)
    with sqlite3.connect(tmp_path / "a.sqlite") as a, sqlite3.connect(tmp_path / "b.sqlite") as b:
        for table in ("analyte", "source", "chem_props", "pathway", "reaction"):
            assert sorted(a.execute(f"SELECT * FROM {table}").fetchall(), key=repr) == sorted(b.execute(f"SELECT * FROM {table}").fetchall(), key=repr)


def resolved_fixture(tmp_path, file_name=FILE_NAME):
    path = tmp_path / file_name
    entries = []
    for accession, genes in [('P1', ['1', '9']), ('P2', ['2', '9'])]:
        entries.append({'primaryAccession': accession, 'uniProtkbId': accession + '_HUMAN',
                        'entryType': 'UniProtKB reviewed (Swiss-Prot)' if accession == 'P1' else 'UniProtKB unreviewed (TrEMBL)',
                        'secondaryAccessions': ['OLD1'] if accession == 'P1' else [],
                        'genes': [{'geneName': {'value': 'SAME'}}],
                        'uniProtKBCrossReferences': [{'database': 'GeneID', 'id': g} for g in genes]})
    with gzip.open(path, 'wt') as out:
        json.dump({'results': entries}, out)
    reader = FixtureReader()
    reader.metadata = copy.deepcopy(reader.metadata)
    reader.metadata['registry_datasets'][1]['files'][0]['path'] = file_name
    reader.metadata['registry_datasets'][1]['files'][0]['sha256'] = hashlib.sha256(path.read_bytes()).hexdigest()
    def dataset_file(name):
        assert name == file_name
        return path
    dataset = SimpleNamespace(snapshot_id='uniprot:human:test', file=dataset_file)
    requested = []
    def resolve(pin):
        requested.append(pin)
        return dataset
    identity = GeneIdentity.from_stage(reader.metadata, SimpleNamespace(resolve=resolve), file_name=file_name)
    assert requested == ['uniprot:human:test']
    return reader, identity


def test_all_human_default_and_explicit_reviewed_file(tmp_path):
    from src.use_cases.ramp.build_sqlite import parser
    args = parser().parse_args(['--stage-id', 'test', '--output', 'test.sqlite'])
    assert args.uniprot_file == FILE_NAME == 'uniprot-human.json.gz'
    _, identity = resolved_fixture(tmp_path)
    key = ('ProteinIdentifier', 'UniProtKB:P2')
    assert identity.groups([key]).input_to_group[key] == ('UniProtKB', 'UniProtKB:P2')
    _, alternate = resolved_fixture(tmp_path, 'uniprot-human-reviewed.json.gz')
    assert alternate.provenance['file'] == 'uniprot-human-reviewed.json.gz'


def test_uniprot_groups_preserve_evidence_and_remap_every_endpoint(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['ProteinIdentifier'][0]['protein_type'] = 'Enzyme'
    reader.data['ProteinIdentifier'] += [node('UniProtKB:OLD1', protein_type='Transporter'), node('UniProtKB:P1-2')]
    reader.data['HmdbMetaboliteProteinAssociationEdge'] += [edge('HMDB:HMDB1', 'UniProtKB:OLD1'),
                                                         edge('HMDB:HMDB1', 'UniProtKB:P1-2')]
    reader.data['GeneIdentifier'] += [node('NCBIGene:9'), node('NCBIGene:UNKNOWN'), node('Symbol:SAME')]
    before = copy.deepcopy(reader.data)
    path = tmp_path / 'resolved.sqlite'
    manifest = run_export(reader, path, gene_identity=identity, progress=lambda _: None)
    assert reader.data == before
    with sqlite3.connect(path) as db:
        ids = dict(db.execute('SELECT sourceId,rampId FROM source'))
        common = ids['uniprot:P1']
        assert {ids[s] for s in ['entrez:1', 'uniprot:OLD1', 'uniprot:P1-2']} == {common}
        assert {ids[s] for s in ['uniprot:P2', 'entrez:9', 'gene_symbol:SAME']} == {common}
        assert ids['entrez:UNKNOWN'] != common
        assert db.execute("SELECT count(*) FROM analyte WHERE type='gene'").fetchone() == (2,)
        assert db.execute('SELECT rampGeneId FROM catalyzed').fetchone() == (common,)
        assert db.execute('SELECT proteinType FROM catalyzed').fetchone() == ('Enzyme; Transporter; Unknown',)
        assert db.execute('SELECT ramp_gene_id,uniprot FROM reaction2protein').fetchone() == (common, 'uniprot:P1')
        assert db.execute('SELECT count(*) FROM analytehaspathway WHERE rampId=?', (common,)).fetchone() == (1,)
        assert db.execute('PRAGMA foreign_key_check').fetchall() == []
        assert db.execute("SELECT dataSource FROM source WHERE sourceId='entrez:1'").fetchone() == ('hmdb',)
        assert db.execute("SELECT sourceId FROM source WHERE dataSource='uniprot' ORDER BY sourceId").fetchall() == [('uniprot:P1',), ('uniprot:P2',)]
    assert manifest['gene_resolution']['counts'] == {
        'GeneIdentifier.resolved': 1, 'GeneIdentifier.ambiguous': 2, 'GeneIdentifier.unmatched': 1,
        'ProteinIdentifier.resolved': 4, 'output_groups': 2,
        'multi_accession_groups': 1, 'largest_canonical_accession_group': 2}
    assert manifest['gene_resolution']['policy_version'] == 2
    for documents in reader.data.values():
        documents.reverse()
    run_export(reader, tmp_path / 'shuffled.sqlite', gene_identity=identity, progress=lambda _: None)
    with sqlite3.connect(path) as a, sqlite3.connect(tmp_path / 'shuffled.sqlite') as b:
        for table in ('source', 'analyte', 'analytehaspathway', 'reaction2protein', 'catalyzed'):
            assert sorted(a.execute(f'SELECT * FROM {table}').fetchall(), key=repr) == sorted(b.execute(f'SELECT * FROM {table}').fetchall(), key=repr)


def test_unused_shared_alias_does_not_merge_proteins(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    # P1/P2 share symbol SAME and GeneID 9 in UniProt, but neither is an input ID.
    path = tmp_path / 'unused.sqlite'
    run_export(reader, path, gene_identity=identity, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        ids = dict(db.execute('SELECT sourceId,rampId FROM source'))
        assert ids['entrez:1'] == ids['uniprot:P1']
        assert ids['uniprot:P1'] != ids['uniprot:P2']
        assert 'gene_symbol:SAME' not in ids


def test_used_alias_collapses_gene_only_matches_and_both_protein_relationships(tmp_path):
    reader, identity = resolved_fixture(tmp_path)
    reader.data['GeneIdentifier'] = [node('NCBIGene:9')]
    reader.data['GenePathwayEdge'] = [edge('NCBIGene:9', 'SMPDB:PW1', details=[{'source': 'HMDB'}])]
    reader.data['RheaProteinReactionEdge'].append(edge('UniProtKB:P2', 'RHEA:1'))
    path = tmp_path / 'collapsed.sqlite'
    run_export(reader, path, gene_identity=identity, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT count(*) FROM analyte WHERE type='gene'").fetchone() == (1,)
        rows = db.execute('SELECT ramp_gene_id,uniprot FROM reaction2protein').fetchall()
        assert len({r[0] for r in rows}) == 1
        assert {r[1] for r in rows} == {'uniprot:P1', 'uniprot:P2'}
    # Canonical IDs without their own graph nodes are also retained by projection.
    reader.data['ProteinIdentifier'] = []
    for collection in ('ProteinPathwayEdge', 'RheaProteinReactionEdge', 'HmdbMetaboliteProteinAssociationEdge'):
        reader.data[collection] = []
    path = tmp_path / 'gene-only.sqlite'
    run_export(reader, path, gene_identity=identity, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT sourceId FROM source WHERE dataSource='uniprot' ORDER BY sourceId").fetchall() == [('uniprot:P1',), ('uniprot:P2',)]


@pytest.mark.parametrize('mismatch', ['snapshot', 'checksum', 'missing'])
def test_gene_resolver_must_match_stage_input_before_publication(tmp_path, mismatch):
    reader, identity = resolved_fixture(tmp_path)
    if mismatch == 'snapshot':
        identity.provenance['snapshot_id'] = 'uniprot:human:other'
    elif mismatch == 'checksum':
        identity.provenance['sha256'] = 'other'
    else:
        reader.metadata['registry_datasets'] = reader.metadata['registry_datasets'][:1]
    path = tmp_path / 'invalid.sqlite'
    with pytest.raises(ValueError):
        run_export(reader, path, gene_identity=identity)
    assert not path.exists()


def test_missing_pinned_resolver_file_fails(tmp_path):
    reader = FixtureReader()
    def missing(name):
        raise FileNotFoundError(name)
    dataset = SimpleNamespace(snapshot_id='uniprot:human:test', file=missing)
    with pytest.raises(FileNotFoundError):
        GeneIdentity.from_stage(reader.metadata, SimpleNamespace(resolve=lambda _: dataset))


@pytest.mark.parametrize("failure", ["missing_member", "missing_endpoint", "changed_graph", "conflicting_class"])
def test_failure_leaves_no_final_artifact(tmp_path, failure):
    reader = FixtureReader()
    if failure == "missing_member":
        reader.data["MetaboliteIdentifier"] = [d for d in reader.data["MetaboliteIdentifier"] if d["id"] != "HMDB:HMDB2"]
    elif failure == "missing_endpoint":
        reader.data["GenePathwayEdge"].append(edge("NCBIGene:MISSING", "SMPDB:PW1"))
    elif failure == "changed_graph":
        def changed():
            raise ValueError("Graph changed")
        reader.verify_unchanged = changed
    else:
        reader.data["MetaboliteClassificationTerm"].append(node("CLASS:other", name="Other", level_name="class", source="HMDB"))
        reader.data["MetaboliteClassificationEdge"].append(edge("HMDB:HMDB1", "CLASS:other"))
    path = tmp_path / "failed.sqlite"
    with pytest.raises((ValueError, sqlite3.IntegrityError)):
        export_sqlite(reader, path, progress=lambda _: None)
    assert not path.exists()
    assert list(tmp_path.iterdir()) == []


def test_existing_file_is_preserved_before_credentials_are_loaded(tmp_path):
    path = tmp_path / "existing.sqlite"
    path.write_bytes(b"keep")
    with pytest.raises(SystemExit):
        main(["--stage-id", "test", "--output", str(path), "--graph-credentials", "missing"])
    assert path.read_bytes() == b"keep"


@pytest.mark.parametrize('existing', [False, True])
def test_overwrite_publishes_validated_database(tmp_path, existing):
    path = tmp_path / 'replace.sqlite'
    if existing:
        path.write_bytes(b'previous export')
    reader = FixtureReader()
    def check_before_publication():
        assert path.read_bytes() == b'previous export' if existing else not path.exists()
    reader.verify_unchanged = check_before_publication
    export_sqlite(reader, path, overwrite=True, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute('PRAGMA integrity_check').fetchone() == ('ok',)
    assert list(tmp_path.iterdir()) == [path]


@pytest.mark.parametrize('failure', ['projection', 'integrity', 'verification', 'publication'])
def test_failed_overwrite_preserves_old_file(tmp_path, monkeypatch, failure):
    from src.use_cases.ramp import sqlite_writer
    reader = FixtureReader()
    path = tmp_path / 'keep.sqlite'
    path.write_bytes(b'previous export')
    def fail(*args):
        raise OSError('simulated failure')
    if failure == 'projection':
        reader.data['GenePathwayEdge'].append(edge('NCBIGene:MISSING', 'SMPDB:PW1'))
    elif failure == 'integrity':
        reader.data['MetaboliteClassificationTerm'].append(node('CLASS:other', name='Other', level_name='class', source='HMDB'))
        reader.data['MetaboliteClassificationEdge'].append(edge('HMDB:HMDB1', 'CLASS:other'))
    elif failure == 'verification':
        reader.verify_unchanged = fail
    else:
        monkeypatch.setattr(sqlite_writer.os, 'replace', fail)
    with pytest.raises((ValueError, OSError, sqlite3.IntegrityError)):
        export_sqlite(reader, path, overwrite=True, progress=lambda _: None)
    assert path.read_bytes() == b'previous export'
    assert list(tmp_path.iterdir()) == [path]


def test_overwrite_rejects_symlink_and_sidecars(tmp_path):
    from src.use_cases.ramp.sqlite_writer import SQLiteWriter
    target = tmp_path / 'target.sqlite'
    target.write_bytes(b'keep')
    link = tmp_path / 'link.sqlite'
    link.symlink_to(target)
    with pytest.raises(ValueError, match='ordinary file'):
        SQLiteWriter(link, overwrite=True)
    journal = tmp_path / 'target.sqlite-wal'
    journal.write_bytes(b'active')
    with pytest.raises(ValueError, match='journal/WAL'):
        SQLiteWriter(target, overwrite=True)
    assert target.read_bytes() == b'keep'


def test_stage_reader_reconstructs_singletons_and_validates_membership():
    reader = StageReader.__new__(StageReader)
    reader.stage_id = "selected"
    reader.stage = {"summary": {"active_identifier_count": 3, "non_singleton_clique_count": 1,
                                "singleton_identifier_count": 1, "clique_count": 2}}

    def execute(query, **kwargs):
        assert kwargs["bind_vars"] == {"s": "selected"}
        if "ActiveIdentifierChunk" in query:
            return [["HMDB:2", "CHEBI:1", "HMDB:1"]]
        if "MemberEdge" in query:
            return [{"harmonized_metabolite_id": "g1", "member_id": i} for i in ("HMDB:1", "CHEBI:1")]
        return [{"id": "g1", "size": 2}]
    reader.db = SimpleNamespace(aql=SimpleNamespace(execute=execute))
    assert reader.groups() == [("CHEBI:1", "HMDB:1"), ("HMDB:2",)]
    reader.stage["summary"]["singleton_identifier_count"] = 0
    with pytest.raises(ValueError, match="totals"):
        reader.groups()


def test_stage_reader_rejects_rebuilt_metadata_even_with_matching_collection_revisions():
    original = {"_key": "old-build", "registry_datasets": [{"snapshot_id": "hmdb:metabolites_xml:5.0"}]}
    stage = {"status": "complete", "summary": {"source_collection_revisions": {"MetaboliteIdentifier": "rev"},
             "source_graph_fingerprint": source_build_fingerprint(original)}}

    class Collection:
        def revision(self):
            return "rev"

        def get(self, key):
            if key == "selected":
                return stage
            if key == "etl_metadata":
                return {"value": {**original, "_key": "new-build"}}
            raise AssertionError(f"Unexpected read: {key}")

    db = SimpleNamespace(collection=lambda _: Collection())
    with pytest.raises(ValueError, match="build metadata differs"):
        StageReader(db, "selected", None)


def test_corrected_chemistry_is_exported_without_mutating_source(tmp_path):
    reader = FixtureReader()
    stage_reader = StageReader.__new__(StageReader)
    path = ["chem_props", {"match": {"source": "HMDB", "source_id": "HMDB:HMDB1"}}, "mw"]
    decision = SimpleNamespace(
        path=path, mode="set", value="20.0", observed_value="18.015", observed_exists=True,
        target={"model_type": "MetaboliteIdentifier", "id": "HMDB:HMDB1"},
        curation_type=METABOLITE_RECORD_PROPERTIES, batch_id="recorded", published_at="2026-10-01", published_by=None,
        source_operation={"operation_id": "correction", "note": "Reviewed mass"},
    )
    stage_reader.decisions = {"MetaboliteIdentifier": {"HMDB:HMDB1": [decision]}}
    stage_reader.schemas = {"MetaboliteIdentifier": {"fields": {
        "chem_props": {"type": "list", "item_type": "object", "fields": {"mw": "str"}},
    }}}
    original_records = reader.records

    def effective_records(collection):
        for doc in original_records(collection):
            yield stage_reader.effective(collection, doc) if collection == "MetaboliteIdentifier" else doc
    reader.records = effective_records
    path = tmp_path / "corrected.sqlite"
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("SELECT mw FROM chem_props").fetchone() == (20.0,)
    assert reader.data["MetaboliteIdentifier"][0]["chem_props"][0]["mw"] == "18.015"


def test_replay_reads_historical_batches_without_latest_manifest():
    kind = METABOLITE_RECORD_PROPERTIES
    operation = {"action": "set_properties", "target": {"kind": "node", "curation_set": "metabolite_harmonization",
                 "model_type": "MetaboliteIdentifier", "id": "HMDB:1"},
                 "decisions": [{"path": ["is_generic_structure"], "mode": "set", "value": True,
                                "observed_exists": False, "observed_value": None}], "note": "reviewed"}
    batch = {"format_version": 2, "curation_type": kind, "curation_batch_id": "old", "operations": [operation]}
    manifest = {"format_version": 2, "curation_type": kind, "revision": 1,
                "batches": [{"batch_id": "old", "sha256": payload_sha256(batch)}]}
    objects = {manifest_key(kind): json.dumps(manifest), batch_key(kind, "old"): json.dumps(batch)}
    storage = SimpleNamespace(bucket="test", read_text=objects.__getitem__)
    original = resolve_curation_type(storage, kind)
    # Latest state no longer contains the old decision; replay must not consult it.
    objects[manifest_key(kind)] = json.dumps({**manifest, "revision": 2, "batches": []})
    assert resolve_curation_type(storage, kind).active_record_property_decisions == []
    del objects[manifest_key(kind)]
    replayed = replay_curation_snapshot(storage, original.metadata())
    assert replayed.active_record_property_decisions[0].value is True
    assert replayed.metadata() == original.metadata()
    objects[batch_key(kind, "old")] = json.dumps({**batch, "operations": []})
    with pytest.raises(ValueError, match="hash mismatch"):
        replay_curation_snapshot(storage, original.metadata())
