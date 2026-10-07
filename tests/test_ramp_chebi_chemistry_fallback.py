import sqlite3
import pytest

from tests.test_ramp_sqlite_export import FixtureReader, export_sqlite


def test_linked_chemical_entity_supplies_missing_identifier_chemistry(tmp_path):
    reader = FixtureReader()
    reader.data['ChebiChemicalEntityMetaboliteIdentifierEdge'] = [
        {'start_id': 'CHEBI:1', 'end_id': 'CHEBI:1', 'source_label': 'Oxidane'},
        {'start_id': 'CHEBI:INACTIVE', 'end_id': 'CHEBI:INACTIVE'},
    ]
    reader.data['ChemicalEntity'] = [
        {'id': 'CHEBI:1', 'name': 'Oxidane', 'smiles': 'O', 'inchi_key': 'KEY-SECOND-PART',
         'inchi': 'InChI=1S/H2O/h1H2', 'formula': 'H2O', 'mass': '18.015', 'monoisotopic_mass': '18.01056'},
        {'id': 'CHEBI:INACTIVE', 'name': 'Inactive', 'formula': 'C2H2'},
        {'id': 'CHEBI:UNLINKED', 'name': 'Unlinked', 'formula': 'C3H3'},
    ]
    path = tmp_path / 'ramp.sqlite'
    manifest = export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select chem_data_source,chem_source_id,iso_smiles,inchi_key_prefix,inchi_key,inchi,mw,monoisotop_mass,common_name,mol_formula from chem_props where chem_source_id='chebi:1'").fetchall() == [
            ('chebi', 'chebi:1', 'O', 'KEY', 'KEY-SECOND-PART', 'InChI=1S/H2O/h1H2', 18.015, 18.01056, 'Oxidane', 'H2O')]
        assert db.execute("select count(*) from chem_props where chem_source_id='hmdb:HMDB1'").fetchone() == (1,)
        assert db.execute("select count(*) from chem_props where chem_source_id in ('chebi:INACTIVE','chebi:UNLINKED')").fetchone() == (0,)
        assert db.execute("select commonName from source where sourceId='chebi:1' and dataSource='chebi'").fetchone() == ('Oxidane',)
        assert db.execute("select count(*) from analyte where type='compound'").fetchone() == (2,)
    assert manifest['chebi_chemistry_fallback'] == {
        'policy': 'Use same-CURIE ChemicalEntity fields only when active MetaboliteIdentifier.chem_props is empty',
        'source_lookup_policy': 'Register a ChEBI bridge as ChEBI only when it supplies retained fallback chemistry',
        'active_bridge_identifiers': 1,
        'added_rows': 1, 'groups_gaining_first_chemistry': 0}


def test_existing_identifier_chemistry_takes_precedence(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][1]['chem_props'] = [{
        'source': 'ChEBI', 'source_id': 'CHEBI:1', 'molecular_formula': 'H2O'}]
    reader.data['ChebiChemicalEntityMetaboliteIdentifierEdge'] = [{'start_id': 'CHEBI:1', 'end_id': 'CHEBI:1'}]
    reader.data['ChemicalEntity'] = [{'id': 'CHEBI:1', 'formula': 'SHOULD_NOT_EXPORT'}]
    path = tmp_path / 'existing.sqlite'
    manifest = export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select mol_formula from chem_props where chem_source_id='chebi:1'").fetchall() == [('H2O',)]
        assert db.execute("select count(*) from source where sourceId='chebi:1' and dataSource='chebi'").fetchone() == (1,)
    assert manifest['chebi_chemistry_fallback']['added_rows'] == 0
    assert manifest['chebi_chemistry_fallback']['active_bridge_identifiers'] == 1


def test_retained_chemistry_registers_reporting_source_without_bridge(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][1]['chem_props'] = [{
        'source': 'ChEBI', 'source_id': 'CHEBI:1', 'molecular_formula': 'H2O',
        'common_name': 'Oxidane'}]
    reader.data['MetaboliteIdentifier'][1]['names'] = [{
        'value': 'Oxidane', 'source': 'ChEBI', 'source_field': 'ChEBI NAME'}]
    reader.data['MetaboliteIdentifier'][1]['sources'].append('ChEBI\t255\t2026-09-09\t2026-09-09')
    path = tmp_path / 'chemistry-source.sqlite'
    export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select count(*) from source where sourceId='chebi:1' and dataSource='chebi'").fetchone() == (1,)
        assert db.execute("select commonName from source where sourceId='chebi:1' and dataSource='chebi'").fetchone() == ('Oxidane',)
        assert db.execute("select Synonym,source from analytesynonym where source='chebi'").fetchall() == [('Oxidane', 'chebi')]


def test_bridge_without_retained_chemistry_does_not_attribute_chebi(tmp_path):
    reader = FixtureReader()
    reader.data['ChebiChemicalEntityMetaboliteIdentifierEdge'] = [
        {'start_id': 'CHEBI:1', 'end_id': 'CHEBI:1', 'source_label': 'Oxidane'}]
    reader.data['ChemicalEntity'] = [{'id': 'CHEBI:1', 'name': 'Oxidane'}]
    path = tmp_path / 'no-chemistry.sqlite'
    manifest = export_sqlite(reader, path, progress=lambda _: None)
    with sqlite3.connect(path) as db:
        assert db.execute("select commonName from source where sourceId='chebi:1' and dataSource='chebi'").fetchall() == []
        assert db.execute("select count(*) from chem_props where chem_source_id='chebi:1'").fetchone() == (0,)
    assert manifest['chebi_chemistry_fallback']['active_bridge_identifiers'] == 1
    assert manifest['chebi_chemistry_fallback']['added_rows'] == 0


def test_fallback_counts_first_chemistry_for_group(tmp_path):
    reader = FixtureReader()
    reader.data['MetaboliteIdentifier'][0]['chem_props'] = []
    reader.data['ChebiChemicalEntityMetaboliteIdentifierEdge'] = [
        {'start_id': 'CHEBI:1', 'end_id': 'CHEBI:1'}]
    reader.data['ChemicalEntity'] = [{'id': 'CHEBI:1', 'formula': 'H2O'}]
    manifest = export_sqlite(reader, tmp_path / 'first.sqlite', progress=lambda _: None)
    assert manifest['chebi_chemistry_fallback']['added_rows'] == 1
    assert manifest['chebi_chemistry_fallback']['groups_gaining_first_chemistry'] == 1


def test_fallback_rejects_link_to_another_chemical_identity(tmp_path):
    reader = FixtureReader()
    reader.data['ChebiChemicalEntityMetaboliteIdentifierEdge'] = [
        {'start_id': 'CHEBI:2', 'end_id': 'CHEBI:1'}]
    with pytest.raises(ValueError, match='changes identifier'):
        export_sqlite(reader, tmp_path / 'bad.sqlite', progress=lambda _: None)
