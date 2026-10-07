from pathlib import Path
import gzip
import json
import zipfile
import rdkit

from src.input_adapters.metabolite_harmonization.chemprops import (
    ChebiMetaboliteChemPropsAdapter,
    HmdbMetaboliteChemPropsAdapter,
    LipidMapsMetaboliteChemPropsAdapter,
    PubchemMetaboliteChemPropsAdapter,
)
from src.interfaces.output_adapter import OutputAdapter
from src.models.metabolite_harmonization import MetaboliteIdentifier
from src.shared.metabolite_generic_structure import classify_generic_structure
from src.shared.metabolite_structure_chemistry import (
    STRUCTURE_CALCULATION_METHOD,
    calculate_metabolite_chem_props_derivatives,
    calculate_smiles_chemistry,
)
from src.shared.record_merger import FieldConflictBehavior


class _ConvertingOutputAdapter(OutputAdapter):
    def store(self, objects, single_source=False, field_conflict_behavior=FieldConflictBehavior.KeepFirst) -> bool:
        return True

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        return True


def _sdf_record(title: str, tags: dict) -> str:
    lines = [
        title,
        "  TEST",
        "",
        "  0  0  0     0  0            999 V2000",
        "M  END",
        "",
    ]
    for key, value in tags.items():
        lines.extend([f"> <{key}>", value, ""])
    lines.append("$$$$")
    return "\n".join(lines) + "\n"


def _sdf_record_with_mol_block(title: str, atom_lines: list[str], bond_lines: list[str], tags: dict) -> str:
    lines = [
        title,
        "  TEST",
        "",
        f"{len(atom_lines):>3}{len(bond_lines):>3}  0     0  0            999 V2000",
        *atom_lines,
        *bond_lines,
        "M  END",
        "",
    ]
    for key, value in tags.items():
        lines.extend([f"> <{key}>", value, ""])
    lines.append("$$$$")
    return "\n".join(lines) + "\n"


def _records(adapter):
    return [record for batch in adapter.get_all() for record in batch]


def test_hmdb_chemprops_adapter_emits_chemprops(tmp_path: Path):
    zip_path = tmp_path / "structures.zip"
    sdf = _sdf_record(
        "HMDB0000001",
        {
            "DATABASE_ID": "HMDB0000001",
            "SMILES": "C[NH+]1C=NC=C1CC(N)C(=O)O",
            "INCHI_KEY": "BRMWTNUJHUMWMS-LURJTMIESA-N",
            "INCHI_IDENTIFIER": "InChI=1S/C7H11N3O2",
            "MOLECULAR_WEIGHT": "169.1811",
            "EXACT_MASS": "169.085126611",
            "GENERIC_NAME": "1-Methylhistidine",
            "FORMULA": "C7H11N3O2",
        },
    )
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    node = _records(HmdbMetaboliteChemPropsAdapter(structures_zip_file=str(zip_path)))[0]

    assert node.id == "HMDB:HMDB0000001"
    assert node.prefix == "HMDB"
    assert node.chem_props[0].source == "HMDB"
    assert node.chem_props[0].source_id == "HMDB:HMDB0000001"
    assert node.chem_props[0].iso_smiles == "C[NH+]1C=NC=C1CC(N)C(=O)O"
    assert node.chem_props[0].inchi_key_prefix == "BRMWTNUJHUMWMS"
    assert node.chem_props[0].derived_inchi_key == "JDHILDINMRGULE-UHFFFAOYSA-O"
    assert node.chem_props[0].derived_inchi_key_prefix == "JDHILDINMRGULE"
    assert node.chem_props[0].derived_inchi_key_input_field == "iso_smiles"
    assert node.chem_props[0].derived_inchi_key_method == "rdkit.Chem.inchi.MolToInchiKey"
    assert node.chem_props[0].derived_inchi_key_method_version == rdkit.__version__
    assert node.chem_props[0].derived_inchi_key_error is None
    assert node.chem_props[0].monoisotopic_mass == "169.085126611"
    assert node.chem_props[0].common_name == "1-Methylhistidine"
    assert node.chem_props[0].molecular_formula == "C7H11N3O2"


def test_hmdb_chemprops_adapter_uses_xml_only_for_accessions_absent_from_sdf(
    tmp_path: Path,
):
    structures_zip = tmp_path / "structures.zip"
    sdf = _sdf_record(
        "HMDB0000001",
        {
            "DATABASE_ID": "HMDB0000001",
            "SMILES": "C",
            "INCHI_KEY": "VNWKTOKETHGBQD-UHFFFAOYSA-N",
            "INCHI_IDENTIFIER": "InChI=1S/CH4/h1H4",
            "MOLECULAR_WEIGHT": "16.043",
            "EXACT_MASS": "16.031300128",
            "GENERIC_NAME": "Methane from SDF",
            "FORMULA": "CH4",
        },
    )
    with zipfile.ZipFile(structures_zip, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    metabolites_zip = tmp_path / "hmdb_metabolites.zip"
    xml = """<hmdb xmlns="http://www.hmdb.ca">
<metabolite>
  <accession>HMDB0000001</accession>
  <name>Conflicting XML value</name>
  <chemical_formula>C2H6O</chemical_formula>
  <average_molecular_weight>46.07</average_molecular_weight>
  <monisotopic_molecular_weight>46.041864812</monisotopic_molecular_weight>
  <smiles>CCO</smiles>
  <inchi>InChI=1S/C2H6O</inchi>
  <inchikey>LFQSCWFLJHTTHZ-UHFFFAOYSA-N</inchikey>
</metabolite>
<metabolite>
  <accession>HMDB0000002</accession>
  <name>Ethanol from XML fallback</name>
  <chemical_formula>C2H6O</chemical_formula>
  <average_molecular_weight>46.07</average_molecular_weight>
  <monisotopic_molecular_weight>46.041864812</monisotopic_molecular_weight>
  <smiles>CCO</smiles>
  <inchi>InChI=1S/C2H6O</inchi>
  <inchikey>LFQSCWFLJHTTHZ-UHFFFAOYSA-N</inchikey>
</metabolite>
<metabolite>
  <accession>HMDB0000003</accession>
  <name>No source chemistry</name>
</metabolite>
<metabolite>
  <accession>HMDB0000002</accession>
  <name>Duplicate XML value that must be ignored</name>
  <chemical_formula>C3H8O</chemical_formula>
</metabolite>
</hmdb>
"""
    with zipfile.ZipFile(metabolites_zip, "w") as archive:
        archive.writestr("hmdb_metabolites.xml", xml)

    records = _records(
        HmdbMetaboliteChemPropsAdapter(
            structures_zip_file=str(structures_zip),
            metabolites_zip_file=str(metabolites_zip),
        )
    )
    nodes = {record.id: record for record in records}

    assert len(records) == 2
    assert set(nodes) == {"HMDB:HMDB0000001", "HMDB:HMDB0000002"}
    sdf_props = nodes["HMDB:HMDB0000001"].chem_props[0]
    assert sdf_props.common_name == "Methane from SDF"
    assert sdf_props.molecular_formula == "CH4"
    assert sdf_props.inchi_key == "VNWKTOKETHGBQD-UHFFFAOYSA-N"

    xml_props = nodes["HMDB:HMDB0000002"].chem_props[0]
    assert xml_props.source == "HMDB"
    assert xml_props.source_id == "HMDB:HMDB0000002"
    assert xml_props.common_name == "Ethanol from XML fallback"
    assert xml_props.molecular_formula == "C2H6O"
    assert xml_props.mw == "46.07"
    assert xml_props.monoisotopic_mass == "46.041864812"
    assert xml_props.iso_smiles == "CCO"
    assert xml_props.inchi == "InChI=1S/C2H6O"
    assert xml_props.inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert xml_props.derived_inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"


def test_chebi_chemprops_adapter_reads_gzipped_sdf(tmp_path: Path):
    gz_path = tmp_path / "chebi_3_stars.sdf.gz"
    sdf = _sdf_record(
        "CHEBI:27596",
        {
            "ChEBI ID": "CHEBI:27596",
            "SMILES": "Cn1cncc1C[C@H](N)C(=O)O",
            "INCHIKEY": "BRMWTNUJHUMWMS-LURJTMIESA-N",
            "INCHI": "InChI=1S/C7H11N3O2",
            "MASS": "169.184",
            "MONOISOTOPIC_MASS": "169.085126611",
            "ChEBI NAME": "1-methyl-L-histidine",
            "FORMULA": "C7H11N3O2",
        },
    )
    with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
        handle.write(sdf)

    node = _records(ChebiMetaboliteChemPropsAdapter(chebi_sdf_file=str(gz_path)))[0]

    assert node.id == "CHEBI:27596"
    assert [name.to_dict() for name in node.names] == [
        {"value": "1-methyl-L-histidine", "source": "ChEBI", "source_field": "ChEBI NAME"}
    ]
    assert node.chem_props[0].source == "ChEBI"
    assert node.chem_props[0].mw == "169.184"
    assert node.chem_props[0].monoisotopic_mass == "169.085126611"
    assert node.chem_props[0].common_name == "1-methyl-L-histidine"
    assert node.chem_props[0].inchi_key == "BRMWTNUJHUMWMS-LURJTMIESA-N"
    assert node.chem_props[0].inchi == "InChI=1S/C7H11N3O2"
    assert node.chem_props[0].molecular_formula == "C7H11N3O2"
    assert node.chem_props[0].derived_inchi_key == "JDHILDINMRGULE-LURJTMIESA-N"


def test_chebi_chemprops_adapter_omits_partial_generic_structure_masses(tmp_path: Path):
    gz_path = tmp_path / "chebi_3_stars.sdf.gz"
    sdf = _sdf_record(
        "CHEBI:21494",
        {
            "ChEBI ID": "CHEBI:21494",
            "SMILES": "*NOC(C)=O",
            "MASS": "74.059",
            "MONOISOTOPIC_MASS": "74.0242",
            "ChEBI NAME": "N-acetoxyarylamine",
            "FORMULA": "C2H4NO2R",
        },
    )
    with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
        handle.write(sdf)

    props = _records(
        ChebiMetaboliteChemPropsAdapter(chebi_sdf_file=str(gz_path))
    )[0].chem_props[0]

    assert props.mw is None
    assert props.monoisotopic_mass is None
    assert props.molecular_formula == "C2H4NO2R"
    assert props.iso_smiles == "*NOC(C)=O"


def test_chebi_chemprops_adapter_preserves_legacy_name_tag(tmp_path: Path):
    gz_path = tmp_path / "legacy_chebi.sdf.gz"
    with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
        handle.write(_sdf_record("CHEBI:1", {
            "ChEBI ID": "CHEBI:1", "ChEBI Name": "Water"}))

    node = _records(ChebiMetaboliteChemPropsAdapter(chebi_sdf_file=str(gz_path)))[0]
    assert [name.to_dict() for name in node.names] == [
        {"value": "Water", "source": "ChEBI", "source_field": "ChEBI Name"}]
    assert node.chem_props[0].common_name == "Water"


def test_lipidmaps_chemprops_adapter_emits_chemprops(tmp_path: Path):
    zip_path = tmp_path / "LMSD.sdf.zip"
    sdf = _sdf_record(
        "LMFA00000001",
        {
            "LM_ID": "LMFA00000001",
            "SMILES": "CCC(=O)C(=O)O",
            "INCHI_KEY": "TYEYBOSBBBHJIV-UHFFFAOYSA-N",
            "INCHI": "InChI=1S/C4H6O3",
            "MASS": "102.09",
            "EXACT_MASS": "102.031695",
            "COMMON_NAME": "3-methyl pyruvic acid",
            "FORMULA": "C4H6O3",
        },
    )
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    node = _records(LipidMapsMetaboliteChemPropsAdapter(sdf_zip_file=str(zip_path)))[0]

    assert node.id == "LIPIDMAPS:LMFA00000001"
    assert node.chem_props[0].source == "LipidMaps"
    assert node.chem_props[0].mw == "102.09"
    assert node.chem_props[0].monoisotopic_mass == "102.031695"
    assert node.chem_props[0].common_name == "3-methyl pyruvic acid"
    assert node.chem_props[0].derived_inchi_key == "TYEYBOSBBBHJIV-UHFFFAOYSA-N"


def test_lipidmaps_chemprops_derives_wildcard_smiles_from_generic_mol_block(tmp_path: Path):
    zip_path = tmp_path / "LMSD.sdf.zip"
    sdf = _sdf_record_with_mol_block(
        "LMGP04050000",
        [
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    1.0000    0.0000    0.0000 R   0  0  0  0  0  0  0  0  0  0  0  0",
        ],
        ["  1  2  1  0  0  0  0"],
        {"LM_ID": "LMGP04050000", "NAME": "LPG"},
    )
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    node = _records(LipidMapsMetaboliteChemPropsAdapter(sdf_zip_file=str(zip_path)))[0]
    props = node.chem_props[0]

    assert "*" in props.canonical_smiles
    assert props.iso_smiles is None
    assert props.derived_inchi_key is None
    assert props.derived_inchi_key_input_field == "sdf_mol_block"
    assert props.derived_inchi_key_error == "query_structure_not_supported"
    assert classify_generic_structure([{"smiles": props.canonical_smiles}]) is True


def test_lipidmaps_chemprops_derives_smiles_and_inchi_key_from_concrete_mol_block(tmp_path: Path):
    zip_path = tmp_path / "LMSD.sdf.zip"
    sdf = _sdf_record_with_mol_block(
        "LMFA00000003",
        [
            "    0.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    1.0000    0.0000    0.0000 C   0  0  0  0  0  0  0  0  0  0  0  0",
            "    2.0000    0.0000    0.0000 O   0  0  0  0  0  0  0  0  0  0  0  0",
        ],
        [
            "  1  2  1  0  0  0  0",
            "  2  3  1  0  0  0  0",
        ],
        {"LM_ID": "LMFA00000003", "NAME": "ethanol"},
    )
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    node = _records(LipidMapsMetaboliteChemPropsAdapter(sdf_zip_file=str(zip_path)))[0]
    props = node.chem_props[0]

    assert props.canonical_smiles == "CCO"
    assert props.derived_inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert props.derived_inchi_key_prefix == "LFQSCWFLJHTTHZ"
    assert props.derived_inchi_key_input_field == "sdf_mol_block"
    assert props.derived_inchi_key_error is None


def test_lipidmaps_tagged_smiles_takes_precedence_over_mol_block(tmp_path: Path):
    zip_path = tmp_path / "LMSD.sdf.zip"
    sdf = _sdf_record_with_mol_block(
        "LMFA00000004",
        ["    0.0000    0.0000    0.0000 R   0  0  0  0  0  0  0  0  0  0  0  0"],
        [],
        {"LM_ID": "LMFA00000004", "SMILES": "CCO"},
    )
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    props = _records(LipidMapsMetaboliteChemPropsAdapter(sdf_zip_file=str(zip_path)))[0].chem_props[0]

    assert props.iso_smiles == "CCO"
    assert props.canonical_smiles is None
    assert props.derived_inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert props.derived_inchi_key_input_field == "iso_smiles"


def test_lipidmaps_zero_atom_mol_block_remains_unknown_without_derivation_error(tmp_path: Path):
    zip_path = tmp_path / "LMSD.sdf.zip"
    sdf = _sdf_record("LMSP0502AB00", {"LM_ID": "LMSP0502AB00"})
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    props = _records(LipidMapsMetaboliteChemPropsAdapter(sdf_zip_file=str(zip_path)))[0].chem_props[0]

    assert props.canonical_smiles is None
    assert props.derived_inchi_key is None
    assert props.derived_inchi_key_input_field is None
    assert props.derived_inchi_key_error is None


def test_lipidmaps_malformed_nonempty_mol_block_records_failure(tmp_path: Path):
    zip_path = tmp_path / "LMSD.sdf.zip"
    sdf = "\n".join([
        "LMFA00000005",
        "  TEST",
        "",
        "  1  0  0     0  0            999 V2000",
        "M  END",
        "",
        "> <LM_ID>",
        "LMFA00000005",
        "",
        "$$$$",
        "",
    ])
    with zipfile.ZipFile(zip_path, "w") as archive:
        archive.writestr("structures.sdf", sdf)

    props = _records(LipidMapsMetaboliteChemPropsAdapter(sdf_zip_file=str(zip_path)))[0].chem_props[0]

    assert props.canonical_smiles is None
    assert props.derived_inchi_key is None
    assert props.derived_inchi_key_input_field == "sdf_mol_block"
    assert props.derived_inchi_key_error == "mol_parse_failed"


def test_pubchem_chemprops_adapter_emits_chemprops(tmp_path: Path):
    tsv_path = tmp_path / "cid_molecular_info.tsv"
    tsv_path.write_text(
        "pubchem_id\tcid\tmonoisotopic_mass\tinchikey\tinchi_key_prefix\tmolecular_formula\tmolecular_weight\tcanonical_smiles\tisomeric_smiles\tinchi\tiupac_name\n"
        "PUBCHEM.COMPOUND:92105\t92105\t169.085126602\tBRMWTNUJHUMWMS-LURJTMIESA-N\tBRMWTNUJHUMWMS\tC7H11N3O2\t169.18\tCN1C=NC=C1CC(N)C(=O)O\tCN1C=NC=C1C[C@H](N)C(=O)O\tInChI=1S/C7H11N3O2\t2-amino-3-(1-methylimidazol-4-yl)propanoic acid\n",
        encoding="utf-8",
    )

    node = _records(PubchemMetaboliteChemPropsAdapter(molecular_info_file=str(tsv_path)))[0]

    assert node.id == "PUBCHEM.COMPOUND:92105"
    assert node.chem_props[0].source == "PubChem"
    assert node.chem_props[0].source_id == "PUBCHEM.COMPOUND:92105"
    assert node.chem_props[0].iso_smiles == "CN1C=NC=C1C[C@H](N)C(=O)O"
    assert node.chem_props[0].canonical_smiles == "CN1C=NC=C1CC(N)C(=O)O"
    assert node.chem_props[0].isomeric_smiles == "CN1C=NC=C1C[C@H](N)C(=O)O"
    assert node.chem_props[0].mw == "169.18"
    assert node.chem_props[0].iupac_name == "2-amino-3-(1-methylimidazol-4-yl)propanoic acid"
    assert node.names == []  # The old pinned TSV has no PubChem Title.
    assert node.chem_props[0].derived_inchi_key == "JDHILDINMRGULE-LURJTMIESA-N"
    assert node.chem_props[0].derived_inchi_key_input_field == "iso_smiles"


def test_pubchem_title_is_a_source_owned_name_when_future_snapshot_provides_it(tmp_path: Path):
    tsv_path = tmp_path / "cid_molecular_info.tsv"
    tsv_path.write_text(
        "pubchem_id\tcid\ttitle\tiupac_name\n"
        "PUBCHEM.COMPOUND:1\t1\tAcetyl-DL-carnitine\t3-acetoxy-4-(trimethylammonio)butanoate\n",
        encoding="utf-8",
    )
    node = _records(PubchemMetaboliteChemPropsAdapter(molecular_info_file=str(tsv_path)))[0]
    assert [(name.value, name.source, name.source_field) for name in node.names] == [
        ("Acetyl-DL-carnitine", "PubChem", "title")]
    assert node.chem_props[0].iupac_name == "3-acetoxy-4-(trimethylammonio)butanoate"


def test_chemprops_are_json_serializable_after_output_conversion(tmp_path: Path):
    tsv_path = tmp_path / "cid_molecular_info.tsv"
    tsv_path.write_text(
        "pubchem_id\tcid\tmonoisotopic_mass\tinchikey\tinchi_key_prefix\tmolecular_formula\tmolecular_weight\tcanonical_smiles\tisomeric_smiles\tinchi\tiupac_name\n"
        "PUBCHEM.COMPOUND:92105\t92105\t169.085126602\tBRMWTNUJHUMWMS-LURJTMIESA-N\tBRMWTNUJHUMWMS\tC7H11N3O2\t169.18\tcanonical\tisomeric\tInChI=1S/C7H11N3O2\tname\n",
        encoding="utf-8",
    )
    adapter = PubchemMetaboliteChemPropsAdapter(molecular_info_file=str(tsv_path))
    output = _ConvertingOutputAdapter()

    converted_groups = output.sort_and_convert_objects(_records(adapter))
    converted_records = [
        record
        for group in converted_groups.values()
        for record in group[0]
    ]

    json.dumps(converted_records)
    node = next(record for record in converted_records if record["id"] == "PUBCHEM.COMPOUND:92105")
    assert node["chem_props"] == [
        {
            "source": "PubChem",
            "source_id": "PUBCHEM.COMPOUND:92105",
            "iso_smiles": "isomeric",
            "canonical_smiles": "canonical",
            "isomeric_smiles": "isomeric",
            "inchi_key_prefix": "BRMWTNUJHUMWMS",
            "inchi_key": "BRMWTNUJHUMWMS-LURJTMIESA-N",
            "derived_inchi_key_prefix": None,
            "derived_inchi_key": None,
            "derived_inchi_key_input_field": "iso_smiles",
            "derived_inchi_key_method": "rdkit.Chem.inchi.MolToInchiKey",
            "derived_inchi_key_method_version": rdkit.__version__,
            "derived_inchi_key_error": "smiles_parse_failed",
            "inchi": "InChI=1S/C7H11N3O2",
            "mw": "169.18",
            "monoisotopic_mass": "169.085126602",
            "common_name": None,
            "iupac_name": "name",
            "molecular_formula": "C7H11N3O2",
            "calculated_mw": None,
            "calculated_monoisotopic_mass": None,
            "structure_components": [],
            "structure_calculation_input_field": "iso_smiles",
            "structure_calculation_method": STRUCTURE_CALCULATION_METHOD,
            "structure_calculation_method_version": rdkit.__version__,
            "structure_calculation_error": "smiles_parse_failed",
        }
    ]


def test_structure_chemistry_preserves_repeated_isotope_aware_components():
    chemistry = calculate_smiles_chemistry("[2H]O[2H].[2H]O[2H].[Cl-]", "iso_smiles")

    assert chemistry.get("structure_calculation_error") is None
    assert len(chemistry["structure_components"]) == 3
    water_components = [
        component
        for component in chemistry["structure_components"]
        if "D2O" in component["molecular_formula"]
    ]
    assert len(water_components) == 2
    assert water_components[0]["monoisotopic_mass"] == water_components[1]["monoisotopic_mass"]
    chloride = next(
        component for component in chemistry["structure_components"]
        if component["formal_charge"] == -1
    )
    assert chloride["molecular_formula"] == "Cl-"


def test_metabolite_chem_props_derivatives_replace_all_structure_managed_fields():
    corrected_smiles = (
        "[2H]C([2H])([2H])C(O)(CCC[C@H](CC#CC(O)(C(F)(F)F)C(F)(F)F)"
        "[C@H]1CC[C@H]2/C(=C/C=C3/C[C@@H](O)C[C@H](O)C3=C)CCC[C@]12C)"
        "C([2H])([2H])[2H]"
    )

    derived = calculate_metabolite_chem_props_derivatives({
        "iso_smiles": corrected_smiles,
        "isomeric_smiles": "CC",
        "canonical_smiles": "C",
    })

    assert derived["structure_calculation_input_field"] == "iso_smiles"
    assert derived["derived_inchi_key_input_field"] == "iso_smiles"
    assert derived["calculated_mw"] == "612.724610668"
    assert derived["calculated_monoisotopic_mass"] == "612.352039684"
    assert derived["derived_inchi_key"] == "VSOWXEHVBCDXAY-CVRMBSQLSA-N"
    assert derived["structure_components"][0]["molecular_formula"] == "C32H38D6F6O4"
    assert derived["structure_calculation_error"] is None
    assert derived["derived_inchi_key_error"] is None


def test_metabolite_chem_props_derivatives_clear_stale_values_without_structure():
    derived = calculate_metabolite_chem_props_derivatives({
        "iso_smiles": None,
        "isomeric_smiles": "",
        "canonical_smiles": None,
    })

    assert derived["calculated_mw"] is None
    assert derived["calculated_monoisotopic_mass"] is None
    assert derived["derived_inchi_key"] is None
    assert derived["structure_components"] == []
    assert derived["structure_calculation_error"] is None


def test_metabolite_chem_props_derivatives_clear_stale_values_for_invalid_structure():
    derived = calculate_metabolite_chem_props_derivatives({
        "iso_smiles": "not-a-smiles",
    })

    assert derived["calculated_mw"] is None
    assert derived["calculated_monoisotopic_mass"] is None
    assert derived["derived_inchi_key"] is None
    assert derived["structure_components"] == []
    assert derived["structure_calculation_error"] == "smiles_parse_failed"
    assert derived["derived_inchi_key_error"] == "smiles_parse_failed"


def test_chemprops_adapters_honor_max_records(tmp_path: Path):
    tsv_path = tmp_path / "cid_molecular_info.tsv"
    tsv_path.write_text(
        "pubchem_id\tcid\tmonoisotopic_mass\tinchikey\tinchi_key_prefix\tmolecular_formula\tmolecular_weight\tcanonical_smiles\tisomeric_smiles\tinchi\tiupac_name\n"
        "PUBCHEM.COMPOUND:1\t1\t1\tA-B-C\tA\tH\t1\tC\tC\tInChI=1\tone\n"
        "PUBCHEM.COMPOUND:2\t2\t2\tD-E-F\tD\tH2\t2\tCC\tCC\tInChI=2\ttwo\n",
        encoding="utf-8",
    )

    records = _records(PubchemMetaboliteChemPropsAdapter(molecular_info_file=str(tsv_path), max_records=1))

    assert [record.id for record in records] == ["PUBCHEM.COMPOUND:1"]


def test_chemprops_without_smiles_do_not_record_a_derivation_attempt(tmp_path: Path):
    gz_path = tmp_path / "chebi_3_stars.sdf.gz"
    sdf = _sdf_record(
        "CHEBI:1",
        {
            "ChEBI ID": "CHEBI:1",
            "InChIKey": "LFQSCWFLJHTTHZ-UHFFFAOYSA-N",
        },
    )
    with gzip.open(gz_path, "wt", encoding="utf-8") as handle:
        handle.write(sdf)

    node = _records(ChebiMetaboliteChemPropsAdapter(chebi_sdf_file=str(gz_path)))[0]

    assert node.chem_props[0].inchi_key == "LFQSCWFLJHTTHZ-UHFFFAOYSA-N"
    assert node.chem_props[0].derived_inchi_key is None
    assert node.chem_props[0].derived_inchi_key_input_field is None
    assert node.chem_props[0].derived_inchi_key_error is None
