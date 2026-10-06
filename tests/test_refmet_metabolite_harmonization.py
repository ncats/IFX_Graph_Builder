from pathlib import Path
import json

from src.input_adapters.metabolite_harmonization.refmet import RefMetMetaboliteEquivalenceAdapter
from src.interfaces.output_adapter import OutputAdapter
from src.models.metabolite_harmonization import MetaboliteIdentifier, MetaboliteIdentifierMappingEdge
from src.shared.metabolite_generic_structure import classify_generic_structure
from src.shared.record_merger import FieldConflictBehavior


class _ConvertingOutputAdapter(OutputAdapter):
    def store(self, objects, single_source=False, field_conflict_behavior=FieldConflictBehavior.KeepFirst) -> bool:
        return True

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        return True


def _write_refmet_csv(path: Path):
    path.write_text(
        " refmet_id,refmet_name,super_class,main_class,sub_class,formula,exactmass,pubchem_cid,chebi_id,hmdb_id,lipidmaps_id,kegg_id,inchi_key\n"
        "RM0108606,Acutumidine,Alkaloids,Alkaloids,Other alkaloids,C18H22ClNO6,383.113567,CID442840,CHEBI:2452,HMDB0015495,LMPK11000003,C10565,SBALNGLYQFMKPR-NQTWQHAWSA-N\n"
        "RM0108637,Adifoline,Alkaloids,Alkaloids,Other alkaloids,C22H20N2O7,424.127053,441972,2488,,,C09020,DJWXVEDJWPDUBQ-DEALGVFLSA-N\n",
        encoding="utf-8",
    )


def _records(adapter):
    return [record for batch in adapter.get_all() for record in batch]


def test_refmet_adapter_emits_nodes_names_and_edges(tmp_path: Path):
    csv_path = tmp_path / "refmet.csv"
    _write_refmet_csv(csv_path)
    adapter = RefMetMetaboliteEquivalenceAdapter(refmet_csv_file=str(csv_path))
    assert adapter.get_datasource_name().value == "RefMet"

    records = _records(adapter)
    nodes = [record for record in records if isinstance(record, MetaboliteIdentifier)]
    edges = [record for record in records if isinstance(record, MetaboliteIdentifierMappingEdge)]
    nodes_by_id = {node.id: node for node in nodes}
    edge_pairs = {(edge.start_node.id, edge.end_node.id) for edge in edges}

    primary = nodes_by_id["REFMET:RM0108606"]
    assert primary.prefix == "REFMET"
    assert primary.names[0].value == "Acutumidine"
    assert primary.names[0].source == "RefMet"
    assert len(primary.chem_props) == 1
    assert primary.chem_props[0].source == "RefMet"
    assert primary.chem_props[0].source_id == "REFMET:RM0108606"
    assert primary.chem_props[0].molecular_formula == "C18H22ClNO6"
    assert primary.chem_props[0].monoisotopic_mass == "383.113567"
    assert primary.chem_props[0].inchi_key == "SBALNGLYQFMKPR-NQTWQHAWSA-N"
    assert primary.chem_props[0].inchi_key_prefix == "SBALNGLYQFMKPR"
    assert classify_generic_structure([{
        "source": primary.chem_props[0].source,
        "source_id": primary.chem_props[0].source_id,
        "formula": primary.chem_props[0].molecular_formula,
    }]) is False
    assert nodes_by_id["PUBCHEM.COMPOUND:442840"].prefix == "PUBCHEM.COMPOUND"
    assert nodes_by_id["PUBCHEM.COMPOUND:442840"].names[0].value == "Acutumidine"
    assert nodes_by_id["PUBCHEM.COMPOUND:442840"].names[0].source == "RefMet"

    expected_node_ids = {
        "REFMET:RM0108606",
        "PUBCHEM.COMPOUND:442840",
        "CHEBI:2452",
        "HMDB:HMDB0015495",
        "LIPIDMAPS:LMPK11000003",
        "KEGG.COMPOUND:C10565",
        "REFMET:RM0108637",
        "PUBCHEM.COMPOUND:441972",
        "CHEBI:2488",
        "KEGG.COMPOUND:C09020",
    }
    assert expected_node_ids <= set(nodes_by_id)
    assert "InChIKey:SBALNGLYQFMKPR-NQTWQHAWSA-N" not in nodes_by_id
    assert "InChIKey:DJWXVEDJWPDUBQ-DEALGVFLSA-N" not in nodes_by_id

    assert ("REFMET:RM0108606", "CHEBI:2452") in edge_pairs
    assert ("REFMET:RM0108606", "PUBCHEM.COMPOUND:442840") in edge_pairs
    assert ("REFMET:RM0108606", "LIPIDMAPS:LMPK11000003") in edge_pairs
    assert ("REFMET:RM0108606", "InChIKey:SBALNGLYQFMKPR-NQTWQHAWSA-N") not in edge_pairs

    details = {
        (edge.start_node.id, edge.end_node.id): edge.details[0]
        for edge in edges
    }
    assert details[("REFMET:RM0108606", "CHEBI:2452")].source == "RefMet"
    assert details[("REFMET:RM0108606", "CHEBI:2452")].source_field == "chebi_id"
    assert details[("REFMET:RM0108606", "CHEBI:2452")].source_id == "REFMET:RM0108606"


def test_refmet_records_are_json_serializable_after_output_conversion(tmp_path: Path):
    csv_path = tmp_path / "refmet.csv"
    _write_refmet_csv(csv_path)
    adapter = RefMetMetaboliteEquivalenceAdapter(refmet_csv_file=str(csv_path), max_records=1)
    output = _ConvertingOutputAdapter()

    converted_groups = output.sort_and_convert_objects(_records(adapter))
    converted_records = [
        record
        for group in converted_groups.values()
        for record in group[0]
    ]

    json.dumps(converted_records)
    primary = next(record for record in converted_records if record["id"] == "REFMET:RM0108606")
    chebi_edge = next(record for record in converted_records if record.get("end_id") == "CHEBI:2452")
    assert primary["prefix"] == "REFMET"
    assert primary["names"] == [
        {"value": "Acutumidine", "source": "RefMet", "source_field": "refmet_name"}
    ]
    assert primary["chem_props"] == [{
        "source": "RefMet",
        "source_id": "REFMET:RM0108606",
        "iso_smiles": None,
        "canonical_smiles": None,
        "isomeric_smiles": None,
        "inchi_key_prefix": "SBALNGLYQFMKPR",
        "inchi_key": "SBALNGLYQFMKPR-NQTWQHAWSA-N",
        "derived_inchi_key_prefix": None,
        "derived_inchi_key": None,
        "derived_inchi_key_input_field": None,
        "derived_inchi_key_method": None,
        "derived_inchi_key_method_version": None,
        "derived_inchi_key_error": None,
        "inchi": None,
        "mw": None,
        "monoisotopic_mass": "383.113567",
        "common_name": None,
        "iupac_name": None,
        "molecular_formula": "C18H22ClNO6",
        "calculated_mw": None,
        "calculated_monoisotopic_mass": None,
        "structure_components": [],
        "structure_calculation_input_field": None,
        "structure_calculation_method": None,
        "structure_calculation_method_version": None,
        "structure_calculation_error": None,
    }]
    assert chebi_edge["details"] == [
        {"source": "RefMet", "source_field": "chebi_id", "source_id": "REFMET:RM0108606"}
    ]


def test_refmet_adapter_honors_max_records(tmp_path: Path):
    csv_path = tmp_path / "refmet.csv"
    _write_refmet_csv(csv_path)
    adapter = RefMetMetaboliteEquivalenceAdapter(refmet_csv_file=str(csv_path), max_records=1)

    records = _records(adapter)
    assert any(isinstance(record, MetaboliteIdentifier) and record.id == "REFMET:RM0108606" for record in records)
    assert not any(isinstance(record, MetaboliteIdentifier) and record.id == "REFMET:RM0108637" for record in records)


def test_refmet_adapter_omits_empty_chem_props(tmp_path: Path):
    csv_path = tmp_path / "refmet.csv"
    csv_path.write_text(
        "refmet_id,refmet_name,formula,exactmass,inchi_key\n"
        "RM0000001,Unknown chemistry,,,\n",
        encoding="utf-8",
    )

    primary = next(
        record
        for record in _records(
            RefMetMetaboliteEquivalenceAdapter(refmet_csv_file=str(csv_path))
        )
        if isinstance(record, MetaboliteIdentifier)
    )

    assert primary.chem_props == []
