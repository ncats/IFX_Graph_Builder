from pathlib import Path
import json
import zipfile

from src.input_adapters.metabolite_harmonization.hmdb import HmdbMetaboliteEquivalenceAdapter
from src.input_adapters.metabolite_harmonization.hmdb_classes_ontology import (
    HmdbMetaboliteClassificationAdapter,
    HmdbOntologyAdapter,
)
from src.interfaces.output_adapter import OutputAdapter
from src.models.metabolite_harmonization import (
    HmdbMetaboliteOntologyEdge,
    HmdbOntologyParentEdge,
    HmdbOntologyTerm,
    MetaboliteClassificationEdge,
    MetaboliteClassificationParentEdge,
    MetaboliteClassificationTerm,
    MetaboliteIdentifier,
    MetaboliteIdentifierMappingEdge,
)
from src.shared.record_merger import FieldConflictBehavior


class _ConvertingOutputAdapter(OutputAdapter):
    def store(self, objects, field_conflict_behavior=FieldConflictBehavior.KeepFirst) -> bool:
        return True

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        return True


def _write_hmdb_zip(path: Path):
    xml = """<?xml version="1.0" encoding="UTF-8"?>
<hmdb xmlns="http://www.hmdb.ca">
<metabolite>
  <accession>HMDB0000001</accession>
  <name>1-Methylhistidine</name>
  <synonyms>
    <synonym>Pi-methylhistidine</synonym>
    <synonym>1-MHis</synonym>
    <synonym>1-MHis</synonym>
  </synonyms>
  <secondary_accessions>
    <accession>HMDB00001</accession>
    <accession>HMDB0004935</accession>
    <accession>HMDB0004935</accession>
  </secondary_accessions>
  <cas_registry_number>332-80-9</cas_registry_number>
  <chemspider_id>83153</chemspider_id>
  <drugbank_id>DB04151</drugbank_id>
  <foodb_id>FDB093588</foodb_id>
  <pubchem_compound_id>92105</pubchem_compound_id>
  <chebi_id>50599</chebi_id>
  <kegg_id>C01152</kegg_id>
  <biocyc_id>CPD-1</biocyc_id>
  <bigg_id>bigg1</bigg_id>
  <metlin_id>3741</metlin_id>
  <taxonomy>
    <super_class>Organic acids and derivatives</super_class>
    <class>Carboxylic acids and derivatives</class>
    <sub_class>Amino acids, peptides, and analogues</sub_class>
  </taxonomy>
  <ontology>
    <root>
      <term>Disposition</term>
      <definition>A concept that describes the origin of a chemical.</definition>
      <parent_id></parent_id>
      <level>1</level>
      <type>parent</type>
      <descendants>
        <descendant>
          <term>Source</term>
          <definition>Natural or synthetic origin of a chemical.</definition>
          <parent_id>7724</parent_id>
          <level>2</level>
          <type>parent</type>
          <synonyms/>
          <descendants>
            <descendant>
              <term>Endogenous</term>
              <definition></definition>
              <parent_id>7735</parent_id>
              <level>3</level>
              <type>child</type>
              <synonyms/>
            </descendant>
          </descendants>
        </descendant>
        <descendant>
          <term>Biological location</term>
          <definition>Physiological origin within an organism.</definition>
          <parent_id>7724</parent_id>
          <level>2</level>
          <type>parent</type>
          <synonyms/>
          <descendants>
            <descendant>
              <term>Biofluid and excreta</term>
              <definition>A liquid, semi-solid or solid material originating in the body.</definition>
              <parent_id>7725</parent_id>
              <level>3</level>
              <type>parent</type>
              <synonyms/>
              <descendants>
                <descendant>
                  <term>Blood</term>
                  <definition>Blood definition.</definition>
                  <parent_id>7731</parent_id>
                  <level>4</level>
                  <type>child</type>
                  <synonyms>
                    <synonym>Whole blood</synonym>
                    <synonym>Whole blood</synonym>
                  </synonyms>
                </descendant>
              </descendants>
            </descendant>
          </descendants>
        </descendant>
        <descendant>
          <term>Pathway</term>
          <definition>Pathway context outside the RaMP ontology categories.</definition>
          <parent_id>7724</parent_id>
          <level>2</level>
          <type>parent</type>
          <synonyms/>
          <descendants>
            <descendant>
              <term>Acetaminophen Metabolism Pathway</term>
              <definition></definition>
              <parent_id>9001</parent_id>
              <level>3</level>
              <type>child</type>
              <synonyms/>
            </descendant>
          </descendants>
        </descendant>
      </descendants>
    </root>
  </ontology>
</metabolite>
<metabolite>
  <accession>HMDB0000002</accession>
  <name>1,3-Diaminopropane</name>
  <synonyms/>
  <secondary_accessions>
    <accession>HMDB00002</accession>
  </secondary_accessions>
  <pubchem_compound_id>92105</pubchem_compound_id>
  <chebi_id></chebi_id>
  <taxonomy>
    <super_class>Organic nitrogen compounds</super_class>
    <class>Organonitrogen compounds</class>
  </taxonomy>
</metabolite>
</hmdb>
"""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("hmdb_metabolites.xml", xml)


def _records(adapter):
    return [record for batch in adapter.get_all() for record in batch]


def test_hmdb_metabolite_harmonization_adapter_emits_nodes_and_equivalence_edges(tmp_path: Path):
    zip_path = tmp_path / "hmdb_metabolites.zip"
    _write_hmdb_zip(zip_path)
    adapter = HmdbMetaboliteEquivalenceAdapter(hmdb_zip_file=str(zip_path))

    records = _records(adapter)
    nodes = [record for record in records if isinstance(record, MetaboliteIdentifier)]
    edges = [record for record in records if isinstance(record, MetaboliteIdentifierMappingEdge)]
    nodes_by_id = {node.id: node for node in nodes}
    edge_pairs = {(edge.start_node.id, edge.end_node.id) for edge in edges}

    assert nodes_by_id["HMDB:HMDB0000001"].names[0].value == "1-Methylhistidine"
    assert nodes_by_id["HMDB:HMDB0000001"].names[0].source == "HMDB"
    assert nodes_by_id["HMDB:HMDB0000001"].prefix == "HMDB"
    assert [synonym.value for synonym in nodes_by_id["HMDB:HMDB0000001"].synonyms] == [
        "Pi-methylhistidine",
        "1-MHis",
    ]
    assert {name.value for node in nodes if node.id == "PUBCHEM.COMPOUND:92105"
            for name in node.names} == {"1-Methylhistidine", "1,3-Diaminopropane"}
    assert all(name.source == "HMDB" for node in nodes if node.id == "PUBCHEM.COMPOUND:92105"
               for name in node.names)
    assert nodes_by_id["PUBCHEM.COMPOUND:92105"].synonyms == []
    assert nodes_by_id["PUBCHEM.COMPOUND:92105"].prefix == "PUBCHEM.COMPOUND"

    expected_node_ids = {
        "HMDB:HMDB0000001",
        "HMDB:HMDB00001",
        "HMDB:HMDB0004935",
        "CAS:332-80-9",
        "ChemSpider:83153",
        "DRUGBANK:DB04151",
        "FoodDB:FDB093588",
        "PUBCHEM.COMPOUND:92105",
        "CHEBI:50599",
        "KEGG.COMPOUND:C01152",
        "BioCyc:CPD-1",
        "BiGG:bigg1",
        "METLIN:3741",
        "HMDB:HMDB0000002",
        "HMDB:HMDB00002",
    }
    assert expected_node_ids <= set(nodes_by_id)

    assert ("HMDB:HMDB0000001", "CHEBI:50599") in edge_pairs
    assert ("HMDB:HMDB0000001", "PUBCHEM.COMPOUND:92105") in edge_pairs
    assert ("HMDB:HMDB0000001", "KEGG.COMPOUND:C01152") in edge_pairs
    assert ("HMDB:HMDB0000001", "HMDB:HMDB0004935") in edge_pairs
    assert ("HMDB:HMDB0000002", "PUBCHEM.COMPOUND:92105") in edge_pairs
    assert ("HMDB:HMDB0000001", "DRUGBANK:DB04151") in edge_pairs
    assert ("HMDB:HMDB0000001", "FoodDB:FDB093588") in edge_pairs
    assert ("HMDB:HMDB0000001", "BiGG:bigg1") in edge_pairs

    details = {
        (edge.start_node.id, edge.end_node.id): edge.details[0]
        for edge in edges
    }
    assert details[("HMDB:HMDB0000001", "CHEBI:50599")].source == "HMDB"
    assert details[("HMDB:HMDB0000001", "CHEBI:50599")].source_field == "chebi_id"
    assert details[("HMDB:HMDB0000001", "CHEBI:50599")].source_id == "HMDB:HMDB0000001"

    assert len([
        edge
        for edge in edges
        if edge.start_node.id == "HMDB:HMDB0000001" and edge.end_node.id == "HMDB:HMDB0004935"
    ]) == 1


def test_hmdb_records_are_json_serializable_after_output_conversion(tmp_path: Path):
    zip_path = tmp_path / "hmdb_metabolites.zip"
    _write_hmdb_zip(zip_path)
    adapter = HmdbMetaboliteEquivalenceAdapter(hmdb_zip_file=str(zip_path), max_records=1)
    output = _ConvertingOutputAdapter()

    converted_groups = output.sort_and_convert_objects(_records(adapter))
    converted_records = [
        record
        for group in converted_groups.values()
        for record in group[0]
    ]

    json.dumps(converted_records)
    hmdb_node = next(record for record in converted_records if record["id"] == "HMDB:HMDB0000001")
    chebi_edge = next(record for record in converted_records if record.get("end_id") == "CHEBI:50599")
    assert hmdb_node["names"] == [
        {"value": "1-Methylhistidine", "source": "HMDB", "source_field": "name"}
    ]
    assert hmdb_node["prefix"] == "HMDB"
    assert chebi_edge["details"] == [
        {"source": "HMDB", "source_field": "chebi_id", "source_id": "HMDB:HMDB0000001"}
    ]


def test_hmdb_adapter_honors_max_records(tmp_path: Path):
    zip_path = tmp_path / "hmdb_metabolites.zip"
    _write_hmdb_zip(zip_path)
    adapter = HmdbMetaboliteEquivalenceAdapter(hmdb_zip_file=str(zip_path), max_records=1)

    records = _records(adapter)
    assert any(isinstance(record, MetaboliteIdentifier) and record.id == "HMDB:HMDB0000001" for record in records)
    assert not any(isinstance(record, MetaboliteIdentifier) and record.id == "HMDB:HMDB0000002" for record in records)


def test_hmdb_metabolite_classification_adapter_emits_terms_hierarchy_and_edges(tmp_path: Path):
    zip_path = tmp_path / "hmdb_metabolites.zip"
    _write_hmdb_zip(zip_path)
    adapter = HmdbMetaboliteClassificationAdapter(hmdb_zip_file=str(zip_path))

    records = _records(adapter)
    terms = [record for record in records if isinstance(record, MetaboliteClassificationTerm)]
    parent_edges = [record for record in records if isinstance(record, MetaboliteClassificationParentEdge)]
    member_edges = [record for record in records if isinstance(record, MetaboliteClassificationEdge)]
    term_by_name = {(term.level_name, term.name): term for term in terms}
    parent_pairs = {(edge.start_node.id, edge.end_node.id) for edge in parent_edges}
    member_pairs = {(edge.start_node.id, edge.end_node.id) for edge in member_edges}

    super_term = term_by_name[("ClassyFire_super_class", "Organic acids and derivatives")]
    class_term = term_by_name[("ClassyFire_class", "Carboxylic acids and derivatives")]
    sub_class_term = term_by_name[("ClassyFire_sub_class", "Amino acids, peptides, and analogues")]

    assert super_term.id == "HMDB.CLASS:ClassyFire_super_class:Organic acids and derivatives"
    assert super_term.source == "HMDB"
    assert (super_term.id, class_term.id) in parent_pairs
    assert (class_term.id, sub_class_term.id) in parent_pairs
    assert ("HMDB:HMDB0000001", super_term.id) not in member_pairs
    assert ("HMDB:HMDB0000001", class_term.id) not in member_pairs
    assert ("HMDB:HMDB0000001", sub_class_term.id) in member_pairs
    assert any(
        edge.start_node.id == "HMDB:HMDB0000002" and edge.end_node.id.endswith(":Organonitrogen compounds")
        for edge in member_edges
    )
    assert member_edges[0].details[0].source == "HMDB"


def test_hmdb_ontology_adapter_emits_terms_parent_edges_and_memberships(tmp_path: Path):
    zip_path = tmp_path / "hmdb_metabolites.zip"
    _write_hmdb_zip(zip_path)
    adapter = HmdbOntologyAdapter(hmdb_zip_file=str(zip_path))

    records = _records(adapter)
    terms = [record for record in records if isinstance(record, HmdbOntologyTerm)]
    parent_edges = [record for record in records if isinstance(record, HmdbOntologyParentEdge)]
    member_edges = [record for record in records if isinstance(record, HmdbMetaboliteOntologyEdge)]
    terms_by_name = {(term.ontology_type, term.name): term for term in terms}
    parent_pairs = {(edge.start_node.id, edge.end_node.id) for edge in parent_edges}
    member_pairs = {(edge.start_node.id, edge.end_node.id) for edge in member_edges}

    source = terms_by_name[("Source", "Source")]
    endogenous = terms_by_name[("Source", "Endogenous")]
    biofluid = terms_by_name[("Biofluid and excreta", "Biofluid and excreta")]
    blood = terms_by_name[("Biofluid and excreta", "Blood")]

    assert source.id == "HMDB.ONTOLOGY:Source:Source"
    assert endogenous.ontology_type == "Source"
    assert blood.synonyms == ["Whole blood"]
    assert (source.id, endogenous.id) in parent_pairs
    assert (biofluid.id, blood.id) in parent_pairs
    assert ("HMDB:HMDB0000001", endogenous.id) in member_pairs
    assert ("HMDB:HMDB0000001", blood.id) in member_pairs
    assert all(
        edge.details[0].source_id == edge.start_node.id
        for edge in member_edges
    )
    assert all(edge.end_node.id != source.id for edge in member_edges)
    assert ("Pathway", "Pathway") not in terms_by_name
    assert ("Pathway", "Acetaminophen Metabolism Pathway") not in terms_by_name
