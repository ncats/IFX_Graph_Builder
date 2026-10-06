from pathlib import Path
import json
import zipfile

from src.input_adapters.metabolite_harmonization.wikipathways import (
    WikiPathwaysMetaboliteEquivalenceAdapter,
    WikiPathwaysPathwayContextAdapter,
)
from src.interfaces.output_adapter import OutputAdapter
from src.models.metabolite_harmonization import (
    GeneIdentifier,
    GenePathwayEdge,
    MetaboliteIdentifier,
    MetaboliteIdentifierMappingEdge,
    MetabolitePathwayEdge,
    ProteinIdentifier,
    ProteinPathwayEdge,
)
from src.shared.record_merger import FieldConflictBehavior


class _ConvertingOutputAdapter(OutputAdapter):
    def store(self, objects, single_source=False, field_conflict_behavior=FieldConflictBehavior.KeepFirst) -> bool:
        return True

    def create_or_truncate_datastore(self, truncate_tables: bool = None) -> bool:
        return True


def _write_wikipathways_zip(path: Path):
    ttl = """
@prefix dc:      <http://purl.org/dc/elements/1.1/> .
@prefix dcterms: <http://purl.org/dc/terms/> .
@prefix obo:     <http://purl.obolibrary.org/obo/> .
@prefix rdf:     <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#> .
@prefix wp:      <http://vocabularies.wikipathways.org/wp#> .

<https://identifiers.org/wikipathways/WP100_r1>
        rdf:type       wp:Pathway ;
        wp:organism   obo:NCBITaxon_9606 .

<https://identifiers.org/kegg.compound/C00051>
        rdf:type            wp:DataNode , wp:Metabolite ;
        rdfs:label          "Glutathione (reduced)" ;
        dc:identifier       <https://identifiers.org/kegg.compound/C00051> ;
        dcterms:identifier  "C00051" ;
        wp:bdbChEBI         <https://identifiers.org/chebi/CHEBI:16856> ;
        wp:bdbChemspider    <https://identifiers.org/chemspider/111188> ;
        wp:bdbHmdb          <https://identifiers.org/hmdb/HMDB0000125> ;
        wp:bdbInChIKey      <https://identifiers.org/inchikey/RWSXRVCMGQZWBV-WDSKDSINSA-N> ;
        wp:bdbKeggCompound  <https://identifiers.org/kegg.compound/C00051> ;
        wp:bdbPubChem       <http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID25246407> ,
                            <http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID124886> ;
        wp:bdbWikidata      <http://www.wikidata.org/entity/Q116907> .

<https://identifiers.org/pubchem.compound/753>
        rdf:type            wp:DataNode , wp:Metabolite ;
        rdfs:label          "glycine" ;
        dc:identifier       <https://identifiers.org/pubchem.compound/753> ;
        wp:bdbPubChem       <http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID753> ;
        wp:bdbChEBI         <https://identifiers.org/chebi/CHEBI:17754> .

<https://identifiers.org/lipidmaps/LMFA01170120>
        rdf:type            wp:DataNode , wp:Metabolite ;
        rdfs:label          "oxaloacetate" ;
        wp:bdbLipidMaps     <https://identifiers.org/lipidmaps/LMFA01170120> ;
        wp:bdbReactome      <https://identifiers.org/reactome/R-HSA-194002> ;
        wp:bdbPubChem       <http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID970> .

<https://identifiers.org/knapsack/C00000128>
        rdf:type            wp:DataNode , wp:Metabolite ;
        rdfs:label          "source oddity" ;
        wp:bdbKeggGlycan    <https://identifiers.org/kegg.glycan/G00170> ;
        wp:bdbKeggDrug      <https://identifiers.org/kegg.drug/D11487> ;
        wp:bdbChEMBL        <https://identifiers.org/chembl.compound/CHEMBL1207858> ;
        wp:bdbPubChemSubstance <https://identifiers.org/pubchem.substance/136368185> ;
        wp:bdbPharmGKB      <https://identifiers.org/pharmgkb.drug/PA162373091> ;
        wp:bdbPidPathway    <https://identifiers.org/pid.pathway/9606> .

<https://identifiers.org/ncbigene/1234>
        rdf:type            wp:DataNode , wp:GeneProduct ;
        wp:bdbEnsembl       <https://identifiers.org/ensembl/ENSG000001234> ;
        wp:bdbUniprot       <https://identifiers.org/uniprot/P99999> .

<https://identifiers.org/uniprot/P12345>
        rdf:type            wp:DataNode , wp:Protein ;
        wp:bdbEntrezGene    <https://identifiers.org/ncbigene/5678> ;
        wp:bdbWikidata      <http://www.wikidata.org/entity/Q12345> .
"""
    with zipfile.ZipFile(path, "w") as archive:
        archive.writestr("wp/WP100.ttl", ttl)
        archive.writestr(
            "wp/WP101.ttl",
            """
@prefix obo:     <http://purl.obolibrary.org/obo/> .
@prefix rdf:     <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#> .
@prefix wp:      <http://vocabularies.wikipathways.org/wp#> .

<https://identifiers.org/wikipathways/WP101_r1>
        rdf:type       wp:Pathway ;
        wp:organism   obo:NCBITaxon_9606 .

<https://identifiers.org/kegg.compound/C00051>
        rdf:type            wp:DataNode , wp:Metabolite ;
        rdfs:label          "reduced glutathione" ;
        wp:bdbChEBI         <https://identifiers.org/chebi/CHEBI:16856> .
""",
        )
        archive.writestr(
            "wp/WP478.ttl",
            """
@prefix obo:     <http://purl.obolibrary.org/obo/> .
@prefix rdf:     <http://www.w3.org/1999/02/22-rdf-syntax-ns#> .
@prefix rdfs:    <http://www.w3.org/2000/01/rdf-schema#> .
@prefix wp:      <http://vocabularies.wikipathways.org/wp#> .

<https://identifiers.org/wikipathways/WP478_r1>
        rdf:type       wp:Pathway ;
        wp:organism   obo:NCBITaxon_4932 .

<https://identifiers.org/cas/9050-36-6>
        rdf:type            wp:DataNode , wp:Metabolite ;
        rdfs:label          "maltodextrin" ;
        wp:bdbChEBI         <https://identifiers.org/chebi/CHEBI:18398> ;
        wp:bdbPubChem       <http://rdf.ncbi.nlm.nih.gov/pubchem/compound/CID62698> .
""",
        )


def _records(adapter):
    return [record for batch in adapter.get_all() for record in batch]


def test_wikipathways_metabolite_adapter_emits_nodes_labels_and_edges(tmp_path: Path):
    zip_path = tmp_path / "wikipathways-rdf-wp.zip"
    _write_wikipathways_zip(zip_path)
    adapter = WikiPathwaysMetaboliteEquivalenceAdapter(rdf_zip_file=str(zip_path))

    records = _records(adapter)
    nodes = [record for record in records if isinstance(record, MetaboliteIdentifier)]
    edges = [record for record in records if isinstance(record, MetaboliteIdentifierMappingEdge)]
    nodes_by_id = {node.id: node for node in nodes}
    edge_pairs = {(edge.start_node.id, edge.end_node.id) for edge in edges}

    kegg_names = [
        node.names[0].value
        for node in nodes
        if node.id == "KEGG.COMPOUND:C00051" and node.names
    ]
    assert "Glutathione (reduced)" in kegg_names
    assert "reduced glutathione" in kegg_names
    assert all(
        node.names[0].source == "WikiPathways"
        for node in nodes
        if node.id == "KEGG.COMPOUND:C00051" and node.names
    )
    assert nodes_by_id["PUBCHEM.COMPOUND:25246407"].names[0].value == "Glutathione (reduced)"
    assert nodes_by_id["PUBCHEM.COMPOUND:25246407"].names[0].source == "WikiPathways"
    assert nodes_by_id["PUBCHEM.COMPOUND:25246407"].prefix == "PUBCHEM.COMPOUND"
    assert nodes_by_id["PUBCHEM.COMPOUND:753"].names[0].value == "glycine"
    assert nodes_by_id["KEGG.COMPOUND:C00051"].prefix == "KEGG.COMPOUND"

    expected_node_ids = {
        "KEGG.COMPOUND:C00051",
        "CHEBI:16856",
        "ChemSpider:111188",
        "HMDB:HMDB0000125",
        "PUBCHEM.COMPOUND:25246407",
        "PUBCHEM.COMPOUND:124886",
        "Wikidata:Q116907",
        "PUBCHEM.COMPOUND:753",
        "CHEBI:17754",
        "LIPIDMAPS:LMFA01170120",
        "Reactome:R-HSA-194002",
        "PUBCHEM.COMPOUND:970",
        "KNApSAcK:C00000128",
        "KEGG.GLYCAN:G00170",
        "KEGG.DRUG:D11487",
        "ChEMBL.COMPOUND:CHEMBL1207858",
        "PUBCHEM.SUBSTANCE:136368185",
        "PharmGKB.DRUG:PA162373091",
        "PID.PATHWAY:9606",
    }
    assert expected_node_ids <= set(nodes_by_id)
    assert "InChIKey:RWSXRVCMGQZWBV-WDSKDSINSA-N" not in nodes_by_id

    assert ("KEGG.COMPOUND:C00051", "CHEBI:16856") in edge_pairs
    assert ("KEGG.COMPOUND:C00051", "PUBCHEM.COMPOUND:25246407") in edge_pairs
    assert ("PUBCHEM.COMPOUND:753", "CHEBI:17754") in edge_pairs
    assert ("LIPIDMAPS:LMFA01170120", "Reactome:R-HSA-194002") in edge_pairs
    assert ("KNApSAcK:C00000128", "PID.PATHWAY:9606") in edge_pairs
    assert ("KEGG.COMPOUND:C00051", "InChIKey:RWSXRVCMGQZWBV-WDSKDSINSA-N") not in edge_pairs
    assert ("KEGG.COMPOUND:C00051", "KEGG.COMPOUND:C00051") not in edge_pairs
    assert "CAS:9050-36-6" not in nodes_by_id
    assert ("CAS:9050-36-6", "CHEBI:18398") not in edge_pairs
    assert ("CAS:9050-36-6", "PUBCHEM.COMPOUND:62698") not in edge_pairs

    details = {
        (edge.start_node.id, edge.end_node.id): edge.details[0]
        for edge in edges
    }
    assert details[("KEGG.COMPOUND:C00051", "CHEBI:16856")].source == "WikiPathways"
    assert details[("KEGG.COMPOUND:C00051", "CHEBI:16856")].source_field == "bdbChEBI"
    assert details[("KEGG.COMPOUND:C00051", "CHEBI:16856")].source_id == "KEGG.COMPOUND:C00051"


def test_wikipathways_records_are_json_serializable_after_output_conversion(tmp_path: Path):
    zip_path = tmp_path / "wikipathways-rdf-wp.zip"
    _write_wikipathways_zip(zip_path)
    adapter = WikiPathwaysMetaboliteEquivalenceAdapter(rdf_zip_file=str(zip_path), max_records=1)
    output = _ConvertingOutputAdapter()

    converted_groups = output.sort_and_convert_objects(_records(adapter))
    converted_records = [
        record
        for group in converted_groups.values()
        for record in group[0]
    ]

    json.dumps(converted_records)
    source_node = next(record for record in converted_records if record["id"] == "KEGG.COMPOUND:C00051")
    chebi_edge = next(record for record in converted_records if record.get("end_id") == "CHEBI:16856")
    assert source_node["names"] == [
        {"value": "Glutathione (reduced)", "source": "WikiPathways", "source_field": "rdfs:label"}
    ]
    assert source_node["prefix"] == "KEGG.COMPOUND"
    assert chebi_edge["details"] == [
        {"source": "WikiPathways", "source_field": "bdbChEBI", "source_id": "KEGG.COMPOUND:C00051"}
    ]


def test_wikipathways_adapter_honors_max_records(tmp_path: Path):
    zip_path = tmp_path / "wikipathways-rdf-wp.zip"
    _write_wikipathways_zip(zip_path)
    adapter = WikiPathwaysMetaboliteEquivalenceAdapter(rdf_zip_file=str(zip_path), max_records=1)

    records = _records(adapter)
    assert any(isinstance(record, MetaboliteIdentifier) and record.id == "KEGG.COMPOUND:C00051" for record in records)
    assert not any(isinstance(record, MetaboliteIdentifier) and record.id == "PUBCHEM.COMPOUND:753" for record in records)


def test_wikipathways_pathway_context_uses_only_source_native_identifiers(tmp_path: Path):
    zip_path = tmp_path / "wikipathways-rdf-wp.zip"
    _write_wikipathways_zip(zip_path)
    adapter = WikiPathwaysPathwayContextAdapter(rdf_zip_file=str(zip_path))

    records = _records(adapter)
    metabolite_nodes = {
        record.id
        for record in records
        if isinstance(record, MetaboliteIdentifier)
    }
    metabolite_pathway_edges = {
        edge.start_node.id
        for edge in records
        if isinstance(edge, MetabolitePathwayEdge)
    }

    expected_primary_ids = {
        "KEGG.COMPOUND:C00051",
        "PUBCHEM.COMPOUND:753",
        "LIPIDMAPS:LMFA01170120",
        "KNApSAcK:C00000128",
    }
    xref_only_ids = {
        "CHEBI:16856",
        "HMDB:HMDB0000125",
        "PUBCHEM.COMPOUND:25246407",
        "Wikidata:Q116907",
        "CHEBI:17754",
        "Reactome:R-HSA-194002",
        "PUBCHEM.COMPOUND:970",
        "KEGG.GLYCAN:G00170",
        "KEGG.DRUG:D11487",
    }

    assert metabolite_nodes == expected_primary_ids
    assert metabolite_pathway_edges == expected_primary_ids
    assert metabolite_nodes.isdisjoint(xref_only_ids)
    assert metabolite_pathway_edges.isdisjoint(xref_only_ids)

    gene_nodes = {record.id for record in records if isinstance(record, GeneIdentifier)}
    gene_pathway_edges = {
        record.start_node.id for record in records if isinstance(record, GenePathwayEdge)
    }
    protein_nodes = {record.id for record in records if isinstance(record, ProteinIdentifier)}
    protein_pathway_edges = {
        record.start_node.id for record in records if isinstance(record, ProteinPathwayEdge)
    }
    assert gene_nodes == {"NCBIGene:1234"}
    assert gene_pathway_edges == {"NCBIGene:1234"}
    assert protein_nodes == {"UniProtKB:P12345"}
    assert protein_pathway_edges == {"UniProtKB:P12345"}

    pathway_detail = next(
        edge.details[0]
        for edge in records
        if isinstance(edge, MetabolitePathwayEdge)
        and edge.start_node.id == "KEGG.COMPOUND:C00051"
    )
    assert pathway_detail.source_id == "KEGG.COMPOUND:C00051"
