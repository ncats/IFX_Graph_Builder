from pathlib import Path
import gzip
import json

from src.input_adapters.metabolite_harmonization.expasy import ExpasyEnzymeClassAdapter
from src.input_adapters.metabolite_harmonization.rhea import RheaReactionAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite_harmonization import (
    MetaboliteIdentifier,
    ProteinIdentifier,
    RheaMetaboliteReactionEdge,
    RheaProteinReactionEdge,
    RheaReaction,
    RheaReactionClass,
    RheaReactionClassParentEdge,
    RheaReactionDirectionEdge,
    RheaReactionReactionClassEdge,
)


def _records(adapter):
    return [record for batch in adapter.get_all() for record in batch]


def test_expasy_enzyme_class_adapter_emits_classes_and_parent_edges(tmp_path: Path):
    enzclass = tmp_path / "enzclass.txt"
    enzclass.write_text(
        "1. -. -.-  Oxidoreductases.\n"
        "1. 1. -.-   Acting on the CH-OH group of donors.\n"
        "1. 1. 1.-    With NAD(+) or NADP(+) as acceptor.\n",
        encoding="utf-8",
    )
    enzyme_dat = tmp_path / "enzyme.dat"
    enzyme_dat.write_text(
        "ID   1.1.1.1\n"
        "DE   alcohol dehydrogenase.\n"
        "//\n"
        "ID   1.1.1.n1\n"
        "DE   non-standard provisional entry.\n"
        "//\n",
        encoding="utf-8",
    )

    records = _records(
        ExpasyEnzymeClassAdapter(
            enzclass_file=str(enzclass),
            enzyme_dat_file=str(enzyme_dat),
        )
    )

    classes = {record.id: record for record in records if isinstance(record, RheaReactionClass)}
    parent_edges = {
        (record.start_node.id, record.end_node.id)
        for record in records
        if isinstance(record, RheaReactionClassParentEdge)
    }

    assert classes["EC:1.1.1.1"].name == "alcohol dehydrogenase."
    assert classes["EC:1.1.1.n1"].name == "non-standard provisional entry."
    assert classes["EC:1.1.1.-"].name == "With NAD(+) or NADP(+) as acceptor."
    assert parent_edges == {
        ("EC:1.1.1.1", "EC:1.1.1.-"),
        ("EC:1.1.1.n1", "EC:1.1.1.-"),
        ("EC:1.1.1.-", "EC:1.1.-.-"),
        ("EC:1.1.-.-", "EC:1.-.-.-"),
    }


def test_rhea_reaction_adapter_emits_active_reactions_and_context_edges(tmp_path: Path):
    rdf = tmp_path / "rhea.rdf"
    rdf.write_text(
        """<?xml version="1.0" encoding="UTF-8"?>
<rdf:RDF xmlns:rdf="http://www.w3.org/1999/02/22-rdf-syntax-ns#"
         xmlns:rdfs="http://www.w3.org/2000/01/rdf-schema#"
         xmlns:rh="http://rdf.rhea-db.org/">
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000">
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/Reaction"/>
    <rh:accession>RHEA:10000</rh:accession>
    <rdfs:label>test reaction</rdfs:label>
    <rh:equation>water = generic</rh:equation>
    <rh:htmlEquation>water = generic</rh:htmlEquation>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
    <rh:isTransport rdf:datatype="http://www.w3.org/2001/XMLSchema#boolean">false</rh:isTransport>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10001">
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/DirectionalReaction"/>
    <rh:accession>RHEA:10001</rh:accession>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10002">
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/DirectionalReaction"/>
    <rh:accession>RHEA:10002</rh:accession>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10003">
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/BidirectionalReaction"/>
    <rh:accession>RHEA:10003</rh:accession>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Approved"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10004">
    <rdfs:subClassOf rdf:resource="http://rdf.rhea-db.org/Reaction"/>
    <rh:accession>RHEA:10004</rh:accession>
    <rh:status rdf:resource="http://rdf.rhea-db.org/Obsolete"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000_L">
    <rh:contains1 rdf:resource="http://rdf.rhea-db.org/Participant_10000_compound_1"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Participant_10000_compound_1">
    <rh:compound rdf:resource="http://rdf.rhea-db.org/Compound_1"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_1">
    <rh:accession>CHEBI:15377</rh:accession>
    <rh:name>H2O</rh:name>
    <rh:htmlName>H&lt;sub&gt;2&lt;/sub&gt;O</rh:htmlName>
    <rh:formula>H2O</rh:formula>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/10000_R">
    <rh:contains2 rdf:resource="http://rdf.rhea-db.org/Participant_10000_compound_2"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Participant_10000_compound_2">
    <rh:compound rdf:resource="http://rdf.rhea-db.org/Compound_2"/>
  </rdf:Description>
  <rdf:Description rdf:about="http://rdf.rhea-db.org/Compound_2">
    <rh:accession>GENERIC:10594</rh:accession>
    <rh:name>[protein]-dithiol</rh:name>
    <rh:htmlName>[protein]-dithiol</rh:htmlName>
  </rdf:Description>
</rdf:RDF>
""",
        encoding="utf-8",
    )
    directions = tmp_path / "rhea-directions.tsv"
    directions.write_text("RHEA_ID_MASTER\tRHEA_ID_LR\tRHEA_ID_RL\tRHEA_ID_BI\n10000\t10001\t10002\t10003\n", encoding="utf-8")
    rhea2ec = tmp_path / "rhea2ec.tsv"
    rhea2ec.write_text("RHEA_ID\tDIRECTION\tMASTER_ID\tID\n10000\tUN\t10000\t1.1.1.1\n", encoding="utf-8")
    sprot = tmp_path / "rhea2uniprot_sprot.tsv"
    sprot.write_text("RHEA_ID\tDIRECTION\tMASTER_ID\tID\n10000\tUN\t10000\tP00001\n10000\tUN\t10000\tNOHUMAN\n", encoding="utf-8")
    trembl = tmp_path / "rhea2uniprot_trembl.tsv.gz"
    with gzip.open(trembl, "wt", encoding="utf-8") as handle:
        handle.write("RHEA_ID\tDIRECTION\tMASTER_ID\tID\n10001\tLR\t10000\tOLDP2\n")
    uniprot = tmp_path / "uniprot-human.json.gz"
    with gzip.open(uniprot, "wt", encoding="utf-8") as handle:
        json.dump(
            {
                "results": [
                    {"primaryAccession": "P00001", "secondaryAccessions": [], "entryType": "UniProtKB reviewed (Swiss-Prot)"},
                    {"primaryAccession": "P00002", "secondaryAccessions": ["OLDP2"], "entryType": "UniProtKB unreviewed (TrEMBL)"},
                ]
            },
            handle,
        )
    rdf_gz = tmp_path / "rhea.rdf.gz"
    with gzip.open(rdf_gz, "wb") as handle:
        handle.write(rdf.read_bytes())

    class RegistryDataset:
        def file(self, name: str) -> Path:
            return {
                "rhea.rdf.gz": rdf_gz,
                "rhea-directions.tsv": directions,
                "rhea2ec.tsv": rhea2ec,
                "rhea2uniprot_sprot.tsv": sprot,
                "rhea2uniprot_trembl.tsv.gz": trembl,
            }[name]

        def version_info(self) -> DatasourceVersionInfo:
            return DatasourceVersionInfo(
                version="142",
                version_date=None,
                download_date=None,
            )

    records = _records(
        RheaReactionAdapter(
            data_source=RegistryDataset(),
            uniprot_human_file=str(uniprot),
        )
    )

    reaction_nodes = {record.id: record for record in records if isinstance(record, RheaReaction)}
    metabolite_nodes = {record.id for record in records if isinstance(record, MetaboliteIdentifier)}
    protein_nodes = {record.id: record for record in records if isinstance(record, ProteinIdentifier)}
    direction_edges = {
        (record.start_node.id, record.end_node.id, record.variant_direction)
        for record in records
        if isinstance(record, RheaReactionDirectionEdge)
    }
    metabolite_edges = [
        record for record in records if isinstance(record, RheaMetaboliteReactionEdge)
    ]
    protein_edges = {
        (record.start_node.id, record.end_node.id, record.source_file)
        for record in records
        if isinstance(record, RheaProteinReactionEdge)
    }
    class_edges = {
        (record.start_node.id, record.end_node.id)
        for record in records
        if isinstance(record, RheaReactionReactionClassEdge)
    }

    assert set(reaction_nodes) == {"RHEA:10000", "RHEA:10001", "RHEA:10002", "RHEA:10003"}
    assert reaction_nodes["RHEA:10002"].direction == "RL"
    assert metabolite_nodes == {"CHEBI:15377", "RHEA.COMP:10594"}
    assert protein_nodes["UniProtKB:P00001"].is_reviewed is True
    assert protein_nodes["UniProtKB:OLDP2"].is_reviewed is False
    assert 'UniProtKB:P00002' not in protein_nodes
    assert all(e.source_id == e.start_node.id for e in records if isinstance(e, RheaProteinReactionEdge))
    assert "UniProtKB:NOHUMAN" not in protein_nodes
    assert direction_edges == {
        ("RHEA:10000", "RHEA:10001", "LR"),
        ("RHEA:10000", "RHEA:10002", "RL"),
        ("RHEA:10000", "RHEA:10003", "BD"),
    }
    assert ("CHEBI:15377", "RHEA:10000", "left", 1) in {
        (edge.start_node.id, edge.end_node.id, edge.side, edge.coefficient)
        for edge in metabolite_edges
    }
    assert ("CHEBI:15377", "RHEA:10002", "right", 1) in {
        (edge.start_node.id, edge.end_node.id, edge.side, edge.coefficient)
        for edge in metabolite_edges
    }
    assert any(edge.name == "H2O" and edge.formula == "H2O" for edge in metabolite_edges)
    assert all(edge.source_id == edge.start_node.id for edge in metabolite_edges)
    assert protein_edges == {
        ("UniProtKB:P00001", "RHEA:10000", "rhea2uniprot_sprot.tsv"),
        ("UniProtKB:OLDP2", "RHEA:10001", "rhea2uniprot_trembl.tsv.gz"),
    }
    assert class_edges == {("RHEA:10000", "EC:1.1.1.1")}
