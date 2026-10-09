import json

import pytest

from src.shared.targetgraph_parser import TargetGraphGeneParser, TargetGraphProteinParser, TargetGraphTranscriptParser
from src.shared.targetgraph_parser import TargetGraphAddtlProteinIDParser
from src.input_adapters.target_graph.gene_transcript_edge import GeneTranscriptEdgeAdapter
from src.input_adapters.target_graph.protein_nodes_and_edges import GeneProteinEdgeAdapter
from src.id_resolvers.target_graph_resolver import (
    TCRDTargetResolver, TargetGraphProteinResolver, TargetGraphResolver,
)
from src.output_adapters.arango_output_adapter import ArangoOutputAdapter


def test_targetgraph_gene_parser_reads_tsv(tmp_path):
    file_path = tmp_path / "gene_ids.tsv"
    file_path.write_text(
        "\t".join(
            [
                "ncats_gene_id",
                "createdAt",
                "updatedAt",
                "consolidated_symbol",
                "consolidated_description",
            ]
        )
        + "\n"
        + "\t".join(
            [
                "IFXGene:ABC1234",
                "2025-07-08T12:25:29.336300",
                "2026-03-15T09:53:28.798010",
                "MT-RNR1",
                "16S ribosomal RNA",
            ]
        )
        + "\n"
    )

    rows = list(TargetGraphGeneParser(file_path=str(file_path)).all_rows())

    assert rows == [
        {
            "ncats_gene_id": "IFXGene:ABC1234",
            "createdAt": "2025-07-08T12:25:29.336300",
            "updatedAt": "2026-03-15T09:53:28.798010",
            "consolidated_symbol": "MT-RNR1",
            "consolidated_description": "16S ribosomal RNA",
        }
    ]


def test_targetgraph_transcript_parser_handles_missing_optional_fields(tmp_path):
    row = {
        "ncats_transcript_id": "IFXTranscript:XYZ9876",
        "ensembl_transcript_id": "ENST00000387314",
        "ensembl_transcript_id_version": "ENST00000387314.1",
        "ensembl_gene_id": "ENSG00000210049",
        "ensembl_transcript_type": "Mt_tRNA",
        "ensembl_trans_length": "71",
        "ensembl_transcript_tsl": "tslNA",
        "ensembl_canonical": "1.0",
        "ensembl_refseq_MANEselect": "",
        "refseq_rna_id": "",
        "refseq_ncbi_id": "",
        "refseq_status": "",
        "createdAt": "2025-07-08T12:49:11.462992",
        "updatedAt": "2025-07-08T12:49:11.462992",
    }

    file_path = tmp_path / "transcript_ids.tsv"
    file_path.write_text("ncats_transcript_id\nIFXTranscript:XYZ9876\n")

    parser = TargetGraphTranscriptParser(file_path=str(file_path))

    assert parser.get_transcript_location(row) is None
    equivalent_ids = parser.get_equivalent_ids(row)

    assert [eq.id_str() for eq in equivalent_ids] == ["ENSEMBL:ENST00000387314"]


def test_targetgraph_gene_location_ignores_ambiguous_multi_strand():
    row = {
        "consolidated_location": "1|16|17",
        "ensembl_strand": "-1.0|1.0",
    }

    location = TargetGraphGeneParser.get_gene_location(row)

    assert location is not None
    assert location.location == "1|16|17"
    assert location.chromosome == 1
    assert location.strand is None


def test_targetgraph_mapping_ratio_ignores_ambiguous_multi_value():
    row = {"Total_Mapping_Ratio": "0.0833333333333333|0.25"}

    assert TargetGraphGeneParser.get_mapping_ratio(row) is None


def test_targetgraph_protein_parser_preserves_exact_uniprot_isoform(tmp_path):
    file_path = tmp_path / "protein_ids.tsv"
    file_path.write_text("ncats_protein_id\nIFXProtein:L7LS7LZ\n")
    parser = TargetGraphProteinParser(file_path=str(file_path))
    row = {
        "ncats_protein_id": "IFXProtein:L7LS7LZ",
        "uniprot_id": "P56856",
        "uniprot_isoform": "P56856-2",
        "consolidated_ensembl_protein_id": "ENSP00000340939.4",
        "consolidated_refseq_protein": "NP_001002026",
        "consolidated_symbol": "CLDN18",
        "uniprot_secondaryAccessions": "",
        "uniprot_uniProtkbId": "",
        "combined_protein_name": "Claudin-18",
    }

    equivalent_ids = [eq.id_str() for eq in parser.get_equivalent_ids(row)]

    assert "UniProtKB:P56856-2" in equivalent_ids
    assert "UniProtKB:P56856" in equivalent_ids


class _ReleaseFile:
    def __init__(self, path):
        self.path = path

    def file(self, name):
        assert name == self.path.name
        return self.path

    def version_info(self):
        return None


def test_release_parent_edges_use_resolved_ifx_genes(tmp_path):
    transcript = tmp_path / "transcript_ids.tsv"
    transcript.write_text(
        "ncats_transcript_id\tensembl_gene_id\tparent_ifx_gene_id\tparent_resolution_source\tparent_confidence\tcreatedAt\tupdatedAt\n"
        "IFXTranscript:T1\tENSG_WRONG\tIFXGene:G1\tensembl_gene_parent\thigh\t2026-10-08T12:00:00\t2026-10-08T12:00:00\n"
    )
    edge = next(GeneTranscriptEdgeAdapter(_ReleaseFile(transcript)).get_all())[0]
    assert edge.start_node.id == "IFXGene:G1"
    assert edge.details[0].confidence == "high"
    serialized = object.__new__(ArangoOutputAdapter).clean_dict(edge, convert_dates=True)
    assert serialized["details"] == [
        {"resolution_source": "ensembl_gene_parent", "confidence": "high"}
    ]
    json.dumps(serialized)

    protein = tmp_path / "protein_ids.tsv"
    protein.write_text(
        "ncats_protein_id\tuniprot_id\tuniprot_NCBI_id\tparent_ifx_gene_id\tparent_resolution_source\tparent_confidence\tcreatedAt\tupdatedAt\n"
        "IFXProtein:P1\tP12345\t808|809\tIFXGene:G2\tuniprot_ncbi_gene\thigh\t2026-10-08T12:00:00\t2026-10-08T12:00:00\n"
        "IFXProtein:P2\t\t\tIFXGene:G3\tensembl_gene\tmedium\t2026-10-08T12:00:00\t2026-10-08T12:00:00\n"
        "IFXProtein:P3\tP34567\t810\t\tunresolved\t\t2026-10-08T12:00:00\t2026-10-08T12:00:00\n"
    )
    adapter = GeneProteinEdgeAdapter(_ReleaseFile(protein))
    edges = next(adapter.get_all())
    assert [(edge.start_node.id, edge.end_node.id) for edge in edges] == [
        ("IFXGene:G2", "IFXProtein:P1"),
        ("IFXGene:G3", "IFXProtein:P2"),
    ]
    assert edges[0].details[0].resolution_source == "uniprot_ncbi_gene"


def test_release_blank_transcript_parent_fails_in_adapter_and_resolver(tmp_path):
    transcript = tmp_path / "transcript_ids.tsv"
    transcript.write_text(
        "ncats_transcript_id\tensembl_gene_id\tparent_ifx_gene_id\tparent_resolution_source\tparent_confidence\tcreatedAt\tupdatedAt\n"
        "IFXTranscript:T1\tENSG00000000001\t\tensembl_gene_parent\thigh\t2026-10-08T12:00:00\t2026-10-08T12:00:00\n"
    )
    with pytest.raises(ValueError, match="no parent_ifx_gene_id"):
        next(GeneTranscriptEdgeAdapter(_ReleaseFile(transcript)).get_all())
    resolver = object.__new__(TCRDTargetResolver)
    resolver.transcript_parser = TargetGraphTranscriptParser(str(transcript))
    with pytest.raises(ValueError, match="no parent_ifx_gene_id"):
        resolver.get_transcript_ids()


def test_release_mapping_splits_pipe_delimited_ids(tmp_path):
    mapping = tmp_path / "uniprot_mapping.csv"
    mapping.write_text("uniprot_id,uniprot_xref_ChEMBL,uniprot_ccds_id\nP12345,CHEMBL1|CHEMBL2,CCDS1.1|CCDS2.2\n")
    parser = TargetGraphAddtlProteinIDParser(str(mapping))
    ids = [value.id_str() for value in parser.get_id_list(next(parser.all_rows()))]
    assert len(ids) == 4
    assert all("|" not in value for value in ids)


def test_tcrd_canonical_collapse_keeps_only_canonical_gene_parent(tmp_path):
    protein = tmp_path / "protein_ids.tsv"
    protein.write_text(
        "ncats_protein_id\tuniprot_id\tis_canonical\tcanonical_ifx_id\tparent_ifx_gene_id\tparent_resolution_source\n"
        "IFXProtein:P1\tP12345\tTrue\t\tIFXGene:G1\tuniprot_ncbi_gene\n"
        "IFXProtein:P2\tP12345\tFalse\tIFXProtein:P1\tIFXGene:G2\tensembl_transcript_parent\n"
        "IFXProtein:P3\t\tFalse\t\tIFXGene:G3\tensembl_gene\n"
    )
    resolver = object.__new__(TCRDTargetResolver)
    resolver.protein_parsers = [TargetGraphProteinParser(str(protein))]
    protein_ids, _, gene_map = resolver.get_protein_ids(reviewed_only=False, collapse_to_canonical=True)
    assert set(protein_ids) == {"IFXProtein:P1"}
    assert gene_map["IFXProtein:P1"] == {"IFXGene:G1"}


def test_protein_resolvers_read_mapping_from_the_same_release(tmp_path, monkeypatch):
    for name, header in (
        ("gene_ids.tsv", "ncats_gene_id"),
        ("transcript_ids.tsv", "ncats_transcript_id"),
        ("protein_ids.tsv", "ncats_protein_id"),
        ("uniprot_mapping.csv", "uniprot_id"),
    ):
        (tmp_path / name).write_text(header + "\n")

    class Bundle:
        def file(self, name):
            path = tmp_path / name
            if not path.is_file():
                raise FileNotFoundError(path)
            return path

    monkeypatch.setattr(TargetGraphResolver, "__init__", lambda self, **kwargs: None)
    bundle = Bundle()
    protein = TargetGraphProteinResolver(data_source=bundle)
    tcrd = TCRDTargetResolver(
        gene_data_source=bundle,
        transcript_data_source=bundle,
        protein_data_source=bundle,
    )
    expected = str(tmp_path / "uniprot_mapping.csv")
    assert protein.parsers[0].additional_id_file_path == expected
    assert tcrd.protein_parsers[0].additional_id_file_path == expected

    (tmp_path / "uniprot_mapping.csv").unlink()
    with pytest.raises(FileNotFoundError, match="uniprot_mapping.csv"):
        TargetGraphProteinResolver(data_source=bundle)
