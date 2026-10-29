import csv
import hashlib

import pytest

from src.input_adapters.disease_harmonizer.associations import (
    ASSERTION_MEMBER, SUMMARY_MEMBER, DiseaseConsumerAssociationAdapter,
)
from src.interfaces.input_adapter import InputAdapter
from src.models.disease import Disease, GeneDiseaseEdge, ProteinDiseaseEdge
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.protein import Protein


class Source:
    snapshot_id = "ifx_harmonizers:diseases:2.3.1"

    def __init__(self, files):
        self.files = files
        self.manifest = {"files": [
            {"path": name, "sha256": hashlib.sha256(path.read_bytes()).hexdigest()}
            for name, path in files.items()
        ]}

    def file(self, name):
        return self.files[name]

    def version_info(self):
        return DatasourceVersionInfo(version="2.3.1", version_date=None, download_date=None)


def _write(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]), delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)
    return path


def _fixture(tmp_path, *, duplicate_link=False):
    assertions = [
        dict(association_id="IFXDA:ONE", ncats_disease_id="IFXDisease:D", association_source="ClinGen",
             evidence_origins="ClinGen", negated="True", subject_type="gene", ncats_gene_id="IFXGene:G",
             ncats_protein_id="", subject_id_raw="HGNC:1", disease_id_raw="MONDO:1",
             disease_mapping_basis="exact", predicate="biolink:gene_associated_with_condition",
             evidence_type="gene_disease_validity", evidence_label="Disputed", evidence_score="",
             publications="https://example.org/assertion", provided_by="ClinGen"),
        dict(association_id="IFXDA:TWO", ncats_disease_id="IFXDisease:D", association_source="Jensen",
             evidence_origins="JensenLab", negated="False", subject_type="protein", ncats_gene_id="",
             ncats_protein_id="IFXProtein:P", subject_id_raw="ENSP0001", disease_id_raw="DOID:1",
             disease_mapping_basis="exact", predicate="biolink:related_to",
             evidence_type="jensen_textmining", evidence_label="z-score", evidence_score="7.2",
             publications="", provided_by="JensenLab"),
    ]
    summary = [dict(ncats_disease_id="IFXDisease:D", ncats_gene_id="IFXGene:G",
                    support_tier="conflicting", support_score_raw="2.5", n_supporting_sources="2",
                    n_independent_origins="2", evidence_origins="ClinGen|JensenLab",
                    association_ids="IFXDA:ONE|IFXDA:TWO")]
    if duplicate_link:
        summary[0]["association_ids"] += "|IFXDA:ONE"
    files = {
        ASSERTION_MEMBER: _write(tmp_path / "assertions.tsv", assertions),
        SUMMARY_MEMBER: _write(tmp_path / "summary.tsv", summary),
    }
    return Source(files)


def test_grouped_assertions_preserve_raw_ids_and_promotable_payload(tmp_path):
    edge = next(DiseaseConsumerAssociationAdapter(_fixture(tmp_path)).get_all())[0]
    assert type(edge) is GeneDiseaseEdge
    assert (edge.start_node.id, edge.end_node.id) == ("IFXGene:G", "IFXDisease:D")
    assert len(edge.details) == 2
    assert edge.disease_harmonizer_support["support_tier"] == "conflicting"
    assert edge.gene_assessments == [{
        "gene_id": "IFXGene:G", "tier": "conflicting", "score": "2.5",
        "origins": ["ClinGen", "JensenLab"],
    }]
    assert "association_ids" not in edge.disease_harmonizer_support
    assert (edge.details[0].source_id, edge.details[0].negated) == ("MONDO:1", True)
    assert edge.details[0].evidence_label == "Disputed"
    assert edge.details[0].evidence_terms == []
    assert edge.details[0].url == "https://example.org/assertion"
    assert edge.details[0].assertion_payload["publications"] == "https://example.org/assertion"
    assert edge.details[1].source_protein_id == "IFXProtein:P"
    assert edge.details[1].source_subject_id == "ENSP0001"
    assert edge.details[1].evidence_label == "z-score"
    assert edge.details[1].evidence_terms == []
    assert edge.details[1].zscore == 7.2
    assert "association_id" not in edge.details[0].assertion_payload
    promoted = InputAdapter._canonicalize_relationship_class(
        edge, Protein(id="IFXProtein:CANONICAL"), Disease(id="IFXDisease:D")
    )
    assert type(promoted) is ProteinDiseaseEdge
    assert promoted.details == edge.details
    assert promoted.gene_assessments == edge.gene_assessments
    assert promoted.disease_harmonizer_support is None


def test_association_adapter_emits_gene_endpoints_before_edges(tmp_path):
    batches = list(DiseaseConsumerAssociationAdapter(
        _fixture(tmp_path), emit_gene_nodes=True
    ).get_all())
    assert [type(item).__name__ for item in batches[0]] == ["Gene"]
    assert batches[0][0].id == "IFXGene:G"
    assert type(batches[1][0]) is GeneDiseaseEdge


def test_duplicate_assertion_link_fails(tmp_path):
    adapter = DiseaseConsumerAssociationAdapter(_fixture(tmp_path, duplicate_link=True))
    with pytest.raises(ValueError, match="missing/repeated assertion"):
        list(adapter.get_all())
