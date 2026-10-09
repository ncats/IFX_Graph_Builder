import csv
import gzip
from pathlib import Path

import pytest
import yaml

from src.input_adapters.drug_central.drug_indication import DrugCentralIndicationAdapter
from src.input_adapters.drug_central.drug_node import DrugNodeAdapter
from src.input_adapters.drug_central.protein_drug_edge import ProteinDrugEdgeAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.disease import Disease, ProteinDiseaseEdge


class FileSource:
    def __init__(self, directory: Path):
        self.directory = directory

    def file(self, name: str) -> Path:
        path = self.directory / name
        if not path.is_file():
            raise FileNotFoundError(path)
        return path

    def version_info(self) -> DatasourceVersionInfo:
        return DatasourceVersionInfo(version="56-export1")


def write_tsv(path: Path, columns: list[str], rows: list[dict[str, str]]) -> None:
    opener = gzip.open if path.suffix == ".gz" else open
    with opener(path, "wt", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerows(rows)


def test_drugcentral_dump_preserves_nodes_activities_and_indication_join(
    tmp_path: Path,
) -> None:
    write_tsv(
        tmp_path / "drugcentral_structures.smiles.tsv",
        ["ID", "DRUG_NAME", "SMILES"],
        [
            {"ID": "1", "DRUG_NAME": "Drug One", "SMILES": "CCO"},
            {"ID": "2", "DRUG_NAME": "Drug Two", "SMILES": ""},
        ],
    )
    activity_columns = [
        "ACTIVITY_ID", "STRUCT_ID", "ACCESSION", "ACT_VALUE", "ACT_TYPE",
        "ACT_COMMENT", "ACT_SOURCE", "MOA_SOURCE", "MOA", "ACTION_TYPE",
        "ACT_PMID", "MOA_PMID", "ORGANISM",
    ]
    first = {
        "ACTIVITY_ID": "10", "STRUCT_ID": "1", "ACCESSION": "P12345|Q99999",
        "ACT_VALUE": "4.5", "ACT_TYPE": "IC50", "ACT_COMMENT": "assay",
        "ACT_SOURCE": "CHEMBL", "MOA_SOURCE": "DRUG LABEL", "MOA": "1",
        "ACTION_TYPE": "INHIBITOR", "ACT_PMID": "123", "MOA_PMID": "456",
        "ORGANISM": "Homo sapiens",
    }
    write_tsv(
        tmp_path / "drugcentral_drug_target_interaction.tsv.gz",
        activity_columns,
        [
            first,
            first | {"ACTIVITY_ID": "11"},  # SQL projection DISTINCT removes this.
            first | {"ACTIVITY_ID": "12", "STRUCT_ID": "2", "ACCESSION": "P12345", "ACT_VALUE": ""},
            first | {"ACTIVITY_ID": "13", "ORGANISM": "Mus musculus"},
        ],
    )
    write_tsv(
        tmp_path / "drugcentral_indications.tsv",
        [
            "STRUCT_ID", "RELATIONSHIP_NAME", "CONCEPT_NAME", "UMLS_CUI",
            "SNOMED_CONCEPTID", "DOID",
        ],
        [
            {"STRUCT_ID": "1", "RELATIONSHIP_NAME": "indication", "CONCEPT_NAME": "Disease One", "UMLS_CUI": "C123", "SNOMED_CONCEPTID": "42", "DOID": "DOID:1"},
            {"STRUCT_ID": "1", "RELATIONSHIP_NAME": "indication", "CONCEPT_NAME": "Disease One", "UMLS_CUI": "C123", "SNOMED_CONCEPTID": "42", "DOID": "DOID:1"},
            {"STRUCT_ID": "2", "RELATIONSHIP_NAME": "indication", "CONCEPT_NAME": "Disease Two", "UMLS_CUI": "", "SNOMED_CONCEPTID": "", "DOID": ""},
        ],
    )
    source = FileSource(tmp_path)

    drugs = [item for batch in DrugNodeAdapter(source).get_all() for item in batch]
    activities = [item for batch in ProteinDrugEdgeAdapter(source).get_all() for item in batch]
    indications = [item for batch in DrugCentralIndicationAdapter(source).get_all() for item in batch]

    assert [item.id for item in drugs] == ["DrugCentral:1", "DrugCentral:2"]
    assert drugs[1].smiles is None
    assert len(activities) == 3
    assert sum(len(item.details) for item in activities) == 3
    first_detail = activities[0].details[0]
    assert first_detail.act_value == 4.5
    assert first_detail.act_pmids == [123]
    assert first_detail.moa_pmid == 456
    assert first_detail.has_moa is True
    assert first_detail.act_source == "CHEMBL"
    assert len([item for item in indications if isinstance(item, Disease)]) == 2
    indication_edges = [item for item in indications if isinstance(item, ProteinDiseaseEdge)]
    assert len(indication_edges) == 3
    assert {item.start_node.id for item in indication_edges} == {
        "UniProtKB:P12345", "UniProtKB:Q99999"
    }
    assert any(item.end_node.id == "UMLS:C123" for item in indication_edges)
    assert any(item.end_node.id.startswith("DrugCentral:INDICATION:") for item in indication_edges)
    assert sum(len(item.details) for item in indication_edges) == 3


def test_drugcentral_dump_missing_required_header_fails(tmp_path: Path) -> None:
    write_tsv(
        tmp_path / "drugcentral_structures.smiles.tsv",
        ["ID", "DRUG_NAME"],
        [{"ID": "1", "DRUG_NAME": "Drug One"}],
    )

    with pytest.raises(ValueError, match="lacks required columns: SMILES"):
        list(DrugNodeAdapter(FileSource(tmp_path)).get_all())


@pytest.mark.parametrize(
    "yaml_name,expected_adapters",
    [
        ("target_graph.yaml", 3),
        ("pharos.yaml", 3),
        ("impatient_pharos.yaml", 2),
        ("pharos_current_tdls_reconstruction.yaml", 2),
    ],
)
def test_pharos_drugcentral_configs_use_pinned_dump(
    yaml_name: str, expected_adapters: int
) -> None:
    path = Path("src/use_cases/pharos") / yaml_name
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    entries = [
        entry
        for entry in config["input_adapters"]
        if "drug_central/" in entry.get("import", "")
    ]
    assert len(entries) == expected_adapters
    for entry in entries:
        assert entry["kwargs"]["data_source"] == {
            "kind": "source_snapshot",
            "snapshot_id": "drugcentral:drug_exports:56-export1",
        }
        assert "credentials" not in entry
