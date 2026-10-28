import csv
import gzip
import json
from pathlib import Path

import pytest
import yaml

from src.input_adapters.chembl.drug_node import DrugNodeAdapter, ProteinDrugEdgeAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.registry_dataset import RegistryDatasetKind
from src.shared.db_credentials import DBCredentials


class SnapshotSource:
    kind = RegistryDatasetKind.SOURCE
    snapshot_id = "chembl:activity_exports:chembl37-export1"

    def __init__(self, directory: Path):
        self.directory = directory

    def file(self, name: str) -> Path:
        path = self.directory / name
        if not path.is_file():
            raise FileNotFoundError(path)
        return path

    def version_info(self) -> DatasourceVersionInfo:
        return DatasourceVersionInfo(version="chembl37-export1")


class ExternalSource:
    kind = RegistryDatasetKind.EXTERNAL

    def version_info(self) -> DatasourceVersionInfo:
        return DatasourceVersionInfo(version="chembl36")


def write_snapshot(
    directory: Path, *, columns=None, value="5.40",
    compound_name="Drug", canonical_smiles="CCO",
) -> SnapshotSource:
    (directory / "chembl_activity_export_contract.json").write_text(
        json.dumps({
            "snapshot_id": SnapshotSource.snapshot_id,
            "pchembl_cutoff": 5,
            "row_count": 1,
            "null_encoding": "backslash_N_with_leading_backslash_escape",
        }),
        encoding="utf-8",
    )
    row = {
        "activity_id": "101", "compound_id": "CHEMBL1", "compound_name": compound_name,
        "canonical_smiles": canonical_smiles, "target_id": "CHEMBL_T1",
        "uniprot_id": "P12345", "action_type": "INHIBITOR",
        "pchembl_value": value, "standard_type": "IC50", "assay_type": "B",
        "pubmed_id": "12345", "src_id": "1",
    }
    columns = columns or list(row)
    with gzip.open(
        directory / "chembl_activities_pchembl5.tsv.gz",
        "wt", encoding="utf-8", newline="",
    ) as handle:
        writer = csv.DictWriter(handle, fieldnames=columns, delimiter="\t")
        writer.writeheader()
        writer.writerow({key: row[key] for key in columns})
    return SnapshotSource(directory)


def test_snapshot_drug_and_activity_edges_share_the_export(tmp_path: Path) -> None:
    source = write_snapshot(tmp_path)
    drugs = [item for batch in DrugNodeAdapter(data_source=source).get_all() for item in batch]
    edges = [item for batch in ProteinDrugEdgeAdapter(data_source=source).get_all() for item in batch]
    assert len(drugs) == len(edges) == 1
    assert drugs[0].id == "CHEMBL.COMPOUND:CHEMBL1"
    assert edges[0].start_node.id == "UniProtKB:P12345"
    detail = edges[0].details[0]
    assert detail.ref_id == 101
    assert detail.act_value == 5.4
    assert detail.act_pmids == [12345]


def test_snapshot_applies_higher_cutoff_and_rejects_lower_cutoff(tmp_path: Path) -> None:
    source = write_snapshot(tmp_path)
    assert len(list(DrugNodeAdapter(data_source=source, pchembl_cutoff=5.2).get_all())[0]) == 1
    assert list(DrugNodeAdapter(data_source=source, pchembl_cutoff=6).get_all()) == [[]]
    source = write_snapshot(tmp_path, value="5.20")
    assert len(list(DrugNodeAdapter(data_source=source, pchembl_cutoff=5.2).get_all())[0]) == 1
    with pytest.raises(ValueError, match="contains only pChEMBL"):
        DrugNodeAdapter(data_source=source, pchembl_cutoff=4)


def test_snapshot_rejects_missing_column_and_below_cutoff(tmp_path: Path) -> None:
    source = write_snapshot(tmp_path, columns=["activity_id", "compound_id"])
    with pytest.raises(ValueError, match="lacks required columns"):
        list(DrugNodeAdapter(data_source=source).get_all())
    source = write_snapshot(tmp_path, value="4.99")
    with pytest.raises(ValueError, match="below export cutoff"):
        list(DrugNodeAdapter(data_source=source).get_all())


def test_historical_external_input_still_accepts_database_credentials() -> None:
    credentials = DBCredentials(
        url="localhost", user="reader", password=None, schema="chembl36"
    )
    adapter = DrugNodeAdapter(data_source=ExternalSource(), credentials=credentials)
    assert adapter.version_info.version == "chembl36"
    assert adapter.credentials.schema == "chembl36"


def test_snapshot_preserves_empty_string_and_sql_null(tmp_path: Path) -> None:
    source = write_snapshot(tmp_path, compound_name="", canonical_smiles="\\N")
    drug = next(iter(DrugNodeAdapter(data_source=source).get_all()))[0]
    assert drug.name == ""
    assert drug.smiles is None
    source = write_snapshot(tmp_path, compound_name="\\\\N")
    drug = next(iter(DrugNodeAdapter(data_source=source).get_all()))[0]
    assert drug.name == "\\N"


@pytest.mark.parametrize("yaml_name", [
    "target_graph.yaml", "pharos.yaml", "impatient_pharos.yaml",
])
def test_active_pharos_chembl_configs_use_pinned_export(yaml_name: str) -> None:
    path = Path("src/use_cases/pharos") / yaml_name
    config = yaml.safe_load(path.read_text(encoding="utf-8"))
    entries = [entry for entry in config["input_adapters"]
               if entry.get("import") == "./src/input_adapters/chembl/drug_node.py"]
    assert len(entries) == 2
    for entry in entries:
        assert entry["kwargs"]["data_source"] == {
            "kind": "source_snapshot",
            "snapshot_id": "chembl:activity_exports:chembl37-export1",
        }
        assert "credentials" not in entry
