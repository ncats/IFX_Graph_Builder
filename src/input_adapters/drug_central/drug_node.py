"""DrugCentral source adapters backed by the pinned Registry export."""

import csv
import gzip
from collections.abc import Generator, Iterable
from pathlib import Path

from src.constants import DataSourceName, Prefix
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.ligand import Ligand
from src.models.node import EquivalentId


class DrugCentralAdapter:
    def __init__(self, data_source):
        self.data_source = data_source
        self.version_info = data_source.version_info()

    def rows(
        self, file_name: str, required_columns: Iterable[str]
    ) -> Generator[dict[str, str], None, None]:
        path: Path = self.data_source.file(file_name)
        opener = gzip.open if path.suffix == ".gz" else open
        with opener(path, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            columns = set(reader.fieldnames or ())
            missing = set(required_columns) - columns
            if missing:
                raise ValueError(
                    f"DrugCentral export {path} lacks required columns: "
                    + ", ".join(sorted(missing))
                )
            for line_number, row in enumerate(reader, start=2):
                if None in row or any(value is None for value in row.values()):
                    raise ValueError(
                        f"Malformed DrugCentral export row in {path} at line {line_number}"
                    )
                yield row


class DrugNodeAdapter(InputAdapter, DrugCentralAdapter):
    def __init__(self, data_source):
        DrugCentralAdapter.__init__(self, data_source)

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.DrugCentral

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[list[Ligand], None, None]:
        drugs: list[Ligand] = []
        seen_ids: set[str] = set()
        for row in self.rows(
            "drugcentral_structures.smiles.tsv", ("ID", "DRUG_NAME", "SMILES")
        ):
            struct_id = row["ID"].strip()
            if not struct_id:
                raise ValueError("DrugCentral structures export contains a blank ID")
            drug_id = EquivalentId(id=struct_id, type=Prefix.DrugCentral).id_str()
            if drug_id in seen_ids:
                continue
            seen_ids.add(drug_id)
            drugs.append(
                Ligand(
                    id=drug_id,
                    name=row["DRUG_NAME"] or None,
                    smiles=row["SMILES"] or None,
                    isDrug=True,
                )
            )
            if len(drugs) >= self.batch_size:
                yield drugs
                drugs = []
        if drugs:
            yield drugs
