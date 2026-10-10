import csv
import re
from typing import Generator, List

from src.constants import DataSourceName, Prefix
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.node import EquivalentId
from src.models.protein import Protein


class AntibodyCountAdapter(InputAdapter):
    batch_size: int = 1000

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.Antibodypedia

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def __init__(self, data_source, identifier_mode: str = "uniprot"):
        super().__init__()
        if identifier_mode not in {"uniprot", "ifx_protein"}:
            raise ValueError(f"Unsupported Antibodypedia identifier mode: {identifier_mode}")
        self.file_path = str(data_source.file())
        self.version_info = data_source.version_info()
        self.identifier_mode = identifier_mode

    def get_all(self) -> Generator[List[Protein], None, None]:
        proteins = []
        with open(self.file_path, mode="r", newline="") as file:
            csv_reader = csv.DictReader(file)
            id_column = "ncats_protein_id" if self.identifier_mode == "ifx_protein" else "uniprot_id"
            required_columns = {id_column, "antibodies"}
            missing = required_columns - set(csv_reader.fieldnames or ())
            if missing:
                raise ValueError(f"Antibodypedia input {self.file_path} lacks columns: {sorted(missing)}")
            for row_number, row in enumerate(csv_reader, start=2):
                antibody_str = row["antibodies"]
                if not antibody_str:
                    continue
                if self.identifier_mode == "ifx_protein":
                    match = re.fullmatch(r"(\d+) antibodies(?: from \d+ providers?)?", antibody_str.strip())
                else:
                    match = re.search(r"\d+", antibody_str)
                if not match:
                    raise ValueError(
                        f"Invalid Antibodypedia count at {self.file_path}:{row_number}: {antibody_str!r}"
                    )
                antibody_count = int(match.group(1) if self.identifier_mode == "ifx_protein" else match.group())
                if antibody_count <= 0:
                    continue
                if self.identifier_mode == "ifx_protein":
                    protein_id = row[id_column].strip()
                    if not re.fullmatch(r"IFXProtein:[A-Za-z0-9]+", protein_id):
                        raise ValueError(
                            f"Invalid Antibodypedia protein ID at {self.file_path}:{row_number}: {protein_id!r}"
                        )
                else:
                    protein_id = EquivalentId(id=row[id_column], type=Prefix.UniProtKB).id_str()
                proteins.append(Protein(id=protein_id, antibody_count=[antibody_count]))

        yield proteins
