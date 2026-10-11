import csv
import gzip
from pathlib import Path
from typing import Generator, List, Optional

from src.constants import DataSourceName, Prefix
from src.input_adapters.flat_file_adapter import FlatFileAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.ppi import PPIEdge
from src.models.protein import Protein


class StringPPIAdapter(FlatFileAdapter):
    version_info: DatasourceVersionInfo

    def __init__(
        self,
        data_source,
        score_cutoff: int = 400,
        max_rows: Optional[int] = None,
    ):
        file_path = str(data_source.file())
        FlatFileAdapter.__init__(self, file_path=file_path)
        self.score_cutoff = int(score_cutoff)
        self.max_rows = max_rows
        self.version_info = data_source.version_info()

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.STRING

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _open_input(self):
        path = Path(self.file_path)
        if path.suffix == ".gz":
            return gzip.open(path, "rt", encoding="utf-8")
        return open(path, "r", encoding="utf-8")

    @staticmethod
    def _strip_taxon_prefix(protein_id: str) -> str:
        if protein_id.startswith("9606."):
            return protein_id.split(".", 1)[1]
        return protein_id

    def get_all(self) -> Generator[List[PPIEdge], None, None]:
        batch: List[PPIEdge] = []
        kept_rows = 0
        with self._open_input() as handle:
            header = handle.readline().strip().split()
            for line in handle:
                if self.max_rows is not None and kept_rows >= self.max_rows:
                    break
                parts = line.strip().split()
                if not parts:
                    continue
                row = dict(zip(header, parts))
                score = int(row["combined_score"])
                if score < self.score_cutoff:
                    continue

                protein1 = self._strip_taxon_prefix(row["protein1"])
                protein2 = self._strip_taxon_prefix(row["protein2"])
                if protein1 == protein2:
                    continue
                protein1, protein2 = sorted((protein1, protein2))

                edge = PPIEdge(
                    start_node=Protein(id=f"{Prefix.ENSEMBL}:{protein1}"),
                    end_node=Protein(id=f"{Prefix.ENSEMBL}:{protein2}"),
                    score=[score],
                )
                batch.append(edge)
                kept_rows += 1
                if len(batch) >= self.batch_size:
                    yield batch
                    batch = []
        yield batch
