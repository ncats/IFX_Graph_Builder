"""Protein–ligand activity assertions from the DrugCentral export."""

from collections import OrderedDict
from collections.abc import Generator

from src.constants import DataSourceName, Prefix
from src.input_adapters.drug_central.drug_node import DrugCentralAdapter
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.ligand import ActivityDetails, Ligand, ProteinLigandEdge
from src.models.node import EquivalentId
from src.models.protein import Protein


_FILE = "drugcentral_drug_target_interaction.tsv.gz"
_PROJECTED_COLUMNS = (
    "STRUCT_ID",
    "ACCESSION",
    "ACT_VALUE",
    "ACT_TYPE",
    "ACT_COMMENT",
    "ACT_SOURCE",
    "MOA_SOURCE",
    "MOA",
    "ACTION_TYPE",
    "ACT_PMID",
    "MOA_PMID",
)
_REQUIRED_COLUMNS = (*_PROJECTED_COLUMNS, "ORGANISM")


class ProteinDrugEdgeAdapter(InputAdapter, DrugCentralAdapter):
    batch_size = 10000

    def __init__(self, data_source):
        DrugCentralAdapter.__init__(self, data_source)

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.DrugCentral

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def get_all(self) -> Generator[list[ProteinLigandEdge], None, None]:
        edges: OrderedDict[str, ProteinLigandEdge] = OrderedDict()
        seen_rows: set[tuple[str, ...]] = set()

        for row in self.rows(_FILE, _REQUIRED_COLUMNS):
            if row["ORGANISM"] != "Homo sapiens":
                continue
            struct_id = row["STRUCT_ID"].strip()
            accession_text = row["ACCESSION"].strip()
            if not struct_id or not accession_text:
                raise ValueError(
                    f"DrugCentral interaction export has blank STRUCT_ID or ACCESSION"
                )
            projected = tuple(row[column] for column in _PROJECTED_COLUMNS)
            if projected in seen_rows:
                continue
            seen_rows.add(projected)

            detail = ActivityDetails(activity_source=DataSourceName.DrugCentral)
            for source_column, destination in (
                ("ACT_TYPE", "act_type"),
                ("ACTION_TYPE", "action_type"),
                ("ACT_SOURCE", "act_source"),
                ("MOA_SOURCE", "moa_source"),
                ("ACT_COMMENT", "comment"),
            ):
                if row[source_column]:
                    setattr(detail, destination, row[source_column])
            if row["ACT_VALUE"]:
                detail.act_value = float(row["ACT_VALUE"])
            if row["ACT_PMID"]:
                detail.act_pmids = [int(row["ACT_PMID"])]
            if row["MOA_PMID"]:
                detail.moa_pmid = int(row["MOA_PMID"])
            if row["MOA"] == "1":
                detail.has_moa = True

            ligand = Ligand(
                id=EquivalentId(id=struct_id, type=Prefix.DrugCentral).id_str()
            )
            for accession in (part.strip() for part in accession_text.split("|")):
                if not accession:
                    raise ValueError(
                        f"DrugCentral interaction has an empty accession token for {struct_id}"
                    )
                protein = Protein(
                    id=EquivalentId(id=accession, type=Prefix.UniProtKB).id_str()
                )
                edge_id = f"{protein.id}|{ligand.id}"
                if edge_id not in edges:
                    edges[edge_id] = ProteinLigandEdge(
                        start_node=protein, end_node=ligand, details=[]
                    )
                edges[edge_id].details.append(detail)

        edge_values = list(edges.values())
        for index in range(0, len(edge_values), self.batch_size):
            yield edge_values[index:index + self.batch_size]
