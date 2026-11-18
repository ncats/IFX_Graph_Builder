from dataclasses import dataclass
from typing import List, Generator

from src.constants import DataSourceName
from src.input_adapters.sqlite_ramp.metabolite_adapter import MetaboliteAdapter
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.metabolite import Metabolite
from src.models.node import Relationship

# This adapter reads from two RaMP sqlite databases and defines relationships between metabolite nodes that share input IDs
# It's used for comparing metabolites between versions of RaMP or changes in RaMP entity resolution code

@dataclass
class SharesMetaboliteIDrelationship(Relationship):
    count: int = 0


class MetaboliteSetEdgeAdapter(InputAdapter):
    name: str = "Metabolite Set Relationship Adapter"

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.Dummy

    def get_version(self) -> DatasourceVersionInfo:
        return DatasourceVersionInfo(
            version="calculated"
        )
    left_adapter: MetaboliteAdapter
    right_adapter: MetaboliteAdapter

    def get_all(self) -> Generator[List[SharesMetaboliteIDrelationship], None, None]:
        left_mets = list(self.left_adapter.get_all())[0]
        right_mets = list(self.right_adapter.get_all())[0]

        reverse_lookup = {}

        for metabolite in left_mets:
            ramp_id = metabolite.id
            equiv_ids = [equiv_id.id_str() for equiv_id in metabolite.xref]
            for id in equiv_ids:
                if id in reverse_lookup:
                    existing_id = reverse_lookup[id]
                    if ramp_id != existing_id:
                        raise Exception(f"this wasn't supposed to happen: both {ramp_id} and {existing_id} mention {id}")
                else:
                    reverse_lookup[id] = ramp_id

        relationship_map = {}
        for metabolite in right_mets:
            ramp_id = metabolite.id
            equiv_ids = [equiv_id.id_str() for equiv_id in metabolite.xref]
            for id in equiv_ids:
                if id not in reverse_lookup:
                    print(f"{id} exists in the 'right' list and not the 'left'. Skipping it...")
                    continue
                else:
                    other_ramp_id = reverse_lookup[id]
                    composite_id = f"{other_ramp_id}|{ramp_id}"
                    if composite_id in relationship_map:
                        val = relationship_map[composite_id]
                        relationship_map[composite_id] = val + 1
                    else:
                        relationship_map[composite_id] = 1

        relationships = []
        for key, count in relationship_map.items():
            start_id, end_id = key.split('|')
            start_met = Metabolite(id=start_id)
            end_met = Metabolite(id=end_id)

            relationships.append(
                SharesMetaboliteIDrelationship(
                    start_node=start_met,
                    end_node=end_met,
                    count=count))

        yield relationships

    def set_left(self, data_source):
        self.left_adapter = MetaboliteAdapter(data_source=data_source)
        return self

    def set_right(self, data_source):
        self.right_adapter = MetaboliteAdapter(data_source=data_source)
        return self
