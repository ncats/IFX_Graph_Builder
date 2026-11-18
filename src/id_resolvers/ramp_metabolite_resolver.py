from typing import List, Dict

from src.input_adapters.sql_adapter import SqliteAdapter
from src.input_adapters.sqlite_ramp.tables import Source as SqliteSource
from src.input_adapters.sqlite_ramp.registry_sqlite import ramp_sqlite_file
from src.interfaces.id_resolver import IdResolver, IdMatch
from src.models.node import Node


class RampMetaboliteIdResolver(IdResolver, SqliteAdapter):
    source_to_ramp: Dict[str, str]
    ramp_to_sources: Dict[str, List[str]]

    @staticmethod
    def _normalize_id(input_id: str) -> str:
        if isinstance(input_id, str) and ':' in input_id:
            prefix, id_part = input_id.split(':', 1)
            return f"{prefix.lower()}:{id_part}"
        return input_id

    def __init__(self, data_source, **kwargs):
        self.data_source = data_source
        SqliteAdapter.__init__(self, str(ramp_sqlite_file(data_source)))
        IdResolver.__init__(self, **kwargs)
        self._build_lookup()

    def _build_lookup(self):
        self.source_to_ramp = {}
        self.ramp_to_sources = {}

        results = self.get_session().query(
            SqliteSource.sourceId,
            SqliteSource.rampId
        ).filter(SqliteSource.geneOrCompound == 'compound').all()

        for source_id, ramp_id in results:
            normalized = self._normalize_id(source_id)
            self.source_to_ramp[normalized] = ramp_id
            if ramp_id not in self.ramp_to_sources:
                self.ramp_to_sources[ramp_id] = []
            if source_id not in self.ramp_to_sources[ramp_id]:
                self.ramp_to_sources[ramp_id].append(source_id)

        # also map ramp IDs to themselves so already-canonical IDs pass through
        for ramp_id in list(self.ramp_to_sources.keys()):
            self.source_to_ramp[self._normalize_id(ramp_id)] = ramp_id

        print(f"RampMetaboliteIdResolver: loaded {len(self.source_to_ramp)} source IDs "
              f"for {len(self.ramp_to_sources)} metabolites")

    def resolve_internal(self, input_nodes: List[Node]) -> Dict[str, List[IdMatch]]:
        result_list = {}
        for node in input_nodes:
            ramp_id = self.source_to_ramp.get(self._normalize_id(node.id))
            if ramp_id is not None:
                match_type = 'exact' if node.id == ramp_id else 'sourceId'
                result_list[node.id] = [
                    IdMatch(
                        input=node.id,
                        match=ramp_id,
                        equivalent_ids=self.ramp_to_sources.get(ramp_id, []),
                        context=[match_type]
                    )
                ]
            else:
                result_list[node.id] = []
        return result_list
