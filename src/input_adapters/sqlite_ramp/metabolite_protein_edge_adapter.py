from typing import List, Generator

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import Catalyzed as SqliteCatalyzed
from src.models.protein import Protein
from src.models.metabolite import MetaboliteProteinEdge, Metabolite


class MetaboliteProteinEdgeAdapter(RaMPSqliteAdapter):
    cached_audit_trail_info = None

    def get_all(self) -> Generator[List[MetaboliteProteinEdge], None, None]:
        results = self.get_session().query(
            SqliteCatalyzed.rampCompoundId,
            SqliteCatalyzed.rampGeneId
        ).all()

        relationships: List[MetaboliteProteinEdge] = [
            MetaboliteProteinEdge(
                start_node=Metabolite(id=row[0]),
                end_node=Protein(id=row[1])
            ) for row in results
        ]
        yield relationships
