from typing import List, Generator

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import ReactionToMetabolite as SqliteReactionToMetabolite
from src.models.metabolite import Metabolite, MetaboliteReactionEdge
from src.models.reaction import Reaction


class MetaboliteReactionEdgeAdapter(RaMPSqliteAdapter):
    def get_all(self) -> Generator[List[MetaboliteReactionEdge], None, None]:
        results = self.get_session().query(
            SqliteReactionToMetabolite.ramp_cmpd_id,
            SqliteReactionToMetabolite.ramp_rxn_id,
            SqliteReactionToMetabolite.substrate_product,
            SqliteReactionToMetabolite.is_cofactor
        ).all()

        relationships: List[MetaboliteReactionEdge] = [
            MetaboliteReactionEdge(
                start_node=Metabolite(id=row[0]),
                end_node=Reaction(id=row[1]),
                substrate_product=row[2],
                is_cofactor=row[3] == 1
            ) for row in results
        ]
        yield relationships

