from typing import List, Generator

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import Reaction as SqliteReaction
from src.models.reaction import Reaction, ReactionDirection


class ReactionAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[Reaction], None, None]:
        results = self.get_session().query(
            SqliteReaction.ramp_rxn_id,
            SqliteReaction.rxn_source_id,
            SqliteReaction.is_transport,
            SqliteReaction.direction,
            SqliteReaction.label,
            SqliteReaction.equation,
            SqliteReaction.html_equation
        ).all()

        reactions: List[Reaction] = [
            Reaction(
                id=row[0],
                source_id=row[1],
                is_transport=row[2],
                direction=ReactionDirection.parse(row[3]),
                label=row[4],
                equation=row[5],
                html_equation=row[6]
            ) for row in results
        ]
        yield reactions

