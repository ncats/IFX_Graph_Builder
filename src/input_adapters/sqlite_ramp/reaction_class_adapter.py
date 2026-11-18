from typing import List, Generator

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import ReactionClass as SqliteReactionClass
from src.models.reaction import ReactionClass


class ReactionClassAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[ReactionClass], None, None]:
        results = self.get_session().query(
            SqliteReactionClass.rxn_class_ec,
            SqliteReactionClass.ec_level,
            SqliteReactionClass.rxn_class
        ).distinct()

        reaction_classes: List[ReactionClass] = [
            ReactionClass(
                id=row[0],
                level=row[1],
                name=row[2]
            ) for row in results
        ]
        yield reaction_classes

