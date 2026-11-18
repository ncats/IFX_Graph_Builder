from typing import Generator, List

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import Pathway as SqlitePathway
from src.models.pathway import Pathway


class PathwayAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[Pathway], None, None]:
        results = self.get_session().query(
            SqlitePathway.pathwayRampId,
            SqlitePathway.sourceId,
            SqlitePathway.type,
            SqlitePathway.pathwayCategory,
            SqlitePathway.pathwayName
        ).all()

        pathways: List[Pathway] = [
            Pathway(
                id=row[0],
                source_id=row[1],
                type=row[2],
                category=row[3],
                name=row[4]
            ) for row in results
        ]
        yield pathways

