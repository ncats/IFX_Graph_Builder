from typing import List, Generator
from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.models.metabolite_class import MetaboliteClass
from src.input_adapters.sqlite_ramp.tables import MetaboliteClass as SqliteMetaboliteClass


class MetaboliteClassAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[MetaboliteClass], None, None]:
        results = self.get_session().query(
            SqliteMetaboliteClass.class_level_name,
            SqliteMetaboliteClass.class_name,
            SqliteMetaboliteClass.source
        ).distinct().all()

        metabolite_classes: List[MetaboliteClass] = [
            MetaboliteClass(
                id=MetaboliteClass.compiled_name(row[0], row[1]),
                level=row[0],
                name=row[1],
                source=row[2]) for row in
            results
        ]
        yield metabolite_classes
