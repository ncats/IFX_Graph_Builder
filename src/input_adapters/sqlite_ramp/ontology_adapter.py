from typing import List, Generator

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import Ontology as SqliteOntology
from src.models.ontology import Ontology


class OntologyAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[Ontology], None, None]:
        results = self.get_session().query(
            SqliteOntology.rampOntologyId,
            SqliteOntology.commonName,
            SqliteOntology.HMDBOntologyType
        ).all()

        ontologies: List[Ontology] = [
            Ontology(
                id=row[0],
                commonName=row[1],
                HMDBOntologyType=row[2]
            ) for row in results
        ]
        yield ontologies

