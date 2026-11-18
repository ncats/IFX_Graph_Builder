from typing import List, Generator

from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import AnalyteHasOntology as SqliteAnalyteHasOntology
from src.models.metabolite import Metabolite
from src.models.ontology import MetaboliteOntologyEdge, Ontology


class MetaboliteOntologyEdgeAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[MetaboliteOntologyEdge], None, None]:
        results = self.get_session().query(
            SqliteAnalyteHasOntology.rampCompoundId,
            SqliteAnalyteHasOntology.rampOntologyId
        ).all()

        analyte_ontology_relationships: List[MetaboliteOntologyEdge] = [
            MetaboliteOntologyEdge(
                start_node=Metabolite(id=row[0]),
                end_node=Ontology(id=row[1])
            ) for row in results
        ]
        yield analyte_ontology_relationships

