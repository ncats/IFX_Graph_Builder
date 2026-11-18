from abc import ABC, abstractmethod
from typing import List, Generator

from src.constants import Prefix
from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import AnalytePathwayEdge as SqliteAnalytePathwayEdge, Analyte as SqliteAnalyte, Source as SqliteSource
from src.models.gene import Gene
from src.models.metabolite import Metabolite
from src.models.node import EquivalentId
from src.models.pathway import GenePathwayEdge, MetabolitePathwayEdge, ProteinPathwayEdge, Pathway
from src.models.protein import Protein


class AnalytePathwayEdgeAdapter(RaMPSqliteAdapter, ABC):

    @abstractmethod
    def get_id_prefix(self) -> str:
        pass

    @abstractmethod
    def get_all(self):
        pass


class MetabolitePathwayEdgeAdapter(AnalytePathwayEdgeAdapter):
    def get_id_prefix(self) -> str:
        return "RAMP_C"

    def get_all(self) -> Generator[List[MetabolitePathwayEdge], None, None]:
        results = self.get_session().query(
            SqliteAnalytePathwayEdge.rampId,
            SqliteAnalytePathwayEdge.pathwayRampId,
            SqliteAnalytePathwayEdge.pathwaySource
        ).filter(SqliteAnalytePathwayEdge.rampId.startswith(self.get_id_prefix())).all()

        yield [
            MetabolitePathwayEdge(
                start_node=Metabolite(id=row[0]),
                end_node=Pathway(id=row[1]),
                source=row[2]
            ) for row in results
        ]


class ProteinPathwayEdgeAdapter(AnalytePathwayEdgeAdapter):
    def get_id_prefix(self) -> str:
        return "RAMP_G"

    def get_all(self) -> Generator[List[ProteinPathwayEdge], None, None]:
        results = self.get_session().query(
            SqliteAnalytePathwayEdge.rampId,
            SqliteAnalytePathwayEdge.pathwayRampId,
            SqliteAnalytePathwayEdge.pathwaySource
        ).filter(SqliteAnalytePathwayEdge.rampId.startswith(self.get_id_prefix())).all()

        yield [
            ProteinPathwayEdge(
                start_node=Protein(id=row[0]),
                end_node=Pathway(id=row[1]),
                source=row[2]
            ) for row in results
        ]


class EnsemblGenePathwayEdgeAdapter(AnalytePathwayEdgeAdapter):
    def get_id_prefix(self) -> str:
        return "RAMP_G"

    def get_all(self) -> Generator[List[GenePathwayEdge], None, None]:
        results = (
            self.get_session()
            .query(
                SqliteAnalytePathwayEdge.pathwayRampId,
                SqliteSource.sourceId,
                SqliteAnalytePathwayEdge.pathwaySource
            )
            .join(
                SqliteAnalyte,
                SqliteAnalyte.rampId == SqliteAnalytePathwayEdge.rampId
            )
            .join(
                SqliteSource,
                SqliteSource.rampId == SqliteAnalyte.rampId
            )
            .filter(
                SqliteAnalyte.type == "gene",
                SqliteSource.IDtype == "ensembl"
            )
            .order_by(SqliteAnalytePathwayEdge.pathwayRampId)
            .all()
        )

        relationships: List[GenePathwayEdge] = []
        for row in results:
            gene_id = row[1].split(":")[1]
            ensembl_gene_id = EquivalentId(id=gene_id, type=Prefix.ENSEMBL).id_str()
            relationships.append(
                GenePathwayEdge(
                    start_node=Gene(id=ensembl_gene_id),
                    end_node=Pathway(id=row[0]),
                    source=row[2]
                )
            )
        yield relationships