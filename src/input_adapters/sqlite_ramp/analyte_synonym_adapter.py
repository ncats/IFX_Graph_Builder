from abc import ABC, abstractmethod
from typing import List, Generator
from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import AnalyteSynonym as SqliteAnalyteSynonym
from src.models.analyte import Synonym, Analyte
from src.models.metabolite import Metabolite
from src.models.protein import Protein


class AnalyteSynonymAdapter(RaMPSqliteAdapter, ABC):
    cls = Analyte

    @abstractmethod
    def get_id_prefix(self) -> str:
        pass

    def get_all(self) -> Generator[List[Analyte], None, None]:
        results = self.get_session().query(
            SqliteAnalyteSynonym.rampId,
            SqliteAnalyteSynonym.Synonym,
            SqliteAnalyteSynonym.source
        ).filter(SqliteAnalyteSynonym.rampId.startswith(self.get_id_prefix())).distinct()

        analyte_dict = {}
        for row in results:
            ramp_id = row[0]
            if ramp_id in analyte_dict:
                analyte_dict[ramp_id].append(row)
            else:
                analyte_dict[ramp_id] = [row]

        objs: List[Analyte] = [
            self.cls(
                id=key,
                synonyms=[
                    Synonym(
                         term=row[1], source=row[2]) for row in synonym_list]
            ) for key, synonym_list in analyte_dict.items()
        ]
        yield objs


class MetaboliteSynonymAdapter(AnalyteSynonymAdapter):
    cls = Metabolite

    def get_id_prefix(self) -> str:
        return "RAMP_C"


class ProteinSynonymAdapter(AnalyteSynonymAdapter):
    cls = Protein

    def get_id_prefix(self) -> str:
        return "RAMP_G"
