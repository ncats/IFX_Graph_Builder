from typing import Generator, List

from src.constants import DataSourceName
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.pounce.experiment import PlatformType
from src.models.pounce.project import ProjectType
from src.models.pounce.stats_summary import AnalyteSummary
from src.shared.arango_adapter import ArangoAdapter
from src.shared.db_credentials import DBCredentials
from src.shared.record_merger import FieldConflictBehavior


class AnalyteSummaryAdapter(InputAdapter, ArangoAdapter):
    field_conflict_behavior = FieldConflictBehavior.KeepLast

    def __init__(self, credentials: DBCredentials, database_name: str):
        ArangoAdapter.__init__(self, credentials=credentials, database_name=database_name)

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.NCATSPounce

    def get_version(self) -> DatasourceVersionInfo:
        return DatasourceVersionInfo()

    def get_all(self) -> Generator[List[AnalyteSummary], None, None]:
        rows = self.runQuery(_analyte_summary_query)
        nodes = []
        for row in rows:
            platform_type = PlatformType.parse(row['platform_type'])
            project_type = ProjectType.parse(row['project_type'])
            if platform_type is None or project_type is None:
                continue
            nodes.append(AnalyteSummary(
                id=f"analyte_summary-{row['platform_type']}-{row['project_type']}",
                platform_type=platform_type,
                project_type=project_type,
                analyte_count=int(row['analyte_count']),
            ))
        yield nodes


_analyte_summary_query = """
FOR e IN `Experiment`
    FILTER e.calculated_properties != null
       AND e.calculated_properties.analyte_count != null
       AND e.calculated_properties.analyte_count > 0
       AND e.platform_type != null
    FOR p IN INBOUND e `ProjectExperimentEdge`
        FILTER p.project_type != null AND LENGTH(p.project_type) > 0
        FOR pt IN p.project_type
            COLLECT platform_type = e.platform_type, project_type = pt
            AGGREGATE analyte_count = SUM(TO_NUMBER(e.calculated_properties.analyte_count))
            RETURN {platform_type, project_type, analyte_count}
"""
