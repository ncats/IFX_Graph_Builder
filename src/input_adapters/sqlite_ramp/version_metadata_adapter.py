from typing import List, Generator

from sqlalchemy import desc
from src.input_adapters.sqlite_ramp.ramp_sqlite_adapter import RaMPSqliteAdapter
from src.input_adapters.sqlite_ramp.tables import DBVersion as SqliteDBVersion, VersionInfo as SqliteVersionInfo
from src.models.version import DatabaseVersion, DataVersion, DatabaseDataVersionEdge


class VersionMetaAdapter(RaMPSqliteAdapter):

    def get_all(self) -> Generator[List[DatabaseVersion], None, None]:
        result = self.get_session().query(
            SqliteDBVersion.ramp_version,
            SqliteDBVersion.load_timestamp,
            SqliteDBVersion.version_notes
        ).order_by(desc(SqliteDBVersion.ramp_version)).first()

        db_version = DatabaseVersion(
            id = result[0],
            timestamp = result[1],
            notes= result[2]
        )

        results = self.get_session().query(
            SqliteVersionInfo.data_source_id,
            SqliteVersionInfo.data_source_name,
            SqliteVersionInfo.data_source_url,
            SqliteVersionInfo.data_source_version
        ).filter(SqliteVersionInfo.status == 'current').all()

        data_versions: List[DataVersion] = [
            DataVersion(
                id=row[0],
                name=row[1],
                url=row[2],
                version=row[3]
            ) for row in results
        ]

        relationships: List[DatabaseDataVersionEdge] = [
            DatabaseDataVersionEdge(
                start_node=db_version,
                end_node=data_version
            ) for data_version in data_versions
        ]
        yield [db_version, *data_versions, *relationships]

