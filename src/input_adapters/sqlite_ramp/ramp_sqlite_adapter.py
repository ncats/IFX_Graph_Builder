from abc import ABC

from src.constants import DataSourceName
from src.input_adapters.sql_adapter import SqliteAdapter
from src.input_adapters.sqlite_ramp.registry_sqlite import ramp_sqlite_file
from src.input_adapters.sqlite_ramp.ramp_version_info import RaMPVersionInfo
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.version import DataVersion, DatabaseVersion


class RaMPSqliteAdapter(SqliteAdapter, InputAdapter, ABC):
    ramp_version_info: RaMPVersionInfo

    def __init__(self, data_source):
        self.data_source = data_source
        super().__init__(str(ramp_sqlite_file(data_source)))
        self.ramp_version_info = RaMPVersionInfo()
        self.ramp_version_info.initialize(self.get_session())

    def get_data_version(self, datakey: str) -> DataVersion:
        requested_data = [data for data in self.ramp_version_info.data_versions if data.id == datakey]
        if len(requested_data) == 0:
            raise Exception(f"datakey not found: {datakey}")
        return requested_data[0]

    def get_database_version(self) -> DatabaseVersion:
        return self.ramp_version_info.db_version

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.RaMP

    def get_version(self) -> DatasourceVersionInfo:
        return self.data_source.version_info()
