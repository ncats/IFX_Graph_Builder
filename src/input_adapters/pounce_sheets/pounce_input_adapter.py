import os
from datetime import datetime
from typing import Generator, List, Union

from src.constants import DataSourceName
from src.input_adapters.pounce_sheets.pounce_node_builder import PounceNodeBuilder
from src.input_adapters.pounce_sheets.pounce_parser import PounceParser
from src.input_adapters.pounce_sheets.validator_loader import load_validators
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.node import Node, Relationship


class PounceInputAdapter(InputAdapter):
    def __init__(self, project_file: str, experiment_files: List[str] = None,
                 stats_results_files: List[str | None] = None, validators_config: str = None):
        self.project_file = project_file
        self.experiment_files = experiment_files or []
        self.stats_results_files = stats_results_files or []
        self.validators_config = validators_config
        self._pounce_parser = PounceParser()
        self._cached_validation_result = None

    # --- Interface methods ---

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.NCATSPounce

    def get_version(self) -> DatasourceVersionInfo:
        data = self.get_validation_data()
        return DatasourceVersionInfo(
            version_date=data.project.date,
            download_date=datetime.fromtimestamp(os.path.getmtime(self.project_file)).date()
        )

    def get_validators(self) -> list:
        if self.validators_config:
            return load_validators(self.validators_config)
        return []

    def _get_validation_result(self):
        if self._cached_validation_result is None:
            self._cached_validation_result = self._pounce_parser.parse_all(
                self.project_file, self.experiment_files, self.stats_results_files
            )
        return self._cached_validation_result

    def get_validation_data(self):
        data, _ = self._get_validation_result()
        return data

    def get_structural_issues(self):
        _, issues = self._get_validation_result()
        return issues

    def get_all(self) -> Generator[List[Union[Node, Relationship]], None, None]:
        builder = PounceNodeBuilder()
        yield from builder.build(
            self.project_file,
            self.experiment_files,
            self.stats_results_files,
            self._pounce_parser,
        )
