from typing import Any, List

import yaml

from src.core.object_factory import create_object_from_config
from src.core.registry_integration import RegistryIntegration, RegistryReference
from src.interfaces.id_resolver import IdResolver
from src.interfaces.input_adapter import InputAdapter
from src.interfaces.output_adapter import OutputAdapter


def _resolve_registry_datasets(
    config_node: Any,
    registry: RegistryIntegration,
    resolved_datasets: dict,
):
    if isinstance(config_node, dict):
        resolved = {}
        for key, value in config_node.items():
            if key == "data_source" or key.endswith("_data_source"):
                reference = RegistryReference.parse(value)
                if reference not in resolved_datasets:
                    resolved_datasets[reference] = registry.resolve(value)
                resolved[key] = resolved_datasets[reference]
            else:
                resolved[key] = _resolve_registry_datasets(value, registry, resolved_datasets)
        return resolved
    if isinstance(config_node, list):
        return [_resolve_registry_datasets(entry, registry, resolved_datasets) for entry in config_node]
    return config_node


def resolve_registry_references(config_dict: dict) -> dict:
    registry_config = config_dict.get("registry")
    if not registry_config:
        return config_dict
    if not isinstance(registry_config, dict):
        raise TypeError("registry configuration must be a mapping")
    registry = RegistryIntegration.connect(registry_config)
    resolved = dict(config_dict)
    resolved_datasets = {}
    for key in ("resolvers", "input_adapters", "post_adapters"):
        if key in resolved:
            resolved[key] = _resolve_registry_datasets(resolved[key], registry, resolved_datasets)
    return resolved

class Config:
    config_dict: {}
    yaml_files: list = []
    def __init__(self, yaml_file):
        self.config_dict = self.load_config_from_yaml(yaml_file)

    def load_config_from_yaml(self, file_path):
        config_dict = self.load_one_yaml(file_path)
        config_dict = self._load_nested_yamls(config_dict)
        config_dict = resolve_registry_references(config_dict)
        print('Configuration loaded from yaml file(s)')
        return config_dict

    def load_one_yaml(self, yaml_file):
        with open(yaml_file, "r") as file:
            config_dict = yaml.safe_load(file)
            self.yaml_files.append(yaml_file)
            return config_dict

    def _load_nested_yamls(self, config_node, parent_key=None):
        if isinstance(config_node, str) and (config_node.endswith(".yaml") or config_node.endswith(".yml")):
            nested_config = self.load_one_yaml(config_node)
            return self._load_nested_yamls(nested_config)
        if isinstance(config_node, dict):
            for key, value in config_node.items():
                if isinstance(value, list):
                    for index, entry in enumerate(value):
                        value[index] = self._load_nested_yamls(entry, parent_key=key)
                if isinstance(value, dict):
                    config_node[key] = self._load_nested_yamls(value, parent_key=key)
                if (
                    isinstance(value, str)
                    and (value.endswith(".yaml") or value.endswith(".yml"))
                    and not (parent_key == "registry" and key == "credentials")
                ):
                    nested_config = self.load_one_yaml(value)
                    config_node[key] = self._load_nested_yamls(nested_config)
        if isinstance(config_node, list):
            if isinstance(config_node, list):
                for index, entry in enumerate(config_node):
                    config_node[index] = self._load_nested_yamls(entry, parent_key=parent_key)
        return config_node

    def __repr__(self):
        return f"{self.__class__.__name__}({self.config_dict})"

    def create_object_list(self, key, required = True):
        objects = []
        if key not in self.config_dict and required == False:
            return objects
        config = self.config_dict[key]
        for c in config:
            obj = create_object_from_config(c)
            objects.append(obj)
        return objects


class ETL_Config(Config):
    resolvers: dict[str, IdResolver] = {}

    def __init__(self, yaml_file):
        self.resolvers = {}
        Config.__init__(self, yaml_file)
        self.create_resolvers()

    def create_resolvers(self):
        if 'resolvers' in self.config_dict:
            for config in self.config_dict['resolvers']:
                key = config['label']
                obj = create_object_from_config(config)
                self.resolvers[key] = obj

    def create_output_adapters(self) -> List[OutputAdapter]:
        return self.create_object_list('output_adapters')

    def create_input_adapters(self) -> List[InputAdapter]:
        return self.create_object_list('input_adapters', False)

    def create_post_adapters(self) -> List[InputAdapter]:
        """Create graph-dependent adapters that run after primary ingestion."""
        return self.create_object_list('post_adapters', False)
