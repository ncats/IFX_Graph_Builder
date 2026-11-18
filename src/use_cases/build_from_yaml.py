from src.core.config import ETL_Config
from src.core.curations import resolve_curation_types
from src.core.etl import ETL
from src.infrastructure.object_storage import (
    AwsAssumeRoleCredentials,
    object_storage_from_credentials,
)
from src.shared.db_credentials import DBCredentials
from src.interfaces.resolver_metadata import resolver_fingerprints_by_type
from src.models.registry_dataset import RegistryDatasetMetadata


def preload_core_model_modules():
    # Preload the model modules involved in ETL inheritance/dataclass composition
    # so YAML-driven adapter imports do not control their first-import order.
    import src.models.node
    import src.models.protein
    import src.models.disease
    import src.models.gwas_trait
    import src.models.panther_class
    import src.models.dto_class
    import src.models.external_link
    import src.models.chebi
    import src.models.metabolite_harmonization


class BuildGraphFromYaml:
    configuration: ETL_Config
    etl: ETL
    yaml_file: str

    def __init__(self, yaml_file: str):
        self.load_yaml(yaml_file)

    def load_yaml(self, yaml_file: str):
        preload_core_model_modules()
        self.yaml_file = yaml_file
        self.configuration = ETL_Config(yaml_file)
        output_adapters = self.configuration.create_output_adapters()
        graph_views = self.configuration.config_dict.get("graph_views", [])
        resolver_fingerprints = resolver_fingerprints_by_type(self.configuration.config_dict.get("resolvers"))
        registry_datasets = _registry_datasets_from_config(self.configuration.config_dict)
        for output_adapter in output_adapters:
            if hasattr(output_adapter, "set_graph_views_metadata"):
                output_adapter.set_graph_views_metadata(graph_views=graph_views, source_yaml=yaml_file)
            if hasattr(output_adapter, "set_resolver_metadata"):
                output_adapter.set_resolver_metadata(
                    resolver_fingerprints_by_type=resolver_fingerprints,
                    source_yaml=yaml_file,
                )
            if hasattr(output_adapter, "set_registry_dataset_metadata"):
                output_adapter.set_registry_dataset_metadata(registry_datasets)
        input_adapters = self.configuration.create_input_adapters()
        post_adapters = self.configuration.create_post_adapters()
        curation_snapshots = self._resolve_curations(output_adapters)
        resolver_map = {}
        for resolver in self.configuration.resolvers.values():
            for type in resolver.types:
                if type in resolver_map:
                    raise Exception(f"It only makes sense to have one resolver for each type. Resolver for {type} already exists", resolver_map[type], resolver)
                resolver_map[type] = resolver

        self.etl = ETL(
            input_adapters=input_adapters,
            post_adapters=post_adapters,
            output_adapters=output_adapters,
            resolver_map=resolver_map,
            curation_snapshots=curation_snapshots,
        )

    def _resolve_curations(self, output_adapters) -> dict:
        config = self.configuration.config_dict.get("curations")
        if not config:
            return {}
        if isinstance(config, list):
            curation_types = config
            credentials_config = self.configuration.config_dict.get("object_storage_credentials")
            allow_missing = False
        elif isinstance(config, dict):
            curation_types = config.get("types") or []
            credentials_config = config.get("credentials")
            allow_missing = bool(config.get("allow_missing", False))
        else:
            raise TypeError("curations must be a list of types or a mapping")
        if not curation_types:
            return {}
        unsupported = [adapter for adapter in output_adapters if not adapter.supports_curations()]
        if unsupported:
            names = ", ".join(type(adapter).__name__ for adapter in unsupported)
            raise RuntimeError(f"Configured curations are not supported by output adapter(s): {names}")
        if not isinstance(credentials_config, dict):
            raise ValueError("curations.credentials must resolve to an object-storage credential mapping")
        credentials = (
            AwsAssumeRoleCredentials.from_yaml(credentials_config)
            if "role_arn" in credentials_config or credentials_config.get("type") == "aws_assume_role"
            else DBCredentials.from_yaml(credentials_config)
        )
        storage = object_storage_from_credentials(credentials)
        return resolve_curation_types(
            storage, curation_types, allow_missing=allow_missing
        )

    def prepare_datastore(self, truncate_tables: bool = True):
        self.etl.create_or_truncate_datastores(truncate_tables=truncate_tables)

    def do_etl(self, do_post_processing = True, clean_edges: bool = True, resume: bool = False,
               rerun_completed: bool = False, run_id: str | None = None):
        self.etl.do_etl(do_post_processing, clean_edges, resume=resume,
                        run_id=run_id or self.yaml_file, rerun_completed=rerun_completed)


def _registry_datasets_from_config(config_node) -> list[dict]:
    datasets_by_identity = {}

    def visit(value, usage: str):
        if isinstance(value, RegistryDatasetMetadata):
            metadata = value.to_metadata()
            identity = (metadata["kind"], metadata["snapshot_id"])
            existing = datasets_by_identity.setdefault(identity, metadata)
            usages = set(existing.get("usages") or [])
            usages.add(usage)
            existing["usages"] = sorted(usages)
            return
        if isinstance(value, dict):
            for entry in value.values():
                visit(entry, usage)
            return
        if isinstance(value, list):
            for entry in value:
                visit(entry, usage)

    for entry in config_node.get("input_adapters", []) or []:
        label = f"adapter:{entry.get('class')}" if isinstance(entry, dict) else "adapter"
        visit(entry, label)
    for entry in config_node.get("post_adapters", []) or []:
        label = f"post_adapter:{entry.get('class')}" if isinstance(entry, dict) else "post_adapter"
        visit(entry, label)
    for entry in config_node.get("resolvers", []) or []:
        label = f"resolver:{entry.get('label') or entry.get('class')}" if isinstance(entry, dict) else "resolver"
        visit(entry, label)

    return sorted(
        datasets_by_identity.values(),
        key=lambda item: (item["kind"], item["snapshot_id"]),
    )
