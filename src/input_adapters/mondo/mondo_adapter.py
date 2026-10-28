import json
import re
from typing import Generator, List, Optional

from src.constants import DataSourceName
from src.input_adapters.flat_file_adapter import FlatFileAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.disease import Disease, MondoDiseaseParentEdge


class MondoBaseAdapter(FlatFileAdapter):
    version_info: Optional[DatasourceVersionInfo] = None

    def __init__(self, data_source):
        file_path = str(data_source.file())
        FlatFileAdapter.__init__(self, file_path=file_path)
        self.version_info = data_source.version_info()
        self._graph = None

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.Mondo

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _get_graph(self) -> dict:
        if self._graph is None:
            with open(self.file_path, "r", encoding="utf-8") as handle:
                payload = json.load(handle)
            graphs = payload.get("graphs") or []
            if not graphs:
                raise ValueError(f"No graphs found in {self.file_path}")
            self._graph = graphs[0]
        return self._graph

    @staticmethod
    def _normalize_mondo_id(value: Optional[str]) -> Optional[str]:
        if not value:
            return None
        if value.startswith("MONDO:"):
            return value
        match = re.search(r"/obo/MONDO_(\d+)$", value)
        if match:
            return f"MONDO:{match.group(1)}"
        return None

    @staticmethod
    def _extract_values(entries) -> List[str]:
        values = []
        for entry in entries or []:
            value = entry.get("val") if isinstance(entry, dict) else entry
            if isinstance(value, str) and value.strip():
                values.append(value.strip())
        return values

    @staticmethod
    def _is_class_node(node: dict) -> bool:
        return node.get("type") == "CLASS"

    @staticmethod
    def _is_deprecated(node: dict) -> bool:
        return (node.get("meta") or {}).get("deprecated") is True

    def _iter_included_disease_nodes(self):
        for node in self._get_graph().get("nodes") or []:
            if not self._is_class_node(node):
                continue
            if self._is_deprecated(node):
                continue
            disease_id = self._normalize_mondo_id(node.get("id"))
            if disease_id is None:
                continue
            yield node

    def _included_disease_ids(self) -> set:
        return {
            self._normalize_mondo_id(node.get("id"))
            for node in self._iter_included_disease_nodes()
        }


class MondoDiseaseAdapter(MondoBaseAdapter):

    def _build_subset_label_map(self) -> dict:
        subset_label_map = {}
        for node in self._get_graph().get("nodes") or []:
            if node.get("type") != "PROPERTY":
                continue
            subset_id = node.get("id")
            if not subset_id:
                continue
            comments = self._extract_values((node.get("meta") or {}).get("comments"))
            if comments:
                subset_label_map[subset_id] = comments[0]
        return subset_label_map

    @staticmethod
    def _subset_fallback_label(subset_uri: str) -> str:
        if "#" in subset_uri:
            return subset_uri.split("#", 1)[1]
        return subset_uri.rstrip("/").rsplit("/", 1)[-1]

    def _resolve_subset_labels(self, subset_values: List[str], subset_label_map: dict) -> List[str]:
        return [
            subset_label_map.get(subset_value, self._subset_fallback_label(subset_value))
            for subset_value in subset_values
        ]

    def _to_disease(self, node: dict, subset_label_map: dict) -> Disease:
        meta = node.get("meta") or {}
        definition = meta.get("definition") or {}
        subset_values = self._extract_values(meta.get("subsets"))

        return Disease(
            id=self._normalize_mondo_id(node.get("id")),
            name=node.get("lbl"),
            type=node.get("type"),
            mondo_description=definition.get("val") if isinstance(definition, dict) else None,
            subsets=self._resolve_subset_labels(subset_values, subset_label_map),
            synonyms=self._extract_values(meta.get("synonyms")),
            comments=self._extract_values(meta.get("comments")),
        )

    def get_all(self) -> Generator[List[Disease], None, None]:
        subset_label_map = self._build_subset_label_map()
        diseases = [self._to_disease(node, subset_label_map) for node in self._iter_included_disease_nodes()]
        yield diseases


class MondoDiseaseParentEdgeAdapter(MondoBaseAdapter):

    def _edge_disease_ids(self, edge: dict) -> tuple[Optional[str], Optional[str]]:
        child_id = self._normalize_mondo_id(edge.get("sub"))
        parent_id = self._normalize_mondo_id(edge.get("obj"))
        return child_id, parent_id

    @staticmethod
    def _is_parent_edge(edge: dict) -> bool:
        return edge.get("pred") == "is_a"

    @staticmethod
    def _to_parent_edge(child_id: str, parent_id: str) -> MondoDiseaseParentEdge:
        return MondoDiseaseParentEdge(
            start_node=Disease(id=child_id),
            end_node=Disease(id=parent_id),
            source="MONDO",
        )

    def get_all(self) -> Generator[List[MondoDiseaseParentEdge], None, None]:
        included_disease_ids = self._included_disease_ids()
        edges: List[MondoDiseaseParentEdge] = []

        for edge in self._get_graph().get("edges") or []:
            if not self._is_parent_edge(edge):
                continue

            child_id, parent_id = self._edge_disease_ids(edge)
            if child_id is None or parent_id is None:
                continue
            if child_id not in included_disease_ids or parent_id not in included_disease_ids:
                continue

            edges.append(self._to_parent_edge(child_id, parent_id))

        yield edges
