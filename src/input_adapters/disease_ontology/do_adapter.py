"""Disease Ontology adapters over the Registry's doid.obo snapshot."""

import re
from typing import Generator, List, Optional

import obonet

from src.constants import DataSourceName
from src.input_adapters.flat_file_adapter import FlatFileAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.disease import Disease, DODiseaseParentEdge


_QUOTED_VALUE = re.compile(r'^"((?:\\.|[^"\\])*)"')
_SYNONYM_SCOPE = re.compile(r'^"(?:\\.|[^"\\])*"\s+(BROAD|EXACT|NARROW|RELATED)\b')
_SYNONYM_SCOPE_ORDER = {"BROAD": 0, "EXACT": 1, "NARROW": 2, "RELATED": 3}


def _obo_quoted_value(value: str) -> str:
    """Extract a definition or synonym, decoding OBO escaped characters."""
    match = _QUOTED_VALUE.match(value)
    if match is None:
        raise ValueError(f"Expected a quoted OBO value: {value[:80]}")
    return re.sub(
        r"\\(.)",
        lambda escaped: "\n" if escaped.group(1) == "n" else escaped.group(1),
        match.group(1),
    )


def _obo_synonym(value: str) -> tuple[int, str]:
    match = _SYNONYM_SCOPE.match(value)
    if match is None:
        raise ValueError(f"Expected an OBO synonym scope: {value[:80]}")
    return _SYNONYM_SCOPE_ORDER[match.group(1)], _obo_quoted_value(value)


class DOBaseAdapter(FlatFileAdapter):
    """Shared OBO reader for Disease Ontology graph and Pharos table adapters."""

    version_info: Optional[DatasourceVersionInfo] = None

    def __init__(self, data_source):
        super().__init__(file_path=str(data_source.file("doid.obo")))
        self.version_info = data_source.version_info()
        self._graph = None

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.DiseaseOntology

    def get_version(self) -> DatasourceVersionInfo:
        return self.version_info

    def _get_graph(self):
        if self._graph is None:
            graph = obonet.read_obo(self.file_path)
            if not any(term_id.startswith("DOID:") for term_id in graph.nodes):
                raise ValueError(f"No DOID terms found in {self.file_path}")
            self._graph = graph
        return self._graph

    def _iter_included_do_nodes(self):
        for term_id, term_data in self._get_graph().nodes(data=True):
            if term_id.startswith("DOID:"):
                yield term_id, term_data

    def _included_do_ids(self) -> set[str]:
        return {term_id for term_id, _ in self._iter_included_do_nodes()}

    def _iter_parent_pairs(self):
        included_do_ids = self._included_do_ids()
        for child_id, parent_id, predicate in self._get_graph().edges(keys=True):
            if (
                predicate == "is_a"
                and child_id in included_do_ids
                and parent_id in included_do_ids
            ):
                yield child_id, parent_id


class DODiseaseAdapter(DOBaseAdapter):

    @staticmethod
    def _to_disease(term_id: str, term_data: dict) -> Disease:
        definition = term_data.get("def")
        comment = term_data.get("comment")
        return Disease(
            id=term_id,
            name=term_data.get("name"),
            type="CLASS",
            do_description=_obo_quoted_value(definition) if definition else None,
            synonyms=[
                text for _, text in sorted(
                    _obo_synonym(value) for value in term_data.get("synonym", [])
                )
            ],
            comments=[comment] if comment else [],
        )

    def get_all(self) -> Generator[List[Disease], None, None]:
        yield [
            self._to_disease(term_id, term_data)
            for term_id, term_data in self._iter_included_do_nodes()
        ]


class DODiseaseParentEdgeAdapter(DOBaseAdapter):

    def get_all(self) -> Generator[List[DODiseaseParentEdge], None, None]:
        yield [
            DODiseaseParentEdge(
                start_node=Disease(id=child_id),
                end_node=Disease(id=parent_id),
                source="DO",
            )
            for child_id, parent_id in self._iter_parent_pairs()
        ]
