import re
from typing import Generator, List

from src.input_adapters.disease_ontology.do_adapter import DOBaseAdapter, DODiseaseAdapter
from src.input_adapters.mondo.mondo_adapter import MondoBaseAdapter, MondoDiseaseAdapter
from src.models.tcrd_disease_ontology import MondoTerm, MondoTermParentEdge, DOTerm, DOTermParentEdge


class MondoTableAdapter(MondoDiseaseAdapter):

    _IDENTIFIERS_ORG_PREFIX_MAP = {
        "mesh": "MESH",
        "medgen": "MEDGEN",
        "snomedct": "SNOMEDCT",
        "omim": "OMIM",
        "orphanet": "ORPHANET",
        "icd10cm": "ICD10CM",
        "icd10": "ICD10",
        "icd9": "ICD9",
        "efo": "EFO",
        "doid": "DOID",
        "ncit": "NCIT",
        "umls": "UMLS",
    }

    @classmethod
    def _normalize_xref_value(cls, value: str) -> str | None:
        if not value:
            return None
        raw = value.strip()
        if not raw:
            return None
        if "://" not in raw and ":" in raw:
            return raw

        obo_match = re.search(r"/obo/([A-Za-z0-9]+)_([^/#]+)$", raw)
        if obo_match:
            return f"{obo_match.group(1).upper()}:{obo_match.group(2)}"

        identifiers_match = re.search(r"identifiers\.org/([^/]+)/([^/?#]+)$", raw)
        if identifiers_match:
            prefix = identifiers_match.group(1).lower()
            mapped_prefix = cls._IDENTIFIERS_ORG_PREFIX_MAP.get(prefix, prefix.upper())
            return f"{mapped_prefix}:{identifiers_match.group(2)}"

        umls_match = re.search(r"/umls/id/([^/?#]+)$", raw, re.IGNORECASE)
        if umls_match:
            return f"UMLS:{umls_match.group(1)}"

        return None

    @classmethod
    def _extract_exact_matches(cls, node: dict) -> list[str]:
        exact_matches = set()
        for entry in (node.get("meta") or {}).get("basicPropertyValues") or []:
            predicate = (entry.get("pred") or "").lower()
            if not predicate.endswith("exactmatch"):
                continue
            normalized = cls._normalize_xref_value(entry.get("val") or "")
            if normalized:
                exact_matches.add(normalized)
        return sorted(exact_matches)

    @classmethod
    def _extract_xrefs(cls, node: dict) -> list[str]:
        xrefs = set()
        for entry in (node.get("meta") or {}).get("xrefs") or []:
            normalized = cls._normalize_xref_value(entry.get("val") if isinstance(entry, dict) else entry)
            if normalized:
                xrefs.add(normalized)
        return sorted(xrefs)

    def _to_term(self, node: dict, subset_label_map: dict) -> MondoTerm:
        disease = self._to_disease(node, subset_label_map)
        comments = disease.comments or []
        return MondoTerm(
            id=disease.id,
            name=disease.name,
            mondo_description=disease.mondo_description,
            comment=comments[0] if comments else None,
            mondo_xrefs=self._extract_xrefs(node),
            exact_matches=self._extract_exact_matches(node),
        )

    def get_all(self) -> Generator[List[MondoTerm], None, None]:
        subset_label_map = self._build_subset_label_map()
        yield [self._to_term(node, subset_label_map) for node in self._iter_included_disease_nodes()]


class MondoTableParentEdgeAdapter(MondoBaseAdapter):

    def get_all(self) -> Generator[List[MondoTermParentEdge], None, None]:
        included_disease_ids = self._included_disease_ids()
        edges: List[MondoTermParentEdge] = []

        for edge in self._get_graph().get("edges") or []:
            if edge.get("pred") != "is_a":
                continue

            child_id = self._normalize_mondo_id(edge.get("sub"))
            parent_id = self._normalize_mondo_id(edge.get("obj"))
            if child_id is None or parent_id is None:
                continue
            if child_id not in included_disease_ids or parent_id not in included_disease_ids:
                continue

            edges.append(MondoTermParentEdge(
                start_node=MondoTerm(id=child_id),
                end_node=MondoTerm(id=parent_id),
            ))

        yield edges


class DOTableAdapter(DODiseaseAdapter):

    def _to_term(self, term_id: str, term_data: dict) -> DOTerm:
        disease = self._to_disease(term_id, term_data)
        return DOTerm(
            id=disease.id,
            name=disease.name,
            do_description=disease.do_description,
        )

    def get_all(self) -> Generator[List[DOTerm], None, None]:
        yield [
            self._to_term(term_id, term_data)
            for term_id, term_data in self._iter_included_do_nodes()
        ]


class DOTableParentEdgeAdapter(DOBaseAdapter):

    def get_all(self) -> Generator[List[DOTermParentEdge], None, None]:
        yield [
            DOTermParentEdge(
                start_node=DOTerm(id=child_id),
                end_node=DOTerm(id=parent_id),
            )
            for child_id, parent_id in self._iter_parent_pairs()
        ]
