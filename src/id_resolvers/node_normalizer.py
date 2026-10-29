import time
from typing import List, Dict, Iterable

import requests
from src.interfaces.id_resolver import IdResolver, IdMatch
from src.models.node import Node
from src.shared.util import yield_per


class TranslatorNodeNormResolver(IdResolver):
    base_url = "https://nodenormalization-sri.renci.org/1.5"
    name = "Translator Node Normalizer"
    conflate_genes_and_proteins: bool = False
    retryable_status_codes = {408, 429, 500, 502, 503, 504}
    node_type_to_biolink_type = {
        "Condition": "biolink:Disease",
        "Disease": "biolink:Disease",
        "Drug": "biolink:Drug",
        "Gene": "biolink:Gene",
        "Ligand": "biolink:SmallMolecule",
        "Metabolite": "biolink:SmallMolecule",
        "Phenotype": "biolink:PhenotypicFeature",
        "PhenotypicFeature": "biolink:PhenotypicFeature",
        "Protein": "biolink:Protein",
        "Transcript": "biolink:Transcript",
    }
    example_ids_by_biolink_type = {
        "biolink:Disease": ["MONDO:0005148", "DOID:9352"],
        "biolink:Drug": ["CHEBI:15365", "DRUGBANK:DB00945"],
        "biolink:Gene": ["NCBIGene:7157", "HGNC:11998"],
        "biolink:PhenotypicFeature": ["HP:0001250", "HP:0004322"],
        "biolink:Protein": ["UniProtKB:P04637", "UniProtKB:P38398"],
        "biolink:SmallMolecule": ["CHEBI:15377", "PUBCHEM.COMPOUND:2244"],
        "biolink:Transcript": ["ENSEMBL:ENST00000269305"],
    }

    def __init__(self,
                 batch_size: int = 50000,
                 request_timeout: int = 120,
                 max_retries: int = 10,
                 retry_backoff_seconds: int = 60,
                 base_url: str | None = None,
                 **kwargs):
        super().__init__(**kwargs)
        if base_url is not None:
            if not isinstance(base_url, str) or not base_url.startswith("https://"):
                raise ValueError("Node Normalizer base_url must be an HTTPS URL")
            self.base_url = base_url.rstrip("/")
        self.name = f"{self.name}: {self.base_url}"
        self.batch_size = batch_size
        self.request_timeout = request_timeout
        self.max_retries = max(1, max_retries)
        self.retry_backoff_seconds = retry_backoff_seconds

    def node_norm_url(self):
        return f"{self.base_url}/get_normalized_nodes"

    def node_norm_prefixes_url(self):
        return f"{self.base_url}/get_curie_prefixes"

    def get_prefix_counts(self) -> List[Dict[str, object]]:
        semantic_types = self._semantic_types_for_registered_types()
        rows = []
        if semantic_types:
            for semantic_type in semantic_types:
                response = requests.get(
                    self.node_norm_prefixes_url(),
                    params={"semantic_type": semantic_type},
                    timeout=self.request_timeout,
                )
                response.raise_for_status()
                rows.extend(self._parse_prefix_counts_payload(response.json()))
        else:
            response = requests.get(
                self.node_norm_prefixes_url(),
                timeout=self.request_timeout,
            )
            response.raise_for_status()
            rows = self._parse_prefix_counts_payload(response.json())
        return self._merge_prefix_counts(rows)

    def _semantic_types_for_registered_types(self) -> List[str]:
        return sorted({
            self.node_type_to_biolink_type[node_type]
            for node_type in self.types
            if node_type in self.node_type_to_biolink_type
        })

    def get_example_ids(self, limit: int = 5) -> List[str]:
        examples = []
        for semantic_type in self._semantic_types_for_registered_types():
            for example_id in self.example_ids_by_biolink_type.get(semantic_type, []):
                if example_id not in examples:
                    examples.append(example_id)
                if len(examples) >= limit:
                    return examples
        return examples

    @staticmethod
    def _merge_prefix_counts(rows: List[Dict[str, object]]) -> List[Dict[str, object]]:
        merged = {}
        unknown_count_prefixes = set()
        for row in rows:
            prefix = str(row.get("prefix") or "").strip()
            if not prefix:
                continue
            count = row.get("count")
            if isinstance(count, int):
                merged[prefix] = merged.get(prefix, 0) + count
            else:
                merged.setdefault(prefix, 0)
                unknown_count_prefixes.add(prefix)
        return sorted(
            [
                {
                    "prefix": prefix,
                    "count": None if prefix in unknown_count_prefixes else count,
                }
                for prefix, count in merged.items()
            ],
            key=lambda row: (
                -(row["count"] if isinstance(row.get("count"), int) else -1),
                row["prefix"].lower(),
            ),
        )

    @staticmethod
    def _parse_prefix_counts_payload(payload) -> List[Dict[str, object]]:
        if isinstance(payload, dict):
            for key in ("curie_prefixes", "prefixes", "data"):
                if key in payload:
                    payload = payload[key]
                    break

        rows = []
        if isinstance(payload, dict):
            for prefix, count in payload.items():
                if isinstance(count, dict) and isinstance(count.get("curie_prefix"), dict):
                    rows.extend(TranslatorNodeNormResolver._parse_prefix_counts_payload(count["curie_prefix"]))
                    continue
                rows.append({"prefix": str(prefix), "count": TranslatorNodeNormResolver._parse_prefix_count(count)})
        elif isinstance(payload, list):
            for item in payload:
                if isinstance(item, dict):
                    prefix = item.get("prefix") or item.get("curie_prefix") or item.get("id")
                    if prefix is None:
                        continue
                    count = item.get("count")
                    rows.append({"prefix": str(prefix), "count": count if isinstance(count, int) else None})
                else:
                    rows.append({"prefix": str(item), "count": None})

        return sorted(
            rows,
            key=lambda row: (
                -(row["count"] if isinstance(row.get("count"), int) else -1),
                row["prefix"].lower(),
            ),
        )

    @staticmethod
    def _parse_prefix_count(value):
        if isinstance(value, int):
            return value
        if isinstance(value, str):
            try:
                return int(value)
            except ValueError:
                return None
        return None

    def resolve_internal(self, input_nodes: List[Node]) -> Dict[str, List[IdMatch]]:

        input_ids = list(set([node.id for node in input_nodes]))

        result_list = {}

        for batch in yield_per(input_ids, self.batch_size):
            response_data = self._post_to_node_normalizer(batch)

            for input_id, results in response_data.items():
                res_obj: List[IdMatch] = []
                if results is not None:
                    res_obj = [
                        IdMatch(
                            input=input_id,
                            match=results['id']['identifier'],
                            equivalent_ids=[equiv_id['identifier'] for equiv_id in results['equivalent_identifiers']]
                        )]
                result_list[input_id] = res_obj
        return result_list

    def _post_to_node_normalizer(self, curies: Iterable[str]):
        batch = list(curies)
        post_body = {
            "curies": batch,
            "conflate": self.conflate_genes_and_proteins,
            "description": False,
            "drug_chemical_conflate": False
        }

        response = self._post_with_retries(batch, post_body)
        print(f"Resolved {len(batch)} from Translator Node Normalizer")
        return response.json()

    def _post_with_retries(self, batch: List[str], post_body: Dict):
        last_exception = None

        for attempt in range(1, self.max_retries + 1):
            try:
                response = requests.post(
                    self.node_norm_url(),
                    json=post_body,
                    timeout=self.request_timeout,
                )
            except requests.exceptions.RequestException as exc:
                last_exception = exc
                if attempt == self.max_retries:
                    break
                self._sleep_before_retry(attempt, batch, f"{type(exc).__name__}: {exc}")
                continue

            if response.status_code == 200:
                return response

            if response.status_code not in self.retryable_status_codes or attempt == self.max_retries:
                self._raise_node_normalizer_error(batch, response=response)

            self._sleep_before_retry(
                attempt,
                batch,
                f"{response.status_code} {response.text[:500]}",
            )

        self._raise_node_normalizer_error(batch, exception=last_exception)

    def _sleep_before_retry(self, attempt: int, batch: List[str], reason: str):
        sleep_seconds = self.retry_backoff_seconds * attempt
        print(
            f"Node Normalizer request failed for {len(batch)} IDs "
            f"(attempt {attempt}/{self.max_retries}): {reason}. "
            f"Retrying in {sleep_seconds} seconds."
        )
        time.sleep(sleep_seconds)

    @staticmethod
    def _sample_batch(batch: List[str]) -> str:
        return ", ".join(batch[:10])

    def _raise_node_normalizer_error(self, batch: List[str], response=None, exception: Exception = None):
        sample = self._sample_batch(batch)
        if response is not None:
            print(f"Failed to resolve IDs: {response.status_code} {response.text}")
            raise RuntimeError(
                f"Translator Node Normalizer failed with HTTP {response.status_code} "
                f"for {len(batch)} IDs after {self.max_retries} attempts. "
                f"Sample input IDs: {sample}"
            )
        raise RuntimeError(
            f"Translator Node Normalizer request failed for {len(batch)} IDs "
            f"after {self.max_retries} attempts. Sample input IDs: {sample}"
        ) from exception
