import pytest
import requests

from src.id_resolvers import node_normalizer
from src.id_resolvers.node_normalizer import TranslatorNodeNormResolver
from src.models.node import Node


class FakeResponse:
    def __init__(self, status_code, payload=None, text=""):
        self.status_code = status_code
        self._payload = payload or {}
        self.text = text

    def json(self):
        return self._payload

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(self.text)


def test_node_normalizer_retries_retryable_http_status(monkeypatch):
    calls = []
    responses = [
        FakeResponse(500, text="temporary service failure"),
        FakeResponse(
            200,
            {
                "MONDO:0000001": {
                    "id": {"identifier": "MONDO:0000001"},
                    "equivalent_identifiers": [{"identifier": "MONDO:0000001"}],
                }
            },
        ),
    ]

    def fake_post(url, json, timeout):
        calls.append((url, json, timeout))
        return responses.pop(0)

    monkeypatch.setattr(node_normalizer.requests, "post", fake_post)
    monkeypatch.setattr(node_normalizer.time, "sleep", lambda seconds: None)

    resolver = TranslatorNodeNormResolver(
        types=["Disease"],
        max_retries=2,
        retry_backoff_seconds=0,
    )

    results = resolver.resolve_internal([Node(id="MONDO:0000001")])

    assert len(calls) == 2
    assert calls[0][2] == 120
    assert results["MONDO:0000001"][0].match == "MONDO:0000001"


def test_node_normalizer_retries_request_exception(monkeypatch):
    calls = []

    def fake_post(url, json, timeout):
        calls.append((url, json, timeout))
        if len(calls) == 1:
            raise requests.exceptions.ConnectionError("connection dropped")
        return FakeResponse(
            200,
            {
                "DOID:1234": {
                    "id": {"identifier": "DOID:1234"},
                    "equivalent_identifiers": [{"identifier": "DOID:1234"}],
                }
            },
        )

    monkeypatch.setattr(node_normalizer.requests, "post", fake_post)
    monkeypatch.setattr(node_normalizer.time, "sleep", lambda seconds: None)

    resolver = TranslatorNodeNormResolver(
        types=["Disease"],
        max_retries=2,
        retry_backoff_seconds=0,
    )

    results = resolver.resolve_internal([Node(id="DOID:1234")])

    assert len(calls) == 2
    assert results["DOID:1234"][0].match == "DOID:1234"


def test_node_normalizer_reports_sample_ids_after_retries(monkeypatch):
    def fake_post(url, json, timeout):
        return FakeResponse(500, text="still broken")

    monkeypatch.setattr(node_normalizer.requests, "post", fake_post)
    monkeypatch.setattr(node_normalizer.time, "sleep", lambda seconds: None)

    resolver = TranslatorNodeNormResolver(
        types=["Disease"],
        max_retries=2,
        retry_backoff_seconds=0,
    )

    with pytest.raises(RuntimeError, match="Sample input IDs: MONDO:0000001"):
        resolver.resolve_internal([Node(id="MONDO:0000001")])


def test_node_normalizer_parses_prefix_count_payloads(monkeypatch):
    calls = []

    def fake_get(url, timeout, params=None):
        calls.append((url, timeout, params))
        if params == {"semantic_type": "biolink:Disease"}:
            return FakeResponse(200, {
                "biolink:Disease": {
                    "curie_prefix": {"MONDO": "3", "DOID": "2"},
                },
            })
        if params == {"semantic_type": "biolink:SmallMolecule"}:
            return FakeResponse(200, {
                "biolink:SmallMolecule": {
                    "curie_prefix": {"CHEBI": "5", "MESH": "1"},
                },
            })
        raise AssertionError((url, params))

    monkeypatch.setattr(node_normalizer.requests, "get", fake_get)

    resolver = TranslatorNodeNormResolver(
        types=["Condition", "Disease", "Ligand"],
        request_timeout=7,
    )

    assert resolver.get_prefix_counts() == [
        {"prefix": "CHEBI", "count": 5},
        {"prefix": "MONDO", "count": 3},
        {"prefix": "DOID", "count": 2},
        {"prefix": "MESH", "count": 1},
    ]
    assert calls == [
        (resolver.node_norm_prefixes_url(), 7, {"semantic_type": "biolink:Disease"}),
        (resolver.node_norm_prefixes_url(), 7, {"semantic_type": "biolink:SmallMolecule"}),
    ]

    assert TranslatorNodeNormResolver._parse_prefix_counts_payload(["MONDO", "DOID"]) == [
        {"prefix": "DOID", "count": None},
        {"prefix": "MONDO", "count": None},
    ]


def test_node_normalizer_returns_examples_for_registered_types():
    resolver = TranslatorNodeNormResolver(
        types=["Condition", "Disease", "Ligand"],
    )

    assert resolver.get_example_ids(limit=4) == [
        "MONDO:0005148",
        "DOID:9352",
        "CHEBI:15377",
        "PUBCHEM.COMPOUND:2244",
    ]


def test_node_normalizer_uses_configured_production_endpoint():
    resolver = TranslatorNodeNormResolver(
        types=["Ligand"],
        base_url="https://nodenorm.transltr.io/1.5/",
    )

    assert resolver.node_norm_url() == "https://nodenorm.transltr.io/1.5/get_normalized_nodes"
    assert resolver.node_norm_prefixes_url() == "https://nodenorm.transltr.io/1.5/get_curie_prefixes"
    assert resolver.name.endswith("https://nodenorm.transltr.io/1.5")
