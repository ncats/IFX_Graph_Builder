import pytest
import requests

from src.constants import DataSourceName
from src.input_adapters.rdas.diseases import RDASRareDiseaseAdapter
from src.models.disease import Disease
from tests.registry_fakes import registry_dataset


class _FakeResponse:
    def __init__(self, payload, status_code=200):
        self.payload = payload
        self.status_code = status_code

    def raise_for_status(self):
        if self.status_code >= 400:
            raise requests.HTTPError(f"{self.status_code} Client Error")

    def json(self):
        return self.payload


def test_rdas_adapter_emits_rare_gard_diseases(monkeypatch, tmp_path):
    calls = []

    def fake_post(url, json, timeout):
        calls.append((url, json, timeout))
        offset = json["variables"]["offset"]
        if offset == 0:
            return _FakeResponse({
                "data": {
                    "diseases": [
                        {"gardId": "GARD:1", "gardName": "GRACILE syndrome"},
                        {"gardId": "GARD:0000005", "gardName": "Achondrogenesis"},
                    ]
                }
            })
        return _FakeResponse({"data": {"diseases": []}})

    monkeypatch.setattr("src.input_adapters.rdas.diseases.requests.post", fake_post)

    adapter = RDASRareDiseaseAdapter(
        graphql_url="https://example.org/graphql",
        batch_size=2,
        data_source=registry_dataset(
            tmp_path, source="rdas", dataset="diseases_graphql", version="test"
        ),
    )
    batches = list(adapter.get_all())
    diseases = [node for batch in batches for node in batch]

    assert adapter.get_datasource_name() == DataSourceName.RDAS
    assert adapter.get_version().version == "test"
    assert adapter.get_version().download_date is None
    assert all(isinstance(disease, Disease) for disease in diseases)
    assert [(disease.id, disease.name, disease.rare_disease) for disease in diseases] == [
        ("GARD:0000001", "GRACILE syndrome", True),
        ("GARD:0000005", "Achondrogenesis", True),
    ]
    assert calls[0][1]["variables"] == {"limit": 2, "offset": 0}
    assert "diseases(limit: $limit, offset: $offset)" in calls[0][1]["query"]


def test_rdas_adapter_dedupes_and_skips_invalid_gard_ids(monkeypatch, tmp_path):
    def fake_post(url, json, timeout):
        return _FakeResponse({
            "data": {
                "diseases": [
                    {"gardId": "GARD:1", "gardName": "First"},
                    {"gardId": "gard:0000001", "gardName": "Duplicate"},
                    {"gardId": "OMIM:1", "gardName": "Wrong prefix"},
                    {"gardId": "", "gardName": "Blank"},
                ]
            }
        })

    monkeypatch.setattr("src.input_adapters.rdas.diseases.requests.post", fake_post)

    adapter = RDASRareDiseaseAdapter(
        batch_size=10,
        data_source=registry_dataset(
            tmp_path, source="rdas", dataset="diseases_graphql", version="test"
        ),
    )
    diseases = [node for batch in adapter.get_all() for node in batch]

    assert len(diseases) == 1
    assert diseases[0].id == "GARD:0000001"
    assert diseases[0].rare_disease is True


def test_rdas_adapter_reports_schema_errors(monkeypatch, tmp_path):
    def fake_post(url, json, timeout):
        return _FakeResponse(
            {"errors": [{"message": 'Cannot query field "gards" on type "Query".'}]},
            status_code=400,
        )

    monkeypatch.setattr("src.input_adapters.rdas.diseases.requests.post", fake_post)
    adapter = RDASRareDiseaseAdapter(
        data_source=registry_dataset(
            tmp_path, source="rdas", dataset="diseases_graphql", version="test"
        )
    )

    with pytest.raises(RuntimeError, match='Cannot query field "gards"'):
        list(adapter.get_all())


def test_rdas_adapter_rejects_missing_disease_page(monkeypatch, tmp_path):
    def fake_post(url, json, timeout):
        return _FakeResponse({"data": {}})

    monkeypatch.setattr("src.input_adapters.rdas.diseases.requests.post", fake_post)
    adapter = RDASRareDiseaseAdapter(
        data_source=registry_dataset(
            tmp_path, source="rdas", dataset="diseases_graphql", version="test"
        )
    )

    with pytest.raises(ValueError, match="missing data.diseases list"):
        list(adapter.get_all())
