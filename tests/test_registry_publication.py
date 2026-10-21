from datetime import date, datetime
from pathlib import Path


from src.core.registry_publication import publish_derived_file, publish_source_file


def test_publish_source_file_uses_standalone_registry_client(monkeypatch, tmp_path: Path):
    path = tmp_path / "records.tsv"
    path.write_text("id\n1\n", encoding="utf-8")
    captured = {}

    class Client:
        def publish_source(self, snapshot_id, **kwargs):  # type: ignore[no-untyped-def]
            captured["snapshot_id"] = snapshot_id
            captured.update(kwargs)
            return object()

    monkeypatch.setattr(
        "src.core.registry_publication.RegistryClient.connect",
        lambda credentials: Client(),
    )

    publish_source_file(
        "example:records:2026-09-11",
        path,
        "registry.yaml",
        capture_method="provider_export",
        source_url="https://example.org/records.tsv",
        metadata={"producer": "test"},
    )

    assert captured["snapshot_id"] == "example:records:2026-09-11"
    assert captured["files"] == {"records.tsv": path}
    assert captured["version_date"] == date(2026, 9, 11)
    assert isinstance(captured["captured_at"], datetime)
    assert captured["captured_at"].tzinfo is not None
    assert captured["capture_method"] == "provider_export"
    assert captured["file_sources"] == {"records.tsv": "https://example.org/records.tsv"}
    assert captured["upstream_urls"] == ("https://example.org/records.tsv",)


def test_publish_derived_file_requires_kind_qualified_inputs(monkeypatch, tmp_path: Path):
    path = tmp_path / "records.tsv"
    path.write_text("id\n1\n", encoding="utf-8")
    captured = {}

    class Derived:
        def publish(self, snapshot_id, **kwargs):  # type: ignore[no-untyped-def]
            captured["snapshot_id"] = snapshot_id
            captured.update(kwargs)
            return object()

    class Client:
        derived = Derived()

    monkeypatch.setattr(
        "src.core.registry_publication.RegistryClient.connect",
        lambda credentials: Client(),
    )

    publish_derived_file(
        "target_graph:gene_ids:2026-09-11",
        path,
        "registry.yaml",
        inputs=["source=ensembl:genes:116", "external=nodenorm:api:2026jul22"],
        producer_release="2.2.0",
        producer_repository="https://github.com/ncats/IFX_Harmonizers",
        producer_revision="a" * 40,
        transform_name="target_harmonizer_export",
        validation={"rows": 1},
    )

    assert captured["snapshot_id"] == "target_graph:gene_ids:2026-09-11"
    assert [ref.kind.value for ref in captured["inputs"]] == [
        "source_snapshot",
        "external_dataset_version",
    ]
    assert captured["producer"].name == "ifx_harmonizers"
    assert captured["transform"] == {
        "name": "target_harmonizer_export",
        "version": "2.2.0",
    }
