from datetime import UTC, date, datetime
from pathlib import Path

import pytest

from src.core.registry_integration import RegistryIntegration, RegistryReference
from src.core.config import Config, resolve_registry_references
from src.models.registry_dataset import RegistryDatasetKind


class Result:
    source = "example"
    dataset = "records"
    version = "1"
    version_date = date(2026, 9, 10)
    download_date = date(2026, 9, 11)
    registered_at = datetime(2026, 9, 12, tzinfo=UTC)
    snapshot_id = "example:records:1"
    manifest_uri = "s3://registry/example/records/1/manifest.yaml"
    local_dir: Path
    manifest = {
        "kind": "source_snapshot",
        "files": [
            {
                "path": "records.tsv",
                "sha256": "a" * 64,
                "size_bytes": 1,
                "storage_uri": "s3://registry/example/records/1/records.tsv",
            }
        ],
    }


class RecordingEndpoint:
    def __init__(self, result: Result):
        self.result = result
        self.calls: list[tuple[str, Path | None]] = []

    def materialize(self, snapshot_id, *, destination):  # type: ignore[no-untyped-def]
        self.calls.append((snapshot_id, Path(destination)))
        return self.result

    def describe(self, snapshot_id):  # type: ignore[no-untyped-def]
        self.calls.append((snapshot_id, None))
        return self.result


class Client:
    def __init__(self, result: Result):
        self.source = RecordingEndpoint(result)
        self.derived = RecordingEndpoint(result)
        self.external = RecordingEndpoint(result)

    def materialize(self, snapshot_id, *, destination):  # type: ignore[no-untyped-def]
        return self.source.materialize(snapshot_id, destination=destination)


def _result(tmp_path: Path) -> Result:
    result = Result()
    result.local_dir = tmp_path
    (tmp_path / "records.tsv").write_text("x", encoding="utf-8")
    return result


def test_source_shorthand_dispatches_only_to_source_client(tmp_path: Path) -> None:
    client = Client(_result(tmp_path))

    resolved = RegistryIntegration(client, tmp_path).resolve("example:records:1")

    assert resolved.kind is RegistryDatasetKind.SOURCE
    assert resolved.file() == tmp_path / "records.tsv"
    assert resolved.version_info().version_date == date(2026, 9, 10)
    assert client.source.calls == [("example:records:1", tmp_path)]
    assert client.derived.calls == []
    assert client.external.calls == []


def test_explicit_derived_dispatches_without_kind_probing(tmp_path: Path) -> None:
    result = _result(tmp_path)
    result.manifest = {**result.manifest, "kind": "derived_snapshot"}
    client = Client(result)

    resolved = RegistryIntegration(client, tmp_path).resolve(
        {"kind": "derived_snapshot", "snapshot_id": "example:records:1"}
    )

    assert resolved.kind is RegistryDatasetKind.DERIVED
    assert client.source.calls == []
    assert client.derived.calls == [("example:records:1", tmp_path)]
    assert client.external.calls == []


def test_explicit_external_describes_without_materializing(tmp_path: Path) -> None:
    result = _result(tmp_path)
    result.manifest = {**result.manifest, "kind": "external_dataset_version"}
    client = Client(result)

    resolved = RegistryIntegration(client, tmp_path).resolve(
        {"kind": "external_dataset_version", "snapshot_id": "example:records:1"}
    )

    assert resolved.kind is RegistryDatasetKind.EXTERNAL
    assert resolved.local_dir is None
    assert resolved.download_date is None
    assert resolved.to_metadata()["registered_at"] == "2026-09-12T00:00:00+00:00"
    assert client.source.calls == []
    assert client.derived.calls == []
    assert client.external.calls == [("example:records:1", None)]
    with pytest.raises(ValueError, match="metadata-only"):
        resolved.file()


def test_reference_requires_explicit_kind_for_non_source() -> None:
    assert RegistryReference.parse("example:records:1").kind.value == "source_snapshot"
    with pytest.raises(ValueError, match="requires kind and snapshot_id"):
        RegistryReference.parse({"snapshot_id": "example:records:1"})
    with pytest.raises(ValueError, match="unsupported fields"):
        RegistryReference.parse(
            {
                "kind": "derived_snapshot",
                "snapshot_id": "example:records:1",
                "fallback": True,
            }
        )


def test_client_errors_propagate_without_legacy_fallback(tmp_path: Path) -> None:
    class FailingClient(Client):
        def materialize(self, snapshot_id, *, destination):  # type: ignore[no-untyped-def]
            raise LookupError(snapshot_id)

    integration = RegistryIntegration(FailingClient(_result(tmp_path)), tmp_path)

    with pytest.raises(LookupError, match="example:records:missing"):
        integration.resolve("example:records:missing")


def test_config_resolves_dataset_arguments_with_one_connection(monkeypatch) -> None:
    calls = []

    class Integration:
        def resolve(self, value):  # type: ignore[no-untyped-def]
            calls.append(value)
            return f"resolved:{value}"

    monkeypatch.setattr(
        RegistryIntegration,
        "connect",
        lambda config: Integration(),
    )
    config = {
        "registry": {"credentials": "registry.yaml"},
        "resolvers": [{
            "kwargs": {
                "data_source": "example:records:1",
                "protein_data_source": "example:mapping:2",
            },
        }],
        "input_adapters": [{
            "kwargs": {
                "data_source": {
                    "kind": "derived_snapshot",
                    "snapshot_id": "example:derived:3",
                },
            },
        }],
        "output_adapters": [{"kwargs": {"data_source": "leave:untouched:1"}}],
    }

    resolved = resolve_registry_references(config)

    assert calls == [
        "example:records:1",
        "example:mapping:2",
        {"kind": "derived_snapshot", "snapshot_id": "example:derived:3"},
    ]
    assert resolved["resolvers"][0]["kwargs"]["data_source"] == "resolved:example:records:1"
    assert resolved["output_adapters"] == config["output_adapters"]


def test_config_preserves_registry_credential_yaml_as_a_path(monkeypatch, tmp_path) -> None:
    credentials_path = tmp_path / "registry_credentials.yaml"
    credentials_path.write_text("bucket: should-not-be-expanded\n", encoding="utf-8")
    config_path = tmp_path / "build.yaml"
    config_path.write_text(
        "registry:\n"
        f"  credentials: {credentials_path}\n"
        "input_adapters: []\n",
        encoding="utf-8",
    )
    captured = []

    class Integration:
        pass

    monkeypatch.setattr(
        RegistryIntegration,
        "connect",
        lambda config: captured.append(config) or Integration(),
    )

    loaded = Config(config_path).config_dict

    assert loaded["registry"]["credentials"] == str(credentials_path)
    assert captured == [{"credentials": str(credentials_path)}]
