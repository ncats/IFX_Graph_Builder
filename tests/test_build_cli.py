from types import SimpleNamespace

import pytest

from src.core.object_factory import create_object_from_config
from src.output_adapters.arango_output_adapter import ArangoOutputAdapter
from src.use_cases import build_cli


def _builder(database_name="from_yaml", adapters=None):
    adapter = object.__new__(ArangoOutputAdapter)
    adapter.database_name = database_name
    return SimpleNamespace(
        yaml_file="src/use_cases/working.yaml",
        etl=SimpleNamespace(output_adapters=[adapter] if adapters is None else adapters),
        prepare_datastore=lambda: None,
    )


def test_confirmation_uses_configured_output_database(monkeypatch):
    builder = _builder("targets_251_legacy_id_files")
    prompts = []
    monkeypatch.setattr("builtins.input", lambda prompt: prompts.append(prompt) or "yes")

    assert build_cli.prepare_primary_builder(builder, resume=False, yes=False)
    assert "'targets_251_legacy_id_files'" in prompts[0]


def test_yaml_imported_arango_adapter_passes_validation(monkeypatch):
    adapter = create_object_from_config({
        "import": "./src/output_adapters/arango_output_adapter.py",
        "class": "ArangoOutputAdapter",
        "kwargs": {"database_name": "targets_251_legacy_id_files"},
        "credentials": {"url": "http://localhost:8529"},
    })
    builder = _builder(adapters=[adapter])
    prompts = []
    monkeypatch.setattr("builtins.input", lambda prompt: prompts.append(prompt) or "no")

    assert not build_cli.prepare_primary_builder(builder, resume=False, yes=False)
    assert "'targets_251_legacy_id_files'" in prompts[0]


def test_resume_error_uses_configured_output_database(monkeypatch, capsys):
    builder = _builder("targets_251_consumer_contract")
    monkeypatch.setattr(build_cli, "_arango_database_exists", lambda adapter: False)

    assert not build_cli.prepare_primary_builder(builder, resume=True, yes=False)
    assert "Cannot resume 'targets_251_consumer_contract'" in capsys.readouterr().out


@pytest.mark.parametrize("adapters", [[], [object()], [object(), object()]])
def test_invalid_output_configuration_fails_before_build(adapters):
    builder = _builder(adapters=adapters)

    with pytest.raises(ValueError, match="src/use_cases/working.yaml must configure exactly one ArangoOutputAdapter"):
        build_cli.prepare_primary_builder(builder, resume=False, yes=True)


@pytest.mark.parametrize("refresh_post_on_resume", [True, False])
def test_resume_can_refresh_post_build_without_resetting_its_run(
    monkeypatch, refresh_post_on_resume,
):
    calls = []

    class Builder:
        def __init__(self, yaml_file):
            self.yaml_file = yaml_file

        def do_etl(self, **kwargs):
            calls.append((self.yaml_file, kwargs))

    monkeypatch.setattr("src.use_cases.build_from_yaml.BuildGraphFromYaml", Builder)
    monkeypatch.setattr(build_cli, "parse_common_build_args", lambda build_name: SimpleNamespace(resume=True, yes=False))
    monkeypatch.setattr(build_cli, "prepare_primary_builder", lambda builder, **kwargs: True)

    build_cli.run_arango_build_cli(
        build_name="pounce",
        primary_yaml="primary.yaml",
        post_yaml="post.yaml",
        refresh_post_on_resume=refresh_post_on_resume,
    )

    assert calls == [
        ("primary.yaml", {"clean_edges": False, "resume": True}),
        ("post.yaml", {
            "resume": True,
            "rerun_completed": refresh_post_on_resume,
            "run_id": "post.yaml:refresh" if refresh_post_on_resume else None,
        }),
    ]
