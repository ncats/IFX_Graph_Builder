from datetime import timedelta
from types import SimpleNamespace

import yaml

from ifx_registry import AuditDisposition, SnapshotRef

from src.use_cases import audit_registry_freshness


def test_collects_all_registry_arguments_and_deduplicates() -> None:
    config = {
        "resolvers": [
            {"kwargs": {"data_source": "uniprot:human:2026_02"}},
        ],
        "input_adapters": [
            {
                "kwargs": {
                    "data_source": {
                        "kind": "derived_snapshot",
                        "snapshot_id": "pubchem:cid_molecular_info:deps-1",
                    },
                    "uniprot_data_source": "uniprot:human:2026_02",
                }
            }
        ],
        "output_adapters": [
            {"kwargs": {"data_source": "ignored:output:1"}},
        ],
    }

    assert audit_registry_freshness.collect_registry_references(config) == (
        SnapshotRef.source("uniprot:human:2026_02"),
        SnapshotRef.derived("pubchem:cid_molecular_info:deps-1"),
    )


def test_collects_pins_from_nested_yaml_but_does_not_open_credentials(
    monkeypatch,
    tmp_path,
) -> None:
    nested = tmp_path / "nested.yaml"
    nested.write_text(
        yaml.safe_dump(
            {
                "kwargs": {
                    "data_source": "reactome:pathways:98",
                    "credentials": "missing-secret.yaml",
                }
            }
        ),
        encoding="utf-8",
    )
    config = {"input_adapters": ["nested.yaml"]}
    monkeypatch.chdir(tmp_path)

    assert audit_registry_freshness.collect_registry_references(config) == (
        SnapshotRef.source("reactome:pathways:98"),
    )


def test_empty_registry_configuration_is_rejected(tmp_path) -> None:
    config_path = tmp_path / "build.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "registry": {},
                "input_adapters": [
                    {"kwargs": {"data_source": "reactome:pathways:97"}}
                ],
            }
        ),
        encoding="utf-8",
    )

    try:
        audit_registry_freshness.audit_build_yaml(config_path)
    except ValueError as error:
        assert "no registry configuration" in str(error)
    else:
        raise AssertionError("empty registry configuration was accepted")


def test_audits_yaml_with_registry_connection_settings(monkeypatch, tmp_path) -> None:
    config_path = tmp_path / "build.yaml"
    config_path.write_text(
        yaml.safe_dump(
            {
                "registry": {
                    "credentials": "registry.yaml",
                    "bucket": "example-bucket",
                },
                "input_adapters": [
                    {"kwargs": {"data_source": "reactome:pathways:97"}}
                ],
            }
        ),
        encoding="utf-8",
    )
    calls = []
    expected_report = SimpleNamespace(entries=(), is_current=True)

    class Client:
        def audit(self, references, *, timeout):  # type: ignore[no-untyped-def]
            calls.append((references, timeout))
            return expected_report

    def connect(credentials, **kwargs):  # type: ignore[no-untyped-def]
        calls.append((credentials, kwargs))
        return Client()

    monkeypatch.setattr(audit_registry_freshness.RegistryAuditClient, "connect", connect)

    report = audit_registry_freshness.audit_build_yaml(
        config_path,
        timeout=timedelta(seconds=12),
    )

    assert report is expected_report
    assert calls == [
        (
            "registry.yaml",
            {
                "bucket": "example-bucket",
                "region": None,
                "prefix": "",
                "source_configuration": None,
                "cure_credentials_file": None,
            },
        ),
        ((SnapshotRef.source("reactome:pathways:97"),), timedelta(seconds=12)),
    ]


def test_main_returns_nonzero_when_build_is_not_current(monkeypatch, capsys) -> None:
    entry = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=SnapshotRef.source("reactome:pathways:97"),
        reason="A checked replacement is registered",
        recommended_snapshot_id="reactome:pathways:98",
        latest_upstream_version=SimpleNamespace(value="98"),
    )
    report = SimpleNamespace(
        roots=(entry.reference,),
        entries=(entry,),
        is_current=False,
    )
    monkeypatch.setattr(
        audit_registry_freshness,
        "audit_build_yaml",
        lambda *args, **kwargs: report,
    )

    assert audit_registry_freshness.main(["ramp.yaml"]) == 1
    output = capsys.readouterr().out
    assert "Update YAML to newest registered pin (1)" in output
    assert "reactome:pathways:97 -> reactome:pathways:98" in output
    assert (
        "Result: 1 known updates · 0 dependent rebuilds · "
            "0 manual sources to confirm · 0 registered with freshness caveats · "
        "0 blocked · 0 fully up to date"
    ) in output


def test_report_groups_actions_and_preserves_entry_order(capsys) -> None:
    def entry(disposition, snapshot_id, **kwargs):  # type: ignore[no-untyped-def]
        defaults = {
            "reason": "explanation",
            "recommended_snapshot_id": None,
            "latest_upstream_version": None,
        }
        return SimpleNamespace(
            disposition=disposition,
            reference=SnapshotRef.source(snapshot_id),
            **(defaults | kwargs),
        )

    report = SimpleNamespace(
        roots=(
            SnapshotRef.source("hmdb:metabolites_xml:5.0"),
            SnapshotRef.source("chebi:three_star_sdf:old"),
            SnapshotRef.source("wikipathways:rdf_wp:old"),
            SnapshotRef.source("reactome:pathways:97"),
        ),
        is_current=False,
        entries=(
            entry(
                AuditDisposition.UNVERIFIABLE,
                "hmdb:metabolites_xml:5.0",
                reason="No automatic upstream check is available.",
            ),
            entry(
                AuditDisposition.REGISTER_SOURCE,
                "chebi:three_star_sdf:old",
                latest_upstream_version=SimpleNamespace(value="new"),
            ),
            entry(
                AuditDisposition.UPDATE_PIN,
                "wikipathways:rdf_wp:old",
                recommended_snapshot_id="wikipathways:rdf_wp:new",
            ),
            entry(AuditDisposition.CURRENT, "reactome:pathways:97"),
        ),
    )

    audit_registry_freshness.print_report(report, "ramp.yaml")

    output = capsys.readouterr().out
    assert output.index("Refresh Registry") < output.index("Update YAML")
    assert output.index("Update YAML") < output.index("Needs manual check")
    assert output.index("Needs manual check") < output.index("Up to date")
    assert "chebi:three_star_sdf:old -> chebi:three_star_sdf:new (not registered yet)" in output
    assert "wikipathways:rdf_wp:old -> wikipathways:rdf_wp:new" in output
    assert "reactome:pathways:97\n" in output
    assert (
        "Result: 2 known updates · 0 dependent rebuilds · "
            "0 manual sources to confirm · 0 registered with freshness caveats · "
        "0 blocked · 1 fully up to date"
    ) in output


def test_registration_only_does_not_claim_yaml_needs_an_update(capsys) -> None:
    entry = SimpleNamespace(
        disposition=AuditDisposition.REGISTER_SOURCE,
        reference=SnapshotRef.source("example:records:2"),
        reason="Version 2 is not registered",
        recommended_snapshot_id=None,
        latest_upstream_version=SimpleNamespace(value="2"),
    )
    report = SimpleNamespace(
        roots=(entry.reference,),
        entries=(entry,),
        is_current=False,
    )

    audit_registry_freshness.print_report(report, "build.yaml")

    output = capsys.readouterr().out
    assert "\nRefresh Registry (1)\n" in output
    assert "Refresh Registry, then update YAML" not in output
    assert "example:records:2 (not registered yet)" in output
    assert "example:records:2 -> example:records:2" not in output


def test_caveated_root_pin_update_keeps_the_action_and_reason(capsys) -> None:
    reference = SnapshotRef.derived("pubchem:cid_molecular_info:deps-old")
    caveat = SimpleNamespace(
        code=SimpleNamespace(value="upstream_check_unverified"),
        message="Live upstream freshness could not be verified for refmet.",
    )
    entry = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=reference,
        reason="A registered derived snapshot matches newest registered inputs.",
        recommended_snapshot_id="pubchem:cid_molecular_info:deps-new",
        latest_upstream_version=None,
        caveats=(caveat,),
    )
    report = SimpleNamespace(roots=(reference,), entries=(entry,), is_current=False)

    audit_registry_freshness.print_report(report, "ramp.yaml")

    output = capsys.readouterr().out
    assert "Update YAML to newest registered pin (1)" in output
    assert "deps-old -> pubchem:cid_molecular_info:deps-new" in output
    assert "A registered derived snapshot matches newest registered inputs." in output
    assert "Caveat: Live upstream freshness could not be verified for refmet." in output


def test_current_upstream_check_caveat_uses_neutral_heading(capsys) -> None:
    caveat = SimpleNamespace(
        code=SimpleNamespace(value="upstream_check_unverified"),
        message="Live upstream freshness could not be verified for refmet.",
    )
    entry = SimpleNamespace(
        disposition=AuditDisposition.CURRENT,
        reference=SnapshotRef.source("refmet:metabolites_csv:sha256-current"),
        reason="The exact pin is registered; live upstream freshness could not be verified.",
        caveats=(caveat,),
    )
    report = SimpleNamespace(roots=(entry.reference,), entries=(entry,), is_current=False)

    audit_registry_freshness.print_report(report, "ramp.yaml")

    output = capsys.readouterr().out
    assert "Registered pin — freshness caveats (1)" in output
    assert "Up to date — manual source check needed" not in output


def test_report_separates_manual_caveats_from_dependent_rebuilds(capsys) -> None:
    manual_ref = SnapshotRef.source("hmdb:metabolites_xml:5.0")
    caveat = SimpleNamespace(
        origin=manual_ref,
        code=SimpleNamespace(value="manual_freshness"),
        message="Confirm that HMDB 5.0 is still current.",
    )
    entries = (
        SimpleNamespace(
            disposition=AuditDisposition.UNVERIFIABLE,
            reference=manual_ref,
            reason="No automatic upstream check is available.",
            caveats=(caveat,),
        ),
        SimpleNamespace(
            disposition=AuditDisposition.REBUILD_DERIVED,
            reference=SnapshotRef.derived("pubchem:compound_cid_set:deps-old"),
            reason="Current inputs produce deps-new; rebuild this dataset.",
            recommended_snapshot_id=None,
            caveats=(caveat,),
        ),
        SimpleNamespace(
            disposition=AuditDisposition.BLOCKED,
            reference=SnapshotRef.derived("pubchem:compound_records:deps-old"),
            reason=(
                "After rebuilding pubchem:compound_cid_set:deps-old, "
                "rebuild this derived dataset"
            ),
            caveats=(caveat,),
        ),
        SimpleNamespace(
            disposition=AuditDisposition.CURRENT,
            reference=SnapshotRef.derived("example:derived:1"),
            reason="Current under automatic checks; confirm HMDB manually.",
            caveats=(caveat,),
        ),
    )
    report = SimpleNamespace(
        roots=(manual_ref, entries[1].reference, entries[3].reference),
        entries=entries,
        is_current=False,
    )

    audit_registry_freshness.print_report(report, "ramp.yaml")

    output = capsys.readouterr().out
    assert "Rebuild in Registry, then update YAML (1)" in output
    assert "Then rebuild dependent datasets (1)" in output
    assert "Needs manual check (1)" in output
    assert "Needs manual source confirmation" not in output
    assert "Up to date — manual source check needed (1)" in output
    assert output.count("Caveat: Confirm that HMDB 5.0 is still current.") == 1
    assert (
        "Result: 1 known updates · 1 dependent rebuilds · "
            "1 manual sources to confirm · 1 registered with freshness caveats · "
        "0 blocked · 0 fully up to date"
    ) in output


def test_transitive_updates_do_not_claim_the_yaml_contains_the_pin(capsys) -> None:
    source = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=SnapshotRef.source("example:records:1"),
        reason="Use version 2",
        recommended_snapshot_id="example:records:2",
        latest_upstream_version=SimpleNamespace(value="2"),
        caveats=(),
    )
    derived = SimpleNamespace(
        disposition=AuditDisposition.REBUILD_DERIVED,
        reference=SnapshotRef.derived("derived:output:deps-old"),
        reason="Current inputs produce deps-new; rebuild this dataset.",
        recommended_snapshot_id=None,
        caveats=(),
        dependencies=(source.reference,),
    )
    report = SimpleNamespace(
        roots=(derived.reference,),
        entries=(source, derived),
        is_current=False,
    )

    audit_registry_freshness.print_report(report, "build.yaml")

    output = capsys.readouterr().out
    assert "Use newest registered inputs when rebuilding (1)" in output
    assert "Rebuild in Registry, then update YAML (1)" in output
    assert "Update YAML to newest registered pin" not in output
    assert "Result: 1 known updates" in output


def test_transitive_updates_are_hidden_without_a_rebuild_action(capsys) -> None:
    root = SnapshotRef.derived("derived:output:deps-old")
    source = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=SnapshotRef.source("example:records:1"),
        reason="Use version 2",
        recommended_snapshot_id="example:records:2",
        latest_upstream_version=None,
        caveats=(),
    )
    derived = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=root,
        reason="A registered derived snapshot matches current inputs.",
        recommended_snapshot_id="derived:output:deps-new",
        latest_upstream_version=None,
        caveats=(),
    )
    report = SimpleNamespace(roots=(root,), entries=(source, derived), is_current=False)

    audit_registry_freshness.print_report(report, "build.yaml")

    output = capsys.readouterr().out
    assert "Update YAML to newest registered pin (1)" in output
    assert "Use newest registered inputs when rebuilding" not in output


def test_transitive_updates_are_scoped_to_their_rebuild_root(capsys) -> None:
    stale_source = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=SnapshotRef.source("stale:records:1"),
        reason="Use version 2",
        recommended_snapshot_id="stale:records:2",
        latest_upstream_version=None,
        caveats=(),
    )
    update_root = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=SnapshotRef.derived("derived:already_built:old"),
        reason="A replacement is registered.",
        recommended_snapshot_id="derived:already_built:new",
        latest_upstream_version=None,
        caveats=(),
        dependencies=(stale_source.reference,),
    )
    rebuild_source = SimpleNamespace(
        disposition=AuditDisposition.UPDATE_PIN,
        reference=SnapshotRef.source("rebuild:records:1"),
        reason="Use version 2",
        recommended_snapshot_id="rebuild:records:2",
        latest_upstream_version=None,
        caveats=(),
    )
    rebuild_root = SimpleNamespace(
        disposition=AuditDisposition.REBUILD_DERIVED,
        reference=SnapshotRef.derived("derived:needs_rebuild:old"),
        reason="Rebuild required.",
        recommended_snapshot_id=None,
        caveats=(),
        dependencies=(rebuild_source.reference,),
    )
    report = SimpleNamespace(
        roots=(update_root.reference, rebuild_root.reference),
        entries=(stale_source, update_root, rebuild_source, rebuild_root),
        is_current=False,
    )

    audit_registry_freshness.print_report(report, "build.yaml")

    output = capsys.readouterr().out
    assert "rebuild:records:1 -> rebuild:records:2" in output
    assert "stale:records:1 -> stale:records:2" not in output


def test_registration_prerequisite_is_a_dependent_rebuild(capsys) -> None:
    entry = SimpleNamespace(
        disposition=AuditDisposition.BLOCKED,
        reference=SnapshotRef.derived("derived:output:deps-old"),
        reason=(
            "After registering example:records:2, rebuild this derived dataset"
        ),
        caveats=(),
    )
    report = SimpleNamespace(
        roots=(entry.reference,),
        entries=(entry,),
        is_current=False,
    )

    audit_registry_freshness.print_report(report, "build.yaml")

    output = capsys.readouterr().out
    assert "Then rebuild dependent datasets (1)" in output
    assert "\nBlocked" not in output
    assert "1 dependent rebuilds" in output


def test_shared_source_report_traces_nested_derived_inputs_and_compares_pins(
    capsys,
) -> None:
    direct = (
        SnapshotRef.source("babel:human_gene_compendium:2026jul22"),
        SnapshotRef.source("uniprot:human:2026_02"),
        SnapshotRef.source("ncbi:gene_summary:2026-06-17"),
    )
    targets = SnapshotRef.derived("ifx_harmonizers:targets:2.4.0")
    baseline = SnapshotRef.derived("ifx_harmonizers:targets:2.3.0")
    other = SnapshotRef.derived("example:other:1")
    entries = (
        SimpleNamespace(reference=direct[0], dependencies=()),
        SimpleNamespace(reference=direct[1], dependencies=()),
        SimpleNamespace(reference=direct[2], dependencies=()),
        SimpleNamespace(
            reference=baseline,
            dependencies=(
                SnapshotRef.source("babel:human_gene_compendium:2026jul22"),
            ),
        ),
        SimpleNamespace(
            reference=targets,
            dependencies=(
                baseline,
                SnapshotRef.source("uniprot:human_idmapping:2026_03"),
                SnapshotRef.source("uniprot:human_uniref100_sparql:2026_03-export1"),
                SnapshotRef.source("ncbi:human_gene_info:2026-10-05"),
            ),
        ),
        SimpleNamespace(
            reference=other,
            dependencies=(SnapshotRef.source("uniprot:human:2026_02"),),
        ),
    )
    report = SimpleNamespace(
        roots=(*direct, targets, other), entries=entries, is_current=True
    )

    comparisons = audit_registry_freshness.compare_shared_source_pins(report)

    assert [item.result for item in comparisons] == [
        "MATCH",  # Exact UniProt pin in the other derived root.
        "MATCH",  # Babel reached through the nested 2.3.0 snapshot.
        "REVIEW",  # The derived release has no NCBI gene summary pin.
        "REVIEW",  # The derived release has no UniProt human pin.
    ]
    printable_report = SimpleNamespace(
        roots=report.roots,
        entries=tuple(
            SimpleNamespace(
                **vars(entry),
                disposition=AuditDisposition.CURRENT,
                caveats=(),
                reason="Current",
            )
            for entry in entries
        ),
        is_current=True,
    )
    audit_registry_freshness.print_report(printable_report, "build.yaml")
    output = capsys.readouterr().out
    assert "Shared source pins (4 comparisons)" in output
    assert "REVIEW (same source, different datasets): ncbi:gene_summary" in output
    assert "YAML: ncbi:gene_summary:2026-06-17" in output
    assert "Derived: no ncbi:gene_summary pin" in output
    assert "Other ncbi inputs:" in output
    assert "ncbi:human_gene_info:2026-10-05" in output
    assert "uniprot:human_idmapping:2026_03" in output
    assert "Review whether these datasets should share a pin" in output
    assert "via ifx_harmonizers:targets:2.4.0" in output
    assert "Pin consistency: 2 matches · 0 mismatches · " in output
    assert "2 same-source dataset reviews · 0 unchecked derived roots" in output


def test_shared_source_mismatch_makes_current_freshness_report_fail(
    monkeypatch, capsys
) -> None:
    direct = SnapshotRef.source("example:records:1")
    derived = SnapshotRef.derived("example:derived:1")
    report = SimpleNamespace(
        roots=(direct, derived),
        entries=(
            SimpleNamespace(
                reference=direct,
                disposition=AuditDisposition.CURRENT,
                dependencies=(),
                caveats=(),
            ),
            SimpleNamespace(
                reference=SnapshotRef.source("example:records:2"),
                disposition=AuditDisposition.CURRENT,
                dependencies=(),
                caveats=(),
            ),
            SimpleNamespace(
                reference=derived,
                disposition=AuditDisposition.CURRENT,
                dependencies=(SnapshotRef.source("example:records:2"),),
                caveats=(),
            ),
        ),
        is_current=True,
    )
    monkeypatch.setattr(
        audit_registry_freshness, "audit_build_yaml", lambda *args, **kwargs: report
    )

    assert audit_registry_freshness.main(["build.yaml"]) == 1
    output = capsys.readouterr().out
    assert "MISMATCH: example:records" in output
    assert "YAML: example:records:1" in output
    assert "Derived: example:records:2" in output


def test_same_source_different_dataset_is_review_not_pin_conflict(
    monkeypatch, capsys
) -> None:
    direct = SnapshotRef.source("example:records:1")
    other = SnapshotRef.source("example:other_records:2")
    derived = SnapshotRef.derived("example:derived:1")
    report = SimpleNamespace(
        roots=(direct, derived),
        entries=(
            SimpleNamespace(
                reference=direct,
                disposition=AuditDisposition.CURRENT,
                dependencies=(),
                caveats=(),
            ),
            SimpleNamespace(
                reference=other,
                disposition=AuditDisposition.CURRENT,
                dependencies=(),
                caveats=(),
            ),
            SimpleNamespace(
                reference=derived,
                disposition=AuditDisposition.CURRENT,
                dependencies=(other,),
                caveats=(),
            ),
        ),
        is_current=True,
    )
    monkeypatch.setattr(
        audit_registry_freshness, "audit_build_yaml", lambda *args, **kwargs: report
    )

    assert audit_registry_freshness.main(["build.yaml"]) == 0
    output = capsys.readouterr().out
    assert "REVIEW (same source, different datasets): example:records" in output
    assert "Derived: no example:records pin" in output
    assert "example:other_records:2" in output


def test_shared_source_report_handles_unavailable_derived_dependencies() -> None:
    direct = SnapshotRef.source("uniprot:human:2026_02")
    derived = SnapshotRef.derived("example:missing:1")
    report = SimpleNamespace(roots=(direct, derived), entries=(), is_current=False)

    assert audit_registry_freshness.compare_shared_source_pins(report) == ()
    assert audit_registry_freshness._unchecked_derived_roots(report) == (derived,)


def test_nested_missing_lineage_is_reported_as_unchecked(capsys) -> None:
    direct = SnapshotRef.source("uniprot:human:2026_02")
    derived = SnapshotRef.derived("example:output:1")
    missing = SnapshotRef.derived("example:missing:1")
    report = SimpleNamespace(
        roots=(direct, derived),
        entries=(
            SimpleNamespace(
                reference=direct,
                disposition=AuditDisposition.CURRENT,
                dependencies=(),
                caveats=(),
            ),
            SimpleNamespace(
                reference=derived,
                disposition=AuditDisposition.BLOCKED,
                dependencies=(missing,),
                reason="Missing dependency",
                caveats=(),
            ),
        ),
        is_current=False,
    )

    audit_registry_freshness.print_report(report, "build.yaml")

    output = capsys.readouterr().out
    assert "Shared source pins cannot be fully checked" in output
    assert "UNCHECKED dependency lineage: example:output:1" in output
    assert "Pin consistency: 0 matches · 0 mismatches · " in output
    assert "0 same-source dataset reviews · 1 unchecked derived roots" in output


def test_unchecked_lineage_prevents_a_successful_audit(monkeypatch) -> None:
    derived = SnapshotRef.derived("example:output:1")
    report = SimpleNamespace(
        roots=(derived,),
        entries=(
            SimpleNamespace(
                reference=derived,
                disposition=AuditDisposition.CURRENT,
                dependencies=(),
                caveats=(),
            ),
        ),
        is_current=True,
    )
    monkeypatch.setattr(
        audit_registry_freshness, "audit_build_yaml", lambda *args, **kwargs: report
    )

    assert audit_registry_freshness.main(["build.yaml"]) == 1
