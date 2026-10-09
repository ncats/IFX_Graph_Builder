"""Audit every IFX Registry pin used by an ODIN build YAML."""

from __future__ import annotations

import argparse
from collections.abc import Mapping
from dataclasses import dataclass
from datetime import timedelta
from pathlib import Path
import sys
from typing import Any, Literal

from ifx_registry import AuditDisposition, RegistryAuditClient, RegistryError, SnapshotRef
import yaml

from src.core.registry_integration import RegistryReference, RegistryReferenceKind


_REFERENCE_FACTORIES = {
    RegistryReferenceKind.SOURCE: SnapshotRef.source,
    RegistryReferenceKind.DERIVED: SnapshotRef.derived,
    RegistryReferenceKind.EXTERNAL: SnapshotRef.external,
}

@dataclass(frozen=True)
class SharedSourcePin:
    derived_root: SnapshotRef
    direct: SnapshotRef
    dependency: SnapshotRef | None
    other_inputs: tuple[SnapshotRef, ...]
    result: Literal["MATCH", "MISMATCH", "REVIEW"]


def collect_registry_references(
    config: Mapping[str, Any],
) -> tuple[SnapshotRef, ...]:
    """Collect the exact Registry roots that ODIN will resolve for a build."""
    references: list[SnapshotRef] = []
    loaded_files: set[Path] = set()

    def visit(node: Any, parent_key: str | None = None) -> None:
        if isinstance(node, Mapping):
            for key, value in node.items():
                if key == "data_source" or str(key).endswith("_data_source"):
                    parsed = RegistryReference.parse(value)
                    references.append(
                        _REFERENCE_FACTORIES[parsed.kind](parsed.snapshot_id)
                    )
                else:
                    visit(value, str(key))
        elif isinstance(node, list):
            for value in node:
                visit(value, parent_key)
        elif (
            isinstance(node, str)
            and (node.endswith(".yaml") or node.endswith(".yml"))
            and not _is_credentials_key(parent_key)
        ):
            # Match Config exactly: nested YAML paths are interpreted from cwd.
            nested_path = Path(node)
            resolved = nested_path.resolve()
            if resolved in loaded_files:
                raise ValueError(f"Nested YAML cycle detected at {nested_path}")
            loaded_files.add(resolved)
            with nested_path.open(encoding="utf-8") as stream:
                nested = yaml.safe_load(stream)
            visit(nested, parent_key)
            loaded_files.remove(resolved)

    # These are the same sections whose Registry references Config resolves.
    for section in ("resolvers", "input_adapters"):
        if section in config:
            visit(config[section], section)
    return tuple(dict.fromkeys(references))


def audit_build_yaml(
    yaml_path: str | Path,
    *,
    timeout: timedelta = timedelta(seconds=30),
) -> Any:
    """Load one build YAML and return Registry's dependency-complete audit."""
    path = Path(yaml_path)
    with path.open(encoding="utf-8") as stream:
        config = yaml.safe_load(stream)
    if not isinstance(config, Mapping):
        raise TypeError(f"Build YAML must contain a mapping: {path}")

    registry_config = config.get("registry")
    if not isinstance(registry_config, Mapping) or not registry_config:
        raise ValueError(f"Build YAML has no registry configuration: {path}")
    references = collect_registry_references(config)
    if not references:
        return None

    client = RegistryAuditClient.connect(
        registry_config.get("credentials"),
        bucket=_optional_string(registry_config.get("bucket")),
        region=_optional_string(registry_config.get("region")),
        prefix=str(registry_config.get("prefix", "")),
        source_configuration=registry_config.get("source_configuration"),
        cure_credentials_file=registry_config.get("cure_credentials"),
    )
    return client.audit(references, timeout=timeout)


def compare_shared_source_pins(report: Any) -> tuple[SharedSourcePin, ...]:
    """Compare direct source pins with source inputs of each derived YAML root."""
    if report is None:
        return ()
    roots = tuple(getattr(report, "roots", ()))
    direct_sources = sorted(
        (ref for ref in roots if ref.kind == RegistryReferenceKind.SOURCE),
        key=lambda ref: ref.snapshot_id,
    )
    entries = {entry.reference: entry for entry in report.entries}
    comparisons: list[SharedSourcePin] = []
    for derived_root in sorted(
        (ref for ref in roots if ref.kind == RegistryReferenceKind.DERIVED),
        key=lambda ref: ref.snapshot_id,
    ):
        dependencies = _transitive_source_dependencies(derived_root, entries)
        for direct in direct_sources:
            shared_source = tuple(
                dependency
                for dependency in dependencies
                if dependency.dataset.source == direct.dataset.source
            )
            if not shared_source:
                continue
            same_dataset = tuple(
                dependency
                for dependency in shared_source
                if dependency.dataset == direct.dataset
            )
            other_inputs = tuple(
                dependency
                for dependency in shared_source
                if dependency.dataset != direct.dataset
            )
            if not same_dataset:
                comparisons.append(
                    SharedSourcePin(
                        derived_root, direct, None, other_inputs, "REVIEW"
                    )
                )
                continue
            for dependency in same_dataset:
                comparisons.append(
                    SharedSourcePin(
                        derived_root,
                        direct,
                        dependency,
                        other_inputs,
                        "MATCH" if direct.version == dependency.version else "MISMATCH",
                    )
                )
    return tuple(comparisons)


def _unchecked_derived_roots(report: Any) -> tuple[SnapshotRef, ...]:
    if report is None:
        return ()
    entries = {entry.reference: entry for entry in report.entries}
    unchecked: set[SnapshotRef] = set()
    for root in getattr(report, "roots", ()):
        if root.kind != RegistryReferenceKind.DERIVED:
            continue
        pending = [root]
        visited: set[SnapshotRef] = set()
        while pending:
            reference = pending.pop()
            if reference in visited:
                continue
            visited.add(reference)
            if reference.kind != RegistryReferenceKind.DERIVED:
                continue
            entry = entries.get(reference)
            dependencies = getattr(entry, "dependencies", ()) if entry else ()
            if not dependencies:
                unchecked.add(root)
                break
            pending.extend(dependencies)
    return tuple(sorted(unchecked, key=lambda ref: ref.snapshot_id))


def _transitive_source_dependencies(
    derived_root: SnapshotRef, entries: Mapping[SnapshotRef, Any]
) -> tuple[SnapshotRef, ...]:
    pending = [derived_root]
    visited: set[SnapshotRef] = set()
    sources: set[SnapshotRef] = set()
    while pending:
        reference = pending.pop()
        if reference in visited:
            continue
        visited.add(reference)
        if reference.kind == RegistryReferenceKind.SOURCE:
            sources.add(reference)
            continue
        entry = entries.get(reference)
        if entry is not None:
            pending.extend(getattr(entry, "dependencies", ()))
    return tuple(sorted(sources, key=lambda ref: ref.snapshot_id))


def print_report(report: Any, yaml_path: str | Path) -> None:
    """Print a dependency-first, operator-oriented audit report."""
    print(f"Registry freshness for {yaml_path}")
    if report is None:
        print("No Registry references found in resolvers or input adapters.")
        return
    roots = frozenset(
        getattr(report, "roots", tuple(entry.reference for entry in report.entries))
    )
    rebuild_dependencies = _rebuild_dependency_references(report.entries)

    groups = (
        (
            "Refresh Registry, then update YAML",
            lambda entry: (
                entry.disposition is AuditDisposition.REGISTER_SOURCE
                and entry.reference in roots
                and _registration_target(entry) is not None
                and _registration_target(entry) != entry.reference.snapshot_id
            ),
            _print_registration_entry,
        ),
        (
            "Refresh Registry",
            lambda entry: (
                entry.disposition is AuditDisposition.REGISTER_SOURCE
                and (
                    entry.reference not in roots
                    or _registration_target(entry) is None
                    or _registration_target(entry) == entry.reference.snapshot_id
                )
            ),
            _print_registration_entry,
        ),
        (
            "Rebuild in Registry, then update YAML",
            lambda entry: (
                entry.disposition is AuditDisposition.REBUILD_DERIVED
                and entry.reference in roots
            ),
            _print_rebuild_entry,
        ),
        (
            "Rebuild in Registry",
            lambda entry: (
                entry.disposition is AuditDisposition.REBUILD_DERIVED
                and entry.reference not in roots
            ),
            _print_rebuild_entry,
        ),
        (
            "Then rebuild dependent datasets",
            _is_dependent_rebuild,
            _print_reason_entry,
        ),
        (
            "Update YAML to newest registered pin",
            lambda entry: (
                entry.disposition is AuditDisposition.UPDATE_PIN
                and entry.reference in roots
            ),
            _print_pin_update_entry,
        ),
        (
            "Use newest registered inputs when rebuilding",
            lambda entry: (
                entry.disposition is AuditDisposition.UPDATE_PIN
                and entry.reference not in roots
                and entry.reference in rebuild_dependencies
            ),
            _print_pin_update_entry,
        ),
        (
            "Needs manual check",
            lambda entry: entry.disposition is AuditDisposition.UNVERIFIABLE,
            _print_reason_entry,
        ),
        (
            "Up to date — manual source check needed",
            lambda entry: (
                entry.disposition is AuditDisposition.CURRENT
                and _has_only_manual_caveats(entry)
            ),
            _print_reason_entry,
        ),
        (
            "Registered pin — freshness caveats",
            lambda entry: (
                entry.disposition is AuditDisposition.CURRENT
                and bool(_caveats(entry))
                and not _has_only_manual_caveats(entry)
            ),
            _print_reason_entry,
        ),
        (
            "Blocked",
            lambda entry: (
                entry.disposition is AuditDisposition.BLOCKED
                and not _is_dependent_rebuild(entry)
            ),
            _print_reason_entry,
        ),
        (
            "Up to date",
            lambda entry: (
                entry.disposition is AuditDisposition.CURRENT
                and not _caveats(entry)
            ),
            _print_current_entry,
        ),
    )
    for heading, predicate, printer in groups:
        entries = [entry for entry in report.entries if predicate(entry)]
        if not entries:
            continue
        print(f"\n{heading} ({len(entries)})")
        for entry in entries:
            printer(entry)

    comparisons = compare_shared_source_pins(report)
    unchecked_roots = _unchecked_derived_roots(report)
    print(f"\nShared source pins ({len(comparisons)} comparisons)")
    if not comparisons:
        if unchecked_roots:
            print("  Shared source pins cannot be fully checked.")
        else:
            print("  No source provider is shared by a direct YAML pin and a derived input.")
    else:
        for comparison in comparisons:
            label = (
                "REVIEW (same source, different datasets)"
                if comparison.result == "REVIEW"
                else comparison.result
            )
            print(f"  {label}: {comparison.direct.dataset}")
            print(f"    YAML: {comparison.direct.snapshot_id}")
            print(
                "    Derived: "
                + (
                    comparison.dependency.snapshot_id
                    if comparison.dependency is not None
                    else f"no {comparison.direct.dataset} pin"
                )
                + f" (via {comparison.derived_root.snapshot_id})"
            )
            if comparison.other_inputs:
                print(f"    Other {comparison.direct.dataset.source} inputs:")
                for other in comparison.other_inputs:
                    print(f"      {other.snapshot_id}")
            if comparison.result == "REVIEW":
                print(
                    "    Review whether these datasets should share a pin "
                    "or remain separately versioned."
                )
    for root in unchecked_roots:
        print(f"  UNCHECKED dependency lineage: {root.snapshot_id}")

    known_update_count = sum(
        entry.disposition
        in {AuditDisposition.REGISTER_SOURCE, AuditDisposition.REBUILD_DERIVED}
        or (
            entry.disposition is AuditDisposition.UPDATE_PIN
            and entry.reference in roots
        )
        for entry in report.entries
    )
    dependent_rebuild_count = sum(_is_dependent_rebuild(entry) for entry in report.entries)
    manual_sources = {
        caveat.origin.snapshot_id
        for entry in report.entries
        for caveat in _caveats(entry)
        if _caveat_code(caveat) == "manual_freshness"
    }
    qualified_current_count = sum(
        entry.disposition is AuditDisposition.CURRENT and bool(_caveats(entry))
        for entry in report.entries
    )
    fully_current_count = sum(
        entry.disposition is AuditDisposition.CURRENT and not _caveats(entry)
        for entry in report.entries
    )
    blocked_count = sum(
        entry.disposition is AuditDisposition.BLOCKED
        and not _is_dependent_rebuild(entry)
        for entry in report.entries
    )
    print(
        "\nResult: "
        f"{known_update_count} known updates · "
        f"{dependent_rebuild_count} dependent rebuilds · "
        f"{len(manual_sources)} manual sources to confirm · "
        f"{qualified_current_count} registered with freshness caveats · "
        f"{blocked_count} blocked · "
        f"{fully_current_count} fully up to date"
    )
    match_count = sum(item.result == "MATCH" for item in comparisons)
    mismatch_count = sum(item.result == "MISMATCH" for item in comparisons)
    review_count = sum(item.result == "REVIEW" for item in comparisons)
    print(
        f"Pin consistency: {match_count} matches · {mismatch_count} mismatches · "
        f"{review_count} same-source dataset reviews · "
        f"{len(unchecked_roots)} unchecked derived roots"
    )


def _print_current_entry(entry: Any) -> None:
    print(f"  {entry.reference.snapshot_id}")


def _print_reason_entry(entry: Any) -> None:
    _print_current_entry(entry)
    print(f"    {entry.reason}")


def _print_pin_update_entry(entry: Any) -> None:
    replacement = entry.recommended_snapshot_id
    if replacement:
        print(f"  {entry.reference.snapshot_id} -> {replacement}")
        print(f"    {entry.reason}")
        _print_caveats(entry)
    else:
        _print_reason_entry(entry)


def _print_registration_entry(entry: Any) -> None:
    replacement = _registration_target(entry)
    if replacement is None:
        _print_reason_entry(entry)
        return
    if replacement == entry.reference.snapshot_id:
        print(f"  {replacement} (not registered yet)")
    else:
        print(f"  {entry.reference.snapshot_id} -> {replacement} (not registered yet)")
    _print_caveats(entry)


def _registration_target(entry: Any) -> str | None:
    upstream = entry.latest_upstream_version
    if upstream is None:
        return None
    return f"{entry.reference.dataset}:{upstream.value}"


def _print_rebuild_entry(entry: Any) -> None:
    if entry.recommended_snapshot_id:
        print(
            f"  {entry.reference.snapshot_id} -> "
            f"{entry.recommended_snapshot_id}"
        )
    else:
        _print_current_entry(entry)
    print(f"    {entry.reason}")
    _print_caveats(entry)


def _print_caveats(entry: Any) -> None:
    for caveat in _caveats(entry):
        print(f"    Caveat: {caveat.message}")


def _caveats(entry: Any) -> tuple[Any, ...]:
    return tuple(getattr(entry, "caveats", ()))


def _caveat_code(caveat: Any) -> str:
    code = caveat.code
    return str(getattr(code, "value", code))


def _has_only_manual_caveats(entry: Any) -> bool:
    caveats = _caveats(entry)
    return bool(caveats) and all(
        _caveat_code(caveat) == "manual_freshness" for caveat in caveats
    )


def _rebuild_dependency_references(entries: tuple[Any, ...]) -> frozenset[SnapshotRef]:
    """Return only the dependency closure of rebuild actions in this report."""
    dependencies = {
        entry.reference: tuple(getattr(entry, "dependencies", ()))
        for entry in entries
    }
    pending = [
        entry.reference
        for entry in entries
        if entry.disposition is AuditDisposition.REBUILD_DERIVED
    ]
    result: set[SnapshotRef] = set()
    while pending:
        reference = pending.pop()
        for dependency in dependencies.get(reference, ()):
            if dependency not in result:
                result.add(dependency)
                pending.append(dependency)
    return frozenset(result)


def _is_dependent_rebuild(entry: Any) -> bool:
    return (
        entry.disposition is AuditDisposition.BLOCKED
        and entry.reason.startswith(("After rebuilding ", "After registering "))
    )


def main(argv: list[str] | None = None) -> int:
    parser = argparse.ArgumentParser(
        description=(
            "Audit all Registry pins in an ODIN build YAML, including transitive "
            "derived-dataset dependencies."
        )
    )
    parser.add_argument("yaml", type=Path, help="ODIN build YAML to inspect")
    parser.add_argument(
        "--timeout-seconds",
        type=float,
        default=30,
        help="timeout for each upstream source check (default: 30)",
    )
    args = parser.parse_args(argv)
    if args.timeout_seconds <= 0:
        parser.error("--timeout-seconds must be positive")

    try:
        report = audit_build_yaml(
            args.yaml,
            timeout=timedelta(seconds=args.timeout_seconds),
        )
    except (OSError, TypeError, ValueError, RegistryError, yaml.YAMLError) as error:
        print(f"Freshness audit failed: {error}", file=sys.stderr)
        return 2

    print_report(report, args.yaml)
    has_pin_mismatch = any(
        comparison.result == "MISMATCH"
        for comparison in compare_shared_source_pins(report)
    )
    has_unchecked_lineage = bool(_unchecked_derived_roots(report))
    return (
        0
        if report is None
        or (report.is_current and not has_pin_mismatch and not has_unchecked_lineage)
        else 1
    )


def _optional_string(value: object) -> str | None:
    return str(value) if value is not None else None


def _is_credentials_key(key: str | None) -> bool:
    return key is not None and "credential" in key.lower()


if __name__ == "__main__":
    raise SystemExit(main())
