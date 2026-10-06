"""Stable identity of an ODIN source build, shared by stages and exporters."""
from src.core.curations import payload_sha256


def source_build_fingerprint(metadata: dict):
    """Hash build identity and pinned inputs, independent of database transport."""
    registry_datasets = sorted([
        {
            "snapshot_id": item.get("snapshot_id"),
            "build_key": item.get("build_key"),
            "publication_fingerprint": item.get("publication_fingerprint"),
        }
        for item in metadata.get("registry_datasets") or []
    ], key=lambda item: (
        item.get("snapshot_id") or "",
        item.get("build_key") or "",
        item.get("publication_fingerprint") or "",
    ))
    identity = {
        "run_key": metadata.get("_key"),
        "run_date": metadata.get("run_date"),
        "source_yaml": metadata.get("source_yaml"),
        "git_commit": (metadata.get("git_info") or {}).get("git_commit"),
        "registry_datasets": registry_datasets,
    }
    if not any(value for value in identity.values()):
        return None
    return payload_sha256(identity)
