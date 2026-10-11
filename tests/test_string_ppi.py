import gzip
from datetime import date

from src.input_adapters.string.string_ppi import StringPPIAdapter
from tests.registry_fakes import registry_dataset


def test_string_ppi_adapter_applies_default_cutoff_and_skips_self_pairs(tmp_path):
    data_path = tmp_path / "9606.protein.links.v12.5.txt.gz"

    with gzip.open(data_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "protein1 protein2 combined_score",
                    "9606.ENSP0001 9606.ENSP0002 399",
                    "9606.ENSP0001 9606.ENSP0001 900",
                    "9606.ENSP0002 9606.ENSP0003 400",
                    "9606.ENSP0003 9606.ENSP0004 700",
                ]
            )
        )

    adapter = StringPPIAdapter(
        data_source=registry_dataset(
            tmp_path,
            data_path.name,
            version="12.5",
            version_date=date(2025, 1, 15),
        ),
    )

    batches = list(adapter.get_all())
    edges = [edge for batch in batches for edge in batch]

    assert len(edges) == 2
    assert [(edge.start_node.id, edge.end_node.id, edge.score) for edge in edges] == [
        ("ENSEMBL:ENSP0002", "ENSEMBL:ENSP0003", [400]),
        ("ENSEMBL:ENSP0003", "ENSEMBL:ENSP0004", [700]),
    ]
    assert all(edge.sources == [] for edge in edges)

    version = adapter.get_version()
    assert version.version == "12.5"
    assert version.version_date.isoformat() == "2025-01-15"


def test_string_ppi_adapter_supports_explicit_cutoff_override(tmp_path):
    data_path = tmp_path / "9606.protein.links.v12.0.txt.gz"

    with gzip.open(data_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "protein1 protein2 combined_score",
                    "9606.ENSP0001 9606.ENSP0002 250",
                    "9606.ENSP0002 9606.ENSP0003 300",
                ]
            ),
        )

    adapter = StringPPIAdapter(
        data_source=registry_dataset(tmp_path, data_path.name),
        score_cutoff=300,
    )

    batches = list(adapter.get_all())
    edges = [edge for batch in batches for edge in batch]

    assert len(edges) == 1
    assert edges[0].start_node.id == "ENSEMBL:ENSP0002"
    assert edges[0].end_node.id == "ENSEMBL:ENSP0003"
    assert edges[0].score == [300]


def test_string_ppi_adapter_honors_max_rows_on_kept_edges(tmp_path):
    data_path = tmp_path / "9606.protein.links.v12.0.txt.gz"

    with gzip.open(data_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "protein1 protein2 combined_score",
                    "9606.ENSP0001 9606.ENSP0002 250",
                    "9606.ENSP0002 9606.ENSP0002 900",
                    "9606.ENSP0003 9606.ENSP0004 400",
                    "9606.ENSP0004 9606.ENSP0005 500",
                    "9606.ENSP0005 9606.ENSP0006 600",
                ]
            ),
        )

    adapter = StringPPIAdapter(
        data_source=registry_dataset(tmp_path, data_path.name),
        max_rows=2,
    )

    batches = list(adapter.get_all())
    edges = [edge for batch in batches for edge in batch]

    assert len(edges) == 2
    assert [(edge.start_node.id, edge.end_node.id, edge.score) for edge in edges] == [
        ("ENSEMBL:ENSP0003", "ENSEMBL:ENSP0004", [400]),
        ("ENSEMBL:ENSP0004", "ENSEMBL:ENSP0005", [500]),
    ]


def test_string_ppi_adapter_canonicalizes_pair_direction(tmp_path):
    data_path = tmp_path / "9606.protein.links.v12.0.txt.gz"

    with gzip.open(data_path, "wt", encoding="utf-8") as handle:
        handle.write(
            "\n".join(
                [
                    "protein1 protein2 combined_score",
                    "9606.ENSP9999 9606.ENSP0001 500",
                ]
            ),
        )

    adapter = StringPPIAdapter(registry_dataset(tmp_path, data_path.name))

    batches = list(adapter.get_all())
    edges = [edge for batch in batches for edge in batch]

    assert len(edges) == 1
    assert edges[0].start_node.id == "ENSEMBL:ENSP0001"
    assert edges[0].end_node.id == "ENSEMBL:ENSP9999"
    assert edges[0].score == [500]
