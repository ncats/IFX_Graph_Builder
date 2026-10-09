from concurrent.futures import ThreadPoolExecutor

from src.constants import Prefix
from src.id_resolvers.sqlite_cache_resolver import MatchingPair, SqliteCacheResolver
from src.id_resolvers.target_graph_resolver import TCRDTargetResolver
from src.models.node import EquivalentId


class _ThreadReusableResolver(SqliteCacheResolver):
    def __init__(self, cache_path, **kwargs):
        self.cache_path = str(cache_path)
        super().__init__(**kwargs)

    def cache_location(self):
        return self.cache_path

    def get_version_info(self) -> str:
        return "thread-reuse-test-v1"

    def matching_ids(self):
        yield MatchingPair("IFXProtein:1", "UniProtKB:P1", "UniProtKB")
        yield MatchingPair("IFXProtein:1", "Symbol:ABC1", "Symbol")


class _ProteinParser:
    def __init__(self, rows):
        self.rows = rows
        self.fieldnames = set()

    def all_rows(self):
        yield from self.rows

    @staticmethod
    def get_uniprot_reviewed(row):
        return row.get("reviewed", True)

    @staticmethod
    def get_uniprot_id(row):
        return row.get("uniprot_id")

    @staticmethod
    def get_is_canonical(row):
        return row.get("is_canonical")

    @staticmethod
    def get_id(row):
        return row["id"]

    @staticmethod
    def get_isoform_id(row):
        return row.get("canonical_ifx_id")

    @staticmethod
    def get_equivalent_ids(row):
        ids = []
        if row.get("uniprot_id"):
            ids.append(EquivalentId(id=row["uniprot_id"], type=Prefix.UniProtKB))
        if row.get("symbol"):
            ids.append(EquivalentId(id=row["symbol"], type=Prefix.Symbol))
        return ids

    @staticmethod
    def get_transcript_ids(row):
        return row.get("transcripts", [])

    @staticmethod
    def get_gene_id(row):
        return row.get("gene_id")

    @staticmethod
    def get_parent_gene_id(row):
        return None


class _TranscriptParser:
    fieldnames = set()
    rows = [
        {"id": "IFXTranscript:1", "ensembl": "ENST000001", "ncbi": "100"},
        {"id": "IFXTranscript:2", "ensembl": "ENST000002", "ncbi": "200"},
    ]

    def all_rows(self):
        yield from self.rows

    @staticmethod
    def get_id(row):
        return row["id"]

    @staticmethod
    def get_associated_ensg_id(row):
        return None

    @staticmethod
    def get_associated_ncbi_id(row):
        return row["ncbi"]

    @staticmethod
    def get_equivalent_ids(row):
        return [EquivalentId(id=row["ensembl"], type=Prefix.ENSEMBL)]


class _GeneParser:
    rows = [
        {"id": "IFXGene:1", "ncbi": "100", "symbol": "GENE1"},
        {"id": "IFXGene:2", "ncbi": "200", "symbol": "GENE2"},
    ]

    def all_rows(self):
        yield from self.rows

    @staticmethod
    def get_id(row):
        return row["id"]

    @staticmethod
    def get_equivalent_ids(row):
        return [
            EquivalentId(id=row["ncbi"], type=Prefix.NCBIGene),
            EquivalentId(id=row["symbol"], type=Prefix.Symbol),
        ]


def _resolver_with_rows(rows):
    resolver = object.__new__(TCRDTargetResolver)
    resolver.protein_parsers = [_ProteinParser(rows)]
    resolver.transcript_parser = _TranscriptParser()
    resolver.gene_parser = _GeneParser()
    resolver.reviewed_only = False
    resolver.collapse_to_canonical = True
    return resolver


def test_sqlite_cache_resolver_can_be_reused_from_worker_thread(tmp_path):
    resolver = _ThreadReusableResolver(tmp_path / "resolver.sqlite", types=["Protein"])

    with ThreadPoolExecutor(max_workers=1) as executor:
        prefix_counts = executor.submit(resolver.get_prefix_counts).result()

    busy_timeout = resolver.connection.execute("PRAGMA busy_timeout").fetchone()[0]

    assert busy_timeout == resolver.sqlite_busy_timeout_ms
    assert {"prefix": "UniProtKB", "count": 1} in prefix_counts
    assert {"prefix": "Symbol", "count": 1} in prefix_counts


def test_tcrd_target_resolver_maps_isoform_aliases_to_canonical_target():
    resolver = _resolver_with_rows([
        {
            "id": "IFXProtein:CANONICAL",
            "is_canonical": True,
            "uniprot_id": "P00001",
            "symbol": "CANON",
            "transcripts": ["ENST000001"],
            "gene_id": "100",
        },
        {
            "id": "IFXProtein:ISOFORM",
            "is_canonical": False,
            "canonical_ifx_id": "IFXProtein:CANONICAL",
            "uniprot_id": "P99999-2",
            "symbol": "ISO",
            "transcripts": ["ENST000002"],
            "gene_id": "200",
        },
    ])

    matches = set(resolver.matching_ids())

    assert MatchingPair("IFXProtein:CANONICAL", "IFXProtein:CANONICAL", "exact") in matches
    assert MatchingPair("IFXProtein:CANONICAL", "IFXProtein:ISOFORM", "isoform") in matches
    assert MatchingPair("IFXProtein:CANONICAL", "UniProtKB:P99999-2", "UniProtKB") in matches
    assert MatchingPair("IFXProtein:CANONICAL", "ENSEMBL:ENST000002", "ENSEMBL") in matches
    assert MatchingPair("IFXProtein:CANONICAL", "NCBIGene:100", "NCBIGene") in matches
    assert MatchingPair("IFXProtein:CANONICAL", "Symbol:GENE1", "Symbol") in matches
    assert MatchingPair("IFXProtein:CANONICAL", "NCBIGene:200", "NCBIGene") not in matches
    assert MatchingPair("IFXProtein:CANONICAL", "Symbol:GENE2", "Symbol") not in matches


def test_tcrd_target_resolver_drops_isoforms_without_canonical_target():
    resolver = _resolver_with_rows([
        {
            "id": "IFXProtein:CANONICAL",
            "is_canonical": True,
            "uniprot_id": "P00001",
        },
        {
            "id": "IFXProtein:ORPHAN",
            "is_canonical": False,
            "uniprot_id": "P99998-2",
            "symbol": "ORPHAN",
            "transcripts": ["ENST000002"],
            "gene_id": "200",
        },
    ])

    matches = set(resolver.matching_ids())

    assert MatchingPair("IFXProtein:CANONICAL", "IFXProtein:ORPHAN", "isoform") not in matches
    assert MatchingPair("IFXProtein:CANONICAL", "UniProtKB:P99998-2", "UniProtKB") not in matches
    assert MatchingPair("IFXProtein:CANONICAL", "ENSEMBL:ENST000002", "ENSEMBL") not in matches
    assert MatchingPair("IFXProtein:CANONICAL", "NCBIGene:200", "NCBIGene") not in matches
