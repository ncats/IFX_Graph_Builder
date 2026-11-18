"""Resolve human gene identifiers from a pinned Ensembl FTP Registry snapshot."""

import csv
import gzip
import re
from collections.abc import Generator
from typing import Any

from src.constants import Prefix
from src.id_resolvers.sqlite_cache_resolver import MatchingPair, SqliteCacheResolver
from src.input_adapters.ensembl.gene import gene_features
from src.interfaces.id_resolver import IdMatch
from src.models.gene import Gene
from src.models.node import EquivalentId


ENTREZ_FILE = "ensembl_entrez.tsv.gz"
RESOLVER_CONTRACT_VERSION = "2"
_HGNC_ACCESSION = re.compile(r"\bAcc:HGNC:(\d+)\b")
_VERSIONED_GENE = re.compile(r"^(ENSEMBL:ENSG\d+)\.\d+$")


class EnsemblGeneResolver(SqliteCacheResolver):
    """Resolve source-asserted Ensembl, HGNC, NCBI Gene, and symbol IDs."""

    name = "Ensembl Gene Resolver"

    def __init__(self, data_source, **kwargs):
        self.data_source = data_source
        # Fail before opening a resolver cache if either required Registry file is absent.
        self.data_source.file("ensembl_homo_sapiens.gff3.gz")
        self.entrez_file = data_source.file(ENTREZ_FILE)
        super().__init__(**kwargs)

    def get_version_info(self) -> str:
        return (
            f"{self.data_source.snapshot_id}:"
            f"ensembl-gene-resolver-{RESOLVER_CONTRACT_VERSION}"
        )

    def resolve_internal(self, input_nodes):
        matches = super().resolve_internal(input_nodes)
        fallback = {}
        for node in input_nodes:
            if node.id in matches:
                continue
            versioned = _VERSIONED_GENE.fullmatch(node.id)
            if versioned:
                fallback[node.id] = versioned.group(1)
        if fallback:
            base_matches = super().resolve_internal([Gene(id=base) for base in set(fallback.values())])
            for source_id, base in fallback.items():
                if base not in base_matches:
                    continue
                matches[source_id] = [
                    IdMatch(
                        input=source_id,
                        match=result.match,
                        equivalent_ids=result.equivalent_ids,
                        context=[*(result.context or []), "version_free_fallback"],
                    )
                    for result in base_matches[base]
                ]
        return matches

    def matching_ids(self) -> Generator[MatchingPair, Any, None]:
        emitted: set[MatchingPair] = set()
        gene_ids = set()
        for gene_id, _, attributes in gene_features(self.data_source, canonical_only=False):
            ensembl_id = _equivalent(gene_id, Prefix.ENSEMBL)
            if ensembl_id is None:
                continue
            gene_ids.add(gene_id)
            pairs = [MatchingPair(id=ensembl_id, match=ensembl_id, type="exact")]
            for hgnc_number in _HGNC_ACCESSION.findall(attributes.get("description", "")):
                pairs.append(MatchingPair(
                    id=ensembl_id,
                    match=_equivalent(hgnc_number, Prefix.HGNC),
                    type=Prefix.HGNC.value,
                ))
            symbol = _equivalent(attributes.get("Name"), Prefix.Symbol)
            if symbol is not None:
                pairs.append(MatchingPair(id=ensembl_id, match=symbol, type=Prefix.Symbol.value))
            for pair in pairs:
                if pair not in emitted:
                    emitted.add(pair)
                    yield pair

        with gzip.open(self.entrez_file, "rt", encoding="utf-8", newline="") as handle:
            reader = csv.DictReader(handle, delimiter="\t")
            required = {"gene_stable_id", "xref", "db_name"}
            if not required.issubset(reader.fieldnames or []):
                raise ValueError(f"Ensembl Entrez file lacks required columns: {self.entrez_file}")
            for row in reader:
                gene_id = (row["gene_stable_id"] or "").strip()
                ncbi_id = (row["xref"] or "").strip()
                if row["db_name"] != "EntrezGene" or gene_id not in gene_ids or not ncbi_id.isdecimal():
                    continue
                pair = MatchingPair(
                    id=_equivalent(gene_id, Prefix.ENSEMBL),
                    match=_equivalent(str(int(ncbi_id)), Prefix.NCBIGene),
                    type=Prefix.NCBIGene.value,
                )
                if pair not in emitted:
                    emitted.add(pair)
                    yield pair


def _equivalent(raw_value: str | None, prefix: Prefix) -> str | None:
    value = (raw_value or "").strip()
    if not value:
        return None
    return EquivalentId(id=value, type=prefix).id_str()
