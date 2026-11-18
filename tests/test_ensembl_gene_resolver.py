import gzip
from pathlib import Path

import pytest

from src.constants import Prefix
from src.id_resolvers.ensembl_gene_resolver import EnsemblGeneResolver
from src.input_adapters.ensembl.gene import GeneAdapter
from src.models.node import EquivalentId
from tests.registry_fakes import registry_dataset


def _id(value: str, prefix: Prefix) -> str:
    return EquivalentId(id=value, type=prefix).id_str()


def _source(tmp_path: Path):
    gff3 = tmp_path / "ensembl_homo_sapiens.gff3.gz"
    with gzip.open(gff3, "wt") as handle:
        handle.write("##gff-version 3\n")
        for chrom, kind, gene_id, name in (
            ("17", "gene", "ENSG00000141510", "TP53"),
            ("1", "ncRNA_gene", "ENSG00000210049", "MT-TF"),
            ("2", "pseudogene", "ENSG00000999999", "PSEUDO"),
            ("CHR_PATCH", "gene", "ENSG00000888888", "PATCH"),
        ):
            handle.write(
                f"{chrom}\tensembl\t{kind}\t100\t200\t.\t+\t.\t"
                f"ID=gene:{gene_id};gene_id={gene_id};Name={name};biotype={kind};"
                "description=some gene [Source:HGNC Symbol%3BAcc:HGNC:11998]\n"
            )
    entrez = tmp_path / "ensembl_entrez.tsv.gz"
    with gzip.open(entrez, "wt") as handle:
        handle.write("gene_stable_id\ttranscript_stable_id\txref\tdb_name\n")
        handle.write("ENSG00000141510\tENST1\t7157\tEntrezGene\n")
        handle.write("ENSG00000141510\tENST2\t7157\tEntrezGene\n")
        handle.write("ENSG00000141510\tENST3\t9999\tEntrezGene\n")
        handle.write("ENSG00000141510\tENST4\tBAD\tEntrezGene\n")
        handle.write("ENSG00000141510\tENST5\t7777\tEntrezGene_trans_name\n")
    return registry_dataset(
        tmp_path, gff3.name, entrez.name,
        source="ensembl", dataset="human_ftp", version="116",
    )


def test_ensembl_resolver_uses_ftp_gene_and_entrez_assertions(tmp_path, monkeypatch):
    source = _source(tmp_path)
    cache_file = tmp_path / "ensembl-resolver.sqlite"
    monkeypatch.setattr(EnsemblGeneResolver, "cache_location", lambda self: str(cache_file))
    resolver = EnsemblGeneResolver(data_source=source, types=["Gene"])
    gene = _id("ENSG00000141510", Prefix.ENSEMBL)
    pairs = {(pair.id, pair.match, pair.type) for pair in resolver.matching_ids()}
    assert (gene, gene, "exact") in pairs
    assert (gene, _id("11998", Prefix.HGNC), Prefix.HGNC.value) in pairs
    assert (gene, _id("TP53", Prefix.Symbol), Prefix.Symbol.value) in pairs
    assert (gene, _id("7157", Prefix.NCBIGene), Prefix.NCBIGene.value) in pairs
    assert (gene, _id("9999", Prefix.NCBIGene), Prefix.NCBIGene.value) in pairs
    assert (gene, _id("7777", Prefix.NCBIGene), Prefix.NCBIGene.value) not in pairs
    for other in ("ENSG00000210049", "ENSG00000999999", "ENSG00000888888"):
        ensembl = _id(other, Prefix.ENSEMBL)
        assert (ensembl, ensembl, "exact") in pairs
    assert resolver.get_version_info() == "ensembl:human_ftp:116:ensembl-gene-resolver-2"


def test_ensembl_gene_adapter_reads_canonical_gene_like_features(tmp_path):
    adapter = GeneAdapter(data_source=_source(tmp_path))
    genes = [gene for batch in adapter.get_all() for gene in batch]
    assert len(genes) == 3
    assert {gene.id for gene in genes} == {
        _id("ENSG00000141510", Prefix.ENSEMBL),
        _id("ENSG00000210049", Prefix.ENSEMBL),
        _id("ENSG00000999999", Prefix.ENSEMBL),
    }
    assert genes[0].symbol == "TP53"
    assert genes[0].biotype == "gene"
    assert genes[0].location.chromosome == "17"
    assert (genes[0].location.start, genes[0].location.end) == (100, 200)
    assert adapter.get_version().version == "116"


def test_ensembl_resolver_requires_both_registry_files(tmp_path):
    source = registry_dataset(
        tmp_path, "ensembl_homo_sapiens.gff3.gz",
        source="ensembl", dataset="human_ftp", version="116",
    )
    with pytest.raises(FileNotFoundError, match="ensembl_homo_sapiens.gff3.gz"):
        EnsemblGeneResolver(data_source=source, types=["Gene"])

    (tmp_path / "ensembl_homo_sapiens.gff3.gz").write_bytes(b"placeholder")
    with pytest.raises(FileNotFoundError, match="ensembl_entrez.tsv.gz"):
        EnsemblGeneResolver(data_source=source, types=["Gene"])
