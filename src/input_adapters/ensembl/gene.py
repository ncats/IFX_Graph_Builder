"""Gene annotations from a pinned Ensembl human FTP GFF3 snapshot."""

import gzip
from typing import Generator, List
from urllib.parse import unquote

from src.constants import DataSourceName, Prefix
from src.interfaces.input_adapter import InputAdapter
from src.models.datasource_version_info import DatasourceVersionInfo
from src.models.gene import Gene, GeneticLocation, Strand
from src.models.node import EquivalentId


GFF3_FILE = "ensembl_homo_sapiens.gff3.gz"
GENE_FEATURE_TYPES = {"gene", "ncRNA_gene", "pseudogene"}
CANONICAL_CHROMOSOMES = {str(number) for number in range(1, 23)} | {"X", "Y", "MT"}


def gene_features(data_source, *, canonical_only=True):
    path = data_source.file(GFF3_FILE)
    with gzip.open(path, "rt", encoding="utf-8") as handle:
        first = handle.readline()
        if not first.startswith("##gff-version 3"):
            raise ValueError(f"Ensembl GFF3 header missing in {path}")
        for line in handle:
            if line.startswith("#"):
                continue
            columns = line.rstrip("\n").split("\t")
            if len(columns) != 9:
                raise ValueError(f"Malformed Ensembl GFF3 row in {path}: {line[:160]}")
            if columns[2] not in GENE_FEATURE_TYPES:
                continue
            if canonical_only and columns[0] not in CANONICAL_CHROMOSOMES:
                continue
            attributes = {
                key: unquote(value)
                for field in columns[8].split(";") if "=" in field
                for key, value in [field.split("=", 1)]
            }
            gene_id = attributes.get("gene_id")
            if not gene_id or attributes.get("ID") != f"gene:{gene_id}":
                raise ValueError(f"Ensembl GFF3 gene identifier is missing or inconsistent: {line[:160]}")
            yield gene_id, columns, attributes


class GeneAdapter(InputAdapter):
    def __init__(self, data_source, canonical_only: bool = True):
        self.data_source = data_source
        self.canonical_only = canonical_only

    def get_name(self) -> str:
        name = super().get_name()
        return name if self.canonical_only else f"{name} [all:{self.data_source.snapshot_id}]"

    def get_all(self) -> Generator[List[Gene], None, None]:
        genes = []
        seen = set()
        for gene_id, columns, attributes in gene_features(
            self.data_source, canonical_only=self.canonical_only
        ):
            if gene_id in seen:
                raise ValueError(f"Duplicate Ensembl GFF3 gene ID: {gene_id}")
            seen.add(gene_id)
            genes.append(Gene(
                id=EquivalentId(id=gene_id, type=Prefix.ENSEMBL).id_str(),
                location=GeneticLocation(
                    chromosome=columns[0],
                    start=int(columns[3]),
                    end=int(columns[4]),
                    strand=Strand.parse(columns[6]),
                ),
                symbol=attributes.get("Name"),
                biotype=attributes.get("biotype"),
            ))
        if not genes:
            raise ValueError(f"Ensembl GFF3 has no selected human gene features: {self.data_source.snapshot_id}")
        yield genes

    def get_datasource_name(self) -> DataSourceName:
        return DataSourceName.ENSEMBL

    def get_version(self) -> DatasourceVersionInfo:
        return self.data_source.version_info()
