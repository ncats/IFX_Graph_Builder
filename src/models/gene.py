from dataclasses import dataclass
from datetime import datetime
from typing import Dict, List, Optional

from src.core.decorators import facets, search
from src.interfaces.simple_enum import SimpleEnum
from src.models.node import Node, Relationship
from src.models.publication import PublicationReference
from src.models.year_score import YearScore


class Strand(SimpleEnum):
    plus1 = "+1"
    minus1 = "-1"

    @classmethod
    def parse(cls, input_value: str):
        if input_value is None or input_value == '':
            return None
        if input_value == '+':
            return Strand.plus1
        if input_value == '-':
            return Strand.minus1
        int_val = int(float(input_value))
        if int_val > 0:
            return Strand.plus1
        return Strand.minus1

@dataclass
class GeneticLocation:
    location: Optional[str] = None
    chromosome: Optional[str] = None
    strand: Optional[Strand] = None
    start: Optional[int] = None
    end: Optional[int] = None

    def to_dict(self) -> Dict[str, str]:
        ret_dict = {}
        if self.location is not None:
            ret_dict['location'] = self.location
        if self.chromosome is not None:
            ret_dict['chromosome'] = self.chromosome
        if self.strand is not None:
            ret_dict['chromosome_strand'] = self.strand.value
        if self.start is not None:
            ret_dict['start_position'] = self.start
        if self.end is not None:
            ret_dict['end_position'] = self.end
        return ret_dict

    @classmethod
    def from_dict(cls, data: dict):
        if data is None:
            return None
        return GeneticLocation(
            location=data.get('location'),
            chromosome=data.get('chromosome'),
            start=data.get('start_position'),
            end=data.get('end_position'),
            strand=Strand.parse(data.get('chromosome_strand')))


@dataclass
class Audited:
    created: Optional[datetime] = None
    updated: Optional[datetime] = None


@dataclass
@facets(category_fields=["symbol"])
@search(text_fields=["symbol", "ncbi_gene_summary"])
class Gene(Audited, Node):
    target_payload: Optional[Dict[str, str]] = None
    location: Optional[GeneticLocation] = None
    publications: Optional[List[PublicationReference]] = None
    ncbi_gene_summary: Optional[str] = None
    pt_score: Optional[List[float]] = None
    pt_score_by_year: Optional[List[YearScore]] = None
    mapping_ratio: Optional[float] = None
    symbol: Optional[str] = None
    biotype: Optional[str] = None
    calculated_properties: Optional[Dict[str, float]] = None
    Name_Provenance: Optional[str] = None
    Location_Provenance: Optional[str] = None
    Ensembl_ID_Provenance: Optional[str] = None
    NCBI_ID_Provenance: Optional[str] = None
    HGNC_ID_Provenance: Optional[str] = None
    Symbol_Provenance: Optional[str] = None

@dataclass
@facets(category_fields=['type'])
@search(text_fields=['symbol', 'biotype'])
class MeasuredGene(Node):
    symbol: Optional[str] = None
    biotype: Optional[str] = None

@dataclass
class MeasuredGeneEdge(Relationship):
    start_node: MeasuredGene
    end_node: Gene
