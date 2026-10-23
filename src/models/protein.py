from dataclasses import dataclass, asdict
from enum import Enum
from typing import Union, Optional, Dict, List

from src.core.decorators import facets, search
from src.interfaces.simple_enum import SimpleEnum
from src.models.analyte import Analyte
from src.models.gene import Audited
from src.models.node import Relationship, Node
from src.models.publication import PublicationReference
from src.models.reaction import Reaction
from src.models.year_score import YearScore


class TDL(SimpleEnum):
    Tclin = 'Tclin'
    Tchem = 'Tchem'
    Tbio = 'Tbio'
    Tdark = 'Tdark'

@dataclass
class TDLMetadata:
    tdl_ligand_count: int
    tdl_drug_count: int
    tdl_go_term_count: int
    tdl_generif_count: int
    tdl_pm_score: float
    tdl_antibody_count: int

    def to_dict(self) -> Dict[str, str]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict):
        fields = ['tdl_ligand_count', 'tdl_drug_count', 'tdl_go_term_count',
                  'tdl_generif_count', 'tdl_pm_score', 'tdl_antibody_count']
        found = False
        for field in fields:
            if field in data:
                found = True
        if not found:
            return None
        return TDLMetadata(
            tdl_ligand_count=data['tdl_ligand_count'],
            tdl_drug_count=data['tdl_drug_count'],
            tdl_go_term_count=data['tdl_go_term_count'],
            tdl_generif_count=data['tdl_generif_count'],
            tdl_pm_score=data['tdl_pm_score'],
            tdl_antibody_count=data['tdl_antibody_count']
        )

class IDGFamily(SimpleEnum):
    Other = 'Other'
    Enzyme = 'Enzyme'
    TranscriptionFactor = 'Transcription Factor'
    Kinase = 'Kinase'
    Transporter = 'Transporter'
    oGPCR = 'oGPCR'
    GPCR = 'GPCR'
    IonChannel = 'Ion Channel'
    Epigenetic = 'Epigenetic'
    NuclearReceptor = 'Nuclear Receptor'
    TFEpigenetic = 'TF-Epigenetic'

    @classmethod
    def parse(cls, tcrd_value: Union[str, Enum]):
        if tcrd_value is None or tcrd_value == '':
            return None
        if tcrd_value == 'IC':
            return IDGFamily.IonChannel
        if tcrd_value == 'TF; Epigenetic':
            return IDGFamily.TFEpigenetic
        if tcrd_value == 'TF':
            return IDGFamily.TranscriptionFactor
        if tcrd_value == 'NR':
            return IDGFamily.NuclearReceptor

        for member in IDGFamily:
            if member.value == tcrd_value:
                return member

        raise Exception(f"Unknown IDG Family: {tcrd_value}")

@dataclass
@facets(
    category_fields=["protein_type", "tdl", "idg_family", "uniprot_reviewed", "canonical_isoform_status"],
    numeric_fields=["pm_score", "pt_score", "antibody_count", "novelty"])
@search(text_fields=["name", "symbol", "preferred_symbol", "description", "gene_name", "uniprot_id", "ensembl_id", "refseq_id", "ncbi_id", "ncbi_gene_summary"])
class Protein(Audited, Analyte):
    target_payload: Optional[Dict[str, str]] = None
    preferred_symbol: Optional[str] = None
    protein_type: Optional[str] = None
    description: Optional[str] = None
    symbol: Optional[str] = None
    antibody_count: Optional[List[int]] = None
    ncbi_id: Optional[str] = None
    ensembl_id: Optional[str] = None
    refseq_id: Optional[str] = None
    uniprot_id: Optional[str] = None
    secondary_uniprot_ids: Optional[List[str]] = None
    sequence: Optional[str] = None
    gene_name: Optional[str] = None
    tdl: Optional[TDL] = None
    tdl_meta: Optional[TDLMetadata] = None
    name: Optional[str] = None
    idg_family: Optional[IDGFamily] = None
    pm_score: Optional[List[float]] = None
    pm_score_by_year: Optional[List[YearScore]] = None
    pt_score: Optional[List[float]] = None
    pt_score_by_year: Optional[List[YearScore]] = None
    patent_family_mentions: Optional[List[str]] = None
    patent_identifier_sources: Optional[List[str]] = None
    novelty: Optional[List[float]] = None
    publications: Optional[List[PublicationReference]] = None
    ncbi_gene_summary: Optional[str] = None
    uniprot_annotationScore: Optional[int] = None
    uniprot_reviewed: Optional[bool] = None
    uniprot_function: Optional[str] = None
    similarity: Optional[str] = None
    Ensembl_ID_Provenance: Optional[str] = None
    RefSeq_ID_Provenance: Optional[str] = None
    Uniprot_ID_Provenance: Optional[str] = None
    ensembl_canonical: Optional[str] = None
    is_canonical: Optional[bool] = None
    uniprot_canonical: Optional[bool] = None
    canonical_isoform_status: Optional[str] = None
    uniprot_entryType: Optional[str] = None
    mapping_ratio: Optional[float] = None
    protein_name_score: Optional[str] = None
    protein_name_method: Optional[str] = None
    uniprot_isoform: Optional[str] = None
    calculated_properties: Optional[Dict[str, float]] = None
    dtoid: Optional[str] = None
    dtoclass: Optional[str] = None


@dataclass
class ProteinReactionEdge(Relationship):
    start_node: Protein
    end_node: Reaction
    is_reviewed: bool = None


@dataclass
@facets(category_fields=['type'])
@search(text_fields=['name', 'type'])
class MeasuredProtein(Node):
    name: Optional[str] = None
    type: Optional[str] = None


@dataclass
class MeasuredProteinEdge(Relationship):
    start_node: MeasuredProtein
    end_node: Protein
