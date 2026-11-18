from dataclasses import dataclass, field
from datetime import date
from typing import Dict, Optional, List, TYPE_CHECKING

from src.core.decorators import facets
from src.interfaces.simple_enum import SimpleEnum
from src.models.node import Node, Relationship
from src.models.pounce.project import AccessLevel
from src.models.pounce.workbook_artifact import WorkbookArtifact


class PlatformType(SimpleEnum):
    BULK_RNASEQ  = "bulk_rnaseq"
    SC_RNASEQ    = "sc_rnaseq"
    ATAC_SEQ     = "atac_seq"
    RASL_SEQ     = "rasl_seq"
    QPCR         = "qpcr"
    CHIP_SEQ     = "chip_seq"
    METABOLOMICS = "metabolomics"
    LIPIDOMICS   = "lipidomics"
    PROTEOMICS   = "proteomics"

if TYPE_CHECKING:
    from src.models.pounce.project import Project, Person
    from src.models.pounce.biosample import Biosample
    from src.models.pounce.biospecimen import Biospecimen


@dataclass
@facets(category_fields=['experiment_type', 'platform_type'])
class Experiment(Node):
    name: Optional[str] = None
    description: Optional[str] = None
    design: Optional[str] = None
    workbook: Optional[WorkbookArtifact] = None
    access: Optional[AccessLevel] = None
    experiment_type: Optional[str] = None
    date: Optional[date] = None
    platform_type: Optional[PlatformType] = None
    platform_name: Optional[str] = None
    platform_provider: Optional[str] = None
    platform_output_type: Optional[str] = None
    public_repo_id: Optional[str] = None
    repo_url: Optional[str] = None
    raw_file_archive_dir: List[str] = field(default_factory=list)
    extraction_protocol: Optional[str] = None
    acquisition_method: Optional[str] = None
    metabolite_identification_description: Optional[str] = None
    experiment_data_file: Optional[str] = None
    attached_files: List[str] = field(default_factory=list)
    calculated_properties: Optional[Dict[str, float]] = None
    geo_submission_date: Optional[str] = None
    curation_notes: Optional[str] = None


@dataclass
class ProjectExperimentEdge(Relationship):
    start_node: "Project" = None
    end_node: "Experiment" = None


@dataclass
class ExperimentPersonEdge(Relationship):
    start_node: "Experiment" = None
    end_node: "Person" = None
    roles: List[str] = field(default_factory=list)


@dataclass
class RunBiosample(Node):
    biological_replicate_number: Optional[int] = None
    technical_replicate_number: Optional[int] = None
    run_order: Optional[int] = None
    batch: Optional[str] = None
    source_file_name: Optional[str] = None
    replicate_label: Optional[str] = None
    geo_sample_accession: Optional[str] = None
    geo_sample_url: Optional[str] = None
    geo_sra_url: Optional[str] = None
    sample_description: Optional[str] = None
    source_name: Optional[str] = None
    library_text: Optional[str] = None


@dataclass
class BiosampleRunBiosampleEdge(Relationship):
    start_node: "Biosample" = None
    end_node: "RunBiosample" = None


@dataclass
class ExperimentBiospecimenEdge(Relationship):
    start_node: "Experiment" = None
    end_node: "Biospecimen" = None
