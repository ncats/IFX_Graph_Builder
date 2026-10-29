from dataclasses import dataclass, field, asdict
from typing import List, Dict, Union, TYPE_CHECKING
from src.interfaces.simple_enum import SimpleEnum, LabeledIntEnum
from src.models.node import Node, Relationship
from src.models.pounce.workbook_artifact import WorkbookArtifact
from datetime import date

if TYPE_CHECKING:
    from src.models.pounce.biosample import Biosample

class AccessLevel(SimpleEnum):
    private = "private"
    ncats = "ncats"
    public = "public"

class LabGroup(LabeledIntEnum):
    CGB = (1, "CGB")
    ETB = (2, "ETB")
    TDB = (3, "TDB")
    ACC = (4, "ACC")
    ADME = (5, "ADME")
    RSC = (6, "RSC")
    IFX = (7, "IFX")
    SCTL = (8, "SCTL")
    TMB = (9, "TMB")

class ProjectType(LabeledIntEnum):
    DRUG_SCREENING                   = (1, "Drug Screening")
    MOLECULAR_EPIDEMIOLOGY           = (2, "Molecular Epidemiology")
    TOXICOLOGY                       = (3, "Toxicology")
    DISEASE_DEVELOPMENT              = (4, "Disease Development/Characterization")
    DRUG_MECHANISM                   = (5, "Drug Mechanism")
    KO_KD_SCREENING                  = (6, "KO/KD screening")
    DRUG_OR_OTHER_EXPOSURE_CHARACTERIZATION = (7, "Drug or other Exposure Characterization")

@dataclass
class Person(Node):
    name: str = None
    email: str = None

    def to_dict(self) -> Dict[str, str]:
        return asdict(self)

    @classmethod
    def from_dict(cls, data: dict):
        return cls(**data)

    @classmethod
    def make(cls, name: str, email: str = None) -> 'Person':
        person_id = email if email else name.replace(' ', '_')
        return cls(id=person_id.lower(), name=name, email=email)


@dataclass
class Project(Node):
    name: str = None
    description: str = None
    date: date = None
    workbook: WorkbookArtifact = None
    lab_groups: List[LabGroup] = field(default_factory=list)
    keywords: List[str] = field(default_factory=list)
    access: AccessLevel = None
    project_type: List[ProjectType] = field(default_factory=list)
    rare_disease_focus: bool = None
    sample_preparation: str = None

    @staticmethod
    def id_from_display_id(display_id: str) -> int:
        try:
            return int(display_id.split('_')[-1])
        except (ValueError, IndexError):
            raise ValueError(f"Invalid project display ID format: {display_id}")

    @staticmethod
    def display_id_from_id(id: Union[str, int]) -> str:
        if isinstance(id, str):
            id = int(id)
        return f"pounce_proj_{id:05d}"

@dataclass
class ProjectPersonEdge(Relationship):
    start_node: Project
    end_node: Person
    roles: List[str] = field(default_factory=list)


@dataclass
class ProjectBiosampleEdge(Relationship):
    start_node: "Project" = None
    end_node: "Biosample" = None
