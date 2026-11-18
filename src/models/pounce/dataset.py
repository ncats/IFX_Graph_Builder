from dataclasses import dataclass, field
from typing import Optional, List, TYPE_CHECKING

from src.core.decorators import facets
from src.models.node import Node, Relationship
from src.models.pounce.project import AccessLevel

if TYPE_CHECKING:
    from src.models.pounce.experiment import Experiment, RunBiosample
    from src.models.gene import Gene
    from src.models.metabolite import Metabolite
    from src.models.protein import Protein
    import pandas as pd


@dataclass
@facets(category_fields=['data_type'])
class Dataset(Node):
    data_type: Optional[str] = None
    access: Optional[AccessLevel] = None
    pre_processing_description: Optional[str] = None
    peri_processing_description: Optional[str] = None
    file_reference: Optional[str] = None
    row_count: Optional[int] = None
    column_count: Optional[int] = None
    gene_id_column: Optional[str] = None
    sample_columns: List[str] = field(default_factory=list)
    reported_data_processing: Optional[str] = None
    data_representation: Optional[str] = None
    source_feature_type_note: Optional[str] = None
    _data_frame: Optional["pd.DataFrame"] = field(default=None, repr=False)


@dataclass
class ExperimentDatasetEdge(Relationship):
    start_node: "Experiment" = None
    end_node: "Dataset" = None


@dataclass
class DatasetRunBiosampleEdge(Relationship):
    start_node: "Dataset" = None
    end_node: "RunBiosample" = None


@dataclass
class DatasetGeneEdge(Relationship):
    start_node: "Dataset" = None
    end_node: "Gene" = None


@dataclass
class DatasetMetaboliteEdge(Relationship):
    start_node: "Dataset" = None
    end_node: "Metabolite" = None


@dataclass
class DatasetProteinEdge(Relationship):
    start_node: "Dataset" = None
    end_node: "Protein" = None
