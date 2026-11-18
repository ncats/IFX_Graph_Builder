"""Per-run single-cell data file references, distinct from numeric datasets."""

from dataclasses import dataclass
from typing import Optional, TYPE_CHECKING

from src.models.node import Node, Relationship

if TYPE_CHECKING:
    from src.models.pounce.experiment import Experiment, RunBiosample


@dataclass
class RawDataFiles(Node):
    run_biosample_id: str = None
    matrix_file: Optional[str] = None
    features_file: Optional[str] = None
    barcodes_file: Optional[str] = None
    seurat_file: Optional[str] = None
    seurat_description: Optional[str] = None
    geo_sample: Optional[str] = None
    feature_count: Optional[int] = None
    barcode_count: Optional[int] = None
    nonzero_count: Optional[int] = None
    supplementary_urls: Optional[str] = None


@dataclass
class ExperimentRawDataFilesEdge(Relationship):
    start_node: "Experiment" = None
    end_node: RawDataFiles = None


@dataclass
class RunBiosampleRawDataFilesEdge(Relationship):
    start_node: "RunBiosample" = None
    end_node: RawDataFiles = None
