"""Container for all parsed metadata from a POUNCE submission.

Data matrices (RawData, PeakData, StatsReadyData) are intentionally excluded —
validators only need the metadata sheets, which keeps validation fast regardless
of data volume.
"""

from dataclasses import dataclass, field
from typing import Dict, List, Optional

from src.input_adapters.pounce_sheets.parsed_classes import (
    ParsedProject,
    ParsedPerson,
    ParsedBiosample,
    ParsedBiospecimen,
    ParsedDemographics,
    ParsedExposure,
    ParsedExperiment,
    ParsedGene,
    ParsedMetab,
    ParsedProtein,
    ParsedProteinDataMeta,
    ParsedPeakDataMeta,
    ParsedRawDataMeta,
    ParsedRunBiosample,
    ParsedStatsResultsMeta,
)


@dataclass
class ParsedPounceData:
    project: Optional[ParsedProject] = None
    people: List[ParsedPerson] = field(default_factory=list)
    param_maps: Dict[str, Dict[str, str]] = field(default_factory=dict)
    workbook_versions: Dict[str, int] = field(default_factory=dict)
    original_sheets: Dict[str, Dict[str, str]] = field(default_factory=dict)
    biosamples: List[ParsedBiosample] = field(default_factory=list)
    biospecimens: List[ParsedBiospecimen] = field(default_factory=list)
    demographics: List[ParsedDemographics] = field(default_factory=list)
    exposures: List[ParsedExposure] = field(default_factory=list)
    experiments: List[ParsedExperiment] = field(default_factory=list)
    genes: List[ParsedGene] = field(default_factory=list)
    proteins: List[ParsedProtein] = field(default_factory=list)
    metabolites: List[ParsedMetab] = field(default_factory=list)
    protein_data_meta: List[ParsedProteinDataMeta] = field(default_factory=list)
    peak_data_meta: List[ParsedPeakDataMeta] = field(default_factory=list)
    raw_data_meta: List[ParsedRawDataMeta] = field(default_factory=list)
    run_biosamples: List[ParsedRunBiosample] = field(default_factory=list)
    stats_results: List[ParsedStatsResultsMeta] = field(default_factory=list)
