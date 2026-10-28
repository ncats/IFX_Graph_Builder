"""Annotated dataclasses for parsed POUNCE Excel data.

Each class maps directly to one (or part of one) Excel sheet.  Every field
carries ``sheet_field`` metadata so that parsers and validators can be
auto-generated from the class definitions alone.

These are *parsed input objects* — NOT Node subclasses.  The
``PounceInputAdapter`` converts them into Nodes/Relationships for the ETL.
"""

from dataclasses import dataclass, field
from datetime import date
from typing import List, Optional

from src.input_adapters.pounce_sheets.sheet_field import sheet_field


# ---------------------------------------------------------------------------
# Project workbook — ProjectMeta (meta sheet, single-row key-value)
# ---------------------------------------------------------------------------

@dataclass
class ParsedProject:
    project_name: Optional[str] = sheet_field(
        key="project_name", sheet="ProjectMeta")
    description: Optional[str] = sheet_field(
        key="description", sheet="ProjectMeta")
    date: Optional[date] = sheet_field(
        key="date", sheet="ProjectMeta", parse="date")
    owner_names: Optional[List[str]] = sheet_field(
        key="owner_name", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    owner_emails: Optional[List[str]] = sheet_field(
        key="owner_email", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    collaborator_names: Optional[List[str]] = sheet_field(
        key="collaborator_name", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    collaborator_emails: Optional[List[str]] = sheet_field(
        key="collaborator_email", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    lab_groups: Optional[List[str]] = sheet_field(
        key="lab_groups", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    organizations: Optional[List[str]] = sheet_field(
        key="Organizations", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    program_tags: Optional[List[str]] = sheet_field(
        key="program_tag", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    project_biosample_tracking: Optional[str] = sheet_field(
        key="biosample_tracking", sheet="ProjectMeta")
    keywords: Optional[List[str]] = sheet_field(
        key="keywords", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    privacy_type: Optional[str] = sheet_field(
        key="privacy_type", sheet="ProjectMeta")
    project_type: Optional[List[str]] = sheet_field(
        key="project_type", sheet="ProjectMeta", parse="string_list",
        default_factory=list)
    rd_tag: Optional[str] = sheet_field(
        key="RD_tag", sheet="ProjectMeta", parse="bool")
    biosample_preparation: Optional[str] = sheet_field(
        key="biosample_preparation", sheet="ProjectMeta")


@dataclass
class ParsedPerson:
    name: Optional[str] = None
    email: Optional[str] = None
    roles: List[str] = field(default_factory=list)


# ---------------------------------------------------------------------------
# Project workbook — BioSampleMeta via BioSampleMap (mapped sheet, multi-row)
# ---------------------------------------------------------------------------

@dataclass
class ParsedBiosample:
    biosample_id: Optional[str] = sheet_field(
        key="biosample_id", sheet="BioSampleMeta")
    biosample_type: Optional[str] = sheet_field(
        key="biosample_type", sheet="BioSampleMeta")
    biosample_tracking: Optional[str] = sheet_field(
        key="biosample_tracking", sheet="BioSampleMeta")
    biosample_batch: Optional[str] = sheet_field(
        key="biosample_batch", sheet="BioSampleMeta")


@dataclass
class ParsedBiospecimen:
    biospecimen_id: Optional[str] = sheet_field(
        key="biospecimen_id", sheet="BioSampleMeta")
    biospecimen_type: Optional[str] = sheet_field(
        key="biospecimen_type", sheet="BioSampleMeta")
    biospecimen_description: Optional[str] = sheet_field(
        key="biospecimen_description", sheet="BioSampleMeta")
    organism_names: Optional[str] = sheet_field(
        key="organism_names", sheet="BioSampleMeta")
    organism_category: Optional[str] = sheet_field(
        key="organism_category", sheet="BioSampleMeta", parse="category")
    disease_names: Optional[List[str]] = sheet_field(
        key="disease_names", sheet="BioSampleMeta", parse="string_list",
        default_factory=list)
    disease_category: Optional[str] = sheet_field(
        key="disease_category", sheet="BioSampleMeta", parse="category")


@dataclass
class ParsedDemographics:
    age: Optional[str] = sheet_field(
        key="age", sheet="BioSampleMeta")
    race: Optional[str] = sheet_field(
        key="race", sheet="BioSampleMeta")
    ethnicity: Optional[str] = sheet_field(
        key="ethnicity", sheet="BioSampleMeta")
    sex: Optional[str] = sheet_field(
        key="sex", sheet="BioSampleMeta")
    demographic_categories: Optional[List[str]] = sheet_field(
        key="demographic{}_category", sheet="BioSampleMeta", indexed=True,
        parse="category", default_factory=list)
    phenotype_categories: Optional[List[str]] = sheet_field(
        key="phenotype{}_category", sheet="BioSampleMeta", indexed=True,
        parse="category", default_factory=list)


@dataclass
class ParsedExposure:
    names: Optional[List[str]] = sheet_field(
        key="exposure{}_names", sheet="BioSampleMeta", indexed=True,
        parse="string_list", default_factory=list)
    type: Optional[str] = sheet_field(
        key="exposure{}_type", sheet="BioSampleMeta", indexed=True)
    category: Optional[str] = sheet_field(
        key="exposure{}_category", sheet="BioSampleMeta", indexed=True,
        parse="category")
    concentration: Optional[str] = sheet_field(
        key="exposure{}_concentration", sheet="BioSampleMeta", indexed=True)
    concentration_unit: Optional[str] = sheet_field(
        key="exposure{}_unit", sheet="BioSampleMeta", indexed=True)
    duration: Optional[str] = sheet_field(
        key="exposure{}_time", sheet="BioSampleMeta", indexed=True)
    duration_unit: Optional[str] = sheet_field(
        key="exposure{}_time_unit", sheet="BioSampleMeta", indexed=True)
    start_time: Optional[str] = sheet_field(
        key="exposure{}_start", sheet="BioSampleMeta", indexed=True)
    end_time: Optional[str] = sheet_field(
        key="exposure{}_end", sheet="BioSampleMeta", indexed=True)
    condition_category: Optional[str] = sheet_field(
        key="condition{}_category", sheet="BioSampleMeta", indexed=True,
        parse="category")
    growth_media: Optional[str] = sheet_field(
        key="growth_media", sheet="BioSampleMeta")
    exposure_index: Optional[int] = None  # which exposure slot (1, 2, ...)


# ---------------------------------------------------------------------------
# Experiment workbook — ExperimentMeta (meta sheet)
# ---------------------------------------------------------------------------

@dataclass
class ParsedExperiment:
    experiment_name: Optional[str] = sheet_field(
        key="experiment_name", sheet="ExperimentMeta")
    associated_project_name: Optional[str] = sheet_field(
        key="associated_project_name", sheet="ExperimentMeta")
    experiment_description: Optional[str] = sheet_field(
        key="experiment_description", sheet="ExperimentMeta")
    experiment_design: Optional[str] = sheet_field(
        key="experiment_design", sheet="ExperimentMeta")
    experiment_type: Optional[str] = sheet_field(
        key="experiment_type", sheet="ExperimentMeta")
    date: Optional[date] = sheet_field(
        key="date", sheet="ExperimentMeta", parse="date")
    lead_data_generator: Optional[str] = sheet_field(
        key="lead_data_generator", sheet="ExperimentMeta")
    lead_data_generator_email: Optional[str] = sheet_field(
        key="lead_data_generator_email", sheet="ExperimentMeta")
    lead_informatician: Optional[str] = sheet_field(
        key="lead_informatician", sheet="ExperimentMeta")
    lead_informatician_email: Optional[str] = sheet_field(
        key="lead_informatician_email", sheet="ExperimentMeta")
    platform_type: Optional[str] = sheet_field(
        key="platform_type", sheet="ExperimentMeta")
    platform_name: Optional[str] = sheet_field(
        key="platform_name", sheet="ExperimentMeta")
    platform_provider: Optional[str] = sheet_field(
        key="platform_provider", sheet="ExperimentMeta")
    platform_output_type: Optional[str] = sheet_field(
        key="platform_output_type", sheet="ExperimentMeta")
    public_repo_id: Optional[str] = sheet_field(
        key="public_repo_id", sheet="ExperimentMeta")
    repo_url: Optional[str] = sheet_field(
        key="repo_url", sheet="ExperimentMeta")
    raw_file_archive_dir: Optional[List[str]] = sheet_field(
        key="raw_file_archive_dir", sheet="ExperimentMeta", parse="string_list",
        default_factory=list)
    extraction_protocol: Optional[str] = sheet_field(
        key="extraction_protocol", sheet="ExperimentMeta")
    acquisition_method: Optional[str] = sheet_field(
        key="acquisition_method", sheet="ExperimentMeta")
    metabolite_identification_description: Optional[str] = sheet_field(
        key="metabolite_identification_description", sheet="ExperimentMeta")
    protein_identification_description: Optional[str] = sheet_field(
        key="protein_identification_description", sheet="ExperimentMeta")
    experiment_data_file: Optional[str] = sheet_field(
        key="experiment_data_file", sheet="ExperimentMeta")
    attached_files: Optional[List[str]] = sheet_field(
        key="attached_file_{}", sheet="ExperimentMeta", indexed=True,
        default_factory=list)
    attached_file_descriptions: Optional[List[str]] = sheet_field(
        key="attached_file_{}_description", sheet="ExperimentMeta", indexed=True,
        default_factory=list)


# ---------------------------------------------------------------------------
# Experiment workbook — GeneMeta via GeneMap (mapped sheet, multi-row)
# ---------------------------------------------------------------------------

@dataclass
class ParsedGene:
    gene_id: Optional[str] = sheet_field(
        key="gene_id", sheet="GeneMeta")
    hgnc_gene_symbol: Optional[str] = sheet_field(
        key="hgnc_gene_symbol", sheet="GeneMeta")
    gene_biotype: Optional[str] = sheet_field(
        key="gene_biotype", sheet="GeneMeta")
    alternate_id: Optional[str] = sheet_field(
        key="alternate_id", sheet="GeneMeta")
    chromosome_location: Optional[str] = sheet_field(
        key="chromosome_location", sheet="GeneMeta")
    strand: Optional[str] = sheet_field(
        key="strand", sheet="GeneMeta")
    pathway_ids: Optional[str] = sheet_field(
        key="pathway_ids", sheet="GeneMeta")
    categories: Optional[List[str]] = sheet_field(
        key="category {}", sheet="GeneMeta", indexed=True, default_factory=list)


@dataclass
class ParsedProtein:
    protein_id: Optional[str] = sheet_field(
        key="protein_id", sheet="ProteinMeta")
    protein_name: Optional[str] = sheet_field(
        key="protein_name", sheet="ProteinMeta")
    gene_name: Optional[str] = sheet_field(
        key="gene_name", sheet="ProteinMeta")
    protein_class: Optional[str] = sheet_field(
        key="protein_class", sheet="ProteinMeta")
    alternate_id: Optional[str] = sheet_field(
        key="alternate_id", sheet="ProteinMeta")
    pathway_ids: Optional[str] = sheet_field(
        key="pathway_ids", sheet="ProteinMeta")
    categories: Optional[List[str]] = sheet_field(
        key="category {}", sheet="ProteinMeta", indexed=True, default_factory=list)


# ---------------------------------------------------------------------------
# Experiment workbook — MetabMeta via MetabMap (mapped sheet, multi-row)
# ---------------------------------------------------------------------------

@dataclass
class ParsedMetab:
    metab_id: Optional[str] = sheet_field(
        key="metab_id", sheet="MetabMeta")
    # In v3, metab_id is the measured feature ID; standardized_metab_id is
    # separately supplied and is the only v3 ID eligible for canonical lookup.
    standardized_metab_id: Optional[str] = sheet_field(
        key="standardized_metab_id", sheet="MetabMeta")
    metab_id_is_feature: bool = False
    source_file: Optional[str] = None
    metab_name: Optional[str] = sheet_field(
        key="metab_name", sheet="MetabMeta")
    metab_chemclass: Optional[str] = sheet_field(
        key="metab_chemclass", sheet="MetabMeta")
    identification_level: Optional[int] = sheet_field(
        key="identification_level", sheet="MetabMeta", parse="int")
    alternate_metab_id: Optional[str] = sheet_field(
        key="alternate_metab_id", sheet="MetabMeta")
    alternate_metab_symbol: Optional[str] = sheet_field(
        key="alternate_metab_symbol", sheet="MetabMeta")
    pathway_ids: Optional[str] = sheet_field(
        key="pathway_ids", sheet="MetabMeta")
    categories: Optional[List[str]] = sheet_field(
        key="category {}", sheet="MetabMeta", indexed=True, default_factory=list)


# ---------------------------------------------------------------------------
# Experiment workbook — PeakDataMeta (meta sheet, single-row key-value)
# ---------------------------------------------------------------------------

@dataclass
class ParsedPeakDataMeta:
    pre_processing_description: Optional[str] = sheet_field(
        key="pre_processing_description", sheet="PeakDataMeta")
    peri_processing_description: Optional[str] = sheet_field(
        key="peri_processing_description", sheet="PeakDataMeta")
    peakdata_tag: Optional[str] = sheet_field(
        key="peakdata_tag", sheet="PeakDataMeta")


# ---------------------------------------------------------------------------
# Experiment workbook — RawDataMeta (meta sheet, single-row key-value)
# ---------------------------------------------------------------------------

@dataclass
class ParsedRawDataMeta:
    pre_processing_description: Optional[str] = sheet_field(
        key="pre_processing_description", sheet="RawDataMeta")
    peakdata_tag: Optional[str] = sheet_field(
        key="peakdata_tag", sheet="RawDataMeta")
    reference_genome: Optional[str] = sheet_field(
        key="reference_genome", sheet="RawDataMeta")
    gene_or_transcript: Optional[str] = sheet_field(
        key="Gene or Transcript", sheet="RawDataMeta")


@dataclass
class ParsedProteinDataMeta:
    pre_processing_description: Optional[str] = sheet_field(
        key="pre_processing_description", sheet="ProteinDataMeta")
    peri_processing_description: Optional[str] = sheet_field(
        key="peri_processing_description", sheet="ProteinDataMeta")
    proteindata_tag: Optional[str] = sheet_field(
        key="proteindata_tag", sheet="ProteinDataMeta")


# ---------------------------------------------------------------------------
# Experiment workbook — RunBioSampleMeta via RunBioSampleMap (mapped sheet, multi-row)
# ---------------------------------------------------------------------------

@dataclass
class ParsedRunBiosample:
    run_biosample_id: Optional[str] = sheet_field(
        key="run_biosample_id", sheet="RunBioSampleMeta")
    biosample_id: Optional[str] = sheet_field(
        key="biosample_id", sheet="RunBioSampleMeta")
    biological_replicate_number: Optional[str] = sheet_field(
        key="biological_replicate_number", sheet="RunBioSampleMeta", parse="int")
    technical_replicate_number: Optional[str] = sheet_field(
        key="technical_replicate_number", sheet="RunBioSampleMeta", parse="int")
    runbiosample_run_order: Optional[str] = sheet_field(
        key="runbiosample_run_order",
        aliases=["biosample_run_order"],
        sheet="RunBioSampleMeta",
        parse="int")
    batch: Optional[str] = sheet_field(
        key="batch", sheet="RunBioSampleMeta")
    runbiosample_file: Optional[str] = sheet_field(
        key="runbiosample_file", sheet="RunBioSampleMeta")


# ---------------------------------------------------------------------------
# StatsResults workbook — StatsResultsMeta (meta sheet)
# ---------------------------------------------------------------------------

@dataclass
class ParsedStatsResultsMeta:
    statsresults_id: Optional[str] = sheet_field(
        key="statsresults_id", sheet="StatsResultsMeta")
    statsresults_name: Optional[str] = sheet_field(
        key="statsresults_name", sheet="StatsResultsMeta")
    stats_description: Optional[str] = sheet_field(
        key="stats_description", sheet="StatsResultsMeta")
    experiment_name: Optional[str] = sheet_field(
        key="experiment_name", sheet="StatsResultsMeta")
    lead_informatician: Optional[str] = sheet_field(
        key="lead_informatician", sheet="StatsResultsMeta")
    lead_informatician_email: Optional[str] = sheet_field(
        key="lead_informatician_email", sheet="StatsResultsMeta")
    pre_processing_description: Optional[str] = sheet_field(
        key="pre_processing_description", sheet="StatsResultsMeta")
    peri_processing_description: Optional[str] = sheet_field(
        key="peri_processing_description", sheet="StatsResultsMeta")
    effect_size: Optional[str] = sheet_field(
        key="ES", sheet="StatsResultsMeta")
    effect_size_pval: Optional[str] = sheet_field(
        key="ESPval", sheet="StatsResultsMeta")
    effect_size_adj_pval: Optional[str] = sheet_field(
        key="ESadjPval", sheet="StatsResultsMeta")
    effect_size_2: Optional[str] = sheet_field(
        key="ES2", sheet="StatsResultsMeta")
    effect_size_prefix: Optional[str] = sheet_field(
        key="ES_1_prefix", sheet="StatsResultsMeta")
    effect_size_pval_prefix: Optional[str] = sheet_field(
        key="ES_1_Pval_prefix", sheet="StatsResultsMeta")
    effect_size_adj_pval_prefix: Optional[str] = sheet_field(
        key="ES_1_adjPval_prefix", sheet="StatsResultsMeta")
    effect_size_2_prefix: Optional[str] = sheet_field(
        key="ES_2_prefix", sheet="StatsResultsMeta")
    effect_size_2_pval: Optional[str] = sheet_field(
        key="ES_2_Pval", sheet="StatsResultsMeta")
    effect_size_2_pval_prefix: Optional[str] = sheet_field(
        key="ES_2_Pval_prefix", sheet="StatsResultsMeta")
    effect_size_2_adj_pval: Optional[str] = sheet_field(
        key="ES_2_adjPval", sheet="StatsResultsMeta")
    effect_size_2_adj_pval_prefix: Optional[str] = sheet_field(
        key="ES_2_adjPval_prefix", sheet="StatsResultsMeta")
    data_analysis_code_link: Optional[str] = sheet_field(
        key="DataAnalysisCode_Link", sheet="StatsResultsMeta")
