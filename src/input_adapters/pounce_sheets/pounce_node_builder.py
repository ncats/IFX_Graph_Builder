import hashlib
import os
import pandas as pd
import re
from dataclasses import dataclass, field
from datetime import datetime
from typing import Any, Dict, Generator, List, Optional, Tuple, Type, Union

from src.constants import Prefix
from src.input_adapters.excel_sheet_adapter import ExcelsheetParser
from src.input_adapters.pounce_sheets.constants import ExperimentWorkbook, ProjectWorkbook, StatsResultsWorkbook


_NA_VALUES = {"NA", "N/A", ""}

from src.input_adapters.pounce_sheets.parsed_classes import ParsedPerson, ParsedProject
from src.input_adapters.pounce_sheets.parsed_pounce_data import ParsedPounceData
from src.input_adapters.pounce_sheets.pounce_parser import PounceParser
from src.input_adapters.pounce_sheets.workbook_format import load_pounce_workbook
from src.models.gene import MeasuredGene, MeasuredGeneEdge, Gene
from src.models.metabolite import MeasuredMetabolite, MeasuredMetaboliteEdge, Metabolite
from src.models.protein import MeasuredProtein, MeasuredProteinEdge, Protein
from src.models.node import Node, Relationship, EquivalentId
from src.models.pounce.biosample import Biosample, BiosampleBiospecimenEdge
from src.models.pounce.biospecimen import Biospecimen
from src.models.pounce.category_value import CategoryValue
from src.models.pounce.dataset import Dataset, ExperimentDatasetEdge, DatasetRunBiosampleEdge, DatasetGeneEdge, DatasetMetaboliteEdge, DatasetProteinEdge
from src.models.pounce.demographics import Demographics
from src.models.pounce.exposure import Exposure
from src.models.pounce.stats_result import StatsResult, ComparisonColumn, ExperimentStatsResultEdge, StatsResultGeneEdge, StatsResultMetaboliteEdge, StatsResultPersonEdge, StatsResultProteinEdge
from src.models.pounce.experiment import Experiment, PlatformType, ProjectExperimentEdge, ExperimentPersonEdge, RunBiosample, BiosampleRunBiosampleEdge, ExperimentBiospecimenEdge
from src.models.pounce.exposure import BiosampleExposureEdge
from src.models.pounce.project import Project, AccessLevel, ProjectType, Person, ProjectPersonEdge, \
    ProjectBiosampleEdge, LabGroup
from src.models.pounce.workbook_artifact import WorkbookArtifact


@dataclass
class PounceContext:
    project: Optional[Project] = None
    biosample_by_original_id: Dict[str, Biosample] = field(default_factory=dict)
    biospecimen_by_biosample_id: Dict[str, Biospecimen] = field(default_factory=dict)
    run_biosample_by_original_id: Dict[str, RunBiosample] = field(default_factory=dict)
    gene_by_raw_id: Dict[str, Any] = field(default_factory=dict)
    protein_by_raw_id: Dict[str, Any] = field(default_factory=dict)
    metabolite_by_raw_id: Dict[str, Any] = field(default_factory=dict)
    experiment_counter: int = 0


def _iter_pairs(experiment_files, stats_files):
    """Yield (exp_file, stats_file_or_None) pairs."""
    if stats_files and len(stats_files) != len(experiment_files):
        raise ValueError(
            "StatsResults workbook slots must match experiment count; use None for an experiment without stats"
        )
    for i, exp_file in enumerate(experiment_files):
        yield exp_file, stats_files[i] if i < len(stats_files) else None


class PounceNodeBuilder:
    @staticmethod
    def _slugify(value: Optional[str], max_len: int = 48) -> str:
        if not value:
            return "unnamed"
        slug = re.sub(r"[^a-z0-9]+", "-", value.strip().lower())
        slug = re.sub(r"-{2,}", "-", slug).strip("-")
        return slug[:max_len].strip("-") or "unnamed"

    @staticmethod
    def _hash_suffix(*parts: Optional[str], length: int = 8) -> str:
        normalized = "|".join("" if part is None else str(part).strip().lower() for part in parts)
        return hashlib.sha1(normalized.encode("utf-8")).hexdigest()[:length]

    @classmethod
    def _project_id(cls, parsed_project: ParsedProject) -> str:
        slug = cls._slugify(parsed_project.project_name)
        date_part = parsed_project.date.isoformat() if parsed_project.date else None
        owner_emails = sorted((parsed_project.owner_emails or []))
        collaborator_emails = sorted((parsed_project.collaborator_emails or []))
        suffix = cls._hash_suffix(
            parsed_project.project_name,
            date_part,
            "|".join(owner_emails),
            "|".join(collaborator_emails),
        )
        return f"pounce_{slug}-{suffix}"

    @classmethod
    def _experiment_id(cls, project_id: str, parsed_experiment) -> str:
        slug = cls._slugify(parsed_experiment.experiment_name)
        date_part = parsed_experiment.date.isoformat() if hasattr(parsed_experiment.date, "isoformat") else parsed_experiment.date
        suffix = cls._hash_suffix(
            project_id,
            parsed_experiment.experiment_name,
            date_part,
            parsed_experiment.platform_type,
            parsed_experiment.platform_name,
            parsed_experiment.platform_output_type,
        )
        return f"{project_id}_{slug}-{suffix}"

    @staticmethod
    def _to_null_if_na(value: Optional[str]) -> Optional[str]:
        """Convert a value to None if it is missing, empty, or 'N/A' (case-insensitive)."""
        if value is None:
            return None
        stripped = value.strip()
        if not stripped or stripped.upper() in _NA_VALUES:
            return None
        return value

    def __init__(self):
        self._ctx = PounceContext()

    def build(self, project_file, experiment_files, stats_results_files, parser) -> Generator:
        print(f"reading: {project_file}")
        project_data, _ = parser.parse_project(load_pounce_workbook(project_file, "project").parser)
        proj_obj = self._project_node(project_data.project, workbook_local_path=project_file)
        person_nodes, person_edges = self._person_nodes_and_edges(proj_obj, project_data.people)
        yield [proj_obj, *person_nodes, *person_edges]

        yield from self._biosample_nodes(proj_obj, project_data)

        for exp_file, stats_file in _iter_pairs(experiment_files, stats_results_files):
            exp_parser = load_pounce_workbook(exp_file, "experiment").parser
            exp_data, _ = parser.parse_experiment(exp_parser)
            stats_parser = load_pounce_workbook(stats_file, "stats").parser if stats_file else None
            stats_data, _ = parser.parse_stats_results(stats_parser) if stats_parser else (None, [])
            yield from self._experiment_nodes(
                proj_obj,
                exp_data,
                exp_parser,
                stats_data,
                stats_parser,
                experiment_file=exp_file,
                stats_file=stats_file,
            )

    # --- Project workbook ---

    @staticmethod
    def _project_node(parsed_project: ParsedProject, workbook_local_path: Optional[str] = None) -> Project:
        return Project(
            id=PounceNodeBuilder._project_id(parsed_project),
            name=parsed_project.project_name,
            description=parsed_project.description,
            date=parsed_project.date,
            workbook=WorkbookArtifact(
                original_filename=os.path.basename(workbook_local_path) if workbook_local_path else None,
                _local_path=workbook_local_path,
            ) if workbook_local_path else None,
            lab_groups=[group for raw in (parsed_project.lab_groups or []) if (group := LabGroup.parse(raw)) is not None],
            access=AccessLevel.parse(parsed_project.privacy_type) if parsed_project.privacy_type else None,
            keywords=parsed_project.keywords or [],
            project_type=[pt for raw in (parsed_project.project_type or []) if (pt := ProjectType.parse(raw)) is not None],
            rare_disease_focus=parsed_project.rd_tag if parsed_project.rd_tag is not None else False,
            sample_preparation=parsed_project.biosample_preparation
        )

    @staticmethod
    def _person_nodes_and_edges(
        proj_obj: Project, people: List[ParsedPerson]
    ) -> Tuple[List[Person], List[ProjectPersonEdge]]:
        nodes = []
        edges = []
        for parsed_person in people:
            person_node = Person.make(parsed_person.name, parsed_person.email)
            nodes.append(person_node)
            edges.append(ProjectPersonEdge(start_node=proj_obj, end_node=person_node, roles=parsed_person.roles))
        return nodes, edges

    def _biosample_nodes(
        self, proj_obj: Project, project_data: ParsedPounceData
    ) -> Generator[List[Union[Node, Relationship]], None, None]:
        biosample_param_map = project_data.param_maps.get(ProjectWorkbook.BiosampleMapSheet.name) or {}
        num_exposure_slots = len(PounceParser._detect_exposure_indices(biosample_param_map))

        biospecimens = {}
        biosamples = []
        exposures = {}
        sample_exposure_edges: List[BiosampleExposureEdge] = []
        project_biosample_edges: List[ProjectBiosampleEdge] = []
        biosample_biospecimen_edges: List[BiosampleBiospecimenEdge] = []

        for i, parsed_bs in enumerate(project_data.biosamples):
            parsed_spec = project_data.biospecimens[i]

            biospecimen_id = parsed_spec.biospecimen_id
            full_biospecimen_id = f"{proj_obj.id}-{biospecimen_id}"
            if full_biospecimen_id not in biospecimens:
                organism_category = self._make_category_value(
                    parsed_spec.organism_category,
                    biosample_param_map.get("organism_category")) if parsed_spec.organism_category else None
                disease_category = self._make_category_value(
                    parsed_spec.disease_category,
                    biosample_param_map.get("disease_category")) if parsed_spec.disease_category else None

                biospecimens[full_biospecimen_id] = Biospecimen(
                    id=full_biospecimen_id,
                    original_id=biospecimen_id,
                    type=parsed_spec.biospecimen_type,
                    description=parsed_spec.biospecimen_description,
                    organism=parsed_spec.organism_names,
                    organism_category=organism_category,
                    disease_category=disease_category,
                    diseases=parsed_spec.disease_names or []
                )

            full_biosample_id = f"{proj_obj.id}-{parsed_bs.biosample_id}-{biospecimen_id}"

            parsed_demo = project_data.demographics[i] if i < len(project_data.demographics) else None
            demographics = self._convert_demographics(parsed_demo, full_biosample_id) if parsed_demo else None

            biosample_obj = Biosample(
                id=full_biosample_id,
                original_id=parsed_bs.biosample_id,
                type=parsed_bs.biosample_type,
                tracking=self._to_null_if_na(parsed_bs.biosample_tracking),
                demographics=demographics
            )
            biosamples.append(biosample_obj)
            self._ctx.biosample_by_original_id[str(parsed_bs.biosample_id)] = biosample_obj
            self._ctx.biospecimen_by_biosample_id[biosample_obj.id] = biospecimens[full_biospecimen_id]

            project_biosample_edges.append(ProjectBiosampleEdge(start_node=proj_obj, end_node=biosample_obj))
            biosample_biospecimen_edges.append(
                BiosampleBiospecimenEdge(start_node=biosample_obj, end_node=biospecimens[full_biospecimen_id])
            )

            if num_exposure_slots > 0:
                row_exposures = project_data.exposures[i * num_exposure_slots:(i + 1) * num_exposure_slots]
                for parsed_exp in row_exposures:
                    exposure_obj = self._convert_exposure(parsed_exp)
                    sample_exposure_edges.append(BiosampleExposureEdge(start_node=biosample_obj, end_node=exposure_obj))
                    if exposure_obj.id not in exposures:
                        exposures[exposure_obj.id] = exposure_obj

        yield list(biospecimens.values())
        yield biosamples
        yield list(exposures.values())
        yield sample_exposure_edges
        yield project_biosample_edges
        yield biosample_biospecimen_edges

    # --- Experiment workbook ---

    def _experiment_nodes(
        self,
        proj_obj: Project,
        exp_data: ParsedPounceData,
        exp_parser: ExcelsheetParser,
        stats_data: Optional[ParsedPounceData] = None,
        stats_parser: Optional[ExcelsheetParser] = None,
        experiment_file: Optional[str] = None,
        stats_file: Optional[str] = None,
    ) -> Generator[List[Union[Node, Relationship]], None, None]:
        print(f"reading: {exp_parser.file_path}")

        parsed = exp_data.experiments[0]

        exp_id = self._experiment_id(proj_obj.id, parsed)

        exp_date = None
        if parsed.date:
            if isinstance(parsed.date, datetime):
                exp_date = parsed.date.date()
            else:
                try:
                    exp_date = datetime.strptime(str(parsed.date), "%Y%m%d").date()
                except ValueError:
                    exp_date = parsed.date

        analyte_count, biosample_count = self._compute_experiment_counts(exp_data, exp_parser, stats_parser)

        experiment_obj = Experiment(
            id=exp_id,
            name=parsed.experiment_name,
            description=parsed.experiment_description,
            design=parsed.experiment_design,
            workbook=WorkbookArtifact(
                original_filename=os.path.basename(experiment_file) if experiment_file else None,
                _local_path=experiment_file,
            ) if experiment_file else None,
            access=proj_obj.access,
            experiment_type=parsed.experiment_type,
            date=exp_date,
            platform_type=PlatformType.parse(parsed.platform_type),
            platform_name=parsed.platform_name,
            platform_provider=parsed.platform_provider,
            platform_output_type=parsed.platform_output_type,
            public_repo_id=PounceNodeBuilder._to_null_if_na(parsed.public_repo_id),
            repo_url=PounceNodeBuilder._to_null_if_na(parsed.repo_url),
            raw_file_archive_dir=parsed.raw_file_archive_dir or [],
            extraction_protocol=parsed.extraction_protocol,
            acquisition_method=parsed.acquisition_method,
            metabolite_identification_description=parsed.metabolite_identification_description,
            experiment_data_file=parsed.experiment_data_file,
            attached_files=parsed.attached_files or [],
            calculated_properties={
                "analyte_count": float(analyte_count),
                "biosample_count": float(biosample_count),
            }
        )

        data_generator = Person.make(parsed.lead_data_generator, parsed.lead_data_generator_email) \
            if parsed.lead_data_generator else None
        informatician = Person.make(parsed.lead_informatician, parsed.lead_informatician_email) \
            if parsed.lead_informatician else None

        persons = [p for p in [data_generator, informatician] if p is not None]
        person_edges = []
        if data_generator:
            person_edges.append(ExperimentPersonEdge(start_node=experiment_obj, end_node=data_generator, roles=["DataGenerator"]))
        if informatician:
            person_edges.append(ExperimentPersonEdge(start_node=experiment_obj, end_node=informatician, roles=["Informatician"]))

        yield [experiment_obj, *persons, ProjectExperimentEdge(start_node=proj_obj, end_node=experiment_obj), *person_edges]

        # Run biosamples
        run_biosamples = []
        rb_edges: List[BiosampleRunBiosampleEdge] = []
        experiment_biospecimen_edges: List[ExperimentBiospecimenEdge] = []
        seen_biospecimen_ids = set()
        for parsed_rb in exp_data.run_biosamples:
            run_biosample_id = parsed_rb.run_biosample_id
            run_biosample_obj = RunBiosample(
                id=f"{proj_obj.id}-{run_biosample_id}",
                biological_replicate_number=self._parse_int(parsed_rb.biological_replicate_number),
                technical_replicate_number=self._parse_int(parsed_rb.technical_replicate_number),
                run_order=self._parse_int(parsed_rb.runbiosample_run_order),
                batch=parsed_rb.batch
            )
            run_biosamples.append(run_biosample_obj)
            self._ctx.run_biosample_by_original_id[str(run_biosample_id)] = run_biosample_obj

            biosample_id = str(parsed_rb.biosample_id)
            if biosample_id in self._ctx.biosample_by_original_id:
                biosample_obj = self._ctx.biosample_by_original_id[biosample_id]
                rb_edges.append(BiosampleRunBiosampleEdge(
                    start_node=biosample_obj,
                    end_node=run_biosample_obj
                ))
                biospecimen_obj = self._ctx.biospecimen_by_biosample_id.get(biosample_obj.id)
                if biospecimen_obj is None:
                    print(f"Warning: Biosample has no biospecimen mapping: {biosample_obj.id}")
                elif biospecimen_obj.id not in seen_biospecimen_ids:
                    seen_biospecimen_ids.add(biospecimen_obj.id)
                    experiment_biospecimen_edges.append(
                        ExperimentBiospecimenEdge(start_node=experiment_obj, end_node=biospecimen_obj)
                    )
            else:
                print(f"Warning: RunBiosample references unknown biosample_id: {biosample_id}")

        yield run_biosamples
        yield rb_edges
        if experiment_biospecimen_edges:
            yield experiment_biospecimen_edges

        sheet_names = exp_parser.sheet_dfs.keys()
        if ExperimentWorkbook.GeneMetaSheet.name in sheet_names:
            yield from self._parse_genes(experiment_obj, exp_parser)
        if ExperimentWorkbook.ProteinMetaSheet.name in sheet_names:
            yield from self._parse_proteins(experiment_obj, exp_parser)
        if ExperimentWorkbook.MetabMetaSheet.name in sheet_names:
            yield from self._parse_metabolites(experiment_obj, exp_parser)
        if ExperimentWorkbook.RawDataSheet.name in sheet_names:
            gene_map = exp_parser.get_parameter_map(ExperimentWorkbook.GeneMapSheet.name)
            analyte_id_col = gene_map.get(ExperimentWorkbook.GeneMapSheet.Key.gene_id)
            yield from self._parse_data_matrix(
                experiment_obj,
                meta_sheet=ExperimentWorkbook.RawDataMetaSheet.name,
                data_sheet=ExperimentWorkbook.RawDataSheet.name,
                analyte_id_col=analyte_id_col,
                data_type="raw data",
                parser=exp_parser
            )
        if ExperimentWorkbook.PeakDataSheet.name in sheet_names:
            metab_map = exp_parser.get_parameter_map(ExperimentWorkbook.MetabMapSheet.name)
            analyte_id_col = metab_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_id)
            yield from self._parse_data_matrix(
                experiment_obj,
                meta_sheet=ExperimentWorkbook.PeakDataMetaSheet.name,
                data_sheet=ExperimentWorkbook.PeakDataSheet.name,
                analyte_id_col=analyte_id_col,
                data_type="raw data",
                parser=exp_parser
            )
        if ExperimentWorkbook.ProteinDataSheet.name in sheet_names:
            protein_map = exp_parser.get_parameter_map(ExperimentWorkbook.ProteinMapSheet.name)
            analyte_id_col = protein_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_id)
            yield from self._parse_data_matrix(
                experiment_obj,
                meta_sheet=ExperimentWorkbook.ProteinDataMetaSheet.name,
                data_sheet=ExperimentWorkbook.ProteinDataSheet.name,
                analyte_id_col=analyte_id_col,
                data_type="raw data",
                parser=exp_parser
            )

        if stats_parser:
            yield from self._stats_results_nodes(experiment_obj, stats_parser, stats_data, stats_file=stats_file)

    def _compute_experiment_counts(
        self, exp_data: ParsedPounceData, exp_parser: ExcelsheetParser,
        stats_parser: Optional[ExcelsheetParser] = None
    ) -> Tuple[int, int]:
        """Return (analyte_count, biosample_count) from actual data matrix contents.

        Counts unique analytes and biosamples across all data matrices (RawData,
        PeakData, StatsReadyData), matching what the database reports after an
        inner-join query.
        """
        analyte_ids: set = set()
        run_biosample_ids: set = set()
        sheet_names = exp_parser.sheet_dfs.keys()

        if ExperimentWorkbook.RawDataSheet.name in sheet_names:
            gene_map = exp_parser.get_parameter_map(ExperimentWorkbook.GeneMapSheet.name)
            analyte_id_col = gene_map.get(ExperimentWorkbook.GeneMapSheet.Key.gene_id)
            df = self._prepare_data_frame(exp_parser, ExperimentWorkbook.RawDataSheet.name, analyte_id_col)
            if df is not None:
                analyte_ids.update(str(i) for i in df.index)
                run_biosample_ids.update(str(c) for c in df.columns)

        if ExperimentWorkbook.ProteinDataSheet.name in sheet_names:
            protein_map = exp_parser.get_parameter_map(ExperimentWorkbook.ProteinMapSheet.name)
            analyte_id_col = protein_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_id)
            df = self._prepare_data_frame(exp_parser, ExperimentWorkbook.ProteinDataSheet.name, analyte_id_col)
            if df is not None:
                analyte_ids.update(str(i) for i in df.index)
                run_biosample_ids.update(str(c) for c in df.columns)

        if ExperimentWorkbook.PeakDataSheet.name in sheet_names:
            metab_map = exp_parser.get_parameter_map(ExperimentWorkbook.MetabMapSheet.name)
            analyte_id_col = metab_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_id)
            df = self._prepare_data_frame(exp_parser, ExperimentWorkbook.PeakDataSheet.name, analyte_id_col)
            if df is not None:
                analyte_ids.update(str(i) for i in df.index)
                run_biosample_ids.update(str(c) for c in df.columns)

        if stats_parser:
            analyte_id_col, _ = self._parse_effect_size_map(stats_parser)
            if StatsResultsWorkbook.StatsReadyDataSheet.name in stats_parser.sheet_dfs:
                df = self._prepare_data_frame(stats_parser, StatsResultsWorkbook.StatsReadyDataSheet.name, analyte_id_col)
                if df is not None:
                    analyte_ids.update(str(i) for i in df.index)
                    run_biosample_ids.update(str(c) for c in df.columns)
            if StatsResultsWorkbook.EffectSizeSheet.name in stats_parser.sheet_dfs:
                df = self._prepare_data_frame(stats_parser, StatsResultsWorkbook.EffectSizeSheet.name, analyte_id_col)
                if df is not None:
                    analyte_ids.update(str(i) for i in df.index)

        rb_to_bs = {
            str(rb.run_biosample_id): str(rb.biosample_id)
            for rb in exp_data.run_biosamples
            if rb.run_biosample_id and rb.biosample_id
        }
        if run_biosample_ids:
            biosample_ids = {rb_to_bs[rb_id] for rb_id in run_biosample_ids if rb_id in rb_to_bs}
        else:
            biosample_ids = set(rb_to_bs.values())

        return len(analyte_ids), len(biosample_ids)

    def _get_metabolite_categories(self, row, column_map: dict) -> List[CategoryValue]:
        categories = []
        cat_num = 1
        while True:
            key = ExperimentWorkbook.MetabMapSheet.Key.category.format(cat_num)
            column_name = column_map.get(key)
            if not column_name:
                break
            value = row.get(column_name)
            if value is not None and str(value).strip() and str(value).strip().lower() not in ('nan', 'na'):
                categories.append(CategoryValue(
                    id=f"{column_name}-{value}",
                    name=column_name, value=str(value).strip())
                )
            cat_num += 1
        return categories

    # --- Analyte parsing ---

    def _parse_genes(self, experiment_obj: Experiment, parser: ExcelsheetParser
                     ) -> Generator[List[Union[Node, Relationship]], None, None]:
        """Parse gene data from GeneMeta sheet."""
        column_map = parser.get_parameter_map(ExperimentWorkbook.GeneMapSheet.name)
        gene_id_col = column_map.get(ExperimentWorkbook.GeneMapSheet.Key.gene_id)
        symbol_col = column_map.get(ExperimentWorkbook.GeneMapSheet.Key.hgnc_gene_symbol)
        biotype_col = column_map.get(ExperimentWorkbook.GeneMapSheet.Key.gene_biotype)

        genes = []
        gene_edges = []
        for _, row in parser.sheet_dfs[ExperimentWorkbook.GeneMetaSheet.name].iterrows():
            gene_id = row.get(gene_id_col) if gene_id_col else None
            if not gene_id:
                continue
            gene_ensembl_id = EquivalentId(id=gene_id, type=Prefix.ENSEMBL).id_str()
            gene_obj = MeasuredGene(
                id=gene_ensembl_id,
                symbol=row.get(symbol_col) if symbol_col else None,
                biotype=row.get(biotype_col) if biotype_col else None
            )
            genes.append(gene_obj)
            self._ctx.gene_by_raw_id[str(gene_id)] = gene_obj

            gene_edges.append(MeasuredGeneEdge(start_node=gene_obj, end_node=Gene(id=gene_ensembl_id)))

        yield genes
        yield gene_edges

    def _parse_metabolites(self, experiment_obj: Experiment, parser: ExcelsheetParser
                           ) -> Generator[List[Union[Node, Relationship]], None, None]:
        """Parse metabolite data from MetabMeta sheet."""
        column_map = parser.get_parameter_map(ExperimentWorkbook.MetabMapSheet.name)
        metab_id_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_id)
        metab_name_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_name)
        metab_type_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.metab_chemclass)
        id_level_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.identification_level)
        alternate_id_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.alternate_metab_id)
        standardized_id_col = column_map.get("standardized_metab_id")
        feature_id_format = getattr(parser, "pounce_version", 2) == 3
        alternate_symbol_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.alternate_metab_symbol)
        pathway_ids_col = column_map.get(ExperimentWorkbook.MetabMapSheet.Key.pathway_ids)

        metabolites = []
        metabolite_edges = []

        for _, row in parser.sheet_dfs[ExperimentWorkbook.MetabMetaSheet.name].iterrows():
            metab_id = row.get(metab_id_col) if metab_id_col else None
            if not metab_id:
                continue

            alternate_ids = self._split_pipe(row.get(alternate_id_col) if alternate_id_col else None)
            standardized_ids = self._split_pipe(row.get(standardized_id_col) if standardized_id_col else None)
            alternate_names = self._split_pipe(row.get(alternate_symbol_col) if alternate_symbol_col else None)
            pathway_ids = self._split_pipe(row.get(pathway_ids_col) if pathway_ids_col else None)

            metab_obj = MeasuredMetabolite(
                id=f"{experiment_obj.id}-{metab_id}",
                name=row.get(metab_name_col) if metab_name_col else None,
                type=row.get(metab_type_col) if metab_type_col else None,
                identification_level=self._parse_int(row.get(id_level_col)) if id_level_col else None,
                alternate_ids=alternate_ids,
                alternate_names=alternate_names,
                pathway_ids=pathway_ids,
                categories=self._get_metabolite_categories(row, column_map)
            )
            metabolites.append(metab_obj)
            self._ctx.metabolite_by_raw_id[str(metab_id)] = metab_obj

            canonical_ids = ([] if feature_id_format else [metab_id]) + standardized_ids + alternate_ids
            for raw_id in dict.fromkeys(canonical_ids):
                metabolite_edges.append(MeasuredMetaboliteEdge(start_node=metab_obj, end_node=Metabolite(id=raw_id)))

        yield metabolites
        yield metabolite_edges

    def _parse_proteins(self, experiment_obj: Experiment, parser: ExcelsheetParser
                        ) -> Generator[List[Union[Node, Relationship]], None, None]:
        """Parse protein data from ProteinMeta sheet."""
        column_map = parser.get_parameter_map(ExperimentWorkbook.ProteinMapSheet.name)
        protein_id_col = column_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_id)
        protein_name_col = column_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_name)
        protein_class_col = column_map.get(ExperimentWorkbook.ProteinMapSheet.Key.protein_class)

        proteins = []
        protein_edges = []
        for _, row in parser.sheet_dfs[ExperimentWorkbook.ProteinMetaSheet.name].iterrows():
            protein_id = row.get(protein_id_col) if protein_id_col else None
            if not protein_id:
                continue

            protein_uniprot_id = EquivalentId(id=protein_id, type=Prefix.UniProtKB).id_str()
            protein_obj = MeasuredProtein(
                id=protein_uniprot_id,
                name=row.get(protein_name_col) if protein_name_col else None,
                type=row.get(protein_class_col) if protein_class_col else None,
            )
            proteins.append(protein_obj)
            self._ctx.protein_by_raw_id[str(protein_id)] = protein_obj
            protein_edges.append(
                MeasuredProteinEdge(start_node=protein_obj, end_node=Protein(id=protein_uniprot_id))
            )

        yield proteins
        yield protein_edges

    # --- Data matrix and stats results parsing ---

    @staticmethod
    def _prepare_data_frame(parser: ExcelsheetParser, data_sheet: str, analyte_id_col: str = None):
        """Clean a data sheet DataFrame: drop empty columns/rows, set analyte ID as index.
        Returns cleaned_df or None if empty."""
        raw_df = parser.sheet_dfs[data_sheet]

        if raw_df.empty or len(raw_df.columns) == 0:
            print(f"Empty data sheet: {data_sheet}")
            return None

        raw_df = raw_df.dropna(axis=1, how='all')
        raw_df = raw_df.dropna(subset=[analyte_id_col])

        if raw_df.empty:
            print(f"No valid rows in data sheet: {data_sheet}")
            return None

        return raw_df.set_index(analyte_id_col)

    def _build_analyte_edges(self, parent_node, df, gene_edge_cls: Type[Relationship],
                             metab_edge_cls: Type[Relationship], protein_edge_cls: Type[Relationship]) -> List[Relationship]:
        """Build edges from a parent node (Dataset or StatsResult) to Gene/Metabolite/Protein nodes."""
        edges = []
        for analyte_id in df.index:
            raw_id = str(analyte_id)
            if raw_id in self._ctx.gene_by_raw_id:
                edges.append(gene_edge_cls(start_node=parent_node, end_node=self._ctx.gene_by_raw_id[raw_id]))
            elif raw_id in self._ctx.metabolite_by_raw_id:
                edges.append(metab_edge_cls(start_node=parent_node, end_node=self._ctx.metabolite_by_raw_id[raw_id]))
            elif raw_id in self._ctx.protein_by_raw_id:
                edges.append(protein_edge_cls(start_node=parent_node, end_node=self._ctx.protein_by_raw_id[raw_id]))
        return edges

    def _qualify_data_frame_ids(self, df):
        """Rename DataFrame index and columns to use qualified entity IDs."""
        index_rename = {}
        for raw_id in df.index:
            raw_id_str = str(raw_id)
            if raw_id_str in self._ctx.metabolite_by_raw_id:
                index_rename[raw_id] = self._ctx.metabolite_by_raw_id[raw_id_str].id
            elif raw_id_str in self._ctx.gene_by_raw_id:
                index_rename[raw_id] = self._ctx.gene_by_raw_id[raw_id_str].id
            elif raw_id_str in self._ctx.protein_by_raw_id:
                index_rename[raw_id] = self._ctx.protein_by_raw_id[raw_id_str].id

        col_rename = {
            col: self._ctx.run_biosample_by_original_id[str(col)].id
            for col in df.columns if str(col) in self._ctx.run_biosample_by_original_id
        }
        if index_rename:
            df = df.rename(index=index_rename)
        if col_rename:
            df = df.rename(columns=col_rename)
        return df

    def _parse_data_matrix(self, experiment_obj: Experiment, meta_sheet: str,
                           data_sheet: str, analyte_id_col: str = None,
                           data_type: str = "raw data",
                           parser: ExcelsheetParser = None
                           ) -> Generator[List[Union[Node, Relationship]], None, None]:
        """Parse a data matrix sheet (RawData, PeakData, or StatsReadyData) into a Dataset node."""
        pre_processing = parser.safe_get_string(meta_sheet, "pre_processing_description")
        peri_processing = parser.safe_get_string(meta_sheet, "peri_processing_description")
        dataset_id_suffix = data_type.replace(" ", "_")
        dataset_id = f"{experiment_obj.id}:{dataset_id_suffix}"

        raw_df = self._prepare_data_frame(parser, data_sheet, analyte_id_col)

        if raw_df is None:
            dataset = Dataset(
                id=dataset_id, data_type=data_type,
                access=experiment_obj.access,
                pre_processing_description=pre_processing, peri_processing_description=peri_processing,
                row_count=0, column_count=0, gene_id_column=analyte_id_col, sample_columns=[]
            )
            yield [dataset, ExperimentDatasetEdge(start_node=experiment_obj, end_node=dataset)]
            return

        self._coerce_numeric_columns(raw_df)

        placeholder = Dataset(
            id=dataset_id, data_type=data_type,
            access=experiment_obj.access,
            pre_processing_description=pre_processing, peri_processing_description=peri_processing,
            row_count=len(raw_df), column_count=len(raw_df.columns),
            gene_id_column=analyte_id_col, sample_columns=list(raw_df.columns)
        )
        run_biosample_edges = [
            DatasetRunBiosampleEdge(start_node=placeholder, end_node=self._ctx.run_biosample_by_original_id[str(col)])
            for col in raw_df.columns if str(col) in self._ctx.run_biosample_by_original_id
        ]
        analyte_edges = self._build_analyte_edges(placeholder, raw_df, DatasetGeneEdge, DatasetMetaboliteEdge, DatasetProteinEdge)

        qualified_df = self._qualify_data_frame_ids(raw_df)

        dataset = Dataset(
            id=dataset_id, data_type=data_type,
            access=experiment_obj.access,
            pre_processing_description=pre_processing, peri_processing_description=peri_processing,
            row_count=len(qualified_df), column_count=len(qualified_df.columns),
            gene_id_column=analyte_id_col, sample_columns=list(qualified_df.columns),
            _data_frame=qualified_df
        )
        for edge in run_biosample_edges + analyte_edges:
            edge.start_node = dataset

        yield [dataset, ExperimentDatasetEdge(start_node=experiment_obj, end_node=dataset)]
        if run_biosample_edges:
            yield run_biosample_edges
        if analyte_edges:
            yield analyte_edges

    def _stats_results_nodes(self, experiment_obj: Experiment, stats_parser: ExcelsheetParser,
                             stats_data: Optional[ParsedPounceData] = None,
                             stats_file: Optional[str] = None
                             ) -> Generator[List[Union[Node, Relationship]], None, None]:
        """Parse StatsReadyData and EffectSize from the stats results workbook."""
        stats_sheet_names = stats_parser.sheet_dfs.keys()

        analyte_id_col = None
        es_column_map: Dict[str, ComparisonColumn] = {}
        if StatsResultsWorkbook.EffectSizeMapSheet.name in stats_sheet_names:
            analyte_id_col, es_column_map = self._parse_effect_size_map(stats_parser)

        if StatsResultsWorkbook.StatsReadyDataSheet.name in stats_sheet_names:
            yield from self._parse_data_matrix(
                experiment_obj,
                meta_sheet=StatsResultsWorkbook.StatsResultsMetaSheet.name,
                data_sheet=StatsResultsWorkbook.StatsReadyDataSheet.name,
                analyte_id_col=analyte_id_col,
                data_type="stats ready data",
                parser=stats_parser
            )

        if StatsResultsWorkbook.EffectSizeSheet.name in stats_sheet_names:
            yield from self._parse_effect_size(
                experiment_obj,
                analyte_id_col,
                es_column_map,
                stats_parser,
                stats_data,
                stats_file=stats_file,
            )

    def _parse_effect_size(self, experiment_obj: Experiment,
                           analyte_id_col: str = None,
                           es_column_map: Dict[str, ComparisonColumn] = None,
                           parser: ExcelsheetParser = None,
                           stats_data: Optional[ParsedPounceData] = None,
                           stats_file: Optional[str] = None
                           ) -> Generator[List[Union[Node, Relationship]], None, None]:
        """Parse EffectSize sheet into a StatsResult node with edges to analytes."""
        if es_column_map is None:
            es_column_map = {}

        stats_result_id = f"{experiment_obj.id}:effect_size"
        meta_sheet = StatsResultsWorkbook.StatsResultsMetaSheet

        meta_kwargs = dict(
            name=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.statsresults_name),
            description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.stats_description),
            pre_processing_description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.pre_processing_description),
            peri_processing_description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.peri_processing_description),
            effect_size_description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.effect_size),
            effect_size_pval_description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.effect_size_pval),
            effect_size_adj_pval_description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.effect_size_adj_pval),
            effect_size_2_description=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.effect_size_2),
            data_analysis_code_link=parser.safe_get_string(meta_sheet.name, meta_sheet.Key.data_analysis_code_link),
        )

        parsed_stats = stats_data.stats_results[0] if stats_data and stats_data.stats_results else None
        informatician = Person.make(parsed_stats.lead_informatician, parsed_stats.lead_informatician_email) \
            if parsed_stats and parsed_stats.lead_informatician else None

        raw_df = self._prepare_data_frame(parser, StatsResultsWorkbook.EffectSizeSheet.name, analyte_id_col)

        if raw_df is None:
            stats_result = StatsResult(
                id=stats_result_id, data_type="effect_size",
                workbook=WorkbookArtifact(
                    original_filename=os.path.basename(stats_file) if stats_file else None,
                    _local_path=stats_file,
                ) if stats_file else None,
                row_count=0, column_count=0,
                analyte_id_column=analyte_id_col,
                comparison_columns=list(es_column_map.values()), **meta_kwargs
            )
            nodes = [stats_result, ExperimentStatsResultEdge(start_node=experiment_obj, end_node=stats_result)]
            if informatician:
                nodes += [informatician, StatsResultPersonEdge(start_node=stats_result, end_node=informatician, roles=["Informatician"])]
            yield nodes
            return

        self._coerce_numeric_columns(raw_df)

        placeholder = StatsResult(id=stats_result_id, data_type="effect_size",
                                  row_count=0, column_count=0, analyte_id_column=analyte_id_col)
        analyte_edges = self._build_analyte_edges(placeholder, raw_df, StatsResultGeneEdge, StatsResultMetaboliteEdge, StatsResultProteinEdge)

        qualified_df = self._qualify_data_frame_ids(raw_df)

        comparison_columns = [
            es_column_map.get(str(col), ComparisonColumn.from_column_name(str(col)))
            for col in qualified_df.columns
        ]

        stats_result = StatsResult(
            id=stats_result_id, data_type="effect_size",
            workbook=WorkbookArtifact(
                original_filename=os.path.basename(stats_file) if stats_file else None,
                _local_path=stats_file,
            ) if stats_file else None,
            row_count=len(qualified_df), column_count=len(qualified_df.columns),
            analyte_id_column=analyte_id_col, comparison_columns=comparison_columns,
            _data_frame=qualified_df, **meta_kwargs
        )
        for edge in analyte_edges:
            edge.start_node = stats_result

        nodes = [stats_result, ExperimentStatsResultEdge(start_node=experiment_obj, end_node=stats_result)]
        if informatician:
            nodes += [informatician, StatsResultPersonEdge(start_node=stats_result, end_node=informatician, roles=["Informatician"])]
        yield nodes
        if analyte_edges:
            yield analyte_edges

    # --- Conversion helpers (ParsedX -> Node) ---

    @staticmethod
    def _make_category_value(value: str, column_name: str = None) -> Optional[CategoryValue]:
        if value is None:
            return None
        col = column_name or "unknown"
        return CategoryValue(id=f"{col}-{value}", name=col, value=value)

    @staticmethod
    def _convert_demographics(parsed, parent_id: str) -> Optional[Demographics]:
        if parsed is None:
            return None
        age = parsed.age
        race = parsed.race
        ethnicity = parsed.ethnicity
        sex = parsed.sex
        categories_raw = parsed.demographic_categories or []
        phenotype_raw = parsed.phenotype_categories or []

        # Parse "column:value" strings into CategoryValue objects
        categories = []
        for cat_str in categories_raw:
            if ":" in cat_str:
                col, val = cat_str.split(":", 1)
                categories.append(CategoryValue(id=f"{col}-{val}", name=col, value=val))

        phenotype_categories = []
        for cat_str in phenotype_raw:
            if ":" in cat_str:
                col, val = cat_str.split(":", 1)
                phenotype_categories.append(CategoryValue(id=f"{col}-{val}", name=col, value=val))

        if all(v is None for v in [age, race, ethnicity, sex]) and not categories and not phenotype_categories:
            return None

        return Demographics(
            id=f"{parent_id}::demographics" if parent_id else None,
            age=age, race=race, ethnicity=ethnicity, sex=sex,
            categories=categories or None,
            phenotype_categories=phenotype_categories or None
        )

    @staticmethod
    def _convert_exposure(parsed) -> Exposure:
        """Convert a ParsedExposure into an Exposure Node."""
        return Exposure(
            id='calculate',
            names=parsed.names or [],
            type=parsed.type,
            category=CategoryValue(
                id=f"exposure_category-{parsed.category}",
                name="exposure_category", value=parsed.category
            ) if parsed.category else None,
            concentration=parsed.concentration,
            concentration_unit=parsed.concentration_unit,
            duration=parsed.duration,
            duration_unit=parsed.duration_unit,
            start_time=parsed.start_time,
            end_time=parsed.end_time,
            growth_media=parsed.growth_media,
            condition=CategoryValue(
                id=f"condition_category-{parsed.condition_category}",
                name="condition_category", value=parsed.condition_category
            ) if parsed.condition_category else None,
        )

    # --- Utilities ---

    @staticmethod
    def _parse_effect_size_map(parser: ExcelsheetParser) -> Tuple[Optional[str], Dict[str, ComparisonColumn]]:
        """Parse EffectSize_Map sheet into (analyte_id_col, {col_name: ComparisonColumn})."""
        sheet_name = StatsResultsWorkbook.EffectSizeMapSheet.name
        df = parser.sheet_dfs.get(sheet_name)
        if df is None or df.empty:
            return None, {}

        analyte_keys = {
            StatsResultsWorkbook.EffectSizeMapSheet.Key.gene_id,
            StatsResultsWorkbook.EffectSizeMapSheet.Key.metabolite_id,
            StatsResultsWorkbook.EffectSizeMapSheet.Key.protein_id,
        }
        notes_col = "Submitter_Notes"
        has_notes = notes_col in df.columns

        analyte_id_col = None
        column_map: Dict[str, ComparisonColumn] = {}

        for _, row in df.iterrows():
            ncats_key = row.get(ExcelsheetParser.KEY_COLUMN, '')
            if not isinstance(ncats_key, str):
                ncats_key = ''
            ncats_key = ncats_key.strip()

            submitter_val = row.get(ExcelsheetParser.MAPPED_VALUE_COLUMN, '')
            if not isinstance(submitter_val, str) or not submitter_val.strip():
                continue
            submitter_val = submitter_val.strip()

            notes = None
            if has_notes:
                raw_notes = row.get(notes_col)
                if isinstance(raw_notes, str) and raw_notes.strip():
                    notes = raw_notes.strip()

            if ncats_key in analyte_keys:
                analyte_id_col = submitter_val
            else:
                parts = [p.strip() for p in submitter_val.split('|') if p.strip()]
                if not parts:
                    continue
                col_name = parts[0]
                properties = parts[1:]
                if len(properties) == 0:
                    properties = None

                col = ComparisonColumn.from_column_name(col_name)
                if getattr(parser, "pounce_version", 2) == 3:
                    v3_role = re.fullmatch(r"(ES|Pval|AdjPval)_result\d+", ncats_key)
                    if v3_role:
                        role = v3_role.group(1)
                        col.stat_type = {
                            "ES": "effect_size",
                            "Pval": "pvalue",
                            "AdjPval": "adjusted_pvalue",
                        }[role]
                        prefix_key = {
                            "ES": "ES_1_prefix",
                            "Pval": "ES_1_Pval_prefix",
                            "AdjPval": "ES_1_adjPval_prefix",
                        }[role]
                        prefix = parser.safe_get_string("StatsResultsMeta", prefix_key)
                        if prefix and col_name.startswith(prefix):
                            remainder = col_name[len(prefix):].lstrip("_")
                            group1, group2 = ComparisonColumn.parse_groups(remainder)
                            col.group1, col.group2 = group1, group2
                            col.comparison = f"{group1}_VS_{group2}" if group1 and group2 else group1
                col.properties = properties
                col.notes = notes
                column_map[col_name] = col

        return analyte_id_col, column_map

    @staticmethod
    def _coerce_numeric_columns(df) -> None:
        """Convert all DataFrame columns to numeric in-place, coercing non-numeric values to NaN."""
        for col in df.columns:
            df[col] = pd.to_numeric(df[col], errors='coerce')

    @staticmethod
    def _split_pipe(value) -> List[str]:
        """Null-safe split of a pipe-delimited cell value. Returns [] for blank/NaN."""
        if not isinstance(value, str) or not value.strip():
            return []
        return [v.strip() for v in value.split('|') if v.strip() and v.strip().upper() not in _NA_VALUES]

    @staticmethod
    def _parse_int(value) -> Optional[int]:
        """Safely parse an integer value."""
        if value is None or value == '' or value == 'NA':
            return None
        try:
            return int(float(value))
        except (ValueError, TypeError):
            return None
