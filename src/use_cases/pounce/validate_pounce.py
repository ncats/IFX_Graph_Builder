import argparse
import os
from src.core.config import ETL_Config, create_object_from_config
from src.input_adapters.pounce_sheets.mapping_coverage import (
    check_metabolite_coverage, check_gene_coverage, check_protein_coverage,
)
from src.input_adapters.pounce_sheets.pounce_parser import (
    _PROJECT_REQUIRED_SHEETS, _EXPERIMENT_RECOGNIZED_SHEETS, _STATS_RECOGNIZED_SHEETS,
)
from src.input_adapters.pounce_sheets.workbook_format import restore_pounce_issue

parser = argparse.ArgumentParser()
parser.add_argument(
    "--show-unmapped",
    action="store_true",
    help="Print the full list of unmapped IDs in coverage output."
)
args = parser.parse_args()
UNMAPPED_PREVIEW_LIMIT = 50

# yaml_file = "./src/use_cases/pounce/test_pounce_validation.yaml"
yaml_file = "./src/use_cases/pounce/pounce.yaml"
# yaml_file = "./src/use_cases/working.yaml"

config = ETL_Config(yaml_file)
resolver_map = {t: r for r in config.resolvers.values() for t in r.types}

pounce_configs = [
    c for c in config.config_dict.get("input_adapters", [])
    if c.get("class") == "PounceInputAdapter"
]


def _assign_content_issue_source_files(all_issues, project_file: str, experiment_files: list[str], stats_files: list[str]) -> None:
    for issue in all_issues:
        if issue.source_file:
            continue

        if issue.entity == "project":
            issue.source_file = project_file
            continue

        if issue.entity == "experiments" and issue.row is not None and 0 <= issue.row < len(experiment_files):
            issue.source_file = experiment_files[issue.row]
            continue

        if issue.entity == "stats_results" and issue.row is not None and 0 <= issue.row < len(stats_files):
            issue.source_file = stats_files[issue.row]
            continue

        if issue.sheet in _PROJECT_REQUIRED_SHEETS:
            issue.source_file = project_file
        elif issue.sheet in _EXPERIMENT_RECOGNIZED_SHEETS and len(experiment_files) == 1:
            issue.source_file = experiment_files[0]
        elif issue.sheet in _STATS_RECOGNIZED_SHEETS and len(stats_files) == 1:
            issue.source_file = stats_files[0]

for chosen in pounce_configs:
    kwargs = chosen.get("kwargs", {})
    project_file = kwargs.get("project_file", "")
    experiment_files = kwargs.get("experiment_files") or []
    stats_files = kwargs.get("stats_results_files") or []

    print(f"\n{'='*60}")
    print(f"Validating: {os.path.basename(project_file)}")
    print(f"{'='*60}")

    adapter = create_object_from_config(chosen)
    parsed_data = adapter.get_validation_data()
    structural_issues = adapter.get_structural_issues()

    content_issues = []
    for validator in adapter.get_validators():
        content_issues.extend(validator.validate(parsed_data))
    _assign_content_issue_source_files(content_issues, project_file, experiment_files, stats_files)
    for issue in content_issues:
        if issue.source_file and issue.sheet:
            restore_pounce_issue(
                issue,
                parsed_data.original_sheets.get(issue.source_file, {}),
                parsed_data.workbook_versions.get(issue.source_file, 2),
            )

    all_issues = structural_issues + content_issues
    errors = [e for e in all_issues if e.severity == "error"]
    warnings = [e for e in all_issues if e.severity == "warning"]

    def _source_label(issue) -> str:
        """Return a display label (basename) for the xlsx file an issue came from."""
        if issue.source_file:
            return os.path.basename(issue.source_file)
        issue_sheet = issue.sheet
        if issue_sheet in _PROJECT_REQUIRED_SHEETS:
            return os.path.basename(project_file) if project_file else "Project workbook"
        if issue_sheet in _EXPERIMENT_RECOGNIZED_SHEETS:
            if len(experiment_files) == 1:
                return os.path.basename(experiment_files[0])
            return "Experiment workbook(s): " + ", ".join(os.path.basename(f) for f in experiment_files)
        if issue_sheet in _STATS_RECOGNIZED_SHEETS:
            if len(stats_files) == 1:
                return os.path.basename(stats_files[0])
            return "Stats results workbook(s): " + ", ".join(os.path.basename(f) for f in stats_files)
        return "(unknown workbook)"

    if all_issues:
        error_count = len(errors)
        warning_count = len(warnings)
        parts = []
        if error_count:
            parts.append(f"{error_count} error(s)")
        if warning_count:
            parts.append(f"{warning_count} warning(s)")
        print(f"\n{' and '.join(parts)} found:")

        # Group: workbook → sheet → issues
        by_workbook: dict = {}
        for item in all_issues:
            wb = _source_label(item)
            sheet = item.sheet or "(no sheet)"
            by_workbook.setdefault(wb, {}).setdefault(sheet, []).append(item)

        for workbook, sheets in by_workbook.items():
            print(f"\n{workbook}:")
            for sheet, group in sheets.items():
                print(f"  {sheet}:")
                for item in group:
                    label = "ERROR" if item.severity == "error" else "WARN "
                    row_info = f" row {item.row}" if item.row is not None else ""
                    print(f"    {label}  {item.field}{row_info}: {item.message}")
    else:
        print("\nValidation passed — no errors or warnings.")

    print(f"\nSummary:")
    print("  Workbook formats:")
    for path, version in parsed_data.workbook_versions.items():
        print(f"    {os.path.basename(path)}: v{version}")
    print(f"  Project: {parsed_data.project.project_name if parsed_data.project else 'N/A'}")
    print(f"  Biosamples: {len(parsed_data.biosamples)}")
    print(f"  Biospecimens: {len(parsed_data.biospecimens)}")
    print(f"  Experiments: {len(parsed_data.experiments)}")
    print(f"  RunBiosamples: {len(parsed_data.run_biosamples)}")
    print(f"  StatsResults: {len(parsed_data.stats_results)}")
    print(f"  People: {len(parsed_data.people)}")

    has_metabolites = bool(parsed_data.metabolites)
    has_genes = bool(parsed_data.genes)
    has_proteins = bool(parsed_data.proteins)

    if has_metabolites or has_genes or has_proteins:
        print(f"\nMapping Coverage:")

        if has_metabolites:
            metab_resolver = resolver_map.get("Metabolite")
            if metab_resolver:
                cov = check_metabolite_coverage(parsed_data.metabolites, metab_resolver)
                print(f"  MetabMeta — {cov.mapped}/{cov.total} ({cov.mapped_pct:.1f}%)"
                      f" metabolites will resolve to canonical Metabolite nodes")
                if cov.unmapped_ids:
                    print(f"    Unmapped: {len(cov.unmapped_ids)}")
                    preview_ids = cov.unmapped_ids[:UNMAPPED_PREVIEW_LIMIT]
                    if args.show_unmapped:
                        print(f"    Unmapped IDs: {', '.join(cov.unmapped_ids)}")
                    else:
                        suffix = "" if len(cov.unmapped_ids) <= UNMAPPED_PREVIEW_LIMIT else ", ..."
                        print(f"    Unmapped preview: {', '.join(preview_ids)}{suffix}")
            else:
                print(f"  MetabMeta — no Metabolite resolver configured")

        if has_genes:
            gene_resolver = resolver_map.get("Gene")
            cov = check_gene_coverage(parsed_data.genes, gene_resolver)
            if cov:
                print(f"  GeneMeta — {cov.mapped}/{cov.total} ({cov.mapped_pct:.1f}%)"
                      f" genes will resolve to canonical Gene nodes")
                if cov.unmapped_ids:
                    print(f"    Unmapped: {len(cov.unmapped_ids)}")
                    preview_ids = cov.unmapped_ids[:UNMAPPED_PREVIEW_LIMIT]
                    if args.show_unmapped:
                        print(f"    Unmapped IDs: {', '.join(cov.unmapped_ids)}")
                    else:
                        suffix = "" if len(cov.unmapped_ids) <= UNMAPPED_PREVIEW_LIMIT else ", ..."
                        print(f"    Unmapped preview: {', '.join(preview_ids)}{suffix}")
            else:
                print(f"  GeneMeta — no Gene resolver configured")

        if has_proteins:
            protein_resolver = resolver_map.get("Protein")
            cov = check_protein_coverage(parsed_data.proteins, protein_resolver)
            if cov:
                print(f"  ProteinMeta — {cov.mapped}/{cov.total} ({cov.mapped_pct:.1f}%)"
                      f" proteins will resolve to canonical Protein nodes")
                if cov.unmapped_ids:
                    print(f"    Unmapped: {len(cov.unmapped_ids)}")
                    preview_ids = cov.unmapped_ids[:UNMAPPED_PREVIEW_LIMIT]
                    if args.show_unmapped:
                        print(f"    Unmapped IDs: {', '.join(cov.unmapped_ids)}")
                    else:
                        suffix = "" if len(cov.unmapped_ids) <= UNMAPPED_PREVIEW_LIMIT else ", ..."
                        print(f"    Unmapped preview: {', '.join(preview_ids)}{suffix}")
            else:
                print(f"  ProteinMeta — no Protein resolver configured")
