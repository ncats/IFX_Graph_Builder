"""POUNCE submission validation routes.

Wired into app.py via:
    app.include_router(_pounce_module.router)
    _pounce_module.set_templates(templates)
    _pounce_module.set_pounce_config(path_to_pounce_yaml)
    _pounce_module.set_smtp_config(path_to_smtp_yaml)   # optional
"""
import os
import shutil
import smtplib
import tempfile
import time
import uuid
from dataclasses import asdict
from email import encoders
from email.mime.base import MIMEBase
from email.mime.multipart import MIMEMultipart
from email.mime.text import MIMEText
from pathlib import Path
from typing import List

import pandas as pd
from fastapi import APIRouter, File, Form, HTTPException, Request, UploadFile
from fastapi.responses import FileResponse, HTMLResponse
from sqlalchemy import text

from src.core.validator import ValidationError
from src.input_adapters.excel_sheet_adapter import ExcelsheetParser
from src.input_adapters.pounce_sheets.mapping_coverage import (
    check_gene_coverage,
    check_metabolite_coverage,
    check_protein_coverage,
)
from src.input_adapters.pounce_sheets.pounce_input_adapter import PounceInputAdapter
from src.input_adapters.pounce_sheets.pounce_parser import (
    PounceParser,
    _EXPERIMENT_RECOGNIZED_SHEETS,
    _PROJECT_REQUIRED_SHEETS,
    _STATS_RECOGNIZED_SHEETS,
)
from src.input_adapters.pounce_sheets.validator_loader import load_validators
from src.input_adapters.pounce_sheets.workbook_format import (
    UnsupportedPounceFormat, normalize_pounce_workbook, restore_pounce_issue,
)
from src.qa_browser.pounce_api_contract import PounceValidationRequest

router = APIRouter(prefix="/pounce")

_templates = None
_VALIDATORS_CONFIG = str(Path(__file__).resolve().parent.parent / "use_cases" / "pounce" / "pounce_validators.yaml")

# Resolver map populated at startup from pounce.yaml.
# Keys are entity type strings (e.g. "Metabolite", "Gene").
# Empty if pounce.yaml is missing or all resolvers fail to load.
_resolver_map: dict = {}
_pounce_config_path: str = ""
_resolvers_loaded: bool = False

# ── Session store ────────────────────────────────────────────────────────────
# Keeps uploaded files alive after validation so the submitter can email them.
# session dict keys:
#   dir          — temp directory path (NOT auto-deleted)
#   created      — unix timestamp for TTL
#   project_name — display name from parsed project or uploaded filename
#   file_paths   — list of all uploaded xlsx paths
#   file_basenames — list of original basenames (for display)
#   summary      — dict: error_count, warning_count, coverage_lines

_sessions: dict = {}
_SESSION_TTL = 3600  # 1 hour

# ── SMTP configuration ───────────────────────────────────────────────────────
# Empty dict means email is disabled.
# Keys: host, port, user, password, from_address, to_address, use_tls

_smtp_config: dict = {}
_get_mysql_db_engine = None
_MYSQL_DB_NAME = "omicsdb_dev2"
_MYSQL_SOURCE_ID = "default"
_public_project_base_url = ""


def set_templates(t):
    global _templates
    _templates = t


def set_smtp_config(path: str):
    """Load SMTP credentials from a YAML file.

    Called once at app startup. If the file is missing or malformed the
    email submission feature is silently disabled.
    """
    global _smtp_config
    if not path or not os.path.exists(path):
        print(f"[pounce] SMTP config not found at {path!r} — email submission disabled")
        return
    try:
        import yaml
        with open(path) as f:
            _smtp_config = yaml.safe_load(f) or {}
        print(f"[pounce] Loaded SMTP config from {path!r}")
    except Exception as e:
        print(f"[pounce] Could not load SMTP config from {path!r}: {e}")


def set_mysql_engine_getter(getter):
    global _get_mysql_db_engine
    _get_mysql_db_engine = getter


def set_public_project_base_url(url: str):
    global _public_project_base_url
    _public_project_base_url = (url or "").rstrip("/")


def set_pounce_config(pounce_yaml_path: str):
    """Store the pounce.yaml path for lazy resolver loading."""
    global _pounce_config_path, _resolver_map, _resolvers_loaded
    _pounce_config_path = pounce_yaml_path or ""
    _resolver_map = {}
    _resolvers_loaded = False
    if not _pounce_config_path or not os.path.exists(_pounce_config_path):
        print(f"[pounce] config not found at {_pounce_config_path!r} — mapping coverage disabled")
    else:
        print(f"[pounce] configured lazy resolver loading from {_pounce_config_path!r}")


def _ensure_resolvers_loaded():
    """Load resolvers from pounce.yaml on first use for mapping coverage checks."""
    global _resolver_map, _resolvers_loaded
    if _resolvers_loaded:
        return
    _resolvers_loaded = True
    if not _pounce_config_path or not os.path.exists(_pounce_config_path):
        return
    try:
        from src.core.config import Config, create_object_from_config
        config = Config(_pounce_config_path)
        for rc in config.config_dict.get("resolvers", []):
            label = rc.get("label", "?")
            try:
                resolver = create_object_from_config(rc)
                for t in getattr(resolver, "types", []):
                    _resolver_map[t] = resolver
                    print(f"[pounce] loaded resolver '{label}' for type '{t}'")
            except Exception as e:
                print(f"[pounce] could not load resolver '{label}': {e}")
    except Exception as e:
        print(f"[pounce] could not load config from {_pounce_config_path!r}: {e}")


# ── Session helpers ───────────────────────────────────────────────────────────

def _cleanup_old_sessions():
    """Remove sessions older than _SESSION_TTL and delete their temp dirs."""
    now = time.time()
    expired = [sid for sid, s in list(_sessions.items()) if now - s["created"] > _SESSION_TTL]
    for sid in expired:
        shutil.rmtree(_sessions[sid]["dir"], ignore_errors=True)
        del _sessions[sid]


# ── File helpers ──────────────────────────────────────────────────────────────

def _save_upload(upload: UploadFile, directory: str, prefix: str = "") -> str:
    safe_name = Path(upload.filename).name
    dest = os.path.join(directory, f"{prefix}{safe_name}" if prefix else safe_name)
    if os.path.exists(dest):
        base, ext = os.path.splitext(safe_name)
        dest = os.path.join(directory, f"{prefix}{base}_1{ext}")
    with open(dest, "wb") as f:
        f.write(upload.file.read())
    return dest


def _write_blob_file(filename: str, content: bytes, directory: str, prefix: str = "") -> str:
    safe_name = Path(filename).name
    dest = os.path.join(directory, f"{prefix}{safe_name}" if prefix else safe_name)
    if os.path.exists(dest):
        base, ext = os.path.splitext(safe_name)
        dest = os.path.join(directory, f"{prefix}{base}_1{ext}")
    with open(dest, "wb") as f:
        f.write(content)
    return dest


def _is_real_upload(upload: UploadFile | None) -> bool:
    return bool(upload and getattr(upload, "filename", ""))


def _build_session(
    project_path: str,
    exp_paths: list[str],
    stats_paths: list[str],
    parsed_data=None,
    submission_label: str = "POUNCE submission",
    file_change_summary: dict | None = None,
    existing_project_id: str | None = None,
) -> tuple[str, dict]:
    session_id = str(uuid.uuid4())
    project_name = (
        parsed_data.project.project_name
        if parsed_data and parsed_data.project
        else os.path.basename(project_path)
    )
    session = {
        "dir": os.path.dirname(project_path),
        "created": time.time(),
        "project_name": project_name,
        "project_path": project_path,
        "exp_paths": exp_paths,
        "stats_paths": stats_paths,
        "file_paths": [project_path] + exp_paths + stats_paths,
        "file_basenames": (
            [os.path.basename(project_path)]
            + [os.path.basename(p) for p in exp_paths]
            + [os.path.basename(p) for p in stats_paths]
        ),
        "summary": None,
        "submission_label": submission_label,
        "file_change_summary": file_change_summary or {},
        "existing_project_id": existing_project_id,
    }
    _sessions[session_id] = session
    return session_id, session


def _summarize_new_submission(project_path: str, exp_paths: list[str], stats_paths: list[str]) -> dict:
    return {
        "added": [os.path.basename(project_path)] + [os.path.basename(p) for p in exp_paths] + [os.path.basename(p) for p in stats_paths],
        "replaced": [],
        "unchanged": [],
    }


def _summarize_edit_submission(
    prior_session: dict,
    project_file: UploadFile,
    experiment_files: list[UploadFile],
    stats_files: list[UploadFile],
    new_experiment_files: list[UploadFile],
    new_stats_files: list[UploadFile],
) -> dict:
    summary = {"added": [], "replaced": [], "unchanged": []}

    current_project_name = os.path.basename(prior_session["project_path"])
    if _is_real_upload(project_file):
        summary["replaced"].append(f"{current_project_name} -> {Path(project_file.filename).name}")
    else:
        summary["unchanged"].append(current_project_name)

    for existing_path, upload in zip(prior_session["exp_paths"], experiment_files):
        existing_name = os.path.basename(existing_path)
        if _is_real_upload(upload):
            summary["replaced"].append(f"{existing_name} -> {Path(upload.filename).name}")
        else:
            summary["unchanged"].append(existing_name)

    for existing_path, upload in zip(prior_session["stats_paths"], stats_files):
        existing_name = os.path.basename(existing_path)
        if _is_real_upload(upload):
            summary["replaced"].append(f"{existing_name} -> {Path(upload.filename).name}")
        else:
            summary["unchanged"].append(existing_name)

    pair_count = max(len(new_experiment_files), len(new_stats_files))
    for i in range(pair_count):
        exp_upload = new_experiment_files[i] if i < len(new_experiment_files) else None
        stats_upload = new_stats_files[i] if i < len(new_stats_files) else None
        if _is_real_upload(exp_upload):
            summary["added"].append(Path(exp_upload.filename).name)
        if _is_real_upload(stats_upload):
            summary["added"].append(Path(stats_upload.filename).name)

    return summary


def _get_mysql_engine():
    if _get_mysql_db_engine is None:
        raise RuntimeError("MySQL access is not configured for the QA browser.")
    return _get_mysql_db_engine(_MYSQL_DB_NAME, source_id=_MYSQL_SOURCE_ID)


def _search_existing_projects(search_term: str = "") -> list[dict]:
    engine = _get_mysql_engine()
    where_sql = ""
    params: dict = {}
    if search_term:
        where_sql = """
        WHERE (
            p.name LIKE :search
            OR p.id LIKE :search
            OR p.description LIKE :search
            OR p.access LIKE :search
            OR COALESCE(lg.lab_groups, '') LIKE :search
            OR COALESCE(kw.keywords, '') LIKE :search
            OR COALESCE(pt.project_types, '') LIKE :search
            OR COALESCE(pp.people, '') LIKE :search
            OR COALESCE(ep.experiment_people, '') LIKE :search
        )
        """
        params["search"] = f"%{search_term}%"

    query = text(f"""
        SELECT
            p.id,
            p.name,
            p.description,
            p.date,
            p.access,
            p.rare_disease_focus,
            COALESCE(lg.lab_groups, '') AS lab_groups,
            COALESCE(kw.keywords, '') AS keywords,
            COALESCE(pt.project_types, '') AS project_types,
            COALESCE(pp.people, '') AS people,
            COALESCE(ep.experiment_people, '') AS experiment_people,
            COUNT(DISTINCT pe.to_id) AS experiment_count,
            COUNT(DISTINCT pb.to_id) AS biosample_count,
            COUNT(DISTINCT CASE WHEN pw.project_id IS NOT NULL THEN p.id END) AS has_project_workbook
        FROM project p
        LEFT JOIN (
            SELECT plg.parent_id, GROUP_CONCAT(lg.label ORDER BY lg.label SEPARATOR ', ') AS lab_groups
            FROM project__lab_groups plg
            JOIN lab_group lg ON lg.id = plg.lab_group_id
            GROUP BY plg.parent_id
        ) lg ON lg.parent_id = p.id
        LEFT JOIN (
            SELECT parent_id, GROUP_CONCAT(value ORDER BY value SEPARATOR ', ') AS keywords
            FROM project__keywords
            GROUP BY parent_id
        ) kw ON kw.parent_id = p.id
        LEFT JOIN (
            SELECT ppt.parent_id, GROUP_CONCAT(pt.label ORDER BY pt.label SEPARATOR ', ') AS project_types
            FROM project__project_type ppt
            JOIN project_type pt ON pt.id = ppt.project_type_id
            GROUP BY ppt.parent_id
        ) pt ON pt.parent_id = p.id
        LEFT JOIN (
            SELECT
                pp.from_id AS project_id,
                GROUP_CONCAT(
                    DISTINCT COALESCE(NULLIF(TRIM(p.name), ''), p.email, p.id)
                    ORDER BY COALESCE(NULLIF(TRIM(p.name), ''), p.email, p.id)
                    SEPARATOR ', '
                ) AS people
            FROM project_to_person pp
            JOIN person p ON p.id = pp.to_id
            GROUP BY pp.from_id
        ) pp ON pp.project_id = p.id
        LEFT JOIN (
            SELECT
                pe.from_id AS project_id,
                GROUP_CONCAT(
                    DISTINCT COALESCE(NULLIF(TRIM(p.name), ''), p.email, p.id)
                    ORDER BY COALESCE(NULLIF(TRIM(p.name), ''), p.email, p.id)
                    SEPARATOR ', '
                ) AS experiment_people
            FROM project_to_experiment pe
            JOIN experiment_to_person ep ON ep.from_id = pe.to_id
            JOIN person p ON p.id = ep.to_id
            GROUP BY pe.from_id
        ) ep ON ep.project_id = p.id
        LEFT JOIN project_to_experiment pe ON pe.from_id = p.id
        LEFT JOIN project_to_biosample pb ON pb.from_id = p.id
        LEFT JOIN project__workbook pw ON pw.project_id = p.id
        {where_sql}
        GROUP BY
            p.id, p.name, p.description, p.date, p.access, p.rare_disease_focus,
            lg.lab_groups, kw.keywords, pt.project_types, pp.people, ep.experiment_people
        ORDER BY p.date DESC, p.name ASC
        LIMIT 100
    """)
    with engine.connect() as conn:
        return [dict(row) for row in conn.execute(query, params).mappings().all()]


def _create_existing_project_session(project_id: str) -> str:
    engine = _get_mysql_engine()
    tmp = tempfile.mkdtemp()
    try:
        with engine.connect() as conn:
            project_row = conn.execute(text("""
                SELECT p.id, p.name, pw.original_filename, pw.content_blob
                FROM project p
                LEFT JOIN project__workbook pw ON pw.project_id = p.id
                WHERE p.id = :project_id
            """), {"project_id": project_id}).mappings().first()

            if not project_row:
                raise RuntimeError(f"Project '{project_id}' not found in {_MYSQL_DB_NAME}.")
            if not project_row["original_filename"] or project_row["content_blob"] is None:
                raise RuntimeError(
                    f"Project '{project_id}' does not have a stored project workbook in {_MYSQL_DB_NAME}."
                )

            pair_rows = conn.execute(text("""
                SELECT
                    e.id AS experiment_id,
                    e.name AS experiment_name,
                    ew.original_filename AS experiment_filename,
                    ew.content_blob AS experiment_blob,
                    sr.id AS stats_result_id,
                    sr.name AS stats_result_name,
                    sw.original_filename AS stats_filename,
                    sw.content_blob AS stats_blob
                FROM project_to_experiment pe
                JOIN experiment e ON e.id = pe.to_id
                LEFT JOIN experiment__workbook ew ON ew.experiment_id = e.id
                LEFT JOIN experiment_to_stats_result es ON es.from_id = e.id
                LEFT JOIN stats_result sr ON sr.id = es.to_id
                LEFT JOIN stats_result__workbook sw ON sw.stats_result_id = sr.id
                WHERE pe.from_id = :project_id
                ORDER BY e.name, e.id, sr.name, sr.id
            """), {"project_id": project_id}).mappings().all()

        project_path = _write_blob_file(project_row["original_filename"], project_row["content_blob"], tmp)
        exp_paths: list[str] = []
        stats_paths: list[str] = []

        for i, row in enumerate(pair_rows):
            if not row["experiment_filename"] or row["experiment_blob"] is None:
                raise RuntimeError(
                    f"Experiment '{row['experiment_id']}' does not have a stored workbook in {_MYSQL_DB_NAME}."
                )
            if not row["stats_result_id"] or not row["stats_filename"] or row["stats_blob"] is None:
                raise RuntimeError(
                    f"Experiment '{row['experiment_id']}' is missing its stored stats workbook in {_MYSQL_DB_NAME}."
                )
            exp_paths.append(_write_blob_file(row["experiment_filename"], row["experiment_blob"], tmp, prefix=f"exp{i}_"))
            stats_paths.append(_write_blob_file(row["stats_filename"], row["stats_blob"], tmp, prefix=f"stats{i}_"))

        session_id, session = _build_session(
            project_path,
            exp_paths,
            stats_paths,
            submission_label="Existing project edit",
            file_change_summary={
                "added": [],
                "replaced": [],
                "unchanged": [os.path.basename(project_path)] + [os.path.basename(p) for p in exp_paths] + [os.path.basename(p) for p in stats_paths],
            },
            existing_project_id=project_id,
        )
        session["project_name"] = project_row["name"] or session["project_name"]
        return session_id
    except Exception:
        shutil.rmtree(tmp, ignore_errors=True)
        raise


def _append_new_pairs(
    experiment_uploads: list[UploadFile],
    stats_uploads: list[UploadFile],
    tmp: str,
    exp_paths: list[str],
    stats_paths: list[str],
) -> None:
    pair_count = max(len(experiment_uploads), len(stats_uploads))
    for i in range(pair_count):
        exp_upload = experiment_uploads[i] if i < len(experiment_uploads) else None
        stats_upload = stats_uploads[i] if i < len(stats_uploads) else None
        has_exp = _is_real_upload(exp_upload)
        has_stats = _is_real_upload(stats_upload)
        if not has_exp and not has_stats:
            continue
        if has_exp != has_stats:
            raise ValueError(
                f"New experiment row {i + 1} must include both an experiment workbook and a StatsResults workbook."
            )
        exp_paths.append(_save_upload(exp_upload, tmp, prefix=f"exp_new{i}_"))
        stats_paths.append(_save_upload(stats_upload, tmp, prefix=f"stats_new{i}_"))


# ── Validation helpers ────────────────────────────────────────────────────────

def _source_label(issue, project_filename: str, exp_filenames: list, stats_filenames: list) -> str:
    if issue.source_file:
        return os.path.basename(issue.source_file)
    sheet = issue.sheet
    if sheet in _PROJECT_REQUIRED_SHEETS:
        return project_filename
    if sheet in _EXPERIMENT_RECOGNIZED_SHEETS:
        return exp_filenames[0] if len(exp_filenames) == 1 else "Experiment workbook(s)"
    if sheet in _STATS_RECOGNIZED_SHEETS:
        return stats_filenames[0] if len(stats_filenames) == 1 else "Stats workbook(s)"
    return "(unknown workbook)"


def _group_issues(all_issues, project_filename, exp_filenames, stats_filenames):
    grouped: dict = {}
    for issue in all_issues:
        wb = _source_label(issue, project_filename, exp_filenames, stats_filenames)
        sheet = issue.sheet or "(no sheet)"
        grouped.setdefault(wb, {}).setdefault(sheet, []).append(issue)
    return [(wb, list(sheets.items())) for wb, sheets in grouped.items()]


def _assign_content_issue_source_files(all_issues: list[ValidationError], project_path: str,
                                       exp_paths: list[str], stats_paths: list[str]) -> None:
    for issue in all_issues:
        if issue.source_file:
            continue

        if issue.entity == "project":
            issue.source_file = project_path
            continue

        if issue.entity == "experiments" and issue.row is not None and 0 <= issue.row < len(exp_paths):
            issue.source_file = exp_paths[issue.row]
            continue

        if issue.entity == "stats_results" and issue.row is not None and 0 <= issue.row < len(stats_paths):
            issue.source_file = stats_paths[issue.row]
            continue

        if issue.sheet in _PROJECT_REQUIRED_SHEETS:
            issue.source_file = project_path
        elif issue.sheet in _EXPERIMENT_RECOGNIZED_SHEETS and len(exp_paths) == 1:
            issue.source_file = exp_paths[0]
        elif issue.sheet in _STATS_RECOGNIZED_SHEETS and len(stats_paths) == 1:
            issue.source_file = stats_paths[0]


def _restore_content_issue_sheets(issues: list[ValidationError], parsed_data) -> None:
    for issue in issues:
        if issue.source_file and issue.sheet:
            original = parsed_data.original_sheets.get(issue.source_file, {})
            version = parsed_data.workbook_versions.get(issue.source_file, 2)
            restore_pounce_issue(issue, original, version)


def _file_version_rows(parsed_data, project_path: str, exp_paths: list[str], stats_paths: list[str]):
    if parsed_data is None:
        return []
    rows = []
    for role, paths in (("Project", [project_path]), ("Experiment", exp_paths),
                        ("Stats Results", stats_paths)):
        for path in paths:
            version = parsed_data.workbook_versions.get(path)
            if version is not None:
                rows.append((role, os.path.basename(path), f"v{version}"))
    return rows


def _validate_parsed_submission(parsed_data, structural_issues, project_path, exp_paths, stats_paths):
    """Apply the same content checks and issue presentation to either input transport."""
    content_issues = []
    for validator in load_validators(_VALIDATORS_CONFIG):
        content_issues.extend(validator.validate(parsed_data))
    _assign_content_issue_source_files(content_issues, project_path, exp_paths, stats_paths)
    _restore_content_issue_sheets(content_issues, parsed_data)
    issues = _collapse_dependent_mapping_issues(structural_issues + content_issues)
    grouped = _group_issues(
        issues, os.path.basename(project_path),
        [os.path.basename(path) for path in exp_paths],
        [os.path.basename(path) for path in stats_paths],
    )
    return issues, grouped


def _sheet_workbook_loader(workbooks):
    """Adapt caller-supplied sheet tables to the existing workbook parser."""
    by_source_id = {workbook.source_id: workbook for workbook in workbooks}

    def load(source_id, role):
        workbook = by_source_id[source_id]
        sheets = {
            sheet.name: pd.DataFrame(
                [["" if value is None else value for value in row] for row in sheet.rows],
                columns=sheet.columns,
            )
            for sheet in workbook.sheets
        }
        parser = ExcelsheetParser.from_sheets(source_id, sheets)
        return normalize_pounce_workbook(parser, role)

    return load


def _collapse_dependent_mapping_issues(all_issues: list[ValidationError]) -> list[ValidationError]:
    """Suppress downstream analyte-ID noise when the root map-sheet key is broken.

    If GeneMap.gene_id or MetabMap.metab_id is invalid, row-level GeneMeta/MetabMeta
    ID failures and downstream matrix cross-reference failures are consequences, not
    independent root causes. Replace that cascade with one explanatory message.
    """
    collapse_specs = [
        {
            "map_sheet": "GeneMap",
            "map_key": "gene_id",
            "meta_sheet": "GeneMeta",
            "meta_field": "gene_id",
            "xref_phrase": "does not exist in GeneMeta",
            "message": (
                "Downstream gene_id checks were skipped because 'GeneMap' maps "
                "'gene_id' incorrectly. Fix the GeneMap gene_id column mapping first; "
                "additional GeneMeta or matrix ID errors may appear after that."
            ),
        },
        {
            "map_sheet": "MetabMap",
            "map_key": "metab_id",
            "meta_sheet": "MetabMeta",
            "meta_field": "metab_id",
            "xref_phrase": "does not exist in MetabMeta",
            "message": (
                "Downstream metab_id checks were skipped because 'MetabMap' maps "
                "'metab_id' incorrectly. Fix the MetabMap metab_id column mapping first; "
                "additional MetabMeta or matrix ID errors may appear after that."
            ),
        },
    ]

    collapsed_issues = list(all_issues)
    for spec in collapse_specs:
        has_root_mapping_issue = any(
            issue.severity == "error"
            and issue.sheet == spec["map_sheet"]
            and (issue.field == spec["map_key"] or issue.column == spec["map_key"])
            for issue in collapsed_issues
        )
        if not has_root_mapping_issue:
            continue

        filtered = []
        suppressed_count = 0
        for issue in collapsed_issues:
            is_meta_fallout = issue.sheet == spec["meta_sheet"] and issue.field == spec["meta_field"]
            is_xref_fallout = spec["xref_phrase"] in issue.message
            if is_meta_fallout or is_xref_fallout:
                suppressed_count += 1
                continue
            filtered.append(issue)

        if suppressed_count:
            filtered.append(ValidationError(
                severity="warning",
                entity="parse",
                field=spec["map_key"],
                message=f"{spec['message']} ({suppressed_count} dependent issue{'s' if suppressed_count != 1 else ''} hidden.)",
                sheet=spec["map_sheet"],
                column=spec["map_key"],
            ))
        collapsed_issues = filtered

    return collapsed_issues


def _compute_coverage(parsed_data) -> list:
    """Return AnalyteCoverage results for analyte types that have a resolver loaded."""
    coverage = []
    if parsed_data is None:
        return coverage
    _ensure_resolvers_loaded()

    if parsed_data.metabolites:
        resolver = _resolver_map.get("Metabolite")
        if resolver:
            try:
                coverage.append(check_metabolite_coverage(parsed_data.metabolites, resolver))
            except Exception as e:
                print(f"[pounce] metabolite coverage check failed: {e}")

    if parsed_data.genes:
        resolver = _resolver_map.get("Gene")
        try:
            cov = check_gene_coverage(parsed_data.genes, resolver)
            if cov:
                coverage.append(cov)
        except Exception as e:
            print(f"[pounce] gene coverage check failed: {e}")

    if parsed_data.proteins:
        resolver = _resolver_map.get("Protein")
        try:
            cov = check_protein_coverage(parsed_data.proteins, resolver)
            if cov:
                coverage.append(cov)
        except Exception as e:
            print(f"[pounce] protein coverage check failed: {e}")

    return coverage


def _compute_summary(
    all_issues: list,
    coverage: list,
    grouped: list | None = None,
    parsed_data=None,
    metab_resolver_missing: bool = False,
    gene_resolver_missing: bool = False,
    protein_resolver_missing: bool = False,
) -> dict:
    """Build a summary dict for inclusion in the submission email."""
    errors   = [i for i in all_issues if i.severity == "error"]
    warnings = [i for i in all_issues if i.severity == "warning"]
    cov_lines = [
        f"  {cov.analyte_type}: {cov.mapped}/{cov.total} ({cov.mapped_pct:.1f}%)"
        for cov in coverage
    ]
    if metab_resolver_missing:
        cov_lines.append("  MetabMeta: resolver not configured, coverage check skipped")
    if gene_resolver_missing:
        cov_lines.append("  GeneMeta: resolver not configured, coverage check skipped")
    if protein_resolver_missing:
        cov_lines.append("  ProteinMeta: resolver not configured, coverage check skipped")

    issue_breakdown = []
    for workbook, sheets in grouped or []:
        workbook_count = sum(len(issues) for _, issues in sheets)
        issue_breakdown.append(f"  {workbook}: {workbook_count} issue{'s' if workbook_count != 1 else ''}")
        for sheet, issues in sheets:
            issue_breakdown.append(f"    - {sheet}: {len(issues)}")

    data_counts = []
    if parsed_data is not None:
        data_counts = [
            f"  Experiments: {len(parsed_data.experiments)}",
            f"  Biospecimens: {len(parsed_data.biospecimens)}",
            f"  Biosamples: {len(parsed_data.biosamples)}",
            f"  Run biosamples: {len(parsed_data.run_biosamples)}",
            f"  Genes: {len(parsed_data.genes)}",
            f"  Proteins: {len(parsed_data.proteins)}",
            f"  Metabolites: {len(parsed_data.metabolites)}",
            f"  Stats datasets: {len(parsed_data.stats_results)}",
        ]

    def _fmt(issue):
        parts = [p for p in [issue.sheet, issue.field,
                              f"row {issue.row}" if issue.row is not None else None]
                 if p]
        loc = " / ".join(parts)
        return f"  [{loc}] {issue.message}" if loc else f"  {issue.message}"

    return {
        "error_count":    len(errors),
        "warning_count":  len(warnings),
        "coverage_lines": cov_lines,
        "error_samples":   [_fmt(i) for i in errors[:5]],
        "warning_samples": [_fmt(i) for i in warnings[:5]],
        "issue_breakdown": issue_breakdown,
        "data_counts": data_counts,
    }


# ── Routes ────────────────────────────────────────────────────────────────────

@router.post("/api/validate")
def validate_sheet_tables(submission: PounceValidationRequest):
    """Validate decoded v2/v3 sheet tables without creating an upload session."""
    project_id = submission.project.source_id
    exp_ids = [workbook.source_id for workbook in submission.experiments]
    stats_ids = [workbook.source_id for workbook in submission.stats_results]
    loader = _sheet_workbook_loader(
        [submission.project, *submission.experiments, *submission.stats_results]
    )
    try:
        parsed_data, structural_issues = PounceParser().parse_all(
            project_id, exp_ids, stats_ids, workbook_loader=loader
        )
    except (UnsupportedPounceFormat, KeyError, ValueError, TypeError) as exc:
        raise HTTPException(status_code=422, detail=str(exc)) from exc

    issues, _ = _validate_parsed_submission(
        parsed_data, structural_issues, project_id, exp_ids, stats_ids
    )
    coverage = _compute_coverage(parsed_data)
    return {
        "schema_version": 1,
        "files": [
            {"role": role, "source_id": source_id,
             "template_version": parsed_data.workbook_versions[source_id]}
            for role, source_ids in (
                ("project", [project_id]),
                ("experiment", exp_ids),
                ("stats_results", stats_ids),
            ) for source_id in source_ids
        ],
        "error_count": sum(issue.severity == "error" for issue in issues),
        "warning_count": sum(issue.severity == "warning" for issue in issues),
        "issues": [
            {**{key: value for key, value in asdict(issue).items() if key != "source_file"},
             "source_id": issue.source_file}
            for issue in issues
        ],
        "coverage": [
            {**asdict(item), "mapped_pct": item.mapped_pct}
            for item in coverage
        ],
    }

@router.get("/validate", response_class=HTMLResponse)
async def submission_home(request: Request):
    return _templates.TemplateResponse(request, "pounce_submission_home.html", {"request": request})


@router.get("/validate/new", response_class=HTMLResponse)
async def upload_form(request: Request):
    return _templates.TemplateResponse(request, "pounce_upload.html", {"request": request})


@router.get("/validate/existing", response_class=HTMLResponse)
async def existing_project_form(request: Request, q: str = ""):
    projects = []
    error = None
    try:
        projects = _search_existing_projects(q.strip())
    except Exception as e:
        error = str(e)

    return _templates.TemplateResponse(request, "pounce_existing_project.html", {
        "request": request,
        "projects": projects,
        "query": q,
        "error": error,
        "mysql_db_name": _MYSQL_DB_NAME,
        "public_project_base_url": _public_project_base_url,
    })


@router.post("/validate/existing/load", response_class=HTMLResponse)
async def load_existing_project(request: Request, project_id: str = Form(...)):
    _cleanup_old_sessions()
    try:
        session_id = _create_existing_project_session(project_id)
        return await edit_form(request, session_id)
    except Exception as e:
        return _templates.TemplateResponse(request, "pounce_existing_project.html", {
            "request": request,
            "projects": [],
            "query": "",
            "error": str(e),
            "mysql_db_name": _MYSQL_DB_NAME,
            "public_project_base_url": _public_project_base_url,
        })


@router.post("/validate", response_class=HTMLResponse)
async def run_validation(
    request: Request,
    project_file:     UploadFile       = File(...),
    experiment_files: List[UploadFile] = File(...),
    stats_files:      List[UploadFile] = File(...),
):
    _cleanup_old_sessions()

    parse_error = None
    parsed_data = None
    grouped     = []
    all_issues  = []
    session_id  = None

    tmp = tempfile.mkdtemp()

    try:
        project_path = _save_upload(project_file, tmp)
        exp_paths    = [_save_upload(f, tmp, prefix=f"exp{i}_")   for i, f in enumerate(experiment_files)]
        stats_paths  = [_save_upload(f, tmp, prefix=f"stats{i}_") for i, f in enumerate(stats_files)]

        adapter = PounceInputAdapter(
            project_file=project_path,
            experiment_files=exp_paths,
            stats_results_files=stats_paths,
            validators_config=_VALIDATORS_CONFIG,
        )
        parsed_data       = adapter.get_validation_data()
        structural_issues = adapter.get_structural_issues()

        all_issues, grouped = _validate_parsed_submission(
            parsed_data, structural_issues, project_path, exp_paths, stats_paths
        )

        # Create session so files survive until the submitter clicks Send.
        session_id, _ = _build_session(
            project_path,
            exp_paths,
            stats_paths,
            parsed_data=parsed_data,
            submission_label="New project submission",
            file_change_summary=_summarize_new_submission(project_path, exp_paths, stats_paths),
        )

    except Exception as e:
        parse_error = str(e)
        shutil.rmtree(tmp, ignore_errors=True)

    errors   = [i for i in all_issues if i.severity == "error"]
    warnings = [i for i in all_issues if i.severity == "warning"]
    coverage = _compute_coverage(parsed_data)

    has_metabolites         = bool(parsed_data and parsed_data.metabolites)
    has_genes               = bool(parsed_data and parsed_data.genes)
    has_proteins            = bool(parsed_data and parsed_data.proteins)
    metab_resolver_missing  = has_metabolites and "Metabolite" not in _resolver_map
    gene_resolver_missing   = has_genes       and "Gene"       not in _resolver_map
    protein_resolver_missing = has_proteins   and "Protein"    not in _resolver_map

    # Fill in session summary now that coverage and resolver state are available.
    if session_id and session_id in _sessions:
        _sessions[session_id]["summary"] = _compute_summary(
            all_issues,
            coverage,
            grouped=grouped,
            parsed_data=parsed_data,
            metab_resolver_missing=metab_resolver_missing,
            gene_resolver_missing=gene_resolver_missing,
            protein_resolver_missing=protein_resolver_missing,
        )

    session_filenames = (
        _sessions[session_id]["file_basenames"] if session_id and session_id in _sessions else []
    )

    return _templates.TemplateResponse(request, "pounce_results.html", {
        "request":                request,
        "parsed_data":            parsed_data,
        "grouped":                grouped,
        "errors":                 errors,
        "warnings":               warnings,
        "parse_error":            parse_error,
        "project_filename":       project_file.filename,
        "coverage":               coverage,
        "metab_resolver_missing": metab_resolver_missing,
        "gene_resolver_missing":  gene_resolver_missing,
        "protein_resolver_missing": protein_resolver_missing,
        "session_id":             session_id,
        "session_filenames":      session_filenames,
        "file_versions":          _file_version_rows(parsed_data, project_path, exp_paths, stats_paths) if parsed_data else [],
        "smtp_configured":        bool(_smtp_config),
    })


@router.get("/download/{session_id}/{filename}")
async def download_file(session_id: str, filename: str):
    """Serve an uploaded file back to the submitter for editing."""
    _cleanup_old_sessions()
    session = _sessions.get(session_id)
    if not session:
        return HTMLResponse("Session expired — please re-upload your files.", status_code=404)
    # Guard against path traversal: only serve filenames actually in the session.
    if filename not in session["file_basenames"]:
        return HTMLResponse("File not found.", status_code=404)
    for fp in session["file_paths"]:
        if os.path.basename(fp) == filename:
            return FileResponse(
                fp,
                filename=filename,
                media_type="application/vnd.openxmlformats-officedocument.spreadsheetml.sheet",
            )
    return HTMLResponse("File not found.", status_code=404)


@router.get("/edit/{session_id}", response_class=HTMLResponse)
async def edit_form(request: Request, session_id: str):
    """Show the upload form pre-populated with the current session's filenames."""
    _cleanup_old_sessions()
    session = _sessions.get(session_id)
    if not session:
        return _templates.TemplateResponse(request, "pounce_submit_confirm.html", {
            "request": request,
            "error": "Session expired — please re-upload your files.",
            "name": "", "email": "", "project_name": "",
        })
    return _templates.TemplateResponse(request, "pounce_edit.html", {
        "request":          request,
        "session_id":       session_id,
        "mode_title":       "Edit Project Files",
        "mode_subtitle":    "Download the current files, replace any of them, and append new experiment and StatsResults workbook pairs if needed.",
        "project_basename": os.path.basename(session["project_path"]),
        "exp_basenames":    [os.path.basename(p) for p in session["exp_paths"]],
        "stats_basenames":  [os.path.basename(p) for p in session["stats_paths"]],
    })


@router.post("/edit/{session_id}", response_class=HTMLResponse)
async def revalidate(
    request:          Request,
    session_id:       str,
    project_file:     UploadFile       = File(...),
    experiment_files: List[UploadFile] = File(default=[]),
    stats_files:      List[UploadFile] = File(default=[]),
    new_experiment_files: List[UploadFile] = File(default=[]),
    new_stats_files:      List[UploadFile] = File(default=[]),
):
    """Re-run validation, replacing only the slots where a new file was uploaded."""
    _cleanup_old_sessions()
    session = _sessions.get(session_id)
    if not session:
        return _templates.TemplateResponse(request, "pounce_submit_confirm.html", {
            "request": request,
            "error": "Session expired — please re-upload your files.",
            "name": "", "email": "", "project_name": "",
        })

    parse_error  = None
    parsed_data  = None
    grouped      = []
    all_issues   = []
    new_sid      = None
    project_path = ""
    tmp          = tempfile.mkdtemp()

    try:
        # Project — replace if a new file was provided, else copy existing.
        if project_file.filename:
            project_path = _save_upload(project_file, tmp)
        else:
            project_path = shutil.copy2(session["project_path"], tmp)

        # Experiments — zip new uploads with existing paths.
        exp_paths = []
        for i, (new_f, existing_p) in enumerate(zip(experiment_files, session["exp_paths"])):
            if new_f.filename:
                exp_paths.append(_save_upload(new_f, tmp, prefix=f"exp{i}_"))
            else:
                exp_paths.append(shutil.copy2(existing_p, tmp))
        # Keep any existing experiments beyond what was submitted (shouldn't happen, but safe).
        for existing_p in session["exp_paths"][len(exp_paths):]:
            exp_paths.append(shutil.copy2(existing_p, tmp))

        # Stats — same pattern.
        stats_paths = []
        for i, (new_f, existing_p) in enumerate(zip(stats_files, session["stats_paths"])):
            if new_f.filename:
                stats_paths.append(_save_upload(new_f, tmp, prefix=f"stats{i}_"))
            else:
                stats_paths.append(shutil.copy2(existing_p, tmp))
        for existing_p in session["stats_paths"][len(stats_paths):]:
            stats_paths.append(shutil.copy2(existing_p, tmp))

        _append_new_pairs(new_experiment_files, new_stats_files, tmp, exp_paths, stats_paths)

        adapter = PounceInputAdapter(
            project_file=project_path,
            experiment_files=exp_paths,
            stats_results_files=stats_paths,
            validators_config=_VALIDATORS_CONFIG,
        )
        parsed_data       = adapter.get_validation_data()
        structural_issues = adapter.get_structural_issues()

        all_issues, grouped = _validate_parsed_submission(
            parsed_data, structural_issues, project_path, exp_paths, stats_paths
        )

        new_sid, _ = _build_session(
            project_path,
            exp_paths,
            stats_paths,
            parsed_data=parsed_data,
            submission_label=session.get("submission_label") or "Updated submission",
            file_change_summary=_summarize_edit_submission(
                session,
                project_file,
                experiment_files,
                stats_files,
                new_experiment_files,
                new_stats_files,
            ),
            existing_project_id=session.get("existing_project_id"),
        )

    except Exception as e:
        parse_error = str(e)
        shutil.rmtree(tmp, ignore_errors=True)

    errors   = [i for i in all_issues if i.severity == "error"]
    warnings = [i for i in all_issues if i.severity == "warning"]
    coverage = _compute_coverage(parsed_data)

    has_metabolites        = bool(parsed_data and parsed_data.metabolites)
    has_genes              = bool(parsed_data and parsed_data.genes)
    has_proteins           = bool(parsed_data and parsed_data.proteins)
    metab_resolver_missing = has_metabolites and "Metabolite" not in _resolver_map
    gene_resolver_missing  = has_genes       and "Gene"       not in _resolver_map
    protein_resolver_missing = has_proteins  and "Protein"    not in _resolver_map

    if new_sid and new_sid in _sessions:
        _sessions[new_sid]["summary"] = _compute_summary(
            all_issues,
            coverage,
            grouped=grouped,
            parsed_data=parsed_data,
            metab_resolver_missing=metab_resolver_missing,
            gene_resolver_missing=gene_resolver_missing,
            protein_resolver_missing=protein_resolver_missing,
        )

    session_filenames = (
        _sessions[new_sid]["file_basenames"] if new_sid and new_sid in _sessions else []
    )

    return _templates.TemplateResponse(request, "pounce_results.html", {
        "request":                request,
        "parsed_data":            parsed_data,
        "grouped":                grouped,
        "errors":                 errors,
        "warnings":               warnings,
        "parse_error":            parse_error,
        "project_filename":       os.path.basename(project_path),
        "coverage":               coverage,
        "metab_resolver_missing": metab_resolver_missing,
        "gene_resolver_missing":  gene_resolver_missing,
        "protein_resolver_missing": protein_resolver_missing,
        "session_id":             new_sid,
        "session_filenames":      session_filenames,
        "file_versions":          _file_version_rows(parsed_data, project_path, exp_paths, stats_paths) if parsed_data else [],
        "smtp_configured":        bool(_smtp_config),
    })


@router.post("/submit/{session_id}", response_class=HTMLResponse)
async def submit_pounce(
    request: Request,
    session_id: str,
    name: str = Form(...),
    email: str = Form(...),
):
    _cleanup_old_sessions()

    session = _sessions.get(session_id)
    if not session:
        return _templates.TemplateResponse(request, "pounce_submit_confirm.html", {
            "request": request,
            "error": "Session expired — please re-upload your files and try again.",
            "name": name,
            "email": email,
            "project_name": "",
        })

    if not _smtp_config:
        return _templates.TemplateResponse(request, "pounce_submit_confirm.html", {
            "request": request,
            "error": "Email delivery is not configured on this server.",
            "name": name,
            "email": email,
            "project_name": session["project_name"],
        })

    # Build email message.
    summary         = session.get("summary") or {}
    error_count     = summary.get("error_count", 0)
    warning_count   = summary.get("warning_count", 0)
    coverage_lines  = summary.get("coverage_lines", [])
    error_samples   = summary.get("error_samples", [])
    warning_samples = summary.get("warning_samples", [])
    issue_breakdown = summary.get("issue_breakdown", [])
    data_counts     = summary.get("data_counts", [])
    file_change_summary = session.get("file_change_summary") or {}
    submission_label = session.get("submission_label") or "POUNCE submission"

    body_lines = [
        "POUNCE Submission",
        "",
        f"Submission type: {submission_label}",
        f"Submitter: {name} <{email}>",
        f"Project:   {session['project_name']}",
    ]
    if session.get("existing_project_id"):
        body_lines.append(f"Existing project ID: {session['existing_project_id']}")

    added_files = file_change_summary.get("added") or []
    replaced_files = file_change_summary.get("replaced") or []
    unchanged_files = file_change_summary.get("unchanged") or []
    if added_files or replaced_files or unchanged_files:
        body_lines.extend([
            "",
            "=== File Change Summary ===",
        ])
        if added_files:
            body_lines.append("  Added:")
            body_lines.extend([f"    - {item}" for item in added_files])
        if replaced_files:
            body_lines.append("  Replaced:")
            body_lines.extend([f"    - {item}" for item in replaced_files])
        if unchanged_files:
            body_lines.append("  Unchanged:")
            body_lines.extend([f"    - {item}" for item in unchanged_files])

    body_lines.extend([
        "",
        "=== Attached Files ===",
    ])
    body_lines.extend([f"  - {filename}" for filename in session.get("file_basenames", [])])
    if data_counts:
        body_lines.append("")
        body_lines.append("=== Parsed Data Snapshot ===")
        body_lines.extend(data_counts)

    body_lines.extend([
        "",
        "=== Validation Summary ===",
        f"  Errors:   {error_count}",
    ])
    if error_samples:
        body_lines.extend(error_samples)
        if error_count > 5:
            body_lines.append(f"  … {error_count - 5} more error(s)")
    body_lines.append(f"  Warnings: {warning_count}")
    if warning_samples:
        body_lines.extend(warning_samples)
        if warning_count > 5:
            body_lines.append(f"  … {warning_count - 5} more warning(s)")
    if issue_breakdown:
        body_lines.append("")
        body_lines.append("=== Issue Breakdown ===")
        body_lines.extend(issue_breakdown)
    if coverage_lines:
        body_lines.append("")
        body_lines.append("=== Mapping Coverage ===")
        body_lines.extend(coverage_lines)

    msg = MIMEMultipart()
    msg["Subject"]  = f"POUNCE Submission: {session['project_name']}"
    msg["From"]     = _smtp_config["from_address"]
    msg["To"]       = _smtp_config["to_address"]
    msg["Cc"]       = email
    msg["Reply-To"] = email
    msg.attach(MIMEText("\n".join(body_lines), "plain"))

    # Attach all xlsx files.
    for file_path in session["file_paths"]:
        if os.path.exists(file_path):
            with open(file_path, "rb") as f:
                part = MIMEBase("application", "octet-stream")
                part.set_payload(f.read())
            encoders.encode_base64(part)
            part.add_header(
                "Content-Disposition",
                "attachment",
                filename=os.path.basename(file_path),
            )
            msg.attach(part)

    # Send via SMTP.
    try:
        host     = _smtp_config.get("host", "localhost")
        port     = _smtp_config.get("port", 587)
        use_tls  = _smtp_config.get("use_tls", True)
        user     = _smtp_config.get("user")
        password = _smtp_config.get("password")
        timeout  = _smtp_config.get("timeout_seconds", 10)

        smtp = smtplib.SMTP(host, port, timeout=timeout)
        if use_tls:
            smtp.starttls()
        if user and password:
            smtp.login(user, password)
        smtp.sendmail(_smtp_config["from_address"], [_smtp_config["to_address"], email], msg.as_string())
        smtp.quit()
    except Exception as e:
        return _templates.TemplateResponse(request, "pounce_submit_confirm.html", {
            "request":      request,
            "error":        f"Failed to send email: {e}",
            "name":         name,
            "email":        email,
            "project_name": session["project_name"],
        })

    # Clean up session immediately after successful send.
    shutil.rmtree(session["dir"], ignore_errors=True)
    del _sessions[session_id]

    return _templates.TemplateResponse(request, "pounce_submit_confirm.html", {
        "request":      request,
        "error":        None,
        "name":         name,
        "email":        email,
        "project_name": session["project_name"],
    })
