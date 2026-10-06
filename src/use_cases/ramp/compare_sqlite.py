"""Read-only, offline comparisons of RaMP SQLite application databases."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import html
import json
import math
import os
from pathlib import Path
import sqlite3
import tempfile


SECTIONS = (
    "Metabolites", "Core coverage", "Source identifiers by namespace", "Source attribution",
    "Pathways by provider", "Chemistry by provider", "Chemistry completeness",
    "Metabolite group sizes", "Relationship integrity", "Table rows",
)
CHEM_FIELDS = ("iso_smiles", "inchi_key", "inchi", "mw", "monoisotop_mass", "common_name", "mol_formula")
GROUP_BINS = ((0, '0'), (1, '1'), (5, '2–5'), (10, '6–10'), (20, '11–20'),
              (50, '21–50'), (100, '51–100'), (500, '101–500'), (1000, '501–1000'), (float('inf'), '1001+'))
GROUP_METRICS = [f'{label} source IDs: groups' for _, label in GROUP_BINS] + ['Largest metabolite group (source IDs)']
MULTI_RAMP_METRICS = (
    "Source IDs linked to multiple metabolite RaMP IDs",
    "Source IDs linked to multiple gene/protein RaMP IDs",
    "Source IDs linked to both metabolites and genes/proteins",
)
METABOLITE_SOURCE_NAMES = {
    'chebi': 'ChEBI', 'hmdb': 'HMDB', 'hmdb_kegg': 'HMDB (KEGG)',
    'lipidmaps': 'LipidMaps', 'pfocr': 'PFOCR', 'pubchem': 'PubChem',
    'reactome': 'Reactome', 'refmet': 'RefMet', 'rhea': 'Rhea',
    'wiki': 'WikiPathways', 'wikipathways_kegg': 'WikiPathways (KEGG)',
}
LINKS = {
    "source": (("rampId", "analyte", "rampId"),),
    "analytesynonym": (("rampId", "analyte", "rampId"),),
    "chem_props": (("ramp_id", "analyte", "rampId"),),
    "metabolite_class": (("ramp_id", "analyte", "rampId"),),
    "analytehaspathway": (("rampId", "analyte", "rampId"), ("pathwayRampId", "pathway", "pathwayRampId")),
    "analytehasontology": (("rampCompoundId", "analyte", "rampId"), ("rampOntologyId", "ontology", "rampOntologyId")),
    "catalyzed": (("rampCompoundId", "analyte", "rampId"), ("rampGeneId", "analyte", "rampId")),
    "reaction2met": (("ramp_cmpd_id", "analyte", "rampId"), ("ramp_rxn_id", "reaction", "ramp_rxn_id")),
    "reaction2protein": (("ramp_gene_id", "analyte", "rampId"), ("ramp_rxn_id", "reaction", "ramp_rxn_id")),
    "reaction_ec_class": (("ramp_rxn_id", "reaction", "ramp_rxn_id"),),
}


def quoted(name):
    return '"' + name.replace('"', '""') + '"'


def metric(value=None, *, state="available", denominator=None):
    return {"value": value, "state": state, "denominator": denominator}


def normalize_identifier(identifier):
    """Normalize only namespace case; accessions and identifier families stay intact."""
    prefix, sep, accession = identifier.strip().partition(":")
    return prefix.casefold() + sep + accession if sep else identifier.strip()


def metabolite_group_metrics(db):
    counts = {label: 0 for _, label in GROUP_BINS}
    largest = 0
    query = "SELECT count(DISTINCT s.sourceId) FROM analyte a LEFT JOIN source s ON s.rampId=a.rampId WHERE a.type='compound' GROUP BY a.rampId"
    for (size,) in db.execute(query):
        largest = max(largest, size)
        label = next(label for upper, label in GROUP_BINS if size <= upper)
        counts[label] += 1
    return {**{f'{label} source IDs: groups': metric(count) for label, count in counts.items()},
            'Largest metabolite group (source IDs)': metric(largest)}


def profile_database(label, path, progress=print):
    path = Path(path).resolve(strict=True)
    progress(f"Profiling {label}: {path.name}")
    db = sqlite3.connect(path.as_uri() + "?mode=ro", uri=True)
    db.row_factory = sqlite3.Row
    report = {"label": label, "path": str(path), "size_bytes": path.stat().st_size,
              "sections": {s: {} for s in SECTIONS}, "notes": [], "source_versions": [],
              "schema": {}, "manifest": {}, "versions": []}
    identity_sets = defaultdict(set)
    try:
        db.execute("PRAGMA query_only=ON")
        db.execute("BEGIN")
        tables = [r[0] for r in db.execute("SELECT name FROM sqlite_master WHERE type='table' AND name NOT LIKE 'sqlite_%'")]
        for table in tables:
            report["schema"][table] = [dict(r) for r in db.execute(f"PRAGMA table_info({quoted(table)})")]

        def has(table, *columns):
            return table in report["schema"] and set(columns) <= {c["name"] for c in report["schema"][table]}

        def scalar(section, name, sql, required):
            report["sections"][section][name] = (
                metric(db.execute(sql).fetchone()[0]) if all(has(t, *cols) for t, cols in required.items())
                else metric(state="unavailable")
            )

        if has("ramp_export_metadata", "key", "value"):
            row = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
            if row:
                report["manifest"] = json.loads(row[0])
        manifest = report["manifest"]
        if manifest.get("status") == "post_processing_pending":
            report["notes"].append("Base-table export: post-processing is pending. Empty derived tables are not data losses.")
        if manifest.get("gene_policy"):
            prefix = "Gene/protein grouping: " if manifest.get('gene_resolution') else "Gene/protein grouping is provisional: "
            report["notes"].append(prefix + manifest["gene_policy"])
        if manifest.get('gene_resolution'):
            report['notes'].append('Gene/protein resolution counts: ' + json.dumps(manifest['gene_resolution']['counts'], sort_keys=True))
        if manifest.get("omitted_tables"):
            report["notes"].append("Intentionally omitted: " + ", ".join(manifest["omitted_tables"]))
        if has("db_version", "ramp_version"):
            columns = [c for c in ("ramp_version", "load_timestamp", "version_notes") if has("db_version", c)]
            report["versions"] = [dict(r) for r in db.execute("SELECT " + ",".join(map(quoted, columns)) + " FROM db_version")]
        if has("version_info", "data_source_id"):
            query = "SELECT * FROM version_info" + (" WHERE status='current'" if has("version_info", "status") else "")
            report["source_versions"] = [dict(r) for r in db.execute(query)]

        for table in tables:
            state = "pending" if table in manifest.get("pending_tables", []) else "available"
            report["sections"]["Table rows"][table] = metric(db.execute(f"SELECT count(*) FROM {quoted(table)}").fetchone()[0], state=state)

        core = "Core coverage"
        for kind, name in (("compound", "Metabolites"), ("gene", "Gene/protein analytes")):
            scalar(core, name, f"SELECT count(DISTINCT rampId) FROM analyte WHERE type='{kind}'", {"analyte": ["rampId", "type"]})
        report["sections"]["Metabolites"]["Total metabolites"] = report["sections"][core]["Metabolites"].copy()
        for table, column, name in (("pathway", "pathwayRampId", "Pathways"), ("reaction", "ramp_rxn_id", "Reactions"),
                                    ("chem_props", "ramp_id", "Metabolites with chemistry"),
                                    ("metabolite_class", "ramp_id", "Metabolites with classification"),
                                    ("analytehasontology", "rampCompoundId", "Metabolites with ontology"),
                                    ("analytesynonym", "rampId", "Analytes with synonyms")):
            scalar(core, name, f"SELECT count(DISTINCT {quoted(column)}) FROM {quoted(table)}", {table: [column]})
        for kind, name in (("compound", "Metabolites mapped to pathways"), ("gene", "Gene/protein analytes mapped to pathways")):
            scalar(core, name, f"SELECT count(DISTINCT e.rampId) FROM analytehaspathway e JOIN analyte a ON a.rampId=e.rampId WHERE a.type='{kind}'",
                   {"analytehaspathway": ["rampId"], "analyte": ["rampId", "type"]})

        if has("source", "sourceId"):
            for row in db.execute("SELECT DISTINCT sourceId FROM source WHERE sourceId IS NOT NULL AND trim(sourceId) != ''"):
                identifier = normalize_identifier(row[0])
                namespace = identifier.split(":", 1)[0] if ":" in identifier else "(no namespace)"
                identity_sets[namespace].add(identifier)
            for namespace, identifiers in sorted(identity_sets.items()):
                report["sections"]["Source identifiers by namespace"][namespace] = metric(len(identifiers))
            report["sections"][core]["Unique source identifiers"] = metric(sum(map(len, identity_sets.values())))
        for kind, label in (("compound", "metabolite"), ("gene", "gene/protein")):
            scalar(core, f"Source IDs linked to multiple {label} RaMP IDs",
                   "SELECT count(*) FROM (SELECT sourceId FROM source "
                   f"WHERE geneOrCompound='{kind}' AND sourceId IS NOT NULL "
                   "GROUP BY sourceId HAVING count(DISTINCT rampId)>1)",
                   {"source": ["sourceId", "rampId", "geneOrCompound"]})
        scalar(core, "Source IDs linked to both metabolites and genes/proteins",
               "SELECT count(*) FROM (SELECT sourceId FROM source "
               "WHERE sourceId IS NOT NULL AND rampId IS NOT NULL "
               "AND geneOrCompound IN ('compound','gene') "
               "GROUP BY sourceId HAVING count(DISTINCT geneOrCompound)=2)",
               {"source": ["sourceId", "rampId", "geneOrCompound"]})

        if has("source", "dataSource", "sourceId", "rampId", "geneOrCompound"):
            for r in db.execute("SELECT dataSource, geneOrCompound, count(*), count(DISTINCT sourceId), count(DISTINCT rampId) FROM source GROUP BY dataSource,geneOrCompound"):
                for suffix, value in zip(("rows", "source IDs", "analytes"), r[2:]):
                    report["sections"]["Source attribution"][f"{r[0] or '(missing)'} / {r[1] or '(missing)'} / {suffix}"] = metric(value)
        if has("source", "dataSource", "rampId", "geneOrCompound"):
            # Count the distinct union before display-name mapping: summing raw
            # provider groups would double-count overlapping metabolite IDs.
            provider_key = "coalesce(nullif(lower(trim(dataSource)),''),'(missing provider)')"
            for raw_name, count in db.execute(
                    f"SELECT {provider_key}, count(DISTINCT rampId) FROM source "
                    f"WHERE geneOrCompound='compound' GROUP BY {provider_key}"):
                display_name = METABOLITE_SOURCE_NAMES.get(raw_name, raw_name)
                report["sections"]["Metabolites"][f"With data from {display_name}"] = metric(count)
        if has("pathway", "type", "pathwayRampId"):
            for r in db.execute("SELECT type,count(DISTINCT pathwayRampId) FROM pathway GROUP BY type"):
                report["sections"]["Pathways by provider"][r[0] or "(missing)"] = metric(r[1])

        if has("chem_props", "chem_data_source", "chem_source_id", "ramp_id"):
            fields = [f for f in CHEM_FIELDS if has("chem_props", f)]
            expressions = [f"count(DISTINCT CASE WHEN {quoted(f)} IS NOT NULL AND trim(CAST({quoted(f)} AS TEXT)) != '' THEN chem_source_id END)" for f in fields]
            query = "SELECT chem_data_source,count(*),count(DISTINCT chem_source_id),count(DISTINCT ramp_id)"
            query += ("," + ",".join(expressions) if expressions else "") + " FROM chem_props GROUP BY chem_data_source"
            for row in db.execute(query):
                src, rows, identities, analytes = row[:4]
                src = src or "(missing)"
                for suffix, value in (("rows", rows), ("source IDs", identities), ("metabolites", analytes)):
                    report["sections"]["Chemistry by provider"][f"{src} / {suffix}"] = metric(value)
                for field in CHEM_FIELDS:
                    value = row[4 + fields.index(field)] if field in fields else None
                    report["sections"]["Chemistry completeness"][f"{src} / {field}"] = metric(value, denominator=identities, state="available" if field in fields else "unavailable")
            registered = {str(r["data_source_id"]).casefold() for r in report["source_versions"]}
            providers = {r[0] for r in db.execute("SELECT DISTINCT chem_data_source FROM chem_props WHERE chem_data_source IS NOT NULL")}
            missing = sorted(p for p in providers if p.casefold() not in registered)
            if missing:
                report["notes"].append("Chemistry providers without current version entries: " + ", ".join(missing))

        if has("source", "rampId", "sourceId") and has("analyte", "rampId", "type"):
            report["sections"]["Metabolite group sizes"] = metabolite_group_metrics(db)

        for table, links in LINKS.items():
            for column, target, key in links:
                scalar("Relationship integrity", f"{table}.{column} → {target}: orphan rows",
                       f"SELECT count(*) FROM {quoted(table)} e WHERE NOT EXISTS (SELECT 1 FROM {quoted(target)} t WHERE t.{quoted(key)}=e.{quoted(column)})",
                       {table: [column], target: [key]})
        # Empty dynamic sections are distinguishable from unavailable schemas.
        report["available_sections"] = {
            "Metabolites": has("source", "dataSource", "rampId", "geneOrCompound"),
            "Source identifiers by namespace": has("source", "sourceId"),
            "Source attribution": has("source", "dataSource", "sourceId", "rampId", "geneOrCompound"),
            "Pathways by provider": has("pathway", "type", "pathwayRampId"),
            "Chemistry by provider": has("chem_props", "chem_data_source", "chem_source_id", "ramp_id"),
            "Chemistry completeness": has("chem_props", "chem_data_source", "chem_source_id", "ramp_id"),
            "Metabolite group sizes": has("source", "sourceId", "rampId") and has("analyte", "rampId", "type"),
        }
    finally:
        db.close()
    return report, identity_sets


def compare(databases, progress=print):
    profiles, sets = [], []
    for label, path in databases:
        profile, identities = profile_database(label, path, progress)
        profiles.append(profile)
        sets.append(identities)
    for profile, identifiers in zip(profiles[1:], sets[1:]):
        profile["source_id_changes"] = {}
        if not profiles[0]["available_sections"]["Source identifiers by namespace"] or not profile["available_sections"]["Source identifiers by namespace"]:
            profile["notes"].append("Source-ID overlap unavailable: source identifiers missing from this database or baseline.")
            continue
        for namespace in sorted(set(sets[0]) | set(identifiers)):
            old, new = sets[0].get(namespace, set()), identifiers.get(namespace, set())
            profile["source_id_changes"][namespace] = {"retained": len(old & new), "added": len(new - old), "removed": len(old - new)}
    return {"format_version": 1, "generated_at": datetime.now(timezone.utc).isoformat(),
            "baseline": profiles[0]["label"], "databases": profiles,
            "methodology": [
                "All inputs are read through read-only SQLite transactions; no graph or network access is performed by this command.",
                "Deltas compare each database with the first input. Increases and decreases are descriptive, not pass/fail judgments.",
                "Source-ID overlap trims outer whitespace and case-folds only the namespace before the first colon; accessions remain case-sensitive. No cross-family reconciliation is performed.",
                "RaMP IDs are used only for within-database joins and counts, never to match entities between releases.",
                "Chemistry completeness counts distinct provider source IDs with at least one nonblank value, divided by all distinct source IDs from that provider. It does not imply every assertion is complete or scientifically valid.",
                "Orphan checks count rows whose endpoint is missing, including NULL endpoints, even where old schemas declare no foreign key.",
                "Missing tables/columns are Unavailable; absent provider groups in an available section are zero. Pending tables show their physical row count without a loss delta.",
                "Group sizes count distinct stored source IDs per compound, not chemistry assertions. They measure export membership, not historical RaMP ID continuity.",
                "Recorded versions may be incomplete. No version is inferred from data presence; coverage checks are not R-package compatibility tests.",
            ]}


def get_metric(profile, section, key):
    if key in profile["sections"][section]:
        return profile["sections"][section][key]
    if section == "Chemistry completeness":
        field = key.rsplit(" / ", 1)[-1]
        if field not in {c["name"] for c in profile["schema"].get("chem_props", [])}:
            return metric(state="unavailable")
    if section == "Table rows" and key in profile["manifest"].get("omitted_tables", []):
        return metric(state="omitted")
    if profile["available_sections"].get(section):
        return metric(0, denominator=0 if section == "Chemistry completeness" else None)
    return metric(state="unavailable")


def comparison_cell(value, baseline, show_delta):
    state, n = value["state"], value["value"]
    if state != "available":
        return {"text": state.title() + (f" ({n:,} rows)" if n is not None else ""), "delta": None, "percent_change": None}
    text = f"{n:,}"
    denominator = value["denominator"]
    if denominator is not None:
        text += f" / {denominator:,}" + (f" ({100*n/denominator:.1f}%)" if denominator else " (—)")
    delta, percent = None, None
    if show_delta and baseline["state"] == "available":
        delta = n - baseline["value"]
        percent = 100 * delta / baseline["value"] if baseline["value"] else None
    return {"text": text, "delta": delta, "percent_change": percent}


def prepare_tables(report):
    profiles = report["databases"]
    tables = {}
    table_databases = {}
    excluded = set(report.get('core_coverage_excluded_labels', []))
    unknown = excluded - {p['label'] for p in profiles}
    if unknown:
        raise ValueError(f'Unknown Core coverage exclusions: {sorted(unknown)}')
    for section in SECTIONS:
        selected = [p for p in profiles if section not in ('Core coverage', 'Metabolites') or p['label'] not in excluded]
        if not selected:
            raise ValueError('Core coverage must include at least one database')
        table_databases[section] = [{'label': p['label'], 'compared_with': selected[i-1]['label'] if i else None}
            for i, p in enumerate(selected)]
        rows = []
        keys = {k for p in selected for k in p["sections"][section]}
        ordered = sorted(keys, key=lambda k: (GROUP_METRICS.index(k) if k in GROUP_METRICS else len(GROUP_METRICS), k)) if section == 'Metabolite group sizes' else sorted(keys)
        if section == 'Metabolites':
            ordered = (['Total metabolites'] if 'Total metabolites' in keys else []) + [
                k for k in ordered if k != 'Total metabolites']
        if section == 'Core coverage':
            ordered = [k for k in ordered if k not in MULTI_RAMP_METRICS]
            after = ordered.index('Unique source identifiers') + 1 if 'Unique source identifiers' in ordered else len(ordered)
            ordered[after:after] = [k for k in MULTI_RAMP_METRICS if k in keys]
        for key in ordered:
            rows.append({"metric": key, "cells": [
                {**get_metric(p, section, key), **comparison_cell(get_metric(p, section, key),
                    get_metric(selected[i-1] if i else selected[0], section, key), i > 0)}
                for i, p in enumerate(selected)
            ]})
        tables[section] = rows
    report["comparison_tables"] = tables
    report['table_databases'] = table_databases
    report['methodology'] = [
        'All table deltas compare each column with the preceding displayed database.'
        if note.startswith('Deltas compare each database') or note.startswith('Core coverage deltas compare') else note
        for note in report['methodology']]
    version_maps = []
    for profile in profiles:
        versions = defaultdict(set)
        for row in profile['source_versions']:
            source = str(row.get('data_source_id') or '(missing)').casefold()
            versions[source].add(str(row.get('data_source_version') or 'Not recorded'))
        version_maps.append(versions)
    report['source_version_table'] = [
        {'source': source, 'versions': ['; '.join(sorted(versions.get(source, {'Not recorded'}))) for versions in version_maps]}
        for source in sorted({source for versions in version_maps for source in versions})]
    return report


def delta_shading(cell):
    """Continuous directional tint; undefined percentages have no magnitude."""
    percent = cell.get('percent_change')
    if cell.get('delta') is None or percent is None or percent == 0 or not math.isfinite(percent):
        return ''
    strength = min(abs(percent), 100) / 100
    endpoint = (173, 222, 188) if percent > 0 else (246, 185, 185)
    rgb = ','.join(str(round(255 + (channel - 255) * strength)) for channel in endpoint)
    return f' style="background-color:rgb({rgb})"'


def render_html(report):
    escape = lambda value: html.escape(str(value))
    profiles = report["databases"]
    parts = ["<!doctype html><html lang='en'><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'>",
             "<title>RaMP database comparison</title><style>body{font:16px/1.5 system-ui,sans-serif;margin:0;color:#182635;background:#f5f7fa}main{max-width:1400px;padding:24px;margin:auto}h1{margin-bottom:4px}h2{margin-top:32px}p{max-width:100ch}.scroll{max-width:100%;overflow-x:auto;background:white;border:1px solid #d2dbe3;border-radius:6px}table{border-collapse:collapse;width:100%;font-size:14px}th,td{padding:9px 12px;border-bottom:1px solid #dde4eb;text-align:right;vertical-align:top;min-width:145px}th:first-child,td:first-child{text-align:left;min-width:240px}th{background:#eaf0f6}tbody tr:nth-child(even){background:#f8fafc}small{display:block;color:#35475a}.delta-legend{margin:16px 0}.delta-legend p{margin-top:6px;font-size:14px}.delta-scale{height:14px;max-width:360px;border:1px solid #d2dbe3;border-radius:4px;background:linear-gradient(to right,rgb(246,185,185),#fff,rgb(173,222,188))}code{overflow-wrap:anywhere}article{background:white;border-left:4px solid #52738c;padding:12px 18px;margin:12px 0}details{margin:12px 0}pre{white-space:pre-wrap;overflow-wrap:anywhere;font-size:12px}nav a{display:inline-block;margin:4px 14px 4px 0;color:#23527c}@media(max-width:600px){main{padding:14px}h1{font-size:26px}}@media print{.scroll{overflow:visible}details{display:block}body{background:white}}</style><main>",
             "<h1>RaMP database comparison</h1>",
             f"<p>All table deltas compare each column with the preceding displayed database. Generated {escape(report['generated_at'])}.</p>",
             "<p>Changes describe coverage, not correctness. RaMP IDs are not matched across releases.</p>",
             '<div class="delta-legend"><div class="delta-scale" aria-hidden="true"></div>'
             '<p>Red: decrease · Green: increase. Stronger shading means a larger percentage change, capped at 100%. '
             'Colors indicate direction, not quality. Undefined percentages are unshaded.</p></div>']
    for p in profiles:
        versions = ", ".join(str(v.get("ramp_version")) for v in p["versions"]) or "Unknown"
        parts.append(f"<article><strong>{escape(p['label'])}</strong> · embedded version {escape(versions)}")
        if p["manifest"].get("stage_id"):
            parts.append(f"<p>Stage: <code>{escape(p['manifest']['stage_id'])}</code></p>")
        parts.extend(f"<p>{escape(note)}</p>" for note in p["notes"])
        parts.append(f"<details><summary>Input file and metadata</summary><code>{escape(p['path'])}</code><p>{p['size_bytes']:,} bytes</p><pre>{escape(json.dumps(p['versions'],indent=2))}</pre></details></article>")
    parts.append("<nav>" + " ".join(f"<a href='#s{i}'>{escape(s)}</a>" for i, s in enumerate(SECTIONS)) + "</nav>")
    for i, section in enumerate(SECTIONS):
        header = "<thead><tr><th>Metric</th>" + "".join(f"<th>{escape(p['label'])}</th>" for p in report['table_databases'][section]) + "</tr></thead>"
        parts.append(f"<h2 id='s{i}'>{escape(section)}</h2>")
        if section == 'Metabolites':
            parts.append('<p>Source rows count distinct metabolite RaMP IDs recorded for each provider in the source table. '
                         'A metabolite can appear under several sources, so source counts do not add up to the total. '
                         'Chemistry recorded only in chem_props is shown separately below and does not count here. '
                         'Changes compare each column with the previous displayed build.</p>')
        if section == 'Core coverage':
            parts.append('<p>Changes compare each column with the preceding displayed database.</p>')
            parts.append('<p>Unique source identifiers ignore namespace capitalization and outer whitespace. '
                         'Each mapping count measures exact stored source identifiers linked to more than one distinct RaMP ID within that entity type. '
                         'Genes and proteins are counted together. Repeated provider rows for the same RaMP ID do not increase the count. '
                         'The cross-type count measures exact stored source identifiers linked to at least one metabolite RaMP ID and at least one gene/protein RaMP ID; '
                         'each identifier is counted once. These counts can overlap.</p>')
        if section == "Chemistry completeness":
            parts.append("<p>Source IDs with a value / all source IDs from that provider. Percentages measure coverage; deltas compare populated-ID counts.</p>")
        if section == 'Metabolite group sizes':
            parts.append('<p>Group size is the number of distinct source identifiers assigned to a compound RaMP ID. '
                         'The ranges count groups; the final row reports the largest size.</p>')
        parts.append("<div class='scroll'><table>" + header + "<tbody>")
        for row in report["comparison_tables"][section]:
            label = escape(row['metric'])
            if section == 'Metabolites' and row['metric'] == 'Total metabolites':
                label = '<strong>' + label + '</strong>'
            parts.append(f"<tr><td>{label}</td>")
            for cell in row["cells"]:
                parts.append("<td" + delta_shading(cell) + ">" + escape(cell["text"]))
                if cell["delta"] is not None:
                    pct = f"{cell['percent_change']:+.1f}%" if cell["percent_change"] is not None else "—"
                    parts.append(f"<small>Δ {cell['delta']:+,} ({pct})</small>")
                parts.append("</td>")
            parts.append("</tr>")
        parts.append("</tbody></table></div>")
    parts.append("<h2>Recorded source versions</h2>")
    parts.append('<p>Not recorded means version metadata is missing; it does not establish that the source is absent.</p>')
    parts.append("<div class='scroll'><table><thead><tr><th>Source</th>" + ''.join(f"<th>{escape(p['label'])}</th>" for p in profiles) + '</tr></thead><tbody>')
    for row in report['source_version_table']:
        parts.append('<tr><td>' + escape(row['source']) + '</td>' + ''.join('<td>' + escape(v) + '</td>' for v in row['versions']) + '</tr>')
    parts.append('</tbody></table></div>')
    parts.append("<h2>Methodology</h2><ul>" + "".join(f"<li>{escape(n)}</li>" for n in report["methodology"]) + "</ul>")
    parts.append("<details><summary>Schema and export manifests</summary><pre>" + escape(json.dumps({p["label"]: {"schema": p["schema"], "manifest": p["manifest"]} for p in profiles}, indent=2)) + "</pre></details></main></html>")
    return "\n".join(parts)


def write_report(report, output):
    output = Path(output).resolve()
    if output.suffix.lower() != ".html":
        raise ValueError("Report output must end in .html")
    paths = (output, output.with_suffix(".json"))
    inputs = [Path(p["path"]).resolve() for p in report["databases"]]
    for path in paths:
        if path in inputs or any(path.exists() and path.samefile(p) for p in inputs):
            raise ValueError("Report output must not overwrite an input database")
    output.parent.mkdir(parents=True, exist_ok=True)
    for path, content in zip(paths, (render_html(report), json.dumps(report, indent=2))):
        fd, name = tempfile.mkstemp(dir=output.parent, suffix=".tmp")
        try:
            with os.fdopen(fd, "w") as handle:
                handle.write(content)
            os.replace(name, path)
        finally:
            Path(name).unlink(missing_ok=True)
    return paths


def main(argv=None):
    parser = argparse.ArgumentParser(description="Compare local RaMP SQLite files read-only; table deltas compare consecutive displayed databases.")
    parser.add_argument("--database", action="append", required=True, metavar="LABEL=PATH")
    parser.add_argument("--output", required=True, type=Path, help="HTML report path; matching JSON is also written")
    parser.add_argument('--core-exclude', action='append', default=[], metavar='LABEL',
                        help='Omit this database label from Core coverage and the Metabolites summary (repeatable)')
    args = parser.parse_args(argv)
    databases = []
    for value in args.database:
        label, sep, path = value.partition("=")
        if not sep or not label.strip() or not path:
            parser.error("Each --database must be LABEL=PATH")
        databases.append((label.strip(), Path(path)))
    if len(databases) < 2 or len({label for label, _ in databases}) != len(databases):
        parser.error("Provide at least two databases with unique labels")
    try:
        unknown = set(args.core_exclude) - {label for label, _ in databases}
        if unknown or len(set(args.core_exclude)) == len(databases):
            raise ValueError('Core exclusions must name supplied databases and leave at least one displayed database')
        report = compare(databases, progress=lambda s: print(s, flush=True))
        report['core_coverage_excluded_labels'] = args.core_exclude
        prepare_tables(report)
        for path in write_report(report, args.output):
            print(path)
    except (OSError, sqlite3.Error, ValueError) as exc:
        parser.exit(1, f"Comparison failed: {exc}\n")


if __name__ == "__main__":
    main()
