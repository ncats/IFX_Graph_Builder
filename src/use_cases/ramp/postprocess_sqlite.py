"""Repeatable post-processing of an existing RaMP SQLite export."""

import argparse
from collections import defaultdict
import json
from pathlib import Path
import sqlite3
import time


ENTITY_STATUS_COLUMNS = ('status_category', 'entity_source_id',
                         'entity_source_name', 'entity_count')
ENTITY_STATUS_DDL = '''CREATE TABLE IF NOT EXISTS entity_status_info (
    status_category VARCHAR(64) NOT NULL,
    entity_source_id VARCHAR(32) NOT NULL,
    entity_source_name VARCHAR(45) NOT NULL,
    entity_count INTEGER NOT NULL
)'''
DB_VERSION_COLUMNS = ('ramp_version', 'load_timestamp', 'version_notes',
                      'met_intersects_json', 'gene_intersects_json',
                      'met_intersects_json_pw_mapped',
                      'gene_intersects_json_pw_mapped', 'db_sql_url')
DB_VERSION_DDL = '''CREATE TABLE IF NOT EXISTS db_version (
    ramp_version VARCHAR(20) NOT NULL,
    load_timestamp DATETIME NOT NULL DEFAULT 'CURRENT_TIMESTAMP',
    version_notes VARCHAR(256),
    met_intersects_json VARCHAR(10000),
    gene_intersects_json VARCHAR(10000),
    met_intersects_json_pw_mapped VARCHAR(10000),
    gene_intersects_json_pw_mapped VARCHAR(10000),
    db_sql_url VARCHAR(256)
)'''
EMPTY_PATHWAY_TABLES = {
    'pathway_duplicates': (
        '''CREATE TABLE IF NOT EXISTS pathway_duplicates (
            pathwayRampId1 varchar(30) not null,
            pathwayRampId2 varchar(30) not null
        )''',
        ('pathwayRampId1', 'pathwayRampId2')),
    'pathway_similarity': (
        '''CREATE TABLE IF NOT EXISTS pathway_similarity (
            pathwayRampId varchar(30) not null primary key,
            analyte_blob blob,
            metabolite_blob blob,
            gene_blob blob,
            metabolite_count integer default 0,
            gene_count integer default 0
        )''',
        ('pathwayRampId', 'analyte_blob', 'metabolite_blob', 'gene_blob',
         'metabolite_count', 'gene_count')),
}
KEGG_SOURCE_ALIASES = {'hmdb_kegg': 'kegg', 'wikipathways_kegg': 'kegg'}
SOURCE_NAMES = {
    'kegg': 'KEGG', 'chebi': 'ChEBI', 'hmdb': 'HMDB',
    'reactome': 'Reactome', 'wiki': 'WikiPathways',
    'lipidmaps': 'LIPIDMAPS', 'rhea': 'Rhea',
    'pfocr': 'PFOCR', 'refmet': 'RefMet', 'pubchem': 'PubChem',
    'uniprot': 'UniProt', 'expasy': 'ExPASy ENZYME',
}
REQUIRED_TABLES = ('analyte', 'source', 'pathway', 'analytehaspathway',
                   'ontology', 'analytehasontology', 'chem_props',
                   'reaction2met', 'reaction2protein',
                   'catalyzed', 'version_info', 'ramp_export_metadata')


def _source_key(value, *, context):
    if not isinstance(value, str) or not value.strip():
        raise ValueError(f'{context} has a missing source identifier')
    return KEGG_SOURCE_ALIASES.get(value.strip().lower(), value.strip().lower())


def entity_status_rows(db):
    """Build legacy category counts from primary SQLite rows and stored source keys."""
    missing = [table for table in REQUIRED_TABLES if not db.execute(
        "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()]
    if missing:
        raise ValueError(f'Entity status requires missing primary tables: {", ".join(missing)}')

    counts = defaultdict(lambda: defaultdict(int))
    for kind, source, count in db.execute('''
        SELECT geneOrCompound,lower(trim(dataSource)),count(DISTINCT rampId)
        FROM source GROUP BY geneOrCompound,lower(trim(dataSource))'''):
        if kind not in ('compound', 'gene'):
            raise ValueError(f'Unsupported source analyte type: {kind!r}')
        category = 'Metabolites' if kind == 'compound' else 'Genes'
        key = _source_key(source, context='source.dataSource')
        if key == 'kegg':
            # Two upstream KEGG spellings can overlap on one RaMP analyte.
            continue
        counts[category][key] = count
    for kind, count in db.execute('''
        SELECT geneOrCompound,count(DISTINCT rampId) FROM source
        WHERE lower(trim(dataSource)) IN ('kegg','hmdb_kegg','wikipathways_kegg')
        GROUP BY geneOrCompound'''):
        counts['Metabolites' if kind == 'compound' else 'Genes']['kegg'] = count

    for source, count in db.execute('''
        SELECT lower(trim(type)),count(DISTINCT pathwayRampId)
        FROM pathway GROUP BY lower(trim(type))'''):
        counts['Pathways'][_source_key(source, context='pathway.type')] = count
    for kind, source, count in db.execute('''
        SELECT a.type,lower(trim(p.pathwaySource)),count(*)
        FROM analytehaspathway p JOIN analyte a ON a.rampId=p.rampId
        GROUP BY a.type,lower(trim(p.pathwaySource))'''):
        if kind not in ('compound', 'gene'):
            raise ValueError(f'Unsupported pathway analyte type: {kind!r}')
        category = ('Metabolite-Pathway Associations' if kind == 'compound'
                    else 'Gene-Pathway Associations')
        counts[category][_source_key(source, context='analytehaspathway.pathwaySource')] += count
    for source, count in db.execute('''
        SELECT lower(trim(chem_data_source)),count(DISTINCT chem_source_id)
        FROM chem_props GROUP BY lower(trim(chem_data_source))'''):
        counts['Chemical Property Records'][_source_key(source, context='chem_props.chem_data_source')] = count

    for table, category in (('reaction2met', 'Metabolite-Reaction Associations'),
                            ('reaction2protein', 'Gene-Reaction Associations')):
        for source, count in db.execute(f'''
            SELECT lower(substr(rxn_source_id,1,instr(rxn_source_id,':')-1)),count(*)
            FROM {table} GROUP BY 1'''):
            counts[category][_source_key(source, context=f'{table}.rxn_source_id')] = count

    # catalyzed has no source column: its present exporter contract is HMDB.
    counts['Metabolite-Gene Associations']['hmdb'] = db.execute('''
        SELECT count(*) FROM (SELECT DISTINCT rampCompoundId,rampGeneId FROM catalyzed)''').fetchone()[0]
    for source, count in db.execute('''
        SELECT source,count(*) FROM (
            SELECT DISTINCT lower(substr(rm.rxn_source_id,1,instr(rm.rxn_source_id,':')-1)) source,
                            rm.ramp_cmpd_id,rp.ramp_gene_id
            FROM reaction2met rm JOIN reaction2protein rp
              ON rm.ramp_rxn_id=rp.ramp_rxn_id
             AND rm.rxn_source_id=rp.rxn_source_id
        ) GROUP BY source'''):
        counts['Metabolite-Gene Associations'][_source_key(source, context='reaction participant source')] += count

    version_names = {}
    for source, name in db.execute("SELECT data_source_id,data_source_name FROM version_info WHERE status='current'"):
        key = _source_key(source, context='version_info.data_source_id')
        if key in version_names and version_names[key] != name:
            raise ValueError(f'Conflicting current display names for {key}')
        version_names[key] = name
    rows = []
    for category, by_source in sorted(counts.items()):
        names_in_category = {}
        for source, count in sorted(by_source.items()):
            name = SOURCE_NAMES.get(source) or version_names.get(source) or source
            if name in names_in_category:
                raise ValueError(f'{category}: {source} and {names_in_category[name]} share display name {name!r}; define an intentional source alias')
            names_in_category[name] = source
            rows.append((category, source, name, count))
    return rows


def update_lookup_counts(db):
    """Fill the two derived lookup counts using the legacy counting rules."""
    db.execute('''CREATE TEMP TABLE _ramp_pathway_counts AS
        SELECT rampId,count(DISTINCT pathwayRampId) AS n
        FROM analytehaspathway WHERE lower(trim(pathwaySource))!='hmdb'
        GROUP BY rampId''')
    db.execute('CREATE UNIQUE INDEX _ramp_pathway_counts_id ON _ramp_pathway_counts(rampId)')
    db.execute('''UPDATE source SET pathwayCount = coalesce((
        SELECT n FROM _ramp_pathway_counts p WHERE p.rampId=source.rampId
    ),0)''')
    db.execute('DROP TABLE _ramp_pathway_counts')
    db.execute('''UPDATE ontology SET metCount = coalesce((
        SELECT count(DISTINCT a.rampCompoundId)
        FROM analytehasontology a WHERE a.rampOntologyId=ontology.rampOntologyId
    ),0)''')


def db_version_intersections(db):
    """Return R-compatible exclusive source intersections from stored identities."""
    mapped = {row[0] for row in db.execute('''
        SELECT DISTINCT a.rampId FROM analytehaspathway a
        JOIN pathway p ON p.pathwayRampId=a.pathwayRampId
        WHERE p.pathwayCategory IS NULL
           OR (lower(trim(p.pathwayCategory)) NOT IN ('smpdb2','smpdb3')
               AND lower(trim(p.pathwayCategory)) NOT LIKE '%:smpdb2'
               AND lower(trim(p.pathwayCategory)) NOT LIKE '%:smpdb3')''')}
    memberships = {'compound': {}, 'gene': {}}
    for rid, kind, source in db.execute('''
        SELECT rampId,geneOrCompound,lower(trim(dataSource))
        FROM source GROUP BY rampId,geneOrCompound,lower(trim(dataSource))'''):
        if kind not in memberships:
            raise ValueError(f'Unsupported source analyte type: {kind!r}')
        if not isinstance(rid, str) or not rid:
            raise ValueError('source contains a missing RaMP ID')
        memberships[kind].setdefault(rid, set()).add(
            _source_key(source, context='source.dataSource'))

    display = dict(SOURCE_NAMES)
    version_names = {}
    for source, name in db.execute(
            "SELECT data_source_id,data_source_name FROM version_info WHERE status='current'"):
        key = _source_key(source, context='version_info.data_source_id')
        name = name.strip() if isinstance(name, str) and name.strip() else key
        if key in version_names and version_names[key] != name:
            raise ValueError(f'Conflicting source display name for {key}')
        version_names[key] = name
        display.setdefault(key, name)
    used_sources = {source for by_id in memberships.values()
                    for sources in by_id.values() for source in sources}
    labels_to_keys = defaultdict(list)
    for source in sorted(used_sources):
        labels_to_keys[display.get(source, source)].append(source)
    collisions = {label: keys for label, keys in labels_to_keys.items() if len(keys) > 1}
    if collisions:
        raise ValueError(f'Source display names are ambiguous for UpSet intersections: {collisions}')
    columns = {}
    summary = {}
    for kind, prefix, global_col, mapped_col in (
            ('compound', 'cmpd', 'met_intersects_json', 'met_intersects_json_pw_mapped'),
            ('gene', 'gene', 'gene_intersects_json', 'gene_intersects_json_pw_mapped')):
        for scope, column in (('global', global_col), ('mapped-to-pathway', mapped_col)):
            counts = defaultdict(int)
            eligible = 0
            for rid, sources in memberships[kind].items():
                if scope == 'mapped-to-pathway' and rid not in mapped:
                    continue
                labels = tuple(sorted({display.get(source, source) for source in sources}))
                if len(labels) != len(sources):
                    raise ValueError(f'{kind} source keys collapse to duplicate display names: {sources}')
                counts[labels] += 1
                eligible += 1
            records = [{'id': f'{prefix}_src_set_{i}', 'sets': list(labels), 'size': count}
                       for i, (labels, count) in enumerate(sorted(counts.items()), 1)]
            if sum(row['size'] for row in records) != eligible:
                raise ValueError(f'{column} intersections do not partition eligible RaMP IDs')
            columns[column] = json.dumps(records, sort_keys=True, separators=(',', ':'))
            summary[column] = {'scope': scope, 'analyte_type': kind,
                               'analytes': eligible, 'intersections': len(records),
                               'sources': sorted({source for labels in counts for source in labels})}
    return columns, summary


def update_db_version(db, manifest, intersections):
    """Create or update one current version row without changing its timestamp."""
    db.execute(DB_VERSION_DDL)
    columns = tuple(row[1] for row in db.execute('PRAGMA table_info(db_version)'))
    if columns != DB_VERSION_COLUMNS:
        raise ValueError(f'Unexpected db_version columns: {columns}')
    existing = db.execute('SELECT rowid,ramp_version,load_timestamp FROM db_version').fetchall()
    if len(existing) > 1:
        raise ValueError(f'Expected one db_version row, found {len(existing)}')
    version, timestamp = manifest['release_version'], manifest['created_at']
    if existing:
        rowid, stored_version, stored_timestamp = existing[0]
        if (stored_version, stored_timestamp) != (version, timestamp):
            raise ValueError('db_version release or timestamp differs from export manifest')
    else:
        db.execute('''INSERT INTO db_version (ramp_version,load_timestamp,version_notes)
            VALUES (?,?,?)''',
            (version, timestamp,
             'Lookup-table diagnostic; not a release database. See ramp_export_metadata.'
             if manifest['scope'] == 'lookup_table_diagnostic' else
             'Base tables; post-processing status is in ramp_export_metadata.'))
        rowid = db.execute('SELECT last_insert_rowid()').fetchone()[0]
    fields = ('met_intersects_json', 'gene_intersects_json',
              'met_intersects_json_pw_mapped', 'gene_intersects_json_pw_mapped')
    db.execute('''UPDATE db_version SET met_intersects_json=?,gene_intersects_json=?,
        met_intersects_json_pw_mapped=?,gene_intersects_json_pw_mapped=? WHERE rowid=?''',
        (*(intersections[field] for field in fields), rowid))
    if db.execute('SELECT count(*) FROM db_version').fetchone()[0] != 1:
        raise ValueError('db_version row count changed during update')


def ensure_empty_pathway_tables(db):
    """Expose the legacy table shape and preserve any already computed rows."""
    counts = {}
    for name, (ddl, expected_columns) in EMPTY_PATHWAY_TABLES.items():
        db.execute(ddl)
        columns = tuple(row[1] for row in db.execute(f'PRAGMA table_info({name})'))
        if columns != expected_columns:
            raise ValueError(f'Unexpected {name} columns: {columns}')
        counts[name] = db.execute(f'SELECT count(*) FROM {name}').fetchone()[0]
    return counts


def postprocess_pathway_similarity(path, *, progress=None):
    """Atomically replace legacy compressed overlaps and exact duplicate pairs."""
    from src.use_cases.ramp.pathway_similarity import populate_temp_tables

    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f'SQLite input must be an existing ordinary file: {path}')
    started = time.perf_counter()
    with sqlite3.connect(path, timeout=30) as db:
        db.execute('PRAGMA busy_timeout=30000')
        mode = db.execute('PRAGMA journal_mode').fetchone()[0].lower()
        if mode not in ('delete', 'truncate', 'persist'):
            raise ValueError(f'Post-processing requires rollback-journal mode, found {mode}')
        db.execute('PRAGMA synchronous=FULL')
        try:
            db.execute('BEGIN IMMEDIATE')
            raw = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
            if raw is None:
                raise ValueError('SQLite has no RaMP export manifest')
            manifest = json.loads(raw[0])
            if manifest.get('scope') not in ('lookup_table_diagnostic', 'base_tables'):
                raise ValueError(f'Unsupported RaMP export scope: {manifest.get("scope")!r}')
            missing = [table for table in ('pathway', 'analytehaspathway') if not db.execute(
                "SELECT 1 FROM sqlite_master WHERE type='table' AND name=?", (table,)).fetchone()]
            if missing:
                raise ValueError(f'Pathway similarity requires missing tables: {", ".join(missing)}')
            ensure_empty_pathway_tables(db)
            summary = populate_temp_tables(db, progress=progress)
            if progress:
                progress('Replacing pathway_similarity and pathway_duplicates atomically')
            db.execute('DELETE FROM pathway_similarity')
            db.execute('DELETE FROM pathway_duplicates')
            db.execute('''INSERT INTO pathway_similarity
                (pathwayRampId,analyte_blob,metabolite_blob,gene_blob,
                 metabolite_count,gene_count)
                SELECT pathwayRampId,analyte_blob,metabolite_blob,gene_blob,
                       metabolite_count,gene_count
                FROM _ramp_similarity_new ORDER BY pathwayRampId''')
            db.execute('''INSERT INTO pathway_duplicates
                (pathwayRampId1,pathwayRampId2)
                SELECT pathwayRampId1,pathwayRampId2 FROM _ramp_duplicates_new
                ORDER BY pathwayRampId1,pathwayRampId2''')
            counts = {name: db.execute(f'SELECT count(*) FROM {name}').fetchone()[0]
                      for name in EMPTY_PATHWAY_TABLES}
            if counts != {'pathway_similarity': summary['rows'],
                          'pathway_duplicates': summary['duplicate_pairs']}:
                raise ValueError(f'Pathway result row counts differ after insertion: {counts}')
            scope_columns = {'analyte': 'analyte_blob', 'metabolite': 'metabolite_blob',
                             'gene': 'gene_blob'}
            for scope, column in scope_columns.items():
                count = db.execute(f'SELECT count(*) FROM pathway_similarity WHERE {column} IS NOT NULL').fetchone()[0]
                if count != summary[scope]['eligible_pathways']:
                    raise ValueError(f'{scope} blob coverage differs from eligible pathways: {count}')
            summary['compressed_bytes'] = db.execute('''SELECT coalesce(sum(
                coalesce(length(analyte_blob),0)+coalesce(length(metabolite_blob),0)
                +coalesce(length(gene_blob),0)),0) FROM pathway_similarity''').fetchone()[0]
            manifest.setdefault('row_counts', {}).update(counts)
            manifest['pending_tables'] = [name for name in manifest.get('pending_tables', [])
                                          if name not in EMPTY_PATHWAY_TABLES]
            included = manifest.get('included_tables')
            for name in EMPTY_PATHWAY_TABLES:
                if included is not None and name not in included:
                    included.append(name)
            manifest.setdefault('post_processing', {})['pathway_tables'] = {
                'status': 'computed',
                'policy': ('Legacy non-HMDB pathway scope (equivalent to excluding SMPDB in this build); '
                           'combined minimum 10 IDs, metabolite/gene minimum 5; exact Jaccard rounded '
                           'to thousandths, positive pairs delta-indexed and zlib-compressed; '
                           'duplicates are identical combined-analyte sets'),
                'summary': summary,
            }
            if (manifest.get('scope') == 'base_tables'
                    and not manifest.get('pending_tables')
                    and not manifest.get('pending_fields')
                    and manifest['post_processing'].get('human_reaction_flags')):
                manifest['status'] = 'complete'
                updated = db.execute(
                    "UPDATE db_version SET version_notes=?",
                    ('Post-processing complete. See ramp_export_metadata.',))
                if updated.rowcount != 1:
                    raise ValueError('Complete RaMP SQLite requires one db_version row')
            db.execute("UPDATE ramp_export_metadata SET value=? WHERE key='manifest'",
                       (json.dumps(manifest, sort_keys=True),))
            if db.execute('PRAGMA foreign_key_check').fetchone():
                raise ValueError('SQLite foreign-key check failed')
            if db.execute('PRAGMA quick_check').fetchone() != ('ok',):
                raise ValueError('SQLite quick check failed')
            db.commit()
        except Exception:
            db.rollback()
            raise
    return {**summary, 'elapsed_seconds': time.perf_counter() - started}


def graph_human_flags(manifest, *, database, graph_credentials,
                      registry_credentials, registry_cache_dir):
    """Read human-scope evidence from the exact graph and UniProt pin."""
    import yaml
    from src.core.graph_build_identity import source_build_fingerprint
    from src.core.registry_integration import RegistryIntegration
    from src.shared.arango_adapter import ArangoAdapter
    from src.shared.db_credentials import DBCredentials
    from src.use_cases.ramp.sqlite_human_reactions import (
        human_primary_accessions, human_reaction_flags,
    )

    credentials = DBCredentials.from_yaml(yaml.safe_load(graph_credentials.read_text()))
    graph = ArangoAdapter(credentials, database).get_db()
    stage_id = manifest['stage_id'].removeprefix('HarmonizationStage:')
    stage = graph.collection('HarmonizationStage').get(stage_id)
    if not stage or stage.get('status') != 'complete':
        raise ValueError(f'Completed stage is missing from graph: {stage_id}')
    if stage.get('_rev') != manifest.get('stage_revision'):
        raise ValueError('Graph stage revision differs from SQLite export')
    metadata = graph.collection('metadata_store').get('etl_metadata')['value']
    if source_build_fingerprint(metadata) != source_build_fingerprint(manifest['graph_build']):
        raise ValueError('Graph input metadata differs from SQLite export')
    relevant = ('RheaReaction', 'RheaMetaboliteReactionEdge', 'RheaProteinReactionEdge',
                'BiologicalRole', 'IsAEdge', 'HasBiologicalRoleEdge')
    for collection in relevant:
        if graph.collection(collection).revision() != manifest['source_revisions'].get(collection):
            raise ValueError(f'{collection} revision differs from SQLite export')
    pin = manifest['gene_resolution']
    registry = RegistryIntegration.connect({'credentials': registry_credentials,
                                            'cache_dir': registry_cache_dir})
    dataset = registry.resolve(pin['snapshot_id'])
    accessions = human_primary_accessions(dataset.file(pin['file']), pin['sha256'])
    flags, summary = human_reaction_flags(graph, accessions)
    for collection in relevant:
        if graph.collection(collection).revision() != manifest['source_revisions'][collection]:
            raise ValueError(f'{collection} changed while deriving human reaction flags')
    return flags, {**summary, 'graph_stage_id': manifest['stage_id'],
                   'policy': 'RaMP legacy reaction scope: ChEBI descendants of five human/drug-metabolite roots via is_a or has_biological_role; all reported ChEBI participants required; any pinned human UniProt primary accession qualifies',
                   'uniprot_snapshot_id': pin['snapshot_id'],
                   'uniprot_file': pin['file'], 'uniprot_sha256': pin['sha256'],
                   'chebi_roots': sorted((
                       'CHEBI:77746', 'CHEBI:85234', 'CHEBI:84087',
                       'CHEBI:76967', 'CHEBI:49103'))}


def postprocess_entity_status(path, *, human_flags=None, human_manifest=None):
    """Replace derived status rows transactionally; preserve the database on failure."""
    path = Path(path)
    if path.is_symlink() or not path.is_file():
        raise ValueError(f'SQLite input must be an existing ordinary file: {path}')
    started = time.perf_counter()
    with sqlite3.connect(path, timeout=30) as db:
        db.execute('PRAGMA busy_timeout=30000')
        mode = db.execute('PRAGMA journal_mode').fetchone()[0].lower()
        if mode not in ('delete', 'truncate', 'persist'):
            raise ValueError(f'Post-processing requires rollback-journal mode, found {mode}')
        db.execute('PRAGMA synchronous=FULL')
        try:
            db.execute('BEGIN IMMEDIATE')
            raw = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
            if raw is None:
                raise ValueError('SQLite has no RaMP export manifest')
            manifest = json.loads(raw[0])
            if manifest.get('scope') not in ('lookup_table_diagnostic', 'base_tables'):
                raise ValueError(f'Unsupported RaMP export scope: {manifest.get("scope")!r}')
            rows = entity_status_rows(db)
            db.execute(ENTITY_STATUS_DDL)
            columns = tuple(row[1] for row in db.execute('PRAGMA table_info(entity_status_info)'))
            if columns != ENTITY_STATUS_COLUMNS:
                raise ValueError(f'Unexpected entity_status_info columns: {columns}')
            db.execute('DELETE FROM entity_status_info')
            db.executemany('INSERT INTO entity_status_info VALUES (?,?,?,?)', rows)
            if db.execute('SELECT count(*) FROM entity_status_info').fetchone()[0] != len(rows):
                raise ValueError('Entity status row count changed during insertion')
            update_lookup_counts(db)
            intersections, intersection_summary = db_version_intersections(db)
            update_db_version(db, manifest, intersections)
            pathway_table_counts = ensure_empty_pathway_tables(db)
            if human_flags is not None:
                if human_manifest is None:
                    raise ValueError('Human reaction flags require provenance metadata')
                stored = {source: rid for rid, source in db.execute(
                    'SELECT ramp_rxn_id,rxn_source_id FROM reaction')}
                missing = sorted(set(stored) - {'rhea:' + key.split(':', 1)[1]
                                               for key in human_flags})
                if missing:
                    raise ValueError(f'Human flags lack SQLite reactions: {missing[:5]}')
                db.executemany('''UPDATE reaction SET has_human_prot=?,only_human_mets=?
                    WHERE ramp_rxn_id=?''', [(*human_flags['RHEA:' + source.split(':', 1)[1]], rid)
                                           for source, rid in stored.items()])
            manifest.setdefault('row_counts', {})['entity_status_info'] = len(rows)
            manifest['row_counts']['db_version'] = 1
            manifest['row_counts'].update(pathway_table_counts)
            manifest['pending_tables'] = [name for name in manifest.get('pending_tables', [])
                                          if name != 'entity_status_info']
            included = manifest.get('included_tables')
            if included is not None and 'entity_status_info' not in included:
                included.append('entity_status_info')
            if included is not None and 'db_version' not in included:
                included.append('db_version')
            for name in EMPTY_PATHWAY_TABLES:
                if included is not None and name not in included:
                    included.append(name)
            manifest.setdefault('post_processing', {})['entity_status_info'] = {
                'policy': 'Counts from stored source-bearing tables; KEGG source aliases combined; catalyzed attributed to HMDB',
                'row_count': len(rows),
            }
            manifest['post_processing']['lookup_counts'] = {
                'policy': 'source.pathwayCount counts distinct non-HMDB pathways per RaMP ID; ontology.metCount counts distinct associated metabolites',
            }
            manifest['post_processing']['db_version_intersections'] = {
                'policy': 'Exact nonempty source memberships from source; KEGG aliases normalized before grouping; pathway scope excludes smpdb2/smpdb3',
                'summary': intersection_summary,
            }
            if manifest['post_processing'].get('pathway_tables', {}).get('status') != 'computed':
                manifest['post_processing']['pathway_tables'] = {
                    'status': 'schema_only' if not any(pathway_table_counts.values()) else 'existing_rows_preserved',
                    'policy': 'Legacy pathway table schemas present; this post-processing pass does not calculate or replace their rows',
                    'row_counts': pathway_table_counts,
                }
            if human_flags is not None:
                manifest['post_processing']['human_reaction_flags'] = human_manifest
            manifest['pending_fields'] = [field for field in manifest.get('pending_fields', [])
                                          if field not in (('source.pathwayCount', 'ontology.metCount',
                                                           'reaction.has_human_prot', 'reaction.only_human_mets',
                                                           'db_version.*intersects*')
                                                           if human_flags is not None else
                                                           ('source.pathwayCount', 'ontology.metCount',
                                                            'db_version.*intersects*'))]
            db.execute("UPDATE ramp_export_metadata SET value=? WHERE key='manifest'",
                       (json.dumps(manifest, sort_keys=True),))
            if db.execute('PRAGMA foreign_key_check').fetchone():
                raise ValueError('SQLite foreign-key check failed')
            if db.execute('PRAGMA quick_check').fetchone() != ('ok',):
                raise ValueError('SQLite quick check failed')
            source_rows = db.execute('SELECT count(*) FROM source').fetchone()[0]
            ontology_rows = db.execute('SELECT count(*) FROM ontology').fetchone()[0]
            reaction_rows = db.execute('SELECT count(*) FROM reaction').fetchone()[0]
            db.commit()
        except Exception:
            db.rollback()
            raise
    return {'rows': len(rows), 'source_rows': source_rows,
            'ontology_rows': ontology_rows,
            'reaction_rows': reaction_rows if human_flags is not None else None,
            'intersections': intersection_summary,
            'elapsed_seconds': time.perf_counter() - started}


def main(argv=None):
    parser = argparse.ArgumentParser(description='Rerun RaMP SQLite post-processing without rebuilding primary tables')
    parser.add_argument('--sqlite', type=Path, required=True, help='Existing RaMP SQLite to update')
    parser.add_argument('--pathway-only', action='store_true',
                        help='Recompute compressed pathway similarity and duplicate pairs only')
    parser.add_argument('--with-human-flags', action='store_true',
                        help='Also derive reaction human flags from the matching graph and pinned UniProt file')
    parser.add_argument('--database', default='metabolite_harmonization')
    parser.add_argument('--graph-credentials', type=Path,
                        default=Path('src/use_cases/secrets/ifxdev_arangodb.yaml'))
    parser.add_argument('--registry-credentials', type=Path,
                        default=Path('src/use_cases/secrets/aws_ifx_registry.yaml'))
    parser.add_argument('--registry-cache-dir', type=Path, default=Path('/var/tmp/ifx-registry-cache'))
    args = parser.parse_args(argv)
    if args.pathway_only:
        if args.with_human_flags:
            parser.error('--pathway-only cannot be combined with --with-human-flags')
        result = postprocess_pathway_similarity(args.sqlite, progress=print)
        print(f"Updated pathway_similarity: {result['rows']:,} rows")
        print(f"Updated pathway_duplicates: {result['duplicate_pairs']:,} pairs")
        print(f"Compressed similarity data: {result['compressed_bytes'] / 1048576:.1f} MiB")
        print(f"Pathway post-processing completed in {result['elapsed_seconds']:.1f}s")
        return 0
    flags = provenance = None
    if args.with_human_flags:
        with sqlite3.connect(f'{args.sqlite.resolve().as_uri()}?mode=ro', uri=True) as db:
            row = db.execute("SELECT value FROM ramp_export_metadata WHERE key='manifest'").fetchone()
        if row is None:
            raise ValueError('SQLite has no RaMP export manifest')
        flags, provenance = graph_human_flags(json.loads(row[0]), database=args.database,
            graph_credentials=args.graph_credentials, registry_credentials=args.registry_credentials,
            registry_cache_dir=args.registry_cache_dir)
    result = postprocess_entity_status(args.sqlite, human_flags=flags, human_manifest=provenance)
    print(f"Updated entity_status_info: {result['rows']:,} rows")
    print(f"Updated source.pathwayCount: {result['source_rows']:,} source rows")
    print(f"Updated ontology.metCount: {result['ontology_rows']:,} ontology rows")
    print('Updated db_version: 4 source-intersection JSON fields')
    if result['reaction_rows'] is not None:
        print(f"Updated reaction human-scope flags: {result['reaction_rows']:,} reactions")
    print(f"Post-processing completed in {result['elapsed_seconds']:.1f}s")
    return 0


if __name__ == '__main__':
    raise SystemExit(main())
