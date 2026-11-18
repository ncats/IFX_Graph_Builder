# Data Update Playbook

## Goal
Provide a repeatable workflow for refreshing existing input files, validating payload drift, and rebuilding a fresh graph/MySQL database with updated source data.

This playbook is for **updating existing sources**, not adding brand new ingest sources. Use `playbooks/ingest_playbook.md` for new-source development.

## Workflow Rules

- Start with payload validation, not code changes.
- Let the user run source acquisition and ETL executions unless they explicitly delegate those runs.
- Prefer a narrow working validation path before any broader rebuild.
- Update the relevant design doc when a refresh exposes a new parsing rule, mapping change, or failure mode.
- End with explicit validation instructions or follow-up rebuild steps for the user.

## Optional Comparison Inputs

- Refreshed raw input files are the primary evidence for payload drift.
- Old Pharos loader code when legacy behavior is relevant, usually under `https://github.com/unmtransinfo/TCRD/tree/master/loaders`
- Current Pharos MySQL in the `pharos319` schema using `src/use_cases/secrets/pharos_credentials.yaml`
- New Pharos MySQL in the `pharos400` schema using `src/use_cases/secrets/pharos_credentials.yaml`
- Graph staging database on `ifxdev` when the refreshed source may already be landing there

## Typical Use Case

- One or more `input_files/auto/<source>/` payloads were refreshed
- Multiple auto-downloaded sources may have changed at once
- The goal is to build a fresh working or production-style database and catch schema, parser, resolver, and downstream conversion problems early

## Checklist

1) **Record what changed**
   - List the sources being refreshed.
   - Record file names, versions, version dates, and download dates.
   - Confirm which files are expected to change and which should remain stable.

2) **Validate that the refresh actually happened**
   - Check file timestamps, sizes, and version metadata files.
   - Compare old vs new row counts or record counts where practical.
   - Watch for suspiciously small files, empty files, truncated archives, or stale version files.

3) **Profile the payload shape before rebuilding**
   - Check column names or top-level fields.
   - Check identifier formats and prefixes.
   - Check representative values, null rates, duplicates, and category drift.
   - Check whether multi-value fields changed delimiter or structure.
   - Check whether any previously expected columns disappeared or renamed.

4) **Review data drift, not just file presence**
   - Look for new categories, evidence codes, relationship types, ontology prefixes, or species labels.
   - Look for changes in cardinality:
     one-to-one becoming one-to-many
     single parent becoming multi-parent
     unique rows becoming duplicate-prone
   - Look for free-text normalization drift:
     capitalization
     punctuation
     embedded suffixes
     identifier version suffixes

5) **Validate assumptions used by adapters and resolvers**
   - Confirm the adapter still reads the right columns/fields.
   - Confirm any resolver input IDs still match the expected identifier family.
   - Confirm joins still line up across node and edge endpoints.
   - Confirm source-specific parsing is still lossless enough for the payload.

6) **Run a small working build first**
   - Start with `src/use_cases/working.yaml` or another intentionally narrow config.
   - Prefer validating a minimal set of affected sources before a full rebuild.
   - For bulk source refreshes, use a working config that exercises the affected source families.
   - Ask the user to run the working ETL and report back with the relevant failure or validation output.

7) **Inspect the working graph directly**
   - Check collection counts.
   - Check key node and edge collections touched by the refresh.
   - Check endpoint integrity:
     edge start IDs resolve to nodes
     edge end IDs resolve to nodes
   - Check representative samples for expected provenance, names, IDs, and details.
   - Check whether post-processing removed dangling edges as expected.
   - Check whether representative raw input records can be traced into the working graph.

8) **Run the matching working MySQL conversion when available**
   - Use `src/use_cases/working_mysql.yaml` when there is a downstream TCRD validation path.
   - Treat every `IntegrityError` as actionable until explained.
   - Compare the refreshed output against both `pharos319` and, when relevant, `pharos400`.
   - Common causes:
     stale non-truncated tables
     duplicate emission across batches
     schema too restrictive for source cardinality
     unresolved graph references

9) **Only then run the broader rebuild**
   - For a broad auto-source refresh, rebuild the full target database only after the working path is clean.
   - Expect this to be a mostly automated pipeline for producing a fresh database from updated inputs.
   - Keep the refresh process reproducible:
     refresh inputs
     validate payload drift
     run working validation
     run full build

10) **Document what broke and what was fixed**
   - Record every adapter, schema, resolver, or converter issue exposed by the refresh.
   - Note whether the issue was:
     bad source assumption
     payload drift
     stale downstream schema
     hidden duplicate previously masked by permissive inserts
   - Update this playbook or source-specific designs when a new failure mode is discovered.
   - Record the exact validation or rebuild commands the user should run next.

## High-Value Checks

- Column names still match adapter expectations
- Data types still match parser expectations
- Identifier families still match resolver expectations
- Distinct IDs still stay distinct after normalization
- Edge endpoints still correspond to real nodes
- Parent-child tables still allow real source cardinality
- Converters dedupe correctly across the whole run, not just one batch
- Fresh rebuild scripts actually truncate when they claim to
- Downstream inserts fail loudly instead of silently dropping bad rows

## Common Failure Modes

- Source file updated but version metadata did not
- Column rename or silent field removal
- New delimiter or changed nested structure inside a text field
- Identifier prefix drift
- Case/punctuation drift causing node-edge mismatches
- New multi-parent or multi-valued relationships exposing overly strict primary keys
- Duplicate associations across batches exposing missing converter dedupe
- Partial working graphs leaving dangling edges before downstream conversion
- Rebuild scripts preserving old rows when a clean rebuild was intended

## Recommended Mindset

- Assume input updates can break previously valid assumptions.
- Validate real payloads before trusting old code paths.
- Prefer narrow working validation before expensive full rebuilds.
- Treat new hard failures as useful signals, especially when permissive insert behavior previously hid bad data.
- Diagnose the root cause before landing defensive guards that might hide the real payload or schema problem.
