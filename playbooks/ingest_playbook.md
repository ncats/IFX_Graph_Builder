# New Data Source Ingest Playbook

## Goal
Provide a repeatable workflow for adding a new data source to the target graph ingest. Start with discovery, document decisions early, validate against real payloads, and only then implement the adapter. Promotion into broader workflows happens only after the working ingest is validated.

## Workflow Rules

- Read this playbook before starting any brand-new ingest source.
- Do discovery before code changes.
- Pause after discovery and propose a short implementation plan for user confirmation.
- Let the user run source acquisition and ETL executions unless they explicitly delegate those runs.
- Start in `src/use_cases/working.yaml`; only promote to `src/use_cases/pharos/target_graph.yaml` after validation.
- Keep source-specific decisions in a design doc under `designs/`.
- End with explicit validation instructions the user can run.

## Optional Comparison Inputs

- Raw input files remain the primary source of truth for payload shape.
- Old Pharos loader code, usually under `https://github.com/unmtransinfo/TCRD/tree/master/loaders`
- Current Pharos MySQL in the `pharos319` schema using `src/use_cases/secrets/pharos_credentials.yaml`
- New Pharos MySQL in the `pharos400` schema using `src/use_cases/secrets/pharos_credentials.yaml`
- Graph staging database on `ifxdev`, especially when a source may already land in the graph built by `build_pharos.py`

## Checklist

1) **Identify and acquire the source data**
   - Prefer an exact, versioned IFX Registry snapshot and pin its `data_source` in the working configuration.
   - If at all possible, get a discrete version and version date for the source.
   - For a source not yet in the Registry, document its acquisition and registration path before implementation.
   - The user should acquire or register required source files and report back when they are available, unless they explicitly delegate that step.
   - Record exact source URLs and snapshot IDs in a design doc.

2) **Capture version metadata**
   - Prefer an official version endpoint if available.
   - If not, use a stable proxy (e.g., `Last-Modified` header).
   - Write version data to a small file so adapters can reuse it.
   - Prefer persisting `version`, `version_date`, and `download_date` during download/prep rather than recomputing them inside the adapter.

3) **Explore and profile the downloaded files**
   - Before adapter implementation, run several iterative passes to inspect real payload shape.
   - Confirm field presence, cardinality, identifiers, edge semantics, and metadata quality.
   - Finalize inclusion/exclusion and mapping decisions based on observed data, not source documentation alone.
   - If the source is Pharos-relevant, compare the raw payload against both current downstream outputs and any existing graph staging output when available.

4) **Write the initial design document**
   - Create a source-specific document under `designs/` as soon as discovery starts.
   - Capture:
     - source URLs and file formats
     - version strategy
     - observed payload shape and identifier families
     - artifact classification when the input is derived: identity map, harmonizer annotation, source-faithful normalized data, or application/presentation output
     - provisional mapping and exclusion decisions
     - open questions or risks to validate during implementation
   - This design doc is the expected written artifact from discovery.

5) **Review previous ETLs for this data**
   - Ask the user about the TCRD ETL, or any later ETLs for the same data.
   - Treat prior ETLs as comparison points, not as the source of truth.
   - Note that the TCRD format is not always a natural fit, but often captures important historical scope.
   - For Pharos-related sources, inspect the old loader implementation in the TCRD repository when it helps explain legacy field choices or filtering.

6) **Check identifier normalization coverage**
   - Before adapter implementation, inspect how the configured resolver path will normalize the source IDs.
   - For Pharos / target_graph disease ingest, check the current Node Normalizer integration in `src/id_resolvers/node_normalizer.py`.
   - When the source offers multiple disease identifier families, profile each candidate family separately (for example `UMLS`, `SNOMEDCT`, `DOID`) rather than assuming the most ontology-like one is best.
   - Query the resolver service metadata when helpful, for example Node Normalizer `GET /get_curie_prefixes`, to confirm accepted prefixes.
   - Measure real coverage on distinct source IDs, not just a few spot checks.
   - Record both:
     - percent of source IDs that resolve at all
     - representative canonical prefixes returned by the resolver
   - Use these findings to choose what raw source ID the adapter should emit and leave canonicalization to the resolver layer whenever possible.
   - For harmonizer-produced or other pre-resolved inputs, read `designs/harmonizer_graph_builder_boundary.md` and determine separately:
     - which columns answer **"who is this?"** and belong in a resolver
     - which columns are genuinely producer-owned scores, methods, confidence, or review decisions and may be ingested through a harmonizer adapter
     - which fields or evidence originated in primary sources and should remain independently ingestible
   - Do not treat Registry registration or an `app_graph`, QA, explorer, or presentation export as evidence that the artifact is an appropriate adapter input.
   - If one proposed adapter would couple selection or versioning of multiple providers, pause and include that ownership change in the implementation plan for explicit approval.

7) **Review data that makes it into TCRD**
   - Currently Pharos uses pharos319.
   - Review the relevant tables and row counts to understand what was ingested previously.
   - Compare previous ingest output against the current raw payload to separate legacy limitations from current source reality.
   - When relevant, also inspect `pharos400` to understand what the newer MySQL path already captures or still misses.

8) **Pause and propose the implementation plan**
   - Summarize the intended adapter scope, node/edge model, resolver dependencies, and validation plan.
   - Keep the first pass intentionally minimal.
   - Get user confirmation before making code changes.

9) **Implement an InputAdapter**
   - Inherit from `src/interfaces/input_adapter.py` (or `FlatFileAdapter`).
   - Implement:
     - `get_all`
     - `get_datasource_name`
     - `get_version` (include `version`, `version_date`, `download_date`).
   - Emit `Node` / `Relationship` models that match the schema.
   - Keep adapters focused on source parsing and structural graph emission.
   - For harmonizer handoffs, keep identity mappings in a resolver and restrict the harmonizer adapter to assertions owned by the harmonizer. Do not re-label bundled primary-source fields as harmonizer provenance.
   - Distinguish preferred-annotation harmonization from evidence selection. Do not consume prefiltered, sampled, capped, ranked, or collapsed evidence edges unless the design explicitly documents every selection rule and demonstrates that provider identity, raw endpoints, row-level evidence, and independent source selection remain available.

10) **Map to the data model**
   - Confirm existing node/edge classes or add new ones in `src/models/`.
   - Use stable IDs and consistent prefixes.
   - Avoid speculative parsing when source text is ambiguous; preserve the source text when parsing would be lossy.
   - Keep source-specific payload that may merge later inside `details` structures instead of flattening it into top-level edge fields.

11) **Wire configuration into YAML**
   - Add the adapter to `src/use_cases/working.yaml` first.
   - Pass file paths and version metadata file paths via `kwargs`.
   - Only after the working ingest is validated, promote the finalized configuration into `src/use_cases/pharos/target_graph.yaml`.

12) **Validate the working ingest**
    - Ask the user to run the working ETL path.
    - Validate that counts, labels, IDs, provenance, and key edge endpoints look correct.
    - Validate that representative input-file records land where expected in the working graph and, when available, in the working MySQL output.
    - When there is a working MySQL validation path such as `src/use_cases/working_mysql.yaml`, compare working MySQL output against `pharos319` before promotion.
    - When relevant, compare the working graph against the graph staging database on `ifxdev` to confirm whether the source is already represented there.
    - Check both row counts and field population:
      - which destination tables received rows
      - which source-specific columns are populated in `pharos319` but still empty in the working MySQL output
      - whether graph data is present in the working graph but not yet mapped into downstream tables

13) **Update the design document**
    - Revise the design doc to reflect what actually ended up in the code:
      - Final field mappings and any decisions that changed during implementation
      - Actual node/edge counts produced
      - Any data quality issues encountered and how they were handled
      - Follow-up scope intentionally deferred from the first pass
    - Include the exact validation steps and comparison points used to accept the ingest.
    - Keep the final design doc concise and end-state focused.
    - Document where the implementation landed, not every exploratory branch or wrong turn from development.
    - If investigation notes are worth keeping, summarize only the conclusions that materially explain the final design.

14) **Update the playbook with reusable lessons learned**
   - After the ingest is working, capture any lessons that would help future ingests in the `Lessons Learned` section of this playbook.
   - Keep these lessons generally applicable across sources and workflows.
   - Prefer patterns, failure modes, and workflow improvements over source-specific facts.
   - Do not duplicate details that belong only in the source design doc.

---

## Lessons Learned

### File format
- **Always survey the file before writing the adapter.** Column names, headers, and ID formats
  frequently differ from documentation or old code. Use `head`, `cut`, and a quick Python profile
  script to confirm shape, cardinality, and value ranges before committing to a design.
- **Headerless files are common.** Pass explicit `fieldnames` to `csv.DictReader` rather than
  relying on a header row.
- **"Human" in the filename doesn't mean human-filtered.** E.g. JensenLab's
  `human_tissue_integrated_full.tsv` contains BTO tissues spanning all organisms. Profile
  the actual gene IDs and tissue ontology prefixes to understand scope.

### Gene / protein identifiers
- **Resolver coverage and graph-edge coverage can differ.** A resolver may
  index IDs that its source adapter does not emit as nodes (for example genes
  on noncanonical Ensembl sequences). Compare matched-ID counts with actual
  edge endpoints in the working graph, and align the resolver and source
  adapter scopes before describing a match as a linked graph entity.
- **Check what ID type the file actually uses.** JensenLab TISSUES uses ENSP (Ensembl protein)
  IDs mixed with gene symbols for non-coding entries — not the ENSG IDs the old TCRD code implied.
  Filter with `startswith("ENS")` to accept both ENSP and ENSG while dropping miRNAs, rRNAs, etc.
- **Don't emit a `Protein` node if there is no extra data to contribute.** When there are no
  `calculated_properties` or other protein-level fields to add, the edge's `start_node` is
  sufficient — a standalone `Protein` node just adds noise.

### Human-specificity: adapters vs. resolvers

- **Adapters are responsible for human filtering; resolvers are not.** A resolver's job is to map
  IDs to canonical equivalents — it doesn't skip unmatched entities, it just returns them as-is.
  Filtering to human-only data must happen in the adapter at emit time, before records enter the
  graph. If a non-human tissue node is emitted, it will land in the graph even without a UBERON
  mapping.
- **`uberon.obo` is multi-species, not human-only.** It covers vertebrates broadly (human, mouse,
  zebrafish, xenopus, etc.). A BTO ID appearing as an xref in the OBO does not guarantee it is
  human. Use the **FMA xref filter**: UBERON/CL terms with an FMA (Foundational Model of Anatomy)
  xref are human-grounded, because FMA is a human-specific ontology. This mirrors the old TCRD
  `getOntologyMap()` approach.
- **The FMA filter belongs in the adapter, not the resolver.** `TissueResolver` should resolve
  whatever it's given — it's a general-purpose tool used across multiple pipelines. Human-specificity
  is a data-source concern. Encode it in `_load_valid_tissue_ids()` (or equivalent) at adapter
  instantiation, not in the resolver logic.

### Tissue / ontology IDs
- **Verify that the ontology prefix appears in `uberon.obo` before adding it to
  `valid_ontologies`.** CLDB (BRENDA Cell Line Database) is a valid prefix in `constants.py` but
  has zero xrefs in UBERON, so adding it to the resolver list does nothing.
- **Watch for typos in source data ontology IDs.** JensenLab has a single erroneous `CLBD:0007212`
  entry (should be `CLDB`) — worth reporting upstream.

### Tau and tissue specificity
- **Tau requires a curated, flat, non-redundant tissue panel.** It breaks down when:
  - The tissue set is ontology-driven with multiple levels of granularity (e.g. "liver",
    "hepatocytes", and "liver parenchymal cell" all present as separate rows).
  - Cell lines and anatomical tissues are mixed in the same dataset.
  - HPM, HPA, and GTEx all use curated atlases where tau is well-defined; JensenLab TISSUES
    does not.

### Versioning
- **Use the `Last-Modified` HTTP header when no official version endpoint exists.** JensenLab
  regenerates files on a weekly Sunday schedule — the `Last-Modified` header captures that publish
  date reliably. Record it in source snapshot metadata so the adapter can use it as `version_date`.
- **Use named parameters for `DatasourceVersionInfo`.** This avoids argument-order regressions when version handling evolves.

### Edge merging and downstream parity
- **If duplicate source rows can collapse onto the same final graph edge, make merge-sensitive edge payload list-valued in the graph.** STRING PPIs needed `score: List[float]` instead of a scalar because reviewed-target collapse could merge multiple source rows onto one final `PPIEdge`; otherwise `KeepLast` silently overwrote earlier scores.
- **A graph edge and a downstream SQL row do not need the same cardinality.** For PPIs, the graph should keep one conceptual undirected edge with merged payload, while the legacy `ncats_ppi` table expects reciprocal directed rows for query-order parity. Model the graph cleanly first, then adapt row cardinality explicitly in the SQL converter.
- **When a legacy downstream table stores one scalar but the graph preserves many values, choose and document an explicit collapse rule.** For STRING PPI export, `max(score)` was used for `ncats_ppi.score`; do not let scalar selection happen implicitly through merge order.
