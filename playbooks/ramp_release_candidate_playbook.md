# RaMP SQLite release-candidate playbook

Use this after a completed, validated harmonization stage has produced a full
SQLite through `src.use_cases.ramp.build_sqlite`. Package the final-version
SQLite locally in `RaMP-DB/db/` for review on a GitHub branch. The database
bytes committed to that branch are the intended release bytes; merging the
branch into `main` is the official release event. Keep the original build so
packaging can be repeated without another graph export.

## 1. Freeze and verify the build

Record the stage ID, graph build revision, pinned source snapshots, and exporter
options from `ramp_export_metadata`; compute the SQLite SHA-256 separately.
Confirm the manifest
has `scope=base_tables`, `status=complete`, no pending tables or fields, computed
pathway tables, and human-reaction flags. Run `PRAGMA integrity_check` and
`PRAGMA foreign_key_check`. `package_release_candidate.py` repeats these checks
and refuses an incomplete build. Never relabel an incomplete diagnostic export.

Keep the comparison report and a fresh RaMP owner review with the candidate.
Compare against released 2.4.0 and 3.0.7 using the same R-package calls and
source IDs, since RaMP IDs are assigned afresh. Record differences that affect
users, including outputs that improve. The review of the October 2026 build is
[`output_files/ramp/ramp_fresh_product_review_20261013.md`](../output_files/ramp/ramp_fresh_product_review_20261013.md).

## 2. Package a local candidate

Choose a version that is not already present in `../RaMP-DB/db/`. For the
current candidate, use `3.0.13`:

```bash
.venv/bin/python -m src.use_cases.ramp.package_release_candidate \
  --sqlite output_files/ramp/ramp-base.sqlite \
  --version 3.0.13 \
  --output ../RaMP-DB/db/RaMP_SQLite_v3.0.13.sqlite.gz
```

The command copies the SQLite, changes `db_version.ramp_version`, every
`version_info.ramp_db_version`, and `ramp_export_metadata.release_version` in
the copy, then verifies the result and compresses it. The original build's
stage ID, source pins, and `created_at` remain unchanged. The manifest marks
the verified version `release_ready=true`; that means ready to release, not
already published. Packaging time and checksums are printed in an external
receipt, so no branch or candidate-specific information enters the SQLite.
The command prints the original, final SQLite, and gzip SHA-256 checksums.
Save that output with the review. It refuses an existing output, so a revised
candidate needs a new version or deliberate removal of an uncommitted local
artifact after review.

`../RaMP-DB/.gitattributes` already tracks `*.sqlite.gz` with Git LFS. Verify
`git -C ../RaMP-DB status --short` shows only the intended new database file
among the work for this candidate. Leave staging and committing to the owner.

## 3. Exercise the exact packaged bytes

Decompress the candidate into temporary storage and check its SQLite integrity,
foreign keys, manifest, and three matching version labels. Use the local file
explicitly with `RaMP:::.RaMP(dbname=path)`; the normal `RaMP()` constructor
discovers remote `main` and cannot select an uncommitted local file. Check
`getCurrentRaMPVersion(db=...)` returns the candidate label. Run representative
updated-vignette and frontend-backed R functions for pathway lookup,
enrichment, reaction, chemical class, chemistry, entity counts, and UpSet data.
Discover the current updated-vignette branch rather than assuming `main` has it.
Record the package, vignette, and frontend revisions used for the test.

The October 2026 owner review found that PFOCR accounted for 2,921 of 3,295
sample pathway rows and 288 of 326 enriched rows for its vignette examples.
For each candidate, check the unfiltered result and the existing enrichment
source-exclusion behavior, then record the product owner's decision about
source visibility, grouping, filtering, and any follow-up UI/package work.
The current frontend applies `dataSourceExclusion` **after**
`runEnrichPathways()`, so the displayed rows are filtered but adjusted
p-values are not recomputed for the reduced source set. Include that behavior
in the decision; a true source-restricted analysis would require an R-package
change and fresh testing.
Do not silently remove PFOCR data during packaging. The same review found six
of eight chemical-class vignette example IDs unmatched in both old and new
builds; refresh or explain that example before publication. These are explicit
review items, not reasons to mislabel an otherwise valid local candidate.

## 4. Handoff and release gate

The owner may commit the exact packaged bytes through Git LFS after reviewing
the report and checksums. Do not relabel or mutate the SQLite after that commit;
merge the same LFS object into `main` to release it. Before merging, decide the PFOCR presentation policy,
update the vignette example, and test the package's normal version discovery
from the intended RaMP-DB branch. The frontend currently initializes
`RaMP(branch='ramp3.0')`; update and test its branch/database selection as part
of deployment. Check that the displayed version is the version of the loaded
SQLite. A local file in `db/` alone does not complete those release steps.
