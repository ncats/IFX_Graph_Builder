# Metabolite curation first-release compilation

Status: completed on 2026-10-01. The first-release manifests are active in S3;
all provisional batch objects were retained for rollback.

## Purpose

The v2 curation batches were developed before their first public release. The
provisional S3 history can therefore be compiled once into a clean active policy
before immutable batch history becomes a compatibility contract.

After this cutover, published batches and manifests are immutable. Later changes
must be additive batches resolved through the normal latest-decision rules.

## Release contents

An ephemeral compiler should read all six metabolite curation streams and write
a portable local tree of batch JSON, manifests, and a migration report. It must
not write to or delete from object storage. Because this is a one-time migration
helper, run it from a temporary location rather than adding it to the repository.

Effective policy streams retain only active decisions:

- equivalence removals;
- record suppressions;
- active expected-clique assertions;
- active record-property overrides, one resolved decision per target/path.

MW adjudications are different. Runtime review deliberately reads their ordered
history so a clique can split into independently reviewed descendants. The
compiler therefore preserves every accept/reopen operation in chronological
order. A reopen command is policy-bearing and must not be discarded.

## Legacy RaMP deny-list attribution

Surviving removals from
`config/curation_mapping_issues_list.txt` are matched to line-level Git blame and
grouped into one batch per responsible commit. Duplicate pairs in the legacy
file are attributed to their newest blamed occurrence.

Valid historical Git authors map the two older batches to John Braisted. Six
later commits contain corrupt raw identity `= <=>`; their GitHub pull-request
associations provide explicit attribution evidence for Keith Kelleher:

- PR 15: commits `4cb9f30` and `6f574e7`;
- PR 18: commit `b858051`;
- PR 20: commits `db071ae`, `63a62d3`, and `c398d8b`.

Each compiled legacy batch keeps:

- original commit SHA and author date;
- raw Git author identity;
- attributed curator and attribution method/evidence;
- source repository and file path;
- separate migration publication time.

Git attribution is shown as attribution evidence in the review UI, rather than
being presented as proof of scientific authorship.

## Verification and cutover

Run the ephemeral read-only compiler only after any in-flight pipeline has
finished. Keep its code and working files outside the repository.

For every stream the command reads the generated objects back through the
production curation resolver and compares the semantic subject-to-decision
policy with the current snapshot. The migration report records old manifest
hashes and batch IDs, new hashes and IDs, active counts, and semantic policy
hashes. Compilation must fail when a configured legacy batch ID is absent or no
active operation is attributed from it. Batch IDs must include every immutable
payload input, including the fixed release publication time, so repeated runs
cannot reuse an object key for different bytes.

The compiler intentionally contains no S3 activation or deletion command. A
future cutover should:

1. Review the local curator assignments and migration report.
2. Upload all new immutable batch objects without replacing existing objects.
3. Read back and verify every uploaded hash and contract.
4. Conditionally replace manifests last, using the old manifest hashes as the
   precondition.
5. Re-run affected pipelines once because batch provenance changes their
   curation fingerprints even when effective cliques are identical.
6. Retain provisional objects until the release is validated; delete them only
   as a separate, explicitly approved cleanup.

The same local JSON tree can later be committed to GitHub without changing its
logical manifest/object-key layout.

## Completed release

Release `metabolite-curations-first-release-2026-10-01` was compiled at
`2026-10-01T17:42:01.770200Z`. Before activation, the generated tree was read
through the production resolver and compared with the prior effective policy
for all six streams. Batch and manifest hashes were then verified locally,
immutable batches were uploaded with create-only preconditions, and manifests
were activated using their previously read ETags. An independent S3 readback
matched all reviewed hashes and semantic policy hashes.

The active manifests now reference 15 first-release batches instead of 38
provisional batches:

- 10 equivalence-denial batches: eight Git-attributed legacy batches, one
  curator-confirmed batch of 63 non-Git decisions recovered from the September
  pre-release compaction, and one batch of 144 later QA decisions;
- one batch each for ChEBI properties, metabolite properties, expected cliques,
  record suppressions, and MW adjudication history.

The release preserves 1,404 effective equivalence removals, 16 ChEBI property
decisions, 98 metabolite property decisions, two expected-clique assertions,
five record suppressions, and the ordered eight-operation MW review history.
The 63 curator-confirmed decisions were absent from the working copy, `HEAD`,
and every historical committed revision of the legacy deny-list file; they
were therefore recorded as pre-release S3 curation decisions rather than
Git-attributed legacy rows.

No old batch object was deleted. Existing pipelines need to be synchronized or
rerun because immutable batch provenance changed their curation fingerprints,
although the effective curation policy is unchanged.
