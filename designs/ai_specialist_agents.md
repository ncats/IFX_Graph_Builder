# AI Specialist Agents

## Purpose

IFX_ODIN uses a small set of project-scoped AI specialists to bring deliberate
perspectives into feature work. The primary agent still owns the plan,
implementation, validation, and final answer. Specialists inspect and advise;
they do not make competing edits.

## Specialists

| Specialist | When it participates | What it contributes |
| --- | --- | --- |
| UX designer | Before user-interface or user-workflow changes | User journey, interaction model, states, accessibility, and acceptance criteria |
| Software architect | Before material changes to ownership, public contracts, persistence, dependency direction, or architectural patterns | Design shape, boundaries, reuse decisions, and tradeoffs |
| Code reviewer | After feature implementation and tests | Concrete findings about correctness, security, concurrency, regressions, maintainability, and test gaps |
| RaMP product reviewer | Manually, after a validated harmonization stage has been exported to SQLite | R-package vignette and frontend use-case comparison against released 2.4.0 and 3.0.7, with a written owner review and evidence figures |

The software architect is selective: routine work that follows an established
design does not require consultation. It also treats composition roots as
wiring boundaries, not
as homes for parsing, credential translation, domain adaptation, or migration
fallback policy. Any temporary compatibility layer should identify its current
consumer, owner, deletion gate, and target release; when a dependency can be
made mandatory, a decisive cutover is preferred over open-ended dual routing.

The architect and code reviewer remain separate to preserve an independent
post-implementation check. The reviewer evaluates the code that was actually
built and should not be anchored to, or repeat, the architect's proposed shape.

The agent definitions live in `.codex/agents/`. The implementation
consultants work read-only so advice remains separate from implementation.
The manually invoked RaMP product reviewer may write only its review report
and supporting artifacts in `output_files/ramp/`.

## Collaboration Flow

1. The primary agent determines which pre-implementation perspectives apply.
2. UX and architecture consultations may run together because neither edits
   files.
3. An implementation-stage UX review renders the affected pages with
   representative data at desktop and mobile widths. If rendering is blocked,
   the review is explicitly incomplete.
4. The primary agent summarizes any material choice for the user. Specialists
   provide a recommended default so routine work does not stall on minor
   questions.
5. After agreement where needed, the primary agent implements and validates the
   feature.
6. The code reviewer examines the completed change. The primary agent addresses
   substantive findings and reports the outcome.

This is selective rather than ceremonial. Documentation-only edits and trivial
mechanical changes do not need a specialist pass. A user can also request or
decline a specialist explicitly for any task.

The RaMP product reviewer is a manual release-assessment role, not an automatic
feature gate. It reads the completed SQLite and released 2.4.0 and 3.0.7
baselines through the RaMP R package. It tests the vignette and functions
called by `../ifx-frontend-library`, using the frontend's example arguments
and visualizations as review guidance, and documents both concerning changes
and outputs that work well. It writes its review to
`output_files/ramp/ramp_product_owner_review.md`, and places supporting
figures and rendered vignettes under `output_files/ramp/ramp_product_review/`.
It uses the unreleased 3.0.12 database only as clearly labeled diagnostic
context where a historical policy change needs explanation, never as a
release baseline. In particular, the RefMet review distinguishes intended
removal of isolated catalog entries from loss of linked or otherwise useful
compounds.
It also checks the chemical-identity concerns raised in
[`RaMP-DB#74`](https://github.com/ncats/RaMP-DB/issues/74), distinguishing
intentional broad pathway-mapping groups from implausible chemistry merges and
testing source-specific chemical-property lookup.
The updated vignette source is tracked on the RaMP-DB `gh-pages` branch rather
than `main`; the reviewer records that branch revision and distinguishes it
from the local execution copy under `../RaMP-DB/notebooks/`.

## Evolving the Agents

These files are team-owned working agreements. Adjust their preferences when a
review repeatedly produces low-value feedback, a missed concern recurs, or the
team develops a clearer design principle. Development feedback from the user is
an important source for those improvements: evaluate whether each comment
reveals a durable specialist preference or workflow. Put durable lessons into
the relevant agent instructions, and keep feature-specific decisions in that
feature's design so the specialist does not become narrowly overfit.

Add another specialist only when it has a distinct decision lens and a
repeatable workflow; avoid creating roles that merely restate the primary
agent's job.
