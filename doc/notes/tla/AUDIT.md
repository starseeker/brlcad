# BObol TLA+ suite audit

Audit date: 2026-09-03

This is the active formal-suite audit.  Current model-checking results and
state counts live only in `baselines/tlc-2.19.json` and generated logs.  This
file records semantic conclusions, not a second copy of verification output.

## Audit result

The structural and semantic findings from the initial reorganization have
been resolved in the suite and its gate:

- Fairness is assigned only to named implementation progress actions.  Input,
  policy, camera, visibility, cancellation, and other environment actions are
  premises, not promises.  `FAIRNESS.md` records the witness classes and the
  runner rejects broad fairness on `Next` or fairness on cataloged environment
  actions.
- Every configuration which disables TLC's built-in deadlock report defines
  and checks `DeadlockOnlyAtTerminal`.  The catalog records the terminal
  predicate, and lint requires complete agreement between catalog, module, and
  configuration.
- Required high-risk actions are cataloged.  The runner reads TLC's action
  coverage report and rejects a required action which is absent or has zero
  coverage.
- Suite-facing outcome contracts use `Ready`, `Constrained`, `Failed`, and
  `Terminal`.  Focused protocols which model ownership or transaction
  completion without classifying a product outcome use `Terminal` alone.
  `TERMINOLOGY.md` qualifies every local `phase` boundary and the overloaded
  revision/profile terms.
- The generated baseline is the sole source for current state counts.  Prose
  documents may describe historical runs but do not restate the active
  baseline.
- Every composition has machine-checked assumption/guarantee identifiers.
  Lint verifies that assumptions belong to cataloged components, guarantees
  belong to the composition, and all named operators exist.

## Formula changes from the semantic review

The audit found formula problems which a directory-only reorganization would
not have exposed:

- Arbitration's old bound made its `Complete` transition unreachable.  The
  configured budget domain now includes the maximum modeled demand.
- Convergence's empty-submission retirement branch was unreachable because
  frame completion suppressed the pending witness.  Empty retirement is now
  explicit, and an open input stream is represented as an external progress
  witness rather than assumed to close.
- Occurrence-control retry was a synthetic state change without corresponding
  in-flight work.  Retryable, successful, and terminal-failure completions are
  now typed actions; bounded repair acquisition owns real work and does not
  assume the final provider attempt succeeds.
- The canonical evidence stamp omitted visibility identity.  Visibility is now
  an independent revision domain and publication is authenticated against all
  six modeled evidence domains.
- The concrete control refinement could classify only ready or failed; a
  safely constrained implementation state would have projected to ready once
  work retired.  Constraint evidence is now independent of failure evidence,
  mutually exclusive, and part of the terminal refinement contract.
- Interaction closure could strand an open gesture when fairness was correctly
  removed from gesture completion.  Closing the input source now retires that
  gesture into the finite handoff protocol.
- Terminal-quality action-level fairness initially referred to partial quality
  actions whose structural counters were constrained only by `Next`.  A shared
  wrapper now completes those actions before both composition and fairness, so
  the promised action is exactly the modeled production transition.
- Asset publication previously used one constrained state for both resource
  limits and unrecoverable provider/geometry failure.  The composition now
  carries exclusive constraint and failure evidence and exposes distinct
  `Constrained` and `Failed` outcomes.
- Durable cache publication previously used a Boolean completion marker and
  could not represent replacement of a valid predecessor.  The focused live-
  publication model and its asset composition now carry a baseline mapping,
  atomic candidate mapping, and denied-commit state; denial is required action
  coverage and must preserve the baseline mapping.
- Host work previously modeled only an orderly producer close whose remaining
  notifications, timers, and frames could drain.  Endpoint loss removes that
  fairness witness.  `CloseEndpoint` now atomically cancels every host/source/
  provider/exact-publication level and exports its terminal postcondition to
  the lifecycle composition contract.
- The unbounded natural-number revision domains could not express the C++
  machine boundary.  `ObolIdentityExhaustion` now uses a finite domain in
  which a requested mutation either publishes a fresh successor or enters a
  closed fail-stop state without committing or reissuing an identity.
- A cross-repository refinement audit found mutable-state caches which treated
  compact hashes as exact evidence.  The formulas already compared their
  abstract records and identities exactly, so no formula change was required:
  Obol GPU keys and BRL-CAD capacity, structural-frontier, compact-source, and
  cross-source presentation credentials now retain exact inputs or receive a
  non-reused token only after exact comparison.
- Aggregate renderer evidence no longer hashes or sums local serials, or uses
  reusable object addresses as source identity.  Obol supplies a non-reused
  assembly identity; exact canonical source tuples produce checked draw,
  timing, resource, and preparation change tokens.  The formulas already use
  exact abstract identities, so this too is a refinement correction rather
  than a formula change.

Three focused models were added because the historical failure families imply
broader risks than the individual incidents:

- `ObolResultAuthentication` independently checks population epoch, demand
  revision, and source-route revision.  It now also models publish,
  exact-current terminal-failure, and retry dispositions, and requires stale
  results to create neither publication nor failure evidence before successor
  work is owned.
- `ObolSharedAssetLease` checks that one consumer cannot cancel a coalesced
  producer while another live consumer still depends on it, including a
  bounded late consumer joining during build, result, or post-cancel states.
- `ObolIdentityExhaustion` checks the fixed-width edge omitted by the natural-
  number control models: no credential is reused, stale evidence cannot
  authenticate, and exhaustion resolves through fail-stop quiescence.

`RISK_COVERAGE.md` maps these and the remaining prospective hazards to formal,
executable, and operational evidence.  `../libbobol_tla_conformance.md` records
the corresponding production audit and closure evidence.

## Deliberate limits

The suite is a control-plane roadmap, not a proof of unbounded C++ behavior.
TLA+ is used for ownership, authentication, publication, cancellation,
bounded progress, and terminal-evidence protocols.  Numeric geometry quality,
allocator behavior, renderer correctness, thread/runtime integration, and
performance still require executable fault injection, image comparisons,
sanitizers, resource measurements, and production traces.

The result-authentication implementation and its exhaustive asynchronous
matrix now refine the strengthened formula.  Explicit generation leases and
six build/result/cancellation lifecycle cases likewise refine
`ObolSharedAssetLease`; secondary consumers receive demand-local replay rather
than another view's payload.  Complete controller transition records now
carry typed endpoints, event/owner identity, and renderer-preparation rank;
their first full graphical trace found and closed the already-modeled
render-in-flight witness seam.  Complete admission stamps are now immutable
six-argument values, with per-domain stale-plan and retained-allocation
checks.  Commit-boundary denial now has executable retained-scene,
presentation, and durable-replacement witnesses plus a prior-mapping formula.
Endpoint teardown now has executable close-during-gesture, claimed-frame,
worker cancellation, and reentrant callback-lifetime witnesses plus an atomic
host-work cancellation formula.  Checked C++ successor primitives and the
finite exhaustion model now close fixed-width revision wrap as an
authentication path.  Remaining data-plane and operational limits are mapped
in `RISK_COVERAGE.md`, so they cannot be mistaken for properties proved by
TLC.
