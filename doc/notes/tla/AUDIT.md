# BObol TLA+ suite audit

Audit date: 2026-09-02

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

Two focused models were added because the historical failure families imply
broader risks than the individual incidents:

- `ObolResultAuthentication` independently checks population epoch, demand
  revision, and source-route revision, and requires stale results to be
  rejected or converted into explicitly owned successor work.
- `ObolSharedAssetLease` checks that one consumer cannot cancel a coalesced
  producer while another live consumer still depends on it, including a
  bounded late consumer joining during build, result, or post-cancel states.

`RISK_COVERAGE.md` maps these and the remaining prospective hazards to formal,
executable, and operational evidence.

## Deliberate limits

The suite is a control-plane roadmap, not a proof of unbounded C++ behavior.
TLA+ is used for ownership, authentication, publication, cancellation,
bounded progress, and terminal-evidence protocols.  Numeric geometry quality,
allocator behavior, renderer correctness, thread/runtime integration, and
performance still require executable fault injection, image comparisons,
sanitizers, resource measurements, and production traces.

The remaining open work is therefore implementation evidence, not an
unrecorded formal assumption.  The highest-priority gaps are the shared-build
consumer-close regression, the combined stale-result identity matrix,
allocator failure at commit boundaries, endpoint loss during interaction,
and a C++ audit of all six evidence-stamp constructors.  They are recorded as
completion work in `../libbobol_active_debt.md` and mapped to their formal and
executable evidence in `RISK_COVERAGE.md`, so they cannot be mistaken for
properties already proved by TLC.
