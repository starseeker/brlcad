# BObol TLA+ terminology and local-state map

Production has no global BObol phase machine.  A model variable named `phase`
is a finite abstraction local to one protocol boundary.  In prose, reviews,
and refinement mappings, qualify it with the boundary below.  The matching
`models.json` `localPhaseMeanings` map is linted exhaustively so a phase-bearing
model cannot be added or renamed without an explicit qualification.

| Module | Qualified meaning of `phase` |
|---|---|
| `ObolInteractionSession` | interaction-session phase |
| `ObolQuietSuccessor` | quiet-successor handoff phase |
| `ObolSpatialProducer` | durable spatial-producer phase |
| `ObolCapacitySearch` | bounded capacity-search phase |
| `ObolCapacityPresentationHandoff` | capacity-result presentation phase |
| `ObolCompletedPassOwnership` | completed-pass disposition phase |
| `ObolPointQualityOwnership` | point-quality ownership phase |
| `ObolPointTerminalEvidence` | point terminal-evidence phase |
| `ObolResidentGrowth` | resident-growth/drain transaction phase |
| `ObolStaticQuality` | static-quality trial phase |
| `ObolSparsePlanIdentity` | sparse retained-plan patch phase |
| `ObolRetainedAllocationPresentation` | retained-allocation presentation phase |
| `ObolTerminalQualityOrdering` | terminal-quality policy phase |
| `ObolTerminalConvergenceComposition` | terminal-convergence work phase |

These values refine control locations and ownership facts, not a shared C++
enum.  A transition between two rows in this table is consequently never
implied merely because both models use a value such as `terminal`.

## Outcome predicates

At suite-facing outcome boundaries:

- `Ready` means the requested result is current and usable;
- `Constrained` means the protocol terminated safely with explicit evidence
  that the requested result could not be fully achieved;
- `Failed` means the protocol terminated with a modeled error outcome; and
- `Terminal` is exactly `Ready \/ Constrained \/ Failed`.

A focused ownership or transaction model may expose only `Terminal` when it
does not have enough facts to classify a product outcome.  It must not invent
`Ready` as a synonym for "no work is currently owned."  Cancellations and
supersessions are terminal for some local lifetimes but are not automatically
ready, constrained, or failed product outcomes.

## Qualified overloaded terms

- **Source profile certification** identifies whether source geometry has a
  usable certified representation.  **Workload profile** is a finite model
  parameter describing a demand shape.  They are unrelated types.
- **Population epoch** identifies the semantic set from which work was built.
  **Demand revision** identifies a requested asset/quality demand.
  **Source-route revision** identifies the provider/cache route authorized to
  produce that result; changing it invalidates an in-flight predecessor even
  when the population slot and demand value are reused.
  **Visibility revision**, **view revision**, **policy revision**, and
  **capacity revision** authenticate their respective domains.  Equality in
  one domain never substitutes for equality in another.
- **Presentation revision** identifies a candidate offered to a frame.
  **Frame acceptance** is the renderer-side acknowledgement for that exact
  presentation.  It is not durable-cache publication.
- **Retained-scene commit** atomically replaces the retained presentation.
  **Durable-cache commit** publishes an immutable reusable asset.  A failure
  before either commit preserves the corresponding prior object.

When a new model needs one of these concepts, reuse the qualified meaning and
type.  A shared TLA+ vocabulary module is warranted only after two models have
identical semantics and types; textual similarity alone is insufficient.
