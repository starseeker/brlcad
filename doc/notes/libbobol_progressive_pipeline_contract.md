# libBObol progressive-presentation contract

Last reviewed: 2026-08-26

This document is the canonical high-level control contract for progressive CAD
presentation.  It exists to keep the implementation understandable while it
serves inputs ranging from one small mesh to hundreds of thousands of distinct
assets and multiple simultaneously visible large meshes.  Workload names such
as Lucy, Hubble, multi-Lucy, and 150k are qualification profiles.  They are not
runtime modes and must not acquire special control paths.

The contract intentionally says less than the implementation.  Numeric screen
error, importance, memory, and frame-cost calculations are pure algorithms
below this boundary.  Mesh construction and rendering are data-plane work.
Neither may create another control owner.

## Canonical pipeline

```text
semantic inventory
        |
        v
immutable asset availability ----+
        |                         |
        +----> view demand -------+----> pure admission plan
                                             |
                                             v
                                  presentation transaction
                                             |
                                             v
                                  completed-frame evidence
                                             |
                                             +----> next plan, if facts changed
```

The implementation has five ledgers and one derived outcome:

1. **Inventory** identifies semantic occurrences and their immutable-asset
   keys for one source population revision.  Discovery is its sole writer.
2. **Availability** records the best complete immutable representation known
   for each asset, plus a named pending or terminal source result.  Bounded
   realization services are its sole writer.
3. **Demand** is a snapshot of occurrence visibility, projection, style, and
   interaction requirements for one view revision.  View evaluation is its
   sole writer.
4. **Admission plan** is a pure, deterministic function of an exact inventory,
   availability, demand, policy, and capacity revision tuple.  It owns all
   representation and quality allocation decisions.
5. **Presentation transaction** applies one plan to a retained candidate frame.
   A complete frame commits atomically; an interrupted or failed frame commits
   nothing and records typed capacity or error evidence.
6. **Outcome** is derived from the ledgers and transaction.  It is never a
   separately mutable phase machine.

No other component may decide that a scene is ready, choose a cut, turn a box
into a mesh, or reopen completed work.

## Revision and evidence rule

Every plan is certified by one revision tuple:

```text
(inventory, availability, view, policy, capacity)
```

A planner invocation records the tuple it consumed.  The same tuple cannot
produce another planning obligation.  A new plan is permitted only when one
component of the tuple changes or when a presentation transaction returns new
typed capacity evidence.  Diagnostic strings, HUD changes, timers, and render
reasons are not evidence.

This rule prevents the observed balancing/refining cycles: a controller cannot
reconsider identical facts merely because a timer fired or a status label was
redrawn.  Numeric planner output must also be deterministic for a tuple.

Planner identity is semantic rather than a bytewise snapshot of caller state.
An input guarded by a disabled feature is absent from the identity key; changing
its dormant raw value cannot invalidate a certificate or reopen planning.  One
canonical key type is used by plan caches, in-flight transaction matching, and
completed certificates.  Diagnostic and display-only fields never enter it.

## Outcomes and progress witnesses

The derived visible-presentation outcome is exactly one of:

- **active**: unfinished work has a named owner and an enabled progress event;
- **ready**: every visible occurrence has the richest pixel-justified admitted
  representation, with no failures or unacknowledged frame;
- **constrained**: all work is terminal, but a typed capacity or policy proof
  prevents richer pixel-justified presentation; or
- **error**: terminal source or presentation failure prevents truthful useful
  coverage.

`terminal` means no enabled foreground presentation work remains.  It is
independent of success.  `viewReady` means
`terminal && outcome in {ready, constrained}`; callers inspect the outcome to
distinguish exact pixel-justified readiness from a proven resource limit.  An
error must never masquerade as readiness.  Background cache construction and
persistence are reported separately and may continue after a terminal visible
presentation.  A successful terminal scene contains no structural boxes.

The retained HUD is an observer of this snapshot and may never retain a
stronger claim across a state transition.  `View ready` is permitted only when
`hasLodState && terminal && viewReady && !backgroundPending && fraction == 1`.
Input handlers publish the inexpensive label/bar transition synchronously when
they begin or end a view interaction; geometry rendering remains independently
queued and coalescible.  A stale terminal label beside an active controller is
therefore a contract violation, even if the following paint would correct it.

Provider-memory evidence constrains availability, not view demand.  For its
exact resident-admission revision, a denial prevents an allocation from naming
an unavailable suffix, but it does not invalidate or block a same-or-coarser
resident presentation.  The allocator retains the unconstrained pixel-demand
cost in its certificate and includes the admission revision in semantic plan
identity, so a capacity-release edge reopens planning without a special
workload path.

A constrained scene may use a truthful aggregate or certified coarse mesh, but
not an unresolved structural fallback.

Every active state identifies one witness from this finite set:

- discovery can publish another inventory result;
- realization can start or complete bounded asset work;
- demand can finish evaluating a new view;
- planning can consume a previously unplanned revision tuple;
- publication can adopt an immutable result;
- presentation can complete or interrupt the current transaction; or
- the host can service a level-triggered owner-thread request.

If no foreground witness is enabled, the visible state is terminal or
defective.  A UI polling timer is not a progress witness.  Background cache
work has the same owner and resource rules, but does not make an already
complete visible frame unready.

A witness is valid only when it advances a finite measure for one exact
revision-bound target.  Renderer-side preparation therefore reports a typed
obligation containing target identity, total and remaining immutable units,
admitted transient bytes, and one of `preparing`, `complete`, `constrained`,
or `failed`.  A retry is justified only when remaining work for the same
target strictly decreased.  An activity counter, elapsed time, a changed GPU
buffer serial, or work completed for a superseded target is diagnostic data,
not a liveness proof.  `ObolPresentationPreparation.tla` models this boundary.

## Workload-independent behavior

All inputs use the same pipeline.  Cardinality, asset size, reuse, visibility,
renderer speed, cache state, and memory pressure alter data and numeric
admission, not state ownership.

- A small affordable scene may have rich availability immediately and receive
  a direct-mesh plan.  This is an availability optimization, not a completion
  shortcut.
- Many small or subpixel occurrences are admitted to bounded aggregates when
  their individual cost is not justified.  They retain semantic identity.
- One large mesh uses spatial PoP pages when a single global prefix cannot
  localize visible demand.
- Several distinct large meshes share bounded byte scheduling fairly; one
  asset cannot monopolize preparation or admission indefinitely.
- Repeated large instances prepare an immutable asset once.  Each occurrence
  still has independent visibility, transform, selection, and demand.
- Cold and warm paths differ only in how availability is produced.  They feed
  the same planner and presentation transaction.
- System GL and OSMesa execute the same plan.  Backend capacity evidence may
  differ, but semantics and terminal classification may not.

An overview or structural box is an early truthful presentation.  It must not
delay immediately available affordable mesh data, and it cannot reappear over
a richer committed representation in the same view revision.

## Resource and execution contracts

Asset preparation is scheduled by admitted transient bytes, not worker count
alone.  Each active preparation owns a finite reservation; completion,
cancellation, and failure release it.  Large-asset fairness and shared-asset
deduplication are required independently of the number of workers.

Steady presentation work and one-time presentation preparation are separate
capacity dimensions.  A completed frame certifies steady draw cost.  Upload,
atlas construction, command packing, or software expansion is admitted by a
bounded preparation obligation and cannot be learned as steady FPS cost.
Preparation which exceeds its reservation returns a typed constraint; it
cannot acquire an unbounded retry merely because some internal allocation or
upload changed.  The representation executor may batch, page, or directly
reuse immutable data, but it must refine to this same obligation contract.

The admission planner owns retained-memory and frame-cost limits.  Producers
may report availability but cannot publish an unadmitted rich representation.
Non-visible resident suffixes may be compacted under pressure without changing
semantic inventory or durable cache identity.

Frame capacity is learned by one bounded search certificate for a frozen
population and revision tuple.  Each candidate is a complete scene-wide,
visual-importance-ordered allocation, not an independently promoted object.
A candidate is not eligible for measurement until every progressive
occurrence stores the cut assigned by that exact allocation certificate and
any renderer-wide ceiling is either removed or proven inert at or above every
active cut.  A pass which changes those occurrence cuts may complete this
handoff; requiring a mechanically unchanged successor would strand the search
behind the obsolete ceiling.  Equal aggregate cost is not an identity proof
for a multi-occurrence allocation because a different visual distribution may
have the same total.  The retained allocator therefore publishes a semantic-
epoch-scoped signature of every occurrence cut and point representation.  A
different numeric budget selecting that same signature reuses the population's
completed measurements immediately; it does not repaint an unchanged view for
another three-frame sample series.  An applied allocation hidden by an effective ceiling may
arm the candidate's sample-frame latch, but only its exact ceiling-free
successor may consume a sample.  Invalid or differently constrained frames do
not consume the candidate's bounded sample allowance.
A completed measurement either consumes one of a fixed number of samples for
that candidate or strictly narrows a known-safe/known-unsafe bracket.  A
terminal search publishes one capacity certificate and cannot reopen without a
new population or revision.  Throughput smoothing may propose a candidate; an
EMA change, timer, or repaint is not evidence and cannot itself advance the
capacity revision.  Likewise, completing a mechanical allocation scan without
selecting a successor population preserves the current capacity revision.  The
revision advances only when the decision names a new population; otherwise the
scan would invalidate the allocation certificate required by its own handoff
and could reopen the same allocation indefinitely.  Search work is bounded by
the candidate bracket, not scene cardinality or elapsed wall time.
`ObolCapacitySearch.tla` models this owner;
`ObolCapacityPresentationHandoff.tla` models the allocation-application,
ceiling-reconciliation, and exact-sample ordering boundary.

At the boundary between a completed mechanical pass and its successor,
capacity calibration and presentation handoff are candidate owners, not
independent cleanup stages.  `ObolCompletedPassOwnership.tla` requires that the
owner be selected before either path mutates evidence.  Resident growth and
coverage precede both.  A clean completed pass with an active presentation or
allocation handoff belongs to that handoff; generic capacity calibration may
run only when no stronger owner owns the boundary.  A current allocation may
also consume its cut/rescan annotations when its already-pending capacity
sample requires a ceiling-free presentation.  That normalization is part of
the immutable owner selection, not a later cancellation of capacity work.
Every effect-producing successor starts with fresh pass annotations, and a
same-population capacity sample does not advance the semantic capacity
revision.  Deactivating a capacity submission after calibration is
insufficient because calibration also installs a retained-reallocation
request.

The model follows growth, coverage, presentation, allocation, cut
presentation, population change, and same-population sample effects across
successive passes.  On 2026-08-27 TLC explored 2,270 generated / 1,804 distinct
states to depth 17 with no invariant or liveness error.  The executable mapping
is `completedPassSelection`, the availability-scheduler successor values, and
the pass-annotation/revision contracts.  The focused coordinator test exhausts
the selector's Boolean input domain and runs 512 deterministic composed traces.

An interrupted presentation has exactly one successor owner.  Strict progress
on the same finite presentation transaction takes priority; otherwise an
unfinished population cursor continues unchanged; only when both are quiescent
may capacity recovery change presentation policy.  A deadline callback cannot
reset a population cursor and request capacity recovery in the same event.
`ObolDeadlineOwnership.tla` models this boundary.

Renderer resource replacement is transactional.  New arrays and all required
companions are validated before retained metadata is committed.  Fixed-function
array state is exact for each command.  An interrupted frame does not advance
the completed-frame revision, update capacity calibration, or release the old
coherent framebuffer.

A presentation certificate contains the complete renderer control vector, not
an approximation of it.  For progressive CAD that vector is at least target
pixel error, integer PoP ceiling, fractional next-cut population, point-proxy
threshold, population identity, view revision, and measured presentation cost.
Saving an integer cut while dropping its fractional page population is a
quality regression, not a valid restore.  Pose, quiet-handoff, and exact-view
history paths must apply this vector atomically.

HUD and faceplate drawing observe an immutable progress snapshot.  Updating
diagnostic presentation may request a presentation-only redraw, but it must not
change a LoD revision, satisfy or create a LoD obligation, or overwrite the
typed reason for a pending transaction.

## Complexity budget

New state is more expensive than new code.  A change that adds a boolean latch,
timer retry, alternate phase, or workload branch must first show why it cannot
be expressed as a ledger value, revision edge, typed obligation, or pure
planner input.  The default answer is to remove or combine ownership.

Every pending state must document:

1. its sole owner;
2. the exact revision tuple it belongs to;
3. the event that enables progress;
4. the edge that retires or supersedes it;
5. its bounded resource reservation; and
6. the formal and executable test which exercises it.

Every policy decision has one owner.  Status and HUD code are pure observers.
Renderer diagnostics are not policy inputs.  Owner-thread work is bounded and
must not copy mesh arrays.  Planning must not be O(scene size times retry count)
or O(scene size times view-window count).

A proposed behavior change must update, as applicable:

- `ObolProgressivePipeline.tla` for ownership, safety, or liveness;
- the independent numeric allocation oracle for quality/resource changes;
- a focused executable state-transition test;
- the relevant graphical matrix profiles; and
- this contract if it changes an owner or invariant.

Image quality, elapsed time, and numeric error are not proven by TLA+.  They
remain release gates measured by images, APNG sequences, perf, memory telemetry,
and cold/warm System GL and OSMesa runs.

## Control grammar and implementation guardrails

The implementation may react to only these control-plane facts:

1. inventory advanced or closed;
2. immutable availability advanced, failed, or was compacted;
3. view or user policy changed;
4. a revision-stamped plan was produced;
5. a presentation transaction completed, was interrupted, or failed; or
6. an admitted resource reservation was released.

Timers and toolkit wakeups may deliver those facts, but are not facts
themselves.  Workload size, asset shape, cache warmth, and renderer backend are
planner inputs, never additional event kinds.

The control implementation follows a functional-core/imperative-shell split:

- ledgers own facts and monotonically advancing revisions;
- the admission planner consumes immutable snapshots and returns a decision
  plus complete successor evidence without mutating its input;
- a resumable execution cursor may track bounded progress through a plan, but
  it contains no quality, capacity, or terminal policy;
- the presentation transaction is the only pixel commit authority; and
- the controller performs effects selected by those values but may not mutate
  evidence fields piecemeal.

This yields reviewable implementation tests in addition to the behavioral
tests above:

- evidence components are private to their ledger/reducer;
- planning is deterministic and preserves unrelated evidence;
- every plan carries the exact revision tuple it consumed;
- reapplying an already planned tuple is an idempotent no-op;
- same-view commits cannot lose coverage or regress quality without a new
  typed capacity constraint;
- a cursor cannot survive superseding inventory, availability, view, or policy
  revisions; and
- every active outcome names exactly one enabled owner or a finite set of
  independent enabled owners.

These are architecture acceptance criteria, not aspirations.  Capacity and
admission evidence is now controller-read-only, and the bounded execution
cursor is a separate allocation-free value.  Plans and cursors now carry the
exact five-domain revision stamp.  Inventory, availability, and capacity edges
advance through one typed coordinator operation; view and policy advance at
their existing sole owners.  Each operation retires the preceding cursor, and
commit rejects a certified plan whose stamp no longer matches.  Thus no cursor
or plan can survive a superseding input, and an unchanged tuple cannot
manufacture another planning obligation.
The focused TLA+ model covers transition ownership and liveness; compile-time
encapsulation and executable transition tests connect that abstraction to C++.
Independent numeric oracles, fuzz/property tests, sanitizers, perf, and image
matrices cover the properties intentionally absent from TLA+.

Camera and input revisions change demand; they do not invalidate a complete
retained presentation.  A new plan starts from the committed representation
and quality.  It may improve that floor immediately, or retain it while richer
data arrives.  It may commit a coarser stable presentation only when the
transaction carries a typed memory, deadline, or explicit user-policy
constraint.  Merely entering a new view epoch, replaying a timer, or beginning
another calibration pass is not such a constraint.  Interactive renderer-only
limits remain reversible and do not rewrite this retained floor.

## Canonical implementation shape

The C++ control plane must make the small abstract model visible in its types.
The target is one value with four parts, not a collection of cooperating phase
machines:

```text
ProgressiveControlState
    EvidenceSnapshot
        InventorySnapshot
        AvailabilitySnapshot
        DemandSnapshot
        PolicySnapshot
        CapacitySnapshot
    optional PlanExecution
    optional PresentationTransaction
    WorkLedger
```

`EvidenceSnapshot` is immutable planner input and carries the exact five-domain
revision stamp.  `PlanExecution` contains only the bounded cursor and immutable
plan it is applying.  `PresentationTransaction` owns the one candidate frame
which may atomically replace the committed frame.  `WorkLedger` has a fixed
entry type for each legitimate producer or owner-thread obligation; an absent
entry means no such work exists.  Outcome, HUD phase, and progress are pure
projections of these four values.

The controller accepts only the finite event alphabet listed above.  Its
reducer returns a complete successor state and a bounded list of effects.  The
imperative shell performs those effects, reports typed results as later
events, and never edits control evidence directly.  A toolkit wakeup may ask
the shell to service an existing obligation, but cannot create an event or
advance a revision merely by firing.

State which remembers a proof is a typed, revision-bound certificate, not a
boolean plus companion fields.  Examples include a completed-frame capacity
certificate, a rejected quality trial, an exact-view quality-history entry,
and a resident-allocation certificate.  State which represents mutually
exclusive progress is a sum type, not several independently set latches.  In
particular, submission pass, retained allocation, publication handoff,
point-proxy trial, static-quality trial, and interaction handoff must each
have one typed owner and one explicit lifecycle.

The retained occurrence-allocation request is one concrete application of
this rule.  It is a finite value with `NONE`, preserve-budget,
recompute-budget, and presentation-reconciliation alternatives.  Only the
last alternative carries a reconciliation budget.  A weaker request cannot
overwrite that completed-frame proof, and resetting the request cannot leave
an orphaned budget.  Reintroducing separate pending, preserve, or budget
fields would recreate invalid states and is prohibited.

The following are not additional control states:

- direct affordable mesh, PoP prefix, spatial page, aggregate point, and
  structural preview are representation choices returned by the planner;
- few-small, many-small, single-large, multi-large, shared-large, cache-cold,
  cache-warm, hardware GL, and software GL are input profiles;
- interactive and quiet are demand/policy values; and
- performance-limited, memory-limited, and source-failed are typed terminal
  evidence used to derive an outcome.

This distinction is the primary complexity guardrail.  A proposed fix which
needs a new workload mode, retry phase, or free boolean is evidence that the
implementation boundary is wrong.  Work stops at that point until the behavior
can be expressed as a planner input, a typed certificate, or one of the finite
obligations above.

## Executable refinement map

Formal and production state must be connected explicitly.  Each production
transition is named by one abstract event and records, behind diagnostic
enablement, its old and new revision stamp, obligation kind, plan identity,
transaction identity, and typed outcome.  An offline trace checker and focused
randomized reducer test must enforce:

- no production transition exists outside the finite event alphabet;
- an unchanged evidence tuple cannot create a second plan;
- a stale plan cannot advance or commit;
- every active obligation has an owner and a strictly advancing progress
  witness expressed as a finite remaining-work measure for one exact target;
- a completed transaction is the only committed-frame mutation;
- same-view coverage and quality do not regress without a new typed capacity
  constraint; and
- after finite input, every fair trace reaches ready, constrained, or error.

The trace is observational and must be cheap or disabled in normal builds.  It
must never be consulted by production policy.  This creates a refinement
mapping from real executions to `ObolProgressivePipeline.tla` without trying to
model triangle arrays, floating-point projection, Qt scheduling, or renderer
internals.

The first production refinement boundary is now executable.  An
allocation-free map projects the remaining concrete controller facts onto the
finite interaction, inventory, availability, publication, planning,
presentation, handoff, compaction, and cache-write obligations, then selects
one owner by fixed precedence.  Convergence consumes this projection rather
than reconstructing separate Boolean unions, and qged records it for every
sample.  Exhaustive C++ coverage checks all 2,097,152 combinations of its 21
concrete inputs for owner/obligation consistency.  The focused
`ObolControlRefinement.tla` model symmetry-reduces structural-frontier and
point-recovery inputs because both map only to the planning obligation; its 20
fact classes are checked with both terminal-error values.  It proves that fair
owner service strictly decreases and eventually retires the finite work set
(4,194,302 generated / 2,097,152 distinct states).  This projection is not a
second scheduler and does not
satisfy the final reducer requirement by itself; each following migration
must move one production effect writer behind the typed boundary and delete
the superseded path.

The first runtime trace acceptance boundary is also executable.  The qged
diagnostic snapshot records the complete five-domain revision tuple, current
retained allocation plan, presentation transaction and required frame, and
committed-frame serial.  A standalone checker rejects sampled revision or
serial regression, more than one nonzero plan for an unchanged tuple,
spontaneous terminal reopening, observer disagreement, and A/B/A control
oscillation without an input or progress witness.  It intentionally accepts
an unchanged state while finite worker/cache work is active and recognizes an
explicit host command before its owner-thread revision synchronization.  The
checker is observational test code and cannot schedule work.  This closes the
sampled-trace gate, not the final requirement above to name and record every
production reducer transition with its finite remaining-work measure.

Numeric allocation remains a separate proof obligation.  Its pure snapshot-to-
plan function is checked against an independent oracle and property tests for
budget safety, deterministic ordering, monotonic error reduction, protected
visual-significance floors, and scene-cardinality independence.  Renderer
transactions retain boundary tests, sanitizers, APNG/image comparison, and
real GUI/perf qualification.  Passing one layer cannot substitute for another.

The refinement map is also the migration boundary.  New reducer variants do
not coexist indefinitely with old latches: each production caller moved behind
the reducer deletes its former write path in the same series.  If a migration
cannot identify the old state it supersedes, it is an additional framework and
does not proceed.

## Maintainability acceptance gate

A change to the progressive controller is reviewable only when all of the
following are true:

1. it names the abstract event and evidence domain it changes;
2. it does not add a workload/backend/cache-specific control branch;
3. new persistent state is a typed certificate or obligation with one owner,
   revision stamp, retirement edge, and resource bound;
4. the reducer remains deterministic and allocation-free on the owner-thread
   hot path except for the bounded plan storage it explicitly returns;
5. the same change deletes the superseded latch, recovery path, or policy
   owner;
6. focused reducer/oracle tests cover the new transition, and graphical tests
   cover its user-visible profile; and
7. the architecture document and TLA+ model change only if the actual contract
   changes, not to rationalize an implementation exception.

Line count alone is not a design metric, but ownership concentration is.  No
translation unit outside the reducer may write evidence or cycle state, and no
renderer/HUD/client callback may decide readiness, quality, or retry policy.
The long-term acceptance criterion is that a maintainer can enumerate every
control state and transition from the types and event reducer without reading
the renderer or workload-specific tests.

The current C++ is not yet at that acceptance point.  The latest 2026-08-25
mechanical audit still found 27 free top-level Boolean latches in the view
coordinator and 491 direct `this->d->lod*` assignments across the controller,
progressive pump, and render completion units.  The latter count includes
legitimate counters, cursors, and numeric presentation controls, so it is an
ownership-concentration indicator rather than a defect count.  Exact revision
stamps make stale work safe, but they
cannot make those combinations understandable.  Until the event reducer and
finite work ledger replace them, these counts are a complexity ratchet: a
change may reduce or encapsulate them, but may not add another free latch or
mutation owner.  A numeric estimator or immutable cache structure is not
control-state debt; preserving the measured data-plane algorithms is an
explicit requirement of this bounded control-plane redesign.

The renderer-preparation counterexample is now closed at the production
boundary.  Obol publishes an exact target, immutable total, monotone completed
work, retained scratch reservation, and typed terminal state.  libBObol
compares the before/after snapshots and grants an unchanged retry only for a
new target, strict same-target progress, or a completion edge.  The former raw
preparation serial is diagnostic-only and has no policy consumer.

The unit named by that certificate is a bounded cancellation atom, not an
arbitrary source-authored range.  PoP cluster ranges may be much larger than a
stable presentation deadline, so Obol splits flat-shaded atlas work into
triangle-aligned ranges no larger than the executor cancellation quantum.
Cancellation inside one range retains earlier committed ranges and leaves the
target `Preparing`; it is not malformed geometry and cannot publish `Failed`.
Only a complete range advances `completedUnits`.  This is the concrete
refinement of `ObolPresentationPreparation.PrepareUnit`.

Retained sparse-plan identity has the same single-owner rule.  A streamed
box-to-mesh replacement may leave a hidden compiled tombstone with the same
stable `InstanceId` as its appended active slot.  Reordering the tombstone
inside its old cut group may update an ID-to-index mapping only when that
mapping still names the moved source slot.  An inactive duplicate cannot steal
the active slot's mapping.  Every cut patch validates the group span, bin
cardinality, and active index before mutation; rejection rebuilds from the
authoritative instance records.  `ObolSparsePlanIdentity.tla` and the Obol
tombstone regression cover this aliasing boundary.

Three further boundaries are now explicit.  Interaction phase, gesture kind,
release timing, and quiet debounce are one typed interaction session rather
than independent latches.  Retained-allocation reuse and in-flight matching use
one canonical semantic input key, so disabled maximum-protected-budget values
cannot manufacture a new plan.  Deadline handling selects one typed successor:
transaction retry, population continuation, or presentation recovery.  The
focused C++ transition matrices and `ObolInteractionSession.tla` and
`ObolDeadlineOwnership.tla` cover those boundaries.

The quiet-restoration counterexample is closed at a narrow functional-core
boundary.  `BObolLodQuietSuccessorReducer` consumes the current presentation,
revision-matched prior-stable, proven-scale, and exact-view certificates and
returns one complete successor and one typed handoff obligation.  Certificate
precedence is fixed and independent of the last transient motion frame: exact
view, retained prior stable, proven scale, then stable one-pixel demand.  A
prior-pose restore requires one current-pose presentation proof before the
allocator may translate it.  Other motion-only renderer limits may survive
only as a bounded allocator handoff; they cannot choose the semantic quiet
target.  The controller no longer performs a later exact-history overwrite or
a separate stable-demand fallback.  The older mouse-release-time pose restore
was also deleted: release retains the responsive motion image through the
debounce, and only the quiet reducer may restore and certify the prior pose.
Thus a pre-quiet miss cannot convert a temporary motion budget into persistent
capacity evidence.

The focused C++ test permutes transient motion target, ceiling, fractional
population, and point-proxy values while holding the semantic certificates
fixed and requires the same successor.  `ObolQuietSuccessor.tla` independently
checks two differently ordered transient schedules, single successor
publication, and eventual terminal behavior.  This is the required pattern
for subsequent control reductions: add one pure owner and delete every older
writer in the same slice.

A 2026-08-25 experiment deliberately tested the tempting shortcut: every
certified restored presentation immediately entered the static-quality trial.
It improved Lucy's post-rotation cut, but a 50k warm OSMesa run failed to reach
its first stable scope within 180 seconds and a 150k cold run failed to settle
after rotation within the same bound.  The experiment was reverted.  This is
not evidence for a many-object exception; it proves that restoration and
authorization of new quality work are distinct events.  The interaction
reducer must first restore one complete certified target, then authorize at
most one bounded successor plan from current demand and capacity evidence.

The bounded capacity-search certificate now owns the retained pixel-demand
endpoint.  It consumes the quiet successor as immutable input and returns only
a revision-stamped bounded successor plan or typed constraint.  Timer expiry,
object count, cache temperature, backend, and workload name are not result
variants.  The remaining migration boundary is to move ordinary coverage
calibration behind this same certificate and delete its older probe/count/
rescan authority.

The exact-current 150k OSMesa trace also separates finite mechanical work from
policy convergence.  Its occurrence cursor advanced monotonically in 2,048-
entry windows, normally taking 5.7--7 ms per window and less than a second for
a continuously scheduled full census.  The long settle time came from repeated
completed-frame calibration and reallocation waves.  This is not justification
for a 150k branch or a larger cursor window.  Complete the certificate
migration above: every same-population successor must consume a finite sample,
narrow its bracket, or publish terminal evidence.

The follow-on exact-current traces exposed two composition failures which the
snapshot refinement map alone could not prevent.  First, capacity sampling and
generic presentation cleanup both selected successors for one completed frame.
Second, an ordinary resident-suffix drain fell through to capacity handling
while the availability ledger still (correctly) prohibited allocation.  The
implementation now gives a capacity sequence exclusive successor ownership and
represents resident growth as one typed ledger phase.  Drain completion can only
select another coalesced drain or one retained allocation.  A later trace found
the complementary pre-drain edge: an ordinary completed pass restarted capacity
calibration while resident growth waited for its cursor to become idle.  One
pure completed-pass successor now gives active drain completion precedence,
then yields to pending resident growth, and permits capacity calibration only
when availability has no claim on the population.

The 2026-08-26 50k replay found the same ownership error at the next boundary.
A finite point-to-triangle recovery allocation could cover the complete
occurrence population while deliberately leaving richer pixel-target demand.
The residual quality annotation was treated as unfinished recovery admission,
and ordinary capacity calibration repeatedly reopened the same all-scene pass.
Recovery now owns capacity after any resident-growth transaction retires.
Residual quality debt is terminal for that recovery only when the current
allocation certificate covers the current population; otherwise admission
continues.  Ordinary capacity search cannot start until the recovery has
either presented its changed cut or retired its proven no-op plan.

`ObolResidentGrowth.tla` is the focused liveness model for that transaction.
It also models the incomplete-coverage boundary exposed by an OSMesa Generic
Twin zoom: coverage cannot guard the resident drain which may make coverage
possible.  A completed drain consumes its growth edge and transfers to either
one ordinary coverage retry or one scene-wide reallocation.  Repeating the
same incomplete coverage pass while the growth phase remains pending is not a
legal transition.
It includes an ordinary capacity-blocked cursor at the growth edge and proves
that capacity cannot restart after resident growth becomes pending.  TLC
checked 525 generated / 377 distinct states to depth 23 with no error.
`ObolPointQualityOwnership.tla` additionally distinguishes adaptive point-cut
calibration, handoff confirmation of an already chosen cut, and retained
triangle recovery.  A recovery frame whose callback is consumed by a stronger
owner retains a level-triggered pump successor and eventually retires.  The
model also covers certified no-op/residual-debt recovery, a point request which
arrives while one finite capacity search drains, atomic rejection of a
one-pixel static trial back to its retained threshold, and replacement of an
adaptive frame by handoff confirmation.  Point recovery and ordinary capacity
calibration remain mutually exclusive.  TLC checked 29 generated / 17 distinct
states to depth 5 with no error.  These models establish transition composition
and liveness, not elapsed time or visual quality; the 50k OSMesa, 150k OSMesa,
and Lucy renderer rows remain the runtime authorities.

`ObolStructuralFrontierOwnership.tla` refines the adjacent renderer-owned
frontier which the point-quality model intentionally did not represent.  An
exact completed frame containing non-terminal structural fallbacks must own
one of two finite successors: another point-distribution frame while its
threshold search can still reduce the frontier, or bounded mesh repair after
that policy is exhausted.  An ownerless nonempty frontier is therefore an
enabled planning obligation, never a terminal or merely diagnostic state.
TLC checked all 30 distinct states to depth 14 for successor availability,
owner/frontier compatibility, no premature readiness, and eventual readiness.
The production abstraction mapping is
`BObolLodControl::Inputs::structuralFrontier`,
`pointSuccessorRequiredForStructuralFrontier`, and the exact structural
histogram branch in `completeRenderTiming()`.  The executable refinement is
the coordinator truth table plus the Lucy, 50k, and 150k OSMesa matrices.

Camera-frame reuse is owned by the typed interaction session through its
release debounce.  Button release is not a demand, capacity, or presentation
event and may not disable reuse or start exact view classification itself.
The quiet-successor transition performs that handoff after the finite debounce
and then schedules exact classification and refinement.  This preserves the
last complete responsive framebuffer while keeping the terminal exact-view
contract unchanged.

## Counterexample classification protocol

Every observed failure is classified before production code changes:

1. **Contract violation:** ownership, revision, progress, atomic commit, or
   terminality is wrong.  Change the finite event/reducer contract, update the
   relevant TLA+ model, and delete the superseded control path.
2. **Numeric policy defect:** the event graph is valid but allocation,
   projection, quality, or capacity arithmetic selects a bad plan.  Change a
   pure function and its independent oracle/property tests; do not add state.
3. **Data-plane defect:** immutable mesh construction, cache data, upload,
   renderer state, or picking is incorrect.  Fix it below the control
   boundary and add sanitizer/boundary/image evidence; do not teach the
   controller about the backend or shape.
4. **Observation defect:** HUD, diagnostic, test timing, or image capture
   reports the wrong state.  Fix the observer without advancing a control
   revision or manufacturing work.

A proposed fix which spans more than one class is split into independently
testable changes.  If the failure cannot be assigned to one class, the
ownership boundary is not yet understood and implementation pauses for a
design update.  This rule prevents a graphical counterexample from becoming
an ad hoc controller regime.

## Required implementation simplification

The implementation predates this single-pipeline contract.  Refactoring must
preserve performance-sensitive algorithms while reducing control structure.
Completed items remain contractual constraints and may not be reintroduced:

- **complete:** the independent imperative state machine and its diagnostic
  phase mirror were deleted; one derived progressive outcome snapshot is the
  public phase/outcome authority;
- **complete:** publication batching and frame acknowledgement are one typed,
  revision-bound presentation transaction with atomic completion;
- **complete functional-core boundary:** capacity, resource pressure, headroom, point-proxy, and
  structural-admission facts now share one allocation-free admission-evidence
  snapshot.  Their transition entry points copy that snapshot and return one
  typed decision plus complete successor evidence.  Static quality predicates
  are part of the same planner, and the focused test verifies determinism and
  preservation of unrelated evidence.  Evidence components are exposed to the
  controller only through const views, and the per-pass execution cursor is a
  physically separate value;
- **complete revision boundary:** every admission plan and bounded cursor
  carries `(inventory, availability, view, policy, capacity)`.  The roughly
  forty free-form production cursor resets were replaced by the five semantic
  revision owners, all five domains now use one allocation-free executable
  revision transition, external policy-epoch synchronization retires the
  cursor through that boundary, and stale certified commits are rejected;
- **complete for coordination bookkeeping:** provider terminality, result
  readiness/age, inventory-delta coalescing, and resident-growth obligation
  now have one availability ledger, while pacing and settled-population
  predicates have one stateless availability scheduler.  The former inventory
  delta, provider pacing, and resident growth policy owners were deleted.
  Completed-pass succession is also pure: availability drain/growth ownership
  precedes ordinary capacity calibration, so those schedulers cannot restart
  the shared cursor behind one another.
  Immutable per-asset availability remains owned by the LoD service/cache as
  required; it must not be copied into the controller ledger;
- **complete presentation-certificate boundary:** integer PoP ceiling and its
  fractional next-cut population are one value with pixel target, point-proxy
  threshold, population identity, view epoch, and measured cost.  Pose restore,
  quiet restore, and exact-view history preserve the complete vector; focused
  executable tests reject a partial restore;
- **complete semantic allocation identity:** certificate reuse, in-flight
  transaction matching, and allocator input comparison use one canonical key;
  guarded-off inputs are normalized out of that identity;
- **complete deadline-successor boundary:** a presentation deadline chooses
  exactly one of finite-transaction retry, monotone population continuation,
  or presentation recovery; it cannot reset a live producer cursor;
- **complete quiet-successor boundary:** one allocation-free reducer selects
  the complete first quiet target from revision-matched certificates; the
  former controller exact-history overwrite and fallback target writer were
  deleted;
- **complete point-quality ownership boundary:** quiet point calibration and
  retained triangle recovery are one mutually exclusive sum type.  Recovery
  owns source admission until its changed cut has a presentation witness; a
  calibration request cannot displace or coexist with it;
- **partially complete capacity boundary:** a revision-bound bounded search now
  owns the retained pixel-demand endpoint, including the one-way steady-to-
  static transition and a fixed sample/candidate bound.  Migrate the remaining
  ordinary probe/calibration branches to that certificate and delete their
  companion latches; an EMA may propose a candidate but may not own retry;
- **remaining:** retain visibility census and projected-demand caching as
  algorithms, not mutable policy authorities; and
- **complete HUD publication boundary:** the render host samples convergence
  independently and publishes only the three LoD progress records through
  `ged_view_lod_progress_sync`.  Periodic progress no longer invokes the
  general faceplate synchronizer or republishes axes, grid, ADC, framebuffer
  composition, and unrelated HUD state.  The focused faceplate test changes
  semantic grid state between calls and proves a progress-only publication
  leaves the retained grid node untouched; the subsequent full sync applies
  the change.

The old wording is retained here only as a review checklist; the implementation
must not recreate any completed owner under a new name.  In particular, there
is no workload-mode dispatcher and no mutable HUD phase.

The next simplification target is therefore:

- replace the remaining overlapping controller latches with the canonical
  `ProgressiveControlState`, typed certificates, finite work ledger, and one
  event reducer;
- extract the remaining controller shell by ledger and effect boundaries,
  without moving or duplicating performance-sensitive numeric algorithms;
- retain source visibility census and projected-demand caches as derived
  algorithms rather than new policy owners;
- preserve the progress-only HUD publication boundary; and
- leave numeric estimators as stateless or revision-keyed helpers; and
- preserve the observational sampled-trace checker while extending the
  canonical reducer to record every abstract production transition.

### Asset production and occurrence presentation

Repeated instances share immutable asset work, not view state.  At most one
producer may build or load a given asset generation, while every occurrence
retains its own transform, visibility, spatial-page demand, active cut, style,
selection, and presentation commit.  Waiting for a coalesced producer is not
pending owner-thread work; the producer result is the only wake edge.

A payload's cache key or resident cut does not certify its retained renderer
binding.  Presentation is complete only when the occurrence is bound to the
deterministic payload part and its prepared whole-mesh geometry or page layers
cover the exact current demand.  Provider results must not copy requested page
ids into `presentedChunks` without that renderer-preparation evidence.

Spatial visibility is view-local presentation state.  A large leaf may have no
page intersecting the current frustum even when its conservative whole-asset
box touches it.  Such an occurrence remains authored-visible and resident but
is culled from the current assembly; its box is neither drawn nor counted as a
convergence obligation.  A later camera epoch restores it through ordinary
page demand without rediscovery.

The temporary whole-target overview is not an LoD leaf.  It retires only after
the producer-certified leaf population has been installed and presented, but
leaf boxes may satisfy that coverage frontier.  Subsequent box-to-mesh repair
is an independent occurrence-presentation transaction.

This is a staged replacement, not a second framework.  Each new ledger or typed
transaction must delete the superseded latches and policy owner in the same
series.  Until replacement is complete, adapters may expose old private fields,
but they may not make decisions.

## Qualification profiles

The same executable contract must pass these independent axes:

| Profile | Occurrences | Unique assets | Largest asset | Reuse | Primary risk |
|---|---:|---:|---:|---:|---|
| Generic Twin | hundreds | hundreds | small | low | direct availability must not be delayed by boxes |
| Hubble | thousands | many | mixed | mixed | visual importance and hierarchy interaction |
| Lucy | one | one | large | none | spatial demand and close-zoom refinement |
| multi-Lucy | several | one | large | high | occurrence demand with one prepared asset |
| distinct multi-large | several | several | large | none | bounded preparation and fairness |
| 50k/150k varied | many | many | mixed | low | aggregate admission and owner-thread bounds |

Each profile is run cold and warm, shaded and wireframe where applicable, on
System GL and OSMesa.  Passing one row cannot authorize a different control
path for another row.

## Formal evidence

`ObolProgressivePipeline.tla` is the bounded model of this contract.  Its two
occurrences and assets are a symmetry reduction, not a cardinality assumption.
It covers few-small, many-small, single-large, distinct multi-large, and
shared-large profiles through finite input, view, and policy epochs.  Inventory,
availability, view, policy, and capacity have independent revisions.  A
bounded plan captures their exact tuple; if any component changes before the
plan finishes, the stale plan can only abort and cannot advance or commit.

On 2026-08-25 TLC explored 2,986,118 generated and 1,427,962 distinct states to
depth 42 with no invariant, deadlock, or liveness error.  It proves bounded
preparation ownership, shared-asset reuse, no direct meshes in large profiles,
atomic non-regressing same-view commits, evidence-gated planning, truthful
terminal outcomes, stale-plan rejection, active progress witnesses, eventual
terminal behavior after finite input, and non-starvation of distinct large
assets.  Its strengthened retained-floor invariant also proves that an input or
view epoch alone cannot authorize a coarser stable mesh commit.  It does not
prove renderer correctness, memory arithmetic, frame rate, or perceptual
quality.

`ObolPresentationPreparation.tla` isolates the renderer/host liveness seam.
On 2026-08-25 TLC explored 1,143 generated and 526 distinct states to depth 12
with no invariant or liveness error.  It proves that a retry requires strict
same-target reduction of a finite remaining-work measure, that active work has
an admitted reservation, that superseded targets cannot commit, and that
finite final input reaches complete, constrained, or failed.  It deliberately
does not estimate upload throughput or choose a representation; those remain
Obol data-plane tests and graphical/perf qualification.

`ObolInteractionSession.tla` isolates input ownership and passed 98 generated /
44 distinct states to depth 10.  `ObolDeadlineOwnership.tla` isolates deadline
successor priority and passed 43 generated / 34 distinct states to depth 5.
`ObolQuietSuccessor.tla` isolates schedule-independent motion-to-quiet
restoration and passed 1,242 generated / 924 distinct states to depth 5.  It
proves that equivalent semantic certificates produce one identical successor
despite different transient motion publication orders, that a prior-pose
restore cannot bypass its current-pose proof, and that finite input terminates.
`ObolCapacitySearch.tla` checks all modeled true capacities for eight ordered
budgets, a monotone map which pairs adjacent budgets onto the same discrete
population, and three bounded samples for every previously unseen population.
Equivalent populations reuse their prior classification and select a midpoint
of the remaining proven bracket; they do not walk adjacent numeric cost units
which cannot change the discrete PoP cuts.
It includes the one-way transition from the preferred steady cadence to the
independent static deadline.  TLC explored 2,918 generated / 2,918 distinct
states to depth 27.  It proves sound strict bracket reduction, immediate reuse
of a previously classified population, non-repetition of measured budgets, at
most one goal transition, single terminal certificate publication, and
eventual completion for a frozen tuple.

`ObolCapacityPresentationHandoff.tla` isolates the exact-presentation boundary
in front of that search.  On 2026-08-27 TLC explored 572 generated / 390
distinct states to depth 10 with no invariant or liveness error.  It covers
changed and already-applied occurrence plans, effective and inert global
ceilings, and an applied occurrence allocation whose protected population
exceeds its own certified budget.  It also models an initially unavailable
assignment and requires availability restriction before the cut can be
applied.  A completed no-op allocation scan preserves the frozen population
and allocation-certificate revisions.  It proves that pre-handoff frames
cannot consume a sample, an effective ceiling is removed before measurement,
and the over-budget state has one finite owner: local representation reduction
followed by either an exact sampled certificate or an explicit bounded
constraint.  The
executable mapping is `claimOverBudgetAllocation`,
`capacitySampleRequiresCeilingFreeHandoff`,
`cadAllocationPlanCutsApplied`, and `capacitySamplePopulationReady`; the
focused coordinator/update-action tests and the 50k/multi-Lucy renderer
replays are the refinement evidence.

A presentation-required handoff retains the recovery budget computed from the
interrupted richer frame even when its corrected renderer frame completes just
beyond the endpoint timer and therefore supplies no new safe extrapolation.
The model's over-budget handoff state represents this persistent evidence: it
may transition only through local reduction or bounded constraint publication,
never by reclassifying the unchanged population as fitting.  The executable
mapping is `BObolLodPresentationPolicy::noteFramePresented`; the focused
timer-edge test and Lucy OSMesa wire zoom lifecycle are the refinement evidence.

`ObolPointQualityOwnership.tla` refines the high-level quality obligation into
mutually exclusive adaptive calibration, handoff confirmation, and triangle-
recovery phases.  TLC explored 29 generated / 17 distinct states to depth 5.
It proves that each active phase has its required frame or producer witness,
an already-active capacity search drains before a waiting point frame,
confirmation cannot enter adaptive recovery logic, rejected one-pixel trials
restore their retained baseline, recovery cannot coexist with calibration,
and finite work reaches a quiescent terminal state.  The focused C++ tests
also check the production effects which schedule the point successor after a
completed allocation pass and the recovery frame after a changed recovery
pass.  These focused models do not prove timing classification, allocation
quality, or renderer performance.

The executable timing refinement must preserve the same distinction.  The
broad `presentationPending` predicate means either an adaptive census or a
confirmation frame is outstanding; only `adaptiveCalibrationPending` may
authorize the finite 400 ms point-census allowance.  Confirmation is an
already-selected presentation handoff and obeys the independently selected
quiet/static-quality deadline without receiving the census extension merely
because it is point work.
Software raster calls are subdivided into 8K-triangle cancellation units so
that the host can enforce that deadline without exposing a partial frame.
This data-plane bound is established by perf and graphical response tests,
not by the TLA+ model.

The formal argument is compositional: `ObolProgressivePipeline.tla` defines the
only high-level event and ledger grammar; each focused model refines one of its
typed obligations without adding a workload regime.  A focused model must name
its abstraction mapping and executable production effect.  It may not expand
the high-level event alphabet merely to describe an implementation latch.  A
counterexample which cannot map back to inventory, availability, demand,
planning, presentation, or resource release is architectural debt and stops
implementation until ownership is clarified.

Every formal action requires an executable production refinement test.  The
deadline model correctly transferred ownership from render replay to a live
population cursor, but the first C++ effect mapping left the old replay latch
set.  A 50k GUI run exposed the resulting closed wait; the focused controller
test now asserts that the producer is actually invoked after that transfer.
Passing TLC without this mapping evidence is not a control-plane acceptance
criterion.
