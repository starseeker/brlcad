# libBObol formal model catalog

Last reviewed: 2026-09-01

The authoritative suite index, proof-flow graph, shared vocabulary, and
verification policy are in `tla/README.md`, `tla/models.json`,
`tla/GLOSSARY.md`, and `tla/verification.md`.  This document explains the
architectural intent and production mapping.  The models are small bounded
model-checking boundaries for ownership, safety, and liveness.  They
are not alternate implementations, workload-specific controllers, or proofs
of visual quality, floating-point projection, memory use, or elapsed time.
All state counts retained in this narrative are dated historical observations;
the sole current result source is `tla/baselines/tlc-2.19.json`.
Representation-only metadata such as cached PCA OBB corners therefore does
not add a state variable: runtime validation proves its geometric contract,
and renderer tests prove equivalent shaded/wire batching.  If proxy shape ever
becomes an independently scheduled allocation decision rather than a pure
projection of the admitted aggregate representation, it must first be added to
the canonical owner/terminal relation.

## Canonical refinement

`ObolProgressivePipeline.tla` is the canonical finite control contract.
`ObolControlRefinement.tla` maps the concrete C++ obligation facts onto that
contract.  `ObolLodComposition.tla` checks the cross-boundary seam from
admission through retained growth, exact presentation, ceiling reconciliation,
capacity sampling, structural repair, point quality, and terminal publication.
`ObolTerminalConvergenceComposition.tla` expands the closed-input tail which
that model formerly treated atomically: an exact visibility delta must finish
its finite census and exact framebuffer classification before reallocation,
distinct static populations may reopen bounded capacity searches, each exact
presentation consumes finite renderer preparation units, over-budget
allocations may consume successively coarser local representations, and
background compaction follows—but cannot reopen—the terminal foreground
population.
`ObolControlLifecycleComposition.tla` closes the orthogonal lifecycle seam
between policy enablement, provider work, view demand, host pumping, exact
frames, semantic-only selection/style revisions, deadline recovery, and
terminal status.  In particular, a semantic-only successor may require an
exact framebuffer without manufacturing an LoD-capacity obligation.
`ObolAssetPublicationComposition.tla` and `ObolCadFrameComposition.tla` close
the two adjacent data-publication seams: immutable producer output through
current-demand presentation/cache completion, and atomic retained mutation
through exact-target frame acceptance.
Every other model isolates one boundary of the same transition relation.  A
focused model may add detail, but it may not introduce another production
owner or control mode.

The production refinement target is the typed state in
`lod_control_private.h`, the evidence and planning policies in
`lod_coordinator_private.h`, and the imperative effect executors in the
`view_controller_*.cpp` units.  Workload labels such as Lucy, Hubble, 50k,
150k, cold, warm, System GL, and OSMesa are test profiles only.

## How to read the proof hierarchy

Focused models prove one building block's abstract contract: legal states,
exclusive ownership, monotonic evidence, finite progress measures, and valid
retirement edges.  Composition models then prove that adjacent building-block
contracts are compatible and that their combined control flow has neither an
ownerless obligation nor an unsupported terminal outcome.  The concrete
refinement model and C++ assertions/tests check that production fields and
effects implement those abstract facts.

This establishes necessary and sufficient control evidence only within the
modeled abstraction.  It does not prove that the abstraction includes every
user requirement, nor does it prove perceptual quality, geometry correctness,
memory consumption, latency, or frame rate.  Those data-plane and quantitative
properties remain the responsibility of executable property tests, renderer
tests, image comparisons, profiling, and the graphical qualification matrix.
Thus a passing focused model is not sufficient evidence for a user-facing
behavior, and a passing composition model is not production clearance.

## Model ownership map

| Model | Contract boundary | Primary implementation |
|---|---|---|
| `ObolProgressivePipeline` | sole owner, finite work ledger, terminal outcome | `lod_control_private.h`, `view_controller_progressive.cpp` |
| `ObolControlRefinement` | concrete facts to finite obligations | `lod_control_private.h`, `QgProgressiveDiagnostics.cpp` |
| `ObolLodComposition` | admission/growth/presentation/capacity/structural/point/static seam ownership and terminality | controller typed policies and `view_controller.cpp`, `view_controller_render.cpp` effect executors |
| `ObolTerminalConvergenceComposition` | finite cross-owner progress after inputs close, including visibility-census-to-exact-frame ordering, exact-target renderer preparation, static-to-capacity re-entry, local reconciliation, retained-ceiling constraint, and subordinate compaction | `BObolLodPlanningObligations`, `BObolLodExactPresentationFrame`, `CadPresentationPreparation`, `BObolLodStaticQualityTrial`, `BObolLodPresentationPolicy`, `BObolLodTerminalQualityReducer`, capacity search and compaction policies, plus their controller effect executors |
| `ObolControlLifecycleComposition` | policy/provider/demand/host/exact-frame composition, semantic-only exact presentation, first-deadline recovery ownership, and derived terminal status | `view_controller.cpp`, `view_controller_progressive.cpp`, `view_controller_host.cpp`, `view_controller_render.cpp` |
| `ObolHostWork` | level-triggered host/controller wakeup, including registered-versus-pending providers and exact/publication pump-to-render ownership | `view_controller_host.cpp`, `view_controller_progressive.cpp`, `window_host.cpp`, and qtcad's exact-frame delivery path |
| `ObolSemanticPresentation` | selection/style revision through an exact successor frame which is neither capacity nor LoD-planning evidence, including mutation during an in-flight predecessor and recovery only when no successor already owns the stale presentation | `BObolLodExactPresentationFrame`, the independent capacity/planning render classifications, `BObolViewController::requestExactCadPresentationRender`, `bobol_display_endpoint_request_presentation_frame`, libged selection publication, and qged manipulator presentation |
| `ObolInteractionSession` | gesture/debounce/quiet lifecycle and interrupted-replay retirement at the quiet boundary | `BObolLodInteractionSession`, `BObolViewController::notePresentationRenderInterrupted` |
| `ObolLodPolicyDisable` | atomic retirement of automatic owners while retaining presentation and provider work | `BObolViewController::synchronizeAutomaticLodControl`, `BObolLodCoordinator::retireAutomaticLodControl` |
| `ObolQuietSuccessor` | schedule-independent retained-detail restoration | `lod_quiet_successor_private.h`, `BObolLodPresentationPolicy` |
| `ObolDeadlineOwnership` | unique successor after an interrupted frame | `view_controller_render.cpp` |
| `ObolCapacitySearch` | finite renderer-capacity bracket and discrete-population reuse | `lod_capacity_search_private.h` |
| `ObolRetainedAllocationPrefix` | fixed-revision budget-to-population monotonicity, exact mixed mesh/point protected-floor certification, zero-cost cut canonicalization, and first-gap termination | `retained_allocation_private.cpp` |
| `ObolRetainedAllocationPresentation` | cut/protection commit through exact presentation, with allocation-invariant unmanaged render cost | controller retained-allocation application, `BObolViewLodState::allocationUnmanagedRenderCost`, and Obol classifier |
| `ObolCadTimingEvidence` | threshold-stamped CAD duration/reuse evidence which non-CAD frames cannot consume | `BObolLodRendererPerformanceEvidence`, `view_controller_render.cpp` |
| `ObolCompletedPassOwnership` | multi-pass growth/coverage/handoff/capacity ownership, effects, annotation and selective-delta lifetime, scene-wide capacity scope, and semantic revision | `BObolLodAvailabilityScheduler`, `BObolLodPresentationPolicy`, `BObolLodRetainedPassAnnotations`, `BObolViewController::beginSceneWideCapacitySubmission`, `view_controller.cpp` |
| `ObolCapacityPresentationHandoff` | resident-availability retry, occurrence allocation to ceiling-free sample, canonical reconciliation budget, and terminal absorption | `BObolLodAvailabilityScheduler`, `BObolLodPresentationPolicy`, `BObolLodCapacityEvidence` |
| `ObolStaticQuality` | bounded quiet-view fidelity improvement | `BObolLodStaticQualityTrial` |
| `ObolPointQualityOwnership` | exclusive point calibration/recovery owner | `BObolLodPointQualityPhase` |
| `ObolPointCalibrationDomains` | responsive/static evidence isolation and finite terminal point refinement | `BObolLodAdmissionEvidence`, point-pressure and structural-preload control |
| `ObolPointTerminalEvidence` | idempotent terminal point witness across unchanged structural censuses | `BObolLodPointProxyEvidence` |
| `ObolTerminalQualityOrdering` | exclusive ceiling-reconciliation, point-quality, static-quality, and stale-certificate recertification successor | `lod_terminal_quality_private.h`, `view_controller_render.cpp` |
| `ObolStructuralFrontierOwnership` | unresolved structural repair, bounded finer-cut preload batches, or capacity-witnessed terminal-proxy successor | `BObolLodStructuralRepair` |
| `ObolLodAdmission` | safe direct admission and bounded large-scene path | `BObolLodAdmissionPlanner` |
| `ObolLodArbitration` | budget-safe visual-importance ordering | `retained_allocation_private.cpp` |
| `ObolLodConvergence` | finite convergence and witnessed constrained endpoints for single-large and many-object scenes | controller convergence snapshot and ledger |
| `ObolSubmissionPass` | bounded submission cursor, rescan debt, level-triggered coverage-census ownership, and the exact visibility-census/frame/allocation sequence | `BObolLodSubmissionPass`, `BObolLodExactPresentationFrame`, `BObolLodPlanningObligations`, source visibility journal, controller progressive pump and host synchronization |
| `ObolResidentGrowth` | coalesced resident suffix drain, coverage transfer, and reallocation | availability ledger and planning obligations |
| `ObolResidentCompaction` | revision-safe background memory reclamation | `view_controller_residency.cpp`, `BLodService` |
| `ObolActiveProducerDemand` | superseding demand for a stable asset producer | `BLodService`, mesh submission action |
| `ObolResultAuthentication` | population-, demand-, and route-authenticated asynchronous result acceptance | `BLodService`, compact-population lookup, mesh result reducer |
| `ObolSharedAssetLease` | multi-view producer lifetime, late join, and final-consumer cancellation | `BLodService` coalescing and consumer subscription lifecycle |
| `ObolAssetPublicationComposition` | producer demand, typed constraint/failure, additive live pages, final hierarchy, and durable cache publication | `BLodService`, source presentation staging, mesh LoD cache |
| `ObolPresentationPreparation` | finite renderer-side preparation | Obol `CadPresentationPreparation`, controller render executor |
| `ObolCadViewPublication` | exact-view preparation/report acceptance | Obol `SoCADAssembly`, libBObol view attachment |
| `ObolCadMutation` | validated retained-scene publication, notification, and resource denial | Obol `CadSceneMutation`, `SoCADAssembly::replaceScene` |
| `ObolCadFrameComposition` | atomic scene mutation through exact-target preparation, host frame claim, and report acceptance | Obol `SoCADAssembly`, renderer preparation/reporting, libBObol host boundary |
| `ObolSparsePlanIdentity` | active-slot ownership when a hidden tombstone shares an `InstanceId` | `BObolCompactInstanceIndex`, sparse PoP cut patching |
| `ObolLodColdPreview` | coverage-safe cold preview, bounded admission, and per-occurrence replacement by a resident progressive binding | serialized source, compact coverage producer, mesh submission action, and view state |
| `ObolSpatialProducer` | cold spatial hierarchy and durable-cache lifecycle | mesh LoD cache/service producer |
| `ObolLiveSpatialPublication` | additive immutable page publication | compact presentation staging and Obol assembly |
| `ObolDeferredAutoview` | one-shot progressive fit ownership | libged Obol draw adapter |
| `ObolOccurrenceControl` | one-occurrence fallback/mesh/repair acknowledgement | source presentation and LoD update action |

`ObolLodColdPreview` supersedes the former `ObolColdPresentation` model.  The
older model was removed because it represented the same boundary with less
complete input, capacity, and coverage ownership.
Its occurrence-local `binding` state is deliberately stronger than a concrete
result-kind tag: `resident` means that the occurrence owns a valid progressive
hierarchy.  A prepared cold preview in the adaptive mesh channel still maps to
`coverage`.  The shared-asset update-action regression checks this refinement
boundary directly; TLC cannot infer the C++ payload classification from an
enum value.
The model additionally distinguishes a whole-prefix preview generation from
the completed spatial generation.  A nonempty private-page binding can become
`resident` only after the spatial generation exists; drawing the same logical
cut from a whole prefix is not an equivalent residency proof.  The 2026-08-31
focused rerun explored all 1,902 generated / 918 distinct states to depth 13
with no invariant, deadlock, or temporal-property failure.

## Change routing

- Change sole ownership, terminality, or progress witnesses in
  `ObolProgressivePipeline` and its focused boundary model.
- Change the concrete-to-abstract field map in `ObolControlRefinement` and its
  C++ per-field plus exhaustive abstract-owner test.
- Change visual-importance arithmetic in the retained-allocation oracle and
  property tests; TLA+ checks only its atomic budget/priority contract.
- Change CAD timing classification, point-threshold reuse, or the capacity
  sample handoff in `ObolCadTimingEvidence`, the typed host timing API, and
  renderer-performance evidence tests.
- Change retained CAD publication in `ObolCadMutation`,
  `ObolCadViewPublication`, or their seam model `ObolCadFrameComposition`, and
  change sparse retained identity ownership in `ObolSparsePlanIdentity`, plus
  Obol boundary tests.
- Change cold/cache ordering in `ObolLodColdPreview`, `ObolSpatialProducer`,
  `ObolLiveSpatialPublication`, or their seam model
  `ObolAssetPublicationComposition`, plus cold/warm GUI and memory tests.
- Change policy-off retirement or the distinction between an automatic
  capacity frame and a presentation-only repaint in `ObolLodPolicyDisable`,
  `ObolControlLifecycleComposition`, the controller off/on/off test, and the
  System GL/OSMesa LoD-off matrix.
- Change provider registration/pending semantics, exact-frame host ownership,
  or terminal-status projection in `ObolHostWork`,
  `ObolControlLifecycleComposition`, and the host/controller tests.
- Change selection/highlight presentation ownership or exact-frame
  supersession in both `ObolSemanticPresentation` and
  `ObolControlLifecycleComposition`, the exact-frame latch and render-intent
  tests, and System GL/OSMesa hierarchy-selection replays.  The focused model
  defines the transaction; the lifecycle model proves that transaction cannot
  re-enter capacity recovery when composed with host and deadline behavior.

`ObolLodPolicyDisable.AutomaticDomains` maps to production state as follows.
This is the review checklist for `retireAutomaticLodControl`; it is not a
second writable ledger.

| Domain | Production state retired atomically |
|---|---|
| service | active LoD generation and its cancellable service work |
| availability | queued-result, resident-growth, and inventory-coalescing deadline evidence |
| submission | source cursor, pass, rescan, intent, delta, and last-submitted stamps |
| planning | planning obligations, retained-pass annotations, and discrete trial permit |
| viewDemand | demand refresh, pose continuity, and view-quality probe state |
| coverage | required/active coverage census and convergence candidates |
| interaction | interaction session and its start certificate |
| presentation | handoff, publication transaction, interrupted replay, exact-frame, and point-admission frame |
| pointQuality | point-quality phase and point-proxy recovery evidence |
| staticQuality | bounded static-quality trial |
| structuralRepair | structural repair and point-relaxation evidence |
| compaction | compaction policy and resident compaction snapshot |
| capacityHost | capacity/headroom/resource-recovery evidence and capacity-labelled host render request |

Database-provider registration and immutable retained presentation are outside
those domains.  Provider *pending work* still participates in the host pump;
an idle registered provider does not.

Do not add a model for a renderer, named asset, cache temperature, or object
count.  If a proposed model cannot identify its production owner, revision
stamp, finite measure, and retirement edge, first simplify the implementation
contract.

## Running TLC

Use the catalog-aware runner so every check first validates and parses the
complete suite and always puts TLC state outside the source tree.  For example:

```text
cmake -DTLA2TOOLS_JAR=/home/cyapp/tla+/tla2tools.jar \
  -DTLA_MODE=check -DTLA_MODEL=ObolCadMutation \
  -DTLA_WORK_DIR=/tmp/obol-tlc-cad-mutation \
  -P misc/CMake/RunTLA.cmake
```

Run the canonical pair for every control-plane ownership change and the
focused models named by the table for the touched boundaries.  Convert every
counterexample into a source invariant and executable regression test; do not
accept the model fix without the implementation guard.

The following are historical verification observations through 2026-09-01.
They explain why guards exist but are not the current baseline.  The sole
machine-checked counts are in `tla/baselines/tlc-2.19.json`:

- `ObolProgressivePipeline`: 2,358,764 generated / 1,095,220
  distinct states, depth 40;
- `ObolControlRefinement`: 475,078 generated / 237,568 distinct states,
  depth 1, covering all 29 concrete input facts independently and every
  combination of the twelve semantically distinct ownership classes;
- `ObolLodComposition`: 70,952 generated / 35,848 distinct states, depth 17,
  including successful and capacity-constrained static-quality successors;
- `ObolTerminalConvergenceComposition`: 1,520 generated / 1,187 distinct
  states, depth 81, including a two-step exact visibility delta and its
  separately prepared exact framebuffer prerequisite, two exact-target
  preparation units per allocation presentation, two static populations which
  each reopen a fresh three-candidate capacity search, three strictly coarser
  reconciliation representations, a retained-ceiling constraint, and two
  background compaction steps which cannot reopen foreground work;
- `ObolControlLifecycleComposition`: 1,092,377 generated / 227,787 distinct
  states, depth 32, including two bounded retire/restart interruptions between
  a current-view demand obligation and its runnable cursor, semantic-only
  mutations before and during an exact frame, distinct capacity and semantic
  deadline-miss successors, and the requirements that every requested
  importance census owns its demand while presentation-only work cannot create
  LoD recovery;
- `ObolAssetPublicationComposition`: 5,782 generated / 1,536 distinct states,
  depth 12;
- `ObolCadFrameComposition`: 4,921 generated / 1,568 distinct states,
  depth 17, with camera and effective renderer-presentation controls sharing
  one exact-frame revision;
- `ObolHostWork`: 78,500 generated / 12,902 distinct states, depth 17, including
  idle registered providers, later provider publication, independent pump
  reasons, exact-frame debt before and after it is attached to a render, and
  batched-publication timer transfer to a coalesced render;
- `ObolSemanticPresentation`: 114 generated / 77 distinct states, depth 12,
  including a style mutation which supersedes an in-flight predecessor frame,
  one bounded lost-owner seam, non-duplicating idle recovery, and the invariant
  that semantic presentation cannot manufacture renderer-capacity or LoD-
  planning evidence;
- `ObolSubmissionPass`: 90 generated / 45 distinct states, depth 11, including
  an active coverage obligation whose cursor was paused by a stronger owner
  and a presentation-only visibility mutation whose exact source delta is
  guaranteed a level-triggered producer.  The exact census and exact
  framebuffer are independent facts: reallocation cannot begin until both are
  current.  A repeated edit supersedes an unpublished predecessor allocation,
  capacity evidence survives exact visibility changes, and quiescent
  visibility input eventually obtains an exact presentation and current
  allocation.  The 2026-09-01 extension also proves that a selective source
  pass cannot consume scene-wide demand debt: publication of new quality debt
  creates a level-triggered demand producer which survives the selective pass
  and is eventually consumed or superseded by a newer semantic revision.  Its
  `FinishSelectivePass` edge also retires selective mode with the cursor; the
  production refinement now does so before starting any successor;
- `ObolCapacitySearch`: 257,534 generated / 187,450 distinct states,
  depth 39;
- `ObolRetainedAllocationPrefix`: 97,643 generated / 55,279 distinct states,
  depth 19, including zero-cost precision upgrades at canonical terminal
  populations;
- `ObolRetainedAllocationPresentation`: 24 generated / 24 distinct states,
  depth 5, including allocation-invariant unmanaged presentation cost;
- `ObolCadTimingEvidence`: 7,726,131 generated / 1,222,796 distinct states,
  depth 12 over the full 64-cut production domain, including every bounded
  accelerated candidate inside a known bracket;
- `ObolCompletedPassOwnership`: 24,418 generated / 16,988 distinct states,
  depth 17;
- `ObolCapacityPresentationHandoff`: 159,114 generated / 17,598 distinct
  states, depth 19, including a stale-allocation resident retry which must
  resolve to drawable data or an explicit bounded constraint, exact empty and
  nonempty visibility censuses, and the requirement that a terminal capacity
  certificate atomically retires any older presentation-reconciliation budget;
- `ObolPresentationPreparation`: 1,143 generated / 526 distinct states,
  depth 12;
- `ObolCadMutation`: 257,805 generated / 128,904 distinct states, depth 5;
- `ObolCadViewPublication`: 1,422 generated / 480 distinct, depth 11;
- `ObolSparsePlanIdentity`: 17 generated / 13 distinct states, depth 7;
- `ObolResidentGrowth`: 2,562 generated / 1,866 distinct states, depth 24;
- `ObolInteractionSession`: 949 generated / 216 distinct, depth 12,
  including completed-motion and deadline-retired gates followed by an
  interrupted presentation replay at the quiet boundary;
- `ObolStructuralFrontierOwnership`: 60,535 generated / 29,266 distinct
  states, depth 386, including bounded provider batches and finite strictly-
  decreasing prefetch remainders;
- `ObolLodPolicyDisable`: 1,081,360 generated / 65,544 distinct states,
  depth 17, including off/on/off transitions from every initial reducer-domain
  combination and disabled camera changes;
- `ObolLodConvergence`: 742 generated / 423 distinct states, depth 25 after
  deadline/memory/presentation evidence became mandatory for a constrained
  terminal outcome and quiet current-view quality debt was required to retain
  a named owner until its successor demand pass completes; and
- `ObolLodAdmission`: 321 generated / 151 distinct states, depth 6, including
  cold PoP-before-discovery ordering and zero-marginal direct replacement; and
- `ObolLiveSpatialPublication`: 481 generated / 254 distinct, depth 11; and
- `ObolTerminalQualityOrdering`: 1,210 generated / 725 distinct states,
  depth 13, with typed static/reconciliation/point/capacity constraints; and
- `ObolPointTerminalEvidence`: 243 generated / 135 distinct states, depth 15.

The 2026-08-30 lifecycle composition pass closed three assumptions which the
focused models had previously proved only in isolation.  First, provider
registration is a capability while provider pending is a work fact; using the
registry as discovery state creates a permanent false owner.  Second, exact
frame debt awaiting the controller pump is distinct from exact debt already
attached to a pending or in-flight render.  Third, terminal/HUD outcome is
always derived from current evidence.  It is never stored and independently
retired.  TLC counterexamples also required a newly published provider result
to supersede an older quality owner before current-view planning resumes.
These are implementation invariants, not optional diagnostics.  Production
now performs the same exclusive PUMP-to-RENDER transfer for exact and batched
publication frames; focused host tests distinguish unfinished progressive
status from runnable pump ownership.  The lifecycle model and the subsequently
extended host model were rerun after that refinement and completed without
errors.  Naming publication in the host model found two additional composition
cases: an
exact request must preserve a still-live publication timer, and any standing
full render may absorb the publication debt.  Production's canonical pump
projection and host tests now encode both rules.

The later Generic Twin 99% stall exposed one more lifecycle seam.  Quiet-view
restoration could request a retained importance census after a pose change,
start its current-view demand pass, and then let a stronger presentation owner
retire that cursor without discharging the level-triggered demand.  The
abstract work ledger correctly refused terminal readiness, but the unchanged-
signature fast path had no transition which could recreate the executable
cursor.  `ObolControlLifecycleComposition` now distinguishes `demandPending`
from `demandCursorActive`, permits two bounded stronger-owner interruptions,
and proves that the pass is restarted and eventually completes.  Production
uses the same guarded transition and retires the importance request only from
that dense pass's completion edge.  The convergence diagnostic exports demand
refresh separately in the exact 29-fact mask, so this failure cannot hide
behind either an inactive `SUBMISSION` alias or the shared `PLANNING`
obligation during a trace audit.

`ObolControlRefinement` now names all 29 production input facts independently.
Facts which happen to refine to the same obligation remain separate so a new
producer cannot disappear behind a claimed symmetry reduction.  TLC checks
each concrete mapping and exhausts owner precedence over the twelve semantically
distinct classes; the C++ test mirrors that per-field/class split.  This avoids
an otherwise accidental full concrete-fact power set without weakening the
contract.

The twenty-ninth fact separates an interrupted closed-transaction replay from
ordinary exact-presentation debt.  The former remains the exclusive highest
owner.  The latter is a downstream effect and cannot preempt the capacity
candidate allocation which will make its requested framebuffer meaningful.
The runtime trace checker independently recomputes this precedence from the
exported facts, so the controller cannot validate its own mistaken owner.

The terminal-quality model now distinguishes ordinary pixel-demand capacity
from an exact selected presentation above its prior allocation certificate.
The latter cannot be absorbed as ownerless merely because the earlier
certificate was smaller: capacity owns a recertification or bounded
successor.  The structural-frontier model likewise makes a terminal proxy an
explicit constrained outcome available only after structural repair has a
capacity-rejection witness.  These focused facts refine the composition
model's existing exact-presentation and constraint-witness predicates; they
do not create workload-specific control modes.

The asset-publication composition check initially found a liveness hole which
the focused producer/page models could not see in isolation: a terminal
constraint belonged to an older demand, but no owner restarted the same
immutable producer for the newer demand.  The contract now resumes that
producer from its retained completed pages, and the executable service test
proves that a new view resubmits and publishes after an exact-demand failure.
Historical failure records may remain for diagnostics; only an exact current-
demand failure may suppress work.  The same composition audit removed an
over-strong assumption from `ObolLiveSpatialPublication`: every page record
must be complete before the final hierarchy/cache marker, but every page need
not have been independently presented first.  Presentation remains bounded
and optional; cache completeness remains authoritative.

The submission-pass boundary was extended after the 150k draw-metadata
migration found a terminal framebuffer coexisting with an active coverage
census and an idle cursor.  Compaction then owned the only reported obligation
but could not pass its coverage gate.  The model now requires every coverage
obligation to have an active/rescan cursor or a level-triggered pump witness
and proves that closed inventory eventually retires both the cursor and
census.  Production maps an active census to foreground inventory work,
re-arms its bounded cursor after a stronger owner retires, and preserves an
exact visibility proof across policy-only quality revisions instead of
creating a redundant census.

The same boundary now distinguishes immutable mesh inventory from effective
presentation visibility.  A hierarchy erase may leave the retained scene
correct while invalidating the LoD visibility denominator.  The model requires
every unconsumed visibility revision to have an active pass, retained rescan,
or host pump witness, and proves that a stable exact delta is eventually
applied without converting it into a complete inventory scan.  It then
requires exactly one allocation successor for the new visibility revision,
while preserving capacity evidence.  An edit arriving during that allocation
supersedes its unpublished result and consumes the newest delta first.  The
liveness property is deliberately conditional on eventually quiescent user
input; an infinite stream of visibility edits cannot promise terminal
convergence.  The production refinement is the source visibility journal, the
post-synchronization host revision check, the typed exact-visibility successor
obligation, and the focused 1/0/1 census regression.

The point-terminal-evidence boundary was added after a warm 150k trace proved
that a rejected finer structural preload correctly recorded a terminal
witness, but an immediately repeated, unchanged structural census erased it.
The focused model makes an unchanged census an idempotent observation and
allows only a semantic view, policy, population, or threshold change to renew
the consumed decision.  Its executable refinement is the no-op structural
seed regression in `test_bobol_lod_coordinator`.

The 2026-08-29 canonical progressive rerun completed with the same 2,358,764
generated / 1,095,220 distinct states and depth 40, including both temporal
branches.  The 2026-09-01 concrete refinement rerun checked all 29 current
facts and all 4,096 combinations of their twelve distinct owner classes in
475,078 generated / 237,568 distinct states.  Both canonical checks reported
no error.
The 2026-08-31 point-quality rerun completed with 66 generated /
36 distinct states and depth 8.  Production refines the point model's
`PhaseHasWitness` predicate with the typed
`PointCalibrationProducerInputs`, an exhaustive 256-case C++ truth table, and
the runtime unwitnessed-presentation violation bit.  That bridge is required:
TLC cannot prove that a host supplied the concrete producer it abstracted.
The host API therefore requires a `BObolPresentationTimingContext` whose
capacity relevance and CAD-execution facts are distinct enum types; omitted
or swapped facts fail to compile.

The timing-evidence model covers every finite safe/unsafe boundary, including
an irreducible population which overloads every candidate, while arbitrary
non-CAD frames change unrelated host timing.  Only a threshold-matched CAD
sample can narrow the bracket.  An overload may select any strictly coarser
candidate up to the known-safe side; this covers both geometric bisection and
the production timing-ratio acceleration without permitting either to cross
the bracket.  The checked-in 64-cut configuration exhausts the production
point-threshold domain and proves the same invariants and temporal property
for wide initial corrections as well as already-bracketed corrections.

The capacity model now makes allocation, exact presentation, and timing
measurement separate phases.  Its liveness proof therefore cannot spend an
invalid timing allowance while occurrence cuts are still being applied or a
temporary renderer ceiling still hides them.  The completed-pass model also
predicts the first capacity sample before the search object exists, so the
exclusive handoff owner removes that ceiling before calibration begins.  It
also models an older population frame barrier which is pending when the next
candidate is selected.  Consuming that barrier publishes a fresh semantic
capacity revision before the candidate plan commits; two allocation plans
cannot share one revision.  The current
32-budget/three-sample/four-candidate-per-goal run also makes approximate
capacity an explicit contract: after one safe population and one rejected
richer population, a goal settles rather than binary-searching the exact
adjacent budget.  Candidate publication is bounded across both ordered goals,
and a sample which misses only the steady target transfers its exact population
as the static safe lower bound.  The allocator's exact protected minimum is
part of the search domain, so no later floor override can undo a candidate and
reopen it.  TLC explored 257,534 generated / 187,450 distinct states to depth
39 with the production three-sample policy.  The expanded
resident-growth run explored 33,630 generated / 22,128 distinct states to
depth 44.

`ObolCapacitySearch` treats the aggregate classification rule and threshold as
part of its abstract immutable revision tuple.  Production refines that
abstraction by storing the exact finite threshold in
`BObolLodCapacitySearchKey`; the executable regression changes only that value
and rejects certificate reuse.  Point-versus-box raster arithmetic remains a
renderer/property-test concern rather than a useful TLA+ state dimension.

The C++ refinement now puts the complete six-domain stamp in every retained
allocation input, cache key, transaction, and result.  Its executable contract
also distinguishes the current certified plan from an older committed plan
which remains visible while replacement planning is in progress.  An
otherwise identical capacity-revision fixture proves that the older plan is
not current, cannot be reused as the new certificate, and is replaced by one
distinct plan serial.  This distinction is essential to interpreting runtime
traces: presentation continuity is state, but it is not planning authority.

The capacity search's monotone `PopulationOf` premise is no longer trusted as
an informal allocator property.  `ObolRetainedAllocationPrefix` models two
budgets executing the production first-gap rule over one immutable importance
order and proves budget safety, exact prefix cost, nested terminal populations,
equal-cost canonical population identity, zero-cost successor inclusion, and
eventual termination.  Its varied-cost fixture deliberately places a costly
upgrade before cheap and zero-cost upgrades: skipping the costly gap would
reproduce the non-nested population defect found in the 50k OSMesa trace,
while stopping before a free precision cut reproduced the ownerless 5k
allocation.  The executable refinement sweeps the mixed-scale retained
allocator fixture over increasing budgets and compares every occurrence cut,
selected cost, and population signature.

The presentation-handoff model carries a certified budget through the final
occurrence allocation and its generic changed-cut completion pass.  Those
mechanical actions preserve the terminal certificate for an unchanged
semantic problem.  The C++ refinement is the distinction between retiring an
in-flight calibration and explicitly invalidating capacity evidence.

The handoff now also distinguishes renderer work from an exact visibility
census.  Zero presented work can become a capacity sample only when the census
is authoritative and empty; a nonempty census must publish actual presentation
work first.  This closes the failure mode in which a populated compact scene
could report a completed zero-work frame and certify an empty terminal view.
Publishing a terminal capacity certificate also clears the predecessor
handoff's reconciliation and active budgets in the same transition.  TLC
found the complementary edge during this change: if reconciliation could be
requested after certificate publication, the stale lower budget could be
rearmed immediately.  The production policy now makes certificate acceptance
the single retirement operation, and the coordinator regression covers that
cross-owner boundary.
The admission model independently permits an already visible PoP payload to
be replaced atomically by affordable direct geometry after late safe-scene
certification, charging only marginal cost.  That is the formal boundary for
the Generic Twin regression in which one leaf otherwise remained progressive.

The composition model prevents those local contracts from relying on mutually
incompatible boundary assumptions.  It checks both empty and nonempty views,
safe and large profiles, every bounded initial capacity, optional independent
and capacity-owned resident growth, effective and absent renderer ceilings,
capacity sampling, structural, point, and static-quality successors, and
persistent proxies.
It proves one semantic owner, budget safety, cache-growth ownership, exact
visibility before timing/readiness, zero-work only for an exact empty view,
constraint-witnessed proxies, non-regressing direct geometry, terminal
absorption, and eventual terminality.  Its 2026-08-30 TLC run explored 70,952
generated / 35,848 distinct states to depth 17 with all temporal branches.

The canonical pair was rerun after the constrained-evidence and completed-pass
ownership changes: `ObolProgressivePipeline` again explored 2,358,764
generated / 1,095,220 distinct states to depth 40 with both temporal branches,
and the current `ObolControlRefinement` explored 475,078 generated / 237,568
distinct states to depth 1 over the exact 29-fact implementation map.  The focused
convergence model connects its abstract
nonempty witness to the runtime `constraintEvidenceMask`; the executable
refinement checker rejects a constrained outcome whose mask is empty.

All completed without errors.  These counts are diagnostic, not acceptance
criteria.  Sanitizers, randomized C++ traces, image/APNG comparison, perf,
and cold/warm System GL and OSMesa qualification remain mandatory for the
data-plane properties outside the models.
