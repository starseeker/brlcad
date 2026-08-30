# libBObol formal model catalog

Last reviewed: 2026-08-29

This is the authoritative index of the drawing stack's TLA+ models.  The
models are small proof boundaries for ownership, safety, and liveness.  They
are not alternate implementations, workload-specific controllers, or proofs
of visual quality, floating-point projection, memory use, or elapsed time.
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

## Model ownership map

| Model | Contract boundary | Primary implementation |
|---|---|---|
| `ObolProgressivePipeline` | sole owner, finite work ledger, terminal outcome | `lod_control_private.h`, `view_controller_progressive.cpp` |
| `ObolControlRefinement` | concrete facts to finite obligations | `lod_control_private.h`, `QgProgressiveDiagnostics.cpp` |
| `ObolLodComposition` | admission/growth/presentation/capacity/structural/point seam ownership and terminality | controller typed policies and `view_controller.cpp`, `view_controller_render.cpp` effect executors |
| `ObolHostWork` | level-triggered host/controller wakeup | `view_controller_host.cpp`, `window_host.cpp` |
| `ObolInteractionSession` | gesture/debounce/quiet lifecycle | `BObolLodInteractionSession` |
| `ObolLodPolicyDisable` | atomic retirement of automatic owners while retaining presentation and provider work | `BObolViewController::synchronizeAutomaticLodControl`, `BObolLodCoordinator::retireAutomaticLodControl` |
| `ObolQuietSuccessor` | schedule-independent retained-detail restoration | `lod_quiet_successor_private.h`, `BObolLodPresentationPolicy` |
| `ObolDeadlineOwnership` | unique successor after an interrupted frame | `view_controller_render.cpp` |
| `ObolCapacitySearch` | finite renderer-capacity bracket and discrete-population reuse | `lod_capacity_search_private.h` |
| `ObolRetainedAllocationPrefix` | fixed-revision budget-to-population monotonicity and first-gap termination | `retained_allocation_private.cpp` |
| `ObolCadTimingEvidence` | threshold-stamped CAD duration/reuse evidence which non-CAD frames cannot consume | `BObolLodRendererPerformanceEvidence`, `view_controller_render.cpp` |
| `ObolCompletedPassOwnership` | multi-pass growth/coverage/handoff/capacity ownership, effects, annotation and selective-delta lifetime, scene-wide capacity scope, and semantic revision | `BObolLodAvailabilityScheduler`, `BObolLodPresentationPolicy`, `BObolLodRetainedPassAnnotations`, `BObolViewController::beginSceneWideCapacitySubmission`, `view_controller.cpp` |
| `ObolCapacityPresentationHandoff` | occurrence allocation to ceiling-free sample, canonical reconciliation budget, and terminal absorption | `BObolLodPresentationPolicy`, `BObolLodCapacityEvidence` |
| `ObolStaticQuality` | bounded quiet-view fidelity improvement | `BObolLodStaticQualityTrial` |
| `ObolPointQualityOwnership` | exclusive point calibration/recovery owner | `BObolLodPointQualityPhase` |
| `ObolPointTerminalEvidence` | idempotent terminal point witness across unchanged structural censuses | `BObolLodPointProxyEvidence` |
| `ObolTerminalQualityOrdering` | exclusive ceiling-reconciliation, point-quality, and static-quality successor | `lod_terminal_quality_private.h`, `view_controller_render.cpp` |
| `ObolStructuralFrontierOwnership` | unresolved structural proxy successor | `BObolLodStructuralRepair` |
| `ObolLodAdmission` | safe direct admission and bounded large-scene path | `BObolLodAdmissionPlanner` |
| `ObolLodArbitration` | budget-safe visual-importance ordering | `retained_allocation_private.cpp` |
| `ObolLodConvergence` | finite convergence and witnessed constrained endpoints for single-large and many-object scenes | controller convergence snapshot and ledger |
| `ObolSubmissionPass` | bounded submission cursor, rescan debt, and level-triggered coverage-census ownership | `BObolLodSubmissionPass`, controller progressive pump |
| `ObolResidentGrowth` | coalesced resident suffix drain, coverage transfer, and reallocation | availability ledger and planning obligations |
| `ObolResidentCompaction` | revision-safe background memory reclamation | `view_controller_residency.cpp`, `BLodService` |
| `ObolActiveProducerDemand` | superseding demand for a stable asset producer | `BLodService`, mesh submission action |
| `ObolAssetPublicationComposition` | producer demand, additive live pages, final hierarchy, and durable cache publication | `BLodService`, source presentation staging, mesh LoD cache |
| `ObolPresentationPreparation` | finite renderer-side preparation | Obol `CadPresentationPreparation`, controller render executor |
| `ObolCadViewPublication` | exact-view preparation/report acceptance | Obol `SoCADAssembly`, libBObol view attachment |
| `ObolCadMutation` | validated retained-scene publication, notification, and resource denial | Obol `CadSceneMutation`, `SoCADAssembly::replaceScene` |
| `ObolCadFrameComposition` | atomic scene mutation through exact-target preparation, host frame claim, and report acceptance | Obol `SoCADAssembly`, renderer preparation/reporting, libBObol host boundary |
| `ObolSparsePlanIdentity` | active-slot ownership when a hidden tombstone shares an `InstanceId` | Obol `CadPlanCache`, sparse PoP cut patching |
| `ObolLodColdPreview` | coverage-safe cold preview and bounded admission | serialized source and compact coverage producer |
| `ObolSpatialProducer` | cold spatial hierarchy and durable-cache lifecycle | mesh LoD cache/service producer |
| `ObolLiveSpatialPublication` | additive immutable page publication | compact presentation staging and Obol assembly |
| `ObolDeferredAutoview` | one-shot progressive fit ownership | libged Obol draw adapter |
| `obol_lod_control` | one-occurrence fallback/mesh/repair acknowledgement | source presentation and LoD update action |

`ObolLodColdPreview` supersedes the former `ObolColdPresentation` model.  The
older model was removed because it represented the same boundary with less
complete input, capacity, and coverage ownership.

## Change routing

- Change sole ownership, terminality, or progress witnesses in
  `ObolProgressivePipeline` and its focused boundary model.
- Change the concrete-to-abstract field map in `ObolControlRefinement` and its
  exhaustive C++ test.
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
  the controller off/on/off test, and the System GL/OSMesa LoD-off matrix.

Do not add a model for a renderer, named asset, cache temperature, or object
count.  If a proposed model cannot identify its production owner, revision
stamp, finite measure, and retirement edge, first simplify the implementation
contract.

## Running TLC

Always put TLC state files outside the source tree.  For example:

```text
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar \
  -workers 1 -metadir /tmp/obol-tlc-cad-mutation \
  -config doc/notes/ObolCadMutation.cfg \
  doc/notes/ObolCadMutation.tla
```

Run the canonical pair for every control-plane ownership change and the
focused models named by the table for the touched boundaries.  Convert every
counterexample into a source invariant and executable regression test; do not
accept the model fix without the implementation guard.

Current focused verification through 2026-08-29:

- `ObolProgressivePipeline`: 2,358,764 generated / 1,095,220
  distinct states, depth 40;
- `ObolControlRefinement`: 4,194,302 generated / 2,097,152
  distinct states, depth 1;
- `ObolLodComposition`: 35,600 generated / 18,136 distinct states, depth 15;
- `ObolAssetPublicationComposition`: 5,782 generated / 1,536 distinct states,
  depth 12;
- `ObolCadFrameComposition`: 8,186 generated / 2,338 distinct states,
  depth 18;
- `ObolHostWork`: 2,007 generated / 714 distinct states, depth 15;
- `ObolSubmissionPass`: 10 generated / 7 distinct states, depth 3, including
  an active coverage obligation whose cursor was paused by a stronger owner;
- `ObolCapacitySearch`: 261,250 generated / 182,647 distinct states,
  depth 60;
- `ObolRetainedAllocationPrefix`: 81,371 generated / 46,576 distinct states,
  depth 19;
- `ObolCadTimingEvidence`: 7,726,131 generated / 1,222,796 distinct states,
  depth 12 over the full 64-cut production domain, including every bounded
  accelerated candidate inside a known bracket;
- `ObolCompletedPassOwnership`: 7,538 generated / 5,684 distinct states,
  depth 17;
- `ObolCapacityPresentationHandoff`: 183,306 generated / 21,438 distinct
  states, depth 20, including exact empty and nonempty visibility censuses;
- `ObolPresentationPreparation`: 1,143 generated / 526 distinct states,
  depth 12;
- `ObolCadMutation`: 257,805 generated / 128,904 distinct states, depth 5;
- `ObolCadViewPublication`: 1,422 generated / 480 distinct, depth 11;
- `ObolSparsePlanIdentity`: 17 generated / 13 distinct states, depth 7;
- `ObolResidentGrowth`: 2,562 generated / 1,866 distinct states, depth 24;
- `ObolInteractionSession`: 468 generated / 132 distinct, depth 11;
- `ObolStructuralFrontierOwnership`: 27,246 generated / 20,938 distinct
  states, depth 386;
- `ObolLodPolicyDisable`: 7,436 generated / 840 distinct states, depth 8,
  including off/on/off transitions from every initial automatic-work
  combination;
- `ObolLodConvergence`: 660 generated / 399 distinct states, depth 27 after
  deadline/memory/presentation evidence became mandatory for a constrained
  terminal outcome; and
- `ObolLodAdmission`: 321 generated / 151 distinct states, depth 6, including
  cold PoP-before-discovery ordering and zero-marginal direct replacement; and
- `ObolLiveSpatialPublication`: 481 generated / 254 distinct, depth 11; and
- `ObolTerminalQualityOrdering`: 245 generated / 170 distinct states,
  depth 11; and
- `ObolPointTerminalEvidence`: 243 generated / 135 distinct states, depth 15.

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

The point-terminal-evidence boundary was added after a warm 150k trace proved
that a rejected finer structural preload correctly recorded a terminal
witness, but an immediately repeated, unchanged structural census erased it.
The focused model makes an unchanged census an idempotent observation and
allows only a semantic view, policy, population, or threshold change to renew
the consumed decision.  Its executable refinement is the no-op structural
seed regression in `test_bobol_lod_coordinator`.

The 2026-08-29 canonical progressive rerun completed with the same 2,358,764
generated / 1,095,220 distinct states and depth 40, including both temporal
branches.  The concrete refinement rerun independently enumerated all
2,097,152 fact combinations (4,194,302 generated transitions), and both checks
reported no error.  The point-quality rerun completed with 29 generated /
17 distinct states and depth 5.  Production refines the point model's
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
32-budget/four-sample/four-candidate-per-goal run also makes approximate
capacity an explicit contract: after one safe population and one rejected
richer population, a goal settles rather than binary-searching the exact
adjacent budget.  Candidate publication is bounded across both ordered goals,
and a sample which misses only the steady target transfers its exact population
as the static safe lower bound.  TLC explored 189,853 generated / 139,661
distinct states to depth 44.  The expanded
resident-growth run explored 33,630 generated / 22,128 distinct states to
depth 44.

`ObolCapacitySearch` treats the aggregate classification rule and threshold as
part of its abstract immutable revision tuple.  Production refines that
abstraction by storing the exact finite threshold in
`BObolLodCapacitySearchKey`; the executable regression changes only that value
and rejects certificate reuse.  Point-versus-box raster arithmetic remains a
renderer/property-test concern rather than a useful TLA+ state dimension.

The C++ refinement now puts the complete five-domain stamp in every retained
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
equal-cost population identity, and eventual termination.  Its varied-cost
fixture deliberately places a costly upgrade before cheap upgrades: skipping
that gap would reproduce the non-nested population defect found in the 50k
OSMesa trace.  The executable refinement sweeps the mixed-scale retained
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
The admission model independently permits an already visible PoP payload to
be replaced atomically by affordable direct geometry after late safe-scene
certification, charging only marginal cost.  That is the formal boundary for
the Generic Twin regression in which one leaf otherwise remained progressive.

The composition model prevents those local contracts from relying on mutually
incompatible boundary assumptions.  It checks both empty and nonempty views,
safe and large profiles, every bounded initial capacity, optional independent
and capacity-owned resident growth, effective and absent renderer ceilings,
capacity sampling, structural and point successors, and persistent proxies.
It proves one semantic owner, budget safety, cache-growth ownership, exact
visibility before timing/readiness, zero-work only for an exact empty view,
constraint-witnessed proxies, non-regressing direct geometry, terminal
absorption, and eventual terminality.  Its 2026-08-29 TLC run explored 35,600
generated / 18,136 distinct states to depth 15 with all temporal branches.

The canonical pair was rerun after the constrained-evidence and completed-pass
ownership changes: `ObolProgressivePipeline` again explored 2,358,764
generated / 1,095,220 distinct states to depth 40 with both temporal branches,
and `ObolControlRefinement` explored 4,194,302 generated / 2,097,152 distinct
states to depth 1.  The focused convergence model connects its abstract
nonempty witness to the runtime `constraintEvidenceMask`; the executable
refinement checker rejects a constrained outcome whose mask is empty.

All completed without errors.  These counts are diagnostic, not acceptance
criteria.  Sanitizers, randomized C++ traces, image/APNG comparison, perf,
and cold/warm System GL and OSMesa qualification remain mandatory for the
data-plane properties outside the models.
