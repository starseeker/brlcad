# libBObol formal model catalog

Last reviewed: 2026-08-28

This is the authoritative index of the drawing stack's TLA+ models.  The
models are small proof boundaries for ownership, safety, and liveness.  They
are not alternate implementations, workload-specific controllers, or proofs
of visual quality, floating-point projection, memory use, or elapsed time.

## Canonical refinement

`ObolProgressivePipeline.tla` is the canonical finite control contract.
`ObolControlRefinement.tla` maps the concrete C++ obligation facts onto that
contract.  Every other model isolates one boundary of the same transition
relation.  A focused model may add detail, but it may not introduce another
production owner or control mode.

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
| `ObolHostWork` | level-triggered host/controller wakeup | `view_controller_host.cpp`, `window_host.cpp` |
| `ObolInteractionSession` | gesture/debounce/quiet lifecycle | `BObolLodInteractionSession` |
| `ObolQuietSuccessor` | schedule-independent retained-detail restoration | `lod_quiet_successor_private.h`, `BObolLodPresentationPolicy` |
| `ObolDeadlineOwnership` | unique successor after an interrupted frame | `view_controller_render.cpp` |
| `ObolCapacitySearch` | finite renderer-capacity bracket and discrete-population reuse | `lod_capacity_search_private.h` |
| `ObolCompletedPassOwnership` | multi-pass growth/coverage/handoff/capacity ownership, effects, annotation lifetime, and semantic revision | `BObolLodAvailabilityScheduler`, `BObolLodPresentationPolicy`, `BObolLodRetainedPassAnnotations`, `view_controller.cpp` |
| `ObolCapacityPresentationHandoff` | occurrence allocation to ceiling-free sample | `BObolLodPresentationPolicy` |
| `ObolStaticQuality` | bounded quiet-view fidelity improvement | `BObolLodStaticQualityTrial` |
| `ObolPointQualityOwnership` | exclusive point calibration/recovery owner | `BObolLodPointQualityPhase` |
| `ObolStructuralFrontierOwnership` | unresolved structural proxy successor | `BObolLodStructuralRepair` |
| `ObolLodAdmission` | safe direct admission and bounded large-scene path | `BObolLodAdmissionPlanner` |
| `ObolLodArbitration` | budget-safe visual-importance ordering | `retained_allocation_private.cpp` |
| `ObolLodConvergence` | finite convergence for single-large and many-object scenes | controller convergence snapshot and ledger |
| `ObolSubmissionPass` | bounded submission cursor and rescan debt | `BObolLodSubmissionPass` |
| `ObolResidentGrowth` | coalesced resident suffix drain, coverage transfer, and reallocation | availability ledger and planning obligations |
| `ObolResidentCompaction` | revision-safe background memory reclamation | `view_controller_residency.cpp`, `BLodService` |
| `ObolActiveProducerDemand` | superseding demand for a stable asset producer | `BLodService`, mesh submission action |
| `ObolPresentationPreparation` | finite renderer-side preparation | Obol `CadPresentationPreparation`, controller render executor |
| `ObolCadViewPublication` | exact-view preparation/report acceptance | Obol `SoCADAssembly`, libBObol view attachment |
| `ObolCadMutation` | validated retained-scene publication, notification, and resource denial | Obol `CadSceneMutation`, `SoCADAssembly::replaceScene` |
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
- Change retained CAD publication in `ObolCadMutation` or
  `ObolCadViewPublication`, and change sparse retained identity ownership in
  `ObolSparsePlanIdentity`, plus Obol boundary tests.
- Change cold/cache ordering in `ObolLodColdPreview`, `ObolSpatialProducer`,
  or `ObolLiveSpatialPublication`, plus cold/warm GUI and memory tests.

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

Current focused verification through 2026-08-28:

- `ObolProgressivePipeline`: 2,358,764 generated / 1,095,220
  distinct states, depth 40;
- `ObolControlRefinement`: 4,194,302 generated / 2,097,152
  distinct states, depth 1;
- `ObolSubmissionPass`: 7 generated / 5 distinct states, depth 4;
- `ObolCapacitySearch`: 4,704 generated / 4,368 distinct states,
  depth 37;
- `ObolCompletedPassOwnership`: 2,334 generated / 1,836 distinct states,
  depth 17;
- `ObolPresentationPreparation`: 1,143 generated / 526 distinct states,
  depth 12;
- `ObolCadMutation`: 257,805 generated / 128,904 distinct states, depth 5;
- `ObolCadViewPublication`: 1,422 generated / 480 distinct, depth 11;
- `ObolSparsePlanIdentity`: 17 generated / 13 distinct states, depth 7;
- `ObolResidentGrowth`: 525 generated / 377 distinct states, depth 23;
- `ObolInteractionSession`: 468 generated / 132 distinct, depth 11; and
- `ObolLiveSpatialPublication`: 343 generated / 212 distinct, depth 11.

The 2026-08-28 canonical progressive rerun completed in 95 seconds with the
same 2,358,764 generated / 1,095,220 distinct states and depth 40, including
both temporal branches.  The point-quality rerun completed with 29 generated /
17 distinct states and depth 5.  Production refines the point model's
`PhaseHasWitness` predicate with the typed
`PointCalibrationProducerInputs`, an exhaustive 256-case C++ truth table, and
the runtime unwitnessed-presentation violation bit.  That bridge is required:
TLC cannot prove that a positional Boolean call site supplied the concrete
producer it abstracted.

The capacity model now makes allocation, exact presentation, and timing
measurement separate phases.  Its liveness proof therefore cannot spend an
invalid timing allowance while occurrence cuts are still being applied or a
temporary renderer ceiling still hides them.  The completed-pass model also
predicts the first capacity sample before the search object exists, so the
exclusive handoff owner removes that ceiling before calibration begins.

All completed without errors.  These counts are diagnostic, not acceptance
criteria.  Sanitizers, randomized C++ traces, image/APNG comparison, perf,
and cold/warm System GL and OSMesa qualification remain mandatory for the
data-plane properties outside the models.
