# libBObol active debt

Last reviewed: 2026-09-01

This is the sole remaining-work list for the Obol drawing stack.  It records
work which is not complete; resolved failure analysis belongs in
`libbobol_engineering_lessons.md`, and executable release rows belong in
`obol_production_readiness.md`.  Historical timing directories and
chronological debugging logs are not design authority.

## Current implementation baseline

The following architecture is implemented and must not be reintroduced as
debt under another name:

- libged owns renderer-neutral draw roots, exact/subtree erase, selection,
  edit scopes, annotations, polygons, faceplate state, and typed scene deltas.
- libBObol owns database realization, compact occurrence inventories, cache
  and worker services, view demand, resource policy, presentation planning,
  and convergence evidence.
- Obol owns retained part/instance storage, renderer planning/execution,
  view-local preparation, picking primitives, and completed-frame evidence.
- Each view owns a distinct `SoCADAssembly` and `CadViewState`.  Immutable,
  admitted `PartGeometry` may be shared; camera-local plans and reports may
  not.
- Part and occurrence identities are distinct fixed-width strong types.
  Persistent mesh topology uses fixed-width values and is validated before it
  becomes renderer-visible.
- `PartGeometryBuilder` is the sole mutable producer form.  Its admission
  validates once and constructs a private, const `PartGeometry` snapshot;
  cached snapshots recover tokens in O(1), and renderers cannot receive a
  mutable alias.  Scene mutation batches have pure preflight functions and are
  atomic with respect to validation failure.  libBObol stages geometry,
  instances, and semantics without touching the live assembly, then commits
  through checked mutation calls.
- Complete-scene replacement constructs its candidate before the no-throw
  publication swap and reports resource denial without changing the live
  scene.  Sparse publication remains proportional to its bounded journal:
  validation and libBObol staging are atomic, but Obol does not clone a 150k
  retained scene to recover from process-wide allocator exhaustion during the
  mechanical commit.  Capacity is reserved from the manifest and sparse
  journal size remains bounded.
- Obol update windows are move-only, nest-safe RAII values.  Manual
  `beginUpdate()`/`endUpdate()` pairing is no longer public.
- One retained assembly presentation never chooses CAD LoD.  It executes the
  producer-certified cut and reports exact cost/resource evidence.
- Progressive mesh policy reads one shared immutable generation snapshot per
  allocation candidate; enrichment or compaction publishes a new generation
  without invalidating an in-flight transaction.  Renderer timing evidence is
  point-threshold stamped, and hosts must classify capacity relevance and
  actual CAD execution with distinct required types.
- Progressive spatial pages are storage/render partitions of one logical CAD
  occurrence.  They inherit its transform, style, selection, picking, and
  semantic path and are not promoted to independent CAD objects.
- Optional PCA-oriented bounds are immutable cache/part metadata, not a LoD
  state.  Monolithic retained parts use them for screen-significant batched
  aggregate boxes in shaded and wire modes; AABB remains the validated
  fallback.  The optional projection pass occurs after cold coverage can be
  published, so it cannot delay the first useful preview.
- Workload names and renderer names are qualification profiles, not control
  modes.  Lucy, Hubble, 50k, 150k, System GL, and OSMesa all use the same
  state transition relation.

Focused Obol CAD tests, libBObol rendering/realization/update tests,
`ged_test_obol_draw_sync`, the pivot guard, and the installed-package consumer
pass after the 2026-08-28 control-contract cleanup.  This is a development
gate, not production clearance.

## P0: finish the control-state reduction

The canonical contract is `libbobol_progressive_pipeline_contract.md`: one
six-domain evidence stamp (inventory, availability, visibility, view, policy,
capacity),
at most one bounded plan cursor, at most one presentation transaction, and a
finite work ledger.  HUD outcome is a projection of those values, never a
second phase machine.

Much of this contract is implemented: typed revision advancement, stale-plan
rejection, finite interaction/static-quality/point-quality states, one
presentation transaction, bounded preparation evidence, and the control
refinement map all have focused tests.  The former 6.1k-line coordinator
header is now a small umbrella over independently compile-checked admission,
capacity, delivery, scene-evidence, presentation, and view-policy boundaries.
Submission source repositioning now atomically resets its source-local plan,
allocation cursor, and visibility census.  Fresh/retired submission passes
atomically consume predecessor rescan state, while deliberate inventory pauses
retain it.  Deadline-safe cuts, retained visual-error bounds, and renderer
timing/upload observations are keyed allocation-free evidence values rather
than independently writable companion scalars.
Terminal renderer-capacity evidence now remains authoritative through both
the final mechanical cut application and subsequent ordinary admission
pumps.  The completed allocation is one typed policy transition: it either
retires only in-flight measurement while preserving a current certificate, or
atomically invalidates the certificate and its budget-limited witness.  A
current terminal certificate caps ordinary planning at its safe budget, so a
throughput estimate cannot immediately reselect the rejected discrete
population and recreate the same search.
Constrained terminal outcomes now expose a typed evidence mask and the runtime
refinement checker rejects an unwitnessed constraint.  This distinguishes an
honest deadline/memory/presentation endpoint from a controller which merely
stops while visual debt remains.  Keep this evidence derived from the owning
policy values; do not add another writable outcome-reason field.
`ObolLodComposition.tla` now closes the formal seam between admission,
retained growth, exact presentation, capacity, structural repair, point
quality, static quality, and terminal publication.
`ObolTerminalConvergenceComposition.tla` closes the previously atomic tail of
that seam.  It permits the production relationships which matter near
termination—an exact visibility census followed by its exact framebuffer
classification before reallocation, a distinct static population reopening a
bounded capacity search,
each exact presentation consuming monotone renderer-preparation units, an
over-budget allocation consuming a strictly coarser local representation, and
quiet compaction after foreground convergence—but requires every re-entry to
consume a finite semantic member.  Its current 2026-08-31 TLC run checked 639
  distinct states to depth 81.  The corresponding executable handoff now
removes the temporary global ceiling after a successful occurrence-local
allocation; only a revision-bound protected-minimum constraint may retain it.
This eliminates the former path which re-walked global PoP ordinals after the
local allocation had already succeeded.
`ObolControlLifecycleComposition.tla` closes the orthogonal policy/provider/
host seam.  It proves that provider registration is not provider work,
policy-off camera bookkeeping cannot arm automatic demand, exact-frame debt
distinguishes work awaiting a request from work already attached to a frame,
the first capacity-relevant hard-deadline miss owns quality recovery and a
fresh exact-frame obligation, and terminal/HUD outcome is derived rather than
stored.  The 2026-09-01 TLC pass checked 1,092,377 generated / 227,787 distinct
lifecycle states to depth 32.  It includes two bounded interruptions which
retire a runnable demand cursor without discharging its level-triggered
importance-census obligation, and semantic-only selection/style mutations
before and during exact presentation.  Those mutations advance their own
presentation revision but preserve every LoD control fact; even an interrupted
semantic-only frame may request only an exact repaint, never capacity recovery.
The subsequently hardened
28-fact refinement map checked 57,344 states: every concrete field
independently and every combination of its ten distinct owner classes.  The
offline checker now requires that complete fact mask, independently projects
its obligations and fixed-precedence owner, and includes it in A/B/A cycle
identity.  Its adversarial gate rejects alias-only cycles, owner/obligation
mismatches, a presentation owner without a typed finite successor, and a
missing concrete mask.  Dense 50k OSMesa tracing then found two refinement
defects which sparse event samples had missed: a consumed capacity candidate
could publish a replacement plan without advancing the capacity domain, and
the runtime validator did not recognize the controller-scoped PUMP-to-RENDER
transition already present in both `ObolControlLifecycleComposition` and
`ObolHostWork`.  Capacity replacement now advances its semantic revision, and
presentation diagnostics export the exact witness source; an unrelated shared
host pump cannot satisfy that check.  The resulting dense trace passes the
refinement checker, although the 30-second 50k OSMesa gate can still expire
while monotone command preparation and worker realization are active.  This
was distinguished from an end-state cycle by an extended full workflow: the
warm 50k OSMesa selection/erase/redraw replay completed in 31.6 seconds, and
its three post-selection convergence waits returned to ownerless readiness in
0.55--0.92 seconds with no obligation or violation.  This
does not close the remaining production refinement debt: not every imperative
effect writer emits a typed transition record yet, and preparation's remaining
unit rank is not yet exported into the dense trace, so the checker still
samples states rather than checking a complete concrete execution against the
transition relation.
Finish that migration at semantic ledger boundaries and delete each superseded
writer in the same change.  Do not add another policy model for the same state,
and do not create a monolithic effect facade around already cohesive reducers.
The 2026-08-30 implementation refinement now gives exact-frame debt distinct
request-required and frame-awaited states, including reattachment to a
coalesced host request after a newer semantic mutation supersedes its target.
Provider registration and provider pending remain separate at the host/status
boundary.  One derived controller-local pump projection covers all reducer
obligations, while service queues and provider work are composed at the host
boundary.  The policy-off transaction asserts the complete modeled retirement
postcondition and also retires the inventory-coalescing deadline.  The
controller-to-ledger projection is now canonical as well: inventory
coalescing, visibility deferral, source deltas, quality probes, resumable
retained allocation, importance census, and resident-admission retry can no
longer keep the pump active while disappearing from convergence diagnostics.
The public diagnostic snapshot now also carries the concrete 28-fact mask;
this distinguishes those aliases without exposing private reducer objects.
Focused
off/on/off, idle-provider, exact-frame, and host-wakeup tests exercise these
bridges.  Requesting an exact or batched-publication frame now transfers the
runnable host level from PUMP to RENDER instead of polling an action which
cannot advance before that frame completes; `BObolProgressiveStatus::hasMore`
still reports the unfinished transaction without becoming a second scheduler.
The HUD audit found no independently writable terminal outcome or controller
phase: its phase and completion are projections of current evidence, while Qt
retains only publication-change observations.
The 2026-08-29 ownership audit retired the last detached lifecycle combinations
relevant to proxy admission: structural relaxation is one four-state value,
the recovery-plan witness lives with point-quality state, renderer-feedback
consumption lives with the interaction session, and the host render/capacity
request is one three-way value.
`ObolAssetPublicationComposition.tla` and `ObolCadFrameComposition.tla` also
close the adjacent producer/cache and retained-scene/frame seams.  The first
composition run found a real ownerless successor: a producer constrained under
an old view demand did not resume for a newer demand.  Production now retries
the current demand, and the focused service test proves a new-view result is
actually published rather than merely clearing the old failure query.  The CAD
frame seam rejects stale target reports across concurrent scene/view changes.
The production shape is nevertheless still concentrated in two places:

- `lod_admission_policy_private.h` now retains the allocation-free evidence
  values and planner declarations; the 71-method numeric planner is compiled
  independently in `lod_admission_policy.cpp`.  Point and structural evidence
  remain in the same 1.6k-line header because their small methods preserve the
  trivially-copyable value contract, but any further extraction must keep the
  header-alone compile and exhaustive policy tests.
- `view_lod_coordinator_state_private.h` is about 1.3k lines.  Submission,
  repair, capacity, presentation, point-quality, and interaction ownership are
  cohesive typed values.  Further changes should remove a measured effect-
  shell responsibility, not add accessors or a generic lifecycle wrapper.
- `view_controller.cpp` is about 6.9k lines after lifetime/camera/host-request,
  exact-picking, and residency extraction.  It remains both reducer caller and
  the LoD policy effect executor.  Exact-view quality history is now an
  independently compiled allocation-free policy rather than another inline
  controller implementation.

Required completion work:

1. Finish inventorying the remaining mutable controller fields by sole owner,
   revision domain, progress witness, and terminal transition.  Delete
   write-only or derivable fields; do not split a keyed certificate back into
   convenient writable scalars.
2. Continue moving nontrivial pure-policy method bodies behind the new
   compile-checked private boundaries.  Keep hot occurrence storage dense and
   allocation-free; this is a responsibility split, not an object hierarchy.
3. Route each remaining direct lifecycle mutation through the typed reducer
   which owns that ledger and returns a complete successor plus bounded
   effects.  Submission cursor positioning, fresh/retired pass transitions,
   structural relaxation, point-recovery completion, interaction feedback,
   and host render requests now satisfy this rule.  Deliberate pause/resume
   operations preserve rescan debt explicitly.  An unmapped transition is
   debt unless it represents a genuinely new semantic event.
4. Extend randomized refinement traces as each remaining effect writer moves
   behind the reducer, and retain the exhaustive focused value tests.  The
   completed-pass boundary now has full selector enumeration plus 512 seeded
   multi-pass traces.  Every nonterminal state must identify a worker, cursor,
   scheduled owner-thread pump, frame acknowledgement, or finite timer which
   can make progress.
5. Re-run the applicable TLC models after ownership/liveness changes and the
   graphical matrix after numeric policy changes.  Formal models do not judge
   visual quality or wall-clock performance.

Acceptance: unchanged evidence cannot reopen planning; an invalid/stale plan
cannot commit; one event cannot select two successor owners; no terminal HUD
state has foreground work; and no nonterminal state is ownerless.

Two 2026-08-28 150k System-GL traces exposed the same missing refinement at
different effect boundaries.  A point-calibration request pauses the active
submission cursor, so that cursor cannot also be counted as the frame producer
which will present the calibration.  Producer classification is now one typed
input value, all 256 input combinations have an executable truth-table test,
and the controller recomputes ownership after deadline effects mutate the
calibration state.  Runtime refinement validation additionally rejects a
`PRESENTATION` owner without a requested frame, independent producer, or
finite publication timer.  Retain those checks while completing the reducer;
an owner label, progressive-pump level, or paused cursor is not a progress
witness.

After immutable generation snapshots removed repeated shared-generation loads
and CAD timing became threshold stamped, the 2026-08-28 warm 150k shaded
replays terminated box-free in about 40.4 seconds on System GL and 36.8 seconds
on OSMesa.  The prior 32/64-pixel balancing cycle was gone: terminal samples
have no owner, obligation, pending foreground/background work, or control
violation.  This resolves the known same-population liveness defect, but not
all scale debt.  Profile true-cold useful-preview latency and cache I/O using
common policy events rather than an object-count regime.

The stricter 2026-08-31 matched-quality workflow selects materially richer
populations.  Its initial 150k System-GL view reaches 20.49M faces through
78,341 meshed occurrences in 220 seconds; OSMesa reaches its constrained
1.17M-face endpoint in 51 seconds.  Both are ownerless, box-free, proxy-free,
and have no prominent-floor violation.  A current 50k `perf` replay confirms
that the former 13.6% Obol scene-mutation snapshot hotspot has been retired;
no current exclusive symbol exceeds 3%.  Treat the remaining 150k latency as
distributed realization/publication throughput debt.  Instrument phase-level
work and test batching or parallel publication before proposing another local
micro-optimization.  `obol_lod_visual_quality.md` owns the exact metrics and
the distinction between matched controls and unsafe whole-scene controls.

The 2026-08-29 draw-metadata migration replay exposed a separate cross-owner
gap: a policy revision could create a coverage census after the exact
camera/source proof was already current, retire its cursor through a stronger
capacity path, and leave compaction waiting on the orphaned census.  Coverage
is now an explicit inventory obligation, its producer is level-triggered, and
policy-only quality changes preserve the existing proof.  The corrected warm
System-GL 50k/150k matrices both pass; the 50k endpoint is fully ownerless,
while 150k presents a terminal box-free framebuffer and truthfully continues
bounded resident-prefix compaction in the background.  Do not weaken that
distinction by treating background reclamation as unfinished visual work.

A subsequent 50k exact-subpath erase exposed a different revision-boundary
defect: compact retained visibility changed without changing immutable mesh
inventory, so the presentation was correct while the convergence denominator
remained at 50,000.  Mesh inventory and effective presentation visibility now
have independent monotonic revisions and bounded journals.  The controller
unions their changed-entry sets, while the host performs a level-triggered
revision check after presentation synchronization.  Exact visibility is also
a first-class revision in retained-allocation plans.  It requires one
successor allocation after its bounded source delta is applied and the
resulting framebuffer has been classified exactly.  The source census alone
is insufficient: a restored occurrence may have lost its retained payload
while hidden, and only the successor frame exposes the structural replacement
work which must precede allocation.  This edge does not
invalidate the renderer-capacity certificate: visibility changes allocation,
not the cost model.  A newer exact edit supersedes an unpublished predecessor
allocation so the latest delta cannot be stranded behind stale work.

The final warm OSMesa 50k workflow updates 50,000 to 49,984 and back in
1.36/1.40 seconds while retaining 1,206,806 faces, zero boxes, and the same
5,130,680-unit certified budget.  The corresponding 150k workflow updates
150,000 to 149,984 and back in 6.05/3.96 seconds while retaining 835,318 faces,
zero boxes, and the same 1,446,651-unit budget within one rounding unit.  Both
terminate ownerless and the six-domain runtime trace passes.  Preserve this
distinction: visibility
is a planning input, not a geometry-inventory or renderer-capacity mutation,
and presentation-only edits still require a progress witness.

A later timing-sensitive 50k replay made this ordering failure deterministic:
erase-time allocation retired three visible mesh payloads, redraw updated the
50,000-occurrence census, and immediate allocation inspected the predecessor
frame.  Sixteen restored occurrences were absent from both the terminal ledger
and every producer.  The controller now requires an exact presentation frame
on every exact visibility delta, and
`BObolLodPlanningObligations::exactVisibilityReallocationReady` is the sole
gate from that frame to reallocation.  The strengthened submission and
terminal-composition models prove the ordering and quiescent liveness.  The
post-fix warm 50k OSMesa workflow at
`/tmp/qged-visibility-prerequisite-50k-osmesa-20260901` restores all 50,000
occurrences in 0.83 seconds, box-free and ownerless.

The pre-typed-host 150k OSMesa report had three quality-floor misses, including
one prominent synthetic occurrence.  The rebuilt host rerun has zero total or
prominent floor misses and zero visual-importance debt, confirming that the
old faceplate/CAD timing mix affected the selected capacity rather than proving
an allocator-ordering failure.  Realistic wheels/blades/hulls must still
demonstrate that the scene-wide importance ordering preserves genuinely
prominent forms.  Treat any future unaffordable residual as numeric
allocation/qualification debt; do not reopen a terminal control search merely
because proven-constrained quality demand remains.

## P0: complete physical responsibility extraction

The initial extraction and the Obol renderer split are real improvements, but
several units still obstruct review:

- `database_source.cpp` is about 17.6k lines.  Compact CAD presentation is now
  an independently compiled 2.1k-line owner, and its bounded copy-on-write
  staging journal remains side-effect-free until Obol accepts the complete
  renderer transaction.  BoT geometry conversion and normal construction are
  also independently compiled worker-safe units.  The next genuine boundary
  is discovery/traversal versus worker realization orchestration; do not move
  hot compact records into per-occurrence objects while extracting it.
- `draw_obol.cpp` is about 12.9k lines after endpoint, geometry, overlay, and
  scene-record extraction.  The libged private mutation vocabulary is now
  consistently named `ged_scene_reducer_*`; separate reducer orchestration
  from source realization/backend effects without reintroducing an adapter
  transaction type.
- Obol's former 5.8k-line `SoCADAssembly.cpp` is now a roughly 2.1k-line Coin
  node/mutation/action surface.  Retained plan maintenance is compiled in
  `CadAssemblyPlan.cpp`, camera-local subpixel classification in
  `CadAssemblyClassification.cpp`, and their state-only private contract in
  `CadAssemblyImpl.h`.  Renderer indirect, instanced, flattened, and
  retained/direct executors remain separate.  Further extraction should occur
  only at a measured retained-cache, picking, or Coin-action boundary.

Do not replace large files with textual `.inc` fragments, mirrored state, or
thin forwarding abstractions.  Independent compilation must enforce one-way
dependencies.  Preserve sparse deltas, bounded 512-record publication,
immutable geometry ownership, and the no-second-150k-vector rule.

Remove superseded APIs and code in the same change which replaces them.
Compatibility with intermediate branch APIs/cache formats is not required.

## P0: graphical release qualification

Run the exact final binaries using true-cold isolated caches and certified
warm caches.  Required cases and controls are enumerated in
`obol_production_readiness.md`; the minimum model set is Generic Twin, Lucy,
multi-Lucy and xpushed multi-Lucy, Hubble, Havoc, NIST BREPs, heterogeneous
50k/150k scenes, Stanford meshes, and independent multi-gigabyte vehicle
models.

For both System GL and OSMesa, cover shaded, wire, hidden-line/evaluated modes
where applicable, LoD on/off, resize and fractional DPR, zoom, rotation,
translation, selection, exact/subpath erase-redraw, and memory turnover.
Compare System GL and OSMesa semantics/images within declared tolerances.
Use APNG and apitrace when diagnosing corruption, flashes, or camera jumps.
Use the matched full-detail methodology and provisional corpus targets in
`obol_lod_visual_quality.md`; libicv SSIM/PHASH and silhouette disagreement are
required evidence, not substitutes for named-feature inspection.

Release evidence must prove:

- a useful bounded initial presentation and truthful cold-work HUD;
- monotone useful convergence or an explicit typed resource constraint;
- no box/mesh cycling, mesh-to-box regression, holes from invalid PoP
  topology, lingering fallback boxes, or premature `View ready`;
- protected prominent shapes meet the visual-significance floor;
- rotation/translation restore proven affordable detail, while zoom starts
  from retained detail and changes demand incrementally;
- invisible/noncritical residency is reclaimed under pressure without OOM;
- input remains responsive and expensive rendering is interruptible; and
- repeated view turnover settles without leaked work or zombie GUI processes.

Run ASan/UBSan across the shared dynamic stack with worker teardown, corrupt
cache input, endpoint replacement, plugin reload, edit cancellation, and rapid
view close.  Run TSan/LSan where the native runtime supports them.  Static
linkage is reserved for explicit static-link tests.

## P0: interaction and editing completion

`qged_editing.md` owns the detailed editing plan.  The remaining production
gate explicitly includes:

- classify every librt primitive operation and audit advertised-operation
  rejection, constraint preservation, error reporting, and command readback;
- complete the reusable manipulator vocabulary (including ARB vertex/edge/
  face interaction) and richer primitive-specific manipulators;
- implement specialized sketch curve/profile interaction using the current
  libqtcad widget only where it remains useful; delete obsolete demos/widgets;
- validate qged and MGED/gsh sessions against the same edit command state,
  including final MGED `sed`/`oed` qualification;
- drive real mouse paths for polygon creation, resize, move, and boolean
  operations, not only equivalent GED commands;
- cover point/rectangle selection and modifiers, hierarchy-scale tree row
  styling, selection/highlight persistence across draw/erase/edit,
  lighting/navigation/axes/grid/framebuffer faceplate control, single/quad
  layouts, resize, and fractional DPR; and
- prove GUI widgets, command readback, retained manipulators, and database
  state agree regardless of which surface made the last change.

The deterministic measurement gate now covers 2D and exact-hit 3D mouse
gestures, distance/angle and degree/radian readback, cancellation, tool
replacement, resize, single/quad layouts, and System GL/OSMesa presentation.
It remains part of fractional-DPR hardware repetition, but is no longer an
unimplemented interaction mode.

The deterministic view-settings gate now covers bidirectional GED/widget
state for ADC, center dot, grid, model/view axes, scale, parameter/FPS text,
framebuffer composition, and the cutting plane in single/quad layouts on both
renderers.  It verifies retained Obol feature presence, command readback,
visible clipping, and exact restoration.  Fractional-DPR repetition and the
interactive in-scene cutting-plane affordance remain; ordinary settings
readback is no longer an unqualified control path.

The existing edit runtime and deterministic GUI replays are regression
substrate, not a substitute for physical-pointer, plugin-lifecycle,
hierarchy-scale, and renderer-backed qualification.

## P1: scale, cache, and geometry completion

- Measure discovery on cold local storage, warm page cache, and slower storage.
  Parallel read-only discovery must be bounded by I/O and memory and publish
  ordinary librt/GED hierarchy data without a parallel object ecosystem.
- Complete xpushed multi-Lucy cold preparation and visibility-turnover tests.
  Shared-instance reuse and thousands of genuinely distinct assets are
  separate workloads and both are required.
- Qualify cold spatial previews as globally representative, budget-aware
  presentations.  A local source page may never replace whole-object coverage.
- Preserve the cross-renderer NIST adaptive-wire regression and the passing
  shaded NIST growth/zoom-out/reclamation/cache-restore lifecycle.  Repeat the
  same constrained-residency qualification on at least one large real BREP;
  the partial Big Boy BoT is not evidence for its original BREP hierarchy.
  The current indexed-face guard rejects the known partial Big Boy tire and
  cache version 3 prevents its replay, while fresh-cache NIST remains green.
  This is only fail-closed containment.  Replace the legacy aggregate-success
  CDT call with one bounded provider contract carrying deadline,
  memory/result limits, cancellation, per-face completion, and typed outcome.
  Both LoD-on and LoD-off shaded drawing must consume that same validated mesh
  before the original Big Boy hierarchy can become a release gate.
- Verify constrained-memory cache reload, corrupt/incompatible cache
  invalidation, background compaction, and cancellation during persistence.
- Qualify simultaneous true-cold qged processes against the same large-asset
  cache.  LMDB transaction publication and the final completeness witness have
  passed the Generic Twin cross-process correctness gate, but independent
  processes still duplicate source classification until one publishes that
  witness.  Measure the memory/CPU cost on a large unique asset and add a
  crash-recoverable per-content generation lease if the stampede can violate
  the working-set contract.
- Measure frame, publication, cache, and resident-memory costs at 50k/150k and
  multi-gigabyte scale.  Preserve essential summaries and shared fixtures;
  remove disposable captures and duplicate large directories after each run.
- Extend the passing heterogeneous 10k OBB cold/warm matrix to no-LoD
  ground-truth comparisons on realistic wheels, blades, booms, and hulls and
  to the 50k/150k pressure range.  Persistence and cross-renderer selection are
  already qualified: the constrained OSMesa pair deterministically retains
  693 OBBs across cold-manifest and warm-manifest runs.  Record the perceptual
  improvement over AABB and the cache-generation cost on real geometry.
- Decide whether a first-cold discovery manifest should be enriched after
  detached PoP characterization.  The live payload already carries its OBB;
  the sealed discovery journal and its cache lifecycle lock deliberately make
  the initial structural cue an AABB.  Any enrichment must remain bounded and
  must not introduce nested draw-cache access or a second occurrence registry.
- If multi-page assets need terminal proxying under measured pressure, add one
  occurrence-level aggregate owner which suppresses/restores all pages
  atomically; never attach the whole-source OBB independently to each page.

## Exit condition

The new stack is production-ready only when the exact final binaries pass the
release matrix on required renderers/platforms, the remaining controller and
source boundaries are understandable and independently testable, real large
models meet the fidelity/responsiveness/memory requirements, and qged, MGED,
gsh, Archer, and rtwizard retain all user-facing drawing and editing behavior.
