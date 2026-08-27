# BRL-CAD Obol drawing architecture

Last reviewed: 2026-08-25

This is the entry point for the production drawing architecture shared by
qged, MGED, gsh, Archer, rtwizard, and other GED clients.  Intermediate branch
APIs and migration designs are not compatibility contracts.  The user-visible
behavior, the installed GED/libBObol APIs, and the state contracts referenced
below are authoritative.

## Documentation map

Read these documents in this order:

1. This overview defines component ownership and the data flow.
2. `libbobol_progressive_pipeline_contract.md` defines the canonical pipeline,
   its sole owners, revision/evidence rule, outcomes, and complexity budget.
3. `libbobol_lod_state_contract.md` defines lower-level progressive-display safety,
   liveness, memory, view, and performance invariants.
4. `ged_draw_api.md`, `libbobol_api_contract.md`, and
   `obol_endpoint_lifecycle.md` define the installed interfaces and lifetime
   rules.
5. `libbobol_platform_threading_matrix.md` defines owner-thread and host
   requirements.
6. `libbobol_active_debt.md` is the only active drawing-system TODO list.
7. `obol_production_readiness.md` is the executable release matrix.
8. `qged_editing.md` is the detailed primitive/sketch editing workstream.
9. `libbobol_engineering_lessons.md` records resolved failures whose causes and
   guards must not be rediscovered.

`ObolHostWork.tla`, `obol_lod_control.tla`, `ObolLodConvergence.tla`,
`ObolLodAdmission.tla`, `ObolLodArbitration.tla`,
`ObolProgressivePipeline.tla`, `ObolInteractionSession.tla`,
`ObolDeadlineOwnership.tla`, and `ObolCapacitySearch.tla` are focused
control-plane models.  They do not model mesh arrays, rendering, numeric
resource policy, or workload cardinality.  `ObolProgressivePipeline` replaces
the earlier mode-composition model: workload labels are model-checking
profiles, never production control modes.

The `brl_obol_*`, `obol_*_coverage`, `obol_draw_perf_debug.txt`, and old TODO
notes are historical design and migration evidence.  Their headers point here
for current authority.  They may explain why a test exists, but they do not
define alternate supported behavior.

## Responsibility boundaries

### libged: semantic CAD intent

Libged owns renderer-neutral user intent and database semantics:

- drawn roots and ancestor/descendant frontier subsumption;
- exact and subtree erase;
- view ownership and independent/shared view scope;
- stable semantic occurrence references;
- selection sets and derived ancestor/impact relations;
- edit-scope acquire/release, conflicts, and database reconciliation;
- typed view features, annotations, polygons, faceplate state, and committed
  scene deltas; and
- database mutation events and client observers.

One semantic reducer commits these operations.  Each commit produces one typed
delta, and one private backend seam applies that delta once.  Libged does not
own realized mesh arrays, mirror the Obol scene, or walk renderer nodes to
repair semantic state.

### libBObol: realization and view policy

LibBObol owns the retained CAD data plane:

- database-source traversal and compact occurrence inventories;
- cold and warm manifest realization;
- bounded source and mesh worker services;
- immutable PoP hierarchy production and cache validation;
- view projection, visibility, screen-error demand, and importance;
- scene-level draw, memory, and interaction budgets;
- resident-prefix growth and compaction;
- page-local immutable presentation of spatial PoP assets without a second
  aggregate mesh copy;
- selection/picking presentation and retained edit features;
- framebuffer composition, faceplate features, and lighting policy; and
- the derived convergence snapshot, atomic presentation transaction, and
  host-work contract.

The progressive control plane is a functional core behind an imperative
controller shell.  Inventory, availability, demand, policy, and capacity facts
are revisioned inputs; the allocation planner returns a complete successor
evidence value and typed decision.  A bounded execution cursor may resume the
mechanical application of a plan, but it does not decide quality or readiness.
The controller is an effect executor, not another policy owner.  The current
implementation enforces the evidence/cursor split and const evidence boundary.
Plans and cursors carry the exact inventory, availability, view, policy, and
capacity revision stamp; typed revision advancement is the sole production
cursor-invalidation path, and stale certified commits are rejected.  See the
canonical pipeline contract for the acceptance gate and change rules.

That boundary is the target architecture, not yet a claim that the complete
controller has reached its final maintainable form.  The current coordinator
still carries overlapping submission, retained-allocation, quality-trial,
point-proxy, interaction-handoff, and recovery latches.  They are historical
implementation debt even when their individual tests pass.  The P0 control-
state reduction in `libbobol_active_debt.md` must replace them with the
canonical evidence snapshot, one bounded plan execution, one presentation
transaction, a finite work ledger, and revision-bound typed certificates.
Workload profiles and renderer backends must not survive that migration as
control modes.

Provider terminality, worker-result readiness/age, inventory-delta coalescing,
and resident-growth obligation have one controller-owned availability ledger.
Per-asset immutable data and residency remain in the LoD service/cache; the
controller ledger records their progress edges rather than copying that asset
state.  One stateless availability scheduler supplies pacing and
settled-population predicates.

Workers produce immutable data.  Coin/Obol nodes, view state, result adoption,
and live presentation remain on the owning presentation thread.  Large arrays
are shared or moved into their final owner; the owner thread must not copy a
mesh merely to publish it.

A spatial PoP page is a storage/rendering partition, not a CAD occurrence.
All presented pages for one occurrence use one active cut and retain the base
occurrence's transform, style, selection, picking, and path semantics.  Warm
cache reads and cold live publication produce the same layer representation.
Unpaged assets use the ordinary one-part cumulative prefix.  A multi-page
asset is never repacked into a monolithic renderer array during refinement or
compaction.

### Obol/Coin: retained presentation

Obol owns retained scene nodes, frame planning, renderer batching, GPU/OSMesa
resources, picking primitives, and exact execution of the cut selected by
libBObol.  `SoCADAssembly` does not independently choose CAD LoD or infer
quality from a cut number.  It consumes producer-certified populations and
reports completed-frame cost and resource evidence.

`SoCADAssembly`'s subpixel classifier also remains semantic: it stores one
logical point record per collapsed occurrence, which preserves exact coverage,
selection, highlighting, picking, and later mesh promotion.  Only the
software renderer may replace that immutable logical stream with a cached,
camera-local screen-bin presentation stream.  Its deterministic representative
priority is selected/hovered/color-emphasized, nearest depth, then stable
instance ID; bins start at native-pixel resolution and grow only enough to
honor the software point cap.  Hardware GL consumes the logical stream
directly.  The physical submitted-point count is diagnostic-only and must
never be used as LoD coverage or selection state.

Discovery publishes an immutable source-profile summary for each completed
compact draw epoch.  That summary is an admission gate only: a small bounded
source population may permit the view allocator to perform a modestly wider
bounded perceptual planning pass and may later permit it to consider direct
full-mesh results,
but projection, subpixel aggregation, visual importance, and scene budgets
remain the sole authority for whether any occurrence receives one.  A large
or unknown profile must retain compact contracts; it must never cause a
source-size-only direct-import storm.

Before that optional scene-level profile exists, source realization already
performs a narrower safety admission: a leaf path derives its outer worker
reservation from the immutable directory record, while a combination reserves
the finite coordinator capacity because its own record cannot describe its
descendant imports.  This is a concurrency/memory guard, not a rendering
policy or a direct-mesh authorization.

Obol may accept application-provided picking or scene extensions through
documented capability seams, but must not depend directly on librt types.

### Clients

qged, MGED, gsh, Archer, and rtwizard create or host display endpoints and
translate toolkit input into GED/view operations.  They do not maintain a
second draw, selection, editing, or LoD authority.  A command, GUI control, or
3D manipulator changes the same state and every client view observes it.

## Private implementation layout

Responsibility-specific translation units are deliberately kept out of unity
builds.  Independent compilation is a dependency test: a private helper must
be declared at the owning boundary rather than becoming visible accidentally
because two implementation files happened to share one compiler invocation.
The first extraction pass is complete, but it is not a claim that every source
is now small: `database_source.cpp` remains about 19k lines and `draw_obol.cpp`
about 13k.  The next split must follow an observed ownership/lifecycle or perf
seam, never line count alone.

- `database_source.cpp` owns traversal and realization orchestration;
  `database_source_cache.cpp`, `database_source_compact_access.cpp`,
  `database_source_interaction.cpp`, `database_source_summary.cpp`, and
  `database_source_types.cpp` own the named subordinate responsibilities.
- `view_controller.cpp` owns LoD demand/allocation policy and controller
  lifetime; `view_controller_progressive.cpp`,
  `view_controller_render.cpp`, `view_controller_scene.cpp`, and
  `view_controller_lighting.cpp` own provider pumping, frame execution,
  forwarding, and lighting respectively.
- `draw_obol.cpp` owns semantic transaction application and progressive source
  realization; `draw_obol_endpoint.cpp`, `draw_obol_geometry.cpp`,
  `draw_obol_overlay.cpp`, and `draw_obol_scene_records.cpp` own the renderer
  boundary surfaces their names describe.
- Obol's `SoCADAssembly` owns retained assembly state and frame planning;
  `CadSoftwareWire` owns the non-GL wire execution path.

The private headers expose only hidden value and lifecycle seams.  They are not
installed APIs.  Further extraction must keep dependency direction one-way and
must not replace compact records with per-occurrence scene objects.  Candidate
seams currently worth measuring are source traversal versus realization
publication, GED transaction reduction versus backend application, and retained
assembly planning versus renderer execution.

For progressive control, physical extraction follows the canonical state
shape rather than the current file boundaries:

- `lod_evidence_private.*` owns the five immutable snapshots and typed
  certificates;
- `lod_planner_private.*` is the pure numeric allocation function and its
  independent test surface;
- `lod_control_private.*` owns the finite event reducer, work ledger, plan
  cursor, and presentation transaction;
- `view_controller_progressive.cpp` is the imperative effect executor and
  provider/host adapter; and
- render units consume committed plans and return transaction evidence only.

These names describe ownership, not a requirement to create empty wrapper
files.  Each extraction must delete old fields and entry points in the same
change.  Temporary mirrored state or a second scheduler is forbidden.

## Identity domains

The asynchronous drawing system keeps these identities separate:

| Domain | Meaning | Owner |
|---|---|---|
| semantic source | drawn database root, mode, and source revision | libged |
| source route | one live producer/controller binding | libBObol |
| source population epoch | one compact occurrence inventory | libBObol |
| immutable asset | content-addressed PoP hierarchy | libBObol LoD service |
| occurrence | path, transform, appearance, and view-local demand | libged/libBObol boundary |
| presentation part | retained geometry and instance records | Obol |
| view/policy epoch | camera and policy snapshot used for a request | libBObol controller |
| frame revision | exact render request captured and acknowledged | host/Obol |

Persistent geometry identity excludes owner-local route and population epochs.
Those epochs authenticate asynchronous publication without fragmenting reusable
cache data.  A dense index is only a hint and must be validated by its fixed-
width population epoch or resolved through the semantic key.

## Progressive data flow

```text
GED command/database event
        |
        v
semantic reducer --committed delta--> private Obol backend
        |                                  |
        |                                  v
        |                         compact source/overview
        |                                  |
        |                         bounded realization workers
        |                                  |
        |                         immutable PoP assets/cache
        |                                  |
        +--------------------------> view demand coordinator
                                           |
                                  owner-thread publication
                                           |
                                      SoCADAssembly
                                           |
                                completed-frame/resource evidence
                                           |
                                  next bounded policy decision
```

Cold startup publishes the cheapest truthful information first: one whole-
target overview, then independently useful leaf coverage, then minimum valid
mesh prefixes.  Richer prefixes extend the same immutable asset.  Warm startup
may skip directly to certified manifest and PoP data, but obeys the same
presentation contract.

Stable does not mean all source triangles are resident.  It means every visible
occurrence has the richest pixel-justified cut affordable under the current
finite frame and memory budgets, all missing work has a named obligation or
terminal reason, and the exact current frame has been acknowledged.

## Hard performance constraints

- No per-leaf libged geometry mirror for a compact draw root.
- No full-scene traversal for a sparse selection, style, erase, or edit delta.
- No GUI-thread mesh copying or rebuilding already resident PoP prefixes.
- No O(scene size times bounded window count) planning pass.
- Source and mesh concurrency is bounded by transient bytes as well as workers.
- Stable framebuffers are retained until a semantic, camera, policy, or
  presentation change actually requests another render.
- Orthographic LoD depends on projected visibility and extent, never depth.
- Rotation and translation begin from retained proven cuts and back off only
  when measured interaction cost requires it.  Zoom retargets from current
  resident levels while input remains active.
- A static view carrying over a coarser-than-one-pixel interaction target uses
  measured frame-time and scene-cost headroom to choose the finest safe
  intermediate target.  Coarse interactive targets are not terminal quality
  tiers; the raster-stable 0.75, 0.5, and 0.25 targets remain separately
  witnessed discrete terminal tiers.
- Non-visible assets may release resident suffixes under memory pressure while
  preserving semantic intent and reusable cache identity.
- Quiet resident demand stops at the allocator-selected presentation cut.
  Finer physical pixel demand remains quality evidence, but cannot repeatedly
  populate suffixes which the current scene allocation cannot present.  Active
  scale interaction may prefetch one bounded transition past a stale allocation
  so a large mesh can refine while zooming; the next allocation remains the
  authority over its visible cut.
- Structural boxes are unresolved/failure fallbacks, not a routine motion LoD.
- Subpixel aggregation is a terminal representation only when exact coverage,
  visibility classification, and the applicable performance or memory limit
  are all proven.
- A point/box classifier threshold is a presentation transaction.  Its exact
  framebuffer acknowledgement releases source admission even when no managed
  mesh exists yet and the frame is therefore ineligible for capacity
  calibration.  Admission may not depend on the population it is responsible
  for creating.
- A hard endpoint abort is capacity evidence outside completed-frame timing
  jitter.  At the irreducible PoP prefix it must advance the reversible point-
  aggregation bracket.  A rejected one-pixel static trial remains a terminal
  capacity witness for its view/capacity epoch; recovery threshold changes do
  not reopen the identical rejected population.
- A scene-wide importance/minimax allocation is one closed-population
  transaction.  Provider-inventory completion does not close it while mesh
  tasks, queued/results publication, or the coalesced resident-suffix drain
  can still enlarge the occurrence population.  Keep the last coherent
  framebuffer during that drain and allocate once on its terminal edge.
- A failed headroom presentation is a retained negative witness keyed by
  view, policy, population, and attempted allowance.  Cancelling and
  forgetting that witness permits the same rich/coarse cycle; only a changed
  population, new epoch, or materially larger proven allowance may retry.
- An interrupted retained traversal closes the presentation transaction.
  Owner-thread scene publication cannot mutate the population while the
  renderer holds a resumable classifier or command-plan cursor; otherwise
  every bounded retry invalidates the work completed by its predecessor.
  The cursor belongs to one exact target and reports finite total/remaining
  work plus its admitted transient reservation.  Only strict same-target
  progress permits another retry; an activity serial does not.  Worker
  execution may continue only through bounded queues.  A completed successor
  frame is the commit edge which reopens publication.

- Structural fallback points are published best-available representations,
  not hidden richer geometry.  Their aggregate point threshold is independent
  of the triangle-prefix recovery controller and remains at its last proven
  value until mesh publication replaces the fallback population.  Lowering
  that threshold would reveal boxes and is a visual regression.  Consequently
  a one-pixel static-quality trial is admissible only when no structural
  fallback occurrences remain.

## Formal and executable proof boundaries

The complete renderer is intentionally not modeled formally.  The focused
models cover three control seams that have produced real failures:

- `obol_lod_control.tla`: fallback identity, payload publication, repair
  admission, and render acknowledgement;
- `ObolHostWork.tla`: level-triggered controller work, idempotent pending
  requests, one-way presentation-to-capacity upgrade, pre-traversal transaction
  claims, and eventual quiescence.  TLC passed the current contract on
  2026-08-25 (2,007 generated / 714 distinct states, depth 15);
- `ObolLodConvergence.tla`: exclusive progress ownership across bounded
  coverage, resident loading, allocation, result handoff, and frame
  acknowledgement for both a large single mesh and a many-occurrence scene.
  It proves eventual quiescence after finite input, including the valid case
  where a constrained terminal allocation retains known physical-quality debt.
- `ObolLodAdmission.tla`: safe-scene direct admission, large-scene PoP-only
  admission, subpixel aggregation, and an explicit constrained outcome.
- `ObolLodArbitration.tla`: protected visual-importance floors, budget-safe
  arbitration, unchanged-epoch cut monotonicity, and terminal constrained
  debt.  Pair it with `test_bobol_retained_allocation_oracle` when changing
  numeric allocation scoring.
- `ObolProgressivePipeline.tla`: the canonical ownership composition.  It
  treats few-small, many-small, single-large, distinct multi-large, and
  shared-large scenes as profiles of one discovery, availability, planning,
  presentation, and outcome pipeline.  It guards against a workload-specific
  control mode, repeated planning of identical evidence, ownerless active
  work, false readiness, and starvation of one of several large assets.
- `ObolPresentationPreparation.tla`: the renderer/host preparation seam.  An
  exact target owns an immutable finite total, monotone completed work, and a
  bounded transient reservation.  Only new-target admission, strict
  same-target progress, completion, constraint, or failure is an abstract
  event; renderer activity counters are outside the control refinement map.
  Obol publishes this certificate for subpixel classification, flat shaded
  atlas construction, and retained-indirect preparation.  libBObol compares
  frame-boundary snapshots and has no raw-serial retry path.
- `ObolLodColdPreview.tla`: the cold serialized-large-mesh seam.  It
  distinguishes a cheap source-order sample from a coverage-certified preview,
  requires the standing overview until that certificate exists, bounds producer
  admission by declared working set, and proves a finite view cannot remain in
  an ownerless or undisplayable intermediate state.  It intentionally cannot
  establish whether a particular proxy looks useful or completes within a
  number of milliseconds; those remain image, perf, and memory qualification.
  The implementation counterpart is an explicitly typed cache callback:
  complete PoP prefixes and source-wide coverage previews are distinct payload
  kinds.  Coverage samples become bounded voxel surface geometry at the
  service/presentation boundary.  A former partial spatial-leaf callback was
  removed because no consumer could present it as whole-object coverage
  without violating that contract.
- `ObolLiveSpatialPublication.tla`: the additive spatial-page seam used by
  both cold production and warm cache reads.  Completed immutable pages are
  adopted as renderer layers without constructing a scene-wide aggregate;
  every visible page uses one common active cut.  The model proves the
  required publication ordering: a page is complete before publication,
  coverage remains through every incomplete page set, cancellation freezes
  new publication and cannot mark the cache complete, and only final complete
  geometry may retire coverage.  It is a control-plane guard; retained-page
  tests and graphical qualification remain the data-plane evidence.
- `ObolActiveProducerDemand.tla`: the stable-asset/changing-demand seam for a
  coalesced producer.  It checks that intermediate and final publications are
  stamped from current demand, an overtaken completion is superseded rather
  than a provider failure, exact-demand failure state is retired on epoch advance,
  and finite input reaches either a current presentation or a genuine current
  terminal failure.

TLC counterexamples must be converted into source-level invariants and focused
tests.  ASan/UBSan, fake-clock coordinator exploration, graphical event replay,
image comparison, perf, and large-model runs cover the data-plane properties
outside those models.

Run the focused LoD models after changing their respective policies:

```
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolLodConvergence.cfg doc/notes/ObolLodConvergence.tla
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolLodAdmission.cfg doc/notes/ObolLodAdmission.tla
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolLodArbitration.cfg doc/notes/ObolLodArbitration.tla
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolProgressivePipeline.cfg \
  doc/notes/ObolProgressivePipeline.tla
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolLodColdPreview.cfg \
  doc/notes/ObolLodColdPreview.tla
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolLiveSpatialPublication.cfg \
  doc/notes/ObolLiveSpatialPublication.tla
java -XX:+UseParallelGC -jar /home/cyapp/tla+/tla2tools.jar -workers 1 \
  -config doc/notes/ObolActiveProducerDemand.cfg \
  doc/notes/ObolActiveProducerDemand.tla
```

## Change rule

New work should refine these boundaries rather than add a parallel owner.  A
change which needs a compatibility draw path, client-side selection repair,
per-occurrence scene duplication, a second LoD policy, or an unbounded owner-
thread pass requires an explicit architecture review.
