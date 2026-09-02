# BRL-CAD Obol drawing architecture

Last reviewed: 2026-08-29

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
7. `obol_lod_visual_quality.md` defines the image-comparison method and the
   user-visible quality contract.
8. `obol_production_readiness.md` is the executable release matrix.
9. `qged_editing.md` is the detailed primitive/sketch editing workstream.
10. `libbobol_engineering_lessons.md` records resolved failures whose causes and
   guards must not be rediscovered.
11. `libbobol_formal_models.md` maps every TLA+ boundary to its sole production
    owner and defines when each model must change.

`ObolHostWork.tla`, `obol_lod_control.tla`, `ObolLodConvergence.tla`,
`ObolLodAdmission.tla`, `ObolLodArbitration.tla`,
`ObolProgressivePipeline.tla`, `ObolLodComposition.tla`,
`ObolAssetPublicationComposition.tla`, `ObolCadFrameComposition.tla`,
`ObolInteractionSession.tla`,
`ObolDeadlineOwnership.tla`, and `ObolCapacitySearch.tla` are focused
control-plane models.  They do not model mesh arrays, rendering, numeric
resource policy, or workload cardinality.  `ObolProgressivePipeline` replaces
the earlier mode-composition model: workload labels are model-checking
profiles, never production control modes.

`ObolLodComposition.tla` is the bounded seam check which the canonical model
formerly described only informally.  It composes admission, retained growth,
exact presentation, capacity, structural repair, and point quality through one
owner and terminal contract.  Passing focused models without this seam check
is not accepted as evidence that their assumptions are mutually compatible.

The adjacent composition models continue that proof boundary through the data
plane without modeling geometry arrays: `ObolAssetPublicationComposition`
checks superseding demand, additive pages, final hierarchy, and durable cache
publication; `ObolCadFrameComposition` checks atomic scene mutation, exact
renderer preparation, host frame ownership, and report acceptance.

`ObolCadViewPublication.tla` isolates the renderer boundary: preparation and
completed-frame evidence carry an exact view identity, view changes cancel
old preparation, and a consumer cannot accept a stale historical report as a
current frame.  Geometry sharing is intentionally outside that model because
it does not transfer presentation ownership.

`ObolCadMutation.tla` isolates retained-scene mutation: invalid or resource-
denied candidates cannot open an update or change the last valid scene, a
valid candidate remains invisible while its RAII update is open, bounded
mechanical sub-batches drain before close, and a stable input eventually
leaves no staged candidate or open update.  The checked model currently
reaches 128,904 distinct states at depth five.  `ObolInteractionSession.tla` also
includes the terminal ownerless-handoff case and reaches 132 distinct states
at depth eleven.

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
controller shell.  Inventory, availability, visibility, demand, policy, and
capacity facts
are revisioned inputs; the allocation planner returns a complete successor
evidence value and typed decision.  A bounded execution cursor may resume the
mechanical application of a plan, but it does not decide quality or readiness.
The controller is an effect executor, not another policy owner.  The current
implementation enforces the evidence/cursor split and const evidence boundary.
Plans and cursors carry the exact inventory, availability, visibility, view,
policy, and capacity revision stamp; typed revision advancement is the sole production
cursor-invalidation path, and stale certified commits are rejected.  See the
canonical pipeline contract for the acceptance gate and change rules.

That boundary is now physically represented by cohesive typed values for
submission, retained allocation, quality trials, point recovery, interaction,
presentation handoff, host render requests, and structural repair.  These
values compose through one finite work-ledger projection; they are not wrapped
in a monolithic lifecycle object.  Such a facade would add a large effect
vocabulary without eliminating an owner, obscure independent concurrency, and
invite copying of scale-sensitive state.  Remaining control debt is semantic-
event trace coverage and further effect-shell extraction, not another phase
machine.  Workload profiles and renderer backends remain qualification inputs,
never control modes.

Provider terminality, worker-result readiness/age, inventory-delta coalescing,
and resident-growth obligation have one controller-owned availability ledger.
Per-asset immutable data and residency remain in the LoD service/cache; the
controller ledger records their progress edges rather than copying that asset
state.  One stateless availability scheduler supplies pacing and
settled-population predicates.

Workers produce mutable `PartGeometryBuilder` values which are validated and
admitted before leaving the producer boundary.  Admission constructs the only
renderer-visible `PartGeometry` form: a non-constructible, non-copyable
snapshot whose members are const.  Moving a builder transfers its large arrays
without copying; admitting an intentional lvalue makes an independent copy.
Cached admitted snapshots recover a lightweight validation token in O(1) and
are never rescanned.  Coin/Obol nodes, view state, result adoption, and live
presentation remain on the owning presentation thread.

Publication is a transaction boundary.  The owner thread stages geometry
tokens, instance records, and semantic records without mutating the live
assembly.  Whole-journal pure validation completes before opening the update
window.  Geometry and removals commit once, while large occurrence/style/cut
journals use bounded mechanical sub-batches under that same observer-visible
RAII update window; this avoids a second 50k/150k scene vector.  A preflight
rejection names both field and reason and leaves the preceding scene unchanged.
Compact presentation bookkeeping uses a bounded copy-on-write overlay and is
adopted only after Obol accepts the journal.  Every mutation result is consumed
through the centralized libBObol publication boundary.  Early returns and
nested batches cannot leave notification or rebuild deferral active.

Complete-scene replacement additionally constructs its retained candidate
before a no-throw swap and reports resource denial without clearing the live
scene.  Sparse publication intentionally does not clone the complete scene:
its strong guarantee covers validation and off-scene staging, while the
bounded mechanical commit relies on manifest-driven capacity reservation.
Treating allocator exhaustion during that owner-thread commit as a recoverable
transaction error would require a second retained-scene population and is not
part of the scale contract.

A spatial PoP page is a storage/rendering partition, not a CAD occurrence.
All presented pages for one occurrence use one active cut and retain the base
occurrence's transform, style, selection, picking, and path semantics.  Warm
cache reads and cold live publication produce the same layer representation.
Unpaged assets use the ordinary one-part cumulative prefix.  A multi-page
asset is never repacked into a monolithic renderer array during refinement or
compaction.

Aggregate proxy shape is immutable asset metadata, not another LoD result or
control state.  Cache generation accumulates PCA moments during its existing
bounded source scan and, after the first coverage preview is eligible for
publication, performs one sequential projection pass.  Cache format 23 stores
the binary-XYZ corners only when the PCA box is materially tighter than the
source AABB.  Cache reads validate finiteness, orthogonality, corner topology,
and conservative axis-aligned extent before admitting the metadata.  A miss,
invalid record, or insignificant improvement simply retains the AABB.

The draw cache is a carrier, not a second geometry authority.  LoD-asset
record format 3 stores optional canonical-asset OBB metadata, and chunked
manifest format 3 stores optional per-occurrence object-coordinate corners
without allocating eight points for every AABB-only occurrence.  Manifest
identity `leaf-v4` and chunk identity `manifest-chunk-v3` deliberately
invalidate the intermediate leaf-semantics records; there is no compatibility
reader.  The
directory-wide draw format remains 5 because mesh LoD and draw records share
the same cache root: independently versioned draw metadata must not erase a
valid multi-gigabyte PoP cache.  A manifest descriptor is published only after
all bounded chunks, and readers stream those chunks without reassembling a
whole-scene vector.  A first-cold manifest may be sealed before detached PoP
characterization produces the OBB.  That session's live payload still carries
the immutable corners, but the earliest structural coverage remains an AABB
unless the hierarchy was already cached.

Warm-start completeness has two independent witnesses.  A complete census
proves every semantic occurrence, transform, material, Boolean role, and bound
is known, but does not imply every terminal representation is resident.  A
mesh record is already sufficient for lazy view-driven PoP realization;
analytic and BREP records without such a source contract are replayed through
the terminal converter only.  Representation completion is published after
that bounded subset succeeds.  A rejected or stale subset clears the census
witness before falling back to ordinary discovery, so a partial replay cannot
be persisted as authoritative.

Whole-target framing uses librt's read-only `rt_display_bounds` traversal.
It evaluates combination operators from transformed primitive bounds without
constructing regions, soltabs, BVHs, or other raytracing preparation.  The
result is conservative for subtraction and intersection and is an immutable
display-extent certificate; it is not a claim that evaluated CSG tightly
occupies every point in the box.

An admitted monolithic `PartGeometry` carries those optional corners beside
the ordinary mesh; it does not replace or duplicate the mesh.  Obol consults
them only when its existing view/load policy has already selected a batched
box representation.  Wire mode emits 12 edges and shaded mode emits 12 lit
triangles from the same transformed corners.  A startup AABB remains a
temporary structural-coverage result, whereas an OBB used for an admitted
mesh is a renderer representation of that same semantic occurrence and does
not reopen discovery, availability, allocation, or convergence.  Spatial
pages deliberately do not each inherit a whole-source OBB: persistent proxying
of a multi-page occurrence requires one group-level aggregate owner which can
atomically suppress or restore all pages.

Not every progressive part requires a service/cache payload.  A source may
publish an immutable, producer-certified progressive part directly, as the
BREP wire producer does.  Such a part participates in the same view demand,
allocation, render-cost, convergence, and diagnostic domains as a service-
backed mesh, but changing its cut updates only the view-local instance record;
it performs no worker, cache, array-copy, or part-publication work.  Direct cut
bindings are authenticated by source routing ID, population epoch, occurrence
key, and immutable geometry revision.  The source LoD delta journal retires
replaced bindings in proportion to changed entries; a lost journal or replaced
population triggers one scan of the view's direct bindings, never a routine
whole-scene scan.

### Obol/Coin: retained presentation

Obol owns retained scene nodes, frame planning, renderer batching, GPU/OSMesa
resources, picking primitives, and exact execution of the cut selected by
libBObol.  `SoCADAssembly` does not independently choose CAD LoD or infer
quality from a cut number.  It consumes producer-certified populations and
reports completed-frame cost and resource evidence.

One `SoCADAssembly` belongs to one view presentation.  Independent views own
independent assembly nodes and camera-local preparation/report state.  Their
assemblies share immutable `PartGeometry` storage rather than sharing a node;
this prevents one view from invalidating another view's classifier, prepared
commands, or completed-frame certificate without copying mesh arrays.

`SoCADAssembly`'s subpixel classifier also remains semantic: it stores one
logical aggregate record per collapsed occurrence, which preserves exact
coverage, selection, highlighting, picking, and later mesh promotion.  A
genuinely tiny footprint becomes a point; a screen-significant collapsed
occurrence remains one batched AABB or producer-certified OBB.  Only the
software renderer may replace the immutable logical point stream with a
cached, camera-local screen-bin presentation stream.  Its deterministic
representative priority is selected/hovered/color-emphasized, nearest depth,
then stable instance ID; bins start at native-pixel resolution and grow only
enough to honor the software point cap.  Hardware GL consumes logical points
directly, and both renderers batch boxes.  Physical submitted primitives are
diagnostic-only and must never be used as LoD coverage or selection state.

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
is now small: `database_source.cpp` remains about 17.4k lines and
`draw_obol.cpp` about 12.9k.  The next split must follow an observed
ownership/lifecycle or performance seam, never line count alone.

- `database_source.cpp` owns traversal and realization orchestration;
  `database_source_cache.cpp`, `database_source_compact_access.cpp`,
  `database_source_interaction.cpp`, `database_source_presentation.cpp`,
  `database_source_summary.cpp`, and `database_source_types.cpp` own the named
  subordinate responsibilities.  `database_source_mesh_geometry.cpp` and
  `cad_normals.cpp` own worker-safe geometry conversion and normal
  construction.  `compact_presentation_staging.cpp` owns the bounded
  copy-on-write retained-scene journal.
- `view_controller.cpp` owns LoD demand/allocation policy;
  `view_controller_host.cpp` owns controller lifetime, camera synchronization,
  and host render requests; `view_controller_progressive.cpp`,
  `view_controller_render.cpp`, `view_controller_scene.cpp`,
  `view_controller_lighting.cpp`, `view_controller_picking.cpp`, and
  `view_controller_residency.cpp` own provider pumping, frame execution,
  forwarding, lighting, exact picking, and memory reclamation respectively.
  `lod_admission_policy.cpp`, `lod_presentation_policy.cpp`, and
  `lod_view_quality_history.cpp` own numeric admission planning, the
  allocation-free handoff, and exact-view history policies.  The formerly
  monolithic coordinator header is now an umbrella over self-contained
  admission, capacity, delivery, scene-evidence, presentation, and view-policy
  boundaries; each focused header is compiled alone by the policy test target.
- `draw_obol.cpp` owns semantic transaction application and progressive source
  realization; `draw_obol_endpoint.cpp`, `draw_obol_geometry.cpp`,
  `draw_obol_overlay.cpp`, and `draw_obol_scene_records.cpp` own the renderer
  boundary surfaces their names describe.
- Obol's presentation-local `SoCADAssembly.cpp` owns the Coin node, retained
  mutation, picking, and action surface while immutable part geometry may be
  shared across views.  `CadAssemblyPlan.cpp` owns retained plan/cache
  maintenance, `CadAssemblyClassification.cpp` owns camera-local subpixel
  classification, and `CadAssemblyImpl.h` is their state-only private seam.
  `CadRendererGLIndirect`, `CadRendererGLInstanced`, `CadRendererGLFlat`, and
  the retained/direct executor unit own distinct GL execution families;
  `CadSoftwareWire` owns the non-GL wire execution path.

The private headers expose only hidden value and lifecycle seams.  They are not
installed APIs.  Further extraction must keep dependency direction one-way and
must not replace compact records with per-occurrence scene objects.  Candidate
seams currently worth measuring are source traversal versus realization
publication, GED transaction reduction versus backend application, and retained
assembly planning versus renderer execution.

For progressive control, physical extraction follows the canonical state
shape rather than workload or renderer names:

- `lod_scene_evidence_private.h` owns coverage, source-profile,
  projected-demand, visibility, resource, and convergence evidence;
- `lod_capacity_policy_private.h` owns frame capacity, the bounded admission
  cursor, and retained-allocation requests;
- `lod_view_policy_private.h` owns camera interaction and static-quality
  policy;
- `lod_admission_policy_private.h` owns point/structural admission and the pure
  admission planner;
- `lod_delivery_policy_private.h` owns availability, presentation, and
  compaction scheduling;
- `lod_presentation_policy.cpp` owns renderer-ceiling handoff and
  reconciliation;
- `lod_control_private.h` owns the finite event reducer, work ledger, plan
  cursor, and typed obligations;
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
  so a large mesh can refine while zooming.  A changed resident-capacity epoch
  may likewise probe the current physical demand past the stale allocation:
  consuming that retry while still clamping demand to the allocation would be a
  no-op and could strand a coarse cut forever.  In both cases the service byte
  governor remains authoritative, and the next scene allocation owns the
  visible presentation cut.
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

The complete renderer is intentionally not modeled formally.  The
authoritative model catalog, implementation ownership map, change routing,
and TLC invocation rules are in `libbobol_formal_models.md`.  Keep detailed
transition commentary with the focused model rather than duplicating it here.

Every counterexample must become a source-level invariant and an executable
regression test.  ASan/UBSan, fake-clock coordinator exploration, graphical
event replay, image comparison, perf, and large-model runs cover visual,
numeric, memory, and timing properties outside the formal abstraction.

## Change rule

New work should refine these boundaries rather than add a parallel owner.  A
change which needs a compatibility draw path, client-side selection repair,
per-occurrence scene duplication, a second LoD policy, or an unbounded owner-
thread pass requires an explicit architecture review.
