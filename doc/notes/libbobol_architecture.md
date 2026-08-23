# BRL-CAD Obol drawing architecture

Last reviewed: 2026-08-22

This is the entry point for the production drawing architecture shared by
qged, MGED, gsh, Archer, rtwizard, and other GED clients.  Intermediate branch
APIs and migration designs are not compatibility contracts.  The user-visible
behavior, the installed GED/libBObol APIs, and the state contracts referenced
below are authoritative.

## Documentation map

Read these documents in this order:

1. This overview defines component ownership and the data flow.
2. `libbobol_lod_state_contract.md` defines progressive-display safety,
   liveness, memory, view, and performance invariants.
3. `ged_draw_api.md`, `libbobol_api_contract.md`, and
   `obol_endpoint_lifecycle.md` define the installed interfaces and lifetime
   rules.
4. `libbobol_platform_threading_matrix.md` defines owner-thread and host
   requirements.
5. `libbobol_active_debt.md` is the only active drawing-system TODO list.
6. `obol_production_readiness.md` is the executable release matrix.
7. `qged_editing.md` is the detailed primitive/sketch editing workstream.
8. `libbobol_engineering_lessons.md` records resolved failures whose causes and
   guards must not be rediscovered.

`ObolHostWork.tla` and `obol_lod_control.tla` are focused control-plane models.
They do not model mesh arrays, rendering, or numeric resource policy.

The older `brl_obol_*`, `obol_*_coverage`, and archived TODO documents are
historical design and migration evidence.  They may explain why a test exists,
but they do not define alternate supported behavior.

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
- selection/picking presentation and retained edit features;
- framebuffer composition, faceplate features, and lighting policy; and
- the controller state machine and host-work contract.

Workers produce immutable data.  Coin/Obol nodes, view state, result adoption,
and live presentation remain on the owning presentation thread.  Large arrays
are shared or moved into their final owner; the owner thread must not copy a
mesh merely to publish it.

### Obol/Coin: retained presentation

Obol owns retained scene nodes, frame planning, renderer batching, GPU/OSMesa
resources, picking primitives, and exact execution of the cut selected by
libBObol.  `SoCADAssembly` does not independently choose CAD LoD or infer
quality from a cut number.  It consumes producer-certified populations and
reports completed-frame cost and resource evidence.

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
must not replace compact records with per-occurrence scene objects.

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
  Worker execution may continue only through bounded queues.  A completed
  successor frame is the commit edge which reopens publication.

- Structural fallback points are published best-available representations,
  not hidden richer geometry.  Their aggregate point threshold is independent
  of the triangle-prefix recovery controller and remains at its last proven
  value until mesh publication replaces the fallback population.  Lowering
  that threshold would reveal boxes and is a visual regression.  Consequently
  a one-pixel static-quality trial is admissible only when no structural
  fallback occurrences remain.

## Formal and executable proof boundaries

The complete renderer is intentionally not modeled formally.  The focused
models cover two concurrency seams that have produced real failures:

- `obol_lod_control.tla`: fallback identity, payload publication, repair
  admission, and render acknowledgement;
- `ObolHostWork.tla`: level-triggered controller work, toolkit notification,
  frame revision capture, and eventual quiescence.

TLC counterexamples must be converted into source-level invariants and focused
tests.  ASan/UBSan, fake-clock coordinator exploration, graphical event replay,
image comparison, perf, and large-model runs cover the data-plane properties
outside those models.

## Change rule

New work should refine these boundaries rather than add a parallel owner.  A
change which needs a compatibility draw path, client-side selection repair,
per-occurrence scene duplication, a second LoD policy, or an unbounded owner-
thread pass requires an explicit architecture review.
