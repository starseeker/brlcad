# libBObol active debt

Last reviewed: 2026-08-27

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
- Progressive spatial pages are storage/render partitions of one logical CAD
  occurrence.  They inherit its transform, style, selection, picking, and
  semantic path and are not promoted to independent CAD objects.
- Workload names and renderer names are qualification profiles, not control
  modes.  Lucy, Hubble, 50k, 150k, System GL, and OSMesa all use the same
  state transition relation.

Focused Obol CAD tests, libBObol rendering/realization/update tests,
`ged_test_obol_draw_sync`, the pivot guard, and the installed-package consumer
pass after the 2026-08-27 publication/API cleanup.  This is a development
gate, not production clearance.

## P0: finish the control-state reduction

The canonical contract is `libbobol_progressive_pipeline_contract.md`: one
five-domain evidence stamp (inventory, availability, view, policy, capacity),
at most one bounded plan cursor, at most one presentation transaction, and a
finite work ledger.  HUD outcome is a projection of those values, never a
second phase machine.

Much of this contract is implemented: typed revision advancement, stale-plan
rejection, finite interaction/static-quality/point-quality states, one
presentation transaction, bounded preparation evidence, and the control
refinement map all have focused tests.  The former 6.1k-line coordinator
header is now a small umbrella over independently compile-checked admission,
capacity, delivery, scene-evidence, presentation, and view-policy boundaries.
The production shape is nevertheless still concentrated in two places:

- `lod_admission_policy_private.h` now retains the allocation-free evidence
  values and planner declarations; the 71-method numeric planner is compiled
  independently in `lod_admission_policy.cpp`.  Point and structural evidence
  remain in the same 1.4k-line header because their small methods preserve the
  trivially-copyable value contract, but any further extraction must keep the
  header-alone compile and exhaustive policy tests.
- `view_lod_coordinator_state_private.h` is about 1.3k lines and still exposes
  many mutable companions around submission and repair work.
- `view_controller.cpp` is about 6.5k lines after lifetime/camera/host-request,
  exact-picking, and residency extraction.  It remains both reducer caller and
  the LoD policy effect executor.  Exact-view quality history is now an
  independently compiled allocation-free policy rather than another inline
  controller implementation.

Required completion work:

1. Inventory every remaining mutable controller field by sole owner, revision
   domain, progress witness, and terminal transition.  Delete write-only or
   derivable fields.
2. Continue moving nontrivial pure-policy method bodies behind the new
   compile-checked private boundaries.  Keep hot occurrence storage dense and
   allocation-free; this is a responsibility split, not an object hierarchy.
3. Route remaining direct lifecycle mutations through one typed reducer which
   returns a complete successor plus bounded effects.  An unmapped transition
   is debt unless it represents a genuinely new semantic event.
4. Add randomized refinement traces and retain the exhaustive focused value
   tests.  Every nonterminal state must identify a worker, cursor, scheduled
   owner-thread pump, frame acknowledgement, or finite timer which can make
   progress.
5. Re-run the applicable TLC models after ownership/liveness changes and the
   graphical matrix after numeric policy changes.  Formal models do not judge
   visual quality or wall-clock performance.

Acceptance: unchanged evidence cannot reopen planning; an invalid/stale plan
cannot commit; one event cannot select two successor owners; no terminal HUD
state has foreground work; and no nonterminal state is ownerless.

## P0: complete physical responsibility extraction

The initial extraction and the Obol renderer split are real improvements, but
several units still obstruct review:

- `database_source.cpp` is about 17.4k lines.  Compact CAD presentation is now
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
  styling, selection/highlight persistence across draw/erase/edit, every
  measurement mode, lighting/navigation/axes/grid/framebuffer faceplate
  control, single/quad layouts, resize, and fractional DPR; and
- prove GUI widgets, command readback, retained manipulators, and database
  state agree regardless of which surface made the last change.

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
- Qualify adaptive BREP wire/shaded tessellation growth and zoom-out
  reclamation on the NIST corpus and at least one large real BREP.
- Verify constrained-memory cache reload, corrupt/incompatible cache
  invalidation, background compaction, and cancellation during persistence.
- Measure frame, publication, cache, and resident-memory costs at 50k/150k and
  multi-gigabyte scale.  Preserve essential summaries and shared fixtures;
  remove disposable captures and duplicate large directories after each run.

## Exit condition

The new stack is production-ready only when the exact final binaries pass the
release matrix on required renderers/platforms, the remaining controller and
source boundaries are understandable and independently testable, real large
models meet the fidelity/responsiveness/memory requirements, and qged, MGED,
gsh, Archer, and rtwizard retain all user-facing drawing and editing behavior.
