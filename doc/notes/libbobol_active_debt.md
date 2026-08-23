# libBObol active debt

Last reviewed: 2026-08-23

This is the authoritative remaining-work list for the Obol drawing stack.
`obol_production_readiness.md` defines the release matrix; `qged_editing.md`
defines the editing workstream; `libbobol_engineering_lessons.md` preserves
resolved failures.  Do not add chronological test logs or alternate designs
here.

## Current baseline

- The semantic drawing API, compact occurrence model, shared-scene/view-local
  split, progressive PoP cache, retained renderer, selection, polygons,
  faceplate, navigation gizmo, and edit-session substrate are implemented.
- The current binary passed `drawing_baseline` plus `bobol_headless`: 86/86
  tests in 408 seconds, with the expected unavailable X11 capability skip.
- Focused TLC models for payload publication, host work, and Lucy/large-scene
  convergence pass.  The convergence model establishes that physical quality
  debt is not itself a quiet residency-work obligation.
- Current OSMesa evidence includes successful Lucy smooth zoom and certified-
  warm 50k shaded plus 150k shaded/wire full lifecycles.  They finish with
  exact occurrence coverage, zero structural boxes, no pending work, and an
  explicit software-performance limit where appropriate.

This is not release clearance.  The source changed after portions of the older
matrix, System GL cannot be exercised on this host, and visual-importance
quality remains insufficiently qualified on real large models.

## P0: visual quality and progressive convergence

- Redesign and qualify the significance frontier under a constrained large
  scene.  It must order optional transitions by projected error reduction,
  footprint/silhouette significance, user emphasis, complete transition cost,
  and the current frame budget.  Prominent wheels, blades, tails, and hulls
  must not remain conspicuously coarse merely because many cheap objects exist.
- Build two complementary small-scene proofs before changing allocator policy:
  a TLA+ arbitration model for coverage priority, affordable prominent-floor
  non-starvation, budget safety, unchanged-epoch monotonicity, and terminal
  constrained debt; and an exact enumerating or SMT oracle for the numeric
  greedy/minimax result.  Images remain the authority for perceptual quality.
- Preserve the resident-allocation boundary.  Quiet residency stops at the
  allocator-selected presentation cut; active scale interaction may prefetch
  one bounded transition past it.  Reintroducing unconstrained quiet prefetch
  recreates the 150k allocation loop; disabling it during scale interaction
  strands Lucy.
- Preserve the closed-population transaction: do not allocate a scene-wide
  population while provider, task, result, publication, residency-drain, or
  interrupted-render work can still change it.  A terminal capacity witness
  must not rearm without a view, policy, population, or capacity edge.
- Establish explicit release thresholds for first useful image, interaction
  latency, static convergence, stable-frame time, peak/resident memory, and
  visual-quality-floor debt.  Record p50/p95/max rather than one favorable run.

## P0: release qualification

- Run the exact-current graphical matrix from truly cold and certified warm
  caches.  Required OSMesa gaps after the latest residency policy are 50k wire
  and the cold scale rows; re-run representative real-model rows if their
  cache or policy path changes.
- Run the complete System GL matrix on a usable X host.  Do not promote older
  System GL captures to current-binary evidence.  Compare semantics and images
  with OSMesa within declared tolerances; use APNG and apitrace for a flash,
  corruption, camera jump, or GL-state mismatch.
- Qualify Generic Twin, Lucy, multi-Lucy/xpush, Hubble, Havoc, NIST BREPs,
  Stanford meshes, and independent multi-gigabyte vehicle models.  Exercise
  shaded/wire/evaluated modes as applicable, LoD on/off, resize, zoom,
  rotate, translate, selection, exact/subpath erase-redraw, and cold/warm
  cache behavior.
- Run ASan/UBSan across the shared dynamic stack with active worker teardown,
  cache corruption, endpoint replacement, plugin reload, edit cancellation,
  and rapid view close.  Run TSan/LSan on a compatible native worker.  A
  sandbox runtime limitation is not a successful race/leak check.

## P0: interaction and editing qualification

- Drive actual Qt controls and mouse paths, not only equivalent GED commands.
  For every interaction, widget state, GED semantic state, retained scene,
  and command readback must agree after GUI- and command-originated changes.
- Complete the qged selection, polygon, measurement, and faceplate matrix:
  point/rectangle selection and modifiers; tree styling and draw/erase under
  selection; polygon creation/manipulation/booleans; all measurement modes;
  lighting/navigation/axes/grid/framebuffer controls; resize and fractional
  DPR in single and quad layouts.
- Complete the primitive/sketch work described in `qged_editing.md`: classify
  every librt operation, audit rejection/readback, finish the reusable
  manipulator vocabulary and specialized sketch curves, and qualify common
  qged/MGED/gsh edit sessions including sed/oed.

## P1: scale, memory, and maintainability

- Measure discovery on cold local storage, warm page cache, and slower storage.
  Keep parallel read-only discovery bounded by I/O throughput and memory, and
  prove ordinary librt/GED hierarchy and tree publication parity.
- Qualify repeated-instance visibility turnover separately from distinct-asset
  stress.  Verify spatial-page reuse, zoom-out compaction, and cache reload
  without boxes or full asset rebuilds.
- Qualify adaptive BREP tessellation/LoD growth and zoom-out reclamation over
  the NIST corpus and at least one large real BREP model.
- Keep extraction ownership-led.  The first extraction pass exists, but
  `database_source.cpp` (about 19k lines) and `draw_obol.cpp` (about 13k) are
  still concentrated.  Extract only at a demonstrated lifecycle or profiling
  seam; preserve compact storage, sparse deltas, bounded owner-thread windows,
  and immutable payload ownership.
- Replace the remaining monolithic static `ged_test_obol_draw_sync` coverage
  with owning-library tests plus a production shared-stack GED integration
  test.  Static linkage belongs only in explicit static-link checks.

## Exit condition

The stack is production-ready only when the exact final binaries pass the
release matrix on required renderers and platforms, no visual correctness or
responsiveness blocker remains, real large-model quality meets the agreed
significance floor, memory recovers after view turnover, and qged/MGED/gsh
editing and interaction remain coherent under lifecycle stress.
