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
- The current binary passed `drawing_baseline` plus `bobol_headless`: 87/87
  tests in 395.5 seconds on 2026-08-23, with the expected unavailable X11
  capability skip.
- Focused TLC models for payload publication, host work, and Lucy/large-scene
  convergence pass.  The convergence model establishes that physical quality
  debt is not itself a quiet residency-work obligation.
- Current OSMesa evidence includes successful Lucy smooth zoom and certified-
  warm 50k shaded plus 150k shaded/wire full lifecycles.  They finish with
  exact occurrence coverage, zero structural boxes, no pending work, and an
  explicit software-performance limit where appropriate.
- Current linked-Obol validation also found and fixed an OSMesa retained-wire
  upload overread when a richer GPU prefix survived a coarser request.  The
  150k wire cold/warm replay is now qualified under the bounded diagnostic
  environment.  Discovery now retains one arc
  population and aliases resolved directory names, coverage avoids a second
  whole-file mapping for large databases, and manifest persistence streams
  bounded 8 MiB occurrence chunks before publishing its descriptor, so an
  interrupted write remains a cache miss.  The optional OSMesa flat-wire atlas
  also declines a failed transient allocation and falls back to retained
  rendering rather than throwing through Coin.  Re-run the cold/warm pair to
  establish the resulting end-to-end memory bound and visual result.  The
  compact-presentation path now publishes retained instance records in bounded
  batches rather than creating a second whole-scene update vector; this fixes
  a traced OSMesa repaint `std::bad_alloc` at 150k scale.  Compact discovery
  reserves its lightweight renderer-update population before the first merge,
  avoiding the former 131k-to-262k retained-instance reallocation, while the
  heavyweight occurrence records remain segmented to preserve bounded partial
  progress.  Under the bounded
  8 GiB diagnostic limit, the requested 150k compact population now stops
  cleanly before records are drained, retaining the useful overview/partial
  scene rather than aborting qged.  That limit is below the retained discovery
  footprint of the 3.4 GiB fixture, not another vector-doubling regression.
  Repeat the exact-current cold/warm pairs on adequately provisioned
  production hardware before treating all 150k OSMesa modes as cleared.

  A fresh short cold replay under a 16 GiB address-space cap now discovered
  all 150001 leaves in 20.6 seconds, admitted 664 mesh payloads, and reported
  neither terminal provider failures nor structural boxes.  The constrained
  cache fallback is 512 MiB rather than 4 GiB so it can coexist with a
  multi-gigabyte database under that cap.  This is a regression checkpoint,
  not the required full cold/warm qualification.  The full wire matrix then
  passed cold in 44 seconds and certified warm in 81 seconds, with 150001
  discovered leaves, zero structural boxes or terminal failures, no pending
  work, and an explicit software-performance limit in both terminal frames.
  The corresponding exact-current shaded matrix also passed cold in 44
  seconds and certified warm in 61 seconds with the same coverage, box,
  failure, pending-work, and camera-contract checks.

This is not release clearance.  The current binary passed System GL smoke
coverage for Generic Twin and Lucy, but the complete System GL model matrix
and visual-importance quality on real large models remain insufficiently
qualified.

## P0: visual quality and progressive convergence

- **Resolved 2026-08-23:** cache-independent deferred-autoview ownership.  The
  2026-08-23 OSMesa
  Generic Twin shaded cold/warm real-canvas pair reached the same 709 mesh
  payloads with zero boxes, failures, or pending work, but its `ae 90 0`
  stable camera differed materially by cache state.  Cold deferred autoview
  and warm cached completion must target the same view context and preserve
  the same caller orientation before release.  The retained evidence bundle
  is `/tmp/qged-generic-smoke-20260823-1947` until the corrected replay is
  captured; do not weaken the camera-contract comparator to hide it.  The
  immediately repeated pair passed, including when timing diagnostics slowed
  execution, so treat this as a timing-sensitive synchronization race rather
  than a deterministic cache-format discrepancy.
  The retained failing telemetry identifies the concrete software-canvas
  asymmetry: cold reaches the requested `ae 90 0` camera after its first
  progressive synchronization, while warm remains at the initial camera.
  `QgSW::need_update` previously requested only a framebuffer repaint after
  `diff_hashes()`, unlike `QgGL`, so it did not copy a GED command's camera
  into the Obol controller.  The source now requests `BV_REFRESH_VIEW` before
  presenting.  A rebuilt qged passed the exact cold/warm OSMesa replay in
  `/tmp/qged-generic-camera-fix-20260823-rerun`: both runs recorded the same
  `ae 90 0` orientation and the matrix reported zero failures.  Retain the
  old failure bundle as a regression specimen and preserve this replay in the
  release evidence until the full GUI matrix supersedes it.
  The first fix exposed the same ownership edge in a fast warm Lucy wire
  presentation: the first background frame could be requested without a
  semantic canvas refresh.  `QgSW::present_frame` now includes
  `BV_REFRESH_VIEW`, which synchronizes the authoritative GED view at the
  presentation boundary without forcing every worker wake to mutate camera
  state.  The validated warm replay in
  `/tmp/qged-lucy-wire-present-fix-20260823` completed in 61 seconds and its
  stable camera exactly matched the preceding cold wire contract.
- `ObolDeferredAutoview.tla` is the companion bounded ownership model for
  that investigation.  It verifies request replacement, same-view framing
  cancellation, orientation preservation, exact-once completion, teardown
  without cross-view writes, and the producer/wake/presentation liveness
  witness.  TLC passed this extension on 2026-08-23 (4,833 generated / 896
  distinct states, depth 11).  It does not validate Qt event ordering or
  floating-point camera values; retain the real qged cold/warm camera replay
  as the implementation-level acceptance test.
- The significance frontier now orders optional transitions by projected error
  reduction, footprint/silhouette significance, user emphasis, complete
  transition cost, and the current frame budget.  `libBObol_lod_update_action`
  exercises one prominent occurrence against sixteen cheap ordinary
  occurrences, proving that the sole affordable richer prefix goes to the
  prominent one.  Still qualify this on constrained real large scenes:
  prominent wheels, blades, tails, and hulls must not remain conspicuously
  coarse merely because many cheap objects exist.
- The bounded protocol proof is in `ObolLodArbitration.tla`: it covers coverage
  priority, affordable prominent-floor non-starvation, budget safety,
  unchanged-epoch monotonicity, and explicit terminal constrained debt.  Run
  TLC with `ObolLodArbitration.cfg` after changing the arbitration protocol.
  The current model passed on 2026-08-23 (102 generated / 72 distinct states,
  depth 12), together with the allocation oracle and LoD update/coordinator
  tests.
- `test_bobol_retained_allocation_oracle` exhaustively checks the bounded
  numeric marginal frontier and deterministic minimax tie handling against an
  independent oracle.  Keep it with the TLA+ model when changing allocator
  policy.  Images remain the authority for perceptual quality.
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

## P0: bounded cold source preparation

- Cold source preparation must be admitted before it calls `rt_db_get_internal`.
  Under a 12 GiB address-space cap, System GL Lucy reached that librt path and
  `bu_calloc` invoked `bu_bomb`; the existing working-set policy permits an
  exceptional oversized task to run alone, which is not a safe guarantee.
  Replace that exception with a source-preparation admission contract that can
  return a constrained, diagnosable result without entering a non-recoverable
  librt allocation path.  It must preserve a useful existing presentation and
  permit retry after a capacity edge.
- The first guard is implemented and verified: the formerly crashing 12 GiB
  System GL Lucy shaded/wire replays now terminate without `bu_bomb`, with one
  constrained occurrence retaining its structural fallback.  Address-space
  refusal is reported as terminal `FALLBACK`, rather than an implementation
  error, so it is diagnosable without poisoning provider-failure accounting.
  It is safety evidence only.  Cold large-mesh preparation still needs a chunked/importable
  source representation so Lucy can provide a useful mesh rather than that
  fallback under this constraint.
- A V5 serialized-source bridge now handles the first safe alternative when
  raw import admission fails.  It borrows the memory-backed database record,
  validates its envelope and indices, converts points in fixed bounded
  blocks, and uses recoverable C++ allocation failure rather than entering
  `rt_db_get_internal`/`bu_bomb`.  It preserves authored winding/culling and
  deliberately declines BoTs with active authored surface normals until their
  optional records are included in the direct-source contract.  The cold
  serialized path is regression-tested after clearing the cache on a reopened
  memory-backed V5 database.  This is a safety/availability bridge, **not**
  completion of the paged source-preparation P0: native point/face arrays and
  the existing global PoP activation structures are still materialized for one
  asset.  Qualify the 12 GiB Lucy route before claiming it improves that case,
  then replace those global arrays with the byte-budgeted page contract below.
- Do not mistake `BObolStagedSourceMesh` for that representation.  The compact
  discovery detail worker deliberately skips authored BoTs, and its staged
  lease is a short handoff window; the later LoD task consequently reimports
  Lucy through `rt_db_get_internal`.  Extend the existing compact producer
  with a bounded authored-BoT page/import path, or persist an equivalent
  coarse PoP prefix before that lease ends.  It must remain view-demanded and
  must not eagerly import every authored BoT during hierarchy discovery.
  The shared V5 serialized-BoT reader now safely exposes fixed vertex/face
  records to the recoverable bridge, but the cache builder still converts them
  into whole native arrays.  The missing seam is therefore a streaming
  cache-builder/source-page contract with an explicit byte budget, not another
  borrowed-view fallback or a longer staged-source lease.
- The page contract must be deliberately narrow and measurable:
  1. a source reader supplies validated indexed point/face records from either
     native arrays or serialized V5 bytes; its native fast path must remain
     pointer-based;
  2. a serialized reader keeps a bounded per-worker decoded vertex/face page
     working set, with no process-lifetime thread-local pages;
  3. cache identity must remain the canonical post-winding point/face stream,
     so a safe source path cannot create a visually identical but redundant
     hierarchy merely because it arrived from serialized bytes;
  4. per-corner authored normals either participate in that reader or decline
     cleanly—never silently alter shaded output; and
  5. the final source lease and all pages are released before quiet residency
     management begins.  Initial source paging need not remove the separate
     global activation/cache-output structures; account for those explicitly
     in admission and reduce them in a subsequent builder pass.
- **Implemented 2026-08-23:** the V5 source reader is now used through PoP
  classification, initial-prefix materialization, spatial chunk construction,
  and cache serialization.  A serialized producer retains one 16,384-record
  point page and face page per active phase/worker; native sources remain
  direct pointer loads.  The test clears the cache, cold-generates from the
  reopened memory-backed V5 record, proves cache-key equality with the native
  path, and checks CW culling.  Active authored surface normals and
  degenerate serialized faces still decline cleanly rather than changing
  shaded semantics.  This removes full native source-point/face arrays from
  the constrained path, but not the global activation/cache-output structures.
- Acceptance for that page contract is a cold V5 Lucy run under the historical
  12 GiB address cap with a mesh preview (not only a box), no `bu_bomb`, a
  bounded measured page working set, stable post-generation source residency,
  and byte-identical cache identity versus the ordinary native source path.
- **Next producer design — spatial PoP leaves:** retain the existing consumer
  contract of one hierarchy with globally stable quantization cuts and
  independently drawable chunks, but stop deriving those chunks from a
  source-wide activation graph.  After the bounded global vertex-bounds scan,
  partition source faces deterministically into fixed spatial leaves.  A leaf
  owns its source-face page list, local vertex remap, per-cut face/point
  counts, bounds, and immutable payload record.  It may be classified and
  persisted independently under a byte reservation; no leaf build may retain
  another leaf's decoded points, indices, or remap.  The hierarchy descriptor
  is published only after its leaf table is complete, while each leaf record
  has its own completeness marker.  View admission may then request visible
  leaves immediately and schedule their refinement independently.

  The first implementation must keep one global quantization schedule and
  use source coordinates at shared boundaries, so neighboring leaves render
  the same snapped position at a given cut.  It must not invent transition
  triangles or an independent per-leaf coordinate system in the first pass.
  That gives crack-free equal-cut rendering and lets the current renderer,
  which already consumes chunk-local records, remain intact.  Neighboring
  unequal-cut stitching and adaptive meshlet error metrics are a subsequent
  enhancement, not a reason to retain global source arrays now.

  Required invariants and tests:
  1. cache identity is the canonical post-winding source stream plus spatial
     partition/version, independent of producer-worker timing;
  2. every valid source face occurs in exactly one leaf; local indices are
     in range and each leaf's terminal record is exact for its face set;
  3. aggregate terminal face count equals the sanitized source count; shared
     vertices may be duplicated across leaves and therefore only have a
     lower-bound aggregate point-count invariant;
  4. a leaf is independently recoverable after interrupted publication; no
     partial leaf is discoverable from a hierarchy descriptor;
  5. cold scheduling can publish a bounded visible leaf mesh without waiting
     for distant leaves, while a warm hierarchy remains view-selective; and
  6. 12 GiB Lucy must prove these facts through System GL visual capture,
     with a bounded peak and no process-wide termination.
- **Producer constraint confirmed 2026-08-23:** chunk serialization no longer
  allocates source-vertex-sized remap tables; each existing page now keeps only
  a leaf-local source-id map.  The page-record writer and its metadata writer
  are now separately callable from a local source-face/activation stream; the
  existing global producer uses that same path and the mesh-cache suite passes.
  Do not enable chunk-only cache publication atop the current
  global builder: a 20 GiB Lucy GUI replay died immediately after global chunk
  construction when page-record writes overlapped the still-live global
  activation/prefix structures.  Page records must instead be emitted by the
  bounded leaf producer, after the prior leaf's scratch is released.  This is
  a peak-memory invariant, not a cache-format preference.
- **Implemented test-gated 2026-08-23:** the serialized source can now run a
  deterministic spatial producer (`BOBOL_LOD_SPATIAL_LEAVES=1`).  It streams
  fixed-width source-face records into anonymous fixed-grid spools, drains one
  bounded page at a time through the local page writer, and publishes only the
  final metadata/marker after all page records succeed.  Its unit fixture
  verifies both a durable completed route and an injected constrained route:
  the completed descriptor has a partition-specific cache key, zero global
  cluster ranges, exact terminal face coverage, independently readable pages,
  survives reopening, and publishes one bounded direct first-page preview;
  the constrained route has nonzero live coverage but remains undiscoverable
  after its generating handle is released.  `ObolSpatialProducer.tla` models
  this producer boundary: direct preview while a cache write is held, atomic
  durable discovery, process-local constrained pages, cancellation, and no
  cache re-entry during a write transaction.  TLC passed its two-leaf bounded
  model on 2026-08-23 (35 generated / 28 distinct states, depth 6), including
  the pre-identity mapped-source bootstrap lifecycle.  It is
  deliberately not the
  default yet: cache identity/versioning, partial live publication after a
  capacity edge, and real Lucy cold/warm qualification remain before it may
  replace the existing useful-prefix path.
  The producer now emits a 4,096-face direct first-page handoff before the
  complete face stream is spooled; durable leaves remain at the ordinary
  65,536-face target.  A System-GL Lucy diagnostic replay reached that handoff
  in 16.1 ms of producer time and the service accepted it (`published=1`),
  eliminating the former cache-lock self-deadlock.  It still missed the
  replay's 1.5-second image checkpoint because task/scene scheduling reached
  the handoff after that checkpoint, and it did not settle inside the
  intentionally short five-second deadline while full classification
  continued.  Treat this as evidence that first-page construction is bounded,
  not as first-useful-image qualification.  Follow-up tracing separated two
  gates: the controller's coverage-only census correctly protects a
  multi-object scene, while serialized generation was hashing and
  bounds-scanning the full source before it could materialize a first page.
  A test-gated sole-source exception now permits the spatial request while
  the census continues, and a nonterminal 4,096-face sampled bootstrap is
  published before global cache identity/bounds preparation.  The 2026-08-23
  System-GL Lucy replay confirms that the service accepts it (`published=1`,
  4,096 scene faces).  An initial whole-model-stride bootstrap took 1.95 s
  because it decoded almost every mapped vertex page; changing it to one
  contiguous source page reduced that work to 9.8 ms.  The remaining
  first-useful-image bottleneck was the roughly 1.3 s structural discovery
  pass, not page caching, mapped point/face decoding, or GUI adoption.  The
  serialized AABB scan is now a bounded parallel reduction only for a large
  individual BoT (at most two scans, four workers each); ordinary many-leaf
  discovery retains its existing outer worker pool.  The repeated real Lucy
  System-GL replay reduced coverage to 845 ms and recorded the accepted
  4,096-face bootstrap by the 0.9-second checkpoint.  Re-qualify cold/warm
  50k and 150k cases before making a broad discovery-performance claim; add a
  visible spatial-frontier representation before enabling the route by
  default.
  A 2026-08-23 20 GiB qged run selected ordinary native preparation, while
  the corresponding 12 GiB process terminated before LoD preparation began;
  neither is evidence for or against the bounded producer.  Do not attribute
  that pre-provider process failure to page classification.  The disposable
  GUI cache bundles were removed after recording this fact.
- **Implemented test-gated constrained publication 2026-08-23:** a page write
  failure can retain a strictly capped 64 MiB immutable live-page set for the
  generated handle.  No name mapping, descriptor, or completeness marker is
  published in that state, so reopening remains a cache miss; the current
  handle can still read its retained pages through the ordinary chunk callback.
  The mesh-cache test injects the write failure and proves nonzero drawable
  coverage, source-limited status, and that the incomplete hierarchy is not
  discoverable after the handle is released.  This is the correct cache
  atomicity boundary.  It still needs view-prioritized first-page selection and
  a structural representation for unbuilt pages before it provides the final
  cold Lucy user experience.
- **Investigated 2026-08-23:** the serialized reader is correct enough to
  classify the full Lucy source under a 20 GiB diagnostic cap: 28,040,274
  faces and 14,020,180 points completed bounds, classification, and attempted
  persistence in 15.9 seconds.  The same route still fails the 12 GiB
  acceptance case before a mesh is published.  This is not an occurrence-route
  rejection or a malformed V5 record.  The current constrained LMDB map retry
  is 512 MiB; a complete Lucy hierarchy exceeds it, and the monolithic cache
  transaction treats persistence failure as source-generation failure.  Do
  not paper over this by raising the virtual map under the address cap or by
  weakening the cache-completeness marker.  Design a bounded cache-publication
  protocol which keeps incomplete records unreachable, allows an already
  materialized minimum prefix to be published as a constrained live asset,
  and records exactly which higher prefixes are unavailable for retry after a
  cache-capacity edge.  Its compatibility with existing header readers must
  be proved by the mesh-cache suite before it replaces the current atomic
  transaction.
- **Implemented but not yet qualified 2026-08-23:** immutable cache records
  now commit independently and the final marker remains the completeness
  witness.  If the durable publication fails, the producer retains its
  view-bounded initial prefix as a live asset, caps that asset's advertised
  hierarchy at the retained cut, and reports it as memory-limited instead of
  publishing a partial name mapping or returning to boxes.  Focused
  mesh-cache and LoD-service tests pass.  Removing a redundant serialized
  face-validation sweep allowed the 20 GiB System GL replay to classify Lucy
  in 4.3 seconds, construct spatial chunks in 13.1 seconds, reject the
  undersized durable cache safely, and visibly retain a detailed stable mesh
  at `/tmp/qged-lucy-paged-fastpath-20g-20260823`.  The 12 GiB replay still
  terminates before classifier completion, so this is **not** the required
  acceptance evidence: reduce the producer's remaining whole-mesh working
  set or add a demand-driven spatial-leaf builder before claiming constrained
  Lucy is fixed.  Keep the source-limited asset contract, but do not use it
  to mask that earlier failure.
- Do not treat a cache-map fallback as the solution.  The cache now has normal,
  constrained, and minimum mapping retries, but its capacity affects reuse,
  not safe raw-BoT materialization.
- Requalify System GL Generic Twin, Lucy, Hubble, and Havoc after this path is
  hardened.  Generic Twin passed cold/warm shaded/wire at 12 GiB; Hubble wire
  materialized payloads but its first-frame checkpoint was empty, while Hubble
  shaded retained a render/refinement-frame liveness witness.  Neither Hubble
  result is release evidence.

  **Progress 2026-08-23:** the Hubble OSMesa shaded cold replay initially
  exposed a test-host liveness gap after zoom: all service queues were empty,
  but the nested `wait_progressive_idle` loop observed an outstanding render
  request without asking the visible endpoint to consume it.  The event driver
  now presents such requests through the same real-canvas path it uses for
  scope barriers; it does not synthesize controller completion.  The rebuilt
  cold replay passed in 17 seconds, including the Hubble prominent-object
  floor assertions.  A separately launched certified-warm replay then passed
  in 19 seconds.  It exposed one validator inconsistency: Hubble's selected
  panel may be depth-occluded, and the comments already allow a zero erase
  delta, but the redraw threshold was still nonzero.  The semantic
  compact-entry/selected-instance checks prove the exact erase/redraw
  transition, so the Hubble redraw threshold now matches the documented
  occlusion rule.  Retain System GL Hubble qualification as open.

## P0: release qualification

- Run the exact-current graphical matrix from truly cold and certified warm
  caches.  The exact-current 50k OSMesa shaded/wire rows are green; the 150k
  150k OSMesa shaded/wire cold/warm rows are green under the bounded
  diagnostic environment.  Re-run them on production hardware scale rows
  before release.
  Re-run representative real-model rows if their cache or policy path changes.
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
- **Validated 2026-08-23:** the dedicated fractional-DPR qged selection replay
  passed on both OSMesa and System GL using the real view-selection palette.
  It covers point select/remove, held drag-rectangle placement in physical
  canvas coordinates, selection-list readback, and selected erase/redraw.
  The runner now accepts either the installed or regression `moss.g` fixture,
  so a space-conscious partial build does not silently lose this coverage.
  This is one matrix row, not clearance for hierarchy-scale selection or the
  remaining controls below.
- Complete the qged selection, polygon, measurement, and faceplate matrix:
  point/rectangle selection and modifiers; tree styling and draw/erase under
  selection; polygon creation/manipulation/booleans; all measurement modes;
  lighting/navigation/axes/grid/framebuffer controls; resize and fractional
  DPR in single and quad layouts.
- Complete the primitive/sketch work described in `qged_editing.md`: the
  librt descriptor contract now classifies all registered operations, but
  rejection/readback behavior, the reusable manipulator vocabulary,
  specialized sketch curves, and common qged/MGED/gsh sessions including
  sed/oed still require qualification.

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
