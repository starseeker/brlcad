# libBObol active debt

Last reviewed: 2026-08-23

This is the only active TODO inventory for the BRL-CAD Obol drawing stack.
Detailed qualification steps are in `obol_production_readiness.md`; primitive
and sketch editing details are in `qged_editing.md`.  Completed migration
history and resolved defects belong in `libbobol_engineering_lessons.md`.

## P0: stabilize the current substrate

- The first behavior-preserving extraction pass is complete: database cache,
  compact access, interaction, summary, and type responsibilities; controller
  lighting, progressive pump, render execution, and scene forwarding;
  libged endpoint, geometry, overlay, and scene-record bridges; and Obol's
  software wire renderer now compile independently.  Review the remaining
  concentrations (`database_source.cpp`, progressive transaction work in
  `draw_obol.cpp`, and `SoCADAssembly` frame planning/execution) only at a real
  ownership or profiling seam.  Preserve compact indexed storage, sparse
  updates, bounded owner-thread windows, and immutable payload ownership.
- Replace the monolithic static `ged_test_obol_draw_sync` white-box executable
  with focused owning-library tests and a normal shared-stack GED integration
  test.  Static libraries should otherwise be exercised only by explicit
  static-link tests.
- Complete the current parallel database-discovery integration and prove
  parity with ordinary librt hierarchy/reference publication, failure cleanup,
  read-only concurrency, memory bounds, mapped-file fallback, and clean
  adoption by GED's normal database index.
- Retain the source-realization/cache static-lifetime barrier and its regression
  which exits with a cache callback active.  Run the full shared BRL-CAD,
  libBObol, and Obol ASan/UBSan stack after the source extraction.

## P0: final production qualification

- The current-source `bobol_headless` and `drawing_baseline` selection passes
  all 86 tests in 389 seconds, including public API/symbol, client, and
  shared-library lifecycle coverage; one unavailable X11 capability is skipped
  as expected.  Repeat after any subsequent production change.
- Pre-admission-change OSMesa evidence covers Generic Twin, Lucy, Havoc, and
  Hubble in shaded/wire cold/certified-warm rows, plus 50k and 150k
  distinct-asset shaded/wire cold/warm rows.  The new joint mesh/point
  allocator has focused coverage: its 50k OSMesa matrix remains green; a fresh
  Lucy cold start now publishes one resident compact mesh at the first stable
  checkpoint (26.3 seconds, 144,976 faces, no point proxy); and the 150k
  terminal allocation removes more than half of the earlier prominent-floor
  violations without boxes.  Re-run the complete OSMesa 16-row and scale
  matrix after this policy change before calling it current qualification.  The
  current host cannot rerun System GL because `/tmp/.X11-unix` is owned by
  `nobody:nogroup`, which prevents Xvfb from establishing a local socket.  Do
  not treat earlier System GL captures as current-binary evidence; rerun that
  matrix on a valid X host.
- Run the organized System GL and OSMesa matrices from truly cold and certified
  warm caches.  A partial cache may never be labeled warm.  Generic Twin
  shaded/wire cold/warm passes on both backends after the current extraction
  and rotation-safe autoview work.  Hubble shaded/wire cold/warm passes on both
  backends with a representative `PANEL_C01` hierarchy target.  The current 50k
  and 150k distinct-asset fixtures pass System GL shaded cold/warm with exact
  terminal coverage and no boxes; the 50k warm perf run has no dominant
  controller/render livelock.  Lucy shaded and wire cold/warm passes on both
  backends with spatial refinement, exact-view recall, zoom-out compaction,
  zero boxes, and explicit software-performance limits.  Havoc shaded/wire
  cold/warm also passes on both backends.  The 50k OSMesa shaded cold/warm rows
  pass in 35/50 seconds, and a fresh-cache 50k OSMesa wire pair passes in
  50/75 seconds with zero terminal boxes, explicit software-performance limits,
  and no stagnant post-rotation dispatches.  The 150k OSMesa shaded cold/warm
  rows now pass in 98/72 seconds after fixing structural-fallback recovery:
  both finish at an aggregate point threshold of 64, with zero boxes and zero
  stagnant dispatches.  The cold/warm terminal populations retain 2,119/2,855
  mesh occurrences and classify roughly 148k subpixel occurrences.  The 150k
  OSMesa wire warm interaction/lifecycle row now also passes in 126 seconds
  with exact 150,001-occurrence classification, zero boxes, no pending work,
  and no terminal renderer-wide PoP ceiling.  It specifically covers zoom,
  pose motion, tree selection, exact subpath erase, and redraw.  Cold wire
  validation and terminal visual-quality consistency remain open: the focused
  closed-population run spent 98 ms and limited worst projected error to about
  18 pixels, while the complete row settled near 53 ms and about 60 pixels.
  Do not promote the row until that unused static allowance is allocated
  consistently by visual importance.
- The 2026-08-23 current-source regression fixes two late continuity defects:
  a completed deadline-safe interactive frame now survives the policy-only
  quiet handoff for the same camera/source-quality domain, and fractional
  pixel targets remain fractional in the submit-action scene ceiling.  Lucy
  shaded passes on both System GL and OSMesa from cold and certified warm
  caches (51/26 and 55/32 seconds); its return view is cut 21 at 0.645 pixels
  for the 0.75-pixel request.  Hubble OSMesa shaded/wire cold/warm now passes
  through an explicit canvas-presentable barrier.  The warm shared-cache 50k
  and 150k OSMesa shaded rows also pass in 24 and 53 seconds, respectively,
  with zero boxes.  The latter ends at 2,066 mesh occurrences plus 147,934
  subpixel points and a 60 ms retained terminal frame.  Repeat cold 50k/150k
  rows and the full shared-stack gate after the final source change.
- The current static-quality retry and logical-occurrence fixes pass the full
  OSMesa graphical matrix again: Generic Twin, Lucy, Havoc, and Hubble in
  shaded/wire cold/certified-warm (16 rows), followed by 50k and 150k
  distinct-asset shaded/wire cold/warm rows.  A first static frame may retry
  once when newly published presentation resources make that frame exceed its
  deadline; only a second unchanged miss becomes a capacity witness.  A source
  with one logical occurrence may never use point aggregation merely because
  its PoP implementation has many internal spatial clusters.  The broad CTest
  sweep also passes every libBObol, drawing, qged, qtcad/Obol, gsh, and
  MGED/LoD row.  Re-run these exact rows after any policy change.
- Qualify Generic Twin, Lucy, multi-Lucy/xpush, Hubble, havoc, NIST BREP,
  Stanford meshes, and the varied 50k/150k fixtures.  Include shaded,
  wireframe, evaluated modes where applicable, LoD auto/off transitions,
  resize, zoom, rotation, translation, selection, and exact/subpath erase and
  redraw.
- Record time to overview, leaf coverage, first useful mesh, ready view, and
  stable return after interaction; p50/p95/max input latency; completed-frame
  time/FPS; resident/peak memory; cache growth; and terminal representation
  counts.
- Compare System GL and OSMesa semantics and images within explicit tolerances.
  Use APNG or apitrace for any remaining flash, corruption, camera jump, or GL
  state disagreement.
- Exercise independent real multi-gigabyte vehicle models.  Synthetic fixtures
  do not prove visual-importance behavior for wheels, blades, hulls, fasteners,
  and mixed-size unique assemblies.

## P0: qged GUI interaction audit

Programmatically drive the actual controls and mouse paths, not only equivalent
GED commands.  Every mode must prove widget state, semantic state, retained
scene state, and command readback agree after either GUI- or command-originated
changes.

- View selection: point, rectangle, append/toggle/subtract, hierarchy impact,
  whole-row tree styling, draw/erase under selection, single-edit promotion,
  multi-selection, resize, HiDPI, and large Hubble selection latency.
- Polygon Create/Modify: mouse creation, constrained resize, movement, point
  editing, append/delete, boolean union/intersection/difference, holes/islands,
  open contours, snapping, style changes, persistence, cancellation, tool
  handoff, single/quad views, and CAD composition.
- Primitive editing: generated and custom controls, CLI/GUI/manipulator
  round-trips, commit/cancel/revert, invalid-input immutability, plugin reload,
  multiple views, selection promotion, and database mutation invalidation.
- Sketch editing: direct viewport curve construction, vertex/segment/topology
  changes, NURB and arc affordances, non-unit planes, linked extrude/revolve
  profiles, persistence, and cancellation.
- Measurement: every measurement mode, mouse creation and adjustment, snap
  modes, unit changes, labels, selection/highlight, resize/HiDPI, command
  readback, cleanup, and operation with progressive CAD still refining.
- Remaining controls: camera/orbit, lighting profiles, nav gizmo styles,
  faceplate/axes/grid/ADC, framebuffer modes, raytrace composition, snapping,
  view settings, export, hierarchy/tree filters, and plugin unload/reload.

## P1: editing completion

- Classify every installed librt primitive edit operation as generated, action,
  custom, or unsupported.  Remove remaining implicit `INFER` behavior before a
  command may drive a retained manipulator.
- Audit every primitive handler for reject-without-mutation and set/readback
  round trips.  Do not hide faulty handlers with whole-primitive rollback in
  the generic hot path.
- Normalize `ft_set_edit_mode` side effects for EBM, DSP, extrude, CLINE, PIPE,
  NMG, and any other transitional handler.
- Add constrained point, plane, axis, and rotation-ring manipulators without
  publishing per-element scene nodes.
- Finish specialized sketch NURB/arc interaction and direct viewport curve
  creation.  Replace extrude/revolve's transitional profile helper with the
  shared sketch adapter.
- Collapse migrated primitive plugins into one registered edit tool once the
  custom adapters no longer need separate lifecycle behavior.
- Complete MGED sed/oed image and interaction qualification on the common GED
  edit-session contract, including alternate occurrence conflicts, erase under
  edit, and exact retained-scene restoration.

## P1: scale and performance

- Preserve the current 150k System GL scale result: cold/warm scripted rows
  complete in 42/39 seconds on this host rather than the former roughly
  218-second post-rotation path, with exact 150,001-leaf coverage, zero boxes,
  and bounded transient mesh work.  Re-profile before changing scale policy.
- Preserve the 50k structural-admission regression: an exact point/box
  classifier frame is a presentation acknowledgement even before any managed
  mesh makes the frame capacity-relevant.  It must release the first mesh wave
  rather than repaint structural coverage indefinitely.
- Preserve the 50k OSMesa wire deadline regression.  A hard-aborted minimum-
  prefix frame must advance small-part aggregation even inside the ordinary
  completed-frame timing deadband, and a rejected one-pixel static trial must
  remain rejected for the current view/capacity epoch.  Point-threshold
  recovery cannot reopen the identical one-pixel population.
- Preserve the 150k OSMesa structural-fallback recovery regression.  Aggregate
  structural points are not hidden richer meshes: lowering their point
  threshold exposes boxes.  Triangle-prefix recovery must retain the proven
  aggregate cut while structural fallbacks remain, and a static one-pixel
  trial is admissible only after the visible population is fully realized.
- Preserve the 150k resident-growth transaction boundary.  A minimax
  occurrence allocation may start only after provider inventory, mesh tasks,
  queued/result publication, and the coalesced residency drain are all closed.
  Treating provider inventory alone as settled made every bounded cache wave
  invalidate and rerun the scene allocation, producing the observed
  coarse/fine cycling.  A hard headroom miss is a negative witness, not a
  cancellation which may immediately re-arm the same population.
- Preserve the allocated-residency boundary.  Quiet views load only through
  the allocator-selected presentation cut; finer physical demand stays
  recorded as quality debt until a view or capacity edge changes that
  allocation.  Active scale interaction alone may prefetch one bounded
  transition beyond it.  Unconstrained quiet prefetch made each 150k cache
  result wave invalidate the closed population and rerun scene allocation,
  while suppressing prefetch entirely stranded Lucy during zoom.  The focused
  TLA+ model and submit-action test cover both sides of this rule.  Repeat the
  complete cold/warm 50k/150k shaded/wire lifecycle on System GL and OSMesa
  after the change.  Focused Lucy and 150k shaded-prefix rows pass, as do the
  full certified-warm OSMesa 50k shaded and 150k shaded/wire lifecycles; 50k
  wire, truly cold rows, and current System GL evidence remain outstanding.
- Preserve the interrupted-presentation transaction boundary.  Once a CAD
  traversal has retained resumable preparation, owner-thread provider, result,
  compaction, and cut publication must remain frozen until the exact successor
  frame completes.  Workers may continue only behind their bounded queues.
  The focused controller regression proves both the freeze and its
  completed-frame release.
- Preserve the static first-presentation retry boundary.  An initial static
  trial can pay immutable-buffer installation, command setup, or software
  warmup not reflected by the preparation serial.  It gets exactly one
  unchanged retry; a second miss is the negative capacity witness that bounds
  the retained allocation.  Do not make every miss retryable, or restore the
  old coarse/fine loop.
- Preserve source-logical occurrence semantics.  A spatially paged single mesh
  is one visual object, not one object per PoP cluster, and cannot be reduced
  to aggregate point coverage.  For multi-occurrence sources, keep using the
  current visible-occurrence census so an otherwise large assembly with only
  one visible part is also not point-aggregated spuriously.
- The large OSMesa matrix now passes at 50k and 150k in both shaded and wire,
  cold and certified warm.  Continue reducing repeated owner-thread
  demand/retarget work, but distinguish it from legitimate software triangle
  setup, lighting, and rasterization.  Re-run the rows after any policy change.
- Reduce the remaining 150k OSMesa convergence variance.  Current warm shaded
  runs reach exact 150,001-occurrence classification with zero boxes and a
  roughly 82--90 ms terminal frame, but setup-to-idle has ranged from about 16
  to 45 seconds depending on the number of resumable preparation/deadline
  slices.  The terminal presentation is viable; repeated closed-population
  allocation and prepared-command transactions are the remaining target.
  A current 50k OSMesa profile assigns about 32% of sampled CPU to software
  rasterization, 31% to libBObol (including roughly 22% inclusive in bounded
  compact submission), and 11% to Obol.  Reduce repeated compact planning and
  static allocation transactions without relaxing the bounded owner-thread
  slices or treating the legitimate software raster cost as a logic defect.
- Completed 2026-08-23: Obol retains the exact one-occurrence logical proxy
  stream for coverage, selection, highlight, picking, and mesh-promotion
  decisions, but the software renderer submits a cached camera-local screen-bin
  stream.  The representative is deterministic (selected/hovered/color
  emphasis, then nearest depth, then stable instance ID); ordinary GL retains
  the full-density stream.  The software path uses one native screen pixel
  where possible and expands bins only to honor its 32,768-point cap.  A 150k
  OSMesa shaded warm interaction/selection/erase/redraw qualification passed
  in 49 seconds with 148,619 logical proxies, 3,457 submitted proxy points,
  zero boxes, convergence fraction 1, and a 97.0 ms terminal frame.  The 50k
  OSMesa shaded cold/warm rows also pass in 51/40 seconds with roughly
  46.7k/46.9k logical proxies, roughly 3.29k submitted points, zero boxes, and
  convergence fraction 1.  The qged scale matrix now rejects
  a high-cardinality OSMesa row unless its physical proxy count is nonzero,
  lower than its logical coverage, and within the cap.  Repeat the broad
  real-model, wire, dense-overlap, selected-color, HiDPI, and resize rows after
  this renderer change before calling it release-qualified.
- A current 150k OSMesa cold shaded profile records no single liveness hot
  loop: libBObol accounts for 34.4% of samples, libosmesa for 11.8%, and Obol
  for 8.3%.  The largest named source hot spot is parallel compact coverage
  bounds (3.2%), alongside expected software triangle rasterization (3.4%).
  Split discovery/coverage from view-presentation measurements before making a
  policy change; this whole-run profile includes both phases.
- Measure the new database-discovery path under cold local storage, warm page
  cache, and slower storage.  Worker count and inflight bytes should respond to
  I/O throughput rather than multiplying contention.
- Validate repeated-instance visibility turnover at 50k/150k-equivalent scale
  separately from unique-asset stress.
- Validate adaptive BREP tessellation growth and zoom-out memory recovery over
  the full NIST corpus and a larger real BREP model.
- Establish release thresholds for time to first useful representation,
  interactive latency, stable convergence, and peak/resident memory.
- Redesign and qualify the significance frontier under a performance-limited
  150k software frame.  The current terminal images are coherent and have no
  structural boxes, but a certified-warm shaded zoom/rotation sequence still
  leaves about 1.2k of 2.45k prominent occurrences below the protected visual
  floor at an 80--97 ms frame.  The immediate defect in that path was fixed:
  a deadline miss supplied a 95-percent strict ceiling and marginal recovery
  applied a second 80-percent reduction.  Do not reintroduce that compounded
  margin.  The remaining problem is not safely solved by promoting every
  point-to-minimum mesh transition or by treating every point-to-protected
  transition as atomic; both focused experiments reduced useful coverage on
  this workload.  Define an explicit, measurable frontier using projected
  error reduction, screen footprint/silhouette, visual emphasis, complete
  transition cost, and frame budget.  Validate it with per-occurrence
  certificates and real wheels/blades/hulls as well as mixed-size 50k/150k
  fixtures before changing allocator priority again.
- Harden the `libBObol_lod_service` test deadline.  A focused run completed in
  9.92 seconds against its exact 10-second timeout after one timeout failure;
  qualification should not depend on scheduler luck.

## P1: physical maintainability after extraction

- Keep each new private translation unit aligned with one owner and lifecycle.
  Do not split merely by line count or create bidirectional private APIs.
- Continue reducing implementation concentration where perf or ownership
  review identifies a real seam.  File size alone is not permission to replace
  compact storage with per-leaf objects or add another publication route.
- Keep internal symbols hidden and private headers uninstalled.  Installed C
  headers must compile independently as C11 and C++17 and match the reviewed
  symbol manifests.

## Platform and optional capability debt

- Windows: native WGL/Qt/Tk lifecycle, fixed-width cache, BREP/PoP zoom,
  sanitizer, memory diagnostics, and the formerly observed transient invalid
  vertex regression.
- macOS: select and qualify a native Coin/Obol context provider, then add host,
  lifecycle, and sanitizer coverage.
- TSan/LSan: run on a compatible native Linux worker when restricted containers
  cannot initialize them.
- Stereo: do not advertise it until a concrete host passes capability,
  mode-output, capture, resize, and mono-restoration tests.
- Appleseed: validate the imgstream tile publication path in an enabled build.
- External custom hosts: require factory capability descriptors before removing
  the current direct-host compatibility seam.

## Definition of done

The drawing stack is ready for production when the exact final binaries pass
the complete lower-level and graphical matrix, no known user-facing correctness
or responsiveness blocker remains, large-scene results meet or improve the
recorded baseline, every constrained terminal state reports its reason, and
the code/documentation tree contains one supported path with no unowned
compatibility or synchronization mechanism.
