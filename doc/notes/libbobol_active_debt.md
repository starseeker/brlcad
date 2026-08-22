# libBObol active debt

Last reviewed: 2026-08-22

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

- Run the complete current-build `bobol_headless`, `drawing_baseline`, public
  API/symbol, client, and shared-library lifecycle gates after the last source
  change.  The 2026-08-22 cleanup build passes the combined 86-test headless/
  drawing selection (one unavailable X11 capability skip); repeat it after any
  subsequent production change.
- Run the organized System GL and OSMesa matrices from truly cold and certified
  warm caches.  A partial cache may never be labeled warm.  Generic Twin
  shaded cold/warm passes on both backends after the current extraction and
  rotation-safe autoview work; the remaining model/mode rows are still open.
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

- Improve the 150k distinct-asset post-rotation refinement path, which has
  terminated correctly but has taken roughly 218 seconds on this host.
  Profile before changing policy; preserve bounded memory and terminal proof.
- Complete the large OSMesa matrix.  Continue reducing repeated owner-thread
  demand/retarget work, but distinguish it from legitimate software triangle
  setup, lighting, and rasterization.
- Measure the new database-discovery path under cold local storage, warm page
  cache, and slower storage.  Worker count and inflight bytes should respond to
  I/O throughput rather than multiplying contention.
- Validate repeated-instance visibility turnover at 50k/150k-equivalent scale
  separately from unique-asset stress.
- Validate adaptive BREP tessellation growth and zoom-out memory recovery over
  the full NIST corpus and a larger real BREP model.
- Establish release thresholds for time to first useful representation,
  interactive latency, stable convergence, and peak/resident memory.

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
