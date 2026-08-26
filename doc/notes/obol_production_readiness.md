# Obol production-readiness matrix

Last reviewed: 2026-08-25

This is the release qualification plan for the BRL-CAD Obol drawing stack.  A
row is green only when it was run against the final binaries and records its
model, renderer, cache state, event script, viewport/DPR, report, images, and
performance evidence.  Historical measurements belong in
`libbobol_engineering_lessons.md`, not here.

## Current evidence snapshot

| Area | Current evidence | Release status |
|---|---|---|
| focused CTest gate | The dependency-complete 2026-08-25 build passes all 28 `bobol_headless` tests in 4.83 s, including renderer, LoD service/coordinator/update, cache, compact ownership, retained allocation, edit manipulator, host, API/symbol, and GED draw-sync contracts.  Independent draw scope is exercised with a real local endpoint across draw/erase/redraw/zap.  The focused interaction gate also passes semantic selection, GED edit, qtcad primitive edit/preview/selection/measurement, qged event replay, polygon/sketch replay, and qged primitive-edit replay in single and quad layouts (13/13). | broader graphical production suite still required |
| control-plane models | Focused publication, host-work, convergence, admission, arbitration, canonical-pipeline, cold-preview, static-quality, renderer-preparation, interaction-session, deadline-ownership, capacity-search, resident-growth, and point-recovery TLC models pass.  `ObolProgressivePipeline` explored 2,986,118 generated / 1,427,962 distinct states to depth 42.  The focused models cover strict finite-work progress, exclusive successor ownership, bounded safe/unsafe search, a typed resident-growth transaction, and deferred point-recovery completion.  Production now uses a typed interaction session, canonical allocation key, exclusive deadline successor, and progress-only HUD publication.  qged reports the five-domain revisions plus plan/transaction/frame identities, and its external checker rejects sampled regression, duplicate-plan, spontaneous-reopen, and unwitnessed A/B/A traces; the focused checker test and exact Generic Twin OSMesa cold/warm run pass. | formal and sampled-runtime boundaries green; remaining ordinary capacity branches and the complete event/effect reducer with per-transition records remain |
| Lucy OSMesa | The exact-current warm shaded lifecycle reaches a terminal box-free 6.0M-face cut with no quality-floor violation.  It still fails the first-useful-image oracle: at 0.2/1.5 s the available 60.6k-face cut is a rectangular slab rather than a recognizable Lucy silhouette.  This is a candidate-quality defect, not permission to weaken the visual test or add a Lucy control path. | warm first-candidate visual quality open; pose/zoom quality and cold/wire rows still require final qualification |
| 50k scale | Cold shaded scale checks pass on System GL and OSMesa in 8.7/10.0 s to terminal CONSTRAINED.  The exact-current full warm OSMesa interaction script passes in 42 s after stale deadline replay ownership was retired: all checkpoints are terminal, box-free, quality-debt-free, and exercise motion, zoom, selection, subpath erase, and redraw. | shaded cold/warm OSMesa liveness green; wire and realistic visual-significance qualification remain |
| 150k scale | The bounded cold crash gate passed under a 16 GiB address-space cap on System GL and OSMesa with all 150,001 leaves discovered and no terminal boxes or failures.  The exact-current warm shaded OSMesa lifecycle at `/tmp/qged-control-terminal-150k-final-20260825` now completes its full interaction/selection script at 100 percent with phase idle, no control owner or obligation, no structural boxes, and a passing corrected report validator.  Exclusive capacity/generic successor selection, the typed resident-growth drain/reallocation transaction, and level-triggered point-recovery completion remove the prior repeated allocation and ownerless 99-percent states. | shaded OSMesa liveness green; cold/wire, production-hardware timing, and realistic visual-significance qualification remain |
| System GL smoke | Generic Twin shaded cold/warm and wire cold/warm passed after the camera-depth change.  The exact-current 2026-08-24 System GL shaded cold/warm replay is retained at `/tmp/qged-generic-system-20260824`; the wire replay is `/tmp/qged-generic-wire-system-debug-20260824`.  Both terminal wire frames contain 709 mesh payloads, zero boxes/failures/pending work, and matching cold/warm camera state.  The exact-current Lucy shaded cold/warm interaction replay at `/tmp/qged-lucy-system-20260824` also passed: both final frames are exact/ready, box-free, and quality-floor-clean with one resident source payload, 7.41M displayed PoP faces, a 230 MB resident mesh set, and one quality-history recall.  Hubble shaded/wire cold/warm now pass after the overview lifecycle repair: the visible whole-target extent remains through early discovery and retires once useful leaf coverage is committed, rather than persisting as a terminal fallback.  A pre-created, empty `BU_DIR_CACHE` Lucy replay now proves separate milestones: first useful payload at 8.25 s (559k faces, zero boxes), terminal idle at 57.7 s; the same cache warm-starts at 0.28 s and idles at 1.76 s.  The GUI driver records payload-ready separately, and the System-GL terminal render-only wake is no longer coalesced indefinitely. | Qualify Lucy wire and complete the System-GL large-model matrix |
| OSMesa Generic Twin | Exact-current 2026-08-24 shaded cold/warm replay at `/tmp/qged-generic-osmesa-20260824` passed.  Both final views are exact/ready with 709 payloads, zero structural boxes/failures/pending work, zero quality-floor debt, and the System-GL-matching camera scale; the software path recorded four/five bounded interrupted frames during cold/warm convergence.  The exact-current wire cold run at `/tmp/qged-generic-wire-osmesa-20260824` and validated-cache warm run at `/tmp/qged-generic-wire-osmesa-warm-20260824` both passed: each has 709 payloads, zero boxes/failures/pending work, and byte-identical final autoview state. | focused regression green; qualify Lucy wire and the larger model matrix |
| spatial Lucy, System GL | Exact-current 2026-08-25 warm shaded retained-page replay passes.  It refines during continuous zoom, reaches its quiet pixel target, compacts on zoom-out, restores the prior view through quality history, and has no boxes or page-level subpixel proxies. | focused retained-page regression green; cold and wire rows remain |
| direct-mesh Generic Twin | Exact-current 2026-08-25 shaded cold/warm replays pass on System GL and OSMesa.  Each terminal frame has 708 direct full-detail payloads, one progressive payload, all 709 BoT occurrences, zero structural boxes, and no pending work.  The overall extent preview remains visible during discovery; it is not counted as a semantic leaf. | focused admission regression green; keep discovery-preview latency distinct from terminal-mesh admission |
| Hubble OSMesa | Exact-current 2026-08-25 shaded cold/warm lifecycles and camera contracts pass after the static-quality restoration change; the prior shaded/wire evidence remains green, with 1,804 warm shaded payloads and zero structural boxes. | focused regression green; rerun wire with the final binary and continue the full model matrix |

The 150k terminal view may validly be performance-limited.  It is not visually
qualified until prominent-object quality floors are demonstrated on realistic
large models, not merely synthetic coverage fixtures.

## Common acceptance criteria

Every graphical row must show:

- correct final camera, extent, draw mode, colors/materials, lighting, and
  retained scene composition;
- exact visible occurrence classification, no unexpected boxes, no empty
  frame, no invalid geometry, and coherent shaded/wire cut semantics;
- no owner-thread pending-without-witness loop, stale result, coordinator
  invariant violation, or unbounded worker/resident allocation;
- semantic selection, tree state, highlights, erase/redraw, and edit
  promotion consistent with command and widget state;
- stable-frame, first-useful-image, interaction-latency, resident/peak-memory,
  cache, and representation telemetry within the row's declared threshold;
- visually inspected checkpoint images.  APNG is required for suspected
  flicker; apitrace is required for System-GL-only corruption or state leaks.

The GUI runner may allow cold asset construction time, but a timeout, crash,
sanitizer finding, empty capture, unexplained terminal box, pending idle loop,
or a missing report fails the row.

## Required models

| Model | Primary purpose |
|---|---|
| Generic Twin | production visual baseline; tail/hull/engine, autoview, lighting |
| Lucy | spatial-page demand, close zoom, smooth refinement, compaction |
| multi-Lucy/xpush | shared-asset reuse and visibility turnover |
| Hubble | deep hierarchy, selection/tree scale, small-component importance |
| Havoc | mixed hierarchy and real-model interaction |
| NIST BREP corpus | wire/shaded BREP and adaptive tessellation behavior |
| Stanford meshes | varied manifold/non-manifold mesh behavior |
| 5k/50k/150k varied fixtures | distinct assets, mixed sizes/colors/regions/hierarchy, aggregate budgets |
| independent multi-gigabyte vehicle | real visual significance, I/O, memory, and unique-mesh behavior |

Synthetic fixtures must include mostly manifold inputs, mixed mesh sizes,
visually prominent components, color/region variation, and deep hierarchy.
Repeated-instance and distinct-asset cases are separate obligations.

## Matrix

### Lower-level/shared stack

- Clean Ninja build of qged, MGED, gsh, Archer/TkObol, rtwizard, and plugins.
- `bobol_contract`, `bobol_headless`, `drawing_baseline`, libbg polygon,
  librt discovery/edit, libbu cache, installed-package consumer, public-header,
  and symbol-manifest tests.
- TLC configurations for occurrence publication, host work, bounded
  Lucy/large-scene convergence, and deferred-autoview ownership.  The latter
  proves control ownership only; real qged cold/warm camera replays remain
  required for Qt/renderer timing and image behavior.
- Shared-stack ASan/UBSan: workers active during teardown, cache corruption,
  endpoint replacement, plugin reload, edit cancellation, and rapid close.
- Native-worker TSan/LSan; do not count a container runtime limitation as a
  successful dynamic analysis run.

### qged controls and interaction

Run each applicable row in single and quad layouts, DPR 1 and fractional DPR,
before/after resize, on System GL and OSMesa.

- Camera: MGED orbit, rotate, translate, center, zoom, smooth round trip,
  autoview, aspect changes, navigation gizmo, and camera/history readback.
  Close orthographic zoom must retain its scene-depth range and never behave
  as an implicit cut; verify this at extreme Hubble zoom on both renderers.
- Selection: point and rectangle, append/toggle/subtract, tree-to-scene and
  scene-to-tree propagation, selection under erase/redraw, and Hubble latency.
- Faceplate: axes/ticks, grid, ADC, HUD progress, lighting profiles,
  framebuffer modes, raytrace overlay/underlay/interlay, snapping, and view
  settings.  The `view cutting` controller is disabled by default and is
  independent of camera clipping; qualify its visible plane affordance,
  faceplate controls, persistence, exact picking, and both renderers.
- Polygon/measurement: direct mouse creation and editing, booleans, snapping,
  persistence, cancellation, command readback, all measurement modes, labels,
  units, and resize.
- Primitive/sketch editing: actual widgets and mouse gestures, CLI/widget/
  manipulator round trips, commit/cancel/revert, invalid-input immutability,
  plugin reload, multiple views, and database mutation invalidation.

### Drawing and LoD

- Modes 0--5 on appropriate moss/rook models, with LoD auto and off.
- Generic Twin shaded/wire cold/warm and LoD auto/off; compare `ae 90 0` and
  oblique views to the main baseline for skin seams, engine, tail, underside,
  lighting, and boxes.
- Lucy cold/warm zoom in/out and rotation; check spatial page coverage,
  crack-free wings, retained history, in-motion refinement, uniform zoom-out
  quality, and resident-memory recovery.
- Hubble shaded/wire selection, resize, erase/redraw, and importance floors.
- NIST BREP shaded/wire, LoD auto/off, adaptive tessellation growth, and
  zoom-out reclamation.
- Multi-Lucy/xpush turnover; bring old and new instances into view while
  retaining shared assets and interaction responsiveness.
- 5k/50k/150k cold and certified-warm shaded/wire runs, including rotation,
  zoom, selection, exact subpath erase/redraw, and retained terminal frames.

### Discovery, performance, and visual significance

- Measure read-only discovery on cold local storage, warm page cache, and
  slower storage.  Record elapsed time, worker concurrency, I/O behavior,
  peak bytes, and tree/GED parity.
- For 50k and a real model, collect perf.  Separate discovery, compact
  planning, worker/cache, upload/preparation, and raster work from diagnostic
  JSON/image compression.
- Capture APNG for any visual cycle; compare terminal and representative active
  frames to the main baseline where available.
- Record time to overview, leaf coverage, first useful mesh, view ready, stable
  return, p50/p95/max input latency, completed-frame time/FPS, resident/peak
  memory, cache growth, and quality-floor debt.
- Evaluate visual-importance results using wheels, blades, tails, hulls, and
  fasteners.  Synthetic terminal coverage alone is not evidence of adequate
  visual quality.

## Evidence retention

Keep one shared copy of reusable large inputs and canonical warm caches.  Keep
the event script, final report, representative checkpoint images, and perf/APNG
only when they explain a result.  Remove redundant transient screenshots,
trace logs, and duplicate generated models after recording the needed evidence.
