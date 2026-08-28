# Obol production-readiness matrix

Last reviewed: 2026-08-26

This is the release qualification plan for the BRL-CAD Obol drawing stack.  A
row is green only when it was run against the final binaries and records its
model, renderer, cache state, event script, viewport/DPR, report, images, and
performance evidence.  Historical measurements belong in
`libbobol_engineering_lessons.md`, not here.

## Current evidence snapshot

| Area | Current evidence | Release status |
|---|---|---|
| focused CTest gate | The dependency-complete 2026-08-26 build passes all 28 `bobol_headless` tests in 6.13 s real / 20.31 process-seconds, including renderer, LoD service/coordinator/update, cache, compact ownership, retained allocation, edit manipulator, host, API/symbol, and GED draw-sync contracts.  The service test proves that a coalesced asset producer may complete against its retained latest demand; the view-state test proves that a superseded result cannot create a terminal occurrence failure, while genuinely stale source/cache data retains its failure semantics.  Independent draw scope is exercised with a real local endpoint across draw/erase/redraw/zap.  The linked Obol CAD suite passes 6/6, including its 131,072-occurrence classifier-reservation regression.  The focused interaction gate also passes semantic selection, GED edit, qtcad primitive edit/preview/selection/measurement, qged event replay, polygon/sketch replay, and qged primitive-edit replay in single and quad layouts (13/13). | broader graphical production suite still required |
| control-plane models | Focused publication, host-work, convergence, admission, arbitration, canonical-pipeline, cold-preview, static-quality, renderer-preparation, interaction-session, deadline-ownership, capacity-search/presentation-handoff, resident-growth, point-recovery, structural-frontier, and active-producer-demand TLC models pass.  `ObolProgressivePipeline` explored 2,986,118 generated / 1,427,962 distinct states to depth 42; the focused capacity-presentation model covers changed/already-applied plans, effective/inert ceilings, over-budget applied allocations, and availability restriction before cut application in 178 generated / 138 distinct states.  The point-quality model distinguishes adaptive calibration, handoff confirmation, and triangle recovery; covers static rejection plus a point request arriving during an active finite capacity search; and passes in 29 generated / 17 distinct states to depth 5.  The structural-frontier model checks 30 distinct states to depth 14 and proves that an exact nonempty fallback frontier always owns an enabled point-frame or mesh-repair successor and cannot publish readiness.  The active-producer model checks 340 distinct states to depth 12, including liveness, and proves that immutable asset work can survive changing demand without converting an overtaken result into a current failure.  Production now uses a typed interaction session, a finite retained-allocation request, canonical allocation identity including resident admission, an exclusive deadline successor, and progress-only HUD publication.  qged reports the five-domain revisions plus plan/transaction/frame identities and retained HUD text; its external checker rejects sampled regression, duplicate-plan, spontaneous reopen, unwitnessed A/B/A traces, and `View ready` unless the snapshot has LoD state, is terminal and view-ready, has no background work, and is exactly 100 percent. | formal and sampled-runtime boundaries green; remaining ordinary capacity branches and the complete event/effect reducer with per-transition records remain |
| Lucy OSMesa | The exact-current 2026-08-26 true-cold shaded lifecycle at `/tmp/qged-lucy-latest-demand-contract-20260826` passes in 76 seconds and separately records globally representative coverage and the first real CAD mesh.  Camera changes made while the immutable hierarchy is being built no longer discard its live pages or final result: the producer resolves page selection against the service-owned latest demand, and service validation uses that same demand.  Two consecutive certified-warm replays after splitting superseded work from stale source failures pass in 63/62 seconds at `/tmp/qged-lucy-superseded-fix-{a,b}-20260826`, each with a terminal ready, box-free mesh and zero occurrence failures.  The typed submission-owner replay at `/tmp/qged-submission-owner-lucy-20260826` passes in 82 seconds with 1.95M presented faces and 0.632-pixel maximum certified error.  After the exact source-delta owner replaced its independent active/source/plan fields, `/tmp/qged-submission-delta-lucy-20260826` also passes in 94 seconds with a box-free 1.265-pixel terminal view.  These lifecycles exercise continuous zoom/refinement, zoom-out recovery, rotation, lighting, hierarchy selection, subpath erase/redraw, camera ownership, and the strict HUD contract.  The exact-current warm wire lifecycle at `/tmp/qged-production-wire-20260826/lucy-osmesa-warm-limit-fix` also passes in 68 s. | cold/warm shaded and warm wire OSMesa interaction rows green; wire cold remains |
| 50k scale | Cold shaded scale checks pass on System GL and OSMesa in 8.7/10.0 s to terminal CONSTRAINED.  A later exact-current run exposed interactive deadline recovery walking down one PoP ordinal per missed frame even when the measured cost ratio already proved that insufficient.  Recovery now operates in render-cost space during motion as well as at rest; the prior-pose deadline floor remains quiet-only.  The warm OSMesa shaded matrix at `/tmp/qged-50k-superseded-fix-20260826` passes in 78 seconds: held-motion recovery reaches its responsive ceiling directly, and the terminal view has 7,566 mesh payloads, 1.52M faces, zero boxes/failures, and no control owner.  After replacing the raw source-plan/census latches with typed owner-thread values, `/tmp/qged-submission-owner-50k-20260826` passes in 72 seconds and terminates ready with 6,633 progressive payloads, 1.32M faces, zero boxes/failures, and no pending work.  The subsequent exact source-delta ownership replay at `/tmp/qged-submission-delta-50k-20260826` passes in 95 seconds with 7,448 progressive payloads, 1.50M faces, 3.000-pixel certified error, and the same box-free terminal contract.  The exact-current warm OSMesa wire lifecycle remains green in 74 s. | warm shaded/wire OSMesa correctness, liveness, and interaction gate green; remaining cold/hardware rows remain |
| 150k scale | The bounded cold crash gate passed under a 16 GiB address-space cap on System GL and OSMesa with all 150,001 leaves discovered and no terminal boxes or failures.  An exact trace then exposed and fixed an ownerless 5,160-box structural frontier.  The warm OSMesa shaded matrix at `/tmp/qged-150k-superseded-fix-20260826` passes in 101 seconds after the deadline-recovery and superseded-result changes, reaching a zero-box, zero-failure terminal performance-limited view with 2,420 mesh payloads and 1.81M faces.  After folding structural-repair activity, frontier, and reservation into one typed transaction, `/tmp/qged-structural-repair-owner-150k-20260826` also passes in 101 seconds and terminates ready with 2,602 progressive payloads, 1.74M faces, 1.000-pixel certified error, zero boxes/failures, and no pending work.  These software-limited endpoints prove liveness/resource behavior, not realistic visual-significance quality. | shaded/wire warm OSMesa liveness green; cold/wire, production-hardware timing, and realistic visual-significance qualification remain |
| System GL smoke | The current shaded Generic Twin cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 12/11 seconds after camera identity became an optional snapshot.  Both exact-view returns recall history and every terminal checkpoint is ready with zero structural boxes and occurrence failures; the final certified errors are 0.236/0.241 pixels.  The 2026-08-27 post-ownership full System GL interaction replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-system-report.json` passed in 9.7 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The prior wire evidence is `/tmp/qged-generic-wire-system-debug-20260824`; both terminal wire frames contain 709 mesh payloads, zero boxes/failures/pending work, and matching cold/warm camera state.  The exact-current Lucy shaded cold/warm interaction replay at `/tmp/qged-lucy-system-20260824` also passed: both final frames are exact/ready, box-free, and quality-floor-clean with one resident source payload, 7.41M displayed PoP faces, a 230 MB resident mesh set, and one quality-history recall.  Hubble shaded/wire cold/warm pass after the overview lifecycle repair. | Qualify Lucy wire and complete the System-GL large-model matrix |
| OSMesa Generic Twin | The current shaded cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 15/14 seconds.  Both exact-view returns recall history and terminate ready with zero boxes/failures; the final certified errors are 0.621/0.676 pixels.  The 2026-08-27 post-ownership full replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-osmesa-report.json` passed in 14.3 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The final scene content and camera match the paired System-GL images.  The 2026-08-26 cross-renderer wire replay at `/tmp/qged-generic-wire-hud-final-20260826` passes cold and warm on System GL and OSMesa after the availability/HUD changes, with no invalid ready-label sample or terminal box. | focused regression green; continue the larger real-model matrix |
| spatial Lucy, System GL | Exact-current 2026-08-25 warm shaded retained-page replay passes.  It refines during continuous zoom, reaches its quiet pixel target, compacts on zoom-out, restores the prior view through quality history, and has no boxes or page-level subpixel proxies. | focused retained-page regression green; cold and wire rows remain |
| direct-mesh Generic Twin | Exact-current 2026-08-25 shaded cold/warm replays pass on System GL and OSMesa.  Each terminal frame has 708 direct full-detail payloads, one progressive payload, all 709 BoT occurrences, zero structural boxes, and no pending work.  The overall extent preview remains visible during discovery; it is not counted as a semantic leaf. | focused admission regression green; keep discovery-preview latency distinct from terminal-mesh admission |
| Hubble OSMesa | Exact-current 2026-08-25 shaded cold/warm lifecycles and camera contracts pass after the static-quality restoration change; the prior shaded/wire evidence remains green, with 1,804 warm shaded payloads and zero structural boxes. | focused regression green; rerun wire with the final binary and continue the full model matrix |
| multi-Lucy/xpush capacity | Exact-current 2026-08-26 warm initial-view runs pass on both renderers after capacity samples were ordered behind ceiling-free occurrence-plan handoff.  The pre-fix timing replay split 2/4 between pixel demand and a premature ordinary-deadline endpoint; after the applied-but-hidden candidate fix, 6/6 early-checkpoint OSMesa runs reached the identical approximately 984.1k-face, 0.457-pixel pixel-demand endpoint in 4.59--7.04 s with no boxes or control work.  Post-recovery-ownership OSMesa replays remain green: multi-Lucy reaches 985,820 faces before a performance-limited zoom endpoint, while xpush reaches 985,808 faces initially and 1.53M after zoom, all box-free and terminal.  System GL reached 3.81M faces at 0.228 pixels in 2.45 s.  The 2026-08-27 full 75-event OSMesa turnover/zoom/rotation replay at `/tmp/qged-multi-lucy-growth-handoff-20260827/formal-contract-report.json` completed in 29.2 s after completed-pass ownership and annotation lifetime were made atomic.  All waits settled terminal/ready with zero failures or control violations; the final view had eight meshes, 319k faces, zero boxes, and no remaining owner or obligation. | warm full-interaction regression green; cold eight-distinct-asset preparation, compaction, and wire rows remain |

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
