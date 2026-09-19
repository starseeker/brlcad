# Obol production-readiness matrix

Reviewed 2026-09-19. This document owns acceptance criteria and evidence status.
The [guide](obol_simplification_guide.md) owns the roadmap and
[active debt](libbobol_active_debt.md) owns implementation work. Production
readiness remains open; neither historical passes nor test counts qualify the
current checkout.

## Current qualification status

| Gate | Status | Evidence needed to close |
|---|---|---|
| S0 reproducible stack | Open: configuration can replace dependencies or lose package discovery | Fresh isolated build, repeated configuration and verified runtime identities |
| S1 supported surface/ownership | Open: viewport mutation boundary narrowed; other families need inventory | Complete caller/entry/owner/contract mapping |
| S2 complete production workflow | Open: focused repairs and framebuffer round trips only | Cold CAD draw, active camera, cancellation/close, reopen/redraw with output and resource evidence |
| S3 remaining consolidation | Open | Delete superseded owners/state and migrate callers |
| S4 resources/lifetime | Open | BREP, contention, measured resource bounds and remaining lifecycle/control cases |
| S5 capability envelope | Open; known visual failures and unqualified native/client/edit rows | All required graphical, geometry, scale, interaction and platform rows |
| S6 release candidate | Open | All required evidence tied to one exact final candidate |

The September 19 assessment ran 14 selected existing CTest rows successfully
and reproduced a raw viewport visibility round-trip failure. Those binaries
predate the transition and are not qualification for its changed source.
Assessment artifacts are under `.build/obol-qualification/20260919-assessment`.
The transition's build/test records belong under
`.build/obol-qualification/20260919-simplification-transition`; consult its
report and the [handoff](obol_session_handoff.md) for exact tested scope.

The rebased baseline's focused build, behavioral and GUI results are recorded
in [main integration](obol_main_integration.md). Its framebuffer transition
checks pass. Annotation visual parity, an intermittent scene-light sweep
assertion, clean native dependency reproduction and the broader acceptance
matrix remain open. These results do not close S0--S6 by themselves.

Prior publication repairs through CXX-SOURCE-089 remain retained regressions.
The unsealed direct-viewport-field experiment is withdrawn in favor of owning
host/stream mutations; it is not recorded as a new conformance closure.
Earlier stage percentages are retired because no stable weighted acceptance
denominator existed. Report named criteria and candidate-specific evidence.

The full [September 19 readiness snapshot](obol_20260919_readiness_history.md)
preserves historical results and counterexamples. Old `/tmp` reports and warm
caches cited there are unavailable; reconstruct them before relying on their
claims. The sole formal baseline is `tla/baselines/tlc-2.19.json`; model state
counts prove only the modeled finite protocol and do not qualify native C++.

## Common acceptance criteria

Every graphical row must show:

- correct final camera, extent, draw mode, colors/materials, lighting, and
  requested normal style across ordinary and spatial-page meshes, and retained
  scene composition;
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

The full database corpus must cover every named solid and combination,
including nested roots. Expensive volumetric and repeated aggregate cases
need declared time bounds, not silent omission. Preserve finite-output and
bare-root mesh/coverage regressions alongside these model workflows.

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

For each required row record source revisions plus dirty changes, dependency
and loaded-library identities, build options, input identity, cache condition,
resource limits, command/replay, result and representative inspected images.
Distinguish native GPU, software System GL/llvmpipe and OSMesa. A missing limit
or unavailable required result leaves the row open. Declare image tolerances
and latency/memory/convergence thresholds before a qualification run; do not
choose them after observing its output.

Store current evidence under a persistent qualification directory with a
manifest. Keep a shared copy of reusable large inputs and certified warm
caches. Retain reports, scripts, meaningful images and diagnostic traces;
remove duplicate transient captures only after recording the useful evidence.
Historical paths remain historical and must not be labeled current passes.
A dependency or implementation change requires re-running affected rows.
