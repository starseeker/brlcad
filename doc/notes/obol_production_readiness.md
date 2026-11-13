# Obol production-readiness matrix

Reviewed 2026-09-20. This document owns acceptance criteria and evidence status.
The [guide](obol_simplification_guide.md) owns the roadmap and
[active debt](libbobol_active_debt.md) owns implementation work. Production
readiness remains open; neither historical passes nor test counts qualify the
current checkout.

## Current qualification status

| Gate | Status | Evidence needed to close |
|---|---|---|
| S0 reproducible stack | Demonstrated: September 19 source-built Linux candidate | Repeat affected checks when source/dependency selection changes; final packaging/platform qualification remains S5/S6 |
| S1 supported surface/ownership | Demonstrated on the September 20 checkout: finite inventory and production caller migration | Re-run contract checks when a live mutation family or public owner changes; final candidate remains S6 |
| S2 complete production workflow | Demonstrated on the September 20 Linux software-rendered checkout | Repeat the named FLOW-01 workflow when its owners change; native GPU/platform and sanitizer lifetime qualification remain S4--S6 |
| S3 remaining consolidation | Demonstrated on the September 20 checkout: three bounded seam closures and final consumer audit | Reopen only for a concrete duplicate writer, mirrored control state or split policy; final candidate remains S6 |
| S4 resources/lifetime | Open | BREP, contention, measured resource bounds and remaining lifecycle/control cases |
| S5 capability envelope | Open; known visual failures and unqualified native/client/edit rows | All required graphical, geometry, scale, interaction and platform rows |
| S6 release candidate | Open | All required evidence tied to one exact final candidate |

The base source-built Linux candidate is recorded under
`.build-main/obol-qualification/20260919-annotation-fills`; focused interaction,
legacy-fallback and immediate-GL evidence for the current binaries is under
`20260919-annotation-interaction`. The S1 owner/caller audit and focused
contract regressions are under `20260920-api-inventory`. The complete headless
and offscreen Qt production workflow is under `20260920-flow-01`. The native stack
uses Obol `a544ae3876bd` plus its captured local patch and clean OSMesa
`bfa90de454a7`; the BRL-CAD build uses the normal build-local install. The exact
source patch, headers, build settings, binaries and loader identities belong
to the candidate manifest. Rebase and preceding repair evidence live in the
[integration record](obol_main_integration.md), rather than defining additional
acceptance gates here.

Nine native suites, four focused drawing checks, 39 affected integration rows
and both library/Qt export tests pass. Both publication sweeps also pass, retaining the allocation and callback-order
assertions. A full BRL-CAD build (1,442 Ninja steps), the isolated installed-package
consumer, repository/license checks and both `annotate` rows pass after the
main integration repairs. Retained protocol-model results belong to the preceding
projection candidate; this numeric/representation change does not alter the
modeled publication protocol.

The repaired behavior includes per-source cancellation, material/update
preservation, overview-to-leaf publication, independent camera projection,
authored stroke width/color/pattern and retained unlit fill geometry. Focused
native pixels exercise fixed VBO, GLSL VBO, instanced, flat-batch and direct
software wire rendering and cover zoom, base-width changes, occurrence color
replacement and same-ID geometry replacement. A diagnostic capability override
now forces the immediate GL fallback on compatibility-profile OSMesa and
System-GL contexts; the same pixel test verifies tier-zero execution, authored
styles and replacement separately from direct software wire.
Software System GL under Xvfb is not native GPU qualification.

Fill-specific native pixels and picking pass in all four native draw modes,
including holes, authored color/replacement, mixed point/fill/stroke geometry,
display-plane projection and same-ID role replacement. Source export retains
area, fill color, semantic background masks and separately authored strokes,
including exact line masks. Mask pixels match the active solid or gradient
background across fixed and GLSL paths, model/display coordinates and selection.
Existing triangle storage, GPU buffers, the unlit point pass and the render
action's frame background own the repair; no additional scene or publication
controller was introduced.

All seven annotation comparisons initially failed the unchanged 0.99
threshold. Follow-up [image attribution](obol_main_integration.md#annotation-image-attribution)
under `20260919-annotation-framing` and `20260919-annotation-framing-1-6`
isolates the old bounds policy, former material lookup, fixed-pixel
display-plane placement and newly retained authored styling. Fitted old cameras
raise frames one through six to 0.961952--0.975646; their model-space annotation
bounds and exact-color pixels align within one pixel. Frame seven passes at
0.999907 under its independently derived old camera. New command assertions
protect representative materials, radial dimension spans, centered complete
framing and display-plane projection. With user authorization, all seven
controls now reflect the stable current render and pass at the unchanged 0.99
threshold. The drawing-image, command-only and Qt export tests pass together.
Production rendering is unchanged from the output-consumer candidate; broader
qualification is retained from the fill manifest, not rerun for these test
additions. A focused follow-up makes snap and measurement use the active
display-plane projection for live viewport traversal and makes the legacy
vlist publisher reject screen annotations instead of treating display offsets
as model units. Full editing/client parity remains open in EDIT-01 and S5. The
former load-sensitive scene-light witness is closed by the later kind-aware
semantic equality repair and deterministic unused-payload regression. S3 is
now demonstrated by the later ownership closures; S4--S6 still require
resource, quality, client and platform evidence.

FLOW-01 adds named headless and graphical workflows. The headless path uses a
real delayed mesh request, changes the camera, closes its production endpoint
with active workers, reopens and redraws, then requires a changed terminal image,
compact source state and empty service queues. The graphical path performs the
same ownership sequence through a real `QgView`; a bounded service-delayed task
isolates widget destruction from small-model planner timing. Both paths measure
Linux worker retirement. Each passed ten consecutive runs, and the combined
ordinary, delayed, stale-source and denial coverage passed after the full build.
This demonstrates S2 on software OSMesa; it is not native GPU, sanitizer,
large-model or cross-platform qualification.

The first S3/OWN-01 consolidation is recorded under
`20260920-own-01-service-work`. The LoD service now owns lock-consistent
service-wide and per-generation work snapshots, the complete transient-idle
predicate and the generation result-work predicate. The latter includes a
consumer's lease on a shared producer even when that generation owns no task.
GED and GUI production waits no longer maintain partial counter lists, and
controller/renderer readiness no longer combines separately locked generation
counters. The controller's duplicate result-pending atomic is also gone; queue
availability comes from the service snapshot, while the independent first-ready
timestamp continues to bound batching latency. Thirteen focused service,
coordinator, API, workflow and progressive rows, six broader controller/model
rows, compact-publication stress and the installed consumer pass; both
production lifecycle rows also pass ten consecutive runs. This closes one
ownership seam while the remaining S3 inventory and S4--S6 qualification
obligations remain open at that checkpoint.

The second S3/OWN-01 consolidation is recorded under
`20260920-own-01-resident-capacity`. Stable resident bytes are directly
published as one scalar by every accounting path. The service samples that
value with its mutex-owned growth reservation and resident limit, and publishes
completed growth before releasing its reservation. Renderer headroom and CPU
pressure no longer assemble capacity policy from three independent calls;
convergence and qged diagnostics use the same contract. Unit arithmetic and
actual publication, compaction, eviction/reload and shutdown transitions pass.
The publication test observes a live full-hierarchy growth reservation while
changing its limit and rejects any zero-capacity handoff gap. The thirteen
focused rows, six controller/model rows, compact-publication stress, installed
consumer and ten repetitions of each production lifecycle path pass. The
remaining S3 inventory and S4--S6 qualification obligations remain open at
that checkpoint.

The third S3/OWN-01 consolidation is recorded under
`20260920-own-01-progress-display`. Qt no longer mirrors libged's decisions
about whether LoD progress is visible, when it is terminal-ready, or which
internal phases are one display transition. The convergence snapshot now owns
one immutable classification consumed by both paths. Libged continues to own
the faceplate's text, colors and geometry; Qt continues to own publication
cadence and scheduling. A direct classification table, ten focused faceplate,
progressive and API rows, the thirteen production-focused rows, six
controller/model rows, the installed consumer and ten repetitions of each
production lifecycle path pass. No annotation pixels or control images changed.

The closing S3 audit is recorded under `20260920-own-01-s3-close`. GED's
service wait and the GUI qualification wait now consume the existing immutable
host-work level rather than rebuilding it from separate controller getters;
service worker/cache quiescence remains a separate service snapshot. The audit
classifies the remaining synchronized state as single-writer boundary transfer
or diagnostics. The affected qged build, ten focused host/API/GED/Qt rows and
ten repetitions of GED cross-run plus both production lifecycle rows pass. S3
is demonstrated for the current inventory. S4--S6 remain open.

The subsequent segment/fill color repair leaves annotation frames one through
five and seven byte-identical to the framing candidate. Frame six alone changes:
authored colors make previously black-on-black segments visible and improve its
SSIM from 0.886535 to 0.887488. Its differently framed control remained
unchanged at that focused checkpoint; this run did not qualify the
candidate-wide matrix.

The subsequent mask repair adds focused native and BRL-CAD source/export
evidence. A continued visual run is byte-identical to the post-color candidate,
so no annotation control image changes for this repair.

The subsequent output-consumer repair gives PostScript and PNG explicit
effective width, exact pattern/factor, color/alpha and annotation-fill behavior.
Semantic masks use each format's declared background. Plot maps patterns to its
five named line modes and explicitly degrades width and filled areas, which the
format cannot represent. Qt object queries retain exact parallel style vectors.
The isolated command regression and focused Qt test pass. A continued
seven-frame run has the same hashes and scores as the mask checkpoint, including
an exact frame-seven match, so that output-only repair did not warrant control
changes. The later camera diagnostic supplies the control attribution. This
closes the known output representation mapping; broader client visual parity
remains part of S5.

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
