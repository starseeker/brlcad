# Obol production-readiness matrix

Last reviewed: 2026-08-23

This is the executable qualification plan for the BRL-CAD Obol drawing stack.
It records required rows and current evidence, not a chronological development
log.  A row is green only when it was produced by the exact final binaries and
its report, images, and performance data satisfy the common criteria below.

## Current release state

The semantic/API redesign and lower-level architecture gates are complete.
The exact current build passes the combined `drawing_baseline` and
`bobol_headless` selection (86 tests in 389 seconds, one unavailable X11
capability skip), including the full NIST BREP corpus, qged polygon and
primitive-edit replay, selection, picking, measurement, framebuffer,
endpoint, and shared-library contracts.  Prior current-architecture evidence
also includes the larger graphical control matrices and System-GL 50k/150k
stress runs.  The complete real-model and platform release matrix below
remains open.

Recently completed focused evidence:

- 2026-08-23 Lucy/large-scene convergence ownership: the focused
  `ObolLodConvergence.tla` model checks two interrupting input epochs for both
  workload profiles, exclusive progress ownership, pending-work witnesses,
  constrained quality debt, and eventual quiet stability.  TLC exhaustively
  generated 562 states (355 distinct) with no safety or liveness failure.  The
  source-level rule derived from the model caps quiet residency at the current
  scene allocation while allowing active scale interaction one bounded
  prefetch transition.  The previously stagnant certified-warm 150k OSMesa
  shaded prefix now completes in 32.2 seconds with exact 150,001-occurrence
  classification, 3,402 mesh payloads, 146,599 logical proxies, 228,211 active
  faces, about 21.2 MiB resident mesh data, zero structural boxes, and no
  pending work.  The matching Lucy smooth-zoom lifecycle completes in 27.1
  seconds, refines during active zoom, returns to 559,494 faces, and compacts
  resident mesh data from about 22.0 to 6.15 MiB while zoomed out before
  restoring it.  The full certified-warm OSMesa lifecycles also pass for 50k
  shaded in 60.6 seconds and 150k shaded/wire in 56.8/86.8 seconds, including
  interaction, selection, exact erase/redraw, exact occurrence coverage, zero
  boxes, and zero pending/stagnant work.  Repeat 50k wire, truly cold rows, and
  the complete System GL matrix before release.

- 2026-08-23 bounded software-proxy presentation: logical subpixel-proxy
  coverage remains per occurrence in `SoCADAssembly`, while the OSMesa
  renderer caches a deterministic camera-local screen-bin stream for physical
  submission.  A current 150k shaded warm interaction, selection, erase, and
  redraw row passed in 49 seconds with 148,619 logical proxies, 3,457 submitted
  proxy points, zero structural boxes, convergence fraction 1, and a 97.0 ms
  terminal frame.  The matching 50k OSMesa shaded cold/warm rows passed in
  51/40 seconds with roughly 46.7k/46.9k logical proxies and roughly 3.29k
  submitted points.  The new qged metric and scale guard keep semantic coverage
  distinct from physical draw cost.  This is a focused renderer qualification;
  rerun the complete OSMesa and System GL real-model/scale matrices, including
  wire, dense overlap, selected color, HiDPI, and resize, before release.

- 2026-08-23 pending allocator requalification: the scene allocator now admits
  modest point-eligible occurrences jointly with mesh cuts, rather than making
  every 12-pixel part a mandatory mesh before visually dominant parts can
  refine.  A one-occurrence scene is deliberately excluded: point aggregation
  cannot reduce its dispatch work and must not hide its only mesh.  A fresh
  OSMesa Lucy cold start reaches its first stable checkpoint in 26.3 seconds
  with one compact resident mesh, 144,976 faces, no proxy, and no compact
  registry gap.  The previous full run found that exact handoff gap before the
  correction; its other 15 rows passed.  The focused 50k OSMesa matrix passes,
  while the whole OSMesa matrix and 150k rows must be repeated before release.

- 2026-08-23 static marginal-capacity correction: a hard-deadline miss already
  records a five-percent-below-failed strict capacity ceiling.  The retained
  marginal allocator was applying its independent 20-percent throughput margin
  to that already conservative ceiling a second time, unnecessarily stranding
  visible detail.  The current source caps the once-margin-derived estimate at
  the strict ceiling directly; `libBObol_lod_coordinator` exercises that
  boundary.  Focused certified-warm 150k OSMesa shaded runs reach an idle,
  zero-box terminal state after initial draw and zoom/rotation (roughly
  34--71 seconds on this host, with 80--98 ms final frames).  A later full
  interaction/selection/erase/redraw replay completed in 55 seconds with
  2,806 mesh payloads and 149,135 subpixel proxies: a useful 150k overview
  is feasible, while a per-occurrence software-raster presentation does not
  meet the interactive target on this host.  This validates liveness and the
  bounded recovery path only: approximately 1.2k of 2.45k prominent
  occurrences can still miss the visual-quality floor.  Do not mark
  visual-significance qualification green until the significance-frontier work
  below has improved a representative real-model result as well as this
  synthetic stress case.

- 2026-08-23 current-source retry/occurrence qualification: after the static
  first-presentation retry and source-logical occurrence fixes, the full OSMesa
  qged matrix passes again for Generic Twin, Lucy, Havoc, and Hubble in
  shaded/wire cold/certified-warm (16 rows), plus 50k and 150k distinct-asset
  shaded/wire cold/warm.  The same source passes every selected libBObol,
  libged drawing, qged interaction/editing, qtcad/Obol, gsh, and MGED/LoD CTest
  row.  The full CTest invocation still has unrelated environment/reference
  failures: loopback-network tests are sandbox-blocked, and the pivot/repository
  guards plus legacy raytrace image references need their own owners.  System
  GL remains unqualified on this host because no usable X socket is available.

- 2026-08-23 exact-current OSMesa qualification: Generic Twin, Lucy, Havoc,
  and Hubble pass fresh-cold and certified-warm shaded/wire qged rows (16
  total), including screenshots, selection/tree transitions, exact erase and
  redraw, smooth zoom, retained-history return, and zero terminal boxes.  The
  50k and 150k distinct-asset rows also pass shaded/wire cold/warm.  The 150k
  final renders are bounded to roughly 34--86 ms and retain exact 150,001
  occurrence coverage, but report a typed software-performance limit and
  nonzero prominent-quality-floor debt; this is not release clearance for
  visual significance.  System GL must be rerun from this binary set on a
  valid X host: the present host's `/tmp/.X11-unix` owner prevents Xvfb from
  creating its local socket.

- 2026-08-23 current-source continuity qualification: a completed CAD work
  record now preserves the richest deadline-safe presentation cut across a
  policy-only interaction-to-quiet handoff, while camera, renderer, and
  source-quality-domain changes invalidate it.  Fractional terminal targets
  remain fractional through submit-action allocation.  Lucy shaded passes
  cold/warm on System GL (51/26 seconds) and OSMesa (55/32 seconds), with a
  0.75-pixel return request realized at cut 21 / 0.645 pixels; Hubble OSMesa
  shaded/wire cold/warm passes after a state-based first-canvas-paint barrier;
  and shared-cache OSMesa shaded scale rows pass at 50k (24 seconds) and 150k
  (53 seconds).  The 150k terminal state has zero boxes, 2,066 mesh
  occurrences, 147,934 subpixel points, and a 60 ms retained frame.  This is
  focused post-change evidence only: repeat the full shared-stack and truly
  cold scale rows before release.

- independently compiled controller frame-execution and libged retained-
  geometry bridge extractions, followed by current System GL and OSMesa
  Generic Twin shaded/wire cold/warm graphical passes (all 709 managed
  occurrences are meshes with zero boxes in the wire rows);
- current System GL and OSMesa Hubble shaded/wire cold/warm passes, including
  tree selection and exact erase/redraw through the representative `PANEL_C01`
  branch;
- current Lucy shaded and wire cold/warm passes on both renderers.  The
  graphical smooth-zoom sequence proves in-motion spatial-prefix growth,
  crack-free close views, bounded wheel dispatch, exact-view history recall,
  zero terminal boxes, and zoom-out resident compaction.  OSMesa wire
  correctly reports a hard static-frame performance limit when the next
  discrete population cannot meet the deadline rather than falsely claiming
  unconstrained pixel-exact convergence;
- current Havoc shaded/wire cold/warm passes on both renderers;
- current 50k distinct-asset System GL shaded cold/warm passes (21/19 seconds
  wall time for the scripted rows), with 50,001 available leaves, exact
  terminal coverage, zero boxes, roughly 15.3k retained mesh payloads plus
  36.6k subpixel occurrences, and a separately passing warm perf capture;
- current 150k distinct-asset System GL shaded cold/warm passes (42/39 seconds
  wall time), with exact 150,001-leaf terminal coverage, zero boxes, roughly
  15--18k retained mesh payloads plus 135--137k subpixel occurrences, and peak
  transient mesh work below 68 MiB;
- current 50k distinct-asset OSMesa shaded cold/warm passes (35/50 seconds),
  and a separately fresh-cache OSMesa wire cold/warm pair (50/75 seconds).
  The wire rows converge after rotation with exact 50,001-occurrence
  classification, zero structural boxes, roughly 1.8k retained mesh
  occurrences plus 48k subpixel occurrences, an explicit software-performance
  limit, and zero stagnant coordinator dispatches;
- current 150k distinct-asset OSMesa shaded cold/warm passes (98/72 seconds).
  Both terminal states have exact occurrence classification, zero boxes, zero
  stagnant dispatches, an aggregate point threshold of 64, and explicit
  software-performance limits.  The cold/warm frames retain 2,119/2,855 mesh
  occurrences plus roughly 148k subpixel occurrences.  A 297-frame APNG
  diagnostic distinguishes one bounded capacity backoff from the former
  coarse/fine alternation, and the final stable screenshots are coherent;
- focused current-build warm revalidation after structural marginal-budget
  certification and the closed interrupted-replay transaction reaches exact
  150,001-occurrence classification with zero boxes.  Clean 150k shaded runs
  completed in 25--26 seconds total (16--18 seconds in the final idle wait),
  retained roughly 1.6k meshes plus 148k subpixel occurrences, and drew the
  terminal OSMesa frame in roughly 82--90 ms.  The corresponding System GL row
  completed in 16 seconds total (8 seconds in the idle wait) with a 25 ms
  terminal render.  A deliberately traced OSMesa run took 53 seconds and a
  clean timeline took 68 seconds total while still finishing with zero boxes;
  this remaining preparation/allocation variance is an open release issue, not
  evidence that the terminal 150k software presentation is intrinsically
  undrawable.  A current 50k OSMesa replay also finished in 15 seconds total,
  with exact coverage, zero boxes, no stagnant dispatches, and an 85 ms
  terminal frame;
- current 150k distinct-asset OSMesa wire warm interaction/lifecycle pass
  (126 seconds), with exact 150,001-occurrence classification, zero structural
  boxes, no pending work, and a ceiling-free terminal local allocation after
  zoom, rotation, selection, exact subpath erase, and redraw.  This proves the
  resident-growth and failed-headroom no-reentry fixes, but is not yet a green
  release row: cold evidence is outstanding and repeated warm runs have used
  materially different portions of the 100 ms static allowance (roughly
  18--60 pixels worst projected error);
- rotation- and aspect-safe autoview framing, exact retained-scene image
  controls, and matching eager/deferred terminal endpoint images;
- repeated warm System GL validation of the level-triggered Qt idle-tail
  witness after fixing a coalesced-update lost wake;
- fixed-width source population result authentication;
- non-destructive indexed cache inspection and bounded cache invalidation;
- parallel read-only database discovery parity on the focused librt fixture;
- cost-weighted progress on an actual 150k compact scene;
- shared BRL-CAD/libBObol/Obol ASan graphical startup and active-worker
  shutdown; and
- normal source-realization, draw-cache, and mesh-cache tests.

## Evidence rules

- Record source and Obol revisions, build configuration, renderer, model
  identity, canonical cache path, cold/warm proof, event script, viewport/DPR,
  elapsed checkpoints, terminal diagnostics, and artifact hashes.
- A timeout, crash, sanitizer finding, invariant violation, empty capture,
  missing report, unexplained terminal box, or manual-only observation fails
  the row.
- Cold uses a new empty cache namespace.  Warm reuses only a successful cold
  completion marker for the same canonical file identity.
- Fixed sleeps are not completion.  Wait on controller, subprocess, frame, or
  quiet-state predicates with an outer timeout.
- Keep ordinary timing at stage granularity.  Per-batch and per-asset logging
  is verbose-only: logging a 50k/150k stream changes the owner-thread schedule
  that the trace is meant to measure.
- Preserve small reports and representative failure images.  Reuse one copy of
  large geometry and delete per-run caches, traces, and duplicate captures once
  findings are recorded.

## Common pass criteria

Every applicable row must prove:

- correct final camera, model extent, draw mode, colors/materials, lighting,
  selection/highlight, faceplate, and framebuffer composition;
- useful overview/coverage during cold startup and no unexplained blank period;
- monotone useful refinement without mesh-to-box regression, zero-work cycling,
  or a coarse restart after ordinary motion;
- no terminal fallback boxes unless a named source failure remains;
- terminal convergence fraction 1 with no queued source/task/result,
  publication, allocation, compaction, timer, or unacknowledged frame work;
- constrained states explicitly report performance, memory, or terminal source
  limits and retain complete visibility classification;
- bounded input/event latency and responsive interruption of expensive frames;
- memory remains within configured resident and transient limits and recovers
  after zoom-out/non-visible compaction; and
- System GL and OSMesa agree semantically and remain within the row's explicit
  image tolerance.

## Model corpus

| Model | Primary purpose |
|---|---|
| `db/moss.g`, rook | modes 0-5, evaluated CSG, labels and editing |
| `db/faa/Generic_Twin.g` | production visual baseline, mixed geometry, tail/hull/engine detail |
| `lucy.g` | one very large mesh, spatial pages, zoom history and memory recovery |
| multi-Lucy and xpush | shared assets, PCA/transform reuse, occurrence turnover |
| Stanford mesh assembly | varied mesh density and hierarchy reuse |
| Hubble | wide hierarchy, many small components, selection/tree/importance |
| `havoc.g`, independent vehicle models | deep real hierarchy and mixed-size visual significance |
| NIST BREP corpus | wire/shaded BREP, adaptive tessellation and cache reuse |
| 5k/50k/150k varied fixtures | unique assets, regions/colors/depth, aggregate budgets |

Synthetic scale fixtures must vary object sizes from roughly 500 to 100,000
faces, include at least 95% manifold inputs, mix regions, colors, hierarchy
depth, shaded and wire modes, and contain visually prominent components which
occupy a meaningful fraction of the overall vehicle extent.

## Lower-level and shared-stack gate

- clean Ninja build of qged, MGED, gsh, Archer/TkObol, rtwizard, and plugins;
- `bobol_contract`, `bobol_headless`, `drawing_baseline`, libbg polygon,
  librt discovery/edit, libbu cache, and installed-package consumer tests;
- installed C11/C++17 header self-containment and public symbol manifests;
- shared-stack ASan/UBSan including active worker teardown, cache corruption,
  endpoint replacement, plugin reload, edit cancellation, and rapid view close;
- TSan/LSan on a compatible native worker;
- focused TLC checks for host work, occurrence publication, and bounded
  Lucy/large-scene convergence ownership.

## qged graphical control matrix

Run every row in single and quad layouts where applicable, System GL and
OSMesa, DPR 1 and fractional DPR, before and after resize.

### Camera and view

- MGED-style orbit, rotate, translate, center, zoom, smooth zoom round trip,
  autoview, preset `ae`, and orthographic/perspective behavior;
- nav gizmo cube and circles/axes styles, mouse clicks/drags, enable/disable,
  and command/widget readback;
- startup, settled, mid-refinement, maximize, fullscreen, restore, minimize,
  hidden/exposed, and quad-layout resizing;
- Studio and MGED lighting profiles, material/culling transitions, and matching
  System GL/OSMesa illumination.

### Selection and hierarchy

- point, rectangle, append, toggle, subtract, select-all/clear, and mode changes;
- exact selected occurrence, active ancestors, impacted descendants, parent and
  grandparent whole-row styling;
- selection before realization; draw, exact erase, subtree erase, ancestor draw
  subsumption, redraw, rename/delete, and edit while selected;
- sole-selection edit geometry versus multi-selection styling only;
- Hubble large rectangle latency and tree scrolling/repaint cost;
- Coin pick and application/librt callback pick parity for mesh-backed targets.

### Polygon editing

- rectangle, square, circle, ellipse, and general polygon mouse creation;
- constrained resize, whole movement, point select/drag/append/delete;
- union, intersection, difference, holes, islands, and open contours;
- grid/line/geometry snapping, fill/edge styles, copy/import, sketch round trip,
  persistence, Escape/button/tool-handoff cancellation, and deletion;
- CLI-to-widget and widget-to-CLI synchronization plus shaded/wire CAD
  composition while progressive drawing is active.

### Primitive and sketch editing

- generated and custom controls for representative primitives;
- command, widget, and retained manipulator changes reflected in all views;
- invalid-value immutability, set/readback, commit, cancel, checkpoint/revert,
  rename/remove, erase under edit, alternate occurrence conflict, and plugin
  unload/reload;
- ARB vertex/edge/face selection and handles, BoT picking surface, constrained
  axes/planes/rings, and selection-driven edit promotion;
- sketch direct curve creation, segment/vertex/topology edits, NURB/arc modes,
  non-unit planes, linked extrude/revolve profiles, and persistence.

### Measurement and remaining tools

- drive every measurement mode by mouse, including creation, adjustment,
  labels, units, snapping, selection/highlight, readback, cancellation, and
  cleanup;
- grid, ADC, model/view axes and ticks, faceplate labels/lines, snap settings,
  view information/settings, export, hierarchy filters, and plugin lifecycle;
- embedded raytrace framebuffer underlay/overlay/interlay/off, subprocess
  completion, resize/rerender, selection/edit composition, and CAD restoration.

## Drawing and LoD matrix

- Modes 0-5 on assigned moss/rook models with LoD auto and off.
- Generic Twin shaded/wire, cold/warm, LoD auto/off, System GL/OSMesa.  Compare
  `ae 90 0` and additional oblique/tail views with the main-branch production
  baseline.  Inspect hull gaps, engine, tail, underside lighting, and boxes.
- NIST BREP modes 0, 1, 2, and 4, cold/warm and LoD auto/off on both renderers.
- Lucy cold/warm smooth zoom in/out and rotation.  Verify page coverage,
  crack-free wings/top/bottom, prefix history, in-motion refinement, uniform
  zoom-out detail, and resident memory recovery.
- Multi-Lucy/xpush view turnover: bring richer and unseen instances into view
  while others leave, preserving shared assets and interaction FPS.
- Hubble shaded/wire resize, selection, erase/redraw, rotation restoration,
  visual-importance floors, and stable subpixel representation.
- LoD-off rows may be slow but must remain correct and interruptible.

## Scale and performance matrix

Run 5k, 50k, and 150k varied fixtures from cold and certified warm caches.
For 50k and one real model capture perf; for observed flicker capture APNG; for
System-GL-only corruption capture apitrace.

Required checkpoints:

1. database directory/reference discovery;
2. whole-target overview;
3. complete leaf/subpixel coverage;
4. first useful meshes;
5. initial ready state;
6. smooth zoom in/out;
7. pure rotation and translation release restoration;
8. visibility turnover;
9. selection and tree update;
10. exact/subpath erase and redraw;
11. final constrained or unconstrained ready state.

Perf review must separate production work from diagnostic JSON, path dumps,
PNG/APNG compression, and capture overhead.  Investigate any O(scene) owner-
thread pass repeated by a bounded allocation window or sparse capacity retry.

## Client matrix

- qged: complete control matrix above.
- MGED: shared LoD, modes, axes/ticks, sed/oed, framebuffer composition,
  selection and resize under TkObol.
- gsh: headless draw/erase/cache/status/LoD/edit semantics and clean shutdown.
- Archer: four-view lifecycle, selection, resize, RtControl framebuffer, and
  database replacement.
- rtwizard: embedded progressive result publication, resize, completion, and
  teardown.

## Platform matrix

Linux System GL and OSMesa are local release gates.  Windows must independently
qualify WGL/Qt/Tk, fixed-width caches, BREP/PoP zoom, transient-vertex safety,
memory reporting, and sanitizer behavior.  macOS requires a selected native
context provider and equivalent host/lifecycle coverage.

## Final audit

After the final fix:

1. rebuild every touched target and rerun its focused group;
2. rerun all lower-level labels and client smokes;
3. execute the full graphical and scale matrix from one recorded binary set;
4. compare against main-branch visual ground truth and previous large-model
   baselines;
5. verify no required row is stale, missing, skipped, manually inferred, or
   using an uncertified cache; and
6. update `libbobol_active_debt.md` with only genuine platform or optional
   capability deferrals.
