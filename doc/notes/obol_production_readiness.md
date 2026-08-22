# Obol production-readiness matrix

Last reviewed: 2026-08-22

This is the executable qualification plan for the BRL-CAD Obol drawing stack.
It records required rows and current evidence, not a chronological development
log.  A row is green only when it was produced by the exact final binaries and
its report, images, and performance data satisfy the common criteria below.

## Current release state

The semantic/API redesign and lower-level architecture gates are complete.
The exact current build passes the combined `drawing_baseline` and
`bobol_headless` selection (86 tests, one capability skip), including the full
NIST BREP corpus, qged polygon and primitive-edit replay, selection, picking,
measurement, framebuffer, endpoint, and shared-library contracts.  Prior
current-architecture evidence also includes the larger graphical control
matrices and System-GL 50k/150k stress runs.  The complete real-model and
platform release matrix below remains open.

Recently completed focused evidence:

- independently compiled controller frame-execution and libged retained-
  geometry bridge extractions, followed by current System GL and OSMesa
  Generic Twin cold/warm graphical passes;
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
- focused TLC checks for host work and occurrence publication.

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
