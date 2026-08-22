# libBObol engineering lessons

Last reviewed: 2026-08-22

This document preserves durable lessons from the Obol drawing migration and
production shake-down.  It is not a chronological status log.  Each entry
records a failure pattern, its architectural cause, the adopted rule, and the
kind of regression that must preserve it.

## Identity and lifetime

### Raw addresses are not durable identity

A compact part-deduplication memo retained raw `PartGeometry *` addresses after
the weakly owned geometry died.  Allocator reuse let unrelated authored
geometry inherit a fallback-box PartId, producing a truthful but unserviceable
box-repair loop.

Rule: pointer fast paths require a live weak witness to the exact object.
Durable identities contain an owner, object identifier, and generation or
revision.  Force allocator reuse in regression tests.

### Dense indices require population authentication

A late LoD result could arrive after a compact registry replacement which
reused the same occurrence key and dense index.  Route and key alone could not
distinguish the old population.

Rule: asynchronous compact requests carry a fixed-width population epoch.
Publication authenticates the epoch before using an index.  The epoch is not
part of persistent content identity, so legitimate cache reuse remains intact.

### Fixed-width persistent topology

Platform-sized counts and aligned typed reads are unsafe cache formats.
Windows-only invalid vertices and crashes were consistent with width/alignment
assumptions that happened to survive Linux testing.

Rule: persistent topology uses deterministic fixed-width integers.  LMDB is
byte storage; copy values into aligned fixed-width objects and validate every
count, range, index, coordinate, normal, and cumulative cut before publication.
Cache incompatibility invalidates old data rather than invoking compatibility
logic.

### Process-wide workers impose destruction ordering

ASan found a source-realization worker reading the draw-cache binding map after
the map's function-local static destructor ran.  The coordinator itself was
still waiting for process-exit destruction.

Rule: all lazy non-trivial globals reachable by process-lifetime workers are
constructed before those workers start.  Reverse destruction must join workers
first.  Regress with an active cache callback across return from `main`, and
test the normal shared-library stack rather than a duplicate static harness.

## State and liveness

### Work levels need progress witnesses

Toolkit notifications are edges and may coalesce; controller obligations are
levels.  Clearing an edge without retaining a timer, scheduled owner-thread
loop, in-flight task, cursor, or frame can strand work.  Treating every repaint
as work hides this defect and wastes frames.

Rule: every nonterminal state names its progress witness.  A completed frame
acknowledges only the revision it captured and cannot clear a newer request.
The focused `ObolHostWork.tla` model and window-host tests guard this contract.

### A zero-work repair is a state-machine bug

Several cycles alternated balancing/refining or repeatedly rendered boxes after
all useful tasks drained.  Causes included stale visibility sources, ordinary
structural visuals mislabeled as LoD fallbacks, an old part remaining bound
behind an unchanged numeric revision, and pressure revisions with no eligible
work.

Rule: only explicitly marked unresolved leaf fallbacks create a box-repair
obligation.  A retry must identify a task, cursor, budget edge, handoff, or
render acknowledgement that can change state.  Otherwise retire it as terminal
or error.  Tests must reject mesh-to-box regression and repeated phase cycling.

### Ready is a semantic proof, not an empty queue

An empty worker queue may still have unpublished results, unacknowledged frames,
partial visibility census, pending resident allocation, or source preparation.
Conversely, a constrained scene can be validly ready without all original
geometry resident.

Rule: readiness combines exact live-source coverage, visible occurrence
classification, affordable cut satisfaction, terminal failure/constraint
reasons, drained publication, and acknowledgement of the current frame.

### Progress percentages must reflect remaining cost

Counting only represented leaves reported roughly 95--99% while expensive mesh
handoff and refinement still remained.  Users interpreted the nearly full bar
as a stalled system.

Rule: progress is phase- and cost-weighted: discovery/coverage, useful
representation, detailed mesh work, and final handoff have separate weight.
The label reports counts and the terminal constraint; it is an estimate, not a
promise based on object count alone.

## Progressive mesh policy

### Cut ordinals have no geometric meaning

The original sixteen lock-step levels produced large population jumps and
encouraged assumptions such as “15 means full.”  Those assumptions spread into
renderer, cache, picking, and scheduler code.

Rule: the producer publishes an admissible cut vector.  Each cut has cumulative
point/face/byte counts, axis precision, conservative object error, and an exact
flag.  Consumers binary-search physical metadata and use a named out-of-range
sentinel for “richest resident.”

### Pixel exact is not source complete

Loading and drawing every source triangle is wasteful when the view cannot
resolve it.  Large stable scenes should spend more than the interactive budget
when useful, but need not retain invisible suffixes.

Rule: stable demand is the richest visible pixel-justified cut affordable
under finite static frame and memory budgets.  Retain the rendered framebuffer
until something changes.  Under memory pressure, compact non-visible suffixes
without deleting semantic intent or disk-cache identity.

### Existing prefixes are the starting point

Zoom and motion originally restarted coarse progression, causing boxes,
flicker, and slow restoration.  PoP data is cumulative; a new cut normally
changes array ranges rather than rebuilding the hierarchy.

Rule: rotation and translation retain the proven cut unless measured FPS
requires a temporary ceiling.  Orthographic pure rotation restores the prior
scale allowance immediately, then resolves newly visible silhouettes.  Zoom
retargets from current resident prefixes and may load suffixes during input.
Only missing data uses a box.

### One global prefix is insufficient for a giant projected mesh

Lucy demonstrated that a single global cut may spend triangles outside the
visible high-detail area and remain coarse where zoomed.  Spatial pages make
subsets selectable while sharing one canonical hierarchy and content identity.

Rule: page/chunk selection is producer metadata over one asset, not independent
semantic mesh objects.  Page boundaries must remain crack-free and active cuts
must be uniform enough to avoid holes or visible section discontinuities.

### Visual importance is scene-wide

Per-object pixel targets alone can leave prominent wheels, blades, tails, or
hulls coarse under aggregate pressure.  Depth is not an orthographic metric.

Rule: allocate a calibrated scene budget by projected visibility, footprint,
screen-error reduction, silhouette/surface significance, selection/edit
emphasis, and the marginal cost of the next cut.  Apply quality floors to
prominent objects and use a priority frontier.  Perspective distance matters
only through projection.

## Camera, view, and presentation

### Autoview has one owner and one certified extent

Fitting partial box unions caused visible center/size jumps and cold/warm
camera differences.  BREP leaf boxes can conservatively overestimate trimmed
semantic bounds.

Rule: the progressive source owns one path-scoped exact extent and may apply at
most one framing transition after its priority coverage record is published.
Conservative overview bounds may render immediately but do not certify final
autoview.  Compare final camera center and size across cold/warm and both
renderers.

A tight AABB also exposed an older lower-level assumption: sizing from only
its largest axis is not rotation invariant.  The display-list system had
usually hidden that defect by contributing sphere-expanded object bounds.
`bv_autoview_bounds` now fits the AABB's bounding sphere and accounts for the
viewport aspect ratio.  This keeps a completed autoview valid across ordinary
rotation and prevents resize-dependent clipping without teaching Obol a
second framing rule.  Image or filled-area tests must encode this complete-
target framing rather than relying on the formerly clipped magnification.

### Viewport and device pixels are distinct

Selection rectangles and framebuffer tests failed at fractional device scale
because logical Qt coordinates were compared with physical render pixels.

Rule: every input/capture path records logical canvas size, device-pixel ratio,
and physical viewport.  Test resize, maximize/fullscreen/restore, minimize,
quad views, and HiDPI on actual widgets.

### Stable repaint is not a semantic render request

Always rendering on Qt exposure hid missing invalidation edges and made static
views consume CPU.  Conversely, retaining the framebuffer without converting
a semantic delta into a presentation request left erase/redraw visually stale.

Rule: expose/repaint blits the matching completed image/FBO.  Semantic,
selection, edit, camera, policy, or result publication explicitly requests a
new frame.  Camera-only refresh with no camera change is passive.

### Every new host-work level needs a retained witness

The Qt idle-tail pass could remove the last progress HUD node and thereby
create a fresh presentation request.  OSMesa repainted synchronously, but the
System GL path posted `update()` and retired its timer.  When Qt coalesced that
update, the render latch and a refinement barrier waiting for the next frame
remained pending forever despite empty worker and result queues.

Rule: after any host callback, timer, paint completion, or HUD synchronization
which can create work, resample the level-triggered host snapshot.  If work is
still present, retain exactly one queued witness until a frame or pump consumes
it.  Edge callbacks are wake hints, not the liveness proof.  Exercise this race
with repeated warm System GL runs as well as the focused host-work model.

## Cache and cold startup

### Cache inspection must be non-destructive

The historical bare cache command erased and regenerated every BoT on the UI
thread.  At 150k objects it looked hung and destroyed valuable warm data.

Rule: status is one read-only indexed snapshot.  Deep validation, bounded
prewarm, database clear, and all-file clear are distinct explicit operations.
Whole-cache clear uses one write transaction; targeted invalidation removes the
name publication without deleting shared content-addressed payloads.

### Warm is a certified state

A nonempty cache is not a warm model.  Partial prior runs caused “warm” tests
to enter cold work and made performance comparisons meaningless.

Rule: a warm test requires a completion marker tied to canonical database path
and file identity plus a complete manifest contract.  Relative and absolute
spellings of the same file share a namespace; a copied file deliberately does
not.

### Cold preview must escape long preparation

Lucy produced a valid minimum prefix after classification but kept the box
visible while spatial pages and persistence continued for many seconds.

Rule: publish an immutable, bounded, view-requested preview as soon as it is
valid.  Do not wait for all private metadata or cache writes.  The final result
extends or replaces it atomically; cancellation and transient bytes remain
bounded.

### Concurrency is also a memory policy

Starting every large mesh on every CPU can exhaust memory.  Pure FIFO can also
starve a large root behind a stream of small tasks.

Rule: outer source work and inner mesh work have separate bounded pools and
working-set reservations.  First-fit admission may bypass a blocked item only
a bounded number of times.  One task larger than the limit may run alone.

## Renderer and lighting

### Renderer equivalence requires explicit state

Coin/GL state leaks, stale light slots, color-material state, and swap behavior
produced flashes, brightness differences, and corrupted-looking frames which
did not appear in OSMesa.

Rule: each backend installs and restores a complete explicit render snapshot.
System GL and OSMesa consume the same material, ambient, light, culling, cut,
and instance semantics.  Deep GL sentinels and image comparisons guard both
paths; apitrace is used when System GL alone is wrong.

### Lighting is policy, not backend accident

Straight-on lighting washed out face-on Lucy while old MGED expectations used
a different headlight model.

Rule: one controller-owned lighting profile feeds all backends.  Studio is the
default and historical MGED remains user-selectable.  Tests compare profile
transitions and clear unused lights.

## GUI and client integration

### Selection must remain semantic and sparse

Repeated path resolution in both GED and Qt paint made a 2,522-occurrence
Hubble rectangle selection take seconds and stall the UI.

Rule: libged is the selection authority.  It retains resolved native IDs and
publishes one bulk relation snapshot.  QgModel caches path hashes and performs
O(1) paint-role lookups; sparse changes notify only affected rows.  Only one
globally selected occurrence receives edit-grade Coin geometry.

### GUI controls and commands share one authority

Early primitive and polygon plugins held private copies and could disagree
with CLI state, other views, or retained manipulators.

Rule: GED/librt owns mutable edit state.  Controls and manipulators read after
typed revision events and submit typed operations to the same session.  Plugin
unload has a synchronous barrier which destroys filters, controls, dialogs,
and retained presenters before unloading the DSO.

### Whole-row tree styling is data, not repaint work

Selection/impact highlighting once depended on expensive per-cell path work.

Rule: expose selection roles from cached model state and let a delegate paint
the complete row.  Never resolve geometry or traverse the scene during paint.

## Testing lessons

- Run actual graphical clients.  Headless scene assertions cannot detect GL
  flashes, camera jumps, wrong lighting, mouse-coordinate errors, or retained
  stale pixels.
- Drive real widgets and mouse gestures.  Equivalent commands do not validate
  event filters, coordinate transforms, focus, control readback, or tool state.
- Wait on semantic completion, not fixed sleeps.  Raytraces, LoD, and worker
  drain vary by renderer and cache state.
- Capture every transient frame with APNG when diagnosing flicker; a final PNG
  can hide a one-frame regression.
- A test must reject empty geometry and terminal boxes explicitly.  Earlier
  tests passed because they checked command/state success without verifying the
  presented result.
- Rebuild every affected executable after shared-layout changes.  A stale test
  binary can mimic a runtime defect.
- Do not diagnose production ownership with a test that loads duplicate static
  and shared copies of the same libraries.  White-box assertions belong in the
  owning library; end-to-end tests use the production shared stack.
- Perf captures must exclude or account for diagnostic JSON, path dumps, and
  PNG/APNG compression, which can contribute material harness-only cost.
