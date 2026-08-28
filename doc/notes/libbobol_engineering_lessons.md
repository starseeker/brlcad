# libBObol engineering lessons

Last reviewed: 2026-08-27

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

### A const shared pointer is not an immutable snapshot

The first admitted-geometry API accepted a `shared_ptr<const PartGeometry>`,
but producers could retain a `shared_ptr<PartGeometry>` to the same allocation.
Validation therefore proved only one instant: a later producer mutation could
invalidate indices, bounds, or progressive intervals while render and picking
threads still consumed the object.  Derived wire geometry repeated the alias
through a shared mutable triangle mesh.

Rule: mutable geometry exists only as `PartGeometryBuilder`.  Admission owns a
copy or move, validates it, and privately constructs a const-member
`PartGeometry`; clients cannot construct or copy that renderer-visible type.
Derived representations retain an admitted parent snapshot, never a mutable
subobject.  Regress both the type traits and lvalue-builder independence, and
search the complete consumer tree for direct construction after API changes.

### Process-wide workers impose destruction ordering

ASan found a source-realization worker reading the draw-cache binding map after
the map's function-local static destructor ran.  The coordinator itself was
still waiting for process-exit destruction.

Rule: all lazy non-trivial globals reachable by process-lifetime workers are
constructed before those workers start.  Reverse destruction must join workers
first.  Regress with an active cache callback across return from `main`, and
test the normal shared-library stack rather than a duplicate static harness.

### External node ABI changes require consumer rebuilds

Replacing `libObol` after adding an `SoCADAssembly` field while retaining old
BRL-CAD test objects made C++ allocate the old class size and the new
constructor write past it.  The resulting double-free looked like renderer
ownership corruption; Valgrind identified the first invalid write as the
constructor boundary.

Rule: after an external Obol header or public node layout changes, install both
header and library and rebuild every allocating consumer before interpreting
runtime failures.  Release qualification starts with a dependency-complete
build.  A shared-library timestamp update alone is not evidence that Ninja
rebuilt targets which were not requested.

## State and liveness

### Planning state must not mutate the live presentation

Compact publication originally updated BRL-CAD's retained part/reference and
layer bookkeeping while it was still constructing an Obol batch.  A late
validation rejection therefore left the renderer unchanged but poisoned the
producer's idea of what had been committed; later passes could skip required
records or retire the wrong parts.

Rule: build sparse presentation deltas in a bounded copy-on-write overlay.
Validate the complete geometry/occurrence/style/cut/removal journal against
the live assembly, commit it under one RAII update window, and only then adopt
the producer bookkeeping.  The overlay stores changed keys and signed
reference deltas, never a second copy of a 50k/150k scene.

### Similar proxy properties are not the same contract

A strict validator briefly required every structural proxy to be eligible for
subpixel point aggregation.  Whole-target startup overviews are structural
wire extents which must remain visible in shaded mode, but deliberately must
not collapse to one point.  The coupling rejected valid wire drawing and made
a ready scene empty.

Rule: `structuralProxy` distinguishes temporary extent geometry from authored
wire; `subpixelProxyEligible` independently authorizes point replacement.  A
subpixel-eligible part must have finite conservative proxy corners.  A
structural-only overview is valid and remains unaggregated.  Exercise both
combinations in Obol validation and the complete GED drawing corpus.

### Unity builds can hide source-boundary defects

Renderer extraction initially left one translation unit using a helper which
was visible only because the build concatenated implementation files.  A
nominal file split therefore did not establish a real dependency boundary.

Rule: responsibility-specific renderer and publication units compile outside
unity builds.  Shared private helpers live in an explicit private header or a
single owning source, and each extracted unit must compile independently.

### Transaction control is not swappable scene payload

A staged `SoCADAssembly::replaceScene()` swapped the retained database as one
no-throw commit, but the database base also contained the live update-depth
counter.  The swap installed the candidate's zero counter while an RAII update
scope was open.  Scope destruction therefore observed no open update and sent
no scene notification: queries saw the replacement while render traversal
could remain unaware of it.

Rule: transaction nesting, publication latches, and producer capacity hints
belong to the live assembly, never to a swappable scene snapshot.  A successful
nonempty replacement emits one completion notification; a rejected replacement
emits none.  Preserve both behaviors with a priority-zero node-sensor test.

### Workload regimes are data, not control modes

Lucy, multi-Lucy, Hubble, Generic Twin, and 150k repeatedly exposed opposite
ends of the same scheduling contract.  Fixing one with a workload-shaped latch
or alternate convergence phase made another oscillate, stall, or discard
useful detail.  The problem was not the number of profiles; it was allowing
cardinality and shape to select different control ownership.

Rule: every profile feeds the same inventory, availability, demand, admission,
and presentation pipeline.  Cardinality, projected importance, reuse, bytes,
and measured frame cost are numeric inputs.  They may alter an allocation but
must not introduce another event kind, phase, or completion rule.  The bounded
pipeline TLA+ model varies workload profiles while preserving one transition
relation; graphical matrices validate the numeric and visual consequences the
model intentionally omits.

### Capacity evidence is not a plan cursor

The historical budget object accumulated frame measurements, quality-floor
proofs, retry latches, requested mode changes, and mutable remaining budgets
for a resumable scene scan.  Grouping related fields reduced drift, but left
the controller able to edit the supposed evidence from many unrelated paths.

Rule: completed-frame capacity evidence is immutable planner input.  A plan is
revision stamped and deterministic.  Its bounded execution cursor owns only
mechanical progress and remaining allowance; it cannot change policy or
survive a superseding revision.  Typed presentation outcomes are the only
normal source of new capacity evidence.  Compile-time encapsulation and
transition tests must enforce this split rather than relying on naming.

### Inactive inputs are not planner identity

A retained-allocation cache compared a raw maximum protected budget even when
the corresponding feature was disabled.  Changing that dormant value created
a different transaction key and permitted a numerically identical plan to run
again.

Rule: define planner identity once as a canonical semantic key.  Normalize
guarded-off inputs out of it, and use the same key for cache reuse, in-flight
matching, and result certificates.  Diagnostic fields and inactive policy
values are never evidence.

### Finite work can still be starved by repeated policy waves

A 150k occurrence cursor was individually bounded, and traced 2,048-entry
windows could scan the entire scene in under a second.  The view still took
tens of seconds because completed frames repeatedly reopened capacity
calibration and scene-wide reallocation.  Optimizing the cursor would not have
fixed the controlling cost.

Rule: capacity calibration is a finite search for one frozen population.  Each
event consumes a bounded sample or strictly narrows a safe/unsafe allocation
bracket, and terminal evidence cannot reopen without a semantic revision.
Throughput smoothing proposes candidates; it does not own retry or revision
edges.  Trace both mechanical progress and the number of policy generations.

### A deadline chooses one successor owner

An interrupted frame could both reset an unfinished population cursor and
advance capacity recovery.  On very large scenes the finite scan therefore
restarted indefinitely even though every individual window made progress.

Rule: strict same-transaction progress wins, otherwise a live population
cursor continues unchanged, and only a quiescent producer may yield to
capacity recovery.  Encode this as one sum-valued decision and exhaustively
test its input matrix; do not compose independent retry booleans.

The abstract deadline model already cleared render ownership when population
won, but the first C++ mapping selected the correct enum without retiring an
older replay latch.  A full 50k GUI run found the resulting pending-without-
witness state.  Formal safety of the abstraction is therefore necessary but
not sufficient: every modeled action needs an executable refinement test which
asserts all production effects, including retirement of the previous owner.

### Diagnostic phase machines become competing authorities

The coordinator maintained an imperative event-driven phase machine alongside
the convergence policy which actually derived readiness from work witnesses.
The machine did not control admission or rendering, but every runtime path had
to reconcile it and public diagnostics exposed both answers.  Roughly twenty
notification sites, nine diagnostic fields, and hundreds of lines of tests
therefore preserved a second answer which could drift without changing the
frame.

Rule: user-visible phase and outcome are pure projections of the canonical
inventory, availability, demand, transaction, and work-witness snapshot.
Diagnostics observe that projection; they never own a phase.  A proposed new
mutable control object must identify the unique fact it owns and delete any
predecessor owner in the same change.

### Publication and acknowledgement are one transaction

Applied-result batching and the required-frame barrier previously had separate
reset, pending, deadline, and completion paths.  A completed CAD frame could
clear one half independently of the other, making stale revisions and
pending-without-witness states difficult to exclude.

Rule: one allocation-free, revision-bound presentation transaction owns the
applied-result batch, deadline/request status, reason set, required render
serial, and atomic retirement.  A frame from another view or policy revision
invalidates the transaction instead of acknowledging it.  Unit tests cover
shared reasons, early and stale frames, serial saturation, unbarriered
publication, and source-failure liveness.

### Work levels need progress witnesses

Toolkit notifications are edges and may coalesce; controller obligations are
levels.  Clearing an edge without retaining a timer, scheduled owner-thread
loop, in-flight task, cursor, or frame can strand work.  Treating every repaint
as work hides this defect and wastes frames.

Rule: every nonterminal state names its progress witness.  Pending render
requests are idempotent levels; a capacity request may upgrade a pending
presentation request once.  The host claims the transaction before traversal,
so work published during traversal is a distinct successor rather than a
duplicate reason which continually advances the serial.  The focused
`ObolHostWork.tla` model and window-host tests guard this contract.

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

### Admission acknowledgement is not capacity calibration

The 50k structural scene completed thousands of exact point/box frames but
submitted no meshes.  The classifier frame correctly contained no managed
mesh cost and was routed around capacity calibration; its admission latch was
incorrectly retired only inside that skipped capacity path.  This made the
absence of meshes prevent the transition whose job was to admit them.

Rule: acknowledge a presentation transaction according to what the exact frame
proved, independently of whether its work is a reusable throughput sample.
Changing ownership of an unchanged effective classifier threshold is not a
framebuffer mutation.  Tests must prove nonzero first-wave submission after
settled structural coverage, not merely that the renderer keeps presenting.

### A removed recovery ceiling does not erase its capacity proof

OSMesa wire correctly stopped a Lucy return view after the next richer static
population missed the hard presentation deadline, but reconciliation removed
the temporary renderer ceiling before convergence status was queried.  The HUD
therefore called a 1.29-pixel result unconstrained against a one-pixel target.

Rule: the rejected static-quality trial is a view/policy/source-epoch capacity
witness, not an implementation detail of the temporary ceiling.  Preserve and
report that typed performance constraint until a genuine camera, policy,
resource, cadence, or population edge invalidates it.  Qualification compares
the current exact frame with its current budget and protected-quality result;
raw work from an earlier view is neither a fidelity ordering nor a valid
capacity predicate after view-local clustering and recalibration.

### A first static frame is not necessarily steady-state capacity evidence

Lucy could complete an approximately 95 ms static frame after a first attempt
had exceeded the 100 ms endpoint deadline while newly published resources were
being installed.  The preparation serial did not capture that renderer-side
first-use cost, so the controller permanently rejected static quality and
forced a much coarser 80 ms marginal allocation even though the stable frame
was valid.

Rule: an explicit static-quality trial gets one unchanged first-use retry.  A
later retry is justified only by a revision-bound preparation obligation whose
finite remaining work strictly decreased.  A changed diagnostic serial is not
that proof: uploads may evict one another or advance work for a superseded
target.  A no-progress interruption is typed capacity evidence.  This rule is
workload and renderer independent and must not become an open-ended retry path.

### Resident demand and presented quality are independent

A warm single-mesh Lucy view had its complete requested PoP suffix resident,
so every provider request was classified as satisfied.  After pose-only
rotation, however, the renderer still submitted through a much coarser global
ceiling learned at the ordinary quiet cadence.  The convergence controller
mistook provider completion for visual completion and never entered the
separate static-quality allowance; timing then decided whether the final image
showed roughly 500k or 2.5M faces.

Rule: quality debt is the union of unsatisfied resident demand and a
presentation ceiling below the active progressive population.  The latter is
actionable even for exactly one progressive occurrence.  Use the reversible
renderer ceiling during input, then evaluate the hidden resident population
under the bounded static-frame contract.  Do not reload or rewrite occurrence
cuts merely to restore data which is already resident.

### Hard aborts cannot use a completed-frame deadband

The first warm 50k OSMesa wire replay repeatedly aborted at about 101 ms
against a 100 ms endpoint deadline.  The small-part controller reused a five-
percent completed-frame jitter deadband, returned the unchanged one-pixel cut,
and requested the identical incomplete frame.  Once that was corrected, a
second legacy path walked the point cut toward safety and then reopened the
same already rejected one-pixel static trial, producing a larger
one-pixel/coarse cycle.

Rule: a hard-aborted irreducible frame is typed endpoint evidence, not a noisy
completed timing sample.  It must advance the remaining reversible population
control by a bounded step.  A rejected static-quality population remains
rejected for its view/capacity epoch; changing the recovery point threshold
does not make the rejected one-pixel population cheaper and cannot invalidate
that proof.  Camera, user-policy, renderer/resource, or cadence epoch changes
may start a new trial.  Regression qualification must include a truly fresh
cold/first-warm pair because later cache and scheduling state can hide this
class of liveness failure.

### Structural points are not hidden richer geometry

The 150k OSMesa shaded replay exposed a second coarse/fine cycle after the
hard-abort fix.  Triangle recovery and static-quality headroom paths reset the
aggregate point threshold to one.  That is valid only when every point stands
in for an already realized mesh.  During streaming population many points are
the best published structural fallback; lowering their threshold reveals
boxes, so the purported quality recovery is actually a visual regression.

Rule: triangle-prefix recovery and point/box aggregation are independent.
Keep the last proven aggregate cut while any structural fallback occurrences
remain.  A one-pixel static-quality trial is allowed only for a fully realized
population.  Mesh publication, rather than a triangle-recovery endpoint,
removes the structural constraint.  Verify this with frame-sequence evidence:
one bounded capacity backoff is legitimate, but a return to a previously
rejected denser population in the same epoch is a control-cycle failure.

### PoP spatial clusters are not scene occurrences

The one-Lucy mesh source carries many spatial PoP clusters so that zoomed views can
load only needed pages.  A point-aggregation guard used that transient cluster
census as if it were a scene occurrence count, which let exact Lucy frames
collapse into points despite there being only one logical object to draw.

Rule: retain a source-wide logical occurrence count separately from the
view-local spatial census.  A count of one vetoes point aggregation.  For
multi-occurrence sources, apply the visible census as well, so an assembly with
only one visible part is not needlessly reduced to a point.  Page/cluster
partitioning is an implementation detail; occurrence identity is a scene
semantic.

### Ready is a semantic proof, not an empty queue

An empty worker queue may still have unpublished results, unacknowledged frames,
partial visibility census, pending resident allocation, or source preparation.
Conversely, a constrained scene can be validly ready without all original
geometry resident.

Rule: readiness combines exact live-source coverage, visible occurrence
classification, affordable cut satisfaction, terminal failure/constraint
reasons, drained publication, and acknowledgement of the current frame.

### A completed frame is capacity evidence even if new work is queued

Lucy could draw a richer in-gesture cut inside the stable 100 ms deadline, but
the next queued submission changed mutable scene status before quiet handoff
consulted it.  The controller then forgot that proof and replaced the displayed
cut with an unnecessarily coarser one.  A second defect rounded a 0.75-pixel
request to one pixel in the scene error-ceiling helper, even though the
provider and retained allocator honored the fractional request.

Rule: associate deadline capacity with the immutable completed-work record,
not later mutable scene status.  The richest completed cut is reusable across
policy-only transitions for the same camera and source-quality domain; camera,
renderer, or source-population changes invalidate it.  Preserve the exact
fractional pixel target at every selection and allocation boundary.  The hard
deadline still clamps an attempted richer replacement below the cut that
actually missed.

### First-useful-image timing starts after the canvas is presentable

An OSMesa Hubble row issued `draw` while the Qt software canvas was still
waiting for its first paint.  The controller later reached an exact zero-box
frame, but the test attributed canvas startup to source realization.

Rule: GUI qualification waits for the active canvas's explicit first-paint
state before measuring draw-to-first-useful-image latency.  This is a
state-based barrier, not a startup sleep; it preserves cold-stream timing once
a user can actually see the canvas.

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

### Resident growth and allocation are one transaction

A compact source can finish database inventory while bounded cache/service
waves are still appending immutable PoP suffixes.  Treating that inventory
edge as a settled scene ran the global minimax allocator in every momentary
queue gap; each late result invalidated its certificate, rebalanced all visible
cuts, and made 150k views alternate coarse and fine for tens of seconds.

Rule: keep the last coherent framebuffer while provider, task, result,
publication, or resident-drain work owns the population.  Run one scene-wide
importance allocation only after all of those owners are quiet.  A hard failed
headroom trial is a negative capacity witness, not a cancellation: remember it
until the view/policy/population changes or materially greater capacity is
proved.  Failures observed before the population closes may constrain that
attempt, but must not poison the final event-driven static-quality pass.

### Physical quality debt is not resident work

Loading every pixel-demanded suffix helped a single large Lucy mesh refine,
but made a performance-limited 150k scene repeatedly invalidate its closed
population: each cache-result wave populated data which the current allocation
could not draw and triggered another global allocation.  Suppressing all such
loads instead stranded large meshes during active zoom.

Rule: a quiet view loads only through the allocator-selected presentation cut.
Keep finer physical demand as explicit quality debt so a later view or capacity
revision can reconsider it.  Active scale interaction may prefetch one bounded
transition past a stale allocation, subject to the normal worker and resident-
memory limits; presentation still stops at the allocated cut.  Model physical
debt separately from admitted resident work, and never report the former alone
as a pending progress owner.

### Resumable rendering requires a closed scene population

A 150k OSMesa frame often spends its deadline preparing the retained point,
instance, atlas, and command record rather than rasterizing it.  The renderer
correctly preserved its cursor, but the ordinary owner-thread pump could apply
another provider or result wave before the successor frame.  Each nominally
resumable slice then described a different scene and discarded the work just
completed, amplifying a bounded preparation into long, schedule-sensitive
convergence.

Rule: after an interrupted CAD traversal, freeze every owner-thread scene
publication path until the exact successor presentation completes.  Keep
workers behind their existing byte/result bounds; do not block or discard
their immutable output.  Release the publication gate on every completed CAD
frame, including presentation-only frames, before capacity-specific early
returns.  Test both halves: a provider must not run while replay is pending and
must run immediately after the completed-frame commit.

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

### Shared retained content needs mutation-aware presentation fanout

GED-wide annotations, features, and polygons live under a shared retained
root, while each graphical view owns a different framebuffer and presentation
callback.  Updating the shared store correctly requested a render from its
controller, but that controller has no graphical host.  System GL therefore
kept an old framebuffer even though an observational OSMesa traversal could
appear correct.  Conversely, waking every view after every successful `view`
command made read-only inspection needlessly redraw large scenes.

Rule: a successful shared-store mutation advances the shared controller's
request serial.  Compare that serial around command execution and fan out one
presentation-only request to graphical views only when it changed.  Store
updates remain renderer-neutral; queries remain passive.

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
valid.  Its bound is the exact per-task share already reserved by the scene
allocator: a minimum PoP prefix is not a free exception.  A cold scene of
distinct large meshes must retain structural coverage for a source whose
indivisible prefix exceeds that share, rather than publishing one costly
preview per worker before a completed frame can calibrate the aggregate.  Do
not wait for all private metadata or cache writes.  The final result extends
or replaces it atomically; cancellation and transient bytes remain bounded.

### Concurrency is also a memory policy

Starting every large mesh on every CPU can exhaust memory.  Pure FIFO can also
starve a large root behind a stream of small tasks.

Rule: outer source work and inner mesh work have separate bounded pools and
working-set reservations.  First-fit admission may bypass a blocked item only
a bounded number of times.  One task larger than the limit may run alone.

### Source size is not screen importance

A tempting cold-start shortcut is to render every small authored BoT directly.
That fails in a large view: thousands of individually tiny meshes can all be
subpixel, and direct admission bypasses the aggregate channel precisely when
it is most useful.

Rule: discovery may classify a whole source population as a conservative
safe-zone candidate using count, source bytes, and reuse facts, but only the
view allocator may choose a representation.  It must project each occurrence
first, aggregate subpixel coverage, and apply the same protected importance
and draw/memory budgets to direct and PoP candidates.  A source profile is a
gate, never a per-object entitlement.

### Semantic population counts are not presentation counts

Generic Twin's complete safe-scene profile contained 2,248 leaf occurrences,
while its compact presentation temporarily contained 2,249 records because
the synthetic whole-target overview is also a renderer record.  Requiring
those counts to be equal permanently disabled direct-terminal admission.  A
second hidden assumption treated all unique primitive leaves as unique mesh
assets and rejected the same small 4.1 MB scene at an arbitrary 2,048 cutoff.

Rule: a completed source profile authenticates its own semantic population.
Do not compare it with a renderer/presentation collection that is permitted to
contain auxiliary records.  Profile limits may conservatively gate a scene,
but exact projection, type, draw-cost, working-set, and resident-memory checks
must remain the per-mesh admission authority.  Graphical tests must report and
assert the selected representation kind; eventual non-box geometry is not
proof that a direct route ran.

## Renderer and lighting

### Compatibility-state allocation is transactional

Under a capped 150k OSMesa replay, `_mesa_PushAttrib` dereferenced a failed
allocation while copying attribute state.  Merely catching later renderer
buffer allocation could not protect this earlier compatibility-stack boundary.

Rule: OSMesa constructs server and client attribute snapshots completely
before publishing a stack entry, cleans nested references on failure, and
reports `GL_OUT_OF_MEMORY`.  Obol verifies that each requested stack depth
actually advanced, pops only successful pushes, and rejects the candidate CAD
frame before scene mutation when the guard is incomplete.  Resource-pressure
tests must cap the address space and exercise both the provider and consumer
side of an allocation boundary.

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

### Tool refresh is not a selection gesture

The polygon editor once treated the combo box's first item as both a display
default and a scene selection.  A paint, checkpoint, or command readback could
therefore select a polygon merely by refreshing the widget, altering command
state and selection appearance without user input.

Rule: the retained selection record is authoritative.  Passive widget refresh
only mirrors the selected record; mouse selection, an explicit combo change,
and a creation/deletion gesture are the only transitions allowed to change
selection.  Creation may establish an initial target only if none exists, and
editor deletion may deliberately advance to a remaining target.

### Whole-row tree styling is data, not repaint work

Selection/impact highlighting once depended on expensive per-cell path work.

Rule: expose selection roles from cached model state and let a delegate paint
the complete row.  Never resolve geometry or traverse the scene during paint.

## Testing lessons

### Classify a counterexample before changing control state

The same coarse Lucy image arose from three different causes during this
work: incomplete renderer preparation, an over-conservative numeric deadline
correction, and overlapping quiet-restoration policy owners.  Treating all
three as "LoD did not refine" invites another retry latch or model-specific
threshold and makes the controller less reviewable.

Rule: first classify an observed failure as a contract, numeric policy,
data-plane, or observation defect.  Contract defects change a finite event or
typed certificate and its formal model.  Numeric defects change pure functions
and independent oracles.  Data-plane defects stay below the controller and
gain boundary/sanitizer/image tests.  Observation defects cannot mutate a LoD
revision.  If a proposed fix spans classes, split it; if the class is unknown,
pause implementation and update the ownership design.

### Deadline ratios must not compound unnamed safety penalties

An exact OSMesa frame completed a fraction of a millisecond beyond a 400 ms
deadline.  Scaling its attempted cost by `deadline / elapsed` was already the
required overload correction; multiplying that result by another unexplained
0.80 collapsed the next stable scene by more than twenty percent and made an
exact zoom return visibly coarse.

Rule: capacity arithmetic is a named pure policy with unit tests at the
boundary.  A deadline ratio and any headroom margin have separate, documented
purposes.  Do not stack arbitrary constants in controller effect code.  Timing
traces are evidence for the pure estimator, not permission to add state.

### A bounded renderer obligation starts before an expired deadline poll

A later CAD assembly may be reached after an earlier scene node has exhausted
the host frame deadline.  If it checks that deadline before publishing and
advancing any finite work, the host sees no progress and either starves the
assembly or mistakes preparation latency for permanent draw capacity.

Rule: publish the exact preparation target first, grant one bounded planning
quantum, and then honor interruption.  Subsequent retries require strict
same-target remaining-work reduction.  This is fair admission of an existing
obligation, not an exception to the frame deadline.

### A broader predicate is not necessarily a common contract

Replacing a single-visible-occurrence condition with a superficially general
"restored certified presentation" predicate improved Lucy, but made 50k warm
startup and 150k post-rotation restoration exceed three-minute test bounds.
The predicate combined two different facts: a target was safe to restore, and
new static-quality work was affordable.  Object count had only hidden that
ownership error.

Rule: generalize behavior by reducing it to distinct typed events, not by
widening an existing condition.  A restored presentation is completed
evidence.  A successor refinement is a new revision-stamped plan requiring
current demand and capacity evidence.  If a change helps one qualification
profile but harms another, revert it and update the common refinement map;
never retain it behind workload-shaped conditions.

### One completion event may select only one control successor

A completed framebuffer simultaneously satisfied a retained presentation
barrier and a capacity sample.  Both effect paths reacted: the capacity owner
requested its next sample while generic presentation cleanup restarted a full
occurrence pass.  At 150k scale the legal individual transitions composed into
an expensive refinement/balancing cycle.

Rule: successor ownership is exclusive.  A stronger typed transaction owns its
whole sequence; weaker observers of the same completion may retire bookkeeping
but may not start work.  Model and test the effect ordering, not merely the set
of valid snapshot owners.

A later multi-Lucy counterexample showed that the implementation still chose
capacity first and asked the handoff whether it superseded that choice only
afterward.  Capacity calibration had already changed its evidence and installed
a retained-allocation request, so deactivating the visible submission did not
cancel the hidden successor.  `completedPassSelection` now chooses one typed
owner from the immutable completion snapshot before either controller runs.
This is the production seam checked by `ObolCompletedPassOwnership.tla`.
Starting the chosen successor also resets the complete pass-annotation
transaction; a stale `changedCut` bit cannot make the next no-op pass look
productive.  A legacy epilogue which advanced the capacity revision after a
no-op terminal pass was removed for the same reason: revisions name semantic
edges, not pump calls.

A subsequent contract audit found that the ceiling-free capacity-sample
exception was still normalized only after capacity calibration.  The later
handoff could therefore run after the supposedly exclusive capacity owner.
Normalization now belongs to the immutable selection, structural/point effects
are gated by that selected owner, and the unselected handoff reducer cannot run.
The multi-pass model covers the successor effect and annotation lifetime, not
just the valid owner set at one snapshot.

### A multi-step availability transaction belongs in one sum type

Resident growth was represented by ledger debt plus a controller drain boolean.
The drain could finish while the debt still correctly prohibited allocation,
and generic budget handling restarted the drain indefinitely.

The first typed repair still omitted one guard from its formal abstraction.
After an OSMesa Generic Twin zoom, useful coverage was incomplete and resident
growth was required.  Coverage yielded to growth, but the growth scheduler
required complete coverage before beginning its drain.  It consequently
replayed the same incomplete coverage pass thousands of times with no task,
result, or framebuffer capable of changing either premise.

Rule: required, active, dirty-active, and allocation-ready are mutually
exclusive phases owned by one ledger.  Coverage is an output of the resident
drain, never its admission guard.  Drain completion has only three successors:
another drain for a coalesced arrival, one scene-wide allocation when coverage
is complete, or one ordinary coverage retry after consuming the growth edge
when coverage remains incomplete.  A new view interrupts the active drain back
to required; it never silently consumes old-view work as a current allocation
proof.  The focused model must include every runtime guard on those successors;
modeling phases without their coverage precondition did not prove this boundary.

### Presented evidence needs a level-triggered retirement path

Point-to-triangle recovery once depended on handling its exact completed frame
in the same callback which retired it.  When a stronger capacity or handoff
barrier owned that callback, recovery remained pending after every task, cursor,
and render request had disappeared.

Rule: preserve the fact that the recovery was presented and let the ordinary
progressive pump retire it idempotently when stronger barriers clear.  A typed
obligation must always name a producer, a pending frame, or a level-triggered
completion witness.

### Coalesced asset work and view demand have different lifetimes

A cold Lucy hierarchy continued correctly across camera changes, but its live
pages and final result carried the request captured when construction began.
The controller rejected those results, leaving only coverage.  Updating the
provider stamp alone exposed the inverse error: the generic service validator
compared the current result with the launch request and converted it into a
terminal `STALE` result.  Old exact-demand failures then survived later policy
epochs and could falsely terminate an otherwise healthy warm view.

Rule: immutable asset production owns a stable asset key; page selection and
result routing own the latest view/policy demand.  The producer must refresh
demand immediately before view-dependent selection, and the service must
validate against that same retained demand.  A later camera change may reject
the completion but may not turn it into a current-demand provider failure.
Failure suppression is exact-demand state and is discarded when its view or
policy epoch advances.  Test the full sequence—launch, multiple demand changes,
intermediate publication, final completion, normalization, owner-thread apply,
and a subsequent warm refinement—not only the provider callback.

Do not overload source validity with scheduling currency.  A stale cache or
changed source metric is a provider failure for an exact demand; a valid result
overtaken by newer demand is `SUPERSEDED` and is simply discarded.  Using one
status for both made normal camera timing capable of terminating convergence.

### Hidden tombstones do not own a stable identity

Sparse box-to-mesh publication retains an old hidden plan slot and appends the
active mesh slot.  Both records deliberately carry the same stable
`InstanceId`.  A later PoP cut-bin swap in the old group used to rewrite the
global ID-to-index map unconditionally, so the hidden tombstone stole the
active replacement's index.  The next cut update combined that old index with
the new group's metadata, underflowed two zero bin counts to `UINT32_MAX`, and
asked the renderer to scan more than eight billion phantom occurrences.

Rule: a duplicate retained record may update an identity index only when the
mapping still names the record's source slot.  Validate group cardinality,
span membership, and source-bin ownership before the first sparse mutation;
fall back to an authoritative rebuild on mismatch.  Test the complete shared
part, replacement tombstone, old-group swap, active-replacement cut sequence.
Checking only final instance count or one ordinary rebind does not exercise the
aliasing failure.

### Cancellation is not invalid geometry

Flat atlas construction formerly returned one Boolean for malformed input and
deadline interruption.  An abort after one committed range therefore changed
a valid `Preparing` certificate into permanent `Failed`.  Separately, one
source-authored PoP range could exceed the entire deadline, so retries never
committed even their first unit.

Rule: preparation outcomes are typed.  Split source ranges into bounded atomic
units, retain complete units across frames, report interruption without
changing validity, and reserve `Failed` for data which cannot become valid on
retry.  A regression must interrupt inside a source range, observe strict
partial progress, then resume the identical target to completion.

### Thread cleanup must retain registry entries through callbacks

The source-realization shutdown contract test exposed two independent races on
2026-08-26.  First, a job published its terminal state before retiring its final
active-work and working-set accounting.  Under concurrent execution this made a
terminal constrained request coexist with one stale active item.  Terminal is
now the source coordinator's commit point: active accounting is retired before
the release-store which publishes completion.

A rarer process-exit crash remained.  Concurrent stress observed three raw
segfaults in 800 pre-fix executions and one in 1,600 executions after the
coordinator ordering fix.  Two debugger captures produced the same main-thread
stack while a realization worker was deliberately still inside its shutdown
cache callback:

```text
___pthread_mutex_lock(invalid cc_storage mutex)
cc_mutex_lock
CoinInternal::StorageRegistry::cleanupThread
CoinInternal::ThreadCleanupTrigger::~ThreadCleanupTrigger
__call_tls_dtors
__run_exit_handlers
```

At the captured failure, the worker threads and coordinator were still alive;
the main thread was running Obol's thread-local cleanup.  Obol's
`StorageRegistry::cleanupThread()` had snapshotted raw `cc_storage *` entries,
released the registry lock, and then invoked arbitrary storage destructors.
Such a callback can unregister and destroy another snapshotted storage before
the later raw pointer is visited.  Snapshotting prevents a registry-lock
deadlock but does not provide object lifetime.  The durable contract is that a
cleanup traversal must pin every entry it may dereference, or serialize entry
destruction with a reentrant registry lifetime lock.  Exception suppression
cannot make dereferencing a freed entry safe.

The original crash report at
`.build/src/libBObol/tests/test_bobol_source_realization-29-bomb.log` is older
and records a distinct `db_clone_dbi` failure; do not conflate it with the two
identical Obol TLS-cleanup stacks above.  Recheck current upstream Obol before
implementing the storage-lifetime repair, then retain a focused thread-exit
regression in Obol and the active-callback process-exit regression in libBObol.

- Run actual graphical clients.  Headless scene assertions cannot detect GL
  flashes, camera jumps, wrong lighting, mouse-coordinate errors, or retained
  stale pixels.
- Drive real widgets and mouse gestures.  Equivalent commands do not validate
  event filters, coordinate transforms, focus, control readback, or tool state.
- Wait on semantic completion, not fixed sleeps.  Raytraces, LoD, and worker
  drain vary by renderer and cache state.
- Keep ordinary timing observational.  Aggregate streamed-source and resident
  mesh timing at the completed realization; per-batch or per-asset logs belong
  behind an explicit verbose switch because synchronous logging changes the
  schedule and can manufacture a large-model regression.
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
