# libBObol engineering lessons

Last reviewed: 2026-09-03

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

### Cache failure is not geometry rejection

The BREP source path once reported `PoP input rejected` whenever cache setup,
cache writing, hierarchy generation, or input validation returned the same
boolean failure.  In a restricted process the cache directory was simply
unwritable, but the message blamed a valid tessellation and the source then
remained a structural box.  This made an environmental acceleration failure
look like a topology defect.

Rule: validate source arrays separately and report their exact invariant
failure.  Treat cache context, generation, write, and reopen as distinct
operations.  A finite tessellation with in-domain indices remains terminal
drawable geometry when optional PoP/cache acceleration fails.  Tests which
exercise cache behavior provide an isolated writable `BU_DIR_CACHE`; tests of
the fallback deliberately make acceleration unavailable and still require a
terminal mesh.

### Multi-face tessellation success means complete source coverage

The legacy BREP CDT could fail several constrained faces, append triangles
from the faces which happened to succeed, and return aggregate success.  A
Big Boy tire consequently supplied 368 allocated points but referenced only
32 of them on one upper plane.  PoP correctly classified that defective input
and made its missing faces visible at every progressive level.  A conventional
BREP edge drawing beside it was not shaded ground truth and made the failure
look like an LoD disagreement.

Rule: a multi-face provider may publish success only when it proves the
required source coverage.  Validate referenced rather than merely allocated
vertices.  Exact BREP boundary vertices are lower-bound coverage witnesses;
the untrimmed NURBS surface box is only an upper envelope because its control
hull may extend beyond trimmed faces.  Invalidate persistent results whenever
that contract changes.  Post-hoc validation is containment, not scheduling:
the provider itself must bound time, memory, result size, and cancellation and
must report per-face completion with a typed outcome.  LoD generation starts
only after that source contract succeeds.

### Progress reports use the producer's finite work rank

Cold BREP preparation could visibly replace boxes one at a time while the HUD
stayed at a derived percentage or appeared to restart discovery after
provisional autoview.  Representation count is not the work queue: one
completed tessellation may contain many renderer records, while one pending
source item may still dominate the remaining time.

Rule: once source planning knows its finite detail-work inventory, publish
saturating completed/total counters and aggregate those counters without a
scene traversal.  The convergence policy owns the phase and maps that exact
rank into the preparation band; the HUD only formats it.  Before the total is
known, say that preparation is active without a percentage or ETA.  Database
discovery retains phase priority only until its leaf census is complete.

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

### Callback removal is a lifetime transaction

A presentation hook copied its callback and raw endpoint pointer under a
mutex, then invoked them after releasing the lock.  Endpoint destruction could
clear the hook and free its state while an already-dispatched callback still
used that pointer.  The frame hook tracked active dispatches, but removal from
inside its own callback waited for itself and deadlocked; a nested callback
could create the same cycle while removing an enclosing dispatch.

Rule: callback registration owns a dispatch state, not just a function and raw
pointer.  External removal closes registration and waits for every acquired
dispatch.  Removal from any active frame on the current callback stack defers
storage deletion until that frame unwinds.  Test blocking external removal,
self-removal, nested cross-removal, and endpoint destruction during dispatch.

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

### Resident progressive data still needs view-binding identity

The first adaptive BREP wire path correctly retained all cut arrays in one
immutable Obol part, but initially left the richest cut active because only
service-backed mesh payloads entered view planning.  After direct cuts became
view-managed, a second hazard remained: an in-place source edit could replace
the immutable part while aggregate cut and convergence counters still
described the old generation.

Rule: availability and view binding are separate.  A fully resident
progressive part needs no worker/cache payload, but it remains a logical
progressive occurrence with requested/active cuts and render cost.  Authenticate
that binding with routing ID, population epoch, occurrence key, and geometry
revision.  Retire its metrics from the source's sparse mutation journal before
planning the replacement; use an authoritative binding scan only when delta
history is unavailable.  The renderer must receive the authored minimum when
no authenticated view binding exists.

### Terminal evidence is not an in-flight transaction

A bounded renderer-capacity search correctly reached a certified Lucy wire
budget, applied its final occurrence allocation, and then restarted from the
same two populations.  The generic changed-cut completion path used a broad
"reset calibration" operation which erased both mechanical frame state and
the terminal certificate.  The search itself was finite; its caller made the
same completed problem look new.

Rule: distinguish retiring an in-flight measurement from invalidating semantic
capacity evidence.  Applying a certified final cut and completing its generic
pass preserve the terminal certificate and budget-limited witness.  Only an
explicit population, view, policy, or resource-capacity invalidation may erase
that evidence.  Model the post-certificate application path, regress the
unchanged-problem rearm directly, and qualify it with the full interaction
lifecycle rather than only the isolated search.

Retaining the certificate is insufficient if ordinary planning ignores it.
A follow-up Lucy trace showed the search finish at the affordable discrete
population, after which the throughput planner immediately selected the
adjacent known-slow population and caused another identical search.  A current
terminal certificate therefore caps ordinary admission as well as preventing
measurement rearm; only semantic invalidation or a separately bounded
structural-coverage transaction may supersede it.

### Exact capacity search is the wrong user-facing objective

The warm 50k OSMesa lifecycle eventually satisfied the finite liveness model,
but exposed 43 scene-allocation plans and 74 capacity frames.  After the last
input event, the view spent roughly 129 seconds alternating between balancing
and refinement.  The search was terminating correctly; it was repeatedly
solving an unnecessarily exact problem as resident availability changed.

Rule: renderer capacity is a guardrail, not the visual objective.  Preserve a
proven-safe visual, enrich from it, and stop a goal after the first rejected
richer population.  Bound cold recovery before any safe result.  A completed
frame classifies every ordered deadline it meets, so changing from steady to
static policy must not discard a population already proven safe for the static
deadline.  Formal liveness must be paired with a bound on candidate
publication and a graphical wall-time qualification; "eventually" alone says
nothing about acceptable settling latency.

### One typed request must not have two scalar truths

The presentation-reconciliation request correctly retained the tighter of
repeated budgets, but its surrounding evidence copied each caller's raw value
into a separate active-budget field.  A weaker reassertion could therefore
leave one transaction reporting a tight request and executing a looser
allocation.  It also restarted the scene cursor even though no semantic input
had changed.  Lucy exposed the split as a pending handoff whose selected cost
matched neither the requested nor certified budget.

Rule: a typed request owns its dependent scalar values.  Mutators report
whether the canonical request changed, consumers read the canonical value
back from the request, and bounded cursors restart only for a real tightening.
Model duplicate requests explicitly and assert that the active allocator
budget always equals the canonical request budget.

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

### A service quantum is not a quality ceiling

The 50k point-relaxation path priced an affordable 10,849-occurrence finer
population, then rejected it because the structural first-wave service limit
was 8,192 occurrences.  That limit exists to bound one provider admission
wave; treating it as a semantic terminal criterion left the view at a
needlessly coarse 16-pixel aggregate threshold.

Rule: price a quality transition using its complete projected population, but
execute it as deterministic bounded batches while retaining the preceding
coherent frame.  Every exact batch completion must strictly decrease the
remaining-occurrence rank.  Zero rank is the only publication edge; no
progress, an invalidated domain, resource denial, or a missed deadline rejects
the private candidate without changing the public cut.

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

### A candidate allocation is not the captured framebuffer

After a pose-only Lucy change, OSMesa correctly retained resident cut 24 at
0.712 normalized presented error.  The new allocator plan described the
coarser affordable cut 23 at 1.103.  The visual-quality harness preferred that
planning maximum and rejected an exact framebuffer which was demonstrably
richer than the unused candidate.

Rule: framebuffer qualification uses error recomputed from the exact presented
cuts and requires an exact completed-work record.  Keep the allocator's
maximum as separate planning evidence for capacity and convergence audits.
Never use a candidate, requested cut, resident suffix, or stale certificate as
a proxy for pixels which the renderer actually presented.  Conversely, a
coarser exact framebuffer cannot borrow the error of a richer candidate; the
two values remain independently named.

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

### Every aggregate producer must preserve the same extent contract

The ordinary view classifier retained screen-significant aggregate bounds,
but the indirect renderer's atlas-admission fallback constructed a separate
center-only point record.  Large objects therefore became large raster points
only in the highest-pressure regime, where ordinary model tests were least
likely to exercise the defect.

Rule: aggregate shape is a shared presentation contract, not a producer-local
shortcut.  Every normal, append, and resource-pressure producer must carry the
conservative world bounds and classify the complete current projection.  A
point is valid only when that projection is fully classifiable and no larger
than five pixels; clipped or larger projections remain boxes.  Exercise the
pressure path with a deliberately undersized atlas and assert physical point,
line, and triangle work, rather than relying on a nonempty framebuffer.

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

### Constrained must name the evidence which stops refinement

A 50k cold OSMesa endpoint could truthfully have no pending workers, results,
or presentation transaction while retaining visible quality debt.  The old
telemetry reported only `performanceLimited`, forcing a reviewer to infer from
budgets and timing whether the allocator had reached a measured endpoint or
had simply stopped early.  Trace logging also changed the scheduling enough to
select a different affordable marginal population, making that inference even
less reliable.

Rule: a constrained terminal outcome carries a derived, nonempty evidence
mask.  Name stable-frame capacity, renderer ceiling, subpixel aggregation,
static-deadline rejection, memory pressure, or a bounded representation limit
directly.  Never treat visual debt, an empty queue, or a coarse image as the
constraint.  Keep the mask diagnostic and derived from the owning policy
values so it cannot become another state machine.  Qualification rejects a
constrained sample without this witness.

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

### CAD timing identity includes its presentation threshold

A 150k System-GL view repeatedly alternated between 32- and 64-pixel point
classification even though its CAD population was unchanged.  The expensive
CAD traversal at one threshold was followed by a cheap faceplate-only repaint;
the controller paired the cheap host duration with the retained CAD work
record and falsely inferred headroom.  The Qt software host already measured
whether CAD executed, but an omitted positional argument silently defaulted
that fact to true.

Rule: CAD work identity, point threshold, duration, upload/replay classifier,
and execution fact form one sample.  A sample is reusable only at its stamped
threshold.  A non-CAD frame may update user-facing host timing but cannot
create or consume CAD capacity evidence.  Hosts supply capacity relevance and
CAD execution through distinct required enum types, never defaulted Boolean
arguments.  Preserve this boundary with `ObolCadTimingEvidence.tla`, focused
evidence tests, and a high-cardinality terminal replay.

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

The same distinction applies after resident reclamation.  A capacity revision
can make a previously denied suffix admissible while the retained scene
allocation still names the older coarse cut.  Capping the retry request at that
stale allocation consumed the revision without requesting any bytes: the
controller then had three known quality debts, zero tasks, and no future event
capable of improving them.

Rule: a resident-admission retry invalidates the old allocation as a physical
residency cap, not as presentation authority.  Probe the current pixel demand
through the ordinary byte-governed service; if admitted, perform one successor
scene allocation, and if denied, retain an explicit bounded constraint.  A
retry token may never retire as a no-op while its requested and drawable cuts
differ.  `ObolCapacityPresentationHandoff` models this as a required/requested/
resolved-or-constrained transaction.

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

### GUI replay owns an isolated desktop contract

Skipping saved-window restoration did not make qged replay hermetic.  Plugins
could still inherit the operator's enabled state, replay shutdown could write
its final geometry back to the operator's QSettings, and a window manager
could acknowledge an older configure request after `QWidget::resize()` had
already reported the new local size.  Autoview then fitted the same model to a
different aspect ratio in managed and full-detail runs even though both drawing
policies were correct.

Rule: test mode uses a process-private temporary QSettings scope before any
client reads settings, and never persists window state.  A scripted top-level
size or state request carries a monotonically increasing generation.  Delayed
guards act only while that exact generation remains current, and a caller may
request a bounded native-window stability interval before issuing camera
commands.  Regress an injected late configure acknowledgement, run on a real
window manager as well as offscreen, and compare the final camera across LoD
policies rather than trusting the immediate resize return value.

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

The first repair compared the shared controller's render-request serial around
command execution.  That serial identifies a coalesced frame transaction, not
content.  Because the shared controller has no endpoint to consume its latch,
the first mutation advanced the serial and every later mutation joined the
same request.  Selection and movement state were correct, but System GL kept
the framebuffer produced before those later changes.

Rule: each shared store owns a monotonic presentation revision which advances
for every actual retained-content mutation, even while one frame request is
latched.  Compare the store revision around command execution and fan out one
presentation-only request to graphical views only when it changed.  Equal
publication and read-only queries preserve the revision.  Render-request
serials remain frame-transaction identity and must never stand in for content
identity.

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

### A mesh-channel preview is not a resident progressive binding

A cold shared-mesh run could publish useful prepared geometry for one
occurrence, complete the shared PoP hierarchy through another occurrence, and
then loop forever trying to refresh the first.  The prepared preview carried
the adaptive mesh result kind, so the submission action incorrectly sent it
through resident-cut retargeting even though it had no progressive hierarchy.
The retarget could never succeed, and shared-asset reuse was never reached.

Rule: classify lifecycle capability from the payload it actually owns, not
from the broad result channel which transported it.  An adaptive mesh payload
without a valid progressive hierarchy is provisional coverage.  Keep it
visible while replacement is prepared, charge only the marginal replacement
cost, and bind it directly to a compatible resident shared hierarchy as soon
as one exists.  Readiness requires the replacement for every occurrence,
including the occurrence which published the original preview.  Exercise this
with cold and warm shared-mesh GUI runs plus a focused update-action test; the
cold-preview TLA+ model's `resident` binding maps to this capability test, not
to the result-kind enum.

### Whole-prefix residency is not spatial-page residency

A cold producer installed a valid unpaged whole-prefix preview while it built
the durable spatial hierarchy.  Once hierarchy metadata named private pages,
`canDrawChunksAtCut()` delegated to the whole-prefix `canDrawCut()` query and
claimed those pages were already resident.  The provider skipped the required
representation transition, then failed when asked to prepare page layers from
the still-unpaged generation.

Rule: representation identity is part of a residency proof.  Empty page
selection denotes whole-leaf accounting only for an actual spatial generation;
a nonempty page selection can never be satisfied by an unpaged preview, even
when that preview can draw the same logical leaf and cut.  A page set with no
active faces at a valid coarse cut is a successful zero-draw presentation, not
a renderer failure.  Test preview-to-spatial replacement under true-cold task
ordering, not only warm page loading.

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

The retained source style and the framebuffer have distinct commit edges.  A
selection mutation therefore requests one exact presentation-only successor;
an older frame which was already in flight cannot acknowledge the newer style
revision.  That successor neither changes the retained LoD allocation nor
contributes renderer-capacity evidence.

### GUI controls and commands share one authority

Early primitive and polygon plugins held private copies and could disagree
with CLI state, other views, or retained manipulators.

Rule: GED/librt owns mutable edit state.  Controls and manipulators read after
typed revision events and submit typed operations to the same session.  Plugin
unload has a synchronous barrier which destroys filters, controls, dialogs,
and retained presenters before unloading the DSO.

### Command refresh needs a durable semantic witness

The qged view-settings replay initially passed on OSMesa but failed on System
GL: a faceplate command changed authoritative libbv state and the rendered
scene, while the open settings panel retained its old checkbox value.  The
post-command checkpoint compared the camera hash and a transient refresh-dirty
latch.  A sufficiently fast renderer could consume that latch before the
comparison, so unrelated later renderer work happened to hide the defect on
the slower backend.

Rule: command-to-widget notification must not depend on paint completion or a
consumable wakeup latch.  Compare the retained monotonic presentation revision
captured before command execution with the value after it.  Keep the dirty
latch for explicit redraws which do not mutate state, but never use it as the
sole proof of a semantic change.  Run the same command/widget replay through
both fast and slow renderers; backend-dependent success is evidence of a
missing ownership or revision boundary.

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

### A selective delta cannot drive a scene-wide successor

A warm Generic Twin zoom completed an exact 14-entry source-delta pass while
the motion-to-stable handoff needed a new allocation over all 709 mesh
occurrences.  The handoff correctly requested that allocation, but the
controller left the selective delta armed.  Each pump reconstructed the same
14-entry plan, skipped its already-satisfied entries, and then observed that
the full allocation was still unapplied.  The framebuffer was useful and the
worker service was idle, but the planning/handoff pair could never terminate.

Rule: a source delta owns one bounded source-local pass.  A scene-wide
retained-allocation successor atomically consumes that delta, clears its pinned
entry plan, resets pass annotations, and begins once from source zero.  Centralize
that transition so reconciliation and local-reduction paths cannot implement
different lifetimes.  `ObolCompletedPassOwnership.tla` now asserts that an
allocation-handoff effect cannot retain a selective delta; the cross-renderer
Generic Twin wire lifecycle is the executable regression.

A later 5k OSMesa run exposed the same ownership error on the result-demand
path.  The selective cursor had completed, but its target set was retained
through an exact frame.  That stale mode claimed planning ownership and made
the still-standing full demand appear selective, so no runnable cursor could
consume it.  Selective mode now retires at the cursor-completion edge, exactly
as `ObolSubmissionPass` specifies.  Both the runtime validator and the offline
trace checker reject a source-delta fact without an active submission cursor.

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

The first implementation of that correction still started the search before
the completed-pass snapshot could say that its first exact sample was pending.
On a large scene the allocator had applied every cut, but a temporary global
ceiling hid the candidate.  Repainting could not change either fact, so the
controller repeatedly presented the same coarse population.  Capacity-sample
ownership is now level-triggered: completed-pass selection predicts the first
sample from the current allocation certificate and requires the ceiling-free
handoff before constructing or advancing the search.

Rule: applying an allocation is not presenting it, and presenting it is not
measuring it.  Keep `ALLOCATING`, `PRESENTING`, and `MEASURING` as distinct
certificate phases.  A population change invalidates the candidate and yields
to the already-pending resident-growth producer.  A hidden candidate yields to
one presentation handoff.  Only an exact current framebuffer may consume a
timing sample or its bounded invalid-sample allowance.  Repainting is never a
successor for an allocation or presentation ownership gap.

### An older frame barrier precedes a newer capacity candidate

Lucy OSMesa reached `ALLOCATING` while an ordinary changed-population frame
barrier was still pending.  The broad "capacity search awaits a sample"
predicate routed completed frames around generic barrier consumption, while
the barrier itself paused submission.  Each owner therefore waited for the
other and qged presented unchanged HUD-only frames indefinitely.  Once the
barrier was consumed, its candidate successor initially reused the previous
allocation plan's semantic revision; the runtime trace checker correctly
rejected the two plan serials sharing one revision.

Rule: completed-frame handling consumes the older typed barrier before
dispatching search-phase work.  If a frozen candidate is already
`ALLOCATING`, barrier completion transfers directly to that exact candidate;
no heuristic may replace it and no generic presentation debt may survive.
The transfer publishes a fresh capacity revision before the successor plan
commits because the older population may already own a plan under the
selection revision.  Preserve this boundary with the focused coordinator
test, `ObolCapacitySearch.tla` plan-revision invariant, and the qged control-
trace uniqueness check.

### A bounded retained cursor keeps its allocation mode until retirement

A timing-perturbed 150k OSMesa replay committed allocation plan 1 and then
applied that plan over several bounded source windows.  Cut application
temporarily made the generic population-settled predicate false.  The
controller consequently changed the still-active cursor back to ordinary
admission and discarded its retained allocation certificate; plan 2 was then
committed for the same inventory, availability, visibility, view, policy, and
capacity revision tuple.  The allocator and per-occurrence cut publication
were deterministic—the cursor-mode transition had destroyed their owner.

Rule: the eligibility predicate may start a retained pass, but it may not
retire one.  A retained pass keeps its mode until the bounded cursor completes
or an explicit semantic invalidation retires that cursor.  Frame publication,
resident-demand bookkeeping, and intermediate population-settlement changes
may pause the pass but cannot reinterpret its unconsumed tail.  This is the
production expression of the complete-population transaction in contract rule
20.  `BObolLodSubmissionIntent::retainedAdmissionForPass()` centralizes that
edge, and the coordinator test covers both active-pass retention and idle-pass
retirement.

Post-fix warm OSMesa replays at
`/tmp/qged-lod-{50k,150k}-retained-pass-pin-20260901` pass the external
six-domain trace in 70 and 157 seconds.  Their plan serial remains unchanged
through every bounded application window; each later plan is paired with a
new availability or capacity revision.  Neither replay leaves a qged process
behind.

A later Hubble OSMesa hierarchy replay exposed the adjacent implementation
boundary.  Selection itself retired its exact presentation frame, but the
selected subpath erase/redraw left triangle recovery active while capacity
search selected a richer static candidate.  The retained-allocation request
encoded that candidate as an ordinary "preserve budget" request, so the
recovery ceiling silently clamped it.  The following exact frame therefore
described a different occurrence population, the certificate correctly
classified it as stale, and the controller restarted the same search.

Rule: a capacity-candidate allocation is a distinct, stronger request kind.
While its certificate is `ALLOCATING`, its exact candidate budget survives
recovery, deadline, throughput, and generic reallocation heuristics.  A
completed frame may narrow or reject it through the bounded search, and a
semantic population invalidation may retire it, but no weaker request may
rewrite it in place.  `ObolPointQualityOwnership` already expresses this
priority; the focused C++ refinement now activates a recovery ceiling while
applying the candidate, and the dual-renderer Hubble hierarchy workflow covers
the complete selection/erase/redraw effect path.

A later 50k trace exposed the observational companion to this rule.  After
the capacity revision advances, the renderer intentionally keeps the preceding
committed allocation visible while it calculates the replacement.  That old
plan is a safe presentation fallback, not a planning certificate for the new
revision.  Reporting the active renderer plan as the current control plan made
the trace checker appear to observe two plans under the new tuple even though
the first was certified under the preceding tuple.

Rule: retained allocation results carry the complete six-domain revision
stamp.  Diagnostics report a nonzero current plan only when that stamp matches
the controller stamp and the plan remains the renderer's active plan.  Report
the older committed plan separately so safe presentation continuity cannot be
confused with current planning authority.  The allocator input key includes
the same stamp, and the focused C++ test advances capacity with otherwise
identical numeric inputs to prove that the fallback becomes non-current and
one distinct replacement is published.

### Exact visibility changes allocation, not renderer capacity

An exact hierarchy erase/redraw initially advanced inventory and capacity even
though neither immutable geometry nor the renderer cost model had changed.
That discarded a useful certificate, drove a 150k scene to a much coarser cut,
and made recovery take tens of seconds.  Merely preserving the certificate was
not sufficient: launching replacement allocation before applying the sparse
visibility journal stranded newly visible entries, while allowing an old
allocation to finish after a second edit could publish a stale successor.

Rule: effective visibility owns a separate revision.  Apply its bounded source
delta first, request and complete one exact framebuffer classification for the
new census, then perform one scene-wide allocation.  The census proves which
occurrences are logically visible; it does not prove how those occurrences
are represented.  A hidden occurrence may have had its mesh payload retired,
so allocation based on the predecessor framebuffer can omit the restored
occurrence and leave no producer capable of repairing it.
Preserve capacity evidence across that edge, and supersede any unpublished
visibility allocation when a newer exact edit arrives.  Terminal liveness is
required once input becomes quiescent; continuous edits retain coherent
presentation but cannot promise a final allocation.

This distinction was missing from both the focused submission abstraction and
the controller.  A 50k erase/redraw replay therefore reached 50,000 logically
visible occurrences but represented only 49,984, with no task, owner, or
obligation.  The corrected formal sequence is `visibility delta -> exact
presentation/classification -> allocation`; the C++ transition uses the same
typed readiness predicate.  Do not merge the exact-census and exact-frame
facts again merely because they often complete in the same GUI event on small
or fast scenes.

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

### Recompute ownership after effects change producer eligibility

A 150k completion path classified an active submission cursor as the producer
of a pending point-calibration frame.  Calibration then paused that same
cursor.  No frame, provider, service result, or finite timer remained, but the
controller kept only the progressive-pump level and waited forever.  A second
variant occurred after a deadline correction: the pre-correction state chose
population continuation, the correction requested calibration and paused the
cursor, and the stale choice suppressed the explicit frame request.

Rule: producer ownership is derived from the post-effect state.  Inputs which
jointly pause or enable a producer belong in one named value, not positional
booleans copied between call sites.  If an effect changes any input, recompute
the successor before publishing it.  Exhaust the small Boolean boundary in an
executable truth-table test; the current point-calibration boundary checks all
256 combinations.

### An owner label is not a progress witness

The abstract control snapshot correctly reported `PRESENTATION`, so the
ordinary ownerless-work invariant remained green even when presentation had
no possible successor.  A scheduled pump did not help: its only candidate
cursor was paused by the presentation obligation itself.

Rule: validate owner and witness separately.  Presentation requires a queued
host frame, an independent producer capable of publishing one, or a finite
publication timer.  Runtime refinement validation now reports an explicit
unwitnessed-presentation violation.  TLA+ proves the abstract
`PhaseHasWitness` property; exhaustive value tests and sampled graphical
reports prove that the concrete-to-abstract mapping supplies a real witness.

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

### Companion fields are hidden transition relations

The submission cursor once reset its source index, entry offset, pinned source
plan, retained-allocation cursor, and visibility count in separate statements
at dozens of call sites.  Pass activity and rescan debt were likewise cleared
in either order.  Each sequence was locally reasonable, but every omitted or
reordered statement defined another reachable state: a new pass could inherit
an old source plan, or a retired pass could retain an unowned rescan.

The same hazard existed in evidence rather than work ownership.  A proven
renderer cut was writable separately from its view and quality-domain key;
retained error bounds were separable from their view/policy/point-threshold
key; GPU elapsed time was separable from its sample serial and retained-upload
replay state.  Partial updates could therefore authorize stale capacity or
quality evidence even though every individual scalar looked valid.

Rule: if values are valid, invalidated, or advanced together, represent them as
one allocation-free value with semantic transitions.  Cursor repositioning,
fresh pass start, pass retirement, evidence replacement, and renderer sample
acceptance are operations, not repeated assignments.  Keep explicit pause and
resume transitions only where debt deliberately survives.  Test key mismatch,
predecessor debt, first-versus-repeated upload, and reset behavior at the value
boundary before exercising the full controller.

### A candidate application is not capacity invalidation

A retained allocation selected by an active capacity search used the generic
changed-cut completion helper.  That helper preserved terminal certificates
but reset every nonterminal search.  The richer population consequently drew
once, lost its static-deadline owner, and ordinary steady planning immediately
coarsened it again.

Rule: applying the exact `ALLOCATING` candidate advances the same frozen search
toward presentation.  It preserves the search key, candidate budget, and
deadline goal while retiring only an obsolete generic frame barrier.  Capacity
is an output epoch and cannot be one of the semantic inputs used to decide
whether that search is current.  Unit-test the complete select, apply, present,
and sample boundary; a terminal-only regression does not cover it.

### Managed presentation is not external allocation cost

After erasing a selected Hubble subpath, the allocator alternated between two
otherwise current populations.  Its external-cost input was reconstructed by
subtracting the full occurrence-managed CAD population from the preceding
frame.  Point aggregation and renderer ceilings can replace that population,
so applying one allocation changed the next allocation's input key and erased
its certificate.

Rule: obtain external presentation cost from retained work outside the
occurrence-local allocation domain.  Applying CAD cuts or point protection
must change total and managed cost together and leave their difference
invariant.  Protect the boundary with the presentation TLA+ invariant, a view
state allocation-application regression, and a selected-path erase replay.

### Current allowance is not the capacity-search domain

The allocator reported its current optional marginal allowance as the maximum
searchable candidate.  A conservative throughput estimate thereby became a
self-certifying upper endpoint: the search could prove that estimate safe but
was forbidden to ask the same deterministic allocator for a richer scene
budget.  Static zoom views stopped visibly coarse despite ample deadline
headroom.

Rule: the finite search domain is the complete currently resident pixel-demand
cost.  The allowance used to build the current population is only the first
candidate.  Later candidates remain bounded by deterministic growth, sample
and candidate limits, persistent deadline ceilings, and the hard render abort.
Do not encode the current answer as the maximum possible question.

### A first hard miss needs a typed successor

A current allocation near the software-renderer deadline could miss before
its first exact frame.  The controller delayed bounded search until a prior
scalar ceiling existed.  Repeated direct corrections could then exhaust their
presentation handoff and leave a nonterminal 99-percent view with no work,
render request, constraint, or capacity owner.

Rule: when an unconstrained current allocation is applied, its first hard miss
is already a valid unsafe observation of that candidate.  Enter the bounded
capacity certificate immediately.  The preferred goal may reject it and
transfer once to the longer static-image deadline; either way the interrupted
frame has exactly one finite successor.  The lifecycle composition model now
asserts that this edge owns quality recovery and a new exact-frame obligation.

### A proof obligation needs an executable refinement witness

The terminal and host-work models already required a presentation obligation
to move from a controller pump to a queued render.  The first dense 50k trace
nevertheless reported an unwitnessed presentation owner because production
recognized only an already queued render, a timer, or an independent producer.
At an adjacent seam, a stale-population capacity result selected a replacement
plan without advancing the capacity revision, allowing two plans to claim the
same six-domain certificate.  Neither failure contradicted the abstract proof:
the missing steps were in the concrete-to-abstract mapping.

Rule: every abstract progress edge needs a named production predicate, a typed
diagnostic source, and an executable adversarial check.  Shared host PUMP is a
presentation witness only when the controller-local presentation projection
is one of its standing reasons.  Every replacement capacity candidate advances
the capacity domain even when another producer owns its mechanical allocation.
Dense traces exclude transaction/frame clocks from cycle identity and reject a
presentation owner whose independently exported witness mask is empty.

### A semantic classifier input must publish its revision at the setter

Hubble's software hidden-line path produced two valid but different retained
allocations under one unchanged six-domain revision tuple.  The allocator was
not nondeterministic: point aggregation had deliberately changed from one to
64 pixels between the plans.  The threshold was present in the allocator's
canonical input key, but its centralized publisher had not advanced the
capacity revision.  The abstract capacity model already treated classifier
rule and threshold as immutable members of its revision; production had failed
to refine that edge.

Rule: a semantic input and its revision are one publication operation.  Point
presentation and discovery thresholds now advance the capacity domain inside
their centralized setters, including the paired setter, so no calibration or
structural-recovery caller can omit the edge.  Runtime traces continue to
reject two nonzero allocation plans for one unchanged revision tuple.  Do not
weaken that check by adding mutable numeric fields only to diagnostics.

The same Hubble trace exposed a neighboring refinement error in the runtime
producer validator.  A complete bounded inventory/submission coverage pass is
a valid producer for the retained importance census; an arbitrary selective
submission is not.  The validator now expresses that distinction instead of
requiring only a demand-refresh producer.

### Specialized pass completion cannot retire broader demand

A later Hubble OSMesa hidden-line trace reached a missing-coverage completion
while a retained importance census was pending.  The coverage path cleared the
level-triggered demand refresh even though its cursor had run as retained or
structural work and had not proved the complete current-view demand.  For one
controller observation the census consequently had no producer; the next pump
happened to recreate a pass, concealing the violation in ordinary final-image
tests.

Rule: physical-demand refresh has one retirement edge: completion of a dense,
ordinary, non-selective demand pass.  Coverage, retained allocation, and
structural repair retire only their own obligations.  A broader demand may
remain paused behind a stronger owner, but it stays level-triggered and visible
to the refinement checker until its own cursor completes.

The same trace also showed why allocation identity must describe the selected
endpoint rather than every attempted option.  Once an allocation already
selects complete pixel demand, enabling a protected-floor trial cannot enrich
it further.  Reuse that canonical endpoint when all other semantic inputs are
identical; issuing a second plan serial under the same revision is not useful
replanning, it is a contract violation.

### A successor certificate retires predecessor control debt

A bounded capacity search could publish an exact terminal certificate while an
older motion/deadline handoff still retained its reconciliation budget.  The
next generic pass then consumed that older, lower budget and replaced the
newly certified population.  This looked like selection instability because a
presentation-only selection frame was often the first frame after the latent
handoff became runnable.

Rule: accepting a terminal certificate is an ownership transfer, not merely a
new fact.  Atomically retire every predecessor handoff budget and prohibit
reconciliation from being rearmed while the certificate is published.  Keep
selection/style revisions outside the capacity key and test that their one
exact presentation frame returns to the same terminal allocation.

### Render cost equality does not imply PoP cut equality

Two adjacent PoP cuts can have identical vertex and primitive counts while the
richer cut improves coordinate precision.  Treating all improvements as
positive-cost marginal decisions stranded those free cuts whenever the numeric
budget was exactly full.

Rule: canonicalize every selected mesh cut through all consecutive equal-cost
successors before publishing it.  Keep only positive-cost upgrades in the
marginal queue.  Formal allocation-prefix checks and executable tests must
cover a saturated budget with a richer zero-cost precision cut.

### A GUI checkpoint is an observation barrier

OSMesa image export necessarily renders the requested state, but a fast
System-GL canvas can still expose the previous front buffer when Qt coalesces a
retained-scene update.  Waiting only for command completion or controller idle
therefore made cross-renderer tests disagree even though both scene graphs were
correct.

Rule: after a semantic scene mutation, drain progressive work and require an
actual presentation before capturing the checkpoint.  Keep passive quality
samples passive; only an explicitly named exact-frame barrier may force the
presentation.

### Native geometry must be stable at the semantic command

A resize replay once observed the requested geometry, waited successfully,
and then received a delayed X11 configure before draw/autoview.  The resulting
camera differed only by the new aspect ratio and looked like a policy-dependent
CAD framing defect.

Rule: require a bounded stable interval before the initial draw and record the
native geometry at the camera-defining command.  A test that proves only the
earlier resize acknowledgement has not proved the autoview precondition.

### Shared layers are not exclusively owned containers

The framebuffer regression asserted that its shared overlay layer had exactly
one child.  Other retained image presentations may legitimately use that same
layer, so the assertion failed even though the migrated framebuffer viewport
was attached exactly once and the controller state was correct.

Rule: validate a feature's owned node by identity and occurrence count.  Do not
infer exclusive ownership from a shared scene-layer child count.

### A replacement is not atomic if its predecessor is removed during staging

The mesh-LoD cache used immutable content-addressed payloads and published a
small name-to-content mapping last, but forced refresh deleted the old mapping
before it began generating the candidate.  A later allocation or cache-write
failure therefore left no partial new hierarchy discoverable—while also making
the complete old hierarchy unreachable.  The usual “publish marker last” rule
had protected first publication but not replacement.

The retained-scene path exposed the notification version of the same mistake.
A denied compact mutation changed no scene data, yet it entered and closed an
outer update scope before rejection; observers saw a new node identity for the
old scene.

Rule: a transaction includes its discoverability and notification identity.
Stage the replacement while the predecessor remains authoritative, perform
resource admission before opening any observable update scope, and commit the
payload/map or scene/notification pair together.  Failure tests must start
with a known-good predecessor and compare that predecessor after denial; an
empty-state test alone cannot expose destructive replacement.

### Skipping zero is still wraparound

Several asynchronous identities incremented an unsigned integer and replaced
zero with one.  That preserved a sentinel but eventually reissued every old
nonzero credential.  A separate trace path reset tokens when diagnostics were
re-enabled, creating the same ABA risk without numeric overflow: a late RAII
scope could close a new frame carrying the recycled token.

Rule: classify every increasing value by authority.  Authentication,
ordering, and ownership credentials use one checked successor and fail-stop
before reuse; clearing a journal does not reset the lifetime of its tokens.
Observation-only totals saturate and must never enter an authorization
comparison.  Test both the representation boundary and lifecycle reset paths,
because either can reissue an apparently current identity.

### An AABB cannot certify an oriented volume

A producer-side precheck required a PCA-oriented proxy to contain every corner
of the geometry's axis-aligned bounds.  That test rejects nearly every useful
rotated box: an AABB contains empty space outside the tighter OBB.  Removing the
precheck without a replacement would have created the opposite failure, where
a structurally valid but non-conservative OBB could reject the whole mesh at
the renderer boundary.

Rule: validate oriented containment against authoritative renderer positions,
not the corners of an axis-aligned envelope.  Keep strict admission as the
default, but let a producer explicitly classify the proxy as optional so only
`InvalidAggregateProxy` falls back to ordinary bounds.  Test both directions:
a tight rotated OBB must survive, and a well-formed OBB which excludes a real
vertex must be discarded without losing the mesh.

### A digest is not an evidence stamp

Several renderer and controller caches folded multiple revisions or a variable
input list into a 64-bit hash and reused work when the hash matched.  The
inputs were valid and the hashes were well distributed, but neither fact makes
equality a proof: a collision can authenticate stale colors, visibility, draw
representation, progressive geometry, capacity evidence, a structural
frontier, or compact-source presentation state.

Rule: authorization compares the complete typed tuple.  Hashes may select a
candidate bucket only when exact equality follows.  For bounded variable input,
retain and compare the canonical vector rather than inventing a stronger claim
for its digest.  This is the implementation analogue of a TLA+ record: omitting
or folding one field changes the state relation, even if the failure is hard to
reproduce statistically.  When an allocation-free scalar must cross a hot
policy boundary, issue a checked non-reused token only after the owner compares
the exact value.  A bounded exact-interner eviction may cause redundant work;
it must never allow a token to name a different retained value.

### Aggregate serials are identities, not arithmetic

Combining per-object monotonic serials by hashing fixed one collision mode but
left probabilistic equality; combining them by addition avoided ordinary
collisions but could saturate while each source remained live.  Using the
object address as the source component added a separate ABA path when an
allocator reused a retired node's storage.  All three are invalid when the
aggregate controls frame classification or measurement deduplication.

Rule: give each source a non-reused lifetime identity, compare the canonical
set of `(source identity, source serial)` records exactly, and issue a checked
non-reused aggregate token only when that set changes.  Retire the current set
when evidence is absent; otherwise `A -> absent -> A` can silently reuse the
old aggregate identity even though no numeric counter wrapped.

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
