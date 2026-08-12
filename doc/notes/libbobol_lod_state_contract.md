# libBObol progressive display state contract

This document is the compact, proof-aided contract for BRL-CAD progressive
drawing.  It supplements the ownership design in `brl_obol_view_lod.txt`.
Implementation details are correct only when they refine this model without
violating its safety, liveness, or cost properties.

The contract is intentionally expressed in terms of state and transitions,
not classes.  qged, MGED, and gsh may drive the state differently, but they
must produce the same visible behavior.

## User-visible objective

For every view:

1. Present the cheapest useful representation immediately.
2. While the view changes, present a coherent mesh cut that meets the
   interaction frame target.  Scale-changing input refines from the current
   cut and publishes newly available suffixes while input continues; it does
   not wait for the quiet transition to restart demand from a proxy.
3. Once the view is quiet, monotonically refine visible geometry toward a
   one-pixel screen-error target, subject to a finite stable frame budget.
4. Declare stability only when every visible mesh-backed occurrence has its
   richest affordable cut and no event is required to make progress.
5. Retain and trim shared data independently from the active draw cut so a
   camera change normally changes prefix counts, not geometry objects.

Boxes are a cold/missing-data fallback, not a normal interaction LoD.  A box
may remain only for a visible occurrence whose minimum mesh prefix is not yet
available or whose source has definitively failed.  OBB is an optional cached
fallback and is not a required progression state.

## State domains

The model distinguishes four identities which must never be conflated:

| State | Identity and owner | Cardinality |
|---|---|---|
| source | BRL-CAD object plus source revision; shared scene | one per drawn source |
| asset | immutable PoP hierarchy and retained prefix; LoD service | one per distinct mesh content |
| occurrence | source path/instance transform and view-local requested/active cut | one per displayed leaf instance |
| presentation | retained part/GPU buffers plus instance records; Obol | one per source/view presentation |

One asset may serve many occurrences.  Each occurrence has an independent
active cut.  A presentation references assets and occurrences; it does not
own BRL-CAD source truth or choose LoD.

For a view `v`, occurrence `o`, and asset `a`:

```text
visible[v,o]          boolean from the renderer-contract frustum test
pixel_extent[v,o]     projected screen extent in pixels
error[a,c]            producer-certified object-space error at cut c
screen_error[v,o,c]   error[a,c] * pixel_extent[v,o] / asset_extent[a]
demand[v,o]           least admissible cut c with screen_error <= target
active[v,o]           currently presented PoP cut, or NONE
resident[a]           richest retained array prefix
resident_target[v,a]  richest visible pixel demand for asset a in view v
drawable[a,c]         resident arrays contain all entries required by cut c
faces[a,c]            cumulative face prefix for cut c
scene_budget[v]       calibrated aggregate submitted-work allowance
render_ceiling[v]     temporary O(1) interaction-only maximum cut, or NONE
```

Cut ordinals have no geometric meaning.  The producer may define any monotone
sequence of cumulative prefixes up to the explicit implementation limit, and
every cut carries its own face count, point count, resident byte estimate,
per-axis quantization precision, conservative object-space error, and exact
flag.  Selection binary-searches that metadata.  No scheduler or renderer may
reconstruct quality with `2^-cut`, assume a fixed cut count, or treat a
particular ordinal as full detail.

Asynchronous request identity and producer-selected demand are deliberately
separate.  `BObolLodRequest::requestedCut` is immutable task-key input.  It is
an admissible cut when a retained hierarchy was available during projection,
or `-1` when the cold submitter could express only projected diameter and
target pixel error.  The provider resolves that physical demand against its
certified hierarchy and reports it in `BObolLodResult::resolvedCut` without
rewriting the copied request.  `geometry.activeCut` is the cut in the returned
presentation and `residentCut` is the retained-array frontier; neither is a
substitute for resolved demand.  Result matching, cancellation, coalescing,
and duplicate suppression use the original request, while presentation,
terminal-state, and refinement decisions use `resolvedCut`.  A local reuse
path which already owns the hierarchy performs the same resolution before it
binds a second occurrence.

The canonical libBObol mesh producer uses 16-bit source quantization and one
axis refinement per cut.  It begins at one bit on every nondegenerate axis,
then increments the axis with the largest remaining object-space error until
all axes are exact.  Thus a fully three-dimensional mesh has 46 cuts; planar
or linear meshes have fewer.  This is a producer policy, not an Obol renderer
assumption.  Obol accepts any producer-certified monotone cut vector within
its explicit limit.  The named unspecified-cut value means “richest resident”
and is outside the valid cut range; no valid ordinal is a sentinel.

In an orthographic view, equal projected extents have equal LoD value
regardless of depth.  In a perspective view, distance matters only through
projection.  No separate depth weighting is allowed.

## Draw extent and autoview contract

Each progressive draw source has one authoritative source-local coverage
extent for `(full draw path, source revision)`.  It is a producer fact, not a
measurement of whichever fallback or PoP cut happens to be visible.

1. Cold coverage computes the extent while enumerating all eligible
   occurrences and publishes it once through the stream's priority lane.
2. The owner thread may fit the camera only after that final priority record
   has been drained and merged.  “Worker queued” is not “scene published.”
3. A warm manifest stores the identical path-scoped extent.  It may declare
   warm coverage authoritative only when that extent and every immutable leaf
   source request are present.
4. Object-name AABB/OBB caches remain object-local.  A transformed occurrence
   union must never overwrite them; the same leaf may be drawn through several
   paths with different transforms and Boolean contexts.
5. During initial progressive draw the provider is the sole camera owner.
   Transaction code must not apply a provisional fit from partial boxes.
6. One request produces at most one framing transition and may never roll back
   to an earlier partial extent.  Cold/warm and System GL/OSMesa runs of the
   same source/view must produce the same final camera.

The initial blank/default frame to exact-model fit can remain visible on a
slow backend.  It is one intentional transition; a partial-fit/exact-fit pair
is a contract violation.  A transient blank presentation is not an idle
endpoint: if source or progressive work is pending, every backend must arm its
provider timer before returning from that paint.  In particular, the software
canvas may not make image availability a prerequisite for the wakeup which
will produce its first image.

## State machine

```text
                         source/result failure
                                |
                                v
  ABSENT --> FALLBACK --> COVERAGE --> INTERACTIVE --> SETTLING --> STABLE
               ^            |             |              |           |
               |            |             +--view event--+           |
               |            +-------source invalidation--------------+
               |                                                    |
               +------------resident eviction below minimum---------+
                                                                    |
                                      quiet + idle + demand snapshot |
                                                                    v
                                                               COMPACTING
                                                                    |
                                                                    +--> STABLE
```

`FALLBACK`
: Source bounds are visible immediately.  Asset loading may be pending.

`COVERAGE`
: Every visible mesh-backed occurrence receives its minimum coherent PoP
  prefix before any already-covered occurrence receives optional refinement.
  A scale-only camera refresh is not cold coverage: it retains the preceding
  cuts and may refine them while its bounded scan verifies current visibility.
  New source/inventory coverage retains the strict minimum-mesh-first rule.

`INTERACTIVE`
: A finite scene budget and an O(1) presentation ceiling bound the next frame.
  The exact per-occurrence allocator catches up in bounded waves.  Finer
  resident arrays are not evicted.  Pose-only input may hold richer worker
  publications to avoid command churn; zoom must apply them in bounded waves
  and request missing cumulative suffixes without ending the interaction.
  Render-cost exhaustion does not prohibit a zoom suffix request: residency
  is admitted by the independent worker working-set and resident-byte limits,
  while the result publishes the currently affordable presentation cut.  A
  completed scale frame may expose one richer population before the next
  input event, using a bounded 10 Hz quality-frame allowance.  This remains
  active while wheel/trackpad events continue; it does not depend on the quiet
  debounce or button release.  A successfully presented quality ceiling is
  retained across later scale epochs.  A deadline miss or a frame slower than
  100 ms lowers only the O(1) render ceiling; it never rewrites a retained
  occurrence to its minimum cut or discards the richer resident prefix.
  System GL motion decisions use the slower of completed endpoint traversal
  and asynchronous CAD GPU time; OSMesa supplies the same evidence through
  completed traversal time.  An aborted Coin traversal has no exact submitted
  work denominator.  It may drive deadline backoff, the global render ceiling,
  or sub-pixel aggregation, but it may not reduce calibrated throughput,
  persistent scene budget, occurrence cuts, or residency.

`SETTLING`
: The interaction cut remains until one coarse frame has actually been
  presented and the quiet debounce expires.  Pose retention is legal only if
  the interaction began from a view already proven ready.  Initial draw,
  partial coverage, and unfinished refinement may never become terminal merely
  because an orthographic camera event was pose-only.  The epoch separately
  records “scale changed at least once” and “the current camera edge changes
  scale”; a later rotation cannot rearm a zoom probe, and a scale round trip is
  tested by its final signature rather than by the presence of wheel events.
  On entry to quiet state, the last proven stable scene allowance replaces the
  transient 60 FPS motion allowance.  A post-release frame which installs a
  new interaction ceiling before debounce must reapply the stable snapshot or
  enter a witnessed handoff; it cannot report ready with that ceiling stuck.

`STABLE`
: Refinement proceeds in bounded publication waves.  Candidate marginal
  transitions are ordered by user emphasis and estimated screen-error
  reduction per added face, but an occurrence may jump across several authored
  cuts when its richest affordable drawable prefix is already resident.  It
  must not pay one frame or worker round trip per cut merely because the
  producer supplied a denser quality schedule.  Stable cuts do not regress
  within an unchanged view epoch.  If
  a newly attempted discrete cut proves too expensive, overload recovery
  begins from the current presentation: minimum coverage is reserved for all
  visible occurrences and the richest already-active cuts which fit are kept
  in screen-priority order.  Recovery must not reset every occurrence to its
  minimum and then replay initial cut walking.

`COMPACTING`
: After a longer quiet interval, aggregate consumer demand may trim resident
  CPU/GPU tails with headroom.  Active cuts remain drawable throughout.  GPU
  atlas maintenance may consolidate compatible live ranges with one bounded
  page-sized device-local scratch buffer, update part offsets atomically, and
  delete emptied pages.  It must never require a second scene-sized atlas or a
  source reload.  Consolidation is permitted only after resident demand has
  been quiet long enough to make the current cuts a stable witness.

Persistent PoP cache records are activation-cut chunks.  Growing a retained
asset after compaction reads only `(resident, requested]`; it must not recreate
a cache-reader copy of `[minimum, resident]`.  The worker constructs and
atomically publishes a replacement immutable renderer snapshot while the old
snapshot remains drawable.  Down-trimming allocates the exact retained prefix
and must not preserve a richer vector's capacity.  Authored corner-normal
meshes currently use the whole-prefix fallback because globally consistent
vertex splitting is not derivable from an unannotated suffix.

Each ordinary coordinate/index stream also carries a nonzero process-local
lineage token.  Equal lineage across immutable geometry generations certifies
that the preceding position, normal, and index values are an exact prefix and
that quantization domains and cut tables are unchanged.  Obol may then keep
the existing GPU prefix, submit only the newly published CPU suffix, and use a
device-local buffer copy when reserved capacity must grow.  A zero or changed
token requires conservative complete replacement.  GPU compaction may retain
a richer prefix than CPU residency until its own quiet maintenance pass; draw
counts still clamp to the active/resident cut, so that retention is invisible.

Only one resident-growth task may be active for a given asset occurrence.
Wheel epochs with successively richer demands coalesce behind that task; its
completion wakes planning, which submits the newest still-missing target.
Serialized cuts must not queue as independent work because they cannot grow
the same asset in parallel and their obsolete result epochs create publication
churn.

A service request owns its duplicate-suppression identity from admission until
the owner thread drains its result, not merely until the worker finishes.  A
completed result waiting in the coalesced publication queue is still active for
submission purposes.  Releasing the identity at worker completion permits a
fast warm-cache result to be resubmitted once per planning wave while GUI
publication is deliberately batched, causing needless work and visible
convergence cycling.

A resident-memory denial is scoped to the resident-admission revision at which
it was observed.  For that unchanged revision it is a terminal capacity fact,
not an actionable unsatisfied occurrence, and an identical scene-budget
retarget must not erase it.  The denial becomes actionable again only after a
capacity/residency revision changes, or after the occurrence's requested
asset, view, or policy identity changes.  This makes a genuinely saturated
scene quiet without preventing progress when eviction or a larger limit later
creates headroom.

The owner-thread coordinator is the authority for these phases.  Client
Booleans are observations, not alternate state machines: named gesture and
cancellation events own their corresponding edges (gesture release retains the
interactive phase through its quiet debounce), a nonzero pending witness
cannot produce `STABLE`, and compaction observed before complete coverage is
canonicalized back to `COVERAGE` or `FALLBACK`.  The raw observation is retained
for invariant diagnostics, while only the canonical state is published.

Completed-frame CPU and GPU pressure are coordinator inputs, not merely HUD
telemetry.  Every unique active CAD assembly is sampled once; all policy and
diagnostic queries consume the retained aggregate in O(1).  A pressure edge
increments a resource revision and requests one compaction pass after coverage
and interaction safety gates.  If the same pressure remains after that
revision is handled, the terminal state is stable and memory-limited rather
than an endless COMPACTING loop.  A later pressure edge creates a new revision
and a new recovery opportunity.  If pressure is first observed while recovery
is disabled, enabling recovery consumes that still-current unhandled revision
once; the system must not wait for pressure to clear and recur.

The System GL atlas has an additional bounded-presentation response.  Exact
preparation orders visible unique parts by projected importance, reserves a
coherent minimum prefix before enrichment, and replaces only the
least-important eligible tail with one-point occurrence records when the
configured atlas cannot hold every minimum.  Those records are a terminal
memory-limited representation, not structural boxes and not unfinished LoD
work.  Selected/highlighted importance participates in the same ordering.
They are legal only while the completed-frame pressure sample, convergence
proof, and idle service state agree; clearing the pressure observation is not
a liveness requirement for an oversized visible working set.

Coverage necessity is part of the coverage policy, not a separate controller
latch.  LoD-disabled or service-detached views are vacuously covered.  When
coverage is required, invalidating its proof immediately gates compaction;
an old planning cursor cannot claim `COMPACTING` concurrently with missing
minimum coverage.  If cancellation, source erase, LoD disable, or an empty
view removes the coverage producer, retiring the associated compaction request
also acknowledges that pressure-recovery revision.  It may not leave an
unwitnessed background obligation.

Minimum-mesh coverage, user-visible convergence, aggregate scene budgeting,
view-demand scheduling, renderer-independent quality calculation, resource-
edge/recovery state, the late calibrated-headroom retry witness, and result-
publication batching are allocation-free coordinator policies.
BObolViewController supplies source, worker, completed-frame, and memory
measurements and executes scheduling decisions, but it does not maintain
parallel latches or reinterpret their one-shot rules.  The coverage policy
accumulates saturating visible/covered counts across bounded compact-index
windows.  A complete non-rescan pass atomically publishes its coverage result
and authoritative visible denominator before retiring those transient
counters.  Camera and inventory changes invalidate that proof; a quality-only
pass may preserve it and avoid a whole-hierarchy recount.

The convergence policy is a scene-pointer-free projection of those retained
proofs and progress witnesses.  It alone decides ready/background status,
memory- and performance-limited terminal states, the public HUD phase, and the
monotonic fraction for one typed view/policy epoch.  A diagnostic query may
sample this value but cannot introduce another scheduler, hierarchy traversal,
or alternate interpretation of readiness.

The budget policy treats projected per-occurrence error as demand and measured
renderer cost as total-scene admission.  It owns the current allowance, the
remaining budget shared by every bounded occurrence window, retained-cut
recovery admission, the one-shot overload witness, and the complete stable
calibration series.  A budget-blocked pass returns one decision describing
probe eligibility, the frame barrier, and whether the next frame is an
unchanged calibration replay or a presentation of newly admitted work.  Frame
completion consumes that witness and either requests the next bounded sample
or restarts submission; the controller has no parallel rescan, probe-count, or
terminal-budget latch.  Its stable growth tiers (4x, 2x, and 1.25x), three-
sample minimum, twelve-sample headroom bound, near-target convergence margin,
10 Hz zoom-quality probe, and measured overload backoff are deterministic
boundary-tested decisions.  Deadline feedback may reduce the current
allowance immediately, but no budget operation rebuilds or discards an
immutable PoP population.

The view-demand policy owns the camera-local distinction between pose and
scale changes, the scale-refresh proof carried through bounded coverage
windows, and the complete interactive quality-probe lifecycle.  One completed
scale frame may arm one coherent next-population probe; refinement/publication
and motion-frame barriers prevent early consumption, and a later resident-
growth completion may rearm exactly one probe.  A successfully presented probe
retains its ceiling across the next scale event only when it met the policy's
single 100 ms/10 Hz responsiveness contract.  Quiet handoff atomically retires
the scale epoch and returns whether stable policy revision must preserve its
demand refresh.  Clearing the last completed-frame barrier explicitly wakes
the host even when the generic progressive latch is already set; a subsequent
scale event therefore cannot starve a now-safe probe by consuming the only
notification edge.  The controller supplies active-cut and frame
observations and applies the returned renderer ceiling, but it has no parallel
pending, active, presented, quality-budget, or scale-refresh latches.

The scalar quality policy converts completed-frame timing into the interactive
pixel error, small-occurrence aggregation threshold, and reversible discrete
PoP ceiling.  It is scene-pointer-free, sanitizes non-finite inputs, applies
repeat overload correction relative to the cut which actually produced the
frame, and rounds toward a known-safe coarser PoP population.  The controller
may supply measurements and apply those values, but it cannot maintain a
second version of their clamps, thresholds, or cut arithmetic.

Every batch of applied immutable results retains its first-unpresented deadline
until the publication policy replaces that timer with exactly one requested-
frame witness.  The first useful mesh and the end of a result stream publish
immediately; otherwise stable and interactive batches use bounded, render-time-
adaptive intervals.  A completed frame atomically retires the whole
accumulated batch, and convergence cannot report ready before that
presentation.

Headroom admission thresholds are strict: demonstrated capacity must exceed
the current allowance by 5 percent and the reusable frame must take less than
120 percent of its target duration.  Stable-presentation handoff and
point-proxy calibration may consume the frame immediately following a terminal
planning pass, so the controller revisits the same one-shot admission when
those barriers clear.  If the interaction ceiling prevented a retained
occurrence from advancing at all, stable handoff converts that
unsatisfied-demand annotation into an explicit post-frame rescan.  The first
ceiling-free frame therefore remains bounded, and no raw demand bit may keep
the pump alive without a render, submission, worker, result, or timer witness.

Resident compaction is admitted by a separate allocation-free coordinator
policy.  Its quiet deadline is necessary but not sufficient: minimum mesh
coverage and every visual-settling witness (coverage/refinement scans,
calibration, the explicit late-headroom replay, stable handoff, result
publication, and batched unpublished results) must be complete.  Partial
service planning continues only when its resident-demand revision still
matches.  LoD-disabled drawing has no mesh-coverage prerequisite.  If active
LoD coverage is incomplete but has no worker, submission, provider, frame, or
timer capable of advancing it, the compaction request is retired; it cannot
keep the progressive pump or HUD pending by itself.

## Transition obligations

Every pending flag has exactly one progress witness:

| Pending condition | Required witness | Completion |
|---|---|---|
| source realization | provider task or source callback | fallback or source data published |
| minimum-mesh coverage | bounded compact-index cursor or full-rescan obligation | complete pass publishes visible/covered proof |
| asset suffix | queued/in-flight/cache task | result applied or definitive failure |
| applied result batch | adaptive timer, then exactly one requested frame | completed frame presents the batch |
| resident retarget | guaranteed requested frame | cut presented |
| budget blocked with measured headroom | guaranteed calibration frame | allowance grows or capacity is revised |
| quiet debounce | host timer/frame pump | stable policy revision |
| compaction | owner-thread compaction pass | retained prefix and demand snapshot agree |

It is invalid to report pending work when all of the following hold:

```text
no queued, delayed, in-flight, result, or cache-write work
no requested frame
no armed timer
no source callback capable of advancing state
```

That is a scheduler deadlock even if the GUI remains responsive.

For a retained asset, refinement is cumulative and budget bounded:

```text
admit(v,o) = richest cut c such that
             active[v,o] < c <= demand[v,o]
             && drawable[asset(o),c]
             && marginal_cost(active[v,o],c) <= remaining_scene_budget
```

The allocator may jump directly to `admit`; it must not jump directly to
`resident` when the view does not demand or the scene cannot afford that cut.
Under contention, fairness ordering may deliberately admit only the next
marginal transition before reconsidering peers.  If no richer cut fits, the
allocator may grow the probed allowance toward calibrated capacity or declare
the cut budget-limited.  It must not retry the same unaffordable transition
forever.

When a probe does not change geometry because the next prefix is discrete,
the next allowance grows from the previous probed allowance, not from the
unchanged active face count:

```text
probe_budget[n+1] > probe_budget[n]
probe_budget[n+1] <= calibrated_capacity
```

The number of unchanged probes is therefore logarithmically bounded.

## Safety invariants

These must hold after every owner-thread transition:

1. `active[v,o] != NONE` implies `visible[v,o]`, except during the bounded
   transaction which removes newly off-frustum bindings.
2. `active[v,o] <= resident[asset(o)]`, or `drawable[asset(o),active[v,o]]`
   for a coordinate-only plateau.
3. Shaded, wireframe, bounds, and picking use the same effective active cut,
   including the temporary interaction ceiling.
4. An occurrence binding references the current source, view, policy, and
   occurrence revisions.  Stale results cannot replace current data.  A
   camera-stale cumulative suffix may be rebased only when asset pointer,
   database revision, source revision, content hash, draw mode, source route,
   and occurrence identity all match; the retained occurrence's old camera
   stamp is not an asset-identity failure.  Rebase retargets metadata to the
   current allocator epoch.  Worker publication may extend residency or admit
   a richer reserved cut, but may never lower an already-richer occurrence
   cut.  The prepared renderer geometry defines the drawable frontier.
5. Stable refinement is monotonic for an unchanged view/policy epoch.
6. All visible mesh-backed occurrences receive minimum coverage before any
   occurrence receives a second optional prefix increment.
7. Off-frustum occurrences issue no draw calls and consume no scene face
   budget.  Their shared asset may remain resident.
8. Selection/highlight changes priority and style, not asset identity.
9. A source edit invalidates affected content; a camera edit does not.
10. Direct GL rendering restores GL state and invalidates Coin state caches
    before returning.
11. Cross-generation GPU prefix reuse requires an equal nonzero producer
    lineage; the ordinary and atlas upload counters must show no cumulative
    CPU re-upload while such a stream grows.
12. Scene render-cost admission governs `active`, never `resident_target`.
    Zoom residency is bounded by worker and resident-memory admission, and a
    prefetched result may not raise `active` above its independently reserved
    presentation cut.
13. Active camera input never rewrites retained occurrence cuts downward.
    Responsiveness pressure changes the renderer-wide effective cut; quiet
    admission and memory maintenance are the only cut/retention authorities.
14. A progressive autoview is fulfilled only from the acknowledged exact
    path-scoped coverage extent.  That certified source-local extent remains
    authoritative after every partial registry merge and final detached
    adoption; a representation-derived union may be used only before the
    certification exists.  Cache state, publication order, and renderer
    cadence cannot alter the final camera.
15. Wanting a richer retained cut is not a progress witness.  A pending phase
    must name a task, result, bounded cursor, timer, or presentation frame which
    can discharge it; otherwise the current performance-limited cut is a stable
    terminal state until a new view, resource, or residency edge arrives.
16. A renderer may report persistent atlas-admission pressure in a terminal
    memory-limited state.  That state is legal only with complete coverage, a
    view-ready convergence proof, no pending visual work, and a coherent
    drawable cut.  Pressure replacement points are permitted only for the
    renderer's least-important eligible tail and only while that same terminal
    proof is true; structural fallback boxes are never permitted.  Clearing
    the observation is not a liveness requirement when the visible working
    set genuinely exceeds the allowance.
17. Duplicate suppression spans worker execution and queued publication.  A
    result waiting for owner-thread drain owns the same request identity as its
    producer; coalescing may replace its payload but may not make that request
    appear absent.
18. A current resident-admission denial is excluded from the actionable
    unsatisfied frontier.  An identical render-budget retarget preserves that
    witness, and a new resident-admission revision reopens it.
19. Worker result identity is immutable.  A cold request with
    `requestedCut == -1` remains byte-for-byte matchable after the provider
    selects a cut; the selected demand is carried only in `resolvedCut`.
    Progressive ready results must use that resolved demand for terminal and
    active-cut decisions.
20. A bounded retained-allocation pass takes one complete population census.
    Its active and minimum render-cost currencies are frozen until every
    bounded source window in that pass has been consumed.  A presentation,
    result wave, or window boundary may reset per-frame work allowance but may
    not rescan the whole population or change the pass currency.  Completing,
    cancelling, or invalidating the pass clears the snapshot.
## Liveness properties

Using `[]` for “always” and `<>` for “eventually,” the implementation must
satisfy:

```text
[] (visible_mesh_backed && minimum_resident
    -> <> (active >= minimum || definitive_failure))

[] (scale_interaction && resident_target > resident
    -> <> (resident increases || memory_limited || definitive_failure))

[] (scale_interaction && richer_drawable && quality_headroom
    -> <> (effective_cut increases before quiet))

[] (pending && affordable_next && host_running
    -> <> presented_next)

[] (quiet && budget_blocked && calibrated_headroom
    -> <> (budget_increased || capacity_revised || next_presented))

[] (view_changed
    -> <> (all_active_bindings_match_current_visibility))

[] (stable
    -> no_pending_work && no_progress_witness && no_affordable_next)
```

Cancellation, result reordering, cache hits, cache misses, and source
replacement must preserve these properties.

## Budget and priority contract

The target demand is per occurrence; admission is per scene.

```text
stable demand:       pixel_error = 1
interactive demand:  adaptive pixel_error >= 1
scene render cost:   sum(cost[asset(o), effective_cut(v,o)])
effective cut:       min(active[v,o], render_ceiling[v]) when ceiling is set
```

The calibrated stable and interactive capacities are separate.  A quiet
readback or a motion presentation interval cannot contaminate the other
capacity.  A new interaction may conservatively seed an uninitialized or
underloaded motion estimate from measured stable geometry throughput; actual
slow motion frames lower it immediately.  A fast, extremely coarse frame is
not a lower throughput proof.  Likewise, a deadline miss while only a tiny
prefix is active may be widget, compositor, or readback overhead; fixed cost
divided by a handful of triangles must not poison either stable or motion
capacity.  Budget reductions affect a future view/motion epoch; they do not
make an unchanged quiet frame visibly regress.

Prepared System GL command replay and an unchanged retained ordinary-VBO
presentation are both reusable calibration witnesses.  The latter is proven
by unchanged cumulative full/suffix upload counters and is required by the
OSMesa fixed-function path.  A reusable frame which demonstrates headroom
after the normal probe series gets one bounded admission retry for the exact
`(view revision, policy revision, current budget)` witness.  The terminal
budget pass must explicitly request this unchanged replay and report it as
pending convergence; a later selection, HUD, checkpoint, or other unrelated
repaint is not a substitute progress event.  The view cannot report `STABLE`
until the witness is consumed or invalidated.  If that retry enlarges the
budget another distinct witness may retry; an unchanged budget cannot repaint
forever.

Stable refinement is monotonic within one camera epoch.  A calibration pass
which changes no visible prefix is a probe, not progress; at most twelve such
probes may occur at one unchanged active population.  This bound lets the
budget span a large discrete next-prefix jump while guaranteeing that
floating-point calibration convergence or inconsistent thresholds cannot
keep an otherwise stable retained view repainting forever.

Presentation timing has two consumers with deliberately different contracts.
The LoD controller receives the short-horizon cadence signal so it can react
within a few frames.  The faceplate receives a separate 750 ms elapsed-time
EMA, rebased at the start of a continuous navigation gesture and displayed to
one decimal place.  Smoothing the faceplate must not delay the controller, and
event-driven idle gaps must not be interpreted as renderer throughput.

Refinement priority is lexicographic:

1. selected occurrence;
2. highlighted occurrence;
3. missing minimum coverage;
4. estimated visible error reduction per added face;
5. projected pixel extent;
6. stable occurrence index.

The analytical PoP score is:

```text
error(c) = producer_error(c) * pixel_extent / asset_extent
benefit(c -> c+1) =
    affected_screen_span * max(0, error(c) - error(c+1))
value = benefit / added_faces
```

This is a constant-time hot-path estimate.  A cached builder-time surface or
silhouette delta may replace it only if measurements show better visual
ordering without per-frame triangle scans.

## Complexity obligations

| Operation | Required bound |
|---|---|
| adaptive interaction ceiling | `O(number of presentations)`, not occurrences |
| one compact scheduling wave | `O(wave_size)` projection plus bounded sorting |
| sparse active-cut update | `O(changed occurrences)` |
| draw submission | `O(visible occurrences + draw groups)` |
| completed-frame resource sample | `O(unique presentations)` once; subsequent policy/HUD queries `O(1)` |
| off-frustum draw cost | zero triangles and no per-asset upload |
| cache/provider work | bounded by distinct demanded assets, not occurrences |
| selection/style update | `O(changed paths/instances)`, no full hierarchy scan |
| stable compaction | background/quiet only, after coverage and every visual-settling witness; never on the input critical path |

Repeated-instance and distinct-asset scale are separate proof obligations.
The former stresses occurrence planning and instancing.  The latter stresses
source import, cache fan-out, PoP construction parallelism, result publication,
GPU upload, and memory pressure.  Passing one does not imply passing the other.

## Verification mapping

| Contract area | Required evidence |
|---|---|
| cold fallback and first mesh | empty-cache GUI timeline and screenshots |
| progressive autoview | exact cold/warm and System GL/OSMesa final-camera equivalence, automated cross-run center/size equality, and no provisional rollback |
| resident retarget | unit test proving no provider/cache work |
| resident suffix/trim | cache test proving suffix-only reads leave no reader prefix, realization test proving exact-capacity trim, and Lucy zoom memory telemetry |
| continuous zoom refinement | cold/warm Lucy held-gesture checkpoints proving resident growth and a richer effective cut before release on System GL and OSMesa |
| discrete prefix progress | unit test with an active cut below a richer resident prefix and an unaffordable direct jump |
| view working-set turnover | close multi-instance occurrence hashes and images across view directions |
| unique asset fan-out | thousands-of-distinct-mesh cold/warm stress, worker/cache telemetry, and perf |
| saturated residency | hard-cap cold/warm stress proving a quiet memory-limited terminal cut, bounded task count, and resubmission only after an admission revision edge |
| queued-result ownership | service test proving duplicate rejection before result drain and readmission after drain; warm stress proving bounded task fan-out |
| scene budget | held-motion and stable face/FPS telemetry |
| tail/silhouette quality | Generic Twin multi-angle image comparison |
| renderer packet semantics | independent exhaustive oracle over sparse ownership, hidden/proxy channel rules, retained ranges, and the full explicit-cut input domain |
| GL state | deep before/after state sentinel on every exercised System GL and OSMesa route, with apitrace for any failure |
| wire parity | shaded/wire active-cut and image matrix |
| selection/edit | hierarchy selection, erase/redraw, promotion/demotion, and picking tests |
| liveness | exhaustive scalar phase/event canonicalization; 512 seeded 96-event fake-clock/fake-service sequences; explicit checkpoint/failure/cancellation-pressure scenarios; and reports rejecting pending-without-witness or stable-with-affordable-next states |

The GUI runner is permitted to use generous wall-clock time for cold asset
construction, but it must not use that generosity to accept a terminal box,
an internally pending idle loop, or a visibly incomplete stable cut.
