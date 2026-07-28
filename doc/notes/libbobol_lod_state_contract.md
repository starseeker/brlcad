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
   interaction frame target.
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
active level.  A presentation references assets and occurrences; it does not
own BRL-CAD source truth or choose LoD.

For a view `v`, occurrence `o`, and asset `a`:

```text
visible[v,o]          boolean from the renderer-contract frustum test
pixel_extent[v,o]     projected screen extent in pixels
demand[v,o]           clamp(ceil(log2(pixel_extent / pixel_error)), 0, 15)
active[v,o]           currently presented PoP level, or NONE
resident[a]           richest retained array prefix
drawable[a,l]         resident arrays contain all entries required by l
faces[a,l]            cumulative face prefix for l
scene_budget[v]       calibrated aggregate submitted-face allowance
render_ceiling[v]     temporary O(1) interaction-only maximum level, or NONE
```

In an orthographic view, equal projected extents have equal LoD value
regardless of depth.  In a perspective view, distance matters only through
projection.  No separate depth weighting is allowed.

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

`INTERACTIVE`
: A finite scene budget and an O(1) presentation ceiling bound the next frame.
  The exact per-occurrence allocator catches up in bounded waves.  Finer
  resident arrays are not evicted.

`SETTLING`
: The interaction cut remains until one coarse frame has actually been
  presented and the quiet debounce expires.

`STABLE`
: Refinement proceeds one drawable prefix step per presented frame.  Candidate
  steps are ordered by user emphasis and estimated screen-error reduction per
  added face.  Stable cuts do not regress within an unchanged view epoch.

`COMPACTING`
: After a longer quiet interval, aggregate consumer demand may trim resident
  CPU/GPU tails with headroom.  Active cuts remain drawable throughout.

## Transition obligations

Every pending flag has exactly one progress witness:

| Pending condition | Required witness | Completion |
|---|---|---|
| source realization | provider task or source callback | fallback or source data published |
| asset suffix | queued/in-flight/cache task | result applied or definitive failure |
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

For a retained asset, refinement is incremental:

```text
next(v,o) = minimum level l such that
            active[v,o] < l <= demand[v,o] and drawable[asset(o),l]
```

The allocator must try `next`, not jump directly to `resident`.  If `next`
does not fit, it may grow the probed allowance toward calibrated capacity or
declare the cut budget-limited.  It must not retry the same unaffordable jump
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
3. Shaded, wireframe, bounds, and picking use the same effective active level,
   including the temporary interaction ceiling.
4. An occurrence binding references the current source, view, policy, and
   occurrence revisions.  Stale results cannot replace current data.
5. Stable refinement is monotonic for an unchanged view/policy epoch.
6. All visible mesh-backed occurrences receive minimum coverage before any
   occurrence receives a second optional prefix increment.
7. Off-frustum occurrences issue no draw calls and consume no scene face
   budget.  Their shared asset may remain resident.
8. Selection/highlight changes priority and style, not asset identity.
9. A source edit invalidates affected content; a camera edit does not.
10. Direct GL rendering restores GL state and invalidates Coin state caches
    before returning.

## Liveness properties

Using `[]` for “always” and `<>` for “eventually,” the implementation must
satisfy:

```text
[] (visible_mesh_backed && minimum_resident
    -> <> (active >= minimum || definitive_failure))

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
interactive demand:  adaptive pixel_error >= 4
scene face cost:     sum(faces[asset(o), effective_level(v,o)])
effective level:     min(active[v,o], render_ceiling[v]) when ceiling is set
```

The calibrated stable and interactive capacities are separate.  A quiet
readback or a motion presentation interval cannot contaminate the other
capacity.  Budget reductions affect a future view/motion epoch; they do not
make an unchanged quiet frame visibly regress.

Refinement priority is lexicographic:

1. selected occurrence;
2. highlighted occurrence;
3. missing minimum coverage;
4. estimated visible error reduction per added face;
5. projected pixel extent;
6. stable occurrence index.

The analytical PoP score is:

```text
error(l) ~= pixel_extent / 2^l
benefit(l -> l+1) =
    affected_screen_span * max(0, error(l) - error(l+1))
value = benefit / added_faces
```

This is a constant-time hot-path estimate.  A cached builder-time surface or
silhouette delta may replace it only if measurements show better visual
ordering without per-frame triangle scans.

## Complexity obligations

| Operation | Required bound |
|---|---|
| interaction emergency coarsen | `O(number of presentations)`, not occurrences |
| one compact scheduling wave | `O(wave_size)` projection plus bounded sorting |
| sparse active-cut update | `O(changed occurrences)` |
| draw submission | `O(visible occurrences + draw groups)` |
| off-frustum draw cost | zero triangles and no per-asset upload |
| cache/provider work | bounded by distinct demanded assets, not occurrences |
| selection/style update | `O(changed paths/instances)`, no full hierarchy scan |
| stable compaction | background/quiet only; never on the input critical path |

Repeated-instance and distinct-asset scale are separate proof obligations.
The former stresses occurrence planning and instancing.  The latter stresses
source import, cache fan-out, PoP construction parallelism, result publication,
GPU upload, and memory pressure.  Passing one does not imply passing the other.

## Verification mapping

| Contract area | Required evidence |
|---|---|
| cold fallback and first mesh | empty-cache GUI timeline and screenshots |
| resident retarget | unit test proving no provider/cache work |
| discrete prefix progress | unit test with active level below a richer resident prefix and an unaffordable direct jump |
| view working-set turnover | close multi-instance occurrence hashes and images across view directions |
| unique asset fan-out | thousands-of-distinct-mesh cold/warm stress, worker/cache telemetry, and perf |
| scene budget | held-motion and stable face/FPS telemetry |
| tail/silhouette quality | Generic Twin multi-angle image comparison |
| GL state | System GL apitrace plus corruption/noise screenshots |
| wire parity | shaded/wire active-level and image matrix |
| selection/edit | hierarchy selection, erase/redraw, promotion/demotion, and picking tests |
| liveness | report rejects pending-without-witness and stable-with-affordable-next states |

The GUI runner is permitted to use generous wall-clock time for cold asset
construction, but it must not use that generosity to accept a terminal box,
an internally pending idle loop, or a visibly incomplete stable cut.
