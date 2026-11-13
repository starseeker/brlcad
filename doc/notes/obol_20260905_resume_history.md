# Obol drawing: September 5 resume history

Archived 2026-09-05 when the simplification roadmap replaced chronological
resume instructions.  This is a frozen evidence snapshot, **not current
instructions, a specification, or a backlog**.  Statements such as “latest”
and “next” below refer to their original context and may be superseded even
within the snapshot.  Do not append new session status here.

Use the [current handoff](obol_session_handoff.md),
[simplification guide](obol_simplification_guide.md), and
[active debt](libbobol_active_debt.md) to resume.  Release qualification and
exact observations remain in [production readiness](obol_production_readiness.md).
These older notes preserve artifact paths, diagnoses and operational details
while keeping them out of the normal reading path.

## Previous handoff

Reviewed 2026-09-05.  **The stack is not production-qualified.** The current
focus is spatial-page visual quality and progressive convergence, followed by
scale and full GUI/editing qualification.  This is a resume snapshot, not a
second backlog: [active debt](libbobol_active_debt.md) owns remaining work.

### Checkout and build facts

| Repository | Reviewed state |
|---|---|
| `/home/cyapp/starseeker/brlcad` | Branch `brlobol6`, `2d73328a92`, with resumed diagnostic, replay-ownership, source-admission and spatial-presentation edits; Ninja build in `.build`, Debug, sanitizers off |
| `../obol` (root `obol` symlink target) | Branch `main`, `a544ae38`, with uncommitted constant-face lighting and indirect atlas-ceiling fixes plus image regressions; `build` is Release, System GL and software enabled |
| `../obol/external/osmesa` | Real submodule, now `bfa90de`; includes repeated-normal lighting optimization |
| `../osmesa` | Same `bfa90de`, with existing `.build` and example tests enabled |

Prior implementation changes are committed; the September 5 resume's
diagnostic/validator, controller/source and documentation edits are uncommitted.  Preserve the
user's untracked `brlobol_rework.txt`: it
proposes a replacement control kernel but is **not an adopted specification**.
Some proposals (formal hierarchy, compositions, refinement tests, typed
ownership) already exist; its service/executor maintainability critique still
deserves evaluation.  Do not launch a wholesale rewrite from that note alone.

**Shared binaries reconciled September 5.** Obol/OSMesa were rebuilt, their
libraries/headers/package explicitly installed into BRL-CAD's `.build`, and
BRL-CAD consumers rebuilt.  Build-tree and installed SHA-256 identities match:
OSMesa `097aac27...`, Obol now `f281b59b...`; qged resolves both from `.build/lib`.
Build/install logs and full hashes are under
`/tmp/qged-obol-resume-20260905/provenance`; the latest atlas fix is under
`/tmp/obol-atlas-ceiling-20260905`, preserving `5400bd13...` under `baseline`.
Earlier lighting logs are in `/tmp/obol-fixed-face-lighting-20260905`, with the
preceding `e02411a8...` library under `baseline`.  The original `dd3ac7b1...`
library remains in `/tmp/obol-flat-face-lighting-20260905/baseline`.  Obol's cached
default install prefix is `/home/cyapp/bext/build/install`, **not** this
BRL-CAD build; do not blindly install to the default.  Changing shared C++
layouts requires rebuilding consumers, including test executables.

### What is implemented and must survive

- GED owns draw/erase frontiers, selection, edit sessions and renderer-neutral
  intent.  libBObol owns realization, immutable PoP/cache/service data, view
  demand and resource policy.  Obol owns retained geometry and view-local
  preparation/rendering.  Spatial pages remain one logical CAD occurrence.
- Six-domain evidence stamps, typed owners, bounded cursors, transactional
  presentation and complete transition tracing are implemented.  TLA+ has
  canonical, composition and component tiers with executable refinement gates.
  It does not establish numeric fidelity, real-time latency or actual memory
  bounds merely by exploring a finite control model.
- Spatial-page authored/flat/smooth normal support and worker-side normal
  generation exist.  Readiness now matches requested normal policy/crease,
  not just the presence of an older drawable mesh.  No GUI-thread expansion
  of millions of corner normals should be reintroduced.
- Fixed-function shaded and wire cut uploads now report exact finite
  preparation work.  Normal-policy publication and service-admission epochs
  have concrete readiness/trace witnesses.  Replaying the same upload must
  not manufacture progress or new capacity evidence.
- Bare-root primitive dispatch was repaired.  The retained headless corpus
  covered 64 databases/16,917 roots in shaded/wire and LoD on/off; this is not
  graphical proof for every root.  HUD workload estimation and generated
  primitive editing are also implemented foundations, not tasks to start over.

### Current evidence and unresolved failures

The [readiness matrix](obol_production_readiness.md), especially its September
5 audit, contains exact observations and qualifications.  Key distinctions:

1. **Observation semantics are repaired; Lucy still fails qualification.**
   Latest paired warm replays are under `/tmp/qged-obol-resume-20260905`.
   Software finishes 217 events in 123.478 s; all 2,505 control transitions
   pass.  Discrete and smooth close zoom still miss the prominent floor
   (normalized errors 3.705 and 3.581).  Zoom-out passes this time; the final
   frame is ready with 2,101,208 faces.  System GL finishes 217 events in
   105.435 s, but smooth close zoom misses the floor (8.709), and its control
   checker rejects an A/B/A at trace indices 714--716, event 139, around
   `lod-static-overscan-resident-growth`.  Frame completion alone is not a
   finite-work witness.  The follow-up fix shares the frame consumer's
   submission/capacity/publication readiness gate with the replay requester.
   `/tmp/qged-obol-static-replay-20260905` completes all 217 events in 131.530 s
   and passes all 1,184 transitions without relaxing the checker.  Its smooth
   close floor still fails (error 5.717).  The earlier desktop timeout remains
   separate and unresolved.

   Disabled compact-registry diagnostics now serialize nulls.  The standalone
   `lod_lucy_contract.jq` names each failed condition, tolerates a startup point
   only before mesh coverage/readiness, and retains the first-useful-image and
   cold/warm mesh deadlines.  Coverage publication and exact mesh presentation
   are separate checkpoints.  Historical reports remain unchanged; replaying
   their new Lucy filter still rejects the genuine quality misses.
2. **Normals are improved, not fully qualified.** Latest software flat/authored/
   smooth/return checkpoints use cut 24 and 2,101,208 faces, at approximately
   326/317/369/308 ms.  The matching System GL run drops from cut 22 to 21 for
   smooth normals, then restores 22; it is not a same-cut comparison.
   Repeat controlled same-cut image/timing comparisons and cross-page
   continuity, zoom-in/zoom-out, cold/warm and both backends.  A local OSMesa
   microbenchmark improves repeated-normal multi-light work roughly 15--17%
   with identical pixels; no end-to-end speedup or performant software GLSL
   conclusion is established.
3. **A desktop frame-delivery stall remains unresolved.**
   `/tmp/qged-lucy-fixed-preparation-20260904` stopped at event 28 awaiting
   `render-preparation-replay`, no workers, no subsequent frame.  The isolated
   X-server replay completes.  Investigate canvas exposure, Qt request delivery
   and endpoint claiming; do not assume another capacity-loop defect.
4. **Scale results are historical.** Generic Twin's latest shaded System GL/
   OSMesa cold/warm matrix passes.  Buddha and the broad feature matrices have
   useful prior evidence.  Current shared changes still need Lucy/multi-Lucy,
   Hubble, heterogeneous 50k/150k and real-vehicle requalification.  Old stress
   caches were format 23 and deleted; current mesh cache format is 24.
5. **Bare giant admission now uses bounded serialized coverage.** The original
   `draw -m1 lucy.s` timeout remains retained.  Streamed v5 BoTs now reserve the
   source coordinator's finite allowance and enter the existing serialized
   coverage census; detail stays behind its own memory governor.  Incomplete
   coverage cannot escape into the unrestricted serial primitive fallback.
   `/tmp/qged-obol-bare-coverage-20260905` passes shaded/wire on both backends
   with the shared warm mesh cache, one logical occurrence, exact ready mesh,
   no boxes/proxies/floor misses, and passing traces, in 4.1--5.7 s.  Each uses
   1,284,926 active mesh faces; wire presents 3,854,778 lines.  Cold giant and
   multiple-root memory/concurrency qualification remain.
6. **Resident page bytes must not suppress presentation repair.** A new
   shortened software run timed out in planning at event 25, repeatedly
   advancing capacity revisions while its selected allocation stayed
   unapplied.  Retain `/tmp/qged-obol-close-perf-20260905` and
   `/tmp/qged-obol-planning-cycle-20260905`; the latter's strict trace passes
   despite the script timeout, so revision changes alone are insufficient
   liveness evidence.  A deterministic regression exposed a retained fast
   path which skipped missing page geometry after a successful same-cut
   demand update.  Submission now requires prepared presentation at the
   admitted cut before skipping the bounded worker repair.  Shaded, wire and
   hidden-line regressions pass with prefetch disabled and zero extra draw
   allowance.  Warm 111-event repeats at
   `/tmp/qged-obol-spatial-repair-20260905/{osmesa,system}` finish in
   56.349/42.194 s and pass 1,472/466 control transitions.  Both are exactly
   ready but still miss the close quality floor (errors 3.581/8.709).  This
   closes the deterministic preparation gap, not every historical timeout.
   The final-binary full software replay at
   `/tmp/qged-obol-spatial-repair-full-20260905` completes 217 events in
   111.814 s and passes all 1,567 control transitions.  Its final 2,101,208-face
   view is ready and ownerless, but the matrix fails the smooth close
   prominent-floor check (error 3.581, events 106/108/110).


7. **Constant face lighting now skips redundant software interpolation.**
   The profiled original library rejected protected cut 29 (cost 4,333,525)
   after a 419.763 ms draw with no preparation, against the unchanged 400 ms
   deadline.  Obol now uses flat shading only for normal-free batches with
   directional lights and an infinite viewer, restoring the caller's shade
   model.  Explicit normals, positional lights and local-viewer lighting
   retain their original path.  A nonempty exact-image regression explicitly
   requires the batch atlas and covers five lighting/state cases.  Seven
   software tests and six System GL tests pass; the software-only cut-cache
   test is correctly skipped on System GL.  Raw Obol executables require
   `OBOL_TEST_RENDER_BACKEND=swrast` to actually select software in this dual
   build.  Nine BRL-CAD focused CTests pass (3.01 s).

   `/tmp/qged-obol-flat-lighting-perf-20260905` completes 111 events in
   62.113 s and passes 840 transitions; close cut 29 presents 3,445,998 faces,
   error 2.859, no floor miss, at 372.390 ms.  The matched original profile is
   `/tmp/qged-obol-close-render-perf-20260905` (54.281 s, close cut 28/error
   3.581).  These close snapshots actually use tier 1 (retained cuts), so
   the batch-only change does not establish a direct close-up speedup.
   Overall replay time is not a same-work speedup measurement.
   The full software row `/tmp/qged-obol-flat-lighting-full-20260905` completes
   217 events in 66.049 s, passes 1,241 transitions and all Lucy quality checks,
   including close error 1.429.  Its matrix **fails the separate HUD check**:
   event 8 has a visible 13-by-7-pixel fill (91 counted pixels), while the
   validator selects it using a newer progress estimate and requires 100.
   The new `lod_hud_contract.sh` repairs that observation boundary by
   selecting a sufficiently large recorded fill.  Its 100-pixel floor remains,
   and gray-only tracks, insufficient colored pixels, missing images and
   invalid geometry fail standalone regressions.
   `/tmp/qged-obol-flat-lighting-system-20260905` completes 217 events in
   94.470 s and passes 860 transitions, but still **fails** close quality
   (error 4.354, cut 27, 1,806,720 faces).  Large preparation slices and other
   liveness/platform/scale debts remain; this is not production qualification.


8. **The HUD contract is repaired; further image/renderer gates remain.**
   The clean follow-up software matrix
   `/tmp/qged-obol-hud-contract-final-20260905` completes 217 events in
   102.957 s and passes all 1,998 transitions and Lucy quality checks.  The HUD
   helper independently verifies 1,365 fill pixels.  The full matrix still
   **fails** the normal-policy comparison: flat/authored report active cut 25,
   authored presents 2,101,208 faces instead of the active 3,183,110, and
   smooth/returned settle at cut 24.  Preserve this evidence; do not waive
   same-cut image validation.  The earlier `hud-contract-full` run overlapped
   a CMake rebuild and was terminated; it is not timing evidence.

   System GL under private Xvfb uses Mesa 25.2.8 llvmpipe (LLVM 20.1.2), not a
   hardware GPU.  `/tmp/qged-obol-system-close-perf-20260905` finishes 111 events
   in 37.904 s and passes 472 transitions, but rejects cut 28 after a 439.327 ms
   preparation-free draw.  The lighting optimization now also covers the
   retained-cut executor, using one shared eligibility predicate.  Both route
   image tests pass against the preserved baseline and current library;
   eight software/seven System GL focused tests and ten BRL-CAD CTests pass.
   `/tmp/qged-obol-fixed-face-system-perf-20260905` finishes 111 events in
   38.790 s and passes 482 transitions, but still rejects cut 28 at 442.349 ms
   and settles cut 27/error 4.354.  **No meaningful System GL speedup is proved.**
   The renderer's 512 MiB batch bound, actual executor selection, large
   preparation latency, and fixed-cut work accounting deserve further
   inspection before changing allocation policy or deadlines.

9. **Fractional preparation retries now retain their exact population.**
   `/tmp/qged-obol-fixed-face-full-osmesa-20260905` times out at event 105
   after 60 seconds of readiness waiting.  Its 1,809 transitions pass the
   strict checker despite 79 preparation targets during that event and no
   pending workers.  Frames are arriving, so this is distinct from the
   desktop lost-frame stall.  The interruption handler preserved the integer
   ceiling but dropped its fractional page subset, completing a cheaper cut
   which could request another fractional candidate.

   The rebuilt controller now retains the fraction on an unchanged-ceiling
   preparation retry.  `/tmp/qged-obol-fractional-retry-20260905` **passes the
   full software matrix**: 217 events in 98.545 s, all 1,960 transitions,
   quality, normal-policy images/coherence, and HUD visibility (1,365 pixels).
   Event 105 actually interrupts a fractional candidate and retains the same
   target through retry, completed preparation, timing replay, and owner
   release.  Close error is 2.859 with no prominent-floor miss; final readiness
   has no owner or host work.  Ten focused CTests pass (2.98 s), and the existing
   preparation model passes TLC (526 distinct states); neither a new isolated
   C++ fractional-handler fixture nor full-suite formal qualification is
   claimed.  The unchanged-deadline System GL quality failure and prior
   normal-policy/planning/desktop reproducers remain separate debts.  See the
   readiness note for the exact transition sequence and binary/test logs.
   The same final binaries at
   `/tmp/qged-obol-fractional-retry-system-20260905` finish all 217 events in
   102.925 s and pass all 916 transitions, with final readiness and no owner or
   host work.  That matrix still **fails** close quality (cut 27/error 4.354).
   Its later normal-image/HUD gates do not run after the quality rejection;
   smooth normals also record a different cut from flat/authored/returned.

10. **Static capacity search now starts in its actual deadline domain.**
    An opt-in llvmpipe programmable replay had completed a richer static frame
    before falling to cut 21/134,950 faces/error 17.417.  Its new search used
    STEADY with both duration fields set to 400 ms, applying an earlier
    258,437-unit preferred-cadence ceiling.  `capacitySearchKey()` now shares
    explicit STATIC entry and the unchanged ordinary preferred duration
    between completed-pass and interrupted-frame paths.  The regression also
    verifies that a real static miss still constrains the candidate.

    `/tmp/qged-obol-static-domain-glsl-close-20260905` completes 111 events in
    30.832 s, passes 430 transitions, and reaches cut 29/3,445,998 faces/error
    2.859 with no floor miss at 376.402 ms.  This uses the existing
    `OBOL_CAD_SOFTWARE_GLSL=1` option; renderer defaults are unchanged.
    The full default System GL row
    `/tmp/qged-obol-static-search-domain-system-20260905` completes 217 events
    in 91.338 s and passes 832 transitions, but still fails close quality at
    cut 27/error 4.354.  Ten focused CTests pass (3.01 s).  The capacity model
    now includes `staticTrial` entry and passes TLC (357,336 distinct states);
    its individual baseline/catalog entry is updated.  Static quality also
    passes TLC (128 states), with all 44 model/config pairs linted.  Keep
    broader renderer/scale qualification separate from this controller fix.
    The matching default OSMesa matrix
    `/tmp/qged-obol-static-search-domain-osmesa-20260905` passes all 217 events
    in 103.681 s, all 1,481 transitions, Lucy quality, normal-policy images,
    and HUD visibility.  Close cut 29 has error 2.859 and no floor miss;
    normal checkpoints all retain cut 24/2,101,208 faces.  Final readiness has
    no owner or host work.
    The full opt-in programmable row
    `/tmp/qged-obol-static-search-domain-glsl-full-20260905` finishes 217 events
    in 76.798 s and passes 914 transitions, settled quality, normal-policy
    images, and HUD visibility.  Its matrix still fails in-gesture smooth-zoom
    refinement: rendered cut bounds 25/23/23/23 do not advance from the
    starting cut, and resident bytes remain unchanged through those captures.
    See `gesture-failure.json`; preserve this next reproducer and the existing
    gate.  The richer quiet close view does not qualify continuous input.

11. **Motion recovery preserves the actual failed-cut witness.** Cost-based
    deadline recovery can jump several cuts; its cheap fallback previously
    prohibited all intervening quality probes.  Both miss paths now supply the
    attempted cut separately.  The coordinator regression permits completed
    intermediate cuts but excludes the actual failed cut; ten focused CTests
    pass (3.39 s).  The instrumented 111-event replay
    `/tmp/qged-obol-motion-miss-witness-glsl-close-20260905` exercises those
    probes and passes 426 transitions in 23.984 s, reaching close cut 29/error
    2.859.  No deadlines, memory limits, renderer defaults or validators change.

    The full default OSMesa row
    `/tmp/qged-obol-motion-miss-witness-osmesa-20260905` **passes** all 217 events
    in 69.244 s, all 1,270 transitions, quality, normal-policy images, smooth
    zoom and HUD (1,358 pixels).  Close cut 29 has error 1.429; normal snapshots
    all retain cut 23/1,284,926 faces.  Default System GL at
    `/tmp/qged-obol-motion-miss-witness-system-20260905` finishes 217 events in
    92.909 s and passes 778 transitions, but still fails close quality at
    cut 27/error 4.354.  It additionally reports no presented mesh at event 32's
    first-mesh checkpoint; preserve this frame-boundary observation.  Later
    image/HUD gates do not run after those failures.

    The full opt-in programmable row
    `/tmp/qged-obol-motion-miss-witness-glsl-full-20260905` finishes 217 events
    in 73.107 s and passes 900 transitions, settled quality, normals and HUD.
    In-gesture cut 27 now exceeds starting cut 26, with resident bytes growing
    from 262,094,432 to 279,277,676 and a 39.479 ms frame.  The matrix still
    fails its atlas-generation reuse/ordinary first-upload condition: those
    counters do not advance.  Both captured frames use indirect tier 6, and
    the atlas reuse counter increments on generation changes.  Audit actual
    incremental presentation before changing the gate; retain the named
    `gesture-failure.json` and `gesture-conditions.json`.  All three final
    views are ready with no owner or host work.  See readiness for exact
    timings, renderer and source observations, and remaining qualification.
    The final cleanup rebuild and ten CTests pass again (9.17 s); the
    `motion-miss-witness-final-{build,tests}` logs are in `/tmp`.  No GUI,
    profiler, TLC or build process remains active.

12. **Indirect ceiling changes now admit their missing atlas prefix.**
    The generation-counter investigation found a real renderer defect:
    occurrence-cut changes could upload a suffix, while renderer-ceiling
    changes only clamped to existing GPU data.  The new regression fails on
    the prior library (128 triangles instead of 38,400).  Both paths now share
    bounded admission and retain exact-preparation fallback on relocation.
    Occurrence and ceiling tests also verify nonempty rich pixels and exact
    rich/coarse/rich restoration without another geometry upload.

    Ten focused System GL and six software Obol tests pass; the final BRL-CAD
    rebuild and eleven CTests pass (3.00 s).  Logs and the prior library are in
    `/tmp/obol-atlas-ceiling-20260905`.  Installed/build Obol is `f281b59b...`.
    `lod_atlas_reuse.jq` now recognizes actual same-generation suffix uploads,
    with 13 adversarial cases.  It accepts the post-upload zoom report and
    still rejects the previous no-upload report.  The first-mesh wait also
    now requires exact presented geometry after delivering queued frames;
    it could previously return after an abort because CPU adoption sufficed.

    `/tmp/qged-obol-atlas-ceiling-glsl-full-20260905` finishes 217 events in
    74.630 s with 909 transitions, settled quality, normals and HUD passing.
    The revised zoom filter passes; the original generation-only witness had
    rejected 24 MB of actual suffix uploads.  The final all-changes opt-in row
    `/tmp/qged-obol-atlas-ceiling-glsl-final-20260905` finishes 217 events in
    69.874 s with 860 passing transitions and close cut 29/error 2.859.
    Its matrix still **fails quiet zoom-out**: cut 20/337,610 faces/error 4.274,
    41.700 ms, scene budget 10,270,264 against active cost 1,067,000, no renderer
    ceiling and no GPU pressure.  The mesh does carry a current memory-denial
    witness (admission revision 1), with resident cut 26 behind displayed cut
    20.  Final readiness has no owner or host work; the zoom-out compaction
    plan is absent as well.  The revised input refinement branch passes.
    This is a resident-denial/allocation and compaction reproducer; do not
    confuse it with the repaired GPU upload or relax the prominent floor.
    Final default System GL at
    `/tmp/qged-obol-atlas-ceiling-system-final-20260905` completes 217 events in
    100.532 s and passes 929 transitions.  First-mesh presentation now passes;
    close cut 27/error 4.354 still fails.  Final OSMesa at
    `/tmp/qged-obol-atlas-ceiling-osmesa-final-20260905` completes 217 events in
    104.708 s and passes 1,460 transitions, settled quality, same-cut normals
    and HUD (1,414 pixels), but fails continuous rendered-cut advancement:
    starting cut 25 versus captured bounds 25/23/25/25 despite resident/cache
    growth.  Its zoom-out compaction does pass.  Both final views are ready
    with no owner/host work.  Retain their distinct `gesture-conditions.json`
    evidence; no GUI, build, profiler or TLC job remains active.

13. **Memory denial no longer hides already prepared resident detail.**
    The failed opt-in zoom-out above has a current memory-denial witness,
    resident cut 26 and active cut 20.  The allocator incorrectly used the
    active cut as the entire denied availability domain.  It now permits
    richer immutable prepared prefixes, checking the demanded channel and
    every populated spatial page, and the retained application passes the
    provider-denial guard without scheduling new source work.  Unavailable
    and unprepared suffixes remain constrained; fixed bounded coverage and
    all memory/render/image limits are preserved.

    The expanded update-action regression fails before the fix and passes
    all three availability/preparation cases afterward.  Twelve CTests pass
    (8.23 s), including service, allocation oracle and coordinator.  Logs are
    `/tmp/obol-resident-denial-{before,after,final-tests}-20260905.log`, with
    sequential build logs under the same prefix.  No formal control model
    changed or new TLC result is claimed.  Current libBObol is `ce8d39fa...`;
    qged remains `2b28af79...` and Obol `f281b59b...`.  Full hashes are in
    `/tmp/obol-resident-denial-20260905/hashes.txt`.

    Both full warm shaded Lucy matrix rows pass every gate:
    `/tmp/qged-obol-resident-denial-glsl-20260905` (217 events, 78.677 s,
    928 transitions) and `/tmp/qged-obol-resident-denial-osmesa-20260905`
    (217 events, 69.911 s, 1,263 transitions).  Opt-in System GL closes at
    renderer cut 29/error 2.859, zooms out to cut 25/error 1.339, and compacts
    327,791,888 to 275,887,428 bytes.  OSMesa closes at cut 29/error 1.429,
    zooms out to cut 24/error 0.814, and compacts 222,923,888 to 95,123,820
    bytes.  Normal comparisons keep a common cut within each row (24 and
    23 respectively).  Final owners/host flags are zero.  The earlier default
    System GL close-floor failure and software cut-25-start input failure
    remain separate qualification debt.  Several diagnostic replays passed
    before this fix too; the deterministic regression is the causal evidence,
    not timing differences between adaptive graphical runs.

    The user has asked whether the implementation process and architecture
    should be reconsidered.  Address that assessment before adopting another
    redesign; no replacement scheduler or rewrite specification is approved.
    Preserve the implemented immutable data plane and existing contracts.

Fresh follow-up checks: seven focused CTests pass (2.62 s), and all four Generic
Twin wire cold/warm System GL/OSMesa matrix rows pass.  The follow-up small
corpus passes 156 shaded/wire checks over 78 roots.  All 44 model/config pairs
pass lint; focused TLC passes terminal convergence, terminal quality ordering
and static quality (1,187/10,232/128 distinct states).  The full canonical run
was stopped before completion; no new full-suite TLC, clean rebuild or platform
qualification is claimed.  Historical software Lucy rows pass; the latest
row retains a continuous-input failure, and repeated cross-regime production
qualification remains open.
The latest submission fix passes nine focused CTests in 2.75 s; its existing
retained-allocation presentation model passes TLC (57 distinct states).

### Recommended next work, in order

1. Preserve the repaired replay ownership and bare-root admission boundaries;
   extend their evidence to cold giant and multiple-root execution.  The strict
   trace checker and import memory bound remain unchanged.
2. Fix actual quality/latency issues at their owning numeric or renderer
   boundary, with same-cut perf/images and a cross-regime regression per fix.
   Start from `/tmp/qged-obol-close-budget-20260905`: the shortened replay
   passes control but misses the close floor.  Its stdout trace records a
   protected allocation followed by richer static work and a cheaper
   reconciliation; the readiness note records exact costs.  Inspection
   confirms one common allocation cut per logical payload: pages have
   independent residency, not independent presentation allocations.  Do not
   change terminal routing based on layer count alone.  Preserve the new
   same-cut page-preparation regression and audit the failed planning trace's
   repeated capacity revisions against real progress witnesses.
   Preserve the constant-face lighting optimization and exact-image test.
   Preserve the repaired HUD visibility contract.  Resolve the new software
   normal-policy comparison failure and continue System GL close-floor work
   from the actual retained-cut executor and latest paired replay.
   Preserve the fractional preparation retry identity and its exercised
   event-105 regression; audit quiet candidate changes as liveness debt,
   rather than accepting a new target signature as progress by itself.
   Preserve explicit STATIC search entry and its independent-deadline-ceiling
   regression.  The opt-in programmable path now provides a close-floor pass;
   qualify its full images and workload costs before changing defaults.
   Preserve the motion failed-cut witness and intermediate-probe regression.
   Preserve indirect ceiling-prefix admission and its image/count regression,
   the suffix-upload witness, and the exact first-mesh barrier.  Start the
   remaining opt-in audit from the latest cheap, coarse quiet zoom-out with
   substantial declared budget headroom; default System GL close quality
   remains a separate numeric/performance debt.
   Close the desktop lost-frame reproducer and measure cancellation latency
   inside large uploads.  Change models only if ownership/liveness changes.
3. Repeat small, giant-single, multiple-giant, and many-distinct-asset cases on
   System GL and OSMesa; cold then certified warm, zoom/rotate/translate,
   selection/erase/redraw, resize/DPR and memory turnover.  Use bounded subset
   LoD-off controls at dangerous scales, not an unbounded 150k full-detail draw.
4. Resume responsibility extraction and operational hardening.  Major units
   remain large: database source 18.3k lines, GED draw 12.9k, view controller
   8.3k, service 7.5k.  Extract independently compiled owners, not `.inc`
   fragments or forwarding wrappers.  Preserve sparse/dense storage, bounded
   publication, immutable sharing, and no second 150k scene vector.
5. Complete the broader production backlog: large-BREP typed/bounded provider,
   real vehicle visual-significance tests, shared-cache cold concurrency,
   Windows corruption/zoom and full shared-stack sanitizers, then qged/MGED/
   gsh/Archer/rtwizard feature qualification.  Original Big Boy BREP remains
   separate from the partial BoT conversion; fail-closed partial-tessellation
   rejection is not a complete BREP solution.

Do not lose the editing workstream: descriptor classification already exists;
runtime legality/rejection/readback still needs comprehensive auditing.  Keep
richer indexed manipulators, specialized sketch arcs/NURBs/profile editing,
MGED `sed`/`oed` lifecycle parity, polygon mouse creation/move/resize/booleans,
measurement modes, tree/selection scaling, and CLI/widget/manipulator agreement
explicit.  [qged_editing.md](qged_editing.md) has the detailed acceptance cases.

### Reproduction and retained assets

Run one heavy GUI/profiler at a time, with private display/settings isolation.
Reuse one model/cache copy, collect perf with matching instrumentation, and
verify processes exit.  Deep diagnostic JSON collection changes timings.

```sh
ctest --test-dir .build --output-on-failure -j2 \
  -R '^(libBObol_(lod_update_action|lod_coordinator|lod_progress_estimator)|qged_lod_control_trace_contract)$'

jq -e -f src/qged/tests/lod_control_trace.jq \
  /tmp/qged-lucy-lighting-reuse-perf-20260904/cases/lucy-osmesa-shaded-swapdefault-warm/report.json

# After rebuilding; choose one reusable output directory, preserving key evidence.
xvfb-run -a -s '-screen 0 1280x1024x24' \
  ./src/qged/tests/qged_gui_matrix.sh --build-dir ./.build \
  --artifact-dir /tmp/qged-lucy-resume \
  --cases lucy --backends osmesa --modes shaded --swap-intervals default \
  --warm-cache /tmp/qged-lucy-control-contract-cold-20260904/caches/lucy-osmesa-shaded-swapdefault/cache \
  --no-baseline --timeout 180
```

The historical `/home/cyapp/brlcad/.build/bin/qged` main baseline is absent.
Reestablish/identify a baseline before claiming new main comparisons.
`src/qged/tests/README.md` explains the lifecycle, visual-quality, resize and
interaction runners; `qged_lod_quality_matrix.sh` supports managed-only and
memory-limited qualification.  A successful event script is not a passing
validator, and a passing control trace is not a passing image.

Reusable inputs exist: `.build/lucy.g` (688 MB), `Happy_Buddha.g` (26 MB),
`Generic_Twin.g`, `unique_mesh_stress.g` (120 MB), 50k (1.20 GB), 150k (3.61 GB),
`stanford_local.g` (6.11 GB), and partial `bigboy.g` (569 MB).  Inspect Stanford
roots for multi-Lucy/xpush before recreating them.  Hubble is at
`/home/cyapp/models/NASA/Hubble/Hubble_Space_Telescope.g`; Havoc and NIST are in
`.build/share/db`.  `/tmp/obol-shared-models` contains other shared fixtures.
The backup `.build` under `/media/cyapp/Backup 2023/Unversioned/brlobol5` is an
additional source if needed, not a reason to duplicate all inputs.

Keep the shared Lucy cache above (~743 MiB), latest normal/perf and desktop
stall reports/images, `.build/testing-artifacts/osmesa-lighting`, and crash/
counterexample records.  `/tmp/obol-current-cache-matrix` is now only 39 MiB of
historical evidence, **not usable warm stress caches**.  Earlier cleanup freed
~40 GiB.  No active qged/debugger/perf/TLC process remained at review.  Prune
redundant captures only after retaining a compact result and reproduction;
other sessions' artifacts and shared assets are not disposable.

### Reading map and guardrails

Start here, then [active debt](libbobol_active_debt.md),
[architecture](libbobol_architecture.md), and the
[pipeline contract](libbobol_progressive_pipeline_contract.md).
[TLA README](tla/README.md), [risk coverage](tla/RISK_COVERAGE.md) and
[conformance audit](libbobol_tla_conformance.md) define the formal boundary;
`tla/models.json` and `tla/baselines/tlc-2.19.json` are the catalog/result
authority.  Use `/home/cyapp/tla+/tla2tools.jar` through
`misc/CMake/RunTLA.cmake`, with generated state outside the source tree.
[Engineering lessons](libbobol_engineering_lessons.md) preserves resolved
failure mechanisms; historical `brl_obol_*` notes are not competing designs.

Classify each new failure as control-contract, numeric-policy, data-plane or
observation before changing code.  Preserve exact identities and progress
witnesses; never fix a stall by relabeling work or declaring a pending normal
layer ready.  Pixel-exact means view-sufficient, not all source triangles.
Zoom starts from retained data; pose-only changes restore proven affordable
detail.  Geometry loading, resident retention and currently presented detail
are separate quantities.  Performance pressure may justify a witnessed
constraint, but a clean control trace cannot excuse missing prominent shapes.

## Previous active-debt resume priorities

Start with `obol_session_handoff.md` for the checkout, artifact, and command
snapshot.  This document remains the detailed remaining-work authority; the
handoff is not a second evolving backlog.  Historical green rows below do not
qualify the September 4 spatial-normal/renderer changes at every scale.

The September 5 resume reconciled shared binaries and repaired disabled
registry/initial-cue observation semantics.  Exact evidence is in the readiness
matrix; these are no longer tasks to repeat.  Current work is:

1. **Finish spatial-page quality and renderer qualification.** Authored,
   flat, and smooth normal layers and their readiness matching are implemented,
   not missing features.  Repeat same-cut images/timing, cross-page continuity,
   cold/warm close zoom, rotation recovery, and zoom-out memory reclamation on
   both backends.  Diagnose the separate desktop requested-frame timeout;
   isolated-display success does not retire it.  Fixed-cut uploads now report
   finite preparation, but cancellation latency within one large upload remains
   unqualified.  The OSMesa microbenchmark win is not proof of a performant
   software GLSL path.  Obol's additional constant-face interpolation
   optimization now passes exact-image regressions and the full software
   Lucy quality/control checks, but the paired System GL close floor still
   fails (error 4.354).  Preserve this software optimization while extending
   cross-regime evidence; differing adaptive workloads do not establish an
   overall speedup.  The HUD observation mismatch is repaired: the validator
   selects a sufficiently large recorded fill and retains its 100-pixel
   visibility floor, with hidden-fill and malformed-input regressions.  The
   subsequent full software replay passes quality/control but fails the normal
   comparison: authored cut 25 presents a cut-24 face population, followed by
   smooth/returned cut 24.  Establish controlled same-cut evidence and audit
   terminal allocation versus actual presentation; do not waive the comparison.
   The close snapshots use the retained-cut renderer, so a batch-only profile
   win is not evidence of a direct close-up speedup.  Preserve the shared
   constant-lighting predicate and separate batch/cut image regressions.
   Exact artifacts are in the readiness matrix.  Preserve the retained-page
   preparation guard and shaded/wire/hidden-line regressions.  The 111-event
   dual-backend zoom repeats pass control but still miss the close quality
   floor.  Audit repeated capacity revisions in the retained failed planning
   reproducer: `/tmp/qged-obol-planning-cycle-20260905` passes the trace filter
   despite its readiness timeout.  A changed revision is not by itself a
   finite-progress witness, and the deterministic repair does not qualify
   every historical stall.  The first full retained-cut software replay additionally
   timed out at event 105 in static fractional preparation.  Its trace passes
   while preparation signatures keep changing; the host is delivering aborted
   frames, so this is distinct from the desktop lost-frame hypothesis.  The
   interruption handler had dropped the fractional selection while preserving
   its integer ceiling.  The correction now preserves both coordinates, and
   the full software matrix passes all 217 events, including an interrupted
   fractional candidate which retains its target through retry and acceptance.
   Its normal-policy/image gates pass too, but one passing adaptive population
   does not retire the earlier differing-cut failure.  Preserve this exercised
   regression while strengthening trace witnesses for repeated candidate
   changes.  A fresh target is not itself progress in a quiet view.
   The opt-in System GL programmable path exposed another repaired boundary:
   a static-quality search started under STEADY with equal duration fields,
   then reused an older cadence ceiling and discarded its rich frame.  The
   shared coordinator key now selects STATIC explicitly, with deterministic
   independent-ceiling coverage and a passing expanded capacity model.  Its
   shortened graphical reproducer now meets the close floor at cut 29; the
   default System GL full replay still fails that floor at cut 27.  Keep
   programmable rendering opt-in: its full replay now passes settled quality,
   normal-policy images and HUD, but fails the continuous-zoom refinement gate.
   That reproducer led to a numeric-policy correction: a cheap recovery hint
   no longer rejects unmeasured intermediate cuts.  The actually failed cut
   remains excluded, with completed-frame regressions and an exercised probe
   trace.  The atlas follow-up found and repaired a renderer gap: raising only
   the renderer ceiling did not attempt to upload a missing GPU prefix.  Both
   cut paths now share bounded admission, with failing-before/passing-after
   image/count regressions.  The validator accepts measured same-generation
   suffix uploads and still rejects the old no-upload reproducer.  Preserve
   those proofs.  The latest opt-in final replay nevertheless fails quiet
   zoom-out: cut 20/337,610 faces/error 4.274 is declared ready at 41.700 ms,
   with a 10,270,264-unit scene budget, no renderer ceiling, and no GPU memory
   pressure; its zoom-out compaction plan is also absent.  The payload has a
   current memory-denial witness and resident cut 26 behind displayed cut 20.
   The follow-up allocation regression now admits already prepared resident
   detail while retaining the denial of unavailable/unprepared suffixes.
   The following full opt-in System GL and default OSMesa rows pass all
   matrix gates, including quality, input, normals, HUD and turnover.  These
   adaptive runs do not establish schedule-independent behavior.  Preserve
   the original failure and new deterministic regression while extending
   qualification; the render deadline remains unchanged.
   The final software row passes settled quality, normals and HUD, but fails
   in-gesture rendered-cut advancement: starting cut 25 versus captured bounds
   25/23/25/25, despite real cache/resident growth.  It does compact correctly
   on zoom-out.  Retain that distinct numeric-policy/input reproducer too.
   The first-mesh wait now checks completed geometry after delivering queued
   work; adoption could previously pass the wait before a frame aborted.
   Exact artifacts and cross-backend rows are in the readiness note.  Qualify
   this input path and cross-scale behavior before changing defaults; a fast,
   coarse frame is not evidence of a better renderer.
2. **Extend repaired replay/source boundaries to cold and concurrent roots.**
   The static replay requester now shares its consumer's producer-readiness
   gate; the follow-up System GL Lucy run passes all 1,184 strict transitions.
   Streamed bare v5 BoTs now use the existing bounded serialized census and
   cannot fall through to an unrestricted primitive import.  Bare Lucy passes
   warm shaded/wire probes on both backends in 4.1--5.7 s, with one logical
   occurrence and exact ready geometry.  Retain those regressions while
   qualifying cold giant, multiple roots, cancellation and memory contention.
   The smooth-close numeric quality miss and desktop frame-delivery timeout
   remain open; passing replay ownership does not retire either.
3. **Requalify both scale extremes before more control changes.** Repeat
   Generic Twin, Buddha, Hubble, Lucy/multi-Lucy, and heterogeneous 50k/150k
   lifecycles using the same binaries.  The old large-scene caches were format
   23 and were removed; current format 24 requires new cold preparation before
   certified warm runs.  Retain one cache per reusable fixture.
4. **Resume structural and editing work against that baseline.** Keep the
   extraction, large-BREP provider, platform/sanitizer, visual-significance,
   and GUI/editing gates below.  `brlobol_rework.txt` is an untracked proposal
   for a replacement control core, not an adopted specification.  Evaluate its
   service/ownership critique against today's implemented reducers and formal
   compositions before choosing another migration.  Do not recreate the
   already implemented model hierarchy or introduce a second scheduler.
