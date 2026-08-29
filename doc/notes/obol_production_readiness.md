# Obol production-readiness matrix

Last reviewed: 2026-08-29

This is the release qualification plan for the BRL-CAD Obol drawing stack.  A
row is green only when it was run against the final binaries and records its
model, renderer, cache state, event script, viewport/DPR, report, images, and
performance evidence.  Historical measurements belong in
`libbobol_engineering_lessons.md`, not here.

## Current evidence snapshot

| Area | Current evidence | Release status |
|---|---|---|
| focused CTest gate | The dependency-complete 2026-08-29 build passes all 27 `bobol_headless` tests in 5.83 seconds, including renderer, LoD service/coordinator/update, cache, compact ownership, retained allocation, edit manipulator, host, API/symbol, and GED draw-sync contracts.  The policy-disable regression additionally proves that off/on/off transitions preserve the current presentation, retire every automatic owner, keep explicit manual generations usable, and cannot be rearmed by a renderer-style change or capacity-labelled repaint.  The focused eight-test LoD/draw-sync/control-trace gate also passes.  The service test proves that a coalesced asset producer may complete against its retained latest demand; the view-state test proves that a superseded result cannot create a terminal occurrence failure, while genuinely stale source/cache data retains its failure semantics.  Independent draw scope is exercised with a real local endpoint across draw/erase/redraw/zap.  The linked Obol CAD suite includes its 131,072-occurrence classifier-reservation regression.  The focused interaction gate also covers semantic selection, GED edit, qtcad primitive edit/preview/selection/measurement, qged event replay, polygon/sketch replay, and qged primitive-edit replay in single and quad layouts. | broader graphical production suite still required |
| qged retained interaction | The final 2026-08-29 binaries pass the focused 12-test qtcad/qged interaction gate and the dual-backend graphical matrices.  `/tmp/qged-selection-presentation-revision-20260829` passes System GL and OSMesa point add/remove, fractional-DPR rectangle placement/commit, selected styling, erase-selected, and redraw-selected.  `/tmp/qged-polygon-visual-presentation-revision-20260829` passes retained polygon styling, selection, movement, resize, zoom, cleanup, and cross-renderer placement; immediate selection/move deltas are 1,583/4,894 pixels on System GL and 1,586/4,921 on OSMesa.  `/tmp/qged-framebuffer-presentation-revision-20260829` passes in 37/37 seconds and combines progressive Generic Twin drawing, hierarchy selection, an ellipsoid edit manipulator, faceplate center marker, external raytrace framebuffer underlay/overlay/interlay, resize, rerender, and framebuffer disable.  Shared feature/polygon stores now publish content revisions independently of their coalesced controller render latch, so every actual shared mutation wakes the attached views while queries and equal publication remain passive. | focused cross-renderer selection/polygon/framebuffer/edit gate green; full control, physical-pointer, and model matrix remains |
| qged measurement | `qged_measure_ui_replay` and its quad-layout counterpart drive the actual palette and canvas through 2D and exact-hit 3D distance/angle gestures, degree/radian readback, cancellation, resize, and tool replacement.  They pass on OSMesa in the ordinary CTest gate and on final System OpenGL and OSMesa binaries in single and quad layouts; retained images are in `/tmp/qged-measure-{system,osmesa}-{single,quad}-final-20260829`.  The replay asserts line-layer point counts and field values, requires visible framebuffer deltas for both measurement modes, and requires cancellation to reproduce the exact baseline (zero changed pixels).  It exposed and fixed three production defects: no endpoint binding for right-click cancel, palette deactivation that failed to detach its semantic filter, and measurement lines hidden by shaded geometry.  Line-layer depth behavior is now explicit, with measurement guides depth-independent and other diagnostic overlays depth-tested by default. | dual-backend, single/quad, 2D/3D measurement gate green; fractional-DPR physical-device repetition remains |
| qged view settings | `qged_view_settings_ui_replay` and its quad-layout counterpart pass on final System OpenGL and OSMesa binaries in single and quad layouts; retained frames are in `/tmp/qged-view-settings-{system,osmesa}-{single,quad}-20260829`.  The replay drives the real settings palette, toolbar, GED commands, endpoint properties, retained Obol features, and framebuffer.  It proves bidirectional ADC/center-dot/grid/model-axes/scale/view-axes/parameter/FPS/framebuffer state, exact cutting-plane command readback, a visible shaded clipping delta, and pixel-exact restoration when clipping is disabled.  Stable replay IDs cover every exercised control.  System GL exposed a host race in which a fast renderer consumed the transient refresh latch before qged's post-command comparison; canvas checkpoints now compare libbv's monotonic frame revision as the durable semantic witness. | dual-backend, single/quad settings and cutting-plane control gate green; fractional-DPR and in-scene plane-affordance qualification remain |
| qged draw modes and resize | `/tmp/qged-resize-moss-hermetic-final-20260829` passes all 20 moss rows: modes 0, 1, 2, 4, and 5, managed LoD and full-detail policy, on System GL and OSMesa.  `/tmp/qged-resize-rook-m3-hermetic-final-20260829` passes the four corresponding mode-3 rows.  Each run resizes during realization and after settling, exercises minimize/maximize/fullscreen restore plus a resize storm and policy round trip, and terminates with an identical cross-policy camera, zero structural boxes, and no progressive work.  GUI replay now uses an application-owned temporary QSettings scope and never reads or writes the operator's geometry or plugin choices.  Top-level resize requests carry generations and an optional native stability barrier, so a delayed window-manager acknowledgement cannot change autoview or let an old retry resurrect a superseded size; the event-player regression injects that race directly. | baseline all-mode, dual-backend resize/policy matrix green; repeat on larger models and fractional-DPR hardware |
| control-plane models | The focused publication, host-work, convergence, admission, arbitration, canonical-pipeline, cold-preview, static-quality, renderer-preparation, interaction-session, deadline-ownership, capacity/handoff, timing-evidence, resident-growth, point-recovery, structural-frontier, and producer-demand models pass.  The canonical pipeline explored 2,358,764 generated / 1,095,220 distinct states; control refinement explored 4,194,302 / 2,097,152.  The new cross-boundary composition model explored 35,600 / 18,136 states to depth 15 and proves that admission, independent/capacity-owned growth, exact visibility presentation, ceiling reconciliation, capacity sampling, structural repair, point quality, and terminal publication share one compatible owner/terminal contract.  The current 32-budget/four-sample/four-candidate-per-goal capacity model explored 189,853 / 139,661 states to depth 44 and proves bounded candidate publication, one-sided settling after a safe result, direct cross-deadline evidence reuse, exact allocation/presentation ownership, and finite completed-presentation failure.  On 2026-08-29 `ObolPointTerminalEvidence` added the missing idempotence contract for a constrained point cut: an unchanged exact structural census cannot erase the only terminal witness after its finer-preload decision was consumed.  TLC explored 243 generated / 135 distinct states to depth 15.  The revised admission model explored 321 / 151 states to depth 6 and proves late safe-scene certification may atomically promote an existing PoP payload by marginal cost.  The revised capacity handoff explored 183,306 / 21,438 states to depth 20 and rejects a zero-work sample unless an exact census proves the scene is empty.  Corresponding coordinator, update-action, and dual-renderer Generic Twin regressions pass.  Runtime refinement checks sampled regression, duplicate plans, spontaneous reopen, unwitnessed presentation/constraints, and invalid readiness. | formal composition and sampled-runtime boundaries green; complete event/effect reducer coverage and per-transition records remain |
| Lucy OSMesa | The exact-current 2026-08-26 true-cold shaded lifecycle at `/tmp/qged-lucy-latest-demand-contract-20260826` passes in 76 seconds and separately records globally representative coverage and the first real CAD mesh.  Camera changes made while the immutable hierarchy is being built no longer discard its live pages or final result: the producer resolves page selection against the service-owned latest demand, and service validation uses that same demand.  Two consecutive certified-warm replays after splitting superseded work from stale source failures pass in 63/62 seconds at `/tmp/qged-lucy-superseded-fix-{a,b}-20260826`, each with a terminal ready, box-free mesh and zero occurrence failures.  The typed submission-owner replay at `/tmp/qged-submission-owner-lucy-20260826` passes in 82 seconds with 1.95M presented faces and 0.632-pixel maximum certified error.  After the exact source-delta owner replaced its independent active/source/plan fields, `/tmp/qged-submission-delta-lucy-20260826` also passes in 94 seconds with a box-free 1.265-pixel terminal view.  These lifecycles exercise continuous zoom/refinement, zoom-out recovery, rotation, lighting, hierarchy selection, subpath erase/redraw, camera ownership, and the strict HUD contract.  The 2026-08-28 exact-current warm wire lifecycle at `/tmp/qged-lucy-wire-capacity-transfer-revision-20260828` passes in 53 seconds after older frame-barrier precedence and capacity-transfer revision ownership were fixed; it completes all 191 events with 861,234 final mesh faces, zero boxes, and a clean external control trace.  Two canonical-handoff replays at `/tmp/qged-lucy-wire-canonical-handoff-{a,b}-20260828.json` first proved request-owned reconciliation budget and terminal-measurement retirement.  A later timing trace showed that ordinary planning could still ignore the retained terminal certificate and reselect the adjacent known-slow PoP population.  The terminal budget is now authoritative for ordinary admission, and two independent full replays at `/tmp/qged-lucy-wire-terminal-planner-clamp-{a,b}-20260828.json` pass all 191 events in 51.4/52.6 seconds.  Events 33 and 147 terminate ready instead of repeating the search; both final frames are box-free with no owner, obligation, pending calibration, or runtime-contract violation.  The final-binary true-cold/certified-warm wire pair at `/tmp/qged-production-lucy-wire-osmesa-20260828` passes in 51/58 seconds.  Both 191-event replays finish box-free, ownerless, and quality-floor-clean; cold presents 559,494 faces within the measured software budget and warm safely admits 1,284,926 faces. | cold/warm shaded and wire OSMesa interaction rows green |
| 50k scale | Cold shaded scale checks pass on System GL and OSMesa in 8.7/10.0 s to terminal CONSTRAINED.  A later exact-current run exposed interactive deadline recovery walking down one PoP ordinal per missed frame even when the measured cost ratio already proved that insufficient.  Recovery now operates in render-cost space during motion as well as at rest; the prior-pose deadline floor remains quiet-only.  The warm OSMesa shaded matrix at `/tmp/qged-50k-superseded-fix-20260826` passes in 78 seconds: held-motion recovery reaches its responsive ceiling directly, and the terminal view has 7,566 mesh payloads, 1.52M faces, zero boxes/failures, and no control owner.  After replacing the raw source-plan/census latches with typed owner-thread values, `/tmp/qged-submission-owner-50k-20260826` passes in 72 seconds and terminates ready with 6,633 progressive payloads, 1.32M faces, zero boxes/failures, and no pending work.  The subsequent exact source-delta ownership replay at `/tmp/qged-submission-delta-50k-20260826` passes in 95 seconds with 7,448 progressive payloads, 1.50M faces, 3.000-pixel certified error, and the same box-free terminal contract.  The exact-current warm OSMesa wire lifecycle remains green in 74 s.  The post-terminal-authority endpoint at `/tmp/qged-50k-terminal-planner-clamp-20260828.json` settles terminal and box-free in 7.0 seconds with 255 meshes, 49,745 subpixel occurrences, no owner/obligation/violation, and matching requested, active, and certified budgets.  Independent final-binary true-cold shaded OSMesa runs at `/tmp/qged-50k-constraint-witness-{final,repeat}-20260828` pass the complete interaction/camera matrix in 48.8/46.0 seconds.  Both select the same 6,935 mesh payloads and 43,065 subpixel occurrences, present 1.57--1.63M faces, have zero quality-floor misses/control violations, and explicitly record stable-budget, subpixel-aggregation, and static-deadline witnesses (`constraintEvidenceMask=13`) rather than relying on an inferred `performanceLimited` label.  The 2026-08-29 final-binary warm shaded OSMesa replay at `/tmp/qged-production-50k-osmesa-shaded-idempotent-20260829` passes the complete scripted interaction in 20 seconds after point-terminal census idempotence was enforced.  It ends box-free with 753 meshes, 49,282 aggregated occurrences, 363,752 faces, and a 126 ms software frame. | warm shaded/wire and repeated true-cold shaded OSMesa correctness, liveness, and interaction gates green; remaining System-GL cold/wire rows remain |
| 150k scale | The bounded cold crash gate passed under a 16 GiB address-space cap on System GL and OSMesa with all 150,001 leaves discovered and no terminal failures.  With the default scheduler (eight active realization workers on this host), the exact-current System-GL replay at `/tmp/qged_150k_timing_fix_report.json` reaches a stable terminal endpoint in 40.4 seconds with 78,341 mesh occurrences, zero visible boxes/failures, and no owner, obligation, background work, or runtime-contract violation.  The rebuilt typed-host OSMesa replay at `/tmp/qged_150k_typed_timing_osmesa_report_retry.json` terminates in 36.8 seconds with 4,097 meshes, 384,755 faces, 145,903 structural occurrences represented by 1,426 point draw records, and the same zero-box/control contract.  It is responsiveness constrained but has zero total/prominent quality-floor violations and zero visual-importance debt.  Threshold-stamped CAD timing prevents the prior 32/64-pixel balancing loop and stops faceplate-only frames from inflating CAD capacity.  The post-terminal-authority endpoint at `/tmp/qged-150k-terminal-planner-clamp-20260828.json` settles terminal and box-free in 29.3 seconds with 686 meshes, 149,344 subpixel occurrences, no owner/obligation/violation, and matching requested, active, and certified budgets.  A later trace found an unchanged structural-distribution seed reopening the already constrained point cut after its finer preload was rejected.  The idempotent reducer fix and formal contract are qualified by `/tmp/qged-production-150k-idempotent-seed-20260829`: the final-binary warm System-GL wire interaction passes in 27 seconds, remains box-free, and each constrained finer candidate is consumed once with no identical terminal replay.  This is liveness/performance evidence rather than final realistic visual-significance clearance. | current shaded System-GL/OSMesa liveness green; true-cold/wire rows and realistic prominent-object quality remain |
| System GL smoke | The current shaded Generic Twin cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 12/11 seconds after camera identity became an optional snapshot.  Both exact-view returns recall history and every terminal checkpoint is ready with zero structural boxes and occurrence failures; the final certified errors are 0.236/0.241 pixels.  The 2026-08-27 post-ownership full System GL interaction replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-system-report.json` passed in 9.7 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The post-terminal-authority replay at `/tmp/qged-generic-system-terminal-planner-clamp-20260828.json` passes in 9.5 seconds and ends ready with 709 meshes, 135,073 faces, 57,102 lines, zero boxes, and no owner/obligation/violation.  The prior wire evidence is `/tmp/qged-generic-wire-system-debug-20260824`; both terminal wire frames contain 709 mesh payloads, zero boxes/failures/pending work, and matching cold/warm camera state.  The exact-current Lucy shaded cold/warm interaction replay at `/tmp/qged-lucy-system-20260824` also passed: both final frames are exact/ready, box-free, and quality-floor-clean with one resident source payload, 7.41M displayed PoP faces, a 230 MB resident mesh set, and one quality-history recall.  The final-binary Lucy true-cold/certified-warm wire pair at `/tmp/qged-production-lucy-wire-system-20260828` passes in 39/16 seconds.  Both runs finish with the same 3,183,110-face cut, 9,549,330 wire lines, no boxes, pending work, control violation, or quality-floor miss, and matching camera state.  Hubble shaded/wire cold/warm pass after the overview lifecycle repair. | Lucy shaded/wire green; complete the remaining System-GL large-model matrix |
| OSMesa Generic Twin | The current shaded cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 15/14 seconds.  Both exact-view returns recall history and terminate ready with zero boxes/failures; the final certified errors are 0.621/0.676 pixels.  The 2026-08-27 post-ownership full replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-osmesa-report.json` passed in 14.3 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The final scene content and camera match the paired System-GL images.  The 2026-08-26 cross-renderer wire replay at `/tmp/qged-generic-wire-hud-final-20260826` passes cold and warm on System GL and OSMesa after the availability/HUD changes, with no invalid ready-label sample or terminal box.  The exact 2026-08-29 cold shaded cross-renderer replay at `/tmp/qged-production-generic-shaded-idempotent-20260829` passes in 11 seconds on System GL and 14 seconds on OSMesa; both end with all 709 meshes, approximately 135k faces, a one-pixel point threshold, and zero structural fallback. | focused regression green; continue the larger real-model matrix |
| spatial Lucy, System GL | Exact-current 2026-08-25 warm shaded retained-page replay passes.  It refines during continuous zoom, reaches its quiet pixel target, compacts on zoom-out, restores the prior view through quality history, and has no boxes or page-level subpixel proxies. | focused retained-page regression green; cold and wire rows remain |
| direct-mesh Generic Twin | The 2026-08-29 cold replays at `/tmp/qged-generic-promotion-{cold,osmesa}` pass on System GL and OSMesa.  Late safe-scene certification now promotes an already visible PoP payload atomically by marginal replacement cost.  Each terminal frame has all 709 BoT occurrences in direct full detail, 185,388 faces, zero progressive/proxy/structural-box payloads, and no pending work.  The overall extent preview remains a discovery-only presentation and is not counted as a semantic leaf. | focused admission regression green; keep discovery-preview latency distinct from terminal-mesh admission |
| concurrent cold cache | On 2026-08-28 two System-GL qged processes opened the same `Generic_Twin.g` and one initially empty `BU_DIR_CACHE` simultaneously.  Both exited successfully with terminal ready, box-free mesh presentations and no occurrence failure.  `mdb_stat` validated all three resulting LMDB environments (117 LoD, 83 draw, and 4 name-map entries), and an independent third process reopened the shared result, loaded from cache, and again converged terminal/ready without a cache write. | cross-process cache correctness green; large-asset duplicate-generation resource pressure remains to qualify |
| Hubble OSMesa | Exact-current 2026-08-25 shaded cold/warm lifecycles and camera contracts pass after the static-quality restoration change, with 1,804 warm shaded payloads and zero structural boxes.  The final-binary cold/warm wire pair at `/tmp/qged-production-hubble-wire-osmesa-20260828` passes in 16/16 seconds.  Both finish terminal, box-free, ownerless, and quality-floor-clean with 1,965 mesh payloads plus 553 subpixel occurrences.  The hierarchy replay retains one selected CAD instance across exact subpath erase/redraw, loads 2,530 tree items without a fallback scan, and applies its path update in at most 436 microseconds.  `/tmp/qged-resize-hubble-body-selection-final-20260829` adds four shaded managed/full-detail resize runs on System GL and OSMesa using the substantial `all.g/BODY` region rather than the old small probe.  Selection survives every window transition and the resize storm; the terminal selected frame is box-free and idle, clearing it changes 4,763--4,925 scene pixels, and all four final cameras match exactly.  The visually inspected frames show the same body panels receiving selected styling on both renderers. | shaded selection/resize and shaded/wire cold/warm focused gates green; complete wire resize and broader physical-pointer modes |
| NIST BREP | The final-binary cold/warm wire lifecycle at `/tmp/qged-nist-resident-wire-final-20260828` passes on System GL and OSMesa.  All four runs select cut 10 with 19,457 lines, refine to cut 11 with 28,892 lines on zoom-in, and return to cut 10/19,457 lines on zoom-out.  Each endpoint is terminal with one direct progressive occurrence, zero structural boxes, and zero LoD-service tasks; cold HUD state remains `Discovering model` until that occurrence arrives.  The focused source-mutation test replaces the resident progressive part through the sparse delta journal and proves its cut and aggregate metrics retire.  The paired shaded OSMesa regression at `/tmp/qged-nist-shaded-regression-20260828` remains box-free and view-responsive. | adaptive wire cold/warm cross-renderer gate green; shaded memory-reclamation and larger real-BREP qualification remain |
| multi-Lucy/xpush capacity | Exact-current 2026-08-26 warm initial-view runs pass on both renderers after capacity samples were ordered behind ceiling-free occurrence-plan handoff.  The pre-fix timing replay split 2/4 between pixel demand and a premature ordinary-deadline endpoint; after the applied-but-hidden candidate fix, 6/6 early-checkpoint OSMesa runs reached the identical approximately 984.1k-face, 0.457-pixel pixel-demand endpoint in 4.59--7.04 s with no boxes or control work.  Post-recovery-ownership OSMesa replays remain green: multi-Lucy reaches 985,820 faces before a performance-limited zoom endpoint, while xpush reaches 985,808 faces initially and 1.53M after zoom, all box-free and terminal.  System GL reached 3.81M faces at 0.228 pixels in 2.45 s.  The 2026-08-27 full 75-event OSMesa turnover/zoom/rotation replay at `/tmp/qged-multi-lucy-growth-handoff-20260827/formal-contract-report.json` completed in 29.2 s after completed-pass ownership and annotation lifetime were made atomic.  All waits settled terminal/ready with zero failures or control violations; the final view had eight meshes, 319k faces, zero boxes, and no remaining owner or obligation. | warm full-interaction regression green; cold eight-distinct-asset preparation, compaction, and wire rows remain |

The 2026-08-29 aggregate-proxy update keeps a genuinely tiny occurrence as one
point but preserves any footprint larger than five pixels as a batched box.
Cache-backed monolithic meshes use a validated PCA OBB when it is at least five
percent tighter by surface measure; absent, invalid, or insignificant metadata
uses the AABB.  Wire mode draws its 12 edges, shaded mode draws 12 lit
triangles, and shaded-with-edges draws both without creating per-occurrence
scene draw calls.
The same classification now applies to indirect-atlas admission pressure;
that path formerly emitted an unconditional point and discarded the source
bounds.  A forced one-MiB-atlas regression proves 128 screen-significant
pressure replacements produce 1,536 shaded triangles (and 1,536 edges when
requested), then collapse to exactly 128 points only after zooming below the
five-pixel ceiling.  Cache round-trip/validation and realization tests preserve
the OBB metadata, and the focused Obol CAD suite passes all 19 contracts on
System GL and OSMesa.  `/tmp/qged-obb-generic-20260829` additionally passes
the final-binary Generic Twin shaded cold/warm matrix on both renderers: all
four terminal frames have 709 direct full-detail meshes and zero AABB, OBB,
generic proxy, progressive, or structural-fallback payloads.  This proves the
OBB path does not defeat direct-to-mesh admission for a modest scene.  A
heterogeneous pressure-driven graphical comparison and group-level multi-page
aggregate ownership remain open qualification items.
The final Generic Twin resize/policy matrix at
`/tmp/qged-production-generic-resize-policy-disable-final-20260829` passes
auto and initially-off shaded cases on System GL in 6/10 seconds and OSMesa in
9/10 seconds.  Every endpoint has zero control owner, obligation, violation,
pending render/progressive work, and structural fallback boxes.  The
off/on/off checkpoints retain meshes rather than restarting at boxes, and
viewport/style repaints while disabled remain presentation-only.
The final point-bracket policy combines its timing-ratio correction with the
known safe/unsafe bounds, so it may skip redundant thresholds but can never
cross the proven coarse side.  `ObolCadTimingEvidence` exhaustively checks all
bounded jumps across the full 64-cut production domain (7,726,131 generated /
1,222,796 distinct states, depth 12), and the focused coordinator regression
covers the measured correction numerically.  The current 150k shaded OSMesa
replay at `/tmp/qged-production-150k-osmesa-threshold-acceleration-20260829`
passes in 49 seconds: its former twelve-frame `32 -> ... -> 64` staircase is
now `32 -> 55.3 -> 64`, and it terminates with 2,075 meshes, 147,964 aggregated
occurrences in 274 renderer records, zero structural fallbacks, and no control
owner or obligation.  Current follow-up replays pass for 50k shaded OSMesa
(`/tmp/qged-production-50k-osmesa-threshold-acceleration-20260829`, 20 s),
Lucy shaded OSMesa
(`/tmp/qged-production-lucy-osmesa-threshold-acceleration-20260829`, 45 s,
340,082 faces, 1.265-pixel maximum error, no quality-floor violation), and
cold shaded Generic Twin on both renderers
(`/tmp/qged-production-generic-shaded-threshold-acceleration-20260829`,
12/14 s).  These runs validate the complete client path in addition to the
renderer-level primitive/accounting checks.

The 2026-08-28 exact-current shaded smokes cover both renderers without
replacing the broader interaction rows above.  Generic Twin reaches a terminal
box-free endpoint in 1.85 seconds on System GL and 4.34 seconds on the rebuilt
OSMesa host,
with 673 meshes.  Lucy reaches a terminal box-free endpoint in 1.75 seconds on
System GL and 2.70 seconds on OSMesa; both present 559,494 faces at a
0.928-pixel maximum certified error.  The 50k fixture reaches a terminal
box-free endpoint in 11.1 seconds on System GL with 27,968 mesh occurrences,
23,006 point occurrences, and 2.27M faces.  OSMesa reaches a constrained
box-free endpoint in 14.1 seconds with 1,094 meshes, 1,431 point draw records,
and 378k faces.  All report zero runtime-contract violations.  The next scale
gates are true-cold large-asset latency and realistic visual-significance
qualification rather than another object-count-specific controller change.

The final-binary Generic Twin wire matrix at
`/tmp/qged-generic-wire-cross-renderer-fix-20260828` passes cold and warm on
System GL in 10/10 seconds and OSMesa in 15/19 seconds.  All four endpoints
contain 709 mesh occurrences, about 135k scene faces, zero structural
fallbacks, zero control work or violations, and exact camera-contract matches.
This row also regresses the completed selective source-delta to scene-wide
handoff: the former OSMesa warm loop repeatedly visited 14 satisfied entries
and could not certify its 709-occurrence allocation.

The 150k terminal view may validly be performance-limited.  It is not visually
qualified until prominent-object quality floors are demonstrated on realistic
large models, not merely synthetic coverage fixtures.

## Common acceptance criteria

Every graphical row must show:

- correct final camera, extent, draw mode, colors/materials, lighting, and
  retained scene composition;
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

Keep one shared copy of reusable large inputs and canonical warm caches.  Keep
the event script, final report, representative checkpoint images, and perf/APNG
only when they explain a result.  Remove redundant transient screenshots,
trace logs, and duplicate generated models after recording the needed evidence.
