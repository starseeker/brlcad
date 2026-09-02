# Obol production-readiness matrix

Last reviewed: 2026-09-01

This is the release qualification plan for the BRL-CAD Obol drawing stack.  A
row is green only when it was run against the final binaries and records its
model, renderer, cache state, event script, viewport/DPR, report, images, and
performance evidence.  Historical measurements belong in
`libbobol_engineering_lessons.md`, not here.

## Current evidence snapshot

| Area | Current evidence | Release status |
|---|---|---|
| focused CTest gate | After a complete current-tree Ninja build on 2026-09-01, all 28 `bobol_headless` tests pass in 9.0 seconds, including renderer, LoD service/coordinator/update, cache, compact ownership, retained allocation, edit manipulator, host, API/symbol, GED draw-sync, and view-command contracts.  The policy-disable regression additionally proves that off/on/off transitions preserve the current presentation, retire every automatic owner, keep explicit manual generations usable, and cannot be rearmed by a renderer-style change or capacity-labelled repaint.  The service test proves that a coalesced asset producer may complete against its retained latest demand; the view-state test proves that a superseded result cannot create a terminal occurrence failure, while genuinely stale source/cache data retains its failure semantics.  Independent draw scope is exercised with a real local endpoint across draw/erase/redraw/zap.  The linked Obol CAD suite includes its 131,072-occurrence classifier-reservation regression.  The current OSMesa/offscreen interaction gate passes all 13 qged event, measurement, settings, polygon, polygon/sketch, primitive-edit, and framebuffer-host tests in 14.7 seconds.  After explicitly relinking the two model-test executables against the current controller ABI, all 19 qtcad Pinewood, Havoc, M35, NIST, and Generic Twin real-model/progressive tests pass in 46 seconds. | broader graphical production suite still required |
| shared clients | Current gsh, MGED, Archer/TkObol, and rtwizard binaries were explicitly rebuilt against the current shared stack.  A fresh private-X-server gate passes all 15 thread-affinity, Tk widget/attach, Archer, gsh draw/ERT/rt-routing/progressive-LoD, MGED host/draw/framebuffer/ERT/progressive-LoD/edit-restore, and rtwizard smoke tests in 135 seconds. | shared-client smoke and core framebuffer/edit routes green; comprehensive interactive client qualification remains |
| qged retained interaction | The final 2026-08-29 binaries pass the focused 12-test qtcad/qged interaction gate and the dual-backend graphical matrices.  `/tmp/qged-selection-presentation-revision-20260829` passes System GL and OSMesa point add/remove, fractional-DPR rectangle placement/commit, selected styling, erase-selected, and redraw-selected.  The 2026-08-31 follow-up at `/tmp/qged-selection-ui-delivery-final-20260831` additionally proves every observed selection frame drains to idle without changing the LoD policy or requesting renderer-capacity evidence.  Single-object ARB/ELL/sketch manipulators and the navigation gizmo now use the typed presentation-only endpoint API, and direct GL explicitly presents the sole exact-frame barrier when Qt would otherwise coalesce it indefinitely.  Full Hubble hierarchy selection/erase/redraw passes on System GL and OSMesa in 9/26 seconds at `/tmp/qged-selection-stability-hubble-{system-final2,final}-20260831`, with the selection frame idle before the next checkpoint and no LoD service work.  A later intermittent quiet-handoff stall was traced to a level-triggered demand refresh whose bounded cursor had been retired by a stronger presentation owner.  The guarded cursor restart now passes a fresh System GL/OSMesa Hubble lifecycle in 11/14 seconds plus two independent OSMesa repetitions at `/tmp/qged-selection-demand-restart*-20260831`; every selected observation preserves the policy revision and returns to IDLE/ready with zero obligations.  After terminal capacity certificates were made to retire predecessor handoff debt atomically, fresh point-selection replays at `/tmp/qged-hubble-point-selection-post-handoff-{system,osmesa}-report.json` pass in 5.2/9.9 seconds: selection creates one exact-presentation obligation, preserves the view/policy revisions, and the following observation is ready with zero control or service debt.  `/tmp/qged-polygon-visual-presentation-revision-20260829` passes retained polygon styling, selection, movement, resize, zoom, cleanup, and cross-renderer placement; immediate selection/move deltas are 1,583/4,894 pixels on System GL and 1,586/4,921 on OSMesa.  `/tmp/qged-framebuffer-presentation-revision-20260829` passes in 37/37 seconds and combines progressive Generic Twin drawing, hierarchy selection, an ellipsoid edit manipulator, faceplate center marker, external raytrace framebuffer underlay/overlay/interlay, resize, rerender, and framebuffer disable.  The current-binary reruns at `/tmp/qged-selection-ui-contract-final-20260901`, `/tmp/qged-polygon-visual-contract-final4-20260901`, and `/tmp/qged-framebuffer-contract-final2-20260901` pass on both renderers; the matching System-GL single/quad widget replays are under `/tmp/qged-system-ui-final-20260901`.  Shared feature/polygon stores now publish content revisions independently of their coalesced controller render latch, so every actual shared mutation wakes the attached views while queries and equal publication remain passive. | focused cross-renderer selection/polygon/framebuffer/edit gate green; full control, physical-pointer, and model matrix remains |
| qged measurement | `qged_measure_ui_replay` and its quad-layout counterpart drive the actual palette and canvas through 2D and exact-hit 3D distance/angle gestures, degree/radian readback, cancellation, resize, and tool replacement.  They pass on OSMesa in the ordinary CTest gate and on final System OpenGL and OSMesa binaries in single and quad layouts; the current System-GL evidence is under `/tmp/qged-system-ui-final-20260901/measure-*`.  The replay asserts line-layer point counts and field values, requires visible framebuffer deltas for both measurement modes, and requires cancellation to reproduce the exact baseline (zero changed pixels).  It exposed and fixed three production defects: no endpoint binding for right-click cancel, palette deactivation that failed to detach its semantic filter, and measurement lines hidden by shaded geometry.  Line-layer depth behavior is now explicit, with measurement guides depth-independent and other diagnostic overlays depth-tested by default. | dual-backend, single/quad, 2D/3D measurement gate green; fractional-DPR physical-device repetition remains |
| qged view settings | `qged_view_settings_ui_replay` and its quad-layout counterpart pass on final System OpenGL and OSMesa binaries in single and quad layouts; the current System-GL evidence is under `/tmp/qged-system-ui-final-20260901/settings-*`.  The replay drives the real settings palette, toolbar, GED commands, endpoint properties, retained Obol features, and framebuffer.  It proves bidirectional ADC/center-dot/grid/model-axes/scale/view-axes/parameter/FPS/framebuffer state, exact cutting-plane command readback, a visible shaded clipping delta, and pixel-exact restoration when clipping is disabled.  Stable replay IDs cover every exercised control.  System GL exposed a host race in which a fast renderer consumed the transient refresh latch before qged's post-command comparison; canvas checkpoints now compare libbv's monotonic frame revision as the durable semantic witness. | dual-backend, single/quad settings and cutting-plane control gate green; fractional-DPR and in-scene plane-affordance qualification remain |
| qged draw modes and resize | `/tmp/qged-resize-moss-all-modes-contract-final2-20260901` passes all 24 moss workflows: modes 0--5, managed LoD and full-detail policy, on System GL and OSMesa.  The apparent mode-5 camera mismatch was a harness race: a delayed native X11 configure changed the window aspect after the resize wait but before draw/autoview.  The runner now requires one second of initial geometry stability and proves the semantic draw batch executes at the requested geometry.  The isolated replacement at `/tmp/qged-resize-moss-m5-framing-final-20260901` passes all four mode-5 rows with identical cross-policy camera size and height.  The current mode-3 rook replacement at `/tmp/qged-resize-rook-m3-contract-final-20260901` passes its four rows in 5--13 seconds.  Each run resizes during realization and after settling, exercises minimize/maximize/fullscreen restore plus a resize storm and policy round trip, and terminates with zero structural boxes and no progressive work.  GUI replay uses an application-owned temporary QSettings scope and never reads or writes the operator's geometry or plugin choices.  Top-level resize requests carry generations and an optional native stability barrier, so a delayed window-manager acknowledgement cannot change autoview or let an old retry resurrect a superseded size; the event-player regression injects that race directly. | baseline all-mode, dual-backend resize/policy matrix green; repeat on larger models and fractional-DPR hardware |
| control-plane models | The focused publication, host-work, convergence, admission, arbitration, canonical-pipeline, cold-preview, static-quality, renderer-preparation, interaction-session, deadline-ownership, capacity/handoff, timing-evidence, resident-growth, point-recovery, structural-frontier, and producer-demand models pass.  The canonical pipeline explored 2,358,764 generated / 1,095,220 distinct states.  The current 29-fact control refinement explored 475,078 / 237,568 states and distinguishes exclusive interrupted replay from ordinary exact-presentation debt, while the lifecycle composition explored 1,092,377 / 227,787 states to depth 32 and proves recovery after two bounded demand-cursor interruptions plus isolation of semantic-only exact frames from capacity recovery.  The cross-boundary composition model explored 35,600 / 18,136 states to depth 15 and proves that admission, independent/capacity-owned growth, exact visibility presentation, ceiling reconciliation, capacity sampling, structural repair, point quality, and terminal publication share one compatible owner/terminal contract.  The terminal-convergence composition explores 1,520 / 1,187 states to depth 81 and requires an exact visibility census and its exact framebuffer classification before reallocation.  The bounded capacity search includes the protected minimum in its immutable numeric domain, freezes capacity-owned allowance inputs, and exports an exact allocation/presentation/three-sample progress rank; its current exploration count is recorded in `libbobol_formal_models.md`.  `ObolCadFrameComposition` now treats camera and effective renderer controls as one exact-frame revision (4,921 / 1,568 states, depth 17).  `ObolHostWork` includes unsubmitted source revisions, timing-driven work opened by frame completion, and provider retirement while an independent exact frame remains; it explores 78,500 / 12,902 states to depth 17.  The 2026-09-01 `ObolSubmissionPass` explores 90 / 45 states to depth 11 and proves that scene-wide demand debt survives a selective source pass until a complete successor consumes it or a newer semantic revision supersedes it.  `ObolLodConvergence` explores 742 / 423 states to depth 25 with the matching quiet-quality owner rule.  On 2026-08-29 `ObolPointTerminalEvidence` added the missing idempotence contract for a constrained point cut: an unchanged exact structural census cannot erase the only terminal witness after its finer-preload decision was consumed.  TLC explored 243 generated / 135 distinct states to depth 15.  The revised admission model explored 321 / 151 states to depth 6 and proves late safe-scene certification may atomically promote an existing PoP payload by marginal cost.  The 2026-09-01 capacity handoff explores 159,114 / 17,598 states to depth 19, requires a resident retry to resolve as drawable or explicitly constrained, rejects a zero-work sample unless an exact census proves the scene is empty, and requires a terminal certificate to retire older reconciliation debt atomically.  Corresponding coordinator, update-action, and dual-renderer Generic Twin regressions pass.  Runtime refinement checks sampled regression, duplicate plans, spontaneous reopen, unwitnessed presentation/constraints, invalid readiness, and six-domain revision monotonicity. | formal composition and sampled-runtime boundaries green; complete event/effect reducer coverage and per-transition records remain |
| Lucy OSMesa | The exact-current 2026-08-26 true-cold shaded lifecycle at `/tmp/qged-lucy-latest-demand-contract-20260826` passes in 76 seconds and separately records globally representative coverage and the first real CAD mesh.  Camera changes made while the immutable hierarchy is being built no longer discard its live pages or final result: the producer resolves page selection against the service-owned latest demand, and service validation uses that same demand.  Two consecutive certified-warm replays after splitting superseded work from stale source failures pass in 63/62 seconds at `/tmp/qged-lucy-superseded-fix-{a,b}-20260826`, each with a terminal ready, box-free mesh and zero occurrence failures.  The typed submission-owner replay at `/tmp/qged-submission-owner-lucy-20260826` passes in 82 seconds with 1.95M presented faces and 0.632-pixel maximum certified error.  After the exact source-delta owner replaced its independent active/source/plan fields, `/tmp/qged-submission-delta-lucy-20260826` also passes in 94 seconds with a box-free 1.265-pixel terminal view.  These lifecycles exercise continuous zoom/refinement, zoom-out recovery, rotation, lighting, hierarchy selection, subpath erase/redraw, camera ownership, and the strict HUD contract.  The 2026-08-28 exact-current warm wire lifecycle at `/tmp/qged-lucy-wire-capacity-transfer-revision-20260828` passes in 53 seconds after older frame-barrier precedence and capacity-transfer revision ownership were fixed; it completes all 191 events with 861,234 final mesh faces, zero boxes, and a clean external control trace.  Two canonical-handoff replays at `/tmp/qged-lucy-wire-canonical-handoff-{a,b}-20260828.json` first proved request-owned reconciliation budget and terminal-measurement retirement.  A later timing trace showed that ordinary planning could still ignore the retained terminal certificate and reselect the adjacent known-slow PoP population.  The terminal budget is now authoritative for ordinary admission, and two independent full replays at `/tmp/qged-lucy-wire-terminal-planner-clamp-{a,b}-20260828.json` pass all 191 events in 51.4/52.6 seconds.  Events 33 and 147 terminate ready instead of repeating the search; both final frames are box-free with no owner, obligation, pending calibration, or runtime-contract violation.  The final-binary true-cold/certified-warm wire pair at `/tmp/qged-production-lucy-wire-osmesa-20260828` passes in 51/58 seconds.  Both 191-event replays finish box-free, ownerless, and quality-floor-clean; cold presents 559,494 faces within the measured software budget and warm safely admits 1,284,926 faces. | cold/warm shaded and wire OSMesa interaction rows green |
| 50k scale | Cold shaded scale checks pass on System GL and OSMesa in 8.7/10.0 s to terminal CONSTRAINED.  A later exact-current run exposed interactive deadline recovery walking down one PoP ordinal per missed frame even when the measured cost ratio already proved that insufficient.  Recovery now operates in render-cost space during motion as well as at rest; the prior-pose deadline floor remains quiet-only.  The warm OSMesa shaded matrix at `/tmp/qged-50k-superseded-fix-20260826` passes in 78 seconds: held-motion recovery reaches its responsive ceiling directly, and the terminal view has 7,566 mesh payloads, 1.52M faces, zero boxes/failures, and no control owner.  After replacing the raw source-plan/census latches with typed owner-thread values, `/tmp/qged-submission-owner-50k-20260826` passes in 72 seconds and terminates ready with 6,633 progressive payloads, 1.32M faces, zero boxes/failures, and no pending work.  The subsequent exact source-delta ownership replay at `/tmp/qged-submission-delta-50k-20260826` passes in 95 seconds with 7,448 progressive payloads, 1.50M faces, 3.000-pixel certified error, and the same box-free terminal contract.  The exact-current warm OSMesa wire lifecycle remains green in 74 s.  The post-terminal-authority endpoint at `/tmp/qged-50k-terminal-planner-clamp-20260828.json` settles terminal and box-free in 7.0 seconds with 255 meshes, 49,745 subpixel occurrences, no owner/obligation/violation, and matching requested, active, and certified budgets.  Independent final-binary true-cold shaded OSMesa runs at `/tmp/qged-50k-constraint-witness-{final,repeat}-20260828` pass the complete interaction/camera matrix in 48.8/46.0 seconds.  Both select the same 6,935 mesh payloads and 43,065 subpixel occurrences, present 1.57--1.63M faces, have zero quality-floor misses/control violations, and explicitly record stable-budget, subpixel-aggregation, and static-deadline witnesses (`constraintEvidenceMask=13`) rather than relying on an inferred `performanceLimited` label.  The 2026-08-29 final-binary warm shaded OSMesa replay at `/tmp/qged-production-50k-osmesa-shaded-idempotent-20260829` passes the complete scripted interaction in 20 seconds after point-terminal census idempotence was enforced.  It ends box-free with 753 meshes, 49,282 aggregated occurrences, 363,752 faces, and a 126 ms software frame. | warm shaded/wire and repeated true-cold shaded OSMesa correctness, liveness, and interaction gates green; remaining System-GL cold/wire rows remain |
| 150k scale | The bounded cold crash gate passed under a 16 GiB address-space cap on System GL and OSMesa with all 150,001 leaves discovered and no terminal failures.  With the default scheduler (eight active realization workers on this host), the exact-current System-GL replay at `/tmp/qged_150k_timing_fix_report.json` reaches a stable terminal endpoint in 40.4 seconds with 78,341 mesh occurrences, zero visible boxes/failures, and no owner, obligation, background work, or runtime-contract violation.  The rebuilt typed-host OSMesa replay at `/tmp/qged_150k_typed_timing_osmesa_report_retry.json` terminates in 36.8 seconds with 4,097 meshes, 384,755 faces, 145,903 structural occurrences represented by 1,426 point draw records, and the same zero-box/control contract.  It is responsiveness constrained but has zero total/prominent quality-floor violations and zero visual-importance debt.  Threshold-stamped CAD timing prevents the prior 32/64-pixel balancing loop and stops faceplate-only frames from inflating CAD capacity.  The post-terminal-authority endpoint at `/tmp/qged-150k-terminal-planner-clamp-20260828.json` settles terminal and box-free in 29.3 seconds with 686 meshes, 149,344 subpixel occurrences, no owner/obligation/violation, and matching requested, active, and certified budgets.  A later trace found an unchanged structural-distribution seed reopening the already constrained point cut after its finer preload was rejected.  The idempotent reducer fix and formal contract are qualified by `/tmp/qged-production-150k-idempotent-seed-20260829`: the final-binary warm System-GL wire interaction passes in 27 seconds, remains box-free, and each constrained finer candidate is consumed once with no identical terminal replay.  This is liveness/performance evidence rather than final realistic visual-significance clearance. | current shaded System-GL/OSMesa liveness green; true-cold/wire rows and realistic prominent-object quality remain |
| System GL smoke | The current shaded Generic Twin cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 12/11 seconds after camera identity became an optional snapshot.  Both exact-view returns recall history and every terminal checkpoint is ready with zero structural boxes and occurrence failures; the final certified errors are 0.236/0.241 pixels.  The 2026-08-27 post-ownership full System GL interaction replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-system-report.json` passed in 9.7 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The post-terminal-authority replay at `/tmp/qged-generic-system-terminal-planner-clamp-20260828.json` passes in 9.5 seconds and ends ready with 709 meshes, 135,073 faces, 57,102 lines, zero boxes, and no owner/obligation/violation.  The prior wire evidence is `/tmp/qged-generic-wire-system-debug-20260824`; both terminal wire frames contain 709 mesh payloads, zero boxes/failures/pending work, and matching cold/warm camera state.  The exact-current Lucy shaded cold/warm interaction replay at `/tmp/qged-lucy-system-20260824` also passed: both final frames are exact/ready, box-free, and quality-floor-clean with one resident source payload, 7.41M displayed PoP faces, a 230 MB resident mesh set, and one quality-history recall.  The final-binary Lucy true-cold/certified-warm wire pair at `/tmp/qged-production-lucy-wire-system-20260828` passes in 39/16 seconds.  Both runs finish with the same 3,183,110-face cut, 9,549,330 wire lines, no boxes, pending work, control violation, or quality-floor miss, and matching camera state.  Hubble shaded/wire cold/warm pass after the overview lifecycle repair. | Lucy shaded/wire green; complete the remaining System-GL large-model matrix |
| OSMesa Generic Twin | The current shaded cold/warm pair at `/tmp/qged-optional-policy-signature-generic-20260826` passes in 15/14 seconds.  Both exact-view returns recall history and terminate ready with zero boxes/failures; the final certified errors are 0.621/0.676 pixels.  The 2026-08-27 post-ownership full replay at `/tmp/qged-generic-growth-handoff-20260827/formal-contract-osmesa-report.json` passed in 14.3 s: all 12 waits were terminal/ready, every wait was box-free with at least 673 meshes, and the final frame had 709 meshes and 135k faces.  The final scene content and camera match the paired System-GL images.  The 2026-08-26 cross-renderer wire replay at `/tmp/qged-generic-wire-hud-final-20260826` passes cold and warm on System GL and OSMesa after the availability/HUD changes, with no invalid ready-label sample or terminal box.  The exact 2026-08-29 cold shaded cross-renderer replay at `/tmp/qged-production-generic-shaded-idempotent-20260829` passes in 11 seconds on System GL and 14 seconds on OSMesa; both end with all 709 meshes, approximately 135k faces, a one-pixel point threshold, and zero structural fallback.  The 2026-08-30 isolated-cold and same-cache hidden-line runs at `/tmp/qged-generic-osmesa-m4-{cold-r12,warm-r13}` pass in 11/13 seconds after cost-matched structural timing and exact-presentation recertification.  Every checkpoint contains all expected meshes with zero boxes/proxies/control violations; libicv SSIM is 0.998978--0.999833 with 0.0997%--0.2804% silhouette disagreement.  The differing cold/warm renderer costs are expected flat-batch versus retained-VBO submitted-work currencies and are consumed only with their paired durations. | compact direct path green across shaded, wire, and hidden-line cold/warm; continue larger real-model matrix |
| spatial Lucy, System GL | Exact-current 2026-08-25 warm shaded retained-page replay passes.  It refines during continuous zoom, reaches its quiet pixel target, compacts on zoom-out, restores the prior view through quality history, and has no boxes or page-level subpixel proxies. | focused retained-page regression green; cold and wire rows remain |
| direct-mesh Generic Twin | The 2026-08-29 cold replays at `/tmp/qged-generic-promotion-{cold,osmesa}` pass on System GL and OSMesa.  Late safe-scene certification now promotes an already visible PoP payload atomically by marginal replacement cost.  Each terminal frame has all 709 BoT occurrences in direct full detail, 185,388 faces, zero progressive/proxy/structural-box payloads, and no pending work.  The overall extent preview remains a discovery-only presentation and is not counted as a semantic leaf. | focused admission regression green; keep discovery-preview latency distinct from terminal-mesh admission |
| matched terminal image quality | The matched full-detail/LoD shaded matrix uses the same qged binary, per-orientation camera, physical canvas, lighting, and renderer, with libicv SSIM/PHASH plus exact and one-pixel-tolerant silhouette disagreement.  Generic Twin is effectively identical to full detail; Hubble, certified-warm Lucy, and NIST BREP satisfy their pixel-target or typed-constrained contracts without missing topology.  The current 4.96M-face heterogeneous 5k control compares at 0.978799--0.990360 SSIM on System GL while managed LoD presents 1.56M--2.26M faces; its tolerant silhouette disagreement is 0.004%--0.031%.  OSMesa presents 988k--1.66M faces at 0.943401--0.981091 SSIM under an explicit 169--333 ms responsiveness constraint.  The partial real-world Big Boy BoT conversion retains 14%--23% of 6.93M instance-expanded faces at 0.968001--0.982677 SSIM on System GL.  Its corrected OSMesa path reaches 0.946401--0.973746 SSIM without boxes, terminal proxies, or protected-floor debt.  A cross-renderer comparison found that non-exact PoP cuts disabled culling but fixed-function OSMesa still used source-oriented one-sided lighting; deriving lighting from the displayed cut's culling contract removed the black boiler/cab/roof patches and is covered by an Obol renderer regression.  Named wheel, running-gear, and boiler/cab crops show that the remaining high exact side-view disagreement is predominantly one-pixel boundary movement rather than a missing prominent component.  The 2026-08-31 wire and shaded-with-edges comparison adds Generic Twin and Hubble on both renderers.  Generic Twin reaches the complete mesh at 0.998978--0.999986 SSIM.  Hubble wireframe has below 0.008% one-pixel-tolerant silhouette disagreement; shaded-with-edges remains below 0.65%, including explicitly responsiveness-limited OSMesa cuts.  Exact metrics, content-digest controls, safe large-control rules, and provisional targets are in `obol_lod_visual_quality.md`. | simple/moderate shaded, wire, and shaded-with-edges; heterogeneous 5k; and partial real-world Big Boy cross-renderer baselines green; full BREP train, multi-large-mesh image comparison, and independent production-vehicle qualification remain |
| 50k/150k visibility mutation | The final warm 50k OSMesa workflow at `/tmp/qged-selection-capacity-50k-osmesa-final-20260831` selects in 0.93 s, erases 16 exact occurrences in 1.36 s, restores them in 1.40 s, and deselects in 0.90 s.  Visibility changes 50,000/49,984/50,000 while the terminal retained scene keeps 1,206,806 faces, zero boxes, and the same 5,130,680-unit certified budget; the complete matrix and camera contract pass.  A later timing-sensitive replay retired three mesh payloads while the subpath was hidden and exposed a missing prerequisite: source visibility became current before the successor framebuffer had classified restored occurrences.  The controller now enforces exact census, then exact presentation/classification, then reallocation.  `/tmp/qged-visibility-prerequisite-50k-osmesa-20260901` restores all 50,000 occurrences in 0.83 s and terminates with 1,343 meshes, 48,657 subpixel occurrences, zero boxes, and no owner or obligation.  The 150k OSMesa qualification at `/tmp/qged-selection-capacity-150k-osmesa-fixed-20260831` selects in 0.65 s, erases in 6.05 s, restores in 3.96 s, and deselects in 0.63 s.  It retains 835,318 faces, zero boxes, and the same 1,446,651-unit budget within one rounding unit across 150,000/149,984/150,000 visible occurrences.  All terminate ownerless with monotonic visibility revisions and pass the six-domain runtime trace.  Exact visibility deltas avoid a full inventory rescan and preserve renderer-capacity evidence. | shaded OSMesa hierarchy selection and exact mutation green at 50k/150k; wire and nested/deep hierarchy rows remain |
| ordinary-model terminal latency | The final Generic Twin cold/warm interaction repeats reach the first terminal mesh presentation in 0.71--0.93 seconds on System GL and 1.27--1.55 seconds on OSMesa.  Each endpoint has 673 active CAD payloads, zero structural boxes, no background/progressive obligation, and a ready convergence witness.  A later intermittent OSMesa run reached 99% after frame timing opened capacity/handoff work without arming another Qt timer.  Frame completion now re-evaluates the level-triggered host-work predicate; `ObolHostWork` checks that transition, and three consecutive untraced OSMesa quality runs terminate in 8--9 seconds with 134k--135k faces and no surviving qged process.  A separate Lucy race showed provider teardown clearing this shared pump after its final mutation left exact-frame debt; teardown now synchronizes every independent obligation, the executable regression passes, and three repeated System-GL Lucy lifecycles settle identically in 4.42--4.50 seconds.  The gate measures draw-to-idle rather than a fixed-delay screenshot.  Its nested-loop presenter forces the real QOpenGLWidget paint path, and a one-command test batch no longer disables its only canvas update.  Structural timing evidence is now an inseparable `(presented cost, duration)` pair, so a cheap box-only frame cannot authorize a richer current scene.  Compact direct scenes use a 400 ms prominent replacement deadline while large censuses retain their 200 ms finite deadline, and an exact presentation above a stale allocation certificate re-enters capacity reconciliation instead of becoming an ownerless 672/673 endpoint. | dual-backend cold/warm latency, exact presentation controls, cost-matched structural admission, and completed-frame wakeup green |
| concurrent cold cache | On 2026-08-28 two System-GL qged processes opened the same `Generic_Twin.g` and one initially empty `BU_DIR_CACHE` simultaneously.  Both exited successfully with terminal ready, box-free mesh presentations and no occurrence failure.  `mdb_stat` validated all three resulting LMDB environments (117 LoD, 83 draw, and 4 name-map entries), and an independent third process reopened the shared result, loaded from cache, and again converged terminal/ready without a cache write. | cross-process cache correctness green; large-asset duplicate-generation resource pressure remains to qualify |
| Hubble OSMesa | Exact-current 2026-08-25 shaded cold/warm lifecycles and camera contracts pass after the static-quality restoration change, with 1,804 warm shaded payloads and zero structural boxes.  The final-binary cold/warm wire pair at `/tmp/qged-production-hubble-wire-osmesa-20260828` passes in 16/16 seconds.  Both finish terminal, box-free, ownerless, and quality-floor-clean with 1,965 mesh payloads plus 553 subpixel occurrences.  The hierarchy replay retains one selected CAD instance across exact subpath erase/redraw, loads 2,530 tree items without a fallback scan, and applies its path update in at most 436 microseconds.  `/tmp/qged-resize-hubble-body-selection-final-20260829` adds four shaded managed/full-detail resize runs on System GL and OSMesa using the substantial `all.g/BODY` region rather than the old small probe.  Selection survives every window transition and the resize storm; the terminal selected frame is box-free and idle, clearing it changes 4,763--4,925 scene pixels, and all four final cameras match exactly.  The visually inspected frames show the same body panels receiving selected styling on both renderers.  The 2026-08-30 selected-subpath erase replay exposed allocator output feeding back through its external-cost input.  `/tmp/qged-hubble-selection-fixed-report.json` now settles after three bounded OSMesa allocation steps, and `/tmp/qged-hubble-selection-fixed-system-report.json` is visually stable from the first 250 ms checkpoint; both end terminal with no owner, obligation, or render request. | shaded selection/resize and shaded/wire cold/warm focused gates green; complete wire resize and broader physical-pointer modes |
| NIST BREP | The final-binary cold/warm wire lifecycle at `/tmp/qged-nist-resident-wire-final-20260828` passes on System GL and OSMesa.  All four runs select cut 10 with 19,457 lines, refine to cut 11 with 28,892 lines on zoom-in, and return to cut 10/19,457 lines on zoom-out.  Each endpoint is terminal with one direct progressive occurrence, zero structural boxes, and zero LoD-service tasks; cold HUD state remains `Discovering model` until that occurrence arrives.  The focused source-mutation test replaces the resident progressive part through the sparse delta journal and proves its cut and aggregate metrics retire.  The paired shaded OSMesa regression at `/tmp/qged-nist-shaded-regression-20260828` remains box-free and view-responsive.  The current resize/policy matrix at `/tmp/qged-resize-nist-contract-final-20260901` passes all eight System-GL/OSMesa, managed/off, shaded/edge workflows in 5--12 seconds; inspection confirms matching shaded silhouettes and the expected richer adaptive edge tessellation when LoD is enabled. | adaptive wire cold/warm and dual-backend shaded/edge resize gates green; shaded zoom-memory reclamation and larger real-BREP qualification remain |
| multi-Lucy/xpush capacity | Exact-current 2026-08-26 warm initial-view runs pass on both renderers after capacity samples were ordered behind ceiling-free occurrence-plan handoff.  The pre-fix timing replay split 2/4 between pixel demand and a premature ordinary-deadline endpoint; after the applied-but-hidden candidate fix, 6/6 early-checkpoint OSMesa runs reached the identical approximately 984.1k-face, 0.457-pixel pixel-demand endpoint in 4.59--7.04 s with no boxes or control work.  Post-recovery-ownership OSMesa replays remain green: multi-Lucy reaches 985,820 faces before a performance-limited zoom endpoint, while xpush reaches 985,808 faces initially and 1.53M after zoom, all box-free and terminal.  System GL reached 3.81M faces at 0.228 pixels in 2.45 s.  The 2026-08-27 full 75-event OSMesa turnover/zoom/rotation replay at `/tmp/qged-multi-lucy-growth-handoff-20260827/formal-contract-report.json` completed in 29.2 s after completed-pass ownership and annotation lifetime were made atomic.  All waits settled terminal/ready with zero failures or control violations; the final view had eight meshes, 319k faces, zero boxes, and no remaining owner or obligation. | warm full-interaction regression green; cold eight-distinct-asset preparation, compaction, and wire rows remain |

The 2026-09-01 post-visibility-contract image qualification supersedes the
older shaded point measurements in the table without weakening the remaining
release rows.  Fresh byte-matched full-detail/LoD comparisons cover Generic
Twin, Hubble, Lucy, heterogeneous 5k, partial Big Boy BoT, and NIST BREP on
System GL and OSMesa; all 48 checkpoints satisfy the exact presentation,
geometric error, topology, proxy, prominent-floor, and constraint-evidence
contracts.  Managed-only multi-Lucy, 50k, and 150k tiers also terminate
box-free and proxy-free with every classified prominent occurrence meshed.
The current 150k System-GL/OSMesa endpoints take 218/84 seconds, present
20.49M/1.12M faces, and peak at 9.27/7.52 GB process RSS while retaining
688/79 MB of controller-accounted mesh data.  Complete reports and realistic
expectations are in `obol_lod_visual_quality.md`; artifacts are under
`/tmp/qged-lod-quality-post-visibility-*-20260901`.

The subsequent terminal-handoff audit found two concrete refinement defects,
not failures of the abstract convergence rule.  A stronger selective pass
could retire the local cursor while leaving scene-wide quality debt, and a
point-threshold change could create a new allocation without advancing its
capacity revision.  The demand obligation is now level-triggered, all
classifier setters publish their semantic revision, and the runtime producer
validator recognizes complete inventory/submission coverage without accepting
an arbitrary selective submission.  The matched Hubble rerun at
`/tmp/qged-lod-quality-hubble-contract-final-20260901` passes modes 0 and 4 on
System GL in 8/9 seconds and OSMesa in 18/27 seconds.  Every endpoint is
terminal and ready with all 2,030 payloads satisfied, zero terminal proxies,
and a clean external control trace; the OSMesa hidden-line return is explicitly
responsiveness constrained rather than silently declared pixel-exact.
Focused point-selection replays at
`/tmp/qged-hubble-point-selection-contract-current-{system,osmesa}-20260901.json`
also pass the external trace.  Selection preserves the view and policy
revisions and the same 1,782/1,782 satisfied population; both renderers are
already idle and ready with zero boxes, proxies, obligations, service work, or
control violations at the first GUI readback and remain so after two seconds.

The 2026-09-01 capacity/HUD audit fixed two further refinement defects before
changing presentation text.  The protected visual floor is now the capacity
search's exact lower bound and capacity-owned allowances remain frozen for the
search key; this removes the out-of-domain candidate/floor alternation.  Exact
frame debt is distinct from exclusive interrupted replay, so a capacity
candidate allocates before its downstream frame.  Finally, a selective source
scope now retires with its completed cursor instead of blocking the broader
result-demand successor.  The final warm 5k OSMesa gate at
`/tmp/qged-lod-hud-unique5k-osmesa-final-20260901` passes in 32 seconds with
987,826 faces, zero boxes/proxies, a terminal ownerless endpoint, and a clean
325-transition external trace.  Its bounded capacity rank reaches 40/40.

The subsequent Hubble OSMesa mode-4 audit found two concrete-to-contract
mapping defects.  A complete pixel-demand endpoint was being assigned a new
plan serial when only its now-irrelevant protected-floor trial changed, and a
missing coverage pass could clear the demand-refresh producer of a still-live
importance census.  Allocation now canonicalizes that complete endpoint, and
only a completed dense ordinary pass may retire physical-demand refresh.  The
exact reproducer at
`/tmp/qged-contract-fixed2-hubble-osmesa-m4-20260901` reaches terminal/ready in
26.3 seconds with all 130 transitions clean.  The focused HUD replay records
capacity search as an exact local rank (for example, `search 13%`, `probe 2/8`)
alongside whole-view progress and names final allocation, publication, and
presentation owners instead of describing every last-percent handoff as
`Improving view 99%`.

The final wire/edge replay at
`/tmp/qged-lod-quality-wire-edge-contract-final-20260901` then passes Generic
Twin and Hubble modes 0 and 4 on System GL and OSMesa: all eight rows terminate
ready in 4--22 seconds with zero structural boxes, terminal proxies, prominent
floor violations, or runtime-contract violations.  All sixteen Generic Twin
views meet their strict safe-direct image targets.  Hubble's software
hidden-line endpoint is explicitly responsiveness constrained; paired-image
inspection retains every major component and the long cable rather than
hiding a topology loss behind aggregate image metrics.

The complete post-contract visual ladder is retained under
`/tmp/qged-lod-quality-post-contract-*-20260901`.  Fresh matched shaded runs
cover Generic Twin, Hubble, NIST BREP, Lucy, 5k distinct meshes, and the
partial Big Boy BoT on both renderers; all 48 views pass their image,
presentation, and external control-trace contracts.  Managed two-L Lucy,
50k, and 150k tiers terminate box/proxy-free with no prominent-floor miss.
The current 150k System-GL/OSMesa initial views finish in 76/67 seconds under a
12 GiB process limit, reach their exact 40/40 capacity rank, and retain all
current prominence candidates.  System GL reduces the preceding terminal
population from 20.49M to 4.51M faces while changing only 0.037 percent of the
one-pixel-tolerant silhouette; the OSMesa change is 0.574 percent.  This is
evidence that the screen-space allocator can remove substantial invisible
work without losing the whole-model envelope.  It does not replace the open
independent production-vehicle significance gate.

The quality harness now emits executable matched and managed assessments.
Safe-direct SSIM is a hard contract.  For larger shaded scenes, SSIM remains
an advisory photometric target while pixel error and one-pixel-tolerant
silhouette are independent geometry gates; coarse normals may change lighting
without moving a silhouette.  Managed-only telemetry can establish a valid
bounded endpoint but is always marked for visual review because it has no safe
whole-scene oracle.

The 2026-08-31 50k selection-stabilization replay exposed a missing render
classification.  Capacity relevance alone could not distinguish an LoD-owned
structural classifier from a selection/style presentation; both are
non-capacity frames when no measurable mesh population is present.  The frame
transaction now carries independent capacity and LoD-planning relevance.
Selection is excluded from both, while the classifier remains planning
relevant.  `ObolSemanticPresentation` proves that its isolated exact successor
can manufacture neither kind of evidence.  The strengthened graphical matrix
waits for selection readiness before erase.  Warm shaded 50k replays pass on
OSMesa and System GL under
`/tmp/qged-50k-selection-planning-{final,system}-20260831`.  The final OSMesa
binary retires the exact selection frame in 0.69 seconds to owner/obligation
zero without changing its 1,195,385-face population.  System GL returns in
0.19 seconds to the same
pre-selection background-compaction owner and unchanged 14,493,432-face
population.  Neither opens a capacity sample or headroom-refinement frame.

The 2026-08-31 selection follow-up also covers a cold shared-mesh publication
edge which the earlier single-object UI probe did not exercise.  A provisional
prepared preview was incorrectly treated as a resident progressive binding,
so one occurrence could keep projection refresh alive after the shared asset
was ready.  The capability-based classifier and direct shared-asset
replacement now pass the focused update-action test, a true-cold eight-instance
Lucy lifecycle, the dual-renderer selection UI matrix, and cold/warm Hubble
hierarchy selection/erase/redraw.  Every terminal observation has zero owner,
obligation, queued task, and in-flight task; all eight Lucy occurrences own
progressive payloads.  Reports are retained under
`/tmp/qged-selection-stabilize-*`; regenerated cache payloads were removed.

The subsequent full Hubble OSMesa hierarchy replay showed that the selection
frame was healthy but selected subpath erase/redraw could still expose a
capacity/triangle-recovery ownership violation.  A stronger `ALLOCATING`
capacity candidate had been represented by the generic preserve-budget request
and was clamped to the recovery ceiling; every exact frame then rejected that
different population as stale and restarted the search.  Capacity candidates
now have a distinct retained-allocation request kind which weaker recovery,
deadline, and reconciliation requests cannot replace.  The focused policy
test exercises the simultaneous recovery ceiling, and
`ObolPointQualityOwnership` passes 66 generated / 36 distinct states to depth
8.  Fresh exact-current warm Hubble workflows pass on System GL and OSMesa at
`/tmp/qged-hubble-post-selection-capacity-owner-final-{system,osmesa}-20260831`.
Selection returns to idle without changing its 620,449-face population, and
erase/redraw/clear-selection all terminate without an owner, obligation, or
capacity-sample loop.

The latest high-cardinality System-GL evidence supersedes the older terminal
populations in the cumulative scale rows above.  The bounded-batch 150k shaded
replay at `/tmp/qged-production-current-150k-batched-20260830` reaches a
whole-model preview during discovery, then terminates pixel-target and box-free
with 91,710 mesh payloads, 61,771 subpixel occurrences, 9.54M triangles, about
642 MB resident mesh data, and no prominent-floor or visual-importance debt.
Its first pose takes 179 seconds and the complete interaction workflow 229
seconds, so tens-of-thousands-of-distinct-mesh throughput remains an
optimization target even though the view is useful much earlier.  The focused
structural-frontier model now covers the bounded batches and their strictly
decreasing exact remainder over 60,535 generated / 29,266 distinct states to
depth 386.

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
OBB path does not defeat direct-to-mesh admission for a modest scene.

The heterogeneous 10k fixture at
`/tmp/obol-current-cache-matrix/unique-obb-10k` mixes mesh size, topology,
color, hierarchy, and baked non-axis-aligned local geometry.  A fresh two-pass
prewarm generated all 10,000 PoP hierarchies (8,349 then 1,651) and retained
1,960 materially tighter PCA boxes in both the asset and draw caches.  Under a
deliberately constrained OSMesa policy, the cold-manifest and warm-manifest
runs select the same 693 OBBs, 2,856 AABBs, 6,442 aggregate points, and nine
mesh payloads, with 29,571 presented faces, no structural fallback, and a
terminal ready state in 1.78/1.83 seconds.  The paired warm System-GL run
instead admits 6,962 mesh payloads and represents only 3,587 subpixel
occurrences as points; it terminates box-free with 920,108 presented faces.
The GUI matrix now derives its smooth-zoom target from database bounds, so
these images cannot pass by examining an empty fixed coordinate.  Repeated
explicit prewarm also skips valid current cache entries instead of resubmitting
the same prefix indefinitely.

The independently versioned draw metadata was then migrated over the preserved
current-format PoP caches instead of erasing 8.8 GiB of valid hierarchy data.
The final System-GL warm interaction replays pass at both pressure scales:
`/tmp/obol-current-cache-matrix/unique_mesh-50k-coverage-contract-warm`
finishes ownerless and box-free in 29 seconds with 34,746 active mesh payloads,
17,906 subpixel occurrences, and 4.59 million presented faces;
`/tmp/obol-current-cache-matrix/unique_mesh-150k-coverage-contract-warm-v2`
passes in 62 seconds with 13,937 resident mesh payloads, 97,087 aggregate
points, 8.86 million presented faces, and zero structural boxes or control
violations.  The 150k final framebuffer is already terminal and usable while
the HUD correctly reports its remaining bounded resident-prefix compaction as
background memory/cache work.  A policy-only revision no longer manufactures
a redundant coverage census; an active current-view census is now explicit
foreground inventory work and is reattached to a bounded cursor if a stronger
owner had paused it.

This clears persistence, cross-renderer semantics, and heterogeneous synthetic
pressure as implementation risks.  Real wheels/blades/booms/hulls still need
perceptual comparison against no-LoD ground truth.  Multi-page terminal proxy
ownership and optional enrichment of a manifest sealed before first-cold PoP
characterization also remain open; the active cold payload already carries
the OBB, but the initial structural cue deliberately remains an AABB.
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

The 2026-09-01 BREP fallback and progress pass adds a current-tree 64/64
`drawing_baseline` run, a 28/28 `bobol_headless` run, and focused NIST PMI7-10
shaded and wire real-model passes.  The shaded real-model test finishes from
an isolated writable cache in about twelve seconds and the wire case in about
six seconds.  A true-cold OSMesa GUI replay exposes the exact four-item source
preparation rank, replaces each temporary box with its completed mesh, and
terminates with no boxes.  The matched managed/control image reaches SSIM
0.998831 and zero one-pixel-tolerant silhouette disagreement while using fewer
faces than the control.  Cache-context failure is now separately diagnosed
and preserves valid tessellation as terminal non-LoD geometry.  This closes
the misleading PoP-rejection and indefinite-box failure modes.

The subsequent constrained-residency lifecycle at
`/tmp/qged-lod-brep-pmi7-10-reclaim-trace15-20260901` (System GL) and
`/tmp/qged-lod-brep-pmi7-10-reclaim-osmesa-20260901` (OSMesa) qualifies the
NIST shaded BREP zoom/reclaim/return/restore path against matched LoD-off
controls.  Both backends present all four targets with zero structural boxes,
terminal proxies, or prominent-floor violations.  The close view and its
cache-backed restore both present 220,637 faces.  System GL records SSIM
0.993568 and OSMesa 0.994099--0.994100; every checkpoint has zero silhouette
disagreement after a one-pixel tolerance.  The deliberately reduced resident
limit produces an honest memory- or responsiveness-constrained terminal
status without losing the pixel-target presentation.  A changed resident
capacity epoch now permits one byte-governed physical-demand retry past its
stale allocation; previously that retry scheduled zero work and could publish
a false coarse terminal state.  NIST adaptive shaded BREP residency is green;
the larger original Big Boy BREP hierarchy remains an open release row.

The HUD now reports finite source preparation and capacity-search ranks rather
than parking at a generic `Improving view 99%`.  A warm OSMesa NIST probe
(`/tmp/qged-lod-progress-probe-report2.json`) showed the source rank advancing
0/4, 1/4, 2/4, and 3/4 before terminal readiness.  The fourth BREP dominates
the remaining wall time, so this is deliberately a finite work rank rather
than a fabricated ETA; the label continues to identify the active producer.

The final current-tree visual spot audit after the background-result handoff
repair is retained under `/tmp/qged-lod-quality-*-20260901` and summarized in
`obol_lod_visual_quality.md`.  Generic Twin and Hubble pass matched control/LoD
wire and shaded comparisons on both backends; NIST BREP and the partial Big Boy
BoT pass their matched System-GL rows; and single/multi-Lucy plus shaded
50k/150k managed endpoints are coherent, terminal, box-free, proxy-free, and
prominent-floor-clean.  A queued cache result can no longer demote an already
resolved semantic revision from terminal to foreground refinement.  The 50k
wire endpoints are correct but explicitly constrained at 105--145 ms per
frame and approximately 2--3 pixels maximum projected error.  Improving
high-cardinality wire submission/raster throughput without changing visible
edge semantics is therefore an optimization gate, not a prerequisite for the
shaded quality result.  The HUD reports `render primitives` from the latest
exact completed-frame record, so direct wire and full-detail channels are no
longer mislabeled or reported as zero.  It falls back to the active PoP face
estimate only while no exact frame record is available.

The next BREP production-quality gate is the larger original Big Boy hierarchy:
complete its shaded source conversion, repeat the now-passing NIST constrained
close/zoom-out/cache-restore lifecycle, and use named train features for an
independent real-model significance review.  Cache-unavailable, invalid-input,
and PoP-construction failures still share a low-level integer return API;
replacing that ambiguity with typed diagnostics is cleanup debt, although valid
BREP tessellation is now preserved correctly in every case.

The September 1 single-wheel probe makes the source-conversion prerequisite
concrete.  The legacy CDT silently aggregated a partial result after five
constrained-edge face failures: only 32 of 368 allocated points were
referenced, and the resulting 128 indices covered one upper plane rather than
the tire's full thickness.  PoP classified that defective source faithfully;
it did not reject or repair it.  The indexed-face provider now fails closed
when referenced output omits the exact BREP boundary extent or escapes the
surface envelope, and shaded BREP cache version 3 invalidates the old partial
payload.  The enhanced `brepdraw` audit shades the same primitive with all
eight faces represented by 1,336 triangles in 0.009 seconds.  A full-train
comparison remains blocked until both display policies share that bounded,
validated source contract, including deadline, memory/result limits,
cancellation, per-face completion, and typed failure.

The guard is selective rather than a global BREP fallback.  All four NIST
PMI7-10 solids pass the provider diagnostic, and the fresh-cache OSMesa
comparison at `/tmp/qged-lod-nist-brep-v3-20260901` completes in 18/24 seconds
for control/LoD.  It presents 178,561 LoD faces with SSIM 0.980645, zero
one-pixel-tolerant silhouette disagreement, and no boxes, proxies,
quality-floor violations, or terminal error.

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
