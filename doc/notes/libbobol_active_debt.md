# libBObol active debt

Reviewed 2026-09-20. This is the sole current backlog. The
[simplification guide](obol_simplification_guide.md) owns the roadmap and
[production readiness](obol_production_readiness.md) owns release acceptance.
Previous closures and reproductions are preserved in the
[September 19 snapshot](obol_20260919_debt_history.md) and
[conformance audit](libbobol_tla_conformance.md). Those historical records do
not add work beyond the obligations below or qualify changed binaries.

## Ordered work and completion evidence

BUILD-01/S0 is demonstrated for the September 19 source-built Linux candidate.
The fresh native install, repeated native/BRL-CAD configuration, default loader
identities, nine native suites and 39 affected BRL-CAD tests pass; the drawing
lifecycle regression also passes after the annotation cancellation repair.
Exact source patches, options and binary identities are retained under
`.build-main/obol-qualification/20260919-native-baseline`. Recheck the affected
evidence after dependency changes; this does not qualify other platforms or
the final release. Use the [integration recipe](obol_main_integration.md#reproducing-the-native-dependency-baseline).

API-01/S1 is demonstrated for the current checkout. The finite
[supported mutation inventory](libbobol_api_contract.md#supported-mutation-inventory)
distinguishes detached construction, live owner publication, worker delivery,
derived output, borrowed inspection and callbacks. The production caller audit
removed the remaining attached-source field, compact-presentation,
wireframe-realization and dead child-attachment bypasses. Public header,
installed-consumer, source publication, selection and GED synchronization
checks are retained under
`.build-main/obol-qualification/20260920-api-inventory`. A new live mutation
family must update that inventory and its acceptance check; final-candidate
requalification remains RELEASE-01.

FLOW-01/S2 is demonstrated for the same checkout. Named headless and graphical
workflows cover cold draw, camera input while work is active, endpoint/widget
close, reopen/redraw, changed terminal output and transient resource release.
The full GED row retains delayed, stale and denied-result coverage. Both focused
lifecycle rows passed ten consecutive runs; the combined focused/broad rows,
full build, package consumer, annotation, repository and license checks pass.
Evidence and limitations are recorded under
`.build-main/obol-qualification/20260920-flow-01`.

| ID | Gate | Remaining work | Closure evidence |
|---|---|---|---|
| LIFE-01 | S4 | Close source/service and borrowed host/widget lifetime gaps; repeat worker-active close, endpoint replacement, plugin reload and cancellation | Shared dynamic-stack ASan/UBSan and supported native TSan/LSan; no late publication, retained leases or lost wakeups |
| CONTROL-01 | S4 | Reproduce and close quiet planning churn and desktop requested-frame stall; investigate the load-sensitive exact-frame witness report | Stable-input termination with finite witnesses, complete transition journal and repeatable adverse scheduling tests |
| BREP-01 | S4/S5 | Replace aggregate-success CDT use with one bounded, typed BREP provider contract used by LoD on/off | Deadline, memory/result limits, cancellation and per-face completion; full original large BREP hierarchy and NIST cases pass |
| CACHE-01 | S4 | Measure simultaneous true-cold large-asset processes sharing a cache, cancellation during persistence, corruption/reload and compaction | Bounded memory/CPU/I/O and complete publication; add a crash-recoverable generation lease only if measured stampede violates bounds |
| QUALITY-01 | S5 | Close Lucy quality floors, close zoom, software cut/lighting issues and named-feature geometry fidelity | Both renderers pass declared image/coverage/normal/importance criteria over the retained failing starting conditions |
| SCALE-01 | S5 | Complete repeated/distinct multi-Lucy/xpush, 50k/150k and independent multi-gigabyte qualification | True-cold/warm discovery and useful-preview latency, per-phase costs, input/frame bounds, peak/resident/cache bytes and turnover recovery |
| EDIT-01 | S5 | Complete actual editing/selection/polygon/faceplate interaction and shared-client parity | Command/widget/manipulator/database agreement, invalid-input preservation, physical mouse paths and lifecycle tests |
| PLATFORM-01 | S5 | Complete required native hardware, DPR and client rows under declared support levels | Current System GL/OSMesa, single/quad, resize/fractional DPR, Windows and other required host results; deferred rows remain explicit |
| RELEASE-01 | S6 | Freeze and qualify the final candidate against every required matrix row | Exact source/dependency/binary manifest, retained reports/images, no open release blocker or unmeasured acceptance threshold |

## Immediate sequence

S0--S3 are demonstrated for this checkout. Begin S4 with the retained LIFE-01
and CONTROL-01 counterexamples. Keep FLOW-01 as the integration guard while
closing source/service and host/widget lifetime under sanitizers, then reproduce
the quiet-planning, requested-frame and load-sensitive exact-frame witnesses
with complete transition journals. A newly demonstrated duplicate writer or
split policy reopens OWN-01 as an S3 regression; cooperation across the declared
boundaries does not.

The first OWN-01 seam is complete. Service quiescence was independently rebuilt
from separately locked counters in GED, qged and test waits, and the headless
lifecycle version omitted active requests. Per-generation controller and
renderer decisions had the same coherence problem and could additionally call
a generation quiet while it still held a lease on another generation's shared
producer. `BObolLodService::workStatus()` and `generationWorkStatus()` now
return lock-consistent snapshots. Their status types own complete-idle and
generation-result-work predicates, including shared-producer leases. Production
callers use those predicates; phase-specific tests keep their explicit narrower
meaning. Generation zero now consistently reports an empty scope. The
controller's duplicate atomic result-pending flag is removed: the service queue
owns result availability, while the controller retains only the independent
first-ready timestamp used to bound publication latency. The control reducer
receives queue availability as an immutable input. The affected build, thirteen
focused rows, six controller/model rows, compact-publication stress, installed
consumer and both production lifecycle rows repeated ten times pass under
`.build-main/obol-qualification/20260920-own-01-service-work`.

The second OWN-01 seam is complete. Stable resident capacity is now published
as one atomic scalar rather than derived from separately sampled total and
reloadable-backing counters. Every realization, compaction, eviction and stop
writer updates it. `BObolLodService::residentCapacityStatus()` samples that
scalar with the mutex-owned growth reservation and limit. A completed growth
publishes exact stable bytes before releasing its reservation, so a concurrent
sample can conservatively count both but cannot expose the same capacity to two
producers. Renderer headroom and pressure policy, convergence status and qged
diagnostics use this contract. Saturating arithmetic and actual retained
publication, compaction, eviction/reload and stop transitions are covered; the
publication test changes the limit while observing a live growth reservation
and rejects any reservation-to-stable zero gap. The affected focused and model
rows, compact-publication stress, installed consumer and both production
lifecycle rows pass under
`.build-main/obol-qualification/20260920-own-01-resident-capacity`.

The third OWN-01 seam is complete. Qt previously reconstructed libged's LoD
progress visibility, terminal-readiness and phase-coalescing policy, then kept
the three results as separate cached fields. `BObolLodConvergenceStatus` now
derives one immutable `BObolLodProgressDisplayStatus` from its complete
snapshot. Libged consumes its visibility and terminal-ready results while
retaining faceplate text, color and geometry; Qt compares the whole value when
deciding whether a state transition needs publication and retains ownership of
cadence and event-loop scheduling. The local Qt enum, two mapping functions,
three split cache fields and duplicate libged readiness/visibility predicates
are gone. Direct classification, faceplate, production workflow, progressive,
controller/model, public/installed API and repeated lifecycle checks pass under
`.build-main/obol-qualification/20260920-own-01-progress-display`. The rendered
annotation contract did not change, so no control image was replaced.

The final S3 audit migrated `view lod service wait` and the GUI qualification
wait from separately sampled controller flags to `BObolHostWorkSnapshot`, while
retaining service quiescence as its separate worker/cache boundary. The
remaining apparent mirrors are deliberate single-writer transfers: the GED/Qt
host synchronizes renderer-neutral view input into an endpoint controller,
controller publication commits renderer limits in one operation, and source
admission receives the owner's point threshold as an input. Diagnostics read
individual fields but do not schedule or mutate from them. The affected qged
build, ten focused host/API/GED/Qt rows, and ten repetitions of GED cross-run
plus both production lifecycle paths pass under
`.build-main/obol-qualification/20260920-own-01-s3-close`. This closes OWN-01
and S3 for the current inventory.

Ordinary production defects discovered by independent behavioral checks should
be fixed at their existing owner while S4 work proceeds. Preserve
all earlier qualified publication and source-lifetime regressions except tests
for explicitly withdrawn, unsealed behavior. New nested-callback combinations
need a supported caller and an acceptance requirement before becoming blockers.

## Boundaries and retained requirements

API-01 distinguishes detached construction, live owning mutations, derived
outputs, borrowed inspection, worker result delivery and application observers
in the [API contract](libbobol_api_contract.md#supported-mutation-inventory).
No universal reentry guard or callback queue has been installed. Audit legacy
callbacks before changing their behavior; preserve committed state on observer
failure and contain exceptions at toolkit/C boundaries.

OWN-01 concerns view policy versus frame effects; source discovery/traversal
versus realization/publication; service task/lease, residency/reservation and
persistence; GED semantic reduction versus backend effects; and retained
presentation versus renderer execution. Some seams already meet the target.
Do not extract solely to reach a file-size goal, add forwarding wrappers or
move hot occurrence data into per-occurrence objects. Keep sparse deltas,
bounded publication, immutable geometry and the no-second-scene-sized-vector
constraint. Remove superseded branch APIs in the change that migrates callers.

CONTROL-01 retains the exact-frame ordering, six-domain identity, finite
obligation and completed-frame journal gates. The September 10
`libBObol_lod_update_action` witness failure appeared under concurrent test load;
predecessor and isolated current runs passed. Reproduce under controlled load
before deciding whether timing assumptions or production ordering are wrong.
The earlier scene-light input failure is closed. Point lights carried an unused,
indeterminate direction payload, and equality compared it as if it were
semantic; a NaN could therefore make an unchanged light unequal to itself and
replace its child during an enablement-only update. Equality now compares only
the fields used by each light kind and preserves an aliased input. A
deterministic unused-NaN regression reaches the original reentrant update and
passes in the complete ASan/UBSan source sequence.
The desktop frame-delivery reproduction stopped at event 28. Quiet planning
cycles must be reproduced with complete traces rather than inferred from HUD
labels. Historical raw `/tmp` evidence is unavailable and must be regenerated.

The expanded direct-primitive color test exposed a separate FLOW-01 failure:
a late whole-target overview replaced the bare root leaf's equal-tier box,
erasing its source request while convergence reported ready. The source
registry now rejects that replacement. The controlled delivery-order test
fails before the repair and passes afterward; ordinary direct draws also pass
20 consecutive repetitions. That counterexample is now also retained as part
of the demonstrated FLOW-01 boundary.

Annotation projection work exposed an OWN-01 defect in the existing per-view
assembly: it wrote private LoD part/placement/cut records into the shared source
array, corrupting later callback compilation. The repair derives view records
from staged presentation state and feeds complete resets through Obol's indexed
replacement reader; sparse updates retain bounded batches. The callback
regression fails before the repair and passes afterward. This removes a duplicate
writer without adding a second full occurrence array. The subsequent native
display-plane representation builds on that boundary; stroke and fill fidelity
still remain open.

QUALITY-01 includes the recorded Lucy close floor (cut 28, 2,619,533 faces,
error 3.5813), System-GL close floor and software cut-25 zoom starting conditions.
Do not retire one using a pass from another starting cut. Numeric, renderer
and resource failures belong to their owners; proven-constrained visual debt
must not reopen an exhausted control search. Named wheels, blades, booms,
tails and hulls need full-detail comparisons and inspected images, alongside
SSIM/PHASH and silhouette metrics from the visual-quality contract.

QUALITY-01 and EDIT-01 also retain the seven upstream annotation image
comparisons. Command/update/color assertions pass and duplicate model-space
anchor translation is repaired. Hide/show exposed an unrelated-source
cancellation defect: an erase retired pending jobs for labels still intended
to be visible. Existing per-source stream admission now owns that cancellation;
the added delivered-stroke regression fails before the repair and passes after.
Direct primitives now use the same material resolution as combination leaves;
legacy `rgb` and canonical `color` edits survive rendering and annotation
regeneration. The delivered-color and round-trip regressions pass and the
cyan leader is restored in the inspected image. Native display-plane projection
now keeps offsets fixed in pixels through independent cameras, zoom, rotation
and resize; ray/bounds tests and GED rectangle/export assertions pass. An
overview refresh now publishes complete leaf identity and selection state
through the existing source owner, and provisional overviews cannot seed
primitive geometry. The green screen-facing leader is visible in the new frame
four. The subsequent style repair retains authored width, color and pattern
runs without duplicating occurrences. Native pixels cover camera changes, live
base width, occurrence color replacement and same-ID part replacement; source
export preserves fractional widths, exact masks, colors and default resets.
The same style test now forces immediate GL on compatibility-profile OSMesa and
System-GL contexts and verifies tier-zero execution separately from direct
software wire. PostScript and PNG now
consume effective widths, masks, colors, alpha and annotation fills; plot uses
its documented named-pattern degradation and cannot represent widths or filled
areas; Qt object queries retain exact parallel style vectors. An isolated
command regression and the Qt export regression cover those mappings. Filled
areas retain librt triangulation in the existing native
triangle channel, draw unlit in every mode and participate in picking and
export. Hole/area, authored fill color/replacement, mixed stroke/fill and
same-ID role replacement checks pass. Semantic fill masks reproduce the active
solid or gradient background in fixed and GLSL rendering and retain an explicit
flag through library/Qt export. A
continued image run leaves frames one through five
and seven byte-identical; frame six restores formerly black-on-black authored
colors and improves from 0.886535 to 0.887488 SSIM. Its control retains a
separate framing difference at that checkpoint.
The [image attribution](obol_main_integration.md#annotation-image-attribution)
isolates frame seven's intentional camera change and the truck's explicit
material colors. A fitted-camera follow-up aligns the model-space annotation
pixels in frames one through six within one pixel of their old controls; the
remaining differences are the independently tested material, display-plane and
authored-style contracts. All seven attributed controls now pass at the
unchanged 0.99 threshold. Preserve those contracts rather than restoring
per-object bounds padding or region-table overrides. Viewport-scoped snap and
measurement now use the same display-plane transform as rendering and export;
ordinary traversal remains viewless, path-local measurement retains stored
offsets, and the legacy vlist publisher explicitly rejects screen annotations.
Cross-client visual parity and fill-boundary interaction semantics remain open.
The [annotation inventory](libbobol_api_contract.md#annotation-path-owners-and-representation-boundaries)
maps display coordinates, stroke style, fills and bounds/picking/export to
their owners; these remain representation work within the existing gates.

BREP-01's indexed-face check rejects the partial Big Boy tire and prevents
reusing its incomplete cache. This is containment, not correct tessellation.
A partial BoT cannot qualify the original BREP hierarchy. Preserve adaptive
wire and shaded growth/zoom-out/reclamation/cache-restore NIST regressions.

SCALE-01 includes cold local storage, warm page cache and slower storage;
globally representative spatial previews; shared-instance reuse versus truly
distinct assets; visibility turnover; and distributed realization/publication
throughput. Extend the heterogeneous OBB comparisons from 10k to realistic
features and 50k/150k pressure. Initial discovery-manifest OBB enrichment is
conditional on measured benefit and must remain bounded without nested cache
access. Occurrence-level aggregation for terminal multi-page proxying is
conditional on measured need; never attach one whole-source proxy per page.

EDIT-01 retains runtime legality/rejection/readback for every advertised edit,
ARB vertex/edge/face and primitive manipulators, specialized sketch interaction,
MGED `sed`/`oed`, polygon create/resize/move/boolean mouse paths, selection
modifiers and hierarchy-scale tree styling. Cover selection/highlight across
draw/erase/edit, actual cutting-plane affordance, plugin lifecycle and
fractional DPR. Existing deterministic measurement and settings round trips
are regression assets; they do not close native pointer/hardware qualification.

PLATFORM-01 follows the [platform matrix](libbobol_platform_threading_matrix.md).
Do not silently waive Windows corruption/zoom/lifecycle debt or expand deferred
macOS hosts. qged, MGED, gsh, Archer and rtwizard share the release behavior.
Programmable rendering remains opt-in until its own acceptance rows pass.

## Exit condition

S1--S3 close structural simplification and are demonstrated for this checkout.
Production readiness requires all
required S4--S6 rows on the exact final candidate. A new counterexample maps to
an existing row; add a row only for a distinct obligation. Report concrete
passed/open criteria, never subjective percentages or test counts as maturity.
