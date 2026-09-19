# libBObol active debt

Reviewed 2026-09-19. This is the sole current backlog. The
[simplification guide](obol_simplification_guide.md) owns the roadmap and
[production readiness](obol_production_readiness.md) owns release acceptance.
Previous closures and reproductions are preserved in the
[September 19 snapshot](obol_20260919_debt_history.md) and
[conformance audit](libbobol_tla_conformance.md). Those historical records do
not add work beyond the obligations below or qualify changed binaries.

## Ordered work and completion evidence

| ID | Gate | Remaining work | Closure evidence |
|---|---|---|---|
| BUILD-01 | S0 | Make a fresh build and repeated configuration select one coherent BRL-CAD/Obol/OSMesa stack; remove dependence on manual restaging | Exact source/options/dirty-patch manifest, loaded-library identity and passing smoke after both build and reconfiguration |
| API-01 | S1 | Finish the live mutation, callback, ownership and lifetime inventory across host, source/service, scene, GED and presentation; audit real callers | Each entry has one owner and contract; deprecated/unsupported use has a migration or removal decision; inventory is finite |
| FLOW-01 | S2 | Qualify cold draw, camera change during loading, cancel/close, reopen/redraw and terminal presentation through actual owners | Deterministic delayed/stale/denied-result tests plus graphical workflow, independent output checks and measured resource release |
| OWN-01 | S3 | Resolve the ownership violations found by API-01, applying FLOW-01's method | Removed duplicate writers/state, independent compilation, production transition tests and migrated callers |
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

Complete BUILD-01 and API-01 sufficiently to make FLOW-01 reproducible and
reviewable. The viewport demonstration removes an unsealed direct-field
observer layer and uses existing host/stream mutation methods. It does not
close API-01 for the rest of libBObol or FLOW-01 for asynchronous CAD drawing.
Do not resume an exhaustive audit of arbitrary raw viewport writes.

Ordinary production defects discovered by independent behavioral checks should
be fixed at their existing owner while this inventory is completed. Preserve
all earlier qualified publication and source-lifetime regressions except tests
for explicitly withdrawn, unsealed behavior. New nested-callback combinations
need a supported caller and an acceptance requirement before becoming blockers.

## Boundaries and retained requirements

API-01 must distinguish detached construction, live owning mutations, derived
outputs, borrowed inspection, worker result delivery and application observers.
The migrated viewport family is listed in the
[API contract](libbobol_api_contract.md#supported-mutation-inventory).
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
The rebased source-publication sweep also failed a scene-light input assertion
while other qualification ran concurrently. The isolated selector and a later
full source sweep pass. Preserve both logs and investigate the scheduling or
shared-test-state dependency before declaring that witness closed.
The desktop frame-delivery reproduction stopped at event 28. Quiet planning
cycles must be reproduced with complete traces rather than inferred from HUD
labels. Historical raw `/tmp` evidence is unavailable and must be regenerated.

QUALITY-01 includes the recorded Lucy close floor (cut 28, 2,619,533 faces,
error 3.5813), System-GL close floor and software cut-25 zoom starting conditions.
Do not retire one using a pass from another starting cut. Numeric, renderer
and resource failures belong to their owners; proven-constrained visual debt
must not reopen an exhausted control search. Named wheels, blades, booms,
tails and hulls need full-detail comparisons and inspected images, alongside
SSIM/PHASH and silhouette metrics from the visual-quality contract.

QUALITY-01 and EDIT-01 also retain the seven upstream annotation image
comparisons. Command/update/color assertions pass and duplicate model-space
anchor translation is repaired, but framing, screen-space behavior and full
text/fill/style parity still require visual qualification. The original image
controls and thresholds are unchanged; see the [integration record](obol_main_integration.md).

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

S1--S3 close structural simplification only. Production readiness requires all
required S4--S6 rows on the exact final candidate. A new counterexample maps to
an existing row; add a row only for a distinct obligation. Report concrete
passed/open criteria, never subjective percentages or test counts as maturity.
