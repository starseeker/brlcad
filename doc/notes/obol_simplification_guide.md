# Obol drawing simplification guide

Reviewed 2026-09-20. This is the roadmap and finish line for the current
BRL-CAD/Obol drawing effort. The system is viable, but is not release-qualified.
Keep the working geometry, storage, rendering and policy foundations; reduce
how applications may mutate live state before extending publication machinery.

## Decision and working method

Use existing owning controllers, stores and hosts for application mutations.
Prepare a complete candidate, validate its identity, commit without callbacks,
then notify observers of committed state. Inventor fields support detached
construction and internal implementation; their existence does not promise
transactional live updates for every assignment, connection or callback order.
The [API contract](libbobol_api_contract.md#supported-mutation-inventory)
records the supported entry points and migration limits.

Application observers inspect committed state. Follow-up mutations belong in
a subsequent owner-thread turn; callbacks must not destroy an active owner or
let exceptions escape into a toolkit/C boundary. Existing defensive handling
remains until its actual callers and tests have been migrated. This is a
contract simplification, not a claim that every legacy entry already rejects
reentry or contains exceptions uniformly.

Start with one complete user workflow through the existing implementation.
Use its failures to select a responsibility boundary. Replace that boundary
only if it cannot be made understandable with fewer writers and less state.
Do not freeze all C++, build a universal controller, generate an implementation
literally from TLA+, or maintain a permanent parallel prototype.

Keep immutable geometry, compact occurrences, sparse semantic deltas, typed
revision identities, bounded publication, reservations and existing reducers.
Workers produce detached results; the owner thread authenticates and publishes
them. Keep requested, resident, prepared and actually presented detail distinct.
An allocation denial must name the prohibited operation and its evidence; it
does not invalidate already affordable resident detail.

## One role per document

| Question | Authority |
|---|---|
| What are we building and when are we done? | This guide |
| Who owns behavior and what is supported? | [Architecture](libbobol_architecture.md), [API](libbobol_api_contract.md), [pipeline](libbobol_progressive_pipeline_contract.md) and their boundary contracts |
| What work remains? | [Active debt](libbobol_active_debt.md), the sole backlog |
| What must pass and what has been measured? | [Production readiness](obol_production_readiness.md) and its visual/editing/platform contracts |
| Where do I resume this checkout? | [Session handoff](obol_session_handoff.md), replaced rather than appended |
| Why does a guard exist? | [Conformance evidence](libbobol_tla_conformance.md), [engineering lessons](libbobol_engineering_lessons.md), and dated archives |

The September 19 handoff, debt and readiness archives preserve previous runs
and counterexamples. They are evidence, not alternate current instructions.
Keep open tasks out of historical closure narratives. Do not create a second
backlog in a handoff or a text proposal. Further history cleanup is optional
unless it resolves an actual conflict or broken reference.

## Capability envelope

The existing required scope is retained. Modes and constrained outcomes are
inputs to the same owners, not separate orchestration systems.

| Dimension | Required coverage |
|---|---|
| Geometry | Small primitives, ordinary BoTs, spatial-page giant meshes, adaptive BREP, repeated and distinct assets, 50k/150k and real multi-gigabyte cases |
| Display policy | Modes 0--5 where applicable, LoD auto/off, authored/flat/smooth normals, mesh/point/box choices, whole-occurrence coverage |
| Resources | True-cold and certified-warm caches, ample/constrained memory, preparation/render deadlines, cancellation, compaction, shared-cache contention |
| View and execution | Quiet/active input, visibility turnover, single/quad panes, resize/fractional DPR, System GL/OSMesa, framebuffer and retained raytrace composition |
| User behavior | Draw, pick, selection/tree state, camera/faceplate, polygon/measurement, primitive/sketch editing and command/widget/manipulator agreement |
| Clients/platforms | qged, MGED, gsh, Archer and rtwizard at the declared support levels in the [platform matrix](libbobol_platform_threading_matrix.md) |

A stable input must reach requested quality or a witnessed constrained result
within declared resource bounds. Infinite changing input need not converge;
it must remain responsive and cancel obsolete work. A fast coarse image alone
does not pass the visual contract. Windows debt and native GPU qualification
are not closed by Linux llvmpipe; deferred macOS hosts remain deferred as the
platform matrix states. Scope changes require an explicit release decision.

## Roadmap and gates

S1--S6 identifiers are retained for existing references. S0 is an early build
prerequisite, rather than postponing dependency coherence to the final release.
Stages measure observable acceptance, not effort. Previous completion
percentages are retired: no defensible weighted denominator was established.
Report each gate as open, demonstrated on a named candidate, or qualified.

| Gate | Deliverable | Exit criterion |
|---|---|---|
| S0: reproducible stack | Recorded BRL-CAD/Obol/OSMesa source state, options and resolved runtime identities | A fresh isolated build and repeated configuration load the intended headers/libraries without manual repairs; smoke passes on the resulting binaries |
| S1: supported mutation and ownership inventory | Finite entry list and sole owners across host, scene, source/service, GED and presentation | Every live application entry has a caller, owner, lifetime, revision effect, callback/failure contract and acceptance check; unsupported paths have a migration or removal decision |
| S2: one complete production workflow | Cold draw, camera change during loading, cancel/close, reopen/redraw and terminal presentation | Actual production owners pass ordinary and delayed/stale/denied-work paths; geometry/images, resource release and convergence agree; replaced paths are deleted |
| S3: consolidate remaining owners | Apply S2's demonstrated method to the S1 inventory | No duplicate policy writers or mirrored control state; independent compilation and transition tests enforce boundaries; all required caller migrations complete |
| S4: resource and lifecycle completion | Typed bounded BREP work, source/service lifetime, cancellation, contention and terminal control | Each admitted operation terminates or cancels, releases its resources and rejects stale publication; unchanged evidence cannot reopen exhausted work; measured bounds pass |
| S5: capability qualification | Required drawing, geometry, performance, visual, editing, client and native-host rows | Every required row passes declared thresholds with retained reports and inspected images; no known blocker is hidden by another passing scenario |
| S6: release candidate | Exact final source/dependency/binary manifest and release decision | All required rows qualify that candidate; changed components trigger the affected checks; no unresolved blocker or unexplained missing evidence |

S1 and S2 are demonstrated on the September 20 checkout. The API contract has
a finite owner/lifetime/revision/failure table, production source/store callers
use those owners, and public/installed consumer plus transition checks cover
the migrated boundary. FLOW-01 exercises cold draw, camera input during real
delayed mesh work, worker-active close, reopen/redraw, changed terminal pixels,
compact mesh presentation and empty transient queues through the production GED
and Qt owners. Treat a newly demonstrated mutation family as an S1
counterexample and a failure of that workflow as an S2 counterexample;
otherwise continue with S3 rather than repeating either audit.

S2 evidence is under
`.build-main/obol-qualification/20260920-flow-01`. Its graphical path uses
offscreen Qt with software OSMesa and therefore does not qualify native GPU or
other platform rows. The source-evidence reducer and earlier counterexample
repairs remain useful production assets. Do not restart completed work merely
because later gates exercise the same owners at larger scale or on more hosts.

S3 proceeded as bounded seam closures. Its first closure makes the LoD
service the sole owner of complete transient-work quiescence and coherent
per-generation work observation. GED and GUI waits consume one lock-consistent
service snapshot instead of maintaining incomplete counter lists. Controller
and renderer readiness consumes one generation snapshot, including the
otherwise invisible interval in which a consumer waits through a shared-
producer lease. The service queue replaces the controller's former writable
result-pending mirror; the independent first-ready timestamp remains as
batch-age evidence. Phase-specific waits retain explicitly narrower
conditions. Evidence is under
`.build-main/obol-qualification/20260920-own-01-service-work`.

The second closure gives resident-capacity policy one service contract. Stable
renderer bytes are published directly instead of reconstructed from separately
changing total and backing counters. The service returns those stable bytes
with its reservation and limit, and exact growth publication precedes
reservation release. Headroom, pressure, convergence and qged consumers no
longer assemble that policy independently. A concurrent full-hierarchy test
changes the limit while a real growth reservation is visible and rejects a
reservation-to-stable gap. Evidence is under
`.build-main/obol-qualification/20260920-own-01-resident-capacity`.

The third closure makes convergence the sole source of progress-display
visibility, terminal readiness and stable display classes. Libged renders that
classification and Qt schedules it; neither rebuilds the other's policy. The
Qt-local class mappings and three split cache fields and libged's duplicate
visibility/readiness predicates are removed. Direct classification and
production faceplate, workflow, progressive, model and package checks are
retained under
`.build-main/obol-qualification/20260920-own-01-progress-display`. This is a
third bounded seam. It did not change pixels, so the existing annotation
controls remain the correct baseline.

The closing inventory sweep moved the remaining GED command and GUI
qualification waits to the existing host-work snapshot instead of sampling
controller flags independently. The apparent mirrors left by the sweep are
single-writer boundary transfers or diagnostic observations, not independent
policy owners. The affected build, focused transitions and repeated production
flows are retained under
`.build-main/obol-qualification/20260920-own-01-s3-close`. S3 is demonstrated
for the current inventory. A new duplicate writer or split policy is an S3
regression; otherwise the active roadmap proceeds through S4--S6.

S1--S3 define the simplification finish line. S4--S6 define production readiness.
Numerical policy and geometry fixes need not wait for every extraction, but
must retain the ownership rules and meet their own measurable release rows.

## How a change is accepted

Before implementation, name the backlog row, supported caller and failure to
fix. Record the current owner and the code/state that will disappear. A change
that only wraps existing owners or relocates interdependent policy does not
close an extraction gate. Line counts are review signals, not acceptance goals.

Use ordinary round trips and independent output checks alongside fault
injection. Check image pixels and geometry, not just cached pointers and fields.
Check memory magnitude, input latency, cancellation and stable convergence,
not just final ownership labels. Use production decision code in deterministic
tests; a separately implemented simulator cannot qualify the GUI's decisions.
Attribute image differences to camera, geometry, material or presentation.
When a control encodes intentionally replaced behavior, retain it and record
an isolated comparison plus independent assertions of the replacement contract
before establishing a new control. Updating images alone does not close a row.
Re-run relevant models for ownership/liveness changes, and renderer/geometry
checks for numeric changes. TLA+ covers its finite abstract protocol; it does
not prove the C++ implementation, visual fidelity or wall-clock bounds.

Each closure records its code boundary, removed state/writers, source and
runtime identities, applicable tests, output evidence and remaining limits.
Keep legacy regression protections until their callers migrate. Reentry that
happens to succeed is not a reason to expand the supported callback contract.
No new scheduler, generic effect vocabulary, scene-sized mirror, per-occurrence
object graph or compatibility layer without a demonstrated requirement.

## When to replace a boundary

If the S2 workflow still needs multiple owners to establish one fact, cannot
state when work ends, or adds more exception paths and mirrored state after
caller migration, stop extending that boundary. Replace its orchestration
behind the existing data/execution interfaces and run the same acceptance
workflow. This is the decision point for selective source/service or
publication replacement. A whole-stack rewrite requires evidence that these
bounded replacements cannot retain the proven domain implementation.
