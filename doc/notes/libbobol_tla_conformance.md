# libBObol TLA+ to C++ conformance audit

Last reviewed: 2026-09-12

This is the implementation-facing conformance ledger for the BObol formal
suite.  `tla/models.json` remains the authoritative proof graph and test
catalog; this document records whether the production owners actually refine
the modeled transitions.  A TLC pass proves the bounded abstract relation,
not the C++ mapping.  A row is conforming only when the formula, production
owner, and executable witness agree.

The audit began with the failure families in
`libbobol_engineering_lessons.md`, then considered adjacent races rather than
limiting its scope to failures already observed.  In particular, address and
dense-slot reuse motivated independent identity domains; stale demand failures
motivated authentication before outcome interpretation; lost wakeups motivated
level-triggered owners; and mutable shared generations motivated explicit
consumer lifetime rather than reliance on pointer ownership.

## Status vocabulary

- **Conforming**: the model-to-C++ mapping is explicit and an executable test
  exercises the refinement boundary.
- **Partial**: the core owner exists, but an identified transition, trace, or
  failure edge lacks executable evidence.
- **Nonconforming**: the current C++ ownership permits a behavior forbidden by
  the model.
- **Outside abstraction**: retain a separate data-plane, fault-injection,
  sanitizer, performance, image, or platform gate.

## Shared terminology

| Formal term | Production meaning | Authoritative C++ representation |
|---|---|---|
| evidence stamp | all semantic inputs which authorize one admission plan | `BObolLodAdmissionRevisionStamp`: inventory, availability, visibility, view, policy, capacity |
| demand | the view-local request which authorizes presentation | `BObolViewEpoch` plus `BObolPolicyEpoch`; neither is asset identity |
| source route | lifetime of the source owner to which a result returns | `BObolSourceRoutingId` |
| source population | lifetime of the compact dense registry addressed by an entry index | `BObolSourcePopulationEpoch` |
| current result | route, population, view, and policy all match their current owners | `BObolLodResultAuthenticationContract` |
| publish | install a useful current result | `BObolLodResultDisposition::PUBLISH` followed by `BObolViewLodState::consumeSourceResult` |
| terminal failure | exact-current provider/cache failure evidence | `RECORD_TERMINAL_FAILURE` and the demand-scoped CAD occurrence failure map |
| retry | no failure evidence; recreate current demand | `RETRY_CURRENT_DEMAND` and the controller demand replay edge |
| supersede | discard an identity-mismatched result before interpreting its outcome | `SUPERSEDE`; no payload or failure-map mutation |
| work owner | the one enabled finite mechanism which can discharge an obligation | `BObolLodControlRefinement::Owner` |
| consumer lease | one live view/generation's claim on an in-flight shared producer or queued result | `BObolSharedProducer::leases`, keyed by service generation and carrying that consumer's newest demand |
| authentication identity | a fixed-width credential whose equality authorizes work, publication, or ownership | checked strong values and `identity_counter_private.h`; exhaustion is fail-stop and never wraps |
| diagnostic counter | an observation total which authorizes no behavior | saturating counter; it may stop increasing but may not be used as an identity |

An occurrence-targeted result must carry nonzero route and population
credentials even when a focused test uses semantic-key routing with
`sourceEntryIndex == UINT32_MAX`.  Empty credentials are not a wildcard.  A
source-wide result has an empty occurrence key and is outside compact
occurrence authentication.

## Conformance by proof boundary

| Proof boundary | Production enforcement | Executable evidence | Status |
|---|---|---|---|
| `ObolProgressivePipeline`: one six-domain stamp, one bounded cursor, one finite owner, derived terminal outcome | `lod_revision_private.h`, `lod_control_private.h`, coordinator evidence values and controller effect executors; exact Obol GPU keys; exact/interned capacity, structural-frontier, and compact-presentation credentials | `test_lod_coordinator.cpp` six-domain revision, owner, witness, typed-outcome, equivalent-population, and collision tests; retained-allocation identity tests; randomized completed-pass and transition-journal tests; focused Obol renderer tests | **Partial**: the source-evidence boundary is conforming after `CXX-SOURCE-001`; the remaining S1 writer/acceptance audit is open. The `CXX-STAMP-001` value boundary remains closed: partial and caller-assembled planning stamps are unrepresentable, every domain independently invalidates stored plans and allocation keys, and mutable-state cache authorization retains exact inputs rather than a folded digest. |
| `ObolControlRefinement`: 29 independent facts refine to nine work classes and one fixed-precedence owner | `BObolLodControlRefinement`, `BObolLodControlTransitionScope`, controller transition journal, and runtime violation projection | per-field and exhaustive owner combinations in `test_lod_coordinator.cpp`; randomized journal test; complete graphical trace checker; result-to-presentation trace oracle | **Conforming** after `CXX-TRACE-001` and `CXX-VTRACE-001`: every retained transition has typed endpoints, the audit sentinel rejects an unmapped writer, and authenticated result publication remains causally paired with its exact completed frame. |
| `ObolResultAuthentication`: independent route, population, and demand authentication before publish/failure/retry | service result normalization, `lod_result_authentication_private.h`, controller publication gate, view-state reducer | direct disposition table plus the 50-case asynchronous controller matrix in `test_lod_update_action.cpp`; service request-matching tests | **Conforming** after `CXX-RP-001`. |
| `ObolIdentityExhaustion`: a fixed-width identity advances without reuse or fail-stops before a semantic mutation | `BObolLodStrongUInt64`, libBObol `identity_counter_private.h`, and Obol `CadIdentityCounter.h` across retained-plan, classifier, preparation, timing, and resource evidence | successor boundary tests in `test_lod_coordinator.cpp` and Obol `test_cad_instance_records.cpp`; cross-session trace-token regression in `test_window_host.cpp` | **Conforming** after `CXX-ABA-001`. |
| `ObolSharedAssetLease`: a surviving or late consumer retains shared build/result ownership; only the final lease cancels | typed shared-producer state, generation leases, producer-aware cancellation, and per-generation payload/replay drains in `BObolLodService` | same-generation coalescing plus six cross-generation build/result/cancellation lifecycle cases in `test_lod_service.cpp` | **Conforming** after `CXX-SL-001`. |
| admission/allocation/capacity and terminal-quality component models | typed coordinator policies and retained-allocation transaction | focused exhaustive/property tests and renderer evidence tests | **Conforming** for modeled control relations; numeric/visual sufficiency remains outside TLA+. |
| retained CAD mutation and exact-frame compositions | private staging, checked Obol mutation, exact target/report revisions | Obol mutation tests, CAD presentation/frame tests, controller endpoint and commit-denial tests | **Conforming** for atomic precommit resource denial after `CXX-TXN-001` and typed endpoint cancellation after `CXX-LIFE-001`. |
| host, policy, interaction, and semantic-presentation lifecycle | typed host request, interaction session, exact-frame and policy retirement ledgers | off/on/off, interrupted frame, callback/unsubscribe, endpoint-loss, and semantic-only presentation tests | **Conforming** after `CXX-LIFE-001` and `CXX-POLICY-001`: endpoint loss retires controller and host owners, callback removal drains or defers every live dispatch, and live policy publication advances its control epoch within the same input operation. |
| cache publication, durability, geometry truth, renderer pixels, resource magnitude | staged payload records, atomic name mapping, provider validation, immutable geometry, renderer implementations | denied replacement, corrupt-cache, geometry, image, GUI, sanitizer, and performance matrices | **Conforming** for modeled atomic name publication; crash consistency, numeric truth, and resource magnitude remain outside the abstraction as noted in `tla/RISK_COVERAGE.md`. |

## Closed vertical finding: result publication

### CXX-RP-001 — stale result could bypass complete authentication

Before this pass, compact route and population were checked separately in the
view state, while view/policy arbitration in the controller allowed an old
mesh to publish as a cold bootstrap when no equal or newer mesh was resident.
That violated `ObolResultAuthentication.StaleResultCannotApply`.  It also
interpreted some provider outcomes before proving identity, so an obsolete
terminal error could become failure evidence for a replacement occurrence.

The production boundary now performs one allocation-free decision over three
independent identity domains.  Any mismatch yields `SUPERSEDE` before provider
status is interpreted, causes no payload or failure-map mutation, and arms one
current-demand replay.  Current `READY` publishes; current `CACHE_MISS`,
`STALE`, and `ERROR` record exact-demand terminal failure; other provider
statuses retry without failure evidence.

The cumulative progressive-mesh rebase remains legal only as a
pre-authentication refinement.  It must prove identical immutable asset
ownership and reconstruct the current page set, counts, and prepared renderer
generation.  Only then may it rewrite the request to current demand and enter
the publisher.  It is not permission to apply a stale result.

`ObolResultAuthentication` now includes the same publish, terminal-failure,
and retry dispositions.  TLC checks that a stale result enables neither
publication nor failure recording and cannot create an accepted outcome.  The
configured finite run explores one bounded retry because production retry
count is deliberately not a formal resource claim.

The C++ table tests every provider status.  The asynchronous matrix crosses
all ten statuses with current identity and independently stale route,
population, view, and policy identities.  Stale rows run while the occurrence
is cold, so a future bootstrap exception cannot silently restore the defect.
The stricter contract also corrected hand-built fixtures which had omitted
route/population credentials that production requests always carry.

## Closed vertical finding: shared producer lifetime

### CXX-SL-001 — coalescing is not a consumer lease

Before this pass, `lod_request_active_key` correctly removed occurrence identity when
`coalesceAssetProducer` is true, and it preserves route and compact-population
identity.  `activeRequestKeyCounts` and `latestActiveRequests` then serialize
the producer and retain its newest demand.  This is request coalescing, not the
multi-consumer lifetime required by `ObolSharedAssetLease`:

- a duplicate submission returns no task ID and records no consumer or
  generation lease;
- the surviving work item and completed result retain only the first task's
  generation;
- `cancelGeneration` cancels that work and erases that generation's queued
  result without consulting another generation which requested the same
  asset; and
- `drainGenerationResults` can deliver the result only to its original
  generation.

Consequently, closing the first view could cancel a build still needed by a
second view, and a late join has no durable claim on a queued result.  The
post-publication `residentMeshConsumerDemands` snapshots protect compaction;
they do not lease an active producer and cannot satisfy this model.

The service now owns one typed producer record per stable active asset key and
one lease per service generation.  Duplicate discovery in both submission and
the earlier demand-update fast path records the lease and that consumer's
newest demand.  Cancelling a generation removes only its lease; worker,
preview, debug-delay, and cache cancellation consult the producer lifetime and
continue while any lease survives.  Pending work and undelivered results
retire only with the final lease.

Only the producer generation receives its demand-specific payload.  Other
generations receive a lightweight, exact-current `SUPERSEDED` replay result,
which the result-authentication reducer maps to `RETRY_CURRENT_DEMAND`.  This
avoids copying large mesh arrays or publishing geometry prepared under another
view's admission certificate; the awakened consumer instead binds the shared
resident asset under its own current demand.  Late joins work during both the
build and queued-result states.

Six executable lifecycle cases cover two generations sharing one delayed
build, first-generation cancellation while the second survives, late join
during build, late join while a result is queued, final-consumer cancellation
during build, and final-consumer cancellation with a queued result.  They
assert one provider call, generation-scoped delivery or replay, and complete
retirement of producer, lease, request, result, and in-flight ownership.
`ObolSharedAssetLease.ResultRetainedForUnsatisfiedConsumer` makes the queued
result rule explicit.  The affected TLC closure passed with baseline counts:
237 generated / 96 distinct lease states, 7,067 / 1,900 asset-publication
states, 712,617 / 356,352 control-refinement states, and 10,167,666 /
4,303,828 canonical pipeline states.

## Closed vertical finding: complete control-transition refinement

### CXX-TRACE-001 — sampled states were not an execution trace

The prior qged observer sampled distinct convergence states at event-loop
intervals.  Its projection and cycle checks were useful, but several reducer
effects could occur between observations.  A passing report therefore proved
only that the states it happened to see were valid; it could not prove that
each concrete transition belonged to the finite refinement relation.

The controller now owns an opt-in bounded journal.  Every retained record has
an immutable before and after endpoint, a contiguous serial, and one event
from the external-input plus nine-owner alphabet.  Owner-thread effect
executors use nested `BObolLodControlTransitionScope` boundaries.  Cohesive
external-input setters suppress their implementation-detail nested scopes so
revision, obligation, and successor wake publication remain one atomic
transition.  Worker result callbacks publish only a typed publication event;
the owner captures its endpoints without reading Coin state from the worker.

The endpoint includes the complete control fact/obligation/owner projection,
all admission revisions, plan and presentation identities, host work, and the
exact renderer-preparation target signature and remaining-unit rank.  When
observable state changes between registered scopes, the journal emits
`unnamed`; it never guesses an owner.  Full buffers report drops rather than
overwriting history.  Neither the journal nor its checker is consulted by
production policy.

The qged checker now rejects an unnamed or unknown event, noncontiguous
serial, missing endpoint field, discontinuous before/after pair, event/owner
mismatch, rank regression, local invariant, and trace truncation.  Its
adversarial shell fixtures cover each new rejection.  A 512-step deterministic
controller trace checks boundary nesting, continuity, and bounded-loss
reporting in C++.

The first complete OSMesa draw/settle/zoom trace exposed a real missing
refinement witness: consuming a queued render cleared `RENDER` before the
exact-frame reducer completed.  `ObolHostWork` already models that interval as
`renderInFlight`; production now publishes the matching typed claimed-frame
host level and presentation witness and retires it only at completed or
interrupted-frame reduction.  The repaired representative trace contains 57
transitions, 15 claimed-frame endpoints, no unnamed events or drops, and
passes the independent refinement checker.

## Closed vertical finding: result-to-presentation trace refinement

### CXX-VTRACE-001 — valid snapshots did not prove publication reached a frame

The complete transition journal closed the gaps between sampled control
states, but the generic checker still treated transaction and presentation
serials as diagnostic clocks.  It did not prove the vertical causal chain from
an asynchronous result-ready notification, through owner-thread publication
and an exact-frame barrier, to retirement by the required completed CAD
render.  In particular, a trace made entirely of locally valid snapshots could
still lose, replace, or retire that barrier at the wrong frame.

`BObolResultPresentationTraceChecker` is a test-side observer over the public
journal.  It rejects missing or dropped prefixes, unnamed events,
noncontiguous serials, discontinuous endpoints, invalid production states,
event/owner mismatches, revision or clock regressions, and disagreement
between the concrete 29-fact mask and an independently recomputed obligation
and owner.  Once a current result arms a publication barrier, the checker
retains its transaction, view, policy, and required-render identities.  Only a
matching completed-render clock may retire it; a typed external-input revision
may instead supersede it.  The exact render-completion clock is now part of the
public trace endpoint and qged's JSON endpoint schema, distinct from the
broader host presentation counter.

The production integration witness submits delayed coarse and refined results
through the real worker service, publishes each immutable payload on the
controller thread, completes both render transactions, drains the bounded
journal, and requires at least one fully retired publication barrier.  Small
synthetic traces independently prove that the oracle rejects lost serials,
discontinuous endpoints, unnamed and wrong-owner events, missing completed
frames, and traces with no vertical flow.

Running that witness exposed three adjacent defects.  Final non-partial result
batches requested a frame but, unlike partial prefixes, did not arm its numeric
publication barrier.  A manually submitted generation was also manufacturing
an automatic demand-refresh successor after result delivery even though
automatic LoD was disabled.  Conversely, the disabled convergence shortcut
reported `terminal` while the manual result's publication or presentation
reducer still owned foreground work.  Every requested result frame now has an
exact successor identity, manual delivery owns only its submitted request, and
the disabled shortcut remains active until any explicit publication/control
debt reaches its frame.  The ordinary no-LoD path is unchanged and immediately
ready.

No formula change was needed.  With `automatic = FALSE`,
`ObolControlLifecycleComposition.CompleteProvider` already creates no demand
pass but does create exact-frame debt; `QueueExactFrame`, `BeginRender`, and
`CompleteCurrentRender` carry that debt to a current presentation.  The new
conformance-catalog entry maps the live C++ trace to exactly those actions.

## Closed vertical finding: complete evidence stamps

### CXX-STAMP-001 — admission stamps permitted partial construction

The canonical TLA+ model already represents both current and in-progress
planning identity as total six-field records.  C++ did not enforce that
shape: `BObolLodAdmissionRevisionStamp` was a default-constructible aggregate,
and its coordinator snapshot and tests assigned companion fields one at a
time.  A future asynchronous constructor could therefore omit visibility or
capture fields on opposite sides of a semantic mutation while still passing
the type checker.

The stamp is now an immutable, allocation-free value requiring typed
inventory, availability, visibility, view, policy, and capacity arguments.
The only empty value is the named `administrative()` sentinel used by
unstamped immediate reducer results and inactive certificates.  Stored plan,
cursor, allocation, constraint, and acceptance values initialize that
sentinel explicitly; asynchronous inputs and results can replace it only by
copying a complete stamp.  The owner-thread coordinator is the sole production
constructor from live revision owners and captures all six in one expression.

All six mutation domains now have executable stale-identity checks for the
admission plan, its resumable cursor, and the retained-allocation input key.
Compile-time checks reject a default stamp and a five-domain constructor while
preserving trivial-copy and fixed-size behavior.  The audit also found three
intentional partial comparisons for renderer-capacity evidence.  They now use
one named `sameRendererCapacityProblem` projection: inventory, availability,
view, and policy define that measured problem; visibility is a replaceable
allocation input and capacity is the measurement's output.  This distinction
prevents an ad hoc partial comparison from masquerading as full plan
authentication.

No formula change was needed: `ObolProgressivePipeline` already requires the
complete stamp at plan start, each planning step, completion, and stale abort.
The affected coordinator, retained-allocation, host, and trace tests pass.

The cross-repository audit subsequently found lower-level Obol GPU caches which
folded several revisions, progressive ranges, or flattened-atlas credentials
into a 64-bit digest and then treated equality as proof that an old upload was
current.  That conflicts with the model even though collision is unlikely:
probability is not an authorization invariant.  Aggregate-proxy uploads now
retain the exact source, attribute, presentation, draw-mode, and software-path
tuple.  Progressive cut buffers retain their exact representation, lineage,
interval, quantization, and packed range vector.  Flattened wire and shaded
atlas entries retain instance, part, cut, source revision, transform, and
source interval.  Hashing remains legal only as a lookup accelerator followed
by exact equality.

The same rule applies above the renderer.  Retained-allocation capacity reuse
now carries a non-reused population identity issued only after the view owner
compares the exact ordered selected population.  Equal diagnostic digests do
not authorize reuse.  Structural-frontier deduplication uses an exact bounded
interner; eviction can cause only a redundant rebuild, never false reuse.
Compact-source and cross-source presentation caches compare exact revision and
membership vectors, including semantic revision, rather than structure/style
digests.  Length-prefixed layered payload keys and exact stable-handle
interning close the corresponding variable-string identity paths.

Completed-frame aggregation follows the same rule.  Obol assigns every
`SoCADAssembly` a process-lifetime identity, including across address reuse.
libBObol canonicalizes the exact `(assembly identity, local serial)` set and
issues a checked non-reused token for draw execution, GPU timing, and GPU
resource observations.  Invalid or absent evidence retires the current set, so
an `A -> absent -> A` source sequence cannot reuse its earlier aggregate token.
Obol likewise replaces the saturating sum of its assembly- and renderer-owned
preparation serials with an exact-pair change identity.

No TLA+ formula changed for these closures.  The formulas already compare
abstract records and identities exactly; the C++ digest shortcuts were invalid
refinements of that relation.  Strong 128-bit immutable geometry identifiers
and persistent mesh content hashes remain deliberate content-address domains,
not mutable evidence stamps.  Their byte-level validity and cache durability
remain executable/operational obligations rather than claims of the bounded
control models.

## Closed vertical finding: transactional publication

### CXX-TXN-001 — denial could notify or erase before commit

The retained-scene path had pure staging and checked mutation calls, but an
outer `SoCADAssembly` update scope was opened before the final resource gate.
An injected denial changed no scene bytes yet still closed that empty scope and
advanced the node notification identity.  Resource admission now occurs after
all compact-presentation preflight and before the first update scope.  The
executable fixture proves that both retained-scene and presentation denial
preserve the prior scene and `SoNode` notification ID, while the succeeding
commit changes both together.

The durable mesh-LoD path had the inverse problem.  Forced refresh deleted the
old name-to-content mapping before candidate generation.  A write, allocation,
or disk failure could therefore leave a complete immutable predecessor on disk
but make it undiscoverable.  Refresh and authored-mesh replacement now build
all content-addressed records first and swap the name mapping only after the
candidate is complete.  Deterministic failure points abort component and
spatial transactions immediately before commit and reject the final mapping
write.  Tests cover both no-prior-object denial and failed replacement of a
known-good hierarchy; the latter must retain the same key and reopenable
payload.  Optional in-process name hints are cleared on post-commit allocation
failure so they cannot mask the authoritative disk mapping.

This finding required a formula correction.  The former
`ObolLiveSpatialPublication` and asset composition represented durability with
only a Boolean marker, so they could not express replacement of an existing
mapping—the exact state lost by the C++ bug.  Both models now carry a baseline
mapping, an atomic candidate publication, and `DenyCacheCommit`.
`DeniedCacheCommitPreservesMapping` proves that denial leaves the baseline
discoverable, while `CacheMarkerRequiresFinalHierarchy` still prevents early
publication.  TLC checks the denial action as required coverage in both the
focused model and its parent composition.

The guarantee is deliberately at the precommit boundary.  Complete-scene
replacement uses a staged no-throw swap.  Sparse retained mutation reserves
for its bounded journal but does not clone a 150k scene to recover from
process-wide allocator exhaustion during the mechanical commit.  Filesystem
crash consistency and kill-at-instruction durability remain separate
operational tests rather than claims of this bounded control model.

## Closed vertical finding: endpoint and callback lifetime

### CXX-LIFE-001 — endpoint loss left work without a host

The display endpoint previously detached its callback and host without
retiring work retained by a borrowed controller.  An open gesture could remain
in its interaction/debounce protocol, a consumed render request could remain a
claimed-frame witness, and worker, coverage, presentation, or capacity owners
could remain live even though no endpoint existed to advance them.  This was
not an orderly producer close: the remaining frame and timer fairness premises
were false after the host disappeared.

Endpoint loss is now one explicit cancellation transaction.  Before clearing
the host levels it resets renderer-derived presentation controls, cancels the
active service generation, retires every automatic-LoD domain, closes the
interaction session, invalidates capacity evidence, and clears pending and
claimed frames plus the shared pump.  The postcondition asserts both the
complete automatic-domain retirement predicate and an empty typed host-work
snapshot.  Immutable resident payloads and user policy survive; binding a new
endpoint starts a fresh coverage/capacity transaction and reprojects any
independent provider work.

The same audit found that presentation callback removal copied a raw callback
and user pointer under a mutex, then invoked it without a dispatch lease.
Removal could return while the callback still used freed endpoint state.
Frame callbacks had a dispatch count, but self-unsubscribe waited for its own
dispatch and deadlocked.  Both callback kinds now use the same reference-
counted close protocol.  External removal waits for all dispatches; removal
from the current or any enclosing nested callback defers deletion until that
dispatch unwinds.  An allocation failure preserves the prior live callback
rather than silently dropping it.

Executable cases close an endpoint around a borrowed controller with an open
gesture, active service generation, progressive work, and claimed capacity
frame, and require zero surviving host flags, control facts, obligations,
owner, generation, or capacity search.  Separate cases prove external
presentation teardown blocks until callback return, self-unsubscribe is safe,
and a nested callback can unsubscribe its outer dispatch without deadlock.
The pre-existing factory race still proves host destruction waits for an
active frame callback.

`ObolInteractionSession.CloseInput` already modeled close during a gesture.
`ObolHostWork` now distinguishes endpoint loss from drainable
`CloseProducers`: `CloseEndpoint` atomically cancels notification, loop, pump,
pending render, claimed render, exact-frame, publication-timer, source, and
provider ownership.  It is required nonzero action coverage and exports
`EndpointClosureRetiresAllWork` to the lifecycle composition contract.  The
focused TLC run passes at 89,521 generated / 12,924 distinct states to depth
17.

## Closed vertical finding: fixed-width identity exhaustion

### CXX-ABA-001 — counter wrap could authenticate retired evidence

The formal revision domains were monotonic naturals, but several production
counters were fixed-width integers.  Most skipped zero after unsigned wrap;
`compact_next_revision` explicitly returned to one.  A sufficiently old
revision, generation, route, handle, lineage, plan serial, or callback token
could therefore compare equal to a new credential and authenticate retired
work.  Skipping the zero sentinel reduced the collision set by one value; it
did not prevent ABA reuse.

All increment sites which create authentication, ordering, or ownership
credentials now use the same checked-successor contract in their owning
repository.  It covers libBObol local and atomic counters, zero-valid
sequences, nonzero identities, public strong LoD epoch values, and Obol's CAD
plan, classifier, preparation, GPU-resource, and timing evidence.  A current
identity may advance to the maximum value; a stored next-identity allocator
reserves that value as its exhaustion marker.  In either representation, an
operation which needs an unavailable successor terminates before it can
commit the associated semantic mutation.  This fail-stop boundary is
intentional: most affected APIs cannot roll back every already-observed
dependent object, and continuing after exhaustion would be less safe than
stopping.  Obol resource teardown no longer resets externally observed upload
or timer allocation domains, so context recreation cannot manufacture an ABA
without numeric overflow.  Diagnostic totals are explicitly separate and
saturate because they authorize no behavior; resource byte/count arithmetic
remains governed by its own saturation and reservation tests.

Aggregate identities are checked at the same boundary.  An assembly object has
a process-lifetime identity rather than relying on its reusable address, and
multi-assembly execution, GPU timing, and resource samples are interned from
their exact canonical source tuples.  Obol's combined preparation token is
also derived from exact component equality rather than a saturating sum.

The audit also found an ABA path without numeric overflow.  Disabling and
re-enabling controller transition tracing reset its token allocator while an
old RAII scope could still exist.  A late destructor could then carry the same
token as a new scope and close the wrong transition frame.  Trace records may
be cleared between sessions, but their token and serial domains now remain
monotonic for the controller lifetime.  A regression leaves an old scope
alive across re-enable and proves that it is rejected rather than accepted as
the current scope.

`ObolIdentityExhaustion` closes the TLA+ abstraction gap with a finite revision
domain.  A requested mutation either atomically publishes a fresh successor
or takes `FailStop` at the boundary; issued revisions are never reused, stale
evidence cannot authenticate, and the halted state is closed and quiescent.
The focused TLC run explores the complete configured graph at 91 generated /
37 distinct states to depth 8.  Pure C++ boundary tests in both repositories
cover zero, the final valid successor, exhausted detection, local/atomic
allocation, and diagnostic saturation without deliberately terminating the
test process.

## Operational model-to-implementation gate

`tla/conformance.json` is the machine-readable bridge from formal actions to
C++ evidence.  It defines one shared observation vector for control ownership,
consumer generation, exact request identity, authenticated result disposition,
terminal lifecycle, semantic presentation debt, render-barrier transaction,
source identity, and evidence validity.  The C++ harness copies those values
from production reducers; it is an observer and never becomes a second
scheduler or source of policy.

The same harness adapts the public `BObolLodConvergenceStatus` and
`BObolHostWorkSnapshot` pair into that vector.  Its public-observation scenario
checks render request/retirement and worker-generation attachment/detachment,
so lifecycle and generation fields cannot remain catalog-only placeholders.

`test_model_conformance.cpp` executes the highest-risk reducer sequences one
formal action at a time and checks the common invariants after every action.
Its deterministic scenarios cover owner discharge; current, failed, retrying,
and late superseded results; exact-set reordering, invalidation, ABA
non-reuse, and numeric exhaustion; coherent retained-replay timing; and early,
incomplete, current, and stale presentation frames.  The semantic-frame steps
use formal action names, while implementation-only barrier steps are identified
as such.  The trace-refinement scenario also mutates a valid vertical witness
to prove each structural rejection, while the asynchronous staged-result test
checks the same oracle against production worker and renderer transitions.
Existing asynchronous tests remain the witnesses for worker
cancellation, callback teardown/reuse, shared producer leases, and full
controller behavior.  The manifest connects both kinds of test to every action
singled out by `models.json.semanticAudit.requiredActions`.

`ValidateBObolConformance.cmake` fails if an observation domain disappears, a
model or action is renamed without updating its C++ scenario, a test source is
missing, a stepwise scenario stops executing its formal action name, or any
required semantic-audit action loses its mapping.  It runs as a focused CTest
and from the TLA+ runner, so the bridge is checked by both implementation and
proof changes.

The runtime gates are deliberately separate:

- `ctest -L bobol_conformance` is the fast ordinary-CI contract;
- `bobol_sanitizer_tests` plus `ctest -L bobol_sanitizer` is the shared
  ASan/UBSan contract and the nightly TSan contract;
- `bobol_graphics_offscreen`, `bobol_graphics_osmesa`, and
  `bobol_graphics_system_gl` keep software/offscreen, OSMesa, and native
  X11/System-GL evidence distinct; and
- `bobol_performance_qualification` records cold quick-profile evidence and
  fails if no database or expected iteration result was actually exercised.

The durable-cache fault matrix now denies opening the backing store, writing
staged payload data, and committing the authoritative name mapping.  Initial
publication remains undiscoverable at every edge, and replacement failure
preserves the prior complete hierarchy.  This checks transactional visibility;
power-loss durability and storage-device behavior remain release/platform
qualification rather than a TLA+ claim.

## S1 ownership inventory: September 5 first pass

This is evidence for the [simplification guide](obol_simplification_guide.md),
not another task backlog.  The audit inspected current declarations and writer
sites in the coordinator/controller, host adapter, database source, source
realization service, LoD service, GED draw adapter and presentation staging.
It also checked the existing tests at those boundaries.  **S1 remains open:**
the table classifies state families, while the remaining per-transition audit
and acceptance-threshold closure remain in active debt.

Physical ownership and formal conformance are separate questions.  A typed
value can enforce valid transitions while several controller callbacks still
choose which transition to invoke.  Conversely, a declaration/reference count
does not establish whether state is redundant.  The low-reference controller
fields inspected here include real identity counters, policy values and
diagnostic outputs; none was justified for deletion merely by that count.

| State family | Current production authority and lifecycle | Evidence and disposition |
|---|---|---|
| Six admission revisions, admission evidence and bounded cursor | `BObolLodCoordinator::advanceAdmissionRevision` and the admission planner; a semantic edge retires a stale cursor, and commit checks the complete stamp | Revision, admission and allocation-oracle tests already exist. Preserve the typed values; source observation/consumption now delegates to `BObolLodSourceEvidence` as detailed below. |
| Observed and submitted source evidence | `BObolLodSourceEvidence::observe` is the sole observation writer; exact source identity plus inventory/visibility journals identify its input. `consume` accepts only the latest observation; service/generation/empty-source retirement resets both handles | `test_lod_source_evidence.cpp` and the real-controller visibility fixture cover publication, stale consumption and retirement. At most two immutable source-signature vectors are retained, sharing after consumption; no occurrence/geometry copy. Host delivery and broader submission effects remain in the other audited families. |
| Submission source plan, visible count, delta and pass | Coordinator cursor/source operations plus `submitLodRequestsIfNeeded`/`submitLodRequests`; source-local progress completes or is superseded, with sparse deltas retaining consumed prefixes | The [writer/successor inventory](#submission-writer-and-successor-inventory) classifies all 51 direct mutation sites and their indirect entries. Production regressions cover input retargeting and demand preservation through presentation repair. The inventory identifies remaining source-supersession proof and physical extraction work. |
| Coverage, visibility census, projected demand and retained continuity | Coverage/census/demand values in the coordinator; source/view changes invalidate the relevant proof, exact deltas update the census, and finite scans establish coverage | Existing coverage, visibility, projected-demand and continuity tests. These are necessary view data and certificates, not duplicate scene ownership. |
| Availability and planning obligations | `BObolLodAvailabilityLedger`, `BObolLodPlanningObligations` and the existing work-ledger projection; result arrival, inventory completion and resource release enable named work | Availability/scheduler, planning-obligation, complete-trace and composed-lifecycle tests exist. Finish auditing their caller-selected successors; do not add a scheduler above them. |
| Static quality | `BObolLodStaticQualityTrial::completeFrame` now selects completed-frame replay, acceptance, rejection, cut advancement and population handoff in independently compiled `lod_static_quality.cpp`. The coordinator applies admission effects; the render callback applies presentation effects | The [static-quality inventory](#static-quality-writer-inventory) records the remaining lifecycle/start/interrupt writers. Production-policy and real two-occurrence cost-provider tests cover `CXX-STATIC-001`; this bounded extraction does not close all controller ownership or graphical qualification. |
| Presentation limits and timing | `BObolViewController::Impl` publication methods write ceiling/point thresholds; presentation/history/deadline values own restoration. Frame completion and coordinator reset/seed operations update throughput evidence | Quality, presentation/history and deadline tests exist. Preserve canonical cost versus actual renderer work. Remaining caller-side numeric decisions belong in their existing policies; timing counters must not become retry identities. |
| Exact frame, interrupted replay, publication and point-quality work | Existing typed exact-frame/replay/presentation/point owners; completion or cancellation retires a specific target, and a retry needs finite preparation progress | `test_cad_presentation_frame_retirement`, exact/interrupted-frame tests and transaction tests. Preserve last valid pixels and bounded preparation. Real cancellation latency and the retained frame-delivery failure remain release work. |
| Host request, claimed frame and capacity claim | `BObolRenderRequestState` plus `consumeRenderRequest`, `retireClaimedRender`, and endpoint retirement under `renderRequestMutex` | Host-work and endpoint-loss tests cover the effect boundary. The two claimed flags have only these writer sites, but still encode a three-way value as two Booleans; include them in S3's cohesive host-state cleanup, not a new phase owner. |
| Source identity, compact registry and journals | `BObolCadSourceState` and `BObolCompactOccurrenceRegistryState`; routing survives source lifetime, population identifies dense storage, and separate inventory/visibility journals carry sparse changes | The [compact-adoption audit](#compact-source-adoption-boundary) closes `CXX-SOURCE-005` and records installation effects. Complete remaining source setters and sparse writers; declaration grouping alone does not prove sole ownership. |
| Detached source production | `BObolSourceRealizationJobPrivate` and coordinator worker queue; item state, cancellation and remaining count survive to terminal adoption/reaping; byte reservations and bounded bypasses control admission | The [source submission audit](#detached-source-submission-and-release-inventory) classifies lifecycle/queue/reservation writers and closes `CXX-SOURCE-002` allocation failure after ownership transfer. Source registry/adoption, pressure/shutdown and wider service audits remain open. |
| Shared producer, queued result and residency | `BObolLodServicePrivate` under its mutex, shared generation leases, resident-demand revisions, compaction targets and growth reservation scope | Shared-producer lifecycle, cancellation, queue-limit, working-set and resident-compaction tests exist. S3 separates task/lease, residency and persistence bodies while preserving these owners; full worker-writer and release-edge audit is still required. |
| GED semantics and deferred effects | Semantic reducer requests are distinct from `ged_obol_progressive_provider_data`, deferred source jobs and command-time autoview data. Jobs delegate production to `BObolSourceRealizationCoordinator`; autoview/refinement flags are still written in several `draw_obol.cpp` routines | Existing GED/deferred-autoview models and client tests are foundations. S3 must separate semantic reduction from deferred realization/autoview effect orchestration. Source appearance snapshots are data, not a second editable authority. |
| Retained presentation staging | `BObolCompactPresentationStaging` owns a temporary change journal; checked Obol mutation accepts the candidate before `commit` publishes the bridge state | Existing mutation-denial and presentation tests. Preserve the independently compiled staging boundary and bounded delta storage; it needs qualification, not another facade. |
| Configuration, caches and observations | User options/deadlines/forced-cut settings, camera/lighting state, picking caches, diagnostic timings/counters and the opt-in trace journal | Classify by consumer. These fields are not additional progressive phases. A field used to authorize work must instead appear in the relevant certificate/obligation row. |

The remaining audit is bounded by these families: complete caller/alias writer
coverage, source/adoption and service release edges, GED deferred effects and
submission/presentation integration, then reconcile every required acceptance
row with a declared threshold and witness.  Do not restart the already mapped
formal hierarchy.  The inventory has already identified a concrete production
counterexample rather than only a file-size concern.

## Compact source adoption boundary

This S1 pass covers detached adoption, index installation and their identity,
journal and resource-retirement effects in `database_source.cpp`.  It does not
close all source setters or sparse mutation writers.

| Fact/resource | Authority and lifetime at adoption |
|---|---|
| Source result identity | `SoBRLDatabaseSource` captures and authenticates `BObolSourceRealizationStamp` before streaming and at final adoption: route, population, database binding, keys/path, source/input revisions, modes and tessellation tolerances. Camera/view policy remains nonbinding for source content. The caller owns producer completion and stream-drain certification. |
| Registry reuse | Only `authoritativeStreamDrained` certifies that the result has already been merged. The count check is a consistency check, not payload identity. Without that certificate, adoption installs the detached index. |
| Routing and dense indices | Index replacement preserves object-lifetime routing; creating a different source object issues a new route. `clearCompactInstanceIndex`, used by installation, advances the population epoch; appends preserve it. Public semantic handles and dense-index cache identities are distinct contracts. |
| Journals and request epoch | Installation invalidates the CAD/inventory journals and republishes the request contract for current source/input revisions. Source/input/database staleness revokes it. Non-consuming sparse journals retain at most 256 batches/65,536 entries each; readers below their floor rescan. Visibility has its own journal. |
| Semantic continuity | `compact_merge_runtime_state` transfers matching runtime presentation state before installation; retained selection/visibility frontiers are reapplied. The previous index is deleted before installation returns. The fallback uses existing population-sized staging; no second persistent registry is added. |
| Stream retirement | Certified reuse hides the temporary overview through the existing sparse CAD update while preserving the live assembly. The source retains a staging stream only while it owns staged bytes; provider claims or invalidation retire that lease. |

The follow-up journal review indexes 39 direct invalidation calls in
`/tmp/obol-ged-stream-denial-20260905/source-journal-callers.json`; its sibling
`source-journal-review.md` records the shared delta helper and metadata effects.
CAD presentation, geometry inventory and visibility have distinct writers and
rescan floors. Existing append, visibility, source-evidence and retained-geometry
tests exercise their ordinary paths. CXX-CONFIG-008 qualifies the twelve public
source-field notification paths; journal retention limits and the remaining caller
aliases still require qualification.

### CXX-SOURCE-005 — matching metadata falsely certifies an adopted payload

The unstreamed adoption path compared count, path, geometry kind, request
presence and transform, then treated the existing registry as the completed
result. A replacement can preserve all those values while changing immutable
geometry or provider inputs. The shortcut could mark the source current while
keeping the old payload and dense-index evidence.

The fix deletes that inferred certificate. The existing explicit drained-stream
path retains its fast adoption; other callers use the existing authoritative
index installation. No new lifecycle state, equality framework or registry copy
is introduced. The public API now states the caller's certificate obligation.

`test_compact_aabb_stream_upgrade` exercises same-path, same-tier and
same-transform replacement after source invalidation. It requires the new
immutable geometry, retained selection, unchanged routing, a new population
epoch, invalidated inventory deltas and the actual retained assembly payload.
Existing tests retain streamed overview retirement, camera independence and
staged-provider lifetime coverage. The [adoption checkpoint](obol_20260919_readiness_history.md#september-5-authoritative-source-adoption)
owns negative/fixed evidence and its limits.

The concrete metadata-to-payload comparison is below the formal source
observation abstraction; the executable counterexample supplies that evidence.
That adoption change needed no model relation or formal baseline change. The
following delivery audit checks GED's caller-owned certificate.

### CXX-SOURCE-006 — producer completion is not successful delivery

GED could stop draining after a memory-denied merge, then adopt a worker which
had already reached COMPLETE. This bypassed the denied merge by installing the
detached registry. Skipping adoption alone preserved the preview but left the
unmet expected population pending forever.

| Evidence | Owner and retirement rule |
|---|---|
| Worker outcome | The source coordinator owns production completion. It does not certify owner-thread delivery. |
| Delivery denial | `ged_obol_fail_stream_delivery` cancels its stream and publishes the source's existing FAILED outcome. Reservation and merge failures share this path; final adoption skips cancelled delivery. `CXX-SOURCE-007` removes the redundant item memory-exhausted flag. |
| Missing population | The source keeps its expected count and retained geometry. Controller inventory checks exclude failed sources from pending discovery while independently retaining live provider work. |
| Current error | Convergence derives failed sources from current rendered roots rather than the last synchronous realization action. Explicit source preparation and failure remain observable with automatic LoD off. |
| View changes and retry | `markStale(STALE_VIEW)` preserves a failed external producer's outcome and diagnostic. A view change does not create a replacement producer. Source replacement or an explicit realization can publish another outcome. |

The GED regression delays its small stream's drain until the real worker is
COMPLETE, then injects allocation denial at the production merge. It requires
unchanged preview geometry/population, retained expected inventory, a current
error without invented source work, successful rendering through policy
off/on/off, and a fresh draw reaching all four authoritative occurrences.
The coordinator regression also preserves unrelated pending preparation while
one source has failed. These extend existing tests and add no control state.

The lifecycle composition now includes provider failure, an error witness that
survives policy changes, and a failed terminal outcome after frame/work
retirement. Required success/failure action coverage and executable regression
mappings are cataloged. A mutant hiding the error with automatic LoD off fails
`ErrorHasWitness`. These mappings are regression evidence, not a claim of
complete stepwise refinement. The failure slice models one provider obligation;
it does not establish composition of a failed source with another live provider.
The coordinator regression checks that pending-work projection in C++.
Full-suite verification and the candidate's
limits are recorded in the [delivery checkpoint](obol_20260919_readiness_history.md#september-5-failed-source-delivery-and-display-transitions).

### CXX-SOURCE-007 — queued output targets a replacement source

GED's early drain checked numeric revisions and representation but omitted the
captured source route and path. Reusing a key and its revision numbers could
therefore merge an old overview into a different owner, or into the same object
after a path change. The later source adoption's path check was too late:
coverage bounds, profile, expected population and geometry were already mutable.
Simply skipping an invalid stream would also leave its queued output counted as
pending delivery.

The initial repair used `ged_obol_deferred_source_matches` as the common admission
predicate for streamed metadata/geometry and terminal adoption. It checked the captured route,
path, source/input revisions, source draw mode and representation after lookup
in the captured primary/view scene. Camera revision remains nonbinding; the
existing immutable launch policy owns any needed CSG successor. The unused
view-revision snapshot is deleted.

An invalid target cancels its existing stream without changing the replacement
source. Cancelled output no longer keeps the queue pending and cannot be adopted;
terminal job retirement releases its storage. Final adoption handles eligible
items individually, so one stale item cannot reject the entire completed batch.
The existing stream cancellation state also replaces the redundant
`memoryExhausted` delivery flag; actual allocation denial still sets the live
source's FAILED outcome and diagnostic.

The worker now observes item-stream cancellation at its existing stage and
completion boundaries. An item cancelled during production reports CANCELLED;
healthy siblings may finish and the batch may report COMPLETE. Job-wide failure
and explicit job cancellation keep their existing rules. Cancellation after an
item has completed production can still retire its delivery independently.

The production GED regression reuses equal numeric stamps in both a new-owner,
same-path case and an in-place, different-path case. Before the fix, occurrence
counts grow from four to five and one to two respectively. After the fix, the
replacement's bounds, geometry, placement, population and realized state remain
unchanged, and preparation retires without error. A controlled real-coordinator
test cancels one blocked item while holding a healthy sibling: the predecessor
reports the cancelled item COMPLETE; the fix reports CANCELLED, leaves the
sibling active, then completes it and releases the worker reservation.

The [routing checkpoint](obol_20260919_readiness_history.md#september-6-source-routing-and-scoped-cancellation)
owns negative/fixed binaries and graphical regression evidence. GED's stale
publication checks are regression mappings to `ObolResultAuthentication`, not
stepwise execution traces. Multi-item source-worker cancellation is checked in
C++; the one-provider lifecycle slice does not establish that composition.
`CXX-SOURCE-008` below moves admission into the source and adds same-lifetime
configuration/population identity. Cancelled storage retention while another
item runs, large-geometry cancellation latency and broader source-writer
qualification remain open.

### CXX-SOURCE-008 — same-owner reconfiguration admits an old stream

The routing repair still accepted a source with the same owner, path and numeric
revisions after an independent compact population replacement. The real GED
regression grows the replacement from four occurrences to five when the old
overview arrives. Database rebinding, representation-key changes and each
tessellation tolerance were also absent from the predicate.

`SoBRLDatabaseSource` now captures and authenticates the immutable
`BObolSourceRealizationStamp`. GED uses that predicate before metadata/geometry
drain, when recording an active CSG successor, and before terminal adoption.
Direct detached adoption requires the launch stamp and verifies the detached
configuration as well. This removes the adapter's parallel field list and
duplicate configuration comparisons. Stream appends preserve population
identity; installing another index invalidates it. Camera revision, view policy
and runtime appearance remain independent, preserving useful source production.

Database binding has a non-reused epoch, advanced solely by
`BObolCadSourceState::setDatabaseBinding`. Its four callers are the public setter,
batched source configuration, detached database initialization and detached
source creation. A new detached template already has a null binding. Rebinding
A to B to A rejects A's old output, while an idempotent assignment preserves
identity. This avoids retaining a borrowed database address after the worker
has released its snapshot-source lease; no extra database lifetime is required.
The existing identity helper fail-stops before counter reuse.

The GED test preserves numeric revisions and modes in nine cases: owner, path,
population, database, database round trip, representation key, and absolute,
relative and normal tolerances. Configuration-only cases retain the previous
population to isolate that boundary. All preserve geometry, placement, bounds,
population and realization status while stale work retires without error.
The update-action test separately rejects a superseded population through direct
unstreamed adoption, while its existing tests preserve same-tier geometry
replacement, semantic selection, valid appends and camera independence.

The [source-identity checkpoint](obol_20260919_readiness_history.md#september-7-source-realization-identity)
owns negative/fixed evidence. These tests refine the existing result
authentication contract; no model relation or accepted TLC baseline changed.
They do not complete the sparse-writer, journal-retention, partial-merge
allocation or broader cancellation/resource audits.

### CXX-PRESENTATION-001 — ordinary fallback leaves an obsolete frame witness

The same integration test exposed a separate presentation defect. After LoD
off/on/off, an ordinary source representation rendered while its cached LoD
assembly remained in the view's required execution set. Every subsequent image
was classified as incomplete because that unused assembly never executed.

`cad_view_lod_assembly_for_action` now retires the old view binding when it
selects ordinary fallback. The existing binding setter invalidates observations
before releasing the assembly; a subsequent frame proves completion against
the remaining bindings. No exact-frame predicate is weakened. The regression
requires real image rendering and terminal error after the full policy sequence.
Broader graphical and host qualification remains open.

### CXX-PRESENTATION-002 — exact zero invalidates a scene aggregate

The graphical delivery replay found that `refreshCadPresentationFrameStatus`
treated a fully culled assembly's exact zero work as a missing observation.
That invalidated the visible assembly's total and reopened an identical frame
repair indefinitely. Primitive and render-cost aggregation now preserve exact
zero contributions; frame execution and requested-control authentication are
unchanged. This changes derived evidence, with no new persistent state.

The real-renderer regression checks a visible/culled pair, an entirely culled
frame and an unexecuted frame. Its negative executable and the two-backend
GUI reproduction are retained in the [aggregation checkpoint](obol_20260919_readiness_history.md#september-5-exact-zero-work-assembly-aggregation).
The unchanged CAD-frame model checks lifecycle ownership; actual work magnitude
is qualified by C++. The following publication-boundary repair closes the
small replay's unnamed transitions; the overview finding below qualifies its
retained extent. Wider graphical
qualification remain open, so these passes do not establish complete refinement.

### CXX-TRACE-001 — source publication crosses an unnamed boundary

The small off/on/off and erase replay exposed three gaps in transition
accounting: `clearViewLodState` closed its cancellation scope before clearing
bindings; source notifications returned before tracing with automatic LoD off;
and GED source/policy changes could reach another controller operation before
their commit was notified. Actual debugger stacks identify each boundary.

The existing clear operation now has one enclosing external-input scope.
Source notification always records the application's committed changes but
publishes automatic planning work only when automatic control is enabled.
GED uses that same notification after shared/local source commits and policy
realization, before view-state invalidation. The post-frontier notification
continues to coalesce through the existing source-evidence owner. There is no
new pending state, event classification, or public transaction interface.

The production controller regression checks clear attribution and rejection of
a later unannounced write after a no-op notification. The existing GED delivery
regression now checks named, invariant-valid transitions across real policy
changes and erase while LoD is off. All four small graphical traces pass the
unchanged strict checker. The [publication checkpoint](obol_20260919_readiness_history.md#september-6-policy-and-erase-publication-boundaries)
records negative/fixed evidence and limits. This closes those reproduced gaps,
not the entire remaining writer audit or a claim of complete stepwise refinement.

### CXX-PRESENTATION-003 — overview geometry and placement lose coherence

The initial GED view-envelope builder discarded the transform accompanying a
shared unit box. The new regression measures actual transformed wire vertices,
then independently checks the retained assembly record. Keeping that transform
exposed two publication gaps: generic realization certified exact bounds without
queueing their overview, and population-sized reservation preceded the priority
overview's merge. Finally, a same-part overview replacement changed coordinate
systems without advancing placement, so the sparse presentation retained its
old normalization matrix.

The initial builder now keeps both returned values. Generic realization reuses
the existing internal coverage builder and priority queue before certifying
completion. Expected leaf count can be recorded without allocating instance
storage; GED commits a queued overview before reserving or merging leaves.
Geometry replacement advances the existing placement revision when its matrix
changes. The shared part and existing control owners remain intact.

The production GED regression captures the certified overview before denying
the next leaf merge, then requires the same geometry and transform, stable
population epoch, actual whole-source wire extent, honest error, real off/on/off
rendering, erase and recovery. A count-only test uses an unallocatable population
without creating an instance index. Six small GUI replays separately exercise
reservation denial, merge denial and success on both backends; their strict
traces and inspected images pass. The [overview checkpoint](obol_20260919_readiness_history.md#september-6-retained-overview-extent-and-delivery)
retains all three counterexamples and binary provenance.

These are C++ geometry, effect-order and regression checks below the lifecycle
model's abstraction. No formal relation or accepted baseline changed. They do
not establish physical OOM recovery, strong exception safety inside merge, or
whole-root coverage before a producer can establish its bounds.

### CXX-PRESENTATION-004 — checkpoint traversal is not presented evidence

The HUD visibility oracle compared current retained features with the previous
System GL framebuffer. The saved erased and redraw-return images have identical
pixels and presentation serial 67, but the latter reports newly published HUD
geometry and a pending render. The 100-pixel assertion therefore selected an
unpresented fill. The corrected selection requires an exact CAD frame and no
outstanding render request, and checks every sufficiently large non-idle
HUD. `CXX-PRESENTATION-006` adds the framebuffer-carried feature identity
needed when interrupted recovery has already retired that render latch. Missing presentation state, malformed geometry, insufficient pixels and
a hidden eligible fill still fail. Pending/incomplete frames do not authorize
an image assertion; the test does not fall back after an eligible image fails.

The same audit found an endpoint mismatch: software checkpoints used image
export and executed a renderer callback even though no presentation completed.
`QgSW::get_presented_frame_image` now returns painted pixels without traversal.
The canvas owns one reference to its displayed image, sharing the completed
fallback's storage during ordinary presentation and preserving distinct
provisional/clear output when necessary. Rebinding clears both image identities;
readback rejects a viewport size or pixel-ratio mismatch. qged uses a common
passive capture helper and requires a matching image at its startup barrier.

The Qt regression observes no renderer callback, unchanged presentation/request
serials, preserved pending work and unchanged captured pixels after mutation.
Software pixels are also compared with QPainter's actual painted orientation
and fractional display scale. Both backend branches pass on private Xvfb.
The [checkpoint evidence](obol_20260919_readiness_history.md#september-7-passive-checkpoints-and-hud-frame-correspondence)
owns the negative/fixed runs and remaining qualification. These are C++ and
pixel checks below the formal abstraction; no model relation, transition
checker or accepted TLC baseline changed.

### CXX-PRESENTATION-005 — terminal error HUD publication

The controller's ERROR phase can contain both a failed source and an unfinished
independent provider. The Qt publication grouping must distinguish that state
from the terminal error after foreground work retires. Otherwise an unchanged
phase bypasses the final sample once periodic publication stops.

`test_qtcad_obol_draw_sync` drives the production helper with its own publication
cache and periodic sampling disabled. A failed source, healthy sibling and
controlled provider establish active ERROR; removing the provider and completing
real controller renders establishes terminal ERROR. The old grouping fails to
publish that transition. The new grouping publishes a full red retained bar,
terminal cap and error label, owns a frame, and preserves feature revisions
without another render request on subsequent pre/post-frame observations.

The test retires preceding camera setup before initializing its cache; the
tested transition uses a bounded number of control/render steps. The GUI HUD
oracle now rejects a missing terminal fill with enabled LoD, and ordinary
denial replays verify actual pixels on both backends. The
[terminal-error checkpoint](obol_20260919_readiness_history.md#september-7-terminal-error-hud-publication)
owns the negative/fixed evidence. This changes a private observation grouping;
no controller lifecycle relation, strict trace predicate or TLC baseline changes.

### CXX-PRESENTATION-006 — CAD execution does not authenticate HUD pixels

The source-identity Generic System wire cold replay reported a yellow HUD but
captured no fill. Its image is byte-identical to the preceding completed frame:
presentation/completion serials stay at 42 while the interrupted count advances
from zero to one (50.01291 ms). The renderer correctly retained its last image.
The checker incorrectly treated exact CAD work and an empty render-request latch
as proof that the current feature store had reached those pixels. Recovery can
retire that latch while a bounded successor remains pending.

Qt now retains the traversed feature-store revision with the completed image.
`renderPending` returns that revision to the GL host on successful completion;
the software traversal captures it at the same boundary. Both hosts publish the
certificate with their actual painted image, preserve it on interrupted fallback,
and invalidate it for controller replacement or provisional/cleared output.
Viewport mismatch rejects an old certificate. These are bounded observation
fields; no scheduler, admission policy or retry predicate consumes them.

qged records current and presented feature revisions. The HUD oracle requires
matching identity before comparing active features to pixels. Unknown or older
provisional/interrupted pixels cannot authenticate current nonterminal features;
missing/malformed identity is an error. With no pending render, a terminal error
must have reached a matching framebuffer. The 100-pixel floor and hidden-fill
failures remain unchanged.

The real Qt regression publishes a newer HUD feature, interrupts rendering,
retires the request, and verifies unchanged pixels plus the old certificate.
Recovery must publish new pixels and the new revision. Passive capture preserves
pending work, and resize rejects old-viewport evidence. Both actual paint paths
pass at fractional scale on private Xvfb. The oracle's old implementation fails
the retained-interrupted-image regression; the new implementation passes all
37 synthetic cases while still rejecting hidden current and terminal fills.

The [frame-evidence checkpoint](obol_20260919_readiness_history.md#september-7-framebuffer-feature-evidence)
owns reports, images, binaries and qualification limits. This observation
boundary lies below the formal models; their accepted relations and TLC baseline
are unchanged. The saved Generic source-identity row remains historical failed
validation, with its false inference explained rather than its pixels altered.

## Detached source submission and release inventory

The audit covers `source_realization.cpp` and its public job/request contract.
Its private item, job, work, worker and coordinator structures are defined and mutated
in that translation unit.  The recorded 23 explicit state, queue, reservation
and active-job/stop mutations are indexed in `source-control-writers.json` under the
[worker-pool checkpoint](obol_20260919_readiness_history.md#september-5-source-worker-pool-construction-and-shutdown).
Constructor initialization, ordinary copied inputs and bypass aging are
classified separately below; the count is not a count of every assignment.

| State/resource family | Writer, identity and enabling evidence | Completion, release and bound |
|---|---|---|
| Request resources and callback lifetime | `submit` validates the whole request vector, prepares item/handle/queue storage, then consumes source/database/callback ownership under the queue lock | Invalid, closed or allocation-failed submission leaves requests unchanged. Admission denial consumes them into a CONSTRAINED job. `BObolSourceRealizationItemPrivate::finish` releases unsuccessful sources before database handles/snapshot files and drops callback context before publishing terminal item state. Successful source/database pairs remain borrowed results until job destruction (`CXX-SOURCE-010`). |
| Item input and production data | `submit` copies client token, draw mode, estimate, stream and callbacks. `source_realize_item` alone writes the detached database/snapshot, realization and manifest results | The owning job plus item index identifies work. Workers retain the job and active callback context. `itemResult` acquires item state before exposing a source, and exposes it only for COMPLETE. Clients retain the job handle for that borrow. GED's publication path requires aggregate COMPLETE. |
| Item outcome | `source_realize_item` publishes RUNNING; `BObolSourceRealizationItemPrivate::finish` alone publishes terminal item states after cleanup. Worker exceptions enter through `source_job_fail` | `CXX-SOURCE-007` makes the worker observe individual stream cancellation before/after its stages. It cancels that item while healthy siblings can finish; a COMPLETE batch may contain CANCELLED items. Job-wide failure still cancels siblings and remains the aggregate failure outcome. Exact geometry/provider completeness is separate qualification. |
| Job cancellation and terminal publication | `source_job_request_cancel` is the cancellation writer for clients, sibling failure and shutdown. Its first atomic request cancels streams; duplicate queue/worker references do not rescan the item vector. `submit` initializes remaining count; the owning worker or queue-retirement path decrements once after cleanup | The final decrement publishes COMPLETE, FAILED or CANCELLED with release/acquire ordering. Admission denial publishes CONSTRAINED without workers. The remaining count and aggregate state are distinct: terminal publication certifies reservation retirement as well as item writes. |
| Queue publication and rollback | `submit` appends a complete batch while holding the queue mutex; failed staging removes only its appended suffix. The worker erases the selected item; cancellation splices cancelled entries into caller-owned cleanup without allocation | Existing entries survive a failed submission. No partial batch is executable. A prepared public handle has no interest until commit, so abort cannot cancel caller-owned streams. No queue copy or new persistent phase is needed. |
| Active count, job and working-set reservation | The worker admits under the same mutex, records its job in its fixed worker slot and increments count/bytes. It executes outside the lock, then clears the slot and reservation before final job publication | Shutdown can reach an active job after it leaves the queue. One reference per worker bounds this storage. Directory-derived admission and eight-bypass aging remain unchanged. Actual memory magnitude and geometry callback latency retain their release gates. |
| Worker-pool lifetime | Construction prepares fixed slots in an RAII candidate after cache initialization. `source_realization_stop` cancels queued and active jobs, wakes and joins workers, then waits for caller-owned queue retirement; both failed construction and normal destruction use it | `CXX-SOURCE-003/004` cover partial-start cleanup/retry, GED error containment and cooperative shutdown with a saturated reservation. `CXX-SOURCE-011` adds a counter for cleanup already detached from the queue, with weak endpoints and a shutdown barrier. Native platforms, real geometry under pressure and full shared-stack sanitizers require their own qualification. |

This inventories the detached-production boundary, not the complete source
registry/adoption, LoD residency/service or GED effect surfaces.  Those stay in
the existing S1/S3 backlog.  The queue and result lifetime are not revision
counters; renderer/view freshness remains the downstream adoption contract.

The current coordinator writer delta is retained as `source-control-writers.json`
in the queued-cancellation artifact. It explicitly includes endpoint installation,
list transfer, the two wakeup edges and the shutdown retirement counter. The
original mutation count above is dated evidence, not the size of this delta.

### Direct realization writer inventory

This audit classifies the remaining terminal-record writers in
`database_source.cpp` by their actual callers. A helper's name does not establish
whether its source is private. Live cached realization now prepares a private
candidate; workers explicitly use the exclusive construction entry.

| Writer family | Owner, identity and trigger | Publication or retirement boundary; remaining acceptance |
|---|---|---|
| Detached template and worker construction | `createDetachedRealizationTemplate` owns its field copy and omits invalidation sensors; `source_realize_item` explicitly calls `bobol_database_source_construct_realization` with its exclusively owned source/database | Intermediate private fields do not require live-scene publication. Stream entries and terminal adoption keep their existing boundaries. Wider constructor/traversal cleanup and actual geometry/snapshot limits remain open. |
| Wire/mesh success and failure | The four cached live wrappers use `BObolPreparedSourcePublication`; `mark_source_realized_current` serves private construction | `CXX-SOURCE-019` publishes candidate geometry, bounds, evaluated roles, owned metadata and child/compiled retirement together. Provider failure retains geometry and publishes FAILED through the status setter; the action's duplicate database-failure write is removed. Public NULL/stream semantics and private worker handoff are retained. |
| External primary geometry | Line, point, triangle and annotation publication; `clearExternalPrimaryGeometry`; primitive line-set, submodel and plot adapters | Source configuration identifies the replacement. `CXX-SOURCE-018` prepares geometry, child/path edits, bounds, terminal metadata and compiled/compact retirement together. Clear retains existing primary shape objects; ordinary replacement preserves auxiliary owners and placement. The submodel adapter replaces its owned leaf population. |
| Compact edit refresh | `refreshCompactObjectGeometry`, called by the scene controller with an edited path and next source revision | `CXX-SOURCE-020` prepares combination edits privately. `CXX-SOURCE-021` prepares only affected leaf occurrences and owner metadata before committing revision, geometry, part ownership and exact bounds. Indexed paths resolve, obsolete profiles are invalidated and source-mesh currency advances. Both retain live sensors. |
| Auxiliary line publication and retirement | `setAuxiliaryLineSet`, empty auxiliary source input, `clearAuxiliaryShapes` | `CXX-SOURCE-022` privately builds named line replacements and publishes placement/children/paths together. Semantic overrides survive; primary geometry and compact handles remain retained. Empty input and clearing use prepared removal. Nested auxiliaries inherit one parent placement. |
| Auxiliary source publication | Nonempty `setAuxiliarySourceLineSet` creates or updates a named child source, copies display/placement state and writes its realized record | `CXX-SOURCE-023` prepares configuration, owned metadata, geometry and child/path edits before committing the complete source. Named routing identity and independent nested owners survive; obsolete bounds/profile/staging and earlier worker stamps are invalidated. Ordinary controller revisions advance at commit before observers, including throwing callbacks. |
| Prototype realization | `realizePrototypeWireframe`, reached through the realization action for diagnostic sources | `CXX-SOURCE-024` builds a privately owned diagnostic shape and publishes it through the existing prepared source boundary. Its geometry, source-local bounds, owned metadata, paths and obsolete compact/compiled retirement agree before callbacks. Auxiliary and nested owners survive. |
| Scene realization effects | `BObolSceneController::realizePending` and `SoBRLRealizeAction` consume the existing source commit | `CXX-REALIZE-001` publishes action progress and each changed source's frame revision before observers. Prepared diagnostics and scoped progress retention preserve the completed prefix on exceptions; rejected preparation and unchanged failure retries do not advance the frame. |
| View realization successor | `BObolViewController::realizePending` contributes prepared effects to the existing source commit | `CXX-REALIZE-002` publishes the standing request with source/action/scene state before observers. Host notification is attempted despite source-observer failure; allocation rejection preserves the preceding request, completed prefixes survive, and host-consumed requests are not recreated. |
| Realization invalidation | `markStale` revokes realization currency and, for source/input/database changes, the source-mesh contract and staged producer lease | `CXX-SOURCE-025` prepares the realization record, exact-bound retirement and owned legacy records before committing. Immutable geometry, compact records and independent nested owners survive; view/draw/tessellation changes retain usable source data. Configuration callers own their separate metadata/display propagation. |
| Configuration and retargeting | Batched source configuration, representation/draw setters, source retargeting and direct field/engine notification | CONFIG-001 closes batched configuration and the separate draw/representation/retarget setters; CONFIG-002 closes their scene effects. CONFIG-003 closes display-name/hierarchy/material-policy setters; CONFIG-004 closes display/placement; CONFIG-005 closes bounds, explicit realization, roles and view-policy scene setters; CONFIG-006 closes explicit database metadata and source/scene material refresh. CONFIG-008 closes direct assignment, connected-field/engine evaluation, delayed delivery and quiet/exception edges for all twelve watched fields. Each qualified boundary preserves its own identity and invalidation rules. Enclosing multi-step callers and traversal/custom cleanup retain their inventory and acceptance rows. |
| Sparse entry-local display and metadata | Seven setters in `database_source_compact_access.cpp`; indexed exact/subtree/object or legacy path matching | `CXX-SPARSE-001` prepares changed records and membership capacity before committing semantics, appearance, memberships, renderer records and compiler/journal effects. Geometry, source stamps and retained rules keep their owners. Retained intent is covered separately by SPARSE-002; enclosing GED operations remain open. |
| Retained presentation and selection | Visibility/highlight/transparency rule setters, highlight clearing, visibility frontier/override replacement and clearing, selected-path replacement/deltas | `CXX-SPARSE-002` prepares retained intent with affected records, commits memberships, renderer/compiler and journal effects before observers, and shares matching with later arrivals. Three separate live reapply writers are removed. Enclosing GED operations remain open. |
| Evaluated-region publication | GED resolves source targets; `setEvaluatedRegionForPath` owns matching shape fields; the scene controller contributes frame effects | `CXX-REGION-001` prepares all matching owned wire/mesh markers or existing compact records, retires compiler/batch evidence and commits each source's scene revision before observers. Completed prefixes and later reentrant edits/replacements survive. Two unused private GED helpers are removed. Other GED metadata/application operations remain open. |
| GED metadata application and material refresh | `applyDatabaseSourceInstanceDrawMetadata` normalizes the existing typed cache record and selects aggregate or compact scope; material sweep retains original source targets | `CXX-GED-001` contributes scene effects to both prepared metadata publishers, checks refresh target owner/material/CAD/node currency, and commits each changed source's frame effect before observers. Duplicate GED conversion and key-by-key refresh are removed. Enclosing direct/deferred draw transactions retain acceptance work. |

The scene-level `publishDatabaseSourceInstance` remains a distinct unqualified composition
of those setters; CONFIG-007 below records its reproduced intermediate-state failure.

This is the terminal-writer classification, not completion of the wider S1
source/service, sparse setter, GED or presentation inventory. Temporary storage
for a replacement may own its new geometry and changed child/path records; it
must not create another persistent occurrence model or copy the entire old scene.

### CXX-SOURCE-017 — failed append leaves unowned references or parent links

`SoBaseList::append` and `insert` acquired a node reference before growing list
storage. Allocation denial left that reference outside the list. `SoChildList`
installed its parent auditor before the same growth, so the failed append could
also leave notifications flowing from a node absent from the child list.
The focused negative tests are retained in `/tmp/obol-child-append-20260907`.

The repair reserves child-list storage before connecting its parent auditor,
then publishes the slot and reference without further allocation. Base-list
append and insertion acquire their references only after storage succeeds.
`SbPList::reserve` preserves length and contents and uses the existing growth
policy. Notification failures occur after a complete append; ordinary append
notification order and existing path indices are unchanged.

The production tests pass five append allocation positions (three preserved,
two committed), covering new and already-shared children, exact references,
parent notification routing after denial, an existing path, retry and throwing
immediate observers. The separate base-list test denies both append and insert
growth. The
[dated evidence](obol_20260919_readiness_history.md#september-7-child-append-ownership)
owns qualification. This does not qualify individual child insertion/removal,
replacement/truncation ordering or complete source replacement.

### CXX-SOURCE-018 — external replacement discards the preceding drawing

The retained line-set probe originally lost the preceding drawing on allocation
denial and exposed incomplete state in 53 of 91 ordinary callbacks. With the
repaired source, denial retains the old geometry and sends no callback; ordinary
replacement produces nine callbacks, all coherent. The final production test
also rejects the preceding libBObol while using the current Obol dependency.
[The dated evidence](obol_20260919_readiness_history.md#september-7-external-geometry-publication)
owns `/tmp/obol-external-publication-20260907` and the qualification artifacts.

| Boundary | Preparation and commit | Qualified behavior |
|---|---|---|
| Nonempty line, point, triangle and annotation | Construct and populate an unpublished shape; prepare source/owned auxiliary metadata, placement and the complete child edit; commit before notifying | Failure preserves the preceding drawing. Bounds and requested/realized stamps agree with the published geometry; nested sources retain their own metadata |
| Empty geometry and clear | Prepare scalar metadata and hold shared geometry while preserving existing primary shape objects; clear arrays and obsolete compact/compiled state in the same commit | Existing primary pointers survive clear. Shared geometry disconnects safely; old handles and compiled children retire with publication |
| Obol child replacement | Reserve lists and parent links; map repeated child occurrences in order; retain affected paths and retired nodes | Commit does not allocate or notify. Retained paths move to their new indices; removed paths truncate before callbacks. Abandoned preparation restores references and parent routing |
| Primitive adapters | Own line/annotation payloads and temporary plot lists; stage submodel leaves and transformed aggregate bounds before publication | Double coordinates survive line conversion; failed plot conversion returns its blocks. Repeated submodel instances remain distinct with database instance IDs enabled or disabled |
| Notification and GED integration | Restore participant flags, attempt remaining callbacks and propagate the first exception; use the record visitor's actual continue/stop convention | Throwing observers leave a complete usable drawing. GED resolves a cleared noncompact source after unrelated records and can update its highlight |

The final external sweep covers **25 scenarios and 2,591 allocation positions**:
2,542 preserve the old drawing and 49 leave complete commits. It includes
source/shape/path observers, explicit/default commands, precision, annotation
text, LoD meshes, clear/empty inputs, shared geometry, compact handles, actual
compiled children, missing placement, disabled fields, retry and throwing
observers. Invalid inputs are rejected without allocating or notifying. The
prepared-child replacement sweep passes **21 positions**, 20 preserved and one
committed, including duplicate/reordered children, path mapping, abandoned
preparation, null rejection and idempotent commit/notification.

The repair also fixes float-only annotation points being erased by an absent
optional precision array, and submodel occurrences sharing an identity when the
borrowed database disables combination instance IDs. It reuses the existing
transaction helpers and removes the superseded external publication paths.
Four primitive-edit replays exercise command/widget/retained-handle agreement
in single/quad layouts on both backends; ten small drawing replays and eight
Generic rows/camera contracts pass.

This qualifies the external publication boundary on the tested owning scene
thread. Connected fields/engines, custom deletion callbacks, remaining live
writers, full allocator/sanitizer coverage, complete traversal teardown and
native hosts remain open. Editing captures establish neither pixel equivalence
nor complete label-layout acceptance. S1--S6 remain incomplete.

### CXX-SOURCE-016 — realization status publishes before owned shape metadata

`setRealizationState` notified each source field as it changed, then rebuilt
compact display state and copied metadata into child shapes. The new production
regression fails without allocation denial: an immediate observer sees a mixed
source/shape record. The expanded case also shows a parent update overwriting
metadata owned by a nested database source. Both negative results are retained
under `/tmp/obol-realization-state-20260907`.

| Owner | Preparation and commit | Failure and bound |
|---|---|---|
| Realization-status setter | Normalize the requested status/revisions; prepare the eight realization fields for the source and each owned legacy wire/mesh shape; copy changed diagnostics before mutation | Failure preserves the old records. An identical request allocates nothing and sends no notification. Failed/unrealized states preserve prior realized stamps and identity; successful explicit/default stamp semantics are unchanged |
| Shape ownership | Visit each legacy node once, retain prepared shapes, and stop at nested database-source owners | Shared children are prepared once; nested source state remains independent. Storage scales with the visited legacy graph and changed diagnostics, with no compact-occurrence mirror or scene-geometry copy |
| Commit and invalidation | Commit every prepared record quietly, invalidate batch eligibility on a transition, and invalidate the compiled assembly before observers | Status-only changes leave compact geometry, display records and memberships alone; placement/material/display setters retain their own duties |
| Field publication | Restore all participants' notification flags, then attempt changed-field notifications and preserve the first exception | The same fixed-size notification guard now serves terminal adoption. Callback failure leaves complete committed values and does not strand later notification attempts |
| Obol field-auditor list | Hold its recursive notification lock through scoped ownership during mutation and snapshot preparation | Allocation-denied snapshot construction previously leaked this lock despite a zero SoDB notification counter; the subsequent owner thread now completes |
| Obol pointer-list storage | Allocate replacement storage before publishing capacity | Failed growth preserves the old length, capacity and elements, so a retry cannot write beyond the old buffer |

The realization-state sweep passes 13 allocation positions, with 11 preserved
and two committed outcomes. It checks immediate source/shape observers, shared
children, nested owners, compact-record preservation, restored flags and retry.
Additional cases exercise explicit revisions, invalid/unrealized requests,
disabled field notifications, throwing observers and child removal by a callback.
The unchanged-request check runs under allocation denial and observes zero
allocations and callbacks.

The full regression originally deadlocked in its later cross-thread notification
check. A retained symbolized GDB stack locates the blocked thread in the recursive
notification lock. The isolated field-auditor case reproduces the hang against
the preceding dependency, while the pointer-list probe rejects the prematurely
published capacity. Both pass with the repair. The
[dated evidence](obol_20260919_readiness_history.md#september-7-realization-state-publication-and-auditor-allocation)
owns final build, test and graphical results.

This closes the status setter used by scene-controller realization updates and
GED failed stream delivery. It does not qualify the complete GED diagnostic
path under persistent physical OOM. `mark_source_realized_current`, direct
realization success/failure writers, traversal/legacy child mutation, other
sparse setters, field connection/engine and delayed-sensor lifetimes, and
throwing deletion/custom cleanup remain open. The other gates retain their
resource, Lucy/scale, editing and platform requirements. No formal relation or
accepted baseline changed; S1--S6 remain incomplete.

### CXX-SOURCE-015 — path callbacks observe partially retired source children

The `014` source/field observer contract did not include direct path sensors.
Both original probes fail against that checkpoint: a path callback sees completed
metadata while an obsolete child is still attached. All four strengthened terminal
scenarios fail with the preceding libBObol. The existing `SoChildList::remove`
notification order is intentional legacy behavior, so terminal adoption now uses
an explicit prepared removal instead of changing that API's callback order.
`/tmp/obol-child-retirement-20260907` retains the reproductions and repairs.

| Owner | Preparation and commit | Failure and lifetime |
|---|---|---|
| `SoChildList::Removal` | Validate ordered removal indices; prepare path truncation/index changes and parent-auditor retirement; retain the parent, removed nodes and affected paths | Preparation and abandonment leave the live graph and reference counts unchanged. Metadata storage scales with removed children and registered path auditors; node graphs/path chains are shared |
| Removal commit | Update paths quietly, preserve parent notification links for surviving occurrences of shared children, then compact the child list once in order | The retained references keep graph nodes alive through path and slot updates; no allocation or change callback is needed during commit |
| Source adoption | Prepare non-auxiliary child removal before the terminal record; commit metadata and child/path state before notification | Placement/auxiliary children survive. Field notification failure does not prevent the path-notification attempt; the first exception is rethrown after publication attempts |
| Immediate notification queue | Finish the existing bounded queue pass after an observer throws, then rethrow the first exception | Remaining immediate callbacks are attempted; queue locks and processing flags retire through scope guards. The existing rescheduling limit is unchanged |
| Sensor/path/child-list retirement | Detach cancels queued node-sensor work; path and child-list destructors remove graph/auditor links without starting change notifications | Destruction cannot dispatch unrelated immediate callbacks through these paths; shared-child notification links are removed from the dying parent |
| Last-reference auditor retirement | Remove the next live sensor from the auditor tree before its dying-reference callback, then find the next live entry | No allocating sensor snapshot. A deletion callback may retire another sensor without leaving a dangling entry to visit |

The source regression now includes direct path observers and an auxiliary child.
Its four terminal sweeps cover 217 allocation positions (205 preserved, 12
committed). The separate prepared-child sweep covers nine positions, including
abandoned/invalid preparation, interleaved removals, path indices, shared children,
retired-node release with eight direct sensors, throwing observers and recovery.
Node-sensor detach and quiet path/group destruction each fail with their preceding
library. Allocation during node release previously aborts; a deletion callback
removing another auditor previously segfaults. Both now pass. The
[dated evidence](obol_20260919_readiness_history.md#september-7-child-path-and-auditor-retirement)
owns final runtime and graphical qualification.

This qualifies the prepared operation and these default node/path/sensor release
paths. Callers must keep the child list unchanged between preparation and commit.
Other source writers and callers of legacy individual removal/truncation still
need their own transaction audit. Throwing deletion callbacks, custom-node
cleanup, connected-field/engine paths outside the twelve CXX-CONFIG-008 source
fields, traversal cleanup, actual resource limits and native/full-sanitizer
qualification remain open.
No formal relation or accepted baseline changed; S1--S6 remain incomplete.

### CXX-SOURCE-014 — terminal adoption exposes mixed realization state

Registry consistency alone did not make completed adoption coherent. The
corrected persistent allocation probe finds six mixed replacements and nine
mixed streamed adoptions against `CXX-SOURCE-013`. A temporarily hidden overview
also kept authored visibility, allowing a later presentation change to revive it.
The four strengthened terminal scenarios fail against that preceding library.
`/tmp/obol-terminal-adoption-20260907` retains the negative and positive evidence.

| Owner | Preparation and commit | Failure and bound |
|---|---|---|
| `PreparedSourceRealization` | Prepare the identity string, desired bounds and change flags; quiet source/terminal fields while committing revisions, status, identity, stale state, bounds and root selection | One fixed fourteen-field record and one identity string; preparation failure restores notification settings without changing the live source |
| Streamed adoption | Reserve overview-retirement indices and hidden membership capacity before changing entries; retire both visible and temporarily hidden authored overviews; publish sparse changes and current request stamps | Existing registry/population identity survives; temporary storage scales with overview count, and previously hidden membership is retained without duplicate insertion or a global sort |
| Replacement adoption | Use the owned-index installer from `013` with field notifications quiet, then transfer staging and commit the prepared terminal record | Failed preparation preserves the source; consumed donor ownership is released. Retry constructs a fresh donor |
| Terminal publication | Commit staging and metadata before ordinary child/path retirement, restore every notification flag, then notify changed fields | Immediate source/field observers see the complete record and current LoD contract. A later notification failure leaves a complete committed source and usable notification machinery |
| Obol notification machinery | `SoSFString` accepts a prepared moved value; scoped notification ownership balances the global counter and lock; scoped field flags clear on unwind; immediate-queue insertion releases its lock if allocation fails | No new field layout or scene copy. Failed propagation defers remaining queued callbacks until a later normal transaction. Callback exceptions propagate; cleanup does not allocate diagnostics |
| GED terminal delivery | Catch adoption allocation failure through the existing failed-delivery path, cancel further delivery and publish availability change | The controlled terminal-preparation fault retains all delivered leaves and exact overview with FAILED status; policy changes and erase remain functional |

The final four terminal sweeps cover 185 allocation positions: 173 preserve the
old record and 12 leave a complete new record. They compare metadata and bounds,
registry/renderer state, revisions, notification settings, child-path retirement,
staging ownership/claim/release and retry. Baseline allocation counting includes
the same immediate observers and child path as injected runs. Failure persists
while unwinding. A separate throwing-observer test checks field/node recovery
and the next writer on another thread. Both notification regressions fail with
the pre-repair Obol library. The seven other installation sweeps and ten merge
sweeps retain their separate registry and prefix contracts. The
[dated evidence](obol_20260919_readiness_history.md#september-7-terminal-adoption-and-notification-recovery)
owns totals and graphical qualification.

The later [child/path retirement audit](#cxx-source-015--path-callbacks-observe-partially-retired-source-children)
addresses direct path observers and default graph/auditor release. At the `014`
checkpoint, unqualified boundaries included direct `SoPath` auditor exceptions
during child-list mutation, connected-field/engine and delayed-sensor paths,
other realization writers such as `setRealizationState` and
`mark_source_realized_current`, or traversal cleanup. GED's diagnostic path is
not a persistent physical-OOM survival guarantee. Actual completed-result bounds,
real snapshot cleanup, native hosts and shared-stack sanitizers remain open.
S1--S6 and the accepted formal baseline are unchanged; no TLA+ relation changed.

### CXX-SOURCE-013 — failed full-index preparation leaks geometry or invalidates live handles

The full-registry setters built a raw replacement index, then applied allocating
visibility and selection updates after installing it. A persistent three-leaf
allocation sweep finds 140 inconsistent positions out of 168: 103 retain
unpublished geometry, 37 publish obsolete bounds, and 30 have inconsistent
presentation (categories overlap). Detached adoption also changed the live
handle-source identity before candidate preparation; allocation zero invalidates
an otherwise unchanged live handle. All nine new installation scenarios fail
against the preceding `CXX-SOURCE-012` library.

| Owner | Preparation and commit | Failure and bound |
|---|---|---|
| Full-registry setters and `realize_walk_data` | Own the candidate with `unique_ptr` from construction through installation | Construction/preparation failure releases the candidate and its shared geometry references; the existing complete candidate remains proportional to the replacement population |
| `compact_prepare_installation` | Transfer runtime occurrence identity/state and rebuild candidate lookup/derived data | Only candidate state changes; the old population, handles, bounds and journals remain valid on preparation failure |
| `compact_prepare_registry_presentation` | Apply current presentation overrides and indexed visibility/selection masks; prepare renderer styles and memberships | Reuses the existing indexed frontier/selection traversal. Two temporary scalar masks replace live reapplication masks; no new scene-sized index or persistent occurrence object |
| `installCompactInstanceIndex` | Transfer the complete index and handle-source identity, advance population/CAD/LoD publication, retire history and publish coherent bounds before field notifications | Preparation has one owner until commit; existing exact whole-source bounds survive replacement. Notification failure cannot leave a half-installed registry |
| Detached replacement | Take the donor index into the same owned candidate and change the live handle identity only when that index commits | Failed preparation preserves the live source; the consumed donor candidate is released. Retry uses a fresh detached result |

`test_compact_publication.cpp` now includes nine complete-installation sweeps:
registry, singleton, detached replacement, active rules, empty source, new paths,
exact source bounds, singleton rules and detached rules. It requires the whole
old or whole new registry, stable old handles and publication identity on failed
preparation, coherent public lookups/renderer records/opacity/memberships,
source aggregates and bounds, release of uncommitted geometry, and successful
retry. It also checks that complete installation reapplies current selection,
including deselection of authored state. Failure remains active while unwinding.
The existing ten merge sweeps remain separate, preserving their prefix oracle.
The [dated evidence](obol_20260919_readiness_history.md#september-7-owned-source-index-installation)
owns allocator totals, focused checks and graphical qualification.

This checkpoint qualifies candidate ownership and registry consistency. The
later [terminal-adoption audit](#cxx-source-014--terminal-adoption-exposes-mixed-realization-state)
strengthens its detached oracle. At the `013` checkpoint, the remaining limits
were: Stream-preserving overview retirement,
child retirement, final realization fields/staging transfer, and traversal
resource cleanup still require their own exception/lifetime audit. In particular,
the detached sweeps admit a complete new registry when a later realization-field
write throws; they do not certify those fields as an atomic terminal result.
Other sparse setters, actual memory magnitude, native hosts and shared-stack
sanitizers remain in S1--S6. No TLA+ relation or accepted baseline changes: this
allocation/ownership repair is below the current formal abstraction.

### CXX-SOURCE-012 — allocation failure exposes an incomplete compact merge

The previous merge appended live records before building all lookup tables and
publishing its journals. Its part helpers could also retain unpublished geometry
or a key without a part. A three-leaf allocation sweep finds 111 inconsistent
positions out of 143: new entries cannot resolve through their handles, or the
view revision still describes the older population. Same-tier source enrichment
can expose a half-written request. An evolving overview can replace the part
before its occurrence changes. The earlier GED fault precedes leaf mutation and
does not cover these boundaries.

| Owner | Preparation and commit | Failure and bound |
|---|---|---|
| `CompactPartPublication` | Prepare content/private-overview identity, map slots and part capacity; commit the immutable binding after the occurrence is prepared | Stack guards erase new uncommitted slots; weak pointer witnesses do not retain abandoned geometry |
| `compact_append_prepared` | Prepare one entry, source presentation rules, instance/member capacity and exact/key/ordered/leaf/name lookup slots; deque insertion is the final allocating operation | Parallel records, lookup values, reference/channel counts and bounds commit together; no scene-sized copy |
| Geometry replacement / source enrichment | Prepare one local entry and required slots; reconstruct at its stable deque address with a checked nonthrowing move constructor | The previous occurrence remains complete if preparation fails; one prepared entry and shared immutable arrays bound staging |
| `mergeCompactOccurrences` | Reserve the input batch's change list and publish every completed occurrence as a prefix on return or exception | Retain population identity, source-contract stamps, overview coverage and current bounds; GED may terminate delivery without losing that prefix |
| `database_source_record_compact_delta` | Prepare sparse indices and the journal node before updating retained counts | Allocation failure clears sparse history and advances its floor, requiring an authoritative refresh; existing 256-batch / 65,536-entry limits remain |
| `setSourceBoundsState` | Assign coherent numeric valid/exact/min/max fields before notifying observers | A notification failure cannot expose mixed bounds; storage is four field pointers and flags |

`SbString` now has defaulted nonthrowing move operations. The registry asserts
nonthrowing entry and instance construction at compile time. This changes no
persistent record layout. Append and overview-to-leaf publication share the
existing presentation-rule application; a leaf inherits active visibility,
selection, highlight and transparency rules before it becomes observable.

`test_compact_publication.cpp` injects persistent scalar/array allocation failure
at every allocation position in ten production merge scenarios: append, empty
source, short paths/shared leaf keys, active presentation rules, geometry
replacement, replacement with rules, source enrichment, overview growth, and
root overview-to-leaf replacement with/without rules. It checks old-or-complete
entry state, public lookup routes, source aggregates and bounds, sparse refresh,
retained renderer records/geometry/memberships, release of unpublished geometry,
and successful retry. All ten scenarios fail against the predecessor library.
The final sweep and graphical results belong to the
[dated evidence](obol_20260919_readiness_history.md#september-7-compact-merge-publication).

GED's additional deterministic fault waits for the small worker to complete,
then fails inside a multi-occurrence merge after its first committed leaf. The
integration requires that leaf and the exact whole-root overview to remain
usable while source delivery becomes FAILED; cancelled delivery cannot install
the completed worker's whole registry. This is separate from the exhaustive
allocator sweep, which keeps failure enabled during stack unwinding.

The allocation/ownership mechanism is below the current formal abstraction; no
TLA+ relation or accepted baseline changes. This closes the tested merge
consistency boundary, not S1/S2 or physical memory qualification. Full source
installation/setter exception safety, other sparse writers, completed-result
magnitude, real snapshot cleanup, cache contention, actual memory pressure,
native clients and shared-stack sanitizers remain in active debt. The bounded
journal's full-refresh fallback does not certify every writer that feeds it.

### CXX-SOURCE-011 — cancelled queue entries wait for unrelated admission

The original queue selected work only when its source-memory reservation fit
and a worker was available. Neither job cancellation nor stream cancellation
removed queued entries. Ten initial regression cases fail behind a full memory
budget or a fully occupied pool: individual stream cancellation, whole-job
cancellation, interest-handle destruction, cancellation before submit and jobs
without streams all leave three queued entries and their resources retained.

The coordinator now owns direct queued retirement through
`bobol_source_realization_cancel_queued`. It checks existing item/job cancellation
facts, splices matching list nodes into a local owner under the queue mutex,
then finishes items and retires their existing remaining counts outside that
mutex. A list replaces the deque so detachment needs no allocation or scene-sized
copy. Worker admission and bounded bypass policy are unchanged. Individual
cancellation preserves healthy siblings; completed results retain their existing
borrowed lifetime.

Streams receive a weak coordinator endpoint only when submission commits. The
binding reports cancellation which preceded registration so submit retires it
after unlocking. A later cancellation closes payload and notifies the endpoint
outside the stream mutex. Job cancellation also notifies after closing streams,
covering streamless requests and shared-stream peers. Item cleanup closes
publication without recursively rescanning the queue for every retired sibling;
admission denial sends one queue notification after its batch closes publication.
The public coordinator's private pointer now shares ownership of its internal
service state; weak endpoints do not extend pool lifetime or recreate the global
service after teardown.

The new `queueRetirements` counter is a resource-lifetime obligation, not an
admission policy: the cancellation helper increments it under the queue mutex
before taking cleanup responsibility, decrements after releasing that ownership,
and shutdown waits for zero after joining workers. A controlled caller destructor
remains active across actual coordinator shutdown. Removing just this wait in an
isolated library makes the shutdown observer fail. A retained completed stream
can still be cancelled after coordinator destruction.

Removing a protected queue root can enable healthy peers immediately. The
initial implementation notified only after cleanup; a destructor waiting for
those peers reproduces a 2.47-second failure. Notification now occurs as soon as
the queue changes, and again when retirement completes for the shutdown barrier.
The final test covers this enabling event and the existing fairness contract.

Sixteen queued cases cover both resource conditions, adding concurrent stream
cancellation and shared-stream peers under job cancellation and admission denial.
They assert source/context release, no cancelled callback execution, unchanged
active reservations, and healthy sibling completion. Both production translation
units now compile independently in normal CMake builds. This exposed a missing
`bu/str.h` include in the stream file that its former unity neighbor supplied.
The failed compile and final normal compile commands are retained. The
[readiness checkpoint](obol_20260919_readiness_history.md#september-7-queued-cancellation-and-admission-wakeup)
owns negative/positive binaries, traces and graphical evidence. No formal
relation or baseline changed. Cancellation work still scales with queued roots
and owned data; large-geometry latency, completed-result bounds, snapshot-file
cleanup and native/shared-stack sanitizer qualification remain open.

### CXX-SOURCE-010 — terminal items retain worker resources behind siblings

The per-item stream fix did not retire detached sources, database handles or
callback contexts. Six controlled warm/cold cases reproduce the next boundary:
a first item completes, cancels or throws while a sibling callback remains
active. Every old terminal item retains callback context; cancelled/failed
sources survive, and cold cancellation/failure plus warm-probe failure retain
their database files. The old result API also exposes a worker-mutated source
before completion, and GED keeps its own raw alias after transferring ownership.

`SourceRealizationItem::finish` now owns terminal item publication. After the
last callback returns or unwinds, it cancels unsuccessful streams, unreferences
the unsuccessful source before closing its database/snapshot resources, drops
the callback context, then release-publishes the outcome. Worker reservations
still retire before aggregate terminal publication. Admission denial and queued
shutdown use the same item cleanup outside the coordinator mutex. No additional
state, resource mirror, scheduler or thread is introduced.

`itemResult` first acquires item state and reads the source pointer only for
COMPLETE. That source/database pair remains valid while the job handle is held,
even if a later explicit cancellation or sibling failure cancels the batch.
This necessary result ownership is distinct from completed callback scratch.
GED clears its transferred pointer, deletes `ownsSource`, and obtains the source
from the completed result at adoption. It still checks stream cancellation,
aggregate completion and the live source-owned launch stamp.

The production regression checks callback/source destruction before sibling
release, database file handles on Linux, preservation of borrowed completed
sources and cold database readback after cancellation, and final handle cleanup.
The shutdown observer checks active and queued source/context destruction while
client job handles remain retained. Admission denial checks the same release
boundary without executing a callback. User cleanup queries coordinator
accounting, exercising release outside its mutex.

The [readiness checkpoint](obol_20260919_readiness_history.md#september-7-terminal-item-resource-retirement)
owns negative/fixed reports and graphical qualification. Public result layout
and signatures are unchanged; non-COMPLETE source exposure is deliberately
removed. Formal cancellation outcomes and the accepted TLC baseline are
unchanged. Completed result retention, queue-wait cancellation, real snapshot
file cleanup, large-geometry magnitude/latency, native hosts and sanitizers
retain separate qualification requirements.

### CXX-SOURCE-009 — cancellation retains stream ownership and accepts late output

The preceding stream cancellation method only set an atomic flag. A controlled
two-source job reaches CANCELLED for one item while its healthy sibling stays
active, yet retains one queued occurrence and an 88-byte triangle import. An
independent weak owner remains live even though that item's worker reservation
has retired. A separate stream test proves that every publication lane still
accepts output after cancellation.

The stream now owns one detachable payload containing queued geometry, priority
coverage, staged imports and the cold/warm persistence journals. Its existing
mutex serializes all payload publication, consumption and cancellation.
Cancellation detaches that owner without allocation, closes subsequent writers,
and destroys the payload after releasing the stream mutex. Fixed progress and
coverage facts remain readable. Already drained geometry and claimed imports
retain their independent owners. Ordinary producer retirement preserves a
completed consumer-owned stream in both the coordinator handle and GED adapter.

The shutdown probe's storage owner queries coordinator accounting on release.
The first implementation deadlocks because shutdown still requests cancellation
under the coordinator mutex. Shutdown now closes admission, copies each active
job under the lock and cancels it outside, then pops and retires queued items
outside the lock. `source_job_retire_item` is the sole decrement/terminal
publication helper for both workers and shutdown; callers first retire queue
ownership or the active reservation. Admission-constrained stream cleanup also
runs outside the queue lock. These are refinements to the dated writer inventory
above, with no new scheduling phase or independent cancellation owner.

The production tests cover late publication, discarded journals, staged-owner
expiry with a healthy sibling, reentrant cleanup, successful staging transfer,
claimed ownership, cooperative shutdown and partial pool construction. The
[readiness checkpoint](obol_20260919_readiness_history.md#september-7-cancelled-stream-ownership)
owns negative/fixed runs, measured cleanup cost and graphical qualification.
Stream/container lifetime and lock mechanics are outside the TLA+ model; the
accepted cancellation outcomes and formal baseline are unchanged. This does
not qualify batch-held detached sources/databases/callback contexts, allocation
failure inside merge, or general cancellation latency and memory bounds.

### CXX-SOURCE-002 — submission allocates after consuming caller ownership

**Conforming at the tested submission transaction boundary.** The old method
claimed that every fallible allocation preceded resource transfer, but it
allocated deque storage and the returned job handle afterward.  An exception
could leave a partial queued batch and no client handle, with caller fields
already cleared.  GED's caller updates its ownership flags only after `submit`
returns, so throwing across this boundary also bypassed that reconciliation.

The deterministic regression injects `std::bad_alloc` at queue growth after
one appended item and at handle allocation for normal and constrained batches.
The instrumented predecessor fails all three cases: caller ownership is lost;
the two normal paths leave queue sizes 2/1 and 3/1 relative to the existing
sentinel work.  The constrained path also cancels the caller's stream before
throwing.  This is controlled exception injection, not a physical OOM test.

Submission now allocates an unattached handle first, stages queue entries under
the publication lock, and rolls back only its own suffix on exception.  All
resource transfers and handle attachment follow successful preparation.
Allocation failure returns an empty handle with caller ownership preserved.
The test verifies unchanged source/database/stream/callback handles, uncancelled
streams, the pre-existing queue and active reservation, no rejected callbacks,
successful retry and final reservation retirement.  Existing admission and
normal completion semantics remain unchanged.

The allocation mechanics are outside the current TLA+ models.  The downstream
lifecycle and deferred-autoview models check their existing control contracts;
they do not prove this exception-safety guarantee.  The executable fault matrix
and inspected commit boundary provide that evidence, as recorded in
[risk coverage](tla/RISK_COVERAGE.md).  Accepted-batch and downstream control
semantics are unchanged.

### CXX-SOURCE-003 — failed worker construction leaks threads and escapes drawing

**Conforming at the tested construction and GED start boundary.** The former
constructor allocated its private state through a raw member pointer.  If
starting a later worker threw, neither the state nor previously started workers
were retired.  The injected failure after the first worker left two process
threads against a baseline of one.  The independent Linux `/proc/self/task`
oracle detects that leak without reading the failed pool's state.

Construction now prepares an RAII candidate, allocates its fixed slots before
starting threads and joins every started thread on exception.  The public
singleton still reports construction failure by exception and can be retried.
The regression checks three failed attempts each return to the original thread
count, then successfully constructs exactly the requested two workers.

A second production regression found that this exception escaped GED's draw
transaction.  The existing deferred-start boundary now catches allocation and
worker-start errors before transferring request ownership.  The strengthened
GED fixture checks retained bounded proxy coverage, no pending source producer,
and successful erase/redraw into its four real occurrences after the fault is
removed.  Its saved negative executable fails at the escaped-exception check.

### CXX-SOURCE-004 — shutdown cannot cancel work already removed from the queue

**Conforming at the tested cooperative shutdown boundary.** The old destructor
set cancellation only on queued jobs.  It could not reach an active job, and
setting a queued job's flag did not notify the stream read by an in-flight
producer.  The shutdown regression retains client interest across actual
coordinator destruction: one active callback reserves the entire allowance,
while a separate two-item job stays queued.  Its cancellation-wait deadline
expired with the old code and the post-destruction observer rejected retirement.

Each fixed worker slot now records its active job under the queue mutex.
The shared stop path cancels jobs in both locations, including their streams,
then joins outside the mutex.  Cancellation visits each immutable job item
vector once even when many queue entries or workers reference the same job.
The regression requires the active callback to observe cancellation, both jobs
to become CANCELLED, no queued callback to execute, and callback contexts to be
released after their interest handles are dropped.  It tests source shutdown
rather than implicitly cancelling local handles before the destructor starts.

The [worker-pool readiness entry](obol_20260919_readiness_history.md#september-5-source-worker-pool-construction-and-shutdown)
owns binaries and results for both findings.  Linux thread enumeration,
injected startup failure and a cooperative callback under controlled admission
pressure do not prove native portability, physical OOM behavior, real geometry
cancellation latency or shared-stack sanitizer cleanliness.  Thread creation
and join mechanics remain outside TLA+; existing downstream control contracts
are unchanged.

## Source notification identity

### CXX-SOURCE-001 — an unkeyed publication flag suppresses distinct inputs

**Conforming at the source-evidence boundary.** The original
`publishPendingLodSourceRevision()` compared current source signatures with the
last submitted snapshot, but its
`lodPendingSourceRevisionPublished` Boolean remembers only that *some* change
was observed.  Before submission consumed that change, the next notification
returned early even if its source revision or evidence domain was different.
`submitLodRequestsIfNeeded()` also used that flag to suppress revision
publication without checking which signature/domain was previously published.

An isolated executable built from the existing visibility-census fixture and
linked to the pre-fix shared libraries demonstrates this sequence without GUI
timing or worker scheduling:

1. Finish the baseline submission.
2. Hide an occurrence and notify the controller.
3. Repeat the notification with no source mutation.
4. Restore visibility and notify before another submission pass.
5. Append a second occurrence and notify, still before submission.

The source reports one changed entry for each real mutation.  Its visibility
journal advances from 3 to 5 and inventory from 2 to 3.  Controller admission
visibility reports `3 -> 4 -> 4 -> 4`; inventory remains `2 -> 2`.  The duplicate
notification is correctly inert, but two distinct semantic changes do not
receive their current admission evidence.  This contradicts the complete
revision contract; the probe does not by itself establish a rendered image
failure or explain an earlier GUI timeout.

Artifacts are in `/tmp/obol-s1-publication-probe-20260905`: instrumented fixture,
`replay.py`, build/run logs, and `binary-identities.txt`.  The compile/link
recipe is `/tmp/obol-s1-probe-commands.txt`.  libBObol is `ce8d39fa...`, Obol
`f281b59b...`, and OSMesa `097aac27...`.  The final probe checks mutation
preconditions and exits 1 at the new assertion.  This retained before-fix probe
did not modify the normal test executable or production sources.

`BObolLodSourceEvidence`, in independently compiled `lod_source_evidence.cpp`,
now owns exact observed and submitted source snapshots.  Notifications and
submission share its `observe` decision; only a successful submission of the
latest observed snapshot may consume it.  Inventory/identity changes publish
inventory, visibility-only changes publish visibility, and duplicate delivery
or consumption publishes neither.  The unkeyed flag and competing controller
publication decisions are deleted.

The submitted baseline survives intervening observations for sparse source
deltas.  The two handles retain at most two immutable vectors of source-contract
signatures and share storage after consumption; they contain no occurrence
records or geometry.  Service/generation retirement releases both.

The production-owner unit test and expanded real-controller visibility fixture
cover two pending visibility changes, both visibility/inventory orderings,
duplicates before/after submission, source replacement/removal/reappearance,
notification-free submission and empty-source retirement.  The owner test also
rejects stale consumption and checks every source-identity field.  The
`ObolSourceEvidence` component and lifecycle/terminal conformance mappings cover
the observation/consumption edge.  The final build, all 14 affected CTests
and the complete 45-model TLC suite pass.  A mutant restoring the old latch
violates `ObservedInputIsDelivered` after two pending inventory observations.
The [dated readiness entry](obol_20260919_readiness_history.md#september-5-source-observation-and-consumption)
owns exact binaries, logs, state counts and graphical results.

All four OSMesa Generic Twin cold/warm shaded/wire rows pass.  The System GL
Generic Twin and both warm Lucy runs pass the strict control trace but retain
full-row HUD, retained-cut or visual-quality failures.  Those failures remain
release work; this focused identity closure does not establish regression-free
graphical qualification or complete S2.  Their causal relationship to this
change has not been established.

## Static-quality writer inventory

This audit covers the existing `BObolLodStaticQualityTrial` value and all its
production mutation callers in the controller/coordinator.  It retains five
private fields: state, sampled ceiling, baseline point threshold, constraint,
and acceptance.  The extraction adds scalar input/decision values, with no
new retained phase, collection, scene copy or scheduler.  The trial remains
`IDLE`, `PROBING`, `RECONCILING`, `REJECTED`, or `ACCEPTED`.

The September 5 lifecycle audit enumerates **19 direct mutation sites**:
four starts, five completion/rejection sites, and ten retirement/restore
sites.  `static-writers.json` in the
[configuration checkpoint](obol_20260919_readiness_history.md#september-5-static-quality-configuration-lifecycle)
records each file, caller and mutation.  The frame-deadline setter retires
through the existing policy-revision writer; it adds no direct trial mutation.

| Transition family | Writer and enabling evidence | Completion, supersession and bound |
|---|---|---|
| Start a trial (4 sites) | `begin` from `advanceStableLodReducerIfReady`, two branches in `recordCompletedRenderTiming`, and `submitLodRequests`: static overscan, accepted one-pixel frame, one-pixel trial, or protected-minimum retry | An existing non-idle trial prevents restart. Cut advancement, exact acceptance or a measured constraint owns completion; starts remain caller-selected policy for S3 extraction. |
| Observe a completed frame | `completeFrame` owns preparation replay, one steady sample per candidate, richer/fractional selection, exact fractional acceptance and predicted denial | Transient preparation is not sustainable capacity. Live production/presentation/interaction owners defer sampling. A selected successor prevents ordinary planning in the same callback. Predictions read fixed-size cut-cost tables; no occurrence traversal enters the trial. |
| Transfer a measured population budget | `BObolLodCoordinator::completeStaticQualityFrame` arms the existing presentation handoff, starts one bounded submission pass, and advances CAPACITY once | Repeated callbacks cannot bypass either active submission or pending handoff. The renderer guard remains until the existing certified allocation handoff; an unavailable single-cut prediction does not create a rejection certificate. |
| Reject an interrupted candidate | `notePresentationRenderInterrupted` calls `reject` after the existing typed terminal-capacity predicate selects a valid stamped constraint | RECONCILING retains the safe presentation until occurrence reconciliation. This interrupt-side selection remains outside the completed-frame extraction. |
| Constrain a protected minimum or saturated allocation | `advanceStableLodReducerIfReady` and `submitLodRequests` call `constrainPresented` with the measured/predicted operation's constraint | REJECTED retains only the narrow revision-bound protected-minimum ceiling exception. Existing constraint tests cover this boundary; selection remains caller policy. |
| Complete occurrence handoff | `submitLodRequests` calls `completeOccurrenceHandoff` only after the presentation policy reports a reconciled allocation | RECONCILING becomes REJECTED once; PROBING releases its guard for one deadline-bounded population frame. The protected-minimum exception retains its certified guard. |
| Retire or restore the value (10 sites) | `setLodService`, `setLodAutoSubmit`, `setLodPolicyRevision`, `advanceLodViewRevision`, `advanceLodPolicyRevision`, `beginLodInteraction`, `setLodFrameRateTargets`, `advanceStableLodReducerIfReady`, `resetRendererPerformanceEvidence`, and `retireAutomaticLodControl` call reset/deactivate/restore | Reset clears the complete value; restored guarded performance evidence re-enters PROBING. Changed policy invalidates; duplicate configuration preserves work. Internal CONTINUE_STATIC_QUALITY preserves the trial. Acceptance/constraint compares the renderer-capacity problem: recording the same frame's CAPACITY edge does not invalidate its own proof. |

Existing `test_static_quality_policy`, capacity-domain, readiness, policy
retirement and occurrence-handoff tests cover the retained lifecycle.
`test_static_quality_completed_frames` now drives the exact compiled policy
which the render callback calls.  `test_static_quality_frame_reconciliation`
also drives the real coordinator effect boundary; the two-progressive result
fixture exercises the production cost provider.  These are focused C++
regressions, not a claim that every GUI event ordering or resource magnitude is
qualified.  Start/interruption policy is inventoried but remains caller-owned
for S3.  The remaining source/service release, GED, host and presentation
writer audit stays in S1's existing backlog.

### CXX-STATIC-002 — unchanged configuration requests another capacity frame

**Conforming at the tested configuration boundary.** Reapplying
`setLodAutoSubmit(TRUE)` to a settled constrained view unconditionally reset
the static trial and synchronized automatic control.  The constrained-source
fixture observed a new render request and obligation mask 32, although its
constraint mask and allocation serial were unchanged.  This counterexample
proves unnecessary work; it does not prove loss of a static rejection
certificate in that particular fixture.

An identical normalized setting now returns before reset or synchronization.
The regression requires no render request or control obligation and preserves
the completed constraint mask and current allocation serial.  Changed settings
retain their existing activation/retirement path.

### CXX-STATIC-003 — changed deadline leaves settled policy evidence current

**Conforming at the tested policy-invalidation boundary.** The lower-level
frame-deadline setter changed endpoint allowances and quality history but did
not advance policy or request work.  After the first fix above, the same
constrained fixture failed at the next assertion: policy remained 3 and no
frame was requested when the stable deadline increased to 500 ms.

A changed deadline now uses the ordinary policy-revision invalidation path,
marks progressive work pending, and requests a capacity frame.  The existing
same-value guard remains.  The regression verifies new policy and frame
authority, then drives normal settling and requires an exact three-occurrence
presentation with more actual triangles than before the deadline change.
It uses controlled traversal cost and real rendering; it is not a latency
benchmark.

Both saved negative builds and final passes are indexed by the
[configuration readiness entry](obol_20260919_readiness_history.md#september-5-static-quality-configuration-lifecycle).
`ObolStaticQuality` maps changed policy to a new input epoch and unchanged
configuration to stuttering.  Numeric deadline changes are outside that
model; `tla/conformance.json` classifies these as production regressions.
No state, model transition, deadline default or validator threshold was added
or relaxed.  These fixes do not close the unrelated graphical reproductions
or complete S1/S2.

### CXX-STATIC-001 — unavailable prediction releases an effective guard

**Conforming at the tested completed-frame decision boundary.** The prior
callback used `singleCadProgressiveCutWithinRenderCost` as an ordinal result.
That provider deliberately returns `-1` for a population with more than one
progressive occurrence.  The old fallback published that value as the
renderer ceiling, removing an effective guard before an occurrence allocation
represented the measured budget.

The first extraction preserved that branch.  Its deterministic probe failed
with `unavailable cost prediction released a renderer guard`; the retained
pre-fix policy, test and executable are under
`/tmp/obol-static-frame-owner-20260905/before-prediction-fix`.
The final production decision selects `HANDOFF_POPULATION`, carrying the
measured current cost and allowed budget to the existing allocator while
preserving PROBING and the renderer guard.  It does not invent a capacity
rejection for missing prediction data.  Known predicted misses still select
RECONCILE_CUT with their explicit candidate constraint.

The former callback's sample, fractional acceptance, rejection and cut-step
decision branches are deleted.  `lod_static_quality.cpp` compiles independently
in the library and policy tests; the callback now supplies scalar evidence and
performs the selected renderer/frame-request effects.

`ObolStaticQuality` covers completed-frame commit, bounded acceptance and
rejected-cut reconciliation.  `ObolCapacityPresentationHandoff` covers the
canonical reconciliation request and subsequent certified guard release.
Their abstract transitions are unchanged.  Neither model encodes the numeric
single-occurrence query or its sentinel: the production regressions supply
that refinement evidence.  `tla/conformance.json` records them as regressions,
not exhaustive stepwise proofs.  The [dated readiness entry](obol_20260919_readiness_history.md#september-5-static-completed-frame-ownership)
owns binary hashes, verification scope and graphical results.  This finding
does not establish the cause of earlier graphical failures or complete S2.

## Submission retarget ownership

### CXX-CURSOR-001 — an untouched pass acquires unnecessary rescan debt

**Conforming at the tested input-retarget boundary.** The Generic Twin System
GL wire trace exposed a planning A/B/A window at unchanged revisions:
SUBMISSION -> SUBMISSION + SUBMISSION_RESCAN -> SUBMISSION, with offsets
0 -> 2048 -> 0.  Rescan completion can be finite; the trace alone did not prove
an infinite loop.  The caller audit found that `submitLodRequestsIfNeeded`
requested a successor rescan for every changed input while a pass was active,
including a newly opened pass which had consumed no predecessor entries.

`test_compact_input_retarget` drives that boundary through the actual automatic
controller and a shared 3,014-occurrence compact fixture, with all occurrences
outside the view to exclude provider work.  After coverage opens a fresh
quality pass, a policy change previously left another active successor at
cursor 0/0.  The retained automatic probe fails before the fix and passes
after it.  Earlier manual-controller probes did not exercise this automatic
lifecycle and are not the counterexample.

| Retarget condition | Required production effect | Executable evidence |
|---|---|---|
| Ordinary cursor at source 0, entry 0 | Rebuild its source-local plan and pass annotations against the new input; preserve any independently owed rescan without inventing one | Policy retarget preserves the fixture's exact-empty census and completes without another pass; scale retarget visits the new census exactly once |
| A source or entry prefix has been consumed | Keep the cursor and request one complete rescan for the predecessor prefix | A policy change during the scale pass consumes the remaining suffix, then visits the full population once: exactly twice the population across the two passes |
| A selective append plan extends or invalidates its scope | Preserve the existing delta-extension/full-rescan rules | Existing append, source replacement and exact visibility regressions remain passing |

The fix is at the existing caller decision.  It introduces no state, revision,
counter, scheduler, renderer mode or scene mirror.  Reinitializing the untouched
source plan is necessary: an old sparse plan must not exclude entries required
by the changed policy.  A consumed prefix continues through the existing
`BObolLodSubmissionPass` owner and bounded action; no new completion mechanism
is added.  Other coverage, demand, result, capacity, presentation and lifecycle
callers of `beginFresh`/`retire` remain in the S1 successor audit.

`ObolSubmissionPass` treats bounded passes atomically; it covers rescan and
demand ownership, while these production tests check concrete offsets and
visit counts.  Its transition relation is unchanged and the catalog now names
the caller and demand-policy owners.  The [dated readiness entry](obol_20260919_readiness_history.md#september-5-submission-retarget-prefix)
owns verification scope, binaries and graphical results.  The original trace
and `/tmp/obol-submission-rescan-20260905/before-auto-fix` retain the failed
boundary; the latter includes the test executable and its old libBObol.

### Submission writer and successor inventory

The September 5 caller audit classifies all **51 direct production mutations**
of `lodSubmissionPass`: 27 in `view_controller.cpp`, four in
`view_controller_progressive.cpp`, 18 in `view_controller_render.cpp`, and two
in `view_lod_coordinator_state_private.h`.  No escaping reference or assignment
to the member was found; writes use `BObolLodSubmissionPass` transitions.
`/tmp/obol-submission-owner-audit-20260905/submission-writers.json` retains each
method, source line, category and source hash.  This is a finite caller
inventory, not proof that every caller composition is correct.

| Caller family | Enabling fact and next owner | Completion, supersession and bound |
|---|---|---|
| Service replacement, generation cancellation, empty source set, automatic-policy retirement | Explicit lifecycle/input edge | Retire the cursor and relevant demand/source state. Provider cancellation and lease release retain their separate lifecycle contracts. |
| Source/view/policy reconciliation | Changed exact source signature or view/policy revision | Source replacement restarts the plan; exact deltas retain their scope; an untouched ordinary plan is rebuilt and a consumed prefix owes one rescan. `CXX-CURSOR-001` checks the production cursor/visit boundary. |
| Deferred discovery rescan | Incomplete inventory pauses a completed selective pass; final inventory enables `beginPendingRescan` | Preserve one owed full pass while paused; consume it once the provider closes. New append deltas extend the current source plan rather than restarting its consumed prefix. |
| Ordinary demand and explicit submission | `demandPassRequired` sees outstanding demand with no stronger owner; a direct public submission explicitly requests a pass | Bounded source/entry progression. Only a complete ordinary, nonselective, nonrepair, nonallocation pass without rescan debt retires demand. Explicit submissions remain available to manual callers. |
| Coverage completion/resumption | A required census is active, or a completed census reports missing coverage/deferred refinement | Missing coverage yields to resident growth or a capacity frame; completed coverage enables one quality successor. Census counters and current-view coverage are distinct from physical demand. |
| Importance census | Quiet-transition importance obligation completes its ordinary census | Retire the obligation before requesting one occurrence allocation; small-source completion has the same successor as bounded coverage. |
| Mechanical/allocation completion | Source index reaches the end, or a retained allocation slice advances | Partial allocation remains active; completed applied cuts yield to their exact frame. Stale unpublished allocation restarts only against changed input. Completed retained allocation cannot spend its allowance twice. |
| Result/resident availability | Authenticated current-demand retry, richer resident prefix, or allocation-invalidating compaction | Results may request demand replay. Resident growth drains bounded suffix work, then consumes its growth edge into coverage or allocation; in-flight producers own the wait. Compaction repair requires actual allocation invalidation. |
| Resident-memory retry | Retryable denial and a new service admission revision, with submission and stronger work quiescent | Record the attempted admission revision before scheduling; completion or semantic invalidation retires its retry annotation. An unchanged admission revision cannot repeatedly retry. |
| Structural/presentation repair | Exact unresolved occurrence frontier, normal-style change, or missing current presentation binding | Existing repair/frontier owners bound structural work. `beginCadPresentationRepairPass` retains ordinary demand; its callers are `setNormalStyle` and exact-payload recovery. Current bindings request an unchanged frame instead of another source mutation. |
| Capacity candidate/reconciliation | Existing capacity decision requests allocation, or a pending candidate/handoff becomes runnable | `beginSceneWideCapacitySubmission` clears selective scope and predecessor annotations. `restartLodCapacitySubmission` also restores an owed coverage census. Sample/candidate progress remains owned by the bounded capacity search. |
| Completed/interrupted presentation | Deadline recovery, consumed handoff frame, retired publication barrier, or confirmed recovery with a stale allocation | A certified current allocation needs no replacement pass. A retired barrier cannot restart planning while capacity owns the successor. Reconciliation keeps its measured budget and guard until occurrence allocation is certified. |
| Static/headroom/point quality | Existing terminal/static/point policy selects a richer tier, measured allowance, structural preload or cheaper recovery | Passes consume a selected budget/threshold/cut/frontier. Search brackets, one-shot headroom, discrete-trial permission and exact frames bound their successors; another callback alone is insufficient. Static start/interrupt and point numeric branches remain caller-side extraction work. |

The indirect capacity entries were also inspected: interrupted-frame decisions,
unmeasurable-frame retirement, completed capacity samples, terminal allocation,
completed-pass calibration, delayed handoff, and exact visibility reallocation.
They funnel through the existing capacity helpers above.  The coordinator's
`completeStaticQualityFrame` starts a pass only for its selected reconciliation
or population-handoff result.  This inventory does not introduce a universal
start/restart dispatcher; the selection belongs with each existing policy.

The related demand and planning values have these writers:

| Value | Request/replacement | Consumption/retirement |
|---|---|---|
| Full-view demand | `refreshForViewRevision`, `refreshForPolicyRevision`, or authenticated result quality debt | `completeDemandRefresh` at ordinary pass completion; service/policy retirement resets it; source coverage and presentation repair preserve it |
| Importance census | Quiet interaction completion arms it together with demand | Completed ordinary census requests allocation; new camera/interaction/source input retires the old importance request |
| Resident-admission retry | A previously untried admission revision with retryable denial | Completed pass, changed input, generation cancellation, or full automatic retirement |
| Exact visibility reallocation | Exact visibility delta, or an allocation's unresolved projection frontier | Consume only after source/presentation/repair/capacity prerequisites yield, then start full allocation; incompatible source input or full retirement resets it |

Existing `test_submission_pass`, `test_view_demand_policy`,
`test_coverage_policy`, planning/availability/capacity policy tests and the
production controller fixtures exercise these owners.  The actual compact
append/visibility, source churn, normal-policy, retained spatial repair and
policy-disable tests cover effects at selected seams.  The composed models
remain separate evidence; the primitive four-state test alone does not prove
all 51 caller effects.

`CXX-DEMAND-002` below closes the pending-demand loss found at source
reconciliation.  The importance request still retires with its old source
domain.  `test_view_controller_source_demand_projection_replacement` now checks
ample and constrained cases using actual software traversal.  Its constrained
case puts retained visible entries beyond the consumed census prefix, then
checks current projection, increasing allocated detail with scale, redistribution
from retained occurrences, actual presented triangles/flat normals and finite
retirement.  A continuous source-transaction journal checks conformance states
and identifies the exact current frame preceding the final allocation.
The [constrained-source readiness entry](obol_20260919_readiness_history.md#september-5-shared-minimum-coverage-under-source-replacement)
records evidence and limits.  This closes that effect proof at the tested
boundary; physical extraction and the remaining state-family audit are open.

### CXX-COVERAGE-001 — shared repair rejects a richer prefix as minimum coverage

**Conforming at the tested shared-prefix repair boundary.** A new visible
occurrence can reuse an already resident progressive asset.  Retained admission
selected its minimum prefix, but structural repair did not: reuse tried the
richer requested/resident cut, rejected its cost, and skipped the occurrence.
The controller could then settle on a terminal box without trying the affordable
minimum mesh.  The regression reproduces this with a 449-unit repair remainder
and a shared minimum costing 74 units.

Structural repair now uses the existing minimum-prefix selection condition.
The shared binding still charges exact per-occurrence cost to the same budget;
provider, residency, preparation and completed-frame limits are unchanged.  No
new owner, state, mode or scheduler was added.  The exact strengthened regression
fails against the saved predecessor library and passes against the fixed one.
It requires three actually rendered meshes and a genuinely constrained allocation,
not merely a nonzero plan certificate or a terminal controller state.

`ObolStructuralFrontierOwnership` abstracts the measured rejection as an outcome;
it does not distinguish numeric prefix costs.  The C++ integration test checks
that refinement boundary.  The model transition relation is unchanged.

### CXX-DEMAND-001 — presentation repair consumes ordinary demand

**Conforming at the tested repair-completion boundary.** A normals change uses
`beginCadPresentationRepairPass` to replace the shared source cursor.  That
pass can refresh an existing CAD payload at its retained cut and return before
updating richer demand.  The caller previously classified only structural
repair and retained allocation as insufficient demand-completion proofs.
Consequently a complete normals repair could clear the ordinary demand flag.

The automatic-controller regression extends the existing 3,014-occurrence
fixture.  A scale change consumes a real demand prefix; a normal-style change
then replaces it with a repair pass.  Before the fix, all repair entries are
visited and demand disappears.  After the fix, repair completion preserves
demand, one complete ordinary successor visits the population, and both demand
and the cursor retire.  All occurrences are outside the view, so the control
counterexample needs no worker tasks, cache writes or GUI timing.

The caller now includes presentation repair in the existing deferred-demand
guard of `completeDemandRefresh`.  No state, revision, pass mode or scheduler
is added.  Existing policy tests cover the guard; the added production test
checks that the normals caller supplies it.  `ObolSubmissionPass` already
preserves demand through presentation completion; its transition relation
is unchanged.  The catalog mapping is regression evidence, not a stepwise
proof that its abstract presentation state is the concrete repair cursor.

`/tmp/obol-submission-owner-audit-20260905/before-fix` retains the failing test
source/executable and matching old libBObol.  The
[dated readiness entry](obol_20260919_readiness_history.md#september-5-demand-preservation-through-presentation-repair)
owns build, model and graphical evidence.  This finding does not establish
the cause of earlier same-cut normals or other graphical failures.

### CXX-DEMAND-002 — source reconciliation discards pending demand

**Conforming at the tested source/demand boundary.** Source reconciliation
used to clear full-view demand whenever planning inputs changed.  An inventory
append or visibility edit during a consumed demand prefix could therefore
discard the obligation before a subsequent normals repair replaced its cursor.
Completing that repair then left no ordinary pass to evaluate the current view.

The automatic-controller regression extends the same compact fixture with
both interleavings.  Inventory append increases the population to 3,015;
visibility mutation hides one entry.  Before the fix, each complete repair
loses demand.  Afterward, repair preserves it and exactly one complete ordinary
successor visits the revised selected population and retires demand/cursor.
The visibility edit occurs during an active nonselective pass and exercises
the full-census fallback.  The separate idle visibility regressions continue
to prove O(delta) behavior.  These fixtures use off-view geometry and prove
control effects, not numeric quality or rendered allocation.

Preserving the flag alone would conflate owed demand with permission to refine
during source coverage.  `BObolLodCoveragePolicy` now replaces its active Boolean
with a private `NONE`/`SOURCE`/`VIEW` census value.  Source coverage retains
minimum-coverage priority while demand remains pending.  A view census may
refine previously covered geometry; camera changes cannot bypass incomplete
source coverage.  Recreating convergence-census keys after a camera change
preserves the already selected owner.  There is no additional independent
latch, scene mirror or scheduler.  The unused coverage `deactivate` entry was
removed.

Demand completion now excludes source coverage as well as repair and retained
allocation.  The deferred-coverage demand census uses the same exclusion,
including presentation repair.  Coordinator tests exercise admission priority,
repeated view invalidation and post-coverage retirement.  The existing one-mesh
camera fixture now drains its demand and threshold-coverage setup before its
unchanged-request and scale-only assertions; those assertions are retained.

`ObolSubmissionPass` now includes initial coverage/demand overlap, preserves the
demand wakeup at coverage completion, and checks that only the demand pass can
consume demand.  A paused rescan is an existing producer obligation while
inventory publication is pending.  The isolated old-completion model violates
`DemandHasProducer`; the corrected component passes.  The model's coverage
state represents source priority, while C++ tests distinguish source and view
censuses.  This remains regression mapping rather than a stepwise refinement
proof of every concrete cursor mode.

`/tmp/obol-source-demand-supersession-20260905` retains both failing production
probes, the test/executable and matching old library, the isolated model
counterexample and implementation patch.  The
[dated readiness entry](obol_20260919_readiness_history.md#september-5-demand-preservation-through-source-coverage)
owns final qualification evidence and its limits.

### CXX-TIMING-001 — completed traversal resamples delivery time

**Conforming at the tested duration-consumption boundary.**
`finishPresentationRenderTiming` accepted a measured traversal duration but
discarded it for exact frames.  It called `completeRenderTiming`, which sampled
the clock again.  Thus delivery and post-traversal classification could change
the duration paired with the completed population.  The automatic deadline
contract regression supplies a 41 ms traversal with 25 ms delivery delay; the
old path reports approximately 66 ms instead of the supplied 41 ms.

`recordCompletedRenderTiming` is now the single consumer.  The existing
start-only completion API obtains elapsed time once; explicit completed-frame
delivery and `renderToImage` preserve the measurement already taken at traversal
completion.  Start identity still qualifies exact-frame retirement, and
incomplete traversal keeps the existing interruption path.  The change adds no
state, clock injection, callback, mode or scheduler.

The existing deadline test checks exact duration equality for late and on-time
frames and retains its incomplete-frame assertion.  `ObolCadTimingEvidence`
models coherent outcome classes independently of unrelated host timing;
preservation of numeric nanoseconds is a C++ regression claim beyond those
classes.  Its transition relation is unchanged.  The
[readiness entry](obol_20260919_readiness_history.md#september-5-completed-traversal-measurement)
owns graphical qualification and artifact identities.  This repair supplies the
controlled-duration effect boundary needed for the pending constrained source
allocation proof; it does not close that proof.

### CXX-TIMING-002 — pose validation omits interrupted frames

The Generic Twin validator inferred responsive rotation from the maximum
completed-render duration. The current and earlier System warm-wire reports
both contain a new interrupted frame during rotation: approximately 50.14 ms
against the 50 ms hard deadline. A subsequent completed frame below 52.5 ms
does not invalidate that miss or its presentation-only recovery ceiling.
The current trace also records `render-deadline` at that boundary.

The matrix now uses `lod_pose_contract.jq`, with the existing 5% completed-frame
tolerance, one-pixel criterion and 95% retention threshold. It distinguishes
new interruption counts within the gesture from historical counts and requires
their measured duration to meet the hard deadline. Evidence after the held
checkpoint cannot justify a ceiling already visible there. Responsive retention
uses exact presented faces or lines; disabled deep face diagnostics previously
made the old population comparison vacuous in both modes.

The standalone contract cases cover valid recovery, missing or decreasing
counters, missing or inconsistent deadlines, stale/later pressure, incorrect
scale classification and lost shaded/wire population. Offline application to
retained reports preserves the original failed logs and records a separate
verdict. This is an observation repair; no production policy writer, threshold,
renderer behavior, formal transition relation or accepted TLC baseline changes.
The [readiness checkpoint](obol_20260919_readiness_history.md#september-7-pose-deadline-observation)
owns qualification and its remaining limits.

## Remaining audit findings

`CXX-SOURCE-001`, `CXX-STATIC-001`, `CXX-CURSOR-001`, `CXX-DEMAND-001`,
`CXX-DEMAND-002`, and `CXX-TIMING-001` are closed at their tested source-evidence,
completed-frame, input-retarget, repair/source demand-completion and
duration-consumption boundaries.  Earlier
named findings retain their focused closure evidence.  S1's physical
ownership/acceptance audit and the recorded graphical failures remain open.
Independent numeric, geometry, image, concurrency, durability, performance and
platform gates remain necessary regardless of the ledger status.

## Change gate

Any change to result routing, shared producer lifetime, revision identity, or
terminal outcome must update together:

1. the focused component model and applicable composition;
2. its production typed value/reducer rather than a parallel latch;
3. the executable table or lifecycle matrix at the real asynchronous seam;
4. `tla/models.json`, `tla/conformance.json`, ownership/test metadata, and this
   ledger; and
5. the focused TLC baseline after review, plus the complete formal-suite gate.

Passing TLC cannot close numeric quality, geometry validity, memory magnitude,
latency, renderer semantics, or filesystem crash consistency.  Those remain
independent release evidence even when their control ownership is modeled.

### CXX-SOURCE-019 — live cached realization exposes private construction

The public wire/mesh APIs and realization action previously ran their builders
on observable sources. A retained public-API probe reported incomplete state in
12/14 wire callbacks, 12/14 shaded-mesh callbacks, 11/48 evaluated-wire callbacks
and 14/48 evaluated-points callbacks. The final production regression also rejects
the preceding libBObol for terminal state/child retirement. These are library
publication defects, independent of renderer timing.

Live calls now create an owned configuration-only template, borrow the database
on the owning thread, and construct geometry privately. The existing prepared
source publication commits the candidate index or children, bounds, realized
record, evaluated role, auxiliary metadata and compiled/path retirement before
notifications. Installation reapplies the live registry's sparse overrides.
A failed provider publishes FAILED coherently while retaining the preceding
geometry and realized revisions; allocation failure before publication leaves
that drawing unchanged. Callback failure leaves a complete replacement.

Workers explicitly call the private construction entry and retain their existing
stream/terminal handoff. Public stream overloads retain their NULL semantics;
stream entries remain incremental publications owned by the stream protocol.
There is no extra worker registry, persistent source mode or copy of the old
scene. The obsolete cache preservation flag and duplicated private terminal
assignments are removed. Private templates omit invalidation sensors and own
partial copies; evaluated adapters scope temporary plot/face-set and shape
ownership. Wider traversal/provider cleanup and resource bounds remain open.

The allocation sweep also reproduced a dependency defect: `SoSFEnum::setEnums`
freed its current mappings before acquiring replacement storage, allowing a
failed field copy to double-free them. It now owns both arrays before replacing
either mapping and handles self-assignment. The final enum test rejects the
preceding Obol and checks both failed allocations and a complete replacement.

The production cases cover live/public/action calls, explicit NULL and non-NULL
streams, repeated assemblies, prior compact/compiled drawing, sparse selection
and appearance overrides, auxiliary/nested owners, placement, disabled fields,
provider failures, path retirement, retry and throwing observers. Evaluated
points use pseudo-random sampling: their observers require complete indexed
triangles and corresponding normals, plus exact terminal/owner state, rather
than equality to an independent sample. Named cases retain the exhaustive
sampled-point sweeps; routine CTest uses three allocation boundaries for each
sampled case while sweeping the other cases completely.

The [dated evidence](obol_20260919_readiness_history.md#september-7-live-cached-realization)
owns final counts, binaries and graphical qualification in
`/tmp/obol-live-cached-realization-20260907`. This closes the live cached
publication boundary, not compact edit revision changes, auxiliary/prototype
publication, invalidation, broader ownership, resource or release qualification.

### CXX-SOURCE-020 — combination edit publication and cache lock recovery

A real repeated-instance edit moved the right occurrence from x=5 to x=10,
but retained the preceding exact maximum x=7 instead of x=12. The combination
refresh also changed the live source revision before building and recreated its
field sensors around that change. A failed preparation could therefore expose
an incomplete edit or lose revision invalidation.

Combination refresh now advances only the private candidate, discards its copied
bound, and uses the shared prepared publication for requested/realized revision,
geometry, bounds, auxiliary owner metadata and compiled/path retirement. The
source-mesh contract receives that same revision. Failed construction preserves
the complete preceding drawing; subsequent fallback/invalidation remains the
caller's responsibility. Stable handles, sparse selection and appearance survive
successful replacement. The live sensor teardown and rollback writes are removed.

The BoT allocation sweep exposed a second production defect: an exception inside
LoD asset cache storage left the draw-cache semaphore held. The next request
blocked in the same semaphore. Store/get now use scoped lock ownership; the
same sweep completes after the repair. The timed-out run, stack and preceding
library are retained. Other draw-cache entry points and full cache/provider
resource cleanup remain open.

The production test covers wire/shaded and BoT request paths, explicit/successor
revisions, disabled notifications, compiled drawings, auxiliary/nested owners,
failed providers, persistent allocation denial, retry and throwing observers.
It verifies actual occurrence placement, bounds, handles, overrides, source-mesh
currency and revision sensor behavior after failures. The
[dated evidence](obol_20260919_readiness_history.md#september-7-combination-edit-publication) owns counts and current GUI results.
The incremental leaf-edit transaction was the next source writer at this checkpoint.

### CXX-SOURCE-021 — incremental leaf edit publication

Leaf refresh previously exposed partially updated revision/geometry to observers.
It also retained a source profile certified before the edit and rejected indexed
selection paths such as `/pair.c/box.s@1`. The regression rejects the preceding
library, and the indexed-path counterexample is retained separately.

Refresh now prepares replacement records only for affected occurrences, new part
ownership, exact bounds and source/owned-shape metadata. It shares the existing
part publication and occurrence commit helpers with streaming replacement. It
commits these facts before notifying, preserves nested source owners, and retains
live sensors and the compiled object for its existing incremental rebuild. Failed
preparation preserves the preceding drawing; a throwing observer sees a complete
edit. Temporary fallback shapes have scoped ownership. Unchanged content reuses
immutable geometry, obsolete source profiles are invalidated, and source-mesh
request currency advances with the edit.

The production tests cover wire/mesh replacements and unchanged geometry,
explicit/successor revisions, allocation denial, missing-object rejection, retry,
throwing observers, disabled notifications, compiled geometry, sparse overrides
and unrelated shared occurrences. Indexed-path editing updates both occurrences
of the database primitive; exact bounds grow and shrink. Preparation allocation
counts agree with 32 and 4,096 unrelated occurrences. This measures counts, not
total bytes or latency; replacement of an extremal bound may scan the index.
The [dated evidence](obol_20260919_readiness_history.md#september-7-incremental-leaf-edit-publication)
owns measured counts and current graphical results. Auxiliary/prototype
publication, invalidation, broader provider/allocator/traversal ownership and
release qualification remain open. No formal relation changed.

### CXX-SOURCE-022 — auxiliary line publication and retirement

The preceding library dispatched 66 callbacks while inserting an auxiliary line
whose geometry and metadata were incomplete. Named replacement also mutated live
fields, and clearing auxiliaries exposed successive partial child populations.
The direct production tests reject that library. Separate regressions expose
missing structural revision advancement and doubled nested placement: a parent
translation of `(10, 20, 30)` produced a minimum of `(20, 40, 60)`.

Named lines now prepare their new geometry privately. The existing primary
publication's placement/child logic is extracted into `PreparedSourceChildren`
and shared with auxiliary replacement. It commits the complete child/path edit
before notifying. Primary geometry, compact handles, unrelated children and
nested owner metadata remain retained. Node names, primitive selection, edit
intent and notification settings survive replacement; obsolete geometry is not
copied into the candidate. Successful replacement retires paths to the old child,
as documented in the API, and ordinary controller calls advance scene structure.

Empty named-line/source input and complete auxiliary clearing use prepared
removal, including repeated child occurrences. A throwing observer does not
interrupt the remaining notification attempts. Nested auxiliary sources retain
placement metadata while inheriting the parent's transform exactly once; the
bounding-box regression verifies both insertion and replacement.

The [dated evidence](obol_20260919_readiness_history.md#september-7-auxiliary-line-publication-and-retirement)
owns allocation counts, binaries and graphical qualification. This closes the
line/child retirement boundary, not the enclosing nonempty auxiliary source's
configuration/terminal transaction, controller revision recovery after observer
exceptions, prototype realization, invalidation or broader resource/traversal
qualification. No formal relation changed.


### CXX-SOURCE-023 — auxiliary source configuration and publication

The preceding library attached a new named source before completing its
configuration and realized record. The current regression rejects it after
observing incomplete source/shape metadata in a 98-callback insertion. A separate
controller regression rejects its post-notification revision update. Allocation
denial during live source construction also exposed an Obol field sensor that
recorded attachment before auditor storage existed; constructor cleanup crashed.

Nonempty auxiliary sources now prepare a detached configuration, the named line,
retained owned-shape metadata and both source/parent child-path edits before
publication. `PreparedScalarValues` is shared with owned-shape preparation;
`PreparedAuxiliaryLine` reuses the existing placement/child publication. This adds
no persistent scene model. Named source routing identity and input sensors
survive updates; nested source owners remain independent. Explicit line display
overrides and per-primitive semantic state survive. Obsolete source bounds,
profiles and staged streams retire at commit. Advancing the existing population
epoch rejects older worker results even when configuration values are unchanged;
the parent source's compact handles remain valid.

A private, nonthrowing controller hook advances ordinary scene structural/frame
revisions after commit and before notifications. Preparation failure leaves both
revisions unchanged. Observer exceptions leave a complete commit and do not stop
the remaining notification attempts. Existing mutation-batch revision deferral
is preserved; reentrant graph mutation and arbitrary custom callbacks are not
qualified by these ordinary-context tests.

Obol now records a field sensor's attachment after `addAuditor` returns. The
regression exercises initial field storage and auditor-list growth, verifies
failed attachment leaves no registration, and checks retry/notification/detach.
The full new-source construction sweep also qualifies constructor unwinding.
`SoField::addAuditor` now rolls back its registration if a custom
`connectionStatusChanged` hook throws; a rejecting-hook regression verifies no
remaining auditor and exactly one delivery after retry. The final regression
rejects the intermediate attachment-only repair. Reentrant auditor-list mutation,
arbitrary hook side effects, connected engines and evaluation/deletion callbacks
remain in the wider dependency audit.

Twenty scenario definitions run directly and through the controller: **40 cases,
3,980 allocation positions, 3,782 preserved and 198 committed outcomes**. Nine
source cases cover insertion, replacement, default commands, display overrides,
disabled fields, alias paths, renamed lines, database rebinding and unchanged
configuration. They check complete recursive scene data, paths, nested owners,
compact handles, sensor recovery, retry and throwing observers. The
[dated evidence](obol_20260919_readiness_history.md#september-8-auxiliary-source-publication)
owns binaries and graphical results. Prototype publication, invalidation, other
writers and full resource/platform qualification remain open. Formal relations,
quality floors and renderer defaults are unchanged.


### CXX-POLICY-001 — live policy published before its control epoch

`ged_view_lod_policy_apply` copied a changed policy directly into the attachment
before `ged_draw_obol_view_lod_policy_changed` cleared view state and advanced the
policy revision. Entering the first controller operation therefore observed new
policy demand under the old epoch. An eager LoD-off draw followed by enable
exposed `NONTERMINAL_WITHOUT_PROGRESS` (mask 128). A settled partial-delivery or
adoption failure instead reopened its terminal state without an input revision.
The implicit `producer_progress` record was evidence of that out-of-order write.

`BObolViewController::setViewLodPolicy` owns the complete live operation. It
sanitizes and compares before mutation, computes the next revision using the
existing exhaustion rule, and enters the existing external-input scope before
publishing the attachment value. It retires obsolete selection/generation when
LoD is enabled, then applies the revision and successor through the existing
policy writer. Master OFF preserves retained geometry even when subordinate
representation switches remain enabled. Equivalent sanitized values and null
input are no-ops. `bv_lod_policy_equal` centralizes the existing value comparison
used by BV, GED and the controller.

GED delegates to this operation; its duplicate clear/revision sequence is
removed. Attachment-only writes remain for initialization/binding or a view
without a controller. Source realization policy and shared-source notification remain in
the GED backend after the coherent view operation. No scheduler, control facts,
trace predicates or formal relations were added or relaxed.

The existing `test_view_controller_policy_disable_retires_automatic_work` now
uses the production operation for live toggles. It checks retained payload
identity, retirement of automatic work, normalized no-op policy values, and one
valid revision-advancing external-input edge for each policy change. The original
eager event sequence and all deferred-delivery cases pass the strict graphical
checker on both renderers, including the four previously failing settled-denial
rows. The SOURCE-022 partial-denial trace which passed had zero retained leaf
payloads and remained nonterminal before re-enable; it never exercised terminal
reopening. The reason that older run lacked a terminal leaf presentation is not
established by this comparison.

The [dated qualification](obol_20260919_readiness_history.md#september-8-live-policy-publication)
owns reports, binaries and broader regression results. This closes the ordinary
live-policy publication boundary. Prototype/configuration writers, general
exception/reentrant-hook behavior, resource magnitudes and release gates remain
in the existing backlog.


### CXX-SOURCE-024 — prototype realization mutates a live drawing

`realizePrototypeWireframe` previously cleared diagnostics, reused the first
line shape, wrote geometry/identity and then updated terminal fields. Immediate
observers saw intermediate records; the preceding library produces 308 callbacks
with an incomplete scene in the new regression. With an auxiliary as the first
line shape, the same lookup could overwrite auxiliary geometry. Bounds also
continued describing the preceding drawing.

The prototype recipe now constructs an explicitly owned private shape. The
existing `BObolPreparedSourcePublication::realize` prepares the replacement and
publishes its geometry, diagnostic identity, source-local bounds, owner metadata,
children/paths and obsolete compact/compiled retirement together. The four-step
revision-dependent square and diagnostic intent remain unchanged. Auxiliary
lines retain their geometry and semantics, and independent nested sources retain
their own records. No additional persistent scene representation or transaction
machinery was introduced. The shared terminal writer replaces duplicate field
assignments, and the shape is retained before child insertion can allocate.

Eight direct/action cases in `test_compact_publication` exercise ordinary,
compact, disabled-field and missing-placement inputs. Every measured allocation
position accepts only the preceding scene or the complete candidate, with
reference/path/handle recovery and retry. Immediate and throwing observers audit
complete recursive scene data. Existing prototype tests exercise stale-input
refresh, picking, snapping, measurement and export; the renderer test exercises
the glyph and color/HUD pixels through OSMesa. The
[dated evidence](obol_20260919_readiness_history.md#september-8-prototype-publication)
owns exact counts and binaries. The subsequent controller effect repair is
recorded below; view successor publication remains separate.

### CXX-REALIZE-001 — observer failure bypasses scene frame publication

The preceding prototype checkpoint reproduced a completed source revision 2
while the scene frame stayed at 1 after a priority-zero observer threw. The
retained old probe exits 1. Source-level atomic publication was insufficient:
action counts, diagnostics and scene revision were updated only after traversal.

The repair attaches a scoped effect to the existing prepared source boundary.
It prepares a failed source's diagnostic before publication, then commits action
counts/diagnostics and advances the scene frame before field/node/path observers.
Prototype, compact wire/mesh and evaluated paths share that boundary. Failure
status publication uses the same commit hook. Repeating an unchanged failure
still reports the attempted failure but creates no new frame revision.

While applying an action, scene getters and summaries read its progress directly.
A scope preserves final counts and moves diagnostics without allocation on every
exit, including a later repository or traversal exception. Each changed source
gets its own frame revision because observers can consume intermediate committed
prefixes. Explicit scene mutation batches retain their documented coalescing.
No worker boundary or additional lifecycle state machine is introduced.

`realization-effects` and `realization-prefix` in `test_compact_publication`
exercise immediate observers, persistent allocation denial, old/complete source
snapshots, later traversal failure, retry, accumulated failure diagnostics,
disabled notifications and explicit batches. The preceding library fails the
new test. The original probe now reports `frame_before=1 frame_after=2` despite
the observer exception. [Dated evidence](obol_20260919_readiness_history.md#september-8-realization-effects)
owns counts and acceptance limits. General reentrant/custom-node/deletion hooks,
repository ownership and large-population resource qualification remain separate
writer/acceptance work. The view-level successor repair is recorded below.

### CXX-REALIZE-002 — observer failure bypasses the view render request

**Repaired September 8.** Previously, `BObolViewController::realizePending`
called `requestLodCapacityRender` only after scene realization returned. With
the scene repair installed, the retained `view_realization_probe.cpp` under
`.build/obol-qualification/20260908-realization-effects` clears the standing
render request, attaches an immediate throwing source observer and invokes the
view API. It prints `threw=1 source_committed=1 frame_before=1 frame_after=2
render_requested=0` and exits 1. That historical probe is retained as failure
evidence; the current binary reports `render_requested=1` and exits 0.

The view now contributes a borrowed prepared effect to the existing synchronous
source publication. Request reasons and policy classification are prepared
before source writes. The commit publishes geometry, action progress, scene
revision and the existing standing request before observers. Request merging
retains stronger current work. Notification attempts both source observers and
the host, preserving the first exception. Completed prefixes retain their
requests, and a request already claimed by the host is not recreated at return.
An unchanged traversal retains the explicit API's ordinary repaint behavior.
No new scheduler, persistent request owner or worker boundary is introduced.

Ordinary render requests share the same preparation/commit/notification helpers.
Previously their kind could change before an allocating reason assignment failed.
The existing trace scope also allowed allocation failure in its destructor to
terminate a completed operation or replace an observer exception. Trace cleanup
now records an explicit dropped transition without throwing; incomplete traces
remain rejected by qualification. Preparation failures still propagate before
the operation starts.

The production `view-realization`, `view-realization-prefix`,
`render-request-publication` and `render-request-trace-unwind` checks cover
immediate field/node and host observers, first-exception retention, persistent
allocation denial, complete prefixes, retry, existing requests, host consumption,
LoD off/auto and tracing off/on. The preceding library fails the new view and
ordinary request checks and aborts the trace-unwind check. The
[dated evidence](obol_20260919_readiness_history.md#september-8-realization-render-requests)
owns exact counts, graphical results and provenance. Resume configuration and
invalidation in the writer inventory, then the wider source/service audit;
large-population resources, general custom hooks and required native hosts remain
separate acceptance work.

### CXX-SOURCE-025 — invalidation publishes before resource and owner retirement

`markStale` wrote `stale` before its reason, terminal record, exact-bound
retirement, source-mesh contract and owned shapes. The retained probe observes a
partial record even without throwing. A throwing observer leaves `stale=1`,
`reason=0`, `status=REALIZED`, exact bounds and an owner still marked current.

The source now shares the existing prepared realization fields and owned-shape
traversal with the status setter. Its fixed source notification record additionally
covers exact-bound retirement. Preparation failure preserves the old state and
staged data; commit publishes every record and retires the applicable producer
lease and compiled/batch currency before callbacks. All notifications are
attempted and the first exception survives. Shared legacy children are visited
once, and traversal stops at independently owned nested sources.

Source/input/database reasons revoke the source-mesh contract and exactness;
view, draw and tessellation reasons preserve usable immutable source data. A
view-only invalidation retains an external producer's terminal failure and
explanation. Compact geometry, population identity and display records remain
unchanged. The unrelated compact-display rebuild is removed from `markStale`;
configuration callbacks and batched configuration retain their own metadata/display
synchronization. This does not qualify those broader writers under exceptions.

The permanent `source-invalidation` check covers seven reason combinations,
current/internal-failed/external-failed initial states, enabled/disabled
notifications, persistent allocation denial, old/complete source and legacy
records, staged-data lifetime, retained compact geometry, shared children,
nested ownership, retry, zero-reason calls and throwing/removing observers. The
preceding library fails the regression; the original probe now observes a complete
record despite its exception. The
[dated evidence](obol_20260919_readiness_history.md#september-8-source-invalidation)
owns exact counts, graphical results and runtime provenance.

### CXX-CONFIG-001 — interrupted configuration leaves mixed fields and lost sensors

**Source configuration boundary repaired; enclosing scene effects remain open.** The original
`configureDatabaseSourceInstanceRepresentation` deleted all twelve field sensors,
then wrote database/identity/path/representation/revision fields individually.
A throwing `instanceKey` observer left a new key with the old path, and subsequent
realization did not recover input-revision detection. The original probe is retained
under `.build/obol-qualification/20260908-source-invalidation`.

Batched configuration and its two wrappers now prepare a detached configuration
candidate, complete owned legacy metadata and changed movable scalar values before
publication. Shared children are visited once; nested sources retain their ownership.
The commit advances the existing database-binding identity, publishes source/shape
fields and applicable exact-bound/source-data retirement, updates compact draw
channels and invalidates compiled presentation before observers. Compact geometry,
occurrence identities and sparse display overrides remain retained. Unchanged
configuration returns without allocation or notification. Aliased arguments are
copied before live writes; field sensors stay attached throughout.

Obol's field touch overload can omit a sensor whose effects were already committed.
That omission travels with this one notification through `SoNotList`; it is not an
ambient suppression flag or sensor detachment. Independent/reentrant field changes
therefore notify normally. Ordinary touch behavior and first-exception retention
remain unchanged. General legacy owner propagation now stops at nested source
boundaries as well, preserving those owners when a retained sensor subsequently fires.

The permanent `source-configuration` regression sweeps allocation failure across
identity and draw-channel changes, realized/failed initial states and enabled/disabled
notifications. It checks complete source/owned metadata, staged-data retirement,
unchanged compact geometry/overrides, worker rejection, retry and no-op behavior.
`committed-field-notification` checks one-notification omission, reentrant edits and
exception recovery. `configuration-notifications` exercises all twelve sensors after
configuration; `configuration-identity` covers database A/B/A, aliased strings,
identity defaults and auxiliary ownership. The preceding physical libBObol fails the
new regression. The original configuration probe now exits 0 with `threw=1
partial_configuration=0 later_input_sensor_recovered=1`.

`setDrawModeState`, `setRepresentationState` and
`retargetDatabaseSourceInstance` now use the same private `ConfigurationPublication`
as batched configuration. Decisions remain in their respective setters: compatible
external draw changes preserve realization and refresh its current identity;
representation changes preserve empty-key semantics; instance-only retargeting with
an unchanged explicit revision retains source data; revision zero retains successor
semantics. Retargeting preserves legacy/material descendant paths with component-boundary
matching, and representation-only refresh preserves those paths exactly. Materials'
multi-value properties and nested sources are untouched. The old live retarget walk
has been removed.

The `configuration-setters` regression covers six changes across realized/failed
and normal/quiet states; `external-draw-configuration` adds four compatible external
draw cases. All 28 sweep persistent allocation rejection and verify preserved or
complete commits, throwing observers, sensor recovery, source stamps, staged-data
retirement and retained geometry. `configuration-path-boundaries` covers slash/prefix
boundaries, aliased inputs, null keys, auxiliary ownership and invalid input. All three
regression groups fail against the preceding physical library and pass with the
repair. Current configuration, setter and path probes exit 0; every setter reports
`threw=1 invalidated_at_exception=1 later_sensor_recovered=1`.

Direct field callbacks and connected-field/engine behavior are qualified separately
by CXX-CONFIG-008. Wider enclosing composition retains its acceptance obligations.
The [dated evidence](obol_20260919_readiness_history.md#september-8-source-configuration-setters)
owns qualification counts and runtime provenance. S1--S6 remain incomplete.

### CXX-CONFIG-008 — direct source fields notify before dependent invalidation

**Repaired and qualified at this boundary.** A direct write to any of the twelve
watched source fields first notified the source container and its observers. The
old field-sensor callback ran afterward, so an observer could see the new source
identity or policy with a current terminal record, old owned-shape metadata and
live resource authority. The retained baseline reproduces that partial state on
a direct `path` assignment.

`SoBRLDatabaseSource::notify` now classifies the triggering field and completes
the same prepared configuration publication before inherited field, node or
delayed observers run. The boundary covers path, instance and representation
identity; representation and draw mode; three tessellation tolerances; BoT LoD
threshold; and source, input and view revisions. A one-shot thread-local marker
identifies notifications emitted by an already committed prepared publication;
the source consumes it before inherited propagation, while a callback's reentrant
raw change remains visible as a new event.

The triggering public field has already committed when preparation begins. The
source therefore establishes an allocation-free fail-safe invalidation first:
the result becomes stale, incompatible source resources and exactness retire,
and compiler/batch currency advances without notifying observers. Preparation
failure cannot leave changed identity authorizing an old result. Retrying the
field event completes the source and owned-shape records. Compatible external
wire/mesh draw changes share one predicate with `setDrawModeState` and preserve
their current geometry.

Obol skips a field's notification path when its container is quiet and the field
has no direct auditors. Twelve private no-op auditors keep that path live; they
do not perform invalidation or own publication state. Disabling an individual
field retains ordinary Obol semantics: the caller must re-enable and touch it to
publish the accumulated change.

The permanent `direct-field-publication` row covers all twelve fields with
immediate and delayed field/node observers, quiet source and field edges, field
and engine connections, callback reentry, observer exception recovery and 87
allocation positions (one preserved, 85 explicit fail-safe, one committed).
The complete state suite passes 80 checks and 805 scenario/sweep rows in 172.74
seconds; the compatible external-draw matrix, GED draw-sync CTest and static-link
probe pass as well. The saved predecessor fails before the source callback can
observe complete invalidation. [Dated evidence](obol_20260919_readiness_history.md#september-10-direct-source-field-publication)
retains logs, exact binaries and hashes.

### CXX-SOURCE-026 — throwing delete callbacks interrupt graph retirement

**Repaired and qualified at the native sensor boundary.** Delete callbacks are
invoked while `SoBase` or an owned field is already being destroyed. A callback
exception previously escaped node and path release before `delete this`, leaving
a reference-count-zero object requiring manual recovery. The same exception from
a field destructor or parent-owned child destruction crossed a `noexcept`
destructor and terminated the process.

`SoDataSensor::invokeDeleteCallback` now contains exceptions from the user hook.
Deletion cannot roll back, so all pending sensors detach and the referenced
object reaches its terminal lifetime state. The public delete-callback contract
documents this rule. Ordinary field, node and path change callbacks retain their
existing propagation behavior.

The retained probe fails against preceding Obol in all four modes: node/path
exit 1 after propagation, and field/child exit 42 through its terminate handler.
All four current modes exit 0. The permanent native test checks one invocation,
automatic detachment and actual destructor completion for node, path,
node-owned-field and parent-owned-child cases. Current Obol passes 1,578 of
1,579 unit tests with one established skip, all 63 integration tests and all
three lifecycle tests.

The complete libBObol compact-publication executable passes 126 direct checks,
21 scenarios and 1,021 pass rows in 309.79 seconds. Eight focused lifecycle,
field and custom-traversal rows also pass. [Dated evidence](obol_20260919_readiness_history.md#september-10-throwing-deletion-callbacks)
retains the baseline/current libraries, exact test binaries, probe, logs and
hashes. Custom-node cleanup, legacy individual remove/truncate callers and the
wider traversal cleanup remain open. This native exception boundary does not
change the formal relations or close resource, platform, sanitizer or graphical
qualification.

### CXX-SOURCE-027 — source clear exposes partially retired realized children

**Repaired and qualified at the public source-clear boundary.**
`SoBRLDatabaseSource::clearRealizedGeometry` previously discarded compact
history, removed its compiled assembly and then called legacy `removeChild` for
each remaining primary child. Immediate source and path observers could see a
partial graph: compact state was already absent while only a prefix of primary
children and paths had retired. A throwing observer could interrupt the loop
and leave that mixed state indefinitely.

The operation now computes all removal indices and prepares one
`SoChildList::Removal` before touching the source. Child order and every audited
path commit quietly; compact history and the compiled assembly then retire
before path and source callbacks run. Allocation failure during preparation
preserves the exact preceding graph and internal state. Notification attempts
continue after an observer exception and rethrow the first failure only after
the complete result is visible. The placement child is always retained,
auxiliary line/source children follow `preserveAuxiliary`, terminal realization
fields remain unchanged and the established return behavior is documented.

The saved predecessor fails the new immediate-observer regression. The focused
current row covers both auxiliary policies, exact retained order, one path per
old child, persistent allocation denial and retry, throwing observers and the
already-clear graph path. Its 24 allocation positions yield 21 preserved old
states and three complete commits. Eight neighboring compiled/external,
retained-presentation, child-lifetime, direct-field and custom-traversal rows
pass. The full compact-publication executable passes 127 direct checks, 21
scenarios and 1,022 pass rows in 310.16 seconds. Prototype and GED draw-sync
CTest rows pass; GED exercises this clear during redraw in all six display
modes. [Dated evidence](obol_20260919_readiness_history.md#september-10-source-realized-geometry-clear)
retains exact baseline/current binaries, dependencies, sources, logs and hashes.

This closes the source clear itself. Custom-node rebuilds, live single-child
detach/reparent sites and the other legacy removal callers remain a finite S1/S3
inventory. Resource limits, native/full-sanitizer platforms and graphical
qualification remain open; the dated completion ranges are unchanged.

### CXX-SOURCE-028 — axes and ADC rebuild expose empty or partial overlay graphs

**Repaired and qualified at the first custom-overlay boundary.**
`SoBRLAxes::rebuildGeometry` and `SoBRLADC::rebuildGeometry` previously removed
all live children before allocating and configuring their replacements. An
immediate node or audited-path observer therefore saw an empty or partial
overlay graph, and allocation failure left the old presentation destroyed.

Both rebuilds now construct complete replacement shapes under local reference
ownership and prepare one `SoChildList::Replacement` before touching the live
graph. Commit publishes child order and path truncation without allocation or
notification; callbacks run only after the complete replacement is visible.
The visible-to-hidden path uses the same transaction. Repeating an already
hidden, empty rebuild allocates and notifies nothing. Preparation failure
preserves the exact old graph, while callback failure leaves the complete new
graph installed, attempts all remaining notifications and then rethrows the
first exception.

The first repaired prototype exposed a native prepared-edit lifetime defect.
Obol temporarily referenced the child-list parent, then released it with
ordinary `unref()` even when a valid caller-owned parent began at reference
count zero. That deleted a newly constructed overlay during its own rebuild.
Prepared removal and replacement now remember whether the parent originally
had an owner and use `unrefNoDelete()` only for the zero-reference case. Native
coverage proves committed replacement, committed removal and abandoned
replacement all preserve that caller-owned parent until its explicit final
release.

The saved libBObol predecessor fails the immediate-observer regression. The
saved preceding Obol combined with repaired libBObol independently fails the
prototype axes-segment assertion. Current focused coverage exercises axes and
ADC replacement and clearing, exact child order, one path per old child,
geometry metadata, persistent allocation denial and retry, throwing observers
and the already-hidden no-op. Its 67 allocation positions yield 63 preserved
old states and four complete commits. Seven neighboring publication/lifetime
rows pass. The full compact-publication executable passes 128 direct checks,
21 scenarios and 1,023 pass rows in 308.17 seconds. Prototype and Qt faceplate
CTest rows pass, as do 1,579 of 1,580 native unit tests with one established
skip, all 63 integration tests and all three lifecycle tests.
[Dated evidence](obol_20260919_readiness_history.md#september-10-axes-and-adc-overlay-rebuild-publication)
retains both negative controls, exact binaries, dependencies, sources, logs and
hashes.

This closes axes and ADC only. Grid and the remaining image, HUD/line, edit,
navigation, viewport, cutting-plane and scene-light builders retain distinct
publication obligations. Live detach/reparent callers and the resource,
platform, sanitizer and graphical gates also remain open; the dated completion
ranges are unchanged.

### CXX-SOURCE-029 — grid rebuild separates derived fields from its HUD graph

**Repaired and qualified at the grid rebuild boundary.**
`SoBRLGrid::rebuildGeometry` previously removed its live HUD child, zeroed three
line-count fields, calculated and published four spacing fields, and then built
and attached replacement geometry. Immediate node, field and audited-path
observers could see values and children from different generations. Allocation
failure after removal permanently discarded the preceding grid.

The method now calculates all derived values and line populations locally and
builds each shape plus the enclosing HUD kit under reference ownership. One
`grid_publish_geometry` helper prepares the complete child replacement and a
seven-field notification set before live mutation. It commits the graph and
derived fields quietly, restores the original notification policy, then drains
path, node and field callbacks. Preparation failure preserves the complete old
state. Callback failure leaves the complete result committed, attempts the
remaining notifications and rethrows the first exception. Hidden clearing uses
the same helper; an already hidden grid with an empty graph and zero counts is a
true no-op.

The saved predecessor fails the new immediate-observer regression. Current
coverage exercises visible replacement and clearing, exact old-child retention,
an audited path, all seven derived fields, nested geometry, persistent
allocation denial and retry, throwing observers and the already-hidden path.
Its 156 allocation positions yield 144 preserved old states and 12 complete
commits. Eight neighboring overlay, source-clear and child-lifetime rows pass.
The full compact-publication executable passes 129 direct checks, 21 scenarios
and 1,024 pass rows in 307.59 seconds. Prototype and Qt faceplate CTest rows
pass. [Dated evidence](obol_20260919_readiness_history.md#september-10-grid-rebuild-publication)
retains exact baseline/current binaries, dependencies, sources, logs and hashes.

This closes `rebuildGeometry`, not the enclosing
`bobol_grid_configure_from_view` writer, whose sequential public input-field
updates retain their own composition obligation. Other custom-node builders,
live detach/reparent callers and the resource/platform/graphical gates remain
open; the dated completion ranges are unchanged.

### CXX-SOURCE-030 — image rebuild clears presentation before stream preparation

**Repaired and qualified for image-plane and viewport-image rebuilds.** Both
methods previously nulled cached texture/face pointers and removed live children
before stream loading or replacement construction. Immediate display, field or
path observers saw an empty or partial presentation. Stream or allocation
failure permanently discarded the preceding image.

Both display nodes now preserve their old presentation through payload load,
fit calculation, textured-quad construction and complete child-replacement
preparation. The shared quad builder returns a reference-owned root and retains
every intermediate node while configuring it, so an exception cannot leak a
partial private graph. One shared image publication helper commits the child,
cached texture/face pointers and two realized-revision fields with notification
disabled, restores policy and then drains path, node and field callbacks.
Observer failure leaves the complete result committed and rethrows only after
the remaining callbacks. Hidden clearing preserves the established realized
revisions; repeated empty clearing does no work.

Payload load refreshes the referenced image source before display commit, so a
dependency notification may expose the complete old presentation. The oracle
accepts complete old or complete new display state and rejects all intermediate
combinations. Source-stream field publication remains a separate obligation.

The saved predecessor fails the immediate-observer regression. Current coverage
exercises image plane and viewport image replacement and clearing, exact
old-child/path ownership, cached pointers, both revisions, dependency preflight,
persistent allocation denial and retry, throwing observers and no-op clearing.
Its 289 allocation positions yield 280 preserved old states and nine complete
commits. Eight neighboring custom-node and child-lifetime rows pass. The full
compact-publication executable passes 130 direct checks, 21 scenarios and 1,025
pass rows in 308.81 seconds. Dedicated image-display, prototype and Qt faceplate
CTest rows pass. [Dated evidence](obol_20260919_readiness_history.md#september-10-image-display-rebuild-publication)
retains exact binaries, dependencies, sources, loader trace, logs and hashes.

This closes the two display rebuild methods and their shared quad builder.
Image-source mutation, edit, navigation, cutting-plane and scene-light builders,
the enclosing grid writer and live detach/reparent sites remain open.
Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-031 — HUD and line-layer rebuilds expose partial child graphs

**Repaired and qualified for HUD-label and line-layer rebuilds.** The HUD-label
method previously cleared its cached label pointer and live children before
allocating its replacement kit and label. The line-layer method removed every
live child, attached a new depth node and then constructed and attached shapes
one at a time. Immediate node or path observers could therefore see an empty or
partially rebuilt graph. Allocation failure discarded the complete old overlay.

The HUD rebuild now constructs a reference-owned kit and label offline, prepares
the complete child replacement, commits the graph and cached pointer, and only
then notifies observers. Its hidden/empty path uses the same prepared replacement
and recognizes an already empty graph with a null pointer as a true no-op. The
line-layer rebuild reserves its owner and child arrays, constructs the depth node
and all usable shapes under reference ownership, and prepares one complete child
replacement. No live child changes before preparation succeeds. Hidden or absent
input clears through the same operation, with an already empty no-op.

The saved predecessor fails the immediate-observer regression. Current coverage
exercises both overlay types under visible replacement and hidden clearing,
exact old-child retention, audited paths, the cached HUD-label pointer, line
shape/point contents, persistent allocation denial and retry, throwing observers
and repeated clearing. Its 142 allocation positions yield 138 preserved old
states and four complete commits. Eight neighboring custom-node and child-
lifetime rows pass. The full compact-publication executable passes 131 direct
checks, 21 scenarios and 1,026 pass rows in 308.60 seconds. Dedicated
image-display, prototype and Qt faceplate CTest rows pass.
[Dated evidence](obol_20260919_readiness_history.md#september-10-hud-and-line-layer-rebuild-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes the HUD-label and line-layer rebuild methods. Image-source mutation,
edit preview/manipulator, navigation, cutting-plane and scene-light builders,
the enclosing grid writer and live detach/reparent sites remain open. Resource,
platform and graphical gates are unchanged.

### CXX-SOURCE-032 — edit manipulator rebuilds expose partial recursive graphs

**Repaired and qualified for the two edit-manipulator rebuild methods.** The
axis and indexed manipulators previously removed their complete live graphs,
attached a light model and then constructed and attached each replacement
subtree one at a time. Immediate node or path observers could see an empty or
partially rebuilt manipulator. Allocation failure discarded the preceding
presentation and could leak an unowned partially constructed subtree.

Both methods now allocate and configure every top-level and nested scene node
under reference ownership outside the live child list. One shared helper
prepares the complete child order, commits it without allocation and notifies
audited paths and the node only after every recursive subtree is complete.
Hidden or empty manipulators clear through the same prepared replacement; an
already empty graph is a true no-op. Precommit failure preserves the exact old
children and paths, while callback failure leaves the complete new graph
installed and drains remaining notifications before rethrowing.

The saved predecessor fails the immediate-observer regression. Current focused
coverage exercises axis and indexed manipulators under visible replacement and
hidden clearing, every old top-level child and path, complete axis/edge/point
and selected/hover/active emphasis subtrees, persistent allocation denial and
retry, throwing observers and repeated clearing. In an isolated process its
449 allocation positions yield 444 preserved old states and five complete
commits. When preceded by the other compact checks, one process-global setup
allocation is absent and the row passes 448 positions with 444 preserved and
four committed. Six neighboring graph/publication rows pass. The full compact-
publication executable passes 132 direct checks, 21 scenarios and 1,027 pass
rows in 307.53 seconds. Dedicated edit-manipulator, image-display, prototype
and Qt faceplate CTest rows pass.
[Dated evidence](obol_20260919_readiness_history.md#september-10-edit-manipulator-rebuild-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes `SoBRLEditManipulator::rebuildGeometry` and
`SoBRLIndexedEditManipulator::rebuildGeometry`. Their public-field and
`setEllipsoidAxes`/`setTopology` composition boundaries remain open, along with
edit-preview, image-source, navigation, cutting-plane, scene-light and live
detach/reparent work. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-033 — edit preview publishes children and status in separate steps

**Repaired and qualified for direct, transformed, failed and cleared preview
output.** `setLineSet`, `setTransformedLineSet` and `clearPreview` previously
removed the live child before replacement allocation, then published realized
revisions, stale state and preview status sequentially. Immediate node, field or
path observers could see mixed generations. Allocation failure lost the old
preview, and transformed construction could leak unowned intermediate nodes.

One shared publication helper now prepares the complete child replacement and
notifications for both realized revisions, status and stale state before live
mutation. Commit installs the graph and four field values with notifications
disabled. It restores policy, then drains path, node and changed-field callbacks,
retaining the first observer exception. Every direct/transformed scene node is
reference-owned during construction. Invalid input atomically publishes an
empty FAILED/stale result while preserving the last realized revision pair;
explicit clear advances both realized revisions and publishes EMPTY/non-stale.
Repeated invalid and clear outcomes are true no-ops.

The saved predecessor fails the immediate-observer regression. Current focused
coverage exercises all four outcomes, exact old-child/path retention, shape and
recursive transform content, all four fields, persistent allocation denial and
retry, throwing observers and repeated no-op behavior. In an isolated process,
its 96 allocation positions yield 76 preserved old states and 20 complete
commits. When preceded by other compact checks, three process-global setup
allocations are absent and it passes 93 positions with 76 preserved and 17
committed. Six neighboring edit/field/child/path rows pass. The full compact-
publication executable passes 133 direct checks, 21 scenarios and 1,028 pass
rows in 306.86 seconds. Prototype, view-store, compact-edit, edit-manipulator
and Qt faceplate CTest rows pass.
[Dated evidence](obol_20260919_readiness_history.md#september-10-edit-preview-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes the three preview output methods. `setEditIntent`,
`markSourceRevision` and `markInputsRevision` remain an input-field composition
boundary. Manipulator setters, image-source mutation, navigation, cutting-plane,
scene-light and live detach/reparent work also remain open. Resource, platform
and graphical gates are unchanged.

### CXX-SOURCE-034 — navigation rebuild clears cached HUD state before construction

**Repaired and qualified for direct navigation-gizmo rebuilds.** The preceding
method cleared both cached node pointers and removed the live HUD before it
allocated any replacement node. It then constructed and attached the HUD,
widget and every recursive face, ring, axis, endpoint and label subtree through
unowned raw pointers. Immediate observers could see an empty graph, and any
construction failure discarded the complete preceding gizmo and could leak the
partially constructed private graph.

`rebuildGeometry` now builds the entire HUD and every nested node under reference
ownership outside the live child list. It prepares one child replacement,
commits the complete graph and both cached pointers, then notifies audited paths
and the gizmo. Hidden clearing uses the same prepared replacement and publishes
null cached pointers before callbacks; an already empty hidden gizmo is a true
no-op. Precommit failure preserves the exact preceding HUD, cached HUD pointer
and path. Callback failure leaves the complete new graph installed and drains
the remaining notifications before rethrowing.

The current test executable forced to load the saved predecessor rejects the
immediate-observer regression. Isolated focused coverage checks visible cube-to-
circles replacement and hidden clearing, the complete 20-child widget graph,
cached HUD identity, retained old HUD and path, all allocation failures, retry,
throwing observers and repeated hidden no-op behavior. It passes 949 allocation
positions: 945 preserve the old state and four reach a complete commit. After
preceding compact checks initialize process-global state, it passes 948
positions with 945 preserved and three committed. Five neighboring graph rows
and the dedicated navigation CTest pass. The full compact-publication executable
passes 134 direct checks, 21 scenarios and 1,029 pass rows in 304.57 seconds.
[Dated evidence](obol_20260919_readiness_history.md#september-10-navigation-gizmo-rebuild-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes only direct `rebuildGeometry` publication. `setCamera`/camera
synchronization and hover/active setters remain enclosing input-plus-rebuild
composition boundaries. Image-source mutation, grid configuration, preview and
manipulator input setters, cutting-plane, scene-light and live detach/reparent
work remain open. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-035 — cutting affordance and scene lights rebuild incrementally

**Repaired and qualified for the two controller-owned output builders.** The
cutting-plane affordance previously removed its current HUD before validating
the plane, camera, viewport and dimensions or allocating projected geometry.
The scene-light rebuild removed every current light and appended replacements
one by one. Its lookup helper also detached and reinserted the light group on
every call, while camera-rig placement did the same for three lights.
Immediate node or path observers could therefore see empty and partial graphs;
allocation failure discarded the complete preceding presentation and could
leak privately constructed nodes.

The affordance now computes and projects its complete line set before building
the reference-owned HUD and six-node widget offline. One prepared replacement
publishes an existing affordance, while a first affordance attaches only after
its complete HUD is present. Invalid or disabled output clears through the same
prepared operation, and an empty disabled affordance is a true no-op. Scene
lights reserve their publication storage, build every typed light under
reference ownership and replace the group once. Empty scene-light output also
has an allocation-free no-op. Camera-rig and scene-light group ordering now use
complete prepared root orders and skip unchanged orders, so paths never observe
temporary detachment.

The current focused executable forced to load the saved predecessor rejects
both immediate-observer rows. Cutting coverage checks replacement and clearing,
the complete HUD widget, 36 projected points, 18 line segments, retained HUD
and path, all allocation failures, retry, throwing observers and disabled no-op
behavior. It passes 147 positions with 145 preserved and two committed. Scene-
light coverage checks directional, spot and point replacements, fields and
enabled state, all old children/paths, clearing, retry, throwing observers and
empty no-op behavior. It passes 80 positions with 69 preserved and 11 committed.
Five neighboring graph rows and the window-host and retained-raytrace CTests
pass. The full compact-publication executable passes 136 direct checks, 21
scenarios and 1,031 pass rows in 305.13 seconds.
[Dated evidence](obol_20260919_readiness_history.md#september-10-cutting-affordance-and-scene-light-rebuild-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes cutting-affordance and scene-light output rebuilding plus their
stable root ordering. `setCuttingPlane`, `setCuttingPlaneEnabled`,
`setSceneLights` and `setSceneLightsEnabled` retain separate input/state/output
composition work. Image-source mutation, grid configuration, preview,
manipulator and navigation input setters, render-environment repair and live
detach/reparent work remain open. Resource, platform and graphical gates are
unchanged.

### CXX-SOURCE-036 — render-environment repair publishes children and root order separately

**Repaired and qualified for render-environment creation/repair and its root
order.** An existing environment repaired a missing depth buffer, ambient
environment, light model and clipping planes through separate insertions. It
then moved the live environment to root index zero through detach/reinsert, and
camera-rig repair published its missing or reordered lights afterward. A root
or environment observer could therefore see a partially repaired environment,
a temporarily absent environment, or a complete environment with an incomplete
root order. Allocation failure could retain any such prefix. When the
environment was absent but a named headlight already existed, the separately
allocated migration headlight also had no owner or parent.

The controller now discovers all retained nodes first and owns every missing
node until publication. It preserves unrelated environment and root extension
children, prepares both live child-list replacements before changing either,
commits both complete orders before notification, and attempts both notification
sets while retaining the first callback exception. A newly created environment
is assembled completely offline before one root publication. The root has one
environment followed by camera/key/fill/rim, or environment/key/fill/rim while
no camera is attached. Repeated pointer occurrences are normalized, and an
already complete graph allocates and notifies nothing.

The current focused executable forced to load the saved CXX-SOURCE-035 library
rejects the immediate-observer row. Current coverage checks absent creation and
in-place repair, all required node types and defaults, unrelated extension
children, every retained root/environment path, every allocation failure,
retry, throwing observers and the complete no-op. It passes 147 allocation
positions with 131 exact preceding states and 16 complete commits. Six
neighboring publication rows plus OSMesa-render, window-host and retained-
raytrace CTests pass. In the full process the same row passes 145 positions with
131 preserved and 14 committed. The exact full executable passes 137 direct
checks, 21 scenarios and 1,032 pass rows in 305.88 seconds. Current compact/
shared/static hashes are `b97e2f64...`, `00375975...` and `db31baa5...`; Obol
and OSMesa remain `21872f7f...` and `097aac27...`.
[Dated evidence](obol_20260919_readiness_history.md#september-10-render-environment-repair-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes the lighting-side live detach/reparent inventory. Framebuffer host
and endpoint migration, feature-store attachment/reordering and Qt faceplate
adaptation remain finite live-child boundaries. Lighting, cutting-plane and
scene-light input setters remain separate composition work, together with image
source, grid, preview, manipulator and navigation input writers. Resource,
platform and graphical gates are unchanged.

### CXX-SOURCE-037 — feature-overlay reordering publishes one complete root order

**Repaired and qualified for pure feature-overlay root reordering.**
`BObolFeatureStore::Impl::reorderOverlayNodes` previously detached every
feature overlay and reattached the sorted records one at a time. Immediate root
and path observers could see truncated prefixes, and allocation or callback
failure could retain one of them. The procedure also relied implicitly on the
attachment helper to avoid duplicating a custom node shared by multiple feature
records.

The store now sorts records, reduces them to the established unique node
occurrences, and constructs the complete successor root offline. Unrelated
endpoint-owned children retain relative order, sorted overlays follow them,
and all occurrences owned by the reorder are normalized. One prepared child-
list replacement commits the final order before notification. Every retained
path remaps directly from its preceding index to its final index. An exact
repeated overlay-setting call remains an allocation-free, notification-free
no-op.

The current focused executable forced to load the saved CXX-SOURCE-036 library
rejects the immediate-observer row. Current coverage uses two typed HUD line
features, two records which share one custom node, and two unrelated root
extensions. It checks final order and identities, every retained path, every
allocation failure, retry, a throwing observer and the exact no-op. The isolated
selector passes 129 allocation positions with 125 exact preceding states and
four complete commits. Six neighboring publication rows plus the feature-store
and Qt faceplate CTests pass. In the full process the row passes 128 positions
with 125 preserved and three committed. The exact full executable passes 138
direct checks, 21 scenarios and 1,033 pass rows in 305.40 seconds. Current
compact/shared/static hashes are `9d403eee...`, `f6a048a4...` and `61e51c22...`;
Obol and OSMesa remain `21872f7f...` and `097aac27...`.
[Dated evidence](obol_20260919_readiness_history.md#september-10-feature-overlay-order-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes only the root reorder step. Feature node replacement/reparenting,
multi-record `setController` migration, and the metadata/node/order composition
in `setOverlayInfo` remain separate view-store boundaries. Framebuffer host,
display-endpoint and Qt faceplate live-child work also remains open. Resource,
platform and graphical gates are unchanged.

### CXX-SOURCE-038 — framebuffer composition publishes viewport and layer roots together

**Repaired and qualified for window-host framebuffer composition changes.**
`BObolWindowHost::setFramebufferComposition` previously removed the retained
viewport from any layer root, wrote its visibility/layer fields, rebuilt its
realized children and then attached it to the destination. Root, field, source
and path observers could see those stages. Allocation or rebuild failure after
removal retained a detached viewport, and a path through its realized image
could retire before the intended destination existed.

The host now prepares one private viewport candidate from the accepted current
stream state. This read does not refresh or notify the shared live image source.
A scalar publication owns the viewport fields; a child replacement owns the
realized graph and cached render-node pointers; three prepared root replacements
own removal, normalization and destination insertion. Every participant
prepares before any live change, all commit before notification, unrelated
layer children keep their order, and notification exceptions do not suppress
the remaining notifications or render request. An exact repeated composition
request allocates and notifies nothing. `SoSFVec2f` joins the scalar publication
types because viewport position and size are part of the committed value.

The current focused executable forced to load the saved CXX-SOURCE-037 library
rejects the immediate-observer row. Current coverage moves a realized viewport
from overlay to underlay among unrelated children in all three roots. It checks
viewport/image-source identity, scalar fields, realized children and cached
texture/face pointers, direct paths and a path through the retired realized
child, every allocation failure, retry, a throwing observer and the exact no-op.
The isolated selector passes 257 allocation positions with 254 exact preceding
states and three complete commits. Six neighboring publication rows plus image-
source, image-display, window-host, headless-window-host and Qt window-host
CTests pass. In the full process the row passes 256 positions with 254 preserved
and two committed. The exact full executable passes 139 direct checks, 21
scenarios and 1,034 pass rows in 305.23 seconds. Current compact/shared/static
hashes are `a7aa017e...`, `f70823e3...` and `0bd218f9...`; Obol and OSMesa remain
`21872f7f...` and `097aac27...`.
[Dated evidence](obol_20260919_readiness_history.md#september-10-framebuffer-composition-publication)
retains exact baseline/current binaries, dependencies, sources, loader traces,
logs and hashes.

This closes composition changes for an already open window-host framebuffer.
Retained-raytrace display-endpoint migration is closed separately by
CXX-SOURCE-039 below. Framebuffer open/close cleanup, viewport and cursor input
setters, feature-store reparenting and Qt faceplate composition remain separate
boundaries. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-039 — raytrace layer migration publishes one retained composition

**Repaired and qualified for retained-raytrace display-endpoint layer changes.**
`composition.rt.layer` previously removed the retained viewport from its old
root, wrote the layer field, rebuilt its realized image, and inserted it into
the destination as separate observable operations. The generic RT-property
tail then restarted the renderer, adding a second presentation publication to
a composition-only change. An immediate observer could interrupt or observe a
detached viewport, incomplete root order, or prematurely retired image path.

The endpoint and window host now share one private three-root successor
preparation. Endpoint migration prepares its scalar viewport publication and
all root replacements before changing live state. It commits the layer field,
the three root orders and endpoint presentation bookkeeping before notifying;
the RT viewport remains first in its destination and unrelated children retain
their order. The layer-only path preserves the accepted image source, realized
children and cached texture/face nodes, and it requests presentation without
restarting ray work. Notification failure does not suppress later prepared
notifications. An exact repeated property request allocates and notifies
nothing.

The current focused executable forced to load the saved CXX-SOURCE-038 library
rejects the immediate-observer row because it sees partial cross-root state.
Current coverage migrates a realized RT viewport from interlay to overlay among
unrelated children in every layer. It checks endpoint policy, viewport/source/
realized-node identity, exact root order, direct and nested retained paths,
every allocation failure, retry, a throwing observer and the exact no-op. The
isolated selector passes 80 allocation positions with 78 exact preceding states
and two complete commits. Six neighboring publication rows, retained-raytrace,
window-host, headless-window-host and Qt window-host tests, and the static link
check pass. The exact full executable passes 140 direct checks, 21 scenarios and
1,035 pass rows in 308.97 seconds. Current compact/shared/static hashes are
`a896c5f9...`, `8e5e26b6...` and `e462056f...`; Obol and OSMesa remain
`21872f7f...` and `097aac27...`. Exact binaries, dependencies, sources, loader
traces, logs and hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-raytrace-layer-composition-publication).

This closes retained-raytrace endpoint layer migration. Endpoint presentation
cleanup, framebuffer open cleanup, viewport and cursor input setters,
feature-store reparenting and Qt faceplate composition remain separate
boundaries. Explicit framebuffer close is closed separately by CXX-SOURCE-040
below. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-040 — framebuffer close publishes attachment and root retirement together

**Repaired and qualified for explicit window-host framebuffer close.**
`BObolWindowHost::closeFramebuffer` previously removed the retained viewport
from its layer roots, released viewport/source ownership, and only then erased
the host attachment. Root and path observers could therefore see no viewport
while the host getters still returned the retiring objects. A throwing callback
could interrupt the sequence and retain that mixed state for the next call.

Close now prepares the complete three-root successor lists before any live
change. After preparation succeeds, it erases the attachment and commits every
root and retained-path update before notification. Viewport/source ownership is
released only after prepared graph notification, and the presentation request
and every cleanup step preserve the first callback failure while continuing the
remaining work. Allocation failure during preparation leaves the exact live
attachment unchanged. Repeating close after retirement allocates and notifies
nothing.

The current focused executable forced to load the saved CXX-SOURCE-039 library
rejects the immediate-observer row. Current coverage closes framebuffers from
off, underlay, interlay and overlay modes among unrelated children in all three
roots. It checks attachment getters/count, root order, direct and nested paths,
retained viewport/source lifetime, every allocation failure, retry, throwing
observers and the exact no-op. The isolated selector passes 54 allocation
positions with 50 exact preceding states and four complete commits. Six
neighboring publication rows, image-source, image-display, window-host,
headless-window-host and Qt window-host tests, and the static link check pass.
The exact full executable passes 141 direct checks, 21 scenarios and 1,036 pass
rows in 308.76 seconds. Current compact/shared/static hashes are `d532d0c4...`,
`d5a6d699...` and `168f5d8e...`; Obol and OSMesa remain `21872f7f...` and
`097aac27...`. Exact binaries, dependencies, sources, loader traces, logs and
hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-framebuffer-close-publication).

This closes explicit framebuffer close. Explicit framebuffer construction/open
is closed separately by CXX-SOURCE-041 below. Window-host destructor-only
cleanup, endpoint presentation cleanup, viewport and cursor input setters,
feature-store reparenting and Qt faceplate composition remain separate
boundaries. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-041 — framebuffer open publishes attachment and root insertion together

**Repaired and qualified for explicit window-host framebuffer construction and
open.** `BObolWindowHost::openFramebuffer` previously reference-owned the new
image source and viewport through raw pointers, inserted the host attachment,
and then called ordinary root insertion. Allocation failure before attachment
could leak those private nodes. Allocation failure in root insertion could
leave the host getters reporting a framebuffer whose viewport was absent from
every composition root.

Source and viewport construction now use scoped node references. The host
builds their complete realized image while private, prepares all three root
successors, and inserts the attachment only after every fallible graph
preparation succeeds. Attachment insertion and root commit then run without an
observable callback between them. The scoped references transfer their exact
host ownership before prepared notifications begin. Notification or
presentation-request failure can report failure only after the complete open
state exists; precommit failure releases private nodes and preserves the exact
absent state. Repeating open for the same framebuffer allocates and notifies
nothing.

The current focused executable forced to load the saved CXX-SOURCE-040 library
rejects an allocation row which retains the attachment without root insertion.
Current coverage opens a memory framebuffer among unrelated children in all
three roots. It checks attachment count/lookups, image source and viewport
realization, destination order, retained extension paths, every allocation
failure, retry, throwing observers and the exact no-op. The isolated selector
passes 458 allocation positions with 266 exact absent states and 192 complete
commits. Six neighboring publication rows, image-source, image-display, window-
host, headless-window-host and Qt window-host tests, and the static link check
pass. The exact full executable passes 142 direct checks, 21 scenarios and
1,037 pass rows in 308.04 seconds. Current compact/shared/static hashes are
`0269e48f...`, `e6ba539e...` and `9f6dc606...`; Obol and OSMesa remain
`21872f7f...` and `097aac27...`. Exact binaries, dependencies, sources, loader
traces, logs and hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-framebuffer-open-publication).

This closes source/viewport construction and attachment/root insertion for
explicit framebuffer open. Neutral-host descriptor/open rollback is closed
separately by CXX-SOURCE-042 below. Viewport and platform-host policy,
constructor/destructor-only cleanup, endpoint presentation cleanup, viewport/
cursor input setters, feature-store reparenting and Qt faceplate composition
remain separate at this checkpoint. Framebuffer-host destruction is closed by
CXX-SOURCE-043 below. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-042 — framebuffer preparation precedes neutral-host policy commit

**Repaired and qualified for neutral-host descriptor/open state with an
unchanged viewport size.** `BObolWindowHost::openFramebuffer` previously called
the virtual host open before validating the stream and before constructing the
image source, viewport, attachment capacity and three-root successor. Any later
failure could therefore leave `isOpen()` and `getDesc()` at the requested
framebuffer policy while the framebuffer attachment remained absent.

The base host now prepares and sanitizes a copied descriptor privately and
commits it by pointer swap only after its prerequisites succeed. Framebuffer
open validates and constructs the stream-backed image graph, reserves the next
attachment slot and prepares all root replacements before invoking host open.
After host open succeeds, attachment insertion, root commit and ownership
transfer cannot allocate; notification and presentation failure report only
after the complete next state exists. The shared private root prerequisite also
lets the headless host create its root and camera before the base descriptor
commit.

The current focused executable forced to load the saved CXX-SOURCE-041 library
rejects a failure row which retains the new host policy without an attachment.
Current coverage starts from both closed and already-open hosts, using equal
preceding and requested viewport dimensions to isolate descriptor/open policy
from the still-separate viewport setter. It checks every descriptor field,
`isOpen()`, controller viewport dimensions, render-request state, attachment and
three-root/path state, every allocation failure, retry, throwing observers and
the exact no-op. The isolated selector passes 603 allocation positions with 439
exact preceding policies and 164 complete commits. Six neighboring publication
rows, image-source, image-display, window-host, headless-window-host and Qt
window-host tests, and the static link check pass. The exact full executable
passes 143 direct checks, 21 scenarios and 1,038 pass rows in 308.71 seconds.
Current compact/shared/static hashes are `ed914796...`, `2eba79a8...` and
`0b00def1...`; Obol and OSMesa remain `21872f7f...` and `097aac27...`. Exact
binaries, dependencies, sources, loader traces, logs and hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-framebuffer-open-policy-publication).

This closes the neutral base-host descriptor/open portion of framebuffer open.
At this checkpoint viewport-size publication, headless context-manager
lifecycle, Qt-owned widget policy, endpoint presentation cleanup, viewport/
cursor input setters, feature-store reparenting and Qt faceplate composition
remain separate boundaries. Framebuffer-host destruction is closed separately
by CXX-SOURCE-043 below, and full-window viewport-size publication by
CXX-SOURCE-044. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-043 — framebuffer host destruction does not run live publication

**Repaired and qualified for base-host framebuffer cleanup with owned and
borrowed controllers.** Base, headless and Qt host destructors previously
invoked their normal live `close()` paths. The owned view-controller destructor
also called live provider removal, camera repair and scene-root publication.
Those operations allocate and notify observers even though destruction has no
successor public state. Because the destructors are implicitly `noexcept`, a
denied allocation or throwing observer reached `std::terminate`.

The base host now contains stream/display detachment failures and invokes a
private no-throw framebuffer cleanup. Each attachment first uses the qualified
close publication. A callback failure occurs after commit and cleanup
continues. If root-successor preparation cannot allocate, the host releases its
manual viewport/source references and leaves the complete graph owned by its
existing roots. Controller ownership then determines final retirement: an
owned controller removes those roots during the same destruction; a borrowed
controller remains valid and keeps the complete image graph. Headless and Qt
destructors invoke this cleanup before their platform resources are destroyed
and contain failures from their platform close operations.

View-controller destruction now releases progressive-provider callbacks,
camera ownership, view-attachment state and retained render composition
directly. It no longer calculates convergence, repairs lighting, constructs
repository membership or requests a frame merely to tear the same objects
down. The embedded scene controller retains its existing allocation-free root
and repository retirement.

The current focused executable forced to load the saved CXX-SOURCE-042 library
exits 42 on the first denied allocation. Current coverage checks exact retained
source/viewport ownership for both owned and borrowed controllers across 18
allocation positions and throwing overlay-root observers. Eight borrowed-path
preparation failures retain the documented complete root-owned fallback; all
other paths detach the graph. The framebuffer open-policy, open, close,
composition, retained raytrace-layer, notification-recovery and graph-
destruction neighbors pass, as do image-source, image-display, window-host,
headless-window-host, Qt window-host, static-link and all seven conformance
checks. The exact full executable passes 144 direct checks, 21 scenarios and
1,039 pass rows in 310.168 seconds. Current compact/shared/static hashes are
`02904d05...`, `4fd8c93c...` and `321acd20...`; Obol and OSMesa remain
`21872f7f...` and `097aac27...`. Exact binaries, dependencies, sources, loader
traces, logs and hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-framebuffer-host-destruction).

This closes framebuffer cleanup and the owned-controller operations reached by
that host path. Exact borrowed-controller root removal when allocation is
unavailable, arbitrary camera/root observer failure, remaining endpoint
cleanup, general sub-viewport/cursor input setters, feature-store reparenting
and Qt faceplate composition remain separate boundaries. Full-window viewport
size is closed separately by CXX-SOURCE-044 below. Resource, platform and
graphical gates are unchanged.

### CXX-SOURCE-044 — viewport size publishes one complete controller state

**Repaired and qualified for full-window viewport-size publication on ordinary
and automatic-LoD controllers.** `BObolViewController::setViewportSize`
previously changed its stored region and the `SoViewport`, then asked the
render manager to detach and reattach its unchanged scene graph before updating
LoD state. The nested LoD request woke the endpoint with `lod-view`; automatic
mode could also wake it with `progressive-work`. Only afterward did the setter
replace the standing reason with `viewport-size`, without another wake when the
request level was already pending. A render-manager or convergence-snapshot
allocation failure could instead leave the three public viewport copies split
or changed under the old LoD revision.

The setter now prepares the final capacity render request before mutation and
updates only the render manager's viewport region. LoD synchronization can
install its immediate view-revision and progressive-work state without
dispatching an intermediate endpoint wake. If the automatic convergence
snapshot cannot allocate, the controller, `SoViewport`, and render-manager
regions all return to the preceding value before the exception escapes. After
successful synchronization, the prepared `viewport-size` request commits and
wakes the endpoint once. Existing private one-argument/no-argument helper
symbols remain as wrappers around the new internal options.

The current focused executable forced to load the exact CXX-SOURCE-043 library
rejects the first success callback because it observes an intermediate state.
Current coverage checks ordinary and automatic-LoD controllers across six
allocation outcomes: four denied automatic snapshots preserve the exact old
state without notification, while two success paths publish the complete new
region, immediate LoD view revision, correct progressive-work level and exact
reason. Both modes also pass retry, allocation-free exact no-op and a throwing
frame callback after commit. The framebuffer open-policy, open, close,
composition, retained raytrace-layer, notification-recovery and graph-
destruction neighbors pass, as do prototype, image-source, image-display,
LoD-update, window-host, headless-window-host, Qt-controller, Qt-window-host,
static-link and all seven conformance checks. The exact full executable passes
145 direct checks, 21 scenarios and 1,040 pass rows in 367.41 seconds. Current
compact/shared/static hashes are `d9ab4dad...`, `fc6c558d...` and
`a85333a8...`; Obol and OSMesa remain `21872f7f...` and `097aac27...`. Exact
binaries, dependencies, sources, loader traces, logs and hashes are retained in
the [dated evidence](obol_20260919_readiness_history.md#september-10-viewport-size-publication).

This closes `setViewportSize()` and its ordinary/automatic-LoD host wakeup.
General sub-viewport publication is closed separately by CXX-SOURCE-045 below.
Arbitrary Coin field observers, platform-host sizing policy, camera/cursor
setters, feature-store reparenting and Qt faceplate composition remain separate
boundaries. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-045 — both viewport setters share one publication transaction

**Repaired and qualified for arbitrary controller sub-viewport publication.**
`BObolViewController::setViewportRegion` retained the pre-CXX-SOURCE-044
sequence: it changed the controller and `SoViewport`, detached and reattached
the render manager's unchanged scene graph, dispatched nested `lod-view` and
automatic `progressive-work` wakes, and only then requested `viewport`. The
same allocation and partial-callback failures therefore remained for callers
which set an origin, partial extent or pixel density directly.

Both public viewport setters now construct their requested region and enter one
private publication transaction. Exact equality retires as a no-op. A real
change prepares the operation's exact render request, installs the controller,
`SoViewport` and render-manager regions, synchronizes the immediate LoD view
revision without an intermediate endpoint wake, then commits and notifies once.
The full-window setter retains its width/height normalization and
`viewport-size` reason; the general setter retains every supplied region field
and uses `viewport`.

The current focused executable forced to load the exact CXX-SOURCE-044 library
rejects the first success callback because it observes an intermediate state.
Current coverage changes window dimensions, sub-viewport origin and extent,
and pixel density in ordinary and automatic-LoD controllers. Six allocation
outcomes pass with four exact old states and two complete commits, plus retry,
allocation-free exact no-op and throwing frame callbacks. The CXX-SOURCE-044
viewport-size selector and framebuffer open-policy, open, close, composition,
retained raytrace-layer, notification-recovery and graph-destruction neighbors
pass, as do prototype, image-source, image-display, LoD-update, window-host,
headless-window-host, Qt-controller, Qt-window-host, static-link and all seven
conformance checks. The exact full executable passes 146 direct checks, 21
scenarios and 1,041 pass rows in 341.95 seconds. Current compact/shared/static
hashes are `82c61f31...`, `77dad89d...` and `7f0abf24...`; Obol and OSMesa
remain `21872f7f...` and `097aac27...`. Exact binaries, dependencies, sources,
loader traces, logs and hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-viewport-region-publication).

This closes `setViewportRegion()` and consolidates both controller viewport
setters. Established-root neutral host composition is closed separately by
CXX-SOURCE-046 below. Arbitrary Coin field observers, first-root and derived-
host policy, camera/cursor setters, feature-store reparenting and Qt faceplate
composition remain separate boundaries. Resource, platform and graphical
gates are unchanged.

### CXX-SOURCE-046 — established-root neutral host open publishes window policy and controller viewport together

**Repaired and qualified for established-root neutral base-host open
composition.**
`BObolWindowHost::open()` previously called the controller's viewport setter
before swapping its copied descriptor and setting the open flag. The viewport
transaction could therefore wake the frame endpoint with the new region, LoD
view revision, progressive-work state and render request while host observers
still saw the preceding descriptor and open state. The exact repeated-open
path also copied its descriptor before discovering that nothing had changed.

Open now compares the normalized request with an already-open host before any
allocation. For this boundary, group-root validation is a no-op because the
host has already established its scene root. A real change privately prepares
the descriptor and controller viewport publication. Only after every fallible
preparation succeeds does it swap the descriptor and set the open flag, then
finish the prepared controller publication and issue its one endpoint wake.
Controller preparation retains the CXX-SOURCE-045 transaction: automatic
convergence snapshot failure restores all three viewport copies before
escaping. Callback failure propagates after the complete host/controller
policy has committed.

The current focused executable forced to load the exact CXX-SOURCE-045 library
rejects its first success callback because that callback sees the preceding
host policy. Current coverage changes every descriptor field and viewport
dimension for closed and already-open hosts after root establishment, using
ordinary and automatic-LoD controllers. Forty allocation outcomes pass: 36
preserve the complete preceding policy without a callback and four publish the
complete successor policy. All four modes also pass retry, exact allocation-
free repeated open and a throwing callback after commit. The two preceding viewport
selectors and framebuffer open-policy, open, close, composition, retained
raytrace-layer, notification-recovery and graph-destruction neighbors pass, as
do prototype, image-source, image-display, LoD-update, window-host, headless-
window-host, Qt-controller, Qt-window-host, static-link and all seven
conformance checks. The exact full executable passes 147 direct checks, 21
scenarios and 1,042 pass rows in 311.13 seconds. Current compact/shared/static
hashes are `0deaa573...`, `5179915b...` and `3169d391...`; Obol and
OSMesa remain `21872f7f...` and `097aac27...`. Exact binaries,
dependencies, sources, loader traces, logs and hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-window-open-viewport-publication).

This closes established-root neutral base-host descriptor/viewport
publication. Host-owned first-root construction is closed separately by
CXX-SOURCE-047 below. Borrowed-controller root establishment, headless context
and camera setup, Qt widget policy, arbitrary camera/cursor setters,
feature-store reparenting and Qt faceplate composition remain separate
boundaries. Resource, platform and graphical gates are unchanged.

### CXX-SOURCE-047 — host-owned root exists before first open can publish

**Repaired and qualified for construction and first open of the neutral
host-owned controller.** A new `BObolWindowHost` previously exposed a
rootless controller. Its first `open()` called the live scene-root setter,
which installed an empty group and requested `scene-root` before the host
descriptor, open flag and viewport transaction committed. An attached endpoint
could observe that mixed state. Because the render-request level was already
pending, the final `viewport-size` request could replace the reason without
another wake. The controller's stored and render-manager regions also started
at 1x1 while `SoViewport` retained Coin's different default until a
successful resize.

The host private constructor now keeps its new controller in temporary unique
ownership, establishes the empty group while the entire host remains
inaccessible, retires the construction-only render/progressive work, and only
then releases the controller into host ownership. If any construction step
throws, temporary ownership destroys the incomplete controller and no host
object becomes visible. The controller private constructor initializes both
`SoViewport` and the render manager from its stored region, establishing the
three-copy invariant for every controller before public construction returns.
The first visible open then uses the CXX-SOURCE-046 descriptor/viewport
transaction and issues one complete `viewport-size` wake.

The current focused executable forced to load the exact CXX-SOURCE-046 library
rejects the first `scene-root` callback because it sees the default closed
host policy. Current ordinary and automatic-LoD first-open coverage passes 20
allocation outcomes: 18 preserve the complete closed state and two commit the
complete open state. Both modes also pass retry, allocation-free exact repeat
and a throwing callback after commit. A separate constructor sweep covers 220
allocation positions; 205 reject construction without exposing an object and
15 produce a complete closed host with one group root, three matching viewport
copies and no pending host work. The CXX-SOURCE-046 host-open selector and
viewport-size, viewport-region, framebuffer open-policy, open, close,
composition, retained raytrace-layer, notification-recovery and graph-
destruction neighbors pass, as do prototype, image-source, image-display,
LoD-update, window-host, headless-window-host, Qt-controller, Qt-window-host,
static-link and all seven conformance checks. The exact full executable passes
148 direct checks, 21 scenarios and 1,043 pass rows in 349.88 seconds. Current
compact/shared/static hashes are `5dd6fcdd...`, `0d7988fe...` and
`8a27709c...`; Obol and OSMesa remain `21872f7f...` and `097aac27...`.
Exact binaries, dependencies, sources, loader traces, logs and hashes are
retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-window-first-open-publication).

This closes construction and first open for the neutral host-owned controller.
Borrowed-controller attachment/root establishment, headless context and camera
setup, Qt widget policy, arbitrary camera/cursor setters, feature-store
reparenting and Qt faceplate composition remain separate boundaries. Resource,
platform and graphical gates are unchanged.

### CXX-SOURCE-048 — default and borrowed controllers satisfy one host attachment invariant

**Repaired and qualified for default-controller construction and neutral host
attachment.** The default `BObolViewController` previously returned without a
scene root. CXX-SOURCE-047 established a root only while constructing a
host-owned controller, leaving a directly constructed or endpoint-owned
controller incomplete. Attaching that controller succeeded without a result,
and its first `open()` published the root independently from the host
descriptor and viewport. Invalid null-root or non-group-root controllers were
also accepted after the host had already closed its current attachment.

Both public controller constructors now call one graph initialization helper.
Temporary unique node owners cover every fallible allocation and commit the
render batch, presentation roots and layer roots only after the complete graph
exists. The default form adds an empty `SoBRLSceneGroup` and retires its
construction-only request, while both forms establish the three-copy 1x1
viewport invariant. `BObolWindowHostPrivate` prepares its descriptor before
constructing a locally owned controller and releases that controller only when
the private state is complete.

`attachController()` returns `SbBool` and validates the candidate group root
before closing or changing the host. Reattaching the same controller does not
close, allocate or notify. Production headless, display-session, display-
endpoint, GED framebuffer and Qt factory callers propagate the result. The C
endpoint constructor uses the common default constructor, contains all C++
exceptions and transfers no ownership on failure. Destruction unsubscribes
from the LoD service before directly cancelling its generation and releasing
resident demand; it does not enter the normal allocating live policy reducer.

The exact CXX-SOURCE-047 library rejects the focused selector because the
default controller lacks a hostable root. Current coverage passes 326
controller construction positions (230 rejected and 96 complete), 233 C
endpoint positions (232 contained rejections and one complete endpoint), two
allocation-free invalid-root cases and 20 ordinary/automatic first-open
positions. Identical attachment and traced service teardown each allocate
zero times; both display modes cover callback failure after commit. Twelve
neighboring publication selectors pass, as do eleven prototype/LoD/image/RT/
host/GED/Qt/qged tests, static link and all seven conformance checks. The exact
full executable passes 149 direct checks, 21 scenarios and 1,044 rows in
360.27 seconds. Current compact/shared/static hashes are `1888ab0e...`,
`922f2260...` and `60a7349c...`; Obol and OSMesa remain `21872f7f...` and
`097aac27...`. Exact binaries, dependencies, sources, loader traces, logs and
hashes are retained in the
[dated evidence](obol_20260919_readiness_history.md#september-10-window-controller-attachment).

This closes default and borrowed-controller root establishment at the neutral
host boundary. Headless context/camera setup, Qt widget lifetime policy,
camera/cursor setters, feature replacement/reparenting, multi-record
controller migration, forced-failure borrowed framebuffer retirement and Qt
faceplate composition remain separate. Resource, platform and graphical gates
are unchanged.

### CXX-SOURCE-049 — headless camera state exists before host lifecycle publication

**Repaired and qualified for default-controller and headless-host camera
ownership.** The default controller previously had a group root but no camera.
Headless `open()` compensated by searching the modeled scene and adding a new
camera there before calling the controller setter. `SoViewport` separately
parents that same camera in its private root and explicitly requires its scene
input to exclude cameras. The adapter therefore owned a scene mutation and
could traverse the same camera through two paths during first open.

The default constructor now builds its orthographic camera under temporary
ownership, then installs it with the already prepared empty scene. Controller,
viewport and render-manager camera pointers agree before construction returns;
the modeled root does not contain that camera and construction-only host work
is cleared. The explicit root/camera constructor retains its caller-selected
semantics.

The headless adapter has no camera builder or scene search. Open and render
validate the attached group root and active camera. A camera-less controller
rejects before descriptor copying and provider binding, including when a
camera is incorrectly present in its modeled scene. The built-in factory runs
the same validation before attachment and contains binding exceptions.

The exact CXX-SOURCE-048 library rejects the focused selector because its
default controller lacks an active camera. Current focused coverage checks the
three camera pointers, private-root membership, modeled-scene exclusion,
successful-open identity and two allocation-free rejection forms. The
strengthened construction selector passes 380 controller positions (259
rejected and 121 complete) and 262 endpoint positions (261 contained
rejections and one complete). Thirteen neighboring selectors, eleven
integration tests, static link and all seven conformance checks pass. The full
executable passes 150 direct checks, 21 scenarios and 1,045 rows in 332.96
seconds. Current compact/shared/static hashes are `da0b8783...`, `9e40d9db...`
and `e6ec6b5c...`; Obol and OSMesa remain `21872f7f...` and `097aac27...`.
Exact evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-10-headless-camera-invariant).

This closes headless camera creation and adoption. Context-provider
publication during headless open, arbitrary live camera replacement, Qt widget
policy and the remaining lifecycle/resource rows stay separate.

### CXX-SOURCE-050 — headless provider and host policy publish as one lifecycle transition

**Repaired and qualified for headless open, provider replacement and close.**
The preceding headless `open()` called the controller's provider setter before
the base host open. Renderer invalidation published a callback before the
action provider, descriptor, viewport and open flag changed. Close published
the inverse partial state, and `setContextManager()` changed the host pointer
before the controller operation could reject allocation.

The controller now prepares the renderer request before mutation, binds the
action provider, commits renderer evidence and progressive work without an
intermediate wake, then publishes the selected final request. Base host open
has matching prepare, commit and notify phases. Headless open composes both
transactions and lets a changed viewport's capacity request subsume the
renderer-performance request. Replacement and close commit both host and
controller policy before notification. The built-in factory defers provider
binding until first open and uses the same host transaction when reasserting an
already open instance.

`SoGLRenderAction` previously removed its registered manager/key and changed
its pointer/key before inserting the replacement into an allocating
`unordered_map`. Obol now registers first. A first manager/key insertion can
fail with the old action intact; manager replacement updates the existing key
without allocation; key replacement retires the preceding key only after the
new entry exists.

The exact CXX-SOURCE-049 host stack rejects the focused selector because its
callback observes a closed host. A second retained stack containing the
controller/host transaction but the preceding Obol library rejects the
allocation-free registered-provider check. Current coverage passes two
cache-key positions, 14 ordinary/automatic open positions and 16 live
replacement/close positions. Every allocation result is the complete old or
new policy, exact no-ops are allocation-free, and throwing callbacks see the
committed state.

Fourteen neighboring selectors, twelve integration tests, three public/install
contract checks, static link, 27 native render-action tests and all seven
conformance checks pass. The full executable passes 151 direct checks, 21
scenarios and 1,046 rows in 356.62 seconds. Current compact/shared/static
hashes are `34a85acc...`, `ad9d0a4f...` and `b4644c79...`; Obol is
`4e36cf66...` and OSMesa remains `097aac27...`. Exact evidence is retained in
the [dated readiness record](obol_20260919_readiness_history.md#september-10-headless-context-provider-publication).

No formal relation or accepted baseline changes. This closes headless provider
publication; live camera replacement, Qt widget policy and the remaining
lifecycle/resource rows stay separate. S1--S6 remain open.

### CXX-SOURCE-051 — feature record, node and root publish as one replacement

**Repaired and qualified for existing-feature replacement and reparenting.**
The feature store previously mutated its live record before rebuilding its
node, detached the old node before attaching its successor, and reordered an
overlay in a third step. Immediate root, path and result observers could see
an empty root or a new record with the preceding graph. Allocation failure
could also leave the new metadata paired with the old node. A custom node
shared by two records could disappear when either record replaced it.

Existing-record publishers now copy and configure a candidate record away
from live state. They retain the candidate node, calculate the effective
ordered children for every old and new attachment root, and prepare all child-
list replacements before mutation. The root edits and record pointer then
commit without allocation; the preceding record releases its node only after
the successor is installed. Root/path notification and presentation requests
run afterward, preserving the committed successor when a callback throws.
One overlay comparator and unique-node helper serve both ordinary reorder and
transactional replacement. Typed publishers, custom nodes, edit previews,
style/geometry rebuilds and non-mesh selection rebuilds share the same
publication primitive.

The exact CXX-SOURCE-050 library rejects the focused selector with a partial-
state diagnostic. Current focused coverage exercises typed replacement,
shared-custom-node replacement and main-root to screen-root reparenting. It
passes 108 forced-allocation outcomes: 98 preserve the complete predecessor
and ten publish the complete successor. Retry, exact reparent no-op, root/path
observers, command-result observation and throwing observers also pass. Nine
neighboring publication selectors, ten feature/GED/Qt integration tests, three
public/install contract checks, static link and all seven conformance checks
pass. The exact full executable passes 152 direct checks, 21 scenarios and
1,050 rows in 356.73 seconds. Current compact/shared/static hashes are
`c838a8d8...`, `8f94eaba...` and `60b37629...`; Obol and OSMesa remain
`4e36cf66...` and `097aac27...`. Exact evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-10-feature-replacement-and-reparenting-publication).

No formal relation or accepted baseline changes. This closes existing feature
replacement/reparenting and `setOverlayInfo()` record/node/root/order
composition. Multi-record controller migration, initial feature insertion,
in-place metadata and mesh-field edits, Qt/lifecycle/resource gates and S1--S6
remain open.

### CXX-SOURCE-052 — feature controller migration publishes all roots and ownership together

**Repaired and qualified for replacement, detachment and adoption of a feature
store controller.** `BObolFeatureStore::setController()` previously detached
and attached each record in map order, changed each attachment pointer
immediately, published the store's controller only after every record moved,
and sorted overlays afterward. Root and path observers could see any migration
prefix paired with the preceding controller. Allocation or callback failure
could leave records split across controllers and roots, with an incomplete
overlay order and no reliable presentation revision.

The store now projects every retained record into its desired controller root
first. One root-order composer preserves unrelated children, reduces aliased
custom nodes to one occurrence, and applies the same overlay ordering used by
single-record replacement. Every affected old and new child-list replacement
prepares before mutation. The root orders, store controller, record attachment
pointers and one presentation revision then commit without allocation. Root
and path observers run afterward, followed by requests to both the preceding
and successor controllers. Every notification is attempted while the first
exception is retained. Rebinding the current controller allocates, changes and
notifies nothing.

The exact CXX-SOURCE-051 library rejects the focused selector because an
observer sees partial roots or ownership. Current coverage uses five records
across model and screen roots, two records which alias one custom node, four
unrelated root extensions and every retained old feature path. Controller
replacement, detachment to no controller and adoption from no controller pass
170 forced-allocation outcomes: 159 preserve the complete predecessor and 11
publish the complete successor. Exact no-op, retry, all old/new root
observers, both controller frame callbacks and exceptions from graph or frame
callbacks also pass. Ten neighboring publication selectors, ten
feature/GED/Qt integration tests, three public/install contract checks, static
link and all seven conformance checks pass. The exact full executable passes
153 direct checks, 21 scenarios and 1,054 rows in 360.67 seconds. Current
compact/shared/static hashes are `d933b9a1...`, `d6c1a1b5...` and
`140a4b1f...`; Obol and OSMesa remain `4e36cf66...` and `097aac27...`.
Exact evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-feature-controller-migration-publication).

No formal relation or accepted baseline changes. This closes multi-record
feature-controller migration. Initial feature insertion, in-place metadata and
mesh-field edits, Qt faceplate/lifetime work, other lifecycle/resource gates
and S1--S6 remain open.

### CXX-SOURCE-053 — Qt canvas validation precedes host identity and ownership changes

**Repaired and qualified for invalid canvas/controller replacement.**
`QgObolWindowHost::setCanvas()` previously destroyed its owned live canvas,
published the candidate canvas pointer and only then asked the base host to
attach the candidate controller. A missing or non-group scene root made the
base host reject that controller after the destructive changes. The Qt host
then exposed the candidate canvas paired with its preceding, and potentially
already destroyed, controller.

The Qt adapter now applies one shared hostability predicate before either
`setCanvas()` or `bindController()` changes canvas, controller or ownership.
`setCanvas()` reports acceptance, and invalid input retains the complete
borrowed or owned predecessor. The same-canvas case is an exact no-op. The
temporary owned-canvas preservation used by factory controller binding is now
scope-bound, so an exception cannot leave later close operations permanently
configured to retain the canvas.

The exact CXX-SOURCE-052 libqtcad rejects the new focused executable because
the borrowed replacement exposes different canvas and controller identities.
Current coverage checks direct controller rejection, borrowed and owned
canvas rejection, retained owned-canvas lifetime, exact no-op and accepted
replacement retirement. The focused executable passes, as do twelve Qt
functional tests, Qt plugin reuse, qged framebuffer-host integration and the
staged public header copy. The System-GL-only host row reaches its declared
environment skip on this seat. Exact focused executable and predecessor/current
libqtcad hashes are `d9d6b275...`, `ef3a5c4e...` and `6709cefc...`;
libBObol, Obol and OSMesa remain `d6c1a1b5...`, `4e36cf66...` and
`097aac27...`. Evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-qt-host-rebinding-validation).

No formal relation or accepted baseline changes. This closes validation and
ownership preservation for rejected Qt canvas/controller inputs. A valid live
rebind interrupted while closing framebuffer attachments or changing canvas
controller resources remains part of the separate forced-failure lifecycle
boundary. Qt faceplate synchronization, widget policy, initial feature
insertion, in-place metadata/mesh-field edits and S1--S6 remain open.

### CXX-SOURCE-054 — Qt fallback faceplate nodes prepare before live publication

**Repaired and qualified for direct Qt grid, axes and ADC node updates.** A Qt
canvas without a GED owner previously configured its attached faceplate nodes
in place. Grid center, spacing, divisions, view inputs and derived geometry,
for example, notified separately. An immediate observer could see a new center
paired with the old spacing and divisions before `rebuildGeometry()` replaced
the derived child graph.

The fallback adapter now constructs each grid, axes or ADC successor away from
the live root. One shared helper composes the root's complete next child list,
prepares it before mutation, commits it once and then attempts both graph
notification and the controller frame request while retaining the first
exception. The same helper handles insertion, replacement and removal. A
single typed child lookup replaces the three duplicate scans. An unchanged
view/frame hash still preserves every existing node and geometry identity.

The exact CXX-SOURCE-053 libqtcad rejects the focused executable because an
immediate sensor observes partially updated grid input fields. The current
executable observes only a complete predecessor or successor and verifies the
final adaptive grid geometry, unchanged-refresh identity, axes/ADC mapping and
removal. Twelve Qt functional tests, Qt plugin reuse, qged framebuffer-host
integration and the graphical GED faceplate row pass. The System-GL-only host
row reaches its declared environment skip on this seat. Exact focused
executable and predecessor/current libqtcad hashes are `0d00e4c0...`,
`6709cefc...` and `9ca01e87...`; libBObol, Obol and OSMesa remain
`d6c1a1b5...`, `4e36cf66...` and `097aac27...`. Evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-qt-fallback-faceplate-node-publication).

No formal relation or accepted baseline changes. This closes each direct Qt
fallback grid/axes/ADC node publication. Whole-faceplate multi-record
composition in the GED adapter, direct callers which configure an already
attached grid, axes/ADC public input setters, widget policy and the remaining
lifecycle/resource gates stay open.

### CXX-SOURCE-055 — initial feature publication commits identity, graph and frame obligation together

**Repaired and qualified.** First publication previously allocated an identity
and inserted the live record and name entries before configuring its payload,
building its retained node or attaching that node. Allocation failure could
therefore consume an unpublished identity or leave a visible partial record.
Even successful root notification preceded the feature presentation revision
and controller render request; failure while preparing that request could
strand committed content whose identical retry was a no-op.

Insertion now builds the complete candidate off-store while advancing only a
local copy of the next identity. It prepares detached map nodes for both
indexes, the complete root child-list replacement and a controller presentation
request before the first live change. The allocation-free commit installs both
indexes, root order, retained ownership, next identity, presentation revision
and request level before graph callbacks. Root and frame notifications are
both attempted while retaining the first exception. The request retains its
exact controller across callback-time store rebinding. Existing-record
publication uses the same prepared request boundary, so its observers also see
the committed presentation obligation.

The exact CXX-SOURCE-054 libBObol rejects the current focused executable at
`feature insertion failure left a partial record, graph or request`. Current
coverage passes 138 allocation positions: 48 preserve the exact empty store
and unpublished identity, while 90 expose only a complete successor. Immediate
graph, frame and command-result observers see the complete record, node, root,
presentation revision and standing request; throwing observers leave that
successor committed. Neighboring feature publication checks, ten feature/GED/
Qt integration tests, three public/install contract tests and the static
consumer pass. The complete compact executable passes 154 direct checks, 21
scenarios and 1,055 PASS rows. Exact focused executable,
predecessor/current shared libBObol, current static libBObol, Obol and OSMesa
hashes are `c933c619...`, `d6c1a1b5...`, `7dbacbea...`, `dc42e322...`,
`4e36cf66...` and `097aac27...`. Evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-initial-feature-publication).

No formal relation or accepted baseline changes. This closes first-record
publication and its render obligation. Whole-faceplate multi-record GED
composition, in-place metadata/mesh-field edits, feature removal/clear,
controller-migration request preparation and the remaining lifecycle/resource
gates stay separate.

### CXX-SOURCE-056 — GED faceplate records publish as one logical composition

**Repaired and qualified for the GED faceplate feature set.** A faceplate
refresh previously called the individual feature publishers and removers in
sequence. Center dot, interactive rectangle, grid, ADC, parameter text, LoD
progress, scale and axes could therefore cross a graph or frame notification
boundary as a mixture of predecessor and successor records.

`BObolFeatureStore::applyPublication` now validates and stages a mixed set of
line, HUD-label and custom-node replacements and removals. Complete candidate
records and retained nodes, new name/record index nodes, every affected final
root order and one controller request prepare before the store changes. Commit
installs the records, indexes, root orders, node ownership, next identity,
presentation revision and request without allocation. Exact unchanged input is
a callback-free no-op. Invalid actions, scopes, kinds, null custom nodes and
duplicate keys reject while retaining the predecessor. Graph, frame and result
callbacks all run after the complete successor is visible, and callback failure
does not roll it back.

The GED adapter accumulates its complete faceplate state before one call to the
new primitive. Against the exact CXX-SOURCE-055 libraries, a saved clean probe
records 18 immediate graph callbacks and observes a partial state before the
final complete state. Current direct coverage passes 373 allocation positions:
231 preserve the exact predecessor and 142 expose only the complete successor.
The real GED selector verifies atomic insertion followed by replacement/removal,
overlay metadata and graph/frame observation. Neighboring feature checks,
eleven feature/GED/Qt integration tests, three public/install contract tests,
the static consumer and seven conformance tests pass. The complete compact
executable passes 155 direct checks, 21 scenarios and 1,056 PASS rows. Exact
direct/GED executables, predecessor/current libBObol and predecessor/current
libged hashes are `42bfc460...`, `f1941c87...`, `7dbacbea...`,
`01c58527...`, `5abb7712...` and `d70e5358...`. Evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-whole-faceplate-ged-composition-publication).

No formal relation or accepted baseline changes. Framebuffer composition and
presentation still occur after this feature transaction, and the public
`ged_view_feature_batch` adapter remains sequential. In-place metadata and mesh
fields, controller-migration request preparation and the remaining lifecycle/
resource gates stay separate.

### CXX-SOURCE-057 — feature removal and clear commit before observation

**Repaired and qualified for the public feature-store APIs.** Handle removal
previously erased its name index before invoking the command-result callback,
then detached its node, erased its record and finally requested presentation.
Any observer or exception could expose a partial retirement. Prefix and scope
removal replayed that sequence for each match, while clear notified each still-
live record and reset command-owner generations only after the record loop.

All handle, name, owned-name, prefix, scope and public clear operations now
prepare removals through the feature publication transaction. Each complete
removal set, final affected root order, result notification list and one
reason-preserving controller request prepare before mutation. The allocation-
free commit retires indexes, records and retained nodes together, advances the
presentation once and clears owner generations inside the clear commit. Graph,
frame and every result callback then see the complete successor; all are
attempted while retaining the first exception. Empty and unmatched operations
do not change presentation state or notify observers.

The focused executable fails against the exact CXX-SOURCE-056 libBObol at
`feature removal measurement exposed a partial publication`. Current handle,
prefix, scope and clear coverage passes 439 allocation positions: 228 preserve
the exact predecessor and 211 expose the complete successor. It also exercises
throwing graph/frame/result callbacks, survivor identity and order, clear-time
generation reset, generation-only clear and unmatched no-ops. Four neighboring
feature selectors, thirteen feature/GED/Qt integration tests, three public/
install contract tests, the static consumer and seven conformance tests pass.
The complete compact executable passes 156 direct checks, 21 scenarios and
1,061 PASS rows. Exact focused executable and predecessor/current shared
libBObol hashes are `b6e41cf6...`, `01c58527...` and `edb42f3e...`. Evidence
is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-feature-removal-and-clear-publication).

No formal relation or accepted baseline changes. At this checkpoint,
destructor-time feature teardown, direct in-place metadata and mesh fields,
generic GED batch composition and the remaining lifecycle/resource gates stayed
separate.

### CXX-SOURCE-058 — feature-controller migration includes both render obligations

**Repaired and qualified.** CXX-SOURCE-052 made migration root orders,
controller ownership and attachment roots one prepared commit, but the previous
detach and next attach render requests were still created afterward. Immediate
root observers therefore saw migrated content before either standing request,
and request-preparation failure could leave the complete graph without its
repaint obligation.

Feature request preparation is now controller-independent and shared by
ordinary feature publication and migration. Both applicable controller
requests and their LoD transition scopes prepare after all affected child-list
replacements and before mutation. Commit publishes root orders, controller and
attachment ownership, presentation revision and both request states without
allocation. Graph notification follows that complete state. Both frame
notifications are attempted while retaining the first exception, including
when graph or the previous-controller frame callback throws. Same-controller
rebinding remains an exact no-op.

The enhanced existing selector fails against the exact CXX-SOURCE-057
libBObol at `feature controller migration exposed partial roots or ownership`.
Current replace, detach and attach sweeps pass 170 allocation positions: 167
preserve their exact predecessors and three expose their complete successors.
They verify standing requests, exact detach/attach reasons, root and path
observation, callback continuation and same-controller no-op behavior. Five
neighboring feature selectors, thirteen feature/GED/Qt integration tests, three
public/install tests, the static consumer and seven conformance tests pass. The
complete compact executable passes 156 direct checks, 21 scenarios and 1,061
PASS rows. Exact focused executable and predecessor/current shared libBObol
hashes are `c5aace3a...`, `edb42f3e...` and `c3d93019...`. Evidence is
retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-feature-controller-migration-request-publication).

No formal relation or accepted baseline changes. At this checkpoint,
destructor-time feature teardown, direct in-place metadata and mesh fields,
generic GED batch composition and the remaining lifecycle/resource gates stayed
separate.

### CXX-SOURCE-059 — direct feature records and mesh fields prepare before publication

**Repaired and qualified.** The remaining retained-mesh edit paths changed
their live feature record and then assigned one or more `SoMField` arrays before
advancing the mesh source revision and requesting presentation. Allocation
failure or an immediate field, node, frame or result observer could therefore
see a new record with old field data, a partially replaced field set, or a
complete edit without its standing render request. `replaceMetadata()` and
`replacePrimitiveMetadata()` likewise assigned directly into the published
record before constructing their result notification.

Obol now provides `SoMField::prepareValueReplacement()`. It validates the
concrete field type, copies the complete successor into a detached same-type
field and suppresses target notification during preparation. `commit()` swaps
the owned value storage, size, capacity and external-storage ownership flag
without allocation or notification; `notify()` restores the target's original
notification state and reports the completed value once. Abandoning an
uncommitted replacement restores notification and leaves the target unchanged.
The target lifetime and no-intervening-mutation requirement is explicit in the
public contract.

The feature store now copies the record and prepares every changed multi-field,
the controller presentation request and the fully constructed command result
before the first live change. One allocation-free commit installs the record,
field arrays, mesh source revision, feature presentation revision and standing
request while retaining the existing mesh node. Field, node, frame and result
notifications then observe the complete successor; all are attempted while the
first exception is retained. Metadata-only changes use the same candidate-
record and prepared-result boundary without manufacturing a visual effect.

The enhanced selector fails against the exact CXX-SOURCE-058 libBObol at
`direct feature mutation exposed partial record, field or request`. Current
indexed-point, selected, highlighted, metadata and primitive-metadata sweeps
pass 529 allocation positions: 144 retain the exact predecessor and 385 expose
the complete successor. They cover field-array growth and shrinkage, retained
feature, field-object and node identity, exact presentation reasons, standing requests,
immediate observers and graph, frame and result callback failure. Six
neighboring feature selectors, thirteen feature/GED/Qt integration tests,
three public/install tests, the static consumer and seven conformance tests
pass. Native Obol passes its focused field replacement test, 1,580 unit tests
with one existing skip and all 63 integration tests. The complete compact
executable passes 157 direct checks, 21 scenarios and 1,067 PASS rows in
357.81 seconds. Evidence is retained in the
[dated readiness record](obol_20260919_readiness_history.md#september-12-direct-feature-record-and-mesh-field-publication).

No formal relation or accepted baseline changes. Externally owned custom-node
mutation, destructor-time feature teardown, generic GED batch composition and
the remaining lifecycle/resource gates stay separate.

### CXX-SOURCE-060 — GED edit transactions publish one complete retained preview

**Repaired and qualified.** `ged_view_edit_transaction_apply()` previously
implemented one logical edit update as a sequence of public operations: ensure
an empty preview, publish its overlay, publish its color, replace its geometry
and semantic revisions, then call a feature-store `touch()`. Immediate graph
observers could see these intermediate generations. A first update advanced the
retained presentation revision three times and depended on an uncomposed direct
record/node notification which did not itself own a frame obligation.

The adapter now finishes primitive transformation and wireframe plotting before
retained mutation, constructs one complete `BObolFeaturePublication`, and calls
the feature store's atomic publication operation once. Geometry, normalized
commands, style, identity, edit intent, source/input revisions, overlay
metadata, retained node and render request therefore become observable as one
successor. Identical live events are exact no-ops. Commit, cancel and discard
retire the preview directly and do not plot stale primitive input. The public
and private feature-store touch entry points were deleted; independent
single-property edit APIs retain their existing single-operation contracts.

The dynamic probe against the exact CXX-SOURCE-059 libged and libBObol observes
three graph callbacks, one partial callback and a retained revision change from
0 to 3. The current libraries expose one complete graph callback, one frame
request and revision 0 to 1. Focused GED coverage proves the same boundary for
caller-supplied points and transformed `rt_db_internal` geometry, replacement,
an exact no-op and terminal retirement. The feature-store batch sweep embeds an
edit preview in mixed replace/remove/insert publication and passes 570 forced-
allocation positions plus graph, frame and result callback failure.

No formal relation or accepted baseline changes. Externally owned custom-node
mutation and revision notification, destructor-time feature teardown, the
general public `ged_view_feature_batch` composition and the remaining
lifecycle/resource gates stay separate.

### CXX-SOURCE-061 — initial custom nodes publish configured overlay state in one commit

**Repaired and qualified.** The qged ARB, ellipsoid, sketch and BoT edit
presenters and the libged navigation gizmo previously published a custom node
and then classified it as an overlay through `setOverlayInfo()`. Immediate
observers could see a retained node without its role, class, lifecycle, order
or source path. ARB, ellipsoid and sketch also attached a default manipulator
before setting its topology or axes, selection domain and session revision.

`BObolFeatureStore` now has a custom-node publication overload which accepts
the complete overlay value. It uses the same prepared insertion/replacement
primitive as other feature publications, so record, node, ordered root,
presentation revision and render request publish together. The five production
callers construct their overlay before publication. ARB, ellipsoid and sketch
configure a detached node before passing it to the store, and a newly published
node relies on the transaction's one frame request. Existing live manipulators
retain their explicit mutable-node update request.

The dynamic probe against the exact CXX-SOURCE-060 libBObol observes two graph
generations, one partial callback and two presentation-revision advances for
the old publish-then-overlay sequence. Current exposes one complete graph
callback, one frame callback and one revision advance. Permanent focused
coverage verifies that graph, frame and result callbacks each observe the
complete custom-node and typed-overlay successor exactly once. Its 75 forced-
allocation positions preserve the empty predecessor 58 times and expose the
complete successor 17 times. The feature-store regression also proves that the
compatible old overload preserves an existing overlay during custom-node
replacement. The complete compact executable passes 158 direct checks, 21
scenarios and 1,068 PASS rows in 316.99 seconds. Fourteen integration/API
tests, seven conformance tests, header self-containment and the static consumer
pass.

No formal relation or accepted baseline changes. Later mutation of an
externally owned published node, including its revision and notification
contract, remains a distinct boundary. Destructor-time feature teardown,
general public `ged_view_feature_batch` composition and wider lifecycle and
resource qualification also remain open.

### CXX-SOURCE-062 — edit-manipulator setters publish one complete node state

**Repaired and qualified at the node-local boundary.** The edit-manipulator
setters changed retained private geometry or public selection fields before
calling `rebuildGeometry()`. Geometry construction can allocate, so failure
could leave new topology, axes, hover, active, domain or visibility state
paired with the preceding children. Immediate field and node observers could
also see that incomplete prefix. `SoBRLIndexedEditManipulator::setTopology()`
was the most severe case: it cleared the live topology before copying the
successor, despite its documented atomic contract.

Both manipulator classes now build owned successor children away from the live
node and prepare the child-list replacement before mutation. One common helper
commits that replacement with the new private geometry and quiet scalar fields,
restores notification policy, then attempts path, node and field callbacks
against the complete successor while retaining the first exception. Indexed
topology builds in a private candidate and swaps only after preparation;
queries and rendering share the candidate's point, edge, face and feature-
position helpers. Exact scalar setter no-ops allocate and notify nothing.

The permanent setter selector covers axis geometry, hover and visibility plus
indexed topology, multi-field selection-domain reset, hover and visibility. It
passes 767 forced-allocation outcomes: 759 preserve the complete predecessor
and eight publish the complete successor. Retry, no-op, retained child paths,
node/path observation and throwing callbacks pass. The same current test
executable forced to load the exact CXX-SOURCE-061 library rejects with
`edit manipulator setter left partial state`. The neighboring rebuild selector
passes 449 positions. The complete compact executable passes 159 direct
checks, 21 scenarios and 1,069 PASS rows in 315.45 seconds. Fourteen
integration/API tests, seven conformance tests, public-header self-containment
and the static consumer pass.

No public or formal contract changed. This closes atomicity inside the two edit
manipulator node classes. Coordinating a published custom node's mutation with
its feature-store revision and frame obligation remains separate, as do the
BoT shared-mesh and navigation-gizmo live-update paths.

### CXX-SOURCE-063 — qged manipulator updates publish one store-visible successor

**Repaired and qualified for the small qged manipulators.** After initial
publication, the ARB, ellipsoid and sketch adapters changed their externally
owned custom node through several setters and then requested a frame. Although
`CXX-SOURCE-062` makes each setter internally atomic, one logical geometry,
selection or pointer update still exposed several complete node generations
without advancing the feature record. The retained feature revision therefore
stayed at one while the graph-visible custom node changed underneath it.

The three adapters now build a fully configured detached manipulator and replace
the named custom-node feature through the `CXX-SOURCE-061` transaction. Node,
topology or axes, selection, hover/active state, overlay ownership, ordered root,
feature revision and frame obligation publish as one successor. ARB and sketch
read authoritative session selection before constructing updated geometry, so
they no longer publish geometry and then selection separately. Stable feature
identity is retained, but plugin state no longer borrows a node pointer: every
use resolves the current record and refreshes its revision. Owner-result
callbacks perform the same resolution after commit, while reentrant graph/frame
work can already reach the committed node by identity.

The production qged replay now requires the ellipsoid, ARB and sketch feature
records to advance after live interaction. Exact preceding plugin binaries fail
at the first ellipsoid check with feature revision one; current passes the full
single-view and quad-view primitive-edit replays. The view-store, edit-
manipulator, public-symbol, GED draw-sync and Qt edit-preview checks also pass.
The evidence and independent verifier are retained under
`.build/obol-qualification/20260912-qged-manipulator-publication`.

No public libBObol API or formal relation changed. Whole-node replacement is
bounded to these small manipulators. The BoT shared-mesh hot path and navigation
gizmo remain separate live-update boundaries.

### CXX-SOURCE-064 — BoT shared-mesh edits publish coordinated surface successors

**Repaired and qualified for interactive BoT point moves.** The BoT adapter
mutated its large shared mesh directly, cleared each already-published wrapper
and requested a frame. Geometry therefore changed without a new feature record,
and selected-face state was absent from that record. Simply replacing wrappers
after the mutation was also insufficient: allocation failure could leave the
predecessor records pointing at the changed shared payload, while per-view
publication callbacks could run before the other views committed.

The feature store now supports one coordinated publication across distinct view
stores. Every record, wrapper, root replacement and frame request prepares
first. An optional allocation-free commit hook then changes shared state, every
store commits, and only afterward do graph, frame and result callbacks run.
Custom mesh publication records carry exact selected and highlighted primitive
state and reject a mismatch with the detached node.

QBot uses this boundary to retain one heavy mesh and topology while constructing
one small detached wrapper per view. After all view stores prepare, the commit
hook changes only the selected one to three points, clears stale normals and
advances the shared source identity. Each view then installs a new wrapper with
the same shared-geometry pointer, the current selected face and the new source
identity. New wrapper identity invalidates its private render cache without a
topology copy.

The permanent two-store test denies 128 allocation positions: 79 preserve both
stores and the shared geometry, while 49 expose the complete successor in both
stores. Immediate root and frame observers see both stores committed. The
existing mixed feature batch still passes 602 positions. The complete compact
executable passes 160 direct checks, 21 scenarios and 1,070 PASS rows in 315.28
seconds. The exact preceding BoT plugin fails the focused GUI replay because it
keeps the old wrapper; current passes that replay and the full single-view and
quad-view primitive-edit replays.

This adds a public coordinated feature-publication operation but no formal
relation. The separate GED point-handle batch and cross-controller feature
composition remain open. The GUI fixture proves retained heavy-geometry
identity, not large-mesh performance. Navigation-gizmo live mutation remains
the next custom-node boundary.

### CXX-SOURCE-065 — navigation-gizmo state publishes as retained successors

**Repaired and qualified for semantic gizmo state.** The libged navigation
plugin retained a borrowed `SoBRLNavigationGizmo *` and changed its public
fields and HUD children directly for style, visibility, hover, press, release
and drag updates. Camera changes also reached the published node through its
field sensor or render traversal. The feature record therefore remained at its
preceding revision while one logical interaction could expose several node
generations and manually request a later frame.

Plugin state now owns one value snapshot and resolves the current custom node
through its stable feature identity. Each presentation change constructs a
fully configured detached gizmo with a fixed camera-orientation snapshot, then
submits one complete custom-node, style, overlay and owner publication. An
allocation-free coordinated-publication hook commits the authoritative plugin
snapshot before the feature store commits; graph, frame and result observers
therefore see matching state and node. Exact presentation no-ops publish
nothing. Gesture bookkeeping which does not affect presentation remains in the
private snapshot. Hosts without a GED view-update callback receive an explicit
endpoint camera synchronization before the camera successor is published.

`SoBRLNavigationGizmo::setCameraOrientationSnapshot()` detaches live camera
tracking and replaces existing orientation-dependent children through the
node's failure-atomic child transaction. Preparation failure restores the
camera, orientation and HUD predecessor; observer failure after commit retains
the complete successor. The geometry builders retain every detached Coin node
through construction, closing the corresponding forced-allocation leaks.

The expanded node selector passes **1,941 forced-allocation positions**: 1,936
preserve a complete predecessor and five expose a complete successor. It also
checks node/path observers, retry, exact hidden no-op behavior and callback
failure. The exact preceding plugin fails all thirteen new integration
assertions for fixed-camera state, exact record advancement and predecessor
stability; current passes. The complete compact executable passes **160 direct
checks, 21 scenarios and 1,070 PASS rows in 335.65 seconds**. Eight adjacent
libBObol, GED and Qt tests pass, as do public-header self-containment, static
linking, the full qged primitive-edit replay and a focused navigation-gizmo
replay in OSMesa single and quad layouts. Evidence is retained under
`.build/obol-qualification/20260912-navigation-gizmo-publication`.

No formal relation changed. Renderer-derived viewport anchoring remains a
render-local layout update rather than semantic feature state. Input-layer
installation and removal are owner-thread lifecycle effects ordered around the
feature publication; no concurrent cross-system transaction is claimed.
Replacement of the controller's camera object while the plugin is completely
idle remains part of the general live-camera boundary. The separate GED point-
handle batch, cross-controller feature composition and wider lifecycle/resource
gates remain open.

### CXX-SOURCE-066 — live camera replacement publishes one controller state

**Repaired and qualified for direct controller camera replacement.** The
preceding controller called `SoViewport::setCamera()` before changing its own
camera pointer, the render manager, the LoD view revision and the render
request. Native viewport replacement itself removed the preceding camera and
notified its paths before inserting the successor. Immediate graph and path
observers could therefore see an empty root or the new viewport camera paired
with the controller's preceding camera; allocation or callback failure could
interrupt the remaining effects.

`SoViewport::prepareCameraReplacement()` now retains both generations and
prepares the complete final child order. Its allocation-free commit changes
the root and viewport camera without notification; notification follows only
after the complete native state is visible. `BObolViewController::setCamera()`
prepares the capacity-level frame request and native replacement, computes the
candidate LoD signature without waking the endpoint, and then commits the
viewport, controller and render-manager camera pointers, LoD revision and
standing request before root/path or frame callbacks. Notification drains the
native observers and endpoint wake while retaining the first exception. Exact
pointer no-ops allocate and notify nothing.

The libged navigation-gizmo plugin now observes the viewport root as the
camera-generation completion edge. It retargets its existing orientation
field sensor before publishing the fixed-orientation gizmo successor. A
retired camera can no longer publish another feature generation, while an
orientation change on the replacement camera still publishes exactly one.

The fresh focused selector passes **181 forced-allocation positions** across
replacement, removal and addition with ordinary and automatic LoD: 174 retain
the complete predecessor and seven expose the complete successor. In the
warmed full matrix it covers 180 positions, with 174 predecessors and six
successors. The exact CXX-SOURCE-065 controller/native pair fails the complete-
state assertion, and the exact preceding gizmo plugin fails the two replacement
and retired-camera assertions against the repaired controller. Current native
Obol passes all 13 viewport tests. The exact full compact executable passes
**161 direct checks, 21 scenarios and 1,071 PASS rows in 357.19 seconds**.
Affected integration tests, public-header compilation, static linking and
exported-symbol checks pass. Exact evidence is retained under
`.build/obol-qualification/20260912-live-camera-replacement`.

No formal relation changed. This boundary covers direct `setCamera()`
replacement and its known retained-graph consumer. The enclosing
`syncCameraFromViewContext()` calculation and its camera, lighting, clipping
and affordance input-field composition remain a separate classified boundary,
as do Qt widget policy, cross-controller feature composition and the wider
lifecycle/resource gates.

### CXX-SOURCE-067 — view-camera input publishes one derived controller state

**Repaired and qualified for `syncCameraFromViewContext()`.** The preceding
operation replaced the projection camera before configuring it, adopted an
initial viewport through another public transaction, updated the tracked
lights and three clip nodes field by field, rebuilt the cutting-plane aid
against the preceding camera, and only then updated the camera, LoD signature
and frame request. Immediate observers could see those intermediate states;
allocation or callback failure could leave a partial successor. Clip-only
changes were also omitted from `changedOut` and the frame obligation, while an
equal cutting-plane input could rebuild retained HUD geometry.

The operation now derives the viewport, camera scalar state, tracked-light
directions, clip planes, affordance geometry and orthographic depth cache in
local candidates. It prepares scalar-field notifications, a detached section
aid, an optional native camera-root replacement, the candidate LoD signature
and one capacity-level render request before an allocation-free commit. That
commit installs all controller, viewport and render-manager copies and every
derived node before notification. Every participant's normal notification
state is restored before the first callback; graph, field and frame observers
are then all attempted, with the first observer exception rethrown after the
complete successor remains visible. Same-projection navigation preserves
camera identity; changing projection replaces it. An exact input is a zero-
allocation, zero-notification no-op and reports `changedOut=FALSE`.

The libged navigation-gizmo integration exercises that enclosing boundary.
Same-projection synchronization publishes one fixed-orientation successor and
retains the camera, an equal synchronization publishes none, projection
change replaces the camera and publishes one successor, and mutation of the
retired generation publishes nothing. This exposed a consumer edge where a
same-orientation camera replacement retargeted the sensor but generic visual
equality suppressed the retained successor. Camera-generation replacement now
forces that single publication without weakening ordinary equal-input
coalescing.

The fresh focused selector passes **1,578 forced-allocation positions** across
initial viewport adoption, camera fields, clip-only enablement and both
projection directions with ordinary and automatic LoD: 1,536 retain the
complete predecessor and 42 expose the complete successor. Immediate root,
camera, affordance, representative field and frame observers see only complete
states; observer exceptions drain; every failed predecessor retries. The
saved CXX-SOURCE-066 controller fails the new complete-state assertion, and
its gizmo plugin fails the same-orientation projection-replacement assertion
against the repaired controller. The focused selector and dynamic plugin pass
under ASan, and the adjacent camera, viewport, lighting, host, LoD, public API,
GED and Qt checks pass. The exact full compact executable passes **162 direct
checks, 21 scenarios and 1,072 PASS rows in 363.92 seconds**. Exact evidence
is retained under `.build/obol-qualification/20260912-view-camera-sync`.

No formal relation changed. This boundary refines one external input into a
single owner-thread controller publication; it does not claim concurrent
cross-thread mutation. The public cutting-plane input setters, other
preview/manipulator/navigation input setters, Qt widget policy and the wider
lifecycle/resource gates remain separate classified work.

### CXX-SOURCE-068 — cutting-plane input publishes intent, clip and aid together

**Repaired and qualified for the public cutting-plane setters.** The preceding
setters changed private intent, the live `SoClipPlane`, the visible section aid
and the frame request in separate operations. Immediate field, node and path
observers could therefore see a mixed generation; allocation or callback
failure could interrupt the successor. Equal plane input rebuilt retained HUD
geometry, and invalid non-finite plane values entered the live controller.

Both setters now use one prepared publication. A detached scalar candidate,
optional section-aid replacement and capacity-level render request are fully
prepared before an allocation-free commit installs the clip fields, aid,
private intent and standing request. Field notification state is restored
before callbacks. Aid, scalar and frame notifications are all attempted, with
the first observer exception rethrown after the complete successor remains
visible. Exact input is an allocation-free, notification-free no-op. Zero-
normal and non-finite planes are rejected without mutation.

The focused selector covers plane changes while disabled and enabled, enable
and disable operations, each with ordinary and automatic LoD. It passes **711
forced-allocation positions**: 636 retain the complete predecessor and 75
expose the complete successor. Immediate root, clip-field, aid and frame
observers see only complete state; every failed predecessor retries; observer
exceptions drain; callback reentry retains the later successor. The exact
CXX-SOURCE-067 library rejects the complete-state assertion. The focused
selector passes under ASan, and the adjacent camera, affordance, environment,
host, GED and Qt checks pass. The exact full compact executable passes **163
direct checks, 21 scenarios and 1,073 PASS rows in 356.77 seconds**. Exact
evidence is retained under
`.build/obol-qualification/20260912-cutting-plane-input`.

No formal relation changed. This boundary refines one owner-thread input
family. The scene-light input setters, other preview/manipulator/navigation
setters, Qt widget policy and the wider lifecycle/resource gates remain
separate classified work.

### CXX-SOURCE-069 — scene-light input publishes values, nodes and policy together

**Repaired and qualified for scene-light content and enablement.** The
preceding content setter copied private input before rebuilding retained
children and did not publish a frame obligation. The enablement setter changed
private intent and each live `SoLight::on` field in sequence before requesting
a frame. GED supplied content and enablement through those two operations.
Immediate graph, field or endpoint observers could see a mixed generation;
allocation or callback failure could retain a partial successor. Equal input
still rebuilt nodes or requested another frame.

Both setters now enter one shared prepared publication. A copied value
snapshot, complete detached light children or retained-node scalar edits, and
one capacity-level render request are ready before an allocation-free commit.
The commit installs private state, all children or fields and the standing
request before callbacks. Enablement preserves light-node identity. Scalar
notification state is restored before callback delivery, all graph, field and
frame notifications are attempted, and callback reentry retains its later
successor. Exact input allocates and notifies nothing.

An additive public overload accepts contents and enablement together. The GED
lighting synchronizer uses it, removing its two-generation scene-light pair.
Existing callers can still publish either content or enablement independently
through the original complete operations.

The focused selector covers installation, replacement, clearing, nonempty and
empty enablement and combined changes, each under ordinary and automatic LoD.
It passes **443 forced-allocation positions**: 392 retain the complete
predecessor and 51 expose the complete successor. Immediate root, group, path,
light-field and frame observers see only complete state; failed predecessors
retry; observer exceptions drain; replacement and enablement reentry retain the
later successor. The saved CXX-SOURCE-068 library lacks the combined public
entry, and its source snapshot records the preceding sequential caller and
setter implementations. Focused ASan and adjacent scene-light, camera,
environment, host, view-store, GED and Qt checks pass. The exact full compact
executable passes **164 direct checks, 21 scenarios and 1,074 PASS rows in
374.97 seconds**. Exact evidence is
retained under `.build/obol-qualification/20260912-scene-light-input`.

No formal relation changed. This boundary refines the scene-light family. The
enclosing camera-rig/profile portion of `ged_view_lighting_sync()`, other
preview/manipulator/navigation setters, Qt widget policy and the wider
lifecycle/resource gates remain separate classified work.

### CXX-SOURCE-070 — complete lighting input publishes one controller generation

**Repaired and qualified for camera-rig and database-light input.** After
CXX-SOURCE-069, `ged_view_lighting_sync()` still published one
`bv_lighting_state` through four camera-rig setters followed by the combined
scene-light setter. Profile, offset, tracking, headlight policy, environment
fields, camera-light fields, database-light children and the frame obligation
could therefore expose several separately observable controller generations.

The controller now has one additive complete-state operation. It prepares
detached environment and camera-light scalar candidates, copied scene-light
intent, replacement children or retained-node enablement edits, and one
capacity render request before mutation. The allocation-free commit installs
all private inputs, scalar fields, children and the standing request before
callbacks. Notification state is restored before delivery, failures drain,
and callback reentry retains the later complete successor. Existing leaf
setters use the same owner. Direction canonicalization is stable across
unrelated setters, so exact inputs allocate and notify nothing; invalid
profile, zero direction and non-finite direction reject without mutation.

The fresh complete-state sweep covers ordinary and automatic LoD and passes
**109 forced-allocation positions**: 98 retain the complete predecessor and 11
expose the complete successor. The warmed full run covers 100 positions (98
predecessors, two successors). Immediate environment, camera-light, scene,
path and frame observers see only complete state; retry, callback failure and
reentry pass. Direct checks also cover each profile, offset, tracking and
headlight leaf plus exact repeats. The saved CXX-SOURCE-069 library lacks the
new complete-state entry. Focused ASan, adjacent publication, LoD, renderer,
host, GED and Qt checks pass. The exact full compact executable passes **165
direct checks, 21 scenarios and 1,075 PASS rows in 347.54 seconds**. Exact
evidence is retained under
`.build/obol-qualification/20260912-lighting-state-input`.

No formal relation changed. This boundary refines one classified controller
input family. Master shading publication, remaining preview/manipulator input
families, Qt widget policy and the wider lifecycle/resource gates remain
separate work.

### CXX-SOURCE-071 — master lighting publishes shading, rig and renderer epoch together

**Repaired and qualified for the master lighting input.**
`setLightingEnabled()` previously assigned the retained `SoLightModel` and
three camera-light enablement fields one at a time, then called renderer
invalidation and a second render request. Immediate field observers could see
a partial rig, and the frame callback received `renderer-performance` while
the standing request had already been renamed `lighting`.

The setter now prepares detached scalar candidates for the light model and
complete camera rig plus a reason-aware renderer invalidation. One
allocation-free commit publishes flat/PHONG shading, profile-dependent rig
enablement, cleared renderer timing/capacity evidence, automatic-LoD work and
one typed `lighting` frame before callbacks. Notification state is restored
before delivery, every field and frame notification is attempted after an
observer exception, and callback reentry retains the later master-lighting
successor. Exact input allocates nothing, preserves timing evidence and does
not request a frame. Unchanged scene-light nodes retain identity. The existing
one-argument private invalidation entry remains in place for its other users.

The focused Studio/MGED, ordinary/automatic sweep passes **120 forced-
allocation positions**: 104 retain the complete predecessor and 16 expose the
complete successor. Immediate light-model, camera-light and frame observers
see only complete state; retry, exact input, callback failure and reentry pass.
The saved CXX-SOURCE-070 library fails the same behavioral selector. Focused
ASan, eight adjacent publication selectors, ten adjacent CTests, public-header
compilation and static linking pass. The exact full compact executable passes
**166 direct checks, 21 scenarios and 1,076 PASS rows in 371.20 seconds**.
Exact evidence is retained under
`.build/obol-qualification/20260912-master-lighting-input`.

No formal relation changed. This boundary refines the renderer-appearance
family. Background, depth-test, transparency, antialiasing, clip-bound,
headlight color/intensity, depth-cue and software-wire inputs remain classified
work, along with the wider lifecycle/resource gates.

### CXX-SOURCE-072 — scalar appearance inputs publish retained and frame state together

**Repaired and qualified for the retained scalar appearance family.** The
background, depth-test, headlight color/intensity and depth-cue setters
previously wrote their private or Coin fields before requesting a frame.
Depth-test and depth-cue also notified a `renderer-performance` invalidation
before replacing it with a differently named request. Immediate field
observers could therefore see incomplete fields, stale renderer evidence or no
corresponding frame. Headlight color/intensity were unnecessarily classified
as LoD-capacity work, and an already matching fog mode prevented depth-cue from
repairing stale fog color or visibility.

The five setters now share one prepared scalar-appearance transaction. A
detached node candidate, complete field notifications and either one
presentation request or one reason-matched renderer invalidation are prepared
before mutation. Background private colors, their derived fog color, depth
test/write, light fields, depth-cue fields, renderer evidence and the frame
level commit before callbacks. Background and light color/intensity preserve
renderer timing and remain presentation-only. Depth-test and depth-cue clear
timing-derived capacity evidence and create automatic-LoD work only when that
policy is active. Exact input allocates nothing; non-finite colors and
intensity do no work; finite light values retain unit clamping; stale fog
derivatives are repaired; callback failures drain; and reentry retains the
later successor.

The seven-operation ordinary/automatic sweep passes **128 forced-allocation
positions**: 96 retain the complete predecessor and 32 expose the complete
successor. Immediate environment, depth-buffer, headlight and frame observers
see only complete state. The saved CXX-SOURCE-071 library fails the same
behavioral selector. Focused ASan, nine adjacent publication selectors, ten
adjacent CTests, public-header compilation and static linking pass. The exact
full compact executable passes **167 direct checks, 21 scenarios and 1,077
PASS rows in 338.31 seconds**. Exact evidence is retained under
`.build/obol-qualification/20260912-scalar-appearance-input`.

No formal relation changed. This boundary refines the renderer-appearance
family. Transparency, antialiasing, clip bounds and software-wire mode remain
classified work, along with the wider lifecycle/resource gates.

### CXX-SOURCE-073 — render-action appearance has one redraw and renderer epoch

**Repaired and qualified for transparency and antialiasing.** Both setters
previously wrote the controller and onscreen/cached-offscreen render actions,
then published a `renderer-performance` invalidation and a second differently
named request. The antialiasing path also called
`SoRenderManager::setAntialiasing()`, which scheduled a toolkit redraw outside
the controller's frame owner. A frame observer could see the changed actions
before the final typed request, and exact private state concealed an externally
stale action.

The setters now prepare one reason-matched renderer invalidation before their
commit. Controller intent and both available render actions then update
through Coin's nonallocating transparency, smoothing and pass-count setters;
renderer timing/capacity evidence and one typed frame commit before the
callback. Antialiasing writes the render action directly with the named
single-pass policy, leaving the controller request as the sole redraw owner.
Exact checks include every available action and repair stale backend state.
Callback failure retains the committed state, and callback reentry retains the
later successor.

The four change/repair operations in ordinary and automatic LoD pass **eight
allocation-free operation states**. The saved CXX-SOURCE-072 library fails the
same complete-state selector. Focused ASan, ten adjacent publication selectors,
ten adjacent CTests, public-header compilation and static linking pass. The
exact full compact executable passes **168 direct checks, 21 scenarios and
1,078 PASS rows in 356.16 seconds**. Exact evidence is retained under
`.build/obol-qualification/20260912-render-action-appearance-input`.

No formal relation changed. This boundary refines the renderer-appearance
family. Clip bounds and software-wire mode remain classified work, along with
the wider lifecycle/resource gates.

### CXX-SOURCE-074 — clip bounds publish derived planes and renderer epoch together

**Repaired and qualified for camera-relative clip bounds.** The setter
previously changed only its two private bounds, notified a separately named
renderer invalidation and then requested a `clip-bounds` frame. The retained
minimum and maximum `SoClipPlane` fields were not recomputed until a later
camera synchronization, so the requested frame traversed the preceding planes.
An equal private value also concealed either stale retained plane.

Camera synchronization and the setter now share one finite camera-relative
plane derivation. The setter derives both world-space planes from the retained
camera and last synchronized horizontal size, prepares detached scalar
candidates and one reason-matched renderer invalidation, then commits the
changed planes, private bounds, cleared timing/capacity evidence and typed frame
before callbacks. Active and disabled clipping retain their enablement; a
one-sided bound change preserves the other plane; exact input allocates and
publishes nothing; stale planes are repaired; invalid or unrepresentable input
does no work; callback failures drain; and callback reentry retains the later
successor. The multi-node scalar helper now lives with the scalar publication
primitive and remains shared by camera synchronization.

The five active/disabled/change/repair operations in ordinary and automatic
LoD pass **104 forced-allocation positions**: 80 retain the complete predecessor
and 24 expose the complete successor. Immediate clip-node, clip-field and frame
observers see only complete state. The saved CXX-SOURCE-073 library fails the
same behavioral selector. Focused ASan, twelve adjacent publication selectors,
ten adjacent CTests, public-header compilation and static linking pass. The
exact full compact executable passes **169 direct checks, 21 scenarios and
1,079 PASS rows in 323.52 seconds**. Exact evidence is retained under
`.build/obol-qualification/20260912-clip-bounds-input`.

No formal relation changed. This boundary refines the renderer-appearance
family. Software-wire mode remains classified work, along with the wider
lifecycle/resource gates.

### CXX-SOURCE-075 — software-wire policy has one complete renderer generation

**Repaired and qualified for software-wire input publication.** The setter
previously kept three runtime copies: controller intent, the LoD wrapper and
the compact render batch. Its equality check consulted only controller intent,
so either render path could remain stale indefinitely. A change also published
a separately named `renderer-performance` callback before replacing the
standing request with `software-wire-mode`; that observer could traverse a
generation whose policy and diagnostic frame owner disagreed.

The setter now normalizes its input, checks every available runtime copy and
prepares one reason-matched renderer invalidation before mutation. Controller,
LoD-wrapper and compact-batch policy then commit with cleared renderer timing
and capacity evidence and one typed frame before the callback. The runtime
policy writes are private scalar operations, so no additional transaction
type is needed. Exact input allocates and publishes nothing; equal controller
intent repairs a stale wrapper or batch; a controller without an LoD wrapper
still updates its compact render path; callback failure retains the committed
generation; and callback reentry retains the later successor.

Seven mode/change/repair/render-path operations in ordinary and automatic LoD
pass **56 forced-allocation positions**: 28 retain the complete predecessor and
28 expose the complete successor. The saved CXX-SOURCE-074 library fails the
same complete-state selector. Focused ASan, thirteen adjacent publication
selectors, ten adjacent CTests, public-header compilation and static linking
pass. The source suite passes **275 PASS rows in 139.29 seconds**. The exact
full compact executable passes **170 direct checks, 21 scenarios and 1,080
PASS rows in 320.72 seconds**. Exact evidence is retained under
`.build/obol-qualification/20260912-software-wire-input`.

No formal relation changed. This boundary closes the classified
renderer-appearance family. Preview, manipulator, navigation, traversal and
lifecycle inputs remain separate classified work, together with the wider
resource and capability gates.

### CXX-SOURCE-076 — controller edit-preview input publishes one retained generation

**Repaired and qualified for the direct controller edit-preview API.** An
existing preview previously received its identity, edit intent and revisions
before line geometry preparation. Allocation failure could therefore retain
new input fields with the preceding realized child and status. Successful
replacement and insertion notified the retained graph before the controller
created its frame obligation; removal similarly notified the root before
requesting its frame.

Replacement now constructs one complete detached preview candidate before
touching the live node. Prepared scalar values and a prepared child-list edit
then preserve the existing preview identity while committing identity, intent,
requested and realized revisions, status, geometry and the capacity-classified
frame together. The prepared field notifications identify the preview's own
stale-marking sensors as already handled, while external field, node and path
observers still run against the complete generation. Initial insertion uses a
prepared root order, and removal uses prepared path/root retirement. All graph,
field and frame notifications drain while retaining the first exception.
Invalid input and missing removal do no work, and callback reentry retains the
later successor.

Insert, replace and remove in ordinary and automatic LoD pass **364 stable
forced-allocation positions**: 350 retain the complete predecessor and 14
expose the complete successor. The selector self-warms one-time Coin
notification state before measurement, so standalone and full-suite counts
match. The saved CXX-SOURCE-075 library fails the same complete-state selector.
Focused ASan, fourteen adjacent publication selectors, ten adjacent CTests,
public-header compilation and static linking pass. The source suite passes
**276 PASS rows in 140.46 seconds**. The exact full compact executable passes
**171 direct checks, 21 scenarios and 1,081 PASS rows in 320.73 seconds**.
Exact evidence is retained under
`.build/obol-qualification/20260912-controller-edit-preview-input`.

No formal relation changed. This closes the direct edit-preview controller
input without reopening the qualified node-local preview or feature-store edit
transaction. Direct line-layer and HUD-label controller inputs, remaining
traversal/lifecycle writers and the wider resource/capability gates remain
separate work.

### CXX-SOURCE-077 — controller overlay inputs publish one complete graph and frame

**Repaired and qualified for the direct line-layer and HUD-label controller
APIs.** Line-layer replacement previously built most geometry off-graph, but
held the candidate without reference ownership and then called the live root's
add, replace or remove operation before requesting a frame. An allocation or
throwing graph observer could leak the candidate or leave a changed graph with
no matching frame obligation. HUD-label replacement was more direct: it wrote
six fields into the existing live node and only then rebuilt allocating HUD
geometry. Failure could retain new input fields with the old label, and all
three insert/replace/remove paths exposed the graph before the frame request.

Both APIs now prepare a complete detached candidate under reference ownership.
Line-layer insertion and replacement prepare the whole root order and preserve
the established replacement-node semantics; removal uses prepared root/path
retirement. HUD insertion uses the same root operation. Existing HUD labels
retain their node identity through one prepared scalar-field and child-list
publication, with the cached `SoHUDLabel` pointer committed in the same quiet
phase. HUD removal also uses prepared retirement. Each operation prepares and
commits one capacity-classified frame before graph or field observation. A
shared notification drain retains the first exception while still attempting
the graph, field and frame callbacks. Invalid or missing input does no work,
and callback reentry retains the later successor.

Line-layer insert, replace and remove in ordinary and automatic LoD pass **348
stable forced-allocation positions**: 322 retain the complete predecessor and
26 expose the complete successor. The corresponding HUD-label operations pass
**538 stable positions**: 490 predecessors and 48 successors. The selectors
self-warm one-time Coin state so standalone and suite counts match. The saved
CXX-SOURCE-076 library fails both complete-state selectors. Focused ASan for
both selectors, fourteen adjacent publication selectors, ten adjacent CTests,
public-header compilation and static linking pass. The source suite passes
**278 PASS rows in 140.33 seconds**. The exact full compact executable passes
**173 direct checks, 21 scenarios and 1,083 PASS rows in 320.43 seconds**.
Exact evidence is retained under
`.build/obol-qualification/20260912-controller-overlay-inputs`.

No formal relation changed. This closes the direct controller overlay input
family without reopening the qualified node-local HUD/line rebuilds or feature-
store overlay transactions. Remaining manipulator/navigation field inputs,
traversal/lifecycle writers and the wider resource/capability gates remain
separate work.

### CXX-SOURCE-078 — grid view input and derived geometry publish one generation

**Repaired and qualified for `bobol_grid_configure_from_view()`.** The helper
previously wrote visibility, snapping, anchor, spacing, division, view,
viewport and color fields into the live grid one at a time, then invoked an
allocating geometry rebuild. Immediate field or node observers could see an
arbitrary prefix of the new input paired with preceding derived spacing,
segment counts and HUD geometry. Allocation failure during the rebuild retained
that mixed live state. The context overload inherited the same boundary through
delegation.

Configuration now copies the target's retained scalar policy into a detached,
reference-owned grid candidate, applies the complete sanitized view input and
performs the already-qualified grid rebuild there. The target child replacement
and every changed scalar field then prepare before mutation. One quiet commit
preserves the target grid identity while publishing input fields, derived
fields and complete visible or empty geometry together. Child and field
notifications drain while retaining the first exception. Invalid input does no
work, and callback reentry retains the later complete configuration.

Visible, hidden and snap-only configurations pass **9,422 stable forced-
allocation positions**: 229 retain the complete predecessor and 9,193 retain
the complete successor. The large successor count is the deliberate sweep of
allocation failures in every post-commit public-field notification; those
failures cannot roll back or interrupt the remaining callbacks. The selector
self-warms one-time Coin state so standalone and suite counts match. The saved
CXX-SOURCE-077 library fails the same complete-state selector. Focused ASan,
seven adjacent publication selectors, ten adjacent CTests, public-header
compilation and static linking pass. The source suite passes **279 PASS rows in
145.02 seconds**. The exact full compact executable passes **174 direct checks,
21 scenarios and 1,084 PASS rows in 323.87 seconds**. Exact evidence is retained
under `.build/obol-qualification/20260912-grid-input-publication`.

No formal relation changed. This closes the enclosing grid view-input
composition without reopening the node-local grid rebuild. Axes/ADC input
composition, image-source mutation, navigation camera/hover input, remaining
traversal/lifecycle writers and the wider resource/capability gates remain
separate work.

### CXX-SOURCE-079 — axes and ADC view inputs publish one retained generation

**Repaired and qualified.** `SoBRLAxes` and `SoBRLADC` had atomic geometry
rebuilds, but retained callers still had no operation which composed their
public input fields with that derived geometry. The only available sequence
wrote the attached fields one at a time and then called `rebuildGeometry()`.
Immediate observers and allocation failure could therefore retain a prefix of
the next input with the preceding shapes.

`bobol_axes_configure_from_view()` and `bobol_adc_configure_from_view()` now
copy retained node policy into a reference-owned detached candidate, apply the
supported `bv_axes_state` or `bv_adc_state` mapping and invoke the already
qualified builder off graph. A shared scalar-field/child-list publisher then
prepares and commits the complete candidate while preserving the attached node
identity. Child and field callbacks drain while retaining the first exception;
invalid pointers do no work, and callback reentry retains the later complete
generation. Grid view configuration now uses the same private publication
primitive. The Qt fallback and GED ADC producer use these canonical mappings
when constructing their detached nodes.

Axes and ADC, each visible and hidden, pass **222 stable forced-allocation
positions**: 200 retain the complete predecessor and 22 retain the complete
successor. The exact CXX-SOURCE-078 library fails because an immediate
observer sees a partial input generation. Focused ASan, eight adjacent
publication selectors, ten adjacent CTests, public-header compilation, both
new exported symbols and static linking pass. The source suite passes **95
direct checks, 21 scenarios and 280 PASS rows in 144.70 seconds**. The exact
full compact executable passes **175 direct checks, 21 scenarios and 1,085
PASS rows in 323.13 seconds**. Exact evidence is retained under
`.build/obol-qualification/20260912-overlay-input-publication`.

No formal relation changed. This closes composition for the fields represented
by the current nodes: axes draw/position/size and ADC draw/center/angle/
distance/colors/width, while preserving ADC crosshair/tick sizing. Rich axes
labels, ticks and style remain a separate capability question rather than a
partially published implementation. Image-source mutation, navigation camera/
hover input, remaining traversal/lifecycle writers and the wider resource/
capability gates remain separate work.

### CXX-SOURCE-080 — navigation gizmo inputs publish one retained generation

**Repaired and qualified.** `SoBRLNavigationGizmo::setHoverPart()`,
`setActivePart()` and `setCamera()` changed live input or camera ownership
before an allocating geometry rebuild. Camera replacement also detached the
old sensor before the successor sensor and HUD were known to be buildable.
Allocation failure or an immediate observer could therefore retain or see a
new input paired with the preceding HUD, or lose the preceding camera
attachment. Camera-snapshot and traversal-driven orientation updates had the
same geometry-publication gap.

A private publication object now copies the gizmo's scalar policy into a
reference-owned detached candidate, stages the requested orientation or
interaction input and builds its complete HUD off graph. Prepared scalar
fields and a prepared child replacement then commit with the private rotation,
HUD and anchor pointers before callbacks. Camera changes additionally prepare
a freshly attached sensor and retained camera before changing target
ownership. The unrealized path changes only private camera state after all
allocation succeeds. Exact no-ops allocate and notify zero times; callback
failures drain while retaining the first exception, and callback reentry keeps
the later complete successor.

Hover, active, realized camera add/replace/remove and unrealized camera
add/replace/remove cover **5,172 stable forced-allocation positions**: 5,162
retain the complete predecessor and 10 retain the complete successor. A fixed
snapshot on an unrealized gizmo is also allocation-free. The exact
CXX-SOURCE-079 library fails because an immediate observer sees a partial input
generation. Focused ASan, nine adjacent publication selectors, ten adjacent
CTests, public-header compilation, the four public method symbols and static
linking pass. The source suite passes **96 direct checks, 21 scenarios and 281
PASS rows in 169.71 seconds**. The exact full compact executable passes **176
direct checks, 21 scenarios and 1,086 PASS rows in 374.37 seconds**. Exact
evidence is retained under
`.build/obol-qualification/20260913-navigation-gizmo-input-publication`.

No formal relation changed. This closes the supported navigation gizmo method
boundary and its fixed-snapshot producer. Direct writes to public style fields,
renderer-derived anchor placement, image-source mutation, remaining
traversal/lifecycle writers and the wider resource/capability gates remain
separately classified work.

### CXX-SOURCE-081 — image-source inputs publish one retained generation

**Repaired and qualified.** `SoBRLImageSource::setStream()`, `setImage()`,
`clearSource()` and `refreshFromStream()` previously released the current
subscriber or owned stream and then wrote source metadata one field at a time.
Allocation failure could reject a replacement after destroying its predecessor,
and immediate field or node observers could see mixed source kind, dimensions,
format, revisions, dirty rectangle and connection state. Reattaching the owned
stream returned by `getStream()` also destroyed that stream before attempting to
reuse its dangling pointer.

The source now prepares a detached scalar candidate and a subscription owner
which contains the prospective stream, subscriber, source kind and dirty-callback
staging state. Subscription, image conversion, stream inspection and all scalar
storage complete before one quiet commit installs the subscription, generation
counters and public fields. The predecessor subscription is deactivated before
retirement, and notification drains while retaining the first callback exception.
Exact stream reattachment, repeated clear/failure states and settled refreshes
allocate and notify zero times. Dirty callbacks attached to a candidate update
only that candidate until commit. Attachment also clears the prior stream's dirty
metadata before applying the new stream information.

Stream add, dirty/clean stream replacement, static-image replacement, clear,
invalid stream, invalid image and refresh cover **174 stable forced-allocation
positions**: 143 retain the complete predecessor and 31 retain the complete
successor. Owned-stream self-reattachment is also an exact no-op. Node, kind,
status, source-revision and data-revision observers see only a complete state;
throwing callbacks drain and reentry retains the later borrowed stream. The exact
CXX-SOURCE-080 library fails because an immediate observer sees a partial input
generation. A before/current compilation probe records the same 928-byte
`SoBRLImageSource` size, so the private ownership change does not alter this
class's C++ object layout.

Focused ASan, nine adjacent publication selectors, ten adjacent CTests,
public-header compilation and static linking pass. The source suite passes **97
direct checks, 21 scenarios and 282 PASS rows in 171.80 seconds**. The exact
full compact executable passes **177 direct checks, 21 scenarios and 1,087 PASS
rows in 373.51 seconds**. Exact evidence is retained under
`.build/obol-qualification/20260913-image-source-input-publication`.

No formal relation changed. This closes the supported image-source method
boundary. Direct writes to public source fields, composition with a display
consumer, remaining traversal/lifecycle writers and the wider resource/
capability gates remain separately classified work.

### CXX-SOURCE-082 — framebuffer cursor input publishes one presentation generation

**Repaired and qualified.** `BObolWindowHost::setFramebufferCursor()`,
`setFramebufferScreenCursor()` and `setFramebufferCursorShape()` previously
wrote cursor visibility, image position and shape directly in sequence, then
requested presentation. Immediate field or node observers could see mixed
cursor values with no standing frame request, and an observer exception could
interrupt the remaining fields or suppress presentation. Equal input still
requested another frame.

All three methods now use one fixed three-field publisher. The presentation
request and its owned diagnostic are prepared first; one allocation-free quiet
commit installs all cursor fields and the standing request before callbacks.
Notification state is restored before delivery, field and endpoint callbacks
drain while retaining the first exception, and callback reentry retains its
later cursor state. Exact input allocates and notifies zero times. Preparation
is limited to the three changed fields and one render-request value.

Cursor enable, screen-cursor disable, custom-shape installation and shape reset
in ordinary and automatic LoD cover **44 stable forced-allocation positions**:
4 retain the complete predecessor and 40 expose the complete successor.
Immediate node, visibility, position, shape and endpoint observers see only a
complete state; failure retry, callback exceptions and reentry pass. The exact
CXX-SOURCE-081 library fails because the first observer sees a partial cursor
generation. Focused ASan, ten adjacent publication selectors, ten adjacent
CTests, the three existing public method symbols, public-header compilation and
static linking pass. The source suite passes **98 direct checks, 21 scenarios
and 283 PASS rows in 155.23 seconds**. The exact full compact executable passes
**178 direct checks, 21 scenarios and 1,088 PASS rows in 335.61 seconds**. Exact
evidence is retained under
`.build/obol-qualification/20260913-framebuffer-cursor-input-publication`.

No formal relation changed. This closes the supported framebuffer cursor-state
methods. Framebuffer reset, viewport/view inputs, source-to-display refresh,
custom bitmap consumption, direct viewport field writes, endpoint cleanup and
the wider resource/capability gates remain separately classified work.

### CXX-SOURCE-083 — framebuffer view input publishes one retained presentation generation

**Repaired and qualified.** `BObolWindowHost::setFramebufferView()` previously
wrote source center and zoom as separate live fields, rebuilt geometry through
an image-source refresh, and requested presentation afterward. Immediate
observers could see mixed transform fields, old geometry and no standing frame
request. Allocation failure could retain changed transform fields without the
matching child graph, while a callback exception could interrupt the rest of
the sequence. A view-only operation also consumed a pending stream generation
as an unrelated source publication.

The window host now stages the transform on the existing detached viewport
candidate. A retained-payload path copies the viewport's accepted texture bytes
and realized revisions, then the shared internal viewport geometry builder
derives texture coordinates and the complete HUD off target. Viewport scalar
fields, child graph, cached texture/face pointers and the presentation request
commit before callbacks. The pending stream generation and all source fields
remain unchanged for flush to publish. The staged builder preserves the
existing public viewport layout and exported member surface.

Nested view publication also exposed a deferred reset in the in-memory texture
helper: `SoTexture2::setImageData()` rewrote an already empty filename, whose
sensor could run after reentrant pixel installation and clear the new image.
Fresh in-memory textures now set their public image field directly, avoiding
that deferred filename event while retaining copied pixel ownership.

Combined pan/zoom, pan-only and zoom-only operations in ordinary and automatic
LoD cover **1,322 stable forced-allocation positions**: 1,308 retain the
complete predecessor and 14 expose the complete successor. Node, field, path,
image-source and endpoint observers see one complete viewport/request state;
retained pixels, texture coordinates, cache pointers and revisions are checked
through failure, callback exception and reentry. Exact input allocates, changes
and notifies zero times. The exact CXX-SOURCE-082 library fails the focused
complete-state selector.

Focused ASan, ten adjacent publication selectors, ten adjacent CTests, the
unchanged public method symbol, absence of new exported helpers, the public
umbrella and static linking pass. The source suite passes **99 direct checks,
21 scenarios and 284 PASS rows in 154.94 seconds**. The exact full compact
executable passes **179 direct checks, 21 scenarios and 1,089 PASS rows in
337.67 seconds**. OSMesa remains `097aac27...`. Exact evidence is retained
under `.build/obol-qualification/20260913-framebuffer-view-input-publication`.

No formal relation changed. This closes the framebuffer view method boundary.
Framebuffer reset and viewport input, source-to-display flush publication,
direct viewport field writes, endpoint cleanup and the wider resource/capability
gates remain separately classified work.

### CXX-SOURCE-084 — framebuffer reset publishes one retained presentation generation

**Repaired and qualified.** `BObolWindowHost::resetFramebuffer()` previously
wrote center, zoom and cursor visibility as separate live fields, refreshed the
image source, rebuilt the viewport, and requested presentation afterward.
Observers and allocation failure could retain mixed reset fields, old geometry
or no standing frame request. The view-state operation also consumed an
unrelated pending stream generation before the display refresh boundary could
publish it.

Reset now performs an exact no-op check before allocation, edits the detached
viewport candidate, and uses the retained-payload transaction introduced by
CXX-SOURCE-083. Center, zoom, cursor visibility, complete HUD geometry, cached
texture/face pointers and the presentation request commit before callbacks.
Accepted texture bytes and realized revisions survive, while pending pixels
and every image-source field remain owned by `flushFramebuffer()`. Reset and
view share the commit/notification protocol rather than duplicating it.

Ordinary and automatic LoD cover **442 stable forced-allocation positions**:
436 retain the complete predecessor and 6 expose the complete successor.
Node, field, path, image-source and endpoint observers see only a complete
state. Exact repeated reset allocates, changes and notifies zero times;
failure retry, callback drain and reentry into a later view pass. The exact
CXX-SOURCE-083 library fails on a partial reset observation. Focused ASan, ten
adjacent publication selectors and ten adjacent CTests pass. The source suite
passes **100 direct checks, 21 scenarios and 285 PASS rows in 155.17 seconds**.
The exact full compact executable passes **180 direct checks, 21 scenarios and
1,090 PASS rows in 337.34 seconds**. Public symbols, the public umbrella and
static linking pass; OSMesa remains `097aac27...`. Exact evidence is retained under
`.build/obol-qualification/20260913-framebuffer-reset-input-publication`.

No formal relation changed. This closes the framebuffer reset method boundary.
Framebuffer viewport input, source-to-display flush publication, direct
viewport field writes, endpoint cleanup and the wider resource/capability
gates remain separately classified work.

### CXX-SOURCE-085 — framebuffer viewport input composes placement, controller size and frame obligation

**Repaired and qualified.** `BObolWindowHost::setFramebufferViewport()`
previously published controller viewport size first, then wrote framebuffer
position and size separately, refreshed source/geometry, and finally requested
presentation. The controller could wake its endpoint while the framebuffer HUD
still had the old placement. Observers and allocation failure could also retain
a split controller/viewport generation or consume pending pixels outside flush.

The host now normalizes the requested rectangle, checks the complete target
before allocation, and prepares retained HUD geometry before beginning the
controller size transition. Position, size, child graph, cached texture/face
pointers, all three controller viewport copies, LoD view state and one
`fb-viewport` request commit before callbacks. A size change uses the
controller's capacity-aware viewport transaction; a position-only change uses
a presentation request. Displayed texture bytes and revisions remain accepted,
while pending source data remains pending for `flushFramebuffer()`.

Combined move/resize, move-only, resize-only and invalid-extent fallback in
ordinary and automatic LoD cover **1,776 stable forced-allocation positions**:
1,756 retain the complete predecessor and 20 expose the complete successor.
Immediate viewport node/field/path, image-source and endpoint observers see
only the composed state, including the three controller region copies and LoD
revision/progressive-work result. Exact input allocates, changes and notifies
zero times; retry, callback drain and reentry pass. The exact CXX-SOURCE-084
library fails on a partial viewport observation. Focused ASan, ten adjacent
publication selectors and ten adjacent CTests pass. The source suite passes
**101 direct checks, 21 scenarios and 286 PASS rows in 157.54 seconds**. The
exact full compact executable passes **181 direct checks, 21 scenarios and
1,091 PASS rows in 338.82 seconds**. Public symbols, the public umbrella and
static linking pass; OSMesa remains `097aac27...`. Exact evidence is retained
under `.build/obol-qualification/20260913-framebuffer-viewport-input-publication`.

No formal relation changed. This closes the framebuffer viewport method
boundary. Source-to-display flush publication, direct viewport field writes,
endpoint cleanup and the wider resource/capability gates remain separately
classified work.

### CXX-SOURCE-086 — framebuffer flush publishes one source-to-display generation

**Repaired and qualified.** `BObolWindowHost::flushFramebuffer()` previously
called `SoBRLImageSource::refreshFromStream()`, then
`SoBRLViewportImage::syncFromSource()`, then requested a frame. Source
observers could therefore run before the new texture, realized revisions and
frame obligation existed. Viewport observers could run before the request, and
allocation failure could retain an accepted source generation with its old
display.

The existing private image-source publisher is now reusable by enclosing
transactions. Flush first performs a nonallocating source/viewport currency
check. A changed call prepares source metadata, reads and validates one stable
stream payload, builds the detached viewport successor and prepares one
`fb-flush` presentation request. Source fields and realized generation,
viewport fields, child graph, texture/face cache pointers, texture bytes and
the standing request commit quietly. Both nodes restore their notification
gates before any callback, then source, viewport/path and endpoint callbacks
are all attempted while retaining the first exception. A post-read stream-info
check rejects a payload whose metadata changed during preparation.

Dirty-pixel flush in ordinary and automatic LoD covers **486 stable
forced-allocation positions**: 476 retain the complete pending predecessor and
10 expose the complete displayed successor. Immediate source,
viewport-field/path and endpoint observers see only the composed generation.
Exact settled and invalid flushes allocate, change and notify zero times;
retry, callback drain and callback-time later-pixel flush retain the last
successor. The exact CXX-SOURCE-085 library fails on its first partial source
observation. Focused ASan, eleven adjacent publication selectors and ten
adjacent CTests pass. The source suite passes **102 direct checks, 21 scenarios
and 287 PASS rows in 156.55 seconds**. The exact full compact executable passes
**182 direct checks, 21 scenarios and 1,092 PASS rows in 338.58 seconds**.
Separate compilation of the three touched implementation units, public
symbols, the public umbrella and static linking pass; OSMesa remains
`097aac27...`. Exact evidence is retained under
`.build/obol-qualification/20260913-framebuffer-flush-input-publication`.

No formal relation changed. This closes the framebuffer source-to-display
flush boundary. Direct viewport field writes, endpoint cleanup and the wider
resource/capability gates remain separately classified work.

### CXX-SOURCE-087 — retained RT image publication no longer exposes an empty or split generation

**Repaired and qualified.** The retained RT endpoint previously attached its
viewport before the first source payload and geometry existed. Restart and
worker presentation then wrote the private stream, published
`SoBRLImageSource::refreshFromStream()` and published
`SoBRLViewportImage::syncFromSource()` as separate operations. A root observer
could see an empty first viewport, while a source observer could see accepted
metadata beside the preceding texture and realized revisions. A callback-time
nested restart could also start the newer worker and then let the older stack
frame overwrite its joinable thread, terminating the process.

RT resources are now constructed under scoped ownership and remain detached
until a complete image is prepared. One internal publisher writes or reuses an
exact staged stream generation, prepares source metadata, validates one stable
payload, builds the viewport successor, prepares first root attachment and an
`rt-restart` capacity request, and commits every live participant before
callbacks. Source and viewport notification gates restore before the first
observer, and all node/path/root/request callbacks drain while retaining the
first exception. Failed restart staging retains its exact pixels and transport
generation for retry; a different stream generation invalidates that reuse.
Worker presentation restores an uncommitted dequeued frame, while a committed
callback failure records its output planes. A restart ownership generation
prevents an older reentrant invocation from starting over its successor.
Numeric and quality RT policy writes now preserve their previous value when a
restart is rejected, and exact repeated values do no work.

Ordinary and automatic LoD restart cover **538 stable forced-allocation
positions**: 474 retain the complete predecessor and 64 expose the complete
successor. Two separately observed first activations expose only a fully
realized viewport/root attachment. Immediate source, viewport-field and old-path
observers see one generation; exact policy input allocates, changes and notifies
zero times; retry, callback drain and callback-time nested restart retain the
last successor. The exact CXX-SOURCE-086 library fails when first attachment
exposes an empty viewport. Focused ASan, ten adjacent publication selectors and
ten adjacent CTests pass. The source suite passes **103 direct checks, 21
scenarios and 288 PASS rows in 165.80 seconds**. The exact full compact
executable passes **183 direct checks, 21 scenarios and 1,093 PASS rows in
351.12 seconds**. Independent `display_endpoint.cpp` compilation, public
symbols, the public umbrella and static linking pass; OSMesa remains
`097aac27...`. Exact evidence is retained under
`.build/obol-qualification/20260913-rt-image-publication`.

No formal relation changed. This closes retained RT first-image and restart
source-to-display publication. Renderer-engine selection, endpoint teardown,
direct public viewport fields and the wider resource/capability gates remain
separately classified work.

### CXX-SOURCE-088 — renderer selection publishes one complete endpoint generation

**Repaired and qualified.** Engine selection previously destroyed an outgoing
RT presentation and assigned the new public policy before renderer-performance
invalidation, graphical render-manager synchronization or the target request
had finished preparing. An allocation failure could therefore leave the new
engine beside the preceding root, or leave the old policy after its RT state
had already been destroyed. Enabling graphics from `NONE` also let Coin throw
while rebinding the render manager after the public engine had changed.

Selection now prepares controller invalidation and any outgoing root
replacement before mutation. RT activation composes with the retained image
publisher: source, payload, viewport, root and request prepare first, then a
precommit hook publishes engine policy, graphical eligibility and performance
invalidation before the image/root/request commit. Leaving RT commits the
target policy and prepared root removal together. A retiring RT state remains
owned through every source, viewport, root and request callback, so a nested
selection can install a distinct successor without an older notification stack
freeing it. Coin render-manager synchronization and non-graphical request
clearing retain explicit retry obligations when their postcommit effects
throw. A graphical frame notification remains deferred until synchronization
succeeds. Repeating the selected policy settles those obligations; a settled
repeat allocates and notifies zero times. A committed image callback failure no
longer prevents the selected RT worker from launching.

Twelve transitions in both ordinary and automatic LoD cover **2,887 stable
forced-allocation positions**: 1,932 retain the complete preceding engine and
955 expose the complete successor. Coverage includes `AUTO`, `HW`, `SW`, `RT`,
`NONE` and `DIAGNOSTIC`, both directions across graphical suppression, RT
activation/removal, settled no-ops, rejected and committed retry, notification
drain, and callback-time RT-to-SW and SW-to-RT replacement. Root and path
observers see only the old or complete target generation. The exact
CXX-SOURCE-087 library exposes a partial generation. Focused ASan and the
adjacent endpoint/host checks pass. The source suite passes **104 direct checks,
21 scenarios and 289 PASS rows in 168.70 seconds**. The exact full compact
executable passes **184 direct checks, 21 scenarios and 1,094 PASS rows in
349.59 seconds**. Independent `display_endpoint.cpp` compilation, public
symbols, the public umbrella and static linking pass; OSMesa remains
`097aac27...`. Exact evidence is retained under
`.build/obol-qualification/20260913-renderer-engine-publication`.

No formal relation changed. This closes endpoint renderer-policy transition
publication. Endpoint teardown, direct public viewport fields and the wider
resource/capability gates remain separately classified work.

### CXX-SOURCE-089 — endpoint destruction is terminal, no-throw and callback-safe

**Repaired and qualified.** Endpoint destruction previously reused live,
fallible publication paths. A throwing presentation observer or allocation
failure during automatic-LoD cancellation could skip later owners, allow an
exception to cross the C destroy entry point or C++ destructor, and leave
render, capacity or progressive work live. Destruction from a synchronous
endpoint callback could also delete the endpoint while its public operation
still held stack state. RT presentation teardown could detach an incomplete
borrowed graph or leave stream ownership split between the endpoint and image
source.

Terminal cleanup now has a dedicated no-throw path. It attempts each
presentation reset independently, makes the active automatic-LoD generation
stale before fallible service cancellation, and then clears every render,
capacity and progressive owner under the controller lock. Endpoint operations
hold an owner-thread scope; callback-time destruction marks the endpoint
terminal immediately and defers deletion until the outer operation drains.
All subsequent endpoint mutation rejects that terminal object. Host teardown
uses a private no-throw detach primitive, preserves an unbound external
context, and severs borrowed controller state even if close fails. The image
source owns an adopted RT stream, and RT retirement restores a complete safe
root for a borrowed controller before releasing its source, viewport and
worker state. The controller, rather than a particular endpoint, owns any
pending render-manager synchronization retry, so a replacement endpoint can
settle it.

The focused matrix covers **30 owned/borrowed engine allocation positions**
and **9 self-contained borrowed-root fallbacks** across `AUTO`, `NONE`,
`DIAGNOSTIC` and `RT`. Throwing observers, active automatic-LoD generations,
source-owned RT streams, callback-time destroy with reentrant mutation, engine
activation callbacks and render-manager retry settlement pass. The exact
CXX-SOURCE-088 library terminates with status 42 when endpoint destruction
propagates terminal cleanup failure. Focused ASan, ten adjacent publication
selectors and ten adjacent CTests pass. The source suite passes **105 direct
checks, 21 scenarios and 290 PASS rows in 170.39 seconds**. The exact full
compact executable passes **185 direct checks, 21 scenarios and 1,095 PASS
rows in 384.03 seconds**. Four implementation units compile independently;
public symbols, the public umbrella and static linking pass. Exact evidence is
retained under
`.build/obol-qualification/20260913-endpoint-destruction`.

No formal relation changed. This closes endpoint terminal teardown. Direct
public viewport fields and the wider resource/capability gates remain
separately classified work.

### CXX-CONFIG-002 — scene configuration effects can be interrupted

**Repaired for scene draw/representation/rename.** The original controller sequenced
source publication and scene effects. A draw observer could interrupt the paired
representation change and frame revision; representation similarly lost its frame.
Rename unindexed the source before retargeting, so an observer could leave a live
source absent from a valid index. Conflict removal and repository retirement could
also run before source preparation succeeded.

Draw and representation now use one prepared source publication with a scene commit
hook. Compatible external wire/mesh geometry remains realized across a draw request's
paired representation update; an explicit representation request retains its separate
invalidation semantics. Frame revisions publish before observers, including the
existing coalescing behavior inside a scene mutation batch.

Rename prepares only affected index entries, path memberships, child/path removal
and repository retirement. Detached hash nodes and required capacity are ready before
commit. The source, graph, indexes, cache ownership and structural/frame revisions
commit without allocation or notification; source and removed-path observers are then
attempted while retaining the first exception. Removed nested owners lose their route,
order and lookup entries. Shared descendants retain their surviving parent; shared
source nodes have one path membership. Cached geometry stays referenced until after
commit, and both the retired source and repository release their leases. No second
scene-sized index copy or catch-only rebuild is used.

An explicit revision change is honored even when the path/key are unchanged. A
normalized unchanged rename remains a no-op. Conflicts which would remove the source
being renamed or leave the conflicting instance reachable elsewhere are rejected
before mutation; ordinary conflict replacement remains supported. Compact occurrence
ownership and distinct source API semantics remain intact.

The `scene-configuration` regression covers nine operations across normal/quiet
notifications and ordinary/batched revisions: draw, representation, rename, conflict,
nested conflict, explicit revision, shared-path rename and compatible external mesh/
wire draw. **36 scenarios sweep 4,868 allocation positions**, with complete graph,
path, ordering/routing, string-key lookup, revision, source-stamp and cache-reference
checks. Ordinary throwing observers inspect indexes immediately; fault sweeps check
allocating string-key lookups after leaving the injection window. Additional tests
cover aliases, normalized no-ops, rejected graph conflicts and reentrant configuration.
The preceding physical library fails both regression groups. The original controller
probe now exits 0 with complete channels, retained indexes and advanced frames for
all three rows. The three touched library translation units also compile independently
with the production warning flags.

The [dated evidence](obol_20260919_readiness_history.md#september-8-scene-configuration-publication)
owns qualification counts and runtime provenance. Other metadata writers, broader
GED/view effects, repository operations outside this retirement path and custom
cleanup callbacks retain their separate acceptance obligations. S1--S6 remain open.

### CXX-CONFIG-003 — metadata setters can interrupt propagation and scene revisions

**Repaired.** Display-name, hierarchy and material-policy setters now use the
prepared source publication and scene commit hook. Source fields, changed owned
metadata, compact effects and compiler/batch invalidation commit before observers.
Preparation failure preserves preceding state; notification failure retains the
complete commit and scene revision. Nested sources retain their ownership, authored
non-database shapes retain their names/materials, and shared children publish once.
Unchanged owned metadata receives no synthetic notification. These changes preserve
realization stamps, staged producer results, geometry and per-occurrence overrides.

Compact occurrences retain their full-path database material color independently
of display policy. Discovery resolves that intrinsic color even under inherited
appearance. Policy changes update RGB while preserving line patterns, widths,
opacity and occurrence semantics. Legacy owned shapes resolve full-path colors
through the existing prefix cache; imported combinations are released if C++
preparation throws. No second occurrence-sized material registry is introduced.
Hierarchy updates the retained renderer parent while preserving each leaf's own
occurrence index and boolean operation. It invalidates compiled compact evidence
because parent identity is outside the existing geometry/style stamps. Display-name
reset restores source-name/path fallback without changing other appearance fields.

The corrected `metadata_probe.cpp` reports complete metadata and advanced scene
frames for all three throwing-observer rows. The original name probe used an authored
non-database diagnostic prototype and incorrectly expected it to be renamed; its
source and failed output remain retained. The corrected probe uses external primary
geometry and fails against the preceding physical library. All three new regression
groups also fail against that library.

Production regressions cover allocation rejection, enabled/quiet notifications,
realized/failed sources, direct/scene/batched calls, aliases, defaults, no-ops,
reentrant callbacks and subsequent input invalidation. A database combination with
subtraction exercises intrinsic colors, legacy/compact material round trips,
renderer records, retained line styling, hierarchy cache retirement and name
fallbacks. The initial focused run caught the missing hierarchy cache invalidation;
its failed output is retained separately from the corrected passing run.
The [dated evidence](obol_20260919_readiness_history.md#september-8-source-metadata-publication)
owns counts, graphical qualification and runtime provenance. S1--S6 remain open.

### CXX-CONFIG-004 — display and placement writers can interrupt their effects

**Closed at the display/placement publication boundary.** `setDisplayState` and
`applyDisplayPatch` now share a prepared display path. A single field mapping filters
source no-ops and applies changed properties to detached source/owned metadata.
Compact membership buffers and a visibility-frontier mask are prepared before commit;
the registry and immutable geometry are not cloned. Numeric style/membership changes,
source invalidation and scene revisions commit before observers. Preparation failure
retains preceding state; observer failure retains the complete commit. Revision
sensors remain attached, and subsequent independent revision changes still invalidate.

Only changed properties propagate. Independent legacy selection, opacity and names,
compact CSG line patterns and intrinsic full-path database colors survive unrelated
updates. Selected/highlighted/normal color variants retain their separate ownership.
Sparse visibility/highlight/transparency overlays keep precedence; registry selection
remains occurrence-owned. Retired whole-target overviews cannot be resurrected.
Visibility journal allocation failure retains the completed state and forces an
authoritative rescan through the existing journal floor. Explicit source/input
revision changes and material-override removal retain their source-invalidation
semantics; ordinary display changes retain realization stamps and producer data.
Authored non-database appearance and independent nested source owners remain intact.

The preceding placement repair prepares source/owned fields, transform nodes and
child/path changes with compact matrices, invalidation and scene revisions. Auxiliary
sources inherit placement metadata without another transform, while ordinary nested
sources keep their own placement. Graph-only repair, aliases, numeric no-ops and
reentrant changes remain covered by the current focused suite. The earlier placement
checkpoint retains its own before/after evidence.

The current `state_probe.cpp` exits 0:

```text
state=display-state threw=1 frame_advanced=1 state_complete=1
state=display-patch threw=1 frame_advanced=1 state_complete=1
state=placement threw=1 frame_advanced=1 state_complete=1
```

Seventy-two display scenarios cover direct/scene/batched callers, ordinary occurrences
and registries, realized/failed sources, ordinary/quiet notification, sparse overlays,
source/input invalidation and persistent allocation failure. Additional production
combination checks exercise CSG line patterns, legacy/compact material precedence,
color round trips, renderer records, unrelated styling and immutable geometry.
Reentrant callbacks, source-revision observer failure, ignored/borrowed arguments,
overlay removal, retired overviews and source no-ops have explicit checks. All three
new regression groups fail against the preceding physical library. The initial
fixture build incorrectly accessed private stamp fields; it now captures public
source revisions. Its failed build and initial test fragment remain retained. The first overlay
assertion treated a whole-target overview as a leaf and consequently hid the
root subtree. The corrected fixture targets explicit leaf paths; both failures
remain separate from qualification.
The [dated evidence](obol_20260919_readiness_history.md#september-9-display-publication) owns counts, graphical qualification and runtime provenance.

This closes these setters, not the remaining direct-writer inventory. Connected
fields/engines, traversal/custom cleanup and broader resource/lifecycle requirements
retain their separate acceptance gates. S1--S6 remain open.

### CXX-CONFIG-005 — remaining scene state setters lose committed revisions

**Closed at these four source/scene setter boundaries.** Bounds and explicit
realization retain their prepared inner source semantics and invoke the scene's
commit hook before observers. Bounds now attempts every changed field notification
and preserves the first exception. Role flags use a fixed scalar publication,
retire compiler/batch evidence before notification and retain unrelated owned state.
View policy prepares all eleven affected fields and owner metadata through the
existing configuration publication, then commits invalidation and scene effects.
Its threshold sensor remains attached; publication suppresses only the already-handled
sensor callback. A subsequent independent threshold change still invalidates.

The old view-policy mutation sequence and role/view-policy full appearance rebuilds
are removed. Geometry, compact semantic/display records, authored auxiliary appearance,
nested source owners, source-production stamps and staged producer data survive.
View-only invalidation retains exact bounds and an external failure's diagnostic.
Role flags retain normalization and do not manufacture realization invalidation.
No-ops allocate nothing at the source boundary; scene batches advance once at end.
Numeric tolerances, nonpositive LoD-scale normalization, invalid status/default
reason handling and borrowed bounds/diagnostic arguments retain their semantics.
Reentrant later changes survive outer observer failure.

The retained `remaining_scene_probe.cpp` now exits 0:

```text
remaining=bounds threw=1 frame_advanced=1 effects_complete=1
remaining=realization threw=1 frame_advanced=1 effects_complete=1
remaining=roles threw=1 frame_advanced=1 effects_complete=1
remaining=view-policy threw=1 frame_advanced=1 effects_complete=1
```

Seventy-two scenarios cover four setters, direct/scene/batched callers, current and
internal/external failed sources, quiet notification and persistent allocation
failure. Dedicated primary-source checks verify compiler effects before throwing
observers; auxiliary sources deliberately use a different rendering path. All three
new regression groups fail against the retained preceding library with verified
loader resolution. The [dated evidence](obol_20260919_readiness_history.md#september-9-scene-state-publication) owns exact counts, fixture corrections, graphical rows and binary provenance.

This does not close enclosing multi-step operations or the wider source/service,
journal and presentation inventory. CXX-CONFIG-008 owns direct source-field
notification acceptance. S1--S6 remain open.

### CXX-CONFIG-006 — database metadata and material refresh can publish partial state

The retained baseline probe interrupts explicit database metadata after region ID, and
compact material refresh after source color validity. Both leave partial records. The
qualified library now produces:

```text
database=metadata threw=1 state_complete=1
database=material-refresh threw=1 state_complete=1 revision=2
```

| Boundary | Prepared state and commit |
|---|---|
| Explicit source database metadata | Normalize invalid color to white and prepare source/owned scalar records and shader storage. Aggregate metadata applies only to matching source paths (including slash/instance normalization); descendants, authored auxiliaries and nested owners retain their records. No full display/owner rebuild runs. |
| Database material refresh | Resolve source and occurrence full paths through one sweep; prepare changed region/color/shader records and strings without copying geometry or the registry. Commit source fields, owned metadata, compact semantics, effective RGB, renderer records and compiler/batch effects before observers. Selection, visibility, line styling and opacity remain intact. |
| Revision/no-change semantics | A repeated material revision avoids lookup/allocation. A compact stamp-only refresh acknowledges the revision and returns 0 without invalidating presentation; a legacy stamp change retains its existing changed return. Explicit metadata no-ops retain borrowed strings without allocation. |
| Scene and bulk callers | Single and bulk callers use the same source publication. Bulk preparation pins its source set and shares path/combination resolution. The scene advances once at the first meaningful source commit; notification failure preserves completed prefixes. Reentrant later material revisions win over the enclosing sweep. Batches defer the visible scene revision until their existing completion boundary. |
| GED metadata application | Remove duplicate noncompact shape writes after the source publication. Direct region writers, sparse metadata setters and enclosing multi-step GED operations remain separately inventoried. |

Twenty-two allocation scenarios and dedicated edge/rendering checks cover these production
paths. Mixed-region colors use librt's independent full-path oracle; shaders and renderer
records retain per-occurrence identity. All three new groups reject the retained preceding
library with verified loader resolution. The
[dated evidence](obol_20260919_readiness_history.md#september-9-database-metadata-and-material-refresh)
owns counts, initial implementation/fixture corrections, graphical regressions and binaries.

Continue sparse metadata/display writers in `database_source_compact_access.cpp` and their
GED region/metadata callers, then the remaining source/service and presentation inventory.
Direct field/engine callbacks, traversal/custom cleanup, complete memory bounds and actual
snapshot-file retirement retain their own acceptance work. This closes CONFIG-006's stated
boundary; **S1--S6 remain incomplete**.

### CXX-SPARSE-001 — sparse entry setters can leave partial records

The preceding library fails allocation probes for display state, selectability, region ID,
region metadata and full metadata. Path traversal and shader/membership allocation can fail
after earlier entries or fields changed. Ordinary material edits also reset authored
transparency and subtraction-line styling through the full display rebuild.

`BObolCompactEntryPublication` now stages only changed scalar/style records, prepares shader
copies and membership capacity, then commits without allocation. Shared intrinsic metadata,
material RGB and numeric-summary helpers serve sparse edits and source-wide material refresh.
The seven covered setters are display state, transparency, selectability, region ID, region
metadata, full metadata and subtraction line style. They retain their existing changed-field
counts and notification behavior. Effective appearance respects retained overlays and the
retired-overview guard; unrelated geometry, placements and source/producer stamps survive.

Compiler invalidation and sparse batch/visibility journals precede callbacks. A journal
allocation denial requests an authoritative rescan after the complete edit. Reentrant later
edits survive outer notification failure. Temporary edited records scale with changed entries;
membership maintenance can scan existing members, and active visibility frontiers retain
their existing transient population mask. This is not a claim that every operation is O(delta).

The three production regression groups cover 28 allocation scenarios, exact/subtree/object
and legacy matching, instance suffixes and prefix boundaries, null/empty/missing queries,
invalid-color normalization, opacity clamping, overlays, retries, no-ops, throwing observers,
reentrant revisions, numeric revisions, exact visibility journals and actual renderer/picking
state. [Dated evidence](obol_20260919_readiness_history.md#september-9-sparse-entry-publication-and-picking-dispatch)
owns qualification and retained baseline failures.

**Follow-up:** SPARSE-002 below closes retained rule/frontier/selection publication;
REGION-001 closes evaluated-region source and scene publication and removes the unused
private region-metadata helper. Other enclosing GED operations, source/service and
presentation ownership retain separate acceptance work; S1--S6 remain open.

### CXX-PICK-001 — direct retained-assembly ray traversal skips picking

The new sparse renderer check found that even its unchanged selectable wire was not picked
when traversing the retained assembly directly. Native reproduction and a debugger trace
confirm that `SoRayPickAction` never dispatched `SoCADAssembly::rayPick`: a `SoNode` subclass
inherits generic pick dispatch unless it registers the ray action explicitly. BRL-CAD's
source wrapper explicitly invokes the assembly's picker and therefore masked this defect.

`SoCADAssembly::initClass` now registers the existing `SoNode::rayPickS` dispatcher. The
native `CadAssemblyPicking.WireRayTraversalHonorsSelectionAndVisibility` regression rejects
the preceding library and exercises translated wire picking, unpickable/hidden suppression
and restoration through normal action traversal. The sparse tests also exercise the derived
BRL-CAD assembly. Numeric picking algorithms and source-wrapper ownership are unchanged.
The same [dated evidence](obol_20260919_readiness_history.md#september-9-sparse-entry-publication-and-picking-dispatch)
records the dependency build, native tests and shared-client qualification.

### CXX-SPARSE-002 — retained intent can disagree with current presentation

The preceding library publishes retained rules or paths before allocating their affected
entry updates. Eight of nine public-operation probes expose partial publication; the ninth
also shows that removing a selected child path clears entries still covered by its retained
parent. Selection delta matching differed from replacement and later streamed arrivals.

The existing `BObolCompactEntryPublication` now prepares retained intent and affected
scalar/style records together. Rule/path vectors, membership capacity and changed-entry
storage are ready before intent, cached overlays, effective appearance, aggregate selection,
renderer records and compiler/journal effects commit. Shared pure helpers resolve both
prepared records and private construction; the three live reapply functions are removed.
Immutable geometry, dense indexes and producer ownership remain in place. Exact rule edits
visit indexed matches; frontier replacement retains its existing transient population mask.
Selection deltas store affected targets and resolve their state against final retained paths,
including remaining parent coverage, object queries and instance suffix normalization.
Unrelated immediate selection survives a delta.

A meaningful retained-intent change now notifies the source even when it has no current
matching occurrence. Identical/no-op and invalid input remain silent. Aggregate selection
field notification and source notification are both attempted after complete publication,
preserving their independent enable flags. Later reentrant edits survive an outer observer
failure. These statements cover the setter boundary; connected-field/engine paths
outside the twelve CXX-CONFIG-008 source fields and enclosing multi-source GED
operations retain their separate acceptance work.

Three regression groups cover 36 allocation scenarios, current and later-arriving entries,
empty sources, quiet notifications, independent expected appearance/revisions, rule
precedence, exact visibility journals, authored restoration, retries, throwing observers,
reentrant edits and overlapping selection paths. The [dated evidence](obol_20260919_readiness_history.md#september-9-retained-presentation-and-selection-publication)
owns counts and qualification. REGION-001 below covers evaluated-region publication.
Other enclosing GED operations and the source/service and presentation inventory remain
open, as do S1--S6.

### CXX-REGION-001 — evaluated-region edits expose partial shapes and miss scene effects

The preceding GED setter writes wire fields individually, updates only the first mesh,
and omits the scene frame effect. An immediate observer sees partial state; even ordinary
completion leaves the second mesh unchanged. The compact route also lacks that enclosing
scene effect. The retained production probe reproduces all three observations.

The source now prepares only changed region integers on matching owned realized shapes.
Shared path matching preserves compact object/subtree and instance-suffix semantics.
Aliases are visited once; auxiliary shapes and nested source owners remain independent.
Unrelated material/region fields, immutable geometry and producer stamps survive.
Compiler/batch retirement and the scene's existing frame effect commit before source or
field observers. Notification flags are restored and all prepared notification attempts
are made before propagating the first observer exception. Compact sources use the existing
sparse publisher with the same scene commit contribution.

The scene operation pins and validates its original targets before editing. Each source
commits independently; allocation or observer failure preserves a completed prefix.
Before applying later targets it checks owner identity, CAD revision and node identity,
so a callback's replacement or later edit, including a quiet edit, supersedes the old
request. The query is copied before callbacks can invalidate borrowed path storage.
Existing explicit mutation batches retain their coalesced frame revision. GED resolves
targets and delegates; the unused private region-metadata setter and evaluated-region
getter are deleted after checking their callers.

Three production groups cover eight source/publication allocation scenarios and two
multi-source prefix sweeps, independent expected fields and retained geometry, throwing
observers, quiet edits, owner replacement, retry/no-op, path boundaries, source picking
readback and compiler retirement/rebuild. Compact checks inspect actual renderer records.
The real GED regression requires every wire/mesh and the frame effect to agree inside
an observer and after an observer exception. Its negative binary links the retained
preceding static GED/BObol archives because these private GED entries are not exported
by the shared library. [Dated evidence](obol_20260919_readiness_history.md#september-9-evaluated-region-source-and-scene-publication)
owns exact counts, corrected test fixtures, graphical rows and binary identities.

GED-001 below covers draw-metadata scene effects and material-refresh target acceptance.
Enclosing direct/deferred GED draw transactions and the remaining source/service and
presentation inventory stay open. Connected fields/engines, custom traversal/cleanup and
measured memory bounds retain separate acceptance work. S1--S6 remain open.

### CXX-GED-001 — metadata application omits scene effects and refresh overwrites later targets

The preceding GED metadata adapter calls prepared source setters without the scene's
frame contribution. A source-field observer sees new metadata with the old frame revision;
quiet compact metadata also leaves that revision unchanged after ordinary completion.
Its material-refresh loop resolves each next key after the preceding source's callbacks,
so it overwrites a later material revision or a replacement owner. The retained production
probe reproduces all four failures.

The scene controller now normalizes the existing `BObolDrawMetadataRecord` once and
contributes its frame effect to the existing aggregate and sparse publishers. A supplied
path selects compact occurrence metadata; an omitted path retains aggregate semantics.
Noncompact sources retain aggregate scope and preserve independent descendant metadata.
GED-001 routed three metadata sites through this boundary; GED-002 below removes the
unreachable direct-leaf site. Valid unchanged metadata
application reports success to GED without advancing the frame, including compact no-ops.
The duplicate GED converter is removed.

Targeted and full-scene material refresh use the existing database sweep. Targeted refresh
preserves each source's database fallback unless an override is supplied. The sweep reuses
one cache while that database remains current and replaces all cached evidence when it
changes; a source without a database cannot suppress a later valid target. Full-scene
refresh and the explicit-database free helper retain their established fallback/error
semantics. The sweep pins original sources and captures material, CAD and node revisions
before callbacks. A scene
acceptance check rejects removed/replaced owners; revision checks preserve later changes,
including quiet metadata edits without a material revision change. Each changed source
contributes its frame effect before observers. The former first-source-only suppression
is removed; explicit mutation batches retain their existing coalescing. Allocation and
observer exceptions preserve every completed source and its effects. Compact stamp-only
refresh keeps its established acknowledgement/no-frame behavior.

GED shares source-target resolution between evaluated-region edits and material refresh.
The direct-leaf adapter pins its owner and copies its key across metadata callbacks before
the following line-style operation; a replaced owner is skipped. This does not establish
combined acceptance of that multi-step draw operation. In particular, later edits to the
same owner and the enclosing transaction's cached source/index references still need audit.

Production regressions add scene/batch aggregate metadata fault sweeps, compact metadata
sweeps with independent expected RGB and renderer readback, and targeted/full-scene
material-prefix sweeps. Observers check each committed source/frame pair. Later material
revisions, quiet metadata changes, replacement/removal, current-source reentrant changes,
throwing observers, per-source database routing, explicit override, missing databases,
retry/no-op and explicit batches have direct coverage. The real GED test checks both
metadata representations, throwing source observers, unchanged application
and later refresh targets. Its negative binary links the retained preceding static
GED/BObol archives. [Dated evidence](obol_20260919_readiness_history.md#september-9-ged-metadata-and-material-refresh-effects)
owns counts, fixture corrections, graphical rows and exact runtime identities.

The growing publication executable exceeded its 180-second CTest timeout. CTest now runs
explicit source and state groups with that same timeout per group. All 71 direct checks
and 21 installation/merge scenarios remain; execution markers verify that each runs once.
The no-argument full run and individual selectors remain available. No application latency
criterion is changed by this test-harness partition.

The follow-up caller audit below replaces the direct-leaf task with deletion of unreachable
code and a reproduced failure in the live scene publisher. All S1--S6 remain incomplete.

### CXX-GED-002 — unreachable leaf redraw machinery

The preceding `ged_obol_apply_draw_paths_to_scene` entry gate and root predicate accept
the same six draw modes. Each mode returns before the leaf walk, so its instance-key
helpers and `ged_obol_direct_apply_leaf_state` cannot execute. The sole caller initializes
default appearance and changes only mode/mixed-mode preservation; its own sole caller
always requests display preservation. Deferred selection, strict fallback, result output
and the alternative display-preservation branch are therefore unused in this helper too.

The redraw helper now accepts the existing path vector and mode, and retains only the
live root-source operation. The separate live deferred adapter remains. The cleanup
removes 438 net lines, including the unreachable metadata/style adapter that GED-001 had
pinned defensively. Its suspected callback-composition issue was not a reachable defect.
The source/caller audit and before/after debugger traces establish reachability; the actual
six-mode redraw probe preserves geometry and line width on both binaries.
[Dated evidence](obol_20260919_readiness_history.md#september-9-unreachable-ged-leaf-redraw-removed)
owns the regression results and retained artifacts.

### CXX-CONFIG-007 — scene source publication exposes intermediate setters

The preceding `BObolSceneController::publishDatabaseSourceInstance` composed individually
prepared setters by executing them in sequence, then attached/moved/indexed the source
and advanced the frame. An existing source's revision observer could see its new revision
with the old line style and old frame revision. A newer line-style edit from that observer
was then overwritten by the outer display setter. The retained `publication_probe` failed
both before and after GED-002; that cleanup changed libged only.

**The reproduced boundary is now repaired.** The source normalizes the existing
`BObolDatabaseSourcePublishState` through setters on a detached scalar candidate. Its
prepared state owner composes owned metadata, compact appearance and optional placement;
the scene prepares group creation, child-list changes and indexes before committing.
Existing rename index preparation and compact style helpers are shared. Database binding,
compiler/journal effects and scene revisions commit before observers; no older setter or
index write follows callbacks. Prepared references retain new/existing sources and parents.
The publish-state declaration now belongs to the source header, preserving the scene's
one-way dependency on source data. Unchanged explicit state returns zero without a frame
change; a pure parent move changes only the scene.

Three production regression groups cover old-or-complete allocation failure, complete
callback observation, quiet fields, batches, retargeting/shared paths, insertion, missing
groups, retry/no-op, pure moves, cycle rejection and callback style/removal/replacement/move.
Actual compact renderer records retain geometry while changing placement, line style,
opacity and full-path materials. Transaction and callback tests reject the retained
preceding library. [Dated evidence](obol_20260919_readiness_history.md#september-9-composed-scene-source-publication)
owns the focused, independent-compilation, GED and graphical results.
The enclosing GED redraw/adoption/deferred callers still require their own current-target
and reference-lifetime audit. S1--S6 remain open.

### CXX-GED-003 — redraw resets and restores current appearance

The root-source redraw saved a scene-wide list of matching appearance records, published
default appearance through `ged_obol_replace_path`, then restored the saved records by key.
All six accepted modes reproduced a line-width callback seeing 1 instead of the preceding
3; its newer edit to 7 was overwritten by the saved 3. Replacing a source under the same
key could also receive the preceding owner's appearance.

The redraw now passes no appearance override. The existing publisher obtains appearance
from the exact current source while constructing its complete candidate. The saved-record
type, collecting/restoring writers and default settings are deleted. An inert root metadata
record containing no material color is also removed; it selected the same material-resolution
branch as a null record. The simplification removes 118 net lines and preserves all six modes.
Twenty-four production GED cases qualify uninterrupted appearance and callbacks that edit,
replace or remove the source. The same test object linked with the retained preceding GED
archive rejects the transient appearance. [Dated evidence](obol_20260919_readiness_history.md#september-9-redraw-appearance-preserved-directly)
owns final builds, focused/graphical checks and before/after probes.

This closes the source-appearance restore writer. It does not close the following group
composition failure, base-source promotion, deferred-proxy clearing, retained-current
marking, broad/narrow occurrence adoption or other deferred target/reference acceptance.

### CXX-GED-004 — group synchronization overwrites a source callback's edit

After publishing the source, the preceding `ged_obol_replace_path` called
`ged_obol_sync_group_state` with its earlier source-state copy. The source callback can
set its group's line width to 7 through `setGroupDisplayState`, after which the enclosing
redraw resets the group to 3. `group_probe` reproduces this in all six modes on the
GED-003 implementation; source appearance and realized geometry remain correct.

**The reproduced source/group boundary is repaired.** A scene-level group-state overload
composes intended group metadata with the existing source, membership, index and revision
publication. Group fields commit and regain their notification policy before source
observers; no saved group write or regrouping follows those observers. The source-only API
remains available. A group-only change advances the frame without changing source state;
invalid metadata targets fail without publication. GED's obsolete group-sync helper is
removed, reducing `draw_obol.cpp` by another 50 net lines.

The standalone group intent/display setters now prepare a retained candidate and their
frame effect before observers. Source and group owners share the existing scalar copy,
commit and notification helpers through a private header; geometry and source-specific
ownership remain in the source implementation. No second notification framework or global
callback prohibition is added.

Four production regression groups cover seven source/group allocation scenarios and four
standalone group scenarios, callbacks in both directions, quiet fields, no-ops, group-only
effects, batches, removal/replacement, exceptions, aliased intent and invalid targets.
The six-mode GED callback regression extends the preceding appearance tests to 30 cases.
[Dated evidence](obol_20260919_readiness_history.md#september-9-source-and-group-publication-composed)
owns baseline rejection, allocation counts, focused/graphical checks and runtime identities.

This qualifies metadata in the composed publisher and the two standalone metadata setters.
The separate group hierarchy entry points are qualified below. GED-005--013 close base
promotion, deferred-proxy clearing, retained-current marking, complete source-record
application, resident-occurrence adoption, live deferred-result acceptance and compact
proxy bootstrap, including the enclosing retained-current source sequence, cross-scene
source copy and core presentation writers. S1--S6 remain open.

### CXX-GED-005 — base-source promotion publishes identity before representation

The live base-source promotion renamed the source instance, notified observers, and then
set its representation before the enclosing full publication. A callback could therefore
find the new key with `REPRESENTATION_DEFAULT`; a later callback edit could also be
overwritten by the remaining writes.

**The transition is repaired.** `publishDatabaseSourceInstance` accepts an optional current
instance key and prepares current-to-next identity, representation, source state, group
state, indexes and frame effects as one publication. It rejects a missing required current
key and an independently owned destination. The GED caller no longer performs a preliminary
rename or representation setter. Its live regression enters the real publication context,
observes the new key and shaded representation together, then changes line width; the same
source remains current with that edit after the outer call.

The adjacent local mesh/vlist factory audit found no source callers and no exported dynamic
symbols. Those factories and their common unsafe helper returned borrowed shapes after
append/move callbacks, but no live contract depended on them. Removing about 320 lines
eliminates that acceptance surface instead of adding machinery to dead code.

The saved preceding GED archive rejects the base-promotion regression with the exact
default-representation observation. Current source and source/group transition allocation
sweeps cover 463 and 493 positions. [Dated evidence](obol_20260919_readiness_history.md#september-10-live-ged-source-transitions)
owns the retained before/after logs and final qualification.

### CXX-GED-006 — deferred proxy state and ownership publish separately

Two live deferred paths exposed partial source records. Base-source replacement first
cleared external primary geometry and marked the source stale before publishing native
ownership. Retained-current marking first added external ownership and then set realization
current. In both cases an observer could see ownership disagree with realization status;
the later call could reach or overwrite a callback-time replacement.

**Both transitions now use complete source publication.** External-to-native role retirement
adds `STALE_SOURCE` inside `publishState`, while retaining the proxy as visible continuity
geometry until native data replaces it. The preliminary clear/stale sequence is removed.
The realization setter has a composed role overload which uses the same prepared notification
envelope as status, revisions, diagnostic and staleness. Retained-current marking calls it
once. Source and owned-shape fields commit with the scene frame before callbacks; notification
restoration occurs once even though the source role field participates.

The source-level regression covers a throwing observer and a later reentrant role edit.
The live GED regression covers callback change, same-key replacement and removal, and checks
the source, retained shape and committed frame. Saved pre-repair libraries reject external
retirement and retained-current marking; current checks pass. [Dated evidence](obol_20260919_readiness_history.md#september-10-live-ged-source-transitions)
owns exact artifacts and timings.

### CXX-GED-007 — source-record application publishes seven field families separately

The private GED source-record bridge represented one logical source state but applied it
through seven scene-controller calls: draw mode, representation, display/revisions,
material policy, realization, roles and view policy. An observer could see any completed
prefix with the old frame state, and later calls could overwrite a callback edit or reach a
same-key replacement. The ordering had a second deterministic defect: it marked a record
current before changing view policy, whose invalidation could leave the requested current
record stale.

**The bridge now performs one complete source publication.** The existing
`BObolDatabaseSourcePublishState` carries an optional explicit realization record. Source
publication applies configuration and view invalidation to its private candidate first,
then applies the authoritative realization status, revisions, diagnostic and stale reason.
Configuration, representation, display, material policy, role ownership, view policy,
owned-shape state and the scene frame commit before observers. The GED bridge preserves
the source's existing appearance while translating its record and calls the current-to-
current scene publisher once.

The source and source/group allocation sweeps now include explicit current realization in
every publication kind. The live GED regression changes all seven former field families
and checks their complete values inside the first callback, then covers a later display
edit, same-key replacement and removal. The identical current test object linked with the
saved preceding GED/BObol archives exits on the first partial callback; the repaired binary
passes all three cases. [Dated evidence](obol_20260919_readiness_history.md#september-10-live-ged-source-transitions)
owns the exact binaries, logs, counts and focused qualification.

### CXX-GED-008 — resident-occurrence adoption notifies before donor retirement

The live broad-source redraw path adopted each narrow source's compact occurrences into
the target and then removed the donor instance in a later controller call. Adoption
mutated the live target incrementally and could notify its observers while the donor edge,
source index and scene revisions still described the old ownership. The GED caller also
cached a raw target pointer across callbacks. Allocation failure during a multi-donor
merge could therefore expose a prefix of the requested handoff.

**Resident data transfer and donor retirement now form one controller transaction.** The
target builds its merged compact index and dependent presentation state on a detached
candidate. The controller prepares every selected donor edge, the complete source index
and affected revisions before mutation. It commits the target and retires the donors
before source or hierarchy notification, including the case where all donor occurrences
already exist in the target. GED groups donor instance keys by target key and lets the
controller resolve and retain the participating objects. Callback edits, same-key
replacement and removal remain current. The old public payload-only adoption method was
removed once its callers used the composed owner.

The transaction regression covers **239 allocation positions**: 236 preserve the old
scene and two complete the publication before notification failure. Edge coverage checks
complete callback state, edit, replacement, removal, throwing observers and invalid
targets. The live GED regression then exercises its natural narrow-to-broad transition,
including compact occurrence and selection survival. Saved predecessor archives export
and call the split operation but lack the composed controller entry point; the current
contract fails to link against that predecessor and the current executable passes. The
seven-check focused qualification is retained under
`.build/obol-qualification/20260910-ged-live-callers`.

### CXX-GED-009 — deferred result publication uses stamped scene acceptance

The live progressive pump merged occurrence batches directly into a borrowed source and
adopted the completed detached source through another source-local operation. It then
published coverage, profile, expected population, scene effects and proxy cleanup in
separate steps. A callback could therefore observe the new source payload with an old
frame or replace/remove the target before a later write through the cached pointer.
Allocation failure from an observer was also indistinguishable from failure before the
requested batch or terminal handoff committed.

**Incremental and terminal delivery now enter through the scene owner.** Each entry point
accepts the launch stamp and resolves the current routing identity before mutation. A
stream batch commits its internally consistent occurrence prefix and every affected frame
revision before source field notification. The full-batch witness remains false for an
internal partial merge and becomes true before notification, allowing the pump to retain
a completed batch when an observer throws. Coverage publication has the same commit
witness. Terminal adoption composes final stream bounds/profile/population, realization
fields, source-local children/indexes and frame effects before callbacks.

GED stores keys and launch stamps across the callback boundary, then resolves the exact
current source again before the next batch or cleanup. Terminal descendant cleanup uses a
pre-callback snapshot of instance key, routing id, parent and ancestors and revalidates
each identity deepest-first. This prevents cleanup from removing a callback-created
replacement or reparented source. Independent primary/view scenes still execute their own
validated cleanup transactions; they are not represented as one cross-scene atomic edit.

The source-level terminal check covers 21 allocation positions. Scene-level terminal
publication covers **135 positions** (84 preserved and 50 complete commits), and streamed
publication covers **88 positions** (39 preserved, 39 retained prefixes and nine complete
batches). Callback edges cover frame/bounds completeness, later edit, same-key replacement,
removal, exceptions and invalid targets. The live GED regression exercises coverage and
terminal callbacks plus stale-result rejection and delivery-denial paths. Seven focused
checks pass in 16.13 seconds under
`.build/obol-qualification/20260910-ged-deferred-acceptance`. The retained symbol comparison
predates GED-008 and is therefore API evidence rather than an immediate-predecessor
behavioral baseline.

### CXX-GED-010 — compact proxy bootstrap publishes one source snapshot

Four live GED bootstrap paths installed a compact occurrence registry directly on an
attached source, then published exact bounds, source profile and realized/current state in
later calls. The view-envelope path also cleared external primary geometry before the
registry install. Each individual write could notify, and the caller retained a borrowed
source pointer across those notification boundaries. An observer could therefore see a
new registry with old realization state or no primary presentation, and a later write
could reach a callback replacement.

**Compact bootstrap is now one stamped controller publication.** The source builds the
complete compact index on a detached template, applies optional producer-certified bounds
and aggregate profile there, and then revalidates the exact routing/configuration stamp.
The existing prepared source publisher commits the index, terminal realization record,
external role, primary-child replacement, scene lookup effects and frame effect before it
notifies. Auxiliary children remain attached. A completion witness is set at commit, so
the GED adapter preserves a complete snapshot when an observer subsequently throws
`std::bad_alloc`; other preparation failures leave the preceding source unchanged.

The old cached-manifest, decoded-manifest, synchronous structural and view-envelope
sequences now call `publishDatabaseSourceInstanceCompactSnapshot` and retain no source
pointer afterward. The view-envelope transaction retires the temporary external line in
the same prepared child order. A deferred root now reads a cached AABB into the compact
snapshot directly instead of creating an external line which it would immediately replace.
For leaf sources which retain an external AABB, the separate exact-bounds setter was
redundant—the external source publication already derives and commits that exact bound
with the line—and has been removed.

The allocation regression covers **263 positions**: 201 preserve the old source and 61
commit the complete snapshot before a notification allocation failure. Callback cases
cover later edit, same-key replacement, removal, exceptions and invalid/stale targets. The
live GED deferred draw observes the compact replacement only after its complete source
state is visible, changes line width, throws `std::bad_alloc`, and proves both the committed
snapshot and later edit survive. Eight focused checks pass in 16.50 seconds under
`.build/obol-qualification/20260910-ged-compact-snapshot`. The retained GED-009 symbol log
lacks the new controller entry point; current libBObol exports it and current libged calls
it. This is an immediate API comparison plus current behavioral evidence, not a full
client or graphical qualification.

### CXX-GED-011 — retained-current source policy publishes before current state

The retained-current state helper repaired under GED-006 was complete by itself, but its
live caller still invoked it after `publishDatabaseSourceInstance`. A view-policy change in
that first publication marked the retained external realization stale and notified before
the helper restored current state. The caller then resolved the instance key again, so a
callback-time same-key replacement with retained geometry could receive that later current
write.

**The enclosing operation now publishes one source record.** When live GED retains current
external geometry, it includes the existing role ownership and the final current
realization record in `BObolDatabaseSourcePublishState`. Configuration and view-policy
invalidation occur on the private candidate before that authoritative state is applied.
The subsequent source lookup and current-state publication are gone, together with the
private helper declaration and implementation.

The live regression enables a view-dependent policy, constructs a current external wire
with a deliberately old policy, and invokes the real GED source ensure path. Its first
callback sees the new policy, current realization, external ownership, owned-shape state
and instance lookup together. A display edit remains current, removal remains absent, and
a same-key replacement deliberately left stale with external geometry is not changed on
unwind. The current test passes all three cases. The same test object linked with the saved
GED-010 archive exits on the partial callback; that archive contains the follow-up helper
symbol while the current library does not.

Four focused checks pass in 23.60 seconds under
`.build/obol-qualification/20260910-ged-retained-sequence`: the seven-case source
publication allocation sweep, source publication edges, full GED draw sync and its CTest
row. This is a focused source/GED boundary, not full client or graphical qualification.

### CXX-GED-012 — cross-scene source copy applies material policy later

The secondary-scene copy path translated a primary source summary into
`BObolDatabaseSourcePublishState`, published the source, and then applied its material
policy through a separate instance lookup. A root observer saw an inserted source with the
default inherited policy before the requested database policy. The later lookup could also
change a same-key replacement installed by that callback.

**Material policy is now part of the source copy transaction.** The summary translator sets
`materialPolicyValid` and `materialPolicy` before calling the existing complete source
publisher. The second setter, result merge and instance lookup are removed.

The live regression creates a database-policy source in the primary scene and invokes the
real reducer draw path against a separate controller. Its root callback sees identity,
representation, revisions, roles and material policy together. A later display edit
survives, removal remains absent and a same-key replacement retains its inherited material
policy. The same test object linked with the saved GED-011 library exits on the partial
insertion callback. Four focused checks pass in 23.90 seconds under
`.build/obol-qualification/20260910-ged-scene-copy`: the seven-case source publication
allocation sweep, source edges, full GED draw sync and its CTest row. This remains a focused
source/GED boundary.

### CXX-GED-013 — presentation commands publish compact and aggregate state separately

GED visibility, highlight and transparency handling applied retained compact rules,
source display state and matching group state through separate calls. A source observer
could therefore see the new compact state with the old aggregate state, or the reverse.
Global highlight clearing also removed retained overrides before publishing the source and
groups. A callback edit, removal or replacement of a later target could be overwritten by
the rest of the original loop.

**The presentation family now has one scene transaction.**
`BObolScenePresentationTransaction` captures every source and group target before the first
publication. Each source prepares its aggregate display fields, retained compact rules,
effective occurrence visibility/highlight/transparency, renderer records, sparse visibility
effects and frame effect, then commits that complete record before observers. Groups use
the existing prepared scalar publisher. The multi-target contract is a consistent prefix:
allocation or observer failure leaves only complete targets committed, while a callback
edit, removal or replacement supersedes pending work for that target. Invalid or duplicate
targets reject the entire request before publication.

The real GED reducer regression covers highlight and global highlight clear on a realized
compact wire source. Its first callback sees the source aggregate, every occurrence and the
frame revision together; a callback line-width edit remains current. The scene regression
also covers simultaneous retained leaf visibility, highlight and transparency, renderer
opacity/visibility, a pending source, a pending group, prevalidation, and **106 injected
allocation positions**: 71 preserve the preceding scene and 36 leave complete prefixes.
The isolated current contract linked with the saved GED-011 archives exits on the split
highlight callback; presentation behavior was unchanged by GED-012. Five focused checks
pass in 16.25 seconds under
`.build/obol-qualification/20260910-ged-presentation`. This is a focused C++ boundary; the
stage estimates and S4--S6 qualification remain unchanged.

### CXX-GED-014 — public source display publishes configuration and appearance in pieces

The database-source display bridge set draw mode, representation and display/material
fields with three controller calls. Its first callback could see only the new draw mode,
and a later call could overwrite a callback edit. Public shape references compounded the
problem by updating the instance key and then its database-path alias, so the same source
could be published twice.

**The bridge now publishes one complete retained source summary.** All matching source
instances are resolved, retained and snapshotted before the first publication. The shared
summary translator fills one `BObolDatabaseSourcePublishState`; requested draw,
representation, display and material changes are applied to that state and the existing
complete source publisher commits it once. Each subsequent target is accepted only while
the retained node and key are still current, so callback changes, removal and replacement
supersede pending work. The public shape path resolves its canonical database path before
publication and no longer performs a second alias write. The public group-reference path
continues to use its single prepared group display publication.

The real GED regression covers public shape and group visibility callbacks, a simultaneous
draw/representation/display/material source update, and two same-path source instances.
Callbacks see one complete record and their edits survive; pending-target edit, same-key
replacement and removal also survive. The same test object linked with the saved GED-013
libraries exits on the overwritten public shape callback. Five focused checks pass in
21.73 seconds under `.build/obol-qualification/20260910-ged-source-display`: the isolated
GED contract, complete scene source publication sweep, adjacent display edges, full GED
draw sync and its CTest row. This closes the public shape/group-reference display row only;
the stage estimates and S4--S6 qualification remain unchanged.

### CXX-GED-015 — deferred selection callback can supersede its source target

Deferred realization validated a captured source and then initialized its compact
selection from the semantic selection service. That initialization notifies observers.
A callback could retarget, replace or remove the source, after which startup continued to
dereference the borrowed pointer and could launch work for an owner which no longer
matched the captured path and representation.

**Deferred startup now retains and revalidates its exact source.** A shared predicate
checks routing identity, path and representation before selection initialization. An
`SoNodeRef` retains the source across the callback. Startup then requires the captured key
to resolve to that same node and repeats the predicate before it captures the realization
stamp, database snapshot and detached worker template. Callback changes therefore
supersede the saved target without stale worker creation.

The real GED reducer regression stages semantic selection and observes the compact source
created by a deferred wire draw. Its selection callback exercises in-place retarget,
same-key replacement and removal. Each edit survives, and one bounded provider pump
reports no source-preparation work. The same test object linked with the saved GED-014
libraries exits with preparation still pending after the retarget. Five focused checks
pass in 21.62 seconds under
`.build/obol-qualification/20260910-ged-deferred-selection`: the isolated contract,
complete scene source publication sweep, adjacent display edges, full GED draw sync and
its CTest row. This closes the deferred-selection initialization row only; the stage
estimates and S4--S6 qualification remain unchanged.

### CXX-GED-016 — streamed subtract style is missing from realization publication

Direct database realization carried a transient dashed flag while building its compact
candidate, but streamed `BObolCompactOccurrence` records carried only their durable
boolean operation. Reconstructing those records therefore published subtractive wire
occurrences with a solid line. A later GED redraw loop called
`setCompactSubtractLineStyle`, but it ran only after eager realization and was already a
no-op there; it neither fixed streaming nor belonged to the realization transaction.

**Subtract style is now intrinsic to compact occurrence construction.** The compact
builder selects the named subtract line pattern when either the direct walk flag or the
durable boolean operation identifies subtraction. Direct, streamed and reconstructed
occurrences therefore publish boolean identity, style and renderer records together. The
redundant post-realization GED writer is removed; the existing compact setter remains the
explicit API for later style overrides.

The direct compact regression observes a newly installed subtract occurrence and sees
its boolean operation and dashed style in the same callback. Existing sparse style
override/allocation coverage continues to pass. The real deferred GED producer now
publishes `/progressive_root.c/box.s` with subtract operation and dashed style. The same
test object linked with the saved GED-015 libraries exits with that occurrence at
`style=0`. Six focused checks pass in 18.54 seconds under
`.build/obol-qualification/20260910-ged-subtract-style`, including the preceding
deferred-selection repair, full GED draw sync and its CTest row. This closes the compact
subtract-style row only; the stage estimates and S4--S6 qualification remain unchanged.

### CXX-GED-017 — database rename publishes source, compact, group and repository state in pieces

The GED rename reducer previously sequenced repository-key rename, compact-path
retargeting, per-source identity changes, saved-state restoration, group updates and a
later realization pass. Each step could notify independently. A callback could therefore
observe old indexes or group paths beside a new source path, and the restoration step
could overwrite a later callback edit. The path helper also treated the renamed object as
the complete source path, so renaming a nested component could collapse a retained path
instead of changing that component in place.

**Database rename is now one prepared scene operation.**
`BObolSceneController::renameDatabaseObject` prepares repository residency for every
affected controller, complete source identity and owned-shape metadata, descendant parent
keys, compact occurrence semantics and renderer hierarchy, group paths, draw-intent paths,
Coin names, source/group indexes and peer revisions. It commits all of them before any
observer runs. Path-derived GED instance keys are declared explicitly; opaque keys remain
stable. The GED bridge only classifies a key as path-derived when it actually embeds the
source path, calls the composed operation once, and realizes the committed scene afterward.

The allocation sweep covers **658 positions**: 649 preserve the preceding scene and eight
fail only after a complete commit. It includes two renamed group occurrences, a path-derived
source, an opaque-key descendant whose path and parent change, compact semantic paths with
stable handles, rebuilt renderer-parent records, same-root and split-root controller peers,
and three independent repositories. Callback edit, removal and exception cases retain the
complete commit; invalid input preserves the scene. The real GED regression proves complete
first-callback state, retains a later display edit and renames a nested database component
without collapsing its source path. The same test object linked with the saved GED-016
libraries exits 1 on the partial callback contract. Ten focused checks pass in **36.84
seconds** under `.build/obol-qualification/20260910-ged-rename`, including adjacent scene,
group, repository, deferred-selection and subtract-style checks, full GED draw sync and its
CTest row. This closes the rename/compact-path row only; the stage estimates and S4--S6
qualification remain unchanged.

### CXX-GED-018 — redraw and erase use live order, global realization and duplicate visibility writers

Redraw formerly traversed the controller's live source array while callbacks could remove
an earlier source and shift the next source into the visited index. Its scoped forms then
called global realization without carrying their view scope, so unrelated work could run
while an accepted source remained pending. Exact erase changed compact occurrence authored
visibility before removing source and group edges in separate notifications. After source
and group removal were composed, its follow-up scan could still select a callback-created or
newly edited retained source. Overview-to-leaf replacement also restored the authored
baseline and discarded that immediate visibility edit. Pure prefix/scope removal started a
global realization pass after its callbacks.

Redraw now captures stable instance keys and exact realization stamps before publication,
including wildcard-mode requests. `realizeDatabaseSourceInstance` accepts only the current
key/stamp pair and confines compact descendant cleanup to that source. Mode and view scope
also filter the accepted visibility owners. The visibility publication uses durable retained compact
overrides, so the requested state survives overview replacement. Exact and prefix erase use
one `BObolSceneRemovalTransaction`. A `BObolSourcePresentationStamp` captures source lifetime
and compact-presentation revision before removal; the later presentation transaction skips
an edited, removed or same-key replacement target. Only compact owners which currently
contain or semantically cover an erased path receive a retained visibility rule. Pure
removal operations no longer call global realization.

The scene removal sweep covers **206 allocation positions** (197 preserve the old scene and
eight complete the commit); the presentation sweep covers **106** (71 preserved and 36
complete prefixes). Realization coverage includes 88 callback histories and eight
source-root cases, 296 descendant allocation positions, six cleanup histories and 248
cleanup allocation positions. The live GED regression covers compact-root retirement,
callback traversal shifts, scoped realization, same-key replacement during visibility,
representation-scoped redraw/erase, a newer retained visibility edit and callback-created
pending work after prefix erase. Nine focused checks pass in **17.22 seconds** under
`.build/obol-qualification/20260910-ged-redraw-erase`, including full GED draw sync and its
CTest row. No immediate predecessor was sealed during the four live reproductions; the saved
GED-016 library lacks the new removal, targeted-realization and presentation-stamp entry
points, so it provides API evidence only. This closes the redraw/erase row at the focused C++
boundary; resource, visual, native-platform and release qualification remain open.

### CXX-GED-019 — mesh-LoD handle and bounds attach through split path writers

The BoT LoD path formerly transferred a raw reader to the first source matching a semantic
path, then published its bounds through a second setter to every source with that path. The
source owned the reader after the first call even when cut loading or bounds inspection later
failed. The bounds check also admitted non-finite values. Source/input/database invalidation
left the reader attached, and the public source setters did not state who owned a rejected
pointer or whether their private state required a scene notification.

Mesh-LoD preparation now keeps the candidate in a local RAII owner until its hierarchy,
initial cut and finite ordered bounds are complete. The GED shape token selects one exact
instance key. `BObolSceneController::adoptDatabaseSourceInstanceMeshLod` accepts that key
only with the realization stamp captured before preparation, then installs the reader and
bounds as one source-resource record. A rejection leaves ownership with the caller; a commit
publishes the new record before destroying the preceding reader. Source, input and database
identity invalidation retire both reader and bounds through the same resource-authority
helper. View-only invalidation retains them. This private service record does not advance a
scene frame; the separately composed mesh-geometry publication owns that observable effect.
The split public source setters, both path-based GED setters and their dead runtime-summary
fields are removed.

The live GED regression covers view-only reuse, non-finite bounds rejection, atomic transfer,
stale-input rejection, retirement/rebuild and exact isolation between wire and shaded
representations. The mesh cache CTest and the source-invalidation/configuration allocation
sweeps also pass. Six focused checks pass in **22.26 seconds** under
`.build/obol-qualification/20260910-ged-mesh-lod`, including full GED draw sync and its
CTest row. The saved GED-016 libraries contain the former split source and GED entry points
and lack the exact controller transfer, so they supply API evidence only. The AddressSanitizer
tree remains blocked by staged Obol headers which predate the required child-list and field
publication APIs; S5 sanitizer qualification remains open. This closes the mesh-LoD
attachment row at the focused C++ boundary and leaves the dated stage estimates unchanged.

### CXX-GED-020 — stream discovery facts and reservation bypass source acceptance

Compact discovery formerly published through three unrelated public source writers. The
expected count grew by a maximum operation rather than remaining an exact fact, the profile
setter accepted internally inconsistent records, and parallel production exposed a
preliminary profile before terminal-empty leaves were resolved. Capacity reservation ran
outside the exact stamped delivery operation and certified the population before allocation,
so a denied reserve could leave partial source state. Source, input and database identity
changes also retained facts discovered for the preceding resource authority.

The occurrence stream now accepts one immutable, internally validated discovery contract.
Its nonzero expected count is exact, its optional profile must agree with that count, and
conflicting or post-cancellation producer writes are rejected. Parallel discovery retains a
producer-local candidate until empty leaves are resolved, then closes the stream with the
final profile and count. `BObolSceneController::certifyDatabaseSourceInstanceCompactStream`
publishes those facts only to the exact current source key and realization stamp. Leaf
capacity reservation is part of that same stamped batch delivery; a denied reserve preserves
the preceding geometry and certified facts. Source, input and database identity invalidation
revokes the contract, while view-only invalidation and a failed delivery retain valid
inventory evidence for the retained preview. Certification and reservation are private
planning state and do not advance a scene frame. The three public source setters are removed.

The direct scene contract covers malformed/conflicting facts, exact certification, stale
stamps, view and resource invalidation, extreme-size certification, denied capacity, retained
geometry and failed-preview evidence. Stream concurrency/cancellation, incremental and
terminal scene publication, source invalidation/configuration allocation sweeps, LoD append
traversal and the live GED denial path also pass. Nine focused checks pass in **26.22
seconds** under `.build/obol-qualification/20260910-ged-stream-contract`, including full GED
draw sync and its CTest row. The saved GED-019 libBObol exposes all three removed setters and
lacks stamped controller certification; the current library has the controller symbol and
none of the split symbols. The AddressSanitizer tree still has stale staged Obol headers, so
S5 sanitizer qualification remains open. This closes the final row in the finite GED/source
presentation writer inventory without changing the dated stage estimates.

### GED presentation writer inventory

This is the finite GED/source presentation sub-inventory for S1/S3. A row is
closed only when its enclosing operation owns target acceptance, retained intent, effective
records and required scene effects.

| Writer family | Current acceptance status |
|---|---|
| Exact/subtree/object visibility, highlight and transparency reducers | **GED-013 qualified:** one scene transaction publishes source aggregate and compact retained/effective state atomically; multiple sources and groups form a callback-safe consistent prefix |
| Global highlight set/clear | **GED-013 qualified:** authored occurrence highlight, retained-rule clearing, source/group state and frame effects share the same transaction |
| Public shape/group-reference display setters | **GED-014 qualified:** canonicalize aliases before one complete source publication; pre-snapshot same-path instances and preserve callback edits, replacement and removal of pending targets; the group branch retains its single prepared publication |
| Redraw and erase realization/publication | **GED-018 qualified:** capture stable source/presentation/realization targets before callbacks; compose exact/prefix removal; retain compact visibility overrides; reject later edits and replacements with exact stamps; keep scoped realization and cleanup on the accepted source |
| Deferred selection initialization (`syncCompactInstanceSelectedPaths`) | **GED-015 qualified:** retain the accepted exact source through selection callbacks and revalidate key, node, routing identity, path and representation before worker capture |
| Compact subtract line style after realization | **GED-016 qualified:** derive default subtract style from durable boolean identity during compact construction and remove the redundant post-realization GED writer |
| Rename and compact-path retargeting | **GED-017 qualified:** prepare repository residency, source identity/hierarchy, compact semantic and renderer records, group paths/names/intents, indexes and peer revisions as one operation; preserve opaque keys and callback successors |
| Mesh-LoD attachment and bounds setters | **GED-019 qualified:** retain the candidate until cut/bounds preparation completes; accept one exact key/stamp; transfer reader and finite bounds together; retire the record on source/input/database invalidation; remove split public/path setters |
| Stream profile, expected-count and reserve setters | **GED-020 qualified:** publish one immutable validated discovery contract through exact stamped source acceptance; reserve inside stamped batch delivery; preserve preceding geometry/facts on denial; revoke facts only when source resource authority changes; remove the three public source setters |

### Group hierarchy writer inventory

These entry points are a finite S1/S3 sub-inventory. An enclosing operation must compose
its membership, indexes and revision effects before notification; calling an individually
qualified helper does not qualify later writes in the caller.

| Entry points | Current acceptance status |
|---|---|
| `ensureGroup` | GROUP-001: complete hierarchy creation, indexes, revision effects and current return target |
| `setGroupDrawIntent`, `setGroupDisplayState` | GED-004: prepared metadata and frame effects |
| `publishDatabaseSourceInstance` | CONFIG-007/GED-004/GED-007/GED-011/GED-012: composed source/membership/metadata/realization publication and qualified live retained/copy callers; GROUP-002 also rejects cycles through prepared descendants |
| `moveDatabaseSourceToGroup`, `moveDatabaseSourceInstanceToGroup` | GROUP-002: prepared membership/index/revision publication, preserved source state, callback acceptance and cycle rejection |
| `renameGroup` | GROUP-003: prepared name/path/index/revision publication, preserved source state and callback acceptance |
| `appendChildToGroup`, `removeChildFromGroup` | GROUP-004: prepared complete subtree membership/index/repository/revision publication, shared-edge survival and callback acceptance |
| `eraseGroupSubpath`, `removeGroup`, `clearGroup` | GROUP-005: composed subtree removal and whole-group clearing, preserved shared ownership, complete path/index/revision effects and callback acceptance |
| `removeDatabaseSource`, `removeDatabaseSourceInstance`, `clearDatabaseSources` | GROUP-006: complete source-edge removal and bulk clearing, selected path identity, duplicate-edge path effects, shared ownership and callback acceptance |
| `moveShapeToGroup`, `removeShape` | GROUP-007: complete parent membership/path/revision publication, retained shape/geometry ownership and callback acceptance |
| `setShapeDrawState`, `setShapeDisplayState`, `setShapePlacementState`, `setShapeSourceState` | GROUP-008: prepared requested scalar fields and frame effects, retained original targets, preserved geometry and callback acceptance |
| `setDatabaseSourceInstanceRealizationState` with role ownership | GED-006: prepared realization record, source roles, owned-shape state and frame effect; later callback edits remain current |
| `subsumeDatabaseSourceInstances` | GED-008: off-scene compact candidate plus composed donor-edge, source-index and revision retirement; retained targets and callback acceptance |
| `mergeDatabaseSourceInstanceCompactOccurrences` | GED-009: stamped target acceptance, consistent-prefix publication, frame effects and full-batch completion witness |
| `adoptDatabaseSourceInstanceRealization` | GED-009: stamped terminal source, stream facts, source-local child/index effects, frame publication and terminal completion witness |
| `setDatabaseSourceInstanceBoundsState` | GED-009 adds a commit witness used by streamed coverage to distinguish preparation from callback failure |
| `publishDatabaseSourceInstanceCompactSnapshot` | GED-010: stamped off-scene compact snapshot, optional certified bounds/profile, terminal state, primary-child retirement, frame effects and commit witness |
| `setSceneRoot`, controller destruction | ROOT-001 qualified: complete root/index/repository/revision publication, callback acceptance, allocation-failure preservation and allocation-free teardown |
| `shareRealizationRepository`, `clearRealizationRepository`, `invalidateRealizationViewVariants`, `renameRealizationObject` | ROOT-001 qualified: prepared membership and cache retirement, affected-key rename, exact per-source residency, retained invocation targets and callback acceptance |

The repository-management entry points retain explicit acceptance rows in this scene
inventory. Their composed root/repository lifetime boundary is qualified below. The
whole-product resource, capability and release matrices remain S4--S6 work.

The adjacent zero-caller GED local mesh/vlist factories were removed under GED-005, the
live source-record bridge is composed under GED-007, resident-occurrence adoption is
composed under GED-008, and live deferred result adapters use stamped scene acceptance
under GED-009. Compact proxy bootstrap is composed under GED-010, and its neighboring
retained-current source sequence is composed under GED-011. Cross-scene source copying is
composed under GED-012, the path/global visibility, highlight and transparency family is
composed under GED-013, public shape/group-reference display is composed under GED-014,
deferred selection target acceptance is composed under GED-015, compact subtract style is
composed under GED-016, database-object rename is composed under GED-017, and redraw/erase
target acceptance and retained visibility are composed under GED-018.
The presentation inventory above names the remaining bounded source/service and enclosing-
operation scope.

### CXX-GROUP-001 — group creation notifies before the hierarchy is complete

The preceding `ensureGroup` appended one group at a time, notifying before indexing it,
then advanced structural/frame revisions after the entire loop. A root observer saw a
partial hierarchy with both a new path and an existing prefix. Its borrowed return pointer
also needed to respect callback removal, rename or replacement.

**The creation boundary is repaired.** `HierarchyPublication` shares path preparation,
retained nodes, child-list replacements and group index preparation with `SourcePublication`.
All missing groups and index entries commit before the scene's revision effect and observer
delivery. Allocation failure before commit preserves the preceding scene. Notification
attempts preserve the first exception; later callback edits remain current. The return
value resolves the canonical requested path after callbacks and can be null. The old live
creation loop and its now-unused `indexSceneGroup` writer are removed.

Two production regression groups cover allocation rejection, existing/new prefixes,
batches, quiet observers, exceptions, recursive no-ops/extensions, removal/replacement,
rename, root replacement and aliased request strings. Both reject the preceding shared
library. [Dated evidence](obol_20260919_readiness_history.md#september-9-complete-group-creation)
owns exact counts, integration/graphical checks and binary identities. Calling other group
mutators from these callbacks exercises creation's acceptance of their completed changes;
it does not qualify those mutators' own notification boundaries.

### CXX-GROUP-002 — source movement publishes intermediate membership

The retained `group_move_probe` observes the preceding
`moveDatabaseSourceInstanceToGroup` with a new destination and an existing one.
Both expose incomplete membership/index/revision state. Creating a destination advances
the frame once before moving the source, then movement advances it again. That implementation
also obtains the source and parent before `ensureGroup` can call observers, but retains the
source only afterward.

**The membership boundary is repaired.** An explicit membership-only constructor reuses
`SourcePublication` preparation, retained targets, child-list/index commit and notification.
It preserves the existing lookup/fallback behavior and source configuration, geometry and
realization identity. The live create/remove/add/reindex sequence is deleted. One complete
structural/frame effect precedes observers, and no later outer mutation overwrites their
edits, moves, removal/replacement, group rename or root replacement.

The cycle concern is now reproduced: the preceding composed source publisher rejected an
existing descendant but accepted a newly prepared group beneath it. `HierarchyPublication`
now checks reachability using prepared child orders where present and live children otherwise.
Both composed source publication and membership-only movement reject existing descendants
and one/two-level new suffixes without changing the scene or revisions.

Three production regression groups cover those six cycle cases, twelve allocation scenarios,
existing/new/root destinations, the path wrapper, no-ops, quiet/throwing observers, batches,
aliased requests and later callback changes. All three reject the preceding shared library.
[Dated evidence](obol_20260919_readiness_history.md#september-9-complete-source-movement)
owns counts and integration/graphical qualification. Other group writers retain their
separate inventory rows; S1--S6 remain open.

### CXX-GROUP-003 — group rename notifies before index and revision effects

**Repaired and qualified at this boundary.** The retained `group_rename_probe` sees partial
path/index/revision state in both callbacks on the preceding implementation. The current
probe sees complete state and one frame increment. Rename shares `GroupPublication` and
prepared indexes, commits all effects before observers and retains original participants.
The old recursive path writer and later index invalidation are removed. Shared groups
prepare once while preserving the previous traversal's final canonical path. Ordinary
scene nodes can remain inside the subtree; rename targets retain BRL scene-group lookup.

Native `SoBase::setName` also lost its old registration on allocation failure. Its private
name dictionaries now prepare destination storage before removing the previous entry,
retire empty registration lists and preserve same-name last-registration semantics.
This enables name assignment as the last allocating step before the prepared scene commit.

`name-registry-transaction`, `group-rename-transaction` and `group-rename-edges` reject the
preceding shared libraries. Production coverage includes native registration/normalization/
clearing, ordinary/shared/nested subtrees, preserved source fields/geometry/realization,
allocation failure, batches, quiet/throwing observers, invalid/colliding names, no-ops,
later/descendant renames, removal/replacement, root replacement and aliased arguments.
The [dated evidence](obol_20260919_readiness_history.md#september-9-complete-group-rename) owns
counts, native suites, independent compilation and integration/graphical qualification.
Other group writers retain their inventory rows; S1--S6 remain open.

### CXX-GROUP-004 — child membership notifies before scene and repository effects

**Repaired and qualified at this boundary.** The retained `membership_probe` exposes
partial scene state on the preceding library for source, group and ordinary-node append
and removal. Direct-source lookups also remain incomplete after the operation. The
current probe observes complete membership, indexes and revisions in all six cases.

`ChildPublication` composes the existing prepared hierarchy and index owners with a
prepared repository update. It indexes nested and auxiliary sources, resolves colliding
paths/instance keys/group keys against the completed graph, and retains shared nodes
still reachable through another edge. Ordinary new insertion prepares the incoming
subtree; collisions and removals inspect remaining reachability without copying the scene.
One structural/frame effect precedes observers. Targets and retired cache values remain
owned through notification, and the publisher performs no later live writes.
The separate child mutation/index sequence and unused direct source-index writer are removed.

Repository seeding previously changed published geometry metadata and could update
reference counts or release old entries before replacement allocation succeeded. Seeding
and release now share a private prepared update: stage affected object counts and cache
entries, commit with map-node transfers, then retire old values after observers. Seeding
shares geometry ownership without reinitializing identity. The common cache-reference
helper also returns its acquired reference if insertion fails.

`group-child-transaction`, `group-child-edges` and `repository-source-update` all reject
the preceding shared library. Coverage includes allocation failure, source/group/plain
and shared children, nested/auxiliary sources, colliding lookup fallback, quiet/throwing
observers, batches, aliased arguments, inverse changes, movement, removal/replacement,
root replacement and cycle rejection. Wire/mesh cache tests cover first seed, replacement,
clearing and shared ownership. [Dated evidence](obol_20260919_readiness_history.md#september-9-complete-child-membership)
owns counts, independent compilation and integration/graphical qualification.
Other removal and root writers retain their separate inventory rows; S1--S6 remain open.

### CXX-GROUP-005 — recursive group removal notifies before complete scene effects

**Repaired and qualified at this boundary.** The retained `removal_probe` observes partial
state for all three entry points with both shared and unshared subtrees on the preceding
library. Each live writer changed children before its revision effect and separately
released repository ownership. Clearing a group also required composition across all its
children, including duplicate edges and sources reachable through several branches.

The existing `ChildPublication` now prepares either a single child edit or a complete
group clear. It retains all affected children, shares the hierarchy/index/repository commit
and resolves surviving source/group lookups against the completed graph. Node and path
observers see complete membership and one structural/frame effect. Old cache values stay
owned through notification; later callback changes remain current. Both path-removal
entry points share one resolver and publisher. The old recursive live removal helper and
separate clear/remove/index/revision sequences are removed.

`group-removal-transaction` and `group-removal-edges` reject the preceding shared library.
Coverage includes allocation failure, retained cache geometry, shared subtrees, duplicate
child edges, nested/auxiliary sources, truncated paths, quiet/throwing observers, batches,
root clearing, no-ops, refill/movement, parent/root replacement and aliased arguments.
Hierarchy-path removal preserves component parsing; indexed clearing preserves its
existing canonical-key lookup and invalid-path behavior. Per-fixture source snapshots
preserve geometry identity; equivalent cross-fixture results compare geometry contents.
[Dated evidence](obol_20260919_readiness_history.md#september-9-complete-group-removal)
owns counts, baseline rejection and integration/graphical qualification.
GROUP-006 below subsequently closes source removal/clearing, and GROUP-008 closes shape
state. Root and repository entry points remain open above.

### CXX-GROUP-006 — source removal notifies before complete scene effects

**Repaired and qualified at this boundary.** The retained `source_removal_probe` observes
partial state in all six preceding-library cases, covering path/instance removal and
bulk clearing with shared and unshared sources. Path removal also retargets a different
node when instance keys collide. The old single-source writer drops surviving lookups
and leaves descendant ownership inconsistent; bulk clearing notifies between parent edits.

The existing `ChildPublication` now prepares the selected source edge or all removable
source edges. Membership, paths, descendant indexes, repository ownership and one
structural/frame effect commit before notification. Ordinary nodes and auxiliary source
subtrees survive bulk clear; explicit auxiliary removal retires its complete subtree.
Shared sources keep remaining edges and lookups, and cache ownership retires only after
the last reachable edge. The return count measures removed physical parent edges.
Path lookup retains its selected source identity through parent resolution instead of
retargeting by a colliding instance key. Single-child edits use prepared removal with
the exact child position; replacement matching cannot distinguish identical siblings.
The obsolete recursive clear helper and incremental source-unindex writers are removed.

`source-removal-transaction` and `source-removal-edges` reject the preceding shared library.
Coverage includes allocation failure, retained shared geometry, duplicate child edges,
shared parent groups, nested/auxiliary sources, path/node observers, quiet/throwing observers,
batches, no-ops, invalid paths, original-target lifetime, colliding keys, callback drain,
reinsertion/replacement, parent/root removal and quiet-policy/argument changes. The outer
operation preserves its result and accepts later callback edits without subsequent writes.
[Dated evidence](obol_20260919_readiness_history.md#september-9-complete-source-removal)
owns exact counts and integration/graphical qualification. GROUP-007/008 below subsequently
close shape membership and scalar state. Root and repository operations remain open;
cross-controller cache ownership is not closed by a shared node within one scene.

### CXX-GROUP-007 — shape membership notifies before complete scene effects

**Repaired and qualified at this boundary.** The retained `shape_probe` exposes partial
state in all four preceding-library cases: mesh/wire movement and removal. Movement
notifies between removing the original edge and adding the destination edge, and both
operations notify before their structural/frame effect. A callback can therefore react
to a detached intermediate shape or obsolete scene revision.

Both entry points now use `ChildPublication`. Movement prepares both parent edits before
commit; removal uses the existing selected-edge transaction. Child paths and one scene
revision effect commit before observers. Shape fields, geometry and unrelated source/cache
ownership remain unchanged. The target stays alive through notification, including when
an observer removes the newly moved shape. Later callback state remains authoritative.
The two direct membership/late-revision sequences are removed without a new state owner.
Existing first-match path lookup, source-subtree exclusion, same-parent/existing-edge
no-ops and invalid-destination behavior remain intact.

`shape-membership-transaction` and `shape-membership-edges` reject the preceding shared
library. Coverage includes mesh/wire shapes, shared geometry, shared nodes and parent
groups, duplicate child positions, allocation failure, node/path observers, quiet/throwing
callbacks, batches, shape style, callback drain/move/reinsertion/replacement, parent/root
removal, aliased requests, colliding paths and target lifetime without an external reference.
Both the new and preceding source-removal fixtures use one independent path-ownership helper.
[Dated evidence](obol_20260919_readiness_history.md#september-9-complete-shape-membership)
owns counts, integration/graphical checks and binary identities. GROUP-008 below closes
the four scalar shape-state setters. Root/repository operations and enclosing GED local
shape callers remain open; membership publication does not qualify their borrowed returns.

### CXX-GROUP-008 — shape state notifies before complete fields and frame effects

**Repaired and qualified at this boundary.** The retained `shape_state_probe` exposes
partial state in all sixteen preceding-library cases, covering four setters, mesh/wire
shapes and observing/editing callbacks. Each setter writes live fields before advancing
the frame revision; a later requested field also overwrites an edit made by the first
callback.

All four entry points now use one `ShapePublication` preparation path and the existing
`PreparedScalarValues`/`PreparedNotifications` machinery. Only requested fields are staged;
geometry arrays, shared geometry and unrelated shape metadata retain their owners. Values
and one frame effect commit before observers. The original shape stays alive through
notification, and no subsequent write overwrites callback state. Structural revisions do
not change. Nullable strings, aliased inputs, existing comparison tolerances, missing-target
results and batch effects are covered. The four typed live writers and their duplicate
mesh/wire dispatch wrappers are removed; the controller is 287 net lines smaller.

`shape-state-transaction` and `shape-state-edges` reject the preceding shared library.
Coverage includes shared shapes/geometry and source caches, allocation failure, node/field/
path observation, quiet fields/nodes, throwing observers, batches, allocation-free tested
no-ops, reentrant updates, removal/movement/replacement, root replacement, notification-policy
changes and target lifetime without an external reference. Preparation failure preserves
the preceding scene; notification failure preserves the complete committed state.
[Dated evidence](obol_20260919_readiness_history.md#september-9-complete-shape-state)
owns exact counts, graphical qualification and binary identities. Root/repository operations
and the enclosing GED local shape callers retain their separate acceptance work. Direct
connected fields/engines and custom cleanup remain inventoried outside these setters.


### CXX-ROOT-001 — root and repository ownership retires before complete scene state

**Repaired and qualified at this boundary.** Four standalone probes against the qualified
GROUP-008 libraries are retained in `.build/obol-qualification/20260909-root-repository`.
They identify five related acceptance failures:

- Replacement and clearing release the old root before publishing root/index/revision
  effects. Both root-deletion callbacks observe obsolete scene state.
- Root replacement has 38 observed allocation positions. Twenty failures preserve the
  preceding scene; eighteen install the next root, release the old cache entry and leave
  revisions unchanged before new cache seeding completes.
- After two controllers share a repository, detaching one root or removing one distinct-
  root source edge evicts geometry still owned through the surviving controller. Both
  same-root and distinct-root detachment reproduce this failure.
- Priming the sibling controller's index before removing a shared-root edge leaves its
  lookup returning the removed source. The unprimed lookup alone masks this defect.
- Destruction calls allocating `setSceneRoot(nullptr)`. Denying its first allocation
  reaches `std::terminate`; the probe's handler verifies `std::bad_alloc` and exits 42.

The ownership boundary must compose root, indexes and cache membership before retirement
can call observers. Source cache refresh is distinct from controller acquisition: realize
also calls `seedSource`, so counting refreshes as new owners is incorrect. Last-owner
release must remain allocation-free during teardown. Shared-parent edits must update or
invalidate every affected controller's lookup before observers, while preserving the
independent single-controller regressions. The working implementation uses native parent
auditor links to identify affected roots without another permanent graph mirror. Child
edits and conflicting-source retirement now share this prepared membership publication.

All four rebuilt reproducers now pass, including all 73 replacement allocation positions
and destruction with no allocation. The shared-controller follow-up reproduces 15 more
failures in rename, group intent, shape state and source display; all 18 follow-up probe
rows now pass, including the three preceding movement successes. One prepared peer-effect
collector composes field and hierarchy changes, preserving each controller's independent
batch and advancing each affected peer once before observers. Source index publication
borrows the enclosing hierarchy instead of maintaining another publisher. Bulk material
refresh prepares the currently accepted target's peers, preserving completed prefixes and
later callback changes; the old owner-only configuration callback is removed.

Three new production checks cover shared allocation failure, nested/shared roots,
independent repositories and batches, quiet/no-op behavior, observer exceptions, peer
replacement/destruction, later edits, source conflicts and per-target material effects.
The realization follow-up reproduces missing peer revisions, obsolete traversal after root
detachment and crashes after source removal or repository replacement. Captured repository,
root and source references now survive callbacks. Prepared peer effects compose with the
existing view effects, and hierarchy/repository changes stop every enclosing native
traversal even inside mutation batches. Completed progress remains observable; a subsequent
call visits the current scene. Two production checks cover 64 callback histories, eight
quiet/shared source-root cases, sole-owner root retirement and 189 allocation positions.
An expanded peer-destruction history also exposed stale entries in native `SoBase::getAuditors`:
its shrinking-list loop skipped retired sensors. Removing from the end fixes this; a native
regression rejects the preceding library.

The source-descendant follow-up treats successful realization child replacement as one
scene publication. Exact replacement order, descendant indexes, repository membership and
peer revisions are prepared before source writes and commit before observers. Controllers
rooted at removed descendants observe complete structural state; retained auxiliary roots
receive a frame effect only when their own metadata changes. The compact-hierarchy cleanup
retains exact source and ancestor nodes, rejects ambiguous keys and revalidates current
identity, parent links and compact ancestry before every removal. Callback-time replacement,
reparenting and root detachment therefore supersede the original cleanup targets.

Two production selectors cover 12 shared/split/independent repository histories, eight
auxiliary cases, one mixed success/failure sequence, six throwing/nonthrowing cleanup
histories, 296 descendant-publication allocation positions and 248 cleanup allocation
positions. The mixed sequence rejects an intermediate implementation that replayed the
preceding source's prepared hierarchy effects after a later source failed. Both selectors
reject the original preceding working libBObol. The complete working suites at that step
contained 38 source checks plus 21 source scenarios and 66 scene checks; 23 rebuilt probe
binaries provided 29 passing rows.

The late-source follow-up applies the same boundary to auxiliary source insertion,
replacement and removal and to external line, point, mesh, annotation, clear,
primitive-wire and submodel publication. Exact child orders, indexes, repository membership
and every affected controller revision prepare before source mutation and commit before
callbacks. The external selector covers seven callback histories and 494 allocation
positions (492 preceding states and two complete commits); the auxiliary selector covers
eight removal histories and 281 allocation positions (268 preceding states and 13 complete
commits). Both reject the saved immediately preceding libBObol. The complete working suites
now contain 40 source checks plus 21 source scenarios and 66 scene checks; six independent
translation units and all 29 rebuilt probe rows pass.

Production integration further distinguishes representation replacement from topology.
An ordinary realization advances affected frames while leaving structural revisions
unchanged. The prepared membership diff promotes the publication to structural when it
adds or removes a nested database source. Attached controller indexes are authoritative:
source and named-group membership changes go through the scene-controller transaction,
while fixture construction may finish raw hierarchy assembly before attachment.

The repository-operation follow-up prepares cache clear, object/view invalidation and
rename before releasing values. Clear and invalidation first extract every affected entry;
rename stages only exact/prefix keys and affected source-owner records, so preparation scales
with the renamed object rather than the full cache. Commit performs no allocation. Residency
counts the union of old/new owners once per source, preventing an overlapping owner from
leaking the final cache lease. The controller retains the invoked repository and stops an
active realization pass before mutation, so a destruction callback can replace the scene's
repository without destroying the active callee or having its later edit overwritten.

Two production selectors cover four direct cache callbacks, three scene repository-
replacement callbacks, allocation-free clear/view invalidation, eight rename allocation
positions and exact last-owner retirement. A third selector expands realization acceptance
to 88 callback histories. The saved preceding library rejects the cache and realization
selectors and segfaults in the isolated scene wrapper selector. Current source publication
passes 42 direct checks plus 21 scenarios in 113.23 seconds; the unchanged 66-check scene
suite passes in 179.04 seconds. After rebuilding static libBObol and GED, synchronization
passes in 7.89 seconds, six translation units compile independently and 23 rebuilt probe
binaries provide 29 passing rows. The renderer smoke test passes with OSMesa still
`097aac27...`.

The final integration evidence records exact results and limits. The saved GROUP-008
library rejects all three integrated root/repository production selectors; the current
library passes them. The complete source suite passes 45 direct checks plus 21 scenarios,
and the state suite passes 66 checks. All 22 focused checks, six independent translation
units, 29 rows from 23 freshly linked probes, the isolated full build, the five-client
target set, GED synchronization and the renderer smoke test pass. Native Obol passes 1,577
unit tests with one existing skip and 63 integration tests.

Fresh graphical qualification provides 24 passing rows across eager drawing, delivery and
resource denial, Generic Twin cold/warm shaded/wire operation and primitive editing on
OSMesa and private-Xvfb System GL. Twenty-four selected originals were inspected. The
focused records and independent verifier are under
`.build/obol-qualification/20260909-root-repository`; graphical reports and images are under
`.build/obol-qualification/20260910-root-integration-2`. This closes ROOT-001. GED local
shape caller acceptance and source/service/presentation ownership remain open, along with
S1--S6 production qualification.
