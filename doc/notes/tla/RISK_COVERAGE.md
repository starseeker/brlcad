# Historical and prospective risk coverage

This matrix was derived from `libbobol_engineering_lessons.md`, the active-debt
and production-readiness records, implementation comments, regression names,
and retained TLC coverage.  Past failures are seeds, not the boundary of the
audit: the “prospective variants” column names adjacent failures which the
same mechanism can plausibly produce.

Status meanings:

- **checked**: a TLA+ invariant/property and an executable C++ witness exist;
- **formal added**: the control risk is modeled, but the listed refinement
  regression must still be made explicit;
- **implementation audit**: TLA+ states the contract, while code inspection,
  fault injection, sanitizers, or stress must establish the data-plane edge;
- **outside TLA+**: retain executable, graphical, or performance evidence.

| Risk family | Historical seed | Prospective variants to preempt | Formal coverage | Executable/refinement status |
|---|---|---|---|---|
| Result identity and ABA reuse | Raw-address reuse, stale dense index, reused occurrence key | Late result with the right route but wrong population; right population but old demand; replacement during completion; revision wrap | Three-domain `ObolResultAuthentication`, `ObolOccurrenceControl`, `ObolSparsePlanIdentity` | **formal added** across independent population, demand, and route revisions; `test_view_controller_compact_direct_result_route` covers a retired compact population, but the combined three-domain C++ overtake matrix remains. Counter wrap is an implementation audit because the formal revisions are monotonic naturals. |
| Complete evidence stamps | Plans and frames accepted with companion fields from different epochs | Visibility-only mutation, policy-only mutation, sequential multi-domain update observed halfway through, inactive input accidentally entering identity | Canonical six-domain stamp in `ObolProgressivePipeline`; `ObolCadFrameComposition`; `ObolSemanticPresentation` | **implementation audit**. Visibility and stale-frame paths have executable witnesses, but every plan/frame constructor must still be shown to populate all six domains exactly once. |
| Shared producer lifetime | Asset demand confused with view demand; shared layers treated as exclusive | One view closes/cancels while another still needs the build; late subscriber joins; final consumer leaves during result publication | `ObolActiveProducerDemand`, new `ObolSharedAssetLease`, `ObolAssetPublicationComposition` | **formal added**, including bounded late join during build/result/post-cancel. `test_asset_producer_coalescing` covers duplicate siblings, but close-one-consumer-during-build and close-final-consumer-during-result remain. |
| Typed completion outcomes | Cache failure treated as invalid geometry; cancellation treated as failure; supersession lost | Retryable transport failure, terminal invalid source, cancellation after result creation, failure for an obsolete demand | `ObolActiveProducerDemand`, `ObolOccurrenceControl`, `ObolPresentationPreparation`, `ObolSpatialProducer`, typed asset composition | **formal added** for success, retry, constraint, failure, cancellation, and supersession; one table-driven C++ result-routing test must still cross every outcome with current and stale demand. |
| Partial publication and transaction failure | Partial BREP success, mutable shared alias, scene visible inside a batch | Allocation/resource failure after transaction open, exception during notification, mixed old/new sparse plan, reentrant observer during close | `ObolCadMutation`, `ObolCadFrameComposition`, `ObolLiveSpatialPublication` | **implementation audit**. The formal contract requires complete preallocation and a no-throw swap; add allocator fault injection before every commit boundary and assert the prior scene/notification pair remains unchanged. |
| Progress owner without producer | Empty cursor, zero-work repair, generic completion erasing a certificate | Cursor selected after its range was retired, wake notification consumed before owner publication, owner recreated without rank decrease | `ObolControlRefinement`, `ObolSubmissionPass`, `ObolCompletedPassOwnership`, corrected `ObolLodConvergence`, corrected `ObolOccurrenceControl` | **checked** by seeded coordinator and empty-drain regressions; required-action coverage now rejects an unreachable empty-tail or repair branch. |
| Starvation and unfair scheduling | Policy waves starved bounded work; one completion selected multiple owners | Two large assets alternate without completion, timer work starved by frames, optional refinement overtakes protected floor | `ObolProgressivePipeline.MultiLargeAssetsDoNotStarve`, `ObolLodArbitration`, `ObolDeadlineOwnership`, `ObolTerminalQualityOrdering` | **checked** at bounded control level. Add long randomized multi-asset stress under worker-count 1 and verify per-asset completion ranks, not just global activity. |
| Capacity and quality evidence | Candidate confused with exact framebuffer, first frame confused with steady state, ceiling erased | Context loss after measurement, equal render cost but richer population, resize/DPR change with retained capacity, non-CAD frames contaminating evidence | `ObolCapacitySearch`, `ObolCadTimingEvidence`, `ObolStaticQuality`, `ObolPointCalibrationDomains`, terminal compositions | **checked** for core transitions; graphical/performance matrices remain authoritative for actual cost. Add endpoint-loss invalidation of all frame-derived capacity certificates. |
| Resource-ledger atomicity | Growth/allocation split, stale allocation causing no-op retry, pressure aggregate lost extents | Overflow in byte/count sums, partial reservation rollback, consumer removal during compaction reservation, eviction racing renewed use | `ObolResidentGrowth`, `ObolResidentCompaction`, `ObolRetainedAllocationPresentation` | **implementation audit** with service stress and saturation tests. Run TSAN for consumer/compaction races and add reservation-failure injection. |
| Cache durability and cold publication | Destructive cache status, late preview, partial hierarchy discovery | Truncated descriptor, disk-full rename, foreign version/endianness, reader during writer crash, stale bounds certifying coverage | `ObolSpatialProducer`, `ObolLiveSpatialPublication`, `ObolLodColdPreview`, asset composition | **checked** for atomic descriptor, cancellation, stale cache, and write failure; crash consistency and filesystem atomicity are **outside TLA+** and need process-kill/fault-injection tests. |
| Interaction and teardown | Host callback lost wake, worker destruction order, callback storage lifetime | Close while gesture remains open, endpoint detach during in-flight frame, unsubscribe from callback, static destruction with queued work | `ObolInteractionSession`, `ObolHostWork`, `ObolControlLifecycleComposition`, `ObolSharedAssetLease` | Callback/unsubscribe and source-realization teardown tests exist. Add close-during-gesture and endpoint-loss-during-exact-frame refinement cases; keep ASan/TSan gates. |
| Determinism and ordering | Deadline selected multiple successors; unordered evidence/cursor conflation | Equal-score tie depends on pointer/hash order, replay differs across process runs, two owners accept one frame | `ObolQuietSuccessor`, `ObolCompletedPassOwnership`, `ObolTerminalQualityOrdering`, `ObolLodArbitration` | Oracle and cross-run tests cover major paths. Require stable identity as the last comparator for every unordered candidate set. |
| Geometry truth and numeric robustness | Partial tessellation published as success, structural points confused with richer geometry | Empty-valid geometry vs invalid, NaN/Inf bounds, degenerate extents, count overflow, mixed primitive coverage | TLA+ models only typed success/failure and structural coverage | **outside TLA+**. Preserve geometry validation, fuzzing, large-count fixtures, and image qualification; never promote a control proof into a geometry-validity claim. |
| Renderer/GUI semantics | Logical/device pixels confused, repaint treated as semantic render, stale native geometry | DPR change mid-frame, selection-only mutation reopening LoD, hidden canvas first-paint, context recreation retaining old report | `ObolSemanticPresentation`, `ObolCadFrameComposition`, lifecycle models | **outside TLA+** for pixel correctness; semantic model plus headless/GUI image matrices and endpoint lifecycle tests are jointly required. |

## New roadmap obligations

The audit makes these near-term C++ checks explicit:

1. Add shared-producer tests where one consumer closes during build and the
   final consumer closes while a result is queued.
2. Add a combined stale-result matrix over population epoch, demand/view
   revision, and source route; exactly one mismatch must force supersession.
3. Add allocator/resource fault injection at retained-scene and presentation
   preparation commit boundaries.
4. Add endpoint-loss and close-during-gesture tests which prove every fair
   frame/timer owner receives a cancellation or failure completion.
5. Audit every six-domain evidence-stamp constructor and reducer write; no
   defaulted companion field may create an unrecorded semantic transition.
6. Keep geometry, numeric, image, performance, and concurrency qualification
   as independent gates; they are refinement evidence, not consequences of
   TLC success.
7. Audit revision increment and comparison code for overflow.  Long-lived
   identity counters must saturate, widen, or deliberately invalidate all
   dependent evidence before reuse; unsigned wrap must never authenticate an
   ABA result.
