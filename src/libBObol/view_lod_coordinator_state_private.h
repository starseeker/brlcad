/*    V I E W _ L O D _ C O O R D I N A T O R _ S T A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_VIEW_LOD_COORDINATOR_STATE_PRIVATE_H
#define LIBBOBOL_VIEW_LOD_COORDINATOR_STATE_PRIVATE_H

#include "BObol/BPickDetail.h"
#include "BObol/BLodService.h"
#include "BObol/BViewController.h"
#include "lod_control_private.h"
#include "lod_coordinator_private.h"
#include "retained_allocation_private.h"

#include <algorithm>
#include <cmath>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <optional>
#include <unordered_map>
#include <utility>
#include <vector>

/* One pinned dense-entry plan for the source currently owned by the bounded
 * submission cursor.  Source identity is the validity witness; a separate
 * boolean can otherwise outlive either the source or its entry vector. */
class BObolLodSubmissionSourcePlan {
public:
    void begin(SoBRLDatabaseSource *source)
    {
	this->sourceValue = source;
	this->entryValues.clear();
    }

    void assign(SoBRLDatabaseSource *source,
	const std::vector<size_t> &entries)
    {
	this->sourceValue = source;
	this->entryValues = entries;
    }

    void reset(void)
    {
	this->sourceValue = NULL;
	this->entryValues.clear();
    }

    bool valid(void) const
    {
	return this->sourceValue != NULL;
    }

    bool validFor(const SoBRLDatabaseSource *source) const
    {
	return source && this->sourceValue == source;
    }

    std::vector<size_t> &entries(void)
    {
	return this->entryValues;
    }

    const std::vector<size_t> &entries(void) const
    {
	return this->entryValues;
    }

    size_t size(void) const
    {
	return this->entryValues.size();
    }

private:
    SoBRLDatabaseSource *sourceValue = NULL;
    std::vector<size_t> entryValues;
};

/* Accumulator for a possibly time-sliced non-compact visibility census.  The
 * count and its source identity are one value so a source change cannot reuse
 * the preceding source's partial count. */
class BObolLodSubmissionVisibleCount {
public:
    void observe(SoBRLDatabaseSource *source, size_t visible)
    {
	if (this->sourceValue != source) {
	    this->sourceValue = source;
	    this->countValue = 0;
	}
	this->countValue = visible > SIZE_MAX - this->countValue ?
	    SIZE_MAX : this->countValue + visible;
    }

    void reset(void)
    {
	this->sourceValue = NULL;
	this->countValue = 0;
    }

    size_t valueFor(const SoBRLDatabaseSource *source) const
    {
	return source && source == this->sourceValue ? this->countValue : 0;
    }

private:
    SoBRLDatabaseSource *sourceValue = NULL;
    size_t countValue = 0;
};

/* Exact source-inventory delta consumed by one bounded submission pass.  A
 * targeted source without a selective entry plan requests a full source scan;
 * a targeted source with a plan requests only those dense entries.  Activity
 * is derived from the target population so it cannot outlive the work it is
 * supposed to describe. */
class BObolLodSubmissionDelta {
public:
    using SourcePlan =
	std::pair<SoBRLDatabaseSource *, std::vector<size_t>>;

    void reset(void)
    {
	this->sourceValues.clear();
	this->planValues.clear();
    }

    bool active(void) const
    {
	return !this->sourceValues.empty();
    }

    bool targets(const SoBRLDatabaseSource *source) const
    {
	return source && std::find(this->sourceValues.begin(),
		this->sourceValues.end(), source) != this->sourceValues.end();
    }

    bool target(SoBRLDatabaseSource *source)
    {
	if (!source || this->targets(source))
	    return false;
	this->sourceValues.push_back(source);
	return true;
    }

    bool targetSelective(SoBRLDatabaseSource *source,
	std::vector<size_t> entries)
    {
	if (!this->target(source))
	    return false;
	this->planValues.emplace_back(source, std::move(entries));
	return true;
    }

    void targetFull(SoBRLDatabaseSource *source)
    {
	(void)this->target(source);
	this->removeSelective(source);
    }

    void replaceSelectivePlans(std::vector<SourcePlan> plans)
    {
	this->reset();
	for (SourcePlan &plan : plans) {
	    std::vector<size_t> *existing =
		this->selectiveEntries(plan.first);
	    if (existing) {
		existing->insert(existing->end(),
		    plan.second.begin(), plan.second.end());
		continue;
	    }
	    (void)this->targetSelective(plan.first, std::move(plan.second));
	}
    }

    std::vector<size_t> *selectiveEntries(
	const SoBRLDatabaseSource *source)
    {
	auto plan = std::find_if(this->planValues.begin(),
		this->planValues.end(), [source](const SourcePlan &entry) {
		    return entry.first == source;
		});
	return plan == this->planValues.end() ? NULL : &plan->second;
    }

    const std::vector<size_t> *selectiveEntries(
	const SoBRLDatabaseSource *source) const
    {
	auto plan = std::find_if(this->planValues.begin(),
		this->planValues.end(), [source](const SourcePlan &entry) {
		    return entry.first == source;
		});
	return plan == this->planValues.end() ? NULL : &plan->second;
    }

    void removeSelective(const SoBRLDatabaseSource *source)
    {
	this->planValues.erase(std::remove_if(this->planValues.begin(),
		this->planValues.end(), [source](const SourcePlan &entry) {
		    return entry.first == source;
		}), this->planValues.end());
    }

    const std::vector<SourcePlan> &plans(void) const
    {
	return this->planValues;
    }

    size_t planCount(void) const
    {
	return this->planValues.size();
    }

private:
    std::vector<SoBRLDatabaseSource *> sourceValues;
    std::vector<SourcePlan> planValues;
};

/**
 * Owner-thread progressive-display state.  This first extraction preserves
 * the exact field layout and algorithms formerly embedded in Impl; it adds
 * no allocation, virtual dispatch, or per-occurrence indirection.
 */
struct BObolLodCoordinator {
    uint64_t residentMeshConsumerId(void) const
    {
	return static_cast<uint64_t>(
	    reinterpret_cast<uintptr_t>(this));
    }

    void invalidateResidentMeshCompactionSnapshot(void)
    {
	if (this->lodService)
	    this->lodService->invalidateResidentMeshConsumer(
		this->residentMeshConsumerId());
    }

    BObolLodAdmissionRevisionStamp admissionRevisionStamp(void) const
    {
	BObolLodAdmissionRevisionStamp stamp;
	stamp.inventory = this->lodAdmissionInventoryRevision;
	stamp.availability = this->lodAdmissionAvailabilityRevision;
	stamp.view = this->lodViewRevision;
	stamp.policy = this->lodPolicyRevision;
	stamp.capacity = this->lodAdmissionCapacityRevision;
	return stamp;
    }

    bool retainedAllocationCertificateCurrent(
	const BObolViewLodState *state) const
    {
	const BObolRetainedAllocationResult &allocation =
	    this->lodRetainedAllocationCertificate;
	return state && allocation.currentPlanSerial(
		this->admissionRevisionStamp(),
		state->activeCadAllocationPlan()) != 0 &&
	    std::isfinite(allocation.pointProxyPixelThreshold) &&
	    std::fabs(allocation.pointProxyPixelThreshold -
		this->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f &&
	    state->cadAllocationPlanCoversCurrentPopulation(
		allocation.allocationPlanSerial, allocation.viewRevision(),
		allocation.policyRevision(),
		allocation.fixedCadPresentationCost);
    }

    bool retainedAllocationRepresentsCurrentView(void) const
    {
	const BObolRetainedAllocationResult &allocation =
	    this->lodRetainedAllocationCertificate;
	if (allocation.selectedPresentationCost > 0)
	    return true;
	/* Retained payloads outside the frustum remain bound so a later camera
	 * change can reuse their immutable PoP data.  Consequently an exact empty
	 * view has a valid zero-cost allocation even though the view state still
	 * owns nonzero resident payloads.  Accept zero only with both halves of the
	 * proof: the allocator found no current demand and the completed dense
	 * visibility census found no visible occurrence.  This keeps a producer
	 * failure which accidentally omitted visible payloads from masquerading as
	 * a completed empty view. */
	return allocation.pixelDemandPresentationCost == 0 &&
	    this->lodCoveragePolicy.hasCompleteVisibleCount() &&
	    this->lodConvergenceCandidateCount() == 0;
    }

    bool retainedAllocationCutsApplied(
	const BObolViewLodState *state) const
    {
	if (!this->retainedAllocationCertificateCurrent(state))
	    return false;
	const BObolRetainedAllocationResult &allocation =
	    this->lodRetainedAllocationCertificate;
	if (!this->retainedAllocationRepresentsCurrentView())
	    return false;
	return state->cadAllocationPlanCutsApplied(
		allocation.allocationPlanSerial, allocation.viewRevision(),
		allocation.policyRevision(),
		allocation.fixedCadPresentationCost);
    }

    bool retainedAllocationPresentationRealized(
	const BObolViewLodState *state) const
    {
	if (!this->retainedAllocationCertificateCurrent(state))
	    return false;
	if (!this->retainedAllocationRepresentsCurrentView())
	    return false;
	const BObolRetainedAllocationResult &allocation =
	    this->lodRetainedAllocationCertificate;
	if (this->lodInteractiveProgressiveCeiling < 0)
	    return this->retainedAllocationCutsApplied(state);
	/* A global cut is occurrence-local only for one progressive payload.
	 * In that domain the framebuffer contains the canonical ceiling cost,
	 * while activeRenderCost() intentionally describes the richer retained
	 * population hidden behind it. */
	if (this->lodInteractiveProgressiveCeiling >= 0 &&
	    state->cadProgressivePayloadCount() == 1) {
	    const size_t ceilingCadCost = state->
		cadRenderCostAtProgressiveCutCeiling(
		    this->lodInteractiveProgressiveCeiling);
	    const size_t presentationCost = BObolLodAdmissionPlanner::
		canonicalSceneCostAtCadCeiling(
		state->activeRenderCost(), state->activeCadRenderCost(),
		ceilingCadCost);
	    return presentationCost == allocation.selectedPresentationCost;
	}
	/* With several independently allocated occurrences, matching aggregate
	 * cost would not prove matching visual distribution.  Only an inert global
	 * ceiling may coexist with the applied occurrence-local plan. */
	return this->retainedAllocationCutsApplied(state) &&
	    state->maximumActiveProgressiveCut() <=
		this->lodInteractiveProgressiveCeiling;
    }

    void advanceAdmissionRevision(BObolLodAdmissionRevisionDomain domain)
    {
	/* A semantic input edge may retarget resident cuts without provider I/O.
	 * Cancel old quiet-memory work before exposing that new demand; compaction
	 * is subordinate to the admission plan, never another presentation owner. */
	this->invalidateResidentMeshCompactionSnapshot();
	const BObolLodAdmissionRevisionStamp next =
	    BObolLodRevisionContract::advance(
		this->admissionRevisionStamp(), domain);
	this->lodAdmissionInventoryRevision = next.inventory;
	this->lodAdmissionAvailabilityRevision = next.availability;
	this->lodViewRevision = next.view;
	this->lodPolicyRevision = next.policy;
	this->lodAdmissionCapacityRevision = next.capacity;
	/* No caller may observe or consume progress certified by the preceding
	 * tuple after this semantic edge. */
	this->lodAdmissionCursor.reset();
    }

    void setAdmissionPolicyRevision(uint64_t revision)
    {
	if (this->lodPolicyRevision.value() != revision)
	    this->invalidateResidentMeshCompactionSnapshot();
	const BObolLodAdmissionRevisionStamp next =
	    BObolLodRevisionContract::setPolicy(
		this->admissionRevisionStamp(), revision);
	this->lodPolicyRevision = next.policy;
	/* An externally supplied policy epoch has the same invalidation contract
	 * as an incremented one.  Keeping this reset here prevents public policy
	 * synchronization from becoming a sixth revision owner. */
	this->lodAdmissionCursor.reset();
    }

    bool commitAdmissionPlan(const BObolLodAdmissionPlan &plan)
    {
	if (BObolLodRevisionContract::planDisposition(
		plan.revisionStamp, this->admissionRevisionStamp()) ==
	    BObolLodRevisionContract::PlanDisposition::STALE)
	    return false;
	this->lodAdmissionEvidence = plan.nextEvidence;
	this->lodAdmissionCursor = plan.nextCursor;
	return true;
    }

    void applyAdmissionEvidenceAction(
	BObolLodAdmissionPlanner::EvidenceAction action)
    {
	this->commitAdmissionPlan(BObolLodAdmissionPlanner::applyEvidenceAction(
	    this->lodAdmissionEvidence, this->lodAdmissionCursor, action));
    }

    void requestCapacityRescan(void)
    {
	this->commitAdmissionPlan(BObolLodAdmissionPlanner::requestCapacityRescan(
	    this->lodAdmissionEvidence, this->lodAdmissionCursor));
    }

    void completeAppliedAllocation(
	const BObolLodCapacityEvidence::CompletedAllocationInputs &inputs)
    {
	this->commitAdmissionPlan(
	    BObolLodAdmissionPlanner::completeAppliedAllocation(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, inputs));
    }

    void requestRetainedRecovery(size_t budget)
    {
	this->commitAdmissionPlan(BObolLodAdmissionPlanner::requestRetainedRecovery(
	    this->lodAdmissionEvidence, this->lodAdmissionCursor, budget));
    }

    void requestRetainedNormalization(size_t budget)
    {
	this->commitAdmissionPlan(
	    BObolLodAdmissionPlanner::requestRetainedNormalization(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, budget));
    }

    void requestRetainedReallocation(bool preserveCurrentBudget = true)
    {
	this->commitAdmissionPlan(
	    BObolLodAdmissionPlanner::requestRetainedReallocation(
		this->lodAdmissionEvidence, this->lodAdmissionCursor,
		preserveCurrentBudget));
    }

    bool resumeCapacityCandidateAllocation(void)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::resumeCapacityCandidateAllocation(
		this->lodAdmissionEvidence, this->lodAdmissionCursor);
	if (!this->commitAdmissionPlan(plan))
	    return false;
	return plan.transitionChanged;
    }

    void requestPresentationReconciliation(size_t budget)
    {
	this->commitAdmissionPlan(
	    BObolLodAdmissionPlanner::requestPresentationReconciliation(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, budget));
    }

    size_t requestCoverageCompletion(size_t activeCost,
	size_t certifiedBudget)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::requestCoverageCompletion(
		this->lodAdmissionEvidence, this->lodAdmissionCursor,
		activeCost, certifiedBudget);
	this->commitAdmissionPlan(plan);
	return plan.transitionValue;
    }

    BObolLodCapacityEvidence::CompletedFrameDecision
    recordDeadlineCapacityMiss(
	const BObolLodCapacityEvidence::DeadlineMissInputs &inputs)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::recordDeadlineCapacityMiss(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, inputs);
	this->commitAdmissionPlan(plan);
	return plan.completedFrameDecision;
    }

    void setRetainedQualityFloor(size_t budget,
	uint64_t populationSignature, size_t activeCost,
	size_t minimumActiveCost)
    {
	this->commitAdmissionPlan(
	    BObolLodAdmissionPlanner::setRetainedQualityFloor(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, budget,
		populationSignature, activeCost, minimumActiveCost));
    }

    bool recordRetainedQualityFloorMiss(void)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::recordRetainedQualityFloorMiss(
		this->lodAdmissionEvidence, this->lodAdmissionCursor);
	this->commitAdmissionPlan(plan);
	return plan.transitionChanged;
    }

    bool rejectRetainedQualityFloor(void)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::rejectRetainedQualityFloor(
		this->lodAdmissionEvidence, this->lodAdmissionCursor);
	this->commitAdmissionPlan(plan);
	return plan.transitionChanged;
    }

    void recordRetainedQualityFloorMet(bool exactProtectedPopulation,
	uint64_t populationSignature, size_t presentedCost)
    {
	this->commitAdmissionPlan(
	    BObolLodAdmissionPlanner::recordRetainedQualityFloorMet(
		this->lodAdmissionEvidence, this->lodAdmissionCursor,
		exactProtectedPopulation, populationSignature, presentedCost));
    }

    bool confirmRetainedRecoveryPresentation(bool onePixelReady)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::confirmRetainedRecoveryPresentation(
		this->lodAdmissionEvidence, this->lodAdmissionCursor,
		onePixelReady);
	this->commitAdmissionPlan(plan);
	return plan.transitionChanged;
    }

    void raiseCapacityBudget(size_t budget)
    {
	this->commitAdmissionPlan(BObolLodAdmissionPlanner::raiseCapacityBudget(
	    this->lodAdmissionEvidence, this->lodAdmissionCursor, budget));
    }

    BObolLodCapacityEvidence::CalibrationDecision finishBlockedCapacityPass(
	const BObolLodCapacityEvidence::CalibrationInputs &inputs)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::finishBlockedCapacityPass(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, inputs);
	this->commitAdmissionPlan(plan);
	return plan.calibrationDecision;
    }

    BObolLodCapacityEvidence::CompletedFrameDecision
    completeCapacityCalibrationFrame(void)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::completeCapacityCalibrationFrame(
		this->lodAdmissionEvidence, this->lodAdmissionCursor);
	this->commitAdmissionPlan(plan);
	return plan.completedFrameDecision;
    }

    BObolLodCapacityEvidence::CompletedFrameDecision
    completeCapacitySearchFrame(
	const BObolLodCapacityEvidence::CompletedFrameInputs &inputs)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::completeCapacitySearchFrame(
		this->lodAdmissionEvidence, this->lodAdmissionCursor, inputs);
	this->commitAdmissionPlan(plan);
	return plan.completedFrameDecision;
    }

    BObolLodCapacityEvidence::CompletedFrameDecision
    retireUnmeasurableCapacityFrame(void)
    {
	const BObolLodAdmissionPlan plan =
	    BObolLodAdmissionPlanner::retireUnmeasurableCapacityFrame(
		this->lodAdmissionEvidence, this->lodAdmissionCursor);
	this->commitAdmissionPlan(plan);
	return plan.completedFrameDecision;
    }

    static constexpr uint64_t defaultInteractivePresentationDeadline(void)
    {
	return 40000000ULL;
    }

    static constexpr uint64_t defaultStablePresentationDeadline(void)
    {
	return 100000000ULL;
    }

    static constexpr uint64_t maximumInteractivePresentationDeadline(void)
    {
	return 100000000ULL;
    }

    static constexpr uint64_t maximumStablePresentationDeadline(void)
    {
	return 250000000ULL;
    }

    static constexpr uint64_t defaultProminentQualityDeadline(void)
    {
	/* A quiet one-shot frame is cancellable between software submission
	 * chunks and is retained without redraw.  This allowance lets a prominent
	 * discrete mesh reach its recognizable neighboring cut without weakening
	 * the ordinary 250 ms stable-interaction endpoint. */
	return 400000000ULL;
    }

    static uint64_t targetPresentationDeadline(float targetFps,
	uint64_t defaultDeadline, uint64_t maximumDeadline)
    {
	if (!std::isfinite(targetFps) || targetFps <= 0.0f)
	    return defaultDeadline;
	const long double frameNanoseconds =
	    1000000000.0L / static_cast<long double>(targetFps);
	if (!std::isfinite(frameNanoseconds) || frameNanoseconds <= 0.0L)
	    return defaultDeadline;
	const uint64_t targetDeadline = frameNanoseconds >=
		static_cast<long double>(UINT64_MAX) ? UINT64_MAX :
	    static_cast<uint64_t>(frameNanoseconds);
	return std::min(maximumDeadline,
	    std::max(defaultDeadline, targetDeadline));
    }

    static float boundedPresentationTargetFps(float targetFps,
	uint64_t maximumDeadline)
    {
	if (!maximumDeadline)
	    return targetFps;
	const long double minimumFps = 1000000000.0L /
	    static_cast<long double>(maximumDeadline);
	return std::max(targetFps, static_cast<float>(minimumFps));
    }

    uint64_t prominentQualityPresentationDeadline(void) const
    {
	return std::max(stablePresentationFrameDeadlineNanoseconds,
	    prominentQualityFrameDeadlineNanoseconds);
    }

    float prominentQualityTargetFps(void) const
    {
	const uint64_t deadline = prominentQualityPresentationDeadline();
	if (!deadline)
	    return lodStableTargetFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(deadline);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return lodStableTargetFps;
	const float bounded = static_cast<float>(deadlineFps);
	return lodStableTargetFps > 0.0f ?
	    std::min(lodStableTargetFps, bounded) : bounded;
    }

    uint64_t staticQualityPresentationDeadline(void) const
    {
	/* Static overscan is the one-shot proof for quality which exceeds the
	 * ordinary refinement cadence.  Requiring the allocator's all-or-nothing
	 * prominent floor to be active here is circular for a discrete mesh: the
	 * floor cannot become active until its complete cut already fits the very
	 * trial this deadline enables.  Give the whole bounded static phase the
	 * longer allowance while the preceding framebuffer remains available.
	 * A rejected rich cut retains this allowance after its bounded handoff:
	 * REJECTED is the terminal certificate that the occurrence-local
	 * population is the richest one proven under this event-driven deadline.
	 * Returning that same population to the ordinary quiet deadline immediately
	 * reopens deadline recovery and makes convergence depend on timing jitter.
	 * A view, policy, renderer, or resource epoch explicitly invalidates the
	 * certificate and restores the ordinary cadence. */
	const BObolLodCapacitySearchCertificate &capacitySearch =
	    lodAdmissionEvidence.capacity().capacitySearch();
	const bool capacitySearchUsesStaticDeadline =
	    capacitySearch.phase() !=
		BObolLodCapacitySearchCertificate::Phase::INACTIVE &&
	    capacitySearch.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC &&
	    capacitySearch.activeTargetNanoseconds() >
		capacitySearch.key().preferredTargetNanoseconds;
	if (lodStaticQualityTrial.usesStaticDeadline() ||
	    lodPointQualityPhase.staticCalibrationPending() ||
	    capacitySearchUsesStaticDeadline)
	    return prominentQualityPresentationDeadline();
	return stablePresentationFrameDeadlineNanoseconds;
    }

    void resetDeadlineSafePresentation(void)
    {
	lodDeadlineSafePresentation.reset();
    }

    void resetRetainedAdmissionQualityProof(void)
    {
	lodRetainedAdmissionQualityEvidence.reset();
	lodRetainedAllocationCertificate = BObolRetainedAllocationResult();
    }

    void resetLodViewQualityHistory(void)
    {
	lodViewQualityHistory.reset();
	resetRetainedAdmissionQualityProof();
	/* A source, renderer, or user-quality reset changes the population or
	 * capacity contract behind a previously measured presentation cut. */
	resetDeadlineSafePresentation();
	lodViewQualityDomainRevision++;
	if (!lodViewQualityDomainRevision)
	    lodViewQualityDomainRevision = 1;
    }

    /* These values are annotations owned by one bounded retained-admission
     * pass.  Reset them as one transaction; clearing only a convenient subset
     * has previously carried stale cut/budget evidence into a successor pass
     * and reopened an otherwise terminal handoff. */
    void resetRetainedPassAnnotations(void)
    {
	lodRetainedPass.reset();
    }

    void retireRetainedRefinementObservation(void)
    {
	lodRetainedPass.retireRefinement();
    }

    void resetRendererPerformanceEvidence(void)
    {
	/* A renderer epoch owns every timing-derived scalar, not just the
	 * exact-camera LRU.  Retained PoP buffers and occurrence cuts are scene
	 * data and deliberately survive this reset; they provide a coherent,
	 * conservative first frame from which the replacement renderer can be
	 * measured. */
	resetLodViewQualityHistory();
	lodInteractiveCalibratedRenderCostPerSecond = 0.0L;
	lodStableCalibratedRenderCostPerSecond = 0.0L;
	lastRenderTimeNanoseconds = 0;
	smoothedRenderTimeNanoseconds = 0;
	lodRendererPerformanceEvidence.reset();
	lodInteractionStartCertificate.reset();
	resetRetainedPassAnnotations();

	/* Capacity ceilings, overload witnesses, and handoff proofs were all
	 * established by the previous renderer.  Start a fresh bounded budget
	 * epoch, but keep the currently applied renderer-only ceiling and point
	 * cut for the first frame.  Their explicit recovery latches will relax
	 * those controls using measurements from the new renderer without
	 * exposing a potentially enormous retained suffix speculatively. */
	this->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY);
	this->commitAdmissionPlan(BObolLodAdmissionPlanner::requestCapacityRescan(
	    this->lodAdmissionEvidence, this->lodAdmissionCursor));
	this->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
	lodPresentationPolicy.reset();
	this->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	this->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_STRUCTURAL_ADMISSION);
	lodPointQualityPhase.reset();
	if (lodPresentationPointProxyPixelThreshold > 1.01f)
	    lodPointQualityPhase.requestCalibration();
	lodStaticQualityTrial.restoreRendererCeiling(
	    lodInteractiveProgressiveCeiling >= 0);
	lodDiscretePopulationTrialPermit.revoke();
	lodInteractionSession.resetCeilingFeedback();
	lodPresentationTransaction.reset();
	lodRefinementNotBeforeMicroseconds = 0;

	/* Delivery-rate telemetry is also renderer-specific.  It is protected
	 * independently because hosts may sample it while an endpoint property
	 * change is being applied on the owner thread. */
	{
	    std::lock_guard<std::mutex> lock(presentationTimingMutex);
	    lastPresentationTimestampNanoseconds = 0;
	    smoothedPresentationIntervalNanoseconds = 0;
	    displayedPresentationIntervalNanoseconds = 0;
	    lastInteractivePresentationTimestampNanoseconds = 0;
	    smoothedInteractivePresentationIntervalNanoseconds = 0;
	}
    }

    /* Retire every automatic-LoD owner for a disabled view policy while
     * retaining the immutable geometry and renderer cuts which produced the
     * current framebuffer.  Those presentation values are a useful stable
     * fallback; resetCadPresentationLimits() belongs to scene/service
     * replacement, not to a policy toggle.
     *
     * Database progressive providers are intentionally absent from this
     * transaction.  They also realize ordinary no-LoD geometry and therefore
     * keep their independent providerPendingCount/wakeup contract. */
    void retireAutomaticLodControl(void)
    {
	if (lodService && lodActiveGeneration != 0)
	    lodService->cancelGeneration(lodActiveGeneration);
	lodActiveGeneration = 0;
	invalidateResidentMeshCompactionSnapshot();
	lodAvailabilityLedger.resetResultQueue();
	lodAvailabilityLedger.resetResidentGrowth();
	rewindLodSubmissionCursor();
	lodSubmissionPass.retire();
	lodSubmissionIntent.reset();
	lodSubmissionDelta.reset();
	lodLastSubmittedViewRevision.reset();
	lodLastSubmittedPolicyRevision.reset();
	lodLastSubmittedSources.clear();
	lodPlanningObligations.reset();
	lodViewDemandPolicy.reset();
	lodDiscretePopulationTrialPermit.revoke();
	lodInteractionSession.reset();
	lodInteractionStartCertificate.reset();
	lodPoseContinuity.reset();
	lodPresentationPolicy.reset();
	lodPresentationTransaction.reset();
	lodInterruptedPresentationReplay.retire();
	lodPointAdmissionFrame.retire();
	lodPointQualityPhase.reset();
	lodStaticQualityTrial.reset();
	lodStructuralRepair.reset();
	resetRetainedPassAnnotations();
	lodCompactionPolicy.retire();
	lodRefinementNotBeforeMicroseconds = 0;
	lodCoveragePolicy.reset();
	lodCoveragePolicy.setRequired(false);
	clearLodConvergenceCandidates();
	resetLodConvergenceFraction();
	applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY);
	applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
	applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_STRUCTURAL_ADMISSION);
	applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::MARK_RESOURCE_RECOVERY_HANDLED);
    }

    float quietAllocationTargetFps(void) const
    {
	/* The ordinary quiet cadence keeps cold streaming and multi-frame
	 * convergence responsive.  Once that phase has produced a terminal exact
	 * frame, an event-driven view may spend up to the separate hard static
	 * presentation deadline on a richer framebuffer which is then retained
	 * without redraw. */
	const uint64_t deadline = staticQualityPresentationDeadline();
	const BObolLodCapacitySearchCertificate &capacitySearch =
	    lodAdmissionEvidence.capacity().capacitySearch();
	const bool capacitySearchUsesStaticDeadline =
	    capacitySearch.phase() !=
		BObolLodCapacitySearchCertificate::Phase::INACTIVE &&
	    capacitySearch.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC;
	if ((!lodStaticQualityTrial.usesStaticDeadline() &&
	     !lodPointQualityPhase.staticCalibrationPending() &&
	     !capacitySearchUsesStaticDeadline) || !deadline)
	    return lodStableTargetFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(deadline);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return lodStableTargetFps;
	const float bounded = static_cast<float>(deadlineFps);
	return lodStableTargetFps > 0.0f ?
	    std::min(lodStableTargetFps, bounded) : bounded;
    }

    float staticQualityTargetFps(void) const
    {
	/* Once the ordinary one-pixel frame is exact and terminal, a static
	 * event-driven view may use the separately configured hard frame deadline
	 * to test a finer pixel target.  The framebuffer is retained afterward,
	 * so requiring every such terminal quality step to meet the ordinary 20 Hz
	 * refinement cadence makes warm-cache convergence depend on one upload or
	 * command-preparation sample. */
	const uint64_t deadline = staticQualityPresentationDeadline();
	if (!deadline)
	    return lodStableTargetFps;
	const long double deadlineFps = 1000000000.0L /
	    static_cast<long double>(deadline);
	if (!std::isfinite(deadlineFps) || deadlineFps <= 0.0L)
	    return lodStableTargetFps;
	const float bounded = static_cast<float>(deadlineFps);
	return lodStableTargetFps > 0.0f ?
	    std::min(lodStableTargetFps, bounded) : bounded;
    }

    uint64_t lastRenderTimeNanoseconds = 0;
    uint64_t smoothedRenderTimeNanoseconds = 0;
    uint64_t lastBackgroundRenderTimeNanoseconds = 0;
    uint64_t lastSceneRenderTimeNanoseconds = 0;
    /* Full CAD plan/atlas construction is a real latency problem, but not a
     * measurement of steady triangle throughput.  Only a prepared replay or
     * an unchanged retained upload population may drive quiet capacity cuts. */
    BObolLodRendererPerformanceEvidence lodRendererPerformanceEvidence;
    uint64_t lastProgressiveAdvanceTimeNanoseconds = 0;
    uint64_t lastLodResultProcessingTimeNanoseconds = 0;
    uint64_t lastProgressiveProviderTimeNanoseconds = 0;
    uint64_t lastLodSubmissionTimeNanoseconds = 0;
    uint64_t lastPresentationSyncTimeNanoseconds = 0;
    uint64_t interactivePresentationFrameDeadlineNanoseconds =
	defaultInteractivePresentationDeadline();
    uint64_t stablePresentationFrameDeadlineNanoseconds =
	defaultStablePresentationDeadline();
    uint64_t prominentQualityFrameDeadlineNanoseconds =
	defaultProminentQualityDeadline();
    uint64_t interruptedPresentationFrameCount = 0;
    uint64_t consecutiveInterruptedPresentationFrames = 0;
    uint64_t lastInterruptedPresentationTimeNanoseconds = 0;
    /* An interrupted CAD traversal owns an exact successor presentation.
     * Keep owner-thread scene publication frozen until that replay completes;
     * otherwise provider/result pumps can invalidate resumable command-plan
     * work between deadline slices.  Workers remain free to fill their
     * bounded queues while the presentation transaction is closed. */
    BObolLodInterruptedPresentationReplay lodInterruptedPresentationReplay;
    uint64_t renderCompletionSerial = 0;
    mutable std::mutex presentationTimingMutex;
    uint64_t presentedFrameSerial = 0;
    uint64_t lastPresentationTimestampNanoseconds = 0;
    uint64_t smoothedPresentationIntervalNanoseconds = 0;
    uint64_t displayedPresentationIntervalNanoseconds = 0;
    uint64_t lastInteractivePresentationTimestampNanoseconds = 0;
    uint64_t smoothedInteractivePresentationIntervalNanoseconds = 0;
    std::vector<BObolProgressiveProviderRecord> progressiveProviders;
    uint64_t progressiveProviderNextToken = 1;
    /* Number of registered non-LoD geometry producers which reported work
     * remaining in the most recent complete provider pass.  Registration is
     * conservatively pending until the first pass proves otherwise. */
    BObolLodAvailabilityLedger lodAvailabilityLedger;
    BObolProgressiveOptions defaultProgressiveOptions;
    BObolLodService *lodService = NULL;
    std::shared_ptr<BObolLodService> managedLodService;
    size_t managedLodWorkerCount = 0;
    uint64_t lodResultSubscriberId = 0;
    SbBool lodAutoSubmit = FALSE;
    uint64_t lodActiveGeneration = 0;
    size_t lodSubmissionSourceIndex = 0;
    size_t lodSubmissionEntryOffset = 0;
    BObolLodSubmissionSourcePlan lodSubmissionSourcePlan;
    /* Expensive stable-view minimax arithmetic is resumable.  The plan owns
     * no scene geometry and publishes nothing until its final epoch-checked
     * owner-thread commit. */
    std::shared_ptr<BObolRetainedAllocationTransaction>
	lodRetainedAllocationTransaction;
    BObolLodSubmissionIntent lodSubmissionIntent;
    BObolLodRetainedQualityEvidence lodRetainedAdmissionQualityEvidence;
    BObolRetainedAllocationResult lodRetainedAllocationCertificate;
    BObolLodSubmissionVisibleCount lodSubmissionVisibleCount;
    /* Large compact scenes are consumed in bounded GUI-thread windows.  A
     * coverage pass proves that every projected leaf has a useful structural
     * presentation before any of those leaves opens a PoP hierarchy.  The
     * following quality pass promotes view-significant leaves under the
     * aggregate budget, preventing early entries from starving later ones. */
    BObolLodCoveragePolicy lodCoveragePolicy;
    /*
     * A camera-scale refresh is not cold coverage.  Existing occurrence cuts
     * remain useful and a zoom must be allowed to retarget them, consume a
     * richer resident prefix, or request the missing cumulative suffix while
     * the bounded visibility scan continues.  New/source inventory coverage
     * keeps the stricter minimum-mesh-first rule below.
     */
    /*
     * The newest completed camera pass proved that every projected mesh
     * occurrence already has a retained mesh presentation.  While input is
     * active, richer worker results may then remain queued without exposing
     * boxes or missing newly visible geometry.
     */
    /* One byte per compact entry is deliberately used instead of a hash set:
     * entry indices are dense and source-stable, so a 150k-occurrence source
     * costs only 150 kB and an exact edit delta updates it in constant time. */
    BObolLodVisibilityCensus lodConvergenceCandidateCensus;
    /* Pose-only motion may reuse the retained assembly without projecting the
     * complete source population.  Keep that missing-proof obligation across
     * progressive pump calls; a function-local flag lets a later empty pump
     * publish the preceding pose's visibility count as if it described the
     * new camera. */
    BObolLodPoseContinuity lodPoseContinuity;
    size_t lodSourceLogicalOccurrenceCount = 0;
    /* Camera/geometric projection evidence is denser than the visibility
     * census but still bounded (roughly a few dozen bytes per occurrence).
     * It belongs to this view, never to the shared database source. */
    BObolLodProjectedDemandCache lodProjectedDemandCache;
    /*
     * Authoritative visible-mesh denominator from the newest completed
     * all-source coverage pass.  Per-source candidate counts are useful while
     * a bounded pass is still being assembled, but its final window used to
     * clear the accumulated coverage counters before convergence retained
     * their total.  A quality-only policy revision preserves this view proof;
     * camera and source-inventory revisions invalidate it.
    */
    mutable BObolLodConvergencePolicy lodConvergencePolicy;
    void clearLodSubmissionSourceState(void)
    {
	lodSubmissionSourcePlan.reset();
	lodRetainedAllocationTransaction.reset();
	lodSubmissionVisibleCount.reset();
    }
    void positionLodSubmissionCursor(size_t sourceIndex)
    {
	lodSubmissionSourceIndex = sourceIndex;
	lodSubmissionEntryOffset = 0;
	clearLodSubmissionSourceState();
    }
    void rewindLodSubmissionCursor(void)
    {
	positionLodSubmissionCursor(0);
    }
    void clearLodConvergenceCandidates(void)
    {
	lodConvergenceCandidateCensus.clear();
	lodPoseContinuity.completeVisibilityCensus();
	lodCoveragePolicy.clearCompleteVisibleCount();
    }
    void resetLodConvergenceFraction(void)
    {
	lodConvergencePolicy.resetFraction();
    }
    static BObolLodVisibilityCensus::SourceKey lodConvergenceSourceKey(
	SoBRLDatabaseSource *source)
    {
	return source ? source->getCompactSourceRoutingId() : 0;
    }
    size_t lodConvergenceCandidateCensusTotal(void) const
    {
	return lodConvergenceCandidateCensus.total();
    }
    void publishExactLodConvergenceCandidateCount(void)
    {
	/* Preserve the coverage proof itself.  Exact deltas mutate only entries
	 * represented by the census; population-changing deltas clear it and run
	 * a complete pass instead. */
	if (lodCoveragePolicy.hasCompleteVisibleCount())
	    lodCoveragePolicy.setCompleteVisibleCount(
		lodConvergenceCandidateCensusTotal());
    }
    void setLodConvergenceCandidateCount(
	SoBRLDatabaseSource *source, size_t count)
    {
	if (!source)
	    return;
	lodConvergenceCandidateCensus.setCount(
	    lodConvergenceSourceKey(source), count);
	publishExactLodConvergenceCandidateCount();
    }
    void beginLodConvergenceCandidateCensus(
	SoBRLDatabaseSource *source, size_t entryCount)
    {
	if (!source)
	    return;
	lodConvergenceCandidateCensus.begin(
	    lodConvergenceSourceKey(source), entryCount);
    }
    void observeLodConvergenceCandidateVisibility(
	SoBRLDatabaseSource *source, size_t entryCount,
	const std::vector<std::pair<size_t, SbBool>> &observations)
    {
	if (!source)
	    return;
	for (const std::pair<size_t, SbBool> &observation : observations) {
	    lodConvergenceCandidateCensus.observe(
		lodConvergenceSourceKey(source), entryCount,
		observation.first, observation.second != FALSE);
	}
	if (lodConvergenceCandidateCensus.complete(
		lodConvergenceSourceKey(source)))
	    publishExactLodConvergenceCandidateCount();
    }
    void completeLodConvergenceCandidateCensus(
	SoBRLDatabaseSource *source)
    {
	if (!source)
	    return;
	/* The dense entry bits are the authority.  A provider result can restart a
	 * bounded source cursor while retaining the coverage transaction; an
	 * action-level accumulated scalar can therefore count a replayed prefix
	 * twice.  Entry observations are idempotent and remain exact across that
	 * restart, as required by later O(delta) visibility edits. */
	lodConvergenceCandidateCensus.finish(
	    lodConvergenceSourceKey(source));
	publishExactLodConvergenceCandidateCount();
    }
    size_t lodConvergenceCandidateCount(void) const
    {
	if (lodCoveragePolicy.hasCompleteVisibleCount())
	    return lodCoveragePolicy.completeVisibleCount();
	return lodConvergenceCandidateCensusTotal();
    }
    bool hasCompleteLodConvergenceCandidateCensus(
	SoBRLDatabaseSource *source) const
    {
	return source && lodConvergenceCandidateCensus.complete(
	    lodConvergenceSourceKey(source));
    }
    bool isVisibleLodConvergenceCandidate(
	SoBRLDatabaseSource *source, size_t entryIndex) const
    {
	return source && lodConvergenceCandidateCensus.visible(
	    lodConvergenceSourceKey(source), entryIndex);
    }
    bool pointProxyAggregationApplicable(void) const
    {
	/* A large mesh may be internally clustered or streamed in multiple
	 * batches, but it remains one logical CAD occurrence.  Point aggregation
	 * can only reduce the draw-call floor of multiple occurrences; applying it
	 * to a single Lucy-like mesh destroys the silhouette it is meant to
	 * preserve.  A source-wide one-occurrence proof therefore vetoes a
	 * transient discovery/cluster census.  Multi-occurrence sources still use
	 * the camera-local visible census, so a mostly off-screen large assembly
	 * does not collapse its one visible feature. */
	if (lodSourceLogicalOccurrenceCount == 1)
	    return false;
	return BObolLodAdmissionPlanner::pointAggregationApplicable(
	    lodConvergenceCandidateCount());
    }
    bool retainedPointProxyAggregationApplicable(void) const
    {
	/* Stable/deadline pressure may only tune a representation which the
	 * current occurrence allocator can actually produce.  Logical occurrence
	 * count alone is sufficient during structural discovery, but it falsely
	 * classified a handful of large, internally chunked meshes as a subpixel
	 * population and repeatedly raised an inert threshold. */
	return BObolLodAdmissionPlanner::retainedPointAggregationApplicable(
	    pointProxyAggregationApplicable(),
	    lodRetainedAllocationCertificate.pointProxyCandidateCount);
    }
    bool deadlinePointProxyAggregationApplicable(void) const
    {
	return BObolLodAdmissionPlanner::deadlinePointAggregationApplicable(
	    retainedPointProxyAggregationApplicable(),
	    lodPointAdmissionFrame.pending(), lodSourceLogicalOccurrenceCount);
    }
    bool pointProxyAggregationApplicableForCameraTransition(void) const
    {
	/* Interaction entry precedes the first changed camera signature.  The
	 * replacement visibility census may therefore be empty even though the
	 * current framebuffer already contains a proven aggregate-point
	 * population.  Preserve that presentation as the transition baseline;
	 * exact frames in the new view will replace its capacity evidence. */
	return BObolLodAdmissionPlanner::
	    pointAggregationApplicableAcrossCameraInvalidation(
		pointProxyAggregationApplicable(),
		lodPresentationPointProxyPixelThreshold);
    }
    void observeLodSourceLogicalOccurrenceCount(size_t count)
    {
	if (count > 0)
	    lodSourceLogicalOccurrenceCount = count;
    }
    bool structuralPointAggregationRequired(
	size_t affordableOccurrenceCount) const
    {
	return BObolLodAdmissionPlanner::structuralPointAggregationRequired(
	    lodSourceLogicalOccurrenceCount, affordableOccurrenceCount);
    }
    bool structuralPointAggregationRequiredForBudget(size_t sceneBudget) const
    {
	return structuralPointAggregationRequired(
	    BObolLodAdmissionPlanner::structuralFirstWaveOccurrenceLimit(
		sceneBudget));
    }
    BObolLodSubmissionPass lodSubmissionPass;
    /* A retained minimax pass selected a level whose PoP suffix is not yet
     * resident.  This is provider work, not a performance-limited quality
     * observation, and must survive the coherent cut presentation barrier. */
    /* The last exact renderer frame observed structural boxes which the
     * predictive point classifier expected to collapse.  One bounded pass
     * must bypass that prediction and obtain their mesh presentations. */
    /* Exact non-terminal box population owned by the current closed repair
     * transaction.  The per-occurrence share is derived only after the
     * admission plan freezes the active-cost baseline which will consume it. */
    BObolLodStructuralRepair lodStructuralRepair;
    /* Accumulated across every bounded window of one scene pass.  Unlike the
     * general refinement-debt bit, this counts only visible box-to-first-mesh
     * admissions rejected by the finite scene allowance. */
    BObolLodRetainedPassAnnotations lodRetainedPass;
    SbBool forceTerminalLodRefinement = FALSE;
    BObolLodViewEpoch lodLastSubmittedViewRevision;
    BObolLodPolicyEpoch lodLastSubmittedPolicyRevision;
    struct LodSourceSnapshot {
	SoBRLDatabaseSource *source = NULL;
	struct db_i *database = NULL;
	BObolLodSourceRoutingId routingId;
	BObolLodInventoryEpoch inventoryRevision;
	SbString databaseId;
	SbString path;
	int drawMode = 0;
	int representationMode = 0;
	SbBool visible = FALSE;
	int lodBotThreshold = 0;
	uint32_t sourceRevision = 0;
	uint32_t inputsRevision = 0;

	bool sameIdentity(const LodSourceSnapshot &other) const
	{
	    return this->database == other.database &&
		this->routingId == other.routingId &&
		this->databaseId == other.databaseId &&
		this->path == other.path &&
		this->drawMode == other.drawMode &&
		this->representationMode == other.representationMode &&
		this->visible == other.visible &&
		this->lodBotThreshold == other.lodBotThreshold &&
		this->sourceRevision == other.sourceRevision &&
		this->inputsRevision == other.inputsRevision;
	}
    };
    std::vector<LodSourceSnapshot> lodLastSubmittedSources;
    BObolLodSubmissionDelta lodSubmissionDelta;
    std::optional<BObolLodViewSnapshot> lodViewSignature;
    BObolLodViewScaleSnapshot lodViewScaleSignature;

    /* Net scale, not the mere presence of a wheel event, decides whether a
     * quiet orthographic view must retarget every occurrence.  Keep the
     * signature and ready-view proof from the start of the continuous input
     * epoch.  A bracketed pose gesture may begin before an unbracketed wheel
     * epoch's debounce expires; in that case these values deliberately keep
     * describing the original stable view. */
    BObolLodInteractionStartCertificate lodInteractionStartCertificate;
    /* A completed, fully covered orthographic pose may verify visibility,
     * selection, and resource pressure, but may not rewrite retained PoP
     * cuts merely because the interaction debounce ended. */
    /* A settled retained camera epoch first records exact projected demand
     * for every visible occurrence, then performs one scene-budgeted
     * importance reallocation.  Keeping this as a census-completion latch
     * prevents stale previous-view metrics and makes the redistribution
     * explicitly one-shot. */
    BObolLodPlanningObligations lodPlanningObligations;
    /* A resident-capacity revision reopens only occurrences whose provider
     * was denied by an older admission epoch.  Keep this distinct from the
     * ordinary unsatisfied-quality frontier: reclaimed bytes must not restart
     * a 150k-entry view census for a handful of denied assets. */
    BObolLodViewEpoch lodViewRevision {1};
    BObolLodPolicyEpoch lodPolicyRevision {1};
    BObolLodViewDemandPolicy lodViewDemandPolicy;
    /* One scale-quality or exact stable-headroom probe may measure one
     * otherwise-unaffordable populated PoP transition across the entire
     * scene.  Individual submit actions are time-sliced and source-local, so
     * uniqueness must live here. */
    BObolLodDiscreteTrialPermit lodDiscretePopulationTrialPermit;
    BObolLodInteractionSession lodInteractionSession;
    /* Renderer-only presentation continuity and motion-to-stable handoff.
     * The policy owns its snapshot/latch proof; this coordinator only applies
     * returned ceiling and point-proxy values to the retained view state. */
    BObolLodPresentationPolicy lodPresentationPolicy;
    /* Exact-view history survives camera epochs, but never source inventory,
     * service, or user quality-contract changes.  The separate domain makes
     * that lifetime explicit and prevents a coincidentally identical camera
     * from accepting evidence belonging to a different scene. */
    BObolLodViewQualityHistory lodViewQualityHistory;
    uint64_t lodViewQualityDomainRevision = 1;
    BObolLodCompactionPolicy lodCompactionPolicy;
    /* Complete-frame resource pressure is retained independently of HUD
     * queries.  A pressure edge requests one safe compaction pass; if the
     * current visible working set itself exceeds a renderer limit, that pass
     * may finish in a stable memory-limited presentation rather than loop. */
    BObolLodAdmissionEvidence lodAdmissionEvidence;
    BObolLodAdmissionCursor lodAdmissionCursor;
    BObolLodInventoryEpoch lodAdmissionInventoryRevision {1};
    BObolLodAvailabilityEpoch lodAdmissionAvailabilityRevision {1};
    BObolLodCapacityEpoch lodAdmissionCapacityRevision {1};
    uint64_t lodResidentAdmissionRevision = 0;
    /* Applied results and their frame acknowledgement are one revision-bound
     * transaction.  A publication frame may include one-time CPU/GPU buffer
     * construction, so it controls refinement pacing without teaching the
     * steady-state capacity estimator that later frames cost the same. */
    BObolLodPresentationTransaction lodPresentationTransaction;
    int64_t lodRefinementNotBeforeMicroseconds = 0;
    float lodTargetPixelError = 1.0f;
    float lodInteractiveTargetFps = 60.0f;
    float lodStableTargetFps = 20.0f;
    /* Event-driven static quality and its terminal capacity witness share one
     * explicit lifecycle.  Internal policy bookkeeping may deactivate a
     * trial while preserving its constraint; genuine camera, service, user
     * policy, and cadence epochs reset the complete value. */
    BObolLodStaticQualityTrial lodStaticQualityTrial;
    long double lodInteractiveCalibratedRenderCostPerSecond = 0.0L;
    long double lodStableCalibratedRenderCostPerSecond = 0.0L;
    void seedInteractiveCalibrationFromStable(void)
    {
	if (!(lodStableCalibratedRenderCostPerSecond > 0.0L))
	    return;
	/* Stable retained rendering is a measured geometry-throughput witness,
	 * but motion also pays camera/cut/update overhead.  Carry half of that
	 * known capacity into a new gesture.  This is deliberately a floor, not
	 * a reset: an established faster interaction estimate is retained, while
	 * a 118-triangle underloaded frame can no longer become a permanent
	 * self-fulfilling capacity ceiling. */
	const long double seed =
	    lodStableCalibratedRenderCostPerSecond * 0.5L;
	lodInteractiveCalibratedRenderCostPerSecond = std::max(
	    lodInteractiveCalibratedRenderCostPerSecond, seed);
    }
    /* A late reusable frame may demonstrate headroom after the ordinary
     * bounded probe series has ended.  Remember the exact camera/policy and
     * budget already retried so each newly enlarged allowance gets one
     * admission pass without permitting an unchanged discrete-cut loop. */
    int lodInteractiveProgressiveCeiling = -1;
    /* Richest renderer-only ceiling which completed inside the hard deadline
     * for the current camera and source-quality domain.  An interaction and
     * its subsequent quiet phase differ only in policy, not in the measured
     * retained data or endpoint capacity.  If a later quiet probe misses,
     * return directly to this proven cut instead of replaying every
     * intervening PoP ordinal. */
    BObolLodDeadlineSafePresentation lodDeadlineSafePresentation;
    /*
     * Current renderer-side small-occurrence aggregation cut.  Interaction
     * raises it immediately when measured frames miss their target.  A quiet
     * performance-limited view may also raise it when every retained mesh is
     * already at its minimum prefix and that irreducible population still
     * exceeds the calibrated scene budget.
     */
    float lodPresentationPointProxyPixelThreshold = 1.0f;
    /* Streaming discovery may temporarily aggregate tiny structural leaf
     * proxies more aggressively than the measured stable-view policy.  It is
     * renderer-only and is retired as soon as the producer inventory settles;
     * keeping a separate value prevents discovery pacing from contaminating
     * the terminal point-calibration bracket. */
    float lodDiscoveryPointProxyPixelThreshold = 1.0f;
    /* A discovery-release or structural-distribution threshold mutation must
     * be presented exactly before mesh admission may consume its classifier.
     * This is intentionally distinct from stable point calibration, which
     * measures an already realized CAD population and never blocks source
     * work. */
    BObolLodPointAdmissionFrame lodPointAdmissionFrame;
    /*
     * A changed quiet aggregation threshold needs one measured presentation
     * before convergence is authoritative.  Keep this distinct from PoP
     * admission: no provider scan or retained-array mutation is required.
     */
    BObolLodPointQualityPhase lodPointQualityPhase;
    /* A completed retained-recovery pass owns one authoritative triangle
     * allocation plan.  Point-threshold calibration may change only the
     * renderer-side small-occurrence classification afterward; it must not
     * submit the same no-op triangle recovery again.  The plan serial plus
     * cadAllocationPlanCoversCurrentPopulation() make this witness expire
     * automatically on a view, policy, occurrence, or resident-mesh change. */
    std::optional<int> lodForcedCut;
    uint64_t maxExactFullDetailFaceCount = 0;
    uint64_t maxExactFullDetailPointCount = 0;
    std::vector<BObolRtPickCache *> rtPickCaches;
    std::vector<SbString> rtPickCachePaths;
    std::vector<struct db_i *> rtPickCacheDatabases;
    std::vector<uint32_t> rtPickCacheSourceRevisions;
    SbBool meshResidencyBudgetEnabled = FALSE;
    size_t maxResidentMeshBytes = 0;
    SbBool meshResidencyEvictDisplayPayloads = TRUE;
    size_t lastMeshBudgetInitialResidentBytes = 0;
    size_t lastMeshBudgetFinalResidentBytes = 0;
    size_t lastMeshBudgetFreedFullDetailBytes = 0;
    size_t lastMeshBudgetFreedDisplayBytes = 0;
    unsigned int lastMeshBudgetVisitedMeshCount = 0;
    unsigned int lastMeshBudgetEvictedFullDetailMeshCount = 0;
    unsigned int lastMeshBudgetEvictedDisplayMeshCount = 0;
    unsigned int lastLodVisitedMeshCount = 0;
    unsigned int lastLodSubmittedTaskCount = 0;
    unsigned int lastLodUpdatedCutCount = 0;
    unsigned int lastLodSkippedMeshCount = 0;
    size_t lastLodResultCount = 0;
    unsigned int lastLodMatchedResultCount = 0;
    unsigned int lastLodAppliedResultCount = 0;
    unsigned int lastLodRejectedResultCount = 0;
    unsigned int lastLodUnmatchedResultCount = 0;
    SbString lastLodDiagnostics = SbString("");
};


#endif /* LIBBOBOL_VIEW_LOD_COORDINATOR_STATE_PRIVATE_H */
