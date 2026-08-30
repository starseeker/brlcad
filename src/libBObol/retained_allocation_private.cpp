/*       R E T A I N E D _ A L L O C A T I O N _ P R I V A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "retained_allocation_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodRealization.h"
#include "BObol/BViewLod.h"
#include "cad_assembly_private.h"
#include "lod_visual_importance_private.h"

#include <Obol/cad/CadProjectedProxy.h>

#include "bu/datetime.h"
#include "bu/log.h"

#include <algorithm>
#include <cmath>
#include <cstring>
#include <cstdint>
#include <limits>
#include <memory>
#include <new>
#include <queue>
#include <string>
#include <unordered_map>
#include <vector>

bool
BObolRetainedAllocationInputKey::operator==(
    const BObolRetainedAllocationInputKey &other) const
{
    return this->externalPresentationCost == other.externalPresentationCost &&
	this->sceneBudget == other.sceneBudget &&
	this->maximumMarginalBudget == other.maximumMarginalBudget &&
	this->maximumProtectedBudget == other.maximumProtectedBudget &&
	this->revisionStamp.same(other.revisionStamp) &&
	this->residentAdmissionRevision == other.residentAdmissionRevision &&
	std::memcmp(&this->pointProxyPixelThreshold,
	    &other.pointProxyPixelThreshold, sizeof(float)) == 0 &&
	this->allowProtectedFloor == other.allowProtectedFloor;
}

BObolRetainedAllocationInputKey
BObolRetainedAllocationInputs::inputKey(void) const
{
    BObolRetainedAllocationInputKey key;
    key.externalPresentationCost = this->externalPresentationCost;
    key.sceneBudget = this->sceneBudget;
    key.maximumMarginalBudget = this->maximumMarginalBudget;
    key.maximumProtectedBudget = this->effectiveMaximumProtectedBudget();
    key.revisionStamp = this->revisionStamp;
    key.residentAdmissionRevision = this->residentAdmissionRevision;
    key.pointProxyPixelThreshold = this->pointProxyPixelThreshold;
    key.allowProtectedFloor = this->allowProtectedFloor;
    return key;
}

BObolRetainedAllocationInputKey
BObolRetainedAllocationResult::inputKey(void) const
{
    BObolRetainedAllocationInputKey key;
    key.externalPresentationCost = this->externalPresentationCost;
    key.sceneBudget = this->requestedSceneBudget;
    key.maximumMarginalBudget = this->maximumMarginalBudget;
    key.maximumProtectedBudget = this->allowProtectedFloor ?
	this->maximumProtectedBudget : 0;
    key.revisionStamp = this->revisionStamp;
    key.residentAdmissionRevision = this->residentAdmissionRevision;
    key.pointProxyPixelThreshold = this->pointProxyPixelThreshold;
    key.allowProtectedFloor = this->allowProtectedFloor;
    return key;
}

bool
bobol_retained_marginal_lower_priority(
    const BObolRetainedMarginalUpgrade &a,
    const BObolRetainedMarginalUpgrade &b)
{
    if (a.visualEmphasis != b.visualEmphasis)
	return a.visualEmphasis < b.visualEmphasis;
    if (a.qualityFloorViolation != b.qualityFloorViolation)
	return !a.qualityFloorViolation;
    if (a.qualityFloorViolation &&
	(a.normalizedError < b.normalizedError ||
	 a.normalizedError > b.normalizedError))
	return a.normalizedError < b.normalizedError;
    if (a.visualBenefitPerCost < b.visualBenefitPerCost ||
	a.visualBenefitPerCost > b.visualBenefitPerCost)
	return a.visualBenefitPerCost < b.visualBenefitPerCost;
    if (a.normalizedError < b.normalizedError ||
	a.normalizedError > b.normalizedError)
	return a.normalizedError < b.normalizedError;
    return a.candidateIndex > b.candidateIndex;
}

class BObolRetainedAllocationTransaction {
public:
    explicit BObolRetainedAllocationTransaction(
	const BObolRetainedAllocationInputs &inputs)
    {
	this->inputKey = inputs.inputKey();
	this->viewState = inputs.viewState;
	this->revisionStamp = inputs.revisionStamp;
	this->residentAdmissionRevision = inputs.residentAdmissionRevision;
	this->cadRevision = inputs.viewState ?
	    inputs.viewState->cadRevision() : 0;
	this->residentDemandRevision = inputs.viewState ?
	    inputs.viewState->residentMeshDemandRevision() : 0;
	this->allocationPlanSerial = inputs.viewState ?
	    inputs.viewState->beginCadAllocationPlan() : 0;
	this->externalPresentationCost = inputs.externalPresentationCost;
	this->fixedCost = inputs.externalPresentationCost;
	this->pixelDemandPresentationCost = inputs.externalPresentationCost;
	this->sceneBudget = inputs.sceneBudget;
	this->maximumMarginalBudget = inputs.maximumMarginalBudget;
	this->allowProtectedFloor = inputs.allowProtectedFloor;
	this->maximumProtectedBudget =
	    inputs.effectiveMaximumProtectedBudget();
	this->pointProxyPixelThreshold = inputs.pointProxyPixelThreshold;
	this->wallStartedMicroseconds = bu_gettime();
	if (!inputs.sources)
	    return;
	this->sources.reserve(inputs.sources->size());
	size_t payloadCount = 0;
	for (SoBRLDatabaseSource *source : *inputs.sources) {
	    SourceSnapshot snapshot;
	    snapshot.source = source;
	    if (source) {
		snapshot.drawMode = source->getEffectiveLodDrawMode();
		snapshot.assembly = this->viewState ?
		    this->viewState->findCadPresentation(source) : NULL;
		if (snapshot.assembly &&
		    this->pointProxyProtectedInstances.find(snapshot.assembly) ==
			this->pointProxyProtectedInstances.end()) {
		    ProtectedSet protection;
		    protection.current = snapshot.assembly->
			pointProxyProtectedInstances();
		    protection.next.reserve(protection.current.size());
		    this->pointProxyProtectedInstances.emplace(
			snapshot.assembly, std::move(protection));
		    this->pointProxyAssemblies.push_back(snapshot.assembly);
		}
		if (this->viewState)
		    /* The marginal allocator uses candidate position as its final
		     * equal-value tie breaker.  A retained allocation may be retried
		     * after a measured frame, so unordered-map iteration here would
		     * select a different but budget-equivalent cut assignment on each
		     * pass.  That churn is a visible mutation and can keep a static
		     * large scene perpetually in its handoff cycle.  Allocation itself
		     * is an intentionally infrequent O(scene) transaction; canonical
		     * occurrence ordering is therefore the correct trade for a stable,
		     * reproducible plan. */
		    this->viewState->findCadPayloads(source, snapshot.payloads);
		payloadCount = snapshot.payloads.size() >
			SIZE_MAX - payloadCount ? SIZE_MAX :
		    payloadCount + snapshot.payloads.size();
	    }
	    this->sources.push_back(std::move(snapshot));
	}
	this->candidates.reserve(payloadCount);
    }

    bool matches(const BObolRetainedAllocationInputs &inputs) const
    {
	if (!inputs.sources || !inputs.viewState ||
	    this->viewState != inputs.viewState ||
	    !this->revisionStamp.same(inputs.revisionStamp) ||
	    this->cadRevision != inputs.viewState->cadRevision() ||
	    this->residentDemandRevision !=
		inputs.viewState->residentMeshDemandRevision() ||
	    !inputs.viewState->isCadAllocationPlanCurrent(
	    this->allocationPlanSerial) ||
	    this->inputKey != inputs.inputKey() ||
	    this->sources.size() != inputs.sources->size())
	    return false;
	for (size_t i = 0; i < this->sources.size(); ++i) {
	    SoBRLDatabaseSource *source = (*inputs.sources)[i];
	    if (this->sources[i].source != source)
		return false;
	    if (source && this->sources[i].drawMode !=
		    source->getEffectiveLodDrawMode())
		return false;
	}
	return true;
    }

    BObolRetainedAllocationStatus advance(
	const BObolRetainedAllocationInputs &inputs,
	uint64_t sliceMicroseconds)
    {
	if (!this->matches(inputs))
	    return BOBOL_RETAINED_ALLOCATION_STALE;
	const int64_t sliceStarted = bu_gettime();
	int64_t phaseStarted = sliceStarted;
	auto transition = [this, &phaseStarted](Phase next) {
	    const int64_t now = bu_gettime();
	    this->phaseMicroseconds[this->phase] +=
		std::max<int64_t>(0, now - phaseStarted);
	    this->phase = next;
	    phaseStarted = now;
	};
	size_t workSinceClockCheck = 0;
	auto shouldYield = [&]() {
	    if (!sliceMicroseconds)
		return false;
	    if (++workSinceClockCheck < 64)
		return false;
	    workSinceClockCheck = 0;
	    return bu_gettime() - sliceStarted >=
		static_cast<int64_t>(sliceMicroseconds);
	};
	auto finishSlice = [this, &phaseStarted]() {
	    const int64_t now = bu_gettime();
	    this->phaseMicroseconds[this->phase] +=
		std::max<int64_t>(0, now - phaseStarted);
	};
	try {
	    for (;;) {
		switch (this->phase) {
		    case DISCOVERY:
			while (this->sourceCursor < this->sources.size()) {
			    SourceSnapshot &source =
				this->sources[this->sourceCursor];
			    while (this->payloadCursor <
				    source.payloads.size()) {
				const BObolViewLodState::CadPayload *payload =
				    source.payloads[this->payloadCursor++];
				this->discover(source, payload);
				if (shouldYield()) {
				    finishSlice();
				    return BOBOL_RETAINED_ALLOCATION_PENDING;
				}
			    }
			    this->payloadCursor = 0;
			    this->sourceCursor++;
			}
			transition(FLOOR);
			break;
		    case PROTECTION_COMPARE:
			while (this->pointProxyAssemblyCursor <
				this->pointProxyAssemblies.size()) {
			    SoCADAssembly *assembly = this->pointProxyAssemblies[
				this->pointProxyAssemblyCursor];
			    ProtectedSet &protection =
				this->pointProxyProtectedInstances[assembly];
			    if (!protection.changed &&
				protection.current.size() != protection.next.size())
				protection.changed = true;
			    while (!protection.changed &&
				protection.compareCursor <
				    protection.current.size()) {
				if (!protection.next.count(protection.current[
					protection.compareCursor++]))
				    protection.changed = true;
				if (shouldYield()) {
				    finishSlice();
				    return BOBOL_RETAINED_ALLOCATION_PENDING;
				}
			    }
			    this->pointProxyAssemblyCursor++;
			}
			transition(COMMIT);
			break;
		    case FLOOR: {
			/* A point representation only reduces scene work when there is a
			 * population to aggregate.  For one terminal occurrence it would
			 * merely hide the sole mesh while retaining the same draw call, and
			 * it can force an unnecessary point-to-mesh republish when the
			 * calibration threshold changes. */
			this->useOptionalPointCandidates =
			    this->hasOptionalPointCandidates &&
			    this->candidates.size() > 1;
			if (this->useOptionalPointCandidates) {
			    this->protectedCandidateCost =
				this->protectedCandidateCost == SIZE_MAX ||
				this->optionalPointCandidateCost == SIZE_MAX ?
				    SIZE_MAX : this->protectedCandidateCost -
				    this->optionalPointCandidateCost;
			    this->protectedFloorSignature ^=
				this->optionalPointCandidateSignature;
			}
			this->protectedFloorCost =
			    this->protectedCandidateCost >
				SIZE_MAX - this->fixedCost ? SIZE_MAX :
			    this->fixedCost + this->protectedCandidateCost;
			this->protectedFloorSignature ^=
			    static_cast<uint64_t>(this->candidates.size()) *
				0x9e3779b97f4a7c15ULL;
			this->protectedFloorSignature ^=
			    static_cast<uint64_t>(this->protectedFloorCost) *
				0x517cc1b727220a95ULL;
			const size_t protectedBudget = std::max(
			    this->sceneBudget, this->maximumProtectedBudget);
			const size_t trialBudget =
			    this->allowProtectedFloor &&
			    this->maximumProtectedBudget > 0 ?
				(protectedBudget > SIZE_MAX -
				    protectedBudget / protectedFloorTrialDivisor ?
				    SIZE_MAX : protectedBudget +
				    protectedBudget / protectedFloorTrialDivisor) :
				protectedBudget;
			this->enforceProtectedFloor = this->allowProtectedFloor &&
			    this->protectedFloorCost <= trialBudget;
			this->effectiveBudget = this->enforceProtectedFloor ?
			    std::max(this->sceneBudget,
				this->protectedFloorCost) :
			    std::max(this->sceneBudget,
				this->maximumMarginalBudget);
			this->allocationBudget = this->enforceProtectedFloor ?
			    this->protectedFloorCost : 0;
			if (this->candidates.empty()) {
			    this->beginBaseline();
			    transition(BASELINE);
			    break;
			}
			/* Capacity search is a bracket over numeric budgets and therefore
			 * requires the allocator's selected populations to be nested.  A
			 * budget-dependent common-error pre-pass could reshuffle several
			 * cuts when the budget changed by one unit, even when the submitted
			 * render cost stayed unchanged.  Start every request from the same
			 * coherent point/minimum/protected floor and let the deterministic
			 * marginal queue extend that population.  The queue's value-per-cost
			 * and weighted-error ordering retains the visual minimax behavior
			 * without creating a second, non-monotonic allocator. */
			this->beginBaseline();
			transition(BASELINE);
			break;
		    }
		    case BASELINE:
			while (this->finalCursor < this->candidates.size()) {
			    const size_t i = this->finalCursor++;
			    const Candidate &candidate = this->candidates[i];
			    this->finalPointProxies[i] =
				this->useOptionalPointCandidates &&
				candidate.pointProxyEligible;
			    this->finalCuts[i] = this->finalPointProxies[i] ?
				candidate.minimumCut :
				(this->enforceProtectedFloor ?
				    candidate.protectedCut : candidate.minimumCut);
			    this->finalCosts[i] = this->finalPointProxies[i] ?
				candidate.pointProxyCost :
				bobol_lod_render_cost_units(
				    this->countsAtCut(candidate, this->finalCuts[i]),
				    candidate.drawMode, 1);
			    this->finalCost = this->finalCosts[i] >
				SIZE_MAX - this->finalCost ? SIZE_MAX :
				this->finalCost + this->finalCosts[i];
			    if (this->finalPointProxies[i])
				this->observePointProxyError(candidate);
			    else {
				this->observeError(candidate, this->finalCuts[i]);
				this->protectPointProxyCandidate(candidate);
			    }
			    /* A point allocation is already realized by the aggregate
			     * presentation.  Stage a mesh cut only when this plan actually
			     * selects a mesh; a later marginal point-to-mesh upgrade stages
			     * its cut in MARGINAL_APPLY. */
			    if (!this->finalPointProxies[i] &&
				!this->stageAllocation(candidate.payload,
				    this->finalCuts[i], candidate.drawMode)) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_FAILED;
			    }
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			this->finalCursor = 0;
			transition(MARGINAL_SEED);
			break;
		    case MARGINAL_SEED:
			if (this->finalCost >= this->effectiveBudget) {
			    transition(PROTECTION_COMPARE);
			    break;
			}
			while (this->finalCursor < this->candidates.size()) {
			    this->queueNextUpgrade(this->finalCursor++);
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			transition(MARGINAL_APPLY);
			break;
		    case MARGINAL_APPLY:
			while (!this->upgrades.empty() &&
			    this->finalCost < this->effectiveBudget) {
			    const BObolRetainedMarginalUpgrade upgrade =
				this->upgrades.top();
			    this->upgrades.pop();
			    if (upgrade.addedCost <=
				    this->effectiveBudget - this->finalCost) {
				this->finalCost += upgrade.addedCost;
				this->finalPointProxies[upgrade.candidateIndex] = FALSE;
				this->finalCuts[upgrade.candidateIndex] =
				    upgrade.nextCut;
				this->finalCosts[upgrade.candidateIndex] =
				    upgrade.nextCost;
				if (!this->stageAllocation(
					this->candidates[upgrade.candidateIndex].payload,
					upgrade.nextCut,
					this->candidates[upgrade.candidateIndex].drawMode)) {
				    finishSlice();
				    return BOBOL_RETAINED_ALLOCATION_FAILED;
				}
				this->protectPointProxyCandidate(
				    this->candidates[upgrade.candidateIndex]);
				this->marginalAccepted = true;
				this->queueNextUpgrade(upgrade.candidateIndex);
			    } else {
				/* The selected population is a prefix of one stable
				 * importance order.  Skipping an unaffordable high-ranked
				 * upgrade lets a lower budget spend its remainder on later
				 * upgrades which a higher budget may then displace.  Stop at
				 * the first gap so a richer budget can only extend this exact
				 * population. */
				while (!this->upgrades.empty())
				    this->upgrades.pop();
			    }
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			if (this->marginalAccepted) {
			    this->finalCursor = 0;
			    this->maximumNormalizedError = 1.0;
			    this->maximumPresentedPixelError =
				this->fixedMaximumPresentedPixelError;
			    this->havePresentedErrorProof =
				this->fixedPresentedErrorProof;
			    transition(PROOF);
			} else {
			    transition(PROTECTION_COMPARE);
			}
			break;
		    case PROOF:
			while (this->finalCursor < this->candidates.size()) {
			    const size_t i = this->finalCursor++;
			    if (this->finalPointProxies[i])
				this->observePointProxyError(this->candidates[i]);
			    else
				this->observeError(
				    this->candidates[i], this->finalCuts[i]);
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			transition(PROTECTION_COMPARE);
			break;
		    case COMMIT:
			if (!this->matches(inputs) ||
			    !this->viewState->commitCadAllocationPlan(
				this->allocationPlanSerial, this->cadRevision,
				this->residentDemandRevision,
				this->viewRevision(), this->policyRevision(),
				this->fixedCadPresentationCost)) {
			    finishSlice();
			    return BOBOL_RETAINED_ALLOCATION_STALE;
			}
			this->commit();
			transition(COMPLETE);
			break;
		    case COMPLETE:
			finishSlice();
			return BOBOL_RETAINED_ALLOCATION_COMPLETE;
		    case PHASE_COUNT:
			finishSlice();
			return BOBOL_RETAINED_ALLOCATION_FAILED;
		}
	    }
	} catch (const std::bad_alloc &) {
	    finishSlice();
	    return BOBOL_RETAINED_ALLOCATION_FAILED;
	}
    }

    void result(BObolRetainedAllocationResult &result) const
    {
	result.maximumNormalizedError = this->candidates.empty() ?
	    (this->havePresentedErrorProof ? 0.0 :
		std::numeric_limits<double>::infinity()) :
	    this->maximumNormalizedError;
	result.maximumProjectedErrorPixels =
	    this->havePresentedErrorProof ?
	    this->maximumPresentedPixelError :
	    std::numeric_limits<double>::infinity();
	result.protectedFloorBudget = this->allocationBudget;
	result.protectedFloorSignature = this->protectedFloorSignature;
	result.selectedPopulationSignature = this->populationSignature();
	result.selectedPresentationCost = this->finalCost;
	result.pixelDemandPresentationCost = this->pixelDemandPresentationCost;
	result.certifiedPresentationBudget = this->effectiveBudget;
	result.allocationPlanSerial = this->allocationPlanSerial;
	result.cadRevision = this->cadRevision;
	result.residentDemandRevision = this->residentDemandRevision;
	result.revisionStamp = this->revisionStamp;
	result.residentAdmissionRevision = this->residentAdmissionRevision;
	result.pointProxyPixelThreshold = this->pointProxyPixelThreshold;
	result.requestedSceneBudget = this->sceneBudget;
	result.externalPresentationCost = this->externalPresentationCost;
	result.fixedCadPresentationCost = this->fixedCadPresentationCost;
	result.maximumMarginalBudget = this->maximumMarginalBudget;
	result.maximumProtectedBudget = this->maximumProtectedBudget;
	result.pointProxyCandidateCount = this->pointProxyCandidateCount;
	for (size_t i = 0; i < this->candidates.size(); ++i) {
	    const Candidate &candidate = this->candidates[i];
	    const bool pointProxy = i < this->finalPointProxies.size() &&
		this->finalPointProxies[i];
	    const int cut = i < this->finalCuts.size() ?
		this->finalCuts[i] : candidate.minimumCut;
	    const double error = pointProxy ? candidate.pointProxyError :
		candidate.mesh.projectedErrorAtCut(cut,
		    candidate.projectedPixelDiameter);
	    const double normalizedError = bobol_lod_normalized_visual_error(
		error, candidate.targetPixelError);
	    if (pointProxy)
		result.selectedPointProxyCount++;
	    const bool prominent = bobol_lod_visual_prominent(
		candidate.visualFootprint);
	    if (prominent)
		result.prominentCandidateCount++;
	    if (prominent && error > candidate.protectedErrorPixels)
		result.prominentQualityFloorViolationCount++;
	    if (std::isfinite(normalizedError))
		result.visualImportanceDebt += candidate.visualFootprint *
		    std::max(0.0, normalizedError - 1.0);
	}
	result.allowProtectedFloor = this->allowProtectedFloor;
    }

    void trace(void) const
    {
	const int64_t discovery = this->phaseMicroseconds[DISCOVERY] +
	    this->phaseMicroseconds[PROTECTION_COMPARE];
	const int64_t marginal = this->phaseMicroseconds[MARGINAL_SEED] +
	    this->phaseMicroseconds[MARGINAL_APPLY] +
	    this->phaseMicroseconds[PROOF];
	int64_t total = 0;
	for (int i = 0; i < PHASE_COUNT; ++i)
	    total += this->phaseMicroseconds[i];
	bu_log("BObol LoD retained allocator phases candidates=%zu "
	       "scan_us=%lld floor_us=%lld baseline_us=%lld "
	       "marginal_us=%lld publish_us=%lld total_us=%lld wall_us=%lld\n",
	       this->candidates.size(),
	       static_cast<long long>(discovery),
	       static_cast<long long>(this->phaseMicroseconds[FLOOR]),
	       static_cast<long long>(this->phaseMicroseconds[BASELINE]),
	       static_cast<long long>(marginal),
	       static_cast<long long>(this->phaseMicroseconds[COMMIT]),
	       static_cast<long long>(total),
	       static_cast<long long>(bu_gettime() -
		   this->wallStartedMicroseconds));
    }

private:
    static void hashPopulationValue(uint64_t &hash, uint64_t value)
    {
	static constexpr uint64_t prime = UINT64_C(1099511628211);
	for (size_t byte = 0; byte < sizeof(value); ++byte) {
	    hash ^= (value >> (byte * 8u)) & UINT64_C(0xff);
	    hash *= prime;
	}
    }

    uint64_t populationSignature(void) const
    {
	static constexpr uint64_t offset = UINT64_C(1469598103934665603);
	uint64_t hash = offset;
	this->hashPopulationValue(hash, this->fixedCadPresentationCost);
	this->hashPopulationValue(hash, this->fixedPointAllocations.size());
	for (const FixedPointAllocation &fixed : this->fixedPointAllocations) {
	    this->hashPopulationValue(hash,
		static_cast<uint64_t>(reinterpret_cast<uintptr_t>(
		    fixed.payload)));
	    this->hashPopulationValue(hash,
		static_cast<uint64_t>(static_cast<uint32_t>(fixed.cut)));
	    this->hashPopulationValue(hash,
		static_cast<uint64_t>(static_cast<uint32_t>(fixed.drawMode)));
	}
	this->hashPopulationValue(hash, this->candidates.size());
	for (size_t i = 0; i < this->candidates.size(); ++i) {
	    const Candidate &candidate = this->candidates[i];
	    this->hashPopulationValue(hash,
		static_cast<uint64_t>(reinterpret_cast<uintptr_t>(
		    candidate.payload)));
	    this->hashPopulationValue(hash,
		static_cast<uint64_t>(static_cast<uint32_t>(
		    this->finalCuts[i])));
	    this->hashPopulationValue(hash,
		static_cast<uint64_t>(static_cast<uint32_t>(
		    candidate.drawMode)));
	    this->hashPopulationValue(hash,
		this->finalPointProxies[i] ? UINT64_C(1) : UINT64_C(0));
	}
	/* Zero denotes an unavailable identity at the policy boundary. */
	return hash ? hash : UINT64_C(1);
    }

    /* A feature below this footprint is normally either small hardware or a
     * thin projected detail.  It still participates in marginal refinement,
     * but does not make the all-or-nothing prominent-feature floor fail for
     * an entire vehicle view. */
    static constexpr size_t protectedFloorTrialDivisor = 4;
    struct Candidate {
	const BObolViewLodState::CadPayload *payload = NULL;
	BObolLodProgressiveMeshSnapshot mesh;
	int minimumCut = -1;
	int maximumCut = -1;
	int drawMode = BOBOL_LOD_DRAW_UNKNOWN;
	SbBool hasNormals = FALSE;
	double projectedPixelDiameter = 0.0;
	double visualFootprint = 0.0;
	double pointProxyError = 0.0;
	size_t pointProxyCost = 0;
	double targetPixelError = 1.0;
	double protectedErrorPixels =
	    std::numeric_limits<double>::infinity();
	unsigned int visualEmphasis = BOBOL_LOD_VISUAL_ORDINARY;
	int protectedCut = -1;
	SoCADAssembly *assembly = NULL;
	Obol::InstanceId instance;
	bool hasInstance = false;
	bool pointProxyEligible = false;
	std::shared_ptr<const
	    BObolViewLodState::CadPayload::ProjectedCutCounts> visibleCounts;
    };

    /* Pixel-exact aggregate points are complete allocation records, but they
     * do not request a mesh cut.  Keep their identity in the population
     * signature without staging a minimum mesh prefix: doing both marks every
     * aggregate point as an unapplied mesh mismatch and expands a sparse
     * retained update back to the complete source population. */
    struct FixedPointAllocation {
	const BObolViewLodState::CadPayload *payload = NULL;
	int cut = -1;
	int drawMode = BOBOL_LOD_DRAW_UNKNOWN;
    };

    struct SourceSnapshot {
	SoBRLDatabaseSource *source = NULL;
	SoCADAssembly *assembly = NULL;
	int drawMode = BOBOL_LOD_DRAW_UNKNOWN;
	std::vector<const BObolViewLodState::CadPayload *> payloads;
    };

    struct ProtectedSet {
	std::vector<Obol::InstanceId> current;
	std::unordered_set<Obol::InstanceId,
	    std::hash<Obol::InstanceId>> next;
	size_t compareCursor = 0;
	bool changed = false;
    };

    struct MarginalUpgradeLess {
	bool operator()(const BObolRetainedMarginalUpgrade &a,
		const BObolRetainedMarginalUpgrade &b) const
	{
	    return bobol_retained_marginal_lower_priority(a, b);
	}
    };

    enum Phase {
	DISCOVERY = 0,
	PROTECTION_COMPARE,
	FLOOR,
	BASELINE,
	MARGINAL_SEED,
	MARGINAL_APPLY,
	PROOF,
	COMMIT,
	COMPLETE,
	PHASE_COUNT
    };

    static double visualFootprint(
	const BObolViewLodState::CadPayload *payload)
    {
	if (!payload)
	    return 0.0;
	return bobol_lod_visual_footprint(payload->projectedPixelArea,
	    payload->projectedPixelPerimeter,
	    payload->projectedPixelDiameter);
    }

    bool pointProxyEligible(const SourceSnapshot &source,
	const BObolViewLodState::CadPayload *payload) const
    {
	if (!source.source || !payload || !payload->progressiveMesh ||
	    !payload->progressiveMesh->isValid() ||
	    !std::isfinite(this->pointProxyPixelThreshold) ||
	    this->pointProxyPixelThreshold <= 0.0f ||
	    !std::isfinite(payload->projectedPixelDiameter) ||
	    payload->projectedPixelDiameter <= 0.0f ||
	    !payload->projectedBoundsContained ||
	    payload->visualEmphasis >= 2 || !source.assembly)
	    return false;
	return payload->projectedPixelDiameter <=
	    this->pointProxyPixelThreshold * 0.75f;
    }

    bool pointProxyIsPixelExact(const SourceSnapshot &source,
	const BObolViewLodState::CadPayload *payload) const
    {
	return this->pointProxyEligible(source, payload) &&
	    visualFootprint(payload) <= 1.0;
    }

    static size_t aggregateProxyCost(
	const BObolViewLodState::CadPayload *payload, int drawMode)
    {
	const SbBool box = payload && payload->projectedPixelDiameter >
	    Obol::CadMaximumPointProxyExtentPixels ? TRUE : FALSE;
	return bobol_lod_aggregate_proxy_render_cost(box, drawMode);
    }

    BObolLodCounts countsAtCut(const Candidate &candidate, int cut) const
    {
	return candidate.visibleCounts ?
	    (*candidate.visibleCounts)[static_cast<size_t>(cut)] :
	    candidate.mesh.hierarchyCountsAtCut(cut, candidate.hasNormals);
    }

    void discover(const SourceSnapshot &source,
	const BObolViewLodState::CadPayload *payload)
    {
	if (!payload || !payload->isValid())
	    return;
	const bool progressive = payload->progressiveMesh &&
	    payload->progressiveMesh->isValid();
	const bool exact = payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL;
	/* Exact geometry is view independent and has no allocator-owned cut.  A
	 * compact source may retire its LoD target inventory after publishing that
	 * geometry, so no later submission pass exists merely to refresh projection
	 * metadata.  Count it as fixed presentation for every view.  Requiring its
	 * old view/policy stamp here manufactures an empty allocation certificate
	 * for a populated exact scene and makes stable handoff retry forever.
	 * Proxies still require current projection evidence because their visual
	 * error is view dependent. */
	if (!progressive &&
	    (exact || (payload->viewRevision == this->viewRevision() &&
		payload->policyRevision == this->policyRevision()))) {
	    this->havePresentedErrorProof =
		this->havePresentedErrorProof || exact;
	    const size_t cost = bobol_lod_render_cost_units(
		payload->counts, source.drawMode, 1);
	    this->fixedCadPresentationCost = cost >
		    SIZE_MAX - this->fixedCadPresentationCost ? SIZE_MAX :
		this->fixedCadPresentationCost + cost;
	    this->fixedCost = cost > SIZE_MAX - this->fixedCost ?
		SIZE_MAX : this->fixedCost + cost;
	    this->pixelDemandPresentationCost = cost >
		    SIZE_MAX - this->pixelDemandPresentationCost ? SIZE_MAX :
		this->pixelDemandPresentationCost + cost;
	    return;
	}
	if (!progressive || payload->viewRevision != this->viewRevision() ||
	    payload->policyRevision != this->policyRevision())
	    return;
	const BObolLodProgressiveMeshSnapshot mesh =
	    payload->progressiveMesh->snapshot();
	if (!mesh.isValid())
	    return;
	const int minimumCut = mesh.minimumCut();
	const int requestedMaximumCut = std::max(minimumCut,
	    std::min(payload->requestedCut,
		mesh.maximumCut()));
	/* A current memory denial proves that the unavailable suffix cannot be
	 * made resident in this admission epoch.  It does not make the retained
	 * active prefix invalid.  Restrict this allocation to that prefix while
	 * retaining the full view demand below for diagnostics and for the next
	 * capacity epoch. */
	const bool currentMemoryDenial = payload->memoryLimited &&
	    payload->residentAdmissionRevision != 0 &&
	    payload->residentAdmissionRevision ==
		this->residentAdmissionRevision &&
	    requestedMaximumCut > payload->activeCut;
	const int maximumCut = currentMemoryDenial ?
	    std::max(minimumCut,
		std::min(requestedMaximumCut, payload->activeCut)) :
	    requestedMaximumCut;
	if (this->pointProxyIsPixelExact(source, payload)) {
	    this->havePresentedErrorProof = true;
	    this->maximumPresentedPixelError = std::max(
		this->maximumPresentedPixelError,
		std::min(static_cast<double>(this->pointProxyPixelThreshold),
		    std::max(0.0, visualFootprint(payload))));
	    BObolLodCounts point;
	    point.pointCount = 1;
	    const size_t cost = bobol_lod_render_cost_units(
		point, source.drawMode, 0);
	    this->fixedCost = cost > SIZE_MAX - this->fixedCost ?
		SIZE_MAX : this->fixedCost + cost;
	    this->pixelDemandPresentationCost = cost >
		    SIZE_MAX - this->pixelDemandPresentationCost ? SIZE_MAX :
		this->pixelDemandPresentationCost + cost;
	    FixedPointAllocation fixed;
	    fixed.payload = payload;
	    fixed.cut = minimumCut;
	    fixed.drawMode = source.drawMode;
	    this->fixedPointAllocations.push_back(std::move(fixed));
	    return;
	}

	/* The source request has already sanitized this value.  In particular, a
	 * fractional target is part of the view contract, not a hint that may be
	 * rounded up to one pixel by the scene-wide allocator. */
	const double target = std::max(
	    static_cast<double>(std::numeric_limits<float>::min()),
	    static_cast<double>(payload->targetPixelError));
	const double diameter = std::max(0.0,
	    static_cast<double>(payload->projectedPixelDiameter));
	const double footprint = visualFootprint(payload);
	Candidate candidate;
	candidate.payload = payload;
	candidate.mesh = mesh;
	candidate.minimumCut = minimumCut;
	candidate.maximumCut = maximumCut;
	candidate.drawMode = source.drawMode;
	candidate.hasNormals = payload->hasNormals;
	candidate.projectedPixelDiameter = diameter;
	candidate.visualFootprint = footprint;
	candidate.pointProxyError = std::max(footprint, diameter);
	candidate.pointProxyCost = aggregateProxyCost(payload, source.drawMode);
	candidate.targetPixelError = target;
	candidate.visualEmphasis = bobol_lod_visual_emphasis(
	    payload->visualEmphasis);
	candidate.assembly = source.assembly;
	candidate.pointProxyEligible =
	    this->pointProxyEligible(source, payload);
	if (source.assembly) {
	    BObolCompactInstanceHandle handle;
	    size_t sourceEntryIndex = static_cast<size_t>(
		payload->sourceEntryIndex);
	    if ((payload->sourcePopulationEpoch !=
		    source.source->getCompactPopulationEpoch() ||
		 payload->sourceEntryIndex == UINT32_MAX) &&
		!source.source->getCompactInstanceIndex(
		    payload->sourceInstanceKey.getString(), sourceEntryIndex))
		sourceEntryIndex = SIZE_MAX;
	    if (sourceEntryIndex <= static_cast<size_t>(
		    std::numeric_limits<int>::max()) &&
		source.source->getCompactInstanceHandle(
		    static_cast<int>(sourceEntryIndex), handle) &&
		handle.isValid()) {
		candidate.instance.w0 = handle.instanceWord0;
		candidate.instance.w1 = handle.instanceWord1;
		candidate.hasInstance = true;
	    }
	}
	candidate.pointProxyEligible = candidate.pointProxyEligible &&
	    candidate.hasInstance;
	if (candidate.pointProxyEligible)
	    this->pointProxyCandidateCount++;
	if (payload->projectedCutCounts &&
	    payload->projectedCutCountsViewRevision == this->viewRevision() &&
	    payload->projectedCutCountsPolicyRevision == this->policyRevision() &&
	    payload->projectedCutCountsMeshRevision == candidate.mesh.revision())
	    candidate.visibleCounts = payload->projectedCutCounts;
	const size_t demandedCost = bobol_lod_render_cost_units(
	    this->countsAtCut(candidate, requestedMaximumCut),
	    candidate.drawMode, 1);
	this->pixelDemandPresentationCost = demandedCost >
		SIZE_MAX - this->pixelDemandPresentationCost ? SIZE_MAX :
	    this->pixelDemandPresentationCost + demandedCost;
	/* Use the same recognizable-feature contract as the production
	 * diagnostics.  Divergent thresholds made the allocator jump over a cut
	 * which the visual oracle considered complete, reject the richer cut on
	 * frame cost, and then fall substantially below recognizable quality.  The
	 * producer request remains the lower bound, and a measured budget may still
	 * reject this atomic candidate floor; the marginal queue then preserves the
	 * same feature priority. */
	const double protectedError = bobol_lod_protected_visual_error(
	    candidate.visualEmphasis, footprint, target);
	candidate.protectedErrorPixels = protectedError;
	candidate.protectedCut = minimumCut;
	if (std::isfinite(protectedError) && diameter > 0.0) {
	    const int protectedCut = candidate.mesh.cutForScreenError(
		diameter, protectedError);
	    if (protectedCut >= 0)
		candidate.protectedCut = std::max(minimumCut,
		    std::min(maximumCut, protectedCut));
	}
	const size_t protectedCost = bobol_lod_render_cost_units(
	    this->countsAtCut(candidate, candidate.protectedCut),
	    candidate.drawMode, 1);
	this->protectedCandidateCost = protectedCost >
	    SIZE_MAX - this->protectedCandidateCost ? SIZE_MAX :
	    this->protectedCandidateCost + protectedCost;
	uint64_t token = static_cast<uint64_t>(
	    reinterpret_cast<uintptr_t>(candidate.payload));
	token ^= static_cast<uint64_t>(
	    static_cast<uint32_t>(candidate.protectedCut)) << 32;
	token ^= static_cast<uint64_t>(
	    static_cast<uint32_t>(candidate.drawMode)) << 16;
	token ^= candidate.hasNormals ? 0xd6e8feb86659fd93ULL : 0;
	token += 0x9e3779b97f4a7c15ULL;
	token = (token ^ (token >> 30)) * 0xbf58476d1ce4e5b9ULL;
	token = (token ^ (token >> 27)) * 0x94d049bb133111ebULL;
	token ^= token >> 31;
	this->protectedFloorSignature ^= token;
	if (candidate.pointProxyEligible) {
	    this->optionalPointCandidateCost = protectedCost >
		SIZE_MAX - this->optionalPointCandidateCost ? SIZE_MAX :
		this->optionalPointCandidateCost + protectedCost;
	    this->optionalPointCandidateSignature ^= token;
	}
	this->candidates.push_back(std::move(candidate));
	this->hasOptionalPointCandidates =
	    this->hasOptionalPointCandidates ||
	    this->candidates.back().pointProxyEligible;
    }

    void beginBaseline(void)
    {
	this->finalCuts.assign(this->candidates.size(), -1);
	this->finalCosts.assign(this->candidates.size(), 0);
	this->finalPointProxies.assign(this->candidates.size(), false);
	this->finalCursor = 0;
	this->finalCost = this->fixedCost;
	this->fixedMaximumPresentedPixelError =
	    this->maximumPresentedPixelError;
	this->fixedPresentedErrorProof = this->havePresentedErrorProof;
	this->maximumNormalizedError = 1.0;
    }

    void observeError(const Candidate &candidate, int cut)
    {
	const double error = candidate.mesh.projectedErrorAtCut(
	    cut, candidate.projectedPixelDiameter);
	if (std::isfinite(error)) {
	    this->havePresentedErrorProof = true;
	    this->maximumPresentedPixelError = std::max(
		this->maximumPresentedPixelError, error);
	}
	const double normalizedError = bobol_lod_normalized_visual_error(
	    error, candidate.targetPixelError);
	if (std::isfinite(normalizedError))
	    this->maximumNormalizedError = std::max(
		this->maximumNormalizedError, normalizedError);
    }

    void observePointProxyError(const Candidate &candidate)
    {
	const double error = candidate.pointProxyError;
	if (std::isfinite(error)) {
	    this->havePresentedErrorProof = true;
	    this->maximumPresentedPixelError = std::max(
		this->maximumPresentedPixelError, error);
	}
	const double normalizedError = bobol_lod_normalized_visual_error(
	    error, candidate.targetPixelError);
	if (std::isfinite(normalizedError))
	    this->maximumNormalizedError = std::max(
		this->maximumNormalizedError, normalizedError);
    }

    void protectPointProxyCandidate(const Candidate &candidate)
    {
	if (candidate.assembly && candidate.hasInstance)
	    this->pointProxyProtectedInstances[candidate.assembly].next.insert(
		candidate.instance);
    }

    void queueNextUpgrade(size_t candidateIndex)
    {
	const Candidate &candidate = this->candidates[candidateIndex];
	if (this->finalPointProxies[candidateIndex]) {
	    const size_t nextCost = bobol_lod_render_cost_units(
		this->countsAtCut(candidate, candidate.minimumCut),
		candidate.drawMode, 1);
	    const double nextError = candidate.mesh.projectedErrorAtCut(
		candidate.minimumCut, candidate.projectedPixelDiameter);
	    BObolRetainedMarginalUpgrade upgrade;
	    upgrade.candidateIndex = candidateIndex;
	    upgrade.nextCut = candidate.minimumCut;
	    upgrade.nextCost = nextCost;
	    upgrade.addedCost = nextCost > this->finalCosts[candidateIndex] ?
		nextCost - this->finalCosts[candidateIndex] : 0;
	    upgrade.visualEmphasis = candidate.visualEmphasis;
	    upgrade.qualityFloorViolation =
		candidate.pointProxyError > candidate.protectedErrorPixels;
	    upgrade.normalizedError = bobol_lod_normalized_visual_error(
		candidate.pointProxyError, candidate.targetPixelError);
	    const double benefit = bobol_lod_marginal_visual_benefit(
		candidate.visualFootprint, candidate.pointProxyError, nextError,
		candidate.targetPixelError);
	    upgrade.visualBenefitPerCost = upgrade.addedCost > 0 ?
		benefit / static_cast<double>(upgrade.addedCost) : benefit + 1.0;
	    this->upgrades.push(upgrade);
	    return;
	}
	const int currentCut = this->finalCuts[candidateIndex];
	if (currentCut < candidate.minimumCut ||
	    currentCut >= candidate.maximumCut)
	    return;
	const double currentError = candidate.mesh.projectedErrorAtCut(
	    currentCut, candidate.projectedPixelDiameter);
	for (int nextCut = currentCut + 1;
	     nextCut <= candidate.maximumCut; ++nextCut) {
	    const size_t nextCost = bobol_lod_render_cost_units(
		this->countsAtCut(candidate, nextCut),
		candidate.drawMode, 1);
	    const double nextError = candidate.mesh.projectedErrorAtCut(
		nextCut, candidate.projectedPixelDiameter);
	    if (nextCost <= this->finalCosts[candidateIndex] &&
		!(nextError < currentError))
		continue;
	    BObolRetainedMarginalUpgrade upgrade;
	    upgrade.candidateIndex = candidateIndex;
	    upgrade.nextCut = nextCut;
	    upgrade.nextCost = nextCost;
	    upgrade.addedCost = nextCost > this->finalCosts[candidateIndex] ?
		nextCost - this->finalCosts[candidateIndex] : 0;
	    upgrade.visualEmphasis = candidate.visualEmphasis;
	    upgrade.qualityFloorViolation =
		currentError > candidate.protectedErrorPixels;
	    upgrade.normalizedError = bobol_lod_normalized_visual_error(
		currentError, candidate.targetPixelError);
	    const double benefit = bobol_lod_marginal_visual_benefit(
		candidate.visualFootprint, currentError, nextError,
		candidate.targetPixelError);
	    upgrade.visualBenefitPerCost = upgrade.addedCost > 0 ?
		benefit / static_cast<double>(upgrade.addedCost) :
		benefit + 1.0;
	    this->upgrades.push(upgrade);
	    return;
	}
    }

    bool stageAllocation(const BObolViewLodState::CadPayload *target,
	int cut, int drawMode)
    {
	return this->viewState && this->viewState->stageCadAllocatedCut(
	    target, cut, this->viewRevision(), this->policyRevision(), drawMode,
	    this->allocationPlanSerial);
    }

    void commit(void)
    {
	for (SoCADAssembly *assembly : this->pointProxyAssemblies) {
	    ProtectedSet &protection =
		this->pointProxyProtectedInstances[assembly];
	    if (protection.changed)
		assembly->adoptPointProxyProtectedInstances(
		    std::move(protection.next));
	}
    }

    BObolViewLodState *viewState = NULL;
    BObolRetainedAllocationInputKey inputKey;
    size_t externalPresentationCost = 0;
    size_t fixedCadPresentationCost = 0;
    size_t pixelDemandPresentationCost = 0;
    uint64_t viewRevision(void) const
    {
	return this->revisionStamp.view.value();
    }

    uint64_t policyRevision(void) const
    {
	return this->revisionStamp.policy.value();
    }

    BObolLodAdmissionRevisionStamp revisionStamp;
    uint64_t residentAdmissionRevision = 0;
    uint64_t cadRevision = 0;
    uint64_t residentDemandRevision = 0;
    uint64_t allocationPlanSerial = 0;
    size_t sceneBudget = 0;
    size_t maximumMarginalBudget = 0;
    bool allowProtectedFloor = false;
    size_t maximumProtectedBudget = 0;
    size_t pointProxyCandidateCount = 0;
    float pointProxyPixelThreshold = 0.0f;
    std::vector<SourceSnapshot> sources;
    size_t sourceCursor = 0;
    size_t payloadCursor = 0;
    std::vector<Candidate> candidates;
    std::vector<FixedPointAllocation> fixedPointAllocations;
    std::unordered_map<SoCADAssembly *, ProtectedSet>
	pointProxyProtectedInstances;
    std::vector<SoCADAssembly *> pointProxyAssemblies;
    size_t pointProxyAssemblyCursor = 0;
    size_t fixedCost = 0;
    size_t protectedCandidateCost = 0;
    size_t optionalPointCandidateCost = 0;
    uint64_t protectedFloorSignature = 0;
    uint64_t optionalPointCandidateSignature = 0;
    double maximumPresentedPixelError = 0.0;
    bool havePresentedErrorProof = false;
    size_t protectedFloorCost = 0;
    bool enforceProtectedFloor = false;
    size_t effectiveBudget = 0;
    size_t allocationBudget = 0;
    bool hasOptionalPointCandidates = false;
    bool useOptionalPointCandidates = false;
    std::vector<int> finalCuts;
    std::vector<size_t> finalCosts;
    std::vector<bool> finalPointProxies;
    size_t finalCursor = 0;
    size_t finalCost = 0;
    double fixedMaximumPresentedPixelError = 0.0;
    bool fixedPresentedErrorProof = false;
    double maximumNormalizedError = 1.0;
    bool marginalAccepted = false;
    std::priority_queue<BObolRetainedMarginalUpgrade,
	std::vector<BObolRetainedMarginalUpgrade>, MarginalUpgradeLess> upgrades;
    Phase phase = DISCOVERY;
    int64_t phaseMicroseconds[PHASE_COUNT] = {0};
    int64_t wallStartedMicroseconds = 0;
};

BObolRetainedAllocationStatus
bobol_retained_allocation_advance(
    std::shared_ptr<BObolRetainedAllocationTransaction> &transaction,
    const BObolRetainedAllocationInputs &inputs,
    uint64_t sliceMicroseconds,
    BObolRetainedAllocationResult &result)
{
    result = BObolRetainedAllocationResult();
    if (!inputs.viewState || inputs.sceneBudget == SIZE_MAX) {
	transaction.reset();
	return BOBOL_RETAINED_ALLOCATION_COMPLETE;
    }
    try {
	if (!transaction || !transaction->matches(inputs))
	    transaction = std::make_shared<
		BObolRetainedAllocationTransaction>(inputs);
    } catch (const std::bad_alloc &) {
	transaction.reset();
	return BOBOL_RETAINED_ALLOCATION_FAILED;
    }
    const BObolRetainedAllocationStatus status =
	transaction->advance(inputs, sliceMicroseconds);
    if (status == BOBOL_RETAINED_ALLOCATION_COMPLETE)
	transaction->result(result);
    else if (status == BOBOL_RETAINED_ALLOCATION_STALE ||
	status == BOBOL_RETAINED_ALLOCATION_FAILED)
	/* A failed or stale transaction can never make forward progress.  In
	 * particular, another allocation-plan serial may invalidate COMMIT without
	 * changing CAD or resident-demand revisions.  Retaining that object makes
	 * later advances retry the same impossible publication forever. */
	transaction.reset();
    return status;
}

void
bobol_retained_allocation_trace(
    const std::shared_ptr<BObolRetainedAllocationTransaction> &transaction)
{
    if (transaction)
	transaction->trace();
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
