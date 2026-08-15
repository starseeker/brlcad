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

class BObolRetainedAllocationTransaction {
public:
    explicit BObolRetainedAllocationTransaction(
	const BObolRetainedAllocationInputs &inputs)
    {
	this->viewState = inputs.viewState;
	this->viewRevision = inputs.viewRevision;
	this->policyRevision = inputs.policyRevision;
	this->cadRevision = inputs.viewState ?
	    inputs.viewState->cadRevision() : 0;
	this->residentDemandRevision = inputs.viewState ?
	    inputs.viewState->residentMeshDemandRevision() : 0;
	this->allocationPlanSerial = inputs.viewState ?
	    inputs.viewState->beginCadAllocationPlan() : 0;
	this->sceneBudget = inputs.sceneBudget;
	this->allowProtectedFloor = inputs.allowProtectedFloor;
	this->maximumProtectedBudget = inputs.maximumProtectedBudget;
	this->progressiveCutCeiling = inputs.progressiveCutCeiling;
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
		    this->viewState->findCadPayloadsUnordered(
			source, snapshot.payloads);
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
	    this->viewRevision != inputs.viewRevision ||
	    this->policyRevision != inputs.policyRevision ||
	    this->cadRevision != inputs.viewState->cadRevision() ||
	    this->residentDemandRevision !=
		inputs.viewState->residentMeshDemandRevision() ||
	    this->sceneBudget != inputs.sceneBudget ||
	    this->allowProtectedFloor != inputs.allowProtectedFloor ||
	    this->maximumProtectedBudget != inputs.maximumProtectedBudget ||
	    this->progressiveCutCeiling != inputs.progressiveCutCeiling ||
	    std::memcmp(&this->pointProxyPixelThreshold,
		&inputs.pointProxyPixelThreshold, sizeof(float)) != 0 ||
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
	auto beginPopulationCost = [this](double targetCeiling) {
	    this->populationCeiling = targetCeiling;
	    this->populationCursor = 0;
	    this->populationCost = this->fixedCost;
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
			transition(PROTECTION_COMPARE);
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
			transition(FLOOR);
			break;
		    case FLOOR: {
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
			    (protectedBudget > SIZE_MAX - protectedBudget / 4 ?
				SIZE_MAX : protectedBudget + protectedBudget / 4) :
			    protectedBudget;
			this->enforceProtectedFloor = this->allowProtectedFloor &&
			    this->protectedFloorCost <= trialBudget;
			this->effectiveBudget = this->enforceProtectedFloor ?
			    std::max(this->sceneBudget,
				this->protectedFloorCost) : this->sceneBudget;
			this->allocationBudget = this->enforceProtectedFloor ?
			    this->protectedFloorCost : 0;
			if (this->candidates.empty()) {
			    this->beginBaseline();
			    transition(BASELINE);
			    break;
			}
			beginPopulationCost(1.0);
			transition(SEARCH_ONE);
			break;
		    }
		    case SEARCH_ONE:
		    case SEARCH_MAXIMUM:
		    case SEARCH_BINARY:
			while (this->populationCursor <
				this->candidates.size()) {
			    const Candidate &candidate =
				this->candidates[this->populationCursor++];
			    const int cut = this->allocatedCut(candidate,
				this->populationCeiling);
			    const size_t cost = bobol_lod_render_cost_units(
				this->countsAtCut(candidate, cut),
				candidate.drawMode, 1);
			    this->populationCost = cost >
				    SIZE_MAX - this->populationCost ? SIZE_MAX :
				this->populationCost + cost;
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			if (this->phase == SEARCH_ONE) {
			    if (this->populationCost <= this->effectiveBudget) {
				this->ceiling = 1.0;
				this->beginBaseline();
				transition(BASELINE);
			    } else {
				beginPopulationCost(this->maximumMinimumError);
				transition(SEARCH_MAXIMUM);
			    }
			    break;
			}
			if (this->phase == SEARCH_MAXIMUM) {
			    if (this->populationCost > this->effectiveBudget) {
				this->ceiling = this->maximumMinimumError;
				this->beginBaseline();
				transition(BASELINE);
			    } else {
				this->lowLog = 0.0;
				this->highLog =
				    std::log2(this->maximumMinimumError);
				this->binaryIteration = 0;
				beginPopulationCost(std::exp2(
				    (this->lowLog + this->highLog) * 0.5));
				transition(SEARCH_BINARY);
			    }
			    break;
			}
			if (this->populationCost <= this->effectiveBudget)
			    this->highLog = std::log2(this->populationCeiling);
			else
			    this->lowLog = std::log2(this->populationCeiling);
			this->binaryIteration++;
			if (this->binaryIteration >= 14) {
			    this->ceiling = std::exp2(this->highLog) *
				(1.0 + 1.0e-9);
			    this->beginBaseline();
			    transition(BASELINE);
			} else {
			    beginPopulationCost(std::exp2(
				(this->lowLog + this->highLog) * 0.5));
			}
			break;
		    case BASELINE:
			while (this->fixedStageCursor <
				this->fixedAllocations.size()) {
			    const FixedAllocation &fixed =
				this->fixedAllocations[this->fixedStageCursor++];
			    if (!this->stageAllocation(
				    fixed.payload, fixed.cut, fixed.drawMode)) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_FAILED;
			    }
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			while (this->finalCursor < this->candidates.size()) {
			    const size_t i = this->finalCursor++;
			    const Candidate &candidate = this->candidates[i];
			    this->finalCuts[i] =
				this->allocatedCut(candidate, this->ceiling);
			    this->finalCosts[i] = bobol_lod_render_cost_units(
				this->countsAtCut(candidate, this->finalCuts[i]),
				candidate.drawMode, 1);
			    this->finalCost = this->finalCosts[i] >
				    SIZE_MAX - this->finalCost ? SIZE_MAX :
				this->finalCost + this->finalCosts[i];
			    this->observeError(candidate, this->finalCuts[i]);
			    if (!this->stageAllocation(candidate.payload,
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
			    transition(COMMIT);
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
			    const MarginalUpgrade upgrade = this->upgrades.top();
			    this->upgrades.pop();
			    if (upgrade.addedCost <=
				    this->effectiveBudget - this->finalCost) {
				this->finalCost += upgrade.addedCost;
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
				this->marginalAccepted = true;
				this->queueNextUpgrade(upgrade.candidateIndex);
			    }
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			if (this->marginalAccepted) {
			    this->finalCursor = 0;
			    this->realizedCeiling = 1.0;
			    this->maximumPresentedPixelError =
				this->fixedMaximumPresentedPixelError;
			    this->havePresentedErrorProof =
				this->fixedPresentedErrorProof;
			    transition(PROOF);
			} else {
			    transition(COMMIT);
			}
			break;
		    case PROOF:
			while (this->finalCursor < this->candidates.size()) {
			    const size_t i = this->finalCursor++;
			    this->observeError(
				this->candidates[i], this->finalCuts[i]);
			    if (shouldYield()) {
				finishSlice();
				return BOBOL_RETAINED_ALLOCATION_PENDING;
			    }
			}
			transition(COMMIT);
			break;
		    case COMMIT:
			if (!this->matches(inputs) ||
			    !this->viewState->commitCadAllocationPlan(
				this->allocationPlanSerial, this->cadRevision,
				this->residentDemandRevision)) {
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
	result.normalizedError = this->candidates.empty() ?
	    std::numeric_limits<double>::infinity() :
	    this->realizedCeiling;
	result.maximumProjectedErrorPixels =
	    this->havePresentedErrorProof ?
	    this->maximumPresentedPixelError :
	    std::numeric_limits<double>::infinity();
	result.protectedFloorBudget = this->allocationBudget;
	result.protectedFloorSignature = this->protectedFloorSignature;
    }

    void trace(void) const
    {
	const int64_t discovery = this->phaseMicroseconds[DISCOVERY] +
	    this->phaseMicroseconds[PROTECTION_COMPARE];
	const int64_t search = this->phaseMicroseconds[SEARCH_ONE] +
	    this->phaseMicroseconds[SEARCH_MAXIMUM] +
	    this->phaseMicroseconds[SEARCH_BINARY];
	const int64_t marginal = this->phaseMicroseconds[MARGINAL_SEED] +
	    this->phaseMicroseconds[MARGINAL_APPLY] +
	    this->phaseMicroseconds[PROOF];
	int64_t total = 0;
	for (int i = 0; i < PHASE_COUNT; ++i)
	    total += this->phaseMicroseconds[i];
	bu_log("BObol LoD retained allocator phases candidates=%zu "
	       "scan_us=%lld floor_us=%lld search_us=%lld baseline_us=%lld "
	       "marginal_us=%lld publish_us=%lld total_us=%lld wall_us=%lld\n",
	       this->candidates.size(),
	       static_cast<long long>(discovery),
	       static_cast<long long>(this->phaseMicroseconds[FLOOR]),
	       static_cast<long long>(search),
	       static_cast<long long>(this->phaseMicroseconds[BASELINE]),
	       static_cast<long long>(marginal),
	       static_cast<long long>(this->phaseMicroseconds[COMMIT]),
	       static_cast<long long>(total),
	       static_cast<long long>(bu_gettime() -
		   this->wallStartedMicroseconds));
    }

private:
    static constexpr double protectedFootprintPixels = 12.0;

    struct Candidate {
	const BObolViewLodState::CadPayload *payload = NULL;
	BObolLodProgressiveMeshPtr mesh;
	int minimumCut = -1;
	int maximumCut = -1;
	int drawMode = BOBOL_LOD_DRAW_UNKNOWN;
	SbBool hasNormals = FALSE;
	double projectedPixelDiameter = 0.0;
	double errorWeight = 1.0;
	int protectedCut = -1;
	std::shared_ptr<const
	    BObolViewLodState::CadPayload::ProjectedCutCounts> visibleCounts;
    };

    struct FixedAllocation {
	const BObolViewLodState::CadPayload *payload = NULL;
	BObolLodProgressiveMeshPtr mesh;
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

    struct MarginalUpgrade {
	size_t candidateIndex = 0;
	int nextCut = -1;
	size_t nextCost = 0;
	size_t addedCost = 0;
	bool qualityFloorViolation = false;
	double weightedError = 0.0;
	double valuePerCost = 0.0;
    };

    struct MarginalUpgradeLess {
	bool operator()(const MarginalUpgrade &a,
		const MarginalUpgrade &b) const
	{
	    if (a.qualityFloorViolation != b.qualityFloorViolation)
		return !a.qualityFloorViolation;
	    if (a.qualityFloorViolation &&
		(a.weightedError < b.weightedError ||
		 a.weightedError > b.weightedError))
		return a.weightedError < b.weightedError;
	    if (a.valuePerCost < b.valuePerCost ||
		a.valuePerCost > b.valuePerCost)
		return a.valuePerCost < b.valuePerCost;
	    if (a.weightedError < b.weightedError ||
		a.weightedError > b.weightedError)
		return a.weightedError < b.weightedError;
	    return a.candidateIndex > b.candidateIndex;
	}
    };

    enum Phase {
	DISCOVERY = 0,
	PROTECTION_COMPARE,
	FLOOR,
	SEARCH_ONE,
	SEARCH_MAXIMUM,
	SEARCH_BINARY,
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
	return std::max(std::sqrt(std::max(0.0,
		static_cast<double>(payload->projectedPixelArea))),
	    std::max(static_cast<double>(
		    payload->projectedPixelPerimeter) * 0.25,
		static_cast<double>(payload->projectedPixelDiameter) * 0.25));
    }

    bool pointProxyCandidate(const SourceSnapshot &source,
	const BObolViewLodState::CadPayload *payload) const
    {
	if (!source.source || !payload || !payload->progressiveMesh ||
	    !payload->progressiveMesh->isValid() ||
	    !std::isfinite(this->pointProxyPixelThreshold) ||
	    this->pointProxyPixelThreshold <= 0.0f ||
	    !std::isfinite(payload->projectedPixelDiameter) ||
	    payload->projectedPixelDiameter <= 0.0f ||
	    !payload->projectedBoundsContained ||
	    payload->visualEmphasis > 0 ||
	    visualFootprint(payload) >= protectedFootprintPixels)
	    return false;
	return payload->projectedPixelDiameter <=
	    this->pointProxyPixelThreshold * 0.75f;
    }

    BObolLodCounts countsAtCut(const Candidate &candidate, int cut) const
    {
	return candidate.visibleCounts ?
	    (*candidate.visibleCounts)[static_cast<size_t>(cut)] :
	    bobol_lod_progressive_counts(
		candidate.mesh, cut, candidate.hasNormals);
    }

    int allocatedCut(const Candidate &candidate, double target) const
    {
	const double targetError = candidate.errorWeight > 0.0 ?
	    target / candidate.errorWeight : target;
	int cut = candidate.mesh->cutForScreenError(
	    candidate.projectedPixelDiameter, targetError);
	if (cut < 0)
	    cut = candidate.minimumCut;
	cut = std::max(candidate.minimumCut,
	    std::min(candidate.maximumCut, cut));
	return this->enforceProtectedFloor ?
	    std::max(cut, candidate.protectedCut) : cut;
    }

    void discover(const SourceSnapshot &source,
	const BObolViewLodState::CadPayload *payload)
    {
	if (!payload || !payload->isValid() ||
	    payload->viewRevision != this->viewRevision ||
	    payload->policyRevision != this->policyRevision)
	    return;
	if (!payload->progressiveMesh ||
	    !payload->progressiveMesh->isValid()) {
	    this->havePresentedErrorProof = true;
	    const size_t cost = bobol_lod_render_cost_units(
		payload->counts, source.drawMode, 1);
	    this->fixedCost = cost > SIZE_MAX - this->fixedCost ?
		SIZE_MAX : this->fixedCost + cost;
	    return;
	}
	const int minimumCut = payload->progressiveMesh->minimumCut();
	const int maximumCut = std::max(minimumCut,
	    std::min(payload->requestedCut,
		payload->progressiveMesh->maximumCut()));
	if (this->pointProxyCandidate(source, payload)) {
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
	    FixedAllocation fixed;
	    fixed.payload = payload;
	    fixed.mesh = payload->progressiveMesh;
	    fixed.cut = minimumCut;
	    fixed.drawMode = source.drawMode;
	    this->fixedAllocations.push_back(std::move(fixed));
	    return;
	}

	const double emphasis = payload->visualEmphasis >= 2 ? 4.0 :
	    (payload->visualEmphasis == 1 ? 2.0 : 1.0);
	const double target = std::max(1.0,
	    static_cast<double>(payload->targetPixelError));
	const double area = std::max(0.0,
	    static_cast<double>(payload->projectedPixelArea));
	const double perimeter = std::max(0.0,
	    static_cast<double>(payload->projectedPixelPerimeter));
	const double diameter = std::max(0.0,
	    static_cast<double>(payload->projectedPixelDiameter));
	const double footprint = std::max(std::sqrt(area),
	    std::max(perimeter * 0.25, diameter * 0.25));
	const double significance = std::max(1.0,
	    std::min(8.0, footprint * 0.50));
	double errorWeight = emphasis * significance / target;
	if (!std::isfinite(errorWeight) || errorWeight <= 0.0)
	    errorWeight = 1.0;
	Candidate candidate;
	candidate.payload = payload;
	candidate.mesh = payload->progressiveMesh;
	candidate.minimumCut = minimumCut;
	candidate.maximumCut = maximumCut;
	candidate.drawMode = source.drawMode;
	candidate.hasNormals = payload->hasNormals;
	candidate.projectedPixelDiameter = diameter;
	candidate.errorWeight = errorWeight;
	if (payload->projectedCutCounts &&
	    payload->projectedCutCountsViewRevision == this->viewRevision &&
	    payload->projectedCutCountsPolicyRevision == this->policyRevision &&
	    payload->projectedCutCountsMeshRevision == candidate.mesh->revision())
	    candidate.visibleCounts = payload->projectedCutCounts;
	const double protectedError = payload->visualEmphasis >= 2 ? 1.0 :
	    (payload->visualEmphasis == 1 ? 2.0 :
		(footprint >= protectedFootprintPixels ? 3.0 :
		    std::numeric_limits<double>::infinity()));
	candidate.protectedCut = minimumCut;
	if (std::isfinite(protectedError) && diameter > 0.0) {
	    const int protectedCut = candidate.mesh->cutForScreenError(
		diameter, protectedError);
	    if (protectedCut >= 0)
		candidate.protectedCut = std::max(minimumCut,
		    std::min(maximumCut, protectedCut));
	}
	if (std::isfinite(protectedError) && source.assembly) {
	    BObolCompactInstanceHandle handle;
	    if (payload->sourceEntryIndex <= static_cast<uint32_t>(
		    std::numeric_limits<int>::max()) &&
		source.source->getCompactInstanceHandle(
		    static_cast<int>(payload->sourceEntryIndex), handle) &&
		handle.isValid()) {
		Obol::InstanceId instance;
		instance.w0 = handle.instanceWord0;
		instance.w1 = handle.instanceWord1;
		this->pointProxyProtectedInstances[source.assembly].next.insert(
		    instance);
	    }
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
	this->candidates.push_back(std::move(candidate));
	const double minimumError =
	    payload->progressiveMesh->projectedErrorAtCut(
		minimumCut, diameter) * errorWeight;
	if (std::isfinite(minimumError))
	    this->maximumMinimumError = std::max(
		this->maximumMinimumError, minimumError);
    }

    void beginBaseline(void)
    {
	this->finalCuts.assign(this->candidates.size(), -1);
	this->finalCosts.assign(this->candidates.size(), 0);
	this->finalCursor = 0;
	this->fixedStageCursor = 0;
	this->finalCost = this->fixedCost;
	this->fixedMaximumPresentedPixelError =
	    this->maximumPresentedPixelError;
	this->fixedPresentedErrorProof = this->havePresentedErrorProof;
	this->realizedCeiling = 1.0;
    }

    void observeError(const Candidate &candidate, int cut)
    {
	const int presentedCut = this->progressiveCutCeiling >= 0 ?
	    std::max(candidate.minimumCut,
		std::min(cut, this->progressiveCutCeiling)) : cut;
	const double error = candidate.mesh->projectedErrorAtCut(
	    presentedCut, candidate.projectedPixelDiameter);
	if (std::isfinite(error)) {
	    this->havePresentedErrorProof = true;
	    this->maximumPresentedPixelError = std::max(
		this->maximumPresentedPixelError, error);
	}
	const double weightedError = error * candidate.errorWeight;
	if (std::isfinite(weightedError))
	    this->realizedCeiling = std::max(
		this->realizedCeiling, weightedError);
    }

    void queueNextUpgrade(size_t candidateIndex)
    {
	const Candidate &candidate = this->candidates[candidateIndex];
	const int currentCut = this->finalCuts[candidateIndex];
	if (currentCut < candidate.minimumCut ||
	    currentCut >= candidate.maximumCut)
	    return;
	const double currentError = candidate.mesh->projectedErrorAtCut(
	    currentCut, candidate.projectedPixelDiameter);
	for (int nextCut = currentCut + 1;
	     nextCut <= candidate.maximumCut; ++nextCut) {
	    const size_t nextCost = bobol_lod_render_cost_units(
		this->countsAtCut(candidate, nextCut),
		candidate.drawMode, 1);
	    const double nextError = candidate.mesh->projectedErrorAtCut(
		nextCut, candidate.projectedPixelDiameter);
	    if (nextCost <= this->finalCosts[candidateIndex] &&
		!(nextError < currentError))
		continue;
	    MarginalUpgrade upgrade;
	    upgrade.candidateIndex = candidateIndex;
	    upgrade.nextCut = nextCut;
	    upgrade.nextCost = nextCost;
	    upgrade.addedCost = nextCost > this->finalCosts[candidateIndex] ?
		nextCost - this->finalCosts[candidateIndex] : 0;
	    upgrade.qualityFloorViolation =
		candidate.protectedCut > currentCut;
	    const double weightedError =
		currentError * candidate.errorWeight;
	    upgrade.weightedError = std::isfinite(weightedError) ?
		weightedError : std::numeric_limits<double>::max();
	    const double benefit = std::max(0.0,
		(currentError - nextError) * candidate.errorWeight);
	    upgrade.valuePerCost = upgrade.addedCost > 0 ?
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
	    target, cut, this->viewRevision, this->policyRevision, drawMode,
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
    uint64_t viewRevision = 0;
    uint64_t policyRevision = 0;
    uint64_t cadRevision = 0;
    uint64_t residentDemandRevision = 0;
    uint64_t allocationPlanSerial = 0;
    size_t sceneBudget = 0;
    bool allowProtectedFloor = false;
    size_t maximumProtectedBudget = 0;
    int progressiveCutCeiling = -1;
    float pointProxyPixelThreshold = 0.0f;
    std::vector<SourceSnapshot> sources;
    size_t sourceCursor = 0;
    size_t payloadCursor = 0;
    std::vector<Candidate> candidates;
    std::vector<FixedAllocation> fixedAllocations;
    size_t fixedStageCursor = 0;
    std::unordered_map<SoCADAssembly *, ProtectedSet>
	pointProxyProtectedInstances;
    std::vector<SoCADAssembly *> pointProxyAssemblies;
    size_t pointProxyAssemblyCursor = 0;
    size_t fixedCost = 0;
    size_t protectedCandidateCost = 0;
    uint64_t protectedFloorSignature = 0;
    double maximumMinimumError = 1.0;
    double maximumPresentedPixelError = 0.0;
    bool havePresentedErrorProof = false;
    size_t protectedFloorCost = 0;
    bool enforceProtectedFloor = false;
    size_t effectiveBudget = 0;
    size_t allocationBudget = 0;
    double populationCeiling = 1.0;
    size_t populationCursor = 0;
    size_t populationCost = 0;
    double lowLog = 0.0;
    double highLog = 0.0;
    int binaryIteration = 0;
    double ceiling = 1.0;
    std::vector<int> finalCuts;
    std::vector<size_t> finalCosts;
    size_t finalCursor = 0;
    size_t finalCost = 0;
    double fixedMaximumPresentedPixelError = 0.0;
    bool fixedPresentedErrorProof = false;
    double realizedCeiling = 1.0;
    bool marginalAccepted = false;
    std::priority_queue<MarginalUpgrade,
	std::vector<MarginalUpgrade>, MarginalUpgradeLess> upgrades;
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
