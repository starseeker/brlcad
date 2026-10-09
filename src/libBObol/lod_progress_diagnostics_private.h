/*   L O D _ P R O G R E S S _ D I A G N O S T I C S _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_PROGRESS_DIAGNOSTICS_PRIVATE_H
#define LIBBOBOL_LOD_PROGRESS_DIAGNOSTICS_PRIVATE_H

#include "common.h"

#include "BObol/BViewController.h"

#include <algorithm>
#include <cstddef>
#include <cstdint>
#include <limits>

/* Diagnostic interpretation of the currently presented proxy population.
 * This is deliberately downstream of scheduling: UI vocabulary must not
 * become another admission policy. */
class BObolLodProxyReasonClassifier {
public:
    struct Inputs {
	size_t presentedSubpixelOccurrenceCount = 0;
	size_t presentedStructuralBoxCount = 0;
	size_t terminalProxyOccurrenceCount = 0;
	size_t terminalOccurrenceFailureCount = 0;
	size_t missingMeshBudgetBlockedCount = 0;
	bool sourcePreparationPending = false;
	bool visibilityPlanningPending = false;
	bool geometryPreparationPending = false;
	bool rendererPreparationPending = false;
	bool frameBudgetLimited = false;
	bool memoryBudgetLimited = false;
    };

    static BObolLodProxyReasonStatus classify(const Inputs &inputs)
    {
	BObolLodProxyReasonStatus result;
	result.intentionalSubpixelOccurrenceCount =
	    inputs.presentedSubpixelOccurrenceCount;

	size_t unassigned = inputs.presentedStructuralBoxCount;
	/* Terminal OBB payloads are rendered structural boxes, not an
	 * additional proxy population.  Remove them from the exact structural
	 * population before assigning the remaining boxes to another cause. */
	result.frameBudgetOccurrenceCount = std::min(
	    unassigned, inputs.terminalProxyOccurrenceCount);
	unassigned -= result.frameBudgetOccurrenceCount;

	result.terminalFailureOccurrenceCount = std::min(
	    unassigned, inputs.terminalOccurrenceFailureCount);
	unassigned -= result.terminalFailureOccurrenceCount;

	const size_t explicitlyBudgetBlocked = std::min(
	    unassigned, inputs.missingMeshBudgetBlockedCount);
	result.frameBudgetOccurrenceCount = saturatingAdd(
	    result.frameBudgetOccurrenceCount, explicitlyBudgetBlocked);
	unassigned -= explicitlyBudgetBlocked;

	if (inputs.sourcePreparationPending) {
	    result.sourcePreparationOccurrenceCount = unassigned;
	} else if (inputs.visibilityPlanningPending) {
	    result.visibilityPlanningOccurrenceCount = unassigned;
	} else if (inputs.geometryPreparationPending) {
	    result.geometryPreparationOccurrenceCount = unassigned;
	} else if (inputs.rendererPreparationPending) {
	    result.rendererPreparationOccurrenceCount = unassigned;
	} else if (inputs.memoryBudgetLimited) {
	    result.memoryBudgetOccurrenceCount = unassigned;
	} else if (inputs.frameBudgetLimited) {
	    result.frameBudgetOccurrenceCount = saturatingAdd(
		result.frameBudgetOccurrenceCount, unassigned);
	} else {
	    result.unclassifiedOccurrenceCount = unassigned;
	}

	if (result.sourcePreparationOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_SOURCE_PREPARATION;
	if (result.visibilityPlanningOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_VISIBILITY_PLANNING;
	if (result.geometryPreparationOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_GEOMETRY_PREPARATION;
	if (result.rendererPreparationOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_RENDERER_PREPARATION;
	if (result.intentionalSubpixelOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_INTENTIONAL_SUBPIXEL;
	if (result.frameBudgetOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_FRAME_BUDGET;
	if (result.memoryBudgetOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_MEMORY_BUDGET;
	if (result.terminalFailureOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_TERMINAL_FAILURE;
	if (result.unclassifiedOccurrenceCount)
	    result.reasonMask |= BOBOL_LOD_PROXY_REASON_UNCLASSIFIED;
	return result;
    }

private:
    static size_t saturatingAdd(size_t left, size_t right)
    {
	const size_t maximum = (std::numeric_limits<size_t>::max)();
	return right > maximum - left ? maximum : left + right;
    }
};

/* Presentation-only episode history.  Milestones are latched from exact
 * frame observations and can neither schedule work nor alter convergence. */
class BObolLodEpisodeTracker {
public:
    struct Inputs {
	uint64_t episodeRevision = 0;
	int64_t episodeStartMicroseconds = 0;
	int64_t observationMicroseconds = 0;
	bool hasLodState = false;
	bool exactPresentation = false;
	bool populationComplete = false;
	bool proxyPresented = false;
	bool meshPresented = false;
	bool stableView = false;
	size_t structuralProxyCount = 0;
    };

    BObolLodEpisodeStatus observe(const Inputs &inputs)
    {
	if (!inputs.hasLodState) {
	    this->reset();
	    return BObolLodEpisodeStatus();
	}
	if (!this->active ||
	    this->episodeRevision != inputs.episodeRevision)
	    this->begin(inputs);

	const uint64_t elapsed = this->elapsedMilliseconds(
	    inputs.observationMicroseconds);
	this->status.elapsedMilliseconds = elapsed;
	if (inputs.exactPresentation && inputs.proxyPresented &&
	    !this->status.firstProxyReached) {
	    this->status.firstProxyReached = TRUE;
	    this->status.firstProxyMilliseconds = elapsed;
	}
	if (inputs.exactPresentation && inputs.meshPresented &&
	    !this->status.firstMeshReached) {
	    this->status.firstMeshReached = TRUE;
	    this->status.firstMeshMilliseconds = elapsed;
	}
	if (inputs.exactPresentation && !this->baselineCaptured)
	    this->maximumStructuralProxyCount = std::max(
		this->maximumStructuralProxyCount,
		inputs.structuralProxyCount);
	if (inputs.exactPresentation && inputs.populationComplete &&
	    this->maximumStructuralProxyCount > 0 && !this->baselineCaptured) {
	    this->baselineCaptured = true;
	    this->status.structuralProxyBaselineCount =
		this->maximumStructuralProxyCount;
	}
	if (inputs.exactPresentation && this->baselineCaptured &&
	    !this->status.halfStructuralProxiesReplaced &&
	    inputs.structuralProxyCount <=
		this->status.structuralProxyBaselineCount / 2) {
	    this->status.halfStructuralProxiesReplaced = TRUE;
	    this->status.halfStructuralProxyReplacementMilliseconds = elapsed;
	}
	if (inputs.stableView && !this->status.stableViewReached) {
	    this->status.stableViewReached = TRUE;
	    this->status.stableViewMilliseconds = elapsed;
	}
	return this->status;
    }

    void reset(void)
    {
	this->active = false;
	this->episodeRevision = 0;
	this->episodeStartMicroseconds = 0;
	this->baselineCaptured = false;
	this->maximumStructuralProxyCount = 0;
	this->status = BObolLodEpisodeStatus();
    }

private:
    void begin(const Inputs &inputs)
    {
	this->active = true;
	this->episodeRevision = inputs.episodeRevision;
	this->episodeStartMicroseconds =
	    inputs.episodeStartMicroseconds > 0 &&
	    inputs.episodeStartMicroseconds <= inputs.observationMicroseconds ?
	    inputs.episodeStartMicroseconds : inputs.observationMicroseconds;
	this->baselineCaptured = false;
	this->maximumStructuralProxyCount = 0;
	this->status = BObolLodEpisodeStatus();
    }

    uint64_t elapsedMilliseconds(int64_t observationMicroseconds) const
    {
	static constexpr uint64_t microsecondsPerMillisecond = 1000;
	if (observationMicroseconds <= this->episodeStartMicroseconds)
	    return 0;
	return static_cast<uint64_t>(observationMicroseconds -
	    this->episodeStartMicroseconds) / microsecondsPerMillisecond;
    }

    bool active = false;
    uint64_t episodeRevision = 0;
    int64_t episodeStartMicroseconds = 0;
    bool baselineCaptured = false;
    size_t maximumStructuralProxyCount = 0;
    BObolLodEpisodeStatus status;
};

#endif /* LIBBOBOL_LOD_PROGRESS_DIAGNOSTICS_PRIVATE_H */
