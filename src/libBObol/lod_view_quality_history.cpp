/*        L O D _ V I E W _ Q U A L I T Y _ H I S T O R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "lod_presentation_policy_private.h"

#include <algorithm>
#include <cmath>
#include <limits>

namespace {

constexpr float minimumTargetPixelError = 0.25f;
constexpr float maximumTargetPixelError = 1.0f;
constexpr float maximumTargetPixelErrorTolerance = 1.01f;
constexpr float minimumTargetPixelErrorTolerance = 0.249f;
constexpr float minimumPointProxyPixelThreshold = 1.0f;
constexpr float maximumPointProxyPixelThreshold = 64.0f;
constexpr float maximumPointProxyPixelThresholdTolerance = 64.01f;
constexpr float fractionTolerance = 1.0001f;
constexpr float controlComparisonTolerance = 1.0e-6f;
constexpr double rankComparisonTolerance = 1.0e-6;
constexpr double errorComparisonTolerance = 1.0e-9;

}

bool
BObolLodViewQualityHistory::remember(const RememberInputs &inputs)
{
    if (!inputs.view.haveCamera || !inputs.domainRevision ||
	!inputs.sceneAvailable || !inputs.exactCompletedFrame ||
	!inputs.terminalPresentationComplete || !inputs.producersSettled ||
	!inputs.presentationDeadlineMet || inputs.resourcePressure ||
	!validQuality(inputs.quality))
	return false;

    size_t index = this->entryCount;
    for (size_t i = 0; i < this->entryCount; ++i) {
	if (this->entries[i].domainRevision == inputs.domainRevision &&
	    this->entries[i].view.same(inputs.view)) {
	    index = i;
	    break;
	}
    }

    /* Throughput and fidelity are independent evidence.  Preserve the
     * largest safe workload while replacing visual controls only with a
     * demonstrably better candidate. */
    if (index < this->entryCount) {
	Entry &current = this->entries[index];
	current.provenRenderCostCapacity = std::max(
	    current.provenRenderCostCapacity,
	    inputs.quality.presentedRenderCost);
	if (!preferFidelityCandidate(inputs.quality, current.quality)) {
	    this->promote(index);
	    return false;
	}
    }

    Entry candidate;
    candidate.valid = true;
    candidate.view = inputs.view;
    candidate.domainRevision = inputs.domainRevision;
    candidate.quality = sanitizeQuality(inputs.quality);
    candidate.provenRenderCostCapacity = index < this->entryCount ?
	this->entries[index].provenRenderCostCapacity :
	candidate.quality.presentedRenderCost;

    if (index == this->entryCount) {
	if (this->entryCount < capacity())
	    this->entryCount++;
	index = this->entryCount - 1;
    }
    for (size_t i = index; i > 0; --i)
	this->entries[i] = this->entries[i - 1];
    this->entries[0] = candidate;
    if (this->rememberCountValue != SIZE_MAX)
	this->rememberCountValue++;
    return true;
}

bool
BObolLodViewQualityHistory::recall(const RecallInputs &inputs,
	Quality &quality)
{
    if (!inputs.view.haveCamera || !inputs.domainRevision ||
	!inputs.sceneAvailable || inputs.resourcePressure)
	return false;
    for (size_t i = 0; i < this->entryCount; ++i) {
	if (!this->entries[i].valid ||
	    this->entries[i].domainRevision != inputs.domainRevision ||
	    !this->entries[i].view.same(inputs.view))
	    continue;
	quality = this->entries[i].quality;
	quality.provenRenderCostCapacity =
	    this->entries[i].provenRenderCostCapacity;
	this->promote(i);
	if (this->recallCountValue != SIZE_MAX)
	    this->recallCountValue++;
	return true;
    }
    return false;
}

void
BObolLodViewQualityHistory::reset(void)
{
    this->entries = {};
    this->entryCount = 0;
    this->rememberCountValue = 0;
    this->recallCountValue = 0;
}

size_t
BObolLodViewQualityHistory::size(void) const
{
    return this->entryCount;
}

size_t
BObolLodViewQualityHistory::rememberCount(void) const
{
    return this->rememberCountValue;
}

size_t
BObolLodViewQualityHistory::recallCount(void) const
{
    return this->recallCountValue;
}

bool
BObolLodViewQualityHistory::validQuality(const Quality &quality)
{
    return std::isfinite(quality.targetPixelError) &&
	quality.targetPixelError >= minimumTargetPixelErrorTolerance &&
	quality.targetPixelError <= maximumTargetPixelErrorTolerance &&
	quality.progressiveCeiling >= -1 &&
	quality.progressiveCeiling < BOBOL_MESH_LOD_CUT_COUNT_MAX &&
	std::isfinite(quality.progressiveNextFraction) &&
	quality.progressiveNextFraction >= 0.0f &&
	quality.progressiveNextFraction <= fractionTolerance &&
	(quality.progressiveCeiling >= 0 ||
	 quality.progressiveNextFraction <= controlComparisonTolerance) &&
	std::isfinite(quality.pointProxyPixelThreshold) &&
	quality.pointProxyPixelThreshold >= minimumPointProxyPixelThreshold &&
	quality.pointProxyPixelThreshold <=
	    maximumPointProxyPixelThresholdTolerance &&
	(std::isfinite(quality.maximumProjectedErrorPixels) ?
	    quality.maximumProjectedErrorPixels >= 0.0 :
	    std::isinf(quality.maximumProjectedErrorPixels) &&
	    quality.maximumProjectedErrorPixels > 0.0) &&
	quality.presentedRenderCost > 0;
}

BObolLodViewQualityHistory::Quality
BObolLodViewQualityHistory::sanitizeQuality(const Quality &quality)
{
    Quality result = quality;
    result.targetPixelError = std::max(minimumTargetPixelError,
	std::min(maximumTargetPixelError, result.targetPixelError));
    result.progressiveCeiling = result.progressiveCeiling < -1 ? -1 :
	std::min<int>(BOBOL_MESH_LOD_CUT_COUNT_MAX - 1,
	    result.progressiveCeiling);
    result.progressiveNextFraction = result.progressiveCeiling < 0 ? 0.0f :
	std::max(0.0f, std::min(1.0f, result.progressiveNextFraction));
    result.pointProxyPixelThreshold = std::max(
	minimumPointProxyPixelThreshold,
	std::min(maximumPointProxyPixelThreshold,
	    result.pointProxyPixelThreshold));
    if (!std::isfinite(result.maximumProjectedErrorPixels) ||
	result.maximumProjectedErrorPixels < 0.0)
	result.maximumProjectedErrorPixels =
	    std::numeric_limits<double>::infinity();
    result.provenRenderCostCapacity = 0;
    return result;
}

double
BObolLodViewQualityHistory::progressiveRank(int ceiling, float nextFraction)
{
    return ceiling < 0 ?
	static_cast<double>(BOBOL_MESH_LOD_CUT_COUNT_MAX) :
	static_cast<double>(ceiling) + static_cast<double>(
	    std::max(0.0f, std::min(1.0f, nextFraction)));
}

bool
BObolLodViewQualityHistory::controlsDominate(const Quality &candidate,
	const Quality &current)
{
    const bool pixelNoWorse = candidate.targetPixelError <=
	current.targetPixelError + controlComparisonTolerance;
    const bool proxyNoWorse = candidate.pointProxyPixelThreshold <=
	current.pointProxyPixelThreshold + controlComparisonTolerance;
    const bool ceilingNoWorse = progressiveRank(
	candidate.progressiveCeiling,
	candidate.progressiveNextFraction) + rankComparisonTolerance >=
	progressiveRank(
	    current.progressiveCeiling, current.progressiveNextFraction);
    const bool strictlyBetter = candidate.targetPixelError <
	    current.targetPixelError - controlComparisonTolerance ||
	candidate.pointProxyPixelThreshold <
	    current.pointProxyPixelThreshold - controlComparisonTolerance ||
	progressiveRank(candidate.progressiveCeiling,
	    candidate.progressiveNextFraction) > progressiveRank(
		current.progressiveCeiling,
		current.progressiveNextFraction) + rankComparisonTolerance;
    return pixelNoWorse && proxyNoWorse && ceilingNoWorse && strictlyBetter;
}

bool
BObolLodViewQualityHistory::controlsEquivalent(const Quality &candidate,
	const Quality &current)
{
    return std::fabs(candidate.targetPixelError -
	    current.targetPixelError) <= controlComparisonTolerance &&
	std::fabs(candidate.pointProxyPixelThreshold -
	    current.pointProxyPixelThreshold) <= controlComparisonTolerance &&
	std::fabs(progressiveRank(candidate.progressiveCeiling,
	    candidate.progressiveNextFraction) - progressiveRank(
		current.progressiveCeiling,
		current.progressiveNextFraction)) <= rankComparisonTolerance;
}

bool
BObolLodViewQualityHistory::preferFidelityCandidate(
	const Quality &candidate, const Quality &current)
{
    const bool candidateHasBound =
	std::isfinite(candidate.maximumProjectedErrorPixels);
    const bool currentHasBound =
	std::isfinite(current.maximumProjectedErrorPixels);
    if (candidateHasBound != currentHasBound)
	return controlsEquivalent(candidate, current) ||
	    controlsDominate(candidate, current);
    if (candidateHasBound && std::fabs(
	    candidate.maximumProjectedErrorPixels -
	    current.maximumProjectedErrorPixels) > errorComparisonTolerance)
	return candidate.maximumProjectedErrorPixels <
	    current.maximumProjectedErrorPixels;
    return controlsDominate(candidate, current);
}

void
BObolLodViewQualityHistory::promote(size_t index)
{
    if (index == 0 || index >= this->entryCount)
	return;
    const Entry entry = this->entries[index];
    for (size_t i = index; i > 0; --i)
	this->entries[i] = this->entries[i - 1];
    this->entries[0] = entry;
}
