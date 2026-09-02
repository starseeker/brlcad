/*       L O D _ Q U I E T _ S U C C E S S O R _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_QUIET_SUCCESSOR_PRIVATE_H
#define LIBBOBOL_LOD_QUIET_SUCCESSOR_PRIVATE_H

#include "common.h"

#include <algorithm>
#include <cmath>
#include <cstddef>
#include <cstdint>
#include <type_traits>

/*
 * Pure selection of the first quiet presentation after camera interaction.
 *
 * Motion limits are transient renderer safety controls.  Their publication
 * timing must not choose the semantic quiet target.  A revision-matched
 * certificate does choose that target; otherwise quiet drawing returns to
 * the one-pixel demand while a bounded retained-allocation handoff translates
 * any motion-only ceiling into occurrence-local cuts.
 *
 * This value deliberately knows nothing about cameras, scenes, clocks,
 * renderers, or workload cardinality.  The presentation policy validates
 * certificate identity before constructing these inputs.
 */
class BObolLodQuietSuccessorReducer {
public:
    enum class Source : uint8_t {
	STABLE_DEMAND = 0,
	PRIOR_STABLE,
	PROVEN_SCALE,
	EXACT_VIEW
    };

    enum class Handoff : uint8_t {
	NONE = 0,
	ALLOCATION,
	PRESENTATION
    };

    struct Target {
	bool valid = false;
	float targetPixelError = 1.0f;
	int progressiveCeiling = -1;
	float progressiveNextFraction = 0.0f;
	float pointProxyPixelThreshold = 1.0f;
	size_t presentedRenderCost = 0;
	size_t provenRenderCostCapacity = 0;
    };

    struct Inputs {
	bool retainedMeshPayloads = false;
	/* Prior stable controls remain a valid starting point after an
	 * orthographic zoom, but only a pose-only change may reuse their exact
	 * occurrence cuts.  Keep those facts separate: conflating them made a
	 * zoom inherit the coarse motion-time point threshold and then visibly
	 * rebuild the stable population. */
	bool priorStableCutsReusable = false;
	Target current;
	Target priorStable;
	Target provenScale;
	Target exactView;
    };

    struct Decision {
	Source source = Source::STABLE_DEMAND;
	Handoff handoff = Handoff::NONE;
	Target target;
	bool apply = false;
	bool clearPresentationLimits = false;

	bool needsHandoff(void) const
	{
	    return this->handoff != Handoff::NONE;
	}
    };

    static Decision reduce(const Inputs &inputs)
    {
	const Target current = sanitize(inputs.current, true);
	Decision decision;
	decision.target = stableDemand();

	if (inputs.exactView.valid) {
	    decision.source = Source::EXACT_VIEW;
	    decision.target = sanitize(inputs.exactView, false);
	} else if (!inputs.retainedMeshPayloads) {
	    decision.clearPresentationLimits = true;
	} else if (inputs.provenScale.valid) {
	    decision.source = Source::PROVEN_SCALE;
	    decision.target = sanitize(inputs.provenScale, false);
	    if (inputs.priorStable.valid) {
		const Target stableControls = sanitize(
		    inputs.priorStable, false);
		/* The current-scale frame proves a drawable presentation and its
		 * cost, not a new user-facing fidelity target.  Preserve its
		 * renderer ceiling while recomputing occurrence demand from the
		 * stable pixel and proxy controls captured before motion. */
		decision.target.targetPixelError =
		    stableControls.targetPixelError;
		decision.target.pointProxyPixelThreshold =
		    stableControls.pointProxyPixelThreshold;
	    }
	    decision.handoff = constrained(decision.target) ?
		Handoff::ALLOCATION : Handoff::NONE;
	} else if (inputs.priorStable.valid) {
	    decision.source = Source::PRIOR_STABLE;
	    decision.target = sanitize(inputs.priorStable, false);
	    /* Pose-only orthographic motion preserves immutable occurrence cuts.
	     * Scale changes preserve only the stable control vector; their new
	     * projected demand still needs one atomic retained allocation. */
	    decision.handoff = inputs.priorStableCutsReusable ?
		Handoff::PRESENTATION : Handoff::ALLOCATION;
	} else {
	    /* Preserve only the renderer safety controls while the allocator
	     * translates them.  The semantic pixel demand is stable and therefore
	     * cannot depend on which intermediate motion frame arrived last. */
	    decision.target.progressiveCeiling = current.progressiveCeiling;
	    decision.target.progressiveNextFraction =
		current.progressiveNextFraction;
	    decision.target.pointProxyPixelThreshold =
		current.pointProxyPixelThreshold;
	    decision.handoff = constrained(current) ?
		Handoff::ALLOCATION : Handoff::NONE;
	}

	decision.apply = !equivalent(current, decision.target);
	return decision;
    }

private:
    static constexpr float MIN_PIXEL_ERROR = 0.25f;
    static constexpr float MAX_PIXEL_ERROR = 64.0f;
    static constexpr float DEFAULT_PIXEL_ERROR = 1.0f;
    static constexpr float MIN_PROXY_PIXEL_THRESHOLD = 1.0f;
    static constexpr float MAX_PROXY_PIXEL_THRESHOLD = 64.0f;
    static constexpr float CONSTRAINED_PROXY_THRESHOLD = 1.01f;
    static constexpr float FLOAT_COMPARISON_EPSILON = 1.0e-6f;

    static Target stableDemand(void)
    {
	Target target;
	target.valid = true;
	return target;
    }

    static Target sanitize(const Target &input, bool allowInvalid)
    {
	Target target = input;
	if (!target.valid && allowInvalid)
	    target = stableDemand();
	target.valid = true;
	target.targetPixelError = std::isfinite(target.targetPixelError) ?
	    std::max(MIN_PIXEL_ERROR,
		std::min(MAX_PIXEL_ERROR, target.targetPixelError)) :
	    DEFAULT_PIXEL_ERROR;
	target.progressiveCeiling = target.progressiveCeiling < -1 ? -1 :
	    target.progressiveCeiling;
	if (target.progressiveCeiling < 0 ||
	    !std::isfinite(target.progressiveNextFraction))
	    target.progressiveNextFraction = 0.0f;
	else
	    target.progressiveNextFraction = std::max(0.0f,
		std::min(1.0f, target.progressiveNextFraction));
	target.pointProxyPixelThreshold =
	    std::isfinite(target.pointProxyPixelThreshold) ?
		std::max(MIN_PROXY_PIXEL_THRESHOLD,
		    std::min(MAX_PROXY_PIXEL_THRESHOLD,
			target.pointProxyPixelThreshold)) :
		MIN_PROXY_PIXEL_THRESHOLD;
	return target;
    }

    static bool constrained(const Target &target)
    {
	return target.progressiveCeiling >= 0 ||
	    target.pointProxyPixelThreshold > CONSTRAINED_PROXY_THRESHOLD;
    }

    static bool equivalent(const Target &left, const Target &right)
    {
	return std::fabs(left.targetPixelError - right.targetPixelError) <=
		FLOAT_COMPARISON_EPSILON &&
	    left.progressiveCeiling == right.progressiveCeiling &&
	    std::fabs(left.progressiveNextFraction -
		right.progressiveNextFraction) <= FLOAT_COMPARISON_EPSILON &&
	    std::fabs(left.pointProxyPixelThreshold -
		right.pointProxyPixelThreshold) <= FLOAT_COMPARISON_EPSILON;
    }
};

static_assert(std::is_trivially_copyable<
	BObolLodQuietSuccessorReducer::Target>::value,
    "quiet presentation targets must remain allocation-free values");
static_assert(std::is_trivially_copyable<
	BObolLodQuietSuccessorReducer::Decision>::value,
    "quiet successor decisions must remain allocation-free values");

#endif /* LIBBOBOL_LOD_QUIET_SUCCESSOR_PRIVATE_H */
