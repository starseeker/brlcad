/*        L O D _ T E R M I N A L _ Q U A L I T Y _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_LOD_TERMINAL_QUALITY_PRIVATE_H
#define LIBBOBOL_LOD_TERMINAL_QUALITY_PRIVATE_H

#include "common.h"

#include <cstddef>
#include <cstdint>
#include <type_traits>

/*
 * Select the sole owner at the quiet terminal-quality boundary.
 *
 * A renderer-wide progressive ceiling hides occurrence-local triangle cuts.
 * Point quality changes which occurrences are independently presented.  A
 * point measurement made while that ceiling is effective therefore cannot
 * describe the requested population.  Reconciliation owns the boundary
 * first, including when the point frame belongs to an enclosing static trial.
 * A static triangle-only trial retains priority when no point presentation is
 * waiting for the ceiling-free population.
 *
 * This is the executable refinement of ObolTerminalQualityOrdering.tla.  It
 * is deliberately allocation-free and independent of controller state.
 */
class BObolLodTerminalQualityReducer {
public:
    enum class Owner : uint8_t {
	NONE = 0,
	STATIC_QUALITY,
	CEILING_RECONCILIATION,
	POINT_QUALITY
    };

    struct Inputs {
	bool stableTerminalContext = false;
	bool exactCompletedFrame = false;
	int progressiveCeiling = -1;
	size_t progressivePayloadCount = 0;
	bool staticPhaseActive = false;
	bool staticPhaseRejected = false;
	bool pointPresentationRequired = false;
    };

    struct Decision {
	Owner owner = Owner::NONE;
	bool globalCeilingDebt = false;
    };

    static bool globalCeilingBlocksPointPresentation(int progressiveCeiling,
	size_t progressivePayloadCount)
    {
	/* One progressive occurrence already has an occurrence-local ordinal.
	 * For every other population, including a stale zero-payload ceiling,
	 * the global control must be retired or translated before point quality
	 * can claim an exact presentation. */
	return progressiveCeiling >= 0 && progressivePayloadCount != 1;
    }

    static Decision reduce(const Inputs &inputs)
    {
	Decision decision;
	decision.globalCeilingDebt =
	    globalCeilingBlocksPointPresentation(
		inputs.progressiveCeiling,
		inputs.progressivePayloadCount);
	if (!inputs.stableTerminalContext || !inputs.exactCompletedFrame)
	    return decision;
	if (decision.globalCeilingDebt) {
	    /* The multi-occurrence ceiling is stronger than both point quality
	     * and an enclosing static trial.  A one-object static ordinal is not
	     * classified as global debt by the predicate above. */
	    decision.owner = Owner::CEILING_RECONCILIATION;
	    return decision;
	}

	if (inputs.pointPresentationRequired)
	    decision.owner = Owner::POINT_QUALITY;
	else if (inputs.staticPhaseActive && !inputs.staticPhaseRejected)
	    decision.owner = Owner::STATIC_QUALITY;
	return decision;
    }
};

static_assert(std::is_trivially_copyable<
	BObolLodTerminalQualityReducer::Inputs>::value,
    "terminal-quality inputs must remain allocation-free values");
static_assert(std::is_trivially_copyable<
	BObolLodTerminalQualityReducer::Decision>::value,
    "terminal-quality decisions must remain allocation-free values");

#endif /* LIBBOBOL_LOD_TERMINAL_QUALITY_PRIVATE_H */
