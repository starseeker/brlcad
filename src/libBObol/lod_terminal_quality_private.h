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
	POINT_QUALITY,
	ALLOCATION_PLANNING,
	ALLOCATION_CAPACITY
    };

    struct Inputs {
	bool stableTerminalContext = false;
	bool exactCompletedFrame = false;
	int progressiveCeiling = -1;
	size_t progressivePayloadCount = 0;
	bool staticPhaseActive = false;
	bool staticPhaseRejected = false;
	/* A completed constrained framebuffer may be the only representation of
	 * a multi-occurrence cut whose local minimum exceeds its hard allowance.
	 * That revision-bound proof makes the renderer ceiling terminal policy,
	 * rather than another reconciliation request. */
	bool retainedCeilingConstraint = false;
	/* A retained constraint is no longer terminal once the current
	 * occurrence-local allocation can encode a complete replacement within
	 * its certified allowance.  This can become true after point aggregation
	 * without changing the camera or geometry revision. */
	bool retainedAllocationReplacesCeiling = false;
	/* Structural fallback repair owns the occurrence population before any
	 * point, allocation, capacity, or static-quality successor.  A retained
	 * quality observation produced by the preceding source pass is subordinate
	 * evidence, not permission to allocate the incomplete population again. */
	bool structuralFrontierPending = false;
	bool pointPresentationRequired = false;
	bool allocationPlanningRequired = false;
	bool allocationCapacityRequired = false;
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

    static bool allocationDemandNeedsCapacity(
	size_t progressivePayloadCount, bool allocationCurrent,
	bool allocationPresentationRealized, size_t selectedPresentationCost,
	size_t pixelDemandPresentationCost, size_t certifiedPresentationBudget,
	bool terminalConstraintRecorded)
    {
	/* A single progressive occurrence uses the renderer's bounded ordinal
	 * staircase.  With several occurrences, only the scene allocator can
	 * decide which richer cut has the greatest screen value.  An exact
	 * presented allocation above its prior certificate is recertification
	 * work, not a reason to abandon ownership: structural publication may
	 * legitimately make the certificate stale before this reducer runs. */
	const bool presentedAllocationNeedsCertification =
	    selectedPresentationCost > certifiedPresentationBudget;
	const bool pixelDemandExceedsSelection =
	    pixelDemandPresentationCost > selectedPresentationCost;
	return progressivePayloadCount > 1 && allocationCurrent &&
	    allocationPresentationRealized && selectedPresentationCost > 0 &&
	    (presentedAllocationNeedsCertification ||
	     pixelDemandExceedsSelection) &&
	    !terminalConstraintRecorded;
    }

    static bool allocationDemandNeedsPlanning(
	size_t progressivePayloadCount, bool allocationCurrent,
	bool allocationPresentationRealized, bool qualityDebt,
	bool terminalConstraintRecorded)
    {
	/* Point refinement changes the occurrence population and invalidates the
	 * preceding allocation certificate.  Before capacity can certify or
	 * constrain that new population, one complete occurrence-local allocation
	 * must make its cuts current and presentable. */
	return progressivePayloadCount > 1 && qualityDebt &&
	    (!allocationCurrent || !allocationPresentationRealized) &&
	    !terminalConstraintRecorded;
    }

    static bool allocationConstraintPresented(bool constraintRecorded,
	bool allocationCurrent, bool allocationPresentationRealized)
    {
	/* Capacity evidence belongs to the exact occurrence population which it
	 * measured.  A point-threshold transition may retain the numeric budget
	 * while invalidating that population; the old certificate cannot suppress
	 * planning for its successor. */
	return constraintRecorded && allocationCurrent &&
	    allocationPresentationRealized;
    }

    static Decision reduce(const Inputs &inputs)
    {
	Decision decision;
	const bool globalCeiling = globalCeilingBlocksPointPresentation(
	    inputs.progressiveCeiling, inputs.progressivePayloadCount);
	const bool effectiveRetainedConstraint =
	    inputs.retainedCeilingConstraint &&
	    !inputs.retainedAllocationReplacesCeiling;
	decision.globalCeilingDebt = globalCeiling &&
	    !effectiveRetainedConstraint;
	if (!inputs.stableTerminalContext || !inputs.exactCompletedFrame)
	    return decision;
	if (inputs.structuralFrontierPending)
	    return decision;
	/* The constraint is a terminal completed-frame proof.  Point and capacity
	 * demand hidden behind it are deliberately constrained, not actionable
	 * successors which may reopen the impossible reconciliation. */
	if (globalCeiling && effectiveRetainedConstraint)
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
	else if (inputs.allocationPlanningRequired)
	    decision.owner = Owner::ALLOCATION_PLANNING;
	else if (inputs.allocationCapacityRequired)
	    decision.owner = Owner::ALLOCATION_CAPACITY;
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
