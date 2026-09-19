/* P R E S E N T A T I O N _ P R E P A R A T I O N _ P R I V A T E . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef LIBBOBOL_PRESENTATION_PREPARATION_PRIVATE_H
#define LIBBOBOL_PRESENTATION_PREPARATION_PRIVATE_H

#include "BObol/BPresentationPreparation.h"

#include <Obol/cad/CadPresentationPreparation.h>

class BObolCadPreparationPolicy {
public:
    static BObolCadPreparationProgress classify(
	const Obol::CadPresentationPreparationSnapshot &before,
	const Obol::CadPresentationPreparationSnapshot &after)
    {
	using State = Obol::CadPresentationPreparationState;
	if (!after.hasTarget())
	    return BOBOL_CAD_PREPARATION_NONE;
	const bool sameTarget = before.hasTarget() &&
	    before.target == after.target;
	/* Validate the finite-work certificate before interpreting its terminal
	 * state.  A newly published or constrained target is not permission to
	 * bypass the same bounds and monotonicity contract. */
	if (after.completedUnits > after.totalUnits ||
		(after.state == State::Complete &&
		 after.completedUnits != after.totalUnits) ||
		(sameTarget &&
		 (before.totalUnits != after.totalUnits ||
		  after.completedUnits < before.completedUnits)))
	    return BOBOL_CAD_PREPARATION_FAILED;
	if (after.state == State::Failed)
	    return BOBOL_CAD_PREPARATION_FAILED;
	if (after.state == State::Constrained)
	    return BOBOL_CAD_PREPARATION_CONSTRAINED;

	if (!sameTarget)
	    return after.state == State::Complete ?
		BOBOL_CAD_PREPARATION_COMPLETED :
		BOBOL_CAD_PREPARATION_STARTED;
	if (after.state == State::Complete &&
		before.state != State::Complete)
	    return BOBOL_CAD_PREPARATION_COMPLETED;
	if (after.completedUnits > before.completedUnits &&
		after.remainingUnits() < before.remainingUnits())
	    return BOBOL_CAD_PREPARATION_ADVANCED;
	return BOBOL_CAD_PREPARATION_NONE;
    }

    static BObolCadPreparationProgress combine(
	BObolCadPreparationProgress left,
	BObolCadPreparationProgress right)
    {
	return rank(right) > rank(left) ? right : left;
    }

    static bool permitsUnchangedRetry(BObolCadPreparationProgress progress)
    {
	return progress == BOBOL_CAD_PREPARATION_STARTED ||
	    progress == BOBOL_CAD_PREPARATION_ADVANCED ||
	    progress == BOBOL_CAD_PREPARATION_COMPLETED;
    }

private:
    static int rank(BObolCadPreparationProgress progress)
    {
	switch (progress) {
	    case BOBOL_CAD_PREPARATION_FAILED: return 5;
	    case BOBOL_CAD_PREPARATION_CONSTRAINED: return 4;
	    case BOBOL_CAD_PREPARATION_ADVANCED: return 3;
	    case BOBOL_CAD_PREPARATION_COMPLETED: return 2;
	    case BOBOL_CAD_PREPARATION_STARTED: return 1;
	    case BOBOL_CAD_PREPARATION_NONE: return 0;
	}
	return 5;
    }
};

#endif /* LIBBOBOL_PRESENTATION_PREPARATION_PRIVATE_H */
