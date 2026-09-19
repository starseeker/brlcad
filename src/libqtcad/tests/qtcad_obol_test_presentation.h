/*      Q T C A D _ O B O L _ T E S T _ P R E S E N T A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#ifndef QTCAD_OBOL_TEST_PRESENTATION_H
#define QTCAD_OBOL_TEST_PRESENTATION_H

#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "qtcad/QgView.h"

#include <QImage>

/* Hidden Qt test views do not receive a normal expose/paint stream.  Present
 * one requested frame with the same transaction ordering as QgSW: claim the
 * request before traversal so a successor published during completion remains
 * pending, and report only the work actually observed by this traversal. */
static inline bool
qtcad_obol_present_requested_frame(QgView &view,
	BObolViewController *controller, QImage &image)
{
    image = QImage();
    if (!controller)
	return false;

    SbBool capacityRelevant = FALSE;
    SbBool planningRelevant = FALSE;
    const SbBool claimed = controller->consumeRenderRequest(NULL,
	&capacityRelevant, &planningRelevant);
    if (!claimed) {
	capacityRelevant = FALSE;
	planningRelevant = FALSE;
    } else if (!controller->isLodPresentationCapacityRelevant()) {
	capacityRelevant = FALSE;
    }

    BObolViewLodState *presentationState = controller->getViewLodState();
    const uint64_t executionBefore = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    if (presentationState)
	presentationState->beginCadPresentationFrame();

    const uint64_t started = controller->beginRenderTiming();
    view.get_viewport_image(image);
    const uint64_t completed = controller->beginRenderTiming();

    const uint64_t executionAfter = presentationState ?
	presentationState->cadPresentationExecutionSerial() : 0;
    if (presentationState)
	presentationState->refreshCadPresentationFrameStatus();
    const BObolCadPreparationProgress preparation = presentationState ?
	presentationState->cadPresentationPreparationProgress() :
	BOBOL_CAD_PREPARATION_NONE;
    const SbBool exact = !presentationState ||
	presentationState->lastCadPresentationFrameExact();
    const BObolPresentationTimingContext timingContext(
	capacityRelevant ? BObolLodCapacityRelevance::RELEVANT :
	    BObolLodCapacityRelevance::EXCLUDED,
	planningRelevant ? BObolLodPlanningRelevance::RELEVANT :
	    BObolLodPlanningRelevance::EXCLUDED,
	executionAfter != executionBefore ?
	    BObolCadPresentationExecution::EXECUTED :
	    BObolCadPresentationExecution::NOT_EXECUTED,
	preparation,
	exact ? BObolCadPresentationCompleteness::EXACT :
	    BObolCadPresentationCompleteness::INCOMPLETE);

    if (image.isNull()) {
	/* The request was claimed, but no framebuffer was produced.  Preserve a
	 * level-triggered presentation obligation instead of certifying work that
	 * was never shown. */
	controller->requestPresentationRender(
	    "qtcad-test-presentation-readback-failed");
	return false;
    }
    if (!exact) {
	controller->notePresentationRenderInterrupted(
	    completed > started ? completed - started : 1, timingContext);
	return false;
    }

    controller->completeRenderTiming(started, timingContext);
    return true;
}

#endif /* QTCAD_OBOL_TEST_PRESENTATION_H */
