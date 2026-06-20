/*                Q G O B O L S N A P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolSnap.cpp */

#include "common.h"

#include "qtcad/QgObolSnap.h"

#include "brlobol/snap_action.h"
#include "brlobol/view_controller.h"
#include "qtcad/QgView.h"

#include <Inventor/SoViewport.h>

QgObolSnapRecord::QgObolSnapRecord(void) :
    point(0.0f, 0.0f, 0.0f),
    path(),
    kind(NONE),
    primitiveIndex(-1),
    distance(0.0f)
{
}

int
qg_obol_snap_point(QgView *display,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	QgObolSnapRecord &record)
{
    record = QgObolSnapRecord();
    if (!display)
	return 0;

    BRLObolViewController *controller = display->obolViewController();
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;

    SoBRLSnapAction snapAction;
    snapAction.setQueryPoint(query);
    snapAction.setTolerance(tolerance > 0.0f ? tolerance : 1.0f);
    snapAction.setEnabledKinds(enabledKinds ? enabledKinds :
	    static_cast<uint32_t>(QgObolSnapRecord::ALL_KINDS));
    snapAction.setPriorityPolicy(SoBRLSnapAction::FEATURE_PRIORITY);
    snapAction.apply(controller->getViewport()->getRoot());

    if (!snapAction.hasCandidate())
	return 0;

    record.point = snapAction.getPoint();
    record.path = snapAction.getPath().getString();
    record.kind = static_cast<int>(snapAction.getKind());
    record.primitiveIndex = snapAction.getPrimitiveIndex();
    record.distance = snapAction.getDistance();
    return 1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
