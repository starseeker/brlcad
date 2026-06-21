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

#include <cmath>

QgObolSnapRecord::QgObolSnapRecord(void) :
    point(0.0f, 0.0f, 0.0f),
    path(),
    editIntentId(),
    editIntentRole(),
    kind(NONE),
    primitiveIndex(-1),
    vertexIndex(-1),
    edgeSlot(-1),
    edgeVertexIndexA(-1),
    edgeVertexIndexB(-1),
    submittedSourceRequestCount(0),
    distance(0.0f),
    sourceFullDetailPending(false)
{
}

static int
qg_obol_snap_kind_priority(int kind)
{
    switch (kind) {
	case QgObolSnapRecord::VERTEX:
	    return 0;
	case QgObolSnapRecord::ENDPOINT:
	    return 1;
	case QgObolSnapRecord::MIDPOINT:
	    return 2;
	case QgObolSnapRecord::EDGE_NEAREST:
	    return 3;
	case QgObolSnapRecord::FACE_NEAREST:
	    return 4;
	case QgObolSnapRecord::CENTER:
	    return 5;
	case QgObolSnapRecord::LINE_NEAREST:
	    return 6;
	case QgObolSnapRecord::CONSTRUCTION_PLANE:
	    return 7;
	case QgObolSnapRecord::NONE:
	default:
	    return 8;
    }
}

static int
qg_obol_snap_record_is_better(const QgObolSnapRecord &candidate,
	const QgObolSnapRecord &current)
{
    const float tieTolerance = 1.0e-6f;

    if (candidate.kind == QgObolSnapRecord::NONE)
	return 0;
    if (current.kind == QgObolSnapRecord::NONE)
	return 1;
    if (candidate.distance < current.distance - tieTolerance)
	return 1;
    if (candidate.distance > current.distance + tieTolerance)
	return 0;

    return qg_obol_snap_kind_priority(candidate.kind) <
	qg_obol_snap_kind_priority(current.kind);
}

static void
qg_obol_snap_record_from_action(const SoBRLSnapAction &snapAction,
	QgObolSnapRecord &record)
{
    record.point = snapAction.getPoint();
    record.path = snapAction.getPath().getString();
    record.editIntentId = snapAction.getEditIntentId().getString();
    record.editIntentRole = snapAction.getEditIntentRole().getString();
    record.kind = static_cast<int>(snapAction.getKind());
    record.primitiveIndex = snapAction.getPrimitiveIndex();
    record.vertexIndex = snapAction.getVertexIndex();
    record.edgeSlot = snapAction.getEdgeSlot();
    record.edgeVertexIndexA = snapAction.getEdgeVertexIndexA();
    record.edgeVertexIndexB = snapAction.getEdgeVertexIndexB();
    record.distance = snapAction.getDistance();
}

static int
qg_obol_snap_consume_source_full_detail(BRLObolViewController *controller,
	SoBRLSnapAction &snapAction)
{
    if (!controller)
	return 0;

    int submitted = 0;
    (void)controller->consumeSnapSourceFullDetail(snapAction, 0,
	    &submitted);
    return submitted;
}

static int
qg_obol_snap_point_with_policy(QgView *display,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
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
    snapAction.setGeometryPolicy(geometryPolicy);
    snapAction.apply(controller->getViewport()->getRoot());

    if (geometryPolicy == SoBRLSnapAction::FULL_DETAIL) {
	record.submittedSourceRequestCount =
	    qg_obol_snap_consume_source_full_detail(controller, snapAction);
	record.sourceFullDetailPending =
	    record.submittedSourceRequestCount > 0;
    }

    if (!snapAction.hasCandidate())
	return 0;

    QgObolSnapRecord candidate;
    qg_obol_snap_record_from_action(snapAction, candidate);
    if (qg_obol_snap_record_is_better(candidate, record))
	record = candidate;
    return 1;
}

int
qg_obol_snap_point(QgView *display,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	QgObolSnapRecord &record)
{
    return qg_obol_snap_point_with_policy(display, query, tolerance,
	    enabledKinds, SoBRLSnapAction::DISPLAY_LEVEL, record);
}

int
qg_obol_snap_point_full_detail(QgView *display,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	QgObolSnapRecord &record)
{
    return qg_obol_snap_point_with_policy(display, query, tolerance,
	    enabledKinds, SoBRLSnapAction::FULL_DETAIL, record);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
