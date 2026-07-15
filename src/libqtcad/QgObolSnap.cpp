/*                Q G O B O L S N A P . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolSnap.cpp */

#include "common.h"

#include "qtcad/QgObolSnap.h"

#include "brlobol/view_controller.h"
#include "brlobol/view_query.h"
#include "qtcad/QgView.h"

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

static void
qg_obol_snap_record(const BRLObolViewSnapRecord &source,
	QgObolSnapRecord &record)
{
    record.point = source.point;
    record.path = source.path.getString();
    record.editIntentId = source.editIntentId.getString();
    record.editIntentRole = source.editIntentRole.getString();
    record.kind = static_cast<int>(source.kind);
    record.primitiveIndex = source.primitiveIndex;
    record.vertexIndex = source.vertexIndex;
    record.edgeSlot = source.edgeSlot;
    record.edgeVertexIndexA = source.edgeVertexIndexA;
    record.edgeVertexIndexB = source.edgeVertexIndexB;
    record.submittedSourceRequestCount = source.submittedSourceRequestCount;
    record.distance = source.distance;
    record.sourceFullDetailPending = source.sourceFullDetailPending ? true : false;
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
    BRLObolViewSnapRecord snap;
    const int found = brlobol_view_snap_point(controller, query, tolerance,
	enabledKinds, geometryPolicy,
	geometryPolicy == SoBRLSnapAction::FULL_DETAIL, snap);
    qg_obol_snap_record(snap, record);
    return found;
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
