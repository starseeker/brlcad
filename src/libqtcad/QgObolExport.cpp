/*                Q G O B O L E X P O R T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolExport.cpp */

#include "common.h"

#include "qtcad/QgObolExport.h"

#include "brlobol/export_action.h"
#include "brlobol/view_controller.h"
#include "qtcad/QgView.h"

#include <Inventor/SoViewport.h>

#include <vector>

QgObolExportTriangleRecord::QgObolExportTriangleRecord(void) :
    path(),
    sourceName(),
    sourceType(),
    editIntentId(),
    editIntentRole(),
    sourceId(0),
    primitiveIndex(-1),
    vertexIndexA(-1),
    vertexIndexB(-1),
    vertexIndexC(-1),
    a(0.0f, 0.0f, 0.0f),
    b(0.0f, 0.0f, 0.0f),
    c(0.0f, 0.0f, 0.0f)
{
}

QgObolExportGeometryRecord::QgObolExportGeometryRecord(void) :
    lineCount(0),
    pointCount(0),
    triangleCount(0),
    submittedSourceRequestCount(0),
    boundsMin(0.0f, 0.0f, 0.0f),
    boundsMax(0.0f, 0.0f, 0.0f),
    boundsValid(false),
    sourceFullDetailPending(false),
    triangles()
{
}

static int
qg_obol_export_consume_source_full_detail(BRLObolViewController *controller,
	SoBRLExportAction &exportAction)
{
    if (!controller)
	return 0;

    int submitted = 0;
    (void)controller->consumeExportSourceFullDetail(exportAction, 0,
	    &submitted);
    return submitted;
}

static void
qg_obol_export_record_from_action(const SoBRLExportAction &exportAction,
	QgObolExportGeometryRecord &record)
{
    record.lineCount = exportAction.getLineCount();
    record.pointCount = exportAction.getPointCount();
    record.triangleCount = exportAction.getTriangleCount();

    const SbBox3f &bounds = exportAction.getBounds();
    if (!bounds.isEmpty()) {
	record.boundsValid = true;
	record.boundsMin = bounds.getMin();
	record.boundsMax = bounds.getMax();
    }

    record.triangles.clear();
    record.triangles.reserve(static_cast<size_t>(record.triangleCount));
    for (int i = 0; i < exportAction.getTriangleCount(); i++) {
	const SoBRLExportAction::TriangleRecord &triangle =
	    exportAction.getTriangle(i);
	QgObolExportTriangleRecord out;
	out.path = triangle.path.getString();
	out.sourceName = triangle.sourceName.getString();
	out.sourceType = triangle.sourceType.getString();
	out.editIntentId = triangle.editIntentId.getString();
	out.editIntentRole = triangle.editIntentRole.getString();
	out.sourceId = triangle.sourceId;
	out.primitiveIndex = triangle.primitiveIndex;
	out.vertexIndexA = triangle.vertexIndexA;
	out.vertexIndexB = triangle.vertexIndexB;
	out.vertexIndexC = triangle.vertexIndexC;
	out.a = triangle.a;
	out.b = triangle.b;
	out.c = triangle.c;
	record.triangles.push_back(out);
    }
}

int
qg_obol_export_geometry_full_detail(QgView *display,
	QgObolExportGeometryRecord &record)
{
    record = QgObolExportGeometryRecord();
    if (!display)
	return 0;

    BRLObolViewController *controller = display->obolViewController();
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;

    SoBRLExportAction exportAction;
    exportAction.setGeometryPolicy(SoBRLExportAction::FULL_DETAIL);
    exportAction.apply(controller->getViewport()->getRoot());
    record.submittedSourceRequestCount =
	qg_obol_export_consume_source_full_detail(controller, exportAction);
    record.sourceFullDetailPending =
	record.submittedSourceRequestCount > 0;

    if (exportAction.getLineCount() <= 0 &&
	    exportAction.getPointCount() <= 0 &&
	    exportAction.getTriangleCount() <= 0)
	return 0;

    qg_obol_export_record_from_action(exportAction, record);
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
