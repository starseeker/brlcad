/*                Q G O B O L E X P O R T . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolExport.cpp */

#include "common.h"

#include "qtcad/QgObolExport.h"

#include "BObol/BExportAction.h"
#include "BObol/BViewController.h"
#include "qtcad/QgView.h"

#include <Inventor/SoViewport.h>

#include <vector>

QgObolExportTriangleRecord::QgObolExportTriangleRecord(void) :
    path(),
    sourceName(),
    sourceType(),
    displayName(),
    geometryName(),
    cacheIdentity(),
    sourceIdentity(),
    cacheIdentityValue(0),
    sourceIdentityValue(0),
    editIntentId(),
    editIntentRole(),
    sourceId(0),
    databaseIntent(false),
    overlayIntent(false),
    hudIntent(false),
    localSource(false),
    sharedSource(false),
    nonDatabaseSource(false),
    drawMode(0),
    recordRole(),
    geometryKind(),
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

QgObolExportObjectLineSummary::QgObolExportObjectLineSummary(void) :
    valid(false),
    pointCount(0),
    segmentCount(0),
    cacheIdentity(),
    sourceIdentity(),
    cacheIdentityValue(0),
    sourceIdentityValue(0)
{
}

QgObolExportObjectPointSummary::QgObolExportObjectPointSummary(void) :
    valid(false),
    pointCount(0),
    cacheIdentity(),
    sourceIdentity(),
    cacheIdentityValue(0),
    sourceIdentityValue(0)
{
}

QgObolExportObjectSurfaceSummary::QgObolExportObjectSurfaceSummary(void) :
    valid(false),
    pointCount(0),
    normalCount(0),
    indexCount(0),
    faceCount(0),
    normalsPerIndex(0),
    materialColorValid(false),
    materialColor(0.0f, 0.0f, 0.0f),
    materialDrawMode(0),
    materialHighlighted(false),
    cacheIdentity(),
    sourceIdentity(),
    cacheIdentityValue(0),
    sourceIdentityValue(0)
{
}

QgObolExportObjectRecord::QgObolExportObjectRecord(void) :
    path(),
    sourceName(),
    sourceType(),
    displayName(),
    geometryName(),
    cacheIdentity(),
    sourceIdentity(),
    cacheIdentityValue(0),
    sourceIdentityValue(0),
    sourceId(0),
    databaseIntent(false),
    overlayIntent(false),
    hudIntent(false),
    localSource(false),
    sharedSource(false),
    nonDatabaseSource(false),
    drawMode(0),
    recordRole(),
    geometryKind(),
    selected(false),
    highlighted(false),
    visible(false),
    colorOverride(false),
    color(0.0f, 0.0f, 0.0f),
    boundsValid(false),
    boundsMin(0.0f, 0.0f, 0.0f),
    boundsMax(0.0f, 0.0f, 0.0f),
    lineCount(0),
    pointPrimitiveCount(0),
    triangleCount(0),
    lineSummary(),
    pointSummary(),
    surfaceSummary(),
    linePoints(),
    lineCommands(),
    points(),
    surfacePoints(),
    surfaceIndices()
{
}

QgObolExportObjectQuery::QgObolExportObjectQuery(void) :
    flags(0),
    glob(),
    drawMode(-1),
    geometryPolicy(QG_OBOL_EXPORT_FULL_DETAIL)
{
}

QgObolExportObjectQueryResult::QgObolExportObjectQueryResult(void) :
    lineCount(0),
    pointCount(0),
    triangleCount(0),
    submittedSourceRequestCount(0),
    sourceFullDetailPending(false),
    objects()
{
}

static int
qg_obol_export_consume_source_full_detail(BObolViewController *controller,
	SoBRLExportAction &exportAction)
{
    if (!controller)
	return 0;

    int submitted = 0;
    (void)controller->consumeExportSourceFullDetail(exportAction, 0,
	    &submitted);
    return submitted;
}

static SoBRLExportAction::GeometryPolicy
qg_obol_export_geometry_policy(int policy)
{
    return policy == QG_OBOL_EXPORT_DISPLAY_LEVEL ?
	SoBRLExportAction::DISPLAY_LEVEL : SoBRLExportAction::FULL_DETAIL;
}

static unsigned int
qg_obol_export_object_query_flags(unsigned int flags)
{
    unsigned int queryFlags = 0;
    if (flags & QG_OBOL_EXPORT_QUERY_VISIBLE_ONLY)
	queryFlags |= SoBRLExportAction::QUERY_VISIBLE_ONLY;
    if (flags & QG_OBOL_EXPORT_QUERY_DATABASE_OBJECTS)
	queryFlags |= SoBRLExportAction::QUERY_DATABASE_OBJECTS;
    if (flags & QG_OBOL_EXPORT_QUERY_VIEW_OBJECTS)
	queryFlags |= SoBRLExportAction::QUERY_VIEW_OBJECTS;
    if (flags & QG_OBOL_EXPORT_QUERY_LOCAL_ONLY)
	queryFlags |= SoBRLExportAction::QUERY_LOCAL_ONLY;
    return queryFlags;
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
	out.displayName = triangle.displayName.getString();
	out.geometryName = triangle.geometryName.getString();
	out.cacheIdentity = triangle.cacheIdentity.getString();
	out.sourceIdentity = triangle.sourceIdentity.getString();
	out.cacheIdentityValue = triangle.cacheIdentityValue;
	out.sourceIdentityValue = triangle.sourceIdentityValue;
	out.editIntentId = triangle.editIntentId.getString();
	out.editIntentRole = triangle.editIntentRole.getString();
	out.sourceId = triangle.sourceId;
	out.databaseIntent = triangle.databaseIntent != 0;
	out.overlayIntent = triangle.overlayIntent != 0;
	out.hudIntent = triangle.hudIntent != 0;
	out.localSource = triangle.localSource != 0;
	out.sharedSource = triangle.sharedSource != 0;
	out.nonDatabaseSource = triangle.nonDatabaseSource != 0;
	out.drawMode = triangle.drawMode;
	out.recordRole = triangle.recordRole.getString();
	out.geometryKind = triangle.geometryKind.getString();
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

static void
qg_obol_export_line_summary_from_action(
	const SoBRLExportAction::ObjectLineSummary &summary,
	QgObolExportObjectLineSummary &out)
{
    out.valid = summary.valid != 0;
    out.pointCount = summary.pointCount;
    out.segmentCount = summary.segmentCount;
    out.cacheIdentity = summary.cacheIdentity.getString();
    out.sourceIdentity = summary.sourceIdentity.getString();
    out.cacheIdentityValue = summary.cacheIdentityValue;
    out.sourceIdentityValue = summary.sourceIdentityValue;
}

static void
qg_obol_export_point_summary_from_action(
	const SoBRLExportAction::ObjectPointSummary &summary,
	QgObolExportObjectPointSummary &out)
{
    out.valid = summary.valid != 0;
    out.pointCount = summary.pointCount;
    out.cacheIdentity = summary.cacheIdentity.getString();
    out.sourceIdentity = summary.sourceIdentity.getString();
    out.cacheIdentityValue = summary.cacheIdentityValue;
    out.sourceIdentityValue = summary.sourceIdentityValue;
}

static void
qg_obol_export_surface_summary_from_action(
	const SoBRLExportAction::ObjectSurfaceSummary &summary,
	QgObolExportObjectSurfaceSummary &out)
{
    out.valid = summary.valid != 0;
    out.pointCount = summary.pointCount;
    out.normalCount = summary.normalCount;
    out.indexCount = summary.indexCount;
    out.faceCount = summary.faceCount;
    out.normalsPerIndex = summary.normalsPerIndex;
    out.materialColorValid = summary.materialColorValid != 0;
    out.materialColor = summary.materialColor;
    out.materialDrawMode = summary.materialDrawMode;
    out.materialHighlighted = summary.materialHighlighted != 0;
    out.cacheIdentity = summary.cacheIdentity.getString();
    out.sourceIdentity = summary.sourceIdentity.getString();
    out.cacheIdentityValue = summary.cacheIdentityValue;
    out.sourceIdentityValue = summary.sourceIdentityValue;
}

static void
qg_obol_export_object_detail_from_action(
	const SoBRLExportAction &exportAction,
	const SoBRLExportAction::ObjectRecord &object,
	QgObolExportObjectRecord &out)
{
    SoBRLExportAction::ObjectLineSummary lineSummary;
    if (exportAction.getObjectRecordLineSummary(object, lineSummary)) {
	qg_obol_export_line_summary_from_action(lineSummary,
		out.lineSummary);
	out.linePoints.reserve(lineSummary.pointCount);
	out.lineCommands.reserve(lineSummary.pointCount);
	for (size_t i = 0; i < lineSummary.pointCount; i++) {
	    SbVec3f point;
	    int command = -1;
	    if (exportAction.getObjectRecordLinePoint(object, i, point)) {
		out.linePoints.push_back(point);
		(void)exportAction.getObjectRecordLineCommand(object, i,
			command);
		out.lineCommands.push_back(command);
	    }
	}
    }

    SoBRLExportAction::ObjectPointSummary pointSummary;
    if (exportAction.getObjectRecordPointSummary(object, pointSummary)) {
	qg_obol_export_point_summary_from_action(pointSummary,
		out.pointSummary);
	out.points.reserve(pointSummary.pointCount);
	for (size_t i = 0; i < pointSummary.pointCount; i++) {
	    SbVec3f point;
	    if (exportAction.getObjectRecordPoint(object, i, point))
		out.points.push_back(point);
	}
    }

    SoBRLExportAction::ObjectSurfaceSummary surfaceSummary;
    if (exportAction.getObjectRecordSurfaceSummary(object, surfaceSummary)) {
	qg_obol_export_surface_summary_from_action(surfaceSummary,
		out.surfaceSummary);
	(void)exportAction.getObjectRecordSurfaceDetail(object,
		&out.surfacePoints, &out.surfaceIndices);
    }
}

static void
qg_obol_export_object_from_action(const SoBRLExportAction &exportAction,
	const SoBRLExportAction::ObjectRecord &object,
	QgObolExportObjectRecord &out)
{
    out.path = object.path.getString();
    out.sourceName = object.sourceName.getString();
    out.sourceType = object.sourceType.getString();
    out.displayName = object.displayName.getString();
    out.geometryName = object.geometryName.getString();
    out.cacheIdentity = object.cacheIdentity.getString();
    out.sourceIdentity = object.sourceIdentity.getString();
    out.cacheIdentityValue = object.cacheIdentityValue;
    out.sourceIdentityValue = object.sourceIdentityValue;
    out.sourceId = object.sourceId;
    out.databaseIntent = object.databaseIntent != 0;
    out.overlayIntent = object.overlayIntent != 0;
    out.hudIntent = object.hudIntent != 0;
    out.localSource = object.localSource != 0;
    out.sharedSource = object.sharedSource != 0;
    out.nonDatabaseSource = object.nonDatabaseSource != 0;
    out.drawMode = object.drawMode;
    out.recordRole = object.recordRole.getString();
    out.geometryKind = object.geometryKind.getString();
    out.selected = object.selected != 0;
    out.highlighted = object.highlighted != 0;
    out.visible = object.visible != 0;
    out.colorOverride = object.colorOverride != 0;
    out.color = object.color;
    if (!object.bounds.isEmpty()) {
	out.boundsValid = true;
	out.boundsMin = object.bounds.getMin();
	out.boundsMax = object.bounds.getMax();
    }
    out.lineCount = object.lineIndices.size();
    out.pointPrimitiveCount = object.pointIndices.size();
    out.triangleCount = object.triangleIndices.size();
    qg_obol_export_object_detail_from_action(exportAction, object, out);
}

static void
qg_obol_export_object_records_from_action(
	const SoBRLExportAction &exportAction,
	const QgObolExportObjectQuery &query,
	QgObolExportObjectQueryResult &record)
{
    record.lineCount = exportAction.getLineCount();
    record.pointCount = exportAction.getPointCount();
    record.triangleCount = exportAction.getTriangleCount();

    std::vector<SoBRLExportAction::ObjectRecord> objects;
    const char *glob = query.glob.empty() ? NULL : query.glob.c_str();
    exportAction.collectObjectRecords(objects,
	    qg_obol_export_object_query_flags(query.flags), glob,
	    query.drawMode);
    record.objects.clear();
    record.objects.reserve(objects.size());
    for (size_t i = 0; i < objects.size(); i++) {
	QgObolExportObjectRecord out;
	qg_obol_export_object_from_action(exportAction, objects[i], out);
	record.objects.push_back(out);
    }
}

int
qg_obol_export_geometry_full_detail(QgView *display,
	QgObolExportGeometryRecord &record)
{
    record = QgObolExportGeometryRecord();
    if (!display)
	return 0;

    BObolViewController *controller = display->obolViewController();
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

int
qg_obol_export_object_records(QgView *display,
	const QgObolExportObjectQuery &query,
	QgObolExportObjectQueryResult &record)
{
    record = QgObolExportObjectQueryResult();
    if (!display)
	return 0;

    BObolViewController *controller = display->obolViewController();
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;

    SoBRLExportAction exportAction;
    exportAction.setGeometryPolicy(qg_obol_export_geometry_policy(
	    query.geometryPolicy));
    exportAction.apply(controller->getViewport()->getRoot());
    record.submittedSourceRequestCount =
	qg_obol_export_consume_source_full_detail(controller, exportAction);
    record.sourceFullDetailPending =
	record.submittedSourceRequestCount > 0;

    qg_obol_export_object_records_from_action(exportAction, query, record);
    return record.objects.empty() ? 0 : 1;
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
