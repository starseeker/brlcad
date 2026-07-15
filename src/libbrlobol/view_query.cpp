/*                   V I E W _ Q U E R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_query.cpp */

#include "common.h"

#include "brlobol/view_query.h"

#include "brlobol/view_controller.h"

#include <Inventor/SbLine.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SoPath.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/details/SoDetail.h>
#include <Inventor/lists/SoPickedPointList.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <float.h>
#include <string>
#include <vector>

BRLObolViewPickRecord::BRLObolViewPickRecord() :
    point(0.0f, 0.0f, 0.0f),
    distance(FLT_MAX)
{
}

BRLObolViewSnapRecord::BRLObolViewSnapRecord() :
    point(0.0f, 0.0f, 0.0f),
    path(),
    editIntentId(),
    editIntentRole(),
    kind(SoBRLSnapAction::NONE),
    primitiveIndex(-1),
    vertexIndex(-1),
    edgeSlot(-1),
    edgeVertexIndexA(-1),
    edgeVertexIndexB(-1),
    submittedSourceRequestCount(0),
    distance(FLT_MAX),
    sourceFullDetailPending(FALSE)
{
}

static const SoBRLPickDetail *
brlobol_view_pick_detail(const SoPickedPoint *pickedPoint)
{
    if (!pickedPoint)
	return NULL;

    const SoDetail *detail = pickedPoint->getDetail();
    if (detail && detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	return static_cast<const SoBRLPickDetail *>(detail);

    SoPath *path = pickedPoint->getPath();
    if (!path)
	return NULL;

    for (int i = path->getLength() - 1; i >= 0; i--) {
	SoNode *node = path->getNode(i);
	if (!node)
	    continue;
	detail = pickedPoint->getDetail(node);
	if (detail && detail->isOfType(SoBRLPickDetail::getClassTypeId()))
	    return static_cast<const SoBRLPickDetail *>(detail);
    }

    return NULL;
}

static BRLObolViewPickRecord
brlobol_view_pick_record(const SoBRLPickDetail &detail,
	const SbVec3f &point,
	float distance)
{
    BRLObolViewPickRecord record;
    record.detail = detail;
    record.point = point;
    record.distance = distance;
    return record;
}

static BRLObolViewPickRecord
brlobol_view_pick_record(const SoPickedPoint *pickedPoint,
	const SoBRLPickDetail *detail,
	const SbVec3f &rayOrigin)
{
    if (!detail)
	return BRLObolViewPickRecord();

    const SbVec3f point = pickedPoint ? pickedPoint->getPoint() :
	SbVec3f(0.0f, 0.0f, 0.0f);
    return brlobol_view_pick_record(*detail, point,
	    (point - rayOrigin).length());
}

static BRLObolViewPickRecord
brlobol_view_pick_record(const BRLObolSourceMeshPickResult &pick)
{
    return brlobol_view_pick_record(pick.detail, pick.point, pick.distance);
}

static BRLObolViewPickRecord
brlobol_view_pick_record(const BRLObolRtPickResult &pick)
{
    return brlobol_view_pick_record(pick.detail, pick.point, pick.distance);
}

static std::string
brlobol_view_pick_normalized_path(const SoBRLPickDetail &detail)
{
    const char *path = detail.getPath().getString();
    if ((!path || !path[0]) && detail.getSourceName().getString())
	path = detail.getSourceName().getString();
    std::string normalized = path ? path : "";
    const size_t start = normalized.find_first_not_of('/');
    return start == std::string::npos ? std::string() :
	normalized.substr(start);
}

static void
brlobol_view_pick_insert(std::vector<BRLObolViewPickRecord> &records,
	const BRLObolViewPickRecord &record,
	bool pickAll)
{
    if (brlobol_view_pick_normalized_path(record.detail).empty())
	return;

    if (pickAll) {
	records.push_back(record);
	return;
    }

    if (records.empty() || record.distance < records[0].distance)
	records.assign(1, record);
}

static bool
brlobol_view_pick_path_already_recorded(
	const std::vector<BRLObolViewPickRecord> &records,
	const BRLObolViewPickRecord &candidate)
{
    const std::string candidatePath =
	brlobol_view_pick_normalized_path(candidate.detail);
    if (candidatePath.empty())
	return false;

    for (const BRLObolViewPickRecord &record : records) {
	if (brlobol_view_pick_normalized_path(record.detail) == candidatePath)
	    return true;
    }
    return false;
}

static bool
brlobol_view_pick_record_nearer(const BRLObolViewPickRecord &a,
	const BRLObolViewPickRecord &b)
{
    return a.distance < b.distance;
}

static void
brlobol_view_pick_display(BRLObolViewController *controller,
	SoRayPickAction &action,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records)
{
    action.setPickAll(pickAll ? TRUE : FALSE);
    action.apply(controller->getViewport()->getRoot());

    const SbVec3f &origin = action.getLine().getPosition();
    if (pickAll) {
	const SoPickedPointList &points = action.getPickedPointList();
	records.reserve(records.size() + static_cast<size_t>(points.getLength()));
	for (int i = 0; i < points.getLength(); i++) {
	    const SoPickedPoint *point = points[i];
	    const SoBRLPickDetail *detail = brlobol_view_pick_detail(point);
	    if (detail)
		brlobol_view_pick_insert(records,
		    brlobol_view_pick_record(point, detail, origin), pickAll);
	}
	return;
    }

    const SoPickedPoint *point = action.getPickedPoint();
    const SoBRLPickDetail *detail = brlobol_view_pick_detail(point);
    if (detail)
	brlobol_view_pick_insert(records,
	    brlobol_view_pick_record(point, detail, origin), pickAll);
}

static int
brlobol_view_pick_source_exact(BRLObolViewController *controller,
	const SbLine &line,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    if (!controller)
	return 0;

    BRLObolSourceMeshPickResult sourcePick;
    int submitted = 0;
    const int hit = controller->pickSourceMeshExactRay(sourcePick,
	line.getPosition(), line.getDirection(), 0, &submitted) > 0 &&
	sourcePick.hit;
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = submitted;
    if (!hit)
	return 0;

    brlobol_view_pick_insert(records, brlobol_view_pick_record(sourcePick),
	pickAll);
    return 1;
}

static void
brlobol_view_pick_rt_exact(BRLObolViewController *controller,
	const SbLine &line,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records)
{
    if (!controller)
	return;

    std::vector<BRLObolRtPickResult> rtPicks;
    controller->pickRtExactRay(rtPicks, line.getPosition(),
	line.getDirection(), pickAll ? TRUE : FALSE);
    for (const BRLObolRtPickResult &pick : rtPicks) {
	BRLObolViewPickRecord record = brlobol_view_pick_record(pick);
	if (brlobol_view_pick_path_already_recorded(records, record))
	    continue;
	brlobol_view_pick_insert(records, record, pickAll);
    }
}

static SbBool
brlobol_view_pick_camera_line(BRLObolViewController *controller,
	int viewportX,
	int viewportY,
	SbLine &line)
{
    if (!controller || !controller->getCamera())
	return FALSE;

    const SbViewportRegion &region = controller->getViewportRegion();
    const SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return FALSE;

    const float nx = std::max(0.0f, std::min(1.0f,
	static_cast<float>(viewportX) / static_cast<float>(size[0])));
    const float ny = std::max(0.0f, std::min(1.0f,
	static_cast<float>(viewportY + 1) / static_cast<float>(size[1])));
    const float aspect = static_cast<float>(size[0]) /
	static_cast<float>(size[1]);
    controller->getCamera()->getViewVolume(aspect).projectPointToLine(
	SbVec2f(nx, ny), line);
    return TRUE;
}

int
brlobol_view_pick_point(BRLObolViewController *controller,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    records.clear();
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return 0;

    const SbViewportRegion &region = controller->getViewportRegion();
    const SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return 0;

    const int viewportX = std::max(0, std::min(x, static_cast<int>(size[0]) - 1));
    const int viewportY = std::max(0, std::min(static_cast<int>(size[1]) - 1 - y,
	static_cast<int>(size[1]) - 1));
    SoRayPickAction action(region);
    action.setPoint(SbVec2s(static_cast<short>(viewportX),
	static_cast<short>(viewportY)));
    action.setRadius(radiusPixels > 0.0f ? radiusPixels : 1.0f);
    brlobol_view_pick_display(controller, action, pickAll, records);

    const SbLine &line = action.getLine();
    brlobol_view_pick_source_exact(controller, line, pickAll, records,
	submittedSourceRequestCount);
    SbLine rtLine = line;
    if (records.empty())
	(void)brlobol_view_pick_camera_line(controller, viewportX, viewportY,
	    rtLine);
    brlobol_view_pick_rt_exact(controller, rtLine, pickAll, records);
    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    brlobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

int
brlobol_view_pick_ray(BRLObolViewController *controller,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	bool pickAll,
	std::vector<BRLObolViewPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    records.clear();
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return 0;

    SoRayPickAction action(controller->getViewportRegion());
    action.setRay(rayOrigin, rayDirection);
    brlobol_view_pick_display(controller, action, pickAll, records);

    SbLine line;
    line.setPosDir(rayOrigin, rayDirection);
    brlobol_view_pick_source_exact(controller, line, pickAll, records,
	submittedSourceRequestCount);
    brlobol_view_pick_rt_exact(controller, line, pickAll, records);
    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    brlobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

static void
brlobol_view_snap_record_from_action(const SoBRLSnapAction &action,
	BRLObolViewSnapRecord &record)
{
    record.point = action.getPoint();
    record.path = action.getPath();
    record.editIntentId = action.getEditIntentId();
    record.editIntentRole = action.getEditIntentRole();
    record.kind = action.getKind();
    record.primitiveIndex = action.getPrimitiveIndex();
    record.vertexIndex = action.getVertexIndex();
    record.edgeSlot = action.getEdgeSlot();
    record.edgeVertexIndexA = action.getEdgeVertexIndexA();
    record.edgeVertexIndexB = action.getEdgeVertexIndexB();
    record.distance = action.getDistance();
}

int
brlobol_view_snap_point(BRLObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
	bool consumeSourceFullDetail,
	BRLObolViewSnapRecord &record)
{
    record = BRLObolViewSnapRecord();
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return 0;

    SoBRLSnapAction action;
    action.setQueryPoint(query);
    action.setTolerance(tolerance > 0.0f ? tolerance : 1.0f);
    action.setEnabledKinds(enabledKinds ? enabledKinds :
	static_cast<uint32_t>(SoBRLSnapAction::ALL_KINDS));
    action.setPriorityPolicy(SoBRLSnapAction::FEATURE_PRIORITY);
    action.setGeometryPolicy(geometryPolicy);
    action.apply(controller->getViewport()->getRoot());

    if (geometryPolicy == SoBRLSnapAction::FULL_DETAIL &&
	consumeSourceFullDetail) {
	(void)controller->consumeSnapSourceFullDetail(action, 0,
	    &record.submittedSourceRequestCount);
	record.sourceFullDetailPending =
	record.submittedSourceRequestCount > 0 ? TRUE : FALSE;
    }

    if (!action.hasCandidate())
	return 0;

    const int submitted = record.submittedSourceRequestCount;
    const SbBool pending = record.sourceFullDetailPending;
    brlobol_view_snap_record_from_action(action, record);
    record.submittedSourceRequestCount = submitted;
    record.sourceFullDetailPending = pending;
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
