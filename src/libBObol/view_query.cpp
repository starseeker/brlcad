/*                   V I E W _ Q U E R Y . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file view_query.cpp */

#include "common.h"

#include "BObol/BViewQuery.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BViewController.h"
#include "BObol/BViewStore.h"

#include <Inventor/SbLine.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SoPath.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/details/SoDetail.h>
#include <Inventor/lists/SoPickedPointList.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>

#include <algorithm>
#include <float.h>
#include <string>
#include <vector>

BObolViewPickRecord::BObolViewPickRecord() :
    point(0.0f, 0.0f, 0.0f),
    distance(FLT_MAX)
{
}

BObolViewSnapRecord::BObolViewSnapRecord() :
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
bobol_view_pick_detail(const SoPickedPoint *pickedPoint)
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

static BObolViewPickRecord
bobol_view_pick_record(const SoBRLPickDetail &detail,
	const SbVec3f &point,
	float distance)
{
    BObolViewPickRecord record;
    record.detail = detail;
    record.point = point;
    record.distance = distance;
    return record;
}

static BObolViewPickRecord
bobol_view_pick_record(const SoPickedPoint *pickedPoint,
	const SoBRLPickDetail *detail,
	const SbVec3f &rayOrigin)
{
    if (!detail)
	return BObolViewPickRecord();

    const SbVec3f point = pickedPoint ? pickedPoint->getPoint() :
	SbVec3f(0.0f, 0.0f, 0.0f);
    return bobol_view_pick_record(*detail, point,
	    (point - rayOrigin).length());
}

static BObolViewPickRecord
bobol_view_pick_record(const BObolSourceMeshPickResult &pick)
{
    return bobol_view_pick_record(pick.detail, pick.point, pick.distance);
}

static BObolViewPickRecord
bobol_view_pick_record(const BObolRtPickResult &pick)
{
    return bobol_view_pick_record(pick.detail, pick.point, pick.distance);
}

static std::string
bobol_view_pick_normalized_path(const SoBRLPickDetail &detail)
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
bobol_view_pick_insert(std::vector<BObolViewPickRecord> &records,
	const BObolViewPickRecord &record,
	bool pickAll)
{
    if (bobol_view_pick_normalized_path(record.detail).empty())
	return;

    if (pickAll) {
	records.push_back(record);
	return;
    }

    if (records.empty() || record.distance < records[0].distance)
	records.assign(1, record);
}

static bool
bobol_view_pick_path_already_recorded(
	const std::vector<BObolViewPickRecord> &records,
	const BObolViewPickRecord &candidate)
{
    const std::string candidatePath =
	bobol_view_pick_normalized_path(candidate.detail);
    if (candidatePath.empty())
	return false;

    for (const BObolViewPickRecord &record : records) {
	if (bobol_view_pick_normalized_path(record.detail) == candidatePath)
	    return true;
    }
    return false;
}

static bool
bobol_view_pick_record_nearer(const BObolViewPickRecord &a,
	const BObolViewPickRecord &b)
{
    return a.distance < b.distance;
}

static void
bobol_view_pick_display_node(SoNode *root,
	SoRayPickAction &action,
	const SbVec3f &worldRayOrigin,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records)
{
    if (!root)
	return;
    action.setPickAll(pickAll ? TRUE : FALSE);
    action.apply(root);

    if (pickAll) {
	const SoPickedPointList &points = action.getPickedPointList();
	records.reserve(records.size() + static_cast<size_t>(points.getLength()));
	for (int i = 0; i < points.getLength(); i++) {
	    const SoPickedPoint *point = points[i];
	    const SoBRLPickDetail *detail = bobol_view_pick_detail(point);
	    if (detail)
		bobol_view_pick_insert(records,
		    bobol_view_pick_record(point, detail, worldRayOrigin),
		    pickAll);
	}
	return;
    }

    const SoPickedPoint *point = action.getPickedPoint();
    const SoBRLPickDetail *detail = bobol_view_pick_detail(point);
    if (detail)
	bobol_view_pick_insert(records,
	    bobol_view_pick_record(point, detail, worldRayOrigin), pickAll);
}

static void
bobol_view_pick_display(BObolViewController *controller,
	SoRayPickAction &action,
	const SbVec3f &worldRayOrigin,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records)
{
    if (!controller || !controller->getViewport())
	return;
    bobol_view_pick_display_node(controller->getViewport()->getRoot(),
	action, worldRayOrigin, pickAll, records);
}

namespace {

struct bobol_feature_point_query {
    SoCamera *camera;
    const SbViewportRegion *region;
    SbVec2s point;
    float radius;
    SbVec3f rayOrigin;
    SbVec3f rayDirection;
    bool pickAll;
    std::vector<BObolViewPickRecord> *records;
};

static int
bobol_view_pick_feature_node(BObolFeatureHandle UNUSED(handle), SoNode *node,
	void *data)
{
    bobol_feature_point_query *query =
	static_cast<bobol_feature_point_query *>(data);
    if (!query || !query->camera || !query->region || !query->records ||
	!node)
	return 1;

    /* Store feature roots are self-contained model-space values attached
     * directly to the view feature group.  Pair each with the active camera
     * in a transient query graph so line/point pixel tolerances get complete
     * camera state without traversing the CAD aggregate. */
    SoSeparator *pickRoot = new SoSeparator;
    pickRoot->ref();
    pickRoot->addChild(query->camera);
    pickRoot->addChild(node);
    SoRayPickAction action(*query->region);
    action.setPoint(query->point);
    action.setRadius(query->radius);
    bobol_view_pick_display_node(pickRoot, action, query->rayOrigin,
	query->pickAll, *query->records);
    pickRoot->unref();
    return 1;
}

static int
bobol_view_pick_features_point(BObolViewController *controller,
	const SbViewportRegion &region, int viewportX, int viewportY,
	float radiusPixels, const SbVec3f &worldRayOrigin,
	const SbVec3f &worldRayDirection, bool pickAll,
	std::vector<BObolViewPickRecord> &records)
{
    if (!controller)
	return 0;
    bobol_feature_point_query query;
    query.camera = controller->getCamera();
    query.region = &region;
    query.point = SbVec2s(static_cast<short>(viewportX),
	static_cast<short>(viewportY));
    query.radius = radiusPixels > 0.0f ? radiusPixels : 1.0f;
    query.rayOrigin = worldRayOrigin;
    query.rayDirection = worldRayDirection;
    query.pickAll = pickAll;
    query.records = &records;
    controller->features().visitNodes(bobol_view_pick_feature_node, &query);
    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    bobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

static int
bobol_view_pick_source_identities(BObolViewController *controller,
	const SbLine &worldLine, bool pickAll,
	std::vector<BObolViewPickRecord> &records)
{
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return 0;

    SoBRLSourceMeshPickAction action;
    action.setRay(worldLine.getPosition(), worldLine.getDirection());
    action.apply(controller->getViewport()->getRoot());
    for (int i = 0;
	 i < action.getSourceBackedFullDetailRequestCount(); ++i) {
	const BObolSourceMeshRequest &request =
	    action.getSourceBackedFullDetailRequest(i);
	if (request.path.getLength() == 0)
	    continue;
	SbVec3f localPoint = request.bounds.isEmpty() ?
	    SbVec3f(0.0f, 0.0f, 0.0f) : request.bounds.getCenter();
	SbVec3f worldPoint;
	request.localToWorld.multVecMatrix(localPoint, worldPoint);
	BObolViewPickRecord record;
	record.point = worldPoint;
	record.distance = (worldPoint - worldLine.getPosition()).length();
	record.detail.setPath(request.path);
	record.detail.setSourceInstanceKey(request.ownerSourceInstanceKey);
	record.detail.setSourceName(request.sourceName);
	record.detail.setSourceType(request.sourceType);
	record.detail.setSourceId(request.sourceId);
	record.detail.setRegionId(request.regionId);
	record.detail.setAirCode(request.airCode);
	record.detail.setMaterialId(request.materialId);
	record.detail.setLos(request.los);
	record.detail.setMaterialColor(request.materialColorValid ? TRUE : FALSE,
	    request.materialColor);
	record.detail.setMaterialShader(request.materialShader);
	record.detail.setEditIntent(request.editIntentId,
	    request.editIntentRole);
	record.detail.setModelPoint(localPoint);
	record.detail.setPrimitive(SoBRLPickDetail::UNKNOWN, -1);
	bobol_view_pick_insert(records, record, pickAll);
    }
    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    bobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

}

static int
bobol_view_pick_source_exact(BObolViewController *controller,
	const SbLine &line,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    if (!controller)
	return 0;

    BObolSourceMeshPickResult sourcePick;
    int submitted = 0;
    const int hit = controller->pickSourceMeshExactRay(sourcePick,
	line.getPosition(), line.getDirection(), 0, &submitted) > 0 &&
	sourcePick.hit;
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = submitted;
    if (!hit)
	return 0;

    bobol_view_pick_insert(records, bobol_view_pick_record(sourcePick),
	pickAll);
    return 1;
}

static void
bobol_view_pick_rt_exact(BObolViewController *controller,
	const SbLine &line,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records)
{
    if (!controller)
	return;

    std::vector<BObolRtPickResult> rtPicks;
    controller->pickRtExactRay(rtPicks, line.getPosition(),
	line.getDirection(), pickAll ? TRUE : FALSE);
    for (const BObolRtPickResult &pick : rtPicks) {
	BObolViewPickRecord record = bobol_view_pick_record(pick);
	if (bobol_view_pick_path_already_recorded(records, record))
	    continue;
	bobol_view_pick_insert(records, record, pickAll);
    }
}

static SbBool
bobol_view_pick_camera_line(BObolViewController *controller,
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
	static_cast<float>(viewportY) / static_cast<float>(size[1])));
    const float aspect = static_cast<float>(size[0]) /
	static_cast<float>(size[1]);
    controller->getCamera()->getViewVolume(aspect).projectPointToLine(
	SbVec2f(nx, ny), line);
    return TRUE;
}

int
bobol_view_pick_display_point(BObolViewController *controller,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records)
{
    records.clear();
    if (!controller || !controller->getViewport() ||
	!controller->getViewport()->getRoot())
	return 0;

    const SbViewportRegion &region = controller->getViewportRegion();
    const SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return 0;

    const int viewportX = std::max(0, std::min(x,
	static_cast<int>(size[0]) - 1));
    const int viewportY = std::max(0, std::min(
	static_cast<int>(size[1]) - 1 - y,
	static_cast<int>(size[1]) - 1));
    SbLine worldLine;
    if (!bobol_view_pick_camera_line(controller, viewportX, viewportY,
	    worldLine))
	return 0;

    /* Explicit view features (edit manipulators, command results, and other
     * plugin nodes) take interaction priority over database geometry.  Query
     * their retained roots directly, avoiding a traversal of the CAD source. */
    if (bobol_view_pick_features_point(controller, region, viewportX,
	    viewportY, radiusPixels, worldLine.getPosition(),
	    worldLine.getDirection(), pickAll,
	    records) > 0 && !pickAll)
	return static_cast<int>(records.size());

    /* Object/path selection is a screen-space identity query, not an exact
     * primitive query.  Prefer the retained occurrence registry so the first
     * click cannot trigger construction of every triangle-pick and librt BVH
     * in a large scene. */
    const int extent = std::max(1, static_cast<int>(std::ceil(
	radiusPixels > 0.0f ? radiusPixels : 1.0f)));
    std::vector<BObolViewPickRecord> boundedRecords;
    const int bounded = bobol_view_pick_rectangle(controller,
	x - extent, y - extent, x + extent, y + extent, boundedRecords);
    if (!boundedRecords.empty()) {
	if (!pickAll && boundedRecords.size() > 1)
	    boundedRecords.resize(1);
	records.insert(records.end(), boundedRecords.begin(),
	    boundedRecords.end());
	if (!pickAll)
	    return static_cast<int>(records.size());
    }

    const int sourceIdentities = bobol_view_pick_source_identities(controller,
	worldLine, pickAll, records);
    if (sourceIdentities > 0 && !pickAll)
	return static_cast<int>(records.size());

    /* Cold/non-mesh database sources may not have occurrence or mesh
     * identities yet.  Their conservative cached target bound is still an
     * immediate object-selection answer and must not initialize exact librt
     * or triangle-pick acceleration state. */
    std::vector<BObolViewPickRecord> sourceRecords;
    const int sourceBounded = bobol_view_pick_source_rectangle(controller,
	x - extent, y - extent, x + extent, y + extent, sourceRecords);
    if (!sourceRecords.empty()) {
	if (!pickAll && sourceRecords.size() > 1)
	    sourceRecords.resize(1);
	records.insert(records.end(), sourceRecords.begin(),
	    sourceRecords.end());
    }
    if (bounded >= 0 || sourceBounded >= 0)
	return static_cast<int>(records.size());

    SoRayPickAction action(region);
    action.setPoint(SbVec2s(static_cast<short>(viewportX),
	static_cast<short>(viewportY)));
    action.setRadius(radiusPixels > 0.0f ? radiusPixels : 1.0f);
    bobol_view_pick_display(controller, action, worldLine.getPosition(),
	pickAll, records);

    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    bobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

int
bobol_view_pick_point(BObolViewController *controller,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records,
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

    const int viewportX = std::max(0, std::min(x,
	static_cast<int>(size[0]) - 1));
    const int viewportY = std::max(0, std::min(
	static_cast<int>(size[1]) - 1 - y,
	static_cast<int>(size[1]) - 1));
    SbLine worldLine;
    if (!bobol_view_pick_camera_line(controller, viewportX, viewportY,
	    worldLine))
	return 0;
    SoRayPickAction action(region);
    action.setPoint(SbVec2s(static_cast<short>(viewportX),
	static_cast<short>(viewportY)));
    action.setRadius(radiusPixels > 0.0f ? radiusPixels : 1.0f);
    bobol_view_pick_display(controller, action, worldLine.getPosition(),
	pickAll, records);

    bobol_view_pick_source_exact(controller, worldLine, pickAll, records,
	submittedSourceRequestCount);
    bobol_view_pick_rt_exact(controller, worldLine, pickAll, records);
    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    bobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

namespace {

struct bobol_rectangle_query_state {
    SbMatrix viewProjection;
    float minimumX;
    float minimumY;
    float maximumX;
    float maximumY;
    bool handled;
    bool sourceBounds;
    std::vector<BObolViewPickRecord> *records;
};

static SoCallbackAction::Response
bobol_rectangle_query_callback(void *data, SoCallbackAction *action,
	const SoNode *node)
{
    bobol_rectangle_query_state *state =
	static_cast<bobol_rectangle_query_state *>(data);
    if (!state || !action || !node || !state->records ||
	!node->isOfType(SoBRLDatabaseSource::getClassTypeId()))
	return SoCallbackAction::CONTINUE;

    const SoBRLDatabaseSource *source =
	static_cast<const SoBRLDatabaseSource *>(node);
    const int count = state->sourceBounds ?
	source->querySourceRectangle(action->getModelMatrix(),
	    state->viewProjection, state->minimumX, state->minimumY,
	    state->maximumX, state->maximumY, *state->records) :
	source->queryCompactRectangle(action->getModelMatrix(),
	    state->viewProjection, state->minimumX, state->minimumY,
	    state->maximumX, state->maximumY, *state->records);
    if (count < 0)
	return SoCallbackAction::CONTINUE;
    state->handled = true;
    return SoCallbackAction::PRUNE;
}

}

static int
bobol_view_pick_rectangle_internal(BObolViewController *controller,
	int x0,
	int y0,
	int x1,
	int y1,
	std::vector<BObolViewPickRecord> &records,
	bool sourceBounds)
{
    records.clear();
    if (!controller || !controller->getCamera() ||
	!controller->getViewport() || !controller->getViewport()->getRoot())
	return -1;

    const SbViewportRegion &region = controller->getViewportRegion();
    const SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return -1;

    const int minimumPixelX = std::max(0, std::min(std::min(x0, x1),
	static_cast<int>(size[0]) - 1));
    const int maximumPixelX = std::max(0, std::min(std::max(x0, x1),
	static_cast<int>(size[0]) - 1));
    const int minimumPixelY = std::max(0, std::min(std::min(y0, y1),
	static_cast<int>(size[1]) - 1));
    const int maximumPixelY = std::max(0, std::min(std::max(y0, y1),
	static_cast<int>(size[1]) - 1));
    const float inverseWidth = 1.0f /
	static_cast<float>(std::max(1, static_cast<int>(size[0]) - 1));
    const float inverseHeight = 1.0f /
	static_cast<float>(std::max(1, static_cast<int>(size[1]) - 1));

    bobol_rectangle_query_state state;
    const float aspect = static_cast<float>(size[0]) /
	static_cast<float>(size[1]);
    state.viewProjection = controller->getCamera()->
	getViewVolume(aspect).getMatrix();
    state.minimumX = 2.0f * minimumPixelX * inverseWidth - 1.0f;
    state.maximumX = 2.0f * maximumPixelX * inverseWidth - 1.0f;
    state.minimumY = 1.0f - 2.0f * maximumPixelY * inverseHeight;
    state.maximumY = 1.0f - 2.0f * minimumPixelY * inverseHeight;
    state.handled = false;
    state.sourceBounds = sourceBounds;
    state.records = &records;

    SoCallbackAction action(region);
    action.addPreCallback(SoBRLDatabaseSource::getClassTypeId(),
	bobol_rectangle_query_callback, &state);
    action.apply(controller->getViewport()->getRoot());
    if (!state.handled) {
	records.clear();
	return -1;
    }

    std::stable_sort(records.begin(), records.end(),
	bobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

int
bobol_view_pick_rectangle(BObolViewController *controller,
	int x0, int y0, int x1, int y1,
	std::vector<BObolViewPickRecord> &records)
{
    return bobol_view_pick_rectangle_internal(controller, x0, y0, x1, y1,
	records, false);
}

int
bobol_view_pick_source_rectangle(BObolViewController *controller,
	int x0, int y0, int x1, int y1,
	std::vector<BObolViewPickRecord> &records)
{
    return bobol_view_pick_rectangle_internal(controller, x0, y0, x1, y1,
	records, true);
}

int
bobol_view_pick_ray(BObolViewController *controller,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	bool pickAll,
	std::vector<BObolViewPickRecord> &records,
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
    bobol_view_pick_display(controller, action, rayOrigin, pickAll, records);

    SbLine line;
    line.setPosDir(rayOrigin, rayDirection);
    bobol_view_pick_source_exact(controller, line, pickAll, records,
	submittedSourceRequestCount);
    bobol_view_pick_rt_exact(controller, line, pickAll, records);
    if (pickAll)
	std::stable_sort(records.begin(), records.end(),
	    bobol_view_pick_record_nearer);
    return static_cast<int>(records.size());
}

static void
bobol_view_snap_record_from_action(const SoBRLSnapAction &action,
	BObolViewSnapRecord &record)
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
bobol_view_snap_point(BObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
	bool consumeSourceFullDetail,
	BObolViewSnapRecord &record)
{
    return bobol_view_snap_point_filtered(controller, query, tolerance,
	    enabledKinds, SoBRLSnapAction::ALL_SOURCES, geometryPolicy,
	    consumeSourceFullDetail, record);
}

int
bobol_view_snap_point_filtered(BObolViewController *controller,
	const SbVec3f &query,
	float tolerance,
	uint32_t enabledKinds,
	uint32_t sourceFilter,
	SoBRLSnapAction::GeometryPolicy geometryPolicy,
	bool consumeSourceFullDetail,
	BObolViewSnapRecord &record)
{
    record = BObolViewSnapRecord();
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
    action.setSourceFilter(sourceFilter);
    BObolPolygonRecord excludedPolygon;
    const BObolPolygonHandle excluded = controller->polygons().snapExclude();
    if (excluded.isValid() &&
	controller->polygons().record(excluded, excludedPolygon))
	action.setExcludedPath(excludedPolygon.name);
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
    bobol_view_snap_record_from_action(action, record);
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
