/*                  Q G O B O L P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolPick.cpp */

#include "common.h"

#include "qtcad/QgObolPick.h"

#include "brlobol/pick_detail.h"
#include "brlobol/view_controller.h"
#include "brlobol/view_store.h"
#include "qtcad/QgView.h"
#include "rt/view.h"
#include "QgLegacyViewContext.h"

#include <Inventor/SbLine.h>
#include <Inventor/SbViewVolume.h>
#include <Inventor/SbVec2f.h>
#include <Inventor/SoPath.h>
#include <Inventor/SoPickedPoint.h>
#include <Inventor/SoViewport.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/lists/SoPickedPointList.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <cmath>
#include <cstdio>
#include <float.h>
#include <string>
#include <unordered_set>

QgObolPickRecord::QgObolPickRecord(void) :
    path(),
    sourceName(),
    sourceType(),
    featureName(),
    materialShader(),
    editIntentId(),
    editIntentRole(),
    point(0.0f, 0.0f, 0.0f),
    modelPoint(0.0f, 0.0f, 0.0f),
    materialColor(1.0f, 1.0f, 1.0f),
    distance(FLT_MAX),
    sourceId(0),
    regionId(0),
    airCode(0),
    materialId(0),
    los(0),
    primitiveKind(UNKNOWN),
    primitiveIndex(-1),
    featurePrimitiveIndex(-1),
    faceVertexIndexA(-1),
    faceVertexIndexB(-1),
    faceVertexIndexC(-1),
    nearestFaceEdgeSlot(-1),
    nearestFaceEdgeVertexIndexA(-1),
    nearestFaceEdgeVertexIndexB(-1),
    nearestFaceVertexSlot(-1),
    nearestFaceVertexIndex(-1),
    materialColorValid(false),
    featurePickResolved(false),
    featurePrimitiveMetadata()
{
}

static void *
qg_obol_pick_view_context(QgView *display)
{
    return display ? qg_legacy_view_to_context(display->view()) : NULL;
}

static std::string
qg_obol_pick_view_scope_name(QgView *display)
{
    void *view_ctx = qg_obol_pick_view_context(display);
    if (!view_ctx)
	return std::string("shared");

    const char *name = rt_view_context_name_get(view_ctx);
    if (name && name[0])
	return std::string(name);

    char fallback[64] = {0};
    std::snprintf(fallback, sizeof(fallback), "%p", view_ctx);
    return std::string(fallback);
}

static BRLObolFeatureOwner
qg_obol_pick_view_owner(QgView *display)
{
    BRLObolFeatureOwner owner;
    owner.ownerToken = qg_obol_pick_view_context(display);
    owner.ownerId = qg_obol_pick_view_scope_name(display).c_str();
    owner.ownerRole = "view";
    return owner;
}

static int
qg_obol_display_extent(QgView *display, int &width, int &height)
{
    width = 0;
    height = 0;
    if (!display)
	return 0;

    const double dpr = display->devicePixelRatioF();
    width = std::max(1, static_cast<int>(std::ceil(display->width() * dpr)));
    height = std::max(1, static_cast<int>(std::ceil(display->height() * dpr)));
    return width > 0 && height > 0;
}

static const SoBRLPickDetail *
qg_obol_brl_pick_detail(const SoPickedPoint *pickedPoint)
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

static bool
qg_obol_feature_pick_from_record(BRLObolViewController *controller,
	const BRLObolFeatureOwner *owner,
	const QgObolPickRecord &record,
	BRLObolFeaturePrimitivePick &pick)
{
    if (!controller || record.primitiveIndex < 0)
	return false;

    const std::string pickedName = record.path.empty() ?
	record.sourceName : record.path;
    if (pickedName.empty())
	return false;

    if (!(owner && controller->features().resolvePrimitivePick(
	    pickedName.c_str(), static_cast<int32_t>(record.primitiveIndex),
	    pick, BRLOBOL_FEATURE_SCOPE_LOCAL, owner)) &&
	    !controller->features().resolvePrimitivePick(pickedName.c_str(),
		static_cast<int32_t>(record.primitiveIndex), pick,
		BRLOBOL_FEATURE_SCOPE_SHARED, NULL))
	return false;

    return true;
}

static void
qg_obol_resolve_feature_pick(BRLObolViewController *controller,
	const BRLObolFeatureOwner *owner,
	QgObolPickRecord &record)
{
    record.featurePickResolved = false;
    record.featureName.clear();
    record.featurePrimitiveIndex = -1;
    record.featurePrimitiveMetadata.clear();

    BRLObolFeaturePrimitivePick pick;
    if (!qg_obol_feature_pick_from_record(controller, owner, record, pick))
	return;

    record.featurePickResolved = true;
    record.featureName = pick.featureName.getString();
    record.featurePrimitiveIndex = static_cast<int>(pick.primitiveIndex);
    for (size_t i = 0; i < pick.metadata.size(); i++) {
	record.featurePrimitiveMetadata.push_back(std::make_pair(
		std::string(pick.metadata[i].key.getString()),
		std::string(pick.metadata[i].value.getString())));
    }
}

static QgObolPickRecord
qg_obol_pick_record(const SoBRLPickDetail *detail,
	const SbVec3f &point,
	float distance)
{
    QgObolPickRecord record;
    record.point = point;
    record.distance = distance;
    if (!detail)
	return record;

    record.path = detail->getPath().getString();
    record.sourceName = detail->getSourceName().getString();
    record.sourceType = detail->getSourceType().getString();
    record.materialShader = detail->getMaterialShader().getString();
    record.editIntentId = detail->getEditIntentId().getString();
    record.editIntentRole = detail->getEditIntentRole().getString();
    record.sourceId = detail->getSourceId();
    record.regionId = detail->getRegionId();
    record.airCode = detail->getAirCode();
    record.materialId = detail->getMaterialId();
    record.los = detail->getLos();
    record.primitiveKind = static_cast<int>(detail->getPrimitiveKind());
    record.primitiveIndex = detail->getPrimitiveIndex();
    record.faceVertexIndexA = detail->getFaceVertexIndexA();
    record.faceVertexIndexB = detail->getFaceVertexIndexB();
    record.faceVertexIndexC = detail->getFaceVertexIndexC();
    record.nearestFaceEdgeSlot = detail->getNearestFaceEdgeSlot();
    record.nearestFaceEdgeVertexIndexA =
	detail->getNearestFaceEdgeVertexIndexA();
    record.nearestFaceEdgeVertexIndexB =
	detail->getNearestFaceEdgeVertexIndexB();
    record.nearestFaceVertexSlot = detail->getNearestFaceVertexSlot();
    record.nearestFaceVertexIndex = detail->getNearestFaceVertexIndex();
    record.modelPoint = detail->getModelPoint();
    record.materialColorValid = detail->hasMaterialColor() ? true : false;
    record.materialColor = detail->getMaterialColor();
    return record;
}

static QgObolPickRecord
qg_obol_pick_record(const SoPickedPoint *pickedPoint,
	const SoBRLPickDetail *detail,
	const SbVec3f *rayOrigin)
{
    SbVec3f point(0.0f, 0.0f, 0.0f);
    float distance = FLT_MAX;

    if (pickedPoint) {
	point = pickedPoint->getPoint();
	if (rayOrigin)
	    distance = (point - *rayOrigin).length();
    }

    return qg_obol_pick_record(detail, point, distance);
}

static QgObolPickRecord
qg_obol_pick_record(const BRLObolSourceMeshPickResult &pick)
{
    return qg_obol_pick_record(&pick.detail, pick.point, pick.distance);
}

static QgObolPickRecord
qg_obol_pick_record(const BRLObolRtPickResult &pick)
{
    return qg_obol_pick_record(&pick.detail, pick.point, pick.distance);
}

static int
qg_obol_apply_feature_states(BRLObolViewController *controller,
	const BRLObolFeatureOwner *owner,
	const std::vector<QgObolPickRecord> &records,
	bool select,
	bool highlight)
{
    if (!controller || records.empty() || (!select && !highlight))
	return 0;

    struct FeaturePrimitiveGroup {
	BRLObolFeatureHandle handle;
	std::vector<int32_t> primitives;
    };
    std::vector<FeaturePrimitiveGroup> groups;

    for (size_t i = 0; i < records.size(); i++) {
	BRLObolFeaturePrimitivePick pick;
	if (!qg_obol_feature_pick_from_record(controller, owner, records[i],
		pick))
	    continue;

	FeaturePrimitiveGroup *group = NULL;
	for (size_t j = 0; j < groups.size(); j++) {
	    if (groups[j].handle.id == pick.handle.id) {
		group = &groups[j];
		break;
	    }
	}
	if (!group) {
	    FeaturePrimitiveGroup next;
	    next.handle = pick.handle;
	    groups.push_back(next);
	    group = &groups.back();
	}
	if (std::find(group->primitives.begin(), group->primitives.end(),
		pick.primitiveIndex) == group->primitives.end())
	    group->primitives.push_back(pick.primitiveIndex);
    }

    int applied = 0;
    for (size_t i = 0; i < groups.size(); i++) {
	if (groups[i].primitives.empty())
	    continue;
	if (select && !controller->features().replaceSelectedPrimitives(
		groups[i].handle, groups[i].primitives))
	    continue;
	if (highlight && !controller->features().replaceHighlightedPrimitives(
		groups[i].handle, groups[i].primitives))
	    continue;
	applied++;
    }

    return applied;
}

static std::string
qg_obol_normalized_pick_path(const std::string &path)
{
    size_t start = 0;
    while (start < path.size() && path[start] == '/')
	start++;
    return path.substr(start);
}

static std::string
qg_obol_pick_record_path(const QgObolPickRecord &record)
{
    if (record.featurePickResolved)
	return qg_obol_normalized_pick_path(record.featureName);
    return qg_obol_normalized_pick_path(record.path);
}

static bool
qg_obol_pick_path_already_recorded(
	const std::vector<QgObolPickRecord> &records,
	const QgObolPickRecord &candidate)
{
    const std::string candidatePath = qg_obol_pick_record_path(candidate);
    if (candidatePath.empty())
	return false;

    for (const QgObolPickRecord &record : records) {
	if (qg_obol_pick_record_path(record) == candidatePath)
	    return true;
    }

    return false;
}

static void
qg_obol_insert_pick_record(std::vector<QgObolPickRecord> &records,
	const QgObolPickRecord &record,
	bool pickAll)
{
    if (pickAll) {
	records.push_back(record);
	return;
    }

    if (records.empty() || record.distance < records[0].distance)
	records.assign(1, record);
}

static bool
qg_obol_pick_record_nearer(const QgObolPickRecord &a,
	const QgObolPickRecord &b)
{
    return a.distance < b.distance;
}

static void
qg_obol_sort_pick_records(std::vector<QgObolPickRecord> &records)
{
    if (records.size() > 1)
	std::stable_sort(records.begin(), records.end(),
		qg_obol_pick_record_nearer);
}

static int
qg_obol_pick_source_full_detail(BRLObolViewController *controller,
	const BRLObolFeatureOwner *owner,
	const SbLine &line,
	bool pickAll,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    if (!controller)
	return 0;

    BRLObolSourceMeshPickResult sourcePick;
    int submitted = 0;
    int added = controller->pickSourceMeshExactRay(sourcePick,
	    line.getPosition(), line.getDirection(), 0, &submitted) > 0 &&
	sourcePick.hit ? 1 : 0;
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = submitted;
    if (added) {
	QgObolPickRecord record = qg_obol_pick_record(sourcePick);
	qg_obol_resolve_feature_pick(controller, owner, record);
	qg_obol_insert_pick_record(records, record, pickAll);
    }

    return added;
}

static SbBool
qg_obol_pick_camera_line(BRLObolViewController *controller,
	int vx,
	int vy,
	SbLine &line)
{
    if (!controller || !controller->getCamera())
	return FALSE;

    const SbViewportRegion &region = controller->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return FALSE;

    const float nx = size[0] > 0 ?
	std::max(0.0f, std::min(1.0f,
	    static_cast<float>(vx) / static_cast<float>(size[0]))) : 0.5f;
    const float ny = size[1] > 0 ?
	std::max(0.0f, std::min(1.0f,
	    static_cast<float>(vy + 1) / static_cast<float>(size[1]))) : 0.5f;
    const float aspect = static_cast<float>(size[0]) /
	static_cast<float>(size[1]);

    SbViewVolume viewVolume = controller->getCamera()->getViewVolume(aspect);
    viewVolume.projectPointToLine(SbVec2f(nx, ny), line);
    return TRUE;
}

static int
qg_obol_pick_rt_exact(BRLObolViewController *controller,
	const BRLObolFeatureOwner *owner,
	const SbLine &line,
	bool pickAll,
	std::vector<QgObolPickRecord> &records)
{
    if (!controller)
	return 0;

    int added = 0;
    std::vector<BRLObolRtPickResult> rtPicks;
    controller->pickRtExactRay(rtPicks, line.getPosition(),
	    line.getDirection(), pickAll ? TRUE : FALSE);
    for (size_t i = 0; i < rtPicks.size(); i++) {
	const BRLObolRtPickResult &rtPick = rtPicks[i];
	QgObolPickRecord record = qg_obol_pick_record(rtPick);
	qg_obol_resolve_feature_pick(controller, owner, record);
	if (qg_obol_pick_path_already_recorded(records, record))
	    continue;

	qg_obol_insert_pick_record(records, record, pickAll);
	added++;
    }

    return added;
}

static int
qg_obol_pick_point_internal(QgView *display,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    records.clear();
    if (!display)
	return 0;

    BRLObolViewController *controller = display->obolViewController();
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;
    BRLObolFeatureOwner owner = qg_obol_pick_view_owner(display);

    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
		static_cast<unsigned int>(height));

    const SbViewportRegion &region = controller->getViewportRegion();
    SbVec2s size = region.getViewportSizePixels();
    if (size[0] <= 0 || size[1] <= 0)
	return 0;

    int vx = x;
    int vy = size[1] - 1 - y;
    if (vx < 0)
	vx = 0;
    if (vy < 0)
	vy = 0;
    if (vx >= size[0])
	vx = size[0] - 1;
    if (vy >= size[1])
	vy = size[1] - 1;

    SoRayPickAction pickAction(region);
    pickAction.setPoint(SbVec2s(static_cast<short>(vx),
	    static_cast<short>(vy)));
    pickAction.setRadius(radiusPixels > 0.0f ? radiusPixels : 1.0f);
    pickAction.setPickAll(pickAll ? TRUE : FALSE);
    pickAction.apply(controller->getViewport()->getRoot());

    const SbLine &line = pickAction.getLine();
    const SbVec3f &rayOrigin = line.getPosition();

    if (pickAll) {
	const SoPickedPointList &pickedPoints =
	    pickAction.getPickedPointList();
	if (pickedPoints.getLength() > 0)
	    records.reserve(records.size() + pickedPoints.getLength());
	for (int i = 0; i < pickedPoints.getLength(); i++) {
	    const SoPickedPoint *pickedPoint = pickedPoints[i];
	    const SoBRLPickDetail *detail =
		qg_obol_brl_pick_detail(pickedPoint);
	    if (detail) {
		QgObolPickRecord record = qg_obol_pick_record(pickedPoint,
			detail, &rayOrigin);
		qg_obol_resolve_feature_pick(controller, &owner, record);
		qg_obol_insert_pick_record(records, record, pickAll);
	    }
	}
    } else {
	const SoPickedPoint *pickedPoint = pickAction.getPickedPoint();
	const SoBRLPickDetail *detail =
	    qg_obol_brl_pick_detail(pickedPoint);
	if (detail) {
	    QgObolPickRecord record = qg_obol_pick_record(pickedPoint,
		    detail, &rayOrigin);
	    qg_obol_resolve_feature_pick(controller, &owner, record);
	    qg_obol_insert_pick_record(records, record, pickAll);
	}
    }

    qg_obol_pick_source_full_detail(controller, &owner, line, pickAll, records,
	    submittedSourceRequestCount);

    {
	SbLine rtLine = line;
	if (records.empty())
	    (void)qg_obol_pick_camera_line(controller, vx, vy, rtLine);
	qg_obol_pick_rt_exact(controller, &owner, rtLine, pickAll, records);
    }

    if (pickAll)
	qg_obol_sort_pick_records(records);

    return static_cast<int>(records.size());
}

int
qg_obol_pick_point(QgView *display,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount)
{
    return qg_obol_pick_point_internal(display, x, y, radiusPixels, pickAll,
	    records, submittedSourceRequestCount);
}

int
qg_obol_pick_ray(QgView *display,
	const SbVec3f &rayOrigin,
	const SbVec3f &rayDirection,
	bool pickAll,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    records.clear();
    if (!display)
	return 0;

    BRLObolViewController *controller = display->obolViewController();
    if (!controller || !controller->getViewport() ||
	    !controller->getViewport()->getRoot())
	return 0;
    BRLObolFeatureOwner owner = qg_obol_pick_view_owner(display);

    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
		static_cast<unsigned int>(height));

    const SbViewportRegion &region = controller->getViewportRegion();
    SoRayPickAction pickAction(region);
    pickAction.setRay(rayOrigin, rayDirection);
    pickAction.setPickAll(pickAll ? TRUE : FALSE);
    pickAction.apply(controller->getViewport()->getRoot());

    SbLine line;
    line.setPosDir(rayOrigin, rayDirection);
    const SbVec3f &origin = line.getPosition();

    if (pickAll) {
	const SoPickedPointList &pickedPoints =
	    pickAction.getPickedPointList();
	if (pickedPoints.getLength() > 0)
	    records.reserve(records.size() + pickedPoints.getLength());
	for (int i = 0; i < pickedPoints.getLength(); i++) {
	    const SoPickedPoint *pickedPoint = pickedPoints[i];
	    const SoBRLPickDetail *detail =
		qg_obol_brl_pick_detail(pickedPoint);
	    if (detail) {
		QgObolPickRecord record = qg_obol_pick_record(pickedPoint,
			detail, &origin);
		qg_obol_resolve_feature_pick(controller, &owner, record);
		qg_obol_insert_pick_record(records, record, pickAll);
	    }
	}
    } else {
	const SoPickedPoint *pickedPoint = pickAction.getPickedPoint();
	const SoBRLPickDetail *detail =
	    qg_obol_brl_pick_detail(pickedPoint);
	if (detail) {
	    QgObolPickRecord record = qg_obol_pick_record(pickedPoint,
		    detail, &origin);
	    qg_obol_resolve_feature_pick(controller, &owner, record);
	    qg_obol_insert_pick_record(records, record, pickAll);
	}
    }

    qg_obol_pick_source_full_detail(controller, &owner, line, pickAll, records,
	    submittedSourceRequestCount);
    qg_obol_pick_rt_exact(controller, &owner, line, pickAll, records);

    if (pickAll)
	qg_obol_sort_pick_records(records);

    return static_cast<int>(records.size());
}

static std::string
qg_obol_pick_unique_key(const QgObolPickRecord &record)
{
    char buffer[128] = {0};
    if (record.featurePickResolved) {
	std::snprintf(buffer, sizeof(buffer), ":feature:%d:%d",
		record.primitiveKind,
		record.featurePrimitiveIndex);
	std::string key = record.featureName;
	key.append(buffer);
	return key;
    }
    std::snprintf(buffer, sizeof(buffer), ":%u:%d:%d:%d",
	    record.sourceId,
	    record.primitiveKind,
	    record.primitiveIndex,
	    record.materialId);
    std::string key = record.path.empty() ? record.sourceName : record.path;
    key.append(buffer);
    return key;
}

int
qg_obol_pick_rect(QgView *display,
	int x0,
	int y0,
	int x1,
	int y1,
	float radiusPixels,
	bool firstOnly,
	std::vector<QgObolPickRecord> &records,
	int *submittedSourceRequestCount)
{
    if (submittedSourceRequestCount)
	*submittedSourceRequestCount = 0;
    records.clear();
    if (!display)
	return 0;

    int minX = std::min(x0, x1);
    int maxX = std::max(x0, x1);
    int minY = std::min(y0, y1);
    int maxY = std::max(y0, y1);

    int displayWidth = 0;
    int displayHeight = 0;
    if (qg_obol_display_extent(display, displayWidth, displayHeight)) {
	minX = std::max(0, std::min(minX, displayWidth - 1));
	maxX = std::max(0, std::min(maxX, displayWidth - 1));
	minY = std::max(0, std::min(minY, displayHeight - 1));
	maxY = std::max(0, std::min(maxY, displayHeight - 1));
    }

    int width = std::max(1, maxX - minX);
    int height = std::max(1, maxY - minY);
    int xSteps = std::max(1, std::min(6, width / 16));
    int ySteps = std::max(1, std::min(6, height / 16));

    std::unordered_set<std::string> seen;
    for (int yi = 0; yi <= ySteps; yi++) {
	int y = minY + (height * yi) / ySteps;
	for (int xi = 0; xi <= xSteps; xi++) {
	    int x = minX + (width * xi) / xSteps;
	    std::vector<QgObolPickRecord> sampledRecords;
	    int submitted = 0;
	    qg_obol_pick_point_internal(display, x, y, radiusPixels,
		    !firstOnly, sampledRecords, &submitted);
	    if (submittedSourceRequestCount)
		*submittedSourceRequestCount += submitted;
	    for (const QgObolPickRecord &record : sampledRecords) {
		std::string key = qg_obol_pick_unique_key(record);
		if (!seen.insert(key).second)
		    continue;
		records.push_back(record);
		if (firstOnly)
		    return static_cast<int>(records.size());
	    }
	}
    }

    return static_cast<int>(records.size());
}

int
qg_obol_pick_apply_feature_state(QgView *display,
	const QgObolPickRecord &record,
	bool select,
	bool highlight)
{
    std::vector<QgObolPickRecord> records;
    records.push_back(record);
    return qg_obol_pick_apply_feature_states(display, records, select,
	    highlight);
}

int
qg_obol_pick_apply_feature_states(QgView *display,
	const std::vector<QgObolPickRecord> &records,
	bool select,
	bool highlight)
{
    if (!display)
	return 0;

    BRLObolViewController *controller = display->obolViewController();
    BRLObolFeatureOwner owner = qg_obol_pick_view_owner(display);
    return qg_obol_apply_feature_states(controller, &owner, records, select,
	    highlight);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
