/*                  Q G O B O L P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolPick.cpp */

#include "common.h"

#include "qtcad/QgObolPick.h"

#include "brlobol/view_controller.h"
#include "brlobol/view_query.h"
#include "brlobol/view_store.h"
#include "bv.h"
#include "qtcad/QgView.h"

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
    return display ? display->viewContext() : NULL;
}

static std::string
qg_obol_pick_view_scope_name(QgView *display)
{
    void *view_ctx = qg_obol_pick_view_context(display);
    if (!view_ctx)
	return std::string("shared");

    const struct bv *view = bv_context_view_const(
	reinterpret_cast<const struct bv_context *>(view_ctx));
    const char *name = bv_name_get(view);
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

    return (owner && controller->features().resolvePrimitivePick(
	pickedName.c_str(), static_cast<int32_t>(record.primitiveIndex), pick,
	BRLOBOL_FEATURE_SCOPE_LOCAL, owner)) ||
	controller->features().resolvePrimitivePick(pickedName.c_str(),
	static_cast<int32_t>(record.primitiveIndex), pick,
	BRLOBOL_FEATURE_SCOPE_SHARED, NULL);
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
qg_obol_pick_record(const BRLObolViewPickRecord &hit)
{
    QgObolPickRecord record;
    const SoBRLPickDetail &detail = hit.detail;
    record.path = detail.getPath().getString();
    record.sourceInstanceKey = detail.getSourceInstanceKey().getString();
    record.sourceName = detail.getSourceName().getString();
    record.sourceType = detail.getSourceType().getString();
    record.materialShader = detail.getMaterialShader().getString();
    record.editIntentId = detail.getEditIntentId().getString();
    record.editIntentRole = detail.getEditIntentRole().getString();
    record.point = hit.point;
    record.modelPoint = detail.getModelPoint();
    record.materialColor = detail.getMaterialColor();
    record.distance = hit.distance;
    record.sourceId = detail.getSourceId();
    record.regionId = detail.getRegionId();
    record.airCode = detail.getAirCode();
    record.materialId = detail.getMaterialId();
    record.los = detail.getLos();
    record.primitiveKind = static_cast<int>(detail.getPrimitiveKind());
    record.primitiveIndex = detail.getPrimitiveIndex();
    record.faceVertexIndexA = detail.getFaceVertexIndexA();
    record.faceVertexIndexB = detail.getFaceVertexIndexB();
    record.faceVertexIndexC = detail.getFaceVertexIndexC();
    record.nearestFaceEdgeSlot = detail.getNearestFaceEdgeSlot();
    record.nearestFaceEdgeVertexIndexA =
	detail.getNearestFaceEdgeVertexIndexA();
    record.nearestFaceEdgeVertexIndexB =
	detail.getNearestFaceEdgeVertexIndexB();
    record.nearestFaceVertexSlot = detail.getNearestFaceVertexSlot();
    record.nearestFaceVertexIndex = detail.getNearestFaceVertexIndex();
    record.materialColorValid = detail.hasMaterialColor() ? true : false;
    return record;
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

    for (const QgObolPickRecord &record : records) {
	BRLObolFeaturePrimitivePick pick;
	if (!qg_obol_feature_pick_from_record(controller, owner, record, pick))
	    continue;

	FeaturePrimitiveGroup *group = NULL;
	for (FeaturePrimitiveGroup &candidate : groups) {
	    if (candidate.handle.id == pick.handle.id) {
		group = &candidate;
		break;
	    }
	}
	if (!group) {
	    groups.push_back(FeaturePrimitiveGroup());
	    group = &groups.back();
	    group->handle = pick.handle;
	}
	if (std::find(group->primitives.begin(), group->primitives.end(),
		pick.primitiveIndex) == group->primitives.end())
	    group->primitives.push_back(pick.primitiveIndex);
    }

    int applied = 0;
    for (const FeaturePrimitiveGroup &group : groups) {
	if (group.primitives.empty())
	    continue;
	if (select && !controller->features().replaceSelectedPrimitives(
		group.handle, group.primitives))
	    continue;
	if (highlight && !controller->features().replaceHighlightedPrimitives(
		group.handle, group.primitives))
	    continue;
	applied++;
    }
    return applied;
}

static void
qg_obol_pick_records(BRLObolViewController *controller,
	const BRLObolFeatureOwner *owner,
	const std::vector<BRLObolViewPickRecord> &hits,
	std::vector<QgObolPickRecord> &records)
{
    records.clear();
    records.reserve(hits.size());
    for (const BRLObolViewPickRecord &hit : hits) {
	QgObolPickRecord record = qg_obol_pick_record(hit);
	qg_obol_resolve_feature_pick(controller, owner, record);
	records.push_back(record);
    }
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
    if (!controller)
	return 0;
    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
	    static_cast<unsigned int>(height));

    std::vector<BRLObolViewPickRecord> hits;
    const int count = brlobol_view_pick_point(controller, x, y,
	radiusPixels, pickAll, hits, submittedSourceRequestCount);
    BRLObolFeatureOwner owner = qg_obol_pick_view_owner(display);
    qg_obol_pick_records(controller, &owner, hits, records);
    return count;
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
    if (!controller)
	return 0;
    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
	    static_cast<unsigned int>(height));

    std::vector<BRLObolViewPickRecord> hits;
    const int count = brlobol_view_pick_ray(controller, rayOrigin,
	rayDirection, pickAll, hits, submittedSourceRequestCount);
    BRLObolFeatureOwner owner = qg_obol_pick_view_owner(display);
    qg_obol_pick_records(controller, &owner, hits, records);
    return count;
}

static std::string
qg_obol_pick_unique_key(const QgObolPickRecord &record)
{
    char buffer[128] = {0};
    if (record.featurePickResolved) {
	std::snprintf(buffer, sizeof(buffer), ":feature:%d:%d",
	    record.primitiveKind, record.featurePrimitiveIndex);
	return record.featureName + buffer;
    }
    std::snprintf(buffer, sizeof(buffer), ":%u:%d:%d:%d", record.sourceId,
	record.primitiveKind, record.primitiveIndex, record.materialId);
    return (record.path.empty() ? record.sourceName : record.path) + buffer;
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
    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height)) {
	minX = std::max(0, std::min(minX, width - 1));
	maxX = std::max(0, std::min(maxX, width - 1));
	minY = std::max(0, std::min(minY, height - 1));
	maxY = std::max(0, std::min(maxY, height - 1));
    }

    const int sampleWidth = std::max(1, maxX - minX);
    const int sampleHeight = std::max(1, maxY - minY);
    const int xSteps = std::max(1, std::min(6, sampleWidth / 16));
    const int ySteps = std::max(1, std::min(6, sampleHeight / 16));
    std::unordered_set<std::string> seen;
    for (int yi = 0; yi <= ySteps; yi++) {
	const int y = minY + (sampleHeight * yi) / ySteps;
	for (int xi = 0; xi <= xSteps; xi++) {
	    const int x = minX + (sampleWidth * xi) / xSteps;
	    std::vector<QgObolPickRecord> sampled;
	    int submitted = 0;
	    qg_obol_pick_point_internal(display, x, y, radiusPixels, !firstOnly,
		sampled, &submitted);
	    if (submittedSourceRequestCount)
		*submittedSourceRequestCount += submitted;
	    for (const QgObolPickRecord &record : sampled) {
		if (!seen.insert(qg_obol_pick_unique_key(record)).second)
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
    std::vector<QgObolPickRecord> records(1, record);
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
