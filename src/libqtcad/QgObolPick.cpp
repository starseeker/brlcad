/*                  Q G O B O L P I C K . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file QgObolPick.cpp */

#include "common.h"

#include "qtcad/QgObolPick.h"

#include "BObol/BViewController.h"
#include "BObol/BViewQuery.h"
#include "bu/vls.h"
#include "ged/draw.h"
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
    featureMetadata(),
    featurePrimitiveMetadata()
{
}

static struct ged_view_context *
qg_obol_pick_view_context(QgView *display)
{
    return display ? ged_view_context_from_bv(display->viewContext()) : NULL;
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
qg_obol_feature_pick_from_record(struct ged_view_context *view_ctx,
	const QgObolPickRecord &record,
	std::string &feature_name,
	int &feature_primitive)
{
    feature_name.clear();
    feature_primitive = -1;
    if (!view_ctx || record.primitiveIndex < 0)
	return false;

    const std::string pickedName = record.path.empty() ?
	record.sourceName : record.path;
    if (pickedName.empty())
	return false;

    struct bu_vls resolved_name = BU_VLS_INIT_ZERO;
    const int resolved = ged_view_feature_pick_primitive_resolve(
	view_ctx, pickedName.c_str(), record.primitiveIndex, 0, 0,
	&resolved_name, &feature_primitive);
    if (resolved)
	feature_name = bu_vls_cstr(&resolved_name);
    bu_vls_free(&resolved_name);
    return resolved ? true : false;
}

static void
qg_obol_resolve_feature_pick(struct ged_view_context *view_ctx,
	QgObolPickRecord &record)
{
    record.featurePickResolved = false;
    record.featureName.clear();
    record.featurePrimitiveIndex = -1;
    record.featureMetadata.clear();
    record.featurePrimitiveMetadata.clear();

    if (!qg_obol_feature_pick_from_record(view_ctx, record,
	record.featureName, record.featurePrimitiveIndex))
	return;

    record.featurePickResolved = true;
    const size_t feature_metadata_count =
	ged_view_feature_metadata_count(view_ctx,
	    record.featureName.c_str());
    for (size_t i = 0; i < feature_metadata_count; i++) {
	struct bu_vls key = BU_VLS_INIT_ZERO;
	struct bu_vls value = BU_VLS_INIT_ZERO;
	if (ged_view_feature_metadata_copy(view_ctx,
	    record.featureName.c_str(), i, &key, &value)) {
	    record.featureMetadata.push_back(std::make_pair(
		std::string(bu_vls_cstr(&key)),
		std::string(bu_vls_cstr(&value))));
	}
	bu_vls_free(&key);
	bu_vls_free(&value);
    }
    const size_t metadata_count =
	ged_view_feature_primitive_metadata_count(view_ctx,
	    record.featureName.c_str(), record.featurePrimitiveIndex);
    for (size_t i = 0; i < metadata_count; i++) {
	struct bu_vls key = BU_VLS_INIT_ZERO;
	struct bu_vls value = BU_VLS_INIT_ZERO;
	if (ged_view_feature_primitive_metadata_copy(view_ctx,
	    record.featureName.c_str(), record.featurePrimitiveIndex, i,
	    &key, &value)) {
	    record.featurePrimitiveMetadata.push_back(std::make_pair(
		std::string(bu_vls_cstr(&key)),
		std::string(bu_vls_cstr(&value))));
	}
	bu_vls_free(&key);
	bu_vls_free(&value);
    }
}

static QgObolPickRecord
qg_obol_pick_record(const BObolViewPickRecord &hit)
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
qg_obol_apply_feature_states(struct ged_view_context *view_ctx,
	const std::vector<QgObolPickRecord> &records,
	bool select,
	bool highlight)
{
    if (!view_ctx || records.empty() || (!select && !highlight))
	return 0;

    struct FeaturePrimitiveGroup {
	std::string name;
	std::vector<int> primitives;
    };
    std::vector<FeaturePrimitiveGroup> groups;

    for (const QgObolPickRecord &record : records) {
	std::string feature_name;
	int feature_primitive = -1;
	if (!qg_obol_feature_pick_from_record(view_ctx, record, feature_name,
		feature_primitive))
	    continue;

	FeaturePrimitiveGroup *group = NULL;
	for (FeaturePrimitiveGroup &candidate : groups) {
	    if (candidate.name == feature_name) {
		group = &candidate;
		break;
	    }
	}
	if (!group) {
	    groups.push_back(FeaturePrimitiveGroup());
	    group = &groups.back();
	    group->name = feature_name;
	}
	if (std::find(group->primitives.begin(), group->primitives.end(),
		feature_primitive) == group->primitives.end())
	    group->primitives.push_back(feature_primitive);
    }

    int applied = 0;
    for (const FeaturePrimitiveGroup &group : groups) {
	if (group.primitives.empty())
	    continue;
	if (select && !ged_view_feature_set_selection(
		view_ctx, group.name.c_str(), group.primitives.data(),
		group.primitives.size()))
	    continue;
	if (highlight &&
	    !ged_view_feature_set_highlights(
		view_ctx, group.name.c_str(), group.primitives.data(),
		group.primitives.size()))
	    continue;
	applied++;
    }
    return applied;
}

static void
qg_obol_pick_records(struct ged_view_context *view_ctx,
	const std::vector<BObolViewPickRecord> &hits,
	std::vector<QgObolPickRecord> &records)
{
    records.clear();
    records.reserve(hits.size());
    for (const BObolViewPickRecord &hit : hits) {
	QgObolPickRecord record = qg_obol_pick_record(hit);
	qg_obol_resolve_feature_pick(view_ctx, record);
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

    BObolViewController *controller = display->obolViewController();
    if (!controller)
	return 0;
    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
	    static_cast<unsigned int>(height));

    std::vector<BObolViewPickRecord> hits;
    const int count = bobol_view_pick_point(controller, x, y,
	radiusPixels, pickAll, hits, submittedSourceRequestCount);
    qg_obol_pick_records(qg_obol_pick_view_context(display), hits, records);
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
qg_obol_pick_display_point(QgView *display,
	int x,
	int y,
	float radiusPixels,
	bool pickAll,
	std::vector<QgObolPickRecord> &records)
{
    records.clear();
    if (!display)
	return 0;

    BObolViewController *controller = display->obolViewController();
    if (!controller)
	return 0;
    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
	    static_cast<unsigned int>(height));

    std::vector<BObolViewPickRecord> hits;
    const int count = bobol_view_pick_display_point(controller, x, y,
	radiusPixels, pickAll, hits);
    qg_obol_pick_records(qg_obol_pick_view_context(display), hits, records);
    return count;
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

    BObolViewController *controller = display->obolViewController();
    if (!controller)
	return 0;
    int width = 0;
    int height = 0;
    if (qg_obol_display_extent(display, width, height))
	controller->setViewportSize(static_cast<unsigned int>(width),
	    static_cast<unsigned int>(height));

    std::vector<BObolViewPickRecord> hits;
    const int count = bobol_view_pick_ray(controller, rayOrigin,
	rayDirection, pickAll, hits, submittedSourceRequestCount);
    qg_obol_pick_records(qg_obol_pick_view_context(display), hits, records);
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

    BObolViewController *controller = display->obolViewController();
    if (controller) {
	controller->setViewportSize(static_cast<unsigned int>(width),
	    static_cast<unsigned int>(height));
	std::vector<BObolViewPickRecord> bounded;
	const int boundedCount = bobol_view_pick_rectangle(controller,
	    minX, minY, maxX, maxY, bounded);
	std::vector<BObolViewPickRecord> sourceBounded;
	const int sourceBoundedCount = bobol_view_pick_source_rectangle(
	    controller, minX, minY, maxX, maxY, sourceBounded);
	bounded.insert(bounded.end(), sourceBounded.begin(),
	    sourceBounded.end());
	if (boundedCount >= 0 || sourceBoundedCount >= 0) {
	    if (firstOnly && bounded.size() > 1)
		bounded.resize(1);
	    qg_obol_pick_records(qg_obol_pick_view_context(display),
		bounded, records);
	    return static_cast<int>(records.size());
	}
    }

    /* A non-CAD scene has no retained occurrence registry.  Keep a bounded
     * display-only fallback for view features and small standalone scenes;
     * exact source/librt work remains exclusive to the explicit ray tool. */
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
	    qg_obol_pick_display_point(display, x, y, radiusPixels,
		!firstOnly, sampled);
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
    return qg_obol_apply_feature_states(qg_obol_pick_view_context(display),
	records, select, highlight);
}

// Local Variables:
// mode: C++
// tab-width: 8
// c-basic-offset: 4
// indent-tabs-mode: t
// c-file-style: "stroustrup"
// End:
// ex: shiftwidth=4 tabstop=8
