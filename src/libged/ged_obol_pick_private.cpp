/*              G E D _ O B O L _ P I C K _ P R I V A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_obol_pick_private.cpp
 *
 * Obol-backed point, rectangle, and snap picking adapter.
 */

#include "common.h"

#include "ged/display_obol_private.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BADC.h"
#include "BObol/BDrawCache.h"
#include "BObol/BGrid.h"
#include "BObol/BInit.h"
#include "BObol/BImageSource.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BLodRealization.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSceneController.h"
#include "BObol/BSceneGroup.h"
#include "BObol/BSnapAction.h"
#include "BObol/BSourceRealization.h"
#include "BObol/BViewportImage.h"
#include "BObol/BVListShape.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewQuery.h"
#include "BObol/BViewStore.h"
#include "bg/line_layer.h"
#include "bg/plane.h"
#include "bg/polygon.h"
#include "bv.h"
#include "bu/hash.h"
#include "bu/malloc.h"
#include "bu/parallel.h"
#include "bu/path.h"
#include "bu/str.h"
#include "bu/datetime.h"
#include "bu/units.h"
#include "bu/vls.h"
#include "ged.h"
#include "ged/db_index.h"
#include "ged/draw.h"
#include "ged/display.h"
#include "ged/view.h"
#include "icv.h"
#include "rt/db5.h"
#include "rt/db_fullpath.h"
#include "rt/search.h"
#include "rt/view.h"
#include "vmath.h"

#include "./ged_bobol_private.hpp"
#include "./draw_obol_bridge_private.hpp"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

#include <algorithm>
#include <Inventor/SbVec2f.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>
#include <Inventor/nodes/SoSeparator.h>
#include <Obol/cad/SoCADAssembly.h>
#include <float.h>
#include <inttypes.h>
#include <limits.h>
#include <math.h>
#include <stdint.h>
#include <stdlib.h>
#include <string.h>
#include <atomic>
#include <memory>
#include <set>
#include <string>
#include <unordered_map>
#include <unordered_set>
#include <utility>
#include <vector>

static const struct bv *
ged_obol_bv_const(const struct ged_view_context *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
}


static BObolViewController *
ged_obol_view_controller_for_context(struct ged_view_context *view_ctx)
{
    return ged_bobol_view_controller(view_ctx);
}


struct ged_obol_pick_candidate {
    ged_obol_pick_candidate(void) :
	distance(FLT_MAX),
	source_id(0),
	primitive_kind(0),
	primitive_index(-1),
	material_id(0),
	nearest_face_vertex_index(-1),
	model_point(0.0f, 0.0f, 0.0f)
    {
	face_vertex_index[0] = -1;
	face_vertex_index[1] = -1;
	face_vertex_index[2] = -1;
    }

    std::string path;
    float distance;
    uint32_t source_id;
    int primitive_kind;
    int primitive_index;
    int material_id;
    int face_vertex_index[3];
    int nearest_face_vertex_index;
    SbVec3f model_point;
};

static std::string
ged_obol_normalized_pick_path(const std::string &path)
{
    size_t start = 0;
    while (start < path.size() && path[start] == '/')
	start++;
    return path.substr(start);
}

static std::string
ged_obol_pick_candidate_key(const ged_obol_pick_candidate &candidate)
{
    char buffer[128] = {0};
    snprintf(buffer, sizeof(buffer), ":%u:%d:%d:%d",
	    candidate.source_id, candidate.primitive_kind,
	    candidate.primitive_index, candidate.material_id);
    return ged_obol_normalized_pick_path(candidate.path) + buffer;
}

static BObolFeatureOwner
ged_obol_pick_view_owner(struct ged_view_context *view_ctx)
{
    BObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    owner.ownerRole = "view";

    const struct bv *view = ged_obol_bv_const(view_ctx);
    const char *name = bv_name_get(view);
    if (name && name[0]) {
	owner.ownerId = name;
	return owner;
    }

    char fallback[64] = {0};
    snprintf(fallback, sizeof(fallback), "%p", (void *)view_ctx);
    owner.ownerId = fallback;
    return owner;
}

static std::string
ged_obol_pick_resolved_path(BObolViewController *controller,
			    const BObolFeatureOwner *owner,
			    const SoBRLPickDetail *detail)
{
    if (!detail)
	return std::string();

    const char *path = detail->getPath().getString();
    const char *source_name = detail->getSourceName().getString();
    const std::string picked_name =
	(path && path[0]) ? std::string(path) :
	(source_name ? std::string(source_name) : std::string());
    if (!controller || picked_name.empty() || detail->getPrimitiveIndex() < 0)
	return picked_name;

    BObolFeaturePrimitivePick pick;
    if ((owner && controller->features().resolvePrimitivePick(
	    picked_name.c_str(), detail->getPrimitiveIndex(), pick,
	    BOBOL_FEATURE_SCOPE_LOCAL, owner)) ||
	    controller->features().resolvePrimitivePick(picked_name.c_str(),
		detail->getPrimitiveIndex(), pick,
		BOBOL_FEATURE_SCOPE_SHARED, NULL))
	return std::string(pick.featureName.getString());

    return picked_name;
}

static ged_obol_pick_candidate
ged_obol_pick_candidate_from_detail(BObolViewController *controller,
				    const BObolFeatureOwner *owner,
				    const SoBRLPickDetail *detail,
				    const SbVec3f &point,
				    float distance)
{
    ged_obol_pick_candidate candidate;
    candidate.path = ged_obol_pick_resolved_path(controller, owner, detail);
    candidate.distance = distance;
    candidate.source_id = detail ? detail->getSourceId() : 0;
    candidate.primitive_kind = detail ?
	static_cast<int>(detail->getPrimitiveKind()) : 0;
    candidate.primitive_index = detail ? detail->getPrimitiveIndex() : -1;
    candidate.material_id = detail ? detail->getMaterialId() : 0;
    if (detail) {
	candidate.face_vertex_index[0] = detail->getFaceVertexIndexA();
	candidate.face_vertex_index[1] = detail->getFaceVertexIndexB();
	candidate.face_vertex_index[2] = detail->getFaceVertexIndexC();
	candidate.nearest_face_vertex_index =
	    detail->getNearestFaceVertexIndex();
	candidate.model_point = detail->getModelPoint();
    } else {
	candidate.model_point = point;
    }
    return candidate;
}

static ged_obol_pick_candidate
ged_obol_pick_candidate_from_record(BObolViewController *controller,
				    const BObolFeatureOwner *owner,
				    const BObolViewPickRecord &record)
{
    return ged_obol_pick_candidate_from_detail(controller, owner,
	&record.detail, record.point, record.distance);
}

static void
ged_obol_pick_insert(std::vector<ged_obol_pick_candidate> &candidates,
		     const ged_obol_pick_candidate &candidate,
		     bool pick_all)
{
    if (candidate.path.empty())
	return;

    if (pick_all) {
	candidates.push_back(candidate);
	return;
    }

    if (candidates.empty() || candidate.distance < candidates[0].distance)
	candidates.assign(1, candidate);
}

static bool
ged_obol_pick_candidate_nearer(const ged_obol_pick_candidate &a,
			       const ged_obol_pick_candidate &b)
{
    return a.distance < b.distance;
}

static void
ged_obol_pick_sort(std::vector<ged_obol_pick_candidate> &candidates)
{
    if (candidates.size() > 1)
	std::stable_sort(candidates.begin(), candidates.end(),
		ged_obol_pick_candidate_nearer);
}

static struct ged_pick_result *
ged_obol_pick_result_from_candidates(
	const std::vector<ged_obol_pick_candidate> &candidates)
{
    struct ged_pick_result *result = ged_pick_result_create();
    if (!result)
	return NULL;

    for (size_t i = 0; i < candidates.size(); i++) {
	struct ged_pick_detail detail = GED_PICK_DETAIL_INIT;
	detail.source_id = candidates[i].source_id;
	detail.primitive_kind = candidates[i].primitive_kind;
	detail.primitive_index = candidates[i].primitive_index;
	detail.material_id = candidates[i].material_id;
	detail.face_vertex_index[0] = candidates[i].face_vertex_index[0];
	detail.face_vertex_index[1] = candidates[i].face_vertex_index[1];
	detail.face_vertex_index[2] = candidates[i].face_vertex_index[2];
	detail.nearest_face_vertex_index =
	    candidates[i].nearest_face_vertex_index;
	detail.model_point[X] = candidates[i].model_point[0];
	detail.model_point[Y] = candidates[i].model_point[1];
	detail.model_point[Z] = candidates[i].model_point[2];
	detail.model_point_valid = 1;
	(void)ged_pick_result_append_detail(result,
		candidates[i].path.c_str(), candidates[i].distance, &detail);
    }

    return result;
}

static int
ged_obol_pick_sync_view_controller(BObolViewController *controller,
				   struct ged_view_context *view_ctx)
{
    if (!controller || !view_ctx)
	return 0;

    const struct bv *view = ged_obol_bv_const(view_ctx);
    const int width = bv_width_get(view);
    const int height = bv_height_get(view);
    if (width > 0 && height > 0)
	controller->setViewportSize(static_cast<unsigned int>(width),
		static_cast<unsigned int>(height));

    controller->syncCameraFromViewContext(view_ctx, TRUE);
    return 1;
}

static int
ged_obol_pick_point_candidates(
	struct ged_view_context *view_ctx,
	int x,
	int y,
	float radius_pixels,
	bool pick_all,
	std::vector<ged_obol_pick_candidate> &candidates)
{
    candidates.clear();
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    ged_obol_pick_sync_view_controller(controller, view_ctx);

    std::vector<BObolViewPickRecord> records;
    (void)bobol_view_pick_point(controller, x, y, radius_pixels,
	pick_all, records, NULL);
    BObolFeatureOwner owner = ged_obol_pick_view_owner(view_ctx);
    candidates.reserve(records.size());
    for (const BObolViewPickRecord &record : records)
	ged_obol_pick_insert(candidates,
	    ged_obol_pick_candidate_from_record(controller, &owner, record),
	    pick_all);

    if (pick_all)
	ged_obol_pick_sort(candidates);

    return static_cast<int>(candidates.size());
}

extern "C" struct ged_pick_result *
ged_draw_obol_view_context_pick_point(struct ged_view_context *view_ctx,
				      int x,
				      int y,
				      int first_only)
{
    std::vector<ged_obol_pick_candidate> candidates;
    (void)ged_obol_pick_point_candidates(view_ctx, x, y, 1.0f,
	    first_only ? false : true, candidates);
    return ged_obol_pick_result_from_candidates(candidates);
}

extern "C" struct ged_pick_result *
ged_draw_obol_view_context_pick_nearest(struct ged_view_context *view_ctx, int x, int y)
{
    return ged_draw_obol_view_context_pick_point(view_ctx, x, y, 1);
}

extern "C" struct ged_pick_result *
ged_draw_obol_view_context_pick_rect(struct ged_view_context *view_ctx,
				     int x0,
				     int y0,
				     int x1,
				     int y1)
{
    std::vector<ged_obol_pick_candidate> candidates;
    std::unordered_set<std::string> seen;
    int min_x = std::min(x0, x1);
    int max_x = std::max(x0, x1);
    int min_y = std::min(y0, y1);
    int max_y = std::max(y0, y1);

    const struct bv *view = ged_obol_bv_const(view_ctx);
    const int width_px = bv_width_get(view);
    const int height_px = bv_height_get(view);
    if (width_px > 0) {
	min_x = std::max(0, std::min(min_x, width_px - 1));
	max_x = std::max(0, std::min(max_x, width_px - 1));
    }
    if (height_px > 0) {
	min_y = std::max(0, std::min(min_y, height_px - 1));
	max_y = std::max(0, std::min(max_y, height_px - 1));
    }

    int width = std::max(1, max_x - min_x);
    int height = std::max(1, max_y - min_y);
    int x_steps = std::max(1, std::min(6, width / 16));
    int y_steps = std::max(1, std::min(6, height / 16));

    for (int yi = 0; yi <= y_steps; yi++) {
	int y = min_y + (height * yi) / y_steps;
	for (int xi = 0; xi <= x_steps; xi++) {
	    int x = min_x + (width * xi) / x_steps;
	    std::vector<ged_obol_pick_candidate> sampled;
	    ged_obol_pick_point_candidates(view_ctx, x, y, 1.0f, true,
		    sampled);
	    for (size_t i = 0; i < sampled.size(); i++) {
		const std::string key = ged_obol_pick_candidate_key(sampled[i]);
		if (!seen.insert(key).second)
		    continue;
		candidates.push_back(sampled[i]);
	    }
	}
    }

    ged_obol_pick_sort(candidates);
    return ged_obol_pick_result_from_candidates(candidates);
}

static uint32_t
ged_obol_snap_kind(enum ged_selection_snap_kind kind)
{
    switch (kind) {
	case GED_SELECTION_SNAP_GRID:
	    return static_cast<uint32_t>(SoBRLSnapAction::GRID);
	case GED_SELECTION_SNAP_ENDPOINT:
	    return static_cast<uint32_t>(SoBRLSnapAction::ENDPOINT);
	default:
	    return 0;
    }
}

extern "C" int
ged_draw_obol_view_context_snap_first_candidate(
	struct ged_view_context *view_ctx,
	const point_t sample,
	enum ged_selection_snap_kind kind,
	point_t candidate)
{
    if (candidate)
	VSETALL(candidate, 0.0);
    if (!view_ctx || !sample || !candidate)
	return 0;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    if (!ged_obol_pick_sync_view_controller(controller, view_ctx))
	return 0;

    const uint32_t enabled_kinds = ged_obol_snap_kind(kind);
    if (!enabled_kinds)
	return 0;

    BObolViewSnapRecord snap;
    if (!bobol_view_snap_point(controller,
	SbVec3f(static_cast<float>(sample[X]), static_cast<float>(sample[Y]),
	    static_cast<float>(sample[Z])), FLT_MAX, enabled_kinds,
	SoBRLSnapAction::FULL_DETAIL, false, snap))
	return 0;

    const SbVec3f &point = snap.point;
    VSET(candidate, point[0], point[1], point[2]);
    return 1;
}

/*
 * Local Variables:
 * mode: C++
 * tab-width: 8
 * c-basic-offset: 4
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
