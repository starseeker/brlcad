/*          G E D _ O B O L _ P O L Y G O N _ P R I V A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_obol_polygon_private.cpp
 *
 * Obol-backed retained view-polygon adapter.
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


static BObolViewController *
ged_obol_shared_view_controller_ensure_for_context(
    struct ged_view_context *view_ctx, int sync_current_scene)
{
    struct ged *gedp = ged_view_context_owner(view_ctx);
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_draw_obol_controller(gedp);
    if (controller)
	return controller;
    if (!ged_draw_obol_scene_controller_ensure_owned(gedp, sync_current_scene))
	return NULL;
    return ged_draw_obol_controller(gedp);
}


static BObolViewController *
ged_obol_shared_view_controller_for_context(struct ged_view_context *view_ctx)
{
    struct ged *gedp = ged_view_context_owner(view_ctx);
    return gedp ? ged_draw_obol_controller(gedp) : NULL;
}


static unsigned char
ged_obol_rgb_channel(float value)
{
    if (value <= 0.0f)
	return 0;
    if (value >= 1.0f)
	return 255;
    return static_cast<unsigned char>(value * 255.0f + 0.5f);
}


static void
ged_obol_rgb_from_color(const SbColor &color, unsigned char rgb[3])
{
    rgb[0] = ged_obol_rgb_channel(color.getValue()[0]);
    rgb[1] = ged_obol_rgb_channel(color.getValue()[1]);
    rgb[2] = ged_obol_rgb_channel(color.getValue()[2]);
}


static BObolPolygonType
ged_obol_polygon_type_from_ged(int type)
{
    switch (type) {
	case GED_VIEW_POLYGON_CIRCLE:
	    return BObolPolygonType::Circle;
	case GED_VIEW_POLYGON_ELLIPSE:
	    return BObolPolygonType::Ellipse;
	case GED_VIEW_POLYGON_RECTANGLE:
	    return BObolPolygonType::Rectangle;
	case GED_VIEW_POLYGON_SQUARE:
	    return BObolPolygonType::Square;
	default:
	    return BObolPolygonType::General;
    }
}

static BObolPolygonUpdate
ged_obol_polygon_update_from_ged(int op)
{
    switch (op) {
	case GED_VIEW_POLYGON_UPDATE_PROPS_ONLY:
	    return BObolPolygonUpdate::PropsOnly;
	case GED_VIEW_POLYGON_UPDATE_PT_SELECT:
	    return BObolPolygonUpdate::PointSelect;
	case GED_VIEW_POLYGON_UPDATE_PT_SELECT_CLEAR:
	    return BObolPolygonUpdate::PointSelectClear;
	case GED_VIEW_POLYGON_UPDATE_PT_MOVE:
	    return BObolPolygonUpdate::PointMove;
	case GED_VIEW_POLYGON_UPDATE_PT_APPEND:
	    return BObolPolygonUpdate::PointAppend;
	case GED_VIEW_POLYGON_UPDATE_PT_DELETE:
	    return BObolPolygonUpdate::PointDelete;
	default:
	    return BObolPolygonUpdate::Default;
    }
}

static BObolFeatureScope
ged_obol_polygon_scope(int local)
{
    return local ? BObolFeatureScope::Local : BObolFeatureScope::Shared;
}

static BObolViewController *
ged_obol_polygon_controller_from_ged_ref(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    if (!view_ctx || !ref.owner || !ref.id || !ref.generation)
	return NULL;

    const int local = (ref.owner & 1u) ? 1 : 0;
    if (ref.owner != ged_view_context_reference_owner(view_ctx, local))
	return NULL;

    BObolViewController *controller = local ?
	ged_obol_view_controller_for_context(view_ctx) :
	ged_obol_shared_view_controller_for_context(view_ctx);
    if (!controller ||
	controller->polygons().referenceGeneration() != ref.generation)
	return NULL;

    BObolPolygonRecord record;
    if (!controller->polygons().record(BObolPolygonHandle(ref.id, 1),
	    record))
	return NULL;
    return controller;
}

static BObolPolygonHandle
ged_obol_polygon_handle_from_ged_ref(ged_view_polygon_ref ref)
{
    return (ref.id && ref.generation) ? BObolPolygonHandle(ref.id, 1) :
	   BObolPolygonHandle();
}

static ged_view_polygon_ref
ged_obol_ged_view_polygon_ref(struct ged_view_context *view_ctx,
			 BObolViewController *controller,
			 BObolPolygonHandle handle,
			 int local)
{
    ged_view_polygon_ref ged_ref = GED_VIEW_POLYGON_REF_NULL_INIT;
    if (!controller || !handle.isValid())
	return ged_ref;
    ged_ref.owner = ged_view_context_reference_owner(view_ctx, local);
    ged_ref.id = handle.id;
    ged_ref.generation = controller->polygons().referenceGeneration();
    if (!ged_ref.owner)
	return GED_VIEW_POLYGON_REF_NULL;
    return ged_ref;
}

static ged_view_polygon_ref
ged_obol_ged_view_polygon_ref_for_handle(struct ged_view_context *view_ctx,
				      BObolViewController *controller,
				      BObolPolygonHandle handle)
{
    BObolPolygonRecord record;
    if (!controller || !controller->polygons().record(handle, record))
	return GED_VIEW_POLYGON_REF_NULL;
    return ged_obol_ged_view_polygon_ref(view_ctx, controller, handle,
	record.scope == BObolFeatureScope::Local ? 1 : 0);
}

static void
ged_obol_polygon_point_from_ged(SbVec3f &dst, const point_t src)
{
    dst = SbVec3f(static_cast<float>(src[X]),
		  static_cast<float>(src[Y]),
		  static_cast<float>(src[Z]));
}

static void
ged_obol_polygon_project_point(point_t dst, struct ged_view_context *view_ctx, const point_t src,
			       plane_t *view_plane)
{
    VMOVE(dst, src);
    if (!view_plane)
	return;
    HSET(*view_plane, 0.0, 0.0, 1.0, src[Z]);

    if (bv_plane_get(view_plane, ged_obol_bv_const(view_ctx)) != 0)
	return;

    fastf_t fx = 0.0;
    fastf_t fy = 0.0;
    bg_plane_closest_pt(&fx, &fy, view_plane, (point_t *)src);

    point_t projected = VINIT_ZERO;
    bg_plane_pt_at(&projected, view_plane, fx, fy);
    VMOVE(dst, projected);
}

static SbColor
ged_obol_polygon_color_from_bu(const struct bu_color *color)
{
    fastf_t rgb[3] = {0.0, 0.0, 0.0};
    if (color)
	bu_color_to_rgb_floats(color, rgb);
    return SbColor(static_cast<float>(rgb[0]),
		   static_cast<float>(rgb[1]),
		   static_cast<float>(rgb[2]));
}

static void
ged_obol_polygon_color_to_bu(struct bu_color *dst, const SbColor &src)
{
    if (!dst)
	return;
    fastf_t rgb[3] = {
	static_cast<fastf_t>(src[0]),
	static_cast<fastf_t>(src[1]),
	static_cast<fastf_t>(src[2])
    };
    bu_color_from_rgb_floats(dst, rgb);
}

static int
ged_obol_polygon_record_to_ged(
    BObolViewController *controller,
    ged_view_polygon_ref ref,
    const BObolPolygonRecord &src,
    struct ged_view_polygon_record *dst)
{
    if (!controller || !dst)
	return 0;

    const char *name = controller->polygons().name(src.handle);
    if (!name)
	return 0;

    memset(dst, 0, sizeof(*dst));
    dst->ref = ref;
    dst->name = name;
    dst->type = static_cast<enum ged_view_polygon_type>(
	static_cast<int>(src.type));
    dst->selected = src.selected ? 1 : 0;
    dst->fill_flag = src.fill ? 1 : 0;
    V2SET(dst->fill_dir, src.fillSlope[0], src.fillSlope[1]);
    dst->fill_delta = src.fillSpacing;
    ged_obol_polygon_color_to_bu(&dst->fill_color, src.fillColor);
    ged_obol_rgb_from_color(src.edgeColor, dst->edge_color);
    dst->curr_contour_i = src.currentContour;
    dst->curr_point_i = src.currentPoint;
    dst->first_contour_open = src.firstContourOpen ? 1 : 0;
    dst->contour_count = src.contourCount;
    dst->point_count = src.pointCount;
    VSET(dst->origin_point, src.originPoint[0], src.originPoint[1],
	 src.originPoint[2]);
    HMOVE(dst->vp, src.viewPlane);
    dst->vZ = src.viewZ;
    dst->sketch_name = controller->polygons().sketchName(src.handle);
    dst->user_data = src.userData;
    return 1;
}

static BObolPolygonHandle
ged_obol_polygon_create_named(
    BObolViewController *controller,
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    int type,
    const point_t input_point)
{
    if (!controller || !name || !name[0] || !input_point)
	return BObolPolygonHandle();

    point_t origin = VINIT_ZERO;
    plane_t view_plane = HINIT_ZERO;
    ged_obol_polygon_project_point(origin, view_ctx, input_point,
				   &view_plane);

    SbVec3f obol_origin;
    ged_obol_polygon_point_from_ged(obol_origin, origin);
    return controller->polygons().create(name,
					 ged_obol_polygon_scope(local),
					 ged_obol_polygon_type_from_ged(type),
					 obol_origin,
					 view_plane,
					 0.0f);
}

extern "C" ged_view_polygon_ref
ged_draw_obol_view_context_polygon_find(struct ged_view_context *view_ctx, const char *name)
{
    if (!name)
	return GED_VIEW_POLYGON_REF_NULL;

    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BObolPolygonHandle local =
	    local_controller->polygons().find(name,
					      BOBOL_FEATURE_SCOPE_LOCAL);
	if (local.isValid())
	    return ged_obol_ged_view_polygon_ref(view_ctx, local_controller, local,
		    1);
    }

    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (!shared_controller)
	return GED_VIEW_POLYGON_REF_NULL;
    return ged_obol_ged_view_polygon_ref(view_ctx, shared_controller,
				    shared_controller->polygons().find(name,
					    BOBOL_FEATURE_SCOPE_SHARED), 0);
}

extern "C" ged_view_polygon_ref
ged_draw_obol_view_context_polygon_find_scoped(
    struct ged_view_context *view_ctx,
    const char *name,
    int local_only)
{
    if (!name)
	return GED_VIEW_POLYGON_REF_NULL;

    if (local_only) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	if (!controller)
	    return GED_VIEW_POLYGON_REF_NULL;
	return ged_obol_ged_view_polygon_ref(view_ctx, controller,
					controller->polygons().find(name, BOBOL_FEATURE_SCOPE_LOCAL), 1);
    }

    return ged_draw_obol_view_context_polygon_find(view_ctx, name);
}

extern "C" ged_view_polygon_ref
ged_draw_obol_view_context_polygon_create(
    struct ged_view_context *view_ctx,
    const char *name,
    int local,
    enum ged_view_polygon_type type,
    const point_t screen_point)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    return ged_obol_ged_view_polygon_ref(view_ctx, controller,
				    ged_obol_polygon_create_named(controller, view_ctx, name, local,
					    type, screen_point), local);
}

extern "C" ged_view_polygon_ref
ged_draw_obol_view_context_polygon_import_sketch(
    const char *name,
    struct db_i *dbip,
    struct directory *dp,
    struct ged_view_context *view_ctx,
    int local)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_scope(view_ctx, local, 1);
    if (!controller || !name)
	return GED_VIEW_POLYGON_REF_NULL;
    return ged_obol_ged_view_polygon_ref(view_ctx, controller,
				    controller->polygons().importSketch(name,
					    ged_obol_polygon_scope(local), dbip, dp), local);
}

extern "C" int
ged_draw_obol_view_context_polygon_set_current(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    long contour_i,
    long point_i)
{
    ged_view_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ged_ref);
    if (!controller)
	return 0;
    return controller->polygons().setCurrent(
	       ged_obol_polygon_handle_from_ged_ref(ged_ref),
	       contour_i, point_i) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_polygon_set_contour_open(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    long contour_i,
    int open)
{
    ged_view_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ged_ref);
    if (!controller)
	return 0;
    return controller->polygons().setContourOpen(
	       ged_obol_polygon_handle_from_ged_ref(ged_ref),
	       contour_i, open ? TRUE : FALSE) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_context_polygon_area(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    fastf_t *area)
{
    if (area)
	*area = 0.0;
    if (!area || !view_ctx)
	return 0;

    ged_view_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ged_ref);
    if (!controller)
	return 0;

    *area = controller->polygons().area(
	ged_obol_polygon_handle_from_ged_ref(ged_ref),
		bv_scale_get(ged_obol_bv_const(view_ctx)));
    return 1;
}

extern "C" int
ged_draw_obol_view_context_polygon_overlap(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    const char *other_name,
    const struct bn_tol *tol,
    int *overlap)
{
    if (overlap)
	*overlap = 0;
    if (!view_ctx || !other_name || !tol || !overlap)
	return 0;

    ged_view_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ged_ref);
    if (!controller)
	return 0;

    BObolPolygonHandle other =
	controller->polygons().find(other_name);
    if (!other.isValid())
	return 0;

    *overlap = controller->polygons().overlaps(
		   ged_obol_polygon_handle_from_ged_ref(ged_ref),
		   other,
		   *tol,
		   bv_scale_get(ged_obol_bv_const(view_ctx))) ? 1 : 0;
    return 1;
}

static size_t
ged_draw_obol_polygon_snap_count(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref exclude,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    BObolPolygonHandle exclude_handle;
    if (ged_obol_polygon_controller_from_ged_ref(view_ctx, exclude) == controller)
	exclude_handle = ged_obol_polygon_handle_from_ged_ref(exclude);
    return controller->polygons().snapCount(exclude_handle);
}

static int
ged_draw_obol_polygon_clear_point_selection(
    struct ged_view_context *view_ctx,
    void *UNUSED(data))
{
    BObolViewController *local = ged_obol_view_controller_for_context(view_ctx);
    BObolViewController *shared =
	ged_obol_shared_view_controller_for_context(view_ctx);
    int changed = 0;
    if (local)
	changed |= local->polygons().clearAllPointSelections() ? 1 : 0;
    if (shared && shared != local)
	changed |= shared->polygons().clearAllPointSelections() ? 1 : 0;
    return changed;
}

static int
ged_draw_obol_polygon_clear_selection(
    struct ged_view_context *view_ctx,
    void *UNUSED(data))
{
    BObolViewController *local = ged_obol_view_controller_for_context(view_ctx);
    BObolViewController *shared =
	ged_obol_shared_view_controller_for_context(view_ctx);
    int changed = 0;
    if (local)
	changed |= local->polygons().clearSelection() ? 1 : 0;
    if (shared && shared != local)
	changed |= shared->polygons().clearSelection() ? 1 : 0;
    return changed;
}

static int
ged_draw_obol_polygon_update(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int utype,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;
    return controller->polygons().update(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_update_screen_pt(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int x,
    int y,
    int utype,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller || !view_ctx)
	return 0;

    point_t model_point = VINIT_ZERO;
    if (!bv_screen_to_model(model_point, ged_obol_bv_const(view_ctx),
				       (fastf_t)x, (fastf_t)y))
	return 0;

    return controller->polygons().updateModelPoint(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       SbVec3f(static_cast<float>(model_point[X]),
		       static_cast<float>(model_point[Y]),
		       static_cast<float>(model_point[Z])),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_update_model_pt(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    const point_t model_point,
    int utype,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller || !model_point)
	return 0;

    return controller->polygons().updateModelPoint(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       SbVec3f(static_cast<float>(model_point[X]),
		       static_cast<float>(model_point[Y]),
		       static_cast<float>(model_point[Z])),
	       ged_obol_polygon_update_from_ged(utype)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_move(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const point_t current_point,
    const point_t previous_point,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller || !current_point || !previous_point)
	return 0;
    return controller->polygons().move(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       SbVec3f(static_cast<float>(current_point[X]),
		       static_cast<float>(current_point[Y]),
		       static_cast<float>(current_point[Z])),
	       SbVec3f(static_cast<float>(previous_point[X]),
		       static_cast<float>(previous_point[Y]),
		       static_cast<float>(previous_point[Z]))) ? 1 : 0;
}

static int
ged_draw_obol_polygon_set_name(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const char *name,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller || !name)
	return 0;
    return controller->polygons().rename(
	       ged_obol_polygon_handle_from_ged_ref(ref), name) ? 1 : 0;
}

static int
ged_draw_obol_polygon_set_selected(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    int selected,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    return controller ? controller->polygons().setSelected(
	ged_obol_polygon_handle_from_ged_ref(ref), selected ? TRUE : FALSE) : 0;
}

static int
ged_draw_obol_polygon_set_visual(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const struct bu_color *edge_color,
    const struct bu_color *fill_color,
    fastf_t fill_slope_x,
    fastf_t fill_slope_y,
    fastf_t fill_density,
    fastf_t vZ,
    int fill_flag,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;

    BObolPolygonVisual visual;
    (void)controller->polygons().visual(
	ged_obol_polygon_handle_from_ged_ref(ref), visual);
    if (edge_color)
	visual.edgeColor = ged_obol_polygon_color_from_bu(edge_color);
    if (fill_color)
	visual.fillColor = ged_obol_polygon_color_from_bu(fill_color);
    visual.fillSlope = SbVec2f(static_cast<float>(fill_slope_x),
			       static_cast<float>(fill_slope_y));
    visual.fillSpacing = static_cast<float>(fill_density);
    visual.viewZ = static_cast<float>(vZ);
    if (fill_flag)
	visual.fillFlags |= BOBOL_POLYGON_FILL_HATCH;
    else
	visual.fillFlags &= ~BOBOL_POLYGON_FILL_HATCH;
    visual.fill = (visual.fillFlags & BOBOL_POLYGON_FILL_HATCH) ?
		  TRUE : FALSE;
    return controller->polygons().setVisual(
	       ged_obol_polygon_handle_from_ged_ref(ref), visual) ? 1 : 0;
}

static int
ged_draw_obol_polygon_set_open(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    int open,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;
    return controller->polygons().setAllContoursOpen(
	       ged_obol_polygon_handle_from_ged_ref(ref),
	       open ? TRUE : FALSE) ? 1 : 0;
}

static int
ged_draw_obol_polygon_close(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *data)
{
    return ged_draw_obol_polygon_set_open(view_ctx, ref, 0, data);
}

static int
ged_draw_obol_polygon_clear_selected_point(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;
    return controller->polygons().clearSelectedPoint(
	       ged_obol_polygon_handle_from_ged_ref(ref)) ? 1 : 0;
}

static int
ged_draw_obol_polygon_remove(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;
    return controller->polygons().remove(
	       ged_obol_polygon_handle_from_ged_ref(ref)) ? 1 : 0;
}

static void *
ged_draw_obol_polygon_user_data(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    return controller ?
	   controller->polygons().userData(
	       ged_obol_polygon_handle_from_ged_ref(ref)) : NULL;
}

static int
ged_draw_obol_polygon_user_data_set(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *user_data,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;
    return controller->polygons().setUserData(
	       ged_obol_polygon_handle_from_ged_ref(ref), user_data) ? 1 : 0;
}

static int
ged_draw_obol_polygon_csg(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref target,
    ged_view_polygon_ref stencil,
    enum bg_polygon_boolean_op op,
    void *UNUSED(data))
{
    BObolViewController *target_controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, target);
    BObolViewController *stencil_controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, stencil);
    if (!target_controller || target_controller != stencil_controller)
	return 0;
    return target_controller->polygons().csg(
	       ged_obol_polygon_handle_from_ged_ref(target),
	       ged_obol_polygon_handle_from_ged_ref(stencil),
	       op) ? 1 : 0;
}

static struct directory *
ged_draw_obol_polygon_export_sketch(
    struct ged_view_context *view_ctx,
    struct db_i *dbip,
    const char *name,
    ged_view_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller || !dbip || !name)
	return NULL;
    return controller->polygons().exportSketch(
	       ged_obol_polygon_handle_from_ged_ref(ref), dbip, name) ?
	   db_lookup(dbip, name, LOOKUP_QUIET) : NULL;
}

static struct directory *
ged_draw_obol_polygon_update_sketch(
    struct ged_view_context *view_ctx,
    struct db_i *dbip,
    ged_view_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller || !dbip)
	return NULL;
    BObolPolygonHandle handle = ged_obol_polygon_handle_from_ged_ref(ref);
    const char *name = controller->polygons().sketchName(handle);
    if (!name || !name[0])
	return NULL;
    return controller->polygons().updateSketch(handle, dbip, name) ?
	db_lookup(dbip, name, LOOKUP_QUIET) : NULL;
}

static int
ged_draw_obol_polygon_sketch_name_set(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const char *name,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;
    return controller->polygons().setSketchName(
	ged_obol_polygon_handle_from_ged_ref(ref), name ? name : "") ? 1 : 0;
}

static int
ged_draw_obol_polygon_snap_exclude_set(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *UNUSED(data))
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return 0;
    if (ged_view_polygon_ref_is_null(ref))
	return controller->polygons().setSnapExclude(
	    BObolPolygonHandle()) ? 1 : 0;
    if (ged_obol_polygon_controller_from_ged_ref(view_ctx, ref) != controller)
	return 0;
    return controller->polygons().setSnapExclude(
	       ged_obol_polygon_handle_from_ged_ref(ref)) ? 1 : 0;
}

extern "C" ged_view_polygon_ref
ged_draw_obol_view_context_polygon_select(
    struct ged_view_context *view_ctx,
    const point_t model_point)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !model_point)
	return GED_VIEW_POLYGON_REF_NULL;
    return ged_obol_ged_view_polygon_ref_for_handle(view_ctx, controller,
	    controller->polygons().selectAtModelPoint(SbVec3f(
		    static_cast<float>(model_point[X]),
		    static_cast<float>(model_point[Y]),
		    static_cast<float>(model_point[Z]))));
}

extern "C" ged_view_polygon_ref
ged_draw_obol_view_context_polygon_dup(
    struct ged_view_context *view_ctx,
    const char *name,
    const char *new_name)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !name || !new_name)
	return GED_VIEW_POLYGON_REF_NULL;
    BObolPolygonHandle src = controller->polygons().find(name);
    return ged_obol_ged_view_polygon_ref_for_handle(view_ctx, controller,
	    controller->polygons().duplicate(src, new_name));
}

extern "C" int
ged_draw_obol_view_polygon_record_get(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    struct ged_view_polygon_record *record)
{
    ged_view_polygon_ref ged_ref = ref;
    BObolViewController *controller =
	ged_obol_polygon_controller_from_ged_ref(view_ctx, ged_ref);
    if (!controller || !record)
	return 0;

    BObolPolygonRecord obol_record;
    if (!controller->polygons().record(
	    ged_obol_polygon_handle_from_ged_ref(ged_ref), obol_record))
	return 0;
    return ged_obol_polygon_record_to_ged(controller, ref,
	    obol_record, record);
}

struct ged_obol_polygon_ged_visit_context {
    struct ged_view_context *view_ctx;
    BObolViewController *controller;
    ged_view_polygon_record_cb callback;
    void *data;
    BObolFeatureScope scope;
};

static int
ged_draw_obol_polygon_visit_ged_cb(const BObolPolygonRecord &obol_record,
				   void *data)
{
    ged_obol_polygon_ged_visit_context *ctx =
	static_cast<ged_obol_polygon_ged_visit_context *>(data);
    if (!ctx || !ctx->controller || !ctx->callback)
	return 0;
    if (obol_record.scope != ctx->scope)
	return 1;

    ged_view_polygon_ref ref =
	ged_obol_ged_view_polygon_ref(ctx->view_ctx, ctx->controller,
		obol_record.handle,
		obol_record.scope == BObolFeatureScope::Local ? 1 : 0);
    struct ged_view_polygon_record record;
    if (!ged_obol_polygon_record_to_ged(ctx->controller, ref, obol_record,
	    &record))
	return 1;
    return ctx->callback(ref, &record, ctx->data);
}

extern "C" void
ged_draw_obol_view_context_polygon_visit_records(
    struct ged_view_context *view_ctx,
    ged_view_polygon_record_cb callback,
    void *data)
{
    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if ((!local_controller && !shared_controller) || !callback)
	return;

    ged_obol_polygon_ged_visit_context ctx;
    ctx.view_ctx = view_ctx;
    ctx.callback = callback;
    ctx.data = data;
    if (local_controller) {
	ctx.controller = local_controller;
	ctx.scope = BObolFeatureScope::Local;
	local_controller->polygons().visitRecords(
		ged_draw_obol_polygon_visit_ged_cb, &ctx);
    }
    if (shared_controller) {
	ctx.controller = shared_controller;
	ctx.scope = BObolFeatureScope::Shared;
	shared_controller->polygons().visitRecords(
		ged_draw_obol_polygon_visit_ged_cb, &ctx);
    }
}

extern "C" size_t
ged_draw_obol_view_context_polygon_snap_count(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref exclude)
{
    return ged_draw_obol_polygon_snap_count(view_ctx,
	    exclude, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_clear_point_selection(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_polygon_clear_point_selection(view_ctx, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_clear_selection(struct ged_view_context *view_ctx)
{
    return ged_draw_obol_polygon_clear_selection(view_ctx, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_snap_exclude_set(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_snap_exclude_set(view_ctx,
	    ref, NULL);
}

extern "C" struct directory *
ged_draw_obol_view_polygon_export_sketch(
    struct ged_view_context *view_ctx,
    struct db_i *dbip,
    const char *name,
    ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_export_sketch(view_ctx, dbip, name,
	    ref, NULL);
}

extern "C" struct directory *
ged_draw_obol_view_polygon_update_sketch(
    struct ged_view_context *view_ctx,
    struct db_i *dbip,
    ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_update_sketch(view_ctx, dbip, ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_sketch_name_set(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const char *name)
{
    return ged_draw_obol_polygon_sketch_name_set(view_ctx, ref, name, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_update(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int op)
{
    return ged_draw_obol_polygon_update(
	    ref, view_ctx, op, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_update_screen_pt(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    int x,
    int y,
    int op)
{
    return ged_draw_obol_polygon_update_screen_pt(
	    ref, view_ctx, x, y, op, NULL);
}

extern "C" int
ged_draw_obol_view_context_polygon_update_model_pt(
    ged_view_polygon_ref ref,
    struct ged_view_context *view_ctx,
    const point_t model_point,
    int op)
{
    return ged_draw_obol_polygon_update_model_pt(
	    ref, view_ctx, model_point, op, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_move(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const point_t current_point,
    const point_t previous_point)
{
    return ged_draw_obol_polygon_move(
	    view_ctx, ref, current_point,
	    previous_point, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_name(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const char *name)
{
    return ged_draw_obol_polygon_set_name(view_ctx, ref, name, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_selected(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    int selected)
{
    return ged_draw_obol_polygon_set_selected(view_ctx, ref, selected, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_visual(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    const struct bu_color *edge_color,
    const struct bu_color *fill_color,
    fastf_t fill_slope_x,
    fastf_t fill_slope_y,
    fastf_t fill_density,
    fastf_t vZ,
    int fill_flag)
{
    return ged_draw_obol_polygon_set_visual(view_ctx, ref, edge_color,
	    fill_color,
	    fill_slope_x, fill_slope_y, fill_density, vZ, fill_flag, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_set_all_contours_open(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    int open)
{
    return ged_draw_obol_polygon_set_open(view_ctx, ref, open, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_close(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_close(view_ctx, ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_clear_selected_point(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_clear_selected_point(view_ctx, ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_remove(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_remove(view_ctx, ref, NULL);
}

extern "C" void *
ged_draw_obol_view_polygon_user_data(struct ged_view_context *view_ctx,
	ged_view_polygon_ref ref)
{
    return ged_draw_obol_polygon_user_data(view_ctx, ref, NULL);
}

extern "C" int
ged_draw_obol_view_polygon_user_data_set(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref ref,
    void *user_data)
{
    return ged_draw_obol_polygon_user_data_set(view_ctx, ref, user_data,
	    NULL);
}

extern "C" int
ged_draw_obol_view_polygon_csg(
    struct ged_view_context *view_ctx,
    ged_view_polygon_ref target,
    ged_view_polygon_ref stencil,
    enum bg_polygon_boolean_op op)
{
    return ged_draw_obol_polygon_csg(view_ctx, target,
	    stencil, op, NULL);
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
