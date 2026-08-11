/*      G E D _ O B O L _ P R E S E N T A T I O N _ P R I V A T E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file ged_obol_presentation_private.cpp
 *
 * Obol-backed retained features, faceplate, lighting, and edit previews.
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


static struct bv *
ged_obol_bv(struct ged_view_context *view_ctx)
{
    return bv_context_view((struct bv_context *)view_ctx);
}


static const struct bv *
ged_obol_bv_const(const struct ged_view_context *view_ctx)
{
    return bv_context_view_const((const struct bv_context *)view_ctx);
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


static SbColor
ged_obol_color_from_rgb(const unsigned char rgb[3])
{
    return SbColor(static_cast<float>(rgb[0]) / 255.0f,
	static_cast<float>(rgb[1]) / 255.0f,
	static_cast<float>(rgb[2]) / 255.0f);
}


static unsigned char
ged_obol_rgb_byte_from_int(int value)
{
    if (value < 0)
	return 0;
    if (value > 255)
	return 255;
    return static_cast<unsigned char>(value);
}


static BObolViewController *
ged_obol_view_controller_for_context(struct ged_view_context *view_ctx)
{
    return ged_bobol_view_controller(view_ctx);
}


static BObolViewController *
ged_obol_view_controller_ensure_for_context(struct ged_view_context *view_ctx,
	int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
	static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) : NULL;
    if (!gedp)
	return NULL;
    BObolViewController *controller = ged_bobol_view_controller(view_ctx);
    if (controller)
	return controller;
    if (!ged_draw_obol_render_endpoint_ensure_for_view(gedp, view_ctx,
	    sync_current_scene))
	return NULL;
    return ged_bobol_view_controller(view_ctx);
}


static BObolViewController *
ged_obol_shared_view_controller_ensure_for_context(
    struct ged_view_context *view_ctx, int sync_current_scene)
{
    struct ged *gedp = view_ctx ?
	static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) : NULL;
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
    struct ged *gedp = view_ctx ?
	static_cast<struct ged *>(ged_view_context_user_data_get(view_ctx)) : NULL;
    return gedp ? ged_draw_obol_controller(gedp) : NULL;
}


BObolFeatureStyle
ged_obol_feature_style_from_ged(
    const struct ged_view_feature_style *style)
{
    BObolFeatureStyle out;
    if (!style)
	return out;

    if (style->visible != -1) {
	out.hasVisible = TRUE;
	out.visible = style->visible ? TRUE : FALSE;
    }
    if (style->selectable != -1) {
	out.hasSelectable = TRUE;
	out.selectable = style->selectable ? TRUE : FALSE;
    }
    if (style->color_valid) {
	out.hasColor = TRUE;
	out.color = ged_obol_color_from_rgb(style->color);
    }
    if (style->line_width >= 0) {
	out.hasLineWidth = TRUE;
	out.lineWidth = style->line_width;
    }
    if (style->line_style >= 0) {
	out.hasLineStyle = TRUE;
	out.lineStyle = style->line_style;
    }
    if (style->arrow != -1) {
	out.hasArrow = TRUE;
	out.arrow = style->arrow ? TRUE : FALSE;
    }
    if (style->arrow_tip_length >= 0.0 || style->arrow_tip_width >= 0.0) {
	out.hasArrowTip = TRUE;
	out.arrowTipLength = style->arrow_tip_length >= 0.0 ?
			     static_cast<float>(style->arrow_tip_length) : 0.0f;
	out.arrowTipWidth = style->arrow_tip_width >= 0.0 ?
			    static_cast<float>(style->arrow_tip_width) : 0.0f;
    }
    return out;
}

int32_t
ged_obol_line_command_from_ged(int command, size_t index)
{
    if (command == GED_DRAW_VIEW_LINE_DRAW ||
	command == BG_GEOMETRY_LINE_DRAW)
	return static_cast<int32_t>(BObolLineCommand::Draw);
    if (command == GED_DRAW_VIEW_LINE_POINT_DRAW ||
	command == BG_GEOMETRY_POINT_DRAW)
	return static_cast<int32_t>(BObolLineCommand::Point);
    if (command == GED_DRAW_VIEW_LINE_MOVE ||
	command == BG_GEOMETRY_LINE_MOVE)
	return static_cast<int32_t>(BObolLineCommand::Move);
    return static_cast<int32_t>(index ? BObolLineCommand::Draw :
				BObolLineCommand::Move);
}

static int
ged_obol_line_command_to_ged(int32_t command)
{
    if (command == static_cast<int32_t>(BObolLineCommand::Draw))
	return GED_DRAW_VIEW_LINE_DRAW;
    if (command == static_cast<int32_t>(BObolLineCommand::Point))
	return GED_DRAW_VIEW_LINE_POINT_DRAW;
    return GED_DRAW_VIEW_LINE_MOVE;
}

std::vector<SbVec3f>
ged_obol_points_from_ged(const point_t *points, size_t point_count)
{
    std::vector<SbVec3f> out;
    if (!points || !point_count)
	return out;

    out.reserve(point_count);
    for (size_t i = 0; i < point_count; i++)
	out.push_back(SbVec3f(
			  static_cast<float>(points[i][X]),
			  static_cast<float>(points[i][Y]),
			  static_cast<float>(points[i][Z])));
    return out;
}

std::vector<int32_t>
ged_obol_commands_from_ged(const int *cmds, size_t point_count)
{
    std::vector<int32_t> out;
    out.reserve(point_count);
    for (size_t i = 0; i < point_count; i++) {
	const int command = cmds ? cmds[i] : -1;
	out.push_back(ged_obol_line_command_from_ged(command, i));
    }
    return out;
}

std::vector<int32_t>
ged_obol_indices_from_ged(const int *indices, size_t index_count)
{
    std::vector<int32_t> out;
    if (!indices || !index_count)
	return out;

    out.reserve(index_count);
    for (size_t i = 0; i < index_count; i++)
	out.push_back(static_cast<int32_t>(indices[i]));
    return out;
}

std::vector<SbVec3f>
ged_obol_vectors_from_ged(const vect_t *vectors, size_t vector_count)
{
    std::vector<SbVec3f> out;
    if (!vectors || !vector_count)
	return out;

    out.reserve(vector_count);
    for (size_t i = 0; i < vector_count; i++)
	out.push_back(SbVec3f(
			  static_cast<float>(vectors[i][X]),
			  static_cast<float>(vectors[i][Y]),
			  static_cast<float>(vectors[i][Z])));
    return out;
}

static BObolFeatureHandle
ged_obol_feature_handle(BObolViewController *controller,
			struct ged_view_context *view_ctx,
			const char *name)
{
    if (!controller || !name)
	return BObolFeatureHandle();

    BObolFeatureOwner owner;
    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";

    BObolFeatureHandle local =
	controller->features().findOwned(name, BOBOL_FEATURE_SCOPE_LOCAL,
					 &owner);
    if (local.isValid())
	return local;

    BObolFeatureHandle shared =
	controller->features().find(name, BOBOL_FEATURE_SCOPE_SHARED);
    if (shared.isValid())
	return shared;

    return controller->features().find(name);
}

struct ged_obol_feature_lookup {
    BObolViewController *controller;
    BObolFeatureHandle handle;
};

static ged_obol_feature_lookup
ged_obol_feature_lookup_for_context(struct ged_view_context *view_ctx, const char *name)
{
    ged_obol_feature_lookup out;
    out.controller = NULL;
    out.handle = BObolFeatureHandle();
    if (!name)
	return out;

    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (local_controller) {
	BObolFeatureHandle handle =
	    ged_obol_feature_handle(local_controller, view_ctx, name);
	if (handle.isValid()) {
	    out.controller = local_controller;
	    out.handle = handle;
	    return out;
	}
    }

    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if (shared_controller && shared_controller != local_controller) {
	BObolFeatureHandle handle =
	    shared_controller->features().find(name,
					       BOBOL_FEATURE_SCOPE_SHARED);
	if (handle.isValid()) {
	    out.controller = shared_controller;
	    out.handle = handle;
	}
    }

    return out;
}

static BObolFeatureScope
ged_obol_feature_scope(int local)
{
    return local ? BObolFeatureScope::Local : BObolFeatureScope::Shared;
}

static BObolFeatureOwner
ged_obol_feature_owner(struct ged_view_context *view_ctx, int local)
{
    BObolFeatureOwner owner;
    if (!local)
	return owner;

    owner.ownerToken = view_ctx;
    owner.ownerId = ged_obol_view_scope_name(view_ctx).c_str();
    owner.ownerRole = "view";
    return owner;
}

BObolOverlayInfo
ged_obol_model_overlay_info(struct ged_view_context *view_ctx,
			    BObolOverlayClass overlay_class,
			    BObolOverlayLifecycle lifecycle,
			    BObolOverlayOrder order,
			    const char *source_path)
{
    BObolOverlayInfo info;
    info.isOverlay = TRUE;
    info.ownerToken = view_ctx;
    info.role = BObolOverlayRole::Model;
    info.overlayClass = overlay_class;
    info.lifecycle = lifecycle;
    info.order = order;
    info.sortOrder = 0;
    info.sourcePath = source_path ? source_path : "";
    return info;
}

int
ged_obol_feature_mark_overlay(BObolViewController *controller,
			      BObolFeatureHandle handle,
			      const BObolOverlayInfo &overlay)
{
    if (!controller || !handle.isValid())
	return 0;
    return controller->features().setOverlayInfo(handle, overlay) ? 1 : 0;
}


static BObolFeatureHandle
ged_obol_publish_line_set(BObolViewController *controller,
			  struct ged_view_context *view_ctx,
			  const char *name,
			  int local,
			  const std::vector<SbVec3f> &points,
			  const std::vector<int32_t> &commands,
			  const BObolFeatureStyle *style)
{
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, local);
    return controller->features().publishLineSet(name,
	    ged_obol_feature_scope(local), points, commands, style,
	    local ? &owner : NULL);
}

static int
ged_obol_remove_feature(BObolViewController *controller,
			struct ged_view_context *view_ctx,
			const char *name,
			int local_mode)
{
    if (!controller || !name)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle;
    if (local_mode > 0) {
	handle = controller->features().findOwned(name,
		 BOBOL_FEATURE_SCOPE_LOCAL, &owner);
    } else if (local_mode == 0) {
	handle = controller->features().find(name,
					     BOBOL_FEATURE_SCOPE_SHARED);
    } else {
	handle = ged_obol_feature_handle(controller, view_ctx, name);
    }

    return handle.isValid() ? (controller->features().remove(handle) ? 1 : 0) : 0;
}

static void
ged_obol_point_from_sb(point_t dst, const SbVec3f &src)
{
    VSET(dst, src[0], src[1], src[2]);
}

BObolLabel
ged_obol_label_from_hud(const struct ged_diagnostic_hud_label &label)
{
    BObolLabel out;
    out.text = label.text ? label.text : "";
    out.point = SbVec3f(
		    static_cast<float>(label.position[0]),
		    static_cast<float>(label.position[1]),
		    0.0f);
    unsigned char rgb[3] = {
	static_cast<unsigned char>(label.color[0] < 0 ? 0 :
	(label.color[0] > 255 ? 255 : label.color[0])),
	static_cast<unsigned char>(label.color[1] < 0 ? 0 :
	(label.color[1] > 255 ? 255 : label.color[1])),
	static_cast<unsigned char>(label.color[2] < 0 ? 0 :
	(label.color[2] > 255 ? 255 : label.color[2]))
    };
    out.hasColor = TRUE;
    out.color = ged_obol_color_from_rgb(rgb);
    out.fontSize = label.font_size > 0.0 ?
		   static_cast<float>(label.font_size) : 12.0f;
    out.sourceId = label.source_id;
    return out;
}

static int
ged_obol_rgb_is_zero(const int rgb[3])
{
    return !rgb || (rgb[0] == 0 && rgb[1] == 0 && rgb[2] == 0);
}

static unsigned char
ged_obol_clamp_color_int(int v)
{
    return static_cast<unsigned char>(v < 0 ? 0 : (v > 255 ? 255 : v));
}

static SbColor
ged_obol_color_from_int_rgb(const int rgb[3],
			    int fallback_r,
			    int fallback_g,
			    int fallback_b)
{
    int r = fallback_r;
    int g = fallback_g;
    int b = fallback_b;
    if (!ged_obol_rgb_is_zero(rgb)) {
	r = rgb[0];
	g = rgb[1];
	b = rgb[2];
    }
    const unsigned char crgb[3] = {
	ged_obol_clamp_color_int(r),
	ged_obol_clamp_color_int(g),
	ged_obol_clamp_color_int(b)
    };
    return ged_obol_color_from_rgb(crgb);
}

static BObolFeatureStyle
ged_obol_faceplate_style(const int rgb[3],
			 int fallback_r,
			 int fallback_g,
			 int fallback_b,
			 int line_width)
{
    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    style.hasColor = TRUE;
    style.color = ged_obol_color_from_int_rgb(rgb, fallback_r, fallback_g,
		  fallback_b);
    style.hasLineWidth = TRUE;
    style.lineWidth = line_width > 0 ? line_width : 1;
    /* Faceplate geometry is specified in screen pixels and must neither
     * participate in nor write the model depth buffer.  Apart from making
     * the overlays camera-independent, the HUD wrapper is important for
     * composite controls such as the LoD progress meter: its narrower fill
     * line occupies the same screen location as the track and would
     * otherwise fail the depth test after the track was drawn. */
    style.hud = TRUE;
    return style;
}

static BObolOverlayInfo
ged_obol_faceplate_overlay_info(struct ged_view_context *view_ctx,
				BObolOverlayOrder order =
				    BObolOverlayOrder::Screen)
{
    BObolOverlayInfo info;
    info.isOverlay = TRUE;
    info.ownerToken = view_ctx;
    info.role = BObolOverlayRole::Screen;
    info.overlayClass = BObolOverlayClass::Faceplate;
    info.lifecycle = BObolOverlayLifecycle::PerView;
    info.order = order;
    info.sortOrder = 0;
    info.sourcePath = "_faceplate";
    return info;
}

/* Map a faceplate overlay coordinate (GED2PM1 space, -1..1 in both axes) to the
 * pixel coordinates an SoHUDKit expects (origin bottom-left, 0..width/height),
 * matching how HUD text labels are placed.  This keeps the yellow faceplate
 * lines screen-locked -- like legacy main's fixed glOrtho(-1,1,-1,1) overlay --
 * instead of routing them through the view->model transform (which made them
 * move and skew with the camera/aspect). */
static void
ged_obol_faceplate_to_pixel(SbVec3f &out,
			    struct ged_view_context *view_ctx,
			    fastf_t x,
			    fastf_t y)
{
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    out = SbVec3f(static_cast<float>((x + 1.0) * 0.5 * width),
		  static_cast<float>((y + 1.0) * 0.5 * height),
		  0.0f);
}

static void
ged_obol_faceplate_append_line(std::vector<SbVec3f> &points,
			       std::vector<int32_t> &commands,
			       struct ged_view_context *view_ctx,
			       fastf_t x1,
			       fastf_t y1,
			       fastf_t x2,
			       fastf_t y2)
{
    SbVec3f a;
    SbVec3f b;
    ged_obol_faceplate_to_pixel(a, view_ctx, x1, y1);
    ged_obol_faceplate_to_pixel(b, view_ctx, x2, y2);

    points.push_back(a);
    commands.push_back(static_cast<int32_t>(BObolLineCommand::Move));
    points.push_back(b);
    commands.push_back(static_cast<int32_t>(BObolLineCommand::Draw));
}

static BObolFeatureHandle
ged_obol_faceplate_publish_lines(BObolViewController *controller,
				 struct ged_view_context *view_ctx,
				 const char *name,
				 const std::vector<SbVec3f> &points,
				 const std::vector<int32_t> &commands,
				 const BObolFeatureStyle &style)
{
    if (!controller || !name || points.empty() || commands.empty())
	return BObolFeatureHandle();

    BObolFeatureHandle handle = ged_obol_publish_line_set(controller,
				  view_ctx, name, 1, points, commands, &style);
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
    return handle;
}

static BObolFeatureHandle
ged_obol_faceplate_publish_hud_labels(BObolViewController *controller,
				      struct ged_view_context *view_ctx,
				      const char *name,
				      const std::vector<BObolLabel> &labels,
				      const BObolFeatureStyle &style)
{
    if (!controller || !name || labels.empty())
	return BObolFeatureHandle();

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle =
	controller->features().publishHudLabels(name,
	    BObolFeatureScope::Local, labels, &style, &owner);
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
    return handle;
}

static void
ged_obol_faceplate_remove(BObolViewController *controller,
			  struct ged_view_context *view_ctx,
			  const char *name)
{
    (void)ged_obol_remove_feature(controller, view_ctx, name, 1);
}

static BObolLabel
ged_obol_faceplate_label(struct ged_view_context *view_ctx,
			 const char *text,
			 fastf_t x,
			 fastf_t y,
			 const int rgb[3],
			 int fallback_r,
			 int fallback_g,
			 int fallback_b,
			 int font_size,
			 int anchor = 0)
{
    BObolLabel label;
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    fastf_t px = x;
    fastf_t py = y;
    if (width > 0 && height > 0) {
	px = (x + 1.0) * 0.5 * (fastf_t)width;
	py = (y + 1.0) * 0.5 * (fastf_t)height;
    }
    label.text = text ? text : "";
    label.point = SbVec3f(static_cast<float>(px),
			  static_cast<float>(py),
			  0.0f);
    label.hasColor = TRUE;
    label.color = ged_obol_color_from_int_rgb(rgb, fallback_r, fallback_g,
		  fallback_b);
    label.fontSize = static_cast<float>(font_size > 0 ? font_size : 20);
    label.anchor = anchor;
    return label;
}

static void
ged_obol_faceplate_sync_center_dot(BObolViewController *controller,
				   struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/center_dot";
    struct bv_other_state state = BV_OTHER_STATE_INIT;
    if (!bv_center_dot_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.gos_draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(4);
    commands.reserve(4);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   -0.01, 0.0, 0.01, 0.0);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   0.0, -0.01, 0.0, 0.01);
    BObolFeatureStyle style = ged_obol_faceplate_style(
				    state.gos_line_color, 255, 255, 0, 1);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx, name,
					   points, commands, style);
}

static void
ged_obol_faceplate_sync_interactive_rect(BObolViewController *controller,
					 struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/interactive_rect";
    struct bv_interactive_rect_state state = BV_INTERACTIVE_RECT_STATE_INIT;
    if (!bv_interactive_rect_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.draw || (ZERO(state.width) && ZERO(state.height))) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    const fastf_t aspect = width > 0 && height > 0 ?
	(fastf_t)width / (fastf_t)height : 1.0;
    const fastf_t x0 = state.x;
    const fastf_t x1 = state.x + state.width;
    const fastf_t y0 = state.y * aspect;
    const fastf_t y1 = (state.y + state.height) * aspect;
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(8);
    commands.reserve(8);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x0, y0, x0, y1);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x0, y1, x1, y1);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x1, y1, x1, y0);
    ged_obol_faceplate_append_line(points, commands, view_ctx, x1, y0, x0, y0);

    BObolFeatureStyle style = ged_obol_faceplate_style(state.color,
	255, 255, 255, state.line_width > 0 ? state.line_width : 1);
    style.lineStyle = state.line_style;
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx, name,
					   points, commands, style);
}

static void
ged_obol_faceplate_sync_grid(BObolViewController *controller,
			     struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/grid";
    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    if (!bv_grid_state_get(&grid, ged_obol_bv_const(view_ctx)) || !grid.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLGrid *node = new SoBRLGrid;
    node->ref();
    node->overlayId = name;
    if (!bobol_grid_configure_from_view_context(node, &grid, view_ctx) ||
	node->getTotalSegmentCount() <= 0) {
	node->unref();
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    BObolFeatureStyle style;
    style.hasVisible = TRUE;
    style.visible = TRUE;
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle =
	controller->features().publishCustomNode(name,
	    BObolFeatureScope::Local, node, &style, &owner);
    node->unref();
    (void)ged_obol_feature_mark_overlay(controller, handle,
					ged_obol_faceplate_overlay_info(view_ctx));
}

static void
ged_obol_faceplate_sync_adc(BObolViewController *controller,
			    struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/adc";
    struct bv_adc_state state = BV_ADC_STATE_INIT;
    if (!bv_adc_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    SoBRLADC *node = new SoBRLADC;
    node->ref();
    node->overlayId = name;
    node->center = SbVec3f(static_cast<float>(state.pos_model[X]),
			  static_cast<float>(state.pos_model[Y]),
			  static_cast<float>(state.pos_model[Z]));
    node->angleDegrees = static_cast<float>(state.a1);
    node->distance = static_cast<float>(state.dst > SMALL_FASTF ?
					state.dst : 1.0);
    node->lineColor = SbColor(state.line_color[0] / 255.0f,
	state.line_color[1] / 255.0f, state.line_color[2] / 255.0f);
    node->tickColor = SbColor(state.tick_color[0] / 255.0f,
	state.tick_color[1] / 255.0f, state.tick_color[2] / 255.0f);
    node->lineWidth = state.line_width > 0 ? state.line_width : 1;
    node->visible = TRUE;
    node->rebuildGeometry();

    BObolFeatureStyle style = ged_obol_faceplate_style(
	state.line_color, 255, 255, 255, state.line_width);
    style.hasSelectable = TRUE;
    style.selectable = FALSE;
    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle =
	controller->features().publishCustomNode(name,
	    BObolFeatureScope::Local, node, &style, &owner);
    node->unref();
    (void)ged_obol_feature_mark_overlay(controller, handle,
	ged_obol_faceplate_overlay_info(view_ctx));
}

static std::vector<std::string>
ged_obol_faceplate_params_parts(struct ged_view_context *view_ctx,
				const struct bv_params_state *params)
{
    std::vector<std::string> parts;
    if (!view_ctx || !params)
	return parts;

    point_t center = VINIT_ZERO;
    mat_t center_mat;
    if (bv_center_mat_get(center_mat, ged_obol_bv_const(view_ctx)))
	MAT_DELTAS_GET_NEG(center, center_mat);
    VSCALE(center, center, bv_base2local_get(ged_obol_bv_const(view_ctx)));

    const char *ustr = bu_units_string(bv_local2base_get(ged_obol_bv_const(view_ctx)));
    if (!ustr)
	ustr = "";

    vect_t aet = VINIT_ZERO;
    (void)bv_aet_get(aet, ged_obol_bv_const(view_ctx));

    struct bu_vls text = BU_VLS_INIT_ZERO;
    if (params->draw_size) {
	bu_vls_sprintf(&text, "size[%s]: %.2f", ustr,
		      bv_size_get(ged_obol_bv_const(view_ctx)) *
		      bv_base2local_get(ged_obol_bv_const(view_ctx)));
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_center) {
	bu_vls_sprintf(&text, "center[%s]: (%.2f, %.2f, %.2f)",
		      ustr, V3ARGS(center));
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_az) {
	bu_vls_sprintf(&text, "az:%.2f", aet[0]);
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_el) {
	bu_vls_sprintf(&text, "el:%.2f", aet[1]);
	parts.push_back(bu_vls_cstr(&text));
    }
    if (params->draw_tw) {
	bu_vls_sprintf(&text, "tw:%.2f", aet[2]);
	parts.push_back(bu_vls_cstr(&text));
    }

    const uint64_t frametime = bv_frametime_get(ged_obol_bv_const(view_ctx));
    if (params->draw_fps && frametime > 0) {
	bu_vls_sprintf(&text, "FPS:%.1f",
		1000000000.0 / (fastf_t)frametime);
	parts.push_back(bu_vls_cstr(&text));
    }
    bu_vls_free(&text);
    return parts;
}

static void
ged_obol_faceplate_sync_params(BObolViewController *controller,
			       struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/params";
    struct bv_params_state params = BV_PARAMS_STATE_INIT;
    if (!bv_params_state_get(&params, ged_obol_bv_const(view_ctx)) ||
	!params.draw) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    if (!params.draw_size && !params.draw_center && !params.draw_az &&
	!params.draw_el && !params.draw_tw && !params.draw_fps) {
	params.draw_size = 1;
	params.draw_center = 1;
	params.draw_az = 1;
	params.draw_el = 1;
	params.draw_tw = 1;
    }

    const int font_size = params.font_size > 0 ? params.font_size : 20;
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    const int available_width = width > 16 ? width - 16 : width;
    const size_t max_chars = available_width > 0 ?
	std::max((size_t)1, (size_t)((fastf_t)available_width /
		(0.58 * (fastf_t)font_size))) : (size_t)80;

    const std::vector<std::string> parts =
	ged_obol_faceplate_params_parts(view_ctx, &params);
    std::vector<std::string> lines;
    std::string line;
    for (size_t i = 0; i < parts.size(); i++) {
	const size_t joined_len = line.size() + (line.empty() ? 0 : 1) +
	    parts[i].size();
	if (!line.empty() && joined_len > max_chars) {
	    lines.push_back(line);
	    line.clear();
	}
	if (!line.empty())
	    line.append(" ");
	line.append(parts[i]);
    }
    if (!line.empty())
	lines.push_back(line);

    std::vector<BObolLabel> labels;
    const fastf_t line_step = height > 0 ?
	(2.0 * (fastf_t)(font_size + 4) / (fastf_t)height) : 0.10;
    for (size_t i = 0; i < lines.size(); i++) {
	int line_font_size = font_size;
	if (available_width > 0 && !lines[i].empty()) {
	    const int fit_size = (int)((fastf_t)available_width /
		(0.58 * (fastf_t)lines[i].size()));
	    line_font_size = std::max(6, std::min(font_size, fit_size));
	}
	labels.push_back(ged_obol_faceplate_label(view_ctx, lines[i].c_str(),
		-0.98, -0.90 + (fastf_t)i * line_step,
		params.color, 255, 255, 0, line_font_size));
    }

    if (labels.empty()) {
	ged_obol_faceplate_remove(controller, view_ctx, name);
	return;
    }

    BObolFeatureStyle style = ged_obol_faceplate_style(
				    params.color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx, name,
	    labels, style);
}

static std::string
ged_obol_lod_compact_count(size_t value)
{
    struct bu_vls text = BU_VLS_INIT_ZERO;
    if (value >= 1000000)
	bu_vls_sprintf(&text, "%.1fM",
	    static_cast<double>(value) / 1000000.0);
    else if (value >= 1000)
	bu_vls_sprintf(&text, "%.1fk",
	    static_cast<double>(value) / 1000.0);
    else
	bu_vls_sprintf(&text, "%zu", value);
    const std::string result = bu_vls_cstr(&text);
    bu_vls_free(&text);
    return result;
}

static void
ged_obol_faceplate_sync_lod_progress(BObolViewController *controller,
				     struct ged_view_context *view_ctx)
{
    static const char track_name[] = "_faceplate/lod_progress_track";
    static const char fill_name[] = "_faceplate/lod_progress_fill";
    static const char label_name[] = "_faceplate/lod_progress_label";
    BObolLodConvergenceStatus status;
    controller->getLodConvergenceStatus(status);

    const bool show =
	status.phase != BOBOL_LOD_CONVERGENCE_IDLE ||
	status.backgroundPending || status.performanceLimited ||
	status.failedSourceCount > 0;
    if (!show) {
	ged_obol_faceplate_remove(controller, view_ctx, track_name);
	ged_obol_faceplate_remove(controller, view_ctx, fill_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name);
	return;
    }

    int color[3] = {96, 220, 255};
    struct bu_vls text = BU_VLS_INIT_ZERO;
    const int percent = static_cast<int>(
	std::floor(std::max(0.0f, std::min(1.0f, status.fraction)) *
	    100.0f + 0.5f));
    const size_t pending = status.pendingTasks > SIZE_MAX - status.inFlight ?
	SIZE_MAX : status.pendingTasks + status.inFlight;
	switch (status.phase) {
	case BOBOL_LOD_CONVERGENCE_DISCOVERING:
	    if (status.expectedLeafCount > status.availableLeafCount) {
		bu_vls_sprintf(&text, "Discovering model %d%%", percent);
		bu_vls_printf(&text, "  %s/%s parts",
		    ged_obol_lod_compact_count(
			std::min(status.availableLeafCount,
			    status.expectedLeafCount)).c_str(),
		    ged_obol_lod_compact_count(
			status.expectedLeafCount).c_str());
	    } else {
		bu_vls_sprintf(&text, "Discovering model");
		if (status.availableLeafCount > 0)
		    bu_vls_printf(&text, "  %s parts found",
			ged_obol_lod_compact_count(
			    status.availableLeafCount).c_str());
	    }
	    break;
	case BOBOL_LOD_CONVERGENCE_INTERACTIVE:
	    color[0] = 255;
	    color[1] = 205;
	    color[2] = 72;
	    bu_vls_sprintf(&text, "Interactive detail  %s triangles",
		ged_obol_lod_compact_count(status.activeFaces).c_str());
	    break;
	case BOBOL_LOD_CONVERGENCE_CALIBRATING:
	    color[0] = 255;
	    color[1] = 170;
	    color[2] = 64;
	    bu_vls_sprintf(&text,
		"Balancing detail for %.0f FPS  %s triangles",
		static_cast<double>(controller->getLodStableTargetFps()),
		ged_obol_lod_compact_count(status.activeFaces).c_str());
	    break;
	case BOBOL_LOD_CONVERGENCE_REFINING:
	    bu_vls_sprintf(&text, "Refining view %d%%", percent);
	    if (status.visibleTargetCount > 0)
		bu_vls_printf(&text, "  %s/%s targets",
		    ged_obol_lod_compact_count(
			std::min(status.satisfiedPayloadCount,
			    status.visibleTargetCount)).c_str(),
		    ged_obol_lod_compact_count(
			status.visibleTargetCount).c_str());
	    if (pending > 0)
		bu_vls_printf(&text, "  %s loading",
		    ged_obol_lod_compact_count(pending).c_str());
	    break;
	case BOBOL_LOD_CONVERGENCE_BACKGROUND:
	    color[0] = 112;
	    color[1] = 235;
	    color[2] = 135;
	    bu_vls_sprintf(&text, "View ready  optimizing memory/cache");
	    if (status.residentMeshBytes > 0)
		bu_vls_printf(&text, "  %.0f MB resident",
		    static_cast<double>(status.residentMeshBytes) /
			(1024.0 * 1024.0));
	    if (status.gpuTrackedBufferBytes > 0)
		bu_vls_printf(&text, "  %.0f MB GPU buffers",
		    static_cast<double>(status.gpuTrackedBufferBytes) /
			(1024.0 * 1024.0));
	    break;
	case BOBOL_LOD_CONVERGENCE_ERROR:
	    color[0] = 255;
	    color[1] = 90;
	    color[2] = 80;
	    bu_vls_sprintf(&text, "View incomplete  %u geometry error%s",
		status.failedSourceCount,
		status.failedSourceCount == 1 ? "" : "s");
	    break;
	case BOBOL_LOD_CONVERGENCE_IDLE:
	default:
	    color[0] = 112;
	    color[1] = 235;
	    color[2] = 135;
	    if (status.gpuMemoryPressure) {
		color[0] = 255;
		color[1] = 170;
		color[2] = 64;
		bu_vls_sprintf(&text,
		    "View ready  GPU-memory-limited  %s triangles",
		    ged_obol_lod_compact_count(status.activeFaces).c_str());
		if (status.gpuPressureProxyCount > 0)
		    bu_vls_printf(&text, "  %s pressure proxies",
			ged_obol_lod_compact_count(
			    status.gpuPressureProxyCount).c_str());
	    } else {
		bu_vls_sprintf(&text,
		    status.performanceLimited ?
			"View ready  responsiveness-limited  %s triangles" :
			"View ready  %s triangles",
		    ged_obol_lod_compact_count(status.activeFaces).c_str());
	    }
	    break;
    }

    const int track_color[3] = {72, 78, 86};
    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(2);
    commands.reserve(2);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
	0.965, -0.58, 0.965, 0.58);
    BObolFeatureStyle track_style = ged_obol_faceplate_style(
	track_color, 72, 78, 86, 9);
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx,
	track_name, points, commands, track_style);

    points.clear();
    commands.clear();
    const fastf_t fill_top =
	-0.58 + 1.16 * std::max(0.0f,
	    std::min(1.0f, status.fraction));
    if (fill_top > -0.58) {
	ged_obol_faceplate_append_line(points, commands, view_ctx,
	    0.965, -0.58, 0.965, fill_top);
	BObolFeatureStyle fill_style = ged_obol_faceplate_style(
	    color, 96, 220, 255, 7);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
	    fill_name, points, commands, fill_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, fill_name);
    }

    std::vector<BObolLabel> labels;
    labels.push_back(ged_obol_faceplate_label(view_ctx,
	bu_vls_cstr(&text), 0.10, -0.70, color,
	96, 220, 255, 12, 0));
    BObolFeatureStyle label_style = ged_obol_faceplate_style(
	color, 96, 220, 255, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
	label_name, labels, label_style);
    bu_vls_free(&text);
}

static void
ged_obol_faceplate_sync_scale(BObolViewController *controller,
			      struct ged_view_context *view_ctx)
{
    static const char line_name[] = "_faceplate/scale";
    static const char label_name[] = "_faceplate/scale_labels";
    struct bv_other_state state = BV_OTHER_STATE_INIT;
    if (!bv_scale_overlay_state_get(&state, ged_obol_bv_const(view_ctx)) ||
	!state.gos_draw) {
	ged_obol_faceplate_remove(controller, view_ctx, line_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name);
	return;
    }

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    points.reserve(6);
    commands.reserve(6);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   -0.5, -0.8, 0.5, -0.8);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   -0.5, -0.79, -0.5, -0.81);
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   0.5, -0.79, 0.5, -0.81);

    BObolFeatureStyle line_style = ged_obol_faceplate_style(
					 state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_lines(controller, view_ctx, line_name,
					   points, commands, line_style);

    struct bu_vls scale = BU_VLS_INIT_ZERO;
    const fastf_t base2local = bv_base2local_get(ged_obol_bv_const(view_ctx));
    const char *unit = !ZERO(base2local) ? bu_units_string(1.0 / base2local) :
		       NULL;
    if (!unit)
	unit = "";
    bu_vls_printf(&scale, "%g%s",
		  bv_size_get(ged_obol_bv_const(view_ctx)) * 0.5 *
		  base2local,
		  unit);
    const int soffset = (int)(strlen(bu_vls_cstr(&scale)) * 0.5);
    std::vector<BObolLabel> labels;
    labels.push_back(ged_obol_faceplate_label(view_ctx, "0", -0.505, -0.78,
		     state.gos_line_color, 255, 255, 0,
		     state.gos_font_size > 0 ? state.gos_font_size : 20));
    labels.push_back(ged_obol_faceplate_label(view_ctx, bu_vls_cstr(&scale),
		     0.5 - (soffset * 0.015), -0.78,
		     state.gos_line_color, 255, 255, 0,
		     state.gos_font_size > 0 ? state.gos_font_size : 20));
    bu_vls_free(&scale);

    BObolFeatureStyle text_style = ged_obol_faceplate_style(
					 state.gos_line_color, 255, 255, 0, 1);
    (void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
	    label_name, labels, text_style);
}

static void
ged_obol_faceplate_append_axis(std::vector<SbVec3f> &points,
			       std::vector<int32_t> &commands,
			       std::vector<BObolLabel> &labels,
			       struct ged_view_context *view_ctx,
			       const mat_t rmat,
			       const struct bv_axes_state *axes,
			       int axis,
			       fastf_t aspect,
			       const char *label_text)
{
    point_t v2 = VINIT_ZERO;
    v2[axis] = axes->axes_size > 0.0 ? axes->axes_size * 0.5 : 0.1;

    point_t rv2;
    point_t rv1;
    MAT4X3PNT(rv2, rmat, v2);
    if (axes->pos_only) {
	VSET(rv1, 0.0, 0.0, 0.0);
    } else {
	VSCALE(rv1, rv2, -1.0);
    }

    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   rv1[X] + axes->axes_pos[X],
				   (rv1[Y] + axes->axes_pos[Y]) * aspect,
				   rv2[X] + axes->axes_pos[X],
				   (rv2[Y] + axes->axes_pos[Y]) * aspect);

    if (axes->label_flag) {
	point_t lv;
	point_t lrv;
	VSET(lv, v2[X] + 0.0078125, v2[Y] + 0.0078125,
	     v2[Z] + 0.0078125);
	MAT4X3PNT(lrv, rmat, lv);
	labels.push_back(ged_obol_faceplate_label(view_ctx, label_text,
			 lrv[X] + axes->axes_pos[X],
			 lrv[Y] + axes->axes_pos[Y],
			 axes->label_color, 255, 255, 0, 20));
    }
}

static void
ged_obol_faceplate_axis_triple_color(int axis, int rgb[3])
{
    if (!rgb)
	return;

    switch (axis) {
	case X:
	    VSET(rgb, 255, 0, 0);
	    break;
	case Y:
	    VSET(rgb, 0, 255, 0);
	    break;
	case Z:
	    VSET(rgb, 0, 0, 255);
	    break;
	default:
	    VSET(rgb, 255, 255, 255);
	    break;
    }
}

static const char *
ged_obol_faceplate_axis_suffix(int axis)
{
    switch (axis) {
	case X:
	    return "/x";
	case Y:
	    return "/y";
	case Z:
	    return "/z";
	default:
	    return "/axis";
    }
}

static void
ged_obol_faceplate_remove_axis_variants(BObolViewController *controller,
					struct ged_view_context *view_ctx,
					const std::string &line_name)
{
    ged_obol_faceplate_remove(controller, view_ctx, line_name.c_str());
    for (int axis = X; axis <= Z; axis++) {
	std::string axis_name = line_name + ged_obol_faceplate_axis_suffix(axis);
	ged_obol_faceplate_remove(controller, view_ctx, axis_name.c_str());
    }
}

static void
ged_obol_faceplate_append_tick_segment(std::vector<SbVec3f> &points,
				       std::vector<int32_t> &commands,
				       struct ged_view_context *view_ctx,
				       const fastf_t axes_pos[3],
				       const point_t t1,
				       const point_t t2,
				       fastf_t aspect)
{
    ged_obol_faceplate_append_line(points, commands, view_ctx,
				   t1[X] + axes_pos[X], (t1[Y] + axes_pos[Y]) * aspect,
				   t2[X] + axes_pos[X], (t2[Y] + axes_pos[Y]) * aspect);
}

static void
ged_obol_faceplate_append_axis_ticks(std::vector<SbVec3f> &tick_points,
				     std::vector<int32_t> &tick_commands,
				     std::vector<SbVec3f> &major_points,
				     std::vector<int32_t> &major_commands,
				     struct ged_view_context *view_ctx,
				     const mat_t rmat,
				     const struct bv_axes_state *axes,
				     fastf_t aspect)
{
    if (!view_ctx || !axes || !axes->tick_enabled ||
	axes->tick_interval <= SMALL_FASTF)
	return;

    const fastf_t view_size = bv_size_get(ged_obol_bv_const(view_ctx));
    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const fastf_t half_axes_size =
	axes->axes_size > 0.0 ? axes->axes_size * 0.5 : 0.1;
    if (view_size <= SMALL_FASTF || width <= 0 ||
	half_axes_size <= SMALL_FASTF)
	return;

    int num_ticks = static_cast<int>(
			view_size / axes->tick_interval * 0.5 * half_axes_size);
    if (num_ticks <= 0)
	return;

    int do_major_only = 0;
    if (axes->tick_threshold > 0 &&
	width <= num_ticks / half_axes_size * axes->tick_threshold * 2) {
	const int ticks_per_major = axes->ticks_per_major > 0 ?
				    axes->ticks_per_major : 1;
	const int num_major_ticks = num_ticks / ticks_per_major;
	if (width <= num_major_ticks / half_axes_size *
	    axes->tick_threshold * 2)
	    return;
	do_major_only = 1;
    }

    const fastf_t interval = axes->tick_interval / view_size * 2.0;
    const fastf_t tlen = axes->tick_length / (fastf_t)width * 2.0;
    const fastf_t maj_tlen =
	axes->tick_major_length / (fastf_t)width * 2.0;

    vect_t xend1 = VINIT_ZERO;
    vect_t xend2 = VINIT_ZERO;
    vect_t yend1 = VINIT_ZERO;
    vect_t yend2 = VINIT_ZERO;
    vect_t zend1 = VINIT_ZERO;
    vect_t zend2 = VINIT_ZERO;
    vect_t maj_xend1 = VINIT_ZERO;
    vect_t maj_xend2 = VINIT_ZERO;
    vect_t maj_yend1 = VINIT_ZERO;
    vect_t maj_yend2 = VINIT_ZERO;
    vect_t maj_zend1 = VINIT_ZERO;
    vect_t maj_zend2 = VINIT_ZERO;
    vect_t rxdir = VINIT_ZERO;
    vect_t neg_rxdir = VINIT_ZERO;
    vect_t rydir = VINIT_ZERO;
    vect_t neg_rydir = VINIT_ZERO;
    vect_t rzdir = VINIT_ZERO;
    vect_t neg_rzdir = VINIT_ZERO;
    vect_t dir = VINIT_ZERO;

    if (!do_major_only) {
	VSET(dir, tlen, 0.0, 0.0);
	MAT4X3PNT(xend1, rmat, dir);
	VSCALE(xend2, xend1, -1.0);
	VSET(dir, 0.0, tlen, 0.0);
	MAT4X3PNT(yend1, rmat, dir);
	VSCALE(yend2, yend1, -1.0);
	VSET(dir, 0.0, 0.0, tlen);
	MAT4X3PNT(zend1, rmat, dir);
	VSCALE(zend2, zend1, -1.0);
    }

    VSET(dir, maj_tlen, 0.0, 0.0);
    MAT4X3PNT(maj_xend1, rmat, dir);
    VSCALE(maj_xend2, maj_xend1, -1.0);
    VSET(dir, 0.0, maj_tlen, 0.0);
    MAT4X3PNT(maj_yend1, rmat, dir);
    VSCALE(maj_yend2, maj_yend1, -1.0);
    VSET(dir, 0.0, 0.0, maj_tlen);
    MAT4X3PNT(maj_zend1, rmat, dir);
    VSCALE(maj_zend2, maj_zend1, -1.0);

    VSET(dir, interval, 0.0, 0.0);
    MAT4X3PNT(rxdir, rmat, dir);
    VSCALE(neg_rxdir, rxdir, -1.0);
    VSET(dir, 0.0, interval, 0.0);
    MAT4X3PNT(rydir, rmat, dir);
    VSCALE(neg_rydir, rydir, -1.0);
    VSET(dir, 0.0, 0.0, interval);
    MAT4X3PNT(rzdir, rmat, dir);
    VSCALE(neg_rzdir, rzdir, -1.0);

    auto append_tick_pair = [&](const vect_t e1, const vect_t e2,
    const vect_t along, int major) {
	point_t t1;
	point_t t2;
	VADD2(t1, e1, along);
	VADD2(t2, e2, along);
	if (major)
	    ged_obol_faceplate_append_tick_segment(major_points,
						   major_commands, view_ctx, axes->axes_pos, t1, t2,
						   aspect);
	else
	    ged_obol_faceplate_append_tick_segment(tick_points,
						   tick_commands, view_ctx, axes->axes_pos, t1, t2,
						   aspect);
    };

    for (int i = 1; i <= num_ticks; i++) {
	const int major = axes->ticks_per_major > 0 ?
			  (i % axes->ticks_per_major == 0) : 0;
	if (!major && do_major_only)
	    continue;

	const vect_t *x1 = major ? &maj_xend1 : &xend1;
	const vect_t *x2 = major ? &maj_xend2 : &xend2;
	const vect_t *y1 = major ? &maj_yend1 : &yend1;
	const vect_t *y2 = major ? &maj_yend2 : &yend2;
	const vect_t *z1 = major ? &maj_zend1 : &zend1;
	const vect_t *z2 = major ? &maj_zend2 : &zend2;

	vect_t tvec;
	VSCALE(tvec, rxdir, i);
	append_tick_pair(*y1, *y2, tvec, major);
	append_tick_pair(*z1, *z2, tvec, major);
	if (!axes->pos_only) {
	    VSCALE(tvec, neg_rxdir, i);
	    append_tick_pair(*y1, *y2, tvec, major);
	    append_tick_pair(*z1, *z2, tvec, major);
	}

	VSCALE(tvec, rydir, i);
	append_tick_pair(*x1, *x2, tvec, major);
	append_tick_pair(*z1, *z2, tvec, major);
	if (!axes->pos_only) {
	    VSCALE(tvec, neg_rydir, i);
	    append_tick_pair(*x1, *x2, tvec, major);
	    append_tick_pair(*z1, *z2, tvec, major);
	}

	VSCALE(tvec, rzdir, i);
	append_tick_pair(*x1, *x2, tvec, major);
	append_tick_pair(*y1, *y2, tvec, major);
	if (!axes->pos_only) {
	    VSCALE(tvec, neg_rzdir, i);
	    append_tick_pair(*x1, *x2, tvec, major);
	    append_tick_pair(*y1, *y2, tvec, major);
	}
    }
}

static void
ged_obol_faceplate_sync_axes_one(BObolViewController *controller,
				 struct ged_view_context *view_ctx,
				 const char *prefix,
				 struct bv_axes_state axes,
				 int position_mode,
				 const mat_t supplied_rotation = NULL)
{
    const int model_axes = position_mode == 2;
    std::string line_name = std::string(prefix) + "/lines";
    std::string label_name = std::string(prefix) + "/labels";
    std::string tick_name = std::string(prefix) + "/ticks";
    std::string major_tick_name = std::string(prefix) + "/major_ticks";
    if (!axes.draw) {
	ged_obol_faceplate_remove_axis_variants(controller, view_ctx, line_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
	return;
    }

    if (axes.axes_size <= 0.0)
	axes.axes_size = model_axes ? 2.0 : 0.2;
    if (ged_obol_rgb_is_zero(axes.axes_color))
	VSET(axes.axes_color, 255, 255, 255);
    if (ged_obol_rgb_is_zero(axes.label_color))
	VSET(axes.label_color, 255, 255, 0);
    if (position_mode == 1 && VNEAR_ZERO(axes.axes_pos, SMALL_FASTF)) {
	VSET(axes.axes_pos, 0.80, -0.80, 0.0);
	axes.pos_only = 1;
	axes.triple_color = 1;
	axes.label_flag = 1;
    }
    if (model_axes)
	axes.label_flag = 1;

    mat_t rmat;
    if (supplied_rotation) {
	MAT_COPY(rmat, supplied_rotation);
    } else if (!bv_rotation_get(rmat, ged_obol_bv_const(view_ctx))) {
	ged_obol_faceplate_remove_axis_variants(controller, view_ctx, line_name);
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
	return;
    }
    if (model_axes) {
	point_t map;
	mat_t model2view;
	if (bv_model2view_get(model2view, ged_obol_bv_const(view_ctx))) {
	    VSCALE(map, axes.axes_pos,
		   bv_local2base_get(ged_obol_bv_const(view_ctx)));
	    MAT4X3PNT(axes.axes_pos, model2view, map);
	}
    }

    const int width = bv_width_get(ged_obol_bv_const(view_ctx));
    const int height = bv_height_get(ged_obol_bv_const(view_ctx));
    const fastf_t aspect = (width > 0 && height > 0) ?
			   (fastf_t)width / (fastf_t)height : 1.0;
    if (!model_axes)
	axes.axes_pos[Y] *= (aspect > SMALL_FASTF) ? 1.0 / aspect : 1.0;

    std::vector<SbVec3f> points;
    std::vector<int32_t> commands;
    std::vector<BObolLabel> labels;
    points.reserve(6);
    commands.reserve(6);
    labels.reserve(3);

    if (axes.triple_color) {
	ged_obol_faceplate_remove(controller, view_ctx, line_name.c_str());
	for (int axis = X; axis <= Z; axis++) {
	    struct bv_axes_state axis_axes = axes;
	    int axis_color[3] = {0, 0, 0};
	    ged_obol_faceplate_axis_triple_color(axis, axis_color);
	    VMOVE(axis_axes.axes_color, axis_color);
	    VMOVE(axis_axes.label_color, axis_color);

	    std::vector<SbVec3f> axis_points;
	    std::vector<int32_t> axis_commands;
	    std::vector<BObolLabel> axis_labels;
	    axis_points.reserve(2);
	    axis_commands.reserve(2);
	    ged_obol_faceplate_append_axis(axis_points, axis_commands,
					   axis_labels, view_ctx, rmat, &axis_axes, axis, aspect,
					   axis == X ? "X" : axis == Y ? "Y" : "Z");

	    BObolFeatureStyle axis_style = ged_obol_faceplate_style(
						 axis_axes.axes_color, 255, 255, 255,
						 axis_axes.line_width);
	    std::string axis_name =
		line_name + ged_obol_faceplate_axis_suffix(axis);
	    (void)ged_obol_faceplate_publish_lines(controller, view_ctx,
						   axis_name.c_str(), axis_points, axis_commands,
						   axis_style);
	    labels.insert(labels.end(), axis_labels.begin(),
			  axis_labels.end());
	}
    } else {
	for (int axis = X; axis <= Z; axis++) {
	    std::string axis_name =
		line_name + ged_obol_faceplate_axis_suffix(axis);
	    ged_obol_faceplate_remove(controller, view_ctx,
				      axis_name.c_str());
	}
	ged_obol_faceplate_append_axis(points, commands, labels, view_ctx,
				       rmat, &axes, X, aspect, "X");
	ged_obol_faceplate_append_axis(points, commands, labels, view_ctx,
				       rmat, &axes, Y, aspect, "Y");
	ged_obol_faceplate_append_axis(points, commands, labels, view_ctx,
				       rmat, &axes, Z, aspect, "Z");

	BObolFeatureStyle line_style = ged_obol_faceplate_style(
					     axes.axes_color, 255, 255, 255, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       line_name.c_str(), points, commands, line_style);
    }

    std::vector<SbVec3f> tick_points;
    std::vector<int32_t> tick_commands;
    std::vector<SbVec3f> major_tick_points;
    std::vector<int32_t> major_tick_commands;
    ged_obol_faceplate_append_axis_ticks(tick_points, tick_commands,
					 major_tick_points, major_tick_commands, view_ctx, rmat, &axes,
					 aspect);
    if (!tick_points.empty()) {
	BObolFeatureStyle tick_style = ged_obol_faceplate_style(
					     axes.tick_color, 255, 255, 0, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       tick_name.c_str(), tick_points, tick_commands, tick_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, tick_name.c_str());
    }
    if (!major_tick_points.empty()) {
	BObolFeatureStyle major_tick_style = ged_obol_faceplate_style(
		axes.tick_major_color, 255, 0, 0, axes.line_width);
	(void)ged_obol_faceplate_publish_lines(controller, view_ctx,
					       major_tick_name.c_str(), major_tick_points,
					       major_tick_commands, major_tick_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx,
				  major_tick_name.c_str());
    }

    if (!labels.empty()) {
	BObolFeatureStyle label_style = ged_obol_faceplate_style(
					      axes.label_color, 255, 255, 0, 1);
	(void)ged_obol_faceplate_publish_hud_labels(controller, view_ctx,
		label_name.c_str(), labels, label_style);
    } else {
	ged_obol_faceplate_remove(controller, view_ctx, label_name.c_str());
    }
}

static void
ged_obol_faceplate_sync_axes(BObolViewController *controller,
			     struct ged_view_context *view_ctx)
{
    struct bv_axes_state axes = BV_AXES_STATE_INIT;
    if (bv_view_axes_state_get(&axes, ged_obol_bv_const(view_ctx)))
	ged_obol_faceplate_sync_axes_one(controller, view_ctx,
					 "_faceplate/view_axes", axes, 1);
    if (bv_model_axes_state_get(&axes, ged_obol_bv_const(view_ctx)))
	ged_obol_faceplate_sync_axes_one(controller, view_ctx,
					 "_faceplate/model_axes", axes, 2);
}

extern "C" int
ged_draw_obol_view_context_hud_axes_replace(
    struct ged_view_context *view_ctx,
    const char *name,
    const struct bv_axes_state *axes,
    const mat_t rotation)
{
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller || !view_ctx || !name || !name[0])
	return 0;

    struct bv_axes_state state = BV_AXES_STATE_INIT;
    if (axes)
	state = *axes;
    ged_obol_faceplate_sync_axes_one(controller, view_ctx, name, state, 0,
	rotation);
    return 1;
}

static void
ged_obol_faceplate_sync_framebuffer(BObolViewController *controller,
			    struct ged_view_context *view_ctx)
{
    static const char name[] = "_faceplate/framebuffer";
    /* Framebuffer images are owned by BObolFramebufferStream.  Remove the
     * retired display-manager snapshot feature if an older host created it. */
    ged_obol_faceplate_remove(controller, view_ctx, name);
}

static int
ged_obol_view_context_faceplate_sync(struct ged *gedp, struct ged_view_context *view_ctx,
	BObolViewController *controller)
{
    if (!gedp || !view_ctx || !controller)
	return BRLCAD_OK;

    ged_obol_faceplate_sync_center_dot(controller, view_ctx);
	ged_obol_faceplate_sync_interactive_rect(controller, view_ctx);
    ged_obol_faceplate_sync_grid(controller, view_ctx);
    ged_obol_faceplate_sync_adc(controller, view_ctx);
    /* Controller render timing is useful performance telemetry, but it is not
     * an observed frame rate.  Hosts publish their completed-presentation
     * cadence separately, including offscreen/headless hosts. */
    const uint64_t presentation_interval =
	controller->getDisplayedPresentationIntervalNanoseconds();
    if (presentation_interval)
	(void)bv_frametime_set(ged_obol_bv(view_ctx), presentation_interval);
    ged_obol_faceplate_sync_params(controller, view_ctx);
    ged_obol_faceplate_sync_lod_progress(controller, view_ctx);
    ged_obol_faceplate_sync_scale(controller, view_ctx);
    ged_obol_faceplate_sync_axes(controller, view_ctx);
    ged_obol_faceplate_sync_framebuffer(controller, view_ctx);

	const int framebuffer_mode =
	bv_framebuffer_mode_get(ged_obol_bv_const(view_ctx));
    const int framebuffer_visible = framebuffer_mode ==
	BV_FRAMEBUFFER_MODE_OVERLAY || framebuffer_mode ==
	BV_FRAMEBUFFER_MODE_UNDERLAY || framebuffer_mode ==
	BV_FRAMEBUFFER_MODE_INTERLAY;
    (void)ged_obol_fbserv_composition_set(gedp, framebuffer_mode);
    if (framebuffer_visible)
	(void)ged_view_framebuffer_present(gedp);

    return BRLCAD_OK;
}

extern "C" int
ged_view_faceplate_sync(struct ged *gedp, struct ged_view_context *view_ctx)
{
    return ged_obol_view_context_faceplate_sync(gedp, view_ctx,
	ged_obol_view_controller_for_context(view_ctx));
}

/* First shader token == "light"? */
static bool
ged_obol_shader_is_light(const char *shader)
{
    if (!shader)
	return false;
    while (*shader == ' ' || *shader == '\t')
	shader++;
    if (bu_strncmp(shader, "light", 5) != 0)
	return false;
    const char c = shader[5];
    return (c == '\0' || c == ' ' || c == '\t' || c == ';');
}

/* Collect in-scene lights from the database's "light"-shader regions.  Done in
 * the GED layer (not the Obol geometry realize walk) so it is independent of the
 * LoD/geometry cache, which otherwise skips the walk on a warm cache.  Positions
 * are model/world-space (region bbox centers); parameters mirror sh_light.c. */
static void
ged_obol_collect_scene_lights(struct ged *gedp,
	std::vector<BObolSceneLightRealization> &out)
{
    out.clear();
    if (!gedp || !gedp->dbip)
	return;
    struct db_i *dbip = gedp->dbip;

    struct directory *dp = RT_DIR_NULL;
    FOR_ALL_DIRECTORY_START(dp, dbip) {
	    if ((dp->d_flags & RT_DIR_HIDDEN) || !(dp->d_flags & RT_DIR_COMB))
		continue;

	    struct rt_db_internal intern;
	    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
		continue;
	    struct rt_comb_internal *comb =
		(struct rt_comb_internal *)intern.idb_ptr;
	    if (!comb || comb->magic != RT_COMB_MAGIC ||
		!ged_obol_shader_is_light(bu_vls_cstr(&comb->shader))) {
		rt_db_free_internal(&intern);
		continue;
	    }

	    /* World-space light position: region geometry bbox center. */
	    struct bu_vls msgs = BU_VLS_INIT_ZERO;
	    point_t rpp_min = VINIT_ZERO;
	    point_t rpp_max = VINIT_ZERO;
	    const char *av[2] = {dp->d_namep, NULL};
	    if (rt_obj_bounds(&msgs, dbip, 1, av, 0, rpp_min, rpp_max) < 0) {
		bu_vls_free(&msgs);
		rt_db_free_internal(&intern);
		continue;
	    }
	    bu_vls_free(&msgs);
	    point_t center;
	    VADD2SCALE(center, rpp_min, rpp_max, 0.5);

	    /* Parameters (defaults mirror sh_light.c): intensity 1, angle 180
	     * (omni), not infinite. */
	    double intensity = 1.0;
	    double angle = 180.0;
	    double fraction = 1.0;
	    double target[3] = {0.0, 0.0, 0.0};
	    int infinite = 0;
	    int exaim = 0;
	    const char *shader = bu_vls_cstr(&comb->shader);
	    struct bu_vls argcopy = BU_VLS_INIT_ZERO;
	    bu_vls_strcpy(&argcopy, shader + 5); /* skip "light" */
	    for (char *s = bu_vls_addr(&argcopy); *s; s++) {
		if (*s == '=' || *s == ';' || *s == ',')
		    *s = ' ';
	    }
	    char *lav[64];
	    size_t lac = bu_argv_from_string(lav, 63, bu_vls_addr(&argcopy));
	    for (size_t k = 0; k < lac; k++) {
		const char *key = lav[k];
		if ((BU_STR_EQUAL(key, "bright") || BU_STR_EQUAL(key, "b") ||
			BU_STR_EQUAL(key, "inten")) && k + 1 < lac)
		    intensity = atof(lav[++k]);
		else if ((BU_STR_EQUAL(key, "angle") || BU_STR_EQUAL(key, "a")) &&
			k + 1 < lac)
		    angle = atof(lav[++k]);
		else if ((BU_STR_EQUAL(key, "fract") || BU_STR_EQUAL(key, "f")) &&
			k + 1 < lac)
		    fraction = atof(lav[++k]);
		else if ((BU_STR_EQUAL(key, "infinite") || BU_STR_EQUAL(key, "i")) &&
			k + 1 < lac)
		    infinite = atoi(lav[++k]);
		else if ((BU_STR_EQUAL(key, "target") || BU_STR_EQUAL(key, "t") ||
			BU_STR_EQUAL(key, "aim") || BU_STR_EQUAL(key, "d") ||
			BU_STR_EQUAL(key, "dir")) && k + 3 < lac) {
		    target[0] = atof(lav[++k]);
		    target[1] = atof(lav[++k]);
		    target[2] = atof(lav[++k]);
		    exaim = 1;
		}
	    }
	    bu_vls_free(&argcopy);

	    vect_t aim = {0.0, 0.0, -1.0};
	    if (exaim)
		VSUB2(aim, target, center);
	    VUNITIZE(aim);

	    BObolSceneLightRealization L;
	    L.name = dp->d_namep ? dp->d_namep : "";
	    if (comb->rgb_valid) {
		L.color = SbColor(comb->rgb[0] / 255.0f, comb->rgb[1] / 255.0f,
		    comb->rgb[2] / 255.0f);
	    } else {
		L.color = SbColor(1.0f, 1.0f, 1.0f);
	    }
	    double norm = intensity * (fraction >= 0.0 ? fraction : 1.0);
	    if (norm < 0.0)
		norm = 0.0;
	    if (norm > 1.0)
		norm = 1.0;
	    L.intensity = (float)norm;
	    L.direction = SbVec3f((float)aim[0], (float)aim[1], (float)aim[2]);
	    if (infinite) {
		L.kind = BOBOL_SCENE_LIGHT_DIRECTIONAL;
	    } else if (angle < 180.0) {
		L.kind = BOBOL_SCENE_LIGHT_SPOT;
		L.position = SbVec3f((float)center[0], (float)center[1],
		    (float)center[2]);
		L.coneAngleDeg = (float)angle;
	    } else {
		L.kind = BOBOL_SCENE_LIGHT_POINT;
		L.position = SbVec3f((float)center[0], (float)center[1],
		    (float)center[2]);
	    }
	    out.push_back(L);
	    rt_db_free_internal(&intern);
    } FOR_ALL_DIRECTORY_END;
}

extern "C" int
ged_view_lighting_sync(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    if (!view_ctx)
	return BRLCAD_ERROR;
    const struct bv *view = ged_obol_bv_const(view_ctx);
    if (!view)
	return BRLCAD_ERROR;
    struct bv_lighting_state lighting;
    if (!bv_lighting_state_get(&lighting, view))
	return BRLCAD_ERROR;

    /* bv_lighting_state is the persistent source of truth; push it to the live
     * Obol controller when one is attached (headless views just keep the bv
     * state).  Select the complete preset before applying its independently
     * adjustable primary direction. */
    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return BRLCAD_OK;

    controller->setLightingProfile(lighting.profile == BV_LIGHTING_MGED ?
	BObolViewController::LIGHTING_MGED :
	BObolViewController::LIGHTING_STUDIO);
    controller->setHeadlightOffset(SbVec3f(
	(float)lighting.headlight_offset[0],
	(float)lighting.headlight_offset[1],
	(float)lighting.headlight_offset[2]));
    controller->setHeadlightCameraTracked(
	lighting.headlight_tracks_camera ? TRUE : FALSE);
    controller->setHeadlightEnabled(
	lighting.headlight_enabled ? TRUE : FALSE);
    /* Supply the in-scene lights from the database (cache-immune) and apply the
     * enable state. */
    std::vector<BObolSceneLightRealization> scene_lights;
    ged_obol_collect_scene_lights(gedp, scene_lights);
    controller->setSceneLights(scene_lights);
    controller->setSceneLightsEnabled(
	lighting.scene_lights_enabled ? TRUE : FALSE);
    return BRLCAD_OK;
}

extern "C" int
ged_view_shading_sync(struct ged *gedp,
	struct ged_view_context *view_ctx)
{
    (void)gedp;
    if (!view_ctx)
	return BRLCAD_ERROR;
    const struct bv *view = ged_obol_bv_const(view_ctx);
    if (!view)
	return BRLCAD_ERROR;
    struct bv_shading_state shading;
    if (!bv_shading_state_get(&shading, view))
	return BRLCAD_ERROR;

    BObolViewController *controller =
	ged_obol_view_controller_for_context(view_ctx);
    if (!controller)
	return BRLCAD_OK;

    BObolViewLodState::NormalStyle style =
	BObolViewLodState::NORMAL_AUTHORED;
    if (shading.normal_style == BV_NORMAL_FLAT)
	style = BObolViewLodState::NORMAL_FLAT;
    else if (shading.normal_style == BV_NORMAL_SMOOTH)
	style = BObolViewLodState::NORMAL_SMOOTH;
    controller->setNormalStyle(style,
	static_cast<float>(shading.normal_crease_angle));
    return BRLCAD_OK;
}

extern "C" int
ged_draw_obol_view_context_feature_store_active(struct ged_view_context *view_ctx)
{
    return ged_obol_view_controller_for_context(view_ctx) ? 1 : 0;
}

extern "C" size_t
ged_draw_obol_view_context_clear(struct ged_view_context *view_ctx, int flags)
{
    if (!(flags & GED_VIEW_CLEAR_VIEW))
	return 0;

    size_t removed = 0;
    if (flags & GED_VIEW_CLEAR_LOCAL) {
	BObolViewController *controller =
	    ged_obol_view_controller_for_context(view_ctx);
	if (!controller)
	    return 0;
	BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
	removed += controller->features().removeScope(
		       BOBOL_FEATURE_SCOPE_LOCAL, &owner);
	removed += controller->polygons().removeScope(
		       BOBOL_FEATURE_SCOPE_LOCAL);
	controller->selection().clear(&owner, BOBOL_SELECTION_ALL);
    } else {
	BObolViewController *controller =
	    ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
	if (!controller)
	    return 0;
	removed += controller->features().removeScope(
		       BOBOL_FEATURE_SCOPE_SHARED, NULL);
	removed += controller->polygons().removeScope(
		       BOBOL_FEATURE_SCOPE_SHARED);
	controller->selection().clear(NULL, BOBOL_SELECTION_ALL);
    }
    return removed;
}

static enum ged_draw_obol_view_feature_kind
ged_obol_view_feature_kind(BObolFeatureKind kind) {
    switch (kind)
    {
	case BObolFeatureKind::Lines:
		    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES;
	case BObolFeatureKind::IndexedLines:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_LINES;
	case BObolFeatureKind::Points:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_POINTS;
	case BObolFeatureKind::Labels:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LABELS;
	case BObolFeatureKind::Arrow:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_ARROW;
	case BObolFeatureKind::Axes:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_AXES;
	case BObolFeatureKind::LineLayer:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINE_LAYER;
	case BObolFeatureKind::EditPreview:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_EDIT_PREVIEW;
	case BObolFeatureKind::IndexedFaceSet:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_INDEXED_FACE_SET;
	case BObolFeatureKind::PolygonOverlay:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_POLYGON_OVERLAY;
	case BObolFeatureKind::HudLabel:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_HUD_LABEL;
	case BObolFeatureKind::CustomNode:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_CUSTOM_NODE;
	case BObolFeatureKind::Unknown:
	default:
	    return GED_DRAW_OBOL_VIEW_FEATURE_KIND_UNKNOWN;
    }
}

struct ged_obol_feature_records_visit {
    struct ged_view_context *view_ctx;
    unsigned int query_flags;
    const char *glob;
    ged_draw_obol_view_feature_record_cb cb;
    void *userdata;
    int count;
};

static BObolFeatureStyle
ged_obol_line_layer_effective_style(const BObolFeatureStyle &parent,
				    const BObolFeatureStyle &layer)
{
    BObolFeatureStyle out = parent;
    if (layer.hasVisible) {
	out.hasVisible = TRUE;
	out.visible = layer.visible;
    }
    if (layer.hasSelectable) {
	out.hasSelectable = TRUE;
	out.selectable = layer.selectable;
    }
    if (layer.hasColor) {
	out.hasColor = TRUE;
	out.color = layer.color;
    }
    if (layer.hasLineWidth) {
	out.hasLineWidth = TRUE;
	out.lineWidth = layer.lineWidth;
    }
    if (layer.hasLineStyle) {
	out.hasLineStyle = TRUE;
	out.lineStyle = layer.lineStyle;
    }
    if (layer.hasArrow) {
	out.hasArrow = TRUE;
	out.arrow = layer.arrow;
    }
    if (layer.hasArrowTip) {
	out.hasArrowTip = TRUE;
	out.arrowTipLength = layer.arrowTipLength;
	out.arrowTipWidth = layer.arrowTipWidth;
    }
    return out;
}

static int
ged_obol_feature_record_visible(const BObolFeatureStyle &parent,
				const BObolFeatureStyle *child)
{
    if (parent.hasVisible && !parent.visible)
	return 0;
    if (child && child->hasVisible && !child->visible)
	return 0;
    return 1;
}

static int
ged_obol_feature_record_glob_matches(const char *glob, const SbString &name)
{
    if (!glob || !glob[0])
	return 1;

    const char *str = name.getString();
    return (str && bu_path_match(glob, str, 0) == 0) ? 1 : 0;
}

static int
ged_obol_feature_record_emit(struct ged_obol_feature_records_visit *ctx,
			     const BObolFeatureRecord &record,
			     const SbString &name,
			     enum ged_draw_obol_view_feature_kind kind,
			     const BObolFeatureStyle &style,
			     int visible,
			     size_t point_count,
			     size_t command_count,
			     size_t index_count,
			     size_t normal_count,
			     size_t label_count,
			     size_t axes_center_count,
			     size_t child_count,
			     const char *line_layer_parent_name,
			     size_t line_layer_index)
{
    if (!ctx || !ctx->cb)
	return 0;

    if ((ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VISIBLE_ONLY) &&
	!visible)
	return 1;

    if (!ged_obol_feature_record_glob_matches(ctx->glob, name))
	return 1;

    struct ged_draw_obol_view_feature_record out;
    memset(&out, 0, sizeof(out));
    out.name = name.getString();
    out.kind = kind;
    out.local = record.scope == BObolFeatureScope::Local ? 1 : 0;
    out.visible = visible;
    out.realized = record.realized ? 1 : 0;
    out.color[0] = 255;
    out.color[1] = 255;
    out.color[2] = 255;
    if (style.hasColor)
	ged_obol_rgb_from_color(style.color, out.color);
    out.line_style = style.hasLineStyle ? style.lineStyle : 0;
    out.line_width = style.hasLineWidth ? style.lineWidth : 1;
    out.point_count = point_count;
    out.command_count = command_count;
    out.index_count = index_count;
    out.normal_count = normal_count;
    out.label_count = label_count;
    out.axes_center_count = axes_center_count;
    out.child_count = child_count;
    out.line_layer_parent_name = line_layer_parent_name;
    out.line_layer_index = line_layer_index;

    ctx->count++;
    return ctx->cb(&out, ctx->userdata);
}

static int
ged_obol_feature_record_visit_cb(const BObolFeatureRecord &record,
				 void *userData)
{
    struct ged_obol_feature_records_visit *ctx =
	(struct ged_obol_feature_records_visit *)userData;
    if (!ctx || !ctx->cb)
	return 0;

    const int wants_db =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(ctx->query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wants_db && !wants_view)
	return 1;

    if (record.kind == BObolFeatureKind::LineLayer &&
	!record.layers.empty()) {
	for (size_t i = 0; i < record.layers.size(); i++) {
	    const BObolLineLayer &layer = record.layers[i];
	    BObolFeatureStyle style =
		ged_obol_line_layer_effective_style(record.style,
						    layer.style);
	    int visible = ged_obol_feature_record_visible(record.style,
			  &layer.style);
	    size_t child_count = layer.points.empty() ? 0 : 1;
	    if (!ged_obol_feature_record_emit(ctx, record,
					      layer.name.getLength() ? layer.name : record.name,
					      GED_DRAW_OBOL_VIEW_FEATURE_KIND_LINES, style, visible,
					      layer.points.size(), layer.commands.size(), 0, 0, 0, 0,
					      child_count, record.name.getString(), i))
		return 0;
	}
	return 1;
    }

    size_t child_count = 0;
    if (!record.layers.empty())
	child_count = record.layers.size();
    else if (!record.labels.empty())
	child_count = record.labels.size();
    else if (!record.axesCenters.empty())
	child_count = record.axesCenters.size();
    else if (!record.points.empty())
	child_count = 1;
    else if (record.kind == BObolFeatureKind::CustomNode &&
	     record.realized)
	child_count = 1;

    return ged_obol_feature_record_emit(ctx, record, record.name,
					ged_obol_view_feature_kind(record.kind), record.style,
					ged_obol_feature_record_visible(record.style, NULL),
					record.points.size(), record.commands.size(),
					record.indices.size(), record.normals.size(), record.labels.size(),
					record.axesCenters.size(), child_count, NULL, 0);
}

extern "C" int
ged_draw_obol_view_context_feature_records_foreach(
    struct ged_view_context *view_ctx,
    unsigned int query_flags,
    const char *glob,
    ged_draw_obol_view_feature_record_cb cb,
    void *userdata)
{
    BObolViewController *local_controller =
	ged_obol_view_controller_for_context(view_ctx);
    BObolViewController *shared_controller =
	ged_obol_shared_view_controller_ensure_for_context(view_ctx, 0);
    if ((!local_controller && !shared_controller) || !cb)
	return 0;

    const int wants_db =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_DB_OBJECTS) ? 1 : 0;
    const int wants_view =
	(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_VIEW_OBJECTS) ? 1 : 0;
    if (wants_db && !wants_view)
	return 0;

    BObolFeatureOwner owner = ged_obol_feature_owner(view_ctx, 1);
    struct ged_obol_feature_records_visit ctx;
    ctx.view_ctx = view_ctx;
    ctx.query_flags = query_flags;
    ctx.glob = glob;
    ctx.cb = cb;
    ctx.userdata = userdata;
    ctx.count = 0;
    if (local_controller) {
	local_controller->features().visitRecords(
		ged_obol_feature_record_visit_cb, &ctx,
		BOBOL_FEATURE_SCOPE_LOCAL, &owner);
    }
    if (!(query_flags & GED_DRAW_VIEW_EXPORT_QUERY_LOCAL_ONLY) &&
	shared_controller) {
	shared_controller->features().visitRecords(
		ged_obol_feature_record_visit_cb, &ctx,
		BOBOL_FEATURE_SCOPE_SHARED, NULL);
    }
    return ctx.count;
}

static BObolViewController *
ged_obol_feature_controller_from_ged_ref(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    if (!view_ctx || !ref.id || !ref.generation)
	return NULL;

    const int local = (ref.owner & 1u) ? 1 : 0;
    if (ref.owner != ged_view_context_reference_owner(view_ctx, local))
	return NULL;

    BObolViewController *controller = local ?
	ged_obol_view_controller_for_context(view_ctx) :
	ged_obol_shared_view_controller_for_context(view_ctx);
    if (!controller ||
	controller->features().referenceGeneration() != ref.generation)
	return NULL;

    BObolFeatureRecord record;
    if (!controller->features().record(BObolFeatureHandle(ref.id, 1),
	    record))
	return NULL;
    return controller;
}

static BObolFeatureHandle
ged_obol_feature_handle_from_ged_ref(ged_view_edit_ref ref)
{
    return (ref.id && ref.generation) ? BObolFeatureHandle(ref.id, 1) :
	   BObolFeatureHandle();
}

static ged_view_edit_ref
ged_obol_ged_feature_ref(struct ged_view_context *view_ctx,
			 BObolViewController *controller,
			 BObolFeatureHandle handle,
			 int local)
{
    ged_view_edit_ref ref = GED_VIEW_EDIT_REF_NULL_INIT;
    if (!controller || !handle.isValid())
	return ref;
    ref.owner = ged_view_context_reference_owner(view_ctx, local);
    ref.id = handle.id;
    ref.generation = controller->features().referenceGeneration();
    if (!ref.owner)
	return GED_VIEW_EDIT_REF_NULL;
    return ref;
}

static BObolOverlayInfo
ged_obol_edit_overlay_info(struct ged_view_context *view_ctx,
			      const void *owner,
			      const char *source_path,
			      int sort_order)
{
    BObolOverlayInfo overlay = ged_obol_model_overlay_info(view_ctx,
				 BObolOverlayClass::EditHandle,
				 BObolOverlayLifecycle::PerTool,
				 BObolOverlayOrder::PostTransparent,
				 source_path);
    overlay.ownerToken = owner ? owner : view_ctx;
    overlay.sortOrder = sort_order;
    return overlay;
}

static std::vector<BObolLabel>
ged_obol_labels_from_ged_feature(
    const struct ged_view_feature_label *labels,
    size_t label_count)
{
    std::vector<BObolLabel> out;
    if (!labels || !label_count)
	return out;

    out.reserve(label_count);
    for (size_t i = 0; i < label_count; i++) {
	BObolLabel label;
	label.text = labels[i].text ? labels[i].text : "";
	label.point = SbVec3f(
			  static_cast<float>(labels[i].point[X]),
			  static_cast<float>(labels[i].point[Y]),
			  static_cast<float>(labels[i].point[Z]));
	if (labels[i].color_valid) {
	    label.hasColor = TRUE;
	    label.color = ged_obol_color_from_rgb(labels[i].color);
	}
	if (labels[i].font_size > 0.0)
	    label.fontSize = static_cast<float>(labels[i].font_size);
	out.push_back(label);
    }
    return out;
}

static unsigned char ged_obol_rgb_byte_from_int(int value);

extern "C" int
ged_draw_obol_view_feature_ref_is_null(ged_view_edit_ref ref)
{
    return (!ref.owner || !ref.id || !ref.generation) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_remove_ref(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    return controller ?
	   (controller->features().remove(
		ged_obol_feature_handle_from_ged_ref(ref)) ? 1 : 0) : 0;
}

extern "C" int
ged_draw_obol_view_context_edit_preview_publish_event(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref feature,
    enum ged_view_edit_preview_event event,
    const char *source_path)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, feature);
    if (!controller)
	return 0;

    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(feature);
    BObolFeatureRecord record;
    if (!controller->features().record(handle, record) ||
	record.kind != BObolFeatureKind::EditPreview)
	return 0;

    if (source_path && source_path[0] && record.overlay.isOverlay) {
	BObolOverlayInfo overlay = record.overlay;
	if (overlay.sourcePath != source_path) {
	    overlay.sourcePath = source_path;
	}
	(void)controller->features().setOverlayInfo(handle, overlay);
    }

    switch (event) {
	case GED_VIEW_EDIT_PREVIEW_COMMIT:
	case GED_VIEW_EDIT_PREVIEW_CANCEL:
	case GED_VIEW_EDIT_PREVIEW_DISCARD:
	    return controller->features().remove(handle) ? 1 : 0;
	case GED_VIEW_EDIT_PREVIEW_BEGIN:
	case GED_VIEW_EDIT_PREVIEW_UPDATE:
	case GED_VIEW_EDIT_PREVIEW_REPLACE_SOURCE:
	default:
	    break;
    }

    {
	return controller->features().touch(handle) ? 1 : 0;
    }
}

extern "C" ged_view_edit_ref
ged_draw_obol_view_context_feature_overlay_ensure(
    struct ged_view_context *view_ctx,
    const char *name,
    const void *owner,
    const char *source_path)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller || !name || !name[0])
	return GED_VIEW_EDIT_REF_NULL;

    BObolFeatureOwner feature_owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle = controller->features().findOwned(name,
				  BOBOL_FEATURE_SCOPE_LOCAL, &feature_owner);

    BObolFeatureRecord record;
    const int needs_publish = !handle.isValid() ||
			      !controller->features().record(handle, record) ||
			      record.kind != BObolFeatureKind::EditPreview;
    if (needs_publish) {
	std::vector<SbVec3f> points;
	std::vector<int32_t> commands;
	handle = controller->features().publishEditPreview(name,
		 source_path && source_path[0] ? source_path : name,
		 name,
		 "edit-handle",
		 points,
		 commands,
		 0,
		 0,
		 &feature_owner);
    }

    if (!handle.isValid())
	return GED_VIEW_EDIT_REF_NULL;

    (void)controller->features().setOverlayInfo(handle,
	    ged_obol_edit_overlay_info(view_ctx, owner, source_path, 0));
    return ged_obol_ged_feature_ref(view_ctx, controller, handle, 1);
}

extern "C" ged_view_edit_ref
ged_draw_obol_view_context_feature_label_ensure(
    struct ged_view_context *view_ctx,
    const char *name,
    const void *owner)
{
    BObolViewController *controller =
	ged_obol_view_controller_ensure_for_context(view_ctx, 1);
    if (!controller || !name || !name[0])
	return GED_VIEW_EDIT_REF_NULL;

    BObolFeatureOwner feature_owner = ged_obol_feature_owner(view_ctx, 1);
    BObolFeatureHandle handle = controller->features().findOwned(name,
				  BOBOL_FEATURE_SCOPE_LOCAL, &feature_owner);

    BObolFeatureRecord record;
    if (!handle.isValid() ||
	!controller->features().record(handle, record) ||
	record.kind != BObolFeatureKind::Labels) {
	std::vector<BObolLabel> labels;
	handle = controller->features().publishLabels(name,
		 BObolFeatureScope::Local, labels, NULL, &feature_owner);
    }

    if (!handle.isValid())
	return GED_VIEW_EDIT_REF_NULL;

    (void)controller->features().setOverlayInfo(handle,
	    ged_obol_edit_overlay_info(view_ctx, owner, name, 1));
    return ged_obol_ged_feature_ref(view_ctx, controller, handle, 1);
}

extern "C" int
ged_draw_obol_view_feature_set_visible(struct ged_view_context *view_ctx,
				       ged_view_edit_ref ref,
				       int visible)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    return controller ?
	   (controller->features().setVisible(
		ged_obol_feature_handle_from_ged_ref(ref),
		visible ? TRUE : FALSE) ? 1 : 0) : 0;
}

extern "C" int
ged_draw_obol_view_feature_set_color(struct ged_view_context *view_ctx,
				     ged_view_edit_ref ref,
				     int r,
				     int g,
				     int b)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;

    const unsigned char rgb[3] = {
	ged_obol_rgb_byte_from_int(r),
	ged_obol_rgb_byte_from_int(g),
	ged_obol_rgb_byte_from_int(b)
    };
    return controller->features().setColor(
	       ged_obol_feature_handle_from_ged_ref(ref),
	       ged_obol_color_from_rgb(rgb)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_touch(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;

    return controller->features().touch(
	ged_obol_feature_handle_from_ged_ref(ref)) ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_labels_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
    const struct ged_view_feature_label *labels,
    size_t label_count)
{
    if (label_count && !labels)
	return 0;

    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;

    BObolFeatureRecord record;
    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(ref);
    if (!controller->features().record(handle, record))
	return 0;

    std::vector<BObolLabel> obol_labels =
	ged_obol_labels_from_ged_feature(labels, label_count);
    if (record.kind == BObolFeatureKind::Labels ||
	record.kind == BObolFeatureKind::HudLabel)
	return controller->features().replaceLabels(handle,
		obol_labels) ? 1 : 0;

    BObolFeatureOwner owner = record.owner;
    BObolFeatureHandle labels_handle =
	controller->features().publishLabels(record.name,
	    record.scope, obol_labels, &record.style,
	    record.scope == BObolFeatureScope::Local ? &owner : NULL);
    return labels_handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_points_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
    enum ged_view_edit_geometry_family family,
    const point_t *points,
    const int *cmds,
    size_t point_count)
{
    if (point_count && !points)
	return 0;

    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;

    BObolFeatureRecord record;
    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(ref);
    if (!controller->features().record(handle, record))
	return 0;

    std::vector<SbVec3f> obol_points =
	ged_obol_points_from_ged(points, point_count);
    std::vector<int32_t> obol_commands =
	ged_obol_commands_from_ged(cmds, point_count);
    if (family == GED_VIEW_EDIT_GEOMETRY_TRANSIENT_PREVIEW ||
	record.kind == BObolFeatureKind::EditPreview) {
	if (record.kind == BObolFeatureKind::EditPreview)
	    return controller->features().replaceEditPreviewGeometry(
		       handle,
		       record.name,
		       obol_points,
		       obol_commands) ? 1 : 0;

	BObolFeatureOwner owner = record.owner;
	BObolFeatureHandle preview_handle =
	    controller->features().publishEditPreview(record.name,
		record.name,
		record.name,
		"edit-handle",
		obol_points,
		obol_commands,
		0,
		0,
		record.scope == BObolFeatureScope::Local ? &owner : NULL);
	return preview_handle.isValid() ? 1 : 0;
    }

    BObolFeatureOwner owner = record.owner;
    BObolFeatureHandle line_handle =
	controller->features().publishLineSet(record.name,
	    record.scope, obol_points, obol_commands, &record.style,
	    record.scope == BObolFeatureScope::Local ? &owner : NULL);
    return line_handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_edit_preview_replace(
    struct ged_view_context *view_ctx,
    ged_view_edit_ref ref,
    const char *source_path,
    const char *edit_intent_id,
    const char *edit_intent_role,
    const point_t *points,
    const int *cmds,
    size_t point_count,
    uint32_t source_revision,
    uint32_t inputs_revision)
{
    if (point_count && !points)
	return 0;

    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    if (!controller)
	return 0;

    BObolFeatureRecord record;
    BObolFeatureHandle handle = ged_obol_feature_handle_from_ged_ref(ref);
    if (!controller->features().record(handle, record))
	return 0;

    std::vector<SbVec3f> obol_points =
	ged_obol_points_from_ged(points, point_count);
    std::vector<int32_t> obol_commands =
	ged_obol_commands_from_ged(cmds, point_count);

    BObolFeatureOwner owner = record.owner;
    const SbString preview_name = record.name;
    const SbString preview_path =
	(source_path && source_path[0]) ? SbString(source_path) :
	record.identity.getLength() > 0 ? record.identity : record.name;
    const SbString intent_id =
	(edit_intent_id && edit_intent_id[0]) ? SbString(edit_intent_id) :
	record.editIntentId.getLength() > 0 ? record.editIntentId : record.name;
    const SbString intent_role =
	(edit_intent_role && edit_intent_role[0]) ? SbString(edit_intent_role) :
	record.editIntentRole.getLength() > 0 ? record.editIntentRole :
	SbString("preview");
    const uint32_t next_source_revision =
	source_revision ? source_revision :
	record.sourceRevision ? record.sourceRevision + 1 : 0;
    const uint32_t next_inputs_revision =
	inputs_revision ? inputs_revision :
	record.inputsRevision ? record.inputsRevision + 1 : 0;

    BObolFeatureHandle preview_handle =
	controller->features().publishEditPreview(preview_name,
	    preview_path,
	    intent_id,
	    intent_role,
	    obol_points,
	    obol_commands,
	    next_source_revision,
	    next_inputs_revision,
	    record.scope == BObolFeatureScope::Local ? &owner : NULL);
    return preview_handle.isValid() ? 1 : 0;
}

extern "C" int
ged_draw_obol_view_feature_clear_geometry(struct ged_view_context *view_ctx,
	ged_view_edit_ref ref)
{
    BObolViewController *controller =
	ged_obol_feature_controller_from_ged_ref(view_ctx, ref);
    return controller ?
	   (controller->features().clearGeometry(
		ged_obol_feature_handle_from_ged_ref(ref)) ? 1 : 0) : 0;
}

extern "C" size_t
ged_draw_obol_view_context_label_count(struct ged_view_context *view_ctx, const char *name)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().labels(lookup.handle, labels))
	return 0;
    return labels.size();
}

extern "C" int
ged_draw_obol_view_context_label_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    struct bu_vls *text,
    point_t point,
    unsigned char rgb[3])
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<BObolLabel> labels;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().labels(lookup.handle, labels) ||
	index >= labels.size())
	return 0;

    const BObolLabel &label = labels[index];
    if (text) {
	bu_vls_trunc(text, 0);
	bu_vls_strcat(text, label.text.getString());
    }
    if (point)
	ged_obol_point_from_sb(point, label.point);
    if (rgb) {
	BObolFeatureStyle style;
	const SbColor color = label.hasColor ? label.color :
			      (lookup.controller->features().style(lookup.handle, style) &&
			       style.hasColor ?
			       style.color : SbColor(1.0f, 1.0f, 1.0f));
	ged_obol_rgb_from_color(color, rgb);
    }
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    point_t **points,
    size_t *point_count)
{
    if (points)
	*points = NULL;
    if (point_count)
	*point_count = 0;
    if (!points || !point_count)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<SbVec3f> obol_points;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().points(lookup.handle,
		obol_points))
	return 0;

    if (obol_points.empty())
	return 1;
    *points = (point_t *)bu_calloc(obol_points.size(), sizeof(point_t),
				   "GED Obol feature points copy");
    for (size_t i = 0; i < obol_points.size(); i++)
	ged_obol_point_from_sb((*points)[i], obol_points[i]);
    *point_count = obol_points.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_indices_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    int **indices,
    size_t *index_count)
{
    if (indices)
	*indices = NULL;
    if (index_count)
	*index_count = 0;
    if (!indices || !index_count)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    std::vector<int32_t> obol_indices;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().indices(lookup.handle,
		obol_indices))
	return 0;

    if (obol_indices.empty())
	return 1;
    *indices = (int *)bu_calloc(obol_indices.size(), sizeof(int),
				"GED Obol feature indices copy");
    for (size_t i = 0; i < obol_indices.size(); i++)
	(*indices)[i] = static_cast<int>(obol_indices[i]);
    *index_count = obol_indices.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_line_command_at(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t index,
    int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    int32_t command = 0;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().lineCommandAt(lookup.handle,
		index, command))
	return 0;
    *out = ged_obol_line_command_to_ged(command);
    return 1;
}

static int
ged_obol_feature_layer_lookup(struct ged_view_context *view_ctx,
			      const char *name,
			      size_t layer_index,
			      BObolLineLayer &layer)
{
    ged_obol_feature_lookup lookup =
	ged_obol_feature_lookup_for_context(view_ctx, name);
    BObolFeatureRecord record;
    if (!lookup.handle.isValid() || !lookup.controller ||
	!lookup.controller->features().record(lookup.handle,
		record) || record.kind != BObolFeatureKind::LineLayer ||
	layer_index >= record.layers.size())
	return 0;

    layer = record.layers[layer_index];
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_layer_points_copy(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t layer_index,
    point_t **points,
    size_t *point_count)
{
    if (points)
	*points = NULL;
    if (point_count)
	*point_count = 0;
    if (!points || !point_count)
	return 0;

    BObolLineLayer layer;
    if (!ged_obol_feature_layer_lookup(view_ctx, name, layer_index, layer))
	return 0;

    if (layer.points.empty())
	return 1;
    *points = (point_t *)bu_calloc(layer.points.size(), sizeof(point_t),
				   "GED Obol feature layer points copy");
    for (size_t i = 0; i < layer.points.size(); i++)
	ged_obol_point_from_sb((*points)[i], layer.points[i]);
    *point_count = layer.points.size();
    return 1;
}

extern "C" int
ged_draw_obol_view_context_feature_layer_line_command_at(
    struct ged_view_context *view_ctx,
    const char *name,
    size_t layer_index,
    size_t point_index,
    int *out)
{
    if (out)
	*out = 0;
    if (!out)
	return 0;

    BObolLineLayer layer;
    if (!ged_obol_feature_layer_lookup(view_ctx, name, layer_index, layer) ||
	point_index >= layer.commands.size())
	return 0;

    *out = ged_obol_line_command_to_ged(layer.commands[point_index]);
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
