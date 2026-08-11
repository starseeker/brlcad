/*             D R A W _ V I E W _ M O V E _ H E L P E R S . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 *
 * This library is distributed in the hope that it will be useful, but
 * WITHOUT ANY WARRANTY; without even the implied warranty of
 * MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
 * Lesser General Public License for more details.
 *
 * You should have received a copy of the GNU Lesser General Public
 * License along with this file; see the file named COPYING for more
 * information.
 */
/**
 * @file draw_view_move_helpers.h
 *
 * Private helpers for retained draw-view pick/move/scale operations in
 * libtclcad.
 */

#ifndef LIBTCLCAD_DRAW_VIEW_MOVE_HELPERS_H
#define LIBTCLCAD_DRAW_VIEW_MOVE_HELPERS_H

#include "common.h"

#include <limits.h>

#include "bu/malloc.h"
#include "vmath.h"

#include "ged/draw.h"
#include "../libged/ged_view_data_line_private.h"

__BEGIN_DECLS

#if defined(__GNUC__)
#  define TCLCAD_DRAW_VIEW_HELPER_STATIC static __attribute__((unused))
#else
#  define TCLCAD_DRAW_VIEW_HELPER_STATIC static
#endif

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_point_count(size_t point_count)
{
    return (point_count > (size_t)INT_MAX) ? 0 : (int)point_count;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_data_lines_points_copy(struct ged_view_context *view_ctx,
					 const char *feature_name,
					 point_t **pts_out)
{
    if (!pts_out)
	return 0;
    *pts_out = NULL;
    struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
    const int staged = strstr(feature_name, "_sdata_") ? 1 : 0;
    if (!ged_view_data_line_state_get(view_ctx, staged, &state))
	return 0;
    int count = _tclcad_draw_view_point_count(state.point_count);
    if (!count || !state.points) {
	ged_view_data_line_state_clear(&state);
	return 0;
    }
    *pts_out = state.points;
    state.points = NULL;
    ged_view_data_line_state_clear(&state);
    return count;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_data_arrows_points_copy(struct ged_view_context *view_ctx,
					  const char *feature_name,
					  point_t **pts_out)
{
    if (!pts_out)
	return 0;
    *pts_out = NULL;
    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return 0;
    struct tclcad_data_arrow_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_arrows :
	&view_state->gv_data_arrows;
    int count = state->gdas_num_points;
    if (count <= 0 || !state->gdas_points)
	return 0;
    *pts_out = (point_t *)bu_calloc((size_t)count, sizeof(point_t),
	"TclCAD draw-view arrow points");
    for (int i = 0; i < count; i++)
	VMOVE((*pts_out)[i], state->gdas_points[i]);
    return count;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_data_axes_centers_copy(struct ged_view_context *view_ctx,
					 const char *feature_name,
					 point_t **pts_out)
{
    if (!pts_out)
	return 0;
    *pts_out = NULL;
    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return 0;
    struct tclcad_data_axes_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_axes :
	&view_state->gv_data_axes;
    int count = state->num_points;
    if (count <= 0 || !state->points)
	return 0;
    *pts_out = (point_t *)bu_calloc((size_t)count, sizeof(point_t),
	"TclCAD draw-view axes points");
    for (int i = 0; i < count; i++)
	VMOVE((*pts_out)[i], state->points[i]);
    return count;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_data_lines_style_read(struct ged_view_context *view_ctx,
					const char *feature_name,
					int color_out[3],
					int *lw_out,
					int *visible_out)
{
    if (color_out) {
	color_out[0] = 255;
	color_out[1] = 255;
	color_out[2] = 0;
    }
    if (lw_out)
	*lw_out = 0;
    if (visible_out)
	*visible_out = 1;

    struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
    if (!ged_view_data_line_state_get(view_ctx,
	    strstr(feature_name, "_sdata_") ? 1 : 0, &state))
	return;
    if (color_out)
	VMOVE(color_out, state.color);
    if (lw_out)
	*lw_out = state.line_width;
    if (visible_out)
	*visible_out = state.draw;
    ged_view_data_line_state_clear(&state);
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_data_arrows_style_read(struct ged_view_context *view_ctx,
					 const char *feature_name,
					 int color_out[3],
					 int *lw_out,
					 int *tip_len_out,
					 int *tip_wid_out,
					 int *visible_out)
{
    if (color_out) {
	color_out[0] = 255;
	color_out[1] = 255;
	color_out[2] = 0;
    }
    if (lw_out)
	*lw_out = 0;
    if (tip_len_out)
	*tip_len_out = 0;
    if (tip_wid_out)
	*tip_wid_out = 0;
    if (visible_out)
	*visible_out = 1;

    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return;
    struct tclcad_data_arrow_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_arrows :
	&view_state->gv_data_arrows;
    if (color_out) {
	color_out[0] = state->gdas_color[0];
	color_out[1] = state->gdas_color[1];
	color_out[2] = state->gdas_color[2];
    }
    if (lw_out)
	*lw_out = state->gdas_line_width;
    if (tip_len_out)
	*tip_len_out = state->gdas_tip_length;
    if (tip_wid_out)
	*tip_wid_out = state->gdas_tip_width;
    if (visible_out)
	*visible_out = state->gdas_draw;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_data_axes_style_read(struct ged_view_context *view_ctx,
				       const char *feature_name,
				       int color_out[3],
				       int *lw_out,
				       int *visible_out)
{
    if (color_out) {
	color_out[0] = 255;
	color_out[1] = 255;
	color_out[2] = 0;
    }
    if (lw_out)
	*lw_out = 0;
    if (visible_out)
	*visible_out = 1;

    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return;
    struct tclcad_data_axes_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_axes :
	&view_state->gv_data_axes;
    if (color_out) {
	color_out[0] = state->color[0];
	color_out[1] = state->color[1];
	color_out[2] = state->color[2];
    }
    if (lw_out)
	*lw_out = state->line_width;
    if (visible_out)
	*visible_out = state->draw;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC fastf_t
_tclcad_draw_view_data_axes_half_size(struct ged_view_context *view_ctx,
	const char *feature_name)
{
    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return 0.0;
    struct tclcad_data_axes_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_axes :
	&view_state->gv_data_axes;
    const struct bv *view = bv_context_view_const(
	(const struct bv_context *)view_ctx);
    const int width = bv_context_width_get((const struct bv_context *)view_ctx);
    const fastf_t scale = bv_size_get(view) /
	(fastf_t)(width > 0 ? width : 512);
    return state->size * 0.5 * scale;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_data_arrows_replace(struct ged_view_context *view_ctx,
				      const char *feature_name,
				      point_t *pts,
				      int npts,
				      int color[3],
				      int line_width,
				      int tip_length,
				      int tip_width,
				      int visible)
{
    if (!view_ctx || !feature_name || !pts || npts < 2)
	return;

    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return;
    struct tclcad_data_arrow_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_arrows :
	&view_state->gv_data_arrows;
    point_t *copy = (point_t *)bu_calloc((size_t)npts, sizeof(point_t),
	"TclCAD moved arrow points");
    for (int i = 0; i < npts; i++)
	VMOVE(copy[i], pts[i]);
    if (state->gdas_points)
	bu_free(state->gdas_points, "TclCAD arrow points");
    state->gdas_points = copy;
    state->gdas_num_points = npts;
    state->gdas_draw = visible;
    state->gdas_line_width = line_width;
    state->gdas_tip_length = tip_length;
    state->gdas_tip_width = tip_width;
    if (color)
	VSET(state->gdas_color, color[0], color[1], color[2]);
    (void)tclcad_data_arrows_publish(view_ctx, feature_name, state);
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_data_lines_replace(struct ged_view_context *view_ctx,
				     const char *feature_name,
				     point_t *pts,
				     int npts,
				     int color[3],
				     int line_width,
				     int visible)
{
    if (!view_ctx || !feature_name || !pts || npts < 2)
	return;

    const int staged = strstr(feature_name, "_sdata_") ? 1 : 0;
    struct ged_view_data_line_state state = GED_VIEW_DATA_LINE_STATE_INIT;
    state.draw = visible;
    if (color)
	VMOVE(state.color, color);
    state.line_width = line_width;
    state.points = (point_t *)bu_calloc((size_t)npts, sizeof(point_t),
	"TclCAD moved line points");
    state.point_count = (size_t)npts;
    for (int i = 0; i < npts; i++)
	VMOVE(state.points[i], pts[i]);
    if (ged_view_data_line_state_replace(view_ctx, staged, &state))
	(void)ged_view_data_line_state_publish(view_ctx, staged);
    ged_view_data_line_state_clear(&state);
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_data_axes_replace(struct ged_view_context *view_ctx,
				    const char *feature_name,
				    point_t *centers,
				    int center_count,
				    fastf_t half_axes_size,
				    int color[3],
				    int line_width,
				    int visible)
{
    if (!view_ctx || !feature_name || !centers || center_count < 1)
	return;

    tclcad_view_state *view_state =
	tclcad_view_tcl_data_from_view_ctx(view_ctx);
    if (!view_state)
	return;
    struct tclcad_data_axes_state *state =
	strstr(feature_name, "_sdata_") ? &view_state->gv_sdata_axes :
	&view_state->gv_data_axes;
    point_t *copy = (point_t *)bu_calloc((size_t)center_count,
	sizeof(point_t), "TclCAD moved axes points");
    for (int i = 0; i < center_count; i++)
	VMOVE(copy[i], centers[i]);
    if (state->points)
	bu_free(state->points, "TclCAD axes points");
    state->points = copy;
    state->num_points = center_count;
    state->draw = visible;
    state->line_width = line_width;
    if (color)
	VSET(state->color, color[0], color[1], color[2]);
    const struct bv *view = bv_context_view_const(
	(const struct bv_context *)view_ctx);
    const int width = bv_context_width_get((const struct bv_context *)view_ctx);
    const fastf_t scale = bv_size_get(view) /
	(fastf_t)(width > 0 ? width : 512);
    state->size = scale > SMALL_FASTF ? half_axes_size * 2.0 / scale : 0.0;
    (void)tclcad_data_axes_publish(view_ctx, feature_name, state);
}

TCLCAD_DRAW_VIEW_HELPER_STATIC size_t
_tclcad_draw_view_data_labels_count(struct ged_view_context *view_ctx,
				    const char *feature_name)
{
    tclcad_label_state *state = tclcad_view_label_state_from_view_ctx(
	view_ctx, strstr(feature_name, "_sdata_") ? 1 : 0);
    return state && state->gdls_num_labels > 0 ?
	(size_t)state->gdls_num_labels : 0;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_data_label_copy(struct ged_view_context *view_ctx,
				  const char *feature_name,
				  size_t index,
				  struct bu_vls *text,
				  point_t point,
				  unsigned char rgb[3])
{
    tclcad_label_state *state = tclcad_view_label_state_from_view_ctx(
	view_ctx, strstr(feature_name, "_sdata_") ? 1 : 0);
    if (!state || index >= (size_t)state->gdls_num_labels)
	return 0;
    if (text)
	bu_vls_strcpy(text, state->gdls_labels[index]);
    if (point)
	VMOVE(point, state->gdls_points[index]);
    if (rgb) {
	rgb[0] = (unsigned char)state->gdls_color[0];
	rgb[1] = (unsigned char)state->gdls_color[1];
	rgb[2] = (unsigned char)state->gdls_color[2];
    }
    return 1;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_data_label_point_set(struct ged_view_context *view_ctx,
				       const char *feature_name,
				       size_t index,
				       const point_t point)
{
    tclcad_label_state *state = tclcad_view_label_state_from_view_ctx(
	view_ctx, strstr(feature_name, "_sdata_") ? 1 : 0);
    if (!state || !point || index >= (size_t)state->gdls_num_labels)
	return 0;
    VMOVE(state->gdls_points[index], point);
    return tclcad_data_labels_publish(view_ctx, state, feature_name);
}

__END_DECLS

#undef TCLCAD_DRAW_VIEW_HELPER_STATIC

#endif /* LIBTCLCAD_DRAW_VIEW_MOVE_HELPERS_H */

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
