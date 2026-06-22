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

#include "../libged/bsg_ged_draw_view_private.h"

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
_tclcad_draw_view_points_copy(struct bsg_view *view,
			      const char *feature_name,
			      point_t **pts_out)
{
    if (!pts_out)
	return 0;
    *pts_out = NULL;
    size_t total = 0;
    if (!ged_draw_view_feature_points_copy(view, feature_name, pts_out, &total) ||
	    !*pts_out)
	return 0;
    int count = _tclcad_draw_view_point_count(total);
    if (!count) {
	bu_free(*pts_out, "TclCAD draw-view points");
	*pts_out = NULL;
    }
    return count;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC int
_tclcad_draw_view_axes_centers_copy(struct bsg_view *view,
				    const char *feature_name,
				    point_t **pts_out)
{
    if (!pts_out)
	return 0;
    *pts_out = NULL;
    size_t total = 0;
    if (!ged_draw_view_feature_axes_centers_copy(view, feature_name,
		pts_out, &total) || !*pts_out)
	return 0;
    int count = _tclcad_draw_view_point_count(total);
    if (!count) {
	bu_free(*pts_out, "TclCAD draw-view axes points");
	*pts_out = NULL;
    }
    return count;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_style_read(struct bsg_view *view,
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

    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    if (!ged_draw_view_feature_style_get(view, feature_name, &style))
	return;

    if (color_out && style.color_valid) {
	color_out[0] = (int)style.color[0];
	color_out[1] = (int)style.color[1];
	color_out[2] = (int)style.color[2];
    }
    if (lw_out)
	*lw_out = style.line_width;
    if (tip_len_out)
	*tip_len_out = (int)style.arrow_tip_length;
    if (tip_wid_out)
	*tip_wid_out = (int)style.arrow_tip_width;
    if (visible_out)
	*visible_out = style.visible;
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_arrows_replace(struct bsg_view *view,
				 const char *feature_name,
				 point_t *pts,
				 int npts,
				 int color[3],
				 int line_width,
				 int tip_length,
				 int tip_width,
				 int visible)
{
    if (!view || !feature_name || !pts || npts < 2)
	return;

    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    style.visible = visible;
    style.color_valid = color ? 1 : 0;
    if (color) {
	style.color[0] = (unsigned char)color[0];
	style.color[1] = (unsigned char)color[1];
	style.color[2] = (unsigned char)color[2];
    }
    style.line_width = line_width;
    style.arrow = 1;
    style.arrow_tip_length = (fastf_t)tip_length;
    style.arrow_tip_width = (fastf_t)tip_width;

    (void)ged_draw_view_tcl_arrows_replace(view, feature_name,
	    (const point_t *)pts, (size_t)npts, &style);
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_lines_replace(struct bsg_view *view,
				const char *feature_name,
				point_t *pts,
				int npts,
				int color[3],
				int line_width,
				int visible)
{
    if (!view || !feature_name || !pts || npts < 2)
	return;

    struct ged_draw_view_line_style style = {{255, 255, 0}, 0};
    if (color) {
	style.color[0] = color[0];
	style.color[1] = color[1];
	style.color[2] = color[2];
    }
    style.line_width = line_width;

    if (ged_draw_view_tcl_lines_replace(view, feature_name,
		(const point_t *)pts, (size_t)npts, &style))
	(void)ged_draw_view_feature_visible_set(view, feature_name, visible);
}

TCLCAD_DRAW_VIEW_HELPER_STATIC void
_tclcad_draw_view_axes_replace(struct bsg_view *view,
			       const char *feature_name,
			       point_t *centers,
			       int center_count,
			       fastf_t half_axes_size,
			       int color[3],
			       int line_width,
			       int visible)
{
    if (!view || !feature_name || !centers || center_count < 1)
	return;

    struct ged_draw_view_feature_style style = GED_DRAW_VIEW_FEATURE_STYLE_INIT;
    style.visible = visible;
    style.color_valid = color ? 1 : 0;
    if (color) {
	style.color[0] = (unsigned char)color[0];
	style.color[1] = (unsigned char)color[1];
	style.color[2] = (unsigned char)color[2];
    }
    style.line_width = line_width;

    (void)ged_draw_view_tcl_axes_replace(view, feature_name,
	    (const point_t *)centers, (size_t)center_count,
	    half_axes_size, &style);
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
