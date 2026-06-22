/*             V I E W _ L E G A C Y _ B S G . C
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
/** @file rt/view_legacy_bsg.c
 *
 * Transitional BSG view adapters for callers that still own bsg_view state.
 */

#include "common.h"

#include <string.h>

#include "bg/plane.h"
#include "bn/qmath.h"
#include "bsg/faceplate.h"
#include "bsg/interaction.h"
#include "bsg/lod.h"
#include "bsg/pick.h"
#include "bsg/polygon.h"
#include "bsg/selection.h"
#include "bsg/snap.h"
#include "bsg/snap_action.h"
#include "bsg/util.h"
#include "bsg/view_set.h"
#include "bsg/view_state.h"
#include "rt/primitives/sketch_legacy_bsg.h"
#include "rt/view_legacy_bsg.h"

struct bsg_mesh_lod *
_rt_mesh_lod_bsg(struct rt_mesh_lod *lod);

void
rt_view_info_from_bsg(struct rt_view_info *info, const struct bsg_view *v)
{
    if (!info)
	return;

    rt_view_info_init(info);
    if (!v)
	return;

    info->width = v->gv_width;
    info->height = v->gv_height;
    info->size = v->gv_size;

    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    if (rt_view_lod_policy_from_bsg(&policy, v)) {
	info->lod.scale = policy.scale;
	info->lod.curve_scale = policy.curve_scale;
	info->lod.point_scale = policy.point_scale;
	info->lod.bot_threshold = policy.bot_threshold;
    }
    rt_view_info_sanitize(info);
}

const char *
rt_view_name_from_bsg(const struct bsg_view *v)
{
    return (v && bu_vls_strlen(&v->gv_name)) ? bu_vls_cstr(&v->gv_name) : NULL;
}

int
rt_view_width_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_width : 0;
}

int
rt_view_height_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_height : 0;
}

fastf_t
rt_view_radius_from_bsg(const struct bsg_view *v)
{
    return v ? v->radius : 1.0;
}

int
rt_view_dimensions_set_bsg(struct bsg_view *v, int width, int height)
{
    if (!v)
	return 0;

    v->gv_width = width;
    v->gv_height = height;
    return 1;
}

int
rt_view_screen_to_view_from_bsg(fastf_t *fx,
				fastf_t *fy,
				struct bsg_view *v,
				fastf_t x,
				fastf_t y)
{
    if (fx)
	*fx = 0.0;
    if (fy)
	*fy = 0.0;
    if (!fx || !fy || !v)
	return 0;

    return bsg_screen_to_view(v, fx, fy, x, y) == 0;
}

int
rt_view_screen_point_from_bsg(point_t point,
			      struct bsg_view *v,
			      fastf_t x,
			      fastf_t y)
{
    if (point)
	VSETALL(point, 0.0);
    if (!point || !v)
	return 0;

    point_t model_point = VINIT_ZERO;
    if (bsg_screen_pt(&model_point, x, y, v))
	return 0;

    VMOVE(point, model_point);
    return 1;
}

int
rt_view_current_point_from_bsg(point_t point, const struct bsg_view *v)
{
    if (!point)
	return 0;

    if (!v) {
	VSETALL(point, 0.0);
	return 0;
    }

    VMOVE(point, v->gv_point);
    return 1;
}

int
rt_view_current_point_set_bsg(struct bsg_view *v, const point_t point)
{
    if (!v || !point)
	return 0;

    VMOVE(v->gv_point, point);
    return 1;
}

int
rt_view_previous_mouse_from_bsg(fastf_t *x, fastf_t *y, const struct bsg_view *v)
{
    if (x)
	*x = 0.0;
    if (y)
	*y = 0.0;
    if (!x || !y || !v)
	return 0;

    *x = v->gv_prevMouseX;
    *y = v->gv_prevMouseY;
    return 1;
}

int
rt_view_previous_mouse_set_bsg(struct bsg_view *v, fastf_t x, fastf_t y)
{
    if (!v)
	return 0;

    v->gv_prevMouseX = x;
    v->gv_prevMouseY = y;
    return 1;
}

int
rt_view_mouse_delta_settings_from_bsg(struct rt_view_mouse_delta_settings *settings,
				      const struct bsg_view *v)
{
    struct rt_view_mouse_delta_settings zero = RT_VIEW_MOUSE_DELTA_SETTINGS_INIT;

    if (!settings)
	return 0;

    *settings = zero;
    if (!v)
	return 0;

    settings->min_delta = v->gv_minMouseDelta;
    settings->max_delta = v->gv_maxMouseDelta;
    settings->rotate_scale = v->gv_rscale;
    settings->scale_scale = v->gv_sscale;
    return 1;
}

int
rt_view_mouse_state_set_bsg(struct bsg_view *v, int x, int y)
{
    if (!v)
	return 0;

    v->gv_prevMouseX = v->gv_mouse_x;
    v->gv_prevMouseY = v->gv_mouse_y;
    v->gv_mouse_x = x;
    v->gv_mouse_y = y;

    point_t current_point = VINIT_ZERO;
    int have_point = rt_view_screen_point_from_bsg(current_point, v,
	    (fastf_t)x, (fastf_t)y);
    VMOVE(v->gv_point, current_point);
    return have_point;
}

int
rt_view_unique_object_name_bsg(struct bu_vls *oname,
			       const char *seed,
			       struct bsg_view *v)
{
    if (!oname || !v)
	return 0;

    bsg_uniq_obj_name(oname, seed, v);
    return 1;
}

static void
_rt_view_interactive_rect_from_bsg_state(struct rt_view_interactive_rect_state *dst,
					 const struct bsg_interactive_rect_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->active = src->active;
    dst->draw = src->draw;
    dst->line_width = src->line_width;
    dst->line_style = src->line_style;
    dst->pos[0] = src->pos[0];
    dst->pos[1] = src->pos[1];
    dst->dim[0] = src->dim[0];
    dst->dim[1] = src->dim[1];
    dst->x = src->x;
    dst->y = src->y;
    dst->width = src->width;
    dst->height = src->height;
    dst->bg[0] = src->bg[0];
    dst->bg[1] = src->bg[1];
    dst->bg[2] = src->bg[2];
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->cdim[0] = src->cdim[0];
    dst->cdim[1] = src->cdim[1];
    dst->aspect = src->aspect;
}

static void
_rt_view_interactive_rect_to_bsg_state(struct bsg_interactive_rect_state *dst,
				       const struct rt_view_interactive_rect_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->active = src->active;
    dst->draw = src->draw;
    dst->line_width = src->line_width;
    dst->line_style = src->line_style;
    dst->pos[0] = src->pos[0];
    dst->pos[1] = src->pos[1];
    dst->dim[0] = src->dim[0];
    dst->dim[1] = src->dim[1];
    dst->x = src->x;
    dst->y = src->y;
    dst->width = src->width;
    dst->height = src->height;
    dst->bg[0] = src->bg[0];
    dst->bg[1] = src->bg[1];
    dst->bg[2] = src->bg[2];
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->cdim[0] = src->cdim[0];
    dst->cdim[1] = src->cdim[1];
    dst->aspect = src->aspect;
}

int
rt_view_interactive_rect_state_from_bsg(struct rt_view_interactive_rect_state *record,
					const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_interactive_rect_state bsg_record = {0};
    if (!bsg_view_interactive_rect_get(v, &bsg_record))
	return 0;

    _rt_view_interactive_rect_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_interactive_rect_state_set_bsg(struct bsg_view *v,
				       const struct rt_view_interactive_rect_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_interactive_rect_state bsg_record = {0};
    _rt_view_interactive_rect_to_bsg_state(&bsg_record, record);
    bsg_view_interactive_rect_set(v, &bsg_record);
    return 1;
}

static void
_rt_view_adc_from_bsg_state(struct rt_view_adc_state *dst,
			    const struct bsg_adc_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->draw = src->draw;
    dst->dv_x = src->dv_x;
    dst->dv_y = src->dv_y;
    dst->dv_a1 = src->dv_a1;
    dst->dv_a2 = src->dv_a2;
    dst->dv_dist = src->dv_dist;
    VMOVE(dst->pos_model, src->pos_model);
    VMOVE(dst->pos_view, src->pos_view);
    VMOVE(dst->pos_grid, src->pos_grid);
    dst->a1 = src->a1;
    dst->a2 = src->a2;
    dst->dst = src->dst;
    dst->anchor_pos = src->anchor_pos;
    dst->anchor_a1 = src->anchor_a1;
    dst->anchor_a2 = src->anchor_a2;
    dst->anchor_dst = src->anchor_dst;
    VMOVE(dst->anchor_pt_a1, src->anchor_pt_a1);
    VMOVE(dst->anchor_pt_a2, src->anchor_pt_a2);
    VMOVE(dst->anchor_pt_dst, src->anchor_pt_dst);
    dst->line_color[0] = src->line_color[0];
    dst->line_color[1] = src->line_color[1];
    dst->line_color[2] = src->line_color[2];
    dst->tick_color[0] = src->tick_color[0];
    dst->tick_color[1] = src->tick_color[1];
    dst->tick_color[2] = src->tick_color[2];
    dst->line_width = src->line_width;
}

static void
_rt_view_adc_to_bsg_state(struct bsg_adc_state *dst,
			  const struct rt_view_adc_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->draw = src->draw;
    dst->dv_x = src->dv_x;
    dst->dv_y = src->dv_y;
    dst->dv_a1 = src->dv_a1;
    dst->dv_a2 = src->dv_a2;
    dst->dv_dist = src->dv_dist;
    VMOVE(dst->pos_model, src->pos_model);
    VMOVE(dst->pos_view, src->pos_view);
    VMOVE(dst->pos_grid, src->pos_grid);
    dst->a1 = src->a1;
    dst->a2 = src->a2;
    dst->dst = src->dst;
    dst->anchor_pos = src->anchor_pos;
    dst->anchor_a1 = src->anchor_a1;
    dst->anchor_a2 = src->anchor_a2;
    dst->anchor_dst = src->anchor_dst;
    VMOVE(dst->anchor_pt_a1, src->anchor_pt_a1);
    VMOVE(dst->anchor_pt_a2, src->anchor_pt_a2);
    VMOVE(dst->anchor_pt_dst, src->anchor_pt_dst);
    dst->line_color[0] = src->line_color[0];
    dst->line_color[1] = src->line_color[1];
    dst->line_color[2] = src->line_color[2];
    dst->tick_color[0] = src->tick_color[0];
    dst->tick_color[1] = src->tick_color[1];
    dst->tick_color[2] = src->tick_color[2];
    dst->line_width = src->line_width;
}

int
rt_view_adc_state_from_bsg(struct rt_view_adc_state *record,
			   const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_adc_state bsg_record = {0};
    if (!bsg_view_adc_get(v, &bsg_record))
	return 0;

    _rt_view_adc_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_adc_state_set_bsg(struct bsg_view *v,
			  const struct rt_view_adc_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_adc_state bsg_record = {0};
    _rt_view_adc_to_bsg_state(&bsg_record, record);
    bsg_view_adc_set(v, &bsg_record);
    return 1;
}

static void
_rt_view_grid_from_bsg_state(struct rt_view_grid_state *dst,
			     const struct bsg_grid_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->rc = src->rc;
    dst->draw = src->draw;
    dst->adaptive = src->adaptive;
    dst->snap = src->snap;
    VMOVE(dst->anchor, src->anchor);
    dst->res_h = src->res_h;
    dst->res_v = src->res_v;
    dst->res_major_h = src->res_major_h;
    dst->res_major_v = src->res_major_v;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
}

static void
_rt_view_grid_to_bsg_state(struct bsg_grid_state *dst,
			   const struct rt_view_grid_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->rc = src->rc;
    dst->draw = src->draw;
    dst->adaptive = src->adaptive;
    dst->snap = src->snap;
    VMOVE(dst->anchor, src->anchor);
    dst->res_h = src->res_h;
    dst->res_v = src->res_v;
    dst->res_major_h = src->res_major_h;
    dst->res_major_v = src->res_major_v;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
}

int
rt_view_grid_state_from_bsg(struct rt_view_grid_state *record,
			    const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_grid_state bsg_record = {0};
    if (!bsg_view_grid_get(v, &bsg_record))
	return 0;

    _rt_view_grid_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_grid_state_set_bsg(struct bsg_view *v,
			   const struct rt_view_grid_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_grid_state bsg_record = {0};
    _rt_view_grid_to_bsg_state(&bsg_record, record);
    bsg_view_grid_set(v, &bsg_record);
    return 1;
}

static void
_rt_view_axes_from_bsg_state(struct rt_view_axes_state *dst,
			     const struct bsg_axes *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->draw = src->draw;
    VMOVE(dst->axes_pos, src->axes_pos);
    dst->axes_size = src->axes_size;
    dst->line_width = src->line_width;
    dst->axes_color[0] = src->axes_color[0];
    dst->axes_color[1] = src->axes_color[1];
    dst->axes_color[2] = src->axes_color[2];
    dst->pos_only = src->pos_only;
    dst->label_flag = src->label_flag;
    dst->label_color[0] = src->label_color[0];
    dst->label_color[1] = src->label_color[1];
    dst->label_color[2] = src->label_color[2];
    dst->triple_color = src->triple_color;
    dst->tick_enabled = src->tick_enabled;
    dst->tick_length = src->tick_length;
    dst->tick_major_length = src->tick_major_length;
    dst->tick_interval = src->tick_interval;
    dst->ticks_per_major = src->ticks_per_major;
    dst->tick_threshold = src->tick_threshold;
    dst->tick_color[0] = src->tick_color[0];
    dst->tick_color[1] = src->tick_color[1];
    dst->tick_color[2] = src->tick_color[2];
    dst->tick_major_color[0] = src->tick_major_color[0];
    dst->tick_major_color[1] = src->tick_major_color[1];
    dst->tick_major_color[2] = src->tick_major_color[2];
}

static void
_rt_view_axes_to_bsg_state(struct bsg_axes *dst,
			   const struct rt_view_axes_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->draw = src->draw;
    VMOVE(dst->axes_pos, src->axes_pos);
    dst->axes_size = src->axes_size;
    dst->line_width = src->line_width;
    dst->axes_color[0] = src->axes_color[0];
    dst->axes_color[1] = src->axes_color[1];
    dst->axes_color[2] = src->axes_color[2];
    dst->pos_only = src->pos_only;
    dst->label_flag = src->label_flag;
    dst->label_color[0] = src->label_color[0];
    dst->label_color[1] = src->label_color[1];
    dst->label_color[2] = src->label_color[2];
    dst->triple_color = src->triple_color;
    dst->tick_enabled = src->tick_enabled;
    dst->tick_length = src->tick_length;
    dst->tick_major_length = src->tick_major_length;
    dst->tick_interval = src->tick_interval;
    dst->ticks_per_major = src->ticks_per_major;
    dst->tick_threshold = src->tick_threshold;
    dst->tick_color[0] = src->tick_color[0];
    dst->tick_color[1] = src->tick_color[1];
    dst->tick_color[2] = src->tick_color[2];
    dst->tick_major_color[0] = src->tick_major_color[0];
    dst->tick_major_color[1] = src->tick_major_color[1];
    dst->tick_major_color[2] = src->tick_major_color[2];
}

int
rt_view_model_axes_state_from_bsg(struct rt_view_axes_state *record,
				  const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_axes bsg_record = {0};
    if (!bsg_view_model_axes_get(v, &bsg_record))
	return 0;

    _rt_view_axes_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_model_axes_state_set_bsg(struct bsg_view *v,
				 const struct rt_view_axes_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_axes bsg_record = {0};
    _rt_view_axes_to_bsg_state(&bsg_record, record);
    bsg_view_model_axes_set(v, &bsg_record);
    return 1;
}

int
rt_view_view_axes_state_from_bsg(struct rt_view_axes_state *record,
				 const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_axes bsg_record = {0};
    if (!bsg_view_view_axes_get(v, &bsg_record))
	return 0;

    _rt_view_axes_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_view_axes_state_set_bsg(struct bsg_view *v,
				const struct rt_view_axes_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_axes bsg_record = {0};
    _rt_view_axes_to_bsg_state(&bsg_record, record);
    bsg_view_view_axes_set(v, &bsg_record);
    return 1;
}

static void
_rt_view_other_from_bsg_state(struct rt_view_other_state *dst,
			      const struct bsg_other_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->gos_draw = src->gos_draw;
    dst->gos_line_color[0] = src->gos_line_color[0];
    dst->gos_line_color[1] = src->gos_line_color[1];
    dst->gos_line_color[2] = src->gos_line_color[2];
    dst->gos_text_color[0] = src->gos_text_color[0];
    dst->gos_text_color[1] = src->gos_text_color[1];
    dst->gos_text_color[2] = src->gos_text_color[2];
    dst->gos_font_size = src->gos_font_size;
}

static void
_rt_view_other_to_bsg_state(struct bsg_other_state *dst,
			    const struct rt_view_other_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->gos_draw = src->gos_draw;
    dst->gos_line_color[0] = src->gos_line_color[0];
    dst->gos_line_color[1] = src->gos_line_color[1];
    dst->gos_line_color[2] = src->gos_line_color[2];
    dst->gos_text_color[0] = src->gos_text_color[0];
    dst->gos_text_color[1] = src->gos_text_color[1];
    dst->gos_text_color[2] = src->gos_text_color[2];
    dst->gos_font_size = src->gos_font_size;
}

int
rt_view_center_dot_state_from_bsg(struct rt_view_other_state *record,
				  const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_other_state bsg_record = {0};
    if (!bsg_view_center_dot_get(v, &bsg_record))
	return 0;

    _rt_view_other_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_center_dot_state_set_bsg(struct bsg_view *v,
				 const struct rt_view_other_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_other_state bsg_record = {0};
    _rt_view_other_to_bsg_state(&bsg_record, record);
    bsg_view_center_dot_set(v, &bsg_record);
    return 1;
}

int
rt_view_scale_overlay_state_from_bsg(struct rt_view_other_state *record,
				     const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_other_state bsg_record = {0};
    if (!bsg_view_scale_state_get(v, &bsg_record))
	return 0;

    _rt_view_other_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_scale_overlay_state_set_bsg(struct bsg_view *v,
				    const struct rt_view_other_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_other_state bsg_record = {0};
    _rt_view_other_to_bsg_state(&bsg_record, record);
    bsg_view_scale_state_set(v, &bsg_record);
    return 1;
}

static void
_rt_view_params_from_bsg_state(struct rt_view_params_state *dst,
			       const struct bsg_params_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->draw = src->draw;
    dst->draw_size = src->draw_size;
    dst->draw_center = src->draw_center;
    dst->draw_az = src->draw_az;
    dst->draw_el = src->draw_el;
    dst->draw_tw = src->draw_tw;
    dst->draw_fps = src->draw_fps;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->font_size = src->font_size;
}

static void
_rt_view_params_to_bsg_state(struct bsg_params_state *dst,
			     const struct rt_view_params_state *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->draw = src->draw;
    dst->draw_size = src->draw_size;
    dst->draw_center = src->draw_center;
    dst->draw_az = src->draw_az;
    dst->draw_el = src->draw_el;
    dst->draw_tw = src->draw_tw;
    dst->draw_fps = src->draw_fps;
    dst->color[0] = src->color[0];
    dst->color[1] = src->color[1];
    dst->color[2] = src->color[2];
    dst->font_size = src->font_size;
}

int
rt_view_params_state_from_bsg(struct rt_view_params_state *record,
			      const struct bsg_view *v)
{
    if (record)
	memset(record, 0, sizeof(*record));
    if (!record || !v)
	return 0;

    struct bsg_params_state bsg_record = {0};
    if (!bsg_view_params_get(v, &bsg_record))
	return 0;

    _rt_view_params_from_bsg_state(record, &bsg_record);
    return 1;
}

int
rt_view_params_state_set_bsg(struct bsg_view *v,
			     const struct rt_view_params_state *record)
{
    if (!v || !record)
	return 0;

    struct bsg_params_state bsg_record = {0};
    _rt_view_params_to_bsg_state(&bsg_record, record);
    bsg_view_params_set(v, &bsg_record);
    return 1;
}

int
rt_view_refresh_request_bsg(struct bsg_view *v, uint32_t flags)
{
    if (!v)
	return 0;

    bsg_view_refresh_request(v, flags);
    return 1;
}

int
rt_view_refresh_dirty_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_dirty(v);
}

uint32_t
rt_view_refresh_consume_bsg(struct bsg_view *v)
{
    return bsg_view_refresh_consume(v);
}

int
rt_view_refresh_complete_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_refresh_complete(v);
    return 1;
}

int
rt_view_refresh_enabled_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_enabled(v);
}

int
rt_view_refresh_enabled_set_bsg(struct bsg_view *v, int enabled)
{
    if (!v)
	return 0;

    bsg_view_refresh_set_enabled(v, enabled);
    return 1;
}

int
rt_view_refresh_suppressed_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_suppressed(v);
}

int
rt_view_refresh_suppress_begin_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_refresh_suppress_begin(v);
    return 1;
}

int
rt_view_refresh_suppress_end_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_refresh_suppress_end(v);
    return 1;
}

int
rt_view_refresh_drawn_count_from_bsg(const struct bsg_view *v)
{
    return bsg_view_refresh_drawn_count(v);
}

int
rt_view_refresh_drawn_count_set_bsg(struct bsg_view *v, int count)
{
    if (!v)
	return 0;

    bsg_view_refresh_set_drawn_count(v, count);
    return 1;
}

int
rt_view_orientation_quat_from_bsg(quat_t orientation, const struct bsg_view *v)
{
    mat_t identity;

    if (!orientation)
	return 0;

    if (!v) {
	MAT_IDN(identity);
	quat_mat2quat(orientation, identity);
	return 0;
    }

    quat_mat2quat(orientation, v->gv_rotation);
    return 1;
}

int
rt_view_aet_from_bsg(vect_t aet, const struct bsg_view *v)
{
    if (!aet)
	return 0;

    if (!v) {
	VSETALL(aet, 0.0);
	return 0;
    }

    VMOVE(aet, v->gv_aet);
    return 1;
}

int
rt_view_aet_set_bsg(struct bsg_view *v, const vect_t aet)
{
    if (!v || !aet)
	return 0;

    bsg_view_set_aet(v, aet);
    return 1;
}

int
rt_view_aet_state_set_bsg(struct bsg_view *v, const vect_t aet)
{
    if (!v || !aet)
	return 0;

    VMOVE(v->gv_aet, aet);
    return 1;
}

int
rt_view_update_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_update(v);
    return 1;
}

int
rt_view_autoview_bsg(struct bsg_view *v, fastf_t scale, int all_view_objs)
{
    if (!v)
	return 0;

    bsg_autoview(v, scale, all_view_objs);
    return 1;
}

int
rt_view_autoview_bounds_bsg(struct bsg_view *v,
			    fastf_t scale,
			    const point_t min,
			    const point_t max)
{
    if (!v || !min || !max)
	return 0;

    bsg_autoview_bounds(v, scale, min, max);
    return 1;
}

int
rt_view_adjust_bsg(struct bsg_view *v,
		   int dx,
		   int dy,
		   point_t keypoint,
		   int mode,
		   unsigned long long flags)
{
    unsigned long long bsg_flags = BSG_IDLE;

    if (!v || !keypoint || flags == RT_VIEW_ADJUST_IDLE)
	return 0;

    if (flags & RT_VIEW_ADJUST_ROT)
	bsg_flags |= BSG_ROT;
    if (flags & RT_VIEW_ADJUST_TRANS)
	bsg_flags |= BSG_TRANS;
    if (flags & RT_VIEW_ADJUST_SCALE)
	bsg_flags |= BSG_SCALE;
    if (flags & RT_VIEW_ADJUST_CENTER)
	bsg_flags |= BSG_CENTER;
    if (flags & RT_VIEW_ADJUST_CON_X)
	bsg_flags |= BSG_CON_X;
    if (flags & RT_VIEW_ADJUST_CON_Y)
	bsg_flags |= BSG_CON_Y;
    if (flags & RT_VIEW_ADJUST_CON_Z)
	bsg_flags |= BSG_CON_Z;
    if (flags & RT_VIEW_ADJUST_CON_GRID)
	bsg_flags |= BSG_CON_GRID;
    if (flags & RT_VIEW_ADJUST_CON_LINES)
	bsg_flags |= BSG_CON_LINES;

    return bsg_adjust(v, dx, dy, keypoint, mode, bsg_flags);
}

unsigned long long
rt_view_hash_bsg(const struct bsg_view *v)
{
    return v ? bsg_hash((struct bsg_view *)v) : 0ULL;
}

int
rt_view_snap_candidates_bsg(struct bsg_view *v,
			    point_t sample,
			    double tol,
			    unsigned long long kinds,
			    struct bsg_snap_result *out)
{
    return bsg_snap_candidates(v, sample, tol, (bsg_snap_kind_mask)kinds, out);
}

int
rt_view_snap_point_2d_bsg(struct bsg_view *v,
			  fastf_t *vx,
			  fastf_t *vy,
			  unsigned long long kinds)
{
    return bsg_snap_point_2d(v, vx, vy, (bsg_snap_kind_mask)kinds);
}

int
rt_view_snap_grid_2d_bsg(struct bsg_view *v, fastf_t *vx, fastf_t *vy)
{
    if (!v || !vx || !vy)
	return 0;

    return bsg_snap_grid_2d(v, vx, vy);
}

static struct bsg_pick_result *
rt_view_pick_result_to_bsg(struct rt_view_pick_result_bsg *result)
{
    return (struct bsg_pick_result *)result;
}

static const struct bsg_pick_result *
rt_view_pick_result_const_to_bsg(const struct rt_view_pick_result_bsg *result)
{
    return (const struct bsg_pick_result *)result;
}

static struct rt_view_pick_result_bsg *
rt_view_pick_result_from_bsg(struct bsg_pick_result *result)
{
    return (struct rt_view_pick_result_bsg *)result;
}

static const char *
rt_view_pick_record_path(const struct bsg_pick_record *record)
{
    if (!record)
	return NULL;

    const char *source_path = bu_vls_cstr(&record->pr_source_path);
    if (source_path && source_path[0])
	return source_path;

    const char *instance_path = bu_vls_cstr(&record->pr_instance_path);
    return (instance_path && instance_path[0]) ? instance_path : NULL;
}

struct rt_view_pick_result_bsg *
rt_view_pick_point_bsg(struct bsg_view *v, int x, int y, int first_only)
{
    return rt_view_pick_result_from_bsg(bsg_pick_point(v, x, y, first_only));
}

struct rt_view_pick_result_bsg *
rt_view_pick_nearest_bsg(struct bsg_view *v, int x, int y)
{
    return rt_view_pick_result_from_bsg(bsg_pick_nearest(v, x, y));
}

struct rt_view_pick_result_bsg *
rt_view_pick_rect_bsg(struct bsg_view *v, int x0, int y0, int x1, int y1)
{
    return rt_view_pick_result_from_bsg(bsg_pick_rect(v, x0, y0, x1, y1));
}

struct rt_view_pick_result_bsg *
rt_view_pick_semantic_path_bsg(struct bsg_view *v, const char *path_pattern)
{
    return rt_view_pick_result_from_bsg(bsg_pick_semantic_path(v, path_pattern));
}

struct rt_view_pick_result_bsg *
rt_view_pick_result_create_bsg(void)
{
    return rt_view_pick_result_from_bsg(bsg_pick_result_create());
}

void
rt_view_pick_result_free_bsg(struct rt_view_pick_result_bsg *result)
{
    bsg_pick_result_free(rt_view_pick_result_to_bsg(result));
}

size_t
rt_view_pick_result_count_bsg(const struct rt_view_pick_result_bsg *result)
{
    return bsg_pick_result_count(rt_view_pick_result_const_to_bsg(result));
}

int
rt_view_pick_result_path_bsg(const struct rt_view_pick_result_bsg *result,
			     size_t index,
			     struct bu_vls *path_out)
{
    if (!path_out)
	return 0;
    bu_vls_trunc(path_out, 0);

    const struct bsg_pick_result *bsg_result =
	rt_view_pick_result_const_to_bsg(result);
    const struct bsg_pick_record *record =
	bsg_pick_result_get(bsg_result, index);
    const char *path = rt_view_pick_record_path(record);
    if (!path || !path[0])
	return 0;

    bu_vls_sprintf(path_out, "%s", path);
    return 1;
}

fastf_t
rt_view_pick_result_hit_dist_bsg(const struct rt_view_pick_result_bsg *result,
				 size_t index)
{
    const struct bsg_pick_result *bsg_result =
	rt_view_pick_result_const_to_bsg(result);
    const struct bsg_pick_record *record =
	bsg_pick_result_get(bsg_result, index);
    return record ? record->pr_hit_dist : -1.0;
}

int
rt_view_pick_result_append_path_bsg(struct rt_view_pick_result_bsg *result,
				    struct bsg_view *v,
				    int screen_x,
				    int screen_y,
				    const char *source_path,
				    fastf_t hit_dist)
{
    if (!result || !source_path || !source_path[0])
	return 0;

    struct bsg_pick_result *bsg_result = rt_view_pick_result_to_bsg(result);
    struct bsg_pick_record *record;
    BU_GET(record, struct bsg_pick_record);
    bu_vls_init(&record->pr_source_path);
    bu_vls_init(&record->pr_instance_path);
    record->pr_scene = (bsg_scene_ref)BSG_SCENE_REF_NULL_INIT;
    record->pr_feature = (bsg_feature_ref)BSG_FEATURE_REF_NULL_INIT;
    record->pr_valid = 1;
    record->pr_view = v;
    record->pr_screen_x = screen_x;
    record->pr_screen_y = screen_y;
    record->pr_primitive_id = -1;
    record->pr_subelement_id = -1;
    record->pr_hit_dist = hit_dist;
    bu_vls_sprintf(&record->pr_source_path, "%s", source_path);
    bu_vls_sprintf(&record->pr_instance_path, "%s", source_path);
    bu_ptbl_ins(&bsg_result->pr_records, (long *)record);
    return 1;
}

int
rt_view_pick_result_append_copy_bsg(struct rt_view_pick_result_bsg *dest,
				    const struct rt_view_pick_result_bsg *src,
				    size_t index,
				    fastf_t hit_dist)
{
    if (!dest || !src)
	return 0;

    const struct bsg_pick_result *bsg_src =
	rt_view_pick_result_const_to_bsg(src);
    const struct bsg_pick_record *src_record =
	bsg_pick_result_get(bsg_src, index);
    if (!src_record)
	return 0;

    struct bsg_pick_result *bsg_dest = rt_view_pick_result_to_bsg(dest);
    struct bsg_pick_record *record;
    BU_GET(record, struct bsg_pick_record);
    bu_vls_init(&record->pr_source_path);
    bu_vls_init(&record->pr_instance_path);
    record->pr_scene = src_record->pr_scene;
    record->pr_feature = src_record->pr_feature;
    record->pr_valid = src_record->pr_valid;
    record->pr_view = src_record->pr_view;
    record->pr_screen_x = src_record->pr_screen_x;
    record->pr_screen_y = src_record->pr_screen_y;
    record->pr_primitive_id = src_record->pr_primitive_id;
    record->pr_subelement_id = src_record->pr_subelement_id;
    record->pr_hit_dist = hit_dist;
    bu_vls_sprintf(&record->pr_source_path, "%s",
	    bu_vls_cstr(&src_record->pr_source_path));
    bu_vls_sprintf(&record->pr_instance_path, "%s",
	    bu_vls_cstr(&src_record->pr_instance_path));
    bu_ptbl_ins(&bsg_dest->pr_records, (long *)record);
    return 1;
}

struct rt_view_pick_result_bsg *
rt_view_pick_result_filter_first_bsg(const struct rt_view_pick_result_bsg *src)
{
    struct rt_view_pick_result_bsg *result = rt_view_pick_result_create_bsg();
    if (!result)
	return NULL;
    if (rt_view_pick_result_count_bsg(src) > 0)
	rt_view_pick_result_append_copy_bsg(result, src, 0,
		rt_view_pick_result_hit_dist_bsg(src, 0));
    return result;
}

static struct bsg_selection *
rt_view_selection_get_bsg(struct bsg_view *v)
{
    return bsg_view_selection(v);
}

int
rt_view_selection_available_bsg(struct bsg_view *v)
{
    return rt_view_selection_get_bsg(v) ? 1 : 0;
}

size_t
rt_view_selection_count_bsg(struct bsg_view *v)
{
    struct bsg_selection *selection = rt_view_selection_get_bsg(v);
    return selection ? bsg_selection_count(selection) : 0;
}

static int
rt_view_selection_apply_pick_result_bsg(struct bsg_view *v,
					const struct bsg_pick_result *result,
					rt_view_selection_path_callback_bsg_t callback,
					void *data)
{
    struct bsg_interaction_result *interactions = result ?
	bsg_interaction_from_pick_result(result) : NULL;
    if (!interactions) {
	struct bsg_selection *selection = rt_view_selection_get_bsg(v);
	if (!selection)
	    return 0;
	bsg_selection_clear(selection);
	return 1;
    }

    if (callback) {
	for (size_t i = 0; i < bsg_interaction_result_count(interactions); i++) {
	    const struct bsg_interaction_record *record =
		bsg_interaction_result_get(interactions, i);
	    const char *path = bsg_interaction_record_path(record);
	    if (path && path[0])
		callback(path, data);
	}
    }

    struct bsg_selection *selection = rt_view_selection_get_bsg(v);
    if (selection)
	bsg_interaction_selection_apply(selection, interactions,
		BSG_INTERACTION_APPLY_SET);
    bsg_interaction_result_free(interactions);
    return selection ? 1 : 0;
}

int
rt_view_selection_set_pick_result_ref_bsg(struct bsg_view *v,
					  const struct rt_view_pick_result_bsg *result,
					  rt_view_selection_path_callback_bsg_t callback,
					  void *data)
{
    return rt_view_selection_apply_pick_result_bsg(v,
	    rt_view_pick_result_const_to_bsg(result), callback, data);
}

int
rt_view_selection_clear_bsg(struct bsg_view *v)
{
    struct bsg_selection *selection = rt_view_selection_get_bsg(v);
    if (!selection)
	return 0;

    bsg_selection_clear(selection);
    return 1;
}

struct bu_ptbl *
rt_view_set_views_bsg(struct bsg_view_set *s)
{
    return bsg_set_views(s);
}

struct bsg_view *
rt_view_set_find_view_bsg(struct bsg_view_set *s, const char *name)
{
    if (!s || !name)
	return NULL;

    return bsg_set_find_view(s, name);
}

void
rt_view_set_init_bsg(struct bsg_view_set *s)
{
    if (!s)
	return;

    bsg_set_init(s);
}

void
rt_view_set_free_bsg(struct bsg_view_set *s)
{
    if (!s)
	return;

    bsg_set_free(s);
}

void
rt_view_init_bsg(struct bsg_view *v, struct bsg_view_set *s)
{
    bsg_init(v, s);
}

void
rt_view_free_bsg(struct bsg_view *v)
{
    bsg_free(v);
}

void
rt_view_set_add_view_bsg(struct bsg_view_set *s, struct bsg_view *v)
{
    bsg_set_add_view(s, v);
}

void
rt_view_set_remove_view_bsg(struct bsg_view_set *s, struct bsg_view *v)
{
    bsg_set_rm_view(s, v);
}

int
rt_view_knobs_reset_bsg(struct bsg_view *v, int category)
{
    if (!v)
	return 0;

    bsg_knobs_reset(&v->k, category);
    return 1;
}

int
rt_view_knob_state_reset_bsg(struct bsg_view_knobs *knobs, int category)
{
    if (!knobs)
	return 0;

    bsg_knobs_reset(knobs, category);
    return 1;
}

unsigned long long
rt_view_knobs_hash_bsg(struct bsg_view *v,
		       struct bu_data_hash_state *state)
{
    return v ? bsg_knobs_hash(&v->k, state) : 0ULL;
}

int
rt_view_knobs_cmd_process_bsg(vect_t *rvec,
			      int *do_rot,
			      vect_t *tvec,
			      int *do_tran,
			      struct bsg_view *v,
			      const char *cmd,
			      fastf_t factor,
			      char origin,
			      int model_flag,
			      int incr_flag)
{
    return bsg_knobs_cmd_process(rvec, do_rot, tvec, do_tran, v, cmd, factor,
	    origin, model_flag, incr_flag);
}

int
rt_view_knobs_translate_bsg(struct bsg_view *v,
			    const vect_t tvec,
			    int model_flag)
{
    if (!v || !tvec)
	return 0;

    bsg_knobs_tran(v, tvec, model_flag);
    return 1;
}

int
rt_view_knobs_rotate_bsg(struct bsg_view *v,
			 const vect_t rvec,
			 char origin,
			 char coords,
			 const matp_t obj_rot,
			 const pointp_t pvt_pt)
{
    if (!v || !rvec)
	return 0;

    bsg_knobs_rot(v, rvec, origin, coords, obj_rot, pvt_pt);
    return 1;
}

int
rt_view_knobs_update_rate_flags_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_update_rate_flags(v);
    return 1;
}

int
rt_view_is_independent_bsg(const struct bsg_view *v)
{
    return bsg_view_is_independent(v);
}

bsg_scene_ref
rt_view_independent_scope_ref_bsg(struct bsg_view *v, int create)
{
    return bsg_view_independent_scope_ref(v, create);
}

void
rt_view_independent_scope_destroy_bsg(struct bsg_view *v)
{
    bsg_view_independent_scope_destroy(v);
}

void
rt_view_init_copy_bsg(struct bsg_view *dest,
		      const struct bsg_view *src,
		      struct bsg_view_set *s)
{
    bsg_view_init_copy(dest, src, s);
}

void
rt_view_tclcad_data_init_bsg(struct bsg_data_tclcad *d)
{
    if (!d)
	return;

    bsg_data_tclcad_init(d);
}

size_t
rt_view_clear_bsg(struct bsg_view *v, int flags)
{
    return bsg_clear(v, flags);
}

fastf_t
rt_view_perspective_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_perspective : 0.0;
}

int
rt_view_perspective_set_bsg(struct bsg_view *v, fastf_t perspective)
{
    if (!v)
	return 0;

    bsg_view_set_perspective(v, perspective);
    return 1;
}

int
rt_view_model2view_from_bsg(mat_t model2view, const struct bsg_view *v)
{
    if (!model2view)
	return 0;

    if (!v) {
	MAT_IDN(model2view);
	return 0;
    }

    MAT_COPY(model2view, v->gv_model2view);
    return 1;
}

int
rt_view_model2view_set_bsg(struct bsg_view *v, const mat_t model2view)
{
    if (!v || !model2view)
	return 0;

    MAT_COPY(v->gv_model2view, model2view);
    return 1;
}

int
rt_view_view2model_from_bsg(mat_t view2model, const struct bsg_view *v)
{
    if (!view2model)
	return 0;

    if (!v) {
	MAT_IDN(view2model);
	return 0;
    }

    MAT_COPY(view2model, v->gv_view2model);
    return 1;
}

int
rt_view_view2model_set_bsg(struct bsg_view *v, const mat_t view2model)
{
    if (!v || !view2model)
	return 0;

    MAT_COPY(v->gv_view2model, view2model);
    return 1;
}

int
rt_view_pmodel2view_from_bsg(mat_t pmodel2view, const struct bsg_view *v)
{
    if (!pmodel2view)
	return 0;

    if (!v) {
	MAT_IDN(pmodel2view);
	return 0;
    }

    MAT_COPY(pmodel2view, v->gv_pmodel2view);
    return 1;
}

int
rt_view_pmodel2view_set_bsg(struct bsg_view *v, const mat_t pmodel2view)
{
    if (!v || !pmodel2view)
	return 0;

    MAT_COPY(v->gv_pmodel2view, pmodel2view);
    return 1;
}

int
rt_view_pmat_from_bsg(mat_t pmat, const struct bsg_view *v)
{
    if (!pmat)
	return 0;

    if (!v) {
	MAT_IDN(pmat);
	return 0;
    }

    MAT_COPY(pmat, v->gv_pmat);
    return 1;
}

int
rt_view_pmat_set_bsg(struct bsg_view *v, const mat_t pmat)
{
    if (!v || !pmat)
	return 0;

    MAT_COPY(v->gv_pmat, pmat);
    return 1;
}

int
rt_view_rotation_from_bsg(mat_t rotation, const struct bsg_view *v)
{
    if (!rotation)
	return 0;

    if (!v) {
	MAT_IDN(rotation);
	return 0;
    }

    MAT_COPY(rotation, v->gv_rotation);
    return 1;
}

int
rt_view_rotation_set_bsg(struct bsg_view *v, const mat_t rotation)
{
    if (!v || !rotation)
	return 0;

    bsg_view_set_rotation(v, rotation);
    return 1;
}

int
rt_view_center_from_bsg(mat_t center, const struct bsg_view *v)
{
    if (!center)
	return 0;

    if (!v) {
	MAT_IDN(center);
	return 0;
    }

    MAT_COPY(center, v->gv_center);
    return 1;
}

int
rt_view_center_vec_set_bsg(struct bsg_view *v, const point_t center)
{
    if (!v || !center)
	return 0;

    bsg_view_set_center_vec(v, center);
    return 1;
}

int
rt_view_plane_from_bsg(plane_t *plane, const struct bsg_view *v)
{
    point_t center;
    vect_t normal;

    if (!plane || !v)
	return -1;

    MAT_DELTAS_GET_NEG(center, v->gv_center);
    VMOVEN(normal, v->gv_rotation + 8, 3);
    VUNITIZE(normal);
    VSCALE(normal, normal, -1.0);

    return bg_plane_pt_nrml(plane, center, normal);
}

int
rt_view_lod_policy_from_bsg(struct rt_view_lod_policy *policy,
			    const struct bsg_view *v)
{
    if (!policy)
	return 0;

    rt_view_lod_policy_init(policy);
    if (!v)
	return 0;

    struct bsg_lod_source_policy_settings bsg_policy;
    memset(&bsg_policy, 0, sizeof(bsg_policy));
    bsg_policy.policy = BSG_LOD_AUTO;
    bsg_policy.scale = 1.0;
    bsg_policy.curve_scale = 1.0;
    bsg_policy.point_scale = 1.0;
    if (!bsg_view_lod_source_policy_get(v, &bsg_policy))
	return 0;

    policy->policy = (int)bsg_policy.policy;
    policy->forced_level = bsg_policy.forced_level;
    policy->mesh_enabled = bsg_policy.mesh_enabled ? 1 : 0;
    policy->csg_enabled = bsg_policy.csg_enabled ? 1 : 0;
    policy->zoom_refresh = bsg_policy.zoom_refresh ? 1 : 0;
    policy->bot_threshold = bsg_policy.bot_threshold;
    policy->scale = bsg_policy.scale;
    policy->curve_scale = bsg_policy.curve_scale;
    policy->point_scale = bsg_policy.point_scale;
    rt_view_lod_policy_sanitize(policy);
    return 1;
}

int
rt_view_lod_policy_apply_bsg(struct bsg_view *v,
			     const struct rt_view_lod_policy *policy)
{
    if (!v || !policy)
	return 0;

    struct rt_view_lod_policy sanitized = *policy;
    rt_view_lod_policy_sanitize(&sanitized);

    struct bsg_lod_source_policy_settings bsg_policy;
    memset(&bsg_policy, 0, sizeof(bsg_policy));
    bsg_policy.policy = (bsg_lod_policy)sanitized.policy;
    bsg_policy.forced_level = sanitized.forced_level;
    bsg_policy.mesh_enabled = sanitized.mesh_enabled ? 1 : 0;
    bsg_policy.csg_enabled = sanitized.csg_enabled ? 1 : 0;
    bsg_policy.zoom_refresh = sanitized.zoom_refresh ? 1 : 0;
    bsg_policy.scale = sanitized.scale;
    bsg_policy.bot_threshold = sanitized.bot_threshold;
    bsg_policy.curve_scale = sanitized.curve_scale;
    bsg_policy.point_scale = sanitized.point_scale;
    bsg_view_lod_source_policy_set(v, &bsg_policy);
    return 1;
}

int
rt_view_lod_policy_copy_bsg(struct bsg_view *dst, const struct bsg_view *src)
{
    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    if (!dst || !rt_view_lod_policy_from_bsg(&policy, src))
	return 0;
    return rt_view_lod_policy_apply_bsg(dst, &policy);
}

void
rt_view_lod_bounds_update_bsg(struct bsg_view *v)
{
    if (v)
	bsg_view_bounds(v);
}

void
rt_view_lod_bounds_callback_set_bsg(struct bsg_view *v)
{
    if (v)
	v->gv_bounds_update = &bsg_view_bounds;
}

int
rt_view_lod_bounds_callback_is_bsg(const struct bsg_view *v)
{
    return (v && v->gv_bounds_update == &bsg_view_bounds) ? 1 : 0;
}

rt_view_bounds_update_callback_bsg_t
rt_view_bounds_update_callback_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_bounds_update : NULL;
}

int
rt_view_bounds_update_callback_set_bsg(struct bsg_view *v,
				       rt_view_bounds_update_callback_bsg_t callback)
{
    if (!v)
	return 0;

    v->gv_bounds_update = callback;
    return 1;
}

int
rt_view_bounds_update_callback_call_bsg(struct bsg_view *v)
{
    rt_view_bounds_update_callback_bsg_t callback =
	rt_view_bounds_update_callback_from_bsg(v);

    if (!callback)
	return 0;

    (*callback)(v);
    return 1;
}

fastf_t
rt_view_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_scale : 1.0;
}

int
rt_view_scale_set_bsg(struct bsg_view *v, fastf_t scale)
{
    if (!v)
	return 0;

    bsg_view_set_scale(v, scale);
    return 1;
}

fastf_t *
rt_view_scale_storage_from_bsg(struct bsg_view *v)
{
    return v ? &v->gv_scale : NULL;
}

int
rt_view_scale_state_set_bsg(struct bsg_view *v,
			    fastf_t scale,
			    fastf_t initial_scale,
			    fastf_t absolute_scale,
			    fastf_t size,
			    fastf_t inverse_size)
{
    if (!v)
	return 0;

    v->gv_scale = scale;
    v->gv_i_scale = initial_scale;
    v->gv_a_scale = absolute_scale;
    v->gv_size = size;
    v->gv_isize = inverse_size;
    return 1;
}

fastf_t
rt_view_initial_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_i_scale : 1.0;
}

int
rt_view_initial_scale_set_bsg(struct bsg_view *v, fastf_t scale)
{
    if (!v)
	return 0;

    v->gv_i_scale = scale;
    return 1;
}

fastf_t
rt_view_absolute_scale_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_a_scale : 0.0;
}

int
rt_view_absolute_scale_set_bsg(struct bsg_view *v, fastf_t scale)
{
    if (!v)
	return 0;

    v->gv_a_scale = scale;
    return 1;
}

fastf_t
rt_view_local2base_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_local2base : 1.0;
}

fastf_t
rt_view_base2local_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_base2local : 1.0;
}

int
rt_view_unit_conversion_set_bsg(struct bsg_view *v,
				fastf_t local2base,
				fastf_t base2local)
{
    if (!v)
	return 0;

    v->gv_local2base = local2base;
    v->gv_base2local = base2local;
    return 1;
}

fastf_t
rt_view_size_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_size : 1.0;
}

int
rt_view_size_set_bsg(struct bsg_view *v, fastf_t size)
{
    if (!v)
	return 0;

    bsg_view_set_size(v, size);
    return 1;
}

fastf_t
rt_view_inverse_size_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_isize : 1.0;
}

int
rt_view_eye_pos_from_bsg(point_t eye_pos, const struct bsg_view *v)
{
    if (!eye_pos)
	return 0;

    if (!v) {
	VSET(eye_pos, 0.0, 0.0, 1.0);
	return 0;
    }

    VMOVE(eye_pos, v->gv_eye_pos);
    return 1;
}

int
rt_view_eye_pos_set_bsg(struct bsg_view *v, const point_t eye_pos)
{
    if (!v || !eye_pos)
	return 0;

    VMOVE(v->gv_eye_pos, eye_pos);
    return 1;
}

int
rt_view_keypoint_from_bsg(point_t keypoint, const struct bsg_view *v)
{
    if (!keypoint)
	return 0;

    if (!v) {
	VSETALL(keypoint, 0.0);
	return 0;
    }

    VMOVE(keypoint, v->gv_keypoint);
    return 1;
}

int
rt_view_keypoint_set_bsg(struct bsg_view *v, const point_t keypoint)
{
    if (!v || !keypoint)
	return 0;

    VMOVE(v->gv_keypoint, keypoint);
    return 1;
}

char
rt_view_rotate_about_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_rotate_about : 'v';
}

int
rt_view_rotate_about_set_bsg(struct bsg_view *v, char rotate_about)
{
    if (!v)
	return 0;

    v->gv_rotate_about = rotate_about;
    return 1;
}

char
rt_view_coord_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_coord : 'v';
}

int
rt_view_coord_set_bsg(struct bsg_view *v, char coord)
{
    if (!v)
	return 0;

    v->gv_coord = coord;
    return 1;
}

int
rt_view_snap_lines_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_lines(v);
}

int
rt_view_snap_lines_set_bsg(struct bsg_view *v, int enabled)
{
    if (!v)
	return 0;

    bsg_view_set_snap_lines(v, enabled);
    return 1;
}

static bsg_feature_ref
rt_view_polygon_ref_to_bsg(rt_view_polygon_ref ref)
{
    bsg_feature_ref bsg_ref = {ref.token, ref.revision};
    return bsg_ref;
}

static rt_view_polygon_ref
rt_view_polygon_ref_from_bsg(bsg_feature_ref ref)
{
    rt_view_polygon_ref rt_ref = {ref.token, ref.revision};
    return rt_ref;
}

rt_view_polygon_ref
rt_view_polygon_create_bsg(struct bsg_view *v, int type, point_t *fp)
{
    if (!v || !fp)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_create_view_polygon_ref(v, type, fp));
}

rt_view_polygon_ref
rt_view_polygon_select_bsg(struct bsg_view *v, point_t *cp)
{
    if (!v || !cp)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_view_select_polygon_ref(v, cp));
}

rt_view_polygon_ref
rt_view_polygon_find_bsg(struct bsg_view *v, const char *name)
{
    if (!v || !name)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_view_polygon_find_ref(v, name));
}

rt_view_polygon_ref
rt_view_polygon_dup_bsg(struct bsg_view *v, const char *name, const char *new_name)
{
    if (!v || !name || !new_name)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_view_polygon_dup_ref(v, name, new_name));
}

static void
rt_view_polygon_record_from_bsg(struct rt_view_polygon_record *record,
				const struct bsg_polygon_record *bsg_record)
{
    if (!record || !bsg_record)
	return;

    memset(record, 0, sizeof(*record));
    record->ref.token = bsg_record->ref.token;
    record->ref.revision = bsg_record->ref.revision;
    record->name = bsg_record->name;
    record->type = bsg_record->type;
    record->fill_flag = bsg_record->fill_flag;
    V2MOVE(record->fill_dir, bsg_record->fill_dir);
    record->fill_delta = bsg_record->fill_delta;
    BU_COLOR_CPY(&record->fill_color, &bsg_record->fill_color);
    record->edge_color[0] = bsg_record->edge_color[0];
    record->edge_color[1] = bsg_record->edge_color[1];
    record->edge_color[2] = bsg_record->edge_color[2];
    record->curr_contour_i = bsg_record->curr_contour_i;
    record->curr_point_i = bsg_record->curr_point_i;
    record->first_contour_open = bsg_record->first_contour_open;
    record->contour_count = bsg_record->contour_count;
    record->point_count = bsg_record->point_count;
    VMOVE(record->origin_point, bsg_record->origin_point);
    HMOVE(record->vp, bsg_record->vp);
    record->vZ = bsg_record->vZ;
    record->user_data = bsg_record->user_data;
}

struct rt_view_polygon_visit_record_state {
    rt_view_polygon_record_callback_bsg_t callback;
    void *data;
};

static int
rt_view_polygon_visit_record_bsg_cb(bsg_feature_ref ref,
				    const struct bsg_polygon_record *bsg_record,
				    void *data)
{
    struct rt_view_polygon_visit_record_state *state =
	(struct rt_view_polygon_visit_record_state *)data;
    if (!state || !state->callback || !bsg_record)
	return 0;

    struct rt_view_polygon_record record;
    rt_view_polygon_record_from_bsg(&record, bsg_record);
    return state->callback(rt_view_polygon_ref_from_bsg(ref), &record, state->data);
}

void
rt_view_polygon_visit_records_bsg(struct bsg_view *v,
				  rt_view_polygon_record_callback_bsg_t callback,
				  void *data)
{
    if (!v || !callback)
	return;

    struct rt_view_polygon_visit_record_state state;
    state.callback = callback;
    state.data = data;
    bsg_view_polygon_visit_records(v, rt_view_polygon_visit_record_bsg_cb,
	    &state);
}

size_t
rt_view_polygon_snap_count_bsg(struct bsg_view *v, rt_view_polygon_ref exclude)
{
    if (!v)
	return 0;

    return bsg_view_polygon_snap_count(v, rt_view_polygon_ref_to_bsg(exclude));
}

int
rt_view_polygon_clear_point_selection_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    return bsg_view_polygon_clear_point_selection(v);
}

int
rt_view_polygon_ref_is_null_bsg(rt_view_polygon_ref ref)
{
    return bsg_polygon_ref_is_null(rt_view_polygon_ref_to_bsg(ref));
}

int
rt_view_polygon_record_get_bsg(rt_view_polygon_ref ref, struct rt_view_polygon_record *record)
{
    if (!record)
	return 0;

    struct bsg_polygon_record bsg_record;
    if (!bsg_polygon_record_get(rt_view_polygon_ref_to_bsg(ref), &bsg_record))
	return 0;

    rt_view_polygon_record_from_bsg(record, &bsg_record);
    return 1;
}

int
rt_view_polygon_update_bsg(rt_view_polygon_ref ref, struct bsg_view *v, int utype)
{
    return bsg_polygon_update(rt_view_polygon_ref_to_bsg(ref), v, utype);
}

int
rt_view_polygon_update_screen_pt_bsg(rt_view_polygon_ref ref, struct bsg_view *v,
				     int x, int y, int utype)
{
    return bsg_polygon_update_screen_pt(rt_view_polygon_ref_to_bsg(ref), v, x, y, utype);
}

int
rt_view_polygon_move_bsg(rt_view_polygon_ref ref, point_t *current_point,
			 point_t *previous_point)
{
    if (!current_point || !previous_point)
	return 0;

    return bsg_polygon_move(rt_view_polygon_ref_to_bsg(ref), current_point, previous_point);
}

int
rt_view_polygon_set_name_bsg(rt_view_polygon_ref ref, const char *name)
{
    if (!name)
	return 0;

    return bsg_polygon_set_name(rt_view_polygon_ref_to_bsg(ref), name);
}

int
rt_view_polygon_set_view_bsg(rt_view_polygon_ref ref, struct bsg_view *v)
{
    return bsg_polygon_set_view(rt_view_polygon_ref_to_bsg(ref), v);
}

int
rt_view_polygon_set_visual_bsg(rt_view_polygon_ref ref,
			       const struct bu_color *edge_color,
			       const struct bu_color *fill_color,
			       fastf_t fill_slope_x,
			       fastf_t fill_slope_y,
			       fastf_t fill_density,
			       fastf_t vZ,
			       int fill_flag)
{
    return bsg_polygon_set_visual(rt_view_polygon_ref_to_bsg(ref), edge_color, fill_color, fill_slope_x,
	    fill_slope_y, fill_density, vZ, fill_flag);
}

int
rt_view_polygon_set_open_bsg(rt_view_polygon_ref ref, int open)
{
    return bsg_polygon_set_open(rt_view_polygon_ref_to_bsg(ref), open);
}

int
rt_view_polygon_close_bsg(rt_view_polygon_ref ref)
{
    return bsg_polygon_close(rt_view_polygon_ref_to_bsg(ref));
}

int
rt_view_polygon_clear_selected_point_bsg(rt_view_polygon_ref ref)
{
    return bsg_polygon_clear_selected_point(rt_view_polygon_ref_to_bsg(ref));
}

int
rt_view_polygon_remove_bsg(rt_view_polygon_ref ref)
{
    return bsg_polygon_remove(rt_view_polygon_ref_to_bsg(ref));
}

void *
rt_view_polygon_user_data_bsg(rt_view_polygon_ref ref)
{
    return bsg_polygon_user_data(rt_view_polygon_ref_to_bsg(ref));
}

int
rt_view_polygon_user_data_set_bsg(rt_view_polygon_ref ref, void *user_data)
{
    return bsg_polygon_user_data_set(rt_view_polygon_ref_to_bsg(ref), user_data);
}

int
rt_view_polygon_csg_bsg(rt_view_polygon_ref target, rt_view_polygon_ref stencil, bg_clip_t op)
{
    return bsg_polygon_csg_ref(rt_view_polygon_ref_to_bsg(target),
	    rt_view_polygon_ref_to_bsg(stencil), op);
}

rt_view_polygon_ref
rt_view_polygon_import_sketch_bsg(const char *name, struct db_i *dbip,
				  struct directory *dp, struct bsg_view *v)
{
    if (!name || !dbip || !dp || !v)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(db_sketch_to_view_polygon_ref(name, dbip, dp, v));
}

struct directory *
rt_view_polygon_export_sketch_bsg(struct db_i *dbip, const char *name,
				  rt_view_polygon_ref ref)
{
    if (!dbip || !name)
	return NULL;

    return db_view_polygon_ref_to_sketch(dbip, name, rt_view_polygon_ref_to_bsg(ref));
}

int
rt_view_polygon_snap_exclude_set_bsg(struct bsg_view *v, rt_view_polygon_ref ref)
{
    return rt_view_snap_exclude_feature_set_bsg(v, rt_view_polygon_ref_to_bsg(ref));
}

int
rt_view_snap_source_flags_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_source_flags(v);
}

int
rt_view_snap_source_flags_set_bsg(struct bsg_view *v, int flags)
{
    if (!v)
	return 0;

    bsg_view_set_snap_source_flags(v, flags);
    return 1;
}

unsigned long long
rt_view_snap_kind_mask_from_bsg(const struct bsg_view *v)
{
    return (unsigned long long)bsg_view_snap_kind_mask(v);
}

int
rt_view_snap_exclude_feature_set_bsg(struct bsg_view *v, bsg_feature_ref ref)
{
    if (!v)
	return 0;

    bsg_view_snap_exclude_feature_set(v, ref);
    return 1;
}

int
rt_view_snap_exclude_feature_clear_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_snap_exclude_feature_clear(v);
    return 1;
}

unsigned long long
rt_view_prepare_tcl_snap_bsg(struct bsg_view *v)
{
    return (unsigned long long)bsg_view_prepare_tcl_snap(v);
}

int
rt_view_center_linesnap_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    bsg_view_center_linesnap(v);
    return 1;
}

int
rt_view_zclip_from_bsg(const struct bsg_view *v)
{
    return bsg_view_zclip(v);
}

int
rt_view_zclip_set_bsg(struct bsg_view *v, int zclip)
{
    if (!v)
	return 0;

    bsg_view_set_zclip(v, zclip);
    return 1;
}

int
rt_view_framebuffer_mode_from_bsg(const struct bsg_view *v)
{
    return bsg_view_framebuffer_mode(v);
}

int
rt_view_framebuffer_mode_set_bsg(struct bsg_view *v, int mode)
{
    if (!v)
	return 0;

    bsg_view_set_framebuffer_mode(v, mode);
    return 1;
}

int
rt_view_cleared_from_bsg(const struct bsg_view *v)
{
    return bsg_view_cleared(v);
}

int
rt_view_cleared_set_bsg(struct bsg_view *v, int cleared)
{
    if (!v)
	return 0;

    bsg_view_set_cleared(v, cleared);
    return 1;
}

int
rt_view_settings_shared_bsg(const struct bsg_view *a, const struct bsg_view *b)
{
    return bsg_view_settings_shared(a, b);
}

double
rt_view_snap_tolerance_factor_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_tolerance_factor(v);
}

int
rt_view_snap_tolerance_factor_set_bsg(struct bsg_view *v, double factor)
{
    if (!v)
	return 0;

    bsg_view_set_snap_tolerance_factor(v, factor);
    return 1;
}

int
rt_mesh_lod_load_view_scene_ref_bsg(struct rt_mesh_lod *lod,
				    bsg_scene_ref visibility_ref,
				    struct bsg_view *v,
				    int reset)
{
    struct bsg_mesh_lod *bsg_lod = _rt_mesh_lod_bsg(lod);
    if (!bsg_lod)
	return -1;

    return bsg_mesh_lod_load_view_scene_ref(bsg_lod, visibility_ref, v, reset);
}

void
rt_mesh_lod_free_scene_ref_bsg(bsg_scene_ref ref)
{
    bsg_mesh_lod_free_scene_ref(ref);
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
