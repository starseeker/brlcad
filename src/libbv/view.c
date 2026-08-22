/*                    B V / V I E W . C
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
/** @file libbv/view.c
 *
 * Scene-graph-free view mechanics.
 */

#include "common.h"

#include <math.h>
#include <string.h>

#include "vmath.h"
#include "bg/plane.h"
#include "bn/mat.h"
#include "bn/qmath.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bv/view.h"

struct bv_context_callback_record {
    bv_context_callback_t callback;
    void *client_data;
};

struct bv_context_owner_data_record {
    const void *key;
    void *data;
};

static void
bv_mat_aet(struct bv *v)
{
    mat_t tmat;
    fastf_t twist;
    fastf_t c_twist;
    fastf_t s_twist;

    if (!v)
	return;

    bn_mat_angles(v->rotation, 270.0 + v->aet[Y], 0.0, 270.0 - v->aet[X]);
    twist = -v->aet[Z] * DEG2RAD;
    c_twist = cos(twist);
    s_twist = sin(twist);
    bn_mat_zrot(tmat, s_twist, c_twist);
    bn_mat_mul2(tmat, v->rotation);
}

struct bv *
bv_create(void)
{
    struct bv *v = NULL;
    BU_ALLOC(v, struct bv);
    bv_init(v);
    return v;
}

void
bv_destroy(struct bv *v)
{
    if (!v)
	return;

    bv_free(v);
    bu_free(v, "bv");
}

void
bv_init(struct bv *v)
{
    if (!v)
	return;

    memset(v, 0, sizeof(*v));
    v->magic = BV_STATE_MAGIC;
    bu_vls_init(&v->name);
    v->scale = BV_DEFAULT_SCALE;
    v->initial_scale = v->scale;
    v->absolute_scale = 0.0;
    v->size = 2.0 * v->scale;
    v->inverse_size = 1.0 / v->size;
    v->local2base = 1.0;
    v->base2local = 1.0;
    v->refresh_enabled = 1;
    bv_lod_policy_init(&v->lod_policy);
    VSET(v->aet, 35.0, 25.0, 0.0);
    VSET(v->eye_pos, 0.0, 0.0, 1.0);
    v->coord_mode = 'v';
    v->rotate_about = 'v';
    v->min_mouse_delta = -20.0;
    v->max_mouse_delta = 20.0;
    v->rotate_scale = 0.4;
    v->scale_scale = 2.0;
    MAT_IDN(v->rotation);
    MAT_IDN(v->center);
    MAT_IDN(v->model2view);
    MAT_IDN(v->pmodel2view);
    MAT_IDN(v->view2model);
    MAT_IDN(v->pmat);
    v->radius = 1.0;
    bv_knobs_reset(&v->knobs, BV_KNOBS_ALL);
    bv_faceplate_defaults(v);
    bv_snap_defaults(v);
    bv_mat_aet(v);
    bv_update(v);
}

void
bv_free(struct bv *v)
{
    if (!v)
	return;

    bu_vls_free(&v->name);
    v->magic = 0;
}

int
bv_is_valid(const struct bv *v)
{
    return (v && v->magic == BV_STATE_MAGIC);
}

void
bv_view_info_init(struct bv_view_info *info)
{
    struct bv_view_info defaults = BV_VIEW_INFO_INIT;
    if (info)
	*info = defaults;
}

void
bv_view_info_sanitize(struct bv_view_info *info)
{
    if (!info)
	return;

    if (info->width <= 0)
	info->width = 1;
    if (info->height <= 0)
	info->height = 1;
    if (info->size <= SMALL_FASTF)
	info->size = 1.0;
    if (info->lod.scale <= SMALL_FASTF)
	info->lod.scale = 1.0;
    if (info->lod.curve_scale <= SMALL_FASTF)
	info->lod.curve_scale = 1.0;
    if (info->lod.point_scale <= SMALL_FASTF)
	info->lod.point_scale = 1.0;
}

void
bv_lod_policy_init(struct bv_lod_policy *policy)
{
    struct bv_lod_policy defaults = BV_LOD_POLICY_INIT;
    if (policy)
	*policy = defaults;
}

void
bv_lod_policy_sanitize(struct bv_lod_policy *policy)
{
    if (!policy)
	return;

    if (policy->scale <= SMALL_FASTF)
	policy->scale = 1.0;
    if (policy->curve_scale <= SMALL_FASTF)
	policy->curve_scale = 1.0;
    if (policy->point_scale <= SMALL_FASTF)
	policy->point_scale = 1.0;
}

int
bv_view_lod_policy_get(struct bv_lod_policy *policy, const struct bv *v)
{
    if (!policy || !bv_is_valid(v))
	return 0;

    *policy = v->lod_policy;
    bv_lod_policy_sanitize(policy);
    return 1;
}

int
bv_view_lod_policy_set(struct bv *v, const struct bv_lod_policy *policy)
{
    struct bv_lod_policy sanitized;

    if (!bv_is_valid(v) || !policy)
	return 0;

    sanitized = *policy;
    bv_lod_policy_sanitize(&sanitized);
    if (v->lod_policy.policy == sanitized.policy &&
	v->lod_policy.forced_level == sanitized.forced_level &&
	v->lod_policy.mesh_enabled == sanitized.mesh_enabled &&
	v->lod_policy.csg_enabled == sanitized.csg_enabled &&
	v->lod_policy.zoom_refresh == sanitized.zoom_refresh &&
	v->lod_policy.bot_threshold == sanitized.bot_threshold &&
	EQUAL(v->lod_policy.scale, sanitized.scale) &&
	EQUAL(v->lod_policy.curve_scale, sanitized.curve_scale) &&
	EQUAL(v->lod_policy.point_scale, sanitized.point_scale))
	return 1;

    v->lod_policy = sanitized;
    v->frame_revision++;
    (void)bv_refresh_request(v, BV_REFRESH_DRAW);
    return 1;
}

static struct bv_lod_settings
bv_view_lod_policy(const struct bv_view_info *info)
{
    struct bv_lod_settings policy = BV_LOD_SETTINGS_INIT;
    if (info)
	policy = info->lod;
    if (policy.curve_scale <= SMALL_FASTF)
	policy.curve_scale = 1.0;
    if (policy.point_scale <= SMALL_FASTF)
	policy.point_scale = 1.0;
    if (policy.scale <= SMALL_FASTF)
	policy.scale = 1.0;
    return policy;
}

fastf_t
bv_view_lod_curve_scale(const struct bv_view_info *info)
{
    return bv_view_lod_policy(info).curve_scale;
}

size_t
bv_view_lod_bot_threshold(const struct bv_view_info *info)
{
    return bv_view_lod_policy(info).bot_threshold;
}

static fastf_t
bv_view_avg_size(const struct bv_view_info *info)
{
    fastf_t view_aspect, x_size, y_size;

    if (!info || info->width <= 0 || info->height <= 0 || info->size <= SMALL_FASTF)
	return 1.0;

    view_aspect = (fastf_t)info->width / info->height;
    x_size = info->size;
    y_size = x_size / view_aspect;

    return (x_size + y_size) / 2.0;
}

fastf_t
bv_view_avg_sample_spacing(const struct bv_view_info *info)
{
    fastf_t avg_view_size, avg_view_samples;

    if (!info || info->width <= 0 || info->height <= 0)
	return 1.0;

    avg_view_size = bv_view_avg_size(info);
    avg_view_samples = (info->width + info->height) / 2.0;

    return avg_view_size / avg_view_samples;
}

fastf_t
bv_view_solid_point_spacing(const struct bv_view_info *info, fastf_t solid_width)
{
    fastf_t radius, avg_view_size, avg_sample_spacing;
    point2d_t p1, p2;

    if (solid_width < SQRT_SMALL_FASTF)
	solid_width = SQRT_SMALL_FASTF;

    avg_view_size = bv_view_avg_size(info);

    radius = solid_width / 4.0;
    if (avg_view_size < solid_width)
	radius = avg_view_size / 4.0;

    p1[Y] = radius;
    p1[X] = 0.0;

    avg_sample_spacing = bv_view_avg_sample_spacing(info);
    if (avg_sample_spacing < radius)
	p2[Y] = radius - avg_sample_spacing;
    else
	p2[Y] = radius;
    p2[X] = sqrt((radius * radius) - (p2[Y] * p2[Y]));

    return DIST_PNT2_PNT2(p1, p2) / bv_view_lod_policy(info).point_scale;
}

int
bv_copy(struct bv *dst, const struct bv *src)
{
    if (!dst || !bv_is_valid(src))
	return 0;

    if (dst == src)
	return 1;

    if (bv_is_valid(dst))
	bv_free(dst);

    bv_init(dst);
    bu_vls_sprintf(&dst->name, "%s", bu_vls_cstr(&src->name));
    dst->user_data = src->user_data;
    dst->width = src->width;
    dst->height = src->height;
    dst->scale = src->scale;
    dst->initial_scale = src->initial_scale;
    dst->absolute_scale = src->absolute_scale;
    dst->size = src->size;
    dst->inverse_size = src->inverse_size;
    dst->perspective = src->perspective;
    dst->local2base = src->local2base;
    dst->base2local = src->base2local;
    dst->frame_revision = src->frame_revision;
    dst->refresh_dirty = src->refresh_dirty;
    dst->refresh_enabled = src->refresh_enabled;
    dst->refresh_suppressed = src->refresh_suppressed;
    dst->refresh_drawn_count = src->refresh_drawn_count;
    dst->frametime = src->frametime;
    dst->zclip = src->zclip;
    dst->framebuffer_mode = src->framebuffer_mode;
    dst->cleared = src->cleared;
    dst->lod_policy = src->lod_policy;
    VMOVE(dst->aet, src->aet);
    VMOVE(dst->eye_pos, src->eye_pos);
    VMOVE(dst->keypoint, src->keypoint);
    VMOVE(dst->current_point, src->current_point);
    dst->mouse_x = src->mouse_x;
    dst->mouse_y = src->mouse_y;
    dst->previous_mouse_x = src->previous_mouse_x;
    dst->previous_mouse_y = src->previous_mouse_y;
    dst->coord_mode = src->coord_mode;
    dst->rotate_about = src->rotate_about;
    dst->min_mouse_delta = src->min_mouse_delta;
    dst->max_mouse_delta = src->max_mouse_delta;
    dst->rotate_scale = src->rotate_scale;
    dst->scale_scale = src->scale_scale;
    MAT_COPY(dst->rotation, src->rotation);
    MAT_COPY(dst->center, src->center);
    MAT_COPY(dst->model2view, src->model2view);
    MAT_COPY(dst->pmodel2view, src->pmodel2view);
    MAT_COPY(dst->view2model, src->view2model);
    MAT_COPY(dst->pmat, src->pmat);
    VMOVE(dst->obb_center, src->obb_center);
    VMOVE(dst->obb_extent1, src->obb_extent1);
    VMOVE(dst->obb_extent2, src->obb_extent2);
    VMOVE(dst->obb_extent3, src->obb_extent3);
    dst->radius = src->radius;
    dst->knobs = src->knobs;
    dst->interactive_rect = src->interactive_rect;
    dst->adc = src->adc;
    dst->grid = src->grid;
    dst->model_axes = src->model_axes;
    dst->view_axes = src->view_axes;
    dst->center_dot = src->center_dot;
    dst->scale_overlay = src->scale_overlay;
    dst->params = src->params;
    dst->background = src->background;
    dst->snap = src->snap;
    return 1;
}

int
bv_name_set(struct bv *v, const char *name)
{
    if (!bv_is_valid(v))
	return 0;

    bu_vls_sprintf(&v->name, "%s", name ? name : "");
    return 1;
}

const char *
bv_name_get(const struct bv *v)
{
    return (bv_is_valid(v) && bu_vls_strlen(&v->name)) ?
	bu_vls_cstr(&v->name) : NULL;
}

int
bv_user_data_set(struct bv *v, void *user_data)
{
    if (!bv_is_valid(v))
	return 0;

    v->user_data = user_data;
    return 1;
}

void *
bv_user_data_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->user_data : NULL;
}

int
bv_dimensions_set(struct bv *v, int width, int height)
{
    if (!bv_is_valid(v))
	return 0;

    v->width = width;
    v->height = height;
    return 1;
}

int
bv_width_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->width : 0;
}

int
bv_height_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->height : 0;
}

fastf_t
bv_scale_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->scale : 1.0;
}

fastf_t *
bv_scale_storage_get(struct bv *v)
{
    return bv_is_valid(v) ? &v->scale : NULL;
}

fastf_t
bv_initial_scale_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->initial_scale : 1.0;
}

int
bv_initial_scale_set(struct bv *v, fastf_t scale)
{
    if (!bv_is_valid(v))
	return 0;

    v->initial_scale = scale;
    v->frame_revision++;
    return 1;
}

fastf_t
bv_absolute_scale_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->absolute_scale : 0.0;
}

int
bv_absolute_scale_set(struct bv *v, fastf_t scale)
{
    if (!bv_is_valid(v))
	return 0;

    v->absolute_scale = scale;
    v->frame_revision++;
    return 1;
}

fastf_t
bv_size_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->size : 0.0;
}

fastf_t
bv_inverse_size_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->inverse_size : 1.0;
}

fastf_t
bv_radius_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->radius : 1.0;
}

fastf_t
bv_perspective_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->perspective : 0.0;
}

int
bv_perspective_set(struct bv *v, fastf_t perspective)
{
    if (!bv_is_valid(v))
	return 0;

    v->perspective = perspective;
    v->frame_revision++;
    return 1;
}

fastf_t
bv_local2base_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->local2base : 1.0;
}

fastf_t
bv_base2local_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->base2local : 1.0;
}

int
bv_unit_conversion_set(struct bv *v, fastf_t local2base, fastf_t base2local)
{
    if (!bv_is_valid(v))
	return 0;

    v->local2base = local2base;
    v->base2local = base2local;
    v->frame_revision++;
    return 1;
}

char
bv_coord_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->coord_mode : 'v';
}

int
bv_coord_set(struct bv *v, char coord)
{
    if (!bv_is_valid(v))
	return 0;

    v->coord_mode = coord;
    v->frame_revision++;
    return 1;
}

char
bv_rotate_about_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->rotate_about : 'v';
}

int
bv_rotate_about_set(struct bv *v, char rotate_about)
{
    if (!bv_is_valid(v))
	return 0;

    v->rotate_about = rotate_about;
    v->frame_revision++;
    return 1;
}

int
bv_refresh_request(struct bv *v, uint32_t flags)
{
    if (!bv_is_valid(v))
	return 0;

    if (v->refresh_enabled && !v->refresh_suppressed)
	v->refresh_dirty |= flags;
    return 1;
}

int
bv_refresh_dirty_get(const struct bv *v)
{
    return (bv_is_valid(v) && v->refresh_dirty) ? 1 : 0;
}

uint32_t
bv_refresh_consume(struct bv *v)
{
    uint32_t dirty;

    if (!bv_is_valid(v))
	return 0;

    dirty = v->refresh_dirty;
    v->refresh_dirty = 0;
    return dirty;
}

int
bv_refresh_complete(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_dirty = 0;
    return 1;
}

int
bv_refresh_enabled_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->refresh_enabled : 0;
}

int
bv_refresh_enabled_set(struct bv *v, int enabled)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_enabled = enabled ? 1 : 0;
    return 1;
}

int
bv_refresh_suppressed_get(const struct bv *v)
{
    return (bv_is_valid(v) && v->refresh_suppressed > 0) ? 1 : 0;
}

int
bv_refresh_suppress_begin(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_suppressed++;
    return 1;
}

int
bv_refresh_suppress_end(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    if (v->refresh_suppressed > 0)
	v->refresh_suppressed--;
    return 1;
}

int
bv_refresh_drawn_count_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->refresh_drawn_count : 0;
}

int
bv_refresh_drawn_count_set(struct bv *v, int count)
{
    if (!bv_is_valid(v))
	return 0;

    v->refresh_drawn_count = count;
    return 1;
}

uint64_t
bv_frame_revision_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->frame_revision : 0;
}

uint64_t
bv_frame_revision_bump(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    v->frame_revision++;
    return v->frame_revision;
}

unsigned long long
bv_hash(const struct bv *v)
{
    struct bu_data_hash_state *state;
    unsigned long long hv;

    if (!bv_is_valid(v))
	return 0ULL;

    state = bu_data_hash_create();
    if (!state)
	return 0ULL;

    bu_data_hash_update(state, &v->width, sizeof(v->width));
    bu_data_hash_update(state, &v->height, sizeof(v->height));
    bu_data_hash_update(state, &v->scale, sizeof(v->scale));
    bu_data_hash_update(state, &v->size, sizeof(v->size));
    bu_data_hash_update(state, &v->perspective, sizeof(v->perspective));
    bu_data_hash_update(state, &v->aet, sizeof(v->aet));
    bu_data_hash_update(state, &v->rotation, sizeof(v->rotation));
    bu_data_hash_update(state, &v->center, sizeof(v->center));
    bu_data_hash_update(state, &v->model2view, sizeof(v->model2view));
    bu_data_hash_update(state, &v->lod_policy.policy, sizeof(v->lod_policy.policy));
    bu_data_hash_update(state, &v->lod_policy.forced_level, sizeof(v->lod_policy.forced_level));
    bu_data_hash_update(state, &v->lod_policy.mesh_enabled, sizeof(v->lod_policy.mesh_enabled));
    bu_data_hash_update(state, &v->lod_policy.csg_enabled, sizeof(v->lod_policy.csg_enabled));
    bu_data_hash_update(state, &v->lod_policy.zoom_refresh, sizeof(v->lod_policy.zoom_refresh));
    bu_data_hash_update(state, &v->lod_policy.bot_threshold, sizeof(v->lod_policy.bot_threshold));
    bu_data_hash_update(state, &v->lod_policy.scale, sizeof(v->lod_policy.scale));
    bu_data_hash_update(state, &v->lod_policy.curve_scale, sizeof(v->lod_policy.curve_scale));
    bu_data_hash_update(state, &v->lod_policy.point_scale, sizeof(v->lod_policy.point_scale));
    hv = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return hv;
}

int
bv_frametime_set(struct bv *v, uint64_t frametime)
{
    if (!bv_is_valid(v))
	return 0;

    v->frametime = frametime;
    return 1;
}

uint64_t
bv_frametime_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->frametime : 0;
}

int
bv_zclip_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->zclip : 0;
}

int
bv_zclip_set(struct bv *v, int zclip)
{
    if (!bv_is_valid(v))
	return 0;

    v->zclip = zclip ? 1 : 0;
    return 1;
}

int
bv_framebuffer_mode_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->framebuffer_mode : 0;
}

int
bv_framebuffer_mode_set(struct bv *v, int mode)
{
    if (!bv_is_valid(v))
	return 0;

    v->framebuffer_mode = mode;
    return 1;
}

int
bv_cleared_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->cleared : 0;
}

int
bv_cleared_set(struct bv *v, int cleared)
{
    if (!bv_is_valid(v))
	return 0;

    v->cleared = cleared ? 1 : 0;
    return 1;
}

void
bv_faceplate_defaults(struct bv *v)
{
    struct bv_interactive_rect_state rect = BV_INTERACTIVE_RECT_STATE_INIT;
    struct bv_adc_state adc = BV_ADC_STATE_INIT;
    struct bv_grid_state grid = BV_GRID_STATE_INIT;
    struct bv_lighting_state lighting = BV_LIGHTING_STATE_INIT;
    struct bv_shading_state shading = BV_SHADING_STATE_INIT;
    struct bv_axes_state model_axes = BV_AXES_STATE_INIT;
    struct bv_axes_state view_axes = BV_AXES_STATE_INIT;
    struct bv_other_state center_dot = BV_OTHER_STATE_INIT;
    struct bv_other_state scale_overlay = BV_OTHER_STATE_INIT;
    struct bv_params_state params = BV_PARAMS_STATE_INIT;

    if (!bv_is_valid(v))
	return;

    grid.res_h = 1.0;
    grid.res_v = 1.0;
    grid.res_major_h = 5;
    grid.res_major_v = 5;
    VSET(grid.color, 255, 255, 255);

    VSET(adc.line_color, 255, 255, 0);
    VSET(adc.tick_color, 255, 255, 255);
    adc.line_width = 1;

    VSET(view_axes.axes_pos, 0.80, -0.80, 0.0);
    view_axes.axes_size = 0.2;
    view_axes.pos_only = 1;
    view_axes.label_flag = 1;
    view_axes.triple_color = 1;
    VSET(view_axes.axes_color, 255, 255, 255);
    VSET(view_axes.label_color, 255, 255, 0);

    model_axes.axes_size = 2.0;
    model_axes.label_flag = 1;
    model_axes.tick_enabled = 1;
    model_axes.tick_length = 4;
    model_axes.tick_major_length = 8;
    model_axes.tick_interval = 100.0;
    model_axes.ticks_per_major = 10;
    model_axes.tick_threshold = 8;
    VSET(model_axes.axes_color, 255, 255, 255);
    VSET(model_axes.label_color, 255, 255, 0);
    VSET(model_axes.tick_color, 255, 255, 0);
    VSET(model_axes.tick_major_color, 255, 0, 0);

    center_dot.gos_font_size = 20;
    VSET(center_dot.gos_line_color, 255, 255, 0);

    scale_overlay.gos_font_size = 20;
    VSET(scale_overlay.gos_line_color, 255, 255, 0);
    VSET(scale_overlay.gos_text_color, 255, 255, 0);

    /* FPS is part of the default faceplate, so the parameter group itself
     * must be visible as well as its individual FPS flag. */
    params.draw = 1;
    params.draw_size = 1;
    params.draw_center = 1;
    params.draw_az = 1;
    params.draw_el = 1;
    params.draw_tw = 1;
    params.draw_fps = 1;
    VSET(params.color, 255, 255, 0);
    params.font_size = 20;

    v->interactive_rect = rect;
    v->adc = adc;
    v->grid = grid;
    v->lighting = lighting;
    v->shading = shading;
    v->model_axes = model_axes;
    v->view_axes = view_axes;
    v->center_dot = center_dot;
    v->scale_overlay = scale_overlay;
    v->params = params;
}

void
bv_snap_defaults(struct bv *v)
{
    struct bv_snap_state snap = BV_SNAP_STATE_INIT;

    if (!bv_is_valid(v))
	return;

    v->snap = snap;
}

#define BV_STATE_GET(_type, _init, _member) \
int \
bv_ ## _member ## _state_get(struct _type *record, const struct bv *v) \
{ \
    if (!record) \
	return 0; \
    if (!bv_is_valid(v)) { \
	struct _type zero = _init; \
	*record = zero; \
	return 0; \
    } \
    *record = v->_member; \
    return 1; \
}

#define BV_STATE_SET(_type, _member) \
int \
bv_ ## _member ## _state_set(struct bv *v, const struct _type *record) \
{ \
    if (!bv_is_valid(v) || !record) \
	return 0; \
    v->_member = *record; \
    return 1; \
}

BV_STATE_GET(bv_adc_state, BV_ADC_STATE_INIT, adc)
BV_STATE_SET(bv_adc_state, adc)

void
bv_adc_model_to_view(struct bv_adc_state *adcs, mat_t model2view, fastf_t amax)
{
    if (!adcs || !model2view)
	return;

    MAT4X3PNT(adcs->pos_view, model2view, adcs->pos_model);
    adcs->dv_x = adcs->pos_view[X] * amax;
    adcs->dv_y = adcs->pos_view[Y] * amax;
}

void
bv_adc_grid_to_view(struct bv_adc_state *adcs, mat_t model2view, fastf_t amax)
{
    point_t model_pt = VINIT_ZERO;
    point_t view_pt;

    if (!adcs || !model2view)
	return;

    MAT4X3PNT(view_pt, model2view, model_pt);
    VADD2(adcs->pos_view, view_pt, adcs->pos_grid);
    adcs->dv_x = adcs->pos_view[X] * amax;
    adcs->dv_y = adcs->pos_view[Y] * amax;
}

void
bv_adc_view_to_grid(struct bv_adc_state *adcs, mat_t model2view)
{
    point_t model_pt = VINIT_ZERO;
    point_t view_pt;

    if (!adcs || !model2view)
	return;

    MAT4X3PNT(view_pt, model2view, model_pt);
    VSUB2(adcs->pos_grid, adcs->pos_view, view_pt);
}

void
bv_adc_reset(struct bv_adc_state *adcs, mat_t view2model, mat_t model2view)
{
    if (!adcs || !view2model || !model2view)
	return;

    adcs->dv_x = adcs->dv_y = 0;
    adcs->dv_a1 = adcs->dv_a2 = 0;
    adcs->dv_dist = 0;

    VSETALL(adcs->pos_view, 0.0);
    MAT4X3PNT(adcs->pos_model, view2model, adcs->pos_view);
    adcs->dst = (adcs->dv_dist * BV_INV_VIEW + 1.0) * M_SQRT1_2;
    adcs->a1 = adcs->a2 = 45.0;
    bv_adc_view_to_grid(adcs, model2view);

    VSETALL(adcs->anchor_pt_a1, 0.0);
    VSETALL(adcs->anchor_pt_a2, 0.0);
    VSETALL(adcs->anchor_pt_dst, 0.0);

    adcs->anchor_pos = 0;
    adcs->anchor_a1 = 0;
    adcs->anchor_a2 = 0;
    adcs->anchor_dst = 0;
}

BV_STATE_GET(bv_grid_state, BV_GRID_STATE_INIT, grid)
BV_STATE_SET(bv_grid_state, grid)
BV_STATE_GET(bv_lighting_state, BV_LIGHTING_STATE_INIT, lighting)
BV_STATE_SET(bv_lighting_state, lighting)
BV_STATE_GET(bv_shading_state, BV_SHADING_STATE_INIT, shading)
BV_STATE_SET(bv_shading_state, shading)
BV_STATE_GET(bv_params_state, BV_PARAMS_STATE_INIT, params)
BV_STATE_SET(bv_params_state, params)
BV_STATE_GET(bv_background_state, BV_BACKGROUND_STATE_INIT, background)
BV_STATE_SET(bv_background_state, background)

#undef BV_STATE_GET
#undef BV_STATE_SET

int
bv_interactive_rect_state_get(struct bv_interactive_rect_state *record, const struct bv *v)
{
    if (!record)
	return 0;
    if (!bv_is_valid(v)) {
	struct bv_interactive_rect_state zero = BV_INTERACTIVE_RECT_STATE_INIT;
	*record = zero;
	return 0;
    }
    *record = v->interactive_rect;
    return 1;
}

int
bv_interactive_rect_state_set(struct bv *v, const struct bv_interactive_rect_state *record)
{
    if (!bv_is_valid(v) || !record)
	return 0;
    v->interactive_rect = *record;
    return 1;
}

int
bv_model_axes_state_get(struct bv_axes_state *record, const struct bv *v)
{
    if (!record)
	return 0;
    if (!bv_is_valid(v)) {
	struct bv_axes_state zero = BV_AXES_STATE_INIT;
	*record = zero;
	return 0;
    }
    *record = v->model_axes;
    return 1;
}

int
bv_model_axes_state_set(struct bv *v, const struct bv_axes_state *record)
{
    if (!bv_is_valid(v) || !record)
	return 0;
    v->model_axes = *record;
    return 1;
}

int
bv_view_axes_state_get(struct bv_axes_state *record, const struct bv *v)
{
    if (!record)
	return 0;
    if (!bv_is_valid(v)) {
	struct bv_axes_state zero = BV_AXES_STATE_INIT;
	*record = zero;
	return 0;
    }
    *record = v->view_axes;
    return 1;
}

int
bv_view_axes_state_set(struct bv *v, const struct bv_axes_state *record)
{
    if (!bv_is_valid(v) || !record)
	return 0;
    v->view_axes = *record;
    return 1;
}

int
bv_center_dot_state_get(struct bv_other_state *record, const struct bv *v)
{
    if (!record)
	return 0;
    if (!bv_is_valid(v)) {
	struct bv_other_state zero = BV_OTHER_STATE_INIT;
	*record = zero;
	return 0;
    }
    *record = v->center_dot;
    return 1;
}

int
bv_center_dot_state_set(struct bv *v, const struct bv_other_state *record)
{
    if (!bv_is_valid(v) || !record)
	return 0;
    v->center_dot = *record;
    return 1;
}

int
bv_scale_overlay_state_get(struct bv_other_state *record, const struct bv *v)
{
    if (!record)
	return 0;
    if (!bv_is_valid(v)) {
	struct bv_other_state zero = BV_OTHER_STATE_INIT;
	*record = zero;
	return 0;
    }
    *record = v->scale_overlay;
    return 1;
}

int
bv_scale_overlay_state_set(struct bv *v, const struct bv_other_state *record)
{
    if (!bv_is_valid(v) || !record)
	return 0;
    v->scale_overlay = *record;
    return 1;
}

int
bv_snap_state_get(struct bv_snap_state *record, const struct bv *v)
{
    if (!record)
	return 0;
    if (!bv_is_valid(v)) {
	struct bv_snap_state zero = BV_SNAP_STATE_INIT;
	*record = zero;
	return 0;
    }
    *record = v->snap;
    return 1;
}

int
bv_snap_state_set(struct bv *v, const struct bv_snap_state *record)
{
    if (!bv_is_valid(v) || !record)
	return 0;

    v->snap = *record;
    v->snap.lines = v->snap.lines ? 1 : 0;
    if (v->snap.tolerance_factor <= 0.0)
	v->snap.tolerance_factor = 1.0;
    return 1;
}

int
bv_snap_lines_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->snap.lines : 0;
}

int
bv_snap_lines_set(struct bv *v, int enabled)
{
    if (!bv_is_valid(v))
	return 0;

    v->snap.lines = enabled ? 1 : 0;
    return 1;
}

int
bv_snap_source_flags_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->snap.source_flags : 0;
}

int
bv_snap_source_flags_set(struct bv *v, int flags)
{
    if (!bv_is_valid(v))
	return 0;

    v->snap.source_flags = flags;
    return 1;
}

unsigned long long
bv_snap_kind_mask_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->snap.kind_mask : 0ULL;
}

int
bv_snap_kind_mask_set(struct bv *v, unsigned long long mask)
{
    if (!bv_is_valid(v))
	return 0;

    v->snap.kind_mask = mask;
    return 1;
}

double
bv_snap_tolerance_factor_get(const struct bv *v)
{
    return bv_is_valid(v) ? v->snap.tolerance_factor : 1.0;
}

int
bv_snap_tolerance_factor_set(struct bv *v, double factor)
{
    if (!bv_is_valid(v))
	return 0;

    v->snap.tolerance_factor = (factor > 0.0) ? factor : 1.0;
    return 1;
}

int
bv_update(struct bv *v)
{
    vect_t work;
    vect_t work1;
    vect_t temp;
    vect_t temp1;

    if (!bv_is_valid(v))
	return 0;

    bn_mat_mul(v->model2view, v->rotation, v->center);
    v->model2view[15] = v->scale;
    bn_mat_inv(v->view2model, v->model2view);

    VSET(work, 0.0, 0.0, 1.0);
    MAT4X3VEC(temp, v->view2model, work);
    VSET(work1, 1.0, 0.0, 0.0);
    MAT4X3VEC(temp1, v->view2model, work1);
    bn_aet_vec(&v->aet[X], &v->aet[Y], &v->aet[Z], temp, temp1,
	    (fastf_t)0.005);

    if ((NEAR_EQUAL(v->aet[Y], 90.0, (fastf_t)0.005) ||
		NEAR_EQUAL(v->aet[Y], -90.0, (fastf_t)0.005)) &&
	    v->aet[X] < 0.0 && !NEAR_ZERO(v->aet[X], (fastf_t)0.005))
	v->aet[X] += 360.0;
    else if (NEAR_ZERO(v->aet[X], (fastf_t)0.005))
	v->aet[X] = 0.0;

    bn_mat_mul(v->pmodel2view, v->pmat, v->model2view);
    v->frame_revision++;
    return 1;
}

int
bv_model2view_get(mat_t model2view, const struct bv *v)
{
    if (!model2view)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(model2view);
	return 0;
    }
    MAT_COPY(model2view, v->model2view);
    return 1;
}

int
bv_model2view_set(struct bv *v, const mat_t model2view)
{
    if (!bv_is_valid(v) || !model2view)
	return 0;

    MAT_COPY(v->model2view, model2view);
    v->frame_revision++;
    return 1;
}

int
bv_view2model_get(mat_t view2model, const struct bv *v)
{
    if (!view2model)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(view2model);
	return 0;
    }
    MAT_COPY(view2model, v->view2model);
    return 1;
}

int
bv_view2model_set(struct bv *v, const mat_t view2model)
{
    if (!bv_is_valid(v) || !view2model)
	return 0;

    MAT_COPY(v->view2model, view2model);
    v->frame_revision++;
    return 1;
}

int
bv_pmodel2view_get(mat_t pmodel2view, const struct bv *v)
{
    if (!pmodel2view)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(pmodel2view);
	return 0;
    }
    MAT_COPY(pmodel2view, v->pmodel2view);
    return 1;
}

int
bv_pmodel2view_set(struct bv *v, const mat_t pmodel2view)
{
    if (!bv_is_valid(v) || !pmodel2view)
	return 0;

    MAT_COPY(v->pmodel2view, pmodel2view);
    v->frame_revision++;
    return 1;
}

int
bv_pmat_get(mat_t pmat, const struct bv *v)
{
    if (!pmat)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(pmat);
	return 0;
    }
    MAT_COPY(pmat, v->pmat);
    return 1;
}

int
bv_pmat_set(struct bv *v, const mat_t pmat)
{
    if (!bv_is_valid(v) || !pmat)
	return 0;

    MAT_COPY(v->pmat, pmat);
    bn_mat_mul(v->pmodel2view, v->pmat, v->model2view);
    v->frame_revision++;
    return 1;
}

int
bv_rotation_get(mat_t rotation, const struct bv *v)
{
    if (!rotation)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(rotation);
	return 0;
    }
    MAT_COPY(rotation, v->rotation);
    return 1;
}

int
bv_rotation_set(struct bv *v, const mat_t rotation)
{
    if (!bv_is_valid(v) || !rotation)
	return 0;

    MAT_COPY(v->rotation, rotation);
    return bv_update(v);
}

int
bv_orientation_quat_get(quat_t orientation, const struct bv *v)
{
    if (!orientation)
	return 0;
    if (!bv_is_valid(v)) {
	mat_t identity;
	MAT_IDN(identity);
	quat_mat2quat(orientation, identity);
	return 0;
    }

    quat_mat2quat(orientation, v->rotation);
    return 1;
}

int
bv_center_mat_get(mat_t center, const struct bv *v)
{
    if (!center)
	return 0;
    if (!bv_is_valid(v)) {
	MAT_IDN(center);
	return 0;
    }
    MAT_COPY(center, v->center);
    return 1;
}

int
bv_center_mat_set(struct bv *v, const mat_t center)
{
    if (!bv_is_valid(v) || !center)
	return 0;

    MAT_COPY(v->center, center);
    v->frame_revision++;
    return 1;
}

int
bv_center_get(point_t center, const struct bv *v)
{
    if (!center)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(center, 0.0);
	return 0;
    }
    MAT_DELTAS_GET_NEG(center, v->center);
    return 1;
}

int
bv_center_set(struct bv *v, const point_t center)
{
    if (!bv_is_valid(v) || !center)
	return 0;

    MAT_DELTAS_VEC_NEG(v->center, center);
    return bv_update(v);
}

int
bv_scale_set(struct bv *v, fastf_t scale)
{
    if (!bv_is_valid(v))
	return 0;

    v->scale = (scale < BV_MIN_SCALE) ? BV_MIN_SCALE : scale;
    v->size = 2.0 * v->scale;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_scale_state_set(struct bv *v,
		   fastf_t scale,
		   fastf_t initial_scale,
		   fastf_t absolute_scale,
		   fastf_t size,
		   fastf_t inverse_size)
{
    if (!bv_is_valid(v))
	return 0;

    v->scale = scale;
    v->initial_scale = initial_scale;
    v->absolute_scale = absolute_scale;
    v->size = size;
    v->inverse_size = inverse_size;
    v->frame_revision++;
    return 1;
}

int
bv_size_set(struct bv *v, fastf_t size)
{
    if (!bv_is_valid(v))
	return 0;

    v->size = (size < BV_MIN_SIZE) ? BV_MIN_SIZE : size;
    v->scale = 0.5 * v->size;
    if (v->scale < BV_MIN_SCALE)
	v->scale = BV_MIN_SCALE;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_aet_get(vect_t aet, const struct bv *v)
{
    if (!aet)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(aet, 0.0);
	return 0;
    }
    VMOVE(aet, v->aet);
    return 1;
}

int
bv_aet_set(struct bv *v, const vect_t aet)
{
    if (!bv_is_valid(v) || !aet)
	return 0;

    VMOVE(v->aet, aet);
    bv_mat_aet(v);
    return bv_update(v);
}

int
bv_aet_state_set(struct bv *v, const vect_t aet)
{
    if (!bv_is_valid(v) || !aet)
	return 0;

    VMOVE(v->aet, aet);
    v->frame_revision++;
    return 1;
}

int
bv_eye_pos_get(point_t eye_pos, const struct bv *v)
{
    if (!eye_pos)
	return 0;
    if (!bv_is_valid(v)) {
	VSET(eye_pos, 0.0, 0.0, 1.0);
	return 0;
    }
    VMOVE(eye_pos, v->eye_pos);
    return 1;
}

int
bv_eye_pos_set(struct bv *v, const point_t eye_pos)
{
    if (!bv_is_valid(v) || !eye_pos)
	return 0;

    VMOVE(v->eye_pos, eye_pos);
    v->frame_revision++;
    return 1;
}

int
bv_keypoint_get(point_t keypoint, const struct bv *v)
{
    if (!keypoint)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(keypoint, 0.0);
	return 0;
    }
    VMOVE(keypoint, v->keypoint);
    return 1;
}

int
bv_keypoint_set(struct bv *v, const point_t keypoint)
{
    if (!bv_is_valid(v) || !keypoint)
	return 0;

    VMOVE(v->keypoint, keypoint);
    v->frame_revision++;
    return 1;
}

int
bv_current_point_get(point_t current_point, const struct bv *v)
{
    if (!current_point)
	return 0;
    if (!bv_is_valid(v)) {
	VSETALL(current_point, 0.0);
	return 0;
    }

    VMOVE(current_point, v->current_point);
    return 1;
}

int
bv_current_point_set(struct bv *v, const point_t current_point)
{
    if (!bv_is_valid(v) || !current_point)
	return 0;

    VMOVE(v->current_point, current_point);
    v->frame_revision++;
    return 1;
}

int
bv_autoview_bounds(struct bv *v, fastf_t scale, const point_t min, const point_t max)
{
    vect_t center = VINIT_ZERO;
    vect_t radial;
    vect_t sqrt_small;
    fastf_t aspect;

    if (!bv_is_valid(v) || !min || !max)
	return 0;

    if (scale < SQRT_SMALL_FASTF)
	scale = 2.0;

    VSETALL(sqrt_small, SQRT_SMALL_FASTF);
    VADD2SCALE(center, max, min, 0.5);
    VSUB2(radial, max, center);
    VMAX(radial, sqrt_small);
    if (VNEAR_ZERO(radial, SQRT_SMALL_FASTF))
	VSETALL(radial, 1.0);

    MAT_IDN(v->center);
    MAT_DELTAS_VEC_NEG(v->center, center);

    /* Autoview is expected to remain valid when the user rotates the view.
     * Fitting only the largest axis of an AABB is not rotation invariant: a
     * later oblique view can project the box diagonal outside the viewport.
     * Historical display-list bounds usually hid this defect by supplying
     * per-object bounding spheres.  Retained scenes correctly publish tight
     * bounds, so establish the sphere contract here where all callers receive
     * it.  Account for the view convention in which size is the horizontal
     * span and the vertical span is size/aspect. */
    v->scale = MAGNITUDE(radial);
    aspect = (v->width > 0 && v->height > 0) ?
	(fastf_t)v->width / (fastf_t)v->height : 1.0;
    if (aspect > 1.0)
	v->scale *= aspect;
    v->size = scale * v->scale;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_screen_to_view(fastf_t *vx, fastf_t *vy, const struct bv *v, fastf_t sx, fastf_t sy)
{
    fastf_t aspect;

    if (vx)
	*vx = 0.0;
    if (vy)
	*vy = 0.0;
    if (!vx || !vy || !bv_is_valid(v) || !v->width || !v->height)
	return 0;

    *vx = sx / (fastf_t)v->width * 2.0 - 1.0;
    aspect = (fastf_t)v->width / (fastf_t)v->height;
    *vy = (sy / (fastf_t)v->height * -2.0 + 1.0) / aspect;
    return 1;
}

int
bv_snap_grid_2d(const struct bv *v, fastf_t *vx, fastf_t *vy)
{
    point_t anchor_view = VINIT_ZERO;
    fastf_t sf = 0.0;
    fastf_t step_h = 0.0;
    fastf_t step_v = 0.0;
    fastf_t qx = 0.0;
    fastf_t qy = 0.0;
    fastf_t ax = 0.0;
    fastf_t ay = 0.0;

    if (!bv_is_valid(v) || !vx || !vy)
	return 0;
    if (ZERO(v->grid.res_h) || ZERO(v->grid.res_v))
	return 0;

    sf = v->scale * v->base2local;
    if (ZERO(sf))
	return 0;

    step_h = v->grid.res_h * v->base2local;
    step_v = v->grid.res_v * v->base2local;
    if (ZERO(step_h) || ZERO(step_v))
	return 0;

    MAT4X3PNT(anchor_view, v->model2view, v->grid.anchor);

    qx = (*vx) * sf;
    qy = (*vy) * sf;
    ax = anchor_view[X] * sf;
    ay = anchor_view[Y] * sf;

    *vx = (ax + round((qx - ax) / step_h) * step_h) / sf;
    *vy = (ay + round((qy - ay) / step_v) * step_v) / sf;
    return 1;
}

int
bv_screen_to_model(point_t model_point, const struct bv *v, fastf_t sx, fastf_t sy)
{
    fastf_t vx = 0.0;
    fastf_t vy = 0.0;
    point_t vpt = VINIT_ZERO;

    if (model_point)
	VSETALL(model_point, 0.0);
    if (!model_point || !bv_screen_to_view(&vx, &vy, v, sx, sy))
	return 0;

    VSET(vpt, vx, vy, 0.0);
    MAT4X3PNT(model_point, v->view2model, vpt);
    return 1;
}

int
bv_plane_get(plane_t *plane, const struct bv *v)
{
    point_t center;
    vect_t normal;

    if (!plane || !bv_is_valid(v))
	return -1;

    MAT_DELTAS_GET_NEG(center, v->center);
    VMOVEN(normal, v->rotation + 8, 3);
    VUNITIZE(normal);
    VSCALE(normal, normal, -1.0);
    return bg_plane_pt_nrml(plane, center, normal);
}

int
bv_previous_mouse_get(fastf_t *x, fastf_t *y, const struct bv *v)
{
    if (x)
	*x = 0.0;
    if (y)
	*y = 0.0;
    if (!x || !y || !bv_is_valid(v))
	return 0;

    *x = v->previous_mouse_x;
    *y = v->previous_mouse_y;
    return 1;
}

int
bv_previous_mouse_set(struct bv *v, fastf_t x, fastf_t y)
{
    if (!bv_is_valid(v))
	return 0;

    v->previous_mouse_x = x;
    v->previous_mouse_y = y;
    return 1;
}

int
bv_mouse_delta_settings_get(struct bv_mouse_delta_settings *settings,
	const struct bv *v)
{
    struct bv_mouse_delta_settings zero = BV_MOUSE_DELTA_SETTINGS_INIT;

    if (!settings)
	return 0;

    *settings = zero;
    if (!bv_is_valid(v))
	return 0;

    settings->min_delta = v->min_mouse_delta;
    settings->max_delta = v->max_mouse_delta;
    settings->rotate_scale = v->rotate_scale;
    settings->scale_scale = v->scale_scale;
    return 1;
}

int
bv_mouse_delta_settings_set(struct bv *v,
	const struct bv_mouse_delta_settings *settings)
{
    if (!bv_is_valid(v) || !settings ||
	!isfinite(settings->min_delta) || !isfinite(settings->max_delta) ||
	!isfinite(settings->rotate_scale) || !isfinite(settings->scale_scale) ||
	settings->min_delta > 0.0 || settings->max_delta < 0.0 ||
	settings->min_delta > settings->max_delta ||
	settings->rotate_scale <= 0.0 || settings->scale_scale <= 0.0)
	return 0;

    v->min_mouse_delta = settings->min_delta;
    v->max_mouse_delta = settings->max_delta;
    v->rotate_scale = settings->rotate_scale;
    v->scale_scale = settings->scale_scale;
    return 1;
}

int
bv_mouse_delta_clamp(fastf_t *dx, fastf_t *dy, const struct bv *v)
{
    struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;

    if (!dx || !dy || !bv_mouse_delta_settings_get(&settings, v))
	return 0;

    if (*dx < settings.min_delta)
	*dx = settings.min_delta;
    else if (*dx > settings.max_delta)
	*dx = settings.max_delta;
    if (*dy < settings.min_delta)
	*dy = settings.min_delta;
    else if (*dy > settings.max_delta)
	*dy = settings.max_delta;
    return 1;
}

static int
bv_adjust_rotation(struct bv *v, fastf_t dx, fastf_t dy,
	const point_t keypoint, fastf_t rotate_scale)
{
    point_t rot_pt;
    point_t new_origin;
    mat_t viewchg;
    mat_t viewchginv;
    point_t new_cent_view;
    point_t new_cent_model;
    mat_t newrot;

    if (!keypoint)
	return 0;
    bn_mat_angles(newrot, dy * rotate_scale, dx * rotate_scale, 0.0);
    MAT4X3PNT(rot_pt, v->model2view, keypoint);
    bn_mat_xform_about_pnt(viewchg, newrot, rot_pt);
    bn_mat_inv(viewchginv, viewchg);
    VSET(new_origin, 0.0, 0.0, 0.0);
    MAT4X3PNT(new_cent_view, viewchginv, new_origin);
    MAT4X3PNT(new_cent_model, v->view2model, new_cent_view);
    MAT_DELTAS_VEC_NEG(v->center, new_cent_model);
    bn_mat_mul2(newrot, v->rotation);
    return bv_update(v);
}

static int
bv_adjust_translation(struct bv *v, fastf_t dx, fastf_t dy)
{
    fastf_t aspect;
    fastf_t fx;
    fastf_t fy;
    vect_t tt;
    point_t delta;
    point_t work;
    point_t vc;
    point_t nvc;

    if (!v->width || !v->height)
	return 0;
    aspect = (fastf_t)v->width / (fastf_t)v->height;
    fx = dx / (fastf_t)v->width * 2.0;
    fy = -dy / (fastf_t)v->height / aspect * 2.0;
    VSET(tt, fx, fy, 0.0);
    MAT4X3PNT(work, v->view2model, tt);
    MAT_DELTAS_GET_NEG(vc, v->center);
    VSUB2(delta, work, vc);
    VSUB2(nvc, vc, delta);
    MAT_DELTAS_VEC_NEG(v->center, nvc);
    return bv_update(v);
}

static int
bv_adjust_scale_factor(struct bv *v, double factor)
{
    if (NEAR_ZERO(factor, SQRT_SMALL_FASTF))
	return 0;
    v->scale /= factor;
    if (v->scale < BV_MIN_SCALE)
	v->scale = BV_MIN_SCALE;
    v->size = 2.0 * v->scale;
    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
    return bv_update(v);
}

int
bv_mouse_delta_adjust(struct bv *v, int dx, int dy, const point_t keypoint,
	unsigned long long flags)
{
    struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;
    fastf_t mouse_dx = (fastf_t)dx;
    fastf_t mouse_dy = (fastf_t)dy;

    if (!bv_is_valid(v) || flags == BV_ADJUST_IDLE ||
	!(flags & (BV_ADJUST_ROT | BV_ADJUST_TRANS | BV_ADJUST_SCALE)) ||
	!bv_mouse_delta_settings_get(&settings, v) ||
	!bv_mouse_delta_clamp(&mouse_dx, &mouse_dy, v))
	return 0;

    if (flags & BV_ADJUST_ROT)
	return bv_adjust_rotation(v, mouse_dx, mouse_dy, keypoint,
	    settings.rotate_scale);
    if (flags & BV_ADJUST_TRANS)
	return bv_adjust_translation(v, mouse_dx, mouse_dy);

    const fastf_t delta = fabs(mouse_dx) > fabs(mouse_dy) ? mouse_dx :
	-mouse_dy;
    if (!v->height || NEAR_ZERO(delta, SQRT_SMALL_FASTF))
	return 0;
    const double factor = 1.0 + (double)delta *
	(double)settings.scale_scale / (double)v->height;
    if (!isfinite(factor) || factor <= SQRT_SMALL_FASTF)
	return 0;
    return bv_adjust_scale_factor(v, factor);
}

int
bv_mouse_state_set(struct bv *v, int x, int y)
{
    if (!bv_is_valid(v))
	return 0;

    v->previous_mouse_x = (fastf_t)v->mouse_x;
    v->previous_mouse_y = (fastf_t)v->mouse_y;
    v->mouse_x = x;
    v->mouse_y = y;
    return bv_screen_to_model(v->current_point, v, (fastf_t)x, (fastf_t)y);
}

int
bv_adjust(struct bv *v, int dx, int dy, const point_t keypoint, int UNUSED(mode), unsigned long long flags)
{
    if (!bv_is_valid(v) || flags == BV_ADJUST_IDLE)
	return 0;

    if (flags & BV_ADJUST_ROT)
	return bv_adjust_rotation(v, (fastf_t)dx, (fastf_t)dy, keypoint, 0.25);

    if (flags & BV_ADJUST_TRANS)
	return bv_adjust_translation(v, (fastf_t)dx, (fastf_t)dy);

    if (flags & BV_ADJUST_SCALE) {
	if (!dx || !dy)
	    return 0;
	return bv_adjust_scale_factor(v, (double)dy / (double)dx);
    }

    if (flags & BV_ADJUST_CENTER) {
	point_t center;
	if (!bv_screen_to_model(center, v, (fastf_t)dx, (fastf_t)dy))
	    return 0;
	MAT_DELTAS_VEC_NEG(v->center, center);
	return bv_update(v);
    }

    return 0;
}

int
bv_knobs_reset(struct bv_knobs *knobs, int category)
{
    if (!knobs)
	return 0;

    if (!category || category == BV_KNOBS_RATE) {
	knobs->rot_model_active = 0;
	VSETALL(knobs->rot_model, 0.0);
	knobs->rot_object_active = 0;
	VSETALL(knobs->rot_object, 0.0);
	knobs->rot_view_active = 0;
	VSETALL(knobs->rot_view, 0.0);
	knobs->trans_model_active = 0;
	VSETALL(knobs->trans_model, 0.0);
	knobs->trans_view_active = 0;
	VSETALL(knobs->trans_view, 0.0);
	knobs->scale_rate_active = 0;
	knobs->scale_rate = 0.0;
    }

    if (!category || category == BV_KNOBS_ABS) {
	VSETALL(knobs->abs_rot_model, 0.0);
	VSETALL(knobs->abs_rot_model_last, 0.0);
	VSETALL(knobs->abs_rot_object, 0.0);
	VSETALL(knobs->abs_rot_object_last, 0.0);
	VSETALL(knobs->abs_rot_view, 0.0);
	VSETALL(knobs->abs_rot_view_last, 0.0);
	VSETALL(knobs->abs_trans_model, 0.0);
	VSETALL(knobs->abs_trans_model_last, 0.0);
	VSETALL(knobs->abs_trans_view, 0.0);
	VSETALL(knobs->abs_trans_view_last, 0.0);
	knobs->abs_scale = 0.0;
    }

    return 1;
}

unsigned long long
bv_knobs_hash(const struct bv_knobs *knobs, struct bu_data_hash_state *state)
{
    int own_state = 0;
    unsigned long long hv;

    if (!knobs)
	return 0ULL;

    if (!state) {
	state = bu_data_hash_create();
	if (!state)
	    return 0ULL;
	own_state = 1;
    }

    bu_data_hash_update(state, &knobs->rot_model, sizeof(knobs->rot_model));
    bu_data_hash_update(state, &knobs->rot_model_active, sizeof(knobs->rot_model_active));
    bu_data_hash_update(state, &knobs->rot_model_origin, sizeof(knobs->rot_model_origin));
    bu_data_hash_update(state, &knobs->rot_object, sizeof(knobs->rot_object));
    bu_data_hash_update(state, &knobs->rot_object_active, sizeof(knobs->rot_object_active));
    bu_data_hash_update(state, &knobs->rot_object_origin, sizeof(knobs->rot_object_origin));
    bu_data_hash_update(state, &knobs->rot_view, sizeof(knobs->rot_view));
    bu_data_hash_update(state, &knobs->rot_view_active, sizeof(knobs->rot_view_active));
    bu_data_hash_update(state, &knobs->rot_view_origin, sizeof(knobs->rot_view_origin));
    bu_data_hash_update(state, &knobs->scale_rate, sizeof(knobs->scale_rate));
    bu_data_hash_update(state, &knobs->scale_rate_active, sizeof(knobs->scale_rate_active));
    bu_data_hash_update(state, &knobs->trans_model, sizeof(knobs->trans_model));
    bu_data_hash_update(state, &knobs->trans_model_active, sizeof(knobs->trans_model_active));
    bu_data_hash_update(state, &knobs->trans_view, sizeof(knobs->trans_view));
    bu_data_hash_update(state, &knobs->trans_view_active, sizeof(knobs->trans_view_active));
    bu_data_hash_update(state, &knobs->abs_rot_model, sizeof(knobs->abs_rot_model));
    bu_data_hash_update(state, &knobs->abs_rot_model_last, sizeof(knobs->abs_rot_model_last));
    bu_data_hash_update(state, &knobs->abs_rot_object, sizeof(knobs->abs_rot_object));
    bu_data_hash_update(state, &knobs->abs_rot_object_last, sizeof(knobs->abs_rot_object_last));
    bu_data_hash_update(state, &knobs->abs_rot_view, sizeof(knobs->abs_rot_view));
    bu_data_hash_update(state, &knobs->abs_rot_view_last, sizeof(knobs->abs_rot_view_last));
    bu_data_hash_update(state, &knobs->abs_scale, sizeof(knobs->abs_scale));
    bu_data_hash_update(state, &knobs->abs_trans_model, sizeof(knobs->abs_trans_model));
    bu_data_hash_update(state, &knobs->abs_trans_model_last, sizeof(knobs->abs_trans_model_last));
    bu_data_hash_update(state, &knobs->abs_trans_view, sizeof(knobs->abs_trans_view));
    bu_data_hash_update(state, &knobs->abs_trans_view_last, sizeof(knobs->abs_trans_view_last));

    if (!own_state)
	return 0ULL;

    hv = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return hv;
}

int
bv_knobs_cmd_process(vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, struct bv *v, const char *cmd, fastf_t factor, char origin, int model_flag, int incr_flag)
{
    struct bv_knobs *k;
    int ind = -1;
    char c;

    if (!bv_is_valid(v) || !cmd || !rvec || !do_rot || !tvec || !do_tran)
	return BRLCAD_ERROR;

    k = &v->knobs;
    c = (cmd[1] == '\0') ? cmd[0] : cmd[1];
    switch (c) {
	case 'x':
	case 'X':
	    ind = X;
	    break;
	case 'y':
	case 'Y':
	    ind = Y;
	    break;
	case 'z':
	case 'Z':
	    ind = Z;
	    break;
	case 'S':
	    ind = 0;
	    break;
	default:
	    return BRLCAD_ERROR;
    }

    if (cmd[1] == '\0') {
	if (cmd[0] == 'x' || cmd[0] == 'y' || cmd[0] == 'z') {
	    fastf_t *coord = model_flag ? &k->rot_model[ind] : &k->rot_view[ind];
	    if (incr_flag)
		*coord += factor;
	    else
		*coord = factor;
	    if (model_flag)
		k->rot_model_origin = origin;
	    else
		k->rot_view_origin = origin;
	    return BRLCAD_OK;
	}
	if (cmd[0] == 'X' || cmd[0] == 'Y' || cmd[0] == 'Z') {
	    fastf_t *coord = model_flag ? &k->trans_model[ind] : &k->trans_view[ind];
	    if (incr_flag)
		*coord += factor;
	    else
		*coord = factor;
	    return BRLCAD_OK;
	}
	if (cmd[0] == 'S') {
	    if (incr_flag)
		k->scale_rate += factor;
	    else
		k->scale_rate = factor;
	    return BRLCAD_OK;
	}
	return BRLCAD_ERROR;
    }

    if (cmd[0] == 'a' && strlen(cmd) == 2) {
	if (cmd[1] == 'x' || cmd[1] == 'y' || cmd[1] == 'z') {
	    fastf_t *abs_rot = model_flag ? k->abs_rot_model : k->abs_rot_view;
	    fastf_t *last = model_flag ? k->abs_rot_model_last : k->abs_rot_view_last;
	    if (incr_flag) {
		abs_rot[ind] += factor;
		(*rvec)[ind] = factor;
	    } else {
		(*rvec)[ind] = factor - last[ind];
		abs_rot[ind] = factor;
	    }
	    if (abs_rot[ind] < -180.0)
		abs_rot[ind] += 360.0;
	    else if (abs_rot[ind] > 180.0)
		abs_rot[ind] -= 360.0;
	    last[ind] = abs_rot[ind];
	    *do_rot = 1;
	    return BRLCAD_OK;
	}
	if (cmd[1] == 'X' || cmd[1] == 'Y' || cmd[1] == 'Z') {
	    fastf_t *abs_tran = model_flag ? k->abs_trans_model : k->abs_trans_view;
	    fastf_t *last = model_flag ? k->abs_trans_model_last : k->abs_trans_view_last;
	    fastf_t scale = (v->scale > SMALL_FASTF) ? v->scale : 1.0;
	    fastf_t sf = factor * v->local2base / scale;
	    if (incr_flag) {
		abs_tran[ind] += sf;
		last[ind] = abs_tran[ind];
		(*tvec)[ind] = factor;
	    } else {
		(*tvec)[ind] = factor - last[ind] * v->scale * v->base2local;
		abs_tran[ind] = sf;
		last[ind] = abs_tran[ind];
	    }
	    *do_tran = 1;
	    return BRLCAD_OK;
	}
	if (cmd[1] == 'S') {
	    if (incr_flag)
		v->absolute_scale += factor;
	    else
		v->absolute_scale = factor;
	    if (-SMALL_FASTF < v->absolute_scale && v->absolute_scale < SMALL_FASTF)
		v->scale = v->initial_scale;
	    else if (v->absolute_scale > 0.0)
		v->scale = v->initial_scale * (1.0 - v->absolute_scale);
	    else
		v->scale = v->initial_scale * (1.0 + (v->absolute_scale * -9.0));
	    if (v->scale < BV_MIN_SCALE)
		v->scale = BV_MIN_SCALE;
	    v->size = 2.0 * v->scale;
	    v->inverse_size = (v->size > 0.0) ? 1.0 / v->size : 0.0;
	    return bv_update(v) ? BRLCAD_OK : BRLCAD_ERROR;
	}
    }

    return BRLCAD_ERROR;
}

int
bv_knobs_translate(struct bv *v, const vect_t tvec, int model_flag)
{
    point_t vc;
    point_t nvc;

    if (!bv_is_valid(v) || !tvec)
	return 0;

    if (model_flag || v->coord_mode == 'o') {
	point_t delta;
	VSCALE(delta, tvec, v->local2base);
	MAT_DELTAS_GET_NEG(vc, v->center);
	VSUB2(nvc, vc, delta);
    } else {
	vect_t tt;
	point_t delta;
	point_t work;
	fastf_t scale = (v->scale > SMALL_FASTF) ? v->scale : 1.0;
	VSCALE(tt, tvec, v->local2base / scale);
	MAT4X3PNT(work, v->view2model, tt);
	MAT_DELTAS_GET_NEG(vc, v->center);
	VSUB2(delta, work, vc);
	VSUB2(nvc, vc, delta);
    }
    MAT_DELTAS_VEC_NEG(v->center, nvc);
    return bv_update(v);
}

int
bv_knobs_rotate(struct bv *v, const vect_t rvec, char origin, char coords, const matp_t obj_rot, const pointp_t pvt_pt)
{
    mat_t base;
    mat_t view_rot;
    point_t pivot_view;
    int recenter = 0;

    if (!bv_is_valid(v) || !rvec)
	return 0;

    MAT_IDN(base);
    bn_mat_angles(base, rvec[X], rvec[Y], rvec[Z]);
    MAT_IDN(view_rot);

    if (coords == 'm') {
	mat_t t1;
	mat_t rvinv;
	bn_mat_inv(rvinv, v->rotation);
	bn_mat_mul(t1, v->rotation, base);
	bn_mat_mul(view_rot, t1, rvinv);
    } else if (coords == 'o' && obj_rot) {
	mat_t obj_inv;
	mat_t model_rot;
	mat_t t1;
	mat_t t2;
	mat_t rvinv;
	bn_mat_inv(obj_inv, obj_rot);
	bn_mat_mul(t2, obj_rot, base);
	bn_mat_mul(model_rot, t2, obj_inv);
	bn_mat_inv(rvinv, v->rotation);
	bn_mat_mul(t1, v->rotation, model_rot);
	bn_mat_mul(view_rot, t1, rvinv);
    } else {
	MAT_COPY(view_rot, base);
    }

    switch (origin) {
	case 'e':
	    VSET(pivot_view, 0.0, 0.0, 1.0);
	    recenter = 1;
	    break;
	case 'm': {
	    point_t mzero = VINIT_ZERO;
	    MAT4X3PNT(pivot_view, v->model2view, mzero);
	    recenter = 1;
	    break;
	}
	case 'k': {
	    point_t mp;
	    if (pvt_pt)
		VMOVE(mp, pvt_pt);
	    else
		VSETALL(mp, 0.0);
	    MAT4X3PNT(pivot_view, v->model2view, mp);
	    recenter = 1;
	    break;
	}
	case 'v':
	default:
	    VSETALL(pivot_view, 0.0);
	    break;
    }

    if (recenter) {
	mat_t about;
	mat_t about_inv;
	point_t new_origin = VINIT_ZERO;
	point_t new_cent_view;
	point_t new_cent_model;
	bn_mat_xform_about_pnt(about, view_rot, pivot_view);
	bn_mat_inv(about_inv, about);
	MAT4X3PNT(new_cent_view, about_inv, new_origin);
	MAT4X3PNT(new_cent_model, v->view2model, new_cent_view);
	MAT_DELTAS_VEC_NEG(v->center, new_cent_model);
    }

    bn_mat_mul2(view_rot, v->rotation);
    return bv_update(v);
}

int
bv_knobs_update_rate_flags(struct bv *v)
{
    struct bv_knobs *k;

    if (!bv_is_valid(v))
	return 0;

    k = &v->knobs;
    k->rot_model_active = (!ZERO(k->rot_model[X]) || !ZERO(k->rot_model[Y]) || !ZERO(k->rot_model[Z]));
    k->trans_model_active = (!ZERO(k->trans_model[X]) || !ZERO(k->trans_model[Y]) || !ZERO(k->trans_model[Z]));
    k->rot_view_active = (!ZERO(k->rot_view[X]) || !ZERO(k->rot_view[Y]) || !ZERO(k->rot_view[Z]));
    k->trans_view_active = (!ZERO(k->trans_view[X]) || !ZERO(k->trans_view[Y]) || !ZERO(k->trans_view[Z]));
    k->scale_rate_active = !ZERO(k->scale_rate);
    return 1;
}

int
bv_knob_values_get(struct bv_knob_values *values, const struct bv *v)
{
    if (!values)
	return 0;

    memset(values, 0, sizeof(*values));
    if (!bv_is_valid(v))
	return 0;

    VMOVE(values->rate_rotation, v->knobs.rot_view);
    VMOVE(values->rate_translation, v->knobs.trans_view);
    values->rate_scale = v->knobs.scale_rate;
    VMOVE(values->absolute_rotation, v->knobs.abs_rot_view);
    VMOVE(values->absolute_translation, v->knobs.abs_trans_view);
    values->absolute_scale = v->absolute_scale;
    return 1;
}

int
bv_knobs_calibrate(struct bv *v)
{
    if (!bv_is_valid(v))
	return 0;

    VSETALL(v->knobs.abs_trans_view, 0.0);
    VSETALL(v->knobs.abs_trans_view_last, 0.0);
    VSETALL(v->knobs.abs_trans_model, 0.0);
    VSETALL(v->knobs.abs_trans_model_last, 0.0);
    return 1;
}

struct bv_context *
bv_context_create(void)
{
    struct bv_context *ctx = NULL;

    BU_ALLOC(ctx, struct bv_context);
    bv_context_init(ctx, NULL);
    return ctx;
}

void
bv_context_destroy(struct bv_context *ctx)
{
    if (!ctx)
	return;

    bv_context_free(ctx);
    bu_free(ctx, "bv_context");
}

void
bv_context_init(struct bv_context *ctx, struct bv *view)
{
    if (!ctx)
	return;

    memset(ctx, 0, sizeof(*ctx));
    ctx->magic = BV_CONTEXT_MAGIC;
    if (view) {
	ctx->view_ptr = view;
	ctx->owns_view = 0;
    } else {
	bv_init(&ctx->view);
	ctx->view_ptr = &ctx->view;
	ctx->owns_view = 1;
    }
    BU_PTBL_INIT(&ctx->callbacks);
    BU_PTBL_INIT(&ctx->owner_data);
}

void
bv_context_free(struct bv_context *ctx)
{
    if (!ctx || ctx->magic != BV_CONTEXT_MAGIC)
	return;

    if (ctx->view_set)
	bv_context_set_remove(ctx->view_set, ctx);
    bv_context_callbacks_clear(ctx);
    bu_ptbl_free(&ctx->callbacks);
    for (size_t i = 0; i < BU_PTBL_LEN(&ctx->owner_data); i++) {
	struct bv_context_owner_data_record *record =
	    (struct bv_context_owner_data_record *)BU_PTBL_GET(
		&ctx->owner_data, i);
	if (record)
	    BU_PUT(record, struct bv_context_owner_data_record);
    }
    bu_ptbl_free(&ctx->owner_data);
    if (ctx->owns_view && ctx->view_ptr == &ctx->view)
	bv_free(&ctx->view);
    ctx->view_ptr = NULL;
    ctx->owns_view = 0;
    ctx->view_set = NULL;
    ctx->magic = 0;
}

int
bv_context_is_valid(const struct bv_context *ctx)
{
    return (ctx && ctx->magic == BV_CONTEXT_MAGIC &&
	    bv_is_valid(ctx->view_ptr));
}

struct bv *
bv_context_view(struct bv_context *ctx)
{
    return bv_context_is_valid(ctx) ? ctx->view_ptr : NULL;
}

const struct bv *
bv_context_view_const(const struct bv_context *ctx)
{
    return bv_context_is_valid(ctx) ? ctx->view_ptr : NULL;
}

int
bv_context_name_set(struct bv_context *ctx, const char *name)
{
    return bv_name_set(bv_context_view(ctx), name);
}

const char *
bv_context_name_get(const struct bv_context *ctx)
{
    return bv_name_get(bv_context_view_const(ctx));
}

int
bv_context_user_data_set(struct bv_context *ctx, void *user_data)
{
    return bv_user_data_set(bv_context_view(ctx), user_data);
}

void *
bv_context_user_data_get(const struct bv_context *ctx)
{
    return bv_user_data_get(bv_context_view_const(ctx));
}

int
bv_context_owner_data_set(struct bv_context *ctx, const void *owner_key,
	void *data)
{
    if (!bv_context_is_valid(ctx) || !owner_key)
	return 0;

    for (size_t i = 0; i < BU_PTBL_LEN(&ctx->owner_data); i++) {
	struct bv_context_owner_data_record *record =
	    (struct bv_context_owner_data_record *)BU_PTBL_GET(
		&ctx->owner_data, i);
	if (!record || record->key != owner_key)
	    continue;
	if (data) {
	    record->data = data;
	    return 1;
	}
	bu_ptbl_rm(&ctx->owner_data, (long *)record);
	BU_PUT(record, struct bv_context_owner_data_record);
	return 1;
    }

    if (!data)
	return 0;

    struct bv_context_owner_data_record *record;
    BU_GET(record, struct bv_context_owner_data_record);
    record->key = owner_key;
    record->data = data;
    bu_ptbl_ins(&ctx->owner_data, (long *)record);
    return 1;
}

void *
bv_context_owner_data_get(const struct bv_context *ctx,
	const void *owner_key)
{
    if (!bv_context_is_valid(ctx) || !owner_key)
	return NULL;
    for (size_t i = 0; i < BU_PTBL_LEN(&ctx->owner_data); i++) {
	const struct bv_context_owner_data_record *record =
	    (const struct bv_context_owner_data_record *)BU_PTBL_GET(
		&ctx->owner_data, i);
	if (record && record->key == owner_key)
	    return record->data;
    }
    return NULL;
}

int
bv_context_dimensions_set(struct bv_context *ctx, int width, int height)
{
    return bv_dimensions_set(bv_context_view(ctx), width, height);
}

int
bv_context_width_get(const struct bv_context *ctx)
{
    return bv_width_get(bv_context_view_const(ctx));
}

int
bv_context_height_get(const struct bv_context *ctx)
{
    return bv_height_get(bv_context_view_const(ctx));
}

int
bv_context_refresh_request(struct bv_context *ctx, uint32_t flags)
{
    return bv_refresh_request(bv_context_view(ctx), flags);
}

int
bv_context_refresh_complete(struct bv_context *ctx)
{
    return bv_refresh_complete(bv_context_view(ctx));
}

int
bv_context_callback_add(struct bv_context *ctx, bv_context_callback_t callback,
	void *client_data)
{
    struct bv_context_callback_record *record = NULL;

    if (!bv_context_is_valid(ctx) || !callback)
	return 0;

    for (size_t i = 0; i < BU_PTBL_LEN(&ctx->callbacks); i++) {
	record = (struct bv_context_callback_record *)BU_PTBL_GET(&ctx->callbacks, i);
	if (record && record->callback == callback &&
		record->client_data == client_data)
	    return 1;
    }

    BU_ALLOC(record, struct bv_context_callback_record);
    record->callback = callback;
    record->client_data = client_data;
    bu_ptbl_ins(&ctx->callbacks, (long *)record);
    return 1;
}

int
bv_context_callback_remove(struct bv_context *ctx,
	bv_context_callback_t callback,
	void *client_data)
{
    int removed = 0;

    if (!ctx || ctx->magic != BV_CONTEXT_MAGIC || !callback)
	return 0;

    for (size_t i = BU_PTBL_LEN(&ctx->callbacks); i > 0; i--) {
	struct bv_context_callback_record *record =
	    (struct bv_context_callback_record *)BU_PTBL_GET(&ctx->callbacks,
		    i - 1);
	if (!record || record->callback != callback ||
		record->client_data != client_data)
	    continue;

	bu_ptbl_rm(&ctx->callbacks, (long *)record);
	BU_PUT(record, struct bv_context_callback_record);
	removed = 1;
    }

    return removed;
}

void
bv_context_callbacks_clear(struct bv_context *ctx)
{
    if (!ctx || ctx->magic != BV_CONTEXT_MAGIC)
	return;

    for (size_t i = 0; i < BU_PTBL_LEN(&ctx->callbacks); i++) {
	struct bv_context_callback_record *record =
	    (struct bv_context_callback_record *)BU_PTBL_GET(&ctx->callbacks, i);
	if (record)
	    BU_PUT(record, struct bv_context_callback_record);
    }
    bu_ptbl_reset(&ctx->callbacks);
}

int
bv_context_notify(struct bv_context *ctx, uint64_t changed_flags)
{
    struct bv_context_callback_record *snapshot = NULL;
    size_t count;

    if (!bv_context_is_valid(ctx))
	return 0;

    count = BU_PTBL_LEN(&ctx->callbacks);
    if (!count)
	return 1;

    snapshot = (struct bv_context_callback_record *)bu_calloc(count,
	    sizeof(struct bv_context_callback_record),
	    "bv_context callback snapshot");
    for (size_t i = 0; i < count; i++) {
	struct bv_context_callback_record *record =
	    (struct bv_context_callback_record *)BU_PTBL_GET(&ctx->callbacks, i);
	if (record)
	    snapshot[i] = *record;
    }

    for (size_t i = 0; i < count; i++) {
	if (snapshot[i].callback)
	    snapshot[i].callback(ctx, changed_flags, snapshot[i].client_data);
    }

    bu_free(snapshot, "bv_context callback snapshot");
    return 1;
}

int
bv_context_update(struct bv_context *ctx, uint64_t changed_flags)
{
    if (!bv_context_is_valid(ctx))
	return 0;

    if (!bv_update(ctx->view_ptr))
	return 0;

    return bv_context_notify(ctx, changed_flags ? changed_flags :
	    BV_CONTEXT_CHANGED_VIEW);
}

int
bv_context_settings_shared(const struct bv_context *a,
	const struct bv_context *b)
{
    return (bv_context_is_valid(a) && bv_context_is_valid(b) &&
	    a->view_set && a->view_set == b->view_set) ? 1 : 0;
}

struct bv_context_set *
bv_context_set_create(void)
{
    struct bv_context_set *s = NULL;
    BU_ALLOC(s, struct bv_context_set);
    bv_context_set_init(s);
    return s;
}

void
bv_context_set_destroy(struct bv_context_set *s)
{
    if (!s)
	return;

    bv_context_set_free(s);
    bu_free(s, "bv_context_set");
}

void
bv_context_set_init(struct bv_context_set *s)
{
    if (!s)
	return;

    memset(s, 0, sizeof(*s));
    s->magic = BV_CONTEXT_SET_MAGIC;
    bv_set_init(&s->bvset);
    BU_PTBL_INIT(&s->views);
    BU_PTBL_INIT(&s->recycle_pool);
}

void
bv_context_set_free(struct bv_context_set *s)
{
    if (!s || s->magic != BV_CONTEXT_SET_MAGIC)
	return;

    for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	struct bv_context *ctx =
	    (struct bv_context *)BU_PTBL_GET(&s->views, i);
	if (bv_context_is_valid(ctx) && ctx->view_set == s)
	    ctx->view_set = NULL;
    }
    bv_set_free(&s->bvset);
    bu_ptbl_free(&s->views);
    bu_ptbl_free(&s->recycle_pool);
    s->magic = 0;
}

struct bu_ptbl *
bv_context_set_views(struct bv_context_set *s)
{
    return (s && s->magic == BV_CONTEXT_SET_MAGIC) ? &s->views : NULL;
}

void *
bv_context_set_recycle_pool(struct bv_context_set *s)
{
    return (s && s->magic == BV_CONTEXT_SET_MAGIC) ? &s->recycle_pool : NULL;
}

struct bv_context *
bv_context_set_find(struct bv_context_set *s, const char *name)
{
    if (!s || s->magic != BV_CONTEXT_SET_MAGIC || !name)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	struct bv_context *ctx =
	    (struct bv_context *)BU_PTBL_GET(&s->views, i);
	const char *vname = bv_name_get(bv_context_view_const(ctx));
	if (vname && BU_STR_EQUAL(vname, name))
	    return ctx;
    }

    return NULL;
}

int
bv_context_set_attach(struct bv_context_set *s, struct bv_context *ctx)
{
    if (!bv_context_is_valid(ctx))
	return 0;
    if (s && s->magic != BV_CONTEXT_SET_MAGIC)
	return 0;

    if (ctx->view_set && ctx->view_set != s)
	bv_context_set_remove(ctx->view_set, ctx);
    ctx->view_set = s;
    return 1;
}

int
bv_context_set_add(struct bv_context_set *s, struct bv_context *ctx)
{
    struct bv *view = bv_context_view(ctx);
    if (!s || s->magic != BV_CONTEXT_SET_MAGIC || !view)
	return 0;

    if (ctx->view_set && ctx->view_set != s)
	bv_context_set_remove(ctx->view_set, ctx);

    bu_ptbl_ins_unique(&s->views, (long *)ctx);
    (void)bv_set_add(&s->bvset, view);
    ctx->view_set = s;
    return 1;
}

int
bv_context_set_remove(struct bv_context_set *s, struct bv_context *ctx)
{
    if (!s || s->magic != BV_CONTEXT_SET_MAGIC)
	return 0;

    if (!ctx) {
	for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	    struct bv_context *sctx =
		(struct bv_context *)BU_PTBL_GET(&s->views, i);
	    if (bv_context_is_valid(sctx) && sctx->view_set == s)
		sctx->view_set = NULL;
	}
	bv_set_free(&s->bvset);
	bv_set_init(&s->bvset);
	bu_ptbl_reset(&s->views);
	return 1;
    }

    if (!bv_context_is_valid(ctx))
	return 0;

    if (ctx->view_set == s)
	ctx->view_set = NULL;
    (void)bv_set_remove(&s->bvset, bv_context_view(ctx));
    (void)bu_ptbl_rm(&s->views, (long *)ctx);
    return 1;
}

struct bv_set *
bv_set_create(void)
{
    struct bv_set *s = NULL;
    BU_ALLOC(s, struct bv_set);
    bv_set_init(s);
    return s;
}

void
bv_set_destroy(struct bv_set *s)
{
    if (!s)
	return;

    bv_set_free(s);
    bu_free(s, "bv_set");
}

void
bv_set_init(struct bv_set *s)
{
    if (!s)
	return;

    memset(s, 0, sizeof(*s));
    s->magic = BV_SET_MAGIC;
    BU_PTBL_INIT(&s->views);
}

void
bv_set_free(struct bv_set *s)
{
    if (!s || s->magic != BV_SET_MAGIC)
	return;

    bu_ptbl_free(&s->views);
    s->magic = 0;
}

int
bv_set_add(struct bv_set *s, struct bv *v)
{
    if (!s || s->magic != BV_SET_MAGIC || !bv_is_valid(v))
	return 0;

    if (v->view_set && v->view_set != s)
	bv_set_remove(v->view_set, v);

    bu_ptbl_ins_unique(&s->views, (long *)v);
    v->view_set = s;
    return 1;
}

int
bv_set_remove(struct bv_set *s, struct bv *v)
{
    if (!s || s->magic != BV_SET_MAGIC || !v)
	return 0;

    if (v->view_set == s)
	v->view_set = NULL;
    (void)bu_ptbl_rm(&s->views, (long *)v);
    return 1;
}

struct bu_ptbl *
bv_set_views(struct bv_set *s)
{
    return (s && s->magic == BV_SET_MAGIC) ? &s->views : NULL;
}

struct bv *
bv_set_find(struct bv_set *s, const char *name)
{
    if (!s || s->magic != BV_SET_MAGIC || !name)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	struct bv *v = (struct bv *)BU_PTBL_GET(&s->views, i);
	const char *vname = bv_name_get(v);
	if (vname && BU_STR_EQUAL(vname, name))
	    return v;
    }

    return NULL;
}

/*
 * Local Variables:
 * mode: C
 * tab-width: 8
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
