/*              V I E W _ C O N T E X T _ N A T I V E . C
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
/** @file librt/view_context_native.c
 *
 * Native RT view context backed by libbv.  This file owns only RT-facing
 * wrapper state; generic view mechanics are delegated to bv_* routines.
 */

#include "common.h"

#include <string.h>

#include "bn/qmath.h"
#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "bu/str.h"
#include "bv.h"
#include "view_context_private.h"

#define RT_VIEW_CONTEXT_NATIVE_MAGIC 0x52545643u
#define RT_VIEW_SET_CONTEXT_NATIVE_MAGIC 0x52545653u

struct rt_view_context_native {
    uint32_t magic;
    struct bv view;
    struct rt_view_set_context_native *view_set;

    void *display_manager;
    void *tclcad_data;
    struct bu_ptbl *callbacks;
    rt_view_context_update_callback_t update_callback;
    void *update_callback_data;

    int edit_matrix_valid;
    mat_t edit_matrix;

    struct rt_view_lod_policy lod_policy;
    int lod_bounds_callback_set;
    int independent_scope_created;
    rt_view_scene_ref scene_root_ref;

    struct rt_view_interactive_rect_state interactive_rect;
    struct rt_view_adc_state adc;
    struct rt_view_grid_state grid;
    struct rt_view_axes_state model_axes;
    struct rt_view_axes_state view_axes;
    struct rt_view_other_state center_dot;
    struct rt_view_other_state scale_overlay;
    struct rt_view_params_state params;

    int snap_lines;
    int snap_source_flags;
    unsigned long long snap_kind_mask;
    double snap_tolerance_factor;
};

struct rt_view_set_context_native {
    uint32_t magic;
    struct bv_set bvset;
    struct bu_ptbl views;
    struct bu_ptbl recycle_pool;
};

static struct rt_view_context_native *
native_ctx(void *ctx)
{
    struct rt_view_context_native *n = (struct rt_view_context_native *)ctx;
    return (n && n->magic == RT_VIEW_CONTEXT_NATIVE_MAGIC) ? n : NULL;
}

static const struct rt_view_context_native *
native_ctx_const(const void *ctx)
{
    const struct rt_view_context_native *n =
	(const struct rt_view_context_native *)ctx;
    return (n && n->magic == RT_VIEW_CONTEXT_NATIVE_MAGIC) ? n : NULL;
}

static struct rt_view_set_context_native *
native_set(void *view_set_ctx)
{
    struct rt_view_set_context_native *s =
	(struct rt_view_set_context_native *)view_set_ctx;
    return (s && s->magic == RT_VIEW_SET_CONTEXT_NATIVE_MAGIC) ? s : NULL;
}

static const struct rt_view_set_context_native *
native_set_const(const void *view_set_ctx)
{
    const struct rt_view_set_context_native *s =
	(const struct rt_view_set_context_native *)view_set_ctx;
    return (s && s->magic == RT_VIEW_SET_CONTEXT_NATIVE_MAGIC) ? s : NULL;
}

static void
native_defaults(struct rt_view_context_native *n)
{
    struct rt_view_lod_policy lod = RT_VIEW_LOD_POLICY_INIT;
    struct rt_view_interactive_rect_state rect =
	RT_VIEW_INTERACTIVE_RECT_STATE_INIT;
    struct rt_view_adc_state adc = RT_VIEW_ADC_STATE_INIT;
    struct rt_view_grid_state grid = RT_VIEW_GRID_STATE_INIT;
    struct rt_view_axes_state axes = RT_VIEW_AXES_STATE_INIT;
    struct rt_view_other_state other = RT_VIEW_OTHER_STATE_INIT;
    struct rt_view_params_state params = RT_VIEW_PARAMS_STATE_INIT;

    if (!n)
	return;

    memset(n, 0, sizeof(*n));
    n->magic = RT_VIEW_CONTEXT_NATIVE_MAGIC;
    bv_init(&n->view);
    n->lod_policy = lod;
    n->interactive_rect = rect;
    n->adc = adc;
    n->grid = grid;
    n->model_axes = axes;
    n->view_axes = axes;
    n->center_dot = other;
    n->scale_overlay = other;
    n->params = params;
    n->snap_kind_mask = RT_VIEW_SNAP_KIND_GRID |
	RT_VIEW_SNAP_KIND_ENDPOINT | RT_VIEW_SNAP_KIND_MIDPOINT |
	RT_VIEW_SNAP_KIND_INTERSECTION | RT_VIEW_SNAP_KIND_PERPENDICULAR |
	RT_VIEW_SNAP_KIND_TANGENT | RT_VIEW_SNAP_KIND_OVERLAY_HANDLE;
    n->snap_tolerance_factor = 10.0;
}

static void
native_free_contents(struct rt_view_context_native *n)
{
    if (!n || n->magic != RT_VIEW_CONTEXT_NATIVE_MAGIC)
	return;

    if (n->callbacks) {
	bu_ptbl_free(n->callbacks);
	BU_PUT(n->callbacks, struct bu_ptbl);
	n->callbacks = NULL;
    }
    if (n->view_set) {
	bu_ptbl_rm(&n->view_set->views, (long *)n);
	(void)bv_set_remove(&n->view_set->bvset, &n->view);
	n->view_set = NULL;
    }
    bv_free(&n->view);
    n->magic = 0;
}

int
_rt_view_context_native_is(const void *ctx)
{
    return native_ctx_const(ctx) ? 1 : 0;
}

int
_rt_view_set_context_native_is(const void *view_set_ctx)
{
    return native_set_const(view_set_ctx) ? 1 : 0;
}

void *
_rt_view_context_native_create(void)
{
    struct rt_view_context_native *n = NULL;
    BU_ALLOC(n, struct rt_view_context_native);
    native_defaults(n);
    return n;
}

void *
_rt_view_context_native_create_with_set(void *view_set_ctx)
{
    struct rt_view_context_native *n =
	(struct rt_view_context_native *)_rt_view_context_native_create();
    _rt_view_context_native_view_set_attach(n, view_set_ctx);
    return n;
}

void *
_rt_view_context_native_create_copy_with_set(const void *src_ctx,
					     void *view_set_ctx)
{
    const struct rt_view_context_native *src = native_ctx_const(src_ctx);
    struct rt_view_context_native *n =
	(struct rt_view_context_native *)_rt_view_context_native_create();

    if (!n)
	return NULL;

    if (src) {
	struct bu_ptbl *callbacks = n->callbacks;
	uint32_t magic = n->magic;
	*n = *src;
	n->magic = magic;
	n->callbacks = callbacks;
	bv_init(&n->view);
	bv_copy(&n->view, &src->view);
	n->view_set = NULL;
    }
    _rt_view_context_native_view_set_attach(n, view_set_ctx);
    return n;
}

void
_rt_view_context_native_free(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return;

    native_free_contents(n);
    bu_free(n, "rt_view_context_native");
}

int
_rt_view_context_native_release_storage(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    native_free_contents(n);
    BU_PUT(n, struct rt_view_context_native);
    return 1;
}

void *
_rt_view_set_context_native_create(void)
{
    struct rt_view_set_context_native *s = NULL;
    BU_ALLOC(s, struct rt_view_set_context_native);
    _rt_view_set_context_native_init(s);
    return s;
}

void
_rt_view_set_context_native_destroy(void *view_set_ctx)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    if (!s)
	return;

    _rt_view_set_context_native_free(s);
    bu_free(s, "rt_view_set_context_native");
}

void
_rt_view_set_context_native_init(void *view_set_ctx)
{
    struct rt_view_set_context_native *s =
	(struct rt_view_set_context_native *)view_set_ctx;
    if (!s)
	return;

    memset(s, 0, sizeof(*s));
    s->magic = RT_VIEW_SET_CONTEXT_NATIVE_MAGIC;
    bv_set_init(&s->bvset);
    BU_PTBL_INIT(&s->views);
    BU_PTBL_INIT(&s->recycle_pool);
}

void
_rt_view_set_context_native_free(void *view_set_ctx)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    if (!s)
	return;

    for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	struct rt_view_context_native *n =
	    (struct rt_view_context_native *)BU_PTBL_GET(&s->views, i);
	if (native_ctx(n))
	    n->view_set = NULL;
    }
    bv_set_free(&s->bvset);
    bu_ptbl_free(&s->views);
    bu_ptbl_free(&s->recycle_pool);
    s->magic = 0;
}

struct bu_ptbl *
_rt_view_set_context_native_views(void *view_set_ctx)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    return s ? &s->views : NULL;
}

void *
_rt_view_set_context_native_recycle_pool(void *view_set_ctx)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    return s ? &s->recycle_pool : NULL;
}

void *
_rt_view_set_context_native_find_view(void *view_set_ctx, const char *name)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    if (!s || !name)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	struct rt_view_context_native *n =
	    (struct rt_view_context_native *)BU_PTBL_GET(&s->views, i);
	const char *vname = _rt_view_context_native_name_get(n);
	if (vname && BU_STR_EQUAL(vname, name))
	    return n;
    }

    return NULL;
}

int
_rt_view_set_context_native_add(void *view_set_ctx, void *view_ctx)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    struct rt_view_context_native *n = native_ctx(view_ctx);
    if (!s || !n)
	return 0;

    if (n->view_set && n->view_set != s)
	_rt_view_set_context_native_remove(n->view_set, n);

    bu_ptbl_ins_unique(&s->views, (long *)n);
    (void)bv_set_add(&s->bvset, &n->view);
    n->view_set = s;
    return 1;
}

int
_rt_view_set_context_native_remove(void *view_set_ctx, void *view_ctx)
{
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    struct rt_view_context_native *n = native_ctx(view_ctx);
    if (!s)
	return 0;

    if (!view_ctx) {
	for (size_t i = 0; i < BU_PTBL_LEN(&s->views); i++) {
	    struct rt_view_context_native *sn =
		(struct rt_view_context_native *)BU_PTBL_GET(&s->views, i);
	    if (native_ctx(sn) && sn->view_set == s)
		sn->view_set = NULL;
	}
	bv_set_free(&s->bvset);
	bv_set_init(&s->bvset);
	bu_ptbl_reset(&s->views);
	return 1;
    }

    if (!n)
	return 0;

    bu_ptbl_rm(&s->views, (long *)n);
    (void)bv_set_remove(&s->bvset, &n->view);
    if (n->view_set == s)
	n->view_set = NULL;
    return 1;
}

int
_rt_view_context_native_view_set_attach(void *ctx, void *view_set_ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    struct rt_view_set_context_native *s = native_set(view_set_ctx);
    if (!n || (view_set_ctx && !s))
	return 0;

    n->view_set = s;
    return 1;
}

void *
_rt_view_context_native_user_data_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_user_data_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_user_data_set(void *ctx, void *user_data)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_user_data_set(n ? &n->view : NULL, user_data);
}

int
_rt_view_context_native_tclcad_data_set(void *ctx, void *tcl_data)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->tclcad_data = tcl_data;
    return 1;
}

int
_rt_view_context_native_callbacks_set(void *ctx, struct bu_ptbl *callbacks)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    if (n->callbacks && n->callbacks != callbacks) {
	bu_ptbl_free(n->callbacks);
	BU_PUT(n->callbacks, struct bu_ptbl);
    }
    n->callbacks = callbacks;
    return 1;
}

int
_rt_view_context_native_update_callback_set(void *ctx,
					   rt_view_context_update_callback_t callback,
					   void *data)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->update_callback = callback;
    n->update_callback_data = data;
    return 1;
}

int
_rt_view_context_native_name_set(void *ctx, const char *name)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_name_set(&n->view, name) : 0;
}

const char *
_rt_view_context_native_name_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? bv_name_get(&n->view) : NULL;
}

int
_rt_view_context_native_width_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.width : 0;
}

int
_rt_view_context_native_height_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.height : 0;
}

int
_rt_view_context_native_dimensions_set(void *ctx, int width, int height)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_dimensions_set(&n->view, width, height) : 0;
}

void *
_rt_view_context_native_display_manager_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->display_manager : NULL;
}

int
_rt_view_context_native_display_manager_set(void *ctx, void *dmp)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->display_manager = dmp;
    return 1;
}

int
_rt_view_context_native_edit_matrix_set(void *ctx, matp_t edit_mat)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !edit_mat)
	return 0;

    MAT_COPY(n->edit_matrix, edit_mat);
    n->edit_matrix_valid = 1;
    return 1;
}

int
_rt_view_context_native_edit_matrix_clear(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    MAT_IDN(n->edit_matrix);
    n->edit_matrix_valid = 0;
    return 1;
}

uint64_t
_rt_view_context_native_frame_revision_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.frame_revision : 0;
}

uint64_t
_rt_view_context_native_frame_revision_bump(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.frame_revision++;
    return n->view.frame_revision;
}

int
_rt_view_context_native_model_matrices_identity(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    MAT_IDN(n->view.model2view);
    MAT_IDN(n->view.view2model);
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_frametime_set(void *ctx, uint64_t frametime)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_frametime_set(n ? &n->view : NULL, frametime);
}

int
_rt_view_context_native_refresh_request(void *ctx, uint32_t flags)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    return bv_refresh_request(&n->view, flags);
}

int
_rt_view_context_native_refresh_dirty_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_refresh_dirty_get(n ? &n->view : NULL);
}

uint32_t
_rt_view_context_native_refresh_consume(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_refresh_consume(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_complete(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_refresh_complete(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_enabled_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_refresh_enabled_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_enabled_set(void *ctx, int enabled)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_refresh_enabled_set(n ? &n->view : NULL, enabled);
}

int
_rt_view_context_native_refresh_suppressed_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_refresh_suppressed_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_suppress_begin(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_refresh_suppress_begin(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_suppress_end(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_refresh_suppress_end(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_drawn_count_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_refresh_drawn_count_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_refresh_drawn_count_set(void *ctx, int count)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_refresh_drawn_count_set(n ? &n->view : NULL, count);
}

fastf_t
_rt_view_context_native_scale_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.scale : 1.0;
}

int
_rt_view_context_native_scale_set(void *ctx, fastf_t scale)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_scale_set(&n->view, scale) : 0;
}

fastf_t *
_rt_view_context_native_scale_storage_get(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? &n->view.scale : NULL;
}

int
_rt_view_context_native_scale_state_set(void *ctx, fastf_t scale,
					fastf_t initial_scale,
					fastf_t absolute_scale,
					fastf_t size,
					fastf_t inverse_size)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.scale = scale;
    n->view.initial_scale = initial_scale;
    n->view.absolute_scale = absolute_scale;
    n->view.size = size;
    n->view.inverse_size = inverse_size;
    return bv_update(&n->view);
}

fastf_t
_rt_view_context_native_initial_scale_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.initial_scale : 1.0;
}

int
_rt_view_context_native_initial_scale_set(void *ctx, fastf_t scale)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.initial_scale = scale;
    n->view.frame_revision++;
    return 1;
}

fastf_t
_rt_view_context_native_absolute_scale_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.absolute_scale : 0.0;
}

int
_rt_view_context_native_absolute_scale_set(void *ctx, fastf_t scale)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.absolute_scale = scale;
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_unit_conversion_set(void *ctx, fastf_t local2base,
					    fastf_t base2local)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.local2base = local2base;
    n->view.base2local = base2local;
    n->view.frame_revision++;
    return 1;
}

fastf_t
_rt_view_context_native_local2base_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.local2base : 1.0;
}

fastf_t
_rt_view_context_native_base2local_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.base2local : 1.0;
}

fastf_t
_rt_view_context_native_size_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.size : 0.0;
}

int
_rt_view_context_native_size_set(void *ctx, fastf_t size)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_size_set(&n->view, size) : 0;
}

fastf_t
_rt_view_context_native_inverse_size_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.inverse_size : 1.0;
}

fastf_t
_rt_view_context_native_perspective_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.perspective : 0.0;
}

int
_rt_view_context_native_perspective_set(void *ctx, fastf_t perspective)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.perspective = perspective;
    n->view.frame_revision++;
    return 1;
}

fastf_t
_rt_view_context_native_radius_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.radius : 1.0;
}

void
_rt_view_context_native_info_get(struct rt_view_info *info, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!info)
	return;

    rt_view_info_init(info);
    if (!n)
	return;

    info->width = n->view.width;
    info->height = n->view.height;
    info->size = n->view.size;
    info->lod.scale = n->lod_policy.scale;
    info->lod.curve_scale = n->lod_policy.curve_scale;
    info->lod.point_scale = n->lod_policy.point_scale;
    info->lod.bot_threshold = n->lod_policy.bot_threshold;
    rt_view_info_sanitize(info);
}

int
_rt_view_context_native_obb_get(const void *ctx, point_t center,
				vect_t extent1, vect_t extent2,
				vect_t extent3)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!n || !center || !extent1 || !extent2 || !extent3)
	return 0;

    VMOVE(center, n->view.obb_center);
    VMOVE(extent1, n->view.obb_extent1);
    VMOVE(extent2, n->view.obb_extent2);
    VMOVE(extent3, n->view.obb_extent3);
    return 1;
}

int
_rt_view_context_native_model2view_get(mat_t model2view, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_model2view_get(model2view, n ? &n->view : NULL);
}

int
_rt_view_context_native_model2view_set(void *ctx, const mat_t model2view)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !model2view)
	return 0;

    MAT_COPY(n->view.model2view, model2view);
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_view2model_get(mat_t view2model, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_view2model_get(view2model, n ? &n->view : NULL);
}

int
_rt_view_context_native_view2model_set(void *ctx, const mat_t view2model)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !view2model)
	return 0;

    MAT_COPY(n->view.view2model, view2model);
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_pmodel2view_get(mat_t pmodel2view, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_pmodel2view_get(pmodel2view, n ? &n->view : NULL);
}

int
_rt_view_context_native_pmodel2view_set(void *ctx, const mat_t pmodel2view)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !pmodel2view)
	return 0;

    MAT_COPY(n->view.pmodel2view, pmodel2view);
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_pmat_get(mat_t pmat, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_pmat_get(pmat, n ? &n->view : NULL);
}

int
_rt_view_context_native_pmat_set(void *ctx, const mat_t pmat)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_pmat_set(&n->view, pmat) : 0;
}

int
_rt_view_context_native_rotation_get(mat_t rotation, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_rotation_get(rotation, n ? &n->view : NULL);
}

int
_rt_view_context_native_rotation_set(void *ctx, const mat_t rotation)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_rotation_set(&n->view, rotation) : 0;
}

int
_rt_view_context_native_center_get(mat_t center, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!center)
	return 0;
    if (!n) {
	MAT_IDN(center);
	return 0;
    }
    MAT_COPY(center, n->view.center);
    return 1;
}

int
_rt_view_context_native_center_mat_set(void *ctx, const mat_t center)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !center)
	return 0;

    MAT_COPY(n->view.center, center);
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_center_set(void *ctx, const point_t center)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_center_set(&n->view, center) : 0;
}

int
_rt_view_context_native_plane_get(plane_t *plane, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_plane_get(plane, n ? &n->view : NULL);
}

int
_rt_view_context_native_orientation_quat_get(quat_t orientation,
					     const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!orientation)
	return 0;
    if (!n) {
	mat_t identity;
	MAT_IDN(identity);
	quat_mat2quat(orientation, identity);
	return 0;
    }
    quat_mat2quat(orientation, n->view.rotation);
    return 1;
}

int
_rt_view_context_native_aet_get(vect_t aet, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_aet_get(aet, n ? &n->view : NULL);
}

int
_rt_view_context_native_aet_set(void *ctx, const vect_t aet)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_aet_set(&n->view, aet) : 0;
}

int
_rt_view_context_native_aet_state_set(void *ctx, const vect_t aet)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_aet_state_set(&n->view, aet) : 0;
}

int
_rt_view_context_native_update(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    int ret;

    if (!n)
	return 0;

    ret = bv_update(&n->view);
    if (ret && n->update_callback)
	n->update_callback(ctx, n->update_callback_data);
    return ret;
}

int
_rt_view_context_native_autoview_bounds(void *ctx, fastf_t scale,
					const point_t min,
					const point_t max)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_autoview_bounds(&n->view, scale, min, max) : 0;
}

int
_rt_view_context_native_adjust(void *ctx, int dx, int dy, point_t keypoint,
			       int mode, unsigned long long flags)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_adjust(&n->view, dx, dy, keypoint, mode, flags) : 0;
}

unsigned long long
_rt_view_context_native_hash(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    struct bu_data_hash_state *state;
    unsigned long long hv;

    if (!n)
	return 0ULL;

    state = bu_data_hash_create();
    if (!state)
	return 0ULL;

    bu_data_hash_update(state, &n->view.width, sizeof(n->view.width));
    bu_data_hash_update(state, &n->view.height, sizeof(n->view.height));
    bu_data_hash_update(state, &n->view.scale, sizeof(n->view.scale));
    bu_data_hash_update(state, &n->view.size, sizeof(n->view.size));
    bu_data_hash_update(state, &n->view.perspective,
	    sizeof(n->view.perspective));
    bu_data_hash_update(state, &n->view.aet, sizeof(n->view.aet));
    bu_data_hash_update(state, &n->view.rotation, sizeof(n->view.rotation));
    bu_data_hash_update(state, &n->view.center, sizeof(n->view.center));
    bu_data_hash_update(state, &n->view.model2view,
	    sizeof(n->view.model2view));
    bu_data_hash_update(state, &n->lod_policy, sizeof(n->lod_policy));
    hv = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return hv;
}

size_t
_rt_view_context_native_clear(void *ctx, int UNUSED(flags))
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    (void)bv_refresh_complete(&n->view);
    return 0;
}

rt_view_scene_ref
_rt_view_context_native_independent_scope_ref(void *ctx, int create)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return rt_view_scene_ref_null();
    if (create)
	n->independent_scope_created = 1;
    return rt_view_scene_ref_null();
}

int
_rt_view_context_native_is_independent(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return (n && n->independent_scope_created) ? 1 : 0;
}

int
_rt_view_context_native_independent_scope_is_null(void *ctx, int create)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 1;
    if (create)
	n->independent_scope_created = 1;
    return n->independent_scope_created ? 0 : 1;
}

void
_rt_view_context_native_independent_scope_destroy(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return;
    n->independent_scope_created = 0;
}

rt_view_scene_ref
_rt_view_context_native_scene_root_ref(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->scene_root_ref : rt_view_scene_ref_null();
}

int
_rt_view_context_native_scene_root_ref_attach(void *ctx,
					     rt_view_scene_ref root_ref)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;
    if (!rt_view_scene_ref_is_null(root_ref) &&
	    rt_view_scene_ref_backend(root_ref) != RT_VIEW_SCENE_BACKEND_OBOL)
	return 0;
    n->scene_root_ref = root_ref;
    return !rt_view_scene_ref_is_null(root_ref) ? 1 : 0;
}

int
_rt_view_context_native_scene_attached(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return (n && !rt_view_scene_ref_is_null(n->scene_root_ref)) ? 1 : 0;
}

int
_rt_view_context_native_measure_candidates(void *ctx, point_t a, point_t b,
					  struct rt_view_measure_result *out)
{
    vect_t delta;
    if (out)
	memset(out, 0, sizeof(*out));
    if (!native_ctx(ctx) || !a || !b || !out)
	return 0;

    VSUB2(delta, b, a);
    out->distance = MAGNITUDE(delta);
    out->projection = out->distance;
    out->normal_alignment = 0.0;
    out->valid = 1;
    return 1;
}

int
_rt_view_context_native_screen_to_view(fastf_t *fx, fastf_t *fy, void *ctx,
				       fastf_t x, fastf_t y)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_screen_to_view(fx, fy, &n->view, x, y) : 0;
}

int
_rt_view_context_native_screen_point_get(point_t point, void *ctx,
					 fastf_t x, fastf_t y)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_screen_to_model(point, &n->view, x, y) : 0;
}

int
_rt_view_context_native_current_point_get(point_t point, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!point)
	return 0;
    if (!n) {
	VSETALL(point, 0.0);
	return 0;
    }
    VMOVE(point, n->view.current_point);
    return 1;
}

int
_rt_view_context_native_current_point_set(void *ctx, const point_t point)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !point)
	return 0;

    VMOVE(n->view.current_point, point);
    n->view.frame_revision++;
    return 1;
}

int
_rt_view_context_native_previous_mouse_get(fastf_t *x, fastf_t *y,
					   const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (x)
	*x = 0.0;
    if (y)
	*y = 0.0;
    if (!x || !y || !n)
	return 0;

    *x = n->view.previous_mouse_x;
    *y = n->view.previous_mouse_y;
    return 1;
}

int
_rt_view_context_native_previous_mouse_set(void *ctx, fastf_t x, fastf_t y)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.previous_mouse_x = x;
    n->view.previous_mouse_y = y;
    return 1;
}

int
_rt_view_context_native_mouse_delta_settings_get(
	struct rt_view_mouse_delta_settings *settings, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    struct rt_view_mouse_delta_settings zero =
	RT_VIEW_MOUSE_DELTA_SETTINGS_INIT;

    if (!settings)
	return 0;

    *settings = zero;
    if (!n)
	return 0;

    settings->min_delta = n->view.min_mouse_delta;
    settings->max_delta = n->view.max_mouse_delta;
    settings->rotate_scale = n->view.rotate_scale;
    settings->scale_scale = n->view.scale_scale;
    return 1;
}

int
_rt_view_context_native_mouse_state_set(void *ctx, int x, int y)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_mouse_state_set(&n->view, x, y) : 0;
}

int
_rt_view_context_native_interactive_rect_state_get(
	struct rt_view_interactive_rect_state *record, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_interactive_rect_state zero =
	    RT_VIEW_INTERACTIVE_RECT_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->interactive_rect;
    return 1;
}

int
_rt_view_context_native_interactive_rect_state_set(
	void *ctx, const struct rt_view_interactive_rect_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->interactive_rect = *record;
    return 1;
}

int
_rt_view_context_native_adc_state_get(struct rt_view_adc_state *record,
				      const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_adc_state zero = RT_VIEW_ADC_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->adc;
    return 1;
}

int
_rt_view_context_native_adc_state_set(void *ctx,
				      const struct rt_view_adc_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->adc = *record;
    return 1;
}

int
_rt_view_context_native_grid_state_get(struct rt_view_grid_state *record,
				       const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_grid_state zero = RT_VIEW_GRID_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->grid;
    return 1;
}

int
_rt_view_context_native_grid_state_set(void *ctx,
				       const struct rt_view_grid_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->grid = *record;
    return 1;
}

int
_rt_view_context_native_model_axes_state_get(
	struct rt_view_axes_state *record, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_axes_state zero = RT_VIEW_AXES_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->model_axes;
    return 1;
}

int
_rt_view_context_native_model_axes_state_set(
	void *ctx, const struct rt_view_axes_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->model_axes = *record;
    return 1;
}

int
_rt_view_context_native_view_axes_state_get(
	struct rt_view_axes_state *record, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_axes_state zero = RT_VIEW_AXES_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->view_axes;
    return 1;
}

int
_rt_view_context_native_view_axes_state_set(
	void *ctx, const struct rt_view_axes_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->view_axes = *record;
    return 1;
}

int
_rt_view_context_native_center_dot_state_get(
	struct rt_view_other_state *record, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_other_state zero = RT_VIEW_OTHER_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->center_dot;
    return 1;
}

int
_rt_view_context_native_center_dot_state_set(
	void *ctx, const struct rt_view_other_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->center_dot = *record;
    return 1;
}

int
_rt_view_context_native_scale_overlay_state_get(
	struct rt_view_other_state *record, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_other_state zero = RT_VIEW_OTHER_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->scale_overlay;
    return 1;
}

int
_rt_view_context_native_scale_overlay_state_set(
	void *ctx, const struct rt_view_other_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->scale_overlay = *record;
    return 1;
}

int
_rt_view_context_native_params_state_get(struct rt_view_params_state *record,
					 const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!record)
	return 0;

    if (!n) {
	struct rt_view_params_state zero = RT_VIEW_PARAMS_STATE_INIT;
	*record = zero;
	return 0;
    }

    *record = n->params;
    return 1;
}

int
_rt_view_context_native_params_state_set(
	void *ctx, const struct rt_view_params_state *record)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !record)
	return 0;

    n->params = *record;
    return 1;
}

int
_rt_view_context_native_eye_pos_get(point_t eye_pos, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_eye_pos_get(eye_pos, n ? &n->view : NULL);
}

int
_rt_view_context_native_eye_pos_set(void *ctx, const point_t eye_pos)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_eye_pos_set(&n->view, eye_pos) : 0;
}

int
_rt_view_context_native_keypoint_get(point_t keypoint, const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_keypoint_get(keypoint, n ? &n->view : NULL);
}

int
_rt_view_context_native_keypoint_set(void *ctx, const point_t keypoint)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_keypoint_set(&n->view, keypoint) : 0;
}

char
_rt_view_context_native_rotate_about_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.rotate_about : 'v';
}

int
_rt_view_context_native_rotate_about_set(void *ctx, char rotate_about)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.rotate_about = rotate_about;
    return 1;
}

char
_rt_view_context_native_coord_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->view.coord_mode : 'v';
}

int
_rt_view_context_native_coord_set(void *ctx, char coord)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->view.coord_mode = coord;
    return 1;
}

int
_rt_view_context_native_zclip_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_zclip_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_zclip_set(void *ctx, int zclip)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_zclip_set(n ? &n->view : NULL, zclip);
}

int
_rt_view_context_native_framebuffer_mode_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_framebuffer_mode_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_framebuffer_mode_set(void *ctx, int mode)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_framebuffer_mode_set(n ? &n->view : NULL, mode);
}

int
_rt_view_context_native_cleared_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return bv_cleared_get(n ? &n->view : NULL);
}

int
_rt_view_context_native_cleared_set(void *ctx, int cleared)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return bv_cleared_set(n ? &n->view : NULL, cleared);
}

int
_rt_view_context_native_settings_shared(const void *a, const void *b)
{
    const struct rt_view_context_native *av = native_ctx_const(a);
    const struct rt_view_context_native *bv = native_ctx_const(b);
    return (av && bv && av->view_set && av->view_set == bv->view_set) ? 1 : 0;
}

int
_rt_view_context_native_snap_lines_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->snap_lines : 0;
}

int
_rt_view_context_native_snap_lines_set(void *ctx, int enabled)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->snap_lines = enabled ? 1 : 0;
    return 1;
}

int
_rt_view_context_native_snap_source_flags_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->snap_source_flags : 0;
}

int
_rt_view_context_native_snap_source_flags_set(void *ctx, int flags)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->snap_source_flags = flags;
    return 1;
}

unsigned long long
_rt_view_context_native_snap_kind_mask_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->snap_kind_mask : 0ULL;
}

int
_rt_view_context_native_snap_exclude_feature_clear(void *ctx)
{
    return native_ctx(ctx) ? 1 : 0;
}

unsigned long long
_rt_view_context_native_prepare_tcl_snap(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0ULL;

    n->snap_source_flags |= RT_VIEW_SNAP_TCL;
    return n->snap_kind_mask;
}

double
_rt_view_context_native_snap_tolerance_factor_get(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->snap_tolerance_factor : 1.0;
}

int
_rt_view_context_native_snap_tolerance_factor_set(void *ctx, double factor)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->snap_tolerance_factor = (factor > 0.0) ? factor : 1.0;
    return 1;
}

static int
native_knob_category(int category)
{
    switch (category) {
	case RT_VIEW_KNOBS_RATE:
	    return BV_KNOBS_RATE;
	case RT_VIEW_KNOBS_ABS:
	    return BV_KNOBS_ABS;
	case RT_VIEW_KNOBS_ALL:
	default:
	    return BV_KNOBS_ALL;
    }
}

static void
native_knobs_from_rt(struct bv_knobs *dst, const struct rt_view_knobs *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    VMOVE(dst->rot_model, src->rot_m);
    dst->rot_model_active = src->rot_m_flag;
    dst->rot_model_origin = src->origin_m;
    dst->rot_model_data = src->rot_m_udata;
    VMOVE(dst->rot_object, src->rot_o);
    dst->rot_object_active = src->rot_o_flag;
    dst->rot_object_origin = src->origin_o;
    dst->rot_object_data = src->rot_o_udata;
    VMOVE(dst->rot_view, src->rot_v);
    dst->rot_view_active = src->rot_v_flag;
    dst->rot_view_origin = src->origin_v;
    dst->rot_view_data = src->rot_v_udata;
    dst->scale_rate = src->sca;
    dst->scale_rate_active = src->sca_flag;
    dst->scale_rate_data = src->sca_udata;
    VMOVE(dst->trans_model, src->tra_m);
    dst->trans_model_active = src->tra_m_flag;
    dst->trans_model_data = src->tra_m_udata;
    VMOVE(dst->trans_view, src->tra_v);
    dst->trans_view_active = src->tra_v_flag;
    dst->trans_view_data = src->tra_v_udata;
    VMOVE(dst->abs_rot_model, src->rot_m_abs);
    VMOVE(dst->abs_rot_model_last, src->rot_m_abs_last);
    VMOVE(dst->abs_rot_object, src->rot_o_abs);
    VMOVE(dst->abs_rot_object_last, src->rot_o_abs_last);
    VMOVE(dst->abs_rot_view, src->rot_v_abs);
    VMOVE(dst->abs_rot_view_last, src->rot_v_abs_last);
    dst->abs_scale = src->sca_abs;
    VMOVE(dst->abs_trans_model, src->tra_m_abs);
    VMOVE(dst->abs_trans_model_last, src->tra_m_abs_last);
    VMOVE(dst->abs_trans_view, src->tra_v_abs);
    VMOVE(dst->abs_trans_view_last, src->tra_v_abs_last);
}

static void
native_knobs_to_rt(struct rt_view_knobs *dst, const struct bv_knobs *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    VMOVE(dst->rot_m, src->rot_model);
    dst->rot_m_flag = src->rot_model_active;
    dst->origin_m = src->rot_model_origin;
    dst->rot_m_udata = src->rot_model_data;
    VMOVE(dst->rot_o, src->rot_object);
    dst->rot_o_flag = src->rot_object_active;
    dst->origin_o = src->rot_object_origin;
    dst->rot_o_udata = src->rot_object_data;
    VMOVE(dst->rot_v, src->rot_view);
    dst->rot_v_flag = src->rot_view_active;
    dst->origin_v = src->rot_view_origin;
    dst->rot_v_udata = src->rot_view_data;
    dst->sca = src->scale_rate;
    dst->sca_flag = src->scale_rate_active;
    dst->sca_udata = src->scale_rate_data;
    VMOVE(dst->tra_m, src->trans_model);
    dst->tra_m_flag = src->trans_model_active;
    dst->tra_m_udata = src->trans_model_data;
    VMOVE(dst->tra_v, src->trans_view);
    dst->tra_v_flag = src->trans_view_active;
    dst->tra_v_udata = src->trans_view_data;
    VMOVE(dst->rot_m_abs, src->abs_rot_model);
    VMOVE(dst->rot_m_abs_last, src->abs_rot_model_last);
    VMOVE(dst->rot_o_abs, src->abs_rot_object);
    VMOVE(dst->rot_o_abs_last, src->abs_rot_object_last);
    VMOVE(dst->rot_v_abs, src->abs_rot_view);
    VMOVE(dst->rot_v_abs_last, src->abs_rot_view_last);
    dst->sca_abs = src->abs_scale;
    VMOVE(dst->tra_m_abs, src->abs_trans_model);
    VMOVE(dst->tra_m_abs_last, src->abs_trans_model_last);
    VMOVE(dst->tra_v_abs, src->abs_trans_view);
    VMOVE(dst->tra_v_abs_last, src->abs_trans_view_last);
}

int
_rt_view_context_native_knobs_reset(void *ctx, int category)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_reset(&n->view.knobs, native_knob_category(category)) : 0;
}

int
_rt_view_context_native_knobs_state_get(struct rt_view_knobs *knobs,
					const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!knobs)
	return 0;

    native_knobs_to_rt(knobs, n ? &n->view.knobs : NULL);
    return n ? 1 : 0;
}

int
_rt_view_context_native_knobs_state_set(void *ctx,
					const struct rt_view_knobs *knobs)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !knobs)
	return 0;

    native_knobs_from_rt(&n->view.knobs, knobs);
    return 1;
}

unsigned long long
_rt_view_context_native_knobs_hash(void *ctx, struct bu_data_hash_state *state)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_hash(&n->view.knobs, state) : 0ULL;
}

int
_rt_view_context_native_knobs_cmd_process(vect_t *rvec, int *do_rot,
					  vect_t *tvec, int *do_tran,
					  void *ctx, const char *cmd,
					  fastf_t factor, char origin,
					  int model_flag, int incr_flag)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_cmd_process(rvec, do_rot, tvec, do_tran, &n->view,
	    cmd, factor, origin, model_flag, incr_flag) : BRLCAD_ERROR;
}

int
_rt_view_context_native_knobs_translate(void *ctx, const vect_t tvec,
					int model_flag)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_translate(&n->view, tvec, model_flag) : 0;
}

int
_rt_view_context_native_knobs_rotate(void *ctx, const vect_t rvec,
				     char origin, char coords,
				     const matp_t obj_rot,
				     const pointp_t pvt_pt)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_rotate(&n->view, rvec, origin, coords, obj_rot,
	    pvt_pt) : 0;
}

int
_rt_view_context_native_knobs_update_rate_flags(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_update_rate_flags(&n->view) : 0;
}

int
_rt_view_context_native_knob_values_get(struct rt_view_knob_values *values,
					const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    struct bv_knob_values bv_values;

    if (!values)
	return 0;

    memset(values, 0, sizeof(*values));
    if (!n || !bv_knob_values_get(&bv_values, &n->view))
	return 0;

    VMOVE(values->rate_rotation, bv_values.rate_rotation);
    VMOVE(values->rate_translation, bv_values.rate_translation);
    values->rate_scale = bv_values.rate_scale;
    VMOVE(values->absolute_rotation, bv_values.absolute_rotation);
    VMOVE(values->absolute_translation, bv_values.absolute_translation);
    values->absolute_scale = bv_values.absolute_scale;
    return 1;
}

int
_rt_view_context_native_knobs_calibrate(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    return n ? bv_knobs_calibrate(&n->view) : 0;
}

int
_rt_view_context_native_lod_policy_get(struct rt_view_lod_policy *policy,
				       const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    if (!policy)
	return 0;

    rt_view_lod_policy_init(policy);
    if (!n)
	return 0;

    *policy = n->lod_policy;
    rt_view_lod_policy_sanitize(policy);
    return 1;
}

int
_rt_view_context_native_lod_policy_apply(void *ctx,
					 const struct rt_view_lod_policy *policy)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n || !policy)
	return 0;

    n->lod_policy = *policy;
    rt_view_lod_policy_sanitize(&n->lod_policy);
    return 1;
}

int
_rt_view_context_native_lod_policy_copy(void *dst_ctx, const void *src_ctx)
{
    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    if (!_rt_view_context_native_lod_policy_get(&policy, src_ctx))
	return 0;
    return _rt_view_context_native_lod_policy_apply(dst_ctx, &policy);
}

int
_rt_view_context_native_lod_bounds_update(void *ctx)
{
    return native_ctx(ctx) ? 1 : 0;
}

int
_rt_view_context_native_lod_bounds_callback_set(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    if (!n)
	return 0;

    n->lod_bounds_callback_set = 1;
    return 1;
}

int
_rt_view_context_native_lod_bounds_callback_is(const void *ctx)
{
    const struct rt_view_context_native *n = native_ctx_const(ctx);
    return n ? n->lod_bounds_callback_set : 0;
}

struct rt_view_native_bounds_update_state {
    int callback_set;
};

void *
_rt_view_context_native_bounds_update_suspend(void *ctx)
{
    struct rt_view_context_native *n = native_ctx(ctx);
    struct rt_view_native_bounds_update_state *state;
    if (!n || !n->lod_bounds_callback_set)
	return NULL;

    BU_GET(state, struct rt_view_native_bounds_update_state);
    state->callback_set = n->lod_bounds_callback_set;
    n->lod_bounds_callback_set = 0;
    return state;
}

int
_rt_view_context_native_bounds_update_restore(void *ctx, void *state_ctx,
					      int UNUSED(refresh_bounds))
{
    struct rt_view_context_native *n = native_ctx(ctx);
    struct rt_view_native_bounds_update_state *state =
	(struct rt_view_native_bounds_update_state *)state_ctx;
    if (!state)
	return 0;

    if (n)
	n->lod_bounds_callback_set = state->callback_set;
    BU_PUT(state, struct rt_view_native_bounds_update_state);
    return n ? 1 : 0;
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
