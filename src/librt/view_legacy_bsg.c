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

#include <math.h>
#include <string.h>

#include "bg/plane.h"
#include "bn/qmath.h"
#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "bu/str.h"
#include "bu/vls.h"
#include "bsg/backend_scene.h"
#include "bsg/export.h"
#include "bsg/faceplate.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"
#include "bsg/interaction.h"
#include "bsg/lod.h"		/* bsg_view_bounds */
#include "bsg/measure.h"
#include "bsg/node.h"
#include "bsg/pick.h"
#include "bsg/polygon.h"
#include "bsg/render.h"
#include "bsg/render_item.h"
#include "bsg/selection.h"
#include "bsg/separator.h"
#include "bsg/snap.h"
#include "bsg/snap_action.h"
#include "bsg/util.h"
#include "bsg/view_set.h"
#include "bsg/view_state.h"
#include "rt/edit_legacy_bsg.h"
#include "rt/primitives/sketch_legacy_bsg.h"
#include "rt/view_legacy_bsg.h"
#include "view_context_private.h"

struct rt_view_line_set_context_bsg {
    bsg_line_set_ref lines;
    bsg_node_ref node;
};

struct rt_view_update_callback_context_bsg {
    rt_view_context_update_callback_t callback;
    void *data;
};

static void
rt_view_context_update_callback_bridge_bsg(struct bsg_view *v, void *data)
{
    struct rt_view_update_callback_context_bsg *state =
	(struct rt_view_update_callback_context_bsg *)data;

    if (state && state->callback)
	state->callback((void *)v, state->data);
}

static void
rt_view_update_callback_bridge_clear_bsg(struct bsg_view *v)
{
    if (!v || v->gv_callback != rt_view_context_update_callback_bridge_bsg)
	return;

    struct rt_view_update_callback_context_bsg *state =
	(struct rt_view_update_callback_context_bsg *)v->gv_clientData;
    if (state)
	BU_PUT(state, struct rt_view_update_callback_context_bsg);
    v->gv_callback = NULL;
    v->gv_clientData = NULL;
}

static bsg_scene_ref
rt_view_scene_ref_to_bsg(rt_view_scene_ref_bsg ref)
{
    bsg_scene_ref bsg_ref = BSG_SCENE_REF_NULL_INIT;
    bsg_ref.opaque = ref.opaque;
    return bsg_ref;
}

static rt_view_scene_ref_bsg
rt_view_scene_ref_from_bsg(bsg_scene_ref ref)
{
    rt_view_scene_ref_bsg rt_ref = RT_VIEW_SCENE_REF_BSG_NULL_INIT;
    rt_ref.opaque = ref.opaque;
    return rt_ref;
}

static bsg_scene_ref
rt_view_neutral_scene_ref_to_bsg(rt_view_scene_ref ref)
{
    bsg_scene_ref bsg_ref = BSG_SCENE_REF_NULL_INIT;
    if (ref.backend != RT_VIEW_SCENE_BACKEND_NONE &&
	ref.backend != RT_VIEW_SCENE_BACKEND_BSG)
	return bsg_ref;
    bsg_ref.opaque = ref.opaque;
    return bsg_ref;
}

static rt_view_scene_ref
rt_view_neutral_scene_ref_from_bsg(bsg_scene_ref ref)
{
    rt_view_scene_ref rt_ref = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_ref.opaque = ref.opaque;
    rt_ref.backend = ref.opaque ? RT_VIEW_SCENE_BACKEND_BSG :
		     RT_VIEW_SCENE_BACKEND_NONE;
    return rt_ref;
}

int
rt_view_context_is_bsg(const void *ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)ctx;
    return (v && v->magic == BSG_VIEW_MAGIC) ? 1 : 0;
}

int
rt_view_context_is_valid(const void *ctx)
{
    return (_rt_view_context_native_is(ctx) || rt_view_context_is_bsg(ctx)) ?
	   1 : 0;
}

int
rt_view_context_is_retained(const void *ctx)
{
    return rt_view_context_is_bsg(ctx);
}

void *
rt_view_context_user_data_from_bsg(const void *ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)ctx;
    return (v && v->magic == BSG_VIEW_MAGIC) ? v->u_data : NULL;
}

void *
rt_view_context_user_data_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_user_data_get(ctx);
    return rt_view_context_user_data_from_bsg(ctx);
}

int
rt_view_context_user_data_set_bsg(void *ctx, void *user_data)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    if (!v || v->magic != BSG_VIEW_MAGIC)
	return 0;

    v->u_data = user_data;
    return 1;
}

int
rt_view_context_user_data_set(void *ctx, void *user_data)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_user_data_set(ctx, user_data);
    return rt_view_context_user_data_set_bsg(ctx, user_data);
}

int
rt_view_context_tclcad_data_set_bsg(void *ctx, void *tcl_data)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    if (!v || v->magic != BSG_VIEW_MAGIC)
	return 0;

    v->gv_tcl = (struct bsg_data_tclcad *)tcl_data;
    return 1;
}

int
rt_view_context_tclcad_data_set(void *ctx, void *tcl_data)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_tclcad_data_set(ctx, tcl_data);
    return rt_view_context_tclcad_data_set_bsg(ctx, tcl_data);
}

void
rt_edit_view_from_bsg(struct rt_edit_view *ev, const void *view_ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)view_ctx;

    if (!ev)
	return;

    ev->gv_scale = rt_view_scale_from_bsg(v);
    ev->gv_base2local = rt_view_base2local_from_bsg(v);
    ev->gv_local2base = rt_view_local2base_from_bsg(v);
    ev->gv_coord = rt_view_coord_from_bsg(v);
    ev->gv_rotate_about = rt_view_rotate_about_from_bsg(v);
    rt_view_rotation_from_bsg(ev->gv_rotation, v);
    rt_view_center_from_bsg(ev->gv_center, v);
    rt_view_model2view_from_bsg(ev->gv_model2view, v);
    rt_view_view2model_from_bsg(ev->gv_view2model, v);
}

void
rt_edit_view_from_context(struct rt_edit_view *ev, const void *view_ctx)
{
    rt_edit_view_from_bsg(ev, view_ctx);
}

struct rt_edit *
rt_edit_create_bsg(struct db_full_path *dfp, struct db_i *dbip,
		   struct bn_tol *tol, const void *view_ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)view_ctx;
    struct rt_edit_view ev;

    if (!v)
	return rt_edit_create(dfp, dbip, tol, NULL);

    rt_edit_view_from_bsg(&ev, v);
    return rt_edit_create(dfp, dbip, tol, &ev);
}

struct rt_edit *
rt_edit_create_context(struct db_full_path *dfp, struct db_i *dbip,
		       struct bn_tol *tol, const void *view_ctx)
{
    return rt_edit_create_bsg(dfp, dbip, tol, view_ctx);
}

int
rt_edit_reinit_bsg(struct rt_edit *s, struct db_full_path *dfp,
		   struct db_i *dbip, struct bn_tol *tol,
		   const void *view_ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)view_ctx;
    struct rt_edit_view ev;

    if (!v)
	return rt_edit_reinit(s, dfp, dbip, tol, NULL);

    rt_edit_view_from_bsg(&ev, v);
    return rt_edit_reinit(s, dfp, dbip, tol, &ev);
}

int
rt_edit_reinit_context(struct rt_edit *s, struct db_full_path *dfp,
		       struct db_i *dbip, struct bn_tol *tol,
		       const void *view_ctx)
{
    return rt_edit_reinit_bsg(s, dfp, dbip, tol, view_ctx);
}

int
rt_edit_knob_cmd_process_bsg(
    struct rt_edit *s,
    vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, int *do_sca,
    const void *view_ctx, const char *cmd, fastf_t f,
    char origin, int incr_flag, void *u_data)
{
    const struct bsg_view *v = (const struct bsg_view *)view_ctx;
    struct rt_edit_view ev;

    if (!v)
	return rt_edit_knob_cmd_process(s, rvec, do_rot, tvec, do_tran,
					do_sca, NULL, cmd, f, origin, incr_flag, u_data);

    rt_edit_view_from_bsg(&ev, v);
    return rt_edit_knob_cmd_process(s, rvec, do_rot, tvec, do_tran,
				    do_sca, &ev, cmd, f, origin, incr_flag, u_data);
}

int
rt_edit_knob_cmd_process_context(
    struct rt_edit *s,
    vect_t *rvec, int *do_rot, vect_t *tvec, int *do_tran, int *do_sca,
    const void *view_ctx, const char *cmd, fastf_t f,
    char origin, int incr_flag, void *u_data)
{
    return rt_edit_knob_cmd_process_bsg(s, rvec, do_rot, tvec, do_tran,
					do_sca, view_ctx, cmd, f, origin, incr_flag, u_data);
}

int
rt_view_context_callbacks_set_bsg(void *ctx, struct bu_ptbl *callbacks)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    if (!v || v->magic != BSG_VIEW_MAGIC)
	return 0;

    if (v->callbacks && v->callbacks != callbacks) {
	bu_ptbl_free(v->callbacks);
	BU_PUT(v->callbacks, struct bu_ptbl);
    }

    v->callbacks = callbacks;
    return 1;
}

int
rt_view_context_callbacks_set(void *ctx, struct bu_ptbl *callbacks)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_callbacks_set(ctx, callbacks);
    return rt_view_context_callbacks_set_bsg(ctx, callbacks);
}

void *
rt_view_context_create_bsg(void)
{
    struct bsg_view *v = NULL;

    BU_ALLOC(v, struct bsg_view);
    bsg_view_init(v, NULL);
    return v;
}

void *
rt_view_context_create(void)
{
    return _rt_view_context_native_create();
}

void *
rt_view_context_create_with_set_bsg(void *view_set_ctx)
{
    struct bsg_view *v = NULL;

    BU_ALLOC(v, struct bsg_view);
    bsg_view_init(v, (struct bsg_view_set *)view_set_ctx);
    return v;
}

void *
rt_view_context_create_with_set(void *view_set_ctx)
{
    return _rt_view_context_native_create_with_set(view_set_ctx);
}

void *
rt_view_context_create_copy_with_set_bsg(const void *src_ctx,
	void *view_set_ctx)
{
    struct bsg_view *v = NULL;

    BU_ALLOC(v, struct bsg_view);
    rt_view_init_copy_bsg(v, (const struct bsg_view *)src_ctx,
			  (struct bsg_view_set *)view_set_ctx);
    return v;
}

static int
rt_view_context_native_copy_from_bsg(void *dst_ctx, const void *src_ctx)
{
    const struct bsg_view *src = (const struct bsg_view *)src_ctx;
    if (!_rt_view_context_native_is(dst_ctx))
	return 0;
    if (!src_ctx)
	return 1;
    if (!src || src->magic != BSG_VIEW_MAGIC)
	return 0;

    (void)_rt_view_context_native_name_set(dst_ctx,
					   bu_vls_cstr(&src->gv_name));
    (void)_rt_view_context_native_user_data_set(dst_ctx, src->u_data);
    (void)_rt_view_context_native_display_manager_set(dst_ctx, src->dmp);
    (void)_rt_view_context_native_dimensions_set(dst_ctx, src->gv_width,
	    src->gv_height);
    (void)_rt_view_context_native_scale_state_set(dst_ctx, src->gv_scale,
	    src->gv_i_scale, src->gv_a_scale, src->gv_size, src->gv_isize);
    (void)_rt_view_context_native_unit_conversion_set(dst_ctx,
	    src->gv_local2base, src->gv_base2local);
    (void)_rt_view_context_native_perspective_set(dst_ctx,
	    src->gv_perspective);
    (void)_rt_view_context_native_aet_state_set(dst_ctx, src->gv_aet);
    (void)_rt_view_context_native_eye_pos_set(dst_ctx, src->gv_eye_pos);
    (void)_rt_view_context_native_keypoint_set(dst_ctx, src->gv_keypoint);
    (void)_rt_view_context_native_rotate_about_set(dst_ctx,
	    src->gv_rotate_about);
    (void)_rt_view_context_native_coord_set(dst_ctx, src->gv_coord);
    (void)_rt_view_context_native_current_point_set(dst_ctx, src->gv_point);
    (void)_rt_view_context_native_previous_mouse_set(dst_ctx,
	    src->gv_prevMouseX, src->gv_prevMouseY);
    (void)_rt_view_context_native_mouse_state_set(dst_ctx, src->gv_mouse_x,
	    src->gv_mouse_y);
    (void)_rt_view_context_native_rotation_set(dst_ctx, src->gv_rotation);
    (void)_rt_view_context_native_center_mat_set(dst_ctx, src->gv_center);
    (void)_rt_view_context_native_model2view_set(dst_ctx, src->gv_model2view);
    (void)_rt_view_context_native_view2model_set(dst_ctx, src->gv_view2model);
    (void)_rt_view_context_native_pmodel2view_set(dst_ctx,
	    src->gv_pmodel2view);
    (void)_rt_view_context_native_pmat_set(dst_ctx, src->gv_pmat);
    return 1;
}

void *
rt_view_context_create_copy_with_set(const void *src_ctx, void *view_set_ctx)
{
    if (_rt_view_context_native_is(src_ctx))
	return _rt_view_context_native_create_copy_with_set(src_ctx,
		view_set_ctx);

    void *ctx = _rt_view_context_native_create_with_set(view_set_ctx);
    if (!rt_view_context_native_copy_from_bsg(ctx, src_ctx)) {
	rt_view_context_free(ctx);
	return NULL;
    }
    return ctx;
}

void
rt_view_context_free_bsg(void *ctx)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return;

    _rt_view_context_scene_adapter_clear(ctx);
    _rt_view_context_feature_adapter_clear(ctx);
    _rt_view_context_polygon_adapter_clear(ctx);
    _rt_view_context_selection_adapter_clear(ctx);
    rt_view_free_bsg(v);
    bu_free(v, "rt_view_context_create_bsg view");
}

void
rt_view_context_free(void *ctx)
{
    if (_rt_view_context_native_is(ctx)) {
	_rt_view_context_scene_adapter_clear(ctx);
	_rt_view_context_feature_adapter_clear(ctx);
	_rt_view_context_polygon_adapter_clear(ctx);
	_rt_view_context_selection_adapter_clear(ctx);
	_rt_view_context_native_free(ctx);
	return;
    }
    rt_view_context_free_bsg(ctx);
}

int
rt_view_context_release_storage_bsg(void *ctx)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return 0;

    _rt_view_context_scene_adapter_clear(ctx);
    _rt_view_context_feature_adapter_clear(ctx);
    _rt_view_context_polygon_adapter_clear(ctx);
    _rt_view_context_selection_adapter_clear(ctx);
    rt_view_update_callback_bridge_clear_bsg(v);
    BU_PUT(v, struct bsg_view);
    return 1;
}

int
rt_view_context_release_storage(void *ctx)
{
    if (_rt_view_context_native_is(ctx)) {
	_rt_view_context_scene_adapter_clear(ctx);
	_rt_view_context_feature_adapter_clear(ctx);
	_rt_view_context_polygon_adapter_clear(ctx);
	_rt_view_context_selection_adapter_clear(ctx);
	return _rt_view_context_native_release_storage(ctx);
    }
    return rt_view_context_release_storage_bsg(ctx);
}

int
rt_view_context_name_set_bsg(void *ctx, const char *name)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v || !name)
	return 0;

    bu_vls_sprintf(&v->gv_name, "%s", name);
    return 1;
}

const char *
rt_view_context_name_from_bsg(const void *ctx)
{
    return rt_view_name_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_name_set(void *ctx, const char *name)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_name_set(ctx, name);
    return rt_view_context_name_set_bsg(ctx, name);
}

const char *
rt_view_context_name_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_name_get(ctx);
    return rt_view_context_name_from_bsg(ctx);
}

int
rt_view_context_width_from_bsg(const void *ctx)
{
    return rt_view_width_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_width_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_width_get(ctx);
    return rt_view_context_width_from_bsg(ctx);
}

int
rt_view_context_height_from_bsg(const void *ctx)
{
    return rt_view_height_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_height_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_height_get(ctx);
    return rt_view_context_height_from_bsg(ctx);
}

int
rt_view_context_dimensions_set_bsg(void *ctx, int width, int height)
{
    return rt_view_dimensions_set_bsg((struct bsg_view *)ctx, width, height);
}

int
rt_view_context_dimensions_set(void *ctx, int width, int height)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_dimensions_set(ctx, width, height);
    return rt_view_context_dimensions_set_bsg(ctx, width, height);
}

void *
rt_view_context_display_manager_from_bsg(const void *ctx)
{
    return rt_view_display_manager_from_bsg((const struct bsg_view *)ctx);
}

void *
rt_view_context_display_manager_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_display_manager_get(ctx);
    return rt_view_context_display_manager_from_bsg(ctx);
}

int
rt_view_context_display_manager_set_bsg(void *ctx, void *dmp)
{
    return rt_view_display_manager_set_bsg((struct bsg_view *)ctx, dmp);
}

int
rt_view_context_display_manager_set(void *ctx, void *dmp)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_display_manager_set(ctx, dmp);
    return rt_view_context_display_manager_set_bsg(ctx, dmp);
}

int
rt_view_context_edit_matrix_set_bsg(void *ctx, matp_t edit_mat)
{
    return rt_view_edit_matrix_set_bsg((struct bsg_view *)ctx, edit_mat);
}

int
rt_view_context_edit_matrix_set(void *ctx, matp_t edit_mat)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_edit_matrix_set(ctx, edit_mat);
    return rt_view_context_edit_matrix_set_bsg(ctx, edit_mat);
}

int
rt_view_context_edit_matrix_clear_bsg(void *ctx)
{
    return rt_view_edit_matrix_clear_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_edit_matrix_clear(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_edit_matrix_clear(ctx);
    return rt_view_context_edit_matrix_clear_bsg(ctx);
}

uint64_t
rt_view_context_frame_revision_from_bsg(const void *ctx)
{
    return rt_view_frame_revision_from_bsg((const struct bsg_view *)ctx);
}

uint64_t
rt_view_context_frame_revision_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_frame_revision_get(ctx);
    return rt_view_context_frame_revision_from_bsg(ctx);
}

uint64_t
rt_view_context_frame_revision_bump(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_frame_revision_bump(ctx);

    struct bsg_view *v = (struct bsg_view *)ctx;
    if (!v)
	return 0;

    v->gv_frame_rev++;
    return v->gv_frame_rev;
}

int
rt_view_context_model_matrices_identity_bsg(void *ctx)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return 0;

    MAT_IDN(v->gv_model2view);
    MAT_IDN(v->gv_view2model);
    return 1;
}

int
rt_view_context_refresh_request_bsg(void *ctx, uint32_t flags)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    if (!v)
	return 0;

    bsg_view_refresh_request(v, flags);
    return 1;
}

int
rt_view_context_refresh_request(void *ctx, uint32_t flags)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_request(ctx, flags);
    return rt_view_context_refresh_request_bsg(ctx, flags);
}

int
rt_view_context_refresh_dirty_from_bsg(const void *ctx)
{
    return rt_view_refresh_dirty_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_refresh_dirty_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_dirty_get(ctx);
    return rt_view_context_refresh_dirty_from_bsg(ctx);
}

uint32_t
rt_view_context_refresh_consume_bsg(void *ctx)
{
    return rt_view_refresh_consume_bsg((struct bsg_view *)ctx);
}

uint32_t
rt_view_context_refresh_consume(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_consume(ctx);
    return rt_view_context_refresh_consume_bsg(ctx);
}

int
rt_view_context_refresh_complete_bsg(void *ctx)
{
    return rt_view_refresh_complete_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_refresh_complete(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_complete(ctx);
    return rt_view_context_refresh_complete_bsg(ctx);
}

int
rt_view_context_refresh_enabled_from_bsg(const void *ctx)
{
    return rt_view_refresh_enabled_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_refresh_enabled_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_enabled_get(ctx);
    return rt_view_context_refresh_enabled_from_bsg(ctx);
}

int
rt_view_context_refresh_enabled_set_bsg(void *ctx, int enabled)
{
    return rt_view_refresh_enabled_set_bsg((struct bsg_view *)ctx, enabled);
}

int
rt_view_context_refresh_enabled_set(void *ctx, int enabled)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_enabled_set(ctx, enabled);
    return rt_view_context_refresh_enabled_set_bsg(ctx, enabled);
}

int
rt_view_context_refresh_suppressed_from_bsg(const void *ctx)
{
    return rt_view_refresh_suppressed_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_refresh_suppressed_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_suppressed_get(ctx);
    return rt_view_context_refresh_suppressed_from_bsg(ctx);
}

int
rt_view_context_refresh_suppress_begin_bsg(void *ctx)
{
    return rt_view_refresh_suppress_begin_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_refresh_suppress_begin(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_suppress_begin(ctx);
    return rt_view_context_refresh_suppress_begin_bsg(ctx);
}

int
rt_view_context_refresh_suppress_end_bsg(void *ctx)
{
    return rt_view_refresh_suppress_end_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_refresh_suppress_end(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_suppress_end(ctx);
    return rt_view_context_refresh_suppress_end_bsg(ctx);
}

int
rt_view_context_refresh_drawn_count_from_bsg(const void *ctx)
{
    return rt_view_refresh_drawn_count_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_refresh_drawn_count_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_drawn_count_get(ctx);
    return rt_view_context_refresh_drawn_count_from_bsg(ctx);
}

int
rt_view_context_refresh_drawn_count_set_bsg(void *ctx, int count)
{
    return rt_view_refresh_drawn_count_set_bsg((struct bsg_view *)ctx, count);
}

int
rt_view_context_refresh_drawn_count_set(void *ctx, int count)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_refresh_drawn_count_set(ctx, count);
    return rt_view_context_refresh_drawn_count_set_bsg(ctx, count);
}

fastf_t
rt_view_context_scale_from_bsg(const void *ctx)
{
    return rt_view_scale_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_scale_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scale_get(ctx);
    return rt_view_context_scale_from_bsg(ctx);
}

fastf_t *
rt_view_context_scale_storage_from_bsg(void *ctx)
{
    return rt_view_scale_storage_from_bsg((struct bsg_view *)ctx);
}

fastf_t *
rt_view_context_scale_storage_get(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scale_storage_get(ctx);
    return rt_view_context_scale_storage_from_bsg(ctx);
}

fastf_t
rt_view_context_size_from_bsg(const void *ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)ctx;
    return v ? v->gv_size : 0.0;
}

fastf_t
rt_view_context_size_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_size_get(ctx);
    return rt_view_context_size_from_bsg(ctx);
}

int
rt_view_context_scale_state_set_bsg(void *ctx,
				    fastf_t scale,
				    fastf_t initial_scale,
				    fastf_t absolute_scale,
				    fastf_t size,
				    fastf_t inverse_size)
{
    return rt_view_scale_state_set_bsg((struct bsg_view *)ctx, scale,
				       initial_scale, absolute_scale, size, inverse_size);
}

int
rt_view_context_scale_state_set(void *ctx,
				fastf_t scale,
				fastf_t initial_scale,
				fastf_t absolute_scale,
				fastf_t size,
				fastf_t inverse_size)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scale_state_set(ctx, scale,
		initial_scale, absolute_scale, size, inverse_size);
    return rt_view_context_scale_state_set_bsg(ctx, scale, initial_scale,
	    absolute_scale, size, inverse_size);
}

fastf_t
rt_view_context_initial_scale_from_bsg(const void *ctx)
{
    return rt_view_initial_scale_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_initial_scale_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_initial_scale_get(ctx);
    return rt_view_context_initial_scale_from_bsg(ctx);
}

int
rt_view_context_initial_scale_set_bsg(void *ctx, fastf_t scale)
{
    return rt_view_initial_scale_set_bsg((struct bsg_view *)ctx, scale);
}

int
rt_view_context_initial_scale_set(void *ctx, fastf_t scale)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_initial_scale_set(ctx, scale);
    return rt_view_context_initial_scale_set_bsg(ctx, scale);
}

fastf_t
rt_view_context_absolute_scale_from_bsg(const void *ctx)
{
    return rt_view_absolute_scale_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_absolute_scale_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_absolute_scale_get(ctx);
    return rt_view_context_absolute_scale_from_bsg(ctx);
}

int
rt_view_context_absolute_scale_set_bsg(void *ctx, fastf_t scale)
{
    return rt_view_absolute_scale_set_bsg((struct bsg_view *)ctx, scale);
}

int
rt_view_context_absolute_scale_set(void *ctx, fastf_t scale)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_absolute_scale_set(ctx, scale);
    return rt_view_context_absolute_scale_set_bsg(ctx, scale);
}

int
rt_view_context_unit_conversion_set_bsg(void *ctx,
					fastf_t local2base,
					fastf_t base2local)
{
    return rt_view_unit_conversion_set_bsg((struct bsg_view *)ctx,
					   local2base, base2local);
}

int
rt_view_context_unit_conversion_set(void *ctx,
				    fastf_t local2base,
				    fastf_t base2local)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_unit_conversion_set(ctx, local2base,
		base2local);
    return rt_view_context_unit_conversion_set_bsg(ctx, local2base,
	    base2local);
}

fastf_t
rt_view_context_local2base_from_bsg(const void *ctx)
{
    return rt_view_local2base_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_local2base_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_local2base_get(ctx);
    return rt_view_context_local2base_from_bsg(ctx);
}

fastf_t
rt_view_context_base2local_from_bsg(const void *ctx)
{
    return rt_view_base2local_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_base2local_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_base2local_get(ctx);
    return rt_view_context_base2local_from_bsg(ctx);
}

int
rt_view_context_frametime_set_bsg(void *ctx, uint64_t frametime)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return 0;

    bsg_view_set_frametime(v, frametime);
    return 1;
}

int
rt_view_context_frametime_set(void *ctx, uint64_t frametime)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_frametime_set(ctx, frametime);
    return rt_view_context_frametime_set_bsg(ctx, frametime);
}

uint64_t
rt_view_context_frametime_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_frametime_get(ctx);
    return bsg_view_frametime((const struct bsg_view *)ctx);
}

int
rt_view_context_lod_bounds_callback_is_bsg(const void *ctx)
{
    return rt_view_lod_bounds_callback_is_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_lod_bounds_callback_is(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_lod_bounds_callback_is(ctx);
    return rt_view_context_lod_bounds_callback_is_bsg(ctx);
}

int
rt_view_context_is_independent_bsg(const void *ctx)
{
    return rt_view_is_independent_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_is_independent(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_is_independent(ctx);
    return rt_view_context_is_independent_bsg(ctx);
}

int
rt_view_context_independent_scope_is_null_bsg(void *ctx, int create)
{
    return rt_view_independent_scope_is_null_bsg((struct bsg_view *)ctx,
	    create);
}

int
rt_view_context_independent_scope_is_null(void *ctx, int create)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_independent_scope_is_null(ctx, create);
    return rt_view_context_independent_scope_is_null_bsg(ctx, create);
}

void
rt_view_context_independent_scope_destroy_bsg(void *ctx)
{
    rt_view_independent_scope_destroy_bsg((struct bsg_view *)ctx);
}

void
rt_view_context_independent_scope_destroy(void *ctx)
{
    if (_rt_view_context_native_is(ctx)) {
	_rt_view_context_native_independent_scope_destroy(ctx);
	return;
    }
    rt_view_context_independent_scope_destroy_bsg(ctx);
}

int
rt_view_set_context_add_bsg(void *view_set_ctx, void *view_ctx)
{
    if (!view_set_ctx || !view_ctx)
	return 0;

    rt_view_set_add_view_bsg((struct bsg_view_set *)view_set_ctx,
			     (struct bsg_view *)view_ctx);
    return 1;
}

int
rt_view_set_context_add(void *view_set_ctx, void *view_ctx)
{
    if (_rt_view_set_context_native_is(view_set_ctx) ||
	_rt_view_context_native_is(view_ctx))
	return _rt_view_set_context_native_add(view_set_ctx, view_ctx);
    return rt_view_set_context_add_bsg(view_set_ctx, view_ctx);
}

int
rt_view_set_context_remove_bsg(void *view_set_ctx, void *view_ctx)
{
    if (!view_set_ctx)
	return 0;

    rt_view_set_remove_view_bsg((struct bsg_view_set *)view_set_ctx,
				(struct bsg_view *)view_ctx);
    return 1;
}

int
rt_view_set_context_remove(void *view_set_ctx, void *view_ctx)
{
    if (_rt_view_set_context_native_is(view_set_ctx) ||
	_rt_view_context_native_is(view_ctx))
	return _rt_view_set_context_native_remove(view_set_ctx, view_ctx);
    return rt_view_set_context_remove_bsg(view_set_ctx, view_ctx);
}

void *
rt_view_line_render_item_create_bsg(void *view_ctx,
				    uint64_t source_id,
				    const char *name,
				    uint64_t geometry_source_identity,
				    uint64_t geometry_revision)
{
    struct bsg_render_item *item = bsg_render_item_create();

    if (!item)
	return NULL;

    item->view = (struct bsg_view *)view_ctx;
    item->source.source_id = source_id;
    item->source.name = name;
    item->geometry.kind = BSG_RENDER_GEOMETRY_LINE_SET;
    item->geometry.source_identity = geometry_source_identity;
    item->geometry.revision = geometry_revision;
    return item;
}

void
rt_view_render_item_free_bsg(void *render_item_ctx)
{
    bsg_render_item_free((struct bsg_render_item *)render_item_ctx);
}

int
rt_view_annotation_curves_add_bsg(void *view_ctx, const char *name)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;
    bsg_separator_ref root;
    bsg_annotation_ref annotation;
    static point_t pts[7] = {
	{0.0, 0.0, 0.0},
	{1.0, 0.0, 0.0},
	{0.0, 1.0, 0.0},
	{1.0, 1.0, 0.0},
	{2.0, 0.0, 0.0},
	{2.0, 1.0, 0.0},
	{3.0, 0.0, 0.0}
    };
    static fastf_t knots[6] = {0.0, 0.0, 0.0, 1.0, 1.0, 1.0};
    static int nurb_controls[3] = {2, 3, 4};
    static int bezier_controls[4] = {3, 4, 5, 6};
    struct bsg_annotation_segment segs[4];
    mat_t model_mat, display_mat;

    if (!v)
	return 0;

    root = bsg_view_scene_separator_ref(v, 1);
    annotation = bsg_annotation_ref_create(v, name ? name : "rt_view_annotation_curves");
    if (bsg_node_ref_is_null(bsg_annotation_ref_as_node(annotation)))
	return 0;

    memset(segs, 0, sizeof(segs));
    segs[0].kind = BSG_ANNOTATION_SEGMENT_LINE;
    segs[0].data.line.start = 0;
    segs[0].data.line.end = 1;
    segs[1].kind = BSG_ANNOTATION_SEGMENT_CARC;
    segs[1].data.carc.start = 0;
    segs[1].data.carc.end = 1;
    segs[1].data.carc.radius = 1.0;
    segs[1].data.carc.center_is_left = 1;
    segs[2].kind = BSG_ANNOTATION_SEGMENT_NURB;
    segs[2].data.nurb.order = 3;
    segs[2].data.nurb.knot_count = 6;
    segs[2].data.nurb.knots = knots;
    segs[2].data.nurb.control_point_count = 3;
    segs[2].data.nurb.control_points = nurb_controls;
    segs[3].kind = BSG_ANNOTATION_SEGMENT_BEZIER;
    segs[3].data.bezier.degree = 3;
    segs[3].data.bezier.control_point_count = 4;
    segs[3].data.bezier.control_points = bezier_controls;

    MAT_IDN(model_mat);
    MAT_IDN(display_mat);
    if (!bsg_annotation_ref_set_record(annotation, "curve annotation",
				       BSG_ANNOTATION_SPACE_MODEL, pts[0], model_mat, display_mat,
				       (const point_t *)pts, 7, segs, 4))
	return 0;

    bsg_separator_ref_append_child(root, bsg_annotation_ref_as_node(annotation));
    return 1;
}

int
rt_view_annotation_display_text_add_bsg(void *view_ctx,
					const char *name,
					const char *text,
					fastf_t size,
					fastf_t rotation)
{
    struct bsg_view *v = (struct bsg_view *)view_ctx;
    bsg_separator_ref root;
    bsg_annotation_ref annotation;
    static point_t pts[1] = {{0.0, 0.0, 0.0}};
    struct bsg_annotation_segment seg;
    mat_t model_mat, display_mat;

    if (!v)
	return 0;

    root = bsg_view_scene_separator_ref(v, 1);
    annotation = bsg_annotation_ref_create(v, name ? name : "rt_view_annotation_text");
    if (bsg_node_ref_is_null(bsg_annotation_ref_as_node(annotation)))
	return 0;

    memset(&seg, 0, sizeof(seg));
    seg.kind = BSG_ANNOTATION_SEGMENT_TEXT;
    seg.data.text.ref_pt = 0;
    seg.data.text.relative_position = BSG_ANNOTATION_TEXT_POS_TOP_RIGHT;
    seg.data.text.text = (char *)(text ? text : "");
    seg.data.text.size = size;
    seg.data.text.rotation = rotation;

    MAT_IDN(model_mat);
    MAT_IDN(display_mat);
    if (!bsg_annotation_ref_set_record(annotation, "display text",
				       BSG_ANNOTATION_SPACE_DISPLAY, pts[0], model_mat, display_mat,
				       (const point_t *)pts, 1, &seg, 1))
	return 0;

    bsg_separator_ref_append_child(root, bsg_annotation_ref_as_node(annotation));
    return 1;
}

void *
rt_view_context_line_set_create_bsg(void *ctx,
				    const char *name,
				    unsigned char r,
				    unsigned char g,
				    unsigned char b)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    struct rt_view_line_set_context_bsg *line_set = NULL;
    bsg_separator_ref root;
    bsg_node_ref root_node;

    if (!v || !name)
	return NULL;

    BU_ALLOC(line_set, struct rt_view_line_set_context_bsg);
    line_set->lines = bsg_line_set_ref_create(v, name);
    line_set->node = bsg_line_set_ref_as_node(line_set->lines);
    if (bsg_node_ref_is_null(line_set->node))
	return line_set;

    root = bsg_view_scene_separator_ref(v, 1);
    root_node = bsg_separator_ref_as_node(root);
    if (!bsg_node_ref_is_null(root_node)) {
	bsg_node_ref_set_color(line_set->node, r, g, b);
	bsg_separator_ref_append_child(root, line_set->node);
    }

    return line_set;
}

void *
rt_view_context_line_set_create(void *ctx,
				const char *name,
				unsigned char r,
				unsigned char g,
				unsigned char b)
{
    if (_rt_view_context_native_is(ctx))
	return NULL;
    return rt_view_context_line_set_create_bsg(ctx, name, r, g, b);
}

int
rt_view_line_set_context_is_null_bsg(const void *line_set_ctx)
{
    const struct rt_view_line_set_context_bsg *line_set =
	(const struct rt_view_line_set_context_bsg *)line_set_ctx;

    return (!line_set || bsg_node_ref_is_null(line_set->node));
}

int
rt_view_line_set_context_is_null(const void *line_set_ctx)
{
    return rt_view_line_set_context_is_null_bsg(line_set_ctx);
}

int
rt_view_line_set_context_set_points_bsg(void *line_set_ctx,
					const point_t *points,
					const int *commands,
					size_t point_count)
{
    struct rt_view_line_set_context_bsg *line_set =
	(struct rt_view_line_set_context_bsg *)line_set_ctx;

    if (rt_view_line_set_context_is_null_bsg(line_set))
	return 0;

    return bsg_line_set_ref_set_points(line_set->lines, points, commands,
				       point_count);
}

int
rt_view_line_set_context_set_points(void *line_set_ctx,
				    const point_t *points,
				    const int *commands,
				    size_t point_count)
{
    return rt_view_line_set_context_set_points_bsg(line_set_ctx, points,
	    commands, point_count);
}

void
rt_view_line_set_context_destroy_bsg(void *line_set_ctx)
{
    struct rt_view_line_set_context_bsg *line_set =
	(struct rt_view_line_set_context_bsg *)line_set_ctx;

    if (!line_set)
	return;

    if (!bsg_node_ref_is_null(line_set->node))
	bsg_node_ref_destroy(line_set->node);
    bu_free(line_set, "RT view BSG line-set context");
}

void
rt_view_line_set_context_destroy(void *line_set_ctx)
{
    rt_view_line_set_context_destroy_bsg(line_set_ctx);
}

void *
rt_view_display_manager_from_bsg(const struct bsg_view *v)
{
    return v ? v->dmp : NULL;
}

int
rt_view_display_manager_set_bsg(struct bsg_view *v, void *dmp)
{
    if (!v)
	return 0;

    v->dmp = dmp;
    return 1;
}

int
rt_view_update_callback_set_bsg(struct bsg_view *v,
				rt_view_update_callback_bsg_t callback,
				void *data)
{
    if (!v)
	return 0;

    rt_view_update_callback_bridge_clear_bsg(v);

    if (callback && !v->callbacks) {
	BU_GET(v->callbacks, struct bu_ptbl);
	bu_ptbl_init(v->callbacks, 8, "rt view BSG callbacks");
    }

    v->gv_callback = callback;
    v->gv_clientData = data;
    return 1;
}

int
rt_view_context_update_callback_set_bsg(void *ctx,
					rt_view_update_callback_bsg_t callback,
					void *data)
{
    return rt_view_update_callback_set_bsg((struct bsg_view *)ctx,
					   callback, data);
}

int
rt_view_context_update_callback_set(void *ctx,
				    rt_view_context_update_callback_t callback,
				    void *data)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_update_callback_set(ctx, callback,
		data);

    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return 0;

    if (!callback) {
	rt_view_update_callback_bridge_clear_bsg(v);
	v->gv_callback = NULL;
	v->gv_clientData = NULL;
	return 1;
    }

    struct rt_view_update_callback_context_bsg *state;
    BU_GET(state, struct rt_view_update_callback_context_bsg);
    state->callback = callback;
    state->data = data;

    return rt_view_update_callback_set_bsg(v,
					   rt_view_context_update_callback_bridge_bsg, state);
}

int
rt_view_edit_matrix_set_bsg(struct bsg_view *v, matp_t edit_mat)
{
    if (!v)
	return 0;

    v->gv_edit_mat = edit_mat;
    return 1;
}

int
rt_view_edit_matrix_clear_bsg(struct bsg_view *v)
{
    return rt_view_edit_matrix_set_bsg(v, NULL);
}

uint64_t
rt_view_frame_revision_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_frame_rev : 0;
}

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

void
rt_view_context_info_from_bsg(struct rt_view_info *info, const void *ctx)
{
    rt_view_info_from_bsg(info, (const struct bsg_view *)ctx);
}

void
rt_view_context_info_get(struct rt_view_info *info, const void *ctx)
{
    if (_rt_view_context_native_is(ctx)) {
	_rt_view_context_native_info_get(info, ctx);
	return;
    }
    rt_view_context_info_from_bsg(info, ctx);
}

int
rt_view_context_obb_from_bsg(const void *ctx,
			     point_t center,
			     vect_t extent1,
			     vect_t extent2,
			     vect_t extent3)
{
    const struct bsg_view *v = (const struct bsg_view *)ctx;

    if (!v || v->magic != BSG_VIEW_MAGIC || !center || !extent1 || !extent2 || !extent3)
	return 0;

    VMOVE(center, v->obb_center);
    VMOVE(extent1, v->obb_extent1);
    VMOVE(extent2, v->obb_extent2);
    VMOVE(extent3, v->obb_extent3);
    return 1;
}

int
rt_view_context_obb_get(const void *ctx,
			point_t center,
			vect_t extent1,
			vect_t extent2,
			vect_t extent3)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_obb_get(ctx, center, extent1, extent2,
					       extent3);
    return rt_view_context_obb_from_bsg(ctx, center, extent1, extent2, extent3);
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

fastf_t
rt_view_context_radius_from_bsg(const void *ctx)
{
    return rt_view_radius_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_radius_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_radius_get(ctx);
    return rt_view_context_radius_from_bsg(ctx);
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
rt_view_context_screen_to_view_from_bsg(fastf_t *fx,
					fastf_t *fy,
					void *ctx,
					fastf_t x,
					fastf_t y)
{
    return rt_view_screen_to_view_from_bsg(fx, fy,
					   (struct bsg_view *)ctx, x, y);
}

int
rt_view_context_screen_to_view(fastf_t *fx,
			       fastf_t *fy,
			       void *ctx,
			       fastf_t x,
			       fastf_t y)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_screen_to_view(fx, fy, ctx, x, y);
    return rt_view_context_screen_to_view_from_bsg(fx, fy, ctx, x, y);
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
rt_view_context_screen_point_from_bsg(point_t point,
				      void *ctx,
				      fastf_t x,
				      fastf_t y)
{
    return rt_view_screen_point_from_bsg(point, (struct bsg_view *)ctx, x, y);
}

int
rt_view_context_screen_point_get(point_t point, void *ctx, fastf_t x, fastf_t y)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_screen_point_get(point, ctx, x, y);
    return rt_view_context_screen_point_from_bsg(point, ctx, x, y);
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
rt_view_context_current_point_from_bsg(point_t point, const void *ctx)
{
    return rt_view_current_point_from_bsg(point,
					  (const struct bsg_view *)ctx);
}

int
rt_view_context_current_point_get(point_t point, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_current_point_get(point, ctx);
    return rt_view_context_current_point_from_bsg(point, ctx);
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
rt_view_context_current_point_set_bsg(void *ctx, const point_t point)
{
    return rt_view_current_point_set_bsg((struct bsg_view *)ctx, point);
}

int
rt_view_context_current_point_set(void *ctx, const point_t point)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_current_point_set(ctx, point);
    return rt_view_context_current_point_set_bsg(ctx, point);
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
rt_view_context_previous_mouse_from_bsg(fastf_t *x, fastf_t *y, const void *ctx)
{
    return rt_view_previous_mouse_from_bsg(x, y, (const struct bsg_view *)ctx);
}

int
rt_view_context_previous_mouse_get(fastf_t *x, fastf_t *y, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_previous_mouse_get(x, y, ctx);
    return rt_view_context_previous_mouse_from_bsg(x, y, ctx);
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
rt_view_context_previous_mouse_set_bsg(void *ctx, fastf_t x, fastf_t y)
{
    return rt_view_previous_mouse_set_bsg((struct bsg_view *)ctx, x, y);
}

int
rt_view_context_previous_mouse_set(void *ctx, fastf_t x, fastf_t y)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_previous_mouse_set(ctx, x, y);
    return rt_view_context_previous_mouse_set_bsg(ctx, x, y);
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
rt_view_context_mouse_delta_settings_from_bsg(struct rt_view_mouse_delta_settings *settings,
	const void *ctx)
{
    return rt_view_mouse_delta_settings_from_bsg(settings,
	    (const struct bsg_view *)ctx);
}

int
rt_view_context_mouse_delta_settings_get(struct rt_view_mouse_delta_settings *settings,
	const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_mouse_delta_settings_get(settings, ctx);
    return rt_view_context_mouse_delta_settings_from_bsg(settings, ctx);
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
rt_view_context_mouse_state_set_bsg(void *ctx, int x, int y)
{
    return rt_view_mouse_state_set_bsg((struct bsg_view *)ctx, x, y);
}

int
rt_view_context_mouse_state_set(void *ctx, int x, int y)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_mouse_state_set(ctx, x, y);
    return rt_view_context_mouse_state_set_bsg(ctx, x, y);
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

int
rt_view_context_unique_object_name_bsg(struct bu_vls *oname,
				       const char *seed,
				       void *ctx)
{
    return rt_view_unique_object_name_bsg(oname, seed,
					  (struct bsg_view *)ctx);
}

int
rt_view_context_unique_object_name(struct bu_vls *oname,
				   const char *seed,
				   void *ctx)
{
    return rt_view_context_unique_object_name_bsg(oname, seed, ctx);
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
rt_view_context_interactive_rect_state_from_bsg(
    struct rt_view_interactive_rect_state *record,
    const void *ctx)
{
    return rt_view_interactive_rect_state_from_bsg(record,
	    (const struct bsg_view *)ctx);
}

int
rt_view_context_interactive_rect_state_get(
    struct rt_view_interactive_rect_state *record,
    const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_interactive_rect_state_get(record, ctx);
    return rt_view_context_interactive_rect_state_from_bsg(record, ctx);
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

int
rt_view_context_interactive_rect_state_set_bsg(
    void *ctx,
    const struct rt_view_interactive_rect_state *record)
{
    return rt_view_interactive_rect_state_set_bsg((struct bsg_view *)ctx,
	    record);
}

int
rt_view_context_interactive_rect_state_set(
    void *ctx,
    const struct rt_view_interactive_rect_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_interactive_rect_state_set(ctx,
		record);
    return rt_view_context_interactive_rect_state_set_bsg(ctx, record);
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
rt_view_context_adc_state_from_bsg(struct rt_view_adc_state *record,
				   const void *ctx)
{
    return rt_view_adc_state_from_bsg(record, (const struct bsg_view *)ctx);
}

int
rt_view_context_adc_state_get(struct rt_view_adc_state *record, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_adc_state_get(record, ctx);
    return rt_view_context_adc_state_from_bsg(record, ctx);
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

int
rt_view_context_adc_state_set_bsg(void *ctx,
				  const struct rt_view_adc_state *record)
{
    return rt_view_adc_state_set_bsg((struct bsg_view *)ctx, record);
}

int
rt_view_context_adc_state_set(void *ctx, const struct rt_view_adc_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_adc_state_set(ctx, record);
    return rt_view_context_adc_state_set_bsg(ctx, record);
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
rt_view_context_grid_state_from_bsg(struct rt_view_grid_state *record,
				    const void *ctx)
{
    return rt_view_grid_state_from_bsg(record, (const struct bsg_view *)ctx);
}

int
rt_view_context_grid_state_get(struct rt_view_grid_state *record, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_grid_state_get(record, ctx);
    return rt_view_context_grid_state_from_bsg(record, ctx);
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

int
rt_view_context_grid_state_set_bsg(void *ctx,
				   const struct rt_view_grid_state *record)
{
    return rt_view_grid_state_set_bsg((struct bsg_view *)ctx, record);
}

int
rt_view_context_grid_state_set(void *ctx, const struct rt_view_grid_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_grid_state_set(ctx, record);
    return rt_view_context_grid_state_set_bsg(ctx, record);
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
rt_view_context_model_axes_state_from_bsg(struct rt_view_axes_state *record,
	const void *ctx)
{
    return rt_view_model_axes_state_from_bsg(record,
	    (const struct bsg_view *)ctx);
}

int
rt_view_context_model_axes_state_get(struct rt_view_axes_state *record,
				     const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_model_axes_state_get(record, ctx);
    return rt_view_context_model_axes_state_from_bsg(record, ctx);
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
rt_view_context_model_axes_state_set_bsg(void *ctx,
	const struct rt_view_axes_state *record)
{
    return rt_view_model_axes_state_set_bsg((struct bsg_view *)ctx, record);
}

int
rt_view_context_model_axes_state_set(void *ctx,
				     const struct rt_view_axes_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_model_axes_state_set(ctx, record);
    return rt_view_context_model_axes_state_set_bsg(ctx, record);
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
rt_view_context_view_axes_state_from_bsg(struct rt_view_axes_state *record,
	const void *ctx)
{
    return rt_view_view_axes_state_from_bsg(record,
					    (const struct bsg_view *)ctx);
}

int
rt_view_context_view_axes_state_get(struct rt_view_axes_state *record,
				    const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_view_axes_state_get(record, ctx);
    return rt_view_context_view_axes_state_from_bsg(record, ctx);
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

int
rt_view_context_view_axes_state_set_bsg(void *ctx,
					const struct rt_view_axes_state *record)
{
    return rt_view_view_axes_state_set_bsg((struct bsg_view *)ctx, record);
}

int
rt_view_context_view_axes_state_set(void *ctx,
				    const struct rt_view_axes_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_view_axes_state_set(ctx, record);
    return rt_view_context_view_axes_state_set_bsg(ctx, record);
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
rt_view_context_center_dot_state_from_bsg(struct rt_view_other_state *record,
	const void *ctx)
{
    return rt_view_center_dot_state_from_bsg(record, (const struct bsg_view *)ctx);
}

int
rt_view_context_center_dot_state_get(struct rt_view_other_state *record,
				     const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_center_dot_state_get(record, ctx);
    return rt_view_context_center_dot_state_from_bsg(record, ctx);
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
rt_view_context_center_dot_state_set_bsg(void *ctx,
	const struct rt_view_other_state *record)
{
    return rt_view_center_dot_state_set_bsg((struct bsg_view *)ctx, record);
}

int
rt_view_context_center_dot_state_set(void *ctx,
				     const struct rt_view_other_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_center_dot_state_set(ctx, record);
    return rt_view_context_center_dot_state_set_bsg(ctx, record);
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
rt_view_context_scale_overlay_state_from_bsg(struct rt_view_other_state *record,
	const void *ctx)
{
    return rt_view_scale_overlay_state_from_bsg(record,
	    (const struct bsg_view *)ctx);
}

int
rt_view_context_scale_overlay_state_get(struct rt_view_other_state *record,
					const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scale_overlay_state_get(record,
		ctx);
    return rt_view_context_scale_overlay_state_from_bsg(record, ctx);
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

int
rt_view_context_scale_overlay_state_set_bsg(void *ctx,
	const struct rt_view_other_state *record)
{
    return rt_view_scale_overlay_state_set_bsg((struct bsg_view *)ctx,
	    record);
}

int
rt_view_context_scale_overlay_state_set(void *ctx,
					const struct rt_view_other_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scale_overlay_state_set(ctx,
		record);
    return rt_view_context_scale_overlay_state_set_bsg(ctx, record);
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
rt_view_context_params_state_from_bsg(struct rt_view_params_state *record,
				      const void *ctx)
{
    return rt_view_params_state_from_bsg(record, (const struct bsg_view *)ctx);
}

int
rt_view_context_params_state_get(struct rt_view_params_state *record,
				 const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_params_state_get(record, ctx);
    return rt_view_context_params_state_from_bsg(record, ctx);
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
rt_view_context_params_state_set_bsg(void *ctx,
				     const struct rt_view_params_state *record)
{
    return rt_view_params_state_set_bsg((struct bsg_view *)ctx, record);
}

int
rt_view_context_params_state_set(void *ctx,
				 const struct rt_view_params_state *record)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_params_state_set(ctx, record);
    return rt_view_context_params_state_set_bsg(ctx, record);
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
rt_view_context_orientation_quat_from_bsg(quat_t orientation, const void *ctx)
{
    return rt_view_orientation_quat_from_bsg(orientation,
	    (const struct bsg_view *)ctx);
}

int
rt_view_context_orientation_quat_get(quat_t orientation, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_orientation_quat_get(orientation, ctx);
    return rt_view_context_orientation_quat_from_bsg(orientation, ctx);
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
rt_view_context_aet_from_bsg(vect_t aet, const void *ctx)
{
    return rt_view_aet_from_bsg(aet, (const struct bsg_view *)ctx);
}

int
rt_view_context_aet_get(vect_t aet, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_aet_get(aet, ctx);
    return rt_view_context_aet_from_bsg(aet, ctx);
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
rt_view_context_aet_set_bsg(void *ctx, const vect_t aet)
{
    return rt_view_aet_set_bsg((struct bsg_view *)ctx, aet);
}

int
rt_view_context_aet_set(void *ctx, const vect_t aet)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_aet_set(ctx, aet);
    return rt_view_context_aet_set_bsg(ctx, aet);
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
rt_view_context_aet_state_set_bsg(void *ctx, const vect_t aet)
{
    return rt_view_aet_state_set_bsg((struct bsg_view *)ctx, aet);
}

int
rt_view_context_aet_state_set(void *ctx, const vect_t aet)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_aet_state_set(ctx, aet);
    return rt_view_context_aet_state_set_bsg(ctx, aet);
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
rt_view_context_update_bsg(void *ctx)
{
    return rt_view_update_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_update(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_update(ctx);
    return rt_view_context_update_bsg(ctx);
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
rt_view_context_autoview_bsg(void *ctx, fastf_t scale, int all_view_objs)
{
    return rt_view_autoview_bsg((struct bsg_view *)ctx, scale,
				all_view_objs);
}

int
rt_view_context_autoview(void *ctx, fastf_t scale, int all_view_objs)
{
    if (_rt_view_context_native_is(ctx))
	return 0;
    return rt_view_context_autoview_bsg(ctx, scale, all_view_objs);
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
rt_view_context_autoview_bounds_bsg(void *ctx,
				    fastf_t scale,
				    const point_t min,
				    const point_t max)
{
    return rt_view_autoview_bounds_bsg((struct bsg_view *)ctx, scale, min, max);
}

int
rt_view_context_autoview_bounds(void *ctx,
				fastf_t scale,
				const point_t min,
				const point_t max)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_autoview_bounds(ctx, scale, min, max);
    return rt_view_context_autoview_bounds_bsg(ctx, scale, min, max);
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

int
rt_view_context_adjust_bsg(void *ctx,
			   int dx,
			   int dy,
			   point_t keypoint,
			   int mode,
			   unsigned long long flags)
{
    return rt_view_adjust_bsg((struct bsg_view *)ctx, dx, dy, keypoint, mode,
			      flags);
}

int
rt_view_context_adjust(void *ctx,
		       int dx,
		       int dy,
		       point_t keypoint,
		       int mode,
		       unsigned long long flags)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_adjust(ctx, dx, dy, keypoint, mode,
					      flags);
    return rt_view_context_adjust_bsg(ctx, dx, dy, keypoint, mode, flags);
}

unsigned long long
rt_view_hash_bsg(const struct bsg_view *v)
{
    return v ? bsg_hash((struct bsg_view *)v) : 0ULL;
}

unsigned long long
rt_view_context_hash_bsg(const void *ctx)
{
    return rt_view_hash_bsg((const struct bsg_view *)ctx);
}

unsigned long long
rt_view_context_hash(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_hash(ctx);
    return rt_view_context_hash_bsg(ctx);
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

static struct bsg_snap_result *
rt_view_snap_result_to_bsg(struct rt_view_snap_result_bsg *result)
{
    return (struct bsg_snap_result *)result;
}

static const struct bsg_snap_result *
rt_view_snap_result_const_to_bsg(const struct rt_view_snap_result_bsg *result)
{
    return (const struct bsg_snap_result *)result;
}

struct rt_view_snap_result_bsg *
rt_view_snap_result_create_bsg(void)
{
    struct bsg_snap_result *result = NULL;

    BU_ALLOC(result, struct bsg_snap_result);
    bsg_snap_result_init(result);
    return (struct rt_view_snap_result_bsg *)result;
}

void
rt_view_snap_result_free_bsg(struct rt_view_snap_result_bsg *result)
{
    struct bsg_snap_result *bsg_result = rt_view_snap_result_to_bsg(result);

    if (!bsg_result)
	return;

    bsg_snap_result_free(bsg_result);
    BU_PUT(bsg_result, struct bsg_snap_result);
}

size_t
rt_view_snap_result_count_bsg(const struct rt_view_snap_result_bsg *result)
{
    const struct bsg_snap_result *bsg_result =
	rt_view_snap_result_const_to_bsg(result);

    return bsg_result ? bsg_result->sr_cnt : 0;
}

int
rt_view_snap_result_point_bsg(const struct rt_view_snap_result_bsg *result,
			      size_t index,
			      point_t point_out)
{
    const struct bsg_snap_result *bsg_result =
	rt_view_snap_result_const_to_bsg(result);

    if (!bsg_result || !point_out || index >= bsg_result->sr_cnt)
	return 0;

    VMOVE(point_out, bsg_result->sr_candidates[index].sc_point);
    return 1;
}

fastf_t
rt_view_snap_result_distance_bsg(const struct rt_view_snap_result_bsg *result,
				 size_t index)
{
    const struct bsg_snap_result *bsg_result =
	rt_view_snap_result_const_to_bsg(result);

    if (!bsg_result || index >= bsg_result->sr_cnt)
	return 0.0;

    return bsg_result->sr_candidates[index].sc_distance;
}

unsigned long long
rt_view_snap_result_kind_bsg(const struct rt_view_snap_result_bsg *result,
			     size_t index)
{
    const struct bsg_snap_result *bsg_result =
	rt_view_snap_result_const_to_bsg(result);

    if (!bsg_result || index >= bsg_result->sr_cnt)
	return 0ULL;

    return (unsigned long long)bsg_result->sr_candidates[index].sc_kind;
}

int
rt_view_snap_result_source_path_bsg(const struct rt_view_snap_result_bsg *result,
				    size_t index,
				    struct bu_vls *path_out)
{
    const struct bsg_snap_result *bsg_result =
	rt_view_snap_result_const_to_bsg(result);

    if (!bsg_result || !path_out || index >= bsg_result->sr_cnt)
	return 0;

    bu_vls_sprintf(path_out, "%s",
		   bu_vls_cstr(&bsg_result->sr_candidates[index].sc_source_path));
    return 1;
}

int
rt_view_snap_candidates_result_bsg(struct bsg_view *v,
				   point_t sample,
				   double tol,
				   unsigned long long kinds,
				   struct rt_view_snap_result_bsg *out)
{
    struct bsg_snap_result *bsg_out = rt_view_snap_result_to_bsg(out);

    if (!bsg_out)
	return 0;

    return bsg_snap_candidates(v, sample, tol, (bsg_snap_kind_mask)kinds,
			       bsg_out);
}

int
rt_view_context_snap_candidates_result_bsg(void *ctx,
	point_t sample,
	double tol,
	unsigned long long kinds,
	struct rt_view_snap_result_bsg *out)
{
    return rt_view_snap_candidates_result_bsg((struct bsg_view *)ctx,
	    sample, tol, kinds, out);
}

static int
rt_view_native_snap_grid_2d(void *ctx, fastf_t *vx, fastf_t *vy)
{
    struct rt_view_grid_state grid = RT_VIEW_GRID_STATE_INIT;
    mat_t model2view;
    point_t anchor_view = VINIT_ZERO;
    fastf_t sf = 0.0;
    fastf_t step_h = 0.0;
    fastf_t step_v = 0.0;
    fastf_t qx = 0.0;
    fastf_t qy = 0.0;
    fastf_t ax = 0.0;
    fastf_t ay = 0.0;

    if (!ctx || !vx || !vy)
	return 0;
    if (!_rt_view_context_native_grid_state_get(&grid, ctx))
	return 0;
    if (ZERO(grid.res_h) || ZERO(grid.res_v))
	return 0;

    sf = _rt_view_context_native_scale_get(ctx) *
	 _rt_view_context_native_base2local_get(ctx);
    if (ZERO(sf))
	return 0;

    step_h = grid.res_h * _rt_view_context_native_base2local_get(ctx);
    step_v = grid.res_v * _rt_view_context_native_base2local_get(ctx);
    if (ZERO(step_h) || ZERO(step_v))
	return 0;

    MAT_IDN(model2view);
    (void)_rt_view_context_native_model2view_get(model2view, ctx);
    MAT4X3PNT(anchor_view, model2view, grid.anchor);

    qx = (*vx) * sf;
    qy = (*vy) * sf;
    ax = anchor_view[X] * sf;
    ay = anchor_view[Y] * sf;

    *vx = (ax + round((qx - ax) / step_h) * step_h) / sf;
    *vy = (ay + round((qy - ay) / step_v) * step_v) / sf;
    return 1;
}

static void
rt_view_native_snap_result_append_grid(struct rt_view_snap_result_bsg *out,
				       const point_t point,
				       fastf_t dist)
{
    struct bsg_snap_result *bsg_out = rt_view_snap_result_to_bsg(out);
    struct bsg_snap_candidate *candidate = NULL;

    if (!bsg_out)
	return;

    bsg_out->sr_candidates = (struct bsg_snap_candidate *)bu_realloc(
				 bsg_out->sr_candidates,
				 (bsg_out->sr_cnt + 1) * sizeof(struct bsg_snap_candidate),
				 "native grid snap candidate");
    candidate = &bsg_out->sr_candidates[bsg_out->sr_cnt++];
    memset(candidate, 0, sizeof(*candidate));
    VMOVE(candidate->sc_point, point);
    candidate->sc_kind = BSG_SNAP_KIND_GRID;
    candidate->sc_distance = dist;
    candidate->sc_feature = (bsg_feature_ref)BSG_FEATURE_REF_NULL_INIT;
    bu_vls_init(&candidate->sc_source_path);
    bu_vls_strcpy(&candidate->sc_source_path, "grid");
}

static int
rt_view_context_snap_candidates_result_native(void *ctx,
	const point_t sample,
	double tol,
	unsigned long long kinds,
	struct rt_view_snap_result_bsg *out)
{
    mat_t model2view;
    mat_t view2model;
    point_t view_pt = VINIT_ZERO;
    point_t snapped = VINIT_ZERO;
    fastf_t vx = 0.0;
    fastf_t vy = 0.0;
    fastf_t dist = 0.0;

    if (!ctx || !sample || !out || !(kinds & RT_VIEW_SNAP_KIND_GRID))
	return 0;

    MAT_IDN(model2view);
    MAT_IDN(view2model);
    (void)_rt_view_context_native_model2view_get(model2view, ctx);
    (void)_rt_view_context_native_view2model_get(view2model, ctx);
    MAT4X3PNT(view_pt, model2view, sample);
    vx = view_pt[X];
    vy = view_pt[Y];
    if (!rt_view_native_snap_grid_2d(ctx, &vx, &vy))
	return 0;

    VSET(view_pt, vx, vy, view_pt[Z]);
    MAT4X3PNT(snapped, view2model, view_pt);
    dist = DIST_PNT_PNT(sample, snapped);
    if (tol > 0.0 && dist > tol)
	return 0;

    rt_view_native_snap_result_append_grid(out, snapped, dist);
    return 1;
}

void *
rt_view_snap_result_context_create(void)
{
    return (void *)rt_view_snap_result_create_bsg();
}

void
rt_view_snap_result_context_free(void *result_ctx)
{
    rt_view_snap_result_free_bsg((struct rt_view_snap_result_bsg *)result_ctx);
}

size_t
rt_view_snap_result_context_count(const void *result_ctx)
{
    return rt_view_snap_result_count_bsg(
	       (const struct rt_view_snap_result_bsg *)result_ctx);
}

int
rt_view_snap_result_context_point(const void *result_ctx,
				  size_t index,
				  point_t point_out)
{
    return rt_view_snap_result_point_bsg(
	       (const struct rt_view_snap_result_bsg *)result_ctx,
	       index, point_out);
}

fastf_t
rt_view_snap_result_context_distance(const void *result_ctx, size_t index)
{
    return rt_view_snap_result_distance_bsg(
	       (const struct rt_view_snap_result_bsg *)result_ctx, index);
}

unsigned long long
rt_view_snap_result_context_kind(const void *result_ctx, size_t index)
{
    return rt_view_snap_result_kind_bsg(
	       (const struct rt_view_snap_result_bsg *)result_ctx, index);
}

int
rt_view_snap_result_context_source_path(const void *result_ctx,
					size_t index,
					struct bu_vls *path_out)
{
    return rt_view_snap_result_source_path_bsg(
	       (const struct rt_view_snap_result_bsg *)result_ctx,
	       index, path_out);
}

int
rt_view_context_snap_candidates_result(void *ctx,
				       point_t sample,
				       double tol,
				       unsigned long long kinds,
				       void *out_ctx)
{
    if (_rt_view_context_native_is(ctx))
	return rt_view_context_snap_candidates_result_native(ctx, sample, tol,
		kinds, (struct rt_view_snap_result_bsg *)out_ctx);
    return rt_view_context_snap_candidates_result_bsg(ctx, sample, tol, kinds,
	    (struct rt_view_snap_result_bsg *)out_ctx);
}

int
rt_view_context_snap_first_candidate(void *ctx,
				     const point_t sample,
				     unsigned long long kinds,
				     point_t candidate)
{
    struct bsg_snap_result sres;
    point_t bsg_sample = VINIT_ZERO;
    int ret = 0;

    if (!candidate)
	return 0;

    VSET(candidate, 0.0, 0.0, 0.0);

    if (!ctx || !sample || !kinds)
	return 0;

    if (_rt_view_context_native_is(ctx)) {
	struct rt_view_snap_result_bsg *native_result =
	    rt_view_snap_result_create_bsg();
	if (rt_view_context_snap_candidates_result_native(ctx,
		sample, 0.0, kinds, native_result) > 0) {
	    ret = rt_view_snap_result_point_bsg(native_result, 0, candidate);
	}
	rt_view_snap_result_free_bsg(native_result);
	return ret;
    }

    VMOVE(bsg_sample, sample);
    bsg_snap_result_init(&sres);
    if (rt_view_snap_candidates_bsg((struct bsg_view *)ctx, bsg_sample, 0.0,
				    kinds, &sres) > 0) {
	VMOVE(candidate, sres.sr_candidates[0].sc_point);
	ret = 1;
    }
    bsg_snap_result_free(&sres);
    return ret;
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
rt_view_context_snap_point_2d_bsg(void *ctx,
				  fastf_t *vx,
				  fastf_t *vy,
				  unsigned long long kinds)
{
    return rt_view_snap_point_2d_bsg((struct bsg_view *)ctx, vx, vy, kinds);
}

int
rt_view_context_snap_point_2d(void *ctx,
			      fastf_t *vx,
			      fastf_t *vy,
			      unsigned long long kinds)
{
    if (_rt_view_context_native_is(ctx))
	return (kinds & RT_VIEW_SNAP_KIND_GRID) ?
	       rt_view_native_snap_grid_2d(ctx, vx, vy) : 0;
    return rt_view_context_snap_point_2d_bsg(ctx, vx, vy, kinds);
}

int
rt_view_snap_grid_2d_bsg(struct bsg_view *v, fastf_t *vx, fastf_t *vy)
{
    if (!v || !vx || !vy)
	return 0;

    return bsg_snap_grid_2d(v, vx, vy);
}

int
rt_view_context_snap_grid_2d_bsg(void *ctx, fastf_t *vx, fastf_t *vy)
{
    return rt_view_snap_grid_2d_bsg((struct bsg_view *)ctx, vx, vy);
}

int
rt_view_context_snap_grid_2d(void *ctx, fastf_t *vx, fastf_t *vy)
{
    if (_rt_view_context_native_is(ctx))
	return rt_view_native_snap_grid_2d(ctx, vx, vy);
    return rt_view_context_snap_grid_2d_bsg(ctx, vx, vy);
}

static void
_rt_view_measure_from_bsg_result(struct rt_view_measure_result *dst,
				 const struct bsg_measure_result *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    dst->distance = src->mr_distance;
    dst->projection = src->mr_projection;
    dst->normal_alignment = src->mr_normal_alignment;
    dst->valid = src->mr_valid;
}

int
rt_view_measure_candidates_bsg(struct bsg_view *v,
			       point_t a,
			       point_t b,
			       struct rt_view_measure_result *out)
{
    struct bsg_measure_result bsg_result = {0.0, 0.0, 0.0, 0};
    int ret = 0;

    if (out)
	memset(out, 0, sizeof(*out));
    if (!out)
	return 0;

    ret = bsg_measure_candidates(v, a, b, &bsg_result);
    _rt_view_measure_from_bsg_result(out, &bsg_result);
    return ret;
}

int
rt_view_context_measure_candidates_bsg(void *ctx,
				       point_t a,
				       point_t b,
				       struct rt_view_measure_result *out)
{
    return rt_view_measure_candidates_bsg((struct bsg_view *)ctx, a, b, out);
}

int
rt_view_context_measure_candidates(void *ctx,
				   point_t a,
				   point_t b,
				   struct rt_view_measure_result *out)
{
    if (out)
	memset(out, 0, sizeof(*out));
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_measure_candidates(ctx, a, b, out);
    return rt_view_context_measure_candidates_bsg(ctx, a, b, out);
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
rt_view_context_pick_point_bsg(void *ctx, int x, int y, int first_only)
{
    return rt_view_pick_point_bsg((struct bsg_view *)ctx, x, y, first_only);
}

void *
rt_view_context_pick_point(void *ctx, int x, int y, int first_only)
{
    if (_rt_view_context_native_is(ctx))
	return NULL;
    return (void *)rt_view_context_pick_point_bsg(ctx, x, y, first_only);
}

struct rt_view_pick_result_bsg *
rt_view_pick_nearest_bsg(struct bsg_view *v, int x, int y)
{
    return rt_view_pick_result_from_bsg(bsg_pick_nearest(v, x, y));
}

struct rt_view_pick_result_bsg *
rt_view_context_pick_nearest_bsg(void *ctx, int x, int y)
{
    return rt_view_pick_nearest_bsg((struct bsg_view *)ctx, x, y);
}

void *
rt_view_context_pick_nearest(void *ctx, int x, int y)
{
    if (_rt_view_context_native_is(ctx))
	return NULL;
    return (void *)rt_view_context_pick_nearest_bsg(ctx, x, y);
}

struct rt_view_pick_result_bsg *
rt_view_pick_rect_bsg(struct bsg_view *v, int x0, int y0, int x1, int y1)
{
    return rt_view_pick_result_from_bsg(bsg_pick_rect(v, x0, y0, x1, y1));
}

struct rt_view_pick_result_bsg *
rt_view_context_pick_rect_bsg(void *ctx, int x0, int y0, int x1, int y1)
{
    return rt_view_pick_rect_bsg((struct bsg_view *)ctx, x0, y0, x1, y1);
}

void *
rt_view_context_pick_rect(void *ctx, int x0, int y0, int x1, int y1)
{
    if (_rt_view_context_native_is(ctx))
	return NULL;
    return (void *)rt_view_context_pick_rect_bsg(ctx, x0, y0, x1, y1);
}

struct rt_view_pick_result_bsg *
rt_view_pick_semantic_path_bsg(struct bsg_view *v, const char *path_pattern)
{
    return rt_view_pick_result_from_bsg(bsg_pick_semantic_path(v, path_pattern));
}

struct rt_view_pick_result_bsg *
rt_view_context_pick_semantic_path_bsg(void *ctx, const char *path_pattern)
{
    return rt_view_pick_semantic_path_bsg((struct bsg_view *)ctx,
					  path_pattern);
}

void *
rt_view_context_pick_semantic_path(void *ctx, const char *path_pattern)
{
    struct rt_view_context_scene_adapter adapter;
    if (rt_view_context_scene_adapter_get(ctx, &adapter) &&
	adapter.pick_semantic_path) {
	void *result = adapter.pick_semantic_path(ctx, path_pattern,
		       adapter.data);
	if (result)
	    return result;
    }

    if (_rt_view_context_native_is(ctx))
	return NULL;

    return (void *)rt_view_context_pick_semantic_path_bsg(ctx, path_pattern);
}

struct rt_view_pick_result_bsg *
rt_view_pick_result_create_bsg(void)
{
    return rt_view_pick_result_from_bsg(bsg_pick_result_create());
}

void *
rt_view_pick_result_context_create_bsg(void)
{
    return (void *)rt_view_pick_result_create_bsg();
}

void *
rt_view_pick_result_context_create(void)
{
    return rt_view_pick_result_context_create_bsg();
}

void
rt_view_pick_result_free_bsg(struct rt_view_pick_result_bsg *result)
{
    bsg_pick_result_free(rt_view_pick_result_to_bsg(result));
}

void
rt_view_pick_result_context_free_bsg(void *result_ctx)
{
    rt_view_pick_result_free_bsg((struct rt_view_pick_result_bsg *)result_ctx);
}

void
rt_view_pick_result_context_free(void *result_ctx)
{
    rt_view_pick_result_context_free_bsg(result_ctx);
}

size_t
rt_view_pick_result_count_bsg(const struct rt_view_pick_result_bsg *result)
{
    return bsg_pick_result_count(rt_view_pick_result_const_to_bsg(result));
}

size_t
rt_view_pick_result_context_count_bsg(const void *result_ctx)
{
    return rt_view_pick_result_count_bsg(
	       (const struct rt_view_pick_result_bsg *)result_ctx);
}

size_t
rt_view_pick_result_context_count(const void *result_ctx)
{
    return rt_view_pick_result_context_count_bsg(result_ctx);
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

int
rt_view_pick_result_context_path_bsg(const void *result_ctx,
				     size_t index,
				     struct bu_vls *path_out)
{
    return rt_view_pick_result_path_bsg(
	       (const struct rt_view_pick_result_bsg *)result_ctx, index,
	       path_out);
}

int
rt_view_pick_result_context_path(const void *result_ctx,
				 size_t index,
				 struct bu_vls *path_out)
{
    return rt_view_pick_result_context_path_bsg(result_ctx, index, path_out);
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

fastf_t
rt_view_pick_result_context_hit_dist_bsg(const void *result_ctx, size_t index)
{
    return rt_view_pick_result_hit_dist_bsg(
	       (const struct rt_view_pick_result_bsg *)result_ctx, index);
}

fastf_t
rt_view_pick_result_context_hit_dist(const void *result_ctx, size_t index)
{
    return rt_view_pick_result_context_hit_dist_bsg(result_ctx, index);
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
rt_view_pick_result_context_append_path_bsg(void *result_ctx,
	void *view_ctx,
	int screen_x,
	int screen_y,
	const char *source_path,
	fastf_t hit_dist)
{
    return rt_view_pick_result_append_path_bsg(
	       (struct rt_view_pick_result_bsg *)result_ctx,
	       (struct bsg_view *)view_ctx, screen_x, screen_y, source_path,
	       hit_dist);
}

int
rt_view_pick_result_context_append_path(void *result_ctx,
					void *view_ctx,
					int screen_x,
					int screen_y,
					const char *source_path,
					fastf_t hit_dist)
{
    return rt_view_pick_result_context_append_path_bsg(result_ctx, view_ctx,
	    screen_x, screen_y, source_path, hit_dist);
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

int
rt_view_pick_result_context_append_copy_bsg(void *dest_ctx,
	const void *src_ctx,
	size_t index,
	fastf_t hit_dist)
{
    return rt_view_pick_result_append_copy_bsg(
	       (struct rt_view_pick_result_bsg *)dest_ctx,
	       (const struct rt_view_pick_result_bsg *)src_ctx, index, hit_dist);
}

int
rt_view_pick_result_context_append_copy(void *dest_ctx,
					const void *src_ctx,
					size_t index,
					fastf_t hit_dist)
{
    return rt_view_pick_result_context_append_copy_bsg(dest_ctx, src_ctx, index,
	    hit_dist);
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

void *
rt_view_pick_result_context_filter_first_bsg(const void *src_ctx)
{
    return (void *)rt_view_pick_result_filter_first_bsg(
	       (const struct rt_view_pick_result_bsg *)src_ctx);
}

void *
rt_view_pick_result_context_filter_first(const void *src_ctx)
{
    return rt_view_pick_result_context_filter_first_bsg(src_ctx);
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

int
rt_view_context_selection_available_bsg(void *ctx)
{
    return rt_view_selection_available_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_selection_available(void *ctx)
{
    struct rt_view_context_selection_adapter adapter;
    if (rt_view_context_selection_adapter_get(ctx, &adapter) &&
	adapter.available)
	return adapter.available(ctx, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_selection_available_bsg(ctx);
}

size_t
rt_view_selection_count_bsg(struct bsg_view *v)
{
    struct bsg_selection *selection = rt_view_selection_get_bsg(v);
    return selection ? bsg_selection_count(selection) : 0;
}

size_t
rt_view_context_selection_count_bsg(void *ctx)
{
    return rt_view_selection_count_bsg((struct bsg_view *)ctx);
}

size_t
rt_view_context_selection_count(void *ctx)
{
    struct rt_view_context_selection_adapter adapter;
    if (rt_view_context_selection_adapter_get(ctx, &adapter) &&
	adapter.count)
	return adapter.count(ctx, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_selection_count_bsg(ctx);
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
rt_view_context_selection_set_pick_result_ref_bsg(
    void *ctx,
    const struct rt_view_pick_result_bsg *result,
    rt_view_selection_path_callback_bsg_t callback,
    void *data)
{
    return rt_view_selection_set_pick_result_ref_bsg((struct bsg_view *)ctx,
	    result, callback, data);
}

int
rt_view_context_selection_set_pick_result_context_bsg(
    void *ctx,
    const void *result_ctx,
    rt_view_selection_path_callback_bsg_t callback,
    void *data)
{
    return rt_view_context_selection_set_pick_result_ref_bsg(ctx,
	    (const struct rt_view_pick_result_bsg *)result_ctx, callback, data);
}

int
rt_view_context_selection_set_pick_result_context(
    void *ctx,
    const void *result_ctx,
    rt_view_selection_path_callback_t callback,
    void *data)
{
    struct rt_view_context_selection_adapter adapter;
    if (rt_view_context_selection_adapter_get(ctx, &adapter) &&
	adapter.set_pick_result_context)
	return adapter.set_pick_result_context(ctx, result_ctx, callback,
					       data, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_selection_set_pick_result_context_bsg(ctx,
	    result_ctx, (rt_view_selection_path_callback_bsg_t)callback, data);
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

int
rt_view_context_selection_clear_bsg(void *ctx)
{
    return rt_view_selection_clear_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_selection_clear(void *ctx)
{
    struct rt_view_context_selection_adapter adapter;
    if (rt_view_context_selection_adapter_get(ctx, &adapter) &&
	adapter.clear)
	return adapter.clear(ctx, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_selection_clear_bsg(ctx);
}

struct bu_ptbl *
rt_view_set_views_bsg(struct bsg_view_set *s)
{
    return bsg_set_views(s);
}

struct bu_ptbl *
rt_view_set_context_views(void *view_set_ctx)
{
    if (_rt_view_set_context_native_is(view_set_ctx))
	return _rt_view_set_context_native_views(view_set_ctx);
    return rt_view_set_views_bsg((struct bsg_view_set *)view_set_ctx);
}

void *
rt_view_set_recycle_pool_bsg(struct bsg_view_set *s)
{
    return bsg_set_fsos(s);
}

void *
rt_view_set_context_recycle_pool(void *view_set_ctx)
{
    if (_rt_view_set_context_native_is(view_set_ctx))
	return _rt_view_set_context_native_recycle_pool(view_set_ctx);
    return rt_view_set_recycle_pool_bsg((struct bsg_view_set *)view_set_ctx);
}

struct bsg_view *
rt_view_set_find_view_bsg(struct bsg_view_set *s, const char *name)
{
    if (!s || !name)
	return NULL;

    return bsg_set_find_view(s, name);
}

void *
rt_view_set_context_find_view(void *view_set_ctx, const char *name)
{
    if (_rt_view_set_context_native_is(view_set_ctx))
	return _rt_view_set_context_native_find_view(view_set_ctx, name);
    return (void *)rt_view_set_find_view_bsg(
	       (struct bsg_view_set *)view_set_ctx, name);
}

void *
rt_view_set_context_create(void)
{
    return _rt_view_set_context_native_create();
}

void
rt_view_set_context_destroy(void *view_set_ctx)
{
    if (!view_set_ctx)
	return;

    if (_rt_view_set_context_native_is(view_set_ctx)) {
	_rt_view_set_context_native_destroy(view_set_ctx);
	return;
    }

    struct bsg_view_set *s = (struct bsg_view_set *)view_set_ctx;
    rt_view_set_context_free(view_set_ctx);
    BU_PUT(s, struct bsg_view_set);
}

void
rt_view_set_init_bsg(struct bsg_view_set *s)
{
    if (!s)
	return;

    bsg_set_init(s);
}

void
rt_view_set_context_init(void *view_set_ctx)
{
    if (_rt_view_set_context_native_is(view_set_ctx)) {
	_rt_view_set_context_native_init(view_set_ctx);
	return;
    }
    rt_view_set_init_bsg((struct bsg_view_set *)view_set_ctx);
}

void
rt_view_set_free_bsg(struct bsg_view_set *s)
{
    if (!s)
	return;

    bsg_set_free(s);
}

void
rt_view_set_context_free(void *view_set_ctx)
{
    if (_rt_view_set_context_native_is(view_set_ctx)) {
	_rt_view_set_context_native_free(view_set_ctx);
	return;
    }
    rt_view_set_free_bsg((struct bsg_view_set *)view_set_ctx);
}

int
rt_view_context_init_bsg(void *ctx, void *view_set_ctx)
{
    if (!ctx)
	return 0;

    rt_view_init_bsg((struct bsg_view *)ctx,
		     (struct bsg_view_set *)view_set_ctx);
    return 1;
}

int
rt_view_context_free_contents_bsg(void *ctx)
{
    if (!ctx)
	return 0;

    rt_view_free_bsg((struct bsg_view *)ctx);
    return 1;
}

int
rt_view_context_view_set_attach_bsg(void *ctx, void *view_set_ctx)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return 0;

    v->vset = (struct bsg_view_set *)view_set_ctx;
    return 1;
}

int
rt_view_context_view_set_attach(void *ctx, void *view_set_ctx)
{
    if (_rt_view_context_native_is(ctx) ||
	_rt_view_set_context_native_is(view_set_ctx))
	return _rt_view_context_native_view_set_attach(ctx, view_set_ctx);
    return rt_view_context_view_set_attach_bsg(ctx, view_set_ctx);
}

void
rt_view_init_bsg(struct bsg_view *v, struct bsg_view_set *s)
{
    bsg_init(v, s);
}

void
rt_view_free_bsg(struct bsg_view *v)
{
    rt_view_update_callback_bridge_clear_bsg(v);
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
rt_view_context_knobs_reset_bsg(void *ctx, int category)
{
    return rt_view_knobs_reset_bsg((struct bsg_view *)ctx, category);
}

int
rt_view_context_knobs_reset(void *ctx, int category)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_reset(ctx, category);
    return rt_view_context_knobs_reset_bsg(ctx, category);
}

static void
rt_view_knobs_from_bsg(struct rt_view_knobs *dst,
		       const struct bsg_view_knobs *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    VMOVE(dst->rot_m, src->rot_m);
    dst->rot_m_flag = src->rot_m_flag;
    dst->origin_m = src->origin_m;
    dst->rot_m_udata = src->rot_m_udata;

    VMOVE(dst->rot_o, src->rot_o);
    dst->rot_o_flag = src->rot_o_flag;
    dst->origin_o = src->origin_o;
    dst->rot_o_udata = src->rot_o_udata;

    VMOVE(dst->rot_v, src->rot_v);
    dst->rot_v_flag = src->rot_v_flag;
    dst->origin_v = src->origin_v;
    dst->rot_v_udata = src->rot_v_udata;

    dst->sca = src->sca;
    dst->sca_flag = src->sca_flag;
    dst->sca_udata = src->sca_udata;

    VMOVE(dst->tra_m, src->tra_m);
    dst->tra_m_flag = src->tra_m_flag;
    dst->tra_m_udata = src->tra_m_udata;

    VMOVE(dst->tra_v, src->tra_v);
    dst->tra_v_flag = src->tra_v_flag;
    dst->tra_v_udata = src->tra_v_udata;

    VMOVE(dst->rot_m_abs, src->rot_m_abs);
    VMOVE(dst->rot_m_abs_last, src->rot_m_abs_last);
    VMOVE(dst->rot_o_abs, src->rot_o_abs);
    VMOVE(dst->rot_o_abs_last, src->rot_o_abs_last);
    VMOVE(dst->rot_v_abs, src->rot_v_abs);
    VMOVE(dst->rot_v_abs_last, src->rot_v_abs_last);
    dst->sca_abs = src->sca_abs;
    VMOVE(dst->tra_m_abs, src->tra_m_abs);
    VMOVE(dst->tra_m_abs_last, src->tra_m_abs_last);
    VMOVE(dst->tra_v_abs, src->tra_v_abs);
    VMOVE(dst->tra_v_abs_last, src->tra_v_abs_last);
}

static void
rt_view_knobs_to_bsg(struct bsg_view_knobs *dst,
		     const struct rt_view_knobs *src)
{
    memset(dst, 0, sizeof(*dst));
    if (!src)
	return;

    VMOVE(dst->rot_m, src->rot_m);
    dst->rot_m_flag = src->rot_m_flag;
    dst->origin_m = src->origin_m;
    dst->rot_m_udata = src->rot_m_udata;

    VMOVE(dst->rot_o, src->rot_o);
    dst->rot_o_flag = src->rot_o_flag;
    dst->origin_o = src->origin_o;
    dst->rot_o_udata = src->rot_o_udata;

    VMOVE(dst->rot_v, src->rot_v);
    dst->rot_v_flag = src->rot_v_flag;
    dst->origin_v = src->origin_v;
    dst->rot_v_udata = src->rot_v_udata;

    dst->sca = src->sca;
    dst->sca_flag = src->sca_flag;
    dst->sca_udata = src->sca_udata;

    VMOVE(dst->tra_m, src->tra_m);
    dst->tra_m_flag = src->tra_m_flag;
    dst->tra_m_udata = src->tra_m_udata;

    VMOVE(dst->tra_v, src->tra_v);
    dst->tra_v_flag = src->tra_v_flag;
    dst->tra_v_udata = src->tra_v_udata;

    VMOVE(dst->rot_m_abs, src->rot_m_abs);
    VMOVE(dst->rot_m_abs_last, src->rot_m_abs_last);
    VMOVE(dst->rot_o_abs, src->rot_o_abs);
    VMOVE(dst->rot_o_abs_last, src->rot_o_abs_last);
    VMOVE(dst->rot_v_abs, src->rot_v_abs);
    VMOVE(dst->rot_v_abs_last, src->rot_v_abs_last);
    dst->sca_abs = src->sca_abs;
    VMOVE(dst->tra_m_abs, src->tra_m_abs);
    VMOVE(dst->tra_m_abs_last, src->tra_m_abs_last);
    VMOVE(dst->tra_v_abs, src->tra_v_abs);
    VMOVE(dst->tra_v_abs_last, src->tra_v_abs_last);
}

int
rt_view_context_knobs_state_from_bsg(struct rt_view_knobs *knobs,
				     const void *ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)ctx;

    if (!knobs)
	return 0;

    rt_view_knobs_from_bsg(knobs, v ? &v->k : NULL);
    return v ? 1 : 0;
}

int
rt_view_context_knobs_state_get(struct rt_view_knobs *knobs,
				const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_state_get(knobs, ctx);
    return rt_view_context_knobs_state_from_bsg(knobs, ctx);
}

int
rt_view_context_knobs_state_set_bsg(void *ctx,
				    const struct rt_view_knobs *knobs)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v || !knobs)
	return 0;

    rt_view_knobs_to_bsg(&v->k, knobs);
    return 1;
}

int
rt_view_context_knobs_state_set(void *ctx,
				const struct rt_view_knobs *knobs)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_state_set(ctx, knobs);
    return rt_view_context_knobs_state_set_bsg(ctx, knobs);
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

unsigned long long
rt_view_context_knobs_hash_bsg(void *ctx,
			       struct bu_data_hash_state *state)
{
    return rt_view_knobs_hash_bsg((struct bsg_view *)ctx, state);
}

unsigned long long
rt_view_context_knobs_hash(void *ctx,
			   struct bu_data_hash_state *state)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_hash(ctx, state);
    return rt_view_context_knobs_hash_bsg(ctx, state);
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
rt_view_context_knobs_cmd_process_bsg(vect_t *rvec,
				      int *do_rot,
				      vect_t *tvec,
				      int *do_tran,
				      void *ctx,
				      const char *cmd,
				      fastf_t factor,
				      char origin,
				      int model_flag,
				      int incr_flag)
{
    return rt_view_knobs_cmd_process_bsg(rvec, do_rot, tvec, do_tran,
					 (struct bsg_view *)ctx, cmd, factor, origin, model_flag, incr_flag);
}

int
rt_view_context_knobs_cmd_process(vect_t *rvec,
				  int *do_rot,
				  vect_t *tvec,
				  int *do_tran,
				  void *ctx,
				  const char *cmd,
				  fastf_t factor,
				  char origin,
				  int model_flag,
				  int incr_flag)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_cmd_process(rvec, do_rot,
		tvec, do_tran, ctx, cmd, factor, origin, model_flag,
		incr_flag);
    return rt_view_context_knobs_cmd_process_bsg(rvec, do_rot, tvec, do_tran,
	    ctx, cmd, factor, origin, model_flag, incr_flag);
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
rt_view_context_knobs_translate_bsg(void *ctx,
				    const vect_t tvec,
				    int model_flag)
{
    return rt_view_knobs_translate_bsg((struct bsg_view *)ctx, tvec,
				       model_flag);
}

int
rt_view_context_knobs_translate(void *ctx,
				const vect_t tvec,
				int model_flag)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_translate(ctx, tvec,
		model_flag);
    return rt_view_context_knobs_translate_bsg(ctx, tvec, model_flag);
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
rt_view_context_knobs_rotate_bsg(void *ctx,
				 const vect_t rvec,
				 char origin,
				 char coords,
				 const matp_t obj_rot,
				 const pointp_t pvt_pt)
{
    return rt_view_knobs_rotate_bsg((struct bsg_view *)ctx, rvec, origin,
				    coords, obj_rot, pvt_pt);
}

int
rt_view_context_knobs_rotate(void *ctx,
			     const vect_t rvec,
			     char origin,
			     char coords,
			     const matp_t obj_rot,
			     const pointp_t pvt_pt)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_rotate(ctx, rvec, origin,
		coords, obj_rot, pvt_pt);
    return rt_view_context_knobs_rotate_bsg(ctx, rvec, origin, coords,
					    obj_rot, pvt_pt);
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
rt_view_context_knobs_update_rate_flags_bsg(void *ctx)
{
    return rt_view_knobs_update_rate_flags_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_knobs_update_rate_flags(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_update_rate_flags(ctx);
    return rt_view_context_knobs_update_rate_flags_bsg(ctx);
}

int
rt_view_context_knob_values_from_bsg(struct rt_view_knob_values *values,
				     const void *ctx)
{
    const struct bsg_view *v = (const struct bsg_view *)ctx;

    if (!values)
	return 0;

    memset(values, 0, sizeof(*values));
    if (!v)
	return 0;

    VMOVE(values->rate_rotation, v->k.rot_v);
    VMOVE(values->rate_translation, v->k.tra_v);
    values->rate_scale = v->k.sca;
    VMOVE(values->absolute_rotation, v->k.rot_v_abs);
    VMOVE(values->absolute_translation, v->k.tra_v_abs);
    values->absolute_scale = rt_view_absolute_scale_from_bsg(v);
    return 1;
}

int
rt_view_context_knob_values_get(struct rt_view_knob_values *values,
				const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knob_values_get(values, ctx);
    return rt_view_context_knob_values_from_bsg(values, ctx);
}

int
rt_view_context_knobs_calibrate_bsg(void *ctx)
{
    struct bsg_view *v = (struct bsg_view *)ctx;

    if (!v)
	return 0;

    VSETALL(v->k.tra_v_abs, 0.0);
    VSETALL(v->k.tra_v_abs_last, 0.0);
    VSETALL(v->k.tra_m_abs, 0.0);
    VSETALL(v->k.tra_m_abs_last, 0.0);
    return 1;
}

int
rt_view_context_knobs_calibrate(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_knobs_calibrate(ctx);
    return rt_view_context_knobs_calibrate_bsg(ctx);
}

int
rt_view_is_independent_bsg(const struct bsg_view *v)
{
    return bsg_view_is_independent(v);
}

rt_view_scene_ref
rt_view_scene_ref_null(void)
{
    rt_view_scene_ref ref = RT_VIEW_SCENE_REF_NULL_INIT;
    return ref;
}

rt_view_scene_ref
rt_view_scene_ref_make(void *opaque, unsigned int backend)
{
    rt_view_scene_ref ref = RT_VIEW_SCENE_REF_NULL_INIT;
    ref.opaque = opaque;
    ref.backend = opaque ? backend : RT_VIEW_SCENE_BACKEND_NONE;
    return ref;
}

int
rt_view_scene_ref_is_null(rt_view_scene_ref ref)
{
    return ref.opaque ? 0 : 1;
}

int
rt_view_scene_ref_equal(rt_view_scene_ref a, rt_view_scene_ref b)
{
    if (rt_view_scene_ref_is_null(a) || rt_view_scene_ref_is_null(b))
	return rt_view_scene_ref_is_null(a) && rt_view_scene_ref_is_null(b);
    if (a.backend != RT_VIEW_SCENE_BACKEND_NONE &&
	b.backend != RT_VIEW_SCENE_BACKEND_NONE &&
	a.backend != b.backend)
	return 0;
    return a.opaque == b.opaque;
}

void *
rt_view_scene_ref_context(rt_view_scene_ref ref)
{
    return ref.opaque;
}

unsigned int
rt_view_scene_ref_backend(rt_view_scene_ref ref)
{
    return ref.opaque ? ref.backend : RT_VIEW_SCENE_BACKEND_NONE;
}

rt_view_scene_ref_bsg
rt_view_scene_ref_null_bsg(void)
{
    return rt_view_scene_ref_from_bsg(bsg_scene_ref_null());
}

int
rt_view_scene_ref_is_null_bsg(rt_view_scene_ref_bsg ref)
{
    return bsg_scene_ref_is_null(rt_view_scene_ref_to_bsg(ref));
}

int
rt_view_scene_ref_equal_bsg(rt_view_scene_ref_bsg a, rt_view_scene_ref_bsg b)
{
    return bsg_scene_ref_equal(rt_view_scene_ref_to_bsg(a),
			       rt_view_scene_ref_to_bsg(b));
}

rt_view_scene_ref_bsg
rt_view_independent_scope_ref_bsg(struct bsg_view *v, int create)
{
    return rt_view_scene_ref_from_bsg(bsg_view_independent_scope_ref(v, create));
}

rt_view_scene_ref
rt_view_context_independent_scope_ref(void *ctx, int create)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_independent_scope_ref(ctx, create);
    return rt_view_neutral_scene_ref_from_bsg(
	       bsg_view_independent_scope_ref((struct bsg_view *)ctx, create));
}


rt_view_scene_ref
rt_view_context_scene_root_ref(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scene_root_ref(ctx);
    return rt_view_neutral_scene_ref_from_bsg(
	       bsg_view_scene_ref((const struct bsg_view *)ctx));
}


int
rt_view_context_scene_root_ref_attach(void *ctx, rt_view_scene_ref root_ref)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scene_root_ref_attach(ctx, root_ref);
    if (!rt_view_scene_ref_is_null(root_ref) &&
	rt_view_scene_ref_backend(root_ref) != RT_VIEW_SCENE_BACKEND_NONE &&
	rt_view_scene_ref_backend(root_ref) != RT_VIEW_SCENE_BACKEND_BSG)
	return 0;
    return bsg_view_scene_ref_attach((struct bsg_view *)ctx,
				     rt_view_neutral_scene_ref_to_bsg(root_ref));
}


int
rt_view_independent_scope_is_null_bsg(struct bsg_view *v, int create)
{
    return rt_view_scene_ref_is_null_bsg(
	       rt_view_independent_scope_ref_bsg(v, create));
}

void
rt_view_independent_scope_destroy_bsg(struct bsg_view *v)
{
    bsg_view_independent_scope_destroy(v);
}

int
rt_view_scene_attached_bsg(const struct bsg_view *v)
{
    if (!v)
	return 0;

    return bsg_view_scene_attached(v);
}

int
rt_view_context_scene_attached_bsg(const void *ctx)
{
    return rt_view_scene_attached_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_scene_attached(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scene_attached(ctx);
    return rt_view_context_scene_attached_bsg(ctx);
}

int
rt_view_scene_anchor_ensure_bsg(struct bsg_view *v)
{
    bsg_separator_ref root;

    if (!v)
	return 0;

    root = bsg_view_scene_separator_ref(v, 1);
    return !bsg_node_ref_is_null(bsg_separator_ref_as_node(root));
}

int
rt_view_context_scene_anchor_ensure_bsg(void *ctx)
{
    return rt_view_scene_anchor_ensure_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_scene_anchor_ensure(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return 0;
    return rt_view_context_scene_anchor_ensure_bsg(ctx);
}

int
rt_view_scene_shared_bsg(const struct bsg_view *a, const struct bsg_view *b)
{
    if (!a || !b)
	return 0;

    return bsg_view_scene_shared(a, b);
}

int
rt_view_context_scene_shared_bsg(const void *a_ctx, const void *b_ctx)
{
    return rt_view_scene_shared_bsg((const struct bsg_view *)a_ctx,
				    (const struct bsg_view *)b_ctx);
}

int
rt_view_context_scene_shared(const void *a_ctx, const void *b_ctx)
{
    if (_rt_view_context_native_is(a_ctx) || _rt_view_context_native_is(b_ctx))
	return 0;
    return rt_view_context_scene_shared_bsg(a_ctx, b_ctx);
}

int
rt_view_visible_render_summary_bsg(struct bsg_view *v,
				   struct rt_view_render_summary_bsg *summary)
{
    struct bsg_render_request *req = NULL;
    struct bsg_render_batch *batch = NULL;
    int ret = 0;
    int nitems = 0;

    if (summary) {
	summary->item_count = 0;
	summary->highlighted_count = 0;
    }

    if (!v || !summary)
	return 0;

    req = bsg_render_request_create(v, NULL);
    batch = bsg_render_batch_create();
    if (!req || !batch)
	goto cleanup;

    bsg_render_request_set_flags(req, BSG_RENDER_FLAG_VISIBLE_ONLY);
    nitems = bsg_render_request_collect(req, batch);
    if (nitems < 0)
	goto cleanup;

    summary->item_count = nitems;
    for (size_t i = 0; i < bsg_render_batch_count(batch); i++) {
	const struct bsg_render_item *item = bsg_render_batch_get(batch, i);
	if (item && item->appearance.highlighted)
	    summary->highlighted_count++;
    }

    ret = 1;

cleanup:
    bsg_render_batch_destroy(batch);
    bsg_render_request_destroy(req);
    return ret;
}

int
rt_view_context_visible_render_summary_bsg(
    void *ctx,
    struct rt_view_render_summary_bsg *summary)
{
    return rt_view_visible_render_summary_bsg((struct bsg_view *)ctx, summary);
}

int
rt_view_context_visible_render_summary(void *ctx,
				       struct rt_view_render_summary *summary)
{
    struct rt_view_render_summary_bsg bsg_summary =
	    RT_VIEW_RENDER_SUMMARY_BSG_INIT;

    if (summary) {
	summary->item_count = 0;
	summary->highlighted_count = 0;
    }

    if (_rt_view_context_native_is(ctx))
	return 0;

    if (!summary ||
	!rt_view_context_visible_render_summary_bsg(ctx, &bsg_summary))
	return 0;

    summary->item_count = bsg_summary.item_count;
    summary->highlighted_count = bsg_summary.highlighted_count;
    return 1;
}

int
rt_view_context_named_line_render_count_bsg(void *ctx, const char *name)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    struct bsg_render_request *req = NULL;
    struct bsg_render_batch *batch = NULL;
    int count = 0;

    if (!v || !name || !bsg_view_scene_attached(v))
	return 0;

    req = bsg_render_request_create(v, NULL);
    batch = bsg_render_batch_create();
    if (!req || !batch)
	goto cleanup;

    bsg_render_request_set_flags(req,
				 BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE);
    if (bsg_render_request_collect(req, batch) < 0)
	goto cleanup;

    for (size_t i = 0; i < bsg_render_batch_count(batch); i++) {
	const struct bsg_render_item *item = bsg_render_batch_get(batch, i);
	if (item && item->name && BU_STR_EQUAL(item->name, name) &&
	    item->geometry.kind == BSG_RENDER_GEOMETRY_LINE_SET)
	    count++;
    }

cleanup:
    bsg_render_batch_destroy(batch);
    bsg_render_request_destroy(req);
    return count;
}

int
rt_view_context_named_line_render_count(void *ctx, const char *name)
{
    if (_rt_view_context_native_is(ctx))
	return 0;
    return rt_view_context_named_line_render_count_bsg(ctx, name);
}

static int
rt_view_path_matches_drawn_prefix_bsg(const char *path, const char *drawn_prefix)
{
    size_t n;

    if (!path || !drawn_prefix)
	return 0;

    n = strlen(drawn_prefix);
    return BU_STR_EQUAL(path, drawn_prefix) ||
	   (strlen(path) > n && bu_strncmp(path, drawn_prefix, n) == 0 && path[n] == '/');
}

static const struct bsg_export_record *
rt_view_find_export_record_bsg(const struct bsg_export_result *result,
			       const char *drawn_prefix)
{
    if (!result || !drawn_prefix)
	return NULL;

    for (size_t i = 0; i < bsg_export_result_count(result); i++) {
	const struct bsg_export_record *rec = bsg_export_result_get(result, i);
	if (rec && rt_view_path_matches_drawn_prefix_bsg(
		bu_vls_cstr(&rec->path), drawn_prefix))
	    return rec;
    }

    return NULL;
}

static const struct bsg_render_item *
rt_view_find_render_item_bsg(const struct bsg_render_batch *batch,
			     const char *drawn_prefix)
{
    if (!batch || !drawn_prefix)
	return NULL;

    for (size_t i = 0; i < bsg_render_batch_count(batch); i++) {
	const struct bsg_render_item *item = bsg_render_batch_get(batch, i);
	if (item && item->source.name &&
	    rt_view_path_matches_drawn_prefix_bsg(item->source.name,
		    drawn_prefix))
	    return item;
    }

    return NULL;
}

int
rt_view_render_export_consistency_bsg(struct bsg_view *v,
				      const char *drawn_prefix,
				      struct rt_view_render_export_consistency_bsg *summary)
{
    struct bsg_export_request ereq;
    struct bsg_export_result *export_result = NULL;
    struct bsg_render_request *req = NULL;
    struct bsg_render_batch *batch = NULL;
    struct bsg_backend_scene *scene = NULL;
    const struct bsg_export_record *export_rec = NULL;
    const struct bsg_render_item *item = NULL;
    int ret = 0;

    if (summary) {
	summary->export_record_found = 0;
	summary->render_item_found = 0;
	summary->backend_node_found = 0;
	summary->export_render_consistent = 0;
	summary->export_backend_consistent = 0;
    }

    if (!v || !drawn_prefix || !summary)
	return 0;

    bsg_export_request_init(&ereq, v);
    ereq.query_flags = BSG_EXPORT_QUERY_VISIBLE_ONLY | BSG_EXPORT_QUERY_DB_OBJECTS;
    ereq.render_flags = BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE;
    export_result = bsg_export_query(&ereq);
    if (!export_result)
	goto cleanup;

    export_rec = rt_view_find_export_record_bsg(export_result, drawn_prefix);
    summary->export_record_found = export_rec ? 1 : 0;

    req = bsg_render_request_create(v, NULL);
    batch = bsg_render_batch_create();
    if (!req || !batch)
	goto cleanup;

    bsg_render_request_set_flags(req,
				 BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE);
    if (bsg_render_request_collect(req, batch) < 0)
	goto cleanup;

    item = rt_view_find_render_item_bsg(batch, drawn_prefix);
    summary->render_item_found = item ? 1 : 0;
    if (export_rec && item) {
	summary->export_render_consistent =
	    export_rec->cache_identity == item->cache_identity &&
	    export_rec->source.source_id == item->source.source_id &&
	    export_rec->geometry_revision == item->geometry.revision &&
	    export_rec->payload_revision == item->payload_revision &&
	    export_rec->draw_mode == item->appearance.draw_mode &&
	    export_rec->visible == item->visible;
    }

    scene = bsg_backend_scene_create();
    if (!scene)
	goto cleanup;
    if (bsg_backend_scene_render_request(v, scene,
					 BSG_RENDER_FLAG_VISIBLE_ONLY | BSG_RENDER_FLAG_PAYLOAD_PREPARE) < 0)
	goto cleanup;

    if (export_rec) {
	const struct bsg_backend_scene_node *node =
	    bsg_backend_scene_find(scene, export_rec->cache_identity);
	summary->backend_node_found = node ? 1 : 0;
	if (node) {
	    summary->export_backend_consistent =
		node->source_identity == export_rec->source.source_id &&
		node->geometry.revision == export_rec->geometry_revision &&
		node->material.draw_mode == export_rec->draw_mode &&
		node->selection.visible == export_rec->visible;
	}
    }

    ret = 1;

cleanup:
    if (scene)
	bsg_backend_scene_destroy(scene);
    bsg_render_batch_destroy(batch);
    bsg_render_request_destroy(req);
    bsg_export_result_free(export_result);
    return ret;
}

int
rt_view_context_render_export_consistency_bsg(
    void *ctx,
    const char *drawn_prefix,
    struct rt_view_render_export_consistency_bsg *summary)
{
    return rt_view_render_export_consistency_bsg((struct bsg_view *)ctx,
	    drawn_prefix, summary);
}

int
rt_view_context_render_export_consistency(
    void *ctx,
    const char *drawn_prefix,
    struct rt_view_render_export_consistency *summary)
{
    struct rt_view_render_export_consistency_bsg bsg_summary =
	    RT_VIEW_RENDER_EXPORT_CONSISTENCY_BSG_INIT;

    if (summary) {
	summary->export_record_found = 0;
	summary->render_item_found = 0;
	summary->backend_node_found = 0;
	summary->export_render_consistent = 0;
	summary->export_backend_consistent = 0;
    }

    struct rt_view_context_scene_adapter adapter;
    if (rt_view_context_scene_adapter_get(ctx, &adapter) &&
	adapter.render_export_consistency &&
	adapter.render_export_consistency(ctx, drawn_prefix, summary,
					  adapter.data))
	return 1;

    if (_rt_view_context_native_is(ctx))
	return 0;

    if (!summary ||
	!rt_view_context_render_export_consistency_bsg(ctx, drawn_prefix,
		&bsg_summary))
	return 0;

    summary->export_record_found = bsg_summary.export_record_found;
    summary->render_item_found = bsg_summary.render_item_found;
    summary->backend_node_found = bsg_summary.backend_node_found;
    summary->export_render_consistent = bsg_summary.export_render_consistent;
    summary->export_backend_consistent = bsg_summary.export_backend_consistent;
    return 1;
}

void
rt_view_init_copy_bsg(struct bsg_view *dest,
		      const struct bsg_view *src,
		      struct bsg_view_set *s)
{
    bsg_view_init_copy(dest, src, s);
}

int
rt_view_context_init_copy_bsg(void *dest_ctx,
			      const void *src_ctx,
			      void *view_set_ctx)
{
    if (!dest_ctx)
	return 0;

    rt_view_init_copy_bsg((struct bsg_view *)dest_ctx,
			  (const struct bsg_view *)src_ctx,
			  (struct bsg_view_set *)view_set_ctx);
    return 1;
}

size_t
rt_view_clear_bsg(struct bsg_view *v, int flags)
{
    return bsg_clear(v, flags);
}

size_t
rt_view_context_clear_bsg(void *ctx, int flags)
{
    return rt_view_clear_bsg((struct bsg_view *)ctx, flags);
}

static int
rt_view_clear_flags_to_bsg(int flags)
{
    int bsg_flags = 0;

    if (flags & RT_VIEW_CLEAR_DB)
	bsg_flags |= RT_VIEW_CLEAR_DB_BSG;
    if (flags & RT_VIEW_CLEAR_VIEW)
	bsg_flags |= RT_VIEW_CLEAR_VIEW_BSG;
    if (flags & RT_VIEW_CLEAR_LOCAL)
	bsg_flags |= RT_VIEW_CLEAR_LOCAL_BSG;

    return bsg_flags;
}

size_t
rt_view_context_clear(void *ctx, int flags)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_clear(ctx, flags);
    return rt_view_context_clear_bsg(ctx, rt_view_clear_flags_to_bsg(flags));
}

fastf_t
rt_view_perspective_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_perspective : 0.0;
}

fastf_t
rt_view_context_perspective_from_bsg(const void *ctx)
{
    return rt_view_perspective_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_perspective_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_perspective_get(ctx);
    return rt_view_context_perspective_from_bsg(ctx);
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
rt_view_context_perspective_set_bsg(void *ctx, fastf_t perspective)
{
    return rt_view_perspective_set_bsg((struct bsg_view *)ctx, perspective);
}

int
rt_view_context_perspective_set(void *ctx, fastf_t perspective)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_perspective_set(ctx, perspective);
    return rt_view_context_perspective_set_bsg(ctx, perspective);
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
rt_view_context_model2view_from_bsg(mat_t model2view, const void *ctx)
{
    return rt_view_model2view_from_bsg(model2view, (const struct bsg_view *)ctx);
}

int
rt_view_context_model2view_get(mat_t model2view, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_model2view_get(model2view, ctx);
    return rt_view_context_model2view_from_bsg(model2view, ctx);
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
rt_view_context_model2view_set_bsg(void *ctx, const mat_t model2view)
{
    return rt_view_model2view_set_bsg((struct bsg_view *)ctx, model2view);
}

int
rt_view_context_model2view_set(void *ctx, const mat_t model2view)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_model2view_set(ctx, model2view);
    return rt_view_context_model2view_set_bsg(ctx, model2view);
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
rt_view_context_view2model_from_bsg(mat_t view2model, const void *ctx)
{
    return rt_view_view2model_from_bsg(view2model,
				       (const struct bsg_view *)ctx);
}

int
rt_view_context_view2model_get(mat_t view2model, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_view2model_get(view2model, ctx);
    return rt_view_context_view2model_from_bsg(view2model, ctx);
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
rt_view_context_view2model_set_bsg(void *ctx, const mat_t view2model)
{
    return rt_view_view2model_set_bsg((struct bsg_view *)ctx, view2model);
}

int
rt_view_context_view2model_set(void *ctx, const mat_t view2model)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_view2model_set(ctx, view2model);
    return rt_view_context_view2model_set_bsg(ctx, view2model);
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
rt_view_context_pmodel2view_from_bsg(mat_t pmodel2view, const void *ctx)
{
    return rt_view_pmodel2view_from_bsg(pmodel2view,
					(const struct bsg_view *)ctx);
}

int
rt_view_context_pmodel2view_get(mat_t pmodel2view, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_pmodel2view_get(pmodel2view, ctx);
    return rt_view_context_pmodel2view_from_bsg(pmodel2view, ctx);
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
rt_view_context_pmodel2view_set_bsg(void *ctx, const mat_t pmodel2view)
{
    return rt_view_pmodel2view_set_bsg((struct bsg_view *)ctx, pmodel2view);
}

int
rt_view_context_pmodel2view_set(void *ctx, const mat_t pmodel2view)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_pmodel2view_set(ctx, pmodel2view);
    return rt_view_context_pmodel2view_set_bsg(ctx, pmodel2view);
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
rt_view_context_pmat_from_bsg(mat_t pmat, const void *ctx)
{
    return rt_view_pmat_from_bsg(pmat, (const struct bsg_view *)ctx);
}

int
rt_view_context_pmat_get(mat_t pmat, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_pmat_get(pmat, ctx);
    return rt_view_context_pmat_from_bsg(pmat, ctx);
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
rt_view_context_pmat_set_bsg(void *ctx, const mat_t pmat)
{
    return rt_view_pmat_set_bsg((struct bsg_view *)ctx, pmat);
}

int
rt_view_context_pmat_set(void *ctx, const mat_t pmat)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_pmat_set(ctx, pmat);
    return rt_view_context_pmat_set_bsg(ctx, pmat);
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
rt_view_context_rotation_from_bsg(mat_t rotation, const void *ctx)
{
    return rt_view_rotation_from_bsg(rotation, (const struct bsg_view *)ctx);
}

int
rt_view_context_rotation_get(mat_t rotation, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_rotation_get(rotation, ctx);
    return rt_view_context_rotation_from_bsg(rotation, ctx);
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
rt_view_context_rotation_set_bsg(void *ctx, const mat_t rotation)
{
    return rt_view_rotation_set_bsg((struct bsg_view *)ctx, rotation);
}

int
rt_view_context_rotation_set(void *ctx, const mat_t rotation)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_rotation_set(ctx, rotation);
    return rt_view_context_rotation_set_bsg(ctx, rotation);
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
rt_view_context_center_from_bsg(mat_t center, const void *ctx)
{
    return rt_view_center_from_bsg(center, (const struct bsg_view *)ctx);
}

int
rt_view_context_center_get(mat_t center, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_center_get(center, ctx);
    return rt_view_context_center_from_bsg(center, ctx);
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
rt_view_context_center_vec_set_bsg(void *ctx, const point_t center)
{
    return rt_view_center_vec_set_bsg((struct bsg_view *)ctx, center);
}

int
rt_view_context_center_set(void *ctx, const point_t center)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_center_set(ctx, center);
    return rt_view_context_center_vec_set_bsg(ctx, center);
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
rt_view_context_plane_from_bsg(plane_t *plane, const void *ctx)
{
    return rt_view_plane_from_bsg(plane, (const struct bsg_view *)ctx);
}

int
rt_view_context_plane_get(plane_t *plane, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_plane_get(plane, ctx);
    return rt_view_context_plane_from_bsg(plane, ctx);
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
rt_view_context_lod_policy_from_bsg(struct rt_view_lod_policy *policy,
				    const void *ctx)
{
    return rt_view_lod_policy_from_bsg(policy, (const struct bsg_view *)ctx);
}

int
rt_view_context_lod_policy_get(struct rt_view_lod_policy *policy,
			       const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_lod_policy_get(policy, ctx);
    return rt_view_context_lod_policy_from_bsg(policy, ctx);
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
rt_view_context_lod_policy_apply_bsg(void *ctx,
				     const struct rt_view_lod_policy *policy)
{
    return rt_view_lod_policy_apply_bsg((struct bsg_view *)ctx, policy);
}

int
rt_view_context_lod_policy_apply(void *ctx,
				 const struct rt_view_lod_policy *policy)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_lod_policy_apply(ctx, policy);
    return rt_view_context_lod_policy_apply_bsg(ctx, policy);
}

int
rt_view_lod_policy_copy_bsg(struct bsg_view *dst, const struct bsg_view *src)
{
    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    if (!dst || !rt_view_lod_policy_from_bsg(&policy, src))
	return 0;
    return rt_view_lod_policy_apply_bsg(dst, &policy);
}

int
rt_view_context_lod_policy_copy_bsg(void *dst_ctx, const void *src_ctx)
{
    return rt_view_lod_policy_copy_bsg((struct bsg_view *)dst_ctx,
				       (const struct bsg_view *)src_ctx);
}

int
rt_view_context_lod_policy_copy(void *dst_ctx, const void *src_ctx)
{
    if (_rt_view_context_native_is(dst_ctx) ||
	_rt_view_context_native_is(src_ctx)) {
	struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;

	if (_rt_view_context_native_is(src_ctx)) {
	    if (!_rt_view_context_native_lod_policy_get(&policy, src_ctx))
		return 0;
	} else {
	    if (!rt_view_context_lod_policy_from_bsg(&policy, src_ctx))
		return 0;
	}

	if (_rt_view_context_native_is(dst_ctx))
	    return _rt_view_context_native_lod_policy_apply(dst_ctx,
		    &policy);
	return rt_view_context_lod_policy_apply_bsg(dst_ctx, &policy);
    }
    return rt_view_context_lod_policy_copy_bsg(dst_ctx, src_ctx);
}

void
rt_view_lod_bounds_update_bsg(struct bsg_view *v)
{
    if (v)
	bsg_view_bounds(v);
}

int
rt_view_context_lod_bounds_update_bsg(void *ctx)
{
    if (!ctx)
	return 0;

    rt_view_lod_bounds_update_bsg((struct bsg_view *)ctx);
    return 1;
}

int
rt_view_context_lod_bounds_update(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_lod_bounds_update(ctx);
    return rt_view_context_lod_bounds_update_bsg(ctx);
}

void
rt_view_lod_bounds_callback_set_bsg(struct bsg_view *v)
{
    if (v)
	v->gv_bounds_update = &bsg_view_bounds;
}

int
rt_view_context_lod_bounds_callback_set_bsg(void *ctx)
{
    if (!ctx)
	return 0;

    rt_view_lod_bounds_callback_set_bsg((struct bsg_view *)ctx);
    return 1;
}

int
rt_view_context_lod_bounds_callback_set(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_lod_bounds_callback_set(ctx);
    return rt_view_context_lod_bounds_callback_set_bsg(ctx);
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

rt_view_bounds_update_callback_bsg_t
rt_view_context_bounds_update_callback_from_bsg(const void *ctx)
{
    return rt_view_bounds_update_callback_from_bsg(
	       (const struct bsg_view *)ctx);
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
rt_view_context_bounds_update_callback_set_bsg(
    void *ctx,
    rt_view_bounds_update_callback_bsg_t callback)
{
    return rt_view_bounds_update_callback_set_bsg((struct bsg_view *)ctx,
	    callback);
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

int
rt_view_context_bounds_update_callback_call_bsg(void *ctx)
{
    return rt_view_bounds_update_callback_call_bsg((struct bsg_view *)ctx);
}

struct rt_view_bounds_update_state {
    rt_view_bounds_update_callback_bsg_t callback;
};

void *
rt_view_context_bounds_update_suspend(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_bounds_update_suspend(ctx);

    rt_view_bounds_update_callback_bsg_t callback =
	rt_view_context_bounds_update_callback_from_bsg(ctx);
    if (!callback)
	return NULL;

    struct rt_view_bounds_update_state *state;
    BU_GET(state, struct rt_view_bounds_update_state);
    state->callback = callback;
    rt_view_context_bounds_update_callback_set_bsg(ctx, NULL);
    return (void *)state;
}

int
rt_view_context_bounds_update_restore(void *ctx, void *state_ctx,
				      int refresh_bounds)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_bounds_update_restore(ctx, state_ctx,
		refresh_bounds);

    struct rt_view_bounds_update_state *state =
	(struct rt_view_bounds_update_state *)state_ctx;
    if (!state)
	return 0;

    int restored = 0;
    if (ctx) {
	restored = rt_view_context_bounds_update_callback_set_bsg(ctx,
		   state->callback);
	if (restored && refresh_bounds)
	    rt_view_context_bounds_update_callback_call_bsg(ctx);
    }

    BU_PUT(state, struct rt_view_bounds_update_state);
    return restored;
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

int
rt_view_context_scale_set_bsg(void *ctx, fastf_t scale)
{
    return rt_view_scale_set_bsg((struct bsg_view *)ctx, scale);
}

int
rt_view_context_scale_set(void *ctx, fastf_t scale)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_scale_set(ctx, scale);
    return rt_view_context_scale_set_bsg(ctx, scale);
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

int
rt_view_context_size_set_bsg(void *ctx, fastf_t size)
{
    return rt_view_size_set_bsg((struct bsg_view *)ctx, size);
}

int
rt_view_context_size_set(void *ctx, fastf_t size)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_size_set(ctx, size);
    return rt_view_context_size_set_bsg(ctx, size);
}

fastf_t
rt_view_inverse_size_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_isize : 1.0;
}

fastf_t
rt_view_context_inverse_size_from_bsg(const void *ctx)
{
    return rt_view_inverse_size_from_bsg((const struct bsg_view *)ctx);
}

fastf_t
rt_view_context_inverse_size_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_inverse_size_get(ctx);
    return rt_view_context_inverse_size_from_bsg(ctx);
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
rt_view_context_eye_pos_from_bsg(point_t eye_pos, const void *ctx)
{
    return rt_view_eye_pos_from_bsg(eye_pos, (const struct bsg_view *)ctx);
}

int
rt_view_context_eye_pos_get(point_t eye_pos, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_eye_pos_get(eye_pos, ctx);
    return rt_view_context_eye_pos_from_bsg(eye_pos, ctx);
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
rt_view_context_eye_pos_set_bsg(void *ctx, const point_t eye_pos)
{
    return rt_view_eye_pos_set_bsg((struct bsg_view *)ctx, eye_pos);
}

int
rt_view_context_eye_pos_set(void *ctx, const point_t eye_pos)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_eye_pos_set(ctx, eye_pos);
    return rt_view_context_eye_pos_set_bsg(ctx, eye_pos);
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
rt_view_context_keypoint_from_bsg(point_t keypoint, const void *ctx)
{
    return rt_view_keypoint_from_bsg(keypoint, (const struct bsg_view *)ctx);
}

int
rt_view_context_keypoint_get(point_t keypoint, const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_keypoint_get(keypoint, ctx);
    return rt_view_context_keypoint_from_bsg(keypoint, ctx);
}

int
rt_view_keypoint_set_bsg(struct bsg_view *v, const point_t keypoint)
{
    if (!v || !keypoint)
	return 0;

    VMOVE(v->gv_keypoint, keypoint);
    return 1;
}

int
rt_view_context_keypoint_set_bsg(void *ctx, const point_t keypoint)
{
    return rt_view_keypoint_set_bsg((struct bsg_view *)ctx, keypoint);
}

int
rt_view_context_keypoint_set(void *ctx, const point_t keypoint)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_keypoint_set(ctx, keypoint);
    return rt_view_context_keypoint_set_bsg(ctx, keypoint);
}

char
rt_view_rotate_about_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_rotate_about : 'v';
}

char
rt_view_context_rotate_about_from_bsg(const void *ctx)
{
    return rt_view_rotate_about_from_bsg((const struct bsg_view *)ctx);
}

char
rt_view_context_rotate_about_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_rotate_about_get(ctx);
    return rt_view_context_rotate_about_from_bsg(ctx);
}

int
rt_view_rotate_about_set_bsg(struct bsg_view *v, char rotate_about)
{
    if (!v)
	return 0;

    v->gv_rotate_about = rotate_about;
    return 1;
}

int
rt_view_context_rotate_about_set_bsg(void *ctx, char rotate_about)
{
    return rt_view_rotate_about_set_bsg((struct bsg_view *)ctx, rotate_about);
}

int
rt_view_context_rotate_about_set(void *ctx, char rotate_about)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_rotate_about_set(ctx, rotate_about);
    return rt_view_context_rotate_about_set_bsg(ctx, rotate_about);
}

char
rt_view_coord_from_bsg(const struct bsg_view *v)
{
    return v ? v->gv_coord : 'v';
}

char
rt_view_context_coord_from_bsg(const void *ctx)
{
    return rt_view_coord_from_bsg((const struct bsg_view *)ctx);
}

char
rt_view_context_coord_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_coord_get(ctx);
    return rt_view_context_coord_from_bsg(ctx);
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
rt_view_context_coord_set_bsg(void *ctx, char coord)
{
    return rt_view_coord_set_bsg((struct bsg_view *)ctx, coord);
}

int
rt_view_context_coord_set(void *ctx, char coord)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_coord_set(ctx, coord);
    return rt_view_context_coord_set_bsg(ctx, coord);
}

int
rt_view_snap_lines_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_lines(v);
}

int
rt_view_context_snap_lines_from_bsg(const void *ctx)
{
    return rt_view_snap_lines_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_snap_lines_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_lines_get(ctx);
    return rt_view_context_snap_lines_from_bsg(ctx);
}

int
rt_view_snap_lines_set_bsg(struct bsg_view *v, int enabled)
{
    if (!v)
	return 0;

    bsg_view_set_snap_lines(v, enabled);
    return 1;
}

int
rt_view_context_snap_lines_set_bsg(void *ctx, int enabled)
{
    return rt_view_snap_lines_set_bsg((struct bsg_view *)ctx, enabled);
}

int
rt_view_context_snap_lines_set(void *ctx, int enabled)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_lines_set(ctx, enabled);
    return rt_view_context_snap_lines_set_bsg(ctx, enabled);
}

static bsg_feature_ref
rt_view_feature_ref_to_bsg(rt_view_feature_ref ref)
{
    bsg_feature_ref bsg_ref = {ref.token, ref.revision};
    return bsg_ref;
}

static rt_view_feature_ref
rt_view_feature_ref_from_bsg(bsg_feature_ref ref)
{
    rt_view_feature_ref rt_ref = {ref.token, ref.revision};
    return rt_ref;
}

static enum bsg_feature_family
rt_view_feature_family_to_bsg(enum rt_view_feature_family family) {
    switch (family)
    {
	case RT_VIEW_FEATURE_TRANSIENT_PREVIEW:
		    return BSG_FEATURE_TRANSIENT_PREVIEW;
	case RT_VIEW_FEATURE_UNKNOWN:
	default:
	    return BSG_FEATURE_UNKNOWN;
    }
}

static bsg_edit_preview_op
rt_view_edit_preview_event_to_bsg(enum rt_view_edit_preview_event event)
{
    switch (event) {
	case RT_VIEW_EDIT_PREVIEW_BEGIN:
	    return BSG_EDIT_PREVIEW_BEGIN;
	case RT_VIEW_EDIT_PREVIEW_UPDATE:
	    return BSG_EDIT_PREVIEW_UPDATE;
	case RT_VIEW_EDIT_PREVIEW_COMMIT:
	    return BSG_EDIT_PREVIEW_COMMIT;
	case RT_VIEW_EDIT_PREVIEW_CANCEL:
	    return BSG_EDIT_PREVIEW_CANCEL;
	case RT_VIEW_EDIT_PREVIEW_REPLACE_SOURCE:
	    return BSG_EDIT_PREVIEW_REPLACE_SOURCE;
	case RT_VIEW_EDIT_PREVIEW_DISCARD:
	    return BSG_EDIT_PREVIEW_DISCARD;
	default:
	    return BSG_EDIT_PREVIEW_UPDATE;
    }
}

struct rt_view_edit_preview_callback_bridge {
    void *preview_ctx;
    struct rt_view_edit_preview_callbacks callbacks;
};

static void
rt_view_edit_preview_bridge_free(void *bridge_ctx)
{
    if (bridge_ctx)
	bu_free(bridge_ctx, "RT view edit preview callback bridge");
}

static uint64_t
rt_view_edit_preview_bridge_revision(void *bridge_ctx)
{
    struct rt_view_edit_preview_callback_bridge *bridge =
	(struct rt_view_edit_preview_callback_bridge *)bridge_ctx;
    if (!bridge || !bridge->callbacks.revision_cb)
	return 0;

    return bridge->callbacks.revision_cb(bridge->preview_ctx);
}

static int
rt_view_edit_preview_bridge_update(void *bridge_ctx, struct bsg_view *UNUSED(v))
{
    struct rt_view_edit_preview_callback_bridge *bridge =
	(struct rt_view_edit_preview_callback_bridge *)bridge_ctx;
    if (!bridge || !bridge->callbacks.update_cb)
	return 0;

    return bridge->callbacks.update_cb(bridge->preview_ctx);
}

static int
rt_view_edit_preview_bridge_pick(void *bridge_ctx, struct bsg_view *UNUSED(v),
				 int x, int y, void *pick_out)
{
    struct rt_view_edit_preview_callback_bridge *bridge =
	(struct rt_view_edit_preview_callback_bridge *)bridge_ctx;
    if (!bridge || !bridge->callbacks.pick_cb)
	return 0;

    return bridge->callbacks.pick_cb(bridge->preview_ctx, x, y, pick_out);
}

static int
rt_view_edit_preview_callbacks_to_bsg(struct bsg_edit_preview_ops *ops,
				      const struct rt_view_edit_preview_callbacks *callbacks,
				      void *preview_ctx)
{
    if (!ops || !callbacks)
	return 0;
    if (!callbacks->revision_cb && !callbacks->update_cb &&
	!callbacks->pick_cb)
	return 1;

    struct rt_view_edit_preview_callback_bridge *bridge =
	(struct rt_view_edit_preview_callback_bridge *)bu_calloc(1,
	    sizeof(struct rt_view_edit_preview_callback_bridge),
	    "RT view edit preview callback bridge");
    bridge->preview_ctx = preview_ctx;
    bridge->callbacks = *callbacks;

    ops->preview_ctx = bridge;
    ops->owns_preview_ctx = 1;
    ops->revision_cb = callbacks->revision_cb ?
		       rt_view_edit_preview_bridge_revision : NULL;
    ops->update_cb = callbacks->update_cb ?
		     rt_view_edit_preview_bridge_update : NULL;
    ops->pick_cb = callbacks->pick_cb ?
		   rt_view_edit_preview_bridge_pick : NULL;
    ops->free_cb = rt_view_edit_preview_bridge_free;
    return 1;
}

int
rt_view_feature_ref_is_null_bsg(rt_view_feature_ref ref)
{
    return bsg_feature_ref_is_null(rt_view_feature_ref_to_bsg(ref));
}

int
rt_view_edit_preview_publish_event_bsg(struct bsg_view *v,
				       rt_view_feature_ref feature,
				       enum rt_view_edit_preview_event event,
				       const char *source_path)
{
    struct bsg_interaction_record *record =
	bsg_interaction_edit_preview_record(v, rt_view_feature_ref_to_bsg(feature),
					    rt_view_edit_preview_event_to_bsg(event), source_path);
    if (!record)
	return 0;

    bsg_interaction_record_free(record);
    return 1;
}

int
rt_view_context_edit_preview_publish_event_bsg(
    void *ctx,
    rt_view_feature_ref feature,
    enum rt_view_edit_preview_event event,
    const char *source_path)
{
    return rt_view_edit_preview_publish_event_bsg((struct bsg_view *)ctx,
	    feature, event, source_path);
}

rt_view_feature_ref
rt_view_feature_overlay_ensure_bsg(struct bsg_view *v,
				   const char *name,
				   const void *owner,
				   void *preview_ctx,
				   const struct rt_view_edit_preview_callbacks *callbacks,
				   const char *source_path)
{
    if (!v || !name)
	return RT_VIEW_FEATURE_REF_NULL;

    bsg_feature_ref ref = bsg_feature_find(v, name);
    if (bsg_feature_ref_is_null(ref)) {
	ref = bsg_feature_create_overlay(v, name, 1 /* local */);
	if (!bsg_feature_ref_is_null(ref)) {
	    bsg_feature_overlay_register_owner(ref, owner,
					       BSG_OVERLAY_ROLE_MODEL,
					       BSG_OVERLAY_CLASS_EDIT_HANDLE,
					       BSG_OVERLAY_LC_PER_TOOL,
					       BSG_OVERLAY_ORDER_POST_TRANSPARENT,
					       NULL, 0);

	    if (callbacks) {
		struct bsg_edit_preview_ops ops = BSG_EDIT_PREVIEW_OPS_INIT;
		if (rt_view_edit_preview_callbacks_to_bsg(&ops, callbacks,
			preview_ctx)) {
		    int attached = bsg_feature_edit_preview_attach(ref,
				   preview_ctx, NULL, &ops);
		    if (!attached && ops.owns_preview_ctx && ops.preview_ctx)
			rt_view_edit_preview_bridge_free(ops.preview_ctx);
		    rt_view_edit_preview_publish_event_bsg(v,
							   rt_view_feature_ref_from_bsg(ref),
							   RT_VIEW_EDIT_PREVIEW_BEGIN, source_path);
		}
	    }
	}
    }

    return rt_view_feature_ref_from_bsg(ref);
}

rt_view_feature_ref
rt_view_context_feature_overlay_ensure_bsg(
    void *ctx,
    const char *name,
    const void *owner,
    void *preview_ctx,
    const struct rt_view_edit_preview_callbacks *callbacks,
    const char *source_path)
{
    return rt_view_feature_overlay_ensure_bsg((struct bsg_view *)ctx, name,
	    owner, preview_ctx, callbacks, source_path);
}

rt_view_feature_ref
rt_view_feature_label_ensure_bsg(struct bsg_view *v,
				 const char *name,
				 const void *owner)
{
    if (!v || !name)
	return RT_VIEW_FEATURE_REF_NULL;

    bsg_feature_ref ref = bsg_feature_find(v, name);
    if (bsg_feature_ref_is_null(ref))
	ref = bsg_feature_create_label(v, name, 1 /* local */);
    if (!bsg_feature_ref_is_null(ref)) {
	bsg_feature_overlay_register_owner(ref, owner,
					   BSG_OVERLAY_ROLE_MODEL,
					   BSG_OVERLAY_CLASS_EDIT_HANDLE,
					   BSG_OVERLAY_LC_PER_TOOL,
					   BSG_OVERLAY_ORDER_POST_TRANSPARENT,
					   NULL, 1);
    }

    return rt_view_feature_ref_from_bsg(ref);
}

rt_view_feature_ref
rt_view_context_feature_label_ensure_bsg(void *ctx,
	const char *name,
	const void *owner)
{
    return rt_view_feature_label_ensure_bsg((struct bsg_view *)ctx, name,
					    owner);
}

int
rt_view_feature_remove_bsg(struct bsg_view *v, const char *name)
{
    return (v && name) ? bsg_feature_remove(v, name) : 0;
}

int
rt_view_context_feature_remove_bsg(void *ctx, const char *name)
{
    return rt_view_feature_remove_bsg((struct bsg_view *)ctx, name);
}

int
rt_view_feature_summary_bsg(struct bsg_view *v,
			    const char *name,
			    struct rt_view_feature_summary_bsg *summary)
{
    if (!summary)
	return 0;

    memset(summary, 0, sizeof(*summary));
    if (!v || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(v, name);
    if (bsg_feature_ref_is_null(ref))
	return 1;

    struct bsg_feature_record record;
    memset(&record, 0, sizeof(record));
    if (!bsg_feature_record_get(ref, &record))
	return 0;

    summary->exists = 1;
    summary->is_overlay = (record.family == BSG_FEATURE_OVERLAY);
    summary->is_label = (record.family == BSG_FEATURE_LABEL);
    summary->is_transient_preview =
	(record.family == BSG_FEATURE_TRANSIENT_PREVIEW);
    summary->visible = record.visible;
    summary->color[0] = record.color[0];
    summary->color[1] = record.color[1];
    summary->color[2] = record.color[2];
    summary->child_count = record.child_count;
    summary->geometry_command_count = record.geometry_command_count;

    return 1;
}

int
rt_view_context_feature_summary_bsg(
    void *ctx,
    const char *name,
    struct rt_view_feature_summary_bsg *summary)
{
    return rt_view_feature_summary_bsg((struct bsg_view *)ctx, name, summary);
}

int
rt_view_context_feature_geometry_summary_bsg(
    void *ctx,
    const char *name,
    struct rt_view_feature_geometry_summary_bsg *summary)
{
    struct bsg_view *v = (struct bsg_view *)ctx;
    point_t *points = NULL;
    int *cmds = NULL;
    size_t point_count = 0;
    int ret = 0;

    if (!summary)
	return 0;

    memset(summary, 0, sizeof(*summary));
    if (!v || !name)
	return 0;

    bsg_feature_ref ref = bsg_feature_find(v, name);
    if (bsg_feature_ref_is_null(ref))
	return 1;

    struct bsg_feature_record record;
    memset(&record, 0, sizeof(record));
    if (!bsg_feature_record_get(ref, &record))
	return 0;

    if (bsg_feature_points_copy(ref, &points, &cmds, &point_count)) {
	summary->exists = 1;
	summary->point_count = point_count;
	summary->command_count = record.geometry_command_count;
	ret = 1;
    }

    if (points)
	bu_free(points, "RT view feature geometry summary points");
    if (cmds)
	bu_free(cmds, "RT view feature geometry summary commands");
    return ret;
}

int
rt_view_context_feature_geometry_summary(
    void *ctx,
    const char *name,
    struct rt_view_feature_geometry_summary *summary)
{
    struct rt_view_feature_geometry_summary_bsg bsg_summary =
	    RT_VIEW_FEATURE_GEOMETRY_SUMMARY_BSG_INIT;

    if (summary) {
	summary->exists = 0;
	summary->point_count = 0;
	summary->command_count = 0;
    }

    if (_rt_view_context_native_is(ctx))
	return 0;

    if (!summary ||
	!rt_view_context_feature_geometry_summary_bsg(ctx, name,
		&bsg_summary))
	return 0;

    summary->exists = bsg_summary.exists;
    summary->point_count = bsg_summary.point_count;
    summary->command_count = bsg_summary.command_count;
    return 1;
}

void
rt_view_feature_set_view_bsg(rt_view_feature_ref ref, struct bsg_view *v)
{
    if (!rt_view_feature_ref_is_null_bsg(ref))
	bsg_feature_set_view(rt_view_feature_ref_to_bsg(ref), v);
}

void
rt_view_feature_set_context_bsg(rt_view_feature_ref ref, void *ctx)
{
    rt_view_feature_set_view_bsg(ref, (struct bsg_view *)ctx);
}

void
rt_view_feature_set_visible_bsg(rt_view_feature_ref ref, int visible)
{
    if (!rt_view_feature_ref_is_null_bsg(ref))
	bsg_feature_set_visible(rt_view_feature_ref_to_bsg(ref), visible);
}

void
rt_view_feature_set_color_bsg(rt_view_feature_ref ref, int r, int g, int b)
{
    if (!rt_view_feature_ref_is_null_bsg(ref))
	bsg_feature_set_color(rt_view_feature_ref_to_bsg(ref), r, g, b);
}

int
rt_view_feature_touch_bsg(rt_view_feature_ref ref)
{
    return rt_view_feature_ref_is_null_bsg(ref) ? 0 :
	   bsg_feature_edit_preview_touch(rt_view_feature_ref_to_bsg(ref));
}

int
rt_view_feature_labels_replace_bsg(rt_view_feature_ref ref,
				   const struct rt_view_feature_label *labels,
				   size_t label_count)
{
    if (rt_view_feature_ref_is_null_bsg(ref))
	return 0;
    if (label_count && !labels)
	return 0;

    struct bsg_feature_label_data *bsg_labels = NULL;
    if (label_count > 0) {
	bsg_labels = (struct bsg_feature_label_data *)bu_calloc(label_count,
		     sizeof(struct bsg_feature_label_data),
		     "RT view feature labels");
	for (size_t i = 0; i < label_count; i++) {
	    struct bsg_feature_label_data init = BSG_FEATURE_LABEL_DATA_INIT;
	    bsg_labels[i] = init;
	    bsg_labels[i].text = labels[i].text;
	    VMOVE(bsg_labels[i].point, labels[i].point);
	    bsg_labels[i].color_valid = labels[i].color_valid;
	    bsg_labels[i].color[0] = labels[i].color[0];
	    bsg_labels[i].color[1] = labels[i].color[1];
	    bsg_labels[i].color[2] = labels[i].color[2];
	    bsg_labels[i].anchor = BSG_ANCHOR_AUTO;
	    bsg_labels[i].font_size = labels[i].font_size;
	}
    }

    int ret = bsg_feature_labels_replace(rt_view_feature_ref_to_bsg(ref),
					 bsg_labels, label_count);
    if (bsg_labels)
	bu_free(bsg_labels, "RT view feature labels");
    return ret;
}

int
rt_view_feature_points_replace_bsg(rt_view_feature_ref ref,
				   enum rt_view_feature_family family,
				   const point_t *points,
				   const int *cmds,
				   size_t point_count)
{
    if (rt_view_feature_ref_is_null_bsg(ref))
	return 0;

    return bsg_feature_points_replace(rt_view_feature_ref_to_bsg(ref),
				      rt_view_feature_family_to_bsg(family),
				      point_count ? points : NULL,
				      point_count ? cmds : NULL,
				      point_count);
}

int
rt_view_feature_clear_geometry_bsg(rt_view_feature_ref ref)
{
    if (rt_view_feature_ref_is_null_bsg(ref))
	return 0;

    return bsg_feature_points_replace(rt_view_feature_ref_to_bsg(ref),
				      BSG_FEATURE_UNKNOWN, NULL, NULL, 0);
}

int
rt_view_feature_ref_is_null(rt_view_feature_ref ref)
{
    return ref.token ? 0 : 1;
}

int
rt_view_context_edit_preview_publish_event(
    void *ctx,
    rt_view_feature_ref feature,
    enum rt_view_edit_preview_event event,
    const char *source_path)
{
    struct rt_view_context_feature_adapter adapter;
    if (rt_view_context_feature_adapter_get(ctx, &adapter) &&
	adapter.edit_preview_publish_event &&
	(!adapter.owns_ref ||
	 adapter.owns_ref(feature, adapter.data))) {
	return adapter.edit_preview_publish_event(ctx, feature, event,
		source_path, adapter.data);
    }

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_edit_preview_publish_event_bsg(ctx, feature, event,
	    source_path);
}

rt_view_feature_ref
rt_view_context_feature_overlay_ensure(
    void *ctx,
    const char *name,
    const void *owner,
    void *preview_ctx,
    const struct rt_view_edit_preview_callbacks *callbacks,
    const char *source_path)
{
    struct rt_view_context_feature_adapter adapter;
    if (rt_view_context_feature_adapter_get(ctx, &adapter) &&
	adapter.overlay_ensure)
	return adapter.overlay_ensure(ctx, name, owner, preview_ctx,
				      callbacks, source_path, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_FEATURE_REF_NULL;

    return rt_view_context_feature_overlay_ensure_bsg(ctx, name, owner,
	    preview_ctx, callbacks, source_path);
}

rt_view_feature_ref
rt_view_context_feature_label_ensure(void *ctx,
				     const char *name,
				     const void *owner)
{
    struct rt_view_context_feature_adapter adapter;
    if (rt_view_context_feature_adapter_get(ctx, &adapter) &&
	adapter.label_ensure)
	return adapter.label_ensure(ctx, name, owner, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_FEATURE_REF_NULL;

    return rt_view_context_feature_label_ensure_bsg(ctx, name, owner);
}

int
rt_view_context_feature_remove(void *ctx, const char *name)
{
    struct rt_view_context_feature_adapter adapter;
    if (rt_view_context_feature_adapter_get(ctx, &adapter) &&
	adapter.remove)
	return adapter.remove(ctx, name, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_feature_remove_bsg(ctx, name);
}

void
rt_view_feature_set_context(rt_view_feature_ref ref, void *ctx)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter)) {
	if (adapter.set_context)
	    (void)adapter.set_context(ref, ctx, adapter.data);
	return;
    }

    if (_rt_view_context_native_is(ctx))
	return;

    rt_view_feature_set_context_bsg(ref, ctx);
}

void
rt_view_feature_set_visible(rt_view_feature_ref ref, int visible)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter)) {
	if (adapter.set_visible)
	    (void)adapter.set_visible(ref, visible, adapter.data);
	return;
    }

    rt_view_feature_set_visible_bsg(ref, visible);
}

void
rt_view_feature_set_color(rt_view_feature_ref ref, int r, int g, int b)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter)) {
	if (adapter.set_color)
	    (void)adapter.set_color(ref, r, g, b, adapter.data);
	return;
    }

    rt_view_feature_set_color_bsg(ref, r, g, b);
}

int
rt_view_feature_touch(rt_view_feature_ref ref)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter))
	return adapter.touch ? adapter.touch(ref, adapter.data) : 0;

    return rt_view_feature_touch_bsg(ref);
}

int
rt_view_feature_labels_replace(rt_view_feature_ref ref,
			       const struct rt_view_feature_label *labels,
			       size_t label_count)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter))
	return adapter.labels_replace ?
	       adapter.labels_replace(ref, labels, label_count,
				      adapter.data) : 0;

    return rt_view_feature_labels_replace_bsg(ref, labels, label_count);
}

int
rt_view_feature_points_replace(rt_view_feature_ref ref,
			       enum rt_view_feature_family family,
			       const point_t *points,
			       const int *cmds,
			       size_t point_count)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter))
	return adapter.points_replace ?
	       adapter.points_replace(ref, family, points, cmds, point_count,
				      adapter.data) : 0;

    return rt_view_feature_points_replace_bsg(ref, family, points, cmds,
	    point_count);
}

int
rt_view_feature_clear_geometry(rt_view_feature_ref ref)
{
    struct rt_view_context_feature_adapter adapter;
    if (_rt_view_feature_adapter_for_ref(ref, &adapter))
	return adapter.clear_geometry ?
	       adapter.clear_geometry(ref, adapter.data) : 0;

    return rt_view_feature_clear_geometry_bsg(ref);
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
rt_view_context_polygon_create_bsg(void *ctx, int type, point_t *fp)
{
    return rt_view_polygon_create_bsg((struct bsg_view *)ctx, type, fp);
}

rt_view_polygon_ref
rt_view_polygon_select_bsg(struct bsg_view *v, point_t *cp)
{
    if (!v || !cp)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_view_select_polygon_ref(v, cp));
}

rt_view_polygon_ref
rt_view_context_polygon_select_bsg(void *ctx, point_t *cp)
{
    return rt_view_polygon_select_bsg((struct bsg_view *)ctx, cp);
}

rt_view_polygon_ref
rt_view_polygon_find_bsg(struct bsg_view *v, const char *name)
{
    if (!v || !name)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_view_polygon_find_ref(v, name));
}

rt_view_polygon_ref
rt_view_context_polygon_find_bsg(void *ctx, const char *name)
{
    return rt_view_polygon_find_bsg((struct bsg_view *)ctx, name);
}

rt_view_polygon_ref
rt_view_polygon_dup_bsg(struct bsg_view *v, const char *name, const char *new_name)
{
    if (!v || !name || !new_name)
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_ref_from_bsg(bsg_view_polygon_dup_ref(v, name, new_name));
}

rt_view_polygon_ref
rt_view_context_polygon_dup_bsg(void *ctx, const char *name, const char *new_name)
{
    return rt_view_polygon_dup_bsg((struct bsg_view *)ctx, name, new_name);
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

void
rt_view_context_polygon_visit_records_bsg(
    void *ctx,
    rt_view_polygon_record_callback_bsg_t callback,
    void *data)
{
    rt_view_polygon_visit_records_bsg((struct bsg_view *)ctx, callback, data);
}

size_t
rt_view_polygon_snap_count_bsg(struct bsg_view *v, rt_view_polygon_ref exclude)
{
    if (!v)
	return 0;

    return bsg_view_polygon_snap_count(v, rt_view_polygon_ref_to_bsg(exclude));
}

size_t
rt_view_context_polygon_snap_count_bsg(void *ctx, rt_view_polygon_ref exclude)
{
    return rt_view_polygon_snap_count_bsg((struct bsg_view *)ctx, exclude);
}

int
rt_view_polygon_clear_point_selection_bsg(struct bsg_view *v)
{
    if (!v)
	return 0;

    return bsg_view_polygon_clear_point_selection(v);
}

int
rt_view_context_polygon_clear_point_selection_bsg(void *ctx)
{
    return rt_view_polygon_clear_point_selection_bsg((struct bsg_view *)ctx);
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
rt_view_polygon_update_context_bsg(rt_view_polygon_ref ref, void *ctx, int utype)
{
    return rt_view_polygon_update_bsg(ref, (struct bsg_view *)ctx, utype);
}

int
rt_view_polygon_update_screen_pt_bsg(rt_view_polygon_ref ref, struct bsg_view *v,
				     int x, int y, int utype)
{
    return bsg_polygon_update_screen_pt(rt_view_polygon_ref_to_bsg(ref), v, x, y, utype);
}

int
rt_view_polygon_update_screen_pt_context_bsg(rt_view_polygon_ref ref, void *ctx,
	int x, int y, int utype)
{
    return rt_view_polygon_update_screen_pt_bsg(ref, (struct bsg_view *)ctx, x,
	    y, utype);
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
rt_view_polygon_set_context_bsg(rt_view_polygon_ref ref, void *ctx)
{
    return rt_view_polygon_set_view_bsg(ref, (struct bsg_view *)ctx);
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

    return db_sketch_to_view_polygon_ref(name, dbip, dp, v);
}

rt_view_polygon_ref
rt_view_polygon_import_sketch_context_bsg(const char *name,
	struct db_i *dbip,
	struct directory *dp,
	void *ctx)
{
    return rt_view_polygon_import_sketch_bsg(name, dbip, dp,
	    (struct bsg_view *)ctx);
}

struct directory *
rt_view_polygon_export_sketch_bsg(struct db_i *dbip, const char *name,
				  rt_view_polygon_ref ref)
{
    if (!dbip || !name)
	return NULL;

    return db_view_polygon_ref_to_sketch(dbip, name, ref);
}

int
rt_view_polygon_snap_exclude_set_bsg(struct bsg_view *v, rt_view_polygon_ref ref)
{
    rt_view_feature_ref feature_ref = {ref.token, ref.revision};
    return rt_view_snap_exclude_feature_set_bsg(v, feature_ref);
}

int
rt_view_context_polygon_snap_exclude_set_bsg(void *ctx, rt_view_polygon_ref ref)
{
    return rt_view_polygon_snap_exclude_set_bsg((struct bsg_view *)ctx, ref);
}

int
rt_view_polygon_ref_is_null(rt_view_polygon_ref ref)
{
    return ref.token ? 0 : 1;
}

int
rt_view_polygon_record_get(rt_view_polygon_ref ref,
			   struct rt_view_polygon_record *record)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.record_get ?
	       adapter.record_get(ref, record, adapter.data) : 0;

    return rt_view_polygon_record_get_bsg(ref, record);
}

rt_view_polygon_ref
rt_view_context_polygon_create(void *ctx, int type, point_t *fp)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.create) {
	rt_view_polygon_ref ref = adapter.create(ctx, type, fp,
				  adapter.data);
	if (!rt_view_polygon_ref_is_null(ref))
	    return ref;
    }

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_context_polygon_create_bsg(ctx, type, fp);
}

rt_view_polygon_ref
rt_view_context_polygon_select(void *ctx, point_t *cp)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.select) {
	rt_view_polygon_ref ref = adapter.select(ctx, cp, adapter.data);
	if (!rt_view_polygon_ref_is_null(ref))
	    return ref;
    }

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_context_polygon_select_bsg(ctx, cp);
}

rt_view_polygon_ref
rt_view_context_polygon_find(void *ctx, const char *name)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.find) {
	rt_view_polygon_ref ref = adapter.find(ctx, name, adapter.data);
	if (!rt_view_polygon_ref_is_null(ref))
	    return ref;
    }

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_context_polygon_find_bsg(ctx, name);
}

rt_view_polygon_ref
rt_view_context_polygon_dup(void *ctx,
			    const char *name,
			    const char *new_name)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.dup) {
	rt_view_polygon_ref ref = adapter.dup(ctx, name, new_name,
					      adapter.data);
	if (!rt_view_polygon_ref_is_null(ref))
	    return ref;
    }

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_context_polygon_dup_bsg(ctx, name, new_name);
}

void
rt_view_context_polygon_visit_records(
    void *ctx,
    rt_view_polygon_record_callback_t callback,
    void *data)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.visit_records) {
	adapter.visit_records(ctx, callback, data, adapter.data);
	return;
    }

    if (_rt_view_context_native_is(ctx))
	return;

    rt_view_context_polygon_visit_records_bsg(ctx, callback, data);
}

size_t
rt_view_context_polygon_snap_count(void *ctx, rt_view_polygon_ref exclude)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.snap_count)
	return adapter.snap_count(ctx, exclude, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_polygon_snap_count_bsg(ctx, exclude);
}

int
rt_view_context_polygon_clear_point_selection(void *ctx)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.clear_point_selection)
	return adapter.clear_point_selection(ctx, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_polygon_clear_point_selection_bsg(ctx);
}

int
rt_view_polygon_update_context(rt_view_polygon_ref ref, void *ctx, int utype)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.update ?
	       adapter.update(ref, ctx, utype, adapter.data) : 0;

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_polygon_update_context_bsg(ref, ctx, utype);
}

int
rt_view_polygon_update_screen_pt_context(rt_view_polygon_ref ref,
	void *ctx,
	int x,
	int y,
	int utype)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.update_screen_pt ?
	       adapter.update_screen_pt(ref, ctx, x, y, utype, adapter.data) :
	       0;

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_polygon_update_screen_pt_context_bsg(ref, ctx, x, y, utype);
}

int
rt_view_polygon_move(rt_view_polygon_ref ref,
		     point_t *current_point,
		     point_t *previous_point)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.move ?
	       adapter.move(ref, current_point, previous_point,
			    adapter.data) : 0;

    return rt_view_polygon_move_bsg(ref, current_point, previous_point);
}

int
rt_view_polygon_set_name(rt_view_polygon_ref ref, const char *name)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.set_name ?
	       adapter.set_name(ref, name, adapter.data) : 0;

    return rt_view_polygon_set_name_bsg(ref, name);
}

int
rt_view_polygon_set_context(rt_view_polygon_ref ref, void *ctx)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.set_context ?
	       adapter.set_context(ref, ctx, adapter.data) : 0;

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_polygon_set_context_bsg(ref, ctx);
}

int
rt_view_polygon_set_visual(rt_view_polygon_ref ref,
			   const struct bu_color *edge_color,
			   const struct bu_color *fill_color,
			   fastf_t fill_slope_x,
			   fastf_t fill_slope_y,
			   fastf_t fill_density,
			   fastf_t vZ,
			   int fill_flag)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.set_visual ?
	       adapter.set_visual(ref, edge_color, fill_color, fill_slope_x,
				  fill_slope_y, fill_density, vZ, fill_flag,
				  adapter.data) : 0;

    return rt_view_polygon_set_visual_bsg(ref, edge_color, fill_color,
					  fill_slope_x, fill_slope_y, fill_density, vZ, fill_flag);
}

int
rt_view_polygon_set_open(rt_view_polygon_ref ref, int open)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.set_open ?
	       adapter.set_open(ref, open, adapter.data) : 0;

    return rt_view_polygon_set_open_bsg(ref, open);
}

int
rt_view_polygon_close(rt_view_polygon_ref ref)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.close ? adapter.close(ref, adapter.data) : 0;

    return rt_view_polygon_close_bsg(ref);
}

int
rt_view_polygon_clear_selected_point(rt_view_polygon_ref ref)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.clear_selected_point ?
	       adapter.clear_selected_point(ref, adapter.data) : 0;

    return rt_view_polygon_clear_selected_point_bsg(ref);
}

int
rt_view_polygon_remove(rt_view_polygon_ref ref)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.remove ? adapter.remove(ref, adapter.data) : 0;

    return rt_view_polygon_remove_bsg(ref);
}

void *
rt_view_polygon_user_data(rt_view_polygon_ref ref)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.user_data ?
	       adapter.user_data(ref, adapter.data) : NULL;

    return rt_view_polygon_user_data_bsg(ref);
}

int
rt_view_polygon_user_data_set(rt_view_polygon_ref ref, void *user_data)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.user_data_set ?
	       adapter.user_data_set(ref, user_data, adapter.data) : 0;

    return rt_view_polygon_user_data_set_bsg(ref, user_data);
}

int
rt_view_polygon_csg(rt_view_polygon_ref target,
		    rt_view_polygon_ref stencil,
		    bg_clip_t op)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(target, &adapter))
	return adapter.csg ?
	       adapter.csg(target, stencil, op, adapter.data) : 0;

    return rt_view_polygon_csg_bsg(target, stencil, op);
}

rt_view_polygon_ref
rt_view_polygon_import_sketch_context(const char *name,
				      struct db_i *dbip,
				      struct directory *dp,
				      void *ctx)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.import_sketch_context) {
	rt_view_polygon_ref ref = adapter.import_sketch_context(name, dbip,
				  dp, ctx, adapter.data);
	if (!rt_view_polygon_ref_is_null(ref))
	    return ref;
    }

    if (_rt_view_context_native_is(ctx))
	return RT_VIEW_POLYGON_REF_NULL;

    return rt_view_polygon_import_sketch_context_bsg(name, dbip, dp, ctx);
}

struct directory *
rt_view_polygon_export_sketch(struct db_i *dbip,
			      const char *name,
			      rt_view_polygon_ref ref)
{
    struct rt_view_context_polygon_adapter adapter;
    if (_rt_view_polygon_adapter_for_ref(ref, &adapter))
	return adapter.export_sketch ?
	       adapter.export_sketch(dbip, name, ref, adapter.data) : NULL;

    return rt_view_polygon_export_sketch_bsg(dbip, name, ref);
}

int
rt_view_context_polygon_snap_exclude_set(void *ctx, rt_view_polygon_ref ref)
{
    struct rt_view_context_polygon_adapter adapter;
    if (rt_view_context_polygon_adapter_get(ctx, &adapter) &&
	adapter.snap_exclude_set)
	return adapter.snap_exclude_set(ctx, ref, adapter.data);

    if (_rt_view_context_native_is(ctx))
	return 0;

    return rt_view_context_polygon_snap_exclude_set_bsg(ctx, ref);
}

int
rt_view_snap_source_flags_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_source_flags(v);
}

int
rt_view_context_snap_source_flags_from_bsg(const void *ctx)
{
    return rt_view_snap_source_flags_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_snap_source_flags_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_source_flags_get(ctx);
    return rt_view_context_snap_source_flags_from_bsg(ctx);
}

int
rt_view_snap_source_flags_set_bsg(struct bsg_view *v, int flags)
{
    if (!v)
	return 0;

    bsg_view_set_snap_source_flags(v, flags);
    return 1;
}

int
rt_view_context_snap_source_flags_set_bsg(void *ctx, int flags)
{
    return rt_view_snap_source_flags_set_bsg((struct bsg_view *)ctx, flags);
}

int
rt_view_context_snap_source_flags_set(void *ctx, int flags)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_source_flags_set(ctx, flags);
    return rt_view_context_snap_source_flags_set_bsg(ctx, flags);
}

unsigned long long
rt_view_snap_kind_mask_from_bsg(const struct bsg_view *v)
{
    return (unsigned long long)bsg_view_snap_kind_mask(v);
}

unsigned long long
rt_view_context_snap_kind_mask_from_bsg(const void *ctx)
{
    return rt_view_snap_kind_mask_from_bsg((const struct bsg_view *)ctx);
}

unsigned long long
rt_view_context_snap_kind_mask_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_kind_mask_get(ctx);
    return rt_view_context_snap_kind_mask_from_bsg(ctx);
}

int
rt_view_snap_exclude_feature_set_bsg(struct bsg_view *v,
				     rt_view_feature_ref ref)
{
    if (!v)
	return 0;

    bsg_view_snap_exclude_feature_set(v, rt_view_feature_ref_to_bsg(ref));
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

int
rt_view_context_snap_exclude_feature_clear_bsg(void *ctx)
{
    return rt_view_snap_exclude_feature_clear_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_snap_exclude_feature_clear(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_exclude_feature_clear(ctx);
    return rt_view_context_snap_exclude_feature_clear_bsg(ctx);
}

unsigned long long
rt_view_prepare_tcl_snap_bsg(struct bsg_view *v)
{
    return (unsigned long long)bsg_view_prepare_tcl_snap(v);
}

unsigned long long
rt_view_context_prepare_tcl_snap_bsg(void *ctx)
{
    return rt_view_prepare_tcl_snap_bsg((struct bsg_view *)ctx);
}

unsigned long long
rt_view_context_prepare_tcl_snap(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_prepare_tcl_snap(ctx);
    return rt_view_context_prepare_tcl_snap_bsg(ctx);
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
rt_view_context_center_linesnap_bsg(void *ctx)
{
    return rt_view_center_linesnap_bsg((struct bsg_view *)ctx);
}

int
rt_view_context_center_linesnap(void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return 0;
    return rt_view_context_center_linesnap_bsg(ctx);
}

int
rt_view_zclip_from_bsg(const struct bsg_view *v)
{
    return bsg_view_zclip(v);
}

int
rt_view_context_zclip_from_bsg(const void *ctx)
{
    return rt_view_zclip_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_zclip_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_zclip_get(ctx);
    return rt_view_context_zclip_from_bsg(ctx);
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
rt_view_context_zclip_set_bsg(void *ctx, int zclip)
{
    return rt_view_zclip_set_bsg((struct bsg_view *)ctx, zclip);
}

int
rt_view_context_zclip_set(void *ctx, int zclip)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_zclip_set(ctx, zclip);
    return rt_view_context_zclip_set_bsg(ctx, zclip);
}

int
rt_view_framebuffer_mode_from_bsg(const struct bsg_view *v)
{
    return bsg_view_framebuffer_mode(v);
}

int
rt_view_context_framebuffer_mode_from_bsg(const void *ctx)
{
    return rt_view_framebuffer_mode_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_framebuffer_mode_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_framebuffer_mode_get(ctx);
    return rt_view_context_framebuffer_mode_from_bsg(ctx);
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
rt_view_context_framebuffer_mode_set_bsg(void *ctx, int mode)
{
    return rt_view_framebuffer_mode_set_bsg((struct bsg_view *)ctx, mode);
}

int
rt_view_context_framebuffer_mode_set(void *ctx, int mode)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_framebuffer_mode_set(ctx, mode);
    return rt_view_context_framebuffer_mode_set_bsg(ctx, mode);
}

int
rt_view_cleared_from_bsg(const struct bsg_view *v)
{
    return bsg_view_cleared(v);
}

int
rt_view_context_cleared_from_bsg(const void *ctx)
{
    return rt_view_cleared_from_bsg((const struct bsg_view *)ctx);
}

int
rt_view_context_cleared_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_cleared_get(ctx);
    return rt_view_context_cleared_from_bsg(ctx);
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
rt_view_context_cleared_set_bsg(void *ctx, int cleared)
{
    return rt_view_cleared_set_bsg((struct bsg_view *)ctx, cleared);
}

int
rt_view_context_cleared_set(void *ctx, int cleared)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_cleared_set(ctx, cleared);
    return rt_view_context_cleared_set_bsg(ctx, cleared);
}

int
rt_view_settings_shared_bsg(const struct bsg_view *a, const struct bsg_view *b)
{
    return bsg_view_settings_shared(a, b);
}

int
rt_view_context_settings_shared_bsg(const void *a, const void *b)
{
    return rt_view_settings_shared_bsg((const struct bsg_view *)a,
				       (const struct bsg_view *)b);
}

int
rt_view_context_settings_shared(const void *a, const void *b)
{
    if (_rt_view_context_native_is(a) || _rt_view_context_native_is(b))
	return _rt_view_context_native_settings_shared(a, b);
    return rt_view_context_settings_shared_bsg(a, b);
}

double
rt_view_snap_tolerance_factor_from_bsg(const struct bsg_view *v)
{
    return bsg_view_snap_tolerance_factor(v);
}

double
rt_view_context_snap_tolerance_factor_from_bsg(const void *ctx)
{
    return rt_view_snap_tolerance_factor_from_bsg(
	       (const struct bsg_view *)ctx);
}

double
rt_view_context_snap_tolerance_factor_get(const void *ctx)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_tolerance_factor_get(ctx);
    return rt_view_context_snap_tolerance_factor_from_bsg(ctx);
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
rt_view_context_snap_tolerance_factor_set_bsg(void *ctx, double factor)
{
    return rt_view_snap_tolerance_factor_set_bsg((struct bsg_view *)ctx,
	    factor);
}

int
rt_view_context_snap_tolerance_factor_set(void *ctx, double factor)
{
    if (_rt_view_context_native_is(ctx))
	return _rt_view_context_native_snap_tolerance_factor_set(ctx,
		factor);
    return rt_view_context_snap_tolerance_factor_set_bsg(ctx, factor);
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
