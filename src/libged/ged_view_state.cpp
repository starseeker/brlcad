/*                   G E D _ V I E W _ S T A T E . C P P
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
/** @file libged/ged_view_state.cpp
 *
 * GED-owned opaque view-context state and façade calls.
 *
 * The lower-level context implementation is still exposed through the RT view
 * API while libbv ownership is being finalized, but GED treats this as its
 * active view-state store rather than a fallback drawing path.
 */

#include "common.h"

#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "dm/obol.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "rt/view.h"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

struct ged_view_state_storage {
    void *active_view_ctx;
    void *view_set_ctx;
    struct bu_ptbl free_views;
};

static struct ged_impl *
ged_view_impl(struct ged *gedp)
{
    return gedp ? gedp->i : NULL;
}

static const struct ged_impl *
ged_view_impl_const(const struct ged *gedp)
{
    return gedp ? gedp->i : NULL;
}

static struct ged_view_state_storage *
ged_view_state(struct ged *gedp)
{
    struct ged_impl *impl = ged_view_impl(gedp);
    return impl ? (struct ged_view_state_storage *)impl->ged_view_state_ctx : NULL;
}

static const struct ged_view_state_storage *
ged_view_state_const(const struct ged *gedp)
{
    const struct ged_impl *impl = ged_view_impl_const(gedp);
    return impl ? (const struct ged_view_state_storage *)impl->ged_view_state_ctx : NULL;
}

static struct ged_view_state_storage *
ged_view_state_create(struct ged *gedp)
{
    struct ged_impl *impl = ged_view_impl(gedp);
    if (!impl)
	return NULL;
    if (impl->ged_view_state_ctx)
	return (struct ged_view_state_storage *)impl->ged_view_state_ctx;

    struct ged_view_state_storage *state =
	(struct ged_view_state_storage *)bu_calloc(1, sizeof(*state),
		"GED view state");
    state->view_set_ctx = rt_view_set_context_create();
    BU_PTBL_INIT(&state->free_views);
    impl->ged_view_state_ctx = (void *)state;
    return state;
}

extern "C" void
ged_view_state_init(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state_create(gedp);
    if (!state || !state->view_set_ctx)
	return;

    void *default_view = rt_view_context_create_with_set(state->view_set_ctx);
    state->active_view_ctx = default_view;
    rt_view_context_name_set(default_view, "default");
    rt_view_set_context_add(state->view_set_ctx, default_view);
    ged_view_context_owned_add(gedp, default_view);
}

extern "C" void
ged_view_state_free(struct ged *gedp)
{
    struct ged_impl *impl = ged_view_impl(gedp);
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!impl || !state)
	return;

    impl->ged_view_state_ctx = NULL;
    state->active_view_ctx = NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&state->free_views); i++) {
	void *view_ctx = (void *)BU_PTBL_GET(&state->free_views, i);
	rt_view_context_free(view_ctx);
    }
    bu_ptbl_free(&state->free_views);
    rt_view_set_context_destroy(state->view_set_ctx);
    bu_free(state, "GED view state");
}

extern "C" GED_EXPORT void *
ged_view_active_ctx(const struct ged *gedp)
{
    const struct ged_view_state_storage *state = ged_view_state_const(gedp);
    return state ? state->active_view_ctx : NULL;
}

extern "C" GED_EXPORT void
ged_view_active_ctx_set(struct ged *gedp, void *view_ctx)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (state)
	state->active_view_ctx = view_ctx;
}

extern "C" GED_EXPORT int
ged_view_context_owned_add(struct ged *gedp, void *view_ctx)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !view_ctx)
	return 0;

    bu_ptbl_ins(&state->free_views, (long *)view_ctx);
    if (!rt_view_context_user_data_get(view_ctx))
	rt_view_context_user_data_set(view_ctx, gedp);
    (void)ged_draw_view_context_obol_scene_adapter_attach(gedp, view_ctx);
    (void)ged_draw_view_context_obol_feature_adapter_attach(gedp, view_ctx);
    (void)ged_draw_view_context_obol_polygon_adapter_attach(gedp, view_ctx);
    (void)ged_draw_view_context_obol_selection_adapter_attach(gedp, view_ctx);
    return 1;
}

extern "C" GED_EXPORT void *
ged_view_set_ctx(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    return state ? state->view_set_ctx : NULL;
}

extern "C" GED_EXPORT struct bu_ptbl *
ged_view_set_views_ctx(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    return state ? rt_view_set_context_views(state->view_set_ctx) : NULL;
}

extern "C" GED_EXPORT void *
ged_view_find_ctx(struct ged *gedp, const char *name)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !name)
	return NULL;

    return rt_view_set_context_find_view(state->view_set_ctx, name);
}

extern "C" GED_EXPORT void
ged_draw_view_context_info_get(struct rt_view_info *view_info,
			       const void *view_ctx)
{
    rt_view_context_info_get(view_info, view_ctx);
}

extern "C" GED_EXPORT int
ged_draw_view_context_lod_policy_get(ged_draw_view_lod_policy *policy,
				     const void *view_ctx)
{
    return rt_view_context_lod_policy_get(policy, view_ctx);
}

extern "C" GED_EXPORT int
ged_draw_view_context_lod_policy_apply(void *view_ctx,
				       const ged_draw_view_lod_policy *policy)
{
    return rt_view_context_lod_policy_apply(view_ctx, policy);
}

extern "C" GED_EXPORT int
ged_draw_view_context_lod_policy_apply_bot_threshold(
	void *view_ctx,
	const ged_draw_view_lod_policy *policy,
	size_t bot_threshold)
{
    if (!policy)
	return 0;

    ged_draw_view_lod_policy override_policy = *policy;
    override_policy.bot_threshold = bot_threshold;
    return ged_draw_view_context_lod_policy_apply(view_ctx, &override_policy);
}

static int
ged_draw_view_snap_kind_mask(enum ged_draw_view_snap_kind kind)
{
    switch (kind) {
	case GED_DRAW_VIEW_SNAP_GRID:
	    return RT_VIEW_SNAP_KIND_GRID;
	case GED_DRAW_VIEW_SNAP_ENDPOINT:
	    return RT_VIEW_SNAP_KIND_ENDPOINT;
	default:
	    return 0;
    }
}

extern "C" GED_EXPORT int
ged_draw_view_context_snap_first_candidate(void *view_ctx,
					   const point_t sample,
					   enum ged_draw_view_snap_kind kind,
					   point_t candidate)
{
    if (!candidate)
	return 0;

    VSET(candidate, 0.0, 0.0, 0.0);

    int snap_kind = ged_draw_view_snap_kind_mask(kind);
    if (!view_ctx || !sample || !snap_kind)
	return 0;

    return rt_view_context_snap_first_candidate(view_ctx, sample, snap_kind,
	    candidate);
}

extern "C" GED_EXPORT void
ged_draw_view_context_lod_bounds_callback_set(void *view_ctx)
{
    rt_view_context_lod_bounds_callback_set(view_ctx);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_scale_get(const void *view_ctx)
{
    return rt_view_context_scale_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_scale_set(void *view_ctx, fastf_t scale)
{
    return rt_view_context_scale_set(view_ctx, scale);
}

extern "C" GED_EXPORT int
ged_view_context_scale_state_set(void *view_ctx,
	fastf_t scale,
	fastf_t initial_scale,
	fastf_t absolute_scale,
	fastf_t size,
	fastf_t inverse_size)
{
    return rt_view_context_scale_state_set(view_ctx, scale,
	    initial_scale, absolute_scale, size, inverse_size);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_initial_scale_get(const void *view_ctx)
{
    return rt_view_context_initial_scale_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_initial_scale_set(void *view_ctx, fastf_t scale)
{
    return rt_view_context_initial_scale_set(view_ctx, scale);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_absolute_scale_get(const void *view_ctx)
{
    return rt_view_context_absolute_scale_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_absolute_scale_set(void *view_ctx, fastf_t scale)
{
    return rt_view_context_absolute_scale_set(view_ctx, scale);
}

extern "C" GED_EXPORT int
ged_view_context_update(void *view_ctx)
{
    return rt_view_context_update(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_is_independent(const void *view_ctx)
{
    return rt_view_context_is_independent(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_independent_scope_is_null(void *view_ctx, int create)
{
    return rt_view_context_independent_scope_is_null(view_ctx, create);
}

extern "C" GED_EXPORT void
ged_view_context_independent_scope_destroy(void *view_ctx)
{
    rt_view_context_independent_scope_destroy(view_ctx);
}

static int
ged_view_clear_flags_to_rt(int flags)
{
    int rt_flags = 0;

    if (flags & GED_VIEW_CLEAR_DB)
	rt_flags |= RT_VIEW_CLEAR_DB;
    if (flags & GED_VIEW_CLEAR_VIEW)
	rt_flags |= RT_VIEW_CLEAR_VIEW;
    if (flags & GED_VIEW_CLEAR_LOCAL)
	rt_flags |= RT_VIEW_CLEAR_LOCAL;

    return rt_flags;
}

extern "C" GED_EXPORT size_t
ged_view_context_clear(void *view_ctx, int flags)
{
    size_t cleared = ged_draw_obol_view_context_clear(view_ctx, flags);
    cleared += rt_view_context_clear(view_ctx, ged_view_clear_flags_to_rt(flags));
    return cleared;
}

extern "C" GED_EXPORT int
ged_view_context_cleared_get(const void *view_ctx)
{
    return rt_view_context_cleared_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_cleared_set(void *view_ctx, int cleared)
{
    return rt_view_context_cleared_set(view_ctx, cleared);
}

extern "C" GED_EXPORT unsigned long long
ged_view_context_hash(const void *view_ctx)
{
    return rt_view_context_hash(view_ctx);
}

extern "C" GED_EXPORT void *
ged_view_context_user_data_get(const void *view_ctx)
{
    return rt_view_context_user_data_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_user_data_set(void *view_ctx, void *user_data)
{
    return rt_view_context_user_data_set(view_ctx, user_data);
}

extern "C" GED_EXPORT int
ged_view_context_tclcad_data_set(void *view_ctx, void *tcl_data)
{
    return rt_view_context_tclcad_data_set(view_ctx, tcl_data);
}

extern "C" GED_EXPORT int
ged_view_context_callbacks_set(void *view_ctx, struct bu_ptbl *callbacks)
{
    return rt_view_context_callbacks_set(view_ctx, callbacks);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_drawn_count_get(const void *view_ctx)
{
    return rt_view_context_refresh_drawn_count_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_drawn_count_set(void *view_ctx, int count)
{
    return rt_view_context_refresh_drawn_count_set(view_ctx, count);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_request(void *view_ctx, uint32_t flags)
{
    return rt_view_context_refresh_request(view_ctx, flags);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_dirty_get(const void *view_ctx)
{
    return rt_view_context_refresh_dirty_get(view_ctx);
}

extern "C" GED_EXPORT uint32_t
ged_view_context_refresh_consume(void *view_ctx)
{
    return rt_view_context_refresh_consume(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_complete(void *view_ctx)
{
    return rt_view_context_refresh_complete(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_enabled_get(const void *view_ctx)
{
    return rt_view_context_refresh_enabled_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_enabled_set(void *view_ctx, int enabled)
{
    return rt_view_context_refresh_enabled_set(view_ctx, enabled);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_suppressed_get(const void *view_ctx)
{
    return rt_view_context_refresh_suppressed_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_suppress_begin(void *view_ctx)
{
    return rt_view_context_refresh_suppress_begin(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_refresh_suppress_end(void *view_ctx)
{
    return rt_view_context_refresh_suppress_end(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_edit_matrix_set(void *view_ctx, matp_t edit_mat)
{
    return rt_view_context_edit_matrix_set(view_ctx, edit_mat);
}

extern "C" GED_EXPORT int
ged_view_context_edit_matrix_clear(void *view_ctx)
{
    return rt_view_context_edit_matrix_clear(view_ctx);
}

extern "C" GED_EXPORT uint64_t
ged_view_context_frame_revision_get(const void *view_ctx)
{
    return rt_view_context_frame_revision_get(view_ctx);
}

extern "C" GED_EXPORT const char *
ged_view_context_name_get(const void *view_ctx)
{
    return rt_view_context_name_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_name_set(void *view_ctx, const char *name)
{
    return rt_view_context_name_set(view_ctx, name);
}

extern "C" GED_EXPORT void *
ged_view_context_create(void)
{
    return rt_view_context_create();
}

extern "C" GED_EXPORT void *
ged_view_context_create_with_set(void *view_set_ctx)
{
    return rt_view_context_create_with_set(view_set_ctx);
}

extern "C" GED_EXPORT void *
ged_view_context_create_copy_with_set(const void *src_view_ctx, void *view_set_ctx)
{
    return rt_view_context_create_copy_with_set(src_view_ctx, view_set_ctx);
}

extern "C" GED_EXPORT void
ged_view_context_free(void *view_ctx)
{
    rt_view_context_free(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_set_context_add(void *view_set_ctx, void *view_ctx)
{
    return rt_view_set_context_add(view_set_ctx, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_update_callback_set(void *view_ctx,
	ged_view_context_update_callback_t callback,
	void *data)
{
    return rt_view_context_update_callback_set(view_ctx, callback, data);
}

extern "C" GED_EXPORT void *
ged_view_context_display_manager_get(const void *view_ctx)
{
    return rt_view_context_display_manager_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_display_manager_set(void *view_ctx, void *dmp)
{
    if (view_ctx && !dmp) {
	void *old_dmp = rt_view_context_display_manager_get(view_ctx);
	void *old_controller = dm_obol_controller((struct dm *)old_dmp);
	if (old_controller) {
	    struct ged *gedp = (struct ged *)rt_view_context_user_data_get(view_ctx);
	    ged_draw_obol_controller_detach_opaque(gedp, old_controller);
	}
    }
    return rt_view_context_display_manager_set(view_ctx, dmp);
}

extern "C" GED_EXPORT int
ged_view_context_width_get(const void *view_ctx)
{
    return rt_view_context_width_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_height_get(const void *view_ctx)
{
    return rt_view_context_height_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_dimensions_set(void *view_ctx, int width, int height)
{
    return rt_view_context_dimensions_set(view_ctx, width, height);
}

extern "C" GED_EXPORT fastf_t *
ged_view_context_scale_storage_get(void *view_ctx)
{
    return rt_view_context_scale_storage_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_framebuffer_mode_get(const void *view_ctx)
{
    return rt_view_context_framebuffer_mode_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_framebuffer_mode_set(void *view_ctx, int mode)
{
    return rt_view_context_framebuffer_mode_set(view_ctx, mode);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_perspective_get(const void *view_ctx)
{
    return rt_view_context_perspective_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_perspective_set(void *view_ctx, fastf_t perspective)
{
    return rt_view_context_perspective_set(view_ctx, perspective);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_size_get(const void *view_ctx)
{
    return rt_view_context_size_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_size_set(void *view_ctx, fastf_t size)
{
    return rt_view_context_size_set(view_ctx, size);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_inverse_size_get(const void *view_ctx)
{
    return rt_view_context_inverse_size_get(view_ctx);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_local2base_get(const void *view_ctx)
{
    return rt_view_context_local2base_get(view_ctx);
}

extern "C" GED_EXPORT fastf_t
ged_view_context_base2local_get(const void *view_ctx)
{
    return rt_view_context_base2local_get(view_ctx);
}

extern "C" GED_EXPORT void
ged_view_context_eye_pos_get(point_t eye_pos, const void *view_ctx)
{
    rt_view_context_eye_pos_get(eye_pos, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_eye_pos_set(void *view_ctx, const point_t eye_pos)
{
    return rt_view_context_eye_pos_set(view_ctx, eye_pos);
}

extern "C" GED_EXPORT void
ged_view_context_keypoint_get(point_t keypoint, const void *view_ctx)
{
    rt_view_context_keypoint_get(keypoint, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_keypoint_set(void *view_ctx, const point_t keypoint)
{
    return rt_view_context_keypoint_set(view_ctx, keypoint);
}

extern "C" GED_EXPORT char
ged_view_context_rotate_about_get(const void *view_ctx)
{
    return rt_view_context_rotate_about_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_rotate_about_set(void *view_ctx, char rotate_about)
{
    return rt_view_context_rotate_about_set(view_ctx, rotate_about);
}

extern "C" GED_EXPORT char
ged_view_context_coord_get(const void *view_ctx)
{
    return rt_view_context_coord_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_coord_set(void *view_ctx, char coord)
{
    return rt_view_context_coord_set(view_ctx, coord);
}

extern "C" GED_EXPORT int
ged_view_context_zclip_get(const void *view_ctx)
{
    return rt_view_context_zclip_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_zclip_set(void *view_ctx, int zclip)
{
    return rt_view_context_zclip_set(view_ctx, zclip);
}

extern "C" GED_EXPORT int
ged_view_context_autoview(void *view_ctx, fastf_t scale, int all_view_objs)
{
    return rt_view_context_autoview(view_ctx, scale, all_view_objs);
}

extern "C" GED_EXPORT int
ged_view_context_autoview_bounds(void *view_ctx, fastf_t scale, const point_t min, const point_t max)
{
    return rt_view_context_autoview_bounds(view_ctx, scale, min, max);
}

extern "C" GED_EXPORT int
ged_view_context_screen_to_view(fastf_t *fx, fastf_t *fy, void *view_ctx, fastf_t x, fastf_t y)
{
    return rt_view_context_screen_to_view(fx, fy, view_ctx, x, y);
}

extern "C" GED_EXPORT int
ged_view_context_screen_point(point_t point, void *view_ctx, fastf_t x, fastf_t y)
{
    return rt_view_context_screen_point_get(point, view_ctx, x, y);
}

extern "C" GED_EXPORT int
ged_view_context_previous_mouse_get(fastf_t *x, fastf_t *y, const void *view_ctx)
{
    return rt_view_context_previous_mouse_get(x, y, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_previous_mouse_set(void *view_ctx, fastf_t x, fastf_t y)
{
    return rt_view_context_previous_mouse_set(view_ctx, x, y);
}

extern "C" GED_EXPORT int
ged_view_context_mouse_delta_settings_get(struct rt_view_mouse_delta_settings *settings, const void *view_ctx)
{
    return rt_view_context_mouse_delta_settings_get(settings, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_mouse_state_set(void *view_ctx, int x, int y)
{
    return rt_view_context_mouse_state_set(view_ctx, x, y);
}

extern "C" GED_EXPORT void
ged_view_context_model2view_get(mat_t model2view, const void *view_ctx)
{
    rt_view_context_model2view_get(model2view, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_model2view_set(void *view_ctx, const mat_t model2view)
{
    return rt_view_context_model2view_set(view_ctx, model2view);
}

extern "C" GED_EXPORT void
ged_view_context_view2model_get(mat_t view2model, const void *view_ctx)
{
    rt_view_context_view2model_get(view2model, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_view2model_set(void *view_ctx, const mat_t view2model)
{
    return rt_view_context_view2model_set(view_ctx, view2model);
}

extern "C" GED_EXPORT void
ged_view_context_pmodel2view_get(mat_t pmodel2view, const void *view_ctx)
{
    rt_view_context_pmodel2view_get(pmodel2view, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_pmodel2view_set(void *view_ctx, const mat_t pmodel2view)
{
    return rt_view_context_pmodel2view_set(view_ctx, pmodel2view);
}

extern "C" GED_EXPORT void
ged_view_context_pmat_get(mat_t pmat, const void *view_ctx)
{
    rt_view_context_pmat_get(pmat, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_pmat_set(void *view_ctx, const mat_t pmat)
{
    return rt_view_context_pmat_set(view_ctx, pmat);
}

extern "C" GED_EXPORT void
ged_view_context_center_get(mat_t center, const void *view_ctx)
{
    rt_view_context_center_get(center, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_center_vec_set(void *view_ctx, const point_t center)
{
    return rt_view_context_center_set(view_ctx, center);
}

extern "C" GED_EXPORT void
ged_view_context_rotation_get(mat_t rotation, const void *view_ctx)
{
    rt_view_context_rotation_get(rotation, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_rotation_set(void *view_ctx, const mat_t rotation)
{
    return rt_view_context_rotation_set(view_ctx, rotation);
}

extern "C" GED_EXPORT void
ged_view_context_aet_get(vect_t aet, const void *view_ctx)
{
    rt_view_context_aet_get(aet, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_aet_set(void *view_ctx, const vect_t aet)
{
    return rt_view_context_aet_set(view_ctx, aet);
}

extern "C" GED_EXPORT int
ged_view_context_aet_state_set(void *view_ctx, const vect_t aet)
{
    return rt_view_context_aet_state_set(view_ctx, aet);
}

extern "C" GED_EXPORT int
ged_view_context_plane_get(plane_t *plane, const void *view_ctx)
{
    return rt_view_context_plane_get(plane, view_ctx);
}

extern "C" GED_EXPORT void
ged_view_context_info_get(struct rt_view_info *info, const void *view_ctx)
{
    rt_view_context_info_get(info, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_orientation_quat_get(quat_t orientation, const void *view_ctx)
{
    return rt_view_context_orientation_quat_get(orientation, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_adc_state_get(struct rt_view_adc_state *adc, const void *view_ctx)
{
    return rt_view_context_adc_state_get(adc, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_adc_state_set(void *view_ctx, const struct rt_view_adc_state *adc)
{
    return rt_view_context_adc_state_set(view_ctx, adc);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_state_get(struct rt_view_knobs *knobs, const void *view_ctx)
{
    return rt_view_context_knobs_state_get(knobs, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_state_set(void *view_ctx, const struct rt_view_knobs *knobs)
{
    return rt_view_context_knobs_state_set(view_ctx, knobs);
}

extern "C" GED_EXPORT unsigned long long
ged_view_context_knobs_hash(void *view_ctx, struct bu_data_hash_state *state)
{
    return rt_view_context_knobs_hash(view_ctx, state);
}

extern "C" GED_EXPORT int
ged_view_context_knob_values_get(struct rt_view_knob_values *values, const void *view_ctx)
{
    return rt_view_context_knob_values_get(values, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_reset(void *view_ctx, int category)
{
    return rt_view_context_knobs_reset(view_ctx, category);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_calibrate(void *view_ctx)
{
    return rt_view_context_knobs_calibrate(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_cmd_process(vect_t *rvec,
	int *do_rot,
	vect_t *tvec,
	int *do_tran,
	void *view_ctx,
	const char *cmd,
	fastf_t factor,
	char origin,
	int model_flag,
	int incr_flag)
{
    return rt_view_context_knobs_cmd_process(rvec, do_rot, tvec, do_tran,
	    view_ctx, cmd, factor, origin, model_flag, incr_flag);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_translate(void *view_ctx, const vect_t tvec, int model_flag)
{
    return rt_view_context_knobs_translate(view_ctx, tvec, model_flag);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_rotate(void *view_ctx,
	const vect_t rvec,
	char origin,
	char coords,
	const matp_t obj_rot,
	const pointp_t pvt_pt)
{
    return rt_view_context_knobs_rotate(view_ctx, rvec, origin, coords,
	    obj_rot, pvt_pt);
}

extern "C" GED_EXPORT int
ged_view_context_knobs_update_rate_flags(void *view_ctx)
{
    return rt_view_context_knobs_update_rate_flags(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_grid_state_get(struct rt_view_grid_state *grid, const void *view_ctx)
{
    return rt_view_context_grid_state_get(grid, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_grid_state_set(void *view_ctx, const struct rt_view_grid_state *grid)
{
    return rt_view_context_grid_state_set(view_ctx, grid);
}

extern "C" GED_EXPORT int
ged_view_context_model_axes_state_get(struct rt_view_axes_state *axes, const void *view_ctx)
{
    return rt_view_context_model_axes_state_get(axes, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_model_axes_state_set(void *view_ctx, const struct rt_view_axes_state *axes)
{
    return rt_view_context_model_axes_state_set(view_ctx, axes);
}

extern "C" GED_EXPORT int
ged_view_context_view_axes_state_get(struct rt_view_axes_state *axes, const void *view_ctx)
{
    return rt_view_context_view_axes_state_get(axes, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_view_axes_state_set(void *view_ctx, const struct rt_view_axes_state *axes)
{
    return rt_view_context_view_axes_state_set(view_ctx, axes);
}

extern "C" GED_EXPORT int
ged_view_context_center_dot_state_get(struct rt_view_other_state *center_dot, const void *view_ctx)
{
    return rt_view_context_center_dot_state_get(center_dot, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_center_dot_state_set(void *view_ctx, const struct rt_view_other_state *center_dot)
{
    return rt_view_context_center_dot_state_set(view_ctx, center_dot);
}

extern "C" GED_EXPORT int
ged_view_context_scale_overlay_state_get(struct rt_view_other_state *scale_state, const void *view_ctx)
{
    return rt_view_context_scale_overlay_state_get(scale_state, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_scale_overlay_state_set(void *view_ctx, const struct rt_view_other_state *scale_state)
{
    return rt_view_context_scale_overlay_state_set(view_ctx, scale_state);
}

extern "C" GED_EXPORT int
ged_view_context_params_state_get(struct rt_view_params_state *params, const void *view_ctx)
{
    return rt_view_context_params_state_get(params, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_params_state_set(void *view_ctx, const struct rt_view_params_state *params)
{
    return rt_view_context_params_state_set(view_ctx, params);
}

extern "C" GED_EXPORT int
ged_view_context_snap_grid_2d(void *view_ctx, fastf_t *vx, fastf_t *vy)
{
    return rt_view_context_snap_grid_2d(view_ctx, vx, vy);
}

extern "C" GED_EXPORT unsigned long long
ged_view_context_prepare_tcl_snap(void *view_ctx)
{
    return rt_view_context_prepare_tcl_snap(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_snap_point_2d(void *view_ctx,
	fastf_t *vx,
	fastf_t *vy,
	unsigned long long kinds)
{
    return rt_view_context_snap_point_2d(view_ctx, vx, vy, kinds);
}

extern "C" GED_EXPORT int
ged_view_context_interactive_rect_state_get(struct rt_view_interactive_rect_state *rect,
	const void *view_ctx)
{
    return rt_view_context_interactive_rect_state_get(rect, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_interactive_rect_state_set(void *view_ctx,
	const struct rt_view_interactive_rect_state *rect)
{
    return rt_view_context_interactive_rect_state_set(view_ctx, rect);
}

extern "C" GED_EXPORT int
ged_view_context_snap_lines_get(const void *view_ctx)
{
    return rt_view_context_snap_lines_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_snap_lines_set(void *view_ctx, int enabled)
{
    return rt_view_context_snap_lines_set(view_ctx, enabled);
}

extern "C" GED_EXPORT int
ged_view_context_center_linesnap(void *view_ctx)
{
    return rt_view_context_center_linesnap(view_ctx);
}

extern "C" GED_EXPORT double
ged_view_context_snap_tolerance_factor_get(const void *view_ctx)
{
    return rt_view_context_snap_tolerance_factor_get(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_snap_tolerance_factor_set(void *view_ctx, double factor)
{
    return rt_view_context_snap_tolerance_factor_set(view_ctx, factor);
}

extern "C" GED_EXPORT int
ged_view_context_unit_conversion_set(void *view_ctx, fastf_t local2base, fastf_t base2local)
{
    return rt_view_context_unit_conversion_set(view_ctx, local2base, base2local);
}

extern "C" GED_EXPORT int
ged_view_context_lod_policy_get(struct rt_view_lod_policy *policy, const void *view_ctx)
{
    return rt_view_context_lod_policy_get(policy, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_lod_policy_apply(void *view_ctx, const struct rt_view_lod_policy *policy)
{
    return rt_view_context_lod_policy_apply(view_ctx, policy);
}

extern "C" GED_EXPORT int
ged_view_context_settings_shared(const void *a, const void *b)
{
    return rt_view_context_settings_shared(a, b);
}
