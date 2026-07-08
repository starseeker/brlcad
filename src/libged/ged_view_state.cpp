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
 * GED-owned opaque view-context state and facade calls.
 *
 * Generic view mechanics are intentionally not mirrored here.  While libbv
 * ownership is being finalized, callers use the lower-level view API directly
 * for scale, matrices, mouse state, knobs, and passive overlay data.  This file
 * owns only GED view registration plus host/policy operations that still have
 * GED behavior attached.
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
ged_view_context_edit_matrix_set(void *view_ctx, matp_t edit_mat)
{
    return rt_view_context_edit_matrix_set(view_ctx, edit_mat);
}

extern "C" GED_EXPORT int
ged_view_context_edit_matrix_clear(void *view_ctx)
{
    return rt_view_context_edit_matrix_clear(view_ctx);
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
