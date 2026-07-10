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

#include "brlobol/view_attachment.h"
#include "bv.h"
#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "ged/draw_obol.h"
#include "ged/view.h"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

struct ged_view_state_storage {
    void *active_view_ctx;
    struct bv_context_set *view_set_ctx;
    struct bu_ptbl free_views;
    struct bu_ptbl host_records;
};

struct ged_view_host_record {
    void *view_ctx;
    struct ged *gedp;
    void *display_manager;
    void *tclcad_data;
    struct bu_ptbl *callbacks;
    ged_view_context_update_callback_t update_callback;
    void *update_callback_data;
    BRLObolViewAttachment *obol_attachment;
};

static struct bu_ptbl ged_view_host_registry;
static int ged_view_host_registry_initialized = 0;

static void
ged_view_host_registry_init(void)
{
    if (ged_view_host_registry_initialized)
	return;

    BU_PTBL_INIT(&ged_view_host_registry);
    ged_view_host_registry_initialized = 1;
}

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
    state->view_set_ctx = bv_context_set_create();
    BU_PTBL_INIT(&state->free_views);
    BU_PTBL_INIT(&state->host_records);
    impl->ged_view_state_ctx = (void *)state;
    return state;
}

static void
ged_view_host_callbacks_free(struct bu_ptbl *callbacks)
{
    if (!callbacks)
	return;

    bu_ptbl_free(callbacks);
    BU_PUT(callbacks, struct bu_ptbl);
}

static struct ged_view_host_record *
ged_view_host_record_find_global(const void *view_ctx)
{
    if (!view_ctx || !ged_view_host_registry_initialized)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&ged_view_host_registry); i++) {
	struct ged_view_host_record *record =
	    (struct ged_view_host_record *)BU_PTBL_GET(&ged_view_host_registry, i);
	if (record && record->view_ctx == view_ctx)
	    return record;
    }

    return NULL;
}

static struct ged_view_host_record *
ged_view_host_record_create(struct ged *gedp, void *view_ctx)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    struct ged_view_host_record *record = NULL;

    if (!state || !view_ctx)
	return NULL;

    record = ged_view_host_record_find_global(view_ctx);
    if (record) {
	if (record->gedp && record->gedp != gedp) {
	    struct ged_view_state_storage *old_state =
		ged_view_state(record->gedp);
	    if (old_state)
		bu_ptbl_rm(&old_state->host_records, (long *)record);
	}
	record->gedp = gedp;
	bu_ptbl_ins_unique(&state->host_records, (long *)record);
	return record;
    }

    record = (struct ged_view_host_record *)bu_calloc(1, sizeof(*record),
	    "GED view host record");
    record->view_ctx = view_ctx;
    record->gedp = gedp;
    record->obol_attachment = new BRLObolViewAttachment;
    record->obol_attachment->ref();

    ged_view_host_registry_init();
    bu_ptbl_ins(&ged_view_host_registry, (long *)record);
    bu_ptbl_ins_unique(&state->host_records, (long *)record);
    return record;
}

static void
ged_view_host_record_destroy(struct ged_view_host_record *record)
{
    if (!record)
	return;

    if (record->gedp) {
	struct ged_view_state_storage *state = ged_view_state(record->gedp);
	if (state)
	    bu_ptbl_rm(&state->host_records, (long *)record);
    }

    if (ged_view_host_registry_initialized) {
	bu_ptbl_rm(&ged_view_host_registry, (long *)record);
	if (!BU_PTBL_LEN(&ged_view_host_registry)) {
	    bu_ptbl_free(&ged_view_host_registry);
	    ged_view_host_registry_initialized = 0;
	}
    }

    ged_view_host_callbacks_free(record->callbacks);
    record->obol_attachment->unref();
    record->obol_attachment = NULL;
    bu_free(record, "GED view host record");
}

static void
ged_view_host_record_destroy_for_view(void *view_ctx)
{
    ged_view_host_record_destroy(ged_view_host_record_find_global(view_ctx));
}

extern "C" void
ged_view_state_init(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state_create(gedp);
    if (!state || !state->view_set_ctx)
	return;

    void *default_view = ged_view_context_create_with_set(state->view_set_ctx);
    state->active_view_ctx = default_view;
    bv_name_set(bv_context_view((struct bv_context *)default_view), "default");
    ged_view_set_context_add(state->view_set_ctx, default_view);
    ged_view_context_owned_add(gedp, default_view);
}

extern "C" void
ged_view_state_free(struct ged *gedp)
{
    struct ged_impl *impl = ged_view_impl(gedp);
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!impl || !state)
	return;

    state->active_view_ctx = NULL;

    while (BU_PTBL_LEN(&state->free_views)) {
	void *view_ctx = (void *)BU_PTBL_GET(&state->free_views, 0);
	bu_ptbl_rm(&state->free_views, (long *)view_ctx);
	ged_view_host_record_destroy_for_view(view_ctx);
	bv_context_destroy((struct bv_context *)view_ctx);
    }
    bu_ptbl_free(&state->free_views);
    while (BU_PTBL_LEN(&state->host_records)) {
	struct ged_view_host_record *record =
	    (struct ged_view_host_record *)BU_PTBL_GET(&state->host_records, 0);
	ged_view_host_record_destroy(record);
    }
    bu_ptbl_free(&state->host_records);
    bv_context_set_destroy(state->view_set_ctx);
    impl->ged_view_state_ctx = NULL;
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

    bu_ptbl_ins_unique(&state->free_views, (long *)view_ctx);
    return ged_view_context_host_attach(gedp, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_host_attach(struct ged *gedp, void *view_ctx)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !view_ctx)
	return 0;

    (void)ged_view_host_record_create(gedp, view_ctx);
    struct bv_context *ctx = (struct bv_context *)view_ctx;
    if (bv_context_view(ctx) && !bv_context_user_data_get(ctx))
	bv_context_user_data_set(ctx, gedp);
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
    return state ? bv_context_set_views(state->view_set_ctx) : NULL;
}

extern "C" GED_EXPORT void *
ged_view_find_ctx(struct ged *gedp, const char *name)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !name)
	return NULL;

    return bv_context_set_find(state->view_set_ctx, name);
}

BRLObolViewAttachment *
ged_view_context_obol_attachment(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? record->obol_attachment : NULL;
}

int
ged_view_context_obol_attachment_bind(void *view_ctx,
				      BRLObolViewAttachment *attachment)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (!record || !attachment)
	return 0;

    if (record->obol_attachment == attachment)
	return 1;

    attachment->copyHostStateFrom(record->obol_attachment);
    attachment->ref();
    record->obol_attachment->unref();
    record->obol_attachment = attachment;
    return 1;
}

int
ged_view_context_obol_attachment_unbind(void *view_ctx,
					BRLObolViewAttachment *attachment)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (!record || !attachment)
	return 0;
    if (record->obol_attachment != attachment)
	return 1;

    BRLObolViewAttachment *replacement = new BRLObolViewAttachment;
    replacement->ref();
    replacement->copyHostStateFrom(record->obol_attachment);
    record->obol_attachment->unref();
    record->obol_attachment = replacement;
    return 1;
}

extern "C" GED_EXPORT void
ged_draw_view_context_info_get(struct bv_view_info *view_info,
			       const void *view_ctx)
{
    if (!view_info)
	return;

    bv_view_info_init(view_info);
    const struct bv *view =
	bv_context_view_const((const struct bv_context *)view_ctx);
    if (!view)
	return;

    view_info->width = bv_width_get(view);
    view_info->height = bv_height_get(view);
    view_info->size = bv_size_get(view);
    bv_view_info_sanitize(view_info);
}

extern "C" GED_EXPORT int
ged_draw_view_context_lod_policy_get(ged_draw_view_lod_policy *policy,
				     const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	if (!policy)
	    return 0;
	record->obol_attachment->getLodPolicy(policy);
	return 1;
    }

    if (!policy || !view_ctx)
	return 0;
    bv_lod_policy_init(policy);
    return 1;
}

extern "C" GED_EXPORT int
ged_draw_view_context_lod_policy_apply(void *view_ctx,
				       const ged_draw_view_lod_policy *policy)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	if (!policy)
	    return 0;
	record->obol_attachment->setLodPolicy(policy);
	return 1;
    }

    return (view_ctx && policy) ? 1 : 0;
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

extern "C" GED_EXPORT int
ged_draw_view_context_lod_bounds_update(void *view_ctx)
{
    return view_ctx ? 1 : 0;
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

    if (!view_ctx || !sample)
	return 0;

    return ged_draw_obol_view_context_snap_first_candidate(view_ctx, sample,
	    kind, candidate);
}

extern "C" GED_EXPORT void
ged_draw_view_context_lod_bounds_callback_set(void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	record->obol_attachment->setLodBoundsCallbackSet(TRUE);
	return;
    }

    (void)view_ctx;
}

extern "C" GED_EXPORT int
ged_draw_view_context_lod_bounds_callback_is(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? (record->obol_attachment->isLodBoundsCallbackSet() ? 1 : 0) : 0;
}

extern "C" GED_EXPORT int
ged_view_context_is_independent(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? (record->obol_attachment->isIndependentScopeCreated() ? 1 : 0) : 0;
}

extern "C" GED_EXPORT int
ged_view_context_independent_scope_is_null(void *view_ctx, int create)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	if (create)
	    record->obol_attachment->setIndependentScopeCreated(TRUE);
	return record->obol_attachment->isIndependentScopeCreated() ? 0 : 1;
    }

    (void)view_ctx;
    (void)create;
    return 1;
}

extern "C" GED_EXPORT void
ged_view_context_independent_scope_destroy(void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	record->obol_attachment->setIndependentScopeCreated(FALSE);
	return;
    }

    (void)view_ctx;
}

extern "C" GED_EXPORT ged_draw_scene_handle
ged_view_context_scene_root_ref(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    (void)view_ctx;
    if (!record || !record->obol_attachment->hasSceneRootToken())
	return ged_draw_scene_handle_null();

    return ged_draw_scene_handle_make(
	    record->obol_attachment->getSceneRootToken(),
	    GED_DRAW_SCENE_BACKEND_OBOL);
}

extern "C" GED_EXPORT int
ged_view_context_scene_root_ref_attach(void *view_ctx,
				       ged_draw_scene_handle root_ref)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	if (!ged_draw_scene_handle_is_null(root_ref) &&
		ged_draw_scene_handle_backend(root_ref) != GED_DRAW_SCENE_BACKEND_OBOL)
	    return 0;
	record->obol_attachment->setSceneRootToken(
		ged_draw_scene_handle_context(root_ref));
	return !ged_draw_scene_handle_is_null(root_ref) ? 1 : 0;
    }

    (void)view_ctx;
    (void)root_ref;
    return 0;
}

extern "C" GED_EXPORT int
ged_view_context_scene_attached(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    (void)view_ctx;
    return record ? (record->obol_attachment->hasSceneRootToken() ? 1 : 0) : 0;
}

extern "C" GED_EXPORT size_t
ged_view_context_clear(void *view_ctx, int flags)
{
    size_t cleared = ged_draw_obol_view_context_clear(view_ctx, flags);
    (void)bv_context_refresh_complete((struct bv_context *)view_ctx);
    return cleared;
}

extern "C" GED_EXPORT void *
ged_view_context_user_data_get(const void *view_ctx)
{
    return bv_context_user_data_get((const struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_user_data_set(void *view_ctx, void *user_data)
{
    return bv_context_user_data_set((struct bv_context *)view_ctx,
	    user_data);
}

extern "C" GED_EXPORT void *
ged_view_context_tclcad_data_get(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? record->tclcad_data : NULL;
}

extern "C" GED_EXPORT int
ged_view_context_tclcad_data_set(void *view_ctx, void *tcl_data)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	record->tclcad_data = tcl_data;
	return 1;
    }

    return 0;
}

extern "C" GED_EXPORT int
ged_view_context_callbacks_set(void *view_ctx, struct bu_ptbl *callbacks)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	if (record->callbacks && record->callbacks != callbacks)
	    ged_view_host_callbacks_free(record->callbacks);
	record->callbacks = callbacks;
	return 1;
    }

    return 0;
}

extern "C" GED_EXPORT void *
ged_view_context_create(void)
{
    return bv_context_create();
}

extern "C" GED_EXPORT void *
ged_view_context_create_with_set(void *view_set_ctx)
{
    struct bv_context *ctx = bv_context_create();
    if (!ctx)
	return NULL;
    if (view_set_ctx &&
	    !bv_context_set_attach((struct bv_context_set *)view_set_ctx, ctx)) {
	bv_context_destroy(ctx);
	return NULL;
    }
    return ctx;
}

extern "C" GED_EXPORT void *
ged_view_context_create_copy_with_set(const void *src_view_ctx, void *view_set_ctx)
{
    struct bv_context *ctx = bv_context_create();
    if (!ctx)
	return NULL;

    const struct bv *src_view =
	bv_context_view_const((const struct bv_context *)src_view_ctx);
    if (src_view)
	bv_copy(bv_context_view(ctx), src_view);
    if (view_set_ctx &&
	    !bv_context_set_attach((struct bv_context_set *)view_set_ctx, ctx)) {
	bv_context_destroy(ctx);
	return NULL;
    }
    return ctx;
}

extern "C" GED_EXPORT void
ged_view_context_free(void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record && record->gedp) {
	struct ged_view_state_storage *state = ged_view_state(record->gedp);
	if (state)
	    bu_ptbl_rm(&state->free_views, (long *)view_ctx);
    }

    ged_view_host_record_destroy(record);
    bv_context_destroy((struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_set_context_add(void *view_set_ctx, void *view_ctx)
{
    return bv_context_set_add((struct bv_context_set *)view_set_ctx,
	    (struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_set_context_remove(void *view_set_ctx, void *view_ctx)
{
    return bv_context_set_remove((struct bv_context_set *)view_set_ctx,
	    (struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_view_set_attach(void *view_ctx, void *view_set_ctx)
{
    return bv_context_set_attach((struct bv_context_set *)view_set_ctx,
	    (struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_update_callback_set(void *view_ctx,
				     ged_view_context_update_callback_t callback,
				     void *data)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	record->update_callback = callback;
	record->update_callback_data = data;
	return 1;
    }

    return 0;
}

extern "C" GED_EXPORT void *
ged_view_context_display_manager_get(const void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? record->display_manager : NULL;
}

extern "C" GED_EXPORT int
ged_view_context_display_manager_set(void *view_ctx, void *dmp)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	void *old_dmp = record->display_manager;
	if (old_dmp && old_dmp != dmp)
	    ged_draw_obol_controller_detach_for_view(record->gedp,
		    view_ctx);
	record->display_manager = dmp;
	return 1;
    }

    return 0;
}

extern "C" GED_EXPORT int
ged_view_context_update(void *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);
    int ret = bv_context_update((struct bv_context *)view_ctx,
	    BV_CONTEXT_CHANGED_VIEW);

    if (ret && record && record->update_callback)
	record->update_callback(view_ctx, record->update_callback_data);

    return ret;
}
