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

#include <cmath>
#include <mutex>
#include <thread>

#include "bn/mat.h"
#include "BObol/BDisplayEndpoint.h"
#include "BObol/BViewController.h"
#include "BObol/BViewAttachment.h"
#include "bv.h"
#include "bu/malloc.h"
#include "bu/assert.h"
#include "bu/ptbl.h"
#include "bu/str.h"
#include "ged/display.h"
#include "ged/display_obol_private.h"
#include "ged/view.h"
#include "./ged_bobol_private.hpp"
#include "./ged_draw_private.h"
#include "./ged_draw_view_private.h"
#include "./ged_private.h"

struct ged_view_state_storage {
    struct ged_view_context *active_view_ctx;
    struct ged_view_set *view_set_ctx;
    struct bu_ptbl free_views;
    struct bu_ptbl host_records;
};

struct ged_view_host_record {
    ged_view_host_record(void) :
	identity(0),
	owner_thread(std::this_thread::get_id()),
	view_ctx(NULL),
	gedp(NULL),
	display_endpoint(NULL),
	owns_display_endpoint(0),
	tclcad_data(NULL),
	callbacks(NULL),
	update_callback(NULL),
	update_callback_data(NULL),
	obol_attachment(NULL)
    {
    }

    uint64_t identity;
    std::thread::id owner_thread;
    struct ged_view_context *view_ctx;
    struct ged *gedp;
    bobol_display_endpoint_t *display_endpoint;
    int owns_display_endpoint;
    void *tclcad_data;
    struct bu_ptbl *callbacks;
    ged_view_context_update_callback_t update_callback;
    void *update_callback_data;
    BObolViewAttachment *obol_attachment;
};

/* Non-owning live-record index used only to resolve value-handle owner ids.
 * Per-GED host_records own the records. */
static struct bu_ptbl ged_view_host_registry;
static int ged_view_host_registry_initialized = 0;
static uint64_t ged_view_host_next_identity = 1;
static std::mutex ged_view_host_registry_mutex;

static void
ged_view_host_registry_init_locked(void)
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
    state->view_set_ctx = (struct ged_view_set *)bv_context_set_create();
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
ged_view_host_record_find_global(const struct ged_view_context *view_ctx)
{
    std::lock_guard<std::mutex> guard(ged_view_host_registry_mutex);
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
ged_view_host_record_find_identity(uint64_t identity)
{
    std::lock_guard<std::mutex> guard(ged_view_host_registry_mutex);
    if (!identity || !ged_view_host_registry_initialized)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&ged_view_host_registry); i++) {
	struct ged_view_host_record *record =
	    (struct ged_view_host_record *)BU_PTBL_GET(&ged_view_host_registry, i);
	if (record && record->identity == identity)
	    return record;
    }

    return NULL;
}

static struct ged_view_host_record *
ged_view_host_record_create(struct ged *gedp, struct ged_view_context *view_ctx)
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

    record = new ged_view_host_record;
    record->view_ctx = view_ctx;
    record->gedp = gedp;
    record->obol_attachment = new BObolViewAttachment;
    record->obol_attachment->ref();

    {
	std::lock_guard<std::mutex> guard(ged_view_host_registry_mutex);
	ged_view_host_registry_init_locked();
	record->identity = ged_view_host_next_identity++;
	if (!record->identity)
	    record->identity = ged_view_host_next_identity++;
	bu_ptbl_ins(&ged_view_host_registry, (long *)record);
    }
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

    {
	std::lock_guard<std::mutex> guard(ged_view_host_registry_mutex);
	if (ged_view_host_registry_initialized) {
	    bu_ptbl_rm(&ged_view_host_registry, (long *)record);
	    if (!BU_PTBL_LEN(&ged_view_host_registry)) {
		bu_ptbl_free(&ged_view_host_registry);
		ged_view_host_registry_initialized = 0;
	    }
	}
    }

    if (record->display_endpoint) {
	(void)bobol_display_endpoint_property_provider_set(
	    record->display_endpoint, NULL, NULL, NULL);
	/* The framebuffer bridge owns retained nodes in the endpoint host.  Close
	 * it while both the host and controller are still valid. */
	if (record->gedp)
	    ged_view_framebuffer_release(record->gedp);
	/* Clear the factory callback before GED releases a controller it owns. */
	bobol_display_endpoint_host_detach(record->display_endpoint);
	ged_draw_obol_controller_detach_for_view(record->gedp,
		record->view_ctx);
	if (record->owns_display_endpoint)
	    bobol_display_endpoint_destroy(record->display_endpoint);
	record->display_endpoint = NULL;
	record->owns_display_endpoint = 0;
    }
    ged_view_host_callbacks_free(record->callbacks);
    record->obol_attachment->unref();
    record->obol_attachment = NULL;
    delete record;
}

static void
ged_view_host_record_destroy_for_view(struct ged_view_context *view_ctx)
{
    ged_view_host_record_destroy(ged_view_host_record_find_global(view_ctx));
}

extern "C" void
ged_view_state_init(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state_create(gedp);
    if (!state || !state->view_set_ctx)
	return;

    struct ged_view_context *default_view =
	ged_view_context_create_with_set(state->view_set_ctx);
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
	struct ged_view_context *view_ctx =
	    (struct ged_view_context *)BU_PTBL_GET(&state->free_views, 0);
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
    bv_context_set_destroy((struct bv_context_set *)state->view_set_ctx);
    impl->ged_view_state_ctx = NULL;
    bu_free(state, "GED view state");
}

extern "C" GED_EXPORT struct ged_view_context *
ged_view_active_ctx(const struct ged *gedp)
{
    const struct ged_view_state_storage *state = ged_view_state_const(gedp);
    return state ? state->active_view_ctx : NULL;
}

extern "C" GED_EXPORT void
ged_view_active_ctx_set(struct ged *gedp, struct ged_view_context *view_ctx)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (state)
	state->active_view_ctx = view_ctx;
}

extern "C" GED_EXPORT int
ged_view_context_owned_add(struct ged *gedp, struct ged_view_context *view_ctx)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !view_ctx)
	return 0;

    bu_ptbl_ins_unique(&state->free_views, (long *)view_ctx);
    return ged_view_context_host_attach(gedp, view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_host_attach(struct ged *gedp, struct ged_view_context *view_ctx)
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

extern "C" GED_EXPORT struct ged_view_set *
ged_view_set_ctx(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    return state ? state->view_set_ctx : NULL;
}

extern "C" GED_EXPORT struct bu_ptbl *
ged_view_set_views_ctx(struct ged *gedp)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    return state ? bv_context_set_views((struct bv_context_set *)state->view_set_ctx) : NULL;
}

extern "C" GED_EXPORT struct ged_view_context *
ged_view_find_ctx(struct ged *gedp, const char *name)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !name)
	return NULL;

    return (struct ged_view_context *)bv_context_set_find(
	    (struct bv_context_set *)state->view_set_ctx, name);
}

BObolViewAttachment *
ged_view_context_obol_attachment(const struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? record->obol_attachment : NULL;
}

uint64_t
ged_view_context_reference_owner(const struct ged_view_context *view_ctx, int local)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (!record || record->identity > (UINT64_MAX >> 1))
	return 0;
    const bool owner_thread =
	record->owner_thread == std::this_thread::get_id();
    BU_ASSERT(owner_thread);
    if (!owner_thread)
	return 0;
    return (record->identity << 1) | (local ? 1u : 0u);
}

struct ged_view_context *
ged_view_context_from_reference_owner(uint64_t owner, int *local)
{
    if (local)
	*local = 0;
    if (!owner)
	return NULL;

    struct ged_view_host_record *record =
	ged_view_host_record_find_identity(owner >> 1);
    if (!record)
	return NULL;
    const bool owner_thread =
	record->owner_thread == std::this_thread::get_id();
    BU_ASSERT(owner_thread);
    if (!owner_thread)
	return NULL;
    if (local)
	*local = (owner & 1u) ? 1 : 0;
    return record->view_ctx;
}

int
ged_view_context_obol_attachment_bind(struct ged_view_context *view_ctx,
				      BObolViewAttachment *attachment)
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
ged_view_context_obol_attachment_unbind(struct ged_view_context *view_ctx,
					BObolViewAttachment *attachment)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (!record || !attachment)
	return 0;
    if (record->obol_attachment != attachment)
	return 1;

    BObolViewAttachment *replacement = new BObolViewAttachment;
    replacement->ref();
    replacement->copyHostStateFrom(record->obol_attachment);
    record->obol_attachment->unref();
    record->obol_attachment = replacement;
    return 1;
}

extern "C" GED_EXPORT void
ged_view_feature_info_get(struct bv_view_info *view_info,
			       const struct ged_view_context *view_ctx)
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
ged_view_lod_policy_get(ged_view_lod_policy *policy,
				     const struct ged_view_context *view_ctx)
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
ged_view_lod_policy_apply(struct ged_view_context *view_ctx,
			  const ged_view_lod_policy *policy)
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
ged_view_lod_policy_apply_bot_threshold(
	struct ged_view_context *view_ctx,
	const ged_view_lod_policy *policy,
	size_t bot_threshold)
{
    if (!policy)
	return 0;

    ged_view_lod_policy override_policy = *policy;
    override_policy.bot_threshold = bot_threshold;
    return ged_view_lod_policy_apply(view_ctx, &override_policy);
}

extern "C" GED_EXPORT int
ged_view_lod_bounds_update(struct ged_view_context *view_ctx)
{
    return view_ctx ? 1 : 0;
}

extern "C" GED_EXPORT int
ged_view_selection_snap(struct ged_view_context *view_ctx,
					   const point_t sample,
					   enum ged_selection_snap_kind kind,
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
ged_draw_source_lod_bounds_callback_set(struct ged_view_context *view_ctx)
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
ged_draw_source_lod_bounds_callback_is(const struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? (record->obol_attachment->isLodBoundsCallbackSet() ? 1 : 0) : 0;
}

extern "C" GED_EXPORT int
ged_view_context_is_independent(const struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? (record->obol_attachment->isIndependentScopeCreated() ? 1 : 0) : 0;
}

extern "C" GED_EXPORT int
ged_view_context_independent_scope_is_null(struct ged_view_context *view_ctx, int create)
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
ged_view_context_independent_scope_destroy(struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    if (record) {
	record->obol_attachment->setIndependentScopeCreated(FALSE);
	return;
    }

    (void)view_ctx;
}

extern "C" ged_draw_scene_handle
ged_view_context_scene_root_ref(const struct ged_view_context *view_ctx)
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

extern "C" int
ged_view_context_scene_root_ref_attach(struct ged_view_context *view_ctx,
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
ged_view_context_scene_attached(const struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    (void)view_ctx;
    return record ? (record->obol_attachment->hasSceneRootToken() ? 1 : 0) : 0;
}

extern "C" GED_EXPORT size_t
ged_view_context_clear(struct ged_view_context *view_ctx, int flags)
{
    size_t cleared = ged_draw_obol_view_context_clear(view_ctx, flags);
    (void)bv_context_refresh_complete((struct bv_context *)view_ctx);
    return cleared;
}

extern "C" GED_EXPORT void *
ged_view_context_user_data_get(const struct ged_view_context *view_ctx)
{
    return bv_context_user_data_get((const struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_user_data_set(struct ged_view_context *view_ctx, void *user_data)
{
    return bv_context_user_data_set((struct bv_context *)view_ctx,
	    user_data);
}

extern "C" GED_EXPORT void *
ged_view_context_tclcad_data_get(const struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ? record->tclcad_data : NULL;
}

extern "C" GED_EXPORT int
ged_view_context_tclcad_data_set(struct ged_view_context *view_ctx, void *tcl_data)
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
ged_view_context_callbacks_set(struct ged_view_context *view_ctx, struct bu_ptbl *callbacks)
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

extern "C" GED_EXPORT struct ged_view_context *
ged_view_context_create(void)
{
    return (struct ged_view_context *)bv_context_create();
}

extern "C" GED_EXPORT struct ged_view_context *
ged_view_context_create_with_set(struct ged_view_set *view_set_ctx)
{
    struct bv_context *ctx = bv_context_create();
    if (!ctx)
	return NULL;
    if (view_set_ctx &&
	    !bv_context_set_attach((struct bv_context_set *)view_set_ctx, ctx)) {
	bv_context_destroy(ctx);
	return NULL;
    }
    return (struct ged_view_context *)ctx;
}

extern "C" GED_EXPORT struct ged_view_context *
ged_view_context_create_copy_with_set(const struct ged_view_context *src_view_ctx,
	struct ged_view_set *view_set_ctx)
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
    return (struct ged_view_context *)ctx;
}

extern "C" GED_EXPORT void
ged_view_context_free(struct ged_view_context *view_ctx)
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

extern "C" GED_EXPORT struct bv_context *
ged_view_context_bv(struct ged_view_context *view_ctx)
{
    return reinterpret_cast<struct bv_context *>(view_ctx);
}

extern "C" GED_EXPORT const struct bv_context *
ged_view_context_bv_const(const struct ged_view_context *view_ctx)
{
    return reinterpret_cast<const struct bv_context *>(view_ctx);
}

extern "C" GED_EXPORT struct ged_view_context *
ged_view_context_from_bv(struct bv_context *view_ctx)
{
    return reinterpret_cast<struct ged_view_context *>(view_ctx);
}

extern "C" GED_EXPORT const struct ged_view_context *
ged_view_context_from_bv_const(const struct bv_context *view_ctx)
{
    return reinterpret_cast<const struct ged_view_context *>(view_ctx);
}

extern "C" GED_EXPORT int
ged_view_set_context_add(struct ged_view_set *view_set_ctx,
	struct ged_view_context *view_ctx)
{
    return bv_context_set_add((struct bv_context_set *)view_set_ctx,
	    (struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_set_context_remove(struct ged_view_set *view_set_ctx,
	struct ged_view_context *view_ctx)
{
    return bv_context_set_remove((struct bv_context_set *)view_set_ctx,
	    (struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_view_set_attach(struct ged_view_context *view_ctx,
	struct ged_view_set *view_set_ctx)
{
    return bv_context_set_attach((struct bv_context_set *)view_set_ctx,
	    (struct bv_context *)view_ctx);
}

extern "C" GED_EXPORT int
ged_view_context_update_callback_set(struct ged_view_context *view_ctx,
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

extern "C" GED_EXPORT ged_display_endpoint_t *
ged_view_context_display_endpoint_get(const struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);

    return record ?
	reinterpret_cast<ged_display_endpoint_t *>(record->display_endpoint) :
	NULL;
}

void
ged_bobol_view_controllers_foreach(struct ged *gedp,
    ged_bobol_view_controller_visit_t callback, void *userdata)
{
    struct ged_view_state_storage *state = ged_view_state(gedp);
    if (!state || !callback)
	return;
    for (size_t i = 0; i < BU_PTBL_LEN(&state->host_records); i++) {
	struct ged_view_host_record *record =
	    (struct ged_view_host_record *)BU_PTBL_GET(&state->host_records, i);
	if (!record || !record->view_ctx || !record->display_endpoint)
	    continue;
	BObolViewController *controller = static_cast<BObolViewController *>(
	    bobol_display_endpoint_controller(record->display_endpoint));
	if (controller && !callback(record->view_ctx, controller, userdata))
	    break;
    }
}

extern "C" GED_EXPORT int
ged_view_context_display_endpoint_ensure(struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);
    if (!record || !record->gedp)
	return 0;
    if (record->display_endpoint)
	return 1;
    return ged_draw_obol_render_endpoint_ensure_for_view(record->gedp,
	view_ctx, 1);
}

/* These adapters own GED view policy even before a display endpoint exists.
 * A null-DM camera context is used by headless clients such as rtwizard. */
static int ged_endpoint_property_get(void *user_data, const char *name,
	struct bv_display_property_value *out);
static int ged_endpoint_property_set(void *user_data, const char *name,
	const struct bv_display_property_value *value);

extern "C" GED_EXPORT int
ged_view_context_display_property_get(const struct ged_view_context *view_ctx,
	const char *name, struct bv_display_property_value *value)
{
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    if (!endpoint)
	return ged_endpoint_property_get((void *)view_ctx, name, value);
    return bobol_display_endpoint_property_get(endpoint, name, value);
}

extern "C" GED_EXPORT int
ged_view_context_display_property_set(struct ged_view_context *view_ctx,
	const char *name,
	const struct bv_display_property_value *value)
{
    bobol_display_endpoint_t *endpoint =
	ged_view_context_obol_endpoint_get(view_ctx);
    if (!endpoint)
	return ged_endpoint_property_set(view_ctx, name, value);
    return bobol_display_endpoint_property_set(endpoint, name, value);
}

static const char *
ged_endpoint_axes_suffix(const char *name, int *model_axes)
{
    static const char model_prefix[] = "view.faceplate.model_axes.";
    static const char view_prefix[] = "view.faceplate.view_axes.";
    if (!name || !model_axes)
	return NULL;
    if (bu_strncmp(name, model_prefix, sizeof(model_prefix) - 1) == 0) {
	*model_axes = 1;
	return name + sizeof(model_prefix) - 1;
    }
    if (bu_strncmp(name, view_prefix, sizeof(view_prefix) - 1) == 0) {
	*model_axes = 0;
	return name + sizeof(view_prefix) - 1;
    }
    return NULL;
}

static int
ged_endpoint_adc_property_get(const struct bv *view, const char *name,
	struct bv_display_property_value *out)
{
    static const char prefix[] = "view.faceplate.adc.";
    if (!name || bu_strncmp(name, prefix, sizeof(prefix) - 1) != 0)
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;

    struct bv_adc_state state = {};
    if (!view || !bv_adc_state_get(&state, view))
	return BV_DISPLAY_PROPERTY_INVALID;
    const char *suffix = name + sizeof(prefix) - 1;
    if (BU_STR_EQUAL(suffix, "visible")) {
	out->bool_value = state.draw ? 1 : 0;
    } else if (BU_STR_EQUAL(suffix, "line_width")) {
	out->uint_value = state.line_width;
    } else if (BU_STR_EQUAL(suffix, "line_color") ||
	BU_STR_EQUAL(suffix, "tick_color")) {
	const int *color = BU_STR_EQUAL(suffix, "line_color") ?
	    state.line_color : state.tick_color;
	for (int i = 0; i < 3; i++) {
	    if (color[i] < 0 || color[i] > 255)
		return BV_DISPLAY_PROPERTY_INVALID;
	    out->color3[i] = color[i] / 255.0;
	}
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_adc_property_set(struct ged_view_context *view_ctx, const char *name,
	const struct bv_display_property_value *value)
{
    static const char prefix[] = "view.faceplate.adc.";
    if (!name || bu_strncmp(name, prefix, sizeof(prefix) - 1) != 0)
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_adc_state state = {};
    if (!view || !value || !bv_adc_state_get(&state, view))
	return BV_DISPLAY_PROPERTY_INVALID;
    const char *suffix = name + sizeof(prefix) - 1;
    if (BU_STR_EQUAL(suffix, "visible")) {
	state.draw = value->bool_value ? 1 : 0;
    } else if (BU_STR_EQUAL(suffix, "line_width")) {
	state.line_width = static_cast<int>(value->uint_value);
    } else if (BU_STR_EQUAL(suffix, "line_color") ||
	BU_STR_EQUAL(suffix, "tick_color")) {
	int *color = BU_STR_EQUAL(suffix, "line_color") ?
	    state.line_color : state.tick_color;
	for (int i = 0; i < 3; i++)
	    color[i] = static_cast<int>(std::lround(value->color3[i] * 255.0));
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }
    if (!bv_adc_state_set(view, &state))
	return BV_DISPLAY_PROPERTY_INVALID;
    (void)ged_view_context_update(view_ctx);
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_axes_property_get(const struct bv *view, const char *name,
	struct bv_display_property_value *out)
{
    int model_axes = 0;
    const char *suffix = ged_endpoint_axes_suffix(name, &model_axes);
    if (!suffix)
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;

    struct bv_axes_state state = {};
    if (!view || !(model_axes ? bv_model_axes_state_get(&state, view) :
	    bv_view_axes_state_get(&state, view)))
	return BV_DISPLAY_PROPERTY_INVALID;

    if (BU_STR_EQUAL(suffix, "visible"))
	out->bool_value = state.draw ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "position.x"))
	out->double_value = state.axes_pos[X];
    else if (BU_STR_EQUAL(suffix, "position.y"))
	out->double_value = state.axes_pos[Y];
    else if (BU_STR_EQUAL(suffix, "position.z"))
	out->double_value = state.axes_pos[Z];
    else if (BU_STR_EQUAL(suffix, "size"))
	out->double_value = state.axes_size;
    else if (BU_STR_EQUAL(suffix, "line_width"))
	out->uint_value = state.line_width;
    else if (BU_STR_EQUAL(suffix, "color")) {
	for (int i = 0; i < 3; i++)
	    out->color3[i] = state.axes_color[i] / 255.0;
    } else if (BU_STR_EQUAL(suffix, "position_only"))
	out->bool_value = state.pos_only ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "labels.visible"))
	out->bool_value = state.label_flag ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "labels.color")) {
	for (int i = 0; i < 3; i++)
	    out->color3[i] = state.label_color[i] / 255.0;
    } else if (BU_STR_EQUAL(suffix, "triple_color"))
	out->bool_value = state.triple_color ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "ticks.visible"))
	out->bool_value = state.tick_enabled ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "ticks.length"))
	out->uint_value = state.tick_length;
    else if (BU_STR_EQUAL(suffix, "ticks.major_length"))
	out->uint_value = state.tick_major_length;
    else if (BU_STR_EQUAL(suffix, "ticks.interval"))
	out->double_value = state.tick_interval;
    else if (BU_STR_EQUAL(suffix, "ticks.per_major"))
	out->uint_value = state.ticks_per_major;
    else if (BU_STR_EQUAL(suffix, "ticks.threshold"))
	out->uint_value = state.tick_threshold;
    else if (BU_STR_EQUAL(suffix, "ticks.color")) {
	for (int i = 0; i < 3; i++)
	    out->color3[i] = state.tick_color[i] / 255.0;
    } else if (BU_STR_EQUAL(suffix, "ticks.major_color")) {
	for (int i = 0; i < 3; i++)
	    out->color3[i] = state.tick_major_color[i] / 255.0;
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_axes_property_set(struct ged_view_context *view_ctx, const char *name,
	const struct bv_display_property_value *value)
{
    int model_axes = 0;
    const char *suffix = ged_endpoint_axes_suffix(name, &model_axes);
    if (!suffix)
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_axes_state state = {};
    if (!view || !value || !(model_axes ?
	bv_model_axes_state_get(&state, view) :
	bv_view_axes_state_get(&state, view)))
	return BV_DISPLAY_PROPERTY_INVALID;

    if (BU_STR_EQUAL(suffix, "visible"))
	state.draw = value->bool_value ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "position.x"))
	state.axes_pos[X] = value->double_value;
    else if (BU_STR_EQUAL(suffix, "position.y"))
	state.axes_pos[Y] = value->double_value;
    else if (BU_STR_EQUAL(suffix, "position.z"))
	state.axes_pos[Z] = value->double_value;
    else if (BU_STR_EQUAL(suffix, "size"))
	state.axes_size = value->double_value;
    else if (BU_STR_EQUAL(suffix, "line_width"))
	state.line_width = static_cast<int>(value->uint_value);
    else if (BU_STR_EQUAL(suffix, "color")) {
	for (int i = 0; i < 3; i++)
	    state.axes_color[i] = static_cast<int>(std::lround(
		value->color3[i] * 255.0));
    } else if (BU_STR_EQUAL(suffix, "position_only"))
	state.pos_only = value->bool_value ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "labels.visible"))
	state.label_flag = value->bool_value ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "labels.color")) {
	for (int i = 0; i < 3; i++)
	    state.label_color[i] = static_cast<int>(std::lround(
		value->color3[i] * 255.0));
    } else if (BU_STR_EQUAL(suffix, "triple_color"))
	state.triple_color = value->bool_value ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "ticks.visible"))
	state.tick_enabled = value->bool_value ? 1 : 0;
    else if (BU_STR_EQUAL(suffix, "ticks.length"))
	state.tick_length = static_cast<int>(value->uint_value);
    else if (BU_STR_EQUAL(suffix, "ticks.major_length"))
	state.tick_major_length = static_cast<int>(value->uint_value);
    else if (BU_STR_EQUAL(suffix, "ticks.interval"))
	state.tick_interval = value->double_value;
    else if (BU_STR_EQUAL(suffix, "ticks.per_major"))
	state.ticks_per_major = static_cast<int>(value->uint_value);
    else if (BU_STR_EQUAL(suffix, "ticks.threshold"))
	state.tick_threshold = static_cast<int>(value->uint_value);
    else if (BU_STR_EQUAL(suffix, "ticks.color")) {
	for (int i = 0; i < 3; i++)
	    state.tick_color[i] = static_cast<int>(std::lround(
		value->color3[i] * 255.0));
    } else if (BU_STR_EQUAL(suffix, "ticks.major_color")) {
	for (int i = 0; i < 3; i++)
	    state.tick_major_color[i] = static_cast<int>(std::lround(
		value->color3[i] * 255.0));
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }

    if (!(model_axes ? bv_model_axes_state_set(view, &state) :
	bv_view_axes_state_set(view, &state)))
	return BV_DISPLAY_PROPERTY_INVALID;
    (void)ged_view_context_update(view_ctx);
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_property_get(void *user_data, const char *name,
	struct bv_display_property_value *out)
{
    struct ged_view_context *view_ctx =
	(struct ged_view_context *)user_data;
    if (!view_ctx || !name || !out)
	return BV_DISPLAY_PROPERTY_INVALID;
    const struct bv *view = bv_context_view_const(
	(const struct bv_context *)view_ctx);
    const int adc_result = ged_endpoint_adc_property_get(view, name, out);
    if (adc_result != BV_DISPLAY_PROPERTY_UNSUPPORTED)
	return adc_result;
    const int axes_result = ged_endpoint_axes_property_get(view, name, out);
    if (axes_result != BV_DISPLAY_PROPERTY_UNSUPPORTED)
	return axes_result;
    if (BU_STR_EQUAL(name, "view.perspective")) {
	out->double_value = bv_perspective_get(view);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.zclip")) {
	out->bool_value = bv_zclip_get(view) ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.navigation.min_delta") ||
	BU_STR_EQUAL(name, "view.navigation.max_delta") ||
	BU_STR_EQUAL(name, "view.navigation.rotate_scale") ||
	BU_STR_EQUAL(name, "view.navigation.scale_scale")) {
	struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;
	if (!bv_mouse_delta_settings_get(&settings, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.navigation.min_delta"))
	    out->double_value = settings.min_delta;
	else if (BU_STR_EQUAL(name, "view.navigation.max_delta"))
	    out->double_value = settings.max_delta;
	else if (BU_STR_EQUAL(name, "view.navigation.rotate_scale"))
	    out->double_value = settings.rotate_scale;
	else
	    out->double_value = settings.scale_scale;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.center_dot.visible")) {
	struct bv_other_state state = {};
	if (!view || !bv_center_dot_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	out->bool_value = state.gos_draw ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.grid.visible")) {
	struct bv_grid_state state = {};
	if (!view || !bv_grid_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	out->bool_value = state.draw ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.grid.adaptive") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.snap") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.anchor.x") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.anchor.y") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.anchor.z") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.resolution.horizontal") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.resolution.vertical") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.major.horizontal") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.major.vertical") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.color")) {
	struct bv_grid_state state = {};
	if (!view || !bv_grid_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.faceplate.grid.adaptive"))
	    out->bool_value = state.adaptive ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.snap"))
	    out->bool_value = state.snap ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.anchor.x"))
	    out->double_value = state.anchor[X];
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.anchor.y"))
	    out->double_value = state.anchor[Y];
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.anchor.z"))
	    out->double_value = state.anchor[Z];
	else if (BU_STR_EQUAL(name,
		"view.faceplate.grid.resolution.horizontal"))
	    out->double_value = state.res_h;
	else if (BU_STR_EQUAL(name,
		"view.faceplate.grid.resolution.vertical"))
	    out->double_value = state.res_v;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.major.horizontal"))
	    out->uint_value = state.res_major_h;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.major.vertical"))
	    out->uint_value = state.res_major_v;
	else {
	    for (int i = 0; i < 3; i++)
		out->color3[i] = state.color[i] / 255.0;
	}
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.interactive.rectangle.visible") ||
	BU_STR_EQUAL(name, "view.interactive.rectangle.line_width") ||
	BU_STR_EQUAL(name, "view.interactive.rectangle.line_style") ||
	BU_STR_EQUAL(name, "view.interactive.rectangle.color")) {
	struct bv_interactive_rect_state state = BV_INTERACTIVE_RECT_STATE_INIT;
	if (!view || !bv_interactive_rect_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.interactive.rectangle.visible"))
	    out->bool_value = state.draw ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.interactive.rectangle.line_width")) {
	    if (state.line_width < 0)
		return BV_DISPLAY_PROPERTY_INVALID;
	    out->uint_value = static_cast<uint64_t>(state.line_width);
	} else if (BU_STR_EQUAL(name, "view.interactive.rectangle.line_style")) {
	    if (state.line_style < 0 || state.line_style > 1)
		return BV_DISPLAY_PROPERTY_INVALID;
	    out->uint_value = static_cast<uint64_t>(state.line_style);
	} else {
	    for (int i = 0; i < 3; i++) {
		if (state.color[i] < 0 || state.color[i] > 255)
		    return BV_DISPLAY_PROPERTY_INVALID;
		out->color3[i] = state.color[i] / 255.0;
	    }
	}
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.scale.visible")) {
	struct bv_other_state state = {};
	if (!view || !bv_scale_overlay_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	out->bool_value = state.gos_draw ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.params.visible") ||
	BU_STR_EQUAL(name, "view.faceplate.params.size") ||
	BU_STR_EQUAL(name, "view.faceplate.params.center") ||
	BU_STR_EQUAL(name, "view.faceplate.params.azimuth") ||
	BU_STR_EQUAL(name, "view.faceplate.params.elevation") ||
	BU_STR_EQUAL(name, "view.faceplate.params.twist") ||
	BU_STR_EQUAL(name, "view.faceplate.params.fps")) {
	struct bv_params_state state = {};
	if (!view || !bv_params_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.faceplate.params.visible"))
	    out->bool_value = state.draw ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.size"))
	    out->bool_value = state.draw_size ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.center"))
	    out->bool_value = state.draw_center ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.azimuth"))
	    out->bool_value = state.draw_az ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.elevation"))
	    out->bool_value = state.draw_el ? 1 : 0;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.twist"))
	    out->bool_value = state.draw_tw ? 1 : 0;
	else
	    out->bool_value = state.draw_fps ? 1 : 0;
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.params.font_size")) {
	struct bv_params_state state = {};
	if (!view || !bv_params_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	out->uint_value = static_cast<uint64_t>(state.font_size);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.center_dot.font_size") ||
	BU_STR_EQUAL(name, "view.faceplate.scale.font_size")) {
	struct bv_other_state state = {};
	const int is_center_dot =
	    BU_STR_EQUAL(name, "view.faceplate.center_dot.font_size");
	if (!view || !(is_center_dot ? bv_center_dot_state_get(&state, view) :
		bv_scale_overlay_state_get(&state, view)))
	    return BV_DISPLAY_PROPERTY_INVALID;
	out->uint_value = static_cast<uint64_t>(state.gos_font_size);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.params.color")) {
	struct bv_params_state state = {};
	if (!view || !bv_params_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	for (int i = 0; i < 3; i++) {
	    if (state.color[i] < 0 || state.color[i] > 255)
		return BV_DISPLAY_PROPERTY_INVALID;
	    out->color3[i] = state.color[i] / 255.0;
	}
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.center_dot.color") ||
	BU_STR_EQUAL(name, "view.faceplate.scale.color")) {
	struct bv_other_state state = {};
	const int is_center_dot =
	    BU_STR_EQUAL(name, "view.faceplate.center_dot.color");
	if (!view || !(is_center_dot ? bv_center_dot_state_get(&state, view) :
		bv_scale_overlay_state_get(&state, view)))
	    return BV_DISPLAY_PROPERTY_INVALID;
	for (int i = 0; i < 3; i++) {
	    if (state.gos_line_color[i] < 0 || state.gos_line_color[i] > 255)
		return BV_DISPLAY_PROPERTY_INVALID;
	    out->color3[i] = state.gos_line_color[i] / 255.0;
	}
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "composition.framebuffer.mode")) {
	const int mode = bv_framebuffer_mode_get(
	    bv_context_view_const((const struct bv_context *)view_ctx));
	if (mode == BV_FRAMEBUFFER_MODE_OFF)
	    out->string_value = "off";
	else if (mode == BV_FRAMEBUFFER_MODE_OVERLAY)
	    out->string_value = "overlay";
	else if (mode == BV_FRAMEBUFFER_MODE_UNDERLAY)
	    out->string_value = "underlay";
	else if (mode == BV_FRAMEBUFFER_MODE_INTERLAY)
	    out->string_value = "interlay";
	else
	    return BV_DISPLAY_PROPERTY_INVALID;
	return BV_DISPLAY_PROPERTY_OK;
    }
    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
}

static int
ged_endpoint_perspective_set(struct ged_view_context *view_ctx, double perspective)
{
    if (!std::isfinite(perspective) || perspective < 0.0 ||
	perspective >= 180.0)
	return BV_DISPLAY_PROPERTY_INVALID;

    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    if (!view || !bv_perspective_set(view, perspective))
	return BV_DISPLAY_PROPERTY_INVALID;

    mat_t pmat;
    if (perspective > SMALL_FASTF) {
	persp_mat(pmat, perspective, (fastf_t)1.0, (fastf_t)0.01,
		  (fastf_t)1.0e10, (fastf_t)1.0);
    } else {
	MAT_IDN(pmat);
    }
    if (!bv_pmat_set(view, pmat))
	return BV_DISPLAY_PROPERTY_INVALID;

    (void)ged_view_context_update(view_ctx);
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_navigation_set(struct ged_view_context *view_ctx, const char *name, double value)
{
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    struct bv_mouse_delta_settings settings = BV_MOUSE_DELTA_SETTINGS_INIT;

    if (!view || !std::isfinite(value) ||
	!bv_mouse_delta_settings_get(&settings, view))
	return BV_DISPLAY_PROPERTY_INVALID;
    if (BU_STR_EQUAL(name, "view.navigation.min_delta"))
	settings.min_delta = value;
    else if (BU_STR_EQUAL(name, "view.navigation.max_delta"))
	settings.max_delta = value;
    else if (BU_STR_EQUAL(name, "view.navigation.rotate_scale"))
	settings.rotate_scale = value;
    else if (BU_STR_EQUAL(name, "view.navigation.scale_scale"))
	settings.scale_scale = value;
    else
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;

    if (!bv_mouse_delta_settings_set(view, &settings))
	return BV_DISPLAY_PROPERTY_INVALID;
    (void)ged_view_context_update(view_ctx);
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_faceplate_font_size_set(struct ged_view_context *view_ctx, const char *name,
	uint64_t font_size)
{
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    if (!view || font_size < 5 || font_size > 96)
	return BV_DISPLAY_PROPERTY_INVALID;

    if (BU_STR_EQUAL(name, "view.faceplate.params.font_size")) {
	struct bv_params_state state = {};
	if (!bv_params_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	state.font_size = static_cast<int>(font_size);
	if (!bv_params_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (BU_STR_EQUAL(name, "view.faceplate.center_dot.font_size") ||
	BU_STR_EQUAL(name, "view.faceplate.scale.font_size")) {
	struct bv_other_state state = {};
	const int is_center_dot =
	    BU_STR_EQUAL(name, "view.faceplate.center_dot.font_size");
	if (!(is_center_dot ? bv_center_dot_state_get(&state, view) :
		bv_scale_overlay_state_get(&state, view)))
	    return BV_DISPLAY_PROPERTY_INVALID;
	state.gos_font_size = static_cast<int>(font_size);
	if (!(is_center_dot ? bv_center_dot_state_set(view, &state) :
		bv_scale_overlay_state_set(view, &state)))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }

    (void)ged_view_context_update(view_ctx);
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_faceplate_color_set(struct ged_view_context *view_ctx, const char *name,
	const double color[3])
{
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);

    if (!view || !color)
	return BV_DISPLAY_PROPERTY_INVALID;
    for (int i = 0; i < 3; i++) {
	if (!std::isfinite(color[i]) || color[i] < 0.0 || color[i] > 1.0)
	    return BV_DISPLAY_PROPERTY_INVALID;
    }

    if (BU_STR_EQUAL(name, "view.faceplate.params.color")) {
	struct bv_params_state state = {};
	if (!bv_params_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	for (int i = 0; i < 3; i++)
	    state.color[i] = static_cast<int>(std::lround(color[i] * 255.0));
	if (!bv_params_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else if (BU_STR_EQUAL(name, "view.faceplate.center_dot.color") ||
	BU_STR_EQUAL(name, "view.faceplate.scale.color")) {
	struct bv_other_state state = {};
	const int is_center_dot =
	    BU_STR_EQUAL(name, "view.faceplate.center_dot.color");
	if (!(is_center_dot ? bv_center_dot_state_get(&state, view) :
		bv_scale_overlay_state_get(&state, view)))
	    return BV_DISPLAY_PROPERTY_INVALID;
	for (int i = 0; i < 3; i++)
	    state.gos_line_color[i] =
		static_cast<int>(std::lround(color[i] * 255.0));
	if (!(is_center_dot ? bv_center_dot_state_set(view, &state) :
		bv_scale_overlay_state_set(view, &state)))
	    return BV_DISPLAY_PROPERTY_INVALID;
    } else {
	return BV_DISPLAY_PROPERTY_UNSUPPORTED;
    }

    (void)ged_view_context_update(view_ctx);
    return BV_DISPLAY_PROPERTY_OK;
}

static int
ged_endpoint_property_set(void *user_data, const char *name,
	const struct bv_display_property_value *value)
{
    struct ged_view_context *view_ctx =
	(struct ged_view_context *)user_data;
    if (!view_ctx || !name || !value)
	return BV_DISPLAY_PROPERTY_INVALID;
    const int adc_result = ged_endpoint_adc_property_set(view_ctx, name,
	value);
    if (adc_result != BV_DISPLAY_PROPERTY_UNSUPPORTED)
	return adc_result;
    const int axes_result = ged_endpoint_axes_property_set(view_ctx, name,
	value);
    if (axes_result != BV_DISPLAY_PROPERTY_UNSUPPORTED)
	return axes_result;
    if (BU_STR_EQUAL(name, "view.perspective"))
	return ged_endpoint_perspective_set(view_ctx, value->double_value);
    if (BU_STR_EQUAL(name, "view.zclip")) {
	if (!bv_zclip_set(bv_context_view((struct bv_context *)view_ctx),
		value->bool_value))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.navigation.min_delta") ||
	BU_STR_EQUAL(name, "view.navigation.max_delta") ||
	BU_STR_EQUAL(name, "view.navigation.rotate_scale") ||
	BU_STR_EQUAL(name, "view.navigation.scale_scale"))
	return ged_endpoint_navigation_set(view_ctx, name, value->double_value);
    if (BU_STR_EQUAL(name, "view.faceplate.params.font_size") ||
	BU_STR_EQUAL(name, "view.faceplate.center_dot.font_size") ||
	BU_STR_EQUAL(name, "view.faceplate.scale.font_size"))
	return ged_endpoint_faceplate_font_size_set(view_ctx, name,
		value->uint_value);
    if (BU_STR_EQUAL(name, "view.faceplate.params.color") ||
	BU_STR_EQUAL(name, "view.faceplate.center_dot.color") ||
	BU_STR_EQUAL(name, "view.faceplate.scale.color"))
	return ged_endpoint_faceplate_color_set(view_ctx, name, value->color3);
    struct bv *view = bv_context_view((struct bv_context *)view_ctx);
    const int enabled = value->bool_value ? 1 : 0;
    if (BU_STR_EQUAL(name, "view.faceplate.center_dot.visible")) {
	struct bv_other_state state = {};
	if (!view || !bv_center_dot_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	state.gos_draw = enabled;
	if (!bv_center_dot_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.grid.visible")) {
	struct bv_grid_state state = {};
	if (!view || !bv_grid_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	state.draw = enabled;
	if (!bv_grid_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.grid.adaptive") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.snap") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.anchor.x") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.anchor.y") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.anchor.z") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.resolution.horizontal") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.resolution.vertical") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.major.horizontal") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.major.vertical") ||
	BU_STR_EQUAL(name, "view.faceplate.grid.color")) {
	struct bv_grid_state state = {};
	if (!view || !bv_grid_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.faceplate.grid.adaptive"))
	    state.adaptive = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.snap"))
	    state.snap = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.anchor.x"))
	    state.anchor[X] = value->double_value;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.anchor.y"))
	    state.anchor[Y] = value->double_value;
	else if (BU_STR_EQUAL(name, "view.faceplate.grid.anchor.z"))
	    state.anchor[Z] = value->double_value;
	else if (BU_STR_EQUAL(name,
		"view.faceplate.grid.resolution.horizontal")) {
	    if (value->double_value <= 0.0)
		return BV_DISPLAY_PROPERTY_INVALID;
	    state.res_h = value->double_value;
	} else if (BU_STR_EQUAL(name,
		"view.faceplate.grid.resolution.vertical")) {
	    if (value->double_value <= 0.0)
		return BV_DISPLAY_PROPERTY_INVALID;
	    state.res_v = value->double_value;
	} else if (BU_STR_EQUAL(name, "view.faceplate.grid.major.horizontal")) {
	    state.res_major_h = static_cast<int>(value->uint_value);
	} else if (BU_STR_EQUAL(name, "view.faceplate.grid.major.vertical")) {
	    state.res_major_v = static_cast<int>(value->uint_value);
	} else {
	    for (int i = 0; i < 3; i++)
		state.color[i] = static_cast<int>(std::lround(
		    value->color3[i] * 255.0));
	}
	if (!bv_grid_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.interactive.rectangle.visible") ||
	BU_STR_EQUAL(name, "view.interactive.rectangle.line_width") ||
	BU_STR_EQUAL(name, "view.interactive.rectangle.line_style") ||
	BU_STR_EQUAL(name, "view.interactive.rectangle.color")) {
	struct bv_interactive_rect_state state = BV_INTERACTIVE_RECT_STATE_INIT;
	if (!view || !bv_interactive_rect_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.interactive.rectangle.visible")) {
	    state.draw = enabled;
	} else if (BU_STR_EQUAL(name,
		"view.interactive.rectangle.line_width")) {
	    state.line_width = static_cast<int>(value->uint_value);
	} else if (BU_STR_EQUAL(name,
		"view.interactive.rectangle.line_style")) {
	    if (value->uint_value > 1)
		return BV_DISPLAY_PROPERTY_INVALID;
	    state.line_style = static_cast<int>(value->uint_value);
	} else {
	    for (int i = 0; i < 3; i++)
		state.color[i] = static_cast<int>(std::lround(
		    value->color3[i] * 255.0));
	}
	if (!bv_interactive_rect_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.scale.visible")) {
	struct bv_other_state state = {};
	if (!view || !bv_scale_overlay_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	state.gos_draw = enabled;
	if (!bv_scale_overlay_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "view.faceplate.params.visible") ||
	BU_STR_EQUAL(name, "view.faceplate.params.size") ||
	BU_STR_EQUAL(name, "view.faceplate.params.center") ||
	BU_STR_EQUAL(name, "view.faceplate.params.azimuth") ||
	BU_STR_EQUAL(name, "view.faceplate.params.elevation") ||
	BU_STR_EQUAL(name, "view.faceplate.params.twist") ||
	BU_STR_EQUAL(name, "view.faceplate.params.fps")) {
	struct bv_params_state state = {};
	if (!view || !bv_params_state_get(&state, view))
	    return BV_DISPLAY_PROPERTY_INVALID;
	if (BU_STR_EQUAL(name, "view.faceplate.params.visible"))
	    state.draw = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.size"))
	    state.draw_size = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.center"))
	    state.draw_center = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.azimuth"))
	    state.draw_az = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.elevation"))
	    state.draw_el = enabled;
	else if (BU_STR_EQUAL(name, "view.faceplate.params.twist"))
	    state.draw_tw = enabled;
	else
	    state.draw_fps = enabled;
	if (!bv_params_state_set(view, &state))
	    return BV_DISPLAY_PROPERTY_INVALID;
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    if (BU_STR_EQUAL(name, "composition.framebuffer.mode")) {
	if (!value->string_value)
	    return BV_DISPLAY_PROPERTY_INVALID;
	int mode = BV_FRAMEBUFFER_MODE_OFF;
	if (BU_STR_EQUAL(value->string_value, "overlay"))
	    mode = BV_FRAMEBUFFER_MODE_OVERLAY;
	else if (BU_STR_EQUAL(value->string_value, "underlay"))
	    mode = BV_FRAMEBUFFER_MODE_UNDERLAY;
	else if (BU_STR_EQUAL(value->string_value, "interlay"))
	    mode = BV_FRAMEBUFFER_MODE_INTERLAY;
	else if (!BU_STR_EQUAL(value->string_value, "off"))
	    return BV_DISPLAY_PROPERTY_INVALID;
	struct bv *mutable_view = bv_context_view(
	    (struct bv_context *)view_ctx);
	const int previous_mode = bv_framebuffer_mode_get(mutable_view);
	if (!bv_framebuffer_mode_set(mutable_view, mode))
	    return BV_DISPLAY_PROPERTY_INVALID;
	struct ged_view_host_record *record =
	    ged_view_host_record_find_global(view_ctx);
	if (record && record->gedp &&
	    ged_obol_fbserv_composition_set(record->gedp, mode) != BRLCAD_OK) {
	    (void)bv_framebuffer_mode_set(mutable_view, previous_mode);
	    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
	}
	(void)ged_view_context_update(view_ctx);
	return BV_DISPLAY_PROPERTY_OK;
    }
    return BV_DISPLAY_PROPERTY_UNSUPPORTED;
}

static void
ged_endpoint_background_initialize(struct ged_view_context *view_ctx,
	bobol_display_endpoint_t *endpoint)
{
    const struct bv *view = bv_context_view_const(
	(const struct bv_context *)view_ctx);
    struct bv_background_state background = BV_BACKGROUND_STATE_INIT;
    if (!endpoint || !view || !bv_background_state_get(&background, view))
	return;

    struct bv_display_property_value value =
	BV_DISPLAY_PROPERTY_VALUE_INIT;
    value.type = BV_DISPLAY_PROPERTY_COLOR3;
    for (int i = 0; i < 3; i++)
	value.color3[i] = background.bottom[i] / 255.0;
    if (bobol_display_endpoint_property_set(endpoint,
	    "controller.background.bottom", &value) !=
	BV_DISPLAY_PROPERTY_OK)
	return;
    for (int i = 0; i < 3; i++)
	value.color3[i] = background.top[i] / 255.0;
    (void)bobol_display_endpoint_property_set(endpoint,
	"controller.background.top", &value);
}

extern "C" GED_EXPORT int
ged_view_context_display_endpoint_set(
	struct ged_view_context *view_ctx,
	ged_display_endpoint_t *display_endpoint,
	int take_ownership)
{
    bobol_display_endpoint_t *endpoint =
	reinterpret_cast<bobol_display_endpoint_t *>(display_endpoint);
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);
    if (!record)
	return 0;

    if (record->display_endpoint == endpoint) {
	record->owns_display_endpoint = take_ownership ? 1 : 0;
	if (endpoint) {
	    (void)bobol_display_endpoint_property_provider_set(endpoint,
		ged_endpoint_property_get, ged_endpoint_property_set, view_ctx);
	    /* bv state remains authoritative even while a view is headless.
	     * Reassert presentation policy when an endpoint is registered so a
	     * controller replacement or delayed GUI attachment cannot silently
	     * restore its own defaults. */
	    (void)ged_view_shading_sync(record->gedp, view_ctx);
	}
	return 1;
    }

    bobol_display_endpoint_t *old_endpoint = record->display_endpoint;
    const int owned_old_endpoint = record->owns_display_endpoint;

    if (endpoint) {
	void *controller = bobol_display_endpoint_controller(endpoint);
	if (!controller || !record->gedp)
	    return 0;
	void *old_controller = old_endpoint ?
	    bobol_display_endpoint_controller(old_endpoint) : NULL;
	const int detached_old = old_controller && old_controller != controller;
	if (detached_old)
	    ged_draw_obol_controller_detach_for_view(record->gedp, view_ctx);
	if (!ged_draw_obol_controller_attach_opaque_for_view(record->gedp,
		view_ctx, controller, 1)) {
	    if (detached_old)
		(void)ged_draw_obol_controller_attach_opaque_for_view(
		    record->gedp, view_ctx, old_controller, 1);
	    return 0;
	}
	} else if (record->gedp) {
	if (old_endpoint)
	    ged_draw_obol_framebuffer_endpoint_detach(record->gedp,
		old_endpoint);
	ged_draw_obol_controller_detach_for_view(record->gedp, view_ctx);
    }

    if (endpoint && old_endpoint)
	ged_draw_obol_framebuffer_endpoint_detach(record->gedp, old_endpoint);
    if (old_endpoint)
	(void)bobol_display_endpoint_property_provider_set(old_endpoint,
	    NULL, NULL, NULL);
    record->display_endpoint = endpoint;
    record->owns_display_endpoint = take_ownership ? 1 : 0;
    if (endpoint) {
	(void)bobol_display_endpoint_property_provider_set(endpoint,
	    ged_endpoint_property_get, ged_endpoint_property_set, view_ctx);
	ged_endpoint_background_initialize(view_ctx, endpoint);
	/* A headless client may configure normal presentation before its
	 * graphical endpoint exists.  Push that retained, per-view policy only
	 * after the new endpoint has become the live controller. */
	(void)ged_view_shading_sync(record->gedp, view_ctx);
    }

    if (old_endpoint && owned_old_endpoint)
	bobol_display_endpoint_destroy(old_endpoint);
    return 1;
}

extern "C" GED_EXPORT int
ged_view_context_update(struct ged_view_context *view_ctx)
{
    struct ged_view_host_record *record =
	ged_view_host_record_find_global(view_ctx);
    int ret = bv_context_update((struct bv_context *)view_ctx,
	    BV_CONTEXT_CHANGED_VIEW);

    if (ret && record && record->update_callback)
	record->update_callback(view_ctx, record->update_callback_data);

    return ret;
}
