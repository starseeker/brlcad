/*             V I E W _ C O N T E X T _ A D A P T E R S . C
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
/** @file librt/view_context_adapters.c
 *
 * Context-bound view adapters that let the neutral RT view API dispatch
 * higher-level scene, feature, polygon, and selection operations to whichever
 * retained scene implementation owns the supplied context or object ref.
 */

#include "common.h"

#include <string.h>

#include "bu/malloc.h"
#include "bu/ptbl.h"
#include "rt/view.h"
#include "view_context_private.h"

struct rt_view_scene_adapter_binding {
    void *ctx;
    struct rt_view_context_scene_adapter adapter;
};

struct rt_view_feature_adapter_binding {
    void *ctx;
    struct rt_view_context_feature_adapter adapter;
};

struct rt_view_polygon_adapter_binding {
    void *ctx;
    struct rt_view_context_polygon_adapter adapter;
};

struct rt_view_selection_adapter_binding {
    void *ctx;
    struct rt_view_context_selection_adapter adapter;
};

static struct bu_ptbl rt_view_scene_adapter_bindings = BU_PTBL_INIT_ZERO;
static int rt_view_scene_adapter_bindings_init = 0;
static struct bu_ptbl rt_view_feature_adapter_bindings = BU_PTBL_INIT_ZERO;
static int rt_view_feature_adapter_bindings_init = 0;
static struct bu_ptbl rt_view_polygon_adapter_bindings = BU_PTBL_INIT_ZERO;
static int rt_view_polygon_adapter_bindings_init = 0;
static struct bu_ptbl rt_view_selection_adapter_bindings = BU_PTBL_INIT_ZERO;
static int rt_view_selection_adapter_bindings_init = 0;

static void
rt_view_scene_adapter_bindings_ensure(void)
{
    if (!rt_view_scene_adapter_bindings_init) {
	BU_PTBL_INIT(&rt_view_scene_adapter_bindings);
	rt_view_scene_adapter_bindings_init = 1;
    }
}

static struct rt_view_scene_adapter_binding *
rt_view_scene_adapter_binding_find(const void *ctx)
{
    if (!ctx || !rt_view_scene_adapter_bindings_init)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&rt_view_scene_adapter_bindings); i++) {
	struct rt_view_scene_adapter_binding *binding =
	    (struct rt_view_scene_adapter_binding *)BU_PTBL_GET(
		&rt_view_scene_adapter_bindings, i);
	if (binding && binding->ctx == ctx)
	    return binding;
    }

    return NULL;
}

void
_rt_view_context_scene_adapter_clear(void *ctx)
{
    struct rt_view_scene_adapter_binding *binding =
	rt_view_scene_adapter_binding_find(ctx);
    if (!binding)
	return;

    (void)bu_ptbl_rm(&rt_view_scene_adapter_bindings, (long *)binding);
    BU_PUT(binding, struct rt_view_scene_adapter_binding);
}

int
rt_view_context_scene_adapter_set(
    void *ctx,
    const struct rt_view_context_scene_adapter *adapter)
{
    if (!ctx)
	return 0;

    if (!adapter || (!adapter->pick_semantic_path &&
		     !adapter->render_export_consistency)) {
	_rt_view_context_scene_adapter_clear(ctx);
	return 1;
    }

    rt_view_scene_adapter_bindings_ensure();
    struct rt_view_scene_adapter_binding *binding =
	rt_view_scene_adapter_binding_find(ctx);
    if (!binding) {
	BU_GET(binding, struct rt_view_scene_adapter_binding);
	memset(binding, 0, sizeof(*binding));
	binding->ctx = ctx;
	bu_ptbl_ins(&rt_view_scene_adapter_bindings, (long *)binding);
    }
    binding->adapter = *adapter;
    return 1;
}

int
rt_view_context_scene_adapter_get(
    void *ctx,
    struct rt_view_context_scene_adapter *adapter)
{
    if (adapter)
	memset(adapter, 0, sizeof(*adapter));
    if (!ctx || !adapter)
	return 0;

    struct rt_view_scene_adapter_binding *binding =
	rt_view_scene_adapter_binding_find(ctx);
    if (!binding)
	return 0;

    *adapter = binding->adapter;
    return 1;
}

static void
rt_view_feature_adapter_bindings_ensure(void)
{
    if (!rt_view_feature_adapter_bindings_init) {
	BU_PTBL_INIT(&rt_view_feature_adapter_bindings);
	rt_view_feature_adapter_bindings_init = 1;
    }
}

static struct rt_view_feature_adapter_binding *
rt_view_feature_adapter_binding_find(const void *ctx)
{
    if (!ctx || !rt_view_feature_adapter_bindings_init)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&rt_view_feature_adapter_bindings);
	 i++) {
	struct rt_view_feature_adapter_binding *binding =
	    (struct rt_view_feature_adapter_binding *)BU_PTBL_GET(
		&rt_view_feature_adapter_bindings, i);
	if (binding && binding->ctx == ctx)
	    return binding;
    }

    return NULL;
}

static int
rt_view_feature_adapter_is_empty(
    const struct rt_view_context_feature_adapter *adapter)
{
    return !adapter || (!adapter->owns_ref &&
			!adapter->edit_preview_publish_event &&
			!adapter->overlay_ensure && !adapter->label_ensure &&
			!adapter->remove && !adapter->set_context &&
			!adapter->set_visible && !adapter->set_color &&
			!adapter->touch && !adapter->labels_replace &&
			!adapter->points_replace && !adapter->clear_geometry);
}

void
_rt_view_context_feature_adapter_clear(void *ctx)
{
    struct rt_view_feature_adapter_binding *binding =
	rt_view_feature_adapter_binding_find(ctx);
    if (!binding)
	return;

    (void)bu_ptbl_rm(&rt_view_feature_adapter_bindings, (long *)binding);
    BU_PUT(binding, struct rt_view_feature_adapter_binding);
}

int
rt_view_context_feature_adapter_set(
    void *ctx,
    const struct rt_view_context_feature_adapter *adapter)
{
    if (!ctx)
	return 0;

    if (rt_view_feature_adapter_is_empty(adapter)) {
	_rt_view_context_feature_adapter_clear(ctx);
	return 1;
    }

    rt_view_feature_adapter_bindings_ensure();
    struct rt_view_feature_adapter_binding *binding =
	rt_view_feature_adapter_binding_find(ctx);
    if (!binding) {
	BU_GET(binding, struct rt_view_feature_adapter_binding);
	memset(binding, 0, sizeof(*binding));
	binding->ctx = ctx;
	bu_ptbl_ins(&rt_view_feature_adapter_bindings, (long *)binding);
    }
    binding->adapter = *adapter;
    return 1;
}

int
rt_view_context_feature_adapter_get(
    void *ctx,
    struct rt_view_context_feature_adapter *adapter)
{
    if (adapter)
	memset(adapter, 0, sizeof(*adapter));
    if (!ctx || !adapter)
	return 0;

    struct rt_view_feature_adapter_binding *binding =
	rt_view_feature_adapter_binding_find(ctx);
    if (!binding)
	return 0;

    *adapter = binding->adapter;
    return 1;
}

int
_rt_view_feature_adapter_for_ref(
    rt_view_feature_ref ref,
    struct rt_view_context_feature_adapter *adapter)
{
    if (adapter)
	memset(adapter, 0, sizeof(*adapter));
    if (!ref.token || !adapter || !rt_view_feature_adapter_bindings_init)
	return 0;

    for (size_t i = 0; i < BU_PTBL_LEN(&rt_view_feature_adapter_bindings);
	 i++) {
	struct rt_view_feature_adapter_binding *binding =
	    (struct rt_view_feature_adapter_binding *)BU_PTBL_GET(
		&rt_view_feature_adapter_bindings, i);
	if (!binding || !binding->adapter.owns_ref)
	    continue;
	if (binding->adapter.owns_ref(ref, binding->adapter.data)) {
	    *adapter = binding->adapter;
	    return 1;
	}
    }

    return 0;
}

static void
rt_view_polygon_adapter_bindings_ensure(void)
{
    if (!rt_view_polygon_adapter_bindings_init) {
	BU_PTBL_INIT(&rt_view_polygon_adapter_bindings);
	rt_view_polygon_adapter_bindings_init = 1;
    }
}

static struct rt_view_polygon_adapter_binding *
rt_view_polygon_adapter_binding_find(const void *ctx)
{
    if (!ctx || !rt_view_polygon_adapter_bindings_init)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&rt_view_polygon_adapter_bindings); i++) {
	struct rt_view_polygon_adapter_binding *binding =
	    (struct rt_view_polygon_adapter_binding *)BU_PTBL_GET(
		&rt_view_polygon_adapter_bindings, i);
	if (binding && binding->ctx == ctx)
	    return binding;
    }

    return NULL;
}

static int
rt_view_polygon_adapter_is_empty(
    const struct rt_view_context_polygon_adapter *adapter)
{
    return !adapter || (!adapter->owns_ref && !adapter->record_get &&
			!adapter->create && !adapter->select && !adapter->find &&
			!adapter->dup && !adapter->visit_records &&
			!adapter->snap_count && !adapter->clear_point_selection &&
			!adapter->update && !adapter->update_screen_pt &&
			!adapter->move && !adapter->set_name && !adapter->set_context &&
			!adapter->set_visual && !adapter->set_open && !adapter->close &&
			!adapter->clear_selected_point && !adapter->remove &&
			!adapter->user_data && !adapter->user_data_set &&
			!adapter->csg && !adapter->import_sketch_context &&
			!adapter->export_sketch && !adapter->snap_exclude_set);
}

void
_rt_view_context_polygon_adapter_clear(void *ctx)
{
    struct rt_view_polygon_adapter_binding *binding =
	rt_view_polygon_adapter_binding_find(ctx);
    if (!binding)
	return;

    (void)bu_ptbl_rm(&rt_view_polygon_adapter_bindings, (long *)binding);
    BU_PUT(binding, struct rt_view_polygon_adapter_binding);
}

int
rt_view_context_polygon_adapter_set(
    void *ctx,
    const struct rt_view_context_polygon_adapter *adapter)
{
    if (!ctx)
	return 0;

    if (rt_view_polygon_adapter_is_empty(adapter)) {
	_rt_view_context_polygon_adapter_clear(ctx);
	return 1;
    }

    rt_view_polygon_adapter_bindings_ensure();
    struct rt_view_polygon_adapter_binding *binding =
	rt_view_polygon_adapter_binding_find(ctx);
    if (!binding) {
	BU_GET(binding, struct rt_view_polygon_adapter_binding);
	memset(binding, 0, sizeof(*binding));
	binding->ctx = ctx;
	bu_ptbl_ins(&rt_view_polygon_adapter_bindings, (long *)binding);
    }
    binding->adapter = *adapter;
    return 1;
}

int
rt_view_context_polygon_adapter_get(
    void *ctx,
    struct rt_view_context_polygon_adapter *adapter)
{
    if (adapter)
	memset(adapter, 0, sizeof(*adapter));
    if (!ctx || !adapter)
	return 0;

    struct rt_view_polygon_adapter_binding *binding =
	rt_view_polygon_adapter_binding_find(ctx);
    if (!binding)
	return 0;

    *adapter = binding->adapter;
    return 1;
}

int
_rt_view_polygon_adapter_for_ref(
    rt_view_polygon_ref ref,
    struct rt_view_context_polygon_adapter *adapter)
{
    if (adapter)
	memset(adapter, 0, sizeof(*adapter));
    if (!ref.token || !adapter || !rt_view_polygon_adapter_bindings_init)
	return 0;

    for (size_t i = 0; i < BU_PTBL_LEN(&rt_view_polygon_adapter_bindings); i++) {
	struct rt_view_polygon_adapter_binding *binding =
	    (struct rt_view_polygon_adapter_binding *)BU_PTBL_GET(
		&rt_view_polygon_adapter_bindings, i);
	if (!binding || !binding->adapter.owns_ref)
	    continue;
	if (binding->adapter.owns_ref(ref, binding->adapter.data)) {
	    *adapter = binding->adapter;
	    return 1;
	}
    }

    return 0;
}

static void
rt_view_selection_adapter_bindings_ensure(void)
{
    if (!rt_view_selection_adapter_bindings_init) {
	BU_PTBL_INIT(&rt_view_selection_adapter_bindings);
	rt_view_selection_adapter_bindings_init = 1;
    }
}

static struct rt_view_selection_adapter_binding *
rt_view_selection_adapter_binding_find(const void *ctx)
{
    if (!ctx || !rt_view_selection_adapter_bindings_init)
	return NULL;

    for (size_t i = 0; i < BU_PTBL_LEN(&rt_view_selection_adapter_bindings);
	 i++) {
	struct rt_view_selection_adapter_binding *binding =
	    (struct rt_view_selection_adapter_binding *)BU_PTBL_GET(
		&rt_view_selection_adapter_bindings, i);
	if (binding && binding->ctx == ctx)
	    return binding;
    }

    return NULL;
}

static int
rt_view_selection_adapter_is_empty(
    const struct rt_view_context_selection_adapter *adapter)
{
    return !adapter || (!adapter->available && !adapter->count &&
			!adapter->set_pick_result_context && !adapter->clear);
}

void
_rt_view_context_selection_adapter_clear(void *ctx)
{
    struct rt_view_selection_adapter_binding *binding =
	rt_view_selection_adapter_binding_find(ctx);
    if (!binding)
	return;

    (void)bu_ptbl_rm(&rt_view_selection_adapter_bindings, (long *)binding);
    BU_PUT(binding, struct rt_view_selection_adapter_binding);
}

int
rt_view_context_selection_adapter_set(
    void *ctx,
    const struct rt_view_context_selection_adapter *adapter)
{
    if (!ctx)
	return 0;

    if (rt_view_selection_adapter_is_empty(adapter)) {
	_rt_view_context_selection_adapter_clear(ctx);
	return 1;
    }

    rt_view_selection_adapter_bindings_ensure();
    struct rt_view_selection_adapter_binding *binding =
	rt_view_selection_adapter_binding_find(ctx);
    if (!binding) {
	BU_GET(binding, struct rt_view_selection_adapter_binding);
	memset(binding, 0, sizeof(*binding));
	binding->ctx = ctx;
	bu_ptbl_ins(&rt_view_selection_adapter_bindings, (long *)binding);
    }
    binding->adapter = *adapter;
    return 1;
}

int
rt_view_context_selection_adapter_get(
    void *ctx,
    struct rt_view_context_selection_adapter *adapter)
{
    if (adapter)
	memset(adapter, 0, sizeof(*adapter));
    if (!ctx || !adapter)
	return 0;

    struct rt_view_selection_adapter_binding *binding =
	rt_view_selection_adapter_binding_find(ctx);
    if (!binding)
	return 0;

    *adapter = binding->adapter;
    return 1;
}
