/*                V I E W _ L E G A C Y _ B S G . C
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
/** @file librt/tests/view_legacy_bsg.c
 *
 * Unit tests for transitional BSG-to-RT view snapshot adapters.
 */

#include "common.h"

#include <math.h>
#include <stdio.h>
#include <string.h>

#include "bg/plane.h"
#include "bn/mat.h"
#include "bn/qmath.h"
#include "bu/app.h"
#include "bu/malloc.h"
#include "bu/str.h"
#include "bsg/feature.h"
#include "bsg/geometry.h"
#include "bsg/measure.h"
#include "bsg/polygon.h"
#include "bsg/scene_builder.h"
#include "bsg/snap.h"
#include "bsg/snap_action.h"
#include "bsg/util.h"
#include "bsg/view_set.h"
#include "bsg/view_state.h"
#include "rt/edit_legacy_bsg.h"
#include "rt/view.h"
#include "rt/view_legacy_bsg.h"

static int
fastf_equal(fastf_t a, fastf_t b)
{
    return fabs(a - b) <= SMALL_FASTF;
}

static int
check_quat(const char *label, const quat_t actual, const quat_t expected)
{
    if (!fastf_equal(actual[0], expected[0]) ||
	    !fastf_equal(actual[1], expected[1]) ||
	    !fastf_equal(actual[2], expected[2]) ||
	    !fastf_equal(actual[3], expected[3])) {
	printf("FAIL: %s quat: got %g,%g,%g,%g expected %g,%g,%g,%g\n",
		label,
		actual[0], actual[1], actual[2], actual[3],
		expected[0], expected[1], expected[2], expected[3]);
	return 1;
    }

    return 0;
}

static int
check_mat(const char *label, const mat_t actual, const mat_t expected)
{
    for (size_t i = 0; i < 16; i++) {
	if (!fastf_equal(actual[i], expected[i])) {
	    printf("FAIL: %s mat[%zu]: got %g expected %g\n",
		    label, i, actual[i], expected[i]);
	    return 1;
	}
    }

    return 0;
}

static int
check_view_frame(const char *label,
		 const struct bsg_view *actual,
		 const struct bsg_view *expected)
{
    if (check_mat(label, actual->gv_model2view, expected->gv_model2view) ||
	    check_mat(label, actual->gv_view2model, expected->gv_view2model) ||
	    check_mat(label, actual->gv_center, expected->gv_center) ||
	    check_mat(label, actual->gv_pmodel2view, expected->gv_pmodel2view)) {
	return 1;
    }

    if (!VNEAR_EQUAL(actual->gv_aet, expected->gv_aet, SMALL_FASTF) ||
	    !fastf_equal(actual->gv_scale, expected->gv_scale) ||
	    !fastf_equal(actual->gv_size, expected->gv_size) ||
	    !fastf_equal(actual->gv_isize, expected->gv_isize)) {
	printf("FAIL: %s frame mismatch\n", label);
	return 1;
    }

    return 0;
}

static int
check_view_info(const char *label,
		const struct rt_view_info *info,
		int width,
		int height,
		fastf_t size,
		fastf_t scale,
		fastf_t curve_scale,
		fastf_t point_scale,
		size_t bot_threshold)
{
    if (!info) {
	printf("FAIL: %s info is NULL\n", label);
	return 1;
    }

    if (info->width != width || info->height != height ||
	    !fastf_equal(info->size, size)) {
	printf("FAIL: %s dimensions/size: got %dx%d size %g\n",
		label, info->width, info->height, info->size);
	return 1;
    }

    if (!fastf_equal(info->lod.scale, scale) ||
	    !fastf_equal(info->lod.curve_scale, curve_scale) ||
	    !fastf_equal(info->lod.point_scale, point_scale) ||
	    info->lod.bot_threshold != bot_threshold) {
	printf("FAIL: %s lod: got scale=%g curve=%g point=%g threshold=%zu\n",
		label, info->lod.scale, info->lod.curve_scale,
		info->lod.point_scale, info->lod.bot_threshold);
	return 1;
    }

    return 0;
}

static int
check_lod_policy(const char *label,
		 const struct rt_view_lod_policy *policy,
		 int mesh_enabled,
		 int csg_enabled,
		 int zoom_refresh,
		 fastf_t scale,
		 fastf_t curve_scale,
		 fastf_t point_scale,
		 size_t bot_threshold)
{
    if (!policy) {
	printf("FAIL: %s policy is NULL\n", label);
	return 1;
    }

    if (policy->mesh_enabled != mesh_enabled ||
	    policy->csg_enabled != csg_enabled ||
	    policy->zoom_refresh != zoom_refresh ||
	    !fastf_equal(policy->scale, scale) ||
	    !fastf_equal(policy->curve_scale, curve_scale) ||
	    !fastf_equal(policy->point_scale, point_scale) ||
	    policy->bot_threshold != bot_threshold) {
	printf("FAIL: %s policy: got enabled=%d/%d refresh=%d scale=%g curve=%g point=%g threshold=%zu\n",
		label, policy->mesh_enabled, policy->csg_enabled,
		policy->zoom_refresh, policy->scale, policy->curve_scale,
		policy->point_scale, policy->bot_threshold);
	return 1;
    }

    return 0;
}

static struct bsg_view *
make_view(const char *name)
{
    struct bsg_view *v;

    BU_ALLOC(v, struct bsg_view);
    bsg_view_init(v, NULL);
    bu_vls_sprintf(&v->gv_name, "%s", name);
    return v;
}

static void
free_view(struct bsg_view *v)
{
    if (!v)
	return;

    bsg_view_free(v);
    bu_free(v, "view legacy bsg test view");
}

static int test_bounds_update_count = 0;
static int test_update_callback_count = 0;
static int test_neutral_update_callback_count = 0;
static void *test_neutral_update_callback_view = NULL;

static void
test_bounds_update_callback(struct bsg_view *UNUSED(v))
{
    test_bounds_update_count++;
}

static void
test_update_callback(struct bsg_view *UNUSED(v), void *data)
{
    int *token = (int *)data;
    test_update_callback_count++;
    if (token)
	(*token)++;
}

static void
test_neutral_update_callback(void *view_ctx, void *data)
{
    int *token = (int *)data;
    test_neutral_update_callback_count++;
    test_neutral_update_callback_view = view_ctx;
    if (token)
	(*token)++;
}

struct test_selection_path_callback_state {
    int count;
    char last_path[64];
};

static void
test_selection_path_callback(const char *path, void *data)
{
    struct test_selection_path_callback_state *state =
	(struct test_selection_path_callback_state *)data;

    if (!state)
	return;
    state->count++;
    bu_strlcpy(state->last_path, path ? path : "", sizeof(state->last_path));
}

struct test_polygon_record_callback_state {
    int count;
    rt_view_polygon_ref last_ref;
};

static int
test_polygon_record_callback(rt_view_polygon_ref ref,
			     const struct rt_view_polygon_record *record,
			     void *data)
{
    struct test_polygon_record_callback_state *state =
	(struct test_polygon_record_callback_state *)data;

    if (!state || !record || rt_view_polygon_ref_is_null(ref))
	return 0;
    state->count++;
    state->last_ref = ref;
    return 1;
}

static bsg_polygon_ref
test_polygon_ref_to_bsg(rt_view_polygon_ref ref)
{
    bsg_polygon_ref bsg_ref = {ref.token, ref.revision};
    return bsg_ref;
}

struct test_feature_preview_callback_state {
    uint64_t revision;
    int pick_calls;
    int last_x;
    int last_y;
    int pick_value;
};

static uint64_t
test_feature_preview_revision(void *data)
{
    struct test_feature_preview_callback_state *state =
	(struct test_feature_preview_callback_state *)data;
    return state ? state->revision : 0;
}

static int
test_feature_preview_pick(void *data, int x, int y, void *pick_out)
{
    struct test_feature_preview_callback_state *state =
	(struct test_feature_preview_callback_state *)data;
    if (!state || !pick_out)
	return 0;

    state->pick_calls++;
    state->last_x = x;
    state->last_y = y;
    *(int *)pick_out = state->pick_value;
    return 1;
}

static bsg_feature_ref
test_feature_ref_to_bsg(rt_view_feature_ref ref)
{
    bsg_feature_ref bsg_ref = {ref.token, ref.revision};
    return bsg_ref;
}

static int
test_bsg_adapter_sanitizes(void)
{
    struct bsg_view *v = make_view("rt_view_info_adapter_sanitize");
    struct bsg_lod_source_policy_settings policy;
    struct rt_view_info info;
    int ret = 0;

    v->gv_width = 0;
    v->gv_height = -64;
    v->gv_size = 0.0;

    if (!bsg_view_lod_source_policy_get(v, &policy)) {
	printf("FAIL: bsg lod policy get\n");
	free_view(v);
	return 1;
    }

    policy.scale = 0.0;
    policy.curve_scale = -3.0;
    policy.point_scale = 0.0;
    policy.bot_threshold = 77;
    bsg_view_lod_source_policy_set(v, &policy);

    rt_view_info_from_bsg(&info, v);
    ret = check_view_info("bsg adapter sanitize", &info, 1, 1, 1.0,
	    1.0, 1.0, 1.0, 77);

    free_view(v);
    return ret;
}

static int
test_bsg_null_adapter(void)
{
    struct rt_view_info info;
    point_t center = VINIT_ZERO;
    vect_t extent1 = VINIT_ZERO;
    vect_t extent2 = VINIT_ZERO;
    vect_t extent3 = VINIT_ZERO;

    info.width = 640;
    info.height = 480;
    info.size = 1000.0;
    info.lod.scale = 4.0;
    info.lod.curve_scale = 5.0;
    info.lod.point_scale = 6.0;
    info.lod.bot_threshold = 7;

    rt_view_info_from_bsg(NULL, NULL);
    rt_view_info_from_bsg(&info, NULL);
    rt_view_context_info_from_bsg(NULL, NULL);
    rt_view_context_info_from_bsg(&info, NULL);
    rt_view_context_info_get(NULL, NULL);
    rt_view_context_info_get(&info, NULL);

    if (rt_view_context_obb_from_bsg(NULL, center, extent1, extent2, extent3) ||
	    rt_view_context_obb_from_bsg(&info, NULL, extent1, extent2, extent3) ||
	    rt_view_context_obb_from_bsg(&info, center, NULL, extent2, extent3) ||
	    rt_view_context_obb_from_bsg(&info, center, extent1, NULL, extent3) ||
	    rt_view_context_obb_from_bsg(&info, center, extent1, extent2, NULL) ||
	    rt_view_context_obb_get(NULL, center, extent1, extent2, extent3) ||
	    rt_view_context_obb_get(&info, NULL, extent1, extent2, extent3) ||
	    rt_view_context_obb_get(&info, center, NULL, extent2, extent3) ||
	    rt_view_context_obb_get(&info, center, extent1, NULL, extent3) ||
	    rt_view_context_obb_get(&info, center, extent1, extent2, NULL)) {
	printf("FAIL: null BSG OBB adapter\n");
	return 1;
    }

    return check_view_info("null bsg adapter", &info, 1, 1, 1.0,
	    1.0, 1.0, 1.0, 0);
}

static int
test_bsg_obb_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_obb_adapter");
    point_t center = VINIT_ZERO;
    vect_t extent1 = VINIT_ZERO;
    vect_t extent2 = VINIT_ZERO;
    vect_t extent3 = VINIT_ZERO;
    point_t expected_center = VINIT_ZERO;
    vect_t expected_extent1 = VINIT_ZERO;
    vect_t expected_extent2 = VINIT_ZERO;
    vect_t expected_extent3 = VINIT_ZERO;
    int ret = 0;

    VSET(expected_center, 1.0, 2.0, 3.0);
    VSET(expected_extent1, 4.0, 0.0, 0.0);
    VSET(expected_extent2, 0.0, 5.0, 0.0);
    VSET(expected_extent3, 0.0, 0.0, 6.0);
    VMOVE(v->obb_center, expected_center);
    VMOVE(v->obb_extent1, expected_extent1);
    VMOVE(v->obb_extent2, expected_extent2);
    VMOVE(v->obb_extent3, expected_extent3);

    if (!rt_view_context_obb_from_bsg(v, center, extent1, extent2, extent3) ||
	    !VNEAR_EQUAL(center, expected_center, SMALL_FASTF) ||
	    !VNEAR_EQUAL(extent1, expected_extent1, SMALL_FASTF) ||
	    !VNEAR_EQUAL(extent2, expected_extent2, SMALL_FASTF) ||
	    !VNEAR_EQUAL(extent3, expected_extent3, SMALL_FASTF)) {
	printf("FAIL: BSG OBB adapter did not copy expected values\n");
	ret = 1;
    }

    VSETALL(center, 0.0);
    VSETALL(extent1, 0.0);
    VSETALL(extent2, 0.0);
    VSETALL(extent3, 0.0);

    if (!rt_view_context_obb_get(v, center, extent1, extent2, extent3) ||
	    !VNEAR_EQUAL(center, expected_center, SMALL_FASTF) ||
	    !VNEAR_EQUAL(extent1, expected_extent1, SMALL_FASTF) ||
	    !VNEAR_EQUAL(extent2, expected_extent2, SMALL_FASTF) ||
	    !VNEAR_EQUAL(extent3, expected_extent3, SMALL_FASTF)) {
	printf("FAIL: neutral OBB adapter did not copy expected values\n");
	ret = 1;
    }

    free_view(v);
    return ret;
}

static int
test_bsg_user_data_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_user_data_adapter");
    struct bsg_data_tclcad tcl_data;
    void *tcl_data_ctx = &tcl_data;
    struct rt_view_info bogus;
    int marker = 42;
    int ret = 0;

    memset(&bogus, 0, sizeof(bogus));
    memset(tcl_data_ctx, 0, sizeof(tcl_data));

    if (rt_view_context_user_data_from_bsg(NULL) ||
	    rt_view_context_user_data_get(NULL) ||
	    rt_view_context_user_data_from_bsg(&bogus) ||
	    rt_view_context_user_data_get(&bogus) ||
	    rt_view_context_user_data_set_bsg(NULL, &marker) ||
	    rt_view_context_user_data_set(NULL, &marker) ||
	    rt_view_context_user_data_set_bsg(&bogus, &marker) ||
	    rt_view_context_user_data_set(&bogus, &marker) ||
	    rt_view_context_tclcad_data_set_bsg(NULL, tcl_data_ctx) ||
	    rt_view_context_tclcad_data_set(NULL, tcl_data_ctx) ||
	    rt_view_context_tclcad_data_set_bsg(&bogus, tcl_data_ctx) ||
	    rt_view_context_tclcad_data_set(&bogus, tcl_data_ctx)) {
	printf("FAIL: null/non-BSG user-data adapter accepted data\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_user_data_set_bsg(v, &marker) ||
	    v->u_data != &marker ||
	    rt_view_context_user_data_from_bsg(v) != &marker ||
	    rt_view_context_user_data_get(v) != &marker) {
	printf("FAIL: BSG user-data adapter did not return expected pointer\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_user_data_set(v, &marker) ||
	    v->u_data != &marker ||
	    rt_view_context_user_data_get(v) != &marker) {
	printf("FAIL: retained context user-data adapter did not return expected pointer\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_tclcad_data_set_bsg(v, tcl_data_ctx) ||
	    !rt_view_context_tclcad_data_set(v, tcl_data_ctx) ||
	    v->gv_tcl != &tcl_data) {
	printf("FAIL: BSG TclCAD data adapter did not bind expected pointer\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_user_data_set_bsg(v, NULL) ||
	    !rt_view_context_user_data_set(v, NULL) ||
	    !rt_view_context_tclcad_data_set_bsg(v, NULL) ||
	    !rt_view_context_tclcad_data_set(v, NULL) ||
	    rt_view_context_user_data_from_bsg(v) ||
	    rt_view_context_user_data_get(v) ||
	    v->gv_tcl) {
	printf("FAIL: BSG user/TclCAD data adapters did not clear pointers\n");
	ret = 1;
    }

cleanup:
    free_view(v);
    return ret;
}

static int
test_bsg_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_info_adapter");
    struct bsg_lod_source_policy_settings policy;
    struct rt_view_info info;
    int ret = 0;

    v->gv_width = 1280;
    v->gv_height = 720;
    v->gv_size = 250.0;

    if (!bsg_view_lod_source_policy_get(v, &policy)) {
	printf("FAIL: bsg lod policy get\n");
	free_view(v);
	return 1;
    }

    policy.scale = 2.0;
    policy.curve_scale = 3.0;
    policy.point_scale = 4.0;
    policy.bot_threshold = 12345;
    bsg_view_lod_source_policy_set(v, &policy);

    rt_view_info_from_bsg(&info, v);
    ret = check_view_info("bsg adapter", &info, 1280, 720, 250.0,
	    2.0, 3.0, 4.0, 12345);
    if (!ret) {
	rt_view_context_info_from_bsg(&info, v);
	ret = check_view_info("bsg context adapter", &info, 1280, 720, 250.0,
		2.0, 3.0, 4.0, 12345);
    }
    if (!ret) {
	rt_view_context_info_get(&info, v);
	ret = check_view_info("retained context adapter", &info, 1280, 720, 250.0,
		2.0, 3.0, 4.0, 12345);
    }
    if (!ret && (!BU_STR_EQUAL(rt_view_name_from_bsg(v), "rt_view_info_adapter") ||
	    !BU_STR_EQUAL(rt_view_context_name_from_bsg(v), "rt_view_info_adapter"))) {
	printf("FAIL: BSG view-name adapter\n");
	ret = 1;
    }
    bu_vls_trunc(&v->gv_name, 0);
    if (!ret && (rt_view_name_from_bsg(v) || rt_view_context_name_from_bsg(v))) {
	printf("FAIL: empty BSG view-name adapter\n");
	ret = 1;
    }
    if (!ret && (rt_view_name_from_bsg(NULL) || rt_view_context_name_from_bsg(NULL))) {
	printf("FAIL: null BSG view-name adapter\n");
	ret = 1;
    }

    free_view(v);
    return ret;
}

static int
test_bsg_orientation_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_orientation_adapter");
    quat_t actual = HINIT_ZERO;
    quat_t expected = HINIT_ZERO;
    mat_t identity;
    int ret = 0;

    if (rt_view_orientation_quat_from_bsg(actual, NULL) ||
	    rt_view_context_orientation_quat_from_bsg(actual, NULL) ||
	    rt_view_context_orientation_quat_get(actual, NULL)) {
	printf("FAIL: null BSG orientation adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(identity);
    quat_mat2quat(expected, identity);
    ret += check_quat("null BSG orientation adapter", actual, expected);

    bn_mat_angles(v->gv_rotation, 12.0, 34.0, 56.0);
    quat_mat2quat(expected, v->gv_rotation);
    if (!rt_view_orientation_quat_from_bsg(actual, v)) {
	printf("FAIL: BSG orientation adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_quat("BSG orientation adapter", actual, expected);
    HSETALL(actual, 0.0);
    if (!rt_view_context_orientation_quat_from_bsg(actual, v)) {
	printf("FAIL: BSG context orientation adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_quat("BSG context orientation adapter", actual, expected);
    HSETALL(actual, 0.0);
    if (!rt_view_context_orientation_quat_get(actual, v)) {
	printf("FAIL: retained context orientation adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_quat("retained context orientation adapter", actual, expected);

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_perspective_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_perspective_adapter");
    int ret = 0;

    if (!fastf_equal(rt_view_perspective_from_bsg(NULL), 0.0) ||
	    !fastf_equal(rt_view_context_perspective_from_bsg(NULL), 0.0) ||
	    !fastf_equal(rt_view_context_perspective_get(NULL), 0.0)) {
	printf("FAIL: null BSG perspective adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_perspective = 37.5;
    if (!fastf_equal(rt_view_perspective_from_bsg(v), 37.5) ||
	    !fastf_equal(rt_view_context_perspective_from_bsg(v), 37.5) ||
	    !fastf_equal(rt_view_context_perspective_get(v), 37.5)) {
	printf("FAIL: BSG perspective adapter got %g\n",
		rt_view_perspective_from_bsg(v));
	ret = 1;
	goto cleanup;
    }

    if (rt_view_perspective_set_bsg(NULL, 12.5) ||
	    rt_view_context_perspective_set_bsg(NULL, 12.5) ||
	    rt_view_context_perspective_set(NULL, 12.5)) {
	printf("FAIL: null BSG perspective set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_perspective_set_bsg(v, 12.5) ||
	    !fastf_equal(rt_view_perspective_from_bsg(v), 12.5) ||
	    !fastf_equal(rt_view_context_perspective_from_bsg(v), 12.5) ||
	    !fastf_equal(rt_view_context_perspective_get(v), 12.5)) {
	printf("FAIL: BSG perspective set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_perspective_set_bsg(v, 22.5) ||
	    !fastf_equal(rt_view_perspective_from_bsg(v), 22.5) ||
	    !fastf_equal(rt_view_context_perspective_from_bsg(v), 22.5) ||
	    !fastf_equal(rt_view_context_perspective_get(v), 22.5)) {
	printf("FAIL: BSG context perspective set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_perspective_set(v, 32.5) ||
	    !fastf_equal(rt_view_context_perspective_get(v), 32.5) ||
	    !fastf_equal(rt_view_perspective_from_bsg(v), 32.5)) {
	printf("FAIL: retained context perspective set adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_display_manager_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_display_manager_adapter");
    int token = 0;
    void *dmp = &token;
    int ret = 0;

    if (rt_view_display_manager_from_bsg(NULL) ||
	    rt_view_context_display_manager_from_bsg(NULL) ||
	    rt_view_context_display_manager_get(NULL)) {
	printf("FAIL: null BSG display-manager adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_display_manager_set_bsg(NULL, dmp) ||
	    rt_view_context_display_manager_set_bsg(NULL, dmp) ||
	    rt_view_context_display_manager_set(NULL, dmp)) {
	printf("FAIL: null BSG display-manager set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_display_manager_set_bsg(v, dmp) ||
	    rt_view_display_manager_from_bsg(v) != dmp ||
	    rt_view_context_display_manager_from_bsg(v) != dmp ||
	    rt_view_context_display_manager_get(v) != dmp ||
	    v->dmp != dmp) {
	printf("FAIL: BSG display-manager set/get adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_display_manager_set_bsg(v, NULL) ||
	    rt_view_display_manager_from_bsg(v) != NULL ||
	    rt_view_context_display_manager_from_bsg(v) != NULL ||
	    rt_view_context_display_manager_get(v) != NULL ||
	    v->dmp != NULL) {
	printf("FAIL: BSG display-manager clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_display_manager_set_bsg(v, dmp) ||
	    rt_view_display_manager_from_bsg(v) != dmp ||
	    rt_view_context_display_manager_from_bsg(v) != dmp ||
	    rt_view_context_display_manager_get(v) != dmp ||
	    v->dmp != dmp) {
	printf("FAIL: BSG context display-manager set/get adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_display_manager_set_bsg(v, NULL) ||
	    rt_view_display_manager_from_bsg(v) != NULL ||
	    rt_view_context_display_manager_from_bsg(v) != NULL ||
	    rt_view_context_display_manager_get(v) != NULL ||
	    v->dmp != NULL) {
	printf("FAIL: BSG context display-manager clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_display_manager_set(v, dmp) ||
	    rt_view_display_manager_from_bsg(v) != dmp ||
	    rt_view_context_display_manager_get(v) != dmp ||
	    v->dmp != dmp) {
	printf("FAIL: retained context display-manager set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_display_manager_set(v, NULL) ||
	    rt_view_display_manager_from_bsg(v) != NULL ||
	    rt_view_context_display_manager_get(v) != NULL ||
	    v->dmp != NULL) {
	printf("FAIL: retained context display-manager clear adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_view_runtime_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_runtime_adapter");
    int token = 0;
    mat_t edit_mat;
    int ret = 0;

    MAT_IDN(edit_mat);

    if (rt_view_update_callback_set_bsg(NULL, test_update_callback, &token) ||
	    rt_view_context_update_callback_set_bsg(NULL,
		test_update_callback, &token) ||
	    rt_view_context_update_callback_set(NULL,
		test_neutral_update_callback, &token)) {
	printf("FAIL: null BSG update callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_edit_matrix_set_bsg(NULL, edit_mat)) {
	printf("FAIL: null BSG edit matrix set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_edit_matrix_set_bsg(NULL, edit_mat)) {
	printf("FAIL: null BSG edit matrix context set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_edit_matrix_set(NULL, edit_mat)) {
	printf("FAIL: null retained edit matrix context set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_edit_matrix_clear_bsg(NULL)) {
	printf("FAIL: null BSG edit matrix clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_edit_matrix_clear_bsg(NULL)) {
	printf("FAIL: null BSG edit matrix context clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_edit_matrix_clear(NULL)) {
	printf("FAIL: null retained edit matrix context clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_frame_revision_from_bsg(NULL)) {
	printf("FAIL: null BSG frame revision adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_frame_revision_from_bsg(NULL)) {
	printf("FAIL: null BSG frame revision context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_frame_revision_get(NULL)) {
	printf("FAIL: null retained frame revision context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_frame_revision_bump(NULL)) {
	printf("FAIL: null retained frame revision bump adapter\n");
	ret = 1;
	goto cleanup;
    }

    test_update_callback_count = 0;
    if (!rt_view_context_update_callback_set_bsg(v, test_update_callback, &token) ||
	    v->gv_callback != test_update_callback ||
	    v->gv_clientData != &token ||
	    !v->callbacks) {
	printf("FAIL: BSG update callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_update(v);
    if (test_update_callback_count != 1 || token != 1) {
	printf("FAIL: BSG update callback adapter did not fire\n");
	ret = 1;
	goto cleanup;
    }

    test_neutral_update_callback_count = 0;
    test_neutral_update_callback_view = NULL;
    token = 0;
    if (!rt_view_context_update_callback_set(v, test_neutral_update_callback,
		&token) ||
	    !v->gv_callback ||
	    v->gv_callback == test_update_callback ||
	    v->gv_clientData == &token ||
	    !v->callbacks) {
	printf("FAIL: retained update callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_update(v);
    if (test_neutral_update_callback_count != 1 ||
	    test_neutral_update_callback_view != v ||
	    token != 1) {
	printf("FAIL: retained update callback adapter did not fire\n");
	ret = 1;
	goto cleanup;
    }

    token = 0;
    test_update_callback_count = 0;
    if (!rt_view_context_update_callback_set_bsg(v, test_update_callback,
		&token) ||
	    v->gv_callback != test_update_callback ||
	    v->gv_clientData != &token) {
	printf("FAIL: BSG update callback set after retained adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_update(v);
    if (test_update_callback_count != 1 || token != 1) {
	printf("FAIL: BSG update callback after retained adapter did not fire\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_update_callback_set(v, NULL, NULL) ||
	    v->gv_callback || v->gv_clientData) {
	printf("FAIL: retained update callback clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_update_callback_set_bsg(v, NULL, NULL) ||
	    v->gv_callback || v->gv_clientData) {
	printf("FAIL: BSG update callback clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_edit_matrix_set_bsg(v, edit_mat) ||
	    v->gv_edit_mat != edit_mat) {
	printf("FAIL: BSG edit matrix set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_edit_matrix_clear_bsg(v) || v->gv_edit_mat) {
	printf("FAIL: BSG edit matrix clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_edit_matrix_set_bsg(v, edit_mat) ||
	    v->gv_edit_mat != edit_mat) {
	printf("FAIL: BSG edit matrix context set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_edit_matrix_clear_bsg(v) || v->gv_edit_mat) {
	printf("FAIL: BSG edit matrix context clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_edit_matrix_set(v, edit_mat) ||
	    v->gv_edit_mat != edit_mat) {
	printf("FAIL: retained edit matrix context set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_edit_matrix_clear(v) || v->gv_edit_mat) {
	printf("FAIL: retained edit matrix context clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_frame_rev = 42;
    if (rt_view_frame_revision_from_bsg(v) != 42) {
	printf("FAIL: BSG frame revision adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_frame_revision_from_bsg(v) != 42) {
	printf("FAIL: BSG frame revision context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_frame_revision_get(v) != 42) {
	printf("FAIL: retained frame revision context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_frame_revision_bump(v) != 43 ||
	    rt_view_context_frame_revision_get(v) != 43) {
	printf("FAIL: retained frame revision bump adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_aet_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_aet_adapter");
    vect_t actual = VINIT_ZERO;
    vect_t expected = VINIT_ZERO;
    int ret = 0;

    VSET(actual, 9.0, 8.0, 7.0);
    if (rt_view_aet_from_bsg(actual, NULL) ||
	    rt_view_context_aet_from_bsg(actual, NULL) ||
	    rt_view_context_aet_get(actual, NULL)) {
	printf("FAIL: null BSG AET adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    VSETALL(expected, 0.0);
    if (!VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: null BSG AET adapter did not return zeros\n");
	ret = 1;
	goto cleanup;
    }

    VSET(v->gv_aet, 12.0, 34.0, 56.0);
    if (!rt_view_aet_from_bsg(actual, v) ||
	    !rt_view_context_aet_from_bsg(actual, v) ||
	    !rt_view_context_aet_get(actual, v)) {
	printf("FAIL: BSG AET adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (!VNEAR_EQUAL(actual, v->gv_aet, SMALL_FASTF)) {
	printf("FAIL: BSG AET adapter did not copy vector\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_aet_set_bsg(NULL, actual) ||
	    rt_view_aet_set_bsg(v, NULL) ||
	    rt_view_context_aet_set_bsg(NULL, actual) ||
	    rt_view_context_aet_set_bsg(v, NULL) ||
	    rt_view_context_aet_set(NULL, actual) ||
	    rt_view_context_aet_set(v, NULL)) {
	printf("FAIL: null BSG AET set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected, 21.0, 22.0, 23.0);
    if (!rt_view_aet_set_bsg(v, expected) ||
	    !rt_view_context_aet_set_bsg(v, expected) ||
	    !rt_view_aet_from_bsg(actual, v) ||
	    !VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: BSG AET set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected, 24.0, 25.0, 26.0);
    if (!rt_view_context_aet_set(v, expected) ||
	    !rt_view_context_aet_get(actual, v) ||
	    !VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: retained context AET set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_aet_state_set_bsg(NULL, expected) ||
	    rt_view_aet_state_set_bsg(v, NULL) ||
	    rt_view_context_aet_state_set_bsg(NULL, expected) ||
	    rt_view_context_aet_state_set_bsg(v, NULL) ||
	    rt_view_context_aet_state_set(NULL, expected) ||
	    rt_view_context_aet_state_set(v, NULL)) {
	printf("FAIL: null BSG AET state set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected, 31.0, 32.0, 33.0);
    if (!rt_view_aet_state_set_bsg(v, expected) ||
	    !rt_view_aet_from_bsg(actual, v) ||
	    !VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: BSG AET state set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected, 41.0, 42.0, 43.0);
    if (!rt_view_context_aet_state_set_bsg(v, expected) ||
	    !rt_view_aet_from_bsg(actual, v) ||
	    !VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: BSG AET context state set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected, 51.0, 52.0, 53.0);
    if (!rt_view_context_aet_state_set(v, expected) ||
	    !rt_view_context_aet_get(actual, v) ||
	    !VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: retained context AET state set adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_update_adapter(void)
{
    struct bsg_view *direct = make_view("rt_view_update_direct");
    struct bsg_view *adapter = make_view("rt_view_update_adapter");
    int ret = 0;

    if (rt_view_update_bsg(NULL) ||
	    rt_view_context_update_bsg(NULL) ||
	    rt_view_context_update(NULL)) {
	printf("FAIL: null BSG update adapter\n");
	ret = 1;
	goto cleanup;
    }

    bn_mat_angles(direct->gv_rotation, 11.0, 22.0, 33.0);
    bn_mat_angles(adapter->gv_rotation, 11.0, 22.0, 33.0);
    MAT_DELTAS(direct->gv_center, 7.0, 8.0, 9.0);
    MAT_DELTAS(adapter->gv_center, 7.0, 8.0, 9.0);
    direct->gv_scale = 4.5;
    adapter->gv_scale = 4.5;
    bn_mat_angles(direct->gv_pmat, 1.0, 2.0, 3.0);
    bn_mat_angles(adapter->gv_pmat, 1.0, 2.0, 3.0);

    bsg_update(direct);
    if (!rt_view_update_bsg(adapter) ||
	    !rt_view_context_update_bsg(adapter) ||
	    !rt_view_context_update(adapter)) {
	printf("FAIL: BSG update adapter should report success\n");
	ret = 1;
	goto cleanup;
    }

    if (check_mat("BSG update adapter model2view",
		adapter->gv_model2view, direct->gv_model2view) ||
	    check_mat("BSG update adapter view2model",
		adapter->gv_view2model, direct->gv_view2model) ||
	    check_mat("BSG update adapter pmodel2view",
		adapter->gv_pmodel2view, direct->gv_pmodel2view)) {
	ret = 1;
	goto cleanup;
    }
    if (!VNEAR_EQUAL(adapter->gv_aet, direct->gv_aet, SMALL_FASTF)) {
	printf("FAIL: BSG update adapter AET: got %g,%g,%g expected %g,%g,%g\n",
		adapter->gv_aet[0], adapter->gv_aet[1], adapter->gv_aet[2],
		direct->gv_aet[0], direct->gv_aet[1], direct->gv_aet[2]);
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(direct);
    free_view(adapter);
    return ret ? 1 : 0;
}

static int
test_bsg_autoview_adapter(void)
{
    struct bsg_view *direct_bounds = make_view("rt_view_autoview_bounds_direct");
    struct bsg_view *adapter_bounds = make_view("rt_view_autoview_bounds_adapter");
    struct bsg_view *context_bounds = make_view("rt_view_autoview_bounds_context");
    struct bsg_view *direct_empty = make_view("rt_view_autoview_empty_direct");
    struct bsg_view *adapter_empty = make_view("rt_view_autoview_empty_adapter");
    point_t min, max;
    int ret = 0;

    VSET(min, -1.0, -2.0, -3.0);
    VSET(max, 5.0, 6.0, 7.0);

    if (rt_view_autoview_bsg(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0) ||
	    rt_view_context_autoview_bsg(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0) ||
	    rt_view_context_autoview(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0) ||
	    rt_view_autoview_bounds_bsg(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max) ||
	    rt_view_context_autoview_bounds_bsg(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max) ||
	    rt_view_context_autoview_bounds(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max) ||
	    rt_view_autoview_bounds_bsg(adapter_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, NULL, max) ||
	    rt_view_autoview_bounds_bsg(adapter_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, NULL) ||
	    rt_view_context_autoview_bounds_bsg(context_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, NULL, max) ||
	    rt_view_context_autoview_bounds_bsg(context_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, NULL) ||
	    rt_view_context_autoview_bounds(context_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, NULL, max) ||
	    rt_view_context_autoview_bounds(context_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, NULL)) {
	printf("FAIL: null BSG autoview adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    bsg_autoview_bounds(direct_bounds, BSG_AUTOVIEW_SCALE_DEFAULT, min, max);
    if (!rt_view_autoview_bounds_bsg(adapter_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max)) {
	printf("FAIL: BSG autoview-bounds adapter should report success\n");
	ret = 1;
	goto cleanup;
    }
    if (check_view_frame("BSG autoview-bounds adapter", adapter_bounds, direct_bounds)) {
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_autoview_bounds_bsg(context_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max)) {
	printf("FAIL: BSG context autoview-bounds adapter should report success\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_autoview_bounds(context_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max)) {
	printf("FAIL: retained context autoview-bounds adapter should report success\n");
	ret = 1;
	goto cleanup;
    }
    if (check_view_frame("BSG context autoview-bounds adapter", context_bounds, direct_bounds)) {
	ret = 1;
	goto cleanup;
    }

    bsg_autoview(direct_empty, BSG_AUTOVIEW_SCALE_DEFAULT, 0);
    if (!rt_view_autoview_bsg(adapter_empty, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0)) {
	printf("FAIL: BSG autoview adapter should report success\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_autoview_bsg(adapter_empty, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0)) {
	printf("FAIL: BSG context autoview adapter should report success\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_autoview(adapter_empty, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0)) {
	printf("FAIL: retained context autoview adapter should report success\n");
	ret = 1;
	goto cleanup;
    }
    if (check_view_frame("BSG autoview adapter", adapter_empty, direct_empty)) {
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(direct_bounds);
    free_view(adapter_bounds);
    free_view(context_bounds);
    free_view(direct_empty);
    free_view(adapter_empty);
    return ret ? 1 : 0;
}

static void
setup_adjust_view(struct bsg_view *v)
{
    v->gv_width = 640;
    v->gv_height = 480;
    v->gv_scale = 12.0;
    v->gv_size = 2.0 * v->gv_scale;
    v->gv_isize = 1.0 / v->gv_size;
    bn_mat_angles(v->gv_rotation, 3.0, 5.0, 7.0);
    MAT_DELTAS(v->gv_center, 2.0, 4.0, 6.0);
    bsg_update(v);
}

static int
check_adjust_equivalence(const char *label,
			 int dx,
			 int dy,
			 const point_t keypoint,
			 unsigned long long bsg_flags,
			 unsigned long long rt_flags)
{
    struct bsg_view *direct = make_view(label);
    struct bsg_view *adapter = make_view(label);
    struct bsg_view *context = make_view(label);
    struct bsg_view *neutral = make_view(label);
    point_t direct_keypoint;
    point_t adapter_keypoint;
    point_t context_keypoint;
    point_t neutral_keypoint;
    int ret = 0;

    setup_adjust_view(direct);
    setup_adjust_view(adapter);
    setup_adjust_view(context);
    setup_adjust_view(neutral);
    VMOVE(direct_keypoint, keypoint);
    VMOVE(adapter_keypoint, keypoint);
    VMOVE(context_keypoint, keypoint);
    VMOVE(neutral_keypoint, keypoint);

    if (!bsg_adjust(direct, dx, dy, direct_keypoint, 0, bsg_flags)) {
	printf("FAIL: direct BSG adjust %s should report success\n", label);
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_adjust_bsg(adapter, dx, dy, adapter_keypoint, 0, rt_flags)) {
	printf("FAIL: BSG adjust adapter %s should report success\n", label);
	ret = 1;
	goto cleanup;
    }
    if (check_view_frame(label, adapter, direct))
	ret = 1;
    if (!rt_view_context_adjust_bsg(context, dx, dy, context_keypoint, 0,
	    rt_flags)) {
	printf("FAIL: BSG context adjust adapter %s should report success\n", label);
	ret = 1;
	goto cleanup;
    }
    if (check_view_frame(label, context, direct))
	ret = 1;
    if (!rt_view_context_adjust(neutral, dx, dy, neutral_keypoint, 0,
	    rt_flags)) {
	printf("FAIL: retained context adjust adapter %s should report success\n", label);
	ret = 1;
	goto cleanup;
    }
    if (check_view_frame(label, neutral, direct))
	ret = 1;

cleanup:
    free_view(direct);
    free_view(adapter);
    free_view(context);
    free_view(neutral);
    return ret ? 1 : 0;
}

static int
test_bsg_adjust_adapter(void)
{
    point_t keypoint;
    point_t origin = VINIT_ZERO;
    struct bsg_view *v = make_view("rt_view_adjust_nulls");
    int ret = 0;

    VSET(keypoint, 1.0, 2.0, 3.0);

    if (rt_view_adjust_bsg(NULL, 1, 2, keypoint, 0, RT_VIEW_ADJUST_SCALE) ||
	    rt_view_context_adjust_bsg(NULL, 1, 2, keypoint, 0,
		RT_VIEW_ADJUST_SCALE) ||
	    rt_view_context_adjust(NULL, 1, 2, keypoint, 0,
		RT_VIEW_ADJUST_SCALE) ||
	    rt_view_adjust_bsg(v, 1, 2, NULL, 0, RT_VIEW_ADJUST_SCALE) ||
	    rt_view_context_adjust_bsg(v, 1, 2, NULL, 0,
		RT_VIEW_ADJUST_SCALE) ||
	    rt_view_context_adjust(v, 1, 2, NULL, 0,
		RT_VIEW_ADJUST_SCALE) ||
	    rt_view_adjust_bsg(v, 1, 2, keypoint, 0, RT_VIEW_ADJUST_IDLE) ||
	    rt_view_context_adjust_bsg(v, 1, 2, keypoint, 0,
		RT_VIEW_ADJUST_IDLE) ||
	    rt_view_context_adjust(v, 1, 2, keypoint, 0,
		RT_VIEW_ADJUST_IDLE)) {
	printf("FAIL: null/idle BSG adjust adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    if (check_adjust_equivalence("BSG adjust scale adapter",
		100, 125, origin, BSG_SCALE, RT_VIEW_ADJUST_SCALE) ||
	    check_adjust_equivalence("BSG adjust translate adapter",
		12, -9, origin, BSG_TRANS, RT_VIEW_ADJUST_TRANS) ||
	    check_adjust_equivalence("BSG adjust rotate adapter",
		7, -5, keypoint, BSG_ROT, RT_VIEW_ADJUST_ROT) ||
	    check_adjust_equivalence("BSG adjust center adapter",
		320, 240, origin, BSG_CENTER, RT_VIEW_ADJUST_CENTER)) {
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_hash_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_hash_adapter");
    unsigned long long direct_hash = 0ULL;
    unsigned long long adapter_hash = 0ULL;
    int ret = 0;

    if (rt_view_hash_bsg(NULL) != 0ULL ||
	    rt_view_context_hash_bsg(NULL) != 0ULL ||
	    rt_view_context_hash(NULL) != 0ULL) {
	printf("FAIL: null BSG hash adapter\n");
	ret = 1;
	goto cleanup;
    }

    setup_adjust_view(v);
    direct_hash = bsg_hash(v);
    adapter_hash = rt_view_hash_bsg(v);
    if (direct_hash != adapter_hash ||
	    direct_hash != rt_view_context_hash_bsg(v) ||
	    direct_hash != rt_view_context_hash(v)) {
	printf("FAIL: BSG hash adapter: got %llu expected %llu\n",
		adapter_hash, direct_hash);
	ret = 1;
	goto cleanup;
    }

    v->gv_width += 7;
    direct_hash = bsg_hash(v);
    adapter_hash = rt_view_hash_bsg(v);
    if (direct_hash != adapter_hash ||
	    direct_hash != rt_view_context_hash_bsg(v) ||
	    direct_hash != rt_view_context_hash(v)) {
	printf("FAIL: BSG hash adapter after mutation: got %llu expected %llu\n",
		adapter_hash, direct_hash);
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_knob_adapter(void)
{
    struct bsg_view *direct = make_view("rt_view_knob_direct");
    struct bsg_view *adapter = make_view("rt_view_knob_adapter");
    vect_t direct_rvec = VINIT_ZERO;
    vect_t adapter_rvec = VINIT_ZERO;
    vect_t direct_tvec = VINIT_ZERO;
    vect_t adapter_tvec = VINIT_ZERO;
    int direct_do_rot = 0;
    int adapter_do_rot = 0;
    int direct_do_tran = 0;
    int adapter_do_tran = 0;
    struct rt_view_knob_values values;
    struct rt_view_knobs neutral_knobs;
    int ret = 0;

    if (rt_view_knobs_reset_bsg(NULL, RT_VIEW_KNOBS_RATE_BSG) ||
	    rt_view_context_knobs_reset_bsg(NULL, RT_VIEW_KNOBS_RATE_BSG) ||
	    rt_view_context_knobs_reset(NULL, RT_VIEW_KNOBS_RATE) ||
	    rt_view_knobs_state_reset(NULL, RT_VIEW_KNOBS_RATE) ||
	    rt_view_knobs_state_hash(NULL, NULL) ||
	    rt_view_context_knobs_state_from_bsg(NULL, adapter) ||
	    rt_view_context_knobs_state_get(NULL, adapter) ||
	    rt_view_context_knobs_state_set_bsg(adapter, NULL) ||
	    rt_view_context_knobs_state_set(adapter, NULL) ||
	    rt_view_knob_state_reset_bsg(NULL, RT_VIEW_KNOBS_RATE_BSG) ||
	    rt_view_knobs_hash_bsg(NULL, NULL) ||
	    rt_view_context_knobs_hash_bsg(NULL, NULL) ||
	    rt_view_context_knobs_hash(NULL, NULL) ||
	    rt_view_knobs_translate_bsg(NULL, direct_tvec, 0) ||
	    rt_view_context_knobs_translate_bsg(NULL, direct_tvec, 0) ||
	    rt_view_context_knobs_translate(NULL, direct_tvec, 0) ||
	    rt_view_knobs_translate_bsg(adapter, NULL, 0) ||
	    rt_view_context_knobs_translate_bsg(adapter, NULL, 0) ||
	    rt_view_context_knobs_translate(adapter, NULL, 0) ||
	    rt_view_knobs_rotate_bsg(NULL, direct_rvec, 'v', 'v', NULL, NULL) ||
	    rt_view_knobs_rotate_bsg(adapter, NULL, 'v', 'v', NULL, NULL) ||
	    rt_view_context_knobs_rotate_bsg(NULL, direct_rvec, 'v', 'v', NULL, NULL) ||
	    rt_view_context_knobs_rotate_bsg(adapter, NULL, 'v', 'v', NULL, NULL) ||
	    rt_view_context_knobs_rotate(NULL, direct_rvec, 'v', 'v', NULL, NULL) ||
	    rt_view_context_knobs_rotate(adapter, NULL, 'v', 'v', NULL, NULL) ||
	    rt_view_knobs_update_rate_flags_bsg(NULL) ||
	    rt_view_context_knobs_update_rate_flags_bsg(NULL) ||
	    rt_view_context_knobs_update_rate_flags(NULL) ||
	    rt_view_context_knob_values_from_bsg(NULL, adapter) ||
	    rt_view_context_knob_values_from_bsg(&values, NULL) ||
	    rt_view_context_knob_values_get(NULL, adapter) ||
	    rt_view_context_knob_values_get(&values, NULL) ||
	    rt_view_context_knobs_calibrate_bsg(NULL) ||
	    rt_view_context_knobs_calibrate(NULL)) {
	printf("FAIL: null BSG knob adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    setup_adjust_view(direct);
    setup_adjust_view(adapter);
    direct->k.sca = 7.0;
    direct->k.rot_v[X] = 3.0;
    direct->k.tra_v[Y] = 4.0;
    direct->k.sca_abs = 9.0;
    direct->k.rot_v_abs[Z] = 8.0;
    direct->k.tra_v_abs[X] = 6.0;
    direct->gv_a_scale = 10.0;
    adapter->k = direct->k;
    adapter->gv_a_scale = direct->gv_a_scale;

    if (!rt_view_context_knobs_state_get(&neutral_knobs, adapter) ||
	    !fastf_equal(neutral_knobs.rot_v[X], 3.0) ||
	    !fastf_equal(neutral_knobs.tra_v[Y], 4.0) ||
	    !fastf_equal(neutral_knobs.sca, 7.0) ||
	    !fastf_equal(neutral_knobs.sca_abs, 9.0) ||
	    rt_view_knobs_state_hash(&neutral_knobs, NULL) != bsg_knobs_hash(&direct->k, NULL)) {
	printf("FAIL: neutral RT knob-state read/hash adapter\n");
	ret = 1;
	goto cleanup;
    }

    neutral_knobs.rot_m[Z] = 5.0;
    neutral_knobs.origin_m = 'm';
    neutral_knobs.tra_m_abs[X] = 2.0;
    if (!rt_view_context_knobs_state_set(adapter, &neutral_knobs) ||
	    !fastf_equal(adapter->k.rot_m[Z], 5.0) ||
	    adapter->k.origin_m != 'm' ||
	    !fastf_equal(adapter->k.tra_m_abs[X], 2.0)) {
	printf("FAIL: neutral RT knob-state write adapter\n");
	ret = 1;
	goto cleanup;
    }

    direct->k = adapter->k;
    bsg_knobs_reset(&direct->k, BSG_KNOBS_RATE);
    if (!rt_view_knobs_state_reset(&neutral_knobs, RT_VIEW_KNOBS_RATE) ||
	    !rt_view_context_knobs_state_set(adapter, &neutral_knobs) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: neutral RT knob-state reset adapter\n");
	ret = 1;
	goto cleanup;
    }

    direct->k.sca = 7.0;
    direct->k.rot_v[X] = 3.0;
    direct->k.tra_v[Y] = 4.0;
    direct->k.sca_abs = 9.0;
    direct->k.rot_v_abs[Z] = 8.0;
    direct->k.tra_v_abs[X] = 6.0;
    adapter->k = direct->k;

    if (!rt_view_context_knob_values_get(&values, adapter) ||
	    !fastf_equal(values.rate_rotation[X], 3.0) ||
	    !fastf_equal(values.rate_translation[Y], 4.0) ||
	    !fastf_equal(values.rate_scale, 7.0) ||
	    !fastf_equal(values.absolute_rotation[Z], 8.0) ||
	    !fastf_equal(values.absolute_translation[X], 6.0) ||
	    !fastf_equal(values.absolute_scale, 10.0)) {
	printf("FAIL: BSG knob values context adapter\n");
	ret = 1;
	goto cleanup;
    }

    bsg_knobs_reset(&direct->k, BSG_KNOBS_RATE);
    if (!rt_view_context_knobs_reset(adapter, RT_VIEW_KNOBS_RATE) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob rate reset context adapter\n");
	ret = 1;
	goto cleanup;
    }

    direct->k.sca_abs = 11.0;
    direct->k.tra_v_abs[X] = 2.0;
    adapter->k = direct->k;
    bsg_knobs_reset(&direct->k, BSG_KNOBS_ABS);
    if (!rt_view_knob_state_reset_bsg(&adapter->k, RT_VIEW_KNOBS_ABS_BSG) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob absolute reset adapter\n");
	ret = 1;
	goto cleanup;
    }

    direct->k.rot_m[Y] = 5.0;
    direct->k.origin_m = 'm';
    adapter->k = direct->k;
    if (rt_view_knobs_hash_bsg(adapter, NULL) != bsg_knobs_hash(&direct->k, NULL) ||
	    rt_view_context_knobs_hash_bsg(adapter, NULL) != bsg_knobs_hash(&direct->k, NULL) ||
	    rt_view_context_knobs_hash(adapter, NULL) != bsg_knobs_hash(&direct->k, NULL)) {
	printf("FAIL: BSG knob hash adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (bsg_knobs_cmd_process(&direct_rvec, &direct_do_rot,
		&direct_tvec, &direct_do_tran, direct, "x", 6.0, 'v', 0, 0) !=
	    rt_view_context_knobs_cmd_process(&adapter_rvec, &adapter_do_rot,
		&adapter_tvec, &adapter_do_tran, adapter, "x", 6.0, 'v', 0, 0) ||
	    direct_do_rot != adapter_do_rot ||
	    direct_do_tran != adapter_do_tran ||
	    !VNEAR_EQUAL(direct_rvec, adapter_rvec, SMALL_FASTF) ||
	    !VNEAR_EQUAL(direct_tvec, adapter_tvec, SMALL_FASTF) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob command-process context adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSETALL(direct_rvec, 0.0);
    VSETALL(adapter_rvec, 0.0);
    VSETALL(direct_tvec, 0.0);
    VSETALL(adapter_tvec, 0.0);
    direct_do_rot = 0;
    adapter_do_rot = 0;
    direct_do_tran = 0;
    adapter_do_tran = 0;
    if (bsg_knobs_cmd_process(&direct_rvec, &direct_do_rot,
		&direct_tvec, &direct_do_tran, direct, "x", 6.0, 'v', 0, 0) !=
	    rt_view_context_knobs_cmd_process_bsg(&adapter_rvec, &adapter_do_rot,
		&adapter_tvec, &adapter_do_tran, adapter, "x", 6.0, 'v', 0, 0) ||
	    direct_do_rot != adapter_do_rot ||
	    direct_do_tran != adapter_do_tran ||
	    !VNEAR_EQUAL(direct_rvec, adapter_rvec, SMALL_FASTF) ||
	    !VNEAR_EQUAL(direct_tvec, adapter_tvec, SMALL_FASTF) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob command-process compatibility context adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSET(direct_tvec, 0.4, -0.2, 0.1);
    VMOVE(adapter_tvec, direct_tvec);
    bsg_knobs_tran(direct, direct_tvec, 0);
    if (!rt_view_context_knobs_translate(adapter, adapter_tvec, 0) ||
	    check_view_frame("BSG knob translate adapter", adapter, direct)) {
	ret = 1;
	goto cleanup;
    }

    VSET(direct_rvec, 2.0, 3.0, 4.0);
    VMOVE(adapter_rvec, direct_rvec);
    bsg_knobs_rot(direct, direct_rvec, 'v', 'v', NULL, NULL);
    if (!rt_view_context_knobs_rotate(adapter, adapter_rvec, 'v', 'v', NULL, NULL) ||
	    check_view_frame("BSG knob rotate adapter", adapter, direct)) {
	ret = 1;
	goto cleanup;
    }

    direct->k.sca = 1.0;
    direct->k.tra_m[Z] = 0.5;
    adapter->k = direct->k;
    bsg_update_rate_flags(direct);
    if (!rt_view_context_knobs_update_rate_flags(adapter) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob rate-flag context adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSET(direct->k.tra_v_abs, 1.0, 2.0, 3.0);
    VSET(direct->k.tra_v_abs_last, 4.0, 5.0, 6.0);
    VSET(direct->k.tra_m_abs, 7.0, 8.0, 9.0);
    VSET(direct->k.tra_m_abs_last, 10.0, 11.0, 12.0);
    adapter->k = direct->k;
    VSETALL(direct->k.tra_v_abs, 0.0);
    VSETALL(direct->k.tra_v_abs_last, 0.0);
    VSETALL(direct->k.tra_m_abs, 0.0);
    VSETALL(direct->k.tra_m_abs_last, 0.0);
    if (!rt_view_context_knobs_calibrate(adapter) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob calibrate context adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(direct);
    free_view(adapter);
    return ret ? 1 : 0;
}

static void
setup_snap_view(struct bsg_view *v)
{
    struct bsg_grid_state grid;

    v->gv_width = 512;
    v->gv_height = 512;
    bsg_view_grid_get(v, &grid);
    grid.res_h = 1.0;
    grid.res_v = 1.0;
    bsg_view_grid_set(v, &grid);
}

static int
test_bsg_snap_adapter(void)
{
    struct bsg_view *direct = make_view("rt_view_snap_direct");
    struct bsg_view *adapter = make_view("rt_view_snap_adapter");
    struct bsg_snap_result direct_res = {0};
    struct bsg_snap_result adapter_res = {0};
    struct rt_view_snap_result_bsg *wrapped_res = NULL;
    point_t sample = VINIT_ZERO;
    point_t first_candidate = VINIT_ZERO;
    point_t wrapped_point = VINIT_ZERO;
    point_t bad_point = VINIT_ZERO;
    struct bu_vls wrapped_path = BU_VLS_INIT_ZERO;
    fastf_t direct_x = 0.1;
    fastf_t direct_y = 0.1;
    fastf_t adapter_x = 0.1;
    fastf_t adapter_y = 0.1;
    fastf_t direct_grid_x = 0.1;
    fastf_t direct_grid_y = 0.1;
    fastf_t adapter_grid_x = 0.1;
    fastf_t adapter_grid_y = 0.1;
    rt_view_feature_ref feature = {1234, 5678};
    bsg_feature_ref direct_feature_seed = {feature.token, feature.revision};
    bsg_feature_ref direct_feature = BSG_FEATURE_REF_NULL_INIT;
    bsg_feature_ref adapter_feature = BSG_FEATURE_REF_NULL_INIT;
    struct rt_view_grid_state context_grid = RT_VIEW_GRID_STATE_INIT;
    int direct_cnt = 0;
    int adapter_cnt = 0;
    int ret = 0;

    VSET(sample, 0.1, 0.1, 0.0);

    if ((RT_VIEW_SNAP_SHARED_BSG |
		RT_VIEW_SNAP_LOCAL_BSG |
		RT_VIEW_SNAP_DB_BSG |
		RT_VIEW_SNAP_VIEW_BSG |
		RT_VIEW_SNAP_TCL_BSG) == 0 ||
	    (RT_VIEW_SNAP_KIND_GRID_BSG |
		RT_VIEW_SNAP_KIND_ENDPOINT_BSG |
		RT_VIEW_SNAP_KIND_MIDPOINT_BSG |
		RT_VIEW_SNAP_KIND_INTERSECTION_BSG |
		RT_VIEW_SNAP_KIND_PERPENDICULAR_BSG |
		RT_VIEW_SNAP_KIND_TANGENT_BSG |
		RT_VIEW_SNAP_KIND_OVERLAY_HANDLE_BSG) == 0) {
	printf("FAIL: BSG snap flag alias coverage\n");
	ret = 1;
	goto cleanup;
    }
    if ((RT_VIEW_SNAP_SHARED |
		RT_VIEW_SNAP_LOCAL |
		RT_VIEW_SNAP_DB |
		RT_VIEW_SNAP_VIEW |
		RT_VIEW_SNAP_TCL) == 0 ||
	    (RT_VIEW_SNAP_KIND_GRID |
		RT_VIEW_SNAP_KIND_ENDPOINT |
		RT_VIEW_SNAP_KIND_MIDPOINT |
		RT_VIEW_SNAP_KIND_INTERSECTION |
		RT_VIEW_SNAP_KIND_PERPENDICULAR |
		RT_VIEW_SNAP_KIND_TANGENT |
		RT_VIEW_SNAP_KIND_OVERLAY_HANDLE) == 0) {
	printf("FAIL: neutral snap flag alias coverage\n");
	ret = 1;
	goto cleanup;
    }

    wrapped_res = rt_view_snap_result_create_bsg();
    if (!wrapped_res ||
	    rt_view_snap_result_count_bsg(NULL) != 0 ||
	    rt_view_snap_result_point_bsg(NULL, 0, bad_point) ||
	    rt_view_snap_result_point_bsg(wrapped_res, 0, NULL) ||
	    rt_view_snap_result_distance_bsg(NULL, 0) != 0.0 ||
	    rt_view_snap_result_kind_bsg(NULL, 0) != 0ULL ||
	    rt_view_snap_result_source_path_bsg(NULL, 0, &wrapped_path) ||
	    rt_view_snap_result_source_path_bsg(wrapped_res, 0, NULL) ||
	    rt_view_snap_candidates_result_bsg(adapter, sample, 0.0,
		RT_VIEW_SNAP_KIND_GRID_BSG, NULL)) {
	printf("FAIL: null/empty BSG snap-result adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_snap_candidates_bsg(NULL, sample, 0.0, RT_VIEW_SNAP_KIND_GRID_BSG,
		&adapter_res) ||
	    rt_view_snap_point_2d_bsg(NULL, &adapter_x, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) ||
	    rt_view_snap_point_2d_bsg(adapter, NULL, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) ||
	    rt_view_snap_point_2d_bsg(adapter, &adapter_x, NULL,
		RT_VIEW_SNAP_KIND_GRID_BSG) ||
	    rt_view_snap_grid_2d_bsg(NULL, &adapter_grid_x, &adapter_grid_y) ||
	    rt_view_snap_grid_2d_bsg(adapter, NULL, &adapter_grid_y) ||
	    rt_view_snap_grid_2d_bsg(adapter, &adapter_grid_x, NULL) ||
	    rt_view_context_snap_grid_2d_bsg(NULL, &adapter_grid_x, &adapter_grid_y) ||
	    rt_view_context_snap_grid_2d_bsg(adapter, NULL, &adapter_grid_y) ||
	    rt_view_context_snap_grid_2d_bsg(adapter, &adapter_grid_x, NULL) ||
	    rt_view_context_snap_point_2d(NULL, &adapter_x, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID) ||
	    rt_view_context_snap_point_2d(adapter, NULL, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID) ||
	    rt_view_context_snap_point_2d(adapter, &adapter_x, NULL,
		RT_VIEW_SNAP_KIND_GRID) ||
	    rt_view_context_snap_first_candidate(NULL, sample,
		RT_VIEW_SNAP_KIND_GRID, first_candidate) ||
	    rt_view_context_snap_first_candidate(adapter, NULL,
		RT_VIEW_SNAP_KIND_GRID, first_candidate) ||
	    rt_view_context_snap_first_candidate(adapter, sample, 0,
		first_candidate) ||
	    rt_view_context_snap_first_candidate(adapter, sample,
		RT_VIEW_SNAP_KIND_GRID, NULL) ||
	    rt_view_context_snap_grid_2d(NULL, &adapter_grid_x, &adapter_grid_y) ||
	    rt_view_context_snap_grid_2d(adapter, NULL, &adapter_grid_y) ||
	    rt_view_context_snap_grid_2d(adapter, &adapter_grid_x, NULL)) {
	printf("FAIL: null BSG snap adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    bsg_snap_result_free(&adapter_res);

    if (rt_view_snap_exclude_feature_set_bsg(NULL, feature) ||
	    rt_view_snap_exclude_feature_clear_bsg(NULL) ||
	    rt_view_context_snap_exclude_feature_clear_bsg(NULL)) {
	printf("FAIL: null BSG snap-exclude adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    setup_snap_view(direct);
    setup_snap_view(adapter);

    if (rt_view_context_grid_state_from_bsg(NULL, adapter) ||
	    rt_view_context_grid_state_from_bsg(&context_grid, NULL) ||
	    !rt_view_context_grid_state_from_bsg(&context_grid, adapter) ||
	    !fastf_equal(context_grid.res_h, 1.0) ||
	    !fastf_equal(context_grid.res_v, 1.0)) {
	printf("FAIL: opaque BSG grid-state context adapter\n");
	ret = 1;
	goto cleanup;
    }
    context_grid.res_h = 2.0;
    context_grid.res_v = 2.0;
    if (rt_view_context_grid_state_set_bsg(NULL, &context_grid) ||
	    rt_view_context_grid_state_set_bsg(adapter, NULL) ||
	    !rt_view_context_grid_state_set_bsg(adapter, &context_grid) ||
	    !rt_view_context_grid_state_from_bsg(&context_grid, adapter) ||
	    !fastf_equal(context_grid.res_h, 2.0) ||
	    !fastf_equal(context_grid.res_v, 2.0)) {
	printf("FAIL: opaque BSG grid-state context set adapter\n");
	ret = 1;
	goto cleanup;
    }
    context_grid.res_h = 1.0;
    context_grid.res_v = 1.0;
    if (!rt_view_context_grid_state_set_bsg(adapter, &context_grid)) {
	printf("FAIL: opaque BSG grid-state context restore adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_context_snap_source_flags_set_bsg(NULL,
		RT_VIEW_SNAP_VIEW_BSG) ||
	    rt_view_context_snap_source_flags_from_bsg(NULL) != 0 ||
	    !rt_view_context_snap_source_flags_set_bsg(adapter,
		RT_VIEW_SNAP_VIEW_BSG) ||
	    rt_view_context_snap_source_flags_from_bsg(adapter) !=
		RT_VIEW_SNAP_VIEW_BSG ||
	    !rt_view_context_snap_source_flags_set_bsg(adapter, 0) ||
	    rt_view_context_snap_source_flags_from_bsg(adapter) != 0) {
	printf("FAIL: opaque BSG snap-source context adapter\n");
	ret = 1;
	goto cleanup;
    }

    direct_cnt = bsg_snap_candidates(direct, sample, 0.0,
	    RT_VIEW_SNAP_KIND_GRID_BSG, &direct_res);
    adapter_cnt = rt_view_snap_candidates_bsg(adapter, sample, 0.0,
	    RT_VIEW_SNAP_KIND_GRID_BSG, &adapter_res);
    if (direct_cnt != adapter_cnt ||
	    direct_res.sr_cnt != adapter_res.sr_cnt ||
	    (direct_res.sr_cnt &&
	     !VNEAR_EQUAL(direct_res.sr_candidates[0].sc_point,
		 adapter_res.sr_candidates[0].sc_point, SMALL_FASTF))) {
	printf("FAIL: BSG snap candidates adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_snap_candidates_result_bsg(adapter, sample, 0.0,
		RT_VIEW_SNAP_KIND_GRID_BSG, wrapped_res) != adapter_cnt ||
	    rt_view_snap_result_count_bsg(wrapped_res) != adapter_res.sr_cnt ||
	    (adapter_res.sr_cnt &&
	     (!rt_view_snap_result_point_bsg(wrapped_res, 0, wrapped_point) ||
	      !VNEAR_EQUAL(wrapped_point, adapter_res.sr_candidates[0].sc_point,
		  SMALL_FASTF) ||
	      !fastf_equal(rt_view_snap_result_distance_bsg(wrapped_res, 0),
		  adapter_res.sr_candidates[0].sc_distance) ||
	      rt_view_snap_result_kind_bsg(wrapped_res, 0) !=
		  (unsigned long long)adapter_res.sr_candidates[0].sc_kind ||
	      !rt_view_snap_result_source_path_bsg(wrapped_res, 0,
		  &wrapped_path) ||
	      !BU_STR_EQUAL(bu_vls_cstr(&wrapped_path),
		  bu_vls_cstr(&adapter_res.sr_candidates[0].sc_source_path)))) ||
	    rt_view_snap_result_point_bsg(wrapped_res, adapter_res.sr_cnt,
		wrapped_point) ||
	    rt_view_snap_result_distance_bsg(wrapped_res, adapter_res.sr_cnt) != 0.0 ||
	    rt_view_snap_result_kind_bsg(wrapped_res, adapter_res.sr_cnt) != 0ULL ||
	    rt_view_snap_result_source_path_bsg(wrapped_res, adapter_res.sr_cnt,
		&wrapped_path)) {
	printf("FAIL: BSG opaque snap-result adapter\n");
	ret = 1;
	goto cleanup;
    }

    if ((adapter_cnt > 0) != rt_view_context_snap_first_candidate(adapter,
		sample, RT_VIEW_SNAP_KIND_GRID, first_candidate) ||
	    (adapter_cnt > 0 &&
	     !VNEAR_EQUAL(first_candidate, adapter_res.sr_candidates[0].sc_point,
		 SMALL_FASTF))) {
	printf("FAIL: neutral snap first-candidate context adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (bsg_snap_point_2d(direct, &direct_x, &direct_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) !=
	    rt_view_snap_point_2d_bsg(adapter, &adapter_x, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) ||
	    !fastf_equal(direct_x, adapter_x) ||
	    !fastf_equal(direct_y, adapter_y)) {
	printf("FAIL: BSG 2D snap adapter\n");
	ret = 1;
	goto cleanup;
    }

    direct_x = direct_y = adapter_x = adapter_y = 0.1;
    if (bsg_snap_point_2d(direct, &direct_x, &direct_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) !=
	    rt_view_context_snap_point_2d_bsg(adapter, &adapter_x, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) ||
	    !fastf_equal(direct_x, adapter_x) ||
	    !fastf_equal(direct_y, adapter_y)) {
	printf("FAIL: opaque BSG 2D snap context adapter\n");
	ret = 1;
	goto cleanup;
    }
    direct_x = direct_y = adapter_x = adapter_y = 0.1;
    if (bsg_snap_point_2d(direct, &direct_x, &direct_y,
		RT_VIEW_SNAP_KIND_GRID_BSG) !=
	    rt_view_context_snap_point_2d(adapter, &adapter_x, &adapter_y,
		RT_VIEW_SNAP_KIND_GRID) ||
	    !fastf_equal(direct_x, adapter_x) ||
	    !fastf_equal(direct_y, adapter_y)) {
	printf("FAIL: neutral 2D snap context adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (bsg_snap_grid_2d(direct, &direct_grid_x, &direct_grid_y) !=
	    rt_view_snap_grid_2d_bsg(adapter, &adapter_grid_x, &adapter_grid_y) ||
	    !fastf_equal(direct_grid_x, adapter_grid_x) ||
	    !fastf_equal(direct_grid_y, adapter_grid_y)) {
	printf("FAIL: BSG 2D grid-snap adapter\n");
	ret = 1;
	goto cleanup;
    }
    direct_grid_x = direct_grid_y = adapter_grid_x = adapter_grid_y = 0.1;
    if (bsg_snap_grid_2d(direct, &direct_grid_x, &direct_grid_y) !=
	    rt_view_context_snap_grid_2d_bsg(adapter, &adapter_grid_x, &adapter_grid_y) ||
	    !fastf_equal(direct_grid_x, adapter_grid_x) ||
	    !fastf_equal(direct_grid_y, adapter_grid_y)) {
	printf("FAIL: opaque BSG 2D grid-snap context adapter\n");
	ret = 1;
	goto cleanup;
    }
    direct_grid_x = direct_grid_y = adapter_grid_x = adapter_grid_y = 0.1;
    if (bsg_snap_grid_2d(direct, &direct_grid_x, &direct_grid_y) !=
	    rt_view_context_snap_grid_2d(adapter, &adapter_grid_x, &adapter_grid_y) ||
	    !fastf_equal(direct_grid_x, adapter_grid_x) ||
	    !fastf_equal(direct_grid_y, adapter_grid_y)) {
	printf("FAIL: neutral 2D grid-snap context adapter\n");
	ret = 1;
	goto cleanup;
    }

    bsg_view_snap_exclude_feature_set(direct, direct_feature_seed);
    if (!rt_view_snap_exclude_feature_set_bsg(adapter, feature) ||
	    !bsg_view_snap_exclude_feature_get(direct, &direct_feature) ||
	    !bsg_view_snap_exclude_feature_get(adapter, &adapter_feature) ||
	    direct_feature.token != adapter_feature.token ||
	    direct_feature.revision != adapter_feature.revision) {
	printf("FAIL: BSG snap-exclude set/get adapter\n");
	ret = 1;
	goto cleanup;
    }

    bsg_view_snap_exclude_feature_clear(direct);
    direct_feature = direct_feature_seed;
    adapter_feature = direct_feature_seed;
    if (!rt_view_snap_exclude_feature_clear_bsg(adapter)) {
	printf("FAIL: BSG snap-exclude clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    int direct_present = bsg_view_snap_exclude_feature_get(direct, &direct_feature);
    int adapter_present = bsg_view_snap_exclude_feature_get(adapter, &adapter_feature);
    if (direct_present != adapter_present ||
	    direct_feature.token != adapter_feature.token ||
	    direct_feature.revision != adapter_feature.revision ||
	    direct_feature.token != 0 ||
	    direct_feature.revision != 0) {
	printf("FAIL: BSG snap-exclude clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_view_snap_exclude_feature_set(direct, direct_feature_seed);
    if (!rt_view_snap_exclude_feature_set_bsg(adapter, feature) ||
	    !rt_view_context_snap_exclude_feature_clear_bsg(adapter)) {
	printf("FAIL: BSG context snap-exclude clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_view_snap_exclude_feature_clear(direct);
    direct_feature = direct_feature_seed;
    adapter_feature = direct_feature_seed;
    direct_present = bsg_view_snap_exclude_feature_get(direct, &direct_feature);
    adapter_present = bsg_view_snap_exclude_feature_get(adapter, &adapter_feature);
    if (direct_present != adapter_present ||
	    direct_feature.token != adapter_feature.token ||
	    direct_feature.revision != adapter_feature.revision ||
	    direct_feature.token != 0 ||
	    direct_feature.revision != 0) {
	printf("FAIL: BSG context snap-exclude clear adapter state\n");
	ret = 1;
	goto cleanup;
    }
    bsg_view_snap_exclude_feature_set(direct, direct_feature_seed);
    if (!rt_view_snap_exclude_feature_set_bsg(adapter, feature) ||
	    !rt_view_context_snap_exclude_feature_clear(adapter)) {
	printf("FAIL: neutral context snap-exclude clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_view_snap_exclude_feature_clear(direct);
    direct_feature = direct_feature_seed;
    adapter_feature = direct_feature_seed;
    direct_present = bsg_view_snap_exclude_feature_get(direct, &direct_feature);
    adapter_present = bsg_view_snap_exclude_feature_get(adapter, &adapter_feature);
    if (direct_present != adapter_present ||
	    direct_feature.token != adapter_feature.token ||
	    direct_feature.revision != adapter_feature.revision ||
	    direct_feature.token != 0 ||
	    direct_feature.revision != 0) {
	printf("FAIL: neutral context snap-exclude clear adapter state\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    bu_vls_free(&wrapped_path);
    bsg_snap_result_free(&direct_res);
    bsg_snap_result_free(&adapter_res);
    rt_view_snap_result_free_bsg(wrapped_res);
    free_view(direct);
    free_view(adapter);
    return ret ? 1 : 0;
}

static int
same_rt_measure_as_bsg(const struct rt_view_measure_result *rt_result,
		       const struct bsg_measure_result *bsg_result)
{
    return rt_result &&
	bsg_result &&
	rt_result->valid == bsg_result->mr_valid &&
	fastf_equal(rt_result->distance, bsg_result->mr_distance) &&
	fastf_equal(rt_result->projection, bsg_result->mr_projection) &&
	fastf_equal(rt_result->normal_alignment,
	    bsg_result->mr_normal_alignment);
}

static int
test_bsg_measure_adapter(void)
{
    struct bsg_view *direct = make_view("rt_view_measure_direct");
    struct bsg_view *adapter = make_view("rt_view_measure_adapter");
    struct bsg_measure_result direct_res = {0.0, 0.0, 0.0, 0};
    struct rt_view_measure_result adapter_res = RT_VIEW_MEASURE_RESULT_INIT;
    point_t a = VINIT_ZERO;
    point_t b = VINIT_ZERO;
    int direct_ret = 0;
    int adapter_ret = 0;
    int ret = 0;

    VSET(a, 0.0, 0.0, 0.0);
    VSET(b, 3.0, 4.0, 0.0);

    if (rt_view_measure_candidates_bsg(adapter, a, b, NULL)) {
	printf("FAIL: null BSG measure adapter result argument\n");
	ret = 1;
	goto cleanup;
    }

    direct_ret = bsg_measure_candidates(direct, a, b, &direct_res);
    adapter_ret = rt_view_measure_candidates_bsg(adapter, a, b, &adapter_res);
    if (direct_ret != adapter_ret ||
	    !same_rt_measure_as_bsg(&adapter_res, &direct_res)) {
	printf("FAIL: BSG measure adapter result mismatch\n");
	ret = 1;
	goto cleanup;
    }

    memset(&direct_res, 0, sizeof(direct_res));
    memset(&adapter_res, 0, sizeof(adapter_res));
    direct_ret = bsg_measure_candidates(NULL, a, b, &direct_res);
    adapter_ret = rt_view_measure_candidates_bsg(NULL, a, b, &adapter_res);
    if (direct_ret != adapter_ret ||
	    !same_rt_measure_as_bsg(&adapter_res, &direct_res)) {
	printf("FAIL: null-view BSG measure adapter result mismatch\n");
	ret = 1;
	goto cleanup;
    }

    VMOVE(b, a);
    adapter_res.distance = 7.0;
    adapter_res.projection = 8.0;
    adapter_res.normal_alignment = 9.0;
    adapter_res.valid = 1;
    adapter_ret = rt_view_measure_candidates_bsg(adapter, a, b, &adapter_res);
    if (adapter_ret ||
	    adapter_res.valid ||
	    !fastf_equal(adapter_res.distance, 0.0) ||
	    !fastf_equal(adapter_res.projection, 0.0) ||
	    !fastf_equal(adapter_res.normal_alignment, 0.0)) {
	printf("FAIL: zero-length BSG measure adapter result\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(direct);
    free_view(adapter);
    return ret ? 1 : 0;
}

static int
test_bsg_feature_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_feature_adapter");
    const int owner_token = 42;
    rt_view_feature_ref null_ref = RT_VIEW_FEATURE_REF_NULL;
    rt_view_feature_ref label_ref = RT_VIEW_FEATURE_REF_NULL;
    rt_view_feature_ref preview_ref = RT_VIEW_FEATURE_REF_NULL;
    bsg_feature_ref bsg_label = BSG_FEATURE_REF_NULL_INIT;
    bsg_feature_ref bsg_preview = BSG_FEATURE_REF_NULL_INIT;
    struct rt_view_feature_summary_bsg summary =
	RT_VIEW_FEATURE_SUMMARY_BSG_INIT;
    struct rt_view_feature_geometry_summary_bsg geom_summary =
	RT_VIEW_FEATURE_GEOMETRY_SUMMARY_BSG_INIT;
    struct rt_view_feature_geometry_summary neutral_geom_summary =
	RT_VIEW_FEATURE_GEOMETRY_SUMMARY_INIT;
    struct bu_vls label_text = BU_VLS_INIT_ZERO;
    point_t label_point = VINIT_ZERO;
    unsigned char label_color[3] = {0, 0, 0};
    point_t preview_points[2] = {VINIT_ZERO, VINIT_ZERO};
    int preview_cmds[2] = {BSG_GEOMETRY_LINE_MOVE, BSG_GEOMETRY_LINE_DRAW};
    point_t *copied_points = NULL;
    int *copied_cmds = NULL;
    size_t copied_count = 0;
    struct test_feature_preview_callback_state callback_state = {37, 0, 0, 0, 29};
    struct rt_view_edit_preview_callbacks callbacks =
	RT_VIEW_EDIT_PREVIEW_CALLBACKS_INIT;
    int pick_value = 0;
    int ret = 0;

    if (!v) {
	printf("FAIL: BSG feature adapter view allocation\n");
	return 1;
    }

    if (!rt_view_feature_ref_is_null(null_ref) ||
	    !rt_view_feature_ref_is_null_bsg(null_ref) ||
	    !rt_view_feature_ref_is_null_bsg(
		rt_view_feature_overlay_ensure_bsg(NULL, "bad", NULL, NULL,
		    NULL, NULL)) ||
	    !rt_view_feature_ref_is_null(
		rt_view_context_feature_overlay_ensure(NULL, "bad", NULL,
		    NULL, NULL, NULL)) ||
	    !rt_view_feature_ref_is_null_bsg(
		rt_view_context_feature_overlay_ensure_bsg(NULL, "bad", NULL,
		    NULL, NULL, NULL)) ||
	    !rt_view_feature_ref_is_null_bsg(
		rt_view_feature_label_ensure_bsg(NULL, "bad", NULL)) ||
	    !rt_view_feature_ref_is_null(
		rt_view_context_feature_label_ensure(NULL, "bad", NULL)) ||
	    !rt_view_feature_ref_is_null_bsg(
		rt_view_context_feature_label_ensure_bsg(NULL, "bad", NULL)) ||
	    rt_view_context_feature_remove(NULL, "bad") ||
	    rt_view_feature_remove_bsg(NULL, "bad") ||
	    rt_view_context_feature_remove_bsg(NULL, "bad") ||
	    rt_view_feature_summary_bsg(NULL, "bad", &summary) ||
	    rt_view_feature_summary_bsg(v, NULL, &summary) ||
	    rt_view_feature_summary_bsg(v, "bad", NULL) ||
	    rt_view_context_feature_summary_bsg(NULL, "bad", &summary) ||
	    rt_view_context_feature_summary_bsg(v, NULL, &summary) ||
	    rt_view_context_feature_summary_bsg(v, "bad", NULL) ||
	    rt_view_context_feature_geometry_summary_bsg(NULL, "bad",
		&geom_summary) ||
	    rt_view_context_feature_geometry_summary_bsg(v, NULL,
		&geom_summary) ||
	    rt_view_context_feature_geometry_summary_bsg(v, "bad", NULL) ||
	    rt_view_context_feature_geometry_summary(NULL, "bad",
		&neutral_geom_summary) ||
	    rt_view_context_feature_geometry_summary(v, NULL,
		&neutral_geom_summary) ||
	    rt_view_context_feature_geometry_summary(v, "bad", NULL) ||
	    rt_view_feature_labels_replace(null_ref, NULL, 0) ||
	    rt_view_feature_labels_replace_bsg(null_ref, NULL, 0) ||
	    rt_view_feature_points_replace(null_ref,
		RT_VIEW_FEATURE_UNKNOWN, NULL, NULL, 0) ||
	    rt_view_feature_points_replace_bsg(null_ref,
		RT_VIEW_FEATURE_UNKNOWN, NULL, NULL, 0) ||
	    rt_view_feature_clear_geometry(null_ref) ||
	    rt_view_feature_clear_geometry_bsg(null_ref) ||
	    rt_view_feature_touch(null_ref) ||
	    rt_view_feature_touch_bsg(null_ref)) {
	printf("FAIL: BSG feature adapter null argument handling\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_feature_geometry_summary_bsg(v, "missing_feature",
	    &geom_summary) ||
	    geom_summary.exists ||
	    geom_summary.point_count != 0 ||
	    geom_summary.command_count != 0) {
	printf("FAIL: BSG missing feature geometry adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_feature_geometry_summary(v, "missing_feature",
	    &neutral_geom_summary) ||
	    neutral_geom_summary.exists ||
	    neutral_geom_summary.point_count != 0 ||
	    neutral_geom_summary.command_count != 0) {
	printf("FAIL: neutral missing feature geometry adapter\n");
	ret = 1;
	goto cleanup;
    }

    label_ref = rt_view_context_feature_label_ensure(v, "feature_label",
	    &owner_token);
    bsg_label = test_feature_ref_to_bsg(label_ref);
    struct bsg_feature_record label_record;
    if (rt_view_feature_ref_is_null(label_ref) ||
	    !bsg_feature_record_get(bsg_label, &label_record) ||
	    label_record.family != BSG_FEATURE_LABEL ||
	    !rt_view_feature_summary_bsg(v, "feature_label", &summary) ||
	    !rt_view_context_feature_summary_bsg(v, "feature_label", &summary) ||
	    !summary.exists ||
	    !summary.is_label ||
	    summary.is_overlay) {
	printf("FAIL: BSG label feature ensure adapter\n");
	ret = 1;
	goto cleanup;
    }

    struct rt_view_feature_label labels[1];
    memset(labels, 0, sizeof(labels));
    labels[0].text = "adapter label";
    VSET(labels[0].point, 1.0, 2.0, 3.0);
    labels[0].color_valid = 1;
    labels[0].color[0] = 10;
    labels[0].color[1] = 20;
    labels[0].color[2] = 30;
    if (!rt_view_feature_labels_replace(label_ref, labels, 1) ||
	    bsg_feature_label_count(bsg_label) != 1 ||
	    !bsg_feature_label_copy(bsg_label, 0, &label_text, label_point,
		label_color) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&label_text), labels[0].text) ||
	    !VNEAR_EQUAL(label_point, labels[0].point, SMALL_FASTF) ||
	    label_color[0] != labels[0].color[0] ||
	    label_color[1] != labels[0].color[1] ||
	    label_color[2] != labels[0].color[2]) {
	printf("FAIL: BSG label feature replace adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_feature_set_context(label_ref, v);
    rt_view_feature_set_visible(label_ref, 0);
    rt_view_feature_set_color(label_ref, 4, 5, 6);
    if (!bsg_feature_record_get(bsg_label, &label_record) ||
	    label_record.visible != 0 ||
	    label_record.color[0] != 4 ||
	    label_record.color[1] != 5 ||
	    label_record.color[2] != 6 ||
	    !rt_view_feature_summary_bsg(v, "feature_label", &summary) ||
	    summary.visible != 0 ||
	    summary.color[0] != 4 ||
	    summary.color[1] != 5 ||
	    summary.color[2] != 6) {
	printf("FAIL: BSG feature style adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_feature_set_visible_bsg(label_ref, 1);
    rt_view_feature_set_color_bsg(label_ref, 7, 8, 9);
    if (!bsg_feature_record_get(bsg_label, &label_record) ||
	    label_record.visible != 1 ||
	    label_record.color[0] != 7 ||
	    label_record.color[1] != 8 ||
	    label_record.color[2] != 9) {
	printf("FAIL: BSG feature style compatibility adapter\n");
	ret = 1;
	goto cleanup;
    }

    callbacks.revision_cb = test_feature_preview_revision;
    callbacks.pick_cb = test_feature_preview_pick;
    preview_ref = rt_view_context_feature_overlay_ensure(v,
	    "feature_preview", &owner_token, &callback_state, &callbacks,
	    "source.s");
    bsg_preview = test_feature_ref_to_bsg(preview_ref);
    struct bsg_feature_record preview_record;
    if (rt_view_feature_ref_is_null(preview_ref) ||
	    !bsg_feature_record_get(bsg_preview, &preview_record) ||
	    preview_record.family != BSG_FEATURE_OVERLAY ||
	    !rt_view_feature_summary_bsg(v, "feature_preview", &summary) ||
	    !summary.exists ||
	    !summary.is_overlay ||
	    summary.is_label ||
	    bsg_feature_edit_preview_revision(bsg_preview) !=
		callback_state.revision ||
	    !bsg_feature_edit_preview_pick(bsg_preview, v, 7, 11,
		&pick_value) ||
	    callback_state.pick_calls != 1 ||
	    callback_state.last_x != 7 ||
	    callback_state.last_y != 11 ||
	    pick_value != callback_state.pick_value ||
	    !rt_view_edit_preview_publish_event_bsg(v, preview_ref,
		RT_VIEW_EDIT_PREVIEW_UPDATE, "source.s") ||
	    !rt_view_context_edit_preview_publish_event(v, preview_ref,
		RT_VIEW_EDIT_PREVIEW_COMMIT, "source.s")) {
	printf("FAIL: BSG edit-preview feature adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSET(preview_points[1], 9.0, 8.0, 7.0);
    if (!rt_view_feature_points_replace(preview_ref,
	    RT_VIEW_FEATURE_TRANSIENT_PREVIEW, preview_points,
	    preview_cmds, 2) ||
	    !bsg_feature_record_get(bsg_preview, &preview_record) ||
	    preview_record.family != BSG_FEATURE_TRANSIENT_PREVIEW ||
	    !rt_view_feature_summary_bsg(v, "feature_preview", &summary) ||
	    !summary.exists ||
	    !summary.is_transient_preview ||
	    summary.geometry_command_count != 2 ||
	    !bsg_feature_points_copy(bsg_preview, &copied_points,
		&copied_cmds, &copied_count) ||
	    copied_count != 2 ||
	    !VNEAR_EQUAL(copied_points[1], preview_points[1], SMALL_FASTF) ||
	    copied_cmds[0] != preview_cmds[0] ||
	    copied_cmds[1] != preview_cmds[1] ||
	    !rt_view_context_feature_geometry_summary_bsg(v,
		"feature_preview", &geom_summary) ||
	    !geom_summary.exists ||
	    geom_summary.point_count != 2 ||
	    geom_summary.command_count != 2 ||
	    !rt_view_context_feature_geometry_summary(v,
		"feature_preview", &neutral_geom_summary) ||
	    !neutral_geom_summary.exists ||
	    neutral_geom_summary.point_count != 2 ||
	    neutral_geom_summary.command_count != 2 ||
	    !rt_view_feature_clear_geometry(preview_ref)) {
	printf("FAIL: BSG feature geometry adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_feature_remove(v, "feature_label") ||
	    !rt_view_context_feature_remove(v, "feature_preview") ||
	    !rt_view_feature_summary_bsg(v, "feature_label", &summary) ||
	    summary.exists ||
	    !rt_view_feature_summary_bsg(v, "feature_preview", &summary) ||
	    summary.exists ||
	    !bsg_feature_ref_is_null(bsg_feature_find(v, "feature_label")) ||
	    !bsg_feature_ref_is_null(bsg_feature_find(v, "feature_preview"))) {
	printf("FAIL: BSG feature remove adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    bu_vls_free(&label_text);
    if (copied_points)
	bu_free(copied_points, "RT view feature adapter copied points");
    if (copied_cmds)
	bu_free(copied_cmds, "RT view feature adapter copied commands");
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_view_scope_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_scope_adapter");
    struct bsg_view *src = make_view("rt_view_init_copy_source");
    struct bsg_view *direct = NULL;
    struct bsg_view *adapter = NULL;
    struct bsg_view_set set = {0};
    struct bsg_view *set_view = NULL;
    struct bsg_view *lifecycle_view = NULL;
    struct bsg_view *free_contents_view = NULL;
    void *release_storage_view = NULL;
    void *neutral_release_storage_view = NULL;
    void *opaque_view = NULL;
    void *neutral_set = NULL;
    struct bu_ptbl *callbacks = NULL;
    struct bu_ptbl *callbacks_replacement = NULL;
    struct bu_ptbl *callbacks_neutral = NULL;
    struct rt_view_pick_result_bsg *pick = NULL;
    struct rt_view_pick_result_bsg *empty_pick = NULL;
    struct rt_view_pick_result_bsg *rt_pick = NULL;
    struct rt_view_pick_result_bsg *rt_first_pick = NULL;
    void *rt_pick_ctx = NULL;
    void *rt_first_pick_ctx = NULL;
    void *neutral_pick_ctx = NULL;
    void *neutral_first_pick_ctx = NULL;
    void *line_set_ctx = NULL;
    void *neutral_line_set_ctx = NULL;
    struct bu_vls rt_pick_path = BU_VLS_INIT_ZERO;
    struct test_selection_path_callback_state callback_state = {0, ""};
    struct rt_view_render_summary_bsg render_summary = RT_VIEW_RENDER_SUMMARY_BSG_INIT;
    struct rt_view_render_summary neutral_render_summary =
	RT_VIEW_RENDER_SUMMARY_INIT;
    struct rt_view_render_export_consistency_bsg consistency =
	RT_VIEW_RENDER_EXPORT_CONSISTENCY_BSG_INIT;
    struct rt_view_render_export_consistency neutral_consistency =
	RT_VIEW_RENDER_EXPORT_CONSISTENCY_INIT;
    rt_view_scene_ref root_neutral = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref direct_root_neutral = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref scope_neutral = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref direct_scope_neutral = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref_bsg scope = RT_VIEW_SCENE_REF_BSG_NULL_INIT;
    rt_view_scene_ref_bsg direct_scope = RT_VIEW_SCENE_REF_BSG_NULL_INIT;
    int set_initialized = 0;
    int ret = 0;

    rt_view_init_bsg(NULL, NULL);
    rt_view_free_bsg(NULL);
    rt_view_set_init_bsg(NULL);
    rt_view_set_context_init(NULL);
    rt_view_set_free_bsg(NULL);
    rt_view_set_context_free(NULL);
    rt_view_set_add_view_bsg(NULL, NULL);
    rt_view_set_remove_view_bsg(NULL, NULL);
    if (rt_view_context_free_contents_bsg(NULL) ||
	    rt_view_context_view_set_attach_bsg(NULL, &set) ||
	    rt_view_context_view_set_attach(NULL, &set)) {
	printf("FAIL: null opaque BSG lifecycle adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_set_views_bsg(NULL)) {
	printf("FAIL: null BSG view-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_context_views(NULL)) {
	printf("FAIL: null neutral view-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_recycle_pool_bsg(NULL) != bsg_set_fsos(NULL)) {
	printf("FAIL: null BSG view-set recycle-pool adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_context_recycle_pool(NULL) != bsg_set_fsos(NULL)) {
	printf("FAIL: null neutral view-set recycle-pool adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_set_context_init(&set);
    set_initialized = 1;
    if (rt_view_set_views_bsg(&set) != bsg_set_views(&set) ||
	    BU_PTBL_LEN(rt_view_set_views_bsg(&set)) != 0) {
	printf("FAIL: empty BSG view-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_context_views(&set) != bsg_set_views(&set) ||
	    BU_PTBL_LEN(rt_view_set_context_views(&set)) != 0) {
	printf("FAIL: empty neutral view-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_recycle_pool_bsg(&set) != bsg_set_fsos(&set)) {
	printf("FAIL: BSG view-set recycle-pool adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_context_recycle_pool(&set) != bsg_set_fsos(&set)) {
	printf("FAIL: neutral view-set recycle-pool adapter\n");
	ret = 1;
	goto cleanup;
    }

    BU_ALLOC(set_view, struct bsg_view);
    bsg_init(set_view, NULL);
    bsg_set_add_view(&set, set_view);
    if (rt_view_set_views_bsg(&set) != bsg_set_views(&set) ||
	    BU_PTBL_LEN(rt_view_set_views_bsg(&set)) != 1 ||
	    (struct bsg_view *)BU_PTBL_GET(rt_view_set_views_bsg(&set), 0) != set_view) {
	printf("FAIL: populated BSG view-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_context_views(&set) != bsg_set_views(&set) ||
	    BU_PTBL_LEN(rt_view_set_context_views(&set)) != 1 ||
	    (struct bsg_view *)BU_PTBL_GET(rt_view_set_context_views(&set), 0) != set_view) {
	printf("FAIL: populated neutral view-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bu_vls_sprintf(&set_view->gv_name, "rt_view_set_find_adapter");
    if (rt_view_set_find_view_bsg(NULL, "rt_view_set_find_adapter") ||
	    rt_view_set_find_view_bsg(&set, NULL) ||
	    rt_view_set_find_view_bsg(&set, "missing") ||
	    rt_view_set_find_view_bsg(&set, "rt_view_set_find_adapter") !=
	    bsg_set_find_view(&set, "rt_view_set_find_adapter")) {
	printf("FAIL: BSG view-set find adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_set_context_find_view(NULL, "rt_view_set_find_adapter") ||
	    rt_view_set_context_find_view(&set, NULL) ||
	    rt_view_set_context_find_view(&set, "missing") ||
	    rt_view_set_context_find_view(&set, "rt_view_set_find_adapter") !=
	    bsg_set_find_view(&set, "rt_view_set_find_adapter")) {
	printf("FAIL: neutral view-set find adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_context_name_set_bsg(NULL, "missing") ||
	    rt_view_context_name_set(NULL, "missing") ||
	    rt_view_context_name_from_bsg(NULL) ||
	    rt_view_context_name_get(NULL) ||
	    rt_view_context_init_bsg(NULL, &set) ||
	    rt_view_context_release_storage_bsg(NULL) ||
	    rt_view_context_release_storage(NULL) ||
	    rt_view_context_callbacks_set_bsg(NULL, NULL) ||
	    rt_view_context_callbacks_set(NULL, NULL) ||
	    rt_view_set_context_add_bsg(NULL, NULL) ||
	    rt_view_set_context_add(NULL, NULL) ||
	    rt_view_set_context_remove_bsg(NULL, NULL) ||
	    rt_view_set_context_remove(NULL, NULL) ||
	    rt_view_context_scale_state_set_bsg(NULL, 1.0, 1.0, 0.0,
		2.0, 0.5) ||
	    rt_view_context_scale_state_set(NULL, 1.0, 1.0, 0.0,
		2.0, 0.5) ||
	    rt_view_context_lod_bounds_callback_set_bsg(NULL) ||
	    rt_view_context_lod_bounds_callback_set(NULL) ||
	    rt_view_context_lod_bounds_callback_is_bsg(NULL) ||
	    rt_view_context_lod_bounds_callback_is(NULL) ||
	    rt_view_context_is_independent_bsg(NULL) ||
	    rt_view_context_is_independent(NULL) ||
	    rt_view_context_is_valid(NULL) ||
	    !rt_view_context_independent_scope_is_null_bsg(NULL, 1) ||
	    !rt_view_context_independent_scope_is_null(NULL, 1) ||
	    !fastf_equal(rt_view_context_size_from_bsg(NULL), 0.0) ||
	    !fastf_equal(rt_view_context_size_get(NULL), 0.0)) {
	printf("FAIL: null opaque BSG view-context setup adapters\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_context_independent_scope_destroy(NULL);
    rt_view_context_free(NULL);

    neutral_set = rt_view_set_context_create();
    opaque_view = rt_view_context_create_with_set(neutral_set);
    if (!neutral_set ||
	    !opaque_view ||
	    !rt_view_context_is_valid(opaque_view) ||
	    rt_view_context_is_bsg(opaque_view) ||
	    rt_view_context_is_retained(opaque_view) ||
	    !rt_view_context_name_set(opaque_view,
		"rt_view_native_lifecycle_adapter") ||
	    !rt_view_set_context_add(neutral_set, opaque_view) ||
	    rt_view_set_context_find_view(neutral_set,
		"rt_view_native_lifecycle_adapter") != opaque_view ||
	    !rt_view_context_scale_state_set(opaque_view, 10.0, 20.0,
		30.0, 40.0, 0.025) ||
	    !fastf_equal(rt_view_context_scale_get(opaque_view), 10.0) ||
	    !fastf_equal(rt_view_context_size_get(opaque_view), 40.0) ||
	    rt_view_context_scene_attached(opaque_view) ||
	    rt_view_context_scene_shared(set_view, opaque_view) ||
	    rt_view_context_is_independent(opaque_view) ||
	    !rt_view_context_independent_scope_is_null(opaque_view, 0) ||
	    !rt_view_context_lod_bounds_callback_set(opaque_view) ||
	    !rt_view_context_lod_bounds_callback_is(opaque_view)) {
	printf("FAIL: native neutral view-context lifecycle\n");
	ret = 1;
	goto cleanup;
    }
    BU_GET(callbacks_neutral, struct bu_ptbl);
    bu_ptbl_init(callbacks_neutral, 8,
	    "rt view native callback table neutral test");
    if (!rt_view_context_callbacks_set(opaque_view, callbacks_neutral)) {
	printf("FAIL: native neutral view-context callback table set\n");
	ret = 1;
	goto cleanup;
    }
    callbacks_neutral = NULL;
    if (!rt_view_context_callbacks_set(opaque_view, NULL)) {
	printf("FAIL: native neutral view-context callback table clear\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_context_free(opaque_view);
    opaque_view = NULL;

    opaque_view = rt_view_context_create_copy_with_set(v, neutral_set);
    if (!opaque_view ||
	    !rt_view_context_is_valid(opaque_view) ||
	    rt_view_context_is_bsg(opaque_view) ||
	    !BU_STR_EQUAL(rt_view_context_name_get(opaque_view),
		bu_vls_cstr(&v->gv_name)) ||
	    !fastf_equal(rt_view_context_scale_get(opaque_view),
		v->gv_scale) ||
	    !fastf_equal(rt_view_context_size_get(opaque_view),
		v->gv_size)) {
	printf("FAIL: native neutral view-context copy from BSG source\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_context_free(opaque_view);
    opaque_view = NULL;
    rt_view_set_context_destroy(neutral_set);
    neutral_set = NULL;

    opaque_view = rt_view_context_create_with_set_bsg(&set);
    if (!opaque_view ||
	    !rt_view_context_is_valid(opaque_view) ||
	    !rt_view_context_is_bsg(opaque_view) ||
	    !rt_view_context_scene_attached_bsg(opaque_view) ||
	    !rt_view_context_scene_attached(opaque_view) ||
	    !rt_view_context_scene_shared_bsg(set_view, opaque_view) ||
	    !rt_view_context_scene_shared(set_view, opaque_view)) {
	printf("FAIL: opaque BSG view-context create-with-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_context_free(opaque_view);
    opaque_view = NULL;

    opaque_view = rt_view_context_create_bsg();
    if (!opaque_view ||
	    !rt_view_context_is_valid(opaque_view) ||
	    !rt_view_context_is_bsg(opaque_view) ||
	    !rt_view_context_name_set_bsg(opaque_view,
		"rt_view_opaque_lifecycle_adapter") ||
	    !BU_STR_EQUAL(rt_view_context_name_from_bsg(opaque_view),
		"rt_view_opaque_lifecycle_adapter") ||
	    !rt_view_context_name_set(opaque_view,
		"rt_view_opaque_lifecycle_adapter_neutral") ||
	    !BU_STR_EQUAL(rt_view_context_name_get(opaque_view),
		"rt_view_opaque_lifecycle_adapter_neutral") ||
	    !rt_view_set_context_add(&set, opaque_view) ||
	    rt_view_set_find_view_bsg(&set,
		"rt_view_opaque_lifecycle_adapter_neutral") != opaque_view ||
	    !rt_view_context_scale_state_set(opaque_view, 10.0, 20.0,
		30.0, 40.0, 0.025) ||
	    !fastf_equal(rt_view_context_scale_get(opaque_view), 10.0) ||
	    !fastf_equal(rt_view_context_size_get(opaque_view), 40.0) ||
	    rt_view_context_is_independent_bsg(opaque_view) ||
	    rt_view_context_is_independent(opaque_view) ||
	    !rt_view_context_independent_scope_is_null_bsg(opaque_view, 0) ||
	    !rt_view_context_independent_scope_is_null(opaque_view, 0)) {
	printf("FAIL: opaque BSG view-context setup adapters\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_lod_bounds_callback_set(opaque_view)) {
	printf("FAIL: opaque BSG view-context LoD callback adapter set\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_lod_bounds_callback_is_bsg(opaque_view) ||
	    !rt_view_context_lod_bounds_callback_is(opaque_view)) {
	printf("FAIL: opaque BSG view-context LoD callback adapter\n");
	ret = 1;
	goto cleanup;
    }
    BU_GET(callbacks, struct bu_ptbl);
    bu_ptbl_init(callbacks, 8, "rt view legacy callback adapter test");
    BU_GET(callbacks_replacement, struct bu_ptbl);
    bu_ptbl_init(callbacks_replacement, 8,
	    "rt view legacy callback adapter replacement test");
    if (!rt_view_context_callbacks_set_bsg(opaque_view, callbacks) ||
	    ((struct bsg_view *)opaque_view)->callbacks != callbacks) {
	printf("FAIL: opaque BSG view-context callback table adapter set\n");
	ret = 1;
	goto cleanup;
    }
    callbacks = NULL;
    if (!rt_view_context_callbacks_set_bsg(opaque_view, callbacks_replacement) ||
	    ((struct bsg_view *)opaque_view)->callbacks != callbacks_replacement) {
	printf("FAIL: opaque BSG view-context callback table adapter replace\n");
	ret = 1;
	goto cleanup;
    }
    callbacks_replacement = NULL;
    if (!rt_view_context_callbacks_set_bsg(opaque_view, NULL) ||
	    ((struct bsg_view *)opaque_view)->callbacks) {
	printf("FAIL: opaque BSG view-context callback table adapter\n");
	ret = 1;
	goto cleanup;
    }
    BU_GET(callbacks_neutral, struct bu_ptbl);
    bu_ptbl_init(callbacks_neutral, 8,
	    "rt view retained callback adapter neutral test");
    if (!rt_view_context_callbacks_set(opaque_view, callbacks_neutral) ||
	    ((struct bsg_view *)opaque_view)->callbacks != callbacks_neutral) {
	printf("FAIL: opaque retained view-context callback table adapter set\n");
	ret = 1;
	goto cleanup;
    }
    callbacks_neutral = NULL;
    if (!rt_view_context_callbacks_set(opaque_view, NULL) ||
	    ((struct bsg_view *)opaque_view)->callbacks) {
	printf("FAIL: opaque retained view-context callback table adapter clear\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_set_context_remove(&set, opaque_view) ||
	    rt_view_set_context_find_view(&set,
		"rt_view_opaque_lifecycle_adapter")) {
	printf("FAIL: opaque neutral view-context remove adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_context_free(opaque_view);
    opaque_view = NULL;

    opaque_view = rt_view_context_create_copy_with_set_bsg(v, &set);
    if (!opaque_view ||
	    !rt_view_context_is_valid(opaque_view) ||
	    !rt_view_context_is_bsg(opaque_view) ||
	    !rt_view_context_scene_attached(opaque_view) ||
	    !rt_view_context_scene_shared(set_view, opaque_view) ||
	    check_view_frame("neutral create-copy-with-set adapter",
		(struct bsg_view *)opaque_view, v)) {
	printf("FAIL: retained neutral view-context create-copy-with-set adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_context_free(opaque_view);
    opaque_view = NULL;

    BU_ALLOC(lifecycle_view, struct bsg_view);
    if (!rt_view_context_init_bsg(lifecycle_view, NULL)) {
	printf("FAIL: opaque BSG view init adapter\n");
	ret = 1;
	goto cleanup;
    }
    bu_vls_sprintf(&lifecycle_view->gv_name, "rt_view_lifecycle_adapter");
    if (lifecycle_view->magic != BSG_VIEW_MAGIC || lifecycle_view->vset) {
	printf("FAIL: BSG view init adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_view_set_attach_bsg(lifecycle_view, &set) ||
	    lifecycle_view->vset != &set ||
	    !rt_view_context_view_set_attach_bsg(lifecycle_view, NULL) ||
	    lifecycle_view->vset) {
	printf("FAIL: BSG view-set attach context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_view_set_attach(lifecycle_view, &set) ||
	    lifecycle_view->vset != &set ||
	    !rt_view_context_view_set_attach(lifecycle_view, NULL) ||
	    lifecycle_view->vset) {
	printf("FAIL: retained view-set attach context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_is_bsg(NULL) ||
	    rt_view_context_is_valid(NULL) ||
	    rt_view_context_is_retained(NULL) ||
	    !rt_view_context_is_valid(lifecycle_view) ||
	    !rt_view_context_is_bsg(lifecycle_view) ||
	    !rt_view_context_is_retained(lifecycle_view)) {
	printf("FAIL: BSG view context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_set_context_add_bsg(&set, lifecycle_view)) {
	printf("FAIL: opaque BSG view-set add adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (lifecycle_view->vset != &set ||
	    rt_view_set_context_find_view(&set, "rt_view_lifecycle_adapter") != lifecycle_view) {
	printf("FAIL: BSG view-set add adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_set_context_remove_bsg(&set, lifecycle_view);
    if (lifecycle_view->vset ||
	    rt_view_set_context_find_view(&set, "rt_view_lifecycle_adapter")) {
	printf("FAIL: BSG view-set remove adapter\n");
	ret = 1;
	goto cleanup;
    }

    BU_ALLOC(free_contents_view, struct bsg_view);
    if (!rt_view_context_init_bsg(free_contents_view, NULL) ||
	    !rt_view_context_name_set_bsg(free_contents_view,
		"rt_view_free_contents_adapter") ||
	    !rt_view_context_free_contents_bsg(free_contents_view)) {
	printf("FAIL: BSG view contents free context adapter\n");
	ret = 1;
	goto cleanup;
    }
    BU_PUT(free_contents_view, struct bsg_view);
    free_contents_view = NULL;

    release_storage_view = rt_view_context_create_bsg();
    if (!release_storage_view ||
	    !rt_view_context_name_set_bsg(release_storage_view,
		"rt_view_release_storage_adapter") ||
	    !rt_view_context_free_contents_bsg(release_storage_view) ||
	    !rt_view_context_release_storage_bsg(release_storage_view)) {
	printf("FAIL: BSG view storage release context adapter\n");
	ret = 1;
	goto cleanup;
    }
    release_storage_view = NULL;

    neutral_release_storage_view = rt_view_context_create();
    if (!neutral_release_storage_view ||
	    !rt_view_context_name_set(neutral_release_storage_view,
		"rt_view_release_storage_adapter_neutral") ||
	    !rt_view_context_release_storage(neutral_release_storage_view)) {
	printf("FAIL: retained view storage release context adapter\n");
	ret = 1;
	goto cleanup;
    }
    neutral_release_storage_view = NULL;

    if (rt_view_selection_available_bsg(NULL) ||
	    rt_view_context_selection_available(NULL) ||
	    rt_view_selection_count_bsg(NULL) != 0 ||
	    rt_view_context_selection_count(NULL) != 0) {
	printf("FAIL: null BSG selection count adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_selection_available_bsg(v) ||
	    !rt_view_context_selection_available(v) ||
	    rt_view_selection_count_bsg(v) != 0 ||
	    rt_view_context_selection_count(v) != 0) {
	printf("FAIL: BSG selection count adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_selection_available_bsg(v) ||
	    rt_view_context_selection_count_bsg(v) != 0 ||
	    !rt_view_context_selection_available(v) ||
	    rt_view_context_selection_count(v) != 0) {
	printf("FAIL: opaque BSG selection count context adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_selection_set_pick_result_ref_bsg(NULL, NULL, NULL, NULL)) {
	printf("FAIL: null opaque selection pick adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_selection_set_pick_result_ref_bsg(NULL, NULL, NULL,
		NULL)) {
	printf("FAIL: null opaque selection pick context adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_selection_set_pick_result_ref_bsg(v, NULL, NULL, NULL) ||
	    rt_view_selection_count_bsg(v) != 0 ||
	    rt_view_context_selection_count(v) != 0) {
	printf("FAIL: opaque selection null-pick clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_selection_set_pick_result_ref_bsg(v, NULL, NULL,
		NULL) ||
	    rt_view_context_selection_count_bsg(v) != 0 ||
	    rt_view_context_selection_count(v) != 0) {
	printf("FAIL: opaque selection null-pick context clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    pick = rt_view_pick_result_create_bsg();
    if (!pick ||
	    !rt_view_pick_result_append_path_bsg(pick, v, 10, 20,
		"/adapter/path", 1.0)) {
	printf("FAIL: opaque selection pick test record setup\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_selection_set_pick_result_ref_bsg(NULL, pick,
		test_selection_path_callback, &callback_state) ||
	    callback_state.count != 1 ||
	    !BU_STR_EQUAL(callback_state.last_path, "/adapter/path")) {
	printf("FAIL: opaque selection pick adapter path callback\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_selection_set_pick_result_ref_bsg(v, pick,
		test_selection_path_callback, &callback_state) ||
	    rt_view_selection_count_bsg(v) != 1 ||
	    rt_view_context_selection_count(v) != 1 ||
	    callback_state.count != 2 ||
	    !BU_STR_EQUAL(callback_state.last_path, "/adapter/path")) {
	printf("FAIL: opaque selection pick adapter apply\n");
	ret = 1;
	goto cleanup;
    }

    empty_pick = rt_view_pick_result_create_bsg();
    if (!empty_pick) {
	printf("FAIL: opaque empty pick setup\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_selection_set_pick_result_ref_bsg(v, empty_pick,
		test_selection_path_callback, &callback_state) ||
	    rt_view_selection_count_bsg(v) != 0 ||
	    rt_view_context_selection_count(v) != 0 ||
	    callback_state.count != 2) {
	printf("FAIL: opaque selection empty-pick clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_selection_clear_bsg(v) ||
	    rt_view_selection_clear_bsg(NULL) ||
	    rt_view_selection_count_bsg(v) != 0 ||
	    rt_view_context_selection_count(v) != 0) {
	printf("FAIL: BSG selection clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_pick_result_count_bsg(NULL) != 0 ||
	    rt_view_pick_result_path_bsg(NULL, 0, &rt_pick_path) ||
	    rt_view_pick_result_hit_dist_bsg(NULL, 0) >= 0.0 ||
	    rt_view_pick_result_append_path_bsg(NULL, v, 10, 20,
		"/adapter/path", 2.0) ||
	    rt_view_pick_result_append_copy_bsg(NULL, NULL, 0, 2.0) ||
	    rt_view_pick_result_context_count_bsg(NULL) != 0 ||
	    rt_view_pick_result_context_path_bsg(NULL, 0, &rt_pick_path) ||
	    rt_view_pick_result_context_hit_dist_bsg(NULL, 0) >= 0.0 ||
	    rt_view_pick_result_context_append_path_bsg(NULL, v, 10, 20,
		"/adapter/path", 2.0) ||
	    rt_view_pick_result_context_append_copy_bsg(NULL, NULL, 0, 2.0) ||
	    rt_view_pick_result_context_count(NULL) != 0 ||
	    rt_view_pick_result_context_path(NULL, 0, &rt_pick_path) ||
	    rt_view_pick_result_context_hit_dist(NULL, 0) >= 0.0 ||
	    rt_view_pick_result_context_append_copy(NULL, NULL, 0, 2.0)) {
	printf("FAIL: null opaque BSG pick adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_pick = rt_view_pick_result_create_bsg();
    if (!rt_pick ||
	    rt_view_pick_result_count_bsg(rt_pick) != 0 ||
	    !rt_view_pick_result_append_path_bsg(rt_pick, v, 30, 40,
		"/adapter/rt-path", 2.5) ||
	    rt_view_pick_result_count_bsg(rt_pick) != 1 ||
	    !rt_view_pick_result_path_bsg(rt_pick, 0, &rt_pick_path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&rt_pick_path), "/adapter/rt-path") ||
	    !fastf_equal(rt_view_pick_result_hit_dist_bsg(rt_pick, 0), 2.5)) {
	printf("FAIL: opaque BSG pick result adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_pick_ctx = rt_view_pick_result_context_create_bsg();
    if (!rt_pick_ctx ||
	    rt_view_pick_result_context_count_bsg(rt_pick_ctx) != 0 ||
	    !rt_view_pick_result_context_append_path_bsg(rt_pick_ctx, v, 31, 41,
		"/adapter/rt-context-path", 3.5) ||
	    rt_view_pick_result_context_count_bsg(rt_pick_ctx) != 1 ||
	    !rt_view_pick_result_context_path_bsg(rt_pick_ctx, 0, &rt_pick_path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&rt_pick_path),
		"/adapter/rt-context-path") ||
	    !fastf_equal(rt_view_pick_result_context_hit_dist_bsg(rt_pick_ctx, 0),
		3.5)) {
	printf("FAIL: opaque BSG pick result context adapter\n");
	ret = 1;
	goto cleanup;
    }

    neutral_pick_ctx = rt_view_pick_result_context_create();
    if (!neutral_pick_ctx ||
	    rt_view_pick_result_context_count(neutral_pick_ctx) != 0 ||
	    !rt_view_pick_result_context_append_copy(neutral_pick_ctx,
		rt_pick_ctx, 0, 4.5) ||
	    rt_view_pick_result_context_count(neutral_pick_ctx) != 1 ||
	    !rt_view_pick_result_context_path(neutral_pick_ctx, 0,
		&rt_pick_path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&rt_pick_path),
		"/adapter/rt-context-path") ||
	    !fastf_equal(rt_view_pick_result_context_hit_dist(neutral_pick_ctx,
		    0), 4.5)) {
	printf("FAIL: opaque retained pick result context adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = rt_view_context_pick_point_bsg(NULL, 0, 0, 0);
    if (rt_first_pick) {
	printf("FAIL: opaque BSG null point-pick context adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = rt_view_context_pick_nearest_bsg(NULL, 0, 0);
    if (rt_first_pick) {
	printf("FAIL: opaque BSG null nearest-pick context adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = rt_view_context_pick_rect_bsg(NULL, 0, 0, 1, 1);
    if (rt_first_pick) {
	printf("FAIL: opaque BSG null rect-pick context adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = NULL;

    rt_first_pick = rt_view_context_pick_point(NULL, 0, 0, 0);
    if (rt_first_pick) {
	printf("FAIL: opaque retained null point-pick context adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = rt_view_context_pick_nearest(NULL, 0, 0);
    if (rt_first_pick) {
	printf("FAIL: opaque retained null nearest-pick context adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = rt_view_context_pick_rect(NULL, 0, 0, 1, 1);
    if (rt_first_pick) {
	printf("FAIL: opaque retained null rect-pick context adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_first_pick = NULL;

    rt_first_pick = rt_view_pick_result_filter_first_bsg(rt_pick);
    if (!rt_first_pick ||
	    rt_view_pick_result_count_bsg(rt_first_pick) != 1 ||
	    !rt_view_pick_result_path_bsg(rt_first_pick, 0, &rt_pick_path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&rt_pick_path), "/adapter/rt-path")) {
	printf("FAIL: opaque BSG pick first-result adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_first_pick_ctx = rt_view_pick_result_context_filter_first_bsg(rt_pick_ctx);
    if (!rt_first_pick_ctx ||
	    rt_view_pick_result_context_count_bsg(rt_first_pick_ctx) != 1 ||
	    !rt_view_pick_result_context_path_bsg(rt_first_pick_ctx, 0,
		&rt_pick_path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&rt_pick_path),
		"/adapter/rt-context-path")) {
	printf("FAIL: opaque BSG pick first-result context adapter\n");
	ret = 1;
	goto cleanup;
    }

    neutral_first_pick_ctx = rt_view_pick_result_context_filter_first(
	    neutral_pick_ctx);
    if (!neutral_first_pick_ctx ||
	    rt_view_pick_result_context_count(neutral_first_pick_ctx) != 1 ||
	    !rt_view_pick_result_context_path(neutral_first_pick_ctx, 0,
		&rt_pick_path) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&rt_pick_path),
		"/adapter/rt-context-path")) {
	printf("FAIL: opaque retained pick first-result context adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_selection_set_pick_result_ref_bsg(NULL, rt_pick,
		test_selection_path_callback, &callback_state) ||
	    !rt_view_selection_set_pick_result_ref_bsg(v, rt_pick,
		test_selection_path_callback, &callback_state) ||
	    rt_view_selection_count_bsg(v) != 1 ||
	    !BU_STR_EQUAL(callback_state.last_path, "/adapter/rt-path")) {
	printf("FAIL: opaque BSG pick selection adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_selection_set_pick_result_ref_bsg(NULL, rt_pick,
		test_selection_path_callback, &callback_state) ||
	    !rt_view_context_selection_set_pick_result_ref_bsg(v, rt_pick,
		test_selection_path_callback, &callback_state) ||
	    rt_view_context_selection_count_bsg(v) != 1 ||
	    !BU_STR_EQUAL(callback_state.last_path, "/adapter/rt-path")) {
	printf("FAIL: opaque BSG pick selection context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_selection_set_pick_result_context_bsg(NULL, rt_pick_ctx,
		test_selection_path_callback, &callback_state) ||
	    !rt_view_context_selection_set_pick_result_context_bsg(v, rt_pick_ctx,
		test_selection_path_callback, &callback_state) ||
	    rt_view_context_selection_count_bsg(v) != 1 ||
	    !BU_STR_EQUAL(callback_state.last_path,
		"/adapter/rt-context-path")) {
	printf("FAIL: opaque BSG pick-result context selection adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_selection_set_pick_result_context(NULL,
		neutral_pick_ctx, test_selection_path_callback, &callback_state) ||
	    !rt_view_context_selection_set_pick_result_context(v,
		neutral_pick_ctx, test_selection_path_callback, &callback_state) ||
	    rt_view_context_selection_count_bsg(v) != 1 ||
	    !BU_STR_EQUAL(callback_state.last_path,
		"/adapter/rt-context-path")) {
	printf("FAIL: opaque retained pick-result context selection adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_selection_clear_bsg(v) ||
	    rt_view_context_selection_clear_bsg(NULL) ||
	    rt_view_context_selection_count_bsg(v) != 0) {
	printf("FAIL: opaque BSG selection clear context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_selection_clear(v) ||
	    rt_view_context_selection_clear(NULL) ||
	    rt_view_context_selection_count_bsg(v) != 0) {
	printf("FAIL: opaque retained selection clear context adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_is_independent_bsg(NULL) != bsg_view_is_independent(NULL)) {
	printf("FAIL: null BSG independent adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_scene_ref_is_null_bsg(rt_view_independent_scope_ref_bsg(NULL, 1))) {
	printf("FAIL: null BSG independent-scope ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_ref_is_null(
		rt_view_context_independent_scope_ref(NULL, 1))) {
	printf("FAIL: null neutral independent-scope ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_ref_is_null(rt_view_context_scene_root_ref(NULL))) {
	printf("FAIL: null neutral scene-root ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scene_ref_context(rt_view_context_scene_root_ref(NULL))) {
	printf("FAIL: null neutral scene-ref context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_independent_scope_is_null_bsg(NULL, 1)) {
	printf("FAIL: null BSG independent-scope null adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_independent_scope_destroy_bsg(NULL);

    if (rt_view_scene_attached_bsg(NULL) ||
	    rt_view_context_scene_attached_bsg(NULL) ||
	    rt_view_context_scene_attached(NULL) ||
	    rt_view_scene_anchor_ensure_bsg(NULL) ||
	    rt_view_context_scene_anchor_ensure_bsg(NULL) ||
	    rt_view_context_scene_anchor_ensure(NULL) ||
	    rt_view_scene_shared_bsg(NULL, v) ||
	    rt_view_context_scene_shared_bsg(NULL, v) ||
	    rt_view_context_scene_shared(NULL, v)) {
	printf("FAIL: null BSG scene-anchor adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (bsg_scene_ref_is_null(bsg_view_scene_ref_ensure(v))) {
	printf("FAIL: BSG view scope adapter test scene root setup\n");
	ret = 1;
	goto cleanup;
    }
    root_neutral = rt_view_context_scene_root_ref(v);
    direct_root_neutral = rt_view_scene_ref_make(
	    bsg_view_scene_ref(v).opaque, RT_VIEW_SCENE_BACKEND_BSG);
    if (rt_view_scene_ref_is_null(root_neutral) ||
	    !rt_view_scene_ref_equal(root_neutral, direct_root_neutral)) {
	printf("FAIL: neutral scene-root ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scene_ref_backend(root_neutral) != RT_VIEW_SCENE_BACKEND_BSG) {
	printf("FAIL: neutral scene-root ref backend adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scene_ref_context(root_neutral) != root_neutral.opaque) {
	printf("FAIL: neutral scene-root ref context adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_scene_root_ref_attach(NULL, root_neutral) ||
	    !rt_view_context_scene_root_ref_attach(v, root_neutral)) {
	printf("FAIL: neutral scene-root ref attach adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_ref_equal(rt_view_context_scene_root_ref(v),
		root_neutral)) {
	printf("FAIL: neutral scene-root ref attach roundtrip\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_attached_bsg(v) ||
	    !rt_view_context_scene_attached_bsg(v) ||
	    !rt_view_context_scene_attached(v)) {
	printf("FAIL: BSG scene-attached adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_anchor_ensure_bsg(v) ||
	    !rt_view_context_scene_anchor_ensure_bsg(v) ||
	    !rt_view_context_scene_anchor_ensure(v) ||
	    !rt_view_scene_attached_bsg(v)) {
	printf("FAIL: BSG scene-anchor ensure adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_shared_bsg(v, v) ||
	    !rt_view_context_scene_shared_bsg(v, v) ||
	    !rt_view_context_scene_shared(v, v)) {
	printf("FAIL: BSG scene-shared adapter\n");
	ret = 1;
	goto cleanup;
    }

    render_summary.item_count = -1;
    render_summary.highlighted_count = -1;
    if (rt_view_visible_render_summary_bsg(NULL, &render_summary) ||
	    rt_view_context_visible_render_summary_bsg(NULL, &render_summary) ||
	    rt_view_context_visible_render_summary(NULL,
		&neutral_render_summary) ||
	    render_summary.item_count != 0 ||
	    render_summary.highlighted_count != 0 ||
	    neutral_render_summary.item_count != 0 ||
	    neutral_render_summary.highlighted_count != 0 ||
	    rt_view_visible_render_summary_bsg(v, NULL) ||
	    rt_view_context_visible_render_summary_bsg(v, NULL) ||
	    rt_view_context_visible_render_summary(v, NULL) ||
	    rt_view_context_named_line_render_count_bsg(NULL, "missing") ||
	    rt_view_context_named_line_render_count_bsg(v, NULL) ||
	    rt_view_context_named_line_render_count(NULL, "missing") ||
	    rt_view_context_named_line_render_count(v, NULL)) {
	printf("FAIL: null BSG visible-render summary adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_visible_render_summary_bsg(v, &render_summary) ||
	    !rt_view_context_visible_render_summary_bsg(v, &render_summary) ||
	    !rt_view_context_visible_render_summary(v,
		&neutral_render_summary) ||
	    render_summary.item_count != 0 ||
	    render_summary.highlighted_count != 0 ||
	    neutral_render_summary.item_count != 0 ||
	    neutral_render_summary.highlighted_count != 0 ||
	    rt_view_context_named_line_render_count_bsg(v, "missing") ||
	    rt_view_context_named_line_render_count(v, "missing")) {
	printf("FAIL: empty BSG visible-render summary adapter\n");
	ret = 1;
	goto cleanup;
    }
    consistency.export_record_found = -1;
    consistency.render_item_found = -1;
    consistency.backend_node_found = -1;
    consistency.export_render_consistent = -1;
    consistency.export_backend_consistent = -1;
    if (rt_view_render_export_consistency_bsg(NULL, "all.g", &consistency) ||
	    rt_view_context_render_export_consistency(NULL, "all.g",
		&neutral_consistency) ||
	    consistency.export_record_found != 0 ||
	    consistency.render_item_found != 0 ||
	    consistency.backend_node_found != 0 ||
	    consistency.export_render_consistent != 0 ||
	    consistency.export_backend_consistent != 0 ||
	    neutral_consistency.export_record_found != 0 ||
	    neutral_consistency.render_item_found != 0 ||
	    neutral_consistency.backend_node_found != 0 ||
	    neutral_consistency.export_render_consistent != 0 ||
	    neutral_consistency.export_backend_consistent != 0 ||
	    rt_view_render_export_consistency_bsg(v, NULL, &consistency) ||
	    rt_view_render_export_consistency_bsg(v, "all.g", NULL) ||
	    rt_view_context_render_export_consistency(v, NULL,
		&neutral_consistency) ||
	    rt_view_context_render_export_consistency(v, "all.g", NULL)) {
	printf("FAIL: null BSG render/export consistency adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_render_export_consistency_bsg(v, "all.g", &consistency) ||
	    !rt_view_context_render_export_consistency(v, "all.g",
		&neutral_consistency) ||
	    consistency.export_record_found != 0 ||
	    consistency.render_item_found != 0 ||
	    consistency.backend_node_found != 0 ||
	    consistency.export_render_consistent != 0 ||
	    consistency.export_backend_consistent != 0 ||
	    neutral_consistency.export_record_found != 0 ||
	    neutral_consistency.render_item_found != 0 ||
	    neutral_consistency.backend_node_found != 0 ||
	    neutral_consistency.export_render_consistent != 0 ||
	    neutral_consistency.export_backend_consistent != 0) {
	printf("FAIL: empty BSG render/export consistency adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_is_independent_bsg(v) != bsg_view_is_independent(v)) {
	printf("FAIL: BSG shared independent adapter mismatch\n");
	ret = 1;
	goto cleanup;
    }

    scope = rt_view_independent_scope_ref_bsg(v, 1);
    direct_scope.opaque = bsg_view_independent_scope_ref(v, 0).opaque;
    if (rt_view_scene_ref_is_null_bsg(scope) ||
	    !rt_view_scene_ref_equal_bsg(scope, direct_scope)) {
	printf("FAIL: BSG independent-scope ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    scope_neutral = rt_view_context_independent_scope_ref(v, 1);
    direct_scope_neutral = rt_view_scene_ref_make(
	    bsg_view_independent_scope_ref(v, 0).opaque,
	    RT_VIEW_SCENE_BACKEND_BSG);
    if (rt_view_scene_ref_is_null(scope_neutral) ||
	    !rt_view_scene_ref_equal(scope_neutral, direct_scope_neutral)) {
	printf("FAIL: neutral independent-scope ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scene_ref_backend(scope_neutral) != RT_VIEW_SCENE_BACKEND_BSG) {
	printf("FAIL: neutral independent-scope ref backend adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_independent_scope_is_null_bsg(v, 0)) {
	printf("FAIL: BSG independent-scope non-null adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_independent_scope_is_null(v, 0)) {
	printf("FAIL: neutral RT independent-scope non-null adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_is_independent_bsg(v) ||
	    rt_view_is_independent_bsg(v) != bsg_view_is_independent(v) ||
	    !rt_view_context_is_independent(v) ||
	    rt_view_context_is_independent(v) != bsg_view_is_independent(v)) {
	printf("FAIL: BSG independent adapter mismatch\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_context_independent_scope_destroy(v);
    if (rt_view_is_independent_bsg(v)) {
	printf("FAIL: BSG independent-scope destroy adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_independent_scope_is_null_bsg(v, 0) ||
	    !rt_view_context_independent_scope_is_null(v, 0)) {
	printf("FAIL: BSG independent-scope destroyed-null adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_view_scene_ref_detach(v);
    if (rt_view_scene_attached_bsg(v) ||
	    rt_view_context_scene_attached(v)) {
	printf("FAIL: BSG detached scene-attached adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_ref_is_null(rt_view_context_scene_root_ref(v))) {
	printf("FAIL: BSG detached neutral scene-root ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scene_anchor_ensure_bsg(v) ||
	    !rt_view_context_scene_anchor_ensure(v) ||
	    !rt_view_scene_attached_bsg(v) ||
	    !rt_view_context_scene_attached(v)) {
	printf("FAIL: BSG detached scene-anchor ensure adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_context_line_set_create_bsg(NULL,
		"rt_view_line_set_context_adapter", 1, 2, 3) ||
	    rt_view_context_line_set_create_bsg(v, NULL, 1, 2, 3) ||
	    rt_view_context_line_set_create(NULL,
		"rt_view_line_set_context_adapter", 1, 2, 3) ||
	    rt_view_context_line_set_create(v, NULL, 1, 2, 3) ||
	    !rt_view_line_set_context_is_null_bsg(NULL) ||
	    !rt_view_line_set_context_is_null(NULL) ||
	    rt_view_line_set_context_set_points_bsg(NULL, NULL, NULL, 0) ||
	    rt_view_line_set_context_set_points(NULL, NULL, NULL, 0)) {
	printf("FAIL: null BSG line-set context adapter\n");
	ret = 1;
	goto cleanup;
    }
    line_set_ctx = rt_view_context_line_set_create_bsg(v,
	    "rt_view_line_set_context_adapter", 11, 22, 33);
    if (!line_set_ctx || rt_view_line_set_context_is_null_bsg(line_set_ctx)) {
	printf("FAIL: BSG line-set context create adapter\n");
	ret = 1;
	goto cleanup;
    }
    {
	point_t line_pts[2] = {{0.0, 0.0, 0.0}, {1.0, 1.0, 0.0}};
	int line_cmds[2] = {BSG_GEOMETRY_LINE_MOVE, BSG_GEOMETRY_LINE_DRAW};
	if (!rt_view_line_set_context_set_points_bsg(line_set_ctx,
		    (const point_t *)line_pts, line_cmds, 2) ||
		rt_view_context_named_line_render_count_bsg(v,
		    "rt_view_line_set_context_adapter") != 1 ||
		rt_view_context_named_line_render_count(v,
		    "rt_view_line_set_context_adapter") != 1 ||
		!rt_view_line_set_context_set_points_bsg(line_set_ctx,
		    NULL, NULL, 0)) {
	    printf("FAIL: BSG line-set context points adapter\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    neutral_line_set_ctx = rt_view_context_line_set_create(v,
	    "rt_view_line_set_context_adapter_neutral", 33, 22, 11);
    if (!neutral_line_set_ctx ||
	    rt_view_line_set_context_is_null(neutral_line_set_ctx)) {
	printf("FAIL: retained line-set context create adapter\n");
	ret = 1;
	goto cleanup;
    }
    {
	point_t line_pts[2] = {{2.0, 0.0, 0.0}, {3.0, 1.0, 0.0}};
	int line_cmds[2] = {BSG_GEOMETRY_LINE_MOVE, BSG_GEOMETRY_LINE_DRAW};
	if (!rt_view_line_set_context_set_points(neutral_line_set_ctx,
		    (const point_t *)line_pts, line_cmds, 2) ||
		rt_view_context_named_line_render_count_bsg(v,
		    "rt_view_line_set_context_adapter_neutral") != 1 ||
		rt_view_context_named_line_render_count(v,
		    "rt_view_line_set_context_adapter_neutral") != 1 ||
		!rt_view_line_set_context_set_points(neutral_line_set_ctx,
		    NULL, NULL, 0)) {
	    printf("FAIL: retained line-set context points adapter\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    if (rt_view_clear_bsg(v, RT_VIEW_CLEAR_VIEW_BSG) !=
	    bsg_clear(v, RT_VIEW_CLEAR_VIEW_BSG) ||
	    rt_view_context_clear_bsg(v, RT_VIEW_CLEAR_VIEW_BSG) !=
	    bsg_clear(v, RT_VIEW_CLEAR_VIEW_BSG) ||
	    rt_view_context_clear(v, RT_VIEW_CLEAR_VIEW) !=
	    bsg_clear(v, RT_VIEW_CLEAR_VIEW_BSG)) {
	printf("FAIL: BSG clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_init_copy_bsg(NULL, src, NULL);
    if (rt_view_context_init_copy_bsg(NULL, src, NULL)) {
	printf("FAIL: null BSG init-copy context adapter\n");
	ret = 1;
	goto cleanup;
    }
    src->gv_width = 321;
    src->gv_height = 123;
    src->gv_prevMouseX = 7;
    src->gv_prevMouseY = 9;
    src->gv_mouse_x = 11;
    src->gv_mouse_y = 13;
    MAT_DELTAS(src->gv_center, -3.0, -4.0, -5.0);
    bn_mat_inv(src->gv_view2model, src->gv_model2view);

    BU_ALLOC(direct, struct bsg_view);
    BU_ALLOC(adapter, struct bsg_view);
    bsg_view_init_copy(direct, src, NULL);
    if (!rt_view_context_init_copy_bsg(adapter, src, NULL)) {
	printf("FAIL: BSG init-copy context adapter failed\n");
	ret = 1;
	goto cleanup;
    }

    if (strcmp(bu_vls_cstr(&direct->gv_name), bu_vls_cstr(&adapter->gv_name)) ||
	    direct->gv_width != adapter->gv_width ||
	    direct->gv_height != adapter->gv_height ||
	    !fastf_equal(direct->gv_prevMouseX, adapter->gv_prevMouseX) ||
	    !fastf_equal(direct->gv_prevMouseY, adapter->gv_prevMouseY) ||
	    !fastf_equal(direct->gv_mouse_x, adapter->gv_mouse_x) ||
	    !fastf_equal(direct->gv_mouse_y, adapter->gv_mouse_y) ||
	    check_view_frame("BSG init-copy adapter", adapter, direct)) {
	printf("FAIL: BSG init-copy adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    rt_view_pick_result_free_bsg(empty_pick);
    rt_view_pick_result_free_bsg(pick);
    rt_view_pick_result_free_bsg(rt_first_pick);
    rt_view_pick_result_free_bsg(rt_pick);
    rt_view_pick_result_context_free_bsg(rt_first_pick_ctx);
    rt_view_pick_result_context_free_bsg(rt_pick_ctx);
    rt_view_pick_result_context_free(neutral_first_pick_ctx);
    rt_view_pick_result_context_free(neutral_pick_ctx);
    rt_view_line_set_context_destroy_bsg(line_set_ctx);
    rt_view_line_set_context_destroy(neutral_line_set_ctx);
    if (callbacks) {
	bu_ptbl_free(callbacks);
	BU_PUT(callbacks, struct bu_ptbl);
    }
    if (callbacks_replacement) {
	bu_ptbl_free(callbacks_replacement);
	BU_PUT(callbacks_replacement, struct bu_ptbl);
    }
    if (callbacks_neutral) {
	bu_ptbl_free(callbacks_neutral);
	BU_PUT(callbacks_neutral, struct bu_ptbl);
    }
    bu_vls_free(&rt_pick_path);
    rt_view_context_free(opaque_view);
    if (neutral_set)
	rt_view_set_context_destroy(neutral_set);
    if (set_initialized && lifecycle_view && lifecycle_view->vset == &set)
	rt_view_set_remove_view_bsg(&set, lifecycle_view);
    if (set_initialized && set_view)
	bsg_set_rm_view(&set, set_view);
    if (set_initialized)
	rt_view_set_context_free(&set);
    if (lifecycle_view) {
	rt_view_free_bsg(lifecycle_view);
	BU_PUT(lifecycle_view, struct bsg_view);
    }
    if (free_contents_view) {
	rt_view_free_bsg(free_contents_view);
	BU_PUT(free_contents_view, struct bsg_view);
    }
    if (release_storage_view) {
	rt_view_context_free_bsg(release_storage_view);
    }
    if (neutral_release_storage_view) {
	rt_view_context_free(neutral_release_storage_view);
    }
    free_view(set_view);
    if (v)
	bsg_view_scene_ref_detach(v);
    free_view(v);
    free_view(src);
    free_view(direct);
    free_view(adapter);
    return ret ? 1 : 0;
}

static int
test_bsg_camera_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_camera_adapter");
    mat_t actual_mat;
    mat_t expected_mat;
    point_t actual_pt;
    point_t expected_pt;
    int ret = 0;

    MAT_ZERO(actual_mat);
    if (rt_view_model2view_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG model2view adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_model2view_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG context model2view adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_model2view_get(actual_mat, NULL)) {
	printf("FAIL: null retained context model2view adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(expected_mat);
    if (memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG model2view adapter did not return identity\n");
	ret = 1;
	goto cleanup;
    }

    MAT_ZERO(actual_mat);
    if (rt_view_view2model_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG view2model adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_view2model_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG context view2model adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_view2model_get(actual_mat, NULL)) {
	printf("FAIL: null retained context view2model adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(expected_mat);
    if (memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG view2model adapter did not return identity\n");
	ret = 1;
	goto cleanup;
    }

    MAT_ZERO(actual_mat);
    if (rt_view_pmodel2view_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG pmodel2view adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_pmodel2view_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG context pmodel2view adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_pmodel2view_get(actual_mat, NULL)) {
	printf("FAIL: null retained context pmodel2view adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(expected_mat);
    if (memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG pmodel2view adapter did not return identity\n");
	ret = 1;
	goto cleanup;
    }

    MAT_ZERO(actual_mat);
    if (rt_view_pmat_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG pmat adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_pmat_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG context pmat adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_pmat_get(actual_mat, NULL)) {
	printf("FAIL: null retained context pmat adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(expected_mat);
    if (memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG pmat adapter did not return identity\n");
	ret = 1;
	goto cleanup;
    }

    MAT_ZERO(actual_mat);
    if (rt_view_rotation_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG rotation adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_rotation_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG context rotation adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_rotation_set_bsg(NULL, actual_mat)) {
	printf("FAIL: null BSG context rotation set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_rotation_get(actual_mat, NULL)) {
	printf("FAIL: null retained context rotation adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_rotation_set(NULL, actual_mat)) {
	printf("FAIL: null retained context rotation set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(expected_mat);
    if (memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG rotation adapter did not return identity\n");
	ret = 1;
	goto cleanup;
    }

    MAT_ZERO(actual_mat);
    if (rt_view_center_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG center adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_center_from_bsg(actual_mat, NULL)) {
	printf("FAIL: null BSG context center adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_center_get(actual_mat, NULL)) {
	printf("FAIL: null retained context center adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_IDN(expected_mat);
    if (memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG center adapter did not return identity\n");
	ret = 1;
	goto cleanup;
    }

    VSETALL(actual_pt, 0.0);
    if (rt_view_eye_pos_from_bsg(actual_pt, NULL)) {
	printf("FAIL: null BSG eye-pos adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_eye_pos_from_bsg(actual_pt, NULL)) {
	printf("FAIL: null BSG context eye-pos adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_eye_pos_get(actual_pt, NULL)) {
	printf("FAIL: null retained context eye-pos adapter should report no source view\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 0.0, 0.0, 1.0);
    if (!VNEAR_EQUAL(actual_pt, expected_pt, SMALL_FASTF)) {
	printf("FAIL: null BSG eye-pos adapter did not return default\n");
	ret = 1;
	goto cleanup;
    }

    bn_mat_angles(v->gv_model2view, 10.0, 20.0, 30.0);
    MAT_DELTAS(v->gv_model2view, 1.0, 2.0, 3.0);
    bn_mat_angles(v->gv_view2model, 30.0, 20.0, 10.0);
    MAT_DELTAS(v->gv_view2model, 3.0, 2.0, 1.0);
    bn_mat_angles(v->gv_pmodel2view, 40.0, 50.0, 60.0);
    MAT_DELTAS(v->gv_pmodel2view, 4.0, 5.0, 6.0);
    bn_mat_angles(v->gv_pmat, 60.0, 50.0, 40.0);
    MAT_DELTAS(v->gv_pmat, 6.0, 5.0, 4.0);
    bn_mat_angles(v->gv_rotation, 15.0, 25.0, 35.0);
    bn_mat_angles(v->gv_center, 4.0, 5.0, 6.0);
    MAT_DELTAS(v->gv_center, 7.0, 8.0, 9.0);
    VSET(v->gv_eye_pos, 3.0, 2.0, 1.0);

    if (!rt_view_model2view_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG model2view adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_model2view, sizeof(mat_t))) {
	printf("FAIL: BSG model2view adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_model2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_model2view, sizeof(mat_t))) {
	printf("FAIL: BSG context model2view adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_model2view_set_bsg(NULL, v->gv_model2view) ||
	    rt_view_model2view_set_bsg(v, NULL) ||
	    rt_view_context_model2view_set_bsg(NULL, v->gv_model2view) ||
	    rt_view_context_model2view_set_bsg(v, NULL) ||
	    rt_view_context_model2view_set(NULL, v->gv_model2view) ||
	    rt_view_context_model2view_set(v, NULL)) {
	printf("FAIL: null BSG model2view set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 11.0, 12.0, 13.0);
    MAT_DELTAS(expected_mat, 14.0, 15.0, 16.0);
    if (!rt_view_model2view_set_bsg(v, expected_mat) ||
	    !rt_view_model2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG model2view set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_model2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context model2view adapter after set\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 17.0, 18.0, 19.0);
    MAT_DELTAS(expected_mat, 20.0, 21.0, 22.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_model2view_set_bsg(v, expected_mat) ||
	    !rt_view_model2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context model2view set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 23.0, 24.0, 25.0);
    MAT_DELTAS(expected_mat, 26.0, 27.0, 28.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_model2view_set(v, expected_mat) ||
	    !rt_view_context_model2view_get(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: retained context model2view set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_view2model_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG view2model adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_view2model, sizeof(mat_t))) {
	printf("FAIL: BSG view2model adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_view2model_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_view2model, sizeof(mat_t))) {
	printf("FAIL: BSG context view2model adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_view2model_set_bsg(NULL, v->gv_view2model) ||
	    rt_view_view2model_set_bsg(v, NULL) ||
	    rt_view_context_view2model_set_bsg(NULL, v->gv_view2model) ||
	    rt_view_context_view2model_set_bsg(v, NULL) ||
	    rt_view_context_view2model_set(NULL, v->gv_view2model) ||
	    rt_view_context_view2model_set(v, NULL)) {
	printf("FAIL: null BSG view2model set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 21.0, 22.0, 23.0);
    MAT_DELTAS(expected_mat, 24.0, 25.0, 26.0);
    if (!rt_view_view2model_set_bsg(v, expected_mat) ||
	    !rt_view_view2model_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG view2model set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_view2model_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context view2model adapter after set\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 27.0, 28.0, 29.0);
    MAT_DELTAS(expected_mat, 30.0, 31.0, 32.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_view2model_set_bsg(v, expected_mat) ||
	    !rt_view_view2model_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context view2model set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 33.0, 34.0, 35.0);
    MAT_DELTAS(expected_mat, 36.0, 37.0, 38.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_view2model_set(v, expected_mat) ||
	    !rt_view_context_view2model_get(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: retained context view2model set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_pmodel2view_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG pmodel2view adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_pmodel2view, sizeof(mat_t))) {
	printf("FAIL: BSG pmodel2view adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmodel2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_pmodel2view, sizeof(mat_t))) {
	printf("FAIL: BSG context pmodel2view adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_pmodel2view_set_bsg(NULL, v->gv_pmodel2view) ||
	    rt_view_pmodel2view_set_bsg(v, NULL) ||
	    rt_view_context_pmodel2view_set_bsg(NULL, v->gv_pmodel2view) ||
	    rt_view_context_pmodel2view_set_bsg(v, NULL) ||
	    rt_view_context_pmodel2view_set(NULL, v->gv_pmodel2view) ||
	    rt_view_context_pmodel2view_set(v, NULL)) {
	printf("FAIL: null BSG pmodel2view set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 31.0, 32.0, 33.0);
    MAT_DELTAS(expected_mat, 34.0, 35.0, 36.0);
    if (!rt_view_pmodel2view_set_bsg(v, expected_mat) ||
	    !rt_view_pmodel2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG pmodel2view set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmodel2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context pmodel2view adapter after set\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 37.0, 38.0, 39.0);
    MAT_DELTAS(expected_mat, 40.0, 41.0, 42.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmodel2view_set_bsg(v, expected_mat) ||
	    !rt_view_pmodel2view_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context pmodel2view set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 43.0, 44.0, 45.0);
    MAT_DELTAS(expected_mat, 46.0, 47.0, 48.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmodel2view_set(v, expected_mat) ||
	    !rt_view_context_pmodel2view_get(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: retained context pmodel2view set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_pmat_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG pmat adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_pmat, sizeof(mat_t))) {
	printf("FAIL: BSG pmat adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmat_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_pmat, sizeof(mat_t))) {
	printf("FAIL: BSG context pmat adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_pmat_set_bsg(NULL, v->gv_pmat) ||
	    rt_view_pmat_set_bsg(v, NULL) ||
	    rt_view_context_pmat_set_bsg(NULL, v->gv_pmat) ||
	    rt_view_context_pmat_set_bsg(v, NULL) ||
	    rt_view_context_pmat_set(NULL, v->gv_pmat) ||
	    rt_view_context_pmat_set(v, NULL)) {
	printf("FAIL: null BSG pmat set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 6.0, 5.0, 4.0);
    MAT_DELTAS(expected_mat, 9.0, 8.0, 7.0);
    if (!rt_view_pmat_set_bsg(v, expected_mat) ||
	    !rt_view_pmat_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG pmat set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmat_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context pmat adapter after set\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 16.0, 15.0, 14.0);
    MAT_DELTAS(expected_mat, 19.0, 18.0, 17.0);
    if (!rt_view_context_pmat_set_bsg(v, expected_mat) ||
	    !rt_view_pmat_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context pmat set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 26.0, 25.0, 24.0);
    MAT_DELTAS(expected_mat, 29.0, 28.0, 27.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_pmat_set(v, expected_mat) ||
	    !rt_view_context_pmat_get(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: retained context pmat set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_rotation_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG rotation adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_rotation, sizeof(mat_t))) {
	printf("FAIL: BSG rotation adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_rotation_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_rotation, sizeof(mat_t))) {
	printf("FAIL: BSG context rotation adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_rotation_set_bsg(NULL, v->gv_rotation) ||
	    rt_view_rotation_set_bsg(v, NULL) ||
	    rt_view_context_rotation_set(NULL, v->gv_rotation) ||
	    rt_view_context_rotation_set(v, NULL)) {
	printf("FAIL: null BSG rotation set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 4.0, 5.0, 6.0);
    if (!rt_view_rotation_set_bsg(v, expected_mat) ||
	    !rt_view_rotation_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG rotation set adapter\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 7.0, 8.0, 9.0);
    if (!rt_view_context_rotation_set_bsg(v, expected_mat) ||
	    !rt_view_rotation_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context rotation set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_rotation_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: BSG context rotation adapter after set\n");
	ret = 1;
	goto cleanup;
    }
    bn_mat_angles(expected_mat, 10.0, 11.0, 12.0);
    MAT_ZERO(actual_mat);
    if (!rt_view_context_rotation_set(v, expected_mat) ||
	    !rt_view_context_rotation_get(actual_mat, v) ||
	    memcmp(actual_mat, expected_mat, sizeof(mat_t))) {
	printf("FAIL: retained context rotation set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_center_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG center adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_center_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_center, sizeof(mat_t))) {
	printf("FAIL: BSG context center adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_center_vec_set_bsg(NULL, v->gv_eye_pos) ||
	    rt_view_center_vec_set_bsg(v, NULL) ||
	    rt_view_context_center_vec_set_bsg(NULL, v->gv_eye_pos) ||
	    rt_view_context_center_vec_set_bsg(v, NULL) ||
	    rt_view_context_center_set(NULL, v->gv_eye_pos) ||
	    rt_view_context_center_set(v, NULL)) {
	printf("FAIL: null BSG center-vector set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 7.0, 8.0, 9.0);
    if (!rt_view_center_vec_set_bsg(v, expected_pt) ||
	    !rt_view_center_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG center-vector set adapter\n");
	ret = 1;
	goto cleanup;
    }
    point_t actual_center = VINIT_ZERO;
    MAT_DELTAS_GET_NEG(actual_center, actual_mat);
    if (!VNEAR_EQUAL(actual_center, expected_pt, SMALL_FASTF)) {
	printf("FAIL: BSG center-vector set adapter did not set expected center\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_center, sizeof(mat_t))) {
	printf("FAIL: BSG center adapter did not copy matrix\n");
	ret = 1;
	goto cleanup;
    }
    MAT_ZERO(actual_mat);
    if (!rt_view_context_center_from_bsg(actual_mat, v) ||
	    memcmp(actual_mat, v->gv_center, sizeof(mat_t))) {
	printf("FAIL: BSG context center adapter after set\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 10.0, 11.0, 12.0);
    if (!rt_view_context_center_vec_set_bsg(v, expected_pt) ||
	    !rt_view_context_center_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG context center-vector set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_DELTAS_GET_NEG(actual_center, actual_mat);
    if (!VNEAR_EQUAL(actual_center, expected_pt, SMALL_FASTF)) {
	printf("FAIL: BSG context center-vector set adapter did not set expected center\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 13.0, 14.0, 15.0);
    if (!rt_view_context_center_set(v, expected_pt) ||
	    !rt_view_context_center_get(actual_mat, v)) {
	printf("FAIL: retained context center set adapter\n");
	ret = 1;
	goto cleanup;
    }
    MAT_DELTAS_GET_NEG(actual_center, actual_mat);
    if (!VNEAR_EQUAL(actual_center, expected_pt, SMALL_FASTF)) {
	printf("FAIL: retained context center set adapter did not set expected center\n");
	ret = 1;
	goto cleanup;
    }

    plane_t actual_plane;
    plane_t expected_plane;
    point_t plane_center;
    vect_t plane_normal;
    if (rt_view_plane_from_bsg(NULL, v) != -1 ||
	    rt_view_plane_from_bsg(&actual_plane, NULL) != -1 ||
	    rt_view_context_plane_from_bsg(NULL, v) != -1 ||
	    rt_view_context_plane_from_bsg(&actual_plane, NULL) != -1 ||
	    rt_view_context_plane_get(NULL, v) != -1 ||
	    rt_view_context_plane_get(&actual_plane, NULL) != -1) {
	printf("FAIL: null BSG view-plane adapter should fail\n");
	ret = 1;
	goto cleanup;
    }
    MAT_DELTAS_GET_NEG(plane_center, v->gv_center);
    VMOVEN(plane_normal, v->gv_rotation + 8, 3);
    VUNITIZE(plane_normal);
    VSCALE(plane_normal, plane_normal, -1.0);
    bg_plane_pt_nrml(&expected_plane, plane_center, plane_normal);
    if (rt_view_plane_from_bsg(&actual_plane, v) ||
	    !VNEAR_EQUAL(actual_plane, expected_plane, SMALL_FASTF)) {
	printf("FAIL: BSG view-plane adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_plane_from_bsg(&actual_plane, v) ||
	    !VNEAR_EQUAL(actual_plane, expected_plane, SMALL_FASTF)) {
	printf("FAIL: BSG context view-plane adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_plane_get(&actual_plane, v) ||
	    !VNEAR_EQUAL(actual_plane, expected_plane, SMALL_FASTF)) {
	printf("FAIL: retained context view-plane adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_eye_pos_from_bsg(actual_pt, v)) {
	printf("FAIL: BSG eye-pos adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (!VNEAR_EQUAL(actual_pt, v->gv_eye_pos, SMALL_FASTF)) {
	printf("FAIL: BSG eye-pos adapter did not copy vector\n");
	ret = 1;
	goto cleanup;
    }
    VSETALL(actual_pt, 0.0);
    if (!rt_view_context_eye_pos_from_bsg(actual_pt, v) ||
	    !VNEAR_EQUAL(actual_pt, v->gv_eye_pos, SMALL_FASTF)) {
	printf("FAIL: BSG context eye-pos adapter did not copy vector\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_eye_pos_get(actual_pt, v) ||
	    !VNEAR_EQUAL(actual_pt, v->gv_eye_pos, SMALL_FASTF)) {
	printf("FAIL: retained context eye-pos adapter did not copy vector\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_eye_pos_set_bsg(NULL, v->gv_eye_pos) ||
	    rt_view_eye_pos_set_bsg(v, NULL) ||
	    rt_view_context_eye_pos_set_bsg(NULL, v->gv_eye_pos) ||
	    rt_view_context_eye_pos_set_bsg(v, NULL) ||
	    rt_view_context_eye_pos_set(NULL, v->gv_eye_pos) ||
	    rt_view_context_eye_pos_set(v, NULL)) {
	printf("FAIL: null BSG eye-pos set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 9.0, 8.0, 7.0);
    if (!rt_view_eye_pos_set_bsg(v, expected_pt) ||
	    !rt_view_eye_pos_from_bsg(actual_pt, v) ||
	    !VNEAR_EQUAL(actual_pt, expected_pt, SMALL_FASTF)) {
	printf("FAIL: BSG eye-pos set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 3.0, 2.0, 1.0);
    if (!rt_view_context_eye_pos_set_bsg(v, expected_pt) ||
	    !rt_view_eye_pos_from_bsg(actual_pt, v) ||
	    !VNEAR_EQUAL(actual_pt, expected_pt, SMALL_FASTF)) {
	printf("FAIL: BSG context eye-pos set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected_pt, 6.0, 5.0, 4.0);
    if (!rt_view_context_eye_pos_set(v, expected_pt) ||
	    !rt_view_context_eye_pos_get(actual_pt, v) ||
	    !VNEAR_EQUAL(actual_pt, expected_pt, SMALL_FASTF)) {
	printf("FAIL: retained context eye-pos set adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_lod_policy_adapter(void)
{
    struct bsg_view *src = make_view("rt_view_lod_policy_source");
    struct bsg_view *dst = make_view("rt_view_lod_policy_destination");
    struct bsg_lod_source_policy_settings bsg_policy;
    struct rt_view_lod_policy policy = RT_VIEW_LOD_POLICY_INIT;
    struct rt_view_lod_policy copied = RT_VIEW_LOD_POLICY_INIT;
    void *bounds_update_state = NULL;
    int ret = 0;

    if (!bsg_view_lod_source_policy_get(src, &bsg_policy)) {
	printf("FAIL: source bsg lod policy get\n");
	ret = 1;
	goto cleanup;
    }

    bsg_policy.mesh_enabled = 1;
    bsg_policy.csg_enabled = 1;
    bsg_policy.zoom_refresh = 1;
    bsg_policy.scale = 2.5;
    bsg_policy.curve_scale = 3.5;
    bsg_policy.point_scale = 4.5;
    bsg_policy.bot_threshold = 2468;
    bsg_view_lod_source_policy_set(src, &bsg_policy);

    if (rt_view_context_lod_policy_from_bsg(NULL, src) ||
	    rt_view_context_lod_policy_from_bsg(&policy, NULL) ||
	    rt_view_context_lod_policy_apply_bsg(NULL, &policy) ||
	    rt_view_context_lod_policy_apply_bsg(dst, NULL) ||
	    rt_view_context_lod_policy_copy_bsg(NULL, src) ||
	    rt_view_context_lod_policy_copy_bsg(dst, NULL) ||
	    rt_view_context_lod_policy_get(NULL, src) ||
	    rt_view_context_lod_policy_get(&policy, NULL) ||
	    rt_view_context_lod_policy_apply(NULL, &policy) ||
	    rt_view_context_lod_policy_apply(dst, NULL) ||
	    rt_view_context_lod_policy_copy(NULL, src) ||
	    rt_view_context_lod_policy_copy(dst, NULL)) {
	printf("FAIL: null rt lod policy context adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_lod_policy_from_bsg(&policy, src)) {
	printf("FAIL: rt lod policy from bsg\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("from bsg", &policy, 1, 1, 1, 2.5, 3.5, 4.5, 2468);
    if (!rt_view_context_lod_policy_from_bsg(&copied, src)) {
	printf("FAIL: rt lod policy context from bsg\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("context from bsg", &copied, 1, 1, 1, 2.5, 3.5, 4.5, 2468);
    if (!rt_view_context_lod_policy_get(&copied, src)) {
	printf("FAIL: rt lod policy neutral context get\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("context neutral get", &copied, 1, 1, 1, 2.5, 3.5, 4.5, 2468);

    policy.mesh_enabled = 0;
    policy.csg_enabled = 1;
    policy.zoom_refresh = 1;
    policy.scale = 0.0;
    policy.curve_scale = 0.0;
    policy.point_scale = 0.0;
    policy.bot_threshold = 1357;
    if (!rt_view_lod_policy_apply_bsg(dst, &policy)) {
	printf("FAIL: rt lod policy apply bsg\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_lod_policy_from_bsg(&copied, dst)) {
	printf("FAIL: rt lod policy reread after apply\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("apply bsg", &copied, 0, 1, 1, 1.0, 1.0, 1.0, 1357);

    policy.bot_threshold = 9753;
    if (!rt_view_context_lod_policy_apply_bsg(dst, &policy)) {
	printf("FAIL: rt lod policy context apply bsg\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_lod_policy_from_bsg(&copied, dst)) {
	printf("FAIL: rt lod policy context reread after apply\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("context apply bsg", &copied, 0, 1, 1, 1.0, 1.0, 1.0, 9753);

    policy.bot_threshold = 6420;
    if (!rt_view_context_lod_policy_apply(dst, &policy)) {
	printf("FAIL: rt lod policy neutral context apply\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_lod_policy_get(&copied, dst)) {
	printf("FAIL: rt lod policy neutral context reread after apply\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("context neutral apply", &copied, 0, 1, 1, 1.0, 1.0, 1.0, 6420);

    if (!rt_view_lod_policy_copy_bsg(dst, src)) {
	printf("FAIL: rt lod policy copy bsg\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_lod_policy_from_bsg(&copied, dst)) {
	printf("FAIL: rt lod policy reread after copy\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("copy bsg", &copied, 1, 1, 1, 2.5, 3.5, 4.5, 2468);

    policy.bot_threshold = 8642;
    if (!rt_view_context_lod_policy_apply_bsg(dst, &policy) ||
	    !rt_view_context_lod_policy_copy_bsg(dst, src) ||
	    !rt_view_context_lod_policy_from_bsg(&copied, dst)) {
	printf("FAIL: rt lod policy context copy bsg\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("context copy bsg", &copied, 1, 1, 1, 2.5, 3.5, 4.5, 2468);
    policy.bot_threshold = 4242;
    if (!rt_view_context_lod_policy_apply(dst, &policy) ||
	    !rt_view_context_lod_policy_copy(dst, src) ||
	    !rt_view_context_lod_policy_get(&copied, dst)) {
	printf("FAIL: rt lod policy neutral context copy\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("context neutral copy", &copied, 1, 1, 1, 2.5, 3.5, 4.5, 2468);

    rt_view_lod_bounds_update_bsg(NULL);
    if (rt_view_context_lod_bounds_update_bsg(NULL) ||
	    rt_view_context_lod_bounds_update(NULL)) {
	printf("FAIL: null BSG context lod bounds update adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_lod_bounds_update_bsg(dst);
    if (!rt_view_context_lod_bounds_update_bsg(dst) ||
	    !rt_view_context_lod_bounds_update(dst)) {
	printf("FAIL: BSG context lod bounds update adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_bounds_update_callback_from_bsg(NULL) ||
	    rt_view_context_bounds_update_callback_from_bsg(NULL)) {
	printf("FAIL: null BSG bounds callback get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_set_bsg(NULL, test_bounds_update_callback) ||
	    rt_view_context_bounds_update_callback_set_bsg(NULL,
		test_bounds_update_callback)) {
	printf("FAIL: null BSG bounds callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_call_bsg(NULL) ||
	    rt_view_context_bounds_update_callback_call_bsg(NULL)) {
	printf("FAIL: null BSG bounds callback call adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_bounds_update_suspend(NULL) ||
	    rt_view_context_bounds_update_restore(NULL, NULL, 1)) {
	printf("FAIL: null retained bounds update suspend/restore adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_from_bsg(dst) ||
	    rt_view_context_bounds_update_callback_from_bsg(dst)) {
	printf("FAIL: initial BSG bounds callback should be null\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_bounds_update_callback_set_bsg(dst, test_bounds_update_callback)) {
	printf("FAIL: BSG bounds callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_from_bsg(dst) != test_bounds_update_callback ||
	    rt_view_context_bounds_update_callback_from_bsg(dst) !=
	    test_bounds_update_callback) {
	printf("FAIL: BSG bounds callback get adapter\n");
	ret = 1;
	goto cleanup;
    }
    test_bounds_update_count = 0;
    if (!rt_view_bounds_update_callback_call_bsg(dst) || test_bounds_update_count != 1) {
	printf("FAIL: BSG bounds callback call adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_bounds_update_callback_set_bsg(dst,
		test_bounds_update_callback) ||
	    rt_view_context_bounds_update_callback_from_bsg(dst) !=
	    test_bounds_update_callback) {
	printf("FAIL: BSG context bounds callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    test_bounds_update_count = 0;
    if (!rt_view_context_bounds_update_callback_call_bsg(dst) ||
	    test_bounds_update_count != 1) {
	printf("FAIL: BSG context bounds callback call adapter\n");
	ret = 1;
	goto cleanup;
    }

    bounds_update_state = rt_view_context_bounds_update_suspend(dst);
    if (!bounds_update_state ||
	    rt_view_context_bounds_update_callback_from_bsg(dst)) {
	printf("FAIL: retained bounds update suspend adapter\n");
	ret = 1;
	goto cleanup;
    }
    test_bounds_update_count = 0;
    if (rt_view_context_bounds_update_callback_call_bsg(dst) ||
	    test_bounds_update_count != 0) {
	printf("FAIL: retained bounds update suspend clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_bounds_update_restore(dst, bounds_update_state, 1) ||
	    rt_view_context_bounds_update_callback_from_bsg(dst) !=
	    test_bounds_update_callback ||
	    test_bounds_update_count != 1) {
	printf("FAIL: retained bounds update restore adapter\n");
	ret = 1;
	goto cleanup;
    }
    bounds_update_state = NULL;

    bounds_update_state = rt_view_context_bounds_update_suspend(dst);
    if (!bounds_update_state ||
	    rt_view_context_bounds_update_callback_from_bsg(dst)) {
	printf("FAIL: retained bounds update suspend-for-null-restore adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_bounds_update_restore(NULL, bounds_update_state, 1) ||
	    rt_view_context_bounds_update_callback_from_bsg(dst)) {
	printf("FAIL: retained bounds update null restore adapter\n");
	ret = 1;
	goto cleanup;
    }
    bounds_update_state = NULL;

    if (!rt_view_bounds_update_callback_set_bsg(dst, NULL) ||
	    rt_view_bounds_update_callback_call_bsg(dst)) {
	printf("FAIL: BSG bounds callback clear adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    if (bounds_update_state)
	rt_view_context_bounds_update_restore(NULL, bounds_update_state, 0);
    free_view(src);
    free_view(dst);
    return ret ? 1 : 0;
}

static void
fill_rt_adc_state(struct rt_view_adc_state *state, int base)
{
    memset(state, 0, sizeof(*state));
    state->draw = 1;
    state->dv_x = base + 1;
    state->dv_y = base + 2;
    state->dv_a1 = base + 3;
    state->dv_a2 = base + 4;
    state->dv_dist = base + 5;
    VSET(state->pos_model, base + 6, base + 7, base + 8);
    VSET(state->pos_view, base + 9, base + 10, base + 11);
    VSET(state->pos_grid, base + 12, base + 13, base + 14);
    state->a1 = base + 15;
    state->a2 = base + 16;
    state->dst = base + 17;
    state->anchor_pos = base + 18;
    state->anchor_a1 = base + 19;
    state->anchor_a2 = base + 20;
    state->anchor_dst = base + 21;
    VSET(state->anchor_pt_a1, base + 22, base + 23, base + 24);
    VSET(state->anchor_pt_a2, base + 25, base + 26, base + 27);
    VSET(state->anchor_pt_dst, base + 28, base + 29, base + 30);
    state->line_color[0] = base + 31;
    state->line_color[1] = base + 32;
    state->line_color[2] = base + 33;
    state->tick_color[0] = base + 34;
    state->tick_color[1] = base + 35;
    state->tick_color[2] = base + 36;
    state->line_width = base + 37;
}

static int
same_rt_adc_as_bsg(const struct rt_view_adc_state *rt_state,
		   const struct bsg_adc_state *bsg_state)
{
    return rt_state->draw == bsg_state->draw &&
	rt_state->dv_x == bsg_state->dv_x &&
	rt_state->dv_y == bsg_state->dv_y &&
	rt_state->dv_a1 == bsg_state->dv_a1 &&
	rt_state->dv_a2 == bsg_state->dv_a2 &&
	rt_state->dv_dist == bsg_state->dv_dist &&
	VNEAR_EQUAL(rt_state->pos_model, bsg_state->pos_model, SMALL_FASTF) &&
	VNEAR_EQUAL(rt_state->pos_view, bsg_state->pos_view, SMALL_FASTF) &&
	VNEAR_EQUAL(rt_state->pos_grid, bsg_state->pos_grid, SMALL_FASTF) &&
	fastf_equal(rt_state->a1, bsg_state->a1) &&
	fastf_equal(rt_state->a2, bsg_state->a2) &&
	fastf_equal(rt_state->dst, bsg_state->dst) &&
	rt_state->anchor_pos == bsg_state->anchor_pos &&
	rt_state->anchor_a1 == bsg_state->anchor_a1 &&
	rt_state->anchor_a2 == bsg_state->anchor_a2 &&
	rt_state->anchor_dst == bsg_state->anchor_dst &&
	VNEAR_EQUAL(rt_state->anchor_pt_a1, bsg_state->anchor_pt_a1, SMALL_FASTF) &&
	VNEAR_EQUAL(rt_state->anchor_pt_a2, bsg_state->anchor_pt_a2, SMALL_FASTF) &&
	VNEAR_EQUAL(rt_state->anchor_pt_dst, bsg_state->anchor_pt_dst, SMALL_FASTF) &&
	rt_state->line_color[0] == bsg_state->line_color[0] &&
	rt_state->line_color[1] == bsg_state->line_color[1] &&
	rt_state->line_color[2] == bsg_state->line_color[2] &&
	rt_state->tick_color[0] == bsg_state->tick_color[0] &&
	rt_state->tick_color[1] == bsg_state->tick_color[1] &&
	rt_state->tick_color[2] == bsg_state->tick_color[2] &&
	rt_state->line_width == bsg_state->line_width;
}

static int
same_rt_adc_state(const struct rt_view_adc_state *a,
		  const struct rt_view_adc_state *b)
{
    struct bsg_adc_state tmp = {0};
    tmp.draw = b->draw;
    tmp.dv_x = b->dv_x;
    tmp.dv_y = b->dv_y;
    tmp.dv_a1 = b->dv_a1;
    tmp.dv_a2 = b->dv_a2;
    tmp.dv_dist = b->dv_dist;
    VMOVE(tmp.pos_model, b->pos_model);
    VMOVE(tmp.pos_view, b->pos_view);
    VMOVE(tmp.pos_grid, b->pos_grid);
    tmp.a1 = b->a1;
    tmp.a2 = b->a2;
    tmp.dst = b->dst;
    tmp.anchor_pos = b->anchor_pos;
    tmp.anchor_a1 = b->anchor_a1;
    tmp.anchor_a2 = b->anchor_a2;
    tmp.anchor_dst = b->anchor_dst;
    VMOVE(tmp.anchor_pt_a1, b->anchor_pt_a1);
    VMOVE(tmp.anchor_pt_a2, b->anchor_pt_a2);
    VMOVE(tmp.anchor_pt_dst, b->anchor_pt_dst);
    tmp.line_color[0] = b->line_color[0];
    tmp.line_color[1] = b->line_color[1];
    tmp.line_color[2] = b->line_color[2];
    tmp.tick_color[0] = b->tick_color[0];
    tmp.tick_color[1] = b->tick_color[1];
    tmp.tick_color[2] = b->tick_color[2];
    tmp.line_width = b->line_width;
    return same_rt_adc_as_bsg(a, &tmp);
}

static void
fill_rt_grid_state(struct rt_view_grid_state *state, int base)
{
    memset(state, 0, sizeof(*state));
    state->rc = base + 1;
    state->draw = 1;
    state->adaptive = 1;
    state->snap = 1;
    VSET(state->anchor, base + 2, base + 3, base + 4);
    state->res_h = base + 5;
    state->res_v = base + 6;
    state->res_major_h = base + 7;
    state->res_major_v = base + 8;
    state->color[0] = base + 9;
    state->color[1] = base + 10;
    state->color[2] = base + 11;
}

static int
same_rt_grid_as_bsg(const struct rt_view_grid_state *rt_state,
		    const struct bsg_grid_state *bsg_state)
{
    return rt_state->rc == bsg_state->rc &&
	rt_state->draw == bsg_state->draw &&
	rt_state->adaptive == bsg_state->adaptive &&
	rt_state->snap == bsg_state->snap &&
	VNEAR_EQUAL(rt_state->anchor, bsg_state->anchor, SMALL_FASTF) &&
	fastf_equal(rt_state->res_h, bsg_state->res_h) &&
	fastf_equal(rt_state->res_v, bsg_state->res_v) &&
	rt_state->res_major_h == bsg_state->res_major_h &&
	rt_state->res_major_v == bsg_state->res_major_v &&
	rt_state->color[0] == bsg_state->color[0] &&
	rt_state->color[1] == bsg_state->color[1] &&
	rt_state->color[2] == bsg_state->color[2];
}

static int
same_rt_grid_state(const struct rt_view_grid_state *a,
		   const struct rt_view_grid_state *b)
{
    return a->rc == b->rc &&
	a->draw == b->draw &&
	a->adaptive == b->adaptive &&
	a->snap == b->snap &&
	VNEAR_EQUAL(a->anchor, b->anchor, SMALL_FASTF) &&
	fastf_equal(a->res_h, b->res_h) &&
	fastf_equal(a->res_v, b->res_v) &&
	a->res_major_h == b->res_major_h &&
	a->res_major_v == b->res_major_v &&
	a->color[0] == b->color[0] &&
	a->color[1] == b->color[1] &&
	a->color[2] == b->color[2];
}

static void
fill_rt_axes_state(struct rt_view_axes_state *state, int base)
{
    memset(state, 0, sizeof(*state));
    state->draw = 1;
    VSET(state->axes_pos, base + 1, base + 2, base + 3);
    state->axes_size = base + 4;
    state->line_width = base + 5;
    state->axes_color[0] = base + 6;
    state->axes_color[1] = base + 7;
    state->axes_color[2] = base + 8;
    state->pos_only = base + 9;
    state->label_flag = base + 10;
    state->label_color[0] = base + 11;
    state->label_color[1] = base + 12;
    state->label_color[2] = base + 13;
    state->triple_color = base + 14;
    state->tick_enabled = base + 15;
    state->tick_length = base + 16;
    state->tick_major_length = base + 17;
    state->tick_interval = base + 18;
    state->ticks_per_major = base + 19;
    state->tick_threshold = base + 20;
    state->tick_color[0] = base + 21;
    state->tick_color[1] = base + 22;
    state->tick_color[2] = base + 23;
    state->tick_major_color[0] = base + 24;
    state->tick_major_color[1] = base + 25;
    state->tick_major_color[2] = base + 26;
}

static int
same_rt_axes_as_bsg(const struct rt_view_axes_state *rt_state,
		    const struct bsg_axes *bsg_state)
{
    return rt_state->draw == bsg_state->draw &&
	VNEAR_EQUAL(rt_state->axes_pos, bsg_state->axes_pos, SMALL_FASTF) &&
	fastf_equal(rt_state->axes_size, bsg_state->axes_size) &&
	rt_state->line_width == bsg_state->line_width &&
	rt_state->axes_color[0] == bsg_state->axes_color[0] &&
	rt_state->axes_color[1] == bsg_state->axes_color[1] &&
	rt_state->axes_color[2] == bsg_state->axes_color[2] &&
	rt_state->pos_only == bsg_state->pos_only &&
	rt_state->label_flag == bsg_state->label_flag &&
	rt_state->label_color[0] == bsg_state->label_color[0] &&
	rt_state->label_color[1] == bsg_state->label_color[1] &&
	rt_state->label_color[2] == bsg_state->label_color[2] &&
	rt_state->triple_color == bsg_state->triple_color &&
	rt_state->tick_enabled == bsg_state->tick_enabled &&
	rt_state->tick_length == bsg_state->tick_length &&
	rt_state->tick_major_length == bsg_state->tick_major_length &&
	fastf_equal(rt_state->tick_interval, bsg_state->tick_interval) &&
	rt_state->ticks_per_major == bsg_state->ticks_per_major &&
	rt_state->tick_threshold == bsg_state->tick_threshold &&
	rt_state->tick_color[0] == bsg_state->tick_color[0] &&
	rt_state->tick_color[1] == bsg_state->tick_color[1] &&
	rt_state->tick_color[2] == bsg_state->tick_color[2] &&
	rt_state->tick_major_color[0] == bsg_state->tick_major_color[0] &&
	rt_state->tick_major_color[1] == bsg_state->tick_major_color[1] &&
	rt_state->tick_major_color[2] == bsg_state->tick_major_color[2];
}

static int
same_rt_axes_state(const struct rt_view_axes_state *a,
		   const struct rt_view_axes_state *b)
{
    return a->draw == b->draw &&
	VNEAR_EQUAL(a->axes_pos, b->axes_pos, SMALL_FASTF) &&
	fastf_equal(a->axes_size, b->axes_size) &&
	a->line_width == b->line_width &&
	a->axes_color[0] == b->axes_color[0] &&
	a->axes_color[1] == b->axes_color[1] &&
	a->axes_color[2] == b->axes_color[2] &&
	a->pos_only == b->pos_only &&
	a->label_flag == b->label_flag &&
	a->label_color[0] == b->label_color[0] &&
	a->label_color[1] == b->label_color[1] &&
	a->label_color[2] == b->label_color[2] &&
	a->triple_color == b->triple_color &&
	a->tick_enabled == b->tick_enabled &&
	a->tick_length == b->tick_length &&
	a->tick_major_length == b->tick_major_length &&
	fastf_equal(a->tick_interval, b->tick_interval) &&
	a->ticks_per_major == b->ticks_per_major &&
	a->tick_threshold == b->tick_threshold &&
	a->tick_color[0] == b->tick_color[0] &&
	a->tick_color[1] == b->tick_color[1] &&
	a->tick_color[2] == b->tick_color[2] &&
	a->tick_major_color[0] == b->tick_major_color[0] &&
	a->tick_major_color[1] == b->tick_major_color[1] &&
	a->tick_major_color[2] == b->tick_major_color[2];
}

static void
fill_rt_other_state(struct rt_view_other_state *state, int base)
{
    memset(state, 0, sizeof(*state));
    state->gos_draw = 1;
    state->gos_line_color[0] = base + 1;
    state->gos_line_color[1] = base + 2;
    state->gos_line_color[2] = base + 3;
    state->gos_text_color[0] = base + 4;
    state->gos_text_color[1] = base + 5;
    state->gos_text_color[2] = base + 6;
    state->gos_font_size = base + 7;
}

static int
same_rt_other_as_bsg(const struct rt_view_other_state *rt_state,
		     const struct bsg_other_state *bsg_state)
{
    return rt_state->gos_draw == bsg_state->gos_draw &&
	rt_state->gos_line_color[0] == bsg_state->gos_line_color[0] &&
	rt_state->gos_line_color[1] == bsg_state->gos_line_color[1] &&
	rt_state->gos_line_color[2] == bsg_state->gos_line_color[2] &&
	rt_state->gos_text_color[0] == bsg_state->gos_text_color[0] &&
	rt_state->gos_text_color[1] == bsg_state->gos_text_color[1] &&
	rt_state->gos_text_color[2] == bsg_state->gos_text_color[2] &&
	rt_state->gos_font_size == bsg_state->gos_font_size;
}

static int
same_rt_other_state(const struct rt_view_other_state *a,
		    const struct rt_view_other_state *b)
{
    return a->gos_draw == b->gos_draw &&
	a->gos_line_color[0] == b->gos_line_color[0] &&
	a->gos_line_color[1] == b->gos_line_color[1] &&
	a->gos_line_color[2] == b->gos_line_color[2] &&
	a->gos_text_color[0] == b->gos_text_color[0] &&
	a->gos_text_color[1] == b->gos_text_color[1] &&
	a->gos_text_color[2] == b->gos_text_color[2] &&
	a->gos_font_size == b->gos_font_size;
}

static void
fill_rt_params_state(struct rt_view_params_state *state, int base)
{
    memset(state, 0, sizeof(*state));
    state->draw = 1;
    state->draw_size = 1;
    state->draw_center = 0;
    state->draw_az = 1;
    state->draw_el = 0;
    state->draw_tw = 1;
    state->draw_fps = 1;
    state->color[0] = base + 1;
    state->color[1] = base + 2;
    state->color[2] = base + 3;
    state->font_size = base + 4;
}

static int
same_rt_params_as_bsg(const struct rt_view_params_state *rt_state,
		      const struct bsg_params_state *bsg_state)
{
    return rt_state->draw == bsg_state->draw &&
	rt_state->draw_size == bsg_state->draw_size &&
	rt_state->draw_center == bsg_state->draw_center &&
	rt_state->draw_az == bsg_state->draw_az &&
	rt_state->draw_el == bsg_state->draw_el &&
	rt_state->draw_tw == bsg_state->draw_tw &&
	rt_state->draw_fps == bsg_state->draw_fps &&
	rt_state->color[0] == bsg_state->color[0] &&
	rt_state->color[1] == bsg_state->color[1] &&
	rt_state->color[2] == bsg_state->color[2] &&
	rt_state->font_size == bsg_state->font_size;
}

static int
same_rt_params_state(const struct rt_view_params_state *a,
		     const struct rt_view_params_state *b)
{
    return a->draw == b->draw &&
	a->draw_size == b->draw_size &&
	a->draw_center == b->draw_center &&
	a->draw_az == b->draw_az &&
	a->draw_el == b->draw_el &&
	a->draw_tw == b->draw_tw &&
	a->draw_fps == b->draw_fps &&
	a->color[0] == b->color[0] &&
	a->color[1] == b->color[1] &&
	a->color[2] == b->color[2] &&
	a->font_size == b->font_size;
}

static int
test_bsg_faceplate_state_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_faceplate_state_adapter");
    struct bsg_adc_state direct_adc = {0};
    struct rt_view_adc_state zero_rt_adc = RT_VIEW_ADC_STATE_INIT;
    struct rt_view_adc_state source_rt_adc = RT_VIEW_ADC_STATE_INIT;
    struct rt_view_adc_state adapter_rt_adc = RT_VIEW_ADC_STATE_INIT;
    struct bsg_grid_state direct_grid = {0};
    struct rt_view_grid_state zero_rt_grid = RT_VIEW_GRID_STATE_INIT;
    struct rt_view_grid_state source_rt_grid = RT_VIEW_GRID_STATE_INIT;
    struct rt_view_grid_state adapter_rt_grid = RT_VIEW_GRID_STATE_INIT;
    struct bsg_axes direct_axes = {0};
    struct rt_view_axes_state zero_rt_axes = RT_VIEW_AXES_STATE_INIT;
    struct rt_view_axes_state source_rt_axes = RT_VIEW_AXES_STATE_INIT;
    struct rt_view_axes_state adapter_rt_axes = RT_VIEW_AXES_STATE_INIT;
    struct bsg_other_state direct_other = {0};
    struct rt_view_other_state zero_rt_other = RT_VIEW_OTHER_STATE_INIT;
    struct rt_view_other_state source_rt_other = RT_VIEW_OTHER_STATE_INIT;
    struct rt_view_other_state adapter_rt_other = RT_VIEW_OTHER_STATE_INIT;
    struct bsg_params_state direct_params = {0};
    struct rt_view_params_state zero_rt_params = RT_VIEW_PARAMS_STATE_INIT;
    struct rt_view_params_state source_rt_params = RT_VIEW_PARAMS_STATE_INIT;
    struct rt_view_params_state adapter_rt_params = RT_VIEW_PARAMS_STATE_INIT;
    int ret = 0;

    fill_rt_adc_state(&source_rt_adc, 500);
    memset(&adapter_rt_adc, 0xff, sizeof(adapter_rt_adc));
    if (rt_view_adc_state_from_bsg(&adapter_rt_adc, NULL) ||
	    rt_view_context_adc_state_from_bsg(&adapter_rt_adc, NULL) ||
	    rt_view_context_adc_state_get(&adapter_rt_adc, NULL) ||
	    !same_rt_adc_state(&adapter_rt_adc, &zero_rt_adc)) {
	printf("FAIL: null RT ADC get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_adc_state_set_bsg(NULL, &source_rt_adc) ||
	    rt_view_adc_state_set_bsg(v, NULL) ||
	    rt_view_adc_state_from_bsg(NULL, v) ||
	    rt_view_context_adc_state_set_bsg(NULL, &source_rt_adc) ||
	    rt_view_context_adc_state_set_bsg(v, NULL) ||
	    rt_view_context_adc_state_from_bsg(NULL, v) ||
	    rt_view_context_adc_state_set(NULL, &source_rt_adc) ||
	    rt_view_context_adc_state_set(v, NULL) ||
	    rt_view_context_adc_state_get(NULL, v)) {
	printf("FAIL: null RT ADC adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_adc_state_set_bsg(v, &source_rt_adc) ||
	    !bsg_view_adc_get(v, &direct_adc) ||
	    !same_rt_adc_as_bsg(&source_rt_adc, &direct_adc) ||
	    !rt_view_adc_state_from_bsg(&adapter_rt_adc, v) ||
	    !same_rt_adc_state(&adapter_rt_adc, &source_rt_adc)) {
	printf("FAIL: RT ADC adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_adc_state(&source_rt_adc, 550);
    memset(&adapter_rt_adc, 0xff, sizeof(adapter_rt_adc));
    if (!rt_view_context_adc_state_set_bsg(v, &source_rt_adc) ||
	    !bsg_view_adc_get(v, &direct_adc) ||
	    !same_rt_adc_as_bsg(&source_rt_adc, &direct_adc) ||
	    !rt_view_context_adc_state_from_bsg(&adapter_rt_adc, v) ||
	    !same_rt_adc_state(&adapter_rt_adc, &source_rt_adc)) {
	printf("FAIL: RT context ADC adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_adc_state(&source_rt_adc, 575);
    memset(&adapter_rt_adc, 0xff, sizeof(adapter_rt_adc));
    if (!rt_view_context_adc_state_set(v, &source_rt_adc) ||
	    !bsg_view_adc_get(v, &direct_adc) ||
	    !same_rt_adc_as_bsg(&source_rt_adc, &direct_adc) ||
	    !rt_view_context_adc_state_get(&adapter_rt_adc, v) ||
	    !same_rt_adc_state(&adapter_rt_adc, &source_rt_adc)) {
	printf("FAIL: retained context RT ADC adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_grid_state(&source_rt_grid, 600);
    memset(&adapter_rt_grid, 0xff, sizeof(adapter_rt_grid));
    if (rt_view_grid_state_from_bsg(&adapter_rt_grid, NULL) ||
	    rt_view_context_grid_state_from_bsg(&adapter_rt_grid, NULL) ||
	    rt_view_context_grid_state_get(&adapter_rt_grid, NULL) ||
	    !same_rt_grid_state(&adapter_rt_grid, &zero_rt_grid)) {
	printf("FAIL: null RT grid get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_grid_state_set_bsg(NULL, &source_rt_grid) ||
	    rt_view_grid_state_set_bsg(v, NULL) ||
	    rt_view_grid_state_from_bsg(NULL, v) ||
	    rt_view_context_grid_state_set_bsg(NULL, &source_rt_grid) ||
	    rt_view_context_grid_state_set_bsg(v, NULL) ||
	    rt_view_context_grid_state_from_bsg(NULL, v) ||
	    rt_view_context_grid_state_set(NULL, &source_rt_grid) ||
	    rt_view_context_grid_state_set(v, NULL) ||
	    rt_view_context_grid_state_get(NULL, v)) {
	printf("FAIL: null RT grid adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_grid_state_set_bsg(v, &source_rt_grid) ||
	    !bsg_view_grid_get(v, &direct_grid) ||
	    !same_rt_grid_as_bsg(&source_rt_grid, &direct_grid) ||
	    !rt_view_grid_state_from_bsg(&adapter_rt_grid, v) ||
	    !same_rt_grid_state(&adapter_rt_grid, &source_rt_grid)) {
	printf("FAIL: RT grid adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_grid_state(&source_rt_grid, 650);
    memset(&adapter_rt_grid, 0xff, sizeof(adapter_rt_grid));
    if (!rt_view_context_grid_state_set(v, &source_rt_grid) ||
	    !bsg_view_grid_get(v, &direct_grid) ||
	    !same_rt_grid_as_bsg(&source_rt_grid, &direct_grid) ||
	    !rt_view_context_grid_state_get(&adapter_rt_grid, v) ||
	    !same_rt_grid_state(&adapter_rt_grid, &source_rt_grid)) {
	printf("FAIL: retained context RT grid adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_axes_state(&source_rt_axes, 1000);
    memset(&adapter_rt_axes, 0xff, sizeof(adapter_rt_axes));
    if (rt_view_model_axes_state_from_bsg(&adapter_rt_axes, NULL) ||
	    rt_view_context_model_axes_state_from_bsg(&adapter_rt_axes, NULL) ||
	    rt_view_context_model_axes_state_get(&adapter_rt_axes, NULL) ||
	    !same_rt_axes_state(&adapter_rt_axes, &zero_rt_axes)) {
	printf("FAIL: null RT model-axes get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_model_axes_state_set_bsg(NULL, &source_rt_axes) ||
	    rt_view_context_model_axes_state_set_bsg(NULL, &source_rt_axes) ||
	    rt_view_model_axes_state_set_bsg(v, NULL) ||
	    rt_view_context_model_axes_state_set_bsg(v, NULL) ||
	    rt_view_context_model_axes_state_from_bsg(NULL, v) ||
	    rt_view_model_axes_state_from_bsg(NULL, v) ||
	    rt_view_context_model_axes_state_set(NULL, &source_rt_axes) ||
	    rt_view_context_model_axes_state_set(v, NULL) ||
	    rt_view_context_model_axes_state_get(NULL, v)) {
	printf("FAIL: null RT model-axes adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_model_axes_state_set_bsg(v, &source_rt_axes) ||
	    !bsg_view_model_axes_get(v, &direct_axes) ||
	    !same_rt_axes_as_bsg(&source_rt_axes, &direct_axes) ||
	    !rt_view_model_axes_state_from_bsg(&adapter_rt_axes, v) ||
	    !same_rt_axes_state(&adapter_rt_axes, &source_rt_axes)) {
	printf("FAIL: RT model-axes adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_axes_state(&source_rt_axes, 1050);
    memset(&adapter_rt_axes, 0xff, sizeof(adapter_rt_axes));
    if (!rt_view_context_model_axes_state_set_bsg(v, &source_rt_axes) ||
	    !bsg_view_model_axes_get(v, &direct_axes) ||
	    !same_rt_axes_as_bsg(&source_rt_axes, &direct_axes) ||
	    !rt_view_context_model_axes_state_from_bsg(&adapter_rt_axes, v) ||
	    !same_rt_axes_state(&adapter_rt_axes, &source_rt_axes)) {
	printf("FAIL: RT model-axes context adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_axes_state(&source_rt_axes, 1075);
    memset(&adapter_rt_axes, 0xff, sizeof(adapter_rt_axes));
    if (!rt_view_context_model_axes_state_set(v, &source_rt_axes) ||
	    !bsg_view_model_axes_get(v, &direct_axes) ||
	    !same_rt_axes_as_bsg(&source_rt_axes, &direct_axes) ||
	    !rt_view_context_model_axes_state_get(&adapter_rt_axes, v) ||
	    !same_rt_axes_state(&adapter_rt_axes, &source_rt_axes)) {
	printf("FAIL: retained context RT model-axes adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_axes_state(&source_rt_axes, 1100);
    memset(&adapter_rt_axes, 0xff, sizeof(adapter_rt_axes));
    if (rt_view_view_axes_state_from_bsg(&adapter_rt_axes, NULL) ||
	    rt_view_context_view_axes_state_from_bsg(&adapter_rt_axes, NULL) ||
	    rt_view_context_view_axes_state_get(&adapter_rt_axes, NULL) ||
	    !same_rt_axes_state(&adapter_rt_axes, &zero_rt_axes)) {
	printf("FAIL: null RT view-axes get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_view_axes_state_set_bsg(NULL, &source_rt_axes) ||
	    rt_view_context_view_axes_state_set_bsg(NULL, &source_rt_axes) ||
	    rt_view_view_axes_state_set_bsg(v, NULL) ||
	    rt_view_context_view_axes_state_set_bsg(v, NULL) ||
	    rt_view_context_view_axes_state_from_bsg(NULL, v) ||
	    rt_view_view_axes_state_from_bsg(NULL, v) ||
	    rt_view_context_view_axes_state_set(NULL, &source_rt_axes) ||
	    rt_view_context_view_axes_state_set(v, NULL) ||
	    rt_view_context_view_axes_state_get(NULL, v)) {
	printf("FAIL: null RT view-axes adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_view_axes_state_set_bsg(v, &source_rt_axes) ||
	    !bsg_view_view_axes_get(v, &direct_axes) ||
	    !same_rt_axes_as_bsg(&source_rt_axes, &direct_axes) ||
	    !rt_view_view_axes_state_from_bsg(&adapter_rt_axes, v) ||
	    !same_rt_axes_state(&adapter_rt_axes, &source_rt_axes)) {
	printf("FAIL: RT view-axes adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_axes_state(&source_rt_axes, 1150);
    memset(&adapter_rt_axes, 0xff, sizeof(adapter_rt_axes));
    if (!rt_view_context_view_axes_state_set_bsg(v, &source_rt_axes) ||
	    !bsg_view_view_axes_get(v, &direct_axes) ||
	    !same_rt_axes_as_bsg(&source_rt_axes, &direct_axes) ||
	    !rt_view_context_view_axes_state_from_bsg(&adapter_rt_axes, v) ||
	    !same_rt_axes_state(&adapter_rt_axes, &source_rt_axes)) {
	printf("FAIL: RT view-axes context adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_axes_state(&source_rt_axes, 1175);
    memset(&adapter_rt_axes, 0xff, sizeof(adapter_rt_axes));
    if (!rt_view_context_view_axes_state_set(v, &source_rt_axes) ||
	    !bsg_view_view_axes_get(v, &direct_axes) ||
	    !same_rt_axes_as_bsg(&source_rt_axes, &direct_axes) ||
	    !rt_view_context_view_axes_state_get(&adapter_rt_axes, v) ||
	    !same_rt_axes_state(&adapter_rt_axes, &source_rt_axes)) {
	printf("FAIL: retained context RT view-axes adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_other_state(&source_rt_other, 700);
    memset(&adapter_rt_other, 0xff, sizeof(adapter_rt_other));
    if (rt_view_center_dot_state_from_bsg(&adapter_rt_other, NULL) ||
	    rt_view_context_center_dot_state_from_bsg(&adapter_rt_other, NULL) ||
	    rt_view_context_center_dot_state_get(&adapter_rt_other, NULL) ||
	    !same_rt_other_state(&adapter_rt_other, &zero_rt_other)) {
	printf("FAIL: null RT center-dot get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_center_dot_state_set_bsg(NULL, &source_rt_other) ||
	    rt_view_context_center_dot_state_set_bsg(NULL, &source_rt_other) ||
	    rt_view_center_dot_state_set_bsg(v, NULL) ||
	    rt_view_context_center_dot_state_set_bsg(v, NULL) ||
	    rt_view_context_center_dot_state_from_bsg(NULL, v) ||
	    rt_view_center_dot_state_from_bsg(NULL, v) ||
	    rt_view_context_center_dot_state_set(NULL, &source_rt_other) ||
	    rt_view_context_center_dot_state_set(v, NULL) ||
	    rt_view_context_center_dot_state_get(NULL, v)) {
	printf("FAIL: null RT center-dot adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_center_dot_state_set_bsg(v, &source_rt_other) ||
	    !bsg_view_center_dot_get(v, &direct_other) ||
	    !same_rt_other_as_bsg(&source_rt_other, &direct_other) ||
	    !rt_view_center_dot_state_from_bsg(&adapter_rt_other, v) ||
	    !same_rt_other_state(&adapter_rt_other, &source_rt_other)) {
	printf("FAIL: RT center-dot adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_other_state(&source_rt_other, 750);
    memset(&adapter_rt_other, 0xff, sizeof(adapter_rt_other));
    if (!rt_view_context_center_dot_state_set_bsg(v, &source_rt_other) ||
	    !bsg_view_center_dot_get(v, &direct_other) ||
	    !same_rt_other_as_bsg(&source_rt_other, &direct_other) ||
	    !rt_view_context_center_dot_state_from_bsg(&adapter_rt_other, v) ||
	    !same_rt_other_state(&adapter_rt_other, &source_rt_other)) {
	printf("FAIL: RT center-dot context adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_other_state(&source_rt_other, 775);
    memset(&adapter_rt_other, 0xff, sizeof(adapter_rt_other));
    if (!rt_view_context_center_dot_state_set(v, &source_rt_other) ||
	    !bsg_view_center_dot_get(v, &direct_other) ||
	    !same_rt_other_as_bsg(&source_rt_other, &direct_other) ||
	    !rt_view_context_center_dot_state_get(&adapter_rt_other, v) ||
	    !same_rt_other_state(&adapter_rt_other, &source_rt_other)) {
	printf("FAIL: retained context RT center-dot adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_other_state(&source_rt_other, 800);
    memset(&adapter_rt_other, 0xff, sizeof(adapter_rt_other));
    if (rt_view_scale_overlay_state_from_bsg(&adapter_rt_other, NULL) ||
	    rt_view_context_scale_overlay_state_from_bsg(&adapter_rt_other, NULL) ||
	    rt_view_context_scale_overlay_state_get(&adapter_rt_other, NULL) ||
	    !same_rt_other_state(&adapter_rt_other, &zero_rt_other)) {
	printf("FAIL: null RT scale overlay get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scale_overlay_state_set_bsg(NULL, &source_rt_other) ||
	    rt_view_context_scale_overlay_state_set_bsg(NULL, &source_rt_other) ||
	    rt_view_scale_overlay_state_set_bsg(v, NULL) ||
	    rt_view_context_scale_overlay_state_set_bsg(v, NULL) ||
	    rt_view_context_scale_overlay_state_from_bsg(NULL, v) ||
	    rt_view_scale_overlay_state_from_bsg(NULL, v) ||
	    rt_view_context_scale_overlay_state_set(NULL, &source_rt_other) ||
	    rt_view_context_scale_overlay_state_set(v, NULL) ||
	    rt_view_context_scale_overlay_state_get(NULL, v)) {
	printf("FAIL: null RT scale overlay adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scale_overlay_state_set_bsg(v, &source_rt_other) ||
	    !bsg_view_scale_state_get(v, &direct_other) ||
	    !same_rt_other_as_bsg(&source_rt_other, &direct_other) ||
	    !rt_view_scale_overlay_state_from_bsg(&adapter_rt_other, v) ||
	    !same_rt_other_state(&adapter_rt_other, &source_rt_other)) {
	printf("FAIL: RT scale overlay adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_other_state(&source_rt_other, 850);
    memset(&adapter_rt_other, 0xff, sizeof(adapter_rt_other));
    if (!rt_view_context_scale_overlay_state_set_bsg(v, &source_rt_other) ||
	    !bsg_view_scale_state_get(v, &direct_other) ||
	    !same_rt_other_as_bsg(&source_rt_other, &direct_other) ||
	    !rt_view_context_scale_overlay_state_from_bsg(&adapter_rt_other, v) ||
	    !same_rt_other_state(&adapter_rt_other, &source_rt_other)) {
	printf("FAIL: RT scale overlay context adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_other_state(&source_rt_other, 875);
    memset(&adapter_rt_other, 0xff, sizeof(adapter_rt_other));
    if (!rt_view_context_scale_overlay_state_set(v, &source_rt_other) ||
	    !bsg_view_scale_state_get(v, &direct_other) ||
	    !same_rt_other_as_bsg(&source_rt_other, &direct_other) ||
	    !rt_view_context_scale_overlay_state_get(&adapter_rt_other, v) ||
	    !same_rt_other_state(&adapter_rt_other, &source_rt_other)) {
	printf("FAIL: retained context RT scale overlay adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_params_state(&source_rt_params, 900);
    memset(&adapter_rt_params, 0xff, sizeof(adapter_rt_params));
    if (rt_view_params_state_from_bsg(&adapter_rt_params, NULL) ||
	    rt_view_context_params_state_from_bsg(&adapter_rt_params, NULL) ||
	    rt_view_context_params_state_get(&adapter_rt_params, NULL) ||
	    !same_rt_params_state(&adapter_rt_params, &zero_rt_params)) {
	printf("FAIL: null RT params get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_params_state_set_bsg(NULL, &source_rt_params) ||
	    rt_view_context_params_state_set_bsg(NULL, &source_rt_params) ||
	    rt_view_params_state_set_bsg(v, NULL) ||
	    rt_view_context_params_state_set_bsg(v, NULL) ||
	    rt_view_context_params_state_from_bsg(NULL, v) ||
	    rt_view_params_state_from_bsg(NULL, v) ||
	    rt_view_context_params_state_set(NULL, &source_rt_params) ||
	    rt_view_context_params_state_set(v, NULL) ||
	    rt_view_context_params_state_get(NULL, v)) {
	printf("FAIL: null RT params adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_params_state_set_bsg(v, &source_rt_params) ||
	    !bsg_view_params_get(v, &direct_params) ||
	    !same_rt_params_as_bsg(&source_rt_params, &direct_params) ||
	    !rt_view_params_state_from_bsg(&adapter_rt_params, v) ||
	    !same_rt_params_state(&adapter_rt_params, &source_rt_params)) {
	printf("FAIL: RT params adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_rt_params_state(&source_rt_params, 950);
    memset(&adapter_rt_params, 0xff, sizeof(adapter_rt_params));
    if (!rt_view_context_params_state_set_bsg(v, &source_rt_params) ||
	    !bsg_view_params_get(v, &direct_params) ||
	    !same_rt_params_as_bsg(&source_rt_params, &direct_params) ||
	    !rt_view_context_params_state_from_bsg(&adapter_rt_params, v) ||
	    !same_rt_params_state(&adapter_rt_params, &source_rt_params)) {
	printf("FAIL: RT params context adapter\n");
	ret = 1;
	goto cleanup;
    }
    fill_rt_params_state(&source_rt_params, 975);
    memset(&adapter_rt_params, 0xff, sizeof(adapter_rt_params));
    if (!rt_view_context_params_state_set(v, &source_rt_params) ||
	    !bsg_view_params_get(v, &direct_params) ||
	    !same_rt_params_as_bsg(&source_rt_params, &direct_params) ||
	    !rt_view_context_params_state_get(&adapter_rt_params, v) ||
	    !same_rt_params_state(&adapter_rt_params, &source_rt_params)) {
	printf("FAIL: retained context RT params adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_scene_ref_adapter_boundary(void)
{
    struct bsg_view *v = make_view("rt_view_lod_bounds_callback");
    struct bsg_view *shared = make_view("rt_view_shared_settings_adapter");
    rt_view_scene_ref neutral_null_ref = rt_view_scene_ref_null();
    rt_view_scene_ref neutral_static_null_ref = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref neutral_obol_ref = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref neutral_obol_ref_copy = RT_VIEW_SCENE_REF_NULL_INIT;
    rt_view_scene_ref neutral_bsg_same_pointer = RT_VIEW_SCENE_REF_NULL_INIT;
    fastf_t *scale_storage = NULL;
    int ret = 0;

    if (!rt_view_scene_ref_is_null(neutral_null_ref) ||
	    !rt_view_scene_ref_equal(neutral_null_ref,
		neutral_static_null_ref)) {
	printf("FAIL: null neutral scene-ref adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_scene_ref_backend(neutral_null_ref) !=
	    RT_VIEW_SCENE_BACKEND_NONE) {
	printf("FAIL: null neutral scene-ref backend adapter\n");
	ret = 1;
	goto cleanup;
    }

    neutral_obol_ref = rt_view_scene_ref_make(v, RT_VIEW_SCENE_BACKEND_OBOL);
    neutral_obol_ref_copy = rt_view_scene_ref_make(v,
	    RT_VIEW_SCENE_BACKEND_OBOL);
    neutral_bsg_same_pointer = rt_view_scene_ref_make(v,
	    RT_VIEW_SCENE_BACKEND_BSG);
    if (rt_view_scene_ref_is_null(neutral_obol_ref) ||
	    rt_view_scene_ref_backend(neutral_obol_ref) !=
		RT_VIEW_SCENE_BACKEND_OBOL ||
	    rt_view_scene_ref_context(neutral_obol_ref) != v ||
	    !rt_view_scene_ref_equal(neutral_obol_ref,
		neutral_obol_ref_copy) ||
	    rt_view_scene_ref_equal(neutral_obol_ref,
		neutral_bsg_same_pointer)) {
	printf("FAIL: tagged neutral scene-ref adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_context_scene_root_ref_attach(v, neutral_obol_ref) != 0) {
	printf("FAIL: BSG scene-root attach accepted non-BSG neutral ref\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_scale_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_scale_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_scale_get(NULL), 1.0)) {
	printf("FAIL: null BSG view scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_scale = 42.0;
    if (!fastf_equal(rt_view_scale_from_bsg(v), 42.0) ||
	    !fastf_equal(rt_view_context_scale_from_bsg(v), 42.0) ||
	    !fastf_equal(rt_view_context_scale_get(v), 42.0)) {
	printf("FAIL: BSG view scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scale_set_bsg(NULL, 25.0) ||
	    rt_view_context_scale_set_bsg(NULL, 25.0) ||
	    rt_view_context_scale_set(NULL, 25.0)) {
	printf("FAIL: null BSG view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scale_set_bsg(v, 25.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 25.0)) {
	printf("FAIL: BSG view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_scale_set_bsg(v, 27.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 27.0) ||
	    !fastf_equal(rt_view_context_scale_get(v), 27.0)) {
	printf("FAIL: BSG context view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_scale_set(v, 29.0) ||
	    !fastf_equal(rt_view_context_scale_get(v), 29.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 29.0)) {
	printf("FAIL: retained context view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scale_storage_from_bsg(NULL)) {
	printf("FAIL: null BSG view scale storage adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_scale_storage_from_bsg(NULL) ||
	    rt_view_context_scale_storage_get(NULL)) {
	printf("FAIL: null BSG context view scale storage adapter\n");
	ret = 1;
	goto cleanup;
    }
    scale_storage = rt_view_scale_storage_from_bsg(v);
    if (scale_storage != &v->gv_scale) {
	printf("FAIL: BSG view scale storage adapter did not return live storage\n");
	ret = 1;
	goto cleanup;
    }
    *scale_storage = 43.0;
    if (!fastf_equal(rt_view_scale_from_bsg(v), 43.0)) {
	printf("FAIL: BSG view scale storage adapter did not update live scale\n");
	ret = 1;
	goto cleanup;
    }
    scale_storage = rt_view_context_scale_storage_from_bsg(v);
    if (scale_storage != &v->gv_scale) {
	printf("FAIL: BSG context view scale storage adapter did not return live storage\n");
	ret = 1;
	goto cleanup;
    }
    *scale_storage = 44.0;
    if (!fastf_equal(rt_view_scale_from_bsg(v), 44.0)) {
	printf("FAIL: BSG context view scale storage adapter did not update live scale\n");
	ret = 1;
	goto cleanup;
    }
    scale_storage = rt_view_context_scale_storage_get(v);
    if (scale_storage != &v->gv_scale) {
	printf("FAIL: retained context view scale storage adapter did not return live storage\n");
	ret = 1;
	goto cleanup;
    }
    *scale_storage = 45.0;
    if (!fastf_equal(rt_view_context_scale_get(v), 45.0)) {
	printf("FAIL: retained context view scale storage adapter did not update live scale\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scale_state_set_bsg(NULL, 1.0, 2.0, 3.0, 4.0, 5.0)) {
	printf("FAIL: null BSG view scale state set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_scale_state_set(NULL, 1.0, 2.0, 3.0, 4.0, 5.0)) {
	printf("FAIL: null retained context view scale state set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scale_state_set_bsg(v, 10.0, 20.0, 30.0, 40.0, 50.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 10.0) ||
	    !fastf_equal(rt_view_initial_scale_from_bsg(v), 20.0) ||
	    !fastf_equal(rt_view_absolute_scale_from_bsg(v), 30.0) ||
	    !fastf_equal(rt_view_size_from_bsg(v), 40.0) ||
	    !fastf_equal(rt_view_inverse_size_from_bsg(v), 50.0)) {
	printf("FAIL: BSG view scale state set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_scale_state_set(v, 11.0, 21.0, 31.0, 41.0, 51.0) ||
	    !fastf_equal(rt_view_context_scale_get(v), 11.0) ||
	    !fastf_equal(rt_view_context_initial_scale_get(v), 21.0) ||
	    !fastf_equal(rt_view_context_absolute_scale_get(v), 31.0) ||
	    !fastf_equal(rt_view_context_size_get(v), 41.0) ||
	    !fastf_equal(rt_view_context_inverse_size_get(v), 51.0)) {
	printf("FAIL: retained context view scale state set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_initial_scale_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_initial_scale_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_initial_scale_get(NULL), 1.0)) {
	printf("FAIL: null BSG initial view scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_i_scale = 500.0;
    if (!fastf_equal(rt_view_initial_scale_from_bsg(v), 500.0) ||
	    !fastf_equal(rt_view_context_initial_scale_from_bsg(v), 500.0) ||
	    !fastf_equal(rt_view_context_initial_scale_get(v), 500.0)) {
	printf("FAIL: BSG initial view scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_initial_scale_set_bsg(NULL, 250.0) ||
	    rt_view_context_initial_scale_set_bsg(NULL, 250.0) ||
	    rt_view_context_initial_scale_set(NULL, 250.0)) {
	printf("FAIL: null BSG initial view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_initial_scale_set_bsg(v, 250.0) ||
	    !fastf_equal(rt_view_initial_scale_from_bsg(v), 250.0)) {
	printf("FAIL: BSG initial view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_initial_scale_set_bsg(v, 275.0) ||
	    !fastf_equal(rt_view_context_initial_scale_from_bsg(v), 275.0) ||
	    !fastf_equal(rt_view_context_initial_scale_get(v), 275.0)) {
	printf("FAIL: BSG context initial view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_initial_scale_set(v, 300.0) ||
	    !fastf_equal(rt_view_context_initial_scale_get(v), 300.0) ||
	    !fastf_equal(rt_view_initial_scale_from_bsg(v), 300.0)) {
	printf("FAIL: retained context initial view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_absolute_scale_from_bsg(NULL), 0.0) ||
	    !fastf_equal(rt_view_context_absolute_scale_from_bsg(NULL), 0.0) ||
	    !fastf_equal(rt_view_context_absolute_scale_get(NULL), 0.0)) {
	printf("FAIL: null BSG absolute view scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_a_scale = -0.25;
    if (!fastf_equal(rt_view_absolute_scale_from_bsg(v), -0.25) ||
	    !fastf_equal(rt_view_context_absolute_scale_from_bsg(v), -0.25) ||
	    !fastf_equal(rt_view_context_absolute_scale_get(v), -0.25)) {
	printf("FAIL: BSG absolute view scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_absolute_scale_set_bsg(NULL, 0.25) ||
	    rt_view_context_absolute_scale_set_bsg(NULL, 0.25) ||
	    rt_view_context_absolute_scale_set(NULL, 0.25)) {
	printf("FAIL: null BSG absolute view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_absolute_scale_set_bsg(v, 0.25) ||
	    !fastf_equal(rt_view_absolute_scale_from_bsg(v), 0.25)) {
	printf("FAIL: BSG absolute view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_absolute_scale_set_bsg(v, 0.75) ||
	    !fastf_equal(rt_view_context_absolute_scale_from_bsg(v), 0.75) ||
	    !fastf_equal(rt_view_context_absolute_scale_get(v), 0.75)) {
	printf("FAIL: BSG context absolute view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_absolute_scale_set(v, 1.25) ||
	    !fastf_equal(rt_view_context_absolute_scale_get(v), 1.25) ||
	    !fastf_equal(rt_view_absolute_scale_from_bsg(v), 1.25)) {
	printf("FAIL: retained context absolute view scale set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_local2base_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_local2base_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_local2base_get(NULL), 1.0) ||
	    !fastf_equal(rt_view_base2local_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_base2local_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_base2local_get(NULL), 1.0)) {
	printf("FAIL: null BSG unit-conversion scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_local2base = 0.5;
    v->gv_base2local = 2.0;
    if (!fastf_equal(rt_view_local2base_from_bsg(v), 0.5) ||
	    !fastf_equal(rt_view_context_local2base_from_bsg(v), 0.5) ||
	    !fastf_equal(rt_view_context_local2base_get(v), 0.5) ||
	    !fastf_equal(rt_view_base2local_from_bsg(v), 2.0) ||
	    !fastf_equal(rt_view_context_base2local_from_bsg(v), 2.0) ||
	    !fastf_equal(rt_view_context_base2local_get(v), 2.0)) {
	printf("FAIL: BSG unit-conversion scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_unit_conversion_set_bsg(NULL, 3.0, 4.0)) {
	printf("FAIL: null BSG unit-conversion set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_unit_conversion_set_bsg(NULL, 3.0, 4.0)) {
	printf("FAIL: null BSG context unit-conversion set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_unit_conversion_set(NULL, 3.0, 4.0)) {
	printf("FAIL: null retained context unit-conversion set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_unit_conversion_set_bsg(v, 3.0, 4.0) ||
	    !fastf_equal(rt_view_local2base_from_bsg(v), 3.0) ||
	    !fastf_equal(rt_view_base2local_from_bsg(v), 4.0)) {
	printf("FAIL: BSG unit-conversion set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_unit_conversion_set_bsg(v, 5.0, 6.0) ||
	    !fastf_equal(rt_view_local2base_from_bsg(v), 5.0) ||
	    !fastf_equal(rt_view_context_local2base_from_bsg(v), 5.0) ||
	    !fastf_equal(rt_view_base2local_from_bsg(v), 6.0) ||
	    !fastf_equal(rt_view_context_base2local_from_bsg(v), 6.0)) {
	printf("FAIL: BSG context unit-conversion set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_unit_conversion_set(v, 7.0, 8.0) ||
	    !fastf_equal(rt_view_context_local2base_get(v), 7.0) ||
	    !fastf_equal(rt_view_context_base2local_get(v), 8.0)) {
	printf("FAIL: retained context unit-conversion set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_context_frametime_set_bsg(NULL, 1000000000) ||
	    rt_view_context_frametime_set(NULL, 1000000000)) {
	printf("FAIL: null BSG context frametime set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_frametime_set_bsg(v, 1000000000) ||
	    bsg_view_frametime(v) != 1000000000) {
	printf("FAIL: BSG context frametime set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_frametime_set(v, 2000000000) ||
	    bsg_view_frametime(v) != 2000000000) {
	printf("FAIL: retained context frametime set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_size_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG view size adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!fastf_equal(rt_view_context_size_from_bsg(NULL), 0.0) ||
	    !fastf_equal(rt_view_context_size_get(NULL), 0.0)) {
	printf("FAIL: null BSG context view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_size = 84.0;
    if (!fastf_equal(rt_view_size_from_bsg(v), 84.0) ||
	    !fastf_equal(rt_view_context_size_get(v), 84.0)) {
	printf("FAIL: BSG view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_size_set_bsg(NULL, 128.0) ||
	    rt_view_context_size_set_bsg(NULL, 128.0) ||
	    rt_view_context_size_set(NULL, 128.0)) {
	printf("FAIL: null BSG view size set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_size_set_bsg(v, 128.0) ||
	    !fastf_equal(rt_view_size_from_bsg(v), 128.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 64.0) ||
	    !fastf_equal(rt_view_inverse_size_from_bsg(v), 1.0 / 128.0)) {
	printf("FAIL: BSG view size set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_size_set_bsg(v, 256.0) ||
	    !fastf_equal(rt_view_context_size_from_bsg(v), 256.0) ||
	    !fastf_equal(rt_view_context_size_get(v), 256.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 128.0) ||
	    !fastf_equal(rt_view_inverse_size_from_bsg(v), 1.0 / 256.0)) {
	printf("FAIL: BSG context view size set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_size_set(v, 512.0) ||
	    !fastf_equal(rt_view_context_size_get(v), 512.0) ||
	    !fastf_equal(rt_view_scale_from_bsg(v), 256.0) ||
	    !fastf_equal(rt_view_context_inverse_size_get(v), 1.0 / 512.0)) {
	printf("FAIL: retained context view size set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_inverse_size_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_inverse_size_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_inverse_size_get(NULL), 1.0)) {
	printf("FAIL: null BSG inverse view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_isize = 1.0 / 84.0;
    if (!fastf_equal(rt_view_inverse_size_from_bsg(v), 1.0 / 84.0) ||
	    !fastf_equal(rt_view_context_inverse_size_from_bsg(v), 1.0 / 84.0) ||
	    !fastf_equal(rt_view_context_inverse_size_get(v), 1.0 / 84.0)) {
	printf("FAIL: BSG inverse view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    point_t keypoint = VINIT_ZERO;
    if (rt_view_keypoint_from_bsg(keypoint, NULL) ||
	    rt_view_context_keypoint_from_bsg(keypoint, NULL) ||
	    rt_view_context_keypoint_get(keypoint, NULL) ||
	    !fastf_equal(keypoint[X], 0.0) ||
	    !fastf_equal(keypoint[Y], 0.0) ||
	    !fastf_equal(keypoint[Z], 0.0)) {
	printf("FAIL: null BSG keypoint adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSET(v->gv_keypoint, 1.0, 2.0, 3.0);
    if (!rt_view_keypoint_from_bsg(keypoint, v) ||
	    !rt_view_context_keypoint_from_bsg(keypoint, v) ||
	    !rt_view_context_keypoint_get(keypoint, v) ||
	    !fastf_equal(keypoint[X], 1.0) ||
	    !fastf_equal(keypoint[Y], 2.0) ||
	    !fastf_equal(keypoint[Z], 3.0)) {
	printf("FAIL: BSG keypoint adapter\n");
	ret = 1;
	goto cleanup;
    }

    point_t set_keypoint;
    VSET(set_keypoint, 4.0, 5.0, 6.0);
    if (rt_view_keypoint_set_bsg(NULL, set_keypoint) ||
	    rt_view_context_keypoint_set_bsg(NULL, set_keypoint) ||
	    rt_view_context_keypoint_set(NULL, set_keypoint) ||
	    rt_view_keypoint_set_bsg(v, NULL) ||
	    rt_view_context_keypoint_set_bsg(v, NULL) ||
	    rt_view_context_keypoint_set(v, NULL)) {
	printf("FAIL: null BSG keypoint set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_keypoint_set_bsg(v, set_keypoint) ||
	    !rt_view_keypoint_from_bsg(keypoint, v) ||
	    !fastf_equal(keypoint[X], 4.0) ||
	    !fastf_equal(keypoint[Y], 5.0) ||
	    !fastf_equal(keypoint[Z], 6.0)) {
	printf("FAIL: BSG keypoint set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(set_keypoint, 7.0, 8.0, 9.0);
    if (!rt_view_context_keypoint_set_bsg(v, set_keypoint) ||
	    !rt_view_context_keypoint_from_bsg(keypoint, v) ||
	    !fastf_equal(keypoint[X], 7.0) ||
	    !fastf_equal(keypoint[Y], 8.0) ||
	    !fastf_equal(keypoint[Z], 9.0)) {
	printf("FAIL: BSG context keypoint set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(set_keypoint, 10.0, 11.0, 12.0);
    if (!rt_view_context_keypoint_set(v, set_keypoint) ||
	    !rt_view_context_keypoint_get(keypoint, v) ||
	    !fastf_equal(keypoint[X], 10.0) ||
	    !fastf_equal(keypoint[Y], 11.0) ||
	    !fastf_equal(keypoint[Z], 12.0)) {
	printf("FAIL: retained context keypoint set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_rotate_about_from_bsg(NULL) != 'v' ||
	    rt_view_context_rotate_about_from_bsg(NULL) != 'v' ||
	    rt_view_context_rotate_about_get(NULL) != 'v') {
	printf("FAIL: null BSG rotate-about adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_rotate_about = 'k';
    if (rt_view_rotate_about_from_bsg(v) != 'k' ||
	    rt_view_context_rotate_about_from_bsg(v) != 'k' ||
	    rt_view_context_rotate_about_get(v) != 'k') {
	printf("FAIL: BSG rotate-about adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_rotate_about_set_bsg(NULL, 'm') ||
	    rt_view_context_rotate_about_set_bsg(NULL, 'm') ||
	    rt_view_context_rotate_about_set(NULL, 'm')) {
	printf("FAIL: null BSG rotate-about set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_rotate_about_set_bsg(v, 'm') ||
	    rt_view_rotate_about_from_bsg(v) != 'm') {
	printf("FAIL: BSG rotate-about set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_rotate_about_set_bsg(v, 'e') ||
	    rt_view_rotate_about_from_bsg(v) != 'e' ||
	    rt_view_context_rotate_about_from_bsg(v) != 'e') {
	printf("FAIL: BSG context rotate-about set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_rotate_about_set(v, 'k') ||
	    rt_view_context_rotate_about_get(v) != 'k') {
	printf("FAIL: retained context rotate-about set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_coord_from_bsg(NULL) != 'v' ||
	    rt_view_context_coord_from_bsg(NULL) != 'v' ||
	    rt_view_context_coord_get(NULL) != 'v') {
	printf("FAIL: null BSG coord adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_coord = 'm';
    if (rt_view_coord_from_bsg(v) != 'm' ||
	    rt_view_context_coord_from_bsg(v) != 'm' ||
	    rt_view_context_coord_get(v) != 'm') {
	printf("FAIL: BSG coord adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_coord_set_bsg(NULL, 'v') ||
	    rt_view_context_coord_set_bsg(NULL, 'v') ||
	    rt_view_context_coord_set(NULL, 'v')) {
	printf("FAIL: null BSG coord set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_coord_set_bsg(v, 'v') ||
	    rt_view_coord_from_bsg(v) != 'v' ||
	    rt_view_context_coord_from_bsg(v) != 'v' ||
	    !rt_view_context_coord_set_bsg(v, 'm') ||
	    rt_view_coord_from_bsg(v) != 'm') {
	printf("FAIL: BSG coord set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_coord_set(v, 'v') ||
	    rt_view_context_coord_get(v) != 'v') {
	printf("FAIL: retained context coord set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_snap_lines_from_bsg(NULL) ||
	    rt_view_context_snap_lines_from_bsg(NULL) ||
	    rt_view_snap_source_flags_from_bsg(NULL) ||
	    rt_view_snap_kind_mask_from_bsg(NULL) ||
	    rt_view_context_snap_kind_mask_from_bsg(NULL) ||
	    rt_view_snap_lines_set_bsg(NULL, 1) ||
	    rt_view_context_snap_lines_set_bsg(NULL, 1) ||
	    rt_view_snap_source_flags_set_bsg(NULL, RT_VIEW_SNAP_TCL_BSG) ||
	    rt_view_context_snap_source_flags_set_bsg(NULL, RT_VIEW_SNAP_TCL_BSG) ||
	    rt_view_context_snap_lines_get(NULL) ||
	    rt_view_context_snap_kind_mask_get(NULL) ||
	    rt_view_context_snap_source_flags_get(NULL) ||
	    rt_view_context_snap_lines_set(NULL, 1) ||
	    rt_view_context_snap_source_flags_set(NULL, RT_VIEW_SNAP_TCL)) {
	printf("FAIL: null BSG snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_lines_set_bsg(v, 1) ||
	    !rt_view_snap_source_flags_set_bsg(v, RT_VIEW_SNAP_TCL_BSG) ||
	    !rt_view_snap_lines_from_bsg(v) ||
	    !rt_view_context_snap_lines_from_bsg(v) ||
	    !rt_view_context_snap_lines_set_bsg(v, 0) ||
	    rt_view_snap_lines_from_bsg(v) ||
	    !rt_view_context_snap_lines_set_bsg(v, 1) ||
	    !rt_view_snap_lines_from_bsg(v) ||
	    rt_view_snap_source_flags_from_bsg(v) != RT_VIEW_SNAP_TCL_BSG ||
	    !(rt_view_snap_kind_mask_from_bsg(v) & RT_VIEW_SNAP_KIND_ENDPOINT_BSG) ||
	    !(rt_view_context_snap_kind_mask_from_bsg(v) & RT_VIEW_SNAP_KIND_ENDPOINT_BSG)) {
	printf("FAIL: BSG snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_snap_lines_set(v, 0) ||
	    rt_view_context_snap_lines_get(v) ||
	    !rt_view_context_snap_lines_set(v, 1) ||
	    !rt_view_context_snap_lines_get(v) ||
	    !rt_view_context_snap_source_flags_set(v, RT_VIEW_SNAP_TCL) ||
	    rt_view_context_snap_source_flags_get(v) != RT_VIEW_SNAP_TCL ||
	    !(rt_view_context_snap_kind_mask_get(v) & RT_VIEW_SNAP_KIND_ENDPOINT)) {
	printf("FAIL: neutral snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }

    point_t polygon_point = VINIT_ZERO;
    if (!rt_view_polygon_ref_is_null_bsg(rt_view_polygon_create_bsg(NULL,
		BSG_POLYGON_GENERAL, &polygon_point)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_create_bsg(v,
		BSG_POLYGON_GENERAL, NULL)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_select_bsg(NULL,
		&polygon_point)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_select_bsg(v, NULL)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_find_bsg(NULL,
		"snap_count_a")) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_find_bsg(v, NULL)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_dup_bsg(NULL,
		"snap_count_a", "snap_count_c")) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_dup_bsg(v,
		NULL, "snap_count_c")) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_polygon_dup_bsg(v,
		"snap_count_a", NULL))) {
	printf("FAIL: null BSG polygon view adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_create_bsg(NULL,
		BSG_POLYGON_GENERAL, &polygon_point)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_create_bsg(v,
		BSG_POLYGON_GENERAL, NULL)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_select_bsg(NULL,
		&polygon_point)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_select_bsg(v,
		NULL)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_find_bsg(NULL,
		"snap_count_a")) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_find_bsg(v,
		NULL)) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_dup_bsg(NULL,
		"snap_count_a", "snap_count_c")) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_dup_bsg(v,
		NULL, "snap_count_c")) ||
	    !rt_view_polygon_ref_is_null_bsg(rt_view_context_polygon_dup_bsg(v,
		"snap_count_a", NULL))) {
	printf("FAIL: null BSG polygon context adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_polygon_ref_is_null(rt_view_context_polygon_create(NULL,
		RT_VIEW_POLYGON_GENERAL, &polygon_point)) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_create(v,
		RT_VIEW_POLYGON_GENERAL, NULL)) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_select(NULL,
		&polygon_point)) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_select(v,
		NULL)) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_find(NULL,
		"snap_count_a")) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_find(v,
		NULL)) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_dup(NULL,
		"snap_count_a", "snap_count_c")) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_dup(v,
		NULL, "snap_count_c")) ||
	    !rt_view_polygon_ref_is_null(rt_view_context_polygon_dup(v,
		"snap_count_a", NULL))) {
	printf("FAIL: null neutral polygon context adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_polygon_ref poly_a = rt_view_polygon_create_bsg(v,
	    BSG_POLYGON_GENERAL, &polygon_point);
    rt_view_polygon_ref poly_b = rt_view_context_polygon_create_bsg(v,
	    BSG_POLYGON_GENERAL, &polygon_point);
    if (rt_view_polygon_ref_is_null_bsg(poly_a) ||
	    rt_view_polygon_ref_is_null_bsg(poly_b) ||
	    !rt_view_polygon_set_name_bsg(poly_a, "snap_count_a") ||
	    !rt_view_polygon_set_name_bsg(poly_b, "snap_count_b") ||
	    rt_view_polygon_snap_count_bsg(NULL, poly_a) ||
	    rt_view_context_polygon_snap_count_bsg(NULL, poly_a) ||
	    rt_view_context_polygon_snap_count(NULL, poly_a) ||
	    rt_view_polygon_snap_count_bsg(v, poly_a) != 1 ||
	    rt_view_context_polygon_snap_count_bsg(v, poly_b) != 1 ||
	    rt_view_context_polygon_snap_exclude_set_bsg(NULL, poly_a) ||
	    !rt_view_context_polygon_snap_exclude_set_bsg(v, poly_a) ||
	    rt_view_context_polygon_snap_exclude_set(NULL, poly_a) ||
	    !rt_view_context_snap_exclude_feature_clear_bsg(v)) {
	printf("FAIL: BSG polygon snap-count adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_polygon_ref poly_neutral = rt_view_context_polygon_create(v,
	    RT_VIEW_POLYGON_GENERAL, &polygon_point);
    if (rt_view_polygon_ref_is_null(poly_neutral) ||
	    !rt_view_polygon_set_name(poly_neutral, "snap_count_neutral") ||
	    rt_view_context_polygon_snap_count(v, poly_neutral) != 2 ||
	    !rt_view_context_polygon_snap_exclude_set(v, poly_neutral) ||
	    !rt_view_context_snap_exclude_feature_clear(v)) {
	printf("FAIL: neutral polygon snap-count adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_polygon_ref found_poly =
	rt_view_context_polygon_find_bsg(v, "snap_count_a");
    rt_view_polygon_ref dup_poly = rt_view_context_polygon_dup_bsg(v,
	    "snap_count_a", "snap_count_c");
    rt_view_polygon_ref selected_poly = rt_view_context_polygon_select_bsg(v,
	    &polygon_point);
    rt_view_polygon_ref found_neutral =
	rt_view_context_polygon_find(v, "snap_count_neutral");
    rt_view_polygon_ref dup_neutral = rt_view_context_polygon_dup(v,
	    "snap_count_neutral", "snap_count_d");
    rt_view_polygon_ref selected_neutral = rt_view_context_polygon_select(v,
	    &polygon_point);
    struct test_polygon_record_callback_state poly_visit_state = {0, {0, 0}};
    rt_view_polygon_visit_records_bsg(NULL, test_polygon_record_callback,
	    &poly_visit_state);
    rt_view_polygon_visit_records_bsg(v, NULL, &poly_visit_state);
    rt_view_context_polygon_visit_records_bsg(NULL, test_polygon_record_callback,
	    &poly_visit_state);
    rt_view_context_polygon_visit_records_bsg(v, NULL, &poly_visit_state);
    rt_view_context_polygon_visit_records_bsg(v, test_polygon_record_callback,
	    &poly_visit_state);
    rt_view_context_polygon_visit_records(NULL, test_polygon_record_callback,
	    &poly_visit_state);
    rt_view_context_polygon_visit_records(v, NULL, &poly_visit_state);
    rt_view_context_polygon_visit_records(v, test_polygon_record_callback,
	    &poly_visit_state);
    if (found_poly.token != poly_a.token ||
	    rt_view_polygon_ref_is_null_bsg(dup_poly) ||
	    rt_view_polygon_ref_is_null_bsg(selected_poly) ||
	    found_neutral.token != poly_neutral.token ||
	    rt_view_polygon_ref_is_null(dup_neutral) ||
	    rt_view_polygon_ref_is_null(selected_neutral) ||
	    poly_visit_state.count < 6) {
	printf("FAIL: BSG polygon view adapter\n");
	ret = 1;
	goto cleanup;
    }

    struct rt_view_polygon_record poly_rec;
    if (!bsg_polygon_set_current(test_polygon_ref_to_bsg(poly_a), 0, 0) ||
	    rt_view_polygon_clear_point_selection_bsg(NULL) ||
	    rt_view_context_polygon_clear_point_selection_bsg(NULL) ||
	    rt_view_context_polygon_clear_point_selection(NULL) ||
	    !rt_view_context_polygon_clear_point_selection_bsg(v) ||
	    !rt_view_polygon_record_get_bsg(poly_a, &poly_rec) ||
	    poly_rec.curr_point_i != -1) {
	printf("FAIL: BSG polygon point-selection clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!bsg_polygon_set_current(test_polygon_ref_to_bsg(poly_a), 0, 0) ||
	    !rt_view_context_polygon_clear_point_selection(v) ||
	    !rt_view_polygon_record_get(poly_a, &poly_rec) ||
	    poly_rec.curr_point_i != -1) {
	printf("FAIL: neutral polygon point-selection clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    struct bu_color edge_color = BU_COLOR_YELLOW;
    struct bu_color fill_color = BU_COLOR_BLUE;
    if (rt_view_polygon_record_get_bsg(poly_a, NULL) ||
	    rt_view_polygon_record_get(poly_a, NULL) ||
	    rt_view_polygon_set_name_bsg(poly_a, NULL) ||
	    rt_view_polygon_set_name(poly_a, NULL) ||
	    rt_view_polygon_move_bsg(poly_a, NULL, &polygon_point) ||
	    rt_view_polygon_move_bsg(poly_a, &polygon_point, NULL) ||
	    rt_view_polygon_move(poly_a, NULL, &polygon_point) ||
	    rt_view_polygon_move(poly_a, &polygon_point, NULL) ||
	    rt_view_polygon_user_data_bsg(RT_VIEW_POLYGON_REF_NULL) ||
	    rt_view_polygon_user_data(RT_VIEW_POLYGON_REF_NULL) ||
	    !rt_view_polygon_set_view_bsg(poly_a, v) ||
	    !rt_view_polygon_set_context_bsg(poly_a, v) ||
	    !rt_view_polygon_set_context(poly_a, v) ||
	    !rt_view_polygon_user_data_set_bsg(poly_a, v) ||
	    rt_view_polygon_user_data_bsg(poly_a) != v ||
	    !rt_view_polygon_user_data_set(poly_a, &ret) ||
	    rt_view_polygon_user_data(poly_a) != &ret) {
	printf("FAIL: BSG polygon object adapter argument handling\n");
	ret = 1;
	goto cleanup;
    }
    (void)rt_view_polygon_set_visual_bsg(poly_a, &edge_color, &fill_color,
	    1.0, 0.0, 5.0, 0.25, 1);
    (void)rt_view_polygon_set_open_bsg(poly_a, 1);
    (void)rt_view_polygon_clear_selected_point_bsg(poly_a);
    (void)rt_view_polygon_update_bsg(poly_a, v, BSG_POLYGON_UPDATE_DEFAULT);
    (void)rt_view_polygon_update_context_bsg(poly_a, v,
	    BSG_POLYGON_UPDATE_DEFAULT);
    (void)rt_view_polygon_update_screen_pt_context_bsg(poly_a, v, 10, 20,
	    BSG_POLYGON_UPDATE_PT_SELECT);
    (void)rt_view_polygon_set_visual(poly_a, &edge_color, &fill_color,
	    0.0, 1.0, 6.0, 0.5, 1);
    (void)rt_view_polygon_set_open(poly_a, 1);
    (void)rt_view_polygon_clear_selected_point(poly_a);
    (void)rt_view_polygon_update_context(poly_a, v, BSG_POLYGON_UPDATE_DEFAULT);
    (void)rt_view_polygon_update_screen_pt_context(poly_a, v, 10, 20,
	    BSG_POLYGON_UPDATE_PT_SELECT);
    if (!rt_view_polygon_record_get(poly_a, &poly_rec) ||
	    poly_rec.fill_flag != 1 ||
	    poly_rec.user_data != &ret) {
	printf("FAIL: BSG polygon object adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_polygon_ref bool_poly = rt_view_polygon_create_bsg(v,
	    BSG_POLYGON_GENERAL, &polygon_point);
    rt_view_polygon_ref bool_poly_neutral = rt_view_context_polygon_create(v,
	    RT_VIEW_POLYGON_GENERAL, &polygon_point);
    if (rt_view_polygon_ref_is_null_bsg(bool_poly) ||
	    rt_view_polygon_ref_is_null(bool_poly_neutral) ||
	    rt_view_polygon_csg_bsg(poly_a, bool_poly, bg_Union) < 0 ||
	    rt_view_polygon_csg(poly_a, bool_poly_neutral, bg_Union) < 0 ||
	    !rt_view_polygon_remove_bsg(bool_poly) ||
	    !rt_view_polygon_remove(bool_poly_neutral) ||
	    !rt_view_polygon_remove(poly_neutral)) {
	printf("FAIL: BSG polygon object CSG/remove adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_polygon_ref_is_null_bsg(
		rt_view_polygon_import_sketch_bsg(NULL, NULL, NULL, v)) ||
	    !rt_view_polygon_ref_is_null_bsg(
		rt_view_polygon_import_sketch_context_bsg(NULL, NULL, NULL, v)) ||
	    !rt_view_polygon_ref_is_null(
		rt_view_polygon_import_sketch_context(NULL, NULL, NULL, v)) ||
	    rt_view_polygon_export_sketch_bsg(NULL, "sketch", poly_a) ||
	    rt_view_polygon_export_sketch(NULL, "sketch", poly_a)) {
	printf("FAIL: BSG polygon sketch adapter argument handling\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_source_flags_set_bsg(v, RT_VIEW_SNAP_DB_BSG) ||
	    !(rt_view_prepare_tcl_snap_bsg(v) & RT_VIEW_SNAP_KIND_ENDPOINT_BSG) ||
	    !(rt_view_context_prepare_tcl_snap_bsg(v) & RT_VIEW_SNAP_KIND_ENDPOINT_BSG) ||
	    rt_view_snap_source_flags_from_bsg(v) != RT_VIEW_SNAP_TCL_BSG ||
	    rt_view_prepare_tcl_snap_bsg(NULL) ||
	    rt_view_context_prepare_tcl_snap_bsg(NULL) ||
	    rt_view_context_prepare_tcl_snap(NULL)) {
	printf("FAIL: BSG Tcl snap preparation adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!(rt_view_context_prepare_tcl_snap(v) & RT_VIEW_SNAP_KIND_ENDPOINT) ||
	    rt_view_context_snap_source_flags_get(v) != RT_VIEW_SNAP_TCL) {
	printf("FAIL: neutral Tcl snap preparation adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_center_linesnap_bsg(NULL) ||
	    rt_view_context_center_linesnap_bsg(NULL) ||
	    rt_view_context_center_linesnap(NULL) ||
	    !rt_view_center_linesnap_bsg(v) ||
	    !rt_view_context_center_linesnap_bsg(v) ||
	    !rt_view_context_center_linesnap(v)) {
	printf("FAIL: BSG line-snap recenter adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_zclip_from_bsg(NULL) ||
	    rt_view_context_zclip_from_bsg(NULL) ||
	    rt_view_context_zclip_get(NULL) ||
	    rt_view_zclip_set_bsg(NULL, 1) ||
	    rt_view_context_zclip_set_bsg(NULL, 1) ||
	    rt_view_context_zclip_set(NULL, 1)) {
	printf("FAIL: null BSG zclip adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_zclip_set_bsg(v, 1) ||
	    rt_view_zclip_from_bsg(v) != 1 ||
	    rt_view_context_zclip_from_bsg(v) != 1 ||
	    rt_view_context_zclip_get(v) != 1 ||
	    !rt_view_context_zclip_set_bsg(v, 0) ||
	    rt_view_zclip_from_bsg(v)) {
	printf("FAIL: BSG zclip adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_zclip_set(v, 1) ||
	    rt_view_context_zclip_get(v) != 1 ||
	    !rt_view_context_zclip_set(v, 0) ||
	    rt_view_context_zclip_get(v)) {
	printf("FAIL: retained context zclip adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_framebuffer_mode_from_bsg(NULL) ||
	    rt_view_context_framebuffer_mode_from_bsg(NULL) ||
	    rt_view_context_framebuffer_mode_get(NULL) ||
	    rt_view_framebuffer_mode_set_bsg(NULL, 2) ||
	    rt_view_context_framebuffer_mode_set_bsg(NULL, 2) ||
	    rt_view_context_framebuffer_mode_set(NULL, 2)) {
	printf("FAIL: null BSG framebuffer mode adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_framebuffer_mode_set_bsg(v, 2) ||
	    rt_view_framebuffer_mode_from_bsg(v) != 2 ||
	    rt_view_context_framebuffer_mode_from_bsg(v) != 2 ||
	    rt_view_context_framebuffer_mode_get(v) != 2 ||
	    !rt_view_context_framebuffer_mode_set(v, 2) ||
	    rt_view_context_framebuffer_mode_get(v) != 2 ||
	    !rt_view_framebuffer_mode_set_bsg(v, 1) ||
	    rt_view_context_framebuffer_mode_from_bsg(v) != 1 ||
	    rt_view_context_framebuffer_mode_get(v) != 1 ||
	    !rt_view_context_framebuffer_mode_set(v, 0) ||
	    rt_view_framebuffer_mode_from_bsg(v) ||
	    rt_view_context_framebuffer_mode_get(v)) {
	printf("FAIL: BSG framebuffer mode adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_cleared_from_bsg(NULL) ||
	    rt_view_context_cleared_from_bsg(NULL) ||
	    rt_view_context_cleared_get(NULL) ||
	    rt_view_cleared_set_bsg(NULL, 1) ||
	    rt_view_context_cleared_set_bsg(NULL, 1) ||
	    rt_view_context_cleared_set(NULL, 1)) {
	printf("FAIL: null BSG cleared-state adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_cleared_set_bsg(v, 1) ||
	    rt_view_cleared_from_bsg(v) != 1 ||
	    rt_view_context_cleared_from_bsg(v) != 1 ||
	    rt_view_context_cleared_get(v) != 1 ||
	    !rt_view_cleared_set_bsg(v, 2) ||
	    rt_view_cleared_from_bsg(v) != 1 ||
	    !rt_view_context_cleared_set(v, 0) ||
	    rt_view_cleared_from_bsg(v)) {
	printf("FAIL: BSG cleared-state adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_settings_shared_bsg(NULL, v) ||
	    rt_view_context_settings_shared_bsg(NULL, v) ||
	    rt_view_context_settings_shared(NULL, v) ||
	    rt_view_settings_shared_bsg(v, NULL) ||
	    rt_view_context_settings_shared_bsg(v, NULL) ||
	    rt_view_context_settings_shared(v, NULL) ||
	    rt_view_settings_shared_bsg(v, shared) ||
	    rt_view_context_settings_shared_bsg(v, shared) ||
	    rt_view_context_settings_shared(v, shared)) {
	printf("FAIL: distinct BSG shared-settings adapter\n");
	ret = 1;
	goto cleanup;
    }
    shared->settings_active = v->settings_active;
    bsg_view_state_sync_from_view(shared->gv_state, shared);
    if (!rt_view_settings_shared_bsg(v, shared) ||
	    !rt_view_context_settings_shared_bsg(v, shared) ||
	    !rt_view_context_settings_shared(v, shared)) {
	printf("FAIL: shared BSG shared-settings adapter\n");
	ret = 1;
	goto cleanup;
    }
    shared->settings_active = shared->settings_local;
    bsg_view_state_sync_from_view(shared->gv_state, shared);

    if (!fastf_equal(rt_view_snap_tolerance_factor_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_snap_tolerance_factor_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_snap_tolerance_factor_get(NULL), 1.0) ||
	    rt_view_snap_tolerance_factor_set_bsg(NULL, 2.0) ||
	    rt_view_context_snap_tolerance_factor_set_bsg(NULL, 2.0) ||
	    rt_view_context_snap_tolerance_factor_set(NULL, 2.0)) {
	printf("FAIL: null BSG snap tolerance factor adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_tolerance_factor_set_bsg(v, 2.5) ||
	    !fastf_equal(rt_view_snap_tolerance_factor_from_bsg(v), 2.5)) {
	printf("FAIL: BSG snap tolerance factor adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_snap_tolerance_factor_set_bsg(v, 3.5) ||
	    !fastf_equal(rt_view_context_snap_tolerance_factor_from_bsg(v), 3.5)) {
	printf("FAIL: BSG context snap tolerance factor adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_snap_tolerance_factor_set(v, 4.5) ||
	    !fastf_equal(rt_view_context_snap_tolerance_factor_get(v), 4.5)) {
	printf("FAIL: neutral snap tolerance factor adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_width_from_bsg(NULL) != 0 ||
	    rt_view_height_from_bsg(NULL) != 0) {
	printf("FAIL: null BSG raw width/height adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_width_from_bsg(NULL) != 0 ||
	    rt_view_context_height_from_bsg(NULL) != 0) {
	printf("FAIL: null BSG context width/height adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_width_get(NULL) != 0 ||
	    rt_view_context_height_get(NULL) != 0) {
	printf("FAIL: null retained context width/height adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_dimensions_set_bsg(NULL, 123, 456)) {
	printf("FAIL: null BSG raw dimensions set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_dimensions_set_bsg(NULL, 123, 456)) {
	printf("FAIL: null BSG context dimensions set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_dimensions_set(NULL, 123, 456)) {
	printf("FAIL: null retained context dimensions set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_dimensions_set_bsg(v, 0, 0)) {
	printf("FAIL: BSG raw dimensions zero set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_width_from_bsg(v) != 0 ||
	    rt_view_height_from_bsg(v) != 0) {
	printf("FAIL: BSG raw width/height zero-preserving adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_dimensions_set_bsg(v, 123, 456)) {
	printf("FAIL: BSG raw dimensions set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_width_from_bsg(v) != 123 ||
	    rt_view_height_from_bsg(v) != 456) {
	printf("FAIL: BSG raw width/height adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_width_from_bsg(v) != 123 ||
	    rt_view_context_height_from_bsg(v) != 456) {
	printf("FAIL: BSG context width/height adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_width_get(v) != 123 ||
	    rt_view_context_height_get(v) != 456) {
	printf("FAIL: retained context width/height adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_context_dimensions_set_bsg(v, 321, 654) ||
	    rt_view_width_from_bsg(v) != 321 ||
	    rt_view_height_from_bsg(v) != 654) {
	printf("FAIL: BSG context dimensions set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_dimensions_set(v, 111, 222) ||
	    rt_view_context_width_get(v) != 111 ||
	    rt_view_context_height_get(v) != 222) {
	printf("FAIL: retained context dimensions set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!fastf_equal(rt_view_radius_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_radius_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_context_radius_get(NULL), 1.0)) {
	printf("FAIL: null BSG view radius adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->radius = 77.0;
    if (!fastf_equal(rt_view_radius_from_bsg(v), 77.0) ||
	    !fastf_equal(rt_view_context_radius_from_bsg(v), 77.0) ||
	    !fastf_equal(rt_view_context_radius_get(v), 77.0)) {
	printf("FAIL: BSG view radius adapter\n");
	ret = 1;
	goto cleanup;
    }

    fastf_t direct_vx = -1.0;
    fastf_t direct_vy = -1.0;
    fastf_t adapter_vx = -1.0;
    fastf_t adapter_vy = -1.0;
    if (rt_view_screen_to_view_from_bsg(&adapter_vx, &adapter_vy,
	    NULL, 60.0, 70.0) ||
	    !fastf_equal(adapter_vx, 0.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null BSG screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    adapter_vx = -1.0;
    adapter_vy = -1.0;
    if (rt_view_context_screen_to_view_from_bsg(&adapter_vx, &adapter_vy,
	    NULL, 60.0, 70.0) ||
	    !fastf_equal(adapter_vx, 0.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null BSG context screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    adapter_vx = -1.0;
    adapter_vy = -1.0;
    if (rt_view_context_screen_to_view(&adapter_vx, &adapter_vy,
	    NULL, 60.0, 70.0) ||
	    !fastf_equal(adapter_vx, 0.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null retained context screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_screen_to_view_from_bsg(NULL, &adapter_vy,
	    v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null output BSG screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    adapter_vy = -1.0;
    if (rt_view_context_screen_to_view_from_bsg(NULL, &adapter_vy,
	    v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null output BSG context screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    adapter_vy = -1.0;
    if (rt_view_context_screen_to_view(NULL, &adapter_vy,
	    v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null output retained context screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (bsg_screen_to_view(v, &direct_vx, &direct_vy, 60.0, 70.0) ||
	    !rt_view_screen_to_view_from_bsg(&adapter_vx, &adapter_vy,
		v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vx, direct_vx) ||
	    !fastf_equal(adapter_vy, direct_vy)) {
	printf("FAIL: BSG screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    adapter_vx = -1.0;
    adapter_vy = -1.0;
    if (!rt_view_context_screen_to_view_from_bsg(&adapter_vx, &adapter_vy,
	    v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vx, direct_vx) ||
	    !fastf_equal(adapter_vy, direct_vy)) {
	printf("FAIL: BSG context screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }
    adapter_vx = -1.0;
    adapter_vy = -1.0;
    if (!rt_view_context_screen_to_view(&adapter_vx, &adapter_vy,
	    v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vx, direct_vx) ||
	    !fastf_equal(adapter_vy, direct_vy)) {
	printf("FAIL: retained context screen-to-view adapter\n");
	ret = 1;
	goto cleanup;
    }

    point_t direct_point = VINIT_ZERO;
    point_t adapter_point = VINIT_ZERO;
    if (rt_view_screen_point_from_bsg(adapter_point, NULL, 60.0, 70.0) ||
	    rt_view_context_screen_point_from_bsg(adapter_point, NULL, 60.0, 70.0) ||
	    rt_view_context_screen_point_get(adapter_point, NULL, 60.0, 70.0) ||
	    !fastf_equal(adapter_point[X], 0.0) ||
	    !fastf_equal(adapter_point[Y], 0.0) ||
	    !fastf_equal(adapter_point[Z], 0.0)) {
	printf("FAIL: null BSG screen-point adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSETALL(adapter_point, -1.0);
    if (bsg_screen_pt(&direct_point, 60.0, 70.0, v) ||
	    !rt_view_screen_point_from_bsg(adapter_point, v, 60.0, 70.0) ||
	    !VNEAR_EQUAL(direct_point, adapter_point, BN_TOL_DIST)) {
	printf("FAIL: BSG screen-point adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSETALL(adapter_point, -1.0);
    if (!rt_view_context_screen_point_from_bsg(adapter_point, v, 60.0, 70.0) ||
	    !VNEAR_EQUAL(direct_point, adapter_point, BN_TOL_DIST)) {
	printf("FAIL: BSG context screen-point adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSETALL(adapter_point, -1.0);
    if (!rt_view_context_screen_point_get(adapter_point, v, 60.0, 70.0) ||
	    !VNEAR_EQUAL(direct_point, adapter_point, BN_TOL_DIST)) {
	printf("FAIL: retained context screen-point adapter\n");
	ret = 1;
	goto cleanup;
    }

    point_t current_point = VINIT_ZERO;
    if (rt_view_current_point_from_bsg(current_point, NULL) ||
	    rt_view_context_current_point_from_bsg(current_point, NULL) ||
	    rt_view_context_current_point_get(current_point, NULL) ||
	    !fastf_equal(current_point[X], 0.0) ||
	    !fastf_equal(current_point[Y], 0.0) ||
	    !fastf_equal(current_point[Z], 0.0) ||
	    rt_view_current_point_from_bsg(NULL, v) ||
	    rt_view_context_current_point_from_bsg(NULL, v) ||
	    rt_view_context_current_point_get(NULL, v)) {
	printf("FAIL: null BSG current-point adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(v->gv_point, 1.0, 2.0, 3.0);
    if (!rt_view_current_point_from_bsg(current_point, v) ||
	    !rt_view_context_current_point_from_bsg(current_point, v) ||
	    !rt_view_context_current_point_get(current_point, v) ||
	    !fastf_equal(current_point[X], 1.0) ||
	    !fastf_equal(current_point[Y], 2.0) ||
	    !fastf_equal(current_point[Z], 3.0)) {
	printf("FAIL: BSG current-point get adapter\n");
	ret = 1;
	goto cleanup;
    }
    point_t updated_point;
    VSET(updated_point, 4.0, 5.0, 6.0);
    if (rt_view_current_point_set_bsg(NULL, updated_point) ||
	    rt_view_current_point_set_bsg(v, NULL) ||
	    rt_view_context_current_point_set_bsg(NULL, updated_point) ||
	    rt_view_context_current_point_set_bsg(v, NULL) ||
	    rt_view_context_current_point_set(NULL, updated_point) ||
	    rt_view_context_current_point_set(v, NULL) ||
	    !rt_view_current_point_set_bsg(v, updated_point) ||
	    !VNEAR_EQUAL(v->gv_point, updated_point, BN_TOL_DIST)) {
	printf("FAIL: BSG current-point set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(updated_point, 7.0, 8.0, 9.0);
    if (!rt_view_context_current_point_set_bsg(v, updated_point) ||
	    !VNEAR_EQUAL(v->gv_point, updated_point, BN_TOL_DIST)) {
	printf("FAIL: BSG context current-point set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(updated_point, 10.0, 11.0, 12.0);
    if (!rt_view_context_current_point_set(v, updated_point) ||
	    !rt_view_context_current_point_get(current_point, v) ||
	    !VNEAR_EQUAL(current_point, updated_point, BN_TOL_DIST)) {
	printf("FAIL: retained context current-point set adapter\n");
	ret = 1;
	goto cleanup;
    }

    fastf_t previous_x = -1.0;
    fastf_t previous_y = -1.0;
    if (rt_view_previous_mouse_from_bsg(&previous_x, &previous_y, NULL) ||
	    !fastf_equal(previous_x, 0.0) || !fastf_equal(previous_y, 0.0) ||
	    rt_view_previous_mouse_from_bsg(NULL, &previous_y, v) ||
	    rt_view_previous_mouse_from_bsg(&previous_x, NULL, v)) {
	printf("FAIL: null BSG previous-mouse adapter\n");
	ret = 1;
	goto cleanup;
    }
    previous_x = -1.0;
    previous_y = -1.0;
    if (rt_view_context_previous_mouse_from_bsg(&previous_x, &previous_y, NULL) ||
	    !fastf_equal(previous_x, 0.0) || !fastf_equal(previous_y, 0.0) ||
	    rt_view_context_previous_mouse_from_bsg(NULL, &previous_y, v) ||
	    rt_view_context_previous_mouse_from_bsg(&previous_x, NULL, v) ||
	    rt_view_context_previous_mouse_get(NULL, &previous_y, v) ||
	    rt_view_context_previous_mouse_get(&previous_x, NULL, v)) {
	printf("FAIL: null BSG context previous-mouse adapter\n");
	ret = 1;
	goto cleanup;
    }
    previous_x = -1.0;
    previous_y = -1.0;
    if (rt_view_context_previous_mouse_get(&previous_x, &previous_y, NULL) ||
	    !fastf_equal(previous_x, 0.0) || !fastf_equal(previous_y, 0.0)) {
	printf("FAIL: null retained context previous-mouse adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_previous_mouse_set_bsg(NULL, 5.5, 6.25) ||
	    !rt_view_previous_mouse_set_bsg(v, 5.5, 6.25) ||
	    !rt_view_previous_mouse_from_bsg(&previous_x, &previous_y, v) ||
	    !fastf_equal(previous_x, 5.5) ||
	    !fastf_equal(previous_y, 6.25) ||
	    !fastf_equal(v->gv_prevMouseX, 5.5) ||
	    !fastf_equal(v->gv_prevMouseY, 6.25)) {
	printf("FAIL: BSG previous-mouse get/set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_previous_mouse_set_bsg(NULL, 7.5, 8.25) ||
	    !rt_view_context_previous_mouse_set_bsg(v, 7.5, 8.25) ||
	    !rt_view_context_previous_mouse_from_bsg(&previous_x, &previous_y, v) ||
	    !fastf_equal(previous_x, 7.5) ||
	    !fastf_equal(previous_y, 8.25) ||
	    !fastf_equal(v->gv_prevMouseX, 7.5) ||
	    !fastf_equal(v->gv_prevMouseY, 8.25)) {
	printf("FAIL: BSG context previous-mouse get/set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_context_previous_mouse_set(NULL, 9.5, 10.25) ||
	    !rt_view_context_previous_mouse_set(v, 9.5, 10.25) ||
	    !rt_view_context_previous_mouse_get(&previous_x, &previous_y, v) ||
	    !fastf_equal(previous_x, 9.5) ||
	    !fastf_equal(previous_y, 10.25) ||
	    !fastf_equal(v->gv_prevMouseX, 9.5) ||
	    !fastf_equal(v->gv_prevMouseY, 10.25)) {
	printf("FAIL: retained context previous-mouse get/set adapter\n");
	ret = 1;
	goto cleanup;
    }

    struct rt_view_mouse_delta_settings mouse_settings = RT_VIEW_MOUSE_DELTA_SETTINGS_INIT;
    if (rt_view_mouse_delta_settings_from_bsg(NULL, v) ||
	    rt_view_context_mouse_delta_settings_from_bsg(NULL, v) ||
	    rt_view_context_mouse_delta_settings_get(NULL, v) ||
	    rt_view_mouse_delta_settings_from_bsg(&mouse_settings, NULL) ||
	    rt_view_context_mouse_delta_settings_from_bsg(&mouse_settings, NULL) ||
	    rt_view_context_mouse_delta_settings_get(&mouse_settings, NULL) ||
	    !fastf_equal(mouse_settings.min_delta, 0.0) ||
	    !fastf_equal(mouse_settings.max_delta, 0.0) ||
	    !fastf_equal(mouse_settings.rotate_scale, 0.0) ||
	    !fastf_equal(mouse_settings.scale_scale, 0.0)) {
	printf("FAIL: null BSG mouse-delta settings adapter\n");
	ret = 1;
	goto cleanup;
    }
    v->gv_minMouseDelta = -3.5;
    v->gv_maxMouseDelta = 9.25;
    v->gv_rscale = 0.125;
    v->gv_sscale = 2.5;
    if (!rt_view_mouse_delta_settings_from_bsg(&mouse_settings, v) ||
	    !rt_view_context_mouse_delta_settings_from_bsg(&mouse_settings, v) ||
	    !rt_view_context_mouse_delta_settings_get(&mouse_settings, v) ||
	    !fastf_equal(mouse_settings.min_delta, -3.5) ||
	    !fastf_equal(mouse_settings.max_delta, 9.25) ||
	    !fastf_equal(mouse_settings.rotate_scale, 0.125) ||
	    !fastf_equal(mouse_settings.scale_scale, 2.5)) {
	printf("FAIL: BSG mouse-delta settings adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_mouse_x = 11;
    v->gv_mouse_y = 13;
    VSETALL(v->gv_point, -1.0);
    VSETALL(direct_point, 0.0);
    if (bsg_screen_pt(&direct_point, 75.0, 85.0, v) ||
	    rt_view_mouse_state_set_bsg(NULL, 75, 85) ||
	    rt_view_context_mouse_state_set_bsg(NULL, 75, 85) ||
	    rt_view_context_mouse_state_set(NULL, 75, 85) ||
	    !rt_view_mouse_state_set_bsg(v, 75, 85) ||
	    !fastf_equal(v->gv_prevMouseX, 11.0) ||
	    !fastf_equal(v->gv_prevMouseY, 13.0) ||
	    v->gv_mouse_x != 75 ||
	    v->gv_mouse_y != 85 ||
	    !VNEAR_EQUAL(v->gv_point, direct_point, BN_TOL_DIST)) {
	printf("FAIL: BSG mouse/current-point state adapter\n");
	ret = 1;
	goto cleanup;
    }
    v->gv_mouse_x = 23;
    v->gv_mouse_y = 29;
    VSETALL(v->gv_point, -1.0);
    VSETALL(direct_point, 0.0);
    if (bsg_screen_pt(&direct_point, 35.0, 45.0, v) ||
	    !rt_view_context_mouse_state_set(v, 35, 45) ||
	    !fastf_equal(v->gv_prevMouseX, 23.0) ||
	    !fastf_equal(v->gv_prevMouseY, 29.0) ||
	    v->gv_mouse_x != 35 ||
	    v->gv_mouse_y != 45 ||
	    !VNEAR_EQUAL(v->gv_point, direct_point, BN_TOL_DIST)) {
	printf("FAIL: retained context mouse/current-point state adapter\n");
	ret = 1;
	goto cleanup;
    }
    v->gv_mouse_x = 17;
    v->gv_mouse_y = 19;
    VSETALL(v->gv_point, -1.0);
    VSETALL(direct_point, 0.0);
    if (bsg_screen_pt(&direct_point, 65.0, 95.0, v) ||
	    !rt_view_context_mouse_state_set_bsg(v, 65, 95) ||
	    !fastf_equal(v->gv_prevMouseX, 17.0) ||
	    !fastf_equal(v->gv_prevMouseY, 19.0) ||
	    v->gv_mouse_x != 65 ||
	    v->gv_mouse_y != 95 ||
	    !VNEAR_EQUAL(v->gv_point, direct_point, BN_TOL_DIST)) {
	printf("FAIL: BSG context mouse/current-point state adapter\n");
	ret = 1;
	goto cleanup;
    }

    struct bu_vls direct_name = BU_VLS_INIT_ZERO;
    struct bu_vls adapter_name = BU_VLS_INIT_ZERO;
    if (rt_view_unique_object_name_bsg(NULL, "view_poly_0", v) ||
	    rt_view_unique_object_name_bsg(&adapter_name, "view_poly_0", NULL) ||
	    rt_view_context_unique_object_name_bsg(NULL, "view_poly_0", v) ||
	    rt_view_context_unique_object_name_bsg(&adapter_name,
		"view_poly_0", NULL) ||
	    rt_view_context_unique_object_name(NULL, "view_poly_0", v) ||
	    rt_view_context_unique_object_name(&adapter_name,
		"view_poly_0", NULL)) {
	printf("FAIL: null BSG unique-object-name adapter\n");
	ret = 1;
	goto cleanup;
    }
    bsg_uniq_obj_name(&direct_name, "view_poly_0", v);
    if (!rt_view_unique_object_name_bsg(&adapter_name, "view_poly_0", v) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&direct_name), bu_vls_cstr(&adapter_name))) {
	printf("FAIL: BSG unique-object-name adapter\n");
	bu_vls_free(&direct_name);
	bu_vls_free(&adapter_name);
	ret = 1;
	goto cleanup;
    }
    bu_vls_trunc(&adapter_name, 0);
    if (!rt_view_context_unique_object_name_bsg(&adapter_name,
	    "view_poly_0", v) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&direct_name), bu_vls_cstr(&adapter_name))) {
	printf("FAIL: BSG context unique-object-name adapter\n");
	bu_vls_free(&direct_name);
	bu_vls_free(&adapter_name);
	ret = 1;
	goto cleanup;
    }
    bu_vls_trunc(&adapter_name, 0);
    if (!rt_view_context_unique_object_name(&adapter_name,
	    "view_poly_0", v) ||
	    !BU_STR_EQUAL(bu_vls_cstr(&direct_name), bu_vls_cstr(&adapter_name))) {
	printf("FAIL: retained context unique-object-name adapter\n");
	bu_vls_free(&direct_name);
	bu_vls_free(&adapter_name);
	ret = 1;
	goto cleanup;
    }
    bu_vls_free(&direct_name);
    bu_vls_free(&adapter_name);

    struct bsg_interactive_rect_state direct_rect;
    struct rt_view_interactive_rect_state neutral_zero_rect = RT_VIEW_INTERACTIVE_RECT_STATE_INIT;
    struct rt_view_interactive_rect_state neutral_adapter_rect;
    struct rt_view_interactive_rect_state neutral_updated_rect = RT_VIEW_INTERACTIVE_RECT_STATE_INIT;
    memset(&neutral_zero_rect, 0, sizeof(neutral_zero_rect));
    memset(&neutral_updated_rect, 0, sizeof(neutral_updated_rect));
    neutral_updated_rect.active = 1;
    neutral_updated_rect.draw = 0;
    neutral_updated_rect.line_width = 7;
    neutral_updated_rect.line_style = 0;
    neutral_updated_rect.pos[0] = 21;
    neutral_updated_rect.pos[1] = 22;
    neutral_updated_rect.dim[0] = 201;
    neutral_updated_rect.dim[1] = 202;
    neutral_updated_rect.x = -0.15;
    neutral_updated_rect.y = 0.45;
    neutral_updated_rect.width = 0.25;
    neutral_updated_rect.height = 0.3;
    neutral_updated_rect.bg[0] = 11;
    neutral_updated_rect.bg[1] = 12;
    neutral_updated_rect.bg[2] = 13;
    neutral_updated_rect.color[0] = 14;
    neutral_updated_rect.color[1] = 15;
    neutral_updated_rect.color[2] = 16;
    neutral_updated_rect.cdim[0] = 500;
    neutral_updated_rect.cdim[1] = 250;
    neutral_updated_rect.aspect = 2.0;
    if (rt_view_interactive_rect_state_from_bsg(&neutral_adapter_rect, NULL) ||
	    rt_view_context_interactive_rect_state_from_bsg(&neutral_adapter_rect, NULL) ||
	    rt_view_context_interactive_rect_state_get(&neutral_adapter_rect, NULL) ||
	    memcmp(&neutral_adapter_rect, &neutral_zero_rect, sizeof(neutral_adapter_rect))) {
	printf("FAIL: null neutral interactive rectangle get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_interactive_rect_state_set_bsg(NULL, &neutral_updated_rect) ||
	    rt_view_context_interactive_rect_state_set_bsg(NULL, &neutral_updated_rect) ||
	    rt_view_context_interactive_rect_state_set(NULL, &neutral_updated_rect) ||
	    rt_view_interactive_rect_state_set_bsg(v, NULL) ||
	    rt_view_context_interactive_rect_state_set_bsg(v, NULL) ||
	    rt_view_context_interactive_rect_state_set(v, NULL) ||
	    rt_view_context_interactive_rect_state_from_bsg(NULL, v) ||
	    rt_view_context_interactive_rect_state_get(NULL, v) ||
	    rt_view_interactive_rect_state_from_bsg(NULL, v)) {
	printf("FAIL: null neutral interactive rectangle adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_interactive_rect_state_set_bsg(v, &neutral_updated_rect) ||
	    !rt_view_interactive_rect_state_from_bsg(&neutral_adapter_rect, v) ||
	    memcmp(&neutral_adapter_rect, &neutral_updated_rect, sizeof(neutral_adapter_rect)) ||
	    !bsg_view_interactive_rect_get(v, &direct_rect) ||
	    direct_rect.active != neutral_updated_rect.active ||
	    direct_rect.draw != neutral_updated_rect.draw ||
	    direct_rect.line_width != neutral_updated_rect.line_width ||
	    direct_rect.line_style != neutral_updated_rect.line_style ||
	    direct_rect.pos[0] != neutral_updated_rect.pos[0] ||
	    direct_rect.pos[1] != neutral_updated_rect.pos[1] ||
	    direct_rect.dim[0] != neutral_updated_rect.dim[0] ||
	    direct_rect.dim[1] != neutral_updated_rect.dim[1] ||
	    !fastf_equal(direct_rect.x, neutral_updated_rect.x) ||
	    !fastf_equal(direct_rect.y, neutral_updated_rect.y) ||
	    !fastf_equal(direct_rect.width, neutral_updated_rect.width) ||
	    !fastf_equal(direct_rect.height, neutral_updated_rect.height) ||
	    direct_rect.bg[0] != neutral_updated_rect.bg[0] ||
	    direct_rect.bg[1] != neutral_updated_rect.bg[1] ||
	    direct_rect.bg[2] != neutral_updated_rect.bg[2] ||
	    direct_rect.color[0] != neutral_updated_rect.color[0] ||
	    direct_rect.color[1] != neutral_updated_rect.color[1] ||
	    direct_rect.color[2] != neutral_updated_rect.color[2] ||
	    direct_rect.cdim[0] != neutral_updated_rect.cdim[0] ||
	    direct_rect.cdim[1] != neutral_updated_rect.cdim[1] ||
	    !fastf_equal(direct_rect.aspect, neutral_updated_rect.aspect)) {
	printf("FAIL: neutral interactive rectangle adapter\n");
	ret = 1;
	goto cleanup;
    }
    neutral_updated_rect.line_width = 11;
    neutral_updated_rect.pos[0] = 31;
    neutral_updated_rect.pos[1] = 32;
    if (!rt_view_context_interactive_rect_state_set_bsg(v, &neutral_updated_rect) ||
	    !rt_view_context_interactive_rect_state_from_bsg(&neutral_adapter_rect, v) ||
	    memcmp(&neutral_adapter_rect, &neutral_updated_rect, sizeof(neutral_adapter_rect))) {
	printf("FAIL: neutral interactive rectangle context adapter\n");
	ret = 1;
	goto cleanup;
    }
    neutral_updated_rect.line_width = 13;
    neutral_updated_rect.pos[0] = 41;
    neutral_updated_rect.pos[1] = 42;
    if (!rt_view_context_interactive_rect_state_set(v, &neutral_updated_rect) ||
	    !rt_view_context_interactive_rect_state_get(&neutral_adapter_rect, v) ||
	    memcmp(&neutral_adapter_rect, &neutral_updated_rect, sizeof(neutral_adapter_rect))) {
	printf("FAIL: retained context interactive rectangle adapter\n");
	ret = 1;
	goto cleanup;
    }

    if ((RT_VIEW_REFRESH_ALL_BSG & (RT_VIEW_REFRESH_VIEW_BSG |
		    RT_VIEW_REFRESH_DRAW_BSG |
		    RT_VIEW_REFRESH_EDIT_BSG |
		    RT_VIEW_REFRESH_FRAMEBUFFER_BSG |
		    RT_VIEW_REFRESH_OVERLAY_BSG |
		    RT_VIEW_REFRESH_FORCE_BSG)) !=
	    (RT_VIEW_REFRESH_VIEW_BSG |
	     RT_VIEW_REFRESH_DRAW_BSG |
	     RT_VIEW_REFRESH_EDIT_BSG |
	     RT_VIEW_REFRESH_FRAMEBUFFER_BSG |
	     RT_VIEW_REFRESH_OVERLAY_BSG |
	     RT_VIEW_REFRESH_FORCE_BSG)) {
	printf("FAIL: BSG refresh flag alias coverage\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_refresh_request_bsg(NULL, RT_VIEW_REFRESH_VIEW_BSG) ||
	    rt_view_context_refresh_request_bsg(NULL, RT_VIEW_REFRESH_VIEW_BSG) ||
	    rt_view_context_refresh_request(NULL, RT_VIEW_REFRESH_VIEW_BSG) ||
	    rt_view_refresh_dirty_from_bsg(NULL) ||
	    rt_view_context_refresh_dirty_from_bsg(NULL) ||
	    rt_view_context_refresh_dirty_get(NULL) ||
	    rt_view_refresh_consume_bsg(NULL) ||
	    rt_view_context_refresh_consume_bsg(NULL) ||
	    rt_view_context_refresh_consume(NULL) ||
	    rt_view_refresh_complete_bsg(NULL) ||
	    rt_view_context_refresh_complete_bsg(NULL) ||
	    rt_view_context_refresh_complete(NULL) ||
	    rt_view_refresh_enabled_from_bsg(NULL) != 1 ||
	    rt_view_context_refresh_enabled_from_bsg(NULL) != 1 ||
	    rt_view_context_refresh_enabled_get(NULL) != 1 ||
	    rt_view_refresh_enabled_set_bsg(NULL, 1) ||
	    rt_view_context_refresh_enabled_set_bsg(NULL, 1) ||
	    rt_view_context_refresh_enabled_set(NULL, 1) ||
	    rt_view_refresh_suppressed_from_bsg(NULL) ||
	    rt_view_context_refresh_suppressed_from_bsg(NULL) ||
	    rt_view_context_refresh_suppressed_get(NULL) ||
	    rt_view_refresh_suppress_begin_bsg(NULL) ||
	    rt_view_context_refresh_suppress_begin_bsg(NULL) ||
	    rt_view_context_refresh_suppress_begin(NULL) ||
	    rt_view_refresh_suppress_end_bsg(NULL) ||
	    rt_view_context_refresh_suppress_end_bsg(NULL) ||
	    rt_view_context_refresh_suppress_end(NULL) ||
	    rt_view_refresh_drawn_count_from_bsg(NULL) ||
	    rt_view_context_refresh_drawn_count_from_bsg(NULL) ||
	    rt_view_context_refresh_drawn_count_get(NULL) ||
	    rt_view_context_refresh_drawn_count_set_bsg(NULL, 1) ||
	    rt_view_context_refresh_drawn_count_set(NULL, 1) ||
	    rt_view_refresh_drawn_count_set_bsg(NULL, 1)) {
	printf("FAIL: null BSG refresh adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_refresh_complete_bsg(v);
    if (!rt_view_refresh_enabled_set_bsg(v, 1) ||
	    !rt_view_refresh_request_bsg(v, RT_VIEW_REFRESH_VIEW_BSG) ||
	    !rt_view_refresh_dirty_from_bsg(v) ||
	    !(rt_view_refresh_consume_bsg(v) & RT_VIEW_REFRESH_VIEW_BSG) ||
	    rt_view_refresh_dirty_from_bsg(v) ||
	    !rt_view_refresh_complete_bsg(v)) {
	printf("FAIL: BSG refresh dirty/consume adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_refresh_complete_bsg(v);
    if (!rt_view_context_refresh_request_bsg(v, RT_VIEW_REFRESH_VIEW_BSG) ||
	    !rt_view_context_refresh_dirty_from_bsg(v) ||
	    !(rt_view_context_refresh_consume_bsg(v) & RT_VIEW_REFRESH_VIEW_BSG) ||
	    rt_view_context_refresh_dirty_from_bsg(v) ||
	    !rt_view_context_refresh_complete_bsg(v)) {
	printf("FAIL: BSG context refresh dirty/consume adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_refresh_request(v, RT_VIEW_REFRESH_VIEW_BSG) ||
	    !rt_view_context_refresh_dirty_get(v) ||
	    !(rt_view_context_refresh_consume(v) & RT_VIEW_REFRESH_VIEW_BSG) ||
	    rt_view_context_refresh_dirty_get(v) ||
	    !rt_view_context_refresh_complete(v)) {
	printf("FAIL: retained context refresh dirty/consume adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_refresh_suppress_begin(v) ||
	    !rt_view_context_refresh_suppressed_get(v) ||
	    !rt_view_context_refresh_request_bsg(v, RT_VIEW_REFRESH_DRAW_BSG) ||
	    rt_view_context_refresh_dirty_get(v) ||
	    rt_view_context_refresh_consume(v) ||
	    !rt_view_context_refresh_suppress_end(v) ||
	    !rt_view_context_refresh_dirty_get(v) ||
	    !rt_view_context_refresh_complete(v) ||
	    rt_view_context_refresh_suppressed_get(v)) {
	printf("FAIL: BSG refresh suppress adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_refresh_enabled_set(v, 0) ||
	    rt_view_context_refresh_enabled_get(v) ||
	    !rt_view_context_refresh_request_bsg(v, RT_VIEW_REFRESH_EDIT_BSG) ||
	    rt_view_context_refresh_dirty_get(v) ||
	    rt_view_context_refresh_consume(v) ||
	    !rt_view_context_refresh_enabled_set(v, 1) ||
	    !rt_view_context_refresh_dirty_get(v) ||
	    !rt_view_context_refresh_complete(v)) {
	printf("FAIL: BSG refresh enabled adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_context_refresh_drawn_count_set(v, 3) ||
	    rt_view_context_refresh_drawn_count_get(v) != 3 ||
	    !rt_view_context_refresh_drawn_count_set(v, -1) ||
	    rt_view_context_refresh_drawn_count_get(v) != 0) {
	printf("FAIL: BSG refresh drawn-count adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_lod_bounds_callback_is_bsg(NULL) ||
	    rt_view_context_lod_bounds_callback_is(NULL) ||
	    rt_view_lod_bounds_callback_is_bsg(v) ||
	    rt_view_context_lod_bounds_callback_is(v)) {
	printf("FAIL: BSG bounds callback initial adapter state\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_lod_bounds_callback_set_bsg(NULL);
    rt_view_lod_bounds_callback_set_bsg(v);
    if (!rt_view_lod_bounds_callback_is_bsg(v) ||
	    !rt_view_context_lod_bounds_callback_is(v)) {
	printf("FAIL: BSG bounds callback set adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    if (shared && shared->settings_local) {
	shared->settings_active = shared->settings_local;
	bsg_view_state_sync_from_view(shared->gv_state, shared);
    }
    free_view(shared);
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_edit_view_adapter(void)
{
    struct bsg_view *v = make_view("rt_edit_view_adapter");
    struct rt_edit_view ev;
    mat_t expected_mat;
    int ret = 0;

    rt_edit_view_from_bsg(&ev, NULL);
    MAT_IDN(expected_mat);
    if (!fastf_equal(ev.gv_scale, 1.0) ||
	    !fastf_equal(ev.gv_base2local, 1.0) ||
	    !fastf_equal(ev.gv_local2base, 1.0) ||
	    ev.gv_coord != 'v' ||
	    ev.gv_rotate_about != 'v' ||
	    memcmp(ev.gv_rotation, expected_mat, sizeof(mat_t)) ||
	    memcmp(ev.gv_center, expected_mat, sizeof(mat_t)) ||
	    memcmp(ev.gv_model2view, expected_mat, sizeof(mat_t)) ||
	    memcmp(ev.gv_view2model, expected_mat, sizeof(mat_t))) {
	printf("FAIL: null BSG edit-view adapter did not return defaults\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_scale = 3.0;
    v->gv_base2local = 4.0;
    v->gv_local2base = 0.25;
    v->gv_coord = 'm';
    v->gv_rotate_about = 'k';
    bn_mat_angles(v->gv_rotation, 11.0, 22.0, 33.0);
    bn_mat_angles(v->gv_center, 44.0, 55.0, 66.0);
    MAT_DELTAS(v->gv_center, 7.0, 8.0, 9.0);
    bn_mat_angles(v->gv_model2view, 10.0, 20.0, 30.0);
    MAT_DELTAS(v->gv_model2view, 1.0, 2.0, 3.0);
    bn_mat_angles(v->gv_view2model, 30.0, 20.0, 10.0);
    MAT_DELTAS(v->gv_view2model, 3.0, 2.0, 1.0);

    rt_edit_view_from_bsg(&ev, v);
    if (!fastf_equal(ev.gv_scale, v->gv_scale) ||
	    !fastf_equal(ev.gv_base2local, v->gv_base2local) ||
	    !fastf_equal(ev.gv_local2base, v->gv_local2base) ||
	    ev.gv_coord != v->gv_coord ||
	    ev.gv_rotate_about != v->gv_rotate_about ||
	    memcmp(ev.gv_rotation, v->gv_rotation, sizeof(mat_t)) ||
	    memcmp(ev.gv_center, v->gv_center, sizeof(mat_t)) ||
	    memcmp(ev.gv_model2view, v->gv_model2view, sizeof(mat_t)) ||
	    memcmp(ev.gv_view2model, v->gv_view2model, sizeof(mat_t))) {
	printf("FAIL: BSG edit-view adapter did not copy expected fields\n");
	ret = 1;
	goto cleanup;
    }

    rt_edit_view_from_context(&ev, v);
    if (!fastf_equal(ev.gv_scale, v->gv_scale) ||
	    !fastf_equal(ev.gv_base2local, v->gv_base2local) ||
	    !fastf_equal(ev.gv_local2base, v->gv_local2base) ||
	    ev.gv_coord != v->gv_coord ||
	    ev.gv_rotate_about != v->gv_rotate_about ||
	    memcmp(ev.gv_rotation, v->gv_rotation, sizeof(mat_t)) ||
	    memcmp(ev.gv_center, v->gv_center, sizeof(mat_t)) ||
	    memcmp(ev.gv_model2view, v->gv_model2view, sizeof(mat_t)) ||
	    memcmp(ev.gv_view2model, v->gv_view2model, sizeof(mat_t))) {
	printf("FAIL: context edit-view adapter did not copy expected fields\n");
	ret = 1;
	goto cleanup;
    }

    struct rt_edit *s = rt_edit_create_context(NULL, NULL, NULL, v);
    if (!s || !s->vp ||
	    !fastf_equal(s->vp->gv_scale, v->gv_scale) ||
	    s->vp->gv_coord != v->gv_coord) {
	printf("FAIL: context edit create adapter did not set expected view\n");
	rt_edit_destroy(s);
	ret = 1;
	goto cleanup;
    }

    if (rt_edit_reinit_context(s, NULL, NULL, NULL, v) != BRLCAD_OK ||
	    !s->vp ||
	    !fastf_equal(s->vp->gv_base2local, v->gv_base2local) ||
	    s->vp->gv_rotate_about != v->gv_rotate_about) {
	printf("FAIL: context edit reinit adapter did not set expected view\n");
	rt_edit_destroy(s);
	ret = 1;
	goto cleanup;
    }

    if (rt_edit_knob_cmd_process_context(NULL, NULL, NULL, NULL, NULL, NULL,
		v, "x", 0.0, 'v', 0, NULL) != BRLCAD_ERROR) {
	printf("FAIL: context edit knob adapter did not reject null arguments\n");
	rt_edit_destroy(s);
	ret = 1;
	goto cleanup;
    }
    rt_edit_destroy(s);

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    if (test_bsg_null_adapter())
	return 1;
    if (test_bsg_obb_adapter())
	return 1;
    if (test_bsg_user_data_adapter())
	return 1;
    if (test_bsg_adapter())
	return 1;
    if (test_bsg_adapter_sanitizes())
	return 1;
    if (test_bsg_orientation_adapter())
	return 1;
    if (test_bsg_aet_adapter())
	return 1;
    if (test_bsg_update_adapter())
	return 1;
    if (test_bsg_autoview_adapter())
	return 1;
    if (test_bsg_adjust_adapter())
	return 1;
    if (test_bsg_hash_adapter())
	return 1;
    if (test_bsg_knob_adapter())
	return 1;
    if (test_bsg_snap_adapter())
	return 1;
    if (test_bsg_measure_adapter())
	return 1;
    if (test_bsg_feature_adapter())
	return 1;
    if (test_bsg_view_scope_adapter())
	return 1;
    if (test_bsg_perspective_adapter())
	return 1;
    if (test_bsg_display_manager_adapter())
	return 1;
    if (test_bsg_view_runtime_adapter())
	return 1;
    if (test_bsg_camera_adapter())
	return 1;
    if (test_bsg_lod_policy_adapter())
	return 1;
    if (test_bsg_faceplate_state_adapter())
	return 1;
    if (test_bsg_scene_ref_adapter_boundary())
	return 1;
    if (test_bsg_edit_view_adapter())
	return 1;

    return 0;
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
