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
#include "bsg/util.h"
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

static void
test_bounds_update_callback(struct bsg_view *UNUSED(v))
{
    test_bounds_update_count++;
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

    info.width = 640;
    info.height = 480;
    info.size = 1000.0;
    info.lod.scale = 4.0;
    info.lod.curve_scale = 5.0;
    info.lod.point_scale = 6.0;
    info.lod.bot_threshold = 7;

    rt_view_info_from_bsg(NULL, NULL);
    rt_view_info_from_bsg(&info, NULL);

    return check_view_info("null bsg adapter", &info, 1, 1, 1.0,
	    1.0, 1.0, 1.0, 0);
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

    if (rt_view_orientation_quat_from_bsg(actual, NULL)) {
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

cleanup:
    free_view(v);
    return ret ? 1 : 0;
}

static int
test_bsg_perspective_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_perspective_adapter");
    int ret = 0;

    if (!fastf_equal(rt_view_perspective_from_bsg(NULL), 0.0)) {
	printf("FAIL: null BSG perspective adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_perspective = 37.5;
    if (!fastf_equal(rt_view_perspective_from_bsg(v), 37.5)) {
	printf("FAIL: BSG perspective adapter got %g\n",
		rt_view_perspective_from_bsg(v));
	ret = 1;
	goto cleanup;
    }

    if (rt_view_perspective_set_bsg(NULL, 12.5)) {
	printf("FAIL: null BSG perspective set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_perspective_set_bsg(v, 12.5) ||
	    !fastf_equal(rt_view_perspective_from_bsg(v), 12.5)) {
	printf("FAIL: BSG perspective set adapter\n");
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
    if (rt_view_aet_from_bsg(actual, NULL)) {
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
    if (!rt_view_aet_from_bsg(actual, v)) {
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
	    rt_view_aet_set_bsg(v, NULL)) {
	printf("FAIL: null BSG AET set adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(expected, 21.0, 22.0, 23.0);
    if (!rt_view_aet_set_bsg(v, expected) ||
	    !rt_view_aet_from_bsg(actual, v) ||
	    !VNEAR_EQUAL(actual, expected, SMALL_FASTF)) {
	printf("FAIL: BSG AET set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_aet_state_set_bsg(NULL, expected) ||
	    rt_view_aet_state_set_bsg(v, NULL)) {
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

cleanup:
    free_view(v);
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
    if (rt_view_model2view_set_bsg(NULL, v->gv_model2view) ||
	    rt_view_model2view_set_bsg(v, NULL)) {
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
    if (rt_view_view2model_set_bsg(NULL, v->gv_view2model) ||
	    rt_view_view2model_set_bsg(v, NULL)) {
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
    if (rt_view_pmodel2view_set_bsg(NULL, v->gv_pmodel2view) ||
	    rt_view_pmodel2view_set_bsg(v, NULL)) {
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
    if (rt_view_pmat_set_bsg(NULL, v->gv_pmat) ||
	    rt_view_pmat_set_bsg(v, NULL)) {
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

    if (rt_view_rotation_set_bsg(NULL, v->gv_rotation) ||
	    rt_view_rotation_set_bsg(v, NULL)) {
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

    if (!rt_view_center_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG center adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_center_vec_set_bsg(NULL, v->gv_eye_pos) ||
	    rt_view_center_vec_set_bsg(v, NULL)) {
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

    plane_t actual_plane;
    plane_t expected_plane;
    point_t plane_center;
    vect_t plane_normal;
    if (rt_view_plane_from_bsg(NULL, v) != -1 ||
	    rt_view_plane_from_bsg(&actual_plane, NULL) != -1) {
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
    if (rt_view_eye_pos_set_bsg(NULL, v->gv_eye_pos) ||
	    rt_view_eye_pos_set_bsg(v, NULL)) {
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

    if (!rt_view_lod_policy_from_bsg(&policy, src)) {
	printf("FAIL: rt lod policy from bsg\n");
	ret = 1;
	goto cleanup;
    }
    ret += check_lod_policy("from bsg", &policy, 1, 1, 1, 2.5, 3.5, 4.5, 2468);

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

    rt_view_lod_bounds_update_bsg(NULL);
    rt_view_lod_bounds_update_bsg(dst);

    if (rt_view_bounds_update_callback_from_bsg(NULL)) {
	printf("FAIL: null BSG bounds callback get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_set_bsg(NULL, test_bounds_update_callback)) {
	printf("FAIL: null BSG bounds callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_call_bsg(NULL)) {
	printf("FAIL: null BSG bounds callback call adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_from_bsg(dst)) {
	printf("FAIL: initial BSG bounds callback should be null\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_bounds_update_callback_set_bsg(dst, test_bounds_update_callback)) {
	printf("FAIL: BSG bounds callback set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_bounds_update_callback_from_bsg(dst) != test_bounds_update_callback) {
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
    if (!rt_view_bounds_update_callback_set_bsg(dst, NULL) ||
	    rt_view_bounds_update_callback_call_bsg(dst)) {
	printf("FAIL: BSG bounds callback clear adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(src);
    free_view(dst);
    return ret ? 1 : 0;
}

static int
test_bsg_mesh_lod_adapter_boundary(void)
{
    struct bsg_view *v = make_view("rt_view_lod_bounds_callback");
    struct bsg_view *shared = make_view("rt_view_shared_settings_adapter");
    bsg_scene_ref null_ref = bsg_scene_ref_null();
    fastf_t *scale_storage = NULL;
    int ret = 0;

    if (rt_mesh_lod_load_view_scene_ref_bsg(NULL, null_ref, v, 0) != -1) {
	printf("FAIL: null mesh lod view-scene-ref load adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_mesh_lod_free_scene_ref_bsg(null_ref);

    if (!fastf_equal(rt_view_scale_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG view scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_scale = 42.0;
    if (!fastf_equal(rt_view_scale_from_bsg(v), 42.0)) {
	printf("FAIL: BSG view scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scale_set_bsg(NULL, 25.0)) {
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
    if (rt_view_scale_storage_from_bsg(NULL)) {
	printf("FAIL: null BSG view scale storage adapter\n");
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
    if (rt_view_scale_state_set_bsg(NULL, 1.0, 2.0, 3.0, 4.0, 5.0)) {
	printf("FAIL: null BSG view scale state set adapter\n");
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

    if (!fastf_equal(rt_view_initial_scale_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG initial view scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_i_scale = 500.0;
    if (!fastf_equal(rt_view_initial_scale_from_bsg(v), 500.0)) {
	printf("FAIL: BSG initial view scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_initial_scale_set_bsg(NULL, 250.0)) {
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

    if (!fastf_equal(rt_view_absolute_scale_from_bsg(NULL), 0.0)) {
	printf("FAIL: null BSG absolute view scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_a_scale = -0.25;
    if (!fastf_equal(rt_view_absolute_scale_from_bsg(v), -0.25)) {
	printf("FAIL: BSG absolute view scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_absolute_scale_set_bsg(NULL, 0.25)) {
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

    if (!fastf_equal(rt_view_local2base_from_bsg(NULL), 1.0) ||
	    !fastf_equal(rt_view_base2local_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG unit-conversion scale adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_local2base = 0.5;
    v->gv_base2local = 2.0;
    if (!fastf_equal(rt_view_local2base_from_bsg(v), 0.5) ||
	    !fastf_equal(rt_view_base2local_from_bsg(v), 2.0)) {
	printf("FAIL: BSG unit-conversion scale adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_unit_conversion_set_bsg(NULL, 3.0, 4.0)) {
	printf("FAIL: null BSG unit-conversion set adapter\n");
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

    if (!fastf_equal(rt_view_size_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_size = 84.0;
    if (!fastf_equal(rt_view_size_from_bsg(v), 84.0)) {
	printf("FAIL: BSG view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_size_set_bsg(NULL, 128.0)) {
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

    if (!fastf_equal(rt_view_inverse_size_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG inverse view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_isize = 1.0 / 84.0;
    if (!fastf_equal(rt_view_inverse_size_from_bsg(v), 1.0 / 84.0)) {
	printf("FAIL: BSG inverse view size adapter\n");
	ret = 1;
	goto cleanup;
    }

    point_t keypoint = VINIT_ZERO;
    if (rt_view_keypoint_from_bsg(keypoint, NULL) ||
	    !fastf_equal(keypoint[X], 0.0) ||
	    !fastf_equal(keypoint[Y], 0.0) ||
	    !fastf_equal(keypoint[Z], 0.0)) {
	printf("FAIL: null BSG keypoint adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSET(v->gv_keypoint, 1.0, 2.0, 3.0);
    if (!rt_view_keypoint_from_bsg(keypoint, v) ||
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
	    rt_view_keypoint_set_bsg(v, NULL)) {
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

    if (rt_view_rotate_about_from_bsg(NULL) != 'v') {
	printf("FAIL: null BSG rotate-about adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_rotate_about = 'k';
    if (rt_view_rotate_about_from_bsg(v) != 'k') {
	printf("FAIL: BSG rotate-about adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_rotate_about_set_bsg(NULL, 'm')) {
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

    if (rt_view_coord_from_bsg(NULL) != 'v') {
	printf("FAIL: null BSG coord adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->gv_coord = 'm';
    if (rt_view_coord_from_bsg(v) != 'm') {
	printf("FAIL: BSG coord adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_coord_set_bsg(NULL, 'v')) {
	printf("FAIL: null BSG coord set adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_coord_set_bsg(v, 'v') ||
	    rt_view_coord_from_bsg(v) != 'v') {
	printf("FAIL: BSG coord set adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_snap_lines_from_bsg(NULL) ||
	    rt_view_snap_source_flags_from_bsg(NULL) ||
	    rt_view_snap_kind_mask_from_bsg(NULL) ||
	    rt_view_snap_lines_set_bsg(NULL, 1) ||
	    rt_view_snap_source_flags_set_bsg(NULL, BSG_SNAP_TCL)) {
	printf("FAIL: null BSG snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_lines_set_bsg(v, 1) ||
	    !rt_view_snap_source_flags_set_bsg(v, BSG_SNAP_TCL) ||
	    !rt_view_snap_lines_from_bsg(v) ||
	    rt_view_snap_source_flags_from_bsg(v) != BSG_SNAP_TCL ||
	    !(rt_view_snap_kind_mask_from_bsg(v) & BSG_SNAP_KIND_ENDPOINT)) {
	printf("FAIL: BSG snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_source_flags_set_bsg(v, BSG_SNAP_DB) ||
	    !(rt_view_prepare_tcl_snap_bsg(v) & BSG_SNAP_KIND_ENDPOINT) ||
	    rt_view_snap_source_flags_from_bsg(v) != BSG_SNAP_TCL ||
	    rt_view_prepare_tcl_snap_bsg(NULL)) {
	printf("FAIL: BSG Tcl snap preparation adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_center_linesnap_bsg(NULL) ||
	    !rt_view_center_linesnap_bsg(v)) {
	printf("FAIL: BSG line-snap recenter adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_zclip_from_bsg(NULL) ||
	    rt_view_zclip_set_bsg(NULL, 1)) {
	printf("FAIL: null BSG zclip adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_zclip_set_bsg(v, 1) ||
	    rt_view_zclip_from_bsg(v) != 1 ||
	    !rt_view_zclip_set_bsg(v, 0) ||
	    rt_view_zclip_from_bsg(v)) {
	printf("FAIL: BSG zclip adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_framebuffer_mode_from_bsg(NULL) ||
	    rt_view_framebuffer_mode_set_bsg(NULL, 2)) {
	printf("FAIL: null BSG framebuffer mode adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_framebuffer_mode_set_bsg(v, 2) ||
	    rt_view_framebuffer_mode_from_bsg(v) != 2 ||
	    !rt_view_framebuffer_mode_set_bsg(v, 1) ||
	    rt_view_framebuffer_mode_from_bsg(v) != 1 ||
	    !rt_view_framebuffer_mode_set_bsg(v, 0) ||
	    rt_view_framebuffer_mode_from_bsg(v)) {
	printf("FAIL: BSG framebuffer mode adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_cleared_from_bsg(NULL) ||
	    rt_view_cleared_set_bsg(NULL, 1)) {
	printf("FAIL: null BSG cleared-state adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_cleared_set_bsg(v, 1) ||
	    rt_view_cleared_from_bsg(v) != 1 ||
	    !rt_view_cleared_set_bsg(v, 2) ||
	    rt_view_cleared_from_bsg(v) != 1 ||
	    !rt_view_cleared_set_bsg(v, 0) ||
	    rt_view_cleared_from_bsg(v)) {
	printf("FAIL: BSG cleared-state adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_settings_shared_bsg(NULL, v) ||
	    rt_view_settings_shared_bsg(v, NULL) ||
	    rt_view_settings_shared_bsg(v, shared)) {
	printf("FAIL: distinct BSG shared-settings adapter\n");
	ret = 1;
	goto cleanup;
    }
    shared->settings_active = v->settings_active;
    bsg_view_state_sync_from_view(shared->gv_state, shared);
    if (!rt_view_settings_shared_bsg(v, shared)) {
	printf("FAIL: shared BSG shared-settings adapter\n");
	ret = 1;
	goto cleanup;
    }
    shared->settings_active = shared->settings_local;
    bsg_view_state_sync_from_view(shared->gv_state, shared);

    if (!fastf_equal(rt_view_snap_tolerance_factor_from_bsg(NULL), 1.0) ||
	    rt_view_snap_tolerance_factor_set_bsg(NULL, 2.0)) {
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

    if (rt_view_width_from_bsg(NULL) != 0 ||
	    rt_view_height_from_bsg(NULL) != 0) {
	printf("FAIL: null BSG raw width/height adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_dimensions_set_bsg(NULL, 123, 456)) {
	printf("FAIL: null BSG raw dimensions set adapter\n");
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

    if (rt_view_lod_bounds_callback_is_bsg(NULL) ||
	    rt_view_lod_bounds_callback_is_bsg(v)) {
	printf("FAIL: BSG bounds callback initial adapter state\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_lod_bounds_callback_set_bsg(NULL);
    rt_view_lod_bounds_callback_set_bsg(v);
    if (!rt_view_lod_bounds_callback_is_bsg(v)) {
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
    if (test_bsg_adapter())
	return 1;
    if (test_bsg_adapter_sanitizes())
	return 1;
    if (test_bsg_orientation_adapter())
	return 1;
    if (test_bsg_aet_adapter())
	return 1;
    if (test_bsg_perspective_adapter())
	return 1;
    if (test_bsg_camera_adapter())
	return 1;
    if (test_bsg_lod_policy_adapter())
	return 1;
    if (test_bsg_mesh_lod_adapter_boundary())
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
