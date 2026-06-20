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

#include "bn/mat.h"
#include "bn/qmath.h"
#include "bu/app.h"
#include "bu/malloc.h"
#include "bsg/util.h"
#include "bsg/view_state.h"
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

    if (!rt_view_center_from_bsg(actual_mat, v)) {
	printf("FAIL: BSG center adapter should report a source view\n");
	ret = 1;
	goto cleanup;
    }
    if (memcmp(actual_mat, v->gv_center, sizeof(mat_t))) {
	printf("FAIL: BSG center adapter did not copy matrix\n");
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

cleanup:
    free_view(src);
    free_view(dst);
    return ret ? 1 : 0;
}

static int
test_bsg_mesh_lod_adapter_boundary(void)
{
    struct bsg_view *v = make_view("rt_view_lod_bounds_callback");
    bsg_scene_ref null_ref = bsg_scene_ref_null();
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
