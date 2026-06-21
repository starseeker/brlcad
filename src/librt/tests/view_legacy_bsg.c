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
#include "bsg/scene_builder.h"
#include "bsg/selection.h"
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
test_bsg_update_adapter(void)
{
    struct bsg_view *direct = make_view("rt_view_update_direct");
    struct bsg_view *adapter = make_view("rt_view_update_adapter");
    int ret = 0;

    if (rt_view_update_bsg(NULL)) {
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
    if (!rt_view_update_bsg(adapter)) {
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
    struct bsg_view *direct_empty = make_view("rt_view_autoview_empty_direct");
    struct bsg_view *adapter_empty = make_view("rt_view_autoview_empty_adapter");
    point_t min, max;
    int ret = 0;

    VSET(min, -1.0, -2.0, -3.0);
    VSET(max, 5.0, 6.0, 7.0);

    if (rt_view_autoview_bsg(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0) ||
	    rt_view_autoview_bounds_bsg(NULL, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, max) ||
	    rt_view_autoview_bounds_bsg(adapter_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, NULL, max) ||
	    rt_view_autoview_bounds_bsg(adapter_bounds, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, min, NULL)) {
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

    bsg_autoview(direct_empty, BSG_AUTOVIEW_SCALE_DEFAULT, 0);
    if (!rt_view_autoview_bsg(adapter_empty, RT_VIEW_AUTOVIEW_SCALE_DEFAULT, 0)) {
	printf("FAIL: BSG autoview adapter should report success\n");
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
    point_t direct_keypoint;
    point_t adapter_keypoint;
    int ret = 0;

    setup_adjust_view(direct);
    setup_adjust_view(adapter);
    VMOVE(direct_keypoint, keypoint);
    VMOVE(adapter_keypoint, keypoint);

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

cleanup:
    free_view(direct);
    free_view(adapter);
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
	    rt_view_adjust_bsg(v, 1, 2, NULL, 0, RT_VIEW_ADJUST_SCALE) ||
	    rt_view_adjust_bsg(v, 1, 2, keypoint, 0, RT_VIEW_ADJUST_IDLE)) {
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

    if (rt_view_hash_bsg(NULL) != 0ULL) {
	printf("FAIL: null BSG hash adapter\n");
	ret = 1;
	goto cleanup;
    }

    setup_adjust_view(v);
    direct_hash = bsg_hash(v);
    adapter_hash = rt_view_hash_bsg(v);
    if (direct_hash != adapter_hash) {
	printf("FAIL: BSG hash adapter: got %llu expected %llu\n",
		adapter_hash, direct_hash);
	ret = 1;
	goto cleanup;
    }

    v->gv_width += 7;
    direct_hash = bsg_hash(v);
    adapter_hash = rt_view_hash_bsg(v);
    if (direct_hash != adapter_hash) {
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
    int ret = 0;

    if (rt_view_knobs_reset_bsg(NULL, RT_VIEW_KNOBS_RATE_BSG) ||
	    rt_view_knob_state_reset_bsg(NULL, RT_VIEW_KNOBS_RATE_BSG) ||
	    rt_view_knobs_hash_bsg(NULL, NULL) ||
	    rt_view_knobs_translate_bsg(NULL, direct_tvec, 0) ||
	    rt_view_knobs_translate_bsg(adapter, NULL, 0) ||
	    rt_view_knobs_rotate_bsg(NULL, direct_rvec, 'v', 'v', NULL, NULL) ||
	    rt_view_knobs_rotate_bsg(adapter, NULL, 'v', 'v', NULL, NULL) ||
	    rt_view_knobs_update_rate_flags_bsg(NULL)) {
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
    adapter->k = direct->k;

    bsg_knobs_reset(&direct->k, BSG_KNOBS_RATE);
    if (!rt_view_knobs_reset_bsg(adapter, RT_VIEW_KNOBS_RATE_BSG) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob rate reset adapter\n");
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
    if (rt_view_knobs_hash_bsg(adapter, NULL) != bsg_knobs_hash(&direct->k, NULL)) {
	printf("FAIL: BSG knob hash adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (bsg_knobs_cmd_process(&direct_rvec, &direct_do_rot,
		&direct_tvec, &direct_do_tran, direct, "x", 6.0, 'v', 0, 0) !=
	    rt_view_knobs_cmd_process_bsg(&adapter_rvec, &adapter_do_rot,
		&adapter_tvec, &adapter_do_tran, adapter, "x", 6.0, 'v', 0, 0) ||
	    direct_do_rot != adapter_do_rot ||
	    direct_do_tran != adapter_do_tran ||
	    !VNEAR_EQUAL(direct_rvec, adapter_rvec, SMALL_FASTF) ||
	    !VNEAR_EQUAL(direct_tvec, adapter_tvec, SMALL_FASTF) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob command-process adapter\n");
	ret = 1;
	goto cleanup;
    }

    VSET(direct_tvec, 0.4, -0.2, 0.1);
    VMOVE(adapter_tvec, direct_tvec);
    bsg_knobs_tran(direct, direct_tvec, 0);
    if (!rt_view_knobs_translate_bsg(adapter, adapter_tvec, 0) ||
	    check_view_frame("BSG knob translate adapter", adapter, direct)) {
	ret = 1;
	goto cleanup;
    }

    VSET(direct_rvec, 2.0, 3.0, 4.0);
    VMOVE(adapter_rvec, direct_rvec);
    bsg_knobs_rot(direct, direct_rvec, 'v', 'v', NULL, NULL);
    if (!rt_view_knobs_rotate_bsg(adapter, adapter_rvec, 'v', 'v', NULL, NULL) ||
	    check_view_frame("BSG knob rotate adapter", adapter, direct)) {
	ret = 1;
	goto cleanup;
    }

    direct->k.sca = 1.0;
    direct->k.tra_m[Z] = 0.5;
    adapter->k = direct->k;
    bsg_update_rate_flags(direct);
    if (!rt_view_knobs_update_rate_flags_bsg(adapter) ||
	    memcmp(&adapter->k, &direct->k, sizeof(direct->k))) {
	printf("FAIL: BSG knob rate-flag adapter\n");
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
    point_t sample = VINIT_ZERO;
    fastf_t direct_x = 0.1;
    fastf_t direct_y = 0.1;
    fastf_t adapter_x = 0.1;
    fastf_t adapter_y = 0.1;
    fastf_t direct_grid_x = 0.1;
    fastf_t direct_grid_y = 0.1;
    fastf_t adapter_grid_x = 0.1;
    fastf_t adapter_grid_y = 0.1;
    bsg_feature_ref feature = {1234, 5678};
    bsg_feature_ref direct_feature = BSG_FEATURE_REF_NULL_INIT;
    bsg_feature_ref adapter_feature = BSG_FEATURE_REF_NULL_INIT;
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
	    rt_view_snap_grid_2d_bsg(adapter, &adapter_grid_x, NULL)) {
	printf("FAIL: null BSG snap adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    bsg_snap_result_free(&adapter_res);

    if (rt_view_snap_exclude_feature_from_bsg(NULL, adapter) ||
	    rt_view_snap_exclude_feature_from_bsg(&adapter_feature, NULL) ||
	    rt_view_snap_exclude_feature_set_bsg(NULL, feature) ||
	    rt_view_snap_exclude_feature_clear_bsg(NULL)) {
	printf("FAIL: null BSG snap-exclude adapter arguments\n");
	ret = 1;
	goto cleanup;
    }

    setup_snap_view(direct);
    setup_snap_view(adapter);

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

    if (bsg_snap_grid_2d(direct, &direct_grid_x, &direct_grid_y) !=
	    rt_view_snap_grid_2d_bsg(adapter, &adapter_grid_x, &adapter_grid_y) ||
	    !fastf_equal(direct_grid_x, adapter_grid_x) ||
	    !fastf_equal(direct_grid_y, adapter_grid_y)) {
	printf("FAIL: BSG 2D grid-snap adapter\n");
	ret = 1;
	goto cleanup;
    }

    bsg_view_snap_exclude_feature_set(direct, feature);
    if (!rt_view_snap_exclude_feature_set_bsg(adapter, feature) ||
	    !bsg_view_snap_exclude_feature_get(direct, &direct_feature) ||
	    !rt_view_snap_exclude_feature_from_bsg(&adapter_feature, adapter) ||
	    direct_feature.token != adapter_feature.token ||
	    direct_feature.revision != adapter_feature.revision) {
	printf("FAIL: BSG snap-exclude set/get adapter\n");
	ret = 1;
	goto cleanup;
    }

    bsg_view_snap_exclude_feature_clear(direct);
    direct_feature = feature;
    adapter_feature = feature;
    if (!rt_view_snap_exclude_feature_clear_bsg(adapter)) {
	printf("FAIL: BSG snap-exclude clear adapter\n");
	ret = 1;
	goto cleanup;
    }
    int direct_present = bsg_view_snap_exclude_feature_get(direct, &direct_feature);
    int adapter_present = rt_view_snap_exclude_feature_from_bsg(&adapter_feature, adapter);
    if (direct_present != adapter_present ||
	    direct_feature.token != adapter_feature.token ||
	    direct_feature.revision != adapter_feature.revision ||
	    direct_feature.token != 0 ||
	    direct_feature.revision != 0) {
	printf("FAIL: BSG snap-exclude clear adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    bsg_snap_result_free(&direct_res);
    bsg_snap_result_free(&adapter_res);
    free_view(direct);
    free_view(adapter);
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
    int set_initialized = 0;
    int ret = 0;

    rt_view_init_bsg(NULL, NULL);
    rt_view_free_bsg(NULL);
    rt_view_set_init_bsg(NULL);
    rt_view_set_free_bsg(NULL);
    rt_view_set_add_view_bsg(NULL, NULL);
    rt_view_set_remove_view_bsg(NULL, NULL);

    if (rt_view_set_views_bsg(NULL)) {
	printf("FAIL: null BSG view-set adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_set_init_bsg(&set);
    set_initialized = 1;
    if (rt_view_set_views_bsg(&set) != bsg_set_views(&set) ||
	    BU_PTBL_LEN(rt_view_set_views_bsg(&set)) != 0) {
	printf("FAIL: empty BSG view-set adapter\n");
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

    BU_ALLOC(lifecycle_view, struct bsg_view);
    rt_view_init_bsg(lifecycle_view, NULL);
    bu_vls_sprintf(&lifecycle_view->gv_name, "rt_view_lifecycle_adapter");
    if (lifecycle_view->magic != BSG_VIEW_MAGIC || lifecycle_view->vset) {
	printf("FAIL: BSG view init adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_set_add_view_bsg(&set, lifecycle_view);
    if (lifecycle_view->vset != &set ||
	    rt_view_set_find_view_bsg(&set, "rt_view_lifecycle_adapter") != lifecycle_view) {
	printf("FAIL: BSG view-set add adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_set_remove_view_bsg(&set, lifecycle_view);
    if (lifecycle_view->vset ||
	    rt_view_set_find_view_bsg(&set, "rt_view_lifecycle_adapter")) {
	printf("FAIL: BSG view-set remove adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_selection_bsg(NULL) ||
	    rt_view_selection_const_bsg(NULL)) {
	printf("FAIL: null BSG selection adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_selection_bsg(v) != bsg_view_selection(v) ||
	    rt_view_selection_const_bsg(v) != bsg_view_selection_const(v) ||
	    bsg_selection_count(rt_view_selection_bsg(v)) != 0) {
	printf("FAIL: BSG selection adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_is_independent_bsg(NULL) != bsg_view_is_independent(NULL)) {
	printf("FAIL: null BSG independent adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!bsg_scene_ref_is_null(rt_view_independent_scope_ref_bsg(NULL, 1))) {
	printf("FAIL: null BSG independent-scope ref adapter\n");
	ret = 1;
	goto cleanup;
    }
    rt_view_independent_scope_destroy_bsg(NULL);

    if (bsg_scene_ref_is_null(bsg_view_scene_ref_ensure(v))) {
	printf("FAIL: BSG view scope adapter test scene root setup\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_is_independent_bsg(v) != bsg_view_is_independent(v)) {
	printf("FAIL: BSG shared independent adapter mismatch\n");
	ret = 1;
	goto cleanup;
    }

    bsg_scene_ref scope = rt_view_independent_scope_ref_bsg(v, 1);
    if (bsg_scene_ref_is_null(scope) ||
	    !bsg_scene_ref_equal(scope, bsg_view_independent_scope_ref(v, 0))) {
	printf("FAIL: BSG independent-scope ref adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_is_independent_bsg(v) ||
	    rt_view_is_independent_bsg(v) != bsg_view_is_independent(v)) {
	printf("FAIL: BSG independent adapter mismatch\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_independent_scope_destroy_bsg(v);
    if (rt_view_is_independent_bsg(v)) {
	printf("FAIL: BSG independent-scope destroy adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (rt_view_clear_bsg(v, RT_VIEW_CLEAR_VIEW_BSG) !=
	    bsg_clear(v, RT_VIEW_CLEAR_VIEW_BSG)) {
	printf("FAIL: BSG clear adapter\n");
	ret = 1;
	goto cleanup;
    }

    rt_view_init_copy_bsg(NULL, src, NULL);
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
    rt_view_init_copy_bsg(adapter, src, NULL);

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
    if (set_initialized && lifecycle_view && lifecycle_view->vset == &set)
	rt_view_set_remove_view_bsg(&set, lifecycle_view);
    if (set_initialized && set_view)
	bsg_set_rm_view(&set, set_view);
    if (set_initialized)
	rt_view_set_free_bsg(&set);
    if (lifecycle_view) {
	rt_view_free_bsg(lifecycle_view);
	BU_PUT(lifecycle_view, struct bsg_view);
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

static void
fill_adc_state(struct bsg_adc_state *state, int base)
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

static void
fill_grid_state(struct bsg_grid_state *state, int base)
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

static void
fill_axes_state(struct bsg_axes *state, int base)
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

static void
fill_other_state(struct bsg_other_state *state, int base)
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

static void
fill_params_state(struct bsg_params_state *state, int base)
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
test_bsg_faceplate_state_adapter(void)
{
    struct bsg_view *v = make_view("rt_view_faceplate_state_adapter");
    struct bsg_adc_state zero_adc = {0};
    struct bsg_adc_state source_adc = {0};
    struct bsg_adc_state direct_adc = {0};
    struct bsg_adc_state adapter_adc = {0};
    struct bsg_grid_state zero_grid = {0};
    struct bsg_grid_state source_grid = {0};
    struct bsg_grid_state direct_grid = {0};
    struct bsg_grid_state adapter_grid = {0};
    struct bsg_axes zero_axes = {0};
    struct bsg_axes source_axes = {0};
    struct bsg_axes direct_axes = {0};
    struct bsg_axes adapter_axes = {0};
    struct bsg_other_state zero_other = {0};
    struct bsg_other_state source_other = {0};
    struct bsg_other_state direct_other = {0};
    struct bsg_other_state adapter_other = {0};
    struct bsg_params_state zero_params = {0};
    struct bsg_params_state source_params = {0};
    struct bsg_params_state direct_params = {0};
    struct bsg_params_state adapter_params = {0};
    int ret = 0;

    fill_adc_state(&source_adc, 100);
    memset(&adapter_adc, 0xff, sizeof(adapter_adc));
    if (rt_view_adc_from_bsg(&adapter_adc, NULL) ||
	    memcmp(&adapter_adc, &zero_adc, sizeof(adapter_adc))) {
	printf("FAIL: null BSG ADC get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_adc_set_bsg(NULL, &source_adc) ||
	    rt_view_adc_set_bsg(v, NULL) ||
	    rt_view_adc_from_bsg(NULL, v)) {
	printf("FAIL: null BSG ADC adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_adc_set_bsg(v, &source_adc) ||
	    !bsg_view_adc_get(v, &direct_adc) ||
	    memcmp(&direct_adc, &source_adc, sizeof(direct_adc)) ||
	    !rt_view_adc_from_bsg(&adapter_adc, v) ||
	    memcmp(&adapter_adc, &source_adc, sizeof(adapter_adc))) {
	printf("FAIL: BSG ADC adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_grid_state(&source_grid, 200);
    memset(&adapter_grid, 0xff, sizeof(adapter_grid));
    if (rt_view_grid_from_bsg(&adapter_grid, NULL) ||
	    memcmp(&adapter_grid, &zero_grid, sizeof(adapter_grid))) {
	printf("FAIL: null BSG grid get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_grid_set_bsg(NULL, &source_grid) ||
	    rt_view_grid_set_bsg(v, NULL) ||
	    rt_view_grid_from_bsg(NULL, v)) {
	printf("FAIL: null BSG grid adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_grid_set_bsg(v, &source_grid) ||
	    !bsg_view_grid_get(v, &direct_grid) ||
	    memcmp(&direct_grid, &source_grid, sizeof(direct_grid)) ||
	    !rt_view_grid_from_bsg(&adapter_grid, v) ||
	    memcmp(&adapter_grid, &source_grid, sizeof(adapter_grid))) {
	printf("FAIL: BSG grid adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_axes_state(&source_axes, 300);
    memset(&adapter_axes, 0xff, sizeof(adapter_axes));
    if (rt_view_model_axes_from_bsg(&adapter_axes, NULL) ||
	    memcmp(&adapter_axes, &zero_axes, sizeof(adapter_axes))) {
	printf("FAIL: null BSG model-axes get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_model_axes_set_bsg(NULL, &source_axes) ||
	    rt_view_model_axes_set_bsg(v, NULL) ||
	    rt_view_model_axes_from_bsg(NULL, v)) {
	printf("FAIL: null BSG model-axes adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_model_axes_set_bsg(v, &source_axes) ||
	    !bsg_view_model_axes_get(v, &direct_axes) ||
	    memcmp(&direct_axes, &source_axes, sizeof(direct_axes)) ||
	    !rt_view_model_axes_from_bsg(&adapter_axes, v) ||
	    memcmp(&adapter_axes, &source_axes, sizeof(adapter_axes))) {
	printf("FAIL: BSG model-axes adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_axes_state(&source_axes, 400);
    memset(&adapter_axes, 0xff, sizeof(adapter_axes));
    if (rt_view_view_axes_from_bsg(&adapter_axes, NULL) ||
	    memcmp(&adapter_axes, &zero_axes, sizeof(adapter_axes))) {
	printf("FAIL: null BSG view-axes get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_view_axes_set_bsg(NULL, &source_axes) ||
	    rt_view_view_axes_set_bsg(v, NULL) ||
	    rt_view_view_axes_from_bsg(NULL, v)) {
	printf("FAIL: null BSG view-axes adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_view_axes_set_bsg(v, &source_axes) ||
	    !bsg_view_view_axes_get(v, &direct_axes) ||
	    memcmp(&direct_axes, &source_axes, sizeof(direct_axes)) ||
	    !rt_view_view_axes_from_bsg(&adapter_axes, v) ||
	    memcmp(&adapter_axes, &source_axes, sizeof(adapter_axes))) {
	printf("FAIL: BSG view-axes adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_other_state(&source_other, 10);
    memset(&adapter_other, 0xff, sizeof(adapter_other));
    if (rt_view_center_dot_from_bsg(&adapter_other, NULL) ||
	    memcmp(&adapter_other, &zero_other, sizeof(adapter_other))) {
	printf("FAIL: null BSG center-dot get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_center_dot_set_bsg(NULL, &source_other) ||
	    rt_view_center_dot_set_bsg(v, NULL) ||
	    rt_view_center_dot_from_bsg(NULL, v)) {
	printf("FAIL: null BSG center-dot adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_center_dot_set_bsg(v, &source_other) ||
	    !bsg_view_center_dot_get(v, &direct_other) ||
	    memcmp(&direct_other, &source_other, sizeof(direct_other)) ||
	    !rt_view_center_dot_from_bsg(&adapter_other, v) ||
	    memcmp(&adapter_other, &source_other, sizeof(adapter_other))) {
	printf("FAIL: BSG center-dot adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_other_state(&source_other, 30);
    memset(&adapter_other, 0xff, sizeof(adapter_other));
    if (rt_view_scale_overlay_from_bsg(&adapter_other, NULL) ||
	    memcmp(&adapter_other, &zero_other, sizeof(adapter_other))) {
	printf("FAIL: null BSG scale overlay get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_scale_overlay_set_bsg(NULL, &source_other) ||
	    rt_view_scale_overlay_set_bsg(v, NULL) ||
	    rt_view_scale_overlay_from_bsg(NULL, v)) {
	printf("FAIL: null BSG scale overlay adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_scale_overlay_set_bsg(v, &source_other) ||
	    !bsg_view_scale_state_get(v, &direct_other) ||
	    memcmp(&direct_other, &source_other, sizeof(direct_other)) ||
	    !rt_view_scale_overlay_from_bsg(&adapter_other, v) ||
	    memcmp(&adapter_other, &source_other, sizeof(adapter_other))) {
	printf("FAIL: BSG scale overlay adapter\n");
	ret = 1;
	goto cleanup;
    }

    fill_params_state(&source_params, 50);
    memset(&adapter_params, 0xff, sizeof(adapter_params));
    if (rt_view_params_from_bsg(&adapter_params, NULL) ||
	    memcmp(&adapter_params, &zero_params, sizeof(adapter_params))) {
	printf("FAIL: null BSG params get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_params_set_bsg(NULL, &source_params) ||
	    rt_view_params_set_bsg(v, NULL) ||
	    rt_view_params_from_bsg(NULL, v)) {
	printf("FAIL: null BSG params adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_params_set_bsg(v, &source_params) ||
	    !bsg_view_params_get(v, &direct_params) ||
	    memcmp(&direct_params, &source_params, sizeof(direct_params)) ||
	    !rt_view_params_from_bsg(&adapter_params, v) ||
	    memcmp(&adapter_params, &source_params, sizeof(adapter_params))) {
	printf("FAIL: BSG params adapter\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    free_view(v);
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
	    rt_view_snap_source_flags_set_bsg(NULL, RT_VIEW_SNAP_TCL_BSG)) {
	printf("FAIL: null BSG snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_lines_set_bsg(v, 1) ||
	    !rt_view_snap_source_flags_set_bsg(v, RT_VIEW_SNAP_TCL_BSG) ||
	    !rt_view_snap_lines_from_bsg(v) ||
	    rt_view_snap_source_flags_from_bsg(v) != RT_VIEW_SNAP_TCL_BSG ||
	    !(rt_view_snap_kind_mask_from_bsg(v) & RT_VIEW_SNAP_KIND_ENDPOINT_BSG)) {
	printf("FAIL: BSG snap policy adapter\n");
	ret = 1;
	goto cleanup;
    }

    if (!rt_view_snap_source_flags_set_bsg(v, RT_VIEW_SNAP_DB_BSG) ||
	    !(rt_view_prepare_tcl_snap_bsg(v) & RT_VIEW_SNAP_KIND_ENDPOINT_BSG) ||
	    rt_view_snap_source_flags_from_bsg(v) != RT_VIEW_SNAP_TCL_BSG ||
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

    if (!fastf_equal(rt_view_radius_from_bsg(NULL), 1.0)) {
	printf("FAIL: null BSG view radius adapter\n");
	ret = 1;
	goto cleanup;
    }

    v->radius = 77.0;
    if (!fastf_equal(rt_view_radius_from_bsg(v), 77.0)) {
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
    if (rt_view_screen_to_view_from_bsg(NULL, &adapter_vy,
	    v, 60.0, 70.0) ||
	    !fastf_equal(adapter_vy, 0.0)) {
	printf("FAIL: null output BSG screen-to-view adapter\n");
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

    point_t direct_point = VINIT_ZERO;
    point_t adapter_point = VINIT_ZERO;
    if (rt_view_screen_point_from_bsg(adapter_point, NULL, 60.0, 70.0) ||
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

    point_t current_point = VINIT_ZERO;
    if (rt_view_current_point_from_bsg(current_point, NULL) ||
	    !fastf_equal(current_point[X], 0.0) ||
	    !fastf_equal(current_point[Y], 0.0) ||
	    !fastf_equal(current_point[Z], 0.0) ||
	    rt_view_current_point_from_bsg(NULL, v)) {
	printf("FAIL: null BSG current-point adapter\n");
	ret = 1;
	goto cleanup;
    }
    VSET(v->gv_point, 1.0, 2.0, 3.0);
    if (!rt_view_current_point_from_bsg(current_point, v) ||
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
	    !rt_view_current_point_set_bsg(v, updated_point) ||
	    !VNEAR_EQUAL(v->gv_point, updated_point, BN_TOL_DIST)) {
	printf("FAIL: BSG current-point set adapter\n");
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

    struct rt_view_mouse_delta_settings mouse_settings = RT_VIEW_MOUSE_DELTA_SETTINGS_INIT;
    if (rt_view_mouse_delta_settings_from_bsg(NULL, v) ||
	    rt_view_mouse_delta_settings_from_bsg(&mouse_settings, NULL) ||
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

    struct bu_vls direct_name = BU_VLS_INIT_ZERO;
    struct bu_vls adapter_name = BU_VLS_INIT_ZERO;
    if (rt_view_unique_object_name_bsg(NULL, "view_poly_0", v) ||
	    rt_view_unique_object_name_bsg(&adapter_name, "view_poly_0", NULL)) {
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
    bu_vls_free(&direct_name);
    bu_vls_free(&adapter_name);

    struct bsg_interactive_rect_state zero_rect = {0};
    struct bsg_interactive_rect_state direct_rect;
    struct bsg_interactive_rect_state adapter_rect;
    struct bsg_interactive_rect_state updated_rect = {0};
    updated_rect.active = 1;
    updated_rect.draw = 1;
    updated_rect.line_width = 3;
    updated_rect.line_style = 1;
    updated_rect.pos[0] = 11;
    updated_rect.pos[1] = 12;
    updated_rect.dim[0] = 101;
    updated_rect.dim[1] = 102;
    updated_rect.x = -0.25;
    updated_rect.y = 0.35;
    updated_rect.width = 0.5;
    updated_rect.height = 0.4;
    updated_rect.bg[0] = 1;
    updated_rect.bg[1] = 2;
    updated_rect.bg[2] = 3;
    updated_rect.color[0] = 4;
    updated_rect.color[1] = 5;
    updated_rect.color[2] = 6;
    updated_rect.cdim[0] = 400;
    updated_rect.cdim[1] = 200;
    updated_rect.aspect = 2.0;
    if (rt_view_interactive_rect_from_bsg(&adapter_rect, NULL) ||
	    memcmp(&adapter_rect, &zero_rect, sizeof(adapter_rect))) {
	printf("FAIL: null BSG interactive rectangle get adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (rt_view_interactive_rect_set_bsg(NULL, &updated_rect) ||
	    rt_view_interactive_rect_set_bsg(v, NULL) ||
	    rt_view_interactive_rect_from_bsg(NULL, v)) {
	printf("FAIL: null BSG interactive rectangle adapter arguments\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_interactive_rect_set_bsg(v, &updated_rect) ||
	    !bsg_view_interactive_rect_get(v, &direct_rect) ||
	    memcmp(&direct_rect, &updated_rect, sizeof(direct_rect)) ||
	    !rt_view_interactive_rect_from_bsg(&adapter_rect, v) ||
	    memcmp(&adapter_rect, &updated_rect, sizeof(adapter_rect))) {
	printf("FAIL: BSG interactive rectangle adapter\n");
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
	    rt_view_refresh_dirty_from_bsg(NULL) ||
	    rt_view_refresh_consume_bsg(NULL) ||
	    rt_view_refresh_complete_bsg(NULL) ||
	    rt_view_refresh_enabled_from_bsg(NULL) != 1 ||
	    rt_view_refresh_enabled_set_bsg(NULL, 1) ||
	    rt_view_refresh_suppressed_from_bsg(NULL) ||
	    rt_view_refresh_suppress_begin_bsg(NULL) ||
	    rt_view_refresh_suppress_end_bsg(NULL) ||
	    rt_view_refresh_drawn_count_from_bsg(NULL) ||
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
    if (!rt_view_refresh_suppress_begin_bsg(v) ||
	    !rt_view_refresh_suppressed_from_bsg(v) ||
	    !rt_view_refresh_request_bsg(v, RT_VIEW_REFRESH_DRAW_BSG) ||
	    rt_view_refresh_dirty_from_bsg(v) ||
	    rt_view_refresh_consume_bsg(v) ||
	    !rt_view_refresh_suppress_end_bsg(v) ||
	    !rt_view_refresh_dirty_from_bsg(v) ||
	    !rt_view_refresh_complete_bsg(v) ||
	    rt_view_refresh_suppressed_from_bsg(v)) {
	printf("FAIL: BSG refresh suppress adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_refresh_enabled_set_bsg(v, 0) ||
	    rt_view_refresh_enabled_from_bsg(v) ||
	    !rt_view_refresh_request_bsg(v, RT_VIEW_REFRESH_EDIT_BSG) ||
	    rt_view_refresh_dirty_from_bsg(v) ||
	    rt_view_refresh_consume_bsg(v) ||
	    !rt_view_refresh_enabled_set_bsg(v, 1) ||
	    !rt_view_refresh_dirty_from_bsg(v) ||
	    !rt_view_refresh_complete_bsg(v)) {
	printf("FAIL: BSG refresh enabled adapter\n");
	ret = 1;
	goto cleanup;
    }
    if (!rt_view_refresh_drawn_count_set_bsg(v, 3) ||
	    rt_view_refresh_drawn_count_from_bsg(v) != 3 ||
	    !rt_view_refresh_drawn_count_set_bsg(v, -1) ||
	    rt_view_refresh_drawn_count_from_bsg(v) != 0) {
	printf("FAIL: BSG refresh drawn-count adapter\n");
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
    if (test_bsg_view_scope_adapter())
	return 1;
    if (test_bsg_perspective_adapter())
	return 1;
    if (test_bsg_camera_adapter())
	return 1;
    if (test_bsg_lod_policy_adapter())
	return 1;
    if (test_bsg_faceplate_state_adapter())
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
