/*                       V I E W _ I N F O . C
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
/** @file librt/tests/view_info.c
 *
 * Unit tests for RT-owned view snapshots.
 */

#include "common.h"

#include <math.h>
#include <stdio.h>

#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "wdb.h"
#include "rt/view.h"
#include "rt/wdb.h"

static int
fastf_equal(fastf_t a, fastf_t b)
{
    return fabs(a - b) <= SMALL_FASTF;
}

static int
fastf_near(fastf_t a, fastf_t b, fastf_t tol)
{
    return fabs(a - b) <= tol;
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
test_defaults(void)
{
    struct rt_view_info info;

    info.width = 23;
    info.height = 42;
    info.size = 99.0;
    info.lod.scale = 7.0;
    info.lod.curve_scale = 8.0;
    info.lod.point_scale = 9.0;
    info.lod.bot_threshold = 100;

    rt_view_info_init(NULL);
    rt_view_info_init(&info);

    return check_view_info("defaults", &info, 1, 1, 1.0, 1.0, 1.0, 1.0, 0);
}

static int
test_sanitize(void)
{
    struct rt_view_info info;

    info.width = 0;
    info.height = -10;
    info.size = 0.0;
    info.lod.scale = 0.0;
    info.lod.curve_scale = -2.0;
    info.lod.point_scale = 0.0;
    info.lod.bot_threshold = 33;

    rt_view_info_sanitize(NULL);
    rt_view_info_sanitize(&info);

    return check_view_info("sanitize", &info, 1, 1, 1.0, 1.0, 1.0, 1.0, 33);
}

static int
test_lod_policy_defaults(void)
{
    struct rt_view_lod_policy policy;

    policy.policy = RT_VIEW_LOD_OFF;
    policy.forced_level = 9;
    policy.mesh_enabled = 1;
    policy.csg_enabled = 1;
    policy.zoom_refresh = 1;
    policy.bot_threshold = 111;
    policy.scale = 7.0;
    policy.curve_scale = 8.0;
    policy.point_scale = 9.0;

    rt_view_lod_policy_init(NULL);
    rt_view_lod_policy_init(&policy);
    if (policy.policy != RT_VIEW_LOD_AUTO ||
	    policy.forced_level != 0 ||
	    policy.mesh_enabled != 0 ||
	    policy.csg_enabled != 0 ||
	    policy.zoom_refresh != 0 ||
	    policy.bot_threshold != 0 ||
	    !fastf_equal(policy.scale, 1.0) ||
	    !fastf_equal(policy.curve_scale, 1.0) ||
	    !fastf_equal(policy.point_scale, 1.0)) {
	printf("FAIL: lod policy defaults\n");
	return 1;
    }

    policy.scale = 0.0;
    policy.curve_scale = -2.0;
    policy.point_scale = 0.0;
    policy.bot_threshold = 222;
    rt_view_lod_policy_sanitize(NULL);
    rt_view_lod_policy_sanitize(&policy);
    if (!fastf_equal(policy.scale, 1.0) ||
	    !fastf_equal(policy.curve_scale, 1.0) ||
	    !fastf_equal(policy.point_scale, 1.0) ||
	    policy.bot_threshold != 222) {
	printf("FAIL: lod policy sanitize\n");
	return 1;
    }

    return 0;
}

static int
test_view_metrics(void)
{
    struct rt_view_info info = RT_VIEW_INFO_INIT;
    struct rt_view_info bad_policy = RT_VIEW_INFO_INIT;
    fastf_t point_spacing;

    info.width = 200;
    info.height = 100;
    info.size = 20.0;
    info.lod.curve_scale = 3.0;
    info.lod.point_scale = 2.0;
    info.lod.bot_threshold = 42;

    if (!fastf_equal(rt_view_lod_curve_scale(&info), 3.0) ||
	    rt_view_lod_bot_threshold(&info) != 42) {
	printf("FAIL: view metric lod policy\n");
	return 1;
    }

    if (!fastf_equal(rt_view_avg_sample_spacing(&info), 0.1)) {
	printf("FAIL: view average sample spacing: got %g\n",
		rt_view_avg_sample_spacing(&info));
	return 1;
    }

    point_spacing = rt_view_solid_point_spacing(&info, 4.0);
    if (!fastf_near(point_spacing, sqrt(0.2) / 2.0, 1.0e-6)) {
	printf("FAIL: view solid point spacing: got %g\n", point_spacing);
	return 1;
    }

    bad_policy.lod.curve_scale = 0.0;
    bad_policy.lod.point_scale = -2.0;
    bad_policy.lod.scale = 0.0;
    bad_policy.lod.bot_threshold = 55;
    if (!fastf_equal(rt_view_lod_curve_scale(&bad_policy), 1.0) ||
	    rt_view_lod_bot_threshold(&bad_policy) != 55) {
	printf("FAIL: view metric bad lod policy\n");
	return 1;
    }

    if (!fastf_equal(rt_view_lod_curve_scale(NULL), 1.0) ||
	    rt_view_lod_bot_threshold(NULL) != 0 ||
	    !fastf_equal(rt_view_avg_sample_spacing(NULL), 1.0)) {
	printf("FAIL: view metric null defaults\n");
	return 1;
    }

    return 0;
}

static int
check_mesh_lod_payload(const char *label, struct rt_mesh_lod *lod,
	size_t max_face_count, size_t max_point_count)
{
    struct rt_mesh_lod_data data;
    struct rt_mesh_lod_info info;

    if (!rt_mesh_lod_data_get(lod, &data) ||
	    data.face_count == 0 || data.face_count > max_face_count ||
	    data.point_count == 0 || data.point_count > max_point_count ||
	    data.point_orig_count == 0 || data.point_orig_count > max_point_count ||
	    !data.faces || !data.points || !data.points_orig) {
	printf("FAIL: %s mesh lod data faces=%zu points=%zu orig=%zu ptrs=%d/%d/%d\n",
		label, data.face_count, data.point_count, data.point_orig_count,
		data.faces ? 1 : 0, data.points ? 1 : 0,
		data.points_orig ? 1 : 0);
	return 1;
    }

    if (!rt_mesh_lod_has_active_data(lod)) {
	printf("FAIL: %s mesh lod active data status\n", label);
	return 1;
    }

    if (!rt_mesh_lod_info_get(lod, &info) ||
	    info.active_level < 0 ||
	    info.face_count != data.face_count ||
	    info.point_count != data.point_count ||
	    info.point_orig_count != data.point_orig_count ||
	    info.normal_count != data.normal_count ||
	    !info.has_faces || !info.has_points || !info.has_original_points ||
	    info.has_normals != (data.normals ? 1 : 0)) {
	printf("FAIL: %s mesh lod info\n", label);
	return 1;
    }

    return 0;
}

static int
test_mesh_lod_api(void)
{
    const char *objname = "rt_view_info_lod_bot";
    const char *mesh_objname = "rt_view_info_lod_mesh_payload";
    const int grid = 12;
    const int vertex_count = (grid + 1) * (grid + 1);
    const int face_count = grid * grid * 2;
    char dbpath[MAXPATHLEN] = {0};
    fastf_t *vertices = NULL;
    int *faces = NULL;
    int ret = 0;
    char cache_dir[MAXPATHLEN] = {0};
    struct db_i *dbip = NULL;
    struct rt_mesh_lod *lod = NULL;
    int memshrink_level = -1;
    struct rt_mesh_lod_cache_status cache_status = RT_MESH_LOD_CACHE_STATUS_INIT;

    bu_dir(cache_dir, MAXPATHLEN, BU_DIR_CURR, "rt_view_info_lod_cache", NULL);
    bu_mkdir(cache_dir);
    bu_setenv("BU_DIR_CACHE", cache_dir, 1);

    vertices = (fastf_t *)bu_calloc((size_t)vertex_count * 3,
	    sizeof(fastf_t), "view_info lod vertices");
    faces = (int *)bu_calloc((size_t)face_count * 3,
	    sizeof(int), "view_info lod faces");
    for (int y = 0; y <= grid; y++) {
	for (int x = 0; x <= grid; x++) {
	    int idx = y * (grid + 1) + x;
	    vertices[idx * 3 + 0] = (fastf_t)x;
	    vertices[idx * 3 + 1] = (fastf_t)y;
	    vertices[idx * 3 + 2] = (fastf_t)((x + y) % 3) * 0.05;
	}
    }
    {
	int fi = 0;
	for (int y = 0; y < grid; y++) {
	    for (int x = 0; x < grid; x++) {
		int v00 = y * (grid + 1) + x;
		int v10 = v00 + 1;
		int v01 = v00 + (grid + 1);
		int v11 = v01 + 1;
		faces[fi++] = v00;
		faces[fi++] = v10;
		faces[fi++] = v11;
		faces[fi++] = v00;
		faces[fi++] = v11;
		faces[fi++] = v01;
	    }
	}
    }

    FILE *fp = bu_temp_file(dbpath, MAXPATHLEN);
    if (!fp) {
	printf("FAIL: mesh lod temp file\n");
	ret = 1;
	goto cleanup;
    }
    fclose(fp);

    dbip = db_create(dbpath, 5);
    if (!dbip) {
	printf("FAIL: mesh lod db_create\n");
	ret = 1;
	goto cleanup;
    }

    struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
    if (!wdbp) {
	printf("FAIL: mesh lod wdb_dbopen\n");
	ret = 1;
	goto cleanup;
    }

    if (mk_bot(wdbp, objname, RT_BOT_SURFACE, RT_BOT_CCW, 0,
		vertex_count, face_count,
		vertices, faces, NULL, NULL) < 0) {
	printf("FAIL: mesh lod mk_bot\n");
	ret = 1;
	goto cleanup;
    }
    {
	point_t center = VINIT_ZERO;
	if (mk_sph(wdbp, mesh_objname, center, 1.0) < 0) {
	    printf("FAIL: mesh lod mk_sph\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    if (db_mesh_lod_status(dbip, objname, &cache_status) != BRLCAD_OK ||
	    !cache_status.directory_found || !cache_status.is_bot ||
	    cache_status.has_cache_key || cache_status.has_cached_payload ||
	    cache_status.stale_cache_entry) {
	printf("FAIL: mesh lod initial cache status\n");
	ret = 1;
	goto cleanup;
    }

    if (db_mesh_lod_refresh(dbip, objname, &cache_status) != BRLCAD_OK ||
	    !cache_status.directory_found || !cache_status.is_bot ||
	    !cache_status.has_cache_key || !cache_status.has_cached_payload ||
	    cache_status.stale_cache_entry ||
	    !cache_status.generated_cache_entry ||
	    !cache_status.cache_key) {
	printf("FAIL: mesh lod refresh status\n");
	ret = 1;
	goto cleanup;
    }

    if (db_mesh_lod_status(dbip, objname, &cache_status) != BRLCAD_OK ||
	    !cache_status.has_cache_key || !cache_status.has_cached_payload ||
	    cache_status.stale_cache_entry ||
	    cache_status.generated_cache_entry ||
	    !cache_status.cache_key) {
	printf("FAIL: mesh lod post-refresh status\n");
	ret = 1;
	goto cleanup;
    }

    lod = db_mesh_lod_get(dbip, objname);
    if (!lod) {
	printf("FAIL: mesh lod get\n");
	ret = 1;
	goto cleanup;
    }

    if (!ret) {
	int got_level_data = 0;
	for (int level = 0; level < 16; level++) {
	    struct rt_mesh_lod_data data;
	    if (rt_mesh_lod_load_level(lod, level, 0) < 0) {
		printf("FAIL: mesh lod load level\n");
		ret = 1;
		break;
	    }
	    if (!rt_mesh_lod_data_get(lod, &data))
		continue;
	    if (data.face_count > (size_t)face_count ||
		    data.point_count > (size_t)vertex_count ||
		    data.point_orig_count > (size_t)vertex_count) {
		printf("FAIL: mesh lod level data counts\n");
		ret = 1;
		break;
	    }
	    got_level_data = 1;
	    memshrink_level = level;
	    break;
	}
	if (!ret && !got_level_data) {
	    printf("FAIL: mesh lod level data\n");
	    ret = 1;
	}
    }

    if (!ret) {
	struct rt_view_info info = RT_VIEW_INFO_INIT;
	info.size = 0.01;
	info.lod.scale = 1.0;
	if (rt_mesh_lod_load_view(lod, &info, 0) < 0) {
	    printf("FAIL: mesh lod load view\n");
	    ret = 1;
	} else if (check_mesh_lod_payload("mesh lod view", lod,
		(size_t)face_count, (size_t)vertex_count)) {
	    ret = 1;
	}
    }

    if (!ret) {
	struct rt_mesh_lod *mesh_lod = NULL;
	if (db_mesh_lod_store_mesh(dbip, mesh_objname,
		(const point_t *)vertices, (size_t)vertex_count, NULL, faces,
		(size_t)face_count, 424242ULL, 0.66, &cache_status) != BRLCAD_OK ||
		!cache_status.directory_found || cache_status.is_bot ||
		!cache_status.has_cache_key || !cache_status.has_cached_payload ||
		cache_status.stale_cache_entry ||
		!cache_status.generated_cache_entry ||
		!cache_status.cache_key) {
	    printf("FAIL: mesh lod store generated mesh status\n");
	    ret = 1;
	} else if (!(mesh_lod = db_mesh_lod_get(dbip, mesh_objname))) {
	    printf("FAIL: mesh lod get generated mesh\n");
	    ret = 1;
	} else {
	    int got_generated_mesh_data = 0;
	    for (int level = 0; level < 16; level++) {
		struct rt_mesh_lod_data data;
		if (rt_mesh_lod_load_level(mesh_lod, level, 0) < 0)
		    continue;
		if (!rt_mesh_lod_data_get(mesh_lod, &data))
		    continue;
		if (check_mesh_lod_payload("mesh lod generated mesh", mesh_lod,
			(size_t)face_count, (size_t)vertex_count)) {
		    ret = 1;
		    break;
		}
		got_generated_mesh_data = 1;
		break;
	    }
	    if (!got_generated_mesh_data) {
		printf("FAIL: mesh lod generated mesh data\n");
		ret = 1;
	    }
	}
	if (mesh_lod)
	    rt_mesh_lod_destroy(mesh_lod);
    }

    if (!ret && memshrink_level >= 0) {
	int active_level;
	struct rt_mesh_lod_data shrink_data;
	struct rt_mesh_lod_info shrink_info = RT_MESH_LOD_INFO_INIT;
	if (rt_mesh_lod_load_level(lod, memshrink_level, 0) != memshrink_level) {
	    printf("FAIL: mesh lod load memshrink level\n");
	    ret = 1;
	} else if ((active_level = rt_mesh_lod_current_level(lod)) < 0) {
	    printf("FAIL: mesh lod current level\n");
	    ret = 1;
	} else if (!rt_mesh_lod_has_active_data(lod)) {
	    printf("FAIL: mesh lod active data status before memshrink\n");
	    ret = 1;
	} else {
	    rt_mesh_lod_memshrink(lod);
	    if (rt_mesh_lod_current_level(lod) != active_level ||
		    rt_mesh_lod_has_active_data(lod) ||
		    rt_mesh_lod_data_get(lod, &shrink_data) ||
		    rt_mesh_lod_info_get(lod, &shrink_info) ||
		    shrink_info.active_level != active_level) {
		printf("FAIL: mesh lod memshrink status\n");
		ret = 1;
	    } else if (rt_mesh_lod_load_level(lod, active_level, 0) != active_level ||
		    check_mesh_lod_payload("mesh lod reload after memshrink", lod,
			(size_t)face_count, (size_t)vertex_count)) {
		printf("FAIL: mesh lod reload after memshrink active=%d\n",
			active_level);
		ret = 1;
	    }
	}
    }

    if (!ret) {
	rt_mesh_lod_destroy(lod);
	lod = NULL;
	db_close(dbip);
	dbip = NULL;

	dbip = db_open(dbpath, DB_OPEN_READWRITE);
	if (!dbip) {
	    printf("FAIL: mesh lod reopen\n");
	    ret = 1;
	    goto cleanup;
	}
	if (db_dirbuild(dbip) < 0) {
	    printf("FAIL: mesh lod reopen dirbuild\n");
	    ret = 1;
	    goto cleanup;
	}

	lod = db_mesh_lod_get(dbip, objname);
	if (!lod) {
	    printf("FAIL: mesh lod get after reopen\n");
	    ret = 1;
	    goto cleanup;
	}
	if (rt_mesh_lod_load_view(lod, NULL, 0) < 0 ||
		check_mesh_lod_payload("mesh lod reopen", lod,
		    (size_t)face_count, (size_t)vertex_count)) {
	    printf("FAIL: mesh lod reopen data\n");
	    ret = 1;
	}
    }

    if (!ret) {
	rt_mesh_lod_destroy(lod);
	lod = NULL;

	if (db_mesh_lod_invalidate(dbip, objname, &cache_status) != BRLCAD_OK ||
		!cache_status.directory_found || !cache_status.is_bot ||
		!cache_status.cleared_cache_entry ||
		!cache_status.cleared_cache_key ||
		cache_status.has_cache_key ||
		cache_status.has_cached_payload ||
		cache_status.stale_cache_entry) {
	    printf("FAIL: mesh lod invalidate status\n");
	    ret = 1;
	    goto cleanup;
	}

	lod = db_mesh_lod_get(dbip, objname);
	if (lod) {
	    printf("FAIL: mesh lod get after invalidate\n");
	    ret = 1;
	    goto cleanup;
	}

	if (db_mesh_lod_status(dbip, objname, &cache_status) != BRLCAD_OK ||
		cache_status.has_cache_key ||
		cache_status.has_cached_payload ||
		cache_status.stale_cache_entry ||
		cache_status.cleared_cache_entry ||
		cache_status.generated_cache_entry) {
	    printf("FAIL: mesh lod status after invalidate\n");
	    ret = 1;
	    goto cleanup;
	}

	if (db_mesh_lod_update(dbip, objname) != BRLCAD_OK) {
	    printf("FAIL: mesh lod compatibility update\n");
	    ret = 1;
	    goto cleanup;
	}
	lod = db_mesh_lod_get(dbip, objname);
	if (!lod || rt_mesh_lod_load_view(lod, NULL, 0) < 0 ||
		check_mesh_lod_payload("mesh lod update after invalidate", lod,
		    (size_t)face_count, (size_t)vertex_count)) {
	    printf("FAIL: mesh lod compatibility update data\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    struct rt_mesh_lod_data null_data;
    struct rt_mesh_lod_detail null_detail;
    struct rt_mesh_lod_info null_info;
    null_info.active_level = 7;
    null_detail.face_count = 42;
    rt_mesh_lod_cache_status_init(&cache_status);
    cache_status.cache_key = 99;
    rt_mesh_lod_cache_status_init(&cache_status);
    rt_mesh_lod_detail_init(&null_detail);
    if (rt_mesh_lod_load_level(NULL, 0, 0) >= 0 ||
	    rt_mesh_lod_load_view(NULL, NULL, 0) >= 0 ||
	    rt_mesh_lod_current_level(NULL) >= 0 ||
	    rt_mesh_lod_has_active_data(NULL) ||
	    rt_mesh_lod_data_get(NULL, &null_data) ||
	    rt_mesh_lod_info_get(NULL, &null_info) ||
	    null_detail.face_count ||
	    rt_mesh_lod_detail_callbacks_set(NULL, NULL, NULL, NULL, NULL) ||
	    null_info.active_level != -1 ||
	    cache_status.cache_key ||
	    db_mesh_lod_status(NULL, objname, &cache_status) != BRLCAD_ERROR ||
	    db_mesh_lod_status(dbip, NULL, &cache_status) != BRLCAD_ERROR ||
	    db_mesh_lod_status(dbip, objname, NULL) != BRLCAD_ERROR ||
	    db_mesh_lod_invalidate(NULL, objname, &cache_status) != BRLCAD_ERROR) {
	printf("FAIL: mesh lod null handling\n");
	ret = 1;
    }
    rt_mesh_lod_memshrink(NULL);

cleanup:
    if (lod)
	rt_mesh_lod_destroy(lod);
    if (dbip) {
	db_mesh_lod_clear(dbip);
	db_close(dbip);
    }
    if (dbpath[0])
	bu_file_delete(dbpath);
    bu_dirclear(cache_dir);
    if (vertices)
	bu_free(vertices, "view_info lod vertices");
    if (faces)
	bu_free(faces, "view_info lod faces");
    return ret;
}

int
main(int argc, char *argv[])
{
    bu_setprogname(argv[0]);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    if (test_defaults())
	return 1;
    if (test_sanitize())
	return 1;
    if (test_lod_policy_defaults())
	return 1;
    if (test_view_metrics())
	return 1;
    if (test_mesh_lod_api())
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
