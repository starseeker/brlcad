/*                 T E S T _ M E S H _ L O D _ C A C H E . C P P
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

#include "common.h"

#include "brlobol/mesh_lod_cache.h"

#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "wdb.h"
#include "rt/wdb.h"

#include <cmath>
#include <cstdio>
#include <cstring>

static int
fastf_equal(fastf_t a, fastf_t b)
{
    return std::fabs(a - b) <= SMALL_FASTF;
}

static int
check_mesh_lod_payload(const char *label,
		       struct BRLObolMeshLod *lod,
		       size_t maxFaceCount,
		       size_t maxPointCount,
		       int requireNormals)
{
    struct BRLObolMeshLodData data;
    struct BRLObolMeshLodInfo info = BRLOBOL_MESH_LOD_INFO_INIT;

    if (!brlobol_mesh_lod_data_get(lod, &data) ||
	data.face_count == 0 || data.face_count > maxFaceCount ||
	data.point_count == 0 || data.point_count > maxPointCount ||
	data.point_orig_count == 0 ||
	data.point_orig_count > maxPointCount ||
	!data.faces || !data.points || !data.points_orig) {
	printf("FAIL: %s mesh lod data faces=%zu points=%zu orig=%zu ptrs=%d/%d/%d\n",
	       label, data.face_count, data.point_count, data.point_orig_count,
	       data.faces ? 1 : 0, data.points ? 1 : 0,
	       data.points_orig ? 1 : 0);
	return 1;
    }

    if (!brlobol_mesh_lod_has_active_data(lod)) {
	printf("FAIL: %s mesh lod active data status\n", label);
	return 1;
    }

    if (!brlobol_mesh_lod_info_get(lod, &info) ||
	info.active_level < 0 ||
	info.face_count != data.face_count ||
	info.point_count != data.point_count ||
	info.point_orig_count != data.point_orig_count ||
	info.normal_count != data.normal_count ||
	!info.has_faces || !info.has_points ||
	!info.has_original_points ||
	info.has_normals != (data.normals ? 1 : 0)) {
	printf("FAIL: %s mesh lod info\n", label);
	return 1;
    }

    if (requireNormals &&
	(!data.normals || data.normal_count != data.face_count * 3 ||
	 !fastf_equal(data.normals[0][Z], 1.0))) {
	printf("FAIL: %s mesh lod normals count=%zu\n",
	       label, data.normal_count);
	return 1;
    }

    return 0;
}

struct mesh_lod_detail_test_data {
    const fastf_t *vertices;
    const int *faces;
    const fastf_t *normals;
    size_t point_count;
    size_t point_orig_count;
    size_t face_count;
    size_t normal_count;
    int setup_count;
    int setup_return;
    int clear_count;
    int free_count;
};

static int
mesh_lod_detail_setup_cb(struct BRLObolMeshLodDetail *detail, void *cbData)
{
    struct mesh_lod_detail_test_data *data =
	    static_cast<struct mesh_lod_detail_test_data *>(cbData);

    if (!detail || !data)
	return -1;

    brlobol_mesh_lod_detail_init(detail);
    detail->faces = data->faces;
    detail->face_count = data->face_count;
    detail->points = reinterpret_cast<const point_t *>(data->vertices);
    detail->point_count = data->point_count;
    detail->points_orig = reinterpret_cast<const point_t *>(data->vertices);
    detail->point_orig_count =
	data->point_orig_count ? data->point_orig_count : data->point_count;
    detail->normals = reinterpret_cast<const vect_t *>(data->normals);
    detail->normal_count = data->normal_count;
    data->setup_count++;
    return data->setup_return;
}

static int
mesh_lod_detail_clear_cb(void *cbData)
{
    struct mesh_lod_detail_test_data *data =
	    static_cast<struct mesh_lod_detail_test_data *>(cbData);

    if (data)
	data->clear_count++;
    return 0;
}

static int
mesh_lod_detail_free_cb(void *cbData)
{
    struct mesh_lod_detail_test_data *data =
	    static_cast<struct mesh_lod_detail_test_data *>(cbData);

    if (data)
	data->free_count++;
    return 0;
}

static int
first_available_level(struct BRLObolMeshLod *lod)
{
    for (int level = 0; level < 16; level++) {
	if (brlobol_mesh_lod_load_level(lod, level, 0) < 0)
	    return -1;
	if (brlobol_mesh_lod_has_active_data(lod))
	    return level;
    }
    return -1;
}

static int
test_detail_callbacks(struct db_i *dbip,
		      const char *name,
		      const fastf_t *vertices,
		      const int *faces,
		      const fastf_t *normals,
		      size_t vertexCount,
		      size_t faceCount)
{
    struct BRLObolMeshLodData detailData;
    struct BRLObolMeshLodInfo detailInfo = BRLOBOL_MESH_LOD_INFO_INIT;
    struct BRLObolMeshLod *lod = brlobol_mesh_lod_get(dbip, name);
    struct mesh_lod_detail_test_data original;
    struct mesh_lod_detail_test_data replacement;
    struct mesh_lod_detail_test_data failing;
    int ret = 0;

    std::memset(&original, 0, sizeof(original));
    std::memset(&replacement, 0, sizeof(replacement));
    std::memset(&failing, 0, sizeof(failing));

    original.vertices = vertices;
    original.faces = faces;
    original.normals = normals;
    original.point_count = vertexCount;
    original.face_count = faceCount;
    original.normal_count = faceCount * 3;

    replacement = original;
    failing = original;
    failing.setup_return = -1;

    if (!lod || first_available_level(lod) < 0 ||
	!brlobol_mesh_lod_detail_callbacks_set(lod,
		mesh_lod_detail_setup_cb, mesh_lod_detail_clear_cb,
		mesh_lod_detail_free_cb, &original) ||
	!brlobol_mesh_lod_detail_callbacks_set(lod,
		mesh_lod_detail_setup_cb, mesh_lod_detail_clear_cb,
		mesh_lod_detail_free_cb, &replacement) ||
	original.free_count != 1 ||
	!brlobol_mesh_lod_has_active_data(lod) ||
	!brlobol_mesh_lod_data_get(lod, &detailData) ||
	replacement.setup_count != 0) {
	printf("FAIL: mesh lod detail callback replacement preserved POP data\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_mesh_lod_load_level(lod, 100, 0) < 0 ||
	!brlobol_mesh_lod_data_get(lod, &detailData) ||
	!brlobol_mesh_lod_info_get(lod, &detailInfo) ||
	detailData.normal_count != faceCount * 3 ||
	detailInfo.normal_count != detailData.normal_count ||
	!detailInfo.has_normals || !detailData.normals ||
	replacement.setup_count != 1) {
	printf("FAIL: mesh lod detail callback full-detail payload\n");
	ret = 1;
	goto cleanup;
    }

    if (!brlobol_mesh_lod_detail_callbacks_set(lod,
	    mesh_lod_detail_setup_cb, mesh_lod_detail_clear_cb,
	    mesh_lod_detail_free_cb, &failing) ||
	replacement.free_count != 1 ||
	brlobol_mesh_lod_current_level(lod) != -1 ||
	brlobol_mesh_lod_has_active_data(lod) ||
	brlobol_mesh_lod_data_get(lod, &detailData)) {
	printf("FAIL: mesh lod detail callback replacement invalidation\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_mesh_lod_load_level(lod, 100, 1) >= 0 ||
	brlobol_mesh_lod_current_level(lod) != -1 ||
	brlobol_mesh_lod_has_active_data(lod) ||
	brlobol_mesh_lod_data_get(lod, &detailData) ||
	failing.setup_count != 1 ||
	failing.clear_count < 1) {
	printf("FAIL: mesh lod detail callback failure clearing\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    if (lod)
	brlobol_mesh_lod_destroy(lod);
    if (!ret && failing.free_count != 1) {
	printf("FAIL: mesh lod detail callback free\n");
	ret = 1;
    }
    return ret;
}

int
main(int argc, char *argv[])
{
    const char *objname = "brlobol_lod_bot";
    const char *meshObjname = "brlobol_lod_mesh_payload";
    const char *invalidBotObjname = "brlobol_invalid_lod_bot";
    const int grid = 12;
    const int vertexCount = (grid + 1) * (grid + 1);
    const int faceCount = grid * grid * 2;
    char dbpath[MAXPATHLEN] = {0};
    char cacheDir[MAXPATHLEN] = {0};
    fastf_t *vertices = NULL;
    fastf_t *detailNormals = NULL;
    int *faces = NULL;
    struct db_i *dbip = NULL;
    struct BRLObolMeshLod *lod = NULL;
    int memshrinkLevel = -1;
    int ret = 0;
    struct BRLObolMeshLodCacheStatus cacheStatus =
	    BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;

    bu_setprogname(argv[0]);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    bu_dir(cacheDir, MAXPATHLEN, BU_DIR_CURR,
	   "brlobol_mesh_lod_cache_test", NULL);
    bu_mkdir(cacheDir);
    bu_setenv("BU_DIR_CACHE", cacheDir, 1);

    vertices = static_cast<fastf_t *>(bu_calloc(
					  static_cast<size_t>(vertexCount) * 3, sizeof(fastf_t),
					  "mesh lod vertices"));
    faces = static_cast<int *>(bu_calloc(
				   static_cast<size_t>(faceCount) * 3, sizeof(int),
				   "mesh lod faces"));
    detailNormals = static_cast<fastf_t *>(bu_calloc(
	    static_cast<size_t>(faceCount) * 3 * 3, sizeof(fastf_t),
	    "mesh lod detail normals"));

    for (int y = 0; y <= grid; y++) {
	for (int x = 0; x <= grid; x++) {
	    int idx = y * (grid + 1) + x;
	    vertices[idx * 3 + X] = static_cast<fastf_t>(x);
	    vertices[idx * 3 + Y] = static_cast<fastf_t>(y);
	    vertices[idx * 3 + Z] =
		static_cast<fastf_t>((x + y) % 3) * 0.05;
	}
    }

    for (int y = 0, fi = 0; y < grid; y++) {
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

    for (size_t ni = 0; ni < static_cast<size_t>(faceCount) * 3; ni++) {
	detailNormals[ni * 3 + X] = 0.0;
	detailNormals[ni * 3 + Y] = 0.0;
	detailNormals[ni * 3 + Z] = 1.0;
    }

    {
	FILE *fp = bu_temp_file(dbpath, MAXPATHLEN);
	if (!fp) {
	    printf("FAIL: mesh lod temp file\n");
	    ret = 1;
	    goto cleanup;
	}
	fclose(fp);
    }

    dbip = db_create(dbpath, 5);
    if (!dbip) {
	printf("FAIL: mesh lod db_create\n");
	ret = 1;
	goto cleanup;
    }

    {
	struct rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
	if (!wdbp) {
	    printf("FAIL: mesh lod wdb_dbopen\n");
	    ret = 1;
	    goto cleanup;
	}

	if (mk_bot(wdbp, objname, RT_BOT_SURFACE, RT_BOT_CCW, 0,
		   vertexCount, faceCount, vertices, faces, NULL, NULL) < 0) {
	    printf("FAIL: mesh lod mk_bot\n");
	    ret = 1;
	    goto cleanup;
	}

	{
	    int invalidFaces[3] = {0, 1, vertexCount};
	    if (mk_bot(wdbp, invalidBotObjname, RT_BOT_SURFACE, RT_BOT_CCW,
		       0, vertexCount, 1, vertices, invalidFaces, NULL, NULL) < 0) {
		printf("FAIL: mesh lod invalid mk_bot\n");
		ret = 1;
		goto cleanup;
	    }
	}

	{
	    point_t center = VINIT_ZERO;
	    if (mk_sph(wdbp, meshObjname, center, 1.0) < 0) {
		printf("FAIL: mesh lod mk_sph\n");
		ret = 1;
		goto cleanup;
	    }
	}
    }

    if (brlobol_mesh_lod_cache_status(dbip, objname, &cacheStatus) !=
	BRLCAD_OK ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	cacheStatus.has_cache_key || cacheStatus.has_cached_payload ||
	cacheStatus.stale_cache_entry) {
	printf("FAIL: mesh lod initial cache status\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_mesh_lod_cache_refresh(dbip, invalidBotObjname,
				       &cacheStatus) != BRLCAD_ERROR ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	cacheStatus.has_cache_key || cacheStatus.has_cached_payload ||
	brlobol_mesh_lod_get(dbip, invalidBotObjname)) {
	printf("FAIL: mesh lod invalid BoT cache rejection\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_mesh_lod_cache_refresh(dbip, objname, &cacheStatus) !=
	BRLCAD_OK ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	!cacheStatus.has_cache_key || !cacheStatus.has_cached_payload ||
	cacheStatus.stale_cache_entry ||
	!cacheStatus.generated_cache_entry || !cacheStatus.cache_key) {
	printf("FAIL: mesh lod refresh status\n");
	ret = 1;
	goto cleanup;
    }

    lod = brlobol_mesh_lod_get(dbip, objname);
    if (!lod) {
	printf("FAIL: mesh lod get\n");
	ret = 1;
	goto cleanup;
    }

    memshrinkLevel = first_available_level(lod);
    if (memshrinkLevel < 0) {
	printf("FAIL: mesh lod level data\n");
	ret = 1;
	goto cleanup;
    }

    {
	struct bv_view_info info = BV_VIEW_INFO_INIT;
	info.size = 0.01;
	if (brlobol_mesh_lod_load_view(lod, &info, 0) < 0 ||
	    check_mesh_lod_payload("mesh lod view", lod,
				   static_cast<size_t>(faceCount),
				   static_cast<size_t>(vertexCount), 0)) {
	    ret = 1;
	    goto cleanup;
	}
    }

    {
	struct BRLObolMeshLod *meshLod = NULL;
	int invalidStoreFaces[3] = {0, vertexCount, 1};
	if (brlobol_mesh_lod_cache_store_mesh(dbip, meshObjname,
					      reinterpret_cast<const point_t *>(vertices),
					      static_cast<size_t>(vertexCount),
					      reinterpret_cast<const vect_t *>(detailNormals),
					      faces, static_cast<size_t>(faceCount), 424242ULL, 0.66,
					      &cacheStatus) != BRLCAD_OK ||
	    !cacheStatus.directory_found || cacheStatus.is_bot ||
	    !cacheStatus.has_cache_key || !cacheStatus.has_cached_payload ||
	    cacheStatus.stale_cache_entry ||
	    !cacheStatus.generated_cache_entry || !cacheStatus.cache_key) {
	    printf("FAIL: mesh lod store generated mesh status\n");
	    ret = 1;
	    goto cleanup;
	}
	meshLod = brlobol_mesh_lod_get(dbip, meshObjname);
	if (!meshLod || first_available_level(meshLod) < 0 ||
	    check_mesh_lod_payload("mesh lod generated mesh", meshLod,
				   static_cast<size_t>(faceCount),
				   static_cast<size_t>(vertexCount), 1)) {
	    printf("FAIL: mesh lod generated mesh data\n");
	    ret = 1;
	}
	if (meshLod)
	    brlobol_mesh_lod_destroy(meshLod);
	if (ret)
	    goto cleanup;

	if (brlobol_mesh_lod_cache_store_mesh(dbip, meshObjname,
					      reinterpret_cast<const point_t *>(vertices),
					      static_cast<size_t>(vertexCount), NULL, invalidStoreFaces, 1,
					      777777ULL, 0.66, &cacheStatus) != BRLCAD_ERROR ||
	    !cacheStatus.has_cache_key ||
	    !cacheStatus.has_cached_payload ||
	    cacheStatus.cache_key != 424242ULL) {
	    printf("FAIL: mesh lod invalid generated mesh store preservation\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    if (test_detail_callbacks(dbip, meshObjname, vertices, faces,
			      detailNormals, static_cast<size_t>(vertexCount),
			      static_cast<size_t>(faceCount))) {
	ret = 1;
	goto cleanup;
    }

    {
	int activeLevel;
	struct BRLObolMeshLodData shrinkData;
	struct BRLObolMeshLodInfo shrinkInfo = BRLOBOL_MESH_LOD_INFO_INIT;
	if (brlobol_mesh_lod_load_level(lod, memshrinkLevel, 0) !=
	    memshrinkLevel ||
	    (activeLevel = brlobol_mesh_lod_current_level(lod)) < 0 ||
	    !brlobol_mesh_lod_has_active_data(lod)) {
	    printf("FAIL: mesh lod memshrink setup\n");
	    ret = 1;
	    goto cleanup;
	}
	brlobol_mesh_lod_memshrink(lod);
	if (brlobol_mesh_lod_current_level(lod) != activeLevel ||
	    brlobol_mesh_lod_has_active_data(lod) ||
	    brlobol_mesh_lod_data_get(lod, &shrinkData) ||
	    brlobol_mesh_lod_info_get(lod, &shrinkInfo) ||
	    shrinkInfo.active_level != activeLevel) {
	    printf("FAIL: mesh lod memshrink status\n");
	    ret = 1;
	    goto cleanup;
	}
	if (brlobol_mesh_lod_load_level(lod, activeLevel, 0) != activeLevel ||
	    check_mesh_lod_payload("mesh lod reload after memshrink", lod,
				   static_cast<size_t>(faceCount),
				   static_cast<size_t>(vertexCount), 0)) {
	    printf("FAIL: mesh lod reload after memshrink\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    brlobol_mesh_lod_destroy(lod);
    lod = NULL;
    db_close(dbip);
    dbip = db_open(dbpath, DB_OPEN_READWRITE);
    if (!dbip || db_dirbuild(dbip) < 0) {
	printf("FAIL: mesh lod reopen\n");
	ret = 1;
	goto cleanup;
    }

    lod = brlobol_mesh_lod_get(dbip, objname);
    if (!lod || brlobol_mesh_lod_load_view(lod, NULL, 0) < 0 ||
	check_mesh_lod_payload("mesh lod reopen", lod,
			       static_cast<size_t>(faceCount),
			       static_cast<size_t>(vertexCount), 0)) {
	printf("FAIL: mesh lod reopen data\n");
	ret = 1;
	goto cleanup;
    }

    brlobol_mesh_lod_destroy(lod);
    lod = NULL;

    if (brlobol_mesh_lod_cache_invalidate(dbip, objname, &cacheStatus) !=
	BRLCAD_OK ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	!cacheStatus.cleared_cache_entry ||
	!cacheStatus.cleared_cache_key ||
	cacheStatus.has_cache_key || cacheStatus.has_cached_payload ||
	cacheStatus.stale_cache_entry) {
	printf("FAIL: mesh lod invalidate status\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_mesh_lod_get(dbip, objname)) {
	printf("FAIL: mesh lod get after invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_mesh_lod_cache_update(dbip, objname) != BRLCAD_OK) {
	printf("FAIL: mesh lod compatibility update\n");
	ret = 1;
	goto cleanup;
    }

    lod = brlobol_mesh_lod_get(dbip, objname);
    if (!lod || brlobol_mesh_lod_load_view(lod, NULL, 0) < 0 ||
	check_mesh_lod_payload("mesh lod update after invalidate", lod,
			       static_cast<size_t>(faceCount),
			       static_cast<size_t>(vertexCount), 0)) {
	printf("FAIL: mesh lod compatibility update data\n");
	ret = 1;
	goto cleanup;
    }

    {
	struct BRLObolMeshLodData nullData;
	struct BRLObolMeshLodDetail nullDetail;
	struct BRLObolMeshLodInfo nullInfo;
	nullInfo.active_level = 7;
	nullDetail.face_count = 42;
	brlobol_mesh_lod_cache_status_init(&cacheStatus);
	cacheStatus.cache_key = 99;
	brlobol_mesh_lod_cache_status_init(&cacheStatus);
	brlobol_mesh_lod_detail_init(&nullDetail);
	if (brlobol_mesh_lod_load_level(NULL, 0, 0) >= 0 ||
	    brlobol_mesh_lod_load_view(NULL, NULL, 0) >= 0 ||
	    brlobol_mesh_lod_current_level(NULL) >= 0 ||
	    brlobol_mesh_lod_has_active_data(NULL) ||
	    brlobol_mesh_lod_data_get(NULL, &nullData) ||
	    brlobol_mesh_lod_info_get(NULL, &nullInfo) ||
	    nullDetail.face_count ||
	    brlobol_mesh_lod_detail_callbacks_set(NULL, NULL, NULL, NULL,
		    NULL) ||
	    nullInfo.active_level != -1 ||
	    cacheStatus.cache_key ||
	    brlobol_mesh_lod_cache_status(NULL, objname,
					  &cacheStatus) != BRLCAD_ERROR ||
	    brlobol_mesh_lod_cache_status(dbip, NULL,
					  &cacheStatus) != BRLCAD_ERROR ||
	    brlobol_mesh_lod_cache_status(dbip, objname,
					  NULL) != BRLCAD_ERROR ||
	    brlobol_mesh_lod_cache_invalidate(NULL, objname,
					      &cacheStatus) != BRLCAD_ERROR) {
	    printf("FAIL: mesh lod null handling\n");
	    ret = 1;
	}
	brlobol_mesh_lod_detail_callbacks_clear(NULL);
	brlobol_mesh_lod_memshrink(NULL);
    }

cleanup:
    if (lod)
	brlobol_mesh_lod_destroy(lod);
    if (dbip) {
	brlobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
    }
    if (dbpath[0])
	bu_file_delete(dbpath);
    bu_dirclear(cacheDir);
    if (vertices)
	bu_free(vertices, "mesh lod vertices");
    if (detailNormals)
	bu_free(detailNormals, "mesh lod detail normals");
    if (faces)
	bu_free(faces, "mesh lod faces");

    return ret ? 1 : 0;
}

/*
 * Local Variables:
 * tab-width: 8
 * mode: C++
 * indent-tabs-mode: t
 * c-file-style: "stroustrup"
 * End:
 * ex: shiftwidth=4 tabstop=8
 */
