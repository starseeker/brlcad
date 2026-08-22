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

#include "BObol/BLodRealization.h"
#include "BObol/BMeshLodCache.h"

#include <Obol/cad/SoCADAssembly.h>

#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "bu/malloc.h"
#include "wdb.h"
#include "rt/wdb.h"

#include <algorithm>
#include <atomic>
#include <cmath>
#include <cstdio>
#include <cstring>
#include <limits>
#include <thread>
#include <vector>

static int
fastf_equal(fastf_t a, fastf_t b)
{
    return std::fabs(a - b) <= SMALL_FASTF;
}

static fastf_t
mesh_pop_expected_snap(fastf_t value, fastf_t minimum, fastf_t maximum,
	uint8_t bits)
{
    if (!(maximum > minimum) || !bits)
	return value;
    const double mask = std::ldexp(1.0,
	BOBOL_MESH_LOD_QUANTIZATION_BITS -
	    std::min<int>(BOBOL_MESH_LOD_QUANTIZATION_BITS, bits));
    const double scaled = std::max(0.0, std::min(65535.0,
	(static_cast<double>(value) - minimum) /
	    (static_cast<double>(maximum) - minimum) * 65535.0));
    const double code = std::floor(scaled);
    const double cell = std::floor(code / mask);
    const double snapped = std::min(65535.0, (cell + 0.5) * mask);
    return static_cast<fastf_t>(
	(snapped / 65535.0) * (maximum - minimum) + minimum);
}

static int
mesh_pop_points_equal(const point_t a, const point_t b)
{
    return fastf_equal(a[X], b[X]) && fastf_equal(a[Y], b[Y]) &&
	fastf_equal(a[Z], b[Z]);
}

struct mesh_lod_suffix_test_data {
    const struct BObolMeshLodHierarchyInfo *hierarchy = NULL;
    int nextCut = -1;
    size_t cumulativePoints = 0;
    size_t streamedPoints = 0;
    size_t streamedFaces = 0;
};

struct mesh_lod_chunk_test_data {
    const struct BObolMeshLodHierarchyInfo *hierarchy = NULL;
    std::vector<uint32_t> seen;
    uint64_t faces = 0;
    uint64_t points = 0;
};

struct mesh_lod_preview_test_data {
    size_t calls = 0;
    size_t faces = 0;
    size_t points = 0;
    int minCut = -1;
    int maxCut = -1;
    bool valid = false;
};

static void
mesh_lod_preview_test_callback(unsigned long long cacheKey,
	const struct BObolMeshLodData *data,
	const struct BObolMeshLodHierarchyInfo *hierarchy,
	void *callbackData)
{
    mesh_lod_preview_test_data *test =
	static_cast<mesh_lod_preview_test_data *>(callbackData);
    if (!test)
	return;
    ++test->calls;
    test->faces = data ? data->face_count : 0;
    test->points = data ? data->point_count : 0;
    test->minCut = hierarchy ? hierarchy->min_cut : -1;
    test->maxCut = hierarchy ? hierarchy->max_cut : -1;
    test->valid = cacheKey && data && hierarchy && data->faces &&
	data->points_orig && data->face_count && data->point_orig_count &&
	hierarchy->min_cut >= 0 &&
	hierarchy->max_cut >= hierarchy->min_cut &&
	data->face_count == hierarchy->cuts[hierarchy->min_cut].face_count &&
	data->point_orig_count ==
	    hierarchy->cuts[hierarchy->min_cut].point_count;
}

static int
mesh_lod_chunk_test_callback(uint32_t chunkId, int cut,
	const point_t *points, size_t pointCount,
	const uint32_t *faces, size_t faceCount,
	const vect_t *normals, size_t normalCount, void *callbackData)
{
    mesh_lod_chunk_test_data *test =
	static_cast<mesh_lod_chunk_test_data *>(callbackData);
    if (!test || !test->hierarchy ||
	chunkId >= test->hierarchy->chunk_count ||
	cut < test->hierarchy->min_cut || cut > test->hierarchy->max_cut)
	return 0;
    const BObolMeshLodChunkInfo &info = test->hierarchy->chunks[chunkId];
    const BObolMeshLodChunkCutInfo &selected = info.cuts[cut];
    if (pointCount != selected.point_count ||
	faceCount != selected.face_count ||
	(pointCount && !points) || (faceCount && !faces) ||
	normalCount != (test->hierarchy->has_normals ? faceCount * 3 : 0) ||
	(normalCount && !normals))
	return 0;
    for (size_t index = 0; index < faceCount * 3; ++index)
	if (faces[index] >= pointCount)
	    return 0;
    for (size_t point = 0; point < pointCount; ++point)
	for (int axis = 0; axis < 3; ++axis)
	    if (!std::isfinite(points[point][axis]) ||
		points[point][axis] < info.bmin[axis] ||
		points[point][axis] > info.bmax[axis])
		return 0;
    test->seen.push_back(chunkId);
    test->faces += faceCount;
    test->points += pointCount;
    return 1;
}

static int
mesh_lod_suffix_test_callback(int cut, const point_t *points,
	size_t pointCount, const uint32_t *faces, size_t faceCount,
	const vect_t *normals, size_t normalCount, void *callbackData)
{
    mesh_lod_suffix_test_data *test =
	static_cast<mesh_lod_suffix_test_data *>(callbackData);
    if (!test || !test->hierarchy || cut != test->nextCut ||
	cut <= 0 || cut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return 0;
    const BObolMeshLodHierarchyInfo &hierarchy = *test->hierarchy;
    if (hierarchy.cuts[cut].point_count < hierarchy.cuts[cut - 1].point_count ||
	hierarchy.cuts[cut].face_count < hierarchy.cuts[cut - 1].face_count ||
	pointCount != hierarchy.cuts[cut].point_count -
	    hierarchy.cuts[cut - 1].point_count ||
	faceCount != hierarchy.cuts[cut].face_count -
	    hierarchy.cuts[cut - 1].face_count ||
	(pointCount && !points) || (faceCount && !faces) ||
	(normalCount && !normals) ||
	(pointCount && reinterpret_cast<uintptr_t>(points) %
	    alignof(fastf_t) != 0) ||
	(faceCount && reinterpret_cast<uintptr_t>(faces) %
	    alignof(uint32_t) != 0) ||
	(normalCount && reinterpret_cast<uintptr_t>(normals) %
	    alignof(fastf_t) != 0) ||
	normalCount != (hierarchy.has_normals ? faceCount * 3 : 0))
	return 0;
    test->cumulativePoints += pointCount;
    for (size_t index = 0; index < faceCount * 3; ++index) {
	if (static_cast<size_t>(faces[index]) >= test->cumulativePoints)
	    return 0;
    }
    test->streamedPoints += pointCount;
    test->streamedFaces += faceCount;
    test->nextCut++;
    return 1;
}

static long double
mesh_signed_six_volume(const struct BObolMeshLodData &data)
{
    if (!data.points_orig || !data.faces || data.point_orig_count == 0 ||
	data.face_count == 0)
	return 0.0L;
    const point_t &origin =
	data.points_orig[static_cast<size_t>(data.faces[0])];
    long double volume = 0.0L;
    for (size_t faceIndex = 0; faceIndex < data.face_count; ++faceIndex) {
	const uint32_t *face = &data.faces[faceIndex * 3];
	long double a[3];
	long double b[3];
	long double c[3];
	for (int axis = 0; axis < 3; ++axis) {
	    a[axis] =
		data.points_orig[static_cast<size_t>(face[0])][axis] -
		origin[axis];
	    b[axis] =
		data.points_orig[static_cast<size_t>(face[1])][axis] -
		origin[axis];
	    c[axis] =
		data.points_orig[static_cast<size_t>(face[2])][axis] -
		origin[axis];
	}
	volume +=
	    a[0] * (b[1] * c[2] - b[2] * c[1]) -
	    a[1] * (b[0] * c[2] - b[2] * c[0]) +
	    a[2] * (b[0] * c[1] - b[1] * c[0]);
    }
    return volume;
}

static int
check_mesh_lod_payload(const char *label,
		       struct BObolMeshLod *lod,
		       size_t maxFaceCount,
		       size_t maxPointCount,
		       int requireNormals)
{
    struct BObolMeshLodData data;
    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;

    if (!bobol_mesh_lod_data_get(lod, &data) ||
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

    if (!bobol_mesh_lod_has_active_data(lod)) {
	printf("FAIL: %s mesh lod active data status\n", label);
	return 1;
    }

    if (!bobol_mesh_lod_info_get(lod, &info) ||
	!bobol_mesh_lod_hierarchy_info_get(lod, &hierarchy) ||
	info.active_cut < 0 ||
	info.face_count != data.face_count ||
	info.point_count != data.point_count ||
	info.point_orig_count != data.point_orig_count ||
	info.normal_count != data.normal_count ||
	!info.has_faces || !info.has_points ||
	!info.has_original_points ||
	info.has_normals != (data.normals ? 1 : 0) ||
	hierarchy.has_normals != info.has_normals ||
	hierarchy.shaded_cull_backfaces !=
	    info.shaded_cull_backfaces) {
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

static int
first_available_cut(struct BObolMeshLod *lod)
{
    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (!bobol_mesh_lod_hierarchy_info_get(lod, &hierarchy))
	return -1;
    for (int cut = hierarchy.min_cut; cut <= hierarchy.max_cut; ++cut) {
	if (bobol_mesh_lod_load_cut(lod, cut, 0) < 0)
	    return -1;
	if (bobol_mesh_lod_has_active_data(lod))
	    return cut;
    }
    return -1;
}

static int
load_terminal_cut(struct BObolMeshLod *lod)
{
    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (!bobol_mesh_lod_hierarchy_info_get(lod, &hierarchy))
	return -1;
    return bobol_mesh_lod_load_cut(lod, hierarchy.max_cut, 0);
}

int
main(int argc, char *argv[])
{
    const char *objname = "bobol_lod_bot";
    const char *duplicateObjname = "bobol_lod_duplicate_bot";
    const char *solidCwObjname = "bobol_lod_solid_cw_bot";
    const char *solidUnorientedObjname =
	"bobol_lod_solid_unoriented_bot";
    const char *brokenUnorientedObjname =
	"bobol_lod_broken_unoriented_bot";
    const char *meshObjname = "bobol_lod_mesh_payload";
    const char *translatedMeshObjname = "bobol_lod_translated_mesh_payload";
    const char *invalidBotObjname = "bobol_invalid_lod_bot";
    const char *degenerateBotObjname = "bobol_degenerate_lod_bot";
    /* Slightly exceed one private chunk so every cache run exercises
     * independent page serialization and selective reads. */
    const int grid = 190;
    const int vertexCount = (grid + 1) * (grid + 1);
    const int faceCount = grid * grid * 2;
    char dbpath[MAXPATHLEN] = {0};
    char cacheDir[MAXPATHLEN] = {0};
    fastf_t *vertices = NULL;
    fastf_t *translatedVertices = NULL;
    fastf_t *detailNormals = NULL;
    int *faces = NULL;
    struct db_i *dbip = NULL;
    struct BObolMeshLod *lod = NULL;
    int memshrinkCut = -1;
    unsigned long long sharedCacheKey = 0;
    int ret = 0;
    struct BObolMeshLodCacheStatus cacheStatus =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    struct BObolMeshLodHierarchyInfo lodHierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;

    bu_setprogname(argv[0]);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    bu_dir(cacheDir, MAXPATHLEN, BU_DIR_CURR,
	   "bobol_mesh_lod_cache_test", NULL);
    bu_dirclear(cacheDir);
    bu_mkdir(cacheDir);
    bu_setenv("BU_DIR_CACHE", cacheDir, 1);

    vertices = static_cast<fastf_t *>(bu_calloc(
					  static_cast<size_t>(vertexCount) * 3, sizeof(fastf_t),
					  "mesh lod vertices"));
    translatedVertices = static_cast<fastf_t *>(bu_calloc(
	static_cast<size_t>(vertexCount) * 3, sizeof(fastf_t),
	"translated mesh lod vertices"));
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
	    translatedVertices[idx * 3 + X] =
		vertices[idx * 3 + X] - 29640.022;
	    translatedVertices[idx * 3 + Y] =
		vertices[idx * 3 + Y] + 5283.2;
	    translatedVertices[idx * 3 + Z] =
		vertices[idx * 3 + Z] + 1300.0;
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
	if (mk_bot(wdbp, duplicateObjname, RT_BOT_SURFACE, RT_BOT_CCW, 0,
		   vertexCount, faceCount, vertices, faces, NULL, NULL) < 0) {
	    printf("FAIL: duplicate mesh lod mk_bot\n");
	    ret = 1;
	    goto cleanup;
	}

	{
	    fastf_t solidVertices[12] = {
		0.0, 0.0, 0.0,
		1.0, 0.0, 0.0,
		0.0, 1.0, 0.0,
		0.0, 0.0, 1.0
	    };
	    /* Clockwise as viewed from the exterior. */
	    int solidFaces[12] = {
		0, 1, 2,
		0, 3, 1,
		0, 2, 3,
		1, 3, 2
	    };
	    if (mk_bot(wdbp, solidCwObjname, RT_BOT_SOLID, RT_BOT_CW, 0,
		       4, 4, solidVertices, solidFaces, NULL, NULL) < 0) {
		printf("FAIL: mesh lod solid CW mk_bot\n");
		ret = 1;
		goto cleanup;
	    }
	    if (mk_bot(wdbp, solidUnorientedObjname, RT_BOT_SOLID,
		       RT_BOT_UNORIENTED, 0, 4, 4, solidVertices, solidFaces,
		       NULL, NULL) < 0) {
		printf("FAIL: mesh lod solid unoriented mk_bot\n");
		ret = 1;
		goto cleanup;
	    }
	    int brokenFaces[12];
	    std::memcpy(brokenFaces, solidFaces, sizeof(brokenFaces));
	    std::swap(brokenFaces[1], brokenFaces[2]);
	    if (mk_bot(wdbp, brokenUnorientedObjname, RT_BOT_SOLID,
		       RT_BOT_UNORIENTED, 0, 4, 4, solidVertices,
		       brokenFaces, NULL, NULL) < 0) {
		printf("FAIL: mesh lod broken unoriented mk_bot\n");
		ret = 1;
		goto cleanup;
	    }
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
	    /* Repeated-index triangles occur in otherwise drawable legacy BoTs.
	     * They have no visible area and must not poison the complete PoP
	     * hierarchy for the valid neighboring faces. */
	    int degenerateFaces[9] = {
		0, 1, grid + 1,
		0, 0, 1,
		0, grid + 1, grid + 2
	    };
	    if (mk_bot(wdbp, degenerateBotObjname, RT_BOT_SURFACE,
		    RT_BOT_CCW, 0, vertexCount, 3, vertices,
		    degenerateFaces, NULL, NULL) < 0) {
		printf("FAIL: mesh lod degenerate mk_bot\n");
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
	    if (mk_sph(wdbp, translatedMeshObjname, center, 1.0) < 0) {
		printf("FAIL: translated mesh lod mk_sph\n");
		ret = 1;
		goto cleanup;
	    }
	}
    }

    if (bobol_mesh_lod_cache_status(dbip, objname, &cacheStatus) !=
	BRLCAD_OK ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	cacheStatus.has_cache_key || cacheStatus.has_cached_payload ||
	cacheStatus.stale_cache_entry) {
	printf("FAIL: mesh lod initial cache status\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_mesh_lod_cache_refresh(dbip, invalidBotObjname,
				       &cacheStatus) != BRLCAD_ERROR ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	cacheStatus.has_cache_key || cacheStatus.has_cached_payload ||
	bobol_mesh_lod_get(dbip, invalidBotObjname)) {
	printf("FAIL: mesh lod invalid BoT cache rejection\n");
	ret = 1;
	goto cleanup;
    }

    {
	struct BObolMeshLod *degenerate = NULL;
	struct BObolMeshLodHierarchyInfo hierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	struct BObolMeshLodData data = {};
	const int refreshStatus = bobol_mesh_lod_cache_refresh(
	    dbip, degenerateBotObjname, &cacheStatus);
	if (refreshStatus != BRLCAD_OK ||
	    !(degenerate = bobol_mesh_lod_get(dbip, degenerateBotObjname)) ||
	    !bobol_mesh_lod_hierarchy_info_get(degenerate, &hierarchy) ||
	    bobol_mesh_lod_load_resident_cut(
		degenerate, hierarchy.max_cut, 0) != hierarchy.max_cut ||
	    !bobol_mesh_lod_data_get(degenerate, &data) ||
	    data.face_count != 2) {
	    printf("FAIL: mesh lod repeated-index face sanitization "
		   "refresh=%d lod=%p min=%d max=%d current=%d "
		   "faces=%zu points=%zu\n", refreshStatus,
		   static_cast<void *>(degenerate), hierarchy.min_cut,
		   hierarchy.max_cut,
		   degenerate ? bobol_mesh_lod_current_cut(degenerate) : -1,
		   data.face_count, data.point_count);
	    if (degenerate)
		bobol_mesh_lod_destroy(degenerate);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(degenerate);
    }

    {
	/* Cold database coverage has already paid to import a potentially huge
	 * BoT in order to publish its leaf bounds.  Cache generation must be able
	 * to consume that caller-owned internal directly without taking
	 * ownership or rereading the database object. */
	struct directory *dp = db_lookup(dbip, objname, LOOKUP_QUIET);
	struct rt_db_internal intern;
	struct BObolMeshLod *opened = NULL;
	struct BObolMeshLodHierarchyInfo openedHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	struct BObolMeshLodData openedData;
	mesh_lod_preview_test_data previewData;
	RT_DB_INTERNAL_INIT(&intern);
	if (!dp || rt_db_get_internal(&intern, dp, dbip, NULL) < 0 ||
	    intern.idb_type != ID_BOT || !intern.idb_ptr ||
	    !(opened = bobol_mesh_lod_cache_refresh_open(
		dbip, objname,
		static_cast<const struct rt_bot_internal *>(intern.idb_ptr),
		&cacheStatus, mesh_lod_preview_test_callback, &previewData)) ||
	    previewData.calls != 1 || !previewData.valid ||
	    !cacheStatus.directory_found || !cacheStatus.is_bot ||
	    !cacheStatus.has_cache_key || !cacheStatus.has_cached_payload ||
	    cacheStatus.stale_cache_entry ||
	    !cacheStatus.generated_cache_entry || !cacheStatus.cache_key ||
	    !bobol_mesh_lod_hierarchy_info_get(
		opened, &openedHierarchy) ||
	    bobol_mesh_lod_current_cut(opened) !=
		openedHierarchy.min_cut ||
	    !bobol_mesh_lod_data_get(opened, &openedData) ||
	    !openedData.face_count || !openedData.point_count) {
	    printf("FAIL: mesh lod staged BoT open refresh status\n");
	    if (opened)
		bobol_mesh_lod_destroy(opened);
	    if (intern.idb_ptr)
		rt_db_free_internal(&intern);
	    ret = 1;
	    goto cleanup;
	}
	sharedCacheKey = cacheStatus.cache_key;
	rt_db_free_internal(&intern);
	if (bobol_mesh_lod_load_resident_cut(
		opened, openedHierarchy.max_cut, 0) != openedHierarchy.max_cut ||
	    !bobol_mesh_lod_data_get(opened, &openedData) ||
	    openedData.face_count != static_cast<size_t>(faceCount) ||
	    openedData.point_count != static_cast<size_t>(vertexCount)) {
	    printf("FAIL: generated open prefix did not append cached suffix\n");
	    bobol_mesh_lod_destroy(opened);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(opened);
    }

    {
	/* A different source name with identical content must reuse the complete
	 * immutable hierarchy.  The old marker-only shortcut returned an
	 * uninitialized state here and made cold generate-and-open fail. */
	struct directory *dp = db_lookup(dbip, duplicateObjname, LOOKUP_QUIET);
	struct rt_db_internal intern;
	struct BObolMeshLod *opened = NULL;
	struct BObolMeshLodHierarchyInfo hierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	struct BObolMeshLodData data;
	RT_DB_INTERNAL_INIT(&intern);
	if (!dp || rt_db_get_internal(&intern, dp, dbip, NULL) < 0 ||
	    intern.idb_type != ID_BOT || !intern.idb_ptr ||
	    !(opened = bobol_mesh_lod_cache_refresh_open(
		dbip, duplicateObjname,
		static_cast<const struct rt_bot_internal *>(intern.idb_ptr),
		&cacheStatus, NULL, NULL)) ||
	    cacheStatus.cache_key != sharedCacheKey ||
	    !bobol_mesh_lod_hierarchy_info_get(opened, &hierarchy) ||
	    bobol_mesh_lod_current_cut(opened) != hierarchy.min_cut ||
	    !bobol_mesh_lod_data_get(opened, &data) ||
	    !data.face_count || !data.point_count) {
	    printf("FAIL: content-identical mesh hierarchy reuse\n");
	    if (opened)
		bobol_mesh_lod_destroy(opened);
	    if (intern.idb_ptr)
		rt_db_free_internal(&intern);
	    ret = 1;
	    goto cleanup;
	}
	rt_db_free_internal(&intern);
	bobol_mesh_lod_destroy(opened);
    }

    if (bobol_mesh_lod_cache_refresh(dbip, objname, &cacheStatus) !=
	BRLCAD_OK ||
	!cacheStatus.directory_found || !cacheStatus.is_bot ||
	!cacheStatus.has_cache_key || !cacheStatus.has_cached_payload ||
	cacheStatus.stale_cache_entry ||
	!cacheStatus.generated_cache_entry || !cacheStatus.cache_key) {
	printf("FAIL: mesh lod refresh status\n");
	ret = 1;
	goto cleanup;
    }

    {
	struct BObolMeshLod *cachedPrefix =
	    bobol_mesh_lod_get_cached_prefix(dbip, cacheStatus.cache_key);
	struct BObolMeshLodData cachedData;
	struct BObolMeshLodHierarchyInfo cachedHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	if (!cachedPrefix ||
	    !bobol_mesh_lod_hierarchy_info_get(
		cachedPrefix, &cachedHierarchy) ||
	    bobol_mesh_lod_load_resident_cut(
		cachedPrefix, cachedHierarchy.max_cut, 0) !=
		cachedHierarchy.max_cut ||
	    !bobol_mesh_lod_data_get(cachedPrefix, &cachedData) ||
	    cachedData.face_count != static_cast<size_t>(faceCount) ||
	    cachedData.point_count != static_cast<size_t>(vertexCount)) {
	    printf("FAIL: mesh lod direct cached-prefix open\n");
	    if (cachedPrefix)
		bobol_mesh_lod_destroy(cachedPrefix);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(cachedPrefix);
    }

    lod = bobol_mesh_lod_get(dbip, objname);
    if (!lod) {
	printf("FAIL: mesh lod get\n");
	ret = 1;
	goto cleanup;
    }

    if (!bobol_mesh_lod_hierarchy_info_get(lod, &lodHierarchy)) {
	printf("FAIL: mesh lod hierarchy metadata\n");
	ret = 1;
	goto cleanup;
    }
    {
	bool validSchedule = lodHierarchy.cut_count ==
		BOBOL_MESH_LOD_CUT_COUNT_MAX &&
	    lodHierarchy.min_cut >= 0 &&
	    lodHierarchy.max_cut == BOBOL_MESH_LOD_CUT_COUNT_MAX - 1;
	uint64_t priorFaces = 0;
	uint64_t priorPoints = 0;
	uint64_t priorBytes = 0;
	double priorError = std::numeric_limits<double>::infinity();
	for (uint32_t cut = 0; validSchedule &&
		cut < lodHierarchy.cut_count; ++cut) {
	    const BObolMeshLodCutInfo &info = lodHierarchy.cuts[cut];
	    int refinedAxes = 0;
	    for (int axis = 0; axis < 3; ++axis) {
		const uint8_t bits = info.quantization_bits[axis];
		if (bits < 1 || bits > BOBOL_MESH_LOD_QUANTIZATION_BITS)
		    validSchedule = false;
		if (cut > 0) {
		    const uint8_t prior = lodHierarchy.cuts[cut - 1].
			quantization_bits[axis];
		    if (bits < prior || bits > prior + 1)
			validSchedule = false;
		    if (bits == prior + 1)
			++refinedAxes;
		}
	    }
	    if ((cut > 0 && refinedAxes != 1) ||
		info.face_count < priorFaces ||
		info.point_count < priorPoints ||
		info.resident_bytes < priorBytes ||
		!std::isfinite(info.object_error) ||
		info.object_error < 0.0 ||
		info.object_error > priorError ||
		(info.exact != 0) !=
		    (cut + 1 == lodHierarchy.cut_count))
		validSchedule = false;
	    priorFaces = info.face_count;
	    priorPoints = info.point_count;
	    priorBytes = info.resident_bytes;
	    priorError = info.object_error;
	}
	const BObolMeshLodCutInfo &terminal =
	    lodHierarchy.cuts[lodHierarchy.max_cut];
	for (int axis = 0; axis < 3; ++axis)
	    if (terminal.quantization_bits[axis] !=
		    BOBOL_MESH_LOD_QUANTIZATION_BITS)
		validSchedule = false;
	if (!validSchedule || terminal.face_count !=
		static_cast<uint64_t>(faceCount) ||
	    terminal.point_count != static_cast<uint64_t>(vertexCount) ||
	    terminal.object_error > 0.0) {
	    printf("FAIL: mesh lod admissible cut schedule contract\n");
	    ret = 1;
	    goto cleanup;
	}

	/* Face activation and vertex presentation are one contract: a face may
	 * remain absent only while at least two of its vertices share a floored
	 * 16-bit prefix cell.  In particular, exact cell boundaries (including
	 * the domain minimum used throughout this grid) must not be snapped by a
	 * different floor/ceil rule. */
	for (int cut = lodHierarchy.min_cut;
		cut <= lodHierarchy.max_cut; ++cut) {
	    struct BObolMeshLodData cutData;
	    if (bobol_mesh_lod_load_cut(lod, cut, 0) != cut ||
		!bobol_mesh_lod_data_get(lod, &cutData) ||
		cutData.face_count != lodHierarchy.cuts[cut].face_count ||
		cutData.point_count != lodHierarchy.cuts[cut].point_count) {
		printf("FAIL: PoP quantization contract cut %d unavailable\n", cut);
		ret = 1;
		goto cleanup;
	    }
	    for (size_t pointIndex = 0;
		    pointIndex < cutData.point_count; ++pointIndex) {
		for (int axis = 0; axis < 3; ++axis) {
		    const fastf_t expected = cut == lodHierarchy.max_cut ?
			cutData.points_orig[pointIndex][axis] :
			mesh_pop_expected_snap(
			    cutData.points_orig[pointIndex][axis],
			    lodHierarchy.quantization_min[axis],
			    lodHierarchy.quantization_max[axis],
			    lodHierarchy.cuts[cut].quantization_bits[axis]);
		    if (!fastf_equal(cutData.points[pointIndex][axis], expected)) {
			printf("FAIL: PoP displayed coordinate disagrees with "
			       "activation prefix (cut=%d point=%zu axis=%d)\n",
			       cut, pointIndex, axis);
			ret = 1;
			goto cleanup;
		    }
		}
	    }
	    if (cut == lodHierarchy.max_cut)
		continue;
	    for (size_t faceIndex = 0; faceIndex < cutData.face_count;
		    ++faceIndex) {
		const point_t &a = cutData.points[cutData.faces[3 * faceIndex]];
		const point_t &b = cutData.points[cutData.faces[3 * faceIndex + 1]];
		const point_t &c = cutData.points[cutData.faces[3 * faceIndex + 2]];
		if (mesh_pop_points_equal(a, b) ||
		    mesh_pop_points_equal(b, c) ||
		    mesh_pop_points_equal(a, c)) {
		    printf("FAIL: PoP activated a face before its prefix vertices "
			   "were distinct (cut=%d face=%zu)\n", cut, faceIndex);
		    ret = 1;
		    goto cleanup;
		}
	    }
	}

	bool validClusters = lodHierarchy.cluster_grid_resolution ==
		BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION &&
	    lodHierarchy.cluster_count > 0 &&
	    lodHierarchy.cluster_count <= BOBOL_MESH_LOD_CLUSTER_COUNT &&
	    lodHierarchy.clusters != NULL;
	uint64_t clusteredIndices = 0;
	std::vector<uint64_t> cutIndices(lodHierarchy.cut_count, 0);
	for (uint32_t cluster = 0; validClusters &&
		cluster < lodHierarchy.cluster_count; ++cluster) {
	    const BObolMeshLodClusterInfo &cell =
		lodHierarchy.clusters[cluster];
	    if (cell.cluster_id >= BOBOL_MESH_LOD_CLUSTER_COUNT ||
		(cluster && cell.cluster_id <=
		 lodHierarchy.clusters[cluster - 1].cluster_id) ||
		!cell.range_count || !cell.ranges)
		validClusters = false;
	    for (uint32_t rangeIndex = 0; validClusters &&
		    rangeIndex < cell.range_count; ++rangeIndex) {
		const BObolMeshLodClusterRange &range =
		    cell.ranges[rangeIndex];
		if (range.activation_cut >= lodHierarchy.cut_count ||
		    range.first_index % 3 || range.index_count < 3 ||
		    range.index_count % 3)
		    validClusters = false;
		else
		    cutIndices[range.activation_cut] += range.index_count;
		clusteredIndices += range.index_count;
	    }
	}
	uint64_t priorCutIndices = 0;
	for (uint32_t cut = 0; validClusters &&
		cut < lodHierarchy.cut_count; ++cut) {
	    const uint64_t cumulative =
		lodHierarchy.cuts[cut].face_count * 3;
	    if (cutIndices[cut] != cumulative - priorCutIndices)
		validClusters = false;
	    priorCutIndices = cumulative;
	}
	if (!validClusters || clusteredIndices !=
		static_cast<uint64_t>(faceCount) * 3) {
	    printf("FAIL: mesh lod clustered range contract\n");
	    ret = 1;
	    goto cleanup;
	}

	bool validChunks = lodHierarchy.chunk_count > 0 &&
	    lodHierarchy.chunks != NULL;
	uint64_t chunkFaces = 0;
	for (uint32_t chunk = 0; validChunks &&
		chunk < lodHierarchy.chunk_count; ++chunk) {
	    const BObolMeshLodChunkInfo &info = lodHierarchy.chunks[chunk];
	    if (info.chunk_id != chunk || info.min_cut < lodHierarchy.min_cut ||
		info.max_cut != lodHierarchy.max_cut)
		validChunks = false;
	    uint32_t priorChunkFaces = 0;
	    uint32_t priorChunkPoints = 0;
	    uint64_t priorChunkBytes = 0;
	    for (uint32_t cut = 0; validChunks &&
		    cut < lodHierarchy.cut_count; ++cut) {
		const BObolMeshLodChunkCutInfo &cutInfo = info.cuts[cut];
		if (cutInfo.face_count < priorChunkFaces ||
		    cutInfo.point_count < priorChunkPoints ||
		    cutInfo.resident_bytes < priorChunkBytes)
		    validChunks = false;
		priorChunkFaces = cutInfo.face_count;
		priorChunkPoints = cutInfo.point_count;
		priorChunkBytes = cutInfo.resident_bytes;
	    }
	    chunkFaces += info.cuts[info.max_cut].face_count;
	}
	std::vector<uint32_t> chunkIds(lodHierarchy.chunk_count);
	for (uint32_t chunk = 0; chunk < lodHierarchy.chunk_count; ++chunk)
	    chunkIds[chunk] = chunk;
	mesh_lod_chunk_test_data chunkData;
	chunkData.hierarchy = &lodHierarchy;
	if (!validChunks || chunkFaces != static_cast<uint64_t>(faceCount) ||
	    !bobol_mesh_lod_read_chunk_prefixes(lod, chunkIds.data(),
		chunkIds.size(), lodHierarchy.max_cut,
		mesh_lod_chunk_test_callback, &chunkData) ||
	    chunkData.seen != chunkIds ||
	    chunkData.faces != static_cast<uint64_t>(faceCount) ||
	    bobol_mesh_lod_read_chunk_prefixes(lod, chunkIds.data(),
		chunkIds.size(), lodHierarchy.max_cut + 1,
		mesh_lod_chunk_test_callback, &chunkData)) {
	    printf("FAIL: mesh lod independent chunk contract\n");
	    ret = 1;
	    goto cleanup;
	}

	/* Private pages are one logical mesh, but they must remain independently
	 * resident all the way through renderer publication.  A richer page
	 * appends a replacement range to the certified packed stream; stable
	 * compaction later drops both unrelated pages and the retired range. */
	if (lodHierarchy.chunk_count < 2) {
	    printf("FAIL: mesh lod chunk realization fixture has one page\n");
	    ret = 1;
	    goto cleanup;
	}
	const uint32_t firstChunk = 0;
	const uint32_t secondChunk = 1;
	int secondCoarseCut = lodHierarchy.chunks[secondChunk].min_cut;
	while (secondCoarseCut < lodHierarchy.max_cut &&
	       lodHierarchy.chunks[secondChunk].cuts[secondCoarseCut].face_count ==
		   lodHierarchy.chunks[secondChunk].cuts[lodHierarchy.max_cut].face_count &&
	       lodHierarchy.chunks[secondChunk].cuts[secondCoarseCut].point_count ==
		   lodHierarchy.chunks[secondChunk].cuts[lodHierarchy.max_cut].point_count)
	    ++secondCoarseCut;
	if (secondCoarseCut >= lodHierarchy.max_cut ||
	    !lodHierarchy.chunks[secondChunk].cuts[secondCoarseCut].face_count) {
	    printf("FAIL: mesh lod chunk fixture has no partial second page\n");
	    ret = 1;
	    goto cleanup;
	}
	BObolLodProgressiveMesh chunkedMesh;
	if (!chunkedMesh.updateChunksFromCache(lod, lodHierarchy,
		{firstChunk}, lodHierarchy.max_cut,
		lodHierarchy.shaded_cull_backfaces ? TRUE : FALSE)) {
	    printf("FAIL: mesh lod first private page realization\n");
	    ret = 1;
	    goto cleanup;
	}
	uint64_t firstRevision = 0;
	const std::shared_ptr<const Obol::PartGeometry> firstGeometry =
	    chunkedMesh.prepareCadGeometry(
		BOBOL_LOD_DRAW_SHADED, &firstRevision);
	const Obol::TriMesh *firstMesh = firstGeometry && firstGeometry->shaded ?
	    &*firstGeometry->shaded : NULL;
	if (!firstMesh || !firstMesh->hasAdaptiveProgressiveClusters() ||
	    firstMesh->progressiveLineage == 0 ||
	    firstMesh->progressiveClusters.size() != 1 ||
	    !chunkedMesh.canDrawChunksAtCut(
		{firstChunk}, lodHierarchy.max_cut)) {
	    printf("FAIL: mesh lod first page was flattened as a global prefix\n");
	    ret = 1;
	    goto cleanup;
	}
	if (!chunkedMesh.updateChunksFromCache(lod, lodHierarchy,
		{secondChunk}, secondCoarseCut,
		lodHierarchy.shaded_cull_backfaces ? TRUE : FALSE)) {
	    printf("FAIL: mesh lod coarse second private page realization\n");
	    ret = 1;
	    goto cleanup;
	}
	uint64_t coarseRevision = 0;
	const std::shared_ptr<const Obol::PartGeometry> coarseGeometry =
	    chunkedMesh.prepareCadGeometry(
		BOBOL_LOD_DRAW_HIDDEN_LINE, &coarseRevision);
	const Obol::TriMesh *coarseMesh =
	    coarseGeometry && coarseGeometry->shaded ?
		&*coarseGeometry->shaded : NULL;
    const Obol::WireRep *coarseWire =
	    coarseGeometry && coarseGeometry->wire ?
		&*coarseGeometry->wire : NULL;
    std::array<BObolLodCounts, BOBOL_MESH_LOD_CUT_COUNT_MAX>
	coarseSelectedCounts;
    int coarseMinimumCut = -1;
    int coarseMaximumDrawableCut = -1;
    const bool coarseCensus = chunkedMesh.drawableCountsAtCuts(
	{firstChunk, secondChunk}, FALSE, coarseSelectedCounts.data(),
	coarseSelectedCounts.size(), &coarseMinimumCut,
	&coarseMaximumDrawableCut);
    int expectedCoarseDrawableCut = lodHierarchy.min_cut - 1;
    for (int cut = lodHierarchy.min_cut; cut <= lodHierarchy.max_cut;
	    ++cut) {
	if (!chunkedMesh.canDrawChunksAtCut(
		{firstChunk, secondChunk}, cut))
	    break;
	expectedCoarseDrawableCut = cut;
    }
    const BObolMeshLodChunkCutInfo &firstCoarseCounts =
	lodHierarchy.chunks[firstChunk].cuts[secondCoarseCut];
    const BObolMeshLodChunkCutInfo &secondCoarseCounts =
	lodHierarchy.chunks[secondChunk].cuts[secondCoarseCut];
    if (!coarseMesh || !coarseWire || coarseRevision <= firstRevision ||
	    !coarseWire->derivesTriangleEdges() ||
	    !coarseWire->triangleEdges ||
	    !coarseWire->segmentPoints.empty() ||
	    coarseWire->segmentCount() != coarseMesh->indices.size() ||
	    coarseMesh->progressiveLineage != firstMesh->progressiveLineage ||
	    coarseWire->progressiveLineage != coarseMesh->progressiveLineage ||
	    coarseMesh->positions.size() <= firstMesh->positions.size() ||
	    coarseMesh->indices.size() <= firstMesh->indices.size() ||
	    !std::equal(firstMesh->positions.begin(), firstMesh->positions.end(),
		coarseMesh->positions.begin()) ||
	    !std::equal(firstMesh->indices.begin(), firstMesh->indices.end(),
		coarseMesh->indices.begin()) ||
	    coarseMesh->progressiveClusters.size() != 2 ||
	    coarseWire->progressiveClusters.size() != 2 ||
	    coarseMesh->progressiveClusters[0].residentCut !=
		lodHierarchy.max_cut ||
	    coarseWire->progressiveClusters[0].residentCut !=
		coarseMesh->progressiveClusters[0].residentCut ||
	    coarseMesh->progressiveClusters[1].residentCut <
		secondCoarseCut ||
	    coarseMesh->progressiveClusters[1].residentCut >=
		lodHierarchy.max_cut ||
	    coarseWire->progressiveClusters[1].residentCut !=
		coarseMesh->progressiveClusters[1].residentCut ||
	    !coarseCensus ||
	    coarseMinimumCut != lodHierarchy.min_cut ||
	    coarseMaximumDrawableCut != expectedCoarseDrawableCut ||
	    coarseSelectedCounts[secondCoarseCut].faceCount !=
		static_cast<uint64_t>(firstCoarseCounts.face_count) +
		secondCoarseCounts.face_count ||
	    coarseSelectedCounts[secondCoarseCut].pointCount !=
		static_cast<uint64_t>(firstCoarseCounts.point_count) +
		secondCoarseCounts.point_count ||
	    chunkedMesh.canDrawChunksAtCut(
		{secondChunk}, lodHierarchy.max_cut)) {
	    printf("FAIL: mesh lod second page did not append independently "
		   "or bulk census changed its frontier (cuts=%d/%d expected=%d)\n",
		   coarseMinimumCut, coarseMaximumDrawableCut,
		   expectedCoarseDrawableCut);
	    ret = 1;
	    goto cleanup;
	}
	if (!chunkedMesh.updateChunksFromCache(lod, lodHierarchy,
		{secondChunk}, lodHierarchy.max_cut,
		lodHierarchy.shaded_cull_backfaces ? TRUE : FALSE)) {
	    printf("FAIL: mesh lod second page refinement\n");
	    ret = 1;
	    goto cleanup;
	}
	uint64_t richRevision = 0;
	const std::shared_ptr<const Obol::PartGeometry> richGeometry =
	    chunkedMesh.prepareCadGeometry(
		BOBOL_LOD_DRAW_HIDDEN_LINE, &richRevision);
	const Obol::TriMesh *richMesh = richGeometry && richGeometry->shaded ?
	    &*richGeometry->shaded : NULL;
	const Obol::WireRep *richWire = richGeometry && richGeometry->wire ?
	    &*richGeometry->wire : NULL;
	if (!richMesh || !richWire || richRevision <= coarseRevision ||
	    !richWire->derivesTriangleEdges() ||
	    !richWire->triangleEdges ||
	    !richWire->segmentPoints.empty() ||
	    richMesh->progressiveLineage != coarseMesh->progressiveLineage ||
	    richWire->progressiveLineage != richMesh->progressiveLineage ||
	    !richWire->hasAdaptiveProgressiveClusters() ||
	    richWire->segmentCount() <= coarseWire->segmentCount() ||
	    richWire->segmentCount() != richMesh->indices.size() ||
	    coarseWire->segmentCount() != coarseMesh->indices.size() ||
	    !std::equal(coarseWire->triangleEdges->indices.begin(),
		coarseWire->triangleEdges->indices.end(),
		richWire->triangleEdges->indices.begin()) ||
	    richMesh->positions.size() <= coarseMesh->positions.size() ||
	    richMesh->indices.size() <= coarseMesh->indices.size() ||
	    !std::equal(coarseMesh->positions.begin(), coarseMesh->positions.end(),
		richMesh->positions.begin()) ||
	    !std::equal(coarseMesh->indices.begin(), coarseMesh->indices.end(),
		richMesh->indices.begin()) ||
	    !chunkedMesh.canDrawChunksAtCut(
		{firstChunk, secondChunk}, lodHierarchy.max_cut)) {
	    printf("FAIL: mesh lod richer page was not a suffix-only publication\n");
	    ret = 1;
	    goto cleanup;
	}

	/* Refining every page is a replacement wave, not a useful append-only
	 * suffix.  It must publish one dense stream under a new lineage; otherwise
	 * each whole-leaf cut retains its preceding mesh as tombstones and repeated
	 * Lucy-scale refinement becomes superlinear in both time and memory. */
	int allPageCoarseCut = lodHierarchy.min_cut;
	for (uint32_t chunk = 0; chunk < lodHierarchy.chunk_count; ++chunk)
	    allPageCoarseCut = std::max(
		allPageCoarseCut, lodHierarchy.chunks[chunk].min_cut);
	while (allPageCoarseCut < lodHierarchy.max_cut) {
	    bool allPopulated = true;
	    bool someRicher = false;
	    for (uint32_t chunk = 0; chunk < lodHierarchy.chunk_count; ++chunk) {
		const BObolMeshLodChunkCutInfo &coarse =
		    lodHierarchy.chunks[chunk].cuts[allPageCoarseCut];
		const BObolMeshLodChunkCutInfo &terminalPage =
		    lodHierarchy.chunks[chunk].cuts[lodHierarchy.max_cut];
		allPopulated = allPopulated && coarse.face_count > 0 &&
		    coarse.point_count > 0;
		someRicher = someRicher ||
		    coarse.face_count < terminalPage.face_count ||
		    coarse.point_count < terminalPage.point_count;
	    }
	    if (allPopulated && someRicher)
		break;
	    ++allPageCoarseCut;
	}
	BObolLodProgressiveMesh broadWaveMesh;
	if (allPageCoarseCut >= lodHierarchy.max_cut ||
	    !broadWaveMesh.updateChunksFromCache(lod, lodHierarchy,
		chunkIds, allPageCoarseCut,
		lodHierarchy.shaded_cull_backfaces ? TRUE : FALSE)) {
	    printf("FAIL: mesh lod broad-wave coarse realization\n");
	    ret = 1;
	    goto cleanup;
	}
	const std::shared_ptr<const Obol::PartGeometry> broadCoarseGeometry =
	    broadWaveMesh.prepareCadGeometry(BOBOL_LOD_DRAW_SHADED, NULL);
	const Obol::TriMesh *broadCoarseMesh =
	    broadCoarseGeometry && broadCoarseGeometry->shaded ?
		&*broadCoarseGeometry->shaded : NULL;
	if (!broadCoarseMesh || !broadCoarseMesh->progressiveLineage ||
	    !broadWaveMesh.updateChunksFromCache(lod, lodHierarchy,
		chunkIds, lodHierarchy.max_cut,
		lodHierarchy.shaded_cull_backfaces ? TRUE : FALSE)) {
	    printf("FAIL: mesh lod broad-wave terminal realization\n");
	    ret = 1;
	    goto cleanup;
	}
	const std::shared_ptr<const Obol::PartGeometry> broadRichGeometry =
	    broadWaveMesh.prepareCadGeometry(BOBOL_LOD_DRAW_SHADED, NULL);
	const Obol::TriMesh *broadRichMesh =
	    broadRichGeometry && broadRichGeometry->shaded ?
		&*broadRichGeometry->shaded : NULL;
	size_t terminalPagePoints = 0;
	for (uint32_t chunk = 0; chunk < lodHierarchy.chunk_count; ++chunk)
	    terminalPagePoints += lodHierarchy.chunks[chunk].
		cuts[lodHierarchy.max_cut].point_count;
	if (!broadRichMesh ||
	    broadRichMesh->progressiveLineage ==
		broadCoarseMesh->progressiveLineage ||
	    broadRichMesh->indices.size() !=
		static_cast<size_t>(faceCount) * 3 ||
	    broadRichMesh->positions.size() != terminalPagePoints) {
	    printf("FAIL: mesh lod broad-wave refinement retained tombstones or "
		   "reused its lineage\n");
	    ret = 1;
	    goto cleanup;
	}
	const size_t richBytes = chunkedMesh.estimateBytes();
	const BObolLodProgressiveMeshTrimPtr pageTrim =
	    chunkedMesh.prepareTrim(
		lodHierarchy.max_cut, {secondChunk});
	if (!pageTrim || !chunkedMesh.commitTrim(pageTrim)) {
	    printf("FAIL: mesh lod page-selective compaction commit\n");
	    ret = 1;
	    goto cleanup;
	}
	std::vector<uint32_t> retainedChunks;
	chunkedMesh.residentChunkIds(retainedChunks);
	uint64_t compactRevision = 0;
	const std::shared_ptr<const Obol::PartGeometry> compactGeometry =
	    chunkedMesh.prepareCadGeometry(
		BOBOL_LOD_DRAW_SHADED, &compactRevision);
	const Obol::TriMesh *compactMesh =
	    compactGeometry && compactGeometry->shaded ?
		&*compactGeometry->shaded : NULL;
	if (!compactMesh || retainedChunks != std::vector<uint32_t>{secondChunk} ||
	    compactRevision <= richRevision ||
	    compactMesh->progressiveLineage == richMesh->progressiveLineage ||
	    compactMesh->progressiveClusters.size() != 1 ||
	    chunkedMesh.estimateBytes() >= richBytes ||
	    chunkedMesh.canDrawChunksAtCut(
		{firstChunk}, lodHierarchy.max_cut) ||
	    !chunkedMesh.canDrawChunksAtCut(
		{secondChunk}, lodHierarchy.max_cut)) {
	    printf("FAIL: mesh lod page compaction retained stale ranges/data\n");
	    ret = 1;
	    goto cleanup;
	}

	double extentSquared = 0.0;
	for (int axis = 0; axis < 3; ++axis) {
	    const double extent = lodHierarchy.quantization_max[axis] -
		lodHierarchy.quantization_min[axis];
	    extentSquared += extent * extent;
	}
	const double extent = std::sqrt(extentSquared);
	const double projectedDiameter = 800.0;
	for (int cut = lodHierarchy.min_cut;
	     extent > 0.0 && cut <= lodHierarchy.max_cut; ++cut) {
	    const double error = lodHierarchy.cuts[cut].object_error;
	    if (!(error > 0.0))
		continue;
	    const double targetPixelError =
		error * projectedDiameter / extent * 1.0000001;
	    const int selected = bobol_mesh_lod_select_cut(
		&lodHierarchy, projectedDiameter, targetPixelError);
	    if (selected < lodHierarchy.min_cut || selected > cut ||
		lodHierarchy.cuts[selected].object_error >
		    error * 1.0000002) {
		printf("FAIL: PoP cut selector violated error bound "
		       "(candidate=%d selected=%d)\n", cut, selected);
		ret = 1;
		goto cleanup;
	    }
	}
	int priorSelected = lodHierarchy.min_cut;
	for (double diameter = 1.0; diameter <= 4096.0; diameter *= 2.0) {
	    const int selected = bobol_mesh_lod_select_cut(
		&lodHierarchy, diameter, 1.0);
	    if (selected < priorSelected || selected > lodHierarchy.max_cut) {
		printf("FAIL: PoP cut selection was not monotonic with projected "
		       "size (%d -> %d)\n", priorSelected, selected);
		ret = 1;
		goto cleanup;
	    }
	    priorSelected = selected;
	}
	if (bobol_mesh_lod_select_cut(NULL, 1.0, 1.0) >= 0 ||
	    bobol_mesh_lod_select_cut(&lodHierarchy, 0.0, 1.0) >= 0 ||
	    bobol_mesh_lod_select_cut(&lodHierarchy, 1.0, 0.0) >= 0) {
	    printf("FAIL: PoP cut selector accepted invalid input\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    {
	/* The production LoD service opens and reads distinct resident handles
	 * concurrently.  LMDB read transactions support that directly; guard
	 * against reintroducing unsafe cache/context state while preserving
	 * parallel warm-cache materialization. */
	const size_t readerCount = 16;
	const size_t iterations = 24;
	std::atomic<size_t> getFailures(0);
	std::atomic<size_t> loadFailures(0);
	std::atomic<size_t> dataFailures(0);
	std::atomic<size_t> countFailures(0);
	std::vector<std::thread> readers;
	readers.reserve(readerCount);
	for (size_t readerIndex = 0; readerIndex < readerCount;
	     ++readerIndex) {
	    readers.emplace_back([&]() {
		for (size_t iteration = 0; iteration < iterations;
		     ++iteration) {
		    struct BObolMeshLod *reader =
			bobol_mesh_lod_get(dbip, objname);
		    struct BObolMeshLodData readerData;
		    if (!reader) {
			getFailures.fetch_add(1, std::memory_order_relaxed);
		    } else {
			if (bobol_mesh_lod_load_cut(
				reader, lodHierarchy.max_cut, 0) !=
				lodHierarchy.max_cut) {
			    loadFailures.fetch_add(1,
				std::memory_order_relaxed);
			} else if (!bobol_mesh_lod_data_get(
				reader, &readerData)) {
			    dataFailures.fetch_add(1,
				std::memory_order_relaxed);
			} else if (readerData.face_count !=
				static_cast<size_t>(faceCount) ||
			    readerData.point_count !=
				static_cast<size_t>(vertexCount)) {
			    countFailures.fetch_add(1,
				std::memory_order_relaxed);
			}
		    }
		    if (reader)
			bobol_mesh_lod_destroy(reader);
		}
	    });
	}
	for (std::thread &reader : readers)
	    reader.join();
	const size_t getFailureCount =
	    getFailures.load(std::memory_order_relaxed);
	const size_t loadFailureCount =
	    loadFailures.load(std::memory_order_relaxed);
	const size_t dataFailureCount =
	    dataFailures.load(std::memory_order_relaxed);
	const size_t countFailureCount =
	    countFailures.load(std::memory_order_relaxed);
	if (getFailureCount || loadFailureCount || dataFailureCount ||
	    countFailureCount) {
	    printf("FAIL: concurrent mesh lod cache reads "
		"(get=%zu load=%zu data=%zu count=%zu)\n",
		getFailureCount, loadFailureCount, dataFailureCount,
		countFailureCount);
	    ret = 1;
	    goto cleanup;
	}
    }

    {
	struct BObolMeshLodData data;
	const int minimumCut =
	    bobol_mesh_lod_load_cut(lod, lodHierarchy.min_cut, 0);
	if (minimumCut < 0 || !bobol_mesh_lod_data_get(lod, &data) ||
	    data.face_count == 0) {
	    printf("FAIL: mesh lod minimum display cut was empty\n");
	    ret = 1;
	    goto cleanup;
	}

	const int terminalCut =
	    bobol_mesh_lod_load_cut(lod, lodHierarchy.max_cut, 0);
	if (terminalCut != lodHierarchy.max_cut ||
	    !bobol_mesh_lod_data_get(lod, &data) ||
	    data.face_count != static_cast<size_t>(faceCount) ||
	    data.point_count != data.point_orig_count ||
	    memcmp(data.points, data.points_orig,
		data.point_count * sizeof(point_t)) != 0) {
	    printf("FAIL: mesh lod terminal PoP cut was not exact and complete "
		   "(cut=%d faces=%zu/%d points=%zu/%zu)\n",
		   terminalCut, data.face_count, faceCount,
		   data.point_count, data.point_orig_count);
	    ret = 1;
	    goto cleanup;
	}
	struct BObolMeshLodInfo surfaceInfo = BOBOL_MESH_LOD_INFO_INIT;
	if (!bobol_mesh_lod_info_get(lod, &surfaceInfo) ||
	    surfaceInfo.shaded_cull_backfaces) {
	    printf("FAIL: open surface was marked safe for backface culling\n");
	    ret = 1;
	    goto cleanup;
	}

	struct BObolMeshLodHierarchyInfo suffixHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	if (!bobol_mesh_lod_hierarchy_info_get(lod, &suffixHierarchy) ||
	    suffixHierarchy.max_cut <= suffixHierarchy.min_cut) {
	    printf("FAIL: mesh lod suffix hierarchy unavailable\n");
	    ret = 1;
	    goto cleanup;
	}
	const int cutBeforeSuffix = bobol_mesh_lod_current_cut(lod);
	bobol_mesh_lod_memshrink(lod);
	mesh_lod_suffix_test_data suffixTest;
	suffixTest.hierarchy = &suffixHierarchy;
	suffixTest.nextCut = suffixHierarchy.min_cut + 1;
	suffixTest.cumulativePoints =
	    suffixHierarchy.cuts[suffixHierarchy.min_cut].point_count;
	if (bobol_mesh_lod_resident_prefix_bytes(lod) != 0 ||
	    !bobol_mesh_lod_read_resident_suffix(lod,
		suffixHierarchy.min_cut, suffixHierarchy.max_cut,
		mesh_lod_suffix_test_callback, &suffixTest) ||
	    suffixTest.nextCut != suffixHierarchy.max_cut + 1 ||
	    suffixTest.streamedPoints !=
		suffixHierarchy.cuts[suffixHierarchy.max_cut].point_count -
		    suffixHierarchy.cuts[suffixHierarchy.min_cut].point_count ||
	    suffixTest.streamedFaces !=
		suffixHierarchy.cuts[suffixHierarchy.max_cut].face_count -
		    suffixHierarchy.cuts[suffixHierarchy.min_cut].face_count ||
	    bobol_mesh_lod_current_cut(lod) != cutBeforeSuffix ||
	    bobol_mesh_lod_resident_prefix_bytes(lod) != 0) {
	    printf("FAIL: mesh lod suffix stream materialized or omitted cache "
		   "records (cuts=%d/%d points=%zu faces=%zu bytes=%zu)\n",
		   suffixTest.nextCut, suffixHierarchy.max_cut + 1,
		   suffixTest.streamedPoints, suffixTest.streamedFaces,
		   bobol_mesh_lod_resident_prefix_bytes(lod));
	    ret = 1;
	    goto cleanup;
	}
    }

    {
	struct BObolMeshLod *translatedLod = NULL;
	struct BObolMeshLodHierarchyInfo baseHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	struct BObolMeshLodHierarchyInfo translatedHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	struct BObolMeshLodData translatedData;
	if (!bobol_mesh_lod_hierarchy_info_get(lod, &baseHierarchy) ||
	    bobol_mesh_lod_cache_store_mesh(dbip, translatedMeshObjname,
		reinterpret_cast<const point_t *>(translatedVertices),
		static_cast<size_t>(vertexCount), NULL, faces,
		static_cast<size_t>(faceCount), 515151ULL, 0,
		&cacheStatus) != BRLCAD_OK ||
	    !(translatedLod = bobol_mesh_lod_get(
		dbip, translatedMeshObjname)) ||
	    !bobol_mesh_lod_hierarchy_info_get(
		translatedLod, &translatedHierarchy)) {
	    printf("FAIL: translated mesh PoP hierarchy setup\n");
	    if (translatedLod)
		bobol_mesh_lod_destroy(translatedLod);
	    ret = 1;
	    goto cleanup;
	}
	for (int cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut) {
	    if (baseHierarchy.cuts[cut].face_count !=
		    translatedHierarchy.cuts[cut].face_count ||
		baseHierarchy.cuts[cut].point_count !=
		    translatedHierarchy.cuts[cut].point_count) {
		printf("FAIL: PoP hierarchy depends on absolute model coordinates "
		       "at cut %d\n", cut);
		ret = 1;
		break;
	    }
	}
	for (int axis = 0; !ret && axis < 3; axis++) {
	    const fastf_t baseExtent =
		baseHierarchy.quantization_max[axis] -
		baseHierarchy.quantization_min[axis];
	    const fastf_t translatedExtent =
		translatedHierarchy.quantization_max[axis] -
		translatedHierarchy.quantization_min[axis];
	    if (!NEAR_EQUAL(baseExtent, translatedExtent, 0.01)) {
		printf("FAIL: translated PoP quantization extent changed on axis %d "
		       "(%.17g != %.17g)\n", axis,
		       baseExtent, translatedExtent);
		ret = 1;
	    }
	}
	if (!ret &&
	    (bobol_mesh_lod_load_cut(
		translatedLod, translatedHierarchy.max_cut, 0) !=
		translatedHierarchy.max_cut ||
	     !bobol_mesh_lod_data_get(translatedLod, &translatedData) ||
	     translatedData.face_count != static_cast<size_t>(faceCount) ||
	     translatedData.point_count !=
		static_cast<size_t>(vertexCount) ||
	     translatedData.point_orig_count != translatedData.point_count ||
	     memcmp(translatedData.points, translatedData.points_orig,
		static_cast<size_t>(vertexCount) * sizeof(point_t)) != 0)) {
	    printf("FAIL: translated terminal PoP cut was not exact\n");
	    ret = 1;
	}
	bobol_mesh_lod_destroy(translatedLod);
	if (ret)
	    goto cleanup;
    }

    {
	struct BObolMeshLod *solidLod = NULL;
	struct BObolMeshLodInfo solidInfo = BOBOL_MESH_LOD_INFO_INIT;
	struct BObolMeshLodHierarchyInfo solidHierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	if (bobol_mesh_lod_cache_refresh(dbip, solidCwObjname,
		&cacheStatus) != BRLCAD_OK ||
	    !(solidLod = bobol_mesh_lod_get(dbip, solidCwObjname)) ||
	    load_terminal_cut(solidLod) < 0 ||
	    !bobol_mesh_lod_info_get(solidLod, &solidInfo) ||
	    !bobol_mesh_lod_hierarchy_info_get(solidLod, &solidHierarchy) ||
	    !solidInfo.shaded_cull_backfaces ||
	    solidInfo.face_count != 4 ||
	    solidHierarchy.cluster_grid_resolution != 0 ||
	    solidHierarchy.cluster_count != 0 || solidHierarchy.clusters ||
	    solidHierarchy.chunk_count != 0 || solidHierarchy.chunks) {
	    printf("FAIL: closed oriented CW BoT culling metadata\n");
	    if (solidLod)
		bobol_mesh_lod_destroy(solidLod);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(solidLod);
    }

    {
	struct BObolMeshLod *solidLod = NULL;
	struct BObolMeshLodInfo solidInfo = BOBOL_MESH_LOD_INFO_INIT;
	struct BObolMeshLodData solidData;
	if (bobol_mesh_lod_cache_refresh(dbip, solidUnorientedObjname,
		&cacheStatus) != BRLCAD_OK ||
	    !(solidLod = bobol_mesh_lod_get(
		dbip, solidUnorientedObjname)) ||
	    load_terminal_cut(solidLod) < 0 ||
	    !bobol_mesh_lod_info_get(solidLod, &solidInfo) ||
	    !bobol_mesh_lod_data_get(solidLod, &solidData) ||
	    solidInfo.shaded_cull_backfaces ||
	    solidInfo.face_count != 4 ||
	    mesh_signed_six_volume(solidData) >= 0.0L) {
	    printf("FAIL: unoriented BoT did not preserve two-sided authored "
		   "winding metadata "
		   "(cull=%d info_faces=%zu data_faces=%zu volume=%.17Lg)\n",
		   solidInfo.shaded_cull_backfaces,
		   solidInfo.face_count, solidData.face_count,
		   mesh_signed_six_volume(solidData));
	    if (solidLod)
		bobol_mesh_lod_destroy(solidLod);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(solidLod);
    }

    {
	struct BObolMeshLod *brokenLod = NULL;
	struct BObolMeshLodInfo brokenInfo = BOBOL_MESH_LOD_INFO_INIT;
	if (bobol_mesh_lod_cache_refresh(dbip, brokenUnorientedObjname,
		&cacheStatus) != BRLCAD_OK ||
	    !(brokenLod = bobol_mesh_lod_get(
		dbip, brokenUnorientedObjname)) ||
	    load_terminal_cut(brokenLod) < 0 ||
	    !bobol_mesh_lod_info_get(brokenLod, &brokenInfo) ||
	    brokenInfo.shaded_cull_backfaces) {
	    printf("FAIL: inconsistently wound unoriented BoT was marked "
		   "safe for backface culling\n");
	    if (brokenLod)
		bobol_mesh_lod_destroy(brokenLod);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(brokenLod);
    }

    memshrinkCut = first_available_cut(lod);
    if (memshrinkCut < 0) {
	printf("FAIL: mesh lod cut data\n");
	ret = 1;
	goto cleanup;
    }

    {
	if (bobol_mesh_lod_load_cut(lod, lodHierarchy.min_cut, 0) < 0 ||
	    check_mesh_lod_payload("mesh lod explicit cut", lod,
				   static_cast<size_t>(faceCount),
				   static_cast<size_t>(vertexCount), 0)) {
	    ret = 1;
	    goto cleanup;
	}
    }

    {
	struct BObolMeshLod *meshLod = NULL;
	int invalidStoreFaces[3] = {0, vertexCount, 1};
	if (bobol_mesh_lod_cache_store_mesh(dbip, meshObjname,
					      reinterpret_cast<const point_t *>(vertices),
					      static_cast<size_t>(vertexCount),
					      reinterpret_cast<const vect_t *>(detailNormals),
					      faces, static_cast<size_t>(faceCount), 424242ULL,
					      0, &cacheStatus) != BRLCAD_OK ||
	    !cacheStatus.directory_found || cacheStatus.is_bot ||
	    !cacheStatus.has_cache_key || !cacheStatus.has_cached_payload ||
	    cacheStatus.stale_cache_entry ||
	    !cacheStatus.generated_cache_entry || !cacheStatus.cache_key) {
	    printf("FAIL: mesh lod store generated mesh status\n");
	    ret = 1;
	    goto cleanup;
	}
	meshLod = bobol_mesh_lod_get(dbip, meshObjname);
	if (!meshLod || first_available_cut(meshLod) < 0 ||
	    check_mesh_lod_payload("mesh lod generated mesh", meshLod,
				   static_cast<size_t>(faceCount),
				   static_cast<size_t>(vertexCount), 1)) {
	    printf("FAIL: mesh lod generated mesh data\n");
	    ret = 1;
	}
	if (meshLod)
	    bobol_mesh_lod_destroy(meshLod);
	if (ret)
	    goto cleanup;

	if (bobol_mesh_lod_cache_store_mesh(dbip, meshObjname,
					      reinterpret_cast<const point_t *>(vertices),
					      static_cast<size_t>(vertexCount), NULL, invalidStoreFaces, 1,
					      777777ULL, 0, &cacheStatus) != BRLCAD_ERROR ||
	    !cacheStatus.has_cache_key ||
	    !cacheStatus.has_cached_payload ||
	    cacheStatus.cache_key != 424242ULL) {
	    printf("FAIL: mesh lod invalid generated mesh store preservation\n");
	    ret = 1;
	    goto cleanup;
	}

	struct BObolMeshLod *firstVariant = NULL;
	struct BObolMeshLod *secondVariant = NULL;
	if (bobol_mesh_lod_cache_store_mesh_variant(dbip, meshObjname,
		reinterpret_cast<const point_t *>(vertices),
		static_cast<size_t>(vertexCount),
		reinterpret_cast<const vect_t *>(detailNormals), faces,
		static_cast<size_t>(faceCount), 525252ULL, 0,
		&cacheStatus) != BRLCAD_OK ||
	    cacheStatus.cache_key != 525252ULL ||
	    !(firstVariant = bobol_mesh_lod_get_cached_prefix(
		dbip, 424242ULL)) ||
	    !(secondVariant = bobol_mesh_lod_get_cached_prefix(
		dbip, 525252ULL)) ||
	    load_terminal_cut(firstVariant) < 0 ||
	    load_terminal_cut(secondVariant) < 0) {
	    printf("FAIL: mesh lod representation variants did not coexist\n");
	    if (firstVariant)
		bobol_mesh_lod_destroy(firstVariant);
	    if (secondVariant)
		bobol_mesh_lod_destroy(secondVariant);
	    ret = 1;
	    goto cleanup;
	}
	bobol_mesh_lod_destroy(firstVariant);
	bobol_mesh_lod_destroy(secondVariant);
    }

    {
	int activeCut;
	struct BObolMeshLodData shrinkData;
	struct BObolMeshLodInfo shrinkInfo = BOBOL_MESH_LOD_INFO_INIT;
	if (bobol_mesh_lod_load_cut(lod, memshrinkCut, 0) !=
	    memshrinkCut ||
	    (activeCut = bobol_mesh_lod_current_cut(lod)) < 0 ||
	    !bobol_mesh_lod_has_active_data(lod)) {
	    printf("FAIL: mesh lod memshrink setup\n");
	    ret = 1;
	    goto cleanup;
	}
	const size_t residentBytesBefore =
	    bobol_mesh_lod_resident_bytes(lod);
	const size_t prefixBytesBefore =
	    bobol_mesh_lod_resident_prefix_bytes(lod);
	bobol_mesh_lod_memshrink(lod);
	if (bobol_mesh_lod_current_cut(lod) != activeCut ||
	    bobol_mesh_lod_has_active_data(lod) ||
	    bobol_mesh_lod_data_get(lod, &shrinkData) ||
	    bobol_mesh_lod_info_get(lod, &shrinkInfo) ||
	    shrinkInfo.active_cut != activeCut ||
	    !prefixBytesBefore ||
	    bobol_mesh_lod_resident_prefix_bytes(lod) != 0 ||
	    bobol_mesh_lod_resident_bytes(lod) >= residentBytesBefore) {
	    printf("FAIL: mesh lod memshrink status\n");
	    ret = 1;
	    goto cleanup;
	}
	if (bobol_mesh_lod_load_cut(lod, activeCut, 0) != activeCut ||
	    check_mesh_lod_payload("mesh lod reload after memshrink", lod,
				   static_cast<size_t>(faceCount),
				   static_cast<size_t>(vertexCount), 0) ||
	    bobol_mesh_lod_resident_prefix_bytes(lod) == 0) {
	    printf("FAIL: mesh lod reload after memshrink\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    bobol_mesh_lod_destroy(lod);
    lod = NULL;
    db_close(dbip);
    dbip = db_open(dbpath, DB_OPEN_READWRITE);
    if (!dbip || db_dirbuild(dbip) < 0) {
	printf("FAIL: mesh lod reopen\n");
	ret = 1;
	goto cleanup;
    }

    lod = bobol_mesh_lod_get(dbip, objname);
    if (!lod || first_available_cut(lod) < 0 ||
	check_mesh_lod_payload("mesh lod reopen", lod,
			       static_cast<size_t>(faceCount),
			       static_cast<size_t>(vertexCount), 0)) {
	printf("FAIL: mesh lod reopen data\n");
	ret = 1;
	goto cleanup;
    }

    bobol_mesh_lod_destroy(lod);
    lod = NULL;

    if (bobol_mesh_lod_cache_invalidate(dbip, objname, &cacheStatus) !=
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

    if (bobol_mesh_lod_get(dbip, objname)) {
	printf("FAIL: mesh lod get after invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_mesh_lod_cache_update(dbip, objname) != BRLCAD_OK) {
	printf("FAIL: mesh lod compatibility update\n");
	ret = 1;
	goto cleanup;
    }

    lod = bobol_mesh_lod_get(dbip, objname);
    if (!lod || first_available_cut(lod) < 0 ||
	check_mesh_lod_payload("mesh lod update after invalidate", lod,
			       static_cast<size_t>(faceCount),
			       static_cast<size_t>(vertexCount), 0)) {
	printf("FAIL: mesh lod compatibility update data\n");
	ret = 1;
	goto cleanup;
    }

    {
	struct BObolMeshLodData nullData;
	struct BObolMeshLodInfo nullInfo;
	nullInfo.active_cut = 7;
	bobol_mesh_lod_cache_status_init(&cacheStatus);
	cacheStatus.cache_key = 99;
	bobol_mesh_lod_cache_status_init(&cacheStatus);
	if (bobol_mesh_lod_load_cut(NULL, 0, 0) >= 0 ||
	    bobol_mesh_lod_current_cut(NULL) >= 0 ||
	    bobol_mesh_lod_has_active_data(NULL) ||
	    bobol_mesh_lod_data_get(NULL, &nullData) ||
	    bobol_mesh_lod_info_get(NULL, &nullInfo) ||
	    nullInfo.active_cut != -1 ||
	    cacheStatus.cache_key ||
	    bobol_mesh_lod_cache_status(NULL, objname,
					  &cacheStatus) != BRLCAD_ERROR ||
	    bobol_mesh_lod_cache_status(dbip, NULL,
					  &cacheStatus) != BRLCAD_ERROR ||
	    bobol_mesh_lod_cache_status(dbip, objname,
					  NULL) != BRLCAD_ERROR ||
	    bobol_mesh_lod_cache_invalidate(NULL, objname,
					      &cacheStatus) != BRLCAD_ERROR) {
	    printf("FAIL: mesh lod null handling\n");
	    ret = 1;
	}
	bobol_mesh_lod_memshrink(NULL);
    }

cleanup:
    if (lod)
	bobol_mesh_lod_destroy(lod);
    if (dbip) {
	bobol_mesh_lod_cache_clear_database(dbip);
	db_close(dbip);
    }
    if (dbpath[0])
	bu_file_delete(dbpath);
    bu_dirclear(cacheDir);
    if (vertices)
	bu_free(vertices, "mesh lod vertices");
    if (translatedVertices)
	bu_free(translatedVertices, "translated mesh lod vertices");
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
