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

#include "BObol/BMeshLodCache.h"

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

struct mesh_lod_suffix_test_data {
    const struct BObolMeshLodHierarchyInfo *hierarchy = NULL;
    int nextCut = -1;
    size_t cumulativePoints = 0;
    size_t streamedPoints = 0;
    size_t streamedFaces = 0;
};

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
    const int grid = 12;
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
	RT_DB_INTERNAL_INIT(&intern);
	if (!dp || rt_db_get_internal(&intern, dp, dbip, NULL) < 0 ||
	    intern.idb_type != ID_BOT || !intern.idb_ptr ||
	    !(opened = bobol_mesh_lod_cache_refresh_from_bot_open(
		dbip, objname,
		static_cast<const struct rt_bot_internal *>(intern.idb_ptr),
		&cacheStatus)) ||
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
	    !(opened = bobol_mesh_lod_cache_refresh_from_bot_open(
		dbip, duplicateObjname,
		static_cast<const struct rt_bot_internal *>(intern.idb_ptr),
		&cacheStatus)) ||
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
	if (bobol_mesh_lod_cache_refresh(dbip, solidCwObjname,
		&cacheStatus) != BRLCAD_OK ||
	    !(solidLod = bobol_mesh_lod_get(dbip, solidCwObjname)) ||
	    load_terminal_cut(solidLod) < 0 ||
	    !bobol_mesh_lod_info_get(solidLod, &solidInfo) ||
	    !solidInfo.shaded_cull_backfaces ||
	    solidInfo.face_count != 4) {
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
