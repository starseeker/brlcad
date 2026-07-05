/*                    T E S T _ D R A W _ C A C H E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 *
 * This library is free software; you can redistribute it and/or
 * modify it under the terms of the GNU Lesser General Public License
 * version 2.1 as published by the Free Software Foundation.
 */

#include "common.h"

#include "brlobol/draw_cache.h"
#include "brlobol/lod_realization.h"

#include "bu/app.h"
#include "bu/env.h"
#include "bu/file.h"
#include "raytrace.h"
#include "rt/db_attr.h"
#include "rt/wdb.h"
#include "wdb.h"

#include <cmath>
#include <cstdio>
#include <cstring>

static const char *path_leaf_name = "path_leaf.s";
static const char *path_region_name = "path_region.r";
static const char *path_top_name = "path_top.c";
static const char *path_region_full_name = "path_top.c/path_region.r";
static const char *path_full_name = "path_top.c/path_region.r/path_leaf.s";

static int
fastf_equal(fastf_t a, fastf_t b)
{
    return std::fabs(a - b) <= SMALL_FASTF;
}

static int
point_equal(const point_t a, const point_t b)
{
    return fastf_equal(a[X], b[X]) && fastf_equal(a[Y], b[Y]) &&
	   fastf_equal(a[Z], b[Z]);
}

static int
check_proxy_point(const char *label,
		  const BRLObolDrawProxyRecord *record,
		  size_t index,
		  const point_t expected)
{
    if (!record || index >= record->pointCount ||
	!point_equal(record->points[index], expected)) {
	printf("FAIL: %s proxy point %zu\n", label, index);
	return 1;
    }
    return 0;
}

static int
make_path_metadata_tree(rt_wdb *wdbp)
{
    point_t path_center;
    unsigned char region_rgb[3] = {10, 20, 30};
    unsigned char top_rgb[3] = {120, 130, 140};
    struct wmember region_members;
    struct wmember top_members;

    if (!wdbp)
	return 1;

    VSET(path_center, 30.0, 0.0, 0.0);
    if (mk_sph(wdbp, path_leaf_name, path_center, 2.0) < 0)
	return 1;

    BU_LIST_INIT(&region_members.l);
    if (!mk_addmember(path_leaf_name, &region_members.l, NULL, WMOP_UNION))
	return 1;
    if (mk_comb(wdbp, path_region_name, &region_members.l, 1,
		"phong", NULL, region_rgb, 88, 5, 66, 101,
		0, 0, 0) < 0)
	return 1;

    BU_LIST_INIT(&top_members.l);
    if (!mk_addmember(path_region_name, &top_members.l, NULL, WMOP_UNION))
	return 1;
    if (mk_comb(wdbp, path_top_name, &top_members.l, 0,
		"plastic", NULL, top_rgb, 0, 0, 0, 0,
		1, 0, 0) < 0)
	return 1;

    return 0;
}

static int
set_test_attributes(db_i *dbip, const char *name)
{
    directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    bu_attribute_value_set avs = BU_AVS_INIT_ZERO;

    if (dp == RT_DIR_NULL)
	return 1;

    bu_avs_init_empty(&avs);
    bu_avs_add(&avs, "region", "yes");
    bu_avs_add(&avs, "region_id", "77");
    bu_avs_add(&avs, "air", "2");
    bu_avs_add(&avs, "los", "100");
    bu_avs_add(&avs, "material_id", "42");
    bu_avs_add(&avs, "inherit", "1");
    bu_avs_add(&avs, "color", "10/20/30");
    bu_avs_add(&avs, "material_name", "DrawCacheMaterial");
    bu_avs_add(&avs, "shader", "plastic {di .7 sp .2}");

    if (db5_update_attributes(dp, &avs, dbip) != 0) {
	bu_avs_free(&avs);
	return 1;
    }

    bu_avs_free(&avs);
    return 0;
}

static int
check_plain_leaf_metadata(const BRLObolDrawMetadataRecord *record)
{
    if (!record || !record->directoryFound || !record->isSolid ||
	record->isComb || record->isRegion || record->hasRegionId ||
	record->hasAircode || record->hasLos || record->hasMaterialId ||
	record->hasColor || record->hasShader) {
	printf("FAIL: plain leaf metadata fields\n");
	return 1;
    }

    return 0;
}

static int
check_metadata(const BRLObolDrawMetadataRecord *record)
{
    if (!record || !record->directoryFound || !record->isSolid ||
	record->isComb || !record->isRegion || record->isHidden ||
	!record->hasRegionId || record->regionId != 77 ||
	!record->hasAircode || record->aircode != 2 ||
	!record->hasLos || record->los != 100 ||
	!record->hasMaterialId || record->materialId != 42 ||
	!record->hasInherit || record->inherit != 1 ||
	!record->hasColor || record->color[0] != 10 ||
	record->color[1] != 20 || record->color[2] != 30 ||
	!record->hasMaterialName ||
	std::strcmp(record->materialName, "DrawCacheMaterial") != 0 ||
	!record->hasShader ||
	std::strcmp(record->shader, "plastic {di .7 sp .2}") != 0) {
	printf("FAIL: draw metadata fields\n");
	return 1;
    }

    return 0;
}

static int
check_path_metadata(const BRLObolDrawMetadataRecord *record)
{
    if (!record || !record->directoryFound || !record->isSolid ||
	record->isComb || !record->isRegion || record->isHidden ||
	!record->hasRegionId || record->regionId != 88 ||
	!record->hasAircode || record->aircode != 5 ||
	!record->hasLos || record->los != 101 ||
	!record->hasMaterialId || record->materialId != 66 ||
	!record->hasInherit || record->inherit != 1 ||
	!record->hasColor || record->color[0] != 120 ||
	record->color[1] != 130 || record->color[2] != 140 ||
	!record->hasShader || !std::strstr(record->shader, "plastic")) {
	printf("FAIL: path draw metadata fields\n");
	return 1;
    }

    return 0;
}

int
main(int argc, char *argv[])
{
    const char *objname = "draw_cache_sph.s";
    char dbpath[MAXPATHLEN] = {0};
    char cacheDir[MAXPATHLEN] = {0};
    db_i *dbip = NULL;
    int ret = 0;
    BRLObolDrawCacheStatus status;
    BRLObolDrawProxyRecord proxy;
    BRLObolDrawMetadataRecord metadata;
    point_t center;
    point_t aabbMin;
    point_t aabbMax;
    point_t obbPoints[8];

    bu_setprogname(argv[0]);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    bu_dir(cacheDir, MAXPATHLEN, BU_DIR_CURR,
	   "brlobol_draw_cache_test", NULL);
    bu_mkdir(cacheDir);
    bu_setenv("BU_DIR_CACHE", cacheDir, 1);

    VSET(center, 1.0, 2.0, 3.0);
    VSET(aabbMin, -3.0, -2.0, -1.0);
    VSET(aabbMax, 5.0, 6.0, 7.0);
    for (size_t i = 0; i < 8; i++)
	VSET(obbPoints[i], (fastf_t)i, (fastf_t)(i + 10),
	     (fastf_t)(i + 20));

    {
	FILE *fp = bu_temp_file(dbpath, MAXPATHLEN);
	if (!fp) {
	    printf("FAIL: draw cache temp file\n");
	    ret = 1;
	    goto cleanup;
	}
	fclose(fp);
    }

    dbip = db_create(dbpath, 5);
    if (!dbip) {
	printf("FAIL: draw cache db_create\n");
	ret = 1;
	goto cleanup;
    }

    {
	rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
	if (!wdbp || mk_sph(wdbp, objname, center, 4.0) < 0 ||
	    make_path_metadata_tree(wdbp)) {
	    printf("FAIL: draw cache mk_sph\n");
	    ret = 1;
	    goto cleanup;
	}
    }

    if (set_test_attributes(dbip, objname)) {
	printf("FAIL: draw cache attributes\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_status(dbip, objname,
					BRLOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || status.hasCachedPayload) {
	printf("FAIL: draw cache initial AABB status\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_refresh(dbip, objname,
					 BRLOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry) {
	printf("FAIL: draw cache AABB refresh\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_AABB,
				     &proxy) != BRLCAD_OK || proxy.kind != BRLOBOL_LOD_PROXY_AABB ||
	proxy.pointCount != 2 ||
	check_proxy_point("AABB min", &proxy, 0, aabbMin) ||
	check_proxy_point("AABB max", &proxy, 1, aabbMax)) {
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_store(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				       obbPoints, 8, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry) {
	printf("FAIL: draw cache OBB store\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK || proxy.kind != BRLOBOL_LOD_PROXY_OBB ||
	proxy.pointCount != 8 ||
	check_proxy_point("OBB first", &proxy, 0, obbPoints[0]) ||
	check_proxy_point("OBB last", &proxy, 7, obbPoints[7])) {
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_metadata_cache_refresh(dbip, objname, &status) !=
	BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry ||
	brlobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_OK ||
	check_metadata(&metadata)) {
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_metadata_cache_refresh(dbip, path_leaf_name, &status) !=
	BRLCAD_OK ||
	brlobol_draw_metadata_cache_get(dbip, path_leaf_name, &metadata) !=
	BRLCAD_OK ||
	check_plain_leaf_metadata(&metadata)) {
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_status(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || status.hasCachedPayload) {
	printf("FAIL: draw cache initial path metadata status\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK ||
	check_path_metadata(&metadata)) {
	ret = 1;
	goto cleanup;
    }

    db_close(dbip);
    dbip = db_open(dbpath, DB_OPEN_READWRITE);
    if (!dbip || db_dirbuild(dbip) < 0) {
	printf("FAIL: draw cache reopen\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK ||
	check_proxy_point("OBB reopen", &proxy, 7, obbPoints[7]) ||
	brlobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_OK ||
	check_metadata(&metadata)) {
	printf("FAIL: draw cache persistent context reopen data\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
	    &metadata) != BRLCAD_OK ||
	check_path_metadata(&metadata)) {
	printf("FAIL: draw cache persistent path metadata reopen data\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_invalidate(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
		&status) != BRLCAD_OK ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK ||
	check_path_metadata(&metadata)) {
	printf("FAIL: draw cache path metadata invalidate/refresh\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_refresh(dbip, path_region_full_name,
	    &status) != BRLCAD_OK ||
	brlobol_draw_path_metadata_cache_get(dbip, path_region_full_name,
		&metadata) != BRLCAD_OK ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK) {
	printf("FAIL: draw cache path metadata object-invalidate setup\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_invalidate_object(dbip,
	    path_leaf_name, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_get(dbip, path_region_full_name,
		&metadata) != BRLCAD_OK) {
	printf("FAIL: draw cache path metadata leaf object invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	brlobol_draw_path_metadata_cache_invalidate_object(dbip,
		path_region_name, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_get(dbip, path_region_full_name,
		&metadata) != BRLCAD_ERROR ||
	brlobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_OK) {
	printf("FAIL: draw cache path metadata ancestor object invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK) {
	printf("FAIL: draw cache path metadata object-invalidate reseed\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_invalidate(dbip, objname,
					    BRLOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_AABB,
				     &proxy) != BRLCAD_ERROR ||
	brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK) {
	printf("FAIL: draw cache targeted invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_cache_clear_database(dbip) != BRLCAD_OK ||
	brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_ERROR ||
	brlobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) !=
	BRLCAD_ERROR) {
	printf("FAIL: draw cache database clear\n");
	ret = 1;
	goto cleanup;
    }

    if (brlobol_draw_proxy_cache_store(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				       obbPoints, 8, &status) != BRLCAD_OK ||
	brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK) {
	printf("FAIL: draw cache reuse after database clear\n");
	ret = 1;
	goto cleanup;
    }

    brlobol_draw_cache_clear_all();
    if (brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_ERROR ||
	brlobol_draw_proxy_cache_store(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				       obbPoints, 8, &status) != BRLCAD_OK ||
	brlobol_draw_proxy_cache_get(dbip, objname, BRLOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK) {
	printf("FAIL: draw cache clear-all context reset\n");
	ret = 1;
	goto cleanup;
    }

    brlobol_draw_cache_status_init(&status);
    status.generatedCacheEntry = 1;
    brlobol_draw_cache_status_init(&status);
    brlobol_draw_proxy_record_init(&proxy);
    proxy.kind = 9;
    brlobol_draw_proxy_record_init(&proxy);
    brlobol_draw_metadata_record_init(&metadata);
    metadata.directoryFound = 1;
    brlobol_draw_metadata_record_init(&metadata);
    if (status.generatedCacheEntry || proxy.kind || metadata.directoryFound ||
	brlobol_draw_proxy_cache_status(NULL, objname,
					BRLOBOL_LOD_PROXY_AABB, &status) != BRLCAD_ERROR ||
	brlobol_draw_proxy_cache_status(dbip, NULL,
					BRLOBOL_LOD_PROXY_AABB, &status) != BRLCAD_ERROR ||
	brlobol_draw_proxy_cache_status(dbip, objname, 999,
					&status) != BRLCAD_ERROR ||
	brlobol_draw_proxy_cache_status(dbip, objname,
					BRLOBOL_LOD_PROXY_AABB, NULL) != BRLCAD_ERROR ||
	brlobol_draw_proxy_cache_store(dbip, objname,
				       BRLOBOL_LOD_PROXY_AABB, NULL, 2, &status) != BRLCAD_ERROR ||
	brlobol_draw_metadata_cache_status(NULL, objname,
					   &status) != BRLCAD_ERROR ||
	brlobol_draw_metadata_cache_status(dbip, NULL,
					   &status) != BRLCAD_ERROR ||
	brlobol_draw_metadata_cache_status(dbip, objname,
					   NULL) != BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_status(NULL, path_full_name,
		&status) != BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_status(dbip, NULL,
		&status) != BRLCAD_ERROR ||
	brlobol_draw_path_metadata_cache_status(dbip, path_full_name,
		NULL) != BRLCAD_ERROR) {
	printf("FAIL: draw cache null handling\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    brlobol_draw_cache_clear_all();
    if (dbip)
	db_close(dbip);
    if (dbpath[0])
	bu_file_delete(dbpath);
    if (cacheDir[0])
	bu_dirclear(cacheDir);

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
