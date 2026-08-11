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

#include "bu/str.h"

#include "BObol/BDrawCache.h"
#include "BObol/BLodRealization.h"

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
static const char *path_table_leaf_name = "path_table_leaf.s";
static const char *path_table_region_name = "path_table_region.r";
static const char *path_table_top_name = "path_table_top.c";
static const char *path_table_full_name =
    "path_table_top.c/path_table_region.r/path_table_leaf.s";
static const char *lod_asset_name = "lod_asset.bot";
static const char *lod_copy_name = "lod_copy.bot";

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
		  const BObolDrawProxyRecord *record,
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
check_manifest(const BObolDrawManifest *manifest)
{
    if (!manifest || manifest->occurrenceCount != 2 ||
	!manifest->coverageBoundsValid ||
	!fastf_equal(manifest->coverageBoundsMin[X], -100.0) ||
	!fastf_equal(manifest->coverageBoundsMin[Y], -200.0) ||
	!fastf_equal(manifest->coverageBoundsMin[Z], -300.0) ||
	!fastf_equal(manifest->coverageBoundsMax[X], 400.0) ||
	!fastf_equal(manifest->coverageBoundsMax[Y], 500.0) ||
	!fastf_equal(manifest->coverageBoundsMax[Z], 600.0) ||
	!manifest->occurrences || !manifest->occurrences[0].path ||
	!manifest->occurrences[0].sourceName ||
	!manifest->occurrences[1].path ||
	!manifest->occurrences[1].sourceName ||
	bu_strcmp(manifest->occurrences[0].path, path_full_name) != 0 ||
	bu_strcmp(manifest->occurrences[0].sourceName, path_leaf_name) != 0 ||
	bu_strcmp(manifest->occurrences[1].path, path_region_full_name) != 0 ||
	bu_strcmp(manifest->occurrences[1].sourceName, path_region_name) != 0 ||
	!manifest->occurrences[0].sourceMeshRequestValid ||
	manifest->occurrences[0].meshAssetKind !=
	    BOBOL_DRAW_CACHE_MESH_ASSET_BREP ||
	manifest->occurrences[0].meshAssetContentHash != 0x123456789abcdef0ULL ||
	!fastf_equal(manifest->occurrences[0].meshAssetTessellationAbsTol,
	    0.125) ||
	!fastf_equal(manifest->occurrences[0].meshAssetTessellationRelTol,
	    0.01) ||
	!fastf_equal(manifest->occurrences[0].meshAssetTessellationNormTol,
	    0.5) ||
	!manifest->occurrences[0].meshAssetPath ||
	!manifest->occurrences[0].meshAssetName ||
	bu_strcmp(manifest->occurrences[0].meshAssetPath,
	    path_full_name) != 0 ||
	bu_strcmp(manifest->occurrences[0].meshAssetName,
	    path_leaf_name) != 0 ||
	manifest->occurrences[0].sourceFaceCount != 123456 ||
	manifest->occurrences[0].sourcePointCount != 65432 ||
	!fastf_equal(
	    manifest->occurrences[0].meshAssetBoundsMin[X], -10.0) ||
	!fastf_equal(
	    manifest->occurrences[0].meshAssetBoundsMax[Z], 60.0) ||
	!fastf_equal(
	    manifest->occurrences[0].meshAssetMatrix[MDX], 2.0) ||
	!fastf_equal(
	    manifest->occurrences[0].meshAssetMatrix[MDZ], 8.0) ||
	manifest->occurrences[1].sourceMeshRequestValid ||
	manifest->occurrences[0].booleanOperation != DB_OP_UNION ||
	manifest->occurrences[1].booleanOperation != DB_OP_SUBTRACT ||
	manifest->occurrences[0].occurrenceIndex != 3 ||
	manifest->occurrences[1].occurrenceIndex != 7 ||
	!manifest->occurrences[0].metadataValid ||
	!manifest->occurrences[0].metadata.hasColor ||
	manifest->occurrences[0].metadata.color[0] != 10 ||
	manifest->occurrences[0].metadata.color[1] != 20 ||
	manifest->occurrences[0].metadata.color[2] != 30 ||
	!fastf_equal(manifest->occurrences[0].localMatrix[MDX], 11.0) ||
	!fastf_equal(manifest->occurrences[1].localMatrix[MDY], 13.0) ||
	!fastf_equal(manifest->occurrences[0].boundsMin[X], -1.0) ||
	!fastf_equal(manifest->occurrences[1].boundsMax[Z], 9.0)) {
	printf("FAIL: draw manifest data\n");
	return 1;
    }
    return 0;
}

static int
make_manifest(BObolDrawManifest *manifest)
{
    if (!manifest)
	return 0;
    bobol_draw_manifest_init(manifest);
    manifest->coverageBoundsValid = 1;
    VSET(manifest->coverageBoundsMin, -100.0, -200.0, -300.0);
    VSET(manifest->coverageBoundsMax, 400.0, 500.0, 600.0);
    manifest->occurrenceCount = 2;
    manifest->occurrences = static_cast<BObolDrawManifestOccurrence *>(
	bu_calloc(manifest->occurrenceCount, sizeof(*manifest->occurrences),
	    "test draw manifest occurrences"));
    if (!manifest->occurrences)
	return 0;

    BObolDrawManifestOccurrence &leaf = manifest->occurrences[0];
    leaf.path = bu_strdup(path_full_name);
    leaf.sourceName = bu_strdup(path_leaf_name);
    leaf.sourceMeshRequestValid = 1;
    leaf.meshAssetKind = BOBOL_DRAW_CACHE_MESH_ASSET_BREP;
    leaf.meshAssetContentHash = 0x123456789abcdef0ULL;
    leaf.meshAssetTessellationAbsTol = 0.125;
    leaf.meshAssetTessellationRelTol = 0.01;
    leaf.meshAssetTessellationNormTol = 0.5;
    leaf.meshAssetPath = bu_strdup(path_full_name);
    leaf.meshAssetName = bu_strdup(path_leaf_name);
    VSET(leaf.meshAssetBoundsMin, -10.0, -20.0, -30.0);
    VSET(leaf.meshAssetBoundsMax, 40.0, 50.0, 60.0);
    MAT_IDN(leaf.meshAssetMatrix);
    MAT_DELTAS(leaf.meshAssetMatrix, 2.0, 4.0, 8.0);
    leaf.sourceFaceCount = 123456;
    leaf.sourcePointCount = 65432;
    MAT_IDN(leaf.localMatrix);
    MAT_DELTAS(leaf.localMatrix, 11.0, 12.0, 13.0);
    VSET(leaf.boundsMin, -1.0, -2.0, -3.0);
    VSET(leaf.boundsMax, 4.0, 5.0, 6.0);
    leaf.booleanOperation = DB_OP_UNION;
    leaf.occurrenceIndex = 3;
    leaf.metadataValid = 1;
    bobol_draw_metadata_record_init(&leaf.metadata);
    leaf.metadata.directoryFound = 1;
    leaf.metadata.hasColor = 1;
    leaf.metadata.color[0] = 10;
    leaf.metadata.color[1] = 20;
    leaf.metadata.color[2] = 30;

    BObolDrawManifestOccurrence &region = manifest->occurrences[1];
    region.path = bu_strdup(path_region_full_name);
    region.sourceName = bu_strdup(path_region_name);
    MAT_IDN(region.localMatrix);
    MAT_DELTAS(region.localMatrix, 11.0, 13.0, 17.0);
    VSET(region.boundsMin, -4.0, -5.0, -6.0);
    VSET(region.boundsMax, 7.0, 8.0, 9.0);
    region.booleanOperation = DB_OP_SUBTRACT;
    region.occurrenceIndex = 7;

    if (!leaf.path || !leaf.sourceName || !leaf.meshAssetPath ||
	!leaf.meshAssetName || !region.path || !region.sourceName) {
	bobol_draw_manifest_free(manifest);
	return 0;
    }
    return 1;
}

static int
make_path_metadata_tree(rt_wdb *wdbp)
{
    point_t path_center;
    point_t table_center;
    unsigned char region_rgb[3] = {10, 20, 30};
    unsigned char top_rgb[3] = {120, 130, 140};
    struct wmember region_members;
    struct wmember top_members;
    struct wmember table_region_members;
    struct wmember table_top_members;

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

    VSET(table_center, 40.0, 0.0, 0.0);
    if (mk_sph(wdbp, path_table_leaf_name, table_center, 2.0) < 0)
	return 1;

    BU_LIST_INIT(&table_region_members.l);
    if (!mk_addmember(path_table_leaf_name, &table_region_members.l, NULL,
	    WMOP_UNION))
	return 1;
    if (mk_comb(wdbp, path_table_region_name, &table_region_members.l, 1,
		NULL, NULL, NULL, 177, 0, 0, 0, 0, 0, 0) < 0)
	return 1;

    BU_LIST_INIT(&table_top_members.l);
    if (!mk_addmember(path_table_region_name, &table_top_members.l, NULL,
	    WMOP_UNION))
	return 1;
    if (mk_comb(wdbp, path_table_top_name, &table_top_members.l, 0,
		NULL, NULL, NULL, 0, 0, 0, 0, 0, 0, 0) < 0)
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
check_plain_leaf_metadata(const BObolDrawMetadataRecord *record)
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
check_metadata(const BObolDrawMetadataRecord *record)
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
	bu_strcmp(record->materialName, "DrawCacheMaterial") != 0 ||
	!record->hasShader ||
	bu_strcmp(record->shader, "plastic {di .7 sp .2}") != 0) {
	printf("FAIL: draw metadata fields\n");
	return 1;
    }

    return 0;
}

static int
check_path_metadata(const BObolDrawMetadataRecord *record)
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

static int
check_path_color_table_metadata(
    db_i *dbip,
    const BObolDrawMetadataRecord *record)
{
    struct region rp;
    unsigned char expected[3] = {0, 0, 0};

    memset(&rp, 0, sizeof(rp));
    rp.reg_regionid = 177;
    db_mater_color_region(dbip, &rp);
    if (!rp.reg_mater.ma_color_valid) {
	printf("FAIL: path draw metadata color table setup\n");
	return 1;
    }
    expected[0] = static_cast<unsigned char>(
		     std::lround(rp.reg_mater.ma_color[0] * 255.0));
    expected[1] = static_cast<unsigned char>(
		     std::lround(rp.reg_mater.ma_color[1] * 255.0));
    expected[2] = static_cast<unsigned char>(
		     std::lround(rp.reg_mater.ma_color[2] * 255.0));

    if (!record || !record->directoryFound || !record->isSolid ||
	record->isComb || !record->isRegion ||
	!record->hasRegionId || record->regionId != 177 ||
	!record->hasColor || record->color[0] != expected[0] ||
	record->color[1] != expected[1] ||
	record->color[2] != expected[2]) {
	printf("FAIL: path draw metadata color table fallback\n");
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
    BObolDrawCacheStatus status;
    BObolDrawProxyRecord proxy;
    BObolDrawLodAssetRecord lodAsset;
    BObolDrawLodAssetRecord loadedLodAsset;
    BObolDrawMetadataRecord metadata;
    BObolDrawManifest manifest;
    BObolDrawManifest loadedManifest;
    BObolDrawManifest manifestDescription;
    point_t center;
    point_t aabbMin;
    point_t aabbMax;
    point_t aabbPoints[2];
    point_t obbPoints[8];

    bu_setprogname(argv[0]);
    bobol_draw_manifest_init(&manifest);
    bobol_draw_manifest_init(&loadedManifest);
    bobol_draw_manifest_init(&manifestDescription);

    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    bu_dir(cacheDir, MAXPATHLEN, BU_DIR_CURR,
	   "bobol_draw_cache_test", NULL);
    bu_mkdir(cacheDir);
    bu_setenv("BU_DIR_CACHE", cacheDir, 1);

    VSET(center, 1.0, 2.0, 3.0);
    VSET(aabbMin, -3.0, -2.0, -1.0);
    VSET(aabbMax, 5.0, 6.0, 7.0);
	VMOVE(aabbPoints[0], aabbMin);
	VMOVE(aabbPoints[1], aabbMax);
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
    db_mater_add(dbip, 177, 177, 41, 42, 43, MATER_NO_ADDR);

    {
	rt_wdb *wdbp = wdb_dbopen(dbip, RT_WDB_TYPE_DB_DISK);
	fastf_t assetVertices[12] = {
	    0.0, 0.0, 0.0,
	    4.0, 0.0, 0.0,
	    0.0, 2.0, 0.0,
	    0.0, 0.0, 1.0
	};
	fastf_t copyVertices[12] = {
	    10.0, -3.0, 2.0,
	    14.0, -3.0, 2.0,
	    10.0, -1.0, 2.0,
	    10.0, -3.0, 3.0
	};
	int botFaces[12] = {
	    0, 1, 2,
	    0, 3, 1,
	    1, 3, 2,
	    2, 3, 0
	};
	if (!wdbp || mk_sph(wdbp, objname, center, 4.0) < 0 ||
	    mk_bot(wdbp, lod_asset_name, RT_BOT_SOLID, RT_BOT_CCW, 0,
		4, 4, assetVertices, botFaces, NULL, NULL) < 0 ||
	    mk_bot(wdbp, lod_copy_name, RT_BOT_SOLID, RT_BOT_CCW, 0,
		4, 4, copyVertices, botFaces, NULL, NULL) < 0 ||
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

    bobol_draw_lod_asset_record_init(&lodAsset);
    bu_strlcpy(lodAsset.assetName, lod_asset_name,
	sizeof(lodAsset.assetName));
    lodAsset.faceCount = 4;
    lodAsset.pointCount = 4;
    VSET(lodAsset.boundsMin, 10.0, -3.0, 2.0);
    VSET(lodAsset.boundsMax, 14.0, -1.0, 3.0);
    VSET(lodAsset.assetBoundsMin, 0.0, 0.0, 0.0);
    VSET(lodAsset.assetBoundsMax, 4.0, 2.0, 1.0);
    MAT_DELTAS(lodAsset.assetToObject, 10.0, -3.0, 2.0);
    if (bobol_draw_lod_asset_cache_store(dbip, lod_copy_name,
	    &lodAsset) != BRLCAD_OK ||
	bobol_draw_lod_asset_cache_get(dbip, lod_copy_name,
	    &loadedLodAsset) != BRLCAD_OK ||
	bu_strcmp(loadedLodAsset.assetName, lod_asset_name) != 0 ||
	loadedLodAsset.faceCount != 4 || loadedLodAsset.pointCount != 4 ||
	!point_equal(loadedLodAsset.boundsMin, lodAsset.boundsMin) ||
	!point_equal(loadedLodAsset.boundsMax, lodAsset.boundsMax) ||
	!point_equal(loadedLodAsset.assetBoundsMin,
	    lodAsset.assetBoundsMin) ||
	!point_equal(loadedLodAsset.assetBoundsMax,
	    lodAsset.assetBoundsMax) ||
	!fastf_equal(loadedLodAsset.assetToObject[MDX], 10.0) ||
	!fastf_equal(loadedLodAsset.assetToObject[MDY], -3.0) ||
	!fastf_equal(loadedLodAsset.assetToObject[MDZ], 2.0)) {
	printf("FAIL: transformed LoD asset cache store/get\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_status(dbip, objname,
					BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || status.hasCachedPayload) {
	printf("FAIL: draw cache initial AABB status\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_refresh(dbip, objname,
					 BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry) {
	printf("FAIL: draw cache AABB refresh\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_AABB,
				     &proxy) != BRLCAD_OK || proxy.kind != BOBOL_LOD_PROXY_AABB ||
	proxy.pointCount != 2 ||
	check_proxy_point("AABB min", &proxy, 0, aabbMin) ||
	check_proxy_point("AABB max", &proxy, 1, aabbMax)) {
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_store(dbip, objname, BOBOL_LOD_PROXY_OBB,
				       obbPoints, 8, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry) {
	printf("FAIL: draw cache OBB store\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK || proxy.kind != BOBOL_LOD_PROXY_OBB ||
	proxy.pointCount != 8 ||
	check_proxy_point("OBB first", &proxy, 0, obbPoints[0]) ||
	check_proxy_point("OBB last", &proxy, 7, obbPoints[7])) {
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_status(dbip, path_top_name,
					BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || status.hasCachedPayload ||
	bobol_draw_proxy_cache_refresh(dbip, path_top_name,
					 BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_ERROR ||
	!status.directoryFound || status.hasCachedPayload ||
	bobol_draw_proxy_cache_store(dbip, path_top_name,
				       BOBOL_LOD_PROXY_AABB, aabbPoints, 2,
				       &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	bobol_draw_proxy_cache_get(dbip, path_top_name,
				     BOBOL_LOD_PROXY_AABB, &proxy) != BRLCAD_OK ||
	proxy.pointCount != 2 ||
	check_proxy_point("comb AABB min", &proxy, 0, aabbMin) ||
	check_proxy_point("comb AABB max", &proxy, 1, aabbMax)) {
	printf("FAIL: draw cache derived combination proxy\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_metadata_cache_refresh(dbip, objname, &status) !=
	BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry ||
	bobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_OK ||
	check_metadata(&metadata)) {
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_metadata_cache_refresh(dbip, path_leaf_name, &status) !=
	BRLCAD_OK ||
	bobol_draw_metadata_cache_get(dbip, path_leaf_name, &metadata) !=
	BRLCAD_OK ||
	check_plain_leaf_metadata(&metadata)) {
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_status(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || status.hasCachedPayload) {
	printf("FAIL: draw cache initial path metadata status\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	!status.generatedCacheEntry ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK ||
	check_path_metadata(&metadata)) {
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_refresh(dbip, path_table_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || !status.hasCachedPayload ||
	bobol_draw_path_metadata_cache_get(dbip, path_table_full_name,
		&metadata) != BRLCAD_OK ||
	check_path_color_table_metadata(dbip, &metadata)) {
	ret = 1;
	goto cleanup;
    }

    if (!make_manifest(&manifest) ||
	bobol_draw_manifest_cache_store(dbip, path_top_name, &manifest) !=
	BRLCAD_OK ||
	bobol_draw_manifest_cache_describe(dbip, path_top_name,
	    &manifestDescription) != BRLCAD_OK ||
	!manifestDescription.coverageBoundsValid ||
	manifestDescription.occurrenceCount != manifest.occurrenceCount ||
	manifestDescription.occurrences != NULL ||
	!point_equal(manifestDescription.coverageBoundsMin,
	    manifest.coverageBoundsMin) ||
	!point_equal(manifestDescription.coverageBoundsMax,
	    manifest.coverageBoundsMax) ||
	bobol_draw_manifest_cache_get(dbip, path_top_name, &loadedManifest) !=
	BRLCAD_OK || check_manifest(&loadedManifest)) {
	printf("FAIL: draw manifest store/describe/get\n");
	ret = 1;
	goto cleanup;
    }
    bobol_draw_manifest_free(&manifest);
    bobol_draw_manifest_free(&loadedManifest);

    db_close(dbip);
    dbip = db_open(dbpath, DB_OPEN_READWRITE);
    if (!dbip || db_dirbuild(dbip) < 0) {
	printf("FAIL: draw cache reopen\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK ||
	check_proxy_point("OBB reopen", &proxy, 7, obbPoints[7]) ||
	bobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_OK ||
	check_metadata(&metadata)) {
	printf("FAIL: draw cache persistent context reopen data\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_get(dbip, path_full_name,
	    &metadata) != BRLCAD_OK ||
	check_path_metadata(&metadata)) {
	printf("FAIL: draw cache persistent path metadata reopen data\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_manifest_cache_get(dbip, path_top_name,
	&loadedManifest) != BRLCAD_OK || check_manifest(&loadedManifest)) {
	printf("FAIL: draw manifest persistent reopen data\n");
	ret = 1;
	goto cleanup;
    }
    bobol_draw_manifest_free(&loadedManifest);
    if (bobol_draw_lod_asset_cache_get(dbip, lod_copy_name,
	    &loadedLodAsset) != BRLCAD_OK ||
	bu_strcmp(loadedLodAsset.assetName, lod_asset_name) != 0 ||
	!fastf_equal(loadedLodAsset.assetToObject[MDX], 10.0)) {
	printf("FAIL: transformed LoD asset cache persistent reopen data\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_invalidate(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
		&status) != BRLCAD_OK ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK ||
	check_path_metadata(&metadata)) {
	printf("FAIL: draw cache path metadata invalidate/refresh\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_refresh(dbip, path_region_full_name,
	    &status) != BRLCAD_OK ||
	bobol_draw_path_metadata_cache_get(dbip, path_region_full_name,
		&metadata) != BRLCAD_OK ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK) {
	printf("FAIL: draw cache path metadata object-invalidate setup\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_invalidate_object(dbip,
	    path_leaf_name, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_get(dbip, path_region_full_name,
		&metadata) != BRLCAD_OK) {
	printf("FAIL: draw cache path metadata leaf object invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	bobol_draw_path_metadata_cache_invalidate_object(dbip,
		path_region_name, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_get(dbip, path_region_full_name,
		&metadata) != BRLCAD_ERROR ||
	bobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_OK) {
	printf("FAIL: draw cache path metadata ancestor object invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_path_metadata_cache_refresh(dbip, path_full_name,
	    &status) != BRLCAD_OK ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) != BRLCAD_OK) {
	printf("FAIL: draw cache path metadata object-invalidate reseed\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_proxy_cache_invalidate(dbip, objname,
					    BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_OK ||
	!status.directoryFound || !status.clearedCacheEntry ||
	status.hasCachedPayload ||
	bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_AABB,
				     &proxy) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK) {
	printf("FAIL: draw cache targeted invalidate\n");
	ret = 1;
	goto cleanup;
    }

    if (bobol_draw_manifest_cache_invalidate_database(dbip) != BRLCAD_OK ||
	bobol_draw_manifest_cache_get(dbip, path_top_name,
		&loadedManifest) != BRLCAD_ERROR ||
	bobol_draw_lod_asset_cache_get(dbip, lod_copy_name,
		&loadedLodAsset) != BRLCAD_ERROR ||
	!make_manifest(&manifest) ||
	bobol_draw_manifest_cache_store(dbip, path_top_name,
		&manifest) != BRLCAD_OK) {
	printf("FAIL: draw manifest database invalidation\n");
	ret = 1;
	goto cleanup;
    }
    bobol_draw_manifest_free(&manifest);

    if (bobol_draw_cache_clear_database(dbip) != BRLCAD_OK ||
	bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_ERROR ||
	bobol_draw_metadata_cache_get(dbip, objname, &metadata) !=
	BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_get(dbip, path_full_name,
		&metadata) !=
	BRLCAD_ERROR ||
	bobol_draw_manifest_cache_get(dbip, path_top_name,
		&loadedManifest) != BRLCAD_ERROR) {
	printf("FAIL: draw cache database clear\n");
	ret = 1;
	goto cleanup;
    }

    if (!make_manifest(&manifest)) {
	printf("FAIL: draw manifest invalid-operation setup\n");
	ret = 1;
	goto cleanup;
    }
    manifest.occurrences[0].booleanOperation = 0;
    if (bobol_draw_manifest_cache_store(dbip, path_top_name, &manifest) !=
	BRLCAD_ERROR) {
	printf("FAIL: draw manifest accepted non-database boolean operation\n");
	ret = 1;
	goto cleanup;
    }
    bobol_draw_manifest_free(&manifest);

    if (!make_manifest(&manifest)) {
	printf("FAIL: draw manifest invalid-coverage setup\n");
	ret = 1;
	goto cleanup;
    }
    manifest.coverageBoundsMin[X] =
	manifest.coverageBoundsMax[X] + 1.0;
    if (bobol_draw_manifest_cache_store(dbip, path_top_name, &manifest) !=
	BRLCAD_ERROR) {
	printf("FAIL: draw manifest accepted invalid coverage bounds\n");
	ret = 1;
	goto cleanup;
    }
    bobol_draw_manifest_free(&manifest);

    if (bobol_draw_proxy_cache_store(dbip, objname, BOBOL_LOD_PROXY_OBB,
				       obbPoints, 8, &status) != BRLCAD_OK ||
	bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK) {
	printf("FAIL: draw cache reuse after database clear\n");
	ret = 1;
	goto cleanup;
    }

    bobol_draw_cache_clear_all();
    if (bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_store(dbip, objname, BOBOL_LOD_PROXY_OBB,
				       obbPoints, 8, &status) != BRLCAD_OK ||
	bobol_draw_proxy_cache_get(dbip, objname, BOBOL_LOD_PROXY_OBB,
				     &proxy) != BRLCAD_OK) {
	printf("FAIL: draw cache clear-all context reset\n");
	ret = 1;
	goto cleanup;
    }

    bobol_draw_cache_status_init(&status);
    status.generatedCacheEntry = 1;
    bobol_draw_cache_status_init(&status);
    bobol_draw_proxy_record_init(&proxy);
    bobol_draw_lod_asset_record_init(&loadedLodAsset);
    loadedLodAsset.faceCount = 1;
    bobol_draw_lod_asset_record_init(&loadedLodAsset);
    proxy.kind = 9;
    bobol_draw_proxy_record_init(&proxy);
    bobol_draw_metadata_record_init(&metadata);
    metadata.directoryFound = 1;
    bobol_draw_metadata_record_init(&metadata);
    if (status.generatedCacheEntry || proxy.kind ||
	loadedLodAsset.faceCount || metadata.directoryFound ||
	bobol_draw_lod_asset_cache_get(NULL, lod_copy_name,
	    &loadedLodAsset) != BRLCAD_ERROR ||
	bobol_draw_lod_asset_cache_get(dbip, NULL,
	    &loadedLodAsset) != BRLCAD_ERROR ||
	bobol_draw_lod_asset_cache_get(dbip, lod_copy_name,
	    NULL) != BRLCAD_ERROR ||
	bobol_draw_lod_asset_cache_store(NULL, lod_copy_name,
	    &lodAsset) != BRLCAD_ERROR ||
	bobol_draw_lod_asset_cache_store(dbip, NULL,
	    &lodAsset) != BRLCAD_ERROR ||
	bobol_draw_lod_asset_cache_store(dbip, lod_copy_name,
	    NULL) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_status(NULL, objname,
					BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_status(dbip, NULL,
					BOBOL_LOD_PROXY_AABB, &status) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_status(dbip, objname, 999,
					&status) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_status(dbip, objname,
					BOBOL_LOD_PROXY_AABB, NULL) != BRLCAD_ERROR ||
	bobol_draw_proxy_cache_store(dbip, objname,
				       BOBOL_LOD_PROXY_AABB, NULL, 2, &status) != BRLCAD_ERROR ||
	bobol_draw_metadata_cache_status(NULL, objname,
					   &status) != BRLCAD_ERROR ||
	bobol_draw_metadata_cache_status(dbip, NULL,
					   &status) != BRLCAD_ERROR ||
	bobol_draw_metadata_cache_status(dbip, objname,
					   NULL) != BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_status(NULL, path_full_name,
		&status) != BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_status(dbip, NULL,
		&status) != BRLCAD_ERROR ||
	bobol_draw_path_metadata_cache_status(dbip, path_full_name,
		NULL) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_get(NULL, path_top_name,
		&loadedManifest) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_get(dbip, NULL,
		&loadedManifest) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_get(dbip, path_top_name,
		NULL) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_describe(NULL, path_top_name,
		&manifestDescription) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_describe(dbip, NULL,
		&manifestDescription) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_describe(dbip, path_top_name,
		NULL) != BRLCAD_ERROR ||
	bobol_draw_manifest_cache_invalidate_database(NULL) != BRLCAD_ERROR) {
	printf("FAIL: draw cache null handling\n");
	ret = 1;
	goto cleanup;
    }

cleanup:
    bobol_draw_manifest_free(&manifest);
    bobol_draw_manifest_free(&loadedManifest);
    bobol_draw_cache_clear_all();
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
