/*             T E S T _ L O D _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/str.h"

#include "brlobol/lod_realization.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbVec3f.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

static BRLObolLodRequest
make_request(void)
{
    BRLObolLodRequest request;

    request.databaseId = "db://unit-test";
    request.databaseRevision = 3;
    request.sourceRevision = 7;
    request.sourceContentHash = 0x1234;
    request.objectPath = "/all.g/mesh.bot";
    request.objectName = "mesh.bot";
    request.viewRevision = 11;
    request.policyRevision = 13;
    request.drawMode = BRLOBOL_LOD_DRAW_SHADED;
    request.lodPolicy = 5;
    request.providerId = "unit-provider";
    request.providerVersion = "1";
    request.qualityTier = BRLOBOL_LOD_QUALITY_FAST_DISPLAY;
    request.bounds = SbBox3f(SbVec3f(-1.0f, -2.0f, -3.0f),
			     SbVec3f(4.0f, 5.0f, 6.0f));
    request.sourceCounts.faceCount = 17;
    request.sourceCounts.pointCount = 19;
    request.sourceCounts.originalPointCount = 23;
    request.sourceCounts.normalCount = 29;
    request.addProviderParam("ratio", "0.66");
    request.addProviderParam("provider-option", "true");

    return request;
}

static int
test_key_determinism(void)
{
    BRLObolLodRequest a = make_request();
    BRLObolLodRequest b = make_request();

    b.providerParams.clear();
    b.addProviderParam("provider-option", "true");
    b.addProviderParam("ratio", "0.66");

    BRLObolLodCacheKey key_a = brlobol_lod_cache_key(a);
    BRLObolLodCacheKey key_b = brlobol_lod_cache_key(b);
    if (!key_a.isValid() || !key_b.isValid() ||
	bu_strcmp(key_a.value.getString(), key_b.value.getString()) != 0) {
	printf("FAIL: LoD cache key not deterministic across parameter order\n");
	return 1;
    }

    b.viewRevision++;
    BRLObolLodCacheKey key_c = brlobol_lod_cache_key(b);
    if (bu_strcmp(key_a.value.getString(), key_c.value.getString()) == 0) {
	printf("FAIL: LoD cache key did not include view revision\n");
	return 1;
    }

    BRLObolLodResult result;
    result.request = a;
    result.cacheKey = key_a;
    if (!brlobol_lod_result_matches_request(result, a)) {
	printf("FAIL: LoD result did not match identical request\n");
	return 1;
    }
    if (brlobol_lod_result_matches_request(result, b)) {
	printf("FAIL: LoD result matched stale request revision\n");
	return 1;
    }

    return 0;
}

static int
test_rt_mesh_result(void)
{
    BRLObolLodRequest request = make_request();
    struct BRLObolMeshLodInfo info = BRLOBOL_MESH_LOD_INFO_INIT;
    struct BRLObolMeshLodCacheStatus status = BRLOBOL_MESH_LOD_CACHE_STATUS_INIT;

    info.active_level = 2;
    info.face_count = 8;
    info.point_count = 9;
    info.point_orig_count = 10;
    info.normal_count = 24;
    info.has_faces = 1;
    info.has_points = 1;
    info.has_original_points = 1;
    info.has_snapped_points = 1;
    info.has_normals = 1;
    VSET(info.bmin, -1.0, -2.0, -3.0);
    VSET(info.bmax, 4.0, 5.0, 6.0);

    status.directory_found = 1;
    status.is_bot = 1;
    status.has_cache_key = 1;
    status.has_cached_payload = 1;
    status.cache_key = 0xfeed;

    BRLObolLodResult result =
	brlobol_lod_result_from_mesh_lod_info(request, info, &status);
    if (result.resultKind != BRLOBOL_LOD_RESULT_MESH ||
	result.qualityTier != BRLOBOL_LOD_QUALITY_FAST_DISPLAY ||
	result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	result.geometry.kind != BRLOBOL_LOD_GEOMETRY_MESH_LOD_CACHE ||
	result.geometry.activeLevel != 2 ||
	result.geometry.providerToken != 0xfeed ||
	!result.geometry.isValid() ||
	result.counts.faceCount != 8 ||
	result.counts.pointCount != 9 ||
	result.counts.originalPointCount != 10 ||
	result.counts.normalCount != 24 ||
	!result.hasSnappedPoints ||
	!result.hasNormals ||
	!brlobol_lod_result_matches_request(result, request)) {
	printf("FAIL: Obol mesh LoD result conversion\n");
	return 1;
    }

    status.stale_cache_entry = 1;
    result = brlobol_lod_result_from_mesh_lod_info(request, info, &status);
    if (!result.stale ||
	result.providerStatus != BRLOBOL_LOD_PROVIDER_STALE ||
	result.diagnostic.getLength() == 0) {
	printf("FAIL: Obol mesh LoD stale status conversion\n");
	return 1;
    }

    return 0;
}

static int
test_rt_mesh_payload_copy(void)
{
    point_t points[4];
    vect_t normals[6];
    int faces[6] = {0, 1, 2, 0, 3, 1};
    struct BRLObolMeshLodData data;
    BRLObolLodMeshPayload payload;

    memset(&data, 0, sizeof(data));
    VSET(points[0], 0.0, 0.0, 0.0);
    VSET(points[1], 2.0, 0.0, 0.0);
    VSET(points[2], 0.0, 2.0, 0.0);
    VSET(points[3], 0.0, 0.0, 2.0);
    for (size_t i = 0; i < 6; i++)
	VSET(normals[i], 0.0, 0.0, 1.0 + (fastf_t)i);
    data.faces = faces;
    data.face_count = 2;
    data.points = points;
    data.point_count = 4;
    data.normals = normals;
    data.normal_count = 6;

    if (!brlobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	!payload.isValid() ||
	payload.points.size() != 4 ||
	payload.normals.size() != 6 ||
	payload.coordIndex.size() != 6 ||
	payload.coordIndex[4] != 3 ||
	fabs((double)payload.normals[5][2] - 6.0) > 1.0e-6 ||
	fabs((double)payload.points[1][0] - 2.0) > 1.0e-6) {
	printf("FAIL: Obol mesh LoD payload copy\n");
	return 1;
    }

    faces[5] = 9;
    if (brlobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	payload.isValid() ||
	!payload.points.empty() ||
	!payload.coordIndex.empty()) {
	printf("FAIL: Obol mesh LoD payload copy accepted invalid index\n");
	return 1;
    }

    faces[5] = 1;
    data.normal_count = 5;
    if (brlobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	payload.isValid() ||
	!payload.points.empty() ||
	!payload.normals.empty() ||
	!payload.coordIndex.empty()) {
	printf("FAIL: Obol mesh LoD payload copy accepted invalid normals\n");
	return 1;
    }

    data.normals = NULL;
    data.normal_count = 0;
    if (!brlobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	!payload.isValid() ||
	!payload.normals.empty()) {
	printf("FAIL: Obol mesh LoD payload copy rejected mesh without normals\n");
	return 1;
    }

    return 0;
}

static int
test_stage_results(void)
{
    BRLObolLodRequest request = make_request();
    std::vector<BRLObolLodDependency> dependencies;
    std::vector<BRLObolLodAttribute> attributes;
    BRLObolLodCounts counts;
    BRLObolLodProxy proxy;

    BRLObolLodDependency dependency;
    dependency.objectPath = "/all.g/child.bot";
    dependency.objectName = "child.bot";
    dependency.sourceRevision = 31;
    dependency.sourceContentHash = 0x77;
    dependency.requiredQualityTier = BRLOBOL_LOD_QUALITY_PROXY;
    dependency.optional = FALSE;
    dependencies.push_back(dependency);

    BRLObolLodResult result =
	brlobol_lod_directory_result(request, dependencies);
    if (result.resultKind != BRLOBOL_LOD_RESULT_DIRECTORY ||
	result.qualityTier != BRLOBOL_LOD_QUALITY_METADATA ||
	result.providerStatus != BRLOBOL_LOD_PROVIDER_READY ||
	result.terminal ||
	result.dependencies.size() != 1 ||
	bu_strcmp(result.dependencies[0].objectPath.getString(),
	       "/all.g/child.bot") != 0 ||
	!brlobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD directory stage result\n");
	return 1;
    }

    BRLObolLodAttribute color;
    color.name = "display.color";
    color.value = "255 0 0";
    attributes.push_back(color);
    result = brlobol_lod_attributes_result(request, attributes);
    if (result.resultKind != BRLOBOL_LOD_RESULT_ATTRIBUTES ||
	result.qualityTier != BRLOBOL_LOD_QUALITY_ATTRIBUTES ||
	result.attributes.size() != 1 ||
	bu_strcmp(result.attributes[0].name.getString(), "display.color") != 0 ||
	!brlobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD attributes stage result\n");
	return 1;
    }

    counts.lineCount = 12;
    result = brlobol_lod_aabb_result(request, request.bounds, &counts);
    if (result.resultKind != BRLOBOL_LOD_RESULT_AABB ||
	result.qualityTier != BRLOBOL_LOD_QUALITY_PROXY ||
	result.proxy.kind != BRLOBOL_LOD_PROXY_AABB ||
	!result.proxy.isValid() ||
	result.bounds.isEmpty() ||
	result.counts.lineCount != 12 ||
	!brlobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD AABB stage result\n");
	return 1;
    }

    proxy.kind = BRLOBOL_LOD_PROXY_OBB;
    proxy.bounds = request.bounds;
    proxy.center.setValue(1.0f, 2.0f, 3.0f);
    proxy.axisX.setValue(1.0f, 0.0f, 0.0f);
    proxy.axisY.setValue(0.0f, 1.0f, 0.0f);
    proxy.axisZ.setValue(0.0f, 0.0f, 1.0f);
    proxy.halfExtents.setValue(4.0f, 5.0f, 6.0f);
    result = brlobol_lod_proxy_result(request, proxy, NULL);
    if (result.resultKind != BRLOBOL_LOD_RESULT_PROXY ||
	result.qualityTier != BRLOBOL_LOD_QUALITY_PROXY ||
	result.proxy.kind != BRLOBOL_LOD_PROXY_OBB ||
	!result.proxy.isValid() ||
	fabs((double)result.proxy.center[0] - 1.0) > 1.0e-6 ||
	fabs((double)result.proxy.halfExtents[2] - 6.0) > 1.0e-6 ||
	!brlobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD proxy stage result\n");
	return 1;
    }

    proxy.clear();
    result = brlobol_lod_proxy_result(request, proxy, NULL);
    if (result.providerStatus != BRLOBOL_LOD_PROVIDER_ERROR ||
	!result.terminal ||
	result.diagnostic.getLength() == 0) {
	printf("FAIL: LoD invalid proxy diagnostic result\n");
	return 1;
    }

    result.clear();
    result.addDependency("/all.g/leaf.bot", "leaf.bot", 41, 0x99,
			 BRLOBOL_LOD_QUALITY_FAST_DISPLAY, TRUE);
    result.addAttribute("draw.mode", "wire");
    if (result.dependencies.size() != 1 ||
	!result.dependencies[0].optional ||
	result.attributes.size() != 1 ||
	bu_strcmp(result.attributes[0].value.getString(), "wire") != 0) {
	printf("FAIL: LoD result append helpers\n");
	return 1;
    }

    return 0;
}

int
main(int argc, char **argv)
{
    bu_setprogname(argv[0]);
    if (argc != 1) {
	printf("Usage: %s\n", argv[0]);
	return 1;
    }

    if (test_key_determinism())
	return 1;
    if (test_rt_mesh_result())
	return 1;
    if (test_rt_mesh_payload_copy())
	return 1;
    if (test_stage_results())
	return 1;

    return 0;
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
