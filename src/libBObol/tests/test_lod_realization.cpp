/*             T E S T _ L O D _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/app.h"

#include "bu/str.h"

#include "BObol/BLodRealization.h"

#include <Obol/cad/SoCADAssembly.h>

#include <Inventor/SbBox.h>
#include <Inventor/SbVec3f.h>

#include <math.h>
#include <stdio.h>
#include <string.h>
#include <vector>

static void
complete_test_hierarchy(struct BObolMeshLodHierarchyInfo &hierarchy)
{
    hierarchy.cut_count = static_cast<uint32_t>(hierarchy.max_cut + 1);
    for (int cut = 0; cut <= hierarchy.max_cut; ++cut) {
	const bool exact = cut == hierarchy.max_cut;
	hierarchy.cuts[cut].object_error =
	    static_cast<double>(hierarchy.max_cut - cut);
	hierarchy.cuts[cut].quantization_bits[X] = exact ? 16 : 15;
	hierarchy.cuts[cut].quantization_bits[Y] = exact ? 16 : 15;
	hierarchy.cuts[cut].quantization_bits[Z] = exact ? 16 : 15;
	hierarchy.cuts[cut].exact = exact ? 1 : 0;
    }
}

static BObolLodRequest
make_request(void)
{
    BObolLodRequest request;

    request.databaseId = "db://unit-test";
    request.databaseRevision = 3;
    request.sourceRevision = 7;
    request.sourceContentHash = 0x1234;
    request.objectPath = "/all.g/mesh.bot";
    request.objectName = "mesh.bot";
    request.viewRevision = 11;
    request.policyRevision = 13;
    request.drawMode = BOBOL_LOD_DRAW_SHADED;
    request.lodPolicy = 5;
    request.providerId = "unit-provider";
    request.providerVersion = "1";
    request.qualityTier = BOBOL_LOD_QUALITY_FAST_DISPLAY;
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
    BObolLodRequest a = make_request();
    BObolLodRequest b = make_request();

    b.providerParams.clear();
    b.addProviderParam("provider-option", "true");
    b.addProviderParam("ratio", "0.66");

    BObolLodCacheKey key_a = bobol_lod_cache_key(a);
    BObolLodCacheKey key_b = bobol_lod_cache_key(b);
    if (!key_a.isValid() || !key_b.isValid() ||
	bu_strcmp(key_a.value.getString(), key_b.value.getString()) != 0) {
	printf("FAIL: LoD cache key not deterministic across parameter order\n");
	return 1;
    }
    if (!bobol_lod_request_keys_equal(a, b)) {
	printf("FAIL: structured LoD key not parameter-order independent\n");
	return 1;
    }

    b.viewRevision++;
    BObolLodCacheKey key_c = bobol_lod_cache_key(b);
    if (bu_strcmp(key_a.value.getString(), key_c.value.getString()) == 0) {
	printf("FAIL: LoD cache key did not include view revision\n");
	return 1;
    }
    if (bobol_lod_request_keys_equal(a, b)) {
	printf("FAIL: structured LoD key ignored view revision\n");
	return 1;
    }

    BObolLodResult result;
    result.request = a;
    result.cacheKey = key_a;
    if (!bobol_lod_result_matches_request(result, a)) {
	printf("FAIL: LoD result did not match identical request\n");
	return 1;
    }
    if (bobol_lod_result_matches_request(result, b)) {
	printf("FAIL: LoD result matched stale request revision\n");
	return 1;
    }

    return 0;
}

static int
test_rt_mesh_result(void)
{
    BObolLodRequest request = make_request();
    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    struct BObolMeshLodCacheStatus status = BOBOL_MESH_LOD_CACHE_STATUS_INIT;

    info.active_cut = 2;
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

    BObolLodResult result =
	bobol_lod_result_from_mesh_lod_info(request, info, &status);
    if (result.resultKind != BOBOL_LOD_RESULT_MESH ||
	result.payloadKind != BOBOL_LOD_PAYLOAD_MESH ||
	!result.payloadIsConsistent() ||
	result.qualityTier != BOBOL_LOD_QUALITY_FAST_DISPLAY ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	result.geometry.kind != BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE ||
	result.geometry.activeCut != 2 ||
	result.resolvedCut != 2 ||
	result.geometry.providerToken != 0xfeed ||
	!result.geometry.isValid() ||
	result.counts.faceCount != 8 ||
	result.counts.pointCount != 9 ||
	result.counts.originalPointCount != 10 ||
	result.counts.normalCount != 24 ||
	!result.hasSnappedPoints ||
	!result.hasNormals ||
	!bobol_lod_result_matches_request(result, request)) {
	printf("FAIL: Obol mesh LoD result conversion\n");
	return 1;
    }

    status.stale_cache_entry = 1;
    result = bobol_lod_result_from_mesh_lod_info(request, info, &status);
    if (!result.stale ||
	result.providerStatus != BOBOL_LOD_PROVIDER_STALE ||
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
    uint32_t faces[6] = {0, 1, 2, 0, 3, 1};
    struct BObolMeshLodData data;
    BObolLodMeshPayload payload;

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

    if (!bobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
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
    if (bobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	payload.isValid() ||
	!payload.points.empty() ||
	!payload.coordIndex.empty()) {
	printf("FAIL: Obol mesh LoD payload copy accepted invalid index\n");
	return 1;
    }

    faces[5] = 1;
    data.normal_count = 5;
    if (bobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	payload.isValid() ||
	!payload.points.empty() ||
	!payload.normals.empty() ||
	!payload.coordIndex.empty()) {
	printf("FAIL: Obol mesh LoD payload copy accepted invalid normals\n");
	return 1;
    }

    data.normals = NULL;
    data.normal_count = 0;
    if (!bobol_lod_mesh_payload_from_mesh_lod_data(payload, data) ||
	!payload.isValid() ||
	!payload.normals.empty()) {
	printf("FAIL: Obol mesh LoD payload copy rejected mesh without normals\n");
	return 1;
    }

    return 0;
}

static int
test_worker_prepared_authored_normals(void)
{
    point_t points[4];
    vect_t normals[6];
    uint32_t faces[6] = {0, 1, 2, 0, 2, 3};
    VSET(points[0], 0.0, 0.0, 0.0);
    VSET(points[1], 1.0, 0.0, 0.0);
    VSET(points[2], 1.0, 1.0, 0.0);
    VSET(points[3], 0.0, 1.0, 0.0);
    for (size_t i = 0; i < 3; ++i)
	VSET(normals[i], 0.0, 0.0, 1.0);
    for (size_t i = 3; i < 6; ++i)
	VSET(normals[i], 0.0, 1.0, 1.0);

    struct BObolMeshLodData data = {};
    data.faces = faces;
    data.face_count = 2;
    data.points = points;
    data.point_count = 4;
    data.points_orig = points;
    data.point_orig_count = 4;
    data.normals = normals;
    data.normal_count = 6;
    VSET(data.bmin, 0.0, 0.0, 0.0);
    VSET(data.bmax, 1.0, 1.0, 0.0);

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    hierarchy.min_cut = 0;
    hierarchy.max_cut = 1;
    hierarchy.resident_cut = 1;
    hierarchy.cuts[0].face_count = 1;
    hierarchy.cuts[1].face_count = 2;
    hierarchy.cuts[0].point_count = 3;
    hierarchy.cuts[1].point_count = 4;
    VSET(hierarchy.quantization_min, 0.0, 0.0, 0.0);
    VSET(hierarchy.quantization_max, 1.0, 1.0, 0.0);
    complete_test_hierarchy(hierarchy);

    BObolLodProgressiveMesh progressive;
    if (!progressive.update(data, hierarchy, 1, FALSE)) {
	printf("FAIL: worker prepared authored-normal fixture update\n");
	return 1;
    }
    uint64_t revision = 0;
    const std::shared_ptr<const Obol::PartGeometry> prepared =
	progressive.prepareCadGeometry(
	    BOBOL_LOD_DRAW_SHADED, &revision);
    uint64_t cachedRevision = 0;
    const std::shared_ptr<const Obol::PartGeometry> cached =
	progressive.prepareCadGeometry(
	    BOBOL_LOD_DRAW_SHADED, &cachedRevision);
    const Obol::TriMesh *mesh =
	prepared && prepared->shaded ? &*prepared->shaded : NULL;
    BObolMeshLodCutInfo coarseCutInfo = {};
    if (!mesh || revision != progressive.revision() ||
	cachedRevision != revision || cached != prepared ||
	mesh->indices.size() != 6 ||
	mesh->normals.size() != mesh->positions.size() ||
	mesh->positions.size() <= 4 ||
	mesh->progressiveLineage != 0 ||
	!mesh->isProgressive() ||
	mesh->indexCountAtCut(0) != 3 ||
	mesh->indexCountAtCut(1) != 6 ||
	mesh->positionCountAtCut(0) != 3 ||
	mesh->positionCountAtCut(1) != mesh->positions.size() ||
	!progressive.cutInfo(0, &coarseCutInfo) ||
	fabs(coarseCutInfo.object_error - 1.0) > 1.0e-12 ||
	progressive.cutForScreenError(100.0, 80.0) != 0 ||
	progressive.cutForScreenError(100.0, 10.0) != 1 ||
	fabs(progressive.projectedErrorAtCut(0, 100.0) -
	    70.71067811865476) > 1.0e-9 ||
	fabs(progressive.projectedErrorAtCut(1, 100.0)) > 1.0e-12) {
	printf("FAIL: authored normals were not canonicalized into a cached "
	       "worker-owned renderer generation "
	       "(mesh=%d revision=%llu/%llu cached_revision=%llu cached=%d "
	       "indices=%zu normals=%zu positions=%zu lineage=%llu "
	       "progressive=%d index_cuts=%zu/%zu position_cuts=%zu/%zu "
	       "cut_info=%d error=%.17g screen_cuts=%d/%d projected=%.17g/%.17g)\n",
	       mesh ? 1 : 0,
	       static_cast<unsigned long long>(revision),
	       static_cast<unsigned long long>(progressive.revision()),
	       static_cast<unsigned long long>(cachedRevision),
	       cached == prepared ? 1 : 0,
	       mesh ? mesh->indices.size() : 0,
	       mesh ? mesh->normals.size() : 0,
	       mesh ? mesh->positions.size() : 0,
	       static_cast<unsigned long long>(
		   mesh ? mesh->progressiveLineage : 0),
	       mesh && mesh->isProgressive() ? 1 : 0,
	       mesh ? mesh->indexCountAtCut(0) : 0,
	       mesh ? mesh->indexCountAtCut(1) : 0,
	       mesh ? mesh->positionCountAtCut(0) : 0,
	       mesh ? mesh->positionCountAtCut(1) : 0,
	       progressive.cutInfo(0, &coarseCutInfo) ? 1 : 0,
	       coarseCutInfo.object_error,
	       progressive.cutForScreenError(100.0, 80.0),
	       progressive.cutForScreenError(100.0, 10.0),
	       progressive.projectedErrorAtCut(0, 100.0),
	       progressive.projectedErrorAtCut(1, 100.0));
	return 1;
    }
    return 0;
}

static int
test_worker_prepared_generation_lifetime(void)
{
    point_t points[4];
    uint32_t faces[6] = {0, 1, 2, 0, 2, 3};
    VSET(points[0], 0.0, 0.0, 0.0);
    VSET(points[1], 1.0, 0.0, 0.0);
    VSET(points[2], 1.0, 1.0, 0.0);
    VSET(points[3], 0.0, 1.0, 0.0);
    struct BObolMeshLodData data = {};
    data.faces = faces;
    data.face_count = 1;
    data.points = points;
    data.point_count = 3;
    data.points_orig = points;
    data.point_orig_count = 3;
    VSET(data.bmin, 0.0, 0.0, 0.0);
    VSET(data.bmax, 1.0, 1.0, 0.0);

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    hierarchy.min_cut = 0;
    hierarchy.max_cut = 1;
    hierarchy.resident_cut = 0;
    hierarchy.cuts[0].face_count = 1;
    hierarchy.cuts[1].face_count = 2;
    hierarchy.cuts[0].point_count = 3;
    hierarchy.cuts[1].point_count = 4;
    VSET(hierarchy.quantization_min, 0.0, 0.0, 0.0);
    VSET(hierarchy.quantization_max, 1.0, 1.0, 0.0);
    complete_test_hierarchy(hierarchy);

    BObolLodProgressiveMesh progressive;
    if (!progressive.update(data, hierarchy, 0, FALSE))
	return 1;
    uint64_t coarseRevision = 0;
    const std::shared_ptr<const Obol::PartGeometry> coarse =
	progressive.prepareCadGeometry(
	    BOBOL_LOD_DRAW_SHADED, &coarseRevision);
    const uint64_t progressiveLineage =
	coarse && coarse->shaded ? coarse->shaded->progressiveLineage : 0;

    data.face_count = 2;
    data.point_count = 4;
    data.point_orig_count = 4;
    hierarchy.resident_cut = 1;
    if (!progressive.update(data, hierarchy, 1, FALSE))
	return 1;
    uint64_t richRevision = 0;
    const std::shared_ptr<const Obol::PartGeometry> rich =
	progressive.prepareCadGeometry(
	    BOBOL_LOD_DRAW_SHADED, &richRevision);
    const size_t richBytes = progressive.estimateBytes();

    if (!coarse || !coarse->shaded || !rich || !rich->shaded ||
	coarseRevision == 0 || richRevision <= coarseRevision ||
	progressiveLineage == 0 ||
	rich->shaded->progressiveLineage != progressiveLineage ||
	coarse->shaded->indices.size() != 3 ||
	rich->shaded->indices.size() != 6 ||
	!progressive.trim(0)) {
	printf("FAIL: immutable progressive generation expansion\n");
	return 1;
    }

    uint64_t trimmedRevision = 0;
    const std::shared_ptr<const Obol::PartGeometry> trimmed =
	progressive.prepareCadGeometry(
	    BOBOL_LOD_DRAW_SHADED, &trimmedRevision);
    const size_t trimmedBytes = progressive.estimateBytes();
    if (!trimmed || !trimmed->shaded ||
	trimmedRevision <= richRevision ||
	trimmed->shaded->progressiveLineage != progressiveLineage ||
	trimmed->shaded->indices.size() != 3 ||
	trimmedBytes >= richBytes ||
	coarse->shaded->indices.size() != 3 ||
	rich->shaded->indices.size() != 6) {
	printf("FAIL: published progressive generation lifetime/memory across "
	       "trim (bytes=%zu/%zu)\n", trimmedBytes, richBytes);
	return 1;
    }
    return 0;
}

static int
test_compaction_preserves_coordinate_only_frontier(void)
{
    point_t points[3];
    uint32_t faces[3] = {0, 1, 2};
    VSET(points[0], 0.0, 0.0, 0.0);
    VSET(points[1], 1.0, 0.0, 0.0);
    VSET(points[2], 0.0, 1.0, 0.0);

    struct BObolMeshLodData data = {};
    data.faces = faces;
    data.face_count = 1;
    data.points = points;
    data.point_count = 3;
    data.points_orig = points;
    data.point_orig_count = 3;
    VSET(data.bmin, 0.0, 0.0, 0.0);
    VSET(data.bmax, 1.0, 1.0, 0.0);

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    hierarchy.min_cut = 0;
    hierarchy.max_cut = 2;
    hierarchy.resident_cut = 2;
    for (int cut = 0; cut <= hierarchy.max_cut; ++cut) {
	hierarchy.cuts[cut].face_count = 1;
	hierarchy.cuts[cut].point_count = 3;
    }
    VSET(hierarchy.quantization_min, 0.0, 0.0, 0.0);
    VSET(hierarchy.quantization_max, 1.0, 1.0, 0.0);
    complete_test_hierarchy(hierarchy);

    BObolLodProgressiveMesh progressive;
    if (!progressive.update(data, hierarchy, 2, FALSE) ||
	!progressive.canDrawCut(2) || !progressive.trim(0)) {
	printf("FAIL: coordinate-only compaction fixture setup\n");
	return 1;
    }

    const std::shared_ptr<const Obol::PartGeometry> trimmed =
	progressive.prepareCadGeometry(BOBOL_LOD_DRAW_SHADED, NULL);
    std::array<BObolLodCounts, BOBOL_MESH_LOD_CUT_COUNT_MAX> populations;
    int minimumCut = -1;
    int maximumDrawableCut = -1;
    const std::vector<uint32_t> wholeLeaf;
    if (progressive.residentCut() != 0 ||
	!progressive.canDrawCut(2) ||
	progressive.cutForScreenError(100.0, 1.0) != 2 ||
	!progressive.drawableCountsAtCuts(wholeLeaf, FALSE,
	    populations.data(), populations.size(), &minimumCut,
	    &maximumDrawableCut) ||
	minimumCut != 0 || maximumDrawableCut != 2 ||
	populations[0].faceCount != 1 ||
	populations[1].faceCount != 1 ||
	populations[2].faceCount != 1 ||
	populations[0].pointCount != 3 ||
	populations[1].pointCount != 3 ||
	populations[2].pointCount != 3 ||
	!trimmed || !trimmed->shaded ||
	trimmed->shaded->progressiveResidentCut != 2 ||
	trimmed->shaded->indices.size() != 3 ||
	trimmed->shaded->positions.size() != 3) {
	printf("FAIL: compaction or bulk census discarded the coordinate-only "
	       "drawable cut frontier (cuts=%d/%d faces=%llu/%llu/%llu)\n",
	       minimumCut, maximumDrawableCut,
	       static_cast<unsigned long long>(populations[0].faceCount),
	       static_cast<unsigned long long>(populations[1].faceCount),
	       static_cast<unsigned long long>(populations[2].faceCount));
	return 1;
    }
    return 0;
}

static int
test_progressive_generation_rejects_corruption(void)
{
    point_t points[3];
    vect_t normals[3];
    uint32_t faces[3] = {0, 1, 2};
    VSET(points[0], 0.0, 0.0, 0.0);
    VSET(points[1], 1.0, 0.0, 0.0);
    VSET(points[2], 0.0, 1.0, 0.0);
    for (size_t i = 0; i < 3; ++i)
	VSET(normals[i], 0.0, 0.0, 1.0);

    struct BObolMeshLodData data = {};
    data.faces = faces;
    data.face_count = 1;
    data.points = points;
    data.point_count = 3;
    data.points_orig = points;
    data.point_orig_count = 3;
    VSET(data.bmin, 0.0, 0.0, 0.0);
    VSET(data.bmax, 1.0, 1.0, 0.0);

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    hierarchy.min_cut = 0;
    hierarchy.max_cut = 0;
    hierarchy.resident_cut = 0;
    for (size_t level = 0; level < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++level) {
	hierarchy.cuts[level].face_count = 1;
	hierarchy.cuts[level].point_count = 3;
    }
    VSET(hierarchy.quantization_min, 0.0, 0.0, 0.0);
    VSET(hierarchy.quantization_max, 1.0, 1.0, 0.0);
    complete_test_hierarchy(hierarchy);

    BObolLodProgressiveMesh progressive;
    if (!progressive.update(data, hierarchy, 0, FALSE)) {
	printf("FAIL: progressive corruption fixture did not initialize\n");
	return 1;
    }
    const uint64_t validRevision = progressive.revision();

    points[2][X] = INFINITY;
    if (progressive.update(data, hierarchy, 0, FALSE) ||
	progressive.revision() != validRevision) {
	printf("FAIL: progressive generation accepted a non-finite point\n");
	return 1;
    }
    points[2][X] = 2.0;
    if (progressive.update(data, hierarchy, 0, FALSE) ||
	progressive.revision() != validRevision) {
	printf("FAIL: progressive generation accepted a point outside its "
	       "quantization domain\n");
	return 1;
    }
    points[2][X] = 0.0;
    faces[2] = 3;
    if (progressive.update(data, hierarchy, 0, FALSE) ||
	progressive.revision() != validRevision) {
	printf("FAIL: progressive generation accepted an invalid index\n");
	return 1;
    }
    faces[2] = 2;
    hierarchy.cuts[0].point_count = 2;
    if (progressive.update(data, hierarchy, 0, FALSE) ||
	progressive.revision() != validRevision) {
	printf("FAIL: progressive generation accepted mismatched hierarchy "
	       "counts\n");
	return 1;
    }
    hierarchy.cuts[0].point_count = 3;
    data.normals = normals;
    data.normal_count = 3;
    normals[1][Y] = NAN;
    if (progressive.update(data, hierarchy, 0, FALSE) ||
	progressive.revision() != validRevision) {
	printf("FAIL: progressive generation accepted a non-finite normal\n");
	return 1;
    }
    return 0;
}

static int
test_stage_results(void)
{
    BObolLodRequest request = make_request();
    std::vector<BObolLodDependency> dependencies;
    std::vector<BObolLodAttribute> attributes;
    BObolLodCounts counts;
    BObolLodProxy proxy;

    BObolLodDependency dependency;
    dependency.objectPath = "/all.g/child.bot";
    dependency.objectName = "child.bot";
    dependency.sourceRevision = 31;
    dependency.sourceContentHash = 0x77;
    dependency.requiredQualityTier = BOBOL_LOD_QUALITY_PROXY;
    dependency.optional = FALSE;
    dependencies.push_back(dependency);

    BObolLodResult result =
	bobol_lod_directory_result(request, dependencies);
    if (result.resultKind != BOBOL_LOD_RESULT_DIRECTORY ||
	result.payloadKind != BOBOL_LOD_PAYLOAD_DIRECTORY ||
	!result.payloadIsConsistent() ||
	result.qualityTier != BOBOL_LOD_QUALITY_METADATA ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	result.terminal ||
	result.dependencies.size() != 1 ||
	bu_strcmp(result.dependencies[0].objectPath.getString(),
	       "/all.g/child.bot") != 0 ||
	!bobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD directory stage result\n");
	return 1;
    }

    BObolLodAttribute color;
    color.name = "display.color";
    color.value = "255 0 0";
    attributes.push_back(color);
    result = bobol_lod_attributes_result(request, attributes);
    if (result.resultKind != BOBOL_LOD_RESULT_ATTRIBUTES ||
	result.payloadKind != BOBOL_LOD_PAYLOAD_ATTRIBUTES ||
	!result.payloadIsConsistent() ||
	result.qualityTier != BOBOL_LOD_QUALITY_ATTRIBUTES ||
	result.attributes.size() != 1 ||
	bu_strcmp(result.attributes[0].name.getString(), "display.color") != 0 ||
	!bobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD attributes stage result\n");
	return 1;
    }

    counts.lineCount = 12;
    result = bobol_lod_aabb_result(request, request.bounds, &counts);
    if (result.resultKind != BOBOL_LOD_RESULT_AABB ||
	result.payloadKind != BOBOL_LOD_PAYLOAD_PROXY ||
	!result.payloadIsConsistent() ||
	result.qualityTier != BOBOL_LOD_QUALITY_PROXY ||
	result.proxy.kind != BOBOL_LOD_PROXY_AABB ||
	!result.proxy.isValid() ||
	result.bounds.isEmpty() ||
	result.counts.lineCount != 12 ||
	!bobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD AABB stage result\n");
	return 1;
    }

    proxy.kind = BOBOL_LOD_PROXY_OBB;
    proxy.bounds = request.bounds;
    proxy.center.setValue(1.0f, 2.0f, 3.0f);
    proxy.axisX.setValue(1.0f, 0.0f, 0.0f);
    proxy.axisY.setValue(0.0f, 1.0f, 0.0f);
    proxy.axisZ.setValue(0.0f, 0.0f, 1.0f);
    proxy.halfExtents.setValue(4.0f, 5.0f, 6.0f);
    result = bobol_lod_proxy_result(request, proxy, NULL);
    if (result.resultKind != BOBOL_LOD_RESULT_PROXY ||
	result.payloadKind != BOBOL_LOD_PAYLOAD_PROXY ||
	!result.payloadIsConsistent() ||
	result.qualityTier != BOBOL_LOD_QUALITY_PROXY ||
	result.proxy.kind != BOBOL_LOD_PROXY_OBB ||
	!result.proxy.isValid() ||
	fabs((double)result.proxy.center[0] - 1.0) > 1.0e-6 ||
	fabs((double)result.proxy.halfExtents[2] - 6.0) > 1.0e-6 ||
	!bobol_lod_result_matches_request(result, request)) {
	printf("FAIL: LoD proxy stage result\n");
	return 1;
    }

    proxy.clear();
    result = bobol_lod_proxy_result(request, proxy, NULL);
    if (result.providerStatus != BOBOL_LOD_PROVIDER_ERROR ||
	result.payloadKind != BOBOL_LOD_PAYLOAD_STATUS ||
	!result.payloadIsConsistent() ||
	!result.terminal ||
	result.diagnostic.getLength() == 0) {
	printf("FAIL: LoD invalid proxy diagnostic result\n");
	return 1;
    }

    result.clear();
    result.addDependency("/all.g/leaf.bot", "leaf.bot", 41, 0x99,
			 BOBOL_LOD_QUALITY_FAST_DISPLAY, TRUE);
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

    /* Dormant shaded attributes are asset residency, not wire-channel work.
     * They must not change the allocator currency of an otherwise identical
     * wire cut, while shaded admission still accounts for their bandwidth. */
    BObolLodCounts normalCostCounts;
    normalCostCounts.faceCount = 100;
    normalCostCounts.pointCount = 75;
    normalCostCounts.normalCount = 300;
    BObolLodCounts noNormalCostCounts = normalCostCounts;
    noNormalCostCounts.normalCount = 0;
    if (bobol_lod_render_cost_units(
	    normalCostCounts, BOBOL_LOD_DRAW_WIRE) !=
	bobol_lod_render_cost_units(
	    noNormalCostCounts, BOBOL_LOD_DRAW_WIRE) ||
	bobol_lod_render_cost_units(
	    normalCostCounts, BOBOL_LOD_DRAW_SHADED) <=
	bobol_lod_render_cost_units(
	    noNormalCostCounts, BOBOL_LOD_DRAW_SHADED)) {
	printf("FAIL: LoD render-cost channel accounting\n");
	return 1;
    }

    if (test_key_determinism())
	return 1;
    if (test_rt_mesh_result())
	return 1;
    if (test_rt_mesh_payload_copy())
	return 1;
    if (test_worker_prepared_authored_normals())
	return 1;
    if (test_worker_prepared_generation_lifetime())
	return 1;
    if (test_compaction_preserves_coordinate_only_frontier())
	return 1;
    if (test_progressive_generation_rejects_corruption())
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
