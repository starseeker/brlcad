/*                L O D _ S E R V I C E . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/env.h"
#include "bu/log.h"
#include "bu/parallel.h"
#include "bu/str.h"
#include "bu/time.h"

#include "BObol/BDrawCache.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshLodCache.h"

#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <iomanip>
#include <limits>
#include <mutex>
#include <new>
#include <set>
#include <sstream>
#include <string.h>
#include <thread>
#include <unordered_map>
#include <vector>

BObolMeshLodProvider::BObolMeshLodProvider(void)
{
    clear();
}

void
BObolMeshLodProvider::clear(void)
{
    service = NULL;
    dbip = NULL;
    bv_view_info_init(&view);
    useView = FALSE;
    refreshMissing = TRUE;
    useForcedLevel = FALSE;
    shrinkAfterCopy = TRUE;
    compactResident = FALSE;
    progressiveDelivery = TRUE;
    initialRefinementFaceBudget = 250000;
    initialRefinementPointBudget = 500000;
    refinementGrowthFactor = 4.0;
    useCurrentDrawLevel = FALSE;
    currentDrawLevel = -1;
    forcedLevel = 0;
    reset = 0;
}

BObolRtSourceFullDetailProvider::BObolRtSourceFullDetailProvider(void)
{
    clear();
}

void
BObolRtSourceFullDetailProvider::clear(void)
{
    dbip = NULL;
    validateSourceMetrics = TRUE;
    maxFullDetailFaceCount = 0;
    maxFullDetailPointCount = 0;
}

BObolRtProxyProvider::BObolRtProxyProvider(void)
{
    clear();
}

void
BObolRtProxyProvider::clear(void)
{
    dbip = NULL;
    proxyKind = BOBOL_LOD_PROXY_AABB;
    useRequestBounds = TRUE;
}

BObolLodTask::BObolLodTask(void)
{
    clear();
}

void
BObolLodTask::clear(void)
{
    generation = 0;
    request.clear();
    dependencies.clear();
    realize = NULL;
    realizeData = NULL;
    realizeDataFree = NULL;
    cacheWrite = NULL;
    cacheWriteData = NULL;
    debugDelayMilliseconds = 0;
    estimatedWorkingSetBytes = 0;
    publishResult = TRUE;
    writeCache = FALSE;
}

void
BObolLodTask::addDependency(uint64_t taskId)
{
    if (taskId != 0)
	dependencies.push_back(taskId);
}

static const char *
lod_request_leaf_name(const char *name)
{
    if (!name)
	return NULL;

    const char *slash = strrchr(name, '/');
    if (slash && slash[1])
	return slash + 1;

    while (*name == '/')
	name++;
    return name[0] ? name : NULL;
}

static const char *
lod_request_object_name(const BObolLodRequest &request)
{
    const char *name = request.objectName.getString();
    const char *leaf = lod_request_leaf_name(name);
    if (leaf)
	return leaf;

    name = request.objectPath.getString();
    return lod_request_leaf_name(name);
}

static BObolLodResult
lod_provider_status_result(const BObolLodRequest &request, int status,
			   const char *diagnostic)
{
    BObolLodResult result;

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.qualityTier = request.qualityTier;
    result.providerStatus = status;
    result.terminal = TRUE;
    result.diagnostic = diagnostic ? diagnostic : "";
    if (status == BOBOL_LOD_PROVIDER_CACHE_MISS ||
	status == BOBOL_LOD_PROVIDER_STALE)
	result.stale = TRUE;

    return result;
}

static SbBool
lod_source_full_detail_exceeds_limits(
    const BObolRtSourceFullDetailProvider *provider,
    uint64_t faceCount, uint64_t pointCount)
{
    if (!provider)
	return FALSE;

    if (provider->maxFullDetailFaceCount != 0 &&
	faceCount > provider->maxFullDetailFaceCount)
	return TRUE;
    if (provider->maxFullDetailPointCount != 0 &&
	pointCount > provider->maxFullDetailPointCount)
	return TRUE;

    return FALSE;
}

static SbBool
lod_request_source_counts_known(const BObolLodRequest &request)
{
    return request.sourceCounts.faceCount != 0 ||
	   request.sourceCounts.pointCount != 0 ? TRUE : FALSE;
}

static BObolLodCounts
lod_counts_from_request(const BObolLodRequest &request)
{
    return request.sourceCounts;
}

static BObolLodResult
lod_aabb_result_from_record(const BObolLodRequest &request,
			    const BObolDrawProxyRecord &record,
			    const char *diagnostic)
{
    BObolLodCounts counts = lod_counts_from_request(request);

    if (record.kind != BOBOL_LOD_PROXY_AABB || record.pointCount != 2)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol AABB draw proxy record is invalid");

    SbBox3f bounds;
    bounds.makeEmpty();
    bounds.extendBy(SbVec3f(static_cast<float>(record.points[0][X]),
			    static_cast<float>(record.points[0][Y]),
			    static_cast<float>(record.points[0][Z])));
    bounds.extendBy(SbVec3f(static_cast<float>(record.points[1][X]),
			    static_cast<float>(record.points[1][Y]),
			    static_cast<float>(record.points[1][Z])));
    BObolLodResult result = bobol_lod_aabb_result(request, bounds,
			      &counts);
    if (diagnostic)
	result.diagnostic = diagnostic;
    return result;
}

static BObolLodResult
lod_aabb_result_from_request(const BObolLodRequest &request,
			     SbBool useRequestBounds)
{
    BObolLodCounts counts = lod_counts_from_request(request);
    if (!useRequestBounds || request.bounds.isEmpty())
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_CACHE_MISS,
					  "Obol AABB draw proxy cache entry unavailable");
    BObolLodResult result = bobol_lod_aabb_result(request,
			      request.bounds, &counts);
    result.diagnostic = "Obol AABB draw proxy using request bounds";
    return result;
}

static BObolLodResult
lod_aabb_result_from_cache_or_db(const BObolLodRequest &request,
				 struct db_i *dbip,
				 SbBool useRequestBounds)
{
    if (!dbip)
	return lod_aabb_result_from_request(request, useRequestBounds);

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol draw proxy provider request has no object name");

    BObolDrawProxyRecord record;
    bobol_draw_proxy_record_init(&record);
    if (bobol_draw_proxy_cache_get(dbip, name, BOBOL_LOD_PROXY_AABB,
				     &record) ==
	BRLCAD_OK)
	return lod_aabb_result_from_record(request, record,
					   "Obol AABB draw proxy loaded from cache");

    if (bobol_draw_proxy_cache_refresh(dbip, name,
					 BOBOL_LOD_PROXY_AABB, NULL) ==
	BRLCAD_OK &&
	bobol_draw_proxy_cache_get(dbip, name, BOBOL_LOD_PROXY_AABB,
				     &record) ==
	BRLCAD_OK)
	return lod_aabb_result_from_record(request, record,
					   "Obol AABB draw proxy generated and cached");

    return lod_aabb_result_from_request(request, useRequestBounds);
}

static SbBool
lod_obb_proxy_from_points(BObolLodProxy &proxy, const point_t *points,
			  size_t pointCount)
{
    if (!points || pointCount != 8)
	return FALSE;

    point_t center;
    VSETALL(center, 0.0);
    for (int i = 0; i < 8; i++)
	VADD2(center, center, points[i]);
    VSCALE(center, center, 1.0 / 8.0);

    vect_t xaxis, yaxis, zaxis;
    VSUB2(xaxis, points[1], points[0]);
    VSUB2(yaxis, points[3], points[0]);
    VSUB2(zaxis, points[4], points[0]);
    const fastf_t xlen = MAGNITUDE(xaxis);
    const fastf_t ylen = MAGNITUDE(yaxis);
    const fastf_t zlen = MAGNITUDE(zaxis);
    if (xlen <= 0.0 && ylen <= 0.0 && zlen <= 0.0)
	return FALSE;

    if (xlen > 0.0)
	VSCALE(xaxis, xaxis, 1.0 / xlen);
    else
	VSET(xaxis, 1.0, 0.0, 0.0);
    if (ylen > 0.0)
	VSCALE(yaxis, yaxis, 1.0 / ylen);
    else
	VSET(yaxis, 0.0, 1.0, 0.0);
    if (zlen > 0.0)
	VSCALE(zaxis, zaxis, 1.0 / zlen);
    else
	VSET(zaxis, 0.0, 0.0, 1.0);

    proxy.clear();
    proxy.kind = BOBOL_LOD_PROXY_OBB;
    proxy.center = SbVec3f(static_cast<float>(center[X]),
			   static_cast<float>(center[Y]),
			   static_cast<float>(center[Z]));
    proxy.axisX = SbVec3f(static_cast<float>(xaxis[X]),
			  static_cast<float>(xaxis[Y]),
			  static_cast<float>(xaxis[Z]));
    proxy.axisY = SbVec3f(static_cast<float>(yaxis[X]),
			  static_cast<float>(yaxis[Y]),
			  static_cast<float>(yaxis[Z]));
    proxy.axisZ = SbVec3f(static_cast<float>(zaxis[X]),
			  static_cast<float>(zaxis[Y]),
			  static_cast<float>(zaxis[Z]));
    proxy.halfExtents = SbVec3f(static_cast<float>(xlen * 0.5),
				static_cast<float>(ylen * 0.5),
				static_cast<float>(zlen * 0.5));
    proxy.bounds.makeEmpty();
    for (int i = 0; i < 8; i++) {
	proxy.bounds.extendBy(SbVec3f(static_cast<float>(points[i][X]),
				      static_cast<float>(points[i][Y]),
				      static_cast<float>(points[i][Z])));
    }

    return proxy.isValid();
}

static BObolLodResult
lod_obb_result_from_request(const BObolLodRequest &request,
			    SbBool useRequestBounds)
{
    if (!useRequestBounds || request.bounds.isEmpty())
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_CACHE_MISS,
				  "Obol OBB draw proxy cache entry unavailable");

    const SbVec3f &minimum = request.bounds.getMin();
    const SbVec3f &maximum = request.bounds.getMax();
    BObolLodProxy proxy;

    proxy.kind = BOBOL_LOD_PROXY_OBB;
    proxy.bounds = request.bounds;
    proxy.center = (minimum + maximum) * 0.5f;
    proxy.halfExtents = (maximum - minimum) * 0.5f;

    BObolLodCounts counts = lod_counts_from_request(request);
    BObolLodResult result = bobol_lod_proxy_result(request, proxy,
						&counts);
    result.diagnostic = "Obol OBB draw proxy using request bounds";
    return result;
}

static BObolLodResult
lod_obb_result_from_cache_or_db(const BObolLodRequest &request,
				struct db_i *dbip)
{
    if (!dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_CACHE_MISS,
					  "Obol OBB draw proxy cache entry unavailable");

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol draw proxy provider request has no object name");

    BObolDrawProxyRecord record;
    bobol_draw_proxy_record_init(&record);
    if (bobol_draw_proxy_cache_get(dbip, name, BOBOL_LOD_PROXY_OBB,
				     &record) !=
	BRLCAD_OK) {
	if (bobol_draw_proxy_cache_refresh(dbip, name,
					     BOBOL_LOD_PROXY_OBB, NULL) !=
	    BRLCAD_OK ||
	    bobol_draw_proxy_cache_get(dbip, name,
					 BOBOL_LOD_PROXY_OBB, &record) !=
	    BRLCAD_OK)
	    return lod_provider_status_result(request,
					      BOBOL_LOD_PROVIDER_CACHE_MISS,
					      "Obol OBB draw proxy cache entry unavailable");
    }

    BObolLodProxy proxy;
    if (!lod_obb_proxy_from_points(proxy, record.points, record.pointCount))
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol OBB draw proxy record is invalid");

    BObolLodCounts counts = lod_counts_from_request(request);
    BObolLodResult result = bobol_lod_proxy_result(request, proxy,
			      &counts);
    result.diagnostic = "Obol OBB draw proxy loaded from cache";
    return result;
}

static SbString
lod_vec3_provider_param(const SbVec3f &value)
{
    std::ostringstream out;

    out << std::setprecision(9)
	<< value[0] << " " << value[1] << " " << value[2];
    return SbString(out.str().c_str());
}

static SbString
lod_bounds_provider_param(const SbBox3f &bounds)
{
    std::ostringstream out;
    const SbVec3f &bmin = bounds.getMin();
    const SbVec3f &bmax = bounds.getMax();

    out << std::setprecision(9)
	<< bmin[0] << " " << bmin[1] << " " << bmin[2] << " "
	<< bmax[0] << " " << bmax[1] << " " << bmax[2];
    return SbString(out.str().c_str());
}

static SbString
lod_float_provider_param(float value)
{
    std::ostringstream out;

    out << std::setprecision(9) << value;
    return SbString(out.str().c_str());
}

static const BObolLodProviderParam *
lod_provider_param(const BObolLodRequest &request, const char *name)
{
    const BObolLodProviderParam *found = NULL;

    if (!name)
	return NULL;
    for (size_t i = 0; i < request.providerParams.size(); i++) {
	if (bu_strcmp(request.providerParams[i].name.getString(), name) != 0)
	    continue;
	if (found)
	    return NULL;
	found = &request.providerParams[i];
    }
    return found;
}

static void
lod_remove_source_query_provider_params(BObolLodRequest &request)
{
    request.providerParams.erase(
	std::remove_if(request.providerParams.begin(),
		       request.providerParams.end(),
    [](const BObolLodProviderParam &param) {
	return bu_strncmp(param.name.getString(), "source_query.",
		       13) == 0;
    }),
    request.providerParams.end());
}

static SbBool
lod_provider_param_has_no_trailing_tokens(std::istringstream &in)
{
    std::string extra;
    return in >> extra ? FALSE : TRUE;
}

static SbBool
lod_parse_float_provider_param(float &value, const SbString &text)
{
    std::istringstream in(text.getString());
    float parsed = 0.0f;
    if (!(in >> parsed))
	return FALSE;
    if (!std::isfinite(parsed) ||
	!lod_provider_param_has_no_trailing_tokens(in))
	return FALSE;
    value = parsed;
    return TRUE;
}

static SbBool
lod_parse_vec3_provider_param(SbVec3f &value, const SbString &text)
{
    std::istringstream in(text.getString());
    float v[3] = {0.0f, 0.0f, 0.0f};
    for (int i = 0; i < 3; i++) {
	if (!(in >> v[i]))
	    return FALSE;
	if (!std::isfinite(v[i]))
	    return FALSE;
    }
    if (!lod_provider_param_has_no_trailing_tokens(in))
	return FALSE;

    value.setValue(v[0], v[1], v[2]);
    return TRUE;
}

static SbBool
lod_parse_bounds_provider_param(SbBox3f &bounds, const SbString &text)
{
    std::istringstream in(text.getString());
    float v[6] = {0.0f, 0.0f, 0.0f, 0.0f, 0.0f, 0.0f};
    for (int i = 0; i < 6; i++) {
	if (!(in >> v[i]))
	    return FALSE;
	if (!std::isfinite(v[i]))
	    return FALSE;
    }
    if (!lod_provider_param_has_no_trailing_tokens(in))
	return FALSE;

    bounds.makeEmpty();
    bounds.extendBy(SbVec3f(v[0], v[1], v[2]));
    bounds.extendBy(SbVec3f(v[3], v[4], v[5]));
    return !bounds.isEmpty();
}

static SbBool
lod_bounds_intersect(const SbBox3f &a, const SbBox3f &b)
{
    if (a.isEmpty() || b.isEmpty())
	return FALSE;

    const SbVec3f amin = a.getMin();
    const SbVec3f amax = a.getMax();
    const SbVec3f bmin = b.getMin();
    const SbVec3f bmax = b.getMax();

    for (int axis = 0; axis < 3; axis++) {
	if (amax[axis] < bmin[axis] || bmax[axis] < amin[axis])
	    return FALSE;
    }
    return TRUE;
}

static SbBool
lod_request_query_space_is_source_local(const BObolLodRequest &request)
{
    const BObolLodProviderParam *spaceParam =
	lod_provider_param(request, "source_query.space");
    return spaceParam &&
	   bu_strcmp(spaceParam->value.getString(), "source_local") == 0 ?
	   TRUE : FALSE;
}

static SbBool
lod_request_snap_query_bounds(const BObolLodRequest &request,
			      SbBox3f &queryBounds)
{
    if (!lod_request_query_space_is_source_local(request))
	return FALSE;

    const BObolLodProviderParam *boundsParam =
	lod_provider_param(request, "source_query.bounds");
    const BObolLodProviderParam *toleranceParam =
	lod_provider_param(request, "source_query.tolerance");
    if (!boundsParam || !toleranceParam)
	return FALSE;

    float tolerance = 0.0f;
    if (!lod_parse_float_provider_param(tolerance, toleranceParam->value) ||
	tolerance < 0.0f)
	return FALSE;

    return lod_parse_bounds_provider_param(queryBounds, boundsParam->value);
}

static SbBool
lod_request_pick_query_ray(const BObolLodRequest &request,
			   SbVec3f &rayOrigin,
			   SbVec3f &rayDirection)
{
    if (!lod_request_query_space_is_source_local(request))
	return FALSE;

    const BObolLodProviderParam *originParam =
	lod_provider_param(request, "source_query.ray.origin");
    const BObolLodProviderParam *directionParam =
	lod_provider_param(request, "source_query.ray.direction");
    if (!originParam || !directionParam)
	return FALSE;

    if (!lod_parse_vec3_provider_param(rayOrigin, originParam->value) ||
	!lod_parse_vec3_provider_param(rayDirection,
				       directionParam->value) ||
	rayDirection.length() <= 0.0f)
	return FALSE;

    rayDirection.normalize();
    return TRUE;
}

static SbBool
lod_request_has_scoped_subset_query(const BObolLodRequest &request)
{
    SbBox3f queryBounds;
    SbVec3f rayOrigin;
    SbVec3f rayDirection;
    const SbBool hasBounds =
	lod_request_snap_query_bounds(request, queryBounds);
    const SbBool hasRay =
	lod_request_pick_query_ray(request, rayOrigin, rayDirection);

    return hasBounds != hasRay ? TRUE : FALSE;
}

static SbBool
lod_ray_intersects_triangle(const SbVec3f &origin,
			    const SbVec3f &direction,
			    const SbVec3f &a,
			    const SbVec3f &b,
			    const SbVec3f &c)
{
    const float epsilon = 1.0e-7f;
    const SbVec3f ab = b - a;
    const SbVec3f ac = c - a;
    const SbVec3f pvec = direction.cross(ac);
    const float det = ab.dot(pvec);
    if (det > -epsilon && det < epsilon)
	return FALSE;

    const float invDet = 1.0f / det;
    const SbVec3f tvec = origin - a;
    const float u = tvec.dot(pvec) * invDet;
    if (u < 0.0f || u > 1.0f)
	return FALSE;

    const SbVec3f qvec = tvec.cross(ab);
    const float v = direction.dot(qvec) * invDet;
    if (v < 0.0f || u + v > 1.0f)
	return FALSE;

    const float t = ac.dot(qvec) * invDet;
    return t >= 0.0f ? TRUE : FALSE;
}

static BObolLodResult
lod_source_full_detail_payload_result(const BObolLodRequest &request,
				      const struct rt_bot_internal *bot)
{
    BObolLodResult result;
    SbBox3f queryBounds;
    SbBool useQueryBounds =
	lod_request_snap_query_bounds(request, queryBounds);
    SbVec3f queryRayOrigin;
    SbVec3f queryRayDirection;
    SbBool useQueryRay =
	lod_request_pick_query_ray(request, queryRayOrigin, queryRayDirection);
    std::vector<size_t> selectedFaces;

    if (useQueryBounds && useQueryRay) {
	useQueryBounds = FALSE;
	useQueryRay = FALSE;
    }

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_FULL_DETAIL;
    result.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = TRUE;

    result.geometry.kind = BOBOL_LOD_GEOMETRY_OBOL_MESH;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = result.cacheKey;
    result.geometry.activeLevel = -1;
    result.geometry.borrowed = FALSE;

    result.bounds.makeEmpty();
    try {
	std::vector<SbVec3f> sourcePoints;
	sourcePoints.reserve(bot->num_vertices);
	for (size_t i = 0; i < bot->num_vertices; i++) {
	    SbVec3f point(static_cast<float>(bot->vertices[i * 3]),
			  static_cast<float>(bot->vertices[i * 3 + 1]),
			  static_cast<float>(bot->vertices[i * 3 + 2]));
	    sourcePoints.push_back(point);
	}

	selectedFaces.reserve(bot->num_faces);
	for (size_t i = 0; i < bot->num_faces; i++) {
	    int ia = bot->faces[i * 3];
	    int ib = bot->faces[i * 3 + 1];
	    int ic = bot->faces[i * 3 + 2];
	    if (ia < 0 || ib < 0 || ic < 0 ||
		static_cast<size_t>(ia) >= bot->num_vertices ||
		static_cast<size_t>(ib) >= bot->num_vertices ||
		static_cast<size_t>(ic) >= bot->num_vertices) {
		result.mesh.clear();
		result.providerStatus = BOBOL_LOD_PROVIDER_ERROR;
		result.diagnostic =
		    "RT source full-detail provider BoT has invalid face indices";
		return result;
	    }

	    if (useQueryBounds) {
		SbBox3f faceBounds;
		faceBounds.makeEmpty();
		faceBounds.extendBy(sourcePoints[static_cast<size_t>(ia)]);
		faceBounds.extendBy(sourcePoints[static_cast<size_t>(ib)]);
		faceBounds.extendBy(sourcePoints[static_cast<size_t>(ic)]);
		if (!lod_bounds_intersect(faceBounds, queryBounds))
		    continue;
	    }
	    if (useQueryRay &&
		!lod_ray_intersects_triangle(queryRayOrigin,
					     queryRayDirection,
					     sourcePoints[static_cast<size_t>(ia)],
					     sourcePoints[static_cast<size_t>(ib)],
					     sourcePoints[static_cast<size_t>(ic)]))
		continue;
	    selectedFaces.push_back(i);
	}

	if (selectedFaces.empty() && (useQueryBounds || useQueryRay)) {
	    result.mesh.clear();
	    result.geometry.clear();
	    result.bounds.makeEmpty();
	    result.counts.clear();
	    result.resultKind = BOBOL_LOD_RESULT_NONE;
	    result.providerStatus = BOBOL_LOD_PROVIDER_FALLBACK;
	    result.diagnostic =
		"RT source full-detail provider scoped query matched no faces";
	    return result;
	}

	if (selectedFaces.empty()) {
	    selectedFaces.reserve(bot->num_faces);
	    for (size_t i = 0; i < bot->num_faces; i++)
		selectedFaces.push_back(i);
	}

	result.counts.faceCount = selectedFaces.size();
	result.mesh.coordIndex.reserve(selectedFaces.size() * 3);
	result.mesh.faceIndex.reserve(selectedFaces.size());
	if (selectedFaces.size() < bot->num_faces) {
	    std::vector<int32_t> sourceToLocal(bot->num_vertices, -1);
	    result.mesh.points.reserve(std::min(bot->num_vertices,
						selectedFaces.size() * 3));
	    result.mesh.vertexIndex.reserve(result.mesh.points.capacity());
	    for (size_t i = 0; i < selectedFaces.size(); i++) {
		size_t faceIndex = selectedFaces[i];
		result.mesh.faceIndex.push_back(static_cast<int32_t>(faceIndex));
		for (size_t j = 0; j < 3; j++) {
		    const int sourceIndex = bot->faces[faceIndex * 3 + j];
		    const size_t sourceSlot = static_cast<size_t>(sourceIndex);
		    int32_t localIndex = sourceToLocal[sourceSlot];
		    if (localIndex < 0) {
			localIndex =
			    static_cast<int32_t>(result.mesh.points.size());
			sourceToLocal[sourceSlot] = localIndex;
			result.mesh.points.push_back(sourcePoints[sourceSlot]);
			result.mesh.vertexIndex.push_back(
			    static_cast<int32_t>(sourceIndex));
		    }
		    result.mesh.coordIndex.push_back(localIndex);
		    result.bounds.extendBy(
			result.mesh.points[static_cast<size_t>(localIndex)]);
		}
	    }
	} else {
	    result.mesh.points.swap(sourcePoints);
	    for (size_t i = 0; i < selectedFaces.size(); i++) {
		size_t faceIndex = selectedFaces[i];
		result.mesh.faceIndex.push_back(static_cast<int32_t>(faceIndex));
		for (size_t j = 0; j < 3; j++) {
		    int idx = bot->faces[faceIndex * 3 + j];
		    result.mesh.coordIndex.push_back(static_cast<int32_t>(idx));
		    result.bounds.extendBy(
			result.mesh.points[static_cast<size_t>(idx)]);
		}
	    }
	}
	result.counts.pointCount = result.mesh.points.size();
    } catch (const std::bad_alloc &) {
	result.mesh.clear();
	result.providerStatus = BOBOL_LOD_PROVIDER_FALLBACK;
	result.diagnostic =
	    "RT source full-detail provider could not allocate BoT payload";
	return result;
    }

    if (!result.mesh.isValid()) {
	result.providerStatus = BOBOL_LOD_PROVIDER_ERROR;
	result.diagnostic =
	    "RT source full-detail provider copied an invalid BoT payload";
    }

    return result;
}

BObolLodResult
bobol_rt_source_full_detail_provider_task(
    const BObolLodRequest &request, void *userData)
{
    BObolRtSourceFullDetailProvider *provider =
	static_cast<BObolRtSourceFullDetailProvider *>(userData);
    if (!provider || !provider->dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider has no database");

    const SbBool scopedSubsetRequest =
	lod_request_has_scoped_subset_query(request);

    if (!scopedSubsetRequest && lod_request_source_counts_known(request) &&
	lod_source_full_detail_exceeds_limits(provider,
		request.sourceCounts.faceCount, request.sourceCounts.pointCount))
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_FALLBACK,
					  "RT source full-detail provider request exceeds full-detail limits");

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider request has no object name");

    struct directory *dp = db_lookup(provider->dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider could not find source object");

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    int internalType = rt_db_get_internal(&intern, dp, provider->dbip, NULL);
    if (internalType < 0)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider could not read source object");

    if (internalType != ID_BOT || intern.idb_type != ID_BOT ||
	intern.idb_ptr == NULL) {
	rt_db_free_internal(&intern);
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider source is not a BoT");
    }

    const struct rt_bot_internal *bot =
	    static_cast<const struct rt_bot_internal *>(intern.idb_ptr);
    if (!bot || !bot->vertices || !bot->faces ||
	bot->num_vertices == 0 || bot->num_faces == 0) {
	rt_db_free_internal(&intern);
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider source BoT has no mesh payload");
    }
    RT_BOT_CK_MAGIC(bot);

    if (!scopedSubsetRequest && lod_source_full_detail_exceeds_limits(provider,
	    bot->num_faces, bot->num_vertices)) {
	rt_db_free_internal(&intern);
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_FALLBACK,
					  "RT source full-detail provider source exceeds full-detail limits");
    }

    if (provider->validateSourceMetrics &&
	((request.sourceCounts.faceCount != 0 &&
	  request.sourceCounts.faceCount != bot->num_faces) ||
	 (request.sourceCounts.pointCount != 0 &&
	  request.sourceCounts.pointCount != bot->num_vertices))) {
	rt_db_free_internal(&intern);
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_STALE,
					  "RT source full-detail provider source metrics changed");
    }

    if (bot->num_vertices >
	static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	bot->num_faces >
	static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	bot->num_faces >
	static_cast<size_t>(std::numeric_limits<size_t>::max() / 3)) {
	rt_db_free_internal(&intern);
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_FALLBACK,
					  "RT source full-detail provider source exceeds copy limits");
    }

    BObolLodResult result =
	lod_source_full_detail_payload_result(request, bot);
    if (result.providerStatus == BOBOL_LOD_PROVIDER_READY &&
	lod_source_full_detail_exceeds_limits(provider,
		result.counts.faceCount, result.counts.pointCount)) {
	result.mesh.clear();
	result.geometry.clear();
	result.bounds.makeEmpty();
	result.counts.clear();
	result.resultKind = BOBOL_LOD_RESULT_NONE;
	result.providerStatus = BOBOL_LOD_PROVIDER_FALLBACK;
	result.diagnostic =
	    "RT source full-detail provider request exceeds full-detail limits";
    }
    rt_db_free_internal(&intern);
    return result;
}

void
bobol_rt_source_full_detail_provider_free(void *userData)
{
    BObolRtSourceFullDetailProvider *provider =
	static_cast<BObolRtSourceFullDetailProvider *>(userData);
    delete provider;
}

SbBool
bobol_lod_rt_source_full_detail_request_from_source_mesh_request(
    BObolLodRequest &request,
    const BObolSourceMeshRequest &sourceRequest,
    const BObolLodRequest *templateRequest)
{
    if (sourceRequest.meshAssetPath.getLength() == 0 &&
	sourceRequest.meshAssetName.getLength() == 0 &&
	sourceRequest.path.getLength() == 0 &&
	sourceRequest.sourceName.getLength() == 0)
	return FALSE;

    if (templateRequest)
	request = *templateRequest;
    else
	request.clear();
    lod_remove_source_query_provider_params(request);

    request.objectPath = sourceRequest.meshAssetPath.getLength() > 0 ?
	sourceRequest.meshAssetPath :
	(sourceRequest.path.getLength() > 0 ?
	 sourceRequest.path : sourceRequest.sourceName);
    request.objectName = sourceRequest.meshAssetName.getLength() > 0 ?
	sourceRequest.meshAssetName : sourceRequest.sourceName;
    if (request.objectName.getLength() == 0) {
	const char *name = lod_request_object_name(request);
	request.objectName = name ? name : "";
    }

    request.providerId = "rt_source_full_detail";
    request.providerVersion = "direct-bot-v1";
    request.qualityTier = BOBOL_LOD_QUALITY_FULL_DETAIL;
    if (request.drawMode == BOBOL_LOD_DRAW_UNKNOWN)
	request.drawMode = BOBOL_LOD_DRAW_SHADED;
    request.bounds = !sourceRequest.meshAssetBounds.isEmpty() ?
	sourceRequest.meshAssetBounds : sourceRequest.bounds;
    request.sourceCounts.clear();
    request.sourceCounts.faceCount = sourceRequest.faceCount;
    request.sourceCounts.pointCount = sourceRequest.pointCount;
    if ((sourceRequest.queryBoundsValid && !sourceRequest.queryBounds.isEmpty()) ||
	sourceRequest.queryRayValid || sourceRequest.queryToleranceValid)
	request.addProviderParam("source_query.space", "source_local");
    if (sourceRequest.queryBoundsValid && !sourceRequest.queryBounds.isEmpty()) {
	request.addProviderParam("source_query.bounds",
				 lod_bounds_provider_param(sourceRequest.queryBounds));
    }
    if (sourceRequest.queryRayValid) {
	request.addProviderParam("source_query.ray.origin",
				 lod_vec3_provider_param(sourceRequest.queryRayOrigin));
	request.addProviderParam("source_query.ray.direction",
				 lod_vec3_provider_param(sourceRequest.queryRayDirection));
    }
    if (sourceRequest.queryToleranceValid) {
	request.addProviderParam("source_query.tolerance",
				 lod_float_provider_param(sourceRequest.queryTolerance));
    }

    return request.objectPath.getLength() > 0 ||
	   request.objectName.getLength() > 0 ? TRUE : FALSE;
}

uint64_t
bobol_lod_submit_rt_source_full_detail_request(
    BObolLodService *service,
    uint64_t generation,
    const BObolSourceMeshRequest &sourceRequest,
    struct db_i *dbip,
    const BObolLodRequest *templateRequest,
    uint64_t maxFullDetailFaceCount,
    uint64_t maxFullDetailPointCount)
{
    if (!service || !dbip)
	return 0;

    BObolRtSourceFullDetailProvider *provider =
	new (std::nothrow) BObolRtSourceFullDetailProvider;
    if (!provider)
	return 0;

    BObolLodTask task;
    task.generation = generation;
    if (!bobol_lod_rt_source_full_detail_request_from_source_mesh_request(
	    task.request, sourceRequest, templateRequest)) {
	delete provider;
	return 0;
    }

    provider->dbip = dbip;
    provider->validateSourceMetrics = TRUE;
    provider->maxFullDetailFaceCount = maxFullDetailFaceCount;
    provider->maxFullDetailPointCount = maxFullDetailPointCount;

    task.realize = bobol_rt_source_full_detail_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = bobol_rt_source_full_detail_provider_free;

    uint64_t taskId = service->submitIfNotActive(task);
    if (taskId == 0)
	bobol_rt_source_full_detail_provider_free(provider);

    return taskId;
}

BObolLodResult
bobol_mesh_lod_provider_task(const BObolLodRequest &request,
			       void *userData)
{
    BObolMeshLodProvider *provider =
	static_cast<BObolMeshLodProvider *>(userData);
    if (!provider || !provider->dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD provider has no database");
    if (provider->service)
	return provider->service->realizeResidentMeshLod(request, *provider);

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD provider request has no object name");

    struct BObolMeshLodCacheStatus status =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_status(provider->dbip, name, &status) != BRLCAD_OK)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD provider could not query cache status");

    if ((!status.has_cache_key || !status.has_cached_payload ||
	 status.stale_cache_entry) && provider->refreshMissing) {
	if (bobol_mesh_lod_cache_refresh(provider->dbip, name, &status) != BRLCAD_OK)
	    return lod_provider_status_result(request,
					      BOBOL_LOD_PROVIDER_CACHE_MISS,
					      "Obol mesh LoD provider could not refresh cache entry");
    }

    struct BObolMeshLod *lod = bobol_mesh_lod_get(provider->dbip, name);
    if (!lod) {
	std::ostringstream diagnostic;
	diagnostic << "Obol mesh LoD provider has no cache payload for "
		   << name << " (cache_key=" << status.cache_key
		   << ", has_key=" << status.has_cache_key
		   << ", has_payload=" << status.has_cached_payload
		   << ", stale=" << status.stale_cache_entry << ")";
	return lod_provider_status_result(request,
					  status.stale_cache_entry ? BOBOL_LOD_PROVIDER_STALE :
					  BOBOL_LOD_PROVIDER_CACHE_MISS,
					  diagnostic.str().c_str());
    }

    int load_ret = provider->useForcedLevel ?
		   bobol_mesh_lod_load_level(lod, provider->forcedLevel, provider->reset) :
		   (request.requestedLevel >= 0 ?
		    bobol_mesh_lod_load_display_level(lod,
			request.requestedLevel, provider->reset) :
		    (provider->useView ?
		     bobol_mesh_lod_load_view(lod, &provider->view, provider->reset) :
		     bobol_mesh_lod_load_view(lod, NULL, provider->reset)));
    if (load_ret < 0) {
	bobol_mesh_lod_destroy(lod);
	return lod_provider_status_result(request,
					  BOBOL_LOD_PROVIDER_CACHE_MISS,
					  "Obol mesh LoD provider could not load a view level");
    }

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    int have_info = bobol_mesh_lod_info_get(lod, &info);
    const char *traceFilter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (traceFilter && traceFilter[0] &&
	((name && strstr(name, traceFilter)) ||
	 (request.objectPath.getLength() > 0 &&
	  strstr(request.objectPath.getString(), traceFilter)))) {
	bu_log("BObol LoD provider trace object=%s request_level=%d "
	       "loaded_level=%d faces=%zu points=%zu have_info=%d "
	       "view_revision=%llu policy_revision=%llu\n",
	       name ? name : "", request.requestedLevel, load_ret,
	       info.face_count, info.point_count, have_info,
	       static_cast<unsigned long long>(request.viewRevision),
	       static_cast<unsigned long long>(request.policyRevision));
    }
    if (!bobol_mesh_lod_has_active_data(lod)) {
	BObolLodResult result =
	    bobol_lod_result_from_mesh_lod_info(request, info, &status);
	bobol_mesh_lod_destroy(lod);
	result.providerStatus = BOBOL_LOD_PROVIDER_CACHE_MISS;
	result.diagnostic = "Obol mesh LoD provider loaded no active mesh data";
	return result;
    }
    if (!have_info) {
	bobol_mesh_lod_destroy(lod);
	return lod_provider_status_result(request,
					  BOBOL_LOD_PROVIDER_CACHE_MISS,
					  "Obol mesh LoD provider loaded no mesh metadata");
    }

    BObolLodResult result =
	bobol_lod_result_from_mesh_lod_info(request, info, &status);
    if (result.providerStatus == BOBOL_LOD_PROVIDER_READY) {
	struct BObolMeshLodData data;
	if (!bobol_mesh_lod_data_get(lod, &data) ||
	    !bobol_lod_mesh_payload_from_mesh_lod_data(result.mesh, data)) {
	    bobol_mesh_lod_destroy(lod);
	    return lod_provider_status_result(request,
					      BOBOL_LOD_PROVIDER_CACHE_MISS,
					      "Obol mesh LoD provider could not copy active mesh payload");
	}
	if (provider->shrinkAfterCopy)
	    bobol_mesh_lod_memshrink(lod);
    }

    bobol_mesh_lod_destroy(lod);
    return result;
}

BObolLodResult
bobol_mesh_lod_cache_provider_task(const BObolLodRequest &request,
				     void *userData)
{
    BObolMeshLodProvider *provider =
	static_cast<BObolMeshLodProvider *>(userData);
    if (!provider || !provider->dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD cache provider has no database");

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD cache provider request has no object name");

    struct BObolMeshLodCacheStatus status =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_status(provider->dbip, name, &status) != BRLCAD_OK)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD cache provider could not query cache status");

    if ((!status.has_cache_key || !status.has_cached_payload ||
	 status.stale_cache_entry) && provider->refreshMissing) {
	if (bobol_mesh_lod_cache_refresh(provider->dbip, name, &status) != BRLCAD_OK)
	    return lod_provider_status_result(request,
					      BOBOL_LOD_PROVIDER_CACHE_MISS,
					      "Obol mesh LoD cache provider could not refresh cache entry");
    }

    BObolLodResult result;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_DIAGNOSTIC;
    result.qualityTier = request.qualityTier;
    result.providerStatus =
	(status.has_cache_key && status.has_cached_payload &&
	 !status.stale_cache_entry) ? BOBOL_LOD_PROVIDER_READY :
	(status.stale_cache_entry ? BOBOL_LOD_PROVIDER_STALE :
	 BOBOL_LOD_PROVIDER_CACHE_MISS);
    result.terminal = TRUE;
    result.geometry.kind = BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.providerToken = status.cache_key;
    result.geometry.cacheKey = result.cacheKey;
    if (result.providerStatus == BOBOL_LOD_PROVIDER_READY)
	result.diagnostic = "Obol mesh LoD cache entry ready";
    else
	result.diagnostic = "Obol mesh LoD cache entry unavailable";
    return result;
}

BObolLodResult
bobol_rt_proxy_provider_task(const BObolLodRequest &request,
			       void *userData)
{
    BObolRtProxyProvider *provider =
	static_cast<BObolRtProxyProvider *>(userData);
    if (!provider)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol draw proxy provider has no provider state");

    if (provider->proxyKind == BOBOL_LOD_PROXY_AABB)
	return lod_aabb_result_from_cache_or_db(request, provider->dbip,
						provider->useRequestBounds);

    if (provider->proxyKind != BOBOL_LOD_PROXY_OBB)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol draw proxy provider has unknown proxy kind");

    BObolLodResult result = lod_obb_result_from_cache_or_db(request,
			      provider->dbip);
    if (result.providerStatus == BOBOL_LOD_PROVIDER_READY)
	return result;

    BObolLodResult fallback = lod_obb_result_from_request(
	request, provider->useRequestBounds);
    return fallback.providerStatus == BOBOL_LOD_PROVIDER_READY ?
	fallback : result;
}

void
bobol_rt_proxy_provider_free(void *userData)
{
    BObolRtProxyProvider *provider =
	static_cast<BObolRtProxyProvider *>(userData);
    delete provider;
}

void
bobol_mesh_lod_provider_free(void *userData)
{
    BObolMeshLodProvider *provider =
	static_cast<BObolMeshLodProvider *>(userData);
    delete provider;
}

struct BObolLodWorkItem {
    uint64_t id;
    BObolLodTask task;
};

struct BObolLodCacheWriteItem {
    BObolLodResult result;
    BObolLodCacheWriteProc write;
    void *writeData;
};

struct BObolLodSubscriber {
    BObolLodSubscriber(void) :
	id(0),
	callback(NULL),
	userData(NULL),
	active(FALSE),
	inFlight(0)
    {
    }

    BObolLodSubscriberId id;
    BObolLodResultReadyCB callback;
    void *userData;
    SbBool active;
    size_t inFlight;
};

static size_t
lod_default_working_set_limit(void)
{
    const size_t mebibyte = 1024ULL * 1024ULL;
    const size_t gibibyte = 1024ULL * mebibyte;
    size_t totalBytes = 0;
    size_t availableBytes = 0;
    const bool haveTotal = bu_mem(BU_MEM_ALL, &totalBytes) >= 0 &&
	totalBytes > 0;
    const bool haveAvailable = bu_mem(BU_MEM_AVAIL, &availableBytes) >= 0 &&
	availableBytes > 0;
    size_t allowance = gibibyte;
    if (haveTotal)
	allowance = std::min(allowance,
	    std::max(mebibyte, totalBytes / 8));
    if (haveAvailable)
	allowance = std::min(allowance,
	    std::max(mebibyte, availableBytes / 4));
    return std::max(mebibyte, allowance);
}

struct BObolGlobalWorkingSetGovernor {
    BObolGlobalWorkingSetGovernor(void) :
	limit(lod_default_working_set_limit()),
	activeBytes(0),
	activeTasks(0),
	peakBytes(0),
	peakTasks(0)
    {
    }

    std::mutex mutex;
    std::condition_variable cv;
    size_t limit;
    size_t activeBytes;
    size_t activeTasks;
    size_t peakBytes;
    size_t peakTasks;
};

static BObolGlobalWorkingSetGovernor &
lod_global_working_set_governor(void)
{
    static BObolGlobalWorkingSetGovernor governor;
    return governor;
}

void
bobol_lod_working_set_acquire(size_t estimatedBytes)
{
    if (!estimatedBytes)
	return;
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::unique_lock<std::mutex> lock(governor.mutex);
    governor.cv.wait(lock, [&]() {
	if (governor.limit == SIZE_MAX)
	    return true;
	/* An exceptional single mesh still has to make progress. */
	if (governor.activeBytes == 0)
	    return true;
	const size_t occupied = std::min(governor.activeBytes,
	    governor.limit);
	return estimatedBytes <= governor.limit - occupied;
    });
    governor.activeBytes =
	estimatedBytes > SIZE_MAX - governor.activeBytes ?
	SIZE_MAX : governor.activeBytes + estimatedBytes;
    governor.activeTasks++;
    governor.peakBytes = std::max(governor.peakBytes,
	governor.activeBytes);
    governor.peakTasks = std::max(governor.peakTasks,
	governor.activeTasks);
}

void
bobol_lod_working_set_release(size_t estimatedBytes)
{
    if (!estimatedBytes)
	return;
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    {
	std::lock_guard<std::mutex> lock(governor.mutex);
	governor.activeBytes =
	    estimatedBytes >= governor.activeBytes ?
	    0 : governor.activeBytes - estimatedBytes;
	if (governor.activeTasks > 0)
	    governor.activeTasks--;
    }
    governor.cv.notify_all();
}

size_t
bobol_lod_working_set_global_limit(void)
{
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::lock_guard<std::mutex> lock(governor.mutex);
    return governor.limit;
}

size_t
bobol_lod_working_set_global_active_bytes(void)
{
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::lock_guard<std::mutex> lock(governor.mutex);
    return governor.activeBytes;
}

size_t
bobol_lod_working_set_global_peak_bytes(void)
{
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::lock_guard<std::mutex> lock(governor.mutex);
    return governor.peakBytes;
}

size_t
bobol_lod_working_set_global_active_tasks(void)
{
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::lock_guard<std::mutex> lock(governor.mutex);
    return governor.activeTasks;
}

size_t
bobol_lod_working_set_global_peak_tasks(void)
{
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::lock_guard<std::mutex> lock(governor.mutex);
    return governor.peakTasks;
}

struct BObolResidentMeshAsset {
    ~BObolResidentMeshAsset()
    {
	if (lod)
	    bobol_mesh_lod_destroy(lod);
	lod = NULL;
    }

    std::mutex mutex;
    struct db_i *dbip = NULL;
    std::string name;
    struct BObolMeshLod *lod = NULL;
    BObolLodProgressiveMeshPtr mesh;
    struct BObolMeshLodCacheStatus status =
	BOBOL_MESH_LOD_CACHE_STATUS_INIT;
};

struct BObolResidentMeshConsumerDemand {
    uint64_t revision = 0;
    std::unordered_map<std::string, int> levels;
};

struct BObolLodServicePrivate {
    explicit BObolLodServicePrivate(BObolLodService *newOwner) :
	owner(newOwner),
	running(FALSE),
	stopping(FALSE),
	cacheWriterStopping(FALSE),
	cacheWriterEnabled(FALSE),
	nextTaskId(1),
	nextSubscriberId(1),
	nextGeneration(0),
	activeGeneration(0),
	maxActiveTasks(1024),
	maxQueuedResults(256),
	maxQueuedCacheWrites(256),
	maxActiveWorkingSetBytes(0),
	activeWorkingSetBytes(0),
	executingTasks(0),
	peakWorkingSetBytes(0),
	peakExecutingTasks(0),
	resultReservations(0),
	cacheWriteReservations(0),
	rejectedTasks(0),
	coalescedResults(0),
	coalescedCacheWrites(0),
	discardedStaleResults(0),
	residentMeshCacheLoads(0),
	residentMeshHits(0),
	residentMeshCompactions(0),
	inFlight(0),
	cacheWriteInFlight(0),
	delayedTasks(0)
    {
	/* This is a concurrent transient-work allowance, not a promise that an
	 * individual mesh can be realized in bounded memory.  The latter needs
	 * streaming/external construction in the provider.  Size the aggregate
	 * allowance from both installed and currently available RAM so adding
	 * CPU workers on a small or already-busy host cannot multiply several
	 * very large topology builds into an avoidable OOM. */
	maxActiveWorkingSetBytes = lod_default_working_set_limit();
    }

    BObolLodService *owner;
    mutable std::mutex mutex;
    std::condition_variable workerCv;
    std::condition_variable cacheWriterCv;
    std::condition_variable subscriberCv;
    std::vector<std::thread> workers;
    std::thread cacheWriter;
    SbBool running;
    SbBool stopping;
    SbBool cacheWriterStopping;
    SbBool cacheWriterEnabled;
    uint64_t nextTaskId;
    BObolLodSubscriberId nextSubscriberId;
    uint64_t nextGeneration;
    uint64_t activeGeneration;
    size_t maxActiveTasks;
    size_t maxQueuedResults;
    size_t maxQueuedCacheWrites;
    size_t maxActiveWorkingSetBytes;
    size_t activeWorkingSetBytes;
    size_t executingTasks;
    size_t peakWorkingSetBytes;
    size_t peakExecutingTasks;
    size_t resultReservations;
    size_t cacheWriteReservations;
    uint64_t rejectedTasks;
    uint64_t coalescedResults;
    uint64_t coalescedCacheWrites;
    uint64_t discardedStaleResults;
    uint64_t residentMeshCacheLoads;
    uint64_t residentMeshHits;
    uint64_t residentMeshCompactions;
    std::deque<BObolLodWorkItem> pending;
    std::deque<BObolLodResult> results;
    std::deque<BObolLodCacheWriteItem> cacheWrites;
    std::vector<BObolLodSubscriber> subscribers;
    std::unordered_map<std::string, size_t> activeRequestKeyCounts;
    std::unordered_map<std::string, std::shared_ptr<BObolResidentMeshAsset>>
	residentMeshes;
    std::unordered_map<uint64_t, BObolResidentMeshConsumerDemand>
	residentMeshConsumerDemands;
    std::set<uint64_t> completed;
    std::unordered_map<uint64_t, uint64_t> taskGenerations;
    std::set<uint64_t> cancelledGenerations;
    std::deque<uint64_t> cancelledGenerationOrder;
    std::unordered_map<uint64_t, size_t> generationTaskCounts;
    size_t inFlight;
    size_t cacheWriteInFlight;
    size_t delayedTasks;
};

static SbBool
lod_generation_cancelled_unlocked(const BObolLodServicePrivate *p,
				  uint64_t generation)
{
    return p->cancelledGenerations.find(generation) !=
	   p->cancelledGenerations.end() ? TRUE : FALSE;
}

static SbBool
lod_generation_cancelled(BObolLodServicePrivate *p, uint64_t generation)
{
    std::lock_guard<std::mutex> lock(p->mutex);
    return lod_generation_cancelled_unlocked(p, generation);
}

static SbBool
lod_generation_cancelled_or_stopping(BObolLodServicePrivate *p,
				     uint64_t generation)
{
    std::lock_guard<std::mutex> lock(p->mutex);
    return p->stopping ||
	   lod_generation_cancelled_unlocked(p, generation) ? TRUE : FALSE;
}

static void
lod_prune_cancelled_generations_unlocked(BObolLodServicePrivate *p)
{
    static const size_t maxHistory = 1024;
    if (!p || p->cancelledGenerationOrder.size() <= maxHistory)
	return;
    for (auto it = p->cancelledGenerationOrder.begin();
	it != p->cancelledGenerationOrder.end() &&
	p->cancelledGenerationOrder.size() > maxHistory;) {
	if (p->generationTaskCounts.find(*it) !=
	    p->generationTaskCounts.end()) {
	    ++it;
	    continue;
	}
	p->cancelledGenerations.erase(*it);
	it = p->cancelledGenerationOrder.erase(it);
    }
}

static void
lod_generation_task_finished_unlocked(BObolLodServicePrivate *p,
	uint64_t generation)
{
    auto found = p->generationTaskCounts.find(generation);
    if (found != p->generationTaskCounts.end()) {
	if (found->second > 1)
	    found->second--;
	else
	    p->generationTaskCounts.erase(found);
    }
    lod_prune_cancelled_generations_unlocked(p);
}

static SbBool
lod_request_has_identity(const BObolLodRequest &request)
{
    return request.databaseId.getLength() > 0 ||
	   request.objectPath.getLength() > 0 ||
	   request.objectName.getLength() > 0 ||
	   request.sourceRevision != 0 ||
	   request.sourceContentHash != 0 ? TRUE : FALSE;
}

static SbString
lod_request_active_key(const BObolLodRequest &request)
{
    /* Camera and policy epochs identify demand observations, not geometry.
     * Coalesce an unchanged source/level demand across those epochs so view
     * motion cannot enqueue the same work repeatedly.  The stable asset key
     * deliberately excludes level, so add the requested draw cut here.
     * Compact
     * occurrences remain distinct consumers: until the service has explicit
     * result fan-out, coalescing siblings here would publish the shared mesh
     * to only one occurrence and strand the others at boxes. */
    SbString key = bobol_lod_geometry_cache_key(request).value;
    key += "|level=";
    key += SbString(std::to_string(request.requestedLevel).c_str());
    if (request.occurrenceKey.getLength() > 0) {
	key += "|occurrence=";
	key += request.occurrenceKey;
    }
    return key;
}

static SbBool
lod_active_request_key_recorded_unlocked(const BObolLodServicePrivate *p,
	const SbString &key)
{
    if (!p || key.getLength() == 0)
	return FALSE;

    return p->activeRequestKeyCounts.find(key.getString()) !=
	   p->activeRequestKeyCounts.end() ? TRUE : FALSE;
}

static void
lod_active_request_key_remove_unlocked(BObolLodServicePrivate *p,
				       const SbString &key)
{
    if (!p || key.getLength() == 0)
	return;

    auto found = p->activeRequestKeyCounts.find(key.getString());
    if (found == p->activeRequestKeyCounts.end())
	return;
    if (found->second > 1)
	found->second--;
    else
	p->activeRequestKeyCounts.erase(found);
}

static BObolLodResult
lod_service_status_result(const BObolLodTask &task, int status,
			  const char *diagnostic)
{
    BObolLodResult result;

    result.generation = task.generation;
    result.request = task.request;
    result.cacheKey = bobol_lod_cache_key(task.request);
    result.qualityTier = task.request.qualityTier;
    result.providerStatus = status;
    result.terminal = TRUE;
    result.diagnostic = diagnostic ? diagnostic : "";
    if (status == BOBOL_LOD_PROVIDER_CANCELLED ||
	status == BOBOL_LOD_PROVIDER_STALE)
	result.stale = TRUE;

    return result;
}

static void
lod_delayed_task_count_add(BObolLodServicePrivate *p, int delta)
{
    std::lock_guard<std::mutex> lock(p->mutex);

    if (delta > 0) {
	p->delayedTasks += (size_t)delta;
	return;
    }

    size_t decrement = (size_t)(-delta);
    p->delayedTasks = decrement > p->delayedTasks ?
		      0 : p->delayedTasks - decrement;
}

static SbBool
lod_wait_for_debug_delay(BObolLodServicePrivate *p,
			 const BObolLodTask &task)
{
    if (task.debugDelayMilliseconds == 0)
	return TRUE;

    lod_delayed_task_count_add(p, 1);

    uint32_t remaining = task.debugDelayMilliseconds;
    while (remaining > 0) {
	if (lod_generation_cancelled_or_stopping(p, task.generation)) {
	    lod_delayed_task_count_add(p, -1);
	    return FALSE;
	}

	uint32_t slice = remaining > 10 ? 10 : remaining;
	std::this_thread::sleep_for(std::chrono::milliseconds(slice));
	remaining -= slice;
    }

    lod_delayed_task_count_add(p, -1);
    return !lod_generation_cancelled_or_stopping(p, task.generation);
}

static void
lod_normalize_result(BObolLodResult &result, const BObolLodTask &task)
{
    if (!result.cacheKey.isValid()) {
	if (lod_request_has_identity(result.request)) {
	    result.cacheKey = bobol_lod_cache_key(result.request);
	} else {
	    result.request = task.request;
	    result.cacheKey = bobol_lod_cache_key(task.request);
	}
    }

    if (!lod_request_has_identity(result.request))
	result.request = task.request;

    if (!bobol_lod_result_matches_request(result, task.request)) {
	result.providerStatus = BOBOL_LOD_PROVIDER_STALE;
	result.stale = TRUE;
	result.terminal = TRUE;
	if (result.diagnostic.getLength() == 0)
	    result.diagnostic = "LoD task returned stale request/result";
    }
}

static SbBool
lod_task_dependencies_ready(const BObolLodServicePrivate *p,
			    const BObolLodTask &task)
{
    if (lod_generation_cancelled_unlocked(p, task.generation))
	return TRUE;

    for (size_t i = 0; i < task.dependencies.size(); i++) {
	if (p->completed.find(task.dependencies[i]) == p->completed.end())
	    return FALSE;
    }

    return TRUE;
}

static size_t
lod_task_estimated_working_set_bytes(const BObolLodTask &task)
{
    if (task.estimatedWorkingSetBytes)
	return task.estimatedWorkingSetBytes;
    const BObolLodCounts &counts = task.request.sourceCounts;
    size_t estimate = static_cast<size_t>(std::min<uint64_t>(
	counts.byteCount, static_cast<uint64_t>(SIZE_MAX)));
    const auto addScaled = [&estimate](uint64_t count, size_t scale) {
	if (!count || estimate == SIZE_MAX)
	    return;
	if (count > static_cast<uint64_t>(SIZE_MAX / scale) ||
	    static_cast<size_t>(count) * scale > SIZE_MAX - estimate)
	    estimate = SIZE_MAX;
	else
	    estimate += static_cast<size_t>(count) * scale;
    };
    /* PoP construction temporarily owns source arrays, sorted topology
     * records, cumulative prefixes, and publication buffers.  These
     * deliberately conservative coefficients are scheduling reservations,
     * not resident-size accounting. */
    addScaled(counts.faceCount, 192);
    addScaled(counts.pointCount, 128);
    return estimate;
}

static SbBool
lod_task_working_set_available(const BObolLodServicePrivate *p,
	const BObolLodTask &task)
{
    if (!p || p->maxActiveWorkingSetBytes == SIZE_MAX)
	return TRUE;
    const size_t estimate = lod_task_estimated_working_set_bytes(task);
    if (!estimate)
	return TRUE;
    /* An exceptional task larger than the cap must still make progress, but
     * it runs alone. */
    if (p->activeWorkingSetBytes == 0)
	return TRUE;
    const size_t occupied = std::min(
	p->activeWorkingSetBytes, p->maxActiveWorkingSetBytes);
    return estimate <= p->maxActiveWorkingSetBytes - occupied ? TRUE : FALSE;
}

static std::deque<BObolLodWorkItem>::iterator
lod_find_ready_task(BObolLodServicePrivate *p)
{
    /* Prefer the coarsest ready task (lowest quality tier) so the cheap proxy /
     * bounding-box stages for the whole scene drain ahead of the expensive mesh
     * stages.  Plain FIFO selection picks an object's mesh task (which became
     * ready as soon as its own proxies finished) before a later object's proxy,
     * so bounding boxes trickle in interleaved with long mesh stalls.  Selecting
     * by tier instead yields a fast "all bounding boxes first, then refine to
     * meshes" frontier.  Ties keep FIFO (submission) order, preserving each
     * an object's explicitly declared dependency order.  The normal display
     * frontier is AABB -> view-selected PoP mesh; OBB remains an optional
     * provider capability rather than a mandatory intermediate stage. */
    std::deque<BObolLodWorkItem>::iterator best = p->pending.end();

    for (std::deque<BObolLodWorkItem>::iterator it = p->pending.begin();
	 it != p->pending.end(); ++it) {
	if (!lod_task_dependencies_ready(p, it->task))
	    continue;
	if (!lod_task_working_set_available(p, it->task))
	    continue;
	if (best == p->pending.end() ||
	    it->task.request.qualityTier < best->task.request.qualityTier)
	    best = it;
    }

    return best;
}

static BObolLodResult
lod_execute_task(BObolLodServicePrivate *p, const BObolLodTask &task)
{
    if (lod_generation_cancelled(p, task.generation))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled");

    if (!lod_wait_for_debug_delay(p, task))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled during debug delay");

    if (!task.realize)
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_ERROR,
					 "LoD task has no realization callback");

    BObolLodResult result = (*task.realize)(task.request, task.realizeData);

    if (lod_generation_cancelled(p, task.generation))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled");

    lod_normalize_result(result, task);
    return result;
}

static void
lod_task_free_realize_data(BObolLodTask &task)
{
    if (task.realizeDataFree && task.realizeData) {
	(*task.realizeDataFree)(task.realizeData);
	task.realizeData = NULL;
	task.realizeDataFree = NULL;
    }
}

struct BObolLodSubscriberCall {
    BObolLodSubscriberId id;
    BObolLodResultReadyCB callback;
    void *userData;
};

/* Result-ready callbacks run on a worker thread.  Track the entire collected
 * dispatch so a callback can also remove a later callback without waiting on
 * the reservation held by this same dispatch. */
static thread_local BObolLodServicePrivate *lod_callback_service = NULL;
static thread_local const std::vector<BObolLodSubscriberCall> *
    lod_callback_dispatch = NULL;

static size_t
lod_callback_dispatch_reservations(const BObolLodServicePrivate *p,
				   BObolLodSubscriberId id)
{
    if (lod_callback_service != p || !lod_callback_dispatch)
	return 0;

    size_t reservations = 0;
    for (size_t i = 0; i < lod_callback_dispatch->size(); i++) {
	if ((*lod_callback_dispatch)[i].id == id)
	    reservations++;
    }

    return reservations;
}

static std::vector<BObolLodSubscriberCall>
lod_collect_result_ready_callbacks(BObolLodServicePrivate *p)
{
    std::vector<BObolLodSubscriberCall> calls;
    std::lock_guard<std::mutex> lock(p->mutex);

    for (size_t i = 0; i < p->subscribers.size(); i++) {
	BObolLodSubscriber &subscriber = p->subscribers[i];
	if (!subscriber.active || !subscriber.callback)
	    continue;

	BObolLodSubscriberCall call;
	call.id = subscriber.id;
	call.callback = subscriber.callback;
	call.userData = subscriber.userData;
	subscriber.inFlight++;
	calls.push_back(call);
    }

    return calls;
}

static void
lod_complete_result_ready_callback(BObolLodServicePrivate *p,
					   BObolLodSubscriberId id)
{
    std::lock_guard<std::mutex> lock(p->mutex);

    for (size_t i = 0; i < p->subscribers.size(); i++) {
	if (p->subscribers[i].id != id)
	    continue;
	if (p->subscribers[i].inFlight > 0)
	    p->subscribers[i].inFlight--;
	break;
    }

    p->subscriberCv.notify_all();
}

static SbBool
lod_result_ready_callback_active(BObolLodServicePrivate *p,
				 BObolLodSubscriberId id)
{
    std::lock_guard<std::mutex> lock(p->mutex);

    for (size_t i = 0; i < p->subscribers.size(); i++) {
	if (p->subscribers[i].id == id)
	    return p->subscribers[i].active;
    }

    return FALSE;
}

static void
lod_notify_result_ready(BObolLodServicePrivate *p)
{
    std::vector<BObolLodSubscriberCall> calls =
	lod_collect_result_ready_callbacks(p);

    for (size_t i = 0; i < calls.size(); i++) {
	BObolLodServicePrivate *previous_service = lod_callback_service;
	const std::vector<BObolLodSubscriberCall> *previous_dispatch =
	    lod_callback_dispatch;
	lod_callback_service = p;
	lod_callback_dispatch = &calls;
	if (calls[i].callback &&
	    lod_result_ready_callback_active(p, calls[i].id))
	    (*calls[i].callback)(p->owner, calls[i].userData);
	lod_callback_service = previous_service;
	lod_callback_dispatch = previous_dispatch;
	lod_complete_result_ready_callback(p, calls[i].id);
    }
}

static std::string
lod_result_slot_key(const BObolLodResult &result)
{
    const BObolLodRequest &request = result.request;
    std::string key = request.databaseId.getString();
    key.push_back('\x1f');
    key += request.occurrenceKey.getLength() > 0 ?
	request.occurrenceKey.getString() :
	(request.objectPath.getLength() > 0 ?
	 request.objectPath.getString() : request.objectName.getString());
    key.push_back('\x1f');
    key += request.providerId.getString();
    char state[160] = {0};
    snprintf(state, sizeof(state),
	"|%d|%d|%d", request.drawMode, result.resultKind,
	result.resultKind == BOBOL_LOD_RESULT_PROXY ? result.proxy.kind : 0);
    key += state;
    return key;
}

static bool
lod_result_supersedes(const BObolLodResult &candidate,
	const BObolLodResult &current)
{
    if (candidate.request.databaseRevision != current.request.databaseRevision)
	return candidate.request.databaseRevision > current.request.databaseRevision;
    if (candidate.request.sourceRevision != current.request.sourceRevision)
	return candidate.request.sourceRevision > current.request.sourceRevision;
    if (candidate.request.viewRevision != current.request.viewRevision)
	return candidate.request.viewRevision > current.request.viewRevision;
    if (candidate.request.policyRevision != current.request.policyRevision)
	return candidate.request.policyRevision > current.request.policyRevision;
    if (candidate.qualityTier != current.qualityTier)
	return candidate.qualityTier > current.qualityTier;
    if (candidate.providerStatus == BOBOL_LOD_PROVIDER_READY &&
	current.providerStatus != BOBOL_LOD_PROVIDER_READY)
	return true;
    if (candidate.providerStatus != BOBOL_LOD_PROVIDER_READY &&
	current.providerStatus == BOBOL_LOD_PROVIDER_READY)
	return false;
    return true;
}

static void
lod_finish_task(BObolLodServicePrivate *p, const BObolLodWorkItem &item,
		BObolLodResult &&result)
{
    SbBool notifyResultReady = FALSE;
    BObolLodResult completedResult = std::move(result);
    completedResult.generation = item.task.generation;
    BObolLodResult cacheResult;
    const bool duplicateForCache = item.task.publishResult &&
	item.task.writeCache && item.task.cacheWrite;
    if (duplicateForCache)
	cacheResult = completedResult;

    {
	std::lock_guard<std::mutex> lock(p->mutex);

	const SbBool discardResult = p->stopping ||
	    lod_generation_cancelled_unlocked(p, item.task.generation);
	if (!discardResult)
	    p->completed.insert(item.id);
	else
	    p->taskGenerations.erase(item.id);
	if (p->inFlight > 0)
	    p->inFlight--;
	if (p->executingTasks > 0)
	    p->executingTasks--;
	const size_t workingBytes =
	    lod_task_estimated_working_set_bytes(item.task);
	p->activeWorkingSetBytes =
	    workingBytes >= p->activeWorkingSetBytes ?
	    0 : p->activeWorkingSetBytes - workingBytes;
	lod_active_request_key_remove_unlocked(p,
					       lod_request_active_key(item.task.request));

	if (item.task.publishResult) {
	    if (p->resultReservations > 0)
		p->resultReservations--;
	    if (discardResult) {
		p->discardedStaleResults++;
	    } else {
		const std::string slot = lod_result_slot_key(completedResult);
		auto existing = std::find_if(p->results.begin(), p->results.end(),
		    [&](const BObolLodResult &queued) {
			return queued.generation == completedResult.generation &&
			       lod_result_slot_key(queued) == slot;
		    });
		if (existing != p->results.end()) {
		    if (lod_result_supersedes(completedResult, *existing))
			*existing = std::move(completedResult);
		    p->coalescedResults++;
		} else {
		    notifyResultReady = p->results.empty() ? TRUE : FALSE;
		    p->results.push_back(std::move(completedResult));
		}
	    }
	}

	if (p->cacheWriterEnabled && item.task.writeCache &&
	    item.task.cacheWrite) {
	    if (p->cacheWriteReservations > 0)
		p->cacheWriteReservations--;
	    if (!discardResult) {
		BObolLodCacheWriteItem writeItem;
		writeItem.result = duplicateForCache ? std::move(cacheResult) :
		    std::move(completedResult);
		writeItem.write = item.task.cacheWrite;
		writeItem.writeData = item.task.cacheWriteData;
		const std::string slot = lod_result_slot_key(writeItem.result);
		auto existing = std::find_if(p->cacheWrites.begin(),
		    p->cacheWrites.end(),
		    [&](const BObolLodCacheWriteItem &queued) {
			return queued.result.generation ==
				writeItem.result.generation &&
			       lod_result_slot_key(queued.result) == slot;
		    });
		if (existing != p->cacheWrites.end()) {
		    if (lod_result_supersedes(writeItem.result,
			    existing->result))
			*existing = std::move(writeItem);
		    p->coalescedCacheWrites++;
		} else {
		    p->cacheWrites.push_back(std::move(writeItem));
		}
	    }
	}
	lod_generation_task_finished_unlocked(p, item.task.generation);
    }

    bobol_lod_working_set_release(
	lod_task_estimated_working_set_bytes(item.task));
    p->workerCv.notify_all();
    p->cacheWriterCv.notify_one();
    if (notifyResultReady)
	lod_notify_result_ready(p);
}

static void
lod_worker_loop(BObolLodServicePrivate *p)
{
    /* Display LoD is background work even when it is CPU intensive.  In
     * particular, a cold topology audit may sort hundreds of millions of
     * half-edges before a richer PoP cut is available.  Keep the host/UI
     * thread eligible to present the already-published coarse proxy.  On
     * platforms without per-thread nice support this is a harmless no-op. */
    bu_nice_set(5);

    for (;;) {
	BObolLodWorkItem item;

	{
	    std::unique_lock<std::mutex> lock(p->mutex);

	    while (!p->stopping && lod_find_ready_task(p) == p->pending.end())
		p->workerCv.wait(lock);

	    if (p->stopping)
		return;

	    std::deque<BObolLodWorkItem>::iterator ready =
		lod_find_ready_task(p);
	    if (ready == p->pending.end())
		continue;
	    item = *ready;
	    p->pending.erase(ready);
	    const size_t workingBytes =
		lod_task_estimated_working_set_bytes(item.task);
	    p->activeWorkingSetBytes =
		workingBytes > SIZE_MAX - p->activeWorkingSetBytes ?
		SIZE_MAX : p->activeWorkingSetBytes + workingBytes;
	    p->executingTasks++;
	    p->peakWorkingSetBytes = std::max(
		p->peakWorkingSetBytes, p->activeWorkingSetBytes);
	    p->peakExecutingTasks = std::max(
		p->peakExecutingTasks, p->executingTasks);
	}

	bobol_lod_working_set_acquire(
	    lod_task_estimated_working_set_bytes(item.task));
	BObolLodResult result = lod_execute_task(p, item.task);
	lod_finish_task(p, item, std::move(result));
	lod_task_free_realize_data(item.task);
    }
}

static void
lod_cache_writer_loop(BObolLodServicePrivate *p)
{
    /* Compression and cache persistence must not delay user input or frame
     * presentation either. */
    bu_nice_set(5);

    for (;;) {
	BObolLodCacheWriteItem item;

	{
	    std::unique_lock<std::mutex> lock(p->mutex);

	    while (!p->cacheWriterStopping && p->cacheWrites.empty())
		p->cacheWriterCv.wait(lock);

	    if (p->cacheWrites.empty() && p->cacheWriterStopping)
		return;

	    if (p->cacheWrites.empty())
		continue;

	    item = std::move(p->cacheWrites.front());
	    p->cacheWrites.pop_front();
	    p->cacheWriteInFlight++;
	}

	/* A cancellation can arrive after this item leaves the queue.  Do not
	 * persist a result that became stale while it was waiting for the cache
	 * writer.  A callback already in progress remains intentionally
	 * non-preemptible. */
	if (item.write && !lod_generation_cancelled_or_stopping(p,
		item.result.generation))
	    (*item.write)(item.result, item.writeData);

	{
	    std::lock_guard<std::mutex> lock(p->mutex);
	    if (p->cacheWriteInFlight > 0)
		p->cacheWriteInFlight--;
	}
	p->cacheWriterCv.notify_all();
    }
}

BObolLodService::BObolLodService(void) :
    p(new BObolLodServicePrivate(this))
{
}

BObolLodService::~BObolLodService(void)
{
    this->stop();
    delete this->p;
    this->p = NULL;
}

SbBool
BObolLodService::start(size_t workerCount, SbBool startCacheWriter)
{
    if (workerCount == 0)
	workerCount = 1;

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (this->p->running)
	    return TRUE;

	this->p->stopping = FALSE;
	this->p->cacheWriterStopping = FALSE;
	this->p->cacheWriterEnabled = startCacheWriter ? TRUE : FALSE;
	this->p->peakWorkingSetBytes = 0;
	this->p->peakExecutingTasks = 0;
	this->p->running = TRUE;
    }

    try {
	if (startCacheWriter)
	    this->p->cacheWriter =
		std::thread(lod_cache_writer_loop, this->p);

	for (size_t i = 0; i < workerCount; i++)
	    this->p->workers.push_back(std::thread(lod_worker_loop, this->p));
    } catch (...) {
	this->stop();
	return FALSE;
    }

    return TRUE;
}

void
BObolLodService::stop(void)
{
    std::vector<std::thread> workers;
    std::thread cacheWriter;
    std::deque<BObolLodWorkItem> pending;
    std::unordered_map<std::string,
	std::shared_ptr<BObolResidentMeshAsset>> residentMeshes;

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (!this->p->running && this->p->workers.empty() &&
	    !this->p->cacheWriter.joinable())
	    return;
	this->p->stopping = TRUE;
    }
    this->p->workerCv.notify_all();

    workers.swap(this->p->workers);
    for (size_t i = 0; i < workers.size(); i++) {
	if (workers[i].joinable())
	    workers[i].join();
    }

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	pending.swap(this->p->pending);
	this->p->inFlight = 0;
	this->p->activeWorkingSetBytes = 0;
	this->p->executingTasks = 0;
	this->p->resultReservations = 0;
	this->p->cacheWriteReservations = 0;
	this->p->activeRequestKeyCounts.clear();
	this->p->cacheWriterStopping = TRUE;
    }
    for (size_t i = 0; i < pending.size(); i++)
	lod_task_free_realize_data(pending[i].task);
    this->p->cacheWriterCv.notify_all();

    if (this->p->cacheWriter.joinable())
	cacheWriter.swap(this->p->cacheWriter);
    if (cacheWriter.joinable())
	cacheWriter.join();

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->cacheWrites.clear();
	this->p->results.clear();
	this->p->completed.clear();
	this->p->taskGenerations.clear();
	this->p->cancelledGenerations.clear();
	this->p->cancelledGenerationOrder.clear();
	this->p->generationTaskCounts.clear();
	this->p->residentMeshConsumerDemands.clear();
	residentMeshes.swap(this->p->residentMeshes);
	this->p->activeGeneration = 0;
	this->p->cacheWriteInFlight = 0;
	this->p->delayedTasks = 0;
	this->p->running = FALSE;
	this->p->stopping = FALSE;
	this->p->cacheWriterStopping = FALSE;
	this->p->cacheWriterEnabled = FALSE;
    }
    /* Destroy cache handles after dropping the service mutex. */
    residentMeshes.clear();
}

SbBool
BObolLodService::isRunning(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->running;
}

static size_t
lod_refinement_growth_budget(size_t current, uint64_t initial,
			     double growth)
{
    const size_t initialBudget = initial >
	static_cast<uint64_t>(std::numeric_limits<size_t>::max()) ?
	std::numeric_limits<size_t>::max() : static_cast<size_t>(initial);
    if (!current)
	return initialBudget;
    if (!std::isfinite(growth) || growth <= 1.0)
	growth = 1.0;
    const long double grown =
	static_cast<long double>(current) * static_cast<long double>(growth);
    const size_t growthBudget =
	grown >= static_cast<long double>(
	    std::numeric_limits<size_t>::max()) ?
	std::numeric_limits<size_t>::max() :
	static_cast<size_t>(std::ceil(grown));
    return std::max(initialBudget, growthBudget);
}

/* Select one bounded presentation step toward the screen-error target.  The
 * target itself remains request.requestedLevel; this helper only prevents a
 * first box->mesh replacement from allocating and uploading an arbitrarily
 * large cumulative prefix. */
static int
lod_progressive_delivery_level(
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int requestedLevel, int currentLevel,
    const BObolMeshLodProvider &provider)
{
    int target = std::max(hierarchy.min_level,
	std::min(hierarchy.max_level, requestedLevel));
    if (!provider.progressiveDelivery || provider.useForcedLevel ||
	provider.reset || currentLevel >= target)
	return target;

    const size_t currentFaces =
	currentLevel >= hierarchy.min_level ?
	hierarchy.face_count[currentLevel] : 0;
    const size_t currentPoints =
	currentLevel >= hierarchy.min_level ?
	hierarchy.point_count[currentLevel] : 0;
    const size_t faceBudget = lod_refinement_growth_budget(currentFaces,
	provider.initialRefinementFaceBudget,
	provider.refinementGrowthFactor);
    const size_t pointBudget = lod_refinement_growth_budget(currentPoints,
	provider.initialRefinementPointBudget,
	provider.refinementGrowthFactor);

    int selected = hierarchy.min_level;
    for (int level = hierarchy.min_level; level <= target; ++level) {
	if (hierarchy.face_count[level] > faceBudget ||
	    hierarchy.point_count[level] > pointBudget)
	    break;
	selected = level;
    }

    /* Sparse hierarchies may jump by more than the nominal growth factor.
     * Always make forward progress by one level; otherwise the controller
     * would keep requesting the same incomplete stage forever. */
    if (currentLevel >= hierarchy.min_level) {
	if (selected <= currentLevel)
	    selected = std::min(target, currentLevel + 1);
	/* Population budgets prevent an oversized initial replacement, but a
	 * later high-resolution target can otherwise skip several already
	 * defined PoP cuts in one presentation (Lucy level 9 -> 12 is an
	 * eighteen-million-face jump).  Publish at most one hierarchy step per
	 * rendered frame.  The retained prefix still appends in place and the
	 * controller immediately schedules the next step. */
	selected = std::min(selected, currentLevel + 1);
    }
    return std::min(selected, target);
}

BObolLodResult
BObolLodService::realizeResidentMeshLod(
    const BObolLodRequest &request,
    const BObolMeshLodProvider &provider)
{
    if (!provider.dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
	    "resident mesh provider has no database");
    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
	    "resident mesh provider request has no object name");

    const BObolLodCacheKey assetKey = bobol_lod_asset_cache_key(request);
    if (!assetKey.isValid())
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
	    "resident mesh provider could not form an asset key");

    std::shared_ptr<BObolResidentMeshAsset> resident;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	auto found = this->p->residentMeshes.find(assetKey.value.getString());
	if (found == this->p->residentMeshes.end()) {
	    resident = std::make_shared<BObolResidentMeshAsset>();
	    resident->dbip = provider.dbip;
	    resident->name = name;
	    resident->mesh = std::make_shared<BObolLodProgressiveMesh>();
	    this->p->residentMeshes[assetKey.value.getString()] = resident;
	} else {
	    resident = found->second;
	}
    }

    std::lock_guard<std::mutex> residentLock(resident->mutex);
    if (resident->dbip != provider.dbip || resident->name != name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_STALE,
	    "resident mesh asset identity collision");

    if (!resident->lod) {
	/* The normal warm path needs no separate status probe:
	 * bobol_mesh_lod_get already resolves the name key and validates the
	 * immutable hierarchy header.  Query/refresh only after that direct
	 * open misses.  At 50k distinct assets this removes two redundant LMDB
	 * transactions and one context lifecycle from every successful task. */
	resident->lod = bobol_mesh_lod_get(provider.dbip, name);
	if (resident->lod) {
	    resident->status.directory_found = 1;
	    resident->status.is_bot = 1;
	    resident->status.has_cache_key = 1;
	    resident->status.has_cached_payload = 1;
	    resident->status.stale_cache_entry = 0;
	    resident->status.cache_key =
		bobol_mesh_lod_cache_key_get(resident->lod);
	} else {
	    if (bobol_mesh_lod_cache_status(provider.dbip, name,
		    &resident->status) != BRLCAD_OK)
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider could not query cache status");
	    if ((!resident->status.has_cache_key ||
		 !resident->status.has_cached_payload ||
		 resident->status.stale_cache_entry) &&
		provider.refreshMissing) {
		const int64_t refreshStarted = bu_gettime();
		if (bobol_mesh_lod_cache_refresh(provider.dbip, name,
			&resident->status) != BRLCAD_OK)
		    return lod_provider_status_result(request,
			BOBOL_LOD_PROVIDER_CACHE_MISS,
			"resident mesh provider could not refresh cache entry");
		if (getenv("BOBOL_DRAW_TIMING"))
		    bu_log("[obol-timing] pop cache: refresh %-24s %8.1f ms\n",
			   name, (bu_gettime() - refreshStarted) / 1000.0);
	    }
	    resident->lod = bobol_mesh_lod_get(provider.dbip, name);
	}
	if (!resident->lod)
	    return lod_provider_status_result(request,
		resident->status.stale_cache_entry ?
		    BOBOL_LOD_PROVIDER_STALE :
		    BOBOL_LOD_PROVIDER_CACHE_MISS,
		"resident mesh provider has no cache payload");
    }

    int requestedLevel = provider.useForcedLevel ?
	provider.forcedLevel : request.requestedLevel;
    if (requestedLevel < 0 && provider.useView) {
	requestedLevel = bobol_mesh_lod_load_view(resident->lod,
	    &provider.view, 0);
	/* load_view is a legacy fallback.  Its arrays will be converted to the
	 * retained exact asset below. */
    }

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (!bobol_mesh_lod_hierarchy_info_get(resident->lod, &hierarchy))
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_CACHE_MISS,
	    "resident mesh provider loaded no hierarchy metadata");

    int currentLevel = bobol_mesh_lod_current_level(resident->lod);
    const int deliveryLevel = provider.useCurrentDrawLevel ?
	provider.currentDrawLevel : currentLevel;
    int drawTarget = requestedLevel;
    if (requestedLevel >= 0)
	drawTarget = lod_progressive_delivery_level(hierarchy, requestedLevel,
	    deliveryLevel, provider);
    int loadTarget = drawTarget;
    if (provider.compactResident && currentLevel >= 0 &&
	drawTarget >= 0 && drawTarget < currentLevel) {
	/* A level is not a bounded unit of memory: one Lucy hierarchy step can
	 * add millions of faces.  Stable reclamation therefore retains exactly
	 * the pixel-demanded prefix. */
	loadTarget = drawTarget;
    } else if (currentLevel >= 0 && drawTarget <= currentLevel &&
	!provider.reset) {
	loadTarget = currentLevel;
    }

    const bool loadNeeded = currentLevel < 0 || loadTarget != currentLevel ||
	provider.reset;
    int residentLevel = currentLevel;
    if (loadNeeded) {
	residentLevel = bobol_mesh_lod_load_resident_level(
	    resident->lod, loadTarget, provider.reset);
	if (residentLevel < 0)
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_CACHE_MISS,
		"resident mesh provider could not load the requested prefix");
	if (currentLevel < residentLevel) {
	    std::lock_guard<std::mutex> lock(this->p->mutex);
	    this->p->residentMeshCacheLoads++;
	}
    } else {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->residentMeshHits++;
    }

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    struct BObolMeshLodData data;
    if (!bobol_mesh_lod_info_get(resident->lod, &info) ||
	!bobol_mesh_lod_data_get(resident->lod, &data))
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_CACHE_MISS,
	    "resident mesh provider loaded no mesh data");

    const uint64_t priorRevision = resident->mesh->revision();
    if (priorRevision == 0 ||
	resident->mesh->residentLevel() != residentLevel) {
	if (!resident->mesh->update(data, hierarchy, residentLevel,
		info.shaded_cull_backfaces ? TRUE : FALSE))
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_ERROR,
		"resident mesh provider could not publish the retained asset");
    }

    int drawLevel = drawTarget;
    if (drawLevel < hierarchy.min_level)
	drawLevel = hierarchy.min_level;
    if (drawLevel > residentLevel)
	drawLevel = residentLevel;

    BObolLodResult result =
	bobol_lod_result_from_mesh_lod_info(request, info, &resident->status);
    result.geometry.cacheKey = assetKey;
    result.geometry.activeLevel = drawLevel;
    result.progressiveMesh = resident->mesh;
    result.residentLevel = residentLevel;
    result.counts.faceCount = resident->mesh->faceCount(drawLevel);
    result.counts.pointCount = resident->mesh->pointCount(drawLevel);
    result.counts.originalPointCount = result.counts.pointCount;
    result.counts.normalCount = info.has_normals ?
	result.counts.faceCount * 3 : 0;
    result.bounds = resident->mesh->bounds();
    result.hasSnappedPoints = FALSE;
    result.hasNormals = info.has_normals ? TRUE : FALSE;
    result.shadedCullBackfaces = resident->mesh->cullBackfaces();
    result.terminal = drawLevel >= std::max(hierarchy.min_level,
	std::min(hierarchy.max_level, requestedLevel)) ? TRUE : FALSE;

    /* Compact CAD occurrences retain the shared asset directly.  Legacy
     * shape consumers still receive a level-local by-value payload until
     * they are migrated to the same retained assembly path. */
    if (request.occurrenceKey.getLength() == 0 &&
	!resident->mesh->copyLevel(result.mesh, drawLevel))
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_ERROR,
	    "resident mesh provider could not materialize legacy payload");

    const char *traceFilter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (traceFilter && traceFilter[0] &&
	(strstr(name, traceFilter) ||
	 (request.objectPath.getLength() > 0 &&
	  strstr(request.objectPath.getString(), traceFilter)))) {
	bu_log("BObol resident LoD trace object=%s request_level=%d "
	       "draw_level=%d resident_level=%d faces=%zu points=%zu "
	       "asset_revision=%llu load=%d compact=%d terminal=%d\n",
	       name, request.requestedLevel, drawLevel, residentLevel,
	       result.counts.faceCount, result.counts.pointCount,
	       static_cast<unsigned long long>(resident->mesh->revision()),
	       loadNeeded ? 1 : 0, provider.compactResident ? 1 : 0,
	       result.terminal ? 1 : 0);
	if (getenv("BOBOL_LOD_TRACE_HIERARCHY")) {
	    for (int level = hierarchy.min_level;
		 level <= hierarchy.max_level; ++level)
		bu_log("BObol resident hierarchy object=%s level=%d "
		       "faces=%zu points=%zu%s\n",
		       name, level, hierarchy.face_count[level],
		       hierarchy.point_count[level],
		       level == hierarchy.max_level ? " terminal" : "");
	}
    }
    return result;
}

size_t
BObolLodService::compactResidentMeshes(
    uint64_t consumerId,
    uint64_t demandRevision,
    const std::vector<BObolLodResidentDemand> &demands)
{
    if (!consumerId)
	return 0;

    BObolResidentMeshConsumerDemand snapshot;
    snapshot.revision = demandRevision;
    for (const BObolLodResidentDemand &demand : demands) {
	if (demand.assetKey.getLength() == 0 || demand.level < 0)
	    continue;
	int &level = snapshot.levels[demand.assetKey.getString()];
	level = std::max(level, demand.level);
    }

    std::unordered_map<std::string, int> aggregate;
    std::vector<std::pair<std::shared_ptr<BObolResidentMeshAsset>, int>>
	work;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	BObolResidentMeshConsumerDemand &current =
	    this->p->residentMeshConsumerDemands[consumerId];
	if (demandRevision < current.revision)
	    return 0;
	current = std::move(snapshot);

	/* Never reclaim a buffer while a worker may still be appending to it.
	 * The committed snapshot remains recorded and the stable controller will
	 * retry its compaction pump once the service is idle. */
	if (!this->p->pending.empty() || this->p->inFlight != 0)
	    return 0;

	for (const auto &consumer : this->p->residentMeshConsumerDemands) {
	    for (const auto &demand : consumer.second.levels) {
		int &level = aggregate[demand.first];
		level = std::max(level, demand.second);
	    }
	}
	for (const auto &demand : aggregate) {
	    const auto resident = this->p->residentMeshes.find(demand.first);
	    if (resident != this->p->residentMeshes.end())
		work.push_back(std::make_pair(resident->second, demand.second));
	}
    }

    size_t compacted = 0;
    for (const auto &item : work) {
	const std::shared_ptr<BObolResidentMeshAsset> &resident = item.first;
	std::lock_guard<std::mutex> residentLock(resident->mutex);
	if (!resident->lod || !resident->mesh)
	    continue;
	struct BObolMeshLodHierarchyInfo hierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	if (!bobol_mesh_lod_hierarchy_info_get(resident->lod, &hierarchy))
	    continue;
	const int currentLevel = bobol_mesh_lod_current_level(resident->lod);
	const int targetLevel = std::min(currentLevel,
	    std::max(hierarchy.min_level,
		std::min(hierarchy.max_level,
		    item.second)));
	if (currentLevel < 0 || targetLevel >= currentLevel)
	    continue;
	const int residentLevel = bobol_mesh_lod_load_resident_level(
	    resident->lod, targetLevel, 0);
	if (residentLevel != targetLevel)
	    continue;
	struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
	struct BObolMeshLodData data;
	if (!bobol_mesh_lod_info_get(resident->lod, &info) ||
	    !bobol_mesh_lod_hierarchy_info_get(resident->lod, &hierarchy) ||
	    !bobol_mesh_lod_data_get(resident->lod, &data) ||
	    !resident->mesh->update(data, hierarchy, residentLevel,
		info.shaded_cull_backfaces ? TRUE : FALSE))
	    continue;
	compacted++;
    }
    if (compacted) {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->residentMeshCompactions += compacted;
    }
    return compacted;
}

void
BObolLodService::releaseResidentMeshConsumer(uint64_t consumerId)
{
    if (!consumerId)
	return;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    this->p->residentMeshConsumerDemands.erase(consumerId);
}

uint64_t
BObolLodService::beginGeneration(void)
{
    std::lock_guard<std::mutex> lock(this->p->mutex);

    this->p->nextGeneration++;
    if (this->p->nextGeneration == 0)
	this->p->nextGeneration++;
    this->p->activeGeneration = this->p->nextGeneration;
    return this->p->activeGeneration;
}

uint64_t
BObolLodService::currentGeneration(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->activeGeneration;
}

void
BObolLodService::cancelGeneration(uint64_t generation)
{
    if (generation == 0)
	return;

    std::deque<BObolLodWorkItem> cancelled;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	if (this->p->cancelledGenerations.insert(generation).second)
	    this->p->cancelledGenerationOrder.push_back(generation);
	if (this->p->activeGeneration == generation)
	    this->p->activeGeneration = 0;

	for (std::deque<BObolLodWorkItem>::iterator it =
		 this->p->pending.begin(); it != this->p->pending.end();) {
	    if (it->task.generation != generation) {
		++it;
		continue;
	    }
	    this->p->completed.insert(it->id);
	    if (this->p->inFlight > 0)
		this->p->inFlight--;
	    lod_generation_task_finished_unlocked(this->p, generation);
	    if (it->task.publishResult && this->p->resultReservations > 0)
		this->p->resultReservations--;
	    if (this->p->cacheWriterEnabled && it->task.writeCache &&
		it->task.cacheWrite && this->p->cacheWriteReservations > 0)
		this->p->cacheWriteReservations--;
	    lod_active_request_key_remove_unlocked(this->p,
		lod_request_active_key(it->task.request));
	    cancelled.push_back(*it);
	    it = this->p->pending.erase(it);
	}
	lod_prune_cancelled_generations_unlocked(this->p);

	for (std::deque<BObolLodResult>::iterator it =
		 this->p->results.begin(); it != this->p->results.end();) {
	    if (it->generation == generation)
		it = this->p->results.erase(it);
	    else
		++it;
	}
	for (std::deque<BObolLodCacheWriteItem>::iterator it =
		 this->p->cacheWrites.begin(); it != this->p->cacheWrites.end();) {
	    if (it->result.generation == generation)
		it = this->p->cacheWrites.erase(it);
	    else
		++it;
	}
	for (auto it = this->p->taskGenerations.begin();
	     it != this->p->taskGenerations.end();) {
	    if (it->second != generation) {
		++it;
		continue;
	    }
	    this->p->completed.erase(it->first);
	    it = this->p->taskGenerations.erase(it);
	}
    }
    for (size_t i = 0; i < cancelled.size(); i++)
	lod_task_free_realize_data(cancelled[i].task);
    this->p->workerCv.notify_all();
    this->p->cacheWriterCv.notify_all();
}

SbBool
BObolLodService::isGenerationCancelled(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_cancelled_unlocked(this->p, generation);
}

void
BObolLodService::setQueueLimits(size_t maxActiveTasks,
	size_t maxQueuedResults, size_t maxQueuedCacheWrites)
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    this->p->maxActiveTasks = maxActiveTasks > 0 ? maxActiveTasks : 1;
    this->p->maxQueuedResults = maxQueuedResults > 0 ? maxQueuedResults : 1;
    this->p->maxQueuedCacheWrites =
	maxQueuedCacheWrites > 0 ? maxQueuedCacheWrites : 1;
}

void
BObolLodService::getQueueLimits(size_t &maxActiveTasks,
	size_t &maxQueuedResults, size_t &maxQueuedCacheWrites) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    maxActiveTasks = this->p->maxActiveTasks;
    maxQueuedResults = this->p->maxQueuedResults;
    maxQueuedCacheWrites = this->p->maxQueuedCacheWrites;
}

void
BObolLodService::setWorkingSetLimit(size_t maxActiveBytes)
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (maxActiveBytes > 0) {
	this->p->maxActiveWorkingSetBytes = maxActiveBytes;
    } else {
	this->p->maxActiveWorkingSetBytes =
	    lod_default_working_set_limit();
    }
    this->p->workerCv.notify_all();
}

size_t
BObolLodService::getWorkingSetLimit(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->maxActiveWorkingSetBytes;
}

size_t
BObolLodService::activeWorkingSetBytesForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->activeWorkingSetBytes;
}

size_t
BObolLodService::executingTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->executingTasks;
}

size_t
BObolLodService::peakWorkingSetBytesForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->peakWorkingSetBytes;
}

size_t
BObolLodService::peakExecutingTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->peakExecutingTasks;
}

size_t
BObolLodService::availableResultTaskCapacity(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const size_t active = this->p->inFlight < this->p->maxActiveTasks ?
	this->p->maxActiveTasks - this->p->inFlight : 0;
    const size_t reserved = this->p->results.size() +
	this->p->resultReservations;
    const size_t results = reserved < this->p->maxQueuedResults ?
	this->p->maxQueuedResults - reserved : 0;
    return std::min(active, results);
}

static uint64_t
lod_service_submit_task(BObolLodServicePrivate *p,
			const BObolLodTask &task,
			SbBool skipActiveDuplicate)
{
    uint64_t id = 0;

    {
	std::lock_guard<std::mutex> lock(p->mutex);
	if (!p->running || p->stopping)
	    return 0;
	if ((task.generation != 0 &&
	    lod_generation_cancelled_unlocked(p, task.generation)) ||
	    p->inFlight >= p->maxActiveTasks ||
	    (task.publishResult &&
	     p->results.size() + p->resultReservations >=
		p->maxQueuedResults) ||
	    (p->cacheWriterEnabled && task.writeCache && task.cacheWrite &&
	     p->cacheWrites.size() + p->cacheWriteInFlight +
		p->cacheWriteReservations >= p->maxQueuedCacheWrites)) {
	    p->rejectedTasks++;
	    return 0;
	}

	BObolLodWorkItem item;
	item.id = p->nextTaskId++;
	if (p->nextTaskId == 0)
	    p->nextTaskId++;
	item.task = task;
	if (item.task.generation == 0) {
	    if (p->activeGeneration == 0) {
		p->nextGeneration++;
		p->activeGeneration = p->nextGeneration;
	    }
	    item.task.generation = p->activeGeneration;
	}

	const SbString activeKey =
	    lod_request_active_key(item.task.request);
	if (skipActiveDuplicate &&
	    lod_active_request_key_recorded_unlocked(p, activeKey))
	    return 0;

	id = item.id;
	p->pending.push_back(item);
	p->activeRequestKeyCounts[activeKey.getString()]++;
	p->generationTaskCounts[item.task.generation]++;
	p->taskGenerations[item.id] = item.task.generation;
	if (item.task.publishResult)
	    p->resultReservations++;
	if (p->cacheWriterEnabled && item.task.writeCache &&
	    item.task.cacheWrite)
	    p->cacheWriteReservations++;
	p->inFlight++;
    }

    p->workerCv.notify_all();
    return id;
}

uint64_t
BObolLodService::submit(const BObolLodTask &task)
{
    return lod_service_submit_task(this->p, task, FALSE);
}

uint64_t
BObolLodService::submitIfNotActive(const BObolLodTask &task)
{
    return lod_service_submit_task(this->p, task, TRUE);
}

SbBool
BObolLodService::hasActiveRequest(
    const BObolLodRequest &request) const
{
    const SbString key = lod_request_active_key(request);

    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_active_request_key_recorded_unlocked(this->p, key);
}

size_t
BObolLodService::drainResults(std::vector<BObolLodResult> &results,
				size_t maxResults)
{
    size_t count = 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);

    while (!this->p->results.empty() &&
	   (maxResults == 0 || count < maxResults)) {
	results.push_back(std::move(this->p->results.front()));
	this->p->results.pop_front();
	count++;
    }

    return count;
}

size_t
BObolLodService::drainMatchingResults(
    std::vector<BObolLodResult> &results,
    const std::vector<BObolLodRequest> &requests,
    size_t maxResults)
{
    if (requests.empty())
	return 0;

    size_t count = 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);

    for (std::deque<BObolLodResult>::iterator it =
	     this->p->results.begin(); it != this->p->results.end();) {
	if (maxResults != 0 && count >= maxResults)
	    break;

	SbBool matched = FALSE;
	for (size_t i = 0; i < requests.size(); i++) {
	    if (bobol_lod_result_matches_request(*it, requests[i])) {
		matched = TRUE;
		break;
	    }
	}
	if (!matched) {
	    ++it;
	    continue;
	}

	results.push_back(std::move(*it));
	it = this->p->results.erase(it);
	count++;
    }

    return count;
}

BObolLodSubscriberId
BObolLodService::subscribeResultReady(BObolLodResultReadyCB callback,
					void *userData)
{
    if (!callback)
	return 0;

    std::lock_guard<std::mutex> lock(this->p->mutex);

    BObolLodSubscriber subscriber;
    subscriber.id = this->p->nextSubscriberId++;
    if (this->p->nextSubscriberId == 0)
	this->p->nextSubscriberId++;
    subscriber.callback = callback;
    subscriber.userData = userData;
    subscriber.active = TRUE;
    this->p->subscribers.push_back(subscriber);

    return subscriber.id;
}

void
BObolLodService::unsubscribeResultReady(BObolLodSubscriberId id)
{
    if (id == 0)
	return;

    std::unique_lock<std::mutex> lock(this->p->mutex);

    for (size_t i = 0; i < this->p->subscribers.size(); i++) {
	if (this->p->subscribers[i].id != id)
	    continue;

	this->p->subscribers[i].active = FALSE;
	const size_t localReservations =
	    lod_callback_dispatch_reservations(this->p, id);
	this->p->subscriberCv.wait(lock, [this, id, localReservations] {
		for (size_t j = 0; j < this->p->subscribers.size(); j++) {
		    if (this->p->subscribers[j].id == id)
			return this->p->subscribers[j].inFlight <=
			    localReservations;
		}
		return true;
	});
	for (size_t j = 0; j < this->p->subscribers.size(); j++) {
	    if (this->p->subscribers[j].id == id) {
		this->p->subscribers.erase(
		    this->p->subscribers.begin() + (long)j);
		break;
	    }
	}
	return;
    }
}

size_t
BObolLodService::inFlightCount(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->inFlight;
}

size_t
BObolLodService::pendingTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->pending.size();
}

size_t
BObolLodService::queuedResultCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->results.size();
}

size_t
BObolLodService::queuedCacheWriteCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->cacheWrites.size() + this->p->cacheWriteInFlight;
}

size_t
BObolLodService::delayedTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->delayedTasks;
}

uint64_t
BObolLodService::rejectedTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->rejectedTasks;
}

uint64_t
BObolLodService::coalescedResultCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->coalescedResults;
}

uint64_t
BObolLodService::coalescedCacheWriteCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->coalescedCacheWrites;
}

uint64_t
BObolLodService::discardedStaleResultCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->discardedStaleResults;
}

size_t
BObolLodService::activeRequestCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->activeRequestKeyCounts.size();
}

size_t
BObolLodService::completedTaskCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->completed.size();
}

size_t
BObolLodService::cancelledGenerationCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->cancelledGenerations.size();
}

size_t
BObolLodService::residentMeshAssetCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshes.size();
}

size_t
BObolLodService::residentMeshBytesForDiagnostics(void) const
{
    std::vector<std::shared_ptr<BObolResidentMeshAsset>> assets;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	assets.reserve(this->p->residentMeshes.size());
	for (const auto &entry : this->p->residentMeshes)
	    assets.push_back(entry.second);
    }

    size_t bytes = 0;
    for (const std::shared_ptr<BObolResidentMeshAsset> &asset : assets) {
	if (!asset)
	    continue;
	/* Diagnostics run on the presentation/UI thread.  A cold cache refresh
	 * owns this mutex while it classifies and persists a very large mesh;
	 * waiting here turned a harmless telemetry sample into a multi-second
	 * input stall.  Stable samples see every idle asset exactly.  An in-flight
	 * sample is intentionally a non-blocking lower bound. */
	std::unique_lock<std::mutex> lock(asset->mutex, std::try_to_lock);
	if (!lock.owns_lock())
	    continue;
	if (!asset->mesh)
	    continue;
	const size_t assetBytes = asset->mesh->estimateBytes();
	if (assetBytes > std::numeric_limits<size_t>::max() - bytes)
	    return std::numeric_limits<size_t>::max();
	bytes += assetBytes;
    }
    return bytes;
}

uint64_t
BObolLodService::residentMeshCacheLoadCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshCacheLoads;
}

uint64_t
BObolLodService::residentMeshHitCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshHits;
}

uint64_t
BObolLodService::residentMeshCompactionCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshCompactions;
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
