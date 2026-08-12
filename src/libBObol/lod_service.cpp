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
#include "bu/datetime.h"

#include "BObol/BDrawCache.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshLodCache.h"

#include "database_source_realization.h"

#include "raytrace.h"
#include "rt/db_io.h"
#include "rt/view.h"

#include <algorithm>
#include <atomic>
#include <chrono>
#include <cmath>
#include <condition_variable>
#include <deque>
#include <iomanip>
#include <list>
#include <limits>
#include <map>
#include <mutex>
#include <new>
#include <set>
#include <sstream>
#include <string.h>
#include <thread>
#include <unordered_map>
#include <unordered_set>
#include <vector>

BObolDatabaseLease::BObolDatabaseLease(struct db_i *source) :
    database(source ? db_clone_dbi(source, NULL) : NULL)
{
}

BObolDatabaseLease::~BObolDatabaseLease(void)
{
    if (this->database)
	db_close_client(this->database, NULL);
    this->database = NULL;
}

std::shared_ptr<BObolDatabaseLease>
BObolDatabaseLease::acquire(struct db_i *database)
{
    if (!database)
	return std::shared_ptr<BObolDatabaseLease>();
    try {
	return std::shared_ptr<BObolDatabaseLease>(
	    new BObolDatabaseLease(database));
    } catch (const std::bad_alloc &) {
	return std::shared_ptr<BObolDatabaseLease>();
    }
}

struct db_i *
BObolDatabaseLease::get(void) const
{
    return this->database;
}

BObolMeshLodProvider::BObolMeshLodProvider(void)
{
    clear();
}

void
BObolMeshLodProvider::clear(void)
{
    service = NULL;
    databaseLease.reset();
    stagedSource.reset();
    meshAssetContentHash = 0;
    generateBrepVariant = FALSE;
    brepTessellationAbsTol = 0.0;
    brepTessellationRelTol = 0.0;
    brepTessellationNormTol = 0.0;
    brepVariantMemoryLimited = FALSE;
    refreshMissing = TRUE;
    useForcedCut = FALSE;
    shrinkAfterCopy = TRUE;
    compactResident = FALSE;
    progressiveDelivery = TRUE;
    initialRefinementFaceBudget = 250000;
    initialRefinementPointBudget = 500000;
    refinementGrowthFactor = 4.0;
    useCurrentDrawCut = FALSE;
    currentDrawCut = -1;
    useDeliveryCutLimit = FALSE;
    deliveryCutLimit = -1;
    usePresentationCutLimit = FALSE;
    presentationCutLimit = -1;
    prefetchCachedTargetOnFirstPublication = FALSE;
    forcedCut = 0;
    reset = 0;
}

SbBool
BObolMeshLodProvider::setDatabase(struct db_i *database)
{
    this->databaseLease = BObolDatabaseLease::acquire(database);
    return this->databaseLease ? TRUE : FALSE;
}

struct db_i *
BObolMeshLodProvider::getDatabase(void) const
{
    return this->databaseLease ? this->databaseLease->get() : NULL;
}

BObolRtSourceFullDetailProvider::BObolRtSourceFullDetailProvider(void)
{
    clear();
}

void
BObolRtSourceFullDetailProvider::clear(void)
{
    databaseLease.reset();
    validateSourceMetrics = TRUE;
    maxFullDetailFaceCount = 0;
    maxFullDetailPointCount = 0;
}

SbBool
BObolRtSourceFullDetailProvider::setDatabase(struct db_i *database)
{
    this->databaseLease = BObolDatabaseLease::acquire(database);
    return this->databaseLease ? TRUE : FALSE;
}

struct db_i *
BObolRtSourceFullDetailProvider::getDatabase(void) const
{
    return this->databaseLease ? this->databaseLease->get() : NULL;
}

BObolRtProxyProvider::BObolRtProxyProvider(void)
{
    clear();
}

void
BObolRtProxyProvider::clear(void)
{
    databaseLease.reset();
    proxyKind = BOBOL_LOD_PROXY_AABB;
    useRequestBounds = TRUE;
}

SbBool
BObolRtProxyProvider::setDatabase(struct db_i *database)
{
    this->databaseLease = BObolDatabaseLease::acquire(database);
    return this->databaseLease ? TRUE : FALSE;
}

struct db_i *
BObolRtProxyProvider::getDatabase(void) const
{
    return this->databaseLease ? this->databaseLease->get() : NULL;
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
    result.geometry.activeCut = -1;
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
    struct db_i *dbip = provider ? provider->getDatabase() : NULL;
    if (!dbip)
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

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "RT source full-detail provider could not find source object");

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    int internalType = rt_db_get_internal(&intern, dp, dbip, NULL);
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

    if (!provider->setDatabase(dbip)) {
	delete provider;
	return 0;
    }
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
    struct db_i *dbip = provider ? provider->getDatabase() : NULL;
    if (!dbip)
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
    if (bobol_mesh_lod_cache_status(dbip, name, &status) != BRLCAD_OK)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD provider could not query cache status");

    if ((!status.has_cache_key || !status.has_cached_payload ||
	 status.stale_cache_entry) && provider->refreshMissing) {
	if (bobol_mesh_lod_cache_refresh(dbip, name, &status) != BRLCAD_OK)
	    return lod_provider_status_result(request,
					      BOBOL_LOD_PROVIDER_CACHE_MISS,
					      "Obol mesh LoD provider could not refresh cache entry");
    }

    struct BObolMeshLod *lod = bobol_mesh_lod_get(dbip, name);
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

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (!bobol_mesh_lod_hierarchy_info_get(lod, &hierarchy)) {
	bobol_mesh_lod_destroy(lod);
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_CACHE_MISS,
	    "Obol mesh LoD provider loaded no hierarchy metadata");
    }
    int requestedCut = provider->useForcedCut ?
	provider->forcedCut : request.requestedCut;
    if (!provider->useForcedCut && request.projectedPixelDiameter > 0.0f &&
	request.targetPixelError > 0.0f)
	requestedCut = bobol_mesh_lod_select_cut(&hierarchy,
	    request.projectedPixelDiameter, request.targetPixelError);
    if (requestedCut < hierarchy.min_cut)
	requestedCut = hierarchy.min_cut;
    if (requestedCut > hierarchy.max_cut)
	requestedCut = hierarchy.max_cut;
    const int load_ret = bobol_mesh_lod_load_cut(
	lod, requestedCut, provider->reset);
    if (load_ret < 0) {
	bobol_mesh_lod_destroy(lod);
	return lod_provider_status_result(request,
					  BOBOL_LOD_PROVIDER_CACHE_MISS,
					  "Obol mesh LoD provider could not load the requested cut");
    }

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    int have_info = bobol_mesh_lod_info_get(lod, &info);
    const char *traceFilter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (traceFilter && traceFilter[0] &&
	((name && strstr(name, traceFilter)) ||
	 (request.objectPath.getLength() > 0 &&
	  strstr(request.objectPath.getString(), traceFilter)))) {
	bu_log("BObol LoD provider trace object=%s request_cut=%d "
	       "loaded_cut=%d faces=%zu points=%zu have_info=%d "
	       "view_revision=%llu policy_revision=%llu\n",
	       name ? name : "", request.requestedCut, load_ret,
	       info.face_count, info.point_count, have_info,
	       static_cast<unsigned long long>(request.viewRevision.value()),
	       static_cast<unsigned long long>(request.policyRevision.value()));
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
    result.resolvedCut = requestedCut;
    {
	BObolLodRequest resolvedIdentity = request;
	resolvedIdentity.requestedCut = requestedCut;
	result.geometry.cacheKey =
	    bobol_lod_geometry_cache_key(resolvedIdentity);
    }
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
    struct db_i *dbip = provider ? provider->getDatabase() : NULL;
    if (!dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD cache provider has no database");

    const char *name = lod_request_object_name(request);
    if (!name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD cache provider request has no object name");

    struct BObolMeshLodCacheStatus status =
	    BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    if (bobol_mesh_lod_cache_status(dbip, name, &status) != BRLCAD_OK)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol mesh LoD cache provider could not query cache status");

    if ((!status.has_cache_key || !status.has_cached_payload ||
	 status.stale_cache_entry) && provider->refreshMissing) {
	if (bobol_mesh_lod_cache_refresh(dbip, name, &status) != BRLCAD_OK)
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
	return lod_aabb_result_from_cache_or_db(request, provider->getDatabase(),
						provider->useRequestBounds);

    if (provider->proxyKind != BOBOL_LOD_PROXY_OBB)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
					  "Obol draw proxy provider has unknown proxy kind");

    BObolLodResult result = lod_obb_result_from_cache_or_db(request,
			      provider->getDatabase());
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

struct BObolLodResultSlotMapKey {
    std::string databaseId;
    std::string occurrence;
    std::string providerId;
    uint64_t generation = 0;
    uint64_t sourceRoutingId = 0;
    int drawMode = 0;
    int resultKind = 0;
    int proxyKind = 0;

    bool operator==(const BObolLodResultSlotMapKey &other) const
    {
	return generation == other.generation &&
	    sourceRoutingId == other.sourceRoutingId &&
	    drawMode == other.drawMode &&
	    resultKind == other.resultKind &&
	    proxyKind == other.proxyKind &&
	    databaseId == other.databaseId &&
	    occurrence == other.occurrence &&
	    providerId == other.providerId;
    }
};

struct BObolLodResultSlotMapKeyHash {
    size_t operator()(const BObolLodResultSlotMapKey &key) const
    {
	size_t hash = std::hash<std::string>()(key.databaseId);
	const auto combine = [&hash](size_t value) {
	    hash ^= value + static_cast<size_t>(0x9e3779b9U) +
		(hash << 6) + (hash >> 2);
	};
	combine(std::hash<std::string>()(key.occurrence));
	combine(std::hash<std::string>()(key.providerId));
	combine(std::hash<uint64_t>()(key.generation));
	combine(std::hash<uint64_t>()(key.sourceRoutingId));
	combine(std::hash<int>()(key.drawMode));
	combine(std::hash<int>()(key.resultKind));
	combine(std::hash<int>()(key.proxyKind));
	return hash;
    }
};

static BObolLodResultSlotMapKey
lod_result_slot_map_key(const BObolLodResult &result)
{
    const BObolLodRequest &request = result.request;
    BObolLodResultSlotMapKey key;
    key.databaseId = request.databaseId.getString();
    if (request.occurrenceKey.getLength() > 0)
	key.occurrence = request.occurrenceKey.getString();
    else if (request.objectPath.getLength() > 0)
	key.occurrence = request.objectPath.getString();
    else
	key.occurrence = request.objectName.getString();
    key.providerId = request.providerId.getString();
    key.generation = result.generation;
    key.sourceRoutingId = request.sourceRoutingId.value();
    key.drawMode = request.drawMode;
    key.resultKind = result.resultKind;
    key.proxyKind = result.resultKind == BOBOL_LOD_RESULT_PROXY ?
	result.proxy.kind : 0;
    return key;
}

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

static size_t
lod_default_resident_mesh_limit(void)
{
    const size_t mebibyte = 1024ULL * 1024ULL;
    const size_t gibibyte = 1024ULL * mebibyte;
    const size_t floor = 256ULL * mebibyte;
    const char *configured = getenv("BOBOL_LOD_RESIDENT_LIMIT_BYTES");
    if (configured && configured[0]) {
	char *end = NULL;
	const unsigned long long value = strtoull(configured, &end, 10);
	if (end && end != configured && *end == '\0') {
	    if (value > static_cast<unsigned long long>(SIZE_MAX))
		return SIZE_MAX;
	    return static_cast<size_t>(value);
	}
    }
    size_t totalBytes = 0;
    const bool haveTotal = bu_mem(BU_MEM_ALL, &totalBytes) >= 0 &&
	totalBytes > 0;

    /*
     * Retained CPU geometry coexists with the database, persistent-cache
     * mappings, transient topology work, GUI state, and a GPU-side copy.
     * Reserve most host memory for those consumers while still allowing a
     * realistic multi-gigabyte vehicle working set on a capable machine.
     */
    size_t allowance = 4 * gibibyte;
    if (haveTotal)
	allowance = std::min(allowance,
	    std::max(floor, totalBytes / 4));
    /* This is the durable resident ceiling, not the transient admission
     * governor above.  Base it on machine capacity so loading the database or
     * constructing cold PoP caches cannot permanently shrink every later
     * view.  Concurrent topology work remains bounded by the separately
     * available-memory-aware working-set governor, while the conservative
     * quarter-RAM share leaves the database, UI, cache mappings, and renderer
     * the majority of host memory. */
    return std::max(floor, allowance);
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
    BObolResidentMeshAsset(void) :
	publishedMinimumCut(-1),
	publishedResidentCut(-1),
	publishedBytes(0),
	publishedBackingPrefixBytes(0),
	useRevision(0),
	orderIndex(SIZE_MAX)
    {
    }

    ~BObolResidentMeshAsset()
    {
	if (lod)
	    bobol_mesh_lod_destroy(lod);
	lod = NULL;
    }

    std::mutex mutex;
    std::string databaseIdentity;
    std::string name;
    struct BObolMeshLod *lod = NULL;
    BObolLodProgressiveMeshPtr mesh;
    struct BObolMeshLodCacheStatus status =
	BOBOL_MESH_LOD_CACHE_STATUS_INIT;
    /* Planner-side summaries avoid retaining every immutable source
     * generation merely to decide whether a stable trim is necessary. */
    std::atomic<int> publishedMinimumCut;
    std::atomic<int> publishedResidentCut;
    /* Total service-owned CPU bytes for this asset: the renderer-neutral
     * immutable prefix plus the opened cache handle and any reloadable cache
     * prefix arrays. */
    std::atomic<size_t> publishedBytes;
    /* Reloadable arrays duplicated temporarily by the cache reader.  Stable
     * maintenance releases these after the immutable generation is
     * published; fixed hierarchy/header state remains in publishedBytes. */
    std::atomic<size_t> publishedBackingPrefixBytes;
    /* An eviction plan is valid only while no later realization has acquired
     * this asset.  This prevents a quiet-view reclamation queued just before
     * renewed input from retiring the asset underneath that new request. */
    std::atomic<uint64_t> useRevision;
    size_t orderIndex;
};

struct BObolResidentMeshDemandValue {
    int cut = -1;
    unsigned int channelMask = 0;
};

struct BObolResidentMeshConsumerDemand {
    uint64_t revision = 0;
    // Resident-asset mutation epoch for which this snapshot was fully
    // compacted.  Zero means the snapshot was recorded while workers were
    // active and must be retried.
    uint64_t residentRevision = 0;
    std::unordered_map<std::string, BObolResidentMeshDemandValue> assets;
    size_t planningCursor = 0;
    size_t planningProjectedResidentBytes = 0;
    uint64_t planningResidentRevision = 0;
    SbBool planning = FALSE;
};

struct BObolResidentMeshCompactionTarget {
    int cut = -1;
    unsigned int channelMask = 0;
    SbBool evict = FALSE;
    uint64_t useRevision = 0;
    uint64_t revision = 0;
};

struct BObolResidentMeshCompactionWork {
    std::string assetKey;
    std::shared_ptr<BObolResidentMeshAsset> resident;
    BObolResidentMeshCompactionTarget target;
    std::vector<uint64_t> consumers;
    size_t estimatedWorkingSetBytes = 0;
};

/*
 * BObolLodService concurrency contract
 * ------------------------------------
 *
 * Lock order:
 *   1. BObolLodServicePrivate::mutex (queue/generation/subscriber/service
 *      residency maps)
 *   2. BObolResidentMeshAsset::mutex (one retained asset)
 *
 * Code which needs both must acquire them in that order.  Expensive provider,
 * cache, mesh-preparation, Coin/presentation, and subscriber callbacks execute
 * with neither lock held.  The callback collector reserves subscribers under
 * the service lock, drops it, invokes callbacks, and reacquires it only to
 * release each reservation.
 *
 * Ownership:
 *   - pending/completed/cache queues and generation tables are pump/service
 *     mutex owned;
 *   - realization callbacks and task-local payloads are worker owned;
 *   - diagnostic byte/count summaries used without the service lock are
 *     atomic;
 *   - resident mesh arrays are guarded by the resident mutex and published as
 *     immutable shared objects;
 *   - Coin nodes and fields are never mutated by this service and remain
 *     presentation-owner-thread only.
 */
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
	/*
	 * Large compact sources are planned in 2048-occurrence quiet-view
	 * windows.  Keeping only 256 result reservations forced every such
	 * window through eight producer/publication barriers and, more
	 * importantly, eight whole-scene update traversals.  Pending tasks and
	 * result handles are lightweight; actual concurrent mesh construction
	 * remains independently bounded by worker count and the byte governor.
	 */
	maxActiveTasks(4096),
	maxQueuedResults(2048),
	maxQueuedCacheWrites(2048),
	maxActiveWorkingSetBytes(0),
	maxResidentMeshBytes(0),
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
	residentMeshEvictions(0),
	residentMeshBytes(0),
	residentMeshBackingBytes(0),
	residentMeshGrowthReservationBytes(0),
	residentMeshRevision(1),
	residentMeshAdmissionRevision(1),
	residentMeshCompactionsInFlight(0),
	residentMeshCompactionResultCount(0),
	residentMeshCompactionResultReservations(0),
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
	maxResidentMeshBytes = lod_default_resident_mesh_limit();
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
    size_t maxResidentMeshBytes;
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
    uint64_t residentMeshEvictions;
    /* Updated only when a retained progressive buffer is published or
     * compacted.  Diagnostics/HUD reads must not walk and try-lock every
     * resident asset on the presentation thread. */
    std::atomic<size_t> residentMeshBytes;
    /* Reloadable cache-reader prefix bytes are part of live diagnostics but
     * not the quiet-state residency target.  The transient working-set
     * governor already bounds their concurrent construction. */
    std::atomic<size_t> residentMeshBackingBytes;
    /* Protected by mutex.  Workers reserve optional stable-prefix growth
     * before loading so independent assets cannot all observe the same free
     * capacity.  Minimum useful prefixes are permitted to exceed the soft
     * target, but are still reserved and therefore constrain richer peers. */
    size_t residentMeshGrowthReservationBytes;
    std::atomic<uint64_t> residentMeshRevision;
    std::atomic<uint64_t> residentMeshAdmissionRevision;
    std::deque<BObolLodWorkItem> pending;
    std::map<int, size_t> pendingQualityCounts;
    std::list<BObolLodResult> results;
    std::list<BObolLodCacheWriteItem> cacheWrites;
    std::unordered_map<BObolLodResultSlotMapKey,
	std::list<BObolLodResult>::iterator,
	BObolLodResultSlotMapKeyHash> resultSlots;
    std::unordered_map<BObolLodResultSlotMapKey,
	std::list<BObolLodCacheWriteItem>::iterator,
	BObolLodResultSlotMapKeyHash> cacheWriteSlots;
    std::vector<BObolLodSubscriber> subscribers;
    std::unordered_map<std::string, size_t> activeRequestKeyCounts;
    /* A completed producer still owns its request identity until the queued
     * presentation result is drained.  This prevents a fast cache hit from
     * being resubmitted while the GUI intentionally coalesces result waves. */
    std::unordered_map<std::string, size_t> queuedResultRequestKeyCounts;
    std::unordered_map<std::string, std::shared_ptr<BObolResidentMeshAsset>>
	residentMeshes;
    /* Append-only while the service is running.  Stable planning advances a
     * bounded cursor through this vector; unordered_map rehashing therefore
     * cannot invalidate a GUI-frame-spanning plan. */
    std::vector<std::pair<std::string,
	std::shared_ptr<BObolResidentMeshAsset>>> residentMeshOrder;
    std::unordered_map<uint64_t, BObolResidentMeshConsumerDemand>
	residentMeshConsumerDemands;
    std::deque<BObolResidentMeshCompactionWork>
	residentMeshCompactionWork;
    std::unordered_set<std::string> residentMeshCompactionQueuedAssets;
    std::unordered_map<std::string, BObolResidentMeshCompactionTarget>
	residentMeshCompactionTargets;
    uint64_t nextResidentMeshCompactionTargetRevision = 1;
    std::unordered_map<uint64_t,
	std::deque<BObolLodResidentCompaction>>
	residentMeshCompactionResults;
    size_t residentMeshCompactionsInFlight;
    size_t residentMeshCompactionResultCount;
    size_t residentMeshCompactionResultReservations;
    std::set<uint64_t> completed;
    std::unordered_map<uint64_t, uint64_t> taskGenerations;
    std::set<uint64_t> cancelledGenerations;
    std::deque<uint64_t> cancelledGenerationOrder;
    std::unordered_map<uint64_t, size_t> generationTaskCounts;
    std::unordered_map<uint64_t, size_t> generationPendingTaskCounts;
    std::unordered_map<uint64_t, size_t> generationExecutingTaskCounts;
    std::unordered_map<uint64_t, size_t> generationDelayedTaskCounts;
    std::unordered_map<uint64_t, size_t> generationResultCounts;
    std::unordered_map<uint64_t, size_t> generationCacheWriteCounts;
    size_t inFlight;
    size_t cacheWriteInFlight;
    size_t delayedTasks;
};

static void
lod_resident_mesh_bytes_replace(std::atomic<size_t> &total,
    size_t priorBytes, size_t currentBytes)
{
    size_t observed = total.load(std::memory_order_relaxed);
    while (true) {
	size_t next = observed;
	if (priorBytes > next)
	    next = 0;
	else
	    next -= priorBytes;
	if (currentBytes > std::numeric_limits<size_t>::max() - next)
	    next = std::numeric_limits<size_t>::max();
	else
	    next += currentBytes;
	if (total.compare_exchange_weak(observed, next,
		std::memory_order_relaxed, std::memory_order_relaxed))
	    return;
    }
}

static void
lod_resident_mesh_revision_advance(std::atomic<uint64_t> &revision)
{
    uint64_t next = revision.fetch_add(1, std::memory_order_relaxed) + 1;
    if (!next)
	revision.store(1, std::memory_order_relaxed);
}

static size_t
lod_resident_stable_bytes(const BObolLodServicePrivate *p)
{
    if (!p)
	return 0;
    const size_t total =
	p->residentMeshBytes.load(std::memory_order_relaxed);
    const size_t backing =
	p->residentMeshBackingBytes.load(std::memory_order_relaxed);
    return backing >= total ? 0 : total - backing;
}

static size_t
lod_resident_cut_stable_bytes(
    const BObolResidentMeshAsset &resident,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int cut)
{
    if (!resident.lod || cut < hierarchy.min_cut ||
	cut > hierarchy.max_cut ||
	cut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return SIZE_MAX;

    const size_t cacheBytes =
	bobol_mesh_lod_resident_bytes(resident.lod);
    const size_t prefixBytes =
	bobol_mesh_lod_resident_prefix_bytes(resident.lod);
    size_t bytes = prefixBytes >= cacheBytes ?
	0 : cacheBytes - prefixBytes;
    const size_t points = hierarchy.cuts[cut].point_count;
    const size_t faces = hierarchy.cuts[cut].face_count;
    const auto addScaled = [&bytes](size_t count, size_t stride) {
	if (bytes == SIZE_MAX || (count && stride > SIZE_MAX / count)) {
	    bytes = SIZE_MAX;
	    return;
	}
	const size_t amount = count * stride;
	bytes = amount > SIZE_MAX - bytes ? SIZE_MAX : bytes + amount;
    };

    if (hierarchy.has_normals) {
	/* Authored corner normals may split every triangle corner into a
	 * distinct renderer vertex.  Count the larger of the source prefix
	 * and that worst-case split, then one 32-bit index per corner. */
	const size_t corners =
	    faces > SIZE_MAX / 3 ? SIZE_MAX : faces * 3;
	const size_t rendererVertices = std::max(points, corners);
	addScaled(rendererVertices, sizeof(SbVec3f) * 2);
	addScaled(corners, sizeof(uint32_t));
    } else {
	addScaled(points, sizeof(SbVec3f));
	addScaled(faces, sizeof(uint32_t) * 3);
    }
    return bytes;
}

class BObolResidentMeshGrowthReservation {
public:
    explicit BObolResidentMeshGrowthReservation(
	BObolLodServicePrivate *service) : p(service)
    {
    }

    ~BObolResidentMeshGrowthReservation()
    {
	release();
    }

    int admit(const BObolResidentMeshAsset &resident,
	const struct BObolMeshLodHierarchyInfo &hierarchy,
	int desiredCut, int publishedCut, size_t priorStableBytes,
	SbBool &limited)
    {
	limited = FALSE;
	if (!p || desiredCut < hierarchy.min_cut)
	    return desiredCut;
	desiredCut = std::min(desiredCut, hierarchy.max_cut);
	if (publishedCut >= desiredCut)
	    return desiredCut;

	const int floorCut = publishedCut >= hierarchy.min_cut ?
	    publishedCut : hierarchy.min_cut;
	int admittedCut = floorCut;
	size_t admittedEstimate =
	    publishedCut >= hierarchy.min_cut ?
		priorStableBytes :
		lod_resident_cut_stable_bytes(
		    resident, hierarchy, floorCut);

	std::lock_guard<std::mutex> lock(p->mutex);
	decisionRevision =
	    p->residentMeshAdmissionRevision.load(
		std::memory_order_relaxed);
	const size_t stableBytes = lod_resident_stable_bytes(p);
	const size_t occupied =
	    p->residentMeshGrowthReservationBytes >
		    SIZE_MAX - stableBytes ?
		SIZE_MAX :
		stableBytes + p->residentMeshGrowthReservationBytes;
	const size_t limit = p->maxResidentMeshBytes;
	for (int cut = floorCut + 1;
	    cut <= desiredCut; ++cut) {
	    const size_t estimate =
		lod_resident_cut_stable_bytes(
		    resident, hierarchy, cut);
	    const size_t growth = estimate > priorStableBytes ?
		estimate - priorStableBytes : 0;
	    const SbBool fits =
		limit == SIZE_MAX ||
		(occupied <= limit && growth <= limit - occupied);
	    if (!fits)
		break;
	    admittedCut = cut;
	    admittedEstimate = estimate;
	}

	/* A first useful mesh is a visual correctness floor.  Reserve it even
	 * when all visible minima collectively exceed the soft target; optional
	 * suffixes observe that overage and are denied. */
	const size_t growth = admittedEstimate > priorStableBytes ?
	    admittedEstimate - priorStableBytes : 0;
	p->residentMeshGrowthReservationBytes =
	    growth > SIZE_MAX - p->residentMeshGrowthReservationBytes ?
		SIZE_MAX :
		p->residentMeshGrowthReservationBytes + growth;
	bytes = growth;
	limited = admittedCut < desiredCut ? TRUE : FALSE;
	return admittedCut;
    }

    uint64_t revision(void) const
    {
	return decisionRevision;
    }

    void release(void)
    {
	if (!p || !bytes)
	    return;
	std::lock_guard<std::mutex> lock(p->mutex);
	p->residentMeshGrowthReservationBytes =
	    bytes >= p->residentMeshGrowthReservationBytes ?
		0 : p->residentMeshGrowthReservationBytes - bytes;
	bytes = 0;
	p->workerCv.notify_all();
    }

private:
    BObolLodServicePrivate *p = NULL;
    size_t bytes = 0;
    uint64_t decisionRevision = 0;
};

static size_t
lod_resident_asset_bytes(const BObolResidentMeshAsset &resident)
{
    const size_t meshBytes =
	resident.mesh ? resident.mesh->estimateBytes() : 0;
    const size_t backingBytes =
	bobol_mesh_lod_resident_bytes(resident.lod);
    return backingBytes > SIZE_MAX - meshBytes ?
	SIZE_MAX : meshBytes + backingBytes;
}

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

static void
lod_generation_count_add_unlocked(
    std::unordered_map<uint64_t, size_t> &counts, uint64_t generation)
{
    counts[generation]++;
}

static void
lod_generation_count_remove_unlocked(
    std::unordered_map<uint64_t, size_t> &counts, uint64_t generation)
{
    const auto found = counts.find(generation);
    if (found == counts.end())
	return;
    if (found->second > 1)
	found->second--;
    else
	counts.erase(found);
}

static size_t
lod_generation_count_unlocked(
    const std::unordered_map<uint64_t, size_t> &counts,
    uint64_t generation)
{
    const auto found = counts.find(generation);
    return found == counts.end() ? 0 : found->second;
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
     * Coalesce all in-flight levels for one occurrence.  Resident PoP growth
     * is serialized by asset, so queuing levels 10, 11, and 12 during a wheel
     * burst cannot execute useful work in parallel; it only loads obsolete
     * prefixes, publishes stale results, and consumes working-set slots.  A
     * completion wakes the bounded planner, which then submits the newest
     * demand if the resident high-water mark is still insufficient.
     *
     * Compact occurrences remain distinct consumers: until the service has
     * explicit result fan-out, coalescing siblings here would publish the
     * shared mesh to only one occurrence and strand the others at boxes. */
    SbString key = bobol_lod_geometry_cache_key(request).value;
    if (request.occurrenceKey.getLength() > 0) {
	key += "|occurrence=";
	key += request.occurrenceKey;
    }
    /*
     * Two live database sources may legitimately request the same database
     * occurrence.  Results are delivered back to a sourceRoutingId, so
     * treating those requests as one active task would strand one source.
     */
    if (request.sourceRoutingId != 0) {
	key += "|route=";
	key += SbString(std::to_string(request.sourceRoutingId.value()).c_str());
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

static SbBool
lod_request_key_recorded_unlocked(const BObolLodServicePrivate *p,
	const SbString &key)
{
    if (lod_active_request_key_recorded_unlocked(p, key))
	return TRUE;
    if (!p || key.getLength() == 0)
	return FALSE;
    return p->queuedResultRequestKeyCounts.find(key.getString()) !=
	p->queuedResultRequestKeyCounts.end() ? TRUE : FALSE;
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

static void
lod_queued_result_request_key_add_unlocked(BObolLodServicePrivate *p,
	const BObolLodRequest &request)
{
    if (!p)
	return;
    const SbString key = lod_request_active_key(request);
    if (key.getLength() > 0)
	p->queuedResultRequestKeyCounts[key.getString()]++;
}

static void
lod_queued_result_request_key_remove_unlocked(BObolLodServicePrivate *p,
	const BObolLodRequest &request)
{
    if (!p)
	return;
    const SbString key = lod_request_active_key(request);
    if (key.getLength() == 0)
	return;
    auto found = p->queuedResultRequestKeyCounts.find(key.getString());
    if (found == p->queuedResultRequestKeyCounts.end())
	return;
    if (found->second > 1)
	found->second--;
    else
	p->queuedResultRequestKeyCounts.erase(found);
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
lod_delayed_task_count_add(BObolLodServicePrivate *p,
	uint64_t generation, int delta)
{
    std::lock_guard<std::mutex> lock(p->mutex);

    if (delta > 0) {
	p->delayedTasks += (size_t)delta;
	lod_generation_count_add_unlocked(
	    p->generationDelayedTaskCounts, generation);
	return;
    }

    size_t decrement = (size_t)(-delta);
    p->delayedTasks = decrement > p->delayedTasks ?
		      0 : p->delayedTasks - decrement;
    lod_generation_count_remove_unlocked(
	p->generationDelayedTaskCounts, generation);
}

static SbBool
lod_wait_for_debug_delay(BObolLodServicePrivate *p,
			 const BObolLodTask &task)
{
    if (task.debugDelayMilliseconds == 0)
	return TRUE;

    lod_delayed_task_count_add(p, task.generation, 1);

    uint32_t remaining = task.debugDelayMilliseconds;
    while (remaining > 0) {
	if (lod_generation_cancelled_or_stopping(p, task.generation)) {
	    lod_delayed_task_count_add(p, task.generation, -1);
	    return FALSE;
	}

	uint32_t slice = remaining > 10 ? 10 : remaining;
	std::this_thread::sleep_for(std::chrono::milliseconds(slice));
	remaining -= slice;
    }

    lod_delayed_task_count_add(p, task.generation, -1);
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
    result.canonicalizePayload();
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

static void
lod_pending_quality_remove(BObolLodServicePrivate *p, int qualityTier)
{
    const auto found = p->pendingQualityCounts.find(qualityTier);
    if (found == p->pendingQualityCounts.end())
	return;
    if (found->second > 1)
	--found->second;
    else
	p->pendingQualityCounts.erase(found);
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

    /*
     * Normal compact waves contain one quality tier and no dependencies.
     * Return their FIFO head in constant time.  The ordered count retains the
     * original global tier preference; mixed/blocked queues take the complete
     * correctness path below.
     */
    if (!p->pending.empty() && !p->pendingQualityCounts.empty() &&
	p->pending.front().task.request.qualityTier ==
	    p->pendingQualityCounts.begin()->first &&
	lod_task_dependencies_ready(p, p->pending.front().task) &&
	lod_task_working_set_available(p, p->pending.front().task))
	return p->pending.begin();

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
    /* The cache writer serializes source PoP data, not renderer presentation
     * snapshots.  Do not let a slow disk queue retain a second strong
     * reference to large GPU-ready arrays. */
    if (duplicateForCache) {
	cacheResult.preparedCadGeometry.reset();
	cacheResult.preparedCadGeometryRevision = 0;
    }

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
	lod_generation_count_remove_unlocked(
	    p->generationExecutingTaskCounts, item.task.generation);
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
		const BObolLodResultSlotMapKey slot =
		    lod_result_slot_map_key(completedResult);
		const auto existing = p->resultSlots.find(slot);
		if (existing != p->resultSlots.end()) {
		    if (lod_result_supersedes(
			    completedResult, *existing->second)) {
			const SbString oldRequestKey = lod_request_active_key(
			    existing->second->request);
			const SbString newRequestKey = lod_request_active_key(
			    completedResult.request);
			if (oldRequestKey != newRequestKey) {
			    lod_queued_result_request_key_remove_unlocked(
				p, existing->second->request);
			    lod_queued_result_request_key_add_unlocked(
				p, completedResult.request);
			}
			*existing->second = std::move(completedResult);
		    }
		    p->coalescedResults++;
		} else {
		    notifyResultReady =
			lod_generation_count_unlocked(
			    p->generationResultCounts,
			    completedResult.generation) == 0 ? TRUE : FALSE;
		    lod_queued_result_request_key_add_unlocked(
			p, completedResult.request);
		    p->results.push_back(std::move(completedResult));
		    p->resultSlots.emplace(slot, std::prev(p->results.end()));
		    lod_generation_count_add_unlocked(
			p->generationResultCounts,
			p->results.back().generation);
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
		const BObolLodResultSlotMapKey slot =
		    lod_result_slot_map_key(writeItem.result);
		const auto existing = p->cacheWriteSlots.find(slot);
		if (existing != p->cacheWriteSlots.end()) {
		    if (lod_result_supersedes(writeItem.result,
			    existing->second->result))
			*existing->second = std::move(writeItem);
		    p->coalescedCacheWrites++;
		} else {
		    p->cacheWrites.push_back(std::move(writeItem));
		    p->cacheWriteSlots.emplace(
			slot, std::prev(p->cacheWrites.end()));
		    lod_generation_count_add_unlocked(
			p->generationCacheWriteCounts,
			p->cacheWrites.back().result.generation);
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

static SbBool
lod_compaction_working_set_available(
    const BObolLodServicePrivate *p,
    size_t estimate)
{
    if (!p || p->maxActiveWorkingSetBytes == SIZE_MAX || !estimate)
	return TRUE;
    if (!p->activeWorkingSetBytes)
	return TRUE;
    if (p->activeWorkingSetBytes >= p->maxActiveWorkingSetBytes)
	return FALSE;
    return estimate <=
	p->maxActiveWorkingSetBytes - p->activeWorkingSetBytes ? TRUE : FALSE;
}

static SbBool
lod_take_resident_compaction_unlocked(
    BObolLodServicePrivate *p,
    BObolResidentMeshCompactionWork &work)
{
    if (!p || p->residentMeshCompactionWork.empty())
	return FALSE;

    BObolResidentMeshCompactionWork &candidate =
	p->residentMeshCompactionWork.front();
    candidate.estimatedWorkingSetBytes =
	candidate.resident ?
	    candidate.resident->publishedBytes.load(
		std::memory_order_relaxed) : 0;
    if (!lod_compaction_working_set_available(
	    p, candidate.estimatedWorkingSetBytes))
	return FALSE;

    candidate.consumers.clear();
    for (const auto &consumer : p->residentMeshConsumerDemands) {
	if (consumer.second.assets.find(candidate.assetKey) !=
	    consumer.second.assets.end())
	    candidate.consumers.push_back(consumer.first);
    }
    const size_t occupied =
	p->residentMeshCompactionResultCount +
	p->residentMeshCompactionResultReservations;
    const size_t reservations = candidate.consumers.size();
    if (reservations > p->maxQueuedResults - std::min(
	    p->maxQueuedResults, occupied) &&
	occupied != 0)
	return FALSE;

    const auto target = p->residentMeshCompactionTargets.find(
	candidate.assetKey);
    if (target == p->residentMeshCompactionTargets.end()) {
	p->residentMeshCompactionQueuedAssets.erase(candidate.assetKey);
	p->residentMeshCompactionWork.pop_front();
	return FALSE;
    }
    candidate.target = target->second;

    work = std::move(candidate);
    p->residentMeshCompactionWork.pop_front();
    p->residentMeshCompactionResultReservations =
	reservations > SIZE_MAX -
	    p->residentMeshCompactionResultReservations ?
	SIZE_MAX :
	p->residentMeshCompactionResultReservations + reservations;
    p->residentMeshCompactionsInFlight++;
    p->executingTasks++;
    p->activeWorkingSetBytes =
	work.estimatedWorkingSetBytes > SIZE_MAX -
	    p->activeWorkingSetBytes ?
	SIZE_MAX :
	p->activeWorkingSetBytes + work.estimatedWorkingSetBytes;
    p->peakWorkingSetBytes = std::max(
	p->peakWorkingSetBytes, p->activeWorkingSetBytes);
    p->peakExecutingTasks = std::max(
	p->peakExecutingTasks, p->executingTasks);
    return TRUE;
}

static BObolLodResidentCompaction
lod_execute_resident_compaction(
    BObolLodServicePrivate *p,
    const BObolResidentMeshCompactionWork &work)
{
    BObolLodResidentCompaction result;
    if (!p || !work.resident || !work.target.revision)
	return result;

    const std::shared_ptr<BObolResidentMeshAsset> &resident = work.resident;
    const auto targetIsCurrent = [&]() {
	const auto found =
	    p->residentMeshCompactionTargets.find(work.assetKey);
	return found != p->residentMeshCompactionTargets.end() &&
	    found->second.revision == work.target.revision &&
	    found->second.useRevision == work.target.useRevision &&
	    found->second.cut == work.target.cut &&
	    found->second.channelMask == work.target.channelMask &&
	    found->second.evict == work.target.evict;
    };

    if (work.target.evict) {
	/* Commit removal under the documented service->asset lock order.  A
	 * provider increments useRevision while holding the service lock, so it
	 * is impossible for an old eviction to win after renewed use. */
	std::unique_lock<std::mutex> serviceLock(p->mutex);
	if (!targetIsCurrent() ||
	    resident->useRevision.load(std::memory_order_relaxed) !=
		work.target.useRevision)
	    return result;
	for (const auto &consumer : p->residentMeshConsumerDemands)
	    if (consumer.second.assets.find(work.assetKey) !=
		    consumer.second.assets.end())
		return result;
	const auto found = p->residentMeshes.find(work.assetKey);
	if (found == p->residentMeshes.end() || found->second != resident)
	    return result;
	std::unique_lock<std::mutex> residentLock(
	    resident->mutex, std::try_to_lock);
	if (!residentLock.owns_lock() || !resident->lod || !resident->mesh)
	    return result;
	const size_t priorBytes =
	    resident->publishedBytes.load(std::memory_order_relaxed);
	const size_t priorBackingBytes =
	    resident->publishedBackingPrefixBytes.load(
		std::memory_order_relaxed);
	p->residentMeshes.erase(found);
	if (resident->orderIndex < p->residentMeshOrder.size()) {
	    auto &ordered = p->residentMeshOrder[resident->orderIndex];
	    if (ordered.first == work.assetKey && ordered.second == resident)
		ordered.second.reset();
	}
	p->residentMeshEvictions++;
	serviceLock.unlock();

	result.priorBytes = priorBytes;
	result.residentBytes = 0;
	resident->publishedBytes.store(0, std::memory_order_relaxed);
	resident->publishedBackingPrefixBytes.store(
	    0, std::memory_order_relaxed);
	resident->publishedMinimumCut.store(-1, std::memory_order_relaxed);
	resident->publishedResidentCut.store(-1, std::memory_order_relaxed);
	resident->mesh.reset();
	if (resident->lod)
	    bobol_mesh_lod_destroy(resident->lod);
	resident->lod = NULL;
	lod_resident_mesh_bytes_replace(p->residentMeshBytes, priorBytes, 0);
	lod_resident_mesh_bytes_replace(
	    p->residentMeshBackingBytes, priorBackingBytes, 0);
	lod_resident_mesh_revision_advance(p->residentMeshRevision);
	if (priorBytes > priorBackingBytes)
	    lod_resident_mesh_revision_advance(
		p->residentMeshAdmissionRevision);
	return result;
    }

    BObolLodProgressiveMeshPtr mesh;
    BObolLodProgressiveMeshTrimPtr preparedTrim;
    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    int preparedTargetCut = -1;
    {
	/* Retain an immutable source generation while constructing the shorter
	 * candidate.  No shared mesh state changes in this phase. */
	std::lock_guard<std::mutex> residentLock(resident->mutex);
	if (!resident->lod || !resident->mesh)
	    return result;
	mesh = resident->mesh;
	const int currentCut = mesh->residentCut();
	const int minimumCut = mesh->minimumCut();
	const int maximumCut = mesh->maximumCut();
	if (currentCut < 0 || minimumCut < 0 ||
	    maximumCut < minimumCut)
	    return result;
	preparedTargetCut = std::min(currentCut,
	    std::max(minimumCut, std::min(maximumCut,
		work.target.cut < 0 ? minimumCut : work.target.cut)));
	if (preparedTargetCut < currentCut) {
	    if (!bobol_mesh_lod_hierarchy_info_get(
		    resident->lod, &hierarchy))
		return result;
	}
    }
    if (preparedTargetCut < mesh->residentCut()) {
	preparedTrim = mesh->prepareTrim(preparedTargetCut);
	if (!preparedTrim)
	    return result;
    }

    /* Publish only after revalidating the exact target and provider-use
     * epoch.  Holding the short service lock through commit makes the check
     * and immutable-generation swap atomic with respect to new requests;
     * the expensive prefix copy above occurs without either service lock. */
    std::unique_lock<std::mutex> serviceLock(p->mutex);
    if (!targetIsCurrent() ||
	resident->useRevision.load(std::memory_order_relaxed) !=
	    work.target.useRevision)
	return result;
    const auto indexed = p->residentMeshes.find(work.assetKey);
    if (indexed == p->residentMeshes.end() || indexed->second != resident)
	return result;
    std::unique_lock<std::mutex> residentLock(
	resident->mutex, std::try_to_lock);
    if (!residentLock.owns_lock() || !resident->lod ||
	resident->mesh != mesh)
	return result;
    const int currentCut = mesh->residentCut();
    const int minimumCut = mesh->minimumCut();
    const int maximumCut = mesh->maximumCut();
    if (currentCut < 0 || minimumCut < 0 ||
	maximumCut < minimumCut)
	return result;
    const int targetCut = std::min(currentCut,
	std::max(minimumCut, std::min(maximumCut,
	    work.target.cut < 0 ? minimumCut : work.target.cut)));
    if (targetCut != preparedTargetCut)
	return result;
    const size_t priorBytes =
	resident->publishedBytes.load(std::memory_order_relaxed);
    const size_t priorBackingBytes =
	resident->publishedBackingPrefixBytes.load(
	    std::memory_order_relaxed);

    if (targetCut >= currentCut) {
	/* Stable demand already matches the immutable prefix.  Release only the
	 * reloadable cache-reader duplicate after the same revision guard. */
	if (!bobol_mesh_lod_resident_prefix_bytes(resident->lod))
	    return result;
	serviceLock.unlock();
	result.priorBytes = priorBytes;
	bobol_mesh_lod_memshrink(resident->lod);
	result.residentBytes = lod_resident_asset_bytes(*resident);
	resident->publishedBytes.store(
	    result.residentBytes, std::memory_order_relaxed);
	resident->publishedBackingPrefixBytes.store(
	    0, std::memory_order_relaxed);
	if (result.residentBytes != priorBytes) {
	    lod_resident_mesh_bytes_replace(
		p->residentMeshBytes, priorBytes, result.residentBytes);
	    lod_resident_mesh_revision_advance(p->residentMeshRevision);
	}
	if (priorBackingBytes)
	    lod_resident_mesh_bytes_replace(
		p->residentMeshBackingBytes, priorBackingBytes, 0);
	return result;
    }

    if (!preparedTrim || !mesh->commitTrim(preparedTrim) ||
	mesh->residentCut() != targetCut)
	return BObolLodResidentCompaction();
    resident->publishedMinimumCut.store(
	hierarchy.min_cut, std::memory_order_relaxed);
    resident->publishedResidentCut.store(
	targetCut, std::memory_order_relaxed);
    serviceLock.unlock();

    result.assetKey = work.assetKey.c_str();
    result.progressiveMesh = mesh;
    result.residentCut = targetCut;
    result.channelMask = work.target.channelMask & 3u;
    /* The replacement immutable generation is self-contained.  Release the
     * duplicate cache-reader prefix before publishing its byte accounting. */
    bobol_mesh_lod_memshrink(resident->lod);
    result.residentBytes = lod_resident_asset_bytes(*resident);
    resident->publishedBytes.store(
	result.residentBytes, std::memory_order_relaxed);
    resident->publishedBackingPrefixBytes.store(
	0, std::memory_order_relaxed);
    lod_resident_mesh_bytes_replace(
	p->residentMeshBytes, priorBytes, result.residentBytes);
    if (priorBackingBytes)
	lod_resident_mesh_bytes_replace(
	    p->residentMeshBackingBytes, priorBackingBytes, 0);
    lod_resident_mesh_revision_advance(p->residentMeshRevision);
    const size_t priorStableBytes = priorBackingBytes >= priorBytes ?
	0 : priorBytes - priorBackingBytes;
    if (result.residentBytes < priorStableBytes)
	lod_resident_mesh_revision_advance(p->residentMeshAdmissionRevision);
    residentLock.unlock();

    if (result.channelMask) {
	const int drawMode = result.channelMask == 3u ?
	    BOBOL_LOD_DRAW_HIDDEN_LINE :
	    (result.channelMask & 1u ?
		BOBOL_LOD_DRAW_WIRE : BOBOL_LOD_DRAW_SHADED);
	result.preparedCadGeometry = mesh->prepareCadGeometry(
	    drawMode, &result.preparedCadGeometryRevision);
    }
    return result;
}

static void
lod_finish_resident_compaction(
    BObolLodServicePrivate *p,
    const BObolResidentMeshCompactionWork &work,
    BObolLodResidentCompaction &&result)
{
    SbBool notifyResultReady = FALSE;
    const bool completed = result.assetKey.getLength() > 0 &&
	result.progressiveMesh && result.residentCut >= 0;
    const bool reclaimedStorage =
	result.priorBytes > result.residentBytes;
    {
	std::lock_guard<std::mutex> lock(p->mutex);
	if (p->residentMeshCompactionsInFlight > 0)
	    p->residentMeshCompactionsInFlight--;
	if (p->executingTasks > 0)
	    p->executingTasks--;
	p->activeWorkingSetBytes =
	    work.estimatedWorkingSetBytes >= p->activeWorkingSetBytes ?
	    0 : p->activeWorkingSetBytes -
		work.estimatedWorkingSetBytes;
	p->residentMeshCompactionResultReservations =
	    work.consumers.size() >=
		p->residentMeshCompactionResultReservations ?
	    0 : p->residentMeshCompactionResultReservations -
		work.consumers.size();
	/* A newer complete demand snapshot may replace this target while the
	 * worker is constructing its immutable candidate.  Do not let the old
	 * completion erase that newer obligation; keep the per-asset queue token
	 * and immediately requeue one work item for the latest target. */
	const auto currentTarget =
	    p->residentMeshCompactionTargets.find(work.assetKey);
	const bool superseded =
	    currentTarget != p->residentMeshCompactionTargets.end() &&
	    currentTarget->second.revision != work.target.revision;
	const auto currentResident = p->residentMeshes.find(work.assetKey);
	if (superseded && !p->stopping &&
	    currentResident != p->residentMeshes.end() &&
	    currentResident->second == work.resident) {
	    BObolResidentMeshCompactionWork successor;
	    successor.assetKey = work.assetKey;
	    successor.resident = work.resident;
	    p->residentMeshCompactionWork.push_back(std::move(successor));
	} else {
	    p->residentMeshCompactionQueuedAssets.erase(work.assetKey);
	    if (currentTarget == p->residentMeshCompactionTargets.end() ||
		currentTarget->second.revision == work.target.revision)
		p->residentMeshCompactionTargets.erase(work.assetKey);
	}

	if ((completed || reclaimedStorage) && !p->stopping)
	    p->residentMeshCompactions++;
	if (completed && !p->stopping) {
	    const size_t resultCountBefore =
		p->residentMeshCompactionResultCount;
	    for (uint64_t consumerId : work.consumers) {
		const auto consumer =
		    p->residentMeshConsumerDemands.find(consumerId);
		if (consumer == p->residentMeshConsumerDemands.end())
		    continue;
		const auto demand =
		    consumer->second.assets.find(work.assetKey);
		if (demand == consumer->second.assets.end() ||
		    demand->second.cut > result.residentCut)
		    continue;
		p->residentMeshCompactionResults[consumerId].push_back(
		    result);
		p->residentMeshCompactionResultCount++;
	    }
	    notifyResultReady =
		resultCountBefore == 0 &&
		p->residentMeshCompactionResultCount > 0 ? TRUE : FALSE;
	}
    }
    p->workerCv.notify_all();
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
	BObolResidentMeshCompactionWork compaction;
	SbBool runCompaction = FALSE;

	{
	    std::unique_lock<std::mutex> lock(p->mutex);

	    std::deque<BObolLodWorkItem>::iterator ready;
	    for (;;) {
		if (p->stopping)
		    return;
		ready = lod_find_ready_task(p);
		if (ready != p->pending.end())
		    break;
		if (lod_take_resident_compaction_unlocked(
			p, compaction)) {
		    runCompaction = TRUE;
		    break;
		}
		p->workerCv.wait(lock);
	    }

	    if (runCompaction) {
		/* Counters and reservations were installed by the take helper. */
	    } else {
	    const int qualityTier = ready->task.request.qualityTier;
		    item = std::move(*ready);
		    p->pending.erase(ready);
		    lod_generation_count_remove_unlocked(
			p->generationPendingTaskCounts,
			item.task.generation);
		    lod_generation_count_add_unlocked(
			p->generationExecutingTaskCounts,
			item.task.generation);
		    lod_pending_quality_remove(p, qualityTier);
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
	}

	if (runCompaction) {
	    bobol_lod_working_set_acquire(
		compaction.estimatedWorkingSetBytes);
	    BObolLodResidentCompaction result =
		lod_execute_resident_compaction(p, compaction);
	    lod_finish_resident_compaction(
		p, compaction, std::move(result));
	    bobol_lod_working_set_release(
		compaction.estimatedWorkingSetBytes);
	    continue;
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

	    p->cacheWriteSlots.erase(
		lod_result_slot_map_key(p->cacheWrites.front().result));
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
	    lod_generation_count_remove_unlocked(
		p->generationCacheWriteCounts, item.result.generation);
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

SbBool
BObolLodService::ensureWorkerCount(size_t workerCount)
{
    if (workerCount == 0)
	workerCount = 1;

    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (!this->p->running || this->p->stopping)
	return FALSE;
    try {
	while (this->p->workers.size() < workerCount)
	    this->p->workers.push_back(
		std::thread(lod_worker_loop, this->p));
    } catch (...) {
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
	this->p->pendingQualityCounts.clear();
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
	this->p->queuedResultRequestKeyCounts.clear();
	this->p->cacheWriteSlots.clear();
	this->p->resultSlots.clear();
	this->p->completed.clear();
	this->p->taskGenerations.clear();
	this->p->cancelledGenerations.clear();
	this->p->cancelledGenerationOrder.clear();
	this->p->generationTaskCounts.clear();
	this->p->generationPendingTaskCounts.clear();
	this->p->generationExecutingTaskCounts.clear();
	this->p->generationDelayedTaskCounts.clear();
	this->p->generationResultCounts.clear();
	this->p->generationCacheWriteCounts.clear();
	this->p->residentMeshConsumerDemands.clear();
	this->p->residentMeshOrder.clear();
	this->p->residentMeshCompactionWork.clear();
	this->p->residentMeshCompactionQueuedAssets.clear();
	this->p->residentMeshCompactionTargets.clear();
	this->p->residentMeshCompactionResults.clear();
	this->p->residentMeshCompactionsInFlight = 0;
	this->p->residentMeshCompactionResultCount = 0;
	this->p->residentMeshCompactionResultReservations = 0;
	residentMeshes.swap(this->p->residentMeshes);
	this->p->residentMeshBytes.store(0, std::memory_order_relaxed);
	this->p->residentMeshBackingBytes.store(
	    0, std::memory_order_relaxed);
	this->p->residentMeshGrowthReservationBytes = 0;
	lod_resident_mesh_revision_advance(
	    this->p->residentMeshAdmissionRevision);
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

size_t
BObolLodService::workerCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->workers.size();
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

/* Select one bounded presentation step toward the producer-resolved
 * screen-error target.  This helper only prevents a first box->mesh
 * replacement from allocating and uploading an arbitrarily large cumulative
 * prefix; it does not mutate asynchronous request identity. */
static int
lod_progressive_delivery_cut(
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int requestedCut, int currentCut,
    const BObolMeshLodProvider &provider)
{
    int target = std::max(hierarchy.min_cut,
	std::min(hierarchy.max_cut, requestedCut));
    if (provider.useDeliveryCutLimit &&
	provider.deliveryCutLimit >= hierarchy.min_cut)
	target = std::min(target, provider.deliveryCutLimit);
    if (!provider.progressiveDelivery || provider.useForcedCut ||
	provider.reset || currentCut >= target)
	return target;

    const size_t currentFaces =
	currentCut >= hierarchy.min_cut ?
	hierarchy.cuts[currentCut].face_count : 0;
    const size_t currentPoints =
	currentCut >= hierarchy.min_cut ?
	hierarchy.cuts[currentCut].point_count : 0;
    const size_t faceBudget = lod_refinement_growth_budget(currentFaces,
	provider.initialRefinementFaceBudget,
	provider.refinementGrowthFactor);
    const size_t pointBudget = lod_refinement_growth_budget(currentPoints,
	provider.initialRefinementPointBudget,
	provider.refinementGrowthFactor);

    int selected = hierarchy.min_cut;
    for (int cut = hierarchy.min_cut; cut <= target; ++cut) {
	if (hierarchy.cuts[cut].face_count > faceBudget ||
	    hierarchy.cuts[cut].point_count > pointBudget)
	    break;
	selected = cut;
    }

    /*
     * Return the richest population that fits the provider's face/point
     * growth allowance.  Nominal cut adjacency is deliberately irrelevant:
     * one cut may add no arrays while the next adds millions.  The budget is
     * the bounded work contract, not "+1 cut".
     */
    if (currentCut >= hierarchy.min_cut && selected <= currentCut)
	selected = std::min(target, currentCut + 1);
    return std::min(selected, target);
}

BObolLodResult
BObolLodService::realizeResidentMeshLod(
    const BObolLodRequest &request,
    const BObolMeshLodProvider &provider)
{
    struct db_i *dbip = provider.getDatabase();
    if (!dbip)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_ERROR,
	    "resident mesh provider has no database");
    const std::string databaseIdentity =
	request.databaseId.getLength() > 0 ?
	    request.databaseId.getString() :
	    (dbip->dbi_filename ? dbip->dbi_filename : "");
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
	    resident->databaseIdentity = databaseIdentity;
	    resident->name = name;
	    resident->mesh = std::make_shared<BObolLodProgressiveMesh>();
	    this->p->residentMeshes[assetKey.value.getString()] = resident;
	    resident->orderIndex = this->p->residentMeshOrder.size();
	    this->p->residentMeshOrder.push_back(std::make_pair(
		std::string(assetKey.value.getString()), resident));
	} else {
	    resident = found->second;
	}
	resident->useRevision.fetch_add(1, std::memory_order_relaxed);
    }

    std::lock_guard<std::mutex> residentLock(resident->mutex);
    if (resident->databaseIdentity != databaseIdentity ||
	resident->name != name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_STALE,
	    "resident mesh asset identity collision");

    bool persistentHierarchyAvailable = resident->lod != NULL;
    if (!resident->lod) {
	/*
	 * Compact requests captured the validated immutable content key while
	 * their source records were streamed.  Open that payload directly and
	 * keep its hierarchy-header snapshot through the first prefix read.
	 * Falling back to the named API preserves the non-compact/legacy path.
	 * At distinct-asset scale this avoids one database lookup, one name-cache
	 * lookup, and one LMDB read transaction per successful warm task.
	 */
	const bool exactVariant = provider.meshAssetContentHash != 0;
	if (exactVariant)
	    resident->lod = bobol_mesh_lod_get_cached_prefix(
		dbip, provider.meshAssetContentHash);
	if (!resident->lod && !exactVariant)
	    resident->lod = bobol_mesh_lod_get_named_cached_prefix(
		dbip, name);
	persistentHierarchyAvailable = resident->lod != NULL;
	if (!resident->lod && exactVariant && provider.generateBrepVariant) {
	    const int64_t variantStarted = bu_gettime();
	    struct bg_tess_tol ttol = BG_TESS_TOL_INIT_TOL;
	    ttol.abs = std::max(0.0, provider.brepTessellationAbsTol);
	    ttol.rel = std::max(0.0, provider.brepTessellationRelTol);
	    ttol.norm = std::max(0.0, provider.brepTessellationNormTol);
	    BObolSourceMeshRequest generatedRequest;
	    std::shared_ptr<BObolStagedSourceMesh> generated =
		bobol_database_brep_staged_mesh_variant(
		    dbip, name, &ttol,
		    provider.meshAssetContentHash,
		    request.sourceRevision.value(), generatedRequest);
	    if (generated && generated->isValid() &&
		bobol_mesh_lod_cache_store_mesh_variant(
		    dbip, name, generated->points,
		    generated->pointCount, generated->normals,
		    generated->faces, generated->faceCount,
		    provider.meshAssetContentHash,
		    generated->shadedCullBackfaces,
		    &resident->status) == BRLCAD_OK) {
		resident->lod = bobol_mesh_lod_get_cached_prefix(
		    dbip, provider.meshAssetContentHash);
	    }
	    if (getenv("BOBOL_DRAW_TIMING"))
		bu_log("BObol BREP variant object=%s content=%llu "
		       "rel_tol=%.17g faces=%zu elapsed_ms=%.3f status=%s\n",
		       name,
		       static_cast<unsigned long long>(
			   provider.meshAssetContentHash),
		       ttol.rel,
		       generated ? generated->faceCount : 0,
		       static_cast<double>(bu_gettime() - variantStarted) /
			   1000.0,
		       resident->lod ? "ready" : "failed");
	}
	if (resident->lod) {
	    struct directory *assetDirectory = db_lookup(dbip, name,
		LOOKUP_QUIET);
	    resident->status.directory_found = assetDirectory ? 1 : 0;
	    resident->status.is_bot = assetDirectory &&
		assetDirectory->d_minor_type == DB5_MINORTYPE_BRLCAD_BOT ? 1 : 0;
	    resident->status.has_cache_key = 1;
	    resident->status.has_cached_payload = 1;
	    resident->status.stale_cache_entry = 0;
	    resident->status.cache_key =
		bobol_mesh_lod_cache_key_get(resident->lod);
	} else {
	    if (bobol_mesh_lod_cache_status(dbip, name,
		    &resident->status) != BRLCAD_OK)
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider could not query cache status");
	    if ((exactVariant || !resident->status.has_cache_key ||
		 !resident->status.has_cached_payload ||
		 resident->status.stale_cache_entry) &&
		provider.refreshMissing) {
		const int64_t refreshStarted = bu_gettime();
		const std::shared_ptr<const BObolStagedSourceMesh> &staged =
		    provider.stagedSource;
		const uint64_t stagedFaceCount = staged ?
		    (staged->faceCount ? staged->faceCount :
		     (staged->bot ? staged->bot->num_faces : 0)) : 0;
		const uint64_t stagedPointCount = staged ?
		    (staged->pointCount ? staged->pointCount :
		     (staged->bot ? staged->bot->num_vertices : 0)) : 0;
		const bool stagedMatches =
		    staged && staged->isValid() &&
		    staged->sourceRevision == request.sourceRevision &&
		    bu_strcmp(staged->assetName.getString(), name) == 0 &&
		    (!request.sourceCounts.faceCount ||
		     request.sourceCounts.faceCount == stagedFaceCount) &&
		    (!request.sourceCounts.pointCount ||
		     request.sourceCounts.pointCount == stagedPointCount);
		int refreshResult = BRLCAD_ERROR;
		if (stagedMatches && staged->bot) {
		    resident->lod =
			bobol_mesh_lod_cache_refresh_from_bot_open(
			    dbip, name, staged->bot,
			    &resident->status);
		    refreshResult = resident->lod ?
			BRLCAD_OK : BRLCAD_ERROR;
		} else if (stagedMatches) {
		    refreshResult = bobol_mesh_lod_cache_store_mesh_variant(
			dbip, name, staged->points,
			staged->pointCount, staged->normals, staged->faces,
			staged->faceCount, staged->contentKey,
			staged->shadedCullBackfaces, &resident->status);
		    if (refreshResult == BRLCAD_OK &&
			resident->status.has_cache_key)
			resident->lod = bobol_mesh_lod_get_cached_prefix(
			    dbip, resident->status.cache_key);
		} else if (!exactVariant) {
		    refreshResult = bobol_mesh_lod_cache_refresh(
			dbip, name, &resident->status);
		} else {
		    return lod_provider_status_result(request,
			BOBOL_LOD_PROVIDER_CACHE_MISS,
			"resident mesh provider has no exact staged variant");
		}
		if (refreshResult != BRLCAD_OK)
		    return lod_provider_status_result(request,
			BOBOL_LOD_PROVIDER_CACHE_MISS,
			"resident mesh provider could not refresh cache entry");
		if (getenv("BOBOL_DRAW_TIMING"))
		    bu_log("[obol-timing] pop cache: refresh %-24s %8.1f ms "
			   "staged=%d\n", name,
			   (bu_gettime() - refreshStarted) / 1000.0,
			   stagedMatches ? 1 : 0);
	    }
	    if (!resident->lod && exactVariant)
		resident->lod = bobol_mesh_lod_get_cached_prefix(
		    dbip, provider.meshAssetContentHash);
	    if (!resident->lod && !exactVariant &&
		resident->status.has_cache_key)
		resident->lod = bobol_mesh_lod_get_cached_prefix(
		    dbip, resident->status.cache_key);
	    if (!resident->lod && !exactVariant)
		resident->lod = bobol_mesh_lod_get_named_cached_prefix(
		    dbip, name);
	}
	if (!resident->lod)
	    return lod_provider_status_result(request,
		resident->status.stale_cache_entry ?
		    BOBOL_LOD_PROVIDER_STALE :
		    BOBOL_LOD_PROVIDER_CACHE_MISS,
		"resident mesh provider has no cache payload");
    }

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (!bobol_mesh_lod_hierarchy_info_get(resident->lod, &hierarchy))
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_CACHE_MISS,
	    "resident mesh provider loaded no hierarchy metadata");
    int requestedCut = provider.useForcedCut ?
	provider.forcedCut : request.requestedCut;
    if (!provider.useForcedCut && request.projectedPixelDiameter > 0.0f &&
	request.targetPixelError > 0.0f)
	requestedCut = bobol_mesh_lod_select_cut(&hierarchy,
	    request.projectedPixelDiameter, request.targetPixelError);
    if (requestedCut < hierarchy.min_cut)
	requestedCut = hierarchy.min_cut;
    if (requestedCut > hierarchy.max_cut)
	requestedCut = hierarchy.max_cut;
    resident->publishedMinimumCut.store(
	hierarchy.min_cut, std::memory_order_relaxed);

    const size_t priorResidentBytes =
	resident->publishedBytes.load(std::memory_order_relaxed);
    const size_t priorBackingBytes =
	resident->publishedBackingPrefixBytes.load(
	    std::memory_order_relaxed);
    const size_t priorStableBytes =
	priorBackingBytes >= priorResidentBytes ?
	    0 : priorResidentBytes - priorBackingBytes;
    int currentCut = bobol_mesh_lod_current_cut(resident->lod);
    const int publishedCut =
	resident->mesh && resident->mesh->isValid() ?
	    resident->mesh->residentCut() : -1;
    const int deliveryCut = provider.useCurrentDrawCut ?
	provider.currentDrawCut :
	(publishedCut >= 0 ? publishedCut : currentCut);
    int residentTarget = requestedCut;
    if (requestedCut >= 0)
	residentTarget = lod_progressive_delivery_cut(hierarchy, requestedCut,
	    deliveryCut, provider);
    if (residentTarget < hierarchy.min_cut)
	residentTarget = hierarchy.min_cut;
    if (residentTarget > hierarchy.max_cut)
	residentTarget = hierarchy.max_cut;
    /*
     * Establish the minimum useful mesh as its own publication.  Besides
     * minimizing time-to-first-content, this prevents early workers from
     * spending the stable byte allowance on optional suffixes before later
     * visible assets have revealed their mandatory coverage floor.  The next
     * bounded view pass may jump directly to the richest face- and
     * byte-admitted target; this is two-phase coverage, not nominal
     * cut-by-cut walking.
     */
    const bool warmFirstPrefetch =
	publishedCut < hierarchy.min_cut &&
	persistentHierarchyAvailable &&
	provider.prefetchCachedTargetOnFirstPublication &&
	!provider.reset && !provider.useForcedCut &&
	requestedCut >= hierarchy.min_cut;
    if (warmFirstPrefetch)
	residentTarget = std::min(hierarchy.max_cut, requestedCut);
    else if (publishedCut < hierarchy.min_cut)
	residentTarget = hierarchy.min_cut;
    /* A warm first task may populate the complete view-demanded immutable
     * prefix, but its publication remains the cheapest useful mesh.  The
     * next presentation pass can select any budget-admitted cut from those
     * already resident arrays without another cache task. */
    int drawTarget = warmFirstPrefetch ? hierarchy.min_cut : residentTarget;
    if (provider.usePresentationCutLimit &&
	provider.presentationCutLimit >= hierarchy.min_cut)
	drawTarget = std::min(drawTarget,
	    provider.presentationCutLimit);
    int loadTarget = residentTarget;
    if (provider.compactResident && publishedCut >= 0 &&
	residentTarget < publishedCut) {
	/* A cut is not a bounded unit of memory: one Lucy hierarchy step can
	 * add millions of faces.  Stable reclamation therefore retains exactly
	 * the pixel-demanded prefix. */
	loadTarget = residentTarget;
    } else if (publishedCut >= 0 && residentTarget <= publishedCut &&
	!provider.reset) {
	loadTarget = publishedCut;
    }

    BObolResidentMeshGrowthReservation growthReservation(this->p);
    SbBool memoryLimited = FALSE;
    loadTarget = growthReservation.admit(*resident, hierarchy,
	loadTarget, publishedCut, priorStableBytes, memoryLimited);
    if (drawTarget > loadTarget)
	drawTarget = loadTarget;

    /*
     * The immutable renderer generation is authoritative residency.  Stable
     * maintenance intentionally drops the cache reader's duplicate prefix;
     * reloading that prefix merely to return an already drawable cut defeats
     * compaction and can produce a load/shrink retry cycle.
     */
    const bool retainedTargetDrawable =
	resident->mesh && resident->mesh->isValid() &&
	resident->mesh->canDrawCut(loadTarget);
    const bool loadNeeded =
	provider.reset || !retainedTargetDrawable ||
	(publishedCut >= 0 && loadTarget != publishedCut);
    int64_t prefixLoadMicroseconds = 0;
    int64_t generationBuildMicroseconds = 0;
    int64_t preparedGeometryMicroseconds = 0;
    int64_t legacyPayloadMicroseconds = 0;
    int residentCut = publishedCut;
    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    const auto populateInfoFromRetainedMesh = [&]() {
	bobol_mesh_lod_info_init(&info);
	info.active_cut = residentCut;
	info.face_count = resident->mesh->faceCount(residentCut);
	info.point_count = resident->mesh->pointCount(residentCut);
	info.point_orig_count = info.point_count;
	info.normal_count = hierarchy.has_normals ?
	    info.face_count * 3 : 0;
	info.has_faces = info.face_count ? 1 : 0;
	info.has_points = info.point_count ? 1 : 0;
	info.has_original_points = info.has_points;
	info.has_snapped_points = 0;
	info.has_normals = hierarchy.has_normals;
	info.shaded_cull_backfaces =
	    hierarchy.shaded_cull_backfaces;
	const SbBox3f bounds = resident->mesh->bounds();
	const SbVec3f minimum = bounds.getMin();
	const SbVec3f maximum = bounds.getMax();
	VSET(info.bmin, minimum[0], minimum[1], minimum[2]);
	VSET(info.bmax, maximum[0], maximum[1], maximum[2]);
    };
    if (loadNeeded) {
	/* Persistent PoP records are already split by activation cut.  Once a
	 * quiet compaction has released the cache reader's duplicate prefix, grow
	 * the immutable renderer generation from only the missing cache suffix.
	 * Corner-normal vertex splitting still needs whole-prefix context and uses
	 * the conservative cumulative fallback. */
	SbBool suffixExtended = FALSE;
	if (!provider.reset && publishedCut >= hierarchy.min_cut &&
	    loadTarget > publishedCut && !hierarchy.has_normals &&
	    resident->mesh && resident->mesh->isValid()) {
	    const int64_t generationBuildStarted = bu_gettime();
	    suffixExtended = resident->mesh->extendFromCache(
		resident->lod, hierarchy, loadTarget,
		hierarchy.shaded_cull_backfaces ? TRUE : FALSE);
	    generationBuildMicroseconds = std::max<int64_t>(
		0, bu_gettime() - generationBuildStarted);
	    if (suffixExtended) {
		residentCut = loadTarget;
		populateInfoFromRetainedMesh();
	    }
	}
	if (!suffixExtended) {
	    const int64_t prefixLoadStarted = bu_gettime();
	    residentCut = bobol_mesh_lod_load_resident_cut(
		resident->lod, loadTarget, provider.reset);
	    prefixLoadMicroseconds = std::max<int64_t>(
		0, bu_gettime() - prefixLoadStarted);
	    if (residentCut < 0)
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_CACHE_MISS,
		    "resident mesh provider could not load the requested prefix");
	    /* hierarchy_info includes the cache handle's current resident cut.
	     * Loading a cheaper prefix changes that value; retaining the snapshot
	     * taken before the load makes an otherwise valid immutable generation
	     * fail its cut-consistency check.  This is common when a shared warm
	     * cache handle last served a richer view and a new consumer first asks
	     * for coverage minimum. */
	    if (!bobol_mesh_lod_hierarchy_info_get(
		    resident->lod, &hierarchy) ||
		hierarchy.resident_cut != residentCut)
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider could not refresh loaded hierarchy metadata");
	    struct BObolMeshLodData data;
	    if (!bobol_mesh_lod_info_get(resident->lod, &info) ||
		!bobol_mesh_lod_data_get(resident->lod, &data))
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_CACHE_MISS,
		    "resident mesh provider loaded no mesh data");
	    const int64_t generationBuildStarted = bu_gettime();
	    if (!resident->mesh->update(data, hierarchy, residentCut,
		    hierarchy.shaded_cull_backfaces ? TRUE : FALSE))
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider could not publish the retained asset");
	    generationBuildMicroseconds = std::max<int64_t>(
		0, bu_gettime() - generationBuildStarted);
	}
	if (publishedCut < residentCut) {
	    std::lock_guard<std::mutex> lock(this->p->mutex);
	    this->p->residentMeshCacheLoads++;
	}
    } else {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->residentMeshHits++;
	populateInfoFromRetainedMesh();
    }
    resident->publishedResidentCut.store(
	residentCut, std::memory_order_relaxed);

    int drawCut = drawTarget;
    if (drawCut < hierarchy.min_cut)
	drawCut = hierarchy.min_cut;
    if (drawCut > residentCut)
	drawCut = residentCut;

    BObolLodResult result =
	bobol_lod_result_from_mesh_lod_info(request, info, &resident->status);
    result.resolvedCut = requestedCut;
    result.geometry.cacheKey = assetKey;
    result.geometry.activeCut = drawCut;
    result.progressiveMesh = resident->mesh;
    result.residentCut = residentCut;
    result.residentAdmissionRevision =
	growthReservation.revision();
    result.counts.faceCount = resident->mesh->faceCount(drawCut);
    result.counts.pointCount = resident->mesh->pointCount(drawCut);
    result.counts.originalPointCount = result.counts.pointCount;
    result.counts.normalCount = info.has_normals ?
	result.counts.faceCount * 3 : 0;
    result.bounds = resident->mesh->bounds();
    result.hasSnappedPoints = FALSE;
    result.hasNormals = hierarchy.has_normals ? TRUE : FALSE;
    result.shadedCullBackfaces = resident->mesh->cullBackfaces();
    result.terminal = drawCut >= std::max(hierarchy.min_cut,
	std::min(hierarchy.max_cut, requestedCut)) ? TRUE : FALSE;
    /* A representation-band admission cap is not a reason to stop walking
     * the admitted band's PoP prefix.  Publish it only on that band's
     * terminal result; treating it as an immediate resident-memory failure
     * would strand a BREP at its first coarse cut. */
    result.memoryLimited = memoryLimited ||
	(provider.brepVariantMemoryLimited && result.terminal);

    /*
     * Prepare the renderer target on this worker for every CAD consumer, not
     * only compact occurrences.  Source-wide database nodes can adopt the
     * same immutable allocation and avoid rebuilding/copying mesh vectors on
     * the GUI thread.  The legacy shape payload remains populated for the
     * explicit SoBRLMeshShape route until that renderer is migrated.
     */
    const int64_t preparedGeometryStarted = bu_gettime();
    result.preparedCadGeometry =
	resident->mesh->prepareCadGeometry(
	    request.drawMode, &result.preparedCadGeometryRevision);
    preparedGeometryMicroseconds = std::max<int64_t>(
	0, bu_gettime() - preparedGeometryStarted);
    /* The renderer generation above is self-contained.  Keeping the cache
     * reader's cumulative arrays duplicates every point/index and previously
     * required one background compaction job per asset merely to release
     * them.  Suffix growth reads its missing ranges directly from persistent
     * storage, so eager release preserves progressive extension while making
     * first publication the final backing-storage cleanup step. */
    if (provider.shrinkAfterCopy && resident->lod)
	bobol_mesh_lod_memshrink(resident->lod);

    const size_t residentBytes = lod_resident_asset_bytes(*resident);
    const size_t backingBytes =
	bobol_mesh_lod_resident_prefix_bytes(resident->lod);
    resident->publishedBytes.store(
	residentBytes, std::memory_order_relaxed);
    resident->publishedBackingPrefixBytes.store(
	backingBytes,
	std::memory_order_relaxed);
    if (residentBytes != priorResidentBytes) {
	lod_resident_mesh_bytes_replace(this->p->residentMeshBytes,
	    priorResidentBytes, residentBytes);
	lod_resident_mesh_revision_advance(
	    this->p->residentMeshRevision);
    }
    if (backingBytes != priorBackingBytes)
	lod_resident_mesh_bytes_replace(
	    this->p->residentMeshBackingBytes,
	    priorBackingBytes, backingBytes);
    const size_t stableBytes = backingBytes >= residentBytes ?
	0 : residentBytes - backingBytes;
    if (stableBytes < priorStableBytes)
	lod_resident_mesh_revision_advance(
	    this->p->residentMeshAdmissionRevision);
    /* Publish exact totals before making this reservation available to a
     * peer.  Otherwise many workers can briefly admit against the gap
     * between estimated release and exact accounting. */
    growthReservation.release();

    if (request.occurrenceKey.getLength() == 0) {
	const int64_t legacyPayloadStarted = bu_gettime();
	if (!resident->mesh->copyCut(result.mesh, drawCut)) {
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_ERROR,
		"resident mesh provider could not materialize legacy payload");
	}
	legacyPayloadMicroseconds = std::max<int64_t>(
	    0, bu_gettime() - legacyPayloadStarted);
    }

    if (loadNeeded && getenv("BOBOL_DRAW_TIMING"))
	bu_log("[obol-timing] resident prefix %-24s cut=%d->%d "
	       "load=%8.1f ms generation=%8.1f ms prepare=%8.1f ms "
	       "legacy=%8.1f ms\n",
	       name, publishedCut, residentCut,
	       prefixLoadMicroseconds / 1000.0,
	       generationBuildMicroseconds / 1000.0,
	       preparedGeometryMicroseconds / 1000.0,
	       legacyPayloadMicroseconds / 1000.0);

    const char *traceFilter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (traceFilter && traceFilter[0] &&
	(strstr(name, traceFilter) ||
	 (request.objectPath.getLength() > 0 &&
	  strstr(request.objectPath.getString(), traceFilter)))) {
	bu_log("BObol resident LoD trace object=%s submitted_cut=%d "
	       "resolved_cut=%d draw_cut=%d resident_cut=%d "
	       "faces=%zu points=%zu "
	       "asset_revision=%llu load=%d compact=%d terminal=%d\n",
	       name, request.requestedCut, result.resolvedCut, drawCut,
	       residentCut,
	       result.counts.faceCount, result.counts.pointCount,
	       static_cast<unsigned long long>(resident->mesh->revision()),
	       loadNeeded ? 1 : 0, provider.compactResident ? 1 : 0,
	       result.terminal ? 1 : 0);
	if (getenv("BOBOL_LOD_TRACE_HIERARCHY")) {
	    for (int cut = hierarchy.min_cut;
		 cut <= hierarchy.max_cut; ++cut)
		bu_log("BObol resident hierarchy object=%s cut=%d "
		       "faces=%zu points=%zu%s\n",
		       name, cut, hierarchy.cuts[cut].face_count,
		       hierarchy.cuts[cut].point_count,
		       cut == hierarchy.max_cut ? " terminal" : "");
	}
    }
    return result;
}

size_t
BObolLodService::scheduleResidentMeshCompaction(
    uint64_t consumerId,
    uint64_t demandRevision,
    const std::vector<BObolLodResidentDemand> &demands,
    SbBool *planningComplete)
{
    if (planningComplete)
	*planningComplete = FALSE;
    if (!consumerId)
	return 0;

    const uint64_t residentRevisionAtEntry =
	this->p->residentMeshRevision.load(std::memory_order_relaxed);
    SbBool continuingPlan = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	const auto current =
	    this->p->residentMeshConsumerDemands.find(consumerId);
	if (current != this->p->residentMeshConsumerDemands.end() &&
	    current->second.revision == demandRevision) {
	    if (current->second.planning) {
		continuingPlan = TRUE;
	    } else if (current->second.residentRevision != 0 &&
		current->second.residentRevision ==
		    residentRevisionAtEntry) {
		if (planningComplete)
		    *planningComplete = TRUE;
		return 0;
	    }
	}
    }

    BObolResidentMeshConsumerDemand snapshot;
    if (!continuingPlan) {
	snapshot.revision = demandRevision;
	for (const BObolLodResidentDemand &demand : demands) {
	    if (demand.assetKey.getLength() == 0 || demand.cut < 0)
		continue;
	    BObolResidentMeshDemandValue &value =
		snapshot.assets[demand.assetKey.getString()];
	    value.cut = std::max(value.cut, demand.cut);
	    value.channelMask |= demand.channelMask & 3u;
	}
	snapshot.planning = TRUE;
	snapshot.planningCursor = 0;
	snapshot.planningProjectedResidentBytes =
	    lod_resident_stable_bytes(this->p);
	snapshot.planningResidentRevision = residentRevisionAtEntry;
    }

    size_t queued = 0;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	BObolResidentMeshConsumerDemand &current =
	    this->p->residentMeshConsumerDemands[consumerId];
	if (demandRevision < current.revision)
	    return 0;
	if (!continuingPlan || current.revision != demandRevision ||
	    !current.planning) {
	    current = std::move(snapshot);
	    current.residentRevision = 0;
	}

	/* Record the complete demand while refinement is active, but wait to
	 * queue trims.  A subsequent stable pump observes the changed resident
	 * revision and constructs the plan. */
	if (!this->p->pending.empty() || this->p->inFlight != 0)
	    return 0;

	/* Planning itself is owner-thread bookkeeping, so bound it just like
	 * result publication.  The append-only resident order makes the cursor
	 * stable while normal workers add later assets. */
	static const size_t planningQuantum = 2048;
	const size_t begin = current.planningCursor;
	const size_t end = std::min(
	    this->p->residentMeshOrder.size(),
	    begin + planningQuantum);
	for (size_t i = begin; i < end; ++i) {
	    const auto &residentEntry = this->p->residentMeshOrder[i];
	    const std::shared_ptr<BObolResidentMeshAsset> &resident =
		residentEntry.second;
	    if (!resident)
		continue;
	    const int residentCut =
		resident->publishedResidentCut.load(
		    std::memory_order_relaxed);
	    const int minimumCut =
		resident->publishedMinimumCut.load(
		    std::memory_order_relaxed);
	    const bool hasReloadableBacking =
		resident->publishedBackingPrefixBytes.load(
		    std::memory_order_relaxed) > 0;
	    if (residentCut < 0 || minimumCut < 0)
		continue;
	    BObolResidentMeshDemandValue aggregate;
	    SbBool demanded = FALSE;
	    for (const auto &consumer :
		    this->p->residentMeshConsumerDemands) {
		const auto demand =
		    consumer.second.assets.find(residentEntry.first);
		if (demand == consumer.second.assets.end())
		    continue;
		demanded = TRUE;
		aggregate.cut = std::max(
		    aggregate.cut, demand->second.cut);
		aggregate.channelMask |=
		    demand->second.channelMask;
	    }
	    const SbBool evict =
		!demanded &&
		this->p->maxResidentMeshBytes != SIZE_MAX &&
		current.planningProjectedResidentBytes >
		    this->p->maxResidentMeshBytes ? TRUE : FALSE;
	    const int targetCut = aggregate.cut < 0 ?
		minimumCut : std::max(minimumCut, aggregate.cut);
	    if (!evict &&
		targetCut >= residentCut && !hasReloadableBacking) {
		this->p->residentMeshCompactionTargets.erase(
		    residentEntry.first);
		continue;
	    }
	    BObolResidentMeshCompactionTarget &target =
		this->p->residentMeshCompactionTargets[
		    residentEntry.first];
	    target.cut = targetCut;
	    target.channelMask = aggregate.channelMask;
	    target.evict = evict;
	    target.useRevision = resident->useRevision.load(
		std::memory_order_relaxed);
	    target.revision =
		this->p->nextResidentMeshCompactionTargetRevision++;
	    if (!target.revision)
		target.revision =
		    this->p->nextResidentMeshCompactionTargetRevision++;
	    if (!this->p->residentMeshCompactionQueuedAssets.insert(
		    residentEntry.first).second)
		continue;
	    if (evict) {
		const size_t bytes = resident->publishedBytes.load(
		    std::memory_order_relaxed);
		const size_t backing =
		    resident->publishedBackingPrefixBytes.load(
			std::memory_order_relaxed);
		const size_t stableBytes =
		    backing >= bytes ? 0 : bytes - backing;
		current.planningProjectedResidentBytes =
		    stableBytes >=
			    current.planningProjectedResidentBytes ?
			0 : current.planningProjectedResidentBytes -
			    stableBytes;
	    }
	    BObolResidentMeshCompactionWork work;
	    work.assetKey = residentEntry.first;
	    work.resident = resident;
	    this->p->residentMeshCompactionWork.push_back(
		std::move(work));
	    queued++;
	}
	current.planningCursor = end;
	if (end >= this->p->residentMeshOrder.size()) {
	    current.planning = FALSE;
	    current.residentRevision = current.planningResidentRevision;
	    if (planningComplete)
		*planningComplete = TRUE;
	}
    }
    if (queued)
	this->p->workerCv.notify_all();
    return queued;
}

size_t
BObolLodService::drainResidentMeshCompactions(
    uint64_t consumerId,
    std::vector<BObolLodResidentCompaction> &results,
    size_t maxResults)
{
    if (!consumerId)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const auto found =
	this->p->residentMeshCompactionResults.find(consumerId);
    if (found == this->p->residentMeshCompactionResults.end())
	return 0;
    std::deque<BObolLodResidentCompaction> &queued = found->second;
    const size_t count = maxResults ?
	std::min(maxResults, queued.size()) : queued.size();
    results.reserve(results.size() + count);
    for (size_t i = 0; i < count; ++i) {
	results.push_back(std::move(queued.front()));
	queued.pop_front();
    }
    this->p->residentMeshCompactionResultCount =
	count >= this->p->residentMeshCompactionResultCount ?
	0 : this->p->residentMeshCompactionResultCount - count;
    if (queued.empty())
	this->p->residentMeshCompactionResults.erase(found);
    this->p->workerCv.notify_all();
    return count;
}

void
BObolLodService::noteResidentMeshUse(const BObolLodCacheKey &assetKey)
{
    if (!assetKey.isValid())
	return;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const auto found = this->p->residentMeshes.find(
	assetKey.value.getString());
    if (found == this->p->residentMeshes.end() || !found->second)
	return;
    found->second->useRevision.fetch_add(1, std::memory_order_relaxed);
}

void
BObolLodService::releaseResidentMeshConsumer(uint64_t consumerId)
{
    if (!consumerId)
	return;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    this->p->residentMeshConsumerDemands.erase(consumerId);
    const auto results =
	this->p->residentMeshCompactionResults.find(consumerId);
    if (results != this->p->residentMeshCompactionResults.end()) {
	const size_t count = results->second.size();
	this->p->residentMeshCompactionResultCount =
	    count >= this->p->residentMeshCompactionResultCount ?
	    0 : this->p->residentMeshCompactionResultCount - count;
	this->p->residentMeshCompactionResults.erase(results);
    }
    this->p->workerCv.notify_all();
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
		    lod_generation_count_remove_unlocked(
			this->p->generationPendingTaskCounts, generation);
		    lod_generation_task_finished_unlocked(this->p, generation);
	    if (it->task.publishResult && this->p->resultReservations > 0)
		this->p->resultReservations--;
	    if (this->p->cacheWriterEnabled && it->task.writeCache &&
		it->task.cacheWrite && this->p->cacheWriteReservations > 0)
		this->p->cacheWriteReservations--;
	    lod_active_request_key_remove_unlocked(this->p,
		lod_request_active_key(it->task.request));
	    lod_pending_quality_remove(
		this->p, it->task.request.qualityTier);
	    cancelled.push_back(std::move(*it));
	    it = this->p->pending.erase(it);
	}
	lod_prune_cancelled_generations_unlocked(this->p);

	for (std::list<BObolLodResult>::iterator it =
		 this->p->results.begin(); it != this->p->results.end();) {
	    if (it->generation == generation) {
		lod_queued_result_request_key_remove_unlocked(
		    this->p, it->request);
		this->p->resultSlots.erase(lod_result_slot_map_key(*it));
		lod_generation_count_remove_unlocked(
		    this->p->generationResultCounts, generation);
		it = this->p->results.erase(it);
	    } else {
		++it;
	    }
	}
	for (std::list<BObolLodCacheWriteItem>::iterator it =
		 this->p->cacheWrites.begin(); it != this->p->cacheWrites.end();) {
	    if (it->result.generation == generation) {
		this->p->cacheWriteSlots.erase(
		    lod_result_slot_map_key(it->result));
		lod_generation_count_remove_unlocked(
		    this->p->generationCacheWriteCounts, generation);
		it = this->p->cacheWrites.erase(it);
	    } else {
		++it;
	    }
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

void
BObolLodService::setResidentMeshLimit(size_t maxResidentBytes)
{
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	const size_t prior = this->p->maxResidentMeshBytes;
	const size_t next = maxResidentBytes > 0 ?
	    maxResidentBytes : lod_default_resident_mesh_limit();
	this->p->maxResidentMeshBytes = next;
	if (next > prior)
	    lod_resident_mesh_revision_advance(
		this->p->residentMeshAdmissionRevision);
    }
    this->p->workerCv.notify_all();
}

size_t
BObolLodService::getResidentMeshLimit(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->maxResidentMeshBytes;
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
lod_service_submit_task_unlocked(BObolLodServicePrivate *p,
				 const BObolLodTask &task,
				 const SbString &activeKey,
				 SbBool skipActiveDuplicate)
{
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

    if (skipActiveDuplicate &&
	lod_request_key_recorded_unlocked(p, activeKey))
	return 0;

    const uint64_t id = item.id;
    const uint64_t generation = item.task.generation;
    const SbBool publishResult = item.task.publishResult;
    const SbBool reserveCacheWrite =
	p->cacheWriterEnabled && item.task.writeCache && item.task.cacheWrite;
    const int qualityTier = item.task.request.qualityTier;
    p->pending.push_back(std::move(item));
    p->pendingQualityCounts[qualityTier]++;
    p->activeRequestKeyCounts[activeKey.getString()]++;
    p->generationTaskCounts[generation]++;
    p->generationPendingTaskCounts[generation]++;
    p->taskGenerations[id] = generation;
    if (publishResult)
	p->resultReservations++;
    if (reserveCacheWrite)
	p->cacheWriteReservations++;
    p->inFlight++;
    return id;
}

static uint64_t
lod_service_submit_task(BObolLodServicePrivate *p,
			const BObolLodTask &task,
			SbBool skipActiveDuplicate)
{
    const SbString activeKey = lod_request_active_key(task.request);
    uint64_t id = 0;
    {
	std::lock_guard<std::mutex> lock(p->mutex);
	id = lod_service_submit_task_unlocked(
	    p, task, activeKey, skipActiveDuplicate);
    }

    if (id)
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

size_t
BObolLodService::submitBatch(
    const std::vector<BObolLodTask> &tasks,
    std::vector<uint64_t> &taskIds,
    SbBool skipActiveDuplicates)
{
    taskIds.assign(tasks.size(), 0);
    if (tasks.empty())
	return 0;

    /*
     * Stable request keys are pure task data.  Build them before taking the
     * shared queue lock so workers can keep completing prior requests while
     * the producer prepares this wave.
     */
    std::vector<SbString> activeKeys;
    activeKeys.reserve(tasks.size());
    for (const BObolLodTask &task : tasks)
	activeKeys.push_back(lod_request_active_key(task.request));

    size_t accepted = 0;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->activeRequestKeyCounts.reserve(
	    this->p->activeRequestKeyCounts.size() + tasks.size());
	this->p->taskGenerations.reserve(
	    this->p->taskGenerations.size() + tasks.size());
	for (size_t i = 0; i < tasks.size(); ++i) {
	    taskIds[i] = lod_service_submit_task_unlocked(
		this->p, tasks[i], activeKeys[i],
		skipActiveDuplicates ? TRUE : FALSE);
	    if (taskIds[i])
		++accepted;
	}
    }

    if (accepted)
	this->p->workerCv.notify_all();
    return accepted;
}

SbBool
BObolLodService::hasActiveRequest(
    const BObolLodRequest &request) const
{
    const SbString key = lod_request_active_key(request);

    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_request_key_recorded_unlocked(this->p, key);
}

static size_t
lod_result_estimated_presentation_bytes(const BObolLodResult &result)
{
    if (result.resultKind != BOBOL_LOD_RESULT_MESH &&
	result.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL)
	return 0;
    const auto saturatingAdd = [](size_t left, size_t right) {
	return right > SIZE_MAX - left ? SIZE_MAX : left + right;
    };
    const auto saturatingMultiply = [](size_t left, size_t right) {
	return left && right > SIZE_MAX / left ? SIZE_MAX : left * right;
    };
    size_t bytes = 0;
    bytes = saturatingAdd(bytes,
	saturatingMultiply(result.mesh.points.size(), sizeof(SbVec3f)));
    bytes = saturatingAdd(bytes,
	saturatingMultiply(result.mesh.normals.size(), sizeof(SbVec3f)));
    bytes = saturatingAdd(bytes,
	saturatingMultiply(result.mesh.coordIndex.size(), sizeof(int32_t)));
    bytes = saturatingAdd(bytes,
	saturatingMultiply(result.mesh.faceIndex.size(), sizeof(int32_t)));
    bytes = saturatingAdd(bytes,
	saturatingMultiply(result.mesh.vertexIndex.size(), sizeof(int32_t)));
    size_t counted = 0;
    counted = saturatingAdd(counted,
	saturatingMultiply(static_cast<size_t>(
	    std::min<uint64_t>(result.counts.pointCount, SIZE_MAX)),
	    sizeof(SbVec3f)));
    counted = saturatingAdd(counted,
	saturatingMultiply(static_cast<size_t>(
	    std::min<uint64_t>(result.counts.faceCount, SIZE_MAX)),
	    3u * sizeof(int32_t)));
    counted = saturatingAdd(counted,
	saturatingMultiply(static_cast<size_t>(
	    std::min<uint64_t>(result.counts.normalCount, SIZE_MAX)),
	    sizeof(SbVec3f)));
    counted = std::max(counted, static_cast<size_t>(
	std::min<uint64_t>(result.counts.byteCount, SIZE_MAX)));
    return std::max(bytes, counted);
}

size_t
BObolLodService::drainResults(std::vector<BObolLodResult> &results,
				size_t maxResults,
				size_t maxEstimatedBytes)
{

    size_t count = 0;
    size_t estimatedBytes = 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);

    while (!this->p->results.empty() &&
	   (maxResults == 0 || count < maxResults)) {
	const size_t frontBytes = lod_result_estimated_presentation_bytes(
	    this->p->results.front());
	if (count && maxEstimatedBytes &&
	    (estimatedBytes >= maxEstimatedBytes ||
	     frontBytes > maxEstimatedBytes - estimatedBytes))
	    break;
	lod_queued_result_request_key_remove_unlocked(
	    this->p, this->p->results.front().request);
	this->p->resultSlots.erase(
	    lod_result_slot_map_key(this->p->results.front()));
	lod_generation_count_remove_unlocked(
	    this->p->generationResultCounts,
	    this->p->results.front().generation);
	results.push_back(std::move(this->p->results.front()));
	this->p->results.pop_front();
	estimatedBytes = frontBytes > SIZE_MAX - estimatedBytes ?
	    SIZE_MAX : estimatedBytes + frontBytes;
	count++;
    }

    return count;
}

size_t
BObolLodService::drainGenerationResults(
    std::vector<BObolLodResult> &results, uint64_t generation,
    size_t maxResults, size_t maxEstimatedBytes)
{
    if (generation == 0)
	return 0;

    size_t count = 0;
    size_t estimatedBytes = 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);

    for (std::list<BObolLodResult>::iterator it =
	     this->p->results.begin(); it != this->p->results.end();) {
	if (maxResults != 0 && count >= maxResults)
	    break;
	if (it->generation != generation) {
	    ++it;
	    continue;
	}

	const size_t resultBytes =
	    lod_result_estimated_presentation_bytes(*it);
	if (count && maxEstimatedBytes &&
	    (estimatedBytes >= maxEstimatedBytes ||
	     resultBytes > maxEstimatedBytes - estimatedBytes))
	    break;

	lod_queued_result_request_key_remove_unlocked(
	    this->p, it->request);
	this->p->resultSlots.erase(lod_result_slot_map_key(*it));
	lod_generation_count_remove_unlocked(
	    this->p->generationResultCounts, generation);
	results.push_back(std::move(*it));
	it = this->p->results.erase(it);
	estimatedBytes = resultBytes > SIZE_MAX - estimatedBytes ?
	    SIZE_MAX : estimatedBytes + resultBytes;
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

    for (std::list<BObolLodResult>::iterator it =
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

	lod_queued_result_request_key_remove_unlocked(
	    this->p, it->request);
	this->p->resultSlots.erase(lod_result_slot_map_key(*it));
	lod_generation_count_remove_unlocked(
	    this->p->generationResultCounts, it->generation);
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
BObolLodService::resultReservationCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->resultReservations;
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

size_t
BObolLodService::activeTaskCountForGeneration(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_count_unlocked(
	this->p->generationTaskCounts, generation);
}

size_t
BObolLodService::pendingTaskCountForGeneration(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_count_unlocked(
	this->p->generationPendingTaskCounts, generation);
}

size_t
BObolLodService::executingTaskCountForGeneration(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_count_unlocked(
	this->p->generationExecutingTaskCounts, generation);
}

size_t
BObolLodService::queuedResultCountForGeneration(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_count_unlocked(
	this->p->generationResultCounts, generation);
}

size_t
BObolLodService::queuedCacheWriteCountForGeneration(
    uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_count_unlocked(
	this->p->generationCacheWriteCounts, generation);
}

size_t
BObolLodService::delayedTaskCountForGeneration(uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_generation_count_unlocked(
	this->p->generationDelayedTaskCounts, generation);
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
    return this->p->residentMeshBytes.load(std::memory_order_relaxed);
}

size_t
BObolLodService::stableResidentMeshBytesForDiagnostics(void) const
{
    return lod_resident_stable_bytes(this->p);
}

size_t
BObolLodService::reservedResidentMeshGrowthBytesForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshGrowthReservationBytes;
}

uint64_t
BObolLodService::residentMeshAdmissionRevision(void) const
{
    return this->p->residentMeshAdmissionRevision.load(
	std::memory_order_relaxed);
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

uint64_t
BObolLodService::residentMeshEvictionCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshEvictions;
}

size_t
BObolLodService::pendingResidentMeshCompactionCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    size_t planning = 0;
    for (const auto &consumer :
	    this->p->residentMeshConsumerDemands)
	if (consumer.second.planning)
	    planning++;
    return this->p->residentMeshCompactionWork.size() +
	this->p->residentMeshCompactionsInFlight + planning;
}

size_t
BObolLodService::queuedResidentMeshCompactionResultCountForDiagnostics(
    uint64_t consumerId) const
{
    if (!consumerId)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const auto found =
	this->p->residentMeshCompactionResults.find(consumerId);
    return found == this->p->residentMeshCompactionResults.end() ?
	0 : found->second.size();
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
