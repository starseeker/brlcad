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
#include "cad_publication_private.h"
#include "draw_cache_private.h"
#include "identity_counter_private.h"
#include "lod_coverage_preview_private.h"
#include "lod_coordinator_private.h"
#include "parallel_budget_private.h"

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

#if !defined(_WIN32)
#  include <sys/resource.h>
#  include <unistd.h>
#endif

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
    generation = 0;
    databaseLease.reset();
    stagedSource.reset();
    useSerializedSpatialSource = FALSE;
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
    initialRefinementCostBudget = 500000;
    refinementGrowthFactor = 4.0;
    useCurrentDrawCut = FALSE;
    currentDrawCut = -1;
    useDeliveryCutLimit = FALSE;
    deliveryCutLimit = -1;
    transientMemoryLimited = FALSE;
    usePresentationCutLimit = FALSE;
    presentationCutLimit = -1;
    presentationAdmissionCertified = FALSE;
    presentationAdmissionViewRevision = 0;
    presentationAdmissionPolicyRevision = 0;
    atomicRepresentationHandoff = FALSE;
    forcedCut = 0;
    resetExisting = FALSE;
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
    dispatchClass = BOBOL_LOD_TASK_DISPATCH_NORMAL;
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
    VSUB2(yaxis, points[2], points[0]);
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

static SbBool
lod_provider_param_enabled(const BObolLodRequest &request, const char *name)
{
    const BObolLodProviderParam *param = lod_provider_param(request, name);
    return param && BU_STR_EQUAL(param->value.getString(), "1") ?
	TRUE : FALSE;
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
    if (!useQueryBounds && !useQueryRay &&
	lod_provider_param_enabled(request,
	    BOBOL_LOD_PREPARED_CAD_ONLY_PARAM)) {
	result.preparedCadGeometry =
	    bobol_database_bot_part_geometry(bot, request.drawMode);
	if (!result.preparedCadGeometry) {
	    result.geometry.clear();
	    result.resultKind = BOBOL_LOD_RESULT_NONE;
	    result.providerStatus = BOBOL_LOD_PROVIDER_ERROR;
	    result.diagnostic =
		"RT source full-detail provider could not prepare CAD geometry";
	    return result;
	}
	result.preparedCadGeometryRevision = 1;
	result.geometry.cacheKey = bobol_lod_geometry_cache_key(request);
	result.cacheKey = bobol_lod_cache_key(request);
	result.bounds = request.bounds;
	result.counts.faceCount = bot->num_faces;
	result.counts.originalPointCount = bot->num_vertices;
	if (result.preparedCadGeometry->shaded) {
	    const Obol::TriMesh &mesh =
		*result.preparedCadGeometry->shaded;
	    result.bounds = mesh.bounds;
	    result.counts.pointCount = mesh.positions.size();
	    result.counts.normalCount = mesh.normals.size();
	    result.hasNormals = mesh.normals.empty() ? FALSE : TRUE;
	} else if (result.preparedCadGeometry->wire) {
	    const Obol::WireRep &wire = *result.preparedCadGeometry->wire;
	    result.bounds = wire.bounds;
	    result.counts.pointCount = bot->num_vertices;
	    result.counts.lineCount = wire.segmentCount();
	}
	result.counts.byteCount =
	    bobol_database_part_geometry_estimate_bytes(
		*result.preparedCadGeometry);
	result.shadedCullBackfaces =
	    result.preparedCadGeometry->shadedCullBackfaces ? TRUE : FALSE;
	return result;
    }
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
	!result.preparedCadGeometry &&
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
	lod, requestedCut, provider->resetExisting ? 1 : 0);
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
    int64_t submittedMicroseconds = 0;
    /* Only submitIfNotActive work has a mutable demand record.  Plain submit
     * may intentionally enqueue several stages with one stable key, and must
     * continue executing each stage against its own request. */
    bool retargetActiveDemand = false;
    /* Service-local admission is reserved when the worker dequeues the task.
     * Preserve the exact charge so a later limit change or an oversized task
     * cannot subtract bytes which this item never reserved. */
    size_t reservedWorkingSetBytes = 0;
    bool exceedsServiceWorkingSetLimit = false;
};

struct BObolLodCacheWriteItem {
    BObolLodResult result;
    BObolLodCacheWriteProc write;
    void *writeData;
};

enum class BObolSharedProducerState : uint8_t {
    BUILDING = 0,
    RESULT
};

struct BObolSharedProducerLease {
    BObolLodRequest demand;
    uint64_t deliveredPublication = 0;
};

struct BObolSharedProducer {
    uint64_t taskId = 0;
    uint64_t producerGeneration = 0;
    uint64_t publication = 0;
    int64_t submittedMicroseconds = 0;
    int64_t executionStartedMicroseconds = 0;
    bool cpuAdmissionWaiting = false;
    bool transientMemoryAdmissionWaiting = false;
    BObolSharedProducerState state = BObolSharedProducerState::BUILDING;
    std::unordered_map<uint64_t, BObolSharedProducerLease> leases;
};

/* Worker-owned progress is published without taking the service mutex.  A
 * cold mesh producer holds its resident-asset mutex for most of realization,
 * so acquiring the outer service lock from a cache callback would invert the
 * service -> resident lock order.  The short record-local transition mutex
 * makes each stage/timing snapshot coherent without entering that hierarchy. */
struct BObolLodProducerProgressRecord {
    uint64_t ownerGeneration = 0;
    std::string activeKey;
    int64_t taskStartedMicroseconds = 0;
    uint64_t queueWaitMicroseconds = 0;
    uint64_t sourceFaceCount = 0;
    uint64_t sourcePointCount = 0;
    uint64_t sourceByteCount = 0;
    std::atomic<bool> active {true};
    std::atomic<int> stage {BOBOL_LOD_PRODUCER_STAGE_NONE};
    std::atomic<uint64_t> completedUnits {0};
    std::atomic<uint64_t> totalUnits {0};
    mutable std::mutex transitionMutex;
    int64_t stageStartedMicroseconds = 0;

    struct Snapshot {
	SbBool active = FALSE;
	int stage = BOBOL_LOD_PRODUCER_STAGE_NONE;
	uint64_t completedUnits = 0;
	uint64_t totalUnits = 0;
	uint64_t queueWaitMicroseconds = 0;
	uint64_t elapsedMicroseconds = 0;
	uint64_t stageElapsedMicroseconds = 0;
	uint64_t sourceFaceCount = 0;
	uint64_t sourcePointCount = 0;
	uint64_t sourceByteCount = 0;
    };

    static uint64_t elapsed(int64_t started, int64_t now)
    {
	return started > 0 && now > started ?
	    static_cast<uint64_t>(now - started) : 0;
    }

    void publish(int nextStage, uint64_t completed, uint64_t total)
    {
	if (nextStage <= BOBOL_LOD_PRODUCER_STAGE_NONE ||
	    nextStage >= BOBOL_LOD_PRODUCER_STAGE_COUNT)
	    return;
	if (total && completed > total)
	    completed = total;
	std::lock_guard<std::mutex> lock(this->transitionMutex);
	const int current = stage.load(std::memory_order_acquire);
	if (current != nextStage) {
	    completedUnits.store(completed, std::memory_order_relaxed);
	    totalUnits.store(total, std::memory_order_relaxed);
	    stageStartedMicroseconds = bu_gettime();
	    stage.store(nextStage, std::memory_order_release);
	    return;
	}
	/* A bounded parallel classifier may finish ranges out of order.  Never
	 * let a late callback move the visible counter backwards. */
	totalUnits.store(total, std::memory_order_relaxed);
	uint64_t observed = completedUnits.load(std::memory_order_relaxed);
	while (observed < completed &&
	       !completedUnits.compare_exchange_weak(observed, completed,
		   std::memory_order_release, std::memory_order_relaxed)) {
	}
    }

    void finish(void)
    {
	std::lock_guard<std::mutex> lock(this->transitionMutex);
	active.store(false, std::memory_order_release);
    }

    Snapshot snapshot(int64_t now) const
    {
	Snapshot result;
	std::lock_guard<std::mutex> lock(this->transitionMutex);
	result.active = active.load(std::memory_order_acquire) ? TRUE : FALSE;
	result.stage = stage.load(std::memory_order_acquire);
	result.totalUnits = totalUnits.load(std::memory_order_relaxed);
	result.completedUnits = std::min(
	    result.totalUnits ? result.totalUnits : UINT64_MAX,
	    completedUnits.load(std::memory_order_acquire));
	result.queueWaitMicroseconds = queueWaitMicroseconds;
	result.elapsedMicroseconds = elapsed(taskStartedMicroseconds, now);
	result.stageElapsedMicroseconds = elapsed(stageStartedMicroseconds, now);
	result.sourceFaceCount = sourceFaceCount;
	result.sourcePointCount = sourcePointCount;
	result.sourceByteCount = sourceByteCount;
	return result;
    }
};

class BObolLodProducerProgressScope {
public:
    explicit BObolLodProducerProgressScope(
	const std::shared_ptr<BObolLodProducerProgressRecord> &value) :
	record(value)
    {
    }

    ~BObolLodProducerProgressScope()
    {
	if (record)
	    record->finish();
    }

private:
    std::shared_ptr<BObolLodProducerProgressRecord> record;
};

struct BObolLodResultSlotMapKey {
    std::string databaseId;
    std::string occurrence;
    std::string providerId;
    uint64_t generation = 0;
    uint64_t sourceRoutingId = 0;
    uint64_t sourcePopulationEpoch = 0;
    int drawMode = 0;
    int resultKind = 0;
    int proxyKind = 0;

    bool operator==(const BObolLodResultSlotMapKey &other) const
    {
	return generation == other.generation &&
	    sourceRoutingId == other.sourceRoutingId &&
	    sourcePopulationEpoch == other.sourcePopulationEpoch &&
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
	combine(std::hash<uint64_t>()(key.sourcePopulationEpoch));
	combine(std::hash<int>()(key.drawMode));
	combine(std::hash<int>()(key.resultKind));
	combine(std::hash<int>()(key.proxyKind));
	return hash;
    }
};

/* An intermediate publication is owned by one active producer, regardless of
 * whether successive callbacks describe temporary coverage, a growing
 * spatial page set, or a global mesh prefix.  Keeping one deferred value per
 * producer bounds lock-contention recovery without retaining every transient
 * geometry snapshot. */
struct BObolDeferredIntermediateKey {
    uint64_t generation = 0;
    std::string activeKey;

    bool operator==(const BObolDeferredIntermediateKey &other) const
    {
	return generation == other.generation && activeKey == other.activeKey;
    }
};

struct BObolDeferredIntermediateKeyHash {
    size_t operator()(const BObolDeferredIntermediateKey &key) const
    {
	size_t hash = std::hash<std::string>()(key.activeKey);
	hash ^= std::hash<uint64_t>()(key.generation) +
	    static_cast<size_t>(0x9e3779b9U) + (hash << 6) + (hash >> 2);
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
    key.sourcePopulationEpoch = request.sourcePopulationEpoch.value();
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
	    std::max(floor, totalBytes / 8));
    size_t availableBytes = 0;
    if (bu_mem(BU_MEM_AVAIL, &availableBytes) >= 0 && availableBytes > 0)
	allowance = std::min(allowance,
	    std::max(floor, availableBytes / 3));
    /* This is the durable resident ceiling, not the transient admission
     * governor above.  Installed capacity supplies the normal share, while
     * the startup available-memory sample prevents concurrent viewers from
     * each claiming that same share.  Concurrent topology work remains
     * bounded by the separate working-set governor.  The conservative
     * one-eighth-RAM capacity share accounts for the fact that one logical
     * prefix may also be referenced by a prepared renderer record and a GL
     * driver allocation.  The available-memory cap prevents a second viewer
     * or another application from independently spending the same physical
     * headroom. */
    return std::max(floor, allowance);
}

/* A retained prefix may coexist with its renderer allocation and, on OSMesa,
 * software graphics storage in the same physical-memory pool.  Reserving at
 * least half of the available-memory snapshot for those consumers is the
 * highest explicit user override we can offer without presenting an OOM knob. */
static constexpr double LOD_MAX_RESIDENT_AVAILABLE_MEMORY_PERCENT = 50.0;

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

SbBool
bobol_lod_working_set_acquire(size_t estimatedBytes)
{
    if (!estimatedBytes)
	return TRUE;
    BObolGlobalWorkingSetGovernor &governor =
	lod_global_working_set_governor();
    std::unique_lock<std::mutex> lock(governor.mutex);
    if (governor.limit != SIZE_MAX && estimatedBytes > governor.limit)
	return FALSE;
    governor.cv.wait(lock, [&]() {
	if (governor.limit == SIZE_MAX)
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
    return TRUE;
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
	orderIndex(SIZE_MAX),
	orientedBoundsPublished(false)
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
    /* A source-limited spatial build owns complete whole-object coverage plus
     * any locally validated pages here.  The result queue only borrows shared
     * immutable geometry; retaining ownership at the asset level prevents a
     * later view request from falling back to an empty page set. */
    std::vector<BObolLodPresentationLayer> limitedSpatialLayers;
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
    /* The draw-asset record is the coverage path's O(1) metadata carrier.
     * Publish the hierarchy OBB there once; never make discovery reopen the
     * much larger PoP payload merely to improve a terminal proxy. */
    bool orientedBoundsPublished;
};

struct BObolResidentMeshDemandValue {
    int cut = -1;
    unsigned int channelMask = 0;
    std::vector<uint32_t> chunkIds;

    bool operator==(const BObolResidentMeshDemandValue &other) const
    {
	return cut == other.cut && channelMask == other.channelMask &&
	    chunkIds == other.chunkIds;
    }
};

struct BObolResidentMeshConsumerDemand {
    uint64_t revision = 0;
    /* assets describe revision only while these values match.  Separating
     * observed revision from snapshot revision makes invalidation O(1)
     * without pretending the prior asset map is current. */
    uint64_t snapshotRevision = 0;
    // Resident-asset mutation epoch for which this snapshot was fully
    // compacted.  Zero means the snapshot was recorded while workers were
    // active and must be retried.
    uint64_t residentRevision = 0;
    std::unordered_map<std::string, BObolResidentMeshDemandValue> assets;
    size_t planningCursor = 0;
    size_t planningProjectedResidentBytes = 0;
    size_t planningCandidateCount = 0;
    size_t completedCandidateCount = 0;
    uint64_t completedPlanRevision = 0;
    uint64_t planningResidentRevision = 0;
    SbBool planning = FALSE;
};

struct BObolResidentMeshCompactionTarget {
    int cut = -1;
    unsigned int channelMask = 0;
    /* Exact stable working set for a chunked asset.  Empty retains the
     * ordinary unchunked-prefix interpretation of cut. */
    std::vector<BObolLodChunkCut> chunkCuts;
    SbBool evict = FALSE;
    uint64_t useRevision = 0;
    uint64_t demandEpoch = 0;
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
 *   3. BObolLodServicePrivate::residentMeshAdmissionMutex (short stable-byte
 *      reservation accounting)
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
 *   - deferred intermediate results are protected by
 *     deferredIntermediateMutex.  No code may acquire the service mutex while
 *     holding that mailbox lock; promotion takes service then mailbox;
 *   - realization callbacks and task-local payloads are worker owned;
 *   - diagnostic byte/count summaries used without the service lock are
 *     atomic;
 *   - resident growth reservations are protected by the admission mutex.
 *     It is a leaf lock: code holding it must not acquire either the service
 *     or a resident-asset mutex;
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
	maxDeferredIntermediateResults(1),
	maxQueuedCacheWrites(2048),
	maxActiveWorkingSetBytes(0),
	maxResidentMeshBytes(0),
	residentMeshLimitPercent(0.0),
	residentMeshLimitBasisBytes(0),
	activeWorkingSetBytes(0),
	executingTasks(0),
	cpuAdmissionWaitingTasks(0),
	transientMemoryAdmissionWaitingTasks(0),
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
	residentMeshStableBytes(0),
	residentMeshGrowthReservationBytes(0),
	residentMeshRevision(1),
	residentMeshAdmissionRevision(1),
	deferredIntermediateResultCount(0),
	deferredIntermediateNotificationPending(false),
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
    mutable std::mutex residentMeshAdmissionMutex;
    mutable std::mutex deferredIntermediateMutex;
    mutable std::mutex producerProgressMutex;
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
    std::atomic<size_t> maxDeferredIntermediateResults;
    size_t maxQueuedCacheWrites;
    size_t maxActiveWorkingSetBytes;
    size_t maxResidentMeshBytes;
    double residentMeshLimitPercent;
    size_t residentMeshLimitBasisBytes;
    size_t activeWorkingSetBytes;
    size_t executingTasks;
    size_t cpuAdmissionWaitingTasks;
    size_t transientMemoryAdmissionWaitingTasks;
    size_t peakWorkingSetBytes;
    size_t peakExecutingTasks;
    size_t resultReservations;
    size_t cacheWriteReservations;
    uint64_t rejectedTasks;
    uint64_t coalescedResults;
    uint64_t coalescedCacheWrites;
    uint64_t discardedStaleResults;
    std::atomic<uint64_t> residentMeshCacheLoads;
    std::atomic<uint64_t> residentMeshHits;
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
    /* Stable immutable renderer bytes are published as one scalar.  Deriving
     * this value from separate total/backing loads allowed a policy reader to
     * observe half of a concurrent accounting replacement. */
    std::atomic<size_t> residentMeshStableBytes;
    /* Protected by residentMeshAdmissionMutex.  Workers reserve optional
     * stable-prefix growth before loading so independent assets cannot all
     * observe the same free capacity.  Minimum useful prefixes are permitted
     * to exceed the soft target, but are still reserved and therefore
     * constrain richer peers. */
    size_t residentMeshGrowthReservationBytes;
    std::atomic<uint64_t> residentMeshRevision;
    std::atomic<uint64_t> residentMeshAdmissionRevision;
    std::deque<BObolLodWorkItem> pending;
    std::map<int, size_t> pendingDispatchCounts;
    std::map<int, size_t> pendingQualityCounts;
    std::list<BObolLodResult> results;
    std::unordered_map<BObolDeferredIntermediateKey, BObolLodResult,
	BObolDeferredIntermediateKeyHash> deferredIntermediateResults;
    std::atomic<size_t> deferredIntermediateResultCount;
    std::atomic<bool> deferredIntermediateNotificationPending;
    std::list<BObolLodCacheWriteItem> cacheWrites;
    std::unordered_map<BObolLodResultSlotMapKey,
	std::list<BObolLodResult>::iterator,
	BObolLodResultSlotMapKeyHash> resultSlots;
    std::unordered_map<BObolLodResultSlotMapKey,
	std::list<BObolLodCacheWriteItem>::iterator,
	BObolLodResultSlotMapKeyHash> cacheWriteSlots;
    std::vector<BObolLodSubscriber> subscribers;
    std::unordered_map<std::string, size_t> activeRequestKeyCounts;
    /* View and policy epochs are demand, not asset identity.  One coalesced
     * cold producer retains the newest request for its stable active key so
     * view-independent page publications do not become stale merely because
     * the camera moved while source preparation continued. */
    std::unordered_map<std::string, BObolLodRequest> latestActiveRequests;
    /* A shared producer has one immutable worker and one explicit lease per
     * interested generation.  Non-owner generations receive a lightweight
     * replay result after each producer publication; the owner alone may
     * receive the demand-specific payload it prepared. */
    std::unordered_map<std::string, BObolSharedProducer> sharedProducers;
    std::unordered_map<uint64_t, std::string> sharedProducerTaskKeys;
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
    uint64_t residentMeshDemandEpoch = 1;
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
    /* Exact enqueue timestamps make queue latency observable without an
     * O(pending task count) scan on every GUI progress sample. */
    std::unordered_map<uint64_t, std::multiset<int64_t>>
	generationPendingTaskTimes;
    std::unordered_map<uint64_t, size_t> generationExecutingTaskCounts;
    std::unordered_map<uint64_t, size_t>
	generationCpuAdmissionWaitingTaskCounts;
    std::unordered_map<uint64_t, size_t>
	generationTransientMemoryAdmissionWaitingTaskCounts;
    std::unordered_map<uint64_t, size_t> generationDelayedTaskCounts;
    std::unordered_map<uint64_t, size_t> generationResultCounts;
    std::unordered_map<uint64_t, size_t> generationCacheWriteCounts;
    uint64_t nextProducerProgressId = 1;
    std::unordered_map<uint64_t,
	std::shared_ptr<BObolLodProducerProgressRecord>> producerProgress;
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

static size_t
lod_resident_asset_stable_bytes(size_t total, size_t backing)
{
    return backing >= total ? 0 : total - backing;
}

static void
lod_resident_mesh_accounting_replace(BObolLodServicePrivate *service,
    size_t priorBytes, size_t currentBytes,
    size_t priorBackingBytes, size_t currentBackingBytes)
{
    if (!service)
	return;
    lod_resident_mesh_bytes_replace(
	service->residentMeshBytes, priorBytes, currentBytes);
    lod_resident_mesh_bytes_replace(
	service->residentMeshBackingBytes,
	priorBackingBytes, currentBackingBytes);
    lod_resident_mesh_bytes_replace(
	service->residentMeshStableBytes,
	lod_resident_asset_stable_bytes(priorBytes, priorBackingBytes),
	lod_resident_asset_stable_bytes(currentBytes, currentBackingBytes));
}

static void
lod_resident_mesh_revision_advance(std::atomic<uint64_t> &revision)
{
    (void)bobol_atomic_identity_advance(revision);
}

static bool
lod_resident_consumer_snapshot_current(
    const BObolResidentMeshConsumerDemand &consumer)
{
    return consumer.revision != 0 &&
	consumer.snapshotRevision == consumer.revision;
}

static void
lod_resident_demand_epoch_advance(BObolLodServicePrivate *p)
{
    if (!p)
	return;
    bobol_identity_advance(p->residentMeshDemandEpoch);
}

static void
lod_discard_resident_compaction_results_unlocked(
    BObolLodServicePrivate *p, uint64_t consumerId)
{
    if (!p || !consumerId)
	return;
    const auto found = p->residentMeshCompactionResults.find(consumerId);
    if (found == p->residentMeshCompactionResults.end())
	return;
    const size_t count = found->second.size();
    p->residentMeshCompactionResultCount =
	count >= p->residentMeshCompactionResultCount ?
	0 : p->residentMeshCompactionResultCount - count;
    p->residentMeshCompactionResults.erase(found);
}

static size_t
lod_resident_stable_bytes(const BObolLodServicePrivate *p)
{
    if (!p)
	return 0;
    return p->residentMeshStableBytes.load(std::memory_order_relaxed);
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
	SbBool &limited,
	const std::vector<uint32_t> *requiredChunks = NULL)
    {
	limited = FALSE;
	if (!p || desiredCut < hierarchy.min_cut)
	    return desiredCut;
	desiredCut = std::min(desiredCut, hierarchy.max_cut);
	if (publishedCut >= desiredCut) {
	    /* This is still an admission decision.  A transient working-set cap
	     * can deliberately present a poorer cut from an already-richer
	     * immutable resident mesh.  Its terminal result must carry the current
	     * capacity epoch, or the owner cannot distinguish a durable denial from
	     * an unstamped result and will resubmit the identical no-load task on
	     * every pump. */
	    decisionRevision = p->residentMeshAdmissionRevision.load(
		std::memory_order_relaxed);
	    return desiredCut;
	}

	const int floorCut = publishedCut >= hierarchy.min_cut ?
	    publishedCut : hierarchy.min_cut;
	int admittedCut = floorCut;
	const auto estimateAtCut = [&](int cut) -> size_t {
	    if (!requiredChunks || requiredChunks->empty() ||
		!hierarchy.chunks)
		return lod_resident_cut_stable_bytes(resident, hierarchy, cut);
	    size_t estimate = 0;
	    for (uint32_t chunk : *requiredChunks) {
		if (chunk >= hierarchy.chunk_count)
		    return SIZE_MAX;
		const uint64_t chunkBytes =
		    hierarchy.chunks[chunk].cuts[cut].resident_bytes;
		const size_t value = chunkBytes > SIZE_MAX ? SIZE_MAX :
		    static_cast<size_t>(chunkBytes);
		estimate = value > SIZE_MAX - estimate ? SIZE_MAX :
		    estimate + value;
	    }
	    return estimate;
	};
	size_t admittedEstimate =
	    publishedCut >= hierarchy.min_cut ?
		priorStableBytes :
		estimateAtCut(floorCut);

	/* Realization owns resident.mutex here.  Admission accounting has its own
	 * leaf lock so a worker never waits for the outer service lock in the
	 * inverse of the documented service -> resident order. */
	std::lock_guard<std::mutex> lock(p->residentMeshAdmissionMutex);
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
	    const size_t estimate = estimateAtCut(cut);
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
	{
	    std::lock_guard<std::mutex> lock(
		p->residentMeshAdmissionMutex);
	    p->residentMeshGrowthReservationBytes =
		bytes >= p->residentMeshGrowthReservationBytes ?
		    0 : p->residentMeshGrowthReservationBytes - bytes;
	    bytes = 0;
	}
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
    size_t bytes = backingBytes > SIZE_MAX - meshBytes ?
	SIZE_MAX : meshBytes + backingBytes;
    std::unordered_set<const Obol::PartGeometry *> countedGeometry;
    for (const BObolLodPresentationLayer &layer :
	 resident.limitedSpatialLayers) {
	if (!layer.geometry ||
	    !countedGeometry.emplace(layer.geometry.get()).second)
	    continue;
	const size_t geometryBytes =
	    bobol_database_part_geometry_estimate_bytes(*layer.geometry);
	bytes = geometryBytes > SIZE_MAX - bytes ?
	    SIZE_MAX : bytes + geometryBytes;
    }
    return bytes;
}

static SbBool
lod_generation_cancelled_unlocked(const BObolLodServicePrivate *p,
				  uint64_t generation)
{
    return p->cancelledGenerations.find(generation) !=
	   p->cancelledGenerations.end() ? TRUE : FALSE;
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
lod_generation_pending_time_add_unlocked(BObolLodServicePrivate *p,
	uint64_t generation, int64_t submittedMicroseconds)
{
    if (!p || !generation || submittedMicroseconds <= 0)
	return;
    p->generationPendingTaskTimes[generation].insert(
	submittedMicroseconds);
}

static void
lod_generation_pending_time_remove_unlocked(BObolLodServicePrivate *p,
	uint64_t generation, int64_t submittedMicroseconds)
{
    if (!p || !generation || submittedMicroseconds <= 0)
	return;
    auto generationTimes = p->generationPendingTaskTimes.find(generation);
    if (generationTimes == p->generationPendingTaskTimes.end())
	return;
    const auto timestamp = generationTimes->second.find(
	submittedMicroseconds);
    if (timestamp != generationTimes->second.end())
	generationTimes->second.erase(timestamp);
    if (generationTimes->second.empty())
	p->generationPendingTaskTimes.erase(generationTimes);
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
     * Compact occurrences remain distinct consumers.  For an explicitly
     * asset-coalesced request, the owner-thread planner binds siblings from
     * the published resident asset (or submits their still-missing spatial
     * pages) after this producer completes.  This serializes expensive cold
     * hierarchy construction without conflating occurrence presentation. */
    SbString key = bobol_lod_geometry_cache_key(request).value;
    if (request.occurrenceKey.getLength() > 0 &&
	!request.coalesceAssetProducer) {
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
    if (request.sourcePopulationEpoch != 0) {
	key += "|population=";
	key += SbString(std::to_string(
	    request.sourcePopulationEpoch.value()).c_str());
    }
    return key;
}

/* Set by the worker around one realize() call.  Direct provider invocations
 * leave these zero and are correctly reported as having no scheduler wait. */
static thread_local int64_t lod_task_submitted_microseconds = 0;
static thread_local int64_t lod_task_started_microseconds = 0;

static std::shared_ptr<BObolLodProducerProgressRecord>
lod_producer_progress_begin(BObolLodServicePrivate *p,
	uint64_t generation, const BObolLodRequest &request)
{
    if (!p || !generation)
	return std::shared_ptr<BObolLodProducerProgressRecord>();
    try {
	std::shared_ptr<BObolLodProducerProgressRecord> record =
	    std::make_shared<BObolLodProducerProgressRecord>();
	record->ownerGeneration = generation;
	record->activeKey = lod_request_active_key(request).getString();
	const int64_t now = std::max<int64_t>(1, bu_gettime());
	record->taskStartedMicroseconds =
	    lod_task_started_microseconds > 0 ?
		lod_task_started_microseconds : now;
	record->queueWaitMicroseconds =
	    lod_task_submitted_microseconds > 0 &&
	    record->taskStartedMicroseconds > lod_task_submitted_microseconds ?
		static_cast<uint64_t>(record->taskStartedMicroseconds -
		    lod_task_submitted_microseconds) : 0;
	record->sourceFaceCount = request.sourceCounts.faceCount;
	record->sourcePointCount = request.sourceCounts.pointCount;
	record->sourceByteCount = request.sourceCounts.byteCount;
	record->publish(BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION, 0, 0);
	std::lock_guard<std::mutex> lock(p->producerProgressMutex);
	const uint64_t id = bobol_nonzero_identity_take(
	    p->nextProducerProgressId);
	p->producerProgress.emplace(id, record);
	return record;
    } catch (const std::bad_alloc &) {
	/* Diagnostics must never make an otherwise admissible mesh fail. */
	return std::shared_ptr<BObolLodProducerProgressRecord>();
    }
}

/* Public stage values are stable diagnostics, not an execution-order
 * encoding.  Later stages were appended for ABI compatibility, so rank them
 * by execution order when choosing the one stage summarized by a progress
 * bar. */
static int
lod_producer_stage_display_priority(int stage)
{
    switch (stage) {
	case BOBOL_LOD_PRODUCER_STAGE_ASSET_SERIALIZATION: return 0;
	case BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP: return 1;
	case BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION: return 2;
	case BOBOL_LOD_PRODUCER_STAGE_BOUNDS_ANALYSIS: return 3;
	case BOBOL_LOD_PRODUCER_STAGE_COVERAGE_PREVIEW: return 4;
	case BOBOL_LOD_PRODUCER_STAGE_SOURCE_HASHING: return 5;
	case BOBOL_LOD_PRODUCER_STAGE_FACE_CLASSIFICATION: return 6;
	case BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION: return 7;
	case BOBOL_LOD_PRODUCER_STAGE_SPATIAL_CONSTRUCTION: return 8;
	case BOBOL_LOD_PRODUCER_STAGE_CACHE_PERSISTENCE: return 9;
	default: return INT_MAX;
    }
}

static void
lod_producer_progress_prune(BObolLodServicePrivate *p)
{
    if (!p)
	return;
    std::lock_guard<std::mutex> lock(p->producerProgressMutex);
    for (auto it = p->producerProgress.begin();
	 it != p->producerProgress.end();) {
	if (!it->second ||
	    !it->second->active.load(std::memory_order_acquire))
	    it = p->producerProgress.erase(it);
	else
	    ++it;
    }
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

static BObolSharedProducer *
lod_shared_producer_unlocked(BObolLodServicePrivate *p,
	const SbString &key)
{
    if (!p || key.getLength() == 0)
	return NULL;
    const auto found = p->sharedProducers.find(key.getString());
    return found == p->sharedProducers.end() ? NULL : &found->second;
}

static const BObolSharedProducer *
lod_shared_producer_unlocked(const BObolLodServicePrivate *p,
	const SbString &key)
{
    if (!p || key.getLength() == 0)
	return NULL;
    const auto found = p->sharedProducers.find(key.getString());
    return found == p->sharedProducers.end() ? NULL : &found->second;
}

static SbBool
lod_request_demand_is_older(const BObolLodRequest &candidate,
	const BObolLodRequest &current)
{
    return candidate.databaseRevision < current.databaseRevision ||
	candidate.sourceRevision < current.sourceRevision ||
	candidate.viewRevision < current.viewRevision ||
	(candidate.viewRevision == current.viewRevision &&
	 candidate.policyRevision < current.policyRevision) ? TRUE : FALSE;
}

static SbBool
lod_shared_producer_record_lease_unlocked(BObolLodServicePrivate *p,
	const SbString &key, const BObolLodRequest &request,
	uint64_t generation, SbBool *replayReady)
{
    BObolSharedProducer *producer = lod_shared_producer_unlocked(p, key);
    if (!producer)
	return FALSE;
    if (!generation)
	generation = producer->producerGeneration;
    if (!generation || lod_generation_cancelled_unlocked(p, generation))
	return FALSE;

    auto found = producer->leases.find(generation);
    if (found != producer->leases.end() &&
	lod_request_demand_is_older(request, found->second.demand))
	return TRUE;

    /* One generation owns one producer payload.  Other occurrences of the
     * same asset are consumers which the next owner-thread pass binds from
     * that resident asset; they must not steal the producer's result slot.
     * The same occurrence, however, may be revisited under a newer camera or
     * policy epoch while the task is still queued, and that demand must
     * replace its obsolete predecessor. */
    if (found != producer->leases.end() &&
	request.occurrenceKey != found->second.demand.occurrenceKey)
	return TRUE;

    if (found == producer->leases.end()) {
	BObolSharedProducerLease lease;
	lease.demand = request;
	found = producer->leases.emplace(generation, std::move(lease)).first;
    } else {
	const SbBool demandChanged = !bobol_lod_request_keys_equal(
	    request, found->second.demand);
	found->second.demand = request;
	/* A queued payload was stamped for the earlier demand.  Preserve it
	 * for ordinary stale-result accounting and add a current replay edge. */
	if (demandChanged && producer->publication)
	    found->second.deliveredPublication = producer->publication - 1;
    }

    if (generation == producer->producerGeneration)
	p->latestActiveRequests[key.getString()] = request;
    if (replayReady && producer->publication &&
	found->second.deliveredPublication < producer->publication)
	*replayReady = TRUE;
    return TRUE;
}

static SbBool
lod_shared_producer_has_leases_unlocked(
	const BObolLodServicePrivate *p, const BObolLodRequest &request)
{
    if (!request.coalesceAssetProducer)
	return FALSE;
    const BObolSharedProducer *producer = lod_shared_producer_unlocked(
	p, lod_request_active_key(request));
    return producer && !producer->leases.empty() ? TRUE : FALSE;
}

static SbBool
lod_producer_cancelled_unlocked(const BObolLodServicePrivate *p,
	uint64_t generation, const BObolLodRequest &request)
{
    if (lod_shared_producer_has_leases_unlocked(p, request))
	return FALSE;
    return lod_generation_cancelled_unlocked(p, generation);
}

static SbBool
lod_producer_cancelled_or_stopping(BObolLodServicePrivate *p,
	uint64_t generation, const BObolLodRequest &request)
{
    std::lock_guard<std::mutex> lock(p->mutex);
    return p->stopping ||
	lod_producer_cancelled_unlocked(p, generation, request) ? TRUE : FALSE;
}

static SbBool
lod_request_key_recorded_unlocked(const BObolLodServicePrivate *p,
	const SbString &key)
{
    if (lod_active_request_key_recorded_unlocked(p, key))
	return TRUE;
    if (!p || key.getLength() == 0)
	return FALSE;
    const BObolSharedProducer *producer =
	lod_shared_producer_unlocked(p, key);
    if (producer && !producer->leases.empty())
	return TRUE;
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
    else {
	p->activeRequestKeyCounts.erase(found);
	p->latestActiveRequests.erase(key.getString());
    }
}

static SbBool
lod_update_active_request_demand_unlocked(BObolLodServicePrivate *p,
	const SbString &key, const BObolLodRequest &request,
	uint64_t generation, SbBool *replayReady = NULL)
{
    if (!p || key.getLength() == 0)
	return FALSE;
    if (lod_shared_producer_record_lease_unlocked(
	    p, key, request, generation, replayReady))
	return TRUE;
    if (p->activeRequestKeyCounts.find(key.getString()) ==
	p->activeRequestKeyCounts.end())
	return FALSE;

    auto found = p->latestActiveRequests.find(key.getString());
    if (found == p->latestActiveRequests.end()) {
	p->latestActiveRequests.emplace(key.getString(), request);
	return TRUE;
    }
    const BObolLodRequest &current = found->second;
    if (lod_request_demand_is_older(request, current))
	return TRUE;
    found->second = request;
    return TRUE;
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
	status == BOBOL_LOD_PROVIDER_STALE ||
	status == BOBOL_LOD_PROVIDER_SUPERSEDED)
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
	if (lod_producer_cancelled_or_stopping(
		p, task.generation, task.request)) {
	    lod_delayed_task_count_add(p, task.generation, -1);
	    return FALSE;
	}

	uint32_t slice = remaining > 10 ? 10 : remaining;
	std::this_thread::sleep_for(std::chrono::milliseconds(slice));
	remaining -= slice;
    }

    lod_delayed_task_count_add(p, task.generation, -1);
    return !lod_producer_cancelled_or_stopping(
	p, task.generation, task.request);
}

static void
lod_normalize_result(BObolLodResult &result,
	const BObolLodRequest &expectedRequest)
{
    if (!result.cacheKey.isValid()) {
	if (lod_request_has_identity(result.request)) {
	    result.cacheKey = bobol_lod_cache_key(result.request);
	} else {
	    result.request = expectedRequest;
	    result.cacheKey = bobol_lod_cache_key(expectedRequest);
	}
    }

    if (!lod_request_has_identity(result.request))
	result.request = expectedRequest;

    if (!bobol_lod_result_matches_request(result, expectedRequest)) {
	result.providerStatus = BOBOL_LOD_PROVIDER_SUPERSEDED;
	result.stale = TRUE;
	result.terminal = TRUE;
	if (result.diagnostic.getLength() == 0)
	    result.diagnostic = "LoD task was superseded by current demand";
    }
    result.canonicalizePayload();
}

static SbBool
lod_task_dependencies_ready(const BObolLodServicePrivate *p,
			    const BObolLodTask &task)
{
    if (lod_producer_cancelled_unlocked(
	    p, task.generation, task.request))
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

/* librt's database import API reports allocation failure with bu_bomb.  That
 * is appropriate for its historical synchronous callers, but a background
 * display worker must never cross that boundary when an explicit address-space
 * cap has already made the conservative import estimate impossible. */
static bool
lod_raw_import_fits_address_space(const BObolLodRequest &request)
{
#if !defined(__linux__)
    (void)request;
    return true;
#else
    struct rlimit limit;
    if (getrlimit(RLIMIT_AS, &limit) != 0 ||
	limit.rlim_cur == RLIM_INFINITY)
	return true;

    const long pageSize = sysconf(_SC_PAGESIZE);
    if (pageSize <= 0)
	return false;
    FILE *statm = fopen("/proc/self/statm", "r");
    if (!statm)
	return false;
    unsigned long long pages = 0;
    const int read = fscanf(statm, "%llu", &pages);
    fclose(statm);
    if (read != 1 || pages >
	static_cast<unsigned long long>(SIZE_MAX / pageSize))
	return false;
    const size_t used = static_cast<size_t>(pages) *
	static_cast<size_t>(pageSize);
    const size_t cap = limit.rlim_cur > SIZE_MAX ? SIZE_MAX :
	static_cast<size_t>(limit.rlim_cur);
    if (used >= cap)
	return false;

    BObolLodTask estimateTask;
    estimateTask.request = request;
    const size_t estimate = lod_task_estimated_working_set_bytes(estimateTask);
    const size_t remaining = cap - used;
    /* Import holds source arrays beside topology and output state.  Retain
     * half of the remaining address space for librt, GUI, and renderer state. */
    return estimate <= remaining / 2;
#endif
}

SbBool
bobol_lod_spatial_source_enabled(const BObolLodRequest &request,
	size_t workingSetLimit)
{
    const char *requested = getenv("BOBOL_LOD_SPATIAL_LEAVES");
    if (requested && requested[0])
	return BU_STR_EQUAL(requested, "0") ? FALSE : TRUE;
    if (workingSetLimit == SIZE_MAX)
	return FALSE;
    BObolLodTask estimateTask;
    estimateTask.request = request;
    return lod_task_estimated_working_set_bytes(estimateTask) >
	workingSetLimit ? TRUE : FALSE;
}

size_t
bobol_lod_spatial_task_working_set_bytes(void)
{
    /* One 64K-face page owns bounded local vertex maps, cumulative page
     * arrays, immutable publication geometry, and at most 64 MiB of live
     * cache spill.  Reserve substantial allocator/hash-table headroom while
     * remaining well below the ordinary 1 GiB service ceiling. */
    static const size_t bytes = 256ULL * 1024ULL * 1024ULL;
    return bytes;
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
    /* Let an oversized task reach the worker only so it can publish the
     * bounded terminal constraint there.  The worker must not invoke its
     * provider or charge an impossible reservation. */
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

static void
lod_pending_dispatch_remove(BObolLodServicePrivate *p, int dispatchClass)
{
    const auto found = p->pendingDispatchCounts.find(dispatchClass);
    if (found == p->pendingDispatchCounts.end())
	return;
    if (found->second > 1)
	--found->second;
    else
	p->pendingDispatchCounts.erase(found);
}

static std::deque<BObolLodWorkItem>::iterator
lod_find_ready_task(BObolLodServicePrivate *p)
{
    /* First-visible previews are latency work: let them pass queued ordinary
     * refinement without changing request identity or pretending they are a
     * coarser quality tier.  Dependencies and working-set admission are still
     * tested before a candidate can run.  Within one dispatch class, prefer
     * the coarsest ready task (lowest quality tier) so the cheap proxy /
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
    if (!p->pending.empty() && !p->pendingDispatchCounts.empty() &&
	p->pending.front().task.dispatchClass ==
	    p->pendingDispatchCounts.rbegin()->first &&
	!p->pendingQualityCounts.empty() &&
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
	    it->task.dispatchClass > best->task.dispatchClass ||
	    (it->task.dispatchClass == best->task.dispatchClass &&
	     it->task.request.qualityTier < best->task.request.qualityTier))
	    best = it;
    }

    return best;
}

/* Provider callbacks intentionally do not receive scheduler-private work-item
 * metadata.  A worker-local timing scope bridges only the enqueue/start times
 * needed when a mesh provider opens its diagnostic progress record. */
class BObolLodTaskTimingScope {
public:
    explicit BObolLodTaskTimingScope(int64_t submittedMicroseconds) :
	priorSubmitted(lod_task_submitted_microseconds),
	priorStarted(lod_task_started_microseconds)
    {
	lod_task_submitted_microseconds = submittedMicroseconds;
	lod_task_started_microseconds = 0;
    }

    void markStarted(void)
    {
	lod_task_started_microseconds = std::max<int64_t>(1, bu_gettime());
    }

    ~BObolLodTaskTimingScope()
    {
	lod_task_submitted_microseconds = priorSubmitted;
	lod_task_started_microseconds = priorStarted;
    }

private:
    int64_t priorSubmitted;
    int64_t priorStarted;
};

/* A service worker has selected this task, but useful provider execution does
 * not begin until the process-wide CPU gate admits it.  Keep that interval
 * explicit: counting it as producer execution hid cross-service contention
 * and made the progress overlay diagnose a stalled mesh build. */
class BObolCpuAdmissionWaitScope {
public:
    BObolCpuAdmissionWaitScope(BObolLodServicePrivate *service,
	const BObolLodWorkItem &work) : p(service), item(work)
    {
	if (!p)
	    return;
	std::lock_guard<std::mutex> lock(p->mutex);
	p->cpuAdmissionWaitingTasks++;
	lod_generation_count_add_unlocked(
	    p->generationCpuAdmissionWaitingTaskCounts,
	    item.task.generation);
	BObolSharedProducer *shared = lod_shared_producer_unlocked(
	    p, lod_request_active_key(item.task.request));
	if (shared && shared->taskId == item.id)
	    shared->cpuAdmissionWaiting = true;
	active = true;
    }

    ~BObolCpuAdmissionWaitScope(void)
    {
	finish(false);
    }

    void admitted(void)
    {
	finish(true);
    }

private:
    void finish(bool admitted)
    {
	if (!active || !p)
	    return;
	std::lock_guard<std::mutex> lock(p->mutex);
	if (p->cpuAdmissionWaitingTasks > 0)
	    p->cpuAdmissionWaitingTasks--;
	lod_generation_count_remove_unlocked(
	    p->generationCpuAdmissionWaitingTaskCounts,
	    item.task.generation);
	BObolSharedProducer *shared = lod_shared_producer_unlocked(
	    p, lod_request_active_key(item.task.request));
	if (shared && shared->taskId == item.id) {
	    shared->cpuAdmissionWaiting = false;
	    if (admitted)
		shared->executionStartedMicroseconds =
		    std::max<int64_t>(1, bu_gettime());
	}
	active = false;
    }

    BObolLodServicePrivate *p = NULL;
    const BObolLodWorkItem &item;
    bool active = false;
};

/* A task which has passed the service-local byte check can still lose the
 * race for the process-wide transient-memory allowance.  Keep that wait
 * distinct from useful provider execution and from pending-queue pressure so
 * a cold view can explain why its first mesh has not started. */
class BObolTransientMemoryAdmissionWaitScope {
public:
    BObolTransientMemoryAdmissionWaitScope(BObolLodServicePrivate *service,
	const BObolLodWorkItem *work, size_t estimatedBytes) :
	p(service), item(work)
    {
	if (!p || !estimatedBytes)
	    return;
	std::lock_guard<std::mutex> lock(p->mutex);
	p->transientMemoryAdmissionWaitingTasks++;
	if (item) {
	    lod_generation_count_add_unlocked(
		p->generationTransientMemoryAdmissionWaitingTaskCounts,
		item->task.generation);
	    BObolSharedProducer *shared = lod_shared_producer_unlocked(
		p, lod_request_active_key(item->task.request));
	    if (shared && shared->taskId == item->id)
		shared->transientMemoryAdmissionWaiting = true;
	}
	active = true;
    }

    ~BObolTransientMemoryAdmissionWaitScope(void)
    {
	finish();
    }

    void admitted(void)
    {
	finish();
    }

private:
    void finish(void)
    {
	if (!active || !p)
	    return;
	std::lock_guard<std::mutex> lock(p->mutex);
	if (p->transientMemoryAdmissionWaitingTasks > 0)
	    p->transientMemoryAdmissionWaitingTasks--;
	if (item) {
	    lod_generation_count_remove_unlocked(
		p->generationTransientMemoryAdmissionWaitingTaskCounts,
		item->task.generation);
	    BObolSharedProducer *shared = lod_shared_producer_unlocked(
		p, lod_request_active_key(item->task.request));
	    if (shared && shared->taskId == item->id)
		shared->transientMemoryAdmissionWaiting = false;
	}
	active = false;
    }

    BObolLodServicePrivate *p = NULL;
    const BObolLodWorkItem *item = NULL;
    bool active = false;
};

class BObolWorkingSetLease {
public:
    ~BObolWorkingSetLease(void)
    {
	if (bytes)
	    bobol_lod_working_set_release(bytes);
    }

    SbBool acquire(size_t estimatedBytes)
    {
	if (!estimatedBytes)
	    return TRUE;
	if (!bobol_lod_working_set_acquire(estimatedBytes))
	    return FALSE;
	bytes = estimatedBytes;
	return TRUE;
    }

private:
    size_t bytes = 0;
};

static BObolSharedProducer *
lod_shared_producer_for_task_unlocked(BObolLodServicePrivate *p,
	uint64_t taskId);

static BObolLodResult
lod_execute_task(BObolLodServicePrivate *p, const BObolLodWorkItem &item,
	BObolParallelBudgetLease &cpuBudget,
	BObolWorkingSetLease &workingSet)
{
    const BObolLodTask &task = item.task;
    BObolLodTaskTimingScope timing(item.submittedMicroseconds);
    if (lod_producer_cancelled_or_stopping(
	    p, task.generation, task.request))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled");

    if (!lod_wait_for_debug_delay(p, task))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled during debug delay");

    if (!task.realize)
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_ERROR,
					 "LoD task has no realization callback");

    /* Every source and LoD path acquires process CPU before transient memory.
     * A reversed pair here can deadlock with source realization, whose outer
     * CPU lease covers coverage/detail workers that acquire transient memory.
     * The cancellable debug delay remains outside both reservations. */
    BObolCpuAdmissionWaitScope cpuWait(p, item);
    cpuBudget.acquireOuter(task.dispatchClass);
    cpuWait.admitted();

    const size_t estimatedBytes = lod_task_estimated_working_set_bytes(task);
    if (item.exceedsServiceWorkingSetLimit)
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_ERROR,
	    "LoD task exceeds service transient working-set limit");
    BObolTransientMemoryAdmissionWaitScope memoryWait(p, &item,
	estimatedBytes);
    const SbBool admitted = workingSet.acquire(estimatedBytes);
    memoryWait.admitted();
    if (!admitted)
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_ERROR,
	    "LoD task exceeds process transient working-set limit");

    timing.markStarted();
    if (lod_producer_cancelled_or_stopping(
	    p, task.generation, task.request))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled");

    /* Queue age is not demand identity.  A bounded 50k-occurrence scan can
     * leave otherwise cheap tasks queued for many seconds while wheel input
     * advances the camera epoch.  Active-key coalescing records the newest
     * demand, so snapshot it immediately before provider execution.  Cold
     * providers then perform their expensive source work once and construct
     * the view-local result for the newest known occurrence demand instead
     * of completing an obsolete request which the owner must discard.
     *
     * Shared producers retain one payload owner per generation.  Sibling
     * occurrences are deliberately replay consumers and therefore never
     * replace that owner's request in the lease table. */
    BObolLodRequest executionRequest = task.request;
    {
	std::lock_guard<std::mutex> lock(p->mutex);
	const SbString activeKey = lod_request_active_key(task.request);
	const BObolSharedProducer *producer =
	    lod_shared_producer_for_task_unlocked(p, item.id);
	if (producer) {
	    const auto lease = producer->leases.find(task.generation);
	    if (lease != producer->leases.end())
		executionRequest = lease->second.demand;
	} else if (item.retargetActiveDemand) {
	    const auto latest = p->latestActiveRequests.find(
		activeKey.getString());
	    if (latest != p->latestActiveRequests.end())
		executionRequest = latest->second;
	}
    }

    BObolLodResult result = (*task.realize)(
	executionRequest, task.realizeData);

    if (lod_producer_cancelled_or_stopping(
	    p, task.generation, task.request))
	return lod_service_status_result(task, BOBOL_LOD_PROVIDER_CANCELLED,
					 "LoD task generation cancelled");

    BObolLodRequest expectedRequest = task.request;
    {
	/* An asset producer may deliberately finish against a newer camera or
	 * policy demand than the one which launched its immutable hierarchy work.
	 * Validate that result against the service-owned latest demand.  If the
	 * camera moved again after the provider's final refresh, the mismatch is
	 * superseded work and the resident hierarchy will satisfy the next pass. */
	std::lock_guard<std::mutex> lock(p->mutex);
	const SbString activeKey = lod_request_active_key(task.request);
	const BObolSharedProducer *producer =
	    lod_shared_producer_unlocked(p, activeKey);
	if (producer) {
	    const auto lease = producer->leases.find(task.generation);
	    if (lease != producer->leases.end())
		expectedRequest = lease->second.demand;
	} else {
	    const auto latest =
		p->latestActiveRequests.find(activeKey.getString());
	    if (latest != p->latestActiveRequests.end())
		expectedRequest = latest->second;
	}
    }
    lod_normalize_result(result, expectedRequest);
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

static size_t
lod_shared_pending_replay_count_unlocked(
    const BObolSharedProducer &producer, uint64_t generation = 0)
{
    size_t count = 0;
    for (const auto &entry : producer.leases) {
	if (generation && entry.first != generation)
	    continue;
	if (entry.second.deliveredPublication < producer.publication)
	    ++count;
    }
    return count;
}

static SbBool
lod_shared_producer_publish_unlocked(BObolSharedProducer &producer,
	SbBool finalPublication, SbBool ownerPayloadQueued)
{
    bobol_identity_advance(producer.publication);
    producer.state = finalPublication ? BObolSharedProducerState::RESULT :
	BObolSharedProducerState::BUILDING;

    if (ownerPayloadQueued) {
	const auto owner = producer.leases.find(
	    producer.producerGeneration);
	if (owner != producer.leases.end())
	    owner->second.deliveredPublication = producer.publication;
    }
    return lod_shared_pending_replay_count_unlocked(producer) ? TRUE : FALSE;
}

enum class BObolIntermediatePublishStatus {
    REJECTED,
    PUBLISHED,
    CAPACITY_BLOCKED
};

static BObolDeferredIntermediateKey
lod_deferred_intermediate_key(uint64_t generation,
    const BObolLodRequest &request)
{
    BObolDeferredIntermediateKey key;
    key.generation = generation;
    key.activeKey = lod_request_active_key(request).getString();
    return key;
}

/* Take ownership without acquiring the service mutex.  Preview callbacks may
 * hold a resident-asset mutex, so waiting for the outer service lock here
 * would invert the documented lock order.  Success means the latest preview
 * is durably retained for promotion, not necessarily that it is already in
 * the presentation result queue. */
static bool
lod_defer_intermediate_result(BObolLodServicePrivate *p,
    uint64_t generation, BObolLodResult &&result)
{
    if (!p || !generation)
	return false;
    BObolDeferredIntermediateKey key = lod_deferred_intermediate_key(
	generation, result.request);
    if (key.activeKey.empty())
	return false;
    result.generation = generation;

    BObolLodResult retired;
    bool accepted = false;
    try {
	std::lock_guard<std::mutex> lock(p->deferredIntermediateMutex);
	const auto found = p->deferredIntermediateResults.find(key);
	if (found != p->deferredIntermediateResults.end()) {
	    if (lod_result_supersedes(result, found->second)) {
		retired = std::move(found->second);
		found->second = std::move(result);
	    }
	    accepted = true;
	} else if (p->deferredIntermediateResults.size() <
		p->maxDeferredIntermediateResults.load(
		    std::memory_order_relaxed)) {
	    p->deferredIntermediateResults.emplace(
		std::move(key), std::move(result));
	    p->deferredIntermediateResultCount.store(
		p->deferredIntermediateResults.size(),
		std::memory_order_release);
	    accepted = true;
	}
    } catch (const std::bad_alloc &) {
	return false;
    }
    return accepted;
}

/* A capacity-blocked value came from an earlier mailbox snapshot.  If its
 * producer published again while that snapshot was being validated, the
 * value now in the mailbox is necessarily newer and must win even when its
 * demand revisions compare equal. */
static void
lod_restore_capacity_blocked_intermediate(BObolLodServicePrivate *p,
    BObolLodResult &&result)
{
    if (!p || !result.generation)
	return;
    BObolDeferredIntermediateKey key = lod_deferred_intermediate_key(
	result.generation, result.request);
    if (key.activeKey.empty())
	return;
    try {
	std::lock_guard<std::mutex> lock(p->deferredIntermediateMutex);
	if (p->deferredIntermediateResults.find(key) !=
		p->deferredIntermediateResults.end() ||
	    p->deferredIntermediateResults.size() >=
		p->maxDeferredIntermediateResults.load(
		    std::memory_order_relaxed))
	    return;
	p->deferredIntermediateResults.emplace(
	    std::move(key), std::move(result));
	p->deferredIntermediateResultCount.store(
	    p->deferredIntermediateResults.size(),
	    std::memory_order_release);
    } catch (const std::bad_alloc &) {
	/* The ordinary final task result remains authoritative. */
    }
}

/* Validate and publish one mailbox value while the service mutex is held.
 * CAPACITY_BLOCKED leaves result intact so the caller can return it to the
 * bounded mailbox after dropping the service lock. */
static BObolIntermediatePublishStatus
lod_publish_intermediate_result_unlocked(BObolLodServicePrivate *p,
    uint64_t generation, BObolLodResult &result, SbBool &notifyResultReady)
{
    const SbString activeKey = lod_request_active_key(result.request);
    BObolSharedProducer *shared = lod_shared_producer_unlocked(p, activeKey);
    if (p->stopping ||
	lod_producer_cancelled_unlocked(p, generation, result.request) ||
	lod_generation_count_unlocked(
	    p->generationExecutingTaskCounts, generation) == 0 ||
	(!lod_active_request_key_recorded_unlocked(p, activeKey) &&
	 !(shared && !shared->leases.empty())))
	return BObolIntermediatePublishStatus::REJECTED;

    if (shared) {
	const auto owner = shared->leases.find(generation);
	if (owner != shared->leases.end())
	    lod_normalize_result(result, owner->second.demand);
	else
	    lod_normalize_result(result, result.request);
    } else {
	lod_normalize_result(result, result.request);
    }
    result.generation = generation;
    const SbBool queueOwnerPayload = !shared ||
	shared->leases.find(generation) != shared->leases.end();

    if (queueOwnerPayload) {
	const BObolLodResultSlotMapKey slot =
	    lod_result_slot_map_key(result);
	const auto existing = p->resultSlots.find(slot);
	if (existing != p->resultSlots.end()) {
	    if (lod_result_supersedes(result, *existing->second)) {
		const SbString oldRequestKey = lod_request_active_key(
		    existing->second->request);
		const SbString newRequestKey = lod_request_active_key(
		    result.request);
		if (oldRequestKey != newRequestKey) {
		    lod_queued_result_request_key_remove_unlocked(
			p, existing->second->request);
		    lod_queued_result_request_key_add_unlocked(
			p, result.request);
		}
		*existing->second = std::move(result);
	    }
	    p->coalescedResults++;
	} else {
	    if (p->results.size() >= p->maxQueuedResults)
		return BObolIntermediatePublishStatus::CAPACITY_BLOCKED;
	    notifyResultReady = lod_generation_count_unlocked(
		p->generationResultCounts, generation) == 0 ? TRUE :
		    notifyResultReady;
	    lod_queued_result_request_key_add_unlocked(p, result.request);
	    p->results.push_back(std::move(result));
	    p->resultSlots.emplace(slot, std::prev(p->results.end()));
	    lod_generation_count_add_unlocked(
		p->generationResultCounts, generation);
	}
    }
    if (shared && lod_shared_producer_publish_unlocked(
	    *shared, FALSE, queueOwnerPayload))
	notifyResultReady = TRUE;
    return BObolIntermediatePublishStatus::PUBLISHED;
}

/* Promote all currently deferred values.  A provider uses the try-lock form;
 * worker and presentation threads use the blocking form only when they hold
 * no resident asset.  Subscriber notification is handed to a safe service
 * worker rather than invoked from a possibly resident-locked callback. */
static bool
lod_promote_deferred_intermediate_results(BObolLodServicePrivate *p,
    bool waitForServiceLock, bool dispatchNotification)
{
    if (!p)
	return false;

    const auto dispatchPending = [p, dispatchNotification]() {
	if (!dispatchNotification ||
	    !p->deferredIntermediateNotificationPending.exchange(
		false, std::memory_order_acq_rel))
	    return;
	lod_notify_result_ready(p);
    };
    if (!p->deferredIntermediateResultCount.load(std::memory_order_acquire)) {
	dispatchPending();
	return true;
    }

    std::unique_lock<std::mutex> serviceLock(p->mutex, std::defer_lock);
    if (waitForServiceLock) {
	serviceLock.lock();
    } else if (!serviceLock.try_lock()) {
	return false;
    }

    std::unordered_map<BObolDeferredIntermediateKey, BObolLodResult,
	BObolDeferredIntermediateKeyHash> pending;
    {
	std::lock_guard<std::mutex> mailboxLock(
	    p->deferredIntermediateMutex);
	pending.swap(p->deferredIntermediateResults);
	p->deferredIntermediateResultCount.store(0,
	    std::memory_order_release);
    }

    SbBool notifyResultReady = FALSE;
    std::vector<BObolLodResult> capacityBlocked;
    capacityBlocked.reserve(pending.size());
    for (auto &entry : pending) {
	BObolIntermediatePublishStatus status =
	    lod_publish_intermediate_result_unlocked(
		p, entry.first.generation, entry.second,
		notifyResultReady);
	if (status == BObolIntermediatePublishStatus::CAPACITY_BLOCKED)
	    capacityBlocked.push_back(std::move(entry.second));
    }
    serviceLock.unlock();

    for (BObolLodResult &blocked : capacityBlocked)
	lod_restore_capacity_blocked_intermediate(p, std::move(blocked));
    if (notifyResultReady)
	p->deferredIntermediateNotificationPending.store(
	    true, std::memory_order_release);
    dispatchPending();
    return true;
}

static BObolLodResult
lod_shared_producer_replay_result(uint64_t generation,
	const BObolLodRequest &request)
{
    BObolLodResult result;
    result.generation = generation;
    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_DIAGNOSTIC;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_SUPERSEDED;
    result.terminal = TRUE;
    result.stale = TRUE;
    result.diagnostic =
	"shared asset producer completed; replay current consumer demand";
    result.canonicalizePayload();
    return result;
}

static std::unordered_map<std::string, BObolSharedProducer>::iterator
lod_shared_producer_erase_unlocked(
    BObolLodServicePrivate *p,
    std::unordered_map<std::string, BObolSharedProducer>::iterator producer)
{
    if (!p || producer == p->sharedProducers.end())
	return producer;
    const uint64_t taskId = producer->second.taskId;
    if (lod_generation_cancelled_unlocked(
	    p, producer->second.producerGeneration)) {
	p->taskGenerations.erase(taskId);
	p->completed.erase(taskId);
    }
    p->sharedProducerTaskKeys.erase(taskId);
    return p->sharedProducers.erase(producer);
}

static void
lod_shared_producer_retire_if_unowned_unlocked(
    BObolLodServicePrivate *p,
    std::unordered_map<std::string, BObolSharedProducer>::iterator producer)
{
    if (!p || producer == p->sharedProducers.end() ||
	!producer->second.leases.empty())
	return;
    (void)lod_shared_producer_erase_unlocked(p, producer);
}

static BObolSharedProducer *
lod_shared_producer_for_task_unlocked(BObolLodServicePrivate *p,
	uint64_t taskId)
{
    if (!p || !taskId)
	return NULL;
    const auto key = p->sharedProducerTaskKeys.find(taskId);
    if (key == p->sharedProducerTaskKeys.end())
	return NULL;
    const auto producer = p->sharedProducers.find(key->second);
    return producer == p->sharedProducers.end() ? NULL : &producer->second;
}

static void
lod_finish_task(BObolLodServicePrivate *p, const BObolLodWorkItem &item,
		BObolLodResult &&result)
{
    /* The provider has returned and therefore owns no resident-asset lock.
     * Give any callback-deferred preview one final promotion opportunity
     * before its authoritative completion is coalesced into the same queue. */
    (void)lod_promote_deferred_intermediate_results(p, true, false);
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

	const SbString activeKey = lod_request_active_key(item.task.request);
	auto sharedEntry = p->sharedProducers.find(activeKey.getString());
	const SbBool sharedTask = sharedEntry != p->sharedProducers.end() &&
	    sharedEntry->second.taskId == item.id ? TRUE : FALSE;
	BObolSharedProducer *shared = sharedTask ?
	    &sharedEntry->second : NULL;
	const SbBool discardResult = p->stopping ||
	    (shared ? shared->leases.empty() :
	     lod_generation_cancelled_unlocked(p, item.task.generation));
	const SbBool queueOwnerPayload = !shared ||
	    shared->leases.find(shared->producerGeneration) !=
		shared->leases.end();
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
	p->activeWorkingSetBytes =
	    item.reservedWorkingSetBytes >= p->activeWorkingSetBytes ?
	    0 : p->activeWorkingSetBytes - item.reservedWorkingSetBytes;
	lod_producer_progress_prune(p);
	lod_active_request_key_remove_unlocked(p, activeKey);

	if (item.task.publishResult) {
	    if (p->resultReservations > 0)
		p->resultReservations--;
	    if (discardResult) {
		p->discardedStaleResults++;
	    } else {
		if (shared && lod_shared_producer_publish_unlocked(
			*shared, TRUE, queueOwnerPayload))
		    notifyResultReady = TRUE;
		if (queueOwnerPayload) {
		    /* The final worker publication supersedes every undrained
		     * preview for this producer/consumer.  Keeping an older slot of
		     * another result kind would make lease retirement depend on drain
		     * order rather than the final publication witness. */
		    if (shared) {
			for (auto queued = p->results.begin();
			     queued != p->results.end();) {
			    if (queued->generation !=
				    shared->producerGeneration ||
				lod_request_active_key(queued->request) !=
				    activeKey) {
				++queued;
				continue;
			    }
			    lod_queued_result_request_key_remove_unlocked(
				p, queued->request);
			    p->resultSlots.erase(
				lod_result_slot_map_key(*queued));
			    lod_generation_count_remove_unlocked(
				p->generationResultCounts,
				queued->generation);
			    queued = p->results.erase(queued);
			}
		    }
		    const BObolLodResultSlotMapKey slot =
			lod_result_slot_map_key(completedResult);
		    const auto existing = p->resultSlots.find(slot);
		    if (existing != p->resultSlots.end()) {
			if (lod_result_supersedes(
				completedResult, *existing->second)) {
			    const SbString oldRequestKey =
				lod_request_active_key(
				    existing->second->request);
			    const SbString newRequestKey =
				lod_request_active_key(
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
			if (lod_generation_count_unlocked(
				p->generationResultCounts,
				completedResult.generation) == 0)
			    notifyResultReady = TRUE;
			lod_queued_result_request_key_add_unlocked(
			    p, completedResult.request);
			p->results.push_back(std::move(completedResult));
			p->resultSlots.emplace(
			    slot, std::prev(p->results.end()));
			lod_generation_count_add_unlocked(
			    p->generationResultCounts,
			    p->results.back().generation);
		    }
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
		if (shared && !shared->leases.empty())
		    writeItem.result.generation = shared->leases.begin()->first;
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
	if (shared && discardResult)
	    lod_shared_producer_retire_if_unowned_unlocked(p, sharedEntry);
    }

    p->workerCv.notify_all();
    p->cacheWriterCv.notify_one();
    if (notifyResultReady ||
	p->deferredIntermediateNotificationPending.exchange(
	    false, std::memory_order_acq_rel))
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
	if (!lod_resident_consumer_snapshot_current(consumer.second))
	    continue;
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
    if (target == p->residentMeshCompactionTargets.end() ||
	target->second.demandEpoch != p->residentMeshDemandEpoch) {
	p->residentMeshCompactionQueuedAssets.erase(candidate.assetKey);
	if (target != p->residentMeshCompactionTargets.end())
	    p->residentMeshCompactionTargets.erase(target);
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
	    found->second.demandEpoch == p->residentMeshDemandEpoch &&
	    found->second.revision == work.target.revision &&
	    found->second.useRevision == work.target.useRevision &&
	    found->second.cut == work.target.cut &&
	    found->second.channelMask == work.target.channelMask &&
	    found->second.chunkCuts == work.target.chunkCuts &&
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
	for (const auto &consumer : p->residentMeshConsumerDemands) {
	    if (!lod_resident_consumer_snapshot_current(consumer.second))
		continue;
	    if (consumer.second.assets.find(work.assetKey) !=
		    consumer.second.assets.end())
		return result;
	}
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
	lod_resident_mesh_accounting_replace(
	    p, priorBytes, 0, priorBackingBytes, 0);
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
    SbBool preparedWorkingSetChanged = FALSE;
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
	if (!work.target.chunkCuts.empty()) {
	    std::vector<BObolLodChunkCut> residentChunkCuts;
	    mesh->residentChunkCuts(residentChunkCuts);
	    preparedWorkingSetChanged =
		residentChunkCuts != work.target.chunkCuts ? TRUE : FALSE;
	    for (const BObolLodChunkCut &chunk : work.target.chunkCuts)
		preparedTargetCut = std::max(preparedTargetCut, chunk.cut);
	} else {
	    preparedTargetCut = std::min(currentCut,
		std::max(minimumCut, std::min(maximumCut,
		    work.target.cut < 0 ? minimumCut : work.target.cut)));
	    preparedWorkingSetChanged =
		preparedTargetCut < currentCut ? TRUE : FALSE;
	}
	if (preparedWorkingSetChanged) {
	    if (!bobol_mesh_lod_hierarchy_info_get(
		    resident->lod, &hierarchy))
		return result;
	}
    }
    if (preparedWorkingSetChanged) {
	preparedTrim = work.target.chunkCuts.empty() ?
	    mesh->prepareTrim(preparedTargetCut) :
	    mesh->prepareTrim(work.target.chunkCuts);
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
    int targetCut = -1;
    if (!work.target.chunkCuts.empty()) {
	for (const BObolLodChunkCut &chunk : work.target.chunkCuts)
	    targetCut = std::max(targetCut, chunk.cut);
    } else {
	targetCut = std::min(currentCut,
	    std::max(minimumCut, std::min(maximumCut,
		work.target.cut < 0 ? minimumCut : work.target.cut)));
    }
    if (targetCut != preparedTargetCut)
	return result;
    const size_t priorBytes =
	resident->publishedBytes.load(std::memory_order_relaxed);
    const size_t priorBackingBytes =
	resident->publishedBackingPrefixBytes.load(
	    std::memory_order_relaxed);

    if (targetCut >= currentCut && !preparedTrim) {
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
	if (result.residentBytes != priorBytes)
	    lod_resident_mesh_revision_advance(p->residentMeshRevision);
	lod_resident_mesh_accounting_replace(
	    p, priorBytes, result.residentBytes, priorBackingBytes, 0);
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
    lod_resident_mesh_accounting_replace(
	p, priorBytes, result.residentBytes, priorBackingBytes, 0);
    lod_resident_mesh_revision_advance(p->residentMeshRevision);
    const size_t priorStableBytes =
	lod_resident_asset_stable_bytes(priorBytes, priorBackingBytes);
    if (result.residentBytes < priorStableBytes)
	lod_resident_mesh_revision_advance(p->residentMeshAdmissionRevision);
    residentLock.unlock();

    if (result.channelMask) {
	const int drawMode = result.channelMask == 3u ?
	    BOBOL_LOD_DRAW_HIDDEN_LINE :
	    (result.channelMask & 1u ?
		BOBOL_LOD_DRAW_WIRE : BOBOL_LOD_DRAW_SHADED);
	if (work.target.chunkCuts.empty()) {
	    result.preparedCadGeometry = mesh->prepareCadGeometry(
		drawMode, &result.preparedCadGeometryRevision);
	}
	/*
	 * Spatial presentation layers are immutable view products, not resident
	 * hierarchy storage.  Rebuilding them here would require choosing one
	 * normal policy for every view sharing this asset, and historically let a
	 * quiet compaction replace smooth or flat shading with authored shading.
	 * The owner keeps any still-drawable layer set across the trim.  Its normal
	 * request prepares a replacement on a worker only when the cut, visible
	 * pages, or presentation policy actually changes.
	 */
    }
    if (getenv("BOBOL_LOD_TRACE_COMPACTION")) {
	std::vector<BObolLodChunkCut> retained;
	mesh->residentChunkCuts(retained);
	SbBool demandedWorkingSetDrawable = TRUE;
	for (const BObolLodChunkCut &chunk : work.target.chunkCuts) {
	    const std::vector<uint32_t> oneChunk = {chunk.chunkId};
	    if (!mesh->canDrawChunksAtCut(oneChunk, chunk.cut))
		demandedWorkingSetDrawable = FALSE;
	}
	size_t preparedFaces = 0;
	if (result.preparedCadGeometry &&
	    result.preparedCadGeometry->shaded) {
	    for (const Obol::ProgressiveTriangleCluster &cluster :
		    result.preparedCadGeometry->shaded->progressiveClusters)
		for (const Obol::ProgressiveTriangleClusterRange &range :
			cluster.ranges)
		    if (range.activationCut <= preparedTargetCut)
			preparedFaces += range.indexCount / 3;
	}
	for (const BObolLodPresentationLayer &layer :
		result.presentationLayers)
	    if (layer.geometry && layer.geometry->shaded)
		preparedFaces += layer.geometry->shaded->indexCountAtCut(
		    static_cast<uint8_t>(std::max(0, layer.activeCut))) / 3;
	bu_log("BObol resident compaction publish asset=%s cut=%d "
	       "target_chunks=%zu retained_chunks=%zu drawable=%d "
	       "revision=%llu prepared_revision=%llu prepared_faces=%zu\n",
	       work.assetKey.c_str(), targetCut, work.target.chunkCuts.size(),
	       retained.size(),
	       work.target.chunkCuts.empty() ? 0 :
		   (demandedWorkingSetDrawable ? 1 : 0),
	       static_cast<unsigned long long>(mesh->revision()),
	       static_cast<unsigned long long>(
		   result.preparedCadGeometryRevision),
	       preparedFaces);
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
		if (!lod_resident_consumer_snapshot_current(
			consumer->second))
		    continue;
		const auto demand =
		    consumer->second.assets.find(work.assetKey);
		if (demand == consumer->second.assets.end())
		    continue;
		const SbBool demandDrawable =
		    !demand->second.chunkIds.empty() && result.progressiveMesh ?
		    result.progressiveMesh->canDrawChunksAtCut(
			demand->second.chunkIds, demand->second.cut) :
		    (demand->second.cut <= result.residentCut ? TRUE : FALSE);
		if (!demandDrawable)
		    continue;
		BObolLodResidentCompaction consumerResult = result;
		consumerResult.consumerDemandRevision =
		    consumer->second.revision;
		p->residentMeshCompactionResults[consumerId].push_back(
		    std::move(consumerResult));
		p->residentMeshCompactionResultCount++;
	    }
	    notifyResultReady =
		resultCountBefore == 0 &&
		p->residentMeshCompactionResultCount > 0 ? TRUE : FALSE;
	}
    }
    p->workerCv.notify_all();
    /* An obsolete asset has no current consumer and therefore produces no
     * presentation result, but freeing its stable bytes advances the
     * admission epoch.  Notify subscribers so memory-limited current assets
     * can consume that capacity edge; otherwise a quiet view can terminate
     * with ample reclaimed headroom and no event capable of starting its
     * sparse retry pass. */
    if (reclaimedStorage)
	notifyResultReady = TRUE;
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
	/* A producer which could not acquire the service lock left its newest
	 * immutable preview in the side mailbox.  An idle peer promotes and
	 * announces it without extending the inverted resident-lock interval. */
	(void)lod_promote_deferred_intermediate_results(p, true, true);
	BObolLodWorkItem item;
	BObolResidentMeshCompactionWork compaction;
	SbBool runCompaction = FALSE;
	SbBool dispatchDeferred = FALSE;

	{
	    std::unique_lock<std::mutex> lock(p->mutex);

	    std::deque<BObolLodWorkItem>::iterator ready = p->pending.end();
	    for (;;) {
		if (p->stopping)
		    return;
		/* A full presentation queue leaves a capacity-blocked preview in
		 * the mailbox.  Do not spin retrying it until a drain wakes us;
		 * pending subscriber notification does not consume queue space. */
		const bool deferredCapacityAvailable =
		    p->results.size() < p->maxQueuedResults;
		if ((deferredCapacityAvailable &&
		     p->deferredIntermediateResultCount.load(
			 std::memory_order_acquire)) ||
		    p->deferredIntermediateNotificationPending.load(
			std::memory_order_acquire)) {
		    dispatchDeferred = TRUE;
		    break;
		}
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

	    if (dispatchDeferred) {
		/* No queue item was selected.  Promote outside the service lock. */
	    } else if (runCompaction) {
		/* Counters and reservations were installed by the take helper. */
	    } else {
		const int qualityTier = ready->task.request.qualityTier;
		const int dispatchClass = ready->task.dispatchClass;
		item = std::move(*ready);
		p->pending.erase(ready);
		lod_generation_count_remove_unlocked(
		    p->generationPendingTaskCounts,
		    item.task.generation);
		lod_generation_pending_time_remove_unlocked(
		    p, item.task.generation,
		    item.submittedMicroseconds);
		lod_generation_count_add_unlocked(
		    p->generationExecutingTaskCounts,
		    item.task.generation);
		lod_pending_quality_remove(p, qualityTier);
		lod_pending_dispatch_remove(p, dispatchClass);
		const size_t workingBytes =
		    lod_task_estimated_working_set_bytes(item.task);
		/* Let a known-over-limit task leave the queue and publish an explicit
		 * terminal constraint below.  Charging it to the active tally would
		 * falsely report an impossible reservation and can block unrelated
		 * affordable work while no provider has started. */
		const bool exceedsServiceLimit =
		    p->maxActiveWorkingSetBytes != SIZE_MAX &&
		    workingBytes > p->maxActiveWorkingSetBytes;
		const size_t accountedBytes = exceedsServiceLimit ? 0 : workingBytes;
		item.exceedsServiceWorkingSetLimit = exceedsServiceLimit;
		item.reservedWorkingSetBytes = accountedBytes;
		p->activeWorkingSetBytes =
		    accountedBytes > SIZE_MAX - p->activeWorkingSetBytes ?
		    SIZE_MAX : p->activeWorkingSetBytes + accountedBytes;
		p->executingTasks++;
		p->peakWorkingSetBytes = std::max(
		    p->peakWorkingSetBytes, p->activeWorkingSetBytes);
		p->peakExecutingTasks = std::max(
		    p->peakExecutingTasks, p->executingTasks);
	    }
	}
	if (dispatchDeferred)
	    continue;

	if (runCompaction) {
	    /* Match source and ordinary LoD execution: never retain transient
	     * memory while waiting for the process CPU gate. */
	    BObolParallelBudgetLease cpuBudget;
	    cpuBudget.acquireOuter();
	    BObolTransientMemoryAdmissionWaitScope memoryWait(
		p, NULL, compaction.estimatedWorkingSetBytes);
	    BObolWorkingSetLease workingSet;
	    const SbBool admitted = workingSet.acquire(
		compaction.estimatedWorkingSetBytes);
	    memoryWait.admitted();
	    BObolLodResidentCompaction result;
	    if (admitted)
		result = lod_execute_resident_compaction(p, compaction);
	    lod_finish_resident_compaction(
		p, compaction, std::move(result));
	    continue;
	}

	/* Keep the transient reservation through result publication and provider
	 * payload cleanup, matching its estimate's complete task lifetime.  CPU
	 * admission only covers provider execution; queue bookkeeping is short and
	 * must not occupy a scarce background execution slot. */
	BObolParallelBudgetLease cpuBudget;
	BObolWorkingSetLease workingSet;
	BObolLodResult result = lod_execute_task(
	    p, item, cpuBudget, workingSet);
	cpuBudget.release();
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
	if (item.write && !lod_producer_cancelled_or_stopping(
		p, item.result.generation, item.result.request)) {
	    BObolParallelBudgetLease cpuBudget;
	    cpuBudget.acquireOuter();
	    if (!lod_producer_cancelled_or_stopping(
		    p, item.result.generation, item.result.request))
		(*item.write)(item.result, item.writeData);
	}

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
	this->p->maxDeferredIntermediateResults.store(
	    std::max<size_t>(1, std::min(
		this->p->maxQueuedResults, workerCount)),
	    std::memory_order_relaxed);
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
	this->p->maxDeferredIntermediateResults.store(
	    std::max<size_t>(1, std::min(
		this->p->maxQueuedResults, this->p->workers.size())),
	    std::memory_order_relaxed);
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
	this->p->deferredIntermediateNotificationPending.store(
	    false, std::memory_order_release);
    }
    this->p->workerCv.notify_all();

    workers.swap(this->p->workers);
    for (size_t i = 0; i < workers.size(); i++) {
	if (workers[i].joinable())
	    workers[i].join();
    }

    {
	std::lock_guard<std::mutex> lock(
	    this->p->deferredIntermediateMutex);
	this->p->deferredIntermediateResults.clear();
	this->p->deferredIntermediateResultCount.store(
	    0, std::memory_order_release);
	this->p->deferredIntermediateNotificationPending.store(
	    false, std::memory_order_release);
    }

    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	pending.swap(this->p->pending);
	this->p->pendingDispatchCounts.clear();
	this->p->pendingQualityCounts.clear();
	this->p->inFlight = 0;
	this->p->activeWorkingSetBytes = 0;
	this->p->executingTasks = 0;
	this->p->cpuAdmissionWaitingTasks = 0;
	this->p->transientMemoryAdmissionWaitingTasks = 0;
	this->p->resultReservations = 0;
	this->p->cacheWriteReservations = 0;
	this->p->activeRequestKeyCounts.clear();
	this->p->latestActiveRequests.clear();
	this->p->sharedProducers.clear();
	this->p->sharedProducerTaskKeys.clear();
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
	this->p->generationPendingTaskTimes.clear();
	this->p->generationExecutingTaskCounts.clear();
	this->p->generationCpuAdmissionWaitingTaskCounts.clear();
	this->p->generationTransientMemoryAdmissionWaitingTaskCounts.clear();
	this->p->generationDelayedTaskCounts.clear();
	this->p->generationResultCounts.clear();
	this->p->generationCacheWriteCounts.clear();
	{
	    std::lock_guard<std::mutex> progressLock(
		this->p->producerProgressMutex);
	    this->p->producerProgress.clear();
	    this->p->nextProducerProgressId = 1;
	}
	this->p->residentMeshConsumerDemands.clear();
	this->p->residentMeshOrder.clear();
	this->p->residentMeshCompactionWork.clear();
	this->p->residentMeshCompactionQueuedAssets.clear();
	this->p->residentMeshCompactionTargets.clear();
	this->p->residentMeshDemandEpoch = 1;
	this->p->residentMeshCompactionResults.clear();
	this->p->residentMeshCompactionsInFlight = 0;
	this->p->residentMeshCompactionResultCount = 0;
	this->p->residentMeshCompactionResultReservations = 0;
	residentMeshes.swap(this->p->residentMeshes);
	this->p->residentMeshBytes.store(0, std::memory_order_relaxed);
	this->p->residentMeshBackingBytes.store(
	    0, std::memory_order_relaxed);
	this->p->residentMeshStableBytes.store(
	    0, std::memory_order_relaxed);
	{
	    std::lock_guard<std::mutex> admissionLock(
		this->p->residentMeshAdmissionMutex);
	    this->p->residentMeshGrowthReservationBytes = 0;
	}
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
lod_refinement_growth_budget(size_t current, size_t initial,
			     double growth)
{
    if (!current)
	return initial;
    if (!std::isfinite(growth) || growth <= 1.0)
	growth = 1.0;
    const long double grown =
	static_cast<long double>(current) * static_cast<long double>(growth);
    const size_t growthBudget =
	grown >= static_cast<long double>(
	    std::numeric_limits<size_t>::max()) ?
	std::numeric_limits<size_t>::max() :
	static_cast<size_t>(std::ceil(grown));
    return std::max(initial, growthBudget);
}

static BObolLodCounts
lod_progressive_delivery_counts(
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    const std::vector<uint32_t> *requiredChunks, int cut)
{
    BObolLodCounts counts;
    if (cut < hierarchy.min_cut || cut > hierarchy.max_cut)
	return counts;

    if (!hierarchy.chunks || !requiredChunks || requiredChunks->empty()) {
	counts.faceCount = hierarchy.cuts[cut].face_count;
	counts.pointCount = hierarchy.cuts[cut].point_count;
    } else {
	for (uint32_t chunk : *requiredChunks) {
	    if (chunk >= hierarchy.chunk_count) {
		counts.clear();
		return counts;
	    }
	    const BObolMeshLodChunkCutInfo &population =
		hierarchy.chunks[chunk].cuts[cut];
	    counts.faceCount = population.face_count >
		    UINT64_MAX - counts.faceCount ?
		UINT64_MAX : counts.faceCount + population.face_count;
	    counts.pointCount = population.point_count >
		    UINT64_MAX - counts.pointCount ?
		UINT64_MAX : counts.pointCount + population.point_count;
	}
    }
    counts.originalPointCount = counts.pointCount;
    /* Preserve the existing delivery-price convention.  Renderer admission
     * accounts for any later corner-normal expansion from exact prepared
     * geometry; this estimate only bounds one progressive worker step. */
    counts.normalCount = hierarchy.has_normals ? counts.pointCount : 0;
    return counts;
}

/* Select one bounded presentation step toward the producer-resolved
 * screen-error target.  This helper only prevents a first box->mesh
 * replacement from allocating and uploading an arbitrarily large cumulative
 * prefix; it does not mutate asynchronous request identity. */
static int
lod_progressive_delivery_cut(
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int requestedCut, int currentCut, int drawMode,
    const BObolMeshLodProvider &provider,
    const std::vector<uint32_t> *requiredChunks)
{
    int target = std::max(hierarchy.min_cut,
	std::min(hierarchy.max_cut, requestedCut));
    if (provider.useDeliveryCutLimit &&
	provider.deliveryCutLimit >= hierarchy.min_cut)
	target = std::min(target, provider.deliveryCutLimit);
    if (!provider.progressiveDelivery || provider.useForcedCut ||
	provider.resetExisting || currentCut >= target)
	return target;

    const BObolLodCounts currentCounts = lod_progressive_delivery_counts(
	hierarchy, requiredChunks, currentCut);
    const size_t currentCost = bobol_lod_render_cost_units(
	currentCounts, drawMode, 1);
    const size_t costBudget = lod_refinement_growth_budget(currentCost,
	provider.initialRefinementCostBudget,
	provider.refinementGrowthFactor);

    int selected = hierarchy.min_cut;
    bool selectedPopulatedCut = false;
    for (int cut = hierarchy.min_cut; cut <= target; ++cut) {
	const BObolLodCounts counts = lod_progressive_delivery_counts(
	    hierarchy, requiredChunks, cut);
	const bool populated = counts.faceCount && counts.pointCount;
	if (!populated) {
	    if (!selectedPopulatedCut)
		selected = cut;
	    continue;
	}
	if (bobol_lod_render_cost_units(counts, drawMode, 1) >
	    costBudget) {
	    /* A page set's first drawable population is indivisible.  The scene
	     * reservation estimates rather than opens every cold hierarchy on the
	     * owner thread; when that estimate is low, publish the smallest useful
	     * prefix and let the next measured frame correct the allowance.
	     * Returning the preceding empty cut turns a valid hierarchy into a
	     * terminal provider failure. */
	    if (!selectedPopulatedCut)
		selected = cut;
	    break;
	}
	selected = cut;
	selectedPopulatedCut = true;
    }

    /*
     * Return the richest population that fits the provider's render-cost
     * growth allowance.  Nominal cut adjacency is deliberately irrelevant:
     * one cut may add no arrays while the next adds millions.  The budget is
     * the bounded work contract, not "+1 cut".
     */
    if (currentCut >= hierarchy.min_cut && selected <= currentCut)
	selected = std::min(target, currentCut + 1);
    return std::min(selected, target);
}

struct BObolColdMeshPreviewContext {
    BObolLodService *service = NULL;
    BObolLodServicePrivate *serviceState = NULL;
    uint64_t generation = 0;
    BObolLodRequest request;
    BObolLodCacheKey assetKey;
    BObolLodProgressiveMeshPtr progressiveMesh;
    /* The submit action partitions one scene-wide first-pass allowance
     * between tasks.  A cold cache callback runs before the ordinary result
     * path, so it must carry that same per-task share explicitly rather than
     * treating a hierarchy minimum as a free preview. */
    size_t renderCostAllowance = 0;
    std::vector<BObolLodPresentationLayer> presentationLayers;
    size_t presentationRenderCost = 0;
    size_t publishedPageCount = 0;
    BObolViewEpoch publishedViewRevision;
    BObolPolicyEpoch publishedPolicyRevision;
    SbString publishedOccurrenceKey;
    uint32_t publishedSourceEntryIndex = UINT32_MAX;
    SbBool spatialLeafProducer = FALSE;
    SbBool spatialCoverageAdmitted = FALSE;
    std::shared_ptr<BObolLodProducerProgressRecord> progress;
};

static void
lod_cold_preview_refresh_demand(BObolColdMeshPreviewContext *context)
{
    if (!context || !context->serviceState)
	return;
    const SbString key = lod_request_active_key(context->request);
    /* Cache callbacks run while the producer owns resident.mutex.  A demand
     * refresh is opportunistic, so never wait for the outer service lock in
     * the inverse direction.  The next page callback (or final result
     * normalization) observes any update missed by this probe. */
    std::unique_lock<std::mutex> lock(
	context->serviceState->mutex, std::try_to_lock);
    if (!lock.owns_lock())
	return;
    const BObolSharedProducer *producer = lod_shared_producer_unlocked(
	context->serviceState, key);
    if (producer) {
	const auto lease = producer->leases.find(context->generation);
	if (lease != producer->leases.end())
	    context->request = lease->second.demand;
	return;
    }
    const auto found = context->serviceState->latestActiveRequests.find(
	key.getString());
    if (found != context->serviceState->latestActiveRequests.end())
	context->request = found->second;
}

static SbBool
lod_cold_preview_demand_unpublished(
    const BObolColdMeshPreviewContext *context)
{
    if (!context)
	return FALSE;
    return context->publishedViewRevision != context->request.viewRevision ||
	context->publishedPolicyRevision != context->request.policyRevision ||
	context->publishedOccurrenceKey != context->request.occurrenceKey ||
	context->publishedSourceEntryIndex != context->request.sourceEntryIndex ?
	    TRUE : FALSE;
}

static void
lod_cold_preview_note_published(BObolColdMeshPreviewContext *context)
{
    if (!context)
	return;
    context->publishedViewRevision = context->request.viewRevision;
    context->publishedPolicyRevision = context->request.policyRevision;
    context->publishedOccurrenceKey = context->request.occurrenceKey;
    context->publishedSourceEntryIndex = context->request.sourceEntryIndex;
}

/* Page callbacks are synchronous on the cold producer worker.  Publishing
 * every completed page would copy an ever-growing descriptor vector and
 * wake the owner thread O(page_count) times.  Preserve immediate first-page
 * feedback, then amortize later publication into bounded waves. */
static constexpr size_t lod_cold_spatial_publication_page_batch = 8;

static void
lod_counts_accumulate(BObolLodCounts &total,
    const BObolLodCounts &addition)
{
    const auto add = [](uint64_t left, uint64_t right) {
	return right > UINT64_MAX - left ? UINT64_MAX : left + right;
    };
    total.faceCount = add(total.faceCount, addition.faceCount);
    total.pointCount = add(total.pointCount, addition.pointCount);
    total.originalPointCount = add(total.originalPointCount,
	addition.originalPointCount);
    total.normalCount = add(total.normalCount, addition.normalCount);
    total.lineCount = add(total.lineCount, addition.lineCount);
    total.byteCount = add(total.byteCount, addition.byteCount);
}

static void
lod_use_limited_spatial_layers(BObolLodResult &result,
    const std::vector<BObolLodPresentationLayer> &layers)
{
    result.presentationLayers = layers;
    result.preparedCadGeometry.reset();
    result.preparedCadGeometryRevision = 0;
    result.counts.clear();
    for (const BObolLodPresentationLayer &layer : result.presentationLayers) {
	result.geometry.activeCut = std::max(
	    result.geometry.activeCut, layer.activeCut);
	result.residentCut = std::max(result.residentCut, layer.activeCut);
	if (layer.geometry)
	    lod_counts_accumulate(result.counts,
		bobol_cad_geometry_counts(*layer.geometry));
    }
}

static int
lod_cold_mesh_preview_cancelled(void *callbackData)
{
    const BObolColdMeshPreviewContext *context =
	static_cast<const BObolColdMeshPreviewContext *>(callbackData);
    if (!context || !context->serviceState || !context->generation)
	return 1;
    /* This callback also runs below resident.mutex.  Missing one cancellation
     * observation is safe and bounded by the producer's next callback; a
     * blocking service-lock acquisition here can instead stop both the
     * producer and its owner-thread progress pump. */
    std::unique_lock<std::mutex> lock(
	context->serviceState->mutex, std::try_to_lock);
    if (!lock.owns_lock())
	return 0;
    return context->serviceState->stopping ||
	lod_producer_cancelled_unlocked(context->serviceState,
	    context->generation, context->request) ? 1 : 0;
}

static void
lod_cold_mesh_progress(int stage, uint64_t completedUnits,
	uint64_t totalUnits, void *callbackData)
{
    BObolColdMeshPreviewContext *context =
	static_cast<BObolColdMeshPreviewContext *>(callbackData);
    if (context && context->progress)
	context->progress->publish(stage, completedUnits, totalUnits);
}

static size_t
lod_cold_preview_render_cost(
    const struct BObolMeshLodHierarchyInfo &hierarchy, int cut,
    int drawMode)
{
    if (cut < hierarchy.min_cut || cut > hierarchy.max_cut ||
	cut < 0 || static_cast<uint32_t>(cut) >= hierarchy.cut_count)
	return SIZE_MAX;

    BObolLodCounts counts;
    counts.faceCount = hierarchy.cuts[cut].face_count;
    counts.pointCount = hierarchy.cuts[cut].point_count;
    counts.normalCount = hierarchy.has_normals ? counts.pointCount : 0;
    return bobol_lod_render_cost_units(counts, drawMode, 1);
}

static void
lod_publish_cold_coverage_preview(
    unsigned long long cacheKey, const struct BObolMeshLodData *data,
    const struct BObolMeshLodHierarchyInfo *hierarchy,
    BObolColdMeshPreviewContext *context)
{
    if (!context || !context->service || !context->generation || !data ||
	!hierarchy || data->faces || data->face_count || !data->points ||
	!data->point_count)
	return;
    const SbBox3f bounds = context->request.bounds.isEmpty() ?
	SbBox3f(SbVec3f(static_cast<float>(data->bmin[X]),
		static_cast<float>(data->bmin[Y]),
		static_cast<float>(data->bmin[Z])),
	    SbVec3f(static_cast<float>(data->bmax[X]),
		static_cast<float>(data->bmax[Y]),
		static_cast<float>(data->bmax[Z]))) : context->request.bounds;
    const int coverageDrawMode = context->spatialLeafProducer ?
	BOBOL_LOD_DRAW_WIRE : context->request.drawMode;
    BObolLodCoveragePreviewBuild preview;
    if (!bobol_lod_build_coverage_preview(*data, bounds, coverageDrawMode,
	    context->renderCostAllowance, preview))
        return;
    std::shared_ptr<const Obol::PartGeometry> geometry =
	bobol_cad_build_geometry(std::move(preview.geometry),
	    "cold coverage preview");
    if (!geometry)
	return;
    const BObolLodCounts coverageCounts = preview.counts;
    const size_t coverageCost = preview.renderCost;

    BObolLodResult result;
    result.request = context->request;
    result.cacheKey = bobol_lod_cache_key(context->request);
    result.geometry.kind = BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
    result.geometry.providerId = context->request.providerId;
    result.geometry.providerVersion = context->request.providerVersion;
    result.geometry.cacheKey = context->assetKey;
    result.geometry.providerToken = cacheKey;
    result.geometry.activeCut = -1;
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = context->request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.bounds = bounds;
    result.counts = coverageCounts;
    result.preparedCadGeometry = geometry;
    /* This preview has no mutable progressive generation.  Nonzero revision
     * identifies an immutable direct PartGeometry handoff without suggesting
     * that it is a PoP resident-cut revision. */
    result.preparedCadGeometryRevision = 1;
    BObolLodPresentationLayer coverageLayer;
    coverageLayer.layerKey = "coverage";
    coverageLayer.geometry = geometry;
    coverageLayer.geometryRevision = 1;
    coverageLayer.coverage = TRUE;
    /* Keep the direct prepared handle for legacy/noncompact consumers, but
     * always carry the semantic layer marker as well.  Without it, ordinary
     * cold previews are indistinguishable from genuine source triangles in
     * convergence diagnostics. */
    result.presentationLayers.push_back(coverageLayer);
    if (context->spatialLeafProducer) {
	context->presentationLayers.clear();
	context->presentationLayers.push_back(coverageLayer);
	context->presentationRenderCost = coverageCost;
	context->publishedPageCount = 0;
	context->spatialCoverageAdmitted = TRUE;
    }
    result.terminal = FALSE;
    std::string diagnostic("cold temporary coverage preview: ");
    if (preview.pointFallback) {
	diagnostic += "points";
    } else {
	diagnostic += std::to_string(preview.cellAxis);
	diagnostic += "^3 voxels";
    }
    result.diagnostic = diagnostic.c_str();
    result.canonicalizePayload();
    if (!result.payloadIsConsistent())
	return;
    const SbBool published = context->service->tryPublishIntermediateResult(
	context->generation, std::move(result));
    if (published)
	lod_cold_preview_note_published(context);
    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	bu_log("[obol-timing] cold coverage preview: source_points=%zu "
	       "representation=%s grid=%zu cost=%zu allowance=%zu "
	       "published=%d\n", data->point_count,
	       preview.pointFallback ? "points" : "voxels",
	       preview.cellAxis, coverageCost, context->renderCostAllowance,
	       published ? 1 : 0);
}

static void
lod_publish_cold_spatial_layers(
    BObolColdMeshPreviewContext *context, unsigned long long cacheKey,
    int activeCut, const char *diagnostic)
{
    if (!context || !context->service || !cacheKey || activeCut < 0 ||
	context->presentationLayers.empty())
	return;

    BObolLodResult result;
    result.request = context->request;
    result.cacheKey = bobol_lod_cache_key(context->request);
    result.geometry.kind = BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
    result.geometry.providerId = context->request.providerId;
    result.geometry.providerVersion = context->request.providerVersion;
    result.geometry.cacheKey = context->assetKey;
    result.geometry.providerToken = cacheKey;
    result.geometry.activeCut = activeCut;
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = context->request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.bounds = context->request.bounds;
    for (const BObolLodPresentationLayer &publishedLayer :
	 context->presentationLayers)
	lod_counts_accumulate(result.counts,
	    bobol_cad_geometry_counts(*publishedLayer.geometry));
    result.presentationLayers = context->presentationLayers;
    result.terminal = FALSE;
    result.diagnostic = diagnostic ? diagnostic :
	"cold validated spatial pages; cache generation continues";
    result.canonicalizePayload();
    if (!result.payloadIsConsistent())
	return;

    if (!context->service->tryPublishIntermediateResult(
	    context->generation, std::move(result)))
	return;
    context->publishedPageCount = context->presentationLayers.size() -
	(context->presentationLayers.front().coverage ? 1u : 0u);
    lod_cold_preview_note_published(context);
}

static int
lod_cold_spatial_latest_cut(const BObolColdMeshPreviewContext *context)
{
    if (!context)
	return -1;
    for (auto layer = context->presentationLayers.rbegin();
	 layer != context->presentationLayers.rend(); ++layer)
	if (!layer->coverage && layer->activeCut >= 0)
	    return layer->activeCut;
    return -1;
}

static void
lod_publish_cold_spatial_page_impl(
    unsigned long long cacheKey,
    const struct BObolMeshLodSpatialPage *page, void *callbackData)
{
    BObolColdMeshPreviewContext *context =
	static_cast<BObolColdMeshPreviewContext *>(callbackData);
    if (!context || !context->service || !context->generation || !page ||
	!cacheKey || !context->spatialLeafProducer ||
	!context->spatialCoverageAdmitted ||
	lod_cold_mesh_preview_cancelled(context) ||
	page->cut < page->hierarchy.min_cut ||
	page->cut > page->hierarchy.max_cut || !page->data.faces ||
	!page->data.points_orig || !page->data.face_count ||
	!page->data.point_orig_count)
	return;

    lod_cold_preview_refresh_demand(context);
    const SbBool demandPublicationPending =
	lod_cold_preview_demand_unpublished(context);

    const size_t remainingAllowance =
	context->presentationRenderCost >= context->renderCostAllowance ? 0 :
	context->renderCostAllowance - context->presentationRenderCost;
    int selectedCut = -1;
    size_t selectedCost = SIZE_MAX;
    for (int cut = page->cut; cut >= page->hierarchy.min_cut; --cut) {
	const size_t candidateCost = lod_cold_preview_render_cost(
	    page->hierarchy, cut, context->request.drawMode);
	if (candidateCost <= remainingAllowance) {
	    selectedCut = cut;
	    selectedCost = candidateCost;
	    break;
	}
    }
    if (selectedCut < page->hierarchy.min_cut || selectedCost == SIZE_MAX) {
	if (demandPublicationPending)
	    lod_publish_cold_spatial_layers(context, cacheKey,
		lod_cold_spatial_latest_cut(context),
		"cold spatial pages rebound to current view");
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] cold spatial page deferred: page=%u "
		   "cut=%d minimum=%d allowance=%zu used=%zu\n",
		   page->page_id, page->cut, page->hierarchy.min_cut,
		   context->renderCostAllowance,
		   context->presentationRenderCost);
	return;
    }

    struct BObolMeshLodData selectedData = page->data;
    struct BObolMeshLodHierarchyInfo selectedHierarchy = page->hierarchy;
    selectedHierarchy.resident_cut = selectedCut;
    selectedData.face_count =
	selectedHierarchy.cuts[selectedCut].face_count;
    selectedData.point_count =
	selectedHierarchy.cuts[selectedCut].point_count;
    selectedData.point_orig_count = selectedData.point_count;
    if (selectedData.face_count > SIZE_MAX / 3u)
	return;
    selectedData.normal_count = selectedData.normals ?
	selectedData.face_count * 3u : 0;

    BObolLodProgressiveMesh pageMesh;
    if (!pageMesh.update(selectedData, selectedHierarchy, selectedCut,
	    selectedHierarchy.shaded_cull_backfaces ? TRUE : FALSE)) {
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] cold spatial page rejected: page=%u "
		   "cut=%d faces=%zu points=%zu\n", page->page_id,
		   selectedCut, selectedData.face_count,
		   selectedData.point_orig_count);
	return;
    }
    uint64_t geometryRevision = 0;
    std::shared_ptr<const Obol::PartGeometry> geometry =
	pageMesh.prepareCadGeometry(
	    context->request.drawMode, &geometryRevision);
    if (!geometry || !geometryRevision)
	return;

    BObolLodPresentationLayer layer;
    std::string layerKey("page:");
    layerKey += std::to_string(page->page_id);
    layer.layerKey = layerKey.c_str();
    layer.geometry = std::move(geometry);
    layer.geometryRevision = geometryRevision;
    layer.activeCut = selectedCut;
    layer.coverage = FALSE;
    context->presentationLayers.push_back(std::move(layer));
    context->presentationRenderCost += selectedCost;

    const size_t pageCount = context->presentationLayers.size() -
	(context->presentationLayers.front().coverage ? 1u : 0u);
    if (!demandPublicationPending && pageCount != 1 &&
	pageCount - context->publishedPageCount <
	    lod_cold_spatial_publication_page_batch)
	return;

    const size_t priorPublishedPageCount = context->publishedPageCount;
    const SbBool demandWasPending =
	lod_cold_preview_demand_unpublished(context);
    lod_publish_cold_spatial_layers(context, cacheKey, selectedCut,
	"cold validated spatial pages; cache generation continues");
    const SbBool published =
	context->publishedPageCount != priorPublishedPageCount ||
	(demandWasPending && !lod_cold_preview_demand_unpublished(context));
    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	bu_log("[obol-timing] cold spatial page: page=%u cut=%d "
	       "faces=%zu points=%zu layers=%zu cost=%zu published=%d\n",
	       page->page_id, selectedCut, selectedData.face_count,
	       selectedData.point_orig_count,
	       context->presentationLayers.size(), selectedCost,
	       published ? 1 : 0);
}

static void
lod_publish_cold_spatial_page(
    unsigned long long cacheKey,
    const struct BObolMeshLodSpatialPage *page, void *callbackData)
{
    try {
	lod_publish_cold_spatial_page_impl(cacheKey, page, callbackData);
    } catch (const std::bad_alloc &) {
	/* Live publication is opportunistic; durable cache generation remains
	 * authoritative under transient memory pressure. */
	return;
    }
}

static void
lod_publish_cold_mesh_preview_impl(
    int previewKind, unsigned long long cacheKey,
    const struct BObolMeshLodData *data,
    const struct BObolMeshLodHierarchyInfo *hierarchy,
    void *callbackData)
{
	BObolColdMeshPreviewContext *context =
	static_cast<BObolColdMeshPreviewContext *>(callbackData);
    if (previewKind == BOBOL_MESH_LOD_PREVIEW_COVERAGE_POINTS) {
	lod_publish_cold_coverage_preview(cacheKey, data, hierarchy, context);
	return;
    }
    if (previewKind != BOBOL_MESH_LOD_PREVIEW_MESH_PREFIX)
	return;
    if (!context || !context->service || !context->generation || !data ||
	!context->progressiveMesh || !hierarchy || hierarchy->min_cut < 0 ||
	!data->faces ||
	!data->points_orig || !data->face_count || !data->point_orig_count)
	return;
    const BObolLodProgressiveMeshPtr &mesh = context->progressiveMesh;
    const int residentPreviewCut =
	hierarchy->resident_cut >= hierarchy->min_cut ?
	hierarchy->resident_cut : hierarchy->min_cut;
    int previewCut = -1;
    size_t previewCost = SIZE_MAX;
    /* The cache producer materializes one globally ordered prefix before it
     * knows the renderer's current capacity.  Do not discard that complete
     * whole-object evidence merely because its richest prepared cut exceeds a
     * software renderer's first-frame allowance: select the richest contained
     * cut the same scene reservation can draw.  This is safe because every
     * lower cut is a prefix of the immutable global PoP order; it is not the
     * source-order spatial bootstrap, which has incomplete coverage. */
    for (int cut = residentPreviewCut; cut >= hierarchy->min_cut; --cut) {
	const size_t candidateCost = lod_cold_preview_render_cost(
	    *hierarchy, cut, context->request.drawMode);
	if (candidateCost <= context->renderCostAllowance) {
	    previewCut = cut;
	    previewCost = candidateCost;
	    break;
	}
    }
	if (previewCut < hierarchy->min_cut || previewCost == SIZE_MAX) {
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] cold PoP preview deferred: resident_cut=%d "
		   "allowance=%zu\n", residentPreviewCut,
		   context->renderCostAllowance);
	return;
    }
    /* The callback borrows arrays sized for residentPreviewCut.  A lower
     * selected cut is still the leading contiguous range of those arrays, but
     * the progressive-mesh validator deliberately requires the data counts to
     * match the active hierarchy cut.  Narrow the borrowed view; do not copy
     * or rebuild the prefix on the worker. */
    struct BObolMeshLodData previewData = *data;
    struct BObolMeshLodHierarchyInfo previewHierarchy = *hierarchy;
    previewHierarchy.resident_cut = previewCut;
    previewData.face_count = hierarchy->cuts[previewCut].face_count;
    previewData.point_count = hierarchy->cuts[previewCut].point_count;
    previewData.point_orig_count = previewData.point_count;
    if (previewData.normals &&
	previewData.face_count > std::numeric_limits<size_t>::max() / 3)
	return;
    previewData.normal_count = previewData.normals ?
	previewData.face_count * 3 : 0;
    if (!previewData.face_count || !previewData.point_count ||
	(previewData.normals && !previewData.normal_count))
	return;

    if (!mesh || !mesh->update(previewData, previewHierarchy, previewCut,
	hierarchy->shaded_cull_backfaces ? TRUE : FALSE)) {
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] cold PoP preview rejected by mesh "
		   "materialization: cut=%d faces=%zu points=%zu\n",
		   previewCut, previewData.face_count,
		   previewData.point_orig_count);
	return;
	}

    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    info.active_cut = previewCut;
    info.face_count = hierarchy->cuts[previewCut].face_count;
    info.point_count = hierarchy->cuts[previewCut].point_count;
    info.point_orig_count = info.point_count;
    info.normal_count = hierarchy->has_normals ? info.face_count * 3 : 0;
    info.has_faces = data->faces && info.face_count ? 1 : 0;
    info.has_points = data->points && info.point_count ? 1 : 0;
    info.has_original_points =
	previewData.points_orig && previewData.point_orig_count ? 1 : 0;
    info.has_snapped_points =
	data->points && data->points_orig != data->points ? 1 : 0;
    info.has_normals =
	previewData.normals && previewData.normal_count ? 1 : 0;
    info.shaded_cull_backfaces = hierarchy->shaded_cull_backfaces;
    VMOVE(info.bmin, data->bmin);
    VMOVE(info.bmax, data->bmax);

    int requestedCut = context->request.requestedCut;
    if (context->request.projectedPixelDiameter > 0.0f &&
	context->request.targetPixelError > 0.0f)
	requestedCut = bobol_mesh_lod_select_cut(
	    hierarchy, context->request.projectedPixelDiameter,
	    context->request.targetPixelError);
    requestedCut = std::max(hierarchy->min_cut,
	std::min(hierarchy->max_cut, requestedCut));

    BObolLodResult result = bobol_lod_result_from_mesh_lod_info(
	context->request, info, NULL);
    result.geometry.providerToken = cacheKey;
    result.geometry.cacheKey = context->assetKey;
    result.geometry.activeCut = previewCut;
    result.resolvedCut = requestedCut;
    result.residentCut = previewCut;
    result.progressiveMesh = mesh;
    result.counts.faceCount = info.face_count;
    result.counts.pointCount = info.point_count;
    result.counts.originalPointCount = info.point_orig_count;
    result.counts.normalCount = info.normal_count;
    result.bounds = mesh->bounds();
    result.hasSnappedPoints = FALSE;
    result.hasNormals = hierarchy->has_normals ? TRUE : FALSE;
    result.shadedCullBackfaces = mesh->cullBackfaces();
    /* This result replaces structural coverage while the authoritative
     * spatial cache is still being built.  It cannot satisfy convergence,
     * even when the cold request's provisional target happens to equal the
     * minimum prefix: the final result supplies the chunk visibility and
     * residency contract. */
    result.terminal = FALSE;
    result.preparedCadGeometry = mesh->prepareCadGeometry(
	context->request.drawMode, &result.preparedCadGeometryRevision);
	if (!result.preparedCadGeometry) {
	if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	    bu_log("[obol-timing] cold PoP preview could not prepare CAD "
		   "geometry: cut=%d\n", previewCut);
	return;
	}
    result.diagnostic =
	"cold minimum PoP prefix; spatial cache generation continues";
    const SbBool published =
	context->service->tryPublishIntermediateResult(
	    context->generation, std::move(result));
    if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	bu_log("[obol-timing] cold PoP preview: cut=%d/%d faces=%zu "
	       "points=%zu cost=%zu published=%d\n", previewCut,
	       residentPreviewCut, info.face_count, info.point_count,
	       previewCost, published ? 1 : 0);
}

static void
lod_publish_cold_mesh_preview(
    int previewKind, unsigned long long cacheKey,
    const struct BObolMeshLodData *data,
    const struct BObolMeshLodHierarchyInfo *hierarchy,
    void *callbackData)
{
    try {
	lod_publish_cold_mesh_preview_impl(
	    previewKind, cacheKey, data, hierarchy, callbackData);
    } catch (const std::bad_alloc &) {
	/* The preview is an early-presentation optimization.  Under memory
	 * pressure, leave its reserved task slot untouched and let the same
	 * worker finish the authoritative spatial cache result. */
	return;
    }
}

/* Return true when no later request needs to retry.  Failure to persist this
 * optional presentation metadata must never turn a usable PoP hierarchy into
 * a provider failure. */
static bool
lod_publish_draw_asset_oriented_bounds(
    struct db_i *dbip, const char *name, const BObolLodRequest &request,
    const struct BObolMeshLodHierarchyInfo &hierarchy)
{
    if (!dbip || !name || !name[0] ||
	hierarchy.oriented_bounds_valid != 1 ||
	!bobol_mesh_lod_oriented_bounds_validate(&hierarchy))
	return true;

    const int maximumCut = hierarchy.max_cut;
    const bool cutValid = maximumCut >= 0 &&
	static_cast<uint32_t>(maximumCut) < hierarchy.cut_count &&
	maximumCut < BOBOL_MESH_LOD_CUT_COUNT_MAX;
    const uint64_t faceCount = request.sourceCounts.faceCount ?
	request.sourceCounts.faceCount :
	(cutValid ? hierarchy.cuts[maximumCut].face_count : 0);
    const uint64_t pointCount = request.sourceCounts.pointCount ?
	request.sourceCounts.pointCount :
	(cutValid ? hierarchy.cuts[maximumCut].point_count : 0);
    if (!faceCount || !pointCount)
	return false;
    return bobol_draw_lod_asset_oriented_bounds_publish(
	dbip, name, faceCount, pointCount, hierarchy.quantization_min,
	hierarchy.quantization_max, hierarchy.oriented_bounds) == BRLCAD_OK;
}

BObolLodResult
BObolLodService::realizeResidentMeshLod(
    const BObolLodRequest &submittedRequest,
    const BObolMeshLodProvider &provider)
{
    /* Asset preparation is stable across camera epochs, but the final chunk
     * selection and result stamp are not.  A cold producer may run long enough
     * for the camera to move several times, so retain a mutable demand copy and
     * refresh it after the immutable hierarchy has been constructed. */
    BObolLodRequest request = submittedRequest;
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
	(void)bobol_atomic_identity_advance(resident->useRevision);
    }

    /* Publish before acquiring the resident lock.  Cache generation and
     * prefix publication deliberately serialize per asset, and on a cold
     * start that wait can otherwise look like an idle or stalled worker. */
    const std::shared_ptr<BObolLodProducerProgressRecord> producerProgress =
	lod_producer_progress_begin(this->p, provider.generation, request);
    BObolLodProducerProgressScope producerProgressScope(producerProgress);
    std::lock_guard<std::mutex> residentLock(resident->mutex);
    if (producerProgress)
	producerProgress->publish(
	    BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP, 0, 1);
    if (resident->databaseIdentity != databaseIdentity ||
	resident->name != name)
	return lod_provider_status_result(request, BOBOL_LOD_PROVIDER_STALE,
	    "resident mesh asset identity collision");

    if (!resident->lod) {
	/*
	 * Compact requests captured the validated immutable content key while
	 * their source records were streamed.  Open that payload directly.
	 * Falling back to the named API preserves the non-compact/legacy path.
	 * At distinct-asset scale this avoids one database lookup, one name-cache
	 * lookup, and one LMDB read transaction per successful warm task.
	 */
	const bool exactVariant = provider.meshAssetContentHash != 0;
	BObolColdMeshPreviewContext previewContext;
	previewContext.service = provider.service;
	previewContext.serviceState = this->p;
	previewContext.generation = provider.generation;
	previewContext.request = request;
	previewContext.assetKey = assetKey;
	previewContext.progressiveMesh = resident->mesh;
	previewContext.progress = producerProgress;
	previewContext.renderCostAllowance =
	    provider.initialRefinementCostBudget;
	struct BObolMeshLodPreviewRequest previewRequest =
	    BOBOL_MESH_LOD_PREVIEW_REQUEST_INIT;
	previewRequest.requested_cut = request.requestedCut;
	previewRequest.projected_pixel_diameter =
	    request.projectedPixelDiameter;
	previewRequest.target_pixel_error = request.targetPixelError;
	previewRequest.cancellation_callback = lod_cold_mesh_preview_cancelled;
	previewRequest.cancellation_data = &previewContext;
	previewRequest.progress_callback = lod_cold_mesh_progress;
	previewRequest.progress_data = &previewContext;
	const bool spatialLeafProducer =
	    provider.useSerializedSpatialSource ? true : false;
	previewRequest.spatial_leaf_producer = spatialLeafProducer ? 1 : 0;
	/* The producer applies the large-source threshold using its authoritative
	 * face count.  Keeping this an explicit request prevents utility/cache API
	 * callers from acquiring presentation behavior they did not ask for. */
	previewRequest.coverage_preview = 1;
	previewContext.spatialLeafProducer =
	    spatialLeafProducer ? TRUE : FALSE;
	if (spatialLeafProducer) {
	    previewRequest.spatial_page_callback =
		lod_publish_cold_spatial_page;
	    previewRequest.spatial_page_data = &previewContext;
	}
	if (!request.bounds.isEmpty()) {
	    const SbVec3f minimum = request.bounds.getMin();
	    const SbVec3f maximum = request.bounds.getMax();
	    for (size_t axis = 0; axis < 3; ++axis) {
		previewRequest.coverage_bmin[axis] = minimum[axis];
		previewRequest.coverage_bmax[axis] = maximum[axis];
	    }
	    previewRequest.coverage_bounds_valid = 1;
	}
	if (exactVariant)
	    resident->lod = bobol_mesh_lod_get_cached_prefix(
		dbip, provider.meshAssetContentHash);
	if (!resident->lod && !exactVariant)
	    resident->lod = bobol_mesh_lod_get_named_cached_prefix(
		dbip, name);
	if (producerProgress)
	    producerProgress->publish(
		BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP, 1, 1);
	if (!resident->lod && exactVariant && provider.generateBrepVariant) {
	    if (producerProgress)
		producerProgress->publish(
		    BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION, 0, 1);
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
	    if (producerProgress)
		producerProgress->publish(
		    BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION, 1, 1);
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
		if (producerProgress)
		    producerProgress->publish(
			BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION, 0, 1);
		if (stagedMatches && staged->bot) {
		    resident->lod =
			bobol_mesh_lod_cache_refresh_open(
			    dbip, name, staged->bot, &resident->status,
			    &previewRequest,
			    lod_publish_cold_mesh_preview,
			    &previewContext);
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
	    if (spatialLeafProducer ||
		!lod_raw_import_fits_address_space(request)) {
		/* Do not cross librt's bu_bomb allocation boundary merely because
		 * a cache is cold.  The explicit spatial producer and an address-space
		 * refusal both require V5's checked serialized source; an unavailable
		 * source remains a diagnosable constrained fallback rather than a
		 * provider error. */
		resident->lod = bobol_mesh_lod_cache_refresh_serialized_open(
		    dbip, name, &resident->status, &previewRequest,
		    lod_publish_cold_mesh_preview, &previewContext);
		if (!resident->lod)
		    return lod_provider_status_result(request,
			BOBOL_LOD_PROVIDER_FALLBACK,
			"resident mesh provider could not admit a checked serialized BoT import");
	    } else {
		resident->lod = bobol_mesh_lod_cache_refresh_open(
		    dbip, name, NULL, &resident->status,
		    &previewRequest,
		    lod_publish_cold_mesh_preview, &previewContext);
	    }
		    refreshResult = resident->lod ?
			BRLCAD_OK : BRLCAD_ERROR;
		} else {
		    return lod_provider_status_result(request,
			BOBOL_LOD_PROVIDER_CACHE_MISS,
			"resident mesh provider has no exact staged variant");
		}
		if (producerProgress &&
		    producerProgress->stage.load(std::memory_order_acquire) ==
			BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION)
		    producerProgress->publish(
			BOBOL_LOD_PRODUCER_STAGE_SOURCE_PREPARATION, 1, 1);
		if (refreshResult != BRLCAD_OK)
		    return lod_provider_status_result(request,
			BOBOL_LOD_PROVIDER_CACHE_MISS,
			"resident mesh provider could not refresh cache entry");
		if (getenv("BOBOL_DRAW_TIMING_VERBOSE"))
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
	    lod_cold_preview_refresh_demand(&previewContext);
	    request = previewContext.request;
	    resident->limitedSpatialLayers =
		std::move(previewContext.presentationLayers);
	}
	if (!resident->lod)
	    return lod_provider_status_result(request,
		resident->status.stale_cache_entry ?
		    BOBOL_LOD_PROVIDER_STALE :
		    BOBOL_LOD_PROVIDER_CACHE_MISS,
		"resident mesh provider has no cache payload");
    }
	if (producerProgress &&
	    producerProgress->stage.load(std::memory_order_acquire) ==
		BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP)
	    producerProgress->publish(
		BOBOL_LOD_PRODUCER_STAGE_CACHE_LOOKUP, 1, 1);

    struct BObolMeshLodHierarchyInfo hierarchy =
	BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
    if (!bobol_mesh_lod_hierarchy_info_get(resident->lod, &hierarchy))
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_CACHE_MISS,
	    "resident mesh provider loaded no hierarchy metadata");
    if (!resident->orientedBoundsPublished)
	resident->orientedBoundsPublished =
	    lod_publish_draw_asset_oriented_bounds(
		dbip, name, request, hierarchy);
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
    const bool sourceLimited = bobol_mesh_lod_source_limited(resident->lod);
    std::vector<uint32_t> requiredChunks = request.requiredChunks;
    const bool chunked = hierarchy.chunks && hierarchy.chunk_count;
    if (chunked && requiredChunks.empty()) {
	if (!request.spatialProjectionValid ||
	    !bobol_lod_visible_chunks(hierarchy, request.localToRoot,
		request.viewProjection, requiredChunks)) {
	    /* Non-view API consumers conservatively request the complete leaf. */
	    requiredChunks.resize(hierarchy.chunk_count);
	    for (uint32_t chunk = 0; chunk < hierarchy.chunk_count; ++chunk)
		requiredChunks[chunk] = chunk;
	}
    }
    for (size_t i = 0; chunked && i < requiredChunks.size(); ++i) {
	if (requiredChunks[i] >= hierarchy.chunk_count ||
	    (i && requiredChunks[i] <= requiredChunks[i - 1]))
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_ERROR,
		"resident mesh provider received an invalid chunk set");
    }
    if (chunked && requiredChunks.empty()) {
	BObolLodResult result;
	result.request = request;
	result.request.requiredChunks.clear();
	result.cacheKey = bobol_lod_cache_key(request);
	result.geometry.kind = BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
	result.geometry.providerId = request.providerId;
	result.geometry.providerVersion = request.providerVersion;
	result.geometry.cacheKey = assetKey;
	result.geometry.providerToken =
	    bobol_mesh_lod_cache_key_get(resident->lod);
	result.resultKind = BOBOL_LOD_RESULT_MESH;
	result.qualityTier = request.qualityTier;
	result.providerStatus = BOBOL_LOD_PROVIDER_READY;
	result.bounds = request.bounds;
	result.resolvedCut = requestedCut;
	result.hasNormals = hierarchy.has_normals ? TRUE : FALSE;
	result.shadedCullBackfaces =
	    hierarchy.shaded_cull_backfaces ? TRUE : FALSE;
	result.terminal = TRUE;
	if (sourceLimited && !resident->limitedSpatialLayers.empty()) {
	    result.residentCut = hierarchy.resident_cut;
	    result.memoryLimited = TRUE;
	    result.diagnostic =
		"bounded spatial coverage retained; no resident page is visible";
	    lod_use_limited_spatial_layers(
		result, resident->limitedSpatialLayers);
	} else if (resident->mesh && resident->mesh->isValid() &&
	    resident->mesh->hasSpatialClusters()) {
	    /* An empty visible-page set is an exact view-local population.  It is
	     * not cancellation: returning CANCELLED asks authentication to replay
	     * the same current demand forever.  Reuse the retained hierarchy as
	     * proof that this is a spatial asset, but publish no renderer layers
	     * and charge no triangles. */
	    result.progressiveMesh = resident->mesh;
	    result.geometry.activeCut = requestedCut;
	    result.residentCut = resident->mesh->residentCut();
	    result.bounds = resident->mesh->bounds();
	    result.counts.clear();
	    result.presentationLayers.clear();
	    result.diagnostic = "view-local spatial page demand is empty";
	} else {
	    /* A cold asset without a retained spatial generation cannot yet prove
	     * the empty-page interpretation.  Make the exceptional condition
	     * terminal for this exact demand rather than reporting cancellation,
	     * whose retry disposition has no state change capable of succeeding. */
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_ERROR,
		"resident mesh provider cannot certify empty spatial demand");
	}
	result.presentationAdmissionCertified =
	    provider.presentationAdmissionCertified;
	result.presentationAdmissionViewRevision =
	    provider.presentationAdmissionViewRevision;
	result.presentationAdmissionPolicyRevision =
	    provider.presentationAdmissionPolicyRevision;
	result.presentationAdmissionCut = result.geometry.activeCut;
	result.canonicalizePayload();
	return result;
    }
    resident->publishedMinimumCut.store(
	hierarchy.min_cut, std::memory_order_relaxed);

    const size_t priorResidentBytes =
	resident->publishedBytes.load(std::memory_order_relaxed);
    const size_t priorBackingBytes =
	resident->publishedBackingPrefixBytes.load(
	    std::memory_order_relaxed);
    const size_t priorStableBytes = lod_resident_asset_stable_bytes(
	priorResidentBytes, priorBackingBytes);
    int currentCut = bobol_mesh_lod_current_cut(resident->lod);
    const int publishedCut =
	resident->mesh && resident->mesh->isValid() ?
	    (chunked ? resident->mesh->residentCutForChunks(requiredChunks) :
		resident->mesh->residentCut()) : -1;
    const int deliveryCut = provider.useCurrentDrawCut ?
	provider.currentDrawCut :
	(publishedCut >= 0 ? publishedCut : currentCut);
    int residentTarget = requestedCut;
    if (requestedCut >= 0)
	residentTarget = lod_progressive_delivery_cut(hierarchy, requestedCut,
	    deliveryCut, request.drawMode, provider,
	    chunked ? &requiredChunks : NULL);
    if (residentTarget < hierarchy.min_cut)
	residentTarget = hierarchy.min_cut;
    if (residentTarget > hierarchy.max_cut)
	residentTarget = hierarchy.max_cut;
    /* Ordinary first publication keeps the bounded useful prefix selected
     * above.  Only replacement of an already visible adaptive representation
     * bypasses that staging and realizes the complete requested handoff. */
    const bool atomicFirstHandoff =
	provider.atomicRepresentationHandoff &&
	publishedCut < hierarchy.min_cut &&
	!provider.resetExisting && !provider.useForcedCut &&
	requestedCut >= hierarchy.min_cut;
    if (atomicFirstHandoff)
	residentTarget = std::min(hierarchy.max_cut, requestedCut);

    /*
     * Changing the visible page set is an atomic presentation transition.
     * In particular, zooming out may add pages which are not resident yet
     * even though the occurrence is already being drawn at a rich cut.  A
     * stale delivery limit must not add those pages at the hierarchy minimum:
     * residentCutForChunks() would then make that minimum the common frontier
     * and the owner thread would briefly replace the whole mesh by a handful
     * of slab-like triangles.
     *
     * Preserve the poorer of the incumbent draw cut and the new physical
     * pixel target.  This still permits a deliberate interactive downgrade,
     * but requires every newly demanded page to reach that downgrade cut
     * before the shared generation is changed.  The old immutable prepared
     * geometry remains drawable until this worker can make the transition in
     * one publication.
     */
    int presentationContinuityCut = hierarchy.min_cut - 1;
    if (!provider.resetExisting && !provider.useForcedCut &&
	provider.useCurrentDrawCut &&
	provider.currentDrawCut >= hierarchy.min_cut) {
	presentationContinuityCut = std::min(
	    provider.currentDrawCut, requestedCut);
	/* A current scene allocation may deliberately choose a coarser cut
	 * while changing the visible spatial page set.  That allocation is the
	 * authority which makes the atomic replacement affordable; continuity
	 * preserves the old framebuffer until publication, not its more
	 * expensive scalar cut in the replacement generation. */
	if (provider.presentationAdmissionCertified &&
	    provider.usePresentationCutLimit &&
	    provider.presentationCutLimit >= hierarchy.min_cut)
	    presentationContinuityCut = std::min(
		presentationContinuityCut,
		provider.presentationCutLimit);
	residentTarget = std::max(
	    residentTarget, presentationContinuityCut);
    }
    /* Publish the exact prefix this task loaded.  First publication uses the
     * scene-admitted useful-prefix allowance; representation replacement uses
     * the atomic handoff above. */
    int drawTarget = residentTarget;
    if (provider.usePresentationCutLimit &&
	provider.presentationCutLimit >= hierarchy.min_cut)
	drawTarget = std::min(drawTarget,
	    provider.presentationCutLimit);
    if (presentationContinuityCut >= hierarchy.min_cut)
	drawTarget = std::max(drawTarget, presentationContinuityCut);
    int loadTarget = residentTarget;
    if (provider.compactResident && publishedCut >= 0 &&
	residentTarget < publishedCut) {
	/* A cut is not a bounded unit of memory: one Lucy hierarchy step can
	 * add millions of faces.  Stable reclamation therefore retains exactly
	 * the pixel-demanded prefix. */
	loadTarget = residentTarget;
    } else if (publishedCut >= 0 && residentTarget <= publishedCut &&
	!provider.resetExisting) {
	loadTarget = publishedCut;
    }

    BObolResidentMeshGrowthReservation growthReservation(this->p);
    SbBool memoryLimited = FALSE;
    loadTarget = growthReservation.admit(*resident, hierarchy,
	loadTarget, publishedCut, priorStableBytes, memoryLimited,
	chunked ? &requiredChunks : NULL);
    if (presentationContinuityCut >= hierarchy.min_cut &&
	loadTarget < presentationContinuityCut) {
	return lod_provider_status_result(request,
	    BOBOL_LOD_PROVIDER_CANCELLED,
	    "resident mesh provider could not admit an atomic page-set transition");
    }
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
	(chunked ? resident->mesh->canDrawChunksAtCut(
	    requiredChunks, loadTarget) :
	    resident->mesh->canDrawCut(loadTarget));
    const bool loadNeeded =
	provider.resetExisting || !retainedTargetDrawable ||
	(publishedCut >= 0 && loadTarget != publishedCut);
	if (producerProgress && loadNeeded)
	    producerProgress->publish(
		BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION, 0,
		1);
    int64_t prefixLoadMicroseconds = 0;
    int64_t generationBuildMicroseconds = 0;
    int64_t preparedGeometryMicroseconds = 0;
    int64_t directPayloadMicroseconds = 0;
    int residentCut = publishedCut;
    struct BObolMeshLodInfo info = BOBOL_MESH_LOD_INFO_INIT;
    const auto populateInfoFromRetainedMesh = [&]() {
	bobol_mesh_lod_info_init(&info);
	info.active_cut = residentCut;
	if (chunked) {
	    uint64_t faces = 0;
	    uint64_t points = 0;
	    for (uint32_t chunk : requiredChunks) {
		faces += hierarchy.chunks[chunk].cuts[residentCut].face_count;
		points += hierarchy.chunks[chunk].cuts[residentCut].point_count;
	    }
	    info.face_count = static_cast<size_t>(std::min<uint64_t>(
		faces, SIZE_MAX));
	    info.point_count = static_cast<size_t>(std::min<uint64_t>(
		points, SIZE_MAX));
	} else {
	    info.face_count = resident->mesh->faceCount(residentCut);
	    info.point_count = resident->mesh->pointCount(residentCut);
	}
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
	if (chunked) {
	    const int64_t generationBuildStarted = bu_gettime();
	    if (!resident->mesh->updateChunksFromCache(
		    resident->lod, hierarchy, requiredChunks, loadTarget,
		    hierarchy.shaded_cull_backfaces ? TRUE : FALSE))
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider could not publish chunk prefixes");
	    generationBuildMicroseconds = std::max<int64_t>(
		0, bu_gettime() - generationBuildStarted);
	    residentCut = resident->mesh->residentCutForChunks(requiredChunks);
	    if (residentCut < hierarchy.min_cut)
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider published an incomplete chunk set");
	    populateInfoFromRetainedMesh();
	}
	/* Persistent PoP records are already split by activation cut.  Once a
	 * quiet compaction has released the cache reader's duplicate prefix, grow
	 * the immutable renderer generation from only the missing cache suffix.
	 * Corner-normal vertex splitting still needs whole-prefix context and uses
	 * the conservative cumulative fallback. */
	SbBool suffixExtended = FALSE;
	if (!chunked && !provider.resetExisting &&
		publishedCut >= hierarchy.min_cut &&
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
	if (!chunked && !suffixExtended) {
	    const int64_t prefixLoadStarted = bu_gettime();
	    residentCut = bobol_mesh_lod_load_resident_cut(
		resident->lod, loadTarget,
		provider.resetExisting ? 1 : 0);
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
	    this->p->residentMeshCacheLoads.fetch_add(
		1, std::memory_order_relaxed);
	}
    } else {
	this->p->residentMeshHits.fetch_add(
	    1, std::memory_order_relaxed);
	populateInfoFromRetainedMesh();
    }
	if (producerProgress && loadNeeded)
	    producerProgress->publish(
		BOBOL_LOD_PRODUCER_STAGE_PREFIX_MATERIALIZATION,
		1, 1);
    resident->publishedResidentCut.store(
	residentCut, std::memory_order_relaxed);

    int drawCut = drawTarget;
    if (drawCut < hierarchy.min_cut)
	drawCut = hierarchy.min_cut;
    if (drawCut > residentCut)
	drawCut = residentCut;

    BObolLodResult result =
	bobol_lod_result_from_mesh_lod_info(request, info, &resident->status);
    result.request.requiredChunks = requiredChunks;
    result.resolvedCut = requestedCut;
    result.geometry.cacheKey = assetKey;
    result.geometry.activeCut = drawCut;
    result.progressiveMesh = resident->mesh;
    result.residentCut = residentCut;
    result.residentAdmissionRevision =
	growthReservation.revision();
    if (chunked) {
	result.counts.clear();
	for (uint32_t chunk : requiredChunks) {
	    result.counts.faceCount +=
		hierarchy.chunks[chunk].cuts[drawCut].face_count;
	    result.counts.pointCount +=
		hierarchy.chunks[chunk].cuts[drawCut].point_count;
	}
    } else {
	result.counts.faceCount = resident->mesh->faceCount(drawCut);
	result.counts.pointCount = resident->mesh->pointCount(drawCut);
    }
    result.counts.originalPointCount = result.counts.pointCount;
    result.counts.normalCount = info.has_normals ?
	result.counts.faceCount * 3 : 0;
    result.bounds = resident->mesh->bounds();
    result.hasSnappedPoints = FALSE;
    result.hasNormals = hierarchy.has_normals ? TRUE : FALSE;
    result.shadedCullBackfaces = resident->mesh->cullBackfaces();
    const bool transientLimitReached = provider.transientMemoryLimited &&
	provider.useDeliveryCutLimit &&
	drawCut >= provider.deliveryCutLimit;
    result.terminal = drawCut >= std::max(hierarchy.min_cut,
	std::min(hierarchy.max_cut, requestedCut)) || transientLimitReached ?
	TRUE : FALSE;
    /* A representation-band admission cap is not a reason to stop walking
     * the admitted band's PoP prefix.  Publish it only on that band's
     * terminal result; treating it as an immediate resident-memory failure
     * would strand a BREP at its first coarse cut. */

    result.memoryLimited = memoryLimited || sourceLimited ||
	provider.transientMemoryLimited ||
	(provider.brepVariantMemoryLimited && result.terminal);
    if (transientLimitReached)
	result.diagnostic =
	    "view detail constrained by transient working-set limit";

    /*
     * Prepare the renderer target on this worker for every CAD consumer, not
     * only compact occurrences.  Source-wide database nodes can adopt the
     * same immutable allocation and avoid rebuilding/copying mesh vectors on
     * the GUI thread.  The vector payload remains populated for direct
     * SoBRLMeshShape clients, which are a supported custom-scene API rather
     * than part of the compact CAD occurrence route.
     */
    const int64_t preparedGeometryStarted = bu_gettime();
    if (chunked) {
	std::vector<BObolLodChunkCut> presentationCuts;
	presentationCuts.reserve(requiredChunks.size());
	for (uint32_t chunkId : requiredChunks)
	    presentationCuts.push_back({chunkId, drawCut});
	if (!resident->mesh->prepareCadPresentationLayers(
		request.drawMode, presentationCuts, request.normalStyle,
		request.normalCreaseAngle,
		result.presentationLayers))
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_ERROR,
		"resident mesh provider could not prepare spatial renderer pages");
	if (result.presentationLayers.empty()) {
	    std::vector<uint32_t> populatedChunks;
	    if (!resident->mesh->populatedChunkIdsAtCut(
		    requiredChunks, drawCut, populatedChunks) ||
		!populatedChunks.empty())
		return lod_provider_status_result(request,
		    BOBOL_LOD_PROVIDER_ERROR,
		    "resident mesh provider omitted populated spatial renderer pages");
	    /* bobol_lod_result_from_mesh_lod_info() conservatively classifies a
	     * zero-count mesh as a cache miss.  For a spatial hierarchy, however,
	     * zero view-local faces at this global cut are an exact drawable state.
	     * The retained progressive generation proves the distinction. */
	    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
	    result.diagnostic =
		"view-local spatial page population is empty at the active cut";
	    result.canonicalizePayload();
	}
    } else {
	result.preparedCadGeometry =
	    resident->mesh->prepareCadGeometry(
		request.drawMode, &result.preparedCadGeometryRevision);
    }
    if (request.normalStyle == BOBOL_LOD_NORMAL_FLAT) {
	result.counts.normalCount = 0;
	result.hasNormals = FALSE;
    } else if (request.normalStyle == BOBOL_LOD_NORMAL_SMOOTH && chunked) {
	uint64_t normalCount = 0;
	for (const BObolLodPresentationLayer &layer :
		result.presentationLayers) {
	    if (!layer.geometry || !layer.geometry->shaded)
		continue;
	    const size_t count = layer.geometry->shaded->normals.size();
	    normalCount = count > UINT64_MAX - normalCount ?
		UINT64_MAX : normalCount + static_cast<uint64_t>(count);
	}
	result.counts.normalCount = normalCount;
	result.hasNormals = normalCount ? TRUE : FALSE;
    }
    preparedGeometryMicroseconds = std::max<int64_t>(
	0, bu_gettime() - preparedGeometryStarted);
    /* A capacity-limited spatial cache contains independently valid local
     * pages but not complete source coverage.  Preserve the coverage layer
     * and every page admitted by this task; replacing them with the resident
     * seed page alone would turn a recognizable whole-object preview into a
     * visibly chopped mesh at task completion. */
    if (sourceLimited && !resident->limitedSpatialLayers.empty()) {
	lod_use_limited_spatial_layers(
	    result, resident->limitedSpatialLayers);
	result.terminal = TRUE;
	result.diagnostic =
	    "bounded spatial coverage and validated pages; durable cache capacity limited";
    }
    result.presentationAdmissionCertified =
	provider.presentationAdmissionCertified;
    result.presentationAdmissionViewRevision =
	provider.presentationAdmissionViewRevision;
    result.presentationAdmissionPolicyRevision =
	provider.presentationAdmissionPolicyRevision;
    result.presentationAdmissionCut = result.geometry.activeCut;
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
    if (residentBytes != priorResidentBytes)
	lod_resident_mesh_revision_advance(
	    this->p->residentMeshRevision);
    lod_resident_mesh_accounting_replace(this->p,
	priorResidentBytes, residentBytes,
	priorBackingBytes, backingBytes);
    const size_t stableBytes =
	lod_resident_asset_stable_bytes(residentBytes, backingBytes);
    if (stableBytes < priorStableBytes)
	lod_resident_mesh_revision_advance(
	    this->p->residentMeshAdmissionRevision);
    /* Publish exact totals before making this reservation available to a
     * peer.  Otherwise many workers can briefly admit against the gap
     * between estimated release and exact accounting. */
    growthReservation.release();

    if (request.occurrenceKey.getLength() == 0) {
	const int64_t directPayloadStarted = bu_gettime();
	if (!resident->mesh->copyCut(result.mesh, drawCut)) {
	    return lod_provider_status_result(request,
		BOBOL_LOD_PROVIDER_ERROR,
		"resident mesh provider could not materialize direct shape payload");
	}
	directPayloadMicroseconds = std::max<int64_t>(
	    0, bu_gettime() - directPayloadStarted);
    }

    if (loadNeeded && getenv("BOBOL_DRAW_TIMING_VERBOSE"))
	bu_log("[obol-timing] resident prefix %-24s cut=%d->%d "
	       "load=%8.1f ms generation=%8.1f ms prepare=%8.1f ms "
	       "direct=%8.1f ms\n",
	       name, publishedCut, residentCut,
	       prefixLoadMicroseconds / 1000.0,
	       generationBuildMicroseconds / 1000.0,
	       preparedGeometryMicroseconds / 1000.0,
	       directPayloadMicroseconds / 1000.0);

    const char *traceFilter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (traceFilter && traceFilter[0] &&
	(strstr(name, traceFilter) ||
	 (request.objectPath.getLength() > 0 &&
	  strstr(request.objectPath.getString(), traceFilter)))) {
	bu_log("BObol resident LoD trace object=%s submitted_cut=%d "
	       "resolved_cut=%d draw_cut=%d resident_cut=%d "
	       "pixels=%.9g target_error=%.9g faces=%zu points=%zu "
	       "memory_limited=%d admission=%llu "
	       "asset_revision=%llu load=%d compact=%d terminal=%d\n",
	       name, request.requestedCut, result.resolvedCut, drawCut,
	       residentCut, request.projectedPixelDiameter,
	       request.targetPixelError,
	       result.counts.faceCount, result.counts.pointCount,
	       result.memoryLimited ? 1 : 0,
	       static_cast<unsigned long long>(
		   result.residentAdmissionRevision),
	       static_cast<unsigned long long>(resident->mesh->revision()),
	       loadNeeded ? 1 : 0, provider.compactResident ? 1 : 0,
	       result.terminal ? 1 : 0);
	if (getenv("BOBOL_LOD_TRACE_HIERARCHY")) {
	    for (int cut = hierarchy.min_cut;
		 cut <= hierarchy.max_cut; ++cut)
		bu_log("BObol resident hierarchy object=%s cut=%d "
		       "faces=%zu points=%zu object_error=%.17g%s\n",
		       name, cut, hierarchy.cuts[cut].face_count,
		       hierarchy.cuts[cut].point_count,
		       hierarchy.cuts[cut].object_error,
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
	    current->second.revision == demandRevision &&
	    lod_resident_consumer_snapshot_current(current->second)) {
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
	snapshot.snapshotRevision = demandRevision;
	for (const BObolLodResidentDemand &demand : demands) {
	    if (demand.assetKey.getLength() == 0 || demand.cut < 0)
		continue;
	    BObolResidentMeshDemandValue &value =
		snapshot.assets[demand.assetKey.getString()];
	    value.cut = std::max(value.cut, demand.cut);
	    value.channelMask |= demand.channelMask & 3u;
	    std::vector<uint32_t> merged;
	    merged.reserve(value.chunkIds.size() + demand.chunkIds.size());
	    std::set_union(value.chunkIds.begin(), value.chunkIds.end(),
		demand.chunkIds.begin(), demand.chunkIds.end(),
		std::back_inserter(merged));
	    value.chunkIds.swap(merged);
	}
	snapshot.planning = TRUE;
	snapshot.planningCursor = 0;
	snapshot.planningProjectedResidentBytes =
	    lod_resident_stable_bytes(this->p);
	snapshot.planningCandidateCount = 0;
	snapshot.completedCandidateCount = 0;
	snapshot.completedPlanRevision = 0;
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
	    const bool snapshotChanged =
		!lod_resident_consumer_snapshot_current(current) ||
		current.revision != demandRevision ||
		current.assets != snapshot.assets;
	    current = std::move(snapshot);
	    current.residentRevision = 0;
	    if (snapshotChanged)
		lod_resident_demand_epoch_advance(this->p);
	}

	/* Record the complete demand while refinement is active, but wait to
	 * queue trims.  A completed provider result remains part of that active
	 * transaction until its owner drains it: its immutable geometry was built
	 * from the preceding demand and may still be rebased to the current one.
	 * Trimming the shared asset between completion and publication made the
	 * rebase validate one generation and publish another.  A subsequent stable
	 * pump observes the changed resident revision and constructs the plan from
	 * the newly installed demand. */
	if (!this->p->pending.empty() || this->p->inFlight != 0 ||
	    !this->p->results.empty())
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
	    std::map<uint32_t, int> aggregateChunkCuts;
	    SbBool demanded = FALSE;
	    for (const auto &consumer :
		    this->p->residentMeshConsumerDemands) {
		if (!lod_resident_consumer_snapshot_current(consumer.second))
		    continue;
		const auto demand =
		    consumer.second.assets.find(residentEntry.first);
		if (demand == consumer.second.assets.end())
		    continue;
		demanded = TRUE;
		aggregate.cut = std::max(
		    aggregate.cut, demand->second.cut);
		aggregate.channelMask |=
		    demand->second.channelMask;
		for (uint32_t chunkId : demand->second.chunkIds) {
		    const auto found = aggregateChunkCuts.find(chunkId);
		    if (found == aggregateChunkCuts.end())
			aggregateChunkCuts.emplace(
			    chunkId, demand->second.cut);
		    else
			found->second = std::max(
			    found->second, demand->second.cut);
		}
	    }
	    const SbBool residentMemoryPressure =
		this->p->maxResidentMeshBytes != SIZE_MAX &&
		current.planningProjectedResidentBytes >
		    this->p->maxResidentMeshBytes ? TRUE : FALSE;
	    const SbBool evict =
		!demanded && residentMemoryPressure ? TRUE : FALSE;
	    const int demandedCut = aggregate.cut < 0 ?
		minimumCut : std::max(minimumCut, aggregate.cut);
	    /* Resident geometry is also the latency cache for a later camera
	     * expansion.  Presentation demand may choose a cheaper active cut,
	     * but only actual capacity pressure authorizes shortening the shared
	     * immutable prefix.  Backing arrays remain independently reclaimable
	     * below because they can be re-read without rebuilding that prefix. */
	    const int targetCut = residentMemoryPressure ?
		demandedCut : residentCut;
	    std::vector<BObolLodChunkCut> residentChunkCuts;
	    if (resident->mesh)
		resident->mesh->residentChunkCuts(residentChunkCuts);
	    std::vector<BObolLodChunkCut> targetChunkCuts;
	    if (!residentChunkCuts.empty()) {
		if (!residentMemoryPressure) {
		    targetChunkCuts = residentChunkCuts;
		} else {
		/* A resident-memory trim is admissible only after every currently
		 * demanded page has reached its requested cut.  Provider growth owns
		 * missing pages and suffixes.  Compacting the intersection while that
		 * growth was pending produced an internally valid but non-drawable
		 * generation, which later replaced a coherent visible frame. */
		SbBool demandedWorkingSetComplete = TRUE;
		for (const auto &desired : aggregateChunkCuts) {
		    const auto residentChunk = std::lower_bound(
			residentChunkCuts.begin(), residentChunkCuts.end(),
			desired.first,
			[](const BObolLodChunkCut &entry, uint32_t id) {
			    return entry.chunkId < id;
			});
		    if (residentChunk == residentChunkCuts.end() ||
			residentChunk->chunkId != desired.first ||
			residentChunk->cut < desired.second) {
			demandedWorkingSetComplete = FALSE;
			break;
		    }
		}
		/* Compaction can only discard or shorten resident pages; a provider
		 * request performs growth.  Under pressure, preserve independently
		 * demanded pages and keep every previously seen page at a recognizable
		 * coverage floor until an entirely undemanded asset must be evicted. */
		if (!demanded) {
		    for (const BObolLodChunkCut &residentChunk :
			residentChunkCuts)
			aggregateChunkCuts[residentChunk.chunkId] = minimumCut;
		} else if (aggregateChunkCuts.empty()) {
		    /* Empty page lists are invalid for chunked consumer demands.  Fail
		     * conservatively instead of turning an integration bug into an
		     * accidental whole-leaf eviction. */
		    for (const BObolLodChunkCut &residentChunk :
			residentChunkCuts)
			aggregateChunkCuts[residentChunk.chunkId] =
			    residentChunk.cut;
		} else {
		    /* The pressure floor is intentionally richer than the mathematical
		     * minimum, whose isolated spatial cells may be slab-like rather than
		     * a recognizable whole-surface approximation. */
		    static const double pressureCoverageReferencePixels = 128.0;
		    static const double coverageTargetPixelError = 1.0;
		    int coverageFloorCut = resident->mesh->cutForScreenError(
			pressureCoverageReferencePixels,
			coverageTargetPixelError);
		    coverageFloorCut = std::max(minimumCut,
			coverageFloorCut < 0 ? minimumCut : coverageFloorCut);
		    for (const BObolLodChunkCut &residentChunk :
			residentChunkCuts) {
			const int retainedCut =
			    std::min(residentChunk.cut, coverageFloorCut);
			aggregateChunkCuts.emplace(
			    residentChunk.chunkId, retainedCut);
		    }
		}
		for (const auto &desired : aggregateChunkCuts) {
		    const auto residentChunk = std::lower_bound(
			residentChunkCuts.begin(), residentChunkCuts.end(),
			desired.first,
			[](const BObolLodChunkCut &entry, uint32_t id) {
			    return entry.chunkId < id;
			});
		    if (residentChunk == residentChunkCuts.end() ||
			residentChunk->chunkId != desired.first)
			continue;
		    BObolLodCounts population;
		    const std::vector<uint32_t> oneChunk = {desired.first};
		    int desiredCut = std::max(minimumCut, desired.second);
		    /* A page may begin after the hierarchy-wide minimum.  Its
		     * coverage floor is the first populated cut, not necessarily
		     * minimumCut. */
		    while (desiredCut <= residentChunk->cut &&
			resident->mesh->hierarchyCountsForChunksAtCut(
			    oneChunk, desiredCut, FALSE, &population) &&
			!population.faceCount)
			desiredCut++;
		    if (desiredCut > residentChunk->cut ||
			!resident->mesh->hierarchyCountsForChunksAtCut(
			    oneChunk, desiredCut, FALSE, &population) ||
			!population.faceCount) {
			/* Keep the page identity and its prior prefix.  A page which
			 * contributes no faces at this common cut is still part of the
			 * coherent spatial working set and may begin at a later cut. */
			targetChunkCuts.push_back(*residentChunk);
			continue;
		    }
		    targetChunkCuts.push_back({desired.first,
			std::min(residentChunk->cut, desiredCut)});
		}
		if (!demandedWorkingSetComplete)
		    targetChunkCuts = residentChunkCuts;
		/* An empty retained set cannot form a drawable immutable generation;
		 * visibility withdrawal should remove the demand instead. */
		if (targetChunkCuts.empty())
		    targetChunkCuts = residentChunkCuts;
		}
	    }
	    const SbBool workingSetChanged = !residentChunkCuts.empty() ?
		residentChunkCuts != targetChunkCuts : targetCut < residentCut;
	    if (getenv("BOBOL_LOD_TRACE_COMPACTION")) {
		BObolLodCounts demandedCounts;
		std::vector<uint32_t> targetChunks;
		targetChunks.reserve(aggregateChunkCuts.size());
		for (const auto &entry : aggregateChunkCuts)
		    targetChunks.push_back(entry.first);
		const SbBool counted = resident->mesh &&
		    !targetChunks.empty() ?
		    resident->mesh->hierarchyCountsForChunksAtCut(
			targetChunks, targetCut, FALSE,
			&demandedCounts) : FALSE;
		bu_log("BObol resident compaction plan asset=%s "
		       "demand_revision=%llu demanded=%d cut=%d "
		       "demand_chunks=%zu resident_chunks=%zu "
		       "target_chunks=%zu faces=%zu counted=%d changed=%d\n",
		       residentEntry.first.c_str(),
		       static_cast<unsigned long long>(demandRevision),
		       demanded ? 1 : 0, targetCut,
		       aggregateChunkCuts.size(), residentChunkCuts.size(),
		       targetChunkCuts.size(), demandedCounts.faceCount,
		       counted ? 1 : 0, workingSetChanged ? 1 : 0);
	    }
	    if (!evict && !workingSetChanged &&
		targetCut >= residentCut && !hasReloadableBacking) {
		this->p->residentMeshCompactionTargets.erase(
		    residentEntry.first);
		continue;
	    }
	    if (current.planningCandidateCount != SIZE_MAX)
		current.planningCandidateCount++;
	    BObolResidentMeshCompactionTarget &target =
		this->p->residentMeshCompactionTargets[
		    residentEntry.first];
	    target.cut = targetCut;
	    target.channelMask = aggregate.channelMask;
	    target.chunkCuts = std::move(targetChunkCuts);
	    target.evict = evict;
	    target.useRevision = resident->useRevision.load(
		std::memory_order_relaxed);
	    target.demandEpoch = this->p->residentMeshDemandEpoch;
	    target.revision = bobol_nonzero_identity_take(
		this->p->nextResidentMeshCompactionTargetRevision);
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
	    current.completedCandidateCount = current.planningCandidateCount;
	    current.completedPlanRevision = current.revision;
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
BObolLodService::invalidateResidentMeshConsumer(uint64_t consumerId)
{
    if (!consumerId)
	return;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const auto consumer =
	this->p->residentMeshConsumerDemands.find(consumerId);
    if (consumer == this->p->residentMeshConsumerDemands.end() ||
	!lod_resident_consumer_snapshot_current(consumer->second))
	return;

    consumer->second.snapshotRevision = 0;
    consumer->second.residentRevision = 0;
    consumer->second.planning = FALSE;
    consumer->second.planningCursor = 0;
    consumer->second.planningProjectedResidentBytes = 0;
    consumer->second.planningCandidateCount = 0;
    consumer->second.completedCandidateCount = 0;
    consumer->second.completedPlanRevision = 0;
    consumer->second.planningResidentRevision = 0;
    lod_resident_demand_epoch_advance(this->p);
    lod_discard_resident_compaction_results_unlocked(this->p, consumerId);
    this->p->workerCv.notify_all();
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
    (void)bobol_atomic_identity_advance(found->second->useRevision);
}

void
BObolLodService::releaseResidentMeshConsumer(uint64_t consumerId)
{
    if (!consumerId)
	return;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (this->p->residentMeshConsumerDemands.erase(consumerId))
	lod_resident_demand_epoch_advance(this->p);
    lod_discard_resident_compaction_results_unlocked(this->p, consumerId);
    this->p->workerCv.notify_all();
}

uint64_t
BObolLodService::beginGeneration(void)
{
    std::lock_guard<std::mutex> lock(this->p->mutex);

    bobol_identity_advance(this->p->nextGeneration);
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

	for (auto producer = this->p->sharedProducers.begin();
	     producer != this->p->sharedProducers.end();) {
	    producer->second.leases.erase(generation);
	    if (!producer->second.leases.empty() ||
		producer->second.state == BObolSharedProducerState::BUILDING) {
		++producer;
		continue;
	    }
	    const uint64_t taskId = producer->second.taskId;
	    this->p->sharedProducerTaskKeys.erase(taskId);
	    this->p->taskGenerations.erase(taskId);
	    this->p->completed.erase(taskId);
	    producer = this->p->sharedProducers.erase(producer);
	}

	for (std::deque<BObolLodWorkItem>::iterator it =
		 this->p->pending.begin(); it != this->p->pending.end();) {
	    BObolSharedProducer *shared =
		lod_shared_producer_for_task_unlocked(this->p, it->id);
	    const SbBool sharedSurvives =
		shared && !shared->leases.empty() ? TRUE : FALSE;
	    const SbBool sharedUnowned =
		shared && shared->leases.empty() ? TRUE : FALSE;
	    if ((it->task.generation != generation && !sharedUnowned) ||
		(it->task.generation == generation && sharedSurvives)) {
		++it;
		continue;
	    }
	    const uint64_t taskId = it->id;
	    const uint64_t taskGeneration = it->task.generation;
	    this->p->completed.insert(taskId);
	    if (this->p->inFlight > 0)
		this->p->inFlight--;
	    lod_generation_count_remove_unlocked(
		this->p->generationPendingTaskCounts, taskGeneration);
	    lod_generation_pending_time_remove_unlocked(
		this->p, taskGeneration, it->submittedMicroseconds);
	    lod_generation_task_finished_unlocked(
		this->p, taskGeneration);
	    if (it->task.publishResult && this->p->resultReservations > 0)
		this->p->resultReservations--;
	    if (this->p->cacheWriterEnabled && it->task.writeCache &&
		it->task.cacheWrite && this->p->cacheWriteReservations > 0)
		this->p->cacheWriteReservations--;
	    lod_active_request_key_remove_unlocked(this->p,
		lod_request_active_key(it->task.request));
	    lod_pending_quality_remove(
		this->p, it->task.request.qualityTier);
	    lod_pending_dispatch_remove(
		this->p, it->task.dispatchClass);
	    cancelled.push_back(std::move(*it));
	    it = this->p->pending.erase(it);
	    const auto producerKey =
		this->p->sharedProducerTaskKeys.find(taskId);
	    if (producerKey != this->p->sharedProducerTaskKeys.end()) {
		this->p->sharedProducers.erase(producerKey->second);
		this->p->sharedProducerTaskKeys.erase(producerKey);
	    }
	    this->p->taskGenerations.erase(taskId);
	    this->p->completed.erase(taskId);
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
	    if (it->result.generation == generation &&
		!lod_shared_producer_has_leases_unlocked(
		    this->p, it->result.request)) {
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
	    BObolSharedProducer *shared =
		lod_shared_producer_for_task_unlocked(this->p, it->first);
	    if (shared && !shared->leases.empty()) {
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
    return this->p->stopping ||
	lod_generation_cancelled_unlocked(this->p, generation);
}

SbBool
BObolLodService::isProducerCancelled(
    uint64_t generation, const BObolLodRequest &request) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->stopping || lod_producer_cancelled_unlocked(
	this->p, generation, request) ? TRUE : FALSE;
}

void
BObolLodService::setQueueLimits(size_t maxActiveTasks,
	size_t maxQueuedResults, size_t maxQueuedCacheWrites)
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    this->p->maxActiveTasks = maxActiveTasks > 0 ? maxActiveTasks : 1;
    this->p->maxQueuedResults = maxQueuedResults > 0 ? maxQueuedResults : 1;
    this->p->maxDeferredIntermediateResults.store(
	std::max<size_t>(1, std::min(this->p->maxQueuedResults,
	    this->p->workers.empty() ? static_cast<size_t>(1) :
		this->p->workers.size())), std::memory_order_relaxed);
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
	std::lock_guard<std::mutex> admissionLock(
	    this->p->residentMeshAdmissionMutex);
	const size_t prior = this->p->maxResidentMeshBytes;
	const size_t next = maxResidentBytes > 0 ?
	    maxResidentBytes : lod_default_resident_mesh_limit();
	this->p->maxResidentMeshBytes = next;
	this->p->residentMeshLimitPercent = 0.0;
	this->p->residentMeshLimitBasisBytes = 0;
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

SbBool
BObolLodService::setResidentMeshAvailableMemoryPercent(
	double availableMemoryPercent)
{
    if (!std::isfinite(availableMemoryPercent) ||
	availableMemoryPercent <= 0.0 ||
	availableMemoryPercent > LOD_MAX_RESIDENT_AVAILABLE_MEMORY_PERCENT)
	return FALSE;

    size_t availableBytes = 0;
    if (bu_mem(BU_MEM_AVAIL, &availableBytes) < 0 || availableBytes == 0)
	return FALSE;

    const long double scaled =
	static_cast<long double>(availableBytes) *
	static_cast<long double>(availableMemoryPercent) / 100.0L;
    const size_t next = scaled >= static_cast<long double>(SIZE_MAX) ?
	SIZE_MAX : std::max<size_t>(1, static_cast<size_t>(scaled));
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	std::lock_guard<std::mutex> admissionLock(
	    this->p->residentMeshAdmissionMutex);
	const size_t prior = this->p->maxResidentMeshBytes;
	this->p->maxResidentMeshBytes = next;
	this->p->residentMeshLimitPercent = availableMemoryPercent;
	this->p->residentMeshLimitBasisBytes = availableBytes;
	if (next > prior)
	    lod_resident_mesh_revision_advance(
		this->p->residentMeshAdmissionRevision);
    }
    this->p->workerCv.notify_all();
    return TRUE;
}

double
BObolLodService::getResidentMeshAvailableMemoryPercent(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshLimitPercent;
}

size_t
BObolLodService::getResidentMeshAvailableMemoryBasisBytes(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentMeshLimitBasisBytes;
}

double
BObolLodService::getMaximumResidentMeshAvailableMemoryPercent(void)
{
    return LOD_MAX_RESIDENT_AVAILABLE_MEMORY_PERCENT;
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
				 SbBool skipActiveDuplicate,
				 SbBool *replayReady)
{
    if (!p->running || p->stopping)
	return 0;

    BObolLodWorkItem item;
    item.id = bobol_nonzero_identity_take(p->nextTaskId);
    item.task = task;
    item.submittedMicroseconds = std::max<int64_t>(1, bu_gettime());
    item.retargetActiveDemand = skipActiveDuplicate ? true : false;
    if (item.task.generation == 0) {
	if (p->activeGeneration == 0) {
	    bobol_identity_advance(p->nextGeneration);
	    p->activeGeneration = p->nextGeneration;
	}
	item.task.generation = p->activeGeneration;
    }


    if (skipActiveDuplicate &&
	lod_request_key_recorded_unlocked(p, activeKey)) {
	(void)lod_update_active_request_demand_unlocked(
	    p, activeKey, task.request, item.task.generation, replayReady);
	return 0;
    }

    if (lod_generation_cancelled_unlocked(
	    p, item.task.generation) ||
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

    const uint64_t id = item.id;
    const uint64_t generation = item.task.generation;
    const SbBool publishResult = item.task.publishResult;
    const SbBool reserveCacheWrite =
	p->cacheWriterEnabled && item.task.writeCache && item.task.cacheWrite;
    const int qualityTier = item.task.request.qualityTier;
    const int dispatchClass = item.task.dispatchClass;
    const int64_t submittedMicroseconds = item.submittedMicroseconds;
    p->pending.push_back(std::move(item));
    p->pendingDispatchCounts[dispatchClass]++;
    p->pendingQualityCounts[qualityTier]++;
    p->activeRequestKeyCounts[activeKey.getString()]++;
    if (skipActiveDuplicate && publishResult &&
	task.request.coalesceAssetProducer &&
	activeKey.getLength() > 0) {
	BObolSharedProducer producer;
	producer.taskId = id;
	producer.producerGeneration = generation;
	producer.submittedMicroseconds = submittedMicroseconds;
	BObolSharedProducerLease lease;
	lease.demand = task.request;
	producer.leases.emplace(generation, std::move(lease));
	p->sharedProducers[activeKey.getString()] = std::move(producer);
	p->sharedProducerTaskKeys[id] = activeKey.getString();
    }
    (void)lod_update_active_request_demand_unlocked(
	p, activeKey, task.request, generation);
    p->generationTaskCounts[generation]++;
    p->generationPendingTaskCounts[generation]++;
    lod_generation_pending_time_add_unlocked(
	p, generation, submittedMicroseconds);
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
    SbBool replayReady = FALSE;
    {
	std::lock_guard<std::mutex> lock(p->mutex);
	id = lod_service_submit_task_unlocked(
	    p, task, activeKey, skipActiveDuplicate, &replayReady);
    }

    if (id)
	p->workerCv.notify_all();
    if (replayReady)
	lod_notify_result_ready(p);
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
BObolLodService::tryPublishIntermediateResult(
    uint64_t generation, BObolLodResult &&result)
{
    if (!generation || !lod_defer_intermediate_result(
	    this->p, generation, std::move(result)))
	return FALSE;

    /* Opportunistically publish without waiting for the outer lock.  If it is
     * busy, the mailbox remains authoritative and an idle worker or the next
     * presentation drain will promote it. */
    (void)lod_promote_deferred_intermediate_results(
	this->p, false, false);
    this->p->workerCv.notify_all();
    return TRUE;
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
    SbBool replayReady = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	this->p->activeRequestKeyCounts.reserve(
	    this->p->activeRequestKeyCounts.size() + tasks.size());
	this->p->taskGenerations.reserve(
	    this->p->taskGenerations.size() + tasks.size());
	for (size_t i = 0; i < tasks.size(); ++i) {
	    taskIds[i] = lod_service_submit_task_unlocked(
		this->p, tasks[i], activeKeys[i],
		skipActiveDuplicates ? TRUE : FALSE, &replayReady);
	    if (taskIds[i])
		++accepted;
	}
    }

    if (accepted)
	this->p->workerCv.notify_all();
    if (replayReady)
	lod_notify_result_ready(this->p);
    return accepted;
}

SbBool
BObolLodService::updateActiveRequestDemand(
    const BObolLodRequest &request, uint64_t generation)
{
    const SbString key = lod_request_active_key(request);

    SbBool replayReady = FALSE;
    SbBool recorded = FALSE;
    {
	std::lock_guard<std::mutex> lock(this->p->mutex);
	recorded = lod_request_key_recorded_unlocked(this->p, key);
	if (recorded)
	    recorded = lod_update_active_request_demand_unlocked(
		this->p, key, request, generation, &replayReady);
    }
    if (replayReady)
	lod_notify_result_ready(this->p);
    return recorded;
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
    if (!result.progressiveMesh && result.preparedCadGeometry)
	bytes = saturatingAdd(bytes,
	    bobol_database_part_geometry_estimate_bytes(
		*result.preparedCadGeometry));
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

static void
lod_shared_owner_result_drained_unlocked(BObolLodServicePrivate *p,
	const BObolLodResult &result)
{
    if (!p || !result.request.coalesceAssetProducer)
	return;
    const SbString key = lod_request_active_key(result.request);
    auto producer = p->sharedProducers.find(key.getString());
    if (producer == p->sharedProducers.end() ||
	producer->second.state != BObolSharedProducerState::RESULT ||
	result.generation != producer->second.producerGeneration)
	return;
    const auto lease = producer->second.leases.find(result.generation);
    if (lease != producer->second.leases.end() &&
	lease->second.deliveredPublication >= producer->second.publication)
	producer->second.leases.erase(lease);
    lod_shared_producer_retire_if_unowned_unlocked(p, producer);
}

enum class BObolSharedReplayFilter : uint8_t {
    ANY = 0,
    GENERATION,
    MATCHING
};

static size_t
lod_drain_shared_replays_unlocked(BObolLodServicePrivate *p,
	std::vector<BObolLodResult> &results, BObolSharedReplayFilter filter,
	uint64_t generation, const std::vector<BObolLodRequest> *requests,
	size_t maxResults)
{
    if (!p || !maxResults)
	return 0;

    size_t count = 0;
    for (auto producer = p->sharedProducers.begin();
	 producer != p->sharedProducers.end() && count < maxResults;) {
	BObolSharedProducer &state = producer->second;
	for (auto lease = state.leases.begin();
	     lease != state.leases.end() && count < maxResults;) {
	    if (lease->second.deliveredPublication >= state.publication ||
		(filter == BObolSharedReplayFilter::GENERATION &&
		 lease->first != generation)) {
		++lease;
		continue;
	    }

	    BObolLodResult replay = lod_shared_producer_replay_result(
		lease->first, lease->second.demand);
	    if (filter == BObolSharedReplayFilter::MATCHING) {
		SbBool matched = FALSE;
		if (requests) {
		    for (const BObolLodRequest &request : *requests) {
			if (bobol_lod_result_matches_request(replay, request)) {
			    matched = TRUE;
			    break;
			}
		    }
		}
		if (!matched) {
		    ++lease;
		    continue;
		}
	    }

	    results.push_back(std::move(replay));
	    lease->second.deliveredPublication = state.publication;
	    ++count;
	    if (state.state == BObolSharedProducerState::RESULT)
		lease = state.leases.erase(lease);
	    else
		++lease;
	}

	if (state.state == BObolSharedProducerState::RESULT &&
	    state.leases.empty()) {
	    producer = lod_shared_producer_erase_unlocked(p, producer);
	} else {
	    ++producer;
	}
    }
    return count;
}

size_t
BObolLodService::drainResults(std::vector<BObolLodResult> &results,
				size_t maxResults,
				size_t maxEstimatedBytes)
{
    this->p->deferredIntermediateNotificationPending.store(
	false, std::memory_order_release);
    (void)lod_promote_deferred_intermediate_results(this->p, true, false);
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
	lod_shared_owner_result_drained_unlocked(
	    this->p, this->p->results.front());
	results.push_back(std::move(this->p->results.front()));
	this->p->results.pop_front();
	estimatedBytes = frontBytes > SIZE_MAX - estimatedBytes ?
	    SIZE_MAX : estimatedBytes + frontBytes;
	count++;
    }

    const size_t replayLimit = maxResults ? maxResults - count : SIZE_MAX;
    count += lod_drain_shared_replays_unlocked(
	this->p, results, BObolSharedReplayFilter::ANY, 0, NULL,
	replayLimit);
    if (count)
	this->p->workerCv.notify_all();

    return count;
}

size_t
BObolLodService::drainGenerationResults(
    std::vector<BObolLodResult> &results, uint64_t generation,
    size_t maxResults, size_t maxEstimatedBytes)
{
    if (generation == 0)
	return 0;

    this->p->deferredIntermediateNotificationPending.store(
	false, std::memory_order_release);
    (void)lod_promote_deferred_intermediate_results(this->p, true, false);
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
	lod_shared_owner_result_drained_unlocked(this->p, *it);
	results.push_back(std::move(*it));
	it = this->p->results.erase(it);
	estimatedBytes = resultBytes > SIZE_MAX - estimatedBytes ?
	    SIZE_MAX : estimatedBytes + resultBytes;
	count++;
    }

    const size_t replayLimit = maxResults ? maxResults - count : SIZE_MAX;
    count += lod_drain_shared_replays_unlocked(
	this->p, results, BObolSharedReplayFilter::GENERATION, generation,
	NULL, replayLimit);
    if (count)
	this->p->workerCv.notify_all();

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

    this->p->deferredIntermediateNotificationPending.store(
	false, std::memory_order_release);
    (void)lod_promote_deferred_intermediate_results(this->p, true, false);
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
	lod_shared_owner_result_drained_unlocked(this->p, *it);
	results.push_back(std::move(*it));
	it = this->p->results.erase(it);
	count++;
    }

    const size_t replayLimit = maxResults ? maxResults - count : SIZE_MAX;
    count += lod_drain_shared_replays_unlocked(
	this->p, results, BObolSharedReplayFilter::MATCHING, 0, &requests,
	replayLimit);
    if (count)
	this->p->workerCv.notify_all();

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
    subscriber.id = bobol_nonzero_identity_take(
	this->p->nextSubscriberId);
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

static size_t
lod_service_queued_result_count_unlocked(const BObolLodServicePrivate *service)
{
    size_t count = service->results.size();
    for (const auto &producer : service->sharedProducers) {
	const size_t pending = lod_shared_pending_replay_count_unlocked(
	    producer.second);
	count = pending > SIZE_MAX - count ? SIZE_MAX : count + pending;
    }
    const size_t deferred = service->deferredIntermediateResultCount.load(
	std::memory_order_acquire);
    count = deferred > SIZE_MAX - count ? SIZE_MAX : count + deferred;
    return count;
}

static size_t
lod_service_active_request_count_unlocked(
    const BObolLodServicePrivate *service)
{
    size_t count = service->activeRequestKeyCounts.size();
    for (const auto &producer : service->sharedProducers) {
	if (!producer.second.leases.empty() &&
	    service->activeRequestKeyCounts.find(producer.first) ==
		service->activeRequestKeyCounts.end())
	    ++count;
    }
    return count;
}

static size_t
lod_service_pending_compaction_count_unlocked(
    const BObolLodServicePrivate *service)
{
    size_t planning = 0;
    for (const auto &consumer : service->residentMeshConsumerDemands)
	if (consumer.second.planning)
	    ++planning;
    return service->residentMeshCompactionWork.size() +
	service->residentMeshCompactionsInFlight + planning;
}

static size_t
lod_service_generation_queued_result_count_unlocked(
    const BObolLodServicePrivate *service, uint64_t generation)
{
    if (!generation)
	return 0;
    size_t count = lod_generation_count_unlocked(
	service->generationResultCounts, generation);
    for (const auto &producer : service->sharedProducers) {
	const size_t pending = lod_shared_pending_replay_count_unlocked(
	    producer.second, generation);
	count = pending > SIZE_MAX - count ? SIZE_MAX : count + pending;
    }
    {
	/* Callers hold the service mutex, preserving service-then-mailbox order. */
	std::lock_guard<std::mutex> lock(
	    service->deferredIntermediateMutex);
	for (const auto &entry : service->deferredIntermediateResults) {
	    if (entry.first.generation != generation)
		continue;
	    count = count == SIZE_MAX ? SIZE_MAX : count + 1;
	}
    }
    return count;
}

static size_t
lod_service_generation_lease_count_unlocked(
    const BObolLodServicePrivate *service, uint64_t generation)
{
    if (!generation)
	return 0;
    size_t count = 0;
    for (const auto &producer : service->sharedProducers)
	if (producer.second.leases.find(generation) !=
	    producer.second.leases.end())
	    ++count;
    return count;
}

size_t
BObolLodService::inFlightCount(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->inFlight;
}

BObolLodServiceWorkStatus
BObolLodService::workStatus(void) const
{
    BObolLodServiceWorkStatus status;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    status.running = this->p->running;
    status.stopping = this->p->stopping;
    status.pendingTasks = this->p->pending.size();
    status.executingTasks = this->p->executingTasks;
    status.cpuAdmissionWaitingTasks = this->p->cpuAdmissionWaitingTasks;
    status.transientMemoryAdmissionWaitingTasks =
	this->p->transientMemoryAdmissionWaitingTasks;
    status.inFlightTasks = this->p->inFlight;
    status.resultReservations = this->p->resultReservations;
    status.cacheWriteReservations = this->p->cacheWriteReservations;
    status.activeRequests =
	lod_service_active_request_count_unlocked(this->p);
    status.queuedResults =
	lod_service_queued_result_count_unlocked(this->p);
    status.queuedCacheWrites =
	this->p->cacheWrites.size() + this->p->cacheWriteInFlight;
    status.delayedTasks = this->p->delayedTasks;
    status.activeWorkingSetBytes = this->p->activeWorkingSetBytes;
    status.pendingResidentMeshCompactions =
	lod_service_pending_compaction_count_unlocked(this->p);
    status.queuedResidentMeshCompactionResults =
	this->p->residentMeshCompactionResultCount;
    status.residentMeshCompactionResultReservations =
	this->p->residentMeshCompactionResultReservations;
    return status;
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
    return lod_service_queued_result_count_unlocked(this->p);
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

BObolLodGenerationWorkStatus
BObolLodService::generationWorkStatus(uint64_t generation) const
{
    BObolLodGenerationWorkStatus status;
    if (!generation)
	return status;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    status.activeTasks = lod_generation_count_unlocked(
	this->p->generationTaskCounts, generation);
    status.pendingTasks = lod_generation_count_unlocked(
	this->p->generationPendingTaskCounts, generation);
    status.executingTasks = lod_generation_count_unlocked(
	this->p->generationExecutingTaskCounts, generation);
    status.cpuAdmissionWaitingTasks = lod_generation_count_unlocked(
	this->p->generationCpuAdmissionWaitingTaskCounts, generation);
    status.transientMemoryAdmissionWaitingTasks =
	lod_generation_count_unlocked(
	    this->p->generationTransientMemoryAdmissionWaitingTaskCounts,
	    generation);
    status.delayedTasks = lod_generation_count_unlocked(
	this->p->generationDelayedTaskCounts, generation);
    status.queuedResults =
	lod_service_generation_queued_result_count_unlocked(
	    this->p, generation);
    status.queuedCacheWrites = lod_generation_count_unlocked(
	this->p->generationCacheWriteCounts, generation);
    status.sharedProducerLeases =
	lod_service_generation_lease_count_unlocked(this->p, generation);
    const int64_t observationMicroseconds = bu_gettime();
    int64_t oldestPendingMicroseconds = 0;
    const auto pendingTimes =
	this->p->generationPendingTaskTimes.find(generation);
    if (pendingTimes != this->p->generationPendingTaskTimes.end() &&
	!pendingTimes->second.empty())
	oldestPendingMicroseconds = *pendingTimes->second.begin();
    /* A shared producer may have been submitted by an older generation.  Its
     * current consumer still needs to see that queue delay even though the
     * physical task belongs to the original owner. */
    for (const auto &producer : this->p->sharedProducers) {
	if (producer.second.state != BObolSharedProducerState::BUILDING ||
	    producer.second.leases.find(generation) ==
		producer.second.leases.end())
	    continue;
	if (producer.second.executionStartedMicroseconds == 0) {
	    const int64_t submitted = producer.second.submittedMicroseconds;
	    if (submitted > 0 && (!oldestPendingMicroseconds ||
		    submitted < oldestPendingMicroseconds))
		oldestPendingMicroseconds = submitted;
	    if (producer.second.cpuAdmissionWaiting &&
		producer.second.producerGeneration != generation &&
		status.cpuAdmissionWaitingTasks != SIZE_MAX)
		status.cpuAdmissionWaitingTasks++;
	}
	if (producer.second.transientMemoryAdmissionWaiting &&
	    producer.second.producerGeneration != generation &&
	    status.transientMemoryAdmissionWaitingTasks != SIZE_MAX)
	    status.transientMemoryAdmissionWaitingTasks++;
    }
    if (oldestPendingMicroseconds > 0 &&
	observationMicroseconds > oldestPendingMicroseconds)
	status.oldestPendingTaskAgeMicroseconds =
	    static_cast<uint64_t>(observationMicroseconds -
		oldestPendingMicroseconds);
    {
	std::lock_guard<std::mutex> progressLock(
	    this->p->producerProgressMutex);
	bool primaryDeterminate = false;
	for (const auto &entry : this->p->producerProgress) {
	    const std::shared_ptr<BObolLodProducerProgressRecord> &record =
		entry.second;
	    if (!record)
		continue;
	    const BObolLodProducerProgressRecord::Snapshot progress =
		record->snapshot(observationMicroseconds);
	    if (!progress.active)
		continue;
	    bool relevant = record->ownerGeneration == generation;
	    if (!relevant && !record->activeKey.empty()) {
		const auto producer = this->p->sharedProducers.find(
		    record->activeKey);
		relevant = producer != this->p->sharedProducers.end() &&
		    producer->second.leases.find(generation) !=
			producer->second.leases.end();
	    }
	    if (!relevant)
		continue;

	    const int stage = progress.stage;
	    if (stage <= BOBOL_LOD_PRODUCER_STAGE_NONE ||
		stage >= BOBOL_LOD_PRODUCER_STAGE_COUNT)
		continue;
	    status.activeProducerCount++;
	    status.maximumProducerQueueWaitMicroseconds = std::max(
		status.maximumProducerQueueWaitMicroseconds,
		progress.queueWaitMicroseconds);
	    status.maximumProducerElapsedMicroseconds = std::max(
		status.maximumProducerElapsedMicroseconds,
		progress.elapsedMicroseconds);
	    const auto add = [](uint64_t left, uint64_t right) {
		return right > UINT64_MAX - left ? UINT64_MAX : left + right;
	    };
	    status.activeProducerSourceFaceCount = add(
		status.activeProducerSourceFaceCount,
		progress.sourceFaceCount);
	    status.activeProducerSourcePointCount = add(
		status.activeProducerSourcePointCount,
		progress.sourcePointCount);
	    status.activeProducerSourceByteCount = add(
		status.activeProducerSourceByteCount,
		progress.sourceByteCount);
	    status.producerStageMask |= 1u << (stage - 1);
	    status.producerStageTaskCounts[stage]++;
	    const uint64_t total = progress.totalUnits;
	    const uint64_t completed = progress.completedUnits;
	    if (status.producerStage == BOBOL_LOD_PRODUCER_STAGE_NONE ||
		lod_producer_stage_display_priority(stage) <
		    lod_producer_stage_display_priority(status.producerStage)) {
		status.producerStage = stage;
		status.producerStageTaskCount = 1;
		status.producerStageCompletedUnits = completed;
		status.producerStageTotalUnits = total;
		status.producerStageElapsedMicroseconds =
		    progress.stageElapsedMicroseconds;
		primaryDeterminate = total != 0;
		continue;
	    }
	    if (stage != status.producerStage)
		continue;
	    status.producerStageTaskCount++;
	    status.producerStageElapsedMicroseconds = std::max(
		status.producerStageElapsedMicroseconds,
		progress.stageElapsedMicroseconds);
	    status.producerStageCompletedUnits = completed >
		    UINT64_MAX - status.producerStageCompletedUnits ? UINT64_MAX :
		    status.producerStageCompletedUnits + completed;
	    if (!primaryDeterminate || !total) {
		status.producerStageTotalUnits = 0;
		primaryDeterminate = false;
	    } else {
		status.producerStageTotalUnits = total >
			UINT64_MAX - status.producerStageTotalUnits ? UINT64_MAX :
		    status.producerStageTotalUnits + total;
	    }
	}
    }
    return status;
}

BObolLodQueueStatus
BObolLodService::generationQueueStatusForDiagnostics(
    uint64_t generation) const
{
    BObolLodQueueStatus status;
    if (!generation)
	return status;

    std::lock_guard<std::mutex> lock(this->p->mutex);
    for (const BObolLodWorkItem &item : this->p->pending) {
	bool relevant = item.task.generation == generation;
	if (!relevant) {
	    const BObolSharedProducer *shared =
		lod_shared_producer_for_task_unlocked(this->p, item.id);
	    relevant = shared && shared->leases.find(generation) !=
		shared->leases.end();
	}
	if (!relevant)
	    continue;
	if (!lod_task_dependencies_ready(this->p, item.task)) {
	    status.dependencyBlockedTasks++;
	} else if (!lod_task_working_set_available(this->p, item.task)) {
	    status.transientMemoryBlockedTasks++;
	} else {
	    status.runnableTasks++;
	}
    }

    status.taskCapacityBlocked =
	this->p->inFlight >= this->p->maxActiveTasks ? TRUE : FALSE;
    const size_t reservedResults =
	this->p->resultReservations > SIZE_MAX - this->p->results.size() ?
	    SIZE_MAX : this->p->results.size() + this->p->resultReservations;
    status.resultCapacityBlocked =
	reservedResults >= this->p->maxQueuedResults ? TRUE : FALSE;
    return status;
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
    return lod_service_generation_queued_result_count_unlocked(
	this->p, generation);
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
    return lod_service_active_request_count_unlocked(this->p);
}

size_t
BObolLodService::sharedProducerCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->sharedProducers.size();
}

size_t
BObolLodService::sharedProducerLeaseCountForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    size_t count = 0;
    for (const auto &producer : this->p->sharedProducers) {
	const size_t leases = producer.second.leases.size();
	count = leases > SIZE_MAX - count ? SIZE_MAX : count + leases;
    }
    return count;
}

size_t
BObolLodService::sharedProducerLeaseCountForGeneration(
    uint64_t generation) const
{
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return lod_service_generation_lease_count_unlocked(
	this->p, generation);
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

BObolLodResidentCapacityStatus
BObolLodService::residentCapacityStatus(void) const
{
    BObolLodResidentCapacityStatus status;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    std::lock_guard<std::mutex> admissionLock(
	this->p->residentMeshAdmissionMutex);
    status.stableResidentBytes =
	this->p->residentMeshStableBytes.load(std::memory_order_relaxed);
    status.reservedGrowthBytes =
	this->p->residentMeshGrowthReservationBytes;
    status.residentLimitBytes = this->p->maxResidentMeshBytes;
    return status;
}

size_t
BObolLodService::stableResidentMeshBytesForDiagnostics(void) const
{
    return lod_resident_stable_bytes(this->p);
}

size_t
BObolLodService::reservedResidentMeshGrowthBytesForDiagnostics(void) const
{
    std::lock_guard<std::mutex> lock(
	this->p->residentMeshAdmissionMutex);
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
    return this->p->residentMeshCacheLoads.load(
	std::memory_order_relaxed);
}

uint64_t
BObolLodService::residentMeshHitCountForDiagnostics(void) const
{
    return this->p->residentMeshHits.load(
	std::memory_order_relaxed);
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
    return lod_service_pending_compaction_count_unlocked(this->p);
}

SbBool
BObolLodService::residentMeshCompactionPlanForDiagnostics(
    uint64_t consumerId, uint64_t *demandRevision,
    size_t *candidateCount) const
{
    if (demandRevision)
	*demandRevision = 0;
    if (candidateCount)
	*candidateCount = 0;
    if (!consumerId)
	return FALSE;

    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (!this->p->residentMeshCompactionWork.empty() ||
	this->p->residentMeshCompactionsInFlight != 0)
	return FALSE;
    const auto queuedResults =
	this->p->residentMeshCompactionResults.find(consumerId);
    if (queuedResults != this->p->residentMeshCompactionResults.end() &&
	!queuedResults->second.empty())
	return FALSE;
    const uint64_t residentRevision =
	this->p->residentMeshRevision.load(std::memory_order_relaxed);
    const auto consumer =
	this->p->residentMeshConsumerDemands.find(consumerId);
    if (consumer == this->p->residentMeshConsumerDemands.end() ||
	consumer->second.planning ||
	!lod_resident_consumer_snapshot_current(consumer->second) ||
	consumer->second.completedPlanRevision != consumer->second.revision ||
	!consumer->second.residentRevision ||
	consumer->second.residentRevision != residentRevision)
	return FALSE;
    if (demandRevision)
	*demandRevision = consumer->second.completedPlanRevision;
    if (candidateCount)
	*candidateCount = consumer->second.completedCandidateCount;
    return TRUE;
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
