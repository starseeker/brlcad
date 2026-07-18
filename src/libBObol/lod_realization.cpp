/*              L O D _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BLodRealization.h"

#include <Inventor/SbVec3f.h>

#include <algorithm>
#include <iomanip>
#include <limits>
#include <sstream>
#include <string>
#include <string.h>

BObolLodDependency::BObolLodDependency(void)
{
    clear();
}

void
BObolLodDependency::clear(void)
{
    objectPath = "";
    objectName = "";
    sourceRevision = 0;
    sourceContentHash = 0;
    requiredQualityTier = BOBOL_LOD_QUALITY_METADATA;
    optional = FALSE;
}

BObolLodAttribute::BObolLodAttribute(void)
{
    clear();
}

void
BObolLodAttribute::clear(void)
{
    name = "";
    value = "";
}

BObolLodCounts::BObolLodCounts(void)
{
    clear();
}

void
BObolLodCounts::clear(void)
{
    faceCount = 0;
    pointCount = 0;
    originalPointCount = 0;
    normalCount = 0;
    lineCount = 0;
    byteCount = 0;
}

BObolLodProxy::BObolLodProxy(void)
{
    clear();
}

void
BObolLodProxy::clear(void)
{
    kind = BOBOL_LOD_PROXY_NONE;
    bounds.makeEmpty();
    center.setValue(0.0f, 0.0f, 0.0f);
    axisX.setValue(1.0f, 0.0f, 0.0f);
    axisY.setValue(0.0f, 1.0f, 0.0f);
    axisZ.setValue(0.0f, 0.0f, 1.0f);
    halfExtents.setValue(0.0f, 0.0f, 0.0f);
    radius = 0.0f;
    geometry.clear();
}

SbBool
BObolLodProxy::isValid(void) const
{
    if (kind == BOBOL_LOD_PROXY_NONE)
	return FALSE;
    if (kind == BOBOL_LOD_PROXY_PROVIDER_TOKEN)
	return geometry.isValid();
    if (kind == BOBOL_LOD_PROXY_SPHERE)
	return radius > 0.0f ? TRUE : FALSE;
    if (kind == BOBOL_LOD_PROXY_AABB)
	return bounds.isEmpty() ? FALSE : TRUE;
    if (kind == BOBOL_LOD_PROXY_OBB)
	return halfExtents[0] >= 0.0f && halfExtents[1] >= 0.0f &&
	       halfExtents[2] >= 0.0f &&
	       (halfExtents[0] > 0.0f || halfExtents[1] > 0.0f ||
		halfExtents[2] > 0.0f) ? TRUE : FALSE;
    return FALSE;
}

BObolLodCacheKey::BObolLodCacheKey(void)
{
    clear();
}

void
BObolLodCacheKey::clear(void)
{
    value = "";
}

SbBool
BObolLodCacheKey::isValid(void) const
{
    return value.getLength() > 0 ? TRUE : FALSE;
}

BObolLodGeometryHandle::BObolLodGeometryHandle(void)
{
    clear();
}

void
BObolLodGeometryHandle::clear(void)
{
    kind = BOBOL_LOD_GEOMETRY_NONE;
    providerId = "";
    providerVersion = "";
    cacheKey.clear();
    providerToken = 0;
    activeLevel = -1;
    borrowed = FALSE;
}

SbBool
BObolLodGeometryHandle::isValid(void) const
{
    return kind != BOBOL_LOD_GEOMETRY_NONE &&
	   (cacheKey.isValid() || providerToken != 0 ||
	    providerId.getLength() > 0) ? TRUE : FALSE;
}

BObolLodMeshPayload::BObolLodMeshPayload(void)
{
    clear();
}

void
BObolLodMeshPayload::clear(void)
{
    points.clear();
    normals.clear();
    coordIndex.clear();
    faceIndex.clear();
    vertexIndex.clear();
}

SbBool
BObolLodMeshPayload::isValid(void) const
{
    const size_t faceCount = coordIndex.size() / 3;
    return !points.empty() && coordIndex.size() >= 3 &&
	   coordIndex.size() % 3 == 0 &&
	   (normals.empty() || normals.size() == coordIndex.size()) &&
	   (faceIndex.empty() || faceIndex.size() == faceCount) &&
	   (vertexIndex.empty() || vertexIndex.size() == points.size()) ?
	   TRUE : FALSE;
}

BObolLodRequest::BObolLodRequest(void)
{
    clear();
}

void
BObolLodRequest::clear(void)
{
    databaseId = "";
    databaseRevision = 0;
    sourceRevision = 0;
    sourceContentHash = 0;
    objectPath = "";
    objectName = "";
    occurrenceKey = "";
    viewRevision = 0;
    policyRevision = 0;
    drawMode = BOBOL_LOD_DRAW_UNKNOWN;
    lodPolicy = 0;
    providerId = "";
    providerVersion = "";
    qualityTier = BOBOL_LOD_QUALITY_METADATA;
    bounds.makeEmpty();
    sourceCounts.clear();
    providerParams.clear();
}

void
BObolLodRequest::addProviderParam(const SbString &name,
				    const SbString &value)
{
    BObolLodProviderParam param;
    param.name = name;
    param.value = value;
    providerParams.push_back(param);
}

BObolLodResult::BObolLodResult(void)
{
    clear();
}

void
BObolLodResult::clear(void)
{
    generation = 0;
    request.clear();
    cacheKey.clear();
    geometry.clear();
    mesh.clear();
    resultKind = BOBOL_LOD_RESULT_NONE;
    qualityTier = BOBOL_LOD_QUALITY_METADATA;
    providerStatus = BOBOL_LOD_PROVIDER_UNKNOWN;
    bounds.makeEmpty();
    counts.clear();
    dependencies.clear();
    attributes.clear();
    proxy.clear();
    estimatedError = -1.0;
    terminal = FALSE;
    fallback = FALSE;
    stale = FALSE;
    hasSnappedPoints = FALSE;
    hasNormals = FALSE;
    diagnostic = "";
}

void
BObolLodResult::addDependency(const SbString &objectPath,
				const SbString &objectName, uint64_t sourceRevision,
				uint64_t sourceContentHash, int requiredQualityTier, SbBool optional)
{
    BObolLodDependency dependency;

    dependency.objectPath = objectPath;
    dependency.objectName = objectName;
    dependency.sourceRevision = sourceRevision;
    dependency.sourceContentHash = sourceContentHash;
    dependency.requiredQualityTier = requiredQualityTier;
    dependency.optional = optional;
    dependencies.push_back(dependency);
}

void
BObolLodResult::addAttribute(const SbString &name, const SbString &value)
{
    BObolLodAttribute attribute;

    attribute.name = name;
    attribute.value = value;
    attributes.push_back(attribute);
}

static void
append_string_field(std::ostringstream &out, const char *name,
		    const SbString &value)
{
    const char *str = value.getString();
    size_t len = str ? strlen(str) : 0;
    out << name << '=' << len << ':' << (str ? str : "") << ';';
}

static void
append_uint_field(std::ostringstream &out, const char *name, uint64_t value)
{
    out << name << '=' << value << ';';
}

static void
append_int_field(std::ostringstream &out, const char *name, int value)
{
    out << name << '=' << value << ';';
}

static void
append_bounds_field(std::ostringstream &out, const char *name,
		    const SbBox3f &bounds)
{
    out << name << '=';
    if (bounds.isEmpty()) {
	out << "empty;";
	return;
    }

    const SbVec3f &bmin = bounds.getMin();
    const SbVec3f &bmax = bounds.getMax();
    out << std::setprecision(9)
	<< bmin[0] << ',' << bmin[1] << ',' << bmin[2] << ','
	<< bmax[0] << ',' << bmax[1] << ',' << bmax[2] << ';';
}

static bool
param_less(const BObolLodProviderParam &a,
	   const BObolLodProviderParam &b)
{
    int name_cmp = bu_strcmp(a.name.getString(), b.name.getString());
    if (name_cmp != 0)
	return name_cmp < 0;
    return bu_strcmp(a.value.getString(), b.value.getString()) < 0;
}

BObolLodCacheKey
bobol_lod_cache_key(const BObolLodRequest &request)
{
    BObolLodCacheKey key;
    std::ostringstream out;

    out << "bobol-lod-v2;";
    append_string_field(out, "database_id", request.databaseId);
    append_uint_field(out, "database_revision", request.databaseRevision);
    append_uint_field(out, "source_revision", request.sourceRevision);
    append_uint_field(out, "source_content_hash", request.sourceContentHash);
    append_string_field(out, "object_path", request.objectPath);
    append_string_field(out, "object_name", request.objectName);
    append_string_field(out, "occurrence_key", request.occurrenceKey);
    append_uint_field(out, "view_revision", request.viewRevision);
    append_uint_field(out, "policy_revision", request.policyRevision);
    append_int_field(out, "draw_mode", request.drawMode);
    append_uint_field(out, "lod_policy", request.lodPolicy);
    append_string_field(out, "provider_id", request.providerId);
    append_string_field(out, "provider_version", request.providerVersion);
    append_int_field(out, "quality_tier", request.qualityTier);
    append_bounds_field(out, "bounds", request.bounds);
    append_uint_field(out, "face_count", request.sourceCounts.faceCount);
    append_uint_field(out, "point_count", request.sourceCounts.pointCount);
    append_uint_field(out, "original_point_count",
		      request.sourceCounts.originalPointCount);
    append_uint_field(out, "normal_count", request.sourceCounts.normalCount);
    append_uint_field(out, "line_count", request.sourceCounts.lineCount);
    append_uint_field(out, "byte_count", request.sourceCounts.byteCount);

    std::vector<BObolLodProviderParam> params = request.providerParams;
    std::sort(params.begin(), params.end(), param_less);
    append_uint_field(out, "provider_param_count",
		      static_cast<uint64_t>(params.size()));
    for (size_t i = 0; i < params.size(); i++) {
	append_string_field(out, "provider_param_name", params[i].name);
	append_string_field(out, "provider_param_value", params[i].value);
    }

    key.value = out.str().c_str();
    return key;
}

SbBool
bobol_lod_mesh_payload_from_mesh_lod_data(BObolLodMeshPayload &payload,
	const struct BObolMeshLodData &data)
{
    payload.clear();

    if (!data.faces || !data.points || data.face_count == 0 ||
	data.point_count == 0)
	return FALSE;

    if (data.point_count >
	static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	data.face_count >
	static_cast<size_t>(std::numeric_limits<int32_t>::max()) / 3)
	return FALSE;

    size_t index_count = data.face_count * 3;
    if ((data.normals && data.normal_count != index_count) ||
	(!data.normals && data.normal_count != 0))
	return FALSE;

    payload.points.reserve(data.point_count);
    for (size_t i = 0; i < data.point_count; i++) {
	payload.points.push_back(SbVec3f(
				     static_cast<float>(data.points[i][X]),
				     static_cast<float>(data.points[i][Y]),
				     static_cast<float>(data.points[i][Z])));
    }

    payload.coordIndex.reserve(index_count);
    for (size_t i = 0; i < index_count; i++) {
	int idx = data.faces[i];
	if (idx < 0 || static_cast<size_t>(idx) >= data.point_count) {
	    payload.clear();
	    return FALSE;
	}
	payload.coordIndex.push_back(static_cast<int32_t>(idx));
    }

    if (data.normals) {
	payload.normals.reserve(index_count);
	for (size_t i = 0; i < index_count; i++) {
	    payload.normals.push_back(SbVec3f(
					  static_cast<float>(data.normals[i][X]),
					  static_cast<float>(data.normals[i][Y]),
					  static_cast<float>(data.normals[i][Z])));
	}
    }

    return payload.isValid();
}

SbBool
bobol_lod_result_matches_request(const BObolLodResult &result,
				   const BObolLodRequest &request)
{
    BObolLodCacheKey expected = bobol_lod_cache_key(request);
    BObolLodCacheKey actual = result.cacheKey.isValid() ?
				result.cacheKey : bobol_lod_cache_key(result.request);

    return bu_strcmp(expected.value.getString(), actual.value.getString()) == 0 ?
	   TRUE : FALSE;
}

BObolLodResult
bobol_lod_result_from_mesh_lod_info(
    const BObolLodRequest &request,
    const struct BObolMeshLodInfo &info,
    const struct BObolMeshLodCacheStatus *status)
{
    BObolLodResult result;

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = BOBOL_LOD_RESULT_MESH;
    result.qualityTier = request.qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.bounds = SbBox3f(
			SbVec3f(static_cast<float>(info.bmin[X]),
				static_cast<float>(info.bmin[Y]),
				static_cast<float>(info.bmin[Z])),
			SbVec3f(static_cast<float>(info.bmax[X]),
				static_cast<float>(info.bmax[Y]),
				static_cast<float>(info.bmax[Z])));
    result.counts.faceCount = info.face_count;
    result.counts.pointCount = info.point_count;
    result.counts.originalPointCount = info.point_orig_count;
    result.counts.normalCount = info.normal_count;
    result.hasSnappedPoints = info.has_snapped_points ? TRUE : FALSE;
    result.hasNormals = info.has_normals ? TRUE : FALSE;

    result.geometry.kind = BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = result.cacheKey;
    result.geometry.activeLevel = info.active_level;
    result.geometry.borrowed = FALSE;

    if (status) {
	result.geometry.providerToken = status->cache_key;
	result.stale = status->stale_cache_entry ? TRUE : FALSE;
	if (status->stale_cache_entry) {
	    result.providerStatus = BOBOL_LOD_PROVIDER_STALE;
	    result.diagnostic = "stale Obol mesh LoD cache entry";
	} else if (!status->has_cache_key || !status->has_cached_payload) {
	    result.providerStatus = BOBOL_LOD_PROVIDER_CACHE_MISS;
	    result.diagnostic = "Obol mesh LoD cache payload unavailable";
	}
    }

    if (!info.has_faces || !info.has_points) {
	result.providerStatus = BOBOL_LOD_PROVIDER_CACHE_MISS;
	result.diagnostic = "Obol mesh LoD result has no active mesh payload";
    }

    return result;
}

static BObolLodResult
stage_result_base(const BObolLodRequest &request, int resultKind,
		  int qualityTier)
{
    BObolLodResult result;

    result.request = request;
    result.cacheKey = bobol_lod_cache_key(request);
    result.resultKind = resultKind;
    result.qualityTier = qualityTier;
    result.providerStatus = BOBOL_LOD_PROVIDER_READY;
    result.terminal = FALSE;

    return result;
}

BObolLodResult
bobol_lod_directory_result(const BObolLodRequest &request,
			     const std::vector<BObolLodDependency> &dependencies)
{
    BObolLodResult result = stage_result_base(request,
			      BOBOL_LOD_RESULT_DIRECTORY, BOBOL_LOD_QUALITY_METADATA);

    result.dependencies = dependencies;
    return result;
}

BObolLodResult
bobol_lod_attributes_result(const BObolLodRequest &request,
			      const std::vector<BObolLodAttribute> &attributes)
{
    BObolLodResult result = stage_result_base(request,
			      BOBOL_LOD_RESULT_ATTRIBUTES, BOBOL_LOD_QUALITY_ATTRIBUTES);

    result.attributes = attributes;
    return result;
}

BObolLodResult
bobol_lod_aabb_result(const BObolLodRequest &request,
			const SbBox3f &bounds, const BObolLodCounts *counts)
{
    BObolLodResult result = stage_result_base(request,
			      BOBOL_LOD_RESULT_AABB, BOBOL_LOD_QUALITY_PROXY);

    result.bounds = bounds;
    result.proxy.kind = BOBOL_LOD_PROXY_AABB;
    result.proxy.bounds = bounds;
    if (counts)
	result.counts = *counts;
    return result;
}

BObolLodResult
bobol_lod_proxy_result(const BObolLodRequest &request,
			 const BObolLodProxy &proxy, const BObolLodCounts *counts)
{
    BObolLodResult result = stage_result_base(request,
			      BOBOL_LOD_RESULT_PROXY, BOBOL_LOD_QUALITY_PROXY);

    result.proxy = proxy;
    result.bounds = proxy.bounds;
    result.geometry = proxy.geometry;
    if (counts)
	result.counts = *counts;

    if (!proxy.isValid()) {
	result.providerStatus = BOBOL_LOD_PROVIDER_ERROR;
	result.terminal = TRUE;
	result.diagnostic = "LoD proxy result has no valid proxy payload";
    }

    return result;
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
