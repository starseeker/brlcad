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
#include <array>
#include <iomanip>
#include <limits>
#include <mutex>
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

BObolLodResidentDemand::BObolLodResidentDemand(void) :
    level(-1)
{
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

struct BObolLodProgressiveMeshPrivate {
    mutable std::mutex mutex;
    BObolLodMeshPayload mesh;
    std::array<size_t, BOBOL_MESH_LOD_LEVEL_COUNT> pointCount = {};
    std::array<size_t, BOBOL_MESH_LOD_LEVEL_COUNT> faceCount = {};
    SbBox3f bounds;
    SbVec3f quantizationMinimum;
    SbVec3f quantizationMaximum;
    int minimumLevel = -1;
    int maximumLevel = -1;
    int residentLevel = -1;
    uint64_t revision = 0;
    SbBool shadedCullBackfaces = FALSE;
};

static int
progressive_level_clamp(const BObolLodProgressiveMeshPrivate *p, int level)
{
    if (!p || p->residentLevel < 0)
	return -1;
    if (level < p->minimumLevel)
	level = p->minimumLevel;
    if (level > p->residentLevel)
	level = p->residentLevel;
    return level;
}

BObolLodProgressiveMesh::BObolLodProgressiveMesh(void) :
    p(new BObolLodProgressiveMeshPrivate)
{
}

BObolLodProgressiveMesh::~BObolLodProgressiveMesh(void)
{
    delete this->p;
    this->p = NULL;
}

SbBool
BObolLodProgressiveMesh::update(
    const struct BObolMeshLodData &data,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int residentLevel,
    SbBool shadedCullBackfaces)
{
    if (!this->p || !data.faces || !data.points_orig ||
	data.face_count == 0 || data.point_orig_count == 0 ||
	residentLevel < hierarchy.min_level ||
	residentLevel > hierarchy.max_level ||
	residentLevel >= BOBOL_MESH_LOD_LEVEL_COUNT)
	return FALSE;

    const size_t indexCount = data.face_count * 3;
    if (data.point_orig_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	data.face_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) / 3 ||
	(data.normals && data.normal_count != indexCount))
	return FALSE;
    for (size_t i = 0; i < indexCount; ++i) {
	if (data.faces[i] < 0 ||
	    static_cast<size_t>(data.faces[i]) >= data.point_orig_count)
	    return FALSE;
    }

    std::lock_guard<std::mutex> lock(this->p->mutex);
    const size_t oldPointCount = this->p->mesh.points.size();
    const size_t oldIndexCount = this->p->mesh.coordIndex.size();
    if (data.point_orig_count < oldPointCount || indexCount < oldIndexCount) {
	/* This is an explicit stable-view trim. */
	this->p->mesh.points.resize(data.point_orig_count);
	this->p->mesh.points.shrink_to_fit();
	this->p->mesh.coordIndex.resize(indexCount);
	this->p->mesh.coordIndex.shrink_to_fit();
	if (!this->p->mesh.normals.empty()) {
	    this->p->mesh.normals.resize(indexCount);
	    this->p->mesh.normals.shrink_to_fit();
	}
    } else {
	this->p->mesh.points.reserve(data.point_orig_count);
	for (size_t i = oldPointCount; i < data.point_orig_count; ++i) {
	    this->p->mesh.points.push_back(SbVec3f(
		static_cast<float>(data.points_orig[i][X]),
		static_cast<float>(data.points_orig[i][Y]),
		static_cast<float>(data.points_orig[i][Z])));
	}
	this->p->mesh.coordIndex.reserve(indexCount);
	for (size_t i = oldIndexCount; i < indexCount; ++i) {
	    this->p->mesh.coordIndex.push_back(
		static_cast<int32_t>(data.faces[i]));
	}
	if (data.normals) {
	    if (this->p->mesh.normals.empty() && oldIndexCount != 0) {
		this->p->mesh.normals.reserve(indexCount);
		for (size_t i = 0; i < oldIndexCount; ++i)
		    this->p->mesh.normals.push_back(SbVec3f(
			static_cast<float>(data.normals[i][X]),
			static_cast<float>(data.normals[i][Y]),
			static_cast<float>(data.normals[i][Z])));
	    }
	    this->p->mesh.normals.reserve(indexCount);
	    for (size_t i = this->p->mesh.normals.size();
		 i < indexCount; ++i)
		this->p->mesh.normals.push_back(SbVec3f(
		    static_cast<float>(data.normals[i][X]),
		    static_cast<float>(data.normals[i][Y]),
		    static_cast<float>(data.normals[i][Z])));
	} else {
	    this->p->mesh.normals.clear();
	}
    }

    if (!this->p->mesh.isValid())
	return FALSE;

    for (int level = 0; level < BOBOL_MESH_LOD_LEVEL_COUNT; ++level) {
	this->p->pointCount[level] = hierarchy.point_count[level];
	this->p->faceCount[level] = hierarchy.face_count[level];
    }
    this->p->minimumLevel = hierarchy.min_level;
    this->p->maximumLevel = hierarchy.max_level;
    this->p->residentLevel = residentLevel;
    this->p->bounds = SbBox3f(
	SbVec3f(static_cast<float>(data.bmin[X]),
		static_cast<float>(data.bmin[Y]),
		static_cast<float>(data.bmin[Z])),
	SbVec3f(static_cast<float>(data.bmax[X]),
		static_cast<float>(data.bmax[Y]),
		static_cast<float>(data.bmax[Z])));
    this->p->quantizationMinimum.setValue(
	static_cast<float>(hierarchy.quantization_min[X]),
	static_cast<float>(hierarchy.quantization_min[Y]),
	static_cast<float>(hierarchy.quantization_min[Z]));
    this->p->quantizationMaximum.setValue(
	static_cast<float>(hierarchy.quantization_max[X]),
	static_cast<float>(hierarchy.quantization_max[Y]),
	static_cast<float>(hierarchy.quantization_max[Z]));
    this->p->shadedCullBackfaces = shadedCullBackfaces;
    this->p->revision++;
    if (this->p->revision == 0)
	this->p->revision++;
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::trim(int residentLevel)
{
    if (!this->p)
	return FALSE;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const int level = progressive_level_clamp(this->p, residentLevel);
    if (level < 0)
	return FALSE;
    const size_t points = this->p->pointCount[level];
    const size_t indices = this->p->faceCount[level] * 3;
    if (points > this->p->mesh.points.size() ||
	indices > this->p->mesh.coordIndex.size())
	return FALSE;
    if (level == this->p->residentLevel)
	return TRUE;
    this->p->mesh.points.resize(points);
    this->p->mesh.points.shrink_to_fit();
    this->p->mesh.coordIndex.resize(indices);
    this->p->mesh.coordIndex.shrink_to_fit();
    if (!this->p->mesh.normals.empty()) {
	this->p->mesh.normals.resize(indices);
	this->p->mesh.normals.shrink_to_fit();
    }
    this->p->residentLevel = level;
    this->p->revision++;
    if (this->p->revision == 0)
	this->p->revision++;
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::copyLevel(
    BObolLodMeshPayload &payload, int requestedLevel) const
{
    payload.clear();
    if (!this->p)
	return FALSE;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const int level = progressive_level_clamp(this->p, requestedLevel);
    if (level < 0)
	return FALSE;
    const size_t points = this->p->pointCount[level];
    const size_t indices = this->p->faceCount[level] * 3;
    if (points > this->p->mesh.points.size() ||
	indices > this->p->mesh.coordIndex.size())
	return FALSE;
    payload.points.assign(this->p->mesh.points.begin(),
	this->p->mesh.points.begin() + points);
    payload.coordIndex.assign(this->p->mesh.coordIndex.begin(),
	this->p->mesh.coordIndex.begin() + indices);
    if (!this->p->mesh.normals.empty())
	payload.normals.assign(this->p->mesh.normals.begin(),
	    this->p->mesh.normals.begin() + indices);
    return payload.isValid();
}

SbBool
BObolLodProgressiveMesh::isValid(void) const
{
    if (!this->p)
	return FALSE;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentLevel >= 0 &&
	this->p->mesh.isValid() ? TRUE : FALSE;
}

SbBool
BObolLodProgressiveMesh::canDrawLevel(int requestedLevel) const
{
    if (!this->p)
	return FALSE;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (requestedLevel < this->p->minimumLevel ||
	requestedLevel > this->p->maximumLevel ||
	requestedLevel >= BOBOL_MESH_LOD_LEVEL_COUNT)
	return FALSE;
    return this->p->pointCount[requestedLevel] <=
	    this->p->mesh.points.size() &&
	this->p->faceCount[requestedLevel] * 3 <=
	    this->p->mesh.coordIndex.size() ? TRUE : FALSE;
}

int
BObolLodProgressiveMesh::minimumLevel(void) const
{
    if (!this->p)
	return -1;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->minimumLevel;
}

int
BObolLodProgressiveMesh::maximumLevel(void) const
{
    if (!this->p)
	return -1;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->maximumLevel;
}

int
BObolLodProgressiveMesh::residentLevel(void) const
{
    if (!this->p)
	return -1;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->residentLevel;
}

uint64_t
BObolLodProgressiveMesh::revision(void) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->revision;
}

size_t
BObolLodProgressiveMesh::pointCount(int requestedLevel) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const int level = progressive_level_clamp(this->p, requestedLevel);
    return level >= 0 ? this->p->pointCount[level] : 0;
}

size_t
BObolLodProgressiveMesh::faceCount(int requestedLevel) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    const int level = progressive_level_clamp(this->p, requestedLevel);
    return level >= 0 ? this->p->faceCount[level] : 0;
}

size_t
BObolLodProgressiveMesh::hierarchyPointCount(int requestedLevel) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (requestedLevel < this->p->minimumLevel)
	requestedLevel = this->p->minimumLevel;
    if (requestedLevel > this->p->maximumLevel)
	requestedLevel = this->p->maximumLevel;
    return requestedLevel >= 0 ?
	this->p->pointCount[requestedLevel] : 0;
}

size_t
BObolLodProgressiveMesh::hierarchyFaceCount(int requestedLevel) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    if (requestedLevel < this->p->minimumLevel)
	requestedLevel = this->p->minimumLevel;
    if (requestedLevel > this->p->maximumLevel)
	requestedLevel = this->p->maximumLevel;
    return requestedLevel >= 0 ?
	this->p->faceCount[requestedLevel] : 0;
}

size_t
BObolLodProgressiveMesh::estimateBytes(void) const
{
    if (!this->p)
	return 0;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->mesh.points.capacity() * sizeof(SbVec3f) +
	this->p->mesh.normals.capacity() * sizeof(SbVec3f) +
	this->p->mesh.coordIndex.capacity() * sizeof(int32_t);
}

SbBox3f
BObolLodProgressiveMesh::bounds(void) const
{
    if (!this->p)
	return SbBox3f();
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->bounds;
}

SbVec3f
BObolLodProgressiveMesh::quantizationMinimum(void) const
{
    if (!this->p)
	return SbVec3f(0.0f, 0.0f, 0.0f);
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->quantizationMinimum;
}

SbVec3f
BObolLodProgressiveMesh::quantizationMaximum(void) const
{
    if (!this->p)
	return SbVec3f(0.0f, 0.0f, 0.0f);
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->quantizationMaximum;
}

SbBool
BObolLodProgressiveMesh::cullBackfaces(void) const
{
    if (!this->p)
	return FALSE;
    std::lock_guard<std::mutex> lock(this->p->mutex);
    return this->p->shadedCullBackfaces;
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
    projectedPixelDiameter = 0.0f;
    targetPixelError = 1.0f;
    requestedLevel = -1;
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
    progressiveMesh.reset();
    residentLevel = -1;
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
    shadedCullBackfaces = FALSE;
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
    /* The selected discrete level determines the provider output.  Raw screen
     * measurements do not: including them would turn sub-pixel camera jitter
     * into distinct active/cache keys for the same payload. */
    append_int_field(out, "requested_level", request.requestedLevel);
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

BObolLodCacheKey
bobol_lod_asset_cache_key(const BObolLodRequest &request)
{
    BObolLodCacheKey key;
    std::ostringstream out;

    out << "bobol-progressive-asset-v1;";
    append_string_field(out, "database_id", request.databaseId);
    append_string_field(out, "object_name", request.objectName);
    if (request.sourceContentHash != 0)
	append_uint_field(out, "source_content_hash",
	    request.sourceContentHash);
    else {
	append_uint_field(out, "database_revision",
	    request.databaseRevision);
	append_uint_field(out, "source_revision", request.sourceRevision);
	/* Database object names are unique.  The occurrence path is a consumer
	 * identity and would duplicate one source asset for every comb leaf. */
	if (request.objectName.getLength() == 0)
	    append_string_field(out, "object_path", request.objectPath);
    }
    append_string_field(out, "provider_id", request.providerId);
    append_string_field(out, "provider_version", request.providerVersion);
    append_int_field(out, "quality_tier", request.qualityTier);

    std::vector<BObolLodProviderParam> params = request.providerParams;
    std::sort(params.begin(), params.end(), param_less);
    append_uint_field(out, "provider_param_count",
		      static_cast<uint64_t>(params.size()));
    for (const BObolLodProviderParam &param : params) {
	append_string_field(out, "provider_param_name", param.name);
	append_string_field(out, "provider_param_value", param.value);
    }

    key.value = out.str().c_str();
    return key;
}

BObolLodCacheKey
bobol_lod_geometry_cache_key(const BObolLodRequest &request)
{
    return bobol_lod_asset_cache_key(request);
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
    result.shadedCullBackfaces =
	info.shaded_cull_backfaces ? TRUE : FALSE;

    result.geometry.kind = BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
    result.geometry.providerId = request.providerId;
    result.geometry.providerVersion = request.providerVersion;
    result.geometry.cacheKey = bobol_lod_geometry_cache_key(request);
    result.geometry.activeLevel = info.active_level;
    result.residentLevel = info.active_level;
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
