/*                L O D _ R E A L I Z A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file brlobol/lod_realization.h */

#ifndef BRLOBOL_LOD_REALIZATION_H
#define BRLOBOL_LOD_REALIZATION_H

#include "brlobol/defines.h"
#include "brlobol/mesh_lod_cache.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <stdint.h>
#include <vector>

enum BRLObolLodDrawMode {
    BRLOBOL_LOD_DRAW_UNKNOWN = 0,
    BRLOBOL_LOD_DRAW_WIRE = 1,
    BRLOBOL_LOD_DRAW_SHADED = 2,
    BRLOBOL_LOD_DRAW_POINTS = 3,
    BRLOBOL_LOD_DRAW_DIAGNOSTIC = 4,
    BRLOBOL_LOD_DRAW_SHADED_BOTS = 5,
    BRLOBOL_LOD_DRAW_HIDDEN_LINE = 6
};

enum BRLObolLodQualityTier {
    BRLOBOL_LOD_QUALITY_METADATA = 0,
    BRLOBOL_LOD_QUALITY_ATTRIBUTES = 1,
    BRLOBOL_LOD_QUALITY_PROXY = 2,
    BRLOBOL_LOD_QUALITY_FAST_DISPLAY = 3,
    BRLOBOL_LOD_QUALITY_PROGRESSIVE = 4,
    BRLOBOL_LOD_QUALITY_FULL_DETAIL = 5
};

enum BRLObolLodResultKind {
    BRLOBOL_LOD_RESULT_NONE = 0,
    BRLOBOL_LOD_RESULT_DIRECTORY = 1,
    BRLOBOL_LOD_RESULT_ATTRIBUTES = 2,
    BRLOBOL_LOD_RESULT_AABB = 3,
    BRLOBOL_LOD_RESULT_PROXY = 4,
    BRLOBOL_LOD_RESULT_MESH = 5,
    BRLOBOL_LOD_RESULT_FULL_DETAIL = 6,
    BRLOBOL_LOD_RESULT_DIAGNOSTIC = 7
};

enum BRLObolLodProviderStatus {
    BRLOBOL_LOD_PROVIDER_UNKNOWN = 0,
    BRLOBOL_LOD_PROVIDER_READY = 1,
    BRLOBOL_LOD_PROVIDER_CACHE_MISS = 2,
    BRLOBOL_LOD_PROVIDER_STALE = 3,
    BRLOBOL_LOD_PROVIDER_RUNNING = 4,
    BRLOBOL_LOD_PROVIDER_TERMINAL = 5,
    BRLOBOL_LOD_PROVIDER_FALLBACK = 6,
    BRLOBOL_LOD_PROVIDER_ERROR = 7,
    BRLOBOL_LOD_PROVIDER_CANCELLED = 8
};

enum BRLObolLodGeometryHandleKind {
    BRLOBOL_LOD_GEOMETRY_NONE = 0,
    BRLOBOL_LOD_GEOMETRY_MESH_LOD_CACHE = 1,
    BRLOBOL_LOD_GEOMETRY_OBOL_MESH = 2,
    BRLOBOL_LOD_GEOMETRY_PROVIDER_TOKEN = 3
};

enum BRLObolLodProxyKind {
    BRLOBOL_LOD_PROXY_NONE = 0,
    BRLOBOL_LOD_PROXY_AABB = 1,
    BRLOBOL_LOD_PROXY_OBB = 2,
    BRLOBOL_LOD_PROXY_SPHERE = 3,
    BRLOBOL_LOD_PROXY_PROVIDER_TOKEN = 4
};

struct BRLOBOL_EXPORT BRLObolLodProviderParam {
    SbString name;
    SbString value;
};

struct BRLOBOL_EXPORT BRLObolLodDependency {
    SbString objectPath;
    SbString objectName;
    uint64_t sourceRevision;
    uint64_t sourceContentHash;
    int requiredQualityTier;
    SbBool optional;

    BRLObolLodDependency(void);
    void clear(void);
};

struct BRLOBOL_EXPORT BRLObolLodAttribute {
    SbString name;
    SbString value;

    BRLObolLodAttribute(void);
    void clear(void);
};

struct BRLOBOL_EXPORT BRLObolLodCounts {
    uint64_t faceCount;
    uint64_t pointCount;
    uint64_t originalPointCount;
    uint64_t normalCount;
    uint64_t lineCount;
    uint64_t byteCount;

    BRLObolLodCounts(void);
    void clear(void);
};

struct BRLOBOL_EXPORT BRLObolLodCacheKey {
    SbString value;

    BRLObolLodCacheKey(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BRLOBOL_EXPORT BRLObolLodGeometryHandle {
    int kind;
    SbString providerId;
    SbString providerVersion;
    BRLObolLodCacheKey cacheKey;
    uint64_t providerToken;
    int activeLevel;
    SbBool borrowed;

    BRLObolLodGeometryHandle(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BRLOBOL_EXPORT BRLObolLodMeshPayload {
    std::vector<SbVec3f> points;
    std::vector<SbVec3f> normals;
    std::vector<int32_t> coordIndex;
    std::vector<int32_t> faceIndex;
    std::vector<int32_t> vertexIndex;

    BRLObolLodMeshPayload(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BRLOBOL_EXPORT BRLObolLodProxy {
    int kind;
    SbBox3f bounds;
    SbVec3f center;
    SbVec3f axisX;
    SbVec3f axisY;
    SbVec3f axisZ;
    SbVec3f halfExtents;
    float radius;
    BRLObolLodGeometryHandle geometry;

    BRLObolLodProxy(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BRLOBOL_EXPORT BRLObolLodRequest {
    SbString databaseId;
    uint64_t databaseRevision;
    uint64_t sourceRevision;
    uint64_t sourceContentHash;
    SbString objectPath;
    SbString objectName;
    uint64_t viewRevision;
    uint64_t policyRevision;
    int drawMode;
    uint32_t lodPolicy;
    SbString providerId;
    SbString providerVersion;
    int qualityTier;
    SbBox3f bounds;
    BRLObolLodCounts sourceCounts;
    std::vector<BRLObolLodProviderParam> providerParams;

    BRLObolLodRequest(void);
    void clear(void);
    void addProviderParam(const SbString &name, const SbString &value);
};

struct BRLOBOL_EXPORT BRLObolLodResult {
    BRLObolLodRequest request;
    BRLObolLodCacheKey cacheKey;
    BRLObolLodGeometryHandle geometry;
    BRLObolLodMeshPayload mesh;
    int resultKind;
    int qualityTier;
    int providerStatus;
    SbBox3f bounds;
    BRLObolLodCounts counts;
    std::vector<BRLObolLodDependency> dependencies;
    std::vector<BRLObolLodAttribute> attributes;
    BRLObolLodProxy proxy;
    double estimatedError;
    SbBool terminal;
    SbBool fallback;
    SbBool stale;
    SbBool hasSnappedPoints;
    SbBool hasNormals;
    SbString diagnostic;

    BRLObolLodResult(void);
    void clear(void);
    void addDependency(const SbString &objectPath, const SbString &objectName,
	uint64_t sourceRevision, uint64_t sourceContentHash,
	int requiredQualityTier, SbBool optional = FALSE);
    void addAttribute(const SbString &name, const SbString &value);
};

BRLOBOL_EXPORT BRLObolLodCacheKey
brlobol_lod_cache_key(const BRLObolLodRequest &request);

BRLOBOL_EXPORT SbBool
brlobol_lod_mesh_payload_from_mesh_lod_data(BRLObolLodMeshPayload &payload,
	const struct BRLObolMeshLodData &data);

BRLOBOL_EXPORT SbBool
brlobol_lod_result_matches_request(const BRLObolLodResult &result,
	const BRLObolLodRequest &request);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_lod_result_from_mesh_lod_info(const BRLObolLodRequest &request,
	const struct BRLObolMeshLodInfo &info,
	const struct BRLObolMeshLodCacheStatus *status);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_lod_directory_result(const BRLObolLodRequest &request,
	const std::vector<BRLObolLodDependency> &dependencies);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_lod_attributes_result(const BRLObolLodRequest &request,
	const std::vector<BRLObolLodAttribute> &attributes);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_lod_aabb_result(const BRLObolLodRequest &request,
	const SbBox3f &bounds, const BRLObolLodCounts *counts);

BRLOBOL_EXPORT BRLObolLodResult
brlobol_lod_proxy_result(const BRLObolLodRequest &request,
	const BRLObolLodProxy &proxy, const BRLObolLodCounts *counts);

#endif /* BRLOBOL_LOD_REALIZATION_H */
