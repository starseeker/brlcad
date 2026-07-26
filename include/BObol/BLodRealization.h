/*             B L O D R E A L I Z A T I O N . H
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */
/** @file BObol/BLodRealization.h */

#ifndef BOBOL_BLODREALIZATION_H
#define BOBOL_BLODREALIZATION_H

#include "BObol/BDefines.h"
#include "BObol/BMeshLodCache.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <memory>
#include <stdint.h>
#include <vector>

enum BObolLodDrawMode {
    BOBOL_LOD_DRAW_UNKNOWN = 0,
    BOBOL_LOD_DRAW_WIRE = 1,
    BOBOL_LOD_DRAW_SHADED = 2,
    BOBOL_LOD_DRAW_POINTS = 3,
    BOBOL_LOD_DRAW_DIAGNOSTIC = 4,
    BOBOL_LOD_DRAW_SHADED_BOTS = 5,
    BOBOL_LOD_DRAW_HIDDEN_LINE = 6
};

enum BObolLodQualityTier {
    BOBOL_LOD_QUALITY_METADATA = 0,
    BOBOL_LOD_QUALITY_ATTRIBUTES = 1,
    BOBOL_LOD_QUALITY_PROXY = 2,
    BOBOL_LOD_QUALITY_FAST_DISPLAY = 3,
    BOBOL_LOD_QUALITY_PROGRESSIVE = 4,
    BOBOL_LOD_QUALITY_FULL_DETAIL = 5
};

enum BObolLodResultKind {
    BOBOL_LOD_RESULT_NONE = 0,
    BOBOL_LOD_RESULT_DIRECTORY = 1,
    BOBOL_LOD_RESULT_ATTRIBUTES = 2,
    BOBOL_LOD_RESULT_AABB = 3,
    BOBOL_LOD_RESULT_PROXY = 4,
    BOBOL_LOD_RESULT_MESH = 5,
    BOBOL_LOD_RESULT_FULL_DETAIL = 6,
    BOBOL_LOD_RESULT_DIAGNOSTIC = 7
};

enum BObolLodProviderStatus {
    BOBOL_LOD_PROVIDER_UNKNOWN = 0,
    BOBOL_LOD_PROVIDER_READY = 1,
    BOBOL_LOD_PROVIDER_CACHE_MISS = 2,
    BOBOL_LOD_PROVIDER_STALE = 3,
    BOBOL_LOD_PROVIDER_RUNNING = 4,
    BOBOL_LOD_PROVIDER_TERMINAL = 5,
    BOBOL_LOD_PROVIDER_FALLBACK = 6,
    BOBOL_LOD_PROVIDER_ERROR = 7,
    BOBOL_LOD_PROVIDER_CANCELLED = 8
};

enum BObolLodGeometryHandleKind {
    BOBOL_LOD_GEOMETRY_NONE = 0,
    BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE = 1,
    BOBOL_LOD_GEOMETRY_OBOL_MESH = 2,
    BOBOL_LOD_GEOMETRY_PROVIDER_TOKEN = 3
};

enum BObolLodProxyKind {
    BOBOL_LOD_PROXY_NONE = 0,
    BOBOL_LOD_PROXY_AABB = 1,
    BOBOL_LOD_PROXY_OBB = 2,
    BOBOL_LOD_PROXY_SPHERE = 3,
    BOBOL_LOD_PROXY_PROVIDER_TOKEN = 4
};

struct BOBOL_EXPORT BObolLodProviderParam {
    SbString name;
    SbString value;
};

struct BOBOL_EXPORT BObolLodDependency {
    SbString objectPath;
    SbString objectName;
    uint64_t sourceRevision;
    uint64_t sourceContentHash;
    int requiredQualityTier;
    SbBool optional;

    BObolLodDependency(void);
    void clear(void);
};

struct BOBOL_EXPORT BObolLodAttribute {
    SbString name;
    SbString value;

    BObolLodAttribute(void);
    void clear(void);
};

struct BOBOL_EXPORT BObolLodCounts {
    uint64_t faceCount;
    uint64_t pointCount;
    uint64_t originalPointCount;
    uint64_t normalCount;
    uint64_t lineCount;
    uint64_t byteCount;

    BObolLodCounts(void);
    void clear(void);
};

struct BOBOL_EXPORT BObolLodCacheKey {
    SbString value;

    BObolLodCacheKey(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolLodResidentDemand {
    SbString assetKey;
    int level;

    BObolLodResidentDemand(void);
};

struct BOBOL_EXPORT BObolLodGeometryHandle {
    int kind;
    SbString providerId;
    SbString providerVersion;
    BObolLodCacheKey cacheKey;
    uint64_t providerToken;
    int activeLevel;
    SbBool borrowed;

    BObolLodGeometryHandle(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolLodMeshPayload {
    std::vector<SbVec3f> points;
    std::vector<SbVec3f> normals;
    std::vector<int32_t> coordIndex;
    std::vector<int32_t> faceIndex;
    std::vector<int32_t> vertexIndex;

    BObolLodMeshPayload(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BObolLodProgressiveMeshPrivate;

/* One thread-safe retained PoP asset shared by every occurrence and view that
 * resolves to the same source geometry.  The resident arrays are exact,
 * activation-ordered cumulative prefixes.  active/draw level deliberately
 * does not live here: each occurrence may draw a different prefix. */
class BOBOL_EXPORT BObolLodProgressiveMesh {
public:
    BObolLodProgressiveMesh(void);
    ~BObolLodProgressiveMesh(void);

    SbBool update(const struct BObolMeshLodData &data,
	const struct BObolMeshLodHierarchyInfo &hierarchy,
	int residentLevel, SbBool shadedCullBackfaces);
    SbBool trim(int residentLevel);
    SbBool copyLevel(BObolLodMeshPayload &payload, int level) const;
    SbBool isValid(void) const;
    /* True when the retained prefix already contains every point and face
     * needed by level.  This may be true above residentLevel() when adjacent
     * PoP levels differ only in coordinate quantization: exact coordinates
     * are retained once and the renderer can select the finer snap without
     * loading or rebuilding geometry. */
    SbBool canDrawLevel(int level) const;
    int minimumLevel(void) const;
    int maximumLevel(void) const;
    int residentLevel(void) const;
    uint64_t revision(void) const;
    size_t pointCount(int level) const;
    size_t faceCount(int level) const;
    size_t estimateBytes(void) const;
    SbBox3f bounds(void) const;
    SbVec3f quantizationMinimum(void) const;
    SbVec3f quantizationMaximum(void) const;
    SbBool cullBackfaces(void) const;

private:
    BObolLodProgressiveMesh(const BObolLodProgressiveMesh &);
    BObolLodProgressiveMesh &operator=(
	const BObolLodProgressiveMesh &);
    BObolLodProgressiveMeshPrivate *p;
};

typedef std::shared_ptr<BObolLodProgressiveMesh>
    BObolLodProgressiveMeshPtr;

struct BOBOL_EXPORT BObolLodProxy {
    int kind;
    SbBox3f bounds;
    SbVec3f center;
    SbVec3f axisX;
    SbVec3f axisY;
    SbVec3f axisZ;
    SbVec3f halfExtents;
    float radius;
    BObolLodGeometryHandle geometry;

    BObolLodProxy(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolLodRequest {
    SbString databaseId;
    uint64_t databaseRevision;
    uint64_t sourceRevision;
    uint64_t sourceContentHash;
    SbString objectPath;
    SbString objectName;
    /* Stable compact-occurrence identity.  Empty for legacy/source-wide
     * requests; populated whenever a result targets one CAD occurrence. */
    SbString occurrenceKey;
    uint64_t viewRevision;
    uint64_t policyRevision;
    int drawMode;
    uint32_t lodPolicy;
    SbString providerId;
    SbString providerVersion;
    int qualityTier;
    /* View-derived display demand.  requestedLevel is a PoP quantization
     * level, never an exact/full-detail request.  A negative value means the
     * producer did not have a projectable view and the provider should use its
     * conservative legacy view estimate. */
    float projectedPixelDiameter;
    float targetPixelError;
    int requestedLevel;
    SbBox3f bounds;
    BObolLodCounts sourceCounts;
    std::vector<BObolLodProviderParam> providerParams;

    BObolLodRequest(void);
    void clear(void);
    void addProviderParam(const SbString &name, const SbString &value);
};

struct BOBOL_EXPORT BObolLodResult {
    uint64_t generation;
    BObolLodRequest request;
    BObolLodCacheKey cacheKey;
    BObolLodGeometryHandle geometry;
    BObolLodMeshPayload mesh;
    BObolLodProgressiveMeshPtr progressiveMesh;
    int residentLevel;
    int resultKind;
    int qualityTier;
    int providerStatus;
    SbBox3f bounds;
    BObolLodCounts counts;
    std::vector<BObolLodDependency> dependencies;
    std::vector<BObolLodAttribute> attributes;
    BObolLodProxy proxy;
    double estimatedError;
    SbBool terminal;
    SbBool fallback;
    SbBool stale;
    SbBool hasSnappedPoints;
    SbBool hasNormals;
    SbBool shadedCullBackfaces;
    SbString diagnostic;

    BObolLodResult(void);
    void clear(void);
    void addDependency(const SbString &objectPath, const SbString &objectName,
	uint64_t sourceRevision, uint64_t sourceContentHash,
	int requiredQualityTier, SbBool optional = FALSE);
    void addAttribute(const SbString &name, const SbString &value);
};

BOBOL_EXPORT BObolLodCacheKey
bobol_lod_cache_key(const BObolLodRequest &request);

/* Stable identity of provider geometry.  Unlike bobol_lod_cache_key this does
 * not contain occurrence or camera epochs, so an unchanged level remains the
 * same display asset across view changes. */
BOBOL_EXPORT BObolLodCacheKey
bobol_lod_geometry_cache_key(const BObolLodRequest &request);

/* Stable source-asset identity.  It excludes occurrence, camera, requested
 * level, and draw mode, so all consumers share one residency high-water mark. */
BOBOL_EXPORT BObolLodCacheKey
bobol_lod_asset_cache_key(const BObolLodRequest &request);

BOBOL_EXPORT SbBool
bobol_lod_mesh_payload_from_mesh_lod_data(BObolLodMeshPayload &payload,
	const struct BObolMeshLodData &data);

BOBOL_EXPORT SbBool
bobol_lod_result_matches_request(const BObolLodResult &result,
	const BObolLodRequest &request);

BOBOL_EXPORT BObolLodResult
bobol_lod_result_from_mesh_lod_info(const BObolLodRequest &request,
	const struct BObolMeshLodInfo &info,
	const struct BObolMeshLodCacheStatus *status);

BOBOL_EXPORT BObolLodResult
bobol_lod_directory_result(const BObolLodRequest &request,
	const std::vector<BObolLodDependency> &dependencies);

BOBOL_EXPORT BObolLodResult
bobol_lod_attributes_result(const BObolLodRequest &request,
	const std::vector<BObolLodAttribute> &attributes);

BOBOL_EXPORT BObolLodResult
bobol_lod_aabb_result(const BObolLodRequest &request,
	const SbBox3f &bounds, const BObolLodCounts *counts);

BOBOL_EXPORT BObolLodResult
bobol_lod_proxy_result(const BObolLodRequest &request,
	const BObolLodProxy &proxy, const BObolLodCounts *counts);

#endif /* BOBOL_BLODREALIZATION_H */
