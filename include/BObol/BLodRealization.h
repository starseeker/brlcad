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
#include "BObol/BLodIdentifiers.h"
#include "BObol/BMeshLodCache.h"

#include <Inventor/SbBox.h>
#include <Inventor/SbMatrix.h>
#include <Inventor/SbString.h>
#include <Inventor/SbVec3f.h>
#include <memory>
#include <stdint.h>
#include <vector>

namespace Obol {
struct PartGeometry;
}

enum BObolLodDrawMode {
    BOBOL_LOD_DRAW_UNKNOWN = 0,
    BOBOL_LOD_DRAW_WIRE = 1,
    BOBOL_LOD_DRAW_SHADED = 2,
    BOBOL_LOD_DRAW_POINTS = 3,
    BOBOL_LOD_DRAW_DIAGNOSTIC = 4,
    BOBOL_LOD_DRAW_SHADED_BOTS = 5,
    BOBOL_LOD_DRAW_HIDDEN_LINE = 6
};

/* Renderer presentation policy captured with a view demand.  This is not
 * source/cache identity: one resident PoP asset may prepare several immutable
 * view-specific normal variants without duplicating persistent geometry. */
enum BObolLodNormalStyle {
    BOBOL_LOD_NORMAL_AUTHORED = 0,
    BOBOL_LOD_NORMAL_FLAT = 1,
    BOBOL_LOD_NORMAL_SMOOTH = 2
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

enum BObolLodPayloadKind {
    BOBOL_LOD_PAYLOAD_NONE = 0,
    BOBOL_LOD_PAYLOAD_DIRECTORY = 1,
    BOBOL_LOD_PAYLOAD_ATTRIBUTES = 2,
    BOBOL_LOD_PAYLOAD_PROXY = 3,
    BOBOL_LOD_PAYLOAD_MESH = 4,
    BOBOL_LOD_PAYLOAD_STATUS = 5
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
    BOBOL_LOD_PROVIDER_CANCELLED = 8,
    /* Valid work which no longer satisfies the service-owned active demand.
     * This is a scheduling outcome, not a source/cache failure. */
    BOBOL_LOD_PROVIDER_SUPERSEDED = 9
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
    BObolSourceEpoch sourceRevision;
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

/*
 * Renderer-neutral scheduling currency, expressed in shaded-triangle
 * equivalents.  It accounts for vertex processing, line expansion, normal
 * bandwidth, and the irreducible per-occurrence planning cost.  The value is
 * deliberately cheap to compute from cached hierarchy counts.
 */
BOBOL_EXPORT size_t
bobol_lod_render_cost_units(const BObolLodCounts &counts, int drawMode,
    size_t occurrenceCount = 1);

/* Physical cost of one renderer-batched aggregate proxy.  Points use the
 * shared one-vertex stream; boxes follow the requested CAD draw channels. */
BOBOL_EXPORT size_t
bobol_lod_aggregate_proxy_render_cost(SbBool box, int drawMode);

struct BOBOL_EXPORT BObolLodCacheKey {
    SbString value;

    BObolLodCacheKey(void);
    void clear(void);
    SbBool isValid(void) const;
};

struct BOBOL_EXPORT BObolLodResidentDemand {
    SbString assetKey;
    int cut;
    /* Renderer channels which retain this asset in the consumer.  Bit 0 is
     * wire geometry and bit 1 is shaded geometry.  Stable compaction uses the
     * aggregate mask to prepare one replacement immutable renderer object on
     * a worker rather than rebuilding it on the presentation thread. */
    unsigned int channelMask;
    /* Sorted union of private spatial subresources needed by this consumer.
     * Empty means the asset is not chunked, never "all chunks". */
    std::vector<uint32_t> chunkIds;

    BObolLodResidentDemand(void);
};

struct BOBOL_EXPORT BObolLodGeometryHandle {
    int kind;
    SbString providerId;
    SbString providerVersion;
    BObolLodCacheKey cacheKey;
    uint64_t providerToken;
    int activeCut;
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
struct BObolLodProgressiveMeshGeneration;
struct BObolLodProgressiveMeshTrim;
typedef std::shared_ptr<const BObolLodProgressiveMeshTrim>
    BObolLodProgressiveMeshTrimPtr;

/* One retained view of immutable progressive-selection metadata.  Allocation
 * transactions use this object to make every decision against the same mesh
 * generation without repeatedly acquiring its atomic shared pointer. */
class BOBOL_EXPORT BObolLodProgressiveMeshSnapshot {
public:
    BObolLodProgressiveMeshSnapshot(void);

    SbBool isValid(void) const;
    uint64_t revision(void) const;
    int minimumCut(void) const;
    int maximumCut(void) const;
    BObolLodCounts hierarchyCountsAtCut(int cut, SbBool hasNormals) const;
    int cutForScreenError(double projectedPixelDiameter,
	double targetPixelError) const;
    double projectedErrorAtCut(int cut,
	double projectedPixelDiameter) const;

private:
    friend class BObolLodProgressiveMesh;
    explicit BObolLodProgressiveMeshSnapshot(
	const std::shared_ptr<const BObolLodProgressiveMeshGeneration> &);
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation;
};

/* One immutable renderer layer owned by an LoD result.  A layer is a
 * presentation subresource of one logical occurrence, never a second CAD
 * occurrence.  layerKey is stable within the source-asset partition and is
 * used by the retained assembly to derive its private part/instance ids. */
struct BOBOL_EXPORT BObolLodPresentationLayer {
    SbString layerKey;
    std::shared_ptr<const Obol::PartGeometry> geometry;
    uint64_t geometryRevision;
    int activeCut;
    SbBool coverage;

    BObolLodPresentationLayer(void);
    SbBool isValid(void) const;
    size_t estimateBytes(void) const;
};

/* Population frontier of one private spatial page.  Records are always
 * sorted by chunkId when exchanged with BObolLodProgressiveMesh. */
struct BOBOL_EXPORT BObolLodChunkCut {
    uint32_t chunkId = 0;
    int cut = -1;

    bool operator==(const BObolLodChunkCut &other) const
    {
	return chunkId == other.chunkId && cut == other.cut;
    }
};

/* One thread-safe retained PoP asset shared by every occurrence and view that
 * resolves to the same source geometry.  A small leaf uses one exact,
 * activation-ordered cumulative prefix.  A large leaf uses independently
 * resident private page prefixes published as immutable renderer layers
 * behind one logical CAD identity.
 * The active draw cut deliberately does not live here: each occurrence may
 * draw a different cut from the pages its view requires. */
class BOBOL_EXPORT BObolLodProgressiveMesh {
public:
    BObolLodProgressiveMesh(void);
    ~BObolLodProgressiveMesh(void);

    SbBool update(const struct BObolMeshLodData &data,
	const struct BObolMeshLodHierarchyInfo &hierarchy,
	int residentCut, SbBool shadedCullBackfaces);
    /* Populate only the selected private spatial pages.  Existing immutable
     * pages are shared into the replacement generation; missing/richer pages
     * are read independently from the cache. */
    SbBool updateChunksFromCache(struct BObolMeshLod *lod,
	const struct BObolMeshLodHierarchyInfo &hierarchy,
	const std::vector<uint32_t> &chunkIds, int residentCut,
	SbBool shadedCullBackfaces);
    /* Extend an existing immutable generation by reading only cache cuts
     * above its published resident frontier.  The old generation remains
     * drawable until the completed replacement is atomically published.
     * Authored corner-normal meshes currently return FALSE and use update(),
     * whose global vertex splitting cannot yet be reconstructed from a suffix
     * alone. */
    SbBool extendFromCache(struct BObolMeshLod *lod,
	const struct BObolMeshLodHierarchyInfo &hierarchy,
	int residentCut, SbBool shadedCullBackfaces);
    /* Build a shorter immutable generation without publishing it.  Stable
     * resident-memory maintenance uses this two-phase form so a newer view or
     * provider use can invalidate an old trim before it changes the shared
     * mesh.  commitTrim succeeds only while the exact source generation used
     * by prepareTrim is still current. */
    BObolLodProgressiveMeshTrimPtr prepareTrim(int residentCut) const;
    BObolLodProgressiveMeshTrimPtr prepareTrim(
	int residentCut, const std::vector<uint32_t> &chunkIds) const;
    /* Dense-repack the exact page/cut working set.  This is the stable-memory
     * operation used after multiple views have contributed different detail
     * demands for different pages of one logical leaf. */
    BObolLodProgressiveMeshTrimPtr prepareTrim(
	const std::vector<BObolLodChunkCut> &chunkCuts) const;
    SbBool commitTrim(const BObolLodProgressiveMeshTrimPtr &trim);
    SbBool trim(int residentCut);
    SbBool copyCut(BObolLodMeshPayload &payload, int cut) const;
    SbBool isValid(void) const;
    /* True when the retained prefix already contains every point and face
     * needed by cut.  This may be true above residentCut() when adjacent
     * PoP cuts differ only in coordinate quantization: exact coordinates
     * are retained once and the renderer can select the finer snap without
     * loading or rebuilding geometry. */
    SbBool canDrawCut(int cut) const;
    SbBool canDrawChunksAtCut(const std::vector<uint32_t> &chunkIds,
	int cut) const;
    SbBool countsForChunksAtCut(const std::vector<uint32_t> &chunkIds,
	int cut, SbBool hasNormals, BObolLodCounts *counts) const;
    /* Snapshot the immutable generation once and report the hierarchy
     * population for every producer cut.  Empty chunkIds means the complete
     * logical leaf; a nonempty sorted set selects private spatial pages.
     * maximumDrawableCut is the richest cumulative prefix which the retained
     * generation can submit coherently.  This bulk form is intended for
     * scene-budget accounting, where independently querying every cut would
     * otherwise acquire the same shared generation dozens of times. */
    SbBool drawableCountsAtCuts(const std::vector<uint32_t> &chunkIds,
	SbBool hasNormals, BObolLodCounts *counts, size_t count,
	int *minimumCut, int *maximumDrawableCut) const;
    SbBool hierarchyCountsForChunksAtCut(
	const std::vector<uint32_t> &chunkIds, int cut,
	SbBool hasNormals, BObolLodCounts *counts) const;
    /* Return the sorted subset which contributes triangles at cut.  Spatial
     * pages with no population at a coarse cut deliberately have no renderer
     * layer; callers use this census to distinguish that valid omission from
     * an incomplete presentation. */
    SbBool populatedChunkIdsAtCut(const std::vector<uint32_t> &chunkIds,
	int cut, std::vector<uint32_t> &populatedChunkIds) const;
    int residentCutForChunks(const std::vector<uint32_t> &chunkIds) const;
    void residentChunkIds(std::vector<uint32_t> &chunkIds) const;
    void residentChunkCuts(std::vector<BObolLodChunkCut> &chunkCuts) const;
    int minimumCut(void) const;
    int maximumCut(void) const;
    int residentCut(void) const;
    uint64_t revision(void) const;
    BObolLodProgressiveMeshSnapshot snapshot(void) const;
    size_t pointCount(int cut) const;
    size_t faceCount(int cut) const;
    /* Immutable hierarchy population, including cuts above the currently
     * resident prefix.  pointCount()/faceCount() retain their drawable-cut
     * semantics and clamp to residentCut(). */
    size_t hierarchyPointCountAtCut(int cut) const;
    size_t hierarchyFaceCountAtCut(int cut) const;
    SbBool hasSpatialClusters(void) const;
    /* Compute the submitted population for every cut when this occurrence
     * straddles the render frustum.  Spatial clusters remain private pieces
     * of one logical CAD leaf; this query is scheduling metadata only.  It
     * returns FALSE for unclustered or wholly-contained occurrences, in
     * which case callers must use the ordinary whole-prefix counts. */
    SbBool visibleCountsAtCuts(const SbMatrix &localToRoot,
	const SbMatrix &viewProjection, SbBool hasNormals,
	BObolLodCounts *counts, size_t count) const;
    SbBool visibleChunkIds(const SbMatrix &localToRoot,
	const SbMatrix &viewProjection,
	std::vector<uint32_t> &chunkIds) const;
    /* Copy immutable selection metadata for one producer-defined cut. */
    SbBool cutInfo(int cut, struct BObolMeshLodCutInfo *info) const;
    int cutForScreenError(double projectedPixelDiameter,
	double targetPixelError) const;
    double projectedErrorAtCut(int cut,
	double projectedPixelDiameter) const;
    size_t estimateBytes(void) const;
    SbBox3f bounds(void) const;
    SbVec3f quantizationMinimum(void) const;
    SbVec3f quantizationMaximum(void) const;
    SbBool cullBackfaces(void) const;
    /* Build one immutable renderer-ready snapshot on the calling worker.
     * The returned geometry is revision tagged and may be adopted directly by
     * SoCADAssembly without copying its large arrays on the GUI thread.
     * Authored corner normals are split into Obol's indexed vertex/normal
     * representation by that worker before publication. */
    std::shared_ptr<const Obol::PartGeometry> prepareCadGeometry(
	int drawMode, uint64_t *preparedRevision = NULL) const;
    /* Publish selected spatial pages without flattening them into one large
     * allocation.  Every returned layer uses activeCut, so neighboring pages
     * retain a coherent quantization boundary while their immutable storage
     * may have different resident high-water marks.  A valid page set may
     * produce an empty layer vector when none of its faces is active at the
     * requested cut; that is a successful zero-draw presentation. */
    SbBool prepareCadPresentationLayers(int drawMode,
	const std::vector<BObolLodChunkCut> &chunkCuts, int normalStyle,
	float normalCreaseAngle,
	std::vector<BObolLodPresentationLayer> &layers) const;

private:
    BObolLodProgressiveMesh(const BObolLodProgressiveMesh &);
    BObolLodProgressiveMesh &operator=(
	const BObolLodProgressiveMesh &);
    BObolLodProgressiveMeshPrivate *p;
};

typedef std::shared_ptr<BObolLodProgressiveMesh>
    BObolLodProgressiveMeshPtr;

/* Construct the complete renderer-cost population for one cumulative PoP
 * cut.  All scene admission paths must use this helper rather than partially
 * initializing BObolLodCounts: omitting authored normals from a predicted cut
 * while charging them to the active payload makes an apparently feasible
 * scene allocation exceed its budget as soon as it is published. */
BOBOL_EXPORT BObolLodCounts
bobol_lod_progressive_counts(
    const BObolLodProgressiveMeshPtr &progressiveMesh, int cut,
    SbBool hasNormals);

/* Resolve the private chunk set intersecting one occurrence's render
 * frustum.  Returns FALSE for an unchunked hierarchy or invalid projection. */
BOBOL_EXPORT SbBool
bobol_lod_visible_chunks(
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    const SbMatrix &localToRoot,
    const SbMatrix &viewProjection,
    std::vector<uint32_t> &chunkIds);

struct BOBOL_EXPORT BObolLodResidentCompaction {
    SbString assetKey;
    BObolLodProgressiveMeshPtr progressiveMesh;
    /* Complete consumer-demand snapshot which authorized this trim.  A view
     * must reject the completion after any occurrence demand changes: quiet
     * memory maintenance is not presentation authority. */
    uint64_t consumerDemandRevision;
    std::shared_ptr<const Obol::PartGeometry> preparedCadGeometry;
    uint64_t preparedCadGeometryRevision;
    /* Chunked assets retain page geometry rather than rebuilding one packed
     * mesh.  The owner filters these immutable handles for each occurrence's
     * current page set and active cut. */
    std::vector<BObolLodPresentationLayer> presentationLayers;
    int residentCut;
    unsigned int channelMask;
    size_t priorBytes;
    size_t residentBytes;

    BObolLodResidentCompaction(void);
};

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
    BObolDatabaseEpoch databaseRevision;
    BObolSourceEpoch sourceRevision;
    uint64_t sourceContentHash;
    SbString objectPath;
    SbString objectName;
    /* Stable compact-occurrence identity.  Empty for legacy/source-wide
     * requests; populated whenever a result targets one CAD occurrence. */
    SbString occurrenceKey;
    /*
     * Owner-thread routing token for a compact database source.  It is not
     * geometry/cache identity and deliberately does not participate in cache
     * keys.  Results use it to reach their retained source directly instead
     * of traversing the scene graph to rediscover an owner already known at
     * submission time.
     */
    BObolSourceRoutingId sourceRoutingId;
    /* Generation of the dense compact-entry population addressed by
     * sourceEntryIndex.  This is owner-thread routing metadata, not geometry
     * or persistent-cache identity.  A compact occurrence request must carry
     * a nonzero value so a late worker result cannot bind to a recycled entry
     * after erase/redraw or edit-driven registry replacement. */
    BObolSourcePopulationEpoch sourcePopulationEpoch;
    /* Multiple compact occurrences may consume one immutable source asset.
     * TRUE lets the service serialize construction/loading of that asset
     * across those occurrences; the owner-thread planner binds or resubmits
     * each occurrence after the producer result is published.  Route and
     * population epochs remain distinct, so separate scene sources cannot
     * consume one another's results. */
    SbBool coalesceAssetProducer;
    /* Fixed-width source-local compact index.  UINT32_MAX means the request
     * is not backed by a compact entry.  The occurrence key remains semantic
     * identity across index rebuilds; owner-thread result publication
     * validates both values against sourcePopulationEpoch before installing
     * the dense binding.  Within that epoch an empty slot is authoritative
     * negative evidence and does not require a string-table fallback. */
    uint32_t sourceEntryIndex;
    BObolViewEpoch viewRevision;
    BObolPolicyEpoch policyRevision;
    /* View-local presentation emphasis captured by the same bounded compact
     * census which projected this request: 0 ordinary, 1 highlighted,
     * 2 selected.  It is scheduling metadata and not cache identity. */
    uint8_t visualEmphasis;
    int drawMode;
    int normalStyle;
    float normalCreaseAngle;
    uint32_t lodPolicy;
    SbString providerId;
    SbString providerVersion;
    int qualityTier;
    /* View-derived display demand.  requestedCut is an admissible PoP cut.
     * A negative value means no producer hierarchy was available when the
     * request was projected.  The provider must select from its certified cut
     * metadata using projectedPixelDiameter and targetPixelError; if the view
     * is not projectable it uses the hierarchy's minimum useful cut. */
    float projectedPixelDiameter;
    /* Projected convex-hull footprint of the occurrence bound.  Diameter is
     * the conservative geometric-error scale; area and perimeter distinguish
     * a genuinely prominent surface/silhouette from a long, mostly empty
     * screen-aligned bounding rectangle. */
    float projectedPixelArea;
    float projectedPixelPerimeter;
    /* TRUE only when the complete conservative bound is inside the renderer
     * clip volume.  SoCADAssembly deliberately will not replace a partly
     * clipped occurrence by one point: doing so can move the visible sliver
     * or make it disappear.  The scheduler must use the same fact when it
     * decides whether a structural proxy is already terminal coverage. */
    SbBool projectedBoundsContained;
    float targetPixelError;
    int requestedCut;
    /* Exact projection state used to resolve a cold chunked hierarchy.  It is
     * request-local scheduling data, not asset identity. */
    SbMatrix localToRoot;
    SbMatrix viewProjection;
    SbBool spatialProjectionValid;
    std::vector<uint32_t> requiredChunks;
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
    std::shared_ptr<const Obol::PartGeometry> preparedCadGeometry;
    uint64_t preparedCadGeometryRevision;
    /* Ordered immutable subresources for one logical CAD occurrence.  Empty
     * retains the ordinary single-geometry representation during migration. */
    std::vector<BObolLodPresentationLayer> presentationLayers;
    /* Producer-resolved display demand.  A cold request may carry
     * request.requestedCut == -1 because no hierarchy existed at submission;
     * the provider selects from certified metadata and reports that cut here
     * without mutating the task identity copied into request. */
    int resolvedCut;
    int residentCut;
    /* Transient owner-budget certificate for geometry.activeCut.  This is
     * intentionally independent of request.viewRevision: coalescing may
     * update the worker request to a newer camera while the provider retains
     * the budget under which its presentation cut was selected. */
    SbBool presentationAdmissionCertified;
    BObolViewEpoch presentationAdmissionViewRevision;
    BObolPolicyEpoch presentationAdmissionPolicyRevision;
    int presentationAdmissionCut;
    /* The service could publish a useful retained mesh, but deliberately
     * withheld a richer suffix to honor its stable resident-byte target.
     * residentAdmissionRevision identifies the capacity epoch which made
     * that decision so an unchanged view cannot immediately retry it. */
    uint64_t residentAdmissionRevision;
    /* Authoritative payload category.  resultKind describes the semantic
     * quality stage; payloadKind makes the active storage alternative
     * explicit and is canonicalized before service publication. */
    int payloadKind;
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
    /* Diagnostic compatibility flag for an unusable/outdated result.  Use
     * providerStatus to distinguish invalid source data from superseded work. */
    SbBool stale;
    SbBool hasSnappedPoints;
    SbBool hasNormals;
    SbBool shadedCullBackfaces;
    SbBool memoryLimited;
    SbString diagnostic;

    BObolLodResult(void);
    void clear(void);
    /* Select the payload category implied by provider status/resultKind and
     * release every inactive large channel. */
    void canonicalizePayload(void);
    SbBool payloadIsConsistent(void) const;
    void addDependency(const SbString &objectPath, const SbString &objectName,
	uint64_t sourceRevision, uint64_t sourceContentHash,
	int requiredQualityTier, SbBool optional = FALSE);
    void addAttribute(const SbString &name, const SbString &value);
};

BOBOL_EXPORT BObolLodCacheKey
bobol_lod_cache_key(const BObolLodRequest &request);

/* Stable identity of provider geometry.  Unlike bobol_lod_cache_key this does
 * not contain occurrence or camera epochs, so an unchanged cut remains the
 * same display asset across view changes. */
BOBOL_EXPORT BObolLodCacheKey
bobol_lod_geometry_cache_key(const BObolLodRequest &request);

/* Stable source-asset identity.  It excludes occurrence, camera, requested
 * cut, and draw mode, so all consumers share one residency high-water mark. */
BOBOL_EXPORT BObolLodCacheKey
bobol_lod_asset_cache_key(const BObolLodRequest &request);

BOBOL_EXPORT SbBool
bobol_lod_mesh_payload_from_mesh_lod_data(BObolLodMeshPayload &payload,
	const struct BObolMeshLodData &data);

BOBOL_EXPORT SbBool
bobol_lod_result_matches_request(const BObolLodResult &result,
	const BObolLodRequest &request);

/**
 * Allocation-free equality of the complete structured request key.
 *
 * Projected diameter and target pixel error are demand diagnostics.  A
 * nonnegative requestedCut is part of provider identity; a cold request may
 * retain -1 while its result reports the producer selection separately in
 * resolvedCut.  Provider parameters are compared as an order-independent
 * multiset, matching cache serialization.
 */
BOBOL_EXPORT SbBool
bobol_lod_request_keys_equal(const BObolLodRequest &lhs,
	const BObolLodRequest &rhs);

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
