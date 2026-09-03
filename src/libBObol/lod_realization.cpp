/*              L O D _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/datetime.h"
#include "bu/parallel.h"
#include "bu/str.h"

#include "BObol/BLodRealization.h"
#include "cad_publication_private.h"
#include "identity_counter_private.h"

#include <Obol/cad/CadGeometryValidation.h>
#include <Obol/cad/CadProjectedProxy.h>
#include <Obol/cad/SoCADAssembly.h>

#include <Inventor/SbVec3f.h>

#include <algorithm>
#include <array>
#include <atomic>
#include <cfloat>
#include <charconv>
#include <cmath>
#include <limits>
#include <mutex>
#include <optional>
#include <string>
#include <string.h>
#include <system_error>
#include <thread>
#include <unordered_map>

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

static size_t
lod_render_cost_add(size_t total, uint64_t count, size_t numerator,
    size_t denominator = 1)
{
    if (total == SIZE_MAX || !count || !numerator)
	return total;
    const uint64_t scaled = count > UINT64_MAX / numerator ?
	UINT64_MAX : count * numerator;
    const uint64_t rounded = scaled / denominator +
	(scaled % denominator != 0 ? 1 : 0);
    const size_t value = rounded > static_cast<uint64_t>(SIZE_MAX) ?
	SIZE_MAX : static_cast<size_t>(rounded);
    return value > SIZE_MAX - total ? SIZE_MAX : total + value;
}

size_t
bobol_lod_render_cost_units(const BObolLodCounts &counts, int drawMode,
    size_t occurrenceCount)
{
    size_t cost = 0;
    /* Retained batching greatly reduces draw calls, but visibility, instance
     * records, selection state, and command construction still impose a real
     * floor for every occurrence. */
    cost = lod_render_cost_add(cost, occurrenceCount, 64);
    cost = lod_render_cost_add(cost, counts.pointCount, 1, 4);
    /* A retained asset may own authored normals while a wire-only view never
     * submits or reads them.  Charge channel work, not dormant asset data. */
    if (drawMode == BOBOL_LOD_DRAW_SHADED ||
	drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	cost = lod_render_cost_add(cost, counts.normalCount, 1, 8);

    const uint64_t lines = counts.lineCount ? counts.lineCount :
	(counts.faceCount > UINT64_MAX / 3 ? UINT64_MAX :
	 counts.faceCount * 3);
    switch (drawMode) {
	case BOBOL_LOD_DRAW_WIRE:
	    cost = lod_render_cost_add(cost, lines, 1, 2);
	    break;
	case BOBOL_LOD_DRAW_HIDDEN_LINE:
	    cost = lod_render_cost_add(cost, counts.faceCount, 1);
	    cost = lod_render_cost_add(cost, lines, 1, 2);
	    break;
	default:
	    cost = lod_render_cost_add(cost, counts.faceCount, 1);
	    break;
    }
    return cost;
}

size_t
bobol_lod_aggregate_proxy_render_cost(SbBool box, int drawMode)
{
    BObolLodCounts counts;
    if (!box) {
	counts.pointCount = 1;
	return bobol_lod_render_cost_units(
	    counts, BOBOL_LOD_DRAW_WIRE, 0);
    }

    const bool knownMode = drawMode == BOBOL_LOD_DRAW_WIRE ||
	drawMode == BOBOL_LOD_DRAW_SHADED ||
	drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    const int effectiveDrawMode = knownMode ?
	drawMode : BOBOL_LOD_DRAW_WIRE;
    const bool wire = effectiveDrawMode == BOBOL_LOD_DRAW_WIRE ||
	effectiveDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    const bool shaded = effectiveDrawMode == BOBOL_LOD_DRAW_SHADED ||
	effectiveDrawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	effectiveDrawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    if (wire) {
	counts.pointCount += Obol::CadAggregateProxyBoxPositionCount;
	counts.lineCount = Obol::CadAggregateProxyBoxLineCount;
    }
    if (shaded) {
	counts.pointCount +=
	    Obol::CadAggregateProxyBoxTrianglePositionCount;
	counts.normalCount +=
	    Obol::CadAggregateProxyBoxTrianglePositionCount;
	counts.faceCount = Obol::CadAggregateProxyBoxTriangleCount;
    }
    return bobol_lod_render_cost_units(counts, effectiveDrawMode, 0);
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
    cut(-1),
    channelMask(0)
{
}

BObolLodResidentCompaction::BObolLodResidentCompaction(void) :
    consumerDemandRevision(0),
    preparedCadGeometryRevision(0),
    residentCut(-1),
    channelMask(0),
    priorBytes(0),
    residentBytes(0)
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
    activeCut = -1;
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

struct BObolLodProgressiveMeshGeneration {
    /* The shaded renderer object is the canonical retained population, not a
     * second conversion layered over an intermediate mesh payload.  Workers
     * build it directly from a cache prefix and the scene adopts this exact
     * immutable allocation. */
    std::shared_ptr<const Obol::PartGeometry> shadedGeometry;
    std::array<size_t, BOBOL_MESH_LOD_CUT_COUNT_MAX> pointCount = {};
    std::array<size_t, BOBOL_MESH_LOD_CUT_COUNT_MAX> faceCount = {};
    std::vector<BObolMeshLodCutInfo> cuts;
    SbBox3f bounds;
    SbVec3f quantizationMinimum;
    SbVec3f quantizationMaximum;
    std::optional<std::array<SbVec3f, 8>> aggregateProxyCorners;
    int minimumCut = -1;
    int maximumCut = -1;
    int residentCut = -1;
    uint64_t revision = 0;
    SbBool shadedCullBackfaces = FALSE;
    struct ChunkMetadata {
	uint32_t chunkId = 0;
	int minimumCut = -1;
	int maximumCut = -1;
	SbBox3f bounds;
	std::array<BObolMeshLodChunkCutInfo,
	    BOBOL_MESH_LOD_CUT_COUNT_MAX> cuts = {};
    };
    struct ResidentChunk {
	uint32_t chunkId = 0;
	int residentCut = -1;
	uint64_t geometryRevision = 0;
	SbBox3f bounds;
	std::vector<SbVec3f> positions;
	std::vector<uint32_t> indices;
	std::vector<SbVec3f> cornerNormals;
	mutable std::mutex preparedMutex;
	mutable std::unordered_map<int,
	    std::weak_ptr<const Obol::PartGeometry>> preparedCadGeometry;
    };
    std::vector<ChunkMetadata> chunks;
    std::vector<std::shared_ptr<const ResidentChunk>> residentChunks;
    /*
     * Renderer-ready channel combinations are immutable generation products.
     * The cache is weak because the scene is the owner once a result is
     * published; dropping a scene must make the large arrays reclaimable.
     */
    mutable std::mutex preparedMutex;
    mutable std::unordered_map<int,
	std::weak_ptr<const Obol::PartGeometry>>
	preparedCadGeometry;
};

static uint64_t
bobol_next_progressive_lineage(void)
{
    static std::atomic<uint64_t> next(1);
    return bobol_atomic_nonzero_identity_take(next);
}

struct BObolLodProgressiveMeshPrivate {
    /*
     * Only writers serialize here.  Readers atomically retain one immutable
     * generation and never wait while a richer PoP prefix is converted.
     */
    std::mutex updateMutex;
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation;
    const uint64_t lineage = bobol_next_progressive_lineage();
};

static std::shared_ptr<const BObolLodProgressiveMeshGeneration>
progressive_generation_load(const BObolLodProgressiveMeshPrivate *p)
{
    if (!p)
	return std::shared_ptr<const BObolLodProgressiveMeshGeneration>();
    return std::atomic_load_explicit(
	&p->generation, std::memory_order_acquire);
}

static void
progressive_generation_store(
    BObolLodProgressiveMeshPrivate *p,
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> &generation)
{
    if (!p)
	return;
    std::atomic_store_explicit(
	&p->generation, generation, std::memory_order_release);
}

static int
progressive_cut_clamp(
    const BObolLodProgressiveMeshGeneration *generation, int cut)
{
    if (!generation || generation->residentCut < 0)
	return -1;
    if (cut < generation->minimumCut)
	cut = generation->minimumCut;
    if (cut > generation->residentCut)
	cut = generation->residentCut;
    return cut;
}

static int
progressive_hierarchy_cut_clamp(
    const BObolLodProgressiveMeshGeneration *generation, int cut)
{
    if (!generation || generation->minimumCut < 0 ||
	generation->maximumCut < generation->minimumCut)
	return -1;
    return std::max(generation->minimumCut,
	std::min(generation->maximumCut, cut));
}

static int
progressive_cut_for_screen_error(
    const BObolLodProgressiveMeshGeneration *generation,
    double projectedPixelDiameter, double targetPixelError)
{
    if (!generation || generation->cuts.empty() ||
	!std::isfinite(projectedPixelDiameter) ||
	!std::isfinite(targetPixelError) || projectedPixelDiameter <= 0.0 ||
	targetPixelError <= 0.0)
	return -1;
    const SbVec3f extent = generation->quantizationMaximum -
	generation->quantizationMinimum;
    const double diagonal = std::sqrt(
	static_cast<double>(extent.sqrLength()));
    if (!(diagonal > 0.0))
	return generation->maximumCut;
    const double maximumObjectError =
	targetPixelError * diagonal / projectedPixelDiameter;
    int low = generation->minimumCut;
    int high = generation->maximumCut;
    while (low < high) {
	const int middle = low + (high - low) / 2;
	if (generation->cuts[static_cast<size_t>(middle)].object_error <=
		maximumObjectError)
	    high = middle;
	else
	    low = middle + 1;
    }
    return low;
}

static double
progressive_projected_error_at_cut(
    const BObolLodProgressiveMeshGeneration *generation, int cut,
    double projectedPixelDiameter)
{
    if (!generation || cut < 0 ||
	static_cast<size_t>(cut) >= generation->cuts.size() ||
	!std::isfinite(projectedPixelDiameter) || projectedPixelDiameter <= 0.0)
	return std::numeric_limits<double>::infinity();
    const SbVec3f extent = generation->quantizationMaximum -
	generation->quantizationMinimum;
    const double diagonal = std::sqrt(
	static_cast<double>(extent.sqrLength()));
    if (!(diagonal > 0.0))
	return 0.0;
    return generation->cuts[static_cast<size_t>(cut)].object_error *
	projectedPixelDiameter / diagonal;
}

static bool
bobol_prepare_authored_corner_normals(
    Obol::TriMesh &mesh, const std::vector<SbVec3f> &cornerNormals);

static bool
bobol_lod_finite_render_coordinate(fastf_t value)
{
    return std::isfinite(value) &&
	std::fabs(value) <= static_cast<fastf_t>(FLT_MAX);
}

static bool
bobol_lod_finite_vector(const fastf_t *value)
{
    return value &&
	bobol_lod_finite_render_coordinate(value[X]) &&
	bobol_lod_finite_render_coordinate(value[Y]) &&
	bobol_lod_finite_render_coordinate(value[Z]);
}

static bool
bobol_lod_valid_bounds(const fastf_t *minimum, const fastf_t *maximum)
{
    return bobol_lod_finite_vector(minimum) &&
	bobol_lod_finite_vector(maximum) &&
	minimum[X] <= maximum[X] &&
	minimum[Y] <= maximum[Y] &&
	minimum[Z] <= maximum[Z];
}

static bool
bobol_lod_point_in_domain(const fastf_t *point,
	const fastf_t *minimum, const fastf_t *maximum)
{
    if (!bobol_lod_finite_vector(point) ||
	!bobol_lod_valid_bounds(minimum, maximum))
	return false;
    for (size_t axis = 0; axis < 3; ++axis) {
	const fastf_t scale = std::max<fastf_t>(1.0,
	    std::max(std::fabs(minimum[axis]),
		std::max(std::fabs(maximum[axis]),
		    std::fabs(maximum[axis] - minimum[axis]))));
	const fastf_t tolerance = 64.0 *
	    std::numeric_limits<fastf_t>::epsilon() * scale;
	if (point[axis] < minimum[axis] - tolerance ||
	    point[axis] > maximum[axis] + tolerance)
	    return false;
    }
    return true;
}

static bool
bobol_lod_oriented_bounds_valid(
    const struct BObolMeshLodHierarchyInfo &hierarchy)
{
    return bobol_mesh_lod_oriented_bounds_validate(&hierarchy) != 0;
}

static std::optional<std::array<SbVec3f, 8>>
bobol_lod_oriented_bounds(
    const struct BObolMeshLodHierarchyInfo &hierarchy)
{
    if (hierarchy.oriented_bounds_valid != 1 ||
	!bobol_lod_oriented_bounds_valid(hierarchy))
	return std::nullopt;
    std::array<SbVec3f, 8> corners;
    for (size_t corner = 0; corner < corners.size(); ++corner)
	corners[corner].setValue(
	    static_cast<float>(hierarchy.oriented_bounds[corner][X]),
	    static_cast<float>(hierarchy.oriented_bounds[corner][Y]),
	    static_cast<float>(hierarchy.oriented_bounds[corner][Z]));

    /* Validate the corner representation after narrowing to renderer floats.
     * Do not use the geometry's axis-aligned bounds as oriented-volume
     * evidence: a useful rotated box normally excludes corners of its AABB.
     * Admission below checks the actual renderer positions and discards this
     * optional optimization if they do not fit. */
    Obol::PartGeometryBuilder proxy;
    proxy.aggregateProxyCorners = corners;
    if (!Obol::cadValidatePartGeometry(proxy))
	return std::nullopt;
    return corners;
}

static void
bobol_lod_set_oriented_bounds(
    struct BObolMeshLodHierarchyInfo &hierarchy,
    const std::optional<std::array<SbVec3f, 8>> &corners)
{
    hierarchy.oriented_bounds_valid = corners ? 1 : 0;
    if (!corners)
	return;
    for (size_t corner = 0; corner < corners->size(); ++corner)
	VSET(hierarchy.oriented_bounds[corner],
	    (*corners)[corner][0], (*corners)[corner][1],
	    (*corners)[corner][2]);
}

static bool
bobol_lod_hierarchy_valid(
    const struct BObolMeshLodData &data,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int residentCut, const char **failureReason = NULL)
{
    const auto fail = [failureReason](const char *reason) {
	if (failureReason)
	    *failureReason = reason;
	return false;
    };
    if (failureReason)
	*failureReason = NULL;
    if (hierarchy.min_cut < 0 ||
	hierarchy.max_cut < hierarchy.min_cut ||
	hierarchy.max_cut >= BOBOL_MESH_LOD_CUT_COUNT_MAX ||
	hierarchy.cut_count !=
	    static_cast<uint32_t>(hierarchy.max_cut + 1) ||
	residentCut < hierarchy.min_cut ||
	residentCut > hierarchy.max_cut ||
	hierarchy.resident_cut != residentCut)
	return fail("cut interval");
    if (!bobol_lod_valid_bounds(data.bmin, data.bmax))
	return fail("mesh bounds");
    if (!bobol_lod_valid_bounds(hierarchy.quantization_min,
	    hierarchy.quantization_max))
	return fail("quantization bounds");
    /* Oriented bounds are an optional proxy optimization.  They are
     * validated immediately before renderer publication, where invalid
     * metadata falls back to the authoritative AABB.  It must not reject an
     * otherwise valid progressive mesh hierarchy. */
    if ((hierarchy.cluster_grid_resolution != 0 ||
	 hierarchy.cluster_count != 0 || hierarchy.clusters != NULL) &&
	(hierarchy.cluster_grid_resolution !=
	     BOBOL_MESH_LOD_CLUSTER_GRID_RESOLUTION ||
	 !hierarchy.cluster_count ||
	 hierarchy.cluster_count > BOBOL_MESH_LOD_CLUSTER_COUNT ||
	 !hierarchy.clusters))
	return fail("cluster layout");

    uint64_t priorPoints = 0;
    uint64_t priorFaces = 0;
    double priorError = std::numeric_limits<double>::infinity();
    for (int cut = 0; cut <= hierarchy.max_cut; ++cut) {
	const uint64_t points = hierarchy.cuts[cut].point_count;
	const uint64_t faces = hierarchy.cuts[cut].face_count;
	if (points < priorPoints || faces < priorFaces ||
	    !std::isfinite(hierarchy.cuts[cut].object_error) ||
	    hierarchy.cuts[cut].object_error < 0.0 ||
	    hierarchy.cuts[cut].object_error > priorError ||
	    points > static_cast<size_t>(
		std::numeric_limits<int32_t>::max()) ||
	    faces > static_cast<size_t>(
		std::numeric_limits<int32_t>::max()) / 3)
	    return fail("cut record");
	priorPoints = points;
	priorFaces = faces;
	priorError = hierarchy.cuts[cut].object_error;
	for (int axis = 0; axis < 3; ++axis) {
	    const uint8_t bits =
		hierarchy.cuts[cut].quantization_bits[axis];
	    if (bits < 1 || bits > BOBOL_MESH_LOD_QUANTIZATION_BITS ||
		(cut > 0 && bits <
		 hierarchy.cuts[cut - 1].quantization_bits[axis]))
		return fail("quantization schedule");
	}
    }
    if (!hierarchy.cuts[hierarchy.max_cut].exact ||
	hierarchy.cuts[hierarchy.max_cut].object_error > 0.0)
	return fail("terminal cut");
    uint64_t clusteredIndices = 0;
    for (uint32_t cluster = 0; cluster < hierarchy.cluster_count;
	 cluster++) {
	const BObolMeshLodClusterInfo &info = hierarchy.clusters[cluster];
	if (info.cluster_id >= BOBOL_MESH_LOD_CLUSTER_COUNT ||
	    (cluster && info.cluster_id <=
	     hierarchy.clusters[cluster - 1].cluster_id) ||
	    !info.range_count || !info.ranges ||
	    !bobol_lod_valid_bounds(info.bmin, info.bmax))
	    return fail("cluster record");
	for (uint32_t rangeIndex = 0; rangeIndex < info.range_count;
	     ++rangeIndex) {
	    const BObolMeshLodClusterRange &range = info.ranges[rangeIndex];
	    if (range.activation_cut >= hierarchy.cut_count ||
		range.first_index % 3 || range.index_count < 3 ||
		range.index_count % 3 ||
		static_cast<uint64_t>(range.first_index) +
		    range.index_count >
		    hierarchy.cuts[range.activation_cut].face_count * 3)
		return fail("cluster range");
	    clusteredIndices += range.index_count;
	}
    }
    if (hierarchy.cluster_count && clusteredIndices !=
	hierarchy.cuts[hierarchy.max_cut].face_count * 3)
	return fail("cluster coverage");
    if (hierarchy.cuts[residentCut].point_count != data.point_orig_count ||
	hierarchy.cuts[residentCut].face_count != data.face_count ||
	!data.point_orig_count || !data.face_count)
	return fail("resident population");
    return true;
}

static std::shared_ptr<BObolLodProgressiveMeshGeneration>
progressive_generation_from_data(
    const struct BObolMeshLodData &data,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int residentCut, SbBool shadedCullBackfaces,
    uint64_t revision, uint64_t progressiveLineage,
    const char **failureReason)
{
    const auto fail = [failureReason](const char *reason) {
	if (failureReason)
	    *failureReason = reason;
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    };
    if (failureReason)
	*failureReason = NULL;
    if (!bobol_lod_hierarchy_valid(
	    data, hierarchy, residentCut, failureReason))
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    std::shared_ptr<BObolLodProgressiveMeshGeneration> generation(
	new BObolLodProgressiveMeshGeneration);
    std::shared_ptr<Obol::PartGeometryBuilder> geometry(
	new Obol::PartGeometryBuilder);
    Obol::TriMesh mesh;
    generation->cuts.assign(hierarchy.cuts,
	hierarchy.cuts + hierarchy.cut_count);
    mesh.progressiveCuts.resize(hierarchy.cut_count);
    for (uint32_t cut = 0; cut < hierarchy.cut_count; ++cut) {
	generation->pointCount[cut] = static_cast<size_t>(
	    hierarchy.cuts[cut].point_count);
	generation->faceCount[cut] = static_cast<size_t>(
	    hierarchy.cuts[cut].face_count);
	mesh.progressiveCuts[cut].quantization = {
	    hierarchy.cuts[cut].quantization_bits[X],
	    hierarchy.cuts[cut].quantization_bits[Y],
	    hierarchy.cuts[cut].quantization_bits[Z]
	};
    }
    generation->minimumCut = hierarchy.min_cut;
    generation->maximumCut = hierarchy.max_cut;
    generation->residentCut = residentCut;
    generation->bounds = SbBox3f(
	SbVec3f(static_cast<float>(data.bmin[X]),
		static_cast<float>(data.bmin[Y]),
		static_cast<float>(data.bmin[Z])),
	SbVec3f(static_cast<float>(data.bmax[X]),
		static_cast<float>(data.bmax[Y]),
		static_cast<float>(data.bmax[Z])));
    generation->quantizationMinimum.setValue(
	static_cast<float>(hierarchy.quantization_min[X]),
	static_cast<float>(hierarchy.quantization_min[Y]),
	static_cast<float>(hierarchy.quantization_min[Z]));
    generation->quantizationMaximum.setValue(
	static_cast<float>(hierarchy.quantization_max[X]),
	static_cast<float>(hierarchy.quantization_max[Y]),
	static_cast<float>(hierarchy.quantization_max[Z]));
    generation->aggregateProxyCorners = bobol_lod_oriented_bounds(hierarchy);
    generation->shadedCullBackfaces = shadedCullBackfaces;
    generation->revision = revision ? revision : 1;

    mesh.positions.resize(data.point_orig_count);
    for (size_t i = 0; i < data.point_orig_count; ++i) {
	if (!bobol_lod_point_in_domain(data.points_orig[i],
		hierarchy.quantization_min,
		hierarchy.quantization_max))
	    return fail("point outside quantization domain");
	mesh.positions[i] = SbVec3f(
	    static_cast<float>(data.points_orig[i][X]),
	    static_cast<float>(data.points_orig[i][Y]),
	    static_cast<float>(data.points_orig[i][Z]));
    }

    if (data.face_count >
	    std::numeric_limits<size_t>::max() / 3)
	return fail("face count overflow");
    const size_t indexCount = data.face_count * 3;
    std::array<size_t, BOBOL_MESH_LOD_CUT_COUNT_MAX> cutIndexCount = {};
    for (uint32_t cut = 0; cut < hierarchy.cut_count; ++cut) {
	const size_t count = generation->faceCount[cut] >
		std::numeric_limits<size_t>::max() / 3 ? indexCount :
	    generation->faceCount[cut] * 3;
	cutIndexCount[cut] = std::min(count, indexCount);
	mesh.progressiveCuts[cut].indexCount = static_cast<uint32_t>(
	    std::min(cutIndexCount[cut],
		static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
    }
    mesh.indices.resize(indexCount);
    uint32_t maximumIndex = 0;
    int completedCut = 0;
    while (completedCut < static_cast<int>(hierarchy.cut_count) &&
	cutIndexCount[completedCut] == 0) {
	mesh.progressiveCuts[completedCut].positionCount = 0;
	++completedCut;
    }
    for (size_t i = 0; i < indexCount; ++i) {
	if (static_cast<size_t>(data.faces[i]) >= mesh.positions.size())
	    return fail("vertex index");
	mesh.indices[i] = static_cast<uint32_t>(data.faces[i]);
	maximumIndex = std::max(maximumIndex, mesh.indices[i]);
	while (completedCut < static_cast<int>(hierarchy.cut_count) &&
	    cutIndexCount[completedCut] <= i + 1) {
	    mesh.progressiveCuts[completedCut].positionCount =
		maximumIndex + 1;
	    ++completedCut;
	}
    }
    if (data.normals) {
	std::vector<SbVec3f> cornerNormals;
	cornerNormals.resize(indexCount);
	for (size_t i = 0; i < indexCount; ++i) {
	    if (!bobol_lod_finite_vector(data.normals[i]))
		return fail("non-finite normal");
	    cornerNormals[i] = SbVec3f(
		static_cast<float>(data.normals[i][X]),
		static_cast<float>(data.normals[i][Y]),
		static_cast<float>(data.normals[i][Z]));
	}
	if (!bobol_prepare_authored_corner_normals(
		mesh, cornerNormals))
	    return fail("corner-normal preparation");
    }

    if (mesh.positions.empty() || mesh.indices.size() < 3 ||
	mesh.indices.size() % 3 != 0 ||
	(!mesh.normals.empty() &&
	 mesh.normals.size() != mesh.positions.size()))
	return fail("renderer mesh arrays");
    mesh.bounds = generation->bounds;
    mesh.progressiveMinimumCut =
	static_cast<uint8_t>(std::max(0, generation->minimumCut));
    int drawableCut = generation->residentCut;
    for (int cut = generation->residentCut + 1;
	 cut <= generation->maximumCut &&
	 cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut) {
	if (generation->faceCount[cut] >
		std::numeric_limits<size_t>::max() / 3 ||
	    generation->faceCount[cut] * 3 > mesh.indices.size() ||
	    generation->pointCount[cut] > data.point_orig_count)
	    break;
	drawableCut = cut;
    }
    mesh.progressiveResidentCut =
	static_cast<uint8_t>(std::max(0, drawableCut));
    mesh.progressiveQuantizationMinimum =
	generation->quantizationMinimum;
    mesh.progressiveQuantizationMaximum =
	generation->quantizationMaximum;
    /* The cache publishes only occupied cells, but these ranges still share
     * one global uniform-grid prefix and one resident cut.  Preserve that
     * distinction from independently resident adaptive chunk pages: treating
     * sparse uniform cells as adaptive changes controller demand semantics
     * and can multiply refinement work in a many-leaf scene. */
    mesh.progressiveClusterGridResolution =
	hierarchy.cluster_grid_resolution;
    mesh.progressiveClusters.resize(hierarchy.cluster_count);
    for (uint32_t cluster = 0; cluster < hierarchy.cluster_count;
	 cluster++) {
	const BObolMeshLodClusterInfo &source = hierarchy.clusters[cluster];
	Obol::ProgressiveTriangleCluster &target =
	    mesh.progressiveClusters[cluster];
	target.residentCut = mesh.progressiveResidentCut;
	target.bounds = SbBox3f(
	    SbVec3f(static_cast<float>(source.bmin[X]),
		static_cast<float>(source.bmin[Y]),
		static_cast<float>(source.bmin[Z])),
	    SbVec3f(static_cast<float>(source.bmax[X]),
		static_cast<float>(source.bmax[Y]),
		static_cast<float>(source.bmax[Z])));
	target.ranges.reserve(source.range_count);
	for (uint32_t range = 0; range < source.range_count; ++range) {
	    target.ranges.push_back({
		source.ranges[range].first_index,
		source.ranges[range].index_count,
		source.ranges[range].activation_cut
	    });
	}
    }
    /* Authored corner-normal canonicalization may globally split and reorder
     * vertices when a richer prefix arrives.  It cannot certify append-only
     * identity yet, so only the ordinary coordinate/index PoP path exposes
     * the lineage optimization to Obol. */
    mesh.progressiveLineage = data.normals ? 0 : progressiveLineage;
    if (data.normals) {
	/* Corner-normal preparation may duplicate vertices and rewrite indices.
	 * Recompute only that uncommon authored-normal case; ordinary PoP assets
	 * populated these cumulative maxima during their single validation/copy
	 * pass above. */
	size_t scannedIndices = 0;
	maximumIndex = 0;
	for (uint32_t cut = 0; cut < hierarchy.cut_count; ++cut) {
	    for (; scannedIndices < cutIndexCount[cut]; ++scannedIndices)
		maximumIndex = std::max(
		    maximumIndex, mesh.indices[scannedIndices]);
	    const size_t positionCount = cutIndexCount[cut] ?
		static_cast<size_t>(maximumIndex) + 1 : 0;
	    mesh.progressiveCuts[cut].positionCount = static_cast<uint32_t>(
		std::min(positionCount,
		    static_cast<size_t>(
			std::numeric_limits<uint32_t>::max())));
	}
    }
    geometry->shaded = std::move(mesh);
    geometry->shadedCullBackfaces =
	shadedCullBackfaces ? true : false;
    geometry->subpixelProxyEligible = true;
    geometry->aggregateProxyCorners = generation->aggregateProxyCorners;
    const std::shared_ptr<const Obol::PartGeometry> completed =
	bobol_cad_build_geometry_with_optional_proxy(
	    std::move(*geometry), "PoP generation");
    if (!completed)
	return fail("CAD geometry admission");
    generation->shadedGeometry = completed;
    return generation;
}

struct BObolChunkReadContext {
    const BObolMeshLodHierarchyInfo *hierarchy = NULL;
    std::vector<std::shared_ptr<const
	BObolLodProgressiveMeshGeneration::ResidentChunk>> *chunks = NULL;
};

static constexpr size_t lod_spatial_parallel_page_threshold = 32;
static constexpr size_t lod_spatial_parallel_worker_limit = 4;

static size_t
lod_spatial_worker_count(size_t pageCount)
{
    if (pageCount < lod_spatial_parallel_page_threshold)
	return 1;
    const int available = bu_avail_cpus();
    return std::min(pageCount, std::min(
	lod_spatial_parallel_worker_limit,
	static_cast<size_t>(std::max(1, available))));
}

static int
bobol_progressive_chunk_read(uint32_t chunkId, int cut,
    const point_t *points, size_t pointCount,
    const uint32_t *faces, size_t faceCount,
    const vect_t *normals, size_t normalCount, void *callbackData)
{
    BObolChunkReadContext *context =
	static_cast<BObolChunkReadContext *>(callbackData);
    if (!context || !context->hierarchy || !context->chunks ||
	chunkId >= context->hierarchy->chunk_count ||
	chunkId >= context->chunks->size())
	return 0;
    const BObolMeshLodChunkInfo &info =
	context->hierarchy->chunks[chunkId];
    if (cut < info.min_cut) {
	(*context->chunks)[chunkId].reset();
	return pointCount == 0 && faceCount == 0 && normalCount == 0;
    }
    if (!points || !faces || !pointCount || !faceCount ||
	faceCount > SIZE_MAX / 3 ||
	normalCount != (context->hierarchy->has_normals ? faceCount * 3 : 0) ||
	(normalCount && !normals))
	return 0;
    std::shared_ptr<BObolLodProgressiveMeshGeneration::ResidentChunk> chunk(
	new BObolLodProgressiveMeshGeneration::ResidentChunk);
    chunk->chunkId = chunkId;
    chunk->residentCut = cut;
    chunk->geometryRevision = bobol_next_progressive_lineage();
    chunk->bounds = SbBox3f(
	SbVec3f(static_cast<float>(info.bmin[X]),
	    static_cast<float>(info.bmin[Y]),
	    static_cast<float>(info.bmin[Z])),
	SbVec3f(static_cast<float>(info.bmax[X]),
	    static_cast<float>(info.bmax[Y]),
	    static_cast<float>(info.bmax[Z])));
    chunk->positions.resize(pointCount);
    for (size_t point = 0; point < pointCount; ++point) {
	if (!bobol_lod_finite_vector(points[point]))
	    return 0;
	chunk->positions[point].setValue(
	    static_cast<float>(points[point][X]),
	    static_cast<float>(points[point][Y]),
	    static_cast<float>(points[point][Z]));
    }
    chunk->indices.assign(faces, faces + faceCount * 3);
    for (uint32_t index : chunk->indices)
	if (index >= pointCount)
	    return 0;
    if (normalCount) {
	chunk->cornerNormals.resize(normalCount);
	for (size_t normal = 0; normal < normalCount; ++normal) {
	    if (!bobol_lod_finite_vector(normals[normal]))
		return 0;
	    chunk->cornerNormals[normal].setValue(
		static_cast<float>(normals[normal][X]),
		static_cast<float>(normals[normal][Y]),
		static_cast<float>(normals[normal][Z]));
	}
    }
    (*context->chunks)[chunkId] = chunk;
    return 1;
}

static std::shared_ptr<BObolLodProgressiveMeshGeneration>
progressive_generation_from_chunks(
    const BObolMeshLodHierarchyInfo &hierarchy,
    const std::vector<std::shared_ptr<const
	BObolLodProgressiveMeshGeneration::ResidentChunk>> &residentChunks,
    SbBool shadedCullBackfaces, uint64_t revision,
    uint64_t progressiveLineage)
{
    if (!hierarchy.chunks || !hierarchy.chunk_count ||
	residentChunks.size() != hierarchy.chunk_count ||
	hierarchy.cut_count == 0 ||
	hierarchy.cut_count > BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    std::shared_ptr<BObolLodProgressiveMeshGeneration> generation(
	new BObolLodProgressiveMeshGeneration);
    generation->minimumCut = hierarchy.min_cut;
    generation->maximumCut = hierarchy.max_cut;
    generation->revision = revision ? revision : 1;
    generation->shadedCullBackfaces = shadedCullBackfaces;
    generation->quantizationMinimum.setValue(
	static_cast<float>(hierarchy.quantization_min[X]),
	static_cast<float>(hierarchy.quantization_min[Y]),
	static_cast<float>(hierarchy.quantization_min[Z]));
    generation->quantizationMaximum.setValue(
	static_cast<float>(hierarchy.quantization_max[X]),
	static_cast<float>(hierarchy.quantization_max[Y]),
	static_cast<float>(hierarchy.quantization_max[Z]));
    generation->bounds = SbBox3f(
	generation->quantizationMinimum, generation->quantizationMaximum);
    generation->aggregateProxyCorners = bobol_lod_oriented_bounds(hierarchy);
    generation->cuts.assign(hierarchy.cuts,
	hierarchy.cuts + hierarchy.cut_count);
    for (uint32_t cut = 0; cut < hierarchy.cut_count; ++cut) {
	generation->pointCount[cut] = static_cast<size_t>(
	    hierarchy.cuts[cut].point_count);
	generation->faceCount[cut] = static_cast<size_t>(
	    hierarchy.cuts[cut].face_count);
    }
    generation->chunks.resize(hierarchy.chunk_count);
    generation->residentChunks = residentChunks;
    for (uint32_t chunk = 0; chunk < hierarchy.chunk_count; ++chunk) {
	const BObolMeshLodChunkInfo &source = hierarchy.chunks[chunk];
	auto &target = generation->chunks[chunk];
	if (source.chunk_id != chunk || source.min_cut < hierarchy.min_cut ||
	    source.max_cut != hierarchy.max_cut)
	    return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	target.chunkId = chunk;
	target.minimumCut = source.min_cut;
	target.maximumCut = source.max_cut;
	target.bounds = SbBox3f(
	    SbVec3f(static_cast<float>(source.bmin[X]),
		static_cast<float>(source.bmin[Y]),
		static_cast<float>(source.bmin[Z])),
	    SbVec3f(static_cast<float>(source.bmax[X]),
		static_cast<float>(source.bmax[Y]),
		static_cast<float>(source.bmax[Z])));
	std::copy(source.cuts,
	    source.cuts + BOBOL_MESH_LOD_CUT_COUNT_MAX,
	    target.cuts.begin());
    }

    /* This scalar is the asset high-water mark, not a promise that every
     * private page is resident to the same cut.  Occurrence/view consumers
     * must use residentCutForChunks() for their demanded page set.  Using the
     * poorest common frontier here stranded a rich visible page whenever an
     * unrelated retained page was deliberately coarser. */
    generation->residentCut = hierarchy.min_cut - 1;
    for (const auto &resident : residentChunks)
	if (resident)
	    generation->residentCut = std::max(
		generation->residentCut, resident->residentCut);

    Obol::TriMesh mesh;
    mesh.bounds = generation->bounds;
    mesh.progressiveMinimumCut = static_cast<uint8_t>(
	std::max(0, generation->minimumCut));
    int adaptiveDrawableCut = generation->minimumCut - 1;
    bool haveResidentChunk = false;
    for (uint32_t chunk = 0; chunk < hierarchy.chunk_count; ++chunk) {
	const auto &resident = residentChunks[chunk];
	if (!resident)
	    continue;
	haveResidentChunk = true;
	/* residentCut records the last population-bearing cut actually read for
	 * this page.  The same arrays may draw later coordinate-only cuts: PoP
	 * quantization refines the resident coordinates without appending another
	 * vertex or triangle.  Advertise that drawable frontier to Obol or a page
	 * whose topology completed early is needlessly rebuilt (and visually
	 * remains snapped coarse) whenever the view asks for finer coordinates. */
	int drawableCut = resident->residentCut;
	for (int candidate = resident->residentCut + 1;
		candidate <= hierarchy.max_cut; ++candidate) {
	    const BObolMeshLodChunkCutInfo &candidateInfo =
		hierarchy.chunks[chunk].cuts[candidate];
	    if (candidateInfo.point_count > resident->positions.size() ||
		static_cast<uint64_t>(candidateInfo.face_count) * 3u >
		    resident->indices.size())
		break;
	    drawableCut = candidate;
	}
	adaptiveDrawableCut = std::max(adaptiveDrawableCut, drawableCut);
    }
    if (!haveResidentChunk)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();

    /*
     * A spatial hierarchy is already partitioned into independently resident
     * renderer pages.  Retain those pages as the generation itself.  Packing
     * them into a second monolithic vector made every view-local refinement
     * copy all previously resident pages and turned close Lucy zooms into a
     * multi-second operation.  The one-page case below remains the ordinary
     * cumulative-prefix representation used by small/direct assets.
     */
    if (hierarchy.chunk_count > 1)
	return generation;

    const auto &resident = residentChunks[0];
    if (!resident || resident->positions.empty() ||
	resident->indices.size() < 3 || resident->indices.size() % 3 != 0)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();

    std::shared_ptr<Obol::PartGeometryBuilder> geometry(new Obol::PartGeometryBuilder);
    mesh.positions = resident->positions;
    mesh.indices = resident->indices;
    if (!resident->cornerNormals.empty() &&
	!bobol_prepare_authored_corner_normals(
	    mesh, resident->cornerNormals))
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    mesh.progressiveResidentCut = static_cast<uint8_t>(
	std::max(0, adaptiveDrawableCut));
    mesh.progressiveQuantizationMinimum = generation->quantizationMinimum;
    mesh.progressiveQuantizationMaximum = generation->quantizationMaximum;
    mesh.progressiveLineage = progressiveLineage;
    mesh.progressiveCuts.resize(hierarchy.cut_count);
    for (uint32_t cut = 0; cut < hierarchy.cut_count; ++cut) {
	mesh.progressiveCuts[cut].indexCount = static_cast<uint32_t>(
	    std::min<uint64_t>(
		static_cast<uint64_t>(
		    hierarchy.chunks[0].cuts[cut].face_count) * 3u,
		mesh.indices.size()));
	mesh.progressiveCuts[cut].positionCount = static_cast<uint32_t>(
	    std::min<uint64_t>(hierarchy.chunks[0].cuts[cut].point_count,
		mesh.positions.size()));
	mesh.progressiveCuts[cut].quantization = {
	    hierarchy.cuts[cut].quantization_bits[X],
	    hierarchy.cuts[cut].quantization_bits[Y],
	    hierarchy.cuts[cut].quantization_bits[Z]
	};
    }
    if (!mesh.normals.empty()) {
	/* Normal canonicalization can split a source position.  Triangle order is
	 * unchanged, so recompute the live vertex frontier for every cut. */
	size_t scanned = 0;
	uint32_t maximum = 0;
	for (Obol::ProgressiveTriangleCut &cut : mesh.progressiveCuts) {
	    const size_t end = std::min<size_t>(
		cut.indexCount, mesh.indices.size());
	    for (; scanned < end; ++scanned)
		maximum = std::max(maximum, mesh.indices[scanned]);
	    cut.positionCount = end ? maximum + 1u : 0u;
	}
    }
    geometry->shaded = std::move(mesh);
    if (!generation->aggregateProxyCorners)
	geometry->conservativeBounds = generation->bounds;
    geometry->shadedCullBackfaces = shadedCullBackfaces ? true : false;
    geometry->subpixelProxyEligible = true;
    geometry->aggregateProxyCorners = generation->aggregateProxyCorners;
    const std::shared_ptr<const Obol::PartGeometry> completed =
	bobol_cad_build_geometry_with_optional_proxy(
	    std::move(*geometry), "spatial PoP generation");
    if (!completed)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    generation->shadedGeometry = completed;
    return generation;
}

static std::shared_ptr<BObolLodProgressiveMeshGeneration>
progressive_generation_prefix(
    const BObolLodProgressiveMeshGeneration &source, int residentCut,
    const std::vector<BObolLodChunkCut> *retainedChunkCuts = NULL)
{
    const int cut = progressive_cut_clamp(&source, residentCut);
    if (cut < 0)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    if (!source.chunks.empty()) {
	if (retainedChunkCuts) {
	    if (retainedChunkCuts->empty())
		return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	    for (size_t i = 0; i < retainedChunkCuts->size(); ++i)
		if ((*retainedChunkCuts)[i].chunkId >= source.chunks.size() ||
		    (*retainedChunkCuts)[i].cut < source.minimumCut ||
		    (*retainedChunkCuts)[i].cut > cut ||
		    (i && (*retainedChunkCuts)[i].chunkId <=
			(*retainedChunkCuts)[i - 1].chunkId))
		    return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	}
	std::vector<BObolMeshLodChunkInfo> chunkInfos(source.chunks.size());
	std::vector<std::shared_ptr<const
	    BObolLodProgressiveMeshGeneration::ResidentChunk>> resident(
	    source.chunks.size());
	for (size_t chunkId = 0; chunkId < source.chunks.size(); ++chunkId) {
	    const auto &metadata = source.chunks[chunkId];
	    BObolMeshLodChunkInfo &info = chunkInfos[chunkId];
	    memset(&info, 0, sizeof(info));
	    info.chunk_id = metadata.chunkId;
	    info.min_cut = metadata.minimumCut;
	    info.max_cut = metadata.maximumCut;
	    const SbVec3f minimum = metadata.bounds.getMin();
	    const SbVec3f maximum = metadata.bounds.getMax();
	    VSET(info.bmin, minimum[0], minimum[1], minimum[2]);
	    VSET(info.bmax, maximum[0], maximum[1], maximum[2]);
	    std::copy(metadata.cuts.begin(), metadata.cuts.end(), info.cuts);
	    int chunkCut = cut;
	    if (retainedChunkCuts) {
		const auto selected = std::lower_bound(
		    retainedChunkCuts->begin(), retainedChunkCuts->end(),
		    static_cast<uint32_t>(chunkId),
		    [](const BObolLodChunkCut &entry, uint32_t id) {
			return entry.chunkId < id;
		    });
		if (selected == retainedChunkCuts->end() ||
		    selected->chunkId != chunkId)
		    continue;
		chunkCut = selected->cut;
	    }
	    if (!metadata.cuts[chunkCut].face_count)
		continue;
	    if (chunkId >= source.residentChunks.size() ||
		!source.residentChunks[chunkId])
		return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	    const auto &prior = source.residentChunks[chunkId];
	    const size_t pointCount = metadata.cuts[chunkCut].point_count;
	    const size_t indexCount =
		static_cast<size_t>(metadata.cuts[chunkCut].face_count) * 3;
	    if (pointCount > prior->positions.size() ||
		indexCount > prior->indices.size() ||
		(!prior->cornerNormals.empty() &&
		 indexCount > prior->cornerNormals.size()))
		return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	    std::shared_ptr<BObolLodProgressiveMeshGeneration::ResidentChunk>
		trimmed(new BObolLodProgressiveMeshGeneration::ResidentChunk);
	    trimmed->chunkId = static_cast<uint32_t>(chunkId);
	    trimmed->residentCut = chunkCut;
	    trimmed->geometryRevision = bobol_next_progressive_lineage();
	    trimmed->bounds = prior->bounds;
	    trimmed->positions.assign(prior->positions.begin(),
		prior->positions.begin() + pointCount);
	    trimmed->indices.assign(prior->indices.begin(),
		prior->indices.begin() + indexCount);
	    if (!prior->cornerNormals.empty())
		trimmed->cornerNormals.assign(prior->cornerNormals.begin(),
		    prior->cornerNormals.begin() + indexCount);
	    resident[chunkId] = trimmed;
	}
	BObolMeshLodHierarchyInfo hierarchy =
	    BOBOL_MESH_LOD_HIERARCHY_INFO_INIT;
	hierarchy.min_cut = source.minimumCut;
	hierarchy.max_cut = source.maximumCut;
	hierarchy.resident_cut = cut;
	hierarchy.cut_count = static_cast<uint32_t>(source.cuts.size());
	hierarchy.has_normals = source.shadedGeometry &&
	    source.shadedGeometry->shaded &&
	    !source.shadedGeometry->shaded->normals.empty();
	hierarchy.shaded_cull_backfaces = source.shadedCullBackfaces;
	VSET(hierarchy.quantization_min,
	    source.quantizationMinimum[0], source.quantizationMinimum[1],
	    source.quantizationMinimum[2]);
	VSET(hierarchy.quantization_max,
	    source.quantizationMaximum[0], source.quantizationMaximum[1],
	    source.quantizationMaximum[2]);
	hierarchy.chunk_count = static_cast<uint32_t>(chunkInfos.size());
	hierarchy.chunks = chunkInfos.data();
	std::copy(source.cuts.begin(), source.cuts.end(), hierarchy.cuts);
	bobol_lod_set_oriented_bounds(hierarchy, source.aggregateProxyCorners);
	return progressive_generation_from_chunks(hierarchy, resident,
	    source.shadedCullBackfaces,
	    bobol_identity_successor_or_terminate(source.revision), 0);
    }
    const size_t indexCount = source.faceCount[cut] * 3;
    const Obol::TriMesh *sourceMesh =
	source.shadedGeometry && source.shadedGeometry->shaded ?
	    &*source.shadedGeometry->shaded : NULL;
    if (!sourceMesh || indexCount > sourceMesh->indices.size())
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    const size_t positionCount =
	sourceMesh->positionCountAtCut(static_cast<uint8_t>(cut));
    if (positionCount > sourceMesh->positions.size() ||
	(!sourceMesh->normals.empty() &&
	 positionCount > sourceMesh->normals.size()))
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();

    std::shared_ptr<BObolLodProgressiveMeshGeneration> generation(
	new BObolLodProgressiveMeshGeneration);
    generation->pointCount = source.pointCount;
    generation->faceCount = source.faceCount;
    generation->cuts = source.cuts;
    generation->bounds = source.bounds;
    generation->quantizationMinimum = source.quantizationMinimum;
    generation->quantizationMaximum = source.quantizationMaximum;
    generation->aggregateProxyCorners = source.aggregateProxyCorners;
    generation->minimumCut = source.minimumCut;
    generation->maximumCut = source.maximumCut;
    generation->residentCut = cut;
    generation->revision =
	bobol_identity_successor_or_terminate(source.revision);
    generation->shadedCullBackfaces = source.shadedCullBackfaces;
    std::shared_ptr<Obol::PartGeometryBuilder> geometry(
	new Obol::PartGeometryBuilder);
    Obol::TriMesh mesh;
    mesh.bounds = sourceMesh->bounds;
    mesh.progressiveCuts = sourceMesh->progressiveCuts;
    mesh.progressiveMinimumCut = sourceMesh->progressiveMinimumCut;
    mesh.progressiveResidentCut = sourceMesh->progressiveResidentCut;
    mesh.progressiveQuantizationMinimum =
	sourceMesh->progressiveQuantizationMinimum;
    mesh.progressiveQuantizationMaximum =
	sourceMesh->progressiveQuantizationMaximum;
    mesh.progressiveLineage = sourceMesh->progressiveLineage;
    mesh.progressiveClusters = sourceMesh->progressiveClusters;
    mesh.progressiveClusterGridResolution =
	sourceMesh->progressiveClusterGridResolution;
    /* Construct directly at the retained size.  Copying the rich vector and
     * then assigning a short range preserves its old capacity, which made a
     * successful quiet compaction report (and actually retain) almost all of
     * the memory it was supposed to release. */
    mesh.positions.assign(
	sourceMesh->positions.begin(),
	sourceMesh->positions.begin() + positionCount);
    if (!sourceMesh->normals.empty())
	mesh.normals.assign(
	    sourceMesh->normals.begin(),
	    sourceMesh->normals.begin() + positionCount);
    mesh.indices.assign(
	sourceMesh->indices.begin(),
	sourceMesh->indices.begin() + indexCount);
    /*
     * residentCut is the persistent-array high-water mark being retained.
     * It is not necessarily the richest cut those arrays can draw: later
     * producer cuts may add only quantization precision.  Preserve that
     * coordinate-only frontier after compaction or a small mesh whose full
     * topology arrived early becomes permanently stuck on the coarse snap
     * associated with its last population-bearing cut.
     */
    int drawableCut = cut;
    for (int candidate = cut + 1;
	 candidate <= source.maximumCut &&
	 static_cast<size_t>(candidate) < sourceMesh->progressiveCuts.size();
	 ++candidate) {
	const Obol::ProgressiveTriangleCut &candidateInfo =
	    sourceMesh->progressiveCuts[static_cast<size_t>(candidate)];
	/* The public AtCut accessors clamp to progressiveResidentCut.  That is
	 * correct for drawing, but not for deciding whether this newly trimmed
	 * snapshot may advertise a later coordinate-only cut: after a second
	 * compaction, clamping makes an actually richer nonresident cut look equal
	 * to the current one.  Inspect the immutable hierarchy metadata directly. */
	if (candidateInfo.positionCount > mesh.positions.size() ||
	    candidateInfo.indexCount > mesh.indices.size())
	    break;
	drawableCut = candidate;
    }
    mesh.progressiveResidentCut =
        static_cast<uint8_t>(std::max(0, drawableCut));
    geometry->shaded = std::move(mesh);
    geometry->shadedCullBackfaces =
	source.shadedCullBackfaces ? true : false;
    geometry->subpixelProxyEligible = true;
    geometry->aggregateProxyCorners = generation->aggregateProxyCorners;
    const std::shared_ptr<const Obol::PartGeometry> completed =
	bobol_cad_build_geometry_with_optional_proxy(
	    std::move(*geometry), "trimmed PoP generation");
    if (!completed)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    generation->shadedGeometry = completed;
    return generation;
}

struct BObolLodProgressiveMeshTrim {
    const BObolLodProgressiveMeshPrivate *owner = NULL;
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> source;
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation;
};

BObolLodProgressiveMeshSnapshot::BObolLodProgressiveMeshSnapshot(void) =
    default;

BObolLodProgressiveMeshSnapshot::BObolLodProgressiveMeshSnapshot(
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> &snapshot) :
    generation(snapshot)
{
}

SbBool
BObolLodProgressiveMeshSnapshot::isValid(void) const
{
    return this->generation && this->generation->minimumCut >= 0 &&
	this->generation->maximumCut >= this->generation->minimumCut ?
	TRUE : FALSE;
}

uint64_t
BObolLodProgressiveMeshSnapshot::revision(void) const
{
    return this->generation ? this->generation->revision : 0;
}

int
BObolLodProgressiveMeshSnapshot::minimumCut(void) const
{
    return this->generation ? this->generation->minimumCut : -1;
}

int
BObolLodProgressiveMeshSnapshot::maximumCut(void) const
{
    return this->generation ? this->generation->maximumCut : -1;
}

BObolLodCounts
BObolLodProgressiveMeshSnapshot::hierarchyCountsAtCut(
    int requestedCut, SbBool hasNormals) const
{
    BObolLodCounts counts;
    const int cut = progressive_hierarchy_cut_clamp(
	this->generation.get(), requestedCut);
    if (cut < 0)
	return counts;
    counts.faceCount = this->generation->faceCount[cut];
    counts.pointCount = this->generation->pointCount[cut];
    counts.originalPointCount = counts.pointCount;
    counts.normalCount = hasNormals ?
	(counts.faceCount > UINT64_MAX / 3 ?
	    UINT64_MAX : counts.faceCount * 3) : 0;
    return counts;
}

int
BObolLodProgressiveMeshSnapshot::cutForScreenError(
    double projectedPixelDiameter, double targetPixelError) const
{
    return progressive_cut_for_screen_error(this->generation.get(),
	projectedPixelDiameter, targetPixelError);
}

double
BObolLodProgressiveMeshSnapshot::projectedErrorAtCut(
    int cut, double projectedPixelDiameter) const
{
    return progressive_projected_error_at_cut(this->generation.get(), cut,
	projectedPixelDiameter);
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
    int residentCut,
    SbBool shadedCullBackfaces)
{
    const auto reject = [](const char *reason) {
	if (getenv("BOBOL_LOD_TRACE_REALIZATION"))
	    bu_log("BObol PoP realization rejected: %s\n", reason);
	return FALSE;
    };
    if (!this->p || !data.faces || !data.points_orig ||
	data.face_count == 0 || data.point_orig_count == 0 ||
	residentCut < hierarchy.min_cut ||
	residentCut > hierarchy.max_cut ||
	residentCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return reject("arguments");

    const char *hierarchyFailure = NULL;
    if (!bobol_lod_hierarchy_valid(
	    data, hierarchy, residentCut, &hierarchyFailure))
	return reject(hierarchyFailure ? hierarchyFailure : "hierarchy");

    if (data.point_orig_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	data.face_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) / 3)
	return reject("renderer address range");
    const size_t indexCount = data.face_count * 3;
    if (data.normals && data.normal_count != indexCount)
	return reject("normal count");

    std::lock_guard<std::mutex> updateLock(this->p->updateMutex);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> prior =
	progressive_generation_load(this->p);
    const uint64_t revision = prior ?
	bobol_identity_successor_or_terminate(prior->revision) : 1;
    const char *failureReason = NULL;
    const std::shared_ptr<BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_from_data(data, hierarchy, residentCut,
	    shadedCullBackfaces, revision, this->p->lineage,
	    &failureReason);
    if (!generation) {
	if (getenv("BOBOL_LOD_TRACE_REALIZATION"))
	    bu_log("BObol PoP realization rejected: %s\n",
		failureReason ? failureReason : "unknown");
	return FALSE;
    }
    progressive_generation_store(this->p, generation);
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::updateChunksFromCache(
    struct BObolMeshLod *lod,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    const std::vector<uint32_t> &chunkIds, int residentCut,
    SbBool shadedCullBackfaces)
{
    const auto fail = [](const char *reason) -> SbBool {
	if (getenv("BOBOL_CHUNK_DEBUG"))
	    bu_log("BObol chunk realization failed: %s\n", reason);
	return FALSE;
    };
    if (!this->p || !lod || !hierarchy.chunks || !hierarchy.chunk_count ||
	chunkIds.empty() || residentCut < hierarchy.min_cut ||
	residentCut > hierarchy.max_cut)
	return fail("arguments");
    for (size_t i = 0; i < chunkIds.size(); ++i)
	if (chunkIds[i] >= hierarchy.chunk_count ||
	    (i && chunkIds[i] <= chunkIds[i - 1]))
	    return fail("chunk set");

    std::lock_guard<std::mutex> updateLock(this->p->updateMutex);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> prior =
	progressive_generation_load(this->p);
    std::vector<std::shared_ptr<const
	BObolLodProgressiveMeshGeneration::ResidentChunk>> resident(
	hierarchy.chunk_count);
    if (prior && prior->residentChunks.size() == resident.size())
	resident = prior->residentChunks;

    const bool verboseTiming = getenv("BOBOL_DRAW_TIMING_VERBOSE") != NULL;
    int64_t readMicroseconds = 0;
    int64_t packMicroseconds = 0;

    std::vector<uint32_t> toRead;
    toRead.reserve(chunkIds.size());
    for (uint32_t chunkId : chunkIds) {
	const BObolMeshLodChunkInfo &info = hierarchy.chunks[chunkId];
	if (residentCut < info.min_cut)
	    continue;
	if (!resident[chunkId] ||
	    info.cuts[residentCut].point_count >
		resident[chunkId]->positions.size() ||
	    static_cast<uint64_t>(info.cuts[residentCut].face_count) * 3u >
		resident[chunkId]->indices.size())
	    toRead.push_back(chunkId);
    }
    if (!toRead.empty()) {
	BObolChunkReadContext context;
	context.hierarchy = &hierarchy;
	context.chunks = &resident;
	const int64_t readStarted = verboseTiming ? bu_gettime() : 0;
	const size_t workerCount = lod_spatial_worker_count(toRead.size());
	SbBool readSucceeded = FALSE;
	if (workerCount == 1) {
	    readSucceeded = bobol_mesh_lod_read_chunk_prefixes(
		lod, toRead.data(), toRead.size(), residentCut,
		bobol_progressive_chunk_read, &context) ? TRUE : FALSE;
	} else {
	    std::vector<BObolMeshLod *> readers(workerCount, NULL);
	    readers[0] = lod;
	    bool peersReady = true;
	    for (size_t worker = 1; worker < workerCount; ++worker) {
		readers[worker] = bobol_mesh_lod_clone_reader(lod);
		if (!readers[worker]) {
		    peersReady = false;
		    break;
		}
	    }
	    bool parallelSetupFailed = !peersReady;

	    if (peersReady) {
		std::atomic<bool> parallelSucceeded(true);
		std::vector<std::thread> workers;
		workers.reserve(workerCount - 1);
		const auto readRange = [&](size_t worker) {
		    const size_t first = toRead.size() * worker / workerCount;
		    const size_t end = toRead.size() * (worker + 1) / workerCount;
		    if (first == end)
			return;
		    if (!bobol_mesh_lod_read_chunk_prefixes(
			    readers[worker], toRead.data() + first, end - first,
			    residentCut, bobol_progressive_chunk_read, &context))
			parallelSucceeded.store(false, std::memory_order_relaxed);
		};
		bool threadsStarted = true;
		try {
		    for (size_t worker = 1; worker < workerCount; ++worker)
			workers.emplace_back(readRange, worker);
		} catch (const std::system_error &) {
		    threadsStarted = false;
		} catch (const std::bad_alloc &) {
		    threadsStarted = false;
		}
		parallelSetupFailed = !threadsStarted;
		if (threadsStarted)
		    readRange(0);
		for (std::thread &worker : workers)
		    worker.join();
		if (threadsStarted)
		    readSucceeded = parallelSucceeded.load(
			std::memory_order_relaxed) ? TRUE : FALSE;
	    }
	    for (size_t worker = 1; worker < readers.size(); ++worker)
		if (readers[worker])
		    bobol_mesh_lod_destroy(readers[worker]);

	    /* Failure to allocate a bounded peer/thread pool is an optimization
	     * failure, not a mesh failure.  Reuse the ordinary serial reader after
	     * every started peer has joined. */
	    if (parallelSetupFailed)
		readSucceeded = bobol_mesh_lod_read_chunk_prefixes(
		    lod, toRead.data(), toRead.size(), residentCut,
		    bobol_progressive_chunk_read, &context) ? TRUE : FALSE;
	}
	if (!readSucceeded)
	    return fail("cache read");
	if (verboseTiming)
	    readMicroseconds = std::max<int64_t>(0, bu_gettime() - readStarted);
    }
    for (uint32_t chunkId : chunkIds) {
	const BObolMeshLodChunkInfo &info = hierarchy.chunks[chunkId];
	if (residentCut >= info.min_cut &&
	    (!resident[chunkId] ||
	     info.cuts[residentCut].point_count >
		 resident[chunkId]->positions.size() ||
	     static_cast<uint64_t>(info.cuts[residentCut].face_count) * 3u >
		 resident[chunkId]->indices.size()))
	    return fail("incomplete cache read");
    }
    const uint64_t revision = prior ?
	bobol_identity_successor_or_terminate(prior->revision) : 1;

    const int64_t packStarted = verboseTiming ? bu_gettime() : 0;
    const std::shared_ptr<BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_from_chunks(hierarchy, resident,
	    shadedCullBackfaces, revision, this->p->lineage);
    if (verboseTiming)
	packMicroseconds = std::max<int64_t>(0, bu_gettime() - packStarted);
    if (!generation)
	return fail("renderer generation");
    progressive_generation_store(this->p, generation);
    if (verboseTiming)
	bu_log("[obol-timing] spatial prefix cut=%d selected=%zu read=%zu "
	       "read_time=%8.1f ms pack=%8.1f ms\n", residentCut,
	       chunkIds.size(), toRead.size(), readMicroseconds / 1000.0,
	       packMicroseconds / 1000.0);
    return TRUE;
}

struct BObolProgressiveSuffixBuilder {
    Obol::TriMesh *mesh = NULL;
    const struct BObolMeshLodHierarchyInfo *hierarchy = NULL;
    int expectedCut = -1;
};

static int
bobol_progressive_suffix_append(int cut, const point_t *points,
	size_t pointCount, const uint32_t *faces, size_t faceCount,
	const vect_t *normals, size_t normalCount, void *callbackData)
{
    BObolProgressiveSuffixBuilder *builder =
	static_cast<BObolProgressiveSuffixBuilder *>(callbackData);
    if (!builder || !builder->mesh || !builder->hierarchy ||
        cut != builder->expectedCut || cut <= 0 ||
        cut >= BOBOL_MESH_LOD_CUT_COUNT_MAX || normals || normalCount)
	return 0;
    Obol::TriMesh &mesh = *builder->mesh;
    const BObolMeshLodHierarchyInfo &hierarchy = *builder->hierarchy;
    const size_t priorPointCount = hierarchy.cuts[cut - 1].point_count;
    const size_t priorFaceCount = hierarchy.cuts[cut - 1].face_count;
    if (hierarchy.cuts[cut].point_count < priorPointCount ||
        hierarchy.cuts[cut].face_count < priorFaceCount ||
        pointCount != hierarchy.cuts[cut].point_count - priorPointCount ||
        faceCount != hierarchy.cuts[cut].face_count - priorFaceCount ||
	mesh.positions.size() != priorPointCount ||
	mesh.indices.size() != priorFaceCount * 3 ||
	(pointCount && !points) || (faceCount && !faces))
	return 0;

    for (size_t i = 0; i < pointCount; ++i) {
	if (!bobol_lod_point_in_domain(points[i],
		hierarchy.quantization_min,
		hierarchy.quantization_max))
	    return 0;
    }
    const size_t finalPointCount = priorPointCount + pointCount;
    for (size_t i = 0; i < faceCount * 3; ++i) {
	if (static_cast<size_t>(faces[i]) >= finalPointCount)
	    return 0;
    }
    for (size_t i = 0; i < pointCount; ++i) {
	mesh.positions.push_back(SbVec3f(
	    static_cast<float>(points[i][X]),
	    static_cast<float>(points[i][Y]),
	    static_cast<float>(points[i][Z])));
    }
    for (size_t i = 0; i < faceCount * 3; ++i) {
	mesh.indices.push_back(static_cast<uint32_t>(faces[i]));
    }
    builder->expectedCut++;
    return 1;
}

SbBool
BObolLodProgressiveMesh::extendFromCache(
    struct BObolMeshLod *lod,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int residentCut,
    SbBool shadedCullBackfaces)
{
    if (!this->p || !lod || hierarchy.has_normals ||
	residentCut < hierarchy.min_cut ||
	residentCut > hierarchy.max_cut ||
	residentCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return FALSE;

    std::lock_guard<std::mutex> updateLock(this->p->updateMutex);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> prior =
	progressive_generation_load(this->p);
    const Obol::TriMesh *priorMesh =
	prior && prior->shadedGeometry && prior->shadedGeometry->shaded ?
	    &*prior->shadedGeometry->shaded : NULL;
    if (!prior || !priorMesh || !priorMesh->normals.empty() ||
	prior->residentCut < hierarchy.min_cut ||
	residentCut <= prior->residentCut ||
	prior->minimumCut != hierarchy.min_cut ||
	prior->maximumCut != hierarchy.max_cut ||
	prior->shadedCullBackfaces !=
	    (shadedCullBackfaces ? TRUE : FALSE))
	return FALSE;
    for (uint32_t cut = 0; cut < hierarchy.cut_count; ++cut) {
	if (prior->pointCount[cut] != hierarchy.cuts[cut].point_count ||
	    prior->faceCount[cut] != hierarchy.cuts[cut].face_count)
	    return FALSE;
    }
    if (hierarchy.cuts[residentCut].point_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	hierarchy.cuts[residentCut].face_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) / 3 ||
	priorMesh->positions.size() !=
	    hierarchy.cuts[prior->residentCut].point_count ||
	priorMesh->indices.size() !=
	    hierarchy.cuts[prior->residentCut].face_count * 3)
	return FALSE;

    std::shared_ptr<BObolLodProgressiveMeshGeneration> generation(
	new BObolLodProgressiveMeshGeneration);
    generation->pointCount = prior->pointCount;
    generation->faceCount = prior->faceCount;
    generation->cuts = prior->cuts;
    generation->bounds = prior->bounds;
    generation->quantizationMinimum = prior->quantizationMinimum;
    generation->quantizationMaximum = prior->quantizationMaximum;
    generation->aggregateProxyCorners = prior->aggregateProxyCorners;
    generation->minimumCut = prior->minimumCut;
    generation->maximumCut = prior->maximumCut;
    generation->residentCut = residentCut;
    generation->revision =
	bobol_identity_successor_or_terminate(prior->revision);
    generation->shadedCullBackfaces = prior->shadedCullBackfaces;

    std::shared_ptr<Obol::PartGeometryBuilder> geometry(
	new Obol::PartGeometryBuilder);
    Obol::TriMesh mesh;
    mesh.bounds = priorMesh->bounds;
    mesh.progressiveCuts = priorMesh->progressiveCuts;
    mesh.progressiveMinimumCut = priorMesh->progressiveMinimumCut;
    mesh.progressiveResidentCut = priorMesh->progressiveResidentCut;
    mesh.progressiveQuantizationMinimum =
	priorMesh->progressiveQuantizationMinimum;
    mesh.progressiveQuantizationMaximum =
	priorMesh->progressiveQuantizationMaximum;
    mesh.progressiveLineage = priorMesh->progressiveLineage;
    mesh.progressiveClusters = priorMesh->progressiveClusters;
    mesh.progressiveClusterGridResolution =
	priorMesh->progressiveClusterGridResolution;
    mesh.positions.reserve(hierarchy.cuts[residentCut].point_count);
    mesh.positions.insert(mesh.positions.end(),
	priorMesh->positions.begin(), priorMesh->positions.end());
    mesh.indices.reserve(hierarchy.cuts[residentCut].face_count * 3);
    mesh.indices.insert(mesh.indices.end(),
	priorMesh->indices.begin(), priorMesh->indices.end());

    BObolProgressiveSuffixBuilder builder;
    builder.mesh = &mesh;
    builder.hierarchy = &hierarchy;
    builder.expectedCut = prior->residentCut + 1;
    if (!bobol_mesh_lod_read_resident_suffix(lod, prior->residentCut,
	    residentCut, bobol_progressive_suffix_append, &builder) ||
	builder.expectedCut != residentCut + 1 ||
	mesh.positions.size() != hierarchy.cuts[residentCut].point_count ||
	mesh.indices.size() != hierarchy.cuts[residentCut].face_count * 3)
	return FALSE;

    size_t scannedIndices =
	std::min<size_t>(priorMesh->progressiveCuts[
	    std::max(0, prior->residentCut)].indexCount,
	    mesh.indices.size());
    uint32_t maximumIndex = priorMesh->progressiveCuts[
	std::max(0, prior->residentCut)].positionCount ?
	priorMesh->progressiveCuts[
	    std::max(0, prior->residentCut)].positionCount - 1 : 0;
    for (int cut = prior->residentCut + 1;
	 cut < static_cast<int>(hierarchy.cut_count); ++cut) {
	const size_t cutIndexCount =
	    hierarchy.cuts[cut].face_count >
		    std::numeric_limits<size_t>::max() / 3 ?
		mesh.indices.size() :
		std::min(hierarchy.cuts[cut].face_count * 3,
		    mesh.indices.size());
	for (; scannedIndices < cutIndexCount; ++scannedIndices)
	    maximumIndex = std::max(maximumIndex,
		mesh.indices[scannedIndices]);
	mesh.progressiveCuts[cut].indexCount = static_cast<uint32_t>(
	    std::min(cutIndexCount,
		static_cast<size_t>(std::numeric_limits<uint32_t>::max())));
	mesh.progressiveCuts[cut].positionCount = cutIndexCount ?
	    maximumIndex + 1 : 0;
    }
    int drawableCut = residentCut;
    for (int cut = residentCut + 1;
	 cut <= hierarchy.max_cut &&
	 cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut) {
	if (hierarchy.cuts[cut].face_count >
		std::numeric_limits<size_t>::max() / 3 ||
	    hierarchy.cuts[cut].face_count * 3 > mesh.indices.size() ||
	    hierarchy.cuts[cut].point_count > mesh.positions.size())
	    break;
	drawableCut = cut;
    }
    mesh.progressiveResidentCut = static_cast<uint8_t>(drawableCut);
    geometry->shaded = std::move(mesh);
    geometry->shadedCullBackfaces =
	shadedCullBackfaces ? true : false;
    geometry->subpixelProxyEligible = true;
    geometry->aggregateProxyCorners = generation->aggregateProxyCorners;
    const std::shared_ptr<const Obol::PartGeometry> completed =
	bobol_cad_build_geometry_with_optional_proxy(
	    std::move(*geometry), "extended PoP generation");
    if (!completed)
	return FALSE;
    generation->shadedGeometry = completed;
    progressive_generation_store(this->p, generation);
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::commitTrim(
    const BObolLodProgressiveMeshTrimPtr &trim)
{
    if (!this->p || !trim || trim->owner != this->p ||
	!trim->source || !trim->generation)
	return FALSE;
    std::lock_guard<std::mutex> updateLock(this->p->updateMutex);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> current =
	progressive_generation_load(this->p);
    if (current != trim->source)
	return FALSE;
    progressive_generation_store(this->p, trim->generation);
    return TRUE;
}

BObolLodProgressiveMeshTrimPtr
BObolLodProgressiveMesh::prepareTrim(int residentCut) const
{
    if (!this->p)
	return BObolLodProgressiveMeshTrimPtr();
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> current =
	progressive_generation_load(this->p);
    const int cut = progressive_cut_clamp(current.get(), residentCut);
    if (cut < 0)
	return BObolLodProgressiveMeshTrimPtr();
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	cut == current->residentCut ? current :
	progressive_generation_prefix(*current, cut);
    if (!generation)
	return BObolLodProgressiveMeshTrimPtr();
    std::shared_ptr<BObolLodProgressiveMeshTrim> trim(
	new BObolLodProgressiveMeshTrim);
    trim->owner = this->p;
    trim->source = current;
    trim->generation = generation;
    return trim;
}

BObolLodProgressiveMeshTrimPtr
BObolLodProgressiveMesh::prepareTrim(
    int residentCut, const std::vector<uint32_t> &chunkIds) const
{
    if (!this->p || chunkIds.empty())
	return BObolLodProgressiveMeshTrimPtr();

    std::vector<BObolLodChunkCut> chunkCuts;
    chunkCuts.reserve(chunkIds.size());
    for (size_t i = 0; i < chunkIds.size(); ++i) {
	if (i && chunkIds[i] <= chunkIds[i - 1])
	    return BObolLodProgressiveMeshTrimPtr();
	chunkCuts.push_back({chunkIds[i], residentCut});
    }
    return this->prepareTrim(chunkCuts);
}

BObolLodProgressiveMeshTrimPtr
BObolLodProgressiveMesh::prepareTrim(
    const std::vector<BObolLodChunkCut> &chunkCuts) const
{
    if (!this->p || chunkCuts.empty())
	return BObolLodProgressiveMeshTrimPtr();
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> current =
	progressive_generation_load(this->p);
    int requestedCut = -1;
    for (size_t i = 0; i < chunkCuts.size(); ++i) {
	if (chunkCuts[i].cut < 0 ||
	    (i && chunkCuts[i].chunkId <= chunkCuts[i - 1].chunkId))
	    return BObolLodProgressiveMeshTrimPtr();
	requestedCut = std::max(requestedCut, chunkCuts[i].cut);
    }
    const int cut = progressive_cut_clamp(current.get(), requestedCut);
    if (cut < 0)
	return BObolLodProgressiveMeshTrimPtr();
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_prefix(*current, cut, &chunkCuts);
    if (!generation)
	return BObolLodProgressiveMeshTrimPtr();
    std::shared_ptr<BObolLodProgressiveMeshTrim> trim(
	new BObolLodProgressiveMeshTrim);
    trim->owner = this->p;
    trim->source = current;
    trim->generation = generation;
    return trim;
}

SbBool
BObolLodProgressiveMesh::trim(int residentCut)
{
    const BObolLodProgressiveMeshTrimPtr trim =
	this->prepareTrim(residentCut);
    return trim && this->commitTrim(trim) ? TRUE : FALSE;
}

SbBool
BObolLodProgressiveMesh::copyCut(
    BObolLodMeshPayload &payload, int requestedCut) const
{
    payload.clear();
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const int cut =
	progressive_cut_clamp(generation.get(), requestedCut);
    if (cut < 0) {
	if (getenv("BOBOL_CHUNK_DEBUG"))
	    bu_log("BObol copyCut has no cut requested=%d generation=%d/%d/%d\n",
		requestedCut, generation ? generation->minimumCut : -99,
		generation ? generation->residentCut : -99,
		generation ? generation->maximumCut : -99);
	return FALSE;
    }

    /* Direct shape clients still receive one conventional payload.  The CAD
     * scene path publishes spatial pages directly and never pays this merge;
     * materialize it only for the explicitly requested compatibility-free
     * direct result. */
    if (!generation->chunks.empty()) {
	if (!canDrawCut(cut))
	    return FALSE;
	for (size_t chunkId = 0; chunkId < generation->chunks.size();
		++chunkId) {
	    const auto &metadata = generation->chunks[chunkId];
	    const size_t pointCount = metadata.cuts[cut].point_count;
	    const size_t indexCount =
		static_cast<size_t>(metadata.cuts[cut].face_count) * 3u;
	    if (!indexCount)
		continue;
	    const auto &resident = generation->residentChunks[chunkId];
	    if (!resident || pointCount > resident->positions.size() ||
		indexCount > resident->indices.size() ||
		payload.points.size() >
		    static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
		pointCount > static_cast<size_t>(
		    std::numeric_limits<int32_t>::max()) - payload.points.size()) {
		payload.clear();
		return FALSE;
	    }
	    const size_t vertexBase = payload.points.size();
	    payload.points.insert(payload.points.end(),
		resident->positions.begin(),
		resident->positions.begin() + pointCount);
	    for (size_t offset = 0; offset < indexCount; ++offset) {
		const uint32_t sourceIndex = resident->indices[offset];
		if (sourceIndex >= pointCount) {
		    payload.clear();
		    return FALSE;
		}
		payload.coordIndex.push_back(static_cast<int32_t>(
		    vertexBase + sourceIndex));
		if (!resident->cornerNormals.empty()) {
		    if (offset >= resident->cornerNormals.size()) {
			payload.clear();
			return FALSE;
		    }
		    payload.normals.push_back(resident->cornerNormals[offset]);
		}
	    }
	}
	return payload.isValid();
    }

    const Obol::TriMesh *mesh =
	generation->shadedGeometry && generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!mesh) {
	if (getenv("BOBOL_CHUNK_DEBUG"))
	    bu_log("BObol copyCut has no shaded mesh cut=%d\n", cut);
	return FALSE;
    }
    if (!generation->chunks.empty() &&
	mesh->hasAdaptiveProgressiveClusters()) {
	payload.points = mesh->positions;
	for (const Obol::ProgressiveTriangleCluster &cluster :
		mesh->progressiveClusters) {
	    for (const Obol::ProgressiveTriangleClusterRange &range :
		    cluster.ranges) {
		if (range.activationCut > cut)
		    continue;
		const uint64_t end = static_cast<uint64_t>(range.firstIndex) +
		    range.indexCount;
		if (range.firstIndex % 3 || range.indexCount % 3 ||
		    end > mesh->indices.size()) {
		    if (getenv("BOBOL_CHUNK_DEBUG"))
			bu_log("BObol copyCut invalid range cut=%d first=%u "
			    "count=%u mesh=%zu\n", cut, range.firstIndex,
			    range.indexCount, mesh->indices.size());
		    payload.clear();
		    return FALSE;
		}
		for (uint64_t indexOffset = range.firstIndex;
			indexOffset < end; ++indexOffset) {
		    const uint32_t index = mesh->indices[indexOffset];
		    if (index >= payload.points.size() || index > INT32_MAX) {
			if (getenv("BOBOL_CHUNK_DEBUG"))
			    bu_log("BObol copyCut invalid index cut=%d index=%u "
				"points=%zu\n", cut, index,
				payload.points.size());
			payload.clear();
			return FALSE;
		    }
		    payload.coordIndex.push_back(static_cast<int32_t>(index));
		    if (!mesh->normals.empty()) {
			if (index >= mesh->normals.size()) {
			    if (getenv("BOBOL_CHUNK_DEBUG"))
				bu_log("BObol copyCut invalid normal cut=%d "
				    "index=%u normals=%zu\n", cut, index,
				    mesh->normals.size());
			    payload.clear();
			    return FALSE;
			}
			payload.normals.push_back(mesh->normals[index]);
		    }
		}
	    }
	}
	const SbBool valid = payload.isValid();
	if (!valid && getenv("BOBOL_CHUNK_DEBUG"))
	    bu_log("BObol copyCut adaptive failed cut=%d chunks=%zu "
		"clusters=%zu positions=%zu indices=%zu output=%zu normals=%zu\n",
		cut, generation->chunks.size(),
		mesh->progressiveClusters.size(), mesh->positions.size(),
		mesh->indices.size(), payload.coordIndex.size(),
		payload.normals.size());
	return valid;
    }
    const size_t points =
	mesh->positionCountAtCut(static_cast<uint8_t>(cut));
    const size_t indices =
	mesh->indexCountAtCut(static_cast<uint8_t>(cut));
    if (!points || points > mesh->positions.size() ||
	indices < 3 || indices > mesh->indices.size()) {
	if (getenv("BOBOL_CHUNK_DEBUG"))
	    bu_log("BObol copyCut prefix failed cut=%d points=%zu/%zu "
		"indices=%zu/%zu resident=%u chunks=%zu\n", cut,
		points, mesh->positions.size(), indices, mesh->indices.size(),
		static_cast<unsigned>(mesh->progressiveResidentCut),
		generation->chunks.size());
	return FALSE;
    }
    payload.points.assign(mesh->positions.begin(),
	mesh->positions.begin() + points);
    payload.coordIndex.reserve(indices);
    if (!mesh->normals.empty())
	payload.normals.reserve(indices);
    for (size_t i = 0; i < indices; ++i) {
	const uint32_t index = mesh->indices[i];
	if (index >= points ||
	    index > static_cast<uint32_t>(
		std::numeric_limits<int32_t>::max())) {
	    if (getenv("BOBOL_CHUNK_DEBUG"))
		bu_log("BObol copyCut prefix invalid index cut=%d "
		    "index=%u points=%zu\n", cut, index, points);
	    payload.clear();
	    return FALSE;
	}
	payload.coordIndex.push_back(static_cast<int32_t>(index));
	if (!mesh->normals.empty()) {
	    if (index >= mesh->normals.size()) {
		if (getenv("BOBOL_CHUNK_DEBUG"))
		    bu_log("BObol copyCut prefix invalid normal cut=%d "
			"index=%u normals=%zu\n", cut, index,
			mesh->normals.size());
		payload.clear();
		return FALSE;
	    }
	    payload.normals.push_back(mesh->normals[index]);
	}
    }
    const SbBool valid = payload.isValid();
    if (!valid && getenv("BOBOL_CHUNK_DEBUG"))
	bu_log("BObol copyCut prefix output invalid cut=%d points=%zu "
	    "indices=%zu normals=%zu\n", cut, payload.points.size(),
	    payload.coordIndex.size(), payload.normals.size());
    return valid;
}

SbBool
BObolLodProgressiveMesh::isValid(void) const
{
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const Obol::TriMesh *mesh =
	generation && generation->shadedGeometry &&
	generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!generation || generation->residentCut < generation->minimumCut)
	return FALSE;
    if (!generation->chunks.empty()) {
	for (const auto &resident : generation->residentChunks)
	    if (resident && !resident->positions.empty() &&
		resident->indices.size() >= 3 &&
		resident->indices.size() % 3 == 0)
		return TRUE;
	return FALSE;
    }
    return mesh && !mesh->positions.empty() && mesh->indices.size() >= 3 &&
	mesh->indices.size() % 3 == 0 &&
	(mesh->normals.empty() ||
	 mesh->normals.size() == mesh->positions.size()) ? TRUE : FALSE;
}

SbBool
BObolLodProgressiveMesh::canDrawCut(int requestedCut) const
{
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation ||
	requestedCut < generation->minimumCut ||
	requestedCut > generation->maximumCut ||
	requestedCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return FALSE;
    const Obol::TriMesh *mesh =
	generation->shadedGeometry && generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!generation->chunks.empty()) {
	/* This whole-leaf query is still used by conservative clients.  Keep it
	 * allocation-free: it is reachable from per-frame admission and a large
	 * scene may contain tens of thousands of private pages. */
	for (size_t chunkId = 0; chunkId < generation->chunks.size(); ++chunkId) {
	    const auto &metadata = generation->chunks[chunkId];
	    if (!metadata.cuts[requestedCut].face_count)
		continue;
	    if (chunkId >= generation->residentChunks.size() ||
		!generation->residentChunks[chunkId])
		return FALSE;
	    const auto &resident = generation->residentChunks[chunkId];
	    if (metadata.cuts[requestedCut].point_count >
		    resident->positions.size() ||
		static_cast<uint64_t>(metadata.cuts[requestedCut].face_count) *
		    3u > resident->indices.size())
		return FALSE;
	}
	return TRUE;
    }
    if (!mesh)
	return FALSE;
    return mesh->isProgressive() &&
	requestedCut <= static_cast<int>(
	    mesh->progressiveResidentCut) ? TRUE : FALSE;
}

SbBool
BObolLodProgressiveMesh::countsForChunksAtCut(
    const std::vector<uint32_t> &chunkIds, int requestedCut,
    SbBool hasNormals, BObolLodCounts *counts) const
{
    if (!counts)
	return FALSE;
    if (!canDrawChunksAtCut(chunkIds, requestedCut)) {
	counts->clear();
	return FALSE;
    }
    return hierarchyCountsForChunksAtCut(
	chunkIds, requestedCut, hasNormals, counts);
}

SbBool
BObolLodProgressiveMesh::drawableCountsAtCuts(
    const std::vector<uint32_t> &chunkIds, SbBool hasNormals,
    BObolLodCounts *counts, size_t count, int *minimumCut,
    int *maximumDrawableCut) const
{
    if (minimumCut)
	*minimumCut = -1;
    if (maximumDrawableCut)
	*maximumDrawableCut = -1;
    if (!counts || count < BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return FALSE;
    for (size_t cut = 0; cut < count; ++cut)
	counts[cut].clear();
    if (!this->p)
	return FALSE;

    /* One atomic shared-generation retain is the point of this API.  Do all
     * population and residency decisions against that exact immutable
     * snapshot so a concurrent suffix publication cannot mix generations. */
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const Obol::TriMesh *mesh =
	generation && generation->shadedGeometry &&
	generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!generation || generation->residentCut < generation->minimumCut ||
	generation->minimumCut < 0 ||
	generation->maximumCut < generation->minimumCut ||
	generation->maximumCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return FALSE;

    if (minimumCut)
	*minimumCut = generation->minimumCut;

    if (generation->chunks.empty()) {
	/* A nonempty selection names private spatial pages.  A whole-prefix
	 * preview can draw the same logical leaf, but it cannot prove residency
	 * for page identifiers belonging to a later durable spatial generation. */
	if (!chunkIds.empty())
	    return FALSE;
	if (!mesh || mesh->positions.empty() || mesh->indices.size() < 3 ||
	    mesh->indices.size() % 3 != 0 ||
	    (!mesh->normals.empty() &&
	     mesh->normals.size() != mesh->positions.size()))
	    return FALSE;
	for (int cut = generation->minimumCut;
		cut <= generation->maximumCut; ++cut) {
	    BObolLodCounts &population = counts[static_cast<size_t>(cut)];
	    population.faceCount = generation->faceCount[cut];
	    population.pointCount = generation->pointCount[cut];
	    population.originalPointCount = population.pointCount;
	    population.normalCount = hasNormals ?
		(population.faceCount > UINT64_MAX / 3 ? UINT64_MAX :
		 population.faceCount * 3) : 0;
	}
	const int drawable = mesh->isProgressive() ?
	    std::min(generation->maximumCut,
		static_cast<int>(mesh->progressiveResidentCut)) : -1;
	if (maximumDrawableCut && drawable >= generation->minimumCut)
	    *maximumDrawableCut = drawable;
	return drawable >= generation->minimumCut ? TRUE : FALSE;
    }

    /* Empty selects the complete logical leaf, matching canDrawCut and the
     * ordinary whole-prefix accounting path.  A nonempty selection must be
     * canonical so results remain deterministic and no page is double
     * charged. */
    if (!chunkIds.empty()) {
	for (size_t i = 0; i < chunkIds.size(); ++i) {
	    if (chunkIds[i] >= generation->chunks.size() ||
		(i && chunkIds[i] <= chunkIds[i - 1]))
		return FALSE;
	}
    }

    int drawable = generation->minimumCut - 1;
    for (int cut = generation->minimumCut;
	    cut <= generation->maximumCut; ++cut) {
	BObolLodCounts &population = counts[static_cast<size_t>(cut)];
	if (chunkIds.empty()) {
	    /* Preserve the ordinary whole-prefix contract exactly.  A producer
	     * may account shared hierarchy vertices differently from a sum of its
	     * private page metadata. */
	    population.faceCount = generation->faceCount[cut];
	    population.pointCount = generation->pointCount[cut];
	}
	bool resident = true;
	const size_t selectedCount = chunkIds.empty() ?
	    generation->chunks.size() : chunkIds.size();
	for (size_t selected = 0; selected < selectedCount; ++selected) {
	    const uint32_t chunkId = chunkIds.empty() ?
		static_cast<uint32_t>(selected) : chunkIds[selected];
	    const BObolMeshLodChunkCutInfo &chunkCut =
		generation->chunks[chunkId].cuts[cut];
	    if (!chunkIds.empty()) {
		population.faceCount = chunkCut.face_count >
			UINT64_MAX - population.faceCount ? UINT64_MAX :
		    population.faceCount + chunkCut.face_count;
		population.pointCount = chunkCut.point_count >
			UINT64_MAX - population.pointCount ? UINT64_MAX :
		    population.pointCount + chunkCut.point_count;
	    }
	    if (!chunkCut.face_count)
		continue;
	    if (chunkId >= generation->residentChunks.size() ||
		!generation->residentChunks[chunkId] ||
		chunkCut.point_count >
		    generation->residentChunks[chunkId]->positions.size() ||
		static_cast<uint64_t>(chunkCut.face_count) * 3u >
		    generation->residentChunks[chunkId]->indices.size())
		resident = false;
	}
	population.originalPointCount = population.pointCount;
	population.normalCount = hasNormals ?
	    (population.faceCount > UINT64_MAX / 3 ? UINT64_MAX :
	     population.faceCount * 3) : 0;
	/* Cumulative page prefixes cannot become drawable again above a missing
	 * resident prefix.  Stop advancing the frontier but retain hierarchy
	 * counts for diagnostics and future capacity planning. */
	if (resident && drawable == cut - 1)
	    drawable = cut;
    }
    if (maximumDrawableCut)
	*maximumDrawableCut = drawable;
    return drawable >= generation->minimumCut ? TRUE : FALSE;
}

SbBool
BObolLodProgressiveMesh::hierarchyCountsForChunksAtCut(
    const std::vector<uint32_t> &chunkIds, int requestedCut,
    SbBool hasNormals, BObolLodCounts *counts) const
{
    if (!counts)
	return FALSE;
    counts->clear();
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || generation->chunks.empty() || chunkIds.empty() ||
	requestedCut < generation->minimumCut ||
	requestedCut > generation->maximumCut)
	return FALSE;
    for (size_t i = 0; i < chunkIds.size(); ++i) {
	const uint32_t chunkId = chunkIds[i];
	if (chunkId >= generation->chunks.size() ||
	    (i && chunkId <= chunkIds[i - 1])) {
	    counts->clear();
	    return FALSE;
	}
	const BObolMeshLodChunkCutInfo &cut =
	    generation->chunks[chunkId].cuts[requestedCut];
	counts->faceCount = cut.face_count >
		UINT64_MAX - counts->faceCount ? UINT64_MAX :
	    counts->faceCount + cut.face_count;
	counts->pointCount = cut.point_count >
		UINT64_MAX - counts->pointCount ? UINT64_MAX :
	    counts->pointCount + cut.point_count;
    }
    counts->originalPointCount = counts->pointCount;
    counts->normalCount = hasNormals ?
	(counts->faceCount > UINT64_MAX / 3 ? UINT64_MAX :
	    counts->faceCount * 3) : 0;
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::canDrawChunksAtCut(
    const std::vector<uint32_t> &chunkIds, int requestedCut) const
{
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation)
	return FALSE;
    if (generation->chunks.empty())
	return FALSE;
    if (requestedCut < generation->minimumCut ||
	requestedCut > generation->maximumCut || chunkIds.empty())
	return FALSE;
    for (size_t i = 0; i < chunkIds.size(); ++i) {
	const uint32_t chunkId = chunkIds[i];
	if (chunkId >= generation->chunks.size() ||
	    (i && chunkId <= chunkIds[i - 1]))
	    return FALSE;
	const auto &metadata = generation->chunks[chunkId];
	if (!metadata.cuts[requestedCut].face_count)
	    continue;
	if (chunkId >= generation->residentChunks.size() ||
	    !generation->residentChunks[chunkId])
	    return FALSE;
	const auto &resident = generation->residentChunks[chunkId];
	if (metadata.cuts[requestedCut].point_count >
		resident->positions.size() ||
	    static_cast<uint64_t>(metadata.cuts[requestedCut].face_count) * 3u >
		resident->indices.size())
	    return FALSE;
    }
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::populatedChunkIdsAtCut(
    const std::vector<uint32_t> &chunkIds, int requestedCut,
    std::vector<uint32_t> &populatedChunkIds) const
{
    populatedChunkIds.clear();
    if (!this->p || chunkIds.empty())
	return FALSE;

    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || generation->chunks.empty() ||
	requestedCut < generation->minimumCut ||
	requestedCut > generation->maximumCut)
	return FALSE;

    populatedChunkIds.reserve(chunkIds.size());
    for (size_t i = 0; i < chunkIds.size(); ++i) {
	const uint32_t chunkId = chunkIds[i];
	if (chunkId >= generation->chunks.size() ||
	    (i && chunkId <= chunkIds[i - 1])) {
	    populatedChunkIds.clear();
	    return FALSE;
	}
	if (generation->chunks[chunkId].cuts[requestedCut].face_count)
	    populatedChunkIds.push_back(chunkId);
    }
    return TRUE;
}

int
BObolLodProgressiveMesh::residentCutForChunks(
    const std::vector<uint32_t> &chunkIds) const
{
    if (!this->p)
	return -1;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation)
	return -1;
    if (generation->chunks.empty())
	return chunkIds.empty() ? generation->residentCut : -1;
    if (chunkIds.empty())
	return -1;
    int frontier = generation->minimumCut - 1;
    for (int cut = generation->minimumCut;
	 cut <= generation->maximumCut; ++cut) {
	bool complete = true;
	for (size_t i = 0; i < chunkIds.size(); ++i) {
	    const uint32_t chunkId = chunkIds[i];
	    if (chunkId >= generation->chunks.size() ||
		(i && chunkId <= chunkIds[i - 1]))
		return -1;
	    if (!generation->chunks[chunkId].cuts[cut].face_count)
		continue;
	    if (chunkId >= generation->residentChunks.size() ||
		!generation->residentChunks[chunkId] ||
		generation->chunks[chunkId].cuts[cut].point_count >
		    generation->residentChunks[chunkId]->positions.size() ||
		static_cast<uint64_t>(
		    generation->chunks[chunkId].cuts[cut].face_count) * 3u >
		    generation->residentChunks[chunkId]->indices.size()) {
		complete = false;
		break;
	    }
	}
	if (!complete)
	    break;
	frontier = cut;
    }
    return frontier;
}

void
BObolLodProgressiveMesh::residentChunkIds(
    std::vector<uint32_t> &chunkIds) const
{
    chunkIds.clear();
    if (!this->p)
	return;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || generation->chunks.empty())
	return;
    chunkIds.reserve(generation->residentChunks.size());
    for (size_t chunk = 0; chunk < generation->residentChunks.size(); ++chunk)
	if (generation->residentChunks[chunk])
	    chunkIds.push_back(static_cast<uint32_t>(chunk));
}

void
BObolLodProgressiveMesh::residentChunkCuts(
    std::vector<BObolLodChunkCut> &chunkCuts) const
{
    chunkCuts.clear();
    if (!this->p)
	return;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || generation->chunks.empty())
	return;
    chunkCuts.reserve(generation->residentChunks.size());
    for (size_t chunk = 0; chunk < generation->residentChunks.size(); ++chunk)
	if (generation->residentChunks[chunk])
	    chunkCuts.push_back({static_cast<uint32_t>(chunk),
		generation->residentChunks[chunk]->residentCut});
}

int
BObolLodProgressiveMesh::minimumCut(void) const
{
    if (!this->p)
	return -1;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->minimumCut : -1;
}

int
BObolLodProgressiveMesh::maximumCut(void) const
{
    if (!this->p)
	return -1;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->maximumCut : -1;
}

int
BObolLodProgressiveMesh::residentCut(void) const
{
    if (!this->p)
	return -1;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->residentCut : -1;
}

uint64_t
BObolLodProgressiveMesh::revision(void) const
{
    if (!this->p)
	return 0;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->revision : 0;
}

BObolLodProgressiveMeshSnapshot
BObolLodProgressiveMesh::snapshot(void) const
{
    return BObolLodProgressiveMeshSnapshot(
	this->p ? progressive_generation_load(this->p) :
	    std::shared_ptr<const BObolLodProgressiveMeshGeneration>());
}

size_t
BObolLodProgressiveMesh::pointCount(int requestedCut) const
{
    if (!this->p)
	return 0;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const int cut =
	progressive_cut_clamp(generation.get(), requestedCut);
    return cut >= 0 ? generation->pointCount[cut] : 0;
}

size_t
BObolLodProgressiveMesh::faceCount(int requestedCut) const
{
    if (!this->p)
	return 0;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const int cut =
	progressive_cut_clamp(generation.get(), requestedCut);
    return cut >= 0 ? generation->faceCount[cut] : 0;
}

size_t
BObolLodProgressiveMesh::hierarchyPointCountAtCut(int requestedCut) const
{
    if (!this->p)
	return 0;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation)
	return 0;
    if (requestedCut < generation->minimumCut)
	requestedCut = generation->minimumCut;
    if (requestedCut > generation->maximumCut)
	requestedCut = generation->maximumCut;
    return requestedCut >= 0 ?
	generation->pointCount[requestedCut] : 0;
}

size_t
BObolLodProgressiveMesh::hierarchyFaceCountAtCut(int requestedCut) const
{
    if (!this->p)
	return 0;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation)
	return 0;
    if (requestedCut < generation->minimumCut)
	requestedCut = generation->minimumCut;
    if (requestedCut > generation->maximumCut)
	requestedCut = generation->maximumCut;
    return requestedCut >= 0 ?
	generation->faceCount[requestedCut] : 0;
}

SbBool
BObolLodProgressiveMesh::hasSpatialClusters(void) const
{
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const Obol::TriMesh *mesh = generation && generation->shadedGeometry &&
	generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    return generation && (!generation->chunks.empty() ||
	(mesh && mesh->hasProgressiveClusters())) ? TRUE : FALSE;
}

static void
bobol_lod_transform_bounds(const SbBox3f &local, const SbMatrix &transform,
    SbBox3f &world)
{
    world.makeEmpty();
    if (local.isEmpty())
	return;
    const SbVec3f minimum = local.getMin();
    const SbVec3f maximum = local.getMax();
    for (unsigned int corner = 0; corner < 8; ++corner) {
	const SbVec3f point(
	    corner & 1u ? maximum[0] : minimum[0],
	    corner & 2u ? maximum[1] : minimum[1],
	    corner & 4u ? maximum[2] : minimum[2]);
	SbVec3f transformed;
	transform.multVecMatrix(point, transformed);
	world.extendBy(transformed);
    }
}

static int
bobol_lod_bounds_frustum_relation(const SbBox3f &bounds,
    const SbMatrix &viewProjection)
{
    if (bounds.isEmpty())
	return 0;
    const SbVec3f minimum = bounds.getMin();
    const SbVec3f maximum = bounds.getMax();
    bool whollyInside = true;
    for (int column = 0; column < 3; ++column) {
	for (int signIndex = 0; signIndex < 2; ++signIndex) {
	    const float sign = signIndex == 0 ? 1.0f : -1.0f;
	    const float a = sign * viewProjection[0][column] +
		viewProjection[0][3];
	    const float b = sign * viewProjection[1][column] +
		viewProjection[1][3];
	    const float c = sign * viewProjection[2][column] +
		viewProjection[2][3];
	    const float d = sign * viewProjection[3][column] +
		viewProjection[3][3];
	    const float positiveX = a < 0.0f ? minimum[0] : maximum[0];
	    const float positiveY = b < 0.0f ? minimum[1] : maximum[1];
	    const float positiveZ = c < 0.0f ? minimum[2] : maximum[2];
	    if (a * positiveX + b * positiveY + c * positiveZ + d < 0.0f)
		return -1;
	    const float negativeX = a < 0.0f ? maximum[0] : minimum[0];
	    const float negativeY = b < 0.0f ? maximum[1] : minimum[1];
	    const float negativeZ = c < 0.0f ? maximum[2] : minimum[2];
	    if (a * negativeX + b * negativeY + c * negativeZ + d < 0.0f)
		whollyInside = false;
	}
    }
    return whollyInside ? 1 : 0;
}

SbBool
bobol_lod_visible_chunks(
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    const SbMatrix &localToRoot, const SbMatrix &viewProjection,
    std::vector<uint32_t> &chunkIds)
{
    chunkIds.clear();
    if (!hierarchy.chunks || !hierarchy.chunk_count)
	return FALSE;
    chunkIds.reserve(hierarchy.chunk_count);
    for (uint32_t chunk = 0; chunk < hierarchy.chunk_count; ++chunk) {
	const BObolMeshLodChunkInfo &info = hierarchy.chunks[chunk];
	if (info.chunk_id != chunk)
	    return FALSE;
	const SbBox3f local(
	    SbVec3f(static_cast<float>(info.bmin[X]),
		static_cast<float>(info.bmin[Y]),
		static_cast<float>(info.bmin[Z])),
	    SbVec3f(static_cast<float>(info.bmax[X]),
		static_cast<float>(info.bmax[Y]),
		static_cast<float>(info.bmax[Z])));
	SbBox3f world;
	bobol_lod_transform_bounds(local, localToRoot, world);
	if (bobol_lod_bounds_frustum_relation(world, viewProjection) >= 0)
	    chunkIds.push_back(chunk);
    }
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::visibleChunkIds(
    const SbMatrix &localToRoot, const SbMatrix &viewProjection,
    std::vector<uint32_t> &chunkIds) const
{
    chunkIds.clear();
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || generation->chunks.empty())
	return FALSE;
    chunkIds.reserve(generation->chunks.size());
    for (const auto &chunk : generation->chunks) {
	SbBox3f world;
	bobol_lod_transform_bounds(chunk.bounds, localToRoot, world);
	if (bobol_lod_bounds_frustum_relation(world, viewProjection) >= 0)
	    chunkIds.push_back(chunk.chunkId);
    }
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::visibleCountsAtCuts(
    const SbMatrix &localToRoot, const SbMatrix &viewProjection,
    SbBool hasNormals, BObolLodCounts *counts, size_t count) const
{
    if (!this->p || !counts || count < BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (generation && !generation->chunks.empty()) {
	/* A complete occurrence uses the ordinary whole-prefix table.  Building
	 * and retaining a second identical per-cut population for every fully
	 * visible chunked leaf made the view census scale with scene size rather
	 * than with the small clipped boundary.  Return a spatial population only
	 * when the occurrence actually straddles the frustum, matching the public
	 * contract and the monolithic-cluster path below. */
	SbBox3f worldBounds;
	bobol_lod_transform_bounds(
	    generation->bounds, localToRoot, worldBounds);
	if (bobol_lod_bounds_frustum_relation(
		worldBounds, viewProjection) != 0)
	    return FALSE;
	for (size_t cut = 0; cut < count; ++cut)
	    counts[cut].clear();
	for (const auto &chunk : generation->chunks) {
	    SbBox3f worldChunk;
	    bobol_lod_transform_bounds(
		chunk.bounds, localToRoot, worldChunk);
	    if (bobol_lod_bounds_frustum_relation(
		    worldChunk, viewProjection) < 0)
		continue;
	    for (size_t cut = 0;
		    cut < count && cut < generation->cuts.size(); ++cut) {
		const uint64_t faces = chunk.cuts[cut].face_count;
		const uint64_t points = chunk.cuts[cut].point_count;
		counts[cut].faceCount = counts[cut].faceCount >
		    UINT64_MAX - faces ? UINT64_MAX :
		    counts[cut].faceCount + faces;
		counts[cut].pointCount = counts[cut].pointCount >
		    UINT64_MAX - points ? UINT64_MAX :
		    counts[cut].pointCount + points;
		counts[cut].originalPointCount = counts[cut].pointCount;
		if (hasNormals) {
		    const uint64_t normals = faces > UINT64_MAX / 3 ?
			UINT64_MAX : faces * 3;
		    counts[cut].normalCount = counts[cut].normalCount >
			UINT64_MAX - normals ? UINT64_MAX :
			counts[cut].normalCount + normals;
		}
	    }
	}
	return TRUE;
    }
    const Obol::TriMesh *mesh = generation && generation->shadedGeometry &&
	generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!mesh || !mesh->hasProgressiveClusters())
	return FALSE;

    SbBox3f worldBounds;
    bobol_lod_transform_bounds(generation->bounds, localToRoot, worldBounds);
    if (bobol_lod_bounds_frustum_relation(
	    worldBounds, viewProjection) != 0)
	return FALSE;

    for (size_t cut = 0; cut < count; ++cut)
	counts[cut].clear();
    std::array<uint64_t, BOBOL_MESH_LOD_CUT_COUNT_MAX> indexDeltas = {};
    for (const Obol::ProgressiveTriangleCluster &cluster :
	 mesh->progressiveClusters) {
	if (cluster.ranges.empty())
	    continue;
	SbBox3f worldCluster;
	bobol_lod_transform_bounds(
	    cluster.bounds, localToRoot, worldCluster);
	if (bobol_lod_bounds_frustum_relation(
		worldCluster, viewProjection) < 0)
	    continue;
	for (const Obol::ProgressiveTriangleClusterRange &range :
	     cluster.ranges) {
	    if (range.activationCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
		continue;
	    const uint64_t prior = indexDeltas[range.activationCut];
	    indexDeltas[range.activationCut] =
		range.indexCount > UINT64_MAX - prior ? UINT64_MAX :
		prior + range.indexCount;
	}
    }
    uint64_t visibleIndices = 0;
    for (size_t cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut) {
	visibleIndices = indexDeltas[cut] > UINT64_MAX - visibleIndices ?
	    UINT64_MAX : visibleIndices + indexDeltas[cut];
	counts[cut].faceCount = visibleIndices / 3;
	/* Partial-cluster draws reference the shared global vertex array, but
	 * their vertex-processing work is one position per submitted index.
	 * Charging that expanded stream matches both the fixed and shader paths
	 * without an O(triangles) unique-index census in the camera hot path. */
	counts[cut].pointCount = visibleIndices;
	counts[cut].originalPointCount = visibleIndices;
	counts[cut].normalCount = hasNormals ? visibleIndices : 0;
    }
    return TRUE;
}

SbBool
BObolLodProgressiveMesh::cutInfo(
    int requestedCut, struct BObolMeshLodCutInfo *info) const
{
    if (!this->p || !info || requestedCut < 0)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation ||
	static_cast<size_t>(requestedCut) >= generation->cuts.size())
	return FALSE;
    *info = generation->cuts[static_cast<size_t>(requestedCut)];
    return TRUE;
}

int
BObolLodProgressiveMesh::cutForScreenError(
    double projectedPixelDiameter, double targetPixelError) const
{
    if (!this->p)
	return -1;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return progressive_cut_for_screen_error(generation.get(),
	projectedPixelDiameter, targetPixelError);
}

double
BObolLodProgressiveMesh::projectedErrorAtCut(
    int cut, double projectedPixelDiameter) const
{
    if (!this->p)
	return std::numeric_limits<double>::infinity();
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return progressive_projected_error_at_cut(generation.get(), cut,
	projectedPixelDiameter);
}

BObolLodCounts
bobol_lod_progressive_counts(
    const BObolLodProgressiveMeshPtr &progressiveMesh, int cut,
    SbBool hasNormals)
{
    BObolLodCounts counts;
    if (!progressiveMesh || !progressiveMesh->isValid())
	return counts;
    counts.faceCount = progressiveMesh->hierarchyFaceCountAtCut(cut);
    counts.pointCount = progressiveMesh->hierarchyPointCountAtCut(cut);
    counts.originalPointCount = counts.pointCount;
    counts.normalCount = hasNormals ?
	(counts.faceCount > UINT64_MAX / 3 ?
	    UINT64_MAX : counts.faceCount * 3) : 0;
    return counts;
}

size_t
BObolLodProgressiveMesh::estimateBytes(void) const
{
    if (!this->p)
	return 0;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (generation && !generation->chunks.empty()) {
	size_t bytes = 0;
	for (const auto &resident : generation->residentChunks) {
	    if (!resident)
		continue;
	    const size_t positions = resident->positions.capacity() *
		sizeof(SbVec3f);
	    const size_t indices = resident->indices.capacity() *
		sizeof(uint32_t);
	    const size_t normals = resident->cornerNormals.capacity() *
		sizeof(SbVec3f);
	    if (positions > SIZE_MAX - bytes ||
		indices > SIZE_MAX - bytes - positions ||
		normals > SIZE_MAX - bytes - positions - indices)
		return SIZE_MAX;
	    bytes += positions + indices + normals;
	}
	return bytes;
    }
    const Obol::TriMesh *mesh =
	generation && generation->shadedGeometry &&
	generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!mesh)
	return 0;
    return mesh->positions.capacity() * sizeof(SbVec3f) +
	mesh->normals.capacity() * sizeof(SbVec3f) +
	mesh->indices.capacity() * sizeof(uint32_t);
}

SbBox3f
BObolLodProgressiveMesh::bounds(void) const
{
    if (!this->p)
	return SbBox3f();
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->bounds : SbBox3f();
}

SbVec3f
BObolLodProgressiveMesh::quantizationMinimum(void) const
{
    if (!this->p)
	return SbVec3f(0.0f, 0.0f, 0.0f);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->quantizationMinimum :
	SbVec3f(0.0f, 0.0f, 0.0f);
}

SbVec3f
BObolLodProgressiveMesh::quantizationMaximum(void) const
{
    if (!this->p)
	return SbVec3f(0.0f, 0.0f, 0.0f);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->quantizationMaximum :
	SbVec3f(0.0f, 0.0f, 0.0f);
}

SbBool
BObolLodProgressiveMesh::cullBackfaces(void) const
{
    if (!this->p)
	return FALSE;
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    return generation ? generation->shadedCullBackfaces : FALSE;
}

struct BObolPreparedCornerKey {
    uint32_t position = 0;
    uint32_t normalX = 0;
    uint32_t normalY = 0;
    uint32_t normalZ = 0;

    bool operator==(const BObolPreparedCornerKey &other) const
    {
	return position == other.position &&
	    normalX == other.normalX &&
	    normalY == other.normalY &&
	    normalZ == other.normalZ;
    }
};

struct BObolPreparedCornerHash {
    size_t operator()(const BObolPreparedCornerKey &key) const
    {
	size_t hash = static_cast<size_t>(key.position);
	const auto combine = [&hash](uint32_t word) {
	    hash ^= static_cast<size_t>(word) +
		static_cast<size_t>(0x9e3779b9u) +
		(hash << 6) + (hash >> 2);
	};
	combine(key.normalX);
	combine(key.normalY);
	combine(key.normalZ);
	return hash;
    }
};

static uint32_t
bobol_prepared_float_bits(float value)
{
    uint32_t bits = 0;
    memcpy(&bits, &value, sizeof(bits));
    return bits;
}

/*
 * BRL-CAD PoP payloads retain authored normals per triangle corner, while
 * Obol binds one normal to each indexed vertex.  Canonicalize that mismatch
 * on the LoD worker which owns this immutable generation.  Publishing the
 * resulting PartGeometry by shared pointer prevents the GUI thread from
 * copying and splitting a large mesh during scene presentation.
 */
static bool
bobol_prepare_authored_corner_normals(
    Obol::TriMesh &mesh, const std::vector<SbVec3f> &cornerNormals)
{
    if (cornerNormals.empty())
	return true;
    if (cornerNormals.size() != mesh.indices.size())
	return false;

    std::unordered_map<BObolPreparedCornerKey, uint32_t,
	BObolPreparedCornerHash> vertexByCorner;
    vertexByCorner.reserve(mesh.indices.size());
    std::vector<SbVec3f> positions;
    std::vector<SbVec3f> normals;
    std::vector<uint32_t> indices;
    positions.reserve(mesh.indices.size());
    normals.reserve(mesh.indices.size());
    indices.reserve(mesh.indices.size());
    for (size_t triangle = 0;
	 triangle + 2 < mesh.indices.size(); triangle += 3) {
	const uint32_t ia = mesh.indices[triangle];
	const uint32_t ib = mesh.indices[triangle + 1];
	const uint32_t ic = mesh.indices[triangle + 2];
	if (ia >= mesh.positions.size() ||
	    ib >= mesh.positions.size() ||
	    ic >= mesh.positions.size())
	    return false;
	SbVec3f faceNormal =
	    (mesh.positions[ib] - mesh.positions[ia]).cross(
		mesh.positions[ic] - mesh.positions[ia]);
	const bool validFaceNormal = faceNormal.sqrLength() > 0.0f;
	if (validFaceNormal)
	    faceNormal.normalize();
	else
	    faceNormal.setValue(0.0f, 0.0f, 1.0f);

	std::array<SbVec3f, 3> triangleNormals;
	float normalDot = 0.0f;
	for (size_t corner = 0; corner < 3; ++corner) {
	    triangleNormals[corner] = cornerNormals[triangle + corner];
	    if (triangleNormals[corner].sqrLength() > 0.0f) {
		triangleNormals[corner].normalize();
		if (validFaceNormal)
		    normalDot += triangleNormals[corner].dot(faceNormal);
	    } else {
		triangleNormals[corner] = faceNormal;
	    }
	}
	if (validFaceNormal && normalDot < 0.0f) {
	    for (SbVec3f &normal : triangleNormals)
		normal.negate();
	}

	for (size_t corner = 0; corner < 3; ++corner) {
	    const size_t sourceCorner = triangle + corner;
	    const uint32_t sourceIndex = mesh.indices[sourceCorner];
	    const SbVec3f &normal = triangleNormals[corner];
	    const BObolPreparedCornerKey key = {
		sourceIndex,
		bobol_prepared_float_bits(normal[0]),
		bobol_prepared_float_bits(normal[1]),
		bobol_prepared_float_bits(normal[2])
	    };
	    auto found = vertexByCorner.find(key);
	    if (found == vertexByCorner.end()) {
		if (positions.size() >
		    static_cast<size_t>(std::numeric_limits<uint32_t>::max()))
		    return false;
		const uint32_t targetIndex =
		    static_cast<uint32_t>(positions.size());
		positions.push_back(mesh.positions[sourceIndex]);
		normals.push_back(normal);
		found = vertexByCorner.emplace(key, targetIndex).first;
	    }
	    indices.push_back(found->second);
	}
    }
    mesh.positions = std::move(positions);
    mesh.normals = std::move(normals);
    mesh.indices = std::move(indices);
    return true;
}

static constexpr uint64_t lod_spatial_page_lineage_salt =
    UINT64_C(0x9e3779b97f4a7c15);

static uint64_t
bobol_spatial_page_lineage(uint64_t assetLineage, uint32_t chunkId)
{
    uint64_t lineage = assetLineage ^
	(lod_spatial_page_lineage_salt + static_cast<uint64_t>(chunkId) +
	 (assetLineage << 6) + (assetLineage >> 2));
    return lineage ? lineage : lod_spatial_page_lineage_salt;
}

static std::shared_ptr<const Obol::PartGeometry>
bobol_prepare_spatial_page_shaded(
    const BObolLodProgressiveMeshGeneration &generation,
    const BObolLodProgressiveMeshGeneration::ResidentChunk &resident,
    uint64_t assetLineage)
{
    static constexpr int shadedChannelKey = 2;
    {
	std::lock_guard<std::mutex> lock(resident.preparedMutex);
	const auto cached = resident.preparedCadGeometry.find(shadedChannelKey);
	if (cached != resident.preparedCadGeometry.end()) {
	    std::shared_ptr<const Obol::PartGeometry> geometry =
		cached->second.lock();
	    if (geometry)
		return geometry;
	}
    }

    if (resident.chunkId >= generation.chunks.size() ||
	resident.positions.empty() || resident.indices.size() < 3 ||
	resident.indices.size() % 3 != 0)
	return std::shared_ptr<const Obol::PartGeometry>();
    const auto &metadata = generation.chunks[resident.chunkId];
    std::shared_ptr<Obol::PartGeometryBuilder> geometry(new Obol::PartGeometryBuilder);
    Obol::TriMesh mesh;
    mesh.positions = resident.positions;
    mesh.indices = resident.indices;
    mesh.bounds = resident.bounds;
    mesh.progressiveMinimumCut = static_cast<uint8_t>(
	std::max(0, metadata.minimumCut));
    int drawableCut = resident.residentCut;
    for (int candidate = resident.residentCut + 1;
	    candidate <= metadata.maximumCut; ++candidate) {
	const auto &candidateInfo = metadata.cuts[candidate];
	if (candidateInfo.point_count > resident.positions.size() ||
	    static_cast<uint64_t>(candidateInfo.face_count) * 3u >
		resident.indices.size())
	    break;
	drawableCut = candidate;
    }
    mesh.progressiveResidentCut = static_cast<uint8_t>(
	std::max(0, drawableCut));
    mesh.progressiveQuantizationMinimum = generation.quantizationMinimum;
    mesh.progressiveQuantizationMaximum = generation.quantizationMaximum;
    mesh.progressiveLineage = bobol_spatial_page_lineage(
	assetLineage, resident.chunkId);
    mesh.progressiveCuts.resize(generation.cuts.size());
    for (size_t cut = 0; cut < generation.cuts.size(); ++cut) {
	mesh.progressiveCuts[cut].indexCount = static_cast<uint32_t>(
	    std::min<uint64_t>(
		static_cast<uint64_t>(metadata.cuts[cut].face_count) * 3u,
		mesh.indices.size()));
	mesh.progressiveCuts[cut].positionCount = static_cast<uint32_t>(
	    std::min<size_t>(metadata.cuts[cut].point_count,
		mesh.positions.size()));
	mesh.progressiveCuts[cut].quantization = {
	    generation.cuts[cut].quantization_bits[X],
	    generation.cuts[cut].quantization_bits[Y],
	    generation.cuts[cut].quantization_bits[Z]
	};
    }
    if (!resident.cornerNormals.empty()) {
	if (!bobol_prepare_authored_corner_normals(
		mesh, resident.cornerNormals))
	    return std::shared_ptr<const Obol::PartGeometry>();
	size_t scanned = 0;
	uint32_t maximum = 0;
	for (Obol::ProgressiveTriangleCut &cut : mesh.progressiveCuts) {
	    const size_t end = std::min<size_t>(
		cut.indexCount, mesh.indices.size());
	    for (; scanned < end; ++scanned)
		maximum = std::max(maximum, mesh.indices[scanned]);
	    cut.positionCount = end ? maximum + 1u : 0u;
	}
	/* Vertex splitting changes the array prefix and cannot inherit an
	 * append-only promise from the source corner-normal stream. */
	mesh.progressiveLineage = 0;
    }
    geometry->shaded = std::move(mesh);
    geometry->conservativeBounds = resident.bounds;
    geometry->shadedCullBackfaces =
	generation.shadedCullBackfaces ? true : false;
    /* A page is a private storage partition of one logical surface, not a
     * semantic object.  Replacing an isolated page by a point can punch a
     * visible pinhole into an otherwise pixel-exact mesh.  Whole-occurrence
     * aggregation remains responsible for subpixel-object proxies. */
    geometry->subpixelProxyEligible = false;
    const std::shared_ptr<const Obol::PartGeometry> completed =
	bobol_cad_build_geometry(std::move(*geometry),
	    "spatial shaded presentation page");
    if (!completed)
	return std::shared_ptr<const Obol::PartGeometry>();
    {
	std::lock_guard<std::mutex> lock(resident.preparedMutex);
	resident.preparedCadGeometry[shadedChannelKey] = completed;
    }
    return completed;
}

static std::shared_ptr<const Obol::PartGeometry>
bobol_prepare_spatial_page_geometry(
    const BObolLodProgressiveMeshGeneration &generation,
    const BObolLodProgressiveMeshGeneration::ResidentChunk &resident,
    int drawMode, uint64_t assetLineage)
{
    const bool wire = drawMode == BOBOL_LOD_DRAW_WIRE ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    const bool shaded = drawMode == BOBOL_LOD_DRAW_SHADED ||
	drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    if (!wire && !shaded)
	return std::shared_ptr<const Obol::PartGeometry>();
    const int channelKey = (wire ? 1 : 0) | (shaded ? 2 : 0);
    if (channelKey == 2)
	return bobol_prepare_spatial_page_shaded(
	    generation, resident, assetLineage);
    {
	std::lock_guard<std::mutex> lock(resident.preparedMutex);
	const auto cached = resident.preparedCadGeometry.find(channelKey);
	if (cached != resident.preparedCadGeometry.end()) {
	    std::shared_ptr<const Obol::PartGeometry> geometry =
		cached->second.lock();
	    if (geometry)
		return geometry;
	}
    }

    const std::shared_ptr<const Obol::PartGeometry> shadedGeometry =
	bobol_prepare_spatial_page_shaded(
	    generation, resident, assetLineage);
    const Obol::TriMesh *sourceMesh =
	shadedGeometry && shadedGeometry->shaded ?
	    &*shadedGeometry->shaded : NULL;
    if (!sourceMesh)
	return std::shared_ptr<const Obol::PartGeometry>();

    std::shared_ptr<Obol::PartGeometryBuilder> geometry(new Obol::PartGeometryBuilder);
    Obol::WireRep wireRep;
    wireRep.bounds = resident.bounds;
    wireRep.triangleEdgeGeometry = shadedGeometry;
    wireRep.triangleEdgeSegmentCount = sourceMesh->indices.size();
    wireRep.progressiveMinimumCut = sourceMesh->progressiveMinimumCut;
    wireRep.progressiveResidentCut = sourceMesh->progressiveResidentCut;
    wireRep.progressiveQuantizationMinimum =
	sourceMesh->progressiveQuantizationMinimum;
    wireRep.progressiveQuantizationMaximum =
	sourceMesh->progressiveQuantizationMaximum;
    wireRep.progressiveLineage = sourceMesh->progressiveLineage;
    wireRep.progressiveCuts.resize(sourceMesh->progressiveCuts.size());
    for (size_t cut = 0; cut < sourceMesh->progressiveCuts.size(); ++cut) {
	wireRep.progressiveCuts[cut].segmentFirst = 0;
	wireRep.progressiveCuts[cut].segmentCount =
	    sourceMesh->progressiveCuts[cut].indexCount;
	wireRep.progressiveCuts[cut].quantization =
	    sourceMesh->progressiveCuts[cut].quantization;
    }
    geometry->wire = std::move(wireRep);
    if (shaded)
	geometry->shaded = *sourceMesh;
    geometry->conservativeBounds = resident.bounds;
    geometry->shadedCullBackfaces =
	generation.shadedCullBackfaces ? true : false;
    geometry->subpixelProxyEligible = false;
    const std::shared_ptr<const Obol::PartGeometry> completed =
	bobol_cad_build_geometry(std::move(*geometry),
	    "spatial wire presentation page");
    if (!completed)
	return std::shared_ptr<const Obol::PartGeometry>();
    {
	std::lock_guard<std::mutex> lock(resident.preparedMutex);
	resident.preparedCadGeometry[channelKey] = completed;
    }
    return completed;
}

SbBool
BObolLodProgressiveMesh::prepareCadPresentationLayers(
    int drawMode, const std::vector<uint32_t> &chunkIds, int activeCut,
    std::vector<BObolLodPresentationLayer> &layers) const
{
    layers.clear();
    const bool traceFailures =
	getenv("BOBOL_LOD_TRACE_PRESENTATION_LAYERS") != NULL;
    const auto traceFailure = [traceFailures, drawMode, activeCut](
	const char *reason, uint32_t chunkId, uint64_t faceCount,
	uint64_t pointCount, size_t residentPoints, size_t residentIndices,
	int residentCut) {
	if (!traceFailures)
	    return;
	bu_log("BObol spatial presentation preparation failure reason=%s "
	       "chunk=%u draw=%d active=%d faces=%llu points=%llu "
	       "resident_points=%zu resident_indices=%zu resident_cut=%d\n",
	       reason, chunkId, drawMode, activeCut,
	       static_cast<unsigned long long>(faceCount),
	       static_cast<unsigned long long>(pointCount), residentPoints,
	       residentIndices, residentCut);
    };
    if (!this->p || chunkIds.empty()) {
	traceFailure("invalid-input", 0, 0, 0, 0, 0, -1);
	return FALSE;
	}
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || generation->chunks.empty()) {
	traceFailure("missing-spatial-generation", 0, 0, 0, 0, 0, -1);
	return FALSE;
	}
    if (activeCut < generation->minimumCut ||
	activeCut > generation->maximumCut) {
	traceFailure("invalid-active-cut", 0, 0, 0, 0, 0, -1);
	return FALSE;
	}
    for (size_t selected = 0; selected < chunkIds.size(); ++selected) {
	const uint32_t chunkId = chunkIds[selected];
	if (chunkId >= generation->chunks.size() ||
	    (selected && chunkId <= chunkIds[selected - 1])) {
	    traceFailure("invalid-page-set", chunkId, 0, 0, 0, 0, -1);
	    return FALSE;
	}
	const auto &cut = generation->chunks[chunkId].cuts[activeCut];
	if (!cut.face_count)
	    continue;
	const auto &resident = generation->residentChunks[chunkId];
	if (!resident || cut.point_count > resident->positions.size() ||
	    static_cast<uint64_t>(cut.face_count) * 3u >
		resident->indices.size()) {
	    traceFailure("incomplete-resident-page", chunkId,
		cut.face_count, cut.point_count,
		resident ? resident->positions.size() : 0,
		resident ? resident->indices.size() : 0,
		resident ? resident->residentCut : -1);
	    return FALSE;
	}
    }
    layers.reserve(chunkIds.size());
    for (uint32_t chunkId : chunkIds) {
	const auto &resident = generation->residentChunks[chunkId];
	if (!resident ||
	    !generation->chunks[chunkId].cuts[activeCut].face_count)
	    continue;
	BObolLodPresentationLayer layer;
	std::string key("page:");
	key += std::to_string(chunkId);
	layer.layerKey = key.c_str();
	layer.geometry = bobol_prepare_spatial_page_geometry(
	    *generation, *resident, drawMode, this->p->lineage);
	if (!layer.geometry) {
	    const auto &cut = generation->chunks[chunkId].cuts[activeCut];
	    traceFailure("renderer-geometry-rejected", chunkId,
		cut.face_count, cut.point_count, resident->positions.size(),
		resident->indices.size(), resident->residentCut);
	    layers.clear();
	    return FALSE;
	}
	layer.geometryRevision = resident->geometryRevision ?
	    resident->geometryRevision : generation->revision;
	layer.activeCut = activeCut;
	layer.coverage = FALSE;
	layers.push_back(std::move(layer));
    }
    /* Spatial pages are private storage, not independent surfaces.  At a
     * coarse global PoP cut every face carried by a view-local page set may
     * legitimately have collapsed into other pages.  An empty selection is
     * therefore a valid zero-draw result, not corrupt renderer data. */
    return TRUE;
}

std::shared_ptr<const Obol::PartGeometry>
BObolLodProgressiveMesh::prepareCadGeometry(
    int drawMode, uint64_t *preparedRevision) const
{
    if (preparedRevision)
	*preparedRevision = 0;
    if (!this->p)
	return std::shared_ptr<const Obol::PartGeometry>();

    const bool wire =
	drawMode == BOBOL_LOD_DRAW_WIRE ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    const bool shaded =
	drawMode == BOBOL_LOD_DRAW_SHADED ||
	drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE;
    if (!wire && !shaded)
	return std::shared_ptr<const Obol::PartGeometry>();
    const int channelKey = (wire ? 1 : 0) | (shaded ? 2 : 0);

    /*
     * Retain one immutable source generation while constructing the renderer
     * object.  A concurrent refinement can publish its replacement without
     * invalidating these arrays or blocking metadata reads on the GUI thread.
     */
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    const Obol::TriMesh *sourceMesh =
	generation && generation->shadedGeometry &&
	generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!generation || generation->residentCut < 0 || !sourceMesh ||
	sourceMesh->positions.empty() || sourceMesh->indices.size() < 3 ||
	sourceMesh->indices.size() % 3 != 0)
	return std::shared_ptr<const Obol::PartGeometry>();

    /*
     * Shaded drawing is the common path.  The worker constructed this exact
     * immutable renderer record while loading the PoP prefix, so publication
     * is one shared-pointer handoff with no array conversion or copy.
     */
    if (shaded && !wire) {
	if (preparedRevision)
	    *preparedRevision = generation->revision;
	return generation->shadedGeometry;
    }

    {
	std::lock_guard<std::mutex> lock(generation->preparedMutex);
	const auto cached =
	    generation->preparedCadGeometry.find(channelKey);
	if (cached != generation->preparedCadGeometry.end()) {
	    std::shared_ptr<const Obol::PartGeometry> geometry =
		cached->second.lock();
	    if (geometry) {
		if (preparedRevision)
		    *preparedRevision = generation->revision;
		return geometry;
	    }
	}
    }

    const int drawableCut = sourceMesh->isProgressive() ?
	static_cast<int>(sourceMesh->progressiveResidentCut) :
	generation->residentCut;

    std::shared_ptr<Obol::PartGeometryBuilder> geometry(
	new Obol::PartGeometryBuilder);
    const SbBox3f bounds(
	generation->quantizationMinimum,
	generation->quantizationMaximum);
    /* Hidden-line requests need both channels.  Wire drawing aliases the
     * immutable indexed triangle stream instead of expanding each face into
     * six copied endpoint positions.  Plain shaded drawing never takes this
     * composition path. */
    if (wire) {
	Obol::WireRep wireRep;
	wireRep.bounds = bounds;
	wireRep.triangleEdgeGeometry = generation->shadedGeometry;
	wireRep.triangleEdgeSegmentCount = sourceMesh->indices.size();
	wireRep.progressiveMinimumCut =
	    static_cast<uint8_t>(std::max(0, generation->minimumCut));
	wireRep.progressiveResidentCut =
	    static_cast<uint8_t>(std::max(0, drawableCut));
	wireRep.progressiveQuantizationMinimum =
	    generation->quantizationMinimum;
	wireRep.progressiveQuantizationMaximum =
	    generation->quantizationMaximum;
	/* Triangle-edge wire records are emitted in cumulative triangle-prefix
	 * order.  They therefore inherit the source PoP stream's append-only
	 * identity; native curve LoD records do not make this promise. */
    wireRep.progressiveLineage = sourceMesh->progressiveLineage;
    wireRep.progressiveCuts.resize(sourceMesh->progressiveCuts.size());
    for (size_t cut = 0; cut < wireRep.progressiveCuts.size(); ++cut) {
	wireRep.progressiveCuts[cut].segmentFirst = 0;
	/* Each triangle index denotes one derived wire segment.  The admitted
	 * source cut is the canonical global prefix.  Reconstructing that count
	 * from hierarchy totals or spatial ranges is unsound: compacted arrays may
	 * expose coordinate-only cuts, and spatial ranges need not form a disjoint
	 * global prefix. */
	wireRep.progressiveCuts[cut].segmentCount =
	    sourceMesh->progressiveCuts[cut].indexCount;
	wireRep.progressiveCuts[cut].quantization =
		sourceMesh->progressiveCuts[cut].quantization;
	}
	if (sourceMesh->hasProgressiveClusters()) {
	    wireRep.progressiveClusters.reserve(
		sourceMesh->progressiveClusters.size());
	    for (const Obol::ProgressiveTriangleCluster &sourceCluster :
		    sourceMesh->progressiveClusters) {
		Obol::ProgressiveWireCluster targetCluster;
		targetCluster.bounds = sourceCluster.bounds;
		targetCluster.residentCut = sourceCluster.residentCut;
		targetCluster.ranges.reserve(sourceCluster.ranges.size());
		for (const Obol::ProgressiveTriangleClusterRange &sourceRange :
			sourceCluster.ranges) {
		    /* Three triangle indices produce three wire segments, so the
		     * source index offset and count are numerically identical to
		     * the corresponding segment offset and count. */
		    targetCluster.ranges.push_back({
			sourceRange.firstIndex, sourceRange.indexCount,
			sourceRange.activationCut});
		}
		wireRep.progressiveClusters.push_back(
		    std::move(targetCluster));
	    }
	    wireRep.progressiveClusterGridResolution =
		sourceMesh->progressiveClusterGridResolution;
	}
	geometry->wire = std::move(wireRep);
    }

    if (shaded) {
	geometry->shaded = *sourceMesh;
	geometry->shadedCullBackfaces =
	    generation->shadedCullBackfaces ? true : false;
    }
    geometry->subpixelProxyEligible = true;
    geometry->aggregateProxyCorners = generation->aggregateProxyCorners;

    std::shared_ptr<const Obol::PartGeometry> prepared =
	bobol_cad_build_geometry_with_optional_proxy(
	    std::move(*geometry), "combined PoP presentation");
    if (!prepared)
	return std::shared_ptr<const Obol::PartGeometry>();
    {
	std::lock_guard<std::mutex> lock(generation->preparedMutex);
	generation->preparedCadGeometry[channelKey] = prepared;
    }
    if (preparedRevision)
	*preparedRevision = generation->revision;
    return prepared;
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
    sourceRoutingId = 0;
    sourcePopulationEpoch = 0;
    coalesceAssetProducer = FALSE;
    sourceEntryIndex = UINT32_MAX;
    viewRevision = 0;
    policyRevision = 0;
    visualEmphasis = 0;
    drawMode = BOBOL_LOD_DRAW_UNKNOWN;
    lodPolicy = 0;
    providerId = "";
    providerVersion = "";
    qualityTier = BOBOL_LOD_QUALITY_METADATA;
    projectedPixelDiameter = 0.0f;
    projectedPixelArea = 0.0f;
    projectedPixelPerimeter = 0.0f;
    projectedBoundsContained = FALSE;
    targetPixelError = 1.0f;
    requestedCut = -1;
    localToRoot.makeIdentity();
    viewProjection.makeIdentity();
    spatialProjectionValid = FALSE;
    requiredChunks.clear();
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

BObolLodPresentationLayer::BObolLodPresentationLayer(void) :
    geometryRevision(0),
    activeCut(-1),
    coverage(FALSE)
{
}

SbBool
BObolLodPresentationLayer::isValid(void) const
{
    return layerKey.getLength() > 0 && geometry ? TRUE : FALSE;
}

size_t
BObolLodPresentationLayer::estimateBytes(void) const
{
    /* PartGeometry is shared immutable storage.  Its owner accounts for its
     * arrays; this result only owns the layer descriptor and its key. */
    return sizeof(*this) + static_cast<size_t>(layerKey.getLength()) + 1u;
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
    preparedCadGeometry.reset();
    preparedCadGeometryRevision = 0;
    presentationLayers.clear();
    resolvedCut = -1;
    residentCut = -1;
    presentationAdmissionCertified = FALSE;
    presentationAdmissionViewRevision = 0;
    presentationAdmissionPolicyRevision = 0;
    presentationAdmissionCut = -1;
    residentAdmissionRevision = 0;
    payloadKind = BOBOL_LOD_PAYLOAD_NONE;
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
    memoryLimited = FALSE;
    diagnostic = "";
}

void
BObolLodResult::canonicalizePayload(void)
{
    if (providerStatus != BOBOL_LOD_PROVIDER_READY &&
	providerStatus != BOBOL_LOD_PROVIDER_FALLBACK) {
	payloadKind = BOBOL_LOD_PAYLOAD_STATUS;
    } else {
	switch (resultKind) {
	    case BOBOL_LOD_RESULT_DIRECTORY:
		payloadKind = BOBOL_LOD_PAYLOAD_DIRECTORY;
		break;
	    case BOBOL_LOD_RESULT_ATTRIBUTES:
		payloadKind = BOBOL_LOD_PAYLOAD_ATTRIBUTES;
		break;
	    case BOBOL_LOD_RESULT_AABB:
	    case BOBOL_LOD_RESULT_PROXY:
		payloadKind = BOBOL_LOD_PAYLOAD_PROXY;
		break;
	    case BOBOL_LOD_RESULT_MESH:
	    case BOBOL_LOD_RESULT_FULL_DETAIL:
		payloadKind = BOBOL_LOD_PAYLOAD_MESH;
		break;
	    case BOBOL_LOD_RESULT_NONE:
		payloadKind = BOBOL_LOD_PAYLOAD_NONE;
		break;
	    case BOBOL_LOD_RESULT_DIAGNOSTIC:
	    default:
		payloadKind = BOBOL_LOD_PAYLOAD_STATUS;
		break;
	}
    }

    if (payloadKind != BOBOL_LOD_PAYLOAD_DIRECTORY)
	dependencies.clear();
    if (payloadKind != BOBOL_LOD_PAYLOAD_ATTRIBUTES)
	attributes.clear();
    if (payloadKind != BOBOL_LOD_PAYLOAD_PROXY)
	proxy.clear();
    if (payloadKind != BOBOL_LOD_PAYLOAD_MESH) {
	mesh.clear();
	progressiveMesh.reset();
	preparedCadGeometry.reset();
	preparedCadGeometryRevision = 0;
	presentationLayers.clear();
	resolvedCut = -1;
	residentCut = -1;
	residentAdmissionRevision = 0;
	hasSnappedPoints = FALSE;
	hasNormals = FALSE;
	shadedCullBackfaces = FALSE;
	memoryLimited = FALSE;
    }
    if (payloadKind != BOBOL_LOD_PAYLOAD_PROXY &&
	payloadKind != BOBOL_LOD_PAYLOAD_MESH)
	geometry.clear();
}

SbBool
BObolLodResult::payloadIsConsistent(void) const
{
    bool layersValid = true;
    std::unordered_set<std::string> layerKeys;
    for (const BObolLodPresentationLayer &layer : presentationLayers) {
	if (!layer.isValid() || !layerKeys.emplace(
		layer.layerKey.getString()).second) {
	    layersValid = false;
	    break;
	}
    }
    const bool hasMesh = mesh.isValid() ||
	(progressiveMesh && progressiveMesh->isValid()) ||
	preparedCadGeometry ||
	(!presentationLayers.empty() && layersValid) ||
	(geometry.isValid() &&
	    (resultKind == BOBOL_LOD_RESULT_MESH ||
	     resultKind == BOBOL_LOD_RESULT_FULL_DETAIL));
    const bool hasProxy = proxy.isValid();
    switch (payloadKind) {
	case BOBOL_LOD_PAYLOAD_NONE:
	    return !hasMesh && !hasProxy && dependencies.empty() &&
		attributes.empty() ? TRUE : FALSE;
	case BOBOL_LOD_PAYLOAD_DIRECTORY:
	    return !hasMesh && !hasProxy && attributes.empty() ? TRUE : FALSE;
	case BOBOL_LOD_PAYLOAD_ATTRIBUTES:
	    return !hasMesh && !hasProxy && dependencies.empty() ? TRUE : FALSE;
	case BOBOL_LOD_PAYLOAD_PROXY:
	    return hasProxy && !hasMesh && dependencies.empty() &&
		attributes.empty() ? TRUE : FALSE;
	case BOBOL_LOD_PAYLOAD_MESH:
	    return hasMesh && !hasProxy && dependencies.empty() &&
		attributes.empty() ? TRUE : FALSE;
	case BOBOL_LOD_PAYLOAD_STATUS:
	    return !hasMesh && !hasProxy && dependencies.empty() &&
		attributes.empty() ? TRUE : FALSE;
	default:
	    return FALSE;
    }
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

/*
 * Geometry asset keys are constructed in the per-occurrence submission hot
 * path.  std::ostringstream performs locale setup and repeatedly allocates
 * while formatting even this integer/string-only key; that cost was visible
 * above projection in 50k-leaf GUI profiles.  Keep the byte format identical
 * while appending directly into one reserved string.
 */
template <typename Integer>
static void
append_integer_string(std::string &out, Integer value)
{
    char buffer[3 * sizeof(Integer) + 4];
    const std::to_chars_result converted =
	std::to_chars(buffer, buffer + sizeof(buffer), value);
    if (converted.ec == std::errc())
	out.append(buffer, converted.ptr);
}

static void
append_string_field(std::string &out, const char *name,
		    const SbString &value)
{
    const char *str = value.getString();
    const size_t len = str ? strlen(str) : 0;
    out.append(name);
    out.push_back('=');
    append_integer_string(out, len);
    out.push_back(':');
    if (len)
	out.append(str, len);
    out.push_back(';');
}

static void
append_uint_field(std::string &out, const char *name, uint64_t value)
{
    out.append(name);
    out.push_back('=');
    append_integer_string(out, value);
    out.push_back(';');
}

static void
append_int_field(std::string &out, const char *name, int value)
{
    out.append(name);
    out.push_back('=');
    append_integer_string(out, value);
    out.push_back(';');
}

static void
append_float_string(std::string &out, float value)
{
    char buffer[64];
    const std::to_chars_result converted = std::to_chars(
	buffer, buffer + sizeof(buffer), value, std::chars_format::general, 9);
    if (converted.ec == std::errc())
	out.append(buffer, converted.ptr);
}

static void
append_bounds_field(std::string &out, const char *name,
		    const SbBox3f &bounds)
{
    out.append(name);
    out.push_back('=');
    if (bounds.isEmpty()) {
	out.append("empty;");
	return;
    }

    const SbVec3f &bmin = bounds.getMin();
    const SbVec3f &bmax = bounds.getMax();
    const float components[6] = {
	bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2]
    };
    for (size_t i = 0; i < 6; ++i) {
	if (i)
	    out.push_back(',');
	append_float_string(out, components[i]);
    }
    out.push_back(';');
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
    std::string out;
    out.reserve(512u +
	request.databaseId.getLength() +
	request.objectPath.getLength() +
	request.objectName.getLength() +
	request.occurrenceKey.getLength() +
	request.providerId.getLength() +
	request.providerVersion.getLength());

    out.append("bobol-lod-v3;");
    append_string_field(out, "database_id", request.databaseId);
    append_uint_field(out, "database_revision",
	request.databaseRevision.value());
    append_uint_field(out, "source_revision", request.sourceRevision.value());
    append_uint_field(out, "source_content_hash", request.sourceContentHash);
    append_string_field(out, "object_path", request.objectPath);
    append_string_field(out, "object_name", request.objectName);
    append_string_field(out, "occurrence_key", request.occurrenceKey);
    append_uint_field(out, "view_revision", request.viewRevision.value());
    append_uint_field(out, "policy_revision", request.policyRevision.value());
    append_int_field(out, "draw_mode", request.drawMode);
    append_uint_field(out, "lod_policy", request.lodPolicy);
    append_string_field(out, "provider_id", request.providerId);
    append_string_field(out, "provider_version", request.providerVersion);
    append_int_field(out, "quality_tier", request.qualityTier);
    /* The selected admissible cut determines the provider output.  Raw screen
     * measurements do not: including them would turn sub-pixel camera jitter
     * into distinct active/cache keys for the same payload. */
    append_int_field(out, "requested_cut", request.requestedCut);
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

    key.value = out.c_str();
    return key;
}

BObolLodCacheKey
bobol_lod_asset_cache_key(const BObolLodRequest &request)
{
    BObolLodCacheKey key;
    std::string out;
    out.reserve(256u +
	request.databaseId.getLength() +
	request.objectName.getLength() +
	request.objectPath.getLength() +
	request.providerId.getLength() +
	request.providerVersion.getLength());

    out.append("bobol-progressive-asset-v1;");
    append_string_field(out, "database_id", request.databaseId);
    append_string_field(out, "object_name", request.objectName);
    if (request.sourceContentHash != 0)
	append_uint_field(out, "source_content_hash",
	    request.sourceContentHash);
    else {
	append_uint_field(out, "database_revision",
	    request.databaseRevision.value());
	append_uint_field(out, "source_revision",
	    request.sourceRevision.value());
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

    key.value = out.c_str();
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
	const uint32_t idx = data.faces[i];
	if (idx > static_cast<uint32_t>(
		std::numeric_limits<int32_t>::max()) ||
	    static_cast<size_t>(idx) >= data.point_count) {
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
bobol_lod_request_keys_equal(const BObolLodRequest &lhs,
			     const BObolLodRequest &rhs)
{
    if (bu_strcmp(lhs.databaseId.getString(),
	    rhs.databaseId.getString()) != 0 ||
	lhs.databaseRevision != rhs.databaseRevision ||
	lhs.sourceRevision != rhs.sourceRevision ||
	lhs.sourceContentHash != rhs.sourceContentHash ||
	bu_strcmp(lhs.objectPath.getString(), rhs.objectPath.getString()) != 0 ||
	bu_strcmp(lhs.objectName.getString(), rhs.objectName.getString()) != 0 ||
	bu_strcmp(lhs.occurrenceKey.getString(),
	    rhs.occurrenceKey.getString()) != 0 ||
	lhs.sourceRoutingId != rhs.sourceRoutingId ||
	lhs.sourcePopulationEpoch != rhs.sourcePopulationEpoch ||
	lhs.viewRevision != rhs.viewRevision ||
	lhs.policyRevision != rhs.policyRevision ||
	lhs.drawMode != rhs.drawMode ||
	lhs.lodPolicy != rhs.lodPolicy ||
	bu_strcmp(lhs.providerId.getString(), rhs.providerId.getString()) != 0 ||
	bu_strcmp(lhs.providerVersion.getString(),
	    rhs.providerVersion.getString()) != 0 ||
	lhs.qualityTier != rhs.qualityTier ||
	lhs.requestedCut != rhs.requestedCut)
	return FALSE;

    if (lhs.bounds.isEmpty() != rhs.bounds.isEmpty())
	return FALSE;
    if (!lhs.bounds.isEmpty()) {
	const SbVec3f lhsBounds[2] = {
	    lhs.bounds.getMin(), lhs.bounds.getMax()
	};
	const SbVec3f rhsBounds[2] = {
	    rhs.bounds.getMin(), rhs.bounds.getMax()
	};
	if (memcmp(lhsBounds, rhsBounds, sizeof(lhsBounds)) != 0)
	    return FALSE;
    }

    const BObolLodCounts &a = lhs.sourceCounts;
    const BObolLodCounts &b = rhs.sourceCounts;
    if (a.faceCount != b.faceCount ||
	a.pointCount != b.pointCount ||
	a.originalPointCount != b.originalPointCount ||
	a.normalCount != b.normalCount ||
	a.lineCount != b.lineCount ||
	a.byteCount != b.byteCount ||
	lhs.providerParams.size() != rhs.providerParams.size())
	return FALSE;

    /*
     * Provider parameter lists are normally tiny.  A counted multiset
     * comparison avoids the sorted vector copy and heap traffic previously
     * paid for every completed result while retaining duplicate semantics.
     */
    for (size_t i = 0; i < lhs.providerParams.size(); ++i) {
	const BObolLodProviderParam &needle = lhs.providerParams[i];
	size_t lhsCount = 0;
	size_t rhsCount = 0;
	for (const BObolLodProviderParam &candidate : lhs.providerParams) {
	    if (bu_strcmp(needle.name.getString(),
		    candidate.name.getString()) == 0 &&
		bu_strcmp(needle.value.getString(),
		    candidate.value.getString()) == 0)
		++lhsCount;
	}
	for (const BObolLodProviderParam &candidate : rhs.providerParams) {
	    if (bu_strcmp(needle.name.getString(),
		    candidate.name.getString()) == 0 &&
		bu_strcmp(needle.value.getString(),
		    candidate.value.getString()) == 0)
		++rhsCount;
	}
	if (lhsCount != rhsCount)
	    return FALSE;
    }

    return TRUE;
}

SbBool
bobol_lod_result_matches_request(const BObolLodResult &result,
				   const BObolLodRequest &request)
{
    return bobol_lod_request_keys_equal(result.request, request);
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
    result.geometry.activeCut = info.active_cut;
    result.resolvedCut = request.requestedCut >= 0 ?
	request.requestedCut : info.active_cut;
    result.residentCut = info.active_cut;
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

    result.canonicalizePayload();
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
    result.canonicalizePayload();
    return result;
}

BObolLodResult
bobol_lod_attributes_result(const BObolLodRequest &request,
			      const std::vector<BObolLodAttribute> &attributes)
{
    BObolLodResult result = stage_result_base(request,
			      BOBOL_LOD_RESULT_ATTRIBUTES, BOBOL_LOD_QUALITY_ATTRIBUTES);

    result.attributes = attributes;
    result.canonicalizePayload();
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
    result.canonicalizePayload();
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

    result.canonicalizePayload();
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
