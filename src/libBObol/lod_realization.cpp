/*              L O D _ R E A L I Z A T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BLodRealization.h"

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
#include <string>
#include <string.h>
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
    int minimumCut = -1;
    int maximumCut = -1;
    int residentCut = -1;
    uint64_t revision = 0;
    SbBool shadedCullBackfaces = FALSE;
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
    uint64_t lineage = next.fetch_add(1, std::memory_order_relaxed);
    if (!lineage)
	lineage = next.fetch_add(1, std::memory_order_relaxed);
    return lineage ? lineage : 1;
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
bobol_lod_hierarchy_valid(
    const struct BObolMeshLodData &data,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int residentCut)
{
    if (hierarchy.min_cut < 0 ||
	hierarchy.max_cut < hierarchy.min_cut ||
	hierarchy.max_cut >= BOBOL_MESH_LOD_CUT_COUNT_MAX ||
	hierarchy.cut_count !=
	    static_cast<uint32_t>(hierarchy.max_cut + 1) ||
	residentCut < hierarchy.min_cut ||
	residentCut > hierarchy.max_cut ||
	hierarchy.resident_cut != residentCut ||
	!bobol_lod_valid_bounds(data.bmin, data.bmax) ||
	!bobol_lod_valid_bounds(hierarchy.quantization_min,
	    hierarchy.quantization_max))
	return false;

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
	    return false;
	priorPoints = points;
	priorFaces = faces;
	priorError = hierarchy.cuts[cut].object_error;
	for (int axis = 0; axis < 3; ++axis) {
	    const uint8_t bits =
		hierarchy.cuts[cut].quantization_bits[axis];
	    if (bits < 1 || bits > BOBOL_MESH_LOD_QUANTIZATION_BITS ||
		(cut > 0 && bits <
		 hierarchy.cuts[cut - 1].quantization_bits[axis]))
		return false;
	}
    }
    if (!hierarchy.cuts[hierarchy.max_cut].exact ||
	hierarchy.cuts[hierarchy.max_cut].object_error > 0.0)
	return false;
    return hierarchy.cuts[residentCut].point_count ==
	data.point_orig_count &&
	hierarchy.cuts[residentCut].face_count == data.face_count &&
	data.point_orig_count > 0 && data.face_count > 0;
}

static std::shared_ptr<BObolLodProgressiveMeshGeneration>
progressive_generation_from_data(
    const struct BObolMeshLodData &data,
    const struct BObolMeshLodHierarchyInfo &hierarchy,
    int residentCut, SbBool shadedCullBackfaces,
    uint64_t revision, uint64_t progressiveLineage)
{
    if (!bobol_lod_hierarchy_valid(data, hierarchy, residentCut))
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    std::shared_ptr<BObolLodProgressiveMeshGeneration> generation(
	new BObolLodProgressiveMeshGeneration);
    std::shared_ptr<Obol::PartGeometry> geometry(
	new Obol::PartGeometry);
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
    generation->shadedCullBackfaces = shadedCullBackfaces;
    generation->revision = revision ? revision : 1;

    mesh.positions.resize(data.point_orig_count);
    for (size_t i = 0; i < data.point_orig_count; ++i) {
	if (!bobol_lod_point_in_domain(data.points_orig[i],
		hierarchy.quantization_min,
		hierarchy.quantization_max))
	    return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	mesh.positions[i] = SbVec3f(
	    static_cast<float>(data.points_orig[i][X]),
	    static_cast<float>(data.points_orig[i][Y]),
	    static_cast<float>(data.points_orig[i][Z]));
    }

    if (data.face_count >
	    std::numeric_limits<size_t>::max() / 3)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
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
	    return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
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
		return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
	    cornerNormals[i] = SbVec3f(
		static_cast<float>(data.normals[i][X]),
		static_cast<float>(data.normals[i][Y]),
		static_cast<float>(data.normals[i][Z]));
	}
	if (!bobol_prepare_authored_corner_normals(
		mesh, cornerNormals))
	    return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
    }

    if (mesh.positions.empty() || mesh.indices.size() < 3 ||
	mesh.indices.size() % 3 != 0 ||
	(!mesh.normals.empty() &&
	 mesh.normals.size() != mesh.positions.size()))
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
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
    generation->shadedGeometry = geometry;
    return generation;
}

static std::shared_ptr<BObolLodProgressiveMeshGeneration>
progressive_generation_prefix(
    const BObolLodProgressiveMeshGeneration &source, int residentCut)
{
    const int cut = progressive_cut_clamp(&source, residentCut);
    if (cut < 0)
	return std::shared_ptr<BObolLodProgressiveMeshGeneration>();
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
    generation->minimumCut = source.minimumCut;
    generation->maximumCut = source.maximumCut;
    generation->residentCut = cut;
    generation->revision = source.revision + 1;
    if (!generation->revision)
	generation->revision = 1;
    generation->shadedCullBackfaces = source.shadedCullBackfaces;
    std::shared_ptr<Obol::PartGeometry> geometry(
	new Obol::PartGeometry);
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
	 candidate <= source.maximumCut; ++candidate) {
	if (sourceMesh->positionCountAtCut(
		static_cast<uint8_t>(candidate)) > mesh.positions.size() ||
	    sourceMesh->indexCountAtCut(
		static_cast<uint8_t>(candidate)) > mesh.indices.size())
	    break;
	drawableCut = candidate;
    }
    mesh.progressiveResidentCut =
        static_cast<uint8_t>(std::max(0, drawableCut));
    geometry->shaded = std::move(mesh);
    geometry->shadedCullBackfaces =
	source.shadedCullBackfaces ? true : false;
    geometry->subpixelProxyEligible = true;
    generation->shadedGeometry = geometry;
    return generation;
}

struct BObolLodProgressiveMeshTrim {
    const BObolLodProgressiveMeshPrivate *owner = NULL;
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> source;
    std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation;
};

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
    if (!this->p || !data.faces || !data.points_orig ||
	data.face_count == 0 || data.point_orig_count == 0 ||
	residentCut < hierarchy.min_cut ||
	residentCut > hierarchy.max_cut ||
	residentCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return FALSE;

    if (!bobol_lod_hierarchy_valid(data, hierarchy, residentCut))
	return FALSE;

    if (data.point_orig_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) ||
	data.face_count >
	    static_cast<size_t>(std::numeric_limits<int32_t>::max()) / 3)
	return FALSE;
    const size_t indexCount = data.face_count * 3;
    if (data.normals && data.normal_count != indexCount)
	return FALSE;

    std::lock_guard<std::mutex> updateLock(this->p->updateMutex);
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> prior =
	progressive_generation_load(this->p);
    uint64_t revision = prior ? prior->revision + 1 : 1;
    if (!revision)
	revision = 1;
    const std::shared_ptr<BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_from_data(data, hierarchy, residentCut,
	    shadedCullBackfaces, revision, this->p->lineage);
    if (!generation)
	return FALSE;
    progressive_generation_store(this->p, generation);
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
    generation->minimumCut = prior->minimumCut;
    generation->maximumCut = prior->maximumCut;
    generation->residentCut = residentCut;
    generation->revision = prior->revision + 1;
    if (!generation->revision)
	generation->revision = 1;
    generation->shadedCullBackfaces = prior->shadedCullBackfaces;

    std::shared_ptr<Obol::PartGeometry> geometry(
	new Obol::PartGeometry);
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
    generation->shadedGeometry = geometry;
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
    if (cut < 0)
	return FALSE;
    const Obol::TriMesh *mesh =
	generation->shadedGeometry && generation->shadedGeometry->shaded ?
	    &*generation->shadedGeometry->shaded : NULL;
    if (!mesh)
	return FALSE;
    const size_t points =
	mesh->positionCountAtCut(static_cast<uint8_t>(cut));
    const size_t indices =
	mesh->indexCountAtCut(static_cast<uint8_t>(cut));
    if (!points || points > mesh->positions.size() ||
	indices < 3 || indices > mesh->indices.size())
	return FALSE;
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
	    payload.clear();
	    return FALSE;
	}
	payload.coordIndex.push_back(static_cast<int32_t>(index));
	if (!mesh->normals.empty()) {
	    if (index >= mesh->normals.size()) {
		payload.clear();
		return FALSE;
	    }
	    payload.normals.push_back(mesh->normals[index]);
	}
    }
    return payload.isValid();
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
    return generation && generation->residentCut >= 0 && mesh &&
	!mesh->positions.empty() && mesh->indices.size() >= 3 &&
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
    if (!mesh)
	return FALSE;
    return mesh->isProgressive() &&
	requestedCut <= static_cast<int>(
	    mesh->progressiveResidentCut) ? TRUE : FALSE;
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
	if (!generation || generation->cuts.empty() ||
	    !std::isfinite(projectedPixelDiameter) ||
	    !std::isfinite(targetPixelError) ||
	    projectedPixelDiameter <= 0.0 || targetPixelError <= 0.0)
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

double
BObolLodProgressiveMesh::projectedErrorAtCut(
    int cut, double projectedPixelDiameter) const
{
    if (!this->p || cut < 0 || !std::isfinite(projectedPixelDiameter) ||
	projectedPixelDiameter <= 0.0)
	return std::numeric_limits<double>::infinity();
    const std::shared_ptr<const BObolLodProgressiveMeshGeneration> generation =
	progressive_generation_load(this->p);
    if (!generation || static_cast<size_t>(cut) >= generation->cuts.size())
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

    std::shared_ptr<Obol::PartGeometry> geometry(
	new Obol::PartGeometry);
    const SbBox3f bounds(
	generation->quantizationMinimum,
	generation->quantizationMaximum);
    /* Hidden-line requests need both channels.  Build the expanded wire
     * channel alongside a copy of the canonical shaded channel.  Plain
     * shaded drawing never takes this composition path. */
    if (wire) {
	Obol::WireRep wireRep;
	wireRep.bounds = bounds;
	if (sourceMesh->indices.size() >
		std::numeric_limits<size_t>::max() / 2)
	    return std::shared_ptr<const Obol::PartGeometry>();
	wireRep.segmentPoints.reserve(sourceMesh->indices.size() * 2);
	wireRep.segmentIds.reserve(sourceMesh->indices.size());
	uint32_t edgeId = 0;
	for (size_t i = 0; i + 2 < sourceMesh->indices.size(); i += 3) {
	    const uint32_t triangle[3] = {
		sourceMesh->indices[i],
		sourceMesh->indices[i + 1],
		sourceMesh->indices[i + 2]
	    };
	    for (size_t edge = 0; edge < 3; ++edge) {
		const uint32_t a = triangle[edge];
		const uint32_t b = triangle[(edge + 1) % 3];
		if (a >= sourceMesh->positions.size() ||
		    b >= sourceMesh->positions.size() ||
		    edgeId == std::numeric_limits<uint32_t>::max())
		    return std::shared_ptr<const Obol::PartGeometry>();
		wireRep.segmentPoints.push_back(sourceMesh->positions[a]);
		wireRep.segmentPoints.push_back(sourceMesh->positions[b]);
		wireRep.segmentIds.push_back(edgeId++);
	    }
	}
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
	    const size_t segments = generation->faceCount[cut] * 3;
	    wireRep.progressiveCuts[cut].segmentFirst = 0;
	    wireRep.progressiveCuts[cut].segmentCount =
		static_cast<uint32_t>(std::min(
		    segments,
		    static_cast<size_t>(
			std::numeric_limits<uint32_t>::max())));
	    wireRep.progressiveCuts[cut].quantization =
		sourceMesh->progressiveCuts[cut].quantization;
	}
	geometry->wire = std::move(wireRep);
    }

    if (shaded) {
	geometry->shaded = *sourceMesh;
	geometry->shadedCullBackfaces =
	    generation->shadedCullBackfaces ? true : false;
    }
    geometry->subpixelProxyEligible = true;

    std::shared_ptr<const Obol::PartGeometry> prepared = geometry;
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
    viewRevision = 0;
    policyRevision = 0;
    drawMode = BOBOL_LOD_DRAW_UNKNOWN;
    lodPolicy = 0;
    providerId = "";
    providerVersion = "";
    qualityTier = BOBOL_LOD_QUALITY_METADATA;
    projectedPixelDiameter = 0.0f;
    projectedPixelArea = 0.0f;
    projectedPixelPerimeter = 0.0f;
    targetPixelError = 1.0f;
    requestedCut = -1;
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
    preparedCadGeometry.reset();
    preparedCadGeometryRevision = 0;
    resolvedCut = -1;
    residentCut = -1;
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
    const bool hasMesh = mesh.isValid() ||
	(progressiveMesh && progressiveMesh->isValid()) ||
	preparedCadGeometry ||
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
