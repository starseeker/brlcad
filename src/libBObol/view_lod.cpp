/*                  V I E W _ L O D . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/str.h"

#include "BObol/BExportAction.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BMeasureAction.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BSnapAction.h"
#include "BObol/BViewLod.h"

#include <Obol/cad/SoCADAssembly.h>
#include <Obol/cad/SoCADViewState.h>

#include <Inventor/actions/SoAction.h>
#include <Inventor/actions/SoCallbackAction.h>
#include <Inventor/actions/SoGLRenderAction.h>
#include <Inventor/actions/SoGetBoundingBoxAction.h>
#include <Inventor/actions/SoPickAction.h>
#include <Inventor/actions/SoRayPickAction.h>
#include <Inventor/misc/SoState.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <cmath>
#include <new>
#include <string.h>
#include <utility>
#include <unordered_set>

SO_ELEMENT_SOURCE(SoBRLViewLodElement);
SO_NODE_SOURCE(SoBRLViewLodGroup);

static std::string
view_lod_string_key(const char *prefix, const SbString &value)
{
    if (!prefix || value.getLength() == 0)
	return std::string();

    std::string key(prefix);
    key += ':';
    key += value.getString();
    return key;
}

static std::string
view_lod_shape_primary_key(const SoBRLMeshShape *shape)
{
    if (!shape)
	return std::string();

    std::string key = view_lod_string_key("identity",
					  shape->sourceIdentity.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("path", shape->sourcePath.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("name", shape->sourceName.getValue());
    if (!key.empty())
	return key;

    char buffer[64] = {0};
    snprintf(buffer, sizeof(buffer), "node:%p", static_cast<const void *>(shape));
    return std::string(buffer);
}

static void
view_lod_append_unique_key(std::vector<std::string> &keys,
			   const std::string &key)
{
    if (key.empty())
	return;

    for (size_t i = 0; i < keys.size(); i++) {
	if (keys[i] == key)
	    return;
    }

    keys.push_back(key);
}

static void
view_lod_shape_keys(std::vector<std::string> &keys,
		    const SoBRLMeshShape *shape)
{
    if (!shape)
	return;

    view_lod_append_unique_key(keys, view_lod_string_key("identity",
			       shape->sourceIdentity.getValue()));
    view_lod_append_unique_key(keys, view_lod_string_key("path",
			       shape->sourcePath.getValue()));
    if (shape->sourcePath.getValue().getLength() > 1 &&
	shape->sourcePath.getValue().getString()[0] == '/')
	view_lod_append_unique_key(keys, std::string("path:") +
				   (shape->sourcePath.getValue().getString() + 1));
    view_lod_append_unique_key(keys, view_lod_string_key("name",
			       shape->sourceName.getValue()));
    view_lod_append_unique_key(keys, view_lod_shape_primary_key(shape));
}

static void
view_lod_result_keys(std::vector<std::string> &keys,
		     const BObolLodResult &result)
{
    /* A compact occurrence may share both object path and object name with
     * siblings.  Its immutable identity is consequently the only valid
     * display binding key. */
    if (result.request.occurrenceKey.getLength() > 0) {
	view_lod_append_unique_key(keys, view_lod_string_key("occurrence",
			       result.request.occurrenceKey));
	return;
    }
    view_lod_append_unique_key(keys, view_lod_string_key("path",
			       result.request.objectPath));
    if (result.request.objectPath.getLength() > 1 &&
	result.request.objectPath.getString()[0] == '/')
	view_lod_append_unique_key(keys, std::string("path:") +
				   (result.request.objectPath.getString() + 1));
    view_lod_append_unique_key(keys, view_lod_string_key("name",
			       result.request.objectName));
}

static const char *
view_lod_leaf_name_from_path(const SbString &path)
{
    const char *p = path.getString();
    if (!p || !p[0])
	return "";
    const char *slash = strrchr(p, '/');
    return (slash && slash[1]) ? slash + 1 : p;
}

static bool
view_lod_paths_equal(const SbString &a, const SbString &b)
{
    const char *ap = a.getString();
    const char *bp = b.getString();
    if (!ap || !bp)
	return false;
    if (ap[0] == '/')
	ap++;
    if (bp[0] == '/')
	bp++;
    return bu_strcmp(ap, bp) == 0;
}

static std::string
view_lod_source_primary_key(const SoBRLDatabaseSource *source)
{
    if (!source)
	return std::string();

    std::string key = view_lod_string_key("instance",
					  source->instanceKey.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("realization",
			      source->realizationIdentity.getValue());
    if (!key.empty())
	return key;

    key = view_lod_string_key("path", source->path.getValue());
    if (!key.empty())
	return key;

    char buffer[64] = {0};
    snprintf(buffer, sizeof(buffer), "source:%p",
	     static_cast<const void *>(source));
    return std::string(buffer);
}

static void
view_lod_source_keys(std::vector<std::string> &keys,
		     const SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    view_lod_append_unique_key(keys, view_lod_string_key("instance",
			       source->instanceKey.getValue()));
    view_lod_append_unique_key(keys, view_lod_string_key("realization",
			       source->realizationIdentity.getValue()));
    view_lod_append_unique_key(keys, view_lod_string_key("path",
			       source->path.getValue()));
    if (source->path.getValue().getLength() > 1 &&
	source->path.getValue().getString()[0] == '/')
	view_lod_append_unique_key(keys, std::string("path:") +
				   (source->path.getValue().getString() + 1));
    view_lod_append_unique_key(keys, view_lod_string_key("name",
			       SbString(view_lod_leaf_name_from_path(
					    source->path.getValue()))));
    view_lod_append_unique_key(keys, view_lod_source_primary_key(source));
}

BObolViewLodState::MeshPayload::MeshPayload(void) :
    sourceContentHash(0),
    resultKind(BOBOL_LOD_RESULT_NONE),
    qualityTier(BOBOL_LOD_QUALITY_METADATA),
    providerStatus(BOBOL_LOD_PROVIDER_UNKNOWN),
    activeCut(-1),
    residentCut(-1),
    requestedCut(-1),
    residentAdmissionRevision(0),
    viewRevision(0),
    policyRevision(0),
    hasSnappedPoints(FALSE),
    hasNormals(FALSE),
    shadedCullBackfaces(FALSE),
    memoryLimited(FALSE)
{
}

SbBool
BObolViewLodState::MeshPayload::isValid(void) const
{
    return this->mesh.isValid() ||
	(this->progressiveMesh && this->progressiveMesh->isValid());
}

size_t
BObolViewLodState::MeshPayload::estimateBytes(void) const
{
    return sizeof(*this) + (this->progressiveMesh ? 0 :
	    this->mesh.points.size() * sizeof(SbVec3f) +
	   this->mesh.normals.size() * sizeof(SbVec3f) +
	   this->mesh.coordIndex.size() * sizeof(int32_t) +
	   this->mesh.faceIndex.size() * sizeof(int32_t) +
	   this->mesh.vertexIndex.size() * sizeof(int32_t));
}

int
BObolViewLodState::MeshPayload::getTriangleCount(void) const
{
    return static_cast<int>(this->mesh.coordIndex.size() / 3);
}

SbBool
BObolViewLodState::MeshPayload::getTriangleVertexIndices(int triangleIndex,
	int &indexA,
	int &indexB,
	int &indexC) const
{
    indexA = -1;
    indexB = -1;
    indexC = -1;
    if (triangleIndex < 0)
	return FALSE;

    size_t base = static_cast<size_t>(triangleIndex) * 3;
    if (base + 2 >= this->mesh.coordIndex.size())
	return FALSE;

    indexA = this->mesh.coordIndex[base];
    indexB = this->mesh.coordIndex[base + 1];
    indexC = this->mesh.coordIndex[base + 2];
    if (indexA < 0 || indexB < 0 || indexC < 0 ||
	static_cast<size_t>(indexA) >= this->mesh.points.size() ||
	static_cast<size_t>(indexB) >= this->mesh.points.size() ||
	static_cast<size_t>(indexC) >= this->mesh.points.size()) {
	indexA = -1;
	indexB = -1;
	indexC = -1;
	return FALSE;
    }

    return TRUE;
}

SbBool
BObolViewLodState::MeshPayload::getTriangle(int triangleIndex,
	SbVec3f &a,
	SbVec3f &b,
	SbVec3f &c) const
{
    int ia = -1;
    int ib = -1;
    int ic = -1;
    if (!this->getTriangleVertexIndices(triangleIndex, ia, ib, ic))
	return FALSE;

    a = this->mesh.points[static_cast<size_t>(ia)];
    b = this->mesh.points[static_cast<size_t>(ib)];
    c = this->mesh.points[static_cast<size_t>(ic)];
    return TRUE;
}

BObolViewLodState::ProxyPayload::ProxyPayload(void) :
    resultKind(BOBOL_LOD_RESULT_NONE),
    qualityTier(BOBOL_LOD_QUALITY_METADATA),
    providerStatus(BOBOL_LOD_PROVIDER_UNKNOWN),
    viewRevision(0),
    policyRevision(0),
    diagnostic("")
{
}

SbBool
BObolViewLodState::ProxyPayload::isValid(void) const
{
    return this->proxy.isValid();
}

size_t
BObolViewLodState::ProxyPayload::estimateBytes(void) const
{
    return sizeof(*this);
}

static BObolLodCounts
view_lod_progressive_counts(
    const BObolLodProgressiveMeshPtr &mesh,
    const std::vector<uint32_t> &chunkIds, int cut, SbBool hasNormals);

BObolViewLodState::CadPayload::CadPayload(void) :
    preparedCadGeometryRevision(0),
    sourceRoutingId(0),
    sourceEntryIndex(UINT32_MAX),
    sourceContentHash(0),
    databaseRevision(0),
    sourceRevision(0),
    resultKind(BOBOL_LOD_RESULT_NONE),
    qualityTier(BOBOL_LOD_QUALITY_METADATA),
    providerStatus(BOBOL_LOD_PROVIDER_UNKNOWN),
    drawMode(BOBOL_LOD_DRAW_UNKNOWN),
    activeCut(-1),
    residentCut(-1),
    requestedCut(-1),
    allocatedCut(-1),
    allocationDrawMode(BOBOL_LOD_DRAW_UNKNOWN),
    allocationViewRevision(0),
    allocationPolicyRevision(0),
    allocationMeshRevision(0),
    allocationPlanSerial(0),
    projectedPixelDiameter(0.0f),
    projectedPixelArea(0.0f),
    projectedPixelPerimeter(0.0f),
    projectedBoundsContained(FALSE),
    targetPixelError(0.0f),
    projectedCutCountsViewRevision(0),
    projectedCutCountsPolicyRevision(0),
    projectedCutCountsMeshRevision(0),
    residentAdmissionRevision(0),
    viewRevision(0),
    policyRevision(0),
    visualEmphasis(0),
    hasSnappedPoints(FALSE),
    hasNormals(FALSE),
    shadedCullBackfaces(FALSE),
    memoryLimited(FALSE),
    ownerState(NULL)
{
}

/* Record the exact partial-frustum population while a bounded demand window
 * already owns the projection work.  The quiet allocator consumes this
 * immutable census later instead of rewalking every clustered hierarchy in
 * one unbounded owner-thread operation.  Null is an intentional result for a
 * fully contained or unclustered occurrence and therefore receives the same
 * epoch stamp; repeated retargets in one epoch do not retry the query.
 *
 * Allocation failure is a conservative degradation, not a failed draw.  The
 * caller retains whole-prefix accounting and the allocator may use its
 * existing exact fallback if memory becomes available. */
static bool
view_lod_update_projected_cut_counts(
    BObolViewLodState::CadPayload *payload,
    const BObolLodRequest &demand)
{
    if (!payload)
	return false;
    const uint64_t viewRevision = demand.viewRevision.value();
    const uint64_t policyRevision = demand.policyRevision.value();
    const uint64_t meshRevision = payload->progressiveMesh ?
	payload->progressiveMesh->revision() : 0;
    if (payload->projectedCutCountsViewRevision == viewRevision &&
	payload->projectedCutCountsPolicyRevision == policyRevision &&
	payload->projectedCutCountsMeshRevision == meshRevision)
	return false;

    const bool hadPartialCounts =
	static_cast<bool>(payload->projectedCutCounts);
    std::shared_ptr<BObolViewLodState::CadPayload::ProjectedCutCounts>
	projectedCounts;
    if (demand.spatialProjectionValid && payload->progressiveMesh &&
	payload->progressiveMesh->isValid() &&
	payload->progressiveMesh->hasSpatialClusters()) {
	try {
	    projectedCounts = std::make_shared<
		BObolViewLodState::CadPayload::ProjectedCutCounts>();
	} catch (const std::bad_alloc &) {
	    projectedCounts.reset();
	}
	if (projectedCounts &&
	    !payload->progressiveMesh->visibleCountsAtCuts(
		demand.localToRoot, demand.viewProjection,
		payload->hasNormals, projectedCounts->data(),
		projectedCounts->size()))
	    projectedCounts.reset();
    }
    payload->projectedCutCounts = projectedCounts;
    payload->projectedCutCountsViewRevision = viewRevision;
    payload->projectedCutCountsPolicyRevision = policyRevision;
    payload->projectedCutCountsMeshRevision = meshRevision;
    return hadPartialCounts || static_cast<bool>(projectedCounts);
}

static BObolLodCounts
view_lod_cad_counts_at_cut(
    const BObolViewLodState::CadPayload *payload, int cut)
{
    BObolLodCounts counts;
    if (!payload || cut < 0)
	return counts;
    if (payload->projectedCutCounts &&
	static_cast<size_t>(cut) < payload->projectedCutCounts->size())
	return (*payload->projectedCutCounts)[static_cast<size_t>(cut)];
    return view_lod_progressive_counts(payload->progressiveMesh,
	payload->requiredChunks, cut, payload->hasNormals);
}

SbBool
BObolViewLodState::CadPayload::isValid(void) const
{
    if (this->providerStatus != BOBOL_LOD_PROVIDER_READY)
	return FALSE;
    if (this->resultKind == BOBOL_LOD_RESULT_MESH ||
	this->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
	return this->mesh.isValid() ||
	    (this->progressiveMesh && this->progressiveMesh->isValid());
    if (this->resultKind == BOBOL_LOD_RESULT_AABB ||
	this->resultKind == BOBOL_LOD_RESULT_PROXY)
	return this->proxy.isValid();
    return FALSE;
}

size_t
BObolViewLodState::CadPayload::estimateBytes(void) const
{
    return (this->progressiveMesh ? 0 :
	    this->mesh.points.size() * sizeof(SbVec3f) +
	   this->mesh.normals.size() * sizeof(SbVec3f) +
	   this->mesh.coordIndex.size() * sizeof(int32_t) +
	   this->mesh.faceIndex.size() * sizeof(int32_t) +
	   this->mesh.vertexIndex.size() * sizeof(int32_t)) +
	   this->requiredChunks.capacity() * sizeof(uint32_t) +
	   this->presentedChunks.capacity() * sizeof(uint32_t) +
	   (this->projectedCutCounts ? sizeof(ProjectedCutCounts) : 0) +
	   (this->renderCostMetrics ? sizeof(RenderCostMetrics) : 0) +
	   sizeof(*this);
}

/* One compact occurrence has one active display representation.  Retaining
 * every completed stage both wastes memory and makes an unordered binding map
 * choose AABB versus OBB versus mesh by accident. */
static int
view_lod_cad_payload_rank(const BObolViewLodState::CadPayload &payload)
{
    int rank = payload.qualityTier * 16;
    switch (payload.resultKind) {
	case BOBOL_LOD_RESULT_FULL_DETAIL:
	    return 512 + rank;
	case BOBOL_LOD_RESULT_MESH:
	    return 448 + rank;
	case BOBOL_LOD_RESULT_PROXY:
	    if (payload.proxy.kind == BOBOL_LOD_PROXY_OBB)
		return 320 + rank;
	    return 256 + rank;
	case BOBOL_LOD_RESULT_AABB:
	    return 192 + rank;
	default:
	    return rank;
    }
}

static size_t
view_lod_saturating_add(size_t lhs, size_t rhs)
{
    return rhs > SIZE_MAX - lhs ? SIZE_MAX : lhs + rhs;
}

static uint64_t
view_lod_saturating_add_u64(uint64_t lhs, uint64_t rhs)
{
    return rhs > UINT64_MAX - lhs ? UINT64_MAX : lhs + rhs;
}

static size_t
view_lod_saturating_subtract(size_t lhs, size_t rhs)
{
    return rhs > lhs ? 0 : lhs - rhs;
}

static unsigned int
view_lod_cad_payload_channel_mask(
    const BObolViewLodState::CadPayload *payload)
{
    if (!payload)
	return 0;
    unsigned int mask = 0;
    if (payload->drawMode == BOBOL_LOD_DRAW_WIRE ||
	payload->drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	mask |= 1u;
    if (payload->drawMode == BOBOL_LOD_DRAW_SHADED ||
	payload->drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	payload->drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	mask |= 2u;
    return mask;
}

static SbBool
view_lod_cad_payload_targets(const BObolViewLodState::CadPayload &payload,
	const std::string &sourceBindingKey, const SbString &occurrenceKey)
{
    return bu_strcmp(payload.sourceBindingKey.getString(),
		     sourceBindingKey.c_str()) == 0 &&
	bu_strcmp(payload.sourceInstanceKey.getString(),
		  occurrenceKey.getString()) == 0 ? TRUE : FALSE;
}

static std::string
view_lod_cad_result_binding_key(const std::string &sourceBindingKey,
	const std::string &resultKey)
{
    std::string key = "source:";
    key += sourceBindingKey;
    key += '|';
    key += resultKey;
    return key;
}

static std::string
view_lod_cad_occurrence_key(const SbString &occurrenceKey)
{
    return occurrenceKey.getLength() > 0 ?
	std::string(occurrenceKey.getString()) : std::string();
}

static bool
view_lod_provider_status_is_terminal_failure(int status)
{
    return status == BOBOL_LOD_PROVIDER_CACHE_MISS ||
	status == BOBOL_LOD_PROVIDER_STALE ||
	status == BOBOL_LOD_PROVIDER_ERROR;
}

static void
view_lod_remove_superseded_cad_payloads(
	std::unordered_map<std::string, BObolViewLodState::CadPayloadPtr>
	&bindings, const std::string &sourceBindingKey,
	const SbString &occurrenceKey,
	const BObolViewLodState::CadPayload &candidate)
{
    const int candidateRank = view_lod_cad_payload_rank(candidate);
    for (auto binding = bindings.begin(); binding != bindings.end();) {
	const BObolViewLodState::CadPayloadPtr &payload = binding->second;
	if (payload && view_lod_cad_payload_targets(*payload,
		sourceBindingKey, occurrenceKey) &&
	    view_lod_cad_payload_rank(*payload) <= candidateRank) {
	    binding = bindings.erase(binding);
	    continue;
	}
	++binding;
    }
}

BObolViewLodState::BObolViewLodState(void) :
    cadFullResyncRevision(1),
    cadBindingsRevision(1),
    normalStyle(NORMAL_AUTHORED),
    normalCreaseAngle(60.0f),
    cadPresentationProgressiveLodCeiling(-1),
    cadPresentationPointProxyPixelThreshold(1.0f),
    cadPresentationDiscoveryPointProxyPixelThreshold(1.0f),
    cadPresentationCameraMotionFrameReuse(FALSE),
    cadPresentationFrameObservationArmed(FALSE),
    cadPresentationFrameStatusValid(FALSE),
    cadLastPresentedPrimitiveCountValid(FALSE),
    cadLastPresentedPrimitiveCount(0),
    cadLastPresentedRenderCostValid(FALSE),
    cadLastPresentedRenderCost(0),
    cadLastPresentationFrameExact(FALSE),
    cadLastSubpixelProxyCount(0),
    cadLastUncollapsedStructuralProxyCount(0),
    cadLastGpuMeasurementValid(FALSE),
    cadLastGpuFaces(0),
    cadLastGpuNanoseconds(0),
    cadLastGpuSerial(0),
    cadLastGpuPointProxyPixelThreshold(1.0f),
    cadLastPreparedReplay(FALSE),
    cadGpuResourceStatusValid(FALSE),
    cadValidPayloadCount(0),
    cadMeshPayloadCountValue(0),
    cadProxyPayloadCountValue(0),
    cadSatisfiedMeshPayloadCount(0),
    cadMemoryLimitedMeshPayloadCount(0),
    cadActiveFaceCount(0),
    cadActiveRenderCost(0),
    cadMinimumActiveRenderCost(0),
    cadDisplayMeshBytes(0),
    cadResidentDemandRevision(1),
    cadActiveAllocationPlanSerial(0),
    cadNextAllocationPlanSerial(0)
{
    memset(this->cadProxyKindCounts, 0,
	sizeof(this->cadProxyKindCounts));
    memset(this->cadProgressiveCutCounts, 0,
	sizeof(this->cadProgressiveCutCounts));
    memset(this->cadProgressiveCeilingRenderCosts, 0,
	sizeof(this->cadProgressiveCeilingRenderCosts));
}

BObolViewLodState::~BObolViewLodState(void)
{
    this->clearCadPresentations();
}

static bool
view_lod_cad_payload_is_mesh(
    const BObolViewLodState::CadPayload *payload)
{
    return payload &&
	(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL);
}

static bool
view_lod_cad_payload_is_satisfied(
    const BObolViewLodState::CadPayload *payload)
{
    return view_lod_cad_payload_is_mesh(payload) &&
	payload->presentedChunks == payload->requiredChunks &&
	(payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL ||
	 payload->requestedCut < 0 ||
	payload->activeCut >= payload->requestedCut);
}

static SbBool
view_lod_progressive_can_draw(
    const BObolLodProgressiveMeshPtr &mesh,
    const std::vector<uint32_t> &chunkIds, int cut)
{
    if (!mesh)
	return FALSE;
    return chunkIds.empty() ? mesh->canDrawCut(cut) :
	mesh->canDrawChunksAtCut(chunkIds, cut);
}

static BObolLodCounts
view_lod_progressive_counts(
    const BObolLodProgressiveMeshPtr &mesh,
    const std::vector<uint32_t> &chunkIds, int cut, SbBool hasNormals)
{
    BObolLodCounts counts;
    if (mesh && !chunkIds.empty() &&
	mesh->countsForChunksAtCut(chunkIds, cut, hasNormals, &counts))
	return counts;
    return bobol_lod_progressive_counts(mesh, cut, hasNormals);
}

static void
view_lod_cad_payload_render_cost_metrics(
    const BObolViewLodState::CadPayload *payload,
    BObolViewLodState::CadPayload::RenderCostMetrics &metrics)
{
    metrics.fill(0);
    if (!payload || !payload->isValid() ||
	!view_lod_cad_payload_is_mesh(payload))
	return;

    BObolLodCounts minimumCounts = payload->counts;
    if (!payload->progressiveMesh || payload->activeCut < 0) {
	const size_t cost = bobol_lod_render_cost_units(
	    payload->counts, payload->drawMode, 1);
	metrics[0] = cost;
	for (int cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut)
	    metrics[static_cast<size_t>(cut) + 1] = cost;
	return;
    }

    std::array<BObolLodCounts, BOBOL_MESH_LOD_CUT_COUNT_MAX> populations;
    int minimumCut = -1;
    int maximumDrawableCut = -1;
    if (!payload->progressiveMesh->drawableCountsAtCuts(
	    payload->requiredChunks, payload->hasNormals,
	    populations.data(), populations.size(), &minimumCut,
	    &maximumDrawableCut) || minimumCut < 0 ||
	maximumDrawableCut < minimumCut) {
	const size_t cost = bobol_lod_render_cost_units(
	    payload->counts, payload->drawMode, 1);
	metrics[0] = cost;
	for (int cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut)
	    metrics[static_cast<size_t>(cut) + 1] = cost;
	return;
    }

    minimumCounts = populations[static_cast<size_t>(minimumCut)];
    metrics[0] = bobol_lod_render_cost_units(
	minimumCounts, payload->drawMode, 1);
    for (int ceiling = 0; ceiling < BOBOL_MESH_LOD_CUT_COUNT_MAX;
	    ++ceiling) {
	int effectiveCut = std::min(payload->activeCut,
	    std::max(0, ceiling));
	effectiveCut = std::max(minimumCut, effectiveCut);
	const BObolLodCounts &counts = effectiveCut <= maximumDrawableCut ?
	    populations[static_cast<size_t>(effectiveCut)] : payload->counts;
	metrics[static_cast<size_t>(ceiling) + 1] =
	    bobol_lod_render_cost_units(counts, payload->drawMode, 1);
    }
}

/*
 * Residency and presentation answer different questions.  activeCut is
 * the cut selected by the aggregate frame budget; requestedCut is the cut
 * justified by the current view.  Preserve the latter only as far as this
 * immutable generation is already resident.  This makes performance-only
 * coarsening a draw-count change, while zoom-out and memory-pressure recovery
 * can still reclaim suffixes without promising data which has not loaded.
 */
static int
view_lod_resident_demand_cut(int activeCut, int requestedCut,
	int residentCut)
{
    if (activeCut < 0)
	return -1;
    if (requestedCut < 0 || residentCut < 0)
	return activeCut;
    return std::max(activeCut, std::min(requestedCut, residentCut));
}

static int
view_lod_cad_resident_demand_cut(
    const BObolViewLodState::CadPayload *payload)
{
    return payload ? view_lod_resident_demand_cut(payload->activeCut,
	payload->requestedCut, payload->residentCut) : -1;
}

void
BObolViewLodState::clearCadPayloadMetrics(void)
{
    this->cadResidentDemandRevision++;
    if (!this->cadResidentDemandRevision)
	this->cadResidentDemandRevision++;
    this->cadValidPayloadCount = 0;
    this->cadMeshPayloadCountValue = 0;
    this->cadProxyPayloadCountValue = 0;
    this->cadSatisfiedMeshPayloadCount = 0;
    this->cadMemoryLimitedMeshPayloadCount = 0;
    this->cadActiveFaceCount = 0;
    this->cadActiveRenderCost = 0;
    this->cadMinimumActiveRenderCost = 0;
    this->cadDisplayMeshBytes = 0;
    memset(this->cadProxyKindCounts, 0,
	sizeof(this->cadProxyKindCounts));
    memset(this->cadProgressiveCutCounts, 0,
	sizeof(this->cadProgressiveCutCounts));
    memset(this->cadProgressiveCeilingRenderCosts, 0,
	sizeof(this->cadProgressiveCeilingRenderCosts));
    this->cadMeshPayloadCountsBySource.clear();
    this->cadUnsatisfiedOccurrencesBySource.clear();
    this->cadPayloadsByAssetKey.clear();
    this->cadResidentDemandStates.clear();
    this->cadResidentDemands.clear();
}

void
BObolViewLodState::addCadResidentDemand(const CadPayload *payload)
{
    const int demandCut = view_lod_cad_resident_demand_cut(payload);
    if (!payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() ||
	payload->cacheKey.getLength() == 0 ||
	demandCut < 0 || demandCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return;
    this->cadResidentDemandRevision++;
    if (!this->cadResidentDemandRevision)
	this->cadResidentDemandRevision++;

    const std::string key = payload->cacheKey.getString();
    const unsigned int channelMask =
	view_lod_cad_payload_channel_mask(payload) & 3u;
    auto inserted = this->cadResidentDemandStates.emplace(
	key, CadResidentDemandState());
    CadResidentDemandState &state = inserted.first->second;
    if (inserted.second) {
	state.demandIndex = this->cadResidentDemands.size();
	BObolLodResidentDemand demand;
	demand.assetKey = payload->cacheKey;
	demand.cut = demandCut;
	demand.channelMask = channelMask;
	demand.chunkIds = payload->requiredChunks;
	this->cadResidentDemands.push_back(demand);
    }
    state.cutCounts[demandCut]++;
    if (demandCut > state.maximumCut) {
	state.maximumCut = demandCut;
	this->cadResidentDemands[state.demandIndex].cut =
	    state.maximumCut;
    }
    state.channelCounts[channelMask]++;
    state.channelMask |= channelMask;
    this->cadResidentDemands[state.demandIndex].channelMask =
	state.channelMask;
    if (!inserted.second) {
	std::vector<uint32_t> &chunks =
	    this->cadResidentDemands[state.demandIndex].chunkIds;
	for (uint32_t chunk : payload->requiredChunks) {
	    size_t &refs = state.chunkCounts[chunk];
	    if (refs++ == 0)
		chunks.insert(std::lower_bound(
		    chunks.begin(), chunks.end(), chunk), chunk);
	}
    } else {
	for (uint32_t chunk : payload->requiredChunks)
	    state.chunkCounts[chunk] = 1;
    }
}

void
BObolViewLodState::removeCadResidentDemand(const CadPayload *payload)
{
    const int demandCut = view_lod_cad_resident_demand_cut(payload);
    if (!payload || payload->cacheKey.getLength() == 0 ||
	demandCut < 0 || demandCut >= BOBOL_MESH_LOD_CUT_COUNT_MAX)
	return;
    this->cadResidentDemandRevision++;
    if (!this->cadResidentDemandRevision)
	this->cadResidentDemandRevision++;
    const auto found = this->cadResidentDemandStates.find(
	payload->cacheKey.getString());
    if (found == this->cadResidentDemandStates.end())
	return;
    CadResidentDemandState &state = found->second;
    const unsigned int channelMask =
	view_lod_cad_payload_channel_mask(payload) & 3u;
    if (state.cutCounts[demandCut])
	state.cutCounts[demandCut]--;
    if (state.channelCounts[channelMask])
	state.channelCounts[channelMask]--;
    std::vector<uint32_t> &chunks =
	this->cadResidentDemands[state.demandIndex].chunkIds;
    for (uint32_t chunk : payload->requiredChunks) {
	const auto foundChunk = state.chunkCounts.find(chunk);
	if (foundChunk == state.chunkCounts.end())
	    continue;
	if (foundChunk->second > 1) {
	    --foundChunk->second;
	    continue;
	}
	state.chunkCounts.erase(foundChunk);
	const auto visible = std::lower_bound(
	    chunks.begin(), chunks.end(), chunk);
	if (visible != chunks.end() && *visible == chunk)
	    chunks.erase(visible);
    }
    state.channelMask = 0;
    for (unsigned int mask = 0; mask < 4; ++mask) {
	if (state.channelCounts[mask])
	    state.channelMask |= mask;
    }
    if (demandCut == state.maximumCut &&
	!state.cutCounts[demandCut]) {
	state.maximumCut = -1;
	for (int cut = demandCut - 1; cut >= 0; --cut) {
	    if (!state.cutCounts[cut])
		continue;
	    state.maximumCut = cut;
	    break;
	}
    }
    if (state.maximumCut >= 0) {
	this->cadResidentDemands[state.demandIndex].cut =
	    state.maximumCut;
	this->cadResidentDemands[state.demandIndex].channelMask =
	    state.channelMask;
	return;
    }

    const size_t removedIndex = state.demandIndex;
    const size_t lastIndex = this->cadResidentDemands.size() - 1;
    if (removedIndex != lastIndex) {
	this->cadResidentDemands[removedIndex] =
	    std::move(this->cadResidentDemands[lastIndex]);
	const auto replacement = this->cadResidentDemandStates.find(
	    this->cadResidentDemands[removedIndex].assetKey.getString());
	if (replacement != this->cadResidentDemandStates.end())
	    replacement->second.demandIndex = removedIndex;
    }
    this->cadResidentDemands.pop_back();
    this->cadResidentDemandStates.erase(found);
}

void
BObolViewLodState::addCadPayloadMetrics(const CadPayload *payload)
{
    if (!payload || !payload->isValid())
	return;

    CadPayload *mutablePayload = const_cast<CadPayload *>(payload);
    CadPayload::RenderCostMetrics renderCostMetrics;
    if (view_lod_cad_payload_is_mesh(payload)) {
	view_lod_cad_payload_render_cost_metrics(
	    payload, renderCostMetrics);
	if (!mutablePayload->renderCostMetrics)
	    mutablePayload->renderCostMetrics =
		std::make_unique<CadPayload::RenderCostMetrics>();
	*mutablePayload->renderCostMetrics = renderCostMetrics;
    } else
	mutablePayload->renderCostMetrics.reset();

    this->addCadResidentDemand(payload);
    if (payload->progressiveMesh &&
	payload->cacheKey.getLength() > 0)
	this->cadPayloadsByAssetKey[payload->cacheKey.getString()].insert(
	    const_cast<CadPayload *>(payload));
    this->cadValidPayloadCount =
	view_lod_saturating_add(this->cadValidPayloadCount, 1);
    this->cadDisplayMeshBytes = view_lod_saturating_add(
	this->cadDisplayMeshBytes, payload->estimateBytes());

    if (payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) {
	this->cadMeshPayloadCountValue =
	    view_lod_saturating_add(this->cadMeshPayloadCountValue, 1);
	if (!view_lod_cad_payload_is_satisfied(payload) &&
	    payload->sourceInstanceKey.getLength() > 0)
	    this->cadUnsatisfiedOccurrencesBySource[
		payload->sourceBindingKey.getString()].insert(
		    payload->sourceInstanceKey.getString());
	size_t &sourceMeshCount = this->cadMeshPayloadCountsBySource[
	    payload->sourceBindingKey.getString()];
	sourceMeshCount =
	    view_lod_saturating_add(sourceMeshCount, 1);
	this->cadActiveFaceCount = view_lod_saturating_add(
	    this->cadActiveFaceCount, payload->counts.faceCount);
	this->cadActiveRenderCost = view_lod_saturating_add(
	    this->cadActiveRenderCost,
	    bobol_lod_render_cost_units(
		payload->counts, payload->drawMode, 1));
	this->cadMinimumActiveRenderCost = view_lod_saturating_add(
	    this->cadMinimumActiveRenderCost,
	    renderCostMetrics[0]);
	if (view_lod_cad_payload_is_satisfied(payload))
	    this->cadSatisfiedMeshPayloadCount = view_lod_saturating_add(
		this->cadSatisfiedMeshPayloadCount, 1);
	if (payload->memoryLimited)
	    this->cadMemoryLimitedMeshPayloadCount =
		view_lod_saturating_add(
		    this->cadMemoryLimitedMeshPayloadCount, 1);
	if (payload->progressiveMesh && payload->activeCut >= 0 &&
	    payload->activeCut <
		static_cast<int>(sizeof(this->cadProgressiveCutCounts) /
		    sizeof(this->cadProgressiveCutCounts[0])))
	    this->cadProgressiveCutCounts[payload->activeCut] =
		view_lod_saturating_add(
		    this->cadProgressiveCutCounts[payload->activeCut], 1);
	for (int cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut)
	    this->cadProgressiveCeilingRenderCosts[cut] =
		view_lod_saturating_add(
		    this->cadProgressiveCeilingRenderCosts[cut],
		    renderCostMetrics[static_cast<size_t>(cut) + 1]);
	return;
    }

    if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	payload->resultKind == BOBOL_LOD_RESULT_PROXY) {
	this->cadProxyPayloadCountValue =
	    view_lod_saturating_add(this->cadProxyPayloadCountValue, 1);
	const int kind = payload->proxy.kind;
	if (kind >= 0 &&
	    kind < static_cast<int>(sizeof(this->cadProxyKindCounts) /
		sizeof(this->cadProxyKindCounts[0])))
	    this->cadProxyKindCounts[kind] =
		view_lod_saturating_add(this->cadProxyKindCounts[kind], 1);
    }
}

void
BObolViewLodState::removeCadPayloadMetrics(const CadPayload *payload)
{
    if (!payload || !payload->isValid())
	return;

    CadPayload::RenderCostMetrics fallbackRenderCostMetrics;
    const CadPayload::RenderCostMetrics *renderCostMetrics =
	payload->renderCostMetrics.get();
    if (view_lod_cad_payload_is_mesh(payload) && !renderCostMetrics) {
	/* This is a defensive path for payloads installed by older/internal
	 * callers.  Ordinary insertion always leaves the exact immutable table
	 * above, so removal performs no shared-generation queries. */
	view_lod_cad_payload_render_cost_metrics(
	    payload, fallbackRenderCostMetrics);
	renderCostMetrics = &fallbackRenderCostMetrics;
    }

    this->removeCadResidentDemand(payload);
    if (payload->cacheKey.getLength() > 0) {
	const auto indexed = this->cadPayloadsByAssetKey.find(
	    payload->cacheKey.getString());
	if (indexed != this->cadPayloadsByAssetKey.end()) {
	    indexed->second.erase(const_cast<CadPayload *>(payload));
	    if (indexed->second.empty())
		this->cadPayloadsByAssetKey.erase(indexed);
	}
    }
    this->cadValidPayloadCount =
	view_lod_saturating_subtract(this->cadValidPayloadCount, 1);
    this->cadDisplayMeshBytes = view_lod_saturating_subtract(
	this->cadDisplayMeshBytes, payload->estimateBytes());

    if (payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) {
	this->cadMeshPayloadCountValue =
	    view_lod_saturating_subtract(this->cadMeshPayloadCountValue, 1);
	if (payload->sourceInstanceKey.getLength() > 0) {
	    const auto unsatisfied =
		this->cadUnsatisfiedOccurrencesBySource.find(
		    payload->sourceBindingKey.getString());
	    if (unsatisfied !=
		    this->cadUnsatisfiedOccurrencesBySource.end()) {
		unsatisfied->second.erase(
		    payload->sourceInstanceKey.getString());
		if (unsatisfied->second.empty())
		    this->cadUnsatisfiedOccurrencesBySource.erase(
			unsatisfied);
	    }
	}
	const auto sourceMeshCount = this->cadMeshPayloadCountsBySource.find(
	    payload->sourceBindingKey.getString());
	if (sourceMeshCount != this->cadMeshPayloadCountsBySource.end()) {
	    sourceMeshCount->second =
		view_lod_saturating_subtract(sourceMeshCount->second, 1);
	    if (!sourceMeshCount->second)
		this->cadMeshPayloadCountsBySource.erase(sourceMeshCount);
	}
	this->cadActiveFaceCount = view_lod_saturating_subtract(
	    this->cadActiveFaceCount, payload->counts.faceCount);
	this->cadActiveRenderCost = view_lod_saturating_subtract(
	    this->cadActiveRenderCost,
	    bobol_lod_render_cost_units(
		payload->counts, payload->drawMode, 1));
	this->cadMinimumActiveRenderCost = view_lod_saturating_subtract(
	    this->cadMinimumActiveRenderCost,
	    renderCostMetrics ? (*renderCostMetrics)[0] : 0);
	if (view_lod_cad_payload_is_satisfied(payload))
	    this->cadSatisfiedMeshPayloadCount = view_lod_saturating_subtract(
		this->cadSatisfiedMeshPayloadCount, 1);
	if (payload->memoryLimited)
	    this->cadMemoryLimitedMeshPayloadCount =
		view_lod_saturating_subtract(
		    this->cadMemoryLimitedMeshPayloadCount, 1);
	if (payload->progressiveMesh && payload->activeCut >= 0 &&
	    payload->activeCut <
		static_cast<int>(sizeof(this->cadProgressiveCutCounts) /
		    sizeof(this->cadProgressiveCutCounts[0])))
	    this->cadProgressiveCutCounts[payload->activeCut] =
		view_lod_saturating_subtract(
		    this->cadProgressiveCutCounts[payload->activeCut], 1);
	for (int cut = 0; cut < BOBOL_MESH_LOD_CUT_COUNT_MAX; ++cut)
	    this->cadProgressiveCeilingRenderCosts[cut] =
		view_lod_saturating_subtract(
		    this->cadProgressiveCeilingRenderCosts[cut],
		    renderCostMetrics ?
			(*renderCostMetrics)[static_cast<size_t>(cut) + 1] : 0);
	return;
    }

    if (payload->resultKind == BOBOL_LOD_RESULT_AABB ||
	payload->resultKind == BOBOL_LOD_RESULT_PROXY) {
	this->cadProxyPayloadCountValue =
	    view_lod_saturating_subtract(
		this->cadProxyPayloadCountValue, 1);
	const int kind = payload->proxy.kind;
	if (kind >= 0 &&
	    kind < static_cast<int>(sizeof(this->cadProxyKindCounts) /
		sizeof(this->cadProxyKindCounts[0])))
	    this->cadProxyKindCounts[kind] = view_lod_saturating_subtract(
		this->cadProxyKindCounts[kind], 1);
    }
}

void
BObolViewLodState::clear(void)
{
    this->clearCadPresentations();
    this->meshBindings.clear();
    this->proxyBindings.clear();
    this->cadSourceBindings.clear();
    this->cadSourceEntryBindings.clear();
    this->cadAssetBindings.clear();
    this->cadBindings.clear();
    this->cadOccurrenceFailures.clear();
    this->cadOccurrenceChanges.clear();
    this->clearCadPayloadMetrics();
    this->noteResidentMeshesChanged("view-state-clear");
}

SoCADAssembly *
BObolViewLodState::findCadPresentation(
    const SoBRLDatabaseSource *source,
    SbString *contentKey) const
{
    if (!source)
	return NULL;
    const std::string key = view_lod_source_primary_key(source);
    const auto found = this->cadPresentations.find(key);
    if (found == this->cadPresentations.end())
	return NULL;
    if (found->second.sourceRoutingId !=
	source->getCompactSourceRoutingId())
	return NULL;
    if (contentKey)
	*contentKey = found->second.contentKey;
    return found->second.assembly;
}

void
BObolViewLodState::setCadPresentation(
    const SoBRLDatabaseSource *source,
    SoCADAssembly *assembly,
    const SbString &contentKey) const
{
    if (!source)
	return;
    const std::string key = view_lod_source_primary_key(source);
    CadPresentation &presentation = this->cadPresentations[key];
    const uint64_t sourceRoutingId = source->getCompactSourceRoutingId();
    if (presentation.sourceRoutingId != sourceRoutingId &&
	presentation.assembly == assembly)
	assembly = NULL;
    if (presentation.assembly != assembly) {
	/* A frame-start observation contains non-owning assembly pointers.  A
	 * presentation replacement may release the last owning reference below,
	 * so retire the observation before changing ownership.  The next render
	 * traversal will arm a new coherent snapshot. */
	this->cadPresentationFrameStartExecutionSerials.clear();
	this->cadPresentationFrameObservationArmed = FALSE;
	if (assembly) {
	    assembly->ref();
	    this->cadPresentationAssemblyUseCounts[assembly]++;
	}
	if (presentation.assembly) {
	    const auto use = this->cadPresentationAssemblyUseCounts.find(
		presentation.assembly);
	    if (use != this->cadPresentationAssemblyUseCounts.end()) {
		if (use->second > 1)
		    use->second--;
		else
		    this->cadPresentationAssemblyUseCounts.erase(use);
	    }
	    presentation.assembly->unref();
	}
	presentation.assembly = assembly;
	this->cadPresentationFrameStatusValid = FALSE;
    }
    presentation.sourceRoutingId = sourceRoutingId;
    presentation.contentKey = contentKey;
    if (presentation.assembly)
	presentation.assembly->progressiveCutCeiling =
	    this->cadPresentationProgressiveLodCeiling;
    if (presentation.assembly)
	presentation.assembly->pointProxyPixelThreshold =
	    std::max(this->cadPresentationPointProxyPixelThreshold,
		this->cadPresentationDiscoveryPointProxyPixelThreshold);
    if (presentation.assembly)
	presentation.assembly->cameraMotionFrameReuse =
	    this->cadPresentationCameraMotionFrameReuse;
    if (!presentation.assembly)
	this->cadPresentations.erase(key);
}

void
BObolViewLodState::clearCadPresentations(void) const
{
    /* Invalidate the non-owning frame-start snapshot before unref'ing any
     * assemblies.  LoD policy changes can clear presentations between a GL
     * traversal and the following HUD/status query; retaining these pointers
     * made that query dereference deleted SoCADAssembly objects. */
    this->cadPresentationFrameStartExecutionSerials.clear();
    this->cadPresentationFrameObservationArmed = FALSE;
    for (auto &binding : this->cadPresentations) {
	if (binding.second.assembly)
	    binding.second.assembly->unref();
	binding.second.assembly = NULL;
    }
    this->cadPresentations.clear();
    this->cadPresentationAssemblyUseCounts.clear();
    this->cadPresentationFrameStatusValid = FALSE;
}

void
BObolViewLodState::setNormalStyle(NormalStyle style, float creaseAngleDegrees)
{
    if (style < NORMAL_AUTHORED || style > NORMAL_SMOOTH)
	style = NORMAL_AUTHORED;
    if (!std::isfinite(creaseAngleDegrees))
	creaseAngleDegrees = 60.0f;
    creaseAngleDegrees = std::max(0.0f,
	std::min(180.0f, creaseAngleDegrees));
    if (this->normalStyle == style &&
	std::fabs(this->normalCreaseAngle - creaseAngleDegrees) <= 1.0e-6f)
	return;
    this->normalStyle = style;
    this->normalCreaseAngle = creaseAngleDegrees;
    /* This is a presentation policy change.  Retain all view-local PoP
     * payloads and invalidate only their renderer assemblies. */
    this->clearCadPresentations();
}

BObolViewLodState::NormalStyle
BObolViewLodState::getNormalStyle(void) const
{
    return this->normalStyle;
}

float
BObolViewLodState::getNormalCreaseAngle(void) const
{
    return this->normalCreaseAngle;
}

SbBool
BObolViewLodState::applyMeshResult(const SoBRLMeshShape *shape,
				     const BObolLodResult &result)
{
    BObolLodResult copy = result;
    return this->applyMeshResultInternal(shape, copy, FALSE);
}

SbBool
BObolViewLodState::applyMeshResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume)
{
    if (!shape ||
	result.resultKind != BOBOL_LOD_RESULT_MESH ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	(!result.mesh.isValid() &&
	 (!result.progressiveMesh || !result.progressiveMesh->isValid())))
	return FALSE;

    MeshPayloadPtr payload(new MeshPayload);
    payload->progressiveMesh = result.progressiveMesh;
    if (consume)
	payload->mesh = std::move(result.mesh);
    else
	payload->mesh = result.mesh;
    payload->sourcePath = shape->sourcePath.getValue();
    payload->sourceName = shape->sourceName.getValue();
    payload->sourceIdentity = shape->sourceIdentity.getValue();
    payload->cacheIdentity = result.cacheKey.value;
    payload->cacheKey = result.geometry.cacheKey.isValid() ?
	result.geometry.cacheKey.value : result.cacheKey.value;
    payload->sourceContentHash = result.request.sourceContentHash;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->activeCut = result.geometry.activeCut;
    payload->residentCut = result.residentCut;
    payload->requestedCut = result.resolvedCut >= 0 ?
	result.resolvedCut : result.request.requestedCut;
    payload->requiredChunks = result.request.requiredChunks;
    payload->residentAdmissionRevision =
	result.residentAdmissionRevision;
    payload->viewRevision = result.request.viewRevision.value();
    payload->policyRevision = result.request.policyRevision.value();
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->hasSnappedPoints = result.hasSnappedPoints;
    payload->hasNormals = result.hasNormals;
    payload->shadedCullBackfaces = result.shadedCullBackfaces;
    payload->memoryLimited = result.memoryLimited;
    payload->diagnostic = result.diagnostic;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++)
	this->meshBindings[keys[i]] = payload;

    return TRUE;
}

SbBool
BObolViewLodState::applyProxyResult(const SoBRLMeshShape *shape,
				      const BObolLodResult &result)
{
    BObolLodResult copy = result;
    return this->applyProxyResultInternal(shape, copy, FALSE);
}

SbBool
BObolViewLodState::applyProxyResultInternal(const SoBRLMeshShape *shape,
	BObolLodResult &result, SbBool consume)
{
    if (!shape ||
	(result.resultKind != BOBOL_LOD_RESULT_AABB &&
	 result.resultKind != BOBOL_LOD_RESULT_PROXY) ||
	result.providerStatus != BOBOL_LOD_PROVIDER_READY ||
	!result.proxy.isValid())
	return FALSE;

    ProxyPayloadPtr payload(new ProxyPayload);
    if (consume)
	payload->proxy = std::move(result.proxy);
    else
	payload->proxy = result.proxy;
    payload->sourcePath = shape->sourcePath.getValue();
    payload->sourceName = shape->sourceName.getValue();
    payload->sourceIdentity = shape->sourceIdentity.getValue();
    payload->cacheIdentity = shape->cacheIdentity.getValue();
    payload->cacheKey = result.cacheKey.value;
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->viewRevision = result.request.viewRevision.value();
    payload->policyRevision = result.request.policyRevision.value();
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->diagnostic = result.diagnostic;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++)
	this->proxyBindings[keys[i]] = payload;

    return TRUE;
}

SbBool
BObolViewLodState::applyDisplayResult(const SoBRLMeshShape *shape,
					const BObolLodResult &result)
{
    if (result.resultKind == BOBOL_LOD_RESULT_MESH ||
	result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
	return this->applyMeshResult(shape, result);

    return this->applyProxyResult(shape, result);
}

SbBool
BObolViewLodState::consumeDisplayResult(const SoBRLMeshShape *shape,
	BObolLodResult &result)
{
    if (result.resultKind == BOBOL_LOD_RESULT_MESH ||
	result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)
	return this->applyMeshResultInternal(shape, result, TRUE);
    return this->applyProxyResultInternal(shape, result, TRUE);
}

SbBool
BObolViewLodState::applySourceResult(
    const SoBRLDatabaseSource *source,
    const BObolLodResult &result)
{
    BObolLodResult copy = result;
    return this->applySourceResultInternal(source, copy, FALSE);
}

SbBool
BObolViewLodState::consumeSourceResult(
    const SoBRLDatabaseSource *source,
    BObolLodResult &result)
{
    return this->applySourceResultInternal(source, result, TRUE);
}

SbBool
BObolViewLodState::applySourceResultInternal(
    const SoBRLDatabaseSource *source,
    BObolLodResult &result,
    SbBool consume)
{
    if (!source)
	return FALSE;

    const std::string sourceBindingKey = view_lod_source_primary_key(source);
    const std::string occurrenceKey =
	view_lod_cad_occurrence_key(result.request.occurrenceKey);
    if (view_lod_provider_status_is_terminal_failure(
	    result.providerStatus)) {
	if (getenv("BOBOL_LOD_TRACE_FAILURES"))
	    bu_log("BObol LoD terminal failure occurrence=%s object=%s "
		   "status=%d cut=%d diagnostic=%s\n",
		   result.request.occurrenceKey.getString(),
		   result.request.objectPath.getString(),
		   result.providerStatus, result.request.requestedCut,
		   result.diagnostic.getString());
	CadOccurrenceFailure &failure =
	    this->cadOccurrenceFailures[sourceBindingKey][occurrenceKey];
	failure.databaseRevision = result.request.databaseRevision.value();
	failure.sourceRevision = result.request.sourceRevision.value();
	failure.sourceContentHash = result.request.sourceContentHash;
	failure.viewRevision = result.request.viewRevision.value();
	failure.policyRevision = result.request.policyRevision.value();
	failure.requestedCut = result.request.requestedCut;
	failure.drawMode = result.request.drawMode;
	failure.qualityTier = result.request.qualityTier;
	failure.providerStatus = result.providerStatus;

	/* A retained lower cut remains useful, but this exact missing suffix is
	 * no longer refinement work.  Keep the payload and remove only its
	 * retry-frontier entry. */
	if (result.request.occurrenceKey.getLength() > 0) {
	    const auto unsatisfied =
		this->cadUnsatisfiedOccurrencesBySource.find(
		    sourceBindingKey);
	    if (unsatisfied !=
		    this->cadUnsatisfiedOccurrencesBySource.end()) {
		unsatisfied->second.erase(occurrenceKey);
		if (unsatisfied->second.empty())
		    this->cadUnsatisfiedOccurrencesBySource.erase(
			unsatisfied);
	    }
	}
	return TRUE;
    }
    if (result.providerStatus != BOBOL_LOD_PROVIDER_READY)
	return FALSE;

    if ((result.resultKind == BOBOL_LOD_RESULT_MESH ||
	 result.resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) &&
	!result.mesh.isValid() &&
	(!result.progressiveMesh || !result.progressiveMesh->isValid()))
	return FALSE;
    if ((result.resultKind == BOBOL_LOD_RESULT_AABB ||
	 result.resultKind == BOBOL_LOD_RESULT_PROXY) &&
	!result.proxy.isValid())
	return FALSE;
    if (result.resultKind != BOBOL_LOD_RESULT_MESH &&
	result.resultKind != BOBOL_LOD_RESULT_FULL_DETAIL &&
	result.resultKind != BOBOL_LOD_RESULT_AABB &&
	result.resultKind != BOBOL_LOD_RESULT_PROXY)
	return FALSE;

    /* Do not let an out-of-order older bootstrap erase a current demand's
     * known failure.  A successful result for the exact demand does prove
     * that the suppression is no longer needed. */
    const auto sourceFailures =
	this->cadOccurrenceFailures.find(sourceBindingKey);
    if (sourceFailures != this->cadOccurrenceFailures.end()) {
	const auto failed = sourceFailures->second.find(occurrenceKey);
	if (failed != sourceFailures->second.end()) {
	    const CadOccurrenceFailure &failure = failed->second;
	    if (failure.databaseRevision ==
		    result.request.databaseRevision &&
		failure.sourceRevision == result.request.sourceRevision &&
		failure.sourceContentHash ==
		    result.request.sourceContentHash &&
		failure.viewRevision == result.request.viewRevision &&
		failure.policyRevision == result.request.policyRevision &&
		failure.requestedCut == result.request.requestedCut &&
		failure.drawMode == result.request.drawMode &&
		failure.qualityTier == result.request.qualityTier) {
		sourceFailures->second.erase(failed);
		if (sourceFailures->second.empty())
		    this->cadOccurrenceFailures.erase(sourceFailures);
	    }
	}
    }

    /*
     * Build the candidate privately, then publish it.  A compact occurrence
     * refines the same authoritative slot many times; when its asset path is
     * unchanged, reuse that slot's CadPayload allocation and move the new
     * small metadata/shared-geometry handles into it.  All alias and asset
     * indexes either resolve the source/occurrence map or already hold the
     * same shared pointer, so pointer identity remains stable while the
     * immutable geometry generation changes atomically between frames.
     */
    CadPayload candidate;
    CadPayload *payload = &candidate;
    payload->progressiveMesh = result.progressiveMesh;
    payload->preparedCadGeometry = result.preparedCadGeometry;
    payload->preparedCadGeometryRevision =
	result.preparedCadGeometryRevision;
    if (consume) {
	payload->mesh = std::move(result.mesh);
	payload->proxy = std::move(result.proxy);
    } else {
	payload->mesh = result.mesh;
	payload->proxy = result.proxy;
    }
    payload->sourcePath = result.request.objectPath.getLength() > 0 ?
	result.request.objectPath : source->path.getValue();
    payload->sourceName = result.request.objectName.getLength() > 0 ?
	result.request.objectName :
	SbString(view_lod_leaf_name_from_path(payload->sourcePath));
    payload->sourceIdentity =
	source->realizationIdentity.getValue().getLength() > 0 ?
	source->realizationIdentity.getValue() : source->path.getValue();
    /* This field intentionally means a compact occurrence key, not the
     * source-node key.  An empty value is a source-wide legacy binding that
     * may fall back to a path; a nonempty value must never alias a sibling. */
    payload->sourceInstanceKey = result.request.occurrenceKey;
    payload->sourceBindingKey = sourceBindingKey.c_str();
    payload->sourceRoutingId = result.request.sourceRoutingId != 0 ?
	result.request.sourceRoutingId.value() :
	source->getCompactSourceRoutingId();
    payload->sourceEntryIndex = result.request.sourceEntryIndex;
    payload->cacheIdentity = result.cacheKey.value;
    payload->cacheKey = result.geometry.cacheKey.isValid() ?
	result.geometry.cacheKey.value : result.cacheKey.value;
    payload->sourceContentHash = result.request.sourceContentHash;
    payload->databaseRevision = result.request.databaseRevision.value();
    payload->sourceRevision = result.request.sourceRevision.value();
    payload->resultKind = result.resultKind;
    payload->qualityTier = result.qualityTier;
    payload->providerStatus = result.providerStatus;
    payload->drawMode = result.request.drawMode;
    payload->activeCut = result.geometry.activeCut;
    payload->residentCut = result.residentCut;
    payload->requestedCut = result.resolvedCut >= 0 ?
	result.resolvedCut : result.request.requestedCut;
    payload->requiredChunks = result.request.requiredChunks;
    payload->presentedChunks = result.request.requiredChunks;
    payload->projectedPixelDiameter =
	result.request.projectedPixelDiameter;
    payload->projectedPixelArea = result.request.projectedPixelArea;
    payload->projectedPixelPerimeter =
	result.request.projectedPixelPerimeter;
    payload->projectedBoundsContained =
	result.request.projectedBoundsContained;
    payload->targetPixelError = result.request.targetPixelError;
    payload->residentAdmissionRevision =
	result.residentAdmissionRevision;
    payload->viewRevision = result.request.viewRevision.value();
    payload->policyRevision = result.request.policyRevision.value();
    payload->visualEmphasis = result.request.visualEmphasis;
    payload->counts = result.counts;
    payload->bounds = result.bounds;
    payload->hasSnappedPoints = result.hasSnappedPoints;
    payload->hasNormals = result.hasNormals;
    payload->shadedCullBackfaces = result.shadedCullBackfaces;
    payload->memoryLimited = result.memoryLimited;
    payload->diagnostic = result.diagnostic;
    (void)view_lod_update_projected_cut_counts(payload, result.request);
    /* Preserve provider accounting for an ordinary whole-prefix result.
     * In particular, it may carry an authored wire line count or byte count
     * which the compact PoP census does not reconstruct.  Only a clipped
     * occurrence needs its provider counts replaced by the retained visible
     * population. */
    if (payload->projectedCutCounts && payload->activeCut >= 0)
	payload->counts = (*payload->projectedCutCounts)[
	    static_cast<size_t>(payload->activeCut)];

    std::unordered_map<std::string, CadPayloadPtr> &sourcePayloads =
	this->cadSourceBindings[sourceBindingKey];
    const auto current = sourcePayloads.find(occurrenceKey);

    /*
     * Provider results publish immutable geometry and residency; they do not
     * own the mutable view demand.  A wheel/trackpad stream can retarget an
     * occurrence several times while one cumulative PoP suffix is in flight.
     * The service deliberately coalesces those requests by asset, so the
     * completing task can carry an older (or merely less exact) requested
     * cut even after the controller rebases its camera epochs.  Replacing the
     * payload wholesale in that state briefly made the older request look
     * satisfied.  No refinement replay was then required until another input
     * event happened to wake the planner.
     *
     * Preserve the owner-thread demand whenever the retained occurrence is
     * at least as current and the result extends the same immutable asset.
     * The newly prepared renderer generation is still installed below.  Its
     * spatial population census is invalidated because the shared generation
     * may have grown; the next bounded demand pass recomputes it from the
     * retained current projection rather than using the worker's stale one.
     */
    if (current != sourcePayloads.end() && current->second) {
	const CadPayload *retained = current->second.get();
	const bool sameAsset = retained->progressiveMesh &&
	    retained->progressiveMesh == payload->progressiveMesh &&
	    retained->databaseRevision == payload->databaseRevision &&
	    retained->sourceContentHash == payload->sourceContentHash &&
	    (payload->sourceContentHash != 0 ||
	     retained->sourceRevision == payload->sourceRevision) &&
	    retained->drawMode == payload->drawMode;
	const bool retainedDemandIsNewer =
	    retained->viewRevision >= payload->viewRevision &&
	    retained->policyRevision >= payload->policyRevision &&
	    (retained->viewRevision > payload->viewRevision ||
	     retained->policyRevision > payload->policyRevision);
	/* The worker may resolve its hierarchy selection more coarsely than a
	 * later coalesced request even after the result has been epoch-rebased.
	 * In an otherwise identical epoch, preserve only that demonstrably richer
	 * owner demand.  Treating all equal-epoch payloads as stale would suppress
	 * a legitimate current provider result after a source/spatial update. */
	const bool retainedCoalescedDemandIsRicher =
	    retained->viewRevision == payload->viewRevision &&
	    retained->policyRevision == payload->policyRevision &&
	    retained->requestedCut > payload->requestedCut;
	if (sameAsset &&
	    (retainedDemandIsNewer || retainedCoalescedDemandIsRicher)) {
	    if (getenv("BOBOL_LOD_TRACE_BUDGET"))
		bu_log("BObol cumulative result retained current demand "
		       "object=%s active=%d/%d requested=%d/%d "
		       "view=%llu/%llu policy=%llu/%llu\n",
		       retained->sourceName.getString(), retained->activeCut,
		       payload->activeCut, retained->requestedCut,
		       payload->requestedCut,
		       static_cast<unsigned long long>(retained->viewRevision),
		       static_cast<unsigned long long>(payload->viewRevision),
		       static_cast<unsigned long long>(retained->policyRevision),
		       static_cast<unsigned long long>(payload->policyRevision));
	    payload->requestedCut = retained->requestedCut;
	    payload->requiredChunks = retained->requiredChunks;
	    payload->projectedPixelDiameter =
		retained->projectedPixelDiameter;
	    payload->projectedPixelArea = retained->projectedPixelArea;
	    payload->projectedPixelPerimeter =
		retained->projectedPixelPerimeter;
	    payload->projectedBoundsContained =
		retained->projectedBoundsContained;
	    payload->targetPixelError = retained->targetPixelError;
	    payload->viewRevision = retained->viewRevision;
	    payload->policyRevision = retained->policyRevision;
	    payload->visualEmphasis = retained->visualEmphasis;
	    payload->projectedCutCounts.reset();
	    payload->projectedCutCountsViewRevision = 0;
	    payload->projectedCutCountsPolicyRevision = 0;
	    payload->projectedCutCountsMeshRevision = 0;
	    /* The old immutable presentation remains authoritative until the
	     * richer generation is committed.  Recomputing its counts against a
	     * concurrently grown/compacted shared asset can observe no common
	     * chunk frontier and transiently publish a zero-cost scene.  Preserve
	     * the exact completed-frame census when the cut is unchanged; for a
	     * genuinely new cut retain the provider's conservative result census.
	     * The next bounded projection pass replaces either with current
	     * view-local counts. */
	    if (payload->activeCut == retained->activeCut)
		payload->counts = retained->counts;
	}
    }
    if (current != sourcePayloads.end() && current->second &&
	view_lod_cad_payload_rank(*current->second) >
	    view_lod_cad_payload_rank(*payload))
	return TRUE;
    CadPayloadPtr publishedPayload;
    if (current != sourcePayloads.end() && current->second) {
	/* Provider results replace renderer generations, not the owner-thread
	 * scene allocation.  Preserve the latter across an asynchronous resident
	 * suffix publication; its epoch checks make a stale allocation inert. */
	payload->allocatedCut = current->second->allocatedCut;
	payload->allocationDrawMode = current->second->allocationDrawMode;
	payload->allocationViewRevision =
	    current->second->allocationViewRevision;
	payload->allocationPolicyRevision =
	    current->second->allocationPolicyRevision;
	payload->allocationMeshRevision =
	    current->second->allocationMeshRevision;
	payload->allocationPlanSerial =
	    current->second->allocationPlanSerial;
	if (current->second->sourceEntryIndex != UINT32_MAX) {
	    const auto indexedSource = this->cadSourceEntryBindings.find(
		current->second->sourceRoutingId);
	    const size_t oldIndex = static_cast<size_t>(
		current->second->sourceEntryIndex);
	    if (indexedSource != this->cadSourceEntryBindings.end() &&
		oldIndex < indexedSource->second.size() &&
		indexedSource->second[oldIndex] == current->second.get())
		indexedSource->second[oldIndex] = NULL;
	}
	this->removeCadPayloadMetrics(current->second.get());
	const bool reuseCompactSlot =
	    result.request.occurrenceKey.getLength() > 0 &&
	    view_lod_paths_equal(current->second->sourcePath,
		payload->sourcePath);
	if (reuseCompactSlot) {
	    publishedPayload = current->second;
	    *publishedPayload = std::move(candidate);
	}
    }
    if (!publishedPayload)
	publishedPayload =
	    std::make_shared<CadPayload>(std::move(candidate));
    payload = publishedPayload.get();
    payload->ownerState = this;

    /*
     * Compact occurrences have one authoritative source/occurrence slot.
     * They are resolved directly by every CAD planner and result consumer;
     * adding the same payload to cadBindings used to allocate another long
     * string/hash record per leaf and retained obsolete result aliases across
     * refinement generations.  Keep aliases only for source-wide legacy
     * results which genuinely need path/name lookup.
     */
    if (payload->sourceInstanceKey.getLength() == 0) {
	view_lod_remove_superseded_cad_payloads(this->cadBindings,
	    sourceBindingKey, payload->sourceInstanceKey, *payload);
	std::vector<std::string> resultKeys;
	view_lod_result_keys(resultKeys, result);
	for (const std::string &resultKey : resultKeys) {
	    const std::string key =
		view_lod_cad_result_binding_key(sourceBindingKey, resultKey);
	    this->cadBindings[key] = publishedPayload;
	}
    }
    sourcePayloads[occurrenceKey] = publishedPayload;
    if (payload->sourceEntryIndex != UINT32_MAX) {
	std::vector<CadPayload *> &indexed =
	    this->cadSourceEntryBindings[payload->sourceRoutingId];
	const size_t entryIndex = static_cast<size_t>(
	    payload->sourceEntryIndex);
	if (indexed.size() <= entryIndex)
	    indexed.resize(entryIndex + 1, NULL);
	indexed[entryIndex] = payload;
    }
    this->addCadPayloadMetrics(payload);
    if (payload->progressiveMesh && payload->progressiveMesh->isValid() &&
	payload->sourcePath.getLength() > 0) {
	std::unordered_map<std::string, CadPayloadPtr> &sourceAssets =
	    this->cadAssetBindings[sourceBindingKey];
	const std::string assetKey = payload->sourcePath.getString();
	const auto resident = sourceAssets.find(assetKey);
	if (resident == sourceAssets.end() || !resident->second ||
	    bu_strcmp(resident->second->cacheKey.getString(),
		payload->cacheKey.getString()) != 0 ||
	    resident->second->residentCut <= payload->residentCut)
	    sourceAssets[assetKey] = publishedPayload;
    }

    /* Preserve the source-wide lookup for a root result.  Entry results must
     * coexist, so binding them under the root keys would make the last
     * completed occurrence replace every sibling. */
	if (result.request.occurrenceKey.getLength() == 0 &&
	(result.request.objectPath.getLength() == 0 ||
	view_lod_paths_equal(source->path.getValue(),
	    result.request.objectPath))) {
	std::vector<std::string> sourceKeys;
	view_lod_source_keys(sourceKeys, source);
	for (const std::string &key : sourceKeys)
	    this->cadBindings[key] = publishedPayload;
    }

    if (payload->sourceInstanceKey.getLength() > 0)
	this->noteCadOccurrenceChanged(
	    sourceBindingKey, payload->sourceInstanceKey);
    else
	this->noteResidentMeshesChanged("source-wide-result");
    return TRUE;
}

const BObolViewLodState::MeshPayload *
BObolViewLodState::findMesh(const SoBRLMeshShape *shape) const
{
    if (this->meshBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, MeshPayloadPtr>::const_iterator it =
	    this->meshBindings.find(keys[i]);
	if (it != this->meshBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::MeshPayload *
BObolViewLodState::findMeshForResult(
    const BObolLodResult &result) const
{
    if (this->meshBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, MeshPayloadPtr>::const_iterator it =
	    this->meshBindings.find(keys[i]);
	if (it != this->meshBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::ProxyPayload *
BObolViewLodState::findProxy(const SoBRLMeshShape *shape) const
{
    if (this->proxyBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_shape_keys(keys, shape);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, ProxyPayloadPtr>::const_iterator it =
	    this->proxyBindings.find(keys[i]);
	if (it != this->proxyBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::ProxyPayload *
BObolViewLodState::findProxyForResult(
    const BObolLodResult &result) const
{
    if (this->proxyBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_result_keys(keys, result);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, ProxyPayloadPtr>::const_iterator it =
	    this->proxyBindings.find(keys[i]);
	if (it != this->proxyBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCad(const SoBRLDatabaseSource *source) const
{
    if (this->cadBindings.empty())
	return NULL;

    std::vector<std::string> keys;
    view_lod_source_keys(keys, source);
    for (size_t i = 0; i < keys.size(); i++) {
	std::unordered_map<std::string, CadPayloadPtr>::const_iterator it =
	    this->cadBindings.find(keys[i]);
	if (it != this->cadBindings.end() && it->second &&
	    it->second->isValid())
	    return it->second.get();
    }

    return NULL;
}

void
BObolViewLodState::findCadPayloadsUnordered(
    const SoBRLDatabaseSource *source,
    std::vector<const CadPayload *> &payloads) const
{
    payloads.clear();
    if (!source || this->cadSourceBindings.empty())
	return;

    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto sourcePayloads = this->cadSourceBindings.find(sourceKey);
    if (sourcePayloads == this->cadSourceBindings.end())
	return;

    payloads.reserve(sourcePayloads->second.size());
    for (const auto &binding : sourcePayloads->second) {
	const CadPayloadPtr &payload = binding.second;
	/*
	 * cadSourceBindings is the authoritative, validated occurrence table.
	 * Payloads enter it only through applySourceResultInternal after complete
	 * validation and are erased rather than invalidated.  Revalidating every
	 * progressive mesh here performs an atomic shared-generation retain for
	 * every occurrence in a full scene synchronization.
	 */
	if (payload)
	    payloads.push_back(payload.get());
    }
}

void
BObolViewLodState::findCadPayloads(
    const SoBRLDatabaseSource *source,
    std::vector<const CadPayload *> &payloads) const
{
    this->findCadPayloadsUnordered(source, payloads);
    std::sort(payloads.begin(), payloads.end(),
	[](const CadPayload *a, const CadPayload *b) {
	    const int instanceOrder = bu_strcmp(a->sourceInstanceKey.getString(),
		b->sourceInstanceKey.getString());
	    if (instanceOrder != 0)
		return instanceOrder < 0;
	    const int pathOrder = bu_strcmp(a->sourcePath.getString(),
		b->sourcePath.getString());
	    if (pathOrder != 0)
		return pathOrder < 0;
	    return bu_strcmp(a->cacheKey.getString(),
		b->cacheKey.getString()) < 0;
	});
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForOccurrence(
    const SoBRLDatabaseSource *source,
    const SbString &occurrenceKey) const
{
    if (!source || occurrenceKey.getLength() == 0 ||
	this->cadSourceBindings.empty())
	return NULL;
    const auto sourcePayloads = this->cadSourceBindings.find(
	view_lod_source_primary_key(source));
    if (sourcePayloads == this->cadSourceBindings.end())
	return NULL;
    const auto payload =
	sourcePayloads->second.find(occurrenceKey.getString());
    return payload != sourcePayloads->second.end() && payload->second ?
	payload->second.get() : NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForSourceEntry(
    const SoBRLDatabaseSource *source, uint32_t sourceEntryIndex,
    const SbString &occurrenceKey) const
{
    if (!source || sourceEntryIndex == UINT32_MAX ||
	occurrenceKey.getLength() == 0 ||
	this->cadSourceEntryBindings.empty())
	return NULL;
    const auto sourceEntries = this->cadSourceEntryBindings.find(
	source->getCompactSourceRoutingId());
    if (sourceEntries == this->cadSourceEntryBindings.end() ||
	static_cast<size_t>(sourceEntryIndex) >= sourceEntries->second.size())
	return NULL;
    const CadPayload *payload = sourceEntries->second[
	static_cast<size_t>(sourceEntryIndex)];
    if (!payload || payload->sourceEntryIndex != sourceEntryIndex ||
	bu_strcmp(payload->sourceInstanceKey.getString(),
	    occurrenceKey.getString()) != 0)
	return NULL;
    return payload;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForAsset(
    const SoBRLDatabaseSource *source,
    const SbString &assetPath) const
{
    if (!source || assetPath.getLength() == 0 ||
	this->cadAssetBindings.empty())
	return NULL;
    const auto sourceAssets = this->cadAssetBindings.find(
	view_lod_source_primary_key(source));
    if (sourceAssets == this->cadAssetBindings.end())
	return NULL;
    const auto payload = sourceAssets->second.find(assetPath.getString());
    return payload != sourceAssets->second.end() && payload->second &&
	payload->second->progressiveMesh &&
	payload->second->progressiveMesh->isValid() ?
	payload->second.get() : NULL;
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForResult(
    const SoBRLDatabaseSource *source,
    const BObolLodResult &result) const
{
    if (!source || this->cadSourceBindings.empty())
	return NULL;

    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto sourcePayloads = this->cadSourceBindings.find(sourceKey);
    if (sourcePayloads == this->cadSourceBindings.end())
	return NULL;
    const std::string occurrenceKey =
	view_lod_cad_occurrence_key(result.request.occurrenceKey);
    const auto payload = sourcePayloads->second.find(occurrenceKey);
    if (payload == sourcePayloads->second.end() || !payload->second ||
	!payload->second->isValid())
	return NULL;

    if (result.request.occurrenceKey.getLength() > 0)
	return payload->second.get();
    if (view_lod_paths_equal(payload->second->sourcePath,
	    result.request.objectPath) ||
	(payload->second->sourceName.getLength() > 0 &&
	 result.request.objectName.getLength() > 0 &&
	 bu_strcmp(payload->second->sourceName.getString(),
	    result.request.objectName.getString()) == 0))
	return payload->second.get();
    return NULL;
}

SbBool
BObolViewLodState::hasCadOccurrenceTerminalFailure(
    const SoBRLDatabaseSource *source,
    const BObolLodRequest &request) const
{
    if (!source)
	return FALSE;
    const auto sourceFailures = this->cadOccurrenceFailures.find(
	view_lod_source_primary_key(source));
    if (sourceFailures == this->cadOccurrenceFailures.end())
	return FALSE;
    const auto failed = sourceFailures->second.find(
	view_lod_cad_occurrence_key(request.occurrenceKey));
    if (failed == sourceFailures->second.end())
	return FALSE;
    const CadOccurrenceFailure &failure = failed->second;
    const SbBool matches = failure.databaseRevision == request.databaseRevision &&
	failure.sourceRevision == request.sourceRevision &&
	failure.sourceContentHash == request.sourceContentHash &&
	failure.viewRevision == request.viewRevision &&
	failure.policyRevision == request.policyRevision &&
	/* A cold traversal intentionally carries -1 until the provider resolves
	 * physical pixel demand against the retained hierarchy.  The failure
	 * stores that resolved ordinal.  Equal view/policy/source epochs make the
	 * unresolved request the same demand; requiring -1 == resolvedCut caused
	 * an otherwise terminal miss to be resubmitted every quiet frame. */
	(request.requestedCut < 0 ||
	 failure.requestedCut == request.requestedCut) &&
	failure.drawMode == request.drawMode &&
	failure.qualityTier == request.qualityTier &&
	view_lod_provider_status_is_terminal_failure(
	    failure.providerStatus) ? TRUE : FALSE;
    if (!matches && getenv("BOBOL_CHUNK_DEBUG"))
	bu_log("BObol terminal failure mismatch db=%llu/%llu source=%llu/%llu "
	    "hash=%llu/%llu view=%llu/%llu policy=%llu/%llu cut=%d/%d "
	    "mode=%d/%d tier=%d/%d status=%d\n",
	    static_cast<unsigned long long>(failure.databaseRevision),
	    static_cast<unsigned long long>(request.databaseRevision.value()),
	    static_cast<unsigned long long>(failure.sourceRevision),
	    static_cast<unsigned long long>(request.sourceRevision.value()),
	    static_cast<unsigned long long>(failure.sourceContentHash),
	    static_cast<unsigned long long>(request.sourceContentHash),
	    static_cast<unsigned long long>(failure.viewRevision),
	    static_cast<unsigned long long>(request.viewRevision.value()),
	    static_cast<unsigned long long>(failure.policyRevision),
	    static_cast<unsigned long long>(request.policyRevision.value()),
	    failure.requestedCut, request.requestedCut,
	    failure.drawMode, request.drawMode,
	    failure.qualityTier, request.qualityTier,
	    failure.providerStatus);
    return matches;
}

size_t
BObolViewLodState::cadOccurrenceTerminalFailureCount(void) const
{
    size_t count = 0;
    for (const auto &source : this->cadOccurrenceFailures)
	count = source.second.size() > SIZE_MAX - count ?
	    SIZE_MAX : count + source.second.size();
    return count;
}

size_t
BObolViewLodState::cadOccurrenceTerminalFailureCountForSource(
    const SoBRLDatabaseSource *source) const
{
    if (!source)
	return 0;
    const auto failures = this->cadOccurrenceFailures.find(
	view_lod_source_primary_key(source));
    return failures == this->cadOccurrenceFailures.end() ?
	0 : failures->second.size();
}

const BObolViewLodState::CadPayload *
BObolViewLodState::findCadForResult(
    const BObolLodResult &result) const
{
    if (this->cadSourceBindings.empty())
	return NULL;

    for (const auto &source : this->cadSourceBindings) {
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->isValid())
		continue;
	if (result.request.occurrenceKey.getLength() > 0) {
	    if (bu_strcmp(payload->sourceInstanceKey.getString(),
		result.request.occurrenceKey.getString()) == 0)
		return payload.get();
	    continue;
	}
	if (view_lod_paths_equal(payload->sourcePath,
		result.request.objectPath) ||
	    (payload->sourceName.getLength() > 0 &&
	     result.request.objectName.getLength() > 0 &&
	     bu_strcmp(payload->sourceName.getString(),
		result.request.objectName.getString()) == 0))
	    return payload.get();
	}
    }

    return NULL;
}

SbBool
BObolViewLodState::retargetMeshPayload(
    const BObolViewLodState::MeshPayload *target,
    int activeCut,
    int requestedCut,
    uint64_t viewRevision,
    uint64_t policyRevision)
{
    if (!target || !target->progressiveMesh ||
	!view_lod_progressive_can_draw(target->progressiveMesh,
	    target->requiredChunks, activeCut) ||
	requestedCut < 0)
	return FALSE;

    MeshPayloadPtr payload;
    for (const auto &binding : this->meshBindings) {
	if (binding.second && binding.second.get() == target) {
	    payload = binding.second;
	    break;
	}
    }
    if (!payload)
	return FALSE;

    payload->activeCut = activeCut;
    payload->requestedCut = requestedCut;
    /* residentCut is the last population-bearing PoP cut published for
     * this payload.  A coordinate-only plateau may make a finer activeCut
     * drawable without changing that population high-water mark, so do not
     * overwrite it from the selected cut or from a newer worker generation. */
    payload->residentAdmissionRevision = 0;
    payload->memoryLimited = FALSE;
    payload->viewRevision = viewRevision;
    payload->policyRevision = policyRevision;
    payload->counts = view_lod_progressive_counts(
	payload->progressiveMesh, payload->requiredChunks,
	activeCut, payload->hasNormals);
    return TRUE;
}

SbBool
BObolViewLodState::retargetCadPayload(
    const BObolViewLodState::CadPayload *target,
    int activeCut,
    const BObolLodRequest &demand)
{
    const int requestedCut = demand.requestedCut;
    const bool spatialPageSetChanged = target &&
	!demand.requiredChunks.empty() &&
	target->presentedChunks != demand.requiredChunks;
    /* The shared progressive asset may already contain a worker-built
     * generation whose result has not reached this view.  The immutable
     * prepared PartGeometry is the authoritative presentation snapshot; its
     * progressive channel frontier also preserves coordinate-only cuts
     * above the last population-bearing residentCut. */
    const auto preparedChannelCut = [target](bool wire) {
	if (!target || !target->progressiveMesh ||
	    !target->preparedCadGeometry)
	    return -1;
	if (wire) {
	    if (!target->preparedCadGeometry->wire)
		return -1;
	    const Obol::WireRep &prepared =
		*target->preparedCadGeometry->wire;
	    return prepared.isProgressive() ?
		static_cast<int>(prepared.progressiveResidentCut) :
		target->progressiveMesh->maximumCut();
	}
	if (!target->preparedCadGeometry->shaded)
	    return -1;
	const Obol::TriMesh &prepared =
	    *target->preparedCadGeometry->shaded;
	return prepared.isProgressive() ?
	    static_cast<int>(prepared.progressiveResidentCut) :
	    target->progressiveMesh->maximumCut();
    };
    /* The immutable progressive asset may be retained across a representation
     * switch.  Validate the channel demanded now, not the channel which last
     * published the payload. */
    const bool needsWire = target &&
	(demand.drawMode == BOBOL_LOD_DRAW_WIRE ||
	 demand.drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE);
    const bool needsShaded = target &&
	(demand.drawMode == BOBOL_LOD_DRAW_SHADED ||
	 demand.drawMode == BOBOL_LOD_DRAW_SHADED_BOTS ||
	 demand.drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE);
    int preparedDrawableCut = target ? target->residentCut : -1;
    if (target && target->preparedCadGeometry) {
	if (needsWire)
	    preparedDrawableCut = preparedChannelCut(true);
	if (needsShaded) {
	    const int shadedCut = preparedChannelCut(false);
	    preparedDrawableCut = needsWire ?
		std::min(preparedDrawableCut, shadedCut) : shadedCut;
	}
    }

    const bool preparedGenerationCurrent = target &&
	target->progressiveMesh &&
	target->preparedCadGeometry &&
	target->preparedCadGeometryRevision ==
	    target->progressiveMesh->revision();
    const bool demandActiveCutDrawable = target &&
	view_lod_progressive_can_draw(target->progressiveMesh,
	    demand.requiredChunks, activeCut);
    const bool canAdoptDemandPageSet = !spatialPageSetChanged ||
	(preparedGenerationCurrent && demandActiveCutDrawable);

    if (!target || !target->progressiveMesh ||
	/* A representation-only retarget does not expose a new cut.  The CAD
	 * assembly owns preparation of the newly selected renderer channel and may
	 * already have cached it on the shared progressive generation even though
	 * this occurrence's historical preparedCadGeometry handle names the prior
	 * channel.  Require a prepared frontier only when the cut itself advances. */
	(activeCut != target->activeCut &&
	 activeCut > preparedDrawableCut) ||
	(activeCut != target->activeCut &&
	 (!canAdoptDemandPageSet || !demandActiveCutDrawable)) ||
	requestedCut < 0)
	return FALSE;

    CadPayloadPtr payload;
    const auto source = this->cadSourceBindings.find(
	target->sourceBindingKey.getString());
    if (source != this->cadSourceBindings.end()) {
	const auto occurrence = source->second.find(
	    view_lod_cad_occurrence_key(target->sourceInstanceKey));
	if (occurrence != source->second.end() &&
	    occurrence->second.get() == target)
	    payload = occurrence->second;
    }
    if (!payload)
	return FALSE;

    /* A scene allocation is valid only for the exact demand epoch and
     * renderer channel it priced.  The controller installs a matching new
     * allocation before its bounded application pass; ordinary camera/policy
     * retargets invalidate an older decision without touching presentation. */
    if (payload->allocationViewRevision != demand.viewRevision.value() ||
	payload->allocationPolicyRevision != demand.policyRevision.value() ||
	payload->allocationDrawMode != demand.drawMode) {
	payload->allocatedCut = -1;
	payload->allocationDrawMode = BOBOL_LOD_DRAW_UNKNOWN;
	payload->allocationViewRevision = 0;
	payload->allocationPolicyRevision = 0;
	payload->allocationMeshRevision = 0;
	payload->allocationPlanSerial = 0;
    }
    const bool projectedPopulationChanged =
	view_lod_update_projected_cut_counts(payload.get(), demand);

    /*
     * A camera-only demand update often changes requestedCut/revisions
     * while deliberately retaining the current drawable prefix.  That is
     * metadata, not a presentation mutation: keep aggregate face/resident
     * metrics and the assembly occurrence journal untouched.  Publishing a
     * false geometry change here made every pose-only motion frame repack
     * hundreds of otherwise identical occurrences.
     */
    if (payload->activeCut == activeCut) {
	/* drawMode is render-cost currency, not only descriptive result
	 * provenance.  A progressive shaded/wire switch needs no provider result,
	 * so update the aggregate metrics in place even when the cut itself is
	 * unchanged. */
	if (payload->drawMode != demand.drawMode) {
	    const BObolLodCounts priorPresentationCounts = payload->counts;
	    this->removeCadPayloadMetrics(payload.get());
	    payload->drawMode = demand.drawMode;
	    payload->requestedCut = requestedCut;
	    payload->requiredChunks = demand.requiredChunks;
	    if (canAdoptDemandPageSet)
		payload->presentedChunks = demand.requiredChunks;
	    payload->projectedPixelDiameter = demand.projectedPixelDiameter;
	    payload->projectedPixelArea = demand.projectedPixelArea;
	    payload->projectedPixelPerimeter = demand.projectedPixelPerimeter;
	    payload->projectedBoundsContained =
		demand.projectedBoundsContained;
	    payload->targetPixelError = demand.targetPixelError;
	    payload->residentAdmissionRevision = 0;
	    payload->memoryLimited = FALSE;
	    payload->viewRevision = demand.viewRevision.value();
	    payload->policyRevision = demand.policyRevision.value();
	    payload->visualEmphasis = demand.visualEmphasis;
	    payload->counts = canAdoptDemandPageSet ?
		view_lod_cad_counts_at_cut(payload.get(), activeCut) :
		priorPresentationCounts;
	    this->addCadPayloadMetrics(payload.get());
	    return TRUE;
	}
	const bool wasSatisfied =
	    view_lod_cad_payload_is_satisfied(payload.get());
	const bool wasMemoryLimited = payload->memoryLimited;
	/* requestedCut contributes to the retained view demand even when the
	 * presented cut is unchanged.  Update that O(1) aggregate in place. */
	const bool populationChanged = projectedPopulationChanged ||
	    payload->requiredChunks != demand.requiredChunks;
	const BObolLodCounts priorPresentationCounts = payload->counts;
	if (populationChanged)
	    this->removeCadPayloadMetrics(payload.get());
	else
	    this->removeCadResidentDemand(payload.get());
	payload->requestedCut = requestedCut;
	payload->requiredChunks = demand.requiredChunks;
	if (canAdoptDemandPageSet)
	    payload->presentedChunks = demand.requiredChunks;
	payload->projectedPixelDiameter = demand.projectedPixelDiameter;
	payload->projectedPixelArea = demand.projectedPixelArea;
	payload->projectedPixelPerimeter = demand.projectedPixelPerimeter;
	payload->projectedBoundsContained = demand.projectedBoundsContained;
	payload->targetPixelError = demand.targetPixelError;
	payload->residentAdmissionRevision = 0;
	payload->memoryLimited = FALSE;
	payload->viewRevision = demand.viewRevision.value();
	payload->policyRevision = demand.policyRevision.value();
	payload->visualEmphasis = demand.visualEmphasis;
	if (populationChanged) {
	    payload->counts = canAdoptDemandPageSet ?
		view_lod_cad_counts_at_cut(payload.get(), activeCut) :
		priorPresentationCounts;
	    this->addCadPayloadMetrics(payload.get());
	} else {
	    const bool isSatisfied =
		view_lod_cad_payload_is_satisfied(payload.get());
	    if (wasSatisfied != isSatisfied) {
		if (isSatisfied)
		    this->cadSatisfiedMeshPayloadCount =
			view_lod_saturating_add(
			    this->cadSatisfiedMeshPayloadCount, 1);
		else
		    this->cadSatisfiedMeshPayloadCount =
			view_lod_saturating_subtract(
			    this->cadSatisfiedMeshPayloadCount, 1);
	    }
	    if (payload->sourceInstanceKey.getLength() > 0) {
		const std::string sourceKey =
		    payload->sourceBindingKey.getString();
		const std::string occurrenceKey =
		    payload->sourceInstanceKey.getString();
		if (isSatisfied) {
		    const auto unsatisfied =
			this->cadUnsatisfiedOccurrencesBySource.find(sourceKey);
		    if (unsatisfied !=
			    this->cadUnsatisfiedOccurrencesBySource.end()) {
			unsatisfied->second.erase(occurrenceKey);
			if (unsatisfied->second.empty())
			    this->cadUnsatisfiedOccurrencesBySource.erase(
				unsatisfied);
		    }
		} else {
		    this->cadUnsatisfiedOccurrencesBySource[sourceKey].insert(
			occurrenceKey);
		}
	    }
	    this->addCadResidentDemand(payload.get());
	}
	if (!populationChanged && wasMemoryLimited)
	    this->cadMemoryLimitedMeshPayloadCount =
		view_lod_saturating_subtract(
		    this->cadMemoryLimitedMeshPayloadCount, 1);
	return TRUE;
    }

    this->removeCadPayloadMetrics(payload.get());
    payload->activeCut = activeCut;
    payload->drawMode = demand.drawMode;
    payload->requestedCut = requestedCut;
    payload->requiredChunks = demand.requiredChunks;
    payload->presentedChunks = demand.requiredChunks;
    payload->projectedPixelDiameter = demand.projectedPixelDiameter;
    payload->projectedPixelArea = demand.projectedPixelArea;
    payload->projectedPixelPerimeter = demand.projectedPixelPerimeter;
    payload->projectedBoundsContained = demand.projectedBoundsContained;
    payload->targetPixelError = demand.targetPixelError;
    payload->residentAdmissionRevision = 0;
    payload->memoryLimited = FALSE;
    payload->viewRevision = demand.viewRevision.value();
    payload->policyRevision = demand.policyRevision.value();
    payload->visualEmphasis = demand.visualEmphasis;
    payload->counts = view_lod_cad_counts_at_cut(payload.get(), activeCut);
    this->addCadPayloadMetrics(payload.get());
    this->noteCadOccurrenceChanged(
	payload->sourceBindingKey.getString(), payload->sourceInstanceKey);
    return TRUE;
}

SbBool
BObolViewLodState::setCadAllocatedCut(
    const BObolViewLodState::CadPayload *target,
    int allocatedCut,
    uint64_t viewRevision,
    uint64_t policyRevision,
    int drawMode)
{
    if (!target || !target->progressiveMesh ||
	!target->progressiveMesh->isValid() || allocatedCut < 0 ||
	allocatedCut < target->progressiveMesh->minimumCut() ||
	allocatedCut > target->progressiveMesh->maximumCut() ||
	viewRevision == 0 || policyRevision == 0)
	return FALSE;

    if (target->ownerState != this)
	return FALSE;

    CadPayload *payload = const_cast<CadPayload *>(target);
    payload->allocatedCut = allocatedCut;
    payload->allocationDrawMode = drawMode;
    payload->allocationViewRevision = viewRevision;
    payload->allocationPolicyRevision = policyRevision;
    payload->allocationMeshRevision = payload->progressiveMesh->revision();
    payload->allocationPlanSerial = 0;
    return TRUE;
}

uint64_t
BObolViewLodState::beginCadAllocationPlan(void)
{
    this->cadNextAllocationPlanSerial++;
    if (!this->cadNextAllocationPlanSerial)
	this->cadNextAllocationPlanSerial++;
    return this->cadNextAllocationPlanSerial;
}

SbBool
BObolViewLodState::stageCadAllocatedCut(
    const BObolViewLodState::CadPayload *target,
    int allocatedCut,
    uint64_t viewRevision,
    uint64_t policyRevision,
    int drawMode,
    uint64_t planSerial)
{
    if (!target || target->ownerState != this || !target->progressiveMesh ||
	!target->progressiveMesh->isValid() || allocatedCut < 0 ||
	allocatedCut < target->progressiveMesh->minimumCut() ||
	allocatedCut > target->progressiveMesh->maximumCut() ||
	viewRevision == 0 || policyRevision == 0 || planSerial == 0 ||
	planSerial != this->cadNextAllocationPlanSerial)
	return FALSE;

    CadPayload *payload = const_cast<CadPayload *>(target);
    payload->allocatedCut = allocatedCut;
    payload->allocationDrawMode = drawMode;
    payload->allocationViewRevision = viewRevision;
    payload->allocationPolicyRevision = policyRevision;
    payload->allocationMeshRevision = payload->progressiveMesh->revision();
    payload->allocationPlanSerial = planSerial;
    return TRUE;
}

SbBool
BObolViewLodState::commitCadAllocationPlan(
    uint64_t planSerial, uint64_t cadRevision,
    uint64_t residentDemandRevision)
{
    if (!planSerial || planSerial != this->cadNextAllocationPlanSerial ||
	cadRevision != this->cadBindingsRevision ||
	residentDemandRevision != this->cadResidentDemandRevision)
	return FALSE;
    this->cadActiveAllocationPlanSerial = planSerial;
    return TRUE;
}

SbBool
BObolViewLodState::isCadAllocationPlanCurrent(uint64_t planSerial) const
{
    return planSerial != 0 &&
	planSerial == this->cadNextAllocationPlanSerial ? TRUE : FALSE;
}

uint64_t
BObolViewLodState::activeCadAllocationPlan(void) const
{
    return this->cadActiveAllocationPlanSerial;
}

SbBool
BObolViewLodState::cadAllocationPlanCoversCurrentPopulation(
    uint64_t planSerial, uint64_t viewRevision,
    uint64_t policyRevision, size_t fixedCadPresentationCost) const
{
    if (!planSerial || !viewRevision || !policyRevision ||
	planSerial != this->cadActiveAllocationPlanSerial)
	return FALSE;

    size_t fixedCost = 0;
    for (const auto &source : this->cadSourceBindings) {
	for (const auto &occurrence : source.second) {
	    const CadPayload *payload = occurrence.second.get();
	    if (!payload || !payload->isValid() ||
		payload->viewRevision != viewRevision ||
		payload->policyRevision != policyRevision)
		continue;
	    if (!payload->progressiveMesh ||
		!payload->progressiveMesh->isValid()) {
		const size_t cost = bobol_lod_render_cost_units(
		    payload->counts, payload->drawMode, 1);
		fixedCost = cost > SIZE_MAX - fixedCost ? SIZE_MAX :
		    fixedCost + cost;
		continue;
	    }
	    if (payload->allocationPlanSerial != planSerial ||
		payload->allocationViewRevision != viewRevision ||
		payload->allocationPolicyRevision != policyRevision ||
		payload->allocationDrawMode != payload->drawMode ||
		payload->allocationMeshRevision !=
		    payload->progressiveMesh->revision() ||
		payload->allocatedCut < payload->progressiveMesh->minimumCut() ||
		payload->allocatedCut > payload->progressiveMesh->maximumCut())
		return FALSE;
	}
    }
    return fixedCost == fixedCadPresentationCost ? TRUE : FALSE;
}

SbBool
BObolViewLodState::removeCadPayload(
    const BObolViewLodState::CadPayload *target)
{
    if (!target)
	return FALSE;
    SbBool removed = FALSE;
    CadPayloadPtr removedPayload;
    const std::string sourceKey = target->sourceBindingKey.getString();
    const std::string occurrenceKey =
	view_lod_cad_occurrence_key(target->sourceInstanceKey);
    auto source = this->cadSourceBindings.find(sourceKey);
    if (source != this->cadSourceBindings.end()) {
	auto occurrence = source->second.find(occurrenceKey);
	if (occurrence != source->second.end() &&
	    occurrence->second.get() == target) {
	    removedPayload = occurrence->second;
	    if (target->sourceEntryIndex != UINT32_MAX) {
		auto indexedSource = this->cadSourceEntryBindings.find(
		    target->sourceRoutingId);
		const size_t entryIndex = static_cast<size_t>(
		    target->sourceEntryIndex);
		if (indexedSource != this->cadSourceEntryBindings.end() &&
		    entryIndex < indexedSource->second.size() &&
		    indexedSource->second[entryIndex] == target)
		    indexedSource->second[entryIndex] = NULL;
	    }
	    this->removeCadPayloadMetrics(occurrence->second.get());
	    source->second.erase(occurrence);
	    if (source->second.empty()) {
		this->cadSourceBindings.erase(source);
		this->cadSourceEntryBindings.erase(
		    removedPayload->sourceRoutingId);
	    }
	    removed = TRUE;
	}
    }

    const CadPayload *removedTarget = removedPayload ?
	removedPayload.get() : target;
    if (removedTarget->sourceInstanceKey.getLength() == 0) {
	for (auto binding = this->cadBindings.begin();
	     binding != this->cadBindings.end();) {
	    if (binding->second && binding->second.get() == removedTarget)
		binding = this->cadBindings.erase(binding);
	    else
		++binding;
	}
    }
    if (removed) {
	if (removedTarget->sourceInstanceKey.getLength() > 0)
	    this->noteCadOccurrenceChanged(
		sourceKey, removedTarget->sourceInstanceKey);
	else {
	    this->noteResidentMeshesChanged("source-wide-remove");
	    this->clearCadPresentations();
	}
    }
    return removed;
}

SbBool
BObolViewLodState::removeMeshPayload(
    const BObolViewLodState::MeshPayload *target)
{
    if (!target)
	return FALSE;
    SbBool removed = FALSE;
    for (auto binding = this->meshBindings.begin();
	 binding != this->meshBindings.end();) {
	if (binding->second && binding->second.get() == target) {
	    binding = this->meshBindings.erase(binding);
	    removed = TRUE;
	} else {
	    ++binding;
	}
    }
    return removed;
}

size_t
BObolViewLodState::bindingCount(void) const
{
    return this->meshBindings.size() + this->proxyBindings.size() +
	   this->cadBindings.size();
}

static std::vector<BObolViewLodState::MeshPayloadPtr>
view_lod_unique_payloads(
    const std::unordered_map<std::string,
    BObolViewLodState::MeshPayloadPtr> &bindings)
{
    std::vector<BObolViewLodState::MeshPayloadPtr> payloads;
    std::unordered_set<const BObolViewLodState::MeshPayload *> seen;
    seen.reserve(bindings.size());
    for (std::unordered_map<std::string,
	 BObolViewLodState::MeshPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (seen.insert(it->second.get()).second)
	    payloads.push_back(it->second);
    }

    return payloads;
}

static std::vector<BObolViewLodState::ProxyPayloadPtr>
view_lod_unique_proxy_payloads(
    const std::unordered_map<std::string,
    BObolViewLodState::ProxyPayloadPtr> &bindings)
{
    std::vector<BObolViewLodState::ProxyPayloadPtr> payloads;
    std::unordered_set<const BObolViewLodState::ProxyPayload *> seen;
    seen.reserve(bindings.size());
    for (std::unordered_map<std::string,
	 BObolViewLodState::ProxyPayloadPtr>::const_iterator it =
	     bindings.begin(); it != bindings.end(); ++it) {
	if (!it->second)
	    continue;
	if (seen.insert(it->second.get()).second)
	    payloads.push_back(it->second);
    }

    return payloads;
}

size_t
BObolViewLodState::payloadCount(void) const
{
    return view_lod_unique_payloads(this->meshBindings).size() +
	   view_lod_unique_proxy_payloads(this->proxyBindings).size() +
	   this->cadPayloadCount();
}

size_t
BObolViewLodState::meshPayloadCount(void) const
{
    return view_lod_unique_payloads(this->meshBindings).size();
}

size_t
BObolViewLodState::proxyPayloadCount(int proxyKind) const
{
    std::vector<ProxyPayloadPtr> payloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    size_t count = 0;
    for (size_t i = 0; i < payloads.size(); i++) {
	if (!payloads[i] || !payloads[i]->isValid())
	    continue;
	if (proxyKind == BOBOL_LOD_PROXY_NONE ||
	    payloads[i]->proxy.kind == proxyKind)
	    count++;
    }

    return count;
}

size_t
BObolViewLodState::cadPayloadCount(void) const
{
    return this->cadValidPayloadCount;
}

size_t
BObolViewLodState::cadPayloadCountForSource(
    const SoBRLDatabaseSource *source) const
{
    if (!source)
	return 0;
    const auto payloads = this->cadSourceBindings.find(
	view_lod_source_primary_key(source));
    return payloads == this->cadSourceBindings.end() ?
	0 : payloads->second.size();
}

size_t
BObolViewLodState::cadMeshPayloadCount(void) const
{
    return this->cadMeshPayloadCountValue;
}

size_t
BObolViewLodState::cadMeshPayloadCountForSource(
    const SoBRLDatabaseSource *source) const
{
    if (!source)
	return 0;
    const auto count = this->cadMeshPayloadCountsBySource.find(
	view_lod_source_primary_key(source));
    return count == this->cadMeshPayloadCountsBySource.end() ?
	0 : count->second;
}

void
BObolViewLodState::unsatisfiedCadOccurrenceKeys(
    const SoBRLDatabaseSource *source,
    uint64_t residentAdmissionRevision,
    std::vector<SbString> &occurrenceKeys) const
{
    occurrenceKeys.clear();
    if (!source)
	return;
    const auto found = this->cadUnsatisfiedOccurrencesBySource.find(
	view_lod_source_primary_key(source));
    if (found == this->cadUnsatisfiedOccurrencesBySource.end())
	return;
    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto sourcePayloads = this->cadSourceBindings.find(sourceKey);
    occurrenceKeys.reserve(found->second.size());
    for (const std::string &key : found->second) {
	/* A resident-memory denial is authoritative only for the capacity
	 * epoch which produced it.  Omitting that occurrence from the sparse
	 * frontier makes hard pressure a stable terminal condition instead of
	 * an idle rescan loop.  Reclamation or a configured-limit increase
	 * advances the service epoch and makes the occurrence actionable again. */
	if (residentAdmissionRevision &&
	    sourcePayloads != this->cadSourceBindings.end()) {
	    const auto payload = sourcePayloads->second.find(key);
	    if (payload != sourcePayloads->second.end() && payload->second &&
		payload->second->memoryLimited &&
		payload->second->residentAdmissionRevision ==
		    residentAdmissionRevision)
		continue;
	}
	occurrenceKeys.push_back(key.c_str());
    }
}

void
BObolViewLodState::retriableMemoryLimitedCadOccurrenceKeys(
    const SoBRLDatabaseSource *source,
    uint64_t residentAdmissionRevision,
    std::vector<SbString> &occurrenceKeys) const
{
    occurrenceKeys.clear();
    if (!source)
	return;
    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto unsatisfied = this->cadUnsatisfiedOccurrencesBySource.find(
	sourceKey);
    const auto sourcePayloads = this->cadSourceBindings.find(sourceKey);
    if (unsatisfied == this->cadUnsatisfiedOccurrencesBySource.end() ||
	sourcePayloads == this->cadSourceBindings.end())
	return;

    occurrenceKeys.reserve(unsatisfied->second.size());
    for (const std::string &key : unsatisfied->second) {
	const auto payload = sourcePayloads->second.find(key);
	if (payload == sourcePayloads->second.end() || !payload->second ||
	    !payload->second->memoryLimited)
	    continue;
	/* A denial from this exact capacity epoch is terminal until the service
	 * reports another genuine reclamation/limit-growth edge. */
	if (residentAdmissionRevision &&
	    payload->second->residentAdmissionRevision ==
		residentAdmissionRevision)
	    continue;
	occurrenceKeys.push_back(key.c_str());
    }
}

size_t
BObolViewLodState::cadProxyPayloadCount(int proxyKind) const
{
    if (proxyKind == BOBOL_LOD_PROXY_NONE)
	return this->cadProxyPayloadCountValue;
    if (proxyKind < 0 ||
	proxyKind >= static_cast<int>(sizeof(this->cadProxyKindCounts) /
	    sizeof(this->cadProxyKindCounts[0])))
	return 0;
    return this->cadProxyKindCounts[proxyKind];
}

void
BObolViewLodState::convergencePayloadCounts(size_t &active,
	size_t &satisfied, size_t &memoryLimited) const
{
    active = 0;
    satisfied = 0;
    memoryLimited = 0;

    const std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->isValid() ||
	    payload->resultKind != BOBOL_LOD_RESULT_MESH)
	    continue;
	active++;
	if (payload->requestedCut < 0 ||
	    payload->activeCut >= payload->requestedCut)
	    satisfied++;
	if (payload->memoryLimited)
	    memoryLimited++;
    }

    active = view_lod_saturating_add(
	active, this->cadMeshPayloadCountValue);
    satisfied = view_lod_saturating_add(
	satisfied, this->cadSatisfiedMeshPayloadCount);
    memoryLimited = view_lod_saturating_add(
	memoryLimited, this->cadMemoryLimitedMeshPayloadCount);
}

size_t
BObolViewLodState::memoryLimitedPayloadCount(void) const
{
    size_t count = this->cadMemoryLimitedMeshPayloadCount;
    const std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads)
	if (payload && payload->isValid() && payload->memoryLimited)
	    count = view_lod_saturating_add(count, 1);
    return count;
}

size_t
BObolViewLodState::activeFaceCount(void) const
{
    size_t faces = 0;
    const std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->isValid() ||
	    payload->resultKind != BOBOL_LOD_RESULT_MESH)
	    continue;
	faces = view_lod_saturating_add(faces, payload->counts.faceCount);
    }

    return view_lod_saturating_add(faces, this->cadActiveFaceCount);
}

size_t
BObolViewLodState::activeRenderCost(void) const
{
    size_t cost = 0;
    const std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->isValid() ||
	    payload->resultKind != BOBOL_LOD_RESULT_MESH)
	    continue;
	cost = view_lod_saturating_add(cost,
	    bobol_lod_render_cost_units(
		payload->counts, BOBOL_LOD_DRAW_SHADED, 1));
    }
    return view_lod_saturating_add(cost, this->cadActiveRenderCost);
}

size_t
BObolViewLodState::activeCadRenderCost(void) const
{
    return this->cadActiveRenderCost;
}

size_t
BObolViewLodState::minimumActiveRenderCost(void) const
{
    size_t cost = 0;
    const std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->isValid() ||
	    payload->resultKind != BOBOL_LOD_RESULT_MESH)
	    continue;
	BObolLodCounts counts = payload->counts;
	if (payload->progressiveMesh &&
	    payload->progressiveMesh->isValid()) {
	    const int minimumCut =
		payload->progressiveMesh->minimumCut();
	    if (minimumCut >= 0 &&
		view_lod_progressive_can_draw(payload->progressiveMesh,
		    payload->requiredChunks, minimumCut)) {
		counts = view_lod_progressive_counts(
		    payload->progressiveMesh, payload->requiredChunks,
		    minimumCut, payload->hasNormals);
	    }
	}
	cost = view_lod_saturating_add(cost,
	    bobol_lod_render_cost_units(
		counts, BOBOL_LOD_DRAW_SHADED, 1));
    }
    return view_lod_saturating_add(
	cost, this->cadMinimumActiveRenderCost);
}

size_t
BObolViewLodState::minimumActiveCadRenderCost(void) const
{
    return this->cadMinimumActiveRenderCost;
}

SbBool
BObolViewLodState::lastCadPresentedPrimitiveCount(size_t &primitives) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    primitives = this->cadLastPresentedPrimitiveCount;
    return this->cadLastPresentedPrimitiveCountValid;
}

SbBool
BObolViewLodState::lastCadPresentedRenderCost(size_t &cost) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    cost = this->cadLastPresentedRenderCost;
    return this->cadLastPresentedRenderCostValid;
}

SbBool
BObolViewLodState::lastCadPresentationFrameExact(void) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    return this->cadLastPresentationFrameExact;
}

SbBool
BObolViewLodState::lastCadPresentationOccurrenceCoverage(
    size_t &subpixelOccurrences, size_t &structuralBoxes) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    subpixelOccurrences = this->cadLastSubpixelProxyCount;
    structuralBoxes = this->cadLastUncollapsedStructuralProxyCount;
    return this->cadLastPresentationFrameExact;
}

SbBool
BObolViewLodState::lastCadStructuralProjectionHistogram(
    CadStructuralProjectionHistogram &histogram) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    histogram = this->cadLastStructuralProjectionHistogram;
    return histogram.exact;
}

SbBool
BObolViewLodState::lastCadGpuMeasurement(
	size_t &faces, uint64_t &nanoseconds, uint64_t &serial,
	float &pointProxyPixelThreshold) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    faces = this->cadLastGpuFaces;
    nanoseconds = this->cadLastGpuNanoseconds;
    serial = this->cadLastGpuSerial;
    pointProxyPixelThreshold = this->cadLastGpuPointProxyPixelThreshold;
    return this->cadLastGpuMeasurementValid;
}

SbBool
BObolViewLodState::hasCadPresentationAssemblies(void) const
{
    for (const auto &entry : this->cadPresentationAssemblyUseCounts) {
	if (entry.first && entry.second && entry.first->instanceCount())
	    return TRUE;
    }
    return FALSE;
}

uint64_t
BObolViewLodState::cadPresentationExecutionSerial(void) const
{
    uint64_t serial = 0;
    for (const auto &entry : this->cadPresentationAssemblyUseCounts) {
	const SoCADAssembly *assembly = entry.first;
	if (!assembly || !entry.second || !assembly->instanceCount())
	    continue;
	const uint64_t value = assembly->renderExecutionSerial();
	serial = value > UINT64_MAX - serial ? UINT64_MAX : serial + value;
    }
    return serial;
}

uint64_t
BObolViewLodState::cadPresentationPreparationSerial(void) const
{
    uint64_t serial = 0;
    for (const auto &entry : this->cadPresentationAssemblyUseCounts) {
	const SoCADAssembly *assembly = entry.first;
	if (!assembly || !entry.second || !assembly->instanceCount())
	    continue;
	const uint64_t value = assembly->renderPreparationSerial();
	serial = value > UINT64_MAX - serial ? UINT64_MAX : serial + value;
    }
    return serial;
}

void
BObolViewLodState::beginCadPresentationFrame(void) const
{
    this->cadPresentationFrameStartExecutionSerials.clear();
    this->cadPresentationFrameStartExecutionSerials.reserve(
	this->cadPresentationAssemblyUseCounts.size());
    for (const auto &entry : this->cadPresentationAssemblyUseCounts) {
	const SoCADAssembly *assembly = entry.first;
	if (!assembly || !entry.second || !assembly->instanceCount())
	    continue;
	this->cadPresentationFrameStartExecutionSerials.emplace(
	    assembly, assembly->renderExecutionSerial());
    }
    this->cadPresentationFrameObservationArmed = TRUE;
    this->cadPresentationFrameStatusValid = FALSE;
}

void
BObolViewLodState::refreshCadPresentationFrameStatus(void) const
{
    this->cadLastPresentedPrimitiveCount = 0;
    this->cadLastPresentedRenderCost = 0;
    this->cadLastGpuFaces = 0;
    this->cadLastGpuNanoseconds = 0;
    this->cadLastGpuSerial = 1469598103934665603ULL;
    this->cadLastGpuPointProxyPixelThreshold = 1.0f;
    this->cadLastPresentationFrameExact = TRUE;
    this->cadLastSubpixelProxyCount = 0;
    this->cadLastUncollapsedStructuralProxyCount = 0;
    this->cadLastStructuralProjectionHistogram =
	CadStructuralProjectionHistogram();
    this->cadLastStructuralProjectionHistogram.exact = TRUE;
    this->cadLastStructuralProjectionHistogram.revision =
	1469598103934665603ULL;
    this->cadGpuResourceStatusValue = CadGpuResourceStatus();

    SbBool haveAssembly = FALSE;
    SbBool presentedValid = TRUE;
    SbBool presentedRenderCostValid = TRUE;
    SbBool gpuMeasurementValid = TRUE;
    SbBool haveGpuPointProxyPixelThreshold = FALSE;
    SbBool preparedReplay = TRUE;
    SbBool resourceSampled = FALSE;
    uint64_t resourceSerial = 1469598103934665603ULL;

    /* lastRenderedWork() belongs to an assembly, not to the host traversal.
     * If Coin stops between assemblies, an unvisited assembly still reports
     * the preceding frame's exact record.  Require every non-empty assembly
     * present at frame start to have entered its renderer during this frame. */
    if (this->cadPresentationFrameObservationArmed) {
	for (const auto &entry :
		this->cadPresentationFrameStartExecutionSerials) {
	    if (!entry.first ||
		entry.first->renderExecutionSerial() == entry.second) {
		this->cadLastPresentationFrameExact = FALSE;
		break;
	    }
	}
    }

    for (const auto &entry : this->cadPresentationAssemblyUseCounts) {
	const SoCADAssembly *assembly = entry.first;
	/* Streaming and edit transactions may retain a presentation object while
	 * its occurrence set is empty.  Such an assembly is deliberately absent
	 * from the frame-start execution snapshot above and contributes no draw
	 * work, resource sample, or exactness requirement here.  Counting it only
	 * at completion made a valid non-empty frame permanently inexact. */
	if (!assembly || !entry.second || !assembly->instanceCount())
	    continue;
	haveAssembly = TRUE;
	const int tier = assembly->lastRenderTier();
	const Obol::CadRenderedWork work = assembly->lastRenderedWork();
	const size_t subpixel = assembly->lastSubpixelProxyCount();
	const size_t structural =
	    assembly->lastUncollapsedStructuralProxyCount();
	const Obol::CadStructuralProxyProjectionHistogram projection =
	    assembly->lastStructuralProxyProjectionHistogram();
	this->cadLastSubpixelProxyCount = view_lod_saturating_add(
	    this->cadLastSubpixelProxyCount, subpixel);
	this->cadLastUncollapsedStructuralProxyCount =
	    view_lod_saturating_add(
		    this->cadLastUncollapsedStructuralProxyCount, structural);
	if (!projection.exact) {
	    this->cadLastStructuralProjectionHistogram.exact = FALSE;
	} else {
	    this->cadLastStructuralProjectionHistogram.visibleCount =
		view_lod_saturating_add(
		    this->cadLastStructuralProjectionHistogram.visibleCount,
		    projection.visibleCount > static_cast<uint64_t>(SIZE_MAX) ?
			SIZE_MAX : static_cast<size_t>(projection.visibleCount));
	    for (size_t bucket = 0;
		    bucket < CadStructuralProjectionHistogram::BucketCount;
		    ++bucket) {
		const uint64_t count = projection.cumulativeCount[bucket];
		this->cadLastStructuralProjectionHistogram.cumulativeCount[bucket] =
		    view_lod_saturating_add(
			this->cadLastStructuralProjectionHistogram.
			    cumulativeCount[bucket],
			count > static_cast<uint64_t>(SIZE_MAX) ?
			    SIZE_MAX : static_cast<size_t>(count));
	    }
	    this->cadLastStructuralProjectionHistogram.revision ^=
		projection.revision;
	    this->cadLastStructuralProjectionHistogram.revision *=
		1099511628211ULL;
	}
	if (!work.exact)
	    this->cadLastPresentationFrameExact = FALSE;
	const uint64_t presented = work.triangleCount >
	    UINT64_MAX - work.lineCount ? UINT64_MAX :
	    work.triangleCount + work.lineCount;
	if (tier < 0 || tier > 6 || !work.exact || !presented) {
	    presentedValid = FALSE;
	} else {
	    this->cadLastPresentedPrimitiveCount =
		presented > static_cast<uint64_t>(
		    SIZE_MAX - this->cadLastPresentedPrimitiveCount) ?
		SIZE_MAX : this->cadLastPresentedPrimitiveCount +
		    static_cast<size_t>(presented);
	}

	const int presentationDrawMode = assembly->drawMode.getValue();
	BObolLodCounts counts;
	int costDrawMode = BOBOL_LOD_DRAW_SHADED;
	if (presentationDrawMode == SoCADAssembly::WIREFRAME) {
	    counts.lineCount = work.lineCount;
	    costDrawMode = BOBOL_LOD_DRAW_WIRE;
	} else if (presentationDrawMode == SoCADAssembly::SHADED) {
	    counts.faceCount = work.triangleCount;
	} else if (presentationDrawMode == SoCADAssembly::SHADED_WITH_EDGES ||
		presentationDrawMode == SoCADAssembly::HIDDEN_LINE) {
	    counts.faceCount = work.triangleCount;
	    counts.lineCount = work.lineCount;
	    costDrawMode = BOBOL_LOD_DRAW_HIDDEN_LINE;
	} else {
	    presentedRenderCostValid = FALSE;
	}
	counts.pointCount = work.positionCount;
	counts.normalCount = work.normalCount;
	if (!work.exact ||
	    (!counts.faceCount && !counts.lineCount && !counts.pointCount)) {
	    presentedRenderCostValid = FALSE;
	} else {
	    const size_t occurrences = work.occurrenceCount >
		static_cast<uint64_t>(SIZE_MAX) ? SIZE_MAX :
		static_cast<size_t>(work.occurrenceCount);
	    this->cadLastPresentedRenderCost = view_lod_saturating_add(
		this->cadLastPresentedRenderCost,
		bobol_lod_render_cost_units(
		    counts, costDrawMode, occurrences));
	}

	if (tier != 6) {
	    gpuMeasurementValid = FALSE;
	} else {
	    const uint64_t timerSerial = assembly->gpuTimerSampleSerial();
	    const uint64_t timerNanoseconds =
		assembly->lastGpuRenderNanoseconds();
	    if (!timerSerial || !timerNanoseconds) {
		gpuMeasurementValid = FALSE;
	    } else {
		const float timerPointThreshold =
		    assembly->lastGpuPointProxyPixelThreshold();
		if (!std::isfinite(timerPointThreshold)) {
		    gpuMeasurementValid = FALSE;
		} else if (!haveGpuPointProxyPixelThreshold) {
		    this->cadLastGpuPointProxyPixelThreshold =
			timerPointThreshold;
		    haveGpuPointProxyPixelThreshold = TRUE;
		} else if (std::fabs(
			this->cadLastGpuPointProxyPixelThreshold -
			timerPointThreshold) > 0.01f) {
		    /* Multiple assemblies share one point-cut policy.  A mixed
		     * asynchronous sample has no single cut which can update it. */
		    gpuMeasurementValid = FALSE;
		}
		const uint64_t triangles =
		    assembly->lastGpuRenderedTriangleCount();
		this->cadLastGpuFaces =
		    triangles > static_cast<uint64_t>(
			SIZE_MAX - this->cadLastGpuFaces) ?
		    SIZE_MAX : this->cadLastGpuFaces +
			static_cast<size_t>(triangles);
		this->cadLastGpuNanoseconds =
		    timerNanoseconds >
			UINT64_MAX - this->cadLastGpuNanoseconds ?
		    UINT64_MAX : this->cadLastGpuNanoseconds +
			timerNanoseconds;
		this->cadLastGpuSerial ^= static_cast<uint64_t>(
		    reinterpret_cast<uintptr_t>(assembly));
		this->cadLastGpuSerial *= 1099511628211ULL;
		this->cadLastGpuSerial ^= timerSerial;
		this->cadLastGpuSerial *= 1099511628211ULL;
	    }
	}

	/* Flattened tiers retain transformed atlas ranges.  Tier 6 must report
	 * exact prepared replay before its timing is reusable calibration. */
	if (tier != 3 && tier != 4 &&
	    (tier != 6 || !assembly->lastRenderUsedPreparedReplay()))
	    preparedReplay = FALSE;

	const Obol::CadGpuResourceSnapshot snapshot =
	    assembly->gpuResourceSnapshot();
	if (!snapshot.frameSerial)
	    continue;
	CadGpuResourceStatus &status = this->cadGpuResourceStatusValue;
	status.trackedBufferBytes = view_lod_saturating_add(
	    status.trackedBufferBytes, snapshot.trackedBufferBytes);
	status.ordinaryPartBufferBytes = view_lod_saturating_add(
	    status.ordinaryPartBufferBytes,
	    snapshot.ordinaryPartBufferBytes);
	status.progressiveCutBufferBytes = view_lod_saturating_add(
	    status.progressiveCutBufferBytes,
	    snapshot.progressiveCutBufferBytes);
	status.progressiveActiveCutBytes = view_lod_saturating_add(
	    status.progressiveActiveCutBytes,
	    snapshot.progressiveActiveCutBytes);
	status.batchBufferBytes = view_lod_saturating_add(
	    status.batchBufferBytes, snapshot.batchBufferBytes);
	status.triangleAtlasAllocatedBytes = view_lod_saturating_add(
	    status.triangleAtlasAllocatedBytes,
	    snapshot.triangleAtlasAllocatedBytes);
	status.triangleAtlasLiveBytes = view_lod_saturating_add(
	    status.triangleAtlasLiveBytes,
	    snapshot.triangleAtlasLiveBytes);
	status.triangleAtlasConfiguredCapacityBytes = view_lod_saturating_add(
	    status.triangleAtlasConfiguredCapacityBytes,
	    snapshot.triangleAtlasBudgetBytes);
	status.triangleAtlasPartCount = view_lod_saturating_add(
	    status.triangleAtlasPartCount,
	    snapshot.triangleAtlasPartCount);
	status.triangleAtlasPageCount = view_lod_saturating_add(
	    status.triangleAtlasPageCount,
	    snapshot.triangleAtlasPageCount);
	status.ordinaryPartFullUploadBytes = view_lod_saturating_add_u64(
	    status.ordinaryPartFullUploadBytes,
	    snapshot.ordinaryPartFullUploadBytes);
	status.ordinaryPartSuffixUploadBytes = view_lod_saturating_add_u64(
	    status.ordinaryPartSuffixUploadBytes,
	    snapshot.ordinaryPartSuffixUploadBytes);
	status.ordinaryPartGpuCopyBytes = view_lod_saturating_add_u64(
	    status.ordinaryPartGpuCopyBytes,
	    snapshot.ordinaryPartGpuCopyBytes);
	status.ordinaryPartLineageReuseCount = view_lod_saturating_add_u64(
	    status.ordinaryPartLineageReuseCount,
	    snapshot.ordinaryPartLineageReuseCount);
	status.ordinaryPartLineageReplacementCount =
	    view_lod_saturating_add_u64(
		status.ordinaryPartLineageReplacementCount,
		snapshot.ordinaryPartLineageReplacementCount);
	status.triangleAtlasFullUploadBytes = view_lod_saturating_add_u64(
	    status.triangleAtlasFullUploadBytes,
	    snapshot.triangleAtlasFullUploadBytes);
	status.triangleAtlasSuffixUploadBytes = view_lod_saturating_add_u64(
	    status.triangleAtlasSuffixUploadBytes,
	    snapshot.triangleAtlasSuffixUploadBytes);
	status.triangleAtlasLineageReuseCount = view_lod_saturating_add_u64(
	    status.triangleAtlasLineageReuseCount,
	    snapshot.triangleAtlasLineageReuseCount);
	status.pressureProxyCount = view_lod_saturating_add(
	    status.pressureProxyCount, snapshot.pressureProxyCount);
	status.progressiveEvictionCount = view_lod_saturating_add_u64(
	    status.progressiveEvictionCount,
	    snapshot.progressiveEvictionCount);
	status.triangleAtlasReclamationCount =
	    view_lod_saturating_add_u64(
		status.triangleAtlasReclamationCount,
		snapshot.triangleAtlasReclamationCount);
	if (snapshot.atlasAdmissionPressure)
	    status.memoryPressure = TRUE;

	/* Hash both identity and local completed-frame serial.  Merely summing
	 * per-context counters can alias when one context advances while another
	 * remains unchanged. */
	resourceSerial ^= static_cast<uint64_t>(
	    reinterpret_cast<uintptr_t>(assembly));
	resourceSerial *= 1099511628211ULL;
	resourceSerial ^= snapshot.frameSerial;
	resourceSerial *= 1099511628211ULL;
	resourceSampled = TRUE;
    }

    if (!haveAssembly || !this->cadLastPresentationFrameExact)
	this->cadLastStructuralProjectionHistogram.exact = FALSE;

    this->cadLastPresentedPrimitiveCountValid =
	haveAssembly && presentedValid;
    this->cadLastPresentedRenderCostValid =
	haveAssembly && presentedRenderCostValid &&
	this->cadLastPresentedRenderCost > 0;
    this->cadLastGpuMeasurementValid =
	haveAssembly && gpuMeasurementValid &&
	haveGpuPointProxyPixelThreshold;
    if (!this->cadLastGpuMeasurementValid)
	this->cadLastGpuSerial = 0;
    else if (!this->cadLastGpuSerial)
	this->cadLastGpuSerial = 1;
    this->cadLastPreparedReplay = haveAssembly && preparedReplay;
    this->cadGpuResourceStatusValid = resourceSampled;
    /* qged diagnostics are JSON and therefore carry this token through a
     * signed 64-bit integer.  It is an opaque change detector, not an
     * arithmetic counter, so reserve the sign bit and keep zero as "none". */
    if (resourceSampled) {
	this->cadGpuResourceStatusValue.sampleSerial =
	    resourceSerial & INT64_MAX;
	if (!this->cadGpuResourceStatusValue.sampleSerial)
	    this->cadGpuResourceStatusValue.sampleSerial = 1;
    }
    this->cadPresentationFrameStatusValid = TRUE;
}

SbBool
BObolViewLodState::cadGpuResourceStatus(
	CadGpuResourceStatus &status) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    status = this->cadGpuResourceStatusValue;
    return this->cadGpuResourceStatusValid;
}

SbBool
BObolViewLodState::lastCadPresentationUsedPreparedReplay(void) const
{
    if (!this->cadPresentationFrameStatusValid)
	this->refreshCadPresentationFrameStatus();
    return this->cadLastPreparedReplay;
}

int
BObolViewLodState::maximumActiveProgressiveCut(void) const
{
    for (int cut =
	    static_cast<int>(sizeof(this->cadProgressiveCutCounts) /
		sizeof(this->cadProgressiveCutCounts[0])) - 1;
	cut >= 0; --cut)
	if (this->cadProgressiveCutCounts[cut])
	    return cut;
    return -1;
}

size_t
BObolViewLodState::cadRenderCostAtProgressiveCutCeiling(int cut) const
{
    if (cut < 0)
	return this->cadActiveRenderCost;
    const int boundedCut = std::min<int>(
	BOBOL_MESH_LOD_CUT_COUNT_MAX - 1, cut);
    return this->cadProgressiveCeilingRenderCosts[boundedCut];
}

int
BObolViewLodState::cadProgressiveCutWithinRenderCost(
    size_t renderCost, int maximumCut) const
{
    if (!renderCost || !this->cadMeshPayloadCountValue)
	return -1;

    int boundedMaximum = maximumCut < 0 ?
	this->maximumActiveProgressiveCut() :
	std::min<int>(BOBOL_MESH_LOD_CUT_COUNT_MAX - 1, maximumCut);
    if (boundedMaximum < 0)
	return -1;
    for (int cut = boundedMaximum; cut >= 0; --cut)
	if (this->cadProgressiveCeilingRenderCosts[cut] <= renderCost)
	    return cut;
    return -1;
}

int
BObolViewLodState::singleCadProgressiveCutWithinRenderCost(
	size_t renderCost) const
{
    if (!renderCost || this->cadMeshPayloadCountValue != 1)
	return -1;

    const CadPayload *payload = NULL;
    for (const auto &source : this->cadSourceBindings) {
	for (const auto &occurrence : source.second) {
	    const CadPayload *candidate = occurrence.second.get();
	    if (!view_lod_cad_payload_is_mesh(candidate))
		continue;
	    if (payload && payload != candidate)
		return -1;
	    payload = candidate;
	}
    }
    if (!payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() || payload->activeCut < 0)
	return -1;

    const int minimumCut = payload->progressiveMesh->minimumCut();
    const int maximumCut = std::min(payload->activeCut,
	payload->progressiveMesh->maximumCut());
    int result = -1;
    for (int cut = minimumCut; cut <= maximumCut; ++cut) {
	if (!view_lod_progressive_can_draw(payload->progressiveMesh,
		payload->requiredChunks, cut))
	    continue;
	const BObolLodCounts counts = view_lod_progressive_counts(
	    payload->progressiveMesh, payload->requiredChunks, cut,
	    payload->hasNormals);
	const size_t cost = bobol_lod_render_cost_units(
	    counts, payload->drawMode, 1);
	if (cost <= renderCost)
	    result = cut;
    }
    return result;
}

void
BObolViewLodState::setCadPresentationProgressiveCutCeiling(int cut) const
{
    cut = cut < 0 ? -1 : std::min<int>(
	Obol::ProgressiveCutLimit - 1, cut);
    this->cadPresentationProgressiveLodCeiling = cut;
    for (const auto &presentation : this->cadPresentations)
	if (presentation.second.assembly)
	    presentation.second.assembly->progressiveCutCeiling = cut;
}

void
BObolViewLodState::setCadPresentationPointProxyPixelThreshold(
    float pixels) const
{
    pixels = std::isfinite(pixels) ?
	std::max(1.0f, std::min(64.0f, pixels)) : 1.0f;
    const float oldEffective = std::max(
	this->cadPresentationPointProxyPixelThreshold,
	this->cadPresentationDiscoveryPointProxyPixelThreshold);
    this->cadPresentationPointProxyPixelThreshold = pixels;
    const float effective = std::max(pixels,
	this->cadPresentationDiscoveryPointProxyPixelThreshold);
    if (std::fabs(effective - oldEffective) <= 1.0e-6f)
	return;
    for (const auto &presentation : this->cadPresentations)
	if (presentation.second.assembly)
	    presentation.second.assembly->pointProxyPixelThreshold = effective;
}

void
BObolViewLodState::setCadPresentationDiscoveryPointProxyPixelThreshold(
    float pixels) const
{
    pixels = std::isfinite(pixels) ?
	std::max(1.0f, std::min(64.0f, pixels)) : 1.0f;
    const float oldEffective = std::max(
	this->cadPresentationPointProxyPixelThreshold,
	this->cadPresentationDiscoveryPointProxyPixelThreshold);
    this->cadPresentationDiscoveryPointProxyPixelThreshold = pixels;
    const float effective = std::max(
	this->cadPresentationPointProxyPixelThreshold, pixels);
    if (std::fabs(effective - oldEffective) <= 1.0e-6f)
	return;
    for (const auto &presentation : this->cadPresentations)
	if (presentation.second.assembly)
	    presentation.second.assembly->pointProxyPixelThreshold = effective;
}

void
BObolViewLodState::setCadPresentationCameraMotionFrameReuse(
    SbBool enabled) const
{
    enabled = enabled ? TRUE : FALSE;
    if (this->cadPresentationCameraMotionFrameReuse == enabled)
	return;
    this->cadPresentationCameraMotionFrameReuse = enabled;
    for (const auto &presentation : this->cadPresentations)
	if (presentation.second.assembly &&
	    presentation.second.assembly->cameraMotionFrameReuse.getValue() !=
		enabled)
	    presentation.second.assembly->cameraMotionFrameReuse = enabled;
}

size_t
BObolViewLodState::estimateDisplayMeshBytes(void) const
{
    std::vector<MeshPayloadPtr> payloads =
	view_lod_unique_payloads(this->meshBindings);
    size_t bytes = 0;
    for (size_t i = 0; i < payloads.size(); i++)
	bytes += payloads[i]->estimateBytes();
    std::vector<ProxyPayloadPtr> proxyPayloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    for (size_t i = 0; i < proxyPayloads.size(); i++)
	bytes += proxyPayloads[i]->estimateBytes();
    return view_lod_saturating_add(bytes, this->cadDisplayMeshBytes);
}

void
BObolViewLodState::residentMeshDemands(
    std::vector<BObolLodResidentDemand> &demands) const
{
    demands = this->cadResidentDemands;
    if (this->meshBindings.empty())
	return;

    struct DemandAggregate {
	int cut = -1;
	unsigned int channelMask = 0;
	std::vector<uint32_t> chunkIds;
    };
    std::unordered_map<std::string, DemandAggregate> richestCuts;
    richestCuts.reserve(
	demands.size() + this->meshBindings.size());
    for (const BObolLodResidentDemand &demand : demands) {
	DemandAggregate &aggregate =
	    richestCuts[demand.assetKey.getString()];
	aggregate.cut = std::max(aggregate.cut, demand.cut);
	aggregate.channelMask |= demand.channelMask;
	std::vector<uint32_t> merged;
	merged.reserve(aggregate.chunkIds.size() + demand.chunkIds.size());
	std::set_union(aggregate.chunkIds.begin(), aggregate.chunkIds.end(),
	    demand.chunkIds.begin(), demand.chunkIds.end(),
	    std::back_inserter(merged));
	aggregate.chunkIds.swap(merged);
    }

    std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    for (const MeshPayloadPtr &payload : meshPayloads) {
	if (!payload || !payload->progressiveMesh ||
	    !payload->progressiveMesh->isValid() ||
	    payload->cacheKey.getLength() == 0)
	    continue;
	const int demandCut = view_lod_resident_demand_cut(
	    payload->activeCut, payload->requestedCut,
	    payload->progressiveMesh->residentCut());
	if (demandCut < 0)
	    continue;
	DemandAggregate &aggregate =
	    richestCuts[payload->cacheKey.getString()];
	aggregate.cut = std::max(aggregate.cut, demandCut);
	std::vector<uint32_t> merged;
	merged.reserve(aggregate.chunkIds.size() +
	    payload->requiredChunks.size());
	std::set_union(aggregate.chunkIds.begin(), aggregate.chunkIds.end(),
	    payload->requiredChunks.begin(), payload->requiredChunks.end(),
	    std::back_inserter(merged));
	aggregate.chunkIds.swap(merged);
    }

    demands.clear();
    demands.reserve(richestCuts.size());
    for (const auto &entry : richestCuts) {
	BObolLodResidentDemand demand;
	demand.assetKey = entry.first.c_str();
	demand.cut = entry.second.cut;
	demand.channelMask = entry.second.channelMask;
	demand.chunkIds = entry.second.chunkIds;
	demands.push_back(demand);
    }
}

uint64_t
BObolViewLodState::residentMeshDemandRevision(void) const
{
    return this->cadResidentDemandRevision;
}

size_t
BObolViewLodState::applyResidentMeshCompaction(
    const BObolLodResidentCompaction &result)
{
    if (result.assetKey.getLength() == 0 || !result.progressiveMesh ||
	!result.progressiveMesh->isValid() || result.residentCut < 0)
	return 0;
    /* Another view may have extended the shared asset after this completion
     * was queued.  Never publish an older prepared generation over that
     * richer replacement; its next stable snapshot will compact again if
     * the aggregate demand permits it. */
    if (result.progressiveMesh->residentCut() != result.residentCut ||
	(result.preparedCadGeometryRevision &&
	 result.preparedCadGeometryRevision !=
	    result.progressiveMesh->revision()))
	return 0;
    const auto indexed = this->cadPayloadsByAssetKey.find(
	result.assetKey.getString());
    if (indexed == this->cadPayloadsByAssetKey.end())
	return 0;

    /*
     * Every occurrence retains the same progressive asset, but each source
     * has its own SoCADAssembly and part namespace.  Update the cheap payload
     * handles for all occurrences, then journal one representative per source
     * so each assembly replaces the shared part exactly once.
     */
    std::unordered_set<std::string> publishedSources;
    size_t published = 0;
    for (CadPayload *payload : indexed->second) {
	if (!payload || payload->progressiveMesh != result.progressiveMesh)
	    continue;
	/*
	 * The resident-demand table is an incremental index over mutable payload
	 * state.  The worker has already changed the shared mesh revision/resident
	 * cut by the time this completion reaches the presentation owner; the
	 * payload's residentCut is therefore the authoritative installed-generation
	 * snapshot needed to remove the old contribution.  Remove it before
	 * adopting the result, otherwise a retired high-water cut can survive in
	 * cadResidentDemands and schedule another, incoherent trim after a zoom.
	 */
	const bool wasSatisfied =
	    view_lod_cad_payload_is_satisfied(payload);
	this->removeCadResidentDemand(payload);
	payload->residentCut = result.residentCut;
	const unsigned int needed =
	    view_lod_cad_payload_channel_mask(payload);
	if (result.preparedCadGeometry &&
	    (result.channelMask & needed) == needed) {
	    payload->preparedCadGeometry = result.preparedCadGeometry;
	    payload->preparedCadGeometryRevision =
		result.preparedCadGeometryRevision;
	    if (view_lod_progressive_can_draw(
		    payload->progressiveMesh, payload->requiredChunks,
		    payload->activeCut))
		payload->presentedChunks = payload->requiredChunks;
	} else {
	    /* Release the retired generation even when a non-authored normal
	     * mode must prepare its replacement later. */
	    payload->preparedCadGeometry.reset();
	    payload->preparedCadGeometryRevision = 0;
	}
	const bool isSatisfied =
	    view_lod_cad_payload_is_satisfied(payload);
	if (wasSatisfied != isSatisfied) {
	    if (isSatisfied)
		this->cadSatisfiedMeshPayloadCount = view_lod_saturating_add(
		    this->cadSatisfiedMeshPayloadCount, 1);
	    else
		this->cadSatisfiedMeshPayloadCount =
		    view_lod_saturating_subtract(
			this->cadSatisfiedMeshPayloadCount, 1);
	}
	if (payload->sourceInstanceKey.getLength() > 0) {
	    const std::string sourceKey =
		payload->sourceBindingKey.getString();
	    const std::string occurrenceKey =
		payload->sourceInstanceKey.getString();
	    if (isSatisfied) {
		const auto unsatisfied =
		    this->cadUnsatisfiedOccurrencesBySource.find(sourceKey);
		if (unsatisfied !=
			this->cadUnsatisfiedOccurrencesBySource.end()) {
		    unsatisfied->second.erase(occurrenceKey);
		    if (unsatisfied->second.empty())
			this->cadUnsatisfiedOccurrencesBySource.erase(
			    unsatisfied);
		}
	    } else {
		this->cadUnsatisfiedOccurrencesBySource[sourceKey].insert(
		    occurrenceKey);
	    }
	}
	this->addCadResidentDemand(payload);
	const char *traceFilter = getenv("BOBOL_LOD_TRACE_OBJECT");
	if (getenv("BOBOL_LOD_TRACE_COMPACTION") && traceFilter &&
	    traceFilter[0] &&
	    (strstr(payload->sourceName.getString(), traceFilter) ||
	     strstr(payload->sourcePath.getString(), traceFilter))) {
	    const SbVec3f qmin = payload->progressiveMesh->quantizationMinimum();
	    const SbVec3f qmax = payload->progressiveMesh->quantizationMaximum();
	    BObolMeshLodCutInfo cutInfo = {};
	    (void)payload->progressiveMesh->cutInfo(
		payload->activeCut, &cutInfo);
	    bu_log("BObol resident compaction adopt object=%s active=%d "
		   "requested=%d pixels=%.9g target=%.9g error=%.9g "
		   "object_error=%.17g q=[%.9g %.9g %.9g]-"
		   "[%.9g %.9g %.9g] revision=%llu chunks=%zu/%zu\n",
		   payload->sourceName.getString(), payload->activeCut,
		   payload->requestedCut, payload->projectedPixelDiameter,
		   payload->targetPixelError,
		   payload->progressiveMesh->projectedErrorAtCut(
		       payload->activeCut, payload->projectedPixelDiameter),
		   cutInfo.object_error,
		   qmin[0], qmin[1], qmin[2], qmax[0], qmax[1], qmax[2],
		   static_cast<unsigned long long>(
		       payload->progressiveMesh->revision()),
		   payload->requiredChunks.size(),
		   payload->presentedChunks.size());
	}

	const std::string sourceKey = payload->sourceBindingKey.getString();
	if (!publishedSources.insert(sourceKey).second)
	    continue;
	if (payload->sourceInstanceKey.getLength() > 0) {
	    this->noteCadOccurrenceChanged(
		sourceKey, payload->sourceInstanceKey);
	} else {
	    this->noteResidentMeshesChanged(
		"legacy-resident-prefix-compaction");
	}
	published++;
    }
    return published;
}

uint64_t
BObolViewLodState::cadRevision(void) const
{
    return this->cadBindingsRevision;
}

void
BObolViewLodState::cadOccurrenceChangesSince(
    const SoBRLDatabaseSource *source, uint64_t revision,
    std::vector<SbString> &occurrenceKeys, SbBool &fullResync) const
{
    occurrenceKeys.clear();
    fullResync = revision < this->cadFullResyncRevision ? TRUE : FALSE;
    if (fullResync || !source)
	return;

    const std::string sourceKey = view_lod_source_primary_key(source);
    const auto changes = this->cadOccurrenceChanges.find(sourceKey);
    if (changes == this->cadOccurrenceChanges.end())
	return;

    std::unordered_set<std::string> seen;
    for (const CadOccurrenceChange &change : changes->second) {
	if (change.revision <= revision ||
	    change.occurrenceKey.getLength() == 0)
	    continue;
	if (seen.insert(change.occurrenceKey.getString()).second)
	    occurrenceKeys.push_back(change.occurrenceKey);
    }
}

void
BObolViewLodState::acknowledgeCadOccurrenceChanges(
    const SoBRLDatabaseSource *source, uint64_t revision) const
{
    if (!source)
	return;
    const std::string sourceKey = view_lod_source_primary_key(source);
    auto changes = this->cadOccurrenceChanges.find(sourceKey);
    if (changes == this->cadOccurrenceChanges.end())
	return;
    std::vector<CadOccurrenceChange> &entries = changes->second;
    entries.erase(std::remove_if(entries.begin(), entries.end(),
	[revision](const CadOccurrenceChange &change) {
	    return change.revision <= revision;
	}), entries.end());
    if (entries.empty())
	this->cadOccurrenceChanges.erase(changes);
}

void
BObolViewLodState::noteCadOccurrenceChanged(
    const std::string &sourceBindingKey, const SbString &occurrenceKey)
{
    if (sourceBindingKey.empty() || occurrenceKey.getLength() == 0) {
	this->noteResidentMeshesChanged("invalid-occurrence-change");
	return;
    }
    this->cadBindingsRevision++;
    if (!this->cadBindingsRevision) {
	this->cadBindingsRevision = 1;
	this->cadFullResyncRevision = 1;
	this->cadOccurrenceChanges.clear();
	return;
    }
    CadOccurrenceChange change;
    change.revision = this->cadBindingsRevision;
    change.occurrenceKey = occurrenceKey;
    this->cadOccurrenceChanges[sourceBindingKey].push_back(
	std::move(change));
}

void
BObolViewLodState::noteResidentMeshesChanged(const char *reason)
{
    this->cadBindingsRevision++;
    if (!this->cadBindingsRevision)
	this->cadBindingsRevision++;
    this->cadFullResyncRevision = this->cadBindingsRevision;
    this->cadOccurrenceChanges.clear();
    if (getenv("BOBOL_COMPACT_PRESENTATION_TIMING"))
	bu_log("BObol compact presentation full-resync reason=%s "
	       "revision=%llu payloads=%zu meshes=%zu\n",
	       reason ? reason : "unknown",
	       static_cast<unsigned long long>(this->cadBindingsRevision),
	       this->cadValidPayloadCount, this->cadMeshPayloadCountValue);
}

size_t
BObolViewLodState::evictDisplayMeshes(unsigned int *evictedMeshCount)
{
    std::vector<MeshPayloadPtr> payloads =
	view_lod_unique_payloads(this->meshBindings);
    size_t bytes = 0;
    for (size_t i = 0; i < payloads.size(); i++)
	bytes += payloads[i]->estimateBytes();
    std::vector<ProxyPayloadPtr> proxyPayloads =
	view_lod_unique_proxy_payloads(this->proxyBindings);
    for (size_t i = 0; i < proxyPayloads.size(); i++)
	bytes += proxyPayloads[i]->estimateBytes();
    size_t cadPayloadCount = 0;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    if (!binding.second || !binding.second->isValid())
		continue;
	    bytes += binding.second->estimateBytes();
	    cadPayloadCount++;
	}

    if (evictedMeshCount)
	*evictedMeshCount = static_cast<unsigned int>(
	    payloads.size() + proxyPayloads.size() + cadPayloadCount);
    this->meshBindings.clear();
    this->proxyBindings.clear();
    this->cadSourceBindings.clear();
    this->cadSourceEntryBindings.clear();
    this->cadAssetBindings.clear();
    this->cadBindings.clear();
    this->clearCadPayloadMetrics();
    this->noteResidentMeshesChanged("display-evict-all");
    return bytes;
}

size_t
BObolViewLodState::evictDisplayMeshPayloads(
    unsigned int *evictedMeshCount)
{
    std::vector<MeshPayloadPtr> meshPayloads =
	view_lod_unique_payloads(this->meshBindings);
    std::vector<CadPayloadPtr> cadMeshPayloads;
    for (const auto &source : this->cadSourceBindings)
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (payload && payload->isValid() &&
		(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
		 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL))
		cadMeshPayloads.push_back(payload);
	}

    size_t bytes = 0;
    for (const MeshPayloadPtr &payload : meshPayloads)
	bytes += payload->estimateBytes();
    for (const CadPayloadPtr &payload : cadMeshPayloads)
	bytes += payload->estimateBytes();
    for (const CadPayloadPtr &payload : cadMeshPayloads)
	this->removeCadPayloadMetrics(payload.get());

    if (evictedMeshCount)
	*evictedMeshCount = static_cast<unsigned int>(
	    meshPayloads.size() + cadMeshPayloads.size());
    this->meshBindings.clear();
    this->cadAssetBindings.clear();
    this->cadSourceEntryBindings.clear();
    for (auto source = this->cadSourceBindings.begin();
	 source != this->cadSourceBindings.end();) {
	for (auto binding = source->second.begin();
	     binding != source->second.end();) {
	    const CadPayloadPtr &payload = binding->second;
	    if (payload &&
		(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
		 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL))
		binding = source->second.erase(binding);
	    else
		++binding;
	}
	if (source->second.empty())
	    source = this->cadSourceBindings.erase(source);
	else
	    ++source;
    }
    for (const auto &source : this->cadSourceBindings) {
	std::vector<CadPayload *> *indexed = NULL;
	for (const auto &binding : source.second) {
	    const CadPayloadPtr &payload = binding.second;
	    if (!payload || !payload->sourceRoutingId ||
		payload->sourceEntryIndex == UINT32_MAX)
		continue;
	    if (!indexed)
		indexed = &this->cadSourceEntryBindings[
		    payload->sourceRoutingId];
	    const size_t entryIndex = static_cast<size_t>(
		payload->sourceEntryIndex);
	    if (indexed->size() <= entryIndex)
		indexed->resize(entryIndex + 1, NULL);
	    (*indexed)[entryIndex] = payload.get();
	}
    }
    for (auto binding = this->cadBindings.begin();
	 binding != this->cadBindings.end();) {
	const CadPayloadPtr &payload = binding->second;
	if (payload &&
	    (payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	     payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL)) {
	    binding = this->cadBindings.erase(binding);
	    continue;
	}
	++binding;
    }
    if (!cadMeshPayloads.empty())
	this->noteResidentMeshesChanged("display-mesh-evict");

    return bytes;
}

SoBRLViewLodElement::~SoBRLViewLodElement(void)
{
}

void
SoBRLViewLodElement::initClass(void)
{
    SO_ELEMENT_INIT_CLASS(SoBRLViewLodElement, SoElement);
}

void
SoBRLViewLodElement::init(SoState *state)
{
    inherited::init(state);
    this->viewState = NULL;
}

void
SoBRLViewLodElement::push(SoState *state)
{
    const SoBRLViewLodElement *prev =
	static_cast<const SoBRLViewLodElement *>(
	    this->getNextInStack());
    this->viewState = prev ? prev->viewState : NULL;
    inherited::push(state);
}

SbBool
SoBRLViewLodElement::matches(const SoElement *element) const
{
    const SoBRLViewLodElement *other =
	static_cast<const SoBRLViewLodElement *>(element);
    return other && this->viewState == other->viewState ? TRUE : FALSE;
}

SoElement *
SoBRLViewLodElement::copyMatchInfo(void) const
{
    SoBRLViewLodElement *element = new SoBRLViewLodElement;
    element->viewState = this->viewState;
    return element;
}

void
SoBRLViewLodElement::set(SoState *state,
			 SoNode *UNUSED(node),
			 const BObolViewLodState *newViewState)
{
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return;

    SoBRLViewLodElement *element =
	static_cast<SoBRLViewLodElement *>(
	    SoElement::getElement(state,
				  SoBRLViewLodElement::getClassStackIndex()));
    element->viewState = newViewState;
}

const BObolViewLodState *
SoBRLViewLodElement::get(SoState *state)
{
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return NULL;

    const SoBRLViewLodElement *element =
	static_cast<const SoBRLViewLodElement *>(
	    SoElement::getConstElement(state,
				       SoBRLViewLodElement::getClassStackIndex()));
    return element ? element->viewState : NULL;
}

SoBRLViewLodGroup::SoBRLViewLodGroup(void) :
    viewState(NULL),
    softwareWireMode(SoCADViewState::SOFTWARE_WIRE_AUTO)
{
    SO_NODE_CONSTRUCTOR(SoBRLViewLodGroup);
}

SoBRLViewLodGroup::~SoBRLViewLodGroup(void)
{
}

void
SoBRLViewLodGroup::initClass(void)
{
    SO_NODE_INIT_CLASS(SoBRLViewLodGroup, SoGroup, "Group");
    SO_ENABLE(SoGLRenderAction, SoBRLViewLodElement);
    SO_ENABLE(SoGetBoundingBoxAction, SoBRLViewLodElement);
    SO_ENABLE(SoRayPickAction, SoBRLViewLodElement);
}

void
SoBRLViewLodGroup::setViewLodState(BObolViewLodState *newViewState)
{
    this->viewState = newViewState;
}

BObolViewLodState *
SoBRLViewLodGroup::getViewLodState(void) const
{
    return this->viewState;
}

void
SoBRLViewLodGroup::setSoftwareWireMode(int mode)
{
    if (mode < SoCADViewState::SOFTWARE_WIRE_AUTO ||
	mode > SoCADViewState::SOFTWARE_WIRE_FAST)
	mode = SoCADViewState::SOFTWARE_WIRE_AUTO;
    this->softwareWireMode = mode;
}

int
SoBRLViewLodGroup::getSoftwareWireMode(void) const
{
    return this->softwareWireMode;
}

SbBool
SoBRLViewLodGroup::pushViewState(SoAction *action)
{
    if (!action || !this->viewState)
	return FALSE;

    SoState *state = action->getState();
    if (!state || !state->isElementEnabled(SoBRLViewLodElement::getClassStackIndex()))
	return FALSE;

    state->push();
    SoBRLViewLodElement::set(state, this, this->viewState);
    if (state->isElementEnabled(
	    SoCADViewStateElement::getClassStackIndex())) {
	Obol::CadViewState cadState = SoCADViewStateElement::get(state);
	cadState.softwareWireMode =
	    static_cast<Obol::CadSoftwareWireMode>(this->softwareWireMode);
	SoCADViewStateElement::set(state, cadState);
    }
    return TRUE;
}

void
SoBRLViewLodGroup::popViewState(SoAction *action, SbBool pushed)
{
    if (!pushed || !action)
	return;

    SoState *state = action->getState();
    if (state)
	state->pop();
}

void
SoBRLViewLodGroup::doAction(SoAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::doAction(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::GLRender(SoGLRenderAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::GLRender(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::callback(SoCallbackAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::callback(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::getBoundingBox(SoGetBoundingBoxAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::getBoundingBox(action);
    this->popViewState(action, pushed);
}

void
SoBRLViewLodGroup::pick(SoPickAction *action)
{
    const SbBool pushed = this->pushViewState(action);
    inherited::pick(action);
    this->popViewState(action, pushed);
}

const BObolViewLodState::MeshPayload *
bobol_view_lod_mesh_for_action(SoAction *action,
				 const SoBRLMeshShape *shape)
{
    if (!action || !shape)
	return NULL;

    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findMesh(shape) : NULL;
}

const BObolViewLodState *
bobol_view_lod_state_for_action(SoAction *action)
{
    return action ? SoBRLViewLodElement::get(action->getState()) : NULL;
}

const BObolViewLodState::ProxyPayload *
bobol_view_lod_proxy_for_action(SoAction *action,
				  const SoBRLMeshShape *shape)
{
    if (!action || !shape)
	return NULL;

    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findProxy(shape) : NULL;
}

const BObolViewLodState::CadPayload *
bobol_view_lod_cad_for_action(SoAction *action,
				const SoBRLDatabaseSource *source)
{
    if (!action || !source)
	return NULL;

    const BObolViewLodState *viewState =
	SoBRLViewLodElement::get(action->getState());
    return viewState ? viewState->findCad(source) : NULL;
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
