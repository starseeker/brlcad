/*        M E S H _ L O D _ S U B M I T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/str.h"

#include "BObol/BDatabaseSource.h"
#include "BObol/BLodService.h"
#include "BObol/BMeshLodCache.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BMeshShape.h"
#include "BObol/BViewLod.h"

#include "raytrace.h"

#include <Inventor/SbString.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoNode.h>

#include <algorithm>
#include <cfloat>
#include <cmath>
#include <errno.h>
#include <stdlib.h>
#include <string.h>
#include <unordered_map>

SO_ACTION_SOURCE(SoBRLMeshLodSubmitAction);

SoBRLMeshLodSubmitAction::SoBRLMeshLodSubmitAction(void) :
    service(NULL),
    dbip(NULL),
    databaseId(""),
    databaseRevision(0),
    viewVolume(),
    useViewVolume(FALSE),
    targetPixelError(1.0f),
    generation(0),
    viewRevision(0),
    policyRevision(0),
    providerId("bobol_mesh_lod"),
    providerVersion(BOBOL_MESH_LOD_PROVIDER_VERSION),
    qualityTier(BOBOL_LOD_QUALITY_FAST_DISPLAY),
    refreshMissing(TRUE),
    reset(0),
    useForcedLevel(FALSE),
    forcedLevel(0),
    requireLodBacked(TRUE),
    allowLevelDowngrade(FALSE),
    viewState(NULL),
    compactEntryFirst(0),
    compactEntryLimit(SIZE_MAX),
    compactEntryNext(0),
    compactEntryTotal(0),
    deferredCompactEntries(FALSE),
    visitedMeshCount(0),
    submittedTaskCount(0),
    updatedCutCount(0),
    skippedMeshCount(0),
    diagnostics("")
{
    SO_ACTION_CONSTRUCTOR(SoBRLMeshLodSubmitAction);
    bv_view_info_init(&this->view);
}

SoBRLMeshLodSubmitAction::~SoBRLMeshLodSubmitAction(void)
{
}

void
SoBRLMeshLodSubmitAction::initClass(void)
{
    SO_ACTION_INIT_CLASS(SoBRLMeshLodSubmitAction, SoAction);
    SO_ACTION_ADD_METHOD(SoNode, SoBRLMeshLodSubmitAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoGroup, SoBRLMeshLodSubmitAction::nodeAction);
    SO_ACTION_ADD_METHOD(SoBRLDatabaseSource,
			 SoBRLMeshLodSubmitAction::databaseSourceAction);
    SO_ACTION_ADD_METHOD(SoBRLMeshShape,
			 SoBRLMeshLodSubmitAction::meshShapeAction);
}

void
SoBRLMeshLodSubmitAction::setService(BObolLodService *newService)
{
    this->service = newService;
}

BObolLodService *
SoBRLMeshLodSubmitAction::getService(void) const
{
    return this->service;
}

void
SoBRLMeshLodSubmitAction::setDatabase(struct db_i *newDbip,
				      const char *newDatabaseId, uint64_t newDatabaseRevision)
{
    this->dbip = newDbip;
    this->databaseId = newDatabaseId ? newDatabaseId : "";
    this->databaseRevision = newDatabaseRevision;
}

void
SoBRLMeshLodSubmitAction::setViewInfo(const struct bv_view_info *info)
{
    if (info)
	this->view = *info;
    else
	bv_view_info_init(&this->view);
    bv_view_info_sanitize(&this->view);
}

const struct bv_view_info &
SoBRLMeshLodSubmitAction::getViewInfo(void) const {
    return this->view;
}

void
SoBRLMeshLodSubmitAction::setViewVolume(const SbViewVolume *volume,
	float newTargetPixelError)
{
    if (volume) {
	this->viewVolume = *volume;
	this->useViewVolume = TRUE;
    } else {
	this->useViewVolume = FALSE;
    }
    this->targetPixelError =
	(std::isfinite(newTargetPixelError) && newTargetPixelError > 0.0f) ?
	newTargetPixelError : 1.0f;
}

void
SoBRLMeshLodSubmitAction::setGeneration(uint64_t newGeneration)
{
    this->generation = newGeneration;
}

uint64_t
SoBRLMeshLodSubmitAction::getGeneration(void) const
{
    return this->generation;
}

void
SoBRLMeshLodSubmitAction::setRevisions(uint64_t newViewRevision,
				       uint64_t newPolicyRevision)
{
    this->viewRevision = newViewRevision;
    this->policyRevision = newPolicyRevision;
}

void
SoBRLMeshLodSubmitAction::setProvider(const char *newProviderId,
				      const char *newProviderVersion)
{
    this->providerId = newProviderId ? newProviderId : "";
    this->providerVersion = newProviderVersion ? newProviderVersion : "";
}

void
SoBRLMeshLodSubmitAction::setQualityTier(int newQualityTier)
{
    this->qualityTier = newQualityTier;
}

void
SoBRLMeshLodSubmitAction::setRefreshMissing(SbBool newRefreshMissing)
{
    this->refreshMissing = newRefreshMissing ? TRUE : FALSE;
}

void
SoBRLMeshLodSubmitAction::setReset(int newReset)
{
    this->reset = newReset;
}

void
SoBRLMeshLodSubmitAction::setForcedLevel(int level)
{
    this->useForcedLevel = TRUE;
    this->forcedLevel = level < 0 ? 0 : level;
}

void
SoBRLMeshLodSubmitAction::clearForcedLevel(void)
{
    this->useForcedLevel = FALSE;
    this->forcedLevel = 0;
}

SbBool
SoBRLMeshLodSubmitAction::hasForcedLevel(void) const
{
    return this->useForcedLevel;
}

int
SoBRLMeshLodSubmitAction::getForcedLevel(void) const
{
    return this->forcedLevel;
}

void
SoBRLMeshLodSubmitAction::setRequireLodBacked(SbBool newRequireLodBacked)
{
    this->requireLodBacked = newRequireLodBacked ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getRequireLodBacked(void) const
{
    return this->requireLodBacked;
}

void
SoBRLMeshLodSubmitAction::setAllowLevelDowngrade(SbBool allow)
{
    this->allowLevelDowngrade = allow ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getAllowLevelDowngrade(void) const
{
    return this->allowLevelDowngrade;
}

void
SoBRLMeshLodSubmitAction::setViewLodState(
    BObolViewLodState *newViewState)
{
    this->viewState = newViewState;
}

const BObolViewLodState *
SoBRLMeshLodSubmitAction::getViewLodState(void) const
{
    return this->viewState;
}

void
SoBRLMeshLodSubmitAction::setCompactEntryRange(size_t first, size_t count)
{
    this->compactEntryFirst = first;
    this->compactEntryLimit = count;
}

size_t
SoBRLMeshLodSubmitAction::getCompactEntryNext(void) const
{
    return this->compactEntryNext;
}

size_t
SoBRLMeshLodSubmitAction::getCompactEntryTotal(void) const
{
    return this->compactEntryTotal;
}

SbBool
SoBRLMeshLodSubmitAction::hasDeferredCompactEntries(void) const
{
    return this->deferredCompactEntries;
}

unsigned int
SoBRLMeshLodSubmitAction::getVisitedMeshCount(void) const
{
    return this->visitedMeshCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getSubmittedTaskCount(void) const
{
    return this->submittedTaskCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getUpdatedCutCount(void) const
{
    return this->updatedCutCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getSkippedMeshCount(void) const
{
    return this->skippedMeshCount;
}

const SbString &
SoBRLMeshLodSubmitAction::getDiagnostics(void) const
{
    return this->diagnostics;
}

void
SoBRLMeshLodSubmitAction::beginTraversal(SoNode *node)
{
    this->visitedMeshCount = 0;
    this->submittedTaskCount = 0;
    this->updatedCutCount = 0;
    this->skippedMeshCount = 0;
    this->diagnostics = "";
    this->compactEntryNext = this->compactEntryFirst;
    this->compactEntryTotal = 0;
    this->deferredCompactEntries = FALSE;
    this->traverse(node);
    this->deferredCompactEntries =
	this->compactEntryNext < this->compactEntryTotal ? TRUE : FALSE;
}

void
SoBRLMeshLodSubmitAction::nodeAction(SoAction *action, SoNode *node)
{
    if (node->isOfType(SoGroup::getClassTypeId()))
	node->doAction(action);
}

static const char *
mesh_lod_source_leaf_name(const SoBRLDatabaseSource *source)
{
    if (!source)
	return "";
    const char *path = source->path.getValue().getString();
    if (!path || !path[0])
	return "";
    const char *slash = strrchr(path, '/');
    return (slash && slash[1]) ? slash + 1 : path;
}

static int
mesh_lod_source_is_threshold_bot(const SoBRLDatabaseSource *source)
{
    if (!source || source->lodBotThreshold.getValue() == 0)
	return 0;

    struct db_i *dbip = source->getDatabase();
    const char *name = mesh_lod_source_leaf_name(source);
    if (!dbip || !name || !name[0])
	return 0;

    struct directory *dp = db_lookup(dbip, name, LOOKUP_QUIET);
    if (dp == RT_DIR_NULL ||
	dp->d_minor_type != DB5_MINORTYPE_BRLCAD_BOT)
	return 0;

    struct rt_db_internal intern;
    RT_DB_INTERNAL_INIT(&intern);
    if (rt_db_get_internal(&intern, dp, dbip, NULL) < 0)
	return 0;

    int ret = 0;
    if (intern.idb_type == ID_BOT && intern.idb_ptr) {
	const struct rt_bot_internal *bot =
	    static_cast<const struct rt_bot_internal *>(intern.idb_ptr);
	if (bot && bot->num_faces >=
	    static_cast<size_t>(source->lodBotThreshold.getValue()))
	    ret = 1;
    }
    rt_db_free_internal(&intern);
    return ret;
}

static void
mesh_lod_make_source_request(const SoBRLDatabaseSource *source,
			     BObolLodRequest &request,
			     const char *databaseId,
			     uint64_t databaseRevision,
			     uint64_t viewRevision,
			     uint64_t policyRevision,
			     int requestDrawMode,
			     const char *providerId,
			     const char *providerVersion,
			     int qualityTier)
{
    request.clear();
    request.databaseId = databaseId ? databaseId : "";
    request.databaseRevision = databaseRevision;
    request.sourceRevision = source ? source->sourceRevision.getValue() : 0;
    request.objectPath = source ? source->path.getValue() : SbString("");
    request.objectName = mesh_lod_source_leaf_name(source);
    request.viewRevision = viewRevision;
    request.policyRevision = policyRevision;
    request.drawMode = requestDrawMode;
    request.lodPolicy = 0;
    request.providerId = providerId ? providerId : "";
    request.providerVersion = providerVersion ? providerVersion : "";
    request.qualityTier = qualityTier;
    if (source) {
	SbBox3f bounds;
	if (source->getSourceBounds(bounds))
	    request.bounds = bounds;
    }
}

/* Convert an occurrence-local bound to a quantized screen-space PoP demand.
 * Returning false means the box is wholly outside the current frustum and no
 * display payload should be loaded for this view. */
static SbBool
mesh_lod_apply_projected_demand(BObolLodRequest &request,
	const SbBox3f &localBounds, const SbMatrix &localToRoot,
	const SbViewVolume &viewVolume, const struct bv_view_info &view,
	float targetPixelError)
{
    if (localBounds.isEmpty())
	return TRUE;

    SbBox3f worldBounds;
    worldBounds.makeEmpty();
    SbVec3f projectedMin(FLT_MAX, FLT_MAX, FLT_MAX);
    SbVec3f projectedMax(-FLT_MAX, -FLT_MAX, -FLT_MAX);
    SbBool haveProjected = FALSE;
    const SbVec3f bmin = localBounds.getMin();
    const SbVec3f bmax = localBounds.getMax();
    SbVec3f worldCorners[8];
    for (int corner = 0; corner < 8; corner++) {
	const SbVec3f local(
	    (corner & 1) ? bmax[0] : bmin[0],
	    (corner & 2) ? bmax[1] : bmin[1],
	    (corner & 4) ? bmax[2] : bmin[2]);
	localToRoot.multVecMatrix(local, worldCorners[corner]);
	worldBounds.extendBy(worldCorners[corner]);
    }
    if (worldBounds.isEmpty() || !viewVolume.intersect(worldBounds))
	return FALSE;

    for (const SbVec3f &world : worldCorners) {
	SbVec3f projected;
	viewVolume.projectToScreen(world, projected);
	if (!std::isfinite(projected[0]) || !std::isfinite(projected[1]))
	    continue;
	haveProjected = TRUE;
	projectedMin[0] = std::min(projectedMin[0], projected[0]);
	projectedMin[1] = std::min(projectedMin[1], projected[1]);
	projectedMax[0] = std::max(projectedMax[0], projected[0]);
	projectedMax[1] = std::max(projectedMax[1], projected[1]);
    }
    if (!haveProjected)
	return TRUE;

    /* Visibility is clipped by the frustum test above, but geometric error
     * must use the full projected extent.  Clipping the extent to the
     * viewport made a very large part that was only partly visible look no
     * larger than the window, selecting a cut several levels too coarse for
     * the visible sliver. */
    const float widthPixels =
	std::max(0.0f, projectedMax[0] - projectedMin[0]) *
	static_cast<float>(std::max(1, view.width));
    const float heightPixels =
	std::max(0.0f, projectedMax[1] - projectedMin[1]) *
	static_cast<float>(std::max(1, view.height));
    const float diameter = std::max(widthPixels, heightPixels);
    const float policyScale =
	(std::isfinite(view.lod.scale) && view.lod.scale > 0.0) ?
	static_cast<float>(view.lod.scale) : 1.0f;
    const float pixelError =
	(std::isfinite(targetPixelError) && targetPixelError > 0.0f) ?
	targetPixelError / policyScale : 1.0f / policyScale;
    int level = 0;
    if (diameter > pixelError)
	level = static_cast<int>(std::ceil(std::log2(diameter / pixelError)));
    /* PoP coordinates are 16-bit, hence levels 0..15.  The cache clamps this
     * to the populated terminal level of the complete hierarchy. */
    level = std::max(0, std::min(15, level));
    request.projectedPixelDiameter = diameter;
    request.targetPixelError = pixelError;
    request.requestedLevel = level;
    return TRUE;
}

static void
mesh_lod_apply_level_hysteresis(BObolLodRequest &request, int activeLevel)
{
    if (activeLevel < 0 || request.requestedLevel < 0 ||
	request.requestedLevel == activeLevel || request.targetPixelError <= 0.0f)
	return;

    const double ratio = static_cast<double>(request.projectedPixelDiameter) /
	static_cast<double>(request.targetPixelError);
    const double upper = std::ldexp(1.0, activeLevel) * 1.25;
    const double lower = activeLevel > 0 ?
	std::ldexp(1.0, activeLevel - 1) * 0.75 : 0.0;
    if ((request.requestedLevel > activeLevel && ratio <= upper) ||
	(request.requestedLevel < activeLevel && ratio >= lower))
	request.requestedLevel = activeLevel;
}

static void
mesh_lod_trace_projected_request(const BObolLodRequest &request,
	const SbBox3f &localBounds, const SbMatrix &localToRoot,
	int activeLevel)
{
    const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (!filter || !filter[0])
	return;
    const char *name = request.objectName.getString();
    const char *path = request.objectPath.getString();
    if ((!name || !strstr(name, filter)) &&
	(!path || !strstr(path, filter)))
	return;

    const SbVec3f bmin = localBounds.isEmpty() ?
	SbVec3f(0.0f, 0.0f, 0.0f) : localBounds.getMin();
    const SbVec3f bmax = localBounds.isEmpty() ?
	SbVec3f(0.0f, 0.0f, 0.0f) : localBounds.getMax();
    const float *m = localToRoot[0];
    bu_log("BObol LoD trace object=%s path=%s occurrence=%s "
	   "pixels=%.9g error=%.9g request_level=%d active_level=%d "
	   "bounds=[%.9g %.9g %.9g]-[%.9g %.9g %.9g] "
	   "matrix=[%.9g %.9g %.9g %.9g; %.9g %.9g %.9g %.9g; "
	   "%.9g %.9g %.9g %.9g; %.9g %.9g %.9g %.9g]\n",
	   name ? name : "", path ? path : "",
	   request.occurrenceKey.getString(),
	   request.projectedPixelDiameter, request.targetPixelError,
	   request.requestedLevel, activeLevel,
	   bmin[0], bmin[1], bmin[2], bmax[0], bmax[1], bmax[2],
	   m[0], m[1], m[2], m[3], m[4], m[5], m[6], m[7],
	   m[8], m[9], m[10], m[11], m[12], m[13], m[14], m[15]);
}

static SbBool
mesh_lod_geometry_key_matches(const SbString &payloadKey,
	const BObolLodRequest &request)
{
    if (payloadKey.getLength() == 0)
	return FALSE;
    const BObolLodCacheKey key = bobol_lod_geometry_cache_key(request);
    return key.isValid() &&
	bu_strcmp(payloadKey.getString(), key.value.getString()) == 0;
}

static SbBool
mesh_lod_cad_payload_is_resident(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request)
{
    return payload && payload->isValid() &&
	(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) &&
	payload->providerStatus == BOBOL_LOD_PROVIDER_READY &&
	mesh_lod_geometry_key_matches(payload->cacheKey, request) &&
	payload->activeLevel == request.requestedLevel;
}

static SbBool
mesh_lod_cad_payload_has_same_source(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request)
{
    if (!payload || payload->requestedLevel < 0)
	return FALSE;
    BObolLodRequest payloadRequest = request;
    payloadRequest.requestedLevel = payload->requestedLevel;
    return mesh_lod_cad_payload_is_resident(payload, payloadRequest);
}

static SbBool
mesh_lod_mesh_payload_is_resident(
	const BObolViewLodState::MeshPayload *payload,
	const BObolLodRequest &request)
{
    return payload && payload->isValid() &&
	payload->resultKind == BOBOL_LOD_RESULT_MESH &&
	payload->providerStatus == BOBOL_LOD_PROVIDER_READY &&
	mesh_lod_geometry_key_matches(payload->cacheKey, request) &&
	payload->activeLevel == request.requestedLevel;
}

static SbBool
mesh_lod_mesh_payload_has_same_source(
	const BObolViewLodState::MeshPayload *payload,
	const BObolLodRequest &request)
{
    if (!payload || payload->requestedLevel < 0)
	return FALSE;
    BObolLodRequest payloadRequest = request;
    payloadRequest.requestedLevel = payload->requestedLevel;
    return mesh_lod_mesh_payload_is_resident(payload, payloadRequest);
}

static uint32_t
mesh_lod_debug_delay_milliseconds(const char *env_name)
{
    const char *delay = getenv(env_name);
    if (!delay || !delay[0])
	return 0;

    errno = 0;
    char *end = NULL;
    unsigned long value = strtoul(delay, &end, 10);
    if (errno != 0 || end == delay || value > UINT32_MAX)
	return 0;

    return (uint32_t)value;
}

void
SoBRLMeshLodSubmitAction::databaseSourceAction(SoAction *action, SoNode *node)
{
    SoBRLMeshLodSubmitAction *submitAction =
	static_cast<SoBRLMeshLodSubmitAction *>(action);
    SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);

    if (!source->hasCompactInstanceIndex()) {
	source->doAction(action);
	return;
    }

    if (source->isCompactOccurrenceRegistry()) {
	const int sourceDrawMode =
	    source->representationMode.getValue() ==
		SoBRLDatabaseSource::REPRESENTATION_HIDDEN_LINE ?
	    BOBOL_LOD_DRAW_HIDDEN_LINE :
	    (source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	    BOBOL_LOD_DRAW_SHADED : BOBOL_LOD_DRAW_WIRE);
	std::unordered_map<std::string,
	    const BObolViewLodState::CadPayload *> activePayloads;
	if (submitAction->viewState) {
	    std::vector<const BObolViewLodState::CadPayload *> payloads;
	    submitAction->viewState->findCadPayloads(source, payloads);
	    activePayloads.reserve(payloads.size());
	    for (const BObolViewLodState::CadPayload *payload : payloads) {
		if (payload && payload->sourceInstanceKey.getLength() > 0)
		    activePayloads[payload->sourceInstanceKey.getString()] =
			payload;
	    }
	}
	const int count = source->getCompactInstanceCount();
	const size_t sourceFirst = submitAction->compactEntryTotal;
	submitAction->compactEntryTotal += static_cast<size_t>(count);
	const size_t rangeFirst = submitAction->compactEntryFirst;
	const size_t rangeLast = submitAction->compactEntryLimit == SIZE_MAX ||
	    rangeFirst > SIZE_MAX - submitAction->compactEntryLimit ? SIZE_MAX :
	    rangeFirst + submitAction->compactEntryLimit;
	const size_t localFirst = rangeFirst > sourceFirst ?
	    std::min(static_cast<size_t>(count), rangeFirst - sourceFirst) : 0;
	const size_t sourceLast = sourceFirst + static_cast<size_t>(count);
	const size_t localLast = rangeLast < sourceLast ?
	    (rangeLast > sourceFirst ? rangeLast - sourceFirst : 0) :
	    static_cast<size_t>(count);
	for (size_t i = localFirst; i < localLast; i++) {
	    BObolCompactInstanceHandle handle;
	    BObolCompactInstanceSummary summary;
	    if (!source->getCompactInstanceHandle(static_cast<int>(i), handle) ||
		!source->getCompactInstanceSummary(handle, summary) ||
		!summary.valid || !summary.visible)
		continue;

	    submitAction->visitedMeshCount++;
	    const SbString target = summary.path.getLength() > 0 ?
		summary.path : source->path.getValue();
	    if (!submitAction->service || !submitAction->service->isRunning()) {
		submitAction->skippedMeshCount++;
		submitAction->appendDiagnostic(target,
		    "LoD service is not running");
		continue;
	    }
	    if (!submitAction->dbip) {
		submitAction->skippedMeshCount++;
		submitAction->appendDiagnostic(target,
		    "LoD submit action has no database");
		continue;
	    }
	    /* Native shaded geometry is already the stable presentation for
	     * analytic primitives such as TGCs.  A source-wide nonzero BoT
	     * threshold does not make those meshes PoP-backed.  Only occurrences
	     * carrying the explicit source-mesh request contract may enter the
	     * cache provider; otherwise a guaranteed cache miss can displace the
	     * useful native presentation with its structural box fallback. */
	    const bool lodEligible =
		summary.lodBacked && summary.sourceMeshRequestValid;
	    const bool meshEligible = lodEligible ||
		(!submitAction->requireLodBacked &&
		 (summary.meshGeometry ||
		  (source->realizationRoleFlags.getValue() &
		   SoBRLDatabaseSource::REALIZATION_ROLE_MESH)));
	    if (!meshEligible) {
		submitAction->skippedMeshCount++;
		continue;
	    }

	    BObolLodRequest request;
	    request.databaseId = submitAction->databaseId;
	    request.databaseRevision = submitAction->databaseRevision;
	    request.sourceRevision = source->sourceRevision.getValue();
	    /* The compact part may currently be only an AABB and will be
	     * replaced by richer presentation data.  Its part identity is not
	     * source-content identity.  Leaving this unset makes the stable
	     * geometry key use database/source revisions and object identity,
	     * which survive box -> PoP presentation replacement. */
	    request.sourceContentHash = 0;
	    request.objectPath = target;
	    request.objectName = summary.sourceName;
	    request.occurrenceKey = summary.sourceInstanceKey;
	    if (request.objectName.getLength() == 0)
		request.objectName = mesh_lod_source_leaf_name(source);
	    request.viewRevision = submitAction->viewRevision;
	    request.policyRevision = submitAction->policyRevision;
	    request.drawMode = sourceDrawMode;
	    request.providerId = submitAction->providerId;
	    request.providerVersion = submitAction->providerVersion;
	    request.qualityTier = submitAction->qualityTier;
	    request.bounds = summary.localBounds;
	    /* Compact summaries already report the complete geometry-to-root
	     * transform, including the source draw matrix.  Applying the source
	     * matrix again corrupts screen bounds (and can make a visible leaf
	     * look tiny or off-screen to the LoD selector). */
	    SbMatrix localToRoot = summary.localToSource;
	    if (submitAction->useViewVolume &&
		!mesh_lod_apply_projected_demand(request, summary.localBounds,
		    localToRoot, submitAction->viewVolume,
		    submitAction->view, submitAction->targetPixelError)) {
		submitAction->skippedMeshCount++;
		continue;
	    }
	    const auto activeIt =
		activePayloads.find(request.occurrenceKey.getString());
	    const BObolViewLodState::CadPayload *activePayload =
		activeIt != activePayloads.end() ? activeIt->second : NULL;
	    if (submitAction->useForcedLevel)
		request.requestedLevel = submitAction->forcedLevel;
	    else if (activePayload)
		mesh_lod_apply_level_hysteresis(request,
		    activePayload->activeLevel);
	    mesh_lod_trace_projected_request(request, summary.localBounds,
		localToRoot, activePayload ? activePayload->activeLevel : -1);

	    /* A camera epoch is not a geometry invalidation.  Keep the current
	     * cut when it already satisfies this demand, and never coarsen it
	     * during motion.  Settled zoom-out may explicitly request a smaller
	     * cut to reclaim display memory. */
	    if (!submitAction->useForcedLevel &&
		mesh_lod_cad_payload_is_resident(activePayload, request)) {
		submitAction->skippedMeshCount++;
		submitAction->appendDiagnostic(target,
		    "current LoD request is already resident");
		continue;
	    }
	    if (!submitAction->useForcedLevel &&
		!submitAction->allowLevelDowngrade && activePayload &&
		activePayload->activeLevel > request.requestedLevel &&
		mesh_lod_cad_payload_has_same_source(activePayload, request)) {
		submitAction->skippedMeshCount++;
		submitAction->appendDiagnostic(target,
		    "finer resident LoD retained during interaction");
		continue;
	    }
	    if (submitAction->reset == 0 && activePayload &&
		mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
		activePayload->activeLevel != request.requestedLevel &&
		submitAction->viewState->retargetCadPayload(activePayload,
		    request.requestedLevel, request.viewRevision,
		    request.policyRevision)) {
		const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (filter && filter[0] &&
		    (strstr(request.objectName.getString(), filter) ||
		     strstr(request.objectPath.getString(), filter)))
		    bu_log("BObol LoD submit trace object=%s level=%d "
			   "retargeted_resident=%d\n",
			   request.objectName.getString(),
			   request.requestedLevel,
			   activePayload->progressiveMesh ?
			       activePayload->progressiveMesh->residentLevel() : -1);
		submitAction->updatedCutCount++;
		continue;
	    }

	const SbBool suppressActiveDuplicate =
		(!submitAction->useForcedLevel && submitAction->reset == 0) ?
		TRUE : FALSE;
	    if (suppressActiveDuplicate &&
		submitAction->service->hasActiveRequest(request)) {
		const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (filter && filter[0] &&
		    (strstr(request.objectName.getString(), filter) ||
		     strstr(request.objectPath.getString(), filter)))
		    bu_log("BObol LoD submit trace object=%s level=%d "
			   "suppressed=active\n",
			   request.objectName.getString(),
			   request.requestedLevel);
		submitAction->skippedMeshCount++;
		continue;
	    }

	    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
	    provider->service = submitAction->service;
	    provider->dbip = submitAction->dbip;
	    provider->view = submitAction->view;
	    provider->useView = TRUE;
	    provider->refreshMissing = submitAction->refreshMissing;
	    provider->useForcedLevel = submitAction->useForcedLevel;
	    provider->forcedLevel = submitAction->forcedLevel;
	    provider->shrinkAfterCopy = FALSE;
		    /* A shared source asset may serve many occurrences at
		     * different levels.  Trimming is therefore a post-generation
		     * aggregate operation, never a leaf-request side effect. */
		    provider->compactResident = FALSE;
	    provider->reset = submitAction->reset;

	    BObolLodTask task;
	    task.generation = submitAction->generation;
	    task.request = request;
	    task.realize = bobol_mesh_lod_provider_task;
	    task.realizeData = provider;
	    task.realizeDataFree = bobol_mesh_lod_provider_free;
	    task.debugDelayMilliseconds = mesh_lod_debug_delay_milliseconds(
		"BOBOL_LOD_TASK_DELAY_MS");
	    const uint64_t taskId = suppressActiveDuplicate ?
		submitAction->service->submitIfNotActive(task) :
		submitAction->service->submit(task);
	    if (!taskId) {
		const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (filter && filter[0] &&
		    (strstr(request.objectName.getString(), filter) ||
		     strstr(request.objectPath.getString(), filter)))
		    bu_log("BObol LoD submit trace object=%s level=%d "
			   "rejected=service\n",
			   request.objectName.getString(),
			   request.requestedLevel);
		bobol_mesh_lod_provider_free(provider);
		submitAction->skippedMeshCount++;
	    } else {
		const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (filter && filter[0] &&
		    (strstr(request.objectName.getString(), filter) ||
		     strstr(request.objectPath.getString(), filter)))
		    bu_log("BObol LoD submit trace object=%s level=%d "
			   "task=%llu\n",
			   request.objectName.getString(),
			   request.requestedLevel,
			   static_cast<unsigned long long>(taskId));
		submitAction->submittedTaskCount++;
	    }
	}
	submitAction->compactEntryNext = std::max(
	    submitAction->compactEntryNext, sourceFirst + localLast);
	/* Auxiliary overlays remain ordinary child nodes and retain their own
	 * scheduling behavior. */
	source->doAction(action);
	return;
    }

    submitAction->visitedMeshCount++;
    SbString target = source->path.getValue();

    if (!submitAction->service || !submitAction->service->isRunning()) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD service is not running");
	source->doAction(action);
	return;
    }

    if (!submitAction->dbip) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD submit action has no database");
	source->doAction(action);
	return;
    }

    const int sourceDrawMode =
	source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ?
	BOBOL_LOD_DRAW_SHADED : BOBOL_LOD_DRAW_WIRE;
    const int roleFlags = source->realizationRoleFlags.getValue();
    const int lodEligible = mesh_lod_source_is_threshold_bot(source);
    const int meshEligible = lodEligible ||
	(!submitAction->requireLodBacked &&
	 ((roleFlags & SoBRLDatabaseSource::REALIZATION_ROLE_MESH) ||
	  source->drawMode.getValue() == SoBRLDatabaseSource::SHADED ||
	  source->representationMode.getValue() ==
	      SoBRLDatabaseSource::REPRESENTATION_SHADED_BOTS ||
	  source->representationMode.getValue() ==
	      SoBRLDatabaseSource::REPRESENTATION_SHADED));

    if (!meshEligible) {
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }

    BObolLodRequest request;
    mesh_lod_make_source_request(source, request,
	submitAction->databaseId.getString(),
	submitAction->databaseRevision,
	submitAction->viewRevision,
	submitAction->policyRevision,
	sourceDrawMode,
	submitAction->providerId.getString(),
	submitAction->providerVersion.getString(),
	submitAction->qualityTier);
    if (submitAction->useViewVolume) {
	SbMatrix localToRoot = source->drawMatrixValid.getValue() ?
	    source->drawMatrix.getValue() : SbMatrix::identity();
	if (!mesh_lod_apply_projected_demand(request, request.bounds,
		localToRoot, submitAction->viewVolume, submitAction->view,
		submitAction->targetPixelError)) {
	    submitAction->skippedMeshCount++;
	    source->doAction(action);
	    return;
	}
    }
    const BObolViewLodState::CadPayload *activePayload =
	submitAction->viewState ? submitAction->viewState->findCad(source) : NULL;
    if (submitAction->useForcedLevel)
	request.requestedLevel = submitAction->forcedLevel;
    else if (activePayload)
	mesh_lod_apply_level_hysteresis(request, activePayload->activeLevel);
    if (!submitAction->useForcedLevel &&
	mesh_lod_cad_payload_is_resident(activePayload, request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current LoD request is already resident");
	source->doAction(action);
	return;
    }
    if (!submitAction->useForcedLevel &&
	!submitAction->allowLevelDowngrade && activePayload &&
	activePayload->activeLevel > request.requestedLevel &&
	mesh_lod_cad_payload_has_same_source(activePayload, request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "finer resident LoD retained during interaction");
	source->doAction(action);
	return;
    }
    if (submitAction->reset == 0 && activePayload &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
	activePayload->activeLevel != request.requestedLevel &&
	submitAction->viewState->retargetCadPayload(activePayload,
	    request.requestedLevel, request.viewRevision,
	    request.policyRevision)) {
	submitAction->updatedCutCount++;
	source->doAction(action);
	return;
    }

    const SbBool suppressActiveDuplicate =
	(!submitAction->useForcedLevel && submitAction->reset == 0) ?
	TRUE : FALSE;
    if (suppressActiveDuplicate &&
	submitAction->service->hasActiveRequest(request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current compact LoD request is already active");
	source->doAction(action);
	return;
    }

    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
    provider->service = submitAction->service;
    provider->dbip = submitAction->dbip;
    provider->view = submitAction->view;
    provider->useView = TRUE;
    provider->refreshMissing = submitAction->refreshMissing;
    provider->useForcedLevel = submitAction->useForcedLevel;
    provider->forcedLevel = submitAction->forcedLevel;
    provider->shrinkAfterCopy = FALSE;
	    provider->compactResident = FALSE;
    provider->reset = submitAction->reset;

    BObolLodTask task;
    task.generation = submitAction->generation;
    task.request = request;
    task.realize = bobol_mesh_lod_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = bobol_mesh_lod_provider_free;
    task.debugDelayMilliseconds =
	mesh_lod_debug_delay_milliseconds("BOBOL_LOD_TASK_DELAY_MS");

    uint64_t taskId = suppressActiveDuplicate ?
		      submitAction->service->submitIfNotActive(task) :
		      submitAction->service->submit(task);
    if (taskId == 0) {
	bobol_mesh_lod_provider_free(provider);
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "LoD service rejected compact mesh task");
    } else {
	submitAction->submittedTaskCount++;
    }

    source->doAction(action);
}

void
SoBRLMeshLodSubmitAction::meshShapeAction(SoAction *action, SoNode *node)
{
    SoBRLMeshLodSubmitAction *submitAction =
	static_cast<SoBRLMeshLodSubmitAction *>(action);
    SoBRLMeshShape *shape = static_cast<SoBRLMeshShape *>(node);

    submitAction->visitedMeshCount++;

    SbString target = shape->sourcePath.getValue().getLength() > 0 ?
		      shape->sourcePath.getValue() : shape->sourceName.getValue();

    if (submitAction->requireLodBacked && !shape->isLodBackedMesh()) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "mesh is not LoD-backed");
	return;
    }

    if (!submitAction->service || !submitAction->service->isRunning()) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD service is not running");
	return;
    }

    if (!submitAction->dbip) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target, "LoD submit action has no database");
	return;
    }

    BObolLodRequest request;
    shape->makeLodRequest(request,
			  submitAction->databaseId.getString(),
			  submitAction->databaseRevision,
			  submitAction->viewRevision,
			  submitAction->policyRevision,
			  BOBOL_LOD_DRAW_SHADED,
			  submitAction->providerId.getString(),
			  submitAction->providerVersion.getString(),
			  submitAction->qualityTier);
    if (submitAction->useViewVolume &&
	!mesh_lod_apply_projected_demand(request, request.bounds,
	    SbMatrix::identity(), submitAction->viewVolume, submitAction->view,
	    submitAction->targetPixelError)) {
	submitAction->skippedMeshCount++;
	return;
    }

    const BObolViewLodState::MeshPayload *viewPayload =
	submitAction->viewState ?
	submitAction->viewState->findMesh(shape) : NULL;
    if (submitAction->useForcedLevel)
	request.requestedLevel = submitAction->forcedLevel;
    else if (viewPayload)
	mesh_lod_apply_level_hysteresis(request, viewPayload->activeLevel);
    if (!submitAction->useForcedLevel &&
	mesh_lod_mesh_payload_is_resident(viewPayload, request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current LoD request is already resident");
	return;
    }
    if (!submitAction->useForcedLevel &&
	!submitAction->allowLevelDowngrade && viewPayload &&
	viewPayload->activeLevel > request.requestedLevel &&
	mesh_lod_mesh_payload_has_same_source(viewPayload, request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "finer resident LoD retained during interaction");
	return;
    }
    if (submitAction->reset == 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request) &&
	viewPayload->activeLevel != request.requestedLevel &&
	request.requestedLevel >= 0 &&
	submitAction->viewState->retargetMeshPayload(viewPayload,
	    request.requestedLevel, request.viewRevision,
	    request.policyRevision)) {
	submitAction->updatedCutCount++;
	return;
    }
    if (!submitAction->useForcedLevel && submitAction->reset == 0 &&
	request.requestedLevel < 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current LoD request is already resident");
	return;
    }

    const SbBool suppressActiveDuplicate =
	(!submitAction->useForcedLevel && submitAction->reset == 0) ?
	TRUE : FALSE;
    if (suppressActiveDuplicate &&
	submitAction->service->hasActiveRequest(request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current LoD request is already active");
	return;
    }

    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
    provider->service = submitAction->service;
    provider->dbip = submitAction->dbip;
    provider->view = submitAction->view;
    provider->useView = TRUE;
    provider->refreshMissing = submitAction->refreshMissing;
    provider->useForcedLevel = submitAction->useForcedLevel;
    provider->forcedLevel = submitAction->forcedLevel;
    provider->shrinkAfterCopy = FALSE;
	    provider->compactResident = FALSE;
    provider->reset = submitAction->reset;

    BObolLodTask task;
    task.generation = submitAction->generation;
    task.request = request;
    task.realize = bobol_mesh_lod_provider_task;
    task.realizeData = provider;
    task.realizeDataFree = bobol_mesh_lod_provider_free;
    task.debugDelayMilliseconds =
	mesh_lod_debug_delay_milliseconds("BOBOL_LOD_TASK_DELAY_MS");

    uint64_t taskId = suppressActiveDuplicate ?
		      submitAction->service->submitIfNotActive(task) :
		      submitAction->service->submit(task);
    if (taskId == 0) {
	bobol_mesh_lod_provider_free(provider);
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       suppressActiveDuplicate &&
				       submitAction->service->hasActiveRequest(request) ?
				       "current LoD request is already active" :
				       "LoD service rejected mesh task");
	return;
    }

    submitAction->submittedTaskCount++;
}

void
SoBRLMeshLodSubmitAction::appendDiagnostic(const SbString &target,
	const char *message)
{
    if (this->diagnostics.getLength() > 0)
	this->diagnostics += "\n";

    this->diagnostics += target.getLength() > 0 ? target : SbString("<unknown>");
    this->diagnostics += ": ";
    this->diagnostics += message ? message : "LoD submit diagnostic";
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
