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
    allowRetainedRefinement(TRUE),
    refinementFaceBudget(SIZE_MAX),
    refinementFaceBudgetUsed(0),
    refinementBudgetBlockedCount(0),
    viewState(NULL),
    compactEntryFirst(0),
    compactEntryLimit(SIZE_MAX),
    compactEntryNext(0),
    compactEntryTotal(0),
    compactEntryPlan(),
    compactEntryPlanSupplied(FALSE),
    submissionTaskLimit(SIZE_MAX),
    retainedSceneFaceBudget(SIZE_MAX),
    retainedSceneFaceBudgetUsed(0),
    retainedRecoveredOccurrences(),
    deferredCompactEntries(FALSE),
    visitedMeshCount(0),
    submittedTaskCount(0),
    updatedCutCount(0),
    pendingRetainedRefinementCount(0),
    skippedMeshCount(0),
    diagnosticCount(0),
    suppressedDiagnosticCount(0),
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
SoBRLMeshLodSubmitAction::setAllowRetainedRefinement(SbBool allow)
{
    this->allowRetainedRefinement = allow ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getAllowRetainedRefinement(void) const
{
    return this->allowRetainedRefinement;
}

void
SoBRLMeshLodSubmitAction::setRefinementFaceBudget(size_t additionalFaces)
{
    this->refinementFaceBudget = additionalFaces;
}

size_t
SoBRLMeshLodSubmitAction::getRefinementFaceBudget(void) const
{
    return this->refinementFaceBudget;
}

size_t
SoBRLMeshLodSubmitAction::getRefinementFaceBudgetUsed(void) const
{
    return this->refinementFaceBudgetUsed;
}

unsigned int
SoBRLMeshLodSubmitAction::getRefinementBudgetBlockedCount(void) const
{
    return this->refinementBudgetBlockedCount;
}

SbBool
SoBRLMeshLodSubmitAction::reserveRefinementFaces(
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int activeLevel, int nextLevel)
{
    if (!progressiveMesh || nextLevel <= activeLevel)
	return TRUE;

    const size_t activeFaces = activeLevel >= 0 ?
	progressiveMesh->hierarchyFaceCount(activeLevel) : 0;
    const size_t nextFaces =
	progressiveMesh->hierarchyFaceCount(nextLevel);
    const size_t additionalFaces =
	nextFaces > activeFaces ? nextFaces - activeFaces : 0;
    if (!additionalFaces)
	return TRUE;

    if (this->refinementFaceBudget != SIZE_MAX &&
	(this->refinementFaceBudgetUsed > this->refinementFaceBudget ||
	 additionalFaces >
	    this->refinementFaceBudget - this->refinementFaceBudgetUsed)) {
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD refinement reservation blocked "
		   "active_level=%d next_level=%d active_faces=%zu "
		   "next_faces=%zu added_faces=%zu used_faces=%zu "
		   "budget_faces=%zu\n",
		   activeLevel, nextLevel, activeFaces, nextFaces,
		   additionalFaces, this->refinementFaceBudgetUsed,
		   this->refinementFaceBudget);
	this->refinementBudgetBlockedCount++;
	return FALSE;
    }

    this->refinementFaceBudgetUsed =
	additionalFaces > SIZE_MAX - this->refinementFaceBudgetUsed ?
	SIZE_MAX : this->refinementFaceBudgetUsed + additionalFaces;
    return TRUE;
}

SbBool
SoBRLMeshLodSubmitAction::reserveInitialFaces(
    uint64_t sourceFaces, size_t &providerFaceAllowance)
{
    providerFaceAllowance = 0;
    if (this->refinementFaceBudget == SIZE_MAX)
	return TRUE;

    /* Four thousand faces is intentionally a reservation, not a demanded
     * cut.  It is large enough to avoid admitting hundreds of unknown sparse
     * hierarchies against a tiny scene budget, while the provider may select
     * a cheaper populated minimum.  A small source reserves its exact full
     * population. */
    static const size_t provisionalFirstCutFaces = 4096;
    const size_t knownFaces =
	sourceFaces > static_cast<uint64_t>(SIZE_MAX) ? SIZE_MAX :
	static_cast<size_t>(sourceFaces);
    const size_t reserveFaces = knownFaces ?
	std::min(knownFaces, provisionalFirstCutFaces) :
	provisionalFirstCutFaces;

    /*
     * The scene face allowance controls optional refinement, not coverage.
     * Once a visible occurrence's minimum PoP mesh is available, falling
     * back to its structural box is both visually disruptive and more
     * expensive semantically than drawing the coherent minimum prefix.  The
     * service's task/result and working-set limits bound cold preparation;
     * consume any remaining refinement allowance here, but never reject the
     * first mesh solely because the learned draw budget is exhausted.
     */
    if (this->refinementFaceBudgetUsed > this->refinementFaceBudget ||
	reserveFaces >
	    this->refinementFaceBudget - this->refinementFaceBudgetUsed) {
	this->refinementFaceBudgetUsed = this->refinementFaceBudget;
	providerFaceAllowance = provisionalFirstCutFaces;
	return TRUE;
    }

    this->refinementFaceBudgetUsed += reserveFaces;
    /* The reservation may use an exact small source count, while a legacy
     * cache's hierarchy metadata can conservatively report a slightly larger
     * populated minimum.  Keep the provider's local first-cut ceiling at the
     * provisional quantum; aggregate admission is still charged by the
     * source count. */
    providerFaceAllowance = std::max(reserveFaces,
	provisionalFirstCutFaces);
    return TRUE;
}

SbBool
SoBRLMeshLodSubmitAction::reserveInitialFaceCost(size_t faceCount)
{
    if (this->refinementFaceBudget != SIZE_MAX &&
	(this->refinementFaceBudgetUsed > this->refinementFaceBudget ||
	 faceCount >
	    this->refinementFaceBudget - this->refinementFaceBudgetUsed)) {
	/* Shared resident assets obey the same coverage floor as provider
	 * results: budget optional enrichment, never box reintroduction. */
	this->refinementFaceBudgetUsed = this->refinementFaceBudget;
	return TRUE;
    }
    this->refinementFaceBudgetUsed =
	faceCount > SIZE_MAX - this->refinementFaceBudgetUsed ?
	SIZE_MAX : this->refinementFaceBudgetUsed + faceCount;
    return TRUE;
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

void
SoBRLMeshLodSubmitAction::setCompactEntryPlan(
    const std::vector<size_t> &entryIndices)
{
    this->compactEntryPlan = entryIndices;
    this->compactEntryPlanSupplied = TRUE;
}

void
SoBRLMeshLodSubmitAction::getCompactEntryPlan(
    std::vector<size_t> &entryIndices) const
{
    entryIndices = this->compactEntryPlan;
}

void
SoBRLMeshLodSubmitAction::setSubmissionTaskLimit(size_t taskCount)
{
    this->submissionTaskLimit = taskCount;
}

void
SoBRLMeshLodSubmitAction::setRetainedSceneFaceBudget(size_t totalFaces)
{
    this->retainedSceneFaceBudget = totalFaces;
}

size_t
SoBRLMeshLodSubmitAction::getRetainedSceneFaceBudgetUsed(void) const
{
    return this->retainedSceneFaceBudgetUsed;
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
SoBRLMeshLodSubmitAction::getPendingRetainedRefinementCount(void) const
{
    return this->pendingRetainedRefinementCount;
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
    this->pendingRetainedRefinementCount = 0;
    this->refinementFaceBudgetUsed = 0;
    this->refinementBudgetBlockedCount = 0;
    this->retainedSceneFaceBudgetUsed = 0;
    this->retainedRecoveredOccurrences.clear();
    this->skippedMeshCount = 0;
    this->diagnosticCount = 0;
    this->suppressedDiagnosticCount = 0;
    this->diagnostics = "";
    this->compactEntryNext = this->compactEntryFirst;
    this->compactEntryTotal = 0;
    this->deferredCompactEntries = FALSE;
    this->traverse(node);
    if (this->suppressedDiagnosticCount > 0) {
	if (this->diagnostics.getLength() > 0)
	    this->diagnostics += "\n";
	SbString summary;
	summary.sprintf("<summary>: %u additional per-occurrence diagnostics "
	    "suppressed", this->suppressedDiagnosticCount);
	this->diagnostics += summary;
    }
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

/* Use the same row-vector view-projection half-spaces as SoCADAssembly's
 * renderer.  SbViewVolume::intersect derives a second set of planes from its
 * parametric volume; at extreme close-up views the two tests can disagree,
 * leaving a structural AABB visible while the LoD scheduler refuses to load
 * its mesh.  The matrix test is the presentation contract: if Obol can draw
 * an instance, the scheduler must consider it visible. */
static SbBool
mesh_lod_box_intersects_render_frustum(const SbBox3f &worldBounds,
	const SbMatrix &viewProjection)
{
    if (worldBounds.isEmpty())
	return TRUE;

    const SbVec3f bmin = worldBounds.getMin();
    const SbVec3f bmax = worldBounds.getMax();
    for (int column = 0; column < 3; column++) {
	for (int signIndex = 0; signIndex < 2; signIndex++) {
	    const float sign = signIndex == 0 ? 1.0f : -1.0f;
	    const float a = sign * viewProjection[0][column] +
		viewProjection[0][3];
	    const float b = sign * viewProjection[1][column] +
		viewProjection[1][3];
	    const float c = sign * viewProjection[2][column] +
		viewProjection[2][3];
	    const float d = sign * viewProjection[3][column] +
		viewProjection[3][3];
	    const float x = a < 0.0f ? bmin[0] : bmax[0];
	    const float y = b < 0.0f ? bmin[1] : bmax[1];
	    const float z = c < 0.0f ? bmin[2] : bmax[2];
	    if (a * x + b * y + c * z + d < 0.0f)
		return FALSE;
	}
    }
    return TRUE;
}

/* Convert an occurrence-local bound to a quantized screen-space PoP demand.
 * Returning false means the box is wholly outside the actual render frustum
 * and no display payload should be loaded for this view. */
static SbBool
mesh_lod_apply_projected_demand(BObolLodRequest &request,
	const SbBox3f &localBounds, const SbMatrix &localToRoot,
	const SbViewVolume &viewVolume, const struct bv_view_info &view,
	float targetPixelError, const SbMatrix *cachedViewProjection = NULL)
{
    if (localBounds.isEmpty())
	return TRUE;

    /* SbViewVolume::projectToScreen rebuilds the double-precision view and
     * projection matrices on every point.  A compact scene evaluates eight
     * corners for every occurrence, so that innocent API use dominated a
     * camera epoch.  One immutable view-projection matrix is sufficient for
     * both renderer-contract frustum testing and normalized projection. */
    const SbMatrix localViewProjection =
	cachedViewProjection ? SbMatrix() : viewVolume.getMatrix();
    const SbMatrix &viewProjection =
	cachedViewProjection ? *cachedViewProjection : localViewProjection;
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
    if (worldBounds.isEmpty() ||
	!mesh_lod_box_intersects_render_frustum(worldBounds, viewProjection))
	return FALSE;

    for (const SbVec3f &world : worldCorners) {
	SbVec3f projected;
	viewProjection.multVecMatrix(world, projected);
	projected *= 0.5f;
	projected += SbVec3f(0.5f, 0.5f, 0.5f);
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
	const BObolLodRequest &request,
	const SbString *cachedGeometryKey = NULL)
{
    if (payloadKey.getLength() == 0)
	return FALSE;
    if (cachedGeometryKey)
	return cachedGeometryKey->getLength() > 0 &&
	    bu_strcmp(payloadKey.getString(),
		cachedGeometryKey->getString()) == 0 ? TRUE : FALSE;
    const BObolLodCacheKey generated = bobol_lod_geometry_cache_key(request);
    return generated.isValid() &&
	bu_strcmp(payloadKey.getString(), generated.value.getString()) == 0;
}

static SbBool
mesh_lod_cad_payload_is_resident(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request,
	const SbString *cachedGeometryKey = NULL)
{
    return payload && payload->isValid() &&
	(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) &&
	payload->providerStatus == BOBOL_LOD_PROVIDER_READY &&
	mesh_lod_geometry_key_matches(payload->cacheKey, request,
	    cachedGeometryKey) &&
	payload->activeLevel == request.requestedLevel;
}

static SbBool
mesh_lod_cad_payload_has_same_source(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request,
	const SbString *cachedGeometryKey = NULL)
{
    if (!payload || payload->requestedLevel < 0)
	return FALSE;
    return payload->isValid() &&
	(payload->resultKind == BOBOL_LOD_RESULT_MESH ||
	 payload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL) &&
	payload->providerStatus == BOBOL_LOD_PROVIDER_READY &&
	mesh_lod_geometry_key_matches(payload->cacheKey, request,
	    cachedGeometryKey);
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

static int
mesh_lod_available_cad_retarget_level(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request,
	SbBool incrementalRefinement)
{
    if (!payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() ||
	request.requestedLevel < 0)
	return request.requestedLevel;
    if (request.requestedLevel <= payload->activeLevel)
	return request.requestedLevel;
    /* residentLevel is the loaded cache prefix, but a finer PoP cut can also
     * be exactly drawable when adjacent levels add only quantization bits and
     * no vertices/faces.  Reuse the richest cut the retained asset can
     * actually draw before scheduling more I/O. */
    for (int level = request.requestedLevel;
	 level > payload->activeLevel; --level) {
	if (!payload->progressiveMesh->canDrawLevel(level))
	    continue;
	if (incrementalRefinement && payload->activeLevel >= 0 &&
	    level > payload->activeLevel + 1) {
	    const int nextLevel = payload->activeLevel + 1;
	    return payload->progressiveMesh->canDrawLevel(nextLevel) ?
		nextLevel : payload->activeLevel;
	}
	return level;
    }
    return payload->activeLevel;
}

static int
mesh_lod_available_mesh_retarget_level(
	const BObolViewLodState::MeshPayload *payload,
	const BObolLodRequest &request,
	SbBool incrementalRefinement)
{
    if (!payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() ||
	request.requestedLevel < 0)
	return request.requestedLevel;
    if (request.requestedLevel <= payload->activeLevel)
	return request.requestedLevel;
    for (int level = request.requestedLevel;
	 level > payload->activeLevel; --level) {
	if (!payload->progressiveMesh->canDrawLevel(level))
	    continue;
	if (incrementalRefinement && payload->activeLevel >= 0 &&
	    level > payload->activeLevel + 1) {
	    const int nextLevel = payload->activeLevel + 1;
	    return payload->progressiveMesh->canDrawLevel(nextLevel) ?
		nextLevel : payload->activeLevel;
	}
	return level;
    }
    return payload->activeLevel;
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
	std::unordered_map<std::string, SbString> geometryKeysByAsset;
	const int count = source->getCompactInstanceCount();
	SbMatrix compactViewProjection;
	const SbMatrix *compactViewProjectionPtr = NULL;
	if (submitAction->useViewVolume) {
	    compactViewProjection = submitAction->viewVolume.getMatrix();
	    compactViewProjectionPtr = &compactViewProjection;
	}
	/* Provider submission is intentionally consumed in bounded windows, but
	 * an over-budget retained scene cannot wait for those windows: camera
	 * changes restart their cursor, so a 50k-leaf scene could repeatedly
	 * inspect the same 2k structural entries while millions of already
	 * resident triangles remained live.
	 *
	 * Recover the displayed population independently.  This scans only
	 * active view bindings, projects them by constant-time occurrence
	 * lookup, and retargets coherent minimum PoP cuts in screen-priority
	 * order.  Minimum mesh coverage is a floor, not an optional budget
	 * allocation: overload may lower every retained draw prefix, but must
	 * never expose a structural box after a mesh has become available. */
	if (submitAction->compactEntryFirst == 0 &&
	    submitAction->retainedSceneFaceBudget != SIZE_MAX &&
	    submitAction->viewState) {
	    struct RetainedCandidate {
		const BObolViewLodState::CadPayload *payload;
		SbString occurrenceKey;
		float projectedPixels;
		int emphasis;
		int minimumLevel;
		int requestedLevel;
		size_t minimumFaces;
	    };
	    std::vector<const BObolViewLodState::CadPayload *> activePayloads;
	    std::vector<RetainedCandidate> retainedCandidates;
	    submitAction->viewState->findCadPayloadsUnordered(
		source, activePayloads);
	    retainedCandidates.reserve(activePayloads.size());
	    for (const BObolViewLodState::CadPayload *payload :
		activePayloads) {
		if (!payload || !payload->isValid())
		    continue;
		BObolCompactLodPlanningSummary summary;
		if (!source->getCompactLodPlanningSummaryForKey(
			payload->sourceInstanceKey.getString(), summary) ||
		    !summary.valid || !summary.visible) {
		    if (submitAction->viewState->removeCadPayload(payload))
			submitAction->updatedCutCount++;
		    continue;
		}
		BObolLodRequest projected;
		projected.bounds = summary.localBounds;
		if (submitAction->useViewVolume &&
		    !mesh_lod_apply_projected_demand(projected,
			summary.localBounds, summary.localToSource,
			submitAction->viewVolume, submitAction->view,
			submitAction->targetPixelError,
			compactViewProjectionPtr)) {
		    if (submitAction->viewState->removeCadPayload(payload))
			submitAction->updatedCutCount++;
		    continue;
		}
		RetainedCandidate candidate;
		candidate.payload = payload;
		candidate.occurrenceKey = summary.sourceInstanceKey;
		candidate.projectedPixels = projected.projectedPixelDiameter;
		candidate.emphasis = summary.selected ? 2 :
		    (summary.highlighted ? 1 : 0);
		candidate.minimumLevel = payload->activeLevel;
		candidate.requestedLevel = projected.requestedLevel;
		candidate.minimumFaces = payload->counts.faceCount;
		if (payload->progressiveMesh) {
		    candidate.minimumLevel =
			payload->progressiveMesh->minimumLevel();
		    if (!payload->progressiveMesh->canDrawLevel(
			    candidate.minimumLevel)) {
			if (submitAction->viewState->removeCadPayload(payload))
			    submitAction->updatedCutCount++;
			continue;
		    }
		    candidate.minimumFaces =
			payload->progressiveMesh->faceCount(
			    candidate.minimumLevel);
		}
		retainedCandidates.push_back(std::move(candidate));
	    }
	    std::sort(retainedCandidates.begin(), retainedCandidates.end(),
		[](const RetainedCandidate &a, const RetainedCandidate &b) {
		if (a.emphasis != b.emphasis)
		    return a.emphasis > b.emphasis;
		if (a.projectedPixels < b.projectedPixels ||
		    a.projectedPixels > b.projectedPixels)
		    return a.projectedPixels > b.projectedPixels;
		return bu_strcmp(a.occurrenceKey.getString(),
		    b.occurrenceKey.getString()) < 0;
		});
	    for (const RetainedCandidate &candidate : retainedCandidates) {
		const size_t remaining =
		    submitAction->retainedSceneFaceBudgetUsed >=
			    submitAction->retainedSceneFaceBudget ?
			0 : submitAction->retainedSceneFaceBudget -
			    submitAction->retainedSceneFaceBudgetUsed;
		const bool exceedsBudget =
		    candidate.minimumFaces > remaining;
		if (exceedsBudget) {
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->refinementBudgetBlockedCount++;
		}
		submitAction->retainedSceneFaceBudgetUsed =
		    candidate.minimumFaces >
			    SIZE_MAX -
				submitAction->retainedSceneFaceBudgetUsed ?
			SIZE_MAX :
			submitAction->retainedSceneFaceBudgetUsed +
			    candidate.minimumFaces;
		submitAction->retainedRecoveredOccurrences.insert(
		    candidate.occurrenceKey.getString());
		if (candidate.payload->activeLevel !=
			candidate.minimumLevel &&
		    candidate.minimumFaces < candidate.payload->counts.faceCount &&
		    submitAction->viewState->retargetCadPayload(
			candidate.payload, candidate.minimumLevel,
			candidate.requestedLevel,
			submitAction->viewRevision,
			submitAction->policyRevision))
		    submitAction->updatedCutCount++;
		if (candidate.minimumLevel < candidate.requestedLevel)
		    submitAction->pendingRetainedRefinementCount++;
	    }
	}
	struct CompactCandidate {
	    size_t index;
	    float projectedPixels;
	    double refinementValuePerFace;
	    int emphasis;
	    bool needsCoverage;
	};
	if (!submitAction->compactEntryPlanSupplied) {
	    std::vector<CompactCandidate> candidates;
	    std::vector<const BObolViewLodState::CadPayload *>
		offscreenPayloads;
	    candidates.reserve(count > 0 ? static_cast<size_t>(count) : 0);
	    for (int candidateIndex = 0; candidateIndex < count;
		candidateIndex++) {
		BObolCompactLodPlanningSummary summary;
		if (!source->getCompactLodPlanningSummary(
			candidateIndex, summary) ||
		    !summary.valid || !summary.visible)
		    continue;
		const bool lodEligible =
		    summary.lodBacked && summary.sourceMeshRequestValid;
		const bool meshEligible = lodEligible ||
		    (!submitAction->requireLodBacked &&
		     (summary.meshGeometry ||
		      (source->realizationRoleFlags.getValue() &
		       SoBRLDatabaseSource::REALIZATION_ROLE_MESH)));
		if (!meshEligible) {
		    if (submitAction->compactEntryFirst != 0)
			continue;
		    submitAction->visitedMeshCount++;
		    submitAction->skippedMeshCount++;
		    continue;
		}

		BObolLodRequest projected;
		projected.bounds = summary.localBounds;
		if (submitAction->useViewVolume &&
		    !mesh_lod_apply_projected_demand(projected,
			summary.localBounds, summary.localToSource,
			submitAction->viewVolume, submitAction->view,
			submitAction->targetPixelError,
			compactViewProjectionPtr)) {
		    const BObolViewLodState::CadPayload *active =
			submitAction->viewState ?
			submitAction->viewState->findCadForOccurrence(
			    source, summary.sourceInstanceKey) : NULL;
		    if (active)
			offscreenPayloads.push_back(active);
		    if (submitAction->compactEntryFirst == 0) {
			submitAction->visitedMeshCount++;
			submitAction->skippedMeshCount++;
		    }
		    continue;
		}
		CompactCandidate candidate;
		candidate.index = static_cast<size_t>(candidateIndex);
		candidate.projectedPixels =
		    projected.projectedPixelDiameter;
		candidate.emphasis = summary.selected ? 2 :
		    (summary.highlighted ? 1 : 0);
		const BObolViewLodState::CadPayload *active =
		    submitAction->viewState ?
		    submitAction->viewState->findCadForOccurrence(
			source, summary.sourceInstanceKey) : NULL;
		candidate.needsCoverage = active == NULL;
		candidate.refinementValuePerFace =
		    static_cast<double>(candidate.projectedPixels);
		/* Rank already-covered leaves by visible PoP error reduction per
		 * additional face.  Raw object size alone lets large, already
		 * smooth skins consume the budget while a cheaper silhouette
		 * feature (such as an aircraft tail) remains spectacularly
		 * coarse.  PoP's quantization bound makes the screen-space error
		 * reduction available analytically; the only per-frame work is
		 * this constant-time score and the sort we already perform.
		 *
		 * The affected-screen-span factor approximates silhouette/area
		 * impact without scanning triangles.  Exact per-level geometric
		 * delta metrics may later be cached by the PoP builder, but are
		 * not required in the hot view path. */
		if (active && active->progressiveMesh &&
		    active->activeLevel >= 0 &&
		    projected.requestedLevel > active->activeLevel) {
		    const int nextLevel = std::min(
			projected.requestedLevel, active->activeLevel + 1);
		    const size_t activeFaces =
			active->progressiveMesh->hierarchyFaceCount(
			    active->activeLevel);
		    const size_t nextFaces =
			active->progressiveMesh->hierarchyFaceCount(nextLevel);
		    const size_t addedFaces = nextFaces > activeFaces ?
			nextFaces - activeFaces : 1;
		    const double currentError = std::ldexp(
			static_cast<double>(candidate.projectedPixels),
			-active->activeLevel);
		    const double nextError = std::ldexp(
			static_cast<double>(candidate.projectedPixels),
			-nextLevel);
		    const double visibleBenefit =
			std::max(1.0,
			    static_cast<double>(candidate.projectedPixels)) *
			std::max(0.0, currentError - nextError);
		    candidate.refinementValuePerFace =
			visibleBenefit / static_cast<double>(addedFaces);
		}
		candidates.push_back(std::move(candidate));
	    }
	    /* Scene coverage is the first progressive objective: give visible
	     * boxes their first useful mesh before enriching already-present
	     * leaves.  Within refinement, spend faces where the next PoP prefix
	     * removes the most estimated screen-space error.  Selection and
	     * highlight remain absolute user-intent priorities. */
	    std::sort(candidates.begin(), candidates.end(),
		[](const CompactCandidate &a, const CompactCandidate &b) {
		if (a.emphasis != b.emphasis)
		    return a.emphasis > b.emphasis;
		if (a.needsCoverage != b.needsCoverage)
		    return a.needsCoverage;
		if (a.refinementValuePerFace <
			b.refinementValuePerFace ||
		    a.refinementValuePerFace >
			b.refinementValuePerFace)
		    return a.refinementValuePerFace >
			b.refinementValuePerFace;
		if (a.projectedPixels < b.projectedPixels ||
		    a.projectedPixels > b.projectedPixels)
		    return a.projectedPixels > b.projectedPixels;
		return a.index < b.index;
		});
	    submitAction->compactEntryPlan.clear();
	    submitAction->compactEntryPlan.reserve(candidates.size());
	    for (const CompactCandidate &candidate : candidates)
		submitAction->compactEntryPlan.push_back(candidate.index);
	    for (const BObolViewLodState::CadPayload *payload :
		offscreenPayloads) {
		if (submitAction->viewState->removeCadPayload(payload))
		    submitAction->updatedCutCount++;
	    }
	}
	const size_t candidateCount = submitAction->compactEntryPlan.size();
	const size_t sourceFirst = submitAction->compactEntryTotal;
	submitAction->compactEntryTotal += candidateCount;
	const size_t rangeFirst = submitAction->compactEntryFirst;
	const size_t rangeLast = submitAction->compactEntryLimit == SIZE_MAX ||
	    rangeFirst > SIZE_MAX - submitAction->compactEntryLimit ? SIZE_MAX :
	    rangeFirst + submitAction->compactEntryLimit;
	const size_t localFirst = rangeFirst > sourceFirst ?
	    std::min(candidateCount, rangeFirst - sourceFirst) : 0;
	const size_t sourceLast = sourceFirst + candidateCount;
	const size_t localLast = rangeLast < sourceLast ?
	    (rangeLast > sourceFirst ? rangeLast - sourceFirst : 0) :
	    candidateCount;
	size_t processedLast = localFirst;
	for (size_t candidateOffset = localFirst;
	    candidateOffset < localLast; candidateOffset++) {
	    processedLast = candidateOffset + 1;
	    const size_t i =
		submitAction->compactEntryPlan[candidateOffset];
	    BObolCompactLodInstanceSummary summary;
	    if (!source->getCompactLodInstanceSummary(
		    static_cast<int>(i), summary))
		continue;
	    if (!summary.valid || !summary.visible)
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
	    request.objectPath = summary.meshAssetPath.getLength() > 0 ?
		summary.meshAssetPath : target;
	    request.objectName = summary.meshAssetName.getLength() > 0 ?
		summary.meshAssetName : summary.sourceName;
	    request.occurrenceKey = summary.sourceInstanceKey;
	    if (request.objectName.getLength() == 0)
		request.objectName = mesh_lod_source_leaf_name(source);
	    request.viewRevision = submitAction->viewRevision;
	    request.policyRevision = submitAction->policyRevision;
	    request.drawMode = sourceDrawMode;
	    request.providerId = submitAction->providerId;
	    request.providerVersion = submitAction->providerVersion;
	    request.qualityTier = submitAction->qualityTier;
	    request.bounds = !summary.meshAssetBounds.isEmpty() ?
		summary.meshAssetBounds : summary.localBounds;
	    request.sourceCounts.faceCount = summary.sourceFaceCount;
	    request.sourceCounts.pointCount = summary.sourcePointCount;
	    request.sourceCounts.originalPointCount = summary.sourcePointCount;
	    /* Compact summaries already report the complete geometry-to-root
	     * transform, including the source draw matrix.  Applying the source
	     * matrix again corrupts screen bounds (and can make a visible leaf
	     * look tiny or off-screen to the LoD selector). */
	    SbMatrix localToRoot = summary.localToSource;
	    if (submitAction->useViewVolume) {
		if (!mesh_lod_apply_projected_demand(
		    request, summary.localBounds, localToRoot,
		    submitAction->viewVolume, submitAction->view,
		    submitAction->targetPixelError,
		    compactViewProjectionPtr)) {
		    submitAction->skippedMeshCount++;
		    continue;
		}
	    }
	    const BObolViewLodState::CadPayload *activePayload =
		submitAction->viewState ?
		submitAction->viewState->findCadForOccurrence(
		    source, request.occurrenceKey) : NULL;
	    /* A completed xpush/PCA adoption may replace its provisional
	     * occurrence key while retaining the exact, distinct source BoT.
	     * Recover that payload only when its asset path is unique in this
	     * source.  Shared instances deliberately remain occurrence-keyed:
	     * choosing one of several identical asset paths would couple their
	     * independent view cuts. */
	    if (!activePayload && count == 1 &&
		request.objectPath.getLength() > 0) {
		activePayload = submitAction->viewState ?
		    submitAction->viewState->findCadForAsset(
			source, request.objectPath) : NULL;
	    }
	    if (submitAction->useForcedLevel)
		request.requestedLevel = submitAction->forcedLevel;
	    else if (activePayload)
		mesh_lod_apply_level_hysteresis(request,
		    activePayload->activeLevel);

	    /* A finite refinement allowance prevents new growth, but by itself
	     * cannot recover when thousands of already-resident minimum cuts
	     * collectively exceed the newly calibrated scene capacity.  During
	     * that recovery pass, retarget retained occurrences in the pinned
	     * screen-priority order.  The coherent minimum prefix is an
	     * unconditional coverage floor even when their aggregate population
	     * exceeds the learned budget. */
	    if (activePayload &&
		submitAction->retainedRecoveredOccurrences.find(
		    request.occurrenceKey.getString()) ==
		    submitAction->retainedRecoveredOccurrences.end() &&
		submitAction->retainedSceneFaceBudget != SIZE_MAX) {
		size_t admittedFaces = activePayload->counts.faceCount;
		int admittedLevel = request.requestedLevel;
		const int targetLevel = request.requestedLevel;
		if (activePayload->progressiveMesh) {
		    const int minimumLevel =
			activePayload->progressiveMesh->minimumLevel();
		    /* Recovery is explicitly coverage-first.  Allocating each
		     * early occurrence's richest affordable target can exhaust
		     * the scene budget before later leaves are visited, making
		     * loaded geometry flash back to structural boxes.  Re-admit
		     * one coherent minimum PoP mesh per occurrence first.  A
		     * following current-view pass spends the remaining allowance
		     * on marginal error reduction in the pinned priority order. */
		    admittedLevel = minimumLevel;
		    admittedFaces =
			activePayload->progressiveMesh->
			    hierarchyFaceCount(admittedLevel);
		}
		const size_t remaining =
		    submitAction->retainedSceneFaceBudgetUsed >=
			    submitAction->retainedSceneFaceBudget ?
			0 : submitAction->retainedSceneFaceBudget -
			    submitAction->retainedSceneFaceBudgetUsed;
		if (admittedFaces > remaining) {
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->refinementBudgetBlockedCount++;
		}
		submitAction->retainedSceneFaceBudgetUsed =
		    admittedFaces >
			    SIZE_MAX -
				submitAction->retainedSceneFaceBudgetUsed ?
			SIZE_MAX :
			submitAction->retainedSceneFaceBudgetUsed +
			    admittedFaces;
		request.requestedLevel = admittedLevel;
		if (admittedLevel < targetLevel)
		    submitAction->pendingRetainedRefinementCount++;
	    }
	    mesh_lod_trace_projected_request(request, summary.localBounds,
		localToRoot, activePayload ? activePayload->activeLevel : -1);
	    const std::string geometryAsset =
		request.objectName.getLength() > 0 ?
		    std::string("name:") + request.objectName.getString() :
		    std::string("path:") + request.objectPath.getString();
	    auto cachedGeometry =
		geometryKeysByAsset.find(geometryAsset);
	    if (cachedGeometry == geometryKeysByAsset.end()) {
		const BObolLodCacheKey generated =
		    bobol_lod_geometry_cache_key(request);
		cachedGeometry = geometryKeysByAsset.emplace(
		    geometryAsset, generated.value).first;
	    }
	    const SbString &requestGeometryKey = cachedGeometry->second;

	    /* Repeated CAD instances share one retained PoP asset but retain
	     * independent view cuts.  Once any occurrence has loaded that asset,
	     * bind another occurrence directly from the resident prefix instead
	     * of scheduling a cache/provider task and result for every instance.
	     * The exact per-occurrence draw population is still charged to the
	     * aggregate scene budget. */
	    if (!activePayload && request.objectPath.getLength() > 0 &&
		submitAction->viewState) {
		const BObolViewLodState::CadPayload *reusable =
		    submitAction->viewState->findCadForAsset(
			source, request.objectPath);
		if (reusable && reusable->progressiveMesh) {
		    int reusedLevel = std::min(request.requestedLevel,
			reusable->progressiveMesh->maximumLevel());
		    const int minimumLevel =
			reusable->progressiveMesh->minimumLevel();
		    /* Retained-admission recovery is coverage-first.  A newly
		     * visible occurrence may share a fully resident rich asset,
		     * but binding that rich prefix before the other visible
		     * occurrences receive their minimum meshes lets one instance
		     * consume the complete calibrated scene budget.  This was
		     * observable as an orthographic view round-trip changing two
		     * visible Lucys into one, despite an identical camera.
		     *
		     * Charge a minimum drawable shared prefix to the total retained
		     * admission budget here.  A later current-view epoch refines it
		     * with the ordinary marginal face allowance. */
		    if (submitAction->retainedSceneFaceBudget != SIZE_MAX)
			reusedLevel = minimumLevel;
		    while (reusedLevel >= minimumLevel &&
			!reusable->progressiveMesh->canDrawLevel(reusedLevel))
			reusedLevel--;
		    if (reusedLevel >= minimumLevel) {
			const size_t reusedFaces =
			    reusable->progressiveMesh->faceCount(reusedLevel);
			SbBool admitted = TRUE;
			if (submitAction->retainedSceneFaceBudget != SIZE_MAX) {
			    const size_t remaining =
				submitAction->retainedSceneFaceBudgetUsed >=
					submitAction->retainedSceneFaceBudget ?
				    0 :
				    submitAction->retainedSceneFaceBudget -
				    submitAction->
					    retainedSceneFaceBudgetUsed;
			    if (reusedFaces > remaining) {
				submitAction->refinementBudgetBlockedCount++;
			    }
			    submitAction->retainedSceneFaceBudgetUsed =
				reusedFaces >
					SIZE_MAX -
					    submitAction->
						retainedSceneFaceBudgetUsed ?
				    SIZE_MAX :
				    submitAction->retainedSceneFaceBudgetUsed +
					reusedFaces;
			} else {
			    admitted = submitAction->reserveInitialFaceCost(
				reusedFaces);
			}
			if (!admitted) {
			    submitAction->pendingRetainedRefinementCount++;
			    submitAction->skippedMeshCount++;
			    continue;
			}
			BObolLodResult reusedResult;
			reusedResult.generation = submitAction->generation;
			reusedResult.request = request;
			reusedResult.cacheKey =
			    bobol_lod_cache_key(request);
			reusedResult.geometry.kind =
			    BOBOL_LOD_GEOMETRY_MESH_LOD_CACHE;
			reusedResult.geometry.providerId = request.providerId;
			reusedResult.geometry.providerVersion =
			    request.providerVersion;
			reusedResult.geometry.cacheKey.value =
			    reusable->cacheKey;
			reusedResult.geometry.activeLevel = reusedLevel;
			reusedResult.geometry.borrowed = TRUE;
			reusedResult.progressiveMesh =
			    reusable->progressiveMesh;
			reusedResult.residentLevel =
			    reusable->progressiveMesh->residentLevel();
			reusedResult.resultKind = BOBOL_LOD_RESULT_MESH;
			reusedResult.qualityTier = std::max(
			    request.qualityTier, reusable->qualityTier);
			reusedResult.providerStatus =
			    BOBOL_LOD_PROVIDER_READY;
			reusedResult.bounds =
			    reusable->progressiveMesh->bounds();
			reusedResult.counts.faceCount = reusedFaces;
			reusedResult.counts.pointCount =
			    reusable->progressiveMesh->pointCount(reusedLevel);
			reusedResult.counts.originalPointCount =
			    reusedResult.counts.pointCount;
			reusedResult.counts.normalCount =
			    reusable->hasNormals ? reusedFaces * 3 : 0;
			reusedResult.hasSnappedPoints = FALSE;
			reusedResult.hasNormals = reusable->hasNormals;
			reusedResult.shadedCullBackfaces =
			    reusable->shadedCullBackfaces;
			reusedResult.terminal =
			    reusedLevel >= request.requestedLevel ?
			    TRUE : FALSE;
			if (submitAction->viewState->applySourceResult(
				source, reusedResult)) {
			    submitAction->updatedCutCount++;
			    activePayload =
				submitAction->viewState->findCadForResult(
				    source, reusedResult);
			    if (!reusedResult.terminal)
				submitAction->
				    pendingRetainedRefinementCount++;
			    continue;
			}
		    }
		}
	    }

	    /* A camera epoch is not a geometry invalidation.  Keep the current
	     * cut when it already satisfies this demand.  The controller may
	     * explicitly allow a cheaper retained draw cut during motion; that
	     * changes only prefix selection and never rereads/rebuilds geometry. */
	    if (!submitAction->useForcedLevel &&
		mesh_lod_cad_payload_is_resident(activePayload, request,
		    &requestGeometryKey)) {
		submitAction->skippedMeshCount++;
		submitAction->appendDiagnostic(target,
		    "current LoD request is already resident");
		continue;
	    }
	    if (!submitAction->useForcedLevel &&
		!submitAction->allowLevelDowngrade && activePayload &&
		activePayload->activeLevel > request.requestedLevel &&
		mesh_lod_cad_payload_has_same_source(activePayload, request,
		    &requestGeometryKey)) {
		submitAction->skippedMeshCount++;
		submitAction->appendDiagnostic(target,
		    "finer resident LoD retained during interaction");
		continue;
	    }
	    if (!submitAction->useForcedLevel &&
		!submitAction->allowRetainedRefinement && activePayload &&
		request.requestedLevel > activePayload->activeLevel &&
		mesh_lod_geometry_key_matches(activePayload->cacheKey, request,
		    &requestGeometryKey)) {
		(void)submitAction->viewState->retargetCadPayload(activePayload,
		    activePayload->activeLevel, request.requestedLevel,
		    request.viewRevision, request.policyRevision);
		submitAction->skippedMeshCount++;
		continue;
	    }
	    const SbBool retainedTargetDrawable =
		!submitAction->useForcedLevel && activePayload &&
		activePayload->progressiveMesh &&
		request.requestedLevel > activePayload->activeLevel &&
		activePayload->progressiveMesh->canDrawLevel(
		    request.requestedLevel);
	    const int retargetLevel =
		mesh_lod_available_cad_retarget_level(activePayload, request,
		    TRUE);
	    if (submitAction->reset == 0 && activePayload &&
		mesh_lod_geometry_key_matches(activePayload->cacheKey, request,
		    &requestGeometryKey) &&
		retargetLevel > activePayload->activeLevel &&
		!submitAction->reserveRefinementFaces(
		    activePayload->progressiveMesh, activePayload->activeLevel,
		    retargetLevel)) {
		(void)submitAction->viewState->retargetCadPayload(activePayload,
		    activePayload->activeLevel, request.requestedLevel,
		    request.viewRevision, request.policyRevision);
		submitAction->pendingRetainedRefinementCount++;
		submitAction->skippedMeshCount++;
		continue;
	    }
	    if (submitAction->reset == 0 && activePayload &&
		mesh_lod_geometry_key_matches(activePayload->cacheKey, request,
		    &requestGeometryKey) &&
		retargetLevel >= 0 &&
		activePayload->activeLevel != retargetLevel &&
		submitAction->viewState->retargetCadPayload(activePayload,
		    retargetLevel, request.requestedLevel, request.viewRevision,
		    request.policyRevision)) {
		const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (filter && filter[0] &&
		    (strstr(request.objectName.getString(), filter) ||
		     strstr(request.objectPath.getString(), filter)))
		    bu_log("BObol LoD submit trace object=%s level=%d "
			   "retargeted_resident=%d\n",
			   request.objectName.getString(),
			   retargetLevel,
			   activePayload->progressiveMesh ?
			       activePayload->progressiveMesh->residentLevel() : -1);
		submitAction->updatedCutCount++;
		if (retainedTargetDrawable &&
		    retargetLevel < request.requestedLevel) {
		    submitAction->pendingRetainedRefinementCount++;
		    continue;
		}
		/* If the full view target is already drawable, this is a pure cut
		 * change and no provider work is needed.  Otherwise keep the richer
		 * resident cut visible and continue below to request only the
		 * still-missing suffix. */
		if (retargetLevel == request.requestedLevel)
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
	    if (!submitAction->useForcedLevel && activePayload &&
		activePayload->progressiveMesh &&
		request.requestedLevel > activePayload->activeLevel &&
		mesh_lod_geometry_key_matches(activePayload->cacheKey, request,
		    &requestGeometryKey)) {
		const int nextLevel = std::min(request.requestedLevel,
		    activePayload->activeLevel + 1);
		if (!submitAction->reserveRefinementFaces(
			activePayload->progressiveMesh,
			activePayload->activeLevel, nextLevel)) {
		    (void)submitAction->viewState->retargetCadPayload(
			activePayload, activePayload->activeLevel,
			request.requestedLevel, request.viewRevision,
			request.policyRevision);
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->skippedMeshCount++;
		    continue;
		}
	    }
	    if (getenv("BOBOL_LOD_TRACE_BUDGET"))
		bu_log("BObol LoD provider admission object=%s occurrence=%s "
		       "active=%d requested=%d key_match=%d budget=%zu\n",
		       request.objectName.getString(),
		       request.occurrenceKey.getString(),
		       activePayload ? activePayload->activeLevel : -1,
		       request.requestedLevel,
		       activePayload && mesh_lod_geometry_key_matches(
			   activePayload->cacheKey, request,
			   &requestGeometryKey) ? 1 : 0,
		       submitAction->refinementFaceBudget);

	    if (submitAction->submittedTaskCount >=
		    submitAction->submissionTaskLimit) {
		/* Do not consume this cursor entry: it still needs provider
		 * work.  The controller will resume here after result capacity
		 * becomes available.  Cheap resident/shared bindings before
		 * this point have already been applied in the same traversal. */
		processedLast = candidateOffset;
		break;
	    }

	    size_t initialFaceAllowance = 0;
	    if (!activePayload &&
		!submitAction->reserveInitialFaces(
		    request.sourceCounts.faceCount, initialFaceAllowance)) {
		submitAction->pendingRetainedRefinementCount++;
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
	    provider->useCurrentDrawLevel = TRUE;
	    provider->currentDrawLevel =
		activePayload ? activePayload->activeLevel : -1;
	    if (initialFaceAllowance) {
		provider->initialRefinementFaceBudget =
		    initialFaceAllowance;
		provider->initialRefinementPointBudget =
		    initialFaceAllowance > SIZE_MAX / 2 ?
			SIZE_MAX : initialFaceAllowance * 2;
	    }
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
	    submitAction->compactEntryNext, sourceFirst + processedLast);
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
    if (!submitAction->useForcedLevel &&
	!submitAction->allowRetainedRefinement && activePayload &&
	request.requestedLevel > activePayload->activeLevel &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request)) {
	(void)submitAction->viewState->retargetCadPayload(activePayload,
	    activePayload->activeLevel, request.requestedLevel,
	    request.viewRevision, request.policyRevision);
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    const SbBool retainedTargetDrawable =
	!submitAction->useForcedLevel && activePayload &&
	activePayload->progressiveMesh &&
	request.requestedLevel > activePayload->activeLevel &&
	activePayload->progressiveMesh->canDrawLevel(request.requestedLevel);
    const int retargetLevel =
	mesh_lod_available_cad_retarget_level(activePayload, request,
	    retainedTargetDrawable);
    if (submitAction->reset == 0 && activePayload &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
	retargetLevel > activePayload->activeLevel &&
	!submitAction->reserveRefinementFaces(
	    activePayload->progressiveMesh, activePayload->activeLevel,
	    retargetLevel)) {
	(void)submitAction->viewState->retargetCadPayload(activePayload,
	    activePayload->activeLevel, request.requestedLevel,
	    request.viewRevision, request.policyRevision);
	submitAction->pendingRetainedRefinementCount++;
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    if (submitAction->reset == 0 && activePayload &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
	retargetLevel >= 0 &&
	activePayload->activeLevel != retargetLevel &&
	submitAction->viewState->retargetCadPayload(activePayload,
	    retargetLevel, request.requestedLevel, request.viewRevision,
	    request.policyRevision)) {
	submitAction->updatedCutCount++;
	if (retainedTargetDrawable && retargetLevel < request.requestedLevel) {
	    submitAction->pendingRetainedRefinementCount++;
	    source->doAction(action);
	    return;
	}
	if (retargetLevel == request.requestedLevel) {
	    source->doAction(action);
	    return;
	}
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
    if (!submitAction->useForcedLevel && activePayload &&
	activePayload->progressiveMesh &&
	request.requestedLevel > activePayload->activeLevel &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request)) {
	const int nextLevel = std::min(request.requestedLevel,
	    activePayload->activeLevel + 1);
	if (!submitAction->reserveRefinementFaces(
		activePayload->progressiveMesh, activePayload->activeLevel,
		nextLevel)) {
	    (void)submitAction->viewState->retargetCadPayload(activePayload,
		activePayload->activeLevel, request.requestedLevel,
		request.viewRevision, request.policyRevision);
	    submitAction->pendingRetainedRefinementCount++;
	    submitAction->skippedMeshCount++;
	    source->doAction(action);
	    return;
	}
    }

    size_t initialFaceAllowance = 0;
    if (!activePayload &&
	!submitAction->reserveInitialFaces(
	    request.sourceCounts.faceCount, initialFaceAllowance)) {
	submitAction->pendingRetainedRefinementCount++;
	submitAction->skippedMeshCount++;
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
    provider->useCurrentDrawLevel = TRUE;
    provider->currentDrawLevel =
	activePayload ? activePayload->activeLevel : -1;
    if (initialFaceAllowance) {
	provider->initialRefinementFaceBudget = initialFaceAllowance;
	provider->initialRefinementPointBudget =
	    initialFaceAllowance > SIZE_MAX / 2 ?
		SIZE_MAX : initialFaceAllowance * 2;
    }
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
    mesh_lod_trace_projected_request(request, request.bounds,
	SbMatrix::identity(), viewPayload ? viewPayload->activeLevel : -1);
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
    if (!submitAction->useForcedLevel &&
	!submitAction->allowRetainedRefinement && viewPayload &&
	request.requestedLevel > viewPayload->activeLevel &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request)) {
	(void)submitAction->viewState->retargetMeshPayload(viewPayload,
	    viewPayload->activeLevel, request.requestedLevel,
	    request.viewRevision, request.policyRevision);
	submitAction->skippedMeshCount++;
	return;
    }
    const SbBool retainedTargetDrawable =
	!submitAction->useForcedLevel && viewPayload &&
	viewPayload->progressiveMesh &&
	request.requestedLevel > viewPayload->activeLevel &&
	viewPayload->progressiveMesh->canDrawLevel(request.requestedLevel);
    const int retargetLevel =
	mesh_lod_available_mesh_retarget_level(viewPayload, request,
	    retainedTargetDrawable);
    if (submitAction->reset == 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request) &&
	retargetLevel > viewPayload->activeLevel &&
	!submitAction->reserveRefinementFaces(
	    viewPayload->progressiveMesh, viewPayload->activeLevel,
	    retargetLevel)) {
	(void)submitAction->viewState->retargetMeshPayload(viewPayload,
	    viewPayload->activeLevel, request.requestedLevel,
	    request.viewRevision, request.policyRevision);
	submitAction->pendingRetainedRefinementCount++;
	submitAction->skippedMeshCount++;
	return;
    }
    if (submitAction->reset == 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request) &&
	viewPayload->activeLevel != retargetLevel &&
	retargetLevel >= 0 &&
	submitAction->viewState->retargetMeshPayload(viewPayload,
	    retargetLevel, request.requestedLevel, request.viewRevision,
	    request.policyRevision)) {
	submitAction->updatedCutCount++;
	if (retainedTargetDrawable && retargetLevel < request.requestedLevel) {
	    submitAction->pendingRetainedRefinementCount++;
	    return;
	}
	if (retargetLevel == request.requestedLevel)
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
    if (!submitAction->useForcedLevel && viewPayload &&
	viewPayload->progressiveMesh &&
	request.requestedLevel > viewPayload->activeLevel &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request)) {
	const int nextLevel = std::min(request.requestedLevel,
	    viewPayload->activeLevel + 1);
	if (!submitAction->reserveRefinementFaces(
		viewPayload->progressiveMesh, viewPayload->activeLevel,
		nextLevel)) {
	    (void)submitAction->viewState->retargetMeshPayload(viewPayload,
		viewPayload->activeLevel, request.requestedLevel,
		request.viewRevision, request.policyRevision);
	    submitAction->pendingRetainedRefinementCount++;
	    submitAction->skippedMeshCount++;
	    return;
	}
    }

    size_t initialFaceAllowance = 0;
    if (!viewPayload &&
	!submitAction->reserveInitialFaces(
	    request.sourceCounts.faceCount, initialFaceAllowance)) {
	submitAction->pendingRetainedRefinementCount++;
	submitAction->skippedMeshCount++;
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
    provider->useCurrentDrawLevel = TRUE;
    provider->currentDrawLevel =
	viewPayload ? viewPayload->activeLevel : -1;
    if (initialFaceAllowance) {
	provider->initialRefinementFaceBudget = initialFaceAllowance;
	provider->initialRefinementPointBudget =
	    initialFaceAllowance > SIZE_MAX / 2 ?
		SIZE_MAX : initialFaceAllowance * 2;
    }
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
    /* Diagnostics are evidence, not a per-leaf event log.  A resident
     * 100,000-leaf scene used to allocate and copy megabytes of identical
     * text on every camera epoch. */
    static const unsigned int diagnosticLimit = 32;
    if (this->diagnosticCount >= diagnosticLimit) {
	this->suppressedDiagnosticCount++;
	return;
    }
    this->diagnosticCount++;
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
