/*        M E S H _ L O D _ S U B M I T _ A C T I O N . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/hash.h"
#include "bu/str.h"
#include "bu/datetime.h"

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
#include <limits>
#include <stdlib.h>
#include <string.h>

SO_ACTION_SOURCE(SoBRLMeshLodSubmitAction);

static size_t
mesh_lod_saturating_scaled_add(size_t total, uint64_t count, size_t scale)
{
    if (total == SIZE_MAX || !count || !scale)
	return total;
    if (count > static_cast<uint64_t>(SIZE_MAX / scale))
	return SIZE_MAX;
    const size_t bytes = static_cast<size_t>(count) * scale;
    return bytes > SIZE_MAX - total ? SIZE_MAX : total + bytes;
}

struct MeshLodBrepVariantDemand {
    uint64_t contentHash = 0;
    double absTol = 0.0;
    double relTol = 0.0;
    double normTol = 0.0;
    size_t estimatedWorkingSetBytes = 0;
    unsigned int band = 0;
    SbBool generate = FALSE;
    SbBool memoryLimited = FALSE;
};

static uint64_t
mesh_lod_brep_variant_hash(uint64_t canonicalHash, unsigned int band)
{
    if (!canonicalHash || !band)
	return canonicalHash;
    struct bu_data_hash_state *state = bu_data_hash_create();
    static const char contract[] = "BObol-BREP-tessellation-band-v1";
    bu_data_hash_update(state, contract, sizeof(contract));
    bu_data_hash_update(state, BOBOL_MESH_LOD_PROVIDER_VERSION,
	strlen(BOBOL_MESH_LOD_PROVIDER_VERSION));
    bu_data_hash_update(state, &canonicalHash, sizeof(canonicalHash));
    bu_data_hash_update(state, &band, sizeof(band));
    uint64_t hash = bu_data_hash_val(state);
    bu_data_hash_destroy(state);
    return hash ? hash : 1;
}

static size_t
mesh_lod_brep_variant_working_set(uint64_t canonicalFaces,
	unsigned int band)
{
    /* Halving a surface tolerance can grow a two-dimensional tessellation by
     * four.  Include detached face-set arrays, corner normals, PoP build
     * scratch, and one immutable publication.  This is deliberately
     * conservative: admission serializes a large refinement before the
     * tessellator can exhaust the host. */
    uint64_t faces = std::max<uint64_t>(1, canonicalFaces);
    for (unsigned int i = 0; i < band; ++i) {
	if (faces > UINT64_MAX / 4)
	    return SIZE_MAX;
	faces *= 4;
    }
    return mesh_lod_saturating_scaled_add(
	64ULL * 1024ULL * 1024ULL, faces, 512);
}

static MeshLodBrepVariantDemand
mesh_lod_brep_variant_demand(
    const BObolSourceMeshRequest &sourceRequest,
    const BObolLodRequest &projectedRequest,
    const BObolViewLodState::CadPayload *activePayload,
    const BObolLodService *service, SbBool allowGeneration)
{
    MeshLodBrepVariantDemand demand;
    demand.contentHash = sourceRequest.meshAssetContentHash;
    demand.absTol = sourceRequest.meshAssetTessellationAbsTol;
    demand.relTol = sourceRequest.meshAssetTessellationRelTol;
    demand.normTol = sourceRequest.meshAssetTessellationNormTol;
    const bool validTolerances =
	std::isfinite(demand.absTol) && demand.absTol >= 0.0 &&
	std::isfinite(demand.relTol) && demand.relTol >= 0.0 &&
	std::isfinite(demand.normTol) && demand.normTol >= 0.0 &&
	(demand.absTol > 0.0 || demand.relTol > 0.0 ||
	 demand.normTol > 0.0);
    const bool brep = BU_STR_EQUAL(
	sourceRequest.sourceType.getString(), "brep") &&
	demand.contentHash && validTolerances;
    if (!brep)
	return demand;

    /* Source tessellation is the reusable canonical representation.  It is
     * intentionally sized to remain pixel-faithful through a normal
     * full-window view; only close-ups beyond that reference generate a new
     * immutable band.  Powers of two prevent smooth wheel zoom from creating
     * a continuum of cache variants. */
    static const double canonicalPixelDiameter = 1024.0;
    const double targetError = std::max(0.125,
	static_cast<double>(projectedRequest.targetPixelError));
    const double ratio = static_cast<double>(
	projectedRequest.projectedPixelDiameter) /
	(canonicalPixelDiameter * targetError);
    unsigned int wantedBand = ratio > 1.0 ?
	static_cast<unsigned int>(std::ceil(std::log2(ratio))) : 0;
    wantedBand = std::min<unsigned int>(wantedBand, 15);

    if (!allowGeneration) {
	/* Interaction never launches a BREP tessellator.  Keep the already
	 * presented representation even while its cheap PoP prefix changes. */
	if (activePayload && activePayload->sourceContentHash)
	    demand.contentHash = activePayload->sourceContentHash;
	return demand;
    }

    unsigned int admittedBand = wantedBand;
    const size_t workingLimit = service ? service->getWorkingSetLimit() :
	SIZE_MAX;
    const size_t generationLimit = workingLimit == SIZE_MAX ? SIZE_MAX :
	workingLimit - workingLimit / 4;
    while (admittedBand > 0 &&
	mesh_lod_brep_variant_working_set(
	    sourceRequest.faceCount, admittedBand) > generationLimit)
	--admittedBand;
    demand.memoryLimited = admittedBand < wantedBand ? TRUE : FALSE;
    demand.band = admittedBand;
    demand.contentHash = mesh_lod_brep_variant_hash(
	sourceRequest.meshAssetContentHash, admittedBand);
    const double scale = std::ldexp(1.0, -static_cast<int>(admittedBand));
    demand.absTol = demand.absTol > 0.0 ? demand.absTol * scale : 0.0;
    demand.relTol *= scale;
    demand.normTol = demand.normTol > 0.0 ? demand.normTol * scale : 0.0;
    demand.estimatedWorkingSetBytes =
	mesh_lod_brep_variant_working_set(
	    sourceRequest.faceCount, admittedBand);
    demand.generate = admittedBand > 0 ? TRUE : FALSE;
    return demand;
}

/*
 * A retained/warm task does not build PoP topology.  It reads one cumulative
 * cache prefix, converts it to an immutable float/index generation, and
 * prepares the requested renderer channels.  Reserve that bounded population
 * instead of charging the cold full-source topology estimate.
 *
 * The coefficients include the cache's double/int prefix, its immutable Obol
 * conversion, the renderer snapshot, vector growth, and a fixed transaction /
 * allocator allowance.  Authored corner normals are intentionally expensive:
 * converting them to Obol's indexed position/normal contract may split every
 * triangle corner.  Wire geometry expands every triangle to three segments.
 */
static size_t
mesh_lod_resident_working_set_estimate(uint64_t pointCount,
	uint64_t faceCount, SbBool hasNormals, int drawMode)
{
    size_t estimate = 8ULL * 1024ULL * 1024ULL;
    estimate = mesh_lod_saturating_scaled_add(estimate, pointCount, 64);
    size_t faceBytes = 48;
    if (hasNormals)
	faceBytes += 320;
    if (drawMode == BOBOL_LOD_DRAW_WIRE ||
	drawMode == BOBOL_LOD_DRAW_HIDDEN_LINE)
	faceBytes += 84;
    return mesh_lod_saturating_scaled_add(estimate, faceCount, faceBytes);
}

static size_t
mesh_lod_resident_task_estimate(
    const BObolLodProgressiveMeshPtr &mesh, int targetCut,
    SbBool hasNormals, int drawMode)
{
    if (!mesh || targetCut < mesh->minimumCut())
	return 0;
    targetCut = std::min(targetCut, mesh->maximumCut());
    return mesh_lod_resident_working_set_estimate(
	mesh->hierarchyPointCountAtCut(targetCut),
	mesh->hierarchyFaceCountAtCut(targetCut), hasNormals, drawMode);
}

static size_t
mesh_lod_warm_source_task_estimate(
    const BObolSourceMeshRequest &sourceRequest,
    const BObolLodRequest &request)
{
    if (!sourceRequest.lodAvailable ||
	sourceRequest.lodActiveCut < 0 ||
	!sourceRequest.lodFaceCount || !sourceRequest.lodPointCount)
	return 0;

    /*
     * The foreground source loaded this cache cut for the same view.  If the
     * request asks beyond it and there is not yet a retained hierarchy, the
     * target population is unknown; fall back to the cold-safe source-count
     * estimator.  Once the first generation is retained, exact hierarchy
     * counts drive every subsequent reservation.
     */
    if (request.requestedCut > sourceRequest.lodActiveCut)
	return 0;
    return mesh_lod_resident_working_set_estimate(
	sourceRequest.lodPointCount, sourceRequest.lodFaceCount,
	sourceRequest.lodHasNormals ? TRUE : FALSE, request.drawMode);
}

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
    useForcedCut(FALSE),
    forcedCut(0),
    requireLodBacked(TRUE),
    allowCutDowngrade(FALSE),
    allowRetainedRefinement(TRUE),
    allowResidentPrefetch(FALSE),
    refinementCutCeiling(-1),
    allowRepresentationRefinement(TRUE),
    preserveMeshCoverage(FALSE),
    refinementCostBudget(SIZE_MAX),
    refinementCostBudgetUsed(0),
    refinementBudgetBlockedCount(0),
    retainedQualityLimitedCount(0),
    retainedAdmissionBlockedCount(0),
    oneOverBudgetRefinementLimit(0),
    oneOverBudgetRefinementUsed(FALSE),
    transitionLimitedRefinement(FALSE),
    viewState(NULL),
    compactEntryFirst(0),
    compactEntryLimit(SIZE_MAX),
    compactEntryNext(0),
    compactEntryTotal(0),
    compactEntryPlan(),
    compactEntryPlanView(NULL),
    compactEntryPlanSupplied(FALSE),
    submissionTaskLimit(SIZE_MAX),
    submissionTimeLimitMicroseconds(0),
    retainedSceneCostBudget(SIZE_MAX),
    retainedSceneCostBudgetUsed(0),
    retainedSceneUpgradeCostBudget(SIZE_MAX),
    retainedSceneUpgradeCostBudgetUsed(0),
    retainedSceneMaximumNormalizedError(
	std::numeric_limits<double>::infinity()),
    retainedRecoveredOccurrences(),
    deferredCompactEntries(FALSE),
    visitedMeshCount(0),
    submittedTaskCount(0),
    updatedCutCount(0),
    pendingRetainedRefinementCount(0),
    pendingResidentRefinementCount(0),
    skippedMeshCount(0),
    visibleMeshCount(0),
    coveredVisibleMeshCount(0),
    diagnosticCount(0),
    suppressedDiagnosticCount(0),
    diagnostics("")
{
    SO_ACTION_CONSTRUCTOR(SoBRLMeshLodSubmitAction);
    bv_view_info_init(&this->view);
}

void
SoBRLMeshLodSubmitAction::setAllowRepresentationRefinement(SbBool allow)
{
    this->allowRepresentationRefinement = allow;
}

SbBool
SoBRLMeshLodSubmitAction::getAllowRepresentationRefinement(void) const
{
    return this->allowRepresentationRefinement;
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
SoBRLMeshLodSubmitAction::setForcedCut(int cut)
{
    this->useForcedCut = TRUE;
    this->forcedCut = cut < 0 ? 0 : cut;
}

void
SoBRLMeshLodSubmitAction::clearForcedCut(void)
{
    this->useForcedCut = FALSE;
    this->forcedCut = 0;
}

SbBool
SoBRLMeshLodSubmitAction::hasForcedCut(void) const
{
    return this->useForcedCut;
}

int
SoBRLMeshLodSubmitAction::getForcedCut(void) const
{
    return this->forcedCut;
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
SoBRLMeshLodSubmitAction::setAllowCutDowngrade(SbBool allow)
{
    this->allowCutDowngrade = allow ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getAllowCutDowngrade(void) const
{
    return this->allowCutDowngrade;
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
SoBRLMeshLodSubmitAction::setAllowResidentPrefetch(SbBool allow)
{
    this->allowResidentPrefetch = allow ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getAllowResidentPrefetch(void) const
{
    return this->allowResidentPrefetch;
}

void
SoBRLMeshLodSubmitAction::setRefinementCutCeiling(int cut)
{
    this->refinementCutCeiling = cut;
}

int
SoBRLMeshLodSubmitAction::getRefinementCutCeiling(void) const
{
    return this->refinementCutCeiling;
}

void
SoBRLMeshLodSubmitAction::setPreserveMeshCoverage(SbBool preserve)
{
    this->preserveMeshCoverage = preserve ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getPreserveMeshCoverage(void) const
{
    return this->preserveMeshCoverage;
}

void
SoBRLMeshLodSubmitAction::setRefinementCostBudget(size_t additionalCost)
{
    this->refinementCostBudget = additionalCost;
}

size_t
SoBRLMeshLodSubmitAction::getRefinementCostBudget(void) const
{
    return this->refinementCostBudget;
}

size_t
SoBRLMeshLodSubmitAction::getRefinementCostBudgetUsed(void) const
{
    return this->refinementCostBudgetUsed;
}

unsigned int
SoBRLMeshLodSubmitAction::getRefinementBudgetBlockedCount(void) const
{
    return this->refinementBudgetBlockedCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getRetainedQualityLimitedCount(void) const
{
    return this->retainedQualityLimitedCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getRetainedAdmissionBlockedCount(void) const
{
    return this->retainedAdmissionBlockedCount;
}

void
SoBRLMeshLodSubmitAction::setOneOverBudgetRefinementLimit(size_t excessCost)
{
    this->oneOverBudgetRefinementLimit = excessCost;
}

size_t
SoBRLMeshLodSubmitAction::getOneOverBudgetRefinementLimit(void) const
{
    return this->oneOverBudgetRefinementLimit;
}

SbBool
SoBRLMeshLodSubmitAction::getOneOverBudgetRefinementUsed(void) const
{
    return this->oneOverBudgetRefinementUsed;
}

void
SoBRLMeshLodSubmitAction::setTransitionLimitedRefinement(SbBool limited)
{
    this->transitionLimitedRefinement = limited ? TRUE : FALSE;
}

SbBool
SoBRLMeshLodSubmitAction::getTransitionLimitedRefinement(void) const
{
    return this->transitionLimitedRefinement;
}

/* Return the next transition which changes the cumulative face population.
 * Quantization-only cuts preceding that transition are free and are folded
 * into it.  If no later population exists, the richest requested snap cut
 * is itself the useful (zero-face-cost) transition. */
static int
mesh_lod_next_population_cut(
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int activeCut, int preferredCut)
{
    if (!progressiveMesh || preferredCut <= activeCut)
	return activeCut;
    preferredCut = std::min(preferredCut,
	progressiveMesh->maximumCut());
    const size_t activeFaces = activeCut >= 0 ?
	progressiveMesh->hierarchyFaceCountAtCut(activeCut) : 0;
    for (int cut = activeCut + 1; cut <= preferredCut; ++cut)
	if (progressiveMesh->hierarchyFaceCountAtCut(cut) > activeFaces)
	    return cut;
    return preferredCut;
}

SbBool
SoBRLMeshLodSubmitAction::reserveRefinementCost(
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int activeCut, int nextCut, int drawMode, SbBool hasNormals)
{
    if (!progressiveMesh || nextCut <= activeCut)
	return TRUE;

    const BObolLodCounts activeCounts = activeCut >= 0 ?
	bobol_lod_progressive_counts(
	    progressiveMesh, activeCut, hasNormals) : BObolLodCounts();
    const BObolLodCounts nextCounts = bobol_lod_progressive_counts(
	progressiveMesh, nextCut, hasNormals);
    const size_t activeCost = bobol_lod_render_cost_units(
	activeCounts, drawMode, 0);
    const size_t nextCost = bobol_lod_render_cost_units(
	nextCounts, drawMode, 0);
    const size_t additionalCost = nextCost > activeCost ?
	nextCost - activeCost : 0;
    if (!additionalCost)
	return TRUE;

    if (this->refinementCostBudget != SIZE_MAX &&
	(this->refinementCostBudgetUsed > this->refinementCostBudget ||
	 additionalCost >
	    this->refinementCostBudget - this->refinementCostBudgetUsed)) {
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD refinement cost reservation blocked "
		   "active_cut=%d next_cut=%d active_faces=%zu "
		   "next_faces=%zu added_cost=%zu used_cost=%zu "
		   "budget_cost=%zu\n",
		   activeCut, nextCut,
		   static_cast<size_t>(activeCounts.faceCount),
		   static_cast<size_t>(nextCounts.faceCount),
		   additionalCost, this->refinementCostBudgetUsed,
		   this->refinementCostBudget);
	this->refinementBudgetBlockedCount++;
	return FALSE;
    }

    this->refinementCostBudgetUsed =
	additionalCost > SIZE_MAX - this->refinementCostBudgetUsed ?
	SIZE_MAX : this->refinementCostBudgetUsed + additionalCost;
    return TRUE;
}

int
SoBRLMeshLodSubmitAction::reserveRefinementCut(
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int activeCut, int preferredCut, int drawMode, SbBool hasNormals)
{
    if (preferredCut <= activeCut)
	return activeCut;

    /* The renderer-wide ceiling is reversible FPS state, not permission for
     * the occurrence allocator to advance invisibly behind it.  Keep an
     * already richer retained cut intact, but admit no new cut above the
     * presentation frontier. */
    if (this->refinementCutCeiling >= 0) {
	if (activeCut >= this->refinementCutCeiling)
	    return activeCut;
	preferredCut = std::min(
	    preferredCut, this->refinementCutCeiling);
	if (preferredCut <= activeCut)
	    return activeCut;
    }

    /*
     * A PoP cut is not necessarily a useful scheduling quantum: adjacent cuts
     * may add no population at all, while another step may add millions of
     * faces.  Select the richest cumulative prefix whose complete marginal
     * cost fits the remaining measured scene allowance.  This preserves the
     * frame budget but avoids manufacturing one provider task, publication,
     * and whole-scene rediscovery pass per cut.
     */
    if (!progressiveMesh)
	return activeCut;
    if (this->transitionLimitedRefinement &&
	this->refinementCostBudget != SIZE_MAX)
	preferredCut = mesh_lod_next_population_cut(
	    progressiveMesh, activeCut, preferredCut);
    const BObolLodCounts activeCounts = activeCut >= 0 ?
	bobol_lod_progressive_counts(
	    progressiveMesh, activeCut, hasNormals) : BObolLodCounts();
    const size_t activeCost = bobol_lod_render_cost_units(
	activeCounts, drawMode, 0);
    if (this->refinementCostBudget == SIZE_MAX) {
	(void)this->reserveRefinementCost(
	    progressiveMesh, activeCut, preferredCut, drawMode,
	    hasNormals);
	return preferredCut;
    }
    const size_t remaining =
	this->refinementCostBudgetUsed >= this->refinementCostBudget ?
	    0 : this->refinementCostBudget -
		this->refinementCostBudgetUsed;
    for (int cut = preferredCut; cut > activeCut; --cut) {
	const BObolLodCounts cutCounts = bobol_lod_progressive_counts(
	    progressiveMesh, cut, hasNormals);
	const size_t cutCost = bobol_lod_render_cost_units(
	    cutCounts, drawMode, 0);
	const size_t additionalCost = cutCost > activeCost ?
	    cutCost - activeCost : 0;
	if (additionalCost > remaining)
	    continue;
	this->refinementCostBudgetUsed =
	    additionalCost >
		    SIZE_MAX - this->refinementCostBudgetUsed ?
		SIZE_MAX :
		this->refinementCostBudgetUsed + additionalCost;
	return cut;
    }

    /* The ordinary predictor cannot learn a discontinuous PoP transition it
     * never renders.  A controller may grant one bounded scene-global trial:
     * select only the next populated prefix, consume all ordinary remainder,
     * and admit no more than the explicitly bounded excess.  Keeping the
     * complete marginal charge in refinementCostBudgetUsed makes subsequent
     * scene-wide accounting conservative even though it exceeds the nominal
     * action budget. */
    if (this->oneOverBudgetRefinementLimit > 0 &&
	!this->oneOverBudgetRefinementUsed) {
	const int trialCut = mesh_lod_next_population_cut(
	    progressiveMesh, activeCut, preferredCut);
	if (trialCut > activeCut) {
	    const BObolLodCounts trialCounts = bobol_lod_progressive_counts(
		progressiveMesh, trialCut, hasNormals);
	    const size_t trialCost = bobol_lod_render_cost_units(
		trialCounts, drawMode, 0);
	    const size_t additionalCost = trialCost > activeCost ?
		trialCost - activeCost : 0;
	    const size_t excessCost = additionalCost > remaining ?
		additionalCost - remaining : 0;
	    if (additionalCost && excessCost > 0 &&
		excessCost <= this->oneOverBudgetRefinementLimit) {
		this->refinementCostBudgetUsed =
		    additionalCost >
			    SIZE_MAX - this->refinementCostBudgetUsed ?
			SIZE_MAX :
			this->refinementCostBudgetUsed + additionalCost;
		this->oneOverBudgetRefinementUsed = TRUE;
		return trialCut;
	    }
	}
    }
    this->refinementBudgetBlockedCount++;
    return activeCut;
}

SbBool
SoBRLMeshLodSubmitAction::reserveInitialCost(
    uint64_t sourceFaces, uint64_t sourcePoints, int drawMode,
    size_t &providerFaceAllowance)
{
    providerFaceAllowance = 0;
    if (this->refinementCostBudget == SIZE_MAX)
	return TRUE;

    /* The exact minimum populated prefix is not known until the provider has
     * read the PoP hierarchy.  Reserve a conservative first-cut estimate here
     * so a warm many-leaf scene can admit a useful wave without either
     * scheduling every mesh at once or walking only a handful of leaves per
     * frame.  The result's exact face count becomes the next frame's measured
     * active population. */
    static const size_t provisionalFirstCutFaces = 512;
    const size_t knownFaces =
	sourceFaces > static_cast<uint64_t>(SIZE_MAX) ? SIZE_MAX :
	static_cast<size_t>(sourceFaces);
    const size_t reserveFaces = knownFaces ?
	std::min(knownFaces, provisionalFirstCutFaces) :
	provisionalFirstCutFaces;
    const size_t knownPoints = sourcePoints > static_cast<uint64_t>(SIZE_MAX) ?
	SIZE_MAX : static_cast<size_t>(sourcePoints);
    BObolLodCounts provisionalCounts;
    provisionalCounts.faceCount = reserveFaces;
    provisionalCounts.pointCount = knownPoints ?
	std::min(knownPoints, reserveFaces * 2) : reserveFaces;
    const size_t reserveCost = bobol_lod_render_cost_units(
	provisionalCounts, drawMode, 1);

    /*
     * Structural leaf proxies are the coverage guarantee.  A minimum PoP
     * mesh is richer data and must obey the same aggregate scene budget as
     * every later prefix.  Exempting first meshes admitted one minimum cut
     * for every visible leaf: thousands of modest meshes then exceeded the
     * calibrated population by orders of magnitude and caused a single
     * multi-hundred-millisecond upload frame.
     *
     * Mark the blocked wave so the controller presents and calibrates what it
     * admitted, then resumes with a larger measured allowance.  A controller
     * may explicitly preserve the minimum visible mesh floor; bounded entry
     * windows and service capacity still prevent an all-at-once cliff.
     */
    if (this->refinementCostBudgetUsed > this->refinementCostBudget ||
	reserveCost >
	    this->refinementCostBudget - this->refinementCostBudgetUsed) {
	this->refinementBudgetBlockedCount++;
	if (!this->preserveMeshCoverage)
	    return FALSE;
    }

    this->refinementCostBudgetUsed =
	reserveCost > SIZE_MAX - this->refinementCostBudgetUsed ?
	SIZE_MAX : this->refinementCostBudgetUsed + reserveCost;
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
SoBRLMeshLodSubmitAction::reserveInitialCost(
    const BObolLodCounts &counts, int drawMode)
{
    const size_t cost = bobol_lod_render_cost_units(counts, drawMode, 1);
    if (this->refinementCostBudget != SIZE_MAX &&
	(this->refinementCostBudgetUsed > this->refinementCostBudget ||
	 cost >
	    this->refinementCostBudget - this->refinementCostBudgetUsed)) {
	this->refinementBudgetBlockedCount++;
	if (!this->preserveMeshCoverage)
	    return FALSE;
    }
    this->refinementCostBudgetUsed =
	cost > SIZE_MAX - this->refinementCostBudgetUsed ?
	SIZE_MAX : this->refinementCostBudgetUsed + cost;
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
    this->compactEntryPlanView = NULL;
    this->compactEntryPlanSupplied = TRUE;
}

void
SoBRLMeshLodSubmitAction::setCompactEntryPlanView(
    const std::vector<size_t> *entryIndices)
{
    this->compactEntryPlanView = entryIndices;
    this->compactEntryPlanSupplied = entryIndices ? TRUE : FALSE;
}

void
SoBRLMeshLodSubmitAction::getCompactEntryPlan(
    std::vector<size_t> &entryIndices) const
{
    entryIndices = this->compactEntryPlanView ?
	*this->compactEntryPlanView : this->compactEntryPlan;
}

void
SoBRLMeshLodSubmitAction::setSubmissionTaskLimit(size_t taskCount)
{
    this->submissionTaskLimit = taskCount;
}

void
SoBRLMeshLodSubmitAction::setSubmissionTimeLimit(uint64_t microseconds)
{
    this->submissionTimeLimitMicroseconds = microseconds;
}

void
SoBRLMeshLodSubmitAction::setRetainedSceneCostBudget(size_t totalCost)
{
    this->retainedSceneCostBudget = totalCost;
}

size_t
SoBRLMeshLodSubmitAction::getRetainedSceneCostBudgetUsed(void) const
{
    return this->retainedSceneCostBudgetUsed;
}

void
SoBRLMeshLodSubmitAction::setRetainedSceneUpgradeCostBudget(
    size_t additionalCost)
{
    this->retainedSceneUpgradeCostBudget = additionalCost;
}

size_t
SoBRLMeshLodSubmitAction::getRetainedSceneUpgradeCostBudgetUsed(void) const
{
    return this->retainedSceneUpgradeCostBudgetUsed;
}

void
SoBRLMeshLodSubmitAction::setRetainedSceneMaximumNormalizedError(double error)
{
    this->retainedSceneMaximumNormalizedError =
	std::isfinite(error) && error > 0.0 ? error :
	std::numeric_limits<double>::infinity();
}

double
SoBRLMeshLodSubmitAction::getRetainedSceneMaximumNormalizedError(void) const
{
    return this->retainedSceneMaximumNormalizedError;
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
SoBRLMeshLodSubmitAction::getPendingResidentRefinementCount(void) const
{
    return this->pendingResidentRefinementCount;
}

unsigned int
SoBRLMeshLodSubmitAction::getSkippedMeshCount(void) const
{
    return this->skippedMeshCount;
}

size_t
SoBRLMeshLodSubmitAction::getVisibleMeshCount(void) const
{
    return this->visibleMeshCount;
}

size_t
SoBRLMeshLodSubmitAction::getCoveredVisibleMeshCount(void) const
{
    return this->coveredVisibleMeshCount;
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
    this->pendingResidentRefinementCount = 0;
    this->refinementCostBudgetUsed = 0;
    this->refinementBudgetBlockedCount = 0;
    this->retainedQualityLimitedCount = 0;
    this->retainedAdmissionBlockedCount = 0;
    this->oneOverBudgetRefinementUsed = FALSE;
    this->retainedSceneCostBudgetUsed = 0;
    this->retainedSceneUpgradeCostBudgetUsed = 0;
    this->retainedRecoveredOccurrences.clear();
    this->skippedMeshCount = 0;
    this->visibleMeshCount = 0;
    this->coveredVisibleMeshCount = 0;
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
struct MeshLodProjectedPoint {
    float x;
    float y;
};

static float
mesh_lod_projected_cross(const MeshLodProjectedPoint &origin,
	const MeshLodProjectedPoint &a, const MeshLodProjectedPoint &b)
{
    return (a.x - origin.x) * (b.y - origin.y) -
	(a.y - origin.y) * (b.x - origin.x);
}

/* The convex hull of eight projected box corners is a constant-time proxy for
 * occupied screen footprint.  It is substantially less misleading than the
 * area of the screen-aligned min/max rectangle, while avoiding any triangle
 * scan in the view-planning path. */
static void
mesh_lod_projected_hull_metrics(MeshLodProjectedPoint *points, int count,
	float &area, float &perimeter)
{
    area = 0.0f;
    perimeter = 0.0f;
    if (!points || count < 2)
	return;
    std::sort(points, points + count,
	[](const MeshLodProjectedPoint &a,
	   const MeshLodProjectedPoint &b) {
	    if (a.x < b.x)
		return true;
	    if (b.x < a.x)
		return false;
	    return a.y < b.y;
	});
    int uniqueCount = 0;
    for (int i = 0; i < count; ++i) {
	if (uniqueCount > 0 &&
	    !(points[i].x < points[uniqueCount - 1].x) &&
	    !(points[uniqueCount - 1].x < points[i].x) &&
	    !(points[i].y < points[uniqueCount - 1].y) &&
	    !(points[uniqueCount - 1].y < points[i].y))
	    continue;
	points[uniqueCount++] = points[i];
    }
    if (uniqueCount < 2)
	return;

    MeshLodProjectedPoint hull[16];
    int hullCount = 0;
    for (int i = 0; i < uniqueCount; ++i) {
	while (hullCount >= 2 &&
	    mesh_lod_projected_cross(hull[hullCount - 2],
		hull[hullCount - 1], points[i]) <= 0.0f)
	    --hullCount;
	hull[hullCount++] = points[i];
    }
    const int lowerCount = hullCount;
    for (int i = uniqueCount - 2; i >= 0; --i) {
	while (hullCount > lowerCount &&
	    mesh_lod_projected_cross(hull[hullCount - 2],
		hull[hullCount - 1], points[i]) <= 0.0f)
	    --hullCount;
	hull[hullCount++] = points[i];
    }
    if (hullCount > 1)
	--hullCount;
    if (hullCount < 2)
	return;

    double twiceArea = 0.0;
    double boundary = 0.0;
    for (int i = 0; i < hullCount; ++i) {
	const MeshLodProjectedPoint &a = hull[i];
	const MeshLodProjectedPoint &b = hull[(i + 1) % hullCount];
	twiceArea += static_cast<double>(a.x) * b.y -
	    static_cast<double>(a.y) * b.x;
	boundary += std::hypot(static_cast<double>(b.x - a.x),
	    static_cast<double>(b.y - a.y));
    }
    area = static_cast<float>(std::fabs(twiceArea) * 0.5);
    perimeter = static_cast<float>(boundary);
}

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
    MeshLodProjectedPoint projectedPoints[8];
    int projectedPointCount = 0;
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
	projectedPoints[projectedPointCount++] = {
	    projected[0] * static_cast<float>(std::max(1, view.width)),
	    projected[1] * static_cast<float>(std::max(1, view.height))
	};
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
    float projectedArea = 0.0f;
    float projectedPerimeter = 0.0f;
    mesh_lod_projected_hull_metrics(projectedPoints, projectedPointCount,
	projectedArea, projectedPerimeter);
    request.projectedPixelDiameter = diameter;
    request.projectedPixelArea = projectedArea;
    request.projectedPixelPerimeter = projectedPerimeter;
    request.targetPixelError = pixelError;
    /* Projection supplies physical demand, not an ordinal estimate.  An
     * existing retained hierarchy resolves it immediately in
     * mesh_lod_apply_cut_hysteresis(); a cold provider resolves the same
     * diameter/error pair from producer-certified cut metadata. */
    request.requestedCut = -1;
    return TRUE;
}

static void
mesh_lod_apply_cut_hysteresis(
    BObolLodRequest &request,
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int activeCut)
{
    if (progressiveMesh) {
	const int selectedCut = progressiveMesh->cutForScreenError(
	    request.projectedPixelDiameter, request.targetPixelError);
	if (selectedCut >= 0)
	    request.requestedCut = selectedCut;
    }
    if (activeCut < 0 || request.requestedCut < 0 ||
	request.requestedCut == activeCut || request.targetPixelError <= 0.0f)
	return;
    if (!progressiveMesh)
	return;
    const double target = static_cast<double>(request.targetPixelError);
    if (request.requestedCut > activeCut) {
	const double activeError = progressiveMesh->projectedErrorAtCut(
	    activeCut, request.projectedPixelDiameter);
	if (activeError <= target * 1.25)
	    request.requestedCut = activeCut;
    } else {
	/* cutForScreenError(target) returns the first cut at or below the
	 * ordinary quality boundary.  Testing that same cut against the stricter
	 * 0.75 release boundary and restoring activeCut on failure creates an
	 * all-or-nothing latch: the selected cut will commonly lie between 0.75
	 * and 1.0 pixels, so a much richer retained prefix can survive every
	 * zoom-out frame.  Lucy, for example, remained at cut 29 after returning
	 * to a scale which originally selected cut 20.
	 *
	 * Ask the hierarchy for the cut which actually satisfies the release
	 * boundary.  It is normally one or two cuts finer than the ordinary
	 * target but can still be substantially coarser than activeCut.  This is
	 * the Schmitt-trigger counterpart of the 1.25 refinement test above: it
	 * prevents boundary chatter without pinning all intervening detail. */
	const int releaseCut = progressiveMesh->cutForScreenError(
	    request.projectedPixelDiameter, target * 0.75);
	if (releaseCut < 0 || releaseCut >= activeCut)
	    request.requestedCut = activeCut;
	else
	    request.requestedCut = std::max(request.requestedCut, releaseCut);
    }
}

static void
mesh_lod_trace_projected_request(const BObolLodRequest &request,
	const SbBox3f &localBounds, const SbMatrix &localToRoot,
	int activeCut)
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
	   "pixels=%.9g error=%.9g request_cut=%d active_cut=%d "
	   "bounds=[%.9g %.9g %.9g]-[%.9g %.9g %.9g] "
	   "matrix=[%.9g %.9g %.9g %.9g; %.9g %.9g %.9g %.9g; "
	   "%.9g %.9g %.9g %.9g; %.9g %.9g %.9g %.9g]\n",
	   name ? name : "", path ? path : "",
	   request.occurrenceKey.getString(),
	   request.projectedPixelDiameter, request.targetPixelError,
	   request.requestedCut, activeCut,
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
	payload->activeCut == request.requestedCut;
}

/*
 * Fast authoritative identity for a view-local CAD binding.
 *
 * The old hot path regenerated the serialized persistent-cache key for every
 * already-resident occurrence on every bounded scene pass.  On a scene with
 * thousands of distinct assets that meant repeated string formatting,
 * sorting, allocation, and hashing even though the occurrence map had
 * already resolved the exact binding.  Retain the source epoch on the
 * payload and compare the small request contract directly.  Persistent cache
 * keys remain the worker/cache boundary; they are not a view-planner data
 * structure.
 */
static SbBool
mesh_lod_cad_payload_matches_asset_epoch(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request)
{
    if (!payload || !payload->isValid() ||
	(payload->resultKind != BOBOL_LOD_RESULT_MESH &&
	 payload->resultKind != BOBOL_LOD_RESULT_FULL_DETAIL) ||
	payload->providerStatus != BOBOL_LOD_PROVIDER_READY ||
	payload->databaseRevision != request.databaseRevision ||
	payload->sourceRevision != request.sourceRevision ||
	payload->sourceContentHash != request.sourceContentHash ||
	payload->qualityTier != request.qualityTier)
	return FALSE;

    /*
     * A retained PoP generation is representation-neutral: shaded and wire
     * presentations select different immutable renderer channels over the
     * same positions/indices and resident prefix.  Do not turn a mode switch
     * into a cache/provider request merely because the channel recorded on
     * the payload differs.  Non-progressive payloads remain mode-specific
     * because they may contain only the requested representation.
     */
    if (payload->drawMode != request.drawMode &&
	!payload->progressiveMesh)
	return FALSE;

    if (request.objectName.getLength() > 0)
	return payload->sourceName.getLength() > 0 &&
	    bu_strcmp(payload->sourceName.getString(),
		request.objectName.getString()) == 0 ? TRUE : FALSE;
    return payload->sourcePath.getLength() > 0 &&
	request.objectPath.getLength() > 0 &&
	bu_strcmp(payload->sourcePath.getString(),
	    request.objectPath.getString()) == 0 ? TRUE : FALSE;
}

static SbBool
mesh_lod_cad_payload_has_same_source(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request,
	const SbString *cachedGeometryKey = NULL)
{
    if (!payload || payload->requestedCut < 0)
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
	payload->activeCut == request.requestedCut;
}

static SbBool
mesh_lod_mesh_payload_has_same_source(
	const BObolViewLodState::MeshPayload *payload,
	const BObolLodRequest &request)
{
    if (!payload || payload->requestedCut < 0)
	return FALSE;
    BObolLodRequest payloadRequest = request;
    payloadRequest.requestedCut = payload->requestedCut;
    return mesh_lod_mesh_payload_is_resident(payload, payloadRequest);
}

template <typename Payload>
static SbBool
mesh_lod_payload_memory_limited_for_epoch(
    const Payload *payload,
    const BObolLodRequest &request,
    const BObolLodService *service)
{
    return payload && service && payload->memoryLimited &&
	payload->residentAdmissionRevision != 0 &&
	payload->residentAdmissionRevision ==
	    service->residentMeshAdmissionRevision() &&
	payload->viewRevision == request.viewRevision &&
	payload->policyRevision == request.policyRevision &&
	request.requestedCut > payload->activeCut ? TRUE : FALSE;
}

static void
mesh_lod_retarget_cad_demand_if_changed(
    BObolViewLodState *state,
    const BObolViewLodState::CadPayload *payload,
    const BObolLodRequest &request)
{
    if (!state || !payload)
	return;
    /* retargetCadPayload also invalidates a previous capacity-denial
     * witness.  That is correct for a new demand, but not for the identical
     * active/requested cut and epoch.  Rewriting identical metadata here
     * made a budget-blocked pass forget a current memory denial and enqueue
     * the same sparse frontier forever. */
    if (payload->requestedCut == request.requestedCut &&
	payload->drawMode == request.drawMode &&
	payload->viewRevision == request.viewRevision &&
	payload->policyRevision == request.policyRevision)
	return;
    (void)state->retargetCadPayload(
	payload, payload->activeCut, request);
}

static void
mesh_lod_retarget_mesh_demand_if_changed(
    BObolViewLodState *state,
    const BObolViewLodState::MeshPayload *payload,
    const BObolLodRequest &request)
{
    if (!state || !payload)
	return;
    if (payload->requestedCut == request.requestedCut &&
	payload->viewRevision == request.viewRevision &&
	payload->policyRevision == request.policyRevision)
	return;
    (void)state->retargetMeshPayload(payload, payload->activeCut,
	request.requestedCut, request.viewRevision.value(),
	request.policyRevision.value());
}

static void
mesh_lod_note_upward_resident_use(
    SoBRLMeshLodSubmitAction *action, const SbString &cacheKey,
    int activeCut, int nextCut)
{
    if (!action || !action->getService() || nextCut <= activeCut ||
	cacheKey.getLength() == 0)
	return;
    BObolLodCacheKey key;
    key.value = cacheKey;
    action->getService()->noteResidentMeshUse(key);
}

/* The scene allocator, rather than an arbitrary PoP-cut step, chooses the
 * exact population cut whose complete added face cost is admitted. */
static int
mesh_lod_bounded_delivery_cut(
	const BObolLodProgressiveMeshPtr &mesh,
	int activeCut, int requestedCut)
{
    (void)mesh;
    (void)activeCut;
    return requestedCut;
}

static int
mesh_lod_available_cad_retarget_cut(
	const BObolViewLodState::CadPayload *payload,
	const BObolLodRequest &request,
	SbBool incrementalRefinement)
{
    if (!payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() ||
	request.requestedCut < 0)
	return request.requestedCut;
    if (request.requestedCut <= payload->activeCut)
	return request.requestedCut;
    /* residentCut is the loaded cache prefix, but a finer PoP cut can also
     * be exactly drawable when adjacent cuts add only quantization bits and
     * no vertices/faces.  Reuse the richest cut the retained asset can
     * actually draw before scheduling more I/O. */
    for (int cut = request.requestedCut;
	 cut > payload->activeCut; --cut) {
	if (!payload->progressiveMesh->canDrawCut(cut))
	    continue;
	const int delivery = incrementalRefinement ?
	    mesh_lod_bounded_delivery_cut(payload->progressiveMesh,
		payload->activeCut, cut) : cut;
	return payload->progressiveMesh->canDrawCut(delivery) ?
	    delivery : payload->activeCut;
    }
    return payload->activeCut;
}

static int
mesh_lod_available_mesh_retarget_cut(
	const BObolViewLodState::MeshPayload *payload,
	const BObolLodRequest &request,
	SbBool incrementalRefinement)
{
    if (!payload || !payload->progressiveMesh ||
	!payload->progressiveMesh->isValid() ||
	request.requestedCut < 0)
	return request.requestedCut;
    if (request.requestedCut <= payload->activeCut)
	return request.requestedCut;
    for (int cut = request.requestedCut;
	 cut > payload->activeCut; --cut) {
	if (!payload->progressiveMesh->canDrawCut(cut))
	    continue;
	const int delivery = incrementalRefinement ?
	    mesh_lod_bounded_delivery_cut(payload->progressiveMesh,
		payload->activeCut, cut) : cut;
	return payload->progressiveMesh->canDrawCut(delivery) ?
	    delivery : payload->activeCut;
    }
    return payload->activeCut;
}

static int
mesh_lod_error_ceiling_cut(const BObolLodRequest &request,
    const BObolLodProgressiveMeshPtr &progressiveMesh,
    int emphasis, int minimumCut, int maximumCut,
    double maximumNormalizedError)
{
    if (minimumCut >= maximumCut ||
	!std::isfinite(maximumNormalizedError) ||
	maximumNormalizedError <= 0.0 ||
	!std::isfinite(request.projectedPixelDiameter) ||
	request.projectedPixelDiameter <= 0.0f)
	return maximumCut;
    const double emphasisWeight = emphasis >= 2 ? 4.0 :
	(emphasis == 1 ? 2.0 : 1.0);
    const double target = std::max(1.0,
	static_cast<double>(request.targetPixelError));
    if (progressiveMesh) {
	const double allowedPixelError =
	    target * maximumNormalizedError / emphasisWeight;
	const int selected = progressiveMesh->cutForScreenError(
	    request.projectedPixelDiameter, allowedPixelError);
	return selected < 0 ? maximumCut :
	    std::max(minimumCut, std::min(maximumCut, selected));
    }
    const double required = std::log2(
	emphasisWeight *
	static_cast<double>(request.projectedPixelDiameter) /
	(target * maximumNormalizedError));
    if (!std::isfinite(required))
	return required > 0.0 ? maximumCut : minimumCut;
    const double cutRequired = required * 3.0;
    const int cut = cutRequired >= static_cast<double>(INT_MAX) ? INT_MAX :
	(cutRequired <= static_cast<double>(INT_MIN) ? INT_MIN :
	 static_cast<int>(std::ceil(cutRequired)));
    return std::max(minimumCut, std::min(maximumCut, cut));
}

static int
mesh_lod_cad_allocated_cut(
    const BObolViewLodState::CadPayload *payload,
    const BObolLodRequest &request)
{
    if (!payload || payload->allocatedCut < 0 ||
	payload->allocationViewRevision != request.viewRevision.value() ||
	payload->allocationPolicyRevision != request.policyRevision.value() ||
	payload->allocationDrawMode != request.drawMode)
	return -1;
    return payload->allocatedCut;
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
	const int sourceDrawMode = source->getEffectiveLodDrawMode();
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
	 * order.  A caller may retain structural proxies as its absolute floor,
	 * or preserve one minimum PoP prefix per visible occurrence and let the
	 * renderer's aggregate point channel bound that mesh floor. */
	if (submitAction->compactEntryFirst == 0 &&
	    submitAction->compactEntryLimit == SIZE_MAX &&
	    submitAction->retainedSceneCostBudget != SIZE_MAX &&
	    submitAction->viewState) {
	    struct RetainedCandidate {
		const BObolViewLodState::CadPayload *payload;
		SbString occurrenceKey;
		float projectedPixels;
		int emphasis;
		int minimumCut;
		int requestedCut;
		BObolLodRequest demand;
		size_t minimumFaces;
		size_t minimumCost;
		int admittedCut;
		size_t admittedFaces;
		size_t admittedCost;
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
		mesh_lod_apply_cut_hysteresis(projected,
		    payload->progressiveMesh, payload->activeCut);
		RetainedCandidate candidate;
		candidate.payload = payload;
		candidate.occurrenceKey = summary.sourceInstanceKey;
		candidate.projectedPixels = projected.projectedPixelDiameter;
		candidate.emphasis = summary.selected ? 2 :
		    (summary.highlighted ? 1 : 0);
		candidate.minimumCut = payload->activeCut;
		candidate.requestedCut = projected.requestedCut;
		projected.viewRevision = submitAction->viewRevision;
		projected.policyRevision = submitAction->policyRevision;
		projected.drawMode = payload->drawMode;
		candidate.demand = projected;
		candidate.minimumFaces = payload->counts.faceCount;
		candidate.minimumCost = bobol_lod_render_cost_units(
		    payload->counts, payload->drawMode, 1);
		candidate.admittedCut = candidate.minimumCut;
		candidate.admittedFaces = candidate.minimumFaces;
		candidate.admittedCost = candidate.minimumCost;
		if (payload->progressiveMesh) {
		    candidate.minimumCut =
			payload->progressiveMesh->minimumCut();
		    if (!payload->progressiveMesh->canDrawCut(
			    candidate.minimumCut)) {
			if (submitAction->viewState->removeCadPayload(payload))
			    submitAction->updatedCutCount++;
			continue;
		    }
		    const BObolLodCounts minimumCounts =
			bobol_lod_progressive_counts(
			    payload->progressiveMesh,
			    candidate.minimumCut,
			    payload->hasNormals);
		    candidate.minimumFaces = minimumCounts.faceCount;
		    candidate.minimumCost = bobol_lod_render_cost_units(
			minimumCounts, payload->drawMode, 1);
		    candidate.admittedCut = candidate.minimumCut;
		    candidate.admittedFaces = candidate.minimumFaces;
		    candidate.admittedCost = candidate.minimumCost;
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
	    /* Reserve coverage for every visible occurrence first, then preserve
	     * the richest already-active prefix which fits the remaining scene
	     * allowance in screen-priority order.  The old recovery path reset
	     * every retained mesh to its minimum prefix and relied on later passes
	     * to walk upward again.  Even a single Lucy therefore flashed from a
	     * cut 7 to cut 0 before eventually returning to cut 6.
	     *
	     * This is still a recovery operation: it never promotes beyond the
	     * existing active cut, never reads or rebuilds geometry, and keeps the
	     * true projected request as the unsatisfied target.  Zero-cost PoP snap
	     * improvements are retained automatically. */
	    size_t minimumSceneCost = 0;
	    for (const RetainedCandidate &candidate : retainedCandidates)
		minimumSceneCost = candidate.minimumCost >
			SIZE_MAX - minimumSceneCost ?
		    SIZE_MAX : minimumSceneCost + candidate.minimumCost;
	    size_t upgradeBudget =
		submitAction->retainedSceneCostBudget > minimumSceneCost ?
		    submitAction->retainedSceneCostBudget - minimumSceneCost : 0;
	    if (minimumSceneCost <= submitAction->retainedSceneCostBudget) {
		for (RetainedCandidate &candidate : retainedCandidates) {
		    if (!candidate.payload->progressiveMesh)
			continue;
		    /* Recovery reallocates the existing resident population; it is
		     * not limited to preserving the cut this occurrence happened to
		     * own in the preceding view.  Newly prominent parts must be able
		     * to reclaim budget from newly insignificant ones without loading
		     * or rebuilding geometry. */
		    int targetCut = std::min(
			candidate.requestedCut,
			std::max(candidate.payload->activeCut,
			    candidate.payload->residentCut));
		    targetCut = mesh_lod_error_ceiling_cut(
			candidate.demand,
			candidate.payload->progressiveMesh,
			candidate.emphasis,
			candidate.minimumCut, targetCut,
			submitAction->retainedSceneMaximumNormalizedError);
		    for (int cut = targetCut;
			 cut >= candidate.minimumCut; --cut) {
			if (!candidate.payload->progressiveMesh->canDrawCut(cut))
			    continue;
			const BObolLodCounts counts =
			    bobol_lod_progressive_counts(
				candidate.payload->progressiveMesh, cut,
				candidate.payload->hasNormals);
			const size_t cost = bobol_lod_render_cost_units(
			    counts, candidate.payload->drawMode, 1);
			const size_t extra = cost > candidate.minimumCost ?
			    cost - candidate.minimumCost : 0;
			if (extra > upgradeBudget)
			    continue;
			candidate.admittedCut = cut;
			candidate.admittedFaces = counts.faceCount;
			candidate.admittedCost = cost;
			upgradeBudget -= extra;
			break;
		    }
		}
	    }
	    for (const RetainedCandidate &candidate : retainedCandidates) {
		const size_t remaining =
		    submitAction->retainedSceneCostBudgetUsed >=
			    submitAction->retainedSceneCostBudget ?
			0 : submitAction->retainedSceneCostBudget -
			    submitAction->retainedSceneCostBudgetUsed;
		const bool exceedsBudget =
		    candidate.admittedCost > remaining;
		if (exceedsBudget) {
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->refinementBudgetBlockedCount++;
		    if (!submitAction->preserveMeshCoverage) {
			if (submitAction->viewState->removeCadPayload(
				candidate.payload))
			    submitAction->updatedCutCount++;
			continue;
		    }
		}
		submitAction->retainedSceneCostBudgetUsed =
		    candidate.admittedCost >
			    SIZE_MAX -
				submitAction->retainedSceneCostBudgetUsed ?
			SIZE_MAX :
			submitAction->retainedSceneCostBudgetUsed +
			    candidate.admittedCost;
		submitAction->retainedRecoveredOccurrences.insert(
		    candidate.occurrenceKey.getString());
		const int priorCut = candidate.payload->activeCut;
		const int retainedCut =
		    priorCut != candidate.admittedCut &&
		    candidate.admittedFaces < candidate.payload->counts.faceCount ?
			candidate.admittedCut : priorCut;
		const SbBool retargeted =
		    submitAction->viewState->retargetCadPayload(
			candidate.payload, retainedCut,
			candidate.demand);
		if (retargeted && priorCut != retainedCut)
		    submitAction->updatedCutCount++;
		if (candidate.admittedCut < candidate.requestedCut)
		    submitAction->pendingRetainedRefinementCount++;
	    }
	}
	struct CompactCandidate {
	    size_t index;
	    float projectedPixels;
	    double projectedErrorPixels;
	    double visualFootprint;
	    double refinementValuePerCost;
	    int emphasis;
	    bool needsCoverage;
	    bool qualityFloorViolation;
	};
	if (!submitAction->compactEntryPlanSupplied) {
	    /* This traversal owns the perceptual order.  Coverage has already
	     * established one mesh per visible occurrence; visit the strongest
	     * deficits first and let each choose its richest budget-fitting
	     * resident cut.  Repeating a whole-scene pass for every nominal PoP
	     * transition made stable convergence proportional to hierarchy depth. */
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
		    const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		    if (filter && filter[0] &&
			strstr(summary.sourceInstanceKey.getString(), filter)) {
			const SbVec3f minimum = summary.localBounds.getMin();
			const SbVec3f maximum = summary.localBounds.getMax();
			SbVec3f origin;
			summary.localToSource.multVecMatrix(
			    SbVec3f(0.0f, 0.0f, 0.0f), origin);
			bu_log("BObol LoD planning trace occurrence=%s "
			       "rejected=projection bounds=(%.9g %.9g %.9g)-"
			       "(%.9g %.9g %.9g) origin=(%.9g %.9g %.9g)\n",
			       summary.sourceInstanceKey.getString(),
			       minimum[0], minimum[1], minimum[2],
			       maximum[0], maximum[1], maximum[2],
			       origin[0], origin[1], origin[2]);
		    }
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
		candidate.visualFootprint = std::max(
		    std::sqrt(std::max(0.0,
			static_cast<double>(projected.projectedPixelArea))),
		    std::max(
			static_cast<double>(projected.projectedPixelPerimeter) *
			    0.25,
			static_cast<double>(candidate.projectedPixels) * 0.25));
		candidate.emphasis = summary.selected ? 2 :
		    (summary.highlighted ? 1 : 0);
		const BObolViewLodState::CadPayload *active =
		    submitAction->viewState ?
		    submitAction->viewState->findCadForOccurrence(
			source, summary.sourceInstanceKey) : NULL;
		candidate.needsCoverage = active == NULL;
		candidate.projectedErrorPixels = active &&
		    active->activeCut >= 0 && active->progressiveMesh ?
		    active->progressiveMesh->projectedErrorAtCut(
			active->activeCut, candidate.projectedPixels) :
		    std::numeric_limits<double>::infinity();
		candidate.qualityFloorViolation = active &&
		    candidate.projectedErrorPixels >
			std::max(1.0, static_cast<double>(
			    projected.targetPixelError) * 3.0);
		candidate.refinementValuePerCost =
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
		 * impact without scanning triangles.  Exact per-cut geometric
		 * delta metrics may later be cached by the PoP builder, but are
		 * not required in the hot view path. */
		if (active && active->progressiveMesh &&
		    active->activeCut >= 0 &&
		    projected.requestedCut > active->activeCut) {
		    const int nextCut = mesh_lod_next_population_cut(
			active->progressiveMesh, active->activeCut,
			projected.requestedCut);
		    const BObolLodCounts activeCounts =
			bobol_lod_progressive_counts(
			    active->progressiveMesh, active->activeCut,
			    active->hasNormals);
		    const BObolLodCounts nextCounts =
			bobol_lod_progressive_counts(
			    active->progressiveMesh, nextCut,
			    active->hasNormals);
		    const size_t activeCost = bobol_lod_render_cost_units(
			activeCounts, sourceDrawMode, 0);
		    const size_t nextCost = bobol_lod_render_cost_units(
			nextCounts, sourceDrawMode, 0);
		    const size_t addedCost = nextCost > activeCost ?
			nextCost - activeCost : 1;
		    const double currentError =
			active->progressiveMesh->projectedErrorAtCut(
			    active->activeCut,
			    candidate.projectedPixels);
		    const double nextError =
			active->progressiveMesh->projectedErrorAtCut(
			    nextCut, candidate.projectedPixels);
		    const double visibleBenefit =
			std::max(1.0, candidate.visualFootprint) *
			std::max(0.0, currentError - nextError);
		    candidate.refinementValuePerCost =
			visibleBenefit / static_cast<double>(addedCost);
		}
		candidates.push_back(std::move(candidate));
	    }
	    /* Scene coverage is the first progressive objective.  For retained
	     * meshes use a minimax quality floor before benefit/cost: this prevents
	     * an expensive but unmistakable wheel, blade, or tail from losing
	     * forever to many cheap low-impact transitions.  Once the conservative
	     * three-pixel floor is met, spend the remainder by marginal projected
	     * footprint/error reduction.  Selection and highlight remain absolute
	     * user-intent priorities. */
	    std::sort(candidates.begin(), candidates.end(),
		[](const CompactCandidate &a, const CompactCandidate &b) {
		if (a.emphasis != b.emphasis)
		    return a.emphasis > b.emphasis;
		if (a.needsCoverage != b.needsCoverage)
		    return a.needsCoverage;
		if (a.qualityFloorViolation != b.qualityFloorViolation)
		    return a.qualityFloorViolation;
		if (a.qualityFloorViolation &&
		    (a.projectedErrorPixels < b.projectedErrorPixels ||
		     a.projectedErrorPixels > b.projectedErrorPixels))
		    return a.projectedErrorPixels > b.projectedErrorPixels;
		if (a.refinementValuePerCost <
			b.refinementValuePerCost ||
		    a.refinementValuePerCost >
			b.refinementValuePerCost)
		    return a.refinementValuePerCost >
			b.refinementValuePerCost;
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
	const std::vector<size_t> &entryPlan =
	    submitAction->compactEntryPlanView ?
		*submitAction->compactEntryPlanView :
		submitAction->compactEntryPlan;
	const size_t candidateCount = entryPlan.size();
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
	const SbBool suppressActiveDuplicate =
	    (!submitAction->useForcedCut && submitAction->reset == 0) ?
	    TRUE : FALSE;
	std::vector<BObolLodTask> pendingTasks;
	const size_t taskCapacity =
	    submitAction->submissionTaskLimit == SIZE_MAX ?
		localLast - localFirst :
		std::min(localLast - localFirst,
		    submitAction->submissionTaskLimit);
	pendingTasks.reserve(taskCapacity);
	const int64_t submissionStarted = bu_gettime();
	for (size_t candidateOffset = localFirst;
	    candidateOffset < localLast; candidateOffset++) {
	    if (candidateOffset > localFirst &&
		submitAction->submissionTimeLimitMicroseconds > 0) {
		const int64_t elapsed = bu_gettime() - submissionStarted;
		if (elapsed >= 0 &&
		    static_cast<uint64_t>(elapsed) >=
			submitAction->submissionTimeLimitMicroseconds) {
		    /* Do not consume this entry.  The controller resumes the
		     * pinned plan at exactly this occurrence next pump. */
		    processedLast = candidateOffset;
		    break;
		}
	    }
	    processedLast = candidateOffset + 1;
	    const size_t i = entryPlan[candidateOffset];
	    BObolCompactLodInstanceSummary summary;
	    if (!source->getCompactLodInstanceSummary(
		    static_cast<int>(i), summary))
		continue;
	    if (!summary.valid || !summary.visible) {
		if ((submitAction->retainedSceneCostBudget != SIZE_MAX ||
		     submitAction->retainedSceneUpgradeCostBudget != SIZE_MAX) &&
		    submitAction->viewState &&
		    summary.sourceInstanceKey.getLength() > 0) {
		    const BObolViewLodState::CadPayload *inactive =
			submitAction->viewState->findCadForOccurrence(
			    source, summary.sourceInstanceKey);
		    if (inactive &&
			submitAction->viewState->removeCadPayload(inactive))
			submitAction->updatedCutCount++;
		}
		continue;
	    }

	    submitAction->visitedMeshCount++;
	    const SbString &target = summary.path.getLength() > 0 ?
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
	    BObolSourceMeshRequest sourceMeshRequest;
	    const SbBool haveSourceMeshRequest =
		source->getCompactSourceMeshRequest(
		    static_cast<int>(i), sourceMeshRequest);

	    BObolLodRequest request;
	    request.databaseId = submitAction->databaseId;
	    request.databaseRevision = submitAction->databaseRevision;
	    request.sourceRevision = source->sourceRevision.getValue();
	    /* Source-backed meshes carry their exact immutable cache identity when
	     * available.  In particular, a generated BREP's tessellation-variant
	     * key is independent of the standing AABB part which the result will
	     * replace. */
	    request.sourceContentHash = summary.sourceContentHash;
	    request.objectPath = summary.meshAssetPath.getLength() > 0 ?
		summary.meshAssetPath : target;
	    request.objectName = summary.meshAssetName.getLength() > 0 ?
		summary.meshAssetName : summary.sourceName;
	    request.occurrenceKey = summary.sourceInstanceKey;
	    request.sourceRoutingId = source->getCompactSourceRoutingId();
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
		    if ((submitAction->retainedSceneCostBudget != SIZE_MAX ||
			 submitAction->retainedSceneUpgradeCostBudget != SIZE_MAX) &&
			submitAction->viewState) {
			const BObolViewLodState::CadPayload *offscreen =
			    submitAction->viewState->findCadForOccurrence(
				source, summary.sourceInstanceKey);
			if (offscreen &&
			    submitAction->viewState->
				removeCadPayload(offscreen))
			    submitAction->updatedCutCount++;
		    }
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
	    submitAction->visibleMeshCount++;
	    if (submitAction->useForcedCut)
		request.requestedCut = submitAction->forcedCut;
	else if (activePayload)
		mesh_lod_apply_cut_hysteresis(request,
		    activePayload->progressiveMesh,
		    activePayload->activeCut);
	    const MeshLodBrepVariantDemand brepVariant =
		haveSourceMeshRequest ? mesh_lod_brep_variant_demand(
		    sourceMeshRequest, request, activePayload,
		    submitAction->service,
		    submitAction->allowRepresentationRefinement) :
		MeshLodBrepVariantDemand();
	    if (brepVariant.contentHash)
		request.sourceContentHash = brepVariant.contentHash;
	    const SbBool terminalFailure =
		submitAction->viewState &&
		submitAction->viewState->
		    hasCadOccurrenceTerminalFailure(source, request);
	    if (activePayload || terminalFailure)
		submitAction->coveredVisibleMeshCount++;
	    if (terminalFailure) {
		/* The source's structural occurrence remains the useful fallback.
		 * A different demand epoch or cut is intentionally not suppressed. */
		submitAction->skippedMeshCount++;
		continue;
	    }
	    if (!submitAction->useForcedCut &&
		mesh_lod_payload_memory_limited_for_epoch(
		    activePayload, request, submitAction->service)) {
		/* The retained minimum/richer prior cut remains the stable
		 * presentation.  Retrying the same suffix in the same capacity
		 * epoch creates an unbounded worker loop without changing a
		 * pixel.  A camera/policy change or reclaimed stable bytes
		 * invalidates this suppression. */
		submitAction->skippedMeshCount++;
		continue;
	    }

	    /* Bounded large-scene recovery reserves the minimum drawable prefix of
	     * every retained occurrence before it reaches this window.  Spend only
	     * the carried upgrade allowance here, preserving the richest part of
	     * the already-active cut which fits.  This avoids the former
	     * minimum-then-refine cycle and makes the pinned plan's visual-priority
	     * order authoritative across all windows. */
	    if (activePayload && activePayload->progressiveMesh &&
		submitAction->retainedSceneUpgradeCostBudget != SIZE_MAX) {
	    const int minimumCut =
		activePayload->progressiveMesh->minimumCut();
	    const int allocatedCut = mesh_lod_cad_allocated_cut(
		activePayload, request);
	    int desiredCut = request.requestedCut;
	    if (allocatedCut >= minimumCut)
		desiredCut = std::min(desiredCut, allocatedCut);
	    else
		desiredCut = mesh_lod_error_ceiling_cut(
		    request, activePayload->progressiveMesh,
		    summary.selected ? 2 :
			(summary.highlighted ? 1 : 0),
		    minimumCut, desiredCut,
		    submitAction->retainedSceneMaximumNormalizedError);
	    desiredCut = std::max(minimumCut, desiredCut);
	    const SbBool sceneQualityLimited =
		desiredCut < request.requestedCut ? TRUE : FALSE;
	    const int availableCut = std::max(activePayload->activeCut,
		activePayload->residentCut);
	    if (desiredCut > availableCut)
		submitAction->pendingResidentRefinementCount++;
	    const int targetCut = std::max(minimumCut,
		std::min(desiredCut, availableCut));
	    const BObolLodCounts minimumCounts =
		    bobol_lod_progressive_counts(
			activePayload->progressiveMesh, minimumCut,
			activePayload->hasNormals);
	    const size_t minimumCost = bobol_lod_render_cost_units(
		    minimumCounts, request.drawMode, 1);
	    const size_t remainingUpgrade =
		    submitAction->retainedSceneUpgradeCostBudgetUsed >=
			submitAction->retainedSceneUpgradeCostBudget ? 0 :
		    submitAction->retainedSceneUpgradeCostBudget -
			submitAction->retainedSceneUpgradeCostBudgetUsed;
	    int admittedCut = minimumCut;
	    BObolLodCounts admittedCounts = minimumCounts;
	    size_t admittedCost = minimumCost;
	    if (!submitAction->allowCutDowngrade) {
		    admittedCut = activePayload->activeCut;
		    admittedCounts = activePayload->counts;
		    admittedCost = bobol_lod_render_cost_units(
			admittedCounts, request.drawMode, 1);
		} else {
		    for (int cut = targetCut; cut >= minimumCut; --cut) {
			if (!activePayload->progressiveMesh->canDrawCut(cut))
			    continue;
			const BObolLodCounts counts =
			    bobol_lod_progressive_counts(
				activePayload->progressiveMesh, cut,
				activePayload->hasNormals);
			const size_t cost = bobol_lod_render_cost_units(
			    counts, request.drawMode, 1);
			const size_t extra = cost > minimumCost ?
			    cost - minimumCost : 0;
			if (extra > remainingUpgrade)
			    continue;
			admittedCut = cut;
			admittedCounts = counts;
			admittedCost = cost;
			break;
		    }
		}
		const size_t upgradeCost = admittedCost > minimumCost ?
		    admittedCost - minimumCost : 0;
		const char *retainedTrace = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (retainedTrace && retainedTrace[0] &&
		    (strstr(request.objectName.getString(), retainedTrace) ||
		     strstr(request.objectPath.getString(), retainedTrace)))
		    bu_log("BObol LoD retained admission object=%s "
			   "minimum_cut=%d target_cut=%d admitted_cut=%d "
			   "active_cut=%d resident_cut=%d requested_cut=%d "
			   "ceiling=%.9g minimum_cost=%zu admitted_cost=%zu "
			   "upgrade_cost=%zu remaining_upgrade=%zu "
			   "used_upgrade=%zu budget_upgrade=%zu\n",
			   request.objectName.getString(), minimumCut,
			   targetCut, admittedCut,
			   activePayload->activeCut,
			   activePayload->residentCut,
			   request.requestedCut,
			   submitAction->retainedSceneMaximumNormalizedError,
			   minimumCost, admittedCost, upgradeCost,
			   remainingUpgrade,
			   submitAction->retainedSceneUpgradeCostBudgetUsed,
			   submitAction->retainedSceneUpgradeCostBudget);
		if (upgradeCost > remainingUpgrade) {
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->refinementBudgetBlockedCount++;
		    submitAction->retainedAdmissionBlockedCount++;
		} else {
		    submitAction->retainedSceneUpgradeCostBudgetUsed =
			upgradeCost > SIZE_MAX -
				submitAction->retainedSceneUpgradeCostBudgetUsed ?
			    SIZE_MAX :
			    submitAction->retainedSceneUpgradeCostBudgetUsed +
				upgradeCost;
		    submitAction->retainedRecoveredOccurrences.insert(
			request.occurrenceKey.getString());
		    const int retainedRequestedCut = std::max(
			admittedCut, request.requestedCut);
		    BObolLodRequest retainedDemand = request;
		    retainedDemand.requestedCut = retainedRequestedCut;
		    const int priorCut = activePayload->activeCut;
		    const SbBool retained =
			submitAction->viewState->retargetCadPayload(
			    activePayload, admittedCut,
			    retainedDemand);
		    if (retained && priorCut != admittedCut)
			submitAction->updatedCutCount++;
		    if (retained && admittedCut < retainedRequestedCut) {
			submitAction->pendingRetainedRefinementCount++;
			/* A minimax ceiling below pixel demand is a real scene-budget
			 * limit even when every discrete allocation selected by that
			 * ceiling fits exactly.  Feed it through the bounded unchanged-
			 * frame calibration path so a conservative cold seed can grow;
			 * that policy still retires the observation after three samples
			 * when the renderer is genuinely saturated.  Missing residency is
			 * kept separate above and starts provider work instead. */
			const SbBool admissionLimited =
			    admittedCut < targetCut ? TRUE : FALSE;
			if (admissionLimited)
			    submitAction->retainedAdmissionBlockedCount++;
			if (sceneQualityLimited)
			    submitAction->retainedQualityLimitedCount++;
			if (admissionLimited || sceneQualityLimited)
			    submitAction->refinementBudgetBlockedCount++;
		    }
		    if (retained) {
			mesh_lod_trace_projected_request(
			    request, summary.localBounds, localToRoot,
			    priorCut);
			submitAction->skippedMeshCount++;
			continue;
		    }
		}
	    }

	    /* A finite refinement allowance prevents new growth, but by itself
	     * cannot recover when thousands of already-resident minimum cuts
	     * collectively exceed newly calibrated scene capacity.  During that
	     * recovery pass, retarget retained occurrences in pinned
	     * screen-priority order.  Entries not admitted retain their
	     * structural leaf proxies. */
	    if (activePayload &&
		submitAction->retainedRecoveredOccurrences.find(
		    request.occurrenceKey.getString()) ==
		    submitAction->retainedRecoveredOccurrences.end() &&
		submitAction->retainedSceneCostBudget != SIZE_MAX) {
		size_t admittedFaces = activePayload->counts.faceCount;
		BObolLodCounts admittedCounts = activePayload->counts;
		int admittedCut = request.requestedCut;
		const int targetCut = request.requestedCut;
		if (activePayload->progressiveMesh) {
		    const int minimumCut =
			activePayload->progressiveMesh->minimumCut();
		    /* Recovery is explicitly coverage-first.  Allocating each
		     * early occurrence's richest affordable target can exhaust
		     * the scene budget before later leaves are visited, making
		     * loaded geometry flash back to structural boxes.  Re-admit
		     * one coherent minimum PoP mesh per occurrence first.  A
		     * following current-view pass spends the remaining allowance
		     * on marginal error reduction in the pinned priority order. */
		    admittedCut = minimumCut;
		    admittedCounts = bobol_lod_progressive_counts(
			activePayload->progressiveMesh, admittedCut,
			activePayload->hasNormals);
		    admittedFaces = admittedCounts.faceCount;
		    /*
		     * If this caller has explicitly prohibited a cut downgrade,
		     * or the minimum prefix saves no submitted faces, retained
		     * admission can account only for what will actually remain
		     * active.  Controller recovery passes enable downgrades;
		     * pose-only passes do not enter this block.  Keeping a richer
		     * snap cut at identical draw cost also avoids a quality loss
		     * with no corresponding FPS benefit.
		     */
		    if (!submitAction->allowCutDowngrade ||
			admittedFaces >= activePayload->counts.faceCount) {
			admittedCut = activePayload->activeCut;
			admittedFaces = activePayload->counts.faceCount;
			admittedCounts = activePayload->counts;
		    }
		}
		const size_t admittedCost = bobol_lod_render_cost_units(
		    admittedCounts, request.drawMode, 1);
		const size_t remaining =
		    submitAction->retainedSceneCostBudgetUsed >=
			    submitAction->retainedSceneCostBudget ?
			0 : submitAction->retainedSceneCostBudget -
			    submitAction->retainedSceneCostBudgetUsed;
		const SbBool admissionBlocked =
		    admittedCost > remaining ? TRUE : FALSE;
		if (admissionBlocked) {
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->refinementBudgetBlockedCount++;
		    if (!submitAction->preserveMeshCoverage) {
			if (submitAction->viewState->
				removeCadPayload(activePayload)) {
			    submitAction->updatedCutCount++;
			    if (submitAction->coveredVisibleMeshCount > 0)
				submitAction->coveredVisibleMeshCount--;
			}
			submitAction->retainedRecoveredOccurrences.insert(
			    request.occurrenceKey.getString());
			continue;
		    }
		}
		submitAction->retainedSceneCostBudgetUsed =
		    admittedCost >
			    SIZE_MAX -
				submitAction->retainedSceneCostBudgetUsed ?
			SIZE_MAX :
			submitAction->retainedSceneCostBudgetUsed +
			    admittedCost;
		submitAction->retainedRecoveredOccurrences.insert(
		    request.occurrenceKey.getString());
		if (activePayload->progressiveMesh) {
		    /*
		     * The admitted prefix is a temporary presentation cut, not a
		     * replacement for the view's desired cut.  Publishing the
		     * minimum as both activeCut and requestedCut made a cold
		     * coverage pass falsely satisfy the occurrence and erase it
		     * from the sparse refinement frontier.  Preserve the true
		     * projected demand so the following pass can refine it.
		     */
		    const int retainedRequestedCut =
			std::max(admittedCut, targetCut);
		    BObolLodRequest retainedDemand = request;
		    retainedDemand.requestedCut = retainedRequestedCut;
		    const int priorCut = activePayload->activeCut;
		    const SbBool retained =
			submitAction->viewState->retargetCadPayload(
			    activePayload, admittedCut,
			    retainedDemand);
		    if (retained && priorCut != admittedCut)
			submitAction->updatedCutCount++;
		    if (retained && !admissionBlocked &&
			admittedCut < retainedRequestedCut)
			submitAction->pendingRetainedRefinementCount++;
		    if (retained) {
			mesh_lod_trace_projected_request(
			    request, summary.localBounds, localToRoot,
			    priorCut);
			submitAction->skippedMeshCount++;
			continue;
		    }
		    /*
		     * A stale/unusable retained binding should not consume this
		     * occurrence forever.  Fall through with the original view
		     * target so the ordinary resident/provider route can repair
		     * it in this same bounded window.
		     */
		}
	    }
	    mesh_lod_trace_projected_request(request, summary.localBounds,
		localToRoot, activePayload ? activePayload->activeCut : -1);

	    const SbBool activeAssetMatches =
		mesh_lod_cad_payload_matches_asset_epoch(
		    activePayload, request);
	    const int allocatedPresentationCut =
		mesh_lod_cad_allocated_cut(activePayload, request);
	    BObolLodRequest presentationDemand = request;
	    if (allocatedPresentationCut >= 0)
		presentationDemand.requestedCut = std::min(
		    presentationDemand.requestedCut,
		    allocatedPresentationCut);
	    if (!submitAction->useForcedCut && activeAssetMatches &&
		activePayload->activeCut == request.requestedCut) {
		mesh_lod_retarget_cad_demand_if_changed(
		    submitAction->viewState, activePayload, request);
		submitAction->skippedMeshCount++;
		continue;
	    }
	    if (!submitAction->useForcedCut && activeAssetMatches &&
		allocatedPresentationCut >= 0 &&
		activePayload->activeCut >=
		    presentationDemand.requestedCut &&
		(!submitAction->allowResidentPrefetch ||
		 !activePayload->progressiveMesh ||
		 activePayload->progressiveMesh->canDrawCut(
		     request.requestedCut))) {
		/* View demand remains recorded at requestedCut, but this occurrence
		 * has reached its exact scene-budget allocation.  Do not let a later
		 * ordinary presentation pass spend another occurrence's share.  Quiet
		 * resident prefetch is independent: when the pixel-demanded suffix is
		 * not drawable yet, continue to the provider while keeping this
		 * allocated cut on screen. */
		/* An allocation stamp says why this occurrence is deliberately
		 * coarser than its physical view demand; it does not make that demand
		 * disappear.  Once the richer suffix is resident, report the remaining
		 * draw-budget debt on every complete policy pass.  The controller then
		 * measures one unchanged presentation and may enlarge/reallocate the
		 * scene allowance.  Silently classifying this as an ordinary skip made
		 * a warm cache settle at the first conservative allocation even when
		 * the retained frame took less than a millisecond (Lucy commonly froze
		 * at cut 14 with cut 20 already resident).
		 *
		 * Missing residency is intentionally excluded by the condition above:
		 * quiet prefetch continues through the provider path before a capacity
		 * probe is requested. */
		if (allocatedPresentationCut < request.requestedCut) {
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->refinementBudgetBlockedCount++;
		}
		mesh_lod_retarget_cad_demand_if_changed(
		    submitAction->viewState, activePayload, request);
		submitAction->skippedMeshCount++;
		continue;
	    }
	    if (!submitAction->useForcedCut &&
		!submitAction->allowCutDowngrade && activeAssetMatches &&
		activePayload->activeCut > request.requestedCut) {
		mesh_lod_retarget_cad_demand_if_changed(
		    submitAction->viewState, activePayload, request);
		submitAction->skippedMeshCount++;
		continue;
	    }
	    if (!submitAction->useForcedCut &&
		!submitAction->allowRetainedRefinement &&
		activeAssetMatches &&
		request.requestedCut > activePayload->activeCut) {
		/*
		 * Pose-only presentation retains the existing draw cut, but the
		 * current projected demand is still the input to the subsequent
		 * scene-wide importance allocation.  Updating this metadata changes
		 * neither the immutable mesh nor the prepared cut.
		 */
		mesh_lod_retarget_cad_demand_if_changed(
		    submitAction->viewState, activePayload, request);
		submitAction->skippedMeshCount++;
		continue;
	    }

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
		if (reusable && reusable->progressiveMesh &&
		    reusable->sourceContentHash == request.sourceContentHash) {
		    /* This occurrence had no hierarchy when its physical view demand
		     * was projected, but the shared resident asset does.  Resolve the
		     * producer-authored cut before attempting the zero-I/O binding;
		     * carrying -1 into prefix arithmetic silently strands repeated
		     * instances at their structural fallback. */
		    mesh_lod_apply_cut_hysteresis(request,
			reusable->progressiveMesh, -1);
		    const int minimumCut =
			reusable->progressiveMesh->minimumCut();
		    if (request.requestedCut < 0)
			request.requestedCut = minimumCut;
		    int reusedCut = std::min(request.requestedCut,
			reusable->progressiveMesh->maximumCut());
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
		    if (submitAction->retainedSceneCostBudget != SIZE_MAX ||
			submitAction->retainedSceneUpgradeCostBudget != SIZE_MAX)
			reusedCut = minimumCut;
		    while (reusedCut >= minimumCut &&
			!reusable->progressiveMesh->canDrawCut(reusedCut))
			reusedCut--;
		    if (reusedCut >= minimumCut) {
			const BObolLodCounts reusedCounts =
			    bobol_lod_progressive_counts(
				reusable->progressiveMesh, reusedCut,
				reusable->hasNormals);
			const size_t reusedFaces = reusedCounts.faceCount;
			const size_t reusedCost = bobol_lod_render_cost_units(
			    reusedCounts, request.drawMode, 1);
			SbBool admitted = TRUE;
			if (submitAction->retainedSceneUpgradeCostBudget !=
				SIZE_MAX) {
			    /* The controller reserved this minimum occurrence in
			     * the irreducible coverage floor.  Its cost is therefore
			     * zero in the carried detail allowance. */
			    admitted = TRUE;
			} else if (submitAction->retainedSceneCostBudget !=
				SIZE_MAX) {
			    const size_t remaining =
				submitAction->retainedSceneCostBudgetUsed >=
					submitAction->retainedSceneCostBudget ?
				    0 :
				    submitAction->retainedSceneCostBudget -
				    submitAction->
					    retainedSceneCostBudgetUsed;
			    if (reusedCost > remaining) {
				submitAction->refinementBudgetBlockedCount++;
				admitted =
				    submitAction->preserveMeshCoverage;
			    }
			    if (admitted)
				submitAction->retainedSceneCostBudgetUsed =
				    reusedCost >
					    SIZE_MAX -
						submitAction->
						    retainedSceneCostBudgetUsed ?
					SIZE_MAX :
					submitAction->
					    retainedSceneCostBudgetUsed +
					    reusedCost;
			} else {
			    admitted = submitAction->reserveInitialCost(
				reusedCounts, request.drawMode);
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
			reusedResult.geometry.activeCut = reusedCut;
			reusedResult.geometry.borrowed = TRUE;
			reusedResult.progressiveMesh =
			    reusable->progressiveMesh;
			reusedResult.resolvedCut = request.requestedCut;
			reusedResult.residentCut =
			    reusable->progressiveMesh->residentCut();
			reusedResult.resultKind = BOBOL_LOD_RESULT_MESH;
			reusedResult.qualityTier = std::max(
			    request.qualityTier, reusable->qualityTier);
			reusedResult.providerStatus =
			    BOBOL_LOD_PROVIDER_READY;
			reusedResult.bounds =
			    reusable->progressiveMesh->bounds();
			reusedResult.counts.faceCount = reusedFaces;
			reusedResult.counts.pointCount =
			    reusable->progressiveMesh->pointCount(reusedCut);
			reusedResult.counts.originalPointCount =
			    reusedResult.counts.pointCount;
			reusedResult.counts.normalCount =
			    reusable->hasNormals ? reusedFaces * 3 : 0;
			reusedResult.hasSnappedPoints = FALSE;
			reusedResult.hasNormals = reusable->hasNormals;
			reusedResult.shadedCullBackfaces =
			    reusable->shadedCullBackfaces;
			reusedResult.terminal =
			    reusedCut >= request.requestedCut ?
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
	    const SbBool retainedTargetDrawable =
		!submitAction->useForcedCut && activePayload &&
		activePayload->progressiveMesh &&
		presentationDemand.requestedCut > activePayload->activeCut &&
		activePayload->progressiveMesh->canDrawCut(
		    presentationDemand.requestedCut);
	    int retargetCut =
		mesh_lod_available_cad_retarget_cut(activePayload,
		    presentationDemand,
		    TRUE);
	    if (submitAction->reset == 0 && activePayload &&
		activeAssetMatches &&
		retargetCut > activePayload->activeCut) {
		retargetCut = submitAction->reserveRefinementCut(
		    activePayload->progressiveMesh,
		    activePayload->activeCut, retargetCut,
		    request.drawMode, activePayload->hasNormals);
		if (retargetCut <= activePayload->activeCut) {
		    mesh_lod_retarget_cad_demand_if_changed(
			submitAction->viewState, activePayload, request);
		    submitAction->pendingRetainedRefinementCount++;
		    submitAction->skippedMeshCount++;
		    continue;
		}
	    }
	    if (submitAction->reset == 0 && activePayload &&
		activeAssetMatches &&
		retargetCut >= 0 &&
		activePayload->activeCut != retargetCut) {
		mesh_lod_note_upward_resident_use(submitAction,
		    activePayload->cacheKey, activePayload->activeCut,
		    retargetCut);
		if (submitAction->viewState->retargetCadPayload(
		    activePayload, retargetCut, request)) {
		const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		if (filter && filter[0] &&
		    (strstr(request.objectName.getString(), filter) ||
		     strstr(request.objectPath.getString(), filter)))
		    bu_log("BObol LoD submit trace object=%s cut=%d "
			   "retargeted_resident=%d\n",
			   request.objectName.getString(),
			   retargetCut,
			   activePayload->progressiveMesh ?
			       activePayload->progressiveMesh->residentCut() : -1);
		submitAction->updatedCutCount++;
		if (retainedTargetDrawable &&
		    retargetCut < request.requestedCut) {
		    submitAction->pendingRetainedRefinementCount++;
		    continue;
		}
		/* If the full view target is already drawable, this is a pure cut
		 * change and no provider work is needed.  Otherwise keep the richer
		 * resident cut visible and continue below to request only the
		 * still-missing suffix. */
		if (retargetCut == request.requestedCut)
		    continue;
		}
	    }

	    int providerDeliveryCut = -1;
	    int providerPresentationCut = -1;
	    if (!submitAction->useForcedCut && activePayload &&
		activePayload->progressiveMesh &&
		request.requestedCut > activePayload->activeCut &&
		activeAssetMatches) {
		providerPresentationCut = activePayload->activeCut;
		if (presentationDemand.requestedCut >
			activePayload->activeCut) {
		    const int preferredCut = mesh_lod_bounded_delivery_cut(
			activePayload->progressiveMesh,
			activePayload->activeCut,
			presentationDemand.requestedCut);
		    providerPresentationCut =
			submitAction->reserveRefinementCut(
			    activePayload->progressiveMesh,
			    activePayload->activeCut, preferredCut,
			    request.drawMode, activePayload->hasNormals);
		}
		if (providerPresentationCut <= activePayload->activeCut) {
		    mesh_lod_retarget_cad_demand_if_changed(
			submitAction->viewState, activePayload, request);
		    submitAction->pendingRetainedRefinementCount++;
		    providerPresentationCut = activePayload->activeCut;
		}
		const SbBool requestedResident =
		    activePayload->progressiveMesh->canDrawCut(
			request.requestedCut);
		if (submitAction->allowResidentPrefetch && !requestedResident)
		    providerDeliveryCut = request.requestedCut;
		else if (providerPresentationCut > activePayload->activeCut)
		    providerDeliveryCut = providerPresentationCut;
		else {
		    submitAction->skippedMeshCount++;
		    continue;
		}
	    }
	    if (getenv("BOBOL_LOD_TRACE_BUDGET"))
		bu_log("BObol LoD provider admission object=%s occurrence=%s "
		       "active=%d requested=%d key_match=%d budget=%zu\n",
		       request.objectName.getString(),
		       request.occurrenceKey.getString(),
		       activePayload ? activePayload->activeCut : -1,
		       request.requestedCut,
		       activeAssetMatches ? 1 : 0,
		       submitAction->refinementCostBudget);

	    /*
	     * An in-flight cumulative-prefix load is coalesced by asset and
	     * occurrence, independently of camera revisions.  Record the newest
	     * view demand before submitBatch() can reject this task as an active
	     * duplicate.  Otherwise a wheel stream can outrun the worker while the
	     * retained occurrence still advertises the older cut.  The completed
	     * suffix is then safely rebased but has no current demand witness to
	     * admit it, so refinement appears only after button-up.  This is a
	     * metadata-only update; it neither journals a presentation mutation nor
	     * traverses the scene.
	     */
	    if (!submitAction->useForcedCut && activePayload &&
		activeAssetMatches &&
		request.requestedCut > activePayload->activeCut)
		mesh_lod_retarget_cad_demand_if_changed(
		    submitAction->viewState, activePayload, request);

	    if (static_cast<size_t>(submitAction->submittedTaskCount) +
		    pendingTasks.size() >=
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
		!submitAction->reserveInitialCost(
		    request.sourceCounts.faceCount,
		    request.sourceCounts.pointCount, request.drawMode,
		    initialFaceAllowance)) {
		submitAction->pendingRetainedRefinementCount++;
		submitAction->skippedMeshCount++;
		continue;
	    }

	    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
	    if (haveSourceMeshRequest &&
		request.sourceContentHash ==
		    sourceMeshRequest.meshAssetContentHash)
		provider->stagedSource = sourceMeshRequest.stagedSource.lock();
	    provider->meshAssetContentHash = request.sourceContentHash;
	    provider->generateBrepVariant = brepVariant.generate;
	    provider->brepTessellationAbsTol = brepVariant.absTol;
	    provider->brepTessellationRelTol = brepVariant.relTol;
	    provider->brepTessellationNormTol = brepVariant.normTol;
	    provider->brepVariantMemoryLimited = brepVariant.memoryLimited;
	    provider->service = submitAction->service;
	    if (!provider->setDatabase(submitAction->dbip)) {
		delete provider;
		submitAction->skippedMeshCount++;
		continue;
	    }
	    provider->refreshMissing = submitAction->refreshMissing;
	    provider->useForcedCut = submitAction->useForcedCut;
	    provider->forcedCut = submitAction->forcedCut;
	    provider->useCurrentDrawCut = TRUE;
	    provider->currentDrawCut =
		activePayload ? activePayload->activeCut : -1;
	    provider->useDeliveryCutLimit =
		providerDeliveryCut >= 0 ? TRUE : FALSE;
	    provider->deliveryCutLimit = providerDeliveryCut;
	    provider->usePresentationCutLimit =
		providerPresentationCut >= 0 && providerDeliveryCut >
		    providerPresentationCut ? TRUE : FALSE;
	    provider->presentationCutLimit = providerPresentationCut;
	    provider->prefetchCachedTargetOnFirstPublication =
		!activePayload ? TRUE : FALSE;
	    /*
	     * The scene allocator has already charged the complete marginal
	     * population of this exact cut.  Do not apply a second per-asset
	     * cut walker in the worker and under-deliver the reserved work.
	     * Cold first coverage, which has no active payload or explicit cut,
	     * retains the provider's independently bounded first-prefix policy.
	     */
	    if (provider->useDeliveryCutLimit)
		provider->progressiveDelivery = FALSE;
	    if (initialFaceAllowance) {
		provider->initialRefinementFaceBudget =
		    initialFaceAllowance;
		provider->initialRefinementPointBudget =
		    initialFaceAllowance > SIZE_MAX / 2 ?
			SIZE_MAX : initialFaceAllowance * 2;
	    }
	    provider->shrinkAfterCopy = TRUE;
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
		    const int workingSetCut = providerDeliveryCut >= 0 ?
			providerDeliveryCut : request.requestedCut;
		    task.estimatedWorkingSetBytes = brepVariant.generate ?
			brepVariant.estimatedWorkingSetBytes : (activePayload ?
			mesh_lod_resident_task_estimate(
			    activePayload->progressiveMesh, workingSetCut,
			    activePayload->hasNormals, request.drawMode) :
			(haveSourceMeshRequest ?
			    mesh_lod_warm_source_task_estimate(
				sourceMeshRequest, request) : 0));
		    pendingTasks.push_back(std::move(task));
	}

	if (!pendingTasks.empty()) {
	    std::vector<uint64_t> taskIds;
	    (void)submitAction->service->submitBatch(
		pendingTasks, taskIds, suppressActiveDuplicate);
	    for (size_t taskIndex = 0;
		    taskIndex < pendingTasks.size(); ++taskIndex) {
		BObolLodTask &task = pendingTasks[taskIndex];
		const uint64_t taskId =
		    taskIndex < taskIds.size() ? taskIds[taskIndex] : 0;
		const BObolLodRequest &request = task.request;
		if (!taskId) {
		    const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		    if (filter && filter[0] &&
			(strstr(request.objectName.getString(), filter) ||
			 strstr(request.objectPath.getString(), filter)))
			bu_log("BObol LoD submit trace object=%s cut=%d "
			       "rejected=service\n",
			       request.objectName.getString(),
			       request.requestedCut);
		    if (task.realizeDataFree && task.realizeData)
			(*task.realizeDataFree)(task.realizeData);
		    task.realizeData = NULL;
		    submitAction->skippedMeshCount++;
		} else {
		    const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
		    if (filter && filter[0] &&
			(strstr(request.objectName.getString(), filter) ||
			 strstr(request.objectPath.getString(), filter)))
			bu_log("BObol LoD submit trace object=%s cut=%d "
			       "task=%llu\n",
			       request.objectName.getString(),
			       request.requestedCut,
			       static_cast<unsigned long long>(taskId));
		    submitAction->submittedTaskCount++;
		}
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

    const int sourceDrawMode = source->getEffectiveLodDrawMode();
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
    if (submitAction->useForcedCut)
	request.requestedCut = submitAction->forcedCut;
    else if (activePayload)
	mesh_lod_apply_cut_hysteresis(request,
	    activePayload->progressiveMesh, activePayload->activeCut);
    if (!submitAction->useForcedCut &&
	mesh_lod_payload_memory_limited_for_epoch(
	    activePayload, request, submitAction->service)) {
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    if (!submitAction->useForcedCut &&
	mesh_lod_cad_payload_is_resident(activePayload, request)) {
	mesh_lod_retarget_cad_demand_if_changed(
	    submitAction->viewState, activePayload, request);
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    if (!submitAction->useForcedCut &&
	!submitAction->allowCutDowngrade && activePayload &&
	activePayload->activeCut > request.requestedCut &&
	mesh_lod_cad_payload_has_same_source(activePayload, request)) {
	mesh_lod_retarget_cad_demand_if_changed(
	    submitAction->viewState, activePayload, request);
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    if (!submitAction->useForcedCut &&
	!submitAction->allowRetainedRefinement && activePayload &&
	request.requestedCut > activePayload->activeCut &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request)) {
	mesh_lod_retarget_cad_demand_if_changed(
	    submitAction->viewState, activePayload, request);
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    const SbBool retainedTargetDrawable =
	!submitAction->useForcedCut && activePayload &&
	activePayload->progressiveMesh &&
	request.requestedCut > activePayload->activeCut &&
	activePayload->progressiveMesh->canDrawCut(request.requestedCut);
    int retargetCut =
	mesh_lod_available_cad_retarget_cut(activePayload, request,
	    retainedTargetDrawable);
    if (submitAction->reset == 0 && activePayload &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
	retargetCut > activePayload->activeCut) {
	retargetCut = submitAction->reserveRefinementCut(
	    activePayload->progressiveMesh, activePayload->activeCut,
	    retargetCut, request.drawMode, activePayload->hasNormals);
	if (retargetCut <= activePayload->activeCut) {
	    mesh_lod_retarget_cad_demand_if_changed(
		submitAction->viewState, activePayload, request);
	    submitAction->pendingRetainedRefinementCount++;
	    submitAction->skippedMeshCount++;
	    source->doAction(action);
	    return;
	}
    }
    if (submitAction->reset == 0 && activePayload &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
	retargetCut >= 0 &&
	activePayload->activeCut != retargetCut) {
	mesh_lod_note_upward_resident_use(submitAction,
	    activePayload->cacheKey, activePayload->activeCut,
	    retargetCut);
	if (submitAction->viewState->retargetCadPayload(
	    activePayload, retargetCut, request)) {
	submitAction->updatedCutCount++;
	if (retainedTargetDrawable && retargetCut < request.requestedCut) {
	    submitAction->pendingRetainedRefinementCount++;
	    source->doAction(action);
	    return;
	}
	if (retargetCut == request.requestedCut) {
	    source->doAction(action);
	    return;
	}
	}
    }

    const SbBool suppressActiveDuplicate =
	(!submitAction->useForcedCut && submitAction->reset == 0) ?
	TRUE : FALSE;
    if (suppressActiveDuplicate &&
	submitAction->service->hasActiveRequest(request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current compact LoD request is already active");
	source->doAction(action);
	return;
    }
    int providerDeliveryCut = -1;
    int providerPresentationCut = -1;
	if (!submitAction->useForcedCut && activePayload &&
	    activePayload->progressiveMesh &&
	request.requestedCut > activePayload->activeCut &&
	mesh_lod_geometry_key_matches(activePayload->cacheKey, request)) {
	const int preferredCut = mesh_lod_bounded_delivery_cut(
	    activePayload->progressiveMesh,
	    activePayload->activeCut, request.requestedCut);
	providerPresentationCut = submitAction->reserveRefinementCut(
	    activePayload->progressiveMesh, activePayload->activeCut,
	    preferredCut, request.drawMode, activePayload->hasNormals);
	if (providerPresentationCut <= activePayload->activeCut) {
	    mesh_lod_retarget_cad_demand_if_changed(
		submitAction->viewState, activePayload, request);
	    submitAction->pendingRetainedRefinementCount++;
	    if (!submitAction->allowResidentPrefetch ||
		activePayload->progressiveMesh->canDrawCut(
		    request.requestedCut)) {
		submitAction->skippedMeshCount++;
		source->doAction(action);
		return;
	    }
	    providerPresentationCut = activePayload->activeCut;
	}
	    providerDeliveryCut = submitAction->allowResidentPrefetch ?
		request.requestedCut : providerPresentationCut;
	}
	if (!submitAction->useForcedCut && activePayload &&
	    mesh_lod_geometry_key_matches(activePayload->cacheKey, request) &&
	    request.requestedCut > activePayload->activeCut)
	    mesh_lod_retarget_cad_demand_if_changed(
		submitAction->viewState, activePayload, request);

	size_t initialFaceAllowance = 0;
    if (!activePayload &&
	!submitAction->reserveInitialCost(
	    request.sourceCounts.faceCount,
	    request.sourceCounts.pointCount, request.drawMode,
	    initialFaceAllowance)) {
	submitAction->pendingRetainedRefinementCount++;
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }

    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
    provider->service = submitAction->service;
    if (!provider->setDatabase(submitAction->dbip)) {
	delete provider;
	submitAction->skippedMeshCount++;
	source->doAction(action);
	return;
    }
    provider->refreshMissing = submitAction->refreshMissing;
    provider->useForcedCut = submitAction->useForcedCut;
    provider->forcedCut = submitAction->forcedCut;
    provider->useCurrentDrawCut = TRUE;
    provider->currentDrawCut =
	activePayload ? activePayload->activeCut : -1;
    provider->useDeliveryCutLimit =
	providerDeliveryCut >= 0 ? TRUE : FALSE;
    provider->deliveryCutLimit = providerDeliveryCut;
    provider->usePresentationCutLimit =
	providerPresentationCut >= 0 && providerDeliveryCut >
	    providerPresentationCut ? TRUE : FALSE;
    provider->presentationCutLimit = providerPresentationCut;
    provider->prefetchCachedTargetOnFirstPublication =
	!activePayload ? TRUE : FALSE;
    if (provider->useDeliveryCutLimit)
	provider->progressiveDelivery = FALSE;
    if (initialFaceAllowance) {
	provider->initialRefinementFaceBudget = initialFaceAllowance;
	provider->initialRefinementPointBudget =
	    initialFaceAllowance > SIZE_MAX / 2 ?
		SIZE_MAX : initialFaceAllowance * 2;
    }
    provider->shrinkAfterCopy = TRUE;
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
	    if (activePayload) {
		const int workingSetCut = providerDeliveryCut >= 0 ?
		    providerDeliveryCut : request.requestedCut;
		task.estimatedWorkingSetBytes =
		    mesh_lod_resident_task_estimate(
			activePayload->progressiveMesh, workingSetCut,
			activePayload->hasNormals, request.drawMode);
	    }

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
    if (submitAction->useForcedCut)
	request.requestedCut = submitAction->forcedCut;
    else if (viewPayload)
	mesh_lod_apply_cut_hysteresis(request,
	    viewPayload->progressiveMesh, viewPayload->activeCut);
    mesh_lod_trace_projected_request(request, request.bounds,
	SbMatrix::identity(), viewPayload ? viewPayload->activeCut : -1);
    if (!submitAction->useForcedCut &&
	mesh_lod_payload_memory_limited_for_epoch(
	    viewPayload, request, submitAction->service)) {
	submitAction->skippedMeshCount++;
	return;
    }
    if (!submitAction->useForcedCut &&
	mesh_lod_mesh_payload_is_resident(viewPayload, request)) {
	submitAction->skippedMeshCount++;
	return;
    }
    if (!submitAction->useForcedCut &&
	!submitAction->allowCutDowngrade && viewPayload &&
	viewPayload->activeCut > request.requestedCut &&
	mesh_lod_mesh_payload_has_same_source(viewPayload, request)) {
	submitAction->skippedMeshCount++;
	return;
    }
    if (!submitAction->useForcedCut &&
	!submitAction->allowRetainedRefinement && viewPayload &&
	request.requestedCut > viewPayload->activeCut &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request)) {
	submitAction->skippedMeshCount++;
	return;
    }
    const SbBool retainedTargetDrawable =
	!submitAction->useForcedCut && viewPayload &&
	viewPayload->progressiveMesh &&
	request.requestedCut > viewPayload->activeCut &&
	viewPayload->progressiveMesh->canDrawCut(request.requestedCut);
    int retargetCut =
	mesh_lod_available_mesh_retarget_cut(viewPayload, request,
	    retainedTargetDrawable);
    if (submitAction->reset == 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request) &&
	retargetCut > viewPayload->activeCut) {
	retargetCut = submitAction->reserveRefinementCut(
	    viewPayload->progressiveMesh, viewPayload->activeCut,
	    retargetCut, request.drawMode, viewPayload->hasNormals);
	if (retargetCut <= viewPayload->activeCut) {
	    (void)submitAction->viewState->retargetMeshPayload(viewPayload,
		viewPayload->activeCut, request.requestedCut,
		request.viewRevision.value(), request.policyRevision.value());
	    submitAction->pendingRetainedRefinementCount++;
	    submitAction->skippedMeshCount++;
	    return;
	}
    }
    if (submitAction->reset == 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request) &&
	viewPayload->activeCut != retargetCut &&
	retargetCut >= 0) {
	mesh_lod_note_upward_resident_use(submitAction,
	    viewPayload->cacheKey, viewPayload->activeCut, retargetCut);
	if (submitAction->viewState->retargetMeshPayload(viewPayload,
	    retargetCut, request.requestedCut, request.viewRevision.value(),
	    request.policyRevision.value())) {
	submitAction->updatedCutCount++;
	if (retainedTargetDrawable && retargetCut < request.requestedCut) {
	    submitAction->pendingRetainedRefinementCount++;
	    return;
	}
	if (retargetCut == request.requestedCut)
	    return;
	}
    }
    if (!submitAction->useForcedCut && submitAction->reset == 0 &&
	request.requestedCut < 0 && viewPayload &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request)) {
	submitAction->skippedMeshCount++;
	return;
    }

    const SbBool suppressActiveDuplicate =
	(!submitAction->useForcedCut && submitAction->reset == 0) ?
	TRUE : FALSE;
    if (suppressActiveDuplicate &&
	submitAction->service->hasActiveRequest(request)) {
	submitAction->skippedMeshCount++;
	submitAction->appendDiagnostic(target,
				       "current LoD request is already active");
	return;
    }
    int providerDeliveryCut = -1;
    int providerPresentationCut = -1;
	if (!submitAction->useForcedCut && viewPayload &&
	    viewPayload->progressiveMesh &&
	request.requestedCut > viewPayload->activeCut &&
	mesh_lod_geometry_key_matches(viewPayload->cacheKey, request)) {
	const int preferredCut = mesh_lod_bounded_delivery_cut(
	    viewPayload->progressiveMesh,
	    viewPayload->activeCut, request.requestedCut);
	providerPresentationCut = submitAction->reserveRefinementCut(
	    viewPayload->progressiveMesh, viewPayload->activeCut,
	    preferredCut, request.drawMode, viewPayload->hasNormals);
	if (providerPresentationCut <= viewPayload->activeCut) {
	    (void)submitAction->viewState->retargetMeshPayload(viewPayload,
		viewPayload->activeCut, request.requestedCut,
		request.viewRevision.value(), request.policyRevision.value());
	    submitAction->pendingRetainedRefinementCount++;
	    if (!submitAction->allowResidentPrefetch ||
		viewPayload->progressiveMesh->canDrawCut(
		    request.requestedCut)) {
		submitAction->skippedMeshCount++;
		return;
	    }
	    providerPresentationCut = viewPayload->activeCut;
	}
	    providerDeliveryCut = submitAction->allowResidentPrefetch ?
		request.requestedCut : providerPresentationCut;
	}
	if (!submitAction->useForcedCut && viewPayload &&
	    mesh_lod_geometry_key_matches(viewPayload->cacheKey, request) &&
	    request.requestedCut > viewPayload->activeCut)
	    mesh_lod_retarget_mesh_demand_if_changed(
		submitAction->viewState, viewPayload, request);

	size_t initialFaceAllowance = 0;
    if (!viewPayload &&
	!submitAction->reserveInitialCost(
	    request.sourceCounts.faceCount,
	    request.sourceCounts.pointCount, request.drawMode,
	    initialFaceAllowance)) {
	submitAction->pendingRetainedRefinementCount++;
	submitAction->skippedMeshCount++;
	return;
    }

    BObolMeshLodProvider *provider = new BObolMeshLodProvider;
    provider->service = submitAction->service;
    if (!provider->setDatabase(submitAction->dbip)) {
	delete provider;
	submitAction->skippedMeshCount++;
	return;
    }
    provider->refreshMissing = submitAction->refreshMissing;
    provider->useForcedCut = submitAction->useForcedCut;
    provider->forcedCut = submitAction->forcedCut;
    provider->useCurrentDrawCut = TRUE;
    provider->currentDrawCut =
	viewPayload ? viewPayload->activeCut : -1;
    provider->useDeliveryCutLimit =
	providerDeliveryCut >= 0 ? TRUE : FALSE;
    provider->deliveryCutLimit = providerDeliveryCut;
    provider->usePresentationCutLimit =
	providerPresentationCut >= 0 && providerDeliveryCut >
	    providerPresentationCut ? TRUE : FALSE;
    provider->presentationCutLimit = providerPresentationCut;
    provider->prefetchCachedTargetOnFirstPublication =
	!viewPayload ? TRUE : FALSE;
    if (provider->useDeliveryCutLimit)
	provider->progressiveDelivery = FALSE;
    if (initialFaceAllowance) {
	provider->initialRefinementFaceBudget = initialFaceAllowance;
	provider->initialRefinementPointBudget =
	    initialFaceAllowance > SIZE_MAX / 2 ?
		SIZE_MAX : initialFaceAllowance * 2;
    }
    provider->shrinkAfterCopy = TRUE;
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
	    const int workingSetCut = providerDeliveryCut >= 0 ?
		providerDeliveryCut : request.requestedCut;
	    if (viewPayload) {
		task.estimatedWorkingSetBytes =
		    mesh_lod_resident_task_estimate(
			viewPayload->progressiveMesh, workingSetCut,
			viewPayload->hasNormals, request.drawMode);
	    } else if (shape->lodAvailable.getValue() &&
		shape->lodActiveCut.getValue() >= workingSetCut &&
		shape->lodFaceCount.getValue() > 0 &&
		shape->lodPointCount.getValue() > 0) {
		task.estimatedWorkingSetBytes =
		    mesh_lod_resident_working_set_estimate(
			shape->lodPointCount.getValue(),
			shape->lodFaceCount.getValue(),
			shape->lodHasNormals.getValue(), request.drawMode);
	    }

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
