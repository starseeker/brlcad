/*                 V I E W _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/log.h"
#include "bu/str.h"
#include "bu/datetime.h"

#include "bv.h"
#include "BObol/BDatabaseSource.h"
#include "BObol/BEditPreview.h"
#include "BObol/BHUDLabelOverlay.h"
#include "BObol/BLineLayerOverlay.h"
#include "BObol/BLodService.h"
#include "BObol/BLodUpdateAction.h"
#include "BObol/BMeshLodSubmitAction.h"
#include "BObol/BPickDetail.h"
#include "BObol/BSceneController.h"
#include "BObol/BViewAttachment.h"
#include "BObol/BViewController.h"
#include "BObol/BViewLod.h"
#include "BObol/BViewStore.h"
#include "cad_assembly_private.h"
#include "lod_coordinator_private.h"
#include "lod_control_private.h"
#include "retained_allocation_private.h"
#include "view_controller_private.h"
#include "raytrace.h"
#include "rt/view.h"

#include <algorithm>
#include <array>
#include <atomic>
#include <chrono>
#include <cmath>
#include <cstdlib>
#include <cstdint>
#include <limits>
#include <memory>
#include <mutex>
#include <queue>
#include <string.h>
#include <type_traits>
#include <unordered_map>
#include <utility>
#include <vector>

#include <Inventor/SbName.h>
#include <Inventor/SbPlane.h>
#include <Inventor/SoDB.h>
#include <Inventor/nodes/SoCamera.h>
#include <Inventor/nodes/SoClipPlane.h>
#include <Inventor/nodes/SoDepthBuffer.h>
#include <Inventor/nodes/SoDirectionalLight.h>
#include <Inventor/nodes/SoEnvironment.h>
#include <Inventor/nodes/SoGroup.h>
#include <Inventor/nodes/SoLight.h>
#include <Inventor/nodes/SoLightModel.h>
#include <Inventor/nodes/SoOrthographicCamera.h>
#include <Inventor/nodes/SoPointLight.h>
#include <Inventor/nodes/SoSpotLight.h>
#include <Inventor/nodes/SoPerspectiveCamera.h>
#include <Inventor/nodes/SoSeparator.h>

std::vector<SoBRLDatabaseSource *> controller_render_database_sources(
    const BObolViewController *controller);
std::vector<SoBRLDatabaseSource *> controller_render_database_source_roots(
    const BObolViewController *controller);

/* A normal live PoP preview is built only after global classification.  For a
 * sole source it therefore has complete object coverage and may replace the
 * structural overview before persistence finishes.  This is deliberately not
 * permission for the serialized spatial bootstrap: its seed page is local
 * cache work only and must remain behind structural coverage. */
static bool
controller_lod_global_preview_requested(size_t sourceMeshRequests)
{
    return sourceMeshRequests == 1;
}

static uint64_t
controller_lod_saturating_add_u64(uint64_t left, uint64_t right)
{
    return right > UINT64_MAX - left ? UINT64_MAX : left + right;
}

/* Decide only whether a complete scene may begin bounded mesh admission in
 * parallel with structural coverage.  Every source-backed mesh participates
 * in one aggregate profile so many individually small roots cannot bypass the
 * large-scene path.  This does not authorize a full-detail cut or reserve any
 * memory/render budget. */
static bool
controller_lod_mesh_first_scene_safe(
    const std::vector<SoBRLDatabaseSource *> &sources)
{
    uint64_t occurrenceCount = 0;
    uint64_t uniqueAssetCount = 0;
    uint64_t encodedSourceBytes = 0;
    uint64_t largestAssetBytes = 0;
    uint64_t meshRequestCount = 0;
    bool haveMeshSource = false;

    for (SoBRLDatabaseSource *source : sources) {
	if (!source)
	    continue;
	const size_t requests = source->getDisplayMeshLodRequestCount();
	if (!requests)
	    continue;
	BObolCompactSourceProfile profile;
	if (!source->getCompactSourceProfile(profile) || !profile.valid)
	    return false;
	haveMeshSource = true;
	occurrenceCount = controller_lod_saturating_add_u64(
	    occurrenceCount, profile.occurrenceCount);
	uniqueAssetCount = controller_lod_saturating_add_u64(
	    uniqueAssetCount, profile.uniqueAssetCount);
	encodedSourceBytes = controller_lod_saturating_add_u64(
	    encodedSourceBytes, profile.encodedSourceBytes);
	largestAssetBytes = std::max(
	    largestAssetBytes, profile.largestAssetBytes);
	meshRequestCount = controller_lod_saturating_add_u64(
	    meshRequestCount, static_cast<uint64_t>(requests));
    }
    return haveMeshSource &&
	BObolLodSourceProfilePolicy::safeMeshFirstPreview(
	    true, occurrenceCount, uniqueAssetCount, encodedSourceBytes,
	    largestAssetBytes, meshRequestCount);
}

/* Large-scene traces can execute thousands of bounded discovery passes before
 * reaching the camera epoch under investigation.  Let a diagnostic opt in to
 * a minimum view revision without changing the ordinary boolean trace knobs.
 * An absent or malformed threshold preserves their historical behavior. */
bool
controller_lod_trace_enabled(const char *name, uint64_t viewRevision)
{
    if (!name || !getenv(name))
	return false;
    const char *minimum = getenv("BOBOL_LOD_TRACE_MIN_VIEW_REVISION");
    if (!minimum || !minimum[0])
	return true;
    char *end = NULL;
    const unsigned long long parsed = strtoull(minimum, &end, 10);
    if (end == minimum || (end && end[0]))
	return true;
    return viewRevision >= static_cast<uint64_t>(parsed);
}
static double
controller_lod_visual_footprint(
    const BObolViewLodState::CadPayload *payload)
{
    if (!payload)
	return 0.0;
    return std::max(
	std::sqrt(std::max(0.0,
	    static_cast<double>(payload->projectedPixelArea))),
	std::max(
	    static_cast<double>(payload->projectedPixelPerimeter) * 0.25,
	    static_cast<double>(payload->projectedPixelDiameter) * 0.25));
}
static double
controller_lod_projected_error(
    const BObolViewLodState::CadPayload *payload, int cut)
{
    if (!payload || cut < 0 ||
	!std::isfinite(payload->projectedPixelDiameter) ||
	payload->projectedPixelDiameter <= 0.0f)
	return std::numeric_limits<double>::infinity();
    const double target = std::max(
	static_cast<double>(std::numeric_limits<float>::min()),
	static_cast<double>(payload->targetPixelError));
    if (!payload->progressiveMesh)
	return std::numeric_limits<double>::infinity();
    return payload->progressiveMesh->projectedErrorAtCut(
	cut, payload->projectedPixelDiameter) / target;
}

static int
controller_lod_log_bucket(double value)
{
    if (!std::isfinite(value))
	return 15;
    if (value <= 1.0)
	return 0;
    return std::max(0, std::min(15,
	static_cast<int>(std::floor(std::log2(value))) + 1));
}


static void
controller_lod_effective_population_cost(
    const BObolViewLodState *viewState, size_t &activeCost,
    size_t &minimumCost)
{
    activeCost = 0;
    minimumCost = 0;
    if (!viewState)
	return;

    /* CadPayload mutation maintains exact occurrence costs incrementally.
     * Reprojecting every retained spatial hierarchy merely to initialize a
     * new budget pass made an otherwise O(1) subpath draw pause for more than
     * 100 ms at 50k leaves, and source aliases could revisit the same retained
     * table.  A new camera epoch deliberately begins with the preceding
     * presentation as its continuity seed; the bounded demand census updates
     * each occurrence to the new visible page set afterward. */
    activeCost = viewState->activeRenderCost();
    minimumCost = viewState->minimumActiveRenderCost();

    /* The retained counters describe occurrence meshes, while the renderer
     * may have replaced many subpixel occurrences with one aggregate point
     * channel.  Prefer its last exact completed-frame CAD work record without
     * losing any non-CAD contribution.  This is a proven presentation cost,
     * not a prediction; live frame deadlines and later census results remain
     * authoritative if the new view costs more. */
    size_t presentedCadCost = 0;
    size_t presentedSubpixelOccurrences = 0;
    size_t presentedStructuralBoxes = 0;
    const SbBool haveCoverageClassification =
	viewState->lastCadPresentationOccurrenceCoverage(
	    presentedSubpixelOccurrences, presentedStructuralBoxes);
    if (viewState->lastCadPresentationFrameExact() &&
	haveCoverageClassification &&
	viewState->lastCadPresentedRenderCost(presentedCadCost)) {
	const size_t retainedCadCost = viewState->activeCadRenderCost();
	const size_t retainedCadMinimum =
	    viewState->minimumActiveCadRenderCost();
	const size_t nonCadCost = activeCost >= retainedCadCost ?
	    activeCost - retainedCadCost : 0;
	const size_t nonCadMinimum = minimumCost >= retainedCadMinimum ?
	    minimumCost - retainedCadMinimum : 0;
	activeCost = presentedCadCost > SIZE_MAX - nonCadCost ?
	    SIZE_MAX : nonCadCost + presentedCadCost;
	/* The exact presentation can contain structural boxes which have not yet
	 * acquired a mesh, and can omit retained meshes replaced by the aggregate
	 * point channel.  Neither retained counter is therefore a valid lower
	 * bound for this frame.  The completed work is a conservative temporary
	 * floor until the next allocation replaces a box or changes a mesh cut. */
	const size_t effectiveCadMinimum =
	    std::min(retainedCadMinimum, presentedCadCost);
	minimumCost = effectiveCadMinimum > SIZE_MAX - nonCadMinimum ?
	    SIZE_MAX : nonCadMinimum + effectiveCadMinimum;
    }
}

/* Resolve retained payload metadata through the compact index captured when
 * its request was submitted.  A matching population epoch authenticates the
 * cheap position directly.  After registry replacement the occurrence string
 * remains authoritative, so fall back to its keyed lookup rather than letting
 * a recycled position affect visibility or quality allocation. */
static bool
controller_lod_planning_summary(SoBRLDatabaseSource *source,
	const BObolViewLodState::CadPayload *payload,
	BObolCompactLodPlanningSummary &summary, size_t *entryIndex = NULL)
{
    summary = BObolCompactLodPlanningSummary();
    if (entryIndex)
	*entryIndex = 0;
    if (!source || !payload ||
	payload->sourceInstanceKey.getLength() == 0)
	return false;

    size_t index = static_cast<size_t>(payload->sourceEntryIndex);
    if (payload->sourceEntryIndex != UINT32_MAX &&
	payload->sourcePopulationEpoch ==
	    source->getCompactPopulationEpoch() &&
	index <= static_cast<size_t>(std::numeric_limits<int>::max()) &&
	source->getCompactLodPlanningSummary(static_cast<int>(index), summary) &&
	summary.valid) {
	if (entryIndex)
	    *entryIndex = index;
	return true;
    }

    if (!source->getCompactInstanceIndex(
	    payload->sourceInstanceKey.getString(), index) ||
	index > static_cast<size_t>(std::numeric_limits<int>::max()) ||
	!source->getCompactLodPlanningSummary(static_cast<int>(index), summary) ||
	!summary.valid)
	return false;
    if (entryIndex)
	*entryIndex = index;
    return true;
}

/* Build the large-scene retained frontier without projecting or sorting every
 * occurrence on the GUI thread.  requestedCut is the last exact projected
 * demand recorded for the payload; requested-active is therefore a
 * conservative logarithmic screen-error estimate.  A fixed bucket table
 * makes construction O(frontier) while ensuring severe quality deficits and
 * explicit user emphasis are visited before cheap, already-smooth details.
 * The bounded submit action recomputes exact projection for each consumed
 * wave, so this retained estimate is only an ordering seed. */
static void
controller_prioritize_lod_frontier(SoBRLDatabaseSource *source,
	const BObolViewLodState *viewState,
	const std::vector<SbString> &occurrenceKeys,
	std::vector<size_t> &entryIndices, bool soleSelectedOccurrence)
{
    static const size_t emphasisBuckets = 3;
    static const size_t errorBuckets = 16;
    static const size_t valueBuckets = 16;
    static const size_t bucketCount =
	emphasisBuckets * errorBuckets * valueBuckets;
    std::array<std::vector<size_t>, bucketCount> buckets;
    entryIndices.clear();
    entryIndices.reserve(occurrenceKeys.size());
    if (!source || !viewState)
	return;

    for (const SbString &key : occurrenceKeys) {
	size_t entryIndex = 0;
	if (!source->getCompactInstanceIndex(key.getString(), entryIndex))
	    continue;
	BObolCompactLodPlanningSummary summary;
	const bool haveSummary = source->getCompactLodPlanningSummary(
	    static_cast<int>(entryIndex), summary) && summary.valid;
	const BObolViewLodState::CadPayload *payload =
	    viewState->findCadForSourceEntry(source,
		static_cast<uint32_t>(entryIndex), key);
	if (!payload)
	    payload = viewState->findCadForOccurrence(source, key);
	const int emphasis = haveSummary && summary.selected &&
	    soleSelectedOccurrence ? 2 :
	    (haveSummary && summary.highlighted ? 1 : 0);
	const int activeCut = payload ? payload->activeCut : -1;
	const int requestedCut = payload ? payload->requestedCut : -1;
	const double retainedError = controller_lod_projected_error(
	    payload, activeCut);
	const int errorBucket = controller_lod_log_bucket(retainedError);

	/* Secondary marginal value uses retained hierarchy counts but no mesh
	 * scan.  The requested cut also supplies a bounded approximation of
	 * projected footprint for ties at the same worst-error tier. */
	double value = static_cast<double>(errorBucket + 1);
	if (payload && payload->progressiveMesh && activeCut >= 0 &&
	    requestedCut > activeCut) {
	    const size_t activeFaces = payload->progressiveMesh->
		hierarchyFaceCountAtCut(activeCut);
	    int nextCut = requestedCut;
	    for (int cut = activeCut + 1; cut <= requestedCut;
		 ++cut) {
		if (payload->progressiveMesh->hierarchyFaceCountAtCut(cut) >
			activeFaces) {
		    nextCut = cut;
		    break;
		}
	    }
	    const size_t nextFaces = payload->progressiveMesh->
		hierarchyFaceCountAtCut(nextCut);
	    const size_t addedFaces = nextFaces > activeFaces ?
		nextFaces - activeFaces : 1;
	    double footprint = controller_lod_visual_footprint(payload);
	    if (footprint <= 0.0)
		footprint = 1.0;
	    double currentError = controller_lod_projected_error(
		payload, activeCut);
	    double nextError = controller_lod_projected_error(
		payload, nextCut);
	    double errorReduction = 0.0;
	    if (std::isfinite(currentError) && std::isfinite(nextError))
		errorReduction = std::max(0.0, currentError - nextError);
	    else if (!std::isfinite(currentError) && std::isfinite(nextError))
		errorReduction = std::numeric_limits<double>::infinity();
	    value = footprint * errorReduction /
		static_cast<double>(addedFaces);
	}
	int valueBucket = 0;
	if (std::isfinite(value) && value > 0.0)
	    valueBucket = std::max(0, std::min(15,
		static_cast<int>(std::floor(std::log2(value))) + 8));
	else if (!std::isfinite(value))
	    valueBucket = 15;
	const size_t bucket =
	    (static_cast<size_t>(emphasis) * errorBuckets +
	     static_cast<size_t>(errorBucket)) * valueBuckets +
	    static_cast<size_t>(valueBucket);
	buckets[bucket].push_back(entryIndex);
    }
    for (size_t bucket = bucketCount; bucket > 0; --bucket)
	entryIndices.insert(entryIndices.end(),
	    buckets[bucket - 1].begin(), buckets[bucket - 1].end());
}

/* Order a retained-cut recovery by what must remain recognizable, not by
 * compact registry insertion order.  requestedCut is the last exact
 * screen-space demand is retained on each view payload, so both the current
 * projected error and the affected screen footprint are available without a
 * fresh scene projection.  Fixed buckets keep this O(N) and avoid the
 * comparison sort which made large-scene recovery a GUI-thread hotspot.
 * Exact projection is still recomputed by each bounded action window. */
static void
controller_prioritize_lod_recovery(SoBRLDatabaseSource *source,
	const BObolViewLodState *viewState,
	std::vector<size_t> &entryIndices, bool soleSelectedOccurrence)
{
    static const size_t emphasisBuckets = 3;
    static const size_t footprintBuckets = 16;
    static const size_t qualityBuckets = 16;
    static const size_t bucketCount =
	emphasisBuckets * footprintBuckets * qualityBuckets;
    std::array<std::vector<size_t>, bucketCount> buckets;
    entryIndices.clear();
    if (!source || !viewState)
	return;

    std::vector<const BObolViewLodState::CadPayload *> payloads;
    viewState->findCadPayloadsUnordered(source, payloads);
    entryIndices.reserve(payloads.size());
    for (const BObolViewLodState::CadPayload *payload : payloads) {
	if (!payload)
	    continue;
	size_t entryIndex = 0;
	BObolCompactLodPlanningSummary summary;
	const bool haveSummary = controller_lod_planning_summary(
	    source, payload, summary, &entryIndex) && summary.visible;
	if (!haveSummary)
	    continue;
	const int emphasis = summary.selected && soleSelectedOccurrence ? 2 :
	    (summary.highlighted ? 1 : 0);
	double visualFootprint = controller_lod_visual_footprint(payload);
	if (visualFootprint <= 0.0)
	    visualFootprint = 1.0;
	double projectedError = controller_lod_projected_error(
	    payload, payload->activeCut);
	const int footprint = controller_lod_log_bucket(visualFootprint);
	const int quality = controller_lod_log_bucket(projectedError);
	const size_t bucket =
	    (static_cast<size_t>(emphasis) * qualityBuckets +
	     static_cast<size_t>(quality)) * footprintBuckets +
	    static_cast<size_t>(footprint);
	buckets[bucket].push_back(entryIndex);
    }
    for (size_t bucket = bucketCount; bucket > 0; --bucket)
	entryIndices.insert(entryIndices.end(),
	    buckets[bucket - 1].begin(), buckets[bucket - 1].end());
}

static SbBool
controller_lod_trace_result(const BObolLodResult &result)
{
    const char *filter = getenv("BOBOL_LOD_TRACE_OBJECT");
    if (!filter || !filter[0])
	return FALSE;
    return (result.request.objectName.getLength() > 0 &&
	    strstr(result.request.objectName.getString(), filter)) ||
	   (result.request.objectPath.getLength() > 0 &&
	    strstr(result.request.objectPath.getString(), filter)) ?
	   TRUE : FALSE;
}

BObolProgressiveOptions::BObolProgressiveOptions(void) :
    /* maxLodResults caps results applied per frame pump.  It is now applied in
     * batched single-scene traversals (see processPendingLodResults) and still
     * bounded by maxLodApplyMicroseconds.  The inner atomic quantum is kept
     * small enough for that time budget to be meaningful; this outer ceiling
     * permits several cheap quanta in one pump. */
    maxLodResults(2048),
    maxLodApplyMicroseconds(4000),
    maxProviders(0),
    /* The microsecond limit is the responsiveness contract.  A 64-item
     * ceiling forced one 16 ms timer turn per batch even when merging a batch
     * took only microseconds, producing a visible many-thousand-leaf replay.
     * Keep a generous safety ceiling and let the measured 4 ms budget govern
     * ordinary pumping. */
    maxProviderItems(4096),
    /* A quiet, fast endpoint may use the whole 8 ms discovery slice.  The
     * provider pacing policy below restores the 4 ms input/software cap from
     * measured presentation state; explicit caller options remain exact. */
    maxProviderMicroseconds(8000),
    forceTerminalLodRefinement(FALSE)
{
}

BObolProgressiveStatus::BObolProgressiveStatus(void)
{
    this->clear();
}

BObolHostWorkSnapshot::BObolHostWorkSnapshot(void) :
    revision(0),
    renderRevision(0),
    flags(BOBOL_HOST_WORK_NONE)
{
}

SbBool
BObolHostWorkSnapshot::pumpPending(void) const
{
    return (this->flags & BOBOL_HOST_WORK_PUMP) ? TRUE : FALSE;
}

SbBool
BObolHostWorkSnapshot::renderPending(void) const
{
    return (this->flags & BOBOL_HOST_WORK_RENDER) ? TRUE : FALSE;
}

SbBool
BObolHostWorkSnapshot::capacitySampleRequested(void) const
{
    return (this->flags & BOBOL_HOST_WORK_CAPACITY_SAMPLE) ? TRUE : FALSE;
}

void
BObolProgressiveStatus::clear(void)
{
    this->providerCount = 0;
    this->providerAdvanced = 0;
    this->lodResultsProcessed = 0;
    this->lodResultsApplied = 0;
    this->submitted = 0;
    this->alreadyCached = 0;
    this->expanded = 0;
    this->existing = 0;
    this->remaining = 0;
    this->proxyPublished = 0;
    this->metadataApplied = 0;
    this->pendingTasks = 0;
    this->inFlight = 0;
    this->queuedResults = 0;
    this->queuedCacheWrites = 0;
    this->changed = 0;
    this->hasMore = 0;
}

BObolLodConvergenceStatus::BObolLodConvergenceStatus(void)
{
    this->clear();
}

void
BObolLodConvergenceStatus::clear(void)
{
    this->phase = BOBOL_LOD_CONVERGENCE_IDLE;
    this->outcome = BOBOL_LOD_PRESENTATION_READY;
    this->controlObligationMask = BOBOL_LOD_CONTROL_OBLIGATION_NONE;
    this->controlOwner = BOBOL_LOD_CONTROL_OWNER_NONE;
    this->controlViolationMask = BOBOL_LOD_CONTROL_VIOLATION_NONE;
    this->viewQualityHistoryEntryCount = 0;
    this->viewQualityHistoryRememberCount = 0;
    this->viewQualityHistoryRecallCount = 0;
    this->inventoryRevision = 0;
    this->availabilityRevision = 0;
    this->viewRevision = 0;
    this->policyRevision = 0;
    this->capacityRevision = 0;
    this->allocationPlanSerial = 0;
    this->presentationTransactionSerial = 0;
    this->presentationRequiredRenderSerial = 0;
    this->presentedFrameSerial = 0;
    this->activeGeneration = 0;
    this->submissionSourceIndex = 0;
    this->submissionEntryOffset = 0;
    this->expectedLeafCount = 0;
    this->availableLeafCount = 0;
    this->visibleTargetCount = 0;
    this->activePayloadCount = 0;
    this->satisfiedPayloadCount = 0;
    this->presentedSubpixelOccurrenceCount = 0;
    this->presentedStructuralBoxCount = 0;
    this->terminalOccurrenceFailureCount = 0;
    this->pendingTasks = 0;
    this->inFlight = 0;
    this->queuedResults = 0;
    this->queuedCacheWrites = 0;
    this->activeFaces = 0;
    this->activeRenderCost = 0;
    this->renderCostBudget = 0;
    this->selectedPresentationCost = 0;
    this->certifiedPresentationBudget = 0;
    this->pixelDemandPresentationCost = 0;
    this->requestedPresentationBudget = 0;
    this->maximumMarginalPresentationBudget = 0;
    this->maximumProtectedPresentationBudget = 0;
    this->pointProxyCandidateCount = 0;
    this->allocationCertificatePlanSerial = 0;
    this->residentMeshBytes = 0;
    this->stableResidentMeshBytes = 0;
    this->reservedResidentMeshGrowthBytes = 0;
    this->residentMeshLimitBytes = 0;
    this->memoryLimitedPayloadCount = 0;
    this->activeWorkingSetBytes = 0;
    this->peakWorkingSetBytes = 0;
    this->residentCompactionCount = 0;
    this->residentCompactionPlanRevision = 0;
    this->residentCompactionCandidateCount = 0;
    this->residentCompactionPlanCurrent = FALSE;
    this->gpuTrackedBufferBytes = 0;
    this->gpuOrdinaryPartBufferBytes = 0;
    this->gpuProgressiveCutBufferBytes = 0;
    this->gpuProgressiveActiveCutBytes = 0;
    this->gpuBatchBufferBytes = 0;
    this->gpuTriangleAtlasAllocatedBytes = 0;
    this->gpuTriangleAtlasLiveBytes = 0;
    this->gpuTriangleAtlasConfiguredCapacityBytes = 0;
    this->gpuTriangleAtlasPartCount = 0;
    this->gpuTriangleAtlasPageCount = 0;
    this->gpuOrdinaryPartFullUploadBytes = 0;
    this->gpuOrdinaryPartSuffixUploadBytes = 0;
    this->gpuOrdinaryPartGpuCopyBytes = 0;
    this->gpuOrdinaryPartLineageReuseCount = 0;
    this->gpuOrdinaryPartLineageReplacementCount = 0;
    this->gpuTriangleAtlasFullUploadBytes = 0;
    this->gpuTriangleAtlasSuffixUploadBytes = 0;
    this->gpuTriangleAtlasLineageReuseCount = 0;
    this->gpuPressureProxyCount = 0;
    this->gpuProgressiveEvictionCount = 0;
    this->gpuTriangleAtlasReclamationCount = 0;
    this->gpuResourceSampleSerial = 0;
    this->fraction = 0.0f;
    this->terminal = TRUE;
    this->terminalError = FALSE;
    this->viewReady = FALSE;
    this->hasLodState = FALSE;
    this->backgroundPending = FALSE;
    this->performanceLimited = FALSE;
    this->memoryLimited = FALSE;
    this->gpuMemoryPressure = FALSE;
    this->refinementFramePending = FALSE;
    this->budgetCalibrationPending = FALSE;
    this->stablePresentationHandoffPending = FALSE;
    this->pointProxyCalibrationPending = FALSE;
    this->pointProxyAdmissionFramePending = FALSE;
    this->stablePointProxyCalibrationPending = FALSE;
    this->pointProxyTriangleRecoveryPending = FALSE;
    this->residentGrowthReallocationPending = FALSE;
    this->publicationFramePending = FALSE;
    this->sourcePreparationPending = FALSE;
    this->sourcePreparationProviderCount = 0;
    this->failedSourceCount = 0;
}

BObolProgressiveProviderRecord::BObolProgressiveProviderRecord(void) :
    token(0),
    callback(NULL),
    userData(NULL),
    userDataFree(NULL)
{
}

static void
controller_append_database_source_unique(
    std::vector<SoBRLDatabaseSource *> &sources,
    SoBRLDatabaseSource *source)
{
    if (!source)
	return;

    const SbString sourceInstance = source->instanceKey.getValue();
    const SbString sourcePath = source->path.getValue();
    const int sourceDrawMode = source->drawMode.getValue();
    struct db_i *sourceDbip = source->getDatabase();

    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *existing = sources[i];
	if (!existing)
	    continue;
	if (existing == source)
	    return;

	const SbString existingInstance = existing->instanceKey.getValue();
	if (sourceInstance.getLength() > 0 &&
	    existingInstance.getLength() > 0 &&
	    bu_strcmp(sourceInstance.getString(),
		   existingInstance.getString()) == 0)
	    return;

	if (sourcePath.getLength() > 0 &&
	    existing->path.getValue().getLength() > 0 &&
	    sourceDbip == existing->getDatabase() &&
	    sourceDrawMode == existing->drawMode.getValue() &&
	    bu_strcmp(sourcePath.getString(),
		   existing->path.getValue().getString()) == 0)
	    return;
    }

    sources.push_back(source);
}

static void
controller_collect_database_sources(SoNode *node,
				    std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
	if (std::find(sources.begin(), sources.end(), source) == sources.end())
	    sources.push_back(source);
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	controller_collect_database_sources(group->getChild(i), sources);
}

static void
controller_collect_database_source_roots(
    SoNode *node,
    std::vector<SoBRLDatabaseSource *> &sources)
{
    if (!node)
	return;

    if (node->isOfType(SoBRLDatabaseSource::getClassTypeId())) {
	SoBRLDatabaseSource *source = static_cast<SoBRLDatabaseSource *>(node);
	controller_append_database_source_unique(sources, source);
	return;
    }

    if (!node->isOfType(SoGroup::getClassTypeId()))
	return;

    SoGroup *group = static_cast<SoGroup *>(node);
    for (int i = 0; i < group->getNumChildren(); i++)
	controller_collect_database_source_roots(group->getChild(i), sources);
}

std::vector<SoBRLDatabaseSource *>
controller_render_database_source_roots(
    const BObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (!controller)
	return sources;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();
    controller_collect_database_source_roots(root, sources);
    if (!sources.empty())
	return sources;

    const int sourceCount = controller->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	controller_append_database_source_unique(sources, source);
    }

    return sources;
}

/* A compact producer publishes its final expected population before all leaf
 * batches reach the owner thread.  Exact inventory deltas are sufficient
 * while that append-only stream is growing; rescanning the complete prefix
 * after every delta makes cold discovery quadratic in the number of batches.
 * The final equality is the producer witness which permits one authoritative
 * all-source coverage pass. */
static bool
controller_lod_compact_inventory_incomplete(
    const std::vector<SoBRLDatabaseSource *> &sources)
{
    for (SoBRLDatabaseSource *source : sources) {
	if (!source)
	    continue;
	const size_t expected = source->getCompactExpectedInstanceCount();
	if (!expected)
	    continue;
	const int currentCount = source->getCompactInstanceCount();
	const size_t current = currentCount > 0 ?
	    static_cast<size_t>(currentCount) : 0;
	if (current < expected)
	    return true;
    }
    return false;
}

std::vector<SoBRLDatabaseSource *>
controller_render_database_sources(const BObolViewController *controller)
{
    std::vector<SoBRLDatabaseSource *> sources;
    if (!controller)
	return sources;

    SoNode *root = controller->getRenderSceneRoot();
    if (!root)
	root = controller->getSceneRoot();
    controller_collect_database_sources(root, sources);
    if (!sources.empty())
	return sources;

    const int sourceCount = controller->getDatabaseSourceCount();
    for (int i = 0; i < sourceCount; i++) {
	SoBRLDatabaseSource *source = controller->getDatabaseSource(i);
	if (source && std::find(sources.begin(), sources.end(), source) ==
	    sources.end())
	    sources.push_back(source);
    }

    return sources;
}

BObolLodPresentationPolicy::Population
controller_lod_presentation_population(const BObolViewLodState *state,
	uint64_t sceneDomainRevision)
{
    BObolLodPresentationPolicy::Population population;
    population.available = state != NULL;
    population.sceneDomainRevision = sceneDomainRevision;
    return population;
}

static int
size_to_int_saturated(size_t value)
{
    return value > static_cast<size_t>(std::numeric_limits<int>::max()) ?
	   std::numeric_limits<int>::max() : static_cast<int>(value);
}

static void
append_controller_lod_diagnostic(SbString &diagnostics,
				 const SbString &target,
				 const char *message)
{
    if (diagnostics.getLength() > 0)
	diagnostics += "\n";

    diagnostics += target.getLength() > 0 ? target : SbString("<unknown>");
    diagnostics += ": ";
    diagnostics += message ? message : "LoD controller diagnostic";
}

const char *
controller_database_id(const struct db_i *dbip)
{
    return dbip && dbip->dbi_filename ? dbip->dbi_filename : "";
}

struct controller_lod_view_signatures {
    BObolLodViewSnapshot view;
    BObolLodViewScaleSnapshot scale;
};

static controller_lod_view_signatures
controller_lod_view_signature(const struct bv_view_info &view,
			      SbBool haveCamera,
			      SoCamera *camera,
			      const SbViewportRegion &region)
{
    controller_lod_view_signatures signatures;
    signatures.view.haveCamera = haveCamera ? 1 : 0;
    signatures.view.width = view.width;
    signatures.view.height = view.height;
    signatures.view.size = view.size;
    signatures.view.lodScale = view.lod.scale;
    signatures.view.curveScale = view.lod.curve_scale;
    signatures.view.pointScale = view.lod.point_scale;
    signatures.view.botThreshold = view.lod.bot_threshold;

    signatures.scale.haveCamera = haveCamera ? 1 : 0;
    signatures.scale.width = view.width;
    signatures.scale.height = view.height;
    signatures.scale.size = view.size;
    signatures.scale.lodScale = view.lod.scale;
    signatures.scale.curveScale = view.lod.curve_scale;
    signatures.scale.pointScale = view.lod.point_scale;
    signatures.scale.botThreshold = view.lod.bot_threshold;
    const SbVec2s viewport = region.getViewportSizePixels();
    signatures.scale.viewportWidth = viewport[0];
    signatures.scale.viewportHeight = viewport[1];

    if (haveCamera && camera) {
	double aspect = controller_aspect_from_region(region);
	if (aspect <= SMALL_FASTF)
	    aspect = 1.0;
	const SbMatrix matrix = camera->getViewVolume(
	    static_cast<float>(aspect)).getMatrix();
	const float *values = matrix[0];
	for (size_t i = 0; i < 16; i++)
	    signatures.view.viewVolumeMatrix[i] = values[i];
	signatures.scale.cameraTypeKey =
	    static_cast<uint64_t>(camera->getTypeId().getKey());
	signatures.scale.aspectRatio = camera->aspectRatio.getValue();
	signatures.scale.focalDistance = camera->focalDistance.getValue();
	if (camera->isOfType(SoOrthographicCamera::getClassTypeId()))
	    signatures.scale.projectionScale =
		static_cast<SoOrthographicCamera *>(camera)->
		    height.getValue();
	else if (camera->isOfType(SoPerspectiveCamera::getClassTypeId()))
	    signatures.scale.projectionScale =
		static_cast<SoPerspectiveCamera *>(camera)->
		    heightAngle.getValue();
    }

    return signatures;
}

struct controller_lod_source_signatures {
    std::vector<BObolLodCoordinator::LodSourceSnapshot> sources;

    bool empty(void) const
    {
	return this->sources.empty();
    }

    bool sameIdentities(
	const std::vector<BObolLodCoordinator::LodSourceSnapshot>
	    &other) const
    {
	if (this->sources.size() != other.size())
	    return false;
	for (size_t i = 0; i < this->sources.size(); ++i) {
	    if (!this->sources[i].sameIdentity(other[i]))
		return false;
	}
	return true;
    }

    bool sameInventories(
	const std::vector<BObolLodCoordinator::LodSourceSnapshot>
	    &other) const
    {
	if (!this->sameIdentities(other))
	    return false;
	for (size_t i = 0; i < this->sources.size(); ++i) {
	    if (this->sources[i].inventoryRevision !=
		    other[i].inventoryRevision)
		return false;
	}
	return true;
    }
};

/*
 * Keep source identity separate from the source's growing display-LoD
 * inventory.  Compact leaf coverage is append-only while a cold/warm
 * realization stream is active: each merge advances displayMeshLodRevision,
 * but it does not replace the source or invalidate the already-consumed
 * prefix of a pinned submission plan.
 *
 * Conflating the two made every arriving box/request batch look like source
 * replacement.  The controller consequently restarted at compact entry zero
 * on every GUI pump and could not begin mesh work until discovery was nearly
 * complete.  Inventory changes instead request/extend a bounded rescan; only
 * an identity/policy change resets the scan.
 */
static controller_lod_source_signatures
controller_lod_source_signature(const BObolViewController *controller)
{
    controller_lod_source_signatures signatures;

    if (!controller)
	return signatures;

    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_sources(controller);
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase())
	    continue;
	const SbBool hasCurrentCompactRequests =
	    source->hasDisplayMeshLodRequests();
	const SbBool hasCurrentRealizedMesh =
	    source->realizationStatus.getValue() ==
		SoBRLDatabaseSource::REALIZED &&
	    !source->needsRealization() &&
	    source->hasRealizedMeshGeometry();
	if (getenv("BOBOL_LOD_TRACE_SOURCE_CONTRACT")) {
	    static std::atomic<unsigned int> traceCount(0);
	    const unsigned int traceIndex = traceCount.fetch_add(1);
	    if (traceIndex < 64)
		bu_log("BObol LoD source contract path=%s compact=%d "
		       "entries=%d requests_current=%d realized_mesh_current=%d "
		       "status=%d needs_realization=%d source_revision=%u "
		       "inputs_revision=%u lod_revision=%llu\n",
		       source->path.getValue().getString(),
		       source->hasCompactInstanceIndex() ? 1 : 0,
		       source->getCompactInstanceCount(),
		       hasCurrentCompactRequests ? 1 : 0,
		       hasCurrentRealizedMesh ? 1 : 0,
		       source->realizationStatus.getValue(),
		       source->needsRealization() ? 1 : 0,
		       source->sourceRevision.getValue(),
		       source->inputsRevision.getValue(),
		       static_cast<unsigned long long>(
			   source->getDisplayMeshLodRevision()));
	}
	/*
	 * A compact request is published only after the stream merge validates
	 * the detached/live source epoch, and hasDisplayMeshLodRequests()
	 * rechecks that epoch.  It is therefore ready for view planning even
	 * while the source-wide authoritative walk is still running.  Requiring
	 * REALIZED here serialized cold mesh work behind complete leaf
	 * discovery and defeated progressive drawing.
	 */
	if (!hasCurrentCompactRequests && !hasCurrentRealizedMesh)
	    continue;

	BObolLodCoordinator::LodSourceSnapshot snapshot;
	snapshot.source = source;
	snapshot.database = source->getDatabase();
	snapshot.routingId.set(source->getCompactSourceRoutingId());
	snapshot.inventoryRevision.set(source->getDisplayMeshLodRevision());
	snapshot.databaseId = controller_database_id(source->getDatabase());
	snapshot.path = source->path.getValue();
	snapshot.drawMode = source->drawMode.getValue();
	snapshot.representationMode = source->representationMode.getValue();
	snapshot.visible = source->visible.getValue();
	snapshot.lodBotThreshold = source->lodBotThreshold.getValue();
	snapshot.sourceRevision = source->sourceRevision.getValue();
	snapshot.inputsRevision = source->inputsRevision.getValue();
	signatures.sources.push_back(std::move(snapshot));
    }

    return signatures;
}

void
BObolViewController::lodResultReadyCB(
    BObolLodService *service, void *userData)
{
    BObolViewController *controller =
	static_cast<BObolViewController *>(userData);
    if (!controller || !service)
	return;

    const uint64_t generation = controller->d->lodActiveGeneration;
    const uint64_t consumerId =
	controller->d->residentMeshConsumerId();
    const bool generationReady = generation != 0 &&
	service->queuedResultCountForGeneration(generation) > 0;
    const bool compactionReady =
	service->queuedResidentMeshCompactionResultCountForDiagnostics(
	    consumerId) > 0;
    if (!generationReady && !compactionReady)
	return;

    if (generationReady) {
	controller->d->lodAvailabilityLedger.noteResultsReady(bu_gettime());
    }
    controller->markProgressiveWorkPending();
}

void
BObolViewController::setLodService(BObolLodService *service)
{
    if (this->d->lodService == service)
	return;

    this->d->resetLodViewQualityHistory();
    this->cancelActiveLodGeneration();
    if (this->d->lodService) {
	this->d->lodService->releaseResidentMeshConsumer(
	    this->d->residentMeshConsumerId());
	if (this->d->lodResultSubscriberId != 0)
	    this->d->lodService->unsubscribeResultReady(
		this->d->lodResultSubscriberId);
    }

    this->d->lodService = service;
    struct bv_lod_policy serviceLodPolicy;
    bv_lod_policy_init(&serviceLodPolicy);
    this->d->viewAttachment->getLodPolicy(&serviceLodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	service && this->d->lodAutoSubmit &&
	serviceLodPolicy.policy != BV_LOD_OFF &&
	serviceLodPolicy.mesh_enabled);
    this->d->lodResultSubscriberId = 0;
    this->d->lodAvailabilityLedger.resetResultQueue();
    this->d->lodActiveGeneration = 0;
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPass.deactivate();
    this->d->lodSubmissionPass.clearRescan();
    this->d->lodSubmissionIntent.reset();
    this->d->lodViewDemandPolicy.reset();
    this->d->lodInteractionStartCertificate.reset();
    this->d->lodPoseContinuity.reset();
    this->d->lodPlanningObligations.retireImportanceCensus();
    this->d->lodAvailabilityLedger.resetResidentGrowth();
    this->d->lodPlanningObligations.retireResidentAdmissionRetry();
    this->d->lodDiscretePopulationTrialPermit.revoke();
    this->d->lodStaticQualityTrial.reset();
    this->d->lodCoveragePolicy.reset();
    if (service)
	this->d->lodCoveragePolicy.activate(true);
    this->d->resetRetainedPassAnnotations();
    this->d->lodStructuralRepair.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR);
    this->d->advanceAdmissionRevision(
        BObolLodAdmissionRevisionDomain::CAPACITY);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
    this->d->lodInteractionSession.clearMotionFrameGate();
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->lodLastSubmittedSources.clear();
    this->d->lodSubmissionDelta.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_RESOURCE_SERVICE);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    this->d->lodResidentAdmissionRevision =
	service ? service->residentMeshAdmissionRevision() : 0;
    this->d->lodCompactionPolicy.resetForServiceChange(
	service != NULL, bu_gettime(), 750000);
    this->d->lodPresentationTransaction.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodInterruptedPresentationReplay.reset();

    if (this->d->lodService)
	this->d->lodResultSubscriberId =
	    this->d->lodService->subscribeResultReady(
		BObolViewController::lodResultReadyCB, this);
}

uint64_t
BObolViewController::beginLodGeneration(void)
{
    if (!this->d->lodService || !this->d->lodService->isRunning())
	return 0;
    this->cancelActiveLodGeneration();
    this->d->lodActiveGeneration =
	this->d->lodService->beginGeneration();
    this->d->lodAvailabilityLedger.resetResultQueue();
    return this->d->lodActiveGeneration;
}

void
BObolViewController::resetDiscoveryPointProxyFloor(SbBool requestFrame)
{
    const float oldEffective = std::max(
	this->d->lodPresentationPointProxyPixelThreshold,
	this->d->lodDiscoveryPointProxyPixelThreshold);
    this->d->lodDiscoveryPointProxyPixelThreshold = 1.0f;
    this->d->lodPointAdmissionFrame.reset();
    BObolViewLodState *viewState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    if (viewState)
	viewState->setCadPresentationDiscoveryPointProxyPixelThreshold(1.0f);
    if (requestFrame && this->d->lodAutoSubmit &&
	this->d->lodPresentationPointProxyPixelThreshold + 1.0e-6f <
	    oldEffective) {
	this->requestRender("lod-discovery-point-reset");
    }
}

void
BObolViewController::cancelActiveLodGeneration(void)
{
    if (this->d->lodService && this->d->lodActiveGeneration != 0) {
	this->d->lodService->cancelGeneration(this->d->lodActiveGeneration);
    }
    this->d->lodActiveGeneration = 0;
    this->d->lodSubmissionSourceIndex = 0;
    this->d->lodSubmissionEntryOffset = 0;
    this->d->clearLodSubmissionPlan();
    this->d->lodSubmissionPass.deactivate();
    this->d->lodSubmissionPass.clearRescan();
    this->d->lodSubmissionIntent.reset();
    this->d->lodViewDemandPolicy.reset();
    this->d->lodDiscretePopulationTrialPermit.revoke();
    this->d->lodAvailabilityLedger.resetResidentGrowth();
    this->d->lodPlanningObligations.retireResidentAdmissionRetry();
    this->d->resetRetainedPassAnnotations();
    this->d->lodStructuralRepair.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR);
    this->d->advanceAdmissionRevision(
        BObolLodAdmissionRevisionDomain::CAPACITY);
    this->d->lodInteractionSession.clearMotionFrameGate();
    this->d->lodAvailabilityLedger.resetResultQueue();
    this->d->lodLastSubmittedViewRevision.reset();
    this->d->lodLastSubmittedPolicyRevision.reset();
    this->d->lodLastSubmittedSources.clear();
    this->d->lodSubmissionDelta.reset();
    this->d->lodPresentationTransaction.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodInterruptedPresentationReplay.reset();
    this->d->lodInteractiveProgressiveCeiling = -1;
    this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_STRUCTURAL_ADMISSION);
    this->d->lodPointQualityPhase.reset();
    this->resetDiscoveryPointProxyFloor(FALSE);
    this->d->lodPresentationPolicy.reset();
    if (this->d->viewAttachment &&
	this->d->viewAttachment->getViewLodState())
	{
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    viewState->setCadPresentationProgressiveCutCeiling(-1);
	    viewState->setCadPresentationPointProxyPixelThreshold(1.0f);
	    viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
}

void
BObolViewController::invalidateDatabaseSourceLodState(void)
{
    this->d->resetLodViewQualityHistory();
    this->cancelActiveLodGeneration();
    /* Source replacement retires every routing/index identity used by the
     * completed view census and projection cache.  cancelActiveLodGeneration
     * deliberately preserves those camera facts for ordinary generation
     * changes, but it also clears the submitted-source snapshots.  Without
     * clearing the source-keyed evidence here, the replacement source is
     * admitted as an initial contract and its visible count is added to the
     * retired source's count (for example, six occurrences become twelve).
     * Begin one fresh authoritative coverage pass for the new inventory. */
    this->d->clearLodConvergenceCandidates();
    this->d->lodProjectedDemandCache.clear();
    this->d->lodCoveragePolicy.activate(true);
    this->d->resetLodConvergenceFraction();
    if (this->d->viewAttachment)
	this->d->viewAttachment->clearViewLodState();
}

BObolLodService *
BObolViewController::getLodService(void) const
{
    return this->d->lodService;
}

BObolLodService *
BObolViewController::ensureManagedLodService(size_t workerCount)
{
    if (workerCount == 0)
	workerCount = 1;
    if (!this->d->managedLodService) {
	this->d->managedLodService =
	    controller_acquire_managed_lod_service(workerCount);
	if (!this->d->managedLodService)
	    return NULL;
    } else if (!this->d->managedLodService->isRunning()) {
	if (!this->d->managedLodService->start(workerCount, TRUE))
	    return NULL;
    } else if (!this->d->managedLodService->ensureWorkerCount(workerCount)) {
	return NULL;
    }
    this->d->managedLodWorkerCount =
	this->d->managedLodService->workerCountForDiagnostics();
    this->setLodService(this->d->managedLodService.get());
    return this->d->managedLodService.get();
}

void
BObolViewController::stopManagedLodService(void)
{
    this->setLodAutoSubmit(FALSE);
    if (this->d->lodService == this->d->managedLodService.get())
	this->setLodService(NULL);
    this->d->managedLodService.reset();
    this->d->managedLodWorkerCount = 0;
}

size_t
BObolViewController::getManagedLodWorkerCount(void) const
{
    if (this->d->managedLodService)
	return this->d->managedLodService->workerCountForDiagnostics();
    return this->d->managedLodWorkerCount;
}

void
BObolViewController::setLodAutoSubmit(SbBool enabled)
{
    const SbBool requested = enabled ? TRUE : FALSE;
    if (this->d->lodAutoSubmit != requested)
	this->d->resetLodViewQualityHistory();
    this->d->lodAutoSubmit = requested;
    this->d->lodStaticQualityTrial.reset();
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    if (this->d->lodAutoSubmit) {
	this->requestRender("lod-auto-submit");
    } else {
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
	this->d->lodAvailabilityLedger.resetResidentGrowth();
	this->d->lodPlanningObligations.retireResidentAdmissionRetry();
	this->d->lodInteractiveProgressiveCeiling = -1;
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_STRUCTURAL_ADMISSION);
	this->d->lodPointQualityPhase.reset();
	this->resetDiscoveryPointProxyFloor(FALSE);
	this->d->lodPresentationPolicy.reset();
	this->d->lodViewDemandPolicy.reset();
	this->d->lodDiscretePopulationTrialPermit.revoke();
	if (this->d->viewAttachment->getViewLodState()) {
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    viewState->setCadPresentationProgressiveCutCeiling(-1);
	    viewState->setCadPresentationPointProxyPixelThreshold(1.0f);
	    viewState->setCadPresentationCameraMotionFrameReuse(FALSE);
	}
    }
}

SbBool
BObolViewController::isLodAutoSubmitEnabled(void) const
{
    return this->d->lodAutoSubmit;
}

void
BObolViewController::setLodForcedCut(int cut)
{
    if (cut < 0)
	cut = 0;

    if (this->d->lodForcedCut && *this->d->lodForcedCut == cut)
	return;

    this->d->resetLodViewQualityHistory();
    this->d->lodForcedCut = cut;
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

void
BObolViewController::clearLodForcedCut(void)
{
    if (!this->d->lodForcedCut)
	return;

    this->d->resetLodViewQualityHistory();
    this->d->lodForcedCut.reset();
    this->advanceLodPolicyRevision();
    this->requestRender("lod-policy");
}

SbBool
BObolViewController::hasLodForcedCut(void) const
{
    return this->d->lodForcedCut.has_value() ? TRUE : FALSE;
}

int
BObolViewController::getLodForcedCut(void) const
{
    return this->d->lodForcedCut.value_or(0);
}

void
BObolViewController::setExactFullDetailBudget(uint64_t maxFaceCount,
						uint64_t maxPointCount)
{
    this->d->maxExactFullDetailFaceCount = maxFaceCount;
    this->d->maxExactFullDetailPointCount = maxPointCount;
}

uint64_t
BObolViewController::getMaxExactFullDetailFaceCount(void) const
{
    return this->d->maxExactFullDetailFaceCount;
}

uint64_t
BObolViewController::getMaxExactFullDetailPointCount(void) const
{
    return this->d->maxExactFullDetailPointCount;
}


SbBool
BObolViewController::hasPendingLodResults(void) const
{
    return this->d->lodAvailabilityLedger.resultsPending() != 0 ? TRUE : FALSE;
}

SbBool
BObolViewController::hasPendingLodSubmissions(void) const
{
    return this->d->lodSubmissionPass.active();
}

SbBool
BObolViewController::hasPendingLodRefinementFrame(void) const
{
    /* A budget calibration probe is just as much a refinement barrier as a
     * newly selected PoP cut.  Reporting idle in the interval between
     * requesting that frame and completeRenderTiming() made warm scenes stop
     * at cache-temperature-dependent coarse cuts. */
    return this->d->lodPresentationTransaction.barrierPending() ||
	this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() ||
	(this->d->lodPresentationPolicy.handoffPending() &&
	 !this->d->lodSubmissionPass.active()) ||
	this->d->lodPointAdmissionFrame.pending() ||
	this->d->lodPointQualityPhase.presentationPending() ||
	this->d->lodAdmissionEvidence.headroom().retryPending();
}

size_t
BObolViewController::processPendingLodResults(size_t maxResults,
	uint64_t maxMicroseconds)
{
    if (!this->d->lodService)
	return 0;

    const auto lod_wakeup_required = [this]() {
	return this->d->lodSubmissionPass.active() ||
	    this->d->lodAvailabilityLedger.resultsPending() != 0 ||
	    this->d->lodAvailabilityLedger.residentGrowthPending() ||
	    this->hasPendingLodRefinementFrame() ||
	    (this->d->lodService &&
	     this->d->lodService->queuedResultCountForGeneration(
		 this->d->lodActiveGeneration) > 0);
    };
    const auto clear_lod_wakeup_if_idle = [this, &lod_wakeup_required]() {
	if (lod_wakeup_required()) {
	    /* A presentation barrier is work even after its last producer has
	     * drained.  In particular, an interrupted software traversal consumes
	     * the render request before it knows whether the frame is complete.
	     * Preserve the pump level until advanceProgressiveWork() can replace
	     * that producer witness with the barrier's explicit successor frame. */
	    this->markProgressiveWorkPending();
	    return;
	}

	/* Merely having a registered provider does not mean that provider has
	 * work.  Keeping the latch set for the lifetime of every GED provider
	 * made a late result-ready callback leave an otherwise idle warm scene
	 * polling forever.  Provider callbacks report their own hasMore state
	 * from advanceProgressiveWork().
	 *
	 * A result-ready callback may race this drain.  Clear first, then recheck
	 * so a concurrent callback cannot lose its frame wakeup. */
	this->clearProgressiveWorkPending();
	if (lod_wakeup_required())
	    this->markProgressiveWorkPending();
    };

    if (!this->hasPendingLodResults() &&
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) == 0) {
	clear_lod_wakeup_if_idle();
	return 0;
    }

    if (maxMicroseconds == 0) {
	(void)this->applyLodResults(this->d->lodService, maxResults);
	clear_lod_wakeup_if_idle();
	return this->d->lastLodResultCount;
    }

    /*
     * Apply drained results in batches to amortize queue synchronization and
     * payload-journal publication.  The time check happens only between these
     * atomic calls, so a 2048-result quiet quantum made the nominal 4 ms
     * contract ineffective (50k GUI traces measured individual applications
     * above 35 ms).  Keep a 256-result atomic quantum, but honor the outer
     * time contract between quanta.  Immutable sparse part publication and
     * adaptive frame batching now bound the deferred renderer update; forcing
     * one timer turn per quantum left a full 2048-result queue idle while the
     * workers were back-pressured.
     *
     * A first useful mesh still yields after one quantum so the user sees it
     * before later generations.  A refinement barrier whose publication
     * deadline is still pending instead belongs to the current coalescing
     * transaction: stopping there would publish one result per frame and
     * serially starve a large result queue.  Explicit callers which pass a
     * zero time budget retain the unbounded drain path above.
     */
    static const size_t apply_batch = 256;
    static const size_t apply_byte_budget = 8u * 1024u * 1024u;
    size_t processed = 0;
    unsigned int matched = 0;
    unsigned int applied = 0;
    unsigned int rejected = 0;
    unsigned int unmatched = 0;
    SbString diagnostics;
    const int64_t started = bu_gettime();
    const bool firstUsefulCoverage =
	this->getActiveLodMeshPayloadCount() == 0;
    while (maxResults == 0 || processed < maxResults) {
	size_t batch = apply_batch;
	if (maxResults != 0) {
	    const size_t remaining = maxResults - processed;
	    if (remaining < batch)
		batch = remaining;
	}
	(void)this->applyLodResults(
	    this->d->lodService, batch, apply_byte_budget);
	if (this->d->lastLodResultCount == 0)
	    break;
	processed += this->d->lastLodResultCount;
	matched += this->d->lastLodMatchedResultCount;
	applied += this->d->lastLodAppliedResultCount;
	rejected += this->d->lastLodRejectedResultCount;
	unmatched += this->d->lastLodUnmatchedResultCount;
	if (this->d->lastLodDiagnostics.getLength() > 0) {
	    if (diagnostics.getLength() > 0)
		diagnostics += "\n";
	    diagnostics += this->d->lastLodDiagnostics;
	}
	const bool awaitingPublicationDeadline =
	    this->d->lodPresentationTransaction.barrierPending() &&
	    this->d->lodPresentationTransaction.publicationAwaitingFrameRequest();
	if (firstUsefulCoverage ||
	    (this->d->lodPresentationTransaction.barrierPending() &&
	     !awaitingPublicationDeadline))
	    break;
	const int64_t now = bu_gettime();
	if (now < started ||
	    static_cast<uint64_t>(now - started) >= maxMicroseconds)
	    break;
    }
    this->d->lastLodResultCount = processed;
    this->d->lastLodMatchedResultCount = matched;
    this->d->lastLodAppliedResultCount = applied;
    this->d->lastLodRejectedResultCount = rejected;
    this->d->lastLodUnmatchedResultCount = unmatched;
    this->d->lastLodDiagnostics = diagnostics;
	clear_lod_wakeup_if_idle();
    return processed;
}

int
BObolViewController::submitLodRequestsIfNeeded(SbBool refreshMissing,
	SbBool resetExisting)
{
    if (!this->d->lodService || !this->d->lodService->isRunning())
	return 0;
    if (!this->d->activeCamera ||
	(!this->getSceneRoot() && !this->getRenderSceneRoot()))
	return 0;

    this->syncLodViewSignature(TRUE);

    const controller_lod_source_signatures signatures =
	controller_lod_source_signature(this);
    if (signatures.empty()) {
	/*
	 * No current source can consume a submission cursor.  This is a valid
	 * transient while a source is replaced, and a terminal state when the
	 * scene has no LoD-backed meshes.  Leaving the old cursor armed creates
	 * a zero-progress GUI pump: submitLodRequestsIfNeeded() returns here on
	 * every turn while hasPendingLodSubmissions() prevents convergence.
	 *
	 * Forget the old source epoch as well.  If a contract reappears after
	 * streaming or editing, even with the same identity text, it must start
	 * one fresh authoritative pass.
	 */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.deactivate();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionDelta.reset();
	/* The empty-contract edge retires the source identity domain, not merely
	 * the current submission cursor.  A Z/redraw may pass through this state
	 * before publishing a replacement source; preserving the old dense census
	 * then adds the replacement population to an unreachable source entry. */
	this->d->lodCoveragePolicy.reset();
	this->d->clearLodConvergenceCandidates();
	this->d->lodProjectedDemandCache.clear();
	this->d->lodAvailabilityLedger.resetResidentGrowth();
	this->d->lodPlanningObligations.retireResidentAdmissionRetry();
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	this->d->lodLastSubmittedSources.clear();
	return 0;
    }

    if (this->d->lodLastSubmittedViewRevision == this->d->lodViewRevision &&
	this->d->lodLastSubmittedPolicyRevision == this->d->lodPolicyRevision &&
	signatures.sameInventories(this->d->lodLastSubmittedSources)) {
	if (this->d->lodSubmissionPass.active()) {
	    /* A pending cursor is not meaningful without an owned generation.
	     * Scene replacement normally clears both together, but a late
	     * resident-growth/coverage edge can legitimately re-arm submission
	     * after cancellation.  The unchanged-signature fast path formerly
	     * returned forever in that state: no task, result, frame, or camera
	     * event remained to create another generation. */
	    if (this->d->lodActiveGeneration == 0) {
		this->d->lodActiveGeneration =
		    this->d->lodService->beginGeneration();
		this->d->lodAvailabilityLedger.resetResultQueue();
		if (this->d->lodActiveGeneration == 0)
		    return 0;
	    }
	    return this->submitLodRequests(this->d->lodService,
		this->d->lodActiveGeneration,
		this->d->lodSubmissionIntent.refreshMissing() ? TRUE : FALSE,
		this->d->lodSubmissionIntent.resetExisting() ? TRUE : FALSE);
	}
	return 0;
    }

    /* A view or LoD-policy epoch does not invalidate source geometry.  Keep
     * useful cold-cache work alive and submit only newly demanded geometry
     * into the same cancellation domain.  In particular, do not restart an
     * in-progress compact-index scan at entry zero on every mouse event.
     * Doing so starves high-index leaves during view interaction and leaves
     * isolated structural boxes on screen.  Finish the current scan using
     * the newest view, then make one complete current-view rescan so entries
     * visited before the change also reach the stable pixel target.  Source
     * replacement paths call cancelActiveLodGeneration explicitly. */
    uint64_t generation = this->d->lodActiveGeneration;
    if (generation == 0)
	generation = this->d->lodService->beginGeneration();
    const bool sourceSetChanged =
	!this->d->lodLastSubmittedSources.empty() &&
	!signatures.sameIdentities(this->d->lodLastSubmittedSources);
    const bool inventoryChanged =
	!this->d->lodLastSubmittedSources.empty() &&
	!signatures.sameInventories(this->d->lodLastSubmittedSources);
    const bool viewOrPolicyChanged =
	this->d->lodLastSubmittedViewRevision != this->d->lodViewRevision ||
	this->d->lodLastSubmittedPolicyRevision !=
	    this->d->lodPolicyRevision;
    /* A contract appearing after the explicit empty-source state has no
     * predecessor against which sourceSetChanged can compare.  Re-arm the
     * authoritative coverage pass here; initial service setup is already
     * active, so this is idempotent before the first observation. */
    if (this->d->lodLastSubmittedSources.empty() &&
	this->d->lodCoveragePolicy.required() &&
	!this->d->lodCoveragePolicy.active())
	this->d->lodCoveragePolicy.activate(true);
    if (sourceSetChanged || inventoryChanged || viewOrPolicyChanged)
	this->d->lodAvailabilityLedger.interruptResidencyDrain();
    /* The first source contract is submitted immediately.  Once that
     * time-to-first-mesh edge exists, allow an append-only producer to build
     * a bounded inventory wave before starting another owner-thread LoD
     * transaction.  All semantic invalidations and the producer's final edge
     * bypass this gate.  Because lodLastSubmittedSources remains unchanged
     * while deferred, the source journal supplies the complete accumulated
     * delta when the deadline expires. */
    const bool deferInventoryDelta =
	!sourceSetChanged && !viewOrPolicyChanged && inventoryChanged &&
	this->d->lodAvailabilityLedger.deferInventoryDelta(
	    true, this->d->lodAvailabilityLedger.providerPendingCount() > 0,
	    this->d->lodSubmissionPass.active() != FALSE,
	    !this->d->lodLastSubmittedSources.empty(),
	    this->d->lodInteractionSession.active() != FALSE, bu_gettime());
    if (deferInventoryDelta) {
	this->markProgressiveWorkPending();
	return 0;
    }
    this->d->lodAvailabilityLedger.commitInventoryDelta();
    if (sourceSetChanged || inventoryChanged || viewOrPolicyChanged)
	this->d->lodPlanningObligations.retireResidentAdmissionRetry();
    bool useSourceDelta = false;
    bool extendedPendingDelta = false;
    bool pendingDeltaNeedsFullRescan = false;
    bool sourceDeltaInvalidatesCoverage = false;
    const bool priorCoverageComplete =
	this->d->lodCoveragePolicy.coverageComplete();
    size_t sourceDeltaFirst = signatures.sources.size();
    /*
     * A source can append another structural batch while an exact delta plan
     * is only partly consumed.  Extend that pinned plan in place.  Replacing
     * it restarts at zero; merely setting rescanPending while leaving delta
     * mode active repeats the old delta forever and never visits the new
     * tail.
     *
     * Keep appended indices at the tail rather than sorting the whole plan:
     * lodSubmissionEntryOffset is a cursor into this exact ordering.  A
     * repeated changed index is harmless and necessary when an entry already
     * consumed by this pass is upgraded again.
     */
    if ((sourceSetChanged || inventoryChanged) &&
	!viewOrPolicyChanged &&
	this->d->lodSubmissionPass.active() &&
	this->d->lodSubmissionDelta.active() &&
	!this->d->lodLastSubmittedSources.empty()) {
	for (size_t currentIndex = 0;
	     currentIndex < signatures.sources.size(); ++currentIndex) {
	    const BObolLodSourceRoutingId routingId =
		signatures.sources[currentIndex].routingId;
	    size_t previousIndex =
		this->d->lodLastSubmittedSources.size();
	    for (size_t candidate = 0;
		 candidate <
		    this->d->lodLastSubmittedSources.size();
		 ++candidate) {
		if (this->d->lodLastSubmittedSources[candidate].routingId ==
			routingId) {
		    previousIndex = candidate;
		    break;
		}
	    }
	    const bool knownSource =
		previousIndex <
		    this->d->lodLastSubmittedSources.size();
	    const bool sameIdentity = knownSource &&
		signatures.sources[currentIndex].sameIdentity(
		    this->d->lodLastSubmittedSources[previousIndex]);
	    const uint64_t previousInventory =
		knownSource ?
		    this->d->lodLastSubmittedSources[previousIndex].
			inventoryRevision.value() : 0;
	    if (sameIdentity && previousInventory ==
		    signatures.sources[currentIndex].
			inventoryRevision.value())
		continue;
	    if (!sameIdentity || !previousInventory) {
		pendingDeltaNeedsFullRescan = true;
		continue;
	    }

	    SoBRLDatabaseSource *changedSource =
		signatures.sources[currentIndex].source;
	    std::vector<size_t> changedEntries;
	    SbBool coverageInvalidated = FALSE;
	    if (!changedSource->getDisplayMeshLodChangedEntries(
		    previousInventory, changedEntries,
		    &coverageInvalidated) ||
		changedEntries.empty()) {
		pendingDeltaNeedsFullRescan = true;
		continue;
	    }
	    if (coverageInvalidated) {
		sourceDeltaInvalidatesCoverage = true;
		pendingDeltaNeedsFullRescan = true;
	    }
	    std::vector<size_t> *existingEntries =
		this->d->lodSubmissionDelta.selectiveEntries(changedSource);
	    if (!existingEntries) {
		/* No selective plan means this source is already being scanned
		 * in full.  Its identity plan is extended from compactCount in
		 * submitLodRequests below. */
		const bool alreadyTargeted =
		    this->d->lodSubmissionDelta.targets(changedSource);
		if (!alreadyTargeted)
		    pendingDeltaNeedsFullRescan = true;
		else
		    extendedPendingDelta = true;
		continue;
	    }
	    existingEntries->insert(existingEntries->end(),
		changedEntries.begin(), changedEntries.end());
	    if (this->d->lodSubmissionSourcePlan.validFor(changedSource))
		this->d->lodSubmissionSourcePlan.entries().insert(
		    this->d->lodSubmissionSourcePlan.entries().end(),
		    changedEntries.begin(), changedEntries.end());
	    extendedPendingDelta = true;
	}
    }
    if ((sourceSetChanged || inventoryChanged) &&
	!viewOrPolicyChanged &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodLastSubmittedSources.empty()) {
	for (size_t currentIndex = 0;
	     currentIndex < signatures.sources.size(); currentIndex++) {
	    const BObolLodSourceRoutingId routingId =
		signatures.sources[currentIndex].routingId;
	    size_t previousIndex =
		this->d->lodLastSubmittedSources.size();
	    for (size_t candidate = 0;
		 candidate <
		    this->d->lodLastSubmittedSources.size();
		 candidate++) {
		if (this->d->lodLastSubmittedSources[candidate].routingId ==
			routingId) {
		    previousIndex = candidate;
		    break;
		}
	    }
	    const bool knownSource =
		previousIndex <
		    this->d->lodLastSubmittedSources.size();
	    const bool sameIdentity = knownSource &&
		signatures.sources[currentIndex].sameIdentity(
		    this->d->lodLastSubmittedSources[previousIndex]);
	    const uint64_t previousInventory =
		knownSource ?
		    this->d->lodLastSubmittedSources[previousIndex].
			inventoryRevision.value() : 0;
	    const bool inventoryDiffers =
		!knownSource ||
		previousInventory !=
		    signatures.sources[currentIndex].
			inventoryRevision.value();
	    if (sameIdentity && !inventoryDiffers)
		continue;

	    SoBRLDatabaseSource *changedSource =
		signatures.sources[currentIndex].source;
	    const bool alreadyTargeted =
		this->d->lodSubmissionDelta.targets(changedSource);
	    sourceDeltaFirst = std::min(sourceDeltaFirst, currentIndex);
	    if (sameIdentity && previousInventory) {
		std::vector<size_t> changedEntries;
		SbBool coverageInvalidated = FALSE;
		if (changedSource->getDisplayMeshLodChangedEntries(
			previousInventory, changedEntries,
			&coverageInvalidated) &&
		    !changedEntries.empty()) {
		    if (coverageInvalidated)
			sourceDeltaInvalidatesCoverage = true;
		    std::vector<size_t> *existingEntries =
			this->d->lodSubmissionDelta.selectiveEntries(
			    changedSource);
		    /* An existing target without a delta plan already requires
		     * a full source scan; a later selective update must not
		     * accidentally narrow it. */
		    if (existingEntries) {
			existingEntries->insert(existingEntries->end(),
			    changedEntries.begin(), changedEntries.end());
			std::sort(existingEntries->begin(),
			    existingEntries->end());
			existingEntries->erase(std::unique(
			    existingEntries->begin(), existingEntries->end()),
			    existingEntries->end());
		    } else if (!alreadyTargeted) {
			(void)this->d->lodSubmissionDelta.targetSelective(
			    changedSource, std::move(changedEntries));
		    }
		} else {
		    sourceDeltaInvalidatesCoverage = true;
		    this->d->lodSubmissionDelta.targetFull(changedSource);
		}
	    } else {
		this->d->lodSubmissionDelta.targetFull(changedSource);
	    }
	}
	useSourceDelta = this->d->lodSubmissionDelta.active();
    } else if ((sourceSetChanged || inventoryChanged) &&
	!this->d->lodSubmissionPass.active()) {
	this->d->lodSubmissionDelta.reset();
    }
    const bool hasExactSourceDelta =
	useSourceDelta || this->d->lodSubmissionDelta.active();
    const bool sourceCoverageInvalidated = sourceSetChanged ||
	(inventoryChanged &&
	 (!hasExactSourceDelta || sourceDeltaInvalidatesCoverage ||
	  pendingDeltaNeedsFullRescan || !priorCoverageComplete));
    if (sourceSetChanged || inventoryChanged) {
	if (sourceSetChanged)
	    this->d->lodProjectedDemandCache.clear();
	this->d->resetLodViewQualityHistory();
	/* A population mutation invalidates the complete-scene denominator.  An
	 * exact visibility/edit delta does not: scanning its changed entries is
	 * sufficient, and re-proving every unchanged occurrence makes a one-path
	 * redraw O(scene size). */
	if (sourceCoverageInvalidated)
	    this->d->lodCoveragePolicy.activate(true);
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	this->d->lodPoseContinuity.clearRetainOccurrenceCuts();
	this->d->lodPlanningObligations.retireImportanceCensus();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR);
	/* A structural append invalidates the old coverage proof, but it does not
	 * invalidate an already observed richer immutable prefix.  Carry that
	 * allocation edge through the new complete coverage pass; clearing it here
	 * made cold convergence depend on whether the last source-inventory batch
	 * raced the last resident-result batch. */
	/* A new occurrence population invalidates both sides of the measured
	 * small-part bracket.  Do not, however, discard the obligation to relax a
	 * coarse presentation-only threshold installed while that population was
	 * streaming.  The final unchanged replay must prove the finest sustainable
	 * cut before convergence may become stable. */
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	const SbBool pointRelaxationRequired =
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f ?
		TRUE : FALSE;
	if (this->d->lodPointQualityPhase.triangleRecoveryPending()) {
	    /* Structural streaming extended the population during the bounded
	     * recovery.  Keep its terminal obligation and its evidence provenance
	     * when restarting the complete-source pass.  Converting a one-pass
	     * handoff normalization into a persistent measured ceiling freezes a
	     * warm cache at whichever tiny prefix happened to arrive first. */
	    if (this->d->lodAdmissionEvidence.capacity().retainedRecoveryCeilingActive())
		this->d->requestRetainedRecovery(
		    this->d->lodAdmissionEvidence.capacity().currentBudget());
	    else
		this->d->requestRetainedNormalization(
		    this->d->lodAdmissionEvidence.capacity().currentBudget());
	} else {
	    if (pointRelaxationRequired)
		this->d->lodPointQualityPhase.requestCalibration();
	    else
		this->d->lodPointQualityPhase.completeCalibration();
	}
	/* Preserve the completed per-entry visibility census for an exact
	 * visibility/edit delta.  The bounded action below updates only its changed
	 * entries and republishes the denominator in O(delta) time.  Population or
	 * identity changes still discard the proof and perform one authoritative
	 * all-source census. */
	if (sourceCoverageInvalidated)
	    this->d->clearLodConvergenceCandidates();
	this->d->resetLodConvergenceFraction();
    }

    /* Source callbacks and owner-thread scene publication may be coalesced.
     * Reconcile the dense convergence census with the authoritative source
     * signature after all ordinary invalidation handling, so an old draw-mode
     * source can never survive merely because no empty scene was observed. */
    this->d->lodConvergenceCandidateCensus.beginSourceSetUpdate();
    bool convergenceSourceDomainChanged = false;
    for (const BObolLodCoordinator::LodSourceSnapshot &source :
	    signatures.sources) {
	if (this->d->lodConvergenceCandidateCensus.retainSource(
		this->d->lodConvergenceSourceKey(source.source)))
	    convergenceSourceDomainChanged = true;
    }
    if (this->d->lodConvergenceCandidateCensus.endSourceSetUpdate())
	convergenceSourceDomainChanged = true;
    if (convergenceSourceDomainChanged) {
	if (this->d->lodCoveragePolicy.hasCompleteVisibleCount())
	    this->d->lodCoveragePolicy.setCompleteVisibleCount(
		this->d->lodConvergenceCandidateCensusTotal());
	this->d->lodCoveragePolicy.activate(true);
	this->d->resetLodConvergenceFraction();
    }
    if (sourceSetChanged || useSourceDelta ||
	!this->d->lodSubmissionPass.active()) {
	this->d->lodSubmissionSourceIndex =
	    useSourceDelta ? sourceDeltaFirst : 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	/*
	 * A population-changing delta is a latency optimization, not a proof that
	 * the complete current view is covered.  In particular, a cold compact
	 * source may publish several large, overlapping append deltas while its
	 * structural registry is still growing.  Declaring coverage from only the
	 * last append made completion timing-dependent: untouched ranges remained
	 * structural boxes until the next camera change.
	 *
	 * Consume such an exact delta first so newly arrived leaves become useful,
	 * then perform one bounded all-source pass.  Pure visibility or in-place
	 * request changes retain the existing complete coverage proof and finish
	 * after their exact entry plan.
	 */
	this->d->lodSubmissionPass.setRescanPending(
	    useSourceDelta && sourceCoverageInvalidated);
	this->d->resetRetainedPassAnnotations();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::INVENTORY);
	if (useSourceDelta)
	    this->d->lodCoveragePolicy.clearPassCounters();
	if (sourceSetChanged || inventoryChanged)
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
    } else if (!extendedPendingDelta) {
	this->d->lodSubmissionPass.requestRescan();
    }
    if (pendingDeltaNeedsFullRescan)
	this->d->lodSubmissionPass.requestRescan();
    this->d->lodSubmissionPass.activate();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    this->d->lodSubmissionIntent.configure(
	refreshMissing != FALSE, resetExisting != FALSE);
    int submitted = this->submitLodRequests(this->d->lodService, generation,
					    refreshMissing, resetExisting);
    if (submitted >= 0) {
	this->d->lodActiveGeneration = generation;
	this->d->lodLastSubmittedViewRevision = this->d->lodViewRevision;
	this->d->lodLastSubmittedPolicyRevision = this->d->lodPolicyRevision;
	this->d->lodLastSubmittedSources = signatures.sources;
    }
    return submitted;
}

int
BObolViewController::submitLodRequests(BObolLodService *service,
	uint64_t generation,
	SbBool refreshMissing,
	SbBool resetExisting)
{
    if (!service)
	service = this->d->lodService;

    this->d->lastLodVisitedMeshCount = 0;
    this->d->lastLodSubmittedTaskCount = 0;
    this->d->lastLodUpdatedCutCount = 0;
    this->d->lastLodSkippedMeshCount = 0;
    this->d->lastLodDiagnostics = "";

    if (!service || !service->isRunning()) {
	this->d->lastLodDiagnostics = "LoD service is not running";
	return -1;
    }
    if (service != this->d->lodService) {
	this->d->lastLodDiagnostics =
	    "LoD submission service is not owned by this controller";
	return -1;
    }
    if (generation == 0) {
	this->d->lastLodDiagnostics =
	    "LoD submission requires an owned generation";
	return -1;
    }
    if (this->d->lodActiveGeneration != 0 &&
	this->d->lodActiveGeneration != generation) {
	this->d->lastLodDiagnostics =
	    "LoD submission generation is not owned by this controller";
	return -1;
    }
    this->d->lodActiveGeneration = generation;

    struct bv_view_info view = BV_VIEW_INFO_INIT;
    if (!this->getViewInfo(&view)) {
	this->d->lastLodDiagnostics = "LoD submission requires an active camera";
	return -1;
    }

    if (!this->getSceneRoot() && !this->getRenderSceneRoot()) {
	this->d->lastLodDiagnostics = "LoD submission requires a scene root";
	return -1;
    }

    /* Source-signature reconciliation in submitLodRequestsIfNeeded() runs
     * before this gate.  A real source/inventory change therefore retires an
     * obsolete sample first, while an unchanged population remains frozen for
     * its presentation-owned measurement. */

    const BObolViewLodState *submissionPresentationState =
	this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
    if (BObolLodAdmissionPlanner::presentationPausesSubmission(
	    this->d->lodPointAdmissionFrame.pending(),
	    this->d->lodPointQualityPhase.presentationPending(),
	    this->d->lodAdmissionEvidence.capacity().rescanAfterFrame(),
	    submissionPresentationState &&
		submissionPresentationState->hasCadPresentationAssemblies())) {
	this->d->lastLodDiagnostics =
	    "LoD submission waits for a presentation measurement";
	return 0;
    }

    if (!this->d->lodSubmissionPass.active()) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->d->lodCoveragePolicy.clearPassCounters();
	this->d->lodSubmissionPass.activate();
	this->d->lodSubmissionIntent.configure(
	    refreshMissing != FALSE, resetExisting != FALSE);
    }

    /* Convert observed aggregate throughput into one total scene population
     * budget.  Per-object projected error remains the quality demand, but it
     * cannot be the admission policy: ten thousand individually modest cuts
     * can collectively be unusable.  The seed makes an all-box scene bounded
     * before any calibration sample exists; later frames learn this exact
     * renderer/draw-mode/viewport's capacity.
     *
     * Quiet views use a lower, still finite FPS target.  Pixel exactness is
     * reached whenever it fits that budget; otherwise the stable state is the
     * richest view-prioritized cut the machine can sustain.  Offline and
     * deterministic callers explicitly request the unbounded policy. */
    std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(this);
    const bool meshFirstSceneSafe =
	controller_lod_mesh_first_scene_safe(sources);
    /* Selection is normally a presentation-only transaction.  The one
     * intentional exception is a sole effective occurrence, which may become
     * a primitive-edit target and needs full Coin interaction geometry.  Use
     * the compact sources' retained selection populations and saturate at
     * two; this remains O(source count), not O(selected leaves). */
    size_t selectedOccurrenceCount = 0;
    for (SoBRLDatabaseSource *selectionSource : sources) {
	const int sourceSelected = selectionSource ?
	    selectionSource->getCompactSelectedInstanceCount() : 0;
	if (sourceSelected <= 0)
	    continue;
	if (sourceSelected > 1 || selectedOccurrenceCount != 0) {
	    selectedOccurrenceCount = 2;
	    break;
	}
	selectedOccurrenceCount = 1;
    }
    double lodAspect = controller_aspect_from_region(this->d->viewportRegion);
    if (lodAspect <= SMALL_FASTF)
	lodAspect = 1.0;
    const SbViewVolume lodViewVolume =
	this->d->activeCamera->getViewVolume(static_cast<float>(lodAspect));
    const SbMatrix lodViewProjection = lodViewVolume.getMatrix();
    /* See the allocation gate below.  Keep this witness outside the
     * admission-plan block: one retained transaction spans many bounded source
     * windows, and its completion/handoff uses the same producer state. */
    const bool retainedProviderInventorySettled =
	this->d->lodAvailabilityLedger.providerPendingCount() == 0 &&
	!controller_lod_compact_inventory_incomplete(sources);
    const bool retainedServiceStreamIdle = !service ||
	(service->activeTaskCountForGeneration(generation) == 0 &&
	 service->queuedResultCountForGeneration(generation) == 0);
    const bool retainedResultDeliveryIdle =
	this->d->lodAvailabilityLedger.resultsPending() == 0 &&
	!this->d->lodPresentationTransaction.publicationPending();
    const bool retainedPopulationSettled =
	BObolLodAvailabilityScheduler::allocationPopulationSettled(
	    retainedProviderInventorySettled, retainedServiceStreamIdle,
	    retainedResultDeliveryIdle,
	    this->d->lodAvailabilityLedger.residentGrowthPending(),
	    this->d->lodAvailabilityLedger.residencyDrainActive());
    const size_t presentationReconciliationBudget =
	this->d->lodPresentationPolicy.handoffReconciliationBudget();
    if (!this->d->lodAdmissionCursor.initialized()) {
	if (presentationReconciliationBudget)
	    this->d->requestPresentationReconciliation(
		presentationReconciliationBudget);
	this->d->lodRetainedPass.clearMissingMeshBudgetBlocked();
	/* A retained minimax allocation is one complete occurrence-plan
	 * transaction, even when result publication or an endpoint deadline
	 * inserts a presentation between bounded planning slices.  Those frames
	 * reset the per-frame allowance, but they must not silently turn the
	 * unconsumed tail of the plan into an ordinary first-come refinement
	 * pass.  Re-arm the same measured-budget allocation until completedPass
	 * retires the mode below. */
	if (this->d->lodSubmissionIntent.retainedAdmission() &&
	    this->d->lodSubmissionPass.active())
	    this->d->requestRetainedReallocation();
	const BObolViewLodState *lodState =
	    this->d->viewAttachment->getViewLodState();
	const size_t activeFaces = lodState ? lodState->activeFaceCount() : 0;
	size_t activeCost = 0;
	size_t minimumActiveCost = 0;
	controller_lod_effective_population_cost(lodState, activeCost,
	    minimumActiveCost);
	const SbBool scaleQualityProbe =
	    this->d->lodInteractionSession.active() &&
	    this->d->lodViewDemandPolicy.qualityBudgetActive() ? TRUE : FALSE;
	/* A zoom quality frame deliberately has a lower cadence target than the
	 * ordinary motion frame: coherent PoP populations are discrete, and the
	 * next useful cut can be just beyond a strict 20/60 Hz allowance.  Use the
	 * same 10 Hz hard floor which decides whether that cut remains visible.
	 * Exact hierarchy costs still prevent an unexpectedly huge cut jump. */
	/* The global progressive-ordinal staircase is valid only for one visible
	 * occurrence, but the static frame allowance is not single-object policy.
	 * A many-part scene spends that allowance through the occurrence-local
	 * importance allocator.  Keeping its target at the ordinary 20 Hz cadence
	 * made an accepted static phase recompute the identical constrained plan
	 * indefinitely instead of using (or rejecting) its bounded 10 Hz budget. */
	const SbBool hardDeadlinePresentation =
	    !this->d->lodInteractionSession.active() &&
	    (this->d->lodStructuralRepair.active() ||
	     this->d->lodStaticQualityTrial.blocksNewTrial() ||
	     this->d->lodPoseContinuity.retainOccurrenceCuts()) ? TRUE : FALSE;
	const float targetFps =
	    scaleQualityProbe ?
		BObolLodViewDemandPolicy::
		    qualityTargetFramesPerSecond() :
		(this->d->lodInteractionSession.active() ?
		    this->d->lodInteractiveTargetFps :
		    (hardDeadlinePresentation ?
			this->d->staticQualityTargetFps() :
			this->d->quietAllocationTargetFps()));
	const long double calibratedCostPerSecond =
	    scaleQualityProbe ?
		this->d->lodStableCalibratedRenderCostPerSecond :
		(this->d->lodInteractionSession.active() ?
		    this->d->lodInteractiveCalibratedRenderCostPerSecond :
		    this->d->lodStableCalibratedRenderCostPerSecond);
    /* The asynchronous GPU query describes the exact work record captured
     * with that query, not necessarily sceneActiveCost after a later result
     * or cut update.  It already contributes a correctly paired
     * cost-per-second sample in completeRenderTiming().  Pairing its duration
     * with the current population here manufactured false overloads and made
     * large scenes oscillate between coarse and rich cuts. */
    const uint64_t observedStableNanoseconds =
	this->d->lodLastRenderWasReusableCadPresentation ?
	    this->d->lastRenderTimeNanoseconds : 0;
	BObolLodCapacityEvidence::Inputs budgetInputs;
	budgetInputs.activeCost = activeCost;
	budgetInputs.minimumActiveCost = minimumActiveCost;
	budgetInputs.targetFps = targetFps;
	budgetInputs.calibratedCostPerSecond = calibratedCostPerSecond;
	budgetInputs.observedStableNanoseconds = observedStableNanoseconds;
	budgetInputs.hardPresentationDeadlineNanoseconds =
	    hardDeadlinePresentation ?
		this->d->staticQualityPresentationDeadline() :
		this->d->stablePresentationFrameDeadlineNanoseconds;
	budgetInputs.lastRenderNanoseconds =
	    this->d->lastRenderTimeNanoseconds;
	budgetInputs.smoothedRenderNanoseconds =
	    this->d->smoothedRenderTimeNanoseconds;
	budgetInputs.interactive = this->d->lodInteractionSession.active() != FALSE;
	budgetInputs.scaleQualityProbe = scaleQualityProbe != FALSE;
	size_t sourceMeshRequestCount = 0;
	size_t sourceLogicalOccurrenceCount = 0;
	for (SoBRLDatabaseSource *source : sources) {
	    if (!source)
		continue;
	    const size_t count = source->getDisplayMeshLodRequestCount();
	    sourceMeshRequestCount = count > SIZE_MAX - sourceMeshRequestCount ?
		SIZE_MAX : sourceMeshRequestCount + count;
	    const int compactCount = source->getCompactInstanceCount();
	    const size_t logicalCount = compactCount > 0 ?
		static_cast<size_t>(compactCount) : count;
	    sourceLogicalOccurrenceCount = logicalCount >
		SIZE_MAX - sourceLogicalOccurrenceCount ? SIZE_MAX :
		sourceLogicalOccurrenceCount + logicalCount;
	}
	this->d->observeLodSourceLogicalOccurrenceCount(
	    sourceLogicalOccurrenceCount);
	budgetInputs.coldSingleOccurrence =
	    activeCost == 0 &&
	    !controller_lod_compact_inventory_incomplete(sources) &&
	    (this->d->lodConvergenceCandidateCount() == 1 ||
	     sourceMeshRequestCount == 1);
	budgetInputs.hardDeadlinePresentation =
	    hardDeadlinePresentation != FALSE;
	budgetInputs.preserveActivePopulation =
	    this->d->lodPoseContinuity.retainOccurrenceCuts();
	budgetInputs.structuralCoverageRepair =
	    this->d->lodStructuralRepair.active();
	budgetInputs.forceTerminal =
	    this->d->forceTerminalLodRefinement != FALSE;
    budgetInputs.releaseCutFloor =
	    this->d->lodInteractionSession.releaseCutFloorActive();
	budgetInputs.stablePresentationHandoff =
		    this->d->lodPresentationPolicy.handoffPending();
    budgetInputs.stablePresentationCostFloor = std::max(
	this->d->lodPresentationPolicy.handoffCostFloor(),
	this->d->lodStaticQualityTrial.acceptedPresentationCostFor(
	    this->d->admissionRevisionStamp()));
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
	    this->d->lodViewRevision.value())) {
	const BObolLodCapacitySearchCertificate &search =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	bu_log("BObol LoD admission input phase=%u candidate=%zu "
	       "samples_remaining=%u current_budget=%zu cursor=%d "
	       "cursor_refinement=%zu cursor_retained=%d "
	       "target_fps=%.3f hard_deadline=%d pose_continuity=%d\n",
	       static_cast<unsigned int>(search.phase()),
	       search.candidateBudget(), search.samplesRemaining(),
	       this->d->lodAdmissionEvidence.capacity().currentBudget(),
	       this->d->lodAdmissionCursor.initialized() ? 1 : 0,
	       this->d->lodAdmissionCursor.refinementRemaining(),
	       this->d->lodAdmissionCursor.retainedAdmission() ? 1 : 0,
	       targetFps, hardDeadlinePresentation ? 1 : 0,
	       this->d->lodPoseContinuity.retainOccurrenceCuts() ? 1 : 0);
    }
    const BObolLodAdmissionPlan admissionPlan =
	BObolLodAdmissionPlanner::plan(
	    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
	    this->d->admissionRevisionStamp(), budgetInputs);
    this->d->commitAdmissionPlan(admissionPlan);
    const BObolLodCapacityEvidence::Decision &budget =
	admissionPlan.capacityDecision;
    if (this->d->lodStructuralRepair.active())
	this->d->lodStructuralRepair.reserveCoverageCost(
	    BObolLodAdmissionPlanner::structuralPerOccurrenceReservation(
		budget.totalBudget, this->d->lodAdmissionCursor.activeCost(),
		this->d->lodStructuralRepair.frontierCount()));
    /* An append-only database producer already publishes useful structural
     * and minimum-mesh deltas.  Its momentarily idle mesh-service queue is
     * not a complete scene population: running the global importance
     * allocator after every append batch serializes discovery behind repeated
     * O(scene) work and makes the HUD alternate BALANCING/REFINING.  Continue
     * the bounded delta/coverage path while the producer is active, then
     * allocate once on its terminal inventory edge. */
	const bool beginRetainedAdmission =
	    budget.retainedAdmission && retainedPopulationSettled;
	/* A sparse unsatisfied-refinement plan is not a retained-recovery plan.
	 * Reusing it lets its first subset consume the complete upgrade allowance;
	 * the later all-occurrence pass then sees zero and normalizes everything
	 * else to minimum.  Give the mode transition an explicit plan epoch and
	 * restart at the first source. */
	if (beginRetainedAdmission &&
	    !this->d->lodSubmissionIntent.retainedAdmission()) {
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionIntent.setRetainedAdmission(true);
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	} else if (!beginRetainedAdmission &&
	    this->d->lodSubmissionIntent.retainedAdmission()) {
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionIntent.setRetainedAdmission(false);
	    this->d->resetRetainedAdmissionQualityProof();
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	}
	if (getenv("BOBOL_LOD_TRACE_BUDGET"))
	    bu_log("BObol LoD frame budget active_faces=%zu active_cost=%zu "
		   "minimum_cost=%zu "
		   "last_render_ms=%.3f smooth_render_ms=%.3f target_fps=%.3f "
		   "calibrated_mcost_s=%.3f total_cost=%zu "
		   "additional_cost=%zu%s\n",
		   activeFaces, activeCost, minimumActiveCost,
		   this->d->lastRenderTimeNanoseconds / 1000000.0,
		   this->d->smoothedRenderTimeNanoseconds / 1000000.0,
		   targetFps,
		   static_cast<double>(
		       calibratedCostPerSecond / 1000000.0L),
		   budget.totalBudget, budget.refinementBudget,
		   budget.totalBudget == SIZE_MAX ? " unbounded" : "");
    }

    /* Population settlement gates the start of a retained allocation above.
     * Once started, that allocation owns one complete bounded traversal.  A
     * suffix becoming resident during the traversal queues a successor; it
     * must not reinterpret the unconsumed tail as an ordinary first-come
     * pass. */
    const bool effectiveRetainedAdmission =
	this->d->lodAdmissionCursor.retainedAdmission() &&
	this->d->lodSubmissionIntent.retainedAdmission();

    /* The ordinary quiet path reaches the raster-stable tier through measured
     * 1.0 -> 0.75/0.5 -> 0.25 pixel admissions.  A deterministic terminal
     * caller has deliberately removed those capacity ceilings, so select the
     * same final physical target directly.  This remains view-dependent and
     * does not imply loading the hierarchy's full topology. */
    const BObolViewLodState *sceneLodState =
	this->d->viewAttachment->getViewLodState();
    const size_t sceneActiveFaces = sceneLodState ?
	sceneLodState->activeFaceCount() : 0;
    /* One admission plan owns a complete population census.  Reuse its frozen
     * currencies across every bounded compact window; rescanning all retained
     * occurrences here made a 150k scene spend roughly three quarters of its
     * CPU in point-proxy cost discovery instead of publishing geometry. */
    const size_t sceneActiveCost =
	this->d->lodAdmissionCursor.activeCost();
    const size_t sceneMinimumActiveCost =
	this->d->lodAdmissionCursor.minimumActiveCost();
    const SbBool retainedAllocationPass =
	effectiveRetainedAdmission &&
	!this->d->lodInteractionSession.active() ? TRUE : FALSE;
    float scenePixelError = this->d->forceTerminalLodRefinement ?
	std::min(this->d->lodTargetPixelError, 0.25f) :
	this->d->lodTargetPixelError;
    if (!this->d->lodForcedCut &&
	!retainedAllocationPass &&
	this->d->lodAdmissionEvidence.capacity().currentBudget() != SIZE_MAX &&
	this->d->lodAdmissionEvidence.capacity().currentBudget() > 0 &&
	sceneActiveCost > this->d->lodAdmissionEvidence.capacity().currentBudget()) {
	const long double over =
	    static_cast<long double>(sceneActiveCost) /
	    static_cast<long double>(this->d->lodAdmissionEvidence.capacity().currentBudget());
	scenePixelError *= static_cast<float>(std::sqrt(over));
    }

	/* Active motion uses the renderer's O(1) aggregate ceiling and leaves
     * occurrence cuts intact while zoom prefetches missing resident suffixes.
     * A complete minimax reallocation is stable-view policy: running it in
     * the first wheel callback blocked input for hundreds of milliseconds at
     * 50k leaves and duplicated work the quiet demand census must perform
     * against the final camera anyway. */
    if (retainedAllocationPass &&
	this->d->lodSubmissionSourceIndex == 0 &&
	this->d->lodSubmissionEntryOffset == 0) {
	size_t maximumProtectedBudget =
	    this->d->lodAdmissionEvidence.capacity().currentBudget();
	size_t maximumMarginalBudget = maximumProtectedBudget;
	if (!this->d->lodInteractionSession.active() &&
	    !this->d->lodStaticQualityTrial.capacityRejected() &&
	!this->d->lodAdmissionEvidence.capacity().retainedQualityFloorRejected() &&
	this->d->lodStableCalibratedRenderCostPerSecond > 0.0L) {
	    const long double hardQualityFps =
		static_cast<long double>(
		    this->d->prominentQualityTargetFps());
	    const long double affordable = hardQualityFps > 0.0L ?
		this->d->lodStableCalibratedRenderCostPerSecond * 0.80L /
		    hardQualityFps : 0.0L;
	    const size_t hardQualityBudget =
		affordable >= static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
		(affordable > 0.0L ? static_cast<size_t>(affordable) : 0);
	    maximumProtectedBudget = std::max(
		maximumProtectedBudget, hardQualityBudget);
	}
	/* The atomic protected floor may be larger than the measured static-frame
	 * allowance even though many valuable marginal improvements fit.  During
	 * the bounded static phase, let the occurrence-local allocator spend that
	 * allowance instead of falling back to the ordinary quiet-frame budget.
	 * After a rejected floor, the same path keeps the floor disabled while the
	 * deadline ceiling prevents a retry of the failed population. */
	if (!this->d->lodInteractionSession.active() &&
	    BObolLodAdmissionPlanner::marginalStaticCapacityAllowed(
		this->d->lodStaticQualityTrial.inProgress(),
		this->d->lodStaticQualityTrial.capacityRejected(),
		this->d->lodAdmissionEvidence.capacity().retainedQualityFloorRejected()) &&
	    this->d->lodStableCalibratedRenderCostPerSecond > 0.0L) {
	    constexpr long double staticMarginalSafetyFactor = 0.80L;
	    const long double hardQualityFps =
		static_cast<long double>(this->d->staticQualityTargetFps());
	    const long double affordable = hardQualityFps > 0.0L ?
		this->d->lodStableCalibratedRenderCostPerSecond *
		    staticMarginalSafetyFactor / hardQualityFps : 0.0L;
	    size_t safeBudget =
		affordable >= static_cast<long double>(SIZE_MAX) ? SIZE_MAX :
		(affordable > 0.0L ? static_cast<size_t>(affordable) : 0);
	    const size_t deadlineCeiling =
		this->d->lodAdmissionEvidence.capacity().
		    staticDeadlineCapacityCeiling();
	    safeBudget = BObolLodAdmissionPlanner::capBudgetAtDeadlineCeiling(
		safeBudget, deadlineCeiling);
	    maximumMarginalBudget = std::max(maximumMarginalBudget, safeBudget);
	}
	size_t protectedFloorBudget = 0;
	uint64_t protectedFloorSignature = 0;
	const int64_t retainedAllocationStarted = bu_gettime();
	BObolViewLodState *retainedViewState =
	    this->d->viewAttachment->getViewLodState();
	BObolRetainedAllocationInputs inputs;
	inputs.sources = &sources;
	inputs.viewState = retainedViewState;
	/* The retained allocator owns CAD occurrence cuts, but its allowance is
	 * learned from whole-scene frame time.  Charge the immutable/non-CAD
	 * portion up front so its certificate and the deadline observer use the
	 * same currency.  Otherwise each recovery omits this constant floor and
	 * must discover it again in another interrupted presentation. */
	const size_t retainedFrameCost = retainedViewState ?
	    retainedViewState->activeRenderCost() : 0;
	const size_t retainedCadCost = retainedViewState ?
	    retainedViewState->activeCadRenderCost() : 0;
	inputs.externalPresentationCost = retainedFrameCost > retainedCadCost ?
	    retainedFrameCost - retainedCadCost : 0;
	inputs.sceneBudget = this->d->lodAdmissionEvidence.capacity().currentBudget();
	inputs.maximumMarginalBudget = maximumMarginalBudget;
	const BObolLodCapacitySearchCertificate &capacitySearch =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	const bool protectedFloorEligibleForSearch =
	    capacitySearch.phase() ==
		BObolLodCapacitySearchCertificate::Phase::INACTIVE ||
	    capacitySearch.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC;
	inputs.allowProtectedFloor = !this->d->lodInteractionSession.active() &&
	    protectedFloorEligibleForSearch &&
	    !this->d->lodStaticQualityTrial.capacityRejected() &&
	    !this->d->lodAdmissionEvidence.capacity().retainedQualityFloorRejected() &&
	    this->d->lodPresentationPolicy.
		handoffReconciliationBudget() == 0;
	inputs.maximumProtectedBudget = inputs.allowProtectedFloor ?
	    maximumProtectedBudget : 0;
	if (presentationReconciliationBudget)
	    inputs.setPresentationReconciliationBudget(
		presentationReconciliationBudget);
	inputs.viewRevision = this->d->lodViewRevision.value();
	inputs.policyRevision = this->d->lodPolicyRevision.value();
	inputs.residentAdmissionRevision = this->d->lodService ?
	    this->d->lodService->residentMeshAdmissionRevision() : 0;
	inputs.pointProxyPixelThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold;
	BObolRetainedAllocationResult allocation;
	const uint64_t allocationSliceMicroseconds = retainedViewState &&
	    retainedViewState->cadMeshPayloadCount() > 4096 &&
	    !this->d->forceTerminalLodRefinement ? 3000 : 0;
	/* A completed allocation is an immutable plan for its exact inputs.  A
	 * measured-frame barrier may request another application pass without
	 * changing any of those inputs.  Re-running minimax in that case is both
	 * O(scene) waste and numerically unsafe near a discrete budget boundary:
	 * a few equal-value cuts can alternate between two equally affordable
	 * plans forever.  Reapply the already committed plan until its occurrence
	 * cuts are unchanged; population, residency, budget, or policy mutations
	 * invalidate this certificate explicitly below. */
	const BObolRetainedAllocationResult &priorAllocation =
	    this->d->lodRetainedAllocationCertificate;
	const bool reuseCommittedAllocation = retainedViewState &&
	    priorAllocation.allocationPlanSerial != 0 &&
	    priorAllocation.allocationPlanSerial ==
		retainedViewState->activeCadAllocationPlan() &&
	    priorAllocation.inputKey() == inputs.inputKey();
	const bool reuseCoveredAllocation = reuseCommittedAllocation &&
	    retainedViewState->cadAllocationPlanCoversCurrentPopulation(
		priorAllocation.allocationPlanSerial,
		inputs.viewRevision, inputs.policyRevision,
		priorAllocation.fixedCadPresentationCost);
	BObolRetainedAllocationStatus allocationStatus;
	if (reuseCoveredAllocation) {
	    allocation = priorAllocation;
	    allocationStatus = BOBOL_RETAINED_ALLOCATION_COMPLETE;
	} else {
	    allocationStatus = bobol_retained_allocation_advance(
		this->d->lodRetainedAllocationTransaction, inputs,
		allocationSliceMicroseconds, allocation);
	}
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value())) {
	    static std::atomic<unsigned int> allocationTraceCount(0);
	    const unsigned int traceIndex = allocationTraceCount.fetch_add(1);
	    if (traceIndex < 256)
		bu_log("BObol LoD retained allocation status=%u reused=%d "
		       "budget=%zu selected=%zu demand=%zu plan=%llu "
		       "active_plan=%llu view=%llu policy=%llu\n",
		       static_cast<unsigned int>(allocationStatus),
		       reuseCoveredAllocation ? 1 : 0,
		       inputs.sceneBudget, allocation.selectedPresentationCost,
		       allocation.pixelDemandPresentationCost,
		       static_cast<unsigned long long>(
			   allocation.allocationPlanSerial),
		       static_cast<unsigned long long>(retainedViewState ?
			   retainedViewState->activeCadAllocationPlan() : 0),
		       static_cast<unsigned long long>(inputs.viewRevision),
		       static_cast<unsigned long long>(inputs.policyRevision));
	}
	if (allocationStatus == BOBOL_RETAINED_ALLOCATION_PENDING ||
	    allocationStatus == BOBOL_RETAINED_ALLOCATION_STALE) {
	    /* Preserve the prior coherent presentation while this unpublished
	     * plan advances.  The ordinary progressive pump owns the timer edge;
	     * no render is needed merely to calculate another slice. */
	    this->d->lodSubmissionPass.activate();
	    this->markProgressiveWorkPending();
	    return 0;
	}
	if (allocationStatus == BOBOL_RETAINED_ALLOCATION_FAILED) {
	    this->d->lodRetainedAllocationTransaction.reset();
	    this->d->lodRetainedAllocationCertificate =
		BObolRetainedAllocationResult();
	    append_controller_lod_diagnostic(
		this->d->lastLodDiagnostics, SbString("retained-allocation"),
		"unable to reserve the bounded stable-view allocation plan; "
		"using the conservative scalar quality floor");
	    this->d->lodRetainedAdmissionMaximumNormalizedError =
		std::numeric_limits<double>::infinity();
	    this->d->lodRetainedAdmissionMaximumProjectedErrorPixels =
		std::numeric_limits<double>::infinity();
	} else {
	    /* A discovery-time pressure sample may raise the point threshold
	     * before the complete retained population is known.  Once the allocator
	     * proves that no occurrence can use that representation, canonicalize
	     * the inert control immediately.  Leaving it active requests calibration
	     * frames which cannot change a draw record and can keep a handful of
	     * large meshes behind the deadline-recovery PoP ceiling indefinitely. */
	    if (allocation.pointProxyCandidateCount == 0 &&
		this->d->lodPresentationPointProxyPixelThreshold > 1.01f) {
		this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
		this->d->lodPointQualityPhase.reset();
		if (retainedViewState)
		    retainedViewState->
			setCadPresentationPointProxyPixelThreshold(1.0f);
		allocation.pointProxyPixelThreshold = 1.0f;
		this->d->lodRetainedAllocationTransaction.reset();
	    }
	    this->d->lodRetainedAllocationCertificate = allocation;
	    this->d->lodRetainedAdmissionMaximumNormalizedError =
		allocation.normalizedError;
	    this->d->lodRetainedAdmissionMaximumProjectedErrorPixels =
		allocation.maximumProjectedErrorPixels;
	    protectedFloorBudget = allocation.protectedFloorBudget;
	    protectedFloorSignature = allocation.protectedFloorSignature;
	    if (!reuseCoveredAllocation &&
		getenv("BOBOL_LOD_TRACE_ALLOCATOR") &&
		this->d->lodRetainedAllocationTransaction)
		bobol_retained_allocation_trace(
		    this->d->lodRetainedAllocationTransaction);
	}
	const int64_t retainedAllocationCompleted = bu_gettime();
	const int64_t retainedAllocationElapsed =
	    retainedAllocationCompleted >= retainedAllocationStarted ?
		retainedAllocationCompleted - retainedAllocationStarted : 0;
	this->d->lodRetainedAdmissionQualityViewRevision =
	    this->d->lodViewRevision.value();
	this->d->lodRetainedAdmissionQualityPolicyRevision =
	    this->d->lodPolicyRevision.value();
	this->d->lodRetainedAdmissionQualityPointProxyPixelThreshold =
	    this->d->lodPresentationPointProxyPixelThreshold;
	this->d->setRetainedQualityFloor(
	    protectedFloorBudget, protectedFloorSignature, sceneActiveCost,
	    sceneMinimumActiveCost);
	if (getenv("BOBOL_LOD_TRACE_BUDGET") ||
	    getenv("BOBOL_LOD_TRACE_ALLOCATOR"))
	    bu_log("BObol LoD retained importance ceiling=%.6f "
		   "budget=%zu selected=%zu certified=%zu "
		   "marginal_limit=%zu protected_limit=%zu reused=%d "
		   "elapsed_us=%lld\n",
		   this->d->lodRetainedAdmissionMaximumNormalizedError,
		   this->d->lodAdmissionEvidence.capacity().currentBudget(),
		   this->d->lodRetainedAllocationCertificate.
		       selectedPresentationCost,
		   this->d->lodRetainedAllocationCertificate.
		       certifiedPresentationBudget,
		   maximumMarginalBudget,
		   inputs.effectiveMaximumProtectedBudget(),
		   reuseCoveredAllocation ? 1 : 0,
		   static_cast<long long>(retainedAllocationElapsed));
    }
    SbBool boundedScenePass = FALSE;
    if (!this->d->lodInteractionSession.active() &&
	this->d->lodPoseContinuity.visibilityCensusDeferred()) {
	/* The motion fast path may have advanced the source cursor to its end
	 * without projecting any occurrences.  Start the first quiet pass from
	 * the beginning and discard any partial coverage counters so its terminal
	 * denominator belongs wholly to this camera. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodCoveragePolicy.clearPassCounters();
	this->d->lodPoseContinuity.completeVisibilityCensus();
    }
    for (size_t i = this->d->lodSubmissionSourceIndex;
	 i < sources.size();) {
	const size_t capacity = service->availableResultTaskCapacity();
	if (!capacity) {
	    break;
	}
	SoBRLDatabaseSource *source = sources[i];
	if (!source) {
	    this->d->lodSubmissionSourceIndex = ++i;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    continue;
	}
	if (this->d->lodSubmissionDelta.active() &&
	    !this->d->lodSubmissionDelta.targets(source)) {
	    this->d->lodSubmissionSourceIndex = ++i;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    continue;
	}
	/* A compact source records its source-backed PoP requests explicitly,
	 * making absence authoritative and cheap to test.  Non-compact sources
	 * may still contain direct SoBRLMeshShape children whose request
	 * metadata is owned by the shape, so let the submit action inspect those
	 * sources.  This preserves the custom scene-node boundary
	 * without rescanning terminal compact analytic meshes on every view. */
	if (source->hasCompactInstanceIndex() &&
	    !source->hasDisplayMeshLodRequests()) {
	    this->d->lodSubmissionSourceIndex = ++i;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    continue;
	}

	struct db_i *dbip = source->getDatabase();
	if (!dbip) {
	    append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
					     source->path.getValue(),
					     "database source has no database for LoD submission");
	    this->d->lodSubmissionSourceIndex = ++i;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    continue;
	}

	BObolViewLodState *viewLodState =
	    this->d->viewAttachment->getViewLodState();
	const size_t sourceMeshRequests =
	    source->getDisplayMeshLodRequestCount();
	/*
	 * Rotation/translation deliberately keeps existing per-occurrence
	 * prefixes.  When every current source-backed mesh occurrence already
	 * has a mesh payload, the renderer's bounded prepared-frame refresh
	 * handles exact frustum visibility and there is no LoD work to discover.
	 * The quiet policy revision performs the new-view pixel-exact pass;
	 * inventory deltas and incomplete coverage always use the normal path.
	 */
	const SbBool scaleInteraction =
	    this->d->lodViewDemandPolicy.scaleChangingInteraction(
		this->d->lodInteractionSession.active() != FALSE) ? TRUE : FALSE;
	const SbBool poseOnlyFullyResident =
	    this->d->lodInteractionSession.active() &&
	    !scaleInteraction &&
	    !this->d->lodSubmissionDelta.active() &&
	    sourceMeshRequests > 0 && viewLodState &&
	    viewLodState->cadMeshPayloadCountForSource(source) >=
		sourceMeshRequests;
	if (poseOnlyFullyResident) {
	    /* Reusing the retained presentation during motion is intentionally
	     * O(1), but resident payload count is not a visibility census.  One
	     * shared mesh payload may represent several visible occurrences, and
	     * point aggregation may leave other visible occurrences without a mesh
	     * payload at all.  Substituting that count for projected visibility let
	     * a pose-only frame publish a false terminal proof (709 resident Generic
	     * Twin payloads versus 1235 visible occurrences) and made exact-view
	     * history quality depend on timing.
	     *
	     * Keep the responsive retained frame, leave coverage active, and let the
	     * first quiet policy pass perform the authoritative camera census before
	     * convergence or quality history can become terminal. */
	    this->d->lodPoseContinuity.deferVisibilityCensus();
	    this->d->lodSubmissionSourceIndex = ++i;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    continue;
	}

	SoBRLMeshLodSubmitAction action;
	action.setService(service);
	action.setDatabase(dbip, controller_database_id(dbip),
			   source->sourceRevision.getValue());
	action.setViewInfo(&view);
	/* Motion-time pixel error is a presentation response, not a residency
	 * target.  During zoom, keep requesting at least the ordinary one-pixel
	 * physical demand so workers can read the needed suffix while input is
	 * still arriving.  The aggregate cost allocator and renderer ceiling below
	 * remain authoritative over which resident cut is exposed in a frame.
	 * Without this split, a retained mesh was merely magnified at the 3-4 pixel
	 * interaction target and then walked through several blocky prefixes only
	 * after button-up. */
	const float residentPixelError = scaleInteraction ?
	    std::min(scenePixelError, 1.0f) : scenePixelError;
	action.setViewVolume(&lodViewVolume, residentPixelError);
	action.setPointProxyPixelThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold);
	action.setAllowTerminalMeshAdmission(
	    meshFirstSceneSafe ? TRUE : FALSE);
	/* Source/inventory discovery already supplied a leaf proxy.  Complete the
	 * useful-coverage census without launching one PoP build per visible leaf;
	 * a fresh budgeted quality pass follows immediately.  A zoom census keeps
	 * existing cuts and is therefore not structural-only.  A sole normal PoP
	 * source can publish its globally classified preview while persistence
	 * continues; multi-leaf scenes still establish their complete structural
	 * frontier first. */
	const SbBool globalColdPreview =
	    (meshFirstSceneSafe ||
	     controller_lod_global_preview_requested(sourceMeshRequests)) ?
		TRUE : FALSE;
	action.setStructuralCoverageOnly(
	    this->d->lodCoveragePolicy.active() &&
	    !this->d->lodViewDemandPolicy.scaleDemandRefreshActive() &&
	    !globalColdPreview);
	action.setStructuralPresentationRepair(
	    this->d->lodStructuralRepair.active() ? TRUE : FALSE);
	if (this->d->lodStructuralRepair.active())
	    action.setStructuralCoverageCostReservation(
		this->d->lodStructuralRepair.coverageCostReservation());
	action.setSelectedOccurrenceCount(selectedOccurrenceCount);
	action.setGeneration(generation);
	action.setRevisions(this->d->lodViewRevision.value(),
	    this->d->lodPolicyRevision.value());
	action.setRefreshMissing(refreshMissing);
    action.setReset(resetExisting);
	action.setViewLodState(viewLodState);
	action.setProjectedDemandCache(&this->d->lodProjectedDemandCache);
	action.setRefinementCutCeiling(
	    this->d->lodInteractionSession.active() ?
		this->d->lodInteractiveProgressiveCeiling : -1);
	/* Coarsening a retained prefix only changes the draw count/snap level;
	 * it does not rebuild or reread the mesh.  Permit it during motion so a
	 * previously settled, expensive cut cannot pin interactive FPS. */
	const SbBool scaleDemandChanged =
	    (!this->d->lodInteractionSession.active() &&
	     !this->d->lodPoseContinuity.retainOccurrenceCuts()) ||
	    scaleInteraction;
	/*
	 * A pose-only interaction already has the assembly-wide render ceiling
	 * as its instantaneous FPS control.  Rewriting hundreds of retained
	 * occurrence cuts underneath that same ceiling cannot reduce submitted
	 * geometry, but it invalidates the prepared command record and turns
	 * rotation into repeated 40-50 ms rebuilds.  Preserve the resident
	 * per-occurrence cuts until zoom, stable view allocation, or later
	 * memory compaction genuinely requires them to change.
	 */
	const SbBool applyRetainedAdmission = retainedAllocationPass;
	/* Keep occurrence cuts monotonic throughout active camera input.  The
	 * renderer-wide ceiling is the O(1), reversible FPS control; rewriting
	 * every retained cut on each wheel event causes visible regression and
	 * invalidates prepared commands.  A quiet scale pass may move directly to
	 * its new pixel target, but a presentation handoff by itself is not
	 * retained-admission authority.  Only explicit recovery or measured
	 * overload may force a broad downgrade. */
	action.setAllowCutDowngrade(
	    BObolLodAvailabilityScheduler::
		presentationCutDowngradeAllowed(
		    this->d->lodInteractionSession.active() != FALSE,
		    this->d->lodInteractionSession.gestureActive() != FALSE,
		    this->d->lodAvailabilityLedger.residencyDrainActive(),
		    scaleDemandChanged != FALSE,
		    applyRetainedAdmission != FALSE));
	/* Pose-only orthographic interaction preserves every existing prefix while
	 * the camera is moving.  Once the view is quiet, however, a newly visible
	 * occurrence still has to reach its current projected demand, and a
	 * measured overload recovery must be able to spend the stable budget after
	 * its coherent cheaper cut is presented.  Treating "preserve existing" as
	 * "forbid every upward retarget" left Hubble at its minimum PoP prefixes
	 * after a rotation even though the resident data and 20 Hz budget could
	 * afford the pre-motion quality.  Downgrade permission above remains the
	 * independent guard which prevents pose changes from needlessly making an
	 * already useful cut coarser. */
	action.setAllowRetainedRefinement(
	    (this->d->forceTerminalLodRefinement ||
	     !this->d->lodInteractionSession.active() || scaleDemandChanged) &&
	    !this->d->lodAvailabilityLedger.residencyDrainActive() &&
	    !this->d->lodPresentationTransaction.barrierPending() ? TRUE : FALSE);
	/* A Schmitt boundary is useful while successive input events perturb the
	 * projection by fractions of a pixel.  It must not redefine the quiet
	 * view's convergence target: the first post-gesture pass recomputes the
	 * exact producer cut and either presents it or records genuine budget /
	 * residency debt. */
	action.setCutHysteresisEnabled(
	    this->d->lodInteractionSession.active());
	/* Once coverage is proven, residency follows pixel demand independently
	 * of the calibrated draw allowance.  This applies during zoom and during
	 * quiet warm-cache convergence: one worker read may append the complete
	 * demanded immutable suffix while the currently affordable cut remains on
	 * screen.  Working-set and resident-byte admission remain authoritative. */
    const SbBool quietCoveredPrefetch =
	!this->d->lodInteractionSession.active() &&
	this->d->lodCoveragePolicy.coverageComplete() ? TRUE : FALSE;
    action.setAllowResidentPrefetch(
	    (scaleInteraction || quietCoveredPrefetch) &&
	    !this->d->lodPresentationTransaction.barrierPending() ? TRUE : FALSE);
    action.setAllowResidentPrefetchPastAllocation(scaleInteraction);
	action.setAllowRepresentationRefinement(
	    scaleDemandChanged &&
	    !this->d->lodAvailabilityLedger.residencyDrainActive() &&
	    !this->d->lodPresentationTransaction.barrierPending() ? TRUE : FALSE);
	/* A calibrated face budget governs richness above the minimum drawable
	 * mesh floor.  Returning a visible resident PoP payload to its box is
	 * both a visual regression and a false economy: the renderer can batch
	 * tiny minimum prefixes into the aggregate point channel without
	 * discarding the useful retained asset. */
	action.setPreserveMeshCoverage(TRUE);
	size_t refinementBudget = this->d->lodAdmissionCursor.refinementRemaining();
	/* Structural proxies occupy the same view while their first useful mesh
	 * replacements are planned, but they are not retained mesh population.
	 * Charging their estimated draw cost against the only initial-mesh
	 * allowance can reduce that allowance to zero and leave a cold view at
	 * boxes forever.  Admit one bounded replacement wave from the current
	 * scene allowance; normal result publication and the next measured frame
	 * immediately return to the ordinary refinement remainder. */
	const SbBool structuralReplacementBootstrap =
	    refinementBudget == 0 &&
	    viewLodState && this->getActiveLodCadPayloadCount() == 0 &&
	    !this->d->lodStructuralRepair.active();
	if (structuralReplacementBootstrap)
	    refinementBudget = std::min(
		this->d->lodAdmissionEvidence.capacity().currentBudget(),
		BObolLodCapacityEvidence::singleOccurrenceBootstrapBudget());
	action.setRefinementCostBudget(refinementBudget);
	if (this->d->lodAdmissionCursor.singleOccurrenceBootstrap())
	    action.setInitialProviderCostBudget(
		this->d->lodAdmissionCursor.refinementRemaining());
	if (this->d->lodDiscretePopulationTrialPermit.available()) {
	    size_t discreteAllowance =
		BObolLodQualityPolicy::discreteTrialOverBudgetAllowance(
		    sceneActiveCost,
		    this->d->lodAdmissionEvidence.capacity().currentBudget());
	    action.setOneOverBudgetRefinementLimit(discreteAllowance);
	}
	/* Every bounded window owns a disjoint part of the pinned occurrence
	 * plan and consumes from this carried scene-wide remainder.  The action
	 * deliberately skips its full priority recovery when a finite window is
	 * configured, avoiding an O(scene size) input stall. */
	if (applyRetainedAdmission)
	    action.setRetainedSceneUpgradeCostBudget(
		this->d->lodAdmissionCursor.retainedAdmissionRemaining());
	if (applyRetainedAdmission)
	    action.setRetainedSceneMaximumNormalizedError(
		this->d->lodRetainedAdmissionMaximumNormalizedError);
	action.setSubmissionTaskLimit(capacity);
	/* Planning and request construction execute on the host thread.  Keep
	 * their per-pump slice below a frame even when a nominal 2k window has
	 * become expensive due to large retained maps.  Deterministic/offline
	 * callers may explicitly remove this interactive deadline. */
	/* The Qt host deliberately spaces ordinary progressive pumps by one
	 * 16 ms event-loop interval.  Applying the same 3 ms planning slice to a
	 * quiet many-leaf view leaves the worker pool empty for most of that
	 * interval: a 50k warm-cache scene then spends more wall time waiting for
	 * the next owner-thread planning slice than reading or publishing meshes.
	 *
	 * Preserve the short motion slice, where input latency is the primary
	 * contract.  A quiet view may spend up to half of one 60 Hz frame filling
	 * the already byte- and result-bounded service queue.  This does not admit
	 * additional memory or rendering work; the service working-set governor,
	 * result reservations, scene budget, and presentation deadline remain the
	 * independent authorities.  Offline callers continue to request an
	 * unbounded deterministic slice explicitly. */
	const uint64_t submissionTimeLimitMicroseconds =
	    this->d->forceTerminalLodRefinement ? 0 :
	    (this->d->lodInteractionSession.active() ? 3000 : 8000);
	action.setSubmissionTimeLimit(submissionTimeLimitMicroseconds);
	/* Never make a complete many-leaf projection/sort one GUI operation.
	 * Initial mesh coverage has no user-visible ordering contract; consume
	 * compact-index order in bounded waves and project each occurrence
	 * against the newest camera when its wave is reached.  Already-retained
	 * payloads are independently re-admitted by screen priority above, so
	 * coarsening still protects the visibly important population.
	 *
	 * This applies to quiet views too.  The former "stable full sort" reached
	 * 66 ms at only 14k leaves and would exceed 200 ms at 50k, making console
	 * input and renewed motion wait behind an LoD planning frame. */
	static const size_t quietCompactWave = 2048;
	static const size_t interactiveCompactWave = 512;
	const size_t compactWave = this->d->lodInteractionSession.active() ?
	    interactiveCompactWave : quietCompactWave;
	const int compactCount = source->getCompactInstanceCount();
	const SbBool boundedLargeCompact =
	    compactCount > static_cast<int>(compactWave) ?
	    TRUE : FALSE;
	if (boundedLargeCompact)
	    boundedScenePass = TRUE;
	const bool usingMemoryAdmissionFrontier =
	    boundedLargeCompact &&
	    this->d->lodPlanningObligations.
		residentAdmissionRetryPending() &&
	    !this->d->lodSubmissionDelta.active() && viewLodState;
	const bool currentViewSourceCensusComplete =
	    this->d->hasCompleteLodConvergenceCandidateCensus(source);
	/* Point-proxy aggregation deliberately leaves subpixel occurrences without
	 * mesh payloads.  Comparing the mesh-payload count with the source request
	 * count can therefore never prove readiness in precisely the large scenes
	 * which need the sparse path: every refinement pass falls back to another
	 * complete source projection.
	 *
	 * The dense per-entry visibility census is the correct authority.  It is
	 * cleared by camera/population revisions, rebuilt by the required full
	 * demand pass, and preserved only across exact source deltas.  Once that
	 * proof and useful-coverage proof are both current, payloads absent from the
	 * unsatisfied set are either satisfied meshes or valid point
	 * representations; structural boxes have their own exact-frame repair
	 * frontier.  A stricter pixel target first activates a coverage/demand pass,
	 * so point occurrences which become mesh-worthy cannot be skipped here. */
	const bool usingUnsatisfiedFrontier =
	    boundedLargeCompact &&
	    !this->d->lodSubmissionDelta.active() &&
	    !this->d->lodStructuralRepair.active() &&
	    !this->d->lodCoveragePolicy.active() &&
	    this->d->lodCoveragePolicy.coverageComplete() &&
	    !this->d->lodAdmissionCursor.retainedAdmission() &&
	    viewLodState && sourceMeshRequests > 0 &&
	    (usingMemoryAdmissionFrontier ||
	     currentViewSourceCensusComplete);
	bool selectiveDeltaPlan = false;
	if (this->d->lodSubmissionDelta.active() &&
	    !this->d->lodSubmissionSourcePlan.validFor(source)) {
	    for (const auto &deltaPlan : this->d->lodSubmissionDelta.plans()) {
		if (deltaPlan.first != source)
		    continue;
		selectiveDeltaPlan = true;
		this->d->lodSubmissionSourcePlan.assign(
		    source, deltaPlan.second);
		break;
	    }
	} else if (this->d->lodSubmissionDelta.active()) {
	    selectiveDeltaPlan =
		this->d->lodSubmissionDelta.selectiveEntries(source) != NULL;
	}
	if (boundedLargeCompact && applyRetainedAdmission &&
	    !this->d->lodSubmissionSourcePlan.validFor(source)) {
	    this->d->lodSubmissionSourcePlan.begin(source);
	    controller_prioritize_lod_recovery(source, viewLodState,
		this->d->lodSubmissionSourcePlan.entries(),
		selectedOccurrenceCount == 1);
	}
	if (usingUnsatisfiedFrontier &&
	    !this->d->lodSubmissionSourcePlan.validFor(source)) {
	    std::vector<SbString> unsatisfiedKeys;
	    const uint64_t admissionRevision = this->d->lodService ?
		this->d->lodService->residentMeshAdmissionRevision() : 0;
	    if (usingMemoryAdmissionFrontier) {
		viewLodState->retriableMemoryLimitedCadOccurrenceKeys(
		    source, admissionRevision, unsatisfiedKeys);
	    } else {
		viewLodState->unsatisfiedCadOccurrenceKeys(
		    source, admissionRevision, unsatisfiedKeys);
	    }
	    this->d->lodSubmissionSourcePlan.begin(source);
	    controller_prioritize_lod_frontier(source, viewLodState,
		unsatisfiedKeys, this->d->lodSubmissionSourcePlan.entries(),
		selectedOccurrenceCount == 1);
	    /* Unsatisfied residency is deliberately retained while an occurrence
	     * leaves the frustum so a later camera pose can reuse its immutable
	     * mesh.  It is not actionable work in the current view, however.  The
	     * exact dense census is already the authority which enabled this sparse
	     * frontier; intersect with it here.  Otherwise an off-screen occurrence
	     * is projected, skipped, and reconstructed as the same one-entry plan
	     * on every GUI pump, preventing a quiet view from ever converging. */
	    this->d->lodSubmissionSourcePlan.entries().erase(
		std::remove_if(
		    this->d->lodSubmissionSourcePlan.entries().begin(),
		    this->d->lodSubmissionSourcePlan.entries().end(),
		    [this, source](size_t entryIndex) {
			return !this->d->isVisibleLodConvergenceCandidate(
			    source, entryIndex);
		    }),
		this->d->lodSubmissionSourcePlan.entries().end());
	}
	/* Transition pacing is an input-time latency tool.  During active zoom it
	 * prevents one exceptional leaf (Lucy is the canonical case) from making
	 * the whole gesture wait for a very large suffix.  Once the view is quiet,
	 * the retained allocator has already established the scene-global,
	 * significance-ordered budget.  Walking every occurrence through one PoP
	 * transition per complete scene pass at that point adds no fairness: it
	 * serializes immutable residency, rebuilds the same many-leaf command
	 * record between tiny result waves, and makes a pose-only restore look like
	 * cold convergence.  Jump directly to each allocated target after input
	 * stops; working-set admission and the hard presentation deadline remain
	 * the independent memory and responsiveness bounds. */
	if (scaleInteraction)
	    action.setTransitionLimitedRefinement(TRUE);
	if (this->d->lodCoveragePolicy.active()) {
	    /*
	     * A bounded index-order window cannot know whether a later window
	     * still contains an uncovered visible leaf.  For cold/source coverage,
	     * count the already-published structural proxy and defer mesh creation
	     * until a complete projected pass proves useful coverage.  A scale
	     * refresh keeps its preceding mesh cuts and refines them in place while
	     * the same scan verifies visibility.
	     */
	    if (!this->d->lodViewDemandPolicy.scaleDemandRefreshActive())
		action.setAllowRetainedRefinement(FALSE);
	    if (boundedLargeCompact)
		this->d->lodCoveragePolicy.markBoundedSource();
	}
	if (boundedLargeCompact &&
	    !this->d->lodSubmissionSourcePlan.validFor(source)) {
	    this->d->lodSubmissionSourcePlan.begin(source);
	    this->d->lodSubmissionSourcePlan.entries().resize(
		static_cast<size_t>(compactCount));
	    for (size_t planIndex = 0;
		 planIndex < this->d->lodSubmissionSourcePlan.size();
		 planIndex++)
		this->d->lodSubmissionSourcePlan.entries()[planIndex] =
		    planIndex;
	} else if (boundedLargeCompact && !usingUnsatisfiedFrontier &&
	    (!this->d->lodSubmissionDelta.active() || !selectiveDeltaPlan) &&
	    this->d->lodSubmissionSourcePlan.validFor(source) &&
	    this->d->lodSubmissionSourcePlan.size() <
		static_cast<size_t>(compactCount)) {
	    /* Structural streaming appends leaves while a pinned scan is in
	     * flight.  Extend the plan tail without restarting the consumed
	     * prefix; restarting at zero on every batch starves late leaves. */
	    const size_t oldCount =
		this->d->lodSubmissionSourcePlan.size();
	    this->d->lodSubmissionSourcePlan.entries().resize(
		static_cast<size_t>(compactCount));
	    for (size_t planIndex = oldCount;
		 planIndex < this->d->lodSubmissionSourcePlan.size();
		 planIndex++)
		this->d->lodSubmissionSourcePlan.entries()[planIndex] =
		    planIndex;
	}
	if (this->d->lodSubmissionSourcePlan.validFor(source))
	    action.setCompactEntryPlanView(
		&this->d->lodSubmissionSourcePlan.entries());
	action.setCompactEntryRange(this->d->lodSubmissionEntryOffset,
	    boundedLargeCompact ? compactWave : SIZE_MAX);
	if (this->d->lodForcedCut)
	    action.setForcedCut(*this->d->lodForcedCut);
	const bool fullViewCandidatePass =
	    !usingUnsatisfiedFrontier &&
	    !this->d->lodSubmissionDelta.active();
	if (fullViewCandidatePass && source->hasCompactInstanceIndex() &&
	    this->d->lodSubmissionEntryOffset == 0)
	    this->d->beginLodConvergenceCandidateCensus(source,
		static_cast<size_t>(std::max(0, compactCount)));
	const int64_t sourceActionStarted = bu_gettime();
	action.apply(source);
	const int64_t sourceActionCompleted = bu_gettime();
	const int64_t sourceActionElapsed =
	    sourceActionCompleted >= sourceActionStarted ?
		sourceActionCompleted - sourceActionStarted : 0;
	if (getenv("BOBOL_LOD_TRACE_SLOW_SUBMISSION") &&
	    sourceActionElapsed >= 10000)
	    bu_log("BObol slow LoD source action path=%s elapsed_us=%lld "
		"entries=%d range_first=%zu range_count=%zu plan=%zu "
		"visited=%u tasks=%u cuts=%u deferred=%d next=%zu "
		"time_limit_us=%llu force_terminal=%d\n",
		source->path.getValue().getString(),
		static_cast<long long>(sourceActionElapsed), compactCount,
		this->d->lodSubmissionEntryOffset,
		boundedLargeCompact ? compactWave : SIZE_MAX,
		this->d->lodSubmissionSourcePlan.size(),
		action.getVisitedMeshCount(), action.getSubmittedTaskCount(),
		action.getUpdatedCutCount(),
		action.hasDeferredCompactEntries() ? 1 : 0,
		action.getCompactEntryNext(),
		static_cast<unsigned long long>(
		    submissionTimeLimitMicroseconds),
		this->d->forceTerminalLodRefinement ? 1 : 0);
	if (action.getOneOverBudgetRefinementUsed())
	    this->d->lodDiscretePopulationTrialPermit.consume();
	if (source->hasCompactInstanceIndex())
	    this->d->observeLodConvergenceCandidateVisibility(source,
		static_cast<size_t>(std::max(0, compactCount)),
		action.getCompactEntryVisibilityObservations());
	if (fullViewCandidatePass) {
	    this->d->lodSubmissionVisibleCount.observe(
		source, action.getVisibleMeshCount());
	}
	if (this->d->lodCoveragePolicy.active()) {
	    const size_t visible = action.getVisibleMeshCount();
	    const size_t covered = action.getCoveredVisibleMeshCount();
	    this->d->lodCoveragePolicy.observe(visible, covered);
	}
	const bool sparseFrontierTrace =
	    getenv("BOBOL_LOD_TRACE_SPARSE_FRONTIER") &&
	    this->d->lodSubmissionSourcePlan.validFor(source) &&
	    this->d->lodSubmissionSourcePlan.size() <= 512;
	if (sparseFrontierTrace || controller_lod_trace_enabled(
		"BOBOL_LOD_TRACE_SUBMISSION",
		this->d->lodViewRevision.value())) {
	    static std::atomic<unsigned int> submissionTraceCount(0);
	    const unsigned int traceIndex =
		submissionTraceCount.fetch_add(1);
	    if (traceIndex < 512)
		bu_log("BObol LoD submission at=%lld path=%s entries=%d "
		       "offset=%zu visited=%u visible=%zu covered=%zu "
		       "tasks=%u cuts=%u skipped=%u deferred=%d next=%zu "
		       "coverage_pass=%d resident_pending=%u retained_pending=%u "
		       "refinement_blocked=%u missing_blocked=%u "
		       "quality_limited=%u admission_blocked=%u plan_size=%zu "
		       "diagnostics=%s\n",
		       static_cast<long long>(bu_gettime()),
		       source->path.getValue().getString(), compactCount,
		       this->d->lodSubmissionEntryOffset,
		       action.getVisitedMeshCount(),
		       action.getVisibleMeshCount(),
		       action.getCoveredVisibleMeshCount(),
		       action.getSubmittedTaskCount(),
		       action.getUpdatedCutCount(),
		       action.getSkippedMeshCount(),
		       action.hasDeferredCompactEntries() ? 1 : 0,
		       action.getCompactEntryNext(),
		       this->d->lodCoveragePolicy.active() ? 1 : 0,
		       action.getPendingResidentRefinementCount(),
		       action.getPendingRetainedRefinementCount(),
		       action.getRefinementBudgetBlockedCount(),
		       action.getMissingMeshBudgetBlockedCount(),
		       action.getRetainedQualityLimitedCount(),
		       action.getRetainedAdmissionBlockedCount(),
		       this->d->lodSubmissionSourcePlan.size(),
		       action.getDiagnostics().getString());
	}
	if (fullViewCandidatePass &&
	    !action.hasDeferredCompactEntries()) {
	    if (source->hasCompactInstanceIndex())
		this->d->completeLodConvergenceCandidateCensus(source);
	    else
		this->d->setLodConvergenceCandidateCount(source,
		    this->d->lodSubmissionVisibleCount.valueFor(source));
	}
	if (!this->d->lodSubmissionSourcePlan.valid() &&
	    source->hasCompactInstanceIndex()) {
	    this->d->lodSubmissionSourcePlan.begin(source);
	    action.getCompactEntryPlan(
		this->d->lodSubmissionSourcePlan.entries());
	}

	this->d->lastLodVisitedMeshCount += action.getVisitedMeshCount();
	this->d->lastLodSubmittedTaskCount += action.getSubmittedTaskCount();
	this->d->lastLodUpdatedCutCount += action.getUpdatedCutCount();
	if (action.getSubmittedTaskCount() > 0 ||
	    action.getUpdatedCutCount() > 0)
	    this->d->lodRetainedPass.noteAdmittedWork();
	/* Accumulate every cut mutation across the complete pass.  Large compact
	 * sources need this across several bounded windows, but the presentation
	 * contract is identical for a small source completed by one action: a
	 * retained allocation which changed a cut must publish that population
	 * before its certificate can release a renderer-wide handoff ceiling. */
	if (action.getUpdatedCutCount() > 0)
	    this->d->lodRetainedPass.noteCutAdvanced();
	if (action.getPendingRetainedRefinementCount() > 0)
	    this->d->lodRetainedPass.noteRefinementPending();
	if (action.getPendingResidentRefinementCount() > 0)
	    this->d->lodRetainedPass.noteResidencyPending();
	if (action.getRefinementBudgetBlockedCount() > 0) {
	    /* A minimax quality ceiling is a terminal allocation for the current
	     * allowance, but it is still the witness which says that allowance
	     * prevented pixel-target convergence.  Route both an unreachable
	     * allocated cut and a deliberately coarser minimax cut through the
	     * same bounded unchanged-frame calibration.  BObolLodCapacityEvidence's
	     * three-sample series and the one-shot headroom witness make a truly
	     * saturated population terminal; suppressing the pure quality case
	     * here instead stranded conservative cold seeds forever (most visibly
	     * a one-leaf Lucy view at its first few thousand faces). */
	    this->d->lodRetainedPass.noteBudgetBlocked();
	}
	if (getenv("BOBOL_LOD_TRACE_BUDGET") &&
	    (action.getRefinementCostBudgetUsed() > 0 ||
	     action.getRefinementBudgetBlockedCount() > 0))
	    bu_log("BObol LoD frame budget source=%s used_faces=%zu "
		   "blocked=%u cuts=%u tasks=%u remaining_before=%zu\n",
		   source->path.getValue().getString(),
		   action.getRefinementCostBudgetUsed(),
		   action.getRefinementBudgetBlockedCount(),
		   action.getUpdatedCutCount(), action.getSubmittedTaskCount(),
		   this->d->lodAdmissionCursor.refinementRemaining());
	this->d->lodAdmissionCursor.consumeRefinement(
	    action.getRefinementCostBudgetUsed());
	if (applyRetainedAdmission) {
	    this->d->lodAdmissionCursor.consumeRetainedAdmission(
		action.getRetainedSceneUpgradeCostBudgetUsed());
	}
	this->d->lastLodSkippedMeshCount += action.getSkippedMeshCount();
	const size_t missingMeshBlocked =
	    action.getMissingMeshBudgetBlockedCount();
	if (missingMeshBlocked > 0) {
	    SbString diagnostic;
	    diagnostic.sprintf("%zu visible mesh requests await the current "
		"scene budget", missingMeshBlocked);
	    append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
		 source->path.getValue(), diagnostic.getString());
	}
	this->d->lodRetainedPass.addMissingMeshBudgetBlocked(
	    missingMeshBlocked);
	if (action.getDiagnostics().getLength() > 0)
	    append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
					     source->path.getValue(),
					     action.getDiagnostics().getString());
	if (action.hasDeferredCompactEntries()) {
	    this->d->lodSubmissionEntryOffset = action.getCompactEntryNext();
	    this->d->lodSubmissionSourceIndex = i;
	    break;
	}
	this->d->lodSubmissionSourceIndex = ++i;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
    }

    const bool completedPass =
	this->d->lodSubmissionSourceIndex >= sources.size();
    if (completedPass)
	this->d->lodPlanningObligations.retireResidentAdmissionRetry();
    const bool completedStructuralRepair =
	completedPass && this->d->lodStructuralRepair.active();
    const size_t completedMissingMeshBudgetBlockedCount = completedPass ?
	this->d->lodRetainedPass.missingMeshBudgetBlockedCount() : 0;
    if (completedStructuralRepair) {
	this->d->lodStructuralRepair.reset();
    }
    if (completedPass)
	this->d->lodRetainedPass.clearMissingMeshBudgetBlocked();
    const bool completedPassRetainedAllocation =
	completedPass &&
	this->d->lodSubmissionIntent.retainedAdmission();
    const bool completedPassChangedCut = completedPass &&
	this->d->lodRetainedPass.cutAdvanced();
    const bool completedPassBudgetBlocked = completedPass &&
	this->d->lodRetainedPass.budgetBlocked();
    const bool completedPassRefinementPending = completedPass &&
	this->d->lodRetainedPass.refinementPending();
    const bool residentWorkPending = completedPass && service &&
	(service->activeTaskCountForGeneration(generation) > 0 ||
	 service->queuedResultCountForGeneration(generation) > 0 ||
	 this->d->lodAvailabilityLedger.resultsPending() > 0);
    const BObolLodAvailabilityScheduler::CompletedPassSuccessor
	completedPassSuccessor =
	    BObolLodAvailabilityScheduler::completedPassSuccessor(
		completedPass,
		this->d->lodAvailabilityLedger.residencyDrainActive(),
		this->d->lodAvailabilityLedger.residentGrowthPending(),
		this->d->lodRetainedPass.residencyPending(),
		residentWorkPending,
		this->d->lodPointQualityPhase.triangleRecoveryPending(),
		this->d->lodPointQualityPhase.presentationPending(),
		completedPassBudgetBlocked);
    if (completedPass && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_PASS", this->d->lodViewRevision.value())) {
	static std::atomic<unsigned int> completedPassTraceCount(0);
	const unsigned int traceIndex =
	    completedPassTraceCount.fetch_add(1);
	if (traceIndex < 512)
		    bu_log("BObol LoD completed pass coverage=%d "
			   "bounded_seen=%d coverage_complete=%d delta=%d "
		   "source_rescan=%d "
		   "handoff=%d refinement_pending=%d residency_pending=%d "
		   "structural_repair=%d cut_advanced=%d "
		   "budget_rescan=%d active_faces=%zu active_cost=%zu "
		   "cost_budget=%zu "
		   "visited=%u cuts=%u\n",
		   this->d->lodCoveragePolicy.active() ? 1 : 0,
		   this->d->lodCoveragePolicy.sawBoundedSource() ? 1 : 0,
		   this->d->lodCoveragePolicy.coverageComplete() ? 1 : 0,
		   this->d->lodSubmissionDelta.active() ? 1 : 0,
		   this->d->lodSubmissionPass.rescanPending() ? 1 : 0,
		   this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
		   this->d->lodRetainedPass.refinementPending() ? 1 : 0,
		   this->d->lodRetainedPass.residencyPending() ? 1 : 0,
		   completedStructuralRepair ? 1 : 0,
		   this->d->lodRetainedPass.cutAdvanced() ? 1 : 0,
		   this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() ? 1 : 0,
		   sceneActiveFaces, sceneActiveCost,
		   this->d->lodAdmissionEvidence.capacity().currentBudget(),
		   this->d->lastLodVisitedMeshCount,
		   this->d->lastLodUpdatedCutCount);
    }
    const BObolLodCoveragePolicy::Completion coverageCompletion =
	this->d->lodCoveragePolicy.completeIfReady(
	    completedPass &&
		!this->d->lodPoseContinuity.visibilityCensusDeferred(),
	    this->d->lodSubmissionPass.rescanPending());
    const bool retainedImportanceCensusCompleted =
	this->d->lodPlanningObligations.importanceCensusPending() &&
	coverageCompletion.completed && !coverageCompletion.missing;
    if (retainedImportanceCensusCompleted) {
	/* The pass which just finished recorded exact projected demand without
	 * discarding the useful retained population.  Reallocate that population
	 * once inside the already measured scene allowance.  This is neither an
	 * overload witness nor a cache/residency operation. */
	this->d->lodPlanningObligations.retireImportanceCensus();
	this->d->requestRetainedReallocation();
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD completed importance census view=%llu "
		   "policy=%llu active_cost=%zu budget=%zu\n",
		   (unsigned long long)this->d->lodViewRevision.value(),
		   (unsigned long long)this->d->lodPolicyRevision.value(),
		   sceneActiveCost,
		   this->d->lodAdmissionEvidence.capacity().currentBudget());
    }
    if (coverageCompletion.completed) {
	/* A targeted delta may address a small compact source.  Its single
	 * synchronous action is already a complete coverage pass; bounded and
	 * unbounded sources publish the same authoritative visible denominator.
	 * Only the bounded case needs the extra capacity/presentation wave below. */
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	if (coverageCompletion.missing) {
	    BObolViewLodState *coverageState =
		this->d->viewAttachment->getViewLodState();
	    if (coverageState)
		coverageState->setCadPresentationCameraMotionFrameReuse(
		    FALSE);
	}
    }
    bool capacityOwnsCompletedPass = false;
    if (coverageCompletion.completed && coverageCompletion.bounded) {
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD completed coverage pass visible=%zu "
		   "covered=%zu missing=%d admitted=%d "
		   "budget_blocked=%d rescan_pending=%d\n",
		   coverageCompletion.visibleCount,
		   coverageCompletion.coveredCount,
		   coverageCompletion.missing ? 1 : 0,
		   this->d->lodRetainedPass.admittedWork() ? 1 : 0,
		   this->d->lodRetainedPass.budgetBlocked() ? 1 : 0,
		   this->d->lodSubmissionPass.rescanPending() ? 1 : 0);
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodRetainedPass.clearAdmittedWork();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	this->d->retireRetainedRefinementObservation();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	if (coverageCompletion.missing) {
	    const bool residentGrowthOwnsSuccessor =
		BObolLodAvailabilityScheduler::residentGrowthOwnsSuccessor(
		    completedPassSuccessor);
	    /*
	     * Capacity for the next promotion wave must be learned from a
	     * presented frame, just like ordinary scene-budget admission.
	     * completeRenderTiming() restarts this coverage pass with the new
	     * allowance.  Until then, pending provider work may continue on the
	     * worker pool without turning the GUI loop into a busy rescan.
	     *
	     * A richer prefix which arrived during this pass is the stronger
	     * successor: it can satisfy the missing occurrence without increasing
	     * the scene budget.  Do not put a capacity-frame barrier in front of
	     * that resident-growth transaction.  The growth scheduler first
	     * rechecks coverage using the new prefix, then drains and reallocates
	     * the settled population.  Reversing those owners leaves each waiting
	     * for the other and alternates an unchanged coverage pass with an
	     * unchanged capacity sample.
	     */
	    this->d->lodSubmissionPass.deactivate();
	    if (residentGrowthOwnsSuccessor) {
		if (completedPassSuccessor ==
			BObolLodAvailabilityScheduler::CompletedPassSuccessor::
			    COMPLETE_RESIDENCY_DRAIN)
		    (void)this->d->lodAvailabilityLedger.
			completeResidencyDrain();
		this->markProgressiveWorkPending();
	    } else {
		this->d->requestCapacityRescan();
		this->requestRender("lod-coverage-admission");
	    }
	} else if (this->d->lodPresentationPolicy.handoffPending() &&
	    sceneLodState &&
	    sceneLodState->activeRenderCost() <=
		sceneLodState->minimumActiveRenderCost()) {
	    /* Coverage has established a useful presentation for every visible
	     * occurrence.  A deadline-created renderer ceiling may
	     * have protected the cold stream, but it also makes every refinement
	     * behind that ceiling look artificially cheap.  End the stale handoff
	     * and measure one ceiling-free minimum frame before spending any
	     * headroom on richer PoP prefixes. */
	    this->d->lodPresentationPolicy.cancelHandoff();
	    this->d->lodInteractiveProgressiveCeiling = -1;
	    sceneLodState->setCadPresentationProgressiveCutCeiling(-1);
	    this->d->lodSubmissionPass.deactivate();
	    this->d->requestCapacityRescan();
	    this->requestRender("lod-coverage-minimum-calibration");
	} else {
	    /* Every projected leaf has a useful structural or mesh presentation. Begin a
	     * fresh bounded pass which may now spend the remaining scene budget
	     * on screen-value-ordered PoP refinement. */
	    this->d->lodSubmissionPass.activate();
	    this->markProgressiveWorkPending();
	}
    } else if (BObolLodCoveragePolicy::needsDeferredQualitySuccessor(
	coverageCompletion.completed, coverageCompletion.missing,
	this->d->lodRetainedPass.refinementPending(),
	completedPassRetainedAllocation ||
	    this->d->lodPointQualityPhase.triangleRecoveryPending() ||
	    this->d->lodPresentationPolicy.handoffPending())) {
	/*
	 * Coverage-first suppression applies to every source, not only compact
	 * sources large enough to need several owner-thread windows.  A small or
	 * single-mesh source completes its coverage census in one action, but an
	 * already-present partial PoP prefix may still have been held back to
	 * preserve the same minimum-mesh-first contract.  Consume that explicit
	 * deferred-quality witness with one fresh pass now that the complete
	 * visible denominator is known.
	 *
	 * This pass begins from the retained prefix and the measured scene budget.
	 * It may retarget resident data immediately or schedule only the missing
	 * suffix; it does not rebuild geometry or return the occurrence to a box.
	 *
	 * A retained-allocation pass, point recovery, or presentation handoff is
	 * already that authoritative successor.  Routing its residual quality
	 * annotation through this branch starts an ordinary pass before the owner
	 * below can consume the allocation proof.  The ordinary pass then asks for
	 * another allocation, creating an allocation -> ordinary -> allocation
	 * loop over the complete scene.
	 */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.setActive(!sources.empty());
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	if (this->d->lodSubmissionPass.active())
	    this->markProgressiveWorkPending();
    } else if (retainedImportanceCensusCompleted) {
	/* Small sources do not take the bounded-coverage successor branch above,
	 * so explicitly start their one-shot importance allocation here. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.setActive(!sources.empty());
	this->d->lodRetainedPass.clearAdmittedWork();
	this->d->retireRetainedRefinementObservation();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	if (this->d->lodSubmissionPass.active())
	    this->markProgressiveWorkPending();
    } else if (completedPass && this->d->lodSubmissionPass.rescanPending()) {
	/*
	 * A rescan requested while consuming an exact source delta is a full
	 * current-view coverage obligation.  Handle it before refinement and
	 * budget barriers: either of those paths may legitimately wait for a
	 * frame, but neither may erase an inventory-coverage obligation.  The
	 * old ordering made cold completion depend on whether the final delta
	 * also changed a retained cut.
	 *
	 * Leaving delta mode armed reconstructs the old selective plan after
	 * clearLodSubmissionPlan() and can loop over that prefix forever while
	 * newly appended leaves remain boxes.
	 */
	const bool compactInventoryIncomplete =
	    controller_lod_compact_inventory_incomplete(sources);
	if (this->d->lodSubmissionDelta.active())
	    this->d->lodSubmissionDelta.reset();
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	if (compactInventoryIncomplete) {
	    /* The exact append delta is complete, but the producer says more
	     * entries are coming.  Preserve accumulated coverage and the final
	     * rescan obligation, then go idle on the LoD cursor until the next
	     * inventory revision.  The progressive source provider owns that wake
	     * edge.  Starting a full scan here repeatedly revisits [0,current) and
	     * dominates cold 50k/150k discovery. */
	    this->d->lodSubmissionPass.deactivate();
	} else {
	    if (this->d->lodCoveragePolicy.active())
		this->d->lodCoveragePolicy.clearPassCounters();
	    this->d->lodSubmissionPass.clearRescan();
	    this->d->lodSubmissionPass.setActive(!sources.empty());
	}
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    } else if (completedPass &&
	this->d->lodRetainedPass.refinementPending() &&
	this->d->lodRetainedPass.cutAdvanced() &&
	!this->d->lodRetainedPass.budgetBlocked()) {
	/* The just-completed pass selected the next resident cut using the
	 * newest view already.  Present it before any requested rescan; the
	 * post-frame submission is itself a full current-view pass.
	 *
	 * A budget-blocked retained allocation is deliberately excluded.  Its
	 * changed cut is also a calibration population: finishBlockedPass()
	 * must retain the post-frame reallocation obligation so measured
	 * headroom is consumed by another minimax pass.  Handling that case here
	 * used to present the better cut but then declare the scene idle with the
	 * larger calibrated allowance unused. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.deactivate();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	this->d->lodRetainedPass.clearAdmittedWork();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
    } else if (completedPassSuccessor ==
	    BObolLodAvailabilityScheduler::CompletedPassSuccessor::
		COMPLETE_RESIDENCY_DRAIN ||
	completedPassSuccessor ==
	    BObolLodAvailabilityScheduler::CompletedPassSuccessor::
		YIELD_TO_RESIDENT_GROWTH) {
	/* A completed pass cannot calibrate capacity while resident growth owns
	 * the shared cursor.  An active drain advances its typed transaction; an
	 * ordinary pass merely yields so the progressive pump can start that
	 * drain.  Both successors preserve the growth debt and leave exactly one
	 * scheduler responsible for the next pass. */
	const bool completeResidencyDrain = completedPassSuccessor ==
	    BObolLodAvailabilityScheduler::CompletedPassSuccessor::
		COMPLETE_RESIDENCY_DRAIN;
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.deactivate();
	this->d->resetRetainedPassAnnotations();
	if (completeResidencyDrain)
	    (void)this->d->lodAvailabilityLedger.completeResidencyDrain();
    } else if (BObolLodAvailabilityScheduler::
	residentRefinementOwnsSuccessor(completedPassSuccessor)) {
	/* A resident suffix is availability work, not another capacity sample.
	 * Reallocating the unchanged occurrence population while its task is in
	 * flight spins the owner thread and cannot make a cut drawable.  Yield to
	 * the result edge when a producer exists; otherwise run one ordinary pass
	 * which can submit the missing request. */
	const bool submitResidentRequests = completedPassSuccessor ==
	    BObolLodAvailabilityScheduler::CompletedPassSuccessor::
		SUBMIT_RESIDENT_REQUESTS;
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.setActive(
	    submitResidentRequests && !sources.empty());
	this->d->resetRetainedPassAnnotations();
	if (this->d->lodSubmissionPass.active()) {
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::AVAILABILITY);
	    this->markProgressiveWorkPending();
	}
    } else if (completedPassSuccessor ==
	BObolLodAvailabilityScheduler::CompletedPassSuccessor::
	    PRESENT_POINT_CALIBRATION) {
	/* A capacity search which was already active may finish before a point
	 * calibration consumes its frame, but it may not reproduce itself behind
	 * that waiting owner.  The completed reallocation is now the exact
	 * population which the point frame must measure. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.deactivate();
	this->d->resetRetainedPassAnnotations();
	this->markProgressiveWorkPending();
	this->requestRender("lod-point-calibration-successor");
    } else if (completedPassSuccessor ==
	BObolLodAvailabilityScheduler::CompletedPassSuccessor::
	    CALIBRATE_CAPACITY) {
	/* A finite scene budget is a stable policy, not a timer to an eventual
	 * unbounded frame.  If this pass admitted a first wave, re-plan only
	 * after that wave renders and calibrates throughput.
	 *
	 * A quiet cut with substantial measured frame headroom gets another
	 * unchanged calibration frame even when this pass admitted nothing.
	 * The stable throughput EMA otherwise receives only one sample at each
	 * population and can stop at an artificially coarse fixed point (for
	 * example Generic Twin's tail at ~50k faces despite a 25 ms software
	 * frame and a 50 ms stable target).  These probes are bounded by the
	 * actual stable frame target and never relax the pixel-error demand. */
	const uint64_t stableTargetNanoseconds =
	    this->d->quietAllocationTargetFps() > 0.0f ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(
		    this->d->quietAllocationTargetFps())) : 0;
	const uint64_t preferredStableTargetNanoseconds =
	    this->d->lodStableTargetFps > 0.0f ?
	    static_cast<uint64_t>(1000000000.0L /
		static_cast<long double>(this->d->lodStableTargetFps)) : 0;
	/* A pose-only orthographic handoff begins with the occurrence cuts from a
	 * previously ready view.  Its first quiet frame is therefore a static
	 * capacity probe, not a cold-start search for the preferred redraw cadence.
	 * Searching the preferred goal first can only coarsen that retained
	 * population and then rebuild it under the independent static deadline.
	 * Use the static goal directly; a missed hard deadline still performs the
	 * ordinary bounded downward search, while a completed frame becomes a safe
	 * lower bound and permits only richer successors. */
	const bool retainedPosePresentation =
	    this->d->lodPoseContinuity.retainOccurrenceCuts();
	const uint64_t capacitySearchPreferredNanoseconds =
	    this->d->lodStaticQualityTrial.usesStaticDeadline() ||
	    retainedPosePresentation ?
		this->d->prominentQualityPresentationDeadline() :
		preferredStableTargetNanoseconds;
	/* Stable pass decisions require a duration for this retained cut.  GPU
	 * query latency is handled by the paired throughput calibration above;
	 * reusing that duration after the occurrence population changes is not a
	 * current-cut deadline measurement. */
	const uint64_t observedStableNanoseconds =
	    this->d->lodLastRenderWasReusableCadPresentation ?
		(this->d->lastRenderTimeNanoseconds ?
		    this->d->lastRenderTimeNanoseconds :
		    this->d->smoothedRenderTimeNanoseconds) : 0;
	/* A probe is useful only while it can establish a larger affordable cut.
	 * Three unchanged samples reject one-time setup noise.  Capacity growth is
	 * then tested by the next distinct admitted population; repeating the same
	 * cut further cannot prove that a richer cut is affordable. */
	/* Never classify a newly reached quiet cut from a single frame.  OSMesa
	 * in particular can pay shader/cache/setup work in the first traversal,
	 * while a compositor stall can perturb System GL.  Three samples remain
	 * a small bounded cost at the stable target and let subsequent ordinary
	 * traversal cost govern the retained budget. */
	BObolLodCapacityEvidence::CalibrationInputs calibrationInputs;
	calibrationInputs.observedNanoseconds = observedStableNanoseconds;
	calibrationInputs.passAdmittedWork =
	    this->d->lodRetainedPass.admittedWork();
	const BObolRetainedAllocationResult &capacityAllocation =
	    this->d->lodRetainedAllocationCertificate;
	const uint64_t activeAllocationPlan = sceneLodState ?
	    sceneLodState->activeCadAllocationPlan() : 0;
	const bool allocationPlanCurrent =
	    capacityAllocation.allocationPlanSerial != 0 &&
	    capacityAllocation.allocationPlanSerial == activeAllocationPlan;
	const bool allocationRevisionCurrent =
	    capacityAllocation.viewRevision == this->d->lodViewRevision.value() &&
	    capacityAllocation.policyRevision == this->d->lodPolicyRevision.value();
	const bool allocationThresholdCurrent =
	    std::isfinite(capacityAllocation.pointProxyPixelThreshold) &&
	    std::fabs(capacityAllocation.pointProxyPixelThreshold -
		this->d->lodPresentationPointProxyPixelThreshold) <= 1.0e-6f;
	const bool allocationCostsValid =
	    capacityAllocation.selectedPresentationCost > 0 &&
	    capacityAllocation.pixelDemandPresentationCost > 0;
	const bool allocationBudgetCurrent =
	    capacityAllocation.requestedSceneBudget ==
		this->d->lodAdmissionEvidence.capacity().currentBudget();
	const bool allocationBudgetCertified =
	    allocationBudgetCurrent &&
	    capacityAllocation.certifiedPresentationBudget >=
		capacityAllocation.selectedPresentationCost;
	const bool allocationPopulationCurrent = sceneLodState &&
	    sceneLodState->cadAllocationPlanCoversCurrentPopulation(
		capacityAllocation.allocationPlanSerial,
		capacityAllocation.viewRevision,
		capacityAllocation.policyRevision,
		capacityAllocation.fixedCadPresentationCost);
	const bool allocationCutsApplied =
	    this->d->retainedAllocationCutsApplied(sceneLodState);
	const bool capacityBudgetFinite =
	    this->d->lodAdmissionEvidence.capacity().currentBudget() != SIZE_MAX;
	/* A renderer-wide ceiling may be the only difference between the richer
	 * retained occurrence cuts and the allocation currently on screen.  The
	 * allocation certificate is realized when that exact ceiling presents its
	 * selected cost; comparing only with the hidden retained cost makes a
	 * completed, drawable plan look perpetually stale.  The resulting
	 * zero-work capacity pass otherwise resubmits forever because no occurrence
	 * mutation can make those two currencies equal while the ceiling remains
	 * authoritative. */
	const bool allocationPresentationRealized =
	    this->d->retainedAllocationPresentationRealized(sceneLodState);
	const size_t allocationPresentationCost =
	    allocationPresentationRealized ?
		capacityAllocation.selectedPresentationCost : sceneActiveCost;
	const bool overBudgetAllocationClaimed =
	    this->d->lodPresentationPolicy.claimOverBudgetAllocation(
		allocationPlanCurrent && allocationRevisionCurrent &&
		    allocationThresholdCurrent && allocationCostsValid &&
		    allocationBudgetCurrent && allocationPopulationCurrent &&
		    capacityBudgetFinite,
		allocationCutsApplied,
		capacityAllocation.selectedPresentationCost,
		capacityAllocation.certifiedPresentationBudget,
		allocationPresentationCost);

	/* Select the effect-producing owner from one immutable completed-pass
	 * snapshot.  In particular, a current occurrence allocation may prove that
	 * a pending capacity sample first needs a ceiling-free handoff.  That
	 * normalization must happen before capacity calibration can alter evidence
	 * or install its own retained-allocation request. */
	BObolLodPresentationPolicy::CompletedPassInputs capacityHandoffInputs;
	capacityHandoffInputs.completed = completedPass != FALSE;
	/* The source cursor has reached its end.  Its mechanical retirement below
	 * is part of this transaction and is not a successor submission. */
	capacityHandoffInputs.submissionPending = false;
	capacityHandoffInputs.rescanAfterFrame =
	    this->d->lodAdmissionEvidence.capacity().rescanAfterFrame();
	capacityHandoffInputs.changedCut = completedPassChangedCut != FALSE;
	const bool handoffServiceQuiescent = !service ||
	    (service->activeTaskCountForGeneration(generation) == 0 &&
	     service->queuedResultCountForGeneration(generation) == 0 &&
	     this->d->lodAvailabilityLedger.resultsPending() == 0);
	const bool handoffPopulationQuiescent =
	    handoffServiceQuiescent && retainedPopulationSettled &&
	    !this->d->lodPresentationTransaction.barrierPending() &&
	    !this->d->lodPresentationTransaction.publicationPending() &&
	    !this->d->lodAvailabilityLedger.residentGrowthPending();
	capacityHandoffInputs.populationQuiescent = handoffPopulationQuiescent;
	const size_t handoffReconciliationBudget =
	    this->d->lodPresentationPolicy.handoffReconciliationBudget();
	const bool allocationCertificateCurrent =
	    this->d->retainedAllocationCertificateCurrent(sceneLodState);
	capacityHandoffInputs.retainedAllocationCompleted =
	    allocationCertificateCurrent;
	capacityHandoffInputs.retainedAllocationCertified = allocationCutsApplied;
	capacityHandoffInputs.presentationLimitsReconciled =
	    BObolLodPresentationPolicy::presentationLimitsReconciled(
		capacityHandoffInputs.retainedAllocationCompleted,
		allocationCutsApplied,
		capacityAllocation.selectedPresentationCost,
		capacityAllocation.certifiedPresentationBudget,
		handoffReconciliationBudget);
	capacityHandoffInputs.retainedRefinementPending =
	    completedPassRefinementPending;
	capacityHandoffInputs.retainedRefinementBudgetBlocked =
	    completedPassBudgetBlocked;
	const BObolLodPresentationPolicy::CompletedPassSelection
	    capacityPassSelection =
		this->d->lodPresentationPolicy.completedPassSelection(
		    capacityHandoffInputs, !overBudgetAllocationClaimed,
		    this->d->lodAdmissionEvidence.capacity().capacitySearch().
			awaitingSample(),
		    this->d->lodInteractiveProgressiveCeiling);
	capacityOwnsCompletedPass = capacityPassSelection.capacityOwns();
	const BObolLodPresentationPolicy::CompletedPassOwner completedPassOwner =
	    capacityPassSelection.owner;
	/* The allocation certificate belongs to the applied occurrence plan, not
	 * to the mechanical pass which happened to create it.  An unchanged
	 * presentation/barrier pass must continue the same bounded search using
	 * that certificate.  Requiring this pass to recompute the allocation made
	 * a terminal search alternate forever with a redundant reallocation. */
	const bool currentCapacityAllocation =
	    allocationPlanCurrent &&
	    allocationRevisionCurrent && allocationThresholdCurrent &&
	    allocationCostsValid &&
	    allocationBudgetCertified && allocationCutsApplied &&
	    capacityBudgetFinite;
	calibrationInputs.allocationCutsApplied = currentCapacityAllocation;
	calibrationInputs.allocationPresentationRealized =
	    allocationPresentationRealized;
	calibrationInputs.activeCost =
	    BObolLodAdmissionPlanner::capacitySearchPresentedCost(
		allocationPresentationRealized, sceneActiveCost,
		capacityAllocation.selectedPresentationCost);
	calibrationInputs.boundedSearch = currentCapacityAllocation;
	calibrationInputs.candidateBudget = currentCapacityAllocation ?
	    capacityAllocation.certifiedPresentationBudget : 0;
	calibrationInputs.populationSignature = currentCapacityAllocation ?
	    capacityAllocation.selectedPopulationSignature : 0;
	calibrationInputs.demandedBudget = currentCapacityAllocation ?
	    capacityAllocation.pixelDemandPresentationCost : 0;
	/* The lower bracket belongs to the certificate's first target.  Seed only
	 * from this exact completed frame when its measured duration proves that
	 * target directly; evidence from a different deadline is not transferable. */
	calibrationInputs.knownSafeBudget =
	    allocationPresentationRealized && observedStableNanoseconds > 0 &&
	    observedStableNanoseconds <= capacitySearchPreferredNanoseconds ?
		calibrationInputs.activeCost : 0;
	calibrationInputs.searchKey =
	    BObolLodCapacitySearchCertificate::keyFor(
		this->d->admissionRevisionStamp(),
		capacitySearchPreferredNanoseconds,
		currentCapacityAllocation &&
			!this->d->lodAdmissionEvidence.resources().anyPressure() ?
		    this->d->prominentQualityPresentationDeadline() :
		    capacitySearchPreferredNanoseconds,
		calibrationInputs.demandedBudget,
		sceneMinimumActiveCost);
	calibrationInputs.searchKey.preferredBudgetCeiling =
	    this->d->lodAdmissionEvidence.capacity().
		deadlineCapacityCeiling();
	calibrationInputs.searchKey.maximumBudgetCeiling =
	    this->d->lodAdmissionEvidence.capacity().
		staticDeadlineCapacityCeiling();
	/* An applied over-budget certificate is already the result of the
	 * capacity candidate named by this pass.  Its typed handoff below owns the
	 * only valid successors: local point aggregation, one bounded static
	 * attempt, or a terminal constraint.  Sending it through the generic
	 * fallback merely requests the same allocation again. */
	const BObolLodCapacityEvidence::CalibrationDecision calibration =
	    completedPassOwner !=
		BObolLodPresentationPolicy::CompletedPassOwner::CAPACITY ?
		    BObolLodCapacityEvidence::CalibrationDecision() :
		this->d->finishBlockedCapacityPass(calibrationInputs);
	/* A structural repair is already ordered by projected importance.  If its
	 * finite scene allowance cannot replace the entire exact box frontier,
	 * preserve the admitted prominent meshes and aggregate the remaining
	 * small-object tail before another pass.  This avoids tens of thousands of
	 * individually affordable frames while retaining immutable resident data
	 * for later views and faster renderers. */
	BObolLodPointProxyEvidence::Decision structuralAggregation;
	const size_t structuralTailBlockedCount =
	    completedStructuralRepair ?
		completedMissingMeshBudgetBlockedCount : 0;
	const float structuralPointThresholdBefore =
	    this->d->lodPresentationPointProxyPixelThreshold;
	if (capacityPassSelection.capacityOwns() &&
	    structuralTailBlockedCount > 0 &&
	    !this->d->lodInteractionSession.active() &&
	    /* Structural coverage is discovered before the occurrence batch is
	     * necessarily installed.  Point aggregation only changes that batch;
	     * arming its completed-frame barrier before one exists blocks the very
	     * source admission which creates it.  Hidden-line and other direct
	     * presentation paths exposed this as an infinite "calibrating
	     * small-part aggregation" loop containing only structural boxes. */
	    sceneLodState && sceneLodState->hasCadPresentationAssemblies() &&
	    this->d->pointProxyAggregationApplicable() &&
	    this->d->structuralPointAggregationRequiredForBudget(
		this->d->lodAdmissionEvidence.capacity().currentBudget())) {
	    const BObolLodAdmissionPlan structuralPlan =
		BObolLodAdmissionPlanner::planPointStructuralCoverageBlocked(
		    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    structuralTailBlockedCount);
	    this->d->commitAdmissionPlan(structuralPlan);
	    structuralAggregation = structuralPlan.pointProxyDecision;
	    if (structuralAggregation.changed) {
		this->d->lodPresentationPointProxyPixelThreshold =
		    structuralAggregation.threshold;
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(
			    structuralAggregation.threshold);
		this->d->lodPointQualityPhase.requestCalibration();
		this->requestRender("lod-structural-tail-aggregation");
	    }
	}
	if (getenv("BOBOL_LOD_TRACE_STRUCTURAL_BUDGET") &&
	    structuralTailBlockedCount > 0)
	    bu_log("BObol LoD structural allocation tail=%zu "
		   "interactive=%d applicable=%d threshold=%.6g "
		   "next=%.6g changed=%d\n",
		   structuralTailBlockedCount,
		   this->d->lodInteractionSession.active() ? 1 : 0,
		   this->d->pointProxyAggregationApplicable() ? 1 : 0,
		   structuralPointThresholdBefore,
		   structuralAggregation.threshold,
		   structuralAggregation.changed ? 1 : 0);
	if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value())) {
	    static std::atomic<unsigned int> passTraceCount(0);
	    const unsigned int traceIndex = passTraceCount.fetch_add(1);
	    if (traceIndex < 256)
		bu_log("BObol LoD completed budget pass "
		       "owner=%u consume_annotations=%d "
		       "active_faces=%zu active_cost=%zu cost_budget=%zu admitted=%d "
		       "retained=%d pose_continuity=%d demanded=%zu "
		       "preferred_ms=%.3f maximum_ms=%.3f "
		       "candidate_reallocation=%d search_active=%d "
		       "sample_frame=%d restart=%d result=%u phase=%u goal=%u "
		       "samples_remaining=%u candidates=%u "
		       "observed_ms=%.3f target_ms=%.3f "
		       "calibrated_mcost_s=%.3f allocation_current=%d "
		       "plan_current=%d revision_current=%d threshold_current=%d "
		       "costs_valid=%d "
		       "budget_certified=%d population_current=%d cuts_applied=%d "
		       "presentation_realized=%d presentation_cost=%zu "
		       "budget_finite=%d selected=%zu pixel_demand=%zu "
		       "certified_budget=%zu requested_budget=%zu "
		       "plan=%llu active_plan=%llu\n",
		       static_cast<unsigned int>(capacityPassSelection.owner),
		       capacityPassSelection.consumePassAnnotations ? 1 : 0,
		       sceneActiveFaces, sceneActiveCost,
		       this->d->lodAdmissionEvidence.capacity().currentBudget(),
		       this->d->lodRetainedPass.admittedWork() ? 1 : 0,
		       completedPassRetainedAllocation ? 1 : 0,
		       retainedPosePresentation ? 1 : 0,
		       calibrationInputs.demandedBudget,
		       capacitySearchPreferredNanoseconds / 1000000.0,
		       calibrationInputs.searchKey.maximumTargetNanoseconds /
			   1000000.0,
		       calibration.candidateReallocation ? 1 : 0,
		       calibration.searchActive ? 1 : 0,
		       calibration.sampleFrame ? 1 : 0,
		       calibration.restartSubmission ? 1 : 0,
		       static_cast<unsigned int>(calibration.searchResult),
		       static_cast<unsigned int>(this->d->lodAdmissionEvidence.
			   capacity().capacitySearch().phase()),
		       static_cast<unsigned int>(this->d->lodAdmissionEvidence.
			   capacity().capacitySearch().goal()),
		       this->d->lodAdmissionEvidence.capacity().capacitySearch().
			   samplesRemaining(),
		       this->d->lodAdmissionEvidence.capacity().capacitySearch().
			   measuredCandidateCount(),
		       observedStableNanoseconds / 1000000.0,
		       stableTargetNanoseconds / 1000000.0,
		       static_cast<double>(
			   this->d->lodStableCalibratedRenderCostPerSecond /
			       1000000.0L),
		       currentCapacityAllocation ? 1 : 0,
		       allocationPlanCurrent ? 1 : 0,
		       allocationRevisionCurrent ? 1 : 0,
		       allocationThresholdCurrent ? 1 : 0,
		       allocationCostsValid ? 1 : 0,
		       allocationBudgetCertified ? 1 : 0,
		       allocationPopulationCurrent ? 1 : 0,
		       allocationCutsApplied ? 1 : 0,
		       allocationPresentationRealized ? 1 : 0,
		       allocationPresentationCost,
		       capacityBudgetFinite ? 1 : 0,
		       capacityAllocation.selectedPresentationCost,
		       capacityAllocation.pixelDemandPresentationCost,
		       capacityAllocation.certifiedPresentationBudget,
		       capacityAllocation.requestedSceneBudget,
		       static_cast<unsigned long long>(
			   capacityAllocation.allocationPlanSerial),
		       static_cast<unsigned long long>(activeAllocationPlan));
	}
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.deactivate();
	if (calibration.restartSubmission) {
	    /* The calibration decision starts a distinct pass.  None of the
	     * completed pass's cut, budget, residency, or work annotations may be
	     * inherited by that successor; doing so makes a clean handoff look as
	     * though every no-op retry changed a cut. */
	    this->d->resetRetainedPassAnnotations();
	    this->d->lodSubmissionPass.activate();
	    this->d->lodSubmissionPass.clearRescan();
	    this->markProgressiveWorkPending();
	}
	if (calibration.requestFrame) {
	    const char *frameReason = calibration.sampleFrame ?
		"lod-budget-sample" : "lod-population-barrier";
	    /* A large compact scan can outlive the render request installed by its
	     * first admitted result.  Publish one standing render level for the
	     * typed frame obligation which must complete before admission resumes. */
	    this->requestRender(frameReason);
	}
	/*
	 * The per-occurrence minimum PoP prefix is a real lower bound on the
	 * triangle allocator.  In a scene with tens of thousands of tiny parts
	 * that floor can remain several times larger than the measured renderer
	 * budget even after every admissible cut has been coarsened.  Do not
	 * declare that 100-300 ms presentation stable merely because another
	 * provider scan cannot lower it.  Raise the existing camera-local point
	 * aggregation threshold and measure the resulting one-call batch.
	 *
	 * This is presentation-only: immutable meshes and their desired cuts
	 * remain resident, so a later view or faster renderer restores them
	 * without cache I/O or geometry rebuilding.
	 */
	if (capacityPassSelection.capacityOwns() &&
	    !structuralAggregation.changed &&
	    !this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() &&
	    this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() &&
	    !this->d->lodInteractionSession.active() &&
	    this->d->retainedPointProxyAggregationApplicable() &&
	    sceneLodState && sceneLodState->hasCadPresentationAssemblies() &&
	    sceneActiveCost <= (sceneLodState ?
		sceneLodState->minimumActiveRenderCost() : sceneActiveCost) &&
	    sceneActiveCost > this->d->lodAdmissionEvidence.capacity().currentBudget() &&
	    stableTargetNanoseconds > 0 &&
	    observedStableNanoseconds >
		static_cast<long double>(stableTargetNanoseconds) * 1.20L) {
	    const BObolLodAdmissionPlan pressurePlan =
		BObolLodAdmissionPlanner::planPointInterrupted(
		    this->d->lodAdmissionEvidence, this->d->lodAdmissionCursor,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    static_cast<uint64_t>(observedStableNanoseconds),
		    this->d->quietAllocationTargetFps());
	    this->d->commitAdmissionPlan(pressurePlan);
	    const BObolLodPointProxyEvidence::Decision &pressure =
		pressurePlan.pointProxyDecision;
	    if (pressure.changed) {
		this->d->lodPresentationPointProxyPixelThreshold =
		    pressure.threshold;
		BObolViewLodState *presentationState =
		    this->d->viewAttachment->getViewLodState();
		if (presentationState)
		    presentationState->
			setCadPresentationPointProxyPixelThreshold(
			    pressure.threshold);
		this->d->lodPointQualityPhase.requestCalibration();
		/* The multiplicative pressure correction intentionally lands on
		 * the safe side of the target.  Its next unchanged frame continues
		 * the bounded bracket search, so terminal quality is the finest cut
		 * which meets the stable FPS contract. */
		this->requestRender("lod-stable-point-calibration");
	    }
	}
	append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
	    SbString(""),
	    this->d->lodPointQualityPhase.presentationPending() ?
		"scene LoD calibrating small-part aggregation" :
	    this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() ?
		(calibration.sampleFrame ?
		    "scene LoD measuring a bounded capacity candidate" :
		    "scene LoD presenting admitted geometry") :
		"scene LoD reached its calibrated face budget");
	this->d->lodRetainedPass.clearAdmittedWork();
	this->d->retireRetainedRefinementObservation();
	/* The capacity revision names a candidate population, not a completed
	 * mechanical scan.  Preserve it while a frame samples that candidate or
	 * while a terminal decision consumes its current allocation certificate.
	 * Advancing it after every no-op blocked pass invalidates the certificate
	 * which the handoff below must consume, reopening the same allocation and
	 * compaction forever. */
	if (calibration.restartSubmission)
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
    } else {
	this->d->lodSubmissionPass.setActive(!completedPass);
	if (completedPass)
	    this->d->lodRetainedPass.clearAdmittedWork();
	if (completedPass)
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    }

    /* A retained-admission allowance belongs to exactly one complete scene
     * pass.  Once that pass changes any cut, its carried remainder has been
     * consumed against the complete population selected by the minimax
     * ceiling.  No coverage/handoff successor is permitted to start another
     * all-occurrence allocation before that population is presented and the
     * policy pass is reset.  Doing so recomputed a fresh full-budget ceiling
     * while passing only the old few-unit remainder to the action; late
     * Hubble occurrences were consequently demoted from valid cut 4 to cut
     * 0. */
    if (completedPass && this->d->lodAdmissionCursor.retainedAdmission() &&
	completedPassChangedCut) {
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.deactivate();
    }

    /* A single-occurrence static-quality allocation is deliberately hidden
     * behind the last completed renderer ceiling.  Changing the occurrence
     * cut therefore need not alter the currently submitted prefix, so the
     * ordinary presentation barrier may have no reason to request a frame.
     * Wake one exact replay explicitly; its completion advances the bounded
     * staircase from the already measured baseline.  Keeping this edge
     * render-serial based also makes GUI-idle detection wait for the actual
     * quality transition instead of reporting ready with richer geometry
     * stranded behind the ceiling. */
    if (completedPass && completedPassRetainedAllocation &&
	this->d->lodStaticQualityTrial.probing() &&
	this->d->lodInteractiveProgressiveCeiling >= 0 && sceneLodState &&
	sceneLodState->maximumActiveProgressiveCut() >
	    this->d->lodInteractiveProgressiveCeiling &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodPresentationTransaction.barrierPending()) {
	this->markProgressiveWorkPending();
	this->requestRender("lod-static-overscan-allocation");
    }

    /*
     * Finish a motion-to-stable handoff only after a complete retained pass
     * made no further occurrence-cut changes.  Until this point the
     * renderer-wide ceiling is the only thing preventing the hidden retained
     * population from appearing in one unbounded frame.
     *
     * Keep the small-part aggregation threshold for the first ceiling-free
     * frame.  Its measured calibration below will relax or tighten that
     * threshold without touching the immutable mesh generations.
     */
    /* A retained-allocation commit proves the owner-thread CAD revision it
     * observed, but it cannot predict a provider result which is still in
     * flight.  Releasing the global deadline ceiling while that stream was
     * growing made each newly installed batch invalidate the measured
     * population, producing an exact constrained -> release -> deadline
     * cycle on 50k warm-cache scenes.  Finish handoff only after both the
     * immutable-result stream and its presentation/coalescing edges are
     * quiet.  The resident-growth policy below supplies the eventual
     * scene-wide retry rather than busy-scanning while workers run. */
    BObolLodPresentationPolicy::CompletedPassInputs handoffInputs;
    handoffInputs.completed = completedPass != FALSE;
    handoffInputs.submissionPending =
	this->d->lodSubmissionPass.active() != FALSE;
    handoffInputs.rescanAfterFrame =
	this->d->lodAdmissionEvidence.capacity().rescanAfterFrame();
    handoffInputs.changedCut = completedPassChangedCut != FALSE;
    const bool handoffServiceQuiescent = !service ||
	(service->activeTaskCountForGeneration(generation) == 0 &&
	 service->queuedResultCountForGeneration(generation) == 0 &&
	 this->d->lodAvailabilityLedger.resultsPending() == 0);
    handoffInputs.populationQuiescent =
	handoffServiceQuiescent && retainedPopulationSettled &&
	!this->d->lodPresentationTransaction.barrierPending() &&
	!this->d->lodPresentationTransaction.publicationPending() &&
	!this->d->lodAvailabilityLedger.residentGrowthPending();
    const size_t handoffReconciliationBudget =
	this->d->lodPresentationPolicy.handoffReconciliationBudget();
    const BObolRetainedAllocationResult &allocationCertificate =
	this->d->lodRetainedAllocationCertificate;
    const bool allocationCertificateCurrent =
	this->d->retainedAllocationCertificateCurrent(sceneLodState);
    const bool handoffAllocationCutsApplied =
	this->d->retainedAllocationCutsApplied(sceneLodState);
    handoffInputs.retainedAllocationCompleted =
	allocationCertificateCurrent;
    handoffInputs.retainedAllocationCertified =
	handoffAllocationCutsApplied;
    handoffInputs.presentationLimitsReconciled =
	BObolLodPresentationPolicy::presentationLimitsReconciled(
	    handoffInputs.retainedAllocationCompleted,
	    handoffAllocationCutsApplied,
	    allocationCertificate.selectedPresentationCost,
	    allocationCertificate.certifiedPresentationBudget,
	    handoffReconciliationBudget);
    handoffInputs.retainedRefinementPending =
	completedPassRefinementPending;
    handoffInputs.retainedRefinementBudgetBlocked =
	completedPassBudgetBlocked;
    BObolLodPresentationPolicy::CompletedPassSelection handoffSelection;
    if (!capacityOwnsCompletedPass) {
	handoffSelection =
	    this->d->lodPresentationPolicy.completedPassSelection(
		handoffInputs, false,
		this->d->lodAdmissionEvidence.capacity().capacitySearch().
		    awaitingSample(),
		this->d->lodInteractiveProgressiveCeiling);
	if (handoffSelection.consumePassAnnotations) {
	    handoffInputs.rescanAfterFrame = false;
	    handoffInputs.changedCut = false;
	}
    }
    BObolLodPresentationPolicy::CompletedPassDecision handoff;
    /* A capacity-owned completion may retire unrelated observation data, but
     * it may not run the handoff reducer after changing capacity evidence.
     * Selection already normalized the one case where the handoff must
     * precede a pending sample. */
    if (!capacityOwnsCompletedPass)
	handoff = this->d->lodPresentationPolicy.completePass(handoffInputs);
    if (completedPass && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_PASS", this->d->lodViewRevision.value())) {
	static std::atomic<unsigned int> handoffTraceCount(0);
	const unsigned int traceIndex = handoffTraceCount.fetch_add(1);
	if (traceIndex < 512)
	    bu_log("BObol LoD handoff pass owner=%u consume_annotations=%d "
		   "pending=%d rescan=%d changed=%d "
	       "reconciled=%d "
	       "allocation=%d certified=%d selected=%zu budget=%zu "
	       "quiescent=%d "
	       "retained=%d blocked=%d presentation_pending=%d "
	       "finish=%d request_allocation=%d request_rescan=%d "
		   "request_local_reduction=%d retire=%d\n",
		   static_cast<unsigned int>(capacityOwnsCompletedPass ?
		       BObolLodPresentationPolicy::CompletedPassOwner::CAPACITY :
		       handoffSelection.owner),
		   handoffSelection.consumePassAnnotations ? 1 : 0,
		   handoffInputs.submissionPending ? 1 : 0,
	       handoffInputs.rescanAfterFrame ? 1 : 0,
	       handoffInputs.changedCut ? 1 : 0,
	       handoffInputs.presentationLimitsReconciled ? 1 : 0,
	       handoffInputs.retainedAllocationCompleted ? 1 : 0,
	       handoffInputs.retainedAllocationCertified ? 1 : 0,
	       allocationCertificate.selectedPresentationCost,
	       allocationCertificate.certifiedPresentationBudget,
	       handoffInputs.populationQuiescent ? 1 : 0,
	       handoffInputs.retainedRefinementPending ? 1 : 0,
	       handoffInputs.retainedRefinementBudgetBlocked ? 1 : 0,
	       this->d->lodPresentationPolicy.handoffPresentationPending() ? 1 : 0,
	       handoff.finishHandoff ? 1 : 0,
	       handoff.requestRetainedAllocation ? 1 : 0,
		   handoff.requestRetainedRescan ? 1 : 0,
		   handoff.requestLocalPresentationReduction ? 1 : 0,
		   handoff.retireRetainedObservation ? 1 : 0);
    }

    const bool handoffResidentPopulationPending =
	!handoffServiceQuiescent ||
	this->d->lodAvailabilityLedger.resultsPending() != 0 ||
	this->d->lodPresentationTransaction.publicationPending();
    if (completedPass && completedPassRetainedAllocation &&
	this->d->lodPresentationPolicy.handoffPending() &&
	handoffResidentPopulationPending) {
	/* Coalesce every late completion/publication into one new allocation
	 * after the stream becomes idle.  This is a liveness witness as well as
	 * a debounce: without it the handoff can remain armed after the last
	 * worker result, while immediately restarting here wastes an O(scene)
	 * pass for every result batch.  A retained cut's own presentation barrier
	 * is not resident growth.  Treating that barrier as a population edge
	 * scheduled another allocation after every coherent frame and replayed a
	 * tiny set of cuts forever. */
	this->d->lodAvailabilityLedger.noteRicherPrefixAvailable();
	this->markProgressiveWorkPending();
    }

    if (handoff.requestRetainedAllocation) {
	/* A targeted source/result delta may be the last pass after the worker
	 * stream becomes quiet.  It cannot release the presentation handoff, but
	 * no external edge remains to request the authoritative allocation.  Start
	 * that transaction now and keep the existing constrained framebuffer on
	 * screen while its bounded owner-thread slices advance. */
	const size_t reconciliationBudget =
	    this->d->lodPresentationPolicy.handoffReconciliationBudget();
	if (reconciliationBudget > 0)
	    this->d->requestPresentationReconciliation(
		reconciliationBudget);
	else
	    this->d->requestRetainedReallocation();
	/* This is another pass in the same camera/policy/source epoch.  The
	 * explicit pending cursor is sufficient to bypass the completed-pass fast
	 * path.  Clearing the epoch witness as well makes the wrapper classify the
	 * already-pending cursor as a view change during submission and append an
	 * unnecessary full rescan after every allocation. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.setActive(!sources.empty());
	/* These latches summarize one complete bounded transaction.  The
	 * allocation request below starts a new transaction in the same epoch;
	 * carrying a preceding cut-advanced bit made an unchanged pass report
	 * changed and prevented handoff completion indefinitely. */
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	if (this->d->lodSubmissionPass.active())
	    this->markProgressiveWorkPending();
    }
    if (handoff.requestLocalPresentationReduction) {
	/* The occurrence allocator has proved that its current point/mesh split
	 * cannot replace the temporary global cut inside the handoff budget.
	 * Coarsen only the independently submitted, visually insignificant tail,
	 * then build another complete importance allocation.  Prominent/selected
	 * occurrences are excluded by the point-proxy classifier, so this retains
	 * their local PoP priority instead of uniformly clipping the whole scene. */
	size_t activePayloads = 0;
	size_t satisfiedPayloads = 0;
	size_t memoryLimitedPayloads = 0;
	if (sceneLodState)
	    sceneLodState->convergencePayloadCounts(activePayloads,
		satisfiedPayloads, memoryLimitedPayloads);
	const size_t unresolved = std::max<size_t>(
	    1, activePayloads > satisfiedPayloads ?
		activePayloads - satisfiedPayloads : activePayloads);
	BObolLodPointProxyEvidence::Decision reduction;
	auto restartHandoffAllocation = [&](size_t retryBudget) {
	    if (retryBudget > 0)
		this->d->requestPresentationReconciliation(retryBudget);
	    else
		this->d->requestRetainedReallocation();
	    this->d->lodSubmissionSourceIndex = 0;
	    this->d->lodSubmissionEntryOffset = 0;
	    this->d->clearLodSubmissionPlan();
	    this->d->lodSubmissionPass.clearRescan();
	    this->d->lodSubmissionPass.setActive(!sources.empty());
	    this->d->resetRetainedPassAnnotations();
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
	    if (this->d->lodSubmissionPass.active())
		this->markProgressiveWorkPending();
	};
	if (this->d->retainedPointProxyAggregationApplicable()) {
	    const BObolLodAdmissionPlan reductionPlan =
		BObolLodAdmissionPlanner::planPointStructuralCoverageBlocked(
		    this->d->lodAdmissionEvidence,
		    this->d->lodAdmissionCursor,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    unresolved);
	    this->d->commitAdmissionPlan(reductionPlan);
	    reduction = reductionPlan.pointProxyDecision;
	}
	if (reduction.changed) {
	    this->d->lodPresentationPointProxyPixelThreshold =
		reduction.threshold;
	    if (sceneLodState)
		sceneLodState->setCadPresentationPointProxyPixelThreshold(
		    reduction.threshold);
	    this->d->resetRetainedAdmissionQualityProof();
	    const size_t retryBudget = handoffReconciliationBudget > 0 ?
		handoffReconciliationBudget :
		allocationCertificate.certifiedPresentationBudget;
	    restartHandoffAllocation(retryBudget);
	} else {
	    /* The ordinary quiet budget is not the final fidelity limit for an
	     * event-driven static view.  If no small-part representation can reduce
	     * the occurrence-local minimum, try that minimum exactly once under the
	     * longer static deadline while retaining the completed constrained
	     * framebuffer. */
	    size_t presentedCost = handoffReconciliationBudget;
	    if (sceneLodState &&
		this->d->lodInteractiveProgressiveCeiling >= 0) {
		const size_t presentedCadCost = sceneLodState->
		    cadRenderCostAtProgressiveCutCeiling(
			this->d->lodInteractiveProgressiveCeiling);
		presentedCost = BObolLodAdmissionPlanner::
		    canonicalSceneCostAtCadCeiling(
			sceneLodState->activeRenderCost(),
			sceneLodState->activeCadRenderCost(),
			presentedCadCost);
	    }
	    const size_t staticBudget = BObolLodQualityPolicy::
		staticLocalMinimumRetryBudget(
		    presentedCost,
		    allocationCertificate.selectedPresentationCost,
		    allocationCertificate.pixelDemandPresentationCost,
		    this->d->lastRenderTimeNanoseconds,
		    this->d->prominentQualityPresentationDeadline());
	    if (!this->d->lodStaticQualityTrial.blocksNewTrial() &&
		staticBudget > handoffReconciliationBudget) {
		this->d->lodStaticQualityTrial.begin(
		    this->d->lodPresentationPointProxyPixelThreshold);
		this->d->lodPresentationPolicy.armHandoff(
		    false, presentedCost, staticBudget);
		restartHandoffAllocation(staticBudget);
		if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
			this->d->lodViewRevision.value()))
		    bu_log("BObol LoD static local-minimum trial "
			   "presented=%zu selected=%zu budget=%zu "
			   "elapsed_ms=%.3f deadline_ms=%.3f\n",
			   presentedCost,
			   allocationCertificate.selectedPresentationCost,
			   staticBudget,
			   this->d->lastRenderTimeNanoseconds / 1000000.0,
			   this->d->prominentQualityPresentationDeadline() /
			       1000000.0);
	    } else {
		/* The bounded static attempt also could not encode a complete local
		 * population.  Preserve the last completed framebuffer and wait for a
		 * real view, resource, or user-policy edge instead of retrying it. */
		append_controller_lod_diagnostic(
		    this->d->lastLodDiagnostics,
		    SbString("retained-allocation"),
		    "protected visible minimum exceeds the measured static frame "
		    "budget; retaining the last completed frame");
		if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
			this->d->lodViewRevision.value()))
		    bu_log("BObol LoD static overscan local reduction saturated "
			   "active=%d rejected=%d ceiling=%d\n",
			   this->d->lodStaticQualityTrial.inProgress() ? 1 : 0,
			   this->d->lodStaticQualityTrial.capacityRejected() ? 1 : 0,
			   this->d->lodInteractiveProgressiveCeiling);
		this->d->lodPresentationPolicy.cancelHandoff();
		BObolLodStaticQualityTrial::Constraint constraint;
		constraint.reason =
		    BObolLodStaticQualityTrial::ConstraintReason::
			PROTECTED_MINIMUM;
		constraint.revisionStamp = this->d->admissionRevisionStamp();
		constraint.committedCeiling =
		    this->d->lodInteractiveProgressiveCeiling;
		constraint.candidateCeiling =
		    this->d->lodInteractiveProgressiveCeiling;
		constraint.candidateCost = std::max<size_t>(1,
		    allocationCertificate.selectedPresentationCost);
		constraint.allowedCost = handoffReconciliationBudget;
		(void)this->d->lodStaticQualityTrial.reject(constraint);
	    }
	}
    }
    /* The policy distinguishes an unsatisfied quality observation from an
     * actual progress witness.  Retiring only the former prevents a
     * performance-limited terminal cut from becoming an endless compaction
     * defer/BACKGROUND loop.  A later result, camera/policy epoch, or capacity
     * edge starts a fresh pass from the retained cut. */
    if (handoff.retireRetainedObservation)
	this->d->lodRetainedPass.clearRefinementPending();
    if (handoff.finishHandoff) {
	BObolViewLodState *presentationState =
	    this->d->viewAttachment->getViewLodState();
	const bool completedStaticReconciliation =
	    this->d->lodStaticQualityTrial.completeReconciliation();
	const bool staticPresentationReconciliationPending =
	    this->d->lodStaticQualityTrial.probing() &&
	    this->d->lodInteractiveProgressiveCeiling >= 0 &&
	    presentationState &&
	    presentationState->maximumActiveProgressiveCut() >
		this->d->lodInteractiveProgressiveCeiling;
	/* Retained occurrence cuts and their immutable PoP suffix are the
	 * continuity baseline across camera epochs.  The former normalization
	 * path rewrote every cut to its minimum merely because a renderer-only
	 * ceiling had protected one handoff frame.  It then walked the same
	 * resident data upward again, producing the conspicuous coarse flash after
	 * every zoom and invalidating otherwise reusable command records.
	 *
	 * Remove only the reversible ceiling and measure the existing cut when the
	 * completed pass actually reconciled the retained occurrences.  The
	 * revision-bound reconciliation budget prevents a deadline recovery from
	 * treating its cheaper constrained frame as permission to retry the richer
	 * cut immediately.
	 *
	 * The endpoint render remains deadline-bounded and double buffered.  An
	 * explicit later capacity/view edge can relax a preserved presentation cut;
	 * memory-pressure compaction remains the sole authority for shrinking the
	 * retained occurrence prefixes. */
	if (completedStaticReconciliation ||
	    !staticPresentationReconciliationPending) {
	    this->d->lodInteractiveProgressiveCeiling = -1;
	    if (presentationState)
		presentationState->setCadPresentationProgressiveCutCeiling(-1);
	}
	/* A rejected richer cut completes the search, but the predecessor is still
	 * a valid event-driven static frame.  Retain that phase after occurrence
	 * reconciliation so unrelated HUD, selection, or screenshot redraws use
	 * the same deadline which proved it.  Camera/resource/capacity transitions
	 * explicitly invalidate the phase elsewhere. */
	if (presentationState)
	    presentationState->setCadPresentationCameraMotionFrameReuse(FALSE);
	if (presentationState &&
	    presentationState->hasCadPresentationAssemblies() &&
	    this->d->lodPresentationPointProxyPixelThreshold > 1.01f)
	    this->d->lodPointQualityPhase.requestHandoffConfirmation();
	else
	    this->d->lodPointQualityPhase.completeCalibration();
	/* The motion-time renderer ceiling may have prevented an otherwise
	 * resident occurrence from advancing at all.  That leaves a real richer
	 * demand but no changed cut, worker, or result which could schedule the
	 * next pass.  Preserve the bounded first ceiling-free presentation, then
	 * turn that demand into the existing post-frame rescan witness.  Leaving
	 * the raw retained-pending annotation armed here creates an unwitnessed
	 * progressive/compaction loop; submitting immediately would skip the
	 * handoff frame and could expose an unbounded population. */
	if (handoff.requestRetainedRescan) {
	    this->d->lodRetainedPass.clearRefinementPending();
	    this->d->requestCapacityRescan();
	}
	this->requestRender(staticPresentationReconciliationPending ?
	    "lod-static-overscan-reconcile" : "lod-stable-handoff");
    }

    /* A retained-prefix recovery can discover that the current occurrence
     * population is already at its minimum cut.  No geometry changes in that
     * case, so no presentation barrier is requested and completeRenderTiming()
     * will not be called merely to retire this planning latch.  Complete the
     * no-op recovery here; recoveries which did change a cut still pass
     * through the render barrier below and are retired by the measured-frame
     * path. */
    if (completedPass && !completedPassChangedCut)
	(void)this->completePointTriangleRecoveryIfReady();
    const bool recoveryPresentationRequired =
	this->d->lodPointQualityPhase.requiresRecoveryPresentation(
	    completedPass, completedPassChangedCut,
	    this->d->lodSubmissionPass.active() != FALSE,
	    this->d->lodPresentationTransaction.barrierPending());
    if (!this->d->lodSubmissionPass.active() &&
	(this->d->lodRetainedPass.cutAdvanced() ||
	 recoveryPresentationRequired) &&
	(boundedScenePass || recoveryPresentationRequired ||
	 this->d->lodRetainedPass.refinementPending() ||
	 this->d->lodPresentationPolicy.handoffPending())) {
	/* Richer prefixes selected by the completed pass are already in memory;
	 * expose that coherent budgeted wave in one completed frame.  This
	 * presentation barrier is also required when
	 * the just-completed pass reached every requested cut: in that case
     * retained-pass refinement pending is false, but the motion-to-stable
	 * handoff still needs one frame with the newly selected cuts before a
	 * second unchanged pass may remove its renderer-wide ceiling.  Gating
	 * this on "pending" stranded the handoff forever at large-scene scale.
	 * Retargeting is metadata/draw-count only; no provider task, cache read,
	 * or geometry rebuild is involved. */
	this->d->retireRetainedRefinementObservation();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->scheduleLodRefinementFrame("lod-cut");
    }
    /* The bounded calibration series may end before its first unchanged,
     * reusable CAD replay.  Request one explicit terminal probe now, while
     * this pass still owns the progress edge.  Its exact view/policy/budget
     * witness is one-shot; no later selection or HUD repaint may substitute
     * for this frame or reopen a view which already reported stable. */
    if (completedPass)
	this->armStableLodHeadroomProbeIfReady();
    if (completedPass && !this->d->lodSubmissionPass.active())
	this->d->lodDiscretePopulationTrialPermit.revoke();
    if (completedPass)
	this->d->lodSubmissionIntent.setRetainedAdmission(false);
    if (this->d->lodSubmissionPass.active())
	this->markProgressiveWorkPending();

    /* Large compact passes are intentionally consumed in small GUI work
	 * slices.  Rebuilding the complete command plan after every 3 ms slice
	 * turns a sub-second metadata scan into minutes of interrupted frames.
	 * The accumulated cut-advanced latch above installs one coherent barrier
	 * after the pass completes. */
    if (this->d->lastLodUpdatedCutCount > 0 &&
	(!boundedScenePass || completedPass)) {
	this->d->lodCompactionPolicy.requestAfter(bu_gettime(), 750000);
	this->requestRender("lod-cut");
    }
    if (completedPass && this->d->lodSubmissionDelta.active() &&
	!this->d->lodSubmissionPass.active() &&
	!this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() &&
	!this->d->lodPresentationTransaction.barrierPending()) {
	this->d->lodSubmissionDelta.reset();
    }

    return size_to_int_saturated(
	       static_cast<size_t>(this->d->lastLodSubmittedTaskCount));
}

int
BObolViewController::applyLodResults(BObolLodService *service,
			       size_t maxResults,
			       size_t maxEstimatedBytes,
			       uint64_t generation)
{
    if (!service)
	service = this->d->lodService;

    this->d->lastLodResultCount = 0;
    this->d->lastLodMatchedResultCount = 0;
    this->d->lastLodAppliedResultCount = 0;
    this->d->lastLodRejectedResultCount = 0;
    this->d->lastLodUnmatchedResultCount = 0;
    this->d->lastLodDiagnostics = "";

    if (!service) {
	this->d->lastLodDiagnostics = "LoD service is not configured";
	return -1;
    }

    const uint64_t drainGeneration = generation != 0 ? generation :
	this->d->lodActiveGeneration;
    if (drainGeneration == 0) {
	this->d->lastLodDiagnostics =
	    "LoD result drain requires an owned generation";
	return 0;
    }

    std::vector<BObolLodResult> drained;
    this->d->lastLodResultCount = service->drainGenerationResults(
	drained, drainGeneration,
	maxResults, maxEstimatedBytes);
    this->d->lodAvailabilityLedger.setResultsPending(
	service->queuedResultCountForGeneration(drainGeneration) > 0);
    if (this->d->lastLodResultCount == 0)
	return 0;

    BObolViewLodState *viewState =
	this->d->viewAttachment->getViewLodState();
    SoBRLLodUpdateAction update;
    update.setViewLodState(viewState);
    /*
     * Compact results already carry both the retained source token and their
     * exact occurrence key.  Resolve the token through the scene controller's
     * structural index and mutate the view payload map directly.  The index
     * is rebuilt only after a genuine structural scene mutation; the steady
     * box->mesh and PoP-cut path is O(1) per result.
     *
     * Keep the action only for non-compact SoBRLMeshShape results.  There is
     * no extra traversal cost on the compact CAD path and no duplicated source
     * geometry or presentation state.
     */
    unsigned int directMatched = 0;
    unsigned int directApplied = 0;
    unsigned int directRejected = 0;
    unsigned int directRetainedPresentationCount = 0;
    SbBool partialRefinementCandidate = FALSE;
    SbBool residentGrowthAdmissionCandidate = FALSE;
    /*
     * GED views render the shared CAD scene plus a view-local scene.  The
     * view controller's sceneController indexes only the latter, so asking it
     * to resolve a shared sourceRoutingId misses and sends every result
     * through a full action traversal.  Build the tiny top-level source route
     * table once for this drained batch.  Source nodes terminate collection,
     * so this does not walk their compact occurrence populations.
     */
    std::unordered_map<uint64_t, SoBRLDatabaseSource *> renderSourceRoutes;
    {
	const std::vector<SoBRLDatabaseSource *> renderSources =
	    controller_render_database_source_roots(this);
	renderSourceRoutes.reserve(renderSources.size());
	for (SoBRLDatabaseSource *source : renderSources) {
	    if (!source)
		continue;
	    const uint64_t routingId = source->getCompactSourceRoutingId();
	    if (routingId)
		renderSourceRoutes[routingId] = source;
	}
    }
    for (size_t i = 0; i < drained.size(); i++) {
	/*
	 * Resolve a compact result's retained source token once.  In particular,
	 * stale-epoch arbitration must not use the source-agnostic generic
	 * lookup: that lookup scans every resident occurrence and made a burst of
	 * 2,048 warm results quadratic on a 50,000-leaf scene.
	 */
	SoBRLDatabaseSource *route = NULL;
	if (drained[i].request.sourceRoutingId != 0 &&
	    drained[i].request.occurrenceKey.getLength() > 0) {
	    const auto renderRoute = renderSourceRoutes.find(
		drained[i].request.sourceRoutingId.value());
	    if (renderRoute != renderSourceRoutes.end())
		route = renderRoute->second;
	    else
		route = this->d->sceneController.findDatabaseSourceRoutingId(
		    drained[i].request.sourceRoutingId.value());
	}
	if (getenv("BOBOL_LOD_TRACE_REJECTIONS") &&
	    drained[i].request.sourceRoutingId != 0 &&
	    drained[i].request.occurrenceKey.getLength() > 0 &&
	    !route) {
	    static std::atomic<unsigned int> routeMissTraceCount(0);
	    const unsigned int traceIndex =
		routeMissTraceCount.fetch_add(1);
	    if (traceIndex < 64)
		bu_log("BObol LoD rejection route-miss object=%s "
		       "occurrence=%s routing=%llu\n",
		       drained[i].request.objectName.getString(),
		       drained[i].request.occurrenceKey.getString(),
		       static_cast<unsigned long long>(
			   drained[i].request.sourceRoutingId.value()));
	}
	const auto findResidentCad = [&]() {
	    if (!viewState)
		return static_cast<const BObolViewLodState::CadPayload *>(NULL);
	    if (drained[i].request.sourceRoutingId != 0) {
		if (!route || drained[i].request.sourcePopulationEpoch == 0 ||
		    drained[i].request.sourcePopulationEpoch.value() !=
			route->getCompactPopulationEpoch())
		    return static_cast<const
			BObolViewLodState::CadPayload *>(NULL);
		return viewState->findCadForResult(route, drained[i]);
	    }
	    return viewState->findCadForResult(drained[i]);
	};
	const BObolViewLodState::CadPayload *residentCadBefore =
	    findResidentCad();
	const BObolViewLodState::MeshPayload *residentMeshBefore = viewState ?
	    viewState->findMeshForResult(drained[i]) : NULL;
	/*
	 * A resident PoP asset is cumulative and independent of the camera
	 * epoch which asked the worker to append its newest suffix.  Continuous
	 * wheel/trackpad input can advance the view several times while one large
	 * suffix is being read.  Rejecting that completed generation merely
	 * because the occurrence metadata has already been retargeted to the
	 * newest view makes refinement possible only after input stops.
	 *
	 * Rebase only an unambiguously identical retained asset: pointer identity
	 * proves that the worker updated the same service-owned progressive mesh,
	 * while the source identities protect edits/replacements.  The retained
	 * occurrence does not have to carry the current camera epoch: continuous
	 * input can advance that metadata faster than a suffix read, which is the
	 * very case this rebase handles.  Preserve a
	 * budget-coarsened active cut.  A stale result may expose only the draw cut
	 * which its provider actually reserved—not its richer resident high-water
	 * mark—and the current renderer ceiling remains authoritative.  The next
	 * bounded view pass handles every other admission decision.
	 */
	const bool incomingProgressiveMesh =
	    drained[i].providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    drained[i].resultKind == BOBOL_LOD_RESULT_MESH &&
	    drained[i].progressiveMesh &&
	    drained[i].progressiveMesh->isValid();
	const bool currentDemandWitness = residentCadBefore &&
	    residentCadBefore->viewRevision ==
		this->d->lodViewRevision.value() &&
	    residentCadBefore->policyRevision ==
		this->d->lodPolicyRevision.value();
	const bool sameCumulativeCadAsset =
	    incomingProgressiveMesh && residentCadBefore &&
	    residentCadBefore->progressiveMesh == drained[i].progressiveMesh &&
	    residentCadBefore->databaseRevision ==
		drained[i].request.databaseRevision.value() &&
	    residentCadBefore->sourceContentHash ==
		drained[i].request.sourceContentHash &&
	    (drained[i].request.sourceContentHash != 0 ||
	     residentCadBefore->sourceRevision ==
		drained[i].request.sourceRevision.value()) &&
	    residentCadBefore->drawMode == drained[i].request.drawMode;
	const bool extendsCumulativeCadAsset =
	    sameCumulativeCadAsset &&
	    (drained[i].residentCut > residentCadBefore->residentCut ||
	     (drained[i].preparedCadGeometryRevision &&
	      drained[i].preparedCadGeometryRevision >
		residentCadBefore->preparedCadGeometryRevision));
	/* A worker completion publishes residency, not a coarsening decision.
	 * The current aggregate allocator may already have selected a richer cut
	 * from the same cumulative asset while this suffix was in flight.  Keep
	 * that cut; camera/policy retargeting is the sole authority which lowers
	 * presentation.  A genuine suffix extension necessarily contains the
	 * older drawable prefix. */
	const bool resultCanPreserveRicherActiveCut =
	    extendsCumulativeCadAsset &&
	    drained[i].geometry.activeCut < residentCadBefore->activeCut &&
	    (residentCadBefore->presentedChunks.empty() ?
		drained[i].progressiveMesh->canDrawCut(
		    residentCadBefore->activeCut) :
		drained[i].progressiveMesh->canDrawChunksAtCut(
		    residentCadBefore->presentedChunks,
		    residentCadBefore->activeCut));
	/* A provider completion owns immutable residency, not the occurrence's
	 * draw budget.  A coalesced worker may return geometry beyond the cut
	 * selected while it was in flight.  The current allocation may admit some
	 * or all of that result, however: clamping to the allocated cut lets a
	 * prepared suffix become visible without either exceeding the scene plan
	 * or waiting for another identical provider request. */
	const int currentAllocatedCut = viewState && residentCadBefore ?
	    viewState->currentCadAllocatedCut(residentCadBefore,
		residentCadBefore->viewRevision,
		residentCadBefore->policyRevision,
		residentCadBefore->drawMode) : -1;
	const bool allocationCoversProviderAdvance = residentCadBefore &&
	    currentAllocatedCut > residentCadBefore->activeCut;
	const int admittedProviderActiveCut = allocationCoversProviderAdvance ?
	    std::min(drained[i].geometry.activeCut,
		currentAllocatedCut) :
	    (residentCadBefore ? residentCadBefore->activeCut : -1);
	const bool providerResultNeedsAllocationClamp =
	    sameCumulativeCadAsset &&
	    drained[i].geometry.activeCut > residentCadBefore->activeCut &&
	    drained[i].request.requiredChunks ==
		residentCadBefore->requiredChunks &&
	    (residentCadBefore->requiredChunks.empty() ?
		drained[i].progressiveMesh->canDrawCut(
		    residentCadBefore->activeCut) :
		drained[i].progressiveMesh->canDrawChunksAtCut(
		    residentCadBefore->requiredChunks,
		    residentCadBefore->activeCut));
	if (resultCanPreserveRicherActiveCut ||
	    providerResultNeedsAllocationClamp) {
	    drained[i].geometry.activeCut = resultCanPreserveRicherActiveCut ?
		residentCadBefore->activeCut : admittedProviderActiveCut;
	    if (residentCadBefore->requiredChunks.empty()) {
		drained[i].counts.faceCount =
		    drained[i].progressiveMesh->faceCount(
			drained[i].geometry.activeCut);
		drained[i].counts.pointCount =
		    drained[i].progressiveMesh->pointCount(
			drained[i].geometry.activeCut);
		drained[i].counts.originalPointCount =
		    drained[i].counts.pointCount;
		drained[i].counts.normalCount = drained[i].hasNormals ?
		    drained[i].counts.faceCount * 3 : 0;
	    } else {
		(void)drained[i].progressiveMesh->countsForChunksAtCut(
		    residentCadBefore->requiredChunks,
		    drained[i].geometry.activeCut, drained[i].hasNormals,
		    &drained[i].counts);
	    }
	    for (BObolLodPresentationLayer &layer :
		 drained[i].presentationLayers)
		layer.activeCut = drained[i].geometry.activeCut;
	}
	const bool staleViewOrPolicy =
	    drained[i].request.viewRevision !=
		this->d->lodViewRevision.value() ||
	    drained[i].request.policyRevision !=
		this->d->lodPolicyRevision.value();
	if (staleViewOrPolicy && incomingProgressiveMesh && residentCadBefore &&
	    getenv("BOBOL_LOD_TRACE_REJECTIONS"))
	    bu_log("BObol LoD cumulative rebase candidate object=%s "
		   "routing=%llu current_demand=%d "
		   "same_mesh=%d database=%llu/%llu source=%llu/%llu "
		   "content=%llu/%llu mode=%d/%d resident=%d/%d\n",
		   drained[i].request.objectName.getString(),
		   static_cast<unsigned long long>(
		       drained[i].request.sourceRoutingId.value()),
		   currentDemandWitness ? 1 : 0,
		   residentCadBefore->progressiveMesh ==
		       drained[i].progressiveMesh ? 1 : 0,
		   static_cast<unsigned long long>(
		       residentCadBefore->databaseRevision),
		   static_cast<unsigned long long>(
		       drained[i].request.databaseRevision.value()),
		   static_cast<unsigned long long>(
		       residentCadBefore->sourceRevision),
		   static_cast<unsigned long long>(
		       drained[i].request.sourceRevision.value()),
		   static_cast<unsigned long long>(
		       residentCadBefore->sourceContentHash),
		   static_cast<unsigned long long>(
		       drained[i].request.sourceContentHash),
		   residentCadBefore->drawMode, drained[i].request.drawMode,
		   residentCadBefore->residentCut, drained[i].residentCut);
	if (staleViewOrPolicy && extendsCumulativeCadAsset) {
	    int activeCut = residentCadBefore->activeCut;
	    const bool waitingForResidentSuffix =
		residentCadBefore->requestedCut >
		    residentCadBefore->residentCut &&
		residentCadBefore->activeCut >=
		    residentCadBefore->residentCut;
	    if (waitingForResidentSuffix) {
		int admittedCut = std::min(
		    residentCadBefore->requestedCut,
		    std::min(drained[i].residentCut,
			drained[i].geometry.activeCut));
		if (this->d->lodInteractiveProgressiveCeiling >= 0)
		    admittedCut = std::min(admittedCut,
			this->d->lodInteractiveProgressiveCeiling);
		activeCut = std::max(activeCut, admittedCut);
	    }
	    if (!(residentCadBefore->presentedChunks.empty() ?
		drained[i].progressiveMesh->canDrawCut(activeCut) :
		drained[i].progressiveMesh->canDrawChunksAtCut(
		    residentCadBefore->presentedChunks, activeCut)))
		activeCut = residentCadBefore->activeCut;

	    /* Rebase one immutable presentation transaction, not only its cut.
	     * The provider task's projection and page census belong to its old
	     * camera.  Publishing those fields under the current epoch made an
	     * equal cut appear satisfied while its HUD error and visible population
	     * still described the prior close view.  A resident compaction may also
	     * complete after the provider result was queued; only the exact prepared
	     * generation may witness the current page set. */
	    const bool preparedGenerationCurrent =
		drained[i].preparedCadGeometry &&
		drained[i].preparedCadGeometryRevision != 0 &&
		drained[i].preparedCadGeometryRevision ==
		    drained[i].progressiveMesh->revision();
	    const bool currentPagesDrawable =
		residentCadBefore->requiredChunks.empty() ?
		drained[i].progressiveMesh->canDrawCut(activeCut) :
		drained[i].progressiveMesh->canDrawChunksAtCut(
		    residentCadBefore->requiredChunks, activeCut);
	    BObolLodCounts currentCounts;
	    bool currentCountsValid = currentPagesDrawable;
	    if (currentCountsValid &&
		residentCadBefore->requiredChunks.empty()) {
		currentCounts.faceCount =
		    drained[i].progressiveMesh->faceCount(activeCut);
		currentCounts.pointCount =
		    drained[i].progressiveMesh->pointCount(activeCut);
		currentCounts.originalPointCount = currentCounts.pointCount;
	    } else if (currentCountsValid) {
		currentCountsValid =
		    drained[i].progressiveMesh->countsForChunksAtCut(
			residentCadBefore->requiredChunks, activeCut,
			drained[i].hasNormals, &currentCounts) ? true : false;
	    }
	    if (preparedGenerationCurrent && currentCountsValid) {
		drained[i].request.viewRevision =
		    this->d->lodViewRevision.value();
		drained[i].request.policyRevision =
		    this->d->lodPolicyRevision.value();
		drained[i].request.visualEmphasis =
		    residentCadBefore->visualEmphasis;
		drained[i].request.projectedPixelDiameter =
		    residentCadBefore->projectedPixelDiameter;
		drained[i].request.projectedPixelArea =
		    residentCadBefore->projectedPixelArea;
		drained[i].request.projectedPixelPerimeter =
		    residentCadBefore->projectedPixelPerimeter;
		drained[i].request.projectedBoundsContained =
		    residentCadBefore->projectedBoundsContained;
		drained[i].request.targetPixelError =
		    residentCadBefore->targetPixelError;
		drained[i].request.requestedCut =
		    residentCadBefore->requestedCut;
		drained[i].request.requiredChunks =
		    residentCadBefore->requiredChunks;
		/* The retained payload does not own the full camera matrices.  The
		 * exact current chunk list and counts above are sufficient for this
		 * publication; invalidate the old task's spatial projection so it
		 * cannot manufacture a stale projected per-cut census. */
		drained[i].request.spatialProjectionValid = FALSE;
		drained[i].resolvedCut = residentCadBefore->requestedCut;
		drained[i].cacheKey = bobol_lod_cache_key(drained[i].request);
		drained[i].geometry.activeCut = activeCut;
		drained[i].counts = currentCounts;
		drained[i].counts.normalCount = drained[i].hasNormals ?
		    drained[i].counts.faceCount * 3 : 0;
		drained[i].terminal =
		    drained[i].resolvedCut < 0 ||
		    activeCut >= drained[i].resolvedCut ? TRUE : FALSE;
		drained[i].stale = FALSE;
	    } else if (getenv("BOBOL_LOD_TRACE_REJECTIONS")) {
		bu_log("BObol LoD cumulative rebase deferred object=%s "
		       "prepared_current=%d pages_drawable=%d counts_valid=%d "
		       "prepared_revision=%llu mesh_revision=%llu\n",
		       drained[i].request.objectName.getString(),
		       preparedGenerationCurrent ? 1 : 0,
		       currentPagesDrawable ? 1 : 0,
		       currentCountsValid ? 1 : 0,
		       static_cast<unsigned long long>(
			   drained[i].preparedCadGeometryRevision),
		       static_cast<unsigned long long>(
			   drained[i].progressiveMesh->revision()));
	    }
	}
	const SbBool traceResult = controller_lod_trace_result(drained[i]);
	if (traceResult) {
	    const BObolViewLodState::CadPayload *resident =
		findResidentCad();
	    bu_log("BObol LoD apply trace object=%s request_cut=%d "
		   "loaded_cut=%d incoming_view=%llu incoming_policy=%llu "
		   "current_view=%llu current_policy=%llu resident_cut=%d "
		   "resident_view=%llu resident_policy=%llu\n",
		   drained[i].request.objectName.getString(),
		   drained[i].resolvedCut,
		   drained[i].geometry.activeCut,
		   static_cast<unsigned long long>(
		       drained[i].request.viewRevision.value()),
		   static_cast<unsigned long long>(
		       drained[i].request.policyRevision.value()),
		   static_cast<unsigned long long>(
		       this->d->lodViewRevision.value()),
		   static_cast<unsigned long long>(
		       this->d->lodPolicyRevision.value()),
		   resident ? resident->activeCut : -1,
		   static_cast<unsigned long long>(
		       resident ? resident->viewRevision : 0),
		   static_cast<unsigned long long>(
		       resident ? resident->policyRevision : 0));
	}
	if (drained[i].request.viewRevision !=
		this->d->lodViewRevision.value() ||
	    drained[i].request.policyRevision !=
		this->d->lodPolicyRevision.value()) {
	    /* A completed mesh from the prior camera remains a valid progressive
	     * bootstrap when this occurrence still has only its box or an older
	     * mesh epoch.  Reject it once equal/newer view data is resident so an
	     * out-of-order completion cannot replay a coarse cut. */
	    const BObolViewLodState::CadPayload *cadPayload =
		findResidentCad();
	    const BObolViewLodState::MeshPayload *meshPayload = viewState ?
		viewState->findMeshForResult(drained[i]) : NULL;
	    const bool cadMesh =
		cadPayload &&
		(cadPayload->resultKind == BOBOL_LOD_RESULT_MESH ||
		 cadPayload->resultKind == BOBOL_LOD_RESULT_FULL_DETAIL);
	    const bool shapeMesh =
		meshPayload &&
		meshPayload->resultKind == BOBOL_LOD_RESULT_MESH;
	    const auto epochNewer = [&drained, i](uint64_t policy,
		    uint64_t view) {
		return policy > drained[i].request.policyRevision ||
		    (policy == drained[i].request.policyRevision &&
		     view > drained[i].request.viewRevision);
	    };
	    const auto epochEqual = [&drained, i](uint64_t policy,
		    uint64_t view) {
		return policy == drained[i].request.policyRevision &&
		    view == drained[i].request.viewRevision;
	    };
	    const bool cadPopulationSupersedes = cadMesh &&
		cadPayload->activeCut >= drained[i].geometry.activeCut &&
		cadPayload->counts.faceCount >= drained[i].counts.faceCount &&
		cadPayload->counts.pointCount >= drained[i].counts.pointCount &&
		cadPayload->presentationLayers.size() >=
		    drained[i].presentationLayers.size();
	    const bool shapePopulationSupersedes = shapeMesh &&
		meshPayload->activeCut >= drained[i].geometry.activeCut &&
		meshPayload->counts.faceCount >= drained[i].counts.faceCount &&
		meshPayload->counts.pointCount >= drained[i].counts.pointCount;
	    const bool residentSupersedes =
		(cadMesh &&
		 (epochNewer(cadPayload->policyRevision,
		      cadPayload->viewRevision) ||
		  (epochEqual(cadPayload->policyRevision,
		       cadPayload->viewRevision) &&
		   cadPopulationSupersedes))) ||
		(shapeMesh &&
		 (epochNewer(meshPayload->policyRevision,
		      meshPayload->viewRevision) ||
		  (epochEqual(meshPayload->policyRevision,
		       meshPayload->viewRevision) &&
		   shapePopulationSupersedes)));
	    const bool incomingMesh =
		drained[i].resultKind == BOBOL_LOD_RESULT_MESH ||
		drained[i].resultKind == BOBOL_LOD_RESULT_FULL_DETAIL;
	    if (!incomingMesh || residentSupersedes) {
		if (traceResult)
		    bu_log("BObol LoD apply trace rejected_by_epoch object=%s\n",
			   drained[i].request.objectName.getString());
		if (getenv("BOBOL_LOD_TRACE_REJECTIONS")) {
		    static std::atomic<unsigned int> rejectionTraceCount(0);
		    const unsigned int traceIndex =
			rejectionTraceCount.fetch_add(1);
		    if (traceIndex < 64)
			bu_log("BObol LoD rejection reason=epoch object=%s "
			       "occurrence=%s incoming_cut=%d requested=%d "
			       "incoming_view=%llu incoming_policy=%llu "
			       "current_view=%llu current_policy=%llu "
			       "resident_cut=%d resident_view=%llu "
			       "resident_policy=%llu supersedes=%d\n",
			       drained[i].request.objectName.getString(),
			       drained[i].request.occurrenceKey.getString(),
			       drained[i].geometry.activeCut,
			       drained[i].resolvedCut,
			       static_cast<unsigned long long>(
				   drained[i].request.viewRevision.value()),
			       static_cast<unsigned long long>(
				   drained[i].request.policyRevision.value()),
			       static_cast<unsigned long long>(
				   this->d->lodViewRevision.value()),
			       static_cast<unsigned long long>(
				   this->d->lodPolicyRevision.value()),
			       cadPayload ? cadPayload->activeCut :
				   (meshPayload ?
				       meshPayload->activeCut : -1),
			       static_cast<unsigned long long>(
				   cadPayload ? cadPayload->viewRevision :
				       (meshPayload ?
					   meshPayload->viewRevision : 0)),
			       static_cast<unsigned long long>(
				   cadPayload ? cadPayload->policyRevision :
				       (meshPayload ?
					   meshPayload->policyRevision : 0)),
			       residentSupersedes ? 1 : 0);
		}
		this->d->lastLodRejectedResultCount++;
		append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
						 drained[i].request.objectPath,
						 "superseded LoD result epoch rejected");
		continue;
	    }
	}
	const int priorActiveCut = residentCadBefore ?
	    residentCadBefore->activeCut : drained[i].geometry.activeCut;
	const int availableActiveCut = std::max(
	    priorActiveCut, drained[i].geometry.activeCut);
	/* Residency and presentation are separate contracts.  A worker may append
	 * the complete requested suffix while preserving the allocator's old
	 * active cut.  Remember that newly usable range only after this result is
	 * accepted; the coalesced scene allocator consumes it below. */
	const SbBool resultOffersRicherResidentPrefix =
	    incomingProgressiveMesh &&
	    drained[i].residentCut > availableActiveCut &&
	    (drained[i].resolvedCut > availableActiveCut ||
	     (residentCadBefore &&
	      residentCadBefore->requestedCut > availableActiveCut)) ?
		TRUE : FALSE;
	const SbBool partialRefinementResult =
	    drained[i].providerStatus == BOBOL_LOD_PROVIDER_READY &&
	    drained[i].resultKind == BOBOL_LOD_RESULT_MESH &&
	    drained[i].geometry.activeCut >= 0 &&
	    drained[i].resolvedCut >
		drained[i].geometry.activeCut &&
	    !drained[i].terminal &&
	    (drained[i].request.sourceRoutingId == 0 ||
	     !residentCadBefore ||
	     residentCadBefore->activeCut <
		 drained[i].geometry.activeCut ||
	     (residentMeshBefore &&
	      residentMeshBefore->activeCut <
		  drained[i].geometry.activeCut)) ?
	    TRUE : FALSE;
	const bool retainPresentationForResidentGrowth =
	    BObolLodAvailabilityScheduler::canRetainPresentation(
		this->d->lodAvailabilityLedger.residencyDrainActive(),
		residentCadBefore != NULL,
		incomingProgressiveMesh && residentCadBefore &&
		    residentCadBefore->progressiveMesh ==
			drained[i].progressiveMesh,
		residentCadBefore &&
		    residentCadBefore->activeCut ==
			drained[i].geometry.activeCut,
		resultOffersRicherResidentPrefix != FALSE);
	if (drained[i].request.sourceRoutingId != 0 &&
	    drained[i].request.occurrenceKey.getLength() > 0) {
	    if (route && route->hasCompactInstanceKey(
		    drained[i].request.occurrenceKey.getString())) {
		directMatched++;
		if (viewState &&
		    viewState->consumeSourceResult(
			route, drained[i])) {
		    directApplied++;
		    if (retainPresentationForResidentGrowth)
			directRetainedPresentationCount++;
		    if (resultOffersRicherResidentPrefix)
			residentGrowthAdmissionCandidate = TRUE;
		    /*
		     * The first useful compact result is normally the minimum
		     * PoP prefix.  Absence of an older payload does not make that
		     * result terminal: publishing it must install the same
		     * render/refinement barrier as every later prefix.  Otherwise
		     * a cold single-mesh scene can go idle at its coarsest cut and
		     * refine only after an unrelated camera or paint event.
		     */
		    if (partialRefinementResult)
			partialRefinementCandidate = TRUE;
		} else {
		    directRejected++;
		    if (getenv("BOBOL_LOD_TRACE_REJECTIONS")) {
			static std::atomic<unsigned int>
			    sourceRejectionTraceCount(0);
			const unsigned int traceIndex =
			    sourceRejectionTraceCount.fetch_add(1);
			if (traceIndex < 64)
			    bu_log("BObol LoD rejection reason=source "
				   "object=%s occurrence=%s cut=%d "
				   "requested=%d progressive=%d valid=%d "
				   "route=%p diagnostic=%s\n",
				   drained[i].request.objectName.getString(),
				   drained[i].request.occurrenceKey.getString(),
				   drained[i].geometry.activeCut,
				   drained[i].resolvedCut,
				   drained[i].progressiveMesh ? 1 : 0,
				   drained[i].progressiveMesh &&
					   drained[i].progressiveMesh->isValid() ?
				       1 : 0,
				   static_cast<void *>(route),
				   drained[i].diagnostic.getString());
		    }
		    append_controller_lod_diagnostic(
			this->d->lastLodDiagnostics,
			drained[i].request.objectPath,
			viewState ?
			    "view-local compact CAD LoD result rejected by source" :
			    "view-local compact CAD LoD update requires a view state");
		}
		continue;
	    }
	    /* A source replacement may legitimately invalidate the token while a
	     * cancellation races a completed result.  Let the generic action
	     * diagnose or match that exceptional stale result by key. */
	}
	if (partialRefinementResult)
	    partialRefinementCandidate = TRUE;
	if (resultOffersRicherResidentPrefix)
	    residentGrowthAdmissionCandidate = TRUE;
	update.addResult(std::move(drained[i]));
    }

    unsigned int rootlessFallbackCount = 0;
    if (update.getResultCount() > 0) {
	SoNode *root = this->getRenderSceneRoot();
	if (!root)
	    root = this->getSceneRoot();
	if (root) {
	    update.apply(root);
	} else {
	    rootlessFallbackCount = static_cast<unsigned int>(
		std::min<size_t>(update.getResultCount(),
		    std::numeric_limits<unsigned int>::max()));
	    append_controller_lod_diagnostic(
		this->d->lastLodDiagnostics, SbString(""),
		"direct-node LoD result application requires a scene root");
	}
    }
    this->d->lastLodMatchedResultCount =
	directMatched + update.getMatchedResultCount();
    this->d->lastLodAppliedResultCount =
	directApplied + update.getAppliedResultCount();
    this->d->lastLodRejectedResultCount +=
	directRejected + update.getRejectedResultCount();
    this->d->lastLodUnmatchedResultCount =
	rootlessFallbackCount + update.getUnmatchedResultCount();

    /* A bounded quiet residency drain deliberately publishes only immutable
     * array growth while retaining every occurrence's current active cut.
     * Do not turn the resulting quality debt into a per-batch presentation
     * barrier.  ResidentGrowthPolicy owns the liveness witness and schedules
     * one aggregate allocation after the complete worker/result wave. */
    const bool retainedPresentationBatch =
	this->d->lastLodAppliedResultCount > 0 &&
	directApplied == directRetainedPresentationCount &&
	update.getAppliedResultCount() == 0;

    /* Provider terminality describes the task snapshot, not necessarily the
     * owner-thread demand which survived its publication.  A coalesced wheel
     * stream or exact-view history recall can leave the newly published
     * occurrence deliberately unsatisfied.  Arm one bounded replay from the
     * authoritative view-state counters so an older terminal result cannot
     * make that current quality debt go idle.  The ordinary allocator remains
     * responsible for deciding whether the replay is drawable,
     * performance-limited, or memory-limited. */
    if (viewState && this->d->lastLodResultCount > 0) {
	size_t activePayloads = 0;
	size_t satisfiedPayloads = 0;
	size_t memoryLimitedPayloads = 0;
	viewState->convergencePayloadCounts(
	    activePayloads, satisfiedPayloads, memoryLimitedPayloads);
	if (!retainedPresentationBatch &&
	    activePayloads > satisfiedPayloads &&
	    activePayloads - satisfiedPayloads > memoryLimitedPayloads)
	    partialRefinementCandidate = TRUE;
    }
    if (update.getResultCount() > 0 &&
	update.getDiagnostics().getLength() > 0) {
	if (this->d->lastLodDiagnostics.getLength() > 0)
	    this->d->lastLodDiagnostics += "\n";
	this->d->lastLodDiagnostics += update.getDiagnostics();
    }
    this->d->lodAvailabilityLedger.setResultsPending(
	service->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) > 0);

    if (this->d->lastLodAppliedResultCount > 0) {
	const size_t applied =
	    static_cast<size_t>(this->d->lastLodAppliedResultCount);
	/* A coalesced resident-growth task may finish after the camera/policy
	 * pass which tried to reserve its next draw cut.  The stale result above
	 * deliberately publishes only its provider-admitted presentation level;
	 * once its richer immutable prefix is installed, invalidate the completed
	 * demand pass so the current aggregate allocator can retarget it without
	 * another cache read, camera event, or quiet transition. */
	if (residentGrowthAdmissionCandidate && this->d->lodAutoSubmit) {
	    this->d->lodAvailabilityLedger.noteRicherPrefixAvailable();
	    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		    this->d->lodViewRevision.value()))
		bu_log("BObol LoD resident growth noted applied=%zu "
		       "interactive=%d scale_changed=%d\n",
		       applied, this->d->lodInteractionSession.active() ? 1 : 0,
		       this->d->lodViewDemandPolicy.interactionScaleChanged() ?
			   1 : 0);
	    /* The missing suffix may complete after this scale epoch has already
	     * spent its one quality probe.  Residency alone is not a visible
	     * improvement, so rearm exactly one bounded probe for the newly
	     * available population.  An ordinary interactive/first-useful result
	     * publishes below; a quiet residency-only wave retains the coherent
	     * frame until its aggregate reallocation consumes this latch. */
	    (void)this->d->lodViewDemandPolicy.rearmAfterResidentGrowth(
		this->d->lodInteractionSession.active() != FALSE);
	    this->markProgressiveWorkPending();
	}
	const int64_t now = bu_gettime();
	/* Residency-only batches do not change the selected presentation.  Their
	 * resident-growth edge above keeps the pump alive and eventually requests
	 * one scene-wide allocation/frame, so recording them as unpresented pixels
	 * would only serialize the worker stream behind repeated whole-scene
	 * renders. */
	if (!retainedPresentationBatch) {
	    /* The queue may already have waited for a complete drain wave.  Preserve
	     * that first-ready timestamp as the publication batch's age origin so a
	     * 250 ms drain bound is not followed by a second 250 ms frame bound. */
	    const int64_t firstReady =
		this->d->lodAvailabilityLedger.firstResultReadyMicroseconds();
	    this->d->lodPresentationTransaction.noteApplied(
		applied, firstReady > 0 ? firstReady : now,
		this->d->lodViewRevision.value(),
		this->d->lodPolicyRevision.value());
	}
	const bool firstUsefulMesh = !retainedPresentationBatch &&
	    this->getActiveLodMeshPayloadCount() <= applied;
	const bool serviceProducer =
	    service->queuedResultCountForGeneration(
		this->d->lodActiveGeneration) > 0 ||
	    service->activeTaskCountForGeneration(
		this->d->lodActiveGeneration) > 0;
    const bool submissionPausedByPresentation =
	BObolLodAdmissionPlanner::presentationPausesSubmission(
	    this->d->lodPointAdmissionFrame.pending(),
	    this->d->lodPointQualityPhase.presentationPending(),
	    this->d->lodAdmissionEvidence.capacity().rescanAfterFrame(),
	    viewState && viewState->hasCadPresentationAssemblies());
	const bool streamIdle =
	    !BObolLodProducerPolicy::canProduceGeometry(
		this->d->lodSubmissionPass.active() != FALSE,
		submissionPausedByPresentation,
		this->d->lodAvailabilityLedger.providerPendingCount() > 0,
		serviceProducer);
	BObolLodPresentationTransaction::Inputs publicationInputs;
	publicationInputs.nowMicroseconds = now;
	publicationInputs.observedRenderNanoseconds = std::max(
	    this->d->lastSceneRenderTimeNanoseconds,
	    this->d->lastRenderTimeNanoseconds);
	publicationInputs.interactive = this->d->lodInteractionSession.active() != FALSE;
	publicationInputs.firstUseful = firstUsefulMesh;
	publicationInputs.streamIdle = streamIdle;
	const BObolLodPresentationTransaction::Decision publicationDecision =
	    this->d->lodPresentationTransaction.decide(publicationInputs);
	const bool publishNow =
	    this->d->lodPresentationTransaction.publicationFramePending();
	if (partialRefinementCandidate &&
	    !this->d->lodPresentationTransaction.barrierPending()) {
	    /*
	     * Hold the next prefix immediately, but allow independent coverage
	     * results to accumulate until the adaptive publication deadline.
	     * scheduleLodRefinementFrame() below supplies the host wakeup once
	     * that deadline is reached.
	     */
	    (void)this->d->lodPresentationTransaction.arm(
		BObolLodPresentationTransaction::REASON_RESULT_PUBLICATION,
		this->d->renderCompletionSerial,
		this->d->lodViewRevision.value(),
		this->d->lodPolicyRevision.value());
	}
	if (publicationDecision.requestFrame) {
	    this->requestRender("lod-result");
	}
	this->d->lodCompactionPolicy.requestAfter(bu_gettime(), 750000);
	(void)this->enforceMeshResidencyBudget();
	if (partialRefinementCandidate && publishNow)
	    this->scheduleLodRefinementFrame("lod-result");
    } else if (partialRefinementCandidate && this->d->lodAutoSubmit) {
	/* A cumulative result may extend the shared resident asset yet be
	 * correctly rejected because its immutable presentation cannot cover
	 * the current page set.  Draining it releases the service's coalescing
	 * key, so immediately restart one current-demand pass.  Waiting for an
	 * applied presentation here leaves no worker, frame barrier, or input
	 * edge capable of resolving the authoritative quality debt. */
	this->d->lodSubmissionSourceIndex = 0;
	this->d->lodSubmissionEntryOffset = 0;
	this->d->clearLodSubmissionPlan();
	this->d->lodSubmissionPass.clearRescan();
	this->d->lodSubmissionPass.activate();
	this->d->lodRetainedPass.noteRefinementPending();
	this->d->lodRetainedPass.clearAdmittedWork();
	this->markProgressiveWorkPending();
	this->requestRender("lod-current-demand-replay");
    }

    return size_to_int_saturated(
	       static_cast<size_t>(this->d->lastLodAppliedResultCount));
}

uint64_t
BObolViewController::getLodViewRevision(void) const
{
    return this->d->lodViewRevision.value();
}

void
BObolViewController::setLodPolicyRevision(uint64_t revision)
{
    uint64_t newRevision = revision == 0 ? 1 : revision;
    if (this->d->lodPolicyRevision.value() == newRevision)
	return;
    this->d->resetLodViewQualityHistory();
    this->d->setAdmissionPolicyRevision(newRevision);
    this->d->lodStaticQualityTrial.reset();
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    /*
     * A quality-policy revision changes the requested PoP cut, not camera
     * visibility.  Retain the last complete current-view denominator until a
     * camera or source-inventory revision invalidates it.  Clearing it here
     * let the quiet transition report a 50k scene converged with a zero
     * visible target while its unsatisfied-frontier fast path deliberately
     * skipped another full visibility scan.
     */
    this->d->resetLodConvergenceFraction();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_RECOVERY_CEILING);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_DEADLINE_CAPACITY_CEILING);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR);
    this->d->lodPointQualityPhase.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    this->requestRender("lod-policy");
}

uint64_t
BObolViewController::getLodPolicyRevision(void) const
{
    return this->d->lodPolicyRevision.value();
}

unsigned int
BObolViewController::getLastLodVisitedMeshCount(void) const
{
    return this->d->lastLodVisitedMeshCount;
}

unsigned int
BObolViewController::getLastLodSubmittedTaskCount(void) const
{
    return this->d->lastLodSubmittedTaskCount;
}

unsigned int
BObolViewController::getLastLodUpdatedCutCount(void) const
{
    return this->d->lastLodUpdatedCutCount;
}

unsigned int
BObolViewController::getLastLodSkippedMeshCount(void) const
{
    return this->d->lastLodSkippedMeshCount;
}

size_t
BObolViewController::getLastLodResultCount(void) const
{
    return this->d->lastLodResultCount;
}

unsigned int
BObolViewController::getLastLodMatchedResultCount(void) const
{
    return this->d->lastLodMatchedResultCount;
}

unsigned int
BObolViewController::getLastLodAppliedResultCount(void) const
{
    return this->d->lastLodAppliedResultCount;
}

unsigned int
BObolViewController::getLastLodRejectedResultCount(void) const
{
    return this->d->lastLodRejectedResultCount;
}

unsigned int
BObolViewController::getLastLodUnmatchedResultCount(void) const
{
    return this->d->lastLodUnmatchedResultCount;
}

const SbString &
BObolViewController::getLastLodDiagnostics(void) const
{
    return this->d->lastLodDiagnostics;
}

size_t
BObolViewController::getActiveLodMeshPayloadCount(void) const
{
    const BObolViewLodState *state =
	this->d->viewAttachment->getViewLodState();
    return state ? state->meshPayloadCount() + state->cadMeshPayloadCount() : 0;
}

size_t
BObolViewController::getActiveLodProxyPayloadCount(int proxyKind) const
{
    const BObolViewLodState *state =
	this->d->viewAttachment->getViewLodState();
    return state ? state->proxyPayloadCount(proxyKind) +
	   state->cadProxyPayloadCount(proxyKind) : 0;
}

size_t
BObolViewController::getActiveLodCadPayloadCount(void) const
{
    return this->d->viewAttachment->getViewLodState() ? this->d->viewAttachment->getViewLodState()->cadPayloadCount() : 0;
}

void
BObolViewController::syncRenderManager(void)
{
    this->d->renderManager->setSceneGraph(this->getRenderRoot());
    this->d->renderManager->setCamera(this->d->activeCamera);
    this->d->renderManager->setViewportRegion(this->d->viewportRegion);
}

void
BObolViewController::advanceLodViewRevision(void)
{
    this->d->advanceAdmissionRevision(
	BObolLodAdmissionRevisionDomain::VIEW);
    this->d->lodInterruptedPresentationReplay.retire();
    this->d->resetRetainedAdmissionQualityProof();
    this->d->lodStaticQualityTrial.reset();
    BObolViewLodState *presentationState =
	this->d->viewAttachment->getViewLodState();
    if (presentationState) {
	presentationState->clearCadOccurrenceTerminalFailures();
	presentationState->setCadPresentationProgressiveCutCeiling(
	    this->d->lodInteractiveProgressiveCeiling);
    }
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    this->d->clearLodConvergenceCandidates();
    this->d->resetLodConvergenceFraction();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_RECOVERY_CEILING);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_DEADLINE_CAPACITY_CEILING);
    this->d->lodPointQualityPhase.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    this->d->lodPresentationPolicy.viewInvalidated();
    this->d->lodPoseContinuity.clearRetainOccurrenceCuts();
    this->d->lodPlanningObligations.retireImportanceCensus();
    this->d->resetDeadlineSafePresentation();
    this->d->lodCoveragePolicy.activate(true);
    this->d->lodViewDemandPolicy.refreshForViewRevision(
	this->d->lodInteractionSession.active() != FALSE);
}

void
BObolViewController::advanceLodPolicyRevision(
    LodPolicyTransition transition)
{
    this->d->advanceAdmissionRevision(
	BObolLodAdmissionRevisionDomain::POLICY);
    BObolViewLodState *presentationState =
	this->d->viewAttachment->getViewLodState();
    if (presentationState)
	presentationState->clearCadOccurrenceTerminalFailures();
    this->d->lodInterruptedPresentationReplay.retire();
    this->d->resetRetainedAdmissionQualityProof();
    /* External policy changes start from the preferred quiet cadence.  An
     * internal static-quality successor changes the requested pixel error
     * only after a completed frame has proved the longer deadline; dropping
     * that proof here made the successor immediately coarsen itself. */
    if (transition != LodPolicyTransition::CONTINUE_STATIC_QUALITY)
	this->d->lodStaticQualityTrial.deactivate();
    struct bv_lod_policy lodPolicy;
    bv_lod_policy_init(&lodPolicy);
    this->d->viewAttachment->getLodPolicy(&lodPolicy);
    this->d->lodCoveragePolicy.setRequired(
	this->d->lodAutoSubmit && this->d->lodService &&
	lodPolicy.policy != BV_LOD_OFF && lodPolicy.mesh_enabled);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    /* A policy epoch may change the target cadence (most importantly the
	 * 60 Hz motion -> 20 Hz quiet transition), refinement authority, and
	 * presentation handoff.  A partially consumed pass is meaningful only for
	 * the inputs which initialized it.  Carrying its retained-admission flag
	 * into the new epoch made a 0.03% budget rounding difference normalize all
	 * 2,500 Hubble occurrences to their minimum prefixes after rotation.  Keep
	 * the calibrated total budget, but recompute the pass decision from the new
	 * epoch's exact inputs.  The revision-stamped cursor rejects the old plan
	 * automatically. */
    /* Quality changes preserve the current view's proven visibility
     * denominator.  Source and camera revisions clear it explicitly. */
    this->d->resetLodConvergenceFraction();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_RETAINED_QUALITY_FLOOR);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CLEAR_DEADLINE_CAPACITY_CEILING);
    this->d->lodPointQualityPhase.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    this->d->lodCoveragePolicy.activate(false);
    this->d->lodViewDemandPolicy.refreshForPolicyRevision(
	transition == LodPolicyTransition::PRESERVE_SCALE_DEMAND,
	this->d->lodInteractionSession.active() != FALSE);
}

void
BObolViewController::beginLodInteraction(void)
{
    if (!this->d->lodAutoSubmit || this->d->lodInteractionSession.gestureActive())
	return;

    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    this->d->lodStaticQualityTrial.deactivate();
    const float previousPixelError = this->d->lodTargetPixelError;
    const bool newInteractionEpoch = !this->d->lodInteractionSession.active();
    int initialQualityFloor = -1;
    if (newInteractionEpoch) {
	BObolLodConvergenceStatus priorStatus;
	this->getLodConvergenceStatus(priorStatus);
	this->d->lodInteractionStartCertificate.capture(
	    this->d->lodViewScaleSignature,
	    this->d->lodViewSignature.has_value(),
	    priorStatus.viewReady,
	    this->d->lodAdmissionEvidence.capacity().currentBudget());
	const BObolViewLodState *snapshotState =
	    this->d->viewAttachment->getViewLodState();
	size_t presentedPrimitives = 0;
	if (priorStatus.viewReady && snapshotState &&
	    this->d->lastRenderTimeNanoseconds > 0 &&
	    this->d->lastRenderTimeNanoseconds <=
		BObolLodViewDemandPolicy::qualityFrameDurationNanoseconds() &&
	    snapshotState->lastCadPresentedPrimitiveCount(
		presentedPrimitives) && presentedPrimitives > 0) {
	    initialQualityFloor =
		snapshotState->maximumActiveProgressiveCut();
	    if (this->d->lodInteractiveProgressiveCeiling >= 0 &&
		initialQualityFloor >= 0)
		initialQualityFloor = std::min(initialQualityFloor,
		    this->d->lodInteractiveProgressiveCeiling);
	}
	this->d->lodPresentationPolicy.capturePrior(
	    this->d->lodTargetPixelError,
	    this->d->lodInteractiveProgressiveCeiling,
	    snapshotState ? snapshotState->
		cadPresentationProgressiveCutNextFraction() : 0.0f,
	    this->d->pointProxyAggregationApplicable() ?
		this->d->lodPresentationPointProxyPixelThreshold : 1.0f,
	    controller_lod_presentation_population(snapshotState,
		this->d->lodViewQualityDomainRevision),
	    this->d->lodViewRevision);
	this->d->seedInteractiveCalibrationFromStable();
    }
    this->d->lodPoseContinuity.clearRetainOccurrenceCuts();
    this->d->lodPlanningObligations.retireImportanceCensus();
    this->d->lodInteractionSession.beginGesture(bu_gettime());
    this->d->lodViewDemandPolicy.beginGesture(newInteractionEpoch);
    if (newInteractionEpoch)
	this->d->lodViewDemandPolicy.seedQualityFloor(initialQualityFloor);
    this->d->lodDiscretePopulationTrialPermit.revoke();
    this->d->lodPresentationPolicy.cancelHandoff();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
    this->d->lodPointQualityPhase.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->advanceAdmissionRevision(
        BObolLodAdmissionRevisionDomain::CAPACITY);
    {
	std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
	/* Rebase the human-facing number on this continuous-demand burst.
	 * Otherwise a long progressive/event-driven interval immediately before
	 * button-down can make a genuinely fluid drag report single-digit FPS
	 * for most of the gesture.  The prior gesture's smoothed interval is
	 * equally stale feedback for the new cut: renderer capacity is retained
	 * separately, while this burst must establish its own delivery cadence. */
	this->d->displayedPresentationIntervalNanoseconds = 0;
	this->d->smoothedInteractivePresentationIntervalNanoseconds = 0;
    }
    /* System GL traversal can return well before the timer-query sample for
     * the submitted CAD batch.  Motion policy must honor the slower of the
     * endpoint traversal and GPU presentation; using CPU time alone left a
     * 30-40 ms retained cut uncapped while targeting 60 FPS.  OSMesa reports
     * its work in the traversal time, so the same maximum is portable. */
    const uint64_t interactionTimingSample = std::max(
	this->d->smoothedRenderTimeNanoseconds,
	this->d->lodLastCadGpuTimeNanoseconds);
    this->d->lodTargetPixelError =
	BObolLodQualityPolicy::interactivePixelError(
	    interactionTimingSample,
	    this->d->lodInteractiveTargetFps);
    if (this->d->pointProxyAggregationApplicableForCameraTransition()) {
	this->d->lodPresentationPointProxyPixelThreshold =
	    BObolLodQualityPolicy::pointProxyThreshold(
		this->d->lodPresentationPointProxyPixelThreshold,
		interactionTimingSample,
		this->d->lodInteractiveTargetFps);
    } else {
	this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    }
    BObolViewLodState *viewState =
	this->d->viewAttachment->getViewLodState();
    if (viewState)
	viewState->setCadPresentationCameraMotionFrameReuse(FALSE);

    /*
     * Wheel and trackpad input is normally an unbracketed camera epoch.  It
     * may already have ratcheted the retained presentation to a measured
     * responsive ceiling before the user presses a mouse button.  Beginning
     * the drag is a continuation of that active interaction, not permission
     * to probe again from the scene maximum.  Preserve the current ceiling
     * as the basis for this correction; the quiet transition below is the
     * only path which may remove it and begin stable refinement.
     */

    const int activeMaximum = viewState ?
	viewState->maximumActiveProgressiveCut() : -1;
    const int responsiveCeiling = BObolLodQualityPolicy::progressiveCeiling(
	    activeMaximum, this->d->lodTargetPixelError,
	    this->d->lodInteractiveProgressiveCeiling);
    /* A fast pose-only gesture needs no ceiling at all: the existing
     * retained cut is already the correct presentation.  Installing a
     * numerically equivalent "maximum active level" ceiling here obscured
     * that contract and could become destructive as a richer suffix arrived
     * during rotation.  Scale-quality floors are introduced only after an
     * actual scale change is observed. */
    if (responsiveCeiling >= 0) {
	if (this->d->lodInteractiveProgressiveCeiling >= 0) {
	    /* A bracketed drag may continue an unbracketed scale epoch.  Keep
	     * that already measured limit and permit only a stricter correction. */
	    if (responsiveCeiling <
		    this->d->lodInteractiveProgressiveCeiling)
		this->d->lodInteractiveProgressiveCeiling =
		    responsiveCeiling;
	} else if (this->d->lodTargetPixelError > 1.01f) {
	    /* A new fast pose gesture whose retained cut already meets the
	     * one-pixel contract needs no numerical ceiling equal to that cut. */
	    this->d->lodInteractiveProgressiveCeiling =
		responsiveCeiling;
	}
    }
    this->d->lodInteractiveCeilingFeedbackRenderSerial =
	this->d->renderCompletionSerial;
    if (viewState)
	viewState->setCadPresentationProgressiveCutCeiling(
	    this->d->lodInteractiveProgressiveCeiling);
    if (viewState)
	viewState->setCadPresentationPointProxyPixelThreshold(
	    this->d->lodPresentationPointProxyPixelThreshold);
    /* The first pointer motion can precede the first camera signature change.
     * Treat the transition to the motion cut as a policy change so the
     * already resident mesh is coarsened before an expensive drag frame. */
    if (fabsf(this->d->lodTargetPixelError - previousPixelError) >
	    std::numeric_limits<float>::epsilon())
	this->advanceLodPolicyRevision();
    this->markProgressiveWorkPending();
    this->requestRender("lod-interaction-begin");
}

void
BObolViewController::endLodInteraction(void)
{
    if (!this->d->lodInteractionSession.gestureActive())
	return;

    (void)this->d->lodInteractionSession.endGesture(bu_gettime());
    this->d->lodViewDemandPolicy.endGesture();
    this->d->lodDiscretePopulationTrialPermit.revoke();
    {
	std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
	this->d->lastInteractivePresentationTimestampNanoseconds = 0;
    }
    /* Release keeps the last responsive motion presentation intact through
     * the debounce.  The quiet-successor reducer is the sole authority which
     * may restore a prior pose and it requires one current-pose presentation
     * proof before occurrence-local allocation.  Restoring here created a
     * second authority: a pre-quiet deadline miss could turn the temporary
     * motion budget into a persistent capacity ceiling. */
    /* A gesture pass may have exhausted the old, smaller allowance.  Force
     * the release pass to derive its admission state from the coherent floor
     * above rather than reusing that stale remainder. */
    this->d->advanceAdmissionRevision(
        BObolLodAdmissionRevisionDomain::CAPACITY);
    /* Keep the interaction epoch through the normal quiet-view debounce.
     * Pose-only presentation may already be restored, but new projection and
     * refinement work still waits so release cannot block on a full software
     * planning frame. */
    this->markProgressiveWorkPending();
    this->requestRender("lod-interaction-end");
}

SbBool
BObolViewController::isLodInteractionActive(void) const
{
    return this->d->lodInteractionSession.active();
}

SbBool
BObolViewController::isLodGestureActive(void) const
{
    return this->d->lodInteractionSession.gestureActive();
}

SbBool
BObolViewController::isLodScaleChangingInteraction(void) const
{
    return this->d->lodViewDemandPolicy.scaleChangingInteraction(
	this->d->lodInteractionSession.active() != FALSE) ? TRUE : FALSE;
}

float
BObolViewController::getLodTargetPixelError(void) const
{
    /* The controller stores the policy's base target.  Projection applies
     * view.lod.scale to that target before selecting a PoP cut, so callers
     * reporting the active quality contract must observe the same physical
     * error rather than the unscaled policy input. */
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    if (this->getViewInfo(&view) && std::isfinite(view.lod.scale) &&
	view.lod.scale > 0.0)
	return this->d->lodTargetPixelError /
	    static_cast<float>(view.lod.scale);
    return this->d->lodTargetPixelError;
}

int
BObolViewController::getLodInteractiveProgressiveCeiling(void) const
{
    return this->d->lodInteractiveProgressiveCeiling;
}

void
BObolViewController::setLodFrameRateTargets(float interactiveFps,
	float stableFps)
{
    if (!std::isfinite(interactiveFps) || interactiveFps <= 0.0f ||
	!std::isfinite(stableFps) || stableFps <= 0.0f)
	return;
    interactiveFps = BObolLodCoordinator::boundedPresentationTargetFps(
	interactiveFps,
	BObolLodCoordinator::maximumInteractivePresentationDeadline());
    stableFps = BObolLodCoordinator::boundedPresentationTargetFps(
	stableFps,
	BObolLodCoordinator::maximumStablePresentationDeadline());
    if (fabsf(this->d->lodInteractiveTargetFps - interactiveFps) <=
	    std::numeric_limits<float>::epsilon() &&
	fabsf(this->d->lodStableTargetFps - stableFps) <=
	    std::numeric_limits<float>::epsilon())
	return;
    this->d->resetLodViewQualityHistory();
    this->d->lodInteractiveTargetFps = interactiveFps;
    this->d->lodStableTargetFps = stableFps;
    /* Admission and traversal must honor the same frame contract.  Keeping
     * the endpoint at its default hard deadline after a user selected a
     * slower fidelity target made the allocator admit the requested detail
     * and the renderer reject it forever.  Preserve the responsive defaults
     * for ordinary targets, and bound how far either user-facing target may
     * relax an endpoint traversal.  Targets below those latency floors are
     * clamped to the effective minimum FPS, so admission cannot again outrun
     * presentation.  An embedding application may still call
     * setPresentationFrameDeadlines() afterward when it needs a stricter
     * contract. */
    this->d->interactivePresentationFrameDeadlineNanoseconds =
	BObolLodCoordinator::targetPresentationDeadline(
	    interactiveFps,
	    BObolLodCoordinator::defaultInteractivePresentationDeadline(),
	    BObolLodCoordinator::maximumInteractivePresentationDeadline());
    this->d->stablePresentationFrameDeadlineNanoseconds =
	BObolLodCoordinator::targetPresentationDeadline(
	    stableFps,
	    BObolLodCoordinator::defaultStablePresentationDeadline(),
	    BObolLodCoordinator::maximumStablePresentationDeadline());
    this->d->prominentQualityFrameDeadlineNanoseconds = std::max(
	this->d->stablePresentationFrameDeadlineNanoseconds,
	BObolLodCoordinator::defaultProminentQualityDeadline());
    this->d->consecutiveInterruptedPresentationFrames = 0;
    this->d->lodStaticQualityTrial.reset();
    this->advanceLodPolicyRevision();
    this->markProgressiveWorkPending();
    this->requestRender("lod-frame-rate");
}

float
BObolViewController::getLodInteractiveTargetFps(void) const
{
    return this->d->lodInteractiveTargetFps;
}

float
BObolViewController::getLodStableTargetFps(void) const
{
    return this->d->lodStableTargetFps;
}

size_t
BObolViewController::getCurrentLodRenderCostBudget(void) const
{
    return this->d->lodAdmissionEvidence.capacity().currentBudget();
}

size_t
BObolViewController::getActiveLodFaceCount(void) const
{
    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    return state ? state->activeFaceCount() : 0;
}

size_t
BObolViewController::getActiveLodRenderCost(void) const
{
    const BObolViewLodState *state = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    return state ? state->activeRenderCost() : 0;
}

void
BObolViewController::getLodConvergenceStatus(
	BObolLodConvergenceStatus &status) const
{
    status.clear();
    status.viewQualityHistoryEntryCount =
	this->d->lodViewQualityHistory.size();
    status.viewQualityHistoryRememberCount =
	this->d->lodViewQualityHistory.rememberCount();
    status.viewQualityHistoryRecallCount =
	this->d->lodViewQualityHistory.recallCount();
    const BObolLodAdmissionRevisionStamp revisionStamp =
	this->d->admissionRevisionStamp();
    status.inventoryRevision = revisionStamp.inventory.value();
    status.availabilityRevision = revisionStamp.availability.value();
    status.viewRevision = revisionStamp.view.value();
    status.policyRevision = revisionStamp.policy.value();
    status.capacityRevision = revisionStamp.capacity.value();
    status.presentationTransactionSerial =
	this->d->lodPresentationTransaction.sequence();
    status.presentationRequiredRenderSerial =
	this->d->lodPresentationTransaction.requiredRenderSerial();
    {
	std::lock_guard<std::mutex> lock(this->d->presentationTimingMutex);
	status.presentedFrameSerial = this->d->presentedFrameSerial;
    }
    status.failedSourceCount = this->getLastFailedSourceCount();
    status.visibleTargetCount =
	this->d->lodConvergenceCandidateCount();
    status.activeFaces = this->getActiveLodFaceCount();
    status.activeRenderCost = this->getActiveLodRenderCost();
    status.renderCostBudget = this->getCurrentLodRenderCostBudget();
    const BObolRetainedAllocationResult &allocation =
	this->d->lodRetainedAllocationCertificate;
    status.selectedPresentationCost = allocation.selectedPresentationCost;
    status.certifiedPresentationBudget =
	allocation.certifiedPresentationBudget;
    status.pixelDemandPresentationCost =
	allocation.pixelDemandPresentationCost;
    status.requestedPresentationBudget = allocation.requestedSceneBudget;
    status.maximumMarginalPresentationBudget =
	allocation.maximumMarginalBudget;
    status.maximumProtectedPresentationBudget =
	allocation.maximumProtectedBudget;
    status.pointProxyCandidateCount = allocation.pointProxyCandidateCount;
    status.allocationCertificatePlanSerial = allocation.allocationPlanSerial;

    const std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(
	    const_cast<BObolViewController *>(this));
    for (const SoBRLDatabaseSource *source : sources) {
	if (!source)
	    continue;
	const size_t available = static_cast<size_t>(
	    std::max(0, source->getCompactInstanceCount()));
	const size_t expected =
	    source->getCompactExpectedInstanceCount();
	const size_t denominator = std::max(available, expected);
	status.availableLeafCount =
	    available > SIZE_MAX - status.availableLeafCount ?
	    SIZE_MAX : status.availableLeafCount + available;
	status.expectedLeafCount =
	    denominator > SIZE_MAX - status.expectedLeafCount ?
	    SIZE_MAX : status.expectedLeafCount + denominator;
    }

    const BObolViewLodState *viewState =
	this->d->viewAttachment->getViewLodState();
    if (viewState) {
	status.allocationPlanSerial = viewState->activeCadAllocationPlan();
	viewState->convergencePayloadCounts(status.activePayloadCount,
	    status.satisfiedPayloadCount,
	    status.memoryLimitedPayloadCount);
	size_t subpixelOccurrences = 0;
	size_t structuralBoxes = 0;
	if (viewState->lastCadPresentationOccurrenceCoverage(
		subpixelOccurrences, structuralBoxes)) {
	    status.presentedSubpixelOccurrenceCount = subpixelOccurrences;
	    status.presentedStructuralBoxCount = structuralBoxes;
	}
	status.terminalOccurrenceFailureCount =
	    viewState->cadOccurrenceTerminalFailureCount();
	BObolViewLodState::CadGpuResourceStatus gpu;
	if (viewState->cadGpuResourceStatus(gpu)) {
	    status.gpuTrackedBufferBytes = gpu.trackedBufferBytes;
	    status.gpuOrdinaryPartBufferBytes =
		gpu.ordinaryPartBufferBytes;
	    status.gpuProgressiveCutBufferBytes =
		gpu.progressiveCutBufferBytes;
	    status.gpuProgressiveActiveCutBytes =
		gpu.progressiveActiveCutBytes;
	    status.gpuBatchBufferBytes = gpu.batchBufferBytes;
	    status.gpuTriangleAtlasAllocatedBytes =
		gpu.triangleAtlasAllocatedBytes;
	    status.gpuTriangleAtlasLiveBytes =
		gpu.triangleAtlasLiveBytes;
	    status.gpuTriangleAtlasConfiguredCapacityBytes =
		gpu.triangleAtlasConfiguredCapacityBytes;
	    status.gpuTriangleAtlasPartCount =
		gpu.triangleAtlasPartCount;
	    status.gpuTriangleAtlasPageCount =
		gpu.triangleAtlasPageCount;
	    status.gpuOrdinaryPartFullUploadBytes =
		gpu.ordinaryPartFullUploadBytes;
	    status.gpuOrdinaryPartSuffixUploadBytes =
		gpu.ordinaryPartSuffixUploadBytes;
	    status.gpuOrdinaryPartGpuCopyBytes =
		gpu.ordinaryPartGpuCopyBytes;
	    status.gpuOrdinaryPartLineageReuseCount =
		gpu.ordinaryPartLineageReuseCount;
	    status.gpuOrdinaryPartLineageReplacementCount =
		gpu.ordinaryPartLineageReplacementCount;
	    status.gpuTriangleAtlasFullUploadBytes =
		gpu.triangleAtlasFullUploadBytes;
	    status.gpuTriangleAtlasSuffixUploadBytes =
		gpu.triangleAtlasSuffixUploadBytes;
	    status.gpuTriangleAtlasLineageReuseCount =
		gpu.triangleAtlasLineageReuseCount;
	    status.gpuPressureProxyCount = gpu.pressureProxyCount;
	    status.gpuProgressiveEvictionCount =
		gpu.progressiveEvictionCount;
	    status.gpuTriangleAtlasReclamationCount =
		gpu.triangleAtlasReclamationCount;
	    status.gpuResourceSampleSerial = gpu.sampleSerial;
	    status.gpuMemoryPressure = gpu.memoryPressure;
	}
    }

    if (this->d->lodService) {
	status.pendingTasks =
	    this->d->lodService->pendingTaskCountForGeneration(
		this->d->lodActiveGeneration);
	const size_t delayed =
	    this->d->lodService->delayedTaskCountForGeneration(
		this->d->lodActiveGeneration);
	status.pendingTasks = delayed > SIZE_MAX - status.pendingTasks ?
	    SIZE_MAX : status.pendingTasks + delayed;
	status.inFlight =
	    this->d->lodService->executingTaskCountForGeneration(
		this->d->lodActiveGeneration);
	status.queuedResults =
	    this->d->lodService->queuedResultCountForGeneration(
		this->d->lodActiveGeneration);
	status.queuedCacheWrites =
	    this->d->lodService->queuedCacheWriteCountForGeneration(
		this->d->lodActiveGeneration);
	status.residentMeshBytes =
	    this->d->lodService->residentMeshBytesForDiagnostics();
	status.stableResidentMeshBytes =
	    this->d->lodService->
		stableResidentMeshBytesForDiagnostics();
	status.reservedResidentMeshGrowthBytes =
	    this->d->lodService->
		reservedResidentMeshGrowthBytesForDiagnostics();
	status.residentMeshLimitBytes =
	    this->d->lodService->getResidentMeshLimit();
	status.activeWorkingSetBytes =
	    this->d->lodService->activeWorkingSetBytesForDiagnostics();
	status.peakWorkingSetBytes =
	    this->d->lodService->peakWorkingSetBytesForDiagnostics();
	status.residentCompactionCount =
	    this->d->lodService->
		residentMeshCompactionCountForDiagnostics();
	status.residentCompactionPlanCurrent =
	    this->d->lodService->residentMeshCompactionPlanForDiagnostics(
		this->d->residentMeshConsumerId(),
		&status.residentCompactionPlanRevision,
		&status.residentCompactionCandidateCount);
    }

    const SbBool structuralPending =
	status.expectedLeafCount > status.availableLeafCount;
    const SbBool sourcePreparationPending =
	this->d->lodAvailabilityLedger.providerPendingCount() > 0 ||
	controller_lod_compact_inventory_incomplete(sources);
    const SbBool structuralDiscovery =
	!structuralPending &&
	status.visibleTargetCount == 0 &&
	status.activePayloadCount == 0 &&
	!this->d->progressiveProviders.empty() &&
	this->hasProgressiveWorkPending();
    const SbBool resultPending =
	this->hasPendingLodResults() || status.queuedResults > 0;
    /* Applied immutable results are not converged until a completed frame has
     * actually presented their batch.  The publication policy owns the
     * corresponding timer-or-frame liveness witness; expose that state here
     * so the HUD cannot report a ready view while richer data is resident but
     * still unpresented. */
    const SbBool publicationPending =
	this->d->lodPresentationTransaction.publicationPending();
    const BObolHostWorkSnapshot hostWork = this->getHostWorkSnapshot();
    /* The render request is the presentation transaction's level-triggered
     * host obligation.  A capacity-relevant request has not yet produced the
     * completed-frame evidence for the current retained controls, even when
     * no source, allocator, or publication cursor remains.  Omitting this
     * edge let waiters observe READY between a fractional-cut request and its
     * frame, then watch the supposedly fixed framebuffer refine during the
     * hold interval.  Presentation-only HUD repaints are intentionally not a
     * LoD obligation. */
    const SbBool capacityFramePending =
	(hostWork.flags & BOBOL_HOST_WORK_CAPACITY_SAMPLE) ? TRUE : FALSE;
    BObolLodControlRefinement::Inputs controlInputs;
    controlInputs.interaction =
	this->d->lodInteractionSession.active() != FALSE;
    controlInputs.inventory = structuralPending || structuralDiscovery;
    controlInputs.availability = sourcePreparationPending != FALSE;
    controlInputs.result = resultPending != FALSE;
    controlInputs.publication = publicationPending != FALSE;
    controlInputs.submission = this->d->lodSubmissionPass.active() != FALSE;
    controlInputs.submissionRescan = this->d->lodSubmissionPass.rescanPending() != FALSE;
    controlInputs.retainedAllocation =
	this->d->lodRetainedPass.refinementPending() ||
	this->d->lodRetainedPass.residencyPending();
    controlInputs.residentGrowth =
	this->d->lodAvailabilityLedger.residentGrowthPending();
    controlInputs.pointTriangleRecovery =
	this->d->lodPointQualityPhase.triangleRecoveryPending();
    controlInputs.structuralFrontier =
	status.presentedStructuralBoxCount >
	    status.terminalOccurrenceFailureCount;
    controlInputs.presentationReplay =
	this->d->lodInterruptedPresentationReplay.pending();
    controlInputs.presentationBarrier =
	this->hasPendingLodRefinementFrame() != FALSE;
    controlInputs.capacityFrame = capacityFramePending != FALSE;
    controlInputs.pointAdmissionFrame =
	this->d->lodPointAdmissionFrame.pending();
    controlInputs.pointCalibration =
	this->d->lodPointQualityPhase.presentationPending();
    controlInputs.capacityCalibration =
	this->d->lodAdmissionEvidence.capacity().rescanAfterFrame();
    controlInputs.headroomProbe =
	this->d->lodAdmissionEvidence.headroom().retryPending();
    controlInputs.handoff =
	this->d->lodPresentationPolicy.handoffPending();
    controlInputs.compaction = this->d->lodCompactionPolicy.pending();
    controlInputs.cacheWrite = status.queuedCacheWrites > 0;
    const BObolLodControlRefinement::Snapshot control =
	BObolLodControlRefinement::evaluate(controlInputs);
    const SbBool calibrationPending = control.calibrationPending() ? TRUE : FALSE;
    const SbBool controlPending = control.controlPending() ? TRUE : FALSE;
    status.refinementFramePending =
	(this->d->lodPresentationTransaction.barrierPending() ||
	 capacityFramePending) ? TRUE : FALSE;
    status.activeGeneration = this->d->lodActiveGeneration;
    status.submissionSourceIndex = this->d->lodSubmissionSourceIndex;
    status.submissionEntryOffset = this->d->lodSubmissionEntryOffset;
    status.budgetCalibrationPending =
	this->d->lodAdmissionEvidence.capacity().rescanAfterFrame() ||
	this->d->lodAdmissionEvidence.headroom().retryPending();
    status.stablePresentationHandoffPending =
	this->d->lodPresentationPolicy.handoffPending() ? TRUE : FALSE;
    status.pointProxyCalibrationPending =
	(this->d->lodPointAdmissionFrame.pending() ||
	 this->d->lodPointQualityPhase.pending()) ? TRUE : FALSE;
    status.pointProxyAdmissionFramePending =
	this->d->lodPointAdmissionFrame.pending() ? TRUE : FALSE;
    status.stablePointProxyCalibrationPending =
	this->d->lodPointQualityPhase.presentationPending() ? TRUE : FALSE;
    status.pointProxyTriangleRecoveryPending =
	this->d->lodPointQualityPhase.triangleRecoveryPending() ? TRUE : FALSE;
    status.residentGrowthReallocationPending =
	this->d->lodAvailabilityLedger.residentGrowthPending() ? TRUE : FALSE;
    status.publicationFramePending = publicationPending;
    status.sourcePreparationPending = sourcePreparationPending;
    status.sourcePreparationProviderCount =
	this->d->lodAvailabilityLedger.providerPendingCount();

    BObolLodConvergencePolicy::Inputs convergenceInputs;
    convergenceInputs.viewEpoch = this->d->lodViewRevision;
    convergenceInputs.policyEpoch = this->d->lodPolicyRevision;
    convergenceInputs.expectedLeafCount = status.expectedLeafCount;
    convergenceInputs.availableLeafCount = status.availableLeafCount;
    convergenceInputs.visibleTargetCount = status.visibleTargetCount;
    convergenceInputs.activePayloadCount = status.activePayloadCount;
    convergenceInputs.satisfiedPayloadCount = status.satisfiedPayloadCount;
    convergenceInputs.presentedSubpixelOccurrenceCount =
	status.presentedSubpixelOccurrenceCount;
    convergenceInputs.presentedStructuralBoxCount =
	status.presentedStructuralBoxCount;
    convergenceInputs.terminalOccurrenceFailureCount =
	status.terminalOccurrenceFailureCount;
    convergenceInputs.memoryLimitedPayloadCount =
	status.memoryLimitedPayloadCount;
    convergenceInputs.pendingTasks = status.pendingTasks;
    convergenceInputs.inFlight = status.inFlight;
    convergenceInputs.queuedCacheWrites = status.queuedCacheWrites;
    convergenceInputs.stableResidentMeshBytes =
	status.stableResidentMeshBytes;
    convergenceInputs.residentMeshLimitBytes =
	status.residentMeshLimitBytes;
    convergenceInputs.gpuTrackedBufferBytes =
	status.gpuTrackedBufferBytes;
    convergenceInputs.failedSourceCount = status.failedSourceCount;
    convergenceInputs.structuralDiscovery = structuralDiscovery != FALSE;
    convergenceInputs.visibilityCensusComplete =
	this->d->lodCoveragePolicy.hasCompleteVisibleCount();
    convergenceInputs.sourcePreparationPending =
	sourcePreparationPending != FALSE;
    convergenceInputs.submissionPending =
	this->d->lodSubmissionPass.active() != FALSE;
    convergenceInputs.resultPending = resultPending != FALSE;
    convergenceInputs.publicationPending = publicationPending != FALSE;
    convergenceInputs.calibrationPending = calibrationPending != FALSE;
    convergenceInputs.controlPending = controlPending != FALSE;
    convergenceInputs.interactive = this->d->lodInteractionSession.active() != FALSE;
    convergenceInputs.compactionPending =
	this->d->lodCompactionPolicy.pending();
    convergenceInputs.progressiveWorkPending =
	this->hasProgressiveWorkPending() != FALSE;
    convergenceInputs.gpuMemoryPressure =
	status.gpuMemoryPressure != FALSE;
    /* Offline terminal mode has deliberately opted out of responsiveness
     * admission.  Do not retain the preceding interactive budget's
     * "limited" diagnosis after its raster-stable target has settled: doing
     * so leaves the progress HUD over deterministic captures even though no
     * work or constraint remains. */
    convergenceInputs.stableBudgetLimited =
	!this->d->forceTerminalLodRefinement &&
	this->d->lodAdmissionEvidence.capacity().stableBudgetLimited();
    convergenceInputs.presentationLimited =
	!this->d->forceTerminalLodRefinement &&
	(this->d->lodInteractiveProgressiveCeiling >= 0 ||
	 this->d->lodPresentationPointProxyPixelThreshold > 1.01f ||
	 /* A rejected static-quality trial is completed-frame proof that the
	  * next richer population missed this view epoch's hard presentation
	  * deadline.  Its temporary renderer ceiling may already have been
	  * reconciled and removed, but that must not make the accepted terminal
	  * cut appear unconstrained to the HUD or qualification harness. */
	 this->d->lodStaticQualityTrial.capacityRejected());
    const BObolLodConvergencePolicy::Decision convergence =
	this->d->lodConvergencePolicy.evaluate(convergenceInputs);

    static_assert(
	static_cast<int>(BObolLodConvergencePolicy::Phase::IDLE) ==
	    BOBOL_LOD_CONVERGENCE_IDLE &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::DISCOVERING) ==
	    BOBOL_LOD_CONVERGENCE_DISCOVERING &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::INTERACTIVE) ==
	    BOBOL_LOD_CONVERGENCE_INTERACTIVE &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::REFINING) ==
	    BOBOL_LOD_CONVERGENCE_REFINING &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::CALIBRATING) ==
	    BOBOL_LOD_CONVERGENCE_CALIBRATING &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::BACKGROUND) ==
	    BOBOL_LOD_CONVERGENCE_BACKGROUND &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::CONVERROR) ==
	    BOBOL_LOD_CONVERGENCE_ERROR &&
	static_cast<int>(BObolLodConvergencePolicy::Phase::PREPARING) ==
	    BOBOL_LOD_CONVERGENCE_PREPARING,
	"public and private LoD convergence phases must agree");
    static_assert(
	static_cast<int>(BObolLodConvergencePolicy::Outcome::ACTIVE) ==
	    BOBOL_LOD_PRESENTATION_ACTIVE &&
	static_cast<int>(BObolLodConvergencePolicy::Outcome::READY) ==
	    BOBOL_LOD_PRESENTATION_READY &&
	static_cast<int>(BObolLodConvergencePolicy::Outcome::CONSTRAINED) ==
	    BOBOL_LOD_PRESENTATION_CONSTRAINED &&
	static_cast<int>(BObolLodConvergencePolicy::Outcome::ERROR) ==
	    BOBOL_LOD_PRESENTATION_ERROR,
	"public and private LoD presentation outcomes must agree");
    static_assert(
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::INTERACTION) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_INTERACTION &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::INVENTORY) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_INVENTORY &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::AVAILABILITY) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_AVAILABILITY &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::PUBLICATION) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_PUBLICATION &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::PLANNING) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_PLANNING &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::PRESENTATION) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_PRESENTATION &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::HANDOFF) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_HANDOFF &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::COMPACTION) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_COMPACTION &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::CACHE_WRITE) ==
	    BOBOL_LOD_CONTROL_OBLIGATION_CACHE_WRITE,
	"public and private LoD work ledgers must agree");
    static_assert(
	static_cast<int>(BObolLodControlRefinement::Owner::NONE) ==
	    BOBOL_LOD_CONTROL_OWNER_NONE &&
	static_cast<int>(BObolLodControlRefinement::Owner::INTERACTION) ==
	    BOBOL_LOD_CONTROL_OWNER_INTERACTION &&
	static_cast<int>(BObolLodControlRefinement::Owner::INVENTORY) ==
	    BOBOL_LOD_CONTROL_OWNER_INVENTORY &&
	static_cast<int>(BObolLodControlRefinement::Owner::AVAILABILITY) ==
	    BOBOL_LOD_CONTROL_OWNER_AVAILABILITY &&
	static_cast<int>(BObolLodControlRefinement::Owner::PUBLICATION) ==
	    BOBOL_LOD_CONTROL_OWNER_PUBLICATION &&
	static_cast<int>(BObolLodControlRefinement::Owner::PLANNING) ==
	    BOBOL_LOD_CONTROL_OWNER_PLANNING &&
	static_cast<int>(BObolLodControlRefinement::Owner::PRESENTATION) ==
	    BOBOL_LOD_CONTROL_OWNER_PRESENTATION &&
	static_cast<int>(BObolLodControlRefinement::Owner::HANDOFF) ==
	    BOBOL_LOD_CONTROL_OWNER_HANDOFF &&
	static_cast<int>(BObolLodControlRefinement::Owner::COMPACTION) ==
	    BOBOL_LOD_CONTROL_OWNER_COMPACTION &&
	static_cast<int>(BObolLodControlRefinement::Owner::CACHE_WRITE) ==
	    BOBOL_LOD_CONTROL_OWNER_CACHE_WRITE,
	"public and private LoD control owners must agree");
    static_assert(
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::OWNERLESS_WORK) ==
	    BOBOL_LOD_CONTROL_VIOLATION_OWNERLESS_WORK &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::TERMINAL_WITH_WORK) ==
	    BOBOL_LOD_CONTROL_VIOLATION_TERMINAL_WITH_WORK &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::INVALID_READINESS) ==
	    BOBOL_LOD_CONTROL_VIOLATION_INVALID_READINESS &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::INVALID_OWNER) ==
	    BOBOL_LOD_CONTROL_VIOLATION_INVALID_OWNER,
	"public and private LoD refinement violations must agree");
    status.phase = static_cast<int>(convergence.phase);
    status.outcome = static_cast<int>(convergence.outcome);
    status.controlObligationMask = control.obligations;
    status.controlOwner = static_cast<int>(control.owner);
    status.controlViolationMask = BObolLodControlRefinement::validate(
	control, convergence.terminal, convergence.viewReady,
	convergence.terminalError);
    status.fraction = convergence.fraction;
    status.terminal = convergence.terminal ? TRUE : FALSE;
    status.terminalError = convergence.terminalError ? TRUE : FALSE;
    status.viewReady = convergence.viewReady ? TRUE : FALSE;
    status.hasLodState = convergence.hasLodState ? TRUE : FALSE;
    status.backgroundPending =
	convergence.backgroundPending ? TRUE : FALSE;
    status.memoryLimited = convergence.memoryLimited ? TRUE : FALSE;
    status.performanceLimited =
	convergence.performanceLimited ? TRUE : FALSE;
}

double
BObolViewController::getCalibratedLodRenderCostPerSecond(void) const
{
    const long double value = this->d->lodInteractionSession.active() ?
	this->d->lodInteractiveCalibratedRenderCostPerSecond :
	this->d->lodStableCalibratedRenderCostPerSecond;
    return value >= static_cast<long double>(
	std::numeric_limits<double>::max()) ?
	std::numeric_limits<double>::max() : static_cast<double>(value);
}

double
BObolViewController::getInteractiveCalibratedLodRenderCostPerSecond(void) const
{
    const long double value =
	this->d->lodInteractiveCalibratedRenderCostPerSecond;
    return value >= static_cast<long double>(
	std::numeric_limits<double>::max()) ?
	std::numeric_limits<double>::max() : static_cast<double>(value);
}

double
BObolViewController::getStableCalibratedLodRenderCostPerSecond(void) const
{
    const long double value =
	this->d->lodStableCalibratedRenderCostPerSecond;
    return value >= static_cast<long double>(
	std::numeric_limits<double>::max()) ?
	std::numeric_limits<double>::max() : static_cast<double>(value);
}

void
BObolViewController::syncLodViewSignature(SbBool advanceOnChange)
{
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    SbBool haveCamera = this->getViewInfo(&view);
    const controller_lod_view_signatures signatures =
	controller_lod_view_signature(
	view, haveCamera, this->d->activeCamera, this->d->viewportRegion);

    const SbBool priorViewSignatureValid =
	this->d->lodViewSignature.has_value() ? TRUE : FALSE;
    const BObolLodViewScaleSnapshot priorScaleSignature =
	this->d->lodViewScaleSignature;
    /* advanceLodViewRevision() clears the view-local visibility census.  Its
     * transient zero must not discard an already proven aggregate point cut;
     * sample this presentation fact before invalidating the census. */
    const SbBool priorPointProxyAggregationApplicable =
	this->d->pointProxyAggregationApplicableForCameraTransition() ?
	    TRUE : FALSE;
    SbBool priorViewReady = FALSE;
    if (this->d->lodAutoSubmit && !this->d->lodInteractionSession.active()) {
	BObolLodConvergenceStatus priorStatus;
	this->getLodConvergenceStatus(priorStatus);
	priorViewReady = priorStatus.viewReady;
    }

    if (this->d->lodViewSignature &&
	this->d->lodViewSignature->same(signatures.view))
	return;

    const SbBool scaleChanged =
	this->d->lodViewSignature &&
	!this->d->lodViewScaleSignature.same(signatures.scale) ?
	    TRUE : FALSE;
    this->d->lodViewSignature = signatures.view;
    this->d->lodViewScaleSignature = signatures.scale;
    if (advanceOnChange) {
	/* Capture the prior scale quality frame before the new camera epoch
	 * clears its probe state.  If it completed within the 10 Hz hard
	 * responsiveness floor, the next wheel event must not immediately undo
	 * that coherent quality step. */
	const int priorQualityCeiling =
	    this->d->lodInteractiveProgressiveCeiling;
	int initialQualityFloor = -1;
	const BObolLodViewDemandPolicy::CameraChangeDecision demandChange =
	    this->d->lodViewDemandPolicy.observeCameraChange(
		scaleChanged != FALSE,
		this->d->lastRenderTimeNanoseconds);
	this->d->lodDiscretePopulationTrialPermit.revoke();
	const SbBool retainPriorQualityCeiling =
	    demandChange.retainPriorQualityCeiling ? TRUE : FALSE;
	this->advanceLodViewRevision();
	/* A camera may be edited through its public Coin fields rather than
	 * through syncCameraFromViewContext().  LoD submission is then the first
	 * place that observes the new signature.  Request the view frame here;
	 * previously isRenderRequested() happened to return true merely because
	 * progressive work existed, masking this missing presentation edge. */
	this->requestRender("lod-view");
	/* Camera bookkeeping is unconditional, but the 150 ms refinement pump
	 * belongs only to an active automatic LoD consumer.  Marking generic
	 * controllers progressive here leaves ordinary retained views
	 * permanently render-pending and suppresses later edge-triggered host
	 * wakeups. */
	if (this->d->lodAutoSubmit) {
	    const SbBool pointProxyAggregationForCameraChange =
		priorPointProxyAggregationApplicable ||
		this->d->pointProxyAggregationApplicable() ? TRUE : FALSE;
	    /*
	     * A bracketed mouse gesture has already entered interaction in
	     * beginLodInteraction(); an unbracketed wheel/trackpad epoch enters
	     * here.  Remember which case this camera change represents before
	     * entering the interaction-session state below.
	     */
	    const SbBool continuingInteractive =
		this->d->lodInteractionSession.active();
	    const SbBool reusePriorScalePresentation =
		!continuingInteractive && scaleChanged && priorViewReady &&
		this->d->lastRenderTimeNanoseconds > 0 &&
		this->d->lastRenderTimeNanoseconds <=
		    BObolLodViewDemandPolicy::
			qualityFrameDurationNanoseconds() ? TRUE : FALSE;
	    if (!continuingInteractive) {
		/* syncLodViewSignature has already installed the new signature by
		 * this point.  Preserve the values sampled above so a sequence of
		 * unbracketed wheel events can prove that it returned exactly to its
		 * starting scale rather than forcing a 50k-leaf stable retarget. */
		this->d->lodInteractionStartCertificate.capture(
		    priorScaleSignature,
		    priorViewSignatureValid != FALSE,
		    priorViewReady,
		    this->d->lodAdmissionEvidence.capacity().currentBudget());
		const BObolViewLodState *snapshotState =
		    this->d->viewAttachment->getViewLodState();
		size_t presentedPrimitives = 0;
		if (priorViewReady && snapshotState &&
		    this->d->lastRenderTimeNanoseconds > 0 &&
		    this->d->lastRenderTimeNanoseconds <=
			BObolLodViewDemandPolicy::
			    qualityFrameDurationNanoseconds() &&
		    snapshotState->lastCadPresentedPrimitiveCount(
			presentedPrimitives) && presentedPrimitives > 0) {
		    initialQualityFloor =
			snapshotState->maximumActiveProgressiveCut();
		    if (priorQualityCeiling >= 0 &&
			initialQualityFloor >= 0)
			initialQualityFloor = std::min(
			    initialQualityFloor, priorQualityCeiling);
		}
		this->d->lodPresentationPolicy.capturePrior(
		    this->d->lodTargetPixelError,
		    this->d->lodInteractiveProgressiveCeiling,
		    snapshotState ? snapshotState->
			cadPresentationProgressiveCutNextFraction() : 0.0f,
		    pointProxyAggregationForCameraChange ?
			this->d->lodPresentationPointProxyPixelThreshold : 1.0f,
		    controller_lod_presentation_population(snapshotState,
			this->d->lodViewQualityDomainRevision),
		    this->d->lodViewRevision);
		this->d->seedInteractiveCalibrationFromStable();
	    }
	    this->d->lodPoseContinuity.clearRetainOccurrenceCuts();
	    this->d->lodPlanningObligations.retireImportanceCensus();
	    this->d->lodViewDemandPolicy.beginCameraInteraction(
		!continuingInteractive, scaleChanged != FALSE);
	    if (!continuingInteractive)
		this->d->lodViewDemandPolicy.seedQualityFloor(
		    initialQualityFloor);
	    const int64_t now = bu_gettime();
	    this->d->lodInteractionSession.observeCameraChange(
		now, this->d->renderCompletionSerial);
	    this->d->lodPresentationPolicy.cancelHandoff();
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_MEASUREMENT);
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
	    this->d->lodRefinementNotBeforeMicroseconds = 0;
	    this->d->advanceAdmissionRevision(
	        BObolLodAdmissionRevisionDomain::CAPACITY);
	    this->d->lodCompactionPolicy.requestAfter(now, 750000);
	    /* Preserve a pending cut's frame gate across a newer camera
	     * signature.  The render ceiling may still respond immediately, but
	     * no path may expose another finer prefix before an intervening
	     * completed frame. */
	    /* A wheel/trackpad camera update is event-driven unless the host
	     * explicitly brackets it as a gesture.  In either case the gap since an
	     * unrelated repaint measures user/test cadence, not renderer capacity;
	     * feeding it to the motion policy can turn a responsive retained cut
	     * into a destructive coarse cut.  Render and GL timer measurements are
	     * the actionable capacity signals. */
	    const uint64_t endpointRenderSample =
		continuingInteractive &&
		this->d->lastRenderTimeNanoseconds > 0 ?
		this->d->lastRenderTimeNanoseconds :
		this->d->smoothedRenderTimeNanoseconds;
	    const uint64_t interactiveRenderSample = std::max(
		endpointRenderSample,
		this->d->lodLastCadGpuTimeNanoseconds);
	    this->d->lodTargetPixelError =
		reusePriorScalePresentation ? 1.0f :
		BObolLodQualityPolicy::interactivePixelError(
		    interactiveRenderSample,
		    this->d->lodInteractiveTargetFps);
	    BObolViewLodState *viewState =
		this->d->viewAttachment->getViewLodState();
	    int entryQualityFloor = -1;
	    if (scaleChanged && !continuingInteractive && viewState &&
		interactiveRenderSample > 0 &&
		interactiveRenderSample <= 66666667ULL) {
		entryQualityFloor = viewState->maximumActiveProgressiveCut();
		if (this->d->lodInteractiveProgressiveCeiling >= 0)
		    entryQualityFloor = std::min(entryQualityFloor,
			this->d->lodInteractiveProgressiveCeiling);
	    }
	    if (viewState)
		viewState->setCadPresentationCameraMotionFrameReuse(
		    TRUE);
	    int responsiveCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    const bool hasNewRenderFeedback =
		this->d->renderCompletionSerial !=
		    this->d->lodInteractiveCeilingFeedbackRenderSerial;
	    /*
	     * cameraMotionFrameReuse deliberately replays the already prepared
	     * occurrence classification.  Its frame therefore cannot validate
	     * a newly computed point-proxy threshold.  Repeatedly feeding those
	     * timings back into the point-proxy quality policy used to grow
	     * an unapplied 1 px cut to 40-64 px during one drag, then expose that
	     * fictitious correction all at once when quiet rendering resumed.
	     *
	     * Seed an unbracketed interaction once from the last actually
	     * presented stable cut.  beginLodInteraction() performs the same
	     * one-time seed for bracketed gestures.  During subsequent motion,
	     * only the progressive draw-count ceiling is safe to ratchet under
	     * prepared replay; quiet, ceiling-free frames calibrate point
	     * aggregation from cuts they really presented.
	     */
	    if (hasNewRenderFeedback && !continuingInteractive &&
		!reusePriorScalePresentation &&
		pointProxyAggregationForCameraChange)
		this->d->lodPresentationPointProxyPixelThreshold =
		    BObolLodQualityPolicy::pointProxyThreshold(
			this->d->lodPresentationPointProxyPixelThreshold,
			interactiveRenderSample,
			this->d->lodInteractiveTargetFps);
	    else if (!pointProxyAggregationForCameraChange) {
		this->d->lodPresentationPointProxyPixelThreshold = 1.0f;
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	    }
	    if (reusePriorScalePresentation) {
		/* The prior exact frame already proved this retained cut and
		 * one-pixel occurrence population below the hard zoom deadline.
		 * Preserve those renderer settings for the first new-scale frame.
		 * Reusing the standing atlas is both faster and visually continuous;
		 * only feedback measured at the new scale may lower it. */
		responsiveCeiling = priorQualityCeiling;
	    } else if (hasNewRenderFeedback) {
		const int activeMaximum = viewState ?
		    viewState->maximumActiveProgressiveCut() : -1;
		responsiveCeiling = BObolLodQualityPolicy::progressiveCeiling(
			activeMaximum, this->d->lodTargetPixelError,
			this->d->lodInteractiveProgressiveCeiling);
		this->d->lodInteractiveCeilingFeedbackRenderSerial =
		    this->d->renderCompletionSerial;
	    }
	    if (entryQualityFloor >= 0)
		responsiveCeiling = std::max(
		    responsiveCeiling, entryQualityFloor);
	    /*
	     * A successful scale-quality frame proves more than the ordinary
	     * motion target: its exact presentation completed inside the 100 ms
	     * interaction deadline.  BObolLodViewDemandPolicy retains that cut as
	     * the quality floor for the rest of the scale epoch and lowers it only
	     * after a real deadline miss.  Honor the retained proof on every wheel
	     * sample, not only on entry to the interaction.  Otherwise the 60 Hz
	     * motion heuristic can ratchet a known-responsive large mesh down one
	     * global PoP ordinal per event (Lucy fell from cut 24 to cut 18 while
	     * zooming out), then make quiet convergence visibly walk back through
	     * all of those already-resident populations.
	     *
	     * This is a renderer-only lower bound.  Occurrence allocation remains
	     * scene-budgeted, resident prefetch remains memory-bounded, and the
	     * completed/aborted-frame pressure paths call noteQualityMiss() before
	     * the next camera sample, so a genuinely unresponsive cut is still
	     * corrected immediately.
	     */
	    const int provenQualityFloor =
		this->d->lodViewDemandPolicy.qualityFloor();
	    if (scaleChanged && provenQualityFloor >= 0)
		responsiveCeiling = std::max(
		    responsiveCeiling, provenQualityFloor);
	    if (retainPriorQualityCeiling && priorQualityCeiling >= 0)
		responsiveCeiling = std::max(
		    responsiveCeiling, priorQualityCeiling);
	    /*
	     * A gesture's first frame is only the initial capacity probe.  If
	     * subsequent measured frames remain slow, ratchet the global
	     * render-only ceiling downward immediately.  Never raise it during
	     * continuous input: that would reintroduce oscillation and a slow
	     * frame.  The normal 150 ms quiet transition clears the ceiling.
	     */
	    if (hasNewRenderFeedback && responsiveCeiling >= 0 &&
		(this->d->lodInteractiveProgressiveCeiling < 0 ||
		 responsiveCeiling <
		    this->d->lodInteractiveProgressiveCeiling))
		this->d->lodInteractiveProgressiveCeiling =
		    responsiveCeiling;
	    if (viewState)
		viewState->setCadPresentationProgressiveCutCeiling(
		    this->d->lodInteractiveProgressiveCeiling);
	    if (viewState)
		viewState->setCadPresentationPointProxyPixelThreshold(
		    this->d->lodPresentationPointProxyPixelThreshold);
	    this->markProgressiveWorkPending();
	}
    }
}
