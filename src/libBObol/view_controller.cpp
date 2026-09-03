/*                 V I E W _ C O N T R O L L E R . C P P
 * BRL-CAD
 *
 * Copyright (c) 2026 United States Government as represented by
 * the U.S. Army Research Laboratory.
 */

#include "common.h"

#include "bu/assert.h"
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
#include "lod_result_authentication_private.h"
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

namespace {

/* Result publication runs on the host thread and yields between 256-result
 * atomic batches.  Give a quiet pump one bounded 16 ms owner-thread slice: a
 * smaller slice serialized a large ready queue behind a much more expensive
 * full-scene render after every few hundred results.  Active input remains
 * interruptible at every batch boundary, and the service's independent byte,
 * task, and result limits continue to bound resource use. */
constexpr uint64_t controller_default_lod_apply_microseconds = 16000;

static_assert(
    static_cast<int>(BObolLodCapacitySearchCertificate::Phase::INACTIVE) ==
	BOBOL_LOD_CAPACITY_SEARCH_INACTIVE &&
    static_cast<int>(BObolLodCapacitySearchCertificate::Phase::ALLOCATING) ==
	BOBOL_LOD_CAPACITY_SEARCH_ALLOCATING &&
    static_cast<int>(BObolLodCapacitySearchCertificate::Phase::PRESENTING) ==
	BOBOL_LOD_CAPACITY_SEARCH_PRESENTING &&
    static_cast<int>(BObolLodCapacitySearchCertificate::Phase::MEASURING) ==
	BOBOL_LOD_CAPACITY_SEARCH_MEASURING &&
    static_cast<int>(BObolLodCapacitySearchCertificate::Phase::TERMINAL) ==
	BOBOL_LOD_CAPACITY_SEARCH_TERMINAL,
    "public and private capacity-search phases must agree");
static_assert(
    static_cast<int>(BObolLodCapacitySearchCertificate::Goal::STEADY) ==
	BOBOL_LOD_CAPACITY_SEARCH_STEADY &&
    static_cast<int>(BObolLodCapacitySearchCertificate::Goal::STATIC) ==
	BOBOL_LOD_CAPACITY_SEARCH_STATIC,
    "public and private capacity-search goals must agree");

}

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
bool
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

bool
controller_lod_adaptive_point_aggregation_allowed(
    const BObolViewController *controller,
    bool staticQualityCapacityRejected)
{
    return BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	controller_lod_mesh_first_scene_safe(
	    controller_render_database_source_roots(controller)),
	staticQualityCapacityRejected);
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

/* Order the cuts which differ from a committed retained allocation by what
 * must remain recognizable, not by compact registry insertion order.
 * requestedCut is the last exact screen-space demand retained on each view
 * payload, so both the current projected error and the affected screen
 * footprint are available without a fresh scene projection.  Fixed buckets
 * keep construction O(retained payloads) and avoid the comparison sort which
 * made large-scene recovery a GUI-thread hotspot.  The returned application
 * plan is sparse: point-aggregated or already-matching occurrences must not
 * make every capacity candidate scan the complete source population again.
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
	const int allocatedCut = viewState->currentCadAllocatedCut(
	    payload, payload->viewRevision, payload->policyRevision,
	    source->getEffectiveLodDrawMode());
	const bool presentationApplied = allocatedCut >= 0 &&
	    viewState->cadAllocatedPresentationApplied(
		payload, payload->viewRevision, payload->policyRevision,
		source->getEffectiveLodDrawMode());
	if (allocatedCut < 0 || presentationApplied)
	    continue;
	const char *recoveryTrace = getenv("BOBOL_LOD_TRACE_RECOVERY_PLAN");
	if (recoveryTrace) {
	    static std::atomic<unsigned int> recoveryTraceCount(0);
	    const unsigned int traceIndex = recoveryTraceCount.fetch_add(1);
	    if (BU_STR_EQUAL(recoveryTrace, "all") || traceIndex < 64)
		bu_log("BObol LoD recovery occurrence path=%s active=%d "
		       "allocated=%d chunks=%zu/%zu prepared=%d "
		       "prepared_revision=%llu mesh_revision=%llu layers=%zu\n",
		       payload->sourcePath.getString(), payload->activeCut,
		       allocatedCut, payload->presentedChunks.size(),
		       payload->requiredChunks.size(),
		       payload->preparedCadGeometry ? 1 : 0,
		       static_cast<unsigned long long>(
			   payload->preparedCadGeometryRevision),
		       static_cast<unsigned long long>(payload->progressiveMesh ?
			   payload->progressiveMesh->revision() : 0),
		       payload->presentationLayers.size());
	}
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
    maxLodApplyMicroseconds(controller_default_lod_apply_microseconds),
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

SbBool
BObolHostWorkSnapshot::frameClaimed(void) const
{
    return (this->flags & BOBOL_HOST_WORK_FRAME_CLAIMED) ? TRUE : FALSE;
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
    this->sourcePreparationCompletedUnits = 0;
    this->sourcePreparationTotalUnits = 0;
    this->sourceAvailabilityChanged = 0;
    this->changed = 0;
    this->hasMore = 0;
}

BObolLodConvergenceStatus::BObolLodConvergenceStatus(void)
{
    this->clear();
}

BObolLodControlTraceState::BObolLodControlTraceState(void) :
    viewRevision(0),
    policyRevision(0),
    interactionActive(FALSE)
{
}

BObolLodControlTransitionRecord::BObolLodControlTransitionRecord(void) :
    serial(0),
    event(BOBOL_LOD_CONTROL_TRANSITION_UNNAMED)
{
}

const char *
bobol_lod_control_transition_event_name(
    BObolLodControlTransitionEvent event)
{
    switch (event) {
	case BOBOL_LOD_CONTROL_TRANSITION_INITIAL:
	    return "initial";
	case BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT:
	    return "external_input";
	case BOBOL_LOD_CONTROL_TRANSITION_INTERACTION:
	    return "interaction";
	case BOBOL_LOD_CONTROL_TRANSITION_INVENTORY:
	    return "inventory";
	case BOBOL_LOD_CONTROL_TRANSITION_AVAILABILITY:
	    return "availability";
	case BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION:
	    return "publication";
	case BOBOL_LOD_CONTROL_TRANSITION_PLANNING:
	    return "planning";
	case BOBOL_LOD_CONTROL_TRANSITION_PRESENTATION:
	    return "presentation";
	case BOBOL_LOD_CONTROL_TRANSITION_HANDOFF:
	    return "handoff";
	case BOBOL_LOD_CONTROL_TRANSITION_COMPACTION:
	    return "compaction";
	case BOBOL_LOD_CONTROL_TRANSITION_CACHE_WRITE:
	    return "cache_write";
	case BOBOL_LOD_CONTROL_TRANSITION_IDLE_SERVICE:
	    return "idle_service";
	case BOBOL_LOD_CONTROL_TRANSITION_UNNAMED:
	default:
	    return "unnamed";
    }
}

static bool
controller_lod_control_trace_state_equal(
    const BObolLodControlTraceState &left,
    const BObolLodControlTraceState &right)
{
    const BObolLodConvergenceStatus &a = left.convergence;
    const BObolLodConvergenceStatus &b = right.convergence;
    return left.hostWork.revision == right.hostWork.revision &&
	left.hostWork.renderRevision == right.hostWork.renderRevision &&
	left.hostWork.flags == right.hostWork.flags &&
	left.viewRevision == right.viewRevision &&
	left.policyRevision == right.policyRevision &&
	left.interactionActive == right.interactionActive &&
	a.phase == b.phase && a.outcome == b.outcome &&
	a.controlFactMask == b.controlFactMask &&
	a.controlObligationMask == b.controlObligationMask &&
	a.controlOwner == b.controlOwner &&
	a.controlViolationMask == b.controlViolationMask &&
	a.controlPresentationWitnessMask ==
	    b.controlPresentationWitnessMask &&
	a.constraintEvidenceMask == b.constraintEvidenceMask &&
	a.inventoryRevision == b.inventoryRevision &&
	a.availabilityRevision == b.availabilityRevision &&
	a.visibilityRevision == b.visibilityRevision &&
	a.viewRevision == b.viewRevision &&
	a.policyRevision == b.policyRevision &&
	a.capacityRevision == b.capacityRevision &&
	a.cadRevision == b.cadRevision &&
	a.residentDemandRevision == b.residentDemandRevision &&
	a.capacitySearchPhase == b.capacitySearchPhase &&
	a.capacitySearchGoal == b.capacitySearchGoal &&
	a.capacitySearchSamplesRemaining ==
	    b.capacitySearchSamplesRemaining &&
	a.capacitySearchMeasuredCandidates ==
	    b.capacitySearchMeasuredCandidates &&
	a.capacitySearchTotalMeasuredCandidates ==
	    b.capacitySearchTotalMeasuredCandidates &&
	a.capacitySearchCandidateLimit == b.capacitySearchCandidateLimit &&
	a.capacitySearchMaximumCandidates ==
	    b.capacitySearchMaximumCandidates &&
	a.capacitySearchSampleLimit == b.capacitySearchSampleLimit &&
	a.capacitySearchCompletedUnits == b.capacitySearchCompletedUnits &&
	a.capacitySearchTotalUnits == b.capacitySearchTotalUnits &&
	a.currentAllocationPlanSerial == b.currentAllocationPlanSerial &&
	a.presentationTransactionSerial == b.presentationTransactionSerial &&
	a.presentationRequiredRenderSerial ==
	    b.presentationRequiredRenderSerial &&
	a.presentedFrameSerial == b.presentedFrameSerial &&
	a.activeGeneration == b.activeGeneration &&
	a.submissionSourceIndex == b.submissionSourceIndex &&
	a.submissionEntryOffset == b.submissionEntryOffset &&
	a.expectedLeafCount == b.expectedLeafCount &&
	a.availableLeafCount == b.availableLeafCount &&
	a.visibleTargetCount == b.visibleTargetCount &&
	a.activePayloadCount == b.activePayloadCount &&
	a.satisfiedPayloadCount == b.satisfiedPayloadCount &&
	a.pendingTasks == b.pendingTasks && a.inFlight == b.inFlight &&
	a.queuedResults == b.queuedResults &&
	a.queuedCacheWrites == b.queuedCacheWrites &&
	a.rendererPreparationTargetSignature ==
	    b.rendererPreparationTargetSignature &&
	a.rendererPreparationTotalUnits == b.rendererPreparationTotalUnits &&
	a.rendererPreparationCompletedUnits ==
	    b.rendererPreparationCompletedUnits &&
	a.rendererPreparationRemainingUnits ==
	    b.rendererPreparationRemainingUnits &&
	a.rendererPreparationReservedBytes ==
	    b.rendererPreparationReservedBytes &&
	a.rendererPreparationTargetCount == b.rendererPreparationTargetCount &&
	a.rendererPreparationPreparingTargetCount ==
	    b.rendererPreparationPreparingTargetCount &&
	a.rendererPreparationConstrainedTargetCount ==
	    b.rendererPreparationConstrainedTargetCount &&
	a.rendererPreparationFailedTargetCount ==
	    b.rendererPreparationFailedTargetCount &&
	a.rendererPreparationInvalidTargetCount ==
	    b.rendererPreparationInvalidTargetCount &&
	a.residentCompactionCount == b.residentCompactionCount &&
	a.residentCompactionCandidateCount ==
	    b.residentCompactionCandidateCount &&
	a.sourcePreparationProviderCount ==
	    b.sourcePreparationProviderCount &&
	a.sourcePreparationCompletedUnits ==
	    b.sourcePreparationCompletedUnits &&
	a.sourcePreparationTotalUnits == b.sourcePreparationTotalUnits &&
	memcmp(&a.fraction, &b.fraction, sizeof(a.fraction)) == 0 &&
	a.terminal == b.terminal &&
	a.terminalError == b.terminalError && a.viewReady == b.viewReady &&
	a.hasLodState == b.hasLodState &&
	a.backgroundPending == b.backgroundPending;
}

static BObolLodControlTransitionEvent
controller_lod_control_owner_event(
    const BObolLodControlTraceState &before,
    const BObolLodControlTraceState &after)
{
    int owner = before.convergence.controlOwner;
    if (owner == BOBOL_LOD_CONTROL_OWNER_NONE)
	owner = after.convergence.controlOwner;
    switch (owner) {
	case BOBOL_LOD_CONTROL_OWNER_INTERACTION:
	    return BOBOL_LOD_CONTROL_TRANSITION_INTERACTION;
	case BOBOL_LOD_CONTROL_OWNER_INVENTORY:
	    return BOBOL_LOD_CONTROL_TRANSITION_INVENTORY;
	case BOBOL_LOD_CONTROL_OWNER_AVAILABILITY:
	    return BOBOL_LOD_CONTROL_TRANSITION_AVAILABILITY;
	case BOBOL_LOD_CONTROL_OWNER_PUBLICATION:
	    return BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION;
	case BOBOL_LOD_CONTROL_OWNER_PLANNING:
	    return BOBOL_LOD_CONTROL_TRANSITION_PLANNING;
	case BOBOL_LOD_CONTROL_OWNER_PRESENTATION:
	    return BOBOL_LOD_CONTROL_TRANSITION_PRESENTATION;
	case BOBOL_LOD_CONTROL_OWNER_HANDOFF:
	    return BOBOL_LOD_CONTROL_TRANSITION_HANDOFF;
	case BOBOL_LOD_CONTROL_OWNER_COMPACTION:
	    return BOBOL_LOD_CONTROL_TRANSITION_COMPACTION;
	case BOBOL_LOD_CONTROL_OWNER_CACHE_WRITE:
	    return BOBOL_LOD_CONTROL_TRANSITION_CACHE_WRITE;
	case BOBOL_LOD_CONTROL_OWNER_NONE:
	default:
	    return BOBOL_LOD_CONTROL_TRANSITION_IDLE_SERVICE;
    }
}

static BObolLodControlTransitionEvent
controller_lod_take_pending_transition_event(std::atomic<int> &pending)
{
    const int value = pending.exchange(BOBOL_LOD_CONTROL_TRANSITION_UNNAMED,
	std::memory_order_acq_rel);
    return value > BOBOL_LOD_CONTROL_TRANSITION_UNNAMED &&
	value <= BOBOL_LOD_CONTROL_TRANSITION_IDLE_SERVICE ?
	static_cast<BObolLodControlTransitionEvent>(value) :
	BOBOL_LOD_CONTROL_TRANSITION_UNNAMED;
}

void
BObolLodConvergenceStatus::clear(void)
{
    this->phase = BOBOL_LOD_CONVERGENCE_IDLE;
    this->outcome = BOBOL_LOD_PRESENTATION_READY;
    this->controlFactMask = 0;
    this->controlObligationMask = BOBOL_LOD_CONTROL_OBLIGATION_NONE;
    this->controlOwner = BOBOL_LOD_CONTROL_OWNER_NONE;
    this->controlViolationMask = BOBOL_LOD_CONTROL_VIOLATION_NONE;
    this->controlPresentationWitnessMask =
	BOBOL_LOD_PRESENTATION_WITNESS_NONE;
    this->constraintEvidenceMask = BOBOL_LOD_CONSTRAINT_NONE;
    this->viewQualityHistoryEntryCount = 0;
    this->viewQualityHistoryRememberCount = 0;
    this->viewQualityHistoryRecallCount = 0;
    this->inventoryRevision = 0;
    this->availabilityRevision = 0;
    this->visibilityRevision = 0;
    this->viewRevision = 0;
    this->policyRevision = 0;
    this->capacityRevision = 0;
    this->cadRevision = 0;
    this->residentDemandRevision = 0;
    this->capacitySearchPhase = 0;
    this->capacitySearchGoal = 0;
    this->capacitySearchSamplesRemaining = 0;
    this->capacitySearchMeasuredCandidates = 0;
    this->capacitySearchTotalMeasuredCandidates = 0;
    this->capacitySearchCandidateLimit = 0;
    this->capacitySearchMaximumCandidates = 0;
    this->capacitySearchSampleLimit = 0;
    this->capacitySearchCompletedUnits = 0;
    this->capacitySearchTotalUnits = 0;
    this->currentAllocationPlanSerial = 0;
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
    this->terminalProxyOccurrenceCount = 0;
    this->terminalOccurrenceFailureCount = 0;
    this->pendingTasks = 0;
    this->inFlight = 0;
    this->queuedResults = 0;
    this->queuedCacheWrites = 0;
    this->rendererPreparationTargetSignature = 0;
    this->rendererPreparationTotalUnits = 0;
    this->rendererPreparationCompletedUnits = 0;
    this->rendererPreparationRemainingUnits = 0;
    this->rendererPreparationReservedBytes = 0;
    this->rendererPreparationTargetCount = 0;
    this->rendererPreparationPreparingTargetCount = 0;
    this->rendererPreparationConstrainedTargetCount = 0;
    this->rendererPreparationFailedTargetCount = 0;
    this->rendererPreparationInvalidTargetCount = 0;
    this->presentedPrimitiveCountValid = FALSE;
    this->presentedPrimitiveCount = 0;
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
    this->reachablePointProxyCandidateCount = 0;
    this->selectedPointProxyCount = 0;
    this->prominentCandidateCount = 0;
    this->prominentQualityFloorViolationCount = 0;
    this->maximumNormalizedVisualError =
	std::numeric_limits<double>::infinity();
    this->visualImportanceDebt = 0.0;
    this->committedAllocationPlanSerial = 0;
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
    this->semanticPresentationFramePending = FALSE;
    this->sourcePreparationPending = FALSE;
    this->sourcePreparationProviderCount = 0;
    this->sourcePreparationCompletedUnits = 0;
    this->sourcePreparationTotalUnits = 0;
    this->failedSourceCount = 0;
}

BObolProgressiveProviderRecord::BObolProgressiveProviderRecord(void) :
    token(0),
    callback(NULL),
    userData(NULL),
    userDataFree(NULL),
    sourcePreparationCompletedUnits(0),
    sourcePreparationTotalUnits(0)
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

    bool samePlanningInputs(
	const std::vector<BObolLodCoordinator::LodSourceSnapshot>
	    &other) const
    {
	return this->sameInventories(other) &&
	    this->sameVisibilityInputs(other);
    }

    bool sameVisibilityInputs(
	const std::vector<BObolLodCoordinator::LodSourceSnapshot>
	    &other) const
    {
	if (!this->sameIdentities(other))
	    return false;
	for (size_t i = 0; i < this->sources.size(); ++i) {
	    if (this->sources[i].visibilityRevision !=
		    other[i].visibilityRevision)
		return false;
	}
	return true;
    }
};

static bool
controller_lod_source_changed_entries(
    const BObolLodCoordinator::LodSourceSnapshot &current,
    const BObolLodCoordinator::LodSourceSnapshot &previous,
    std::vector<size_t> &entryIndices, SbBool &coverageInvalidated)
{
    entryIndices.clear();
    coverageInvalidated = FALSE;
    if (!current.source || !current.sameIdentity(previous))
	return false;

    if (current.inventoryRevision != previous.inventoryRevision) {
	std::vector<size_t> inventoryEntries;
	if (!previous.inventoryRevision.value() ||
	    !current.source->getDisplayMeshLodChangedEntries(
		previous.inventoryRevision.value(), inventoryEntries,
		&coverageInvalidated))
	    return false;
	entryIndices.insert(entryIndices.end(), inventoryEntries.begin(),
	    inventoryEntries.end());
    }
    if (current.visibilityRevision != previous.visibilityRevision) {
	std::vector<size_t> visibilityEntries;
	if (!previous.visibilityRevision ||
	    !current.source->getDisplayMeshLodVisibilityChangedEntries(
		previous.visibilityRevision, visibilityEntries))
	    return false;
	entryIndices.insert(entryIndices.end(), visibilityEntries.begin(),
	    visibilityEntries.end());
    }
    if (entryIndices.empty())
	return false;
    std::sort(entryIndices.begin(), entryIndices.end());
    entryIndices.erase(std::unique(entryIndices.begin(), entryIndices.end()),
	entryIndices.end());
    return true;
}

/* Keep every consumer of source-population counts on the same domain as the
 * submission planner.  Ordinary analytic/wire sources may have a compact
 * scene index for picking and tree synchronization without owning any mesh
 * LoD work; counting those entries in convergence creates a denominator no
 * LoD action can satisfy. */
static bool
controller_lod_source_has_planning_contract(
    const SoBRLDatabaseSource *source)
{
    if (!source || !source->getDatabase())
	return false;
    return source->hasDisplayLodTargets() ||
	(source->realizationStatus.getValue() ==
	     SoBRLDatabaseSource::REALIZED &&
	 !source->needsRealization() && source->hasRealizedMeshGeometry());
}

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
    BObolViewLodState *viewState = controller->getViewLodState();
    for (size_t i = 0; i < sources.size(); i++) {
	SoBRLDatabaseSource *source = sources[i];
	if (!source || !source->getDatabase())
	    continue;
	if (viewState)
	    (void)viewState->synchronizeResidentCadProgressiveSource(source);
	const SbBool hasCurrentCompactTargets =
	    source->hasDisplayLodTargets();
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
		       "entries=%d targets_current=%d requests_current=%d "
		       "resident_progressive=%zu realized_mesh_current=%d "
		       "status=%d needs_realization=%d source_revision=%u "
		       "inputs_revision=%u lod_revision=%llu\n",
		       source->path.getValue().getString(),
		       source->hasCompactInstanceIndex() ? 1 : 0,
		       source->getCompactInstanceCount(),
		       hasCurrentCompactTargets ? 1 : 0,
		       source->hasDisplayMeshLodRequests() ? 1 : 0,
		       source->getDisplayResidentProgressiveGeometryCount(),
		       hasCurrentRealizedMesh ? 1 : 0,
		       source->realizationStatus.getValue(),
		       source->needsRealization() ? 1 : 0,
		       source->sourceRevision.getValue(),
		       source->inputsRevision.getValue(),
		       static_cast<unsigned long long>(
			   source->getDisplayMeshLodRevision()));
	}
	/*
	 * A provider request is published only after the stream merge validates
	 * the detached/live source epoch.  Resident progressive parts are already
	 * immutable and require no provider epoch.  Both are therefore ready for
	 * view planning while source-wide realization is still running.  Requiring
	 * terminal realization here serializes cold work behind complete leaf
	 * discovery and strands native progressive wire at its minimum cut.
	 */
	if (!controller_lod_source_has_planning_contract(source))
	    continue;

	BObolLodCoordinator::LodSourceSnapshot snapshot;
	snapshot.source = source;
	snapshot.database = source->getDatabase();
	snapshot.routingId.set(source->getCompactSourceRoutingId());
	snapshot.inventoryRevision.set(source->getDisplayMeshLodRevision());
	snapshot.visibilityRevision =
	    source->getDisplayMeshLodVisibilityRevision();
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

    /* A result-ready notification may already be in flight when a service is
     * replaced.  The old service owns no successor transition. */
    if (service != controller->d->lodService)
	return;

    const uint64_t generation = controller->d->lodActiveGeneration;
    const uint64_t consumerId =
	controller->d->residentMeshConsumerId();
    const bool generationReady = generation != 0 &&
	service->queuedResultCountForGeneration(generation) > 0;

    /* An explicitly begun generation remains a valid manual workflow when
     * auto-submit is disabled.  Resident compaction is automatic policy work,
     * however, and must not resurrect the pump after policy retirement. */
    const bool compactionReady =
	controller->automaticLodControlEnabled() &&
	service->queuedResidentMeshCompactionResultCountForDiagnostics(
	    consumerId) > 0;
    /* Reclaiming an asset which no consumer still demands deliberately
     * produces no drawable compaction result.  It does, however, advance the
     * admission epoch and may unblock a memory-limited current asset.  Keep
     * this comparison atomic: result-ready callbacks execute on a worker
     * thread while the presentation owner advances the observation cursor. */
    const bool residentCapacityReady =
	controller->automaticLodControlEnabled() &&
	service->residentMeshAdmissionRevision() !=
	    controller->d->lodResidentAdmissionRevision.load(
		std::memory_order_acquire);
    if (!generationReady && !compactionReady && !residentCapacityReady)
	return;

    if (generationReady) {
	controller->d->lodAvailabilityLedger.noteResultsReady(bu_gettime());
    }
    controller->d->lodControlPendingExternalEvent.store(
	BOBOL_LOD_CONTROL_TRANSITION_PUBLICATION,
	std::memory_order_release);
    controller->publishProgressiveWorkPending();
}

void
BObolViewController::setLodService(BObolLodService *service)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
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
    const bool automaticLod =
	this->automaticLodControlEnabled() != FALSE;
    this->d->lodCoveragePolicy.setRequired(
	service && automaticLod);
    this->d->lodResultSubscriberId = 0;
    this->d->lodAvailabilityLedger.resetResultQueue();
    this->d->lodActiveGeneration = 0;
    this->d->rewindLodSubmissionCursor();
    this->d->lodSubmissionPass.retire();
    this->d->lodSubmissionIntent.reset();
    this->d->lodViewDemandPolicy.reset();
    this->d->lodInteractionStartCertificate.reset();
    this->d->lodRetainedViewContinuity.reset();
	this->d->lodPlanningObligations.reset();
    this->d->lodAvailabilityLedger.resetResidentGrowth();
    this->d->lodAvailabilityLedger.commitInventoryDelta();
    this->d->lodDiscretePopulationTrialPermit.revoke();
    this->d->lodStaticQualityTrial.reset();
    this->d->lodCoveragePolicy.reset();
    if (service && automaticLod)
	this->d->lodCoveragePolicy.activate(true);
    this->d->resetRetainedPassAnnotations();
    this->d->lodStructuralRepair.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
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
    this->d->lodResidentAdmissionRevision.store(
	service ? service->residentMeshAdmissionRevision() : 0,
	std::memory_order_release);
    this->d->lodResidentAdmissionRetryRevision = 0;
    this->d->lodCompactionPolicy.resetForServiceChange(
	service != NULL, bu_gettime(), 750000);
    this->d->lodPresentationTransaction.reset();
    this->d->lodRefinementNotBeforeMicroseconds = 0;
    this->d->lodInterruptedPresentationReplay.retire();
    this->d->lodExactPresentationFrame.reset();

    if (this->d->lodService)
	this->d->lodResultSubscriberId =
	    this->d->lodService->subscribeResultReady(
		BObolViewController::lodResultReadyCB, this);

    (void)this->synchronizeAutomaticLodControl();
}

uint64_t
BObolViewController::beginLodGeneration(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
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
    this->d->publishCadDiscoveryPointProxyThreshold(
	BObolViewController::Impl::POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
    this->d->lodPointAdmissionFrame.retire();
    if (requestFrame && this->automaticLodControlEnabled() &&
	this->d->lodPresentationPointProxyPixelThreshold + 1.0e-6f <
	    oldEffective) {
	this->requestLodCapacityRender("lod-discovery-point-reset");
    }
}

void
BObolViewController::cancelActiveLodGeneration(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    if (this->d->lodService && this->d->lodActiveGeneration != 0) {
	this->d->lodService->cancelGeneration(this->d->lodActiveGeneration);
    }
    this->d->lodActiveGeneration = 0;
    this->d->rewindLodSubmissionCursor();
    this->d->lodSubmissionPass.retire();
    this->d->lodSubmissionIntent.reset();
    this->d->lodViewDemandPolicy.reset();
    this->d->lodDiscretePopulationTrialPermit.revoke();
    this->d->lodAvailabilityLedger.resetResidentGrowth();
    this->d->lodAvailabilityLedger.commitInventoryDelta();
    this->d->lodPlanningObligations.retireResidentAdmissionRetry();
    this->d->resetRetainedPassAnnotations();
    this->d->lodStructuralRepair.reset();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
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
    this->d->lodInterruptedPresentationReplay.retire();
    this->d->lodExactPresentationFrame.reset();
    this->d->resetCadPresentationLimits();
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
	    this->d->publishCadCameraMotionFrameReuse(FALSE);
	}
}

SbBool
BObolViewController::lodViewPolicyEnabled(void) const
{
    if (!this->d->viewAttachment)
	return FALSE;

    struct bv_lod_policy policy;
    bv_lod_policy_init(&policy);
    this->d->viewAttachment->getLodPolicy(&policy);
    return policy.policy != BV_LOD_OFF && policy.mesh_enabled ? TRUE : FALSE;
}

SbBool
BObolViewController::automaticLodControlEnabled(void) const
{
    return this->d->lodAutoSubmit && this->lodViewPolicyEnabled() ?
	TRUE : FALSE;
}

SbBool
BObolViewController::synchronizeAutomaticLodControl(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    const bool coverageWasRequired =
	this->d->lodCoveragePolicy.required();
    const bool automatic = this->automaticLodControlEnabled() != FALSE;
    const bool coverageRequired = automatic && this->d->lodService;
    this->d->lodCoveragePolicy.setRequired(coverageRequired);
    if (!automatic) {
	this->retireAutomaticLodControl();
	return FALSE;
    }
    if (coverageRequired && !coverageWasRequired) {
	this->d->clearLodConvergenceCandidates();
	this->d->lodCoveragePolicy.activate(true);
    }
    return TRUE;
}

void
BObolViewController::retireAutomaticLodControl(void)
{
    BObolLodControlTransitionScope controlTransition(this);
    this->d->retireAutomaticLodControl();
    this->retireLodCapacityRenderRequest();
    BU_ASSERT(this->d->automaticLodControlRetired() &&
	!this->getHostWorkSnapshot().capacitySampleRequested());
    BObolViewLodState *viewState = this->d->viewAttachment ?
	this->d->viewAttachment->getViewLodState() : NULL;
    if (viewState)
	viewState->setCadPresentationCameraMotionFrameReuse(FALSE);

    /* The shared pump also owns database realization.  Clear its old LoD
     * edge only when no independent provider has reported work remaining. */
    if (this->d->lodAvailabilityLedger.providerPendingCount() == 0)
	this->clearProgressiveWorkPending();
}

void
BObolViewController::invalidateDatabaseSourceLodState(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
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
    this->d->lodCoveragePolicy.reset();
    if (this->synchronizeAutomaticLodControl() && this->d->lodService)
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    const SbBool requested = enabled ? TRUE : FALSE;
    if (this->d->lodAutoSubmit != requested)
	this->d->resetLodViewQualityHistory();
    this->d->lodAutoSubmit = requested;
    this->d->lodStaticQualityTrial.reset();
    if (this->synchronizeAutomaticLodControl()) {
	this->requestLodCapacityRender("lod-auto-submit");
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    if (cut < 0)
	cut = 0;

    if (this->d->lodForcedCut && *this->d->lodForcedCut == cut)
	return;

    this->d->resetLodViewQualityHistory();
    this->d->lodForcedCut = cut;
    this->advanceLodPolicyRevision();
    this->requestLodCapacityRender("lod-policy");
}

void
BObolViewController::clearLodForcedCut(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    if (!this->d->lodForcedCut)
	return;

    this->d->resetLodViewQualityHistory();
    this->d->lodForcedCut.reset();
    this->advanceLodPolicyRevision();
    this->requestLodCapacityRender("lod-policy");
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
    /* Report only work whose next transition is a completed presentation.
     * Allocation-stage handoff and capacity certificates are planning debt;
     * classifying them as frame debt repaints an unchanged population and can
     * hide a missing producer.  Budget, point, headroom, and presentation-
     * stage handoff probes do require a frame and remain level triggered. */
    return this->d->lodPresentationFramePending() ? TRUE : FALSE;
}

size_t
BObolViewController::processPendingLodResults(size_t maxResults,
	uint64_t maxMicroseconds)
{
    BObolLodControlTransitionScope controlTransition(this);
    if (!this->d->lodService)
	return 0;

    if (!this->hasPendingLodResults() &&
	this->d->lodService->queuedResultCountForGeneration(
	    this->d->lodActiveGeneration) == 0) {
	this->synchronizeProgressiveWorkPending();
	return 0;
    }

    if (maxMicroseconds == 0) {
	(void)this->applyLodResults(this->d->lodService, maxResults);
	this->synchronizeProgressiveWorkPending();
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
    this->synchronizeProgressiveWorkPending();
    return processed;
}

int
BObolViewController::submitLodRequestsIfNeeded(SbBool refreshMissing,
	SbBool resetExisting)
{
    BObolLodControlTransitionScope controlTransition(this);
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.retire();
	this->d->lodSubmissionDelta.reset();
	this->d->lodStructuralRepair.reset();
	/* The empty-contract edge retires the source identity domain, not merely
	 * the current submission cursor.  A Z/redraw may pass through this state
	 * before publishing a replacement source; preserving the old dense census
	 * then adds the replacement population to an unreachable source entry. */
	this->d->lodCoveragePolicy.reset();
	this->d->clearLodConvergenceCandidates();
	this->d->lodProjectedDemandCache.clear();
	this->d->lodAvailabilityLedger.resetResidentGrowth();
	this->d->lodAvailabilityLedger.commitInventoryDelta();
	this->d->lodPlanningObligations.reset();
	this->d->lodViewDemandPolicy.clearDemandRefresh();
	this->d->lodLastSubmittedSources.clear();
	return 0;
    }

    if (this->d->lodLastSubmittedViewRevision == this->d->lodViewRevision &&
	this->d->lodLastSubmittedPolicyRevision == this->d->lodPolicyRevision &&
	signatures.samePlanningInputs(this->d->lodLastSubmittedSources)) {
	if (!this->d->lodSubmissionPass.active() &&
	    this->d->lodSubmissionPass.rescanPending()) {
	    const std::vector<SoBRLDatabaseSource *> sources =
		controller_render_database_source_roots(this);
	    const bool inventorySettled =
		this->d->lodAvailabilityLedger.providerPendingCount() == 0 &&
		!controller_lod_compact_inventory_incomplete(sources);
	    if (inventorySettled &&
		this->d->lodSubmissionPass.beginPendingRescan()) {
		this->d->rewindLodSubmissionCursor();
		this->d->lodSubmissionDelta.reset();
		this->d->resetRetainedPassAnnotations();
		if (this->d->lodCoveragePolicy.active())
		    this->d->lodCoveragePolicy.clearPassCounters();
	    }
	}
	BObolLodViewDemandPolicy::DemandPassInputs demandPassInputs;
	demandPassInputs.submissionActive =
	    this->d->lodSubmissionPass.active();
	demandPassInputs.rescanPending =
	    this->d->lodSubmissionPass.rescanPending();
	demandPassInputs.selectivePass = this->d->lodSubmissionDelta.active();
	    demandPassInputs.structuralRepair =
		this->d->lodStructuralRepair.active();
	    demandPassInputs.retainedAllocation =
		this->d->lodSubmissionIntent.retainedAdmission();
	    demandPassInputs.strongerOwnerPending =
		this->d->lodAdmissionEvidence.capacity().
		    capacityTransactionPending() ||
		this->d->lodPresentationTransaction.barrierPending() ||
		this->d->lodPresentationTransaction.publicationPending() ||
		this->d->lodAvailabilityLedger.residentGrowthPending() ||
		this->d->lodPresentationPolicy.handoffPending() ||
		this->d->lodExactPresentationFrame.pending() ||
		this->d->lodPointQualityPhase.pending();
	if (this->d->lodViewDemandPolicy.demandPassRequired(demandPassInputs)) {
	    /* A demand refresh is a complete ordinary current-view pass.  It may
	     * outlive the policy-revision edge which first armed it when a stronger
	     * presentation transaction retires that pass.  Recreate the runnable
	     * cursor here instead of waiting for an unrelated camera/source change.
	     * Do not manufacture a new evidence epoch: the existing policy revision
	     * is precisely the domain this pass must prove. */
	    this->d->rewindLodSubmissionCursor();
	    this->d->resetRetainedPassAnnotations();
	    this->d->lodSubmissionIntent.configure(
		refreshMissing != FALSE, resetExisting != FALSE);
	    this->d->lodSubmissionPass.beginFresh();
	}
	if (this->d->lodPlanningObligations.
		exactVisibilityReallocationReady(
		    this->d->lodSubmissionPass.active(),
		    this->d->lodExactPresentationFrame.pending(),
		    this->d->lodStructuralRepair.active(),
		    this->d->lodAdmissionEvidence.capacity().
			capacityTransactionPending(),
		    this->d->lodPresentationTransaction.barrierPending())) {
	    /* The selective source pass and its exact successor frame must publish
	     * and classify the new hidden-instance set before minimax examines the
	     * population.  Run exactly one complete allocation afterward; starting
	     * it on the revision edge sees the predecessor visibility state and can
	     * strand restored entries. */
	    this->d->lodPlanningObligations.
		retireExactVisibilityReallocation();
	    this->d->requestRetainedReallocation();
	    this->beginSceneWideCapacitySubmission();
	}
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
    const bool visibilityChanged =
	!this->d->lodLastSubmittedSources.empty() &&
	!signatures.sameVisibilityInputs(this->d->lodLastSubmittedSources);
    const bool planningInputsChanged =
	!this->d->lodLastSubmittedSources.empty() &&
	!signatures.samePlanningInputs(this->d->lodLastSubmittedSources);
    const bool viewOrPolicyChanged =
	this->d->lodLastSubmittedViewRevision != this->d->lodViewRevision ||
	this->d->lodLastSubmittedPolicyRevision !=
	    this->d->lodPolicyRevision;
    const bool visibilitySupersedesRetainedAllocation = visibilityChanged &&
	!sourceSetChanged && !inventoryChanged && !viewOrPolicyChanged &&
	this->d->lodSubmissionPass.active() &&
	this->d->lodSubmissionIntent.retainedAdmission();
    if (visibilitySupersedesRetainedAllocation) {
	/* A retained allocation publishes only at its final epoch-checked commit.
	 * Supersede that unpublished plan so the sparse visibility journal is
	 * consumed first.  The successor allocation then observes the newest
	 * hidden-instance set; allowing the old plan to finish can omit entries
	 * restored by a second edit. */
	this->d->lodSubmissionPass.retire();
	this->d->lodSubmissionIntent.setRetainedAdmission(false);
	this->d->lodRetainedAllocationTransaction.reset();
	this->d->resetRetainedPassAnnotations();
    }
    /* A contract appearing after the explicit empty-source state has no
     * predecessor against which sourceSetChanged can compare.  Re-arm the
     * authoritative coverage pass here; initial service setup is already
     * active, so this is idempotent before the first observation. */
    if (this->d->lodLastSubmittedSources.empty() &&
	this->d->lodCoveragePolicy.required() &&
	!this->d->lodCoveragePolicy.active())
	this->d->lodCoveragePolicy.activate(true);
    if (sourceSetChanged || planningInputsChanged || viewOrPolicyChanged)
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
	!visibilityChanged &&
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
    if (sourceSetChanged || planningInputsChanged || viewOrPolicyChanged)
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
    if ((sourceSetChanged || planningInputsChanged) &&
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
	    const BObolLodCoordinator::LodSourceSnapshot &currentSource =
		signatures.sources[currentIndex];
	    const BObolLodCoordinator::LodSourceSnapshot *previousSource =
		knownSource ?
		    &this->d->lodLastSubmittedSources[previousIndex] : NULL;
	    if (sameIdentity && previousSource &&
		currentSource.inventoryRevision ==
		    previousSource->inventoryRevision &&
		currentSource.visibilityRevision ==
		    previousSource->visibilityRevision)
		continue;
	    if (!sameIdentity || !previousSource) {
		pendingDeltaNeedsFullRescan = true;
		continue;
	    }

	    SoBRLDatabaseSource *changedSource = currentSource.source;
	    std::vector<size_t> changedEntries;
	    SbBool coverageInvalidated = FALSE;
	    if (!controller_lod_source_changed_entries(
		    currentSource, *previousSource, changedEntries,
		    coverageInvalidated) ||
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
    if ((sourceSetChanged || planningInputsChanged) &&
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
	    const BObolLodCoordinator::LodSourceSnapshot &currentSource =
		signatures.sources[currentIndex];
	    const BObolLodCoordinator::LodSourceSnapshot *previousSource =
		knownSource ?
		    &this->d->lodLastSubmittedSources[previousIndex] : NULL;
	    const bool planningDiffers = !previousSource ||
		currentSource.inventoryRevision !=
		    previousSource->inventoryRevision ||
		currentSource.visibilityRevision !=
		    previousSource->visibilityRevision;
	    if (sameIdentity && !planningDiffers)
		continue;

	    SoBRLDatabaseSource *changedSource = currentSource.source;
	    const bool alreadyTargeted =
		this->d->lodSubmissionDelta.targets(changedSource);
	    sourceDeltaFirst = std::min(sourceDeltaFirst, currentIndex);
	    if (sameIdentity && previousSource) {
		std::vector<size_t> changedEntries;
		SbBool coverageInvalidated = FALSE;
		if (controller_lod_source_changed_entries(
			currentSource, *previousSource, changedEntries,
			coverageInvalidated) &&
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
    } else if ((sourceSetChanged || planningInputsChanged) &&
	!this->d->lodSubmissionPass.active()) {
	this->d->lodSubmissionDelta.reset();
    }
    const bool hasExactSourceDelta =
	useSourceDelta || this->d->lodSubmissionDelta.active();
    const bool exactVisibilityDelta = visibilityChanged &&
	!sourceSetChanged && !inventoryChanged && !viewOrPolicyChanged &&
	hasExactSourceDelta && !sourceDeltaInvalidatesCoverage &&
	!pendingDeltaNeedsFullRescan && priorCoverageComplete;
    const bool sourceCoverageInvalidated = sourceSetChanged ||
	(planningInputsChanged &&
	 (!hasExactSourceDelta || sourceDeltaInvalidatesCoverage ||
	  pendingDeltaNeedsFullRescan || !priorCoverageComplete));
    if (this->d->lodStructuralRepair.pointRelaxationPending() &&
	BObolLodStructuralRepair::pointRelaxationDomainChanged(
	    viewOrPolicyChanged, visibilityChanged,
	    sourceCoverageInvalidated)) {
	/* A private candidate belongs to one exact occurrence and projection
	 * domain.  Keep immutable resident results produced by its own preload,
	 * including their inventory revisions, but cancel it when that semantic
	 * domain changes underneath the transaction. */
	this->d->lodStructuralRepair.cancelPointRelaxation();
    }
    if (sourceSetChanged || planningInputsChanged) {
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
	this->d->lodRetainedViewContinuity.clearHandoff();
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
	    if (pointRelaxationRequired) {
		/* The source population changed, so the retained point cut first
		 * needs a classifier frame for the new structural occurrence set.
		 * Stable timing of the preceding population is not transferable. */
		this->requestStructuralPointAdmissionFrame(
		    "lod-source-point-reclassification");
	    } else {
		this->d->lodPointQualityPhase.completeCalibration();
	    }
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
	this->d->positionLodSubmissionCursor(
	    useSourceDelta ? sourceDeltaFirst : 0);
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
	if (exactVisibilityDelta) {
	    /* Visibility changes the occurrence allocation, not the immutable
	     * scene or the renderer's measured cost capacity.  The changed compact
	     * instance set must nevertheless reach one exact framebuffer before
	     * minimax can use its point/box classification.  An erase-time
	     * allocation may have retired meshes which become visible again on
	     * redraw; allocating directly from the selective census sees neither
	     * those missing meshes nor their structural presentation and can leave
	     * a nonterminal view with no producer.
	     *
	     * Keep a terminal capacity certificate, require the presentation
	     * prerequisite, and reallocate its budget only after that frame.  A
	     * missing or stale delta journal never reaches this branch and retains
	     * the conservative full invalidation below. */
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::VISIBILITY);
	    this->d->lodRetainedViewContinuity.beginStableMutation(
		this->getActiveLodMeshPayloadCount() > 0);
	    this->d->requireExactPresentationFrame();
	    this->d->lodPlanningObligations.
		requestExactVisibilityReallocation();
	} else {
	    this->d->lodPlanningObligations.
		retireExactVisibilityReallocation();
	    this->d->applyAdmissionEvidenceAction(
		BObolLodAdmissionPlanner::EvidenceAction::
		    INVALIDATE_CAPACITY_MEASUREMENT);
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::INVENTORY);
	}
	if (useSourceDelta)
	    this->d->lodCoveragePolicy.clearPassCounters();
	if (sourceSetChanged || planningInputsChanged)
	    this->d->applyAdmissionEvidenceAction(
	        BObolLodAdmissionPlanner::EvidenceAction::RESET_CAPACITY_OVERLOAD);
    } else if (!extendedPendingDelta) {
	this->d->lodSubmissionPass.requestRescan();
    }
    if (pendingDeltaNeedsFullRescan)
	this->d->lodSubmissionPass.requestRescan();
    this->d->lodSubmissionPass.activate();
    if (!exactVisibilityDelta)
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
    BObolLodControlTransitionScope controlTransition(this);
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
	    this->d->lodAdmissionEvidence.capacity().capacityAllocationPending(),
	    this->d->lodAdmissionEvidence.capacity().presentationFramePending(),
	    submissionPresentationState &&
		submissionPresentationState->hasCadPresentationAssemblies())) {
	/* This is an ordinary producer barrier, not a submission failure.  The
	 * convergence snapshot and HUD carry its live status; diagnostics are
	 * reserved for actionable failures or terminal constraints. */
	return 0;
    }

    if (!this->d->lodSubmissionPass.active()) {
	this->d->rewindLodSubmissionCursor();
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	this->d->lodCoveragePolicy.clearPassCounters();
	this->d->lodSubmissionPass.beginFresh();
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
    const bool adaptivePointAggregationAllowed =
	BObolLodAdmissionPlanner::adaptivePointAggregationAllowed(
	    meshFirstSceneSafe,
	    this->d->lodStaticQualityTrial.capacityRejected());
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
    const bool capacitySearchOwnsAllocation =
	this->d->lodAdmissionEvidence.capacity().capacitySearch().
	    awaitingSample();
    const size_t presentationReconciliationBudget =
	this->d->lodPresentationPolicy.allocationReconciliationBudget(
	    capacitySearchOwnsAllocation);
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
	 * next useful cut can be just beyond either ordinary cadence.  Use the
	 * same 10 Hz hard floor which decides whether that cut remains visible.
	 * Exact hierarchy costs still prevent an unexpectedly huge cut jump. */
	/* The global progressive-ordinal staircase is valid only for one visible
	 * occurrence, but the static frame allowance is not single-object policy.
	 * A many-part scene spends that allowance through the occurrence-local
	 * importance allocator.  Keeping its target at the ordinary quiet cadence
	 * made an accepted static phase recompute the identical constrained plan
	 * indefinitely instead of using (or rejecting) its bounded static budget. */
	const bool retainedStaticPresentation =
	    this->d->lodRetainedViewContinuity.startCapacityAtStatic();
	const SbBool hardDeadlinePresentation =
	    !this->d->lodInteractionSession.active() &&
	    (this->d->lodStructuralRepair.active() ||
	     this->d->lodStaticQualityTrial.blocksNewTrial() ||
	     retainedStaticPresentation) ? TRUE : FALSE;
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
	    this->d->lodRendererPerformanceEvidence.cadPresentationNanosecondsAt(
		this->d->lodPresentationPointProxyPixelThreshold);
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
	size_t sourceLodTargetCount = 0;
	size_t sourceLogicalOccurrenceCount = 0;
	for (SoBRLDatabaseSource *source : sources) {
	    if (!source)
		continue;
	    const size_t count = source->getDisplayLodTargetCount();
	    sourceLodTargetCount = count > SIZE_MAX - sourceLodTargetCount ?
		SIZE_MAX : sourceLodTargetCount + count;
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
	     sourceLodTargetCount == 1);
	budgetInputs.hardDeadlinePresentation =
	    hardDeadlinePresentation != FALSE;
	budgetInputs.preserveActivePopulation =
	    this->d->lodRetainedViewContinuity.retainOccurrenceCuts() ||
	    retainedStaticPresentation;
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
	std::max(
	    this->d->lodStaticQualityTrial.acceptedPresentationCostFor(
		this->d->admissionRevisionStamp()),
	    this->d->lodStaticQualityTrial.constrainedPresentationBudgetFor(
		this->d->admissionRevisionStamp())));
    if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
	    this->d->lodViewRevision.value())) {
	const BObolLodCapacitySearchCertificate &search =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	bu_log("BObol LoD admission input phase=%u candidate=%zu "
	       "samples_remaining=%u current_budget=%zu cursor=%d "
	       "cursor_refinement=%zu cursor_retained=%d "
	       "target_fps=%.3f hard_deadline=%d retain_cuts=%d "
	       "retained_static=%d\n",
	       static_cast<unsigned int>(search.phase()),
	       search.candidateBudget(), search.samplesRemaining(),
	       this->d->lodAdmissionEvidence.capacity().currentBudget(),
	       this->d->lodAdmissionCursor.initialized() ? 1 : 0,
	       this->d->lodAdmissionCursor.refinementRemaining(),
	       this->d->lodAdmissionCursor.retainedAdmission() ? 1 : 0,
	       targetFps, hardDeadlinePresentation ? 1 : 0,
	       this->d->lodRetainedViewContinuity.retainOccurrenceCuts() ? 1 : 0,
	       retainedStaticPresentation ? 1 : 0);
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
	const bool retainedAdmissionRequested =
	    budget.retainedAdmission && retainedPopulationSettled &&
	    !this->d->lodCoveragePolicy.demandCensusRequired();
	const bool retainAdmissionMode =
	    this->d->lodSubmissionIntent.retainedAdmissionForPass(
		retainedAdmissionRequested,
		this->d->lodSubmissionPass.active());
	/* A sparse unsatisfied-refinement plan is not a retained-recovery plan.
	 * Reusing it lets its first subset consume the complete upgrade allowance;
	 * the later all-occurrence pass then sees zero and normalizes everything
	 * else to minimum.  Give an idle mode transition an explicit plan epoch
	 * and restart at the first source.  Once that retained cursor starts, its
	 * mode is immutable until completion or explicit semantic invalidation. */
	if (retainAdmissionMode &&
	    !this->d->lodSubmissionIntent.retainedAdmission()) {
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionIntent.setRetainedAdmission(true);
	} else if (!retainAdmissionMode &&
	    this->d->lodSubmissionIntent.retainedAdmission()) {
	    this->d->rewindLodSubmissionCursor();
	    this->d->lodSubmissionIntent.setRetainedAdmission(false);
	    this->d->resetRetainedAdmissionQualityProof();
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
	const BObolLodCapacitySearchCertificate &capacitySearch =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	/* A bounded search owns a frozen numeric population domain from its first
	 * candidate through its terminal certificate.  Timing EMA updates are
	 * observations used by that search, not permission to manufacture another
	 * allocation under the same admission revision. */
	const bool capacitySearchOwnsAllowances =
	    capacitySearch.phase() !=
		BObolLodCapacitySearchCertificate::Phase::INACTIVE;
	const size_t activeProtectedFloorBudget =
	    this->d->lodAdmissionEvidence.capacity().
		retainedQualityFloorBudget();
	size_t maximumProtectedBudget =
	    std::max(this->d->lodAdmissionEvidence.capacity().currentBudget(),
		activeProtectedFloorBudget);
	size_t maximumMarginalBudget = maximumProtectedBudget;
	if (!this->d->lodInteractionSession.active() &&
	    !capacitySearchOwnsAllowances &&
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
	    /* A complete mesh-first profile has already bounded source count,
	     * encoded bytes, largest asset, and aggregate working set.  Let its
	     * ordinary marginal tail use the same one-shot static allowance as an
	     * atomic prominent floor.  Restricting that allowance to prominent
	     * candidates can leave a handful of inexpensive ordinary cuts below
	     * pixel demand with no legal successor, even though their completed
	     * frame fits the retained static-image contract. */
	    if (meshFirstSceneSafe)
		maximumMarginalBudget = std::max(
		    maximumMarginalBudget, hardQualityBudget);
	}
	/* The atomic protected floor may be larger than the measured static-frame
	 * allowance even though many valuable marginal improvements fit.  During
	 * the bounded static phase, let the occurrence-local allocator spend that
	 * allowance instead of falling back to the ordinary quiet-frame budget.
	 * After a rejected floor, the same path keeps the floor disabled while the
	 * deadline ceiling prevents a retry of the failed population. */
	if (!this->d->lodInteractionSession.active() &&
	    !capacitySearchOwnsAllowances &&
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
	uint64_t protectedFloorIdentity = 0;
	const int64_t retainedAllocationStarted = bu_gettime();
	BObolViewLodState *retainedViewState =
	    this->d->viewAttachment->getViewLodState();
	BObolRetainedAllocationInputs inputs;
	inputs.sources = &sources;
	inputs.viewState = retainedViewState;
	/* Charge retained work outside occurrence-local allocation, including
	 * compact progressive cuts selected by the bounded source pass.  Deriving
	 * this value from a completed frame feeds point aggregation and renderer
	 * ceilings back into the next allocation key.  View state updates total and
	 * managed occurrence cost together, so their difference stays invariant as
	 * this plan changes active cuts.  Compact-cut mutations invalidate coverage,
	 * preventing the certificate from outliving a real external-cost change. */
	inputs.externalPresentationCost = retainedViewState ?
	    retainedViewState->allocationUnmanagedRenderCost() : 0;
	inputs.sceneBudget = this->d->lodAdmissionEvidence.capacity().currentBudget();
	inputs.maximumMarginalBudget = maximumMarginalBudget;
	const bool protectedFloorEligibleForSearch =
	    capacitySearch.phase() ==
		BObolLodCapacitySearchCertificate::Phase::INACTIVE ||
	    capacitySearch.goal() ==
		BObolLodCapacitySearchCertificate::Goal::STATIC;
	inputs.allowProtectedFloor = !this->d->lodInteractionSession.active() &&
	    protectedFloorEligibleForSearch &&
	    !this->d->lodAdmissionEvidence.capacity().retainedQualityFloorRejected() &&
	    this->d->lodPresentationPolicy.
		handoffReconciliationBudget() == 0;
	inputs.maximumProtectedBudget = inputs.allowProtectedFloor ?
	    maximumProtectedBudget : 0;
	if (presentationReconciliationBudget)
	    inputs.setPresentationReconciliationBudget(
		presentationReconciliationBudget);
	inputs.revisionStamp = this->d->admissionRevisionStamp();
    inputs.residentAdmissionRevision =
	bobol_retained_allocation_resident_admission_revision(
	    retainedViewState, this->d->lodService ?
		this->d->lodService->residentMeshAdmissionRevision() : 0);
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
	const bool reuseCommittedSerial = retainedViewState &&
	    priorAllocation.allocationPlanSerial != 0 &&
	    priorAllocation.allocationPlanSerial ==
		retainedViewState->activeCadAllocationPlan();
	const bool reuseAllocationInputs =
	    priorAllocation.inputKey() == inputs.inputKey() ||
	    priorAllocation.pixelDemandInputEquivalent(inputs);
	const bool reuseCommittedAllocation = reuseCommittedSerial &&
	    reuseAllocationInputs;
	const bool allocationPopulationCurrent = retainedViewState &&
	    priorAllocation.allocationPlanSerial != 0 &&
	    retainedViewState->cadAllocationPlanCoversCurrentPopulation(
		priorAllocation.allocationPlanSerial,
		inputs.viewRevision(), inputs.policyRevision(),
		priorAllocation.fixedCadPresentationCost);
	const bool reuseCoveredAllocation = reuseCommittedAllocation &&
	    allocationPopulationCurrent;
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
		       "serial_current=%d inputs_current=%d population_current=%d "
		       "budget=%zu/%zu marginal=%zu/%zu protected=%zu/%zu "
		       "allow_protected=%d/%d external=%zu/%zu "
		       "selected=%zu demand=%zu unresolved=%zu "
		       "point_threshold=%.3f/%.3f resident_admission=%llu/%llu "
		       "plan=%llu active_plan=%llu "
		       "search_phase=%u search_goal=%u floor_rejected=%d "
		       "trial_rejected=%d reconciliation=%zu "
		       "view=%llu policy=%llu\n",
		       static_cast<unsigned int>(allocationStatus),
		       reuseCoveredAllocation ? 1 : 0,
		       reuseCommittedSerial ? 1 : 0,
		       reuseAllocationInputs ? 1 : 0,
		       allocationPopulationCurrent ? 1 : 0,
		       inputs.sceneBudget, priorAllocation.requestedSceneBudget,
		       inputs.maximumMarginalBudget,
		       priorAllocation.maximumMarginalBudget,
		       inputs.maximumProtectedBudget,
		       priorAllocation.maximumProtectedBudget,
		       inputs.allowProtectedFloor ? 1 : 0,
		       priorAllocation.allowProtectedFloor ? 1 : 0,
		       inputs.externalPresentationCost,
		       priorAllocation.externalPresentationCost,
		       allocation.selectedPresentationCost,
		       allocation.pixelDemandPresentationCost,
		       allocation.unresolvedViewDependentPayloadCount,
		       static_cast<double>(inputs.pointProxyPixelThreshold),
		       static_cast<double>(
			   priorAllocation.pointProxyPixelThreshold),
		       static_cast<unsigned long long>(
			   inputs.residentAdmissionRevision),
		       static_cast<unsigned long long>(
			   priorAllocation.residentAdmissionRevision),
		       static_cast<unsigned long long>(
			   allocation.allocationPlanSerial),
		       static_cast<unsigned long long>(retainedViewState ?
			   retainedViewState->activeCadAllocationPlan() : 0),
		       static_cast<unsigned int>(capacitySearch.phase()),
		       static_cast<unsigned int>(capacitySearch.goal()),
		       this->d->lodAdmissionEvidence.capacity().
			   retainedQualityFloorRejected() ? 1 : 0,
		       this->d->lodStaticQualityTrial.capacityRejected() ? 1 : 0,
		       presentationReconciliationBudget,
		       static_cast<unsigned long long>(inputs.viewRevision()),
		       static_cast<unsigned long long>(inputs.policyRevision()));
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
	} else if (allocation.unresolvedViewDependentPayloadCount > 0) {
	    /* Allocation is downstream of the current-view projection census.  A
	     * policy edge can invalidate those stamps without removing the useful
	     * retained meshes.  Keep the typed refresh frontier across bounded
	     * source windows, but not as a valid presentation certificate.  Accepting
	     * the partial population would starve the source pass; discarding it
	     * would rescan the scene merely to rediscover the same stale records.
	     *
	     * The current retained pass consumes projectionRefreshPlans below.  Its
	     * result cannot itself be the terminal allocation: it was calculated
	     * before those payloads acquired current projection evidence.  Preserve
	     * one explicit successor obligation so the completed refresh pass is
	     * followed by exactly one allocation over the refreshed population.
	     * Without this edge, a large visibility-only redraw can retire the
	     * mechanical pass with every queue empty while the controller still has
	     * no current terminal-quality certificate. */
	    this->d->lodRetainedAllocationCertificate = allocation;
	    /* COMPLETE freezes a transaction's discovery result.  The payload
	     * refresh does not change the allocator input key, so retaining this
	     * object would replay the same refresh request instead of discovering
	     * the newly current projections.  The typed frontier above remains the
	     * durable request; retire only its completed discovery transaction. */
	    this->d->lodRetainedAllocationTransaction.reset();
	    this->d->lodPlanningObligations.
		requestExactVisibilityReallocation();
	} else {
	    /* A discovery-time pressure sample may raise the point threshold
	     * before the complete retained population is known.  Once the allocator
	     * and an exact renderer census prove that no occurrence can use that
	     * representation, canonicalize the inert control immediately.  The
	     * allocation alone sees only mesh-backed candidates; resetting while
	     * structural fallbacks are aggregated re-exposes those boxes and
	     * restarts the same structural-distribution transaction indefinitely. */
	    BObolViewLodState::CadStructuralProjectionHistogram
		structuralProjection;
	    const bool exactStructuralProjection = retainedViewState &&
		retainedViewState->lastCadStructuralProjectionHistogram(
		    structuralProjection);
	    const size_t terminalFailureCount = retainedViewState ?
		retainedViewState->cadOccurrenceTerminalFailureCount() : 0;
	    if (BObolLodAdmissionPlanner::pointProxyThresholdInert(
		    allocation.pointProxyCandidateCount,
		    exactStructuralProjection,
		    structuralProjection.visibleCount, terminalFailureCount) &&
		this->d->lodPresentationPointProxyPixelThreshold > 1.01f) {
		this->d->publishCadPointProxyThreshold(
		    BObolViewController::Impl::
			POINT_PROXY_PIXEL_THRESHOLD_MINIMUM);
		this->d->applyAdmissionEvidenceAction(
		    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
		this->d->lodPointQualityPhase.reset();
		allocation.pointProxyPixelThreshold = 1.0f;
		this->d->lodRetainedAllocationTransaction.reset();
	    }
	    this->d->lodRetainedAllocationCertificate = allocation;
	    protectedFloorBudget = allocation.protectedFloorBudget;
	    protectedFloorIdentity = allocation.protectedFloorIdentity;
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
	const bool allocationSucceeded =
	    allocationStatus != BOBOL_RETAINED_ALLOCATION_FAILED &&
	    allocation.unresolvedViewDependentPayloadCount == 0;
	this->d->lodRetainedAdmissionQualityEvidence.remember(
	    allocationSucceeded ? allocation.maximumNormalizedError :
		std::numeric_limits<double>::infinity(),
	    allocationSucceeded ? allocation.maximumProjectedErrorPixels :
		std::numeric_limits<double>::infinity(),
	    this->d->lodViewRevision.value(), this->d->lodPolicyRevision.value(),
	    this->d->lodPresentationPointProxyPixelThreshold);
	this->d->setRetainedQualityFloor(
	    protectedFloorBudget, protectedFloorIdentity, sceneActiveCost,
	    sceneMinimumActiveCost);
	if (getenv("BOBOL_LOD_TRACE_BUDGET") ||
	    getenv("BOBOL_LOD_TRACE_ALLOCATOR"))
	    bu_log("BObol LoD retained importance ceiling=%.6f "
		   "budget=%zu selected=%zu certified=%zu "
		   "marginal_limit=%zu protected_limit=%zu reused=%d "
		   "elapsed_us=%lld\n",
		   this->d->lodRetainedAdmissionQualityEvidence.
		       maximumNormalizedError(),
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
	this->d->lodRetainedViewContinuity.visibilityCensusDeferred()) {
	/* The motion fast path may have advanced the source cursor to its end
	 * without projecting any occurrences.  Start the first quiet pass from
	 * the beginning and discard any partial coverage counters so its terminal
	 * denominator belongs wholly to this camera. */
	this->d->rewindLodSubmissionCursor();
	this->d->lodCoveragePolicy.clearPassCounters();
	this->d->lodRetainedViewContinuity.completeVisibilityCensus();
    }
    for (size_t i = this->d->lodSubmissionSourceIndex;
	 i < sources.size();) {
	const size_t capacity = service->availableResultTaskCapacity();
	if (!capacity) {
	    break;
	}
	SoBRLDatabaseSource *source = sources[i];
	if (!source) {
	    this->d->positionLodSubmissionCursor(++i);
	    continue;
	}
	if (this->d->lodSubmissionDelta.active() &&
	    !this->d->lodSubmissionDelta.targets(source)) {
	    this->d->positionLodSubmissionCursor(++i);
	    continue;
	}
	/* A compact source records its source-backed PoP requests explicitly,
	 * making absence authoritative and cheap to test.  Non-compact sources
	 * may still contain direct SoBRLMeshShape children whose request
	 * metadata is owned by the shape, so let the submit action inspect those
	 * sources.  This preserves the custom scene-node boundary
	 * without rescanning terminal compact analytic meshes on every view. */
	if (source->hasCompactInstanceIndex() &&
	    !source->hasDisplayLodTargets()) {
	    this->d->positionLodSubmissionCursor(++i);
	    continue;
	}

	struct db_i *dbip = source->getDatabase();
	if (!dbip) {
	    append_controller_lod_diagnostic(this->d->lastLodDiagnostics,
					     source->path.getValue(),
					     "database source has no database for LoD submission");
	    this->d->positionLodSubmissionCursor(++i);
	    continue;
	}

	BObolViewLodState *viewLodState =
	    this->d->viewAttachment->getViewLodState();
	const size_t sourceMeshRequests =
	    source->getDisplayMeshLodRequestCount();
	const size_t sourceLodTargets = source->getDisplayLodTargetCount();
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
	    this->d->lodRetainedViewContinuity.deferVisibilityCensus();
	    this->d->positionLodSubmissionCursor(++i);
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
	     controller_lod_global_preview_requested(sourceLodTargets)) ?
		TRUE : FALSE;
	size_t presentedSubpixelOccurrences = 0;
	size_t presentedStructuralBoxes = 0;
	const bool exactStructuralPresentation = viewLodState &&
	    viewLodState->lastCadPresentationFrameExact() &&
	    viewLodState->lastCadPresentationOccurrenceCoverage(
		presentedSubpixelOccurrences, presentedStructuralBoxes);
	const size_t terminalOccurrenceFailures = viewLodState ?
	    viewLodState->cadOccurrenceTerminalFailureCount() : 0;
	const bool deferStructuralFallbackAdmission =
	    BObolLodAdmissionPlanner::structuralFallbackAdmissionDeferred(
		!this->d->lodInteractionSession.active(),
		this->d->lodCoveragePolicy.coverageComplete(),
		exactStructuralPresentation, presentedStructuralBoxes,
		terminalOccurrenceFailures, globalColdPreview,
		this->d->lodStructuralRepair.active());
	const SbBool structuralCoverageOnly =
	    BObolLodAdmissionPlanner::structuralCoverageOnly(
		this->d->forceTerminalLodRefinement != FALSE,
		this->d->lodCoveragePolicy.active(),
		this->d->lodViewDemandPolicy.demandRefreshActive(),
		globalColdPreview != FALSE,
		deferStructuralFallbackAdmission) ?
		TRUE : FALSE;
	action.setStructuralCoverageOnly(structuralCoverageOnly);
	action.setStructuralPresentationRepair(
	    this->d->lodStructuralRepair.active() ? TRUE : FALSE);
	action.setStructuralTerminalProxy(
	    this->d->lodStructuralRepair.terminalProxy() ? TRUE : FALSE);
	action.setPointRelaxationPreload(
	    this->d->lodStructuralRepair.active() &&
	    this->d->lodStructuralRepair.pointRelaxationPending() ?
		TRUE : FALSE);
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
	     !this->d->lodRetainedViewContinuity.retainOccurrenceCuts()) ||
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
	 * after a rotation even though the resident data and quiet budget could
	 * afford the pre-motion quality.  Downgrade permission above remains the
	 * independent guard which prevents pose changes from needlessly making an
	 * already useful cut coarser. */
	const SbBool retainedRefinementAllowed =
	    (this->d->forceTerminalLodRefinement ||
	     !this->d->lodInteractionSession.active() || scaleDemandChanged) &&
	    !this->d->lodAvailabilityLedger.residencyDrainActive() &&
	    !this->d->lodPresentationTransaction.barrierPending() ? TRUE : FALSE;
	action.setAllowRetainedRefinement(retainedRefinementAllowed);
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
	const bool residentAdmissionRetry =
	    this->d->lodPlanningObligations.residentAdmissionRetryPending();
	action.setAllowResidentPrefetchPastAllocation(
	    BObolLodAvailabilityScheduler::
		residentPrefetchPastAllocationAllowed(
		    scaleInteraction != FALSE, residentAdmissionRetry));
	const SbBool representationRefinementAllowed =
	    (scaleDemandChanged || residentAdmissionRetry) &&
	    !this->d->lodAvailabilityLedger.residencyDrainActive() &&
	    !this->d->lodPresentationTransaction.barrierPending() ? TRUE : FALSE;
	action.setAllowRepresentationRefinement(
	    representationRefinementAllowed);
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
	    viewLodState && this->getActiveLodMeshPayloadCount() == 0 &&
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
		this->d->lodRetainedAdmissionQualityEvidence.
		    maximumNormalizedError());
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
	 * contract.  A quiet view may spend a bounded 8 ms slice filling
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
	    residentAdmissionRetry &&
	    !this->d->lodSubmissionDelta.active() && viewLodState;
	const bool currentViewSourceCensusComplete =
	    this->d->hasCompleteLodConvergenceCandidateCensus(source);
	const bool currentViewIsExactlyEmpty =
	    currentViewSourceCensusComplete &&
	    this->d->lodConvergenceCandidateCount() == 0;
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
	    (boundedLargeCompact || usingMemoryAdmissionFrontier) &&
	    !this->d->lodViewDemandPolicy.demandRefreshActive() &&
	    !this->d->lodCoveragePolicy.demandCensusRequired() &&
	    !this->d->lodSubmissionDelta.active() &&
	    !this->d->lodStructuralRepair.active() &&
	    !this->d->lodCoveragePolicy.active() &&
	    this->d->lodCoveragePolicy.coverageComplete() &&
	    !this->d->lodAdmissionCursor.retainedAdmission() &&
	    viewLodState &&
	    (viewLodState->hasCadPresentationAssemblies() ||
	     currentViewIsExactlyEmpty) &&
	    sourceMeshRequests > 0 &&
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
	bool projectionRefreshPlan = false;
	if (applyRetainedAdmission &&
	    !this->d->lodSubmissionSourcePlan.validFor(source)) {
	    const BObolRetainedAllocationResult &retainedAllocation =
		this->d->lodRetainedAllocationCertificate;
	    for (const BObolRetainedProjectionRefreshPlan &refresh :
		    retainedAllocation.projectionRefreshPlans) {
		if (refresh.source != source)
		    continue;
		if (!refresh.denseRefreshRequired) {
		    this->d->lodSubmissionSourcePlan.begin(source);
		    std::vector<size_t> &entries =
			this->d->lodSubmissionSourcePlan.entries();
		    entries.reserve(refresh.compactEntryIndices.size());
		    for (uint32_t entry : refresh.compactEntryIndices)
			entries.push_back(static_cast<size_t>(entry));
		    projectionRefreshPlan = true;
		}
		break;
	    }
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
	    if (usingMemoryAdmissionFrontier &&
		getenv("BOBOL_LOD_TRACE_COMPACTION"))
		bu_log("BObol resident admission frontier source=%s "
		       "revision=%llu keys=%zu entries=%zu\n",
		       source->path.getValue().getString(),
		       static_cast<unsigned long long>(admissionRevision),
		       unsatisfiedKeys.size(),
		       this->d->lodSubmissionSourcePlan.entries().size());
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
	    if (!this->d->lodViewDemandPolicy.demandRefreshActive())
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
	} else if (boundedLargeCompact && !applyRetainedAdmission &&
	    !usingUnsatisfiedFrontier &&
	    (!this->d->lodSubmissionDelta.active() || !selectiveDeltaPlan) &&
	    this->d->lodSubmissionSourcePlan.validFor(source) &&
	    this->d->lodSubmissionSourcePlan.size() <
		static_cast<size_t>(compactCount)) {
	    /* Structural streaming appends leaves while a dense pinned scan is in
	     * flight.  Extend that plan tail without restarting the consumed
	     * prefix; restarting at zero on every batch starves late leaves.  A
	     * retained-allocation plan is intentionally sparse and complete: its
	     * omitted point-presented and already-matching occurrences must never
	     * be reintroduced as a synthetic dense tail. */
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
	if (this->d->lodSubmissionEntryOffset == 0 &&
	    controller_lod_trace_enabled("BOBOL_LOD_TRACE_PASS",
		this->d->lodViewRevision.value()))
	    bu_log("BObol LoD pass start path=%s retained=%d intent=%d "
		   "capacity_allocation=%d point_recovery=%d "
		   "unsatisfied=%d delta=%d selective=%d plan=%zu\n",
		   source->path.getValue().getString(),
		   retainedAllocationPass ? 1 : 0,
		   this->d->lodSubmissionIntent.retainedAdmission() ? 1 : 0,
		   this->d->lodAdmissionEvidence.capacity().
		       capacityAllocationPending() ? 1 : 0,
		   this->d->lodPointQualityPhase.triangleRecoveryPending() ? 1 : 0,
		   usingUnsatisfiedFrontier ? 1 : 0,
		   this->d->lodSubmissionDelta.active() ? 1 : 0,
		   selectiveDeltaPlan ? 1 : 0,
		   this->d->lodSubmissionSourcePlan.validFor(source) ?
		       this->d->lodSubmissionSourcePlan.size() : 0);
	action.setCompactEntryRange(this->d->lodSubmissionEntryOffset,
	    boundedLargeCompact ? compactWave : SIZE_MAX);
	if (this->d->lodForcedCut)
	    action.setForcedCut(*this->d->lodForcedCut);
	/* Only a dense all-entry projection may replace the authoritative
	 * visibility census.  Retained allocation deliberately uses a sparse plan
	 * containing only occurrences whose cuts differ from its certificate.  If
	 * that plan resets and then completes the dense census, every omitted
	 * already-matching or point-presented occurrence becomes falsely invisible;
	 * an empty sparse plan can consequently prove an occupied view empty. */
	const bool denseViewCandidatePass =
	    !retainedAllocationPass &&
	    !projectionRefreshPlan &&
	    !usingUnsatisfiedFrontier &&
	    !this->d->lodSubmissionDelta.active();
	if (denseViewCandidatePass && source->hasCompactInstanceIndex() &&
	    this->d->lodSubmissionEntryOffset == 0)
	    this->d->beginLodConvergenceCandidateCensus(source,
		static_cast<size_t>(std::max(0, compactCount)));
	const int64_t sourceActionStarted = bu_gettime();
	action.apply(source);
	if (action.getDeferredMeshDemandCount() > 0)
	    this->d->lodCoveragePolicy.noteDemandDeferred();
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
	if (denseViewCandidatePass) {
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
	if (denseViewCandidatePass &&
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
	if (usingMemoryAdmissionFrontier &&
	    getenv("BOBOL_LOD_TRACE_COMPACTION"))
	    bu_log("BObol resident admission action source=%s visited=%u "
		   "skipped=%u tasks=%u cuts=%u budget_used=%zu "
		   "budget_blocked=%u pending_refinement=%u "
		   "allow_refinement=%d allow_representation=%d\n",
		   source->path.getValue().getString(),
		   action.getVisitedMeshCount(), action.getSkippedMeshCount(),
		   action.getSubmittedTaskCount(), action.getUpdatedCutCount(),
		   action.getRefinementCostBudgetUsed(),
		   action.getRefinementBudgetBlockedCount(),
		   action.getPendingRetainedRefinementCount(),
		   retainedRefinementAllowed ? 1 : 0,
		   representationRefinementAllowed ? 1 : 0);
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
	this->d->positionLodSubmissionCursor(++i);
    }

    const bool completedPass =
	this->d->lodSubmissionSourceIndex >= sources.size();
    const bool completedSelectivePass = completedPass &&
	this->d->lodSubmissionDelta.active();

    /* A physical-demand policy change cannot use the sparse frontier built
     * for its predecessor epoch: payloads which satisfied the old target are
     * precisely the entries which may need a richer cut now.  The gate above
     * therefore forces one dense, bounded pass.  Retire the level-triggered
     * obligation only after that pass has consumed every source; a selective,
     * structural-repair, or retained-allocation pass is not such a proof. */
    const bool completedViewDemandRefresh =
	this->d->lodViewDemandPolicy.completeDemandRefresh(
	completedPass,
	completedSelectivePass,
	this->d->lodSubmissionPass.rescanPending(),
	this->d->lodStructuralRepair.active(),
	this->d->lodSubmissionIntent.retainedAdmission());
    if (completedPass)
	this->d->lodPlanningObligations.retireResidentAdmissionRetry();
    const bool completedStructuralRepair =
	completedPass && this->d->lodStructuralRepair.active();
    const bool completedTerminalProxyRepair = completedStructuralRepair &&
	this->d->lodStructuralRepair.terminalProxy();
    const bool completedPointRelaxation = completedStructuralRepair &&
	this->d->lodStructuralRepair.pointRelaxationPending();
    const size_t completedMissingMeshBudgetBlockedCount = completedPass ?
	this->d->lodRetainedPass.missingMeshBudgetBlockedCount() : 0;
    const bool completedPassRetainedAllocation =
	completedPass &&
	this->d->lodSubmissionIntent.retainedAdmission();
    const bool completedPassChangedCut = completedPass &&
	this->d->lodRetainedPass.cutAdvanced();
    /* Point-threshold relaxation is a private preload transaction.  Its
     * existing cut may hit the current scene allowance while the old coherent
     * point presentation remains on screen; that annotation belongs to the
     * candidate's later framebuffer, not to the current capacity certificate.
     * Letting it enter CALIBRATE_CAPACITY resets a terminal static proof
     * before the candidate threshold is committed and alternates the same
     * coarse/fine populations indefinitely. */
    const bool completedPassBudgetBlocked = completedPass &&
	!completedPointRelaxation && this->d->lodRetainedPass.budgetBlocked();
    const bool completedPassRefinementPending = completedPass &&
	!completedPointRelaxation &&
	this->d->lodRetainedPass.refinementPending();
    const bool completedCapacityAllocation = completedPass &&
	this->d->lodAdmissionEvidence.capacity().capacityAllocationPending();
    const bool residentWorkPending = completedPass && service &&
	(service->activeTaskCountForGeneration(generation) > 0 ||
	 service->queuedResultCountForGeneration(generation) > 0 ||
	 this->d->lodAvailabilityLedger.resultsPending() > 0);
    if (completedPointRelaxation) {
	if (completedMissingMeshBudgetBlockedCount > 0) {
	    /* The current point cut remains a coherent, fully presented image.
	     * If the sparse preload cannot fit its finite scene allowance, retain
	     * that image and record the capacity boundary rather than exposing
	     * boxes or replaying the same unsuccessful preload forever. */
	    const BObolLodAdmissionPlan structuralLimitPlan =
		BObolLodAdmissionPlanner::settlePointAtStructuralLimit(
		    this->d->lodAdmissionEvidence,
		    this->d->lodAdmissionCursor,
		    BObolLodPointCalibrationGoal::RESPONSIVE,
		    this->d->lodPresentationPointProxyPixelThreshold);
	    this->d->commitAdmissionPlan(structuralLimitPlan);
	    this->d->lodStructuralRepair.reset();
	} else {
	    /* Keep the candidate threshold private until an exact presentation
	     * proves that every structural occurrence it would reveal has been
	     * replaced.  Worker results may arrive after this submission pass, so
	     * the completed framebuffer—not task completion—is the commit edge. */
	    this->d->lodStructuralRepair.completePointRelaxationAdmission();
	    this->requestLodCapacityRender("lod-point-prefetch-publication");
	}
    } else if (completedStructuralRepair) {
	if (completedMissingMeshBudgetBlockedCount == 0 &&
	    !completedTerminalProxyRepair) {
	    /* Structural first-wave thresholds bound provider admission; they are
	     * not renderer-unsafe timing samples.  Once every selected fallback has
	     * acquired a mesh, discard that startup bracket so completed mesh
	     * frames alone decide how far the point population may relax.  A
	     * terminal-proxy pass preserves its point cut: resetting it would expose
	     * the smaller structural occurrences which that same capacity witness
	     * deliberately left in the aggregate channel. */
	    this->d->applyAdmissionEvidenceAction(
		BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
	}
	this->d->lodStructuralRepair.reset();
    }
    if (completedPass)
	this->d->lodRetainedPass.clearMissingMeshBudgetBlocked();
    if (completedPointRelaxation)
	/* Private preload annotations describe a hidden candidate, not the
	 * displayed point population.  Its publication/frame transaction owns the
	 * only successor after this mechanical pass. */
	this->d->resetRetainedPassAnnotations();
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
		completedPassBudgetBlocked || completedCapacityAllocation);
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
		   completedSelectivePass ? 1 : 0,
		   this->d->lodSubmissionPass.rescanPending() ? 1 : 0,
		   this->d->lodPresentationPolicy.handoffPending() ? 1 : 0,
		   this->d->lodRetainedPass.refinementPending() ? 1 : 0,
		   this->d->lodRetainedPass.residencyPending() ? 1 : 0,
		   completedStructuralRepair ? 1 : 0,
		   this->d->lodRetainedPass.cutAdvanced() ? 1 : 0,
		   this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ? 1 : 0,
		   sceneActiveFaces, sceneActiveCost,
		   this->d->lodAdmissionEvidence.capacity().currentBudget(),
		   this->d->lastLodVisitedMeshCount,
		   this->d->lastLodUpdatedCutCount);
    }
    const BObolLodCoveragePolicy::Completion coverageCompletion =
	this->d->lodCoveragePolicy.completeIfReady(
	    completedPass &&
		!this->d->lodRetainedViewContinuity.visibilityCensusDeferred(),
	    this->d->lodSubmissionPass.rescanPending());
    const bool completedDemandCensus =
	completedPass && !coverageCompletion.completed &&
	!completedPassRetainedAllocation && !completedStructuralRepair &&
	!completedSelectivePass &&
	!this->d->lodSubmissionPass.rescanPending() &&
	this->d->lodCoveragePolicy.completeDemandCensus();
    if (completedDemandCensus && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_PASS", this->d->lodViewRevision.value()))
	bu_log("BObol LoD completed dense demand census view=%llu "
	       "policy=%llu\n",
	       static_cast<unsigned long long>(
		   this->d->lodViewRevision.value()),
	       static_cast<unsigned long long>(
		   this->d->lodPolicyRevision.value()));
    const bool retainedImportanceCensusCompleted =
	this->d->lodPlanningObligations.importanceCensusPending() &&
	((coverageCompletion.completed && !coverageCompletion.missing) ||
	 completedViewDemandRefresh);
    if (completedSelectivePass) {
	/* The selective target set is the cursor's mode, not a successor-frame
	 * obligation.  Once that cursor has consumed every target, later exact
	 * presentation or capacity work owns its own state.  Keeping the target
	 * set alive makes an ordinary demand refresh appear selective and prevents
	 * its required scene-wide pass from ever starting. */
	this->d->lodSubmissionDelta.reset();
    }
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
	 * Only the bounded case needs the extra capacity/presentation wave below.
	 *
	 * Coverage completion must not retire physical-demand refresh.  The same
	 * cursor may have run as retained allocation or structural repair, neither
	 * of which proves a current-view importance census.  completeDemandRefresh()
	 * above is the sole retirement edge for that level-triggered obligation;
	 * clearing it here left the importance census producerless until another
	 * pump happened to recreate the ordinary pass. */
	if (coverageCompletion.missing) {
	    this->d->publishCadCameraMotionFrameReuse(FALSE);
	}
    }
    bool capacityOwnsCompletedPass = false;
    bool capacitySuccessorDeferredByHandoff = false;
    BObolLodCapacityEvidence::CalibrationInputs deferredCalibrationInputs;
    const auto applyCapacityCalibrationDecision =
	[this](const BObolLodCapacityEvidence::CalibrationDecision &decision) {
	    if (decision.searchTerminal && decision.searchResult ==
		    BObolLodCapacitySearchCertificate::Result::CERTIFIED) {
		/* A bounded pass may consume the final exact sample instead of the
		 * render-completion path.  Retire the same superseded handoff at
		 * either refinement boundary. */
		this->d->lodPresentationPolicy.acceptCapacityCertificate();
	    }
	    if (decision.restartSubmission) {
		if (decision.candidateReallocation)
		    (void)this->d->lodAvailabilityLedger.
			transferResidentGrowthToCapacityCandidate();
		/* A capacity candidate starts a distinct pass.  No completed-pass
		 * annotation may cross this ownership boundary. */
		this->beginSceneWideCapacitySubmission();
		this->markProgressiveWorkPending();
	    }
	    if (decision.requestFrame) {
		const char *frameReason = decision.sampleFrame ?
		    "lod-budget-sample" : "lod-population-barrier";
		this->requestLodCapacityRender(frameReason);
	    }
	};
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodRetainedPass.clearAdmittedWork();
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::CLEAR_CAPACITY_LIMIT);
	this->d->retireRetainedRefinementObservation();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
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
	    this->d->lodSubmissionPass.retire();
	    if (residentGrowthOwnsSuccessor) {
		if (completedPassSuccessor ==
			BObolLodAvailabilityScheduler::CompletedPassSuccessor::
			    COMPLETE_RESIDENCY_DRAIN)
		    (void)this->d->lodAvailabilityLedger.
			completeResidencyDrain();
		this->markProgressiveWorkPending();
	    } else {
		this->d->requestCapacityRescan();
		this->requestLodCapacityRender("lod-coverage-admission");
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
	    this->d->publishCadProgressiveCeiling(-1);
	    this->d->lodSubmissionPass.retire();
	    this->d->requestCapacityRescan();
	    this->requestLodCapacityRender("lod-coverage-minimum-calibration");
	} else {
	    /* Every projected leaf has a useful structural or mesh presentation. Begin a
	     * fresh bounded pass which may now spend the remaining scene budget
	     * on screen-value-ordered PoP refinement. */
	    this->d->lodSubmissionPass.beginFresh();
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh(!sources.empty());
	this->d->resetRetainedPassAnnotations();
	this->d->advanceAdmissionRevision(
	    BObolLodAdmissionRevisionDomain::CAPACITY);
	if (this->d->lodSubmissionPass.active())
	    this->markProgressiveWorkPending();
    } else if (retainedImportanceCensusCompleted) {
	/* Small sources do not take the bounded-coverage successor branch above,
	 * so explicitly start their one-shot importance allocation here. */
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh(!sources.empty());
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
	 * Leaving delta mode armed reconstructs the old selective plan after the
	 * source-local cursor state is cleared and can loop over that prefix
	 * forever while newly appended leaves remain boxes.
	 */
	const bool compactInventoryIncomplete =
	    controller_lod_compact_inventory_incomplete(sources);
	this->d->rewindLodSubmissionCursor();
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
	    this->d->lodSubmissionPass.beginFresh(!sources.empty());
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
	 * larger calibrated allowance unused.
	 *
	 * This branch also applies the final cut selected by a terminal bounded
	 * capacity certificate.  Retire only the mechanical measurement in that
	 * case.  Clearing the terminal evidence and its budget-limited witness
	 * makes the next blocked pass restart the identical search, visibly
	 * alternating between the same coarse and rich Lucy populations. */
	const BObolRetainedAllocationResult &completedAllocation =
	    this->d->lodRetainedAllocationCertificate;
	BObolLodCapacityEvidence::CompletedAllocationInputs completion;
	completion.revisionStamp = this->d->admissionRevisionStamp();
	completion.requestedSceneBudget =
	    completedAllocation.requestedSceneBudget;
	completion.certifiedPresentationBudget =
	    completedAllocation.certifiedPresentationBudget;
	completion.allocationCertificateCurrent =
	    completedPassRetainedAllocation &&
	    this->d->retainedAllocationCertificateCurrent(sceneLodState);
	completion.allocationCutsApplied =
	    this->d->retainedAllocationCutsApplied(sceneLodState);
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.retire();
	this->d->lodRetainedPass.clearAdmittedWork();
	this->d->completeAppliedAllocation(completion);
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.retire();
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh(
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.retire();
	this->d->resetRetainedPassAnnotations();
	this->markProgressiveWorkPending();
	this->requestLodCapacityRender("lod-point-calibration-successor");
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
	/* A ready handoff begins with a coherent retained presentation.  Pose-only
	 * orthographic changes can preserve its exact occurrence cuts; zoom still
	 * needs a new pixel-demand allocation.  Neither case needs to discard that
	 * useful presentation merely to rediscover the preferred redraw cadence.
	 * Start its bounded search at the event-driven static deadline.  A measured
	 * miss may still coarsen it, while a safe frame establishes a lower bound
	 * from which only richer candidates may be explored. */
	const bool retainedStaticPresentation =
	    this->d->lodRetainedViewContinuity.startCapacityAtStatic();
	const BObolLodCapacitySearchCertificate &capacitySearch =
	    this->d->lodAdmissionEvidence.capacity().capacitySearch();
	const bool capacitySearchActive = capacitySearch.awaitingSample();
	const uint64_t capacitySearchPreferredNanoseconds =
	    capacitySearchActive ?
		capacitySearch.key().preferredTargetNanoseconds :
	    (this->d->lodStaticQualityTrial.usesStaticDeadline() ?
		this->d->prominentQualityPresentationDeadline() :
		preferredStableTargetNanoseconds);
	/* Stable pass decisions require a duration for this retained cut.  GPU
	 * query latency is handled by the paired throughput calibration above;
	 * reusing that duration after the occurrence population changes is not a
	 * current-cut deadline measurement. */
	const uint64_t observedStableNanoseconds =
	    this->d->lodRendererPerformanceEvidence.cadPresentationNanosecondsAt(
		this->d->lodPresentationPointProxyPixelThreshold);
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
	    capacityAllocation.revisionStamp.same(
		this->d->admissionRevisionStamp());
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
		capacityAllocation.viewRevision(),
		capacityAllocation.policyRevision(),
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
	size_t completedSubpixelOccurrences = 0;
	size_t completedStructuralBoxes = 0;
	const bool completedOccurrenceCoverageExact = sceneLodState &&
	    sceneLodState->lastCadPresentationFrameExact() &&
	    sceneLodState->lastCadPresentationOccurrenceCoverage(
		completedSubpixelOccurrences, completedStructuralBoxes);
	const size_t completedTerminalFailures = sceneLodState ?
	    sceneLodState->cadOccurrenceTerminalFailureCount() : 0;
	/* The projection histogram includes structural source records already
	 * represented by the exact frame's subpixel aggregate.  Those records are
	 * valid terminal presentations at the current view.  Only uncollapsed
	 * fallbacks in the framebuffer belong to structural repair. */
	const bool structuralFrontierPending =
	    BObolLodPresentationPolicy::nonterminalStructuralFrontier(
		completedOccurrenceCoverageExact, completedStructuralBoxes,
		completedTerminalFailures);
	const bool overBudgetAllocationClaimed =
	    !structuralFrontierPending &&
	    !this->d->lodStaticQualityTrial.capacityRejected() &&
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
	capacityHandoffInputs.capacityTransactionPending =
	    this->d->lodAdmissionEvidence.capacity().capacityTransactionPending();
	capacityHandoffInputs.structuralFrontierPending =
	    structuralFrontierPending;
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
	const bool pendingCapacitySample = BObolLodPresentationPolicy::
	    capacitySamplePending(
		this->d->lodAdmissionEvidence.capacity().capacitySearch().
		    awaitingSample(),
		currentCapacityAllocation, allocationPresentationRealized);
	const BObolLodPresentationPolicy::CompletedPassSelection
	    capacityPassSelection =
	this->d->lodPresentationPolicy.completedPassSelection(
	    capacityHandoffInputs, !overBudgetAllocationClaimed,
	    pendingCapacitySample,
	    pendingCapacitySample && !allocationPresentationRealized,
	    this->d->lodInteractiveProgressiveCeiling);
	if (capacityPassSelection.structuralFrontierOwns()) {
	    /* Any capacity sample describes a mesh/point occurrence population.
	     * Structural-only records make that population inexact, so retire only
	     * its measurement transaction.  The current budget remains available
	     * to the finite structural successor selected below. */
	    this->d->applyAdmissionEvidenceAction(
		BObolLodAdmissionPlanner::EvidenceAction::
		    INVALIDATE_CAPACITY_MEASUREMENT);
	}
	capacityOwnsCompletedPass = capacityPassSelection.capacityOwns();
	const BObolLodPresentationPolicy::CompletedPassOwner completedPassOwner =
	    capacityPassSelection.owner;
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
	calibrationInputs.populationDigest = currentCapacityAllocation ?
	    capacityAllocation.selectedPopulationDigest : 0;
	calibrationInputs.populationIdentity = currentCapacityAllocation ?
	    capacityAllocation.selectedPopulationIdentity :
	    0;
	calibrationInputs.populationMinimumBudget = currentCapacityAllocation ?
	    capacityAllocation.selectedPresentationCost : 0;
	calibrationInputs.nextDistinctPopulationBudget =
	    currentCapacityAllocation ?
		capacityAllocation.nextDistinctPresentationBudget : 0;
	calibrationInputs.nextDistinctPopulationBudgetKnown =
	    currentCapacityAllocation &&
	    capacityAllocation.nextDistinctPresentationBudgetKnown;
	calibrationInputs.maximumCandidateBudget = currentCapacityAllocation ?
	    (capacitySearchActive ?
		capacitySearch.key().maximumCandidateBudget :
		capacityAllocation.maximumCapacitySearchBudget()) : 0;
	const uint64_t prospectiveMaximumTargetNanoseconds =
	    currentCapacityAllocation &&
		!this->d->lodAdmissionEvidence.resources().anyPressure() ?
	    this->d->prominentQualityPresentationDeadline() :
	    capacitySearchPreferredNanoseconds;
	const uint64_t capacitySearchInitialTargetNanoseconds =
	    capacitySearchActive ? capacitySearch.activeTargetNanoseconds() :
	    (retainedStaticPresentation ? prospectiveMaximumTargetNanoseconds :
		capacitySearchPreferredNanoseconds);
	/* The lower bracket belongs to the certificate's first target.  Seed only
	 * from this exact completed frame when its measured duration proves that
	 * target directly; evidence from a different deadline is not transferable. */
	calibrationInputs.knownSafeBudget =
	    allocationPresentationRealized && observedStableNanoseconds > 0 &&
	    observedStableNanoseconds <= capacitySearchInitialTargetNanoseconds ?
		calibrationInputs.activeCost : 0;
	if (capacitySearchActive) {
	    /* Candidate budgets are coordinates in one frozen search domain.  The
	     * current retained allocation was built with that candidate as its
	     * marginal allowance, so deriving a new domain from the allocation
	     * would shrink maximumCandidateBudget after every miss.  That rekeys
	     * the certificate, loses its steady bracket, and can certify before the
	     * independently bounded static goal is ever attempted. */
	    calibrationInputs.searchKey = capacitySearch.key();
	} else {
	    const size_t capacityMinimumBudget =
		this->d->lodAdmissionEvidence.capacity().
		    capacitySearchMinimumBudget(sceneMinimumActiveCost,
			capacityAllocation.protectedFloorBudget,
			capacityAllocation.protectedFloorIdentity);
	    calibrationInputs.searchKey =
		BObolLodCapacitySearchCertificate::keyFor(
		    this->d->admissionRevisionStamp(),
		    capacitySearchPreferredNanoseconds,
		    prospectiveMaximumTargetNanoseconds,
			this->d->lodPresentationPointProxyPixelThreshold,
			calibrationInputs.maximumCandidateBudget,
			capacityMinimumBudget, retainedStaticPresentation);
	    calibrationInputs.searchKey.preferredBudgetCeiling =
		this->d->lodAdmissionEvidence.capacity().
		    deadlineCapacityCeiling();
	    calibrationInputs.searchKey.maximumBudgetCeiling =
		this->d->lodAdmissionEvidence.capacity().
		    staticDeadlineCapacityCeiling();
	}
	capacitySuccessorDeferredByHandoff =
	    capacityPassSelection.deferredCapacitySuccessor;
	if (capacitySuccessorDeferredByHandoff)
	    deferredCalibrationInputs = calibrationInputs;
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
	    completedStructuralRepair && !completedPointRelaxation ?
		completedMissingMeshBudgetBlockedCount : 0;
	const float structuralPointThresholdBefore =
	    this->d->lodPresentationPointProxyPixelThreshold;
	if (capacityPassSelection.structuralFrontierOwns() &&
	    structuralTailBlockedCount > 0 &&
	    !this->d->lodInteractionSession.active() &&
	    /* Structural coverage is discovered before the occurrence batch is
	     * necessarily installed.  Point aggregation only changes that batch;
	     * arming its completed-frame barrier before one exists blocks the very
	     * source admission which creates it.  Hidden-line and other direct
	     * presentation paths exposed this as an infinite "calibrating
	     * small-part aggregation" loop containing only structural boxes. */
	    sceneLodState && sceneLodState->hasCadPresentationAssemblies() &&
	    adaptivePointAggregationAllowed &&
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
		this->d->publishCadPointProxyThreshold(
		    structuralAggregation.threshold);
		this->requestStructuralPointAdmissionFrame(
		    "lod-structural-tail-aggregation");
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
		       "retained=%d retained_static=%d demanded=%zu "
		       "preferred_ms=%.3f maximum_ms=%.3f "
		       "candidate_reallocation=%d search_active=%d "
		       "sample_frame=%d restart=%d result=%u phase=%u goal=%u "
		       "samples_remaining=%u candidates=%u total_candidates=%u "
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
		       retainedStaticPresentation ? 1 : 0,
		       calibrationInputs.maximumCandidateBudget,
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
		       this->d->lodAdmissionEvidence.capacity().capacitySearch().
			   totalMeasuredCandidateCount(),
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.retire();
	applyCapacityCalibrationDecision(calibration);
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
	    !this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() &&
	    this->d->lodAdmissionEvidence.capacity().stableBudgetLimited() &&
	    !this->d->lodInteractionSession.active() &&
	    adaptivePointAggregationAllowed &&
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
		    BObolLodPointCalibrationGoal::RESPONSIVE,
		    this->d->lodPresentationPointProxyPixelThreshold,
		    static_cast<uint64_t>(observedStableNanoseconds),
		    this->d->pointCalibrationTargetFps(
			BObolLodPointCalibrationGoal::RESPONSIVE));
	    this->d->commitAdmissionPlan(pressurePlan);
	    const BObolLodPointProxyEvidence::Decision &pressure =
		pressurePlan.pointProxyDecision;
	    if (pressure.changed) {
		this->d->publishCadPointProxyThreshold(pressure.threshold);
		this->d->lodPointQualityPhase.requestCalibration();
		/* The multiplicative pressure correction intentionally lands on
		 * the safe side of the target.  Its next unchanged frame continues
		 * the bounded bracket search, so terminal quality is the finest cut
		 * which meets the stable FPS contract. */
		this->requestLodCapacityRender("lod-stable-point-calibration");
	    }
	}
	/* Capacity progress is typed control state, not a diagnostic.  Publishing
	 * it here left the final report saying "measuring" after a later frame had
	 * certified the search.  QgProgressiveDiagnostics derives the current HUD
	 * label from the convergence snapshot instead. */
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
	/* Completing an ordinary or handoff-owned mechanical pass is not a new
	 * capacity problem.  In particular, a presentation reconciliation may
	 * finish immediately after a bounded search reaches TERMINAL; clearing the
	 * limit here leaves the terminal certificate in place but lets the generic
	 * throughput planner exceed it.  Only the explicit view, policy,
	 * inventory, availability, or coverage transitions invalidate capacity
	 * evidence. */
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
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.retire();
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
	this->requestLodCapacityRender("lod-static-overscan-allocation");
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
    handoffInputs.capacityTransactionPending =
	this->d->lodAdmissionEvidence.capacity().capacityTransactionPending();
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
    size_t handoffSubpixelOccurrences = 0;
    size_t handoffStructuralBoxes = 0;
    const bool handoffOccurrenceCoverageExact = sceneLodState &&
	sceneLodState->lastCadPresentationFrameExact() &&
	sceneLodState->lastCadPresentationOccurrenceCoverage(
	    handoffSubpixelOccurrences, handoffStructuralBoxes);
    const size_t handoffTerminalFailures = sceneLodState ?
	sceneLodState->cadOccurrenceTerminalFailureCount() : 0;
    handoffInputs.structuralFrontierPending =
	BObolLodPresentationPolicy::nonterminalStructuralFrontier(
	    handoffOccurrenceCoverageExact, handoffStructuralBoxes,
	    handoffTerminalFailures);
    BObolLodPresentationPolicy::CompletedPassSelection handoffSelection;
    if (!capacityOwnsCompletedPass) {
	const bool handoffCapacitySamplePending =
	    BObolLodPresentationPolicy::capacitySamplePending(
		this->d->lodAdmissionEvidence.capacity().capacitySearch().
		    awaitingSample(),
		allocationCertificateCurrent && handoffAllocationCutsApplied,
		this->d->retainedAllocationPresentationRealized(sceneLodState));
	handoffSelection =
	    this->d->lodPresentationPolicy.completedPassSelection(
		handoffInputs, false,
		handoffCapacitySamplePending,
		handoffCapacitySamplePending &&
		    !this->d->retainedAllocationPresentationRealized(
			sceneLodState),
		this->d->lodInteractiveProgressiveCeiling);
	if (handoffSelection.consumePassAnnotations) {
	    handoffInputs.capacityTransactionPending = false;
	    handoffInputs.changedCut = false;
	}
    }
    BObolLodPresentationPolicy::CompletedPassDecision handoff;
    /* A capacity-owned completion may retire unrelated observation data, but
     * it may not run the handoff reducer after changing capacity evidence.
     * Selection already normalized the one case where the handoff must
     * precede a pending sample.  Structural ownership still runs completePass:
     * that operation consumes the ordinary pass's richer-cut observation so
     * armStableLodHeadroomProbeIfReady() can start the finite structural
     * successor below. */
    if (!capacityOwnsCompletedPass)
	handoff = this->d->lodPresentationPolicy.completePass(handoffInputs);
    if (completedPass && controller_lod_trace_enabled(
	    "BOBOL_LOD_TRACE_PASS", this->d->lodViewRevision.value())) {
	static std::atomic<unsigned int> handoffTraceCount(0);
	const unsigned int traceIndex = handoffTraceCount.fetch_add(1);
	if (traceIndex < 512)
	    bu_log("BObol LoD handoff pass owner=%u consume_annotations=%d "
		   "pending=%d rescan=%d changed=%d "
	       "reconciled=%d structural=%d "
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
	       handoffInputs.capacityTransactionPending ? 1 : 0,
	       handoffInputs.changedCut ? 1 : 0,
	       handoffInputs.presentationLimitsReconciled ? 1 : 0,
	       handoffInputs.structuralFrontierPending ? 1 : 0,
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
    const bool handoffNeedsIndependentResidentGrowth =
	BObolLodAvailabilityScheduler::needsIndependentResidentGrowth(
	    handoffResidentPopulationPending,
	    this->d->lodAdmissionEvidence.capacity().capacitySearch().
		awaitingSample());
    if (completedPass && completedPassRetainedAllocation &&
	this->d->lodPresentationPolicy.handoffPending() &&
	handoffNeedsIndependentResidentGrowth) {
	/* Coalesce every late completion/publication into one new allocation
	 * after the stream becomes idle.  This is a liveness witness as well as
	 * a debounce: without it the handoff can remain armed after the last
	 * worker result, while immediately restarting here wastes an O(scene)
	 * pass for every result batch.  A capacity candidate's suffix producer and
	 * a retained cut's own presentation barrier are not resident growth.
	 * Treating either as a separate population edge
	 * scheduled another allocation after every coherent frame and replayed a
	 * tiny set of cuts forever. */
	this->d->lodAvailabilityLedger.noteRicherPrefixAvailable();
	this->markProgressiveWorkPending();
    }

    const auto beginSceneWideHandoffAllocation =
	[this, &sources](size_t reconciliationBudget) {
	    if (reconciliationBudget > 0)
		this->d->requestPresentationReconciliation(
		    reconciliationBudget);
	    else
		this->d->requestRetainedReallocation();

	    this->beginSceneWideCapacitySubmission(
		sources.empty() ? FALSE : TRUE);
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
	    if (this->d->lodSubmissionPass.active())
		this->markProgressiveWorkPending();
	};

    if (handoff.requestRetainedAllocation) {
	/* A targeted source/result delta may be the last pass after the worker
	 * stream becomes quiet.  It cannot release the presentation handoff, but
	 * no external edge remains to request the authoritative allocation.  Start
	 * that transaction now and keep the existing constrained framebuffer on
	 * screen while its bounded owner-thread slices advance. */
	const size_t reconciliationBudget =
	    this->d->lodPresentationPolicy.handoffReconciliationBudget();
	beginSceneWideHandoffAllocation(reconciliationBudget);
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
	if (adaptivePointAggregationAllowed &&
	    this->d->retainedPointProxyAggregationApplicable()) {
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
	    this->d->publishCadPointProxyThreshold(reduction.threshold);
	    this->d->resetRetainedAdmissionQualityProof();
	    const size_t retryBudget = handoffReconciliationBudget > 0 ?
		handoffReconciliationBudget :
		allocationCertificate.certifiedPresentationBudget;
	    beginSceneWideHandoffAllocation(retryBudget);
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
	    const uint64_t currentCadSampleNanoseconds =
		this->d->lodRendererPerformanceEvidence.
		    cadPresentationNanosecondsAt(
			this->d->lodPresentationPointProxyPixelThreshold);
	    const size_t staticBudget = BObolLodQualityPolicy::
		staticLocalMinimumRetryBudget(
		    presentedCost,
		    allocationCertificate.selectedPresentationCost,
		    allocationCertificate.pixelDemandPresentationCost,
		    currentCadSampleNanoseconds,
		    this->d->prominentQualityPresentationDeadline());
	    if (!this->d->lodStaticQualityTrial.blocksNewTrial() &&
		staticBudget > handoffReconciliationBudget) {
		this->d->lodStaticQualityTrial.begin(
		    this->d->lodPresentationPointProxyPixelThreshold);
		this->d->lodPresentationPolicy.armHandoff(
		    false, presentedCost, staticBudget);
		beginSceneWideHandoffAllocation(staticBudget);
		if (controller_lod_trace_enabled("BOBOL_LOD_TRACE_HEADROOM",
			this->d->lodViewRevision.value()))
		    bu_log("BObol LoD static local-minimum trial "
			   "presented=%zu selected=%zu budget=%zu "
			   "elapsed_ms=%.3f deadline_ms=%.3f\n",
			   presentedCost,
			   allocationCertificate.selectedPresentationCost,
			   staticBudget,
			   currentCadSampleNanoseconds / 1000000.0,
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
		(void)this->d->lodStaticQualityTrial.
		    constrainPresented(constraint);
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
	const BObolLodStaticQualityTrial::HandoffCompletion staticHandoff =
	    this->d->lodStaticQualityTrial.completeOccurrenceHandoff(
		this->d->admissionRevisionStamp(),
		this->d->lodInteractiveProgressiveCeiling);
	if (staticHandoff.completedRejectedConstraint) {
	    const size_t constrainedBudget = this->d->lodStaticQualityTrial.
		constrainedPresentationBudgetFor(
		    this->d->admissionRevisionStamp());
	    this->d->acceptStaticPresentationConstraint(constrainedBudget);
	}
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
	if (staticHandoff.releaseRendererCeiling) {
	    this->d->publishCadProgressiveCeiling(-1);
	}
	const size_t progressivePayloadCount = presentationState ?
	    presentationState->cadProgressivePayloadCount() : 0;
	if (BObolLodPresentationPolicy::
		capacityProducerRequiredAfterHandoff(
		    capacitySuccessorDeferredByHandoff, handoff.finishHandoff,
		    progressivePayloadCount)) {
	    /* The handoff consumed the allocation pass because its exact sample
	     * first required a ceiling-free presentation.  Recreate that one
	     * displaced producer now; the search begins in PRESENTING and consumes
	     * the same frame requested below, so no second owner or repaint is
	     * introduced. */
	    const BObolLodCapacityEvidence::CalibrationDecision deferred =
		this->d->finishBlockedCapacityPass(deferredCalibrationInputs);
	    applyCapacityCalibrationDecision(deferred);
	    if (deferred.restartSubmission)
		this->d->advanceAdmissionRevision(
		    BObolLodAdmissionRevisionDomain::CAPACITY);
	}
	/* A rejected richer cut completes the search, but the predecessor is still
	 * a valid event-driven static frame.  Retain that phase after occurrence
	 * reconciliation so unrelated HUD, selection, or screenshot redraws use
	 * the same deadline which proved it.  Camera/resource/capacity transitions
	 * explicitly invalidate the phase elsewhere. */
	if (presentationState)
	    this->d->publishCadCameraMotionFrameReuse(FALSE);
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
	/* The handoff population was a reversible safety presentation, not the
	 * bounded search's successor allocation.  If the search still owns an
	 * ALLOCATING candidate, create that producer explicitly now.  Treating the
	 * certificate's broad "active" predicate as a render request only repaints
	 * the safe population forever without applying the candidate. */
	(void)this->scheduleCapacityCandidateAllocationIfReady();
	this->requestLodCapacityRender(staticHandoff.measureCeilingFreeCandidate ?
	    "lod-static-ceiling-free-candidate" : "lod-stable-handoff");
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
    const BObolRetainedAllocationResult &completedAllocation =
	this->d->lodRetainedAllocationCertificate;
    const bool pointProtectionPresentationRequired =
	completedPass && completedPassRetainedAllocation && sceneLodState &&
	completedAllocation.pointProxyProtectionChanged &&
	(!sceneLodState->cadPointProxyProtectionClassified() ||
	 !sceneLodState->lastCadPresentationFrameExact());
    if (!this->d->lodSubmissionPass.active() &&
	(this->d->lodRetainedPass.cutAdvanced() ||
	 recoveryPresentationRequired ||
	 pointProtectionPresentationRequired) &&
	(boundedScenePass || recoveryPresentationRequired ||
	 pointProtectionPresentationRequired ||
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
	/* Applying a certified occurrence allocation does not change the capacity
	 * problem which produced it.  Keep that plan current until the requested
	 * frame presents it; the frame is the measurement witness for the same
	 * candidate.  An ordinary cut mutation has no such complete-population
	 * certificate and must still invalidate prior capacity evidence. */
	if (!this->d->retainedAllocationCutsApplied(sceneLodState))
	    this->d->advanceAdmissionRevision(
		BObolLodAdmissionRevisionDomain::CAPACITY);
	this->scheduleLodRefinementFrame(
	    pointProtectionPresentationRequired ?
		"lod-point-protection" : "lod-cut");
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
	this->requestLodCapacityRender("lod-cut");
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
    BObolLodControlTransitionScope controlTransition(this);
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
    SbBool currentDemandReplayCandidate = FALSE;
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
	const int authorizedProviderCut =
	    BObolLodDeliveryPolicy::authorizedPresentationCut(
		currentAllocatedCut,
		drained[i].presentationAdmissionCertified,
		drained[i].presentationAdmissionViewRevision.value(),
		drained[i].presentationAdmissionPolicyRevision.value(),
		drained[i].presentationAdmissionCut,
		this->d->lodViewRevision.value(),
		this->d->lodPolicyRevision.value(),
		drained[i].geometry.activeCut);
	const bool allocationCoversProviderAdvance = residentCadBefore &&
	    authorizedProviderCut > residentCadBefore->activeCut;
	const int admittedProviderActiveCut = allocationCoversProviderAdvance ?
	    std::min(drained[i].geometry.activeCut,
		authorizedProviderCut) :
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
	/* A reusable immutable generation may be rebased above only after its
	 * current page set, counts, and prepared renderer generation have all been
	 * reconstructed.  Every result which still carries an obsolete compact
	 * route, population, or demand at this point is superseded.  Publishing a
	 * stale mesh as a first-view bootstrap used to bypass that rule and could
	 * make an older camera's projection visible until unrelated input woke the
	 * planner. */
	if (drained[i].request.occurrenceKey.getLength() > 0) {
	    BObolLodResultAuthenticationContext authenticationContext;
	    if (route) {
		authenticationContext.sourceRoutingId =
		    route->getCompactSourceRoutingId();
		authenticationContext.sourcePopulationEpoch =
		    route->getCompactPopulationEpoch();
	    }
	    authenticationContext.viewRevision = this->d->lodViewRevision;
	    authenticationContext.policyRevision = this->d->lodPolicyRevision;
	    const BObolLodResultAuthentication authentication =
		BObolLodResultAuthenticationContract::evaluate(
		    drained[i].request, drained[i].providerStatus,
		    authenticationContext);
	    if (!authentication.identityCurrent()) {
		drained[i].providerStatus = BOBOL_LOD_PROVIDER_SUPERSEDED;
		drained[i].stale = TRUE;
		drained[i].terminal = TRUE;
		if (drained[i].diagnostic.getLength() == 0)
		    drained[i].diagnostic =
			"LoD result identity was superseded before publication";
		if (traceResult || getenv("BOBOL_LOD_TRACE_REJECTIONS"))
		    bu_log("BObol LoD rejection reason=authentication "
			   "object=%s occurrence=%s mismatch_mask=%u\n",
			   drained[i].request.objectName.getString(),
			   drained[i].request.occurrenceKey.getString(),
			   static_cast<unsigned int>(
			       authentication.identityMismatchMask));
		directRejected++;
		currentDemandReplayCandidate = TRUE;
		append_controller_lod_diagnostic(
		    this->d->lastLodDiagnostics,
		    drained[i].request.objectPath,
		    "superseded LoD result identity rejected");
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
	const bool resultNeedsIndependentResidentGrowth =
	    BObolLodAvailabilityScheduler::needsIndependentResidentGrowth(
		resultOffersRicherResidentPrefix != FALSE,
		this->d->lodAdmissionEvidence.capacity().capacitySearch().
		    awaitingSample());
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
		const BObolViewLodState::SourceResultDisposition disposition =
		    viewState ? viewState->consumeSourceResult(
			route, drained[i]) :
		    BObolViewLodState::SourceResultDisposition::
			RETRY_CURRENT_DEMAND;
		if (disposition == BObolViewLodState::
			SourceResultDisposition::ACCEPTED) {
		    directApplied++;
		    if (retainPresentationForResidentGrowth)
			directRetainedPresentationCount++;
		    if (resultNeedsIndependentResidentGrowth)
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
		    if (disposition == BObolViewLodState::
			    SourceResultDisposition::RETRY_CURRENT_DEMAND)
			currentDemandReplayCandidate = TRUE;
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
	if (resultNeedsIndependentResidentGrowth)
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
    if (update.getCurrentDemandRetryResultCount() > 0)
	currentDemandReplayCandidate = TRUE;

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
	const size_t unsatisfiedPayloads = activePayloads > satisfiedPayloads ?
	    activePayloads - satisfiedPayloads : 0;
	const bool actionableQualityDebt =
	    unsatisfiedPayloads > memoryLimitedPayloads;
	if (actionableQualityDebt) {
	    /* Current-view quality debt is semantic, level-triggered work.  A
	     * selective result/delta, capacity candidate, or presentation barrier
	     * may temporarily own the shared cursor, but none of those mechanical
	     * transactions is allowed to erase the complete demand pass.  The
	     * demand policy retires this obligation only after an ordinary pass has
	     * consumed every current source. */
	    this->d->lodViewDemandPolicy.requestDemandRefresh();
	    this->markProgressiveWorkPending();
	}
	const bool retainedPublicationSuccessorOwned =
	    residentGrowthAdmissionCandidate ||
	    this->d->lodAdmissionEvidence.capacity().
		capacityTransactionPending();
	if (BObolLodAvailabilityScheduler::
		retainedPublicationRequiresDemandReplay(
		    retainedPresentationBatch, actionableQualityDebt,
		    retainedPublicationSuccessorOwned)) {
	    /* A retained publication which neither presents the newly resident
	     * population nor creates a resident-growth transaction has not
	     * discharged the current view demand.  Recreate its finite source
	     * cursor now; otherwise the HUD can remain in REFINING with no task,
	     * frame, or planning owner. */
	    currentDemandReplayCandidate = TRUE;
	} else if (!retainedPresentationBatch && actionableQualityDebt) {
	    partialRefinementCandidate = TRUE;
	}
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

    const auto restartCurrentDemand = [this]() {
	this->d->lodViewDemandPolicy.requestDemandRefresh();
	this->d->rewindLodSubmissionCursor();
	this->d->lodSubmissionPass.beginFresh();
	this->d->lodRetainedPass.noteRefinementPending();
	this->d->lodRetainedPass.clearAdmittedWork();
	this->markProgressiveWorkPending();
    };
    if (currentDemandReplayCandidate && this->d->lodAutoSubmit) {
	restartCurrentDemand();
	this->requestLodCapacityRender("lod-current-demand-replay");
    }

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
	    this->d->lodAdmissionEvidence.capacity().capacityAllocationPending(),
	    this->d->lodAdmissionEvidence.capacity().presentationFramePending(),
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
	    this->requestLodCapacityRender("lod-result");
	}
	this->d->lodCompactionPolicy.requestAfter(bu_gettime(), 750000);
	(void)this->enforceMeshResidencyBudget();
	if (partialRefinementCandidate && publishNow)
	    this->scheduleLodRefinementFrame("lod-result");
    } else if (partialRefinementCandidate && this->d->lodAutoSubmit &&
	!currentDemandReplayCandidate) {
	/* A cumulative result may extend the shared resident asset yet be
	 * overtaken by the compact source population or current view demand.
	 * Draining it releases the service's coalescing key, so immediately
	 * restart one authoritative pass.  Waiting for an applied presentation
	 * here leaves no worker, frame barrier, or input edge capable of resolving
	 * the current occurrence. */
	restartCurrentDemand();
	this->requestLodCapacityRender("lod-current-demand-replay");
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    uint64_t newRevision = revision == 0 ? 1 : revision;
    if (this->d->lodPolicyRevision.value() == newRevision)
	return;
    this->d->resetLodViewQualityHistory();
    this->d->setAdmissionPolicyRevision(newRevision);
    this->d->lodStaticQualityTrial.reset();
    if (!this->synchronizeAutomaticLodControl()) {
	this->requestPresentationRender("lod-policy");
	return;
    }
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
    this->requestLodCapacityRender("lod-policy");
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
    const BObolViewLodState *state =
	this->d->viewAttachment->getViewLodState();
    if (!state)
	return 0;
    const size_t payloadCount = state->cadPayloadCount();
    const size_t residentProgressiveCount =
	state->residentCadProgressiveCount();
    return residentProgressiveCount > SIZE_MAX - payloadCount ? SIZE_MAX :
	payloadCount + residentProgressiveCount;
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    const bool automaticLod =
	this->automaticLodControlEnabled() != FALSE;
    this->d->advanceAdmissionRevision(
	BObolLodAdmissionRevisionDomain::VIEW);
    this->d->lodInterruptedPresentationReplay.retire();
    BObolViewLodState *presentationState =
	this->d->viewAttachment->getViewLodState();
    if (automaticLod && presentationState &&
	presentationState->hasCadPresentationAssemblies())
	this->d->requireExactPresentationFrame();
    else
	this->d->lodExactPresentationFrame.reset();
    this->d->resetRetainedAdmissionQualityProof();
    this->d->lodStaticQualityTrial.reset();
    if (presentationState) {
	presentationState->clearCadOccurrenceTerminalFailures();
	this->d->publishCadProgressiveCeiling(
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
    this->d->lodRetainedViewContinuity.clearHandoff();
    this->d->lodPlanningObligations.retireImportanceCensus();
    this->d->resetDeadlineSafePresentation();
    if (automaticLod) {
	this->d->lodCoveragePolicy.activate(true);
	this->d->lodViewDemandPolicy.refreshForViewRevision(
	    this->d->lodInteractionSession.active() != FALSE);
    } else {
	/* Camera bookkeeping is shared by all retained views, but a controller
	 * with automatic LoD disabled has no automatic consumer.  Some reset
	 * helpers above also publish exact-presentation or capacity debt when
	 * their renderer values change.  Retire the complete automatic domain
	 * transaction after bookkeeping rather than trying to duplicate its list
	 * here.  The caller's ordinary camera repaint remains authoritative. */
	this->retireAutomaticLodControl();
    }
}

void
BObolViewController::advanceLodPolicyRevision(
    LodPolicyTransition transition)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    this->d->advanceAdmissionRevision(
	BObolLodAdmissionRevisionDomain::POLICY);
    BObolViewLodState *presentationState =
	this->d->viewAttachment->getViewLodState();
    if (presentationState)
	presentationState->clearCadOccurrenceTerminalFailures();
    this->d->lodInterruptedPresentationReplay.retire();
    if (presentationState && presentationState->hasCadPresentationAssemblies())
	this->d->requireExactPresentationFrame();
    else
	this->d->lodExactPresentationFrame.reset();
    this->d->resetRetainedAdmissionQualityProof();
    /* External policy changes start from the preferred quiet cadence.  An
     * internal static-quality successor changes the requested pixel error
     * only after a completed frame has proved the longer deadline; dropping
     * that proof here made the successor immediately coarsen itself. */
    if (transition != LodPolicyTransition::CONTINUE_STATIC_QUALITY)
	this->d->lodStaticQualityTrial.deactivate();
    if (!this->synchronizeAutomaticLodControl()) {
	return;
    }
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);

    /* A policy epoch may change the target cadence (most importantly the
	 * motion-to-quiet transition), refinement authority, and
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
    /* A policy-only transition keeps the camera and source population.  Its
     * existing exact visibility/coverage proof therefore remains current;
     * the policy revision itself starts the quality submission pass.  Marking
     * coverage active here created a second census owner which some terminal
     * capacity paths could retire without completing.  Compaction would then
     * wait forever on that owner even though no census cursor remained. */
    this->d->lodViewDemandPolicy.refreshForPolicyRevision(
	transition != LodPolicyTransition::ORDINARY,
	this->d->lodInteractionSession.active() != FALSE);
}

void
BObolViewController::beginLodInteraction(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    if (!this->automaticLodControlEnabled() ||
	this->d->lodInteractionSession.gestureActive())
	return;

    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::CANCEL_HEADROOM_RETRY);
    this->d->lodStaticQualityTrial.deactivate();
    const float previousPixelError = this->d->lodTargetPixelError;
    const bool newInteractionEpoch = !this->d->lodInteractionSession.active();
    const uint64_t priorCadSampleNanoseconds =
	this->d->lodRendererPerformanceEvidence.cadPresentationNanosecondsAt(
	    this->d->lodPresentationPointProxyPixelThreshold);
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
	    priorCadSampleNanoseconds > 0 &&
	    priorCadSampleNanoseconds <=
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
    this->d->lodRetainedViewContinuity.clearHandoff();
    this->d->lodPlanningObligations.retireImportanceCensus();
    this->d->lodInteractionSession.beginGesture(bu_gettime());
    this->d->lodViewDemandPolicy.beginGesture(newInteractionEpoch);
    if (newInteractionEpoch)
	this->d->lodViewDemandPolicy.seedQualityFloor(initialQualityFloor);
    this->d->lodDiscretePopulationTrialPermit.revoke();
    this->d->lodPresentationPolicy.cancelHandoff();
    this->d->applyAdmissionEvidenceAction(
        BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
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
     * retained cut uncapped despite missing the motion target.  OSMesa reports
     * its work in the traversal time, so the same maximum is portable. */
    const uint64_t interactionTimingSample = std::max(
	this->d->smoothedRenderTimeNanoseconds,
	this->d->lodRendererPerformanceEvidence.lastGpuTimeNanoseconds());
    this->d->lodTargetPixelError =
	BObolLodQualityPolicy::interactivePixelError(
	    interactionTimingSample,
	    this->d->lodInteractiveTargetFps);
    float pointProxyThreshold =
	this->d->lodPresentationPointProxyPixelThreshold;
    if (this->d->pointProxyAggregationApplicableForCameraTransition()) {
	pointProxyThreshold = BObolLodQualityPolicy::pointProxyThreshold(
	    pointProxyThreshold, interactionTimingSample,
	    this->d->lodInteractiveTargetFps);
    } else {
	pointProxyThreshold = BObolViewController::Impl::
	    POINT_PROXY_PIXEL_THRESHOLD_MINIMUM;
	this->d->applyAdmissionEvidenceAction(
	    BObolLodAdmissionPlanner::EvidenceAction::RESET_POINT_PROXY);
    }
	BObolViewLodState *viewState =
	    this->d->viewAttachment->getViewLodState();
	if (viewState)
	    this->d->publishCadCameraMotionFrameReuse(FALSE);

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
    int presentationCeiling = this->d->lodInteractiveProgressiveCeiling;
    if (responsiveCeiling >= 0) {
	if (this->d->lodInteractiveProgressiveCeiling >= 0) {
	    /* A bracketed drag may continue an unbracketed scale epoch.  Keep
	     * that already measured limit and permit only a stricter correction. */
	    if (responsiveCeiling <
		    this->d->lodInteractiveProgressiveCeiling)
		presentationCeiling = responsiveCeiling;
	} else if (this->d->lodTargetPixelError > 1.01f) {
	    /* A new fast pose gesture whose retained cut already meets the
	     * one-pixel contract needs no numerical ceiling equal to that cut. */
	    presentationCeiling = responsiveCeiling;
	}
    }
    this->d->lodInteractionSession.noteCeilingFeedback(
	this->d->renderCompletionSerial);
    this->d->publishCadPresentationLimits(
	presentationCeiling, 0.0f, pointProxyThreshold);
    /* The first pointer motion can precede the first camera signature change.
     * Treat the transition to the motion cut as a policy change so the
     * already resident mesh is coarsened before an expensive drag frame. */
    if (fabsf(this->d->lodTargetPixelError - previousPixelError) >
	    std::numeric_limits<float>::epsilon())
	this->advanceLodPolicyRevision();
    this->markProgressiveWorkPending();
    this->requestLodCapacityRender("lod-interaction-begin");
}

void
BObolViewController::endLodInteraction(void)
{
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
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
    this->requestLodCapacityRender("lod-interaction-end");
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
    this->requestLodCapacityRender("lod-frame-rate");
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
    size_t activeCost = 0;
    size_t minimumCost = 0;
    controller_lod_effective_population_cost(
	state, activeCost, minimumCost);
    return activeCost;
}

BObolLodControlTraceState
BObolViewController::captureLodControlTraceState(void) const
{
    BObolLodControlTraceState state;
    this->getLodConvergenceStatus(state.convergence);
    state.hostWork = this->getHostWorkSnapshot();
    state.viewRevision = this->getLodViewRevision();
    state.policyRevision = this->getLodPolicyRevision();
    state.interactionActive = this->isLodInteractionActive();
    return state;
}

void
BObolViewController::recordLodControlTransition(
    BObolLodControlTransitionEvent event,
    const BObolLodControlTraceState &before,
    const BObolLodControlTraceState &after, SbBool force)
{
    if (!this->d->lodControlTransitionTracing ||
	(!force && controller_lod_control_trace_state_equal(before, after)))
	return;
    if (this->d->lodControlTransitionRecords.size() >=
	this->d->lodControlTransitionRecordLimit) {
	bobol_saturating_counter_advance(
	    this->d->lodControlTransitionDropped);
	return;
    }
    BObolLodControlTransitionRecord record;
    record.serial = bobol_nonzero_identity_take(
	this->d->lodControlTransitionNextSerial);
    record.event = event;
    record.before = before;
    record.after = after;
    this->d->lodControlTransitionRecords.push_back(std::move(record));
}

uint64_t
BObolViewController::beginLodControlTransition(
    BObolLodControlTransitionEvent event, SbBool ownerEvent)
{
    if (!this->d->lodControlTransitionTracing)
	return 0;

    const BObolLodControlTraceState state =
	this->captureLodControlTraceState();
    if (!this->d->lodControlTransitionFrames.empty()) {
	const Impl::LodControlTransitionFrame &parent =
	    this->d->lodControlTransitionFrames.back();
	if (parent.suppressed ||
	    (!parent.ownerEvent && parent.event ==
		BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT)) {
	    Impl::LodControlTransitionFrame frame;
	    frame.token = bobol_nonzero_identity_take(
		this->d->lodControlTransitionNextToken);
	    frame.suppressed = TRUE;
	    frame.state = state;
	    this->d->lodControlTransitionFrames.push_back(frame);
	    return frame.token;
	}
    }
    if (this->d->lodControlTransitionFrames.empty()) {
	const BObolLodControlTransitionEvent pendingEvent =
	    controller_lod_take_pending_transition_event(
		this->d->lodControlPendingExternalEvent);
	if (this->d->lodControlTransitionHasEndpoint &&
	    !controller_lod_control_trace_state_equal(
		this->d->lodControlTransitionEndpoint, state)) {
	    this->recordLodControlTransition(
		pendingEvent,
		this->d->lodControlTransitionEndpoint, state);
	}
    } else {
	Impl::LodControlTransitionFrame &parent =
	    this->d->lodControlTransitionFrames.back();
	const BObolLodControlTransitionEvent parentEvent = parent.ownerEvent ?
	    controller_lod_control_owner_event(parent.state, state) :
	    parent.event;
	this->recordLodControlTransition(parentEvent, parent.state, state);
	parent.state = state;
    }

    Impl::LodControlTransitionFrame frame;
    frame.token = bobol_nonzero_identity_take(
	this->d->lodControlTransitionNextToken);
    frame.event = event;
    frame.ownerEvent = ownerEvent;
    frame.state = state;
    this->d->lodControlTransitionFrames.push_back(frame);
    return frame.token;
}

void
BObolViewController::endLodControlTransition(uint64_t token)
{
    if (!this->d->lodControlTransitionTracing || token == 0 ||
	this->d->lodControlTransitionFrames.empty())
	return;
    Impl::LodControlTransitionFrame frame =
	this->d->lodControlTransitionFrames.back();
    if (frame.token != token) {
	/* A mismatched scope is itself outside the named reducer relation.  Do
	 * not discard the stack: the owning scopes may still close correctly. */
	const BObolLodControlTraceState state =
	    this->captureLodControlTraceState();
	this->recordLodControlTransition(
	    BOBOL_LOD_CONTROL_TRANSITION_UNNAMED, frame.state, state, TRUE);
	return;
    }
    if (frame.suppressed) {
	this->d->lodControlTransitionFrames.pop_back();
	return;
    }

    const BObolLodControlTraceState state =
	this->captureLodControlTraceState();
    const BObolLodControlTransitionEvent event = frame.ownerEvent ?
	controller_lod_control_owner_event(frame.state, state) : frame.event;
    this->recordLodControlTransition(event, frame.state, state);
    this->d->lodControlTransitionFrames.pop_back();
    if (this->d->lodControlTransitionFrames.empty()) {
	this->d->lodControlTransitionEndpoint = state;
	this->d->lodControlTransitionHasEndpoint = TRUE;
    } else {
	this->d->lodControlTransitionFrames.back().state = state;
    }
}

void
BObolViewController::setLodControlTransitionTracing(SbBool enabled,
    size_t recordLimit)
{
    if (!enabled) {
	this->d->lodControlTransitionTracing = FALSE;
	this->d->lodControlTransitionFrames.clear();
	this->d->lodControlTransitionRecords.clear();
	this->d->lodControlTransitionHasEndpoint = FALSE;
	this->d->lodControlPendingExternalEvent.store(
	    BOBOL_LOD_CONTROL_TRANSITION_UNNAMED,
	    std::memory_order_release);
	return;
    }
    if (this->d->lodControlTransitionTracing) {
	this->d->lodControlTransitionRecordLimit = std::max<size_t>(1,
	    recordLimit);
	return;
    }

    this->d->lodControlTransitionRecordLimit = std::max<size_t>(1,
	recordLimit);
    this->d->lodControlTransitionDropped = 0;
    this->d->lodControlTransitionFrames.clear();
    this->d->lodControlTransitionRecords.clear();
    this->d->lodControlPendingExternalEvent.store(
	BOBOL_LOD_CONTROL_TRANSITION_UNNAMED,
	std::memory_order_release);
    this->d->lodControlTransitionTracing = TRUE;
    this->d->lodControlTransitionEndpoint =
	this->captureLodControlTraceState();
    this->d->lodControlTransitionHasEndpoint = TRUE;
    this->recordLodControlTransition(BOBOL_LOD_CONTROL_TRANSITION_INITIAL,
	this->d->lodControlTransitionEndpoint,
	this->d->lodControlTransitionEndpoint, TRUE);
}

SbBool
BObolViewController::isLodControlTransitionTracing(void) const
{
    return this->d->lodControlTransitionTracing;
}

size_t
BObolViewController::drainLodControlTransitions(
    std::vector<BObolLodControlTransitionRecord> &records)
{
    if (this->d->lodControlTransitionTracing &&
	this->d->lodControlTransitionFrames.empty()) {
	const BObolLodControlTraceState state =
	    this->captureLodControlTraceState();
	const BObolLodControlTransitionEvent pendingEvent =
	    controller_lod_take_pending_transition_event(
		this->d->lodControlPendingExternalEvent);
	if (this->d->lodControlTransitionHasEndpoint &&
	    !controller_lod_control_trace_state_equal(
		this->d->lodControlTransitionEndpoint, state)) {
	    this->recordLodControlTransition(
		pendingEvent,
		this->d->lodControlTransitionEndpoint, state);
	}
	this->d->lodControlTransitionEndpoint = state;
	this->d->lodControlTransitionHasEndpoint = TRUE;
    }
    const size_t count = this->d->lodControlTransitionRecords.size();
    records.insert(records.end(),
	std::make_move_iterator(
	    this->d->lodControlTransitionRecords.begin()),
	std::make_move_iterator(
	    this->d->lodControlTransitionRecords.end()));
    this->d->lodControlTransitionRecords.clear();
    return count;
}

uint64_t
BObolViewController::getDroppedLodControlTransitionCount(void) const
{
    return this->d->lodControlTransitionDropped;
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
    status.inventoryRevision = revisionStamp.inventory().value();
    status.availabilityRevision = revisionStamp.availability().value();
    status.visibilityRevision = revisionStamp.visibility().value();
    status.viewRevision = revisionStamp.view().value();
    status.policyRevision = revisionStamp.policy().value();
    status.capacityRevision = revisionStamp.capacity().value();
    const BObolLodCapacitySearchCertificate &capacitySearch =
	this->d->lodAdmissionEvidence.capacity().capacitySearch();
    status.capacitySearchPhase = static_cast<int>(capacitySearch.phase());
    status.capacitySearchGoal = static_cast<int>(capacitySearch.goal());
    status.capacitySearchSamplesRemaining = capacitySearch.samplesRemaining();
    status.capacitySearchMeasuredCandidates =
	capacitySearch.measuredCandidateCount();
    status.capacitySearchTotalMeasuredCandidates =
	capacitySearch.totalMeasuredCandidateCount();
    status.capacitySearchCandidateLimit =
	BObolLodCapacitySearchCertificate::candidateLimit();
    status.capacitySearchMaximumCandidates =
	capacitySearch.maximumCandidateCount();
    status.capacitySearchSampleLimit =
	BObolLodCapacitySearchCertificate::sampleLimit();
    status.capacitySearchCompletedUnits =
	capacitySearch.progressCompletedUnits();
    status.capacitySearchTotalUnits = capacitySearch.progressTotalUnits();
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
    status.reachablePointProxyCandidateCount =
	allocation.reachablePointProxyCandidateCount;
    status.selectedPointProxyCount = allocation.selectedPointProxyCount;
    status.prominentCandidateCount = allocation.prominentCandidateCount;
    status.prominentQualityFloorViolationCount =
	allocation.prominentQualityFloorViolationCount;
    status.maximumNormalizedVisualError =
	allocation.maximumNormalizedError;
    status.visualImportanceDebt = allocation.visualImportanceDebt;
    status.committedAllocationPlanSerial = allocation.allocationPlanSerial;

    const std::vector<SoBRLDatabaseSource *> sources =
	controller_render_database_source_roots(
	    const_cast<BObolViewController *>(this));
    for (const SoBRLDatabaseSource *source : sources) {
	if (!controller_lod_source_has_planning_contract(source))
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
	status.presentedPrimitiveCountValid =
	    viewState->lastCadPresentedPrimitiveCount(
		status.presentedPrimitiveCount);
	status.cadRevision = viewState->cadRevision();
	status.residentDemandRevision =
	    viewState->residentMeshDemandRevision();
	status.currentAllocationPlanSerial = allocation.currentPlanSerial(
	    revisionStamp, viewState->activeCadAllocationPlan());
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
	status.terminalProxyOccurrenceCount =
	    viewState->cadProxyPayloadCount(BOBOL_LOD_PROXY_OBB);
	const BObolViewLodState::CadPresentationPreparationStatus preparation =
	    viewState->cadPresentationPreparationStatus();
	status.rendererPreparationTargetSignature =
	    preparation.targetSignature;
	status.rendererPreparationTotalUnits = preparation.totalUnits;
	status.rendererPreparationCompletedUnits = preparation.completedUnits;
	status.rendererPreparationRemainingUnits = preparation.remainingUnits;
	status.rendererPreparationReservedBytes = preparation.reservedBytes;
	status.rendererPreparationTargetCount = preparation.targetCount;
	status.rendererPreparationPreparingTargetCount =
	    preparation.preparingTargetCount;
	status.rendererPreparationConstrainedTargetCount =
	    preparation.constrainedTargetCount;
	status.rendererPreparationFailedTargetCount =
	    preparation.failedTargetCount;
	status.rendererPreparationInvalidTargetCount =
	    preparation.invalidTargetCount;
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
	this->d->lodAvailabilityLedger.providerPendingCount() > 0 &&
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
    BObolLodControlRefinement::Inputs controlInputs =
	this->d->lodControllerControlInputs();
    /* Coverage is an inventory proof, not compaction state.  Keep its active
     * census in the finite work ledger even when an older retained framebuffer
     * still supplies a useful visual.  Otherwise that stale proof can make a
     * new-camera census look terminal while compaction waits on it forever. */
    controlInputs.inventory = controlInputs.inventory ||
	structuralPending || structuralDiscovery;
    controlInputs.availability = sourcePreparationPending != FALSE;
    controlInputs.result = controlInputs.result || resultPending != FALSE;
    controlInputs.publication =
	controlInputs.publication || publicationPending != FALSE;
    controlInputs.structuralFrontier = controlInputs.structuralFrontier ||
	status.presentedStructuralBoxCount >
	    status.terminalOccurrenceFailureCount;
    controlInputs.capacityFrame = capacityFramePending != FALSE;
    controlInputs.cacheWrite = status.queuedCacheWrites > 0;
    const BObolLodControlRefinement::Snapshot rawControl =
	BObolLodControlRefinement::evaluate(controlInputs);

    /* Result delivery may be background work when the current view is already
     * completely represented.  Exclude the result edge from the convergence
     * ledger here and let BObolLodConvergencePolicy classify it from the
     * complete presentation counts.  Other control facts remain authoritative
     * foreground obligations. */
    BObolLodControlRefinement::Inputs nonResultControlInputs = controlInputs;
    nonResultControlInputs.result = false;
    const BObolLodControlRefinement::Snapshot nonResultControl =
	BObolLodControlRefinement::evaluate(nonResultControlInputs);

    BObolLodAdmissionPlanner::PointCalibrationProducerInputs
	pointProducerInputs;
    pointProducerInputs.submissionPending =
	this->d->lodSubmissionPass.active() != FALSE;
    pointProducerInputs.discoveryCalibrationPending =
	this->d->lodPointAdmissionFrame.pending();
    pointProducerInputs.stableCalibrationPending =
	this->d->lodPointQualityPhase.presentationPending();
    pointProducerInputs.capacityAllocationPending =
	this->d->lodAdmissionEvidence.capacity().capacityAllocationPending();
    pointProducerInputs.capacitySamplePending =
	this->d->lodAdmissionEvidence.capacity().presentationFramePending();
    const BObolViewLodState *witnessPresentationState =
	this->d->viewAttachment ?
	    this->d->viewAttachment->getViewLodState() : NULL;
    pointProducerInputs.stablePresentationAvailable =
	witnessPresentationState &&
	witnessPresentationState->hasCadPresentationAssemblies();
    pointProducerInputs.providerPending =
	this->d->lodAvailabilityLedger.providerPendingCount() > 0;
    pointProducerInputs.servicePending = status.pendingTasks > 0 ||
	status.inFlight > 0 || status.queuedResults > 0;
    pointProducerInputs.publicationAwaitingFrameRequest =
	this->d->lodPresentationTransaction.
	    publicationAwaitingFrameRequest();
    const bool pointProducerOwnsFrame =
	BObolLodAdmissionPlanner::pointProducerOwnsCalibrationFrame(
	    pointProducerInputs);
    /* The composed lifecycle model permits a presentation obligation to move
     * from PUMP to RENDER.  Qualify the shared pump with the controller-local
     * presentation projection: a provider-only wakeup is not evidence that
     * this owner can advance. */
    BObolLodControlRefinement::PresentationProgress presentationProgress;
    presentationProgress.renderPending = hostWork.renderPending();
    presentationProgress.controllerPumpPending = hostWork.pumpPending() &&
	this->d->lodControllerPresentationPumpPending(
	    hostWork.renderPending());
    presentationProgress.finiteTimerPending =
	this->d->lodPresentationTransaction.
	    publicationAwaitingFrameRequest();
    presentationProgress.independentProducerPending =
	this->d->lodPointQualityPhase.presentationPending() &&
	pointProducerOwnsFrame;
    presentationProgress.claimedFramePending =
	hostWork.frameClaimed() != FALSE;
    const bool presentationProgressWitness =
	presentationProgress.witnessed();
    const SbBool calibrationPending =
	rawControl.calibrationPending() ? TRUE : FALSE;
    status.refinementFramePending =
	(this->d->lodPresentationTransaction.barrierPending() ||
	 capacityFramePending) ? TRUE : FALSE;
    status.activeGeneration = this->d->lodActiveGeneration;
    status.submissionSourceIndex = this->d->lodSubmissionSourceIndex;
    status.submissionEntryOffset = this->d->lodSubmissionEntryOffset;
    status.budgetCalibrationPending =
	this->d->lodAdmissionEvidence.capacity().capacityTransactionPending() ||
	this->d->lodAdmissionEvidence.headroom().retryPending();
    status.stablePresentationHandoffPending =
	this->d->lodPresentationPolicy.handoffPending() ? TRUE : FALSE;
    status.pointProxyCalibrationPending =
	(this->d->lodPointAdmissionFrame.pending() ||
	 this->d->lodPointQualityPhase.pending() ||
	 this->d->lodStructuralRepair.pointRelaxationPending()) ? TRUE : FALSE;
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
    for (const BObolProgressiveProviderRecord &provider :
	 this->d->progressiveProviders) {
	status.sourcePreparationCompletedUnits =
	    controller_lod_saturating_add_u64(
		status.sourcePreparationCompletedUnits,
		provider.sourcePreparationCompletedUnits);
	status.sourcePreparationTotalUnits = controller_lod_saturating_add_u64(
	    status.sourcePreparationTotalUnits,
	    provider.sourcePreparationTotalUnits);
    }

    BObolLodConvergencePolicy::Inputs convergenceInputs;
    convergenceInputs.viewEpoch = this->d->lodViewRevision;
    convergenceInputs.policyEpoch = this->d->lodPolicyRevision;
    convergenceInputs.enabled = this->automaticLodControlEnabled() != FALSE;
    convergenceInputs.expectedLeafCount = status.expectedLeafCount;
    convergenceInputs.availableLeafCount = status.availableLeafCount;
    convergenceInputs.visibleTargetCount = status.visibleTargetCount;
    convergenceInputs.activePayloadCount = status.activePayloadCount;
    convergenceInputs.satisfiedPayloadCount = status.satisfiedPayloadCount;
    convergenceInputs.presentedSubpixelOccurrenceCount =
	status.presentedSubpixelOccurrenceCount;
    convergenceInputs.presentedStructuralBoxCount =
	status.presentedStructuralBoxCount;
    convergenceInputs.terminalProxyOccurrenceCount =
	status.terminalProxyOccurrenceCount;
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
    convergenceInputs.sourcePreparationCompletedUnits =
	status.sourcePreparationCompletedUnits;
    convergenceInputs.sourcePreparationTotalUnits =
	status.sourcePreparationTotalUnits;
    convergenceInputs.submissionPending =
	this->d->lodSubmissionPass.active() != FALSE ||
	this->d->lodViewDemandPolicy.demandRefreshActive();
    convergenceInputs.resultPending = resultPending != FALSE;
    convergenceInputs.publicationPending = publicationPending != FALSE;
    convergenceInputs.calibrationPending = calibrationPending != FALSE;
    convergenceInputs.capacitySearchCompletedUnits =
	status.capacitySearchCompletedUnits;
    convergenceInputs.capacitySearchTotalUnits =
	status.capacitySearchTotalUnits;
    /* The refinement ledger is the authoritative foreground-work projection.
     * Individual fields above select the user-facing phase and progress
     * denominator; this aggregate guard prevents a newly added work fact from
     * being accidentally omitted from terminal readiness. */
    convergenceInputs.controlPending = nonResultControl.foregroundPending();
    convergenceInputs.interactive = this->d->lodInteractionSession.active() != FALSE;
    convergenceInputs.compactionPending =
	this->d->lodCompactionPolicy.pending();
    /* Resident reclamation is a level-triggered producer input.  The worker
     * publishes the new admission epoch before its queued owner-thread wake
     * can run.  Include that atomic mismatch directly in readiness so a
     * waiter cannot observe a transient READY frame between reclamation and
     * the sparse retry pass. */
    const bool residentAdmissionPending =
	this->automaticLodControlEnabled() && this->d->lodService &&
	this->d->lodService->residentMeshAdmissionRevision() !=
	    this->d->lodResidentAdmissionRevision.load(
		std::memory_order_acquire);
    convergenceInputs.progressiveWorkPending =
	this->hasProgressiveWorkPending() != FALSE ||
	residentAdmissionPending;
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
    const bool allocationPresentationRealized = viewState &&
	this->d->retainedAllocationPresentationRealized(viewState);
    convergenceInputs.pixelDemandPresentationProven =
	status.currentAllocationPlanSerial != 0 &&
	status.visibleTargetCount > 0 &&
	allocationPresentationRealized && allocation.selectsPixelDemand();
    const bool allocationConstraintPresented =
	convergenceInputs.stableBudgetLimited &&
	allocationPresentationRealized;
    const bool staticDeadlineConstraintPresented =
	this->d->lodStaticQualityTrial.rejectedFor(revisionStamp);
    convergenceInputs.constrainedPresentationProven =
	allocationConstraintPresented || staticDeadlineConstraintPresented ||
	status.terminalProxyOccurrenceCount > 0;
    convergenceInputs.presentationLimited =
	!this->d->forceTerminalLodRefinement &&
	(this->d->lodInteractiveProgressiveCeiling >= 0 ||
	 this->d->lodPresentationPointProxyPixelThreshold > 1.01f ||
	 /* A rejected static-quality trial is completed-frame proof that the
	  * next richer population missed this view epoch's hard presentation
	  * deadline.  Its temporary renderer ceiling may already have been
	  * reconciled and removed, but that must not make the accepted terminal
	  * cut appear unconstrained to the HUD or qualification harness. */
	 this->d->lodStaticQualityTrial.capacityRejected() ||
	 status.terminalProxyOccurrenceCount > 0);
    const BObolLodConvergencePolicy::Decision convergence =
	this->d->lodConvergencePolicy.evaluate(convergenceInputs);

    /* A terminal decision proves that an outstanding result cannot be needed
     * by the current presentation.  Report only the remaining foreground
     * control facts in that state; queuedResults/backgroundPending retain the
     * background-delivery diagnostic. */
    const BObolLodControlRefinement::Inputs &effectiveControlInputs =
	convergence.terminal ? nonResultControlInputs : controlInputs;
    const BObolLodControlRefinement::Snapshot effectiveControl =
	convergence.terminal ? nonResultControl : rawControl;

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
	    BOBOL_LOD_CONTROL_VIOLATION_INVALID_OWNER &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::UNWITNESSED_PRESENTATION) ==
	    BOBOL_LOD_CONTROL_VIOLATION_UNWITNESSED_PRESENTATION &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::UNWITNESSED_CONSTRAINT) ==
	    BOBOL_LOD_CONTROL_VIOLATION_UNWITNESSED_CONSTRAINT &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::UNWITNESSED_PLANNING) ==
	    BOBOL_LOD_CONTROL_VIOLATION_UNWITNESSED_PLANNING &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Violation::
		NONTERMINAL_WITHOUT_PROGRESS) ==
	    BOBOL_LOD_CONTROL_VIOLATION_NONTERMINAL_WITHOUT_PROGRESS,
	"public and private LoD refinement violations must agree");
    static_assert(
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::PresentationWitness::RENDER) ==
	    BOBOL_LOD_PRESENTATION_WITNESS_RENDER &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::PresentationWitness::CONTROLLER_PUMP) ==
	    BOBOL_LOD_PRESENTATION_WITNESS_CONTROLLER_PUMP &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::PresentationWitness::TIMER) ==
	    BOBOL_LOD_PRESENTATION_WITNESS_TIMER &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::PresentationWitness::
		INDEPENDENT_PRODUCER) ==
	    BOBOL_LOD_PRESENTATION_WITNESS_INDEPENDENT_PRODUCER &&
	BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::PresentationWitness::CLAIMED_FRAME) ==
	    BOBOL_LOD_PRESENTATION_WITNESS_CLAIMED_FRAME,
	"public and private presentation witnesses must agree");
    status.phase = static_cast<int>(convergence.phase);
    status.outcome = static_cast<int>(convergence.outcome);
    status.controlFactMask =
	BObolLodControlRefinement::factMask(effectiveControlInputs);
    status.controlObligationMask = effectiveControl.obligations;
    status.controlOwner = static_cast<int>(effectiveControl.owner);

    if (this->d->lodAdmissionEvidence.capacity().stableBudgetLimited())
	status.constraintEvidenceMask |= BOBOL_LOD_CONSTRAINT_STABLE_BUDGET;
    if (this->d->lodInteractiveProgressiveCeiling >= 0)
	status.constraintEvidenceMask |=
	    BOBOL_LOD_CONSTRAINT_PROGRESSIVE_CEILING;
    if (this->d->lodPresentationPointProxyPixelThreshold > 1.01f)
	status.constraintEvidenceMask |=
	    BOBOL_LOD_CONSTRAINT_SUBPIXEL_AGGREGATION;
    if (this->d->lodStaticQualityTrial.capacityRejected())
	status.constraintEvidenceMask |= BOBOL_LOD_CONSTRAINT_STATIC_DEADLINE;
    if (convergence.memoryLimited)
	status.constraintEvidenceMask |= BOBOL_LOD_CONSTRAINT_MEMORY;
    if (status.terminalProxyOccurrenceCount > 0)
	status.constraintEvidenceMask |= BOBOL_LOD_CONSTRAINT_TERMINAL_PROXY;
    status.controlViolationMask = BObolLodControlRefinement::validate(
	effectiveControl, convergence.terminal, convergence.viewReady,
	convergence.terminalError, presentationProgressWitness,
	convergence.performanceLimited,
	status.constraintEvidenceMask != BOBOL_LOD_CONSTRAINT_NONE);
    const bool externalProgressWitness = status.pendingTasks > 0 ||
	status.inFlight > 0 || status.queuedResults > 0 ||
	hostWork.renderPending();
    status.controlViolationMask |=
	BObolLodControlRefinement::validateLiveness(
	    effectiveControl, convergence.terminal, convergence.hasLodState,
	    externalProgressWitness);
    status.controlPresentationWitnessMask =
	presentationProgress.witnessMask();
    status.controlViolationMask |=
	BObolLodControlRefinement::validateProducers(effectiveControlInputs);
    status.fraction = convergence.fraction;
    status.terminal = convergence.terminal ? TRUE : FALSE;
    status.terminalError = convergence.terminalError ? TRUE : FALSE;
    status.viewReady = convergence.viewReady ? TRUE : FALSE;
    status.hasLodState = convergence.hasLodState ? TRUE : FALSE;
    status.backgroundPending =
	convergence.backgroundPending ? TRUE : FALSE;
    status.semanticPresentationFramePending =
	effectiveControl.owner == BObolLodControlRefinement::Owner::PRESENTATION &&
	effectiveControl.obligations == BObolLodControlRefinement::bit(
	    BObolLodControlRefinement::Work::PRESENTATION) &&
	!status.budgetCalibrationPending &&
	!status.pointProxyCalibrationPending &&
	!status.publicationFramePending &&
	!status.backgroundPending ? TRUE : FALSE;
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
    BObolLodControlTransitionScope controlTransition(
	this, BOBOL_LOD_CONTROL_TRANSITION_EXTERNAL_INPUT);
    struct bv_view_info view = BV_VIEW_INFO_INIT;
    SbBool haveCamera = this->getViewInfo(&view);
    const controller_lod_view_signatures signatures =
	controller_lod_view_signature(
	view, haveCamera, this->d->activeCamera, this->d->viewportRegion);

    const SbBool priorViewSignatureValid =
	this->d->lodViewSignature.has_value() ? TRUE : FALSE;
    const BObolLodViewScaleSnapshot priorScaleSignature =
	this->d->lodViewScaleSignature;
    const uint64_t priorCadSampleNanoseconds =
	this->d->lodRendererPerformanceEvidence.cadPresentationNanosecondsAt(
	    this->d->lodPresentationPointProxyPixelThreshold);
    /* advanceLodViewRevision() clears the view-local visibility census.  Its
     * transient zero must not discard an already proven aggregate point cut;
     * sample this presentation fact before invalidating the census. */
    const SbBool priorPointProxyAggregationApplicable =
	this->d->pointProxyAggregationApplicableForCameraTransition() ?
	    TRUE : FALSE;
    SbBool priorViewReady = FALSE;
    if (this->automaticLodControlEnabled() &&
	!this->d->lodInteractionSession.active()) {
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
		priorCadSampleNanoseconds);
	this->d->lodDiscretePopulationTrialPermit.revoke();
	const SbBool retainPriorQualityCeiling =
	    demandChange.retainPriorQualityCeiling ? TRUE : FALSE;
	this->advanceLodViewRevision();
	/* A camera may be edited through its public Coin fields rather than
	 * through syncCameraFromViewContext().  LoD submission is then the first
	 * place that observes the new signature.  Request the view frame here;
	 * previously isRenderRequested() happened to return true merely because
	 * progressive work existed, masking this missing presentation edge. */
	this->requestLodCapacityRender("lod-view");
	/* Camera bookkeeping is unconditional, but the 150 ms refinement pump
	 * belongs only to an active automatic LoD consumer.  Marking generic
	 * controllers progressive here leaves ordinary retained views
	 * permanently render-pending and suppresses later edge-triggered host
	 * wakeups. */
	if (this->automaticLodControlEnabled()) {
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
		priorCadSampleNanoseconds > 0 &&
		priorCadSampleNanoseconds <=
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
		    priorCadSampleNanoseconds > 0 &&
		    priorCadSampleNanoseconds <=
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
	    this->d->lodRetainedViewContinuity.clearHandoff();
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
	        BObolLodAdmissionPlanner::EvidenceAction::INVALIDATE_CAPACITY_MEASUREMENT);
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
		continuingInteractive ? priorCadSampleNanoseconds : 0;
	    const uint64_t interactiveRenderSample = std::max(
		endpointRenderSample,
		this->d->lodRendererPerformanceEvidence.
		    lastGpuTimeNanoseconds());
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
		this->d->publishCadCameraMotionFrameReuse(TRUE);
	    int responsiveCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    const bool hasNewRenderFeedback =
		this->d->lodInteractionSession.hasNewCeilingFeedback(
		    this->d->renderCompletionSerial);
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
	    float pointProxyThreshold =
		this->d->lodPresentationPointProxyPixelThreshold;
	    if (hasNewRenderFeedback && !continuingInteractive &&
		!reusePriorScalePresentation &&
		pointProxyAggregationForCameraChange)
		pointProxyThreshold = BObolLodQualityPolicy::pointProxyThreshold(
			pointProxyThreshold,
			interactiveRenderSample,
			this->d->lodInteractiveTargetFps);
	    else if (!pointProxyAggregationForCameraChange) {
		pointProxyThreshold = BObolViewController::Impl::
		    POINT_PROXY_PIXEL_THRESHOLD_MINIMUM;
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
		this->d->lodInteractionSession.noteCeilingFeedback(
		    this->d->renderCompletionSerial);
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
	     * sample, not only on entry to the interaction.  Otherwise the motion
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
	    int presentationCeiling =
		this->d->lodInteractiveProgressiveCeiling;
	    if (hasNewRenderFeedback && responsiveCeiling >= 0 &&
		(this->d->lodInteractiveProgressiveCeiling < 0 ||
		 responsiveCeiling <
		    this->d->lodInteractiveProgressiveCeiling))
		presentationCeiling = responsiveCeiling;
	    this->d->publishCadPresentationLimits(
		presentationCeiling, 0.0f, pointProxyThreshold);
	    this->markProgressiveWorkPending();
	}
    }
}
